#!/usr/bin/env python
# Muyuan Chen 2017-03
import numpy as np
from EMAN2 import *
from EMAN2_utils import *
from time import time
import json
from scipy.signal import argrelextrema
from multiprocessing import pool, Manager
try:
	import matplotlib
	matplotlib.use("AGG")
	import matplotlib.pyplot as plt
	pltcolors=["k","b","g","r","m","c","darkblue","darkgreen","darkred","darkmagenta","darkcyan","0.5"]
except:
	print "Matplotlib not available, some output will not be generated"


#### This is a numpy version of EMAN2Ctf.compute_1d, and most code is copied from that function and converted to python. 
#### It can take a numpy array of defocus values and return a numpy matrix of 1D CTF curves
def calc_ctf(defocus, bxsz=256, voltage=300, cs=4.7, apix=1. ,ampcnt=0.):


	b2=bxsz/2
	ds=1.0/(apix*bxsz)
	ns=min(int(floor(.25/ds)),bxsz/2)

	ctfout=np.zeros(b2)
	lbda = 12.2639 / np.sqrt(voltage * 1000.0 + 0.97845 * voltage * voltage)

	g1=np.pi/2.0*cs*1.0e7*pow(lbda,3.0);  
	g2=np.pi*lbda*defocus*10000.0;		 
	acac=np.arccos(ampcnt/100.0);				 

	s=np.arange(b2, dtype=float)*ds
	gam=-g1*(s**4)+np.asarray(np.dot(np.asmatrix(g2).T, np.matrix(s**2)))
	ctfout = (np.cos(gam-acac))**2

	return ctfout


#### Calculate matching score of a set of power spectrum curves (curve) and a set of 1D CTF curves (allctf)
#### and return the a matrix in which each row corresponds to a ctf curve and each column represent a poser spectrum
#### It does the background subtraction in the same way as e2ctf
def calc_all_scr(curve, allctf, zeros, bxsz):
	allscr=np.zeros((len(allctf), len(curve)))+np.inf

	for i,cf in enumerate(allctf):
		zz=zeros[1][zeros[0]==i]

		z0=zz[0]
		z1=np.minimum(zz[-1], bxsz/4)

		if z1-z0<10: continue
		bg=np.array([np.interp(np.arange(z0, z1), zz, p[zz]) for p in curve])
		bsub=curve[:,z0:z1]-bg
		scr=-np.dot(bsub,cf[z0:z1])#/(z1-z0)
		allscr[i]=scr
	return allscr



#### function to get back the index of the allctf array from a defocus value..
def get_ctf_index(x, defmin, defstep):
	return np.array((np.round(x,2)-defmin)/defstep, dtype=int)

#### calculate the defocus value for one tilt 
def calc_defocus_onetlt(args):
	
	#### unpack inputs
	imgi, options=args
	
	stepx=options.stepx
	stepy=options.stepy
	ntry=options.ntry
	bxsz=options.tilesize
	atoum=float(10*1000)
	allctf=options.allctf
	zeros=options.zeros
	tlts=options.tlts
	defrg=options.defrg
	
	
	ang=tlts[imgi]
	rawimg=EMData(options.alifile, 0, False, Region(0,0,imgi,options.nx,options.ny,1))
	apix=rawimg["apix_x"]
	centerx=rawimg["nx"]/2.
	tilex=((np.arange(stepx+1,dtype=float))/(stepx)*(rawimg["nx"]-bxsz)).astype(int)+bxsz/2
	scrtlt=[]
	
	alldefs=np.zeros(ntry)
	for tries in range(ntry):
		psall=[]
		pxlst=np.zeros(stepx)
		
		for ix in range(stepx):
			px=np.random.randint(tilex[ix], tilex[ix+1])
			pxlst[ix]=px
			rrd=[]
			for iy in range(stepy*2): #### make sure we actually have stepy tiles total...
				if len(rrd)>=stepy: break
				py=np.random.randint(rawimg["ny"]-bxsz-2)+bxsz/2-1
				cc=rawimg.get_clip(Region(px-bxsz/2, py-bxsz/2, bxsz, bxsz))
				#### throw away mostly empty tiles
				if cc["sigma"]<rawimg["sigma"]*.5:
					continue
				cc.do_fft_inplace()
				rd=np.array(cc.calc_radial_dist(bxsz/2, 0,1,0))
				rd=np.log10(rd)
				if np.max(rd[1:])<1: 
					print rd
					break
				rd-=np.min(rd)
				rd/=np.max(rd[1:])
				rd[rd>1]=1
				rrd.append(rd)
			psall.append(rrd)

		psavgx=np.array([np.mean(p, axis=0) for p in psall])

		allscr=calc_all_scr(psavgx, allctf, zeros, bxsz)

		

		xdif=centerx-pxlst
		dfdf=apix/atoum*np.sin(ang/180.*np.pi)* xdif
		
		dfs=np.repeat(defrg.copy()[:,None],len(pxlst), axis=1)
		dfspx=dfs+dfdf
		idx=get_ctf_index(dfspx, options.defmin, options.defstep)
		exclude=np.logical_or(np.min(idx, axis=1)<0, np.max(idx, axis=1)>len(allctf)-1)
		idx1=np.clip(idx,0, len(allctf)-1)
		scrs=np.array([allscr[idx1[:, i],i] for i in range(len(pxlst))]).T
		scr= np.mean(scrs, axis=1)
		scr[exclude]=np.max(scr)+.01
		scrtlt.append(scr)
		finaldef=defrg[np.argmin(scr)]
		alldefs[tries]=finaldef
	
	#si=np.argmin(np.mean(np.array(scrtlt), axis=0))
	
	return [imgi, ang, np.mean(alldefs), np.std(alldefs)]


def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--tltfile", type=str,help="IMOD .tlt file", default=None)
	parser.add_argument("--alifile", type=str,help="IMOD .ali file", default=None)
	parser.add_argument("--output", type=str,help="Output text file contains defocus values", default="defocus_output.txt")
	parser.add_argument("--tardef", type=float,help="Target defocus when taking the image (in positive um). The program will search tardef+-defrange for the best defocus", default=5)
	parser.add_argument("--defrange", type=float,help="Search range of defocus. default is 3", default=3)
	parser.add_argument("--tilesize", type=int,help="Size of tile to calculate FFT, default is 256", default=256)
	parser.add_argument("--voltage", type=int,help="Voltage of microscope in kV", default=300)
	parser.add_argument("--cs", type=int,help="Cs of microscope in kV", default=4.7)
	parser.add_argument("--apix", type=float,help="A per pixel", default=-1)
	parser.add_argument("--stepx", type=int,help="Number of tiles to generate on x-axis", default=20)
	parser.add_argument("--stepy", type=int,help="Number of tiles to generate on y-axis", default=20)
	parser.add_argument("--ntry", type=int,help="Measure the CTF n times for each tilt and use the mean defocus value. default is 10", default=10)
	parser.add_argument("--threads", type=int,help="Number of threads to use", default=5)
	parser.add_argument("--ptclin", type=str,help="2D particles input", default=None)
	parser.add_argument("--ptclout", type=str,help="2D particles output", default=None)
	
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	if options.ptclin:
		pname=options.ptclin
		e=EMData(pname, 0, True)
		num=EMUtil.get_image_count(pname)
		bxsz=e["nx"]
		apix=e["apix_x"]
		ctf=e["ctf"]
		if options.voltage>0:
			ctf.voltage=options.voltage
		if options.cs>=0:
			ctf.cs=options.cs
		print "{:d} particles, Voltage={:.1f}kV, Cs={:.1f}".format(num, ctf.voltage, ctf.cs)
		
		if options.ptclout==None:
			fpname=options.ptclin[:options.ptclin.rfind('.')]+"__ctf_flip.hdf"
		else:
			fpname=options.ptclout
		print "Saving output to {}".format(fpname)
		for i in range(num):
			e=EMData(pname, i)
			fft1=e.do_fft()
			defocus=e["defocus"]
			ctf.defocus=defocus
			flipim=fft1.copy()
			ctf.compute_2d_complex(flipim,Ctf.CtfType.CTF_SIGN)
			fft1.mult(flipim)
			out=fft1.do_ift()
			out["ctf"]=ctf
			out["apix_x"] = ctf.apix
			out["apix_y"] = ctf.apix
			out["apix_z"] = ctf.apix
			
			out.write_image(fpname,i)
		
	elif options.tltfile:
	
		time0=time()
		
		tlts=np.loadtxt(options.tltfile)
		print "Read {} tilts from {} to {}.".format(len(tlts), np.min(tlts), np.max(tlts))
		options.tlts=tlts
		
		#### Set up the parameters 
		rawimg=EMData(options.alifile, 0, True)
		options.defmin,options.defmax=options.tardef-options.defrange,options.tardef+options.defrange
		options.nx=rawimg["nx"]
		options.ny=rawimg["ny"]
		
		if options.apix<0:
			apix=rawimg["apix_x"]
			print "Apix from .ali file header: {}".format(apix)
		else:
			apix=options.apix
			
		
		
		options.defstep=.01
		options.defrg=np.arange(options.defmin, options.defmax, options.defstep)
		
		
		#### Generate all possible CTF curves first
		options.allctf=calc_ctf(options.defrg, options.tilesize, voltage=options.voltage, cs=options.cs, apix=apix)
		options.zeros=np.array(argrelextrema(options.allctf, np.less, axis=1))
		
		print "Start working on {} threads...".format(options.threads)

		pl=pool.Pool(options.threads)
		ret=pl.map_async(calc_defocus_onetlt, [(i, options) for i in range(len(tlts))])
		pl.close()
		pl.join()
		
		tltdefs=ret.get()
		print "idx\ttilt angle\tdefocus mean\tdefocue std"
		for t in tltdefs:
			print "{}\t{:.2f}\t{:.2f}\t{:.2f}".format(t[0],t[1],t[2],t[3])
			
			
		info="Input file: {}\nVoltage: {:.1f}\nCs: {:.1f}\n\nidx,tilt angle,defocus mean,std".format(options.alifile,options.voltage, options.cs)
		np.savetxt(options.output,tltdefs, fmt="%.2f", header=info)
		
		print "Done. Output saved to {}. Total time: {}.".format(options.output, time()-time0)
		
		
		try:
			d= np.array(tltdefs)
			plt.figure(figsize=(10,6))
			plt.errorbar(d[:,1], d[:,2], d[:,3], fmt='.-')
			plt.xlim(np.min(d[:,1])-1,np.max(d[:,1])+1)
			plt.xlabel("tilt degree")
			plt.ylabel("defocus/ um")
			plt.savefig(options.output[:options.output.rfind('.')]+"_plot.png")
		except:
			pass
		
	else:
		print "no input..."
	
	E2end(logid)
	
def run(cmd):
	print cmd
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()
	