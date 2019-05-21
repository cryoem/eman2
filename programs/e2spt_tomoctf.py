#!/usr/bin/env python
# Muyuan Chen 2017-03
from __future__ import print_function
from __future__ import division
from past.utils import old_div
from builtins import range
import numpy as np
from EMAN2 import *
import json
from scipy.signal import argrelextrema
from multiprocessing import pool, Manager

def calc_ctf(defocus, bxsz=256, voltage=300, cs=4.7, apix=1. ,phase=0.):
	ds=1.0/(apix*bxsz)
	ctf=EMAN2Ctf()
	ctf.from_dict({"defocus":1.0,"voltage":voltage,"bfactor":0.,"cs":cs,"ampcont":0,"apix":apix,"dfdiff":0,"defang":0})
	ctf.set_phase(phase*np.pi/180.)
	cv=[]
	try: l=len(defocus)
	except: defocus=[defocus]
	
	for df in defocus:
		ctf.defocus=df
		cv.append(ctf.compute_1d(bxsz,ds,Ctf.CtfType.CTF_INTEN))
	return np.array(cv)


def get_xf_pos(tpm, pk):
	### first project the point to xy plane
	xf0=Transform({"type":"xyz","xtilt":float(tpm[4]),"ytilt":float(tpm[3])})
	p0=[pk[0], pk[1], pk[2]]
	p1=xf0.transform(p0)#).astype(int)
	
	### undo the 2D alignment
	xf1=Transform({"type":"2d","tx":tpm[0], "ty":tpm[1],"alpha":tpm[2]})
	p2=xf1.transform([p1[0], p1[1]])
	

	return [p2[0], p2[1], p1[2]]


#### Calculate matching score of a set of power spectrum curves (curve) and a set of 1D CTF curves (allctf)
#### and return the a matrix in which each row corresponds to a ctf curve and each column represent a poser spectrum
#### It does the background subtraction in the same way as e2ctf
def calc_all_scr(curve, allctf, zlist, bxsz, exclude=[]):
	#print curve
	allscr=np.zeros((len(allctf), len(curve)))+np.inf

	for i,cf in enumerate(allctf):
		if i in exclude: continue
		zz=zlist[i]
		if len(zz)==0:
			allscr[i]=1
			continue

		z0=zz[0]
		z1=np.minimum(zz[-1], int(bxsz/2*.7))

		if z1-z0<bxsz/10: continue

		bg=np.array([np.interp(np.arange(z0, z1), zz, p[zz]) for p in curve])

		
		bsub=curve.copy()
		bsub[:, z0:z1]-=bg
		for iz in range(1, len(zz)):
			if zz[iz-1]>=z1: break
			c=bsub[:,zz[iz-1]:zz[iz]]
			mx=np.max(c, axis=1)
			m0=(mx<=0)
			mx[m0]=1
			mx=1./mx
			mx[m0]=0
			c*=mx[:, None]

		bsub=bsub[:, z0:z1]
		bsub[bsub<0]=0
		scr=-np.dot(bsub,cf[z0:z1])/(np.sum(cf[z0:z1]))
#		 scr=np.mean((bsub-cf[z0:z1])**2, axis=1)
		allscr[i]=scr
	return allscr



def main():
	
	usage=" "
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_pos_argument(name="tiltseries",help="Specify tiltseries you want to apply CTF correction.", default="", guitype='filebox', browser="EMTiltseriesTable(withmodal=True,multiselect=True)", row=0, col=0,rowspan=1, colspan=3, mode="model")
	parser.add_argument("--alltiltseries", action="store_true",help="Use all tilt series in the folder. Acceptable file extensions include hdf, mrc, mrcs, st.", default=False,guitype='boolbox',row=1, col=0, rowspan=1, colspan=1,mode="model")

	parser.add_header(name="orblock1", help='Just a visual separation', title="Options", row=2, col=1, rowspan=1, colspan=1, mode="model")
	parser.add_argument("--dfrange", type=str,help="Search range of defocus (start, end, step). default is 2., 7, 0.1", default="2.,7.,0.1", guitype='strbox',row=4, col=0,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--psrange", type=str,help="phase shift range (start, end, step). default is 0, 5, 5", default="0,5,5", guitype='strbox',row=4, col=1,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--tilesize", type=int,help="Size of tile to calculate FFT, default is 256", default=256, guitype='intbox',row=4, col=2,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--voltage", type=int,help="Voltage of microscope in kV", default=300, guitype='intbox',row=6, col=0,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--cs", type=float,help="Cs of microscope in kV", default=2.7, guitype='floatbox',row=6, col=1,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--nref", type=int,help="Using N tilt images near the center tilt to estimate the range of defocus for all images. Default is 15", default=15, guitype='intbox',row=6, col=2,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--stepx", type=int,help="Number of tiles to generate on x-axis (different defocus)", default=20, guitype='intbox',row=8, col=0,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--stepy", type=int,help="Number of tiles to generate on y-axis (same defocus)", default=40, guitype='intbox',row=8, col=1,rowspan=1, colspan=1, mode="model")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	
	(options, args) = parser.parse_args()
	logid=E2init(sys.argv)
	
	#### deal with multiple inputs
	if options.alltiltseries:
		fld="tiltseries/"
		args=[fld+f for f in os.listdir(fld) if (
			f.endswith(".hdf") or f.endswith(".mrc") or f.endswith(".mrcs") or f.endswith(".st"))]
	
	if len(args)==1:
		print("Reading tilt series {}...".format(args[0]))
	else:
		print("Processing {} tilt series in sequence..".format(len(args)))
		cmd=sys.argv
		opt=' '.join([s for s in cmd if s.startswith("-")])
		opt=opt.replace("--alltiltseries","")
		for a in args:
			run("{} {} {}".format(cmd[0], a, opt))
		E2end(logid)
		return
	
	tfile=args[0]
	try:
		js=js_open_dict(info_name(tfile))
		tltparams=np.array(js["tlt_params"])
		js=None
	except:
		print("Cannot find alignment parameters. Please make sure tomogram reconstruction is done for selected tilt series")
		return

  
	img=EMData(tfile,0)
	apix=img["apix_x"]
	upix=apix/10./1000.
	nx, ny=img["nx"], img["ny"]
	if img["nz"]>1:
		imgs=[img.get_clip(Region(0, 0, i, img["nx"], img["ny"], 1)).copy() for i in range(img["nz"])]
	else:
		imgs=EMData.read_images(tfile)
		
	for m in imgs: m.process_inplace("normalize")
	nz=len(imgs)
	
	box=options.tilesize
	sz=min(nx, ny)-box
	print("Reading {}, size {} x {} x {}, apix {:.2f}".format(tfile, nx, ny, nz, apix))


	drg=[float(i) for i in options.dfrange.split(',')]
	defrg=np.arange(drg[0],drg[1],drg[2])
	defstep=drg[2]
	
	prg=[float(i) for i in options.psrange.split(',')]
	pshift=np.arange(prg[0], prg[1], prg[2])
	
	print("Search defocus from {:.1f} to {:.1f}, step {:.2f}".format(drg[0], drg[1], drg[2]))
	print("Search phase shift from {:.1f} to {:.1f}, step {:.1f}".format(prg[0], prg[1], prg[2]))
	
	allctf=[]
	zlist=[]
	for ps in pshift:
		ctf=calc_ctf(defrg, box, voltage=options.voltage, cs=options.cs, apix=apix, phase=ps)
		zeros=np.array(argrelextrema(ctf, np.less, axis=1))
		zl=[zeros[1, zeros[0,:]==i].astype(int) for i in range(len(ctf))]
		allctf.append(ctf)
		zlist.append(zl)


	npt=options.stepy
	nstep=options.stepx
	powerspecs=[]
	print("Generating power spectrum of tiles from tilt series...")
	for imgi in range(nz):
		rawimg=imgs[imgi]
		tpm=tltparams[imgi].copy()
		xrg=old_div(sz,2)/np.cos(tpm[3]/180.*np.pi)*.9
		allrd=[]
		pzus=[]
		xstep=old_div(xrg,float(nstep))
		for xpos in np.arange(-xrg, xrg+1, xstep):

			pts=np.zeros((npt, 3))
			pts[:,0]=xpos
			pts[:,1]=(np.random.rand(npt)-.5)*sz*.9

			ptsxf=np.array([get_xf_pos(tpm, p) for p in pts])
			pz=ptsxf[:,2].copy()
			ptsxf=ptsxf[:,:2]
			ptsxf+=[nx//2, ny//2]
			
			ptsxf=ptsxf[np.sum((ptsxf>[nx, ny]) + (ptsxf<0), axis=1)==0, :]
			if len(ptsxf)<npt*.8:
				continue
			
			rds=[]
			for p in ptsxf:
				tile=rawimg.get_clip(Region(p[0]-box//2, p[1]-box//2, box, box))
				if tile["sigma"]<1e-3 or tile["mean"] != tile["mean_nonzero"]:  # strong criteria to exclude edges
					#rds.append(np.zeros(box//2))
					continue
				tile.do_fft_inplace()
				rd=np.array(tile.calc_radial_dist(box//2, 0,1,0))
				rd=np.log10(rd)
				rd-=np.min(rd)
				rd/=np.max(rd[1:])
				rd[rd>1]=1
				rds.append(rd)
			
			if len(rds) == 0: rd = np.zeros(box//2)
			else: rd=np.mean(rds, axis=0)
			allrd.append(rd)
			pzus.append(np.mean(pz*upix))

		allrd=np.array(allrd)
		powerspecs.append([allrd, pzus])
		
		print("{}, {:.2f}".format(imgi, tpm[3]))
		
		
	tltsrt=np.argsort(abs(tltparams[:,3]))
	nref=options.nref
	tltrefs=tltsrt[:nref]
	print("Using first {} images near center tilt to estimate defocus range...".format(nref))
	dfs=[]
	exclude=[]
	ctfparams=np.zeros((len(tltparams), 2))
	for it in tltsrt:
		allrd, pzus=powerspecs[it]
		
		allscr=[]
		for ic, ctf in enumerate(allctf):
			scr=calc_all_scr(allrd, ctf, zlist[ic], box, exclude)
			idxsft=np.round(old_div(-np.array(pzus),(defstep))).astype(int)
			stilt=np.zeros(len(defrg))+np.inf
			for i,df in enumerate(defrg):
				idx=idxsft+i
				outb=(idx<0)+(idx>=len(defrg))
				idx[outb]=0
		
				s=scr[idx, np.arange(scr.shape[1])].copy()
				sinf=np.isinf(s)
				s[outb]=np.max(s)
				sval=s[np.isinf(s)==0]

				stilt[i]=old_div(np.sum(sval),len(s))
			
			allscr.append(stilt)
			
		allscr=np.array(allscr)
		#print(allscr)
		amp, df= np.array(np.where(allscr==np.min(allscr))).T[0]
		print("ID {}, angle {:.1f}, defocus {:.2f}, phase shift {:.0f}".format(it, tltparams[it,3],  defrg[df], pshift[amp]))
		dfs.append(defrg[df])
		ctfparams[it]=[pshift[amp], defrg[df]]

		#### get enough references
		if len(dfs)==nref: 
			print("In the first {} image, defocus mean {:.2f}, std {:.2f}".format(nref, np.mean(dfs), np.std(dfs)))
			#if (np.mean(dfs)+np.std(dfs)*3>defrg[int(len(defrg)*.9)]):
			if (np.std(dfs)>1):
				print ("Warning: variance of defocus estimation is larger than normal. Something may be wrong...")
				#return
			
			searchrg=min(max(np.std(dfs)*3, 3), 1.5)
			dfsel=abs(defrg-np.mean(dfs))<searchrg
			
			
			print("We will search defocus range {:.2f} to {:.2f} for the rest images.".format(np.min(defrg[dfsel]), np.max(defrg[dfsel])))
			
			exclude=np.where(dfsel==False)[0]
			
			
	js=js_open_dict(info_name(tfile))
	js["defocus"]=ctfparams[:,1].tolist()
	js["phase"]=ctfparams[:,0].tolist()
	js["cs"]=float(options.cs)
	js["voltage"]=float(options.voltage)
	js=None
	

	
	E2end(logid)
	
def run(cmd):
	print(cmd)
	launch_childprocess(cmd)
	
	
if __name__ == '__main__':
	main()

