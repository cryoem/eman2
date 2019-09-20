#!/usr/bin/env python
# Author: Jesus Galaz, 11/01/2012; last update sep/2019
# Copyright (c) 2011 Baylor College of Medicine
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or

# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#
from __future__ import print_function
from __future__ import division
from past.utils import old_div
from builtins import range
import os
from EMAN2 import *
from EMAN2_utils import *
#from time import time


#import matplotlib
#matplotlib.use('Agg',warn=False)		 

#import matplotlib.pyplot as plt
import collections
import sys
import numpy		 
import math
import e2ctf
	 
def main():
	
	progname = os.path.basename(sys.argv[0])
	usage = """Program to fit CTF (find defocus and estimate highest resolution conent) in raw images of tiltseries by periodogram averaging of 
	tiles in strips parallel to the tilt axis."""
		
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--ampcont",type=float, default=0.05,help="""Default=0.05. Amplitude contrast to use for CTF correction phase flipping.""")
	
	parser.add_argument("--apix",type=float, default=None,help="""Default=whatever is on the header of the images. Sampling of the images in angstroms/pixel.""")	
	
	parser.add_argument("--bfactor",type=int, default=1000,help="""Default=1000. Bfactor ("temperature factor") to use.""")

	parser.add_argument("--cs", type=float, default=2.7,help="""Default=2.7. Cs of the microscope with which the images were collected.""")
	
	parser.add_argument("--defocus", type=float, default=None,help="""Default=None. Target defocus at the tilt axis. If not provided, a "global defocus" value will be estimated automatially.""")
	
	parser.add_argument("--debug", action='store_true', default=False,help="""Default=False. Shows plots as they are being produced -requires user feedback to proceed.""")

	parser.add_argument("--input", type=str, default=None, help="""Default=None. Single 2D image or image stack to calculate CTF for.""")
				
	parser.add_argument("--path",type=str, default='tomoctfraw',help="""Directory to store results in. The default is a numbered series of directories containing the prefix 'tomoctfraw'; for example, sptctf_02 will be the directory by default if 'sptctf_01' already exists.""")

	parser.add_argument("--phaseplate", action="store_true", default=False, help="""Default=False. Supply this if the data were collected with hole-free phase plate.""")
	
	parser.add_argument("--plots", action='store_true', default=False,help="""Default=False. Turn this option on to generate plots of the background-subtracted 1D power spectrum overlaid with the best CTF fit. Running on a cluster or via ssh remotely might not support plotting.""")

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	parser.add_argument("--savefft", action="store_true", default=False, help="""Default=False. Saves the average of the ffts of the tiles for an image or its strips.""")

	parser.add_argument("--tilesize",type=int, default=256,help="""Tile size to use for strips when --autofit is provided.""")
		
	parser.add_argument("--tiltangle", type=float, default=None, help="""Default=None. Single 2D image or image stack to calculate CTF for.""")

	parser.add_argument("--tltfile", type=str, default=None, help="""File containing a list of tilt angles corresponding to the tilt angles of images 0 to n.""")
	
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness.")

	parser.add_argument("--voltage", type=int, default=300,help="""Default=300. Voltage of the microscope with which the images where collected.""")
	
	(options, args) = parser.parse_args()		
	

	#'''
	#Make sure input file is EMAN2 readable
	#'''
	options=checkinput(options)
	n=EMUtil.get_image_count(options.input)
	hdr=EMData(options.input,0,True)
	tomo=False
	
	#'''
	#if --input has more than one slice in z, it is likely a tiltseries
	#'''
	nz=hdr['nz']	
	if nz>n:
		tomo=True
		n=nz
	
	#'''
	#Check whether paramters are valid, sensible, workable; if yes, load images and tiltangles into a numbered dictionary
	#'''
	check_healthy_params(options,n)
	angles=get_angles(options,n)

	dirname = os.path.dirname( options.input )
	options.path = dirname + "/" + options.path
	imgs_unstacked = {}
	if tomo:
		options = makepath( options, 'tomoctfraw')
		stem,extension = os.path.splitext(os.path.basename(options.input))
		out = stem + "_unstacked.hdf"

		cmdunstack = "e2proc2d.py " + options.input + ' ' + out + " --unstacking"
		imgs_unstacked_tmp = [ out.replace('.hdf', '-' + str(i+1).zfill(len(str(nz)))) + '.hdf' for i in range(nz) ]

		runcmd(options,cmdunstack)
		for img in imgs_unstacked_tmp:
			out_final =  options.path + '/' + img.replace('unstacked-','unstacked_')
			os.rename( img, out_final )
			num = int(img.split('unstacked-')[-1].replace('.hdf','')) - 1 				#indexes should start from 0
			imgs_unstacked.update({ num : out_final })

	else:
		imgs_unstacked.update({0:options.input})

	#'''
	#Log current run of the program ONLY if input parameters are healthy and images and angles have been loaded without issues
	#'''
	logger = E2init(sys.argv, options.ppid)

	global_ctfs,bad_indexes=get_global_defocuses(options,imgs_unstacked,n,angles,dirname)

	E2end(logger)
	return
	

def get_global_defocuses(options,imgs_unstacked,n,angles,dirname):

	global_defoci = {}
	global_ctfs = {}
	bad_indexes = []
	
	for i in range(n):
		angle=angles[i]
		#img = EMData(options.input,i)
		img_file = imgs_unstacked[i]
		img = EMData(img_file,0)

		apix = img['apix_x']
		if options.apix:
			apix = options.apix
		else:
			options.apix = apix

		nx=img['nx']
		ny=img['ny']
		
		#'''
		#Fit the "global defocus" first, averaging tiles over the entire image, without caring about the defocus gradient direction or the tilt angle
		#'''
		defocusmin = 0.5
		defocusmax = 9.0

		if options.defocus:
			defocusmin = options.defocus - 1.5
			defocusmax = options.defocus + 1.5

		try:
			global_ctf = fit_global_ctf(options,img_file,img,apix,defocusmin,defocusmax)
			global_ctfs.update({i:global_ctf})

			global_defoci.update({i:global_ctf.defocus})
			
			maxres = math.sqrt(global_ctf.bfactor/6.0)
			print("\n(e2tomo_ctfraw)(main) for img {}/{}, global_defocus={}, maxres={}".format(i+1,n,global_ctf.defocus,maxres))

		except:
			print("\nWARNING: img=img_file has a poor fit; retrying")
			if global_defoci:
				try:
					global_defocus_avg = numpy.mean(global_defoci.values())
					global_defocus_std = numpy.std(global_defoci.values())
					defocusmin = global_defocus_avg - 2*global_defocus_std
					defocusmax = global_defocus_avg + 2*global_defocus_std
					print("\nretrying fit with constrained defocusmin={}, defocusmax={}".format(defocusmin,defocusmax))
					global_ctf = fit_global_ctf(options,img_file,img,apix,defocusmin,defocusmax)
					global_ctfs.update({i:global_ctf})
					global_defoci.update({i:global_ctf.defocus})

				except:
					print("\nERROR: skipping img={} due to poor fit; I'll try ONE more time after going over all images".format(img_file))
					bad_indexes.append(i)
					global_ctfs.update({i:None})
			else:
				bad_indexes.append(i)
				continue


	for j in bad_indexes:
		angle=angles[j]
		#img = EMData(options.input,i)
		img_file = imgs_unstacked[j]
		img = EMData(img_file,0)

		apix = img['apix_x']
		if options.apix:
			apix = options.apix
		else:
			options.apix = apix

		nx=img['nx']
		ny=img['ny']

		try:
			global_defocus_avg = numpy.mean(global_defoci)
			global_defocus_std = numpy.std(global_defoci)
			
			defocusmin = global_defocus_avg - 1*global_defocus_std
			defocusmax = global_defocus_avg + 1*global_defocus_std

			print("\nretrying fit for bad image={} with HYPER-constrained defocusmin={}, defocusmax={}".format(img_file,defocusmin,defocusmax))
			
			global_ctf = fit_global_ctf(options,img_file,img,apix,defocusmin,defocusmax)
			global_ctfs.update({j:global_ctf})
			global_defoci.update({j:global_ctf.defocus})

		except:
			print("\nERROR: skipping img={} PERMANENTLY for global defocus fitting due to poor fit; I'll try ONE more time after going over all images".format(img_file))
			global_ctfs.update({j:None})
			#bad_indexes.append(i)

		#except:
		#	print("\nWARNING!!!\n!!!\n!!!: failed to fit CTF for input={}, index={}; the image is likely garbage".format(options.input,i))

			#ctfer(options,i,angle,n)	

	globald_f = 'global_defoci.txt'
	if n > 1: 
		globald_f = options.path + '/' + globald_f
	else:
		globald_f =  dirname + '/global_defoci.txt'

	
	
	#with open(globald_f,'w') as gdf:
	#	lines=[ str(i) + '\t' + str(global_defoci[i]) + '\n' for i in range(len(global_defoci))]
	#	gdf.writelines(lines)
	
	write_txt( global_defoci.keys(), global_defoci.values(), globald_f )
	
	return global_ctfs,bad_indexes


def fit_global_ctf(options,img_file,img,apix,defocusmin,defocusmax):
	#grid_coords=tile_grid(nx,ny,options.tilesize):
	print("\n(fit_global_ctf) entered function")

	verbose=False
	if options.verbose > 5.0:
		verbose=True
	

	tiles = get_tiles(img,options.tilesize,True)
	tiles_fft_avg = incoherent_sum_from_imglist(tiles,checkcorners=True,verbose=verbose)

	if options.savefft:
		extension = os.path.splitext(os.path.basename( img_file ))[-1]
		print("\n(get_global_defocus) extension={}".format(extension))
		out_fft = img_file.replace(extension,'_fft_tile_avg.hdf')
		print("\n(get_global_defocus) out_fft={}".format(out_fft))
		if out_fft == img_file:
			print("\nERROR: risk of overwritting input image")
			sys.exit(1)
		tiles_fft_avg.write_image(out_fft,0)
	
	tiles_fft_avg.process_inplace('filter.highpass.gauss',{"cutoff_pixels":3})
	tiles_fft_avg.set_value_at(0,0,0.0)
	
	print("\n(fit_global_ctf) will call fit_defocus")

	ctf = fit_defocus(options,img_file,tiles_fft_avg, options.voltage, options.cs, options.ampcont, options.phaseplate, apix, options.defocus, defocusmin, defocusmax, 0.1)
	print("\n(fit_global_ctf) leaving function")

	return ctf


def get_angles( options, n ):
	
	angles = {}
	
	if n >1:
		with open( options.tltfile, 'r' ) as f:
			lines = [ float(line.replace('\t','').replace('\n','')) for line in f.readlines()]
			angles = { i : lines[i] for i in range(0, len(lines) ) }

		if angles:
			if len(angles) != n:
				print("\nERROR: the number of angles na={} does not match the number of images in the tiltseries nt={}".format(len(angles),n))
				sys.exit(1)

			if options.verbose > 9:
				print("\n(e2tomo_ctfraw.py)(getangles) angles={}".format(angles) )
			return angles
		else:
			print("\nERROR: failed to read angles from --tltfile={}".format(options.tltfile))
			sys.exit(1)
	elif n==1:
		angles.update({0:options.tiltangle})
		return angles

	return


def check_healthy_params(options,n):
	if n==1 and options.tiltangle == None:
		print("\nERROR: for single images, --tiltangle is required; --tiltangle={}".format(options.tiltangle))
		sys.exit(1)
	elif n>1 and not options.tltfile:
		print("\nERROR: for a full tiltseries, --tltfile is required")
		sys.exit(1)
	elif n>1 and options.tltfile:
		with open(options.tltfile,'r') as f:
			lines=f.readlines()
			nlines=len(lines)
			if len(lines) < n:
				print("\nERROR: insufficient lines n={} in tltfile={} for input={} with n={} images.".format(nlines,options.tltfile,options.input,n))
				sys.exit(1)
			elif len(lines) > n:
				print("\nWARNING: too many lines n={} in tltfile={} for input={} with n={} images. I will attempt to clean --tltfile in case some parasitic lines are confounding the matter".format(nlines,options.tltfile,options.input,n))
				try:
					options.tltfile = remove_blank_lines(options.tltfile,True)
				except:
					print("\nERROR: coudl not clean --tltfile. Please verify that there are the same number of lines as tilt images in --input")
					sys.exit(1)
	return


def fit_defocus(options,img_file,fft_avg,voltage,cs,ampcont,phaseplate,apix,target_defocus,defocusmin=0.5,defocusmax=9.0,defocusstep=0.1):
	"""Fit defocus from an incoherent average of tiles 
	Author: Jesus Montoya, jgalaz@gmail.com, 09/2019
	"""
	print("\n(fit_defocus) entering function")

	tilesize = fft_avg['ny']
	
	fft_bg = fft_avg.process("math.nonconvex")

	fft_avg_1d = fft_avg.calc_radial_dist(old_div(fft_avg.get_ysize(),2),0.0,1.0,1)	# note that this handles the ri2inten averages properly
	
	#print("\n(e2tomo_ctfraw)(fit_defocus) fft_avg_1d={}".format(fft_avg_1d))

	ctf = EMAN2Ctf()

	#print("\n(e2tomo_ctfraw)(fit_defocus) apix={}".format(apix))
	#print("\n(e2tomo_ctfraw)(fit_defocus) tilesize={}".format(tilesize))

	# Compute 1-D curve and background
	ds = old_div(1.0,( apix * tilesize ))
	#print("\n(e2tomo_ctfraw)(fit_defocus) ds={}".format(ds))
	
	bg_1d = e2ctf.low_bg_curve(fft_avg_1d,ds)
	#print("\n(e2tomo_ctfraw)(fit_defocus) bg_1d={}".format(bg_1d))
	
	dfhint=( defocusmin, defocusmax, defocusstep )
	ctf = e2ctf.ctf_fit( fft_avg_1d, bg_1d, bg_1d, fft_avg, fft_bg, float(voltage), float(cs), float(ampcont), phaseplate, float(apix), bgadj=1, autohp=True, dfhint=dfhint)		
	#print("\n(e2tomo_ctfraw)(fit_defocus) ctf={}".format(ctf))

	#bg_sub = numpy.array(fft_avg_1d) - numpy.array(bg_1d)

	bg_sub,bg_1d = e2ctf.calc_1dfrom2d(ctf, fft_avg, fft_bg)

	fit = ctf_get_fit_curve(ctf, ds, fft_avg_1d, bg_1d, bg_sub)
	fit/=max(fit)

	
	#ctf = adjust_ctf_bg(ctf,fft_avg_1d)

	#maxres,bfactor = ctf_fit_bfactor_and_maxres( list(bg_sub),ds,ctf)
	#print("\ndefocus={}, FIRST trial".format(ctf.defocus))


	ctf2 = EMAN2Ctf()
	dfhint2=(max(dfhint[0],ctf.defocus-0.1),min(dfhint[1],ctf.defocus+0.1),min(old_div(dfhint[2],2.0),0.01))
	ctf2 = e2ctf.ctf_fit( fft_avg_1d, bg_1d, bg_1d, fft_avg, fft_bg, float(voltage), float(cs), float(ampcont), phaseplate, float(apix), bgadj=1, autohp=True, dfhint=dfhint2 )		
	bg_sub2,bg_1d2 = e2ctf.calc_1dfrom2d(ctf2, fft_avg, fft_bg)
	
	#bg_sub2,bg_1d2 = e2ctf.calc_1dfrom2d(ctf2, fft_avg, fft_bg)
	fit2 = ctf_get_fit_curve(ctf2, ds, fft_avg_1d, bg_1d2, bg_sub2)
	fit2/=max(fit2)

	#maxres,bfactor = ctf_fit_bfactor_and_maxres( list(bg_sub2),ds,ctf)
	#print("\nmaxres calc={}, bfactor={}, defocus={}, SECOND trial".format(maxres,ctf2.bfactor,ctf2.defocus))

	bg_sub_final = numpy.array(bg_sub2) - numpy.array(bg_1d2)
	bg_sub_final/=max(bg_sub_final)
	import matplotlib.pyplot as plt
	
	if options.plots or options.debug:
		r = len(ctf.background)
		s = numpy.arange(0,ds*r,ds)
		
		if options.debug:
			plt.plot(s, fit2, 'k')
			plt.plot(s, bg_sub_final, 'b')
			plt.show()

		if options.plots:
			extension = os.path.splitext(os.path.basename( img_file ))[-1]

			print("\n(fit_defocus) should save plots")
			ps_file = img_file.replace(extension,'_1d_ps.txt')
			if ps_file == img_file:
				print("\nERROR: risk of overwritting input image")
				sys.exit(1)
			write_txt(s,bg_sub_final,ps_file)
			#with open(options.path + '/' + pfile,'w') as f:
			#	lines = [ str(s[i])+'\t'+str(bg_sub_final[i])+'\n' for i in range(len(s-1))]
			#	lines.append( str(s[-1]) + '\t' + str(bg_sub_final[-1]) )	#the last element does not need a 'return' or 'line break'
			
			fit_file = img_file.replace(extension,'_1d_fit.txt')
			if fit_file == img_file:
				print("\nERROR: risk of overwritting input image")
				sys.exit(1)
			#with open(options.path + '/' + pfile,'w') as f:
			#	lines = [ str(s[i])+'\t'+str(bg_sub_final[i])+'\n' for i in range(len(s-1))]
			#	lines.append( str(s[-1]) + '\t' + str(bg_sub_final[-1]) )	#the last element does not need a 'return' or 'line break'
			write_txt(s,fit2,fit_file)

	print("\n(fit_defocus) leaving function")

	return ctf


def write_txt(xdata,ydata,filename):
	print("\n(write_txt) writing plot txt file for f={}".format(filename))
	with open(filename,'w') as f:
		lines = [ str(xdata[i])+'\t'+str(ydata[i])+'\n' for i in range(len(xdata)) ]
		f.writelines(lines)

	return



def ctf_get_fit_curve(ctf, ds, fft_avg_1d, bg_1d, bg_sub):
	"""From e2ctf.py"""
	r = len(ctf.background)
	s = numpy.arange(0,ds*r,ds)
	
	fit = numpy.array(ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP))		# The fit curve
	fit = fit*fit			# squared

	# auto-amplitude for b-factor adjustment
	rto,nrto=0,0
	for i in range(int(old_div(.04,ds))+1,min(int(old_div(0.15,ds)),len(s)-1)):
		if bg_sub[i]>0 :
			rto+=fit[i]
			nrto+=fabs(bg_sub[i])
	if nrto==0 : rto=1.0
	else : rto/=nrto
	fit=[old_div(fit[i],rto) for i in range(len(s))]

	return fit


'''
def adjust_ctf_bg(ctf,fg_1d):
	"""Smooths the background based on the values of the foreground near the CTF zeroes and puts the
	smoothed background into the CTF object. From e2ctf.py"""
	
	ds=ctf.dsbg
	#print "\n (bgAdj) ds is", ds
	#print "\n (bgAdj) ctf is", ctf
	bg_1d=list(ctf.background)
	#print "\n (bgAdj) bg_1d is", bg_1d

	xyd=XYData()

	# Find the minimum value near the origin, which we'll use as a zero (though it likely should not be)
	mv=(fg_1d[1],1)
	fz=int(old_div(ctf.zero(0),(ds*2)))
	
	for lz in range(1,fz):
		mv=min(mv,(fg_1d[lz],lz))

	xyd.insort(mv[1],mv[0])

	# now we add all of the zero locations to our XYData object
	for i in range(100):
		z=int(old_div(ctf.zero(i),ds))
		if z>=len(bg_1d)-1: break
		if fg_1d[z-1]<fg_1d[z] and fg_1d[z-1]<fg_1d[z+1]: mv=(z-1,fg_1d[z-1])
		elif fg_1d[z]<fg_1d[z+1] : mv=(z,fg_1d[z])
		else : mv=(z+1,fg_1d[z+1])
		xyd.insort(mv[0],mv[1])

	# new background is interpolated XYData
	ctf.background=[xyd.get_yatx_smooth(i,1) for i in range(len(bg_1d))]

	# if our first point (between the origin and the first 0) is too high, we readjust it once
	bs=[fg_1d[i]-ctf.background[i] for i in range(fz)]
	
	#print "bs first is", bs
	
	if min(bs)<0 :
		mv=(bs[0],fg_1d[0],0)
		for i in range(1,fz): mv=min(mv,(bs[i],fg_1d[i],i))
		xyd.set_x(0,mv[2])
		xyd.set_y(0,mv[1])
		
		ctf.background=[xyd.get_yatx_smooth(i,1) for i in range(len(bg_1d))]
		
		#bs2=[fg_1d[i]-ctf.background[i] for i in xrange(fz)]
	return ctf
'''


if '__main__' == __name__:
	main()

