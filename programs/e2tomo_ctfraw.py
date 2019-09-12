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
	
	parser.add_argument("--ampcont",type=float,default=0.05,help="""Default=0.05. Amplitude contrast to use for CTF correction phase flipping.""")
	
	parser.add_argument("--apix",type=float,default=None,help="""Default=whatever is on the header of the images. Sampling of the images in angstroms/pixel.""")	
	
	parser.add_argument("--bfactor",type=int,default=1000,help="""Default=1000. Bfactor ("temperature factor") to use.""")

	parser.add_argument("--cs", type=float,default=2.7,help="""Default=2.7. Cs of the microscope with which the images were collected.""")
	
	parser.add_argument("--defocus", type=float,default=None,help="""Default=None. Target defocus at the tilt axis. If not provided, a "global defocus" value will be estimated automatially.""")

	parser.add_argument("--input", type=str, default=None, help="""Default=None. Single 2D image or image stack to calculate CTF for.""")
				
	parser.add_argument("--path",type=str,default='sptctf',help="""Directory to store results in. The default is a numbered series of directories containing the prefix 'sptctfraw'; for example, sptctf_02 will be the directory by default if 'sptctf_01' already exists.""")

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness.")

	parser.add_argument("--tilesize",type=int,default=256,help="""Tile size to use for strips when --autofit is provided.""")
		
	parser.add_argument("--tiltangle", type=str, default=None, help="""Default=None. Single 2D image or image stack to calculate CTF for.""")

	parser.add_argument("--tltfile",default=None,type=str,help="""File containing a list of tilt angles corresponding to the tilt angles of images 0 to n.""")
	
	parser.add_argument("--voltage", type=int,default=300,help="""Default=300. Voltage of the microscope with which the images where collected.""")
	
	(options, args) = parser.parse_args()		
	

	'''
	Make sure input file and paramters are valid/workable
	'''
	options=checkinput(options)
	n=EMUtil.get_image_count()
	check_healthy_input(options,n)
	angles=get_angles(options,n)


	'''
	Log current run of the program ONLY if input parameters are healthy
	'''
	logger = E2init(sys.argv, options.ppid)
	
	for i in range(n):
		angle=angles[i]
		img = EMData(options.input,i)
		nx=img['nx']
		ny=img['ny']
		
		defocusmin = 0.5
		defocusmax = 9.0

		global_defocus = get_global_defocus(options,img,defocusmin,defocusmax)

	if options.defocus:
		defocusmin = options.defocus - 1.5
		defocusmax = options.defocus + 1.5
		#ctfer(options,i,angle,n)	

	E2end(logger)
	return


#def ctfer(options,index,angle,n):
#
#	coords=tile_grid(nx,ny,tilesize,overlap=True,pad=False)
#	ntiles-len(coords)
#	for coords
#
#	return


def get_global_defocus(options,img,defocusmin,defocusmax):
	#grid_coords=tile_grid(nx,ny,options.tilesize):
	tiles = get_tiles(img,tilesize,True)
	tiles_fft_avg = incoherent_sum_from_imglist(tiles,scale=False,checkcorners=True)
	global_defocus = None
	#global_defocus = fit_defocus(tiles_fft_avg, options.voltage, options.cs, options.ampcont, options.apix, options.defocus, defocusmin, defocusmax, 0.1)

	return global_defocus


def get_angles( options ):
	
	angles = {}
	
	if n >1:
		with open( options.tltfile, 'r' ) as f:
			lines = [ int(line.replace('\t','').replace('\n','')) for line in f.readlines()]
			angles = { i : lines[i] for i in range(0, len(lines) ) }

		if angles:
			if options.verbose > 9:
				print("\n(e2tomo_ctfraw.py)(getangles) angles={}".format(angles) )
			return angles
		else:
			print("\nERROR: failed to read angles from --tltfile={}".format(options.tltfile))
			sys.exit(1)
	elif n==1:
		angles.update({0:int(options.tiltangle)})
		return angles

	return


def check_healthy_input(options,n):
	if n==1 and not options.tiltangle:
		print("\nERROR: for single images, --tiltangle is required")
		sys.exit(1)
	elif not options.tltfile:
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

'''
def fit_defocus(fft_avg,voltage,cs,ampcont,apix,target_defocus,defocusmin=0.5,defocusmax=9.0,defocusstep=0.1):
	"""Fit defocus from an incoherent average of tiles 
	Author: Jesus Montoya, jgalaz@gmail.com, 09/2019
	"""
	tilesize = fft_avg['ny']
	
	fftbg = fft_avg.process("math.nonconvex")
	fft1d = fft_avg.calc_radial_dist(old_div(fft_avg.get_ysize(),2),0.0,1.0,1)	# note that this handles the ri2inten averages properly
	
	ctf = EMAN2Ctf()
	try:

		# Compute 1-D curve and background
		ds = old_div(1.0,( apix * tilesize ))
		bg_1d = e2ctf.low_bg_curve(fft1d,ds)
		
		#initial fit, background adjustment, refine fit, final background adjustment
		#ctf = e2ctf.ctf_fit(fft1d,bg_1d,bg_1d,ffta,fftbg,options.voltage,options.cs,options.ampcont,apix,1,dfhint=( defocusmingradient, defocusmaxgradient, step ))

		ctf = e2ctf.ctf_fit( fft1d, bg_1d, bg_1d, ffta, fftbg, options.voltage, options.cs, options.ampcont, apix, 1, dfhint=( defocusmin, defocusmax, defocusstep ) )
		ctf = bgAdj(ctf,fft1d)
		return ctf.defocus
	except:
		print("ctf fit failed! first try")
		print("len fft1d is", len(fft1d))
		print("ffta is", ffta)
		ctf = EMAN2Ctf()
			
		try:
			#ctf = e2ctf.ctf_fit(fft1d,ctf.background,ctf.background,ffta,fftbg,options.voltage,options.cs,options.ampcont,apix,1,dfhint=( defocusmingradient, defocusmaxgradient, step))
			ctf = e2ctf.ctf_fit(fft1d,ctf.background,ctf.background,ffta,fftbg,options.voltage,options.cs,options.ampcont,apix,1,dfhint=( defocusmin, defocusmax, defocusstep ))	
			ctf = bgAdj(ctf,fft1d)

			#if options.astigmatism: 
			#	e2ctf.ctf_fit_stig(ffta,fftbg,ctf)
		except:
			print("ctf fit failed! second adjustment try")
			print("len fft1d is", len(fft1d))
			print("ffta is", ffta)
			return None

	return 
'''

def adjust_ctf_bg(ctf,fg_1d):
	"""Smooths the background based on the values of the foreground near the CTF zeroes and puts the
	smoothed background into the CTF object. From e2ctf.py"""
	
	ds=ctf.dsbg
	#print "\n (bgAdj) ds is", ds
	ctf=ctf
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

if '__main__' == __name__:
	main()

