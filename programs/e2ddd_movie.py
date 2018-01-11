#!/usr/bin/env python
from __future__ import print_function

#
# Author: Steven Ludtke, 02/12/2013 (sludtke@bcm.edu). Updated on 08/28/16.
# Modified by James Michael Bell, 03/27/2017 (jmbell@bcm.edu)
# Copyright (c) 2000-2013 Baylor College of Medicine
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

from EMAN2 import *
from numpy import *
import pprint
import sys
import os
from sys import argv
from time import sleep,time
import threading
import Queue
import numpy as np
from sklearn import linear_model
from scipy import optimize

def main():

	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <ddd_movie_stack>

	This program will do various processing operations on "movies" recorded on direct detection cameras. It
	is primarily used to do whole-frame alignment of movies using all-vs-all CCFs with a global optimization
	strategy. Several outputs including different frame subsets are produced, as well as a text file with the
	translation vector map.

	See e2ddd_particles for per-particle alignment.

	Note: We have found the following to work on DE64 images:
	e2ddd_movie.py <movies> --de64 --dark <dark_frames> --gain <gain_frames> --gain_darkcorrected --reverse_gain --invert_gain

	Note: For multi-image files in MRC format, use the .mrcs extension. Do not use .mrc, as it will handle input stack as a 3D volume.
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_pos_argument(name="movies",help="List the movies to align.", default="", guitype='filebox', browser="EMMovieDataTable(withmodal=True,multiselect=True)",  row=0, col=0,rowspan=1, colspan=3, mode="align")

	parser.add_header(name="orblock1", help='Just a visual separation', title="Dark/Gain Correction", row=2, col=0, rowspan=2, colspan=3, mode="align")

	#parser.add_header(name="orblock2", help='Just a visual separation', title="- CHOOSE FROM -", row=3, col=0, rowspan=1, colspan=3, mode="align")

	parser.add_argument("--dark",type=str,default=None,help="Perform dark image correction using the specified image file",guitype='filebox',browser="EMMovieDataTable(withmodal=True,multiselect=False)", row=4, col=0, rowspan=1, colspan=3, mode="align")
	parser.add_argument("--rotate_dark",  default = "0", type=str, choices=["0","90","180","270"], help="Rotate dark reference by 0, 90, 180, or 270 degrees. Default is 0. Transformation order is rotate then reverse.",guitype='combobox', choicelist='["0","90","180","270"]', row=5, col=0, rowspan=1, colspan=1, mode='align')
	parser.add_argument("--reverse_dark", default=False, help="Flip dark reference along y axis. Default is False. Transformation order is rotate then reverse.",action="store_true",guitype='boolbox', row=5, col=1, rowspan=1, colspan=1, mode='align')

	parser.add_argument("--gain",type=str,default=None,help="Perform gain image correction using the specified image file",guitype='filebox',browser="EMMovieDataTable(withmodal=True,multiselect=False)", row=6, col=0, rowspan=1, colspan=3, mode="align")
	parser.add_argument("--k2", default=False, help="Perform gain image correction on gain images from a Gatan K2. Note, these are the reciprocal of typical DDD gain images.",action="store_true",guitype='boolbox', row=7, col=0, rowspan=1, colspan=1, mode='align')
	parser.add_argument("--rotate_gain", default = 0, type=str, choices=["0","90","180","270"], help="Rotate gain reference by 0, 90, 180, or 270 degrees. Default is 0. Transformation order is rotate then reverse.",guitype='combobox', choicelist='["0","90","180","270"]', row=7, col=1, rowspan=1, colspan=1, mode='align')
	parser.add_argument("--reverse_gain", default=False, help="Flip gain reference along y axis (about x axis). Default is False. Transformation order is rotate then reverse.",action="store_true",guitype='boolbox', row=7, col=2, rowspan=1, colspan=1, mode='align')
	parser.add_argument("--de64", default=False, help="Perform gain image correction on DE64 data. Note, these should not be normalized.",action="store_true",guitype='boolbox', row=8, col=0, rowspan=1, colspan=1, mode='align')
	parser.add_argument("--gain_darkcorrected", default=False, help="Do not dark correct gain image. False by default.",action="store_true",guitype='boolbox', row=8, col=1, rowspan=1, colspan=1, mode='align')
	parser.add_argument("--invert_gain", default=False, help="Use reciprocal of input gain image",action="store_true",guitype='boolbox', row=8, col=2, rowspan=1, colspan=1, mode='align')

	parser.add_argument("--bad_columns", type=str, help="Comma separated list of camera defect columns",default="")
	parser.add_argument("--bad_rows", type=str, help="Comma separated list of camera defect rows",default="")

	#parser.add_header(name="orblock3", help='Just a visual separation', title="- OR -", row=6, col=0, rowspan=1, colspan=3, mode="align")

	parser.add_header(name="orblock4", help='Just a visual separation', title="Output: ", row=10, col=0, rowspan=2, colspan=1, mode="align")

	parser.add_argument("--align_frames", action="store_true",help="Perform whole-frame alignment of the input stacks",default=False, guitype='boolbox', row=11, col=1, rowspan=1, colspan=1, mode='align[True]')

	parser.add_argument("--allali", default=False, help="Average of all aligned frames.",action="store_true", guitype='boolbox', row=12, col=0, rowspan=1, colspan=1, mode='align[True]')
	parser.add_argument("--noali", default=False, help="Average of non-aligned frames.",action="store_true", guitype='boolbox', row=12, col=1, rowspan=1, colspan=1, mode='align')
	parser.add_argument("--rangeali", default="", help="Average frames 'n1-n2'",type=str, guitype='strbox', row=12, col=2, rowspan=1, colspan=1, mode='align')
	parser.add_argument("--goodali", default=False, help="Average of good aligned frames.",action="store_true", guitype='boolbox', row=13, col=0, rowspan=1, colspan=1, mode='align')
	parser.add_argument("--bestali", default=False, help="Average of best aligned frames.",action="store_true", guitype='boolbox', row=13, col=1, rowspan=1, colspan=1, mode='align')
	parser.add_argument("--ali4to14", default=False, help="Average of frames from 4 to 14.",action="store_true",guitype='boolbox', row=13, col=2, rowspan=1, colspan=1, mode='align')

	parser.add_header(name="orblock5", help='Just a visual separation', title="Optional parameters: ", row=14, col=0, rowspan=2, colspan=3, mode="align")

	parser.add_argument("--step",type=str,default="0,1",help="Specify <first>,<step>,[last]. Processes only a subset of the input data. ie- 0,2 would process all even particles. Same step used for all input files. [last] is exclusive. Default= 0,1",guitype='strbox', row=17, col=0, rowspan=1, colspan=1, mode="align")
	parser.add_argument("--fixbadpixels",action="store_true",default=False,help="Tries to identify bad pixels in the dark/gain reference, and fills images in with sane values instead", guitype='boolbox', row=17, col=1, rowspan=1, colspan=1, mode='align[True]')
	parser.add_argument("--normaxes",action="store_true",default=False,help="Tries to erase vertical/horizontal line artifacts in Fourier space by replacing them with the mean of their neighboring values.",guitype='boolbox', row=17, col=2, rowspan=1, colspan=1, mode='align')
	parser.add_argument("--highdose", default=False, help="Use this flag when aligning high dose data (where features in each frame can be distinguished visually).",action="store_true",guitype='boolbox', row=18, col=0, rowspan=1, colspan=1,mode='align')
	parser.add_argument("--phaseplate", default=False, help="Use this flag when aligning phase plate frames.",action="store_true",guitype='boolbox', row=18, col=1, rowspan=1, colspan=1,mode='align')

	parser.add_argument("--frames",action="store_true",default=False,help="Save the dark/gain corrected frames. Note that frames will be overwritten if identical --suffix is already present.", guitype='boolbox', row=19, col=0, rowspan=1, colspan=1, mode='align')
	parser.add_argument("--ext",default="hdf",type=str, choices=["hdf","mrcs","mrc"],help="Save frames with this extension. Default is 'hdf'.", guitype='strbox', row=19, col=1, rowspan=1, colspan=1, mode='align')
	parser.add_argument("--suffix",type=str,default="proc",help="Specify a unique suffix for output frames. Default is 'proc'. Note that the output of --frames will be overwritten if identical suffix is already present.",guitype='strbox', row=19, col=2, rowspan=1, colspan=1, mode="align")

	parser.add_argument("--round", choices=["float","int"],help="If float (default), apply subpixel frame shifts. If integer, use integer shifts.",default="float",guitype='combobox', choicelist='["float","integer"]', row=18, col=2, rowspan=1, colspan=1, mode='align')
	parser.add_argument("--threads", default=4,type=int,help="Number of threads to run in parallel on a single computer when multi-computer parallelism isn't useful", guitype='intbox', row=20, col=0, rowspan=1, colspan=2, mode="align")

	parser.add_header(name="orblock6", help='Just a visual separation', title="Alignment optimization: ", row=22, col=0, rowspan=2, colspan=3, mode="align")

	parser.add_argument("--optbox", type=int,help="Box size to use during alignment optimization. Default is 256.",default=256, guitype='intbox', row=24, col=0, rowspan=1, colspan=1, mode="align")
	parser.add_argument("--optstep", type=int,help="Step size to use during alignment optimization. Default is 224.",default=224,  guitype='intbox', row=24, col=1, rowspan=1, colspan=1, mode="align")
	parser.add_argument("--optalpha", type=float,help="Penalization to apply during robust regression. Default is 0.1. If 0.0, unpenalized least squares will be performed (i.e., no trajectory smoothing).",default=0.1, guitype='floatbox', row=24, col=2, rowspan=1, colspan=1, mode="align")

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--debug", default=False, action="store_true", help="run with debugging output")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)

	(options, args) = parser.parse_args()

	if len(args)<1:
		print(usage)
		parser.error("Specify input DDD stack to be processed.")

	if options.frames == False and options.noali == False and options.align_frames == False:
		print("No outputs specified. See --frames, --noali, or --align_frames. Exiting.") 
		sys.exit(1)

	if options.align_frames == True:
		if options.allali == False and options.rangeali == False and options.goodali == False and options.bestali == False and options.ali4to14 == False:
			print("No post alignment outputs specified. Try with --allali, --rangeali, --goodali, --bestali, or --ali4to14. Exiting.")
			sys.exit(1)

	if options.bad_columns == "": options.bad_columns = []
	else: 
		try: options.bad_columns = [int(c) for c in options.bad_columns.split(",")]
		except:
			print("Error: --bad_columns contains nonnumeric input.")
			sys.exit(1)

	if options.bad_rows == "": options.bad_rows = []
	else: 
		try: options.bad_rows = [int(r) for r in options.bad_rows.split(",")]
		except:
			print("Error: --bad_rows contains nonnumeric input.")
			sys.exit(1)

	# try: os.mkdir("micrographs")
	# except: pass

	pid=E2init(sys.argv)

	if options.dark :
		if options.verbose: print("Loading Dark Reference")
		if "e2ddd_darkref" in options.dark:
			dark = EMData(options.dark,-1)
		else:
			if options.dark[-4:].lower() in (".mrc") :
				dark_hdr = EMData(options.dark,0,True)
				nx = dark_hdr["nx"]
				ny = dark_hdr["ny"]
				nd = dark_hdr["nz"]
				dark=EMData(options.dark,0,False,Region(0,0,0,nx,ny,1))
			else:
				nd=EMUtil.get_image_count(options.dark)
				dark = EMData(options.dark,0)
			if nd>1:
				sigd=dark.copy()
				sigd.to_zero()
				a=Averagers.get("mean",{"sigma":sigd,"ignore0":1})
				print("Summing dark")
				for i in xrange(0,nd):
					if options.verbose:
						sys.stdout.write("({}/{})   \r".format(i+1,nd))
						sys.stdout.flush()
					if options.dark[-4:].lower() in (".mrc") :
						t=EMData(options.dark,0,False,Region(0,0,i,nx,ny,1))
					else:
						t=EMData(options.dark,i)
					t.process_inplace("threshold.clampminmax",{"minval":0,"maxval":t["mean"]+t["sigma"]*3.5,"tozero":1})
					a.add_image(t)
				dark=a.finish()
				sigd.write_image(options.dark.rsplit(".",1)[0]+"_sig.hdf")
				if options.fixbadpixels:
					sigd.process_inplace("threshold.binary",{"value":sigd["sigma"]/10.0}) # Theoretically a "perfect" pixel would have zero sigma, but in reality, the opposite is true
					dark.mult(sigd)
				dark.write_image(options.dark.rsplit(".",1)[0]+"_sum.hdf")
			#else: dark.mult(1.0/99.0)
			dark.process_inplace("threshold.clampminmax.nsigma",{"nsigma":3.0})
			dark2=dark.process("normalize.unitlen")
	else : dark=None

	if options.gain :
		if options.verbose: print("Loading Gain Reference")
		if "e2ddd_gainref" in options.gain:
			gain = EMData(options.gain,-1)
		else:
			if options.k2: gain=EMData(options.gain)
			else:
				if options.gain[-4:].lower() in (".mrc") :
					gain_hdr = EMData(options.gain,0,True)
					nx = gain_hdr["nx"]
					ny = gain_hdr["ny"]
					nd = gain_hdr["nz"]
					gain=EMData(options.gain,0,False,Region(0,0,0,nx,ny,1))
				else:
					nd=EMUtil.get_image_count(options.gain)
					gain = EMData(options.gain,0)
				if nd>1:
					sigg=gain.copy()
					sigg.to_zero()
					a=Averagers.get("mean",{"sigma":sigg,"ignore0":1})
					print("Summing gain")
					for i in xrange(0,nd):
						if options.verbose:
							sys.stdout.write("({}/{})   \r".format(i+1,nd))
							sys.stdout.flush()
						if options.dark[-4:].lower() in (".mrc") :
							t=EMData(options.gain,0,False,Region(0,0,i,nx,ny,1))
						else:
							t=EMData(options.gain,i)
						#t.process_inplace("threshold.clampminmax.nsigma",{"nsigma":4.0,"tozero":1})
						t.process_inplace("threshold.clampminmax",{"minval":0,"maxval":t["mean"]+t["sigma"]*3.5,"tozero":1})
						a.add_image(t)
					gain=a.finish()
					sigg.write_image(options.gain.rsplit(".",1)[0]+"_sig.hdf")
					if options.fixbadpixels:
						sigg.process_inplace("threshold.binary",{"value":sigg["sigma"]/10.0}) # Theoretically a "perfect" pixel would have zero sigma, but in reality, the opposite is true
						if dark!=None : sigg.mult(sigd)
						gain.mult(sigg)
					gain.write_image(options.gain.rsplit(".",1)[0]+"_sum.hdf")
				if options.de64:
					gain.process_inplace( "threshold.clampminmax", { "minval" : gain[ 'mean' ] - 8.0 * gain[ 'sigma' ], "maxval" : gain[ 'mean' ] + 8.0 * gain[ 'sigma' ], "tomean" : True } )
				else:
					gain.process_inplace("math.reciprocal",{"zero_to":0.0})
					#gain.mult(1.0/99.0)
					#gain.process_inplace("threshold.clampminmax.nsigma",{"nsigma":3.0})

			if dark!=None and options.gain_darkcorrected == False: gain.sub(dark) # dark correct the gain-reference

			if options.de64:
				mean_val = gain["mean"]
				if mean_val <= 0.: mean_val=1.
				gain.process_inplace("threshold.belowtominval",{"minval":0.01,"newval":mean_val})

			gain.mult(1.0/gain["mean"])

			if options.invert_gain: gain.process_inplace("math.reciprocal")
	#elif options.gaink2 :
	#	gain=EMData(options.gaink2)
	else : gain=None

	if options.rotate_gain and gain != None:
		tf = Transform({"type":"2d","alpha":int(options.rotate_gain)})
		gain.process_inplace("xform",{"transform":tf})

	if options.reverse_gain: gain.process_inplace("xform.reverse",{"axis":"y"})

	if options.rotate_dark and dark != None:
		tf = Transform({"type":"2d","alpha":int(options.rotate_dark)})
		dark.process_inplace("xform",{"transform":tf})

	if options.reverse_dark: dark.process_inplace("xform.reverse",{"axis":"y"})

	if gain or dark:
		try: os.mkdir("movies")
		except: pass

	if gain:
		gainname="movies/e2ddd_gainref.hdf"
		gain.write_image(gainname,-1)
		gainid=EMUtil.get_image_count(gainname)-1
		gain["filename"]=gainname
		gain["fileid"]=gainid

	if dark:
		darkname="movies/e2ddd_darkref.hdf"
		dark.write_image(darkname,-1)
		darkid=EMUtil.get_image_count(darkname)-1
		dark["filename"]=darkname
		dark["fileid"]=darkid

	step = options.step.split(",")

	if len(step) == 3 : last = int(step[2])
	else :	last = -1

	first = int(step[0])
	step  = int(step[1])

	if options.verbose : print("Range = {} - {}, Step = {}".format(first, last, step))

	# the user may provide multiple movies to process at once
	for fsp in args:
		if options.verbose : print("Processing {}".format(fsp))

		n = EMUtil.get_image_count(fsp)

		if n < 3 :
			hdr = EMData(fsp, 0, True)

			if hdr["nz"] < 2 :
				print("ERROR: {} has only {} images. Min 3 required.".format(fsp, n))
				continue

			n = hdr["nz"]

		if last <= 0 : flast = n
		else : flast = last

		if flast > n : flast = n

		process_movie(fsp, dark, gain, first, flast, step, options)

	E2end(pid)

def process_movie(fsp,dark,gain,first,flast,step,options):

		if fsp[-4:].lower() in (".mrc"):
			hdr=EMData(fsp,0,True)			# read header
			nx,ny=hdr["nx"],hdr["ny"]

		# bgsub and gain correct the stack
		outim=[]
		nfs = 0
		t = time()
		for ii in xrange(first,flast,step):
			if options.verbose:
				sys.stdout.write(" {}/{}   \r".format(ii-first+1,flast-first+1))
				sys.stdout.flush()

			#if fsp[-4:].lower() in (".mrc","mrcs") :
			if fsp[-4:].lower() in (".mrc") :
				im=EMData(fsp,0,False,Region(0,0,ii,nx,ny,1))
			else: im=EMData(fsp,ii)

			if dark!=None : im.sub(dark)
			if gain!=None : im.mult(gain)
			#im.process_inplace("threshold.clampminmax",{"minval":0,"maxval":im["mean"]+im["sigma"]*3.5,"tozero":1})
			if options.de64: im.process_inplace( "threshold.clampminmax", { "minval" : im[ 'minimum' ], "maxval" : im[ 'mean' ] + 8.0 * im[ 'sigma' ], "tomean" : True } )
			if options.fixbadpixels : im.process_inplace("threshold.outlier.localmean",{"sigma":3.5,"fix_zero":1}) # fixes clear outliers as well as values which were exactly zero

			#im.process_inplace("threshold.clampminmax.nsigma",{"nsigma":3.0})
#			im.mult(-1.0)
			#if options.normalize: im.process_inplace("normalize.edgemean") 

			if options.bad_rows != [] or options.bad_columns != []:
				im = im.process("math.xybadlines",{"rows":options.bad_rows,"cols":options.bad_columns})

			#if options.frames : im.write_image(outname[:-4]+"_corr.hdf",ii-first)
			if options.frames:
				if options.ext == "mrc": outname="{}_{}.{}".format(fsp.rsplit(".",1)[0],options.suffix,"mrcs") #Output contents vary with options
				else: outname="{}_{}.{}".format(fsp.rsplit(".",1)[0],options.suffix,options.ext) #Output contents vary with options
				im.write_image(outname,ii-first)
			outim.append(im)

		if options.ext == "mrc": os.rename(outname,outname.replace(".mrcs",".mrc"))

		if options.noali:
			mgdirname = "micrographs_noali"
			try: os.mkdir(mgdirname)
			except: pass
			alioutname="{}/{}.hdf".format(mgdirname,base_name(fsp))
			out=qsum(outim)
			#write out the unaligned average movie
			out.write_image(alioutname,0)

		t1 = time()-t
		print("{:.1f} s".format(time()-t))

		if options.align_frames :

			start = time()

			nx=outim[0]["nx"]
			ny=outim[0]["ny"]

			md = min(nx,ny)
			if md <= 2048:
				print("This program does not facilitate alignment of movies with frames smaller than 2048x2048 pixels.")
				sys.exit(1)

			n=len(outim)
			nx=outim[0]["nx"]
			ny=outim[0]["ny"]

			print("{} frames read {} x {}".format(n,nx,ny))

			ccfs=Queue.Queue(0)

			# prepare image data (outim) by clipping and FFT'ing all tiles (this is threaded as well)
			immx=[0]*n
			thds = []
			for i in range(n):
				thd = threading.Thread(target=split_fft,args=(options,outim[i],i,options.optbox,options.optstep,ccfs))
				thds.append(thd)
			sys.stdout.write("\rPrecompute  /{} FFTs".format(len(thds)))
			t0=time()

			thrtolaunch=0
			while thrtolaunch<len(thds) or threading.active_count()>1:
				if thrtolaunch<len(thds) :
					while (threading.active_count()==options.threads ) : sleep(.01)
					#if options.verbose :
					#	sys.stdout.write("\rPrecompute {}/{} FFTs {}".format(thrtolaunch+1,len(thds),threading.active_count()))
					#	sys.stdout.flush()
					thds[thrtolaunch].start()
					thrtolaunch+=1
				else: sleep(0.5)

				while not ccfs.empty():
					i,d=ccfs.get()
					immx[i]=d

			for th in thds: th.join()
			print()

			# create threads
			thds=[]
			peak_locs=Queue.Queue(0)
			i=-1
			for ima in range(n-1):
				for imb in range(ima+1,n):
					if options.verbose>3: i+=1		# if i>0 then it will write pre-processed CCF images to disk for debugging
					thds.append(threading.Thread(target=calc_ccf_wrapper,args=(options,(ima,imb),options.optbox,options.optstep,immx[ima],immx[imb],ccfs,peak_locs,i,fsp)))

			print("{:1.1f} s\nCompute {} ccfs".format(time()-t0,len(thds)))
			t0=time()

			# here we run the threads and save the results, no actual alignment done here
			csum2={}

			thrtolaunch=0
			while thrtolaunch<len(thds) or threading.active_count()>1:
				# If we haven't launched all threads yet, then we wait for an empty slot, and launch another
				# note that it's ok that we wait here forever, since there can't be new results if an existing
				# thread hasn't finished.
				if thrtolaunch<len(thds) :
					while (threading.active_count()==options.threads ) : sleep(.01)
					#if options.verbose : print "Starting thread {}/{}".format(thrtolaunch,len(thds))
					thds[thrtolaunch].start()
					thrtolaunch+=1
				else:
					sleep(0.5)

				while not ccfs.empty():
					i,d=ccfs.get()
					csum2[i]=d

				if options.verbose:
					sys.stdout.write("\r  {}/{} ({})".format(thrtolaunch,len(thds),threading.active_count()))
					sys.stdout.flush()

			for th in thds: th.join()
			print()

			avgr=Averagers.get("minmax",{"max":0})
			avgr.add_image_list(csum2.values())
			csum=avgr.finish()

			#####
			# Alignment code
			#####

			# array of x,y locations of each frame, all relative to the last frame in the series, which will always have 0,0 shift
			locs=[0]*(n*2) # we store the value for the last frame as well as a conveience

			print("{:1.1f} s\nAlign {} frames".format(time()-t0,n))
			t0=time()

			peak_locs = {p[0]:p[1] for p in peak_locs.queue}

			if options.debug and options.verbose == 9:
				print("PEAK LOCATIONS:")
				for l in peak_locs.keys():
					print(peak_locs[l])

			# if options.ccweight:
			# 	# normalize ccpeak values
			# 	vals = []
			# 	for ima,(i,j) in enumerate(sorted(peak_locs.keys())):
			# 		for imb in range(i,j):
			# 			try:
			# 				vals.append(peak_locs[(i,j)][-1])
			# 			except:
			# 				pass
			# 	ccmean = np.mean(vals)
			# 	ccstd = np.std(vals)

			m = n*(n-1)/2
			bx = np.ones(m)
			by = np.ones(m)
			A = np.zeros([m,n]) # coefficient matrix
			for ima,(i,j) in enumerate(sorted(peak_locs.keys())):
				for imb in range(i,j):
					try:
						bx[ima] = peak_locs[(i,j)][0]
						by[ima] = peak_locs[(i,j)][1]
						A[ima,imb] = 1
						#A[ima,imb] = float(n-np.fabs(i-j))/n
						#A[ima,imb] = np.exp(1-peak_locs[(i,j)][3])
						#A[ima,imb] = sqrt(float(n-fabs(i-j))/n)
					except:
						pass # CCF peak was not found
			b = np.c_[bx,by]
			A = np.asmatrix(A)
			b = np.asmatrix(b)

			# remove all zero rows from A and corresponding entries in b
			z = np.argwhere(np.all(A==0,axis=1))
			A = np.delete(A,z,axis=0)
			b = np.delete(b,z,axis=0)

			regr = linear_model.Ridge(alpha=options.optalpha,normalize=True,fit_intercept=True)
			regr.fit(A,b)

			traj = regr.predict(np.tri(n))
			#shifts = regr.predict(np.eye(n))-options.optbox/2

			traj -= traj[0]

			if options.round == "int": traj = np.round(traj,0)#.astype(np.int8)

			locs = traj.ravel()
			quals=[0]*n # quality of each frame based on its correlation peak summed over all images
			cen=options.optbox/2 #csum2[(0,1)]["nx"]/2
			for i in xrange(n-1):
				for j in xrange(i+1,n):
					val=csum2[(i,j)].sget_value_at_interp(int(cen+locs[j*2]-locs[i*2]),int(cen+locs[j*2+1]-locs[i*2+1]))*sqrt(float(n-fabs(i-j))/n)
					quals[i]+=val
					quals[j]+=val

			print("{:1.1f} s".format(time()-t0,n))

			runtime = time()-start
			print("Runtime: {:.1f} s".format(runtime))

			# if options.plot:
			# 	import matplotlib.pyplot as plt
			# 	fig,ax = plt.subplots(1,3,figsize=(12,3))
			# 	ax[0].plot(traj[:,0],traj[:,1],c='b',alpha=0.5)
			# 	ax[0].scatter(traj[:,0],traj[:,1],c='b',alpha=0.5)
			# 	ax[0].set_title("Trajectory (x/y pixels)")
			# 	ax[1].set_title("Quality (cumulative pairwise ccf value)")
			# 	ax[1].plot(quals,'k')
			# 	for k in peak_locs.keys():
			# 		try:
			# 			p = peak_locs[k]
			# 			ax[2].scatter(p[0],p[1])
			# 		except: pass
			# 	ax[2].set_title("CCF Peak Coordinates")

			print("{:1.1f} s\nShift images".format(time()-t0))
			for i,im in enumerate(outim):
				if options.round == "int":
					dx = int(round(locs[i*2],0))
					dy = int(round(locs[i*2+1],0))
					im.translate(dx,dy,0)
				else: # float by default
					dx = float(locs[i*2])
					dy = float(locs[i*2+1])
					im.translate(dx,dy,0)

			if options.normaxes:
				for f in outim:
					f.process_inplace("filter.xyaxes0",{"neighbor":1})

			if options.allali:
				mgdirname = "micrographs_allali"
				try: os.mkdir(mgdirname)
				except: pass
				alioutname="{}/{}.hdf".format(mgdirname,base_name(fsp))
				out=qsum(outim)
				out.write_image(alioutname,0)

			# write translations and qualities
			db=js_open_dict(info_name(fsp))
			db["movieali_trans"]=locs
			db["movieali_qual"]=quals
			db["movie_name"]=fsp
			if gain:
				db["gain_name"]=gain["filename"]
				db["gain_id"]=gain["fileid"]
			if dark:
				db["dark_name"]=dark["filename"]
				db["dark_id"]=dark["fileid"]
			db["runtime"]=runtime
			db["precision"]=options.round
			db["optbox"]=options.optbox
			db["optstep"]=options.optstep
			db["optalpha"]=options.optalpha
			db.close()

			out=open("{}_info.txt".format(fsp[:-4]),"w")
			out.write("#i,dx,dy,dr,rel dr,qual\n")
			for i in range(1,n):
				dx,dy = traj[i]
				dxlast,dylast = traj[i-1]
				dr = hypot(dx,dy)
				reldr = hypot(dx-dxlast,dy-dylast)
				out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(i,dx,dy,dr,reldr,quals[i]))

			if options.goodali:
				mgdirname = "micrographs_goodali"
				try: os.mkdir(mgdirname)
				except: pass
				alioutname="{}/{}.hdf".format(mgdirname,base_name(fsp))
				thr=(max(quals[1:])-min(quals))*0.4+min(quals)	# max correlation cutoff for inclusion
				best=[im for i,im in enumerate(outim) if quals[i]>thr]
				out=qsum(best)
				print("Keeping {}/{} frames".format(len(best),len(outim)))
				out.write_image(alioutname,0)

			if options.bestali:
				mgdirname = "micrographs_bestali"
				try: os.mkdir(mgdirname)
				except: pass
				alioutname="{}/{}.hdf".format(mgdirname,base_name(fsp))
				thr=(max(quals[1:])-min(quals))*0.6+min(quals)	# max correlation cutoff for inclusion
				best=[im for i,im in enumerate(outim) if quals[i]>thr]
				out=qsum(best)
				print("Keeping {}/{} frames".format(len(best),len(outim)))
				out.write_image(alioutname,0)

			if options.ali4to14:
				mgdirname = "micrographs_4-14"
				try: os.mkdir(mgdirname)
				except: pass
				alioutname="{}/{}.hdf".format(mgdirname,base_name(fsp))
				# skip the first 4 frames then keep 10
				out=qsum(outim[4:14])
				out.write_image(alioutname,0)

			if len(options.rangeali)>0:
				try: rng=[int(i) for i in options.rangeali.split("-")]
				except:
					print("Error: please specify --rangeali as X-Y where X and Y are inclusive starting with 0")
					sys.exit(1)
				mgdirname = "micrographs_{}-{}".format(rng[0],rng[1])
				try: os.mkdir(mgdirname)
				except: pass
				alioutname="{}/{}.hdf".format(mgdirname,base_name(fsp))
				out=qsum(outim[rng[0]:rng[1]+1])
				out.write_image(alioutname,0)

			print("Done")

# CCF calculation
def calc_ccf_wrapper(options,N,box,step,dataa,datab,out,locs,ii,fsp):

	for i in range(len(dataa)):
		c=dataa[i].calc_ccf(datab[i],fp_flag.CIRCULANT,True)
		try: csum.add(c)
		except: csum=c

	if options.debug:
		if fsp[-5:] == ".mrcs":
			ggg = "{}-ccf_imgs.hdf".format(fsp.replace(".mrcs",""))
		elif fsp[-4:] == ".hdf":
			ggg = "{}-ccf_imgs.hdf".format(fsp.replace(".hdf",""))
		elif fsp[:-4] == ".mrc":
			ggg = "{}-ccf_imgs.hdf".format(fsp.replace(".mrc",""))
		elif fsp[-4:] == ".tif":
			ggg = "{}-ccf_imgs.hdf".format(fsp.replace(".tif",""))
		csum.process("normalize.edgemean").write_image(ggg,ii)

	xx = np.linspace(0,box,box)
	yy = np.linspace(0,box,box)
	xx,yy = np.meshgrid(xx,yy)

	popt,ccpeakval = bimodal_peak_model(options,csum)
	if popt == None:
		csum = from_numpy(np.zeros((box,box))).process("normalize.edgemean")
		locs.put((N,[]))
		out.put((N,csum))
	else:
		#if ii>=0: csum.process("normalize.edgemean").write_image("ccf_models.hdf",ii)
		if options.phaseplate:
			ncc = csum.numpy().copy()
			pcsum = neighbormean_origin(ncc)
			csum = from_numpy(pcsum)
			locs.put((N,[popt[0],popt[1],ccpeakval,csum["maximum"]]))
		else:
			cc_model = correlation_peak_model((xx,yy),popt[0],popt[1],popt[2],popt[3]).reshape(box,box)
			csum = from_numpy(cc_model)
			locs.put((N,[popt[0],popt[1],popt[2],popt[3],popt[4],popt[5],ccpeakval,csum["maximum"]]))
		out.put((N,csum))
	if ii>=0 and options.debug:
		if fsp[-5:] == ".mrcs":
			fff = "{}-ccf_models.hdf".format(fsp.replace(".mrcs",""))
		elif fsp[-4:] == ".hdf":
			fff = "{}-ccf_models.hdf".format(fsp.replace(".hdf",""))
		elif fsp[:-4] == ".mrc":
			fff = "{}-ccf_models.hdf".format(fsp.replace(".mrc",""))
		elif fsp[-4:] == ".tif":
			fff = "{}-ccf_models.hdf".format(fsp.replace(".tif",""))
		csum.process("normalize.edgemean").write_image(fff,ii)

# preprocess regions by normalizing and doing FFT
def split_fft(options,img,i,box,step,out):
	lst=[]
	# if min(img["nx"],img["ny"]) > 6000:
	# 	img.process_inplace("math.fft.resample",{"n":1.5})
	nx = img["nx"]
	ny = img["ny"]
	#proc = img.process("filter.highpass.gauss",{"cutoff_pixels":2})
	# if options.binning == -1:
	# 	if img["nx"] < 5000: pass
	# 	else:
	# 		#if options.debug: print("FFT resample by 6")
	# 		img.process_inplace("math.fft.resample",{"n":4})
	# elif options.binning >= 0:
	# 	#if options.debug: print("FFT resample by {}".format(options.binning))
	# 	img.process_inplace("math.fft.resample",{"n":options.binning})
	#img.process_inplace("filter.highpass.gauss",{"cutoff_pixels":2})
	#img.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.4})
	#proc.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.3})
	#patchid = 0
	for dx in range(box/2,nx-box,step):
		for dy in range(box/2,ny-box,step):
			clp = img.get_clip(Region(dx,dy,box,box))
			#if options.normalize: clp.process_inplace("normalize.edgemean")
			#if box >= 512:
			#	if not options.tomo:
			#		if img["apix_x"] == 1.0: # likely an image with incorrect header or simulated data
			#			clp.process_inplace("filter.highpass.gauss",{"cutoff_pixels":2})
			#		else: 
			#	clp.process_inplace("filter.highpass.gauss",{"cutoff_abs":0.01})
				#clp.process_inplace("math.fft.resample",{"n":1.})
				#if options.tomo:
				#	clp.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.3})
			clp.process_inplace("normalize.edgemean")
			#if options.debug:
			#	clp["patch_id"] = patchid
			#	clp["frame_id"] = i
			#	clp.write_image("patch_{}.hdf".format(str(patchid).zfill(4)),i)
			#patchid += 1
			lst.append(clp.do_fft())
	out.put((i,lst))

def correlation_peak_model(x_y, xo, yo, sigma, amp):
	x, y = x_y
	if sigma <= 0: return np.zeros_like(x)
	xo = float(xo)
	yo = float(yo)
	g = amp*np.exp(-(((x-xo)**2)+((y-yo)**2))/(2.*sigma**2))
	return g.ravel()

def fixedbg_peak_model(x_y, sigma, amp):
	x, y = x_y
	if sigma <= 0: return np.zeros_like(x)
	xo = float(len(x)/2)
	yo = float(len(y)/2)
	g = amp*np.exp(-(((x-xo)**2)+((y-yo)**2))/(2.*sigma**2))
	return g.ravel()

def twod_bimodal(x_y,x1,y1,sig1,amp1,sig2,amp2):
	x, y = x_y
	if sig1 <= 0 or sig2 <= 0: return np.zeros_like(x)
	#correlation_peak = correlation_peak_model((x,y),x1,y1,sig1,amp1)
	cp = amp1*np.exp(-(((x-x1)**2+(y-y1)**2))/(2.*sig1**2))
	#fixedbg_peak = fixedbg_peak_model((x,y),sig2,amp2)
	fp = amp2*np.exp(-(((x-len(x)/2.)**2+(y-len(y)/2.)**2))/(2.*sig2**2))
	return cp.ravel() + fp.ravel() # + noise

# def edgemean(a,xc,yc):
# 	# mean of values around perimeter of image
# 	return (np.sum(a[:,0])+np.sum(a[0,1:])+np.sum(a[1:,-1])+np.sum(a[-1,1:-1]))/(4*(a.shape[0]-1))

def neighbormean_origin(a): # replace origin pixel with mean of surrounding pixels
	proc = a.copy()
	if a.shape[0] % 2 == 0:
		ac = proc.shape[0]/2
		r = proc[ac-2:ac+2,ac-2:ac+2].copy()
		rc = r.shape[0]/2
		r[rc,rc] = np.nan
		r[rc,rc-1] = np.nan
		r[rc-1,rc] = np.nan
		r[rc-1,rc-1] = np.nan
		nm = np.nanmean(r)
		proc[ac,ac] = nm
		proc[ac,ac-1] = nm
		proc[ac-1,ac] = nm
		proc[ac-1,ac-1] = nm
	else:
		ac = proc.shape[0]/2
		r = proc[ac-2:ac+1,ac-2:ac+1].copy()
		plt.imshow(r)
		rc = r.shape[0]/2
		r[rc,rc] = np.nan
		proc[ac,ac] = np.nanmean(r)
	return proc

def find_com(ccf): # faster alternative to gaussian fitting...less robust in theory.
	thresh = (np.ones(ccf.shape) * np.mean(ccf))+2.5*np.std(ccf)
	m = np.greater(ccf,thresh) * 1.0
	m = m / np.sum(m)
	# marginal distributions
	dx = np.sum(m, 1)
	dy = np.sum(m, 0)
	# expected values
	cx = np.sum(dx * np.arange(ccf.shape[0]))
	cy = np.sum(dy * np.arange(ccf.shape[1]))
	return cx-ccf.shape[0]/2,cy-ccf.shape[1]/2

def bimodal_peak_model(options,ccf):
	nxx = ccf["nx"]
	bs = int(nxx/4)

	xx = np.linspace(0,bs,bs)
	yy = np.linspace(0,bs,bs)
	xx,yy = np.meshgrid(xx,yy)

	r = Region(nxx/2-bs/2,nxx/2-bs/2,bs,bs)

	ccfreg = ccf.get_clip(r)
	ncc = ccfreg.numpy().copy()

	x1 = int(bs/2.)
	y1 = int(bs/2.)
	a1 = 1000.0
	s1 = 10.0
	a2 = 20000.0
	s2 = 0.6

	if options.phaseplate:
		initial_guess = [x1,y1,s1,a1]
		# drop the origin (replace with perimeter mean)
		#mval = edgemean(ncc)
		#ix,iy = ncc.shape
		pncc = neighbormean_origin(ncc)
		x1,y1 = find_com(pncc)
		x1 = x1+nxx/2
		y1 = y1+nxx/2
		# if options.debug:
		# 	try:
		# 		ii = EMUtil.get_image_count("tmp.hdf")
		# 		from_numpy(pncc).write_image("tmp.hdf",ii)`
		# 	except:
		# 		from_numpy(pncc).write_image("tmp.hdf",0)
		return [x1,y1],ccf.sget_value_at_interp(x1,y1)
		# try: # run optimization
		# 	bds = [(-np.inf, -np.inf,  0.01, 0.01),(np.inf, np.inf, 100.0, 100000.0)]
		# 	#bds = [(-bs/2, -bs/2,  0.0, 0.0),(bs/2, bs/2, 100.0, 100000.0)]
		# 	popt,pcov=optimize.curve_fit(correlation_peak_model,(xx,yy),ncc.ravel(),p0=initial_guess,bounds=bds,method='dogbox',max_nfev=50)#,xtol=0.1)#,ftol=0.0001,gtol=0.0001)
		# except:
		# 	return None, -1
	elif options.highdose: # only useful for extremely high contrast frames
		yc,xc = np.where(ncc==ncc.max())
		popt = [float(xc[0]+nxx/2),float(yc[0]+nxx/2),ncc.max(),1.,0.,0.]
		return popt,ccf.sget_value_at_interp(popt[0],popt[1])
	else:
		initial_guess = [x1,y1,s1,a1,s2,a2]
		bds = [(-bs/2, -bs/2, 0.01, 0.01, 0.6, 0.01),(bs/2, bs/2, 100.0, 20000.0, 2.5, 100000.0)]
		try:
			popt,pcov=optimize.curve_fit(twod_bimodal,(xx,yy),ncc.ravel(),p0=initial_guess,bounds=bds,method="dogbox",max_nfev=100,xtol=1e-6,ftol=1e-8,loss='linear')
		except RuntimeError:
			return None,-1 

		popt = [p for p in popt]
		popt[0] = popt[0] + nxx/2 - bs/2
		popt[1] = popt[1] + nxx/2 - bs/2
		popt[2] = np.abs(popt[2])
		return popt,ccf.sget_value_at_interp(popt[0],popt[1])

def qsum(imlist):
	avg=Averagers.get("mean")
	avg.add_image_list(imlist)
	return avg.finish()

def qual(locs,ccfs):
	"""computes the quality of the current alignment. Passed a dictionary of CCF images keyed by (i,j) tuple and
	an (x0,y0,x1,y1,...)  shift array. Smaller numbers are better since that's what the simplex does"""
	nrg=0.0
	cen=ccfs[(0,1)]["nx"]/2
	n=len(locs)/2
	for i in xrange(n-1):
		for j in xrange(i+1,n):
			penalty = sqrt(float(n-fabs(i-j))/n)**2 # This is a recognition that we will tend to get better correlation with near neighbors in the sequence
			locx = int(cen+locs[j*2]-locs[i*2])
			locy = int(cen+locs[j*2+1]-locs[i*2+1])
			nrg-=ccfs[(i,j)].sget_value_at_interp(locx,locy)*penalty
	return nrg

if __name__ == "__main__":
	main()


	#parser.add_argument("--normalize",action="store_true",default=False,help="Apply edgenormalization to input images after dark/gain. Do not use this option when aligning frames with MotionCor2.", guitype='boolbox', row=13, col=0, rowspan=1, colspan=1, mode='align')
	#parser.add_argument("--gaink2",type=str,default=None,help="Perform gain image correction. Gatan K2 gain images are the reciprocal of DDD gain images.",guitype='filebox',browser="EMMovieDataTable(withmodal=True,multiselect=False)", row=7, col=0, rowspan=1, colspan=3, mode="align")
	#parser.add_argument("--plot", default=False,help="Display a plot of the movie trajectory after alignment",action="store_true")
	#parser.add_argument("--simpleavg", action="store_true",help="Will save a simple average of the dark/gain corrected frames (no alignment or weighting)",default=False)
	#parser.add_argument("--avgs", action="store_true",help="Testing",default=False)
	#parser.add_argument("--movie", type=int,help="Display an n-frame averaged 'movie' of the stack, specify number of frames to average",default=0)
	#parser.add_argument("--ccweight", action="store_true",help="Supply coefficient matrix with cross correlation peak values rather than 1s.",default=False)
	#parser.add_argument("--optfsc", default=False, help="Specify whether to compute FSC during alignment optimization. Default is False.",action="store_true")
	#parser.add_argument("--falcon", default=False, help="Use this flag to optimize alignment for falcon detector data.",action="store_true",guitype='boolbox', row=5, col=1, rowspan=1, colspan=1)
	#parser.add_argument("--binning", type=int,help="Bin images by this factor by resampling in Fourier space. Default (-1) will choose based on input box size.",default=-1)
#	parser.add_argument("--align_frames_tree", action="store_true",help="Perform whole-frame alignment of the stack hierarchically",default=False)
#	parser.add_argument("--align_frames_countmode", action="store_true",help="Perform whole-frame alignment of frames collected in counting mode",default=False)
	#parser.add_argument("--save_aligned", action="store_true",help="Save dark/gain corrected and optionally aligned stack",default=False, guitype='boolbox', row=14, col=0, rowspan=1, colspan=1, mode='align[True]')



# def calc_incoherent_pws(frames,bs=2048):
# 	mx = np.arange(bs+50,frames[0]['nx']-bs+50,bs)
# 	my = np.arange(bs+50,frames[1]['ny']-bs+50,bs)
# 	regions = {}
# 	for i in xrange(len(frames)): regions[i] = [[x,y] for y in my for x in mx]
# 	ips = Averagers.get('mean')
# 	for i in xrange(len(frames)):
# 		img = frames[i].copy()
# 		frame_avg = Averagers.get('mean')
# 		for r in regions[i]:
# 			reg = frames[i].get_clip(Region(r[0],r[1],bs,bs))
# 			reg.process_inplace("normalize.unitlen")
# 			reg.do_fft_inplace()
# 			reg.ri2inten()
# 			frame_avg.add_image(reg)
# 		ips.add_image(frame_avg.finish())
# 	ips = ips.finish()
# 	ips.process_inplace("math.sqrt")
# 	ips.process_inplace('normalize.edgemean')
# 	ips_ra = ips.process('math.rotationalaverage')
# 	#ips = ips-ips_ra
# 	#ips.process_inplace("math.rotationalaverage")
# 	return ips

# def calc_coherent_pws(frames,bs=2048):
# 	mx = np.arange(bs+50,frames[0]['nx']-bs+50,bs)
# 	my = np.arange(bs+50,frames[1]['ny']-bs+50,bs)
# 	regions = {}
# 	for i in xrange(len(frames)):
# 		regions[i] = [[x,y] for y in my for x in mx]
# 	stacks = {}
# 	for ir in xrange(len(regions[0])):
# 		stacks[ir] = [regions[i][ir] for i in xrange(len(frames))]
# 	cps = Averagers.get('mean')
# 	for s in xrange(len(stacks)):
# 		stack_avg = Averagers.get('mean')
# 		for i,r in enumerate(stacks[s]):
# 			stack_avg.add_image(frames[i].copy().get_clip(Region(r[0],r[1],bs,bs)))
# 		avg = stack_avg.finish()
# 		avg.process_inplace('normalize.unitlen')
# 		avg.do_fft_inplace()
# 		avg.ri2inten()
# 		cps.add_image(avg)
# 	cps = cps.finish()
# 	cps.process_inplace('math.sqrt')
# 	cps.process_inplace('normalize.edgemean')
# 	return cps


		# A simple average

		# if options.simpleavg :
		# 	if options.verbose : print("Simple average")
		# 	avgr=Averagers.get("mean")
		# 	for i in xrange(len(outim)):						# only use the first second for the unweighted average
		# 		if options.verbose:
		# 			sys.stdout.write(" {}/{}   \r".format(i+1,len(outim)))
		# 			sys.stdout.flush()
		# 		avgr.add_image(outim[i])
		# 	print("")

		# 	av=avgr.finish()
		# 	if first!=1 or flast!=-1 : av.write_image(outname[:-4]+"_{}-{}_mean.hdf".format(first,flast),0)
		# 	else: av.write_image(outname[:-4]+"_mean.hdf",0)


			#out=sum(outim[5:15])	# FSC with the earlier frames instead of whole average
			# compute fsc between each aligned frame and the average
			# we tile this for better curves, since we don't need the detail
			#fscq=[0]*n
			# if options.optfsc:
			# 	for i in range(n):
			# 		rgnc=0
			# 		for x in range(64,out["nx"]-192,64):
			# 			for y in range(64,out["ny"]-192,64):
			# 				rgnc+=1.0
			# 				cmpto=out.get_clip(Region(x,y,64,64))
			# 				cscen=outim[i].get_clip(Region(x,y,64,64))
			# 				s,f=calcfsc(cmpto,cscen)
			# 				f=array(f)
			# 				try: fs+=f
			# 				except: fs=f
			# 		fs/=rgnc
			# 		fs=list(fs)
			# 		fscq.append(qsum(fs[2:24]))

			# 		Util.save_data(s[1],s[1]-s[0],fs[1:-1],"{}_fsc_{:02d}.txt".format(outname[:-4],i))

# # CCF calculation
# def calc_ccf(N,box,step,dataa,datab,out):
# 	for i in range(len(dataa)):
# 		c=dataa[i].calc_ccf(datab[i],fp_flag.CIRCULANT,True)
# 		try: csum.add(c)
# 		except: csum=c
# 	return csum
# #	csum.process_inplace("normalize.edgemean")
# #	csum.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.15})
# 	#out.put((N,csum))

			# for i,k in enumerate(sorted(csum2.keys())):
			# 	im=csum2[k]
			# #	norm=im[BOX/2,BOX/2]/csum[BOX/2,BOX/2]
			# #	norm=im.get_clip(Region(BOX/2-5,BOX/2-5,11,11))["mean"]/csum.get_clip(Region(BOX/2-5,BOX/2-5,11,11))["mean"]
			# #	im.write_image("aa.hdf",i)

			# # This has been disabled since it eliminates the peak for zero shift. Instead we try the zero/zero elimination hack
			# 	norm=1.0
			# 	im.sub(csum*norm)

			# 	# set the 0,0 peak to the average of neighboring pixels to reduce fixed pattern noise issues (this worked poorly)
			# 	# im[BOX/2,BOX/2]=(im[BOX/2-1,BOX/2]+im[BOX/2+1,BOX/2]+im[BOX/2,BOX/2+1]+im[BOX/2,BOX/2-1])
			# #	s=im.process("math.sub.optimal",{"ref":csum,"ctfweight":0})

			# #	im.write_image("a.hdf",i+1)
			# 	# This is critical. Without this, after filtering we get too many false peaks
			# 	thr=im["mean"]+im["sigma"]*1.5
			# 	im.process_inplace("threshold.belowtozero",{"minval":thr})


			# # we start with a heavy filter, optimize, then repeat for successively less filtration
			# for scale in [0.02,0.04,0.07,0.1,0.5]:
			# 	csum3={k:csum2[k].process("filter.lowpass.gauss",{"cutoff_abs":scale}) for k in csum2.keys()}

			# 	incr=[16]*len(locs)
			# 	incr[-1]=incr[-2]=4	# if step is zero for last 2, it gets stuck as an outlier, so we just make the starting step smaller
			# 	simp=Simplex(qual,locs,incr,data=csum3)
			# 	locs=simp.minimize(maxiters=int(100/scale),epsilon=.01)[0]
			# 	locs=[int(floor(i*10+.5))/10.0 for i in locs]
			# 	print locs
			# 	if options.verbose > 7:
			# 		out=file("{}_path_{:02d}.txt".format(outname[:-4],int(1.0/scale)),"w")
			# 		for i in xrange(0,len(locs),2): out.write("%f\t%f\n"%(locs[i],locs[i+1]))

			# # compute the quality of each frame
			# quals=[0]*n			# quality of each frame based on its correlation peak summed over all images
			# cen=csum2[(0,1)]["nx"]/2
			# for i in xrange(n-1):
			# 	for j in xrange(i+1,n):
			# 		val=csum2[(i,j)].sget_value_at_interp(int(cen+locs[j*2]-locs[i*2]),int(cen+locs[j*2+1]-locs[i*2+1]))*sqrt(float(n-fabs(i-j))/n)
			# 		quals[i]+=val
			# 		quals[j]+=val


		# show a little movie of 5 averaged frames

		# if options.movie>0 :
		# 	mov=[]
		# 	for i in xrange(options.movie+1,len(outim)):
		# 		im=sum(outim[i-options.movie-1:i])
		# 		#im.process_inplace("filter.lowpass.gauss",{"cutoff_freq":.02})
		# 		mov.append(im)

		# 	display(mov)

			#mov2=[]
			#for i in xrange(0,len(outim)-10,2):
				#im=sum(outim[i+5:i+10])-sum(outim[i:i+5])
				#mov2.append(im)

			#display(mov2)

			#mov=[i.get_clip(Region(1000,500,2048,2048)) for i in mov]
			#s=sum(mov)
#			fsc=[i.calc_fourier_shell_correlation(s)[1025:2050] for i in mov]
#			plot(fsc)

			#csum=sum(csum2.values())
			#csum.mult(1.0/len(csum2))
			#csum.process_inplace("normalize.edgemean")
			#display(csum)
			#csum.write_image("a.hdf",0)

# def qual(locs,ccfs):
# 	"""computes the quality of the current alignment. Passed a dictionary of CCF images keyed by (i,j) tuple and
# 	an (x0,y0,x1,y1,...)  shift array. Smaller numbers are better since that's what the simplex does"""
# 	nrg=0.0
# 	cen=ccfs[(0,1)]["nx"]/2
# 	n=len(locs)/2
# 	for i in xrange(n-1):
# 		for j in xrange(i+1,n):
# #			nrg-=ccfs[(i,j)].sget_value_at_interp(int(cen+locs[j*2]-locs[i*2]),int(cen+locs[j*2+1]-locs[i*2+1]))
# 			# This is a recognition that we will tend to get better correlation with near neighbors in the sequence
# 			nrg-=ccfs[(i,j)].sget_value_at_interp(int(cen+locs[j*2]-locs[i*2]),int(cen+locs[j*2+1]-locs[i*2+1]))*sqrt(float(n-fabs(i-j))/n)
# 	return nrg

		# Generates different possibilites for resolution-weighted, but unaligned, averages

		# xy=XYData()
		# xy.set_size(2)
		# xy.set_x(0,0)
		# xy.set_y(0,1.0)
		# xy.set_x(1,0.707)
		# xy.set_y(1,0.0)
# 		if options.avgs :
# 			if options.verbose : print "Weighted average"
# 			normim=EMData(nx/2+1,ny)
# 			avgr=Averagers.get("weightedfourier",{"normimage":normim})
# 			for i in xrange(min(len(outim),25)):						# only use the first second for the unweighted average
# 				if options.verbose:
# 					print " {}/{}   \r".format(i+1,len(outim)),
# 					sys.stdout.flush()
# 				xy.set_y(1,1.0)					# no weighting
# 				outim[i]["avg_weight"]=xy
# 				avgr.add_image(outim[i])
# 			print ""

# 			av=avgr.finish()
# 			av.write_image(outname[:-4]+"_a.hdf",0)
# #			display(normim)

# 			# linear weighting with shifting 0 cutoff

# 			xy.set_y(1,0.0)
# 			for i in xrange(len(outim)):
# 				if options.verbose:
# 					print " {}/{}   \r".format(i+1,len(outim)),
# 					sys.stdout.flush()
# 				xy.set_x(1,0.025+0.8*(len(outim)-i)/len(outim))
# 				outim[i]["avg_weight"]=xy
# 				avgr.add_image(outim[i])
# 			print ""

# 			av=avgr.finish()
# 			av.write_image(outname[:-4]+"_b.hdf",0)

# 			# exponential falloff with shifting width

# 			xy.set_size(64)
# 			for j in xrange(64): xy.set_x(j,0.8*j/64.0)
# 			for i in xrange(len(outim)):
# 				if options.verbose:
# 					print " {}/{}   \r".format(i+1,len(outim)),
# 					sys.stdout.flush()
# 				for j in xrange(64) : xy.set_y(j,exp(-j/(3.0+48.0*(len(outim)-i)/float(len(outim)))))
# #				plot(xy)
# 				outim[i]["avg_weight"]=xy
# 				avgr.add_image(outim[i])
# 			print ""

# 			av=avgr.finish()
# 			av.write_image(outname[:-4]+"_c.hdf",0)

# def calcfsc(map1,map2):
# 	fsc=map1.calc_fourier_shell_correlation(map2)
# 	third=len(fsc)/3
# 	xaxis=fsc[0:third]
# 	fsc=fsc[third:2*third]

# 	return(xaxis,fsc)

# preprocess regions by normalizing and doing FFT
# def split_fft(img,i,box,step,out):
# 	lst=[]
# 	nx = img["nx"]
# 	ny = img["ny"]
# 	for dx in range(box/2,nx-box,step):
# 		for dy in range(box/2,ny-box,step):
# 			lst.append(img.get_clip(Region(dx,dy,box,box)).process("normalize.edgemean").do_fft())
# 	out.put((i,lst))

# def align(s1,s2,guess=(0,0),localrange=192,verbose=0):
# 	"""Aligns a pair of images, and returns a (dx,dy,Z) tuple. Z is the Z-score of the best peak, not a shift.
# 	The search will be limited to a region of +-localrange/2 about the guess, a (dx,dy) tuple. Resulting dx,dy
# 	is relative to the initial guess. guess and return both indicate the shift required to bring s2 in register
# 	with s1"""

# 	# reduce region used for alignment a bit (perhaps a lot for superresolution imaging

# 	guess=(int(guess[0]),int(guess[1]))
# 	if localrange<5 : localrange=192
# 	newbx=good_boxsize(min(s1["nx"],s1["ny"],4096)*0.8,larger=False)
# 	s1a=s1.get_clip(Region((s1["nx"]-newbx)/2,(s1["ny"]-newbx)/2,newbx,newbx))
# 	s2a=s2.get_clip(Region((s2["nx"]-newbx)/2-guess[0],(s2["ny"]-newbx)/2-guess[1],newbx,newbx))

# #	s1a.process_inplace("math.xystripefix",{"xlen":200,"ylen":200})
# 	s1a.process_inplace("filter.xyaxes0",{"neighbor":1})
# #	s1a.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.05})
# #	s1a.process_inplace("threshold.compress",{"value":0,"range":s1a["sigma"]/2.0})
# 	s1a.process_inplace("filter.highpass.gauss",{"cutoff_abs":.002})

# #	s2a.process_inplace("math.xystripefix",{"xlen":200,"ylen":200})
# #	s2a.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.05})
# 	s2a.process_inplace("filter.xyaxes0",{"neighbor":1})
# 	s2a.process_inplace("filter.highpass.gauss",{"cutoff_abs":.002})

# 	tot=s1a.calc_ccf(s2a)
# 	tot.process_inplace("xform.phaseorigin.tocenter")
# 	tot.process_inplace("normalize.edgemean")

# 	if verbose>2 :
# 		display((s1a,s2a,tot))

# 	if verbose>3 : display((s1a,s2a,tot),force_2d=True)

# 	dx,dy=(tot["nx"]/2-int(guess[0]),tot["ny"]/2-int(guess[1]))					# the 'false peak' should always be at the origin, ie - no translation
# 	mn=(tot[dx-2,dy-2]+tot[dx+2,dy+2]+tot[dx-2,dy+2]+tot[dx+2,dy-2])/4.0
# #	tot[dx,dy]=mn
# 	for x in xrange(dx-1,dx+2):
# 		for y in xrange(dy-1,dy+2):
# 			tot[x,y]=mn		# exclude from COM
# 			pass

# 	# first pass to have a better chance at finding the first peak, using a lot of blurring

# 	tot2=tot.get_clip(Region(tot["nx"]/2-96,tot["ny"]/2-96,192,192))
# 	tot2.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.04})		# This is an empirical value. Started with 0.04 which also seemed to be blurring out high-res features.
# 	tot2=tot2.get_clip(Region(tot2["nx"]/2-localrange/2,tot2["ny"]/2-localrange/2,localrange,localrange))
# 	dx1,dy1,dz=tot2.calc_max_location()
# 	dx1-=localrange/2
# 	dy1-=localrange/2

# 	# second pass with less blurring to fine tune it
# 	tot=tot.get_clip(Region(tot["nx"]/2-12+dx1,tot["ny"]/2-12+dy1,24,24))
# 	tot.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.12})		# This is an empirical value. Started with 0.04 which also seemed to be blurring out high-res features.
# 	dx,dy,dz = tot.calc_max_location()

# 	dev = tot["sigma"]

# 	if dev == 0.0 :
# 		dev  = 1.0
# 		print "Warning: sigma is zero in 'align' in iterative step for guess (", guess[0], ",", guess[1], ")."

# 	zscore = tot[dx,dy] / dev		# a rough Z score for the peak
# 	dx -= 12
# 	dy -= 12

# 	tot.write_image("tot.hdf",-1)

# 	if verbose>1: print "{},{} + {},{}".format(dx1,dy1,dx,dy)
# 	if verbose>2: display(tot)

# 	return dx1+dx+guess[0],dy1+dy+guess[1],zscore

# def align_subpixel(s1,s2,guess=(0,0),localrange=192,verbose=0):
# 	"""Aligns a pair of images to 1/4 pixel precision, and returns a (dx,dy,Z) tuple. Z is the Z-score of the best peak, not a shift.
# 	The search will be limited to a region of +-localrange/2 about the guess, a (dx,dy) tuple. Resulting dx,dy
# 	is relative to the initial guess. guess and return both indicate the shift required to bring s2 in register
# 	with s1"""

# 	# reduce region used for alignment a bit (perhaps a lot for superresolution imaging
# 	guess=(int(guess[0]*2.0),int(guess[1]*2.0))
# 	localrange*=2
# 	if localrange<5 : localrange=192*2
# 	newbx=good_boxsize(min(s1["nx"],s1["ny"],2048)*0.8,larger=False)
# 	newbx*=2
# 	s1a=s1.get_clip(Region((s1["nx"]-newbx)/2,(s1["ny"]-newbx)/2,newbx,newbx))
# 	s1a.scale(2)
# 	s2a=s2.get_clip(Region((s2["nx"]-newbx)/2-guess[0],(s2["ny"]-newbx)/2-guess[1],newbx,newbx))
# 	s2a.scale(2)

# #	s1a.process_inplace("math.xystripefix",{"xlen":200,"ylen":200})
# 	s1a.process_inplace("filter.xyaxes0",{"neighbor":1})
# #	s1a.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.05})
# #	s1a.process_inplace("threshold.compress",{"value":0,"range":s1a["sigma"]/2.0})
# 	s1a.process_inplace("filter.highpass.gauss",{"cutoff_abs":.002})

# #	s2a.process_inplace("math.xystripefix",{"xlen":200,"ylen":200})
# #	s2a.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.05})
# 	s2a.process_inplace("filter.xyaxes0",{"neighbor":1})
# 	s2a.process_inplace("filter.highpass.gauss",{"cutoff_abs":.002})

# 	tot=s1a.calc_ccf(s2a)
# 	tot.process_inplace("xform.phaseorigin.tocenter")
# 	tot.process_inplace("normalize.edgemean")

# 	if verbose>3 : display((s1a,s2a,tot),force_2d=True)

# 	mn=tot["mean"]
# 	dx,dy=(tot["nx"]/2,tot["ny"]/2)					# the 'false peak' should always be at the origin, ie - no translation
# 	for x in xrange(dx-2,dx+3):
# 		for y in xrange(dy-2,dy+3):
# 			tot[x,y]=mn		# exclude from COM
# #			pass

# 	# first pass to have a better chance at finding the first peak, using a lot of blurring
# 	tot2=tot.get_clip(Region(tot["nx"]/2-localrange/2,tot["ny"]/2-localrange/2,localrange,localrange))
# 	tot2.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.04})	# This is an empirical value. Started with 0.04 which also seemed to be blurring out high-res features.
# 	dx1,dy1,dz=tot2.calc_max_location()
# 	dx1-=localrange/2
# 	dy1-=localrange/2

# 	# second pass with less blurring to fine tune it
# 	tot=tot.get_clip(Region(tot["nx"]/2-24+dx1,tot["ny"]/2-24+dy1,48,48))
# 	tot.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.08})		# This is an empirical value. Started with 0.04 which also seemed to be blurring out high-res features.
# 	dx,dy,dz = tot.calc_max_location()

# 	dev = tot["sigma"]

# 	if dev == 0.0 :
# 		dev  = 1.0
# 		print "Warning: sigma is zero in 'align_subpixel' in iterative step for guess (", guess[0], ",", guess[1], ")."

# 	zscore = tot[dx,dy] / dev 	# a rough Z score for the peak
# 	dx -= 24
# 	dy -= 24

# 	if verbose>1: print "{},{} + {},{}".format(dx1,dy1,dx,dy)
# 	if verbose>2: display(tot)

# 	return (dx1+dx+guess[0])/2.0,(dy1+dy+guess[1])/2.0,zscore
