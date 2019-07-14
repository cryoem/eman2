#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

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

from past.utils import old_div
from future import standard_library
standard_library.install_aliases()
from builtins import range
from EMAN2 import *
from numpy import *
import pprint
import sys
import os
from sys import argv
from time import sleep,time,ctime
import threading
import queue
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
	parser.add_pos_argument(name="movies",help="List the tilt modies to align. Note: You must specify images in tilt angle order (negative to positive).", default="", guitype='filebox', browser="EMMovieDataTable(withmodal=True,multiselect=True)",  row=0, col=0,rowspan=1, colspan=3, mode="tomo")

	parser.add_header(name="orblock1", help='Just a visual separation', title="Tomography", row=1, col=0, rowspan=1, colspan=1, mode="tomo")

	#parser.add_argument("--rawtlt", default="", help="Specify a text file containing tilt angles that correspond to the input movies in alphabetical/numerical order.",guitype='filebox',browser="EMMovieRefsTable(withmodal=True,multiselect=False)", row=3, col=0, rowspan=1, colspan=3,mode='tomo')
	parser.add_argument("--tiltseries_name", default="", help="Specify a name for the tilt series to be generated from the input movies.",guitype='strbox', row=4, col=0, rowspan=1, colspan=1,mode='tomo')
	parser.add_argument("--tomo", default=False, help="Process input movies as tilts from a tomogram.",action="store_true", guitype='boolbox', row=4, col=2, rowspan=1, colspan=1, mode='tomo[True]')

	parser.add_header(name="orblock1", help='Just a visual separation', title="Corrections", row=5, col=0, rowspan=1, colspan=1, mode="align,tomo")

	parser.add_argument("--dark",type=str,default="",help="Perform dark image correction using the specified image file",guitype='filebox',browser="EMMovieRefsTable(withmodal=True,multiselect=True)", row=7, col=0, rowspan=1, colspan=3, mode="align,tomo,refs")
	parser.add_argument("--gain",type=str,default="",help="Perform gain image correction using the specified image file",guitype='filebox',browser="EMMovieRefsTable(withmodal=True,multiselect=True)", row=8, col=0, rowspan=1, colspan=3, mode="align,tomo,refs")

	parser.add_argument("--ref_label",type=str,default="",help="Optional: Specify a label for the averaged dark and gain references when using multiple, individual frames.\nA labeled will be written as movierefs/dark_<label>.hdf and movierefs/gain_<label>.hdf.\nNote: This option is ignored when using a single reference image/stack.",guitype='strbox', row=9, col=0, rowspan=1, colspan=3, mode="refs")
	parser.add_argument("--rotate_dark",  default = "0", type=str, choices=["0","90","180","270"], help="Rotate dark reference by 0, 90, 180, or 270 degrees. Default is 0. Transformation order is rotate then reverse.",guitype='combobox', choicelist='["0","90","180","270"]', row=10, col=0, rowspan=1, colspan=1, mode="refs")
	parser.add_argument("--reverse_dark", default=False, help="Flip dark reference along y axis. Default is False. Transformation order is rotate then reverse.",action="store_true",guitype='boolbox', row=10, col=1, rowspan=1, colspan=1, mode="refs")
	parser.add_argument("--k2", default=False, help="Perform gain image correction on gain images from a Gatan K2. Note, these are the reciprocal of typical DDD gain images.",action="store_true",guitype='boolbox', row=11, col=0, rowspan=1, colspan=1, mode="refs")
	parser.add_argument("--de64", default=False, help="Perform gain image correction on DE64 data. Note, these should not be normalized.",action="store_true",guitype='boolbox', row=11, col=1, rowspan=1, colspan=1, mode="refs")
	parser.add_argument("--rotate_gain", default = 0, type=str, choices=["0","90","180","270"], help="Rotate gain reference by 0, 90, 180, or 270 degrees. Default is 0. Transformation order is rotate then reverse.",guitype='combobox', choicelist='["0","90","180","270"]', row=12, col=0, rowspan=1, colspan=1, mode="refs")
	parser.add_argument("--reverse_gain", default=False, help="Flip gain reference along y axis (about x axis). Default is False. Transformation order is rotate then reverse.",action="store_true",guitype='boolbox', row=12, col=1, rowspan=1, colspan=1, mode="refs")
	parser.add_argument("--gain_darkcorrected", default=False, help="Do not dark correct gain image. False by default.",action="store_true",guitype='boolbox', row=13, col=0, rowspan=1, colspan=1, mode="refs")
	parser.add_argument("--invert_gain", default=False, help="Use reciprocal of input gain image",action="store_true",guitype='boolbox', row=13, col=1, rowspan=1, colspan=1, mode="refs")

	parser.add_argument("--bad_columns", type=str, help="Comma separated list of camera defect columns",default="", guitype='strbox', row=14, col=0, rowspan=1, colspan=3, mode="align,tomo")
	parser.add_argument("--bad_rows", type=str, help="Comma separated list of camera defect rows",default="", guitype='strbox', row=15, col=0, rowspan=1, colspan=3, mode="align,tomo")

	parser.add_header(name="orblock4", help='Just a visual separation', title="Alignment: ", row=16, col=0, rowspan=2, colspan=1, mode="align,tomo")

	parser.add_argument("--align_frames", action="store_true",help="Perform whole-frame alignment of the input stacks",default=False, guitype='boolbox', row=18, col=0, rowspan=1, colspan=1, mode='align[True],tomo[False]')
	parser.add_argument("--realign", action="store_true",help="Align frames using previous alignment parameters.",default=False, guitype='boolbox', row=18, col=1, rowspan=1, colspan=1, mode='align[False],tomo[False]')

	parser.add_argument("--round", choices=["float","int"],help="If float (default), apply subpixel frame shifts. If integer, use integer shifts.",default="float",guitype='combobox', choicelist='["float","integer"]', row=18, col=2, rowspan=1, colspan=1, mode="align,tomo")

	parser.add_argument("--noali", default=False, help="Average of non-aligned frames.",action="store_true", guitype='boolbox', row=19, col=0, rowspan=1, colspan=1, mode="align,tomo[True]")
	parser.add_argument("--allali", default=False, help="Average of all aligned frames.",action="store_true", guitype='boolbox', row=19, col=1, rowspan=1, colspan=1, mode='align[True],tomo[False]')
	parser.add_argument("--rangeali", default="", help="Average frames 'n1-n2'",type=str, guitype='strbox', row=19, col=2, rowspan=1, colspan=1, mode="align,tomo")
	parser.add_argument("--goodali", default=False, help="Average of good aligned frames.",action="store_true", guitype='boolbox', row=20, col=0, rowspan=1, colspan=1, mode="align,tomo")
	parser.add_argument("--bestali", default=False, help="Average of best aligned frames.",action="store_true", guitype='boolbox', row=20, col=1, rowspan=1, colspan=1, mode="align,tomo")
	# parser.add_argument("--ali4to14", default=False, help="Average of frames from 4 to 14.",action="store_true",guitype='boolbox', row=13, col=2, rowspan=1, colspan=1, mode="align,tomo")

	parser.add_argument("--groupby", type=int,help="Combine every N frames using a moving window approach.",default=1,  guitype='intbox', row=20, col=2, rowspan=1, colspan=1, mode="align,tomo")

	parser.add_header(name="orblock6", help='Just a visual separation', title="Optimization: ", row=21, col=0, rowspan=2, colspan=3, mode="align,tomo")

	parser.add_argument("--optbox", type=int,help="Box size to use during alignment optimization. Default is 512.",default=512, guitype='intbox', row=23, col=0, rowspan=1, colspan=1, mode="align,tomo")
	parser.add_argument("--optstep", type=int,help="Step size to use during alignment optimization. Default is 448.",default=448,  guitype='intbox', row=23, col=1, rowspan=1, colspan=1, mode="align,tomo")
	parser.add_argument("--optalpha", type=float,help="Penalization to apply during robust regression. Default is 0.1. If 0.0, unpenalized least squares will be performed (i.e., no trajectory smoothing).",default=0.1, guitype='floatbox', row=23, col=2, rowspan=1, colspan=1, mode="align,tomo")
	parser.add_argument("--optccf",default="robust",type=str, choices=["robust","centerofmass","ccfmax"],help="Use this approach to determine relative frame translations.\nNote: 'robust' utilizes a bimodal Gaussian to robustly determine CCF peaks between pairs of frames in the presence of a fixed background.", guitype='combobox', row=24, col=0, rowspan=1, colspan=2, mode='align["robust"],tomo["robust"]',choicelist='["robust","centerofmass","ccfmax"]')

	parser.add_header(name="orblock5", help='Just a visual separation', title="Optional: ", row=25, col=0, rowspan=2, colspan=3, mode="align,tomo")

	parser.add_argument("--step",type=str,default="0,1",help="Specify <first>,<step>,[last]. Processes only a subset of the input data. ie- 0,2 would process all even particles. Same step used for all input files. [last] is exclusive. Default= 0,1",guitype='strbox', row=27, col=0, rowspan=1, colspan=1, mode="align,tomo")
	parser.add_argument("--frames", default=False, help="Write corrected stack of frames.",action="store_true", guitype='boolbox', row=27, col=1, rowspan=1, colspan=1, mode="align,tomo")
	parser.add_argument("--ext",default="hdf",type=str, choices=["hdf","mrcs","mrc"],help="Save frames with this extension. Default is 'hdf'.", guitype='strbox', row=27, col=2, rowspan=1, colspan=1, mode="align,tomo")

	parser.add_argument("--threads", default=4,type=int,help="Number of threads to run in parallel. The default is 4, and our alignment routine requires 2+ threads. Using more threads will result in faster processing times.", guitype='intbox', row=28, col=0, rowspan=1, colspan=1, mode="align,tomo")

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=4, help="verbose level [0-9], higher number means higher level of verboseness",guitype="intbox",row=28,col=1,rowspan=1,colspan=1,mode="align,tomo")
	parser.add_argument("--debug", default=False, action="store_true", help="run with debugging output")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)

	parser.add_argument("--fixbadpixels",action="store_true",default=False,help="Tries to identify bad pixels in the dark/gain reference, and fills images in with sane values instead", guitype='boolbox', row=17, col=1, rowspan=1, colspan=1, mode='align[True]')
	#parser.add_argument("--normaxes",action="store_true",default=False,help="Tries to erase vertical/horizontal line artifacts in Fourier space by replacing them with the mean of their neighboring values.",guitype='boolbox', row=17, col=2, rowspan=1, colspan=1, mode='align')
	#parser.add_argument("--frames",action="store_true",default=False,help="Save the dark/gain corrected frames. Note that frames will be overwritten if identical --suffix is already present.", guitype='boolbox', row=19, col=0, rowspan=1, colspan=1, mode="align,tomo")
	parser.add_argument("--suffix",type=str,default="proc",help="Specify a unique suffix for output frames. Default is 'proc'. Note that the output of --frames will be overwritten if identical suffix is already present.",guitype='strbox', row=28, col=2, rowspan=1, colspan=1, mode="align,tomo")
	#parser.add_argument("--highdose", default=False, help="Use this flag when aligning high dose data (where features in each frame can be distinguished visually).",action="store_true",guitype='boolbox', row=18, col=0, rowspan=1, colspan=1,mode='align')
	parser.add_argument("--phaseplate", default=False, help="Use this flag when aligning phase plate frames.",action="store_true",guitype='boolbox', row=18, col=1, rowspan=1, colspan=1,mode='align')
	# parser.add_argument('--import_movies', action="store_true",default=False,help="List the references to import into 'movies' directory without additional processing.", default="", guitype='boolbox', row=0, col=0,rowspan=1, colspan=3, mode="import[True],default[False]")

	(options, args) = parser.parse_args()

	if len(args)<1:
		print(usage)
		parser.error("Specify input DDD stack to be processed.")

	if options.frames == False and options.noali == False and options.align_frames == False:
		print("No outputs specified. See --frames, --noali, or --align_frames. Exiting.")
		sys.exit(1)

	if options.align_frames == True:
		if options.allali == False and options.rangeali == False and options.goodali == False and options.bestali == False:# and options.ali4to14 == False:
			print("No post alignment outputs specified. Try with --allali, --rangeali, --goodali, --bestali. Exiting.")
			sys.exit(1)

	# moviesdir = os.path.join(".","movies")
	# if not os.access(moviesdir, os.R_OK):
	# 	os.mkdir(moviesdir)
		#print("Error: Could not locate movie data. Please import movie stacks using e2import.py.")
		#sys.exit(1)

	if options.bad_columns == "":
		options.bad_columns = []
	else:
		try:
			options.bad_columns = [int(c) for c in options.bad_columns.split(",")]
		except:
			print("Error: --bad_columns contains nonnumeric input.")
			sys.exit(1)

	if options.bad_rows == "":
		options.bad_rows = []
	else:
		try:
			options.bad_rows = [int(r) for r in options.bad_rows.split(",")]
		except:
			print("Error: --bad_rows contains nonnumeric input.")
			sys.exit(1)

	if options.align_frames and options.realign:
		print("Error: Running --align_frames and --realign would remove any existing alignment.")
		print("If you wish to do so, simply run with --align_frames only. Otherwise, use --realign.")
		sys.exit(1)

	if options.dark != "" or options.gain != "":
		refsdir = os.path.join(".","movierefs")
		if not os.access(refsdir, os.R_OK):
			os.mkdir(refsdir)

	if len(options.dark.split(",")) > 1:
		if options.ref_label != "":
			newfile = "movierefs/gain_{}.lst".format(options.ref_label)
		else:
			count = 0
			for f in os.listdir("movierefs"):
				if "dark{}.hdf".format(count) in f: count += 1
			newfile = "movierefs/dark_{}.lst".format(count)
		run("e2proclst.py {} --create {}".format(options.dark.replace(","," "),newfile))
		options.dark = newfile

	if len(options.gain.split(",")) > 1:
		if options.ref_label != "":
			newfile = "movierefs/gain_{}.lst".format(options.ref_label)
		else:
			count = 0
			for f in os.listdir("movierefs"):
				if "gain{}.hdf".format(count) in f:
					count += 1
			newfile = "movierefs/gain_{}.lst".format(count)
		run("e2proclst.py {} --create {}".format(options.gain.replace(","," "),newfile))
		options.gain = newfile

	pid=E2init(sys.argv)

	if options.noali or options.allali or options.goodali or options.bestali or options.rangeali:
		if options.tomo:
			options.outdir = os.path.join(".","tiltseries")
		else:
			options.outdir = os.path.join(".","micrographs")
		if not os.access(options.outdir, os.R_OK):
			os.mkdir(options.outdir)

	step = options.step.split(",")

	if len(step) == 3 : last = int(step[2])
	else :	last = -1

	first = int(step[0])
	step  = int(step[1])

	if options.verbose : print("Range = {} - {}, Step = {}".format(first, last, step))

	if options.dark != "":
		print("Loading Dark Reference")
		if options.dark[-4:].lower() in (".mrc") :
			dark_hdr = EMData(options.dark,0,True)
			nx = dark_hdr["nx"]
			ny = dark_hdr["ny"]
			nd = dark_hdr["nz"]
			if nd == 1: dark = EMData(options.dark)
			else: dark=EMData(options.dark,0,False,Region(0,0,0,nx,ny,1))
		else:
			nd=EMUtil.get_image_count(options.dark)
			dark = EMData(options.dark,0)
			nx = dark["nx"]
			ny = dark["ny"]
		if nd>1:
			sigd=dark.copy()
			sigd.to_zero()
			a=Averagers.get("mean",{"sigma":sigd,"ignore0":1})
			print("Summing Dark Frames")
			for i in range(0,nd):
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
			#if options.debug: sigd.write_image(options.dark.rsplit(".",1)[0]+"_sig.hdf")
			if options.fixbadpixels:
				sigd.process_inplace("threshold.binary",{"value":old_div(sigd["sigma"],10.0)}) # Theoretically a "perfect" pixel would have zero sigma, but in reality, the opposite is true
				dark.mult(sigd) # mask non-varying pixels in dark reference (set to zero)
		#if options.debug: dark.write_image(options.dark.rsplit(".",1)[0]+"_sum.hdf")
		#else: dark.mult(1.0/99.0)
		dark.process_inplace("threshold.clampminmax.nsigma",{"nsigma":3.0})
		#dark2=dark.process("normalize.unitlen")

		if options.rotate_dark!="0" and dark != None:
			if options.rotate_dark == "90":
				dark.process_inplace("xform.transpose")
				dark.process_inplace("xform.reverse",{"axis":"y"})
			elif options.rotate_dark == "180":
				dark.process_inplace("xform.transpose")
				dark.process_inplace("xform.reverse",{"axis":"y"})
				dark.process_inplace("xform.transpose")
				dark.process_inplace("xform.reverse",{"axis":"y"})
			else:
				dark.process_inplace("xform.transpose")
				dark.process_inplace("xform.reverse",{"axis":"x"}) #reversing over X axis for 270 deg rotate
			# nrot = int(options.rotate_dark)/90
			# for i in range(nrot):
			# 	dark.process_inplace("xform.transpose")
			# 	dark.process_inplace("xform.reverse")#,{"axis":"x"})
			#dark.set_size(maxdim,maxdim)
			#dark.process_inplace("xform.scale",{"clip":maxdim,"scale":1})
			#dark.rotate(0,0,int(options.rotate_dark))

		if options.reverse_dark: dark.process_inplace("xform.reverse",{"axis":"y"})

		options.dark = "movierefs/{}.hdf".format(base_name(options.dark,nodir=True))
		dark.write_image(options.dark,0)

	else : dark=None

	if options.gain != "":
		print("Loading Gain Reference")
		if options.k2: gain=EMData(options.gain)
		else:
			if options.gain[-4:].lower() in (".mrc") :
				gain_hdr = EMData(options.gain,0,True)
				nx = gain_hdr["nx"]
				ny = gain_hdr["ny"]
				nd = gain_hdr["nz"]
				if nd == 1: gain = EMData(options.gain)
				else: gain=EMData(options.gain,0,False,Region(0,0,0,nx,ny,1))
			else:
				nd=EMUtil.get_image_count(options.gain)
				gain = EMData(options.gain,0)
				nx = gain["nx"]
				ny = gain["ny"]
			if nd>1:
				sigg=gain.copy()
				sigg.to_zero()
				a=Averagers.get("mean",{"sigma":sigg,"ignore0":1})
				print("Summing Gain Frames")
				for i in range(0,nd):
					if options.verbose:
						sys.stdout.write("({}/{})   \r".format(i+1,nd))
						sys.stdout.flush()
					if options.gain != "" and options.gain[-4:].lower() in (".mrc") :
						t=EMData(options.gain,0,False,Region(0,0,i,nx,ny,1))
					else:
						t=EMData(options.gain,i)
					#t.process_inplace("threshold.clampminmax.nsigma",{"nsigma":4.0,"tozero":1})
					t.process_inplace("threshold.clampminmax",{"minval":0,"maxval":t["mean"]+t["sigma"]*3.5,"tozero":1})
					a.add_image(t)
				gain=a.finish()
				#if options.debug: sigg.write_image(options.gain.rsplit(".",1)[0]+"_sig.hdf")
				if options.fixbadpixels:
					sigg.process_inplace("threshold.binary",{"value":old_div(sigg["sigma"],10.0)}) # Theoretically a "perfect" pixel would have zero sigma, but in reality, the opposite is true
					if dark!="":
						try: sigg.mult(sigd) # set bad pixels identified in dark reference to 0 in gain reference
						except: pass # exception: dark has only 1 frame
					gain.mult(sigg) # set bad pixels to 0 in gain reference.
			#if options.debug: gain.write_image(options.gain.rsplit(".",1)[0]+"_sum.hdf")
			if options.de64:
				gain.process_inplace( "threshold.clampminmax", { "minval" : gain[ 'mean' ] - 8.0 * gain[ 'sigma' ], "maxval" : gain[ 'mean' ] + 8.0 * gain[ 'sigma' ], "tomean" : True } )
			else:
				gain.process_inplace("math.reciprocal",{"zero_to":0.0})
			#gain.mult(1.0/99.0)
			#gain.process_inplace("threshold.clampminmax.nsigma",{"nsigma":3.0})

		if options.dark!="" and options.gain_darkcorrected == False: gain.sub(dark) # dark correct the gain-reference

		if options.de64:
			mean_val = gain["mean"]
			if mean_val <= 0.: mean_val=1.
			gain.process_inplace("threshold.belowtominval",{"minval":0.01,"newval":mean_val})

		gain.mult(old_div(1.0,gain["mean"]))

		if options.invert_gain: gain.process_inplace("math.reciprocal",{"zero_to":0.0})

		if options.rotate_gain!=0 and gain != None:
			if options.rotate_gain == 90:
				gain.process_inplace("xform.transpose")
				gain.process_inplace("xform.reverse",{"axis":"y"})
			elif options.rotate_gain == 180:
				gain.process_inplace("xform.transpose")
				gain.process_inplace("xform.reverse",{"axis":"y"})
				gain.process_inplace("xform.transpose")
				gain.process_inplace("xform.reverse",{"axis":"y"})
			else:
				gain.process_inplace("xform.transpose")
				gain.process_inplace("xform.reverse",{"axis":"x"}) #reversing over X axis for 270 deg rotate
			# nrot = int(options.rotate_gain)/90
			# for i in range(nrot):
			# 	gain.process_inplace("xform.transpose")
			# 	gain.process_inplace("xform.reverse")
			# dim = [gain["nx"],gain["ny"]]
			# maxdim = max(dim[0],dim[1])
			# gain.set_size(maxdim,maxdim)
			# gain.rotate(0,0,int(options.rotate_gain))
			# if options.rotate_gain == 0 or options.rotate_gain == 180:
			# 	gain.set_size(dim[0],dim[1])
			# else:
			# 	gain.set_size(dim[1],dim[0])

		if options.reverse_gain:
			gain.process_inplace("xform.reverse",{"axis":"y"})

		options.gain = "movierefs/{}.hdf".format(base_name(options.gain,nodir=True))
		gain.write_image(options.gain,0)

	else: gain=None

	if options.tomo:
		#with open(options.rawtlt) as tlt:
		#	angles = [a for a in np.loadtxt(tlt)]
		db=js_open_dict(info_name(options.tomo_name,nodir=True))
		#db["rawtlt_source"] = options.rawtlt
		#db["tilt_angles"] = angles
		if gain:
			db["ddd_gainref"] = options.gain
			# if options.rotate_gain:
			# 	db["ddd_rotate_gain"] = options.rotate_gain
			# if options.reverse_gain:
			# 	db["ddd_reverse_gain"] = options.reverse_gain
			# if options.invert_gain:
			# 	db["ddd_invert_gain"] = options.invert_gain
			# if options.gain_darkcorrected:
			# 	db["ddd_gain_darkcorrected"] = options.gain_darkcorrected
		if dark:
			db["ddd_darkref"] = options.dark
		# 	if options.rotate_dark:
		# 		db["ddd_rotate_dark"]=options.rotate_dark
		# 	if options.reverse_dark:
		# 		db["ddd_reverse_dark"]=options.reverse_dark
		# if options.de64:
		# 	db["ddd_de64"] = options.de64
		# if options.k2:
		# 	db["ddd_k2"] = options.k2
		if options.bad_rows:
			db["ddd_bad_rows"] = options.bad_rows
		if options.bad_columns:
			db["ddd_bad_columns"] = options.bad_columns
		#db["ddd_bad_pixel_file"] =
		#if options.fixbadpixels:
		#	db["ddd_fixbadpixels"] = options.fixbadpixels
		#for a in angles:
		#	db[a] = {}
		db.close()

	# the user may provide multiple movies to process at once
	for idx,fsp in enumerate(sorted(args)):
		print("Processing {}".format(base_name(fsp,nodir=True)))

		if options.tomo:
			# write reference image info to corresponding movie info.json files
			#angle = angles[idx]
			db=js_open_dict(info_name(options.tomo_name,nodir=True))
			#db[angle] = {"data_source":fsp}
		else:
			db=js_open_dict(info_name(fsp,nodir=True))
			db["data_source"]=fsp
			if gain:
				db["ddd_gainref"] = options.gain
				# if options.rotate_gain:
				# 	db["ddd_rotate_gain"] = options.rotate_gain
				# if options.reverse_gain:
				# 	db["ddd_reverse_gain"] = options.reverse_gain
				# if options.invert_gain:
				# 	db["ddd_invert_gain"] = options.invert_gain
				# if options.gain_darkcorrected:
				# 	db["ddd_gain_darkcorrected"] = options.gain_darkcorrected
			if dark:
				db["ddd_darkref"] = options.dark
			# 	if options.rotate_dark:
			# 		db["ddd_rotate_dark"]=options.rotate_dark
			# 	if options.reverse_dark:
			# 		db["ddd_reverse_dark"]=options.reverse_dark
			# if options.de64:
			# 	db["ddd_de64"] = options.de64
			# if options.k2:
			# 	db["ddd_k2"] = options.k2
			if options.bad_rows:
				db["ddd_bad_rows"] = options.bad_rows
			if options.bad_columns:
				db["ddd_bad_columns"] = options.bad_columns
			#if options.fixbadpixels:
			#	db["ddd_fixbadpixels"] = options.fixbadpixels
		db.close()

		hdr = EMData(fsp, 0, True)
		try: n = EMUtil.get_image_count(fsp)
		except: n = hdr["nz"]
		if n < 3 :
			print("ERROR: {} has only {} images. Min 3 required.".format(fsp, n))
			continue

		if last <= 0 : flast = n
		else : flast = last

		if flast > n : flast = n

		process_movie(options, fsp, dark, gain, first, flast, step, idx)

	print("Done")
	E2end(pid)


def process_movie(options,fsp,dark,gain,first,flast,step,idx):
	cwd = os.getcwd()

	# format outname
	if options.frames: outname="{}/{}_{}".format(cwd,base_name(fsp,nodir=True),options.suffix) #Output contents vary with options
	else: outname="{}/{}".format(cwd,base_name(fsp,nodir=True))

	if options.groupby > 1: outname = "{}_group{}".format(outname,options.groupby)

	if options.ext == "mrc": outname = "{}.mrcs".format(outname)
	else: outname = "{}.{}".format(outname,options.ext)

	# prepare to read file
	if fsp[-4:].lower() in (".mrc"):
		hdr=EMData(fsp,0,True)			# read header
		nx,ny=hdr["nx"],hdr["ny"]

	# bgsub and gain correct the stack
	outim=[]
	nfs = 0
	t = time()
	for ii in range(first,flast,step):
		if options.verbose:
			sys.stdout.write(" {}/{}   \r".format(ii-first+1,flast-first+1))
			sys.stdout.flush()

		if fsp[-4:].lower() in (".mrc") :
		#if fsp[-4:].lower() in (".mrc") :
			im=EMData(fsp,0,False,Region(0,0,ii,nx,ny,1))
		else: im=EMData(fsp,ii)

		if dark!=None : im.sub(dark)
		if gain!=None : im.mult(gain)
		#im.process_inplace("threshold.clampminmax",{"minval":0,"maxval":im["mean"]+im["sigma"]*3.5,"tozero":1})
		if options.de64: im.process_inplace( "threshold.clampminmax", { "minval" : im[ 'minimum' ], "maxval" : im[ 'mean' ] + 8.0 * im[ 'sigma' ], "tomean" : True } )
		#if options.fixbadpixels : im.process_inplace("threshold.outlier.localmean",{"sigma":3.5,"fix_zero":1}) # fixes clear outliers as well as values which were exactly zero

		#im.process_inplace("threshold.clampminmax.nsigma",{"nsigma":3.0})
#			im.mult(-1.0)
		#if options.normalize: im.process_inplace("normalize.edgemean")

		if options.bad_rows != [] or options.bad_columns != []:
			im = im.process("math.xybadlines",{"rows":options.bad_rows,"cols":options.bad_columns})

		outim.append(im)

	nfs_read = len(outim)

	if options.noali:
		out=qsum(outim)
		if options.tomo:
			alioutname = os.path.join(".","tiltseries","{}__noali.hdf".format(base_name(options.tomo_name,nodir=True)))
			out.write_image(alioutname,idx) #write out the unaligned average movie
		else:
			alioutname = os.path.join(".","micrographs","{}__noali.hdf".format(base_name(fsp,nodir=True)))
			out.write_image(alioutname,0) #write out the unaligned average movie

	# group frames by moving window
	if options.groupby > 1:
		print("Grouping frames")
		grouped = []
		for i in range(len(outim)):
			avgr = Averagers.get("mean")
			avgr.add_image_list(outim[i:i+options.groupby])
			avg = avgr.finish()
			if options.frames: avg.write_image(outname,i)
			grouped.append(avg)
		outim = grouped

	# # group frames
	# if options.groupby > 1:
	# 	print("Grouping frames")
	# 	grouped = []
	# 	for i in range(0,len(outim)-options.groupby,options.groupby):
	# 		avgr = Averagers.get("mean")
	# 		avgr.add_image_list(outim[i:i+options.groupby])
	# 		grouped.append(avg)

	# 		if options.frames: avg.write_image(outname,ii-first)
	# 	outim = grouped

	if options.frames and options.ext == "mrc":
		os.rename(outname,outname.replace(".mrcs",".mrc"))

	t1 = time()-t
	print("{:.1f} s".format(time()-t))

	if options.align_frames and not options.realign:

		start = time()

		n=len(outim)
		nx=outim[0]["nx"]
		ny=outim[0]["ny"]

		md = min(nx,ny)
		if md <= 2048:
			print("This program does not facilitate alignment of movies with frames smaller than 2048x2048 pixels.")
			sys.exit(1)

		print("{} frames read ({} x {}). Grouped by {}.".format(nfs_read,nx,ny,options.groupby,n))

		ccfs=queue.Queue(0)

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
		peak_locs=queue.Queue(0)
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
		avgr.add_image_list(list(csum2.values()))
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
			for l in list(peak_locs.keys()):
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

		m = old_div(n*(n-1),2)
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
		cen=old_div(options.optbox,2) #csum2[(0,1)]["nx"]/2
		for i in range(n-1):
			for j in range(i+1,n):
				val=csum2[(i,j)].sget_value_at_interp(int(cen+locs[j*2]-locs[i*2]),int(cen+locs[j*2+1]-locs[i*2+1]))*sqrt(old_div(float(n-fabs(i-j)),n))
				quals[i]+=val
				quals[j]+=val

		print("{:1.1f} s".format(time()-t0,n))

		runtime = time()-start
		print("Runtime: {:.1f} s".format(runtime))

		# print("{:1.1f} s\nShift images".format(time()-t0))
		# for i,im in enumerate(outim):
		# 	if options.round == "int":
		# 		dx = int(round(locs[i*2],0))
		# 		dy = int(round(locs[i*2+1],0))
		# 		im.translate(dx,dy,0)
		# 	else: # float by default
		# 		dx = float(locs[i*2])
		# 		dy = float(locs[i*2+1])
		# 		im.translate(dx,dy,0)

		# locs = traj.ravel()
		# quals=[0]*n # quality of each frame based on its correlation peak summed over all images
		# cen=options.optbox/2 #csum2[(0,1)]["nx"]/2
		# for i in xrange(n-1):
		# 	for j in xrange(i+1,n):
		# 		val=csum2[(i,j)].sget_value_at_interp(int(cen+locs[j*2]-locs[i*2]),int(cen+locs[j*2+1]-locs[i*2+1]))*sqrt(float(n-fabs(i-j))/n)
		# 		quals[i]+=val
		# 		quals[j]+=val

		# print("{:1.1f} s".format(time()-t0,n))

		# runtime = time()-start
		# print("Runtime: {:.1f} s".format(runtime))

		drs = []
		reldrs = []
		#io = "micrographs/{}_info.txt".format(base_name(fsp,nodir=True))
		#out=open(io,"w")
		#out.write("#i,dx,dy,dr,rel dr,qual\n")
		for i in range(1,n):
			dx,dy = traj[i]
			dxlast,dylast = traj[i-1]
			dr = hypot(dx,dy)
			drs.append(dr)
			reldr = hypot(dx-dxlast,dy-dylast)
			reldrs.append(reldr)
			#out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(i,dx,dy,dr,reldr,quals[i]))
		#out.close

		# store alignment parameters
		if options.tomo:
			db=js_open_dict(info_name(options.tomo_name,nodir=True))
			db[idx]["ddd_alignment_trans"]=[i for i in locs]
			db[idx]["ddd_alignment_qual"]=[q for q in quals]
			db[idx]["ddd_alignment_time"]=runtime
			db[idx]["ddd_alignment_precision"]=options.round
			db[idx]["ddd_alignment_optbox"]=options.optbox
			db[idx]["ddd_alignment_optstep"]=options.optstep
			db[idx]["ddd_alignment_optalpha"]=options.optalpha
		else:
			db=js_open_dict(info_name(fsp,nodir=True))
			db["ddd_alignment_trans"]=[i for i in locs]
			db["ddd_alignment_qual"]=[q for q in quals]
			db["ddd_alignment_time"]=runtime
			db["ddd_alignment_precision"]=options.round
			db["ddd_alignment_optbox"]=options.optbox
			db["ddd_alignment_optstep"]=options.optstep
			db["ddd_alignment_optalpha"]=options.optalpha
		db.close()

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
		#	plt.show()

	if options.align_frames or options.realign:
		try:
		# Load previous/current alignment params (or input translations) (BOX FORMAT, tab separated values):
			if options.tomo:
				db=js_open_dict(info_name(options.tomo_name,nodir=True))
				locs = db[idx]["ddd_alignment_trans"]
				db.close()
			else:
				db=js_open_dict(info_name(fsp,nodir=True))
				locs = db["ddd_alignment_trans"]
				db.close()

			# shift frames
			print("{:1.1f} s\nShift images".format(time()-t0))
			alioutname = os.path.join(".","micrographs","{}__movieali.hdf".format(base_name(fsp,nodir=True)))
			for im in range(len(outim)):
				if options.round == "int":
					dx = int(round(locs[im*2],0))
					dy = int(round(locs[im*2+1],0))
					outim[im].translate(dx,dy,0)
				else: # float by default
					dx = float(locs[im*2])
					dy = float(locs[im*2+1])
					from fundamentals import fshift
					outim[im] = fshift(outim[im],dx,dy)
					outim[im].write_image(alioutname,im) #write out the unaligned average movie
					#im.translate(dx,dy,0)
				if options.debug or options.verbose > 5:
					print("{}\t{}\t{}".format(im,dx,dy))

			#if options.normaxes:
			#	for f in outim:
			#		f.process_inplace("filter.xyaxes0",{"neighbor":1})
				# or try padding before averaging and clip result to original box size?

			if options.allali:
				out=qsum(outim)
				if options.tomo:
					alioutname = os.path.join(".","tiltseries","{}__allali.hdf".format(base_name(options.tomo_name,nodir=True)))
					out.write_image(alioutname,idx) #write out the unaligned average movie
				else:
					alioutname = os.path.join(".","micrographs","{}__allali.hdf".format(base_name(fsp,nodir=True)))
					out.write_image(alioutname,0) #write out the unaligned average movie

			if options.goodali:
				thr=(max(quals[1:])-min(quals))*0.4+min(quals)	# max correlation cutoff for inclusion
				best=[im for i,im in enumerate(outim) if quals[i]>thr]
				out=qsum(best)
				print("Keeping {}/{} frames".format(len(best),len(outim)))
				if options.tomo:
					alioutname = os.path.join(".","tiltseries","{}__goodali.hdf".format(base_name(options.tomo_name,nodir=True)))
					out.write_image(alioutname,idx) #write out the unaligned average movie
				else:
					alioutname = os.path.join(".","micrographs","{}__goodali.hdf".format(base_name(fsp,nodir=True)))
					out.write_image(alioutname,0) #write out the unaligned average movie

			if options.bestali:
				thr=(max(quals[1:])-min(quals))*0.6+min(quals)	# max correlation cutoff for inclusion
				best=[im for i,im in enumerate(outim) if quals[i]>thr]
				out=qsum(best)
				print("Keeping {}/{} frames".format(len(best),len(outim)))
				if options.tomo:
					alioutname = os.path.join(".","tiltseries","{}__bestali.hdf".format(base_name(options.tomo_name,nodir=True)))
					out.write_image(alioutname,idx) #write out the unaligned average movie
				else:
					alioutname = os.path.join(".","micrographs","{}__bestali.hdf".format(base_name(fsp,nodir=True)))
					out.write_image(alioutname,0) #write out the unaligned average movie

			# if options.ali4to14:
			# 	out=qsum(outim[4:14]) # skip the first 4 frames then keep 10
			# 	if options.tomo:
			# 		alioutname = os.path.join(".","tiltseries","{}__4-14.hdf".format(base_name(fsp)))
			# 		out.write_image(alioutname,idx) #write out the unaligned average movie
			# 	else:
			# 		alioutname = os.path.join(".","micrographs","{}__4-14.hdf".format(base_name(fsp)))
			# 		out.write_image(alioutname,0) #write out the unaligned average movie

			if len(options.rangeali)>0:
				rng=[int(i) for i in options.rangeali.split("-")]
				out=qsum(outim[rng[0]:rng[1]+1])
				if options.tomo:
					alioutname = os.path.join(".","tiltseries","{}__{}.hdf".format(base_name(options.tomo_name,nodir=True),"-".join(rng)))
					out.write_image(alioutname,idx) #write out the unaligned average movie
				else:
					alioutname = os.path.join(".","micrographs","{}__{}.hdf".format(base_name(fsp,nodir=True),"-".join(rng)))
					out.write_image(alioutname,0) #write out the unaligned average movie

		except:
			print("Error: Could not find prior alignment for {}. Exiting".format(fsp,nodir=True))

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

	# img.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.45})

	for dx in range(old_div(box,2),nx-box,step):
		for dy in range(old_div(box,2),ny-box,step):
			clp = img.get_clip(Region(dx,dy,box,box))

			#clp.process_inplace("math.fft.resample",{"n":4})
			#clp.process_inplace("math.meanshrink",{"n":2.})
			#clp.process_inplace("filter.highpass.gauss",{"cutoff_pixels":2})
			#if options.normalize: clp.process_inplace("normalize.edgemean")
			#if box >= 512:
			#	if not options.tomo:
			#		if img["apix_x"] == 1.0: # likely an image with incorrect header or simulated data
			#			clp.process_inplace("filter.highpass.gauss",{"cutoff_pixels":2})
			#		else:
			#clp.process_inplace("mask.soft",{"outer_radius":(nx-box/2)/2,"width":box/8})
			#clp.process_inplace("filter.xyaxes0",{"neighbor":1})
			#clp.process_inplace("filter.highpass.gauss",{"cutoff_pixels":2})
			#clp.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.45})
			clp.process_inplace("normalize")
			#if options.debug:
			#	clp["patch_id"] = patchid
			#	clp["frame_id"] = i
			#	clp.write_image("patch_{}.hdf".format(str(patchid).zfill(4)),i)
			#patchid += 1
			clp.do_fft_inplace()
			lst.append(clp)
	out.put((i,lst))

def correlation_peak_model(x_y, xo, yo, sigma, amp):
	x, y = x_y
	if sigma <= 0: return np.ones_like(x)*np.inf
	xo = float(xo)
	yo = float(yo)
	g = amp*np.exp(old_div(-(((x-xo)**2)+((y-yo)**2)),(2.*sigma**2)))
	return g.ravel()

def fixedbg_peak_model(x_y, sigma, amp):
	x, y = x_y
	if sigma <= 0: return np.ones_like(x)*np.inf
	xo = float(old_div(len(x),2))
	yo = float(old_div(len(y),2))
	g = amp*np.exp(old_div(-(((x-xo)**2)+((y-yo)**2)),(2.*sigma**2)))
	return g.ravel()

def twod_bimodal(x_y,x1,y1,sig1,amp1,sig2,amp2,offset):
	x, y = x_y
	if sig1 <= 0 or sig2 <= 0: return np.ones_like(x)*np.inf #np.zeros_like(x)
	#correlation_peak = correlation_peak_model((x,y),x1,y1,sig1,amp1)
	cp = amp1*np.exp(old_div(-(((x-x1)**2+(y-y1)**2)),(2.*sig1**2)))
	#fixedbg_peak = fixedbg_peak_model((x,y),sig2,amp2)
	fp = amp2*np.exp(old_div(-(((x-old_div(len(x),2.))**2+(y-old_div(len(y),2.))**2)),(2.*sig2**2)))
	return offset + cp.ravel() + fp.ravel() # + noise

# def edgemean(a,xc,yc):
# 	# mean of values around perimeter of image
# 	return (np.sum(a[:,0])+np.sum(a[0,1:])+np.sum(a[1:,-1])+np.sum(a[-1,1:-1]))/(4*(a.shape[0]-1))

def neighbormean_origin(a): # replace origin pixel with mean of surrounding pixels
	proc = a.copy()
	if a.shape[0] % 2 == 0:
		ac = old_div(proc.shape[0],2)
		r = proc[ac-2:ac+2,ac-2:ac+2].copy()
		rc = old_div(r.shape[0],2)
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
		ac = old_div(proc.shape[0],2)
		r = proc[ac-2:ac+1,ac-2:ac+1].copy()
		plt.imshow(r)
		rc = old_div(r.shape[0],2)
		r[rc,rc] = np.nan
		proc[ac,ac] = np.nanmean(r)
	return proc

def find_com(ccf): # faster alternative to gaussian fitting...less robust in theory.
	thresh = (np.ones(ccf.shape) * np.mean(ccf))+2.5*np.std(ccf)
	m = np.greater(ccf,thresh) * 1.0
	m = old_div(m, np.sum(m))
	# marginal distributions
	dx = np.sum(m, 1)
	dy = np.sum(m, 0)
	# expected values
	cx = np.sum(dx * np.arange(ccf.shape[0]))
	cy = np.sum(dy * np.arange(ccf.shape[1]))
	return cx-old_div(ccf.shape[0],2),cy-old_div(ccf.shape[1],2)

def bimodal_peak_model(options,ccf):
	nxx = ccf["nx"]
	bs = nxx/4

	xx = np.linspace(0,bs,bs)
	yy = np.linspace(0,bs,bs)
	xx,yy = np.meshgrid(xx,yy)

	r = Region(old_div(nxx,2)-old_div(bs,2),old_div(nxx,2)-old_div(bs,2),bs,bs)

	ccfreg = ccf.get_clip(r)
	ncc = ccfreg.numpy().copy()

	x1 = float(bs)/2.
	y1 = float(bs)/2.
	a1 = ncc.mean()
	s1 = 10.0
	a2 = ncc.max()
	s2 = 0.6
	offset=ncc.mean()

	if options.optccf == "centerofmass":
		#initial_guess = [x1,y1,s1,a1,s2,a2]
		# drop the origin (replace with perimeter mean)
		#mval = edgemean(ncc)
		#ix,iy = ncc.shape
		pncc = neighbormean_origin(ncc)
		x1,y1 = find_com(pncc)
		x1 = x1+old_div(nxx,2)
		y1 = y1+old_div(nxx,2)
		# if options.debug:
		# 	try:
		# 		ii = EMUtil.get_image_count("tmp.hdf")
		# 		from_numpy(pncc).write_image("tmp.hdf",ii)`
		# 	except:
		# 		from_numpy(pncc).write_image("tmp.hdf",0)
		return [x1,y1,s1,a1,s2,a2,offset],ccf.sget_value_at_interp(x1,y1)
		# try: # run optimization
		# 	bds = [(-np.inf, -np.inf,  0.01, 0.01),(np.inf, np.inf, 100.0, 100000.0)]
		# 	#bds = [(-bs/2, -bs/2,  0.0, 0.0),(bs/2, bs/2, 100.0, 100000.0)]
		# 	popt,pcov=optimize.curve_fit(correlation_peak_model,(xx,yy),ncc.ravel(),p0=initial_guess,bounds=bds,method='dogbox',max_nfev=50)#,xtol=0.1)#,ftol=0.0001,gtol=0.0001)
		# except:
		# 	return None, -1
	elif options.optccf == "ccfmax": # only useful for extremely high contrast frames
		yc,xc = np.where(ncc==ncc.max())
		popt = [float(xc[0]+old_div(nxx,2)),float(yc[0]+old_div(nxx,2)),ncc.max(),1.,0.,0.]
		return popt,ccf.sget_value_at_interp(popt[0],popt[1])
	elif options.optccf == "robust":
		initial_guess = [x1,y1,s1,a1,s2,a2,offset]
		bds = [(-bs/2, -bs/2, 0, 0, 0.6, 0.0,ncc.min()),(bs/2, bs/2, 100.0, 20000.0, 2.5, 100000.0,ncc.max())]
		try:
			popt,pcov=optimize.curve_fit(twod_bimodal,(xx,yy),ncc.ravel(),p0=initial_guess,bounds=bds)
			#popt,pcov=optimize.curve_fit(twod_bimodal,(xx,yy),ncc.ravel(),p0=initial_guess,bounds=bds,method="dogbox",max_nfev=250,xtol=1e-3,ftol=1e-6,loss='linear')
		except:
			try:
				x1,y1 = find_com(ncc)
				popt = [x1,y1,0.,0.,0.,0.,0.]
			except:
				return None,-1

		popt = [p for p in popt]
		popt[0] = popt[0] + old_div(nxx,2) - old_div(bs,2)
		popt[1] = popt[1] + old_div(nxx,2) - old_div(bs,2)
		popt[2] = np.abs(popt[2])

		return popt,ccf.sget_value_at_interp(popt[0],popt[1])
	# elif options.optccf == "robust":
	# 	initial_guess = [x1,y1,s1,a1]
	# 	bds = [(-bs/2, -bs/2, 0.6, ncc.min()),(bs/2, bs/2, 100.0, ncc.max())]
	# 	try:
	# 		popt,pcov=optimize.curve_fit(correlation_peak_model,(xx,yy),ncc.ravel(),p0=initial_guess,bounds=bds,method="dogbox",max_nfev=250,xtol=1e-3,ftol=1e-6,loss='linear')
	# 	except RuntimeError:
	# 		return None,-1

	# 	popt = [p for p in popt]
	# 	popt[0] = popt[0] + nxx/2 - bs/2
	# 	popt[1] = popt[1] + nxx/2 - bs/2
	# 	popt[2] = np.abs(popt[2])
	# 	popt.extend([s2,a2])
	# 	return popt,ccf.sget_value_at_interp(popt[0],popt[1])

def qsum(imlist):
	avg=Averagers.get("mean")
	avg.add_image_list(imlist)
	return avg.finish()

def qual(locs,ccfs):
	"""computes the quality of the current alignment. Passed a dictionary of CCF images keyed by (i,j) tuple and
	an (x0,y0,x1,y1,...)  shift array. Smaller numbers are better since that's what the simplex does"""
	nrg=0.0
	cen=old_div(ccfs[(0,1)]["nx"],2)
	n=old_div(len(locs),2)
	for i in range(n-1):
		for j in range(i+1,n):
			penalty = sqrt(old_div(float(n-fabs(i-j)),n))**2 # This is a recognition that we will tend to get better correlation with near neighbors in the sequence
			locx = int(cen+locs[j*2]-locs[i*2])
			locy = int(cen+locs[j*2+1]-locs[i*2+1])
			nrg-=ccfs[(i,j)].sget_value_at_interp(locx,locy)*penalty
	return nrg

def run(command):
	"Mostly here for debugging, allows you to control how commands are executed (os.system is normal)"

	print("{}: {}".format(ctime(time()),command))
	ret=launch_childprocess(command)

	# We put the exit here since this is what we'd do in every case anyway. Saves replication of error detection code above.
	if ret !=0 :
		print("Error running: ",command)
		sys.exit(1)

if __name__ == "__main__":
	main()
