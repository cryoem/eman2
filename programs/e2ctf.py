#!/usr/bin/env python

#
# Author: Steven Ludtke, 10/29/2008 (sludtke@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
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

# e2ctf.py  10/29/2008 Steven Ludtke
# This is a program for determining CTF parameters and (optionally) phase flipping images

from EMAN2 import *
from EMAN2db import db_open_dict, db_close_dict, db_check_dict, db_list_dicts
from optparse import OptionParser
from math import *
import os
import sys
import weakref
import traceback
from numpy import array,arange
import numpy
import threading

from Simplex import Simplex

debug=False
logid=None

sfcurve=None		# This will store a global structure factor curve if specified
sfcurve2=None
envelopes=[]		# simplex minimizer needs to use a global at the moment

def imcount(fsp):
	try: return EMUtil.get_image_count(fsp)
	except: return 0

def get_df(itm):
	return itm[1].defocus

def main():
	global debug,logid
	progname = os.path.basename(sys.argv[0])

	usage = """prog [options] <input stack/image> ...
Various CTF-related operations on images, including automatic fitting.Input particles should be unmasked and unfiltered. A
minimum of ~20 percent padding around the particles is required for background extraction, even if this brings the edge
of another particle into the box in some cases. Particles should be reasonably well centered. Can also optionally phase
flip and Wiener filter particles. Wiener filtration comes after phase-flipping, so if phase flipping is performed Wiener
filtered particles will also be phase-flipped. Note that both operations are performed on oversampled images if specified,
but this will produce phase-flipped images which are irreversable, so, while oversampling can be useful for fitting, it
is not recommended for phase-flipping.

Increasing padding during the particle picking process will improve the accuracy of phase-flipping, particularly for
images far from focus.

NOTE: This program should be run from the project directory, not from within the particles directory
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_pos_argument(name="particles",help="List the file to process with e2ctf here.", default="", guitype='filebox', browser="EMCTFParticlesTable(withmodal=True,multiselect=True)",  filecheck=False, row=0, col=0,rowspan=1, colspan=2, mode='autofit,tuning,genoutp,gensf')
#	parser.add_header(name="ctfheader", help='Options below this label are specific to e2ctflassaverage3d', title="### e2ctf options ###", default=None, row=1, col=0, rowspan=1, colspan=2, mode="autofit,tuning,genoutp,gensf")

	parser.add_argument("--allparticles",action="store_true",help="Will process all particle stacks stored in the particles subdirectory (no list of files required)",default=False, guitype='boolbox',row=1, col=0, mode='autofit,tuning,genoutp,gensf')
	parser.add_argument("--onlynew",action="store_true",help="Will skip any files for which __ctf_flip files already exist.",default=False)
	parser.add_argument("--sortdefocus",action="store_true",help="Sorts the micrographs in order by defocus",default=False,guitype='boolbox',row=3,col=1, mode='tuning')
	parser.add_argument("--minptcl",type=int,help="Files with fewer than the specified number of particles will be skipped",default=0,guitype='intbox', row=2, col=0, mode='autofit,tuning,genoutp,gensf')
	parser.add_argument("--minqual",type=int,help="Files with a quality value lower than specified will be skipped",default=0,guitype='intbox', row=2, col=1, mode='autofit,tuning,genoutp,gensf')
	parser.add_argument("--chunk",type=str,help="<chunksize>,<nchunk>. Will process files in groups of chunksize, and process the <nchunk>th group. eg - 100,3 will read files 300-399 ",default=None,guitype='strbox',row=1,col=1, mode='autofit,tuning,genoutp,gensf')
	parser.add_argument("--gui",action="store_true",help="Start the GUI for interactive fitting",default=False, guitype='boolbox', row=3, col=0, rowspan=1, colspan=1, mode="tuning[True]")
	parser.add_argument("--autofit",action="store_true",help="Runs automated CTF fitting on the input images",default=False, guitype='boolbox', row=8, col=0, rowspan=1, colspan=1, mode='autofit[True]')
	parser.add_argument("--wholeimage",action="store_true",help="Display an additional curve using the whole micrograph, not just particles.",default=False,guitype='boolbox',row=8,col=2, mode='autofit')
	parser.add_argument("--highdensity",action="store_true",help="If particles are very close together, this will interfere with SSNR estimation. If set uses an alternative strategy, but may over-estimate SSNR.",default=False, guitype='boolbox', row=9, col=2, rowspan=1, colspan=1, mode='autofit[False]')
	parser.add_argument("--zerook",action="store_true",help="Normally particles with zero value on the edge are considered to be bad. This overrides that behavior, primarily for simulated data.",default=False)
	parser.add_argument("--astigmatism",action="store_true",help="Includes astigmatism in automatic fitting",default=False, guitype='boolbox', row=8, col=1, rowspan=1, colspan=1, mode='autofit[False]')
	parser.add_argument("--curdefocushint",action="store_true",help="Rather than doing the defocus from scratch, use existing values in the project as a starting point",default=False, guitype='boolbox', row=9, col=0, rowspan=1, colspan=1, mode='autofit[True]')
	parser.add_argument("--curdefocusfix",action="store_true",help="Fixes the defocus at the current determined value (if any) (+-.001 um), but recomputes SSNR, etc.",default=False, guitype='boolbox', row=9, col=1, rowspan=1, colspan=1, mode='autofit[False]')
	parser.add_argument("--bgmask",type=int,help="Background is computed using a soft mask of the center/edge of each particle with the specified radius. Default radius is boxsize/2.6.",default=0)
	parser.add_argument("--fixnegbg",action="store_true",help="Will perform a final background correction to avoid slight negative values near zeroes")
	parser.add_argument("--computesf",action="store_true",help="Will determine the structure factor*envelope for the aggregate set of images", default=False, nosharedb=True, guitype='boolbox', row=9, col=0, rowspan=1, colspan=1, mode="gensf[True]")
	parser.add_argument("--apix",type=float,help="Angstroms per pixel for all images",default=0, guitype='floatbox', row=4, col=0, rowspan=1, colspan=1, mode="autofit['self.pm().getAPIX()']")
	parser.add_argument("--voltage",type=float,help="Microscope voltage in KV",default=0, guitype='floatbox', row=4, col=1, rowspan=1, colspan=1, mode="autofit['self.pm().getVoltage()']")
	parser.add_argument("--cs",type=float,help="Microscope Cs (spherical aberation)",default=0, guitype='floatbox', row=5, col=0, rowspan=1, colspan=1, mode="autofit['self.pm().getCS()']")
	parser.add_argument("--ac",type=float,help="Amplitude contrast (percentage, default=10)",default=10, guitype='floatbox', row=5, col=1, rowspan=1, colspan=1, mode='autofit')
	parser.add_argument("--defocusmin",type=float,help="Minimum autofit defocus",default=0.6, guitype='floatbox', row=6, col=0, rowspan=1, colspan=1, mode="autofit[0.6]")
	parser.add_argument("--defocusmax",type=float,help="Maximum autofit defocus",default=4, guitype='floatbox', row=6, col=1, rowspan=1, colspan=1, mode='autofit[4.0]')
	parser.add_argument("--constbfactor",type=float,help="Set B-factor to fixed specified value, negative value autofits",default=-1.0, guitype='floatbox', row=12, col=0, rowspan=1, colspan=1, mode='autofit[-1.0],tuning[-1.0],genoutp[-1.0]')
	parser.add_argument("--autohp",action="store_true",help="Automatic high pass filter of the SNR only to remove initial sharp peak, phase-flipped data is not directly affected (default false)",default=False, guitype='boolbox', row=7, col=0, rowspan=1, colspan=1, mode='autofit[True]')
	parser.add_argument("--invert",action="store_true",help="Invert the contrast of the particles in output files (default false)",default=False, guitype='boolbox', row=5, col=1, rowspan=1, colspan=1, mode='genoutp')
	parser.add_argument("--nonorm",action="store_true",help="Suppress per image real-space normalization",default=False)
	parser.add_argument("--nosmooth",action="store_true",help="Disable smoothing of the background (running-average of the log with adjustment at the zeroes of the CTF)",default=False, guitype='boolbox', row=7, col=1, rowspan=1, colspan=1, mode='autofit')
	parser.add_argument("--refinebysnr",action="store_true",help="Refines the defocus value by looking at the high resolution smoothed SNR. Requires good starting defocus. Important: also replaces the SNR with a smoothed version.",default=False, guitype='boolbox', row=3, col=0, rowspan=1, colspan=1, mode='genoutp')
	parser.add_argument("--phaseflip",action="store_true",help="Perform phase flipping after CTF determination and writes to specified file.",default=False, guitype='boolbox', row=4, col=0, rowspan=1, colspan=1, mode='genoutp[True]')
	parser.add_argument("--phaseflipproc",help="If specified _proc particles will be generated. Typical = filter.lowpass.gauss:cutoff_freq=.07",default=None, guitype='strbox', row=6, col=0, rowspan=1, colspan=3, mode='genoutp["filter.lowpass.gauss:cutoff_freq=.07"]')
	parser.add_argument("--phaseflipproc2",help="If specified _proc particles will be generated. Typical = filter.highpass.gauss:cutoff_freq=.005",default=None, guitype='strbox', row=7, col=0, rowspan=1, colspan=3, mode='genoutp["filter.highpass.gauss:cutoff_freq=.005"]')
	parser.add_argument("--phaseflipproc3",help="If specified _proc particles will be generated. Typical = math.meanshrink:n=2",default=None, guitype='strbox', row=8, col=0, rowspan=1, colspan=3, mode='genoutp["math.meanshrink:n=2"]')
	parser.add_argument("--phasefliphp",action="store_true",help="Perform phase flipping with auto-high pass filter",default=False, guitype='boolbox', row=5, col=0, rowspan=1, colspan=1, mode='genoutp')
	parser.add_argument("--wiener",action="store_true",help="Wiener filter (optionally phaseflipped) particles.",default=False, guitype='boolbox', row=4, col=1, rowspan=1, colspan=1, mode='genoutp[True]')
#	parser.add_argument("--virtualout",type=str,help="Make a virtual stack copy of the input images with CTF parameters stored in the header. BDB only.",default=None)
	parser.add_argument("--storeparm",action="store_true",help="Output files will include CTF info. CTF parameters are used from the database, rather than values that may be present in the input image header. Critical to use this when generating output !",default=False,guitype='boolbox', row=3, col=1, rowspan=1, colspan=1, mode='genoutp[True]')
	parser.add_argument("--oversamp",type=int,help="Oversampling factor",default=1, guitype='intbox', row=3, col=0, rowspan=1, colspan=2, mode='autofit')
	parser.add_argument("--classify",type=int,help="Highly experimental ! Subclassify particles (hopefully by defocus) into n groups.",default=0)
	parser.add_argument("--sf",type=str,help="The name of a file containing a structure factor curve. Specify 'none' to use the built in generic structure factor. Default=auto",default="auto",guitype='strbox',nosharedb=True,returnNone=True,row=12,col=1,rowspan=1,colspan=1, mode='autofit,tuning')
	parser.add_argument("--parallel", default=None, help="parallelism argument. This program supports only thread:<n>")
	parser.add_argument("--threads", default=1,type=int,help="Number of threads to run in parallel on a single computer when multi-computer parallelism isn't useful",guitype='intbox', row=12, col=2, rowspan=1, colspan=1, mode='autofit[1]')
	parser.add_argument("--debug",action="store_true",default=False)
	parser.add_argument("--dbds",type=str,default=None,help="Obsolete option for old e2workflow. Present only to provide warning messages.")
	parser.add_argument("--source_image",type=str,default=None,help="Filters particles only with matching ptcl_source_image parameters in the header")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()

	if options.dbds!=None :
		print "--dbds no longer supported, as this was part of the retired e2workflow interface. Exiting."
		sys.exit(1)

	if options.threads : nthreads=options.threads
	elif options.parallel!=None :
		if options.parallel[:7]!="thread:":
			print "ERROR: only thread:<n> parallelism supported by this program. It is i/o limited."
			sys.exit(1)
		nthreads=int(options.parallel[7:])
	else: nthreads=1

	if options.allparticles:
		args=["particles/"+i for i in os.listdir("particles") if "__" not in i and i[0]!="." and ".hed" not in i ]
		args.sort()
		if options.verbose : print "%d particle stacks identified"%len(args)

	if options.chunk!=None:
		print sys.argv
		ninchunk,nchunk=options.chunk.split(",")
		ninchunk=int(ninchunk)
		nchunk=int(nchunk)
		if options.verbose : print "{} stacks with chunks of {}".format(len(args),ninchunk)
		args=args[ninchunk*nchunk:ninchunk*(nchunk+1)]
		if options.verbose: print "{} stacks in specified chunk".format(len(args))
		nthreads=1		# no threading with chunks

	if options.onlynew:
		print "%d files to process"%len(args)
		dl=os.listdir("particles")
		args=[i for i in args if not info_name(i)+"__ctf_flip.hdf" in dl]
		if options.verbose: print "{} stacks after onlynew filter".format(len(args))


	if len(args)<1 : parser.error("Input image required")
	if options.autofit:
		if options.voltage==0 : parser.error("Please specify voltage")
		if options.cs==0 : parser.error("Please specify Cs")
	if options.apix==0 : print "Using A/pix from header"

	debug=options.debug
	img_sets=None

	global sfcurve,sfcurve2
	init_sfcurve(options.sf)

	logid=E2init(sys.argv, options.ppid)

 #	if options.oversamp>1 : options.apix/=float(options.oversamp)


	#db_project=db_open_dict("bdb:project")
	#db_parms=db_open_dict("bdb:e2ctf.parms")
#	db_misc=db_open_dict("bdb:e2ctf.misc")

	# remove any files that don't have enough particles from the list
	if options.minptcl>0 :
		args=[i for i in args if imcount(i)>=options.minptcl]
		if options.verbose: print "{} stacks after minptcl filter".format(len(args))


	# remove files with quality too low
	if options.minqual>0 :
		outargs=[]
		for i in args:
			try:
				js=js_open_dict(info_name(i))
				if js["quality"]>=options.minqual : outargs.append(i)
				js.close()
			except:
#				traceback.print_exc()
				print "Unknown quality for {}, including it".format(info_name(i))
				outargs.append(i)

		args=outargs

		if options.verbose: print "{} stacks after quality filter".format(len(args))

	options.filenames = args

	### Power spectrum and CTF fitting
	if nthreads>1:
		print "Fitting in parallel with ",nthreads," threads"
		chunksize=int(ceil(float(len(args))/nthreads))
#		print " ".join(sys.argv+["--chunk={},{}".format(chunksize,0)])
		threads=[threading.Thread(target=os.system,args=[" ".join(sys.argv+["--chunk={},{}".format(chunksize,i)])]) for i in xrange(nthreads)]
		for t in threads: t.start()
		for t in threads: t.join()
		print "Parallel fitting complete"
	else:
		img_sets=None
		if options.autofit:
			img_sets=pspec_and_ctf_fit(options,debug) # converted to a function so to work with the workflow

			if options.constbfactor>0:
				for i in img_sets: i[1].bfactor=options.constbfactor

	### GUI - user can update CTF parameters interactively
	if options.gui :
		if img_sets==None : img_sets = get_gui_arg_img_sets(options.filenames)

		if options.constbfactor>0:
			for i in img_sets: i[1].bfactor=options.constbfactor

		if options.sortdefocus:
			img_sets.sort(key=get_df)

		if len(img_sets) == 0:
			E2end(logid)
			exit(1)
		from emapplication import EMApp
		app=EMApp()
		gui=GUIctf(app,img_sets,options.autohp,options.nosmooth)
		gui.show_guis()
		app.exec_()

#		print "done execution"

	### Refine defocus and smooth SNR
	if options.refinebysnr:
		refine_and_smoothsnr(options,sfcurve2,debug=False)

	### Process input files
	if debug : print "Phase flipping / Wiener filtration"
	# write wiener filtered and/or phase flipped particle data to the local database
	if options.phaseflip or options.wiener or options.phasefliphp or options.phaseflipproc or options.storeparm: # only put this if statement here to make the program flow obvious
		write_e2ctf_output(options) # converted to a function so to work with the workflow

	if options.computesf :
		img_sets = get_gui_arg_img_sets(options.filenames)
		if options.constbfactor>0:
			for i in img_sets: i[1].bfactor=options.constbfactor
		print "Recomputing structure factor"
		envelope=compute_envelope(img_sets)

#		print envelope

		#db_close_dict("bdb:e2ctf.misc")

		out=file("strucfac.txt","w")
		for i in envelope: out.write("%f\t%f\n"%(i[0],i[1]))
		out.close()


	E2end(logid)

def init_sfcurve(opt):

	global sfcurve,sfcurve2,hasgoodsf

	sfcurve2=None
	if opt != None and opt.lower()!="none" and opt.lower()!="auto" :
		sfcurve2=XYData()
		sfcurve2.read_file(opt)
		hasgoodsf=True
		# We used to pass in log10(sf), but no longer do this
		#for i in range(sfcurve2.get_size()):
			#v=sfcurve2.get_y(i)
			#if v<=0 :
				#print "Warning values <=0 found in structure factor file. Please remove."
				#sfcurve2.set_y(i,-10.0)
			#else : sfcurve2.set_y(i,log10(v))
		#sfcurve2.update()

	elif opt != None and opt.lower()=="auto":
		try:
			sfcurve2=XYData()		# this is really slow and stupid
			sfcurve2.read_file("strucfac.txt")

			sfcurve2.update()
			hasgoodsf=True
		except : sfcurve2=None

	# this empirical curve was generated by manually combining a groel solution scattering curve with an alpha-crystallin solution scattering curve, and 15 PDB models containing different folds
	cv=[2.80906, 1.97088, 0.71626, 0.44646, 0.17816, -0.03370, -0.27300, -0.43296, -0.47462, -0.49176, -0.51401, -0.50851, -0.50085, -0.51879, -0.50726, -0.44237, -0.46572, -0.41184, -0.37315, -0.36693, -0.38623, -0.36812, -0.38944, -0.44176, -0.49944, -0.59203, -0.67172, -0.70637, -0.75822, -0.82767, -0.80866, -0.79560, -0.79147, -0.77391, -0.75435, -0.74013, -0.71295, -0.67304, -0.63188, -0.59686, -0.56459, -0.53561, -0.52926, -0.51478, -0.52703, -0.54996, -0.56983, -0.59393, -0.61916, -0.64065, -0.65594, -0.66507, -0.67619, -0.69587, -0.72263, -0.74979, -0.77228, -0.79427, -0.81728, -0.84210, -0.86782, -0.88952, -0.90666, -0.92398, -0.93935, -0.95353, -0.96825, -0.98245, -0.99630, -1.00828, -1.01905, -1.02951,-1.04015, -1.04975, -1.05807, -1.06691, -1.07601, -1.08674, -1.09222, -1.09494, -1.09815, -1.10561, -1.11427, -1.11832, -1.11867, -1.11744, -1.12003, -1.12583, -1.13025, -1.13495, -1.13707, -1.13804, -1.14301, -1.14933, -1.14846, -1.14018, -1.12828, -1.11983, -1.12223]

	sfcurve=XYData()
	for i,j in enumerate(cv):
		sfcurve.set_x(i,i/200.0+.002)
		sfcurve.set_y(i,pow(10.0,cv[i]))

	if sfcurve2==None:
		print "No  structure factor found, using default internal structure factor. If fitting results are poor, consider rerunning --autofit once structure factor has been computed."
		sfcurve2=sfcurve
		hasgoodsf=False

def get_gui_arg_img_sets(filenames):
	'''
	returns the img_sets list required to intialized the GUI correctly. Each set has
	filename,EMAN2CTF,im_1d,bg_1d,im_2d,bg_2d,qual,bg_1d_low
	'''

	img_sets = []
	for fsp in filenames:
		try:
			js_parms=js_open_dict(info_name(fsp))
			img_set=js_parms["ctf"]
		except:
			print "Warning, you must run auto-fit before running the GUI. No parameters for ",info_name(fsp)
#			traceback.print_exc()
			continue
		try:
			img_set.append(js_parms["quality"])
		except:
			img_set.append(5)

		img_set.append(low_bg_curve(img_set[2],img_set[0].dsbg))
		try: img_set.append(js_parms["ctf_microbox"])
		except: pass
		img_sets.append([fsp]+img_set)

	return img_sets

def write_e2ctf_output(options):
	"write wiener filtered and/or phase flipped particle data to the local database"
	global logid

	if options.phaseflip or options.wiener or options.phasefliphp or options.storeparm:
		for i,filename in enumerate(options.filenames):
			name=base_name(filename)
			if debug: print "Processing ",filename

			try: im=EMData(filename,0,True)
			except:
				print "Error processing {}. Does not appear to be an image stack. Skipping.".format(filename)
				continue

			if options.phaseflip: phaseout="particles/{}__ctf_flip.hdf".format(name)
			else: phaseout=None

			if options.phasefliphp: phasehpout="particles/{}__ctf_flip_hp.hdf".format(name)
			else: phasehpout=None

			if options.wiener:
				if options.autohp: wienerout="particles/{}__ctf_wiener_hp.hdf".format(name)
				else: wienerout="particles/{}__ctf_wiener.hdf".format(name)
			else : wienerout=None

			phaseprocout=None
			if options.phaseflipproc!=None:
				phaseprocout=["particles/{}__ctf_flip_proc.hdf".format(name),parsemodopt(options.phaseflipproc)]

				if options.phaseflipproc2!=None:
					phaseprocout.append(parsemodopt(options.phaseflipproc2))

				if options.phaseflipproc3!=None:
					phaseprocout.append(parsemodopt(options.phaseflipproc3))

			try:
				js=js_open_dict(info_name(filename))
				ctf=js["ctf"][0]		# EMAN2CTF object from disk
				js.close()
			except:
				print "No CTF parameters found in {}, skipping {}.".format(info_name(filename),filename)
				continue
			if options.constbfactor>0: ctf.bfactor=options.constbfactor

			if phaseout : print "Phase image out: ",phaseout,"\t",
			if phaseprocout : print "Processed phase image out: ",phaseprocout[0],"\t",
			if phasehpout : print "Phase-hp image out: ",phasehpout,"\t",
			if wienerout : print "Wiener image out: ",wienerout,
			print "  defocus=",ctf.defocus

			process_stack(filename,phaseout,phasehpout,wienerout,phaseprocout,not options.nonorm,options.oversamp,ctf,invert=options.invert,storeparm=options.storeparm,source_image=options.source_image,zero_ok=options.zerook)

			if logid : E2progress(logid,float(i+1)/len(options.filenames))

def compute_envelope(img_sets,smax=.06):
		"""This computes the intensity of the background subtracted power spectrum around each CTF maximum for
		all data sets, then attempts to merge/normalize the results into a single aggregate curve. This curve
		should be proportional to the structure factor * average envelope function."""

#		global envelopes # needs to be a global for the Simplex minimizer
		# envelopes is essentially a cache of information that could be useful at later stages of refinement

		envelope=[]
		for i in img_sets:
			envelopes.append(ctf_env_points(i[2],i[3],i[1]))
#			envelope.extend(ctf_env_points(i[2],i[3],i[1]))

		# we use a simplex minimizer to try to rescale the individual sets to match as best they can
		scales=[1.0]*(len(img_sets)-1)
		if (len(img_sets)>3) :
			incr=[0.2]*len(img_sets)
			simp=Simplex(env_cmp,scales,incr,data=envelopes)
			scales=simp.minimize(maxiters=8000)[0]
		scales.insert(0,1.0)	# we keep this scale factor fixed

#		print envelopes
#		print scales

		# apply the final rescaling
		envelope=[]
		for i in range(len(scales)):
			cur=envelopes[i]
			for j in range(len(cur)):
				envelope.append([cur[j][0],cur[j][1]*scales[i]])

		envelope.sort()

		out=file("strucfac.allpoints.txt","w")
		for i in envelope: out.write("%f\t%f\n"%(i[0],i[1]))
		out.close()

		# this averages together all of the points at the same spatial frequency
		# not perfect, since different spatial frequencies may have contributions from
		# different curves
		sf=[]
		last=envelope[0]
		n=1
		for i in range(1,len(envelope)):
			if envelope[i][0]==last[0] :
				last[1]+=envelope[i][1]
				n+=1
			else :
				sf.append((last[0],last[1]/n))
				last=envelope[i]
				n=1
		envelope=sf

#		envelope=[i for i in envelope if i[1]>0]	# filter out all negative peak values
		# at smax we transition from the computed curve to the empirical curve, by default this is at ~16 A
		# we find this point, and also simultaneously write the data-based portion of the structure factor to disk
		out=file("strucfac.fromdata.txt","w")
		for i,j in enumerate(envelope):
			if j[0]>=smax :break
			out.write("{}\t{}\n".format(j[0],j[1]))

		sc=j[1]/sfact(j[0])
		print "\nTransitioning to internal structure factor at %1.1f A with scale factor %1.2f"%(1/j[0],sc)
		ds=envelope[i][0]-envelope[i-1][0]
		envelope=envelope[:i]
		s=j[0]
		while s<.293 :
			envelope.append((s,sfact(s)*sc))
			s+=ds
		envelope.append((1.0,sfact(.293)*sc))

		return envelope



def fixnegbg(bg_1d,im_1d,ds):
	"""This will make sure that we don't have a lot of negative values in the background
	subtracted curve near the CTF zeros. This would likely be due to the exluded volume
	of the particle from the solvent"""

	start=int(1.0/(25.0*ds))	# We don't worry about slight negatives below 1/40 A/pix
	start=max(1,start)

	# Find the worst negative peak
	ratio=1.0
	for i in range(start,len(bg_1d)/2):
		if (im_1d[i]+im_1d[i+1]+im_1d[i-1])/(bg_1d[i]+bg_1d[i+1]+bg_1d[i-1])<ratio :
			ratio=(im_1d[i]+im_1d[i+1]+im_1d[i-1])/(bg_1d[i]+bg_1d[i+1]+bg_1d[i-1])

	# return the corrected background
	print "BG correction ratio %1.4f"%ratio
	return [i*ratio for i in bg_1d]

def pspec_and_ctf_fit(options,debug=False):
	"""Power spectrum and CTF fitting. Returns an 'image sets' list. Each item in this list contains
	filename,EMAN2CTF,im_1d,bg_1d,im_2d,bg_2d,qual,bg_1d_low,micro_1d/None"""
	global logid
	img_sets=[]
	#db_parms=db_open_dict("bdb:e2ctf.parms")
	#db_im2d=db_open_dict("bdb:e2ctf.im2d")
	#db_bg2d=db_open_dict("bdb:e2ctf.bg2d")

	for i,filename in enumerate(options.filenames):
		name=base_name(filename)
		try : js_parms=js_open_dict(info_name(filename))
		except :
			print "ERROR: Cannot open {} for metadata storage. Exiting.".format(info_name(filename))
			sys.exit(1)

		# compute the power spectra
		if options.verbose or debug : print "Processing ",filename
		apix=options.apix
		if apix<=0 : apix=EMData(filename,0,1)["apix_x"]

		# After this, PS contains a list of (im_1d,bg_1d,im_2d,bg_2d,bg_1d_low) tuples. If classify is <2 then this list will have only 1 tuple in it
		if options.classify>1 : ps=split_powspec_with_bg(filename,options.source_image,radius=options.bgmask,edgenorm=not options.nonorm,oversamp=options.oversamp,apix=apix,nclasses=options.classify,zero_ok=options.zerook)
		else: ps=list((powspec_with_bg(filename,options.source_image,radius=options.bgmask,edgenorm=not options.nonorm,oversamp=options.oversamp,apix=apix,zero_ok=options.zerook,wholeimage=options.wholeimage,highdensity=options.highdensity),))
		# im_1d,bg_1d,im_2d,bg_2d,bg_1d_low,micro_1d/none
		if ps==None :
			print "Error fitting CTF on ",filename
			continue
		try: ds=1.0/(apix*ps[0][2].get_ysize())
		except:
			print "Error fitting CTF (ds) on ",filename
			continue
		
		for j,p in enumerate(ps):
			try: im_1d,bg_1d,im_2d,bg_2d,bg_1d_low,micro_1d=p
			except:
				im_1d,bg_1d,im_2d,bg_2d,bg_1d_low=p
				micro_1d=None
			if not options.nosmooth : bg_1d=smooth_bg(bg_1d,ds)
			if options.fixnegbg :
				bg_1d=fixnegbg(bg_1d,im_1d,ds)		# This insures that we don't have unreasonable negative values

			if debug: Util.save_data(0,ds,bg_1d,"ctf.bgb4.txt")

			# Fit the CTF parameters
			if debug : print "Fit CTF"
			if options.curdefocushint or options.curdefocusfix:
				try:
					ctf=js_parms["ctf"][0]
					curdf=ctf.defocus
					curdfdiff=ctf.dfdiff
					curdfang=ctf.dfang
					if options.curdefocushint: dfhint=(curdf-0.1,curdf+0.1)
					else: dfhint=(curdf-.001,curdf+.001)
					print "Using existing defocus as hint :",dfhint
				except :
					try:
						ctf=js_parms["ctf_frame"][1]
						curdf=ctf.defocus
						curdfdiff=ctf.dfdiff
						curdfang=ctf.dfang
						if options.curdefocushint: dfhint=(curdf-0.1,curdf+0.1)
						else: dfhint=(curdf-.001,curdf+.001)
						print "Using existing defocus from frame as hint :",dfhint
					except:
						dfhint=None
						print "No existing defocus to start with"
			else: dfhint=(options.defocusmin,options.defocusmax)
			ctf=ctf_fit(im_1d,bg_1d,bg_1d_low,im_2d,bg_2d,options.voltage,options.cs,options.ac,apix,bgadj=not options.nosmooth,autohp=options.autohp,dfhint=dfhint,highdensity=options.highdensity,verbose=options.verbose)
			if options.astigmatism and not options.curdefocusfix : ctf_fit_stig(im_2d,bg_2d,ctf,verbose=1)
			elif options.astigmatism:
				ctf.dfdiff=curdfdiff
				ctf.dfang=curdfang

			im_1d,bg_1d=calc_1dfrom2d(ctf,im_2d,bg_2d)
			if options.constbfactor>0 : ctf.bfactor=options.constbfactor
			else: ctf.bfactor=ctf_fit_bfactor(list(array(im_1d)-array(bg_1d)),ds,ctf)


			if debug:
				Util.save_data(0,ds,im_1d,"ctf.fg.txt")
				Util.save_data(0,ds,bg_1d,"ctf.bg.txt")
				Util.save_data(0,ds,ctf.snr,"ctf.snr.txt")

			try : qual=js_parms["quality"]
			except :
				qual=5
				js_parms["quality"]=5
			if j==0: img_sets.append([filename,ctf,im_1d,bg_1d,im_2d,bg_2d,qual,bg_1d_low,micro_1d])
			else: img_sets.append([filename+"_"+str(j),ctf,im_1d,bg_1d,im_2d,bg_2d,qual,bg_1d_low,micro_1d])

		# store the results back in the database. We omit the filename, quality and bg_1d_low (which can be easily recomputed)
		if img_sets[-1][-1]==None: js_parms.delete("ctf_microbox")
		else: js_parms["ctf_microbox"]=img_sets[-1][-1]
		js_parms["ctf"]=img_sets[-1][1:-3]
		js_parms.close()

		if logid : E2progress(logid,float(i+1)/len(options.filenames))

	project_db = js_open_dict("info/project.json")
	try: project_db.update({ "global.microscope_voltage":options.voltage, "global.microscope_cs":options.cs, "global.apix":apix })
	except:
		print "ERROR: apix not found. This probably means that no CTF curves were sucessfully fit !"

	return img_sets

def refine_and_smoothsnr(options,strfact,debug=False):
	"""This will refine already determined defocus values by maximizing high-resolution smoothed
	SNR and also replace the stored SNR curve with the smoothed version. Note that this REQUIRES
	that a structure factor be loaded and present"""
	global logid
	img_sets=[]
	#db_parms=db_open_dict("bdb:e2ctf.parms")
	#db_parms_old=db_open_dict("bdb:e2ctf.parms.presmooth")

	olddf=[]
	newdf=[]
	skipped=0
	for i,filename in enumerate(options.filenames):
		name=base_name(filename)
		js_parms=js_open_dict(info_name(filename))

		if debug : print "Processing ",filename

		try:
			orig=js_parms["ctf"]
			ctf=orig[0]
		except:
			print "Error: no fit data for {}. Skipping.".format(name)
			skipped+=1
			continue

		if ctf.dfdiff!=0 :
			print "Skipping {}. SSNR based refinement not available for astigmatic images.".format(name)
			skipped+=1
			continue

		olddf.append(ctf.defocus)
		im_1d=orig[1]
		bg_1d=orig[2]

		ds=ctf.dsbg
		r=len(ctf.background)
		s=[ds*ii for ii in range(r)]

		# Restore SNR to original unfiltered version
		ctf.snr=[snr_safe(im_1d[i],bg_1d[i]) for i in range(len(im_1d))]

		# Tune the defocus to maximize high res snr
		if debug : print "Fit Defocus"
		best=(0,olddf[-1])
		for df in [olddf[-1]+ddf/1000.0 for ddf in xrange(-100,101)]:
			ctf.defocus=df
			ssnr=ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_SNR_SMOOTH,strfact)		# The smoothed curve
#			print len( ssnr),len(s)
#			ssnr=[a*b for a,b in enumerate(ssnr)]	# this would impose a r weighting to make high res agreement more important
			qual=sum(ssnr[int(.08/ds):len(s)-2])
			best=max(best,(qual,df))

		newdf.append(best[1])
		if options.verbose : print "%s: %1.4f -> %1.4f"%(filename,olddf[-1],best[1])

		ctf.defocus=best[1]
		ctf.snr=ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_SNR_SMOOTH,strfact)
		js_parms["ctf"]=[ctf]+orig[1:]

		if logid : E2progress(logid,float(i+1)/len(options.filenames))

	if skipped>0 : print "Warning: %d files skipped"%skipped
	if len(olddf)>0 : print "Mean defocus adjustment : %1.4f um"%((sum([fabs(olddf[ii]-newdf[ii]) for ii in xrange(len(olddf))]))/len(olddf))


def env_cmp(sca,envelopes):
#	global envelopes
	env=envelopes
	total=[(j[0],j[1]) for j in env[0]]		# the first set is unscaled to 'anchor' the results

	# env contains a list of lists of points. Each outermost list represents a data set
	# sca is an array of scale factors. We set the first one to 1.0 to avoid arbitrary data scaling
	for i,ii in enumerate(env[1:]):
		for j in ii:
			total.append((j[0],j[1]*sca[i]))

	total.sort()

	ret=0
	for i in range(2,len(total)-2):
		if total[i][1] :
			ret+=((total[i-2][1]/total[i][1])**2+(total[i-1][1]/total[i][1])**2+(total[i+1][1]/total[i][1])**2+(total[i+2][1]/total[i][1])**2)
#			ret+=fabs(total[i-2][1]-total[i][1])+fabs(total[i-1][1]-total[i][1])+fabs(total[i+1][1]-total[i][1])+fabs(total[i+2][1]-total[i][1])

	for i in sca:
		if i<0 : ret*=1.5

	#ret=0
	#for i in range(1,len(total)):
		#if total[i][1] :
			#ret+=fabs((total[i-1][1]/total[i][1])-1.0)/(total[i][0]-total[i-1][0]+.0005)

	return ret

def process_stack(stackfile,phaseflip=None,phasehp=None,wiener=None,phaseproc=None,edgenorm=True,oversamp=1,default_ctf=None,invert=False,storeparm=False,source_image=None,zero_ok=False):
	"""Will phase-flip and/or Wiener filter particles in a file based on their stored CTF parameters.
	phaseflip should be the path for writing the phase-flipped particles
	wiener should be the path for writing the Wiener filtered (and possibly phase-flipped) particles
	oversamp will oversample as part of the processing, ostensibly permitting phase-flipping on a wider range of defocus values
	"""

	im=EMData(stackfile,0)
	ys=im.get_ysize()*oversamp
	ys2=im.get_ysize()
	n=EMUtil.get_image_count(stackfile)
	lctf=None
	#db_parms=db_open_dict("bdb:e2ctf.parms")

	#if virtualout:
		#vin=db_open_dict(stackfile)
		#vout=db_open_dict(virtualout)

	name=base_name(stackfile)
	try :
		js_parms=js_open_dict(info_name(stackfile))
		ctf=js_parms["ctf"]
	except :
		print "ERROR: Cannot find CTF parameters in {}. Skipping.".format(info_name(filename))
		return

	if phasehp:
		p1d=js_parms["ctf"][1]
		for c in xrange(2,len(p1d)):
			if p1d[c-1]<p1d[c] : break
		c-=1
		trg=p1d[c-1]
		hpfilt=[1.0 for i in range(int(ys*1.5))]
		hpfilt[0]=0.0
		for i in xrange(1,c):
			try :hpfilt[i]=sqrt(trg/p1d[i])
			except: hpfilt[i]=0.0

		oscor=2.0*len(p1d)/ys
		if oscor!=floor(oscor) : print "Warning, incompatible oversampling from earlier results %d vs %d"%(len(p1d),ys/2)

		oscor=int(oscor)
#		print hpfilt[:c+4]
		hpfilt=[0.0]+[sum(hpfilt[i+1:i+1+oscor])/oscor for i in range(0,len(hpfilt)-1,oscor)]	# this downsamples the correction curve

#		print hpfilt[:c+4]

	js_parms.close()

	for i in range(n):
		if source_image!=None :
			im1=EMData()
			im1.read_image(stackfile,i,True)
			try:
				if im1["ptcl_source_image"]!=source_image : continue
			except:
				print "Image %d doesn't have the ptcl_source_image parameter. Skipping."%i
				continue

		try: im1 = EMData(stackfile,i)
		except: im1.to_zero()
		im1.process_inplace("mask.zeroedgefill",{"nonzero":1})		# This tries to deal with particles that were boxed off the edge of the micrograph

		# If we detected a zero edge, we mark the particle as bad
		if not zero_ok and im1.has_attr("hadzeroedge") and im1["hadzeroedge"]!=0:
			print "Particle outside of micrograph detected, marking as bad ({},{})".format(stackfile,i)
			js=js_open_dict(info_name(stackfile))
			try:
				s=js["sets"]
			except:
				js["sets"]={"bad_particles":[i]}
			else:
				try: s["bad_particles"].append(i)
				except : s["bad_particles"]=[i]
				js["sets"]=s

		try: ctf=im1["ctf"]
		except : ctf=default_ctf
		if storeparm :
			if stackfile[-4:].lower()!=".hdf" and stackfile[:4].lower()!="bdb:" :
				if i==0: print "Warning, --storeparm option ignored. Input file must be HDF or BDB for this option to work."
			else :
				ctf=default_ctf		# otherwise we're stuck with the values in the file forever
				im1["ctf"]=ctf
				im1.write_image(stackfile,i,EMUtil.ImageType.IMAGE_UNKNOWN,True)
		if type(ctf)==EMAN1Ctf : ctf=default_ctf	# EMAN1 ctf needs a structure factor for this to work


		if edgenorm : im1.process_inplace("normalize.edgemean")
		if oversamp>1 :
			im1.clip_inplace(Region(-(ys2*(oversamp-1)/2),-(ys2*(oversamp-1)/2),ys,ys))
#			print -(ys2*(oversamp-1)/2),-(ys2*(oversamp-1)/2),ys,ys
#		print i
		fft1=im1.do_fft()

		if phaseflip or phasehp or phaseproc:
			if not lctf or not lctf.equal(ctf):
				flipim=fft1.copy()
				ctf.compute_2d_complex(flipim,Ctf.CtfType.CTF_SIGN)
#				if i==0: flipim.write_image("flip.mrc")
			fft1.mult(flipim)
			out=fft1.do_ift()
			out["ctf"]=ctf
			out["apix_x"] = ctf.apix
			out["apix_y"] = ctf.apix
			out["apix_z"] = ctf.apix
			out.clip_inplace(Region(int(ys2*(oversamp-1)/2.0),int(ys2*(oversamp-1)/2.0),ys2,ys2))
			if invert: out.mult(-1.0)
			if edgenorm: out.process("normalize.edgemean")
			if phaseflip: out.write_image(phaseflip,i)

			if phaseproc!=None:
				out2=out.copy()				# processor may or may not be in Fourier space
				out2["ctf"]=ctf
				out2["apix_x"] = ctf.apix
				out2["apix_y"] = ctf.apix
				out2["apix_z"] = ctf.apix
				# we take a sequence of processor option 2-tuples
				for op in phaseproc[1:]:
					out2.process_inplace(op[0],op[1])
#				out2.clip_inplace(Region(int(ys2*(oversamp-1)/2.0),int(ys2*(oversamp-1)/2.0),ys2,ys2))

#				print fft2.get_ysize(),len(hpfilt)

				if edgenorm: out2.process_inplace("normalize.edgemean")
				out2.write_image(phaseproc[0],i)

			if phasehp:
				fft2=fft1.copy()
				fft2.process_inplace("filter.radialtable",{"table":hpfilt})
				out=fft2.do_ift()
				out["ctf"]=ctf
				out["apix_x"] = ctf.apix
				out["apix_y"] = ctf.apix
				out["apix_z"] = ctf.apix
				out.clip_inplace(Region(int(ys2*(oversamp-1)/2.0),int(ys2*(oversamp-1)/2.0),ys2,ys2))

#				print fft2.get_ysize(),len(hpfilt)

				if edgenorm: out.process_inplace("normalize.edgemean")
				if invert: out.mult(-1.0)
				#process_inplace("filter.highpass.autopeak")
				out.write_image(phasehp,i)



		if wiener :
			if not lctf or not lctf.equal(ctf):
				wienerim=fft1.copy()
				ctf.compute_2d_complex(wienerim,Ctf.CtfType.CTF_WIENER_FILTER)
#				print wienerim.get_attr_dict()
#				display(wienerim)
#				print ctf.to_string()
#				plot(ctf.background)
#				plot(ctf.snr)
#				plot(ctf.compute_1d(ys,Ctf.CtfType.CTF_WIENER_FILTER))
			fft1.mult(wienerim)
			out=fft1.do_ift()
			out["ctf"]=ctf
			out["apix_x"] = ctf.apix
			out["apix_y"] = ctf.apix
			out["apix_z"] = ctf.apix
			out.clip_inplace(Region(int(ys2*(oversamp-1)/2.0),int(ys2*(oversamp-1)/2.0),ys2,ys2))
			if invert : out.mult(-1.0)
			out.process("normalize.edgemean")
			try: out.write_image(wiener,i)
			except:
				print wiener,i
				try: out.write_image(wiener,i)
				except:
					print "!!! ",wiener,i
					out.write_image("error.hed",-1)

		#if virtualout:
			#im1["data_path"]=vin.get_data_path(i)
			#vout[vout["maxrec"]+1]=im1

		lctf=ctf

	#if wiener and wiener[:4]=="bdb:" : db_close_dict(wiener)
	#if phaseflip and phaseflip[:4]=="bdb:" : db_close_dict(phaseflip)

	db_close_dict(stackfile)	# this is safe even for non bdb: files

	return

def powspec(stackfile,source_image=None,mask=None,edgenorm=True):
	"""This routine will read the images from the specified file, optionally edgenormalize,
	optionally apply a mask then compute the average
	2-D power spectrum for the stack. If source_image is provided, then it
	will only read particles with ptcl_source_image set to the specified value.

	Results returned as a 2-D FFT intensity/0 image"""

	n=EMUtil.get_image_count(stackfile)

	for i in range(n):
		if source_image!=None :
			im1=EMData()
			im1.read_image(stackfile,i,True)
			try:
				if im1["ptcl_source_image"]!=source_image : continue
			except:
				print "Image %d doesn't have the ptcl_source_image parameter. Skipping."%i
				continue

		im=EMData(stackfile,i)
		if edgenorm : im.process_inplace("normalize.edgemean")
		if mask : im*=mask
		imf=im.do_fft()
		imf.ri2inten()
		if i==0: av=imf
		else: av+=imf

	av/=(float(n)*av.get_ysize()*av.get_ysize())
	av.set_value_at(0,0,0.0)
#	av.process_inplace("xform.fourierorigin.tocenter")

	av.set_complex(1)
	av.set_attr("is_intensity", 1)

	#db_close_dict(stackfile)		# safe for non bdb urls
	return av

masks={}		# mask cache for background/foreground masking
def powspec_with_bg(stackfile,source_image=None,radius=0,edgenorm=True,oversamp=1,apix=2,ptclns=None,zero_ok=False,wholeimage=False,highdensity=False):
	"""This routine will read the images from the specified file, optionally edgenormalize,
	then apply a gaussian mask with the specified radius then compute the average 2-D power
	spectrum for the stack. It will also compute the average 2-D power spectrum using 1-mask + edge
	apotization to get an appoximate 'background' power spectrum. 2-D results returned as a 2-D FFT
	intensity/0 image. 1-D results returned as a list of floats. If source_image is provided, then it
	will only read particles with ptcl_source_image set to the specified value.

	At present, highdensity does nothing here, that is implemented in ctf_fit

	returns a 5-tuple with spectra for (1d particle,1d background,2d particle,2d background,1d background non-convex,1d foreground from wholeimage or None)
	"""

	global masks

	im=EMData(stackfile,0)
	ys=im.get_ysize()*oversamp
	ys2=im.get_ysize()
	if radius<=0 : radius=ys2/2.6
	n=EMUtil.get_image_count(stackfile)
	nn=0
	ds=1.0/(apix*ys)	# oversampled ds

	# set up the inner and outer Gaussian masks
	try:
		mask1,ratio1,mask2,ratio2=masks[(ys,radius)]
	except:
		# Mask 1 is an "inner" Gaussian to extract primarily the particle from the middle of the image
		mask1=EMData(ys2,ys2,1)
		mask1.to_one()
		mask1.process_inplace("mask.gaussian",{"outer_radius":radius,"exponent":4.0})
		# Mask 2 is the 'inverse' (1.0-val) of mask1, with the addition of a soft outer edge to reduce periodic boundary condition issues
		mask2=mask1.copy()*-1+1
#		mask1.process_inplace("mask.decayedge2d",{"width":4})
		mask2.process_inplace("mask.decayedge2d",{"width":4})
		mask1.clip_inplace(Region(-(ys2*(oversamp-1)/2),-(ys2*(oversamp-1)/2),ys,ys))
		mask2.clip_inplace(Region(-(ys2*(oversamp-1)/2),-(ys2*(oversamp-1)/2),ys,ys))

		# ratio1,2 give us info about how much of the image the mask covers for normalization purposes
		ratio1=mask1.get_attr("square_sum")/(ys*ys)	#/1.035
		ratio2=mask2.get_attr("square_sum")/(ys*ys)
		masks[(ys,radius)]=(mask1,ratio1,mask2,ratio2)
#		display((mask1,mask2))

	av1,av2=None,None
	for i in range(n):
		if ptclns!=None and i not in ptclns :continue
		im1 = EMData()
		if source_image!=None :
			im1.read_image(stackfile,i,True)
			try:
				if im1["ptcl_source_image"]!=source_image : continue
			except:
				print "Image %d doesn't have the ptcl_source_image parameter. Skipping."%i
				continue

		im1.read_image(stackfile,i)

		# Images with flat edges due to boxing too close to the edge can adversely impact the power spectrum
		if not zero_ok :
			im1.process_inplace("mask.zeroedgefill",{"nonzero":1})		# This tries to deal with particles that were boxed off the edge of the micrograph
			if im1.has_attr("hadzeroedge") and im1["hadzeroedge"]!=0:
				print "Skipped particle with bad edge ({}:{})".format(stackfile,i)
				continue

		nn+=1
#		im1=EMData(stackfile,i)

		if edgenorm : im1.process_inplace("normalize.edgemean")
		if oversamp>1 :
#			print Region(-(ys2*(oversamp-1)/2),-(ys2*(oversamp-1)/2),ys,ys)
			im1.clip_inplace(Region(-(ys2*(oversamp-1)/2),-(ys2*(oversamp-1)/2),ys,ys))

		im2=im1.copy()

#		print im2.get_size(), im1.get_size()

		# now we compute power spectra for the 2 regions defined by the masks
		# av1/2 contain the incoherent averaged power spectra (intensity average)
		im1*=mask1
		imf=im1.do_fft()
		imf.ri2inten()
		if av1==None:
			av1=EMData(imf["nx"],imf["ny"],imf["nz"])	# we make a new object to avoid copying the header of imf
			av1.set_complex(True)
			av1.to_zero()
		av1+=imf

		im2*=mask2
		imf=im2.do_fft()
		imf.ri2inten()
		if av2==None:
			av2=EMData(imf["nx"],imf["ny"],imf["nz"])
			av2.set_complex(True)
			av2.to_zero()
		av2+=imf

	if nn==0 : return None

	# normalize the 2 curves
	av1/=(float(nn)*av1.get_ysize()*av1.get_ysize()*ratio1)
	av1.set_value_at(0,0,0.0)
	av1.set_complex(1)
	av1["is_intensity"]=1
	av1["ptcl_repr"]=nn

	av2/=(float(nn)*av2.get_ysize()*av2.get_ysize()*ratio2)
	av2.set_value_at(0,0,0.0)
	av2.set_complex(1)
	av2["is_intensity"]=1
	av2["ptcl_repr"]=nn

	# These now should represent the FG+BG and BG curves
	av1_1d=av1.calc_radial_dist(av1.get_ysize()/2,0.0,1.0,1)
	av2_1d=av2.calc_radial_dist(av2.get_ysize()/2,0.0,1.0,1)

	# in this mode we box out "particles" from the original micrograph to get the smoothest curve for defocus fitting
	# must not be use for (for example) structure factor
	av3_1d=None
	if wholeimage:
		try: micro=EMData("micrographs/{}.hdf".format(base_name(stackfile)))
		except:
			print "Error: --wholeimage specified, but could not find ","micrographs/{}.hdf".format(base_name(stackfile))
			sys.exit(1)

		bs=av1["ny"]
		av3=av1.copy()
		av3.to_zero()
		nrg=0
		# overlapping regions
		for y in xrange(bs/2,micro["ny"]-bs,bs/2):
			for x in xrange(bs/2,micro["nx"]-bs,bs/2):
				im1=micro.get_clip(Region(x,y,bs,bs))
				im1.process_inplace("normalize.edgemean")
				im1*=mask1
				imf=im1.do_fft()
				imf.ri2inten()
				av3+=imf
				nrg+=1

		print nrg
		av3/=(float(nrg)*av1.get_ysize()*av1.get_ysize()*ratio1)
		av3.set_value_at(0,0,0.0)
		av3.set_complex(1)
		av3["is_intensity"]=1
		av3_1d=av3.calc_radial_dist(av3["ny"]/2,0.0,1.0,1)


	# added to make the 2D display look better. Should have no other impact at the time
	# it's being added, though autofitting may rely on it in future, so it shouldn't be removed --steve (8/3/11)
	av1.process_inplace("math.sqrt")
	av1["is_intensity"]=0

	av2.process_inplace("math.sqrt")
	av2["is_intensity"]=0

	# This is a new addition (2/4/10) to prevent negative BG subtracted curves near the origin
	maxpix=int(0.04/ds)               # we do this up to ~25 A
	avsnr=0
	avc=0
	for i in xrange(maxpix):
		if av1_1d[i]>av2_1d[i] :
			avsnr+=(av1_1d[i]-av2_1d[i])/av2_1d[i]
			avc+=1

	if avc==0 :
		print "Failed to readjust background in {}. Returning what I can..."
		return (av1_1d,av2_1d,av1,av2,low_bg_curve(av2_1d,ds),av3_1d)

	avsnr/=avc

	for i in xrange(maxpix) :
		if av2_1d[i]>av1_1d[i] : av2_1d[i]=av1_1d[i]/(avsnr+1)          # we set the background equal to the foreground if it's too big

	#db_close_dict(stackfile)	# safe for non-bdb files

	return (av1_1d,av2_1d,av1,av2,low_bg_curve(av2_1d,ds),av3_1d)

# img_sets : filename,ctf,im_1d,bg_1d,im_2d,bg_2d,qual,bg_1d_low
def split_powspec_with_bg(stackfile,source_image=None,radius=0,edgenorm=True,oversamp=1,apix=2,nclasses=2,zero_ok=False):
	"""Same as powspec_with_bg, but also subclassifies the data into groups based on differences in the 1-D background subtracted power spectrum.
Rather than returning a single tuple, returns a list of nclasses tuples.
	"""

	global masks

	im=EMData(stackfile,0)
	ys=im.get_ysize()*oversamp
	ys2=im.get_ysize()
	if radius<=0 : radius=ys2/2.6
	n=EMUtil.get_image_count(stackfile)
	nn=0
	ds=1.0/(apix*ys)	# oversampled ds

	# set up the inner and outer Gaussian masks
	try:
		mask1,ratio1,mask2,ratio2=masks[(ys,radius)]
	except:
		mask1=EMData(ys2,ys2,1)
		mask1.to_one()
		mask1.process_inplace("mask.gaussian",{"outer_radius":radius,"exponent":4.0})
		mask2=mask1.copy()*-1+1
#		mask1.process_inplace("mask.decayedge2d",{"width":4})
		mask2.process_inplace("mask.decayedge2d",{"width":4})
		mask1.clip_inplace(Region(-(ys2*(oversamp-1)/2),-(ys2*(oversamp-1)/2),ys,ys))
		mask2.clip_inplace(Region(-(ys2*(oversamp-1)/2),-(ys2*(oversamp-1)/2),ys,ys))
		ratio1=mask1.get_attr("square_sum")/(ys*ys)	#/1.035
		ratio2=mask2.get_attr("square_sum")/(ys*ys)
		masks[(ys,radius)]=(mask1,ratio1,mask2,ratio2)
#		display((mask1,mask2))

	av_1d_n=[]
	for i in range(n):
		im1 = EMData()
		if source_image!=None :
			im1.read_image(stackfile,i,True)
			try:
				if im1["ptcl_source_image"]!=source_image : continue
			except:
				print "Image %d doesn't have the ptcl_source_image parameter. Skipping."%i
				continue

		im1.read_image(stackfile,i)
		nn+=1
#		im1=EMData(stackfile,i)

		if edgenorm : im1.process_inplace("normalize.edgemean")
		if oversamp>1 :
#			print Region(-(ys2*(oversamp-1)/2),-(ys2*(oversamp-1)/2),ys,ys)
			im1.clip_inplace(Region(-(ys2*(oversamp-1)/2),-(ys2*(oversamp-1)/2),ys,ys))

		im2=im1.copy()

#		print im2.get_size(), im1.get_size()

		im1*=mask1
		imf=im1.do_fft()
		imf.ri2inten()
		if i==0:
			av1=EMData(imf["nx"],imf["ny"],imf["nz"])	# we make a new object to avoid copying the header of imf
			av1.set_complex(True)
			av1.to_zero()

		av_1d=imf.calc_radial_dist(av1.get_ysize()/2,0.0,1.0,1)
		av_1d_n.append(av_1d)
		av1+=imf

		im2*=mask2
		imf=im2.do_fft()
		imf.ri2inten()
		if i==0:
			av2=EMData(imf["nx"],imf["ny"],imf["nz"])
			av2.set_complex(True)
			av2.to_zero()
		av2+=imf


	av1/=(float(nn)*av1.get_ysize()*av1.get_ysize()*ratio1)
	av1.set_value_at(0,0,0.0)
	av1.set_complex(1)
	av1["is_intensity"]=1
	av1["ptcl_repr"]=nn

	av2/=(float(nn)*av2.get_ysize()*av2.get_ysize()*ratio2)
	av2.set_value_at(0,0,0.0)
	av2.set_complex(1)
	av2["is_intensity"]=1
	av2["ptcl_repr"]=nn

	av1_1d=av1.calc_radial_dist(av1.get_ysize()/2,0.0,1.0,1)
	av2_1d=av2.calc_radial_dist(av2.get_ysize()/2,0.0,1.0,1)

	n0=int(0.05/ds)		# N at 20 A. We ignore low resolution information due to interference of structure factor
	bg=EMData(len(av2_1d)-n0,1,1)
	for i in xrange(n0,len(av2_1d)): bg[i-n0]=av2_1d[i]

	nvec=10
	nxl=len(av2_1d)-n0
	mask=EMData(nxl,1,1)
	mask.to_one()
#	pca=Analyzers.get("pca_large",{"mask":mask,"nvec":nvec})

	# now BG subtract each particle 1-D average
	for j in range(n):
		# converts each list of radial values into an image
		im=EMData(nxl,1,1)
		for i in xrange(n0,len(av2_1d)): im[i-n0]=av_1d_n[j][i]/(av1["ny"]*av1["ny"]*ratio1)
#		im.write_image("presub.hdf",j)
		im.sub(bg)
		sm=0
		for k in xrange(nxl*3/4,nxl): sm+=im[k]
		sm/=nxl-nxl*3/4
#		print sm
		im.sub(sm)
#		im.write_image("postsub.hdf",j)
		av_1d_n[j]=im
#		pca.insert_image(im)

#	result=pca.analyze()		# this gives us a basis for decomposing and classifying the images
#	for j in range(nvec): result[j].process_inplace("normalize.unitlen")

	# Instead of using a MSA data-based subspace, let's use the spectrum of different CTF curves instead
	ctf=EMAN2Ctf()
	ctf.bfactor=50.0		# These parameters do not need to be accurate, they are just used to make a nice set of overlapping CTF shaped curves
	ctf.ampcont=10.0
	ctf.voltage=200.0
	ctf.cs=2.0
	ctf.apix=apix
	ctf.dsbg=ds
	dfs=[5.0,4.2,3.5,2.8,2.1,1.7,1.4,1.1,0.9,0.6]		# a nice set of overlapping defocus values for classification
	nvec=len(dfs)
	result=[]
	for i in xrange(nvec):
		ctf.defocus=dfs[i]
		vec=ctf.compute_1d(ys*2,ds,Ctf.CtfType.CTF_AMP,None)  # note that size is not the returned size, but 2*the returned size
		result.append(EMData(nxl,1,1))
		for j in xrange(n0,len(av2_1d)): result[-1][j-n0]=vec[j]
		result[-1].write_image("cvec.hdf",-1)

	# project into the subspace
	prj_1d=[]
	for j in range(n):
		im=EMData(nvec,1,1)
		for k in range(nvec): im[k]=av_1d_n[j].cmp("ccc",result[k])
		im.write_image("postsub.proj.hdf",j)
		prj_1d.append(im)

	#out=file("projs.txt","w")
	#for y in range(n):
		#for x in range(nvec):
			#out.write("{}\t".format(prj_1d[y][x]))
		         #out.write("\n")

	# kmeans classification
	km=Analyzers.get("kmeans",{"ncls":nclasses,"mininclass":5})
	for im in prj_1d: km.insert_image(im)
	cls=km.analyze()

	plists={}
	for i,im in enumerate(prj_1d):
		try: plists[im["class_id"]].append(i)
		except: plists[im["class_id"]]=[i]

#	print plists

	return [powspec_with_bg(stackfile,source_image,radius,edgenorm,oversamp,apix,plists[p]) for p in plists]

	sys.exit(1)
	# added to make the 2D display look better. Should have no other impact at the time
	# it's being added, though autofitting may rely on it in future, so it shouldn't be removed --steve (8/3/11)
	av1.process_inplace("math.sqrt")
	av1["is_intensity"]=0

	av2.process_inplace("math.sqrt")
	av2["is_intensity"]=0

	# This is a new addition (2/4/10) to prevent negative BG subtracted curves near the origin
	maxpix=int(apix*ys2/25.0)               # we do this up to ~80 A
	avsnr=0
	avc=0
	for i in xrange(maxpix):
		if av1_1d[i]>av2_1d[i] :
			avsnr+=(av1_1d[i]-av2_1d[i])/av2_1d[i]
			avc+=1
	avsnr/=avc

	for i in xrange(maxpix) :
		if av2_1d[i]>av1_1d[i] : av2_1d[i]=av1_1d[i]/(avsnr+1)          # we set the background equal to the foreground if it's too big

	#db_close_dict(stackfile)	# safe for non-bdb files

	return (av1_1d,av2_1d,av1,av2,low_bg_curve(av2_1d,ds))


def bgedge2d(stackfile,width):
	"""This routine will read the images from the specified file, and compute the average
	2-D power spectrum computed using boxes taken from the edge of the image. Returns the
	1-D power spectrum as a list of floats. This is not presently used in e2ctf since it
	produces a heavily downsampled background curve, and is provided only for experimentation."""

	n=EMUtil.get_image_count(stackfile)
	av=None

	for i in range(n):
		im=EMData(stackfile,i)

		xs=im.get_xsize()		# x size of image
		xst=int(floor(xs/ceil(xs/width)))	# step to use so we cover xs with width sized blocks

		# Build a list of all boxes around the edge
		boxl=[]
		for x in range(0,xs-xst/2,xst):
			boxl.append((x,0))
			boxl.append((x,xs-xst))
		for y in range(xst,xs-3*xst/2,xst):
			boxl.append((0,y))
			boxl.append((xs-xst,y))

		for b in boxl:
			r=im.get_clip(Region(b[0],b[1],width,width))
			imf=r.do_fft()
			imf.ri2inten()
			if av : av+=imf
			else: av=imf

	av/=(n*len(boxl)*width*width)
	av.set_value_at(0,0,0.0)

	av.set_complex(1)
	av.set_attr("is_intensity", 1)

	#db_close_dict(stackfile)
	return av

def smooth_bg(curve,ds):
	"""Smooths a background curve by doing a running average of the log of the curve, ignoring the first few points"""

	first=int(.02/ds)	# start at 1/50 1/A
	if first<2 : first=2

	return curve[:first]+[pow(curve[i-1]*curve[i]*curve[i+1],.33333) for i in range(first,len(curve)-2)]+[curve[-2],curve[-1]]
#	return curve[:first]+[pow(curve[i-2]*curve[i-1]*curve[i]*curve[i+1]*curve[i+2],.2) for i in range(first,len(curve)-2)]+[curve[-2],curve[-1]]

### Note that the following function has been replaced by a C++ function in Util
def smooth_by_ctf(curve,ds,ctf):
	"""Smooths a curve based on an (already fit) CTF object, starting with the first zero. Assumes that locally the value should be the sum of
	an offset and the ctf * a constant. This does make the curve have the general appearance of the CTF, locally.
	Takes the curve, with the ds value for the curve and a CTF object, returns a new curve."""

	curvea=array(curve)
	bf=ctf.bfactor
	ctf.bfactor=50
	ccurv=array(ctf.compute_1d(len(curve)*2,ds,Ctf.CtfType.CTF_AMP))**2
	ctf.bfactor=bf
#	z=int(zero(1,ctf.voltage,ctf.cs,ctf.defocus,ctf.ampcont)/ds+.5)				# location of first zero in terms of n
	z=ctf.zero(0)
	z=max(4,z)

	ret=curve[:z]
	for i in xrange(z,len(curve)-4):
		lc =curvea[i-4:i+5]			# data to be smoothed
		lcf=array(ccurv[i-4:i+5])	# reference CTF curve
		lcf-=lcf.mean()				# so we can do a dot product
		lcf/=(lcf**2).sum()			# unit length
		ret.append(lc.mean()+lcf[4]*lcf.dot(lc))
		print i*ds,lcf.dot(lc)

	ret+=curve[-4:]
	Util.save_data(0,ds,curve,"a.txt")
	Util.save_data(0,ds,ret,"b.txt")
	Util.save_data(0,ds,list(ccurv),"c.txt")
	return ret

def ctf_fit_bfactor(curve,ds,ctf):
	"""Determines a decent B-factor for the BG subtracted 'curve' (returns the new value). 'ctf' must represent a pretty accurate fit (other than B-factor)"""

	bf=ctf.bfactor
	ctf.bfactor=50
	ccurv=ctf.compute_1d(len(curve)*2,ds,Ctf.CtfType.CTF_AMP)
	ccurv=[f*f for f in ccurv]
	ctf.bfactor=bf

	# We compute a windowed normalized running correlation coefficient between the data and the fit CTF
	# window size depends on curve length
	if len(curve)<64 : wdw=4
	elif len(curve)<128 : wdw=6
	elif len(curve)<256 : wdw=8
	else : wdw=10
	sim=Util.windowdot(curve,ccurv,wdw,1)

	#for i in xrange(len(curve)):
		#print i*ds,a[i],curve[i],ccurv[i]

	risethr=0.75*max(sim[int(0.04/ds):])
	# find the last point where the curve rises above 0.75 its max value
	for i in xrange(len(curve)-wdw,int(0.04/ds),-1):
		if sim[i]>risethr : break

	# now find the first place where it falls below 0.1
	for i in xrange(i,len(curve)-wdw):
		if sim[i]<0.1 : break

	maxres=1.0/(i*ds)
	print "maxres ",maxres," B -> ",maxres*maxres*6.0

	return maxres*maxres*6.0


def least_square(data):
	"simple linear regression for y=mx+b on a list of (x,y) points. Use the C routine if you need speed."
	sum,sum_x,sum_y,sum_xx,sum_xy=0,0,0,0,0
	for d in data:
		y=d[1]

		sum_x+=d[0]
		sum_xx+=d[0]*d[0]
		sum_y+=y
		sum_xy+=d[0]*y
		sum+=1.0

	denom=sum*sum_xx-sum_x*sum_x
	if denom==0 : denom=.00001

	m=(sum*sum_xy-sum_x*sum_y)/denom
	b=(sum_xx*sum_y-sum_x*sum_xy)/denom

	return(m,b)

def snr_safe(s,n) :
	if s<=0 or n<=0 or s-n<0 : return 0.001		# this used to return 0, but I think it may have been causing some snr weighting problems
	return (s-n)/n

def sfact2(s):
	"""This will return a structure factor for the current protein. If no specific structure factor was available or
	specified, it will default to the same curve as sfact()"""

	global sfcurve2
	if sfcurve2==None : init_sfcurve("auto")

	return sfcurve2.get_yatx(s,1)

def sfact(s):
	"""This will return a curve shaped something like the structure factor of a typical protein. It is not designed to be
	highly accurate, but be good enough for approximate B-factor estimation"""

	global sfcurve
	if sfcurve==None : init_sfcurve("auto")

	return sfcurve.get_yatx(s,1)

	#if s<.004 : return 0
##	if s<.006 : return 0
##	if s>.2934 : s=.2934
	#if s>.21 : s=.21

	## this curve should be pretty valid on the range 0.004 - 0.2934, we limit it a bit more to prevent distractions from the sharp peak
	#return pow(10.0,3.6717 - 364.58 * s + 15597 * s**2 - 4.0678e+05 * s**3 + 6.7098e+06 * s**4 - 7.0735e+07 * s**5 + 4.7839e+08 * s**6 - 2.0574e+09 * s**7 +5.4288e+09 * s**8 - 8.0065e+09 * s**9 + 5.0518e+09 * s**10)

def low_bg_curve(bg_1d,ds):
	"Computes a 1-D curve that follows the minima of a background function as best possible, focusing on the lower resolution range"

	# This code has been replaced by an equivalent C++ function
	#ret=bg_1d[:]

	## force the curve to be non-convex, after a point
	#cont=True
	#while cont:
		#cont=False
		#for i in range(3,len(ret)-1):
			#qt = (ret[i-1]+ret[i+1])/2.0
			#if ret[i]>qt :
				#ret[i] = qt
				#cont = True

	#return ret

	return Util.nonconvex(bg_1d,3);

def elambda(V):
	"""returns relativistic electron wavelength. V in KV. Wavelength in A"""
	return 12.3/sqrt(1000*V+0.97845*V*V)

# This isn't really right, a correct version is now implemented in the EMAN2Ctf object
#def zero(N,V,Cs,Z,AC):
	#"""Return the spatial frequency of the order N zero of the CTF with Voltage V, Cs in mm, Defocus Z (positive underfocus) and amplitude contrast AC (0-100)"""
	#acshift=-acos(sqrt(1-AC*AC/10000.0))	# phase shift in radians due to amplitude contrast
	#gamma=pi*N+acshift
	#l=elambda(V)
##	print acshift,gamma,l,Cs*gamma*l*l*l+5.0*l*l*pi*Z*Z
	#return sqrt(-Z/(1000.0*Cs*l*l)+sqrt(Cs*gamma*l*l*l+5.0*l*l*pi*Z*Z)/(3963.327*Cs*l*l*l))

def calc_1dfrom2d(ctf,fg2d,bg2d):
	"""Computes adjusted 1-D power spectra from 2-D FG and BG data. Needs ctf object, since 1-D average uses astigmatism compensation.
returns (fg1d,bg1d)"""
	ds=ctf.dsbg
#	r=len(ctf.background)
#	s=[ds*i for i in range(r)]

	fg=array(ctf.compute_1d_fromimage(len(ctf.background)*2, ds, fg2d))
	bg=array(ctf.compute_1d_fromimage(len(ctf.background)*2, ds, bg2d))
	bgc=bg.copy()
#			bg=smooth_by_ctf(bg,ds,ctf)					# still need to work on this routine, similar concept to SNR smoothing, but has some problems
#	bglow=low_bg_curve(bg,ds)

	#fg2d.write_image("a.hdf")
	#bg2d.write_image("b.hdf")
	#Util.save_data(0,ds,list(fg),"av1.txt")
	#Util.save_data(0,ds,list(bg),"av2.txt")


	# here we adjust the background to make the zeroes zero
	n=2
#	lz=int(zero(1,ctf.voltage,ctf.cs,ctf.defocus,ctf.ampcont)/ds+.5)
	lz=int(ctf.zero(0)/ds+.5)
	while 1:
#		z=int(zero(n,ctf.voltage,ctf.cs,ctf.defocus,ctf.ampcont)/ds+.5)
		z=int(ctf.zero(n-1)/ds+.5)
		if z>len(ctf.background)-2 or z-lz<2: break
		d1=min(fg[lz-1:lz+2]-bg[lz-1:lz+2])
		d2=min(fg[ z-1: z+2]-bg[ z-1: z+2])
		#d1=min(fg[lz]-bg[lz],fg[lz-1]-bg[lz-1],fg[lz+1]-bg[lz+1])
		#d2=min(fg[z]-bg[z],fg[z-1]-bg[z-1],fg[z+1]-bg[z+1])
#				print lz,d1,z,d2
		for x in xrange(lz,z):
			bg[x]+=(z-x)/float(z-lz)*d1+(x-lz)/float(z-lz)*d2
		if n==2 : lwd,lwz=d1,lz
		lz=z
		n+=1

	# deal with the points from the origin to the first zero
	try :
		bg[:lwz]+=lwd
		bg[:lwz]=numpy.minimum(bg[:lwz],fg[:lwz])	# this makes sure the bg subtracted curve isn't negative before the first zero
	except:
		print "ERROR in flattening background. This should only occur if the defocus is either too close or too far from focus with the current box-size and sampling."
		return (list(fg),list(bg))


	# deal with the points from where the zeroes got too close together, just to make a smooth curve
	for x in xrange(lz,len(ctf.background)) :
		bg[x]+=fg[x-2:x+3].mean()-bgc[x-2:x+3].mean()

	return (list(fg),list(bg))


def ctf_fit_stig(im_2d,bg_2d,ctf,verbose=1):
	"""Refines the astigmatism parameters given a good initial fit. Modifies CTF object in-place !!!  No return value."""

	## we start with some astigmatism or the process may not converge well
	#if ctf.dfdiff==0 : ctf.dfdiff=0.1

	bgsub=im_2d-bg_2d
	bgcp=bgsub.copy()
	#sim=Simplex(ctf_stig_cmp,[ctf.dfdiff,ctf.dfang],[0.01,15.0],data=(bgsub,bgcp,ctf))
	#oparm=sim.minimize(epsilon=.00000001,monitor=0)
	#print oparm

	# Give a little arbitrary astigmatism for the angular search
	if ctf.dfdiff==0 : ctf.dfdiff=ctf.defocus/20.0

	oldb=ctf.bfactor
	# coarse angular alignment
	besta=(1.0e15,0)
	ctf.bfactor=500
	for ang in xrange(0,180,15):
		v=ctf_stig_cmp((ctf.defocus+ctf.dfdiff/2.0,ctf.defocus-ctf.dfdiff/2.0,ang),(bgsub,bgcp,ctf))
		besta=min(besta,(v,ang))
	ctf.dfang=besta[1]
	print "best angle:", besta

	# Use a simplex minimizer to find the final fit
	# we minimize using defocusU and defocusV rather than defocus & dfdfiff
	ctf.bfactor=200
	sim=Simplex(ctf_stig_cmp,[ctf.defocus+ctf.dfdiff/2.0,ctf.defocus-ctf.dfdiff/2.0,ctf.dfang],[0.01,0.01,5.0],data=(bgsub,bgcp,ctf))
	oparm=sim.minimize(epsilon=.00000001,monitor=0)
	dfmaj,dfmin,ctf.dfang=oparm[0]		# final fit result
	print "Coarse refine: defocus={:1.4f} dfdiff={:1.5f} dfang={:3.2f} defocusU={:1.4f} defocusV={:1.4f}".format((dfmaj+dfmin)/2.0,(dfmaj-dfmin),ctf.dfang,dfmaj,dfmin)

	# Use a simplex minimizer to refine the local neighborhood
	ctf.bfactor=80
	sim=Simplex(ctf_stig_cmp,oparm[0],[0.005,2.0,.005],data=(bgsub,bgcp,ctf))
	oparm=sim.minimize(epsilon=.00000001,monitor=0)
	dfmaj,dfmin,ctf.dfang=oparm[0]		# final fit result
	print "  Fine refine: defocus={:1.4f} dfdiff={:1.5f} dfang={:3.2f} defocusU={:1.4f} defocusV={:1.4f}".format((dfmaj+dfmin)/2.0,(dfmaj-dfmin),ctf.dfang,dfmaj,dfmin)

	ctf.bfactor=oldb
	ctf.defocus=(dfmaj+dfmin)/2.0
	ctf.dfdiff=dfmaj-dfmin


	## extract points at maxima for B-factor estimation
	#ds=ctf.dsbg
	#fg,bg=calc_1dfrom2d(ctf,im_2d,bg_2d)
	#bglow=low_bg_curve(bg,ds)

	#n=2
	#lz=int(zero(1,ctf.voltage,ctf.cs,ctf.defocus,ctf.ampcont)/ds+.5)
	#while 1:
		#z=int(zero(n,ctf.voltage,ctf.cs,ctf.defocus,ctf.ampcont)/ds+.5)
		#if z>len(ctf.background)-2 or n>30: break
		#m=(z+lz)/2
		#print m*ds,fg[m]-bg[m]
		#lz=z
		#n+=1

	#print oparm

	## Print out smoothed fit values
	#smooth_by_ctf(list(array(fg)-array(bg)),ds,ctf)

	#ctf.bfactor=200
	## Fit dfdiff and defocus simultaneously by exhaustive search
	#outim=EMData(21,21)
	#bestd=(1.0e15,0,0)
	#dfcen=ctf.defocus
	#dfstep=ctf.defocus/100.0
	#dfdstep=ctf.dfdiff/10.0
	#for defocus in xrange(-10,11):
		#for dfdiff in xrange(0,21):
			#v=ctf_stig_cmp((dfdiff*dfstep,ctf.dfang,dfcen+defocus*dfstep),(bgsub,bgcp,ctf))
			#bestd=min(bestd,(v,dfdiff*dfstep,dfcen+defocus*dfstep))
			#outim[defocus+10,dfdiff]=v
	#print "best dfdiff/defocus ",bestd
	#ctf.dfdiff=bestd[1]
	#ctf.defocus=bestd[2]
	#outim.update()
	#outim.write_image("tst.hdf")

	## fine angular alignment
	#for ang in xrange(besta[1]-14,besta[1]+15):
		#v=ctf_stig_cmp((ctf.dfdiff,ang,ctf.defocus),(bgsub,bgcp,ctf))
		#besta=min(besta,(v,ang))
	#ctf.dfang=besta[1]
	#print "best angle:", besta
	#ctf.bfactor=oldb

def ctf_stig_cmp(parms,data):
	"""energy function for fitting astigmatism, parms is (dfmaj, dfmin, dfang) instead of using the internal defocus/dfdiff"""

	bgsub,bgcp,ctf=data

	dfmaj,dfmin,ctf.dfang=parms
	ctf.defocus=(dfmaj+dfmin)/2.0
	ctf.dfdiff=dfmaj-dfmin

	ctf.compute_2d_complex(bgcp,Ctf.CtfType.CTF_FITREF,None)

	#bgcp.write_image("a.hdf",-1)
	#bgsub.write_image("b.hdf",0)

	penalty=0.0
	if ctf.dfdiff<0 : penalty-=ctf.dfdiff
	#if ctf.dfdiff>ctf.defocus : penalty+=ctf.dfdiff-ctf.defocus
	if dfmaj<0 : penalty-=dfmaj
	if dfmin<0 : penalty-=dfmin

#	print parms,ctf.defocus,ctf.dfdiff,bgcp.cmp("dot",bgsub,{"normalize":1}),penalty
	return bgcp.cmp("dot",bgsub,{"normalize":1})+penalty

def ctf_fit(im_1d,bg_1d,bg_1d_low,im_2d,bg_2d,voltage,cs,ac,apix,bgadj=0,autohp=False,dfhint=None,highdensity=False,verbose=1):
	"""Determines CTF parameters given power spectra produced by powspec_with_bg()
	The bgadj option will result in adjusting the bg_1d curve to better match the zeroes
	of the CTF (in which case bg_1d is modified in place)."""
	# defocus estimation
	global debug
	global sfcurve
#	from bisect import insort_left

	if dfhint==None : dfhint = (0.5,5.5,0.05)
	elif isinstance(dfhint,float) : dfhint=(dfhint-.2, dfhint+.2,0.02)
	else: dfhint=(dfhint[0],dfhint[1],min(0.05,(dfhint[1]-dfhint[0])/5))

	ys=im_2d.get_ysize()
	ds=1.0/(apix*ys)
	if ac<0 or ac>100 :
		print "Invalid %%AC, defaulting to 10"
		ac=10.0

	curve=[im_1d[i]-bg_1d[i] for i in xrange(len(im_1d))]

	ctf=EMAN2Ctf()
	ctf.from_dict({"defocus":1.0,"voltage":voltage,"bfactor":150.0,"cs":cs,"ampcont":ac,"apix":apix,"dsbg":ds,"background":bg_1d,"dfdiff":0,"defang":0})

	if len(curve)<64 : wdw=4
	elif len(curve)<128 : wdw=6
	elif len(curve)<256 : wdw=8
	else : wdw=10

	best=(-1,1.0)

	# we might get better results if we readjusted the bg correction after each try, but it's reaaaaly slow
	# so we compute an unadjusted bg one time, then do one local correction once the defocus is close
	im=im_1d
	bg=ctf.compute_1d_fromimage(len(ctf.background)*2, ds, bg_2d)

	for rng in (0,1):
		# second pass is +-0.1 unless the original hint range was narrower
		if rng==1: dfhint=(max(dfhint[0],ctf.defocus-0.1),min(dfhint[1],ctf.defocus+0.1),min(dfhint[2]/2.0,0.005))

		curve=[im[i]-bg[i] for i in xrange(len(im_1d))]
		for df in arange(dfhint[0],dfhint[1],dfhint[2]):
			ctf.defocus=df
			ccurv=ctf.compute_1d(len(curve)*2,ds,Ctf.CtfType.CTF_AMP)
			ccurv=[sfact2(ds*i)*ccurv[i]**2 for i in range(len(ccurv))]		# squared * structure factor

			## Recompute the background assuming the defocus is correct
			#im,bg=calc_1dfrom2d(ctf,im_2d,bg_2d)
			#curve=[im[i]-bg[i] for i in xrange(len(im_1d))]

			sim=array(Util.windowdot(curve,ccurv,wdw,1))
			qual=sim[int(ctf.zero(0)/ds):int(ctf.zero(5)/ds)].mean()
			if qual>best[0]: best=(qual,df)
#			print df,sum(sim),qual

		ctf.defocus=best[1]
		print "Best defocus: {:1.03f}".format(best[1])

		# determine a good B-factor now that the defocus is pretty good
		ctf.bfactor=ctf_fit_bfactor(curve,ds,ctf)

		im,bg=calc_1dfrom2d(ctf,im_2d,bg_2d)

	if bgadj :
		for i in xrange(len(bg)): bg_1d[i]=bg[i]		# overwrite the input background with our final adjusted curve
		bglow=low_bg_curve(bg,ds)
		for i in xrange(len(bg)): bg_1d_low[i]=bglow[i]
	ctf.background=bg_1d

	ctf.snr=[snr_safe(im[i],bg[i]) for i in range(len(im_1d))]

	return ctf

#def ctf_fit(im_1d,bg_1d,bg_1d_low,im_2d,bg_2d,voltage,cs,ac,apix,bgadj=0,autohp=False,dfhint=None,verbose=1):
	#"""Determines CTF parameters given power spectra produced by powspec_with_bg()
	#The bgadj option will result in adjusting the bg_1d curve to better match the zeroes
	#of the CTF (in which case bg_1d is modified in place)."""
	## defocus estimation
	#global debug
	#global sfcurve
##	from bisect import insort_left

	#ys=im_2d.get_ysize()
	#ds=1.0/(apix*ys)
	#if ac<0 or ac>100 :
		#print "Invalid %%AC, defaulting to 10"
		#ac=10.0

	#ctf=EMAN2Ctf()
	#ctf.from_dict({"defocus":1.0,"voltage":voltage,"bfactor":500.0,"cs":cs,"ampcont":ac,"apix":apix,"dsbg":ds,"background":bg_1d})

	#sf = [sfact(i*ds) for i in range(ys)]

	#bgsub=[im_1d[s]-bg_1d[s] for s in range(len(im_1d))]	# background subtracted curve, good until we readjust the background
	#bglowsub=[im_1d[s]-bg_1d_low[s] for s in range(len(im_1d))]	# background subtracted curve, using a non-convex version of the background

	#s1=min(int(.167/ds),ys/3-4)

##	plot(bglowsub)
	## We find the first minimum (also <1/50 max) after the value has fallen to 1/5 max
	#if isinstance(dfhint,float) :
##		s0=max(2,int(zero(1,voltage,cs,dfhint,ac)/ds-3))
		#s0=max(2,int(ctf.zero(0)/ds-3))
		#while (bglowsub[s0]>bglowsub[s0+1] or bglowsub[s0]>bglowsub[s0-1]) and s0<s1: s0+=1	# look for a minimum in the data curve
	#elif isinstance(dfhint,tuple) :
##		s0=max(2,int(zero(1,voltage,cs,(dfhint[0]+dfhint[1])/2.0,ac)/ds-3))
		#s0=max(2,int(ctf.zero(0)/ds-3))
		#while (bglowsub[s0]>bglowsub[s0+1] or bglowsub[s0]>bglowsub[s0-1]) and s0<s1: s0+=1	# look for a minimum in the data curve
	#else :
		#maxsub=max(bglowsub[int(.01/ds):])
		#for i in range(len(bglowsub)-1,0,-1):
			#if bglowsub[i]>maxsub/5.0 : break
		#s0=max(int(.01/ds),i)
		#while (bglowsub[s0]>bglowsub[s0+1] or bglowsub[s0+1]>maxsub/50.0) and s0<s1: s0+=1	# look for a minimum in the data curve
		#if s0==s1 : s0=int(.03/ds)
	#if verbose: print "Minimum at 1/%1.1f 1/A (%1.4f), highest s considered 1/%1.1f 1/A (%1.4f)"%(1.0/(s0*ds),s0*ds,1.0/(s1*ds),s1*ds)

	## we now have a range from the first minimia to ~6 A. Adjust it to include a bit of data closer to the origin
	#for i in range(s0,s0*2/3,-1):
		#if bglowsub[i]>bglowsub[i-1] : break
	#s0m=s0									# save for a defocus estimate
	#s0=i

	#if verbose: print "s0 adjusted to 1/%1.1f 1/A"%(1.0/(s0*ds))

	#if s1<=s0 :
		#print "Error: couldn't fit this data set due to inability to locate appropriate minima"
		#return ctf	# This way we can still proceed to the GUI... a bit of a bad hack...

	#if debug:
		#dfout=file("ctf.df.txt","w")

##	dfbest1=(0,-1.0e20)
	#dfbest1=[]

	## This implies that we already have a general range for the defocus, so we limit the search
	#if isinstance(dfhint,float) :
		#dfhint=int(dfhint*20.0)
		#rng=range(dfhint-3,dfhint+4)
		#if verbose: print "Using defocus hint %1.2f (%1.2f - %1.2f)"%(dfhint/20.0,rng[0]/20.0,rng[-1]/20.0)
	#elif isinstance(dfhint,tuple) :
		#rng=(int(dfhint[0]*20.0),int(dfhint[1]*20.0))
		#if verbose: print "Using limited defocus range (%1.2f - %1.2f)"%(rng[0]/20.0,rng[-1]/20.0)
	#else:
		#rng=range(3,128)
		#if verbose: print "(%s) Using full defocus range (%1.2f - %1.2f)"%(str(dfhint),rng[0]/20.0,rng[-1]/20.0)


	## This loop tries to find the best few possible defocuses
	#for dfi in rng:			# loop over defocus
		#df=dfi/20.0
		#ctf.defocus=df
##		ctf.ampcont=ac
		#cc=ctf.compute_1d(ys,ds,Ctf.CtfType.CTF_AMP)
		#norm=0
		#for fz in range(len(cc)):
			#if cc[fz]<0 : break

		#tot,tota,totb=0,0,0
		#zroa,zrob=0,0
		#for s in range(s0,s1):
			## This is basicaly a dot product of the reference curve vs the simulation
			#a=(cc[s]**2)				# ctf curve ^2
			#b=bglowsub[s]					# bg subtracted intensity
			#tot+=a*b
			#tota+=a*a
			#totb+=b*b

			## This computes the mean value near zeroes
			##if a<.1 and s<ys/3:
				##zroa+=(im_1d[s]-bg_1d[s])*s
				##zrob+=1.0

		#tot/=sqrt(tota*totb)	# correct if we want a normalized dot product
##		tot/=tota				# funnny 'normalization' which seems to help bring out the more likely peaks
##		if zrob==0 : continue
##		tot*=-(zroa/zrob)

		#dfbest1.append((tot,df))		# we keep all of the results, then process them at the end
##		if tot>dfbest1[1] : dfbest1=(df,tot)
		#try :dfout.write("%1.2f\t%g\n"%(df,tot))
		#except : pass

	## Now we find the best few defocus choices
	#dfbest1a=[]						# keep only peaks, and sort them
	#for i in range(1,len(dfbest1)-1):
		#if dfbest1[i]>dfbest1[i-1] and dfbest1[i]>dfbest1[i+1] : dfbest1a.append(dfbest1[i])		# keep only local peaks

	#if len(dfbest1a)==0 : dfbest1a=[max(dfbest1)]

	#dfbest1a.sort()
	#dfbest1a=dfbest1a[-5:]		# keep at most the best 5 peaks

	#if verbose: print "Initial defocus possibilities: ",
	#for i in dfbest1a: print i[1],
	#print

	#if len(dfbest1a)==1 :
		#best=[[0.0,[dfbest1a[0][1],300.0]]]
	#else :
		## Next, we use a simplex minimizer to try for the best CTF parameters for each defocus
		#best=[]
		#for b1a in dfbest1a:
			## our parameter set is (defocus,bfactor)
			#parm=[b1a[1],500.0]

##			print "Initial guess : ",parm
			#sim=Simplex(ctf_cmp_2,parm,[.02,20.0],data=(ctf,bglowsub,s0,s1,ds,parm[0],rng))
			#oparm=sim.minimize(epsilon=.0000001,monitor=0)
##			print "Optimized : ",oparm
			#best.append((oparm[1],oparm[0]))


		## best now contains quality,(df,bfac) for each optimized answer
		#best.sort()
##		if verbose: print "Best value is df=%1.3f  B=%1.1f"%(best[0][1][0],best[0][1][1])

	#ctf.defocus=best[0][1][0]
	#ctf.bfactor=best[0][1][1]
##	cc=ctf.compute_1d(ys,ds,Ctf.CtfType.CTF_AMP)
##	Util.save_data(0,ds,cc,"ctf.ctf.txt")
##	print best[0]

	#if bgadj:
		#bg2=bg_1d[:]
		#if verbose: print "BG adjustment using df=%1.3f  B=%1.1f"%(best[0][1][0],best[0][1][1])
		#for i in range(6):
			## now we try to construct a better background based on the CTF zeroes being zero
			#df=best[0][1][0]
			#ctf.defocus=best[0][1][0]
			#ctf.bfactor=best[0][1][1]
			#cc=ctf.compute_1d(ys,ds,Ctf.CtfType.CTF_AMP)
##			bg2=bg_1d[:]
			#last=0,1.0
			#for x in range(1,len(bg2)-1) :
				#if cc[x]*cc[x+1]<0 :
					#if x<len(bg2)/2 :
						## Below 1/2 nyquist, we work harder to avoid negative values
						#cur=(x,min(im_1d[x]/bg2[x],im_1d[x+1]/bg2[x+1],im_1d[x-1]/bg2[x-1],im_1d[x+2]/bg2[x+2]))
					#else:
						## we search the two points 'at' the zero for the minimum
						## if we're too aggressive about this, we will end up exaggerating the high
						## resolution SNR
						#cur=(x,min(im_1d[x]/bg2[x],im_1d[x+1]/bg2[x+1]))

					## once we have a pair of zeros we adjust the background values between
					#for xx in range(last[0],cur[0]):
						#w=(xx-last[0])/float(cur[0]-last[0])
						#bg_1d[xx]=bg2[xx]*(cur[1]*w+last[1]*(1.0-w))

					#last=cur
			## cover the area from the last zero crossing to the end of the curve
			#for xx in range(last[0],len(bg2)):
				#bg_1d[xx]=bg2[xx]*last[1]

	##	s0=int(.04/ds)+1
	##	s1=min(int(0.15/ds),len(bg_1d)-1)

			## rerun the simplex with the new background
			#bgsub=[im_1d[s]-bg_1d[s] for s in range(len(im_1d))]
			#bg_1d_low=low_bg_curve(bg_1d,ds)
			#bglowsub=[im_1d[s]-bg_1d_low[s] for s in range(len(im_1d))]	# background subtracted curve, using a non-convex version of the background
			#best[0][1][1]=500.0		# restart the fit with B=200.0
			#try:
				#if hasgoodsf: sim=Simplex(ctf_cmp_a,best[0][1],[.02,20.0],data=(ctf,bgsub,s0,s1,ds,best[0][1][0],rng))
				#else: sim=Simplex(ctf_cmp,best[0][1],[.02,20.0],data=(ctf,bgsub,s0,s1,ds,best[0][1][0],rng))
				#oparm=sim.minimize(epsilon=.00000001,monitor=0)
				#if fabs(df-oparm[0][0])/oparm[0][0]<.001:
					#best[0]=(oparm[1],oparm[0])
					#break
				#best[0]=(oparm[1],oparm[0])
				#if verbose: print "After BG correction, value is df=%1.3f  B=%1.1f"%(best[0][1][0],best[0][1][1])
			#except:
				#print "Serious fitting error. Retaining bad fit parameters so manual fitting is possible. This usually indicates a data problem."
				#break
	#else:
		#if verbose: print "Best value is df=%1.3f  B=%1.1f"%(best[0][1][0],best[0][1][1])
		## rerun the simplex with the new background
		#best[0][1][1]=500.0		# restart the fit with B=200.0
		#if hasgoodsf : sim=Simplex(ctf_cmp_a,best[0][1],[.02,20.0],data=(ctf,bgsub,s0,s1,ds,best[0][1][0],rng))
		#else : sim=Simplex(ctf_cmp,best[0][1],[.02,20.0],data=(ctf,bgsub,s0,s1,ds,best[0][1][0],rng))
		#oparm=sim.minimize(epsilon=.0000001,monitor=0)
		#best[0]=(oparm[1],oparm[0])
		#if verbose: print "After BG correction, value is df=%1.3f  B=%1.1f"%(best[0][1][0],best[0][1][1])

	#snr=[snr_safe(im_1d[i],bg_1d[i]) for i in range(len(im_1d))]

	## This will dramatically reduce the intensity of the initial sharp peak found in almost all single particle data
	## this applies to the SNR curve only, downweighting the importance of this section of the spectrum without actually
	## removing the information by filtering the image data. It will, of course also impact Wiener filters.
	#if autohp:
		#for x in range(2,len(snr)-2):
			#if snr[x]>snr[x+1] and snr[x+1]<snr[x+2] : break	# we find the first minimum

		#snr1max=max(snr[1:x])				# find the intensity of the first peak
		#snr2max=max(snr[x+2:len(snr)/2])		# find the next highest snr peak

		#for xx in range(1,x+1): snr[xx]*=0.5*snr2max/snr1max		# scale the initial peak to 50% of the next highest peak


	## store the final results
	#ctf.snr=snr
	#ctf.defocus=best[0][1][0]
	#ctf.bfactor=best[0][1][1]

	#if verbose : print "Best DF = %1.3f   B-factor = %1.0f   avSNR = %1.4f"%(ctf.defocus,ctf.bfactor,(sum(ctf.snr)/float(len(ctf.snr))))


##	Util.save_data(0,ds,ctf.snr,"ctf.snr")
##	Util.save_data(0,ds,bg_1d,"ctf.bg1d")


	#return ctf

def ctf_cmp(parms,data):
	"""This function is a quality metric for a set of CTF parameters vs. data"""
	ctf,bgsub,s0,s1,ds,dforig,rng=data
	ctf.defocus=parms[0]
	ctf.bfactor=parms[1]
	cc=ctf.compute_1d(len(bgsub)*2,ds,Ctf.CtfType.CTF_AMP)	# this is for the error calculation
	ctf.defocus=parms[0]
#	ctf.bfactor=400.0
	ctf.bfactor=parms[1]
#	c2=ctf.compute_1d(len(bgsub)*2,ds,Ctf.CtfType.CTF_AMP)	# this is for the error calculation

	# now compute the error
	# This is an interesting similarity metric. It divides one curve by the other. In theory, this should be flat, with the mean
	# value equal to the 'amplitude' of the curve. However, of course, near zeroes there are singularities, so this tends to focus
	# on making these singularities as symmetric as possible. Other variances are
	# also common. So, what we compute is the standard deviation of this value about the mean, trying to make it as flat as possible
	global sfcurve
	if sfcurve.get_size()!=99 : s0/=2
	a,b,c=0,0,0
	mx=0
#	tp=[]
#	tp2=[]
	for i in range(s0,s1):
		v=sfact2(i*ds)*cc[i]**2
#		tp.append(v)

		if cc[i]**2>.001:
			a+=bgsub[i]/v
#			tp2.append(bgsub[i]/v)
			b+=(bgsub[i]/v)**2
			c+=1
			mx=i
#		else : tp2.append(0)
#	plot(tp2,tp,bgsub)

	mean=a/c
	sig=b/c-mean*mean

	er=sig/fabs(mean)
#	print er

#	er*=(1.0+300.0*(parms[0]-dforig)**4)		# This is a weight which biases the defocus towards the initial value
	er*=1.0+(parms[1]-200)/20000.0+exp(-(parms[1]-50.0)/30.0)		# This is a bias towards small B-factors and to prevent negative B-factors
	er*=max(1.0,1.0+parms[0]*20.0-rng[1])**2		# penalty for being outside range (high)
	er*=max(1.0,1.0+rng[0]-parms[0]*20.0)**2		# penalty for being outside range (low)

	#out=file("dbg","a")
	#out.write("%f\t%f\t%f\t%f\t%f\t%d\n"%(parms[0],parms[1],er,mean,sig,mx))

	return er

def ctf_cmp_a(parms,data):
	"""This function is a quality metric for a set of CTF parameters vs. data. It is a replacement for ctf_cmp"""
	ctf,bgsub,s0,s1,ds,dforig,rng=data
	ctf.defocus=parms[0]
	ctf.bfactor=parms[1]
	cc=ctf.compute_1d(len(bgsub)*2,ds,Ctf.CtfType.CTF_AMP)	# this is for the error calculation
	s0=s0*3/2			# This should be roughlythe position of the first zero
	s0a=max(3,s0/4)		# This is a lower s0 which will include some of the low resolution structure factor region

	# make a complete CTF curve over 0,s1 range
	cc=[sfact2(i*ds)*cc[i]**2 for i in range(0,s1)]
	bs=[max(0,f) for f in bgsub[0:s1]]				# clip negative values

	# normalize (adjust amplitude)
	s1=sum(cc[s0a:])
	s2=sum(bs[s0a:])
	cc=[f*s2/s1 for f in cc]

#	plot(cc,bs)

	er=0.0
#	for i,f in enumerate(cc): er+=fabs(f-bs[i])
	for i in range(s0,len(cc)):
		if i>=s0a and i<s0 : er+=fabs(cc[i]-bs[i])
		if i>=s0 : er+=(cc[i]-bs[i])**2.0

	er*=max(1.0,1.0+parms[0]*20.0-rng[-1])**2		# penalty for being outside range (high)
	er*=max(1.0,1.0+rng[0]-parms[0]*20.0)**2		# penalty for being outside range (low)


	return er


def ctf_cmp_2(parms,data):
	"""This function is a quality metric for a set of CTF parameters vs. data
	This version is used for rough fitting of the defocus. ctf_cmp is used to fine-tune the values. """
	ctf,bgsub,s0,s1,ds,dforig,rng=data
	ctf.defocus=parms[0]
	ctf.bfactor=0.0
	#c2=ctf.compute_1d(len(bgsub)*2,ds,Ctf.CtfType.CTF_AMP)	# we use this for amplitude scaling
	#ctf.bfactor=parms[1]
	cc=ctf.compute_1d(len(bgsub)*2,ds,Ctf.CtfType.CTF_AMP)	# this is for the error calculation

	## compute an optimal amplitude for these parameters
	#a,b=0,0
	#for i in range(s0,s1):
		#if fabs(c2[i])>.8 and bgsub[i]>0 :
			#a+=sfact(i*ds)*cc[i]*cc[i]
			#b+=bgsub[i]
	#if b>0 : norm=a/b
	#else :
		#norm=1.0

#	Util.save_data(0,ds,cc,"a.txt")
#	Util.save_data(0,ds,bgsub,"b.txt")

	# now compute the error
	global sfcurve
#	if sfcurve.get_size()!=99 : s0/=2
	#er=0
	#for i in range(s0,s1):
		#v=sfact(i*ds)*cc[i]*cc[i]*parms[2]
		#er+=(v-bgsub[i])**2

#	out=file("xxx","w")
	a,b,c=0,0,0
	for i in range(s0,s1):
		v=sfact(i*ds)*cc[i]*cc[i]		# structure factor * CTF intensity
		a+=v*bgsub[i]					#
		b+=v*v
		c+=bgsub[i]*bgsub[i]
#		if v>.001 : out.write("%f\t%f\n"%(i*ds,bgsub[i]/v))
#	print s0,s1,a,b,c
	er=1.0-a/sqrt(b*c)

	er1=er

#	print er,(parms[0]-dforig)**2,parms[1]

	er*=(1.0+300.0*(parms[0]-dforig)**4)		# This is a weight which biases the defocus towards the initial value
	er*=max(1.0,1.0+parms[0]*20.0-rng[-1])**2		# penalty for being outside range (high)
	er*=max(1.0,1.0+rng[0]-parms[0]*20.0)**2		# penalty for being outside range (low)
	er*=1.0+(parms[1]-200)/20000.0+exp(-(parms[1]-50.0)/30.0)		# This is a bias towards small B-factors and to prevent negative B-factors
#	er*=(1.0+fabs(parms[1]-200.0)/100000.0)		# This is a bias towards small B-factors
#	print "%1.3g\t%1.3g\t%1.3g\t%1.3g\t"%(er,er1,er/er1,parms[0]),(1.0+300.0*(parms[0]-dforig)**4),max(1.0,1.0+parms[0]*20.0-rng[-1])**2,max(1.0,1.0+rng[0]-parms[0]*20.0)**2,1.0+(parms[1]-200)/20000.0+exp(-(parms[1]-50.0)/30.0)

	#out=file("dbg","a")
	#out.write("%f\t%f\t%f\n"%(parms[0],sqrt(parms[1]),er))

	return er

def ctf_env_points(im_1d,bg_1d,ctf) :
	"""This will return a list of x,y points corresponding to the CTF corrected power spectrum near the maxima"""
	ys=len(bg_1d)
	ds=ctf.dsbg
	cc=ctf.compute_1d(ys,ds,Ctf.CtfType.CTF_AMP)
	ret=[]

	lo=1
	for i in range(1,len(cc)-1):
		if (lo or fabs(cc[i])>0.5) and im_1d[i]-bg_1d[i]>0 :
			ret.append((i*ds,(im_1d[i]-bg_1d[i])/cc[i]**2))
#			ret.append((i*ds,(im_1d[i]-bg_1d[i])/sfact(i*ds)))		# this version removes the structure factor (in theory)
		if lo and fabs(cc[i])>0.5 : lo=0							# this gets the low frequencies before the first maximum

	return ret

#def ctf_env_points(im_1d,bg_1d,ctf) :
#	"""This will return a list of x,y points corresponding to the maxima of the ctf in the background
#	subtracted power spectrum"""
#	ys=len(bg_1d)
#	ds=ctf.dsbg
#	cc=ctf.compute_1d(ys,ds,Ctf.CtfType.CTF_AMP)
#	ret=[]
#
#	for i in range(1,len(cc)-1):
#		if cc[i-1]<cc[i] and cc[i]>cc[i+1] and im_1d[i]-bg_1d[i]>0 :
#			ret.append((i*ds,(im_1d[i]-bg_1d[i])))
##			ret.append((i*ds,(im_1d[i]-bg_1d[i])/sfact(i*ds)))		# this version removes the structure factor (in theory)
#
#	return ret

try:
	from PyQt4 import QtCore, QtGui, QtOpenGL
	from PyQt4.QtCore import Qt
	from OpenGL import GL,GLUT
	from emshape import *
	from valslider import ValSlider,CheckBox
except:
	print "Warning: PyQt4 must be installed to use the --gui option"
	class dummy:
		pass
	class QWidget:
		"A dummy class for use when Qt not installed"
		def __init__(self,parent):
			print "Qt4 has not been loaded"
	class QListWidget:
		"A dummy class"
		def __init__(self,parent):
			print "Qt4 has not been loaded"
	QtGui=dummy()
	QtGui.QWidget=QWidget
	QtGui.QListWidget=QListWidget

def notzero(x):
	if x==0 : return 1.0
	return x

class MyListWidget(QtGui.QListWidget):
	"""Exactly like a normal list widget but intercepts a few keyboard events"""

	def keyPressEvent(self,event):

		if event.key() in (Qt.Key_Up,Qt.Key_Down) :
			QtGui.QListWidget.keyPressEvent(self,event)
			return

		self.emit(QtCore.SIGNAL("keypress"),event)
#		event.key()==Qt.Key_I


class GUIctf(QtGui.QWidget):
	def __init__(self,application,data,autohp=True,nosmooth=False,highdensity=False):
		"""Implements the CTF fitting dialog using various EMImage and EMPlot2D widgets
		'data' is a list of (filename,EMAN2CTF,im_1d,bg_1d,im_2d,bg_2d,qual,bg_1d_low)
		"""
		try:
			from emimage2d import EMImage2DWidget
		except:
			print "Cannot import EMAN image GUI objects (EMImage2DWidget)"
			sys.exit(1)
		try:
			from emplot2d import EMPlot2DWidget
		except:
			print "Cannot import EMAN plot GUI objects (is matplotlib installed?)"
			sys.exit(1)

		self.app = weakref.ref(application)
		self.autohp=autohp
		self.nosmooth=nosmooth
		self.highdensity=highdensity

		QtGui.QWidget.__init__(self,None)
		self.setWindowIcon(QtGui.QIcon(get_image_directory() + "ctf.png"))

		self.data=data
		self.curset=0
		self.plotmode=0

		self.guiim=EMImage2DWidget(application=self.app())
		self.guiiminit = True # a flag that's used to auto resize the first time the gui's set_data function is called
		self.guiplot=EMPlot2DWidget(application=self.app())
		self.guirealim=EMImage2DWidget(application=self.app())	# This will show the original particle images
		self.flipim=None

		self.guirealim.connect(self.guirealim,QtCore.SIGNAL("keypress"),self.realimgkey)
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mousedown"),self.imgmousedown)
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mousedrag"),self.imgmousedrag)
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mouseup")  ,self.imgmouseup)
		self.guiplot.connect(self.guiplot,QtCore.SIGNAL("mousedown"),self.plotmousedown)



		self.guiim.mmode="app"

		# This object is itself a widget we need to set up
		self.hbl = QtGui.QHBoxLayout(self)
		self.hbl.setMargin(0)
		self.hbl.setSpacing(6)
		self.hbl.setObjectName("hbl")

		# plot list and plot mode combobox
		self.vbl2 = QtGui.QVBoxLayout()
		self.setlist=MyListWidget(self)
		self.setlist.setSizePolicy(QtGui.QSizePolicy.Preferred,QtGui.QSizePolicy.Expanding)
		self.vbl2.addWidget(self.setlist)

		self.splotmode=QtGui.QComboBox(self)
		self.splotmode.addItem("Bgsub & fit")
		self.splotmode.addItem("Ptcl & BG power")
		self.splotmode.addItem("SNR")
		self.splotmode.addItem("Smoothed SNR")
		self.splotmode.addItem("Integrated SNR")
		self.splotmode.addItem("Total CTF")
		self.splotmode.addItem("Total SNR")
		self.splotmode.addItem("All SNR (peak)")
		self.splotmode.addItem("SNR Scaling")
		self.splotmode.addItem("SNR vs Defocus")
		self.splotmode.addItem("<debug>")
		self.vbl2.addWidget(self.splotmode)
		self.hbl.addLayout(self.vbl2)

		# ValSliders for CTF parameters
		self.vbl = QtGui.QVBoxLayout()
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		self.hbl.addLayout(self.vbl)

		#self.samp = ValSlider(self,(0,5.0),"Amp:",0)
		#self.vbl.addWidget(self.samp)

		self.imginfo=QtGui.QLabel("Info",self)
		self.vbl.addWidget(self.imginfo)

		self.sdefocus=ValSlider(self,(0,5),"Defocus:",0,90)
		self.vbl.addWidget(self.sdefocus)

		self.sbfactor=ValSlider(self,(0,1600),"B factor:",0,90)
		self.vbl.addWidget(self.sbfactor)

		self.sdfdiff=ValSlider(self,(0,1),"DF Diff:",0,90)
		self.vbl.addWidget(self.sdfdiff)

		self.sdfang=ValSlider(self,(0,180),"Df Angle:",0,90)
		self.vbl.addWidget(self.sdfang)

		self.sampcont=ValSlider(self,(0,100),"% AC",0,90)
		self.vbl.addWidget(self.sampcont)

#		self.sapix=ValSlider(self,(.2,10),"A/Pix:",2,90)
#		self.vbl.addWidget(self.sapix)

		self.svoltage=ValSlider(self,(0,500),"Voltage (kV):",0,90)
		self.vbl.addWidget(self.svoltage)

		self.scs=ValSlider(self,(0,5),"Cs (mm):",0,90)
		self.vbl.addWidget(self.scs)

		self.squality=ValSlider(self,(0,9),"Quality (0-9):",0,90)
		self.squality.setIntonly(True)
		self.vbl.addWidget(self.squality)


		self.hbl_buttons = QtGui.QHBoxLayout()
		self.saveparms = QtGui.QPushButton("Save parms")
		self.recallparms = QtGui.QPushButton("Recall")
		self.refit = QtGui.QPushButton("Refit")
		self.show2dfit = CheckBox(label="Show 2D Sim:",value=False)
		self.showzerorings = CheckBox(label="Show Zeroes:",value=False)
		self.output = QtGui.QPushButton("Output")
		self.hbl_buttons.addWidget(self.refit)
		self.hbl_buttons.addWidget(self.saveparms)
		self.hbl_buttons.addWidget(self.recallparms)
		self.hbl_buttons2 = QtGui.QHBoxLayout()
		self.hbl_buttons2.addWidget(self.show2dfit)
		self.hbl_buttons2.addWidget(self.showzerorings)
		self.hbl_buttons2.addWidget(self.output)
		self.vbl.addLayout(self.hbl_buttons)
		self.vbl.addLayout(self.hbl_buttons2)

		QtCore.QObject.connect(self.sdefocus, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.sbfactor, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.sdfdiff, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.sdfang, QtCore.SIGNAL("valueChanged"), self.newCTF)
#		QtCore.QObject.connect(self.sapix, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.sampcont, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.svoltage, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.scs, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.squality, QtCore.SIGNAL("valueChanged"), self.newQual)
		QtCore.QObject.connect(self.showzerorings, QtCore.SIGNAL("valueChanged"), self.update_plot)
		QtCore.QObject.connect(self.setlist,QtCore.SIGNAL("currentRowChanged(int)"),self.newSet)
		QtCore.QObject.connect(self.setlist,QtCore.SIGNAL("keypress"),self.listkey)
		QtCore.QObject.connect(self.splotmode,QtCore.SIGNAL("currentIndexChanged(int)"),self.newPlotMode)

	   	QtCore.QObject.connect(self.saveparms,QtCore.SIGNAL("clicked(bool)"),self.on_save_params)
		QtCore.QObject.connect(self.recallparms,QtCore.SIGNAL("clicked(bool)"),self.on_recall_params)
		QtCore.QObject.connect(self.refit,QtCore.SIGNAL("clicked(bool)"),self.on_refit)
		QtCore.QObject.connect(self.output,QtCore.SIGNAL("clicked(bool)"),self.on_output)

		self.neednewps=False
		self.update_data()

		self.resize(720,380) # figured these values out by printing the width and height in resize event


		E2loadappwin("e2ctf","main",self)
		E2loadappwin("e2ctf","image",self.guiim.qt_parent)
		E2loadappwin("e2ctf","realimage",self.guirealim.qt_parent)
		E2loadappwin("e2ctf","plot",self.guiplot.qt_parent)

		self.setWindowTitle("CTF")

	def listkey(self,event):

		if event.key()>=Qt.Key_0 and event.key()<=Qt.Key_9 :
			q=int(event.key())-Qt.Key_0
			self.squality.setValue(q)
		elif event.key() == Qt.Key_Left:
			self.sdefocus.setValue(self.sdefocus.getValue()-0.01)
		elif event.key() == Qt.Key_Right:
			self.sdefocus.setValue(self.sdefocus.getValue()+0.01)
		elif event.key() == Qt.Key_B:
			self.sbfactor.setValue(100.0)
		elif event.key()==Qt.Key_S :
			print "Save Parms ",str(self.setlist.item(self.curset).text())
			self.on_save_params()
		elif event.key()==Qt.Key_R :
			self.on_recall_params()


	def on_save_params(self):

		if len(self.setlist.selectedItems()) == 0: return

		val = self.curset
		name = base_name(str(self.setlist.item(val).text()))

		try: js_parms = js_open_dict(info_name(str(self.setlist.item(val).text())))
		except:
			print "Error writing CTF parameters for {}".format(name)
			return

#		if not db_check_dict(name):
#			print "error, the db doesn't exist:",name
#

		tmp=js_parms["ctf"]
		tmp[0]=self.data[val][1]	# EMAN2CTF object
		js_parms["ctf"]=tmp

		js_parms["quality"]=self.data[val][6]

	def on_recall_params(self):
		if len(self.setlist.selectedItems()) == 0: return

		val = self.curset
		name = base_name(str(self.setlist.item(val).text()))

		try: js_parms = js_open_dict(info_name(str(self.setlist.item(val).text())))
		except:
			print "Error reading CTF parameters for {}".format(name)
			return

		self.data[val][1]=js_parms["ctf"][0]
		self.data[val][6]=js_parms["quality"]
		self.newSet(self.curset)

#	def get_output_params(self):

	def on_refit(self):
		# self.data[n] contains filename,EMAN2CTF,im_1d,bg_1d,im_2d,bg_2d,qual
		tmp=list(self.data[self.curset])

		dfdiff=self.sdfdiff.value
		dfang=self.sdfang.value
		ctf=ctf_fit(tmp[2],tmp[3],tmp[7],tmp[4],tmp[5],tmp[1].voltage,tmp[1].cs,tmp[1].ampcont,tmp[1].apix,bgadj=not self.nosmooth,autohp=self.autohp,dfhint=self.sdefocus.value,highdensity=self.highdensity)
		ctf.dfdiff=dfdiff
		ctf.dfang=dfang
		if ctf.dfdiff!=0 :
			ctf_fit_stig(tmp[4],tmp[5],ctf,True)

		tmp[2],tmp[3]=calc_1dfrom2d(ctf,tmp[4],tmp[5])			# update 1-D curves (impacted if astigmatism changed)

		ctf.bfactor=ctf_fit_bfactor(list(array(tmp[2])-array(tmp[3])),ctf.dsbg,ctf)
		ctf.snr=[snr_safe(tmp[2][i],tmp[3][i]) for i in range(len(tmp[2]))]

		val = self.curset
		name = base_name(str(self.setlist.item(val).text()))

		try: js_parms = js_open_dict(info_name(str(self.setlist.item(val).text())))
		except:
			print "Error writing CTF parameters for {}".format(name)
			return

#		if not db_check_dict(name):
#			print "error, the db doesn't exist:",name
#

		tmp=js_parms["ctf"]
		tmp[0]=ctf	# EMAN2CTF object
		js_parms["ctf"]=tmp

		self.data[self.curset][1]=ctf
		self.newSet(self.curset)

	def on_output(self):
		from emsprworkflow import E2CTFOutputTaskGeneral

		n = self.setlist.count()
		names = [str(self.setlist.item(i).text()) for i in xrange(0,n)]

		self.form = E2CTFOutputTaskGeneral()
		self.form.set_names(names)
		self.form.run_form()

	def show_guis(self):
		if self.guiim != None:
			self.app().show_specific(self.guiim)
		if self.guiplot != None:
			self.app().show_specific(self.guiplot)
		if self.guirealim != None:
			self.app().show_specific(self.guirealim)

		self.show()

	def closeEvent(self,event):
#		QtGui.QWidget.closeEvent(self,event)
#		self.app.app.closeAllWindows()
		E2saveappwin("e2ctf","main",self)

		if self.guiim != None:
			E2saveappwin("e2ctf","image",self.guiim.qt_parent)
			self.app().close_specific(self.guiim)
			self.guiim = None
		if self.guiplot != None:
			E2saveappwin("e2ctf","plot",self.guiplot.qt_parent)
			self.app().close_specific(self.guiplot)
		if self.guirealim != None:
			E2saveappwin("e2ctf","realimage",self.guirealim.qt_parent)
			self.app().close_specific(self.guirealim)

		event.accept()
		self.app().close_specific(self)
		self.emit(QtCore.SIGNAL("module_closed")) # this signal is important when e2ctf is being used by a program running its own event loop

	def newData(self,data):
		self.data=data
		self.update_data()

	def update_data(self):
		"""This will make sure the various widgets properly show the current data sets"""
		self.setlist.clear()
		for i,j in enumerate(self.data):
			self.setlist.addItem(j[0])
		self.setlist.setCurrentRow(self.curset)

	def update_plot(self):
		if self.guiplot == None: return # it's closed/not visible
		val=self.curset
		ctf=self.data[val][1]
		ds=ctf.dsbg
		r=len(ctf.background)
		s=[ds*i for i in range(r)]

		if r==0 or len(s)<2 :
			print "Trying to plot bad data (set %s): %s"%(str(val),str(ctf))
			return

		# With astigmatism support, the curves may need an update
		# recompute self.data[val][2] and [3] (FG and BG 1-D curves) if necessary
		if self.neednewps :
			self.data[val][2],self.data[val][3]=calc_1dfrom2d(ctf,self.data[val][4],self.data[val][5])
			self.data[val][7]=low_bg_curve(self.data[val][3],ds)
			ctf.snr=[snr_safe(self.data[val][2][i],self.data[val][3][i]) for i in range(r)]

		# This updates the image circles
		fit=ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP)
		shp={}

		# Old way of finding zeroes from the 1-d curve dynamically
		#nz=0
		#for i in range(1,len(fit)):
			#if fit[i-1]*fit[i]<=0.0:
				#nz+=1
				#shp["z%d"%i]=EMShape(("circle",0.0,0.0,1.0/nz,r,r,i,1.0))
##				if nz==1: print ("circle",0.0,0.0,1.0/nz,r,r,i,1.0)

		if self.show2dfit.getValue() and self.guiim != None and self.flipim != None:
			ctf.compute_2d_complex(self.flipim,Ctf.CtfType.CTF_FITREF)
			self.flipim.mult(self.data[val][4]["sigma"]/self.flipim["sigma"])
			self.flipim.update()
			self.guiim.set_data([self.data[val][4],self.flipim])
#			self.guiim.set_data(self.data[val][4])
		else:
			self.guiim.set_data(self.data[val][4])

		# We draw the first 5 zeroes with computed zero locations
		if self.showzerorings.getValue():
			for i in range(1,10):
				if ctf.dfdiff>0 :
					#z1=zero(i,ctf.voltage,ctf.cs,ctf.defocus-ctf.dfdiff/2,ctf.ampcont)/ctf.dsbg
					#z2=zero(i,ctf.voltage,ctf.cs,ctf.defocus+ctf.dfdiff/2,ctf.ampcont)/ctf.dsbg
					d=ctf.defocus
					ctf.defocus=d-ctf.dfdiff/2
					z1=ctf.zero(i-1)/ctf.dsbg
					ctf.defocus=d+ctf.dfdiff/2
					z2=ctf.zero(i-1)/ctf.dsbg
					ctf.defocus=d
					if z2>len(s) : break
					shp["z%d"%i]=EMShape(("ellipse",0,0,.75,r,r,z2,z1,ctf.dfang,1.0))
				else:
	#				z=zero(i,ctf.voltage,ctf.cs,ctf.defocus,ctf.ampcont)/ctf.dsbg
					z=ctf.zero(i-1)/ctf.dsbg
					if z>len(s) : break
					shp["z%d"%i]=EMShape(("circle",0,0,.75,r,r,z,1.0))

		self.guiim.del_shapes()
		if len(shp)>0 : self.guiim.add_shapes(shp)
		self.guiim.updateGL()

		if self.plotmode==1:
			self.guiplot.set_data((s,self.data[val][2]),"fg",True,True,color=1)
			self.guiplot.set_data((s,self.data[val][7]),"bg(concave)",quiet=True,color=0,linetype=2)
			if len(self.data[val])>8:
				self.guiplot.set_data((s,self.data[val][8]),"Micrograph",quiet=True,color=2,linetype=2)
			self.guiplot.set_data((s,self.data[val][3]),"bg",color=0)
			self.guiplot.setAxisParms("s (1/"+ "$\AA$" +")","Intensity (a.u)")
		elif self.plotmode==0:
			bgsub=[self.data[val][2][i]-self.data[val][3][i] for i in range(len(self.data[val][2]))]
			self.guiplot.set_data((s,bgsub),"fg-bg",True,True,color=0)

			lowbg=self.data[val][7]
			lowbgsub=[self.data[val][2][i]-lowbg[i] for i in range(len(self.data[val][2]))]
			self.guiplot.set_data((s,lowbgsub),"fg-lowbg",False,True,color=0,linetype=2)

			# This means we have a whole micrograph curve
			if len(self.data[val])>8:
				bgsub2=[self.data[val][8][i]-self.data[val][3][i] for i in range(len(self.data[val][2]))]
				self.guiplot.set_data((s,bgsub2),"micro-bg",False,True,color=2,linetype=2)

			fit=ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP)		# The fit curve
			fit=[sfact2(s[i])*fit[i]**2 for i in range(len(s))]		# squared * structure factor

			# auto-amplitude for b-factor adjustment
			rto,nrto=0,0
			for i in range(int(.04/ds)+1,min(int(0.15/ds),len(s)-1)):
				if bgsub[i]>0 :
					#rto+=fit[i]**2/fabs(bgsub[i])
					#nrto+=fit[i]
					#rto+=fit[i]**2
					#nrto+=bgsub[i]**2
					rto+=fit[i]
					nrto+=fabs(bgsub[i])
			if nrto==0 : rto=1.0
			else : rto/=nrto
			fit=[fit[i]/rto for i in range(len(s))]

#			print ctf_cmp((self.sdefocus.value,self.sbfactor.value,rto),(ctf,bgsub,int(.04/ds)+1,min(int(0.15/ds),len(s)-1),ds,self.sdefocus.value))

			self.guiplot.set_data((s,fit),"fit",color=1)
			self.guiplot.setAxisParms("s (1/"+ "$\AA$" + ")","Intensity (a.u)")
		elif self.plotmode==2:
			snr=ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_SNR)		# The snr curve
			self.guiplot.set_data((s,snr[:len(s)]),"snr",True)
			self.guiplot.setAxisParms("s (1/"+ "$\AA$" +")","SNR (intensity ratio)")
		elif self.plotmode==3:
			global sfcurve2
			snr=ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_SNR)		# The snr curve
			self.guiplot.set_data((s,snr[:len(s)]),"SNR",True,True,color=0)
			ssnr=ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_SNR_SMOOTH,sfcurve2)		# The smoothed curve
			self.guiplot.set_data((s,ssnr[:len(s)]),"Smoothed SNR",color=1)

			# Also display the original unfiltered SNR, which may differ from the CTF stored SNR
			#osnr=[snr_safe(self.data[val][2][i],self.data[val][3][i]) for i in range(len(s))]
			#self.guiplot.set_data((s,osnr),"Original SNR",color=2)

			# we can also optionally display the computed SNR used in the smoothing process
			#csnr=ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP)		# The fit curve
			#csnr=[sfact2(s[i])*csnr[i]**2/notzero(self.data[val][3][i]) for i in range(len(s))]		# squared * structure factor/background
			#self.guiplot.set_data((s,csnr[:len(s)]),"Computed SNR",color=2)

			self.guiplot.setAxisParms("s (1/"+ "$\AA$" +")","SNR (intensity ratio)")
		elif self.plotmode==4:
			snr=ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_SNR)		# The snr curve
			for i in range(1,len(snr)): snr[i]=snr[i]*i+snr[i-1]			# integrate SNR*s
#			for i in range(1,len(snr)): snr[i]/=snr[-1]				# normalize
			for i in range(1,len(snr)): snr[i]/=len(snr)			# this way the relative quality of images can be compared
			self.guiplot.set_data((s,snr[:len(s)]),"snr",True)
			self.guiplot.setAxisParms("s (1/"+ "$\AA$" +")","Integrated SNR")
		# Total CTF
		elif self.plotmode==5:
			inten=[fabs(i) for i in ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP)]		# The snr curve
			self.guiplot.set_data((s,inten[:len(s)]),"single",True,quiet=True)
			all=[0 for i in inten]
			for st in self.data:
#				print st
				inten=[fabs(i) for i in st[1].compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP)]
				for i in range(len(all)): all[i]+=inten[i]
			self.guiplot.set_data((s,all[:len(s)]),"total",False,True)
			self.guiplot.setAxisParms("s (1/"+ "$\AA$" +")","CTF Sum")
		# Total SNR
		elif self.plotmode==6:
			inten=[fabs(i) for i in ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_SNR)]		# The snr curve
#			self.guiplot.set_data((s,inten[:len(s)]),"single",True)
			allt=[0 for i in inten]
			allw=[0 for i in inten]
			for st in self.data:
#				print st
				inten=st[1].compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_SNR)
				for i in range(len(allt)):
					allt[i]+=inten[i]
					allw[i]+=inten[i]*st[4]["ptcl_repr"]/((i+1.0)*pi)  # *number of particles / number of pixels
			self.guiplot.set_data((s,allt[:len(s)]),"Total",True,True)

			# This is a sort of estimated 3-D SNR assuming each pixel contributes to 2 Fourier voxels
			# and not taking symmetry or oversampling into account. Useful for VERY rough estimates only
			self.guiplot.set_data((s,allw[:len(s)]),"Total 3D /(sym*oversamp)",False,True)

			self.guiplot.setAxisParms("s (1/"+ "$\AA$" +")","SNR Sum")
		# All SNR
		elif self.plotmode==7:
			self.guiplot.set_data(None,None,True,True)		# erase existing data quietly
			for st in self.data:
				inten=st[1].compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_SNR)
				inten=[(s[i],inten[i]) for i in xrange(1,len(inten)-1) if inten[i]>inten[i+1] and inten[i]>inten[i-1]]
				ls,li=zip(*inten)	# this confusing idiom is unzipping the list
#				print len(ls),len(li)
				self.guiplot.set_data((ls,li),st[0],quiet=True,linetype=-1,symtype=0,symsize=2)

			self.guiplot.setAxisParms("s (1/"+ "$\AA$" +")","SNR")
#			self.guiplot.updateGL()
		# All SNR vs defocus
		elif self.plotmode==9:
			dflist=[d[1].defocus for d in self.data]
			snrlist=[sum(d[1].snr)/len(d[1].snr) for d in self.data]
			self.guiplot.set_data((dflist,snrlist),"df vs snr",replace=True,linetype=-1,symtype=0,symsize=2)		# erase existing data quietly

			self.guiplot.setAxisParms("Defocus (um)","Mean SNR")
#			self.guiplot.updateGL()

		# Debug
		elif self.plotmode==10:
#			bgsub=[self.data[val][2][i]-self.data[val][3][i] for i in range(len(self.data[val][2]))]
#			self.guiplot.set_data("fg-bg",(s,bgsub),True,True)

			#fit=ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP)		# The fit curve
			#fit=[sfact2(s[i])*fit[i]**2 for i in range(len(s))]		# squared * a generic structure factor

			#fit=[bgsub[i]/fit[i] for i in range(len(s))]
			#for i in range(len(fit)) :
				#if fabs(fit[i])>10000 : fit[i]=0

			#a=[i**2 for i in bgsub[33:]]
			#b=[i**2 for i in fit[33:]]
			#r=sqrt(sum(a)*sum(b))
			#fit=[bgsub[i]*fit[i]/r for i in range(len(s))]

			#self.guiplot.set_data((s,fit),"fit",True)
			#self.guiplot.setAxisParms("s (1/A)","Intensity (a.u)")

			#bgsub=[self.data[val][2][i]-self.data[val][3][i] for i in range(len(self.data[val][2]))]
			#self.guiplot.set_data("fg-bg",(s,bgsub),True,True)

			#fit=[bgsub[i]/sfact(s[i]) for i in range(len(s))]		# squared * a generic structure factor

			#self.guiplot.set_data("fit",(s,fit))

			bgsub=[self.data[val][2][i]-self.data[val][3][i] for i in range(len(self.data[val][2]))]
			self.guiplot.set_data((s,bgsub),"fg-bg",True,True,color=0)

			fit=ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP)		# The fit curve
			fit=[sfact2(s[i])*fit[i]**2 for i in range(len(s))]		# squared * structure factor

			# auto-amplitude for b-factor adjustment
			rto,nrto=0,0
			for i in range(int(.04/ds)+1,min(int(0.15/ds),len(s)-1)):
				if bgsub[i]>0 :
					rto+=fit[i]
					nrto+=fabs(bgsub[i])
			if nrto==0 : rto=1.0
			else : rto/=nrto
			fit=[fit[i]/rto for i in range(len(s))]
			self.guiplot.set_data((s,fit),"fit",color=1)

			if len(bgsub)<64 : wdw=4
			elif len(bgsub)<128 : wdw=6
			elif len(bgsub)<256 : wdw=8
			else : wdw=10
			sim=Util.windowdot(bgsub,fit,wdw,1)
			self.guiplot.set_data((s,sim),"Local Sim",color=2)

			print sum(sim),sum(sim[int(.04/ds):int(.12/ds)])

#			print ctf_cmp((self.sdefocus.value,self.sbfactor.value,rto),(ctf,bgsub,int(.04/ds)+1,min(int(0.15/ds),len(s)-1),ds,self.sdefocus.value))

			self.guiplot.setAxisParms("s (1/"+ "$\AA$" + ")","Intensity (a.u)")


		# SNR Scaling
		elif self.plotmode==8:
			# SNR computed from fit and background
			fit=ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP)		# The fit curve
			fit=[sfact2(s[i])*fit[i]**2/max(.001,self.data[val][3][i]) for i in range(len(s))]		# squared * structure factor/background

			# SNR from data
			snr=ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_SNR)

			self.guiplot.set_data((fit,snr),"SNR plot",True,True)
			self.guiplot.setAxisParms("SNR (fit)","SNR (meas)")


	def newSet(self,val):
		"called when a new data set is selected from the list"
		self.curset=val

		self.sdefocus.setValue(self.data[val][1].defocus,True)
		self.sbfactor.setValue(self.data[val][1].bfactor,True)
#		self.sapix.setValue(self.data[val][1].apix)
		self.sampcont.setValue(self.data[val][1].ampcont,True)
		self.svoltage.setValue(self.data[val][1].voltage,True)
		self.scs.setValue(self.data[val][1].cs,True)
		self.sdfdiff.setValue(self.data[val][1].dfdiff,True)
		self.sdfang.setValue(self.data[val][1].dfang,True)
		self.squality.setValue(self.data[val][6],True)
		if self.sdfdiff.value!=0 : self.neednewps=True

		try : ptcl=str(self.data[val][4]["ptcl_repr"])
		except: ptcl="?"
		try:
			ctf=self.data[val][1]
			ssnr="%1.4f"%(sum(ctf.snr)/float(len(ctf.snr)))
		except:
			ssnr="?"
		self.imginfo.setText("%s particles     SNR = %s"%(ptcl,ssnr))

		if self.guiim != None:
#			print self.data
			self.flipim=self.data[val][4].copy()
#			self.guiim.set_data([self.data[val][4],self.flipim])
			self.guiim.set_data(self.data[val][4])
			if self.guiiminit:
				self.guiim.optimally_resize()
				self.guiiminit = False
			self.guiim.updateGL()
		self.update_plot()

#		print "self.data[val]=",self.data[val][0].split('#')[-1]


		self.guiim.qt_parent.setWindowTitle("e2ctf - 2D FFT - "+self.data[val][0].split('#')[-1])
		self.guirealim.qt_parent.setWindowTitle("e2ctf - "+self.data[val][0].split('#')[-1])
		self.guiplot.qt_parent.setWindowTitle("e2ctf - Plot - "+self.data[val][0].split('#')[-1])

		n=EMUtil.get_image_count(self.data[val][0])
		if n>1:
			self.ptcldata=EMData.read_images(self.data[val][0],range(0,min(20,n)))
			im=sum(self.ptcldata)
			im.mult(4.0/len(self.ptcldata))	# 4 compensatess for noise averaging
			self.ptcldata.insert(0,im)
			self.guirealim.set_data(self.ptcldata)
		else : self.guirealim.set_data([EMData()])

	def newPlotMode(self,mode):
		self.plotmode=mode
		self.update_plot()

	def newCTF(self) :
		if self.data[self.curset][1].dfdiff!=0 or self.data[self.curset][1].dfdiff!=self.sdfdiff.value:
			if self.data[self.curset][1].dfdiff!=self.sdfdiff.value or self.data[self.curset][1].dfang!=self.sdfang.value or self.data[self.curset][1].defocus!=self.sdefocus.value :
				self.neednewps=True
		self.data[self.curset][1].defocus=self.sdefocus.value
		self.data[self.curset][1].bfactor=self.sbfactor.value
		self.data[self.curset][1].dfdiff=self.sdfdiff.value
		self.data[self.curset][1].dfang=self.sdfang.value
#		self.data[self.curset][1].apix=self.sapix.value
		self.data[self.curset][1].ampcont=self.sampcont.value
		self.data[self.curset][1].voltage=self.svoltage.value
		self.data[self.curset][1].cs=self.scs.value
		self.update_plot()

	def newQual(self):
		self.data[self.curset][6]=int(self.squality.value)

		val = self.curset
		name = base_name(str(self.setlist.item(val).text()))

		try: js_parms = js_open_dict(info_name(str(self.setlist.item(val).text())))
		except:
			print "Error writing CTF parameters for {}".format(name)
			return

		js_parms["quality"]=self.data[val][6]

	def realimgkey(self,event):
		"""Keypress in the image display window"""

		if event.key()==Qt.Key_I:			# if user presses I in this window we invert the stack on disk
			fsp=self.data[self.curset][0]
			n=EMUtil.get_image_count(fsp)
			print "Inverting images in %s"%fsp
			for i in xrange(n):
				img=EMData(fsp,i)
				img.mult(-1.0)
				img.write_image(fsp,i)

			self.ptcldata=EMData.read_images(fsp,range(0,20))
			self.guirealim.set_data(self.ptcldata)


	def imgmousedown(self,event) :
		m=self.guiim.scr_to_img((event.x(),event.y()))
		#self.guiim.add_shape("cen",["rect",.9,.9,.4,x0,y0,x0+2,y0+2,1.0])

	def imgmousedrag(self,event) :
		m=self.guiim.scr_to_img((event.x(),event.y()))

		# box deletion when shift held down
		#if event.modifiers()&Qt.ShiftModifier:
			#for i,j in enumerate(self.boxes):

	def imgmouseup(self,event) :
		m=self.guiim.scr_to_img((event.x(),event.y()))

	def plotmousedown(self,event) :
		m=self.guiim.scr_to_img((event.x(),event.y()))

	def run(self):
		"""If you make your own application outside of this object, you are free to use
		your own local app.exec_(). This is a convenience for ctf-only programs."""
		self.app.exec_()

#		E2saveappwin("boxer","imagegeom",self.guiim)
#		try:
#			E2setappval("boxer","imcontrol",self.guiim.inspector.isVisible())
#			if self.guiim.inspector.isVisible() : E2saveappwin("boxer","imcontrolgeom",self.guiim.inspector)
#		except : E2setappval("boxer","imcontrol",False)

		return


if __name__ == "__main__":
	main()
