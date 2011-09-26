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
from OpenGL import GL,GLUT
from math import *
import os
import sys
import weakref

from Simplex import Simplex

debug=False
logid=None

sfcurve=None		# This will store a global structure factor curve if specified
sfcurve2=None
envelopes=[]		# simplex minimizer needs to use a global at the moment

def main():
	global debug,logid
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <input stack/image> ...
	
Various CTF-related operations on images, including automatic fitting. Note that automatic fitting is limited to defocuses
less than ~5 microns. Input particles should be unmasked and unfiltered. A minimum of ~20% padding around the
particles is required for background extraction, even if this brings the edge of another particle into the box in some cases.
Particles should be reasonably well centered. Can also optionally phase flip and Wiener filter particles. Wiener filtration comes
after phase-flipping, so if phase flipping is performed Wiener filtered particles will also be phase-flipped. Note that both
operations are performed on oversampled images if specified, but this will produce phase-flipped images which are irreversable,
so, while oversampling can be useful for fitting, it is not recommended for phase-flipping.

Increasing padding during the particle picking process will improve the accuracy of phase-flipping, particularly for
images far from focus."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_pos_argument(name="particles",help="List the file to process with e2ctf here.", default="", guitype='filebox', positional=True, row=0, col=0,rowspan=1, colspan=2, mode='autofit,tuning,genoutp,gensf')
	parser.add_header(name="ctfheader", help='Options below this label are specific to e2ctflassaverage3d', title="### e2ctf options ###", default=None, guitype='filebox', row=1, col=0, rowspan=1, colspan=2, mode="autofit,tuning,genoutp,gensf")
	parser.add_argument("--allparticles",action="store_true",help="Will process all particle sets stored in BDB in the particles subdirectory witho",default=False)
	parser.add_argument("--minptcl",type=int,help="Files with fewer than the specified number of particles will be skipped",default=0)
	parser.add_argument("--gui",action="store_true",help="Start the GUI for interactive fitting",default=False, guidefault=True, guitype='boolbox', row=2, col=0, rowspan=1, colspan=1, mode="tuning")
	parser.add_argument("--autofit",action="store_true",help="Runs automated CTF fitting on the input images",default=False, guidefault=True, guitype='boolbox', row=6, col=0, rowspan=1, colspan=1, mode='autofit')
	parser.add_argument("--curdefocushint",action="store_true",help="Rather than doing the defocus from scratch, use existing values in the project as a starting point",default=False, guidefault=True, guitype='boolbox', row=6, col=1, rowspan=1, colspan=1, mode='autofit')
	parser.add_argument("--refinebysnr",action="store_true",help="Refines the defocus value by looking at the high resolution smoothed SNR. Important: also replaces the SNR with a smoothed version.",default=False)
	parser.add_argument("--bgmask",type=int,help="Background is computed using a soft mask of the center/edge of each particle with the specified radius. Default radius is boxsize/2.6.",default=0)
	parser.add_argument("--fixnegbg",action="store_true",help="Will perform a final background correction to avoid slight negative values near zeroes")
	parser.add_argument("--computesf",action="store_true",help="Will determine the structure factor*envelope for the aggregate set of images", guidefault=True, guitype='boolbox', row=3, col=0, rowspan=1, colspan=1, mode="gensf")
	parser.add_argument("--apix",type=float,help="Angstroms per pixel for all images",default=0, guitype='intbox', row=3, col=0, rowspan=1, colspan=1, mode='autofit')
	parser.add_argument("--voltage",type=float,help="Microscope voltage in KV",default=0, guitype='intbox', row=3, col=1, rowspan=1, colspan=1, mode='autofit')
	parser.add_argument("--cs",type=float,help="Microscope Cs (spherical aberation)",default=0, guitype='intbox', row=4, col=0, rowspan=1, colspan=1, mode='autofit')
	parser.add_argument("--ac",type=float,help="Amplitude contrast (percentage, default=10)",default=10, guitype='intbox', row=4, col=1, rowspan=1, colspan=1, mode='autofit')
	parser.add_argument("--dfmax",type=float,help="Maximum defocus for autofitting (default=6.4)",default=None)
	parser.add_argument("--autohp",action="store_true",help="Automatic high pass filter of the SNR only to remove initial sharp peak, phase-flipped data is not directly affected (default false)",default=False, guidefault=True, guitype='boolbox', row=5, col=0, rowspan=1, colspan=1, mode='autofit')
	parser.add_argument("--invert",action="store_true",help="Invert the contrast of the particles in output files (default false)",default=False)
	parser.add_argument("--nonorm",action="store_true",help="Suppress per image real-space normalization",default=False)
	parser.add_argument("--nosmooth",action="store_true",help="Disable smoothing of the background (running-average of the log with adjustment at the zeroes of the CTF)",default=False, guitype='boolbox', row=5, col=1, rowspan=1, colspan=1, mode='autofit')
	parser.add_argument("--phaseflip",action="store_true",help="Perform phase flipping after CTF determination and writes to specified file.",default=False, guidefault=True, guitype='boolbox', row=3, col=0, rowspan=1, colspan=1, mode='genoutp')
	parser.add_argument("--phasefliphp",action="store_true",help="Perform phase flipping with auto-high pass filter",default=False, guitype='boolbox', row=4, col=0, rowspan=1, colspan=1, mode='genoutp')
	parser.add_argument("--wiener",action="store_true",help="Wiener filter (optionally phaseflipped) particles.",default=False, guidefault=True, guitype='boolbox', row=3, col=1, rowspan=1, colspan=1, mode='genoutp')
	parser.add_argument("--virtualout",type=str,help="Make a virtual stack copy of the input images with CTF parameters stored in the header. BDB only.",default=None)
	parser.add_argument("--storeparm",action="store_true",help="Write the CTF parameters back to the header of the input images. BDB and HDF only.",default=False)
	parser.add_argument("--oversamp",type=int,help="Oversampling factor",default=1, guitype='intbox', row=2, col=0, rowspan=1, colspan=2, mode='autofit')
	parser.add_argument("--sf",type=str,help="The name of a file containing a structure factor curve. Specify 'none' to use the built in generic structure factor. Default=auto",default="auto")
	parser.add_argument("--debug",action="store_true",default=False)
	parser.add_argument("--dbds",type=str,default=None,help="Data base dictionary storage, used by the workflow for storing which files have been filtered. You can ignore this argument")
	parser.add_argument("--source_image",type=str,default=None,help="Filters particles only with matching ptcl_source_image parameters in the header")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()
	
	if options.allparticles:
		args=["bdb:particles#"+i for i in db_list_dicts("bdb:particles") if "_ctf_" not in i]
		args.sort()
		if options.verbose : print "%d particle files identified"%len(args)
		
	if len(args)<1 : parser.error("Input image required")
	if options.autofit:
		if options.voltage==0 : parser.error("Please specify voltage")
		if options.cs==0 : parser.error("Please specify Cs")
	if options.apix==0 : print "Using A/pix from header"
		
	debug=options.debug
	img_sets=None

	global sfcurve,sfcurve2
	init_sfcurve(options.sf)

	logid=E2init(sys.argv)

 #	if options.oversamp>1 : options.apix/=float(options.oversamp)


	db_project=db_open_dict("bdb:project")
	db_parms=db_open_dict("bdb:e2ctf.parms")
#	db_misc=db_open_dict("bdb:e2ctf.misc")

	# remove any files that don't have enough particles from the list
	if options.minptcl>0 :
		args=[i for i in args if EMUtil.get_image_count(i)>=options.minptcl]

	options.filenames = args
		
	### Power spectrum and CTF fitting
	if options.autofit:
		img_sets=pspec_and_ctf_fit(options,debug) # converted to a function so to work with the workflow
		

	### GUI - user can update CTF parameters interactively
	if options.gui :
		img_sets = get_gui_arg_img_sets(options.filenames)
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
	if options.phaseflip or options.wiener or options.phasefliphp or options.virtualout or options.storeparm: # only put this if statement here to make the program flow obvious
		write_e2ctf_output(options) # converted to a function so to work with the workflow

	if options.computesf :
		img_sets = get_gui_arg_img_sets(options.filenames)
		print "Recomputing structure factor"
		envelope=compute_envelope(img_sets)
		
		db_misc=db_open_dict("bdb:e2ctf.misc")
		db_misc["strucfac"]=envelope
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
			db_misc=db_open_dict("bdb:e2ctf.misc",True)
			m=db_misc["strucfac"]
			print "Using previously generated structure factor from bdb:e2ctf.misc"
			sfcurve2=XYData()		# this is really slow and stupid
			for i,j in enumerate(m):
				sfcurve2.set_x(i,j[0])
				sfcurve2.set_y(i,j[1])
			
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
	returns the img_sets list required to intialized the GUI correctly
	'''
	
	# This is a REAL mess. Extremely non-transparent !
	img_sets = []
	if db_check_dict("bdb:e2ctf.parms"):
		db_parms=db_open_dict("bdb:e2ctf.parms",ro=True)
		db_im2d=db_open_dict("bdb:e2ctf.im2d")
		db_bg2d=db_open_dict("bdb:e2ctf.bg2d")
	else: return img_sets
	for fsp in filenames:
		name = get_file_tag(fsp)
		if name[:3]=="00i" : continue
		if not db_parms.has_key(name):
			print "Warning, you must first run auto fit before running the gui - there are no parameters for",name
			continue
		img_set = db_parms[name]
		ctf=EMAN2Ctf()
		ctf.from_string(img_set[0]) # convert to ctf object seeing as it's a string
		img_set = list(img_set)
		img_set[0]=ctf
		img_set.insert(-1,db_im2d[name])
		img_set.insert(-1,db_bg2d[name])
		img_set.append(low_bg_curve(img_set[2],ctf.dsbg))
		img_sets.append([fsp]+img_set)
		
	return img_sets

def write_e2ctf_output(options):
	"write wiener filtered and/or phase flipped particle data to the local database"
	global logid

	if options.phaseflip or options.wiener or options.phasefliphp :
		db_parms=db_open_dict("bdb:e2ctf.parms")
		for i,filename in enumerate(options.filenames):
			name=get_file_tag(filename)
			if debug: print "Processing ",filename

			if options.phaseflip: phaseout="bdb:particles#"+name+"_ctf_flip"
			else: phaseout=None
		
			if options.phasefliphp: phasehpout="bdb:particles#"+name+"_ctf_flip_hp"
			else: phasehpout=None
		
			if options.wiener:
				if options.autohp: wienerout="bdb:particles#"+name+"_ctf_wiener_hp"
				else: wienerout="bdb:particles#"+name+"_ctf_wiener"
			else : wienerout=None

			#if options.phaseflip: phaseout=name+"_ctf_flip.hed"
			#else: phaseout=None
		
			#if options.wiener: wienerout=name+"_ctf_wiener.hed"
			#else : wienerout=None
			
			if phaseout : print "Phase image out: ",phaseout,"\t",
			if phasehpout : print "Phase-hp image out: ",phasehpout,"\t",
			if wienerout : print "Wiener image out: ",wienerout,
			print ""
			ctf=EMAN2Ctf()
			ctf.from_string(db_parms[name][0])
			process_stack(filename,phaseout,phasehpout,wienerout,not options.nonorm,options.oversamp,ctf,invert=options.invert,virtualout=options.virtualout,storeparm=options.storeparm,source_image=options.source_image)
			if options.dbds != None:
				pdb = db_open_dict("bdb:project")
				dbds = pdb.get(options.dbds,dfl={})
				data_entry = {}
				file_key = get_file_tag(filename)
				if dbds.has_key(file_key):
					data_entry = dbds[file_key]
				
				if phaseout:
					data_entry["Phase flipped"] = phaseout
				if phasehpout:
					data_entry["Phase flipped-hp"] = phasehpout
				if wienerout:
					data_entry["Wiener filtered"] = wienerout
					
				dbds[file_key] = data_entry
				pdb[options.dbds] = dbds
					
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
			scales=simp.minimize(maxiters=5000)[0]
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
		for i,j in enumerate(envelope):
			if j[0]>=smax :break
		
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
	filename,EMAN2CTF,im_1d,bg_1d,im_2d,bg_2d,qual"""
	global logid
	img_sets=[]
	db_parms=db_open_dict("bdb:e2ctf.parms")
	db_im2d=db_open_dict("bdb:e2ctf.im2d")
	db_bg2d=db_open_dict("bdb:e2ctf.bg2d")

	for i,filename in enumerate(options.filenames):
		name=get_file_tag(filename)

		# compute the power spectra
		if debug : print "Processing ",filename
		apix=options.apix
		if apix<=0 : apix=EMData(filename,0,1)["apix_x"] 
		im_1d,bg_1d,im_2d,bg_2d,bg_1d_low=powspec_with_bg(filename,options.source_image,radius=options.bgmask,edgenorm=not options.nonorm,oversamp=options.oversamp,apix=apix)
		ds=1.0/(apix*im_2d.get_ysize())
		if not options.nosmooth : bg_1d=smooth_bg(bg_1d,ds)
		if options.fixnegbg : 
			bg_1d=fixnegbg(bg_1d,im_1d,ds)		# This insures that we don't have unreasonable negative values

		if debug: Util.save_data(0,ds,bg_1d,"ctf.bgb4.txt")
		
		# Fit the CTF parameters
		if debug : print "Fit CTF"
		if options.dfmax != None : dfhint=(.15,options.dfmax)
		elif options.curdefocushint :
			try:
				ctf=EMAN2Ctf()
				ctf.from_string(db_parms[name][0])
				curdf=ctf.defocus
				dfhint=(curdf-0.4,curdf+0.4)
				print "Using existing defocus as hint :",dfhint
			except : 
				dfhint=None
				print "No existing defocus to start with"
		else: dfhint=None
		ctf=ctf_fit(im_1d,bg_1d,bg_1d_low,im_2d,bg_2d,options.voltage,options.cs,options.ac,apix,bgadj=not options.nosmooth,autohp=options.autohp,dfhint=dfhint,verbose=options.verbose)
		try : qual=db_parms[name][3]
		except: qual=5
		db_parms[name]=[ctf.to_string(),im_1d,bg_1d,qual]		# the 5 is a default quality value
		db_im2d[name]=im_2d
		db_bg2d[name]=bg_2d
		
		if debug:
			Util.save_data(0,ds,im_1d,"ctf.fg.txt")
			Util.save_data(0,ds,bg_1d,"ctf.bg.txt")
			Util.save_data(0,ds,ctf.snr,"ctf.snr.txt")
			
		img_sets.append([filename,ctf,im_1d,bg_1d,im_2d,bg_2d,5,bg_1d_low])
		if logid : E2progress(logid,float(i+1)/len(options.filenames))
		
	project_db = db_open_dict("bdb:project")
	project_db["global.microscope_voltage"] = options.voltage
	project_db["global.microscope_cs"] = options.cs
	project_db["global.apix"] = apix
	
	#db_close_dict("bdb:project")
	#db_close_dict("bdb:e2ctf.parms")
	
	return img_sets

def refine_and_smoothsnr(options,strfact,debug=False):
	"""This will refine already determined defocus values by maximizing high-resolution smoothed
	SNR and also replace the stored SNR curve with the smoothed version. Note that this REQUIRES
	that a structure factor be loaded and present"""
	global logid
	img_sets=[]
	db_parms=db_open_dict("bdb:e2ctf.parms")
	db_parms_old=db_open_dict("bdb:e2ctf.parms.presmooth")

	olddf=[]
	newdf=[]
	skipped=0
	for i,filename in enumerate(options.filenames):
		name=get_file_tag(filename)

		if debug : print "Processing ",filename
		
		ctf=EMAN2Ctf()
		orig=db_parms[name]

		#if db_parms_old.has_key(name) :
			#if options.verbose: print filename," already smoothed. Skipping."
			#skipped+=1
			#continue
		
		db_parms_old[name]=orig
		try: ctf.from_string(orig[0])
		except :
			print "no old data for %s. Skipping"%name
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
		db_parms[name]=[ctf.to_string()]+orig[1:]
		
		if logid : E2progress(logid,float(i+1)/len(options.filenames))
		
	if skipped>0 : print "Warning: %d files skipped, because they were already smoothed"%skipped
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
	
def process_stack(stackfile,phaseflip=None,phasehp=None,wiener=None,edgenorm=True,oversamp=1,default_ctf=None,invert=False,virtualout=None,storeparm=False,source_image=None):
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
	db_parms=db_open_dict("bdb:e2ctf.parms")

	if virtualout: 
		vin=db_open_dict(stackfile)
		vout=db_open_dict(virtualout)
	
	name = get_file_tag(stackfile)
	if phasehp:
		p1d=db_parms[name][1]
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

	
	for i in range(n):
		if source_image!=None :
			im1=EMData()
			im1.read_image(stackfile,i,True)
			try:
				if im1["ptcl_source_image"]!=source_image : continue
			except:
				print "Image %d doesn't have the ptcl_source_image parameter. Skipping."%i
				continue
		
		im1 = EMData(stackfile,i)
		try: ctf=im1["ctf"]
		except : ctf=default_ctf
		if storeparm : 
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
			
		if phaseflip or phasehp:
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
			out.process("normalize.edgemean")
			if phaseflip: out.write_image(phaseflip,i)

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

				out.process_inplace("normalize.edgemean")
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

		if virtualout:
			im1["data_path"]=vin.get_data_path(i)
			vout[vout["maxrec"]+1]=im1

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
	
	db_close_dict(stackfile)		# safe for non bdb urls
	return av

masks={}		# mask cache for background/foreground masking
def powspec_with_bg(stackfile,source_image=None,radius=0,edgenorm=True,oversamp=1,apix=2):
	"""This routine will read the images from the specified file, optionally edgenormalize,
	then apply a gaussian mask with the specified radius then compute the average 2-D power 
	spectrum for the stack. It will also compute the average 2-D power spectrum using 1-mask + edge 
	apotization to get an appoximate 'background' power spectrum. 2-D results returned as a 2-D FFT 
	intensity/0 image. 1-D results returned as a list of floats. If source_image is provided, then it
	will only read particles with ptcl_source_image set to the specified value.
	
	returns a 5-tuple with spectra for (1d particle,1d background,2d particle,2d background,1d background non-convex)
	"""
	
	global masks
	
	im = EMData()
	im.read_image(stackfile,0)
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
			im1.clip_inplace(Region(-(ys2*(oversamp-1)/2),-(ys2*(oversamp-1)/2),ys,ys))
		
		im2=im1.copy()

#		print im2.get_size(), im1.get_size()

		im1*=mask1
		imf=im1.do_fft()
		imf.ri2inten()
		if i==0: av1=imf
		else: av1+=imf
	
		im2*=mask2
		imf=im2.do_fft()
		imf.ri2inten()
		if i==0: av2=imf
		else: av2+=imf
		
	
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

	# added to make the 2D display look better. Should have no other impact at the time
	# it's being added, though autofitting may rely on it in future, so it shouldn't be removed --steve (8/3/11)
	av1.process_inplace("math.sqrt")
	av1["is_intensity"]=0

	av2.process_inplace("math.sqrt")
	av2["is_intensity"]=0

	# This is a new addition (2/4/10) to prevent negative BG subtracted curves near the origin
	maxpix=int(apix*ys2/80.0)		# we do this up to ~80 A
	for i in xrange(maxpix) :
		if av2_1d[i]>av1_1d[i] : av2_1d[i]=av1_1d[i]		# we set the background equal to the foreground if it's too big

	db_close_dict(stackfile)	# safe for non-bdb files
	
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
	
	db_close_dict(stackfile)
	return av

def smooth_bg(curve,ds):
	"""Smooths a background curve by doing a running average of the log of the curve, ignoring the first few points"""
	
	first=int(.02/ds)	# start at 1/50 1/A
	if first<2 : first=2

	return curve[:first]+[pow(curve[i-1]*curve[i]*curve[i+1],.33333) for i in range(first,len(curve)-2)]+[curve[-2],curve[-1]]
#	return curve[:first]+[pow(curve[i-2]*curve[i-1]*curve[i]*curve[i+1]*curve[i+2],.2) for i in range(first,len(curve)-2)]+[curve[-2],curve[-1]]

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
	if s<=0 or n<=0 : return 0.0
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
	
	ret=bg_1d[:]
	
	# force the curve to be non-convex, after a point
	cont=True
	while cont:
		cont=False
		for i in range(3,len(ret)-1):
			if ret[i]*2.0>ret[i-1]+ret[i+1] : 
				ret[i]=(ret[i-1]+ret[i+1])/2.0
				cont=True
	
	return ret

def ctf_fit(im_1d,bg_1d,bg_1d_low,im_2d,bg_2d,voltage,cs,ac,apix,bgadj=0,autohp=False,dfhint=None,verbose=1):
	"""Determines CTF parameters given power spectra produced by powspec_with_bg()
	The bgadj option will result in adjusting the bg_1d curve to better match the zeroes
	of the CTF (in which case bg_1d is modified in place)."""
	# defocus estimation
	global debug
	global sfcurve
#	from bisect import insort_left
	
	ys=im_2d.get_ysize()
	ds=1.0/(apix*ys)
	if ac<0 or ac>100 :
		print "Invalid %%AC, defaulting to 10"
		ac=10.0
	
	ctf=EMAN2Ctf()
	ctf.from_dict({"defocus":1.0,"voltage":voltage,"bfactor":500.0,"cs":cs,"ampcont":ac,"apix":apix,"dsbg":ds,"background":bg_1d})
	
	sf = [sfact(i*ds) for i in range(ys)]

	bgsub=[im_1d[s]-bg_1d[s] for s in range(len(im_1d))]	# background subtracted curve, good until we readjust the background
	bglowsub=[im_1d[s]-bg_1d_low[s] for s in range(len(im_1d))]	# background subtracted curve, using a non-convex version of the background 

	s1=min(int(.167/ds),ys/3-4)

#	plot(bglowsub)
	# We find the first minimum (also <1/50 max) after the value has fallen to 1/5 max
	maxsub=max(bglowsub[4:])
	for i in range(len(bglowsub)-1,0,-1):
		if bglowsub[i]>maxsub/5.0 : break
	s0=max(int(.005/ds),i)
	while bglowsub[s0]>bglowsub[s0+1] or bglowsub[s0+1]>maxsub/50.0: s0+=1	# look for a minimum in the data curve
	if verbose: print "Minimum at 1/%1.1f 1/A (%1.4f), highest s considered 1/%1.1f 1/A (%1.4f)"%(1.0/(s0*ds),s0*ds,1.0/(s1*ds),s1*ds)
	
	# we now have a range from the first minimia to ~6 A. Adjust it to include a bit of data closer to the origin
	for i in range(s0,s0*2/3,-1):
		if bglowsub[i]>bglowsub[i-1] : break
	s0m=s0									# save for a defocus estimate
	s0=i
	
	if verbose: print "s0 adjusted to 1/%1.1f 1/A"%(1.0/(s0*ds))

	if s1<=s0 :
		print "Error: couldn't fit this data set due to inability to locate appropriate minima"
		return ctf	# This way we can still proceed to the GUI... a bit of a bad hack...
	
	if debug:
		dfout=file("ctf.df.txt","w")

#	dfbest1=(0,-1.0e20)
	dfbest1=[]
	
	# This implies that we already have a general range for the defocus, so we limit the search
	if isinstance(dfhint,float) :
		dfhint=int(dfhint*20.0)
		rng=range(dfhint-3,dfhint+4)
		if verbose: print "Using defocus hint %1.2f (%1.2f - %1.2f)"%(dfhint/20.0,rng[0]/20.0,rng[-1]/20.0)
	elif isinstance(dfhint,tuple) :
		rng=(int(dfhint[0]*20.0),int(dfhint[1]*20.0))
		if verbose: print "Using limited defocus range (%1.2f - %1.2f)"%(rng[0]/20.0,rng[-1]/20.0)
	else:
		rng=range(3,128)
		if verbose: print "(%s) Using full defocus range (%1.2f - %1.2f)"%(str(dfhint),rng[0]/20.0,rng[-1]/20.0)


	# This loop tries to find the best few possible defocuses
	for dfi in rng:			# loop over defocus
		df=dfi/20.0
		ctf.defocus=df
#		ctf.ampcont=ac
		cc=ctf.compute_1d(ys,ds,Ctf.CtfType.CTF_AMP)
		norm=0
		for fz in range(len(cc)): 
			if cc[fz]<0 : break
	
		tot,tota,totb=0,0,0
		zroa,zrob=0,0
		for s in range(s0,s1):
			# This is basicaly a dot product of the reference curve vs the simulation
			a=(cc[s]**2)				# ctf curve ^2
			b=bglowsub[s]					# bg subtracted intensity
			tot+=a*b
			tota+=a*a
			totb+=b*b
			
			# This computes the mean value near zeroes
			#if a<.1 and s<ys/3: 
				#zroa+=(im_1d[s]-bg_1d[s])*s
				#zrob+=1.0

		tot/=sqrt(tota*totb)	# correct if we want a normalized dot product
#		tot/=tota				# funnny 'normalization' which seems to help bring out the more likely peaks
#		if zrob==0 : continue
#		tot*=-(zroa/zrob)
		
		dfbest1.append((tot,df))		# we keep all of the results, then process them at the end
#		if tot>dfbest1[1] : dfbest1=(df,tot)
		try :dfout.write("%1.2f\t%g\n"%(df,tot))
		except : pass
	
	# Now we find the best few defocus choices
	dfbest1a=[]						# keep only peaks, and sort them
	for i in range(1,len(dfbest1)-1):
		if dfbest1[i]>dfbest1[i-1] and dfbest1[i]>dfbest1[i+1] : dfbest1a.append(dfbest1[i])		# keep only local peaks
	
	if len(dfbest1a)==0 : dfbest1a=[max(dfbest1)]
	
	dfbest1a.sort()
	dfbest1a=dfbest1a[-5:]		# keep at most the best 5 peaks
	
	if verbose: print "Initial defocus possibilities: ",
	for i in dfbest1a: print i[1],
	print

	if len(dfbest1a)==1 :
		best=[[0.0,[dfbest1a[0][1],300.0]]]
	else :
		# Next, we use a simplex minimizer to try for the best CTF parameters for each defocus
		best=[]
		for b1a in dfbest1a:
			# our parameter set is (defocus,bfactor)
			parm=[b1a[1],500.0]

	#		print "Initial guess : ",parm
			sim=Simplex(ctf_cmp_2,parm,[.02,20.0],data=(ctf,bglowsub,s0,s1,ds,parm[0],rng))
			oparm=sim.minimize(epsilon=.0000001,monitor=0)
	#		print "Optimized : ",oparm
			best.append((oparm[1],oparm[0]))			


		# best now contains quality,(df,bfac) for each optimized answer
		best.sort()
#		if verbose: print "Best value is df=%1.3f  B=%1.1f"%(best[0][1][0],best[0][1][1])

	ctf.defocus=best[0][1][0]
	ctf.bfactor=best[0][1][1]
#	cc=ctf.compute_1d(ys,ds,Ctf.CtfType.CTF_AMP)
#	Util.save_data(0,ds,cc,"ctf.ctf.txt")
#	print best[0]

	if bgadj:
		bg2=bg_1d[:]
		if verbose: print "BG adjustment using df=%1.3f  B=%1.1f"%(best[0][1][0],best[0][1][1])
		for i in range(6):
			# now we try to construct a better background based on the CTF zeroes being zero
			df=best[0][1][0]
			ctf.defocus=best[0][1][0]
			ctf.bfactor=best[0][1][1]
			cc=ctf.compute_1d(ys,ds,Ctf.CtfType.CTF_AMP)
#			bg2=bg_1d[:]
			last=0,1.0
			for x in range(1,len(bg2)-1) : 
				if cc[x]*cc[x+1]<0 :
					if x<len(bg2)/2 :
						# Below 1/2 nyquist, we work harder to avoid negative values
						cur=(x,min(im_1d[x]/bg2[x],im_1d[x+1]/bg2[x+1],im_1d[x-1]/bg2[x-1],im_1d[x+2]/bg2[x+2]))						
					else:
						# we search the two points 'at' the zero for the minimum
						# if we're too aggressive about this, we will end up exaggerating the high
						# resolution SNR
						cur=(x,min(im_1d[x]/bg2[x],im_1d[x+1]/bg2[x+1]))

					# once we have a pair of zeros we adjust the background values between
					for xx in range(last[0],cur[0]):
						w=(xx-last[0])/float(cur[0]-last[0])
						bg_1d[xx]=bg2[xx]*(cur[1]*w+last[1]*(1.0-w))

					last=cur
			# cover the area from the last zero crossing to the end of the curve
			for xx in range(last[0],len(bg2)):
				bg_1d[xx]=bg2[xx]*last[1]

	#	s0=int(.04/ds)+1
	#	s1=min(int(0.15/ds),len(bg_1d)-1)

			# rerun the simplex with the new background
			bgsub=[im_1d[s]-bg_1d[s] for s in range(len(im_1d))]
			bg_1d_low=low_bg_curve(bg_1d,ds)
			bglowsub=[im_1d[s]-bg_1d_low[s] for s in range(len(im_1d))]	# background subtracted curve, using a non-convex version of the background 
			best[0][1][1]=500.0		# restart the fit with B=200.0
			if hasgoodsf: sim=Simplex(ctf_cmp_a,best[0][1],[.02,20.0],data=(ctf,bgsub,s0,s1,ds,best[0][1][0],rng))
			else: sim=Simplex(ctf_cmp,best[0][1],[.02,20.0],data=(ctf,bgsub,s0,s1,ds,best[0][1][0],rng))
			oparm=sim.minimize(epsilon=.00000001,monitor=0)
			if fabs(df-oparm[0][0])/oparm[0][0]<.001:
				best[0]=(oparm[1],oparm[0])
				break
			best[0]=(oparm[1],oparm[0])
			if verbose: print "After BG correction, value is df=%1.3f  B=%1.1f"%(best[0][1][0],best[0][1][1])
	else:
		if verbose: print "Best value is df=%1.3f  B=%1.1f"%(best[0][1][0],best[0][1][1])
		# rerun the simplex with the new background
		best[0][1][1]=500.0		# restart the fit with B=200.0
		if hasgoodsf : sim=Simplex(ctf_cmp_a,best[0][1],[.02,20.0],data=(ctf,bgsub,s0,s1,ds,best[0][1][0],rng))
		else : sim=Simplex(ctf_cmp,best[0][1],[.02,20.0],data=(ctf,bgsub,s0,s1,ds,best[0][1][0],rng))
		oparm=sim.minimize(epsilon=.0000001,monitor=0)
		best[0]=(oparm[1],oparm[0])
		if verbose: print "After BG correction, value is df=%1.3f  B=%1.1f"%(best[0][1][0],best[0][1][1])

	snr=[snr_safe(im_1d[i],bg_1d[i]) for i in range(len(im_1d))]
	
	# This will dramatically reduce the intensity of the initial sharp peak found in almost all single particle data
	# this applies to the SNR curve only, downweighting the importance of this section of the spectrum without actually
	# removing the information by filtering the image data. It will, of course also impact Wiener filters.
	if autohp:
		for x in range(2,len(snr)-2):
			if snr[x]>snr[x+1] and snr[x+1]<snr[x+2] : break	# we find the first minimum
		
		snr1max=max(snr[1:x])				# find the intensity of the first peak
		snr2max=max(snr[x+2:len(snr)/2])		# find the next highest snr peak

		for xx in range(1,x+1): snr[xx]*=0.5*snr2max/snr1max		# scale the initial peak to 50% of the next highest peak

	
	# store the final results
	ctf.snr=snr
	ctf.defocus=best[0][1][0]
	ctf.bfactor=best[0][1][1]

	if verbose : print "Best DF = %1.3f   B-factor = %1.0f   avSNR = %1.4f"%(ctf.defocus,ctf.bfactor,(sum(ctf.snr)/float(len(ctf.snr))))
	

	Util.save_data(0,ds,ctf.snr,"ctf.snr")
	Util.save_data(0,ds,bg_1d,"ctf.bg1d")

	return ctf

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

	er*=max(1.0,1.0+parms[0]*20.0-rng[1])**2		# penalty for being outside range (high)
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
		v=sfact(i*ds)*cc[i]*cc[i]
		a+=v*bgsub[i]
		b+=v*v
		c+=bgsub[i]*bgsub[i]
#		if v>.001 : out.write("%f\t%f\n"%(i*ds,bgsub[i]/v))
#	print s0,s1,a,b,c
	er=1.0-a/sqrt(b*c)


#	print er,(parms[0]-dforig)**2,parms[1]
	
	er*=(1.0+300.0*(parms[0]-dforig)**4)		# This is a weight which biases the defocus towards the initial value
	er*=max(1.0,1.0+parms[0]*20.0-rng[1])**2		# penalty for being outside range (high)
	er*=max(1.0,1.0+rng[0]-parms[0]*20.0)**2		# penalty for being outside range (low)
	er*=1.0+(parms[1]-200)/20000.0+exp(-(parms[1]-50.0)/30.0)		# This is a bias towards small B-factors and to prevent negative B-factors
#	er*=(1.0+fabs(parms[1]-200.0)/100000.0)		# This is a bias towards small B-factors
	
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
	from emshape import *
	from valslider import ValSlider
except:
	print "Warning: PyQt4 must be installed to use the --gui option"
	class dummy:
		pass
	class QWidget:
		"A dummy class for use when Qt not installed"
		def __init__(self,parent):
			print "Qt4 has not been loaded"
	QtGui=dummy()
	QtGui.QWidget=QWidget
		
def notzero(x):
	if x==0 : return 1.0
	return x
		
class GUIctf(QtGui.QWidget):
	def __init__(self,application,data,autohp=True,nosmooth=False):
		"""Implements the CTF fitting dialog using various EMImage and EMPlot2D widgets
		'data' is a list of (filename,ctf,im_1d,bg_1d,im_2d,bg_2d)
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
		
		QtGui.QWidget.__init__(self,None)
		self.setWindowIcon(QtGui.QIcon(get_image_directory() + "ctf.png"))
		
		self.data=data
		self.curset=0
		self.plotmode=0
		
		self.guiim=EMImage2DWidget(application=self.app())
		self.guiiminit = True # a flag that's used to auto resize the first time the gui's set_data function is called
		self.guiplot=EMPlot2DWidget(application=self.app())
		self.guirealim=EMImage2DWidget(application=self.app())	# This will show the original particle images
		
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
		self.setlist=QtGui.QListWidget(self)
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
		
		self.squality=ValSlider(self,(0,10),"Quality (0-10):",0,90)
		self.squality.setIntonly(True)
		self.vbl.addWidget(self.squality)
		
		
		self.hbl_buttons = QtGui.QHBoxLayout()
		self.saveparms = QtGui.QPushButton("Save parms")
		self.recallparms = QtGui.QPushButton("Recall")
		self.refit = QtGui.QPushButton("Refit")
		self.output = QtGui.QPushButton("Output")
		self.hbl_buttons.addWidget(self.refit)
		self.hbl_buttons.addWidget(self.saveparms)
		self.hbl_buttons.addWidget(self.recallparms)
		self.hbl_buttons2 = QtGui.QHBoxLayout()
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
		QtCore.QObject.connect(self.setlist,QtCore.SIGNAL("currentRowChanged(int)"),self.newSet)
		QtCore.QObject.connect(self.splotmode,QtCore.SIGNAL("currentIndexChanged(int)"),self.newPlotMode)

	   	QtCore.QObject.connect(self.saveparms,QtCore.SIGNAL("clicked(bool)"),self.on_save_params)
		QtCore.QObject.connect(self.recallparms,QtCore.SIGNAL("clicked(bool)"),self.on_recall_params)
		QtCore.QObject.connect(self.refit,QtCore.SIGNAL("clicked(bool)"),self.on_refit)
		QtCore.QObject.connect(self.output,QtCore.SIGNAL("clicked(bool)"),self.on_output)
		
		self.update_data()
		
		self.resize(720,380) # figured these values out by printing the width and height in resize event
		
		
		self.setWindowTitle("CTF")

	def on_save_params(self):
		
		if len(self.setlist.selectedItems()) == 0: return
			
		val = self.curset
		name = str(self.setlist.item(val).text())
		name = get_file_tag(name)
		
#		if not db_check_dict(name):
#			print "error, the db doesn't exist:",name
#			
		db_parms=db_open_dict("bdb:e2ctf.parms")
		ctf = self.data[val][1].to_string()

		tmp=db_parms[name]
		tmp[0]=ctf
		tmp[3]=self.data[val][6]
		db_parms[name] = tmp

	def on_recall_params(self):
		if len(self.setlist.selectedItems()) == 0: return
			
		val = self.curset
		name = str(self.setlist.item(val).text())
		name = get_file_tag(name)
		
		db_parms=db_open_dict("bdb:e2ctf.parms")
		if not db_parms.has_key(name):
			print "error, ctf parameters do not exist for:",name
		
		ctf=EMAN2Ctf()
		ctf.from_string(db_parms[name][0])
		
		self.data[val][1]=ctf
		self.data[val][6]=db_parms[name][3]
		self.newSet(self.curset)
	
#	def get_output_params(self):

	def on_refit(self):
		# self.data[n] contains filename,EMAN2CTF,im_1d,bg_1d,im_2d,bg_2d,qual
		tmp=list(self.data[self.curset])
		ctf=ctf_fit(tmp[2],tmp[3],tmp[7],tmp[4],tmp[5],tmp[1].voltage,tmp[1].cs,tmp[1].ampcont,tmp[1].apix,bgadj=not self.nosmooth,autohp=self.autohp,dfhint=self.sdefocus.value)

		# refit will also write parameters back to the db
		db_parms=db_open_dict("bdb:e2ctf.parms")
		name=get_file_tag(tmp[0])
		tmp=db_parms[name]
		tmp[0]=ctf.to_string()		# the 5 is a default quality value
		db_parms[name]=tmp

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
		if self.guiim != None:
			self.app().close_specific(self.guiim)
			self.guiim = None 
		if self.guiplot != None:
			self.app().close_specific(self.guiplot)
		if self.guirealim != None:
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
		ds=self.data[val][1].dsbg
		r=len(ctf.background)
		s=[ds*i for i in range(r)]
		
		if r==0 or len(s)<2 :
			print "Trying to plot bad data (set %s): %s"%(str(val),str(ctf))
			return
		
		# This updates the image circles
		fit=ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP)
		shp={}
		nz=0
		for i in range(1,len(fit)):
			if fit[i-1]*fit[i]<=0.0: 
				nz+=1
				shp["z%d"%i]=EMShape(("circle",0.0,0.0,1.0/nz,r,r,i,1.0))
#				if nz==1: print ("circle",0.0,0.0,1.0/nz,r,r,i,1.0)
		
		self.guiim.del_shapes()
		self.guiim.add_shapes(shp)
		self.guiim.updateGL()
		
		if self.plotmode==1:
			self.guiplot.set_data((s,self.data[val][2]),"fg",True,True,color=1)
			self.guiplot.set_data((s,self.data[val][7]),"bg(concave)",quiet=True,color=0,linetype=2)
			self.guiplot.set_data((s,self.data[val][3]),"bg",color=0)
			self.guiplot.setAxisParms("s (1/"+ u"\u212B" +")","Intensity (a.u)")
		elif self.plotmode==0: 
			bgsub=[self.data[val][2][i]-self.data[val][3][i] for i in range(len(self.data[val][2]))]
			self.guiplot.set_data((s,bgsub),"fg-bg",True,True,color=0)
			
			lowbg=self.data[val][7]
			lowbgsub=[self.data[val][2][i]-lowbg[i] for i in range(len(self.data[val][2]))]
			self.guiplot.set_data((s,lowbgsub),"fg-lowbg",False,True,color=0,linetype=2)

			
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
			self.guiplot.setAxisParms("s (1/"+ u"\u212B" + ")","Intensity (a.u)")
		elif self.plotmode==2:
			snr=ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_SNR)		# The snr curve
			self.guiplot.set_data((s,snr[:len(s)]),"snr",True)
			self.guiplot.setAxisParms("s (1/"+ u"\u212B" +")","SNR (intensity ratio)")
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

			self.guiplot.setAxisParms("s (1/"+ u"\u212B" +")","SNR (intensity ratio)")
		elif self.plotmode==4:
			snr=ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_SNR)		# The snr curve
			for i in range(1,len(snr)): snr[i]=snr[i]*i+snr[i-1]			# integrate SNR*s
#			for i in range(1,len(snr)): snr[i]/=snr[-1]				# normalize
			for i in range(1,len(snr)): snr[i]/=len(snr)			# this way the relative quality of images can be compared
			self.guiplot.set_data((s,snr[:len(s)]),"snr",True)
			self.guiplot.setAxisParms("s (1/"+ u"\u212B" +")","Integrated SNR")
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
			self.guiplot.setAxisParms("s (1/"+ u"\u212B" +")","CTF Sum")
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

			self.guiplot.setAxisParms("s (1/"+ u"\u212B" +")","SNR Sum")
		# All SNR
		elif self.plotmode==7:
			self.guiplot.set_data(None,None,True,True)		# erase existing data quietly
			for st in self.data:
				inten=st[1].compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_SNR)
				inten=[(s[i],inten[i]) for i in xrange(1,len(inten)-1) if inten[i]>inten[i+1] and inten[i]>inten[i-1]]
				ls,li=zip(*inten)	# this confusing idiom is unzipping the list
#				print len(ls),len(li)
				self.guiplot.set_data((ls,li),st[0],quiet=True,linetype=-1,symtype=0,symsize=2)
				
			self.guiplot.setAxisParms("s (1/"+ u"\u212B" +")","SNR")
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
			
			self.guiplot.set_data((s,self.data[val][2]),"fg",True,True,color=1)
			self.guiplot.set_data((s,self.data[val][3]),"bg",color=0)
			self.guiplot.set_data((s,low_bg_curve(self.data[val][3],ds)),"lowbg",color=2)
			self.guiplot.setAxisParms("s (1/"+ u"\u212B" +")","Intensity (a.u)")
			
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
			self.guiim.set_data(self.data[val][4])
			if self.guiiminit:
				self.guiim.optimally_resize()
				self.guiiminit = False
			self.guiim.updateGL()
		self.update_plot()
		
		self.ptcldata=EMData.read_images(self.data[val][0],range(0,20))
		self.guirealim.set_data(self.ptcldata)

	def newPlotMode(self,mode):
		self.plotmode=mode
		self.update_plot()

	def newCTF(self) :
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
		name = str(self.setlist.item(val).text())
		name = get_file_tag(name)

		db_parms=db_open_dict("bdb:e2ctf.parms")
		tmp=db_parms[name]
		tmp[3]=self.data[val][6]
		db_parms[name] = tmp

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
