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
# This is a program for determining CTF parameters

import global_def
from global_def import *

from EMAN2 import *
from sparx import *
from optparse import OptionParser
from OpenGL import GL,GLUT
from numpy import *
import time
import os
import sys

from Simplex import Simplex

debug=False
logid=None

sfcurve=None		# This will store a global structure factor curve if specified
envelopes=[]		# simplex minimizer needs to use a global at the moment
import weakref

def main():
	global debug,logid
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <input stack/image> ...
	
Various CTF-related operations on images, including automatic fitting. Note that automatic fitting is limited to 5 microns
underfocus at most. Input particles should be unmasked and unfiltered. A minimum of ~20% padding around the
particles is required for background extraction, even if this brings the edge of another particle into the box in some cases.
Particles should be reasonably well centered. Can also optionally phase flip and Wiener filter particles. Wiener filtration comes
after phase-flipping, so if phase flipping is performed Wiener filtered particles will also be phase-flipped. Note that both
operations are performed on oversampled images if specified (though final real-space images are clipped back to their original
size. Increasing padding during the particle picking process will improve the accuracy of phase-flipping, particularly for
images far from focus."""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--gui",action="store_true",help="Start the GUI for interactive fitting",default=False)
	parser.add_option("--auto_fit",action="store_true",help="Runs automated CTF fitting on the input images",default=False)
	parser.add_option("--bgmask",type="int",help="Compute the background power spectrum from the edge of the image, specify a mask radius in pixels which would largely mask out the particles. Default is boxsize/2.",default=0)
	parser.add_option("--apix",type="float",help="Angstroms per pixel for all images",default=0)
	parser.add_option("--voltage",type="float",help="Microscope voltage in KV",default=0)
	parser.add_option("--cs",type="float",help="Microscope Cs (spherical aberation)",default=0)
	parser.add_option("--ac",type="float",help="Amplitude contrast (percentage, default=10)",default=10)
	parser.add_option("--autohp",action="store_true",help="Automatic high pass filter of the SNR only to remove initial sharp peak, phase-flipped data is not directly affected (default false)",default=False)
	#parser.add_option("--invert",action="store_true",help="Invert the contrast of the particles in output files (default false)",default=False)
	parser.add_option("--nonorm",action="store_true",help="Suppress per image real-space normalization",default=False)
	parser.add_option("--nosmooth",action="store_true",help="Disable smoothing of the background (running-average of the log with adjustment at the zeroes of the CTF)",default=False)
	#parser.add_option("--phaseflip",action="store_true",help="Perform phase flipping after CTF determination and writes to specified file.",default=False)
	#parser.add_option("--wiener",action="store_true",help="Wiener filter (optionally phaseflipped) particles.",default=False)
	parser.add_option("--oversamp",type="int",help="Oversampling factor",default=1)
	parser.add_option("--sf",type="string",help="The name of a file containing a structure factor curve. Can improve B-factor determination.",default=None)
	parser.add_option("--debug",action="store_true",default=False)
	
	(options, args) = parser.parse_args()

	if len(args)<1 : parser.error("Input image required")
	
	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()

	if options.auto_fit:
		if options.voltage==0 : parser.error("Please specify voltage")
		if options.cs==0 : parser.error("Please specify Cs")
	if options.apix==0 : print "Using A/pix from header"
		
	debug=options.debug

	global sfcurve
	if options.sf :
		sfcurve=XYData()
		sfcurve.read_file(options.sf)

	logid=E2init(sys.argv)

#	if options.oversamp>1 : options.apix/=float(options.oversamp)

	db_project=db_open_dict("bdb:project")
	db_parms=db_open_dict("bdb:e2ctf.parms")
	db_misc=db_open_dict("bdb:e2ctf.misc")

	options.filenames = args
	### Power spectrum and CTF fitting
	if options.auto_fit:
		img_sets=pspec_and_ctf_fit(options,debug) # converted to a function so to work with the workflow
		
		### This computes the intensity of the background subtracted power spectrum at each CTF maximum for all sets
		global envelopes # needs to be a global for the Simplex minimizer
		# envelopes is essentially a cache of information that could be useful at later stages of refinement
		# as according to Steven Ludtke
		for i in img_sets:
			envelopes.append(ctf_env_points(i[2],i[3],i[1]))

		# we use a simplex minimizer to try to rescale the individual sets to match as best they can
		scales=[1.0]*len(img_sets)
		if (len(img_sets)>3) :
			incr=[0.2]*len(img_sets)
			simp=Simplex(env_cmp,scales,incr)
			scales=simp.minimize(maxiters=1000)[0]
	#		print scales
			print " "

		# apply the final rescaling
		envelope=[]
		for i in range(len(scales)):
			cur=envelopes[i]
			for j in range(len(cur)):
				envelope.append((cur[j][0],cur[j][1]*scales[i]))

		envelope.sort()
		envelope=[i for i in envelope if i[1]>0]	# filter out all negative peak values

		db_misc=db_open_dict("bdb:e2ctf.misc")
		db_misc["envelope"]=envelope
		#db_close_dict("bdb:e2ctf.misc")

		#out=file("envelope.txt","w")
		#for i in envelope: out.write("%f\t%f\n"%(i[0],i[1]))
		#out.close()

	### GUI - user can update CTF parameters interactively
	if options.gui :
		img_sets = get_gui_arg_img_sets(options.filenames)
		if len(img_sets) == 0:
			E2end(logid)
			exit(1)
		from emapplication import EMStandAloneApplication
		app=EMStandAloneApplication()
		gui=GUIctfModule(app,img_sets)
		gui.show_guis()
		app.exec_()

		print "done execution"

	### Process input files
	#if debug : print "Phase flipping / Wiener filtration"
	# write wiener filtered and/or phase flipped particle data to the local database
	#if options.phaseflip or options.wiener: # only put this if statement here to make the program flow obvious
	#	write_e2ctf_output(options) # converted to a function so to work with the workflow

	E2end(logid)

def get_gui_arg_img_sets(filenames):
	'''
	returns the img_sets list required to intialized the GUI correctly
	'''
	
	img_sets = []
	if db_check_dict("bdb:e2ctf.parms"):
		db_parms=db_open_dict("bdb:e2ctf.parms",ro=True)
	else: return img_sets
	for file in filenames:
		name = get_file_tag(file)
		if not db_parms.has_key(name):
			print "error, you must first run auto fit before running the gui - there are no parameters for",name
			return []
		img_set = db_parms[name]
		ctf=EMAN2Ctf()
		ctf.from_string(img_set[0]) # convert to ctf object seeing as it's a string
		img_set[0] = ctf
		actual = [file]
		actual.extend(img_set)
		img_sets.append(actual)
		
	return img_sets

def write_e2ctf_output(options):
	"write wiener filtered and/or phase flipped particle data to the local database"
	global logid
	
	if options.phaseflip or options.wiener:
		db_parms=db_open_dict("bdb:e2ctf.parms")
		for i,filename in enumerate(options.filenames):
			name=get_file_tag(filename)
			if debug: print "Processing ",filename

			if options.phaseflip: phaseout="bdb:particles#"+name+"_ctf_flip"
			else: phaseout=None
		
			if options.wiener:
				if options.autohp: wienerout="bdb:particles#"+name+"_ctf_wiener_hp"
				else: wienerout="bdb:particles#"+name+"_ctf_wiener"
			else : wienerout=None

			#if options.phaseflip: phaseout=name+"_ctf_flip.hed"
			#else: phaseout=None
		
			#if options.wiener: wienerout=name+"_ctf_wiener.hed"
			#else : wienerout=None
			
			if phaseout : print "Phase image out: ",phaseout,"\t",
			if wienerout : print "Wiener image out: ",wienerout,
			print ""
			ctf=EMAN2Ctf()
			ctf.from_string(db_parms[name][0])
			process_stack(filename,phaseout,wienerout,not options.nonorm,options.oversamp,ctf,invert=options.invert)
			if logid : E2progress(logid,float(i+1)/len(options.filenames))
			
			
	
def pspec_and_ctf_fit(options,debug=False):
	"Power spectrum and CTF fitting"
	global logid
	img_sets=[]

	db_parms=db_open_dict("bdb:e2ctf.parms")

	for i,filename in enumerate(options.filenames):
		name=get_file_tag(filename)

		# compute the power spectra
		if debug : print "Processing ",filename
		apix=options.apix
		if apix<=0 : apix=EMData(filename,0,1)["apix_x"] 
		im_1d,bg_1d,im_2d,bg_2d=powspec_with_bg(filename,radius=options.bgmask,edgenorm=not options.nonorm,oversamp=options.oversamp)
		ds=1.0/(apix*im_2d.get_ysize())
		if not options.nosmooth : bg_1d=smooth_bg(bg_1d,ds)

		Util.save_data(0,ds,bg_1d,"ctf.bgb4.txt")

		# Fit the CTF parameters
		if debug : print "Fit CTF"
		ctf=ctf_fit(im_1d,bg_1d,im_2d,bg_2d,options.voltage,options.cs,options.ac,apix,bgadj=not options.nosmooth,autohp=options.autohp)
		db_parms[name]=[ctf.to_string(),im_1d,bg_1d,im_2d,bg_2d]

		if debug:
			Util.save_data(0,ds,im_1d,"ctf.fg.txt")
			Util.save_data(0,ds,bg_1d,"ctf.bg.txt")
			Util.save_data(0,ds,ctf.snr,"ctf.snr.txt")
			
		img_sets.append((filename,ctf,im_1d,bg_1d,im_2d,bg_2d))
		if logid : E2progress(logid,float(i+1)/len(options.filenames))
		
	project_db = db_open_dict("bdb:project")
	project_db["global.microscope_voltage"] = options.voltage
	project_db["global.microscope_cs"] = options.cs
	project_db["global.apix"] = apix
	
	#db_close_dict("bdb:project")
	#db_close_dict("bdb:e2ctf.parms")
	
	return img_sets

def env_cmp(sca,data):
	global envelopes
	env=envelopes
	total=[]
	for i,ii in enumerate(env):
		for j in ii:
			total.append((j[0],j[1]*sca[i]))
	
	total.sort()
	
	ret=0
	for i in range(2,len(total)-2):
		if total[i][1] :
			ret+=((total[i-2][1]-total[i][1])**2+(total[i-1][1]-total[i][1])**2+(total[i+1][1]-total[i][1])**2+(total[i+2][1]-total[i][1])**2)*i
#			ret+=fabs(total[i-2][1]-total[i][1])+fabs(total[i-1][1]-total[i][1])+fabs(total[i+1][1]-total[i][1])+fabs(total[i+2][1]-total[i][1])

	#ret=0
	#for i in range(1,len(total)):
		#if total[i][1] :
			#ret+=fabs((total[i-1][1]/total[i][1])-1.0)/(total[i][0]-total[i-1][0]+.0005)

	return ret
	
def process_stack(stackfile,phaseflip=None,wiener=None,edgenorm=True,oversamp=1,default_ctf=None,invert=False):
	"""Will phase-flip and/or Wiener filter particles in a file based on their stored CTF parameters.
	phaseflip should be the path for writing the phase-flipped particles
	wiener should be the path for writing the Wiener filtered (and possibly phase-flipped) particles
	oversamp will oversample as part of the processing, ostensibly permitting phase-flipping on a wider range of defocus values
	"""
	
	im=EMData()
	im.read_image(stackfile,0) # can't use the constructor if bdb terminology is being used
	ys=im.get_ysize()*oversamp
	ys2=im.get_ysize()
	n=EMUtil.get_image_count(stackfile)
	lctf=None
	
	
	for i in range(n):
		im1 = EMData()
		im1.read_image(stackfile,i)
		try: ctf=im1["ctf"]
		except : ctf=default_ctf
		if type(ctf)==EMAN1Ctf : ctf=default_ctf	# EMAN1 ctf needs a structure factor for this to work
		
		if edgenorm : im1.process_inplace("normalize.edgemean")
		if oversamp>1 :
			im1.clip_inplace(Region(-(ys2*(oversamp-1)/2),-(ys2*(oversamp-1)/2),ys,ys))
		
		fft1=im1.do_fft()
			
		if phaseflip :
			if not lctf or not lctf.equal(ctf):
				flipim=fft1.copy()
				ctf.compute_2d_complex(flipim,Ctf.CtfType.CTF_SIGN)
			fft1.mult(flipim)
			out=fft1.do_ift()
			out["ctf"]=ctf
			out.clip_inplace(Region(int(ys2*(oversamp-1)/2.0),int(ys2*(oversamp-1)/2.0),ys2,ys2))
			if invert: out.mult(-1.0)
			out.process("normalize.edgemean")
			out.write_image(phaseflip,i)

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
		lctf=ctf

	#if wiener and wiener[:4]=="bdb:" : db_close_dict(wiener)
	#if phaseflip and phaseflip[:4]=="bdb:" : db_close_dict(phaseflip)
	
	return

def powspec(stackfile,mask=None,edgenorm=True,):
	"""This routine will read the images from the specified file, optionally edgenormalize,
	optionally apply a mask then compute the average
	2-D power spectrum for the stack. Results returned as a 2-D FFT intensity/0 image"""
	
	n=EMUtil.get_image_count(stackfile)
	
	for i in range(n):
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
	return av

masks={}		# mask cache for background/foreground masking
def powspec_with_bg(stackfile,radius=0,edgenorm=True,oversamp=1):
	"""This routine will read the images from the specified file, optionally edgenormalize,
	then apply a gaussian mask with the specified radius then compute the average 2-D power 
	spectrum for the stack. It will also compute the average 2-D power spectrum using 1-mask + edge 
	apotization to get an appoximate 'background' power spectrum. 2-D results returned as a 2-D FFT 
	intensity/0 image. 1-D results returned as a list of floats.
	
	returns a 4-tuple with spectra for (1d particle,1d background,2d particle,2d background)
	"""
	
	global masks
	
	im = EMData()
	im.read_image(stackfile,0)
	ys=im.get_ysize()*oversamp
	ys2=im.get_ysize()
	n=EMUtil.get_image_count(stackfile)
	if(radius == 0):  radius = ys2//2 - ys2//4
	
	# set up the inner and outer Gaussian masks
	try:
		mask1,ratio1,mask2,ratio2=masks[(ys,radius)]
	except:
		mask1=EMData(ys2,ys2,1)
		mask1.to_one()
		#mask1.process_inplace("mask.gaussian",{"outer_radius":radius})
		mask1.process_inplace("mask.gaussian",{"inner_radius":radius,"outer_radius":ys2//10})
		mask2=mask1.copy()*-1+1
#		mask1.process_inplace("mask.decayedge2d",{"width":4})
		#mask2.process_inplace("mask.decayedge2d",{"width":4})
		mask1.clip_inplace(Region(-(ys2*(oversamp-1)/2),-(ys2*(oversamp-1)/2),ys,ys))
		mask2.clip_inplace(Region(-(ys2*(oversamp-1)/2),-(ys2*(oversamp-1)/2),ys,ys))
		ratio1=mask1.get_attr("square_sum")/(ys*ys)	#/1.035
		ratio2=mask2.get_attr("square_sum")/(ys*ys)
		masks[(ys,radius)]=(mask1,ratio1,mask2,ratio2)
		print  "  RATIOS  ", ratio1, ratio2,"    ",radius
		#mask1.write_image("mask1.hdf",0)
		#mask2.write_image("mask2.hdf",0)
	pav1 = model_blank(ys2,ys2)
	pva1 = model_blank(ys2,ys2)
	pav2 = model_blank(ys2,ys2)
	pva2 = model_blank(ys2,ys2)
	pf = float(ys2*ys2)**2/4.0
	fofo = []
	for i in range(n):
		im1 = EMData()
		im1.read_image(stackfile,i)
#		im1=EMData(stackfile,i)
		#info(im1)
		#if edgenorm : im1.process_inplace("normalize.edgemean")
		if oversamp>1 :
			im1.clip_inplace(Region(-(ys2*(oversamp-1)/2),-(ys2*(oversamp-1)/2),ys,ys))
		im2=im1.copy()

		im1*=mask1
		pp = periodogram(im1)*pf/ratio1
		fofo.append(rot_avg_table(pp))
		Util.add_img(pav1, pp)
		Util.add_img2(pva1, pp)
		imf=im1.do_fft()
		#print_col(imf,0)
		imf.ri2inten()
		#info(imf)
		imf.set_complex(0)
		#print_col(imf,0)
		#imf.write_image("imf.hdf",0)
		#pop = periodogram(im1)
		#print_col(pop,90)
		#pop.write_image("pop.hdf",0)
		if i==0: av1=imf
		else: av1+=imf
	
		im2*=mask2
		pp = periodogram(im2)*pf/ratio2
		Util.add_img(pav2, pp)
		Util.add_img2(pva2, pp)
		imf=im2.do_fft()
		imf.ri2inten()
		if i==0: av2=imf
		else: av2+=imf
	from sys import exit
	
	av1/=(float(n)*av1.get_ysize()*av1.get_ysize()*ratio1)
	av1.set_value_at(0,0,0.0)
	av1.set_complex(1)
	av1.set_attr("is_intensity", 1)

	av2/=(float(n)*av2.get_ysize()*av2.get_ysize()*ratio2)
	av2.set_value_at(0,0,0.0)
	av2.set_complex(1)
	av2.set_attr("is_intensity", 1)

	av1_1d=av1.calc_radial_dist(av1.get_ysize()/2,0.0,1.0,1)
	av2_1d=av2.calc_radial_dist(av2.get_ysize()/2,0.0,1.0,1)
	pav1 /= n
	pva1 = square_root((pva1 - pav1*pav1*n)/(n-1)/n)
	pav2 /= n
	pva2 = square_root((pva2 - pav2*pav2*n)/(n-1)/n)
	pav1.write_image("av1.hdf",0)
	pav2.write_image("av2.hdf",0)
	pva1.write_image("va1.hdf",0)
	pva2.write_image("va2.hdf",0)
	write_text_file([range(ys2//2+1), rot_avg_table(pav1), rot_avg_table(pva1), rot_avg_table(pav2), rot_avg_table(pva2)],"pwsa.txt")
	#write_text_file([range(ys2//2+1), fofo[0], fofo[1], fofo[2], fofo[3], fofo[4], fofo[5]],"ipwsa.txt")
	return (av1_1d,av2_1d,av1,av2)


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
	return av

def smooth_bg(curve,ds):
	"""Smooths a background curve by doing a running average of the log of the curve, ignoring the first few points"""
	
	first=int(.02/ds)	# start at 1/50 1/A
	if first<2 : first=2

	return curve[:first]+[pow(curve[i-1]*curve[i]*curve[i+1],.33333) for i in range(first,len(curve)-2)]+[curve[-2],curve[-1]]
#	return curve[:first]+[pow(curve[i-2]*curve[i-1]*curve[i]*curve[i+1]*curve[i+2],.2) for i in range(first,len(curve)-2)]+[curve[-2],curve[-1]]

def least_square(data,dolog=0):
	"simple linear regression for y=mx+b on a list of (x,y) points. Use the C routine if you need speed."
	sum,sum_x,sum_y,sum_xx,sum_xy=0,0,0,0,0
	for d in data:
		if dolog : y=log10(d[1])
		else : y=d[1]

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
	return s/n-1.0

def sfact_(ss):
	"""This will return a curve shaped something like the structure factor of a typical protein. It is not designed to be
	highly accurate, but be good enough for approximate B-factor estimation"""
	s = 2*ss
	global sfcurve
	if sfcurve:			# Replace the generic structure factor with a user-provided one
		v=sfcurve.get_yatx(s)
		return max(0.001,v)
	if s<.004 : return 0
	if s>.2934 : s=.2934
	return pow(10.0,3.6717 - 364.58 * s + 15597 * s**2 - 4.0678e+05 * s**3 + 6.7098e+06 * s**4 - 7.0735e+07 * s**5 + 4.7839e+08 * s**6 - 2.0574e+09 * s**7 +5.4288e+09 * s**8 - 8.0065e+09 * s**9 + 5.0518e+09 * s**10)


def sfact(s, specimen = "ribosome", mode = "original"):
	mode = "nono"
	if(specimen == "ribosome"):
		if(mode == "original"):
			cof = poly1d(
			[  1.74299751e+11,  -6.33552567e+11,	 1.03272555e+12  ,-9.96735886e+11  ,
			   6.33473586e+11,  -2.78968406e+11,	  8.72441041e+10 , -1.95433202e+10 ,
			   3.12563862e+09,  -3.52114681e+08,	  2.72836330e+07 , -1.40398308e+06 ,
			   4.58118110e+04,  -9.16255616e+02,	  8.49640066e+00]
			)
		else:
			cof = poly1d(
			[  1.73873409e+11 , -6.31788788e+11 ,  1.02940787e+12 , -9.92979623e+11,
			   6.30617550e+11 , -2.77427902e+11 ,  8.66371847e+10 , -1.93662061e+10,
			   3.08723052e+09 , -3.45958145e+08 ,  2.65619992e+07 , -1.34260182e+06,
			   4.19721112e+04 , -7.35986734e+02 ,  3.27651774e+00]
			)
	elif(specimen == "groel"):
		if(mode == "original"):
			cof = poly1d(
			[  1.94998954e+11,  -7.15063399e+11 ,  1.17442894e+12 , -1.13967832e+12,
			   7.25788915e+11,  -3.18555091e+11 ,  9.84631356e+10 , -2.15159505e+10,
			   3.28840050e+09,  -3.42603807e+08 ,  2.32975708e+07 , -9.69959192e+05,
			   2.32759627e+04,  -3.90343922e+02 ,  3.47891555e+00]
			)
		else:
			cof = poly1d(
			[  1.93687379e+11,  -7.10019987e+11,   1.16564309e+12 , -1.13050123e+12,
			   7.19376677e+11,  -3.15388629e+11,   9.73246889e+10 , -2.12134190e+10,
			   3.22868925e+09,  -3.33878805e+08,   2.23619961e+07 , -8.96805906e+05,
			   1.90387356e+04,  -2.02856142e+02,  -1.79285394e+00]
			)


	return  exp(polyval(cof, s))


def ctf_fit(im_1d,bg_1d,im_2d,bg_2d,voltage,cs,ac,apix,bgadj=0,autohp=False):
	"""Determines CTF parameters given power spectra produced by powspec_with_bg()
	The bgadj option will result in adjusting the bg_1d curve to better match the zeroes
	of the CTF (in which case bg_1d is modified in place)."""
	# defocus estimation
	global debug
	
	ys = im_2d.get_ysize()
	ds = 1.0/(apix*ys)
	
	ctf=EMAN2Ctf()
	ctf.from_dict({"defocus":1.0,"voltage":voltage,"bfactor":0.0,"cs":cs,"ampcont":ac,"apix":apix,"dsbg":ds,"background":bg_1d})
	
	sf = sfact([i*ds for i in xrange(ys)], "ribosome","nono")
	#print "  SF ",sf
	
	if debug: dfout=file("ctf.df.txt","w")
	dfbest1=(0,-1.0e20)
	for dfi in range(5,128):			# loop over defocus
		ac=10
		df=dfi/20.0
		ctf.defocus=df
		ctf.ampcont=ac
		cc=ctf.compute_1d(ys,ds,Ctf.CtfType.CTF_AMP)
		st=.04/ds
		norm=0
		for fz in range(len(cc)): 
			if cc[fz]<0 : break
	
		tot,totr=0,0
		for s in range(int(st),ys/2): 
			tot+=(cc[s]**2)*(im_1d[s]-bg_1d[s])
			totr+=cc[s]**4
		#for s in range(int(ys/2)): tot+=(cc[s*ctf.CTFOS]**2)*ps1d[-1][s]/norm
		#for s in range(int(fz/ctf.CTFOS),ys/2): tot+=(cc[s*ctf.CTFOS]**2)*ps1d[-1][s]
		#for s in range(int(fz/ctf.CTFOS),ys/2): tot+=(cc[s*ctf.CTFOS]**2)*snr[s]
		tot/=sqrt(totr)
		#tot/=totr
		if tot>dfbest1[1] : dfbest1=(df,tot)
		try :dfout.write("%1.2f\t%g\n"%(df,tot))
		except : pass
	
	
	
	
	#out=file("bg1d2.txt","w")
	#for a,b in enumerate(bg2): out.write("%1.4f\t%1.5f\n"%(a*ds,b))
	#out.close()

	dfbest=dfbest1
	for dfi in range(-10,10):			# loop over defocus
		df=dfi/100.0+dfbest1[0]
		ctf.defocus=df
		cc=ctf.compute_1d(ys,ds,Ctf.CtfType.CTF_AMP)
		st=.04/ds
		norm=0
		for fz in range(len(cc)): 
			#norm+=cc[fz]**2
			if cc[fz]<0 : break

		tot,totr=0,0
		for s in range(int(st),ys/2): 
			tot+=(cc[s]**2)*(im_1d[s]-bg_1d[s])
			totr+=cc[s]**4

		tot/=sqrt(totr)
		if tot>dfbest[1] : 
			dfbest=(df,tot)
		if debug : dfout.write("%1.2f\t%g\n"%(df,tot))

	ctf.defocus=dfbest[0]
	cc=ctf.compute_1d(ys,ds,Ctf.CtfType.CTF_AMP)
	Util.save_data(0,ds,cc,"ctf.ctf.txt")

	if bgadj:
		# now we try to construct a better background based on the CTF zeroes being zero
		bg2=bg_1d[:]
		last=0,1.0
		for x in xrange(1,len(bg2)-1) : 
			if cc[x]*cc[x+1]<0 :
				# we search +-1 point from the zero for the minimum
				cur=(x,min(im_1d[x]/bg_1d[x],im_1d[x-1]/bg_1d[x-1],im_1d[x+1]/bg_1d[x+1]))
				# once we have a pair of zeros we adjust the background values between
				for xx in xrange(last[0],cur[0]):
					w=(xx-last[0])/float(cur[0]-last[0])
					bg_1d[xx]=bg2[xx]*(cur[1]*w+last[1]*(1.0-w))
#					print xx,"\t",(cur[1]*w+last[1]*(1.0-w)) #,"\t",cur[1],last[1]
				last=cur
		# cover the area from the last zero crossing to the end of the curve
		for xx in xrange(last[0],len(bg2)):
			bg_1d[xx]=bg2[xx]*last[1]


	snr=[snr_safe(im_1d[i],bg_1d[i]) for i in xrange(len(im_1d))]
	# This will dramatically reduce the intensity of the initial sharp peak found in almost all single particle data
	# this applies to the SNR curve only, downweighting the importance of this section of the spectrum without actually
	# removing the information by filtering the image data. It will, of course also impact Wiener filters.
	if autohp:
		for x in range(2,len(snr)-2):
			if snr[x]>snr[x+1] and snr[x+1]<snr[x+2] : break	# we find the first minimum

		snr1max=max(snr[1:x])				# find the intensity of the first peak
		snr2max=max(snr[x+2:len(snr)/2])		# find the next highest snr peak
		qtmp = 0.5*snr2max/snr1max
		for xx in range(1,x+1): snr[xx] *= qtmp		# scale the initial peak to 50% of the next highest peak

	# store the final results
	ctf.snr=snr
	ctf.defocus=dfbest[0]
	
	# This smooths the SNR curve
	#snr=ctf.compute_1d(ys,ds,Ctf.CtfType.CTF_SNR_SMOOTH)
	#snr=snr[:len(ctf.snr)]
	#ctf.snr=snr

	# Last, we try to get a decent B-factor
	# This is a quick hack and not very efficiently coded
	bfs=[0.0,50.0,100.0,200.0,400.0,600.0,800.0,1200.0,1800.0,2500.0,4000.0]
	best=(0,0)
	s0=int(.05/ds)+1
	s1=min(int(0.14/ds),len(bg_1d)-1)
	print  "  FREQ RANGE",s0,s1
	for b in range(1,len(bfs)-1):
		ctf.bfactor=bfs[b]
		cc=ctf.compute_1d(ys,ds,Ctf.CtfType.CTF_AMP)
		sf = sfact([ds*i for i in xrange(len(cc))], "ribosome","original")
		cc=[sf[i]*cc[i]**2 for i in xrange(len(cc))]

		# adjust the amplitude to match well
		a0,a1=0,0
		for s in range(s0,s1):
			a0 += cc[s]
			a1 += fabs(im_1d[s]-bg_1d[s])
		if a1==0 : a1=1.0
		a0/=a1
		cc=[i/a0 for i in cc]

		er=0
		# compute the error
		for s in range(s0,len(bg_1d)-1):
			er+=(cc[s]-im_1d[s]+bg_1d[s])**2

		if best[0]==0 or er<best[0] : best=(er,b)

	# Stupid replication here, in a hurry
	bb=best[1]
	best=(best[0],bfs[best[1]])
	for b in xrange(20):
		ctf.bfactor=bfs[bb-1]*(1.0-b/20.0)+bfs[bb+1]*(b/20.0)
		cc=ctf.compute_1d(ys,ds,Ctf.CtfType.CTF_AMP)
		sf = sfact([i*ds for i in xrange(len(cc))], "ribosome","original")
		cc=[sf[i]*cc[i]**2 for i in xrange(len(cc))]
		# adjust the amplitude to match well
		a0,a1=0,0
		for s in range(s0,s1): 
			a0+=cc[s]
			a1+=fabs(im_1d[s]-bg_1d[s])
		if a1==0 : a1=1.0
		a0/=a1
		cc=[i/a0 for i in cc]
		
		er=0
		# compute the error
		for s in xrange(s0,len(bg_1d)-1):
			er+=(cc[s]-im_1d[s]+bg_1d[s])**2

		if best[0]==0 or er<best[0] :
			best=(er,ctf.bfactor)
			print "  BEST  ",best
			#for s in xrange(s0,len(bg_1d)-1):
			#	print  s,cc[s],im_1d[s]-bg_1d[s],sfact(ds*s)
	#for i in xrange(3*len(cc)):
	#	print  ds*i,"   ",sfact(ds*i)

#		print bfs[b],best


	ctf.bfactor=best[1]

	if 1 : print "Best DF = %1.3f   B-factor = %1.0f"%(dfbest[0],ctf.bfactor)


	return ctf

def ctf_env_points(im_1d,bg_1d,ctf) :
	"""This will return a list of x,y points corresponding to the maxima of the ctf in the background
	subtracted power spectrum"""
	ys=len(bg_1d)
	ds=ctf.dsbg
	cc=ctf.compute_1d(ys,ds,Ctf.CtfType.CTF_AMP)
	ret=[]
	
	for i in range(1,len(cc)-1):
		if cc[i-1]<cc[i] and cc[i]>cc[i+1] and im_1d[i]-bg_1d[i]>0 :
			ret.append((i*ds,(im_1d[i]-bg_1d[i])))
#			ret.append((i*ds,(im_1d[i]-bg_1d[i])/sfact(i*ds)))		# this version removes the structure factor (in theory)
		
	return ret

try:
	from PyQt4 import QtCore, QtGui, QtOpenGL
	from PyQt4.QtCore import Qt
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
	
from emapplication import EMQtWidgetModule

class GUIctfModule(EMQtWidgetModule):
	def __init__(self,application,data):
		self.guictf = GUIctf(application,data,self)
		EMQtWidgetModule.__init__(self,self.guictf)
		self.application = weakref.ref(application)
		self.setWindowTitle("CTF")
		
		#GLUT.glutInit([]) # this should be taken care of now

		
	def get_desktop_hint(self):
		return "inspector"
	
	def show_guis(self):
		self.guictf.show_guis()
		
class GUIctf(QtGui.QWidget):
	def __init__(self,application,data,module=None):
		"""Implements the CTF fitting dialog using various EMImage and EMPlot2D widgets
		'data' is a list of (filename,ctf,im_1d,bg_1d,im_2d,bg_2d)
		"""
		try:
			from emimage import EMImageModule
			from emimage2d import EMImage2DModule
		except:
			print "Cannot import EMAN image GUI objects (emimage,etc.)"
			sys.exit(1)
		try: 
			from emplot2d import EMPlot2DModule
		except:
			print "Cannot import EMAN plot GUI objects (is matplotlib installed?)"
			sys.exit(1)
		
		self.app = weakref.ref(application)
		self.module = weakref.ref(module)
		
		QtGui.QWidget.__init__(self,None)
		self.setWindowIcon(QtGui.QIcon(get_image_directory() + "ctf.png"))
		
		self.data=data
		self.curset=0
		self.plotmode=0
		
		self.guiim=EMImage2DModule(application=self.app())
		self.guiplot=EMPlot2DModule(application=self.app())
		
		# FIXME these "emitters" might be unnecessary.
		im_qt_target = self.app().get_qt_emitter(self.guiim)
		plot_qt_target = self.app().get_qt_emitter(self.guiplot)
		
		im_qt_target.connect(im_qt_target,QtCore.SIGNAL("mousedown"),self.imgmousedown)
		im_qt_target.connect(im_qt_target,QtCore.SIGNAL("mousedrag"),self.imgmousedrag)
		im_qt_target.connect(im_qt_target,QtCore.SIGNAL("mouseup")  ,self.imgmouseup)
		plot_qt_target.connect(plot_qt_target,QtCore.SIGNAL("mousedown"),self.plotmousedown)
		
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
		
		self.sdefocus=ValSlider(self,(0,5),"Defocus:",0,90)
		self.vbl.addWidget(self.sdefocus)
		
		self.sbfactor=ValSlider(self,(0,1600),"B factor:",0,90)
		self.vbl.addWidget(self.sbfactor)
		
		self.sampcont=ValSlider(self,(0,100),"% AC",0,90)
		self.vbl.addWidget(self.sampcont)
		
#		self.sapix=ValSlider(self,(.2,10),"A/Pix:",2,90)
#		self.vbl.addWidget(self.sapix)
		
		self.svoltage=ValSlider(self,(0,500),"Voltage (kV):",0,90)
		self.vbl.addWidget(self.svoltage)
		
		self.scs=ValSlider(self,(0,5),"Cs (mm):",0,90)
		self.vbl.addWidget(self.scs)
		self.hbl_buttons = QtGui.QHBoxLayout()
		self.saveparms = QtGui.QPushButton("Save parms")
		self.recallparms = QtGui.QPushButton("Recall")
		self.output = QtGui.QPushButton("Output")
		self.hbl_buttons.addWidget(self.saveparms)
		self.hbl_buttons.addWidget(self.recallparms)
		self.hbl_buttons2 = QtGui.QHBoxLayout()
		self.hbl_buttons2.addWidget(self.output)
		self.vbl.addLayout(self.hbl_buttons)
		self.vbl.addLayout(self.hbl_buttons2)
		
		QtCore.QObject.connect(self.sdefocus, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.sbfactor, QtCore.SIGNAL("valueChanged"), self.newCTF)
#		QtCore.QObject.connect(self.sapix, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.sampcont, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.svoltage, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.scs, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.setlist,QtCore.SIGNAL("currentRowChanged(int)"),self.newSet)
		QtCore.QObject.connect(self.splotmode,QtCore.SIGNAL("currentIndexChanged(int)"),self.newPlotMode)

	   	QtCore.QObject.connect(self.saveparms,QtCore.SIGNAL("clicked(bool)"),self.on_save_params)
		QtCore.QObject.connect(self.recallparms,QtCore.SIGNAL("clicked(bool)"),self.on_recall_params)
		QtCore.QObject.connect(self.output,QtCore.SIGNAL("clicked(bool)"),self.on_output)
		
		self.update_data()
		
		self.update_data()
		self.resize(460,380) # figured these values out by printing the width and height in resize event
		
#	def resizeEvent(self,event):
#		print self.width(),self.height()
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
		output = []
		for i,val in enumerate(self.data[val]):
			# ignore i == 0 it's just the filename
			if i > 1:
				output.append(val)
			elif i == 1:
				output.append(ctf)

		db_parms[name] = output

	def on_recall_params(self):
		if len(self.setlist.selectedItems()) == 0: return
			
		val = self.curset
		name = str(self.setlist.item(val).text())
		data = [name]
		name = get_file_tag(name)
		
		db_parms=db_open_dict("bdb:e2ctf.parms")
		if not db_parms.has_key(name):
			print "error, ctf parameters do not exist for:",name
#			
		
		data.extend(db_parms[name])
		ctf=EMAN2Ctf()
		ctf.from_string(data[1])
		data[1] = ctf
		
		self.data[val] = data
		self.newSet(self.curset)
	
#	def get_output_params(self):
	
	def on_output(self):
		from emsprworkflow import E2CTFOutputTaskGeneral
		self.form = E2CTFOutputTaskGeneral(self.app())
		self.form.run_form()
	
	def show_guis(self):
		if self.guiim != None:
			self.app().show_specific(self.guiim)
		if self.guiplot != None:
			self.app().show_specific(self.guiplot)
		
		self.app().show_specific(self.module())
	def closeEvent(self,event):
#		QtGui.QWidget.closeEvent(self,event)
#		self.app.app.closeAllWindows()
		if self.guiim != None:
			self.app().close_specific(self.guiim)
			self.guiim = None 
		if self.guiplot != None:
			self.app().close_specific(self.guiplot)
		event.accept()
		if self.module != None:
			self.app().close_specific(self.module())
			self.module().emit(QtCore.SIGNAL("module_closed")) # this signal is important when e2ctf is being used by a program running its own event loop

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
		val=self.curset
		ctf=self.data[val][1]
		ds=self.data[val][1].dsbg
		s=[ds*i for i in xrange(len(ctf.background))]
		if self.plotmode==1:
			self.guiplot.set_data((s,self.data[val][2]),"fg",True,True)
			self.guiplot.set_data((s,self.data[val][3]),"bg")
			self.guiplot.setAxisParms("s (1/A)","Intensity (a.u)")
		elif self.plotmode==0: 
			bgsub=[self.data[val][2][i]-self.data[val][3][i] for i in range(len(self.data[val][2]))]
			self.guiplot.set_data((s,bgsub),"fg-bg",True,True)
			
			fit=ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP)		# The fit curve
			sf = sfact(s, "ribosome","nono")
			fit=[sf[i]*fit[i]**2 for i in xrange(len(s))]		# squared * a generic structure factor

			# auto-amplitude for b-factor adjustment
			rto,nrto=0,0
			for i in xrange(int(.04/ds)+1,min(int(0.15/ds),len(s)-1)): 
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

			self.guiplot.set_data((s,fit),"fit")
			self.guiplot.setAxisParms("s (1/A)","Intensity (a.u)")
		elif self.plotmode==2:
			snr=ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_SNR)		# The snr curve
			self.guiplot.set_data((s,snr[:len(s)]),"snr",True)
			self.guiplot.setAxisParms("s (1/A)","SNR (intensity ratio)")
		elif self.plotmode==3:
			snr=ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_SNR)		# The snr curve
			self.guiplot.set_data((s,snr[:len(s)]),"snr",True)
			ssnr=ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_SNR_SMOOTH)		# The fit curve
			self.guiplot.set_data((s,ssnr[:len(s)]),"ssnr")
			self.guiplot.setAxisParms("s (1/A)","SNR (intensity ratio)")
		elif self.plotmode==4:
			snr=ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_SNR)		# The snr curve
			for i in range(1,len(snr)): snr[i]=snr[i]*i+snr[i-1]			# integrate SNR*s
#			for i in range(1,len(snr)): snr[i]/=snr[-1]				# normalize
			for i in range(1,len(snr)): snr[i]/=len(snr)			# this way the relative quality of images can be compared
			self.guiplot.set_data((s,snr[:len(s)]),"snr",True)
			self.guiplot.setAxisParms("s (1/A)","Integrated SNR")
		elif self.plotmode==5:
			inten=[fabs(i) for i in ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP)]		# The snr curve
			self.guiplot.set_data((s,inten[:len(s)]),"single",True)
			all=[0 for i in inten]
			for st in self.data:
				print st
				inten=[fabs(i) for i in st[1].compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP)]
				for i in range(len(all)): all[i]+=inten[i]
			self.guiplot.set_data((s,all[:len(s)]),"total")
						
			#bgsub=[self.data[val][2][i]-self.data[val][3][i] for i in range(len(self.data[val][2]))]
			#self.guiplot.set_data("fg-bg",(s,bgsub),True,True)
			
			#fit=[bgsub[i]/sfact(s[i]) for i in range(len(s))]		# squared * a generic structure factor

			#self.guiplot.set_data("fit",(s,fit))

	def newSet(self,val):
		"called when a new data set is selected from the list"
		self.curset=val

		self.sdefocus.setValue(self.data[val][1].defocus,True)
		self.sbfactor.setValue(self.data[val][1].bfactor,True)
#		self.sapix.setValue(self.data[val][1].apix)
		self.sampcont.setValue(self.data[val][1].ampcont,True)
		self.svoltage.setValue(self.data[val][1].voltage,True)
		self.scs.setValue(self.data[val][1].cs,True)
		
		self.guiim.set_data(self.data[val][4])
		self.update_plot()

	def newPlotMode(self,mode):
		self.plotmode=mode
		self.update_plot()

	def newCTF(self) :
		self.data[self.curset][1].defocus=self.sdefocus.value
		self.data[self.curset][1].bfactor=self.sbfactor.value
#		self.data[self.curset][1].apix=self.sapix.value
		self.data[self.curset][1].ampcont=self.sampcont.value
		self.data[self.curset][1].voltage=self.svoltage.value
		self.data[self.curset][1].cs=self.scs.value
		self.update_plot()

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
