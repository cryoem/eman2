#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu),
# David Woolford 2007-2008 (woolford@bcm.edu)
# Copyright (c) 2000-2008 Baylor College of Medicine
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

# initial version of make3d

from EMAN2 import *
from EMAN2db import db_open_dict
from copy import deepcopy
from math import ceil
import os
import sys
import math
import random
import traceback
from numpy import array

def get_usage():
	progname = os.path.basename(sys.argv[0])
	usage = progname + """ [options]
	Reconstructs 3D volumes using a set of 2D images. Euler angles are extracted from the 2D image headers and symmetry is imposed. Several reconstruction methods are available (see e2help.py reconstructors) - the fourier method is the default and recommended reconstructor.

	WARNING - the icosahedral symmetry axes are aligned differently in eman1 and eman2 - that means
	using this script to reconstruct icosahedral class averages from eman1 will give very bad results.
	Contact the developers for the simple solution.

	A simple example of usage is:

	e2make3d.py --input=classes.img --sym=c3 --output=recon.mrc --pad=128

	Because there are several default settings, this is more or less equivalent to:

	e2make3d.py --input=classes.img --sym=c3 --output=recon.mrc --pad=128 --keep=1 --recon=fourier --iter=3

	Because the padding is always done using zeroes it is critical your data (or preferably
	the edge pixels of your data) have mean zero. This will be done automatically when the data is read, but without
	any amplitude scaling. If you wish to do additional preprocessing, you can use the --preprocess
	modifier to perform additional preprocessing, eg - 

	e2make3d.py --input=classes.img --sym=c3 --output=recon.mrc --pad=128 --preprocess=normalize.edgemean

	You can add as many --preprocess arguments as you like, which are applied in
	the order in which they are specified, before padding occurs.

	If you specify a value of the keep parameter that is not one, i.e. --keep=0.9, it
	allows for the automatic exclusion of those images that agree poorly with the rest
	of the data set.

	If constructing large volumes use the --lowmem option.
	"""
	return usage

def print_usage():

	usage = get_usage()
	print "usage " + usage;
	print "Please run '" + progname + " -h' for detailed options"

def main():
	parser = EMArgumentParser(usage=get_usage())

	parser.add_argument("--output", default="threed.hdf", help="Output reconstructed volume file name.")
	parser.add_argument("--input", default=None, help="The input projections. Project should usually have the xform.projection header attribute, which is used for slice insertion")
	parser.add_argument("--input_model", default=None, help="If the class-averages have the model_id parameter (produced by e2refinemulti.py), this will use only class-averages with the specified model_id for the reconstruction.")
	parser.add_argument("--tlt", help="An imod tlt file containing alignment angles. If specified slices will be inserted using these angles in the IMOD convention", type=str, default=None)
	parser.add_argument("--sym", dest="sym", default="c1", help="Set the symmetry; if no value is given then the model is assumed to have no symmetry.\nChoices are: i, c, d, tet, icos, or oct.")

	parser.add_argument("--pad", metavar="x or x,y", default=None,type=str, help="Will zero-pad images to the specifed size (x,y) or (x,x) prior to reconstruction. If not specified no padding occurs.")
	parser.add_argument("--padvol", metavar="x or x,y,z", default=None,type=str, help="Defines the dimensions (x,y,z) or (x,x,x) of the reconstructed volume. If ommitted, implied value based on padded 2D images is used.")
	parser.add_argument("--outsize", metavar="x or x,y,z", default=None, type=str, help="Defines the dimensions (x,y,z) or (x,x,x) of the final volume written to disk, if ommitted, size will be based on unpadded input size")

	parser.add_argument("--recon", dest="recon_type", default="fourier", help="Reconstructor to use see e2help.py reconstructors -v. Default is fourier:mode=gauss_2")
	parser.add_argument("--keep", type=float, dest="keep", help="The fraction of slices to keep, based on quality scores (1.0 = use all slices). See keepsig.",default=1.0)
	parser.add_argument("--keepsig",action="store_true",default=False, dest="keepsig", help="If set, keep will be interpreted as a standard deviation coefficient instead of as a percentage.")
	parser.add_argument("--keepabs",action="store_true",default=False, dest="keepabs", help="If set, keep will refer to the absolute quality of the class-average, not a local quality relative to other similar sized classes.")
	parser.add_argument("--no_wt", action="store_true", dest="no_wt", default=False, help="This argument turns automatic weighting off causing all images to be weighted by 1. If this argument is not specified images inserted into the reconstructed volume are weighted by the number of particles that contributed to them (i.e. as in class averages), which is extracted from the image header (as the ptcl_repr attribute).")
	parser.add_argument("--iter", type=int, dest="iter", default=2, help="Set the number of iterations (default is 2). Iterative reconstruction improves the overall normalization of the 2D images as they are inserted into the reconstructed volume, and allows for the exclusion of the poorer quality images.")
	parser.add_argument("--force", "-f",dest="force",default=False, action="store_true",help="deprecated")

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--nofilecheck",action="store_true",help="Turns file checking off in the check functionality - used by e2refine.py.",default=False)
	parser.add_argument("--check","-c",action="store_true",help="Performs a command line argument check only.",default=False)

	parser.add_argument("--lowmem",action="store_true",help="Causes images to be loaded as needed to reduce memory usage at the cost of time.",default=False)
	parser.add_argument("--preprocess", metavar="processor_name(param1=value1:param2=value2)", type=str, action="append", help="preprocessor to be applied to the projections prior to 3D insertion. There can be more than one preprocessor and they are applied in the order in which they are specifed. Applied before padding occurs. See e2help.py processors for a complete list of available processors.")
	parser.add_argument("--setsf",type=str,help="Force the structure factor to match a 'known' curve prior to postprocessing (<filename>, auto or none). default=none",default="none")
	parser.add_argument("--postprocess", metavar="processor_name(param1=value1:param2=value2)", type=str, action="append", help="postprocessor to be applied to the 3D volume once the reconstruction is completed. There can be more than one postprocessor, and they are applied in the order in which they are specified. See e2help.py processors for a complete list of available processors.")
	parser.add_argument("--apix",metavar="A/pix",type=float,help="A/pix value for output, overrides automatic values",default=None)

	parser.add_argument("--start", default=None,type=str, help="This is a starting model for FFT reconstruction")
	parser.add_argument("--startweight", default=1.0,type=float, help="This is the starting model weight")
	# Database Metadata storage
#	parser.add_argument("--dbls",type=str,default=None,help="data base list storage, used by the workflow. You can ignore this argument.")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

	if options.nofilecheck: options.check = True

	# make sure that the user has atleast specified an input file
#	if len(args) <1:
#		parser.error("No input file specified")
#	options.input = args[0]

	if (not options.keep and not options.keepsig):
		print "Warning, neither the keep nor the keepsig argument was specified. Setting keep=1 (keeping 100% of inserted slices)"
		options.keep=1

	# Checking is disabled for now. ie - check will always pass
	if options.check: exit(0)

	if options.input_model!=None : options.input_model=int(options.input_model)

	print "e2make3d.py"
	logger=E2init(sys.argv,options.ppid)

	# get basic image parameters
	tmp=EMData()
	tmp.read_image(options.input,0,True)
	if not options.no_wt :
		try:
			n=1
			while tmp["ptcl_repr"]==0 :
				tmp.read_image(options.input,n,True)
				n+=1
		except: pass

	nx=tmp["nx"]
	ny=tmp["ny"]
	nslice=tmp["nz"]
	if options.apix!=None : apix=options.apix
	else : apix=tmp["apix_x"]


	if options.verbose>0: print "Image dimensions %d x %d"%(nx,ny)

	# parse the padding options, to make sure we have a 2 or 3 tuple for each
	try :
		if options.pad==None : options.pad=(max(nx,ny),max(nx,ny))
		elif "," in options.pad :
			s=options.pad.split(",")
			options.pad=(int(s[0]),int(s[1]))
		else : options.pad=(int(options.pad),int(options.pad))
		if options.verbose>0 : print "Pad to %d x %d"%(options.pad[0],options.pad[1])
	except:
		print "Couldn't parse pad option :",options.pad
		exit(1)

	try :
		if options.padvol==None : padvol=(max(options.pad),max(options.pad),max(options.pad))
		elif "," in options.padvol :
			s=options.padvol.split(",")
			padvol=(int(s[0]),int(s[1]),int(s[2]))
		else : padvol=(int(options.padvol),int(options.padvol),int(options.padvol))
		if options.verbose>0 : print "Padded volume to reconstruct %d x %d x %d"%(padvol[0],padvol[1],padvol[2])
	except:
		print "Couldn't parse padvol option :",options.padvol
		exit(1)

	try :
		if options.outsize==None : outsize=(max(nx,ny),max(nx,ny),max(nx,ny))
		elif "," in options.outsize :
			s=options.outsize.split(",")
			outsize=(int(s[0]),int(s[1]),int(s[2]))
		else : outsize=(int(options.outsize),int(options.outsize),int(options.outsize))
		if options.verbose>0 : print "Final output volume %d x %d x %d"%(outsize[0],outsize[1],outsize[2])
	except:
		print "Couldn't parse outsize option :",options.outsize
		exit(1)

	data=initialize_data(options.input,options.input_model,options.tlt,options.pad,options.no_wt,options.lowmem,options.preprocess)

	# Get the reconstructor and initialize it correctly
	a = parsemodopt(options.recon_type)
	a[1]["size"]=padvol
	a[1]["sym"] = options.sym
	a[1]["verbose"] = options.verbose - 1
	recon=Reconstructors.get(a[0], a[1])

	start=None
	if options.start :
		start=EMData(options.start,0)
		if start["nx"]!=padvol[0] or start["ny"]!=padvol[1] or start["nz"]!=padvol[2] :
			start.clip_inplace(Region((start["nx"]-padvol[0])/2,(start["ny"]-padvol[1])/2,(start["nz"]-padvol[2])/2, padvol[0],padvol[1],padvol[2]))
		start.do_fft_inplace()

	#########################################################
	# The actual reconstruction
	output=reconstruct(data,recon,options.preprocess,options.pad,options.iter,options.keep,options.keepsig,start,options.startweight,options.verbose-1,options.keepabs)
	#
	########################################################3

	# clip to the requested final dimensions
	if output["nx"]!=outsize[0] or output["ny"]!=outsize[1] or output["nz"]!=outsize[2] :
		output.clip_inplace(Region((output["nx"]-outsize[0])/2,(output["ny"]-outsize[1])/2,(output["nz"]-outsize[2])/2, outsize[0],outsize[1],outsize[2]))

	if options.apix!=None : apix=options.apix
	output["apix_x"]=apix
	output["apix_y"]=apix
	output["apix_z"]=apix

	# Structure factor setting
	if options.setsf.lower() != "none" :
		if options.setsf.lower()!="auto" :
			try:
				sfcurve=XYData()
				sfcurve.read_file(options.setsf)
				for i in range(sfcurve.get_size()):
					v=sfcurve.get_y(i)
					if v<=0 :
						print "Warning values <=0 found in structure factor file. Please remove."
				sfcurve.update()
			except:
				print "ERROR: Specified structure factor ({}) not found.".format(options.setsf)
				sys.exit(1)
		else:
			try:
				sfcurve=XYData()
				sfcurve.read_file("strucfac.txt")
				for i in range(sfcurve.get_size()):
					v=sfcurve.get_y(i)
					if v<=0 :
						print "Warning values <=0 found in structure factor file. Please remove."
				sfcurve.update()
			except : sfcurve=None

		if sfcurve==None:
			print "ERROR : Structure factor read failed. Not applying structure factor"
		else:
			output.process_inplace("filter.setstrucfac",{"apix":apix,"strucfac":sfcurve})

	if options.postprocess != None:
		for p in options.postprocess:
			try:
				(processorname, param_dict) = parsemodopt(p)
				if not param_dict : param_dict={}
				output.process_inplace(str(processorname), param_dict)
			except:
				print "warning - application of the post processor",p," failed. Continuing anyway"

#	output.process_inplace("normalize.circlemean")

	# just remove the output file
	if file_exists(options.output):
		remove_file(options.output)

	# write the reconstruction to disk
	output.write_image(options.output,0)
	if options.verbose>0:
			print "Output File: "+options.output

	E2end(logger)

	print "Exiting"

def initialize_data(inputfile,inputmodel,tltfile,pad,no_weights,lowmem,preprocess):
	"""returned list will contain dictionaries containing metadata about each image, and optionally the image itself
	xform - a Transform object
	weight - a floating point weight
	filename - image source filename, not required if data is set
	filenum - image number within filename
	data - an EMData object, if not present or None, filename/num will be used
	nx,ny - size of a single image in a stack file"""

	n_input=EMUtil.get_image_count(inputfile)
	nx,ny,nslice= gimme_image_dimensions3D(inputfile)
	if n_input==1 and nslice>1 and tltfile==None : raise Exception,"Require tlt file to work with volumetric stacks"
	print n_input," input images"

	data=[]

	# The TLT file will override anything stored in the image itself, implies no_weights
	if tltfile:
		f=file(tltfile,'r')
		lines=f.readlines()
		for i,line in enumerate(lines):
			elem={"xform":Transform({"type":"eman","az":90,"alt":float(line),"phi":-90}),"weight":1.0}
			elem["filename"]=inputfile
			if nslice>1 :
				elem["filenum"]=0
				elem["fileslice"]=i
				elem["nx"]=nx
				elem["ny"]=ny
			else :
				elem["fileslice"]=-1
				elem["filenum"]=i
			if not lowmem :
				if nslice>1 : elem["data"]=get_processed_image(inputfile,0,i,preprocess,pad,nx,ny)
				else : elem["data"]=get_processed_image(inputfile,i,-1,preprocess,pad,nx,ny)
			data.append(elem)

		if len(data)!=n_input and len(data)!=nslice :
			print "Error : %d images and only %d lines in .tlt file"%(n_input,len(data))
			exit(1)
	else :
		tmp=EMData()
		for i in xrange(n_input):
			if lowmem: tmp.read_image(inputfile,i,True)
			else : tmp=get_processed_image(inputfile,i,-1,preprocess,pad)

			# these rely only on the header
			try: elem={"xform":tmp["xform.projection"]}
			except : continue
				#raise Exception,"Image %d doesn't have orientation information in its header"%i

			# skip any particles targeted at a different model
			if inputmodel != None and tmp["model_id"]!=inputmodel : continue

			if no_weights: elem["weight"]=1.0
			else :
				try: elem["weight"]=float(tmp["ptcl_repr"])
				except: elem["weight"]=1.0
				# This is bad if you have actual empty classes...
				#if elem["weight"]<=0 :
					#print "Warning, weight %1.2f on particle %d. Setting to 1.0"%(elem["weight"],i)
					#elem["weight"]=1.0

			elem["filename"]=inputfile
			elem["filenum"]=i
			elem["fileslice"]=-1
			if not lowmem:
				elem["data"]=tmp
				tmp=EMData()

			data.append(elem)

		print "Using %d images"%len(data)
	return data

def get_processed_image(filename,nim,nsl,preprocess,pad,nx=0,ny=0):
	"""This function will read an image from disk and preprocess it. nim is the number of the image in
	the file. nsl is the number of a slice within the image if it is 3-D, pass -1 to avoid Region reading.
	nx and ny are the x and y dimensions of a single image on disk, required only for slice reading.
	preprocess takes a list of command-line processor strings. pad is a 2-tuple with
	the dimensions the image should be zero-padded to."""

	if nsl>=0 : ret=EMData(filename,nim,False,Region(0,0,nsl,nx,ny,1))
	else : ret=EMData(filename,nim)

	ret-=ret.get_edge_mean()			# It is critical that the edges be zero, even if we aren't changing the scaling

	# Apply any preprocessing (like normalize.edgemean)
	if isinstance(preprocess,str) : preprocess=[preprocess]
	if preprocess != None and len(preprocess)>0:
		for p in preprocess:
			(processorname, param_dict) = parsemodopt(p)
			if not param_dict : param_dict={}
			ret.process_inplace(str(processorname), param_dict)

	if pad[0]!=ret.get_xsize() or pad[1]!=ret.get_ysize() :
		ret.clip_inplace(Region((ret.get_xsize()-pad[0])/2,(ret.get_ysize()-pad[1])/2, pad[0], pad[1]))

	return ret

def reconstruct(data,recon,preprocess,pad,niter=2,keep=1.0,keepsig=False,start=None,startweight=0,verbose=0,keepabs=False):
	"""Do an actual reconstruction using an already allocated reconstructor, and a data list as produced
	by initialize_data(). preprocess is a list of processor strings to be applied to the image data in the
	event that it hasn't already been read into the data array. start and startweight are optional parameters
	to seed the reconstruction with some non-zero values."""

	for it in xrange(niter) :
		excluded=[]
		included=[]
		output=None		# deletes the results from the previous iteration if any

		if verbose>0: print "Initializing the reconstructor ..."

		if start : recon.setup(start,startweight)
		else : recon.setup()

		if verbose>0:print "Inserting Slices (%d)"%len(data)

		ptcl=0
		for i,elem in enumerate(data):
			if not elem.has_key("norm") : elem["norm"]=1.0

			# If the image is below the quality cutoff, skip it
			try:
				if (keepabs and elem["reconstruct_absqual"]<qcutoff) or (not keepabs and elem["reconstruct_qual"]<qcutoff) or elem["weight"]==0 :
					if verbose>0 : print i," *  (%1.3g)"%(elem["reconstruct_qual"])
					if it==niter-1 :
						if elem["fileslice"]<0 : excluded.append(elem["filenum"])
						else : excluded.append(elem["fileslice"])
					continue
			except: pass

			if elem["fileslice"]<0 : included.append(elem["filenum"])
			else : included.append(elem["fileslice"])

			# get the image to insert
			try:
				img=elem["data"]
				if img["sigma"]==0 : contine
				if not img.has_attr("reconstruct_preproc") :
					img=recon.preprocess_slice(img,elem["xform"])
					elem["data"]=img		# cache this for use in later iterations
				elem["data"].mult(elem["norm"])
				elem["norm"]=1.0
			except:
	#			print traceback.print_exc()
				if elem["fileslice"]>=0 : img=get_processed_image(elem["filename"],elem["filenum"],elem["fileslice"],preprocess,pad,elem["nx"],elem["ny"])
				else : img=get_processed_image(elem["filename"],elem["filenum"],-1,preprocess,pad)
				if img["sigma"]==0 : continue
				img=recon.preprocess_slice(img,elem["xform"])	# no caching here, with the lowmem option
				img.mult(elem["norm"])
	#		img["n"]=i
	#		if i==7 : display(img)

			rd=elem["xform"].get_rotation("eman")
			if verbose>0 : print " %d/%d\r"%(i,len(data)),
			sys.stdout.flush()
	#		print "%d.\t%6.2f  %6.2f  %6.2f    %6.2g\t%6.4g\t%6.4g"%(i,rd["az"],rd["alt"],rd["phi"],elem["weight"],img["mean"],img["sigma"])
			ptcl+=elem["weight"]

			# Actual slice insertion into the volume
	#		if i==len(data)-1 : display(img)
			recon.insert_slice(img,elem["xform"],elem["weight"])

		if it!=niter-1:
			if verbose>0: print "\t   az      alt    phi  \t tweight      norm     absqual    weight"
			for i,elem in enumerate(data):

				try:
					if (keepabs and elem["reconstruct_absqual"]<qcutoff) or (not keepabs and elem["reconstruct_qual"]<qcutoff) or elem["weight"]==0:
						if verbose>0 : print i," *"
#						continue
				except: pass

				# get the image to insert
				try:
					img=elem["data"]
				except:
					if elem["fileslice"]>=0 : img=get_processed_image(elem["filename"],elem["filenum"],elem["fileslice"],preprocess,pad,elem["nx"],elem["ny"])
					else : img=get_processed_image(elem["filename"],elem["filenum"],-1,preprocess,pad)
					img=recon.preprocess_slice(img,elem["xform"])	# no caching here, with the lowmem option
					img.mult(elem["norm"])

				rd=elem["xform"].get_rotation("eman")
				if verbose>0 : print "%d.\t%6.2f  %6.2f  %6.2f\t"%(i,rd["az"],rd["alt"],rd["phi"]),

				if img["sigma"]==0 :
					elem["reconstruct_weight"]=0
					elem["reconstruct_norm"]=0
					elem["reconstruct_absqual"]=0
					if verbose>0 : print ""
					continue

				dosub=True
				try:
					if (keepabs and elem["reconstruct_absqual"]<qcutoff) or (not keepabs and elem["reconstruct_qual"]<qcutoff) or elem["weight"]==0: dosub=False
				except: pass

				recon.determine_slice_agreement(img,elem["xform"],elem["weight"],dosub)

				# These are the parameters returned by determine_slice_agreement. Since the images may be reloaded, we cache them in the data dictionary
				for i in ("reconstruct_weight","reconstruct_norm","reconstruct_absqual"):
					elem[i]=img[i]
					if verbose>0 : print "%8.3g  "%elem[i],
				elem["norm"]*=img["reconstruct_norm"]

				if verbose>0 : print "%d"%int(elem["weight"])

			# Convert absolute qualities to relative qualities by local averaging vs classes with similar numbers of particles
			squal=sorted([[data[i]["weight"],data[i]["reconstruct_absqual"],0,i,data[i]] for i in xrange(len(data))])

			# compute a running average of qualities (sorted by weight), then measure each value vs. the local average
			qlist=[]
			for i in range(len(squal)):
				# skip averages with weight=0
				if squal[i][0]==0 :
					squal[i][4]["reconstruct_qual"]=-1
					continue
				sub=[squal[j][1] for j in xrange(max(0,i-10),min(i+10,len(squal)))]
				squal[i][2]=sum(sub)/len(sub)
				try:
					squal[i][4]["reconstruct_qual"]=squal[i][1]/squal[i][2]
					qlist.append(squal[i][1]/squal[i][2])
				except :
					traceback.print_exc()
					print "##############  ",sub, squal[i], i,len(squal)
					squal[i][4]["reconstruct_qual"]=-1
					qlist.append(-1)
#			plot(qlist)

			# we also get some statistics on the absolute qualities in case the user doesn't want local averaging
			sq=array(sorted([data[i]["reconstruct_absqual"] for i in xrange(len(data))]))
			aqmean=sq.mean()
			aqsigma=sq.std()

			# set exclusion thresholds
			if keepabs:
				if keepsig:
					qcutoff=aqmean-aqsigma*keepabs
					if verbose>0: print "Absolute Quality: mean=%1.3f sigma=%1.3f  ->  cutoff = %1.3f (%d ptcl)"%(aqmean,aqsigma,qcutoff,ptcl)
				else:
					qcutoff=sq[-int(keep*len(qlist))]
					if verbose>0: print "Absolute Quality: min=%1.3f max=%1.3f  ->  cutoff = %1.3f (%d ptcl)"%(sq[0],sq[-1],qcutoff,ptcl)
			else:
				if keepsig:
					qmean=sum(qlist)/len(qlist)
					qsigma=sqrt(sum([i*i for i in qlist])/len(qlist)-qmean**2)
					qcutoff=qmean-qsigma*keep
					if verbose>0: print "Quality: mean=%1.3f sigma=%1.3f  ->  cutoff = %1.3f (%d ptcl)"%(qmean,qsigma,qcutoff,ptcl)
				else:
					qlist.sort()
					qcutoff=qlist[-int(keep*len(qlist))]
					if verbose>0: print "Quality: min=%1.3f max=%1.3f  ->  cutoff = %1.3f (%d ptcl)"%(qlist[0],qlist[-1],qcutoff,ptcl)

		output = recon.finish(True)

	try:
		output.set_attr("ptcl_repr",ptcl)
		if len(included)>0 : output.set_attr("threed_ptcl_idxs",included)
		if len(excluded)>0 : output.set_attr("threed_excl_ptcl_idxs",excluded)
		output.set_attr("threed_ptcl_src",data[0]["filename"])
	except:
		print "Warning, error setting reconstruction attributes"
#	progress += 10
#	E2progress(logid,float(progress)/total_progress)
	if verbose>0 : print "Finished Reconstruction"

	return output

if __name__=="__main__":
	main()
