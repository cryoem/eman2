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
import numpy as np

def get_usage():
	progname = os.path.basename(sys.argv[0])
	usage = progname + """ [options]
	Reconstructs 3D volumes using a set of 2D images. Euler angles are extracted from the 2D image headers and symmetry is imposed.
	This is an optimized version with fewer options than e2make3d, but can run in parallel using threads, and has a number of other changes.

	A simple example of usage is:

	e2make3d.py --input=classes.img --sym=c3 --output=recon.mrc --pad=128

	Because there are several default settings, this is more or less equivalent to:

	e2make3d.py --input=classes.img --sym=c3 --output=recon.mrc --pad=128 --keep=1 --recon=fourier --iter=3

	Because the padding is always done using zeroes it is best if your data (or preferably
	the edge pixels of your data) have mean zero. If you are unsure whether your data are
	appropriately normalized you can add the --preprocess flag

	e2make3d.py --input=classes.img --sym=c3 --output=recon.mrc --pad=128 --preprocess=normalize.edgemean

	You can add as many --preprocess arguments as you like, which are applied in
	the order in which they are specified, before padding occurs.

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

	parser.add_argument("--fillangle", type=float, dest="fillangle", help="An angular range used for both alt & az over which the projection should be averaged. Generally the angular step used when making projections.",default=0)
	parser.add_argument("--pad", metavar="x or x,y", default=None,type=str, help="Will zero-pad images to the specifed size (x,y) or (x,x) prior to reconstruction. If not specified no padding occurs.")
	parser.add_argument("--padvol", metavar="x or x,y,z", default=None,type=str, help="Defines the dimensions (x,y,z) or (x,x,x) of the reconstructed volume. If ommitted, implied value based on padded 2D images is used.")
	parser.add_argument("--outsize", metavar="x or x,y,z", default=None, type=str, help="Defines the dimensions (x,y,z) or (x,x,x) of the final volume written to disk, if ommitted, size will be based on unpadded input size")
	parser.add_argument("--savenorm", default=None, type=str, help="If set, will save the normalization volume showing Fourier space filling to the specified file")

	parser.add_argument("--keep", type=float, dest="keep", help="The fraction of slices to keep, based on quality scores (1.0 = use all slices). See keepsig.",default=1.0)
	parser.add_argument("--keepsig",action="store_true",default=False, dest="keepsig", help="If set, keep will be interpreted as a standard deviation coefficient instead of as a percentage.")
	parser.add_argument("--keepabs",action="store_true",default=False, dest="keepabs", help="If set, keep will refer to the absolute quality of the class-average, not a local quality relative to other similar sized classes.")
	parser.add_argument("--no_wt", action="store_true", dest="no_wt", default=False, help="This argument turns automatic weighting off causing all images to be weighted by 1. If this argument is not specified images inserted into the reconstructed volume are weighted by the number of particles that contributed to them (i.e. as in class averages), which is extracted from the image header (as the ptcl_repr attribute).")
	parser.add_argument("--mode", type=str, default="gauss_2", help="Fourier reconstruction 'mode' to use. The default should not normally be changed. default='gauss_2'")

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	parser.add_argument("--threads", default=4,type=int,help="Number of threads to run in parallel on a single computer. This is the only parallelism supported by e2make3dpar", guitype='intbox', row=24, col=2, rowspan=1, colspan=1, mode="refinement")
	parser.add_argument("--preprocess", metavar="processor_name(param1=value1:param2=value2)", type=str, action="append", help="preprocessor to be applied to the projections prior to 3D insertion. There can be more than one preprocessor and they are applied in the order in which they are specifed. Applied before padding occurs. See e2help.py processors for a complete list of available processors.")
	parser.add_argument("--setsf",type=str,help="Force the structure factor to match a 'known' curve prior to postprocessing (<filename>, auto or none). default=none",default="none")
	parser.add_argument("--postprocess", metavar="processor_name(param1=value1:param2=value2)", type=str, action="append", help="postprocessor to be applied to the 3D volume once the reconstruction is completed. There can be more than one postprocessor, and they are applied in the order in which they are specified. See e2help.py processors for a complete list of available processors.")
	parser.add_argument("--apix",metavar="A/pix",type=float,help="A/pix value for output, overrides automatic values",default=None)

	# Database Metadata storage
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()

#	options.fillangle=options.fillangle*pi/180.0
	# make sure that the user has atleast specified an input file
#	if len(args) <1:
#		parser.error("No input file specified")
#	options.input = args[0]

	if (not options.keep and not options.keepsig):
		print "Warning, neither the keep nor the keepsig argument was specified. Setting keep=1 (keeping 100% of inserted slices)"
		options.keep=1

	if options.input_model!=None : options.input_model=int(options.input_model)

	print "e2make3dpar.py"
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

	data=initialize_data(options.input,options.input_model,options.tlt,options.pad,options.no_wt,options.preprocess)

	# Filter out averages/images which aren't good enough
	if options.verbose: print "Filtering data, %d images ->"%len(data),

	# we have an absolute threshold
	if options.keepabs: thr=options.keep
	else :
		quals=np.array(sorted([i["quality"] for i in data]))
		if np.std(quals)==0 : thr=max(quals)+1.0
		else :
			if options.keepsig:
				thr=np.mean(quals)+options.keep*np.std(quals)
			else:
				try:
					thr=quals[int(floor(len(quals)*options.keep))]
				except: thr=max(quals)+1.0

	excluded=[max(i["fileslice"],i["filenum"]) for i in data if i["quality"]>thr]
	included=[max(i["fileslice"],i["filenum"]) for i in data if i["quality"]<=thr]
	data=[i for i in data if i["quality"]<=thr]
	ptclcount=sum([i["weight"] for i in data])

	if options.verbose: print "After filter, %d images"%len(data)

	# Get the reconstructor and initialize it correctly
	a = {"size":padvol,"sym":options.sym,"mode":options.mode,"verbose":options.verbose-1}
	if options.savenorm!=None : a["savenorm"]=options.savenorm
	recon=Reconstructors.get("fourier", a)

	#########################################################
	# The actual reconstruction

	threads=[threading.Thread(target=reconstruct,args=(data[i::options.threads],recon,options.preprocess,options.pad,
			options.fillangle,options.verbose-1)) for i in xrange(options.threads)]

	recon.setup()
	for i,t in enumerate(threads):
		if options.verbose>1: print "started thread ",i
		t.start()

	for t in threads: t.join()

	output = recon.finish(True)

	if options.verbose>0 : print "Finished Reconstruction"

	try:
		output.set_attr("ptcl_repr",ptclcount)
		if len(included)>0 : output.set_attr("threed_ptcl_idxs",included)
		if len(excluded)>0 : output.set_attr("threed_excl_ptcl_idxs",excluded)
		output.set_attr("threed_ptcl_src",data[0]["filename"])
	except:
		print "Warning, error setting reconstruction attributes"
#	progress += 10
#	E2progress(logid,float(progress)/total_progress)

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
			# need to be really careful about the corners
#			for i in range(sfcurve.get_size()):
#				if sfcurve.get_x(i)>1.0/(2.0*apix) : sfcurve.set_y(i,0.0)
			output.process_inplace("filter.setstrucfac",{"apix":apix,"strucfac":sfcurve})
#			output.process_inplace("filter.lowpass.tophat",{"apix":apix,"cutoff_abs":0.49})


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

def initialize_data(inputfile,inputmodel,tltfile,pad,no_weights,preprocess):
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
			#if not lowmem :
				#if nslice>1 : elem["data"]=get_processed_image(inputfile,0,i,preprocess,pad,nx,ny)
				#else : elem["data"]=get_processed_image(inputfile,i,-1,preprocess,pad,nx,ny)
			elem["quality"]=1.0
			data.append(elem)

		if len(data)!=n_input and len(data)!=nslice :
			print "Error : %d images and only %d lines in .tlt file"%(n_input,len(data))
			exit(1)
	else :
		tmp=EMData()
		for i in xrange(n_input):
			tmp.read_image(inputfile,i,True)
			#else : tmp=get_processed_image(inputfile,i,-1,preprocess,pad)

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

			try: elem["quality"]=float(tmp["class_qual"])
			except:
				try: elem["quality"]=1.0/(elem["weight"]+.00001)
				except: elem["quality"]=1.0
			elem["filename"]=inputfile
			elem["filenum"]=i
			elem["fileslice"]=-1
			#if not lowmem:
				#elem["data"]=tmp
				#tmp=EMData()

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

def reconstruct(data,recon,preprocess,pad,fillangle,verbose=0):
	"""Do an actual reconstruction using an already allocated reconstructor, and a data list as produced
	by initialize_data(). preprocess is a list of processor strings to be applied to the image data in the
	event that it hasn't already been read into the data array. start and startweight are optional parameters
	to seed the reconstruction with some non-zero values."""

	output=None		# deletes the results from the previous iteration if any

	if verbose>0:print "Inserting Slices (%d)"%len(data)

	astep=atan2(1.0,max(pad)/2.0)*180./pi
#	astep=atan2(1.0,max(pad)/2.0)/1.5*180./pi		# experimental smaller step size
	den=floor(fillangle/astep)
	if den>9 :
		den=9
		if verbose>0 : print "Note: Reducing oversampling in make3dpar for speed, this will make higher resolution 'smearing' less effective"
	if den<=1:
		fillangle=0
		if verbose: print "No filling"
	else:
		astep=fillangle/(den-1)-.00001
		if verbose: print "Filling %dx%d, %1.2f deg  %1.3f step"%(den,den,fillangle,astep)

	ptcl=0
	for i,elem in enumerate(data):
		# get the image to insert
		try:
			img=elem["data"]
			if img["sigma"]==0 : contine
			if not img.has_attr("reconstruct_preproc") :
				img=recon.preprocess_slice(img,elem["xform"])
				elem["data"]=img		# cache this for use in later iterations
		except:
#			print traceback.print_exc()
			if elem["fileslice"]>=0 : img=get_processed_image(elem["filename"],elem["filenum"],elem["fileslice"],preprocess,pad,elem["nx"],elem["ny"])
			else : img=get_processed_image(elem["filename"],elem["filenum"],-1,preprocess,pad)
			if img["sigma"]==0 : continue
			img=recon.preprocess_slice(img,elem["xform"])	# no caching here, with the lowmem option
#		img["n"]=i
#		if i==7 : display(img)

		rd=elem["xform"].get_rotation("eman")
		if verbose>0 : print " %d/%d\r"%(i,len(data)),
		sys.stdout.flush()
#		print "%d.\t%6.2f  %6.2f  %6.2f    %6.2g\t%6.4g\t%6.4g"%(i,rd["az"],rd["alt"],rd["phi"],elem["weight"],img["mean"],img["sigma"])
		ptcl+=elem["weight"]

		# Actual slice insertion into the volume
#		if i==len(data)-1 : display(img)
		if fillangle<=0:
			recon.insert_slice(img,elem["xform"],elem["weight"])
		else:
			xf=elem["xform"].get_rotation("eman")
			alt,az=xf["alt"],xf["az"]
			for dalt in np.arange(-fillangle/2.0,fillangle/2.0,astep):
				for daz in np.arange(-fillangle/2.0,fillangle/2.0,astep):
					weightmod=exp(-(dalt**2+daz**2)/(fillangle/4.0)**2)
					newxf=Transform({"type":"eman","alt":alt+dalt,"az":az+daz})
#					print i,elem["filenum"],newxf
					recon.insert_slice(img,newxf,elem["weight"]*weightmod)


	return

if __name__=="__main__":
	main()
