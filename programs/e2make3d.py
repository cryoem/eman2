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
from optparse import OptionParser
from copy import deepcopy
from math import ceil
import os
import sys
import math
import random
HEADER_ONLY=True
HEADER_AND_DATA=False


CONTRIBUTING_PARTICLE_LIMIT = 1000000

def get_usage():
	progname = os.path.basename(sys.argv[0])
	usage = progname + """ <inputfile> [options]
	Reconstructs 3D volumes using a set of 2D images. Euler angles are extracted from the 2D image headers and symmetry is imposed. Several reconstruction methods are available (see e2help.py reconstructors) - the fourier method is the default and recommended reconstructor.
	
	WARNING - the icosahedral symmetry axes are aligned dirrently in eman1 and eman2 - that means
	using this script to reconstruct icosahedral class averages from eman1 will give very bad results.
	Contact the developers for the simple solution. This problem will probably disappear - eman2
	will probably adopt the same symmetric axes as eman1, there needs to be some discussion before that
	occurs. For the time being assume there is a problem.
	
	A simple example of usage is:
	
	e2make3d.py classes.img --sym=c3 --out=recon.mrc --pad=128
	
	Because there are several default settings, this is more or less equivalent to:
	
	e2make3d.py classes.img --sym=c3 --out=recon.mrc --pad=128 --keep=1 --recon=fourier --iter=3
	
	Because the padding is always done using zeroes it is best if your data (or preferably
	the edge pixels of your data) have mean zero. If you are unsure whether your data are 
	appropriately normalized you can add the --preprocess flag
	
	e2make3d.py classes.img --sym=c3 --out=recon.mrc --pad=128 --preprocess=normalize.edgemean
	
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
	parser=OptionParser(usage=get_usage())
	parser.add_option("--output", default="threed.hdf", help="Output reconstructed volume file name.")
	parser.add_option("--input", default=None, help="The input projections. Project should usually have the xform.projection header attribute, which is used for slice insertion")

	parser.add_option("--sym", dest="sym", default="c1", help="Set the symmetry; if no value is given then the model is assumed to have no symmetry.\nChoices are: i, c, d, tet, icos, or oct.")
	parser.add_option("--recon", dest="recon_type", default="fourier", help="Reconstructor to use see e2help.py reconstructors -v. Default is fourier using mode 2.")
	parser.add_option("--verbose", "-v",dest="verbose",default=False, action="store_true",help="Toggle verbose mode - prints extra infromation to the command line while executing.")
	parser.add_option("--keep", type=float, dest="keep", help="The percentage of slices to keep, based on quality scores.")
	parser.add_option("--keepsig",action="store_true",default=False, dest="keepsig", help="The standard deviation alternative to the --keep argument.")
	
	parser.add_option("--no_wt", action="store_true", dest="no_wt", default=False, help="This argument turns automatic weighting off causing all images to be weighted by 1. If this argument is not specified images inserted into the reconstructed volume are weighted by the number of particles that contributed to them (i.e. as in class averages), which is extracted from the image header (as the ptcl_repr attribute).")
	parser.add_option("--iter", type=int, dest="iter", default=3, help="Set the number of iterations (default is 3). Iterative reconstruction improves the overall normalization of the 2D images as they are inserted into the reconstructed volume, and allows for the exclusion of the poorer quality images.")
	parser.add_option("--force", "-f",dest="force",default=False, action="store_true",help="Force overwrite the output file if it exists.")
	parser.add_option("--nofilecheck",action="store_true",help="Turns file checking off in the check functionality - used by e2refine.py.",default=False)
	parser.add_option("--check","-c",action="store_true",help="Performs a command line argument check only.",default=False)
	parser.add_option("--lowmem",action="store_true",help="Causes images to be loaded as the are inserted into the 3D volume, as opposed to having them all read and stored in memory for the duration of the program.",default=False)
	parser.add_option("--tlt", help="An imod tlt containing alignment angles. If specified slices will be inserted using these angles", type="string", default=None)
	
	parser.add_option("--preprocess", metavar="processor_name(param1=value1:param2=value2)", type="string", action="append", help="preprocessor to be applied to the projections prior to 3D insertion. There can be more than one preprocessor and they are applied in the order in which they are specifed. Applied before padding occurs. See e2help.py processors for a complete list of available processors.")
	parser.add_option("--postprocess", metavar="processor_name(param1=value1:param2=value2)", type="string", action="append", help="postprocessor to be applied to the 3D volume once the reconstruction is completed. There can be more than one postprocessor, and they are applied in the order in which they are specified. See e2help.py processors for a complete list of available processors.")
	
	parser.add_option("--pad", metavar="m or m,n", default=None,type="string", help="This can be a single value or two values. If a single value is specified (m) the input images are padded with zeroes uniformly in both directions so that the dimensions are mxm. If two values are specified (m,n) the images are padded with zeroes such that the dimension are mxn (in x and y, respectively). Padding occurs after preprocessing.")
	
	parser.add_option("--start", default=None,type="string", help="This is a starting model for FFT reconstruction")
	parser.add_option("--smw", default=1.0,type="float", help="This is the starting model weight")
	# Database Metadata storage
	parser.add_option("--dbls",type="string",default=None,help="data base list storage, used by the workflow. You can ignore this argument.")
	
	(options, args) = parser.parse_args()

	if options.nofilecheck: options.check = True
	
	# make sure that the user has atleast specified an input file
#	if len(args) <1:
#		parser.error("No input file specified")
#	options.input = args[0]
	
	if (not options.keep and not options.keepsig):
		print "Warning, neither the keep nor the keepsig argument was specified. Setting keep=1 (keeping 100% of inserted slices)"
		options.keep=1
		
	if (options.check): options.verbose = True # turn verbose on if the user is only checking...
	
	if not options.nofilecheck: options.nx,options.ny,options.nz = gimme_image_dimensions3D(options.input)
	if options.tlt:# the check function makes everyone I assume about the tlt aspect true
		options.no_wt = True
		angles = None
		f=file(options.tlt,'r')
		lines=f.readlines()
		angles=[]
		for line in lines:
			angles.append(float(line))
			
		options.angles = angles
		
	# Check for potential errors be examining the contents of the options
	error = check(options,True)
		 
	
	# If there were errors then exit
	if (options.verbose):
		if (error):
			print "e2make3d.py command line arguments test.... FAILED"
		else:
			print "e2make3d.py command line arguments test.... PASSED"
	if error : exit(1)
	# e2refine is sensitive to the return code, so it's important that 0 is returned here
	if options.check: exit(0)
	

	print "e2make3d.py"
	logger=E2init(sys.argv)
	
	# just remove the output file - if the user didn't specify force then the error should have been found in the check function
	if file_exists(options.output):
		if options.force:
			remove_file(options.output)
		else: 
			print "foo, this error should have already been found in check"
			exit(1)
	
	# store image dimensions early on - it saves lots of time at other stages in the code
	image_0 = get_processed_image(options,0,True)
	(xsize, ysize )= image_0.get_xsize(),image_0.get_ysize()
	options.xsize = xsize
	options.ysize = ysize
	options.zsize = xsize
	if ( ysize < xsize ): options.zsize = ysize
	if ysize == 1:
		options.ndim = 1
	else:
		options.ndim = 2
	
	# if we don't have to be careful about memory, store all the images in memory now
	if (options.lowmem == False):
		images = []
		if options.tlt: total_images = options.nz
		else: total_images=EMUtil.get_image_count(options.input)
		force_read = True
		for i in range(0,total_images):
			images.append(get_processed_image(options,i,force_read))
		
		options.images = images
		
	

	# Some reconstructors (only one atm) have special sampling properties, this complicates
	# post depadding of the reconstruction (if padding was used). 
	options.suppress_x_clip = False
	options.suppress_y_clip = False
	options.suppress_z_clip = False
	
	# Now do the reconstruction
	(recon_type, parms) = (parsemodopt(options.recon_type))
	if recon_type=="fourier":
		output=fourier_reconstruction(options,logger)
		# If the user is specifying reconstructor sampling, then suppress
		# depadding after reconstruction in the particular direction
		if "xsample" in parms: options.suppress_x_clip = True
		if "ysample" in parms: options.suppress_y_clip = True
		if "zsample" in parms: options.suppress_z_clip = True
			
	elif recon_type=="baldwinwoolford":
		output=bw_reconstruction(options)
	elif recon_type=="back_projection":
		output=back_projection_reconstruction(options)
	elif recon_type=="fourier2D":
		output=fourier2D_reconstruction(options,logger)
	else:
		# this point should never be reached
		sys.stderr.write("%s reconstuctor is not supported" %options.recon_type)
		exit(1)

	# post depad the reconstruction if padding was used.
	pad_dims = parse_pad(options.pad)
	if pad_dims != None:
		if options.ndim == 2:
			if len(pad_dims) == 1:
				xpad = pad_dims[0]
				ypad = pad_dims[0]
			elif len(pad_dims) == 2:
				xpad = pad_dims[0]
				ypad = pad_dims[1]
			
			if options.suppress_x_clip: 
				startx = 0
				lengthx = output.get_xsize()
			else:
				startx = (xpad-options.nx)/2
				lengthx = options.nx
				
			if options.suppress_y_clip: 
				starty = 0
				lengthy = output.get_ysize()
			else:
				starty = (ypad-options.ny)/2
				lengthy = options.ny
				
			if options.suppress_z_clip: 
				startz = 0
				lengthz = output.get_zsize()
			else:
				startz = (xpad-options.nx)/2
				lengthz = options.nx
			output.clip_inplace(Region(startx,starty,startz,lengthx,lengthy,lengthz))
		elif options.ndim == 1:
			pad = pad_dims[0]
			output.clip_inplace(Region((pad-options.nx)/2,(pad-options.nx)/2,options.nx,options.nx))

	# apply any post processing operations
	if db_check_dict("bdb:project"):
		# this is a temporary workaround to get things working in the workflow
		db = db_open_dict("bdb:project")
		apix = db.get("global.apix",dfl=1)
		output.set_attr("apix_x",apix)
		output.set_attr("apix_y",apix)
		output.set_attr("apix_z",apix)
	
	if options.postprocess == None:
		options.postprocess = ["normalize.circlemean"]
	if options.postprocess != None:
		for p in options.postprocess:
			try:
				(processorname, param_dict) = parsemodopt(p)
				if not param_dict : param_dict={}
				output.process_inplace(str(processorname), param_dict)
			except:
				print "warning - application of the post processor",p," failed. Continuing anyway"
					
	# write the reconstruction to disk
	output.write_image(options.output,0)
	if options.verbose:
			print "Output File: "+options.output

	if options.dbls:
		pdb = db_open_dict("bdb:project")
		tmp_d = pdb.get(options.dbls, dfl={})
		if isinstance(tmp_d,list): # this was added June 2009 - it's a back compatibility measure. We used to store these things as lists, but now it's dicts
			d = {}
			for name in tmp_d:
				s = {}
				s["Original Data"] = name
				d[name] = s
			
			tmp_d = d
		s = {}
		s["Original Data"] = options.output
		tmp_d[options.output] = s
		pdb[options.dbls] = tmp_d

	E2end(logger)
	
	exit(0)
	print "Exiting"

# this function will either read a particular image from disk, or return the image if it's stored in the options
def get_processed_image(options,i, force_read_from_disk=False):
	if (options.lowmem or force_read_from_disk):
		d = EMData()
		if options.tlt:
			roi=Region(0,0,i,options.nx,options.ny,1)
			d.read_image(options.input,0, HEADER_AND_DATA, roi)
			d.set_attr("xform.projection",Transform({"type":"eman","az":90,"alt":options.angles[i],"phi":90}))
		else:
			d.read_image(options.input, i)
		
		if options.preprocess != None:
			for p in options.preprocess:
				(processorname, param_dict) = parsemodopt(p)
				if not param_dict : param_dict={}
				d.process_inplace(str(processorname), param_dict)
				
		pad_dims = parse_pad(options.pad)
		if pad_dims != None:
			if len(pad_dims) == 1:
				xpad = pad_dims[0]
				ypad = pad_dims[0]
			elif len(pad_dims) == 2:
				xpad = pad_dims[0]
				ypad = pad_dims[1]
			
			ndim = d.get_ndim()
			if ndim == 2:
				d.clip_inplace(Region((d.get_xsize()-xpad)/2,(d.get_ysize()-ypad)/2, xpad, ypad))
			elif ndim == 1:
				d.clip_inplace(Region((d.get_xsize()-xpad)/2, xpad))
			#else should never happen
		return d
	else:
		return options.images[i]

def gimme_image_dimensions2D_consider_pad(options):
	
	if options.ndim == 2:
		pad_dims = parse_pad(options.pad)
		if pad_dims != None:
			if len(pad_dims) == 1:
				xsize = pad_dims[0]
				ysize = pad_dims[0]
			elif len(pad_dims) == 2:
				xsize = pad_dims[0]
				ysize = pad_dims[1]
				
			return (xsize, ysize)
			
		else:
			return (options.xsize,options.ysize)
	elif options.ndim == 1:
		pad_dims = parse_pad(options.pad)
		if pad_dims != None:
			return (pad_dims[0],1)
		else:
			return (options.xsize,1)
	#else should never happen

def parse_pad(padstring):
	if padstring == None: return None
	else :
		tmp = padstring.split(',')
		ret = []
		for i in tmp:
			ret.append(int(i))
			
		return ret


def bw_reconstruction(options):
	# Get the reconstructor and initialize it correctly
	a = parsemodopt(options.recon_type)
	recon=Reconstructors.get(a[0], a[1])
	params = recon.get_params()
	(xsize, ysize ) = gimme_image_dimensions2D_consider_pad( options )
	params["x_in"] = xsize;
	params["y_in"] = ysize;
	params["sym"] = options.sym
	recon.insert_params(params)
	recon.setup()
	
	if options.tlt: total_images = options.nz
	else: total_images=EMUtil.get_image_count(options.input)
	
	
	if (options.verbose):
		print "Accruing densities"

	for i in xrange(0,total_images):
		d = EMData()
		d.read_image(options.input, i, True)
		t = d.get_attr("xform.projection")
		recon.insert_slice_weights(t)
	
	if (options.verbose):
		print "Inserting slices"

	for i in xrange(0,total_images):
		
		d = get_processed_image(options,i)
		
		weight = float (d.get_attr("ptcl_repr"))
		d.mult(weight)
	
		t = d.get_attr("xform.projection")
		r = t.get_params("eman")
		failure = recon.insert_slice(d,t)

	return recon.finish()

def back_projection_reconstruction(options):
	
	if options.verbose:
		print "Initializing the reconstructor"

	a = parsemodopt(options.recon_type)

	recon=Reconstructors.get(a[0], a[1])

	params = recon.get_params()
	(xsize, ysize) = gimme_image_dimensions2D_consider_pad( options );
	if ( xsize != ysize ):
		print "Error, back space projection currently only works for images with uniform dimensions"
		exit(1)
		
	print xsize
	params["size"] = xsize;
	params["sym"] = options.sym
	recon.insert_params(params)
	recon.setup()

	if options.tlt: total_images = options.nz
	else: total_images=EMUtil.get_image_count(options.input)
	
	for i in xrange(total_images):
		d = get_processed_image(options,i)
		
		try:
			num_img=d.get_attr("ptcl_repr")
		except:
			num_img = 1 
		if ( num_img<=0 and options.no_wt == False):
			continue
		else:
			num_img = 1
		
		weight = 1
		if ( options.no_wt == False ):
			weight = float (num_img)
		
		param = {}
		param["weight"] = weight
		recon.insert_params(param)
		
		t = d.get_attr("xform.projection")
		r = t.get_params("eman")
		recon.insert_slice(d, t)

		if options.verbose:
			print "%2d/%d  %3d\t%5.1f  %5.1f  %5.1f\t\t%6.2g %6.2g" %(
					(i+1,total_images,num_img,
					r["alt"], r["az"],r["phi"],
					d.get_attr("maximum"),d.get_attr("minimum")))
	
	return recon.finish()
		
def fourier2D_reconstruction(options):
	a = parsemodopt(options.recon_type)
	recon=Reconstructors.get(a[0], a[1])
	(xsize, ysize) = gimme_image_dimensions2D_consider_pad( options )
	if (ysize != 1):
		print "Error, dimensions are greater than 1"
		exit(1)
		
	params = recon.get_params()
	params["nx"] = xsize;
	params["sym"] = options.sym
	
	recon.insert_params(params)
	recon.setup()
	
	if options.tlt: raise NotImplementedException("The 2D reconstruction does not work the the tlt option")
	
	total_images=EMUtil.get_image_count(options.input)
	
	for i in xrange(0,total_images):
		image = get_processed_image(options,i)
	
		transform = Transform({"type":"2d","alpha":image.get_attr("euler_alt")})
		r = transform.get_params("2d")
		failure = recon.insert_slice(image,transform)
			
		if (options.verbose):
			sys.stdout.write( "%2d/%d  %3d\t%5.1f  %5.1f  %5.1f\t\t%6.2f %6.2f" %
							(i+1,total_images, image.get_attr("IMAGIC.imgnum"),
							r["alpha"],image.get_attr("maximum"),image.get_attr("minimum")))
				
			if ( failure ):
				sys.stdout.write( " X" )
			
			sys.stdout.write("\n")


	if (options.verbose):
		print "Inverting 3D Fourier volume to generate the real space reconstruction"
	output = recon.finish()

	return output
	
def fourier_reconstruction(options,logid):
	if (options.verbose):
		print "Initializing the reconstructor ..."
	
	# Get the reconstructor and initialize it correctly
	a = parsemodopt(options.recon_type)
	recon=Reconstructors.get(a[0], a[1])
	params = recon.get_params()
	(xsize, ysize) = gimme_image_dimensions2D_consider_pad( options )
	params["x_in"] = xsize;
	params["y_in"] = ysize;
	params["sym"] = options.sym
	params["quiet"] = not options.verbose
	
	if options.start:
		start = EMData(options.start)
		x = start.get_xsize()
		y = start.get_ysize()
		z = start.get_zsize()
		pad_dims = parse_pad(options.pad)
		if pad_dims != None:
			if len(pad_dims) == 1:
				xpad = pad_dims[0]
				ypad = pad_dims[0]
				zpad = pad_dims[0]
			else:
				# non uniform needs some thinking
				raise
				#xpad = pad_dims[0]
				#ypad = pad_dims[1]
				#zpad = pad_dims[0]
			
			start.clip_inplace(Region((x-xpad)/2,(y-ypad)/2,(z-zpad)/2, xpad, ypad,zpad))

		params["start_model"] = start
		params["start_model_weight"] = options.smw
		
	
	recon.insert_params(params)
	recon.setup()
	
	if options.tlt: total_images = options.nz
	else: 
		try: total_images=EMUtil.get_image_count(options.input)
		except: total_images=len(options.images)

	if (options.verbose):
		print "Inserting Slices"
	
	total_progress = (options.iter+1)*total_images + (options.iter)*total_images+10 
	progress = 0
	for j in xrange(0,options.iter+1):
		removed = 0;
		
		if ( j > 0 ):
			tot = 0
			if options.verbose: print "Determining slice agreement"
			ptcl_repr =[]
			for i in range(total_images):
				image = get_processed_image(options,i)
				
				if image.has_attr("ptcl_repr"):
					num_img=image.get_attr("ptcl_repr")
				else:
					num_img = 0
				ptcl_repr.append(num_img) 
				if image.has_attr("ptcl_repr"):
					if (num_img<=0 and options.no_wt == False):
						continue
					else:
						num_img = 1
				
				if ( options.no_wt == False ):
					weight = float (num_img)
					weight_params = {"weight":weight}
					recon.insert_params(weight_params)
				
				t = image.get_attr("xform.projection")
				recon.determine_slice_agreement(image,t,num_img)
				progress += 1
				E2progress(logid,float(progress)/total_progress)
				tot += 1
				if options.verbose:
					sys.stdout.write(".")
					sys.stdout.flush()
	
			if ( options.keep != 1.0 or options.keepsig == True  ):
				idx = 0
				fsc_scores = []
				hard = 0.0
				for i in range(total_images):
					if (ptcl_repr[i]<=0 and options.no_wt == False): continue
					else:
						fsc_scores.append(-recon.get_score(idx))
						idx += 1
					
					
				if ( options.keepsig ):
					a = Util.get_stats_cstyle(fsc_scores)
					mean = a["mean"]
					std_dev = a["std_dev"]
					hard  = -(mean + options.keepsig*std_dev)
				else:
					b = deepcopy(fsc_scores)
					b.sort()
					# The ceil reflects a conservative policy. If the user specified keep=0.93
					# and there were 10 particles, then they would all be kept. If floor were
					# used instead of ceil, the last particle would be thrown away (in the
					# class average)
					idx = int(ceil(options.keep*len(b))-1)
					hard  = -b[idx]
				param = {}
				param["hard"] = hard
				recon.insert_params(param)
				if options.verbose: print "using new hard parameter %f" %hard
			else:
				param = {}
				param["hard"] = 0.0
				recon.insert_params(param)
				if options.verbose: print "using hard parameter 0.0"
			if options.verbose: print " Done"
			
		idx = 0	
		contributing_ptcls = 0 # calculate how many particles went into the final reconstruction
		for i in range(total_images):
			#print i
			image = get_processed_image(options,i)
			
			if image.has_attr("ptcl_repr"):
				if (image.get_attr("ptcl_repr")<=0 and options.no_wt == False):
					continue
				
			if ( options.no_wt == False ):
				weight = float (image.get_attr("ptcl_repr"))
				weight_params = {"weight":weight}
				recon.insert_params(weight_params)

			t = image.get_attr("xform.projection")
			r = t.get_params("eman")
			failure = recon.insert_slice(image,t)
			
			progress += 1
			E2progress(logid,float(progress)/total_progress)
			ptcl_repr = 0
			if image.has_attr("ptcl_repr"):
				ptcl_repr = image.get_attr("ptcl_repr")
			if not failure: contributing_ptcls += ptcl_repr
			if (options.verbose):
				sys.stdout.write( "%2d/%d  %3d\t%5.1f  %5.1f  %5.1f\t\t%6.2f %6.2f" %
								(i,total_images, ptcl_repr,
								r["az"],r["alt"],r["phi"],
								image.get_attr("maximum"),image.get_attr("minimum")))
				if ( j > 0):
					sys.stdout.write("\t%f %f" %(recon.get_norm(idx), recon.get_score(idx) ))
					
				if ( failure ):
					sys.stdout.write( " X" )
					removed += 1
				
				sys.stdout.write("\n")

			idx += 1
		if options.verbose:	
			print "Iteration %d excluded %d images " %(j,removed)

	if (options.verbose):
		print "Inverting 3D Fourier volume to generate the real space reconstruction"
	output = recon.finish()
	output.set_attr("contributing_ptcls",contributing_ptcls)
	progress += 10
	E2progress(logid,float(progress)/total_progress)
	if (options.verbose):
		print "Finished Reconstruction"

	return output
	
def check(options,verbose=False):
	
	error = False
	
	if ( not options.sym ):
		if verbose:
			print  "Error: you must specify the sym argument"
		error = True
	
	if ( options.iter < 0 ):
		if verbose:
			print  "Error, --iter must be greater than or equal to 0"
		error = True
	
	if ( options.recon_type == None or options.recon_type == ""):
		if (verbose):
			print "Error: the --recon argument must be specified"
		error = True
	else:
		if ( check_eman2_type(options.recon_type,Reconstructors,"Reconstructor") == False ):
			error = True
	
	
	if options.keepsig and not options.keep:
		if verbose:
			print "Error: the --keepsig can only be supplied in tandem with the --keep= argument"
		error = True
	if (options.keep and not options.keepsig and ( options.keep > 1 or options.keep <= 0)):
		if verbose:
			print "Error: the --keep option is a percentage expressed as a fraction - it must be between 0 and 1. To switch to the sigma based alternative use the keepsig argument."
		error = True
	
	if ( options.nofilecheck == False ):
		if not file_exists(options.input):
			if verbose:
				print  "Error: the 2D image data file (%s) does not exist, cannot run e2make3d.py" %(options.input)
			error = True
		else:
			pad_dims = parse_pad(options.pad)
			if ( pad_dims != None ):
				if (len(pad_dims) > 2 or len(pad_dims)==0):
					if verbose:
						print  "Error, could not interpret the --pad argument, got",pad_dims
					error = True
				else:
					image_0 = get_processed_image(options,0,True)
					(xsize, ysize )= image_0.get_xsize(),image_0.get_ysize()
					if ( len(pad_dims) == 1 ) :
						pad = pad_dims[0]
						if pad < xsize or pad < ysize :
							if verbose:
								print  "Error, you specified a padding size (%d) that was smaller than the image dimensions (%dx%d)"%(pad,xsize,ysize)
							error = True;
					elif ( len(pad_dims) == 2 ) :
						xpad = pad_dims[0]
						ypad = pad_dims[1]
						if xpad < xsize :
							if verbose:
								print  "Error, you specified a x padding size (%d) that was smaller than the x dimension of the image (%d)"%(xpad,xsize)
							error = True
						if ypad < ysize :
							if verbose:
								print  "Error, you specified a y padding size (%d) that was smaller than the y dimension of the image (%d)"%(ypad,ysize)
							error = True
		
		if ( file_exists(options.output )):
			if ( not options.force ):
				if verbose:
					print  "Error: The output file exists and will not be written over. Use -f to overwrite"
				error = True
		
		
		if options.tlt != None:
			angles = None
			try:
				f=file(options.tlt,'r')
				lines=f.readlines()
				angles=[]
				for line in lines:
					angles.append(float(line))
					
			except:
				print "Error: the tlt file (%s) is un-interpretable." %options.tlt
				error = True
				
			if angles != None:
				# we assume that when the user is supplying a tlt file that the input is a 3D
				# image (such as an ali file). This could be different but I am rushed
				nx,ny,nz = gimme_image_dimensions3D(options.input)
				if nz != len(angles):
					print "Error: the tlt file (%s) has %d angles, but the image has %d slices. These two numbers need to be the same" %(options.tlt,len(angles),nz)
					error = True
			
	
	# if weighting is being used, this code checks to make sure atleast one image has an attribute "ptcl_repr" that is
	# greater than zero. If this is not the case, the insertion code will think there is nothing to insert...
	if ( options.no_wt == False and not options.nofilecheck):
		total_images=EMUtil.get_image_count(options.input)
		ptcl_repr = False;
		# A sanity test
		if not options.tlt:
			for i in xrange(0,total_images):
				image = EMData()
				read_header_only = True
				image.read_image(options.input, i, read_header_only)
				num_img=image.get_attr("ptcl_repr") 
				if (num_img > 0):
					ptcl_repr = True;
					break
		else:
			ptcl_repr = True
		
		if (ptcl_repr == False):
			print "Error - no image ptcl_repr attribute encountered that was greater than 0 - nothing done"
			print "Specify --no_wt to override this behaviour and give all inserted slices an equal weighting"
			error = True
	#if options.sym:
		#options.sym=options.sym.lower()
		#if (options.sym[0] in ["c","d", "h"]):
			#if not(options.sym[1:].isdigit()):
				#if verbose:
					#print "Error: %s is an invalid symmetry type. You must specify the --sym argument"%options.sym
				#error = True
		#else :
			#if not (options.sym in ["tet","oct","icos"]):
				#if verbose:
					#print "Error: %s is an invalid symmetry type. You must specify the --sym argument"%options.sym
				#error = True
	
	if options.preprocess != None:
		for p in options.preprocess:
			if ( check_eman2_type(p,Processors,"Processor") == False ):
				error = True
	if options.postprocess != None:
		for p in options.postprocess:
			if ( check_eman2_type(p,Processors,"Processor") == False ):
				error = True
	
	return error

if __name__=="__main__":
	main()
