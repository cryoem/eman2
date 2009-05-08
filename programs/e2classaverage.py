#!/usr/bin/env python

#
# Author: David Woolford, 9/7/2007 (woolford@bcm.edu)
# Copyright (c) 2000-2007 Baylor College of Medicine
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

from EMAN2 import *
from optparse import OptionParser
import math
from copy import deepcopy
import os
import sys

READ_HEADER_ONLY = True

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog <input particles> <class mx> <output> [options]

	Produces class averages.
	Can perform iterative alignment.
	Provides bootstrapping functionality.
	"""
		
	parser = OptionParser(usage=usage,version=EMANVERSION)
	
	parser.add_option("--iter", type="int", help="The number of iterations to perform. Default is 0.", default=0)
	parser.add_option("--ref", type="string", help="Reference image. If specified, the metadata in this image is used to assign euler angles to the generated classes. This is typically the projections that were used for the classification.", default=None)
	parser.add_option("--align",type="string",help="This is the aligner used to align particles to the previous class average. Default is None.", default=None)
	parser.add_option("--aligncmp",type="string",help="The comparitor used for the --align aligner. Default is dot.",default="phase")
	parser.add_option("--ralign",type="string",help="This is the second stage aligner used to refine the first alignment. This is usually the \'refine\' aligner.", default=None)
	parser.add_option("--raligncmp",type="string",help="The comparitor used by the second stage aligner.",default="phase")
	parser.add_option("--averager",type="string",help="The type of averager used to produce the class average.",default="mean")
	parser.add_option("--cmp",type="string",help="The comparitor used to generate quality scores for the purpose of particle exclusion in classes, strongly linked to the keep argument.", default="phase")
	parser.add_option("--keep",type="float",help="The fraction of particles to keep in each class.",default=1.0)
	parser.add_option("--keepsig", action="store_true", help="Causes the keep argument to be interpreted in standard deviations.",default=False)
	parser.add_option("--verbose","-v",action="store_true",help="Print useful information while the program is running. Default is off.",default=False)
	parser.add_option("--force", "-f",dest="force",default=False, action="store_true",help="Force overwrite the output file if it exists.")
	parser.add_option("--debug","-d",action="store_true",help="Print debugging infromation while the program is running. Default is off.",default=False)
	parser.add_option("--nofilecheck",action="store_true",help="Turns file checking off in the check functionality - used by e2refine.py.",default=False)
	parser.add_option("--check","-c",action="store_true",help="Performs a command line argument check only.",default=False)
	parser.add_option("--bootstrap",action="store_true",help="Bootstraps iterative alignment by using the first particle in each class to seed the iterative alignment. Only works if the number of iterations is greater than 0.")
	parser.add_option("--resultmx",type="string",help="Specify an output image to store the result matrix. This contains 5 images where row is particle number. Rows in the first image contain the class numbers and in the second image consist of 1s or 0s indicating whether or not the particle was included in the class. The corresponding rows in the third, fourth and fifth images are the refined x, y and angle (respectively) used in the final alignment, these are updated and accurate, even if the particle was excluded from the class.", default=None)
	parser.add_option("--normproc",type="string",help="The normalization processor. Default is normalize.edgemean. If you want to turn this option off specify \'None\'", default="normalize.edgemean")
	parser.add_option("--usefilt", dest="usefilt", default=None, help="Specify a particle data file that has been low pass or Wiener filtered. Has a one to one correspondence with your particle data. If specified will be used to align particles to the running class average, however the original particle will be used to generate the actual class average")
	parser.add_option("--idxcache", default=False, action="store_true", help="Stores the indices of all particles in the given class in the Python list in the e2classaverage.indices database")
	parser.add_option("--dbpath", help="Use in conjunction with --idxcache to specify a directory where the database entries should be stored, e.g. \"refine_01\" ", default=".")
	parser.add_option("--odd", default=False, help="Used by EMAN2 when running eotests. Includes only odd numbered particles in class averages.", action="store_true")
	parser.add_option("--even", default=False, help="Used by EMAN2 when running eotests. Includes only even numbered particles in class averages.", action="store_true")
	
	
	(options, args) = parser.parse_args()
	
	if len(args)<3 : parser.error("Input, classification matix, and output files required")
	
	if (options.check): options.verbose = True # turn verbose on if the user is only checking...
		
	if (options.nofilecheck == False):
		options.classifyfile=args[1]
		options.datafile = args[0]
		options.outfile = args[2]
	
	error = check(options,True)
	
	if (options.verbose):
		if (error):
			print "e2classaverage.py command line arguments test.... FAILED"
		else:
			print "e2classaverage.py command line arguments test.... PASSED"
	
	# returning a different error code is currently important to e2refine.py - returning 0 tells e2refine.py that it has enough
	# information to execute this script
	if error : exit(1)
	if options.check: exit(0)
	
	logger=E2init(sys.argv)
	
	# just remove the file - if the user didn't specify force then the error should have been found in the check function
	if os.path.exists(options.outfile):
		if options.force:
			remove_file(options.outfile)
	
	(num_classes, num_part ) = gimme_image_dimensions2D(args[1]);

	
	if (options.verbose):
		print "Classifications per particle %d, particles %d" %(num_classes, num_part)
	
	# classes contains the classifications - row is particle number, column data contains class numbers (could be greater than 1)
	classes = EMData()
	classes.read_image(args[1], 0)
	class_max = int(classes.get_attr("maximum"))
	class_min = int(classes.get_attr("minimum"))
	
	# NOT SURE IF THiS IS ALREADY IN THE CHECK FUNCTION FIXME 
	# double check that the argument reference image makes sense
	if (options.ref):
		if not os.path.exists(options.ref) and not db_check_dict(options.ref):
			parser.error("File %s does not exist" %options.ref)
			
		num_ref= EMUtil.get_image_count(options.ref)
		if ( class_max > num_ref ):
			print "Error, the classification matrix refers to a class number (%d) that is beyond the number of images (%d) in the reference image (%s)." %(class_max,num_ref,options.ref)
			exit(1)
	
	# double check that the number of particles in the particle image matches the rows in the classification matrix (image)
	num_part_check =  EMUtil.get_image_count(args[0])
	if ( num_part != num_part_check ):
		print "Error, the number of rows (%d) in the classification matrix (image) does not match the number of particles (%d) in the input image." %(num_part,num_part_check)
		exit(1)
	
	# weights contains the weighting of the classification scheme stored in the EMData object "classes" - above
	# dx contains the x translation of the alignment
	# dy contains the y translation of the alignment
	# da contains is the azimuthal rotation of the alignment
	# dflip contains is the mirror alignment
	# row is particle number, column data contains weights - rows should add to 1, but this is not checked.
	weights, dx, dy, da, dflip = EMData(),EMData(),EMData(),EMData(),EMData()
	if options.bootstrap:
		weights.set_size(classes.get_xsize(),classes.get_ysize())
		weights.to_zero()
		dx.set_size(classes.get_xsize(),classes.get_ysize())
		dx.to_zero()
		dy.set_size(classes.get_xsize(),classes.get_ysize())
		dy.to_zero()
		da.set_size(classes.get_xsize(),classes.get_ysize())
		da.to_zero()
		dflip.set_size(classes.get_xsize(),classes.get_ysize())
		dflip.to_zero()
	else:
		if EMUtil.get_image_count(args[1]) != 6:
			print "error, the classification matrix is the wrong size, it needs to contain one image for the classes, weights, dx, dy, da, and dflip. You can bypass this requirement if you supply the bootstrap argument"
			sys.exit(1)
			
		else:
			weights.read_image(args[1], 1)
			dx.read_image(args[1],2)
			dy.read_image(args[1],3)
			da.read_image(args[1],4)
			dflip.read_image(args[1],5)
	
	# empty space for storing x-flipping flags (0s or 1s, if a 1 is stored it will be used at the necessary point to flip prior to adding to the average)
	# sometimes this will not be used at all (it depends on whether or not the aligners that the user has specified do flipping and set flip flags)
#	dflip = EMData(da.get_xsize(),da.get_ysize())
#	dflip.to_zero()
	
	options.norm = parsemodopt(options.normproc)
	
	if (options.iter > 0 or options.bootstrap):
		set_aligner_params_in_options(options)
	
	if options.idxcache:
		try:
			name = "bdb:"+options.dbpath+"#class_indices"
			if options.even: name += "_even"
			elif options.odd: name += "_odd"
			db_name = numbered_bdb(name)
			class_db = db_open_dict(db_name)
		except:
			print "error with db terminology: can't call numbered_bdb with this argument:", "bdb:"+options.dbdir+"#class_indices"
			sys.exit(1)
	# do one class at a time
	for cl in range(class_min,class_max+1):
		E2progress(logger,float(cl-class_min)/(class_max+1-class_min))
		class_cache = ClassCache(args,options,num_part,num_classes,classes)
		class_cache.update_cache(cl)
		
		if (options.verbose): ndata = []
			
		if options.idxcache:
			class_indices = class_cache.class_indices

		if (options.iter > 0 or options.verbose):
			ccache = class_cache.ccache
		
		# this should work, if it doesn't it should have been caught by check function
		averager_parms=parsemodopt(options.averager)
		average = None
		image_cache = class_cache.image_cache
		filtered_image_cache = class_cache.filtered_image_cache
		if ( not options.bootstrap ):
			if options.verbose: print "generating the original class average using alignment parameters in the classification matrix"
			# generate the first class average by applying the transformations, adding and finally normalizing...		
			average,np = get_basic_average(averager_parms,class_cache,da,dflip,dx,dy,weights,options)
		else:
			# generate a bootstrapped initial average. Do this 'inductively' by aligning the 2nd image to the first, then adding. Continue until done...
			if options.verbose: print "bootstrapping the original class average"
			average,np = get_boot_strapped_average(class_cache,options,args)
		
		if options.idxcache:
			class_db[str(cl)] = class_indices
				
		if np == 0:
			if options.verbose:
				print "Class",cl,"...no particles"
			# FIXME
			continue
		
		if options.verbose: print "done original average"
	
		if (options.iter > 0):
			options.cull = True
			if ( options.keepsig == False and options.keep == 1.0 ) : options.cull = False
		
		for it in range(0,options.iter):
			itfrac=it/float(options.iter)
			# do alignment
			align_to_average_and_store_metadata(class_cache,options,average,dx,dy,da,dflip,weights)
#			
			# get the culling threshold
			if options.cull:
				qual_scores = []
				for d in ccache:
					p = d[0]
					c = d[1]
					qual_scores.append(weights.get(c,p))
				
				if ( options.keepsig ):
					a = Util.get_stats_cstyle(qual_scores)
					mean = a["mean"]
					std_dev = a["std_dev"]
					cullthresh = mean + (5.0*(1.0-itfrac)+itfrac*options.keep)*std_dev
				else:
					b = deepcopy(qual_scores)
					b.sort()
					# The ceil reflects a conservative policy. If the user specified keep=0.93
					# and there were 10 particles, then they would all be kept. If floor were
					# used instead of ceil, the last particle would be thrown away (in the
					# class average)
					idx = int(math.ceil(((1.0-itfrac)+itfrac*options.keep)*len(b))-1)
					cullthresh = b[idx]
			
			#finally average with culling
			average,np = get_basic_average_with_cull(averager_parms,class_cache,options,weights,da,dx,dy,dflip,ndata,cullthresh)			
			if options.verbose: print "done iteration",it

		if average == None: continue
		# Align to the reference and extract euler data if it was specified
		if ( options.ref):
			e = EMData()
			if options.iter > 0:
				e.read_image(options.ref, cl)
				
				average_to_ref = align(e,average,options)
				averager=Averagers.get(averager_parms[0], averager_parms[1]) # need an empty averager
				fine_transform = average_to_ref.get_attr("xform.align2d")
				fine_transform.invert()
				
				#average,np = get_basic_average_with_tweak(averager_parms,class_cache,options,weights,da,dx,dy,dflip,fine_transform)
				#avg = None
				np = 0 # number of particles in the average
				for d in ccache:
					p = d[0]
					c = d[1]
					if (options.cull ):
						if ( weights.get(c,p) == 0 ) : continue
					
#					if (options.lowmem):
	   	   	   	   	image = image_cache[p]
					#image = EMData(args[0],p)
#					else:
#						image = images[p].copy()
					
	#				if str(options.normproc) != "None": image.process_inplace(options.norm[0],options.norm[1])
	
					t = Transform({"type":"2d","alpha":da.get(c,p)})
					t.set_trans(dx.get(c,p),dy.get(c,p))
					if dflip.get(c,p) != 0: t.set_mirror(True)
					ct = fine_transform*t
					#print "final transform is"
					#ct.printme()
					#fine_transform.printme()
					rslt = image.process("math.transform",{"transform":ct})
					
					np += 1
					averager.add_image(rslt)
					
					# Don't forget to store the alignment parameters
					params = ct.get_params("2d")
					dx.set(c,p, params["tx"])
					dy.set(c,p, params["ty"])
					da.set(c,p, params["alpha"])
					dflip.set(c,p, params["mirror"])
				
				if np == 0:
					if (options.verbose):
						print "Class",cl,"...no particles on iteration",it
					continue
#				
				average = averager.finish()
#				average.transform(fine_transform)
				
				if str(options.normproc) != "None": average.process_inplace(options.norm[0],options.norm[1])
				#average.process_inplace("xform.centerofmass") this shouldn't be necessary if we aligned to the projection
				average.process_inplace("mask.sharp",{"outer_radius":ta.get_xsize()/2})
			else:
				e.read_image(options.ref, cl, READ_HEADER_ONLY)

			average.set_attr("xform.projection", e.get_attr("xform.projection"))
			average.set_attr("projection_image",options.ref)
			average.set_attr("projection_image_idx",cl)
			#average.set_attr("euler_az",e.get_attr("euler_az"))
			#average.set_attr("euler_phi",e.get_attr("euler_phi"))
			if options.verbose:
				edata = []
				t = e.get_attr("xform.projection")
				rot = t.get_rotation("eman")
				edata.append(rot["alt"])
				edata.append(rot["az"])
				edata.append(rot["phi"])
		
		# now write to disk
		average.set_attr("ptcl_repr",np)
		if np == 0:
			a = average.copy() # get all the attributes
			a.to_zero()
			a.write_image(args[2],-1)
		else:
			average.write_image(args[2],-1)
			
		if options.verbose:
			sys.stdout.write( "Class %d: particles..%d" %(cl,len(ccache)) )
			for t in range(0,options.iter):
				sys.stdout.write("..%d" %ndata[t] )
			if ( options.ref  ):
				sys.stdout.write(" : Eulers..")
				for t in edata:
					sys.stdout.write(" %f" %t)
			sys.stdout.write("\n")
		
	if (options.resultmx != None ):
		if os.path.exists(options.resultmx):
			remove_file(options.resultmx) #oooh not sure about this!
		
		# note the order is important!
		classes.write_image(options.resultmx,-1)
		weights.write_image(options.resultmx,-1)
		dx.write_image(options.resultmx,-1)
		dy.write_image(options.resultmx,-1)
		da.write_image(options.resultmx,-1)
		dflip.write_image(options.resultmx,-1)
		
	if options.idxcache:
			name = "bdb:"+options.dbpath+"#class_indices"
			if options.even: name += "_even"
			elif options.odd: name += "_odd"
			db_name = numbered_bdb(name)
			#db_close_dict(db_name)
		
	
	E2end(logger)
	
	
def align_to_average_and_store_metadata(class_cache,options,average,dx,dy,da,dflip,weights):
	if options.usefilt: image_cache = class_cache.filtered_image_cache
	else: image_cache = class_cache.image_cache
	for p,image in image_cache.items():
		
		c = class_cache.idx_cache[p]
		ta = align(average,image,options)
		t = ta.get_attr("xform.align2d")
		t.invert() # because we aligned the average to the image, and not the other way around
		params = t.get_params("2d")
		
		# store the refined translational and rotational values
		dx.set(c,p, params["tx"])
		dy.set(c,p, params["ty"])
		da.set(c,p, params["alpha"])
		dflip.set(c,p, params["mirror"])
		
		# store the quality score on top of the weights, seeing as the starting weights are no longer required
		if (options.cull): # but only if we need to
			weights.set(c,p, ta.cmp(options.cmp[0],image,options.cmp[1]))

def get_basic_average_with_tweak(averager_parms,class_cache,options,weights,da,dx,dy,dflip,fine_transform):
	averager=Averagers.get(averager_parms[0], averager_parms[1]) # need an empty averager
	np = 0 # number of particles in the average
	for p,image in class_cache.image_cache.items():
		c = class_cache.idx_cache[p]
		if (options.cull ):
			if ( weights.get(c,p) == 0 ) : continue
		
		t = Transform({"type":"2d","alpha":da.get(c,p)})
		t.set_trans(dx.get(c,p),dy.get(c,p))
		if dflip.get(c,p) != 0: t.set_mirror(True)
		ct = fine_transform*t

		rslt = image.process("math.transform",{"transform":ct})
		
		np += 1
		averager.add_image(rslt)
		
		# Don't forget to store the alignment parameters
		params = ct.get_params("2d")
		dx.set(c,p, params["tx"])
		dy.set(c,p, params["ty"])
		da.set(c,p, params["alpha"])
		dflip.set(c,p, params["mirror"])
	
	if np == 0:
		if (options.verbose):
			print "Class",cl,"...no particles on iteration",it
		average = None
	else:	
		average = averager.finish()
		#				average.transform(fine_transform)
		
		if str(options.normproc) != "None": average.process_inplace(options.norm[0],options.norm[1])
		#average.process_inplace("xform.centerofmass") this shouldn't be necessary if we aligned to the projection
		average.process_inplace("mask.sharp",{"outer_radius":ta.get_xsize()/2})
		
	return average,np

def get_basic_average_with_cull(averager_parms,class_cache,options,weights,da,dx,dy,dflip,ndata,cullthresh):
	averager=Averagers.get(averager_parms[0], averager_parms[1]) # need an empty averager
	#avg = None
	np = 0 # number of particles in the average
	for p,image in class_cache.image_cache.items():
		c = class_cache.idx_cache[p]
		if (options.cull ):
			# if we're culling then find the bad particles here...
			# Also store whether or not culling occured by placing 1s and 0s in the weights object
			# This is in case users ever want to know which particles were excluded. This information
			# is written to disk if the resultmx option is specfied
			# FIXME setting the weights matrix to 1 and 0 should probably only occur if we're at the
			# last iteration
			if ( weights.get(c,p) >= cullthresh ) :
				weights.set(c,p,0)
				continue
			else: weights.set(c,p,1)
		else: weights.set(c,p,1)

		t = Transform({"type":"2d","alpha":da.get(c,p)})
		t.set_trans(dx.get(c,p),dy.get(c,p))
		if dflip.get(c,p) != 0: t.set_mirror(True)

		rslt = image.process("math.transform",{"transform":t})
	
		np += 1
		averager.add_image(rslt)
		
	if options.verbose:
		ndata.append(np)

	if np == 0:
		if (options.verbose):
			print "Class",cl,"...no particles on iteration",it
		# FIXME
		average = None
	else:

		average = averager.finish()
		#should this be centeracf?
		if str(options.normproc) != "None": average.process_inplace(options.norm[0],options.norm[1])
		average.process_inplace("xform.centerofmass")
		average.process_inplace("mask.sharp",{"outer_radius":average.get_xsize()/2})
	
	return average,np
			

def get_boot_strapped_average(class_cache,options,args):
	average = None
	np = 0
	
	for p,image in class_cache.image_cache.items():
		c = class_cache.idx_cache[p]

				
		if (average == None):
			average = EMData()
			average.read_image(args[0],p)

			if str(options.normproc) != "None": average.process_inplace(options.norm[0],options.norm[1])
			average.process_inplace("mask.sharp",{"outer_radius":average.get_xsize()/2})
			#average.write_image("ave_.hdf",-1)
			np = 1
		else:
			# there are a lot of ways to do and there is probably
			# an improvement on what happens here, but it suffices
			if not options.usefilt:
				ta = align(average,image,options)
				t = ta.get_attr("xform.align2d")
				t.invert()
				ta = image.copy()
				ta.process_inplace("math.transform",{"transform":(t)})
			else:
				filt_image = class_cache.filtered_image_cache[p]
					
				ta = align(average,filt_image,options)
				t = ta.get_attr("xform.align2d")
				t.invert()
				ta = image.process("math.transform",{"transform":(t)})
				ta.process_inplace("mask.sharp",{"outer_radius":ta.get_xsize()/2})				
				np += 1
				average.add(ta) # now add the image
	
	if np != 0:
		if str(options.normproc) != "None": average.process_inplace(options.norm[0],options.norm[1])
		average.process_inplace("xform.centeracf")
		average.process_inplace("mask.sharp",{"outer_radius":average.get_xsize()/2})
	
	return average,np

def get_basic_average(averager_parms,class_cache,da,dflip,dx,dy,weights,options):
	averager=Averagers.get(averager_parms[0], averager_parms[1])
	# do the initial average, based on the program inputs
	weightsum = 0 # used to normalize the average
	np = 0 # number of particles in the average
	
	for p,image in class_cache.image_cache.items():
		c = class_cache.idx_cache[p]

		t = Transform({"type":"2d","alpha":da.get(c,p),"mirror":int(dflip(c,p))})
		t.set_trans(dx.get(c,p),dy.get(c,p))
		
		rslt = image.process("math.transform",{"transform":t})
				
		np += 1
		weight = weights(c,p)
		weightsum += weight
		rslt.mult(weight)
		
		# Add the image to the averager
		averager.add_image(rslt)
	
	if np == 0 or weightsum == 0:
		if options.verbose:
			print "Class",cl,"...no particles"
		# FIXME
		
		# write blank image? Write meta data?
		return None
		
	average = averager.finish()
	average.mult(float(np)) # Undo the division of np by the averager - this was incorrect because the particles were weighted.
	average.mult(1.0/weightsum) # Do the correct division
	average.process_inplace("xform.centeracf")
	
	return average,np

class ClassCache:
	
	'''
	Stores useful information at each iteration
	'''
	def __init__(self,args,options,num_part,num_classes,class_data):
		self.args = args
		self.options = options
		self.num_part = num_part
		self.num_classes = num_classes
		self.class_data = class_data
		self.image_cache = {}
		self.filtered_image_cache= {}
		self.idx_cache = {}
		self.class_indices = []
		self.ccache = []
		
	def update_cache(self,class_number):
		
		for p in xrange(0,self.num_part):
			if self.options.odd and p % 2 == 0: continue # enforce even/odd constraints
			if self.options.even and p % 2 == 1: continue # enforce even/odd constraints
			for c in xrange(0,self.num_classes):
				if self.class_data.get(c,p) == class_number:
					
					if self.options.idxcache: self.class_indices.append(p)
					# cache the hit if necessary - this is if there is more than one iteration, and/or the user has specifed verbose. In the
					# latter case the class cache is used to print information
					if (self.options.iter > 0 or self.options.verbose): self.ccache.append((p,c))
					
					self.idx_cache[p] = c
					
					image = EMData()
					image.read_image(self.args[0],p)
					if str(self.options.normproc) != "None": image.process_inplace(self.options.norm[0],self.options.norm[1])
					self.image_cache[p] = image
					
					if self.options.usefilt:
						filt_image = EMData()
						filt_image.read_image(options.usefilt,p)
						if str(self.options.normproc) != "None": filt_image.process_inplace(self.options.norm[0],self.options.norm[1])
						self.filtered_image_cache[p] = filt_image

						

def set_aligner_params_in_options(options):
	'''
	Call this before calling align
	'''
	options.align=parsemodopt(options.align)
	options.alicmp=parsemodopt(options.aligncmp)
	# note the parsing of the options.ralign parameters is left for later
	options.alircmp=parsemodopt(options.raligncmp)
	options.cmp=parsemodopt(options.cmp)

def align(this,to,options):
	# Align the particle to the average
	ta=this.align(options.align[0],to,options.align[1],options.alicmp[0],options.alicmp[1])
	
	if ( options.ralign != None ): # potentially employ refine alignment
		
		refineparms=parsemodopt(options.ralign)
		# this parameters I think west best with the refine aligner, but they're not well tested
		# i.e. I need to do some rigorous testing before I can claim anything
		#refineparms[1]["az"] = ta.get_attr_default("align.az",0)-1
		#refineparms[1]["dx"] = ta.get_attr_default("align.dx",0)-1
		#refineparms[1]["dy"] = ta.get_attr_default("align.dy",0)-1
		#refineparms[1]["mode"] = 0
		#refineparms[1]["stepx"] = 2
		#refineparms[1]["stepy"] = 2
		#refineparms[1]["stepaz"] = 5
		
		refineparms[1]["xform.align2d"] = ta.get_attr("xform.align2d")
		#refineparms[1]["dx"] = ta.get_attr_default("align.dx",0)
		#refineparms[1]["dy"] = ta.get_attr_default("align.dy",0)
		#flip = ta.get_attr_default("align.flip",0)
		
		#refine_this = this
		#if flip:
			#refine_this = this.process("xform.flip",{"axis":"x"})
		
		ta = this.align(refineparms[0],to,refineparms[1],options.alircmp[0],options.alircmp[1])
		
		#ta.set_attr("align.flip",flip)
	return ta
		
def check(options, verbose=False):

	error = False
	if ( options.nofilecheck == False ):
		
		if os.path.exists(options.outfile):
			if not options.force:
				error = True
				if (verbose):
					print "Error: output file %s exists, force not specified, will not overwrite, exiting" %options.outfile
		
		if not os.path.exists(options.classifyfile) and not db_check_dict(options.classifyfile):
			error = True
			if (verbose):
				print "Error: the file expected to contain the classification matrix (%s) was not found, cannot run e2classaverage.py" %(options.classifyfile)
		
		if not file_exists(options.datafile) and not db_check_dict(options.datafile):
			error = True
			if (verbose):
				print "Error:  failed to find the particle data (%s)" %options.datafile
		else:
			if (options.usefilt != None):
				if not file_exists(options.usefilt):
					error = True
					if verbose: print "Error: failed to find usefilt file %s" %options.usefilt
				
				n1 = EMUtil.get_image_count(options.usefilt)
				n2 = EMUtil.get_image_count(options.datafile)
				if n1 != n2:
					if verbose: print "Error, the number of images in the starting particle set:",n2,"does not match the number in the usefilt set:",n1
					error = True
					
				read_header_only=True
				img1 = EMData()
				img1.read_image(options.datafile,0,read_header_only)
				img2 = EMData()
				img2.read_image(options.usefilt,0,read_header_only)
				
				nx1 = img1.get_attr("nx") 
				nx2 = img2.get_attr("nx") 
				
				ny1 = img1.get_attr("ny") 
				ny2 = img2.get_attr("ny") 
				
				if nx1 != nx2 or ny1 != ny2:
					error = True
					if verbose: print "Error, the dimensions of particle data (%i x %i) and the usefilt data (%i x %i) do not match" %(nx1,ny1,nx2,ny2)
				
		
		if os.path.exists(options.classifyfile) and os.path.exists(options.datafile):
			(xsize, ysize ) = gimme_image_dimensions2D(options.classifyfile);
			numimg = EMUtil.get_image_count(options.datafile)
			if ( numimg != ysize ):
				error = True
				if (verbose):
					print "Error - the number of rows (%d) in the classification matrix image %s does not match the number of images (%d) in %s" %(ysize, options.classifyfile,numimg,options.datafile)
				
			
		if options.ref != None and not os.path.exists(options.ref) and not db_check_dict(options.ref):
			print "Error: the file expected to contain the reference images (%s) does not exist" %(options.ref)
			error = True
		elif options.ref and (os.path.exists(options.datafile) or db_check_dict(options.datafile)):
			(xsize, ysize ) = gimme_image_dimensions2D(options.datafile);
			(pxsize, pysize ) = gimme_image_dimensions2D(options.ref);
			if ( xsize != pxsize ):
				error = True
				if (verbose):
					print "Error - the dimensions of the reference and particle images do not match"

	if (options.iter > 1 or options.bootstrap) and options.align == None:
		print "Error: you must specify the align argument"
		error = True
		
	if ( options.even and options.odd ):
		print "Error, the even and odd arguments are mutually exclusive"
		error = True
	#if ( options.keep and options.keepsig ):
		#error = True
		#if ( verbose ):
			#print "Error: --keep and --keepsig are mutually exclusive"
	
	if ( options.keep > 1 or options.keep <= 0) and not options.keepsig :
		error = True
		if (verbose):
			print "The --keep option is a percentage expressed as a fraction - it must be between 0 and 1"
	
	if ( options.iter < 0 ):
		error = True
		if (verbose):
			print "Error, --iter must be greater than or equal to 0 - you specified %d" %(options.iter)
		
	if ( check_eman2_type(options.averager,Averagers,"Averager") == False ):
		if (verbose):
			print "Unknown averager",options.averager
		error = True
	
	if ( options.iter > 0 ):
		
		if ( check_eman2_type(options.cmp,Cmps,"Comparitor") == False ):
			error = True
			
		if (options.align == None):
			print "If --classiter is greater than zero, the -align argument must be specified"
			error = True
			
		if ( check_eman2_type(options.align,Aligners,"Aligner") == False ):
			error = True

		if ( check_eman2_type(options.aligncmp,Cmps,"Comparitor") == False ):
			error = True
		
		if ( options.ralign != None ):
			if ( check_eman2_type(options.ralign,Aligners,"Aligner") == False ):
				error = True
				
			if ( check_eman2_type(options.raligncmp,Cmps,"Comparitor") == False ):
				error = True
		if ( str(options.normproc) != "None" ):
			if ( check_eman2_type(options.normproc,Processors,"Processor") == False ):
				error = True

	return error
	
if __name__ == "__main__":
    main()
