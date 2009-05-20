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

from EMAN2db import EMTask

class EMGenClassAverages:
	'''
	A utility class that knows how to break the command line into EMClassAveTasks
	Use like this:
	job = EMGenClassAverages(options,args,logger)
	job.execute()
	'''
	def __init__(self,options,args,logger=None):
		'''
		@param options the options produced by (options, args) = parser.parse_args()
		@param args the options produced by (options, args) = parser.parse_args()
		assumes you have already called the check function.
		'''
		self.options = options
		self.args = args
		self.logger = logger
		
		self.__task_options = None # will eventually be the options parameter of the EMDCClassAverageTask
		
	def __get_task_options(self,options):
		'''
		Get the options required by each task as a dict
		@param options is always self.options - the initialization argument. Could be changed.
		'''
		if self.__task_options == None:
			d = {}
			d["iter"] = options.iter
			d["align"] = parsemodopt(options.align)
			d["aligncmp"] = parsemodopt(options.aligncmp)
			d["cmp"] = parsemodopt(options.cmp)
			d["averager"] = parsemodopt(options.averager)
			
			if hasattr(options,"ralign") and options.ralign != None: 
				d["ralign"] = parsemodopt(options.ralign)
				d["raligncmp"] = parsemodopt(options.raligncmp)  # raligncmp must be specified if using ralign
			else: 
				d["ralign"] = None
				d["raligncmp"] = None
			
			d["keep"] = options.keep
			
			if hasattr(options,"keepsig") and options.keepsig != None: d["keepsig"] = options.keepsig 
			else: d["keepsig"] = False
			
			if hasattr(options,"bootstrap") and options.bootstrap != None: d["bootstrap"] = options.bootstrap 
			else: d["bootstrap"] = False
			
			if hasattr(options,"normproc") and options.normproc != None: d["normproc"] = parsemodopt(options.normproc)
			else: d["normproc"] = None
			
			if hasattr(options,"verbose") and options.verbose != None: d["verbose"] = options.verbose
			else: d["verbose"] = None
			
			self.__task_options = d
			
		return self.__task_options
	
	def __get_class_data(self,class_number,options):
		'''
		Return a dict and a list with useful class data in them
		The dict's keys are the particle index (image row), and dict's values are the class number (image col)
		The list contains the particle indices
		self.__init_memory must be called before this function is called
		@param class_number the class number
		@param options is always self.options - the initialization argument. Could be changed.
		'''
		ptcl_indices = [] # this stores a particle index, for each partice in the class
		dcol_idx_cache = {}# a dictionary version of the col_idx_cache 
		for p in xrange(0,self.num_part):
			if options.odd and p % 2 == 0: continue # enforce even/odd constraints
			if options.even and p % 2 == 1: continue # enforce even/odd constraints
			for c in xrange(0,self.num_classes):
				if self.classes.get(c,p) == class_number:
					
					ptcl_indices.append(p)
					
					dcol_idx_cache[p] = c
		
		return ptcl_indices,dcol_idx_cache
		
	def __get_alignment_data(self,dcol_idx_cache):
		'''
		Get the alignment data out of the input alignment matrix
		Returns a dict, the keys are the particle index, the values are Transforms storing the initial alignment
		@param dcol_idx_cachhe - the dictionary returned by a call to self.__get_class_data
		'''
		transforms = {}
		for ptcl_idx,col in dcol_idx_cache.items():
			t = Transform({"type":"2d","alpha":self.da.get(col,ptcl_idx),"mirror":int(self.dflip(col,ptcl_idx))})
			t.set_trans(self.dx.get(col,ptcl_idx),self.dy.get(col,ptcl_idx))
			transforms[ptcl_idx] = t
		return transforms
	
	
	def __store_alignment_data(self,alis,dcol_idx_cache):
		'''
		Store alignment data in internal images which can later be written to disk
		alis and dcol_idx_cache should have the same keys. This is assumed
		@param dcol_idx_cachhe - the dictionary returned by a call to self.__get_class_data
		@param alis - a dictionary, keys are particle indices, values are Transforms storing 2D alignments
		'''
		for ptcl_idx,col in dcol_idx_cache.items():
			transform = alis[ptcl_idx]
			params = transform.get_params("2d")
			self.dx.set(col,ptcl_idx, params["tx"])
			self.dy.set(col,ptcl_idx, params["ty"])
			self.da.set(col,ptcl_idx, params["alpha"])
			self.dflip.set(col,ptcl_idx, params["mirror"])
			
	def __store_weight_data(self,weights,dcol_idx_cache):
		'''
		Store weight data in an internal image which can later be written to disk
		weights and dcol_idx_cache should have the same keys. This is assumed
		@param dcol_idx_cachhe - the dictionary returned by a call to self.__get_class_data
		@param weights - a dictionary, keys are particle indices, values are single values (ints, floats, bools are fine)
		'''
		for ptcl_idx,col in dcol_idx_cache.items():
			val = weights[ptcl_idx]
			self.weights.set(col,ptcl_idx, float(val))
	
	def __get_weight_data(self,dcol_idx_cache):
		'''
		Get the particle weights out of the input image matrix
		Returns a dict, the keys are the particle index, the values are weights
		@param dcol_idx_cachhe - the dictionary returned by a call to self.__get_class_data
		'''
		weights = {}
		for ptcl_idx,col in dcol_idx_cache.items():
			weight = self.weights.get(col,ptcl_idx)
			weights[ptcl_idx] = weight
		return weights				
	
	def __init_memory(self,args,options):
		'''
		Called internally to read the alignment and classification images into memory
		Also opens a local database for storing classification data, if the options has the idxcache attribute 
		@param args - invariably self.args - see initializer
		@param options - invariably self.options - see initializer
		'''
		# classes contains the classifications - row is particle number, column data contains class numbers (could be greater than 1)
		self.classes = EMData()
		self.classes.read_image(args[1], 0)
		(self.num_classes, self.num_part ) = gimme_image_dimensions2D(args[1]);
		self.class_max = int(self.classes.get_attr("maximum"))
		self.class_min = int(self.classes.get_attr("minimum"))
	
		# weights contains the weighting of the classification scheme stored in the EMData object "classes" - above
		# dx contains the x translation of the alignment
		# dy contains the y translation of the alignment
		# da contains is the azimuthal rotation of the alignment
		# dflip contains is the mirror alignment
		# row is particle number, column data contains weights - rows should add to 1, but this is not checked.
		self.weights, self.dx, self.dy, self.da, self.dflip = EMData(),EMData(),EMData(),EMData(),EMData()
		if options.bootstrap:
			for image in [self.weights,self.dx,self.dy,self.da,self.dflip]:
				image.set_size(self.classes.get_xsize(),self.classes.get_ysize())
				image.to_zero()
			
		elif EMUtil.get_image_count(args[1]) != 6:
			print "error, the classification matrix is the wrong size, it needs to contain one image for the classes, weights, dx, dy, da, and dflip. You can bypass this requirement if you supply the bootstrap argument"
			sys.exit(1)
		else:
			self.weights.read_image(args[1], 1)
			self.dx.read_image(args[1],2)
			self.dy.read_image(args[1],3)
			self.da.read_image(args[1],4)
			self.dflip.read_image(args[1],5)
			
		if options.idxcache:
			try:
				name = "bdb:"+options.dbpath+"#class_indices"
				if options.even: name += "_even"
				elif options.odd: name += "_odd"
				self.db_name = numbered_bdb(name)
				self.class_db = db_open_dict(self.db_name)
			except:
				print "error with db terminology: can't call numbered_bdb with this argument:", "bdb:"+options.dbpath+"#class_indices"
				sys.exit(1)
				
	def execute(self):
		'''
		Make every class average
		Write results to disk
		The main function
		'''
		
		if hasattr(self.options,"parallel") and self.options.parallel != None:
			
			self.task_customers = []
			self.tids = []
			self.__init_memory(self.args, self.options)
			for class_idx in xrange(self.class_min,self.class_max+1):
				ptcl_indices,dcol_idx_cache =  self.__get_class_data(class_idx, self.options)
				if self.options.idxcache:
					self.class_db[str(class_idx)] = ptcl_indices

				init_alis = None
				init_weights = None
				if not self.options.bootstrap:
					init_alis =  self.__get_alignment_data(dcol_idx_cache)
					init_weights = self.__get_weight_data(dcol_idx_cache)
#					for t in init_alis.values(): print t
#					for t in init_weights.values(): print t
					
				data = {}
				data["input"] = ("cache",self.args[0],ptcl_indices)
				if self.options.usefilt:
					data["usefilt"] = ("cache",options.usefilt,ptcl_indices)
				if init_alis != None:
					data["init_alis"] = init_alis
				if init_weights != None:
					data["weights"] = init_weights
				data["class_idx"] = class_idx 
				
				if hasattr(self.options,"ref") and self.options.ref != None:
					data["ref"] = ("cache",self.options.ref,class_idx)
				
				task = EMClassAveTaskDC(data=data,options=self.__get_task_options(self.options))
				
				from EMAN2PAR import EMTaskCustomer
				etc=EMTaskCustomer("dc:localhost:9990")
				#print "Est %d CPUs"%etc.cpu_est()
				tid=etc.send_task(task)
				#print "Task submitted tid=",tid
				
				self.task_customers.append(etc)
				self.tids.append(tid)

			
			classes = range(self.class_min,self.class_max+1)
			while 1:
				if len(self.task_customers) == 0: break
				print len(self.task_customers),"tasks left in main loop"
				for i in xrange(len(self.task_customers)-1,-1,-1):
					task_customer = self.task_customers[i]
					tid = self.tids[i] 
					st=task_customer.check_task((tid,))[0]
					if st==100:
						
						self.task_customers.pop(i)
						self.tids.pop(i)

						rslts = task_customer.get_results(tid)
						
						
						self.__write_class_data(rslts[1])
#						for image in rslts[1]["averages"]:
#							if image != None: image.write_image("dc_average.hdf",-1)
#							
#						if hasattr(self.options,"ref") and self.options.ref != None:
#							classes.remove(rslts[1]["class_idx"])
#							a = EMData(self.options.ref,rslts[1]["class_idx"])
#							a.write_image("dc_average.hdf",-1)		
						if self.logger != None:
							E2progress(self.logger,1.0-len(self.task_customers)/(self.class_max+1-self.class_min))
				
				time.sleep(5)
				
		else: # process in the context of the current program
		
			self.__init_memory(self.args, self.options)
			for class_idx in xrange(self.class_min,self.class_max+1):
				ptcl_indices,dcol_idx_cache =  self.__get_class_data(class_idx, self.options)
				if self.options.idxcache:
					self.class_db[str(class_idx)] = ptcl_indices

				init_alis = None
				init_weights = None
				if not self.options.bootstrap:
					init_alis =  self.__get_alignment_data(dcol_idx_cache)
					init_weights = self.__get_weight_data(dcol_idx_cache)
#					for t in init_alis.values(): print t
#					for t in init_weights.values(): print t
					
				data = {}
				data["input"] = (self.args[0],ptcl_indices)
				if self.options.usefilt:
					data["usefilt"] = (options.usefilt,ptcl_indices)
				if init_alis != None:
					data["init_alis"] = init_alis
				if init_weights != None:
					data["weights"] = init_weights
				data["class_idx"] = class_idx 
				
				if hasattr(self.options,"ref") and self.options.ref != None:
					data["ref"] = (self.options.ref,class_idx)
				
				task = EMClassAveTask(data=data,options=self.__get_task_options(self.options))
				
				print "generating class",class_idx
				task.execute()
				rslt = task.get_return_data()
				self.__write_class_data(rslt)
				
		self.__finalize_writing()
		
	def __finalize_writing(self):
		'''
		Write the results matrix to disk
		'''
		if hasattr(self.options, "resultmx") and self.options.resultmx != None:
			# note the order is important!
			self.classes.write_image(self.options.resultmx,-1)
			self.weights.write_image(self.options.resultmx,-1)
			self.dx.write_image(self.options.resultmx,-1)
			self.dy.write_image(self.options.resultmx,-1)
			self.da.write_image(self.options.resultmx,-1)
			self.dflip.write_image(self.options.resultmx,-1)
			

	def __write_class_data(self,rslts):
		'''
		Write class image to disk
		Store alignment and inclusion metadata in internal images
		@param rslts a dictionary that was returned by an EMClassAveTask or an EMClassAveTaskDC 
		'''
		average = rslts["final_average"]
		if average != None:
			if hasattr(self.options,"ref") and self.options.ref != None:
				average.set_attr("projection_image",self.options.ref)
			average.write_image(self.args[2],-1)
			ptcl_indices,dcol_idx_cache =  self.__get_class_data(rslts["class_idx"], self.options)
			
			final_alis = rslts["final_alis"]
			if  final_alis != None: self.__store_alignment_data(final_alis,dcol_idx_cache)
			final_incs = rslts["final_inclusions"]
			if final_incs != None:
				self.__store_weight_data(final_incs,dcol_idx_cache)
	
	
class EMClassAveTask(EMTask):
	'''
	This is the class average task used to generate class average sequentially in a single program
	Its build to be generic - you can inherit from it and redefine your own init_memory function,
	As in the case of the EMClassAveTaskDC
	'''

	def __init__(self,command="e2classaverage",data=None,options=None):
		EMTask.__init__(self,command,data,options)
		self.averages = [] # will eventually store a list of averages.
		self.all_alis = [] # will eventually contain all of the alignment parameters of each iteration
		self.all_inclusions = [] # will eventually contain dicts of flags for each image (True or False) indicating inclusion or nor
		self.final_average = None
		self.final_alis = None # will be a dictionary of 
		self.class_idx = data["class_idx"] # so it's easy to tell the calling function which class this is
		# options should have these keys:
		# iter - the total number of iterations. 0 is fine
		# align - the main aligner, a list of two strings
		# alligncmp - the main align cmp - a list of two strings
		# ralign - the refine aligner, a list of two string. May be None which turns it off
		# raligncmp - the refinealigncmp - a list of two strings. Needs to specified if ralign is not None
		# averager - the averager - a list of two strings
		# cmp - the final cmp - a list of two strings
		# keep - keep argument, interpreted as a percentage or a number of sigmas depending on the keepsig argument
		# keepsig - if True turns the keep into a sigma based threshold. May be None, False, or unspecified
		# normproc - A normalization processor - a list of two strings. May be None or unspecified	
		# debug - True, False, None, or unspecified - If True you get extra information returned from get_return_data
		# verbose - True, False, None,0,1, or unspecified - If True extra information is printed to stdout
		# bootstrap - True, False, None, or unspecified - If True original average is generated using a bootstrapping approache
	
		# data should have
		# input - a Task-style cached list of input images
		# usefilt - a Task-style cached list of usefilt images
		# ref - a Task-style cached single image
		# class_idx - The class index, a number
		
	def init_memory(self):
		'''
		Reads images into memory, normalize them if the options dictate it
		'''
		raw_data_name=self.data["input"][0]
		self.ptcl_indices = self.data["input"][1]
		
		self.images = {}
		if self.data.has_key("usefilt") and self.data["usefilt"] != None:
			self.usefilt_images = {}
			usefilt_data_name=self.data["usefilt"][0]
		else:
			self.usefilt_images = None
		
		norm = None
		if self.options.has_key("normproc") and self.options["normproc"] != None:
			norm = self.options["normproc"] # a list of length two
		self.norm = norm
			
		for ptcl_idx in self.ptcl_indices:
			image = EMData(raw_data_name,ptcl_idx)
			if norm != None:
				image.process_inplace(norm[0],norm[1])
			self.images[ptcl_idx] = image
			
			if self.usefilt_images != None:
				usefilt_image = EMData(usefilt_data_name,ptcl_idx)
				if norm != None:
					usefilt_image.process_inplace(norm[0],norm[1])
				self.usefilt_images[ptcl_idx] = usefilt_image
		# set a flag that can be used to decide if quality scores need to be
		# evaluated and used for culling the worst, or least similar, particles
		# from the class
		self.culling = False
		if (self.options["iter"] > 0):
			self.culling = True
			if (self.options.has_key("keepsig") and self.options["keepsig"] == False) and self.options["keep"] == 1.0 : self.culling=False
		
		# the reference is used to perform a final alignment, if specified
		self.ref = None
		if self.data.has_key('ref') and self.data['ref'] != None:
			ref_data_name=self.data["ref"][0]
			idx = self.data["ref"][1]
			self.ref = EMData(ref_data_name,idx)
		
		self.verbose = 0
		if self.options.has_key("verbose") and self.options["verbose"] != None:
			self.verbose = self.options["verbose"]
		
	def execute(self):
		'''
		Called to perform class averaging 
		May boot strap the original average, iteratively refines averages, aligns final average to ref 
		'''
		self.init_memory()
		
		if self.verbose:
			print "################## Class",self.data["class_idx"] 
		
	   	if self.options.has_key("bootstrap") and self.options["bootstrap"] == True:
	   		if self.verbose:
	   			print "Bootstrapping the initial class average, #ptcls ", len(self.ptcl_indices)
	   		average,alis = self.__get_bootstrapped_average()
	   	else:
	   		if self.verbose:
	   			print "Generating initial class average using input alignment parameters, #ptcls ", len(self.ptcl_indices)
	   		average = self.__get_init_average_from_ali()
	   		alis = self.data["init_alis"]

	   	self.all_alis.append(alis)
	   	self.averages.append(average)
	   	inclusions = None
	   	
	   	for i in xrange(0,self.options["iter"]):
	   		if i != 0 and False:
	   			alis = self.__align_all(self.averages[-1])
	   			self.all_alis.append(alis)
	   		
	   		threshold = None
	   		sims = None
	   		if self.culling:
	   			sims = self.__cmp_using_ali(self.averages[-1],self.all_alis[-1])
	   			itfrac = float(i+1)/float(self.options["iter"])
	   			threshold = self.__get_cull_threshold(sims,itfrac)
	   		
	   		#print sims,threshold
	   		
	   		average, inclusions = self.__get_average_with_culling(self.all_alis[-1],sims,threshold)
	   		self.all_inclusions.append(inclusions)
	   		
	   		if self.verbose:
	   			vals = inclusions.values()
	   			kept=0
	   			for v in vals:
	   				kept += v
	   			print "Calculating iterative class average ",i+1, "#ptcls",kept
	   		#print inclusions
	   		self.averages.append(average)
	   		
	   	if self.ref != None:
	   		if self.verbose:
	   			print "Aligning final average to reference image"
	   		ref_to_average = self.__align(self.ref,average)
			ref_ali = ref_to_average.get_attr("xform.align2d")
			ref_ali.invert()
			
			ref_alis = {}
			for ptcl_idx,ali in self.all_alis[-1].items():
				ref_alis[ptcl_idx] = ref_ali * ali # this is correct it is easily verified
			self.all_alis.append(ref_alis)
			average = self.__get_average(ref_alis,inclusions)
			self.all_inclusions.append(inclusions) # just for completeness, yes redundant, but yes generically reasonable
			average.set_attr("xform.projection", self.ref.get_attr("xform.projection"))
#			average.set_attr("projection_image",options.ref) # Have to look into this
			average.set_attr("projection_image_idx",self.class_idx)
			self.averages.append(average)
		
		if self.verbose:
			print "****************** Produced",len(self.averages),"averages for class",self.data["class_idx"]
	
	def get_return_data(self):
		'''
		What you call after you call execute
		Returns a dictionary with these keys/values:
		--
		final_average: the final average
		final_alis: a dictionary, keys are particle index, vals are Transforms storing final alignment parameters
		final_inclusions: a dictionary, keys are particle index, vals are True or False to indiciate inclusion in the final average
		--
		If you specified the debug flag you will also get this keys/values:
		---
		averages: the list of averages as they were generated
		all_alis: the list of alignment dictionaries, one corresponding to each average
		all_inclusions: a list of inclusion dictionaries, one for each average, except the first
		'''
		d = {}
		
		if self.options.has_key("debug") and self.options["debug"] == True:
			# if you specify the debug options you get all of the averages back
			d["averages"] = self.averages
			# All of the alignment data
			d["all_alis"] = self.all_alis
			# All of the inclusion information
			d["all_inclusions"] = self.all_inclusions
			
		d["class_idx"] = self.class_idx
		if len(self.averages) != 0:	d["final_average"] = self.averages[-1]
		else: d["final_average"] = None
		
		if len(self.all_alis) != 0:	d["final_alis"] = self.all_alis[-1]
		else: d["final_alis"] = None
		
		if len(self.all_inclusions) != 0: d["final_inclusions"] = self.all_inclusions[-1]
		else: d["final_inclusions"] = None
		
		return d

	def __get_cull_threshold(self,sim_values,itfrac):
		'''
		Get the threshold used to remove particles from the class average
		@param sim_values a dictionary, keys are particle indices (but not used), values are similarity/quality scores
		@param itfrac is the iteration fraction. i.e. if there are 4 iterations, and it's the second iteration, itfrac is 0.5
		@return the threshold one can use to remove particles from the class average
		'''
		
#		def get_cull_threshold(class_cache,weights,options,itfrac):
		qual_scores = sim_values.values()
		
		if (self.options.has_key("keepsig") and self.options["keepsig"] == True):
			a = Util.get_stats_cstyle(qual_scores)
			mean = a["mean"]
			std_dev = a["std_dev"]
			cullthresh = mean + ((self.options["keep"]+1)*(1.0-itfrac)+itfrac*self.options["keep"])*std_dev
		else:
			qual_scores.sort()
			# The ceil reflects a conservative policy. If the user specified keep=0.93
			# and there were 10 particles, then they would all be kept. If floor were
			# used instead of ceil, the last particle would be thrown away (in the
			# class average)
			wf = (1.0-itfrac)+itfrac*self.options["keep"] # weighting function
			idx = int(math.ceil(wf*len(qual_scores)))-1
			#print idx,"was the idx",itfrac*self.options["keep"],wf
			if idx >= len(qual_scores): idx = len(qual_scores)-1
			#raise RuntimeError("Cull thresh idx was too large, something is wrong with the formula")
			cullthresh = qual_scores[idx]
		
		return cullthresh
		
	def __align(self,this,to):
		'''
		Align one image to another, returning the aligned result, using alignment
		parameters stored in self.options
		@param this the image that will be aligned
		@param to the image that this will be aligned to
		@return the aligned image
		'''
		# Align the particle to "to"
		options = self.options
		aligned=this.align(options["align"][0],to,options["align"][1],options["aligncmp"][0],options["aligncmp"][1])
		
		if options.has_key("ralign") and options["ralign"] != None: # potentially employ refine alignment
			refine_parms=options["ralign"][1]
			# this parameters I think west best with the refine aligner, but they're not well tested
			# i.e. I need to do some rigorous testing before I can claim anything
			#refineparms[1]["az"] = ta.get_attr_default("align.az",0)-1
			#refineparms[1]["dx"] = ta.get_attr_default("align.dx",0)-1
			#refineparms[1]["dy"] = ta.get_attr_default("align.dy",0)-1
			#refineparms[1]["mode"] = 0
			#refineparms[1]["stepx"] = 2
			#refineparms[1]["stepy"] = 2
			#refineparms[1]["stepaz"] = 5
			
			refine_parms["xform.align2d"] = aligned.get_attr("xform.align2d")
			aligned = this.align(options["ralign"][0],to,refine_parms,options["raligncmp"][0],options["raligncmp"][1])
			
		return aligned
	
	def __cmp_using_ali(self,average,alis):
		'''
		compare every image to the given average, first transforming using the given alis
		@param average the averaged image, that which you want to compare the raw particles to
		@param alis a dict, key is particle index, values are alignment transforms
		@return a dict, keys are particle index, values are the result of the comparison operation
		'''
		sims = {} # stores similarity score, a measure of quality
			
		for ptcl_idx,ali in alis.items():
			if self.usefilt_images == None:
				image = self.images[ptcl_idx]
			else:
				image = self.usefilt_images[ptcl_idx]
			
			aligned = image.process("math.transform",{"transform":ali})
			
			sims[ptcl_idx] = aligned.cmp(self.options["cmp"][0],average,self.options["cmp"][1])
				
		return sims
	
	def __align_all(self,average):
		'''
		Align all images to the given image
		@param average most probably averaged image, that which you want to align the raw particles to
		@return a dict, keys are particle index, values are the Transforms storing alignment parameters
		'''
		alis = {} # stores Transforms that have alignment parameters in them
		
		for ptcl_idx in self.ptcl_indices:
			if self.usefilt_images == None:
				image = self.images[ptcl_idx]
			else:
				images = self.usefilt_images[ptcl_idx]
				
			aligned = self.__align(average,image)
			ali_parms = aligned.get_attr("xform.align2d")
			ali_parms.invert()
			alis[ptcl_idx] = ali_parms
				
		return alis
	
	def __get_bootstrapped_average(self):
		'''
		Get the bootstrapped average
		You have to call this after calling self.init_memory
		@return the boot strapped average
		'''
		average = None
		np = 0
		alis = {} # holds the transformation matrices that define the alignment 
		for ptcl_idx in self.ptcl_indices:
			if self.usefilt_images == None:
				image = self.images[ptcl_idx]
			else:
				images = self.usefilt_images[ptcl_idx]
					
			if (average == None):
				average = image.copy()
				alis[ptcl_idx] = Transform() # the identity
				np = 1
			else:

				aligned = self.__align(average,image)
				ali_parms = aligned.get_attr("xform.align2d")
				ali_parms.invert()
				alis[ptcl_idx] = ali_parms
				aligned = image.process("math.transform",{"transform":(ali_parms)})
								
				np += 1
				average.add(aligned) # now add the image
	
		if np != 0:
			if self.norm != None: average.process_inplace(self.norm[0],self.norm[1])
			average.process_inplace("xform.centeracf")
			average.process_inplace("mask.sharp",{"outer_radius":average.get_xsize()/2})
			average.set_attr("ptcl_repr",np)
		return average,alis

	def __get_init_average_from_ali(self):
		'''
		Get the initial average using initial alignment parameters
		You have to call this after calling self.init_memory
		@return the initial average
		'''
		averager_parms =  self.options["averager"]
		averager=Averagers.get(averager_parms[0], averager_parms[1])
		
		weightsum = 0 # used to normalize the average
		np = 0 # number of particles in the average
		
		for ptcl_idx,ali in self.data["init_alis"].items():
			image = self.images[ptcl_idx]			
			rslt = image.process("math.transform",{"transform":ali})
					
			np += 1
			weight = self.data["weights"][ptcl_idx]
			weightsum += weight
			rslt.mult(weight)
			
			# Add the image to the averager
			averager.add_image(rslt)
		
		if np == 0 or weightsum == 0:
#			if options.verbose:
#				print "Class",cl,"...no particles"
			# write blank image? Write meta data?
			return None
			
		average = averager.finish()
		average.mult(float(np)) # Undo the division of np by the averager - this was incorrect because the particles were weighted.
		average.mult(1.0/weightsum) # Do the correct division
		if self.norm != None: average.process_inplace(self.norm[0],self.norm[1])
		average.process_inplace("xform.centeracf")
		average.process_inplace("mask.sharp",{"outer_radius":average.get_xsize()/2})
		average.set_attr("ptcl_repr",np)
		
		return average
		
	def __get_average_with_culling(self,alis,sims,cullthresh):
		'''
		get an average using the given alignment parameters, similarity scores, and cull threshold
		cull threshold can be None, in which case no attention is payed to the sims argument
		@param alis a dict, keys are particle indices, values are Transforms storing alignment parameters
		@param sims a dict, keys are parictle indices, values are numbers. This argument is ignored if cullthresh is None
		@param cullthresh, a threshold value used to keep particles out of the class average as based on sims
		@return the average
		'''
		averager_parms =  self.options["averager"]
		averager=Averagers.get(averager_parms[0], averager_parms[1])
		
		inclusion = {} # stores a flag, True or False, indicating whether the particle was included in the class average
		np = 0 # number of particles
		for ptcl_idx,ali in alis.items():
			if ( cullthresh != None ):
				if ( sims[ptcl_idx] > cullthresh ) :
					inclusion[ptcl_idx] = False
					continue
				
			inclusion[ptcl_idx] = True
			image = self.images[ptcl_idx]			
			rslt = image.process("math.transform",{"transform":ali})
		
			np += 1
			averager.add_image(rslt)
	
		if np == 0:
			average = None
		else:
			average = averager.finish()
			if self.norm != None: average.process_inplace(self.norm[0],self.norm[1])
			#should this be centeracf?
			#average.process_inplace("xform.centerofmass", {"int_shift_only":1})
			average.process_inplace("xform.centeracf")
			average.process_inplace("mask.sharp",{"outer_radius":average.get_xsize()/2})
			average.set_attr("ptcl_repr",np)
		return average,inclusion
	
	def __get_average(self,alis,inclusions):
		'''
		get an average using the given alignment parameters and inclusion flags
		inclusions may be None, in which case every particle is included
		@param alis a dict, keys are particle indices, values are Transforms storing alignment parameters
		@param inclusions a dict, keys are parictle indices, values are True or False to indicate inclusion.
		@return the average
		'''
		
		averager_parms =  self.options["averager"]
		averager=Averagers.get(averager_parms[0], averager_parms[1])
		
		np = 0 # number of particles
		for ptcl_idx,ali in alis.items():
			if inclusions != None and inclusions[ptcl_idx] == False: continue
			
			image = self.images[ptcl_idx]			
			rslt = image.process("math.transform",{"transform":ali})
		
			np += 1
			averager.add_image(rslt)
			
		if np == 0:
			average = None
		else:
			average = averager.finish()
			if self.norm != None: average.process_inplace(self.norm[0],self.norm[1])
			#should this be centeracf?
			#average.process_inplace("xform.centerofmass")
			average.process_inplace("xform.centeracf")
			average.process_inplace("mask.sharp",{"outer_radius":average.get_xsize()/2})
			average.set_attr("ptcl_repr",np)
		return average
		
		return average,inclusion
		
class EMClassAveTaskDC(EMClassAveTask):
	'''
	This is the class average task used by the distributing computing parallelism framework
	'''
	
	def __init__(self,command="e2classaverage",data=None,options=None):
		EMClassAveTask.__init__(self,command,data,options)
	
	def init_memory(self):
		'''
		Reads images of the class into memory, normalize them if the options dictate it
		Can't be private because of class inheritance issues
		'''
		cache_name=self.data["input"][1]
		image_cache=db_open_dict(cache_name)
		self.ptcl_indices = self.data["input"][2]
		
		self.images = {}
		if self.data.has_key("usefilt") and self.data["usefilt"] != None:
			self.usefilt_images = {}
			usefilt_cache_name=self.data["usefilt"][1]
			usefilt_image_cache=db_open_dict(usefilt_cache_name)
		else:
			self.usefilt_images = None
		
		norm = None
		if self.options.has_key("normproc") and self.options["normproc"] != None:
			norm = self.options["normproc"] # a list of length two
		self.norm = norm
			
		for ptcl_idx in self.ptcl_indices:
			image = image_cache[ptcl_idx]
			if norm != None:
				image.process_inplace(norm[0],norm[1])
			self.images[ptcl_idx] = image
			
			if self.usefilt_images != None:
				usefilt_image = usefilt_image_cache[ptcl_idx]
				if norm != None:
					usefilt_image.process_inplace(norm[0],norm[1])
				self.usefilt_images[ptcl_idx] = usefilt_image
		# set a flag that can be used to decide if quality scores need to be
		# evaluated and used for culling the worst, or least similar, particles
		# from the class
		self.culling = False
		if (self.options["iter"] > 0):
			self.culling = True
			if (self.options.has_key("keepsig") and self.options["keepsig"] == False) and self.options["keep"] == 1.0 : self.culling=False
		
		# the reference is used to perform a final alignment, if specified
		self.ref = None
		if self.data.has_key('ref') and self.data['ref'] != None:
			cache_name=self.data["ref"][1]
			ref_cache=db_open_dict(cache_name)
			self.ref = ref_cache[self.data["ref"][2]]
		
		self.verbose = 0
		if self.options.has_key("verbose") and self.options["verbose"] != None:
			self.verbose = self.options["verbose"]
			
	

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
	parser.add_option("--parallel", default=None, help="parallelism argument")

	
	(options, args) = parser.parse_args()
	
	if len(args)<3 : parser.error("Input, classification matix, and output files required")
	
	if (options.check): options.verbose = True # turn verbose on if the user is only checking...
		
	if (options.nofilecheck == False):
		options.classifyfile=args[1]
		options.datafile = args[0]
		options.outfile = args[2]
	
	error = check(options,True,args)
	
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
	
	class_gen = EMGenClassAverages(options,args,logger)
	class_gen.execute()
	
	E2end(logger)
#	
#	# just remove the file - if the user didn't specify force then the error should have been found in the check function
#	if os.path.exists(options.outfile):
#		if options.force:
#			remove_file(options.outfile)
#	
#	(num_classes, num_part ) = gimme_image_dimensions2D(args[1]);
#
#	
#	if (options.verbose):
#		print "Classifications per particle %d, particles %d" %(num_classes, num_part)
#	
#	# classes contains the classifications - row is particle number, column data contains class numbers (could be greater than 1)
#	classes = EMData()
#	classes.read_image(args[1], 0)
#	class_max = int(classes.get_attr("maximum"))
#	class_min = int(classes.get_attr("minimum"))
#	
#	# weights contains the weighting of the classification scheme stored in the EMData object "classes" - above
#	# dx contains the x translation of the alignment
#	# dy contains the y translation of the alignment
#	# da contains is the azimuthal rotation of the alignment
#	# dflip contains is the mirror alignment
#	# row is particle number, column data contains weights - rows should add to 1, but this is not checked.
#	weights, dx, dy, da, dflip = EMData(),EMData(),EMData(),EMData(),EMData()
#	if options.bootstrap:
#		weights.set_size(classes.get_xsize(),classes.get_ysize())
#		weights.to_zero()
#		dx.set_size(classes.get_xsize(),classes.get_ysize())
#		dx.to_zero()
#		dy.set_size(classes.get_xsize(),classes.get_ysize())
#		dy.to_zero()
#		da.set_size(classes.get_xsize(),classes.get_ysize())
#		da.to_zero()
#		dflip.set_size(classes.get_xsize(),classes.get_ysize())
#		dflip.to_zero()
#	else:
#		if EMUtil.get_image_count(args[1]) != 6:
#			print "error, the classification matrix is the wrong size, it needs to contain one image for the classes, weights, dx, dy, da, and dflip. You can bypass this requirement if you supply the bootstrap argument"
#			sys.exit(1)
#			
#		else:
#			weights.read_image(args[1], 1)
#			dx.read_image(args[1],2)
#			dy.read_image(args[1],3)
#			da.read_image(args[1],4)
#			dflip.read_image(args[1],5)
#	
#	# empty space for storing x-flipping flags (0s or 1s, if a 1 is stored it will be used at the necessary point to flip prior to adding to the average)
#	# sometimes this will not be used at all (it depends on whether or not the aligners that the user has specified do flipping and set flip flags)
#	# dflip = EMData(da.get_xsize(),da.get_ysize())
#	# dflip.to_zero()
#	
#	options.norm = parsemodopt(options.normproc)
#	
#	if (options.iter > 0 or options.bootstrap):
#		set_aligner_params_in_options(options)
#	
#	if options.idxcache:
#		try:
#			name = "bdb:"+options.dbpath+"#class_indices"
#			if options.even: name += "_even"
#			elif options.odd: name += "_odd"
#			db_name = numbered_bdb(name)
#			class_db = db_open_dict(db_name)
#		except:
#			print "error with db terminology: can't call numbered_bdb with this argument:", "bdb:"+options.dbdir+"#class_indices"
#			sys.exit(1)
#	# do one class at a time
#	for cl in xrange(class_min,class_max+1):
#		print cl
#		E2progress(logger,float(cl-class_min)/(class_max+1-class_min))
#		class_cache = ClassCache(args,options,num_part,num_classes,classes)
#		class_cache.update_cache(cl)
#		
#		ndata = []
#			
#		if options.idxcache:
#			class_indices = class_cache.class_indices
#
#		if (options.iter > 0 or options.verbose):
#			ccache = class_cache.ccache
#		
#		# this should work, if it doesn't it should have been caught by check function
#		averager_parms=parsemodopt(options.averager)
#		average = None
#		image_cache = class_cache.image_cache
#		filtered_image_cache = class_cache.filtered_image_cache
#		if ( not options.bootstrap ):
#			if options.verbose: print "generating the original class average using alignment parameters in the classification matrix"
#			# generate the first class average by applying the transformations, adding and finally normalizing...		
#			average,np = get_basic_average(averager_parms,class_cache,da,dflip,dx,dy,weights,options)
#		else:
#			# generate a bootstrapped initial average. Do this 'inductively' by aligning the 2nd image to the first, then adding. Continue until done...
#			if options.verbose: print "bootstrapping the original class average"
#			average,np = get_boot_strapped_average(class_cache,options,args)
#		
#		if np == 0:
#			if options.verbose:
#				print "Class",cl,"...no particles"
#			# FIXME
#			continue
#		
#		if options.idxcache:
#			class_db[str(cl)] = class_indices
#				
#		
#		
#		if options.verbose: print "done original average"
#	
#		if (options.iter > 0):
#			options.cull = True
#			if ( options.keepsig == False and options.keep == 1.0 ) : options.cull = False
#		
#		for it in range(0,options.iter):
#			if average == None:
#				print "error, the average was None"
#				break
#			itfrac=(it+1)/float(options.iter)
#			# do alignment
#			align_to_average_and_store_metadata(class_cache,options,average,dx,dy,da,dflip,weights)
##			
#			# get the culling threshold
#			cullthresh = None
#			if options.cull:
#				cullthresh = get_cull_threshold(class_cache,weights,options,itfrac)
#							
#			#finally average with culling
#			average,np = get_basic_average_with_cull(averager_parms,class_cache,options,weights,da,dx,dy,dflip,ndata,cullthresh)			
#			if options.verbose: print "done iteration",it
#
#		if average == None: 
#			continue
#		# Align to the reference and extract euler data if it was specified
#		if ( options.ref):
#			e = EMData()
#			if options.iter > 0:
#				e.read_image(options.ref, cl)
#				
#				average_to_ref = align(e,average,options)
#				fine_transform = average_to_ref.get_attr("xform.align2d")
#				fine_transform.invert()
#				
#				average,np = get_basic_average_with_tweak(averager_parms,class_cache,options,weights,da,dx,dy,dflip,fine_transform)
#			else:
#				e.read_image(options.ref, cl, READ_HEADER_ONLY)
#
#			average.set_attr("xform.projection", e.get_attr("xform.projection"))
#			average.set_attr("projection_image",options.ref)
#			average.set_attr("projection_image_idx",cl)
#			
#
#			if options.verbose:
#				edata = []
#				t = e.get_attr("xform.projection")
#				rot = t.get_rotation("eman")
#				edata.append(rot["alt"])
#				edata.append(rot["az"])
#				edata.append(rot["phi"])
#		
#		# now write to disk
#		average.set_attr("ptcl_repr",np)
#		if np == 0:
#			a = average.copy() # get all the attributes
#			a.to_zero()
#			a.write_image(args[2],-1)
#		else:
#			average.write_image(args[2],-1)
#			
#		if options.verbose:
#			sys.stdout.write( "Class %d: particles..%d" %(cl,len(ccache)) )
#			for t in range(0,options.iter):
#				sys.stdout.write("..%d" %ndata[t] )
#			if ( options.ref  ):
#				sys.stdout.write(" : Eulers..")
#				for t in edata:
#					sys.stdout.write(" %f" %t)
#			sys.stdout.write("\n")
#		
#	if (options.resultmx != None ):
#		if os.path.exists(options.resultmx):
#			remove_file(options.resultmx) #oooh not sure about this!
#		
#		# note the order is important!
#		classes.write_image(options.resultmx,-1)
#		weights.write_image(options.resultmx,-1)
#		dx.write_image(options.resultmx,-1)
#		dy.write_image(options.resultmx,-1)
#		da.write_image(options.resultmx,-1)
#		dflip.write_image(options.resultmx,-1)
#		
#	if options.idxcache:
#			name = "bdb:"+options.dbpath+"#class_indices"
#			if options.even: name += "_even"
#			elif options.odd: name += "_odd"
#			db_name = numbered_bdb(name)
#			#db_close_dict(db_name)
#
#	E2end(logger)
#
#def get_cull_threshold(class_cache,weights,options,itfrac):
#	qual_scores = []
#	for d in class_cache.ccache:
#		p = d[0]
#		c = d[1]
#		qual_scores.append(weights.get(c,p))
#	
#	if ( options.keepsig ):
#		a = Util.get_stats_cstyle(qual_scores)
#		mean = a["mean"]
#		std_dev = a["std_dev"]
#		cullthresh = mean + ((options.keep+1)*(1.0-itfrac)+itfrac*options.keep)*std_dev
#	else:
#		b = deepcopy(qual_scores)
#		b.sort()
#		# The ceil reflects a conservative policy. If the user specified keep=0.93
#		# and there were 10 particles, then they would all be kept. If floor were
#		# used instead of ceil, the last particle would be thrown away (in the
#		# class average)
#		wf = (1.0-itfrac)+itfrac*options.keep # weighting function
#		idx = int(math.ceil(wf*len(b)))-1
#		#print idx,"was the idx",itfrac*self.options["keep"],wf
#		if idx >= len(qual_scores): idx = len(b)-1
#			#raise RuntimeError("Cull thresh idx was too large, something is wrong with the formula")
##		cullthresh = b[idx]
##		idx = int(math.ceil(((1.0-itfrac)+itfrac*options.keep)*(len(b)-1)))
##		if idx >= len(b):
##			raise RuntimeError("Cull thresh idx was too large, something is wrong with the formula")
#		cullthresh = b[idx]
#	
#	return cullthresh
#	
#def align_to_average_and_store_metadata(class_cache,options,average,dx,dy,da,dflip,weights):
#	if options.usefilt: image_cache = class_cache.filtered_image_cache
#	else: image_cache = class_cache.image_cache
#	for p,image in image_cache.items():
#		
#		c = class_cache.idx_cache[p]
#		ta = align(average,image,options)
#		t = ta.get_attr("xform.align2d")
#		t.invert() # because we aligned the average to the image, and not the other way around
#		params = t.get_params("2d")
#		
#		# store the refined translational and rotational values
#		dx.set(c,p, params["tx"])
#		dy.set(c,p, params["ty"])
#		da.set(c,p, params["alpha"])
#		dflip.set(c,p, params["mirror"])
#		
#		# store the quality score on top of the weights, seeing as the starting weights are no longer required
#		if (options.cull): # but only if we need to
#			weights.set(c,p, ta.cmp(options.cmp[0],image,options.cmp[1]))
#
#def get_basic_average_with_tweak(averager_parms,class_cache,options,weights,da,dx,dy,dflip,fine_transform):
#	averager=Averagers.get(averager_parms[0], averager_parms[1]) # need an empty averager
#	np = 0 # number of particles in the average
#	for p,image in class_cache.image_cache.items():
#		c = class_cache.idx_cache[p]
#		if (options.cull ):
#			if ( weights.get(c,p) == 0 ) : continue
#		
#		t = Transform({"type":"2d","alpha":da.get(c,p)})
#		t.set_trans(dx.get(c,p),dy.get(c,p))
#		if dflip.get(c,p) != 0: t.set_mirror(True)
#		ct = fine_transform*t
#
#		rslt = image.process("math.transform",{"transform":ct})
#		
#		np += 1
#		averager.add_image(rslt)
#		
#		# Don't forget to store the alignment parameters
#		params = ct.get_params("2d")
#		dx.set(c,p, params["tx"])
#		dy.set(c,p, params["ty"])
#		da.set(c,p, params["alpha"])
#		dflip.set(c,p, params["mirror"])
#	
#	if np == 0:
#		if (options.verbose):
#			print "Class",cl,"...no particles on iteration",it
#		average = None
#	else:	
#		average = averager.finish()
#		#				average.transform(fine_transform)
#		
#		if str(options.normproc) != "None": average.process_inplace(options.norm[0],options.norm[1])
#		#average.process_inplace("xform.centerofmass") this shouldn't be necessary if we aligned to the projection
#		average.process_inplace("mask.sharp",{"outer_radius":average.get_xsize()/2})
#		average.set_attr("ptcl_repr",np)
#		
#	return average,np
#
#def get_basic_average_with_cull(averager_parms,class_cache,options,weights,da,dx,dy,dflip,ndata,cullthresh):
#	averager=Averagers.get(averager_parms[0], averager_parms[1]) # need an empty averager
#	#avg = None
#	np = 0 # number of particles in the average
#	for p,image in class_cache.image_cache.items():
#		c = class_cache.idx_cache[p]
#		if (options.cull ):
#			# if we're culling then find the bad particles here...
#			# Also store whether or not culling occured by placing 1s and 0s in the weights object
#			# This is in case users ever want to know which particles were excluded. This information
#			# is written to disk if the resultmx option is specfied
#			# FIXME setting the weights matrix to 1 and 0 should probably only occur if we're at the
#			# last iteration
#			if ( weights.get(c,p) > cullthresh ) :
#				weights.set(c,p,0)
#				continue
#			else: weights.set(c,p,1)
#		else: weights.set(c,p,1)
#
#		t = Transform({"type":"2d","alpha":da.get(c,p)})
#		t.set_trans(dx.get(c,p),dy.get(c,p))
#		if dflip.get(c,p) != 0: t.set_mirror(True)
#
#		rslt = image.process("math.transform",{"transform":t})
#	
#		np += 1
#		averager.add_image(rslt)
#		
#	if options.verbose:
#		ndata.append(np)
#
#	if np == 0:
##		if (options.verbose):
##			print "Class",cl,"...no particles on iteration",it
#		# FIXME
#		average = None
#	else:
#
#		average = averager.finish()
#		#should this be centeracf?
#		if str(options.normproc) != "None": average.process_inplace(options.norm[0],options.norm[1])
#		#average.process_inplace("xform.centerofmass",{"int_shift_only":1})
#		average.process_inplace("xform.centeracf")
#		average.process_inplace("mask.sharp",{"outer_radius":average.get_xsize()/2})
#		average.set_attr("ptcl_repr",np)
#	
#	return average,np
#			
#
#def get_boot_strapped_average(class_cache,options,args):
#	average = None
#	np = 0
#	
#	for p,image in class_cache.image_cache.items():
#		c = class_cache.idx_cache[p]
#
#				
#		if (average == None):
#			average = class_cache.image_cache[p].copy()
#			#average.read_image(args[0],p)
#
#			if str(options.normproc) != "None": average.process_inplace(options.norm[0],options.norm[1])
#			#average.process_inplace("mask.sharp",{"outer_radius":average.get_xsize()/2})
#			#average.write_image("ave_.hdf",-1)
#			np = 1
#		else:
#			# there are a lot of ways to do and there is probably
#			# an improvement on what happens here, but it suffices
#			if options.usefilt:
#				image = class_cache.filtered_image_cache[p]
#			else:
#				image = class_cache.image_cache[p]
#				
#			ta = align(average,image,options)
#			t = ta.get_attr("xform.align2d")
#			t.invert()
#			ta = image.process("math.transform",{"transform":(t)})
#
#			#ta.process_inplace("mask.sharp",{"outer_radius":ta.get_xsize()/2})				
#			np += 1
#			average.add(ta) # now add the image
#	
#	if np != 0:
#		if str(options.normproc) != "None": average.process_inplace(options.norm[0],options.norm[1])
#		average.process_inplace("xform.centeracf")
#		average.process_inplace("mask.sharp",{"outer_radius":average.get_xsize()/2})
#		average.set_attr("ptcl_repr",np)
#	
#	return average,np
#
#def get_basic_average(averager_parms,class_cache,da,dflip,dx,dy,weights,options):
#	averager=Averagers.get(averager_parms[0], averager_parms[1])
#	# do the initial average, based on the program inputs
#	weightsum = 0 # used to normalize the average
#	np = 0 # number of particles in the average
#	
#	for p,image in class_cache.image_cache.items():
#		c = class_cache.idx_cache[p]
#
#		t = Transform({"type":"2d","alpha":da.get(c,p),"mirror":int(dflip(c,p))})
#		t.set_trans(dx.get(c,p),dy.get(c,p))
#		
#		rslt = image.process("math.transform",{"transform":t})
#				
#		np += 1
#		weight = weights(c,p)
#		weightsum += weight
#		rslt.mult(weight)
#		
#		# Add the image to the averager
#		averager.add_image(rslt)
#	
#	if np == 0 or weightsum == 0:
##		if options.verbose:
##			print "Class",cl,"...no particles"
#		# FIXME
#		
#		# write blank image? Write meta data?
#		return None,0
#		
#	average = averager.finish()
#	average.mult(float(np)) # Undo the division of np by the averager - this was incorrect because the particles were weighted.
#	average.mult(1.0/weightsum) # Do the correct division
#	if str(options.normproc) != "None": average.process_inplace(options.norm[0],options.norm[1])
#	average.process_inplace("xform.centeracf")
#	average.process_inplace("mask.sharp",{"outer_radius":average.get_xsize()/2})
#	average.set_attr("ptcl_repr",np)
#	
#	return average,np
#
#class ClassCache:
#	
#	'''
#	Stores useful information during each iteration
#	'''
#	def __init__(self,args,options,num_part,num_classes,class_data):
#		self.args = args
#		self.options = options
#		self.num_part = num_part
#		self.num_classes = num_classes
#		self.class_data = class_data
#		self.image_cache = {}
#		self.filtered_image_cache= {}
#		self.idx_cache = {}
#		self.class_indices = []
#		self.ccache = []
#		
#	def update_cache(self,class_number):
#		
#		for p in xrange(0,self.num_part):
#			if self.options.odd and p % 2 == 0: continue # enforce even/odd constraints
#			if self.options.even and p % 2 == 1: continue # enforce even/odd constraints
#			for c in xrange(0,self.num_classes):
#				if self.class_data.get(c,p) == class_number:
#					
#					if self.options.idxcache: self.class_indices.append(p)
#					# cache the hit if necessary - this is if there is more than one iteration, and/or the user has specifed verbose. In the
#					# latter case the class cache is used to print information
#					if (self.options.iter > 0 or self.options.verbose): self.ccache.append((p,c))
#					
#					self.idx_cache[p] = c
#					
#					image = EMData()
#					image.read_image(self.args[0],p)
#					if str(self.options.normproc) != "None": image.process_inplace(self.options.norm[0],self.options.norm[1])
#					self.image_cache[p] = image
#					
#					if self.options.usefilt:
#						filt_image = EMData()
#						filt_image.read_image(options.usefilt,p)
#						if str(self.options.normproc) != "None": filt_image.process_inplace(self.options.norm[0],self.options.norm[1])
#						self.filtered_image_cache[p] = filt_image
#
#						
#
#def set_aligner_params_in_options(options):
#	'''
#	Call this before calling align
#	'''
#	options.align=parsemodopt(options.align)
#	options.alicmp=parsemodopt(options.aligncmp)
#	# note the parsing of the options.ralign parameters is left for later
#	options.alircmp=parsemodopt(options.raligncmp)
#	options.cmp=parsemodopt(options.cmp)
#
#def align(this,to,options):
#	# Align the particle to the average
#	ta=this.align(options.align[0],to,options.align[1],options.alicmp[0],options.alicmp[1])
#	
#	if ( options.ralign != None ): # potentially employ refine alignment
#		
#		refineparms=parsemodopt(options.ralign)
#		# this parameters I think west best with the refine aligner, but they're not well tested
#		# i.e. I need to do some rigorous testing before I can claim anything
#		#refineparms[1]["az"] = ta.get_attr_default("align.az",0)-1
#		#refineparms[1]["dx"] = ta.get_attr_default("align.dx",0)-1
#		#refineparms[1]["dy"] = ta.get_attr_default("align.dy",0)-1
#		#refineparms[1]["mode"] = 0
#		#refineparms[1]["stepx"] = 2
#		#refineparms[1]["stepy"] = 2
#		#refineparms[1]["stepaz"] = 5
#		
#		refineparms[1]["xform.align2d"] = ta.get_attr("xform.align2d")
#		ta = this.align(refineparms[0],to,refineparms[1],options.alircmp[0],options.alircmp[1])
#		
#	return ta
#
#		
def check(options,verbose=False,args=[]):
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
				
			
		if options.ref != None and not file_exists(options.ref):
			print "Error: the file expected to contain the reference images (%s) does not exist" %(options.ref)
			error = True
		elif options.ref and (os.path.exists(options.datafile) or db_check_dict(options.datafile)):
			(xsize, ysize ) = gimme_image_dimensions2D(options.datafile);
			(pxsize, pysize ) = gimme_image_dimensions2D(options.ref);
			if ( xsize != pxsize ):
				error = True
				if (verbose):
					print "Error - the dimensions of the reference and particle images do not match"
		elif len(args) > 1 and options.ref:
				# classes contains the classifications - row is particle number, column data contains class numbers (could be greater than 1)
			classes = EMData()
			classes.read_image(args[1], 0,True)
			class_max = int(classes["maximum"])
			num_ref= EMUtil.get_image_count(options.ref)
			if ( class_max > num_ref ):
				print "Error, the classification matrix refers to a class number (%d) that is beyond the number of images (%d) in the reference image (%s)." %(class_max,num_ref,options.ref)
			
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
