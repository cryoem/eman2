#!/usr/bin/env python
#
# Author:  David Woolford 06/02/2008 (woolford@bcm.edu)
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

from EMAN2 import *
from time import time
from math import ceil
from copy import copy


class EMProjectDB:
	"""
	It's implemented as a singleton
	"""
	
	class __impl:
		""" Implementation of the singleton interface """

		def __init__(self):
			DB = EMAN2db.EMAN2DB.open_db(".")
			DB.open_dict("boxer_cache")
			self.project_db = DB.boxer_cache
			self.memory_db = {}
			#self.project_db = shelve.open('.eman2project_db','c',-1,True)

		#def __del__(self):
			#print "closing project_db"
			#self.project_db.close()
	# storage for the instance reference
		#def getAutoBoxers(self):
			#autoboxers = []
			#for i in self.project_db.items():
				#if i[0][0:10] == "autoboxer_":
					#ab = self.project_db[i[0]]
					#swarmAutoBoxer = SwarmAutoBoxer(None)
					#swarmAutoBoxer.become(ab)
					#autoboxers.append(
			
		def close(self):
			#self.project_db.close()
			DB = EMAN2db.EMAN2DB.open_db(".")
			DB.close_dict("boxer_cache")
			
		def set_key_entry(self,key,entry):
			self.project_db[key]= entry
			
		def get_key_entry(self,key):
			return self.project_db[key]
			
		def set_key_entry_in_memory(self,key,entry):
			self.memory_db[key] = entry
		
		def get_key_entry_in_memory(self,key):
			try:
				
				return self.memory_db[key]
			except:
				return self.project_db[key]
			#self.project_db.sync()
			
	__instance = None
	
	def __init__(self):
		""" Create singleton instance """
		# Check whether we already have an instance
		if EMProjectDB.__instance is None:
			# Create and remember instance
			EMProjectDB.__instance = EMProjectDB.__impl()
	
	def __getattr__(self, attr):
		""" Delegate access to implementation """
		return getattr(self.__instance.project_db, attr)

	def __setattr__(self, attr, value):
		""" Delegate access to implementation """
		return setattr(self.__instance, attr, value)
	
	def set_key_entry(self,key,entry):
		return self.__instance.set_key_entry(key,entry)
	
	def get_key_entry(self,key):
		return self.__instance.get_key_entry(key)
	
	def get_key_entry_in_memory(self,key):
		return self.__instance.get_key_entry_in_memory(key)
	
	def set_key_entry_in_memory(self,key,entry):
		return self.__instance.set_key_entry_in_memory(key,entry)
	
#	def cache_key_entry_from_memory_to_disk(self,key,entry):
#		return self.__instance.set_key_entry_in_memory(key,entry)
	
class Box:
	CENTERACF = "centeracf"
	#CENTERALIGNINT = "cenlignint"
	CENTEROFMASS = "centerofmass"
	CENTERPROPAGATE = "centerpropagate"
	CENTERMETHODS = [CENTERACF,CENTERPROPAGATE,CENTEROFMASS]
	def become(self,trimbox):
		'''
		This is like a copy constructor
		'''
		self.xcorner = trimbox.xcorner			# the xcorner - bottom left
		self.ycorner = trimbox.ycorner			# the ycorner - bottom left
		self.xsize = trimbox.xsize				# the xsize of the box
		self.ysize = trimbox.ysize				# the ysize of the box
		self.isref = trimbox.isref				# a flag that can be used to tell if the box is being used as a reference
		self.changed = trimbox.changed			# a flag signalling the box has changed and display needs updatin
		self.isanchor = trimbox.isanchor		# a flag signalling the box has changed and display needs updatin
		self.TS = trimbox.TS
		self.image_name = trimbox.image_name
		self.ismanual = trimbox.ismanual
		
		self.moved = trimbox.moved
		self.origxcorner = trimbox.origxcorner
		self.origycorner = trimbox.origycorner
		
	def __init__(self,xcorner=-1,ycorner=-1,xsize=-1,ysize=-1,isref=0,correlationscore=0,image_name=None):
		self.xcorner = xcorner			# the xcorner - bottom left
		self.ycorner = ycorner			# the ycorner - bottom left
		self.xsize = xsize				# the xsize of the box
		self.ysize = ysize				# the ysize of the box
		self.isref = isref				# a flag that can be used to tell if the box is being used as a reference
		self.correlationscore = correlationscore	# the correlation score
		self.ismanual = False			# a flag to store whether or this box was manually added by the user and was not a reference. It's just a plain box
		
		self.optprofile = None			# a correlation worst-case profile, used for selective auto boxing
		self.changed = False			# a flag signalling the box has changed and display needs updating
		self.corx = -1			# stores the x coordinate of the correlation peak
		self.cory = -1			# stores the y coordinate of the correlation peak
		self.shape = None		# stores the shape used by the image2d widget
		self.image = None 		# stores the image itself, an emdata object
		self.r = 0.4			# RGB red
		self.g = 0.9			# RGB green
		self.b = 0.4			# RGB blue
		self.rorig = 0.4			# RGB red
		self.gorig = 0.9			# RGB green
		self.borig = 0.4			# RGB blue
		self.footprint = None	# stores the image footprint as an emdata object
		self.group = None		# stores a group, typically an int
		self.footprintshrink = 1
		self.isanchor = True		# A flag used by AutoBoxer routines that - if set to true the box will not be included in the generation the template) - This is specific to the SwarmPS autoboxer
		self.TS = None
		
		self.isdummy = False # this can be used to avoid parameters updates - i.e. when the user interactively changes parameters forcefully
		self.image_name = image_name
		
		self.moved = False
		self.origxcorner = -1
		self.origycorner = -1
		
	
	def set_image_name(self,image_name):
		self.image_name = image_name
		
	def get_image_name(self):
		return self.image_name
	
	def move(self,dx,dy):
		if self.moved == False:
			self.origxcorner = self.xcorner
			self.origycorner = self.ycorner
			self.moved = True
		
		self.xcorner += dx
		self.ycorner += dy
		
		if not self.ismanual:
			movedboxes = get_idd_key_entry_in_memory(self.image_name,"moved_boxes") # this may potentially be None
			
			if movedboxes == None:
				movedboxes = []
			
			found = False
			for data in movedboxes:
				if data[0] == self.origxcorner and data[1] == self.origycorner:
					data[2] = self.xcorner
					data[3] = self.ycorner
					found = True
					break
					
			if not found:
				movedboxes.append([self.origxcorner,self.origycorner,self.xcorner,self.xcorner])
			
			
			set_idd_key_entry_in_memory(self.image_name,"moved_boxes",movedboxes)
		
		self.changed = False
		self.update_box_image()

	def update_position_from_data(self,movedboxes):
		
		# 0.0064 = 0.08*0.08, then divided by two to make in terms of the radius
		# ... so the proximity limit is 8% of the radius
		sq_size_limit = (self.xsize**2 + self.ysize**2)*0.0032

		for data in movedboxes:
			sq_dif = (data[0] - self.xcorner)**2 + (data[1] - self.ycorner)**2
			if sq_dif < sq_size_limit:
				self.origxcorner = data[0]
				self.origycorner = data[1]
				self.xcorner = data[2]
				self.ycorner = data[3]
				self.changed = True
				self.moved = True
				self.update_box_image()
				return 1
		
		return 0
	
	def update_box_image(self,norm=True):
		bic = BigImageCache()
		image = bic.get_image(self.image_name)
		#print "getting region",self.xcorner,self.ycorner,self.xsize,self.ysize
		if self.xcorner + self.xsize > image.get_xsize(): self.xcorner = image.get_xsize()-self.xsize
		if self.ycorner + self.ysize > image.get_ysize(): self.ycorner = image.get_ysize()-self.ysize
		
		self.image = image.get_clip(Region(self.xcorner,self.ycorner,self.xsize,self.ysize))
		if norm:
			self.image.process_inplace("normalize.edgemean")
		
		# make sure there are no out of date footprints hanging around
		self.footprint = None

	def change_box_size(self,box_size):
		'''
			Changes the box_size if it really is different to that which is currently stored
			Release references to self.image and self.footprint so that they are generated next time
			they are asked for	
		'''
		if box_size != self.xsize or box_size != self.ysize:
			self.xcorner -= (box_size-self.xsize)/2
			self.xsize = box_size
			self.ycorner -= (box_size-self.ysize)/2
			self.ysize = box_size
			
			self.image = None
			self.footprint = None
			
	def get_box_image(self,norm=True,force=False):
		bic = BigImageCache()
		image = bic.get_image(self.image_name)
		
		if self.image == None or force:
			if image == None:
				print 'error, need to specify the image argument when first calling get_box_image'
			self.update_box_image(norm)
		return self.image
	
	def get_small_box_image(self,autoboxer):
		'''
		gets a shrunken version of the box by asking the database if the shrunken (entire) image
		exists and then clipping out from it
		'''
		image = self.get_small_image(autoboxer)
		if image == None:
			return None
		else:
			shrink = autoboxer.get_best_shrink()
			return image.get_clip(Region(int(self.xcorner/shrink),int(self.ycorner/shrink),int(self.xsize/shrink),int(self.ysize/shrink)))
		
	def get_small_image(self,autoboxer):
		
		#cache = CoarsenedFlattenedImageCache()
		cache = SubsamplerCache
		return  cache.get_image(self.image_name, autoboxer.get_params_mediator())
	
	
	def get_flcf_image(self,autoboxer):
		
		cache = FLCFImageCache()
		return cache.get_image(self.image_name,autoboxer.get_params_mediator())

	
	def get_foot_print(self,shrink=1):
		if self.footprint == None or shrink != self.footprintshrink:
			self.footprintshrink = shrink
			if self.image == None:
				print "error, you can not make a footprint if there is no image"
				exit(1)
			if shrink == 1:
				self.footprint = self.image.make_footprint()
			else :
				self.footprint = self.image.process("math.meanshrink",{"n":shrink}).make_footprint()
				
		return self.footprint
			
	def center(self,method,extrasomething,low_res=False):
		'''
		Ask the box to center itself using one of the available methods (as stored in Box.CENTERMETHODS
		extrasomething has to be an AutoBoxer if using CENTEROFMASS or CENTERACF (it's asked for get_best_shrink)
		extrasomething has to be the template image and have the same dimensions as the results of get_box_image if using 
		CENTERPROPOGATE
		The low_res argument has no effect on the CENTERPROPOGATE option
		'''
		if method not in Box.CENTERMETHODS:
			print "error, you called center using an unknown method:",method
			return 0
		
			
		if method == Box.CENTEROFMASS:
			if low_res == True:
				# WONT WORK there is no self.autoboxer
				image = self.get_small_box_image(self.autoboxer.get_template_radius(),self.autoboxer.get_best_shrink())
				ali = image.calc_center_of_mass()
				dx = -int((ali[0]+0.5-image.get_xsize()/2))*extrasomething.get_best_shrink()
				dy = -int((ali[1]+0.5-image.get_ysize()/2))*extrasomething.get_best_shrink()
			else:
				image = self.get_box_image()
				ali = image.calc_center_of_mass()
				dx = -int((ali[0]+0.5-image.get_xsize()/2))
				dy = -int((ali[1]+0.5-image.get_ysize()/2))

		elif method == Box.CENTERACF:
			if low_res == True:
				# WONT WORK there is no self.autoboxer
				image = self.get_small_box_image(self.autoboxer.get_template_radius(),self.autoboxer.get_best_shrink())
				ccf  = image.calc_ccf(None)
				trans = ccf.calc_max_location_wrap(-1,-1,-1)
				dx = trans[0]/2*extrasomething.get_best_shrink()
				dy = trans[1]/2*extrasomething.get_best_shrink()
			else:
				image = self.get_box_image()
				ccf  = image.calc_ccf(None)
				trans = ccf.calc_max_location_wrap(-1,-1,-1)
				dx = trans[0]/2
				dy = trans[1]/2
		
		elif method == Box.CENTERPROPAGATE:
			template = extrasomething
			image =self.get_box_image()
			ccf  = image.calc_ccf(template)
			#sig = image.calc_fast_sigma_image(None)
			#ccf.div(sig)
			trans = ccf.calc_max_location_wrap(-1,-1,-1)
			dx = trans[0]
			dy = trans[1]
			
		#print "here we are",dx,dy
		
		self.xcorner += dx
		self.ycorner += dy
				
		# have to calculate offsets here
		if low_res == True and not method == Box.CENTERPROPAGATE:
			self.correct_resolution_centering(extrasomething.get_best_shrink(),False)
	
		self.update_box_image()
		self.changed = True
		
		return 1
		
	def correct_resolution_centering(self,shrink,update=True):
		
		nx = self.get_box_image().get_xsize()
		smallx = int(nx)/shrink
		ny = self.get_box_image().get_ysize()
		smally = int(ny)/shrink
			
		difx = int(shrink*int(smallx/2.0+0.5)-int(nx/2.0+0.5))
		dify = int(shrink*int(smally/2.0+0.5)-int(ny/2.0+0.5))
		self.xcorner += difx
		self.ycorner += dify
		
		print "correction",difx,dify
		
		if update and (difx != 0 or dify != 0):
			self.update_box_image()
			self.changed = True
			
	def update_params(self,autoboxer,center=False,force=False):
		'''
		Updates internally stored parameters, currently works only for SwarmAutoBoxer, but
		have attempted to lay basic framework if in future we use a different autoboxer which
		requires its own parameters
		'''
		if self.isdummy:
			return 0
		
		correlation = self.get_flcf_image(autoboxer)
		if correlation == None:
			print 'error, can not update the parameters of a Box because the Boxable has no correlation image'
			return 0
		
		if isinstance(autoboxer,SwarmAutoBoxer):
			shrink = autoboxer.get_best_shrink()
			invshrink = 1/shrink
	
			# the central coordinates of the box in terms of the shrunken correlation image
			x = (self.xcorner+self.xsize/2.0)*invshrink
			y = (self.ycorner+self.ysize/2.0)*invshrink
			
			#the search radius is used in correlation space - it limits the radial distance
			# up to which 'profile' data can be accrued
			# it is currently half the box_size in terms of the correlation image's dimensions
			searchradius = autoboxer.get_search_radius()
		
			peak_location = BoxingTools.find_radial_max(correlation,int(x),int(y),searchradius )
			peak_location2 = BoxingTools.find_radial_max(correlation,peak_location[0],peak_location[1],searchradius )
			if (peak_location != peak_location2):
				# this represents a troubling condition
				# setting box.correlationscore is the flag that other functions can act on in order to exclude
				# this box from consideration
				self.correlationscore = None
				if not force :
					#print "Error, peak location unrefined"
					return 0
		
			# store the peak location
			self.corx = peak_location[0]
			self.cory = peak_location[1]
		
			# store the correlation value at the correlation max
			self.correlationscore = correlation.get(self.corx,self.cory)
		
			# store the profile
			self.optprofile = BoxingTools.get_min_delta_profile(correlation,self.corx,self.cory, searchradius )
			# center on the correlation peak
			if (center):
				self.xcorner = self.corx*shrink-self.xsize/2.0
				self.ycorner = self.cory*shrink-self.ysize/2.0
				self.changed = True
			
			return 1
			
		else:
			print 'error, the autoboxer you are using is not currently known by the Box class'
			return 0
		
class TrimBox:
	'''
	A trimmed down version of a box
	'''
	def __init__(self,box):
		self.xcorner = box.xcorner			# the xcorner - bottom left
		self.ycorner = box.ycorner			# the ycorner - bottom left
		self.xsize = box.xsize				# the xsize of the box
		self.ysize = box.ysize				# the ysize of the box
		self.isref = box.isref				# a flag that can be used to tell if the box is being used as a reference
		self.changed = box.changed			# a flag signalling the box has changed and display needs updatin
		self.isanchor = box.isanchor		# a flag signalling the box has changed and display needs updatin
		self.TS = box.TS					# a time stamp flag
		self.image_name = box.image_name
		self.moved = box.moved
		self.origxcorner = box.origxcorner
		self.origycorner = box.origycorner
		self.ismanual = box.ismanual
		
	def change_box_size(self,box_size):
		'''
		Fixme, this is just a copy of Box.change_box_size
		'''
		if box_size != self.xsize or box_size != self.ysize:
			self.xcorner -= (box_size-self.xsize)/2
			self.xsize = box_size
			self.ycorner -= (box_size-self.ysize)/2
			self.ysize = box_size
			self.changed = True

class Cache:
	'''
	Provides a cache of static size (as defined by self.maxsize)
	As the cache grows objects are popped off the end of the self.cache tuple
	
	get the cache via get_cache - iterate through it to find your object
	add to the cache via add_to_cache
	reset the size of the cache using set_max_size
	'''
	def __init__(self):
		self.maxsize = 10
		self.cache = []

	def set_max_size(self,size):
		'''
		Will resize the cache if it is current larger than the new maxsize
		'''
		if len(self.cache) > size:
			self.cache = self.cache[0:size]

		self.maxsize = size
		
	def add_to_cache(self,object):
		'''
		Add 
		'''
		oldcache = self.cache
		self.cache = [object]
		self.cache.extend(oldcache)
		if len(self.cache) > self.maxsize:
			self.cache.pop(self.maxsize)
			if len(self.cache) != self.maxsize:
				print "error, the caching mechanism is not working correctly"
				
	def get_cache(self):
		return self.cache

	
class ExclusionImageCache:
	'''
	A singleton - used for caching exclusion images
	'''
	class __impl(Cache):
		""" Implementation of the singleton interface """

		def __init__(self):
			Cache.__init__(self)
			
		def get_image(self,image_name,xsize,ysize):
			excImage = None
			# first see if the object is already stored
			for object in self.get_cache():
				if object.get_input_image_name() == image_name:
					excImage = object
					break;
				
				
			if excImage == None:
				#if we make it here the cfimage is not cached
				#print "had to cache a cf image for",image_name
				excImage = ExclusionImage(image_name)
				self.add_to_cache(excImage)
			#else: print "found a cached cf image for",image_name
				
			
			image = excImage.get_image_carefully(xsize,ysize)
			if image != None:
				return image
			else:
				print "there was an error getting the image in ExclusionImageCache"
				return None
		
		
	# storage for the instance reference	
	__instance = None

	def __init__(self):
		""" Create singleton instance """
		# Check whether we already have an instance
		if ExclusionImageCache.__instance is None:
			# Create and remember instance
			ExclusionImageCache.__instance = ExclusionImageCache.__impl()
	
	def __getattr__(self, attr):
		""" Delegate access to implementation """
		return getattr(self.__instance, attr)

	def __setattr__(self, attr, value):
		""" Delegate access to implementation """
		return setattr(self.__instance, attr, value)
	
class ExclusionImage:
	def __init__(self,image_name):
		self.image_name = image_name
		self.image = None
		
		try:
			# we may have the image already on disk, if so parse it
			# the image on disk is likely up to date but not necessarily so
			self.image = get_idd_key_entry(self.image_name,"exclusion_image")
			#self.image = EMData(self.ouputimage_name)
			#print "I read the image",self.ouputimage_name
		except:
			#print "could not read", self.ouputimage_name 
			pass
		
	def get_input_image_name(self):
		return self.image_name
		
	def get_xsize(self):
		if self.image != None: return self.image.get_xsize()
		else: return 0
		
	def get_ysize(self):
		if self.image != None: return self.image.get_ysize()
		else: return 0
	
	def __update_image(self,xsize,ysize):
		'''
		Updates the image using the function arguments
		If they match current parameters than nothing happens - the correct image is already cached
		'''
		
		if self.image == None:
			self.image = EMData(xsize,ysize)
		else:
			# if the image already exists then we must retain the information in it by scaling and resizing it
			oldxsize = self.get_xsize()
			oldysize = self.get_ysize()
			r = Region( (oldxsize-xsize)/2, (oldysize-ysize)/2,xsize,ysize )
			#print "clipping to",(oldxsize-xsize)/2, (oldysize-ysize)/2,xsize,ysize
			scale = float(xsize)/float(oldxsize)
			
			# the order in which you clip and scale is dependent on whether or not scale is > 1
			if scale > 1:
				# if it's greater than one than clip (enlargen the image) first
				self.image.clip_inplace(r)
				# then scale the pixels
				self.image.scale(float(xsize)/float(oldxsize))
			else:
				# if it's less than one scale first so that we retain the maximum amount of the pixel information
				self.image.scale(float(xsize)/float(oldxsize))
				self.image.clip_inplace(r)
			
			#set_idd_key_entry(self.image_name,"exclusion_image",self.image) # not sure if this is necessary
			
		#else:
			#print "doing nothing to currently stored small image in CoarsenedFlattenedImage"
			
	def get_image(self):
		'''
		Should only be called if you know the stored image is up to date
		'''
		return self.image
	
	
	def get_image_carefully(self,xsize,ysize):
		
		if self.image == None or not self.query_params_match(xsize,ysize):
			#print "regenerating cf image"
			self.__update_image(xsize,ysize)
		#else: print "cf image is up to date"
		
		return self.get_image()
	
	def query_params_match(self,flattenradius,shrink):
		try:
			if xsize != self.get_xsize() or ysize != self.get_ysize():
				return False
			else: return True
		except: return False # exception will be thrown if self.smallimage = None

class SincBlackmanSubsamplerCache:
	'''
	A singleton - usef for caching coarsened-flattened images
	'''
	class __impl(Cache):
		""" Implementation of the singleton interface """

		def __init__(self):
			Cache.__init__(self)
			
		def get_image(self,image_name,params_mediator):
			ps_image = None
			# first see if the object is already stored
			for object in self.get_cache():
				if object.get_input_image_name() == image_name:
					ps_image = object
					break;
				
				
			if ps_image == None:
				#if we make it here the cfimage is not cached
				#print "had to cache a cf image for",image_name
				ps_image = SincBlackmanSubsampledImage(image_name)
				self.add_to_cache(ps_image)
			#else: print "found a cached cf image for",image_name
				
			
			image = ps_image.get_image_carefully(params_mediator)
			if image != None:
				return image
			else:
				print "there was an error getting the image in CoarsenedFlattenedImageCache"
				return None
		
	# storage for the instance reference	
	__instance = None

	def __init__(self):
		""" Create singleton instance """
		# Check whether we already have an instance
		if SincBlackmanSubsamplerCache.__instance is None:
			# Create and remember instance
			SincBlackmanSubsamplerCache.__instance = SincBlackmanSubsamplerCache.__impl()
	
	def __getattr__(self, attr):
		""" Delegate access to implementation """
		return getattr(self.__instance, attr)

	def __setattr__(self, attr, value):
		""" Delegate access to implementation """
		return setattr(self.__instance, attr, value)

class SincBlackmanSubsampledImage:
	'''
	Pawel Penczek optimal subsampling approach is encapsulated in this class
	'''
	def __init__(self,image_name):
		self.smallimage = None		# a small copy of an image which has had its background flattened
		self.image_name = image_name
		
		try:
			# we may have the image already on disk, if so parse it
			# the image on disk is likely up to date but not necessarily so
			self.smallimage = get_idd_key_entry(self.image_name,"subsampled_image")
			#print "I read the image",self.ouputimage_name
		except:
			#print "could not read", self.ouputimage_name 
			pass
		
	def get_input_image_name(self):
		return self.image_name
		
	def get_creation_ts(self):
		return	self.smallimage.get_attr("creation_time_stamp")
	
	def get_subsample_rate(self):
		return self.smallimage.get_attr("subsample_rate")
	
	def get_frequency_cutoff(self):
		return self.smallimage.get_attr("frequency_cutoff")
	
	def get_template_min(self):
		return self.smallimage.get_attr("template_min")
	
	def __update_image(self,params_mediator):
		'''
		Updates the image using the parameters that are deduced from the params_mediator
		'''
		subsample_rate = params_mediator.get_best_shrink()
		template_min = params_mediator.get_template_min()
		frequency_cutoff = params_mediator.get_frequency_cutoff()
		bic = BigImageCache()
		image = bic.get_image(self.image_name)
		
		sb = Util.sincBlackman(template_min, frequency_cutoff,1999) # 1999 taken directly from util_sparx.h

		self.smallimage = image.downsample(sb,1.0/subsample_rate)

		self.smallimage.set_attr("subsample_rate",subsample_rate)
		self.smallimage.set_attr("frequency_cutoff",frequency_cutoff)
		self.smallimage.set_attr("template_min",template_min)
		self.smallimage.set_attr("creation_time_stamp",gm_time_string())

		set_idd_key_entry(self.image_name,"subsampled_image",self.smallimage)
				
		#else:
			#print "doing nothing to currently stored small image in CoarsenedFlattenedImage"
			
	def get_image(self):
		'''
		Should only be called if you know the stored image is up to date
		'''
		return self.smallimage
	
	def get_image_carefully(self,params_mediator):
		
		
		if self.smallimage == None or not self.query_params_match(params_mediator):
			#print "regenerating sb image"
			self.__update_image(params_mediator)
		#else: print "sb image is up to date"
		
		return self.get_image()
	
	def query_params_match(self,params_mediator):
		subsample_rate = params_mediator.get_best_shrink()
		template_min = params_mediator.get_template_min()
		frequency_cutoff = params_mediator.get_frequency_cutoff()
		try:
			if subsample_rate != self.get_subsample_rate() or template_min != self.get_template_min() or frequency_cutoff != self.get_frequency_cutoff():
				return False
			else: return True
		except: return False # exception will be thrown if self.smallimage = None

class ImageProcParamsMediator:
	'''
	A mediator class - coordinates the requests of the various image processing tasks, as embodied 
	in the Cache classes, which are currently retrieved from the SwarmAutoBoxer
	'''
	def __init__(self,parent):
		if not isinstance(parent,AutoBoxer):
			print "Error, the ImageProcParamsMediator can only be instantiated with an AutoBoxer as its parent"
			exit(1)
		
		self.parent = parent
		
	def get_best_shrink(self):
		return self.parent.get_best_shrink()
	
	def get_template_radius(self):
		return self.parent.get_template_radius()
	
	def get_template_object(self):
		return self.parent.get_template_object()

	def get_template_min(self):
		return 20
	
	def get_frequency_cutoff(self):
		return 0.12

class CoarsenedFlattenedImageCache:
	'''
	A singleton - usef for caching coarsened-flattened images
	'''
	class __impl(Cache):
		""" Implementation of the singleton interface """

		def __init__(self):
			Cache.__init__(self)
			
		def get_image(self,image_name,params_mediator):
			cfImage = None
			# first see if the object is already stored
			for object in self.get_cache():
				if object.get_input_image_name() == image_name:
					cfImage = object
					break;
				
				
			if cfImage == None:
				#if we make it here the cfimage is not cached
				#print "had to cache a cf image for",image_name
				cfImage = CoarsenedFlattenedImage(image_name)
				self.add_to_cache(cfImage)
			#else: print "found a cached cf image for",image_name
				
			
			image = cfImage.get_image_carefully(params_mediator)
			if image != None:
				return image
			else:
				print "there was an error getting the image in CoarsenedFlattened"
				return None

	# storage for the instance reference	
	__instance = None

	def __init__(self):
		""" Create singleton instance """
		# Check whether we already have an instance
		if CoarsenedFlattenedImageCache.__instance is None:
			# Create and remember instance
			CoarsenedFlattenedImageCache.__instance = CoarsenedFlattenedImageCache.__impl()
	
	def __getattr__(self, attr):
		""" Delegate access to implementation """
		return getattr(self.__instance, attr)

	def __setattr__(self, attr, value):
		""" Delegate access to implementation """
		return setattr(self.__instance, attr, value)

class CoarsenedFlattenedImage:
	def __init__(self,image_name):
		self.smallimage = None		# a small copy of an image which has had its background flattened
		self.image_name = image_name
		
		try:
			# we may have the image already on disk, if so parse it
			# the image on disk is likely up to date but not necessarily so
			self.smallimage = get_idd_key_entry(self.image_name,"coarse_flat_image")
			#print "I read the image",self.ouputimage_name
		except:
			#print "could not read", self.ouputimage_name 
			pass
		
	def get_input_image_name(self):
		return self.image_name
		
	def get_creation_ts(self):
		return	self.smallimage.get_attr("creation_time_stamp")
	
	def get_flatten_radius(self):
		return self.smallimage.get_attr("flatten_radius")
	
	def get_shrink_factor(self):
		return self.smallimage.get_attr("shrink_factor")
	
	def __update_image(self,params_mediator):
		'''
		Updates the image using the function arguments
		If they match current parameters than nothing happens - the correct image is already cached
		'''
		bic = BigImageCache()
		image = bic.get_image(self.image_name)

		flattenradius = params_mediator.get_template_radius()
		shrink =  params_mediator.get_best_shrink()

		self.smallimage = image.process("math.meanshrink",{"n":shrink})
		self.smallimage.process_inplace("filter.flattenbackground",{"radius":flattenradius})
			
		self.smallimage.set_attr("flatten_radius",flattenradius)
		self.smallimage.set_attr("shrink_factor",shrink)
		self.smallimage.set_attr("creation_time_stamp",gm_time_string())
		
		set_idd_key_entry(self.image_name,"coarse_flat_image",self.smallimage)
				
		#else:
			#print "doing nothing to currently stored small image in CoarsenedFlattenedImage"
			
	def get_image(self):
		'''
		Should only be called if you know the stored image is up to date
		'''
		return self.smallimage
	
	
	def get_image_carefully(self,params_mediator):

		if self.smallimage == None or not self.query_params_match(params_mediator):
			#print "regenerating cf image"
			self.__update_image(params_mediator)
		#else: print "cf image is up to date"
		
		return self.get_image()
	
	def query_params_match(self,params_mediator):
		try:
			if self.get_flatten_radius() != params_mediator.get_template_radius() or self.get_shrink_factor() != params_mediator.get_best_shrink():
				return False
			else: return True
		except: return False # exception will be thrown if self.smallimage = None


class InverseSigmaImageCache:
	class __impl(Cache):
		""" A cache for storing big images """

		def __init__(self):
			Cache.__init__(self)

	
		def get_image(self,image_name,params_mediator,forceupdate=False):
			# this loop takes care of things if the image is cached
			object = None
			for sigmaImage in self.get_cache():
				if sigmaImage.get_image_name() == image_name:
					object = sigmaImage
					break
				
			
			# if we make it here the image is not cached
			if object == None:
				#print "I am generating a big image in the cache for",image_name
				object = InverseSigmaImage(image_name)
				self.add_to_cache(object)
			#else: print "I am returning a big image from the cache"
				
			return object.get_image_carefully(params_mediator,forceupdate)
			
			
	# storage for the instance reference
	__instance = None

	def __init__(self):
		""" Create singleton instance """
		# Check whether we already have an instance
		if InverseSigmaImageCache.__instance is None:
			# Create and remember instance
			InverseSigmaImageCache.__instance = InverseSigmaImageCache.__impl()
	
	def __getattr__(self, attr):
		""" Delegate access to implementation """
		return getattr(self.__instance, attr)

	def __setattr__(self, attr, value):
		""" Delegate access to implementation """
		return setattr(self.__instance, attr, value)
	
	
class InverseSigmaImage:
	def __init__(self,image_name):
		self.image_name = image_name
		self.image = None
		
		#try:
			## we may have the image already on disk, if so parse it
			## the image on disk is likely up to date but not necessarily so
			#self.image = get_idd_key_entry(self.image_name,"coarse_sigma_image")
			##print "I read the image",self.ouputimage_name
		#except:
			##print "could not read", self.ouputimage_name 
			#pass
		
	
	def get_image_name(self):
		return self.image_name
	
	def get_flatten_radius(self):
		return self.image.get_attr("flatten_radius")
	
	def get_shrink_factor(self):
		return self.image.get_attr("shrink_factor")
	
	def __update_image(self,params_mediator):
		flattenradius = params_mediator.get_template_radius()
		shrinkfactor = params_mediator.get_best_shrink()
		
		cache = SubsamplerCache
		image = cache.get_image(self.image_name,params_mediator)
		tmp = EMData(flattenradius*2,flattenradius*2)
		tmp.process_inplace("testimage.circlesphere")
		self.image = image.calc_fast_sigma_image(tmp)
		self.image.set_attr("flatten_radius",flattenradius)
		self.image.set_attr("shrink_factor",shrinkfactor)
		
		self.image.process_inplace("math.invert.carefully",{"zero_to":1.0})
		
		#set_idd_key_entry(self.image_name,"coarse_sigma_image",self.image)
		
		return self.image
	
	def get_image_carefully(self,params_mediator,forceupdate=False):

		
		action = False
		if forceupdate == True: action = True
		elif self.image == None: action = True
		elif params_mediator.get_template_radius() != self.get_flatten_radius() or params_mediator.get_best_shrink() != self.get_shrink_factor(): action = True
		
		if action: self.__update_image(params_mediator)
		
		return self.image

class BinaryCircleImageCache:
	'''
	A cache of white binary circles of a certain radius - used to make erasing fast
	Only ever stores these images in memory
	'''
	class __impl(Cache):
		""" A cache for storing big images """

		def __init__(self):
			Cache.__init__(self)

	
		def get_image(self,circleradius):
			# this loop takes care of things if the image is cached
			object = None
			for circleImage in self.get_cache():
				if circleImage.get_circle_radius() == circleradius :
					object = circleImage
					break

			# if we make it here the image is not cached
			if object == None:
				#print "I am generating a big image in the cache for",image_name
				object = BinaryCircleImage(circleradius)
				self.add_to_cache(object)
			#else: print "I am returning a big image from the cache"
				
			return object.get_image()
			
	# storage for the instance reference
	__instance = None

	def __init__(self):
		""" Create singleton instance """
		# Check whether we already have an instance
		if BinaryCircleImageCache.__instance is None:
			# Create and remember instance
			BinaryCircleImageCache.__instance = BinaryCircleImageCache.__impl()
	
	def __getattr__(self, attr):
		""" Delegate access to implementation """
		return getattr(self.__instance, attr)

	def __setattr__(self, attr, value):
		""" Delegate access to implementation """
		return setattr(self.__instance, attr, value)
		

class BinaryCircleImage:
	def __init__(self,circleradius):
		self.image = EMData(2*circleradius+1,2*circleradius+1)
		self.image.process_inplace("testimage.circlesphere")
		self.image.set_attr("circle_radius",circleradius)

	def get_circle_radius(self):
		return self.image.get_attr("circle_radius")
	
	def get_image(self):
		return self.image

class BigImageCache:
	class __impl(Cache):
		""" A cache for storing big images """

		def __init__(self):
			Cache.__init__(self)
			self.set_max_size(1)

	
		def get_image(self,image_name):
			# this loop takes care of things if the image is cached
			object = None
			for bigImage in self.get_cache():
				if bigImage.get_image_name() == image_name:
					object = bigImage
					break
				
			
			# if we make it here the image is not cached
			if object == None:
				#print "I am generating a big image in the cache for",image_name
				object = BigImage(image_name)
				self.add_to_cache(object)
			#else: print "I am returning a big image from the cache"
				
			return object.get_image()
			
	# storage for the instance reference
	__instance = None

	def __init__(self):
		""" Create singleton instance """
		# Check whether we already have an instance
		if BigImageCache.__instance is None:
			# Create and remember instance
			BigImageCache.__instance = BigImageCache.__impl()
	
	def __getattr__(self, attr):
		""" Delegate access to implementation """
		return getattr(self.__instance, attr)

	def __setattr__(self, attr, value):
		""" Delegate access to implementation """
		return setattr(self.__instance, attr, value)
	
class BigImage:
	def __init__(self,image_name):
		self.image_name = image_name
		self.image = None
	
	def get_image_name(self):
		return self.image_name
	
	def get_image(self):
		if self.image == None:
			self.image = EMData(self.image_name)
			self.image.process_inplace("normalize.edgemean") # this seams to be the normal behavior
			
		return self.image
		

#SubsamplerCache = CoarsenedFlattenedImageCache()
SubsamplerCache = SincBlackmanSubsamplerCache()
class FLCFImageCache:
	'''
	A singleton - used for caching flcf images
	'''
	class __impl(Cache):
		""" Implementation of the singleton interface """

		def __init__(self):
			Cache.__init__(self)


		def get_image(self,image_name,params_mediator):
			
			flcfImage = None
			# first see if the object is already stored
			for object in self.get_cache():
				if object.get_input_image_name() == image_name:
					flcfImage = object
					break;
				
				
			if flcfImage == None:
				# if we make it here the flcfimage is not cached
				#print "generated flcf for",image_name
				flcfImage = FLCFImage(image_name)
				self.add_to_cache(flcfImage)
			#else: print "found flcf for",image_name
				
			
			image = flcfImage.get_image_carefully(params_mediator)
			if image != None:
				return image
			else:
				print "there was an error getting the image"
				return None
		
		# storage for the instance reference	
	__instance = None

	def __init__(self):
		""" Create singleton instance """
		# Check whether we already have an instance
		if FLCFImageCache.__instance is None:
			# Create and remember instance
			FLCFImageCache.__instance = FLCFImageCache.__impl()
	
	def __getattr__(self, attr):
		""" Delegate access to implementation """
		return getattr(self.__instance, attr)

	def __setattr__(self, attr, value):
		""" Delegate access to implementation """
		return setattr(self.__instance, attr, value)



class FLCFImage:
	def __init__(self,image_name):
		self.flcfimage = None	# this is the flcf image
		self.image_name=image_name # we must store this it's used externally to determine if the FLCFImage is cached
		
		try: # try to read the image from disk - it may already exist and save us lots of time
			self.flcfimage = get_idd_key_entry(self.image_name,"flcf_image")
		except:
			# the image doesn't exist, that's okay
			pass
		
	def query_params_match(self,params_mediator):
		flatten_radius = params_mediator.get_template_radius()
		shrink_factor = params_mediator.get_best_shrink()
		template = params_mediator.get_template_object()
		#print cfimage.get_attr("creation_time_stamp"), self.get_cfi_ts(),"template",template.get_template_ts(),self.get_template_ts()
		try:
			if flatten_radius != self.get_sigma_image_flatten_radius() or shrink_factor != self.get_sigma_image_shrink_factor() or template.get_template_ts() != self.get_template_ts():
				#print "params did not match"
				return False
			else: return True
		except: return False
	
	def get_input_image_name(self):
		return self.image_name
	
	def get_output_image_name(self):
		return self.outputimage_name
	
	def get_sigma_image_flatten_radius(self):
		return self.flcfimage.get_attr("get_sigma_image_flatten_radius")
	
	def get_sigma_image_shrink_factor(self):
		return self.flcfimage.get_attr("get_sigma_image_shrink_factor")
	
	def get_template_ts(self):
		'''
		get template time stamp
		'''
		return self.flcfimage.get_attr("template_time_stamp")
		
	def get_image(self):
		'''
		Returns the currently stored flcfimage
		'''
		return self.flcfimage
	
	def __update_image(self,params_mediator):

		sicache = InverseSigmaImageCache()
		inv_sigma_image = sicache.get_image(self.image_name,params_mediator)
		
		cfcache = SubsamplerCache
		cfimage = cfcache.get_image(self.image_name,params_mediator)
		
		flatten_radius = params_mediator.get_template_radius()
		shrink_factor = params_mediator.get_best_shrink()
		template = params_mediator.get_template_object()
		
		self.flcfimage = cfimage.calc_ccf( template.get_template() )
		self.flcfimage.mult(inv_sigma_image)
		self.flcfimage.process_inplace("xform.phaseorigin.tocenter")
		self.flcfimage.set_attr("template_time_stamp",template.get_template_ts())
		self.flcfimage.set_attr("get_sigma_image_shrink_factor",shrink_factor)
		self.flcfimage.set_attr("get_sigma_image_flatten_radius",flatten_radius)
		
		set_idd_key_entry(self.image_name,"flcf_image",self.flcfimage)
		
	def get_image(self):
		return self.flcfimage
	
	def get_image_carefully(self,params_mediator):
		'''
		Checks to see if the arguments are the right types
		Then checks to see if the currently store correlation image is up to date
		If it's not up to date or it doesn't exist then it is (re)generated
		Then it is returned
		'''
		template = params_mediator.get_template_object()
		if not isinstance(template,SwarmTemplate) and not isinstance(template,TrimSwarmTemplate):
			print "you can't call genFLCF on an object that is not a SwarmTemplate"
			return None
		
		action = False
		if self.flcfimage != None:
			if not self.query_params_match(params_mediator):
				action = True
		else: action = True
		
		if action:
			#print "generating correlation image"
			self.__update_image(params_mediator)
		#else: print "returning cached correlation image"
		return self.get_image()

class Boxable:
	UNERASE = 'Unerase'
	ERASE = 'Erase'
	
	GREAT = '4 - Excellent'
	GOOD = '3'
	AVERAGE = '2 - Average'
	POOR = '1'
	EXCLUDE = '0 - Exclude'
	QUALITY_META_DATA = [EXCLUDE,POOR,AVERAGE,GOOD,GREAT]
	QUALITY_META_DATA_MAP = {}
	QUALITY_META_DATA_MAP[EXCLUDE] = 0
	QUALITY_META_DATA_MAP[POOR] = 1
	QUALITY_META_DATA_MAP[AVERAGE] = 2
	QUALITY_META_DATA_MAP[GOOD] = 3
	QUALITY_META_DATA_MAP[GREAT] =  4
	
	def __init__(self,image_name,parent,autoboxer=None):
		self.parent = parent		# keep track of the parent in case we ever need it
		self.boxes = []				# a list of boxes
		self.refboxes = []			# a list of boxes
		
		# a list of deleted auto boxes - should be cleared everytime an autobox occurs, and store interactively deleted boxes
		# This is because when the user unerases things the design of e2boxer demands that boxes in unerased regions should 
		# return...
		self.deleted_auto_boxes = [] 
		
		self.box_size = -1			#  the box_size
		self.image_name = image_name
		
		self.fpshrink = -1
		self.exclusionimage = None
		self.template = None
		self.refcache = []
		self.template_ts = -1 # a template time stamp, used to avoid unecessarily regenerating the template in self.auto_box
		self.autoboxer_state_ts = -1 # and autoboxer time stamp, used to avoid unecessary autoboxing, and to force autoboxing when appropriate
		self.autoboxer_id = -1 # Stores the unique ID of the autoboxer - this is to facilitate having many autoboxers in the project data base
		
		self.autoboxer = autoboxer
		self.__frozen = False
		self.__quality = Boxable.QUALITY_META_DATA_MAP[Boxable.AVERAGE] # this makes it the number, not the string
		
		try:
			eicache = ExclusionImageCache()
			#print "reading exclusion image",excimage_name
			self.exclusionimage = eicache.get_image(self.image_name,self.get_small_image().get_xsize(),self.get_small_image().get_ysize())
		except: pass
		
		self.get_db_stamps()
		self.reload_boxes()
		self.get_frozen_from_db()	
		self.get_quality_from_db()
		try:
			self.check_store_image_tag_db()	
		except: pass
		if autoboxer == None: self.get_auto_boxer_from_db()
	
	def reload_boxes(self):
		'''
		You might call this if the box size was changed in the main interface, or in the constructor
		'''
		self.boxes = []
		self.deleted_auto_boxes = []
		try: 
			self.get_auto_selected_from_db(True)
		except: pass # this probably means there is a project_db but it doesn't store any autoboxing results from this image
		#print "got auto, now have", len(self.boxes)
		try:
			self.get_manual_boxes_from_db()
		except: pass
		#print "got manual, now have", len(self.boxes)
		try:
			self.get_references_from_db()	
		except: pass
		#print "got references, now have", len(self.boxes)
		
	def resize_boxes(self,boxes,box_size):
		if boxes == None:
			return
		for box in boxes:
			box.change_box_size(box_size)
		
	def resize_moved_box_data(self,data,box_size,oldbox_size):
		adj = (box_size-oldbox_size)/2
		for d in data:
			d[0] -= adj
			d[1] -= adj
			d[2] -= adj
			d[3] -= adj
	
			
	def set_box_size(self,box_size):
		action = True
		if self.box_size != box_size and self.box_size != -1:
			# if you call set box size and there are boxes with this box size then you should
			# really be calling change box size
			if len(self.boxes) != 0:
				print "warning, in calling set box size, change box size is being called"
				self.change_box_size(box_size)
				action = False
	
		if action:
			self.box_size = box_size
	
	def change_box_size(self,box_size):
		self.box_size = box_size
		self.deleted_auto_boxes = []
		oldxsize  = None
		autoboxes = get_idd_key_entry(self.image_name,"auto_boxes")
		if oldxsize == None and autoboxes != None and len(autoboxes) != 0:
			oldxsize = autoboxes[0].xsize
		self.resize_boxes(autoboxes,box_size)
		set_idd_key_entry(self.image_name,"auto_boxes",autoboxes)
		
		manboxes = get_idd_key_entry(self.image_name,"manual_boxes")
		if oldxsize == None and manboxes != None and len(manboxes) != 0:
			oldxsize = manboxes[0].xsize
		self.resize_boxes(manboxes,box_size)
		set_idd_key_entry(self.image_name,"manual_boxes",manboxes)
		
		refboxes = get_idd_key_entry(self.image_name,"reference_boxes")
		if oldxsize == None and refboxes != None and len(refboxes) != 0:
			oldxsize = refboxes[0].xsize
		self.resize_boxes(refboxes,box_size)
		set_idd_key_entry(self.image_name,"reference_boxes",refboxes)
		
		
		movedboxes =  get_idd_key_entry_in_memory(self.image_name,"moved_boxes")
		if movedboxes == None: return
		if oldxsize == None and movedboxes != None and len(movedboxes) != 0:
			print "warning, changing box sizes, old movement information has been lost"
			movedboxes = []
		else:
			self.resize_moved_box_data(movedboxes,box_size,oldxsize)
		set_idd_key_entry_in_memory(self.image_name,"moved_boxes",movedboxes)
	
	def clear_and_cache(self,keepmanual=False):
		self.boxes = []
		self.deleted_auto_boxes = []
		
		project_db = EMProjectDB()
		try:
			data = project_db[self.get_dd_key()]
		except:
			data = {}
			
		if data == None:
			data = {}
		
		data["reference_boxes"] = []
		data["auto_boxes"] = []
		
		
		if keepmanual:
			self.get_manual_boxes_from_db()
		else :
			data["manual_boxes"] = []
	
	
		project_db.set_key_entry(self.get_dd_key(),data)
		
	def clear_and_reload_images(self):
		self.boxes = []
		self.deleted_auto_boxes = []
		try: 
			self.get_auto_selected_from_db(True)
		except: pass # this probably means there is a project_db but it doesn't store any autoboxing results from this image
		try:
			self.get_manual_boxes_from_db()	
		except: pass
		try:
			self.get_references_from_db()	
		except: pass
		
	def cache_exc_to_db(self):
		if self.exclusionimage != None:
			set_idd_key_entry(self.image_name,"exclusion_image",self.exclusionimage)
			merge_idd_key_entry_memory_to_disk(self.image_name,"moved_boxes")

	def center(self,method):
		if method == Box.CENTERACF or method == Box.CENTEROFMASS:
			extrasomething = self.autoboxer
		elif method == Box.CENTERPROPAGATE:
			extrasomething = self.autoboxer.get_high_res_template_image()
		else:
			print "error, the method you specified is unsupported in Boxable:",method
		
		for box in self.boxes:
			if not box.center(method,extrasomething,False):
				print "there was an error boxing"
				return 0
				
		return 1

	def set_stamps(self,autoboxer_state_ts,template_ts,autoboxer_unique_id):
		'''
		A convenience function for setting all three important time stamps/ids at once
		
		'''
		self.set_autoboxer_state_ts(autoboxer_state_ts)
		self.set_template_ts(template_ts)
		self.set_autoboxer_id(autoboxer_unique_id)

	def get_auto_boxer_state_ts(self):
		return self.autoboxer_state_ts
	
	def set_autoboxer_state_ts(self,autoboxer_state_ts):
		self.autoboxer_state_ts = autoboxer_state_ts
	
	def get_template_ts(self):
		return self.template_ts
	
	def set_template_ts(self,template_ts):
		self.template_ts = template_ts

	def get_autoboxer_id(self):
		return self.autoboxer_id
	
	def set_autoboxer_id(self,autoboxer_unique_id):
		self.autoboxer_id = autoboxer_unique_id
	
	def write_to_db(self):
		'''
		Writes fundamentally important information to the database
		'''
		project_db = EMProjectDB()
		
		try:
			data = project_db[self.get_dd_key()]
		except:
			data = {}
		
		if data == None:
			data = {}
			
		data["autoboxer_state_TS"] = self.get_auto_boxer_state_ts()
		data["template_TS"] = self.get_template_ts()
		data["frozen_state"] = self.__frozen
		data["quality"] = self.__quality
		data["autoboxer_unique_id"] = self.get_autoboxer_id()
		
		project_db.set_key_entry(self.get_dd_key(),data)
	
	def set_key_entryToIDD(self,key,object):
		set_idd_key_entry(self.image_name,key,object)
		
	def get_image_name(self):
		return self.image_name
	
	def set_autoboxer(self,autoboxer):
		self.autoboxer = autoboxer
	
	def extend_boxes(self,boxes):
		self.boxes.extend(boxes)
		
	def is_frozen(self):
		return self.__frozen
	
	def is_interactive(self):
		return not (self.is_frozen() or self.is_excluded())
	
	def toggle_frozen(self):
		self.__frozen = not self.__frozen
		
	def set_frozen(self,frozen):
		self.__frozen = frozen

	def get_dd_key(self):
		'''
		Get Database Dictionary Key
		'''
		return get_idd_key(self.image_name)
	
	def get_auto_boxer_from_db(self):
		
		project_db = EMProjectDB()
		data = project_db[self.get_dd_key()]
		self.autoboxer_id = data["autoboxer_unique_id"]
		trim_autoboxer = project_db[self.autoboxer_id]
		self.autoboxer = SwarmAutoBoxer(self.parent)
		self.autoboxer.become(trim_autoboxer)
		
	def set_quality(self,quality):
		if quality not in Boxable.QUALITY_META_DATA:
			print "error",quality,"is not in",Boxable.QUALITY_META_DATA,"can't do anything"
			return
		else:
			self.__quality =  Boxable.QUALITY_META_DATA_MAP[quality]
	
	def is_excluded(self):
		return self.__quality == Boxable.QUALITY_META_DATA_MAP[Boxable.EXCLUDE]
	
	def get_quality(self):
		return self.__quality
		
	def check_store_image_tag_db(self):
		#print "in cheack and Store image tag in db"
		project_db = EMProjectDB()
		data = project_db[self.get_dd_key()]
		newimagetag = get_file_tag(self.image_name)
		try:
			oldimagetag = data["image_tag"]
		except:
			#print "stored image tag for first time"
			data["image_tag"] = newimagetag
			project_db.set_key_entry(self.get_dd_key(),data)
			return
		
		if oldimagetag != newimagetag:
			print "warning with respect to",self.image_name,"- you are using information in the database that was generated using an image of type",newimagetag,"on an image of type",oldimagetag,". This will potentially cause problems if the images are not equivalent. Suggest renaming the image or boxing it in a separate directory"
			return
		else:
			#print "storing image tag"
			data["image_tag"] = newimagetag
			project_db.set_key_entry(self.get_dd_key(),data)

		
	def get_quality_from_db(self):
		self.__quality = get_idd_key_entry(self.image_name,"quality")
		
		if self.__quality == None:
			self.__quality = Boxable.QUALITY_META_DATA_MAP[Boxable.AVERAGE]
		
	def get_frozen_from_db(self):
		self.__frozen = get_idd_key_entry(self.image_name,"frozen_state")
		if self.__frozen == None:
			self.__frozen = False
		
	def get_db_stamps(self):
		self.autoboxer_state_ts = get_idd_key_entry(self.image_name,"autoboxer_state_TS")
		self.template_ts = get_idd_key_entry(self.image_name,"template_ts")
		self.autoboxer_id = get_idd_key_entry(self.image_name,"autoboxer_unique_id")
	
	def append_stored_auto_boxes(self,trimboxes):
		'''
		Sometimes this functionality is needed when the currently stored auto-selected
		boxes are not removed prior to the execution of autoboxing - in this case 
		the autoboxes stored in the database need to include what was already stored...
		'''
		self.deleted_auto_boxes = []
		for box in self.boxes:
			if not (box.ismanual or box.isref):
				trimboxes.append(TrimBox(box))
	
	def get_auto_selected_from_db(self,forcereadall=False):	
		
		self.deleted_auto_boxes = []
		
		trimboxes = get_idd_key_entry(self.image_name,"auto_boxes")
		if trimboxes == None or len(trimboxes) == 0:
			return 0
		
		movedboxes = get_idd_key_entry_in_memory(self.image_name,"moved_boxes") # this may potentially be None
		
		for trimbox in trimboxes:
			if trimbox.ismanual or trimbox.isref:
				print "error, the box was manual or it was a reference"
				continue
			if trimbox.changed or forcereadall:
				box = Box()
				
				# had to do conversion stuff so pickle would work
				box.become(trimbox)
				box.set_image_name(self.image_name)
				if forcereadall:
					box.changed = True
				if box.isref and not forcereadall:
					continue;
				elif box.isref:
					box.rorig = 0			# RGB red
					box.gorig = 0			# RGB green
					box.borig = 0			# RGB blue
					box.r = 0
					box.g = 0
					box.b = 0
				
				# The box may have a more correct centering
				if movedboxes != None and len(movedboxes) != 0:
					box.update_position_from_data(movedboxes)
				self.boxes.append(box)
		
		# Sometimes an exclusion area is added after the autoboxing has occured, in which case
		# autoboxes in the db will be in the excluded area and hence we have to make sure they 
		# are not included
		self.update_excluded_boxes()	
		
		#print "Added",len(self.boxes)-a," autoboxes"
	
	def get_references_from_db(self):
		#debug
		#a = len(self.boxes)
		
		refboxes = get_idd_key_entry(self.image_name,"reference_boxes")
		if refboxes == None or len(refboxes) == 0:
			return 0
		
		movedboxes = get_idd_key_entry_in_memory(self.image_name,"moved_boxes")# movedboxes is potentially None
		
		for trimbox in refboxes:
			box = Box()
			
			# had to do conversion stuff so pickle would work
			box.become(trimbox)
			box.set_image_name(self.image_name)
			box.changed = True
			box.rorig = 0			# RGB red
			box.gorig = 0			# RGB green
			box.borig = 0			# RGB blue
			box.r = 0
			box.g = 0
			box.b = 0
			
			# The box may have a more correct centering
			if movedboxes != None and len(movedboxes) != 0:
				box.update_position_from_data(movedboxes)
			self.boxes.append(box)

		#print "Added",len(self.boxes)-a," references"
	def get_manual_boxes_from_db(self):
		#a = len(self.boxes) # debug		
		manualboxes = get_idd_key_entry(self.image_name,"manual_boxes")
		if manualboxes == None or len(manualboxes) == 0:
			return 0
	
		for trimbox in manualboxes:
			box = Box()
			
			# had to do conversion stuff so pickle would work
			box.become(trimbox)
			box.set_image_name(self.image_name)
			box.changed = True
			box.rorig = 1			# RGB red
			box.gorig = 1			# RGB green
			box.borig = 1			# RGB blue
			box.r = 1
			box.g = 1
			box.b = 1
			self.boxes.append(box)
		#print "Added",len(self.boxes)-a," manual boxes"

	def get_coord_file_name(self):
		return strip_file_tag(self.image_name)+".box"
		
	def write_coord_file(self,box_size=-1,force=False,verbose=True):
		'''
		If box_size is -1 then the current box_size is used to write output
		If force is True then the old file with the same name is written over (as opposed to backed up)
		
		'''
		if len(self.boxes) == 0:
			print "no boxes to write, doing nothing. Image name is",self.image_name
		else:
			boxname = self.get_coord_file_name()
			if file_exists(boxname):
				if not force:
					f=file(boxname,'r')
					boxname_backup =  strip_file_tag(self.image_name)+str(time()) + ".box.bak"
					print "warning, found box name",boxname,"- am renaming it to", boxname_backup, "- use force to overwrite this behavior"
					fbak=file(boxname_backup,'w')
					fbak.writelines(f.readlines())
					fbak.close()
					f.close()
				else:
					remove_file(boxname)
				
			f=file(boxname,'w')
			
			if verbose: print "writing",self.num_boxes(),"box coordinates to file",boxname
			
			for box in self.boxes:
				if box_size != -1:
					# FOO - this will not work if the box dimensions are not equal...
					origbox_size = box.xsize
					if origbox_size != box.ysize:
						print "error, uniform box dimensions are not supported"
						return
					box.change_box_size(box_size)
						
				f.write(str(int(box.xcorner))+'\t'+str(int(box.ycorner))+'\t'+str(box.xsize)+'\t'+str(box.ysize)+'\n')
				
				if box_size != -1:
					# change it back
					box.change_box_size(origbox_size)

			f.close()
			
	def get_image_file_name(self,imageformat="hdf"):
		return strip_file_tag(self.image_name)+"."+imageformat

	def write_box_images(self,box_size=-1,force=False,imageformat="hdf",verbose=True):
		'''
		If box_size is -1 then the current box_size is used to write output
		If force is True then output is written over (if it already exists) - else an error is printed and nothing happens
		
		'''
		if len(self.boxes) == 0:
			print "no boxes to write, doing nothing. Image name is",self.image_name
		else:
			image_name = self.get_image_file_name(imageformat)
			
			if file_exists(image_name):
				if not force:
					print "warning, file already exists - ", image_name, " doing nothing. Use force to override this behavior"
					return
				else:
					remove_file(image_name)
			
			if verbose:	print "writing",self.num_boxes(),"boxed images to", image_name
			
			for box in self.boxes:
				if box_size != -1:
					# FOO - this will not work if the box dimensions are not equal...
					origbox_size = box.xsize
					if origbox_size != box.ysize:
						print "error, uniform box dimensions are not supported"
						return
					box.change_box_size(box_size)
						
				image = box.get_box_image()
				
				image.set_attr("original_x_corner",box.xcorner)
				image.set_attr("original_y_corner",box.ycorner)
				image.set_attr("original_x_size",box.xsize)
				image.set_attr("original_y_size",box.ysize)
				
				image.set_attr("originating_image_name",self.get_image_name())
				
				image.write_image(image_name,-1)
				
				if box_size != -1:
					box.change_box_size(origbox_size)

	def move_box(self,box,dx,dy,box_num):
		if box.ismanual:
			self.move_manual_box(box,dx,dy)
		
		box.move(dx,dy)

	def add_box(self,box):	
		if not isinstance(box,Box):
			print "You can not add a box to this box set if it is not of type Box"
			return;

		box.boxingobj = self
		
		if box.isref:
			box.rorig = 0			# RGB red
			box.gorig = 0			# RGB green
			box.borig = 0			# RGB blue
			box.r = 0
			box.g = 0
			box.b = 0
		else: # box.ismanual = True
			box.rorig = 1			# RGB red
			box.gorig = 1			# RGB green
			box.borig = 1			# RGB blue
			box.r = 1
			box.g = 1
			box.b = 1
			self.cache_manual_box(box)
		
		#print "adding box",box.xcorner,box.ycorner,box.xsize,box.ysize
		self.boxes.append(box)
		self.refboxes.append(box)
	
	def cache_manual_box(self,box):
		manualboxes = get_idd_key_entry(self.image_name,"manual_boxes")
		
		if manualboxes == None:
			manualboxes = []
		manualboxes.append(TrimBox(box))
		set_idd_key_entry(self.image_name,"manual_boxes",manualboxes)

	def delete_manual_box(self,box):
		
		manualboxes = get_idd_key_entry(self.image_name,"manual_boxes")
		
		if manualboxes == None or len(manualboxes) == 0:
			print "error, you can't move a manual box if there are none!"
			return
		
		found = False
		for j,b in enumerate(manualboxes):
			if b.xcorner == box.xcorner and b.ycorner == box.ycorner:
				manualboxes.pop(j)
				set_idd_key_entry(self.image_name,"manual_boxes",manualboxes)
				found = True
				break
		
		if not found:
			print "error, couldn't find the manual box you tried to delete, nothing happened"
			return
		
	
	def move_manual_box(self,box,dx,dy):
		manualboxes = get_idd_key_entry(self.image_name,"manual_boxes")
		
		if manualboxes == None or len(manualboxes) == 0:
			print "error, you can't move a manual box if there are none!"
			return
		
		found = False
		for j,b in enumerate(manualboxes):
			if b.xcorner == box.xcorner and b.ycorner == box.ycorner:
				b.xcorner += dx
				b.ycorner += dy
				found = True
				set_idd_key_entry(self.image_name,"manual_boxes",manualboxes)
				break
		
		if not found:
			print "error, couldn't find the manual box you tried to move, nothing happened"
			return
		
	
	def delete_box(self,i):
		tmp = self.boxes.pop(i)
		self.deleted_auto_boxes.append(tmp)
		if tmp.ismanual:
			self.delete_manual_box(tmp)
		#yuck, this is horribly inefficient
		for j,box in enumerate(self.refboxes):
			if box.isref and box.TS == tmp.TS:
				self.refboxes.pop(j)
				return True
		
		
		return False
	
	def delete_non_refs(self,updatedisplay=True):
		boxestodelete = []
		self.deleted_auto_boxes = []
		n = len(self.boxes)
		for m in range(n-1,-1,-1):
			box = self.boxes[m]
			if box.isref == False and box.ismanual == False:
				self.delete_box(m)
				boxestodelete.append(m)

		if updatedisplay:
			self.parent.delete_display_boxes(boxestodelete)
	
	def add_non_refs(self,boxes):
		'''
		Add boxes that are stored in eman1 format
		box[0] = xnorner, box[1] = ycorner, box[2] = xsize, box[3] = ysize
		'''
		for box in boxes:
			b = Box(box[0],box[1],box[2],box[3])
			b.set_image_name(self.image_name)
			b.isref = False
			b.changed = True
			self.boxes.append(b)

	def num_boxes(self):
		return len(self.boxes)
	
	def update_box_size(self,box_size):
		'''
		Updates only the box size and corner coordinates
		Switches the changed flag to True to trigger redisplay (but the calling function
		is responsible for knowing and testing for this)
		'''
		# do nothing if it's the same size as what we already have
		if  box_size == self.box_size: return
		
		for box in self.boxes:
			if box.xsize != box_size:
				box.xcorner -= (box_size-box.xsize)/2
				box.xsize = box_size
				box.changed = True
			if box.ysize != box_size:
				box.ycorner -= (box_size-box.ysize)/2
				box.ysize = box_size
				box.changed = True
			
			box.image = None
			box.footprint = None

		self.fprink = -1
		self.flattenimager = -1
		self.box_size = box_size
		
	def get_footprint_shrink(self):
		if self.fpshrink == -1:
			shrink = 1
			tn = self.box_size/2
			while ( tn >= 32 ):
				tn /= 2
				shrink *= 2
			self.fpshrink = shrink
		
		return self.fpshrink
		
	def get_best_shrink(self):
		'''
		FIXME - there should probably be a more well established framework for doing this
		At the moment it is possible that the self.autBoxer is actually None, which isn't good.	
		'''
		if self.autoboxer != None:
			return self.autoboxer.get_best_shrink()
		else:
			print 'warning, there is no autoboxer set, am not sure how to shrink, returning 1 as the shrink factor'
			return 1
		
	def get_correlation_image(self,autoboxer):
		
		cache = FLCFImageCache()
		correlation = cache.get_image(self.image_name,autoboxer.get_params_mediator())
		
		
		# the template time stamp may not be necessary ultimately, seeing as the correlation image
		# (which has its own template time stamp) is tied to this Boxable by the image_name itself
		template = autoboxer.get_template_object()
		template_ts = get_idd_key_entry(self.image_name,"template_ts")
		if template_ts != template.get_template_ts():	
			set_idd_key_entry(self.image_name,"template_ts",template.get_template_ts())
			self.template_ts = template.get_template_ts() # this shouldn't be necessary in future
			
		return correlation
	
	def get_small_image(self):
		#cache = CoarsenedFlattenedImageCache()
		cache = SubsamplerCache
		return  cache.get_image(self.image_name,self.autoboxer.get_params_mediator())
	
	def update_included_boxes(self):
		added_boxes = []
		added_ref_boxes = []
		invshrink = 1.0/self.get_best_shrink()
		exclusionimage = self.get_exclusion_image()
		n = len(self.deleted_auto_boxes)
		for i in range(n-1,-1,-1):
			box = self.deleted_auto_boxes[i]
			x = int((box.xcorner+box.xsize/2.0)*invshrink)
			y = int((box.ycorner+box.ysize/2.0)*invshrink)
			
			if ( exclusionimage.get(x,y) == 0):
				box.changed = True
				added_boxes.append(box)
				self.boxes.append(box)
				self.deleted_auto_boxes.pop(i)

				if box.isref: added_ref_boxes.append(box)
				
		return [added_boxes,added_ref_boxes]
	
	def update_excluded_boxes(self, useinternal=True,exclusionimage= None):
		'''
		
		'''
		if useinternal:
			exclusionimage = self.get_exclusion_image()

		lostboxes = []
			
		invshrink = 1.0/self.get_best_shrink()
		n = len(self.boxes)
		refs = []
		for i in range(n-1,-1,-1):
			box = self.boxes[i]
			x = int((box.xcorner+box.xsize/2.0)*invshrink)
			y = int((box.ycorner+box.ysize/2.0)*invshrink)
			
			if ( exclusionimage.get(x,y) != 0):
				lostboxes.append(i)

				box = self.boxes.pop(i)
				self.deleted_auto_boxes.append(box)
				if box.isref: refs.append(box)
	
		return [lostboxes,refs]
	
	def add_exclusion_particle(self,box):
		
		xx = box.xcorner+box.xsize/2-1
		yy = box.ycorner+box.ysize/2-1
		
		self.add_exclusion_area(None,xx,yy,box.xsize/2)
	
	def add_exclusion_area(self,UNUSEDtype,x,y,radius,flag=ERASE):
		'''
		UNUSEDtype was meant to be a flag for adding other exclusion areas like squares
		At the moment only circular exclusion areas can be written
		'''
		#print "Add exclusion area using",self.get_best_shrink()
		xx = int(x/self.get_best_shrink())
		yy = int(y/self.get_best_shrink())
		
		rr = int(radius/self.get_best_shrink())
		bciCache = BinaryCircleImageCache()
		mask = bciCache.get_image(rr)
		
		
		if flag == Boxable.ERASE:
			val = 0.1
		elif flag == Boxable.UNERASE:
			val = 0
		else:
			print "error - unknow flag:",flag,"doing nothing"
			return
		
		BoxingTools.set_region(self.get_exclusion_image(),mask,xx,yy,val)
		
	
	def get_exclusion_image(self,force=False):
		
		if self.exclusionimage == None or force:
			
			si = self.get_small_image()
			eicache = ExclusionImageCache()
			self.exclusionimage = eicache.get_image(self.image_name,si.get_xsize(),si.get_ysize())
	
		return self.exclusionimage
	
	def classify(self):
		
		# accrue all params
		n = self.autoboxer.optprofileradius+1
		for box in self.boxes:
			# set the force flag true or else the optprofile won't be set when the peak is 'faulty'
			box.update_params(self.autoboxer,False,True)
		
		v = []
		for box in self.boxes:
			#b = copy(box.optprofile[0:n])
			b = copy(box.optprofile)
			#for a in b: 
				#a = box[6]-a
			#print b
			v.append(b)
			
		cl = BoxingTools.classify(v,4)
		self.parent.update_box_colors(cl)
	
	def gen_ref_images(self):
		tmpimage = "tmpparticles.img"
		self.parent.writeBoxesTo(tmpimage)
		
		self.process = QtCore.QProcess()

		program = QtCore.QString("e2refine2d.py")
		args = QtCore.QStringList()
		args.append("--input="+tmpimage)
		args.append("--ncls=25")
		
		QtCore.QObject.connect(self.process, QtCore.SIGNAL("finished(int)"), self.process_finished)
		QtCore.QObject.connect(self.process, QtCore.SIGNAL("started()"), self.process_start)
		print self.process.start(program,args)

	def process_start(self):
		print "received process start signal"
		
	def box_selected(self,event,lc):
		#print "selected",lc[0]
		for box in self.boxes:
			if box.group == lc[0]:
				box.r = 1
				box.g = 1
				box.b = 1
				box.changed = True
			elif box.r == 1 and box.g == 1 and box.b == 1:
				box.r = box.rorig
				box.g = box.gorig
				box.b = box.borig
				box.changed = True
		self.imagemx2.setSelected(lc[0])
		self.parent.box_display_update()
		
	def process_finished(self,int):
		try:
			from emimage import EMImage
		except:
			print "Cannot import EMAN image GUI objects (emimage,etc.)"
			sys.exit(1)
		
		e = EMData().read_images("classes.init.hdf")
		self.imagemx2p = EMImage(e)
		self.imagemx2 = self.imagemx2p.child
		self.imagemx2.setmmode("app")
		QtCore.QObject.connect(self.imagemx2,QtCore.SIGNAL("mousedown"),self.box_selected)
		self.imagemx2p.show()
		
		ef = []
		for image in e:
			image.process_inplace("normalize.edgemean")
			if self.get_best_shrink() != 1:
				image = image.process("math.meanshrink",{"n":self.get_footprint_shrink()})	
			ef.append(image.make_footprint())
		
		for box in self.boxes:
			best = -1
			group = -1
			for i,g in enumerate(ef): 
				s = box.get_foot_print(self.get_footprint_shrink()).cmp("optvariance",g,{"matchfilt":1,"matchamp":1})
				# REMEMBER - cmp returns values that have potentially been negated - a smaller value is better
				if best == -1 or s < best:
					group = i
					best = s
			
			box.group = group
					
		#print scores
		
		print "received finish signal"

def merge_idd_key_entry_memory_to_disk(image_name,key):
	dbkey = get_idd_key(image_name)
	
	project_db = EMProjectDB()
	try:
		data_in_memory = project_db.get_key_entry_in_memory(dbkey)
	except:
#		print "cannot merge data from memory, it's not there :",image_name
		return
	
	try: data_in_memory[key]
	except:
#		print "data in memory didn't have key moved boxes"
		return
	
	try:
		data = project_db[dbkey]
	except:
		data = {}
		
	if data == None: data = {}
	
	try: data[key] = data_in_memory[key]
	except:
#		print "data in memory didn't have key moved boxes"
		return
	#project_db[self.get_dd_key()] = data
		#project_db.sync()
	project_db.set_key_entry(dbkey,data)
	
def get_idd_key(image_name):
	return strip_file_tag(image_name)+"_DD"

def set_idd_key_entry_in_memory(image_name,key,object):
	dbkey = get_idd_key(image_name)
	
	project_db = EMProjectDB()
	try:
		data = project_db.get_key_entry_in_memory(dbkey)
	except:
#		print "there was no data in memory so i am setting a new one",image_name,key
		data = None
		
	if data == None: data = {}
	
	data[key] = object
	project_db.set_key_entry_in_memory(dbkey,data)
	

def set_idd_key_entry(image_name,key,object):
	'''
	write a key/object pair to the Image Database Dictionary associat
	'''
	dbkey = get_idd_key(image_name)
	
	project_db = EMProjectDB()
	try:
		data = project_db[dbkey]
	except:
		data = {}
		
	if data == None: data = {}
	
	data[key] = object
	#project_db[self.get_dd_key()] = data
		#project_db.sync()
	project_db.set_key_entry(dbkey,data)

def get_idd_key_entry_in_memory(image_name,key):
	
	dbkey = get_idd_key(image_name)
	project_db = EMProjectDB()
	try:
		data = project_db.get_key_entry_in_memory(dbkey)
		if data == None: return None
		else: return data[key]
	except: return None
	
	
def get_idd_key_entry(image_name,key):
	dbkey = get_idd_key(image_name)
	
	project_db = EMProjectDB()
	try:
		data = project_db[dbkey]
		if data == None: return None
		else: return data[key] 
	except: return None


class AutoBoxer:
	'''
	FIXME - ALL COMMENTS IN THIS CLASS ARE LIABLE TO BE OUTDATED
	Base class design for auto boxers
	'''
	def __init__(self):
		self.version = 1.0

	def get_template(self):
		'''This should return a single template which is an EMData object'''
		raise Exception
	
	def name(self):
		'''
		Every autoboxer should return a unique name
		'''
		raise Exception
	
	def add_reference(self,box):
		'''
		add a reference box - the box should be in the format of a Box object, see above
		Returns 0 if there is a problem, returns 1 if it's all good
		Adds a reference to a list
		'''
		raise Exception
	
	def remove_reference(self,box):
		'''
		Remove a reference box - the box should in the format of a Box object, see above
		Pops a reference from a list
		'''
		raise Exception
	
	def reference_moved(self,ref):
		'''
		If a reference was moved interactively in the interface this is the function that should be called
		'''
		raise Exception

	def get_template(self):
		'''
		Return the template that is being used. Returns None if there is not template
		'''
		raise Exception

	def set_box_size(self,box_size):
		'''
		Hard set the box_size. Note that nothing is done to the reference boxes. It is
		assumed whichever part of the program calls this function also updates the Box objects
		independently (which implicitly affects the boxes stored internally in the AutoBoxer
		class, because it only ever stores programmatic references)
		'''
		raise Exception
	
	def auto_box(self,correlation,boxes=[],exclusion=None):
		'''
		The main auto_box routine. The calling program should pass in its own correlation map (EMData), and optionally
		an exclusion map of ones and zeros (0 means include, non zero means exclude). Also a list of boxes that will
		not be removed prior to the auto_boxing (and hence probably have associated excluded regions in the exlcusion
		image, but this is beside the point), The model of use here is that
		the calling program should get the current template from the AutoBoxer to generate the correlation map. The
		calling program should be able to cache the correlation map, so should be able to detect if there's been
		a template update by asking for the current set of references (get_references) and cross checking against a list of its own.
		@Returns a list of Boxes
		'''
		raise Exception

class TrimSwarmTemplate:
	'''
	used from writing a template to the database
	'''
	def __init__(self,swarmTemplate):
		self.refboxes = []		# this will eventually be a list of Box objects
		for ref in swarmTemplate.refboxes:
			self.refboxes.append(TrimBox(ref))
		self.template = swarmTemplate.template	# an EMData object that is the template
		self.template_ts = swarmTemplate.template_ts 	# a time stamp that records when the template was generate

class SwarmTemplate:
	def __init__(self,autoboxer):
		self.refboxes = []		# this will eventually be a list of Box objects
		self.template = None	# an EMData object that is the template
		self.template_ts = -1 	# a time stamp that records when the template was generate
		self.autoboxer = autoboxer
	
	def become(self,trimSwarmTemplate):
		if not isinstance(trimSwarmTemplate,TrimSwarmTemplate):
			print "error, can't become anything other than a TrimSwarmTemplate"
			return 0
		else:
			
			self.template = trimSwarmTemplate.template	# an EMData object that is the template
			self.template_ts = trimSwarmTemplate.template_ts 	# a time stamp that records when the template was generate
			self.refboxes = []		# this will eventually be a list of Box objects
			for ref in trimSwarmTemplate.refboxes:
				b = Box()
				b.become(ref)
				self.refboxes.append(b)
			
	def get_template(self):
		return self.template
	
	def get_image(self):
		# FIXME - get_template should be replaced with this function
		return self.template
	
	def get_template_ts(self):
		return self.template.get_attr("template_time_stamp")
	
	def get_references(self):
		#print "asking template for references, there are ",len(self.refboxes)
		return self.refboxes
	
	def append_reference(self,ref):
		if isinstance(ref,Box):
			self.refboxes.append(ref)
		else:
			print "error, can't append that reference, it's not of type Box"
	
	def change_box_size(self,box_size):
		if len(self.refboxes) == 0: return
		elif box_size == self.refboxes[0].xsize: return
		else:
			for box in self.refboxes:
				box.change_box_size(box_size)
	
	def remove_reference(self,box):
		'''
		Returns 1 if the reference was removed
		Returns 0 if it wasn't found
		'''
		if not isinstance(box,Box):
			print "error, can't remove a reference that isn't a box"
			return 0
		
		for j,tmp in enumerate(self.refboxes):
			if box.xcorner == tmp.xcorner and box.ycorner == tmp.ycorner:
				t = self.refboxes.pop(j)
				return 1	
				
		return 0
	
	def gen_template(self):
		'''
		Returns 0 if there are errors
		Return 1 if not
		'''
		# you can only generate a template if there are references
		if len(self.refboxes) <= 0: 
			print 'error, cant call private function gen_template when there are no refboxes, this is an internal error'
			return 0
		
		images_copy = []
		for ref in self.refboxes:
			# some references can be excluded from the template generation procedure, this is flagged
			# by the isanchor flag
			if ref.isanchor == False:
				continue
			image = ref.get_small_box_image(self.autoboxer)
			images_copy.append(image)
		if len(images_copy) == 0:
			print 'error, you have probably set references that all have the isanchor flag set to false, which exluded them all from the template making process'
			print 'can not proceed without references to create template'
			return 0
			
		ave = images_copy[0].copy()
		
		for i in range(1,len(images_copy)):
			#ta = images_copy[i].align("rotate_translate",ave,{},"dot",{"normalize":1})
			ave.add(images_copy[i])
		
		#ave.write_image("prealigned.hdf")
		ave.mult(1.0/len(images_copy))
		ave.process_inplace("math.radialaverage")
		ave.process_inplace("xform.centeracf")
		ave.process_inplace("mask.sharp",{'outer_radius':ave.get_xsize()/2})
		
		#for image in images_copy:
		#	image.write_image("aligned_refs.img",-1)
		#
		#ave.write_image("aligned_refs.img",-1)
		#
		#black = EMData(image.get_xsize(),image.get_ysize())
		#black.to_zero()
		#black.write_image("aligned_refs.img",-1)
		
		#ave.write_image("ave.hdf")
		shrink = self.autoboxer.get_best_shrink()
		# 4 is a magic number
		for n in range(0,4):
			t = []
			for idx,i in enumerate(images_copy):
				ta = i.align("translational",ave,{},"dot",{"normalize":1})
				#t.append(ta)
				
				# FIXME - make it so that a newly clipped portion of the original image
				# is used as the 'aligned' image, to avoid zeroing effects at the edges
				# The problem with this approach is one of persistence. If the box has no associated
				# boxingobj then it will fail. The box will have no boxing obj using the persistent database 
				# approach of e2boxer
				dx = ta.get_attr("align.dx")
				dy = ta.get_attr("align.dy")
				box = self.refboxes[idx]
				size = ta.get_xsize()
				image = box.get_small_image(self.autoboxer)
				a = image.get_clip(Region(int(box.xcorner/shrink-dx),int(box.ycorner/shrink-dy),size,size))
				a.process_inplace("normalize.edgemean")
				
				t.append(a)
				
		
			ave = t[0].copy()
			for i in range(1,len(images_copy)):
				ave.add(t[i])
				
			ave.mult(1.0/len(t))
			ave.process_inplace("math.radialaverage")
			ave.process_inplace("xform.centeracf")
			# edge normalize here SL before
			ave.process_inplace("mask.sharp",{'outer_radius':ave.get_xsize()/2})
			# or normalize and no mask
		
		#debug, un-comment to see the aligned refs and the final template
		#for image in t:
		#	image.write_image("aligned_refs.img",-1)
		
		#ave.write_image("aligned_refs.img",-1)
		
		#black = EMData(image.get_xsize(),image.get_ysize())
		#black.to_zero()
		#black.write_image("aligned_refs.img",-1)
		#END uncomment block
		self.template = ave
		
		self.template_ts = gm_time_string()
		self.template.set_attr("template_time_stamp",self.template_ts)
		return 1
	
	

class TrimSwarmAutoBoxer:
	def __init__(self,swarmAutoBoxer):
			
		self.box_size = swarmAutoBoxer.box_size
		self.shrink = swarmAutoBoxer.shrink
		self.templatedimmin = swarmAutoBoxer.templatedimmin
		
		self.optthreshold = swarmAutoBoxer.optthreshold
		self.optprofile = copy(swarmAutoBoxer.optprofile)
		self.optprofileradius = swarmAutoBoxer.optprofileradius
		self.selection_mode = swarmAutoBoxer.selection_mode
		self.template_ts = swarmAutoBoxer.template_ts
		self.state_ts = swarmAutoBoxer.state_ts
		self.mode = swarmAutoBoxer.mode
		self.refupdate = swarmAutoBoxer.refupdate
		self.regressiveflag = swarmAutoBoxer.regressiveflag
		
		self.template = TrimSwarmTemplate(swarmAutoBoxer.template)
		self.creationTS = swarmAutoBoxer.creationTS
		self.convenienceString = swarmAutoBoxer.convenienceString
	
	def get_unique_stamp(self):
		return "autoboxer_" + self.creationTS
	
	def set_convenience_name(self,string):
		self.convenienceString = string
		#print "set db string",self.get_unique_stamp()
	
	def get_convenience_name(self):
		return self.convenienceString
	
	def set_creation_ts(self,TS):
		''' 
		This function is used when creating a fresh copy of this autoboxer and placing it
		in the database in GUIbox
		'''
		#print "setting creationTS",TS
		self.creationTS = TS
		
	def get_creation_ts(self):
		return self.creationTS
	
	def get_state_ts(self):
		return self.state_ts
	
	def get_template_ts(self):
		return self.template.get_template_ts()
	
	def write_to_db(self):
		data = {}
		data["autoboxer"] = self
		data["convenience_name"] = self.get_convenience_name()
		
		project_db = EMProjectDB()
		project_db.set_key_entry(self.get_unique_stamp(),data)


class SwarmAutoBoxer(AutoBoxer):
	'''
	This is an autoboxer that encapsulates the boxing approach first developed in SwarmPS
	'''
	THRESHOLD = "Threshold"
	SELECTIVE = "Selective"
	MORESELECTIVE = "More Selective"
	DYNAPIX = 1
	ANCHOREDDYNAPIX = 2
	USERDRIVEN = 3
	ANCHOREDUSERDRIVEN = 4
	COMMANDLINE = 5
	def __init__(self,parent):
		AutoBoxer.__init__(self)

		self.template = SwarmTemplate(self)	# an EMData object that is the template
		self.box_size = -1		# stores the global box_size, this is the value being used by boxer in the main interface
		self.shrink = -1
		
		# more privately stuff
		self.templatedimmin = 20  # the smallest amount the template can be shrunken to. Will attempt to get as close to as possible. This is an important part of speeding things up.
		self.optthreshold = -1	# the correlation threshold, used to as the basis of finding local maxima
		self.optprofile = []	# the optimum correlation profile used as the basis of auto selection
		self.optprofileradius = -1 # the optimum radius - used to choose which part of the optprofile is used as the basis of selection
		self.selection_mode = SwarmAutoBoxer.SELECTIVE	# the autobox method - see EMData::BoxingTools for more details
		self.cmp_mode = BoxingTools.CmpMode.SWARM_RATIO
		BoxingTools.set_mode(self.cmp_mode)
		self.__shrink = -1
		
		self.template_ts = -1 # a template time stamp to 
		self.state_ts = -1 # the time stamp that records when the current template and parameters are completely up to date
		
		
		self.mode = SwarmAutoBoxer.DYNAPIX
		self.refupdate = False # this is a flag used when self.mode is USERDRIVEN
		self.permissablemodes = [SwarmAutoBoxer.DYNAPIX,SwarmAutoBoxer.ANCHOREDDYNAPIX,SwarmAutoBoxer.USERDRIVEN,SwarmAutoBoxer.ANCHOREDUSERDRIVEN,SwarmAutoBoxer.COMMANDLINE]
		self.permissablecmp_modes = [BoxingTools.CmpMode.SWARM_RATIO,BoxingTools.CmpMode.SWARM_DIFFERENCE, BoxingTools.CmpMode.SWARM_AVERAGE_RATIO]  # the permissiable peak profile comparitor modes - for convenience when double
		self.permissableselection_modes = [SwarmAutoBoxer.THRESHOLD,SwarmAutoBoxer.SELECTIVE,SwarmAutoBoxer.MORESELECTIVE]  # the permissiable selection modes - for convenience when double checking the calling program is setting the selectionmode explicitly (through set_selection_mode )
		self.regressiveflag = False	# flags a force removal of non references in the Boxable in auto_box
		
		self.dummybox = None
		self.parent = parent
		
		self.creationTS = gm_time_string()
		self.convenienceString = ""
		self.set_convenience_name(self.get_creation_ts()) # this string is the string that users will use to name this autoboxer in the GUIboxCtrl
		
		self.image_proc_params_mediator = ImageProcParamsMediator(self)
	
	def get_params_mediator(self):
		return self.image_proc_params_mediator

	def dynapix_on(self):
		return (self.mode == SwarmAutoBoxer.DYNAPIX or self.mode == SwarmAutoBoxer.ANCHOREDDYNAPIX)
	
	def anchor_on(self):
		return (self.mode == SwarmAutoBoxer.ANCHOREDUSERDRIVEN or self.mode == SwarmAutoBoxer.ANCHOREDDYNAPIX)

	def become(self,trimSwarmAutoBoxer):			
		self.box_size = trimSwarmAutoBoxer.box_size
		self.shrink = trimSwarmAutoBoxer.shrink
		self.templatedimmin = trimSwarmAutoBoxer.templatedimmin
		
		self.optthreshold = trimSwarmAutoBoxer.optthreshold
		self.optprofile = copy(trimSwarmAutoBoxer.optprofile)
		self.optprofileradius = trimSwarmAutoBoxer.optprofileradius
		self.selection_mode = trimSwarmAutoBoxer.selection_mode
		self.template_ts = trimSwarmAutoBoxer.template_ts
		self.state_ts = trimSwarmAutoBoxer.state_ts
		self.creationTS = trimSwarmAutoBoxer.creationTS
		self.mode = trimSwarmAutoBoxer.mode
		self.refupdate = trimSwarmAutoBoxer.refupdate
		self.regressiveflag = trimSwarmAutoBoxer.regressiveflag
		self.convenienceString = trimSwarmAutoBoxer.convenienceString
		self.template = SwarmTemplate(self)
		self.template.become(trimSwarmAutoBoxer.template)
		
		# Things that only the SwarmAutoBoxer (not necessary for the TrimSwarmAutoBoxer to do this in its constructor
		self.__update_ref_params()
		
		
	def set_convenience_name(self,string):
		self.convenienceString = string
		
	def get_convenience_name(self):
		return self.convenienceString
	
	def get_boxable(self):
		return self.parent.get_boxable()
	
	def get_state_ts(self):
		return self.state_ts
	
	def get_template_ts(self):
		return self.template.get_template_ts()
	
	def get_creation_ts(self):
		return self.creationTS
	
	def get_template(self):
		return self.template
	
	def set_dummy_box(self,box):
		if not box==None and not box.isdummy:
			print "you can never set a dummy box unless the isdummy flag is true"
			return 0
		
		if box != None:
			for i,ref in enumerate(self.get_ref_boxes()):
				if ref.isdummy:
					self.remove_reference(ref)
				break
			self.add_reference(box)
		else:
			if self.dummybox != None:
				self.remove_reference(self.dummybox)
			
		self.dummybox = box
	
	def set_mode_explicit(self,mode):
		if mode in self.permissablemodes:
			self.mode = mode
		else:
			print "error, that mode:", mode, "was not in the list of permissable modes"
			exit(1)
	def set_mode(self,dynapix,anchortemplate):
		if dynapix:
			if anchortemplate: self.mode = SwarmAutoBoxer.ANCHOREDDYNAPIX
			else: self.mode = SwarmAutoBoxer.DYNAPIX
		else:
			if anchortemplate: self.mode = SwarmAutoBoxer.ANCHOREDUSERDRIVEN
			else: self.mode = SwarmAutoBoxer.USERDRIVEN

	def set_cmp_mode(self,cmp_mode):
		if cmp_mode in self.permissablecmp_modes:
			if self.cmp_mode != cmp_mode:
				self.cmp_mode = cmp_mode
				BoxingTools.set_mode(self.cmp_mode)
				if not self.__full_update(): return 0
				self.regressiveflag = True
				if self.mode == SwarmAutoBoxer.DYNAPIX or self.mode == SwarmAutoBoxer.ANCHOREDDYNAPIX:
					self.auto_box(self.get_boxable())
				elif self.mode == SwarmAutoBoxer.COMMANDLINE:
					print "warning, haven't double check SwarmAutoBoxer.COMMANDLINE scenario in set_cmp_mode"
				return 1
			else:
				print "warning, attempted to set the cmp_mode to that which was already stored, no action taken"
				return 0
		else:
			print "the peak profile comparitor mode you specified:", cmp_mode,"was not recognized, no action was taken"
			return 0
	def get_mode(self):
		return self.mode
	
	def set_selection_mode(self,selection_mode):
		if selection_mode in self.permissableselection_modes:
			if self.selection_mode != selection_mode:
				self.selection_mode = selection_mode
				self.__plot_update()
				self.state_ts = gm_time_string()
				self.regressiveflag = True
				if self.mode == SwarmAutoBoxer.DYNAPIX or self.mode == SwarmAutoBoxer.ANCHOREDDYNAPIX:
					self.auto_box(self.get_boxable())
				elif self.mode == SwarmAutoBoxer.COMMANDLINE:
					print "warning, haven't double check SwarmAutoBoxer.COMMANDLINE scenario in set_selection_mode"
				return 1
			else:
				print "warning, attempted to set the selection_mode to that which was already stored, no action taken"
				return 0
		else:
			print "the selection mode you specified:", selection_mode,"was not recognized, no action was taken"
			return 0

	def name(self):
		return 'swarmautoboxer'

	def add_reference(self,boxes):
		'''
		 add a reference box - the box should be in the format of a Box, see above):
		'''
		isanchor = False
		try:
			for box in boxes:
				self.template.append_reference(box)
				if box.isanchor: isanchor = True
		except:
			box = boxes	
			isanchor = box.isanchor
			
		
		if isinstance(box,Box):
			if box.xsize != box.ysize:
				print 'error, support for uneven box dimensions is not currently implemented'
				return 0
		
			# store the box_size if we don't have one already
			if self.box_size == -1 and not box.isdummy:
				self.box_size = box.xsize
			# do a sanity check, this shouldn't happen if the program is managing everything carefully
			elif self.box_size != box.xsize and not box.isdummy:
				print 'error, the currently stored box size does not match the box_size of the reference that was just added'
				return 0
			
			self.template.append_reference(box)
			# update the data base
			self.write_specific_references_to_db(self.get_boxable().get_image_name())
			
			if self.mode == SwarmAutoBoxer.DYNAPIX:
				if not box.isanchor and not box.isdummy:
					print 'the box flag is internally inconsistent when using pure dynapix'
					return 0
				if not self.__full_update() : return 0
				self.auto_box(self.get_boxable())
			elif self.mode == SwarmAutoBoxer.ANCHOREDDYNAPIX:
				if box.isanchor and not box.isdummy:
					print 'the box flag is internally inconsistent when anchoring'
					return 0
				box.update_params(self)
				if not self.__accrue_opt_params() :
					self.state_ts = gm_time_string()
					print "there is a problem with the references"
					return 0
				self.state_ts = gm_time_string()
				self.auto_box(self.get_boxable())
			elif self.mode == SwarmAutoBoxer.USERDRIVEN:
				self.refupdate = True
				self.state_ts = -1
				self.template_ts = -1
			elif self.mode == SwarmAutoBoxer.ANCHOREDUSERDRIVEN:
				box.update_params(self)
				if not self.__accrue_opt_params() :
					self.state_ts = gm_time_string()
					print "there is a problem with the references"
					return 0
				self.state_ts = gm_time_string()
			else:
				print 'error, unknown mode in SwarmAutoBoxer'
				return 0
		
			return 1
	
			
		else:
			print "error, you cannot add a reference to the AutoBoxer if it is not in the format of a Box object"
			return 0
	
	
	def remove_reference(self,boxes):
		'''
		Should potentially be called remove_references
		This is somewhat of a hack - the function should be called remove_references because
		it works if the calling function supplies a list of boxes...
		'''
		
		try:
			isanchor = False
			for box in boxes:
				self.template.remove_reference(box)
				if box.isanchor: isanchor = True
		except:
			# there is only one box
			box = boxes
			self.template.remove_reference(box)
			isanchor = box.isanchor
		
		# update the data base
		self.write_specific_references_to_db(self.get_boxable().get_image_name())
			
		if len(self.template.refboxes) == 0:
			self.__reset()
			return 2
			
		if self.mode == SwarmAutoBoxer.DYNAPIX or self.mode == SwarmAutoBoxer.ANCHOREDDYNAPIX:
			if isanchor:
				if not self.__full_update(): return 0
				self.regressiveflag = True
				self.auto_box(self.get_boxable())
			else:
				if not self.__accrue_opt_params() :
					self.state_ts = gm_time_string()
					print "there is a problem with the references"
					return 0
				self.state_ts = gm_time_string()
				self.regressiveflag = True
				self.auto_box(self.get_boxable())
			return 1
		elif self.mode == SwarmAutoBoxer.USERDRIVEN or self.mode == SwarmAutoBoxer.ANCHOREDUSERDRIVEN:
			if isanchor:
				self.refupdate = True
				self.state_ts = -1
				self.template_ts = -1
			else:
				box.update_params(self)
				if not self.__accrue_opt_params() :
					self.state_ts = gm_time_string()
					print "there is a problem with the references"
					return 0
				self.state_ts = gm_time_string()
				self.regressiveflag = True
				
			return 1
		
		return 0
	
	def get_template_radius(self):
		'''
		Returns what will be or is the template radius
		'''
		return int(self.box_size/2/self.get_best_shrink())
	
	def reference_moved(self,box):
		'''
		If a reference was moved interactively in the interface this is the function that should be called
		The return value is whether or not autoboxing occured and hence whether or not display should be updated
		A -1 is returned if an error occured
		'''
		if self.mode == SwarmAutoBoxer.DYNAPIX or self.mode == SwarmAutoBoxer.ANCHOREDDYNAPIX:
			if box.isanchor:
				if not self.__full_update() : return 0
				self.regressiveflag = True
				self.auto_box(self.get_boxable())
			else:
				box.update_params(self)
				if not self.__accrue_opt_params() :
					self.state_ts = gm_time_string()
					print "there is a problem with the references"
					return 0
				self.state_ts = gm_time_string()
				self.regressiveflag = True
				self.auto_box(self.get_boxable())
			
			# update the data base
			self.write_specific_references_to_db(self.get_boxable().get_image_name())
			return 1
		elif self.mode == SwarmAutoBoxer.USERDRIVEN or self.mode == SwarmAutoBoxer.ANCHOREDUSERDRIVEN:
			if box.isanchor:
				self.refupdate = True
				self.state_ts = -1
				self.template_ts = -1
			else:
				box.update_params(self)
				if not self.__accrue_opt_params() :
					self.state_ts = gm_time_string()
					print "there is a problem with the references"
					return 0
				self.state_ts = gm_time_string()
				self.regressiveflag = True
				
			# update the data base
			self.write_specific_references_to_db(self.get_boxable().get_image_name())
			return 0
		else:
			print 'error, unknown mode in SwarmAutoBoxer'
			return -1
		
	def get_template_object(self):
		return self.template
		
	def get_template(self):
		if self.refupdate:
			if not self.__full_update(): return None
			self.refupdate = False
			
		if self.template == None:
			print 'error, you have either asked for the template without setting a reference, or you have added a reference and not set the refupdate flag'
			return None
		
		return self.template
		
	def get_box_size(self):
		return self.box_size
		
	def set_box_size(self,box_size,image_names):
		if (box_size < 6 ):
			print 'error, a hard limit of 6 for the box size is currently enforced. Email developers if this is a problem'
			return
		if self.box_size == box_size:	return
		
		# FIXME - how should we deal with this?
		# update the data base
		#self.write_specific_references_to_db(self.boxable.get_image_name())
		
		self.box_size = box_size
		self.get_best_shrink(True)
		
		# changing the box size of all assocated Boxables should change the box size of the references
		self.template.change_box_size(box_size)
		
		project_db = EMProjectDB()
		for image_name in image_names:
			found = False
			try:
				data = project_db[get_idd_key(image_name)]
				#trim_autoboxer = project_db[data["autoboxer_unique_id"]]
				if data != None:found = True
			except: pass
			
			
			if found:
				try:
					if data["autoboxer_unique_id"] == self.get_unique_stamp():
						#print "we have an associated image"
						boxable = Boxable(image_name,None,self)
						boxable.change_box_size(box_size)
						boxable.get_exclusion_image(True)
				except: pass
			#else:
				#print "image is not associated",image_name
			
		# make sure the shrink value is updated - use the force flag to do it
		
		if self.mode == SwarmAutoBoxer.DYNAPIX or self.mode == SwarmAutoBoxer.ANCHOREDDYNAPIX:
			if not self.__full_update(): 
				print "box size change failed, can't full update"
				return
			self.auto_box(self.get_boxable())
		elif self.mode == SwarmAutoBoxer.USERDRIVEN or self.mode == SwarmAutoBoxer.ANCHOREDUSERDRIVEN :
			self.refupdate = True
			self.state_ts = -1
			self.template_ts = -1
		else:
			print 'error, unknown mode in SwarmAutoBoxer'
	
	def get_search_radius(self):
		return int(0.75*(self.box_size)/float(self.get_best_shrink()))
	
	def get_constraining_radius(self):
		return int(0.5*(self.box_size)/float(self.get_best_shrink()))
	
	def get_best_shrink(self,force=True):	
		if self.box_size == -1:
			print "error - the box_size is currently -1 - I can't figure out the best value to shrink by"	
			return -1
			
		if self.shrink == -1 or force:	
			self.shrink = ceil(float(self.box_size)/float(self.templatedimmin))	
			
		return self.shrink
		
	def auto_box(self,boxable,updatedisplay=True,force=False):
		# If it's user driven then the user has selected a bunch of references and then hit 'autobox'.
		# In which case we do a complete reference update, which generates the template and the
		# best autoboxing parameters
		
		# this is fine - if a boxable is frozen this is more or less a flag for the autoboxer not
		# to autobox it...
		#print "in autobox"
		if boxable.is_excluded():
			print "Image is excluded, doing nothing"
			return 0
		
		if boxable.is_frozen():
			print "Image is frozen, maintaining current state"
			return 0
		# the project_db is used at various points in this function
		project_db = EMProjectDB()
		
		if len(self.get_ref_boxes()) == 0:
			boxable.clear_and_cache(True)
			# FIXME - debug/double check this functionality
			boxable.set_stamps(self.get_state_ts(),-1,self.get_unique_stamp())
			boxable.write_to_db()

			self.write_to_db(True)
			
			self.parent.autoboxer_db_changed()
			#if boxable.get_autoboxer_id() != self.get_creation_ts():
				#boxable.set_autoboxer_id(self.get_unique_stamp())
				#self.parent.autoboxer_db_changed()

			return 1

		# ref update should only be toggled if we are in user driven mode
		if self.refupdate:
			if not self.__full_update(): return 0
			self.refupdate = False

		correlation = boxable.get_correlation_image(self)

		autoboxer_state_ts = boxable.autoboxer_state_ts
		# auto boxing will occur if the time stamp of the AutoBoxer is not the
		# same as the time stamp cached by the Boxable. -1 means it's the first time.
		# Auto boxing will also occur if the autoboxer id of the boxable does not match the creation time stamp of this object
		# the creation time stamp of this object is persistent, even it was recovered from the data base
		
		if autoboxer_state_ts == -1 or autoboxer_state_ts != self.state_ts or boxable.get_autoboxer_id() != self.get_unique_stamp() or force:
			
			if self.mode == SwarmAutoBoxer.DYNAPIX or self.mode == SwarmAutoBoxer.USERDRIVEN or self.regressiveflag or self.mode == SwarmAutoBoxer.COMMANDLINE:
				# we must clear all non-refs if we're using dynapix
				boxable.delete_non_refs(updatedisplay)
				self.regressiveflag = False
			
			exclusion = boxable.get_exclusion_image().copy()
			self.__paint_excluded_box_areas(exclusion,boxable.boxes)

			boxes = self.__auto_box(correlation,boxable,exclusion)
			print "autoboxed",len(boxes)

			trimboxes = []
			# if anchoring is on we need to store the auto boxes that were not removed 
			# This is because the boxable didn't remove it's autoboxed particles prior to autoboxing
			# And we need to make the database complete for persistence
			boxable.append_stored_auto_boxes(trimboxes)
			for box in boxes:
				#if box.ismanual: print "I added a manual box"
				t = TrimBox(box)
				trimboxes.append(t)
			
			boxable.set_key_entryToIDD("auto_boxes",trimboxes)
			
			boxable.set_stamps(self.get_state_ts(),self.get_template_ts(),self.get_unique_stamp())
			boxable.write_to_db()
			
			if self.mode != SwarmAutoBoxer.COMMANDLINE:
				self.write_to_db(True)
		
			# These two lines are a hack because it may not need to happen (redundancy)
			boxable.set_autoboxer_id(self.get_unique_stamp())
			if not (self.mode != SwarmAutoBoxer.COMMANDLINE or self.parent==None):
				self.parent.autoboxer_db_changed()
				#self.__plot_update() # this is inefficent because it may have already happened FIXME
				
			self.write_specific_references_to_db(boxable.get_image_name())
			#print "set boxer id",self.get_creation_ts()
			
			# finally tell the boxable to read the results from the database
			# the False argument tells the boxable not to re-read those
			# boxes stored by the call to append_stored_auto_boxes above
			boxable.get_auto_selected_from_db(False)
				
			return 1

		else: print 'no auto boxing was necessary, up-2-date' # DEBUG
	def get_unique_stamp(self):
		return "autoboxer_"+self.get_creation_ts()
	
	def write_to_db(self,writecurrent=False):
		'''
		Writes this object to the DB using its time stamp
		If writecurrent is True then this SwarmAutoBoxer is also written to the DB as the "currentautobxer",
		meaning it is the most recently added autoboxer
		'''
		project_db = EMProjectDB()
		autoboxer_db_string = self.get_unique_stamp()
		
		trimself = TrimSwarmAutoBoxer(self)
		data = {}	
		data["autoboxer"] = trimself
		data["convenience_name"] = self.get_convenience_name()
		
		project_db.set_key_entry(autoboxer_db_string,data)
		
		if writecurrent:
			project_db.set_key_entry("current_autoboxer",trimself)
		
	def write_specific_references_to_db(self,image_name):
		'''
		Writes the references originating in the image given by image_name
		to the database.
		
		This is called in auto_box 
		'''
		
		refs = self.template.get_references()
		refs_to_write = []
		for ref in refs:
			if ref.get_image_name() == image_name:
				refs_to_write.append(TrimBox(ref))

		set_idd_key_entry(image_name,"reference_boxes",refs_to_write)
	
	def __reset(self):
		#self.box_size = -1
		self.state_ts = -1
		self.template_ts = -1

	def __auto_box(self,correlation,boxable,exclusion=None):
		'''
		Does the autoboxing. Returns a list of Boxes
		'''
		if not isinstance(correlation,EMData):
			print 'error, cannot autobox, the correlation argument is not an EMData object'
			return 0
			
			#print "using opt radius",self.radius, "which has value",tmp,"shrink was",self.shrink
		if self.selection_mode == SwarmAutoBoxer.THRESHOLD:
			mode = 0
		elif self.selection_mode == SwarmAutoBoxer.SELECTIVE:
			mode = 1
		elif self.selection_mode == SwarmAutoBoxer.MORESELECTIVE:
			mode = 2
		
		shrink = self.get_best_shrink()
		# Warning, this search radius value should be the same as the one used by the BoxSets that contributed the reference boxes
		# to this AutoBoxer object. There should be one place/function in the code where both parties access this value
		searchradius = self.get_search_radius()
		soln = BoxingTools.auto_correlation_pick(correlation,self.optthreshold,searchradius,self.optprofile,exclusion,self.optprofileradius,mode)


		template = self.get_high_res_template_image()
		template.write_image("template.hdf")
		# This is what should be written to the database
		boxes = []
		
		for b in soln:
			x = b[0]
			y = b[1]
			xx = int(x*shrink)
			yy = int(y*shrink)
			box = Box(xx-self.box_size/2,yy-self.box_size/2,self.box_size,self.box_size,0)
			box.set_image_name(boxable.get_image_name())
			box.correlationscore = correlation.get(x,y)
			box.corx = b[0]
			box.cory = b[1]
			box.changed = True
#			box.correct_resolution_centering(self.get_best_shrink(),False)
			box.center(Box.CENTERPROPAGATE,template,False)
			boxes.append(box)
	
		return boxes
		
	def get_high_res_template_image(self):
		t = self.get_template_object() # this is the template object
		template = t.get_template() # this is the image
		template.write_image('template_low_res.hdf')
		template = template.copy()
		newx = self.box_size
		newy = self.box_size
		oldx = template.get_xsize()
		oldy = template.get_ysize()
		
		scale = float(newx)/float(oldx)
		
		new_center = [newx/2.0,newy/2.0]
		scale_center = [scale*oldx/2.0,scale*oldy/2.0]

		template.clip_inplace(Region((oldx-newx)/2,(oldy-newy)/2,newx,newy))
		
		template.scale(scale)
		# sometimes centers could be off... FIXME double check
		template.translate(new_center[0]-scale_center[0],new_center[1]-scale_center[1],0)
		template.process_inplace("xform.centeracf") # fixes a big problem
		return template
		
	
	def __full_update(self):
		'''
		Forces a template update, then updates all correlation images
		that the references come from, then does a parameters update.
		This is like a completeness function - this needs to happen for internal
		consistency. It is motivated by the thought that if references
		come from many images, then the correlation images of each of the images
		needs to be updated in order for the correlation parameters to be generated
		consistently (from a correlation image generated by a universal template)
		'''
		
		if not self.template.gen_template():
			print 'error, couldnt generate template'
			return 0
		
		# First tell all references' associated boxing objects to be open to the prospect 
		# if update their correlation images
		
		self.__update_ref_params()

		# parameters should be updated now
		# it's important that the BoxingObjext.update_correlation updated the parameters stored in the boxes
		if not self.__accrue_opt_params(): return 0
		
		self.state_ts = gm_time_string()

		return 1
	
	def __update_ref_params(self):
		for ref in self.get_ref_boxes():
			ref.update_params(self)
			
	def get_ref_boxes(self):
		return self.template.get_references()
	
	def __accrue_opt_params(self):
		'''
		A function for accruing the parameters of the SwarmPSAutoBoxer autoboxing technique
		returns True if optimal parameters were accrued
		return False if not
		'''

		# To determine the threshold from what we've got, iterate through all of the reference
		# boxes and use the lowest correlation score as the correlation threshold
		#print 'current params are, using a total of',len(self.refboxes),'references'
		#print 'threshod:',self.optthreshold
		#print 'profile:',self.optprofile
		#print 'optrad:',self.optprofileradius
		
		if self.dummybox == None:
			found = False
			for i,box in enumerate(self.get_ref_boxes()):
				if box.correlationscore == None:
					# this is an error which probably means that the box, as created by the user, has a strong correlation maximum next to it which is disrupting the auto parameters
					# this is mostly an error for dwoolfords attention
					# for the time being just ignoring it  probably suffices
					# FIXME
					#print "continuing on faulty"
					continue
				if found == False:
					self.optthreshold = box.correlationscore
					found = True
				else:	
					if box.correlationscore < self.optthreshold: self.optthreshold = box.correlationscore
	
			# catch the circumstance where for some strange reason things just didn't work
			# probably the user has some strange data and the rotational template isn't responding normally. 
			# correlation peaks aren't where the user thinks they are.
			if not found:
				print 'error, there were no parameter data that I could inspect. I cant make the optimal parameters'
				return False
			
			# Iterate through the reference boxes and accrue what you can think of
			# as the worst case scenario, in terms of correlation profiles
			
			
			found = False
			for i,box in enumerate(self.get_ref_boxes()):
				if box.correlationscore == None:
					##print "continuing on faulty" - this was already printed above
					continue
				
				#print i,box.optprofile
				if found == False:
					self.optprofile = copy(box.optprofile)
					n = len(self.optprofile)
					found = True
				else:
					profile = box.optprofile
					for j in range(0,n):
						if profile[j] < self.optprofile[j]: self.optprofile[j] = profile[j]
		else:
			self.optprofile = self.dummybox.optprofile
			self.optthreshold = self.dummybox.correlationscore
		
	
		# determine the point in the profile where the drop in correlation score is the greatest, store it in radius
		self.optprofileradius = -1
		tmp = self.optprofile[0]
		for i in range(1,self.get_constraining_radius()):
			# the tmp > 0 is a
			if self.optprofile[i] > tmp and tmp > 0:
				tmp = self.optprofile[i]
				self.optprofileradius = i
		
		self.__plot_update()
		#print 'NOW THEY ARE'
		#print 'threshod:',self.optthreshold
		#print 'profile:',self.optprofile			
		#print 'optrad:',self.optprofileradius
		return True
	
	def __plot_update(self):
		prof = [] # self.selfmod == SwarmAutoBoxer.THRESHOLD (nothing is the right setting in this case)
		if self.selection_mode == SwarmAutoBoxer.SELECTIVE:
			prof = [self.optprofileradius]
		elif self.selection_mode == SwarmAutoBoxer.MORESELECTIVE:
			for i in range(0,self.optprofileradius+1): prof.append(i)
		
		try:
			self.parent.opt_params_updated(self.optthreshold,self.optprofile,prof)
		except: pass
	
	def __paint_excluded_box_areas(self,exclusionimage,boxes):
	
		searchradius = self.get_search_radius()

		for box in boxes:
			# xx and yy are the centers of the image, but in real image coordinates
			xx = box.xcorner + box.xsize/2
			yy = box.ycorner + box.ysize/2
			# shrink them to the small correlation image coordinates
			xx /= self.get_best_shrink()
			yy /= self.get_best_shrink()
			# Set a positive circle into the exclusionimage
			BoxingTools.set_radial_non_zero(exclusionimage,int(xx),int(yy),searchradius)
			