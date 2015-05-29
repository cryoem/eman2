#!/usr/bin/env python

#
# Authors: James Michael Bell & Muyuan Chen, 05/28/2015
# Copyright (c) 2015 Baylor College of Medicine
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

from EMAN2 import *
import os
import sys
import numpy as np
from e2boxinfo_rescale import boxes

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <ddd_movie_stack>
	
	Align the frames of a DDD movie stack.
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--dark",type=str,default=None,help="Perform dark image correction using the specified image file")
	parser.add_argument("--gain",type=str,default=None,help="Perform gain image correction using the specified image file")
	parser.add_argument("--gaink2",type=str,default=None,help="Perform gain image correction. Gatan K2 gain images are the reciprocal of DDD gain images.")
	parser.add_argument("--fixbadpixels",action="store_true",default=False,help="Tries to identify bad pixels in the dark/gain reference, and fills images in with sane values instead")
	parser.add_argument("--parallel", default=None, help="parallelism argument. This program supports only thread:<n>")
	parser.add_argument("--threads", default=1,type=int,help="Number of threads to run in parallel on a single computer when multi-computer parallelism isn't useful", guitype='intbox', row=24, col=2, rowspan=1, colspan=1, mode="refinement[4]")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()
	
	if len(args)<1:
		print usage
		parser.error("Specify input DDD stack")
	
	if options.parallel!=None :
		if options.parallel[:7]!="thread:":
			print "ERROR: only thread:<n> parallelism supported by this program. It is i/o limited."
			sys.exit(1)
		threads=int(options.parallel[7:])
	else: threads=1
	
	if options.threads>1: threads=max(threads,options.threads)
	if threads>1: print "Sorry, limited to one thread at the moment."

	pid=E2init(sys.argv)
	
	if options.dark : 
		nd=EMUtil.get_image_count(options.dark)
		dark=EMData(options.dark,0)
		if nd>1:
			sigd=dark.copy()
			sigd.to_zero()
			a=Averagers.get("mean",{"sigma":sigd,"ignore0":1})
			print "Summing dark"
			for i in xrange(0,nd):
				if options.verbose:
					print " {}/{}   \r".format(i+1,nd),
					sys.stdout.flush()
				t=EMData(options.dark,i)
				t.process_inplace("threshold.clampminmax",{"minval":0,"maxval":t["mean"]+t["sigma"]*3.5,"tozero":1})
				a.add_image(t)
			dark=a.finish()
			sigd.write_image(options.dark.rsplit(".",1)[0]+"_sig.hdf")
			if options.fixbadpixels:
				sigd.process_inplace("threshold.binary",{"value":sigd["sigma"]/10.0})		# Theoretically a "perfect" pixel would have zero sigma, but in reality, the opposite is true
				dark.mult(sigd)
			dark.write_image(options.dark.rsplit(".",1)[0]+"_sum.hdf")
		dark.process_inplace("threshold.clampminmax.nsigma",{"nsigma":3.0})
		dark2=dark.process("normalize.unitlen")
	else: 
		dark=None
	
	if options.gain : 
		nd=EMUtil.get_image_count(options.gain)
		gain=EMData(options.gain,0)
		if nd>1:
			sigg=gain.copy()
			sigg.to_zero()
			a=Averagers.get("mean",{"sigma":sigg,"ignore0":1})
			print "Summing gain"
			for i in xrange(0,nd):
				if options.verbose:
					print " {}/{}   \r".format(i+1,nd),
					sys.stdout.flush()
				t=EMData(options.gain,i)
				t.process_inplace("threshold.clampminmax",{"minval":0,"maxval":t["mean"]+t["sigma"]*3.5,"tozero":1})
				a.add_image(t)
			gain=a.finish()
			sigg.write_image(options.gain.rsplit(".",1)[0]+"_sig.hdf")
			if options.fixbadpixels:
				sigg.process_inplace("threshold.binary",{"value":sigg["sigma"]/10.0})		# Theoretically a "perfect" pixel would have zero sigma, but in reality, the opposite is true
				if dark!=None: 
					sigg.mult(sigd)
				gain.mult(sigg)
			gain.write_image(options.gain.rsplit(".",1)[0]+"_sum.hdf")
		if dark!=None: 
			gain.sub(dark)												# dark correct the gain-reference
		gain.mult(1.0/gain["mean"])									# normalize so gain reference on average multiplies by 1.0
		gain.process_inplace("math.reciprocal",{"zero_to":0.0})		# setting zero values to zero helps identify bad pixels
	elif options.gaink2 :
		gain=EMData(options.gaink2)
	else : gain=None
	
	for movie in args:
		if options.verbose : print "Processing", movie	
		##################################
		aligner = MovieModeAligner(movie,dark,gain)
		aligner.optimize()
		aligner.write()
		##################################
		
	E2end(pid)

class MovieModeAligner:
	
	"""Class to hold information for optimized alignment of DDD cameras."""
	
	def __init__(movie, dark, gain, boxsize=512, transform=None):
		# set path and metadata parameters
		self.path = movie
		self.hdr = EMData(movie,0,True).get_attr_dict()
		self.hdr['nimg'] = EMUtil.get_image_count(movie)
		# set references
		self.dark = dark
		self.gain = gain
		# perform background subtraction
		self._remove_background()
		# calculate regions to be averaged based on size of data
		self._initialize_regions(boxsize)
		# set initial transformations
		self._initialize_transforms(transform)
		# set incoherent and initial coherent power spectra
		self._set_ips()
		self._set_cps()
		# set initial cost to be minimized via optimization
		self._cost = np.inf
		self._optimized = False
	
	def _initialize_regions(self,boxsize):
		self.boxsize = boxsize
		self.regions = {}
		for i in xrange(self.hdr['nimg']):
			self.regions[i] = []
			for x in xrange(self.hdr['nx'] / boxsize - 1):
				for y in xrange(self.hdr['ny'] / boxsize - 1):
					r = Region(x*self.boxsize+self.boxsize/2,y*self.boxsize+self.boxsize/2,self.boxsize,self.boxsize)
					self.regions[i].append(r)
		self.nregions = len(self.regions)
	
	def _initialize_transforms(self,initial):
		if initial == None:
			initial = Transform({"type":"eman","tx":0.0,"ty":0.0})
		self._transforms = [initial for i in self.hdr['nimg']]
		self.optimal_transforms = self._transforms
	
	def _set_ips(self): # FUNCTION WITH PARAMETER TO BE OPTIMIZED 
		"""function to compute the 2D incoherent power spectrum"""
		self.ips = boxes = box = EMData(self.boxsize,self.boxsize)
		for i in xrange(self.hdr['nimg']):
			for r in self.regions:
				box.read_image_c(self.path,i,False,r)
				boxes += box
			boxes /= self.nregions
			boxes.process_inplace("normalize.edgemean")
			boxes.do_fft_inplace()
			self.ips.ri2inten()
			self.ips += boxes
		self.ips /= self.hdr['nimg']
		self.ips.process_inplace('math.rotationalaverage')
	
	def _set_cps(self): # TARGET FUNCTION FOR OPTIMIZATION
		"""function to compute the 2D coherent power spectrum"""
		self.cps = boxes = box = EMData(self.boxsize,self.boxsize)
		for i in xrange(self.hdr['nimg']):
			for r in self.regions:
				box.read_image_c(self.path,i,False,r)
				box.process_inplace("normalize.edgemean")
				box.do_fft_inplace()
				box.ri2inten()
				boxes += box
			boxes /= self.nregions
			self.cps += boxes
		self.cps /= self.hdr['nimg']
		self.cps.process_inplace('math.rotationalaverage')
	
	def _remove_background(self):
		"""function to subtract background noise from power spectra"""
		print('Background subtraction not implemented')
	
	def _update_frame_params(self,imgnum,transform):
		self._transforms[imgnum] = transform
		for region in self.regions:
			origin = region.get_origin()
			region.set_origin(origin + [transform['tx'],transform['ty']])
	
	def _update_cost(self, transforms):
		"""
		Our cost function is the dot product of the incoherent and coherent 2D power spectra
		@param transforms: list of transform objects, one for each frame in the movie
		"""
		for i,transform in enumerate(proposed_transforms):
			if self.transforms[i] != transform:
				self._update_frame_params(i,transform)
		self._set_ips()
		new_cost = -1*self.ips.dot(self.cps)
		if new_cost < self.cost:
			self._cost = -1*self.ips.dot(self.cps)
			self.optimal_transforms = self.transforms
	
	def optimize(self):
		"""Optimization of objective function"""
		if self._optimized: return
		else: self._optimized = True
		print('Not implemented yet.')
		return
	
	def write(self,name=None):
		"""Writes aligned results to disk"""
		print("Writing not yet implemented")
	
	def get_transforms(self):
		return self.transforms
	
	def get_data(self):
		return EMData(self.path)
	
	def get_header(self):
		return self.hdr

if __name__ == "__main__":
	main()