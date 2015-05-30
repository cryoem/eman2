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

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """e2moviealigner.py [options] <ddd_movie_stack>
	
	Determines the optimal whole-frame alignment of a DDD movie. It can be used 
	to generate the affine transforms for each image and can perform the actual 
	alignment of the movie frames according to the transformations.
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--dark",type=str,default=None,help="Perform dark image correction using the specified image file")
	parser.add_argument("--gain",type=str,default=None,help="Perform gain image correction using the specified image file")
	parser.add_argument("--gaink2",type=str,default=None,help="Perform gain image correction. Gatan K2 gain images are the reciprocal of DDD gain images.")
	parser.add_argument("--step",type=str,default="1,1",help="Specify <first>,<step>,[last]. Processes only a subset of the input data. ie- 0,2 would process all even particles. Same step used for all input files. [last] is exclusive. Default= 1,1 (first image skipped)")
	parser.add_argument("--fixbadpixels",action="store_true",dest="fbp",default=False,help="Tries to identify bad pixels in the dark/gain reference, and fills images in with sane values instead")
	parser.add_argument("--frames",action="store_true",default=False,help="Save the dark/gain corrected frames")
	parser.add_argument("--normalize",action="store_true",default=False,help="Apply edgenormalization to input images after dark/gain")
	parser.add_argument("--simpleavg", action="store_true",help="Will save a simple average of the dark/gain corrected frames (no alignment or weighting)",default=False)
	parser.add_argument("--movie", type=int,help="Display an n-frame averaged 'movie' of the stack, specify number of frames to average",default=0)	
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
	
	if options.dark: dark = MovieModeAlignment.dark_correct(options)
	else: dark=None
	
	if options.gain: gain = MovieModeAlignment.gain_correct(options)
	elif options.gaink2: gain = EMData(options.gaink2)
	else: gain = None
	
	pid=E2init(sys.argv)
	
	for fname in args:
		if options.verbose: print "Processing", fname
		# perform background subtraction
		bgsub = MovieModeAlignment.background_subtract(options,fname,dark,gain)
		# simply average frames before alignment
		if options.simpleavg: MovieModeAlignment.simple_average(bgsub)
		# actual alignment of movie frames
		alignment = MovieModeAlignment(bgsub)
		alignment.optimize()
		alignment.write()
		# movie mode viewing
		if options.movie: alignment.show_movie(options.movie)
	
	E2end(pid)

class MovieModeAlignment:
	
	"""
	Class to hold information for optimized alignment of DDD cameras.
	"""
	
	def __init__(self, path, boxsize=512, transforms=None):
		"""
		Initialization method for MovieModeAlignment objects.
		@param path 		:	File location and name.
		@param dark 		:	Dark reference movie.
		@param gain 		:	Gain reference movie.
		@param boxsize  	:	(optional) Size of boxes used to compute average power spectra.
		@param transforms	: 	(optional) A list of Transform objects.
		""" 
		# set path and metadata parameters
		self.path = path
		self.hdr = EMData(path,0).get_attr_dict()
		#self.hdr = hdr
		if path[-4:].lower() in (".mrc"): self.hdr['nimg'] = self.hdr['nz']
		else: self.hdr['nimg'] = EMUtil.get_image_count(path)
		self.outfile = path.rsplit(".",1)[0]+"_proc.hdf"		
		# calculate regions and initialize transforms
		self._initialize_params(boxsize,transforms)
		# set incoherent and initial coherent power spectra
		self._set_ips()
		self._set_cps()
		# set initial cost to be minimized via optimization
		self._cost = sys.float_info.max
		self._optimized = False
	
	def _initialize_params(self,boxsize,transforms):
		"""
		A purely organizational private method to keep the initialization code clean and readable.
		This function takes care of initializing the variables (and allocating the 
		subsequent memory) that will be used through the lifetime of the MovieModeAligner.
		
		@param boxsize		:	Size of boxes used to compute average power spectra.
		@param transforms	:	A list of Transform objects.
		"""
		self._boxsize = boxsize
		self._regions = {}
		for i in xrange(self.hdr['nimg']):
			self._regions[i] = []
			for x in xrange(self.hdr['nx'] / boxsize - 1):
				for y in xrange(self.hdr['ny'] / boxsize - 1):
					r = Region(x*self._boxsize+self._boxsize/2,y*self._boxsize+self._boxsize/2,self._boxsize,self._boxsize)
					self._regions[i].append(r)
		self.nregions = len(self._regions)
		if not transforms:
			t = Transform({"type":"eman","tx":0.0,"ty":0.0})
			self._transforms = [t for i in xrange(self.hdr['nimg'])]
		else: self._transforms = transforms
		self.optimal_transforms = self._transforms
		self._rbox = EMData(self._boxsize,self._boxsize)
		self._cbox = EMData(self._boxsize,self._boxsize).do_fft()
		self._rboxes = EMData(self._boxsize,self._boxsize)
		self._cboxes = EMData(self._boxsize,self._boxsize).do_fft()
		self._ips = EMData(self._boxsize,self._boxsize).do_fft()
		self._cps = EMData(self._boxsize,self._boxsize).do_fft()

	def _set_ips(self):
		"""
		Private method to compute and store the 2D incoherent power spectrum.
		Regions are updated by the _update_frame_params method, which
		is called by the _update_cost method. 
		"""
		for i in xrange(self.hdr['nimg']):
			for r in self._regions[i]:
				self._rbox.read_image_c(self.path,i,False,r)
				self._rboxes += self._rbox
			self._rboxes /= self.nregions
			self._rboxes.process_inplace("normalize.edgemean")
			self._cboxes = self._rboxes.do_fft()
			self._rboxes.to_zero()
			self._ips += self._cboxes
		self._ips /= self.hdr['nimg']
		self._ips.process_inplace('math.rotationalaverage')
	
	def _set_cps(self):
		"""
		Private method to compute and store the 2D coherent power spectrum. 
		Regions are updated by the _update_frame_params method, which
		is called by the _update_cost method. Since this function
		represents our target, it is inadvisable to run this method
		more than once in a single instance.
		"""
		for i in xrange(self.hdr['nimg']):
			for r in self._regions[i]:
				self._rbox.read_image_c(self.path,i,False,r)
				self._rbox.process_inplace("normalize.edgemean")
				self._cbox = self._rbox.do_fft()
				self._cbox.ri2inten()
				self._cboxes += self._cbox
			self._cboxes /= self.nregions
			self._cps += self._cboxes
			self._cboxes.to_zero()
		self._cps /= self.hdr['nimg']
		self._cps.process_inplace('math.rotationalaverage')
	
	def _update_frame_params(self,imgnum,transform):
		"""
		Private method to update a single image according to an affine transformation. 
		Note that transformations should by be in 2D.
		
		@param transform:	An EMAN Transform object
		"""
		self._transforms[imgnum] = transform
		for region in self._regions[imgnum]:
			origin = region.get_origin()
			region.set_origin(origin + [transform['tx'],transform['ty']])
	
	def _update_cost(self, transforms):
		"""
		Private method to update the cost associated with a particular set of frame alignments.
		Our cost function is the dot product of the incoherent and coherent 2D power spectra.
		The optimizer needs to pass a list of transformations which will then be applied. If
		the alignment is improved, the transforms supplied will be stored in the 
		optimal_transforms variable for later access.
		
		@param transforms: 	List of proposed EMAN Transform objects, one for each frame in movie.
		"""
		for i,t in enumerate(transforms):
			if self._transforms[i] != t:
				self._update_frame_params(i,t)
		self._set_ips()
		cost = -1*self._ips.dot(self._cps)
		if cost < self._cost:
			self._cost = cost
			self.optimal_transforms = self._transforms
	
	def optimize(self):
		"""
		Method to perform optimization
		"""
		if self._optimized: print("Optimal alignment already determined.")
		else:
			self._optimized = True
			test = self._transforms
			self._update_cost(test)
			print('Optimizer not yet implemented. Running test with initial transforms.')
		return
	
	def write(self,name=None):
		"""
		Method to write aligned results to disk
		
		@param name:	file name to write aligned movie stack
		"""
		if not name: name=self.outfile
		if not self._optimized:
			print("Warning: Saving non-optimal alignment.\nRun the optimize method to determine best frame translations.")
		for i in xrange(self.hdr['nimg']):
			im = EMData(self.path,i)
			im.transform(self.transforms[i])
			im.write_image_c(name,i)
	
	def get_transforms(self): return self.optimal_transforms
	def get_data(self): return EMData(self.path)
	def get_header(self): return self.hdr
	
	def show_movie(self,f2avg=5):
		"""
		Method to display 'movie' of averaged frames
		
		@param f2avg: number of frames to average. Default is 5.
		"""
		mov=[]
		for i in xrange(f2avg+1,self.hdr['nimg']):
			im=sum(outim[i-f2avg-1:i])
			mov.append(im)
		display(mov)
	
	@classmethod
	def dark_correct(cls,options):
		"""
		Class method to dark correct a DDD movie stack according to a dark reference.
		
		@param options:	"argparse" options from e2ddd_movie.py
		"""
		if path[-4:].lower() in (".mrc"): nd = self.hdr['nz']
		else: nd = EMUtil.get_image_count(options.dark)
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
		return dark
	
	@classmethod
	def gain_correct(cls,options):
		"""
		Class method to gain correct a DDD movie stack.
		
		@param options:	"argparse" options from e2ddd_movie.py
		"""
		if path[-4:].lower() in (".mrc"): nd = self.hdr['nz']
		else: nd = EMUtil.get_image_count(options.gain)
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
				if dark!=None : sigg.mult(sigd)
				gain.mult(sigg)
			gain.write_image(options.gain.rsplit(".",1)[0]+"_sum.hdf")
		if dark!=None : gain.sub(dark)								# dark correct the gain-reference
		gain.mult(1.0/gain["mean"])									# normalize so gain reference on average multiplies by 1.0
		gain.process_inplace("math.reciprocal",{"zero_to":0.0})		# setting zero values to zero helps identify bad pixels
		return gain
	
	@classmethod
	def background_subtract(cls,options,path,dark,gain,outfile=None):
		"""
		Class method to background subtract and gain correct a DDD movie
		
		@param options	:	"argparse" options from e2ddd_movie.py
		@param path		:	path of movie to be background subtracted
		@param dark		:	dark reference image
		@param gain		:	gain reference image
		@param outfile	:	(optional) name of the file to be written with background subtracted movie
		"""
		if not outfile: outfile = path.rsplit(".",1)[0]+"_bgsub.hdf"
		step = options.step.split(",")
		if len(step) == 3: last = int(step[2])
		else: last = -1
		first = int(step[0])
		step  = int(step[1])
		if options.verbose : print "Range = {} - {}, Step = {}".format(first, last, step)
		for i in xrange(first,last,step):
			if path[-4:].lower() in (".mrc"):
				r = Region(0,0,i,nx,ny,1)
				im=EMData(path,0,False,r)
			else:
				im=EMData(path,i)
			if dark: im.sub(dark)
			if gain: im.mult(gain)
			im.process_inplace("threshold.clampminmax",{"minval":0,"maxval":im["mean"]+im["sigma"]*3.5,"tozero":1})
			if options.fixbadpixels: im.process_inplace("threshold.outlier.localmean",{"sigma":3.5,"fix_zero":1})		# fixes clear outliers as well as values which were exactly zero
			if options.normalize: im.process_inplace("normalize.edgemean")
			if options.frames: im.write_image(outname[:-4]+"_corr.hdf",i-first)
			im.write_image_c(outfile,i)
		return outfile
	
	@classmethod
	def simple_average(cls, path):
		"""
		Class method to compute a simple averge of all frames in a DDD movie stack.
		
		@param path	:	file location DDD movie stack
		"""
		hdr = EMData(path,0,True).get_attr_dict()
		nx = hdr['nx']
		ny = hdr['ny']
		if path[-4:].lower() in (".mrc"):
			mrc = True
			nimg = self.hdr['nz']
			r = Region(0,0,i,nx,ny,1)
		else:
			nimg = EMUtil.get_image_count(path)
		avgr=Averagers.get("mean")
		for i in xrange(nimg):
			if mrc:
				r.set_origin(0,0,i)
				avgr.add_image(EMData(path,i,False,r))
			else:
				avgr.add_image(EMData(path,i))
		av=avgr.finish()
		av.write_image(outname[:-4]+"_mean.hdf",0)

if __name__ == "__main__":
	main()