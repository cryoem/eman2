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
from Simplex import Simplex
import itertools as it
import numpy as np
import os
import sys
if 'DISPLAY' in os.environ:
	import matplotlib.pyplot as plt
else:
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """e2moviealigner.py [options] <ddd_movie_stack>
	
	Determines the optimal whole-frame alignment of a DDD movie. It can be used 
	to generate the affine transforms for each image and can perform the actual 
	alignment of the movie frames according to the transformations.
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--path",type=str,default=None,help="Specify the path to the DDD movie you wish to align")
	parser.add_argument("--dark",type=str,default=None,help="Perform dark image correction using the specified image file")
	parser.add_argument("--gain",type=str,default=None,help="Perform gain image correction using the specified image file")
	parser.add_argument("--gaink2",type=str,default=None,help="Perform gain image correction. Gatan K2 gain images are the reciprocal of DDD gain images.")
	parser.add_argument("--boxsize", type=int, help="Set the boxsize used to compute power spectra across movie frames",default=256)
	parser.add_argument("--min", type=int, help="Set the minimum translation in pixels",default=-4)
	parser.add_argument("--max", type=int, help="Set the maximum translation in pixels",default=4)
	parser.add_argument("--step",type=str,default="0,1",help="Specify <first>,<step>,[last]. Processes only a subset of the input data. ie- 0,2 would process all even particles. Same step used for all input files. [last] is exclusive. Default= 0,1 (first image skipped)")
	parser.add_argument("--fixbadpixels",action="store_true",default=False,help="Tries to identify bad pixels in the dark/gain reference, and fills images in with sane values instead")
	parser.add_argument("--frames",action="store_true",default=False,help="Save the dark/gain corrected frames")
	parser.add_argument("--normalize",action="store_true",default=False,help="Apply edgenormalization to input images after dark/gain")
	parser.add_argument("--simpleavg", action="store_true",help="Will save a simple average of the dark/gain corrected frames (no alignment or weighting)",default=False)
	parser.add_argument("--movie", type=int,help="Display an n-frame averaged 'movie' of the stack, specify number of frames to average",default=0)	
	parser.add_argument("--parallel", default=None, help="parallelism argument. This program supports only thread:<n>")
	parser.add_argument("--threads", default=1,type=int,help="Number of threads to run in parallel on a single computer when multi-computer parallelism isn't useful", guitype='intbox', row=24, col=2, rowspan=1, colspan=1, mode="refinement[4]")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()
	
	if len(args)<1 and not options.path:
		print usage
		parser.error("You must specify an input DDD stack.")
	
	if options.path: args.append(options.path) # will eventually convert to nargs=+ or nargs=*
	
	if options.parallel!=None :
		if options.parallel[:7]!="thread:":
			print "ERROR: only thread:<n> parallelism supported by this program. It is i/o limited."
			sys.exit(1)
		threads=int(options.parallel[7:])
	else: threads=1
	
	if options.threads>1: threads=max(threads,options.threads)
	if threads>1: print "Sorry, limited to one thread at the moment."
	
	if options.boxsize: bs = options.boxsize
	hdr = EMData(args[0],0,True).get_attr_dict()
	if (hdr['nx'] / bs - 1) < 3 or (hdr['ny'] / bs - 1) < 3:
		print("You will need to use a smaller box size with your data.")
		sys.exit(1)
	
	if options.min: 
		if options.min > 0: min = -options.min
		else: min = options.min
	
	if options.max:
		if options.max < 0: max = -options.max
		elif options.max < min:
			print("The parameter '--max' must be greater than --min")
			sys.exit(1)
		else: max = options.max
	
	pid=E2init(sys.argv)
	
	for fname in args:
		if options.verbose: print "Processing", fname
		options.path = fname
		if options.dark: dark = MovieModeAlignment.dark_correct(options)
		else: dark=None
		if options.gain: gain = MovieModeAlignment.gain_correct(options,dark)
		elif options.gaink2: gain = EMData(options.gaink2)
		else: gain = None
		bgsub = MovieModeAlignment.background_subtract(options,dark,gain)
		if options.simpleavg: savg = MovieModeAlignment.simple_average(bgsub,min=min,max=max)
		if options.verbose: print("Creating movie aligner object")
		alignment = MovieModeAlignment(fname,bgsub,bs)
		if options.verbose: print("Optimizing movie frame alignment")
		alignment.optimize()
		if options.verbose: print("Plotting derived alignment data")
		alignment.plot_energies()
		alignment.plot_translations()
		if options.verbose: print("Writing aligned frames to disk")
		alignment.write()
		if options.movie: alignment.show_movie(options.movie)
	
	E2end(pid)

class MovieModeAlignment:
	
	"""
	MovieModeAlignment: Class to hold information for optimized alignment of DDD cameras.
	"""
	
	def __init__(self, fname, bgsub, boxsize=256, transforms=None, min=-4, max=4):
		"""
		Initialization method for MovieModeAlignment objects.
		
		@param str fname 		:	File location original data
		@param str bgsub		:	File location of background subtracted data
		@param int boxsize  	:	(optional) Size of boxes used to compute average power spectra. Default is 512.
		@param list transforms	: 	(optional) A list of Transform objects.
		@param float min	:	(optional) The minimum alignment translation in pixels. Default is -4
		@param float max	:	(optional) The maximum alignment translation in pixels. Default is 4
		"""
		self.orig = fname
		self.hdr = EMData(fname,0).get_attr_dict()
		self.path = bgsub
		if fname[-4:].lower() in (".mrc"): self.hdr['nimg'] = self.hdr['nz']
		else: self.hdr['nimg'] = EMUtil.get_image_count(fname)
		self.outfile = fname.rsplit(".",1)[0]+"_align.hdf"
		self.dir = os.path.dirname(os.path.abspath(fname))
		self._initialize_params(boxsize,transforms,min,max)
		self._computed_objective = False
		self._calc_incoherent_power_spectrum()
		self.write_incoherent_power_spectrum()
		self._calc_coherent_power_spectrum()
		os.remove(self.dir+'/coherent.hdf') # a hack to be dealt with later
		self.write_coherent_power_spectrum()
		self._energies = [sys.float_info.max]
		self._all_energies = []
		self._optimized = False
	
	def _initialize_params(self,boxsize,transforms,min,max):
		"""
		A purely organizational private method to keep the initialization code clean and readable.
		This function takes care of initializing the variables (and allocating the 
		subsequent memory) that will be used through the lifetime of the MovieModeAligner.
		
		@param int boxsize		:	Size of boxes used to compute average power spectra.
		@param list transforms	:	A list of Transform objects.
		"""
		self._boxsize = boxsize
		self._min = min
		self._max = max
		self._regions = {}
		mx = xrange(1,self.hdr['nx'] / boxsize - 1,1)
		my = xrange(1,self.hdr['ny'] / boxsize - 1,1)
		for i in xrange(self.hdr['nimg']):
			self._regions[i] = [Region(x*boxsize+boxsize/2,y*boxsize+boxsize/2,boxsize,boxsize) for y in my for x in mx]
		self._nregions = len(self._regions)
		self._stacks = {}
		for ir in xrange(self._nregions):
			self._stacks[ir] = [self._regions[i][ir] for i in xrange(self.hdr['nimg'])]
		self._nstacks = len(self._stacks)
		if transforms: self._transforms = transforms
		else: self._transforms = np.repeat(Transform({"type":"eman","tx":0.0,"ty":0.0}),self.hdr['nimg']).tolist()
		self.optimal_transforms = self._transforms
		self._cboxes = EMData(self._boxsize,self._boxsize).do_fft()
		self._ips = EMData(self._boxsize,self._boxsize).do_fft()
		self._cps = EMData(self._boxsize,self._boxsize).do_fft()
		self._cpsflag = False
		self._ipsflag = False
	
	def _calc_incoherent_power_spectrum(self):
		"""
		Private method to compute and store the 2D incoherent power spectrum.
		This 2D function represents our objective and is only computed once
		during object initialiation.
		"""
		if self._computed_objective:
			print("Incoherent power spectrum has been computed. It is attainable via the get_incoherent_power_spectrum method.")
		else:
			self._ips.to_zero()
			for i in xrange(self.hdr['nimg']):
				img = EMData(self.path,i)
				for r in self._regions[i]: # get region
					box = img.get_clip(r)
					box.process_inplace("normalize.edgemean")
					box.do_fft_inplace() # fft region
					box.ri2inten() # convert to intensities
					self._cboxes += box # sum across regions
				self._cboxes /= self._nregions
				self._ips += self._cboxes # sum across frames
				self._cboxes.to_zero()
			self._ips /= self.hdr['nimg']
			self._ips.process_inplace('math.rotationalaverage') # smooth
			self._computed_objective = True
	
	def _calc_coherent_power_spectrum(self):
		"""
		Private method to compute and store the 2D coherent power spectrum. 
		Regions are updated by the _update_frame_params method, which
		is called by the _update_energy method.
		"""
		self._cps.to_zero()
		self._cboxes.to_zero()
		b = self._cboxes.do_ift()
		for s in xrange(self._nstacks):
			for i,r in enumerate(self._stacks[s]):
				img = EMData(self.orig,i,False,r)
				b += img
			b /= self.hdr['nimg'] # average each region across all movie frames
			b.process_inplace("normalize.edgemean")
			self._cboxes = b.do_fft() # fft
			b.to_zero()
			self._cboxes.ri2inten() # ri2inten
			self._cps += self._cboxes
			self._cboxes.to_zero()
		self._cps /= self._nregions # average
		self._cps.process_inplace('math.rotationalaverage') # smooth
	
	def _update_frame_params(self,i,t):
		"""
		Private method to update a single image according to an affine transformation. 
		Note that transformations should by be in 2D.
		
		@param int i		: 	An integer specifying which transformation to update
		@param Transform t	:	An EMAN Transform object
		"""
		self._transforms[i] = t  # update transform
		newregions = [] # update regions
		for r in self._regions[i]:
			x =  [a+b for a,b in zip(r.get_origin()[:2],t.get_trans_2d())]
			newregions.append(Region(x[0],x[1],self._boxsize,self._boxsize))
		self._regions[i] = newregions
		for ir in xrange(self._nregions): # update stacks
			self._stacks[ir][i] = self._regions[i][ir]
	
	def _update_energy(self):
		"""
		Private method to update the energy associated with a particular set of frame alignments.
		Our energy function is the dot product of the incoherent and coherent 2D power spectra.
		The optimizer needs to pass a list of transformations which will then be applied. If
		the alignment is improved, the transforms supplied will be stored in the 
		optimal_transforms variable for later access.
		
		@param list transforms	: 	List of proposed EMAN Transform objects, one for each frame in movie.
		"""
		self._calc_coherent_power_spectrum()
		self.write_coherent_power_spectrum(imnum=-1)
		energy = EMData.cmp(self._ips,'dot',self._cps,{'normalize':1})
		if energy < self._energies[-1]:
			self._energies.append(energy)
			self.optimal_transforms = self._transforms
	
	def optimize(self, epsilon = 0.001, maxiters = 250, verbose = 1):
		"""
		Method to perform optimization of movie alignment.
		
		@param float epsilon	: the learning rate for the simplex optimizer
		@param int maxiters		: the maximum number of iterations to be computed by the simplex optimizer
		@param int verbose		: 1 or 0 specifying whether (or not) simplex optimization steps will be printed
		"""
		if verbose != 0: verbose = 1
		if self._optimized: 
			print("Optimal alignment already determined.")
			return
		nm=2*self.hdr['nimg']
		guess=[np.random.randint(self._min,self._max) for i in range(nm)]
		sm=Simplex(self._compares,guess,[5]*nm,data=self)
		mn=sm.minimize(epsilon = epsilon, maxiters = maxiters, monitor = verbose) 
		print("\n\nBest Parameters: {}\nError: {}\nIterations: {}\n".format(mn[0],mn[1],mn[2]))
		self._optimized = True
	
	@staticmethod
	def _compares(vec,data):
		for vi in range(0,len(vec),2):
			t = Transform({'type':'eman','tx':vec[vi],'ty':vec[vi+1]})
			data._update_frame_params(vi/2,t)
		data._update_energy()
		return np.log(1+data.get_energy())
		
	def write(self,name=None):
		"""
		Method to write aligned results to disk
		
		@param str name: (optional) file name to write aligned movie stack
		"""
		if not name: name=self.outfile
		elif name[:-4] != '.hdf': name = name[:-4] + '.hdf'  # force HDF
		if not self._optimized: print("Warning: Saving non-optimal alignment. Run the optimize method to determine best frame translations.")
		for i in xrange(self.hdr['nimg']):
			im = EMData(self.orig,i)
			im.transform(self._transforms[i])
			im.write_image_c(name,i)

	def write_coherent_power_spectrum(self,name='coherent.hdf',imnum=None):
		"""
		Method to write coherent power spectrum to current directory.
		
		@param str name	: name of coherent power spectrum to be written to disk
		@param int num	: image slice into which the coherent power spectrum will be saved
		"""
		if name[:-4] != '.hdf': name = name[:-4] + '.hdf'  # force HDF
		rcps = self._cps.do_ift()
		if num and self._cpsflag: rcps.write_image(self.dir+'/'+name,num)
		else: rcps.write_image(self.dir+'/'+name)
		self._cpsflag = True
	
	def write_incoherent_power_spectrum(self,name='incoherent.hdf',num=None):
		"""
		Method to write incoherent power spectrum to current directory.
		
		@param str name	: name of incoherent power spectrum to be written to disk
		@param int num	: image slice into which the incoherent power spectrum will be saved
		"""
		if name[:-4] != '.hdf': name = name[:-4] + '.hdf'  # force HDF
		rips = self._ips.do_ift()
		if num and self._ipsflag: rips.write_image(self.dir+'/'+name,num)
		else: rips.write_image(self.dir+'/'+name)
		self._ipsflag = True
	
	def get_transforms(self):
		return self.optimal_transforms
	
	def get_data(self):
		return EMData(self.path)
	
	def get_header(self): 
		return self.hdr
	
	def get_energies(self):
		return self._energies[1:]

	def get_energy(self): 
		return self._energies[-1]
	
	def get_incoherent_power_spectrum(self): 
		return self._ips
	
	def get_coherent_power_spectrum(self): 
		return self._cps
	
	def show_power_spectra(self):
		print("Displaying incoherent (1) and coherent (2) power spectra.")
		display([self._ips,self._cps])

	def show_bgsub(self): 
		display([EMData(self.bgsub,i) for i in xrange(self.hdr['nimg'])])
	
	def show_orig(self): 
		display([EMData(self.orig,i) for i in xrange(self.hdr['nimg'])])
	
	def show_movie(self,f2avg=5):
		"""
		Method to display 'movie' of averaged frames
		
		@param int f2avg : 	number of frames to average. Default is 5.
		"""
		mov=[]
		for i in xrange(f2avg+1,self.hdr['nimg']):
			outim = EMData(self.path,i)
			im=sum(outim[i-f2avg-1:i])
			mov.append(im)
		display(mov)
	
	def plot_energies(self,fname=None,savefig=True,showfig=False):
		"""
		Method to plot the energies from the optimization in 1D.
		"""
		if len(self._energies) <= 1: 
			print("You must first optimize the movie alignment before running this plotting routine.")
			return
		if not fname: fname = self.dir + '/' + 'energy.png'
		fig = plt.figure()
		ax = fig.add_subplot(1,1,1)
		ax.plot(self._energies[1:])
		ax.set_title('DDD Movie Alignment: Energy')
		if savefig: plt.savefig(fname)
		if showfig: plt.show()
		return fig
	
	def plot_translations(self,fname=None,savefig=True,showfig=False):
		"""
		Method to display optimized, labeled whole-frame translation vectors in 2D.
		"""
		if not fname: fname = self.dir + '/' + 'trans.png'
		trans2d = []
		for i in xrange(self.hdr['nimg']-1):
			ti = self._transforms[i].get_trans_2d()
			tf = self._transforms[i+1].get_trans_2d()
			trans2d.append([ti[0],ti[1],tf[0]-ti[0],tf[1]-ti[1]])
		trans2d = np.array(trans2d)
		X,Y,U,V = zip(*trans2d)
		plt.figure()
		ax = plt.gca()
		ax.quiver(X,Y,U,V,angles='xy',scale_units='xy',scale=1)
		ax.set_xlim([np.min(X),np.max(U)])
		ax.set_ylim([np.min(Y),np.max(V)])
		ax.set_title('DDD Movie Alignment: Frame Motion')
		if savefig: plt.savefig(fname)
		plt.draw()
		if showfig: plt.show()
	
	@classmethod
	def dark_correct(cls,options):
		"""
		Class method to dark correct a DDD movie stack according to a dark reference.
		
		@param namespace options	:	"argparse" options from e2ddd_movie.py
		"""
		hdr = EMData(options.dark,0,True).get_attr_dict()
		if options.path[-4:].lower() in (".mrc"): nd = hdr['nz']
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
	def gain_correct(cls,options,dark):
		"""
		Class method to gain correct a DDD movie stack.
		
		@param namespace options	:	"argparse" options from e2ddd_movie.py
		"""
		hdr = EMData(options.gain,0,True).get_attr_dict()
		if options.path[-4:].lower() in (".mrc"): nd = hdr['nz']
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
	def background_subtract(cls,options,dark,gain,outfile=None):
		"""
		Class method to background subtract and gain correct a DDD movie
		
		@param namespace options	:	"argparse" options from e2ddd_movie.py
		@param str path			:	path of movie to be background subtracted
		@param str dark			:	dark reference image
		@param str gain			: 	gain reference image
		@param str outfile		:	(optional) name of the file to be written with background subtracted movie
		"""
		hdr = EMData(options.path,0,True).get_attr_dict()
		if options.path[-4:].lower() in (".mrc"): nd = hdr['nz']
		else: nd = EMUtil.get_image_count(options.path)
		if not outfile: outfile = options.path.rsplit(".",1)[0]+"_bgsub.hdf"
		step = options.step.split(",")
		if len(step) == 3: last = int(step[2])
		else: last = nd
		first = int(step[0])
		step  = int(step[1])
		if options.verbose > 5: print("Range = {} - {}, Step = {}".format(first, last, step))
		for i in xrange(first,last,step):
			if options.path[-4:].lower() in (".mrc"):
				r = Region(0,0,i,nx,ny,1)
				im=EMData(path,0,False,r)
			else: im=EMData(options.path,i)
			if dark: im.sub(dark)
			if gain: im.mult(gain)
			im.process_inplace("threshold.clampminmax",{"minval":0,"maxval":im["mean"]+im["sigma"]*3.5,"tozero":1})
			if options.fixbadpixels: im.process_inplace("threshold.outlier.localmean",{"sigma":3.5,"fix_zero":1})		# fixes clear outliers as well as values which were exactly zero
			if options.normalize: im.process_inplace("normalize.edgemean")
			if options.frames: im.write_image(outname[:-4]+"_corr.hdf",i-first)
			im.write_image(outfile,i)
		return outfile
	
	@classmethod
	def simple_average(cls, path):
		"""
		Class method to compute a simple averge of all frames in a DDD movie stack.
		
		@param str path:	file location DDD movie stack
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
		return av

if __name__ == "__main__":
	main()