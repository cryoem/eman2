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
from Anneal import SimulatedAnnealer
import numpy as np
import sys
import os
import matplotlib
if 'DISPLAY' in os.environ: # user has a display server running
	import matplotlib.pyplot as plt
else: # user is logged in via SSH without X-Forwarding (no display server)
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """e2moviealigner.py [options] <ddd_movie_stack> <ddd_movie_stack> ... <ddd_movie_stack>

	Determines the optimal whole-frame alignment of a DDD movie. It can be used
	to determine the proper x and y translations for each movie frame and can 
	perform the actual alignment according to the translations.
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--dark",type=str,default=None,help="Perform dark image correction using the specified image file")
	parser.add_argument("--gain",type=str,default=None,help="Perform gain image correction using the specified image file")
	parser.add_argument("--gaink2",type=str,default=None,help="Perform gain image correction. Gatan K2 gain images are the reciprocal of DDD gain images.")
	parser.add_argument("--boxsize", type=int, help="Set the boxsize used to compute power spectra across movie frames",default=512)
	parser.add_argument("--min", type=int, help="Set the minimum translation in pixels",default=-50.0)
	parser.add_argument("--max", type=int, help="Set the maximum translation in pixels",default=50.0)
	parser.add_argument("--step",type=str,default="0,1",help="Specify <first>,<step>,[last]. Processes only a subset of the input data. ie- 0,2 would process all even particles. Same step used for all input files. [last] is exclusive. Default= 0,1 (first image skipped)")
	parser.add_argument("--tmax", type=float, help="Set the maximum achievable temperature for simulated annealing. Default is 25000.0.",default=25000.0)
	parser.add_argument("--tmin", type=float, help="Set the minimum achievable temperature for simulated annealing. Default is 2.5.",default=2.5)
	parser.add_argument("--steps", type=int, help="Set the number of steps to run simulated annealing. Default is 50000.",default=50000)
	parser.add_argument("--updates", type=int, help="Set the number of times to update the temperature when running simulated annealing. Default is 100.",default=100)
	parser.add_argument("--nopresearch",action="store_true",default=False,help="D not run a search of the parameter space before annealing to determine 'best' max and min temperatures as well as the annealing duration.")
	parser.add_argument("--presteps",type=int,default=2000,help="The number of steps to run the exploratory pre-annealing phase.")
	parser.add_argument("--premins",type=float,default=1.0,help="The number of minutes to run each phase of the exploratory phase before annealing.")
	parser.add_argument("--maxiters", type=int, help="Set the maximum number of iterations for the simplex minimizer to run before stopping. Default is 250.",default=250)
	parser.add_argument("--epsilon", type=float, help="Set the learning rate for the simplex minimizer. Smaller is better, but takes longer. Default is 0.001.",default=0.001)
	parser.add_argument("--fixbadpixels",action="store_true",default=False,help="Tries to identify bad pixels in the dark/gain reference, and fills images in with sane values instead")
	parser.add_argument("--frames",action="store_true",default=False,help="Save the dark/gain corrected frames")
	parser.add_argument("--normalize",action="store_true",default=False,help="Apply edgenormalization to input images after dark/gain")
	parser.add_argument("--simpleavg", action="store_true",help="Will save a simple average of the dark/gain corrected frames (no alignment or weighting)",default=False)
	parser.add_argument("--compareavg",action="store_true",help="Display averaged movie frames from before and after alignment for a visual comparison.", default=False)
	parser.add_argument("--movie", type=int,help="Display an n-frame averaged 'movie' of the stack, specify number of frames to average",default=0)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()

	if len(args)<1:
		print usage
		parser.error("You must specify at least one input DDD stack.")

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

		if options.dark:
			if options.verbose: print("Performing dark reference correction")
			dark = MovieModeAligner.dark_correct(options)
		else: dark=None

		if options.gain:
			if options.verbose: print("Performing gain reference correction")
			gain = MovieModeAligner.gain_correct(options,dark)
		elif options.gaink2:
			if options.verbose: print("Performing K2 gain reference correction")
			gain = EMData(options.gaink2)
		else: gain = None

		if options.gain or options.dark:
			if options.verbose: print("Performing background subtraction")
			bgsub = MovieModeAligner.background_subtract(options,dark,gain)
		else:
			if options.verbose: print("No dark or gain references supplied. Skipping background subtraction.")
			bgsub = None

		if options.simpleavg or options.compareavg:
			if options.verbose: print("Generating unaligned average")
			if bgsub: savg = MovieModeAligner.simple_average(bgsub)
			else: savg = MovieModeAligner.simple_average(fname)

		if options.verbose: print("Creating movie aligner object")
		alignment = MovieModeAligner(fname,bgsub,bs,options.min,options.max)

		if options.verbose: print("Optimizing movie frame alignment")
		alignment.optimize(options)

		if options.verbose: print("Plotting derived alignment data")
		alignment.plot_energies(fname=fname[:-4]+'_energies.png')
		alignment.plot_translations(fname=fname[:-4]+'_translations.png')

		if options.verbose: print("Writing aligned frames to disk")
		alignment.write()

		if options.simpleavg or options.compareavg:
			if options.verbose: print("Generating simple aligned average")
			aligned_movie_file = alignment.get_aligned_filename()
			savg2 = MovieModeAligner.simple_average(aligned_movie_file)

		if options.compareavg and len(args) < 2:
			if options.verbose: print("Showing averaged movie frames before and after optimized alignment.")
			display([savg,savg2])

		if options.movie and len(args) < 2:
			if options.verbose: print("Displaying aligned movie")
			alignment.show_movie(options.movie)

	E2end(pid)

class MovieModeAligner:

	"""
	MovieModeAligner: Class to hold information for optimized alignment of DDD cameras.
	"""

	def __init__(self, fname, bgsub=None, boxsize=512, min_shift=-4, max_shift=4, transforms=None):
		"""
		Initialization method for MovieModeAligner objects.

		@param str fname 		:	File location original data
		@param str bgsub		:	File location of background subtracted data
		@param int boxsize  	:	(optional) Size of boxes used to compute average power spectra. Default is 512.
		@param list transforms	: 	(optional) A list of Transform objects.
		@param float min	:	(optional) The minimum alignment translation in pixels. Default is -4
		@param float max	:	(optional) The maximum alignment translation in pixels. Default is 4
		"""
		self.orig = fname
		self.hdr = EMData(fname,0).get_attr_dict()
		if bgsub == None: self.path = fname
		else: self.path = bgsub
		if fname[-4:].lower() in (".mrc"): self.hdr['nimg'] = self.hdr['nz']
		else: self.hdr['nimg'] = EMUtil.get_image_count(fname)
		self.outfile = fname.rsplit(".",1)[0]+"_align.hdf"
		self._dir = os.path.dirname(os.path.abspath(fname))
		self._initialize_params(boxsize,transforms,min_shift,max_shift)
		self._computed_objective = False
		if os.path.isfile(self._dir+'/'+fname[:-4]+'_incoherent.hdf'): os.remove(self._dir+'/'+fname[:-4]+'_incoherent.hdf')
		self._calc_incoherent_power_spectrum()
		self.write_incoherent_power_spectrum(name=fname[:-4]+'_incoherent.hdf')
		if os.path.isfile(self._dir+'/'+fname[:-4]+'_coherent.hdf'): os.remove(self._dir+'/'+fname[:-4]+'_coherent.hdf')
		self._calc_coherent_power_spectrum()
		self.write_coherent_power_spectrum(name=fname[:-4]+'_coherent.hdf')
		self.write_coherent_power_spectrum()
		self._energies = [sys.float_info.max]
		self._all_energies = []
		self._optimized = False

	def _initialize_params(self,boxsize,transforms,min_shift,max_shift):
		"""
		A purely organizational private method to keep the initialization code clean and readable.
		This function takes care of initializing the variables (and allocating the
		subsequent memory) that will be used through the lifetime of the MovieModeAligner.

		@param int boxsize		:	Size of boxes used to compute average power spectra.
		@param list transforms	:	A list of Transform objects.
		"""
		self._boxsize = boxsize
		self._min = min_shift
		self._max = max_shift
		self._regions = {}
		mx = xrange(1,self.hdr['nx'] / boxsize - 1,1)
		my = xrange(1,self.hdr['ny'] / boxsize - 1,1)
		for i in xrange(self.hdr['nimg']):
			self._regions[i] = [Region(x*boxsize+boxsize/2,y*boxsize+boxsize/2,boxsize,boxsize) for y in my for x in mx]
		self._nregions = len(self._regions)
		self._origins = [r.get_origin() for r in self._regions[0]]
		self._norigins = len(self._origins)
		self._stacks = {}
		for ir in xrange(self._norigins): #self._nregions):
			self._stacks[ir] = [self._regions[i][ir] for i in xrange(self.hdr['nimg'])]
		self._nstacks = len(self._stacks)
		if transforms: self._transforms = transforms
		else: self._transforms = np.repeat(Transform({"type":"eman","tx":0.0,"ty":0.0}), self.hdr['nimg']).tolist()
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
		if self._computed_objective: print("Incoherent power spectrum has been computed. It is attainable via the get_incoherent_power_spectrum method.")
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
		for ir in xrange(self._norigins): #self._nregions): # update stacks
			x = np.add(self._origins[ir][:2],t.get_trans_2d())
			newregions.append(Region(x[0],x[1],self._boxsize,self._boxsize))
		self._regions[i] = newregions
		for ir in xrange(self._norigins): #self._nregions): # update stacks
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
		energy = EMData.cmp(self._ips,'dot',self._cps,{'normalize':1})
		if energy < self._energies[-1]:
			self.write_coherent_power_spectrum(num=-1)
			self._energies.append(energy)
			self.optimal_transforms = self._transforms

	def optimize(self, options):
		"""
		Method to perform optimization of movie alignment.
		
		@param namespace options	:	"argparse" options from e2ddd_movie.py
		"""
		if self._optimized:
			print("Optimal alignment already determined.")
			return

		if options.verbose: print("Starting coarse-grained alignment")
		cs = CoarseSearch(self, tmax=options.tmax, tmin=options.tmin, steps=options.steps, updates=options.updates)
		if not options.nopresearch:
			schedule = cs.auto(options.premins,options.presteps)
			cs.set_schedule(schedule)
		state, energy = cs.anneal()
		if options.verbose: print("\n\nBest Parameters: {}\nEnergy: {}\nIterations: {}\n".format(state,energy,iters))

		if options.verbose: print("Starting fine-grained alignment")
		fs = FineSearch(self,state)
		result, error, iters = fs.minimize(options.epsilon,options.maxiters,monitor=1)

		if options.verbose: print("Initializing simplex minimizer")
		if options.verbose: print("\n\nBest Parameters: {}\nError: {}\nIterations: {}\n".format(result,error,iters))
		self._optimized = True

	@staticmethod
	def _compares(vec,aligner):
		"""
		Energy functon which also updates parameters
		
		@param list options vectors		:	translations [x1,y1,x2,y2,...] of the movie frames
		@param MovieModeAligner aligner	:	aligner object
		"""
		for vi in range(0,len(vec),2):
			t = Transform({'type':'eman','tx':vec[vi],'ty':vec[vi+1]})
			aligner._update_frame_params(vi/2,t)
		aligner._update_energy()
		return np.log(1+aligner.get_energy())

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

	def write_coherent_power_spectrum(self,name=None,num=None):
		"""
		Method to write coherent power spectrum to current directory.

		@param str name	: name of coherent power spectrum to be written to disk
		@param int num	: image slice into which the coherent power spectrum will be saved
		"""
		if not name: name = self.orig[:-4] + '_coherent.hdf'
		rcps = self._cps.do_ift()
		if num and self._cpsflag: rcps.write_image(self._dir+'/'+name,num)
		else: rcps.write_image(self._dir+'/'+name)
		self._cpsflag = True

	def write_incoherent_power_spectrum(self,name=None,num=None):
		"""
		Method to write incoherent power spectrum to current directory.

		@param str name	: name of incoherent power spectrum to be written to disk
		@param int num	: image slice into which the incoherent power spectrum will be saved
		"""
		if not name: name = self.orig[:-4] + '_coherent.hdf'
		rips = self._ips.do_ift()
		if num and self._ipsflag: rips.write_image(self._dir+'/'+name,num)
		else: rips.write_image(self._dir+'/'+name)
		self._ipsflag = True

	def get_transforms(self): return self.optimal_transforms
	def get_data(self): return EMData(self.path)
	def get_header(self): return self.hdr
	def get_aligned_filename(self): return self.outfile
	def get_energies(self): return self._energies[1:]
	def get_energy(self): return self._energies[-1]
	def get_incoherent_power_spectrum(self): return self._ips
	def get_coherent_power_spectrum(self): return self._cps

	def show_power_spectra(self): display([self._ips,self._cps])
	def show_bgsub(self): display([EMData(self.bgsub,i) for i in xrange(self.hdr['nimg'])])
	def show_orig(self): display([EMData(self.orig,i) for i in xrange(self.hdr['nimg'])])

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
		if not fname: fname = self._dir + '/' + 'energy.png'
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
		if not fname: fname = self._dir + '/' + 'trans.png'
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
		plt.draw()
		if savefig: plt.savefig(fname)
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
			if options.verbose > 2: print "Summing dark"
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
				sigd.process_inplace("threshold.binary",{"value":sigd["sigma"]/10.0})  # Theoretically a "perfect" pixel would have zero sigma, but in reality, the opposite is true
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
			if options.verbose > 2: print "Summing gain"
			for i in xrange(0,nd):
				if options.verbose:
					print " {}/{}   \r".format(i+1,nd),
					sys.stdout.flush()
				t=EMData(options.gain,i)
				t.process_inplace("threshold.clampminmax",{"minval":0,"maxval":t["mean"]+t["sigma"]*3.5,"tozero":1})
				a.add_image(t)
			gain=a.finish()
			sigg.write_image(options.gain.rsplit(".",1)[0]+"_sig.hdf")
			if options.fixbadpixels: sigg.process_inplace("threshold.binary",{"value":sigg["sigma"]/10.0})  # Theoretically a "perfect" pixel would have zero sigma, but in reality, the opposite is true
			if dark!=None:
				sigd=EMData(options.dark.rsplit(".",1)[0]+"_sig.hdf",0,False)
				sigg.mult(sigd)
			gain.mult(sigg)
			gain.write_image(options.gain.rsplit(".",1)[0]+"_sum.hdf")
		if dark!=None : gain.sub(dark)  # dark correct the gain-reference
		gain.mult(1.0/gain["mean"])  # normalize so gain reference on average multiplies by 1.0
		gain.process_inplace("math.reciprocal",{"zero_to":0.0})  # setting zero values to zero helps identify bad pixels
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
		if not outfile: outfile = options.path[:-4] + '_corr.hdf'
		#if options.verbose > 5: print("Range = {} - {}, Step = {}".format(first, last, step))
		for i in xrange(first,last,step):
			if options.path[-4:].lower() in (".mrc"):
				r = Region(0,0,i,nx,ny,1)
				im=EMData(path,0,False,r)
			else: im=EMData(options.path,i)
			if dark: im.sub(dark)
			if gain: im.mult(gain)
			im.process_inplace("threshold.clampminmax",{"minval":0,"maxval":im["mean"]+im["sigma"]*3.5,"tozero":1})
			if options.fixbadpixels: im.process_inplace("threshold.outlier.localmean",{"sigma":3.5,"fix_zero":1})  # fixes clear outliers as well as values which were exactly zero
			if options.normalize: im.process_inplace("normalize.edgemean")
			if options.frames: im.write_image(outfile,i-first)
			im.write_image(outfile,i)
		return outfile

	@classmethod
	def simple_average(cls, path, outfile=None):
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
			mrc = False
			nimg = EMUtil.get_image_count(path)
		avgr=Averagers.get("mean")
		for i in xrange(nimg):
			if mrc:
				r.set_origin(0,0,i)
				avgr.add_image(EMData(path,i,False,r))
			else: avgr.add_image(EMData(path,i))
		av=avgr.finish()
		if not outfile: outfile = path[:-4] + '_mean.hdf'
		av.write_image(outfile,0)
		return av

class FineSearch(Simplex):

	"""
	Simplex minimizer to finely search the annealed minimum yielded by CoarseSearch.
	"""

	def __init__(self, aligner, init, increments=None, kR=-1.0, kE=2.0, kC=0.5):
		self.aligner = aligner
		if increments == None: increments = [5] * len(init)
		self.kR = kR
		self.kE = kE
		self.kC = kC
		super(FineSearch, self).__init__(aligner._compares,init,increments,kR,kE,kC,aligner)

class CoarseSearch(SimulatedAnnealer):

	"""
	Simulated Annealer to coarsely search the translational alignment parameter space
	"""

	def __init__(self,aligner,state=None,tmax=25000.0,tmin=2.5,steps=50000,updates=100):
		if state == None: self.state = [s for trans in aligner._transforms for s in trans.get_trans_2d()]
		super(CoarseSearch, self).__init__(self.state)
		self.Tmax = tmax
		self.Tmin = tmin
		self.steps = steps
		self.updates = updates
		self.aligner = aligner
		self.slen = len(self.state)
		self.srange = range(0,self.slen-1,2)
		self.count = 0
		self.copy_strategy = 'slice'
		self.save_state_on_exit = False

	def move(self,scale=25.0,incr=1):
		self.state[self.count] = self.state[self.count]+scale*(2*np.random.random()-1)
		self.count = (self.count+incr) % self.slen

	def energy(self):
		for vi in self.srange:
			t = Transform({'type':'eman','tx':self.state[vi],'ty':self.state[vi+1]})
			self.aligner._update_frame_params(vi/2,t)
		self.aligner._update_energy()
		return self.aligner.get_energy()

if __name__ == "__main__":
	main()
