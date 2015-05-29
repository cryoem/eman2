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
		orig_file = movie.get_attr('source_path')
		outfile = name.rsplit(".",1)[0]+"_proc.hdf"
		
		##################################
		aligner = MovieModeAligner(movie)
		aligner.optimize()
		aligner.write(outfile)
		##################################
		
	E2end(pid)

class MovieModeAligner:
	
	"""
	Class to hold information for optimized alignment of DDD cameras.
	"""
	
	def __init__(frames, dark=None,):
		self.data = frames
		self.nimgs = frames.get_zsize()
		self.shape = [frames.get_xsize(),frames.get_ysize()]
 
 		t = Transform({"type":"eman","tx":0.0,"ty":0.0})
		self.transforms = [t for i in self.nimgs]
		self.opt_trans = self.transforms
		
		self.optim = np.inf
		self.cost = self.compute_cost()
	
	def compute_cost(self):
		
		"""
		Our cost function is the dot product of the 
		incoherent and coherent power spectra
		"""
		
		# incoherent power spectrum
		frames = EMData(self.shape[0],self.shape[1])
		for i in self.nimgs:
			r = Region(0,0,i,self.shape[0],self.shape[1],1)
			frame = self.data.get_clip(r)
			frames += frame.do_fft() # FFT each frame
		frames /= self.nimgs # take average of FFTs
		ips = frames.rotavg()
		
		# coherent power spectrum
		frames = EMData(self.shape[0],self.shape[1])
		for i in self.nimgs:
			r = Region(0,0,i,self.shape[0],self.shape[1],1)
			frame = self.data.get_clip(r)
			frame.transform(self.transforms[i])
			frames += frame
		frames /= self.nimgs # average frames
		fft_frames = frames.do_fft() # compute FFT of average
		cps = fft_frames.rotavg()
		
		return -np.dot(ips,cps) # I prefer to minimize 'cost'
	
	def optimize(self):
		"""Optimization of objective function"""
		
		print('Not implemented yet.')

		# Muyuan, you'll need to update the transforms attribute
		
		self.cost = self.compute_cost()
		
		if self.cost < self.optim:
			self.optim = self.cost
			self.opt_trans = self.tansforms
	
	def write(self,fname):
		"""Writes aligned results to disk"""
		
		for i in self.nimgs:
			r = Region(0,0,i,self.shape[0],self.shape[1],1)
			frame = frames.get_clip(r)
			frame.transform(self.transforms[i])
			frame.write_image(fname,i)
	
if __name__ == "__main__":
	main()