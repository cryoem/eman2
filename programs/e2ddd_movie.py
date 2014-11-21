#!/usr/bin/env python

#
# Author: Steven Ludtke, 02/12/2013 (sludtke@bcm.edu)
# Copyright (c) 2000-2013 Baylor College of Medicine
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

import pprint
from EMAN2 import *
import sys
from numpy import *
import numpy.linalg as LA

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <ddd_movie_stack>

	This (still experimental) program will do various processing operations on "movies" recorded on direct detection cameras.
	It is primarily used to do whole-frame alignment of movies. This program uses a binary tree approach, first computing the
	alignment of the 1st half vs the second half, then gradually reducing the set size. 

	See e2ddd_particles for per-particle alignment.
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--align_frames", action="store_true",help="Perform whole-frame alignment of the stack",default=False)
#	parser.add_argument("--align_frames_tree", action="store_true",help="Perform whole-frame alignment of the stack hierarchically",default=False)
#	parser.add_argument("--align_frames_countmode", action="store_true",help="Perform whole-frame alignment of frames collected in counting mode",default=False)
	parser.add_argument("--save_aligned", action="store_true",help="Save dark/gain corrected and optionally aligned stack",default=False)
	parser.add_argument("--dark",type=str,default=None,help="Perform dark image correction using the specified image file")
	parser.add_argument("--gain",type=str,default=None,help="Perform gain image correction using the specified image file")
	parser.add_argument("--gaink2",type=str,default=None,help="Perform gain image correction. Gatan K2 gain images are the reciprocal of DDD gain images.")
	parser.add_argument("--step",type=str,default="1,1",help="Specify <first>,<step>,[last]. Processes only a subset of the input data. ie- 0,2 would process all even particles. Same step used for all input files. [last] is exclusive. Default= 1,1 (first image skipped)")
	parser.add_argument("--fixbadpixels",action="store_true",default=False,help="Tries to identify bad pixels in the dark/gain reference, and fills images in with sane values instead")
	parser.add_argument("--frames",action="store_true",default=False,help="Save the dark/gain corrected frames")
	parser.add_argument("--normalize",action="store_true",default=False,help="Apply edgenormalization to input images after dark/gain")
	parser.add_argument("--movie", type=int,help="Display an n-frame averaged 'movie' of the stack, specify number of frames to average",default=0)
	parser.add_argument("--simpleavg", action="store_true",help="Will save a simple average of the dark/gain corrected frames (no alignment or weighting)",default=False)
	parser.add_argument("--avgs", action="store_true",help="Testing",default=False)
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
	
	if options.threads>1 : 
		threads=max(threads,options.threads)

	#if threads>1 : print "Running with {} threads".format(threads)
	if threads>1 : print "Sorry, limited to one thread at the moment."

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
		#else: dark.mult(1.0/99.0)
		dark.process_inplace("threshold.clampminmax.nsigma",{"nsigma":3.0})
		dark2=dark.process("normalize.unitlen")
	else : dark=None
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
				#t.process_inplace("threshold.clampminmax.nsigma",{"nsigma":4.0,"tozero":1})
				t.process_inplace("threshold.clampminmax",{"minval":0,"maxval":t["mean"]+t["sigma"]*3.5,"tozero":1})
				a.add_image(t)
			gain=a.finish()
			sigg.write_image(options.gain.rsplit(".",1)[0]+"_sig.hdf")
			if options.fixbadpixels:
				sigg.process_inplace("threshold.binary",{"value":sigg["sigma"]/10.0})		# Theoretically a "perfect" pixel would have zero sigma, but in reality, the opposite is true
				if dark!=None : sigg.mult(sigd)
				gain.mult(sigg)
			gain.write_image(options.gain.rsplit(".",1)[0]+"_sum.hdf")
		#else: gain.mult(1.0/99.0)
#		gain.process_inplace("threshold.clampminmax.nsigma",{"nsigma":3.0})

		
		if dark!=None : gain.sub(dark)												# dark correct the gain-reference
		gain.mult(1.0/gain["mean"])									# normalize so gain reference on average multiplies by 1.0
		gain.process_inplace("math.reciprocal",{"zero_to":0.0})		# setting zero values to zero helps identify bad pixels
	elif options.gaink2 :
		gain=EMData(options.gaink2)
	else : gain=None

	#try: display((dark,gain,sigd,sigg))
	#except: display((dark,gain))

	step = options.step.split(",")

	if len(step) == 3 :
		last = int(step[2])
	else :
		last = -1

	first = int(step[0])
	step  = int(step[1])

	if options.verbose : print "Range = {} - {}, Step = {}".format(first, last, step)

	# the user may provide multiple movies to process at once

	for fsp in args:
		if options.verbose : print "Processing", fsp
				
		n = EMUtil.get_image_count(fsp)

		if n < 3 : 
			hdr = EMData(fsp, 0, True)

			if hdr["nz"] < 2 :
				print "ERROR: {} has only {} images. Min 3 required.".format(fsp, n)
				continue

			n = hdr["nz"]

		if last <= 0 :
			flast = n
		else :
			flast = last

		if flast > n :
			flast = n
		
		process_movie(fsp, dark, gain, first, flast, step, options)
		
	E2end(pid)

def process_movie(fsp,dark,gain,first,flast,step,options):
		outname=fsp.rsplit(".",1)[0]+"_proc.hdf"		# always output to an HDF file. Output contents vary with options

		#if fsp[-4:].lower() in (".mrc","mrcs") :
		if fsp[-4:].lower() in (".mrc") :
			hdr=EMData(fsp,0,True)			# read header
			nx,ny=hdr["nx"],hdr["ny"]
		
		# bgsub and gain correct the stack

		outim=[]
		for ii in xrange(first,flast,step):
			if options.verbose:
				print " {}/{}   \r".format(ii-first+1,flast-first+1),
				sys.stdout.flush()

			#if fsp[-4:].lower() in (".mrc","mrcs") :
			if fsp[-4:].lower() in (".mrc") :
				im=EMData(fsp,0,False,Region(0,0,ii,nx,ny,1))
			else: im=EMData(fsp,ii)
			
			if dark!=None : im.sub(dark)
			if gain!=None : im.mult(gain)
			im.process_inplace("threshold.clampminmax",{"minval":0,"maxval":im["mean"]+im["sigma"]*3.5,"tozero":1})
			if options.fixbadpixels : im.process_inplace("threshold.outlier.localmean",{"sigma":3.5,"fix_zero":1})		# fixes clear outliers as well as values which were exactly zero

			#im.process_inplace("threshold.clampminmax.nsigma",{"nsigma":3.0})
#			im.mult(-1.0)
			if options.normalize : im.process_inplace("normalize.edgemean")
			
			if options.frames : im.write_image(outname[:-4]+"_corr.hdf",ii-first)
			outim.append(im)
			#im.write_image(outname,ii-first)

		nx=outim[0]["nx"]
		ny=outim[0]["ny"]

		# show a little movie of 5 averaged frames

		if options.movie>0 :
			mov=[]
			for i in xrange(options.movie+1,len(outim)):
				im=sum(outim[i-options.movie-1:i])
				#im.write_image("movie%d.hdf"%(i/5-1),0)
				#im.process_inplace("filter.lowpass.gauss",{"cutoff_freq":.02})
				mov.append(im)
				
			display(mov)

			#mov2=[]
			#for i in xrange(0,len(outim)-10,2):
				#im=sum(outim[i+5:i+10])-sum(outim[i:i+5])
				#mov2.append(im)
				
			#display(mov2)
			
			#mov=[i.get_clip(Region(1000,500,2048,2048)) for i in mov]
			#s=sum(mov)
#			fsc=[i.calc_fourier_shell_correlation(s)[1025:2050] for i in mov]
#			plot(fsc)
		
		# A simple average

		if options.simpleavg :
			if options.verbose : print "Simple average"
			avgr=Averagers.get("mean")
			for i in xrange(len(outim)):						# only use the first second for the unweighted average
				if options.verbose:
					print " {}/{}   \r".format(i+1,len(outim)),
					sys.stdout.flush()
				avgr.add_image(outim[i])
			print ""

			av=avgr.finish()
			av.write_image(outname[:-4]+"_mean.hdf",0)

		# Generates different possibilites for resolution-weighted, but unaligned, averages

		xy=XYData()
		xy.set_size(2)
		xy.set_x(0,0)
		xy.set_y(0,1.0)
		xy.set_x(1,0.707)
		xy.set_y(1,0.0)
		if options.avgs :
			if options.verbose : print "Weighted average"
			normim=EMData(nx/2+1,ny)
			avgr=Averagers.get("weightedfourier",{"normimage":normim})
			for i in xrange(min(len(outim),25)):						# only use the first second for the unweighted average
				if options.verbose:
					print " {}/{}   \r".format(i+1,len(outim)),
					sys.stdout.flush()
				xy.set_y(1,1.0)					# no weighting
				outim[i]["avg_weight"]=xy
				avgr.add_image(outim[i])
			print ""

			av=avgr.finish()
			av.write_image(outname[:-4]+"_a.hdf",0)
#			display(normim)

			# linear weighting with shifting 0 cutoff

			xy.set_y(1,0.0)
			for i in xrange(len(outim)):
				if options.verbose:
					print " {}/{}   \r".format(i+1,len(outim)),
					sys.stdout.flush()
				xy.set_x(1,0.025+0.8*(len(outim)-i)/len(outim))
				outim[i]["avg_weight"]=xy
				avgr.add_image(outim[i])
			print ""

			av=avgr.finish()
			av.write_image(outname[:-4]+"_b.hdf",0)

			# exponential falloff with shifting width

			xy.set_size(64)
			for j in xrange(64): xy.set_x(j,0.8*j/64.0)
			for i in xrange(len(outim)):
				if options.verbose:
					print " {}/{}   \r".format(i+1,len(outim)),
					sys.stdout.flush()
				for j in xrange(64) : xy.set_y(j,exp(-j/(3.0+48.0*(len(outim)-i)/float(len(outim)))))
#				plot(xy)
				outim[i]["avg_weight"]=xy
				avgr.add_image(outim[i])
			print ""

			av=avgr.finish()
			av.write_image(outname[:-4]+"_c.hdf",0)

		# we iterate the alignment process several times

		#if options.align_frames:
			#outim2=[]
##			for im in outim: im.process_inplace("threshold.clampminmax.nsigma",{"nsigma":4,"tomean":True})
##			for im in outim: im.process_inplace("threshold.clampminmax.nsigma",{"nsigma":4})
			#av=sum(outim[-5:])
##			av=outim[-1].copy()
##			av.mult(1.0/len(outim))
			#fav=[]
			#for it in xrange(2):

				#for im in outim:
					#dx,dy,z=align(im,av,verbose=options.verbose)
					#im2=im.process("xform",{"transform":Transform({"type":"2d","tx":-dx,"ty":-dy})})
					#if options.verbose==1 : print "{}, {}".format(dx,dy)
					#outim2.append(im2)

				#print "-----"
				
				#av=sum(outim2)
				#av.mult(1.0/len(outim))
				#fav.append(av)
				#if it!=2 : outim2=[]
							
			#av.write_image(outname[:-4]+"_aliavg.hdf",0)
			#if options.save_aligned:
				#for i,im in enumerate(outim2): im.write_image(outname[:-4]+"_align.hdf",i)
			#if options.verbose>1 : display(fav,True)

		# we iterate the alignment process several times

		if options.align_frames :
			outim2 = []
			
			print len(outim)
			
			aliavg = sum(outim)				# we start with a simple average of all frames
			
			for it in xrange(3) :
				step = len(outim)		# coarsest search aligns the first 1/2 of the images against the second, step=step/2 each cycle
				xali = XYData()		# this will contain the alignments which are hierarchically estimated and improved
				yali = XYData()		# x is time in both cases, y is x or y

				while step > 1 :
					step /= 2
					i0 = 0

					while i0 < len(outim) :
						i1 = min(i0+step, len(outim))
						av0 = sum(outim[i0:i1])

						tloc=(i0+i1-1)/2.0		# the "time" of the current average
						lrange=hypot(xali.get_yatx_smooth(i1,1)-xali.get_yatx_smooth(i0,1),yali.get_yatx_smooth(i1,1)-yali.get_yatx_smooth(i0,1))*1.5
						if lrange<8 : lrange=8		
						
						guess=(xali.get_yatx_smooth(tloc,1),yali.get_yatx_smooth(tloc,1))
						if xali.get_size()>1 and guess[0]<2 and guess[1]<2 : 
							i0+=step
							continue				# if the predicted shift is too small, then we won't get it right anyway, so we just interpolate
						
##						print step,i0,xali.get_yatx_smooth(tloc,1),yali.get_yatx_smooth(tloc,1),lrange,
	#					dx,dy,Z=align_subpixel(av0,av1,guess=alignments[i1+step/2]-alignments[i0+step/2],localrange=LA.norm(alignments[i1+step-1]-alignments[i0]))

						if step==len(outim)/2 :
							dx,dy,Z=align(aliavg,av0,guess=(0,0),localrange=192,verbose=options.verbose-1)
						else:
							dx,dy,Z=align(aliavg,av0,guess=guess,localrange=lrange,verbose=options.verbose-1)
##						print dx,dy,Z			
						
						xali.insort(tloc,dx)
						yali.insort(tloc,dy)
											
						i0+=step
					
					# possible sometimes to have multiple values for the same x (img #), average in these cases

					xali.dedupx()
					yali.dedupx()
					
					# Smoothing
					# we should have all possible x-values at this point, so we just do a very simplistic smoothing

					for i in xrange(xali.get_size()-2):
						xali.set_y(i+1,(xali.get_y(i)+xali.get_y(i+1)*2.0+xali.get_y(i+2))/4.0)
						yali.set_y(i+1,(yali.get_y(i)+yali.get_y(i+1)*2.0+yali.get_y(i+2))/4.0)
					
##					print ["%6.1f"%i for i in xali.get_xlist()]
##					print ["%6.2f"%i for i in xali.get_ylist()]
##					print ["%6.2f"%i for i in yali.get_ylist()]
					
				outim2=[outim[i].get_clip(Region(-xali.get_yatx_smooth(i,1),-yali.get_yatx_smooth(i,1),outim[i]["nx"],outim[i]["ny"])) for i in xrange(len(outim))]
				
				aliavg=sum(outim2)
				aliavg.mult(1.0/len(outim2))
			
				if options.verbose>2 : 
					out=file("align%d.txt"%it,"w")
					for i in xrange(xali.get_size()):
						out.write("%1.2f\t%1.2f\n"%(xali.get_y(i),yali.get_y(i)));
					out=file("alignsm%d.txt"%it,"w")
					for i in xrange(len(outim)):
						out.write("%1.2f\t%1.2f\n"%(xali.get_yatx_smooth(i,1),yali.get_yatx_smooth(i,1)));
					xali.write_file("alignx%d.txt"%it)
					yali.write_file("aligny%d.txt"%it)

			aliavg.write_image(outname[:-4]+"_aliavg.hdf",0)
			if options.save_aligned:
				for i,im in enumerate(outim2): im.write_image(outname[:-4]+"_align.hdf",i,IMAGE_HDF, False, None, EM_USHORT)
				
			if options.verbose>2 : 
				t=sum(outim)
				t.mult(1.0/len(outim2))
				t=t.get_clip(Region(500,500,3072,3072))
				aliavg=aliavg.get_clip(Region(500,500,3072,3072))
				display([t,aliavg],True)

		# we iterate the alignment process several times

		#if options.align_frames_countmode:
			#outim2=[]
##			for im in outim: im.process_inplace("threshold.clampminmax.nsigma",{"nsigma":6,"tomean":True})	# nsigma was normally 4, but for K2 images even 6 may not be enough
			#av=sum(outim)
##			av.mult(1.0/len(outim))
			#fav=[]
			#for it in xrange(2):		# K2 data seems to converge pretty much immediately in my tests

				#for im in outim:
					#if it==0 : av.sub(im)
					#dx,dy,z=align(im,av,verbose=options.verbose)
					#im2=im.process("xform",{"transform":Transform({"type":"2d","tx":-dx,"ty":-dy})})
					#if options.verbose==1 : print "{}, {}".format(dx,dy)
					#if it==0: av.add(im2)
					#outim2.append(im2)

				#print "-----"
				
				#av=sum(outim2)
				#av.mult(1.0/len(outim))
				#fav.append(av)
				#if it!=2 : outim2=[]
							
			#av.write_image(outname[:-4]+"_aliavg.hdf",0)
			#if options.save_aligned:
				#for i,im in enumerate(outim2): im.write_image(outname[:-4]+"_align.hdf",i)
			#if options.verbose>2 : display(fav,True)

def align(s1,s2,guess=(0,0),localrange=192,verbose=0):
	"""Aligns a pair of images, and returns a (dx,dy,Z) tuple. Z is the Z-score of the best peak, not a shift.
	The search will be limited to a region of +-localrange/2 about the guess, a (dx,dy) tuple. Resulting dx,dy
	is relative to the initial guess. guess and return both indicate the shift required to bring s2 in register
	with s1"""

	# reduce region used for alignment a bit (perhaps a lot for superresolution imaging

	guess=(int(guess[0]),int(guess[1]))
	if localrange<5 : localrange=192
	newbx=good_boxsize(min(s1["nx"],s1["ny"],4096)*0.8,larger=False)
	s1a=s1.get_clip(Region((s1["nx"]-newbx)/2,(s1["ny"]-newbx)/2,newbx,newbx))
	s2a=s2.get_clip(Region((s2["nx"]-newbx)/2-guess[0],(s2["ny"]-newbx)/2-guess[1],newbx,newbx))

#	s1a.process_inplace("math.xystripefix",{"xlen":200,"ylen":200})
	s1a.process_inplace("filter.xyaxes0")
#	s1a.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.05})
#	s1a.process_inplace("threshold.compress",{"value":0,"range":s1a["sigma"]/2.0})
	s1a.process_inplace("filter.highpass.gauss",{"cutoff_abs":.002})
	
#	s2a.process_inplace("math.xystripefix",{"xlen":200,"ylen":200})
#	s2a.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.05})
	s2a.process_inplace("filter.xyaxes0")
	s2a.process_inplace("filter.highpass.gauss",{"cutoff_abs":.002})

	tot=s1a.calc_ccf(s2a)
	tot.process_inplace("xform.phaseorigin.tocenter")
	tot.process_inplace("normalize.edgemean")
	
	if verbose>2 : 
		display((s1a,s2a,tot))
#		s1a.write_image("s1a.hdf",0)
#		s2a.write_image("s2a.hdf",0)
#		tot.write_image("stot.hdf",0)
			
	if verbose>3 : display((s1a,s2a,tot),force_2d=True)
	
	dx,dy=(tot["nx"]/2-int(guess[0]),tot["ny"]/2-int(guess[1]))					# the 'false peak' should always be at the origin, ie - no translation
	mn=(tot[dx-2,dy-2]+tot[dx+2,dy+2]+tot[dx-2,dy+2]+tot[dx+2,dy-2])/4.0
#	tot[dx,dy]=mn
	for x in xrange(dx-1,dx+2):
		for y in xrange(dy-1,dy+2):
			tot[x,y]=mn		# exclude from COM
			pass

	# first pass to have a better chance at finding the first peak, using a lot of blurring

	tot2=tot.get_clip(Region(tot["nx"]/2-96,tot["ny"]/2-96,192,192))
	tot2.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.04})		# This is an empirical value. Started with 0.04 which also seemed to be blurring out high-res features.
	tot2=tot2.get_clip(Region(tot2["nx"]/2-localrange/2,tot2["ny"]/2-localrange/2,localrange,localrange))
	dx1,dy1,dz=tot2.calc_max_location()
	dx1-=localrange/2
	dy1-=localrange/2

	# second pass with less blurring to fine tune it

	tot=tot.get_clip(Region(tot["nx"]/2-12+dx1,tot["ny"]/2-12+dy1,24,24))
	tot.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.12})		# This is an empirical value. Started with 0.04 which also seemed to be blurring out high-res features.
	dx,dy,dz = tot.calc_max_location()

	dev = tot["sigma"]

	if dev == 0.0 :
		dev  = 1.0
		print "Warning: sigma is zero in 'align' in iterative step for guess (", guess[0], ",", guess[1], ")."

	zscore = tot[dx,dy] / dev		# a rough Z score for the peak
	dx -= 12
	dy -= 12
	
	tot.write_image("tot.hdf",-1)
	print 
	#while hypot(dx-tot["nx"]/2,dy-tot["ny"]/2)>64 :
		#tot[dx,dy]=0
		#dx,dy,dz=tot.calc_max_location()

	if verbose>1: print "{},{} + {},{}".format(dx1,dy1,dx,dy)
	if verbose>2: display(tot)
	
	return dx1+dx+guess[0],dy1+dy+guess[1],zscore
		
	#cl=tot.get_clip(Region(dx-8,dy-8,17,17))
	#cm=cl.calc_center_of_mass(0)
	#return dx+cm[0]-8-256,dy+cm[1]-8-256

def align_subpixel(s1,s2,guess=(0,0),localrange=192,verbose=0):
	"""Aligns a pair of images to 1/4 pixel precision, and returns a (dx,dy,Z) tuple. Z is the Z-score of the best peak, not a shift.
	The search will be limited to a region of +-localrange/2 about the guess, a (dx,dy) tuple. Resulting dx,dy
	is relative to the initial guess. guess and return both indicate the shift required to bring s2 in register
	with s1"""

	# reduce region used for alignment a bit (perhaps a lot for superresolution imaging

	guess=(int(guess[0]*2.0),int(guess[1]*2.0))
	localrange*=2
	if localrange<5 : localrange=192*2
	newbx=good_boxsize(min(s1["nx"],s1["ny"],2048)*0.8,larger=False)
	newbx*=2
	s1a=s1.get_clip(Region((s1["nx"]-newbx)/2,(s1["ny"]-newbx)/2,newbx,newbx))
	s1a.scale(2)
	s2a=s2.get_clip(Region((s2["nx"]-newbx)/2-guess[0],(s2["ny"]-newbx)/2-guess[1],newbx,newbx))
	s2a.scale(2)

#	s1a.process_inplace("math.xystripefix",{"xlen":200,"ylen":200})
	s1a.process_inplace("filter.xyaxes0")
#	s1a.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.05})
#	s1a.process_inplace("threshold.compress",{"value":0,"range":s1a["sigma"]/2.0})
	s1a.process_inplace("filter.highpass.gauss",{"cutoff_abs":.002})
	
#	s2a.process_inplace("math.xystripefix",{"xlen":200,"ylen":200})
#	s2a.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.05})
	s2a.process_inplace("filter.xyaxes0")
	s2a.process_inplace("filter.highpass.gauss",{"cutoff_abs":.002})
	
	tot=s1a.calc_ccf(s2a)
	tot.process_inplace("xform.phaseorigin.tocenter")
	tot.process_inplace("normalize.edgemean")
	
	#if verbose>1 : 
		#s1a.write_image("s1a.hdf",0)
		#s2a.write_image("s2a.hdf",0)
		#tot.write_image("stot.hdf",0)
			
	if verbose>3 : display((s1a,s2a,tot),force_2d=True)
	
	mn=tot["mean"]
	dx,dy=(tot["nx"]/2,tot["ny"]/2)					# the 'false peak' should always be at the origin, ie - no translation
	for x in xrange(dx-2,dx+3):
		for y in xrange(dy-2,dy+3):
			tot[x,y]=mn		# exclude from COM
#			pass
			
	# first pass to have a better chance at finding the first peak, using a lot of blurring

	tot2=tot.get_clip(Region(tot["nx"]/2-localrange/2,tot["ny"]/2-localrange/2,localrange,localrange))
	tot2.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.04})	# This is an empirical value. Started with 0.04 which also seemed to be blurring out high-res features.
	dx1,dy1,dz=tot2.calc_max_location()
	dx1-=localrange/2
	dy1-=localrange/2

	# second pass with less blurring to fine tune it

	tot=tot.get_clip(Region(tot["nx"]/2-24+dx1,tot["ny"]/2-24+dy1,48,48))
	tot.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.08})		# This is an empirical value. Started with 0.04 which also seemed to be blurring out high-res features.
	dx,dy,dz = tot.calc_max_location()

	dev = tot["sigma"]

	if dev == 0.0 :
		dev  = 1.0
		print "Warning: sigma is zero in 'align_subpixel' in iterative step for guess (", guess[0], ",", guess[1], ")."

	zscore = tot[dx,dy] / dev 	# a rough Z score for the peak
	dx -= 24
	dy -= 24

	#while hypot(dx-tot["nx"]/2,dy-tot["ny"]/2)>64 :
		#tot[dx,dy]=0
		#dx,dy,dz=tot.calc_max_location()

	if verbose>1: print "{},{} + {},{}".format(dx1,dy1,dx,dy)
	if verbose>2: display(tot)
	
	return (dx1+dx+guess[0])/2.0,(dy1+dy+guess[1])/2.0,zscore
		
	#cl=tot.get_clip(Region(dx-8,dy-8,17,17))
	#cm=cl.calc_center_of_mass(0)
	#return dx+cm[0]-8-256,dy+cm[1]-8-256


if __name__ == "__main__":
	main()
