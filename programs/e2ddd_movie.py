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
#

import pprint
from EMAN2 import *
import sys


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <ddd_movie_stack>

	This program will take an image stack from movie mode on a DirectElectron DDD camera and process it in various ways.
	The input stack should be <dark ref> <gain ref> <img 1> ...
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--align_frames", action="store_true",help="Perform whole-frame alignment of the stack",default=False)
	parser.add_argument("--dark",type=str,default=None,help="Perform dark image correction using the specified image file")
	parser.add_argument("--gain",type=str,default=None,help="Perform gain image correction using the specified image file")
	parser.add_argument("--step",type=str,default="1,1",help="Specify <first>,<step>,[last]. Processes only a subset of the input data. ie- 0,2 would process all even particles. Same step used for all input files. [last] is exclusive. Default= 1,1 (first image skipped)")
	parser.add_argument("--movie", action="store_true",help="Display a 5-frame averaged 'movie' of the frames",default=False)
	parser.add_argument("--parallel", default=None, help="parallelism argument. This program supports only thread:<n>")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()
	
	if len(args)<1:
		print usage
		parser.error("Specify input DDD stack")

	if options.dark : 
		nd=EMUtil.get_image_count(options.dark)
		dark=EMData(options.dark,0)
		if nd>1:
			sigd=dark.copy()
			sigd.to_zero()
			a=Averagers.get("mean",{"sigma":sigd,"ignore0":1})
			print "Summing dark"
			for i in xrange(0,nd):
				t=EMData(options.dark,i)
				t.process_inplace("threshold.clampminmax",{"minval":0,"maxval":t["mean"]+t["sigma"]*3.0,"tozero":1})
				a.add_image(t)
			dark=a.finish()
			dark.write_image(options.dark.rsplit(".",1)[0]+"_sum.hdf")
			sigd.write_image(options.dark.rsplit(".",1)[0]+"_sig.hdf")
		#else: dark.mult(1.0/99.0)
		dark.process_inplace("threshold.clampminmax.nsigma",{"nsigma":2.0})
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
				t=EMData(options.gain,i)
				#t.process_inplace("threshold.clampminmax.nsigma",{"nsigma":4.0,"tozero":1})
				t.process_inplace("threshold.clampminmax",{"minval":0,"maxval":t["mean"]+t["sigma"]*3.0,"tozero":1})
				a.add_image(t)
			gain=a.finish()
			gain.write_image(options.gain.rsplit(".",1)[0]+"_sum.hdf")
			sigg.write_image(options.gain.rsplit(".",1)[0]+"_sig.hdf")
		#else: gain.mult(1.0/99.0)
		gain.process_inplace("threshold.clampminmax.nsigma",{"nsigma":2.0})
	else : gain=None
	if dark!=None and gain!=None : gain.sub(dark)												# dark correct the gain-reference
	if gain!=None : 
		gain.mult(1.0/gain["mean"])									# normalize so gain reference on average multiplies by 1.0
		gain.process_inplace("math.reciprocal",{"zero_to":1.0})		 
	
	pid=E2init(sys.argv)

	#try: display((dark,gain,sigd,sigg))
	#except: display((dark,gain))

	step=options.step.split(",")
	if len(step)==3 : last=int(step[2])
	else: last=-1
	first=int(step[0])
	step=int(step[1])
	if options.verbose: print "Range={} - {}, Step={}".format(first,last,step)

	# the user may provide multiple movies to process at once
	for fsp in args:
		if options.verbose : print "Processing ",fsp
		outname=fsp.rsplit(".",1)[0]+"_proc.hdf"		# always output to an HDF file. Output contents vary with options
		
		n=EMUtil.get_image_count(fsp)
		if n<3 : 
			print "ERROR: {} has only {} images. Min 3 required.".format(fsp,n)
			continue
		
		# bgsub and gain correct the stack
		outim=[]
		if last<=0 : flast=n
		else : flast=last
		for ii in xrange(first,flast,step):
			if options.verbose:
				print " {}/{}   \r".format(ii-first+1,flast-first+1),
				sys.stdout.flush()

			im=EMData(fsp,ii)
			
			if dark!=None : im.sub(dark)
			if gain!=None : im.mult(gain)

			#im.process_inplace("threshold.clampminmax.nsigma",{"nsigma":3.0})
			im.process_inplace("threshold.clampminmax",{"minval":0,"maxval":im["mean"]+im["sigma"]*2.0,"tozero":1})		# TODO - not sure if 2 is really safe here, even on the high end
#			im.mult(-1.0)
			
			outim.append(im)
			#im.write_image(outname,ii-first)


		# show a little movie of 5 averaged frames
		if options.movie :
			mov=[]
			for i in xrange(5,len(outim),5):
				im=sum(outim[i-5:i])
	#			im.process_inplace("normalize.edgemean")
				#im.write_image("movie%d.hdf"%(i/5-1),0)
				#im.process_inplace("filter.lowpass.gauss",{"cutoff_freq":.02})
				mov.append(im)
			display(mov)
		
		# we iterate the process several times
		outim2=[]
		av=sum(outim)
		av.mult(1.0/len(outim))
		fav=[av]
		for it in xrange(4):

			for im in outim:
				dx,dy=zonealign(im,av)
				im2=im.process("xform",{"transform":Transform({"type":"2d","tx":dx,"ty":dy})})
				print "{}, {}".format(dx,dy)
				outim2.append(im2)

			print "-----"
			
			av=sum(outim2)
			av.mult(1.0/len(outim))
			fav.append(av)
			outim2=[]
			
		display(fav)

		# show CCF between first and last frame
		#cf=mov[0].calc_ccf(mov[-1])
		#cf.process_inplace("xform.phaseorigin.tocenter")
		#display(cf)

		## save 10-frame averages without alignment
		#im=sum(outim[:10])
		#im.process_inplace("normalize.edgemean")
		#im.write_image("sum0-10.hdf",0)
		#im=sum(outim[10:20])
		#im.process_inplace("normalize.edgemean")
		#im.write_image("sum10-20.hdf",0)
		#im=sum(outim[20:30])
		#im.process_inplace("normalize.edgemean")
		#im.write_image("sum20-30.hdf",0)
	
	

		#try:
			#dot1=s1.cmp("ccc",dark2,{"negative":0})
			#dot2=s2.cmp("ccc",dark2,{"negative":0})

			##s1.sub(dark2*dot1)
			##s2.sub(dark2*dot2)

			#print dot1,dot2
		#except:
			#print "no dark"

		# alignment
		#ni=len(outim)		# number of input images in movie
		
		#s1=sum(outim[:ni/4])
		#s1.process_inplace("normalize.edgemean")
		#s2=sum(outim[ni*3/4:])
		#s2.process_inplace("normalize.edgemean")
		#dx,dy=align(s1,s2)
		#print "half vs half: ",dx,dy
		
		#dxn=dx/(ni/2.0)		# dx per n
		#dyn=dy/(ni/2.0)
		
		#s1=sum(outim[:ni/4])
		#s1.process_inplace("normalize.edgemean")
		#s2=sum(outim[ni/4:ni/2])
		#s2.process_inplace("normalize.edgemean")
		#dx,dy=align(s1,s2,(dxn*ni/4.0,dyn*ni/4.0))
		#print dx,dy
		
		#s1=sum(outim[ni/4:ni/2])
		#s1.process_inplace("normalize.edgemean")
		#s2=sum(outim[ni/2:ni*3/4])
		#s2.process_inplace("normalize.edgemean")
		#dx,dy=align(s1,s2,(dxn*ni/4.0,dyn*ni/4.0))
		#print dx,dy
		
		#s1=sum(outim[ni/2:ni*3/4])
		#s1.process_inplace("normalize.edgemean")
		#s2=sum(outim[ni*3/4:])
		#s2.process_inplace("normalize.edgemean")
		#dx,dy=align(s1,s2,(dxn*ni/4.0,dyn*ni/4.0))
		#print dx,dy
		
			
		E2end(pid)

def zonealign(s1,s2):
	s1a=s1.copy()
	s1a.process_inplace("math.xystripefix",{"xlen":200,"ylen":200})
	s1a.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.02})
	
	s2a=s2.copy()
	s2a.process_inplace("math.xystripefix",{"xlen":200,"ylen":200})
	s1a.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.02})
	
	tot=None
	for x in range(256,s1["nx"]-512,512):
		for y in range(256,s1["ny"]-512,512):
			s1b=s1a.get_clip(Region(x,y,512,512))
			s2b=s2a.get_clip(Region(x,y,512,512))
			
			c12=s1b.calc_ccf(s2b)
			c12.process_inplace("xform.phaseorigin.tocenter")
			c12.process_inplace("normalize.edgemean")
						
			cm=c12.calc_center_of_mass(0)
			try: tot.add(cm)
			except: tot=c12
			
	dx,dy=(c12["nx"]/2,c12["ny"]/2)					# the 'false peak' should always be at the origin, ie - no translation
	for x in xrange(dx-2,dx+3):
		for y in xrange(dy-2,dy+3):
			tot[x,y]=-1.0		# exclude from COM


	dx,dy,dz=tot.calc_max_location()
	while hypot(dx-256,dy-256)>50 :
		tot[dx,dy]=0
		dx,dy,dz=tot.calc_max_location()
		
		
	cl=tot.get_clip(Region(dx-8,dy-8,17,17))
	cm=cl.calc_center_of_mass(0)
	return dx+cm[0]-8-256,dy+cm[1]-8-256

def align(s1,s2,hint=None):

	s11=s1.get_clip(Region(1024,500,2048,2048))
#	s11.process_inplace("math.addsignoise",{"noise":0.5})
#	s11.process_inplace("normalize.local",{"radius":6,"threshold":0})
#	s11.process_inplace("math.xystripefix",{"xlen":200,"ylen":200})
	#s11.process_inplace("threshold.compress",{"value":s11["mean"],"range":s11["sigma"]})
	s11.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.02})
	s21=s2.get_clip(Region(1024,500,2048,2048))
#	s21.process_inplace("math.addsignoise",{"noise":0.5})
#	s21.process_inplace("normalize.local",{"radius":6,"threshold":0})
#	s21.process_inplace("math.xystripefix",{"xlen":200,"ylen":200})
	#s21.process_inplace("threshold.compress",{"value":s21["mean"],"range":s21["sigma"]})
	s21.process_inplace("filter.lowpass.gauss",{"cutoff_abs":.02})

	c12=s11.calc_ccf(s21)
	m=c12["minimum"]
#	for x in xrange(c12["nx"]): c12[x,0]=m
	c12.process_inplace("normalize.edgemean")
	c12.process_inplace("xform.phaseorigin.tocenter")

	# This peak is the false peak caused by the imperfect dark noise subtraction
	# we want to wipe this peak out
	#dx,dy,dz=tuple(c12.calc_max_location())		# dz is obviously 0
	dx,dy=(c12["nx"]/2,c12["ny"]/2)					# the 'false peak' should always be at the origin, ie - no translation
	newval=(c12[dx-3,dy]+c12[dx+3,dy]+c12[dx,dy-3]+c12[dx,dy+3])/8		# /4 would be the average, we intentionally downweight it
	for x in xrange(dx-2,dx+3):
		for y in xrange(dy-2,dy+3):
			c12[x,y]=newval
	#c12[dx-1,dy]=(c12[dx-1,dy-1]+c12[dx-1,dy+1])/2.0
	#c12[dx+1,dy]=(c12[dx+1,dy-1]+c12[dx+1,dy+1])/2.0
	#c12[dx,dy+1]=(c12[dx+1,dy+1]+c12[dx-1,dy+1])/2.0
	#c12[dx,dy-1]=(c12[dx+1,dy-1]+c12[dx-1,dy-1])/2.0
	#c12[dx,dy]=(c12[dx-1,dy]+c12[dx+1,dy]+c12[dx,dy+1]+c12[dx,dy-1])/4.0
	
	display((s11,s21,c12))
#	display(c12)
	
	if hint!=None:
		cl=c12.get_clip(Region(1024+hint[0]-4,1024+hint[1]-4,9,9))
		dx,dy,dz=cl.calc_max_location()
		cl=c12.get_clip(Region(1024+dx-3,1024+dy-3,7,7))
		cm=cl.calc_center_of_mass(0)
		return dx+cm[0]-3,dy+cm[1]-3
	else:
		dx,dy,dz=c12.calc_max_location()
		cl=c12.get_clip(Region(dx-8,dy-8,17,17))
		cm=cl.calc_center_of_mass(0)
		return dx+cm[0]-8-1024,dy+cm[1]-8-1024

if __name__ == "__main__":
	main()
