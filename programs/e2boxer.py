#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
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

# e2boxer.py  07/27/2004  Steven Ludtke
# This program is used to box out particles from micrographs/CCD frames

from EMAN2 import *
from optparse import OptionParser
from emshape import EMShape
from math import *
import time
import os
import sys

pl=()

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <image>
	
Automatic and manual particle selection. This version is specifically aimed at square boxes
for single particle analysis."""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--gui",action="store_true",help="Start the GUI for interactive boxing",default=False)
	parser.add_option("--boxsize","-B",type="int",help="Box size in pixels",default=-1)
	parser.add_option("--dbin","-D",type="string",help="Filename to read an existing box database from",default=None)
	parser.add_option("--auto","-A",type="string",action="append",help="Autobox using specified method: circle, ref, grid, pspec",default=[])
	parser.add_option("--threshold","-T",type="float",help="(auto:ref) Threshold for keeping particles. 0-4, 0 excludes all, 4 keeps all.",default=2.0)
	parser.add_option("--refptcl","-R",type="string",help="(auto:ref) A stack of reference images. Must have the same scale as the image being boxed.",default=None)
	parser.add_option("--nretest",type="int",help="(auto:ref) Number of reference images (starting with the first) to use in the final test for particle quality.",default=-1)
	parser.add_option("--retestlist",type="string",help="(auto:ref) Comma separated list of image numbers for retest cycle",default="")
	parser.add_option("--ptclsize","-P",type="int",help="(auto:circle) Approximate size (diameter) of the particle in pixels. Not required if reference particles are provided.",default=0)
	parser.add_option("--overlap",type="int",help="(auto:grid) number of pixels of overlap between boxes. May be negative.")
#	parser.add_option("--refvol","-V",type="string",help="A 3D model to use as a reference for autoboxing",default=None)
#	parser.add_option("--sym","-S",type="string",help="Symmetry of the 3D model",default=None)
	parser.add_option("--farfocus",type="string",help="filename or 'next', name of an aligned far from focus image for preliminary boxing",default=None)
	parser.add_option("--dbout",type="string",help="filename to write EMAN1 style box database file to",default=None)
	parser.add_option("--norm",action="store_true",help="Edgenormalize boxed particles",default=False)
	parser.add_option("--ptclout",type="string",help="filename to write boxed out particles to",default=None)
	parser.add_option("--savealiref",action="store_true",help="Stores intermediate aligned particle images in boxali.hdf. Mainly for debugging.",default=False)
	
	(options, args) = parser.parse_args()
	if len(args)<1 : parser.error("Input image required")

	logid=E2init(sys.argv)
	
	if (options.farfocus==None): initial=args[0]
	elif (options.farfocus=="next") : 
		initial=str(int(args[0].split('.')[0])+1)+'.'+args[0].split('.')[1]
		if not os.access(initial,os.R_OK) :
			print "No such file: ",initial
			sys.exit(1)
		print "Using initial file: ",initial
	else: initial=options.farfocus
	
	if not options.dbout and not options.ptclout : print "WARNING : No output files specified"
	
	# read the image in, though it is likely to get destroyed/modified later
	image=EMData()
	image.read_image(initial)
	image.set_attr("datatype",7)
	
	# Store this for later use, even if the image is later corrupted
	image_size=(image.get_xsize(),image.get_ysize())
	
	# read in the reference particles if provided, these could be used
	# by different autoboxing algorithms
	refptcl=None
	if options.refptcl :
		print options.refptcl
		refptcl=EMData.read_images(options.refptcl)
		refbox=refptcl[0].get_xsize()
		print "%d reference particles read (%d x %d)"%(len(refptcl),refbox,refbox)
	
	boxes=[]
	boxthr=1.0
	
	# Read box database from a file
	if options.dbin:
		# x,y,xsize,ysize,quality,changed
		boxes=[[int(j) for j in i.split()] for i in file(options.dbin,"r").readlines() if i[0]!="#"]	# this reads the whole box db file
		
		for i in boxes: 
			try: i[4]=1.0		# all qualities set to 1
			except: i.append(1.0)
			i.append(1)		# 'changed' flag initially set to 1

		if options.boxsize<5 : options.boxsize=boxes[0][2]
		else :
			for i in boxes:
				i[0]-=(options.boxsize-i[2])/2
				i[1]-=(options.boxsize-i[2])/2
				i[2]=options.boxsize
				i[3]=options.boxsize

	# we need to know how big to make the boxes. If nothing is specified, but
	# reference particles are, then we use the reference particle size
	if options.boxsize<5 :
		if options.refptcl : options.boxsize=refptcl[0].get_xsize()
		elif options.ptclsize : 
			pass
			options.boxsize=options.ptclsize
#			options.boxsize=good_boxsize(options.ptclsize*1.2)
		else : parser.error("Please specify a box size")
	else:
		if not options.boxsize in good_box_sizes:
			print "Note: EMAN2 processing would be more efficient with a boxsize of %d"%good_boxsize(options.boxsize)
	
	try: options.retestlist=[int(i) for i in options.retestlist.split(',')]
	except: options.retestlist=[]
			
	shrinkfactor=int(ceil(options.boxsize/16))
	if "pspec" in options.auto : shrinkfactor/=2
	print "Shrink factor = ",shrinkfactor
	#shrinkfactor=int(ceil(image.get_ysize()/1024.0))
	#if options.boxsize/shrinkfactor<12 : shrinkfactor/=2
	
	image.process_inplace("normalize")
	shrink=image
	shrink.mean_shrink(shrinkfactor)		# shrunken original image
	
	# This confusing line insures the shrunken images have even dimensions
	# since odd FFTs haven't been fixed yet
	if shrink.get_xsize()&1 or shrink.get_ysize()&1 :
		shrink=shrink.get_clip(Region(0,0,(shrink.get_xsize()|1)^1,(shrink.get_ysize()|1)^1))

	# now we try to clean up long range density variations in the image
	filtrad=options.boxsize*2/shrinkfactor
	flt=EMData()
	flt.set_size(shrink.get_xsize()+filtrad*4,shrink.get_ysize()+filtrad*4,1)
	flt.to_one()
	flt.process_inplace("mask.sharp",{"outer_radius":filtrad})
	flt/=(float(flt.get_attr("mean"))*flt.get_xsize()*flt.get_ysize())
	flt.process_inplace("xform.phaseorigin")
	a=shrink.get_clip(Region(-filtrad*2,-filtrad*2,shrink.get_xsize()+filtrad*4,shrink.get_ysize()+filtrad*4))
	a.process_inplace("mask.zeroedgefill")
#	a.write_image("q0.hdf",0)
#	flt.write_image("q0.hdf",1)
	a=a.convolute(flt)
#	a.write_image("q0.hdf",2)
	a=a.get_clip(Region(filtrad*2,filtrad*2,shrink.get_xsize(),shrink.get_ysize()))
#	shrink.write_image("q1.hdf",0)
#	a.write_image("q1.hdf",1)
	shrink-=a
#	shrink.write_image("q1.hdf",2)
	a=None
	
	shrink2=shrink.copy()
	shrink2.process_inplace("math.squared")
#	image=EMData()
#	image.read_image(args[0])
#	shrink.write_image("e.mrc")
#	shrink2.write_image("f.mrc")
		
	if len(options.auto)>0 : print "Autobox mode ",options.auto[0]
	
	
	# Reference-based automatic particle picking
	if "ref" in options.auto:
		if not refptcl: error_exit("Reference particles required")
		
		print "Prepare references"
		# refptcls will contain shrunken normalized reference particles
		refptcls=[]
		for n,i in enumerate(refptcl):
			# first a circular mask
			i.process_inplace("normalize.circlemean")
			i.process_inplace("mask.sharp",{"outer_radius":i.get_xsize()/2-1})
						
			ic=i.copy()
			refptcls.append(ic)
			ic.mean_shrink(shrinkfactor)
			# make the unmasked portion mean -> 0
			ic-=float(ic.get_attr("mean_nonzero"))
			ic.process_inplace("normalize.unitlen")
#			ic.write_image("scaled_refs.hdf",-1)

		# prepare a mask to use for local sigma calculaton
		circle=shrink.copy_head()
		circle.to_one()
		circle.process_inplace("mask.sharp",{"outer_radius":options.boxsize/(shrinkfactor*2)-1})
		circle/=(float(circle.get_attr("mean"))*circle.get_xsize()*circle.get_ysize())
		
		ccfmean=shrink.calc_ccf(circle,fp_flag.CIRCULANT)
#		circle.write_image("z0a.hdf")
#		shrink2.write_image("z0b.hdf")
		ccfsig=shrink2.calc_ccf(circle,fp_flag.CIRCULANT)
		ccfmean.process_inplace("math.squared")
		ccfsig-=ccfmean		# ccfsig is now pointwise standard deviation of local mean
		ccfsig.process_inplace("math.sqrt")
#		shrink.write_image("z0.hdf")
#		ccfsig.write_image("z1.hdf")
		
		print "Locating possible particles"
		xs=shrink.get_xsize()
		ys=shrink.get_ysize()
		pks=[]
		for n,i in enumerate(refptcls):
			j=i.get_clip(Region(-(xs-i.get_xsize())/2,-(ys-i.get_ysize())/2,xs,ys))
#			j.write_image("0.%0d.hdf"%n)
			ccfone=shrink.calc_ccf(j,fp_flag.CIRCULANT)
#			ccfone.write_image("a.%0d.hdf"%n)
			ccfone/=ccfsig
#			ccfone.write_image("b.%0d.hdf"%n)
			sig=float(ccfone.get_attr("sigma"))
			ccfone.process_inplace("mask.onlypeaks",{"npeaks":0})
#			ccfone.write_image("c.%0d.hdf"%n)
			pk=ccfone.calc_highest_locations(sig*3.5)
			for m,p in enumerate(pk):
				pk[m]=(-p.value,n,p.x,p.y)
			pks+=pk
			
		pks.sort()		# an ordered list of the best particle locations
		print "%d putative particles located"%len(pks)
				
		# this will produce a new list excluding any lower valued boxes within
		# 1/2 a box size of a higher one. It also rescales the boxes.
		# (ok, you could do this with clever syntax, but this is more readable)

		# first we prepare a grid table of putative box locations for speed
		grid={}
		for n,i in enumerate(pks):
			x=int(floor(i[2]*shrinkfactor/options.boxsize))
			y=int(floor(i[3]*shrinkfactor/options.boxsize))
			try: grid[(x,y)].append(n)
			except: grid[(x,y)]=[n]

		
		goodpks=[]
		bf=options.boxsize/(shrinkfactor*2)
		for n,i in enumerate(pks):
			if i[2]<bf or i[3]<bf or i[2]>xs-bf-1 or i[3]>ys-bf-1 : continue

			# local is a list of putative peaks near the current peak
                        x=int(floor(i[2]*shrinkfactor/options.boxsize))
                        y=int(floor(i[3]*shrinkfactor/options.boxsize))
			local=[]
			for xx in range(x-1,x+2):
				for yy in range(y-1,y+2):
					try: local+=grid[(xx,yy)]
					except: pass
			local=filter(lambda x: x>n,local)

			for nn in local:
				ii=pks[nn]
				if hypot(i[2]-ii[2],i[3]-ii[3])<bf*3/2 : break
			else: goodpks.append([i[0],i[1],i[2]*shrinkfactor-options.boxsize/2,i[3]*shrinkfactor-options.boxsize/2])
		
		print "%d putative particles after local exclusion"%len(goodpks)
		
		print "refine particle locations"
		# This will optimize the center location of each particle and improve
		# the similarity calculation
		goodpks2=[]
		n=0
		for i in goodpks:
			b=EMData()
			
			# on the first pass, rather than just using the best reference determined using the shrunken images
			# we also try a subset of the references to see which gives the best alignment/match
			# refns is the list of ref image #'s to test against
			if options.nretest>0 :
				refns=[i[1]]+range(options.nretest)	# list of the ref numbers of the particles to use in the realignment
			else :
				refns=[i[1]]+options.retestlist
			
			# read in the area where we think a particle exists
			try: b.read_image(initial,0,0,Region(i[2],i[3],options.boxsize,options.boxsize))
			except: continue
			b.process_inplace("normalize.edgemean")
			b.process_inplace("eman1.filter.lowpass.gaussian",{"lowpass":.1})
#			ba=refptcl[i[1]].align("rotate_translate",b,{},"SqEuclidean")
			
			# we iterate over each reference then find the best, eventually recentering ( i[2-3] ) and
			# updating the reference number for the second pass ( i[1] )
			tsts=[]
			for j in refns: 
				ba=b.align("rotate_translate",refptcl[j],{},"optvariance",{"matchfilt":1})
				tsts.append([ba.get_attr("align.score"),j,ba.get_attr("align.dx"),ba.get_attr("align.dy"),ba.get_attr("align.az")])
			tsts.sort()
#			if tsts[0][1]!=i[1] : print i[1]," -> ",tsts[0][1],"    %f,%f  %f"%(tsts[0][2],tsts[0][3],tsts[0][4])
			i[1]=tsts[0][1]
			i[2]-= cos(tsts[0][4])*tsts[0][2]+sin(tsts[0][4])*tsts[0][3]
			i[3]-=-sin(tsts[0][4])*tsts[0][2]+cos(tsts[0][4])*tsts[0][3]

# this code can be used to test alignment
#			b.write_image("cmp.hdf",-1)
#			refptcl[i[1]].write_image("cmp.hdf",-1)
#			ba.write_image("cmp.hdf",-1)
#			try: 
#				b.read_image(args[0],0,0,Region(i[2],i[3],options.boxsize,options.boxsize))
#				b.write_image("cmp.hdf",-1)
#			except: pass
			
			# now we refine this by doing a second pass with the best reference
			try: b.read_image(args[0],0,0,Region(i[2],i[3],options.boxsize,options.boxsize))
			except: continue
			b.process_inplace("normalize.edgemean")
			b.process_inplace("eman1.filter.lowpass.gaussian",{"lowpass":.1})
#			ba=refptcl[i[1]].align("rotate_translate",b,{},"SqEuclidean")
			ba=b.align("rotate_translate",refptcl[i[1]],{},"optvariance",{"matchfilt":1})
			dx=ba.get_attr("align.dx")
			dy=ba.get_attr("align.dy")
			da=ba.get_attr("rotational")
			i[2]-= cos(da)*dx+sin(da)*dy
			i[3]-=-sin(da)*dx+cos(da)*dy
			if hypot(dx,dy)>12.0 or ba.get_attr("ovcmp_m")<=0: 
				print '****'
				continue
			
#			refptcl[i[1]].write_image("at.hdf",-1)
#			ba.write_image("at.hdf",-1)

			# for testing, write out intermediate images 
#			b.write_image("t2.hdf",-1)
#			cc=refptcl[i[1]].rot_scale_trans2D(-da,1.0,-(cos(da)*dx+sin(da)*dy),-(-sin(da)*dx+cos(da)*dy))
#			cc-=ba.get_attr("ovcmp_b")
#			cc/=ba.get_attr("ovcmp_m")
#			cc.write_image("t2.hdf",-1)
#			b-=cc
#			b.write_image("t2.hdf",-1)
			
			# now we calculate the particle quality using the final
			# alignment parameters
			try: b.read_image(args[0],0,0,Region(i[2],i[3],options.boxsize,options.boxsize))
			except: continue
			b.process_inplace("normalize.edgemean")
#			b.process_inplace("eman1.filter.lowpass.gaussian",{"lowpass":.05})
#			print "%d ROT %f"%(n*2,da)
			rr=refptcl[i[1]].rot_scale_trans2D(da,1.0,0,0)
			rr.process_inplace("normalize")
			
#			b.cmp("optvariance",rr,{"keepzero":1})
#			b*=b.get_attr("ovcmp_m")
#			b+=b.get_attr("ovcmp_b")
#			rr.write_image("a.hdf",-1)
#			b.write_image("a.hdf",-1)

#			score=rr.cmp("quadmindot",b,{"normalize":1})+1.0			# This is 1.0-normalized dot product, ie 0 is best 2 is worst
#			score=rr.cmp("phase",b,{})+rr.cmp("optvariance",b,{"radweight":1,"matchamp":1})/rr.get_xsize()
			score=sqrt(rr.cmp("optvariance",b,{"matchfilt":1}))
#			score=sqrt(b.cmp("optvariance",rr,{"matchfilt":1}))
#			score=b.get_attr("ovcmp_m")*b.get_attr("sigma")
#			if (score<=0) : continue

			# now record the fixed up location
#			goodpks2.append((ba.get_attr("align.score")*ba.get_attr("ovcmp_m"),i[2],i[3],i[1],ba.get_attr("ovcmp_m"),n))
			goodpks2.append((score,i[2],i[3],i[1],ba.get_attr("ovcmp_m"),ba.get_attr("ovcmp_b"),n,score))
			print "%d\t%1.2f\t%1.2f\t%1.1f\t%1.4f\t%1.6f"%(n,ba.get_attr("align.dx"),ba.get_attr("align.dy"),ba.get_attr("rotational")*180.0/pi,ba.get_attr("ovcmp_m"),goodpks2[-1][0])
#			ba.write_image("ttt.hdf",-1)
			n+=1
#			display([b,ba,refptcl[i[1]]])
						
		goodpks2.sort()
		
		# now we do 1-D k-means to split the data into 3 groups
		pl=[(goodpks2[0][0],0),((goodpks2[0][0]+goodpks2[-1][0])/2.0,0),(goodpks2[-1][0],0)]
		
		for i in range(20):
			pl2=[(0,0),(0,0),(0,0)]
			for j in goodpks2:
				if j[0]<(pl[0][0]+pl[1][0])/2.0 : pl2[0]=(pl2[0][0]+j[0],pl2[0][1]+1.0)
				elif j[0]>(pl[1][0]+pl[2][0])/2.0 : pl2[2]=(pl2[2][0]+j[0],pl2[2][1]+1.0)
				else : pl2[1]=(pl2[1][0]+j[0],pl2[1][1]+1.0)
#			pl=[(pl2[0][0]/pl2[0][1],pl2[0][1]),(pl2[1][0]/pl2[1][1],pl2[1][1]),(pl2[2][0]/pl2[2][1],pl2[2][1])]
			# complexity here is necessary to deal with empty sets in the k-means
			if pl2[0][1]==0 : pl[0]=(goodpks2[0][0],0)
			else : pl[0]=(pl2[0][0]/pl2[0][1],pl2[0][1])
			if pl2[1][1]==0 : pl[1]=((goodpks2[0][0]+goodpks2[-1][0])/2.0,0)
			else : pl[1]=(pl2[1][0]/pl2[1][1],pl2[1][1])
			if pl2[2][1]==0 : pl[2]=(goodpks2[-1][0],0)
			else : (pl2[2][0]/pl2[2][1],pl2[2][1])
			print pl
		
		
		
		if   (options.threshold<1.0) : thr=pl[0][0]*options.threshold+goodpks2[0][0]*(1.0-options.threshold)
		elif (options.threshold<2.0) : thr=pl[1][0]*(options.threshold-1.0)+pl[0][0]*(2.0-options.threshold)
		elif (options.threshold<3.0) : thr=pl[2][0]*(options.threshold-2.0)+pl[1][0]*(3.0-options.threshold)
		elif (options.threshold<4.0) : thr=goodpks2[-1][0]*(options.threshold-3.0)+pl[2][0]*(4.0-options.threshold)
		else : thr=goodpks2[-1][0]
		
		print "Threshold : ",thr
		boxthr=thr
		
		# put results in standard 'box' array so GUI can modify if desired
		for i in goodpks2:
			boxes.append([i[1],i[2],options.boxsize,options.boxsize,i[0],1])		# x,y,xsize,ysize,quality,changed
			
	# 		out=open("box.stats","w")
# 		for i in goodpks2: out.write("%f\n"%i[0])
# 		out.close()
	
	if "pspec" in options.auto:
		shrink.write_image("shr.mrc")
		for y in range(0,image_size[1]-options.boxsize,options.boxsize/2):
			for x in range(0,image_size[0]-options.boxsize,options.boxsize/2):
				boxes.append([x,y,options.boxsize,options.boxsize,0.0,1])

		for b in boxes:
			cl=shrink.get_clip(Region(b[0]/shrinkfactor,b[1]/shrinkfactor,options.boxsize/shrinkfactor,options.boxsize/shrinkfactor))
			f=cl.do_fft()
			r=f.calc_radial_dist(options.boxsize/shrinkfactor/2-2,2.0,1.0,0)
			plot(r,1)
	
	if "grid" in options.auto:
		try:
			dx=-options.overlap
			if dx+options.boxsize<=0 : dx=0.0
			dy=dx
		except:
			dy=(image_size[1]%options.boxsize)*options.boxsize/image_size[1]-1
			dx=(image_size[0]%options.boxsize)*options.boxsize/image_size[0]-1
			if dy<=0 : dy=((image_size[1]-1)%options.boxsize)*options.boxsize/image_size[1]-1
			if dx<=0 : dx=((image_size[0]-1)%options.boxsize)*options.boxsize/image_size[0]-1
		
#		print image_size,dx,dy,options.boxsize
		for y in range(options.boxsize/2,image_size[1]-options.boxsize,dy+options.boxsize):
			for x in range(options.boxsize/2,image_size[0]-options.boxsize,dx+options.boxsize):
				boxes.append([x,y,options.boxsize,options.boxsize,0.0,1])
				
	if "circle" in options.auto:
		shrinksq=shrink.copy()
		shrinksq*=shrinksq			# shrunken original image squared
		
		# outer and inner ring mask
		outer=EMData()
		sbox=int(options.boxsize/shrinkfactor)
		outer.set_size(shrink.get_xsize(),shrink.get_ysize(),1)
		outer.to_one()
		inner=outer.copy()
		
		outer.process_inplace("mask.sharp",{"inner_radius":sbox*2/5,"outer_radius":sbox/2})
		inner.process_inplace("mask.sharp",{"outer_radius":sbox*2/5})
		
		outer.write_image("b_outer.hdf")
		inner.write_image("b_inner.hdf")

		ccf1=shrinksq.calc_ccf(inner,fp_flag.CIRCULANT)
		ccf2=shrinksq.calc_ccf(outer,fp_flag.CIRCULANT)
		
		ccf1.write_image("b_ccf1.hdf")
		ccf2.write_image("b_ccf2.hdf")
		
		ccf1/=ccf2
		
		ccf1.write_image("b_result.hdf")
	
	E2end(logid)

	# invoke the GUI if requested
	if options.gui:
		gui=GUIbox(args[0],boxes,boxthr,options.boxsize)
		gui.run()

	if options.dbout:
		# Write EMAN1 style box database
		out=open(options.dbout,"w")
		n=0
		for i in boxes:
			if i[4]>boxthr : continue
			out.write("%d\t%d\t%d\t%d\t-3\n"%(i[0],i[1],i[2],i[3]))		
		out.close()

	if options.ptclout:
		# write boxed particles
		n=0
		b=EMData()
		for i in boxes:
			if i[4]>boxthr : continue
			try: b.read_image(args[0],0,0,Region(i[0],i[1],i[2],i[3]))
			except: continue
			if options.norm: b.process_inplace("normalize.edgemean")
#			print n,i
#			print i[4]
			b.write_image(options.ptclout,n)
			n+=1
			if options.savealiref:
				refptcl[i[3]].write_image("boxali.hdf",-1)
				b.write_image("boxali.hdf",-1)
		print "Wrote %d/%d particles to %s"%(n,len(boxes),options.ptclout)
				

try:
	from PyQt4 import QtCore, QtGui, QtOpenGL
	from PyQt4.QtCore import Qt
	from valslider import ValSlider
except:
	print "Warning: PyQt4 must be installed to use the --gui option"
	class dummy:
		pass
	class QWidget:
		"A dummy class for use when Qt not installed"
		def __init__(self,parent):
			print "Qt4 has not been loaded"
	QtGui=dummy()
	QtGui.QWidget=QWidget

class GUIbox:
	def __init__(self,imagefsp,boxes,thr,boxsize=-1):
		"""Implements the 'boxer' GUI. image is the entire image, and boxes and thr specify current boxes
		to begin with. Modified boxes/thr are returned. 
		
		'boxes' is a list of box locations:
		[x0,y0,xsize,ysize,quality,changed]
		quality is used in conjuntion with 'thr' to decide which
		boxes are 'real'. 'changed' is used by the GUI to decide when
		redrawing is necessary (should be set to 1 initially)."""
		try:
			from emimage import EMImage,get_app
		except:
			print "Cannot import EMAN image GUI objects (emimage,etc.)"
			sys.exit(1)
		
		if len(boxes)>0 and boxsize==-1: self.boxsize=boxes[0][2]
		elif boxsize==-1: self.boxsize=128
		else: self.boxsize=boxsize
		
		self.app=get_app()
		self.image=EMData()					# the image to be boxed
		self.image.read_image(imagefsp)
		self.boxes=boxes					# the list of box locations
		self.threshold=thr					# Threshold to decide which boxes to use
		self.ptcl=[]						# list of actual boxed out EMImages
		
		self.moving=None					# Used during a user box drag
		
		self.guiim=EMImage(self.image)		# widget for displaying large image
		self.guimx=EMImage(self.ptcl)		# widget for displaying matrix of smaller images
		
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mousedown"),self.mousedown)
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mousedrag"),self.mousedrag)
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mouseup")  ,self.mouseup  )
		self.guimx.connect(self.guimx,QtCore.SIGNAL("mousedown"),self.boxsel)
		
		self.guimx.mmode="app"
		self.guictl=GUIboxPanel(self)
		
		try:
			E2loadappwin("boxer","imagegeom",self.guiim)
			E2loadappwin("boxer","matrixgeom",self.guimx)
			E2loadappwin("boxer","controlgeom",self.guictl)
			self.guimx.nperrow=E2getappval("boxer","matrixnperrow")
			
			if E2getappval("boxer","imcontrol") :
				self.guiim.showInspector(True)
				E2loadappwin("boxer","imcontrolgeom",self.guiim.inspector)
			if E2getappval("boxer","mxcontrol") :
				self.guimx.showInspector(True)
				E2loadappwin("boxer","mxcontrolgeom",self.guimx.inspector)
		except:
			pass
		
		self.guiim.show()
		self.guimx.show()
		self.guictl.show()
		
		self.boxupdate()

	def boxsel(self,event,lc):
		im=lc[0]
		self.guiim.setActive(im,.9,.9,.4)
		self.guimx.setSelected(im)
		try: self.guiim.scrollTo(self.boxes[im][0]+self.boxes[im][2]/2,self.boxes[im][1]+self.boxes[im][3]/2)
		except: print "boxsel() scrolling error"

	def mousedown(self,event) :
		m=self.guiim.scr2img((event.x(),event.y()))
		print event.x(),event.y(),m[0],m[1]
		for i,j in enumerate(self.boxes):
			if m[0]<j[0] or m[0]>j[0]+j[2] or m[1]<j[1] or m[1]>j[1]+j[3] : continue
			break
		else: 
			if event.modifiers()&Qt.ShiftModifier : return
			# If we get here, we need to make a new box
			self.boxes.append([m[0]-self.boxsize/2,m[1]-self.boxsize/2,self.boxsize,self.boxsize,0,1])
			i=len(self.boxes)-1
			self.boxupdate()
		# an existing box
		if event.modifiers()&Qt.ShiftModifier :
			# with shift, we delete
			self.delbox(i)
			return
		self.moving=(i,m)
		self.guiim.setActive(i,.9,.9,.4)
		x0=self.boxes[i][0]+self.boxes[i][2]/2-1
		y0=self.boxes[i][1]+self.boxes[i][3]/2-1
		self.guiim.addShape("cen",EMShape(["rect",.9,.9,.4,x0,y0,x0+2,y0+2,1.0]))
		if not self.guimx.isVisible(i) : self.guimx.scrollTo(i,yonly=1)
		self.guimx.setSelected(i)
		
	def mousedrag(self,event) :
		m=self.guiim.scr2img((event.x(),event.y()))
		
		# box deletion when shift held down
		if event.modifiers()&Qt.ShiftModifier:
			for i,j in enumerate(self.boxes):
				if m[0]<j[0] or m[0]>j[0]+j[2] or m[1]<j[1] or m[1]>j[1]+j[3] : continue
				break
			else: return
			self.delbox(i)
			return
			
		if self.moving:
			self.boxes[self.moving[0]][0]+=m[0]-self.moving[1][0]
			self.boxes[self.moving[0]][1]+=m[1]-self.moving[1][1]
			self.boxes[self.moving[0]][-1]=1
			self.moving=(self.moving[0],m)
			
			# update the center marker
			i=self.moving[0]
			x0=self.boxes[i][0]+self.boxes[i][2]/2-1
			y0=self.boxes[i][1]+self.boxes[i][3]/2-1
			self.guiim.addShape("cen",EMShape(["rect",.9,.9,.4,x0,y0,x0+2,y0+2,1.0]))
			
			self.boxupdate()
		
	def mouseup(self,event) :
		m=self.guiim.scr2img((event.x(),event.y()))
		self.moving=None
	
	def delbox(self,i):
		"Deletes the numbered box completely"
		sh=self.guiim.shapes
		k=sh.keys()
		k.sort()
		del sh[i]
		for j in k:
			if isinstance(j,int):
				if j>i :
					sh[j-1]=sh[j]
					del sh[j]
		self.guiim.delShapes()
		self.guiim.addShapes(sh)
		self.guiim.setActive(None,.9,.9,.4)
		
		del(self.boxes[i])
		del(self.ptcl[i])
		self.guimx.setData(self.ptcl)
	
	def boxupdate(self,force=False):
		goodboxes=[i for i in self.boxes if i[4]<=self.threshold]
		if len(self.ptcl)>len(goodboxes) :
			del self.ptcl[len(goodboxes):]
			self.guiim.delShapes()

		ns={}
		for j,i in enumerate(goodboxes):
			if not i[5] and not force: continue
			
			self.boxes[j][-1]=0
			im=self.image.get_clip(Region(i[0],i[1],i[2],i[3]))
			ns[j]=EMShape(["rect",.4,.9,.4,i[0],i[1],i[0]+i[2],i[1]+i[3],2.0])
			if j>=len(self.ptcl) : self.ptcl.append(im)
			else : self.ptcl[j]=im

		self.guiim.addShapes(ns)
		self.guimx.setData(self.ptcl)
		
		self.guictl.nboxesChanged(len(self.ptcl))
#		self.emit(QtCore.SIGNAL("nboxes"),len(self.ptcl))

	def run(self):
		"""If you make your own application outside of this object, you are free to use
		your own local app.exec_(). This is a convenience for boxer-only programs."""
		self.app.exec_()
		
		E2saveappwin("boxer","imagegeom",self.guiim)
		E2saveappwin("boxer","matrixgeom",self.guimx)
		E2saveappwin("boxer","controlgeom",self.guictl)
		E2setappval("boxer","matrixnperrow",self.guimx.nperrow)
		try:
			E2setappval("boxer","imcontrol",self.guiim.inspector.isVisible())
			if self.guiim.inspector.isVisible() : E2saveappwin("boxer","imcontrolgeom",self.guiim.inspector)
		except : E2setappval("boxer","imcontrol",False)
		try:
			E2setappval("boxer","mxcontrol",self.guimx.inspector.isVisible())
			if self.guimx.inspector.isVisible() : E2saveappwin("boxer","mxcontrolgeom",self.guimx.inspector)
		except : E2setappval("boxer","mxcontrol",False)
		
		return (self.boxes,self.threshold)

class GUIboxPanel(QtGui.QWidget):
	def __init__(self,target) :
		
		QtGui.QWidget.__init__(self,None)
		self.target=target
		
		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		
		self.info = QtGui.QLabel("%d Boxes"%len(target.boxes),self)
		self.vbl.addWidget(self.info)

		self.thr = ValSlider(self,(0.0,3.0),"Threshold:")
		self.thr.setValue(target.threshold)
		self.vbl.addWidget(self.thr)

		self.hbl1=QtGui.QHBoxLayout()
		self.hbl1.setMargin(0)
		self.hbl1.setSpacing(2)
		self.vbl.addLayout(self.hbl1)
		
		self.lblbs=QtGui.QLabel("Box Size:",self)
		self.hbl1.addWidget(self.lblbs)
		
		self.bs = QtGui.QLineEdit(str(target.boxsize),self)
		self.hbl1.addWidget(self.bs)

		self.hbl2=QtGui.QHBoxLayout()
		self.hbl2.setMargin(0)
		self.hbl2.setSpacing(2)
		self.vbl.addLayout(self.hbl2)

		self.done=QtGui.QPushButton("Done")
		self.vbl.addWidget(self.done)
		
		self.connect(self.bs,QtCore.SIGNAL("editingFinished()"),self.newBoxSize)
		self.connect(self.thr,QtCore.SIGNAL("valueChanged"),self.newThresh)
		self.connect(self.done,QtCore.SIGNAL("clicked(bool)"),self.target.app.quit)
#		self.target.connect(self.target,QtCore.SIGNAL("nboxes"),self.nboxesChanged)
		
	def nboxesChanged(self,n):
		self.info.setText("%d Boxes"%n)
	
	def newBoxSize(self):
		try:
			v=int(self.bs.text())
			if v<12 : raise Exception
		except:
			self.bs.setText(str(self.target.boxsize))
			return
		
		for i in self.target.boxes:
			if i[2]!=v or i[3]!=v : i[5]=1
			else: continue
			i[0]-=(v-i[2])/2
			i[1]-=(v-i[3])/2
			i[2]=v
			i[3]=v
		
		self.target.boxsize=v
		self.target.boxupdate()
		
	def newThresh(self,val):
		self.target.threshold=val
		self.target.boxupdate(True)


if __name__ == "__main__":
	main()
