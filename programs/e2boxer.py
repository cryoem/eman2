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
from time import *
import os
import sys

from copy import *


from emglplot import *

pl=()

THRESHOLD = "Threshold"
SELECTIVE = "Selective"
MORESELECTIVE = "More Selective"

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <image>
	
Automatic and manual particle selection. This version is specifically aimed at square boxes
for single particle analysis."""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--gui",action="store_true",help="Start the GUI for interactive boxing",default=False)
	parser.add_option("--boxsize","-B",type="int",help="Box size in pixels",default=-1)
	parser.add_option("--dbin","-D",type="string",help="Filename to read an existing box database from",default=None)
	parser.add_option("--auto","-A",type="string",action="append",help="Autobox using specified method: ref, grid",default=[])
	parser.add_option("--overlap",type="int",help="(auto:grid) number of pixels of overlap between boxes. May be negative.")
	parser.add_option("--refptcl","-R",type="string",help="(auto:ref) A stack of reference images. Must have the same scale as the image being boxed.",default=None)
	parser.add_option("--nretest",type="int",help="(auto:ref) Number of reference images (starting with the first) to use in the final test for particle quality.",default=-1)
	parser.add_option("--retestlist",type="string",help="(auto:ref) Comma separated list of image numbers for retest cycle",default="")
	parser.add_option("--farfocus",type="string",help="filename or 'next', name of an aligned far from focus image for preliminary boxing",default=None)
	parser.add_option("--dbout",type="string",help="filename to write EMAN1 style box database file to",default=None)
	parser.add_option("--norm",action="store_true",help="Edgenormalize boxed particles",default=False)
	parser.add_option("--ptclout",type="string",help="filename to write boxed out particles to",default=None)
	
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
	# what is this attribute for?
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
		else : parser.error("Please specify a box size")
	else:
		if not options.boxsize in good_box_sizes:
			print "Note: EMAN2 processing would be more efficient with a boxsize of %d"%good_boxsize(options.boxsize)
	
	try: options.retestlist=[int(i) for i in options.retestlist.split(',')]
	except: options.retestlist=[]
			
	shrinkfactor=int(ceil(options.boxsize/16))
	print "Shrink factor = ",shrinkfactor
	#shrinkfactor=int(ceil(image.get_ysize()/1024.0))
	#if options.boxsize/shrinkfactor<12 : shrinkfactor/=2
	
	image.process_inplace("normalize.edgemean")
	shrink=image
	shrink.process_inplace("math.meanshrink",{"n":shrinkfactor})		# shrunken original image
	
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
	flt.process_inplace("xform.phaseorigin.tocorner")
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
			ic.process_inplace("math.meanshrink",{"n":shrinkfactor})
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
	
	E2end(logid)

	# invoke the GUI if requested
	if options.gui:
		gui=GUIbox(args[0],boxes,boxthr,options.boxsize)
		gui.run()
		
	print "finished running"

	if options.dbout:
		boxes = gui.getboxes()
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

class Box:
	def __init__(self,xcorner=-1,ycorner=-1,xsize=-1,ysize=-1,isref=0,correlationscore=0):
		self.xcorner = xcorner			# the xcorner - bottom left
		self.ycorner = ycorner			# the ycorner - bottom left
		self.xsize = xsize				# the xsize of the box
		self.ysize = ysize				# the ysize of the box
		self.isref = isref				# a flag that can be used to tell if the box is being used as a reference
		self.correlationscore = correlationscore	# the correlation score
		
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
		
		# Experimental things that will probably be needed
		
		self.otheridx1 = -1		# stores a unique idx, for BoxSet convenience and efficiency
		self.otheridx2 = -1		# stores a unique idx, for BoxSet convenience and efficiency
		
	
	def updateBoxImage(self,image,norm=True):
		#print "getting region",self.xcorner,self.ycorner,self.xsize,self.ysize
		self.image = image.get_clip(Region(self.xcorner,self.ycorner,self.xsize,self.ysize))
		if norm:
			self.image.process_inplace("normalize.edgemean")
		
		# make sure there are no out of date footprints hanging around
		self.footprint = None
		
			
	def getBoxImage(self,image,norm=True,force=False):
		if self.image == None or force:
			self.updateBoxImage(image,norm)
		return self.image
			

	def getFootPrint(self,shrink=1):
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
			
class BoxSet2:
	def __init__(self,image,parent=None):
		self.image = image			# the image containing the boxes
		self.parent = parent		# keep track of the parent in case we ever need it
		self.boxes = []				# a list of boxes
		self.refboxes = []			# a list of boxes
		self.shrink = -1			# the amount by which the subject image is shrunken before the correlation image is generated. -1 means recalculate the shrink factor. This is an important parameter that speeds autoboxing significantly
		self.boxsize = -1			#  the boxsize
		self.smallimage = None		# a small copy of the image which has had its background flattened
		self.flattenimager = -1		# the r value used to run the flatten image processor
		self.searchradius = -1		# search radius in the correlation image. Important parameter
		
		self.optprofile = None		# An optimum correlation profile, used as the basis of selective autoboxing
		self.optthreshold = None		# The correlation threshold, used for autoboxing
		self.optprofileradius = None	# A point in the optprofile that is used for selective picking
		self.autoboxmethod = SELECTIVE # The method of autoboxing
		self.templateupdate = True	# When turned off, prevents the correlation map from being updated - this is useful the act of adding a ref only affects the optimum parameters, not the correlation map
		self.fpshrink = -1
		self.exclusionimage = None
		
		
	def addbox(self,box):
		if not isinstance(box,Box):
			print "You can not add a box to this box set if it is not of type Box"
			return;
		
		box.isref = True # make sure it knows that it's a reference box
		box.otheridx1 = len(self.boxes) # store this, it's used when the user deletes the box, for efficiency
		box.otheridx2 = len(self.refboxes) # store this, it's used when the user deletes the box, for efficiency
		
		box.rorig = 0			# RGB red
		box.gorig = 0			# RGB green
		box.borig = 0			# RGB blue
		box.r = 0
		box.g = 0
		box.b = 0
		
		#print "adding box",box.xcorner,box.ycorner,box.xsize,box.ysize
		self.boxes.append(box)
		self.refboxes.append(box)
	
	def delbox(self,i):
		tmp = self.boxes.pop(i)
		#yuck, this is horribly inefficient
		for j,box in enumerate(self.refboxes):
			if box.isref and box.xcorner == tmp.xcorner and box.ycorner == tmp.ycorner:
				self.refboxes.pop(j)
				return True
			
		return False
	
	def addnonrefs(self,boxes):
		'''
		Add boxes that are stored in eman1 format
		box[0] = xnorner, box[1] = ycorner, box[2] = xsize, box[3] = ysize
		'''
		for box in boxes:
			b = Box(box[0],box[1],box[2],box[3])
			b.isref = False
			b.changed = True
			self.boxes.append(b)

	def numboxes(self):
		return len(self.boxes)
	
	def updateboxsize(self,boxsize):
		'''
		Updates only the box size and corner coordinates
		Switches the changed flag to True to trigger redisplay (but the calling function
		is responsible for knowing and testing for this)
		'''
		for box in self.boxes:
			if box.xsize != boxsize:
				box.xcorner -= (boxsize-box.xsize)/2
				box.xsize = boxsize
				box.changed = True
			if box.ysize != boxsize:
				box.ycorner -= (boxsize-box.ysize)/2
				box.ysize = boxsize
				box.changed = True

		self.boxsize = boxsize
		
	def getfootprintshrink(self):
		if self.fpshrink == -1:
			shrink = 1
			tn = self.boxsize/2
			while ( tn >= 32 ):
				tn /= 2
				shrink *= 2
			self.fpshrink = shrink
		
		return self.fpshrink
		
	def getbestshrink(self):
		if self.image == None or self.boxsize == -1:
			print "error - either the image is not set, or the boxsize is not set"
			exit(1)
		
		if self.shrink == -1:
			shrink = 1
			inx = self.image.get_xsize()/2
			iny = self.image.get_ysize()/2
			tn = self.boxsize/2
			while ( inx >= 512 and iny >= 512 and tn >= 16 ):
				inx /= 2
				iny /= 2
				tn /= 2
				shrink *= 2
		
			self.shrink=shrink
		
		return self.shrink
	
	def updatetemplate(self,boxsize=-1):
		'''
		update template implicitly updates the correlation image too.
		It's generally called when you know that the reference images have changed
		'''
		#print 'in update template'
		# You can turn off template updating - but it only makes sense if there is
		# not currently a store correlation map. You would do this is you were
		# adding refs to changes the autoboxing parameters only, as opposed to doing that plus
		# updating the correlation map
		if self.templateupdate == False: return
		
		# Warning - read the error statement I am printing
		if boxsize != -1:
			self.boxsize = boxsize
		elif self.boxsize == -1:
			print "error, the first time you call update you must specify the box size"
			exit(1)
		
		
		# Update the template internally, this makes self.template
		self.genrotalignedaverage()
		if self.template == None:
			print "Error, something went wrong with the template generation"
			exit(1)
		
		# newr is the parameter that will potentially be used to run the flattenbackround processor
		newr = self.template.get_xsize()/2.0
		# now we only recalculate the small copy of the subject image if necessary
		if self.smallimage == None or newr != self.flattenimager:
			self.flattenimager = newr
			#FIXME - we could avoid a deep copy by writing the meanshrink processor
			# i.e. section = self.image.process("math.meanshrink",{"n":self.getbestshrink()}
			
			if (self.getbestshrink() != 1):
				self.smallimage = self.image.process("math.meanshrink",{"n":self.getbestshrink()})
			else: self.smallimage = self.image.copy()
			
			self.smallimage.process_inplace("filter.flattenbackground",{"radius":self.flattenimager})

		self.correlation = self.smallimage.calc_flcf( self.template )
		
		# this may not be necessary if we ever want to be completely efficient
		self.correlation.process_inplace("xform.phaseorigin.tocenter")
		self.correlation.write_image("tttttt.hdf")
	
	def genrotalignedaverage(self):
		# you can only generate a template if there are references
		if len(self.refboxes) <= 0: 
			print 'error, cant call genrotalignedaverage when there are no refboxes'
			exit(1)
		
		images_copy = []
		for i in self.refboxes:
			image = i.getBoxImage(self.image)
			if (self.getbestshrink() != 1):
				e = image.process("math.meanshrink",{"n":self.getbestshrink()})
			else : e = image.copy()
			images_copy.append(e)
			
		ave = images_copy[0].copy()
		
		for i in range(1,len(images_copy)):
			#ta = images_copy[i].align("rotate_translate",ave,{},"dot",{"normalize":1})
			ave.add(images_copy[i])
			
		#ave.write_image("prealigned.hdf")
		ave.mult(1.0/len(images_copy))
		ave.process_inplace("math.radialaverage")
		ave.process_inplace("normalize.edgemean")
		ave.process_inplace("mask.sharp",{'outer_radius':ave.get_xsize()/2})
		#ave.write_image("ave.hdf")
		
		# 5 is a magic number
		for n in range(0,5):
			t = []
			for i in images_copy:
				#FIXME - make it so that a newly clipped portion of the original image
				# is used as the 'aligned' image, to avoid zeroing effects at the edges
				ta = i.align("translational",ave,{},"dot",{"normalize":1})
				t.append(ta)
		
			ave = t[0].copy()
			for i in range(1,len(images_copy)):
				ave.add(t[i])
				
			ave.mult(1.0/len(t))
			ave.process_inplace("math.radialaverage")
			ave.process_inplace("normalize.edgemean")
			ave.process_inplace("mask.sharp",{'outer_radius':ave.get_xsize()/2})
		
		self.template = ave
		self.template.write_image("template.hdf")
	
	def addExclusionArea(self, type,x,y,radius):
		
		xx = int(x/self.getbestshrink())
		yy = int(y/self.getbestshrink())
		
		rr = int(radius/self.getbestshrink())
		rrs = rr**2
		#print xx,yy,rr
		
		# this does implicit initialization
		self.getExclusionImage()
		
		ny = self.correlation.get_ysize()
		nx = self.correlation.get_xsize()
		for j in range(-rr,rr):
			for i in range(-rr,rr):
				if (i**2 + j**2) > rrs: continue
				jj = j+yy
				ii = i+xx
				if jj >= ny or jj < 0:continue
				if ii >= nx or ii < 0:continue
				
				self.exclusionimage.set(ii,jj,0)
				self.excl2.set(ii,jj,0.1)
	
	def getExclusionImage(self):
		if self.exclusionimage == None:
			self.exclusionimage = EMData(self.correlation.get_xsize(),self.correlation.get_ysize())
			self.exclusionimage.to_one()
			self.excl2 = EMData(self.correlation.get_xsize(),self.correlation.get_ysize())
			self.excl2.to_zero()
		return self.exclusionimage
	
	def updateefficiency(self,efficiency):
		'''
		paints black circles in the efficiency - which is a binary EMData object
		useful for making things efficient, should probably be called updateexclusions
		'''
		
		oldsearchr = self.searchradius
		self.searchradius = int(0.5*(self.boxsize)/float(self.getbestshrink()))
		
		if self.searchradius != oldsearchr and oldsearchr != -1:
			print "warning, the search radius changed or. Take note david"
		
		for box in self.boxes:
			xx = box.xcorner + box.xsize/2
			yy = box.ycorner + box.ysize/2
			xx /= self.getbestshrink()
			yy /= self.getbestshrink()
			
			BoxingTools.set_radial_zero(efficiency,int(xx),int(yy),self.searchradius)
			
	def accrueparams(self,boxes,center=True):
		if (self.correlation == None):
			print "Error, can't accrue params if now correlation map exists"
			exit(1)

		invshrink = 1.0/self.getbestshrink()
		
		#print "accruing params with a total",len(boxes),"boxes"
		for box in boxes:
			
			# the central coordinates of the box in terms of the shrunken correlation image
			x = (box.xcorner+box.xsize/2.0)*invshrink
			y = (box.ycorner+box.ysize/2.0)*invshrink
			
			#the search radius is used in correlation space - it limits the radial distance
			# up to which 'profile' data can be accrued
			# it is currently half the boxsize in terms of the correlation image's dimensions
			self.searchradius = int(0.5*(self.boxsize)/float(self.getbestshrink()))
			
			peak_location = BoxingTools.find_radial_max(self.correlation,int(x),int(y), self.searchradius )
			peak_location2 = BoxingTools.find_radial_max(self.correlation,peak_location[0],peak_location[1],self.searchradius )
			if (peak_location != peak_location2):
				# this represents a troubling condition
				box.correlationscore = None
				print "Error, peak location unrefined"
				continue
			
			# store the peak location
			box.corx = peak_location[0]
			box.cory = peak_location[1]
			
			# store the correlation value at the correlation max
			box.correlationscore = self.correlation.get(box.corx,box.cory)
			
			# store the profile
			box.optprofile = BoxingTools.get_min_delta_profile(self.correlation,box.corx,box.cory, self.searchradius )
			
			# center on the correlation peak
			if (center):
				box.xcorner = box.corx*self.getbestshrink()-box.xsize/2.0
				box.ycorner = box.cory*self.getbestshrink()-box.ysize/2.0
				box.changed = True
			
			#l = self.searchradius 
			#im = self.correlation.get_clip(Region(box[7]-l/2.0,box[8]-l/2.0,l,l))
			#im.write_image(self.outfile,-1)
			
	def autobox(self,optprofile=None,thr=None):
		if (self.correlation == None):
			print "Error, can't autobox of the correlation image doesn't exist"
			exit(1)
		
		if len(self.refboxes) == 0 :
			print "Error, can't autobox if there are no selected reference boxes"
			exit(1)
			
		efficiency = self.getExclusionImage().copy()
		self.updateefficiency(efficiency)
		efficiency.write_image("efficiency.mrc")
		self.parent.guiim.setOtherData(self.excl2,self.getbestshrink(),True)
		
		# we must accrue the the parameters of the reference images if the optprofile
		# or thr value has not been set in the function arguments
		if ( optprofile == None or thr == None ):
			self.accrueparams(self.refboxes)

		
		# Determine what will be the optimum correlation threshold
		if thr == None:
			# If we must determine the threshold from what we've got, iterate through all of the reference
			# boxes and use the lowest correlation score as the correlation threshold
			found = False
			# POTENTIAL PROBLEM - i have hard coded the use of self.refboxes here for the optimum parameter generation
			# we may want to be able to use other box sets for generating parameters
			for i,box in enumerate(self.refboxes):
				if box.correlationscore == None:
					print "continuing on faulty"
					continue
				if found == False:
					self.optthreshold = box.correlationscore
					found = True
				else:	
					if box.correlationscore < self.optthreshold: self.optthreshold = box.correlationscore
		else:
			self.optthreshold = thr
				
		# Determine what to use as the optimum profile
		if optprofile == None:
			# Iterate through the reference boxes and accrue what you can think of
			# as the worst case scenario, in terms of correlation profiles
			found = False
			# POTENTIAL PROBLEM - i have hard coded the use of self.refboxes here for the optimum parameter generation
			# we may want to be able to use other box sets for generating parameters
			for i,box in enumerate(self.refboxes):
				if box.correlationscore == None:
					##print "continuing on faulty" - this was already printed above
					continue
				if found == False:
					self.optprofile = box.optprofile
					n = len(self.optprofile)
					found = True
				else:
					profile = box.optprofile
					for j in range(0,n):
						if profile[j] < self.optprofile[j]: self.optprofile[j] = profile[j]
			# Tell the parent to update the data it is displaying to the user in a GUI somewhere
			self.parent.updatedata(self.optprofile,self.optthreshold)
		else:
			self.optprofile=optprofile
	
	
		# determine the point in the profile where the drop in correlation score is the greatest, store it in radius
		self.optprofileradius = 0
		tmp = self.optprofile[0]
		for i in range(1,len(self.optprofile)):
			if self.optprofile[i] > tmp and tmp > 0:
				tmp = self.optprofile[i]
				self.optprofileradius = i
			
		
		#print self.optprofile
		#print "using opt radius",self.radius, "which has value",tmp,"shrink was",self.shrink
		if self.autoboxmethod == THRESHOLD:
			mode = 0
		elif self.autoboxmethod == SELECTIVE:
			mode = 1
		elif self.autoboxmethod == MORESELECTIVE:
			mode = 2
		
		soln = BoxingTools.auto_correlation_pick(self.correlation,self.optthreshold,self.searchradius,self.optprofile,efficiency,self.optprofileradius,mode)

		for b in soln:
			x = b[0]
			y = b[1]
			xx = int(x*self.shrink)
			yy = int(y*self.shrink)
			box = Box(xx-self.boxsize/2,yy-self.boxsize/2,self.boxsize,self.boxsize,0)
			box.correlationscore =  self.correlation.get(x,y)
			box.corx = b[0]
			box.cory = b[1]
			box.changed = True
			self.boxes.append(box)
		self.parent.boxupdate()
		
	def classify(self):
		v = []
		# accrue all params
		self.accrueparams(self.boxes)
		
		for box in self.boxes:
			b = copy(box.optprofile[0:self.radius])
			b.sort()
			#for a in b:
				#a = box[6]-a
			#print b
			v.append(b)
			
		cl = BoxingTools.classify(v,4)
		self.parent.updateboxcolors(cl)
	
	def gen_ref_images(self):
		tmpimage = "tmpparticles.img"
		self.parent.write_boxes_to(tmpimage)
		
		self.process = QtCore.QProcess()

		program = QtCore.QString("e2refine2d.py")
		args = QtCore.QStringList()
		args.append("--input="+tmpimage)
		args.append("--ncls=15")
		
		QtCore.QObject.connect(self.process, QtCore.SIGNAL("finished(int)"), self.process_finished)
		QtCore.QObject.connect(self.process, QtCore.SIGNAL("started()"), self.process_start)
		print self.process.start(program,args)

	def process_start(self):
		print "received process start signal"
		
	def boxsel(self,event,lc):
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
		self.parent.boxupdate()
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
		QtCore.QObject.connect(self.imagemx2,QtCore.SIGNAL("mousedown"),self.boxsel)
		self.imagemx2p.show()
		
		ef = []
		for image in e:
			image.process_inplace("normalize.edgemean")
			if self.getbestshrink() != 1:
				image = image.process("math.meanshrink",{"n":self.getfootprintshrink()})	
			ef.append(image.make_footprint())
		
		for box in self.boxes:
			best = -1
			group = -1
			for i,g in enumerate(ef): 
				s = box.getFootPrint(self.getfootprintshrink()).cmp("optvariance",g,{"matchfilt":1,"matchamp":1})
				# REMEMBER - cmp returns values that have potentially been negated - a smaller value is better
				if best == -1 or s < best:
					group = i
					best = s
			
			box.group = group
					
		
		#print scores
		
		print "received finish signal"

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
		
		self.eraseradius = 50
		
		self.app=get_app()
		self.image=EMData()					# the image to be boxed
		self.image.read_image(imagefsp)
		self.imagefsp = imagefsp
		#if abs(self.image.get_attr("mean")) < 1:
			#print "adding 10"
			#self.image.add(10)
		self.boxset2 = BoxSet2(self.image,self)
		self.boxset2.addnonrefs(boxes)
		self.threshold=thr					# Threshold to decide which boxes to use
		self.ptcl=[]						# list of actual boxed out EMImages
		self.boxm = None
		self.moving=None					# Used during a user box drag
		
		self.guiimp=EMImage(self.image)		# widget for displaying large image
		self.guiim=self.guiimp.child
		self.guimxp=EMImage(self.ptcl)		# widget for displaying matrix of smaller images
		self.guimx=self.guimxp.child
		
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mousedown"),self.mousedown)
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mousedrag"),self.mousedrag)
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mouseup")  ,self.mouseup  )
		self.guiim.connect(self.guiim,QtCore.SIGNAL("keypress"),self.keypress)
		self.guimx.connect(self.guimx,QtCore.SIGNAL("removeshape"),self.removeshape)
		self.guimx.connect(self.guimx,QtCore.SIGNAL("mousedown"),self.boxsel)
		self.guimx.connect(self.guimx,QtCore.SIGNAL("mousedrag"),self.boxmove)
		self.guimx.connect(self.guimx,QtCore.SIGNAL("mouseup"),self.boxrelease)
		
		self.mmode = 0
		self.guiim.setmmode(0)
		self.guimx.setmmode("app")
		self.ppc = 1.0
		self.mouseclicks = 0
		self.guictl=GUIboxPanel(self)
		
		
		self.ap = False
		
		try:
			E2loadappwin("boxer","imagegeom",self.guiimp)
			E2loadappwin("boxer","matrixgeom",self.guimxp)
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
		
		self.guiimp.show()
		self.guimxp.show()
		self.guictl.show()
		
		self.boxupdate()
	def write_boxes_to(self,imagename,norm=True):
		boxes = self.getboxes()
		
		# Write EMAN1 style box database
		n = 0
		for box in boxes:
			image = box.getBoxImage(self.image)
#			print n,i
#			print i[4]
			image.write_image(imagename,n)
			n += 1

	def boxmove(self,event,scale):
		dx = (self.boxm[0] - event.x())/scale
		dy = (event.y() - self.boxm[1])/scale
		self.movebox(self.boxm[2],dx,dy)
		self.boxm[0] = event.x()
		self.boxm[1] = event.y()
		self.updateAllImageDisplay()
	def boxrelease(self,event):
		if self.getboxes()[self.boxm[2]].isref :
			self.autobox()
		self.boxm = None

	def boxsel(self,event,lc):
		im=lc[0]
		self.boxm = [event.x(),event.y(),im]
		self.guiim.setActive(im,.9,.9,.4)
		self.guimx.setSelected(im)
		boxes = self.getboxes()
		self.guiim.registerScrollMotion(boxes[im].xcorner+boxes[im].xsize/2,boxes[im].ycorner+boxes[im].ysize/2)
		try:
			#self.guiim.scrollTo(boxes[im].xcorner+boxes[im].xsize/2,boxes[im].ycorner+boxes[im].ysize/2)
			pass
			
		except: print "boxsel() scrolling error"

	def getboxes(self):
		return self.boxset2.boxes
	
	def mousedown(self,event) :
		if self.mmode == 0:
			m=self.guiim.scr2img((event.x(),event.y()))
			collision = False
			# basic strategy is to detect if there was a collision with any other box
			# do this by incrementing through all the box sets
			boxes = self.getboxes()
			
			boxnum = self.collisiondetect(m,boxes)
			if boxnum == -1:
				#if we make it here, that means the user has clicked on an area that is not in any box
				if event.modifiers()&Qt.ShiftModifier : return # the user tried to delete nothing
				# If we get here, we need to make a new box
				box = Box(m[0]-self.boxsize/2,m[1]-self.boxsize/2,self.boxsize,self.boxsize,True)
				box.changed = True
				self.boxset2.addbox(box)
				self.mouseclicks += 1
				self.boxupdate()
				if self.ap:	self.autobox()
				self.updateppc()
				self.updateAllImageDisplay()
				boxnum = len(self.getboxes())-1
				x0=boxes[boxnum].xcorner+boxes[boxnum].xsize/2-1
				y0=boxes[boxnum].ycorner+boxes[boxnum].ysize/2-1
				self.guiim.addShape("cen",EMShape(["rect",.9,.9,.4,x0,y0,x0+2,y0+2,1.0]))
				if not self.guimx.isVisible(boxnum) : self.guimx.scrollTo(boxnum,yonly=1)
				self.guimx.setSelected(boxnum)
				return
				
				
			self.mouseclicks += 1
			# Deleting a box happens here
			if event.modifiers()&Qt.ShiftModifier :
				# with shift, we delete
				self.mouseclicks += 1
				self.delbox(boxnum)
				self.updateAllImageDisplay()
				self.updateppc()
				return 
			
			
			# if we make it here than the we're moving a box
			self.moving=[boxes[boxnum],m,boxnum]
			self.guiim.setActive(boxnum,.9,.9,.4)
				
			x0=boxes[boxnum].xcorner+boxes[boxnum].xsize/2-1
			y0=boxes[boxnum].ycorner+boxes[boxnum].ysize/2-1
			self.guiim.addShape("cen",EMShape(["rect",.9,.9,.4,x0,y0,x0+2,y0+2,1.0]))
			if not self.guimx.isVisible(boxnum) : self.guimx.scrollTo(boxnum,yonly=1)
			self.guimx.setSelected(boxnum)
			self.updateAllImageDisplay()
		elif self.mmode == 1:
			m=self.guiim.scr2img((event.x(),event.y()))
			s = EMShape(["circle",.9,.9,.9,m[0],m[1],self.eraseradius,5])
			self.boxset2.addExclusionArea("circle",m[0],m[1],self.eraseradius)
			self.guiim.addShape("eraser",EMShape(["circle",.9,.9,.9,m[0],m[1],self.eraseradius,5]))
			
	def updateppc(self):
		self.guictl.ppcChanged(len(self.getboxes())/float(self.mouseclicks))
	
	def collisiondetect(self,m,boxes):
			
		for boxnum,box in enumerate(boxes):
			if m[0]<box.xcorner or m[0]>(box.xcorner +box.xsize) or m[1]<box.ycorner or m[1]>(box.ycorner +box.ysize) :
				# no collision
				continue
			# if we make it here there has been a collision, the box already exists
			return boxnum
		
		return -1
	
	def mousedrag(self,event) :
		if self.mmode == 0:
			m=self.guiim.scr2img((event.x(),event.y()))
			
			if event.modifiers()&Qt.ShiftModifier:
				boxnum = self.collisiondetect(m,self.getboxes())
				
				if ( boxnum != -1):
					self.delbox(boxnum)
				return
				
			if self.moving:
				# self.moving[0] is the box, self.moving[1] are the mouse coordinates
				box = self.moving[0]
				# the old m in in self.moving[2]
				oldm = self.moving[1]
				
				self.movebox(self.moving[2],m[0]-oldm[0],m[1]-oldm[1])
				self.moving[1] = m
				self.updateAllImageDisplay()
		elif self.mmode == 1:
			m=self.guiim.scr2img((event.x(),event.y()))
			self.guiim.addShape("eraser",EMShape(["circle",.9,.9,.9,m[0],m[1],self.eraseradius,3]))
			self.boxset2.addExclusionArea("circle",m[0],m[1],self.eraseradius)
			self.updateImageDisplay()
			
	def movebox(self,boxnum,dx,dy):
		box = self.getboxes()[boxnum]
		box.xcorner += dx
		box.ycorner += dy
		box.updateBoxImage(self.image)
			# we have to update the reference also
		self.ptcl[boxnum] = box.getBoxImage(self.image)
			
		x0=box.xcorner+box.xsize/2-1
		y0=box.ycorner+box.ysize/2-1
		self.guiim.addShape("cen",EMShape(["rect",.9,.9,.4,x0,y0,x0+2,y0+2,1.0]))
		box.shape = EMShape(["rect",box.r,box.g,box.b,box.xcorner,box.ycorner,box.xcorner+box.xsize,box.ycorner+box.ysize,2.0])
		self.guiim.addShape(boxnum,box.shape)
		
	def updateboxcolors(self,classify):
		sh=self.guiim.getShapes()
		for i in classify.items():
			color = BoxingTools.get_color(i[1])
			
			sh[int(i[0])].shape[1] = color[0]
			sh[int(i[0])].shape[2] = color[1]
			sh[int(i[0])].shape[3] = color[2]
		
		self.updateImageDisplay()

	def deletenonrefs(self):
		#print event.x(),event.y(),m[0],m[1]
		# do collision detection
		boxes = self.getboxes()
		for m in range(self.boxset2.numboxes(),0,-1):
			box = boxes[m-1]
			if box.isref == False:
				self.delbox(m-1)

			#print "boxset now has",len(boxes)
	
	def mouseup(self,event) :
		if self.mmode == 0:
			m=self.guiim.scr2img((event.x(),event.y()))
			if self.moving != None:
				box = self.moving[0]
				if box.isref and self.ap: self.autobox()
			
			self.moving=None
		elif self.mmode == 1:
			self.guiim.addShape("eraser",EMShape(["circle",0,0,0,0,0,0,0.1]))
			self.updateImageDisplay()
	
	def keypress(self,event):
		if event.key() == Qt.Key_Tab:
			pass
			#self.rendcor = not self.rendcor
			#if ( self.correlation != None ):
				#if self.rendcor == True:
					#self.image.insert_clip(self.correlation,self.correlationcoords)
				#else:
					#self.image.insert_clip(self.correlationsection,self.correlationcoords)
		else: pass
	
	def updateAllImageDisplay(self):
		self.updateImageDisplay()
		self.updateMXDisplay()
	
	def updateImageDisplay(self):
		self.guiim.updateGL()
		
	def updateMXDisplay(self):
		self.guimx.updateGL()
	
	def removeshape(self,boxnum):
		print "removing shape",boxnum
		sh=self.guiim.getShapes()
		k=sh.keys()
		k.sort()
		del sh[int(boxnum)]
		for j in k:
			if isinstance(j,int):
				if j>boxnum :
					sh[j-1]=sh[j]
					del sh[j]
		self.guiim.delShapes()
		self.guiim.addShapes(sh)
		self.guiim.setActive(None,.9,.9,.4)
		
	def delbox(self,boxnum):
		"""
		Deletes the numbered box completely
		"""
		sh=self.guiim.getShapes()
		k=sh.keys()
		k.sort()
		del sh[int(boxnum)]
		for j in k:
			if isinstance(j,int):
				if j>boxnum :
					sh[j-1]=sh[j]
					del sh[j]
		self.guiim.delShapes()
		self.guiim.addShapes(sh)
		self.guiim.setActive(None,.9,.9,.4)

		self.ptcl.pop(boxnum)
		
		if self.boxset2.delbox(boxnum) and self.ap: self.autobox()
			
		
		#del(self.guimx.data[globalboxnum])
		#self.guimx.updateGL()
		#self.guimx.setData(self.ptcl)
	

	def updateboxsize(self,boxsize):
		for boxset in self.boxsets:
			boxset.updateboxsize(boxsize)
		self.boxsize=boxsize
		self.boxupdate()
	
	def boxupdate(self,force=False):
		
		ns = {}
		idx = 0
			# get the boxes
		boxes =self.getboxes()
		for j,box in enumerate(boxes):
	
			if not box.changed and not force:
				idx += 1
				continue
			
			box.changed=False
		
			im=box.getBoxImage(self.image)
			box.shape = EMShape(["rect",box.r,box.g,box.b,box.xcorner,box.ycorner,box.xcorner+box.xsize,box.ycorner+box.ysize,2.0])
			if not box.isref:
				box.shape.isanimated = True
				box.shape.blend = 0
			ns[idx]=box.shape
			if idx>=len(self.ptcl) : self.ptcl.append(im)
			else : self.ptcl[idx]=im
			idx += 1
		

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
		
		return (self.getboxes(),self.threshold)
	
	def getrefimages(self):
		images = []
		global_box_num = 0
		# basic strategy is to detect if there was a collision with any other box
		# do this by incrementing through all the box sets
		for boxset in self.boxsets:
			# get the boxes
			boxes = boxset.boxes
			
			#print event.x(),event.y(),m[0],m[1]
			# do collision detection
			for boxnum,box in enumerate(boxes):
				if box[4] == 1:
					images.append(self.ptcl[global_box_num])
				
				global_box_num += 1
			#print "there are",len(images),"refs"
			return images
	def updatetemplate(self):
		self.boxset2.updatetemplate(self.boxsize)
			
	def autoboxbutton(self):
		self.autobox()
		
	def autobox(self):
		if self.boxset2.templateupdate:
			self.deletenonrefs()
		self.updatetemplate()
		self.boxset2.autobox()

	def dynapick(self):
		self.ap = not self.ap
		
	def done(self):
		self.app.quit
		
	def trydata(self,data,thr):
		self.deletenonrefs()
		
		self.boxset2.autobox()
		
	def updatedata(self,data,thresh):
		self.guictl.updatedata(data,thresh)
		
	def setautobox(self,s):
		self.boxset2.autoboxmethod = s
	def nocupdate(self,bool):
		self.boxset2.templateupdate = not bool	
	def classify(self,bool):
		self.boxset2.gen_ref_images()
		
	def erasetoggled(self,bool):
		# for the time being there are only two mouse modes
		self.mmode = bool
		
	def updateEraseRad(self,rad):
		self.eraseradius = rad

class GUIboxPanel(QtGui.QWidget):
	def __init__(self,target) :
		
		QtGui.QWidget.__init__(self,None)
		self.target=target
		
		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		
		self.plothbl = QtGui.QHBoxLayout()
		
		self.window = EMGLPlotWidget(self)
		self.window.setInit()
		self.window.resize(100,100)
		self.window2=EMParentWin(self.window)
		self.window2.resize(100,100)
		
		self.plothbl.addWidget(self.window2)
		
		self.plotbuttonvbl = QtGui.QVBoxLayout()
		
		self.trythat=QtGui.QPushButton("Try That")
		self.plotbuttonvbl.addWidget(self.trythat)
		
		self.reset=QtGui.QPushButton("Reset")
		self.plotbuttonvbl.addWidget(self.reset)
		
		self.plothbl.addLayout(self.plotbuttonvbl)
		
		self.vbl2 = QtGui.QVBoxLayout()
		
		self.vbl2.addLayout(self.plothbl)
		
		self.thr = ValSlider(self,(0.0,3.0),"Threshold:")
		self.thr.setValue(target.threshold)
		self.vbl2.addWidget(self.thr)
		
		self.interbox = QtGui.QGroupBox("Interactive Parameters")
		self.interbox.setLayout(self.vbl2)
		self.vbl.addWidget(self.interbox)
		
		self.thrbut = QtGui.QRadioButton(THRESHOLD)
		self.selbut = QtGui.QRadioButton(SELECTIVE)
		self.selbut.setChecked(True)
		self.morselbut = QtGui.QRadioButton(MORESELECTIVE)
		
		self.methodhbox = QtGui.QHBoxLayout()
		self.methodhbox.addWidget(self.thrbut)
		self.methodhbox.addWidget(self.selbut)
		self.methodhbox.addWidget(self.morselbut)
		
		self.groupbox = QtGui.QGroupBox("Auto Box Method")
		self.groupbox.setLayout(self.methodhbox)
		
		self.vbl.addWidget(self.groupbox)
		
		#self.vbl.addLayout(self.groupbox)
		
		self.infohbl = QtGui.QHBoxLayout()
		self.info = QtGui.QLabel("%d Boxes"%len(target.getboxes()),self)
		self.ppc = QtGui.QLabel("%f particles per click"%0,self)
		self.infohbl.addWidget(self.info)
		self.infohbl.addWidget(self.ppc)
		
		
		self.statsbox = QtGui.QGroupBox("Stats")
		self.statsbox.setLayout(self.infohbl)
		self.vbl.addWidget(self.statsbox)
		
		self.hbl1=QtGui.QHBoxLayout()
		self.hbl1.setMargin(0)
		self.hbl1.setSpacing(2)
		#self.vbl.addLayout(self.hbl1)
		
		self.lblbs=QtGui.QLabel("Box Size:",self)
		self.hbl1.addWidget(self.lblbs)
		
		self.bs = QtGui.QLineEdit(str(target.boxsize),self)
		self.hbl1.addWidget(self.bs)
		self.vbl.addLayout(self.hbl1)
		
		self.hbl3=QtGui.QHBoxLayout()
		self.dynapick = QtGui.QRadioButton("Dynapix")
		self.hbl3.addWidget(self.dynapick)
		self.nocpick = QtGui.QRadioButton("No correlation update")
		self.hbl3.addWidget(self.nocpick)
		
		self.vbl.addLayout(self.hbl3)
		
		self.autobox=QtGui.QPushButton("Auto Box")
		self.vbl.addWidget(self.autobox)


		self.hbl2=QtGui.QHBoxLayout()
		self.hbl2.setMargin(2)
		self.hbl2.setSpacing(6)
		#self.vbl.addLayout(self.hbl1)
		
		
		#try:
		self.erasepic = QtGui.QIcon("/home/d.woolford/erase.png");
		self.erase=QtGui.QPushButton(self.erasepic,"Erase")
		self.erase.setCheckable(1)
		self.hbl2.addWidget(self.erase)
		
		self.eraseradtext=QtGui.QLabel("Circle radius",self)
		self.hbl2.addWidget(self.eraseradtext)
		
		self.eraserad = QtGui.QLineEdit(str(self.target.eraseradius),self)
		self.hbl2.addWidget(self.eraserad)
		self.eraserad.setEnabled(False)
		
		self.vbl.addLayout(self.hbl2)


		self.done=QtGui.QPushButton("Done")
		self.vbl.addWidget(self.done)

		self.classifybut=QtGui.QPushButton("Classify")
		self.vbl.addWidget(self.classifybut)
		#except: pass
		
		self.connect(self.bs,QtCore.SIGNAL("editingFinished()"),self.newBoxSize)
		self.connect(self.eraserad,QtCore.SIGNAL("editingFinished()"),self.updateEraseRad)
		self.connect(self.thr,QtCore.SIGNAL("valueChanged"),self.newThresh)
		self.connect(self.done,QtCore.SIGNAL("clicked(bool)"),self.target.app.quit)
		self.connect(self.classifybut,QtCore.SIGNAL("clicked(bool)"),self.target.classify)
		self.connect(self.autobox,QtCore.SIGNAL("clicked(bool)"),self.target.autoboxbutton)
		self.connect(self.dynapick,QtCore.SIGNAL("clicked(bool)"),self.dynapickd)
		self.connect(self.trythat,QtCore.SIGNAL("clicked(bool)"),self.trythatd)
		self.connect(self.thrbut, QtCore.SIGNAL("clicked(bool)"), self.gboxclick)
		self.connect(self.selbut, QtCore.SIGNAL("clicked(bool)"), self.gboxclick)
		self.connect(self.morselbut, QtCore.SIGNAL("clicked(bool)"), self.gboxclick)
		self.connect(self.nocpick, QtCore.SIGNAL("clicked(bool)"), self.target.nocupdate)
		self.connect(self.erase, QtCore.SIGNAL("clicked(bool)"), self.erasetoggled)
#		self.target.connect(self.target,QtCore.SIGNAL("nboxes"),self.nboxesChanged)
	
	def erasetoggled(self,bool):
		self.eraserad.setEnabled(bool)
		self.target.erasetoggled(bool)
	
	def dynapickd(self,bool):
		if bool == True:
			self.autobox.setEnabled(False)
		else:
			self.autobox.setEnabled(True)
		self.target.dynapick()
	
	def gboxclick(self,bool):
		if self.thrbut.isChecked():
			s = self.thrbut.text()
		elif self.selbut.isChecked():
			s = self.selbut.text()
		elif self.morselbut.isChecked():
			s = self.morselbut.text()
		else:
			print "Bug intercepted in e2boxer.py. Please email the development team."
			
		self.target.setautobox(str(s))
	
	def trythatd(self):
		self.target.trydata(self.window.getData(),float(self.thr.getValue()))
	
	def updatedata(self,data,thresh):
		#print data
		self.window.setData(data)
		self.thr.setValue(thresh,True)
		self.resize(self.width(),self.height())
		#self.window.resizeGL(self.window.width(),self.window.height())
		#self.window.updateGL()
	def nboxesChanged(self,n):
		self.info.setText("%d Boxes"%n)
		
	def ppcChanged(self,f):
		self.ppc.setText("%f ppc"%f)
	
	def updateEraseRad(self):
		v = int(self.eraserad.text())
		if ( v < 1 ): raise Exception
		self.target.updateEraseRad(v)
	
	def newBoxSize(self):
		try:
			v=int(self.bs.text())
			if v<12 : raise Exception
		except:
			self.bs.setText(str(self.target.boxsize))
			return
		
		self.target.updateboxsize(v)
		
	def newThresh(self,val):
		#print "new threshold"
		self.trythatd()


if __name__ == "__main__":
	main()
