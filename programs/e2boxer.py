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
from emimagemx import EMImageMX
from math import *
from time import *
import os
import sys

from copy import *

from emglplot import *

from shelve import open

from time import time

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
	
	
	image.process_inplace("normalize.edgemean")

		
	if len(options.auto)>0 : print "Autobox mode ",options.auto[0]
	
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
		self.isanchor = True		# A flag used by AutoBoxer routines that - if set to true the box will not be included in the generation the template) - This is specific to the SwarmPS autoboxer
		self.boxingobj = None

	def updateBoxImage(self,image,norm=True):
		#print "getting region",self.xcorner,self.ycorner,self.xsize,self.ysize
		self.image = image.get_clip(Region(self.xcorner,self.ycorner,self.xsize,self.ysize))
		if norm:
			self.image.process_inplace("normalize.edgemean")
		
		# make sure there are no out of date footprints hanging around
		self.footprint = None
		
			
	def getBoxImage(self,image=None,norm=True,force=False):
		if self.image == None or force:
			if image == None:
				print 'error, need to specify the image argument when first calling getBoxImage'
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
			
	def updateParams(self,autoBoxer,center=True):
		'''
		Updates internally stored parameters, currently works only for SwarmAutoBoxer, but
		have attempted to lay basic framework if in future we use a different autoBoxer which
		requires its own parameters
		'''
		
		correlation = self.boxingobj.getCorrelation()
		if correlation == None:
			print 'error, can not update the parameters of a Box because the Boxable has no correlation image'
			return 0
		
		if isinstance(autoBoxer,SwarmAutoBoxer):
			shrink = autoBoxer.getBestShrink()
			invshrink = 1/shrink
	
			# the central coordinates of the box in terms of the shrunken correlation image
			x = (self.xcorner+self.xsize/2.0)*invshrink
			y = (self.ycorner+self.ysize/2.0)*invshrink
			
			#the search radius is used in correlation space - it limits the radial distance
			# up to which 'profile' data can be accrued
			# it is currently half the boxsize in terms of the correlation image's dimensions
			searchradius = autoBoxer.getSearchRadius()
		
			peak_location = BoxingTools.find_radial_max(correlation,int(x),int(y),searchradius )
			peak_location2 = BoxingTools.find_radial_max(correlation,peak_location[0],peak_location[1],searchradius )
			if (peak_location != peak_location2):
				# this represents a troubling condition
				# setting box.correlationscore is the flag that other functions can act on in order to exclude
				# this box from consideration
				self.correlationscore = None
				print "Error, peak location unrefined"
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
			print 'error, the autoBoxer you are using is not currently known by the Box class'
			return 0
class Reference(Box):
	'''
	A reference is a box, but with extra capabilities. It knows the Boxable from which it originated.
	It can therefore tell 
	'''
	def __init__(self):
		Box.__init__(self)
		self.boxingobj = None
	

class Boxable:
	def __init__(self,image,parent=None,autoBoxer=None):
		self.image = image			# the image containing the boxes
		self.parent = parent		# keep track of the parent in case we ever need it
		self.boxes = []				# a list of boxes
		self.refboxes = []			# a list of boxes
		self.boxsize = -1			#  the boxsize
		self.smallimage = None		# a small copy of the image which has had its background flattened
		self.flattenimager = -1		# the r value used to run the flatten image processor
		
		self.fpshrink = -1
		self.exclusionimage = None
		self.template = None
		self.correlation = None
		self.refcache = []
		self.allowcorrelationupdate = False	# a temporary flag that can be used by externally objects, for instance a reference box, which is forcing an update for example
		self.templateTS = -1 # a template time stamp, used to avoid unecessarily regenerating the template in self.autoBox
		self.autoBoxerTS = -1 # and autoBoxer time stamp, used to avoid unecessary autoboxing, and to force autoboxing when appropriate
		
		self.autoBoxer = autoBoxer
	
	def setAutoBoxer(self,autoBoxer):
		self.autoBoxer = autoBoxer
		
	def getCorrelation(self):
		return self.correlation
	
	def addbox(self,box):
		if not isinstance(box,Box):
			print "You can not add a box to this box set if it is not of type Box"
			return;
		
		box.isref = True # make sure it knows that it's a reference box
		box.boxingobj = self
		box.otheridx1 = len(self.boxes) # store this, it's used when the user deletes the box, for exclusion
		box.otheridx2 = len(self.refboxes) # store this, it's used when the user deletes the box, for exclusion
		
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
	
	def updateBoxSize(self,boxsize):
		'''
		Updates only the box size and corner coordinates
		Switches the changed flag to True to trigger redisplay (but the calling function
		is responsible for knowing and testing for this)
		'''
		# do nothing if it's the same size as what we already have
		if  boxsize == self.boxsize: return
		
		for box in self.boxes:
			if box.xsize != boxsize:
				box.xcorner -= (boxsize-box.xsize)/2
				box.xsize = boxsize
				box.changed = True
			if box.ysize != boxsize:
				box.ycorner -= (boxsize-box.ysize)/2
				box.ysize = boxsize
				box.changed = True
			
			box.image = None
			box.footprint = None

		self.fprink = -1
		self.flattenimager = -1
		self.boxsize = boxsize
		self.smallimage = None
		self.correlation = None
		
	def getfootprintshrink(self):
		if self.fpshrink == -1:
			shrink = 1
			tn = self.boxsize/2
			while ( tn >= 32 ):
				tn /= 2
				shrink *= 2
			self.fpshrink = shrink
		
		return self.fpshrink
		
	def getBestShrink(self):
		'''
		FIXME - there should probably be a more well established framework for doing this
		At the moment it is possible that the self.autBoxer is actually None, which isn't good.	
		'''
		if self.autoBoxer != None:
			return self.autoBoxer.getBestShrink()
		else:
			print 'warning, there is not autoboxer set, am not sure how to shrink, returning 1 as the shrink factor'
			return 1
		
	def updateCorrelation(self,template,templateTS,autoBoxer):
		'''
		A function that will update the correlation image if the correlationupdate flag is set to true
		Useful if a template has been updated somewhere, yet many references originate from this BoxingOject -
		All the references will call this function, but we only need to act the first time it happens
		
		Be warned - whoever opens the gate by setting self.allowcorrelationupdate to True should set it 
		to False once their done with the specialized operation
		
		'''
		if self.allowcorrelationupdate:
			self.templateTS = templateTS # Time Stamp, used for efficiency in autoBox to save an unecessary correlation update
			self.__genCorrelation(template,autoBoxer)
			
			# I made a conscientious decision to leave the responsibility of turning this flag off
			# to that of the calling program/function. This uncommented line is left only for documentation purposes
			#self.allowcorrelationupdate = False

	def __genCorrelation(self,template,autoBoxer,forceupdate=False):
		'''
		The force update flag is only meant to be used if the box size has changed - this changes 
		the shrink factor, and also affects the background flattening process.
		'''
		self.template = template
		newr = self.template.get_xsize()/2.0
		# now we only recalculate the small copy of the subject image if necessary
		if self.smallimage == None or newr != self.flattenimager or forceupdate:
			self.flattenimager = newr
			#FIXME - we could avoid a deep copy by writing the meanshrink processor
			# i.e. section = self.image.process("math.meanshrink",{"n":self.getBestShrink()}
			
			if (self.getBestShrink() != 1):
				self.smallimage = self.image.process("math.meanshrink",{"n":autoBoxer.getBestShrink()})
			else: self.smallimage = self.image.copy()
			
			self.smallimage.process_inplace("filter.flattenbackground",{"radius":self.flattenimager})
			print "generated shrunken image" # DEBUG

		self.correlation = self.smallimage.calc_flcf( self.template )
		
		# this may not be necessary if we ever want to be completely efficient
		self.correlation.process_inplace("xform.phaseorigin.tocenter")
		#self.correlation.write_image("tttttt.hdf")
		print "generated correlation image" # DEBUG
		
		
	#def autobox(self,autoBoxer,optprofile=None,thr=None):
		#'''
		#FIXME - act on the argument optprofile and thr
		#'''
		## sometimes the template would have been updated transparently via a reference update in the
		## main interface.
		
		#print 'called autobox' # DEBUG
		
		#if self.templateTS == -1 or self.correlation == None or self.templateTS != autoBoxer.templateTS:
			#print 'regenerating the correlation image' # DEBUG
			#template = autoBoxer.getTemplate()
			#self.templateTS = autoBoxer.templateTS
			#if template != None: self.__genCorrelation(template,autoBoxer)
			#else:
				#print 'error, cant ask the autoBoxer for its template, it doesnt have one'
				#return 0
		#else : print 'using cached correlation image' # DEBUG

		## auto boxing will only ever occur if the time stamp of the AutoBoxer is not the
		## same as the time stamp cached by this object. -1 means it's the first time.
		#if self.autoBoxerTS == -1 or self.autoBoxerTS != autoBoxer.stateTS:
			#print 'auto boxing' # DEBUG
			#exclusion = self.getExclusionImage().copy()
			
			## paint white exlclusion circles around the already selected and reference boxes
			## this will save the autoBox time and prevent redundant boxing
			#self.updateExclusion(exclusion)
		
			#boxes = autoBoxer.autoBox(self.correlation,self.boxes,exclusion)
			#self.autoBoxerTS = autoBoxer.stateTS

			## autoBoxer should return 0 if there was a problem
			#if boxes != 0:
				#for box in boxes: self.boxes.append(box)

			#self.parent.boxupdate()
		#else: print 'no auto boxing was necessary, up-2-date' # DEBUG

	def updateExcludedBoxes(self):
		lostboxes = []
		
		invshrink = 1.0/self.getBestShrink()
		exc = self.getExclusionImage()
		n = len(self.boxes)
		for i in range(n-1,-1,-1):
			box = self.boxes[i]
			x = int((box.xcorner+box.xsize/2.0)*invshrink)
			y = int((box.ycorner+box.ysize/2.0)*invshrink)
			
			if ( exc.get(x,y) != 0):
				lostboxes.append(i)
	
		return lostboxes
	
	def addExclusionArea(self, type,x,y,radius):
		
		xx = int(x/self.getBestShrink())
		yy = int(y/self.getBestShrink())
		
		rr = int(radius/self.getBestShrink())
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
				
				self.exclusionimage.set(ii,jj,0.1)
				
		# FIXME - is this the best place for this?
		self.parent.guiim.setOtherData(self.getExclusionImage(),self.getBestShrink(),True)
	
	def getExclusionImage(self):
		if self.exclusionimage == None:
			self.exclusionimage = EMData(self.correlation.get_xsize(),self.correlation.get_ysize())
			self.exclusionimage.to_zero()
		return self.exclusionimage
	
	def updateExclusion(self,exclusion):
		'''
		paints black circles in the exclusion - which is a binary EMData object
		useful for making things efficient, should probably be called updateExclusions
		'''
		
		if self.autoBoxer != None:
			searchradius = self.autoBoxer.getSearchRadius()
		else:
			print 'warning, there is no autoboxer. Am unsure of how to update excluded zones - doing nothing. This may result in strange autoboxing results'
			return
		
		for box in self.boxes:
			xx = box.xcorner + box.xsize/2
			yy = box.ycorner + box.ysize/2
			xx /= self.getBestShrink()
			yy /= self.getBestShrink()
			
			BoxingTools.set_radial_non_zero(exclusion,int(xx),int(yy),searchradius)
			
		
	#def classify(self):
		#v = []
		## accrue all params
		#for box in self.boxes:
		#self.accrueparams(self.boxes)
		
		#for box in self.boxes:
			#b = copy(box.optprofile[0:self.radius])
			#b.sort()
			##for a in b:
				##a = box[6]-a
			##print b
			#v.append(b)
			
		#cl = BoxingTools.classify(v,4)
		#self.parent.updateboxcolors(cl)
	
	def genRefImages(self):
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
			if self.getBestShrink() != 1:
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

class AutoBoxer:
	'''
	Base class design for auto boxers
	'''
	def __init__(self,boxable=None):
		self.version = 1.0

	def getTemplate(self):
		'''This should return a single template which is an EMData object'''
		raise Exception
	
	def name(self):
		'''
		Every autoboxer should return a unique name
		'''
		raise Exception
	
	def addReference(self,box):
		'''
		add a reference box - the box should be in the format of a Box object, see above
		Returns 0 if there is a problem, returns 1 if it's all good
		Adds a reference to a list
		'''
		raise Exception
	
	def removeReference(self,box):
		'''
		Remove a reference box - the box should in the format of a Box object, see above
		Pops a reference from a list
		'''
		raise Exception
	
	def referenceMoved(self,ref):
		'''
		If a reference was moved interactively in the interface this is the function that should be called
		'''
		raise Exception

	def getTemplate(self):
		'''
		Return the template that is being used. Returns None if there is not template
		'''
		raise Exception

	def setBoxSize(self,boxsize):
		'''
		Hard set the boxsize. Note that nothing is done to the reference boxes. It is
		assumed whichever part of the program calls this function also updates the Box objects
		independently (which implicitly affects the boxes stored internally in the AutoBoxer
		class, because it only ever stores programmatic references)
		'''
		raise Exception
	
	def autoBox(self,correlation,boxes=[],exclusion=None):
		'''
		The main autoBox routine. The calling program should pass in its own correlation map (EMData), and optionally
		an exclusion map of ones and zeros (0 means include, non zero means exclude). Also a list of boxes that will
		not be removed prior to the autoBoxing (and hence probably have associated excluded regions in the exlcusion
		image, but this is beside the point), The model of use here is that
		the calling program should get the current template from the AutoBoxer to generate the correlation map. The
		calling program should be able to cache the correlation map, so should be able to detect if there's been
		a template update by asking for the current set of references (getReferences) and cross checking against a list of its own.
		@Returns a list of Boxes
		'''
		raise Exception
	
class SwarmAutoBoxer(AutoBoxer):
	'''
	This is an autoboxer that encapsulates the boxing approach first developed in SwarmPS
	'''
	DYNAPIX = 1
	ANCHOREDDYNAPIX = 2
	USERDRIVEN = 3
	def __init__(self,boxable):
		AutoBoxer.__init__(self,boxable)
		self.refboxes = []		# this will eventually be a list of Box objects
		self.template = None	# an EMData object that is the template
		self.boxsize = -1		# stores the global boxsize, this is the value being used by boxer in the main interface
		self.shrink = -1
		
		# more privately stuff
		self.templatedimmin = 24  # the smallest amount the template can be shrunken to. Will attempt to get as close to as possible. This is an important part of speeding things up.
		self.optthreshold = -1	# the correlation threshold, used to as the basis of finding local maxima
		self.optprofile = []	# the optimum correlation profile used as the basis of auto selection
		self.optprofileradius = -1 # the optimum radius - used to choose which part of the optprofile is used as the basis of selection
		self.autoboxmethod = SELECTIVE	# the autobox method - see EMData::BoxingTools for more details
		self.__shrink = -1
		
		self.templateTS = -1 # a template time stamp to 
		self.stateTS = -1 # the time stamp that records when the current template and parameters are completely up to date
		
		self.boxable = boxable
		
		self.mode = SwarmAutoBoxer.DYNAPIX
		self.refupdate = False # this is a flag used when self.mode is USERDRIVEN
		self.permissablemodes = [SwarmAutoBoxer.DYNAPIX,SwarmAutoBoxer.ANCHOREDDYNAPIX,SwarmAutoBoxer.USERDRIVEN]  # if another mode is added you would have to find all places where self.mode is used to make decisions and alter
	
	def setBoxable(self,boxable):
		self.boxable = boxable
	
	def setMode(self,mode):
		if mode in self.permissablemodes: self.mode = mode
		else: print 'error, unknow mode:',mode,'- no action taken'
		
	def name(self):
		return 'swarmautoboxer'

	def addReference(self,box):
		'''
		 add a reference box - the box should be in the format of a Box, see above):
		'''
		if isinstance(box,Box):
			if box.xsize != box.ysize:
				print 'error, support for uneven box dimensions is not currently implemented'
				return 0
		
			# store the boxsize if we don't have one already
			if self.boxsize == -1:
				self.boxsize = box.xsize
			# do a sanity check, this shouldn't happen if the program is managing everything carefully
			elif self.boxsize != box.xsize:
				print 'error, the currently stored box size does not match the boxsize of the reference that was just added'
				return 0
			
			self.refboxes.append(box)
		
			if self.mode == SwarmAutoBoxer.DYNAPIX:
				if not box.isanchor:
					print 'the box flag is internally inconsistent when using pure dynapix'
					return 0
				self.__fullUpdate()
				self.autoBox(self.boxable)
			elif self.mode == SwarmAutoBoxer.ANCHOREDDYNAPIX:
				if box.isanchor:
					print 'the box flag is internally inconsistent when anchoring'
					return 0
				box.updateParams(self)
				self.__accrueOptParams()
				self.stateTS = time()
				self.autoBox(self.boxable)
			elif self.mode == SwarmAutoBoxer.USERDRIVEN:
				self.refupdate = True
				self.stateTS = -1
				self.templateTS = -1
			else:
				print 'error, unknown mode in SwarmAutoBoxer'
				return 0
			
			return 1
			
		else:
			print "error, you cannot add a reference to the AutoBoxer if it is not in the format of a Box object"
			return 0
	
	def removeReference(self,box):
		for j,tmp in enumerate(self.refboxes):
			if box.isref and box.xcorner == tmp.xcorner and box.ycorner == tmp.ycorner:
				tmp = self.refboxes.pop(j)
				if len(self.refboxes) == 0:
					self.boxsize = -1
					
				if self.mode == SwarmAutoBoxer.DYNAPIX:
					if not box.isanchor:
						print 'the box flag is internally inconsistent when using pure dynapix'
						return 0
					self.__fullUpdate()
					self.autoBox(self.boxable)
				elif self.mode == SwarmAutoBoxer.ANCHOREDDYNAPIX:
					if box.isanchor:
						print 'the box flag is internally inconsistent when anchoring'
						return 0
					box.updateParams(self)
					self.__accrueOptParams()
					self.stateTS = time()
					self.autoBox(self.boxable)
				elif self.mode == SwarmAutoBoxer.USERDRIVEN:
					self.refupdate = True
					self.stateTS = -1
					self.templateTS = -1
				else:
					print 'error, unknown mode in SwarmAutoBoxer'
				
				return True
		
		return False
	
	def referenceMoved(self,ref):
		'''
		If a reference was moved interactively in the interface this is the function that should be called
		'''
		if self.mode == SwarmAutoBoxer.DYNAPIX:
			if not box.isanchor:
				print 'the box flag is internally inconsistent when using pure dynapix'
				return 0
			self.__fullUpdate()
			self.autoBox(self.boxable)
		elif self.mode == SwarmAutoBoxer.ANCHOREDDYNAPIX:
			if box.isanchor:
				print 'the box flag is internally inconsistent when anchoring'
				return 0
			box.updateParams(self)
			self.__accrueOptParams()
			self.stateTS = time()
			self.autoBox(self.boxable)
		elif self.mode == SwarmAutoBoxer.USERDRIVEN:
			self.refupdate = True
			self.stateTS = -1
			self.templateTS = -1
		else:
			print 'error, unknown mode in SwarmAutoBoxer'
		
	def getTemplate(self):
		if self.mode == SwarmAutoBoxer.USERDRIVEN:
			# If it's user driven then the user has selected a bunch of references and then hit 'autobox'.
			# In which case we do a complete reference update, which generates the template and the
			# best autoboxing parameters
			if len(self.refboxes) == 0:
				print 'error, cant get template if there are no references'
				return None
			
			if self.refupdate:
				self.__fullUpdate()
				self.refupdate = False
				
		elif self.mode == SwarmAutoBoxer.DYNAPIX or self.mode == SwarmAutoBoxer.ANCHOREDDYNAPIX:
			
			if self.template == None:
				print 'error, you have either asked for the template without setting a reference, or you have added a reference and not set the refupdate flag'
				return None
		
		return self.template
		
	def setBoxSize(self,boxsize):
		if (boxsize < 6 ):
			print 'error, a hard limit of 6 for the box size is currently enforced. Email developers if this is a problem'
			return
		if self.boxsize == boxsize:	return
		
		self.boxsize = boxsize
		# make sure the shrink value is updated - use the force flag to do it
		self.getBestShrink(True)
		
		if self.mode == SwarmAutoBoxer.DYNAPIX or self.mode == SwarmAutoBoxer.ANCHOREDDYNAPIX:
			# update references
			self.__fullUpdate()
			self.autoBox(self.boxable)
		elif self.mode == SwarmAutoBoxer.USERDRIVEN:
			self.refupdate = True
			self.stateTS = -1
			self.templateTS = -1
		else:
			print 'error, unknown mode in SwarmAutoBoxer'
	
	def getSearchRadius(self):
		return int(0.5*(self.boxsize)/float(self.getBestShrink()))
	
	
	def getBestShrink(self,force=True):	
		if self.boxsize == -1:	
			print "error - the boxsize is currently -1 - I can't figure out the best value to shrink by"	
			return -1
			
		if self.shrink == -1 or force:	
			self.shrink = ceil(float(self.boxsize)/float(self.templatedimmin))	
			
		return self.shrink
		
	def autoBox(self,boxable):
		if self.mode == SwarmAutoBoxer.USERDRIVEN:
			# If it's user driven then the user has selected a bunch of references and then hit 'autobox'.
			# In which case we do a complete reference update, which generates the template and the
			# best autoboxing parameters
			if len(self.refboxes) == 0:
				print 'error, cant get template if there are no references'
				return 0
			
			if self.refupdate:
				self.__fullUpdate()
				self.refupdate = False
				
				boxable.allowcorrelationupdate = True
				boxable.updateCorrelation(self.template,self.templateTS,self)
				boxable.allowcorrelationupdate = False
				
				correlation = boxable.correlation
				
				# FIXME, this is a copy paste from below
				exclusion = boxable.getExclusionImage().copy()
				# paint white exlclusion circles around the already selected and reference boxes
				# this will save the autoBox time and prevent redundant boxing
				boxable.updateExclusion(exclusion)
				
				boxes = self.__autoBox(correlation,boxable.boxes,exclusion)
				boxable.autoBoxerTS = self.stateTS
				
				# This shouldn't happen in the Database instance
				if boxes != 0:
					for box in boxes: boxable.boxes.append(box)
					
				return 1
					
				# end FIXME
			else: 
				# the user has not triggered any update flags
				# this probably means the user has hit the autobox button twice in a row, or somesuch
				return 0
		elif self.mode == SwarmAutoBoxer.DYNAPIX or self.mode == SwarmAutoBoxer.ANCHOREDDYNAPIX:
			templateTS = boxable.templateTS
			correlation = boxable.correlation
		
			if templateTS == -1 or correlation == None or self.templateTS != templateTS:
				print 'regenerating the correlation image' # DEBUG
				boxable.allowcorrelationupdate = True
				if self.template != None:
					boxable.allowcorrelationupdate = True
					print 'told the boxable to update its template'
					boxable.updateCorrelation(self.template,self.templateTS,self)
					boxable.allowcorrelationupdate = False
					
					correlation = boxable.correlation
				else:
					print 'error, cant ask the autoBoxer for its template, it doesnt have one'
					return 0
			else : print 'using cached correlation image' # DEBUG

			autoBoxerTS = boxable.autoBoxerTS
			# auto boxing will only ever occur if the time stamp of the AutoBoxer is not the
			# same as the time stamp cached by the Boxable. -1 means it's the first time.
			if autoBoxerTS == -1 or autoBoxerTS != self.stateTS:
				print 'auto boxing' # DEBUG
				exclusion = boxable.getExclusionImage().copy()
				
				if self.mode == SwarmAutoBoxer.ANCHOREDDYNAPIX:
					# we should remove non references here
					pass
					
				# paint white exlclusion circles around the already selected and reference boxes
				# this will save the autoBox time and prevent redundant boxing
				boxable.updateExclusion(exclusion)
			
				boxes = self.__autoBox(correlation,boxable.boxes,exclusion)
				boxable.autoBoxerTS = self.stateTS
	
				# This shouldn't happen in the Database instance
				if boxes != 0:
					for box in boxes: boxable.boxes.append(box)
					
				return 1

			else: print 'no auto boxing was necessary, up-2-date' # DEBUG
			
	def __autoBox(self,correlation,boxes=[],exclusion=None):
		'''
		Does the autoboxing. Returns a list of Boxes
		'''
		if not isinstance(correlation,EMData):
			print 'error, cannot autobox, the correlation argument is not an EMData object'
			return 0
			
			#print "using opt radius",self.radius, "which has value",tmp,"shrink was",self.shrink
		if self.autoboxmethod == THRESHOLD:
			mode = 0
		elif self.autoboxmethod == SELECTIVE:
			mode = 1
		elif self.autoboxmethod == MORESELECTIVE:
			mode = 2
		
		shrink = self.getBestShrink()
		# Warning, this search radius value should be the same as the one used by the BoxSets that contributed the reference boxes
		# to this AutoBoxer object. There should be one place/function in the code where both parties access this value
		searchradius = self.getSearchRadius()
		
		soln = BoxingTools.auto_correlation_pick(correlation,self.optthreshold,searchradius,self.optprofile,exclusion,self.optprofileradius,mode)

		# This is what should be written to the database
		boxes = []
		
		for b in soln:
			x = b[0]
			y = b[1]
			xx = int(x*shrink)
			yy = int(y*shrink)
			box = Box(xx-self.boxsize/2,yy-self.boxsize/2,self.boxsize,self.boxsize,0)
			box.correlationscore = correlation.get(x,y)
			box.corx = b[0]
			box.cory = b[1]
			box.changed = True
			boxes.append(box)
	
		return boxes
		
	#  PRIVATE
	#  SWARMAUTOBOXER PRIVATE FUNCTIONS
	#  PRIVATE
	def __genTemplate(self):
		'''
		Returns 0 if there are errors
		Return 1 if not
		'''
		# you can only generate a template if there are references
		if len(self.refboxes) <= 0: 
			print 'error, cant call private function genTemplate when there are no refboxes, this is an internal error'
			return 0
		
		if self.boxable.image == None:
			print 'error, need the big image to be set when generating template so the references can be extracted from it, if need be'
			return 0
		
		print "using",len(self.refboxes),"references"
		images_copy = []
		for ref in self.refboxes:
			# some references can be excluded from the template generation procedure, this is flagged
			# by the isanchor flag
			if ref.isanchor == False:
				continue
			image = ref.getBoxImage(self.boxable.image)
			if (self.getBestShrink() != 1):
				e = image.process("math.meanshrink",{"n":self.getBestShrink()})
			else : e = image.copy()
			images_copy.append(e)
			
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
			ave.process_inplace("xform.centeracf")
			ave.process_inplace("mask.sharp",{'outer_radius':ave.get_xsize()/2})
		
		# debug, un-comment to see the aligned refs and the final template
		#for image in t:
			#image.write_image("aligned_refs.img",-1)
		
		#ave.write_image("aligned_refs.img",-1)
		
		#black = EMData(image.get_xsize(),image.get_ysize())
		#black.to_zero()
		#black.write_image("aligned_refs.img",-1)
		
		self.template = ave
		self.templateTS = time()
		return 1
	
	def __fullUpdate(self):
		'''
		Forces a template update, then updates all correlation images
		that the references come from, then does a parameters update.
		This is like a completeness function - this needs to happen for internal
		consistency. It is motivated by the thought that if references
		come from many images, then the correlation images of each of the images
		needs to be updated in order for the correlation parameters to be generated
		consistently (from a correlation image generated by a universal template)
		'''
		
		if not self.__genTemplate():
			print 'error, couldnt generate template'
			return 0
		
		# First tell all references' associated boxing objects to be open to the prospect 
		# if update their correlation images
		for ref in self.refboxes:
			ref.boxingobj.allowcorrelationupdate = True
			
		# Now iterate through and force the correlation update
		# the boxing Boxables 
		for ref in self.refboxes:
			ref.boxingobj.updateCorrelation(self.template,self.templateTS,self)
			ref.boxingobj.allowcorrelationupdate = False
			
		for ref in self.refboxes:
			ref.updateParams(self)
	
		# parameters should be updated now
		# it's important that the BoxingObjext.updateCorrelation updated the parameters stored in the boxes
		self.__accrueOptParams()
		
		self.stateTS = time()

	
	def __accrueOptParams(self):
		'''
		A function for accruing the parameters of the SwarmPSAutoBoxer autoboxing technique
		returns True if optimal parameters were accrued
		return False if not
		'''

		# To determine the threshold from what we've got, iterate through all of the reference
		# boxes and use the lowest correlation score as the correlation threshold
		found = False
		for i,box in enumerate(self.refboxes):
			if box.correlationscore == None:
				# this is an error which probably means that the box, as created by the user, has a strong correlation maximum next to it which is disrupting the auto parameters
				# this is mostly an error for dwoolfords attention
				# for the time being just ignoring it  probably suffices
				# FIXME
				print "continuing on faulty"
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
		
	
		# determine the point in the profile where the drop in correlation score is the greatest, store it in radius
		self.optprofileradius = -1
		tmp = self.optprofile[0]
		for i in range(1,len(self.optprofile)):
			# the tmp > 0 is a
			if self.optprofile[i] > tmp and tmp > 0:
				tmp = self.optprofile[i]
				self.optprofileradius = i
				
		return True
		
	
class GUIbox:
	def __init__(self,imagefsp,boxes,thr,boxsize=-1):
		"""Implements the 'boxer' GUI. image is the entire image, and boxes and thr specify current boxes
		to begin with. Modified boxes/thr are returned. 
		
		'boxes' is a list of box locations:
		[x0,y0,xsize,ysize,quality,changed]
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
		
		self.eraseradius = 2*boxsize
		
		self.dynapixp=get_app()
		self.image=EMData().read_images(imagefsp,[0])[0]					# original image to be boxed
		self.imagefsp = imagefsp
		self.currentimage = 0
		self.boxsetcache = None
		self.itshrink = -1 # image thumb shrink
		self.imagethumbs = None # image thumbs
		
		self.boxable = Boxable(self.image,self,None)
		self.boxable.addnonrefs(boxes)
		self.boxable.boxsize = boxsize
		
		self.autoBoxer = SwarmAutoBoxer(self.boxable)
		self.autoBoxer.boxsize = boxsize
		
		self.boxable.setAutoBoxer(self.autoBoxer)
		
		
		self.threshold=thr					# Threshold to decide which boxes to use
		self.ptcl=[]						# list of actual boxed out EMImages
		self.boxm = None
		self.moving=None					# Used during a user box drag
		
		self.guiimp=EMImage(self.image)		# widget for displaying large image
		self.guiim=self.guiimp.child
		self.guimxp= None # widget for displaying matrix of smaller images
		self.guimx=EMImageMX()	
		
		self.guimxitp = None
		self.genImageThumbnailsWidget(self.imagefsp)
			
		
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mousedown"),self.mousedown)
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mousedrag"),self.mousedrag)
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mouseup")  ,self.mouseup  )
		self.guiim.connect(self.guiim,QtCore.SIGNAL("keypress"),self.keypress)
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mousewheel"),self.mousewheel)
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mousemove"),self.mousemove)
		#self.guimx.connect(self.guimx,QtCore.SIGNAL("removeshape"),self.removeshape)
		self.guimx.connect(self.guimx,QtCore.SIGNAL("mousedown"),self.boxsel)
		self.guimx.connect(self.guimx,QtCore.SIGNAL("mousedrag"),self.boxmove)
		self.guimx.connect(self.guimx,QtCore.SIGNAL("mouseup"),self.boxrelease)
		
		
		self.mmode = 0
		self.guiim.setmmode(0)
		self.guimx.setmmode("app")
		self.ppc = 1.0
		self.mouseclicks = 0
		self.movingeraser = False
		self.guictl=GUIboxPanel(self)
		
		self.dynapix = True
		self.anchortemplate = False
		
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
		if self.guimxitp != None:
			self.guimxitp.show()
			self.guimxit.connect(self.guimxit,QtCore.SIGNAL("mousedown"),self.imagesel)
			
		self.guictl.show()
		
		self.boxupdate()
		
		#debug
		
		#f = file('autoparams.eman2','r')
		#from pickle import load
		#self.autoBoxer = load(f)
			
	def setPtclMxData(self,data=None):
		'''
		Call this to set the Ptcl Mx data 
		'''
		if data != None and len(data) != 0:
			self.guimx.setData(data)
			if self.guimxp == None:
				self.guimxp = EMParentWin(self.guimx)
				self.guimxp.show()
	
	def imagesel(self,event,lc):
		im=lc[0]
		if im != self.currentimage:
			oldboxset  = self.boxsetcache[self.currentimage]
			self.guimxit.setSelected(im)
			self.image = EMData().read_images(self.imagefsp,[im])[0]
			self.image.process_inplace("normalize.edgemean")
			self.guiim.setData(self.image)
			if self.boxsetcache[im] == None:
				self.boxsetcache[im] = Boxable(self.image,self,self.autoBoxer)
				self.boxsetcache[im].boxsize = self.boxsize
				
			self.boxable = self.boxsetcache[im]
			self.autoBoxer.setBoxable(self.boxable)
			if self.dynapix: self.boxable.autobox(self.autoBoxer)

			for box in self.boxable.boxes: box.changed = True
			self.currentimage = im
			self.ptcl = []
			self.guiim.delShapes()
			self.boxupdate()
			self.updateAllImageDisplay()
	
	def genImageThumbnailsWidget(self,imagename):
		'''
		Generates image thumbnails for a single image name
		if there is only one image in the image file on disk
		this function returns 0 and does nothing.
		Else thumbnails are generated, and image matrix is initialized (but not shown),
		and 1 is returned
		'''
		# warning, the imagename is invalid something could go wrong here
		try: n = EMUtil.get_image_count(imagename)
		except: 
			# warning - bad hacking going on
			print "the image name ", imagename, "probably doesn't exist"
			raise Exception
		
		nim = EMUtil.get_image_count(imagename)
		if (nim == 1): return 0
		
		n = self.getImageThumbShrink()
		self.imagethumbs = []
		
		for i in range(0,nim):
			if i == 0 and self.boxsetcache == None:
				self.boxsetcache = []
				self.boxsetcache.append(self.boxable)
				image = self.image
			else:
				self.boxsetcache.append(None)
				image = EMData().read_images(imagename,[i])[0]
			# could save some time here by using self.image when i == 0
			
			i = image.process("math.meanshrink",{"n":n})
			i.process_inplace("normalize.edgemean")
			self.imagethumbs.append(i)
		
		
		self.guimxit=EMImageMX()		# widget for displaying image thumbs
		self.guimxit.setData(self.imagethumbs)
		self.guimxit.setmmode("app")
		self.guimxitp = EMParentWin(self.guimxit)
		return 1
		
	def getImageThumbShrink(self):
		if self.itshrink == -1:
			if self.image == None:
				print "error - the image is not set, I need it to calculate the image thumb shrink"
				exit(1)
			shrink = 1
			inx = self.image.get_xsize()/2
			iny = self.image.get_ysize()/2
			while ( inx >= 128 and iny >= 128):
				inx /= 2
				iny /= 2
				shrink *= 2
		
			self.itshrink=shrink
		
		return self.itshrink
		
	def writeBoxesTo(self,imagename,norm=True):
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
		return self.boxable.boxes
	
	def mousemove(self,event):
		if self.mmode == 1 and event.modifiers()&Qt.ShiftModifier:
			m=self.guiim.scr2img((event.x(),event.y()))
			self.guiim.addShape("eraser",EMShape(["circle",.1,.1,.1,m[0],m[1],self.eraseradius,3]))
			self.updateImageDisplay()
			self.movingeraser = True
		else:
			if (self.movingeraser):
				self.guiim.addShape("eraser",EMShape(["circle",0,0,0,0,0,0,0.1]))
				self.updateImageDisplay()
				self.movingeraser = False
	def mousewheel(self,event):
		if self.mmode == 1 and event.modifiers()&Qt.ShiftModifier:
			self.guictl.adjustEraseRad(event.delta())
			m=self.guiim.scr2img((event.x(),event.y()))
			self.guiim.addShape("eraser",EMShape(["circle",.1,.1,.1,m[0],m[1],self.eraseradius,3]))
			self.updateImageDisplay()
	
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
				
				# If we get here, we need to add a new reference 
				box = Box(m[0]-self.boxsize/2,m[1]-self.boxsize/2,self.boxsize,self.boxsize,True)
				
				if self.anchortemplate:	box.isanchor = False
				else: box.isanchor = True # this is the default behaviour 
				box.changed = True # this is so image2D nows to repaint the shape
				self.boxable.addbox(box)
				
				
				# autoBoxer will autobox depending on the state of its mode
				self.autoBoxer.addReference(box)
				
				boxnum = len(self.getboxes())-1
				x0=boxes[boxnum].xcorner+boxes[boxnum].xsize/2-1
				y0=boxes[boxnum].ycorner+boxes[boxnum].ysize/2-1
				self.guiim.addShape("cen",EMShape(["rect",.9,.9,.4,x0,y0,x0+2,y0+2,1.0]))
				if not self.guimx.isVisible(boxnum) : self.guimx.scrollTo(boxnum,yonly=1)
				self.guimx.setSelected(boxnum)
				
				self.boxupdate()
				self.updateAllImageDisplay()
				self.mouseclicks += 1
				self.updateppc()
				return
				
				
			self.mouseclicks += 1
			# Deleting a box happens here
			if event.modifiers()&Qt.ShiftModifier :
				# with shift, we delete
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
			#self.boxable.addExclusionArea("circle",m[0],m[1],self.eraseradius)
			self.guiim.addShape("eraser",EMShape(["circle",.9,.9,.9,m[0],m[1],self.eraseradius,3]))
			self.updateAllImageDisplay()
			
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
			self.boxable.addExclusionArea("circle",m[0],m[1],self.eraseradius)
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
		for m in range(self.boxable.numboxes(),0,-1):
			box = boxes[m-1]
			if box.isref == False:
				self.delbox(m-1)

			#print "boxset now has",len(boxes)
	
	def mouseup(self,event) :
		if self.mmode == 0:
			m=self.guiim.scr2img((event.x(),event.y()))
			if self.moving != None:
				box = self.moving[0]
				if box.isref and self.dynapix:
					self.autoBoxer.referenceMoved(box)
					self.autobox()
			
			self.moving=None
		elif self.mmode == 1:
			self.guiim.addShape("eraser",EMShape(["circle",0,0,0,0,0,0,0.1]))
			
			lostboxes = self.boxable.updateExcludedBoxes()
			for n in lostboxes:
				self.delbox(n)
			self.mouseclicks += 1
			self.updateppc()
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
	
	#def removeshape(self,boxnum):
		#print "removing shape",boxnum
		#sh=self.guiim.getShapes()
		#k=sh.keys()
		#k.sort()
		#del sh[int(boxnum)]
		#for j in k:
			#if isinstance(j,int):
				#if j>boxnum :
					#sh[j-1]=sh[j]
					#del sh[j]
		#self.guiim.delShapes()
		#self.guiim.addShapes(sh)
		#self.guiim.setActive(None,.9,.9,.4)
		
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

		box = self.getboxes()[boxnum]
		self.ptcl.pop(boxnum)
		
		if self.boxable.delbox(boxnum) and self.dynapix:
			print 'told autoboxer to remove reference'
			self.autoBoxer.removeReference(box)

	def updateBoxSize(self,boxsize):
		if boxsize != self.boxsize:
			self.deletenonrefs()
			self.boxsize = boxsize
			self.boxable.updateBoxSize(boxsize)
			self.autoBoxer.setBoxSize(boxsize)
			
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
		self.setPtclMxData(self.ptcl)
		
		self.guictl.nboxesChanged(len(self.ptcl))
#		self.emit(QtCore.SIGNAL("nboxes"),len(self.ptcl))

	def run(self):
		"""If you make your own application outside of this object, you are free to use
		your own local app.exec_(). This is a convenience for boxer-only programs."""
		self.dynapixp.exec_()
		
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
	

	def autoboxbutton(self):
		print 'autobox button was pressed, this is currently disabled'

	def toggleDynapix(self,bool):
		if self.autoBoxer != None:
			if bool == True:
				print 'setting mode DYNAPIX'
				self.autoBoxer.setMode(SwarmAutoBoxer.DYNAPIX)
			else:
				print 'setting mode USERDRIVEN'
				self.autoBoxer.setMode(SwarmAutoBoxer.USERDRIVEN)
		
		self.dynapix = bool
		print 'self.dynapix is', self.dynapix
		
	def done(self):
		self.dynapixp.quit
		
	def trydata(self,data,thr):
		print 'trydata was pressed, this feature is currently disabled'
		
	def updatedata(self,data,thresh):
		self.guictl.updatedata(data,thresh)
		
	def setautobox(self,s):
		self.boxable.autoboxmethod = s
		
	def nocupdate(self,bool):
		self.anchortemplate = bool
		
		if bool:
			self.autoBoxer.setMode(SwarmAutoBoxer.ANCHOREDDYNAPIX)
		else:
			if self.dyanapix:
				self.autoBoxer.setMode(SwarmAutoBoxer.DYNAPIX)	
			else:
				self.autoBoxer.setMode(SwarmAutoBoxer.USERDRIVEN)	
	def classify(self,bool):
		self.boxable.genRefImages()
		
	def erasetoggled(self,bool):
		# for the time being there are only two mouse modes
		self.mmode = bool
		
	def updateEraseRad(self,rad):
		self.eraseradius = rad

	def quit(self):
		self.dynapixp.quit()
		
	def saveparams(self):
		name = QtGui.QFileDialog.getSaveFileName(None,'Save Autoboxing Params')
		print name
		if str(name) != '':
			f = file(str(name),'w')
		else: return
		
		from pickle import dump
		dump(self.autoBoxer,f,-1) # -1 forces the use of the HIGHEST_PROTOCOL, which is presumable for efficient
		f.close()
		
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
		self.dynapick.setChecked(True)
		self.hbl3.addWidget(self.dynapick)
		self.nocpick = QtGui.QRadioButton("Anchor")
		self.hbl3.addWidget(self.nocpick)
		
		self.vbl.addLayout(self.hbl3)
		
		self.autobox=QtGui.QPushButton("Auto Box")
		self.vbl.addWidget(self.autobox)


		self.hbl2=QtGui.QHBoxLayout()
		self.hbl2.setMargin(2)
		self.hbl2.setSpacing(6)
		#self.vbl.addLayout(self.hbl1)
		
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
		
		self.saveparams = QtGui.QPushButton("Save Params")
		self.vbl.addWidget(self.saveparams)

		self.classifybut=QtGui.QPushButton("Classify")
		self.vbl.addWidget(self.classifybut)
		
		self.connect(self.bs,QtCore.SIGNAL("editingFinished()"),self.newBoxSize)
		self.connect(self.eraserad,QtCore.SIGNAL("editingFinished()"),self.updateEraseRad)
		self.connect(self.thr,QtCore.SIGNAL("valueChanged"),self.newThresh)
		self.connect(self.done,QtCore.SIGNAL("clicked(bool)"),self.target.quit)
		self.connect(self.saveparams,QtCore.SIGNAL("clicked(bool)"),self.target.saveparams)
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
		self.target.guiim.setMouseTracking(bool)
		self.target.erasetoggled(bool)
	
	def dynapickd(self,bool):
		if bool == True:
			self.autobox.setEnabled(False)
		else:
			self.autobox.setEnabled(True)
		self.target.toggleDynapix(bool)
	
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
	
	def adjustEraseRad(self,delta):
		v = float(self.eraserad.text())
		if delta > 0:
			v = 1.1*v
		if delta < 0:
			v = 0.9*v
			
		self.eraserad.setText(str(int(v)))
		# this makes sure the target updates itself 
		# there may be a better approach, seeing as
		# the target called this function
		self.updateEraseRad()
		
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
		
		self.target.updateBoxSize(v)
		
	def newThresh(self,val):
		#print "new threshold"
		self.trythatd()


if __name__ == "__main__":
	main()
