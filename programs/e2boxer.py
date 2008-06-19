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
from pyemtbx.boxertools import *
from optparse import OptionParser
from emshape import EMShape
from emimagemx import EMImageMX
from math import *
from time import *
import os
import sys
import signal
from copy import *

from emglplot import *

from shelve import open

from time import time,sleep

from sys import getrefcount

if os.name == 'nt':
	def kill(pid):
		"""kill function for Win32"""
		import win32api
		handle = win32api.OpenProcess(1, 0, pid)
		return (0 != win32api.TerminateProcess(handle, 0))
else:
	from os import kill

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
	parser.add_option("--auto","-A",type="string",action="append",help="Autobox using specified method: ref, grid, db",default=[])
	parser.add_option("--writedb",action="store_true",help="Write data box files",default=False)
	parser.add_option("--writeboximages",action="store_true",help="Write data box files",default=False)
	parser.add_option("--force","-f",action="store_true",help="Force overwrites old files",default=False)
	parser.add_option("--overlap",type="int",help="(auto:grid) number of pixels of overlap between boxes. May be negative.")
	parser.add_option("--refptcl","-R",type="string",help="(auto:ref) A stack of reference images. Must have the same scale as the image being boxed.",default=None)
	parser.add_option("--nretest",type="int",help="(auto:ref) Number of reference images (starting with the first) to use in the final test for particle quality.",default=-1)
	parser.add_option("--retestlist",type="string",help="(auto:ref) Comma separated list of image numbers for retest cycle",default="")
	parser.add_option("--farfocus",type="string",help="filename or 'next', name of an aligned far from focus image for preliminary boxing",default=None)
	parser.add_option("--dbout",type="string",help="filename to write EMAN1 style box database file to",default=None)
	parser.add_option("--norm",action="store_true",help="Edgenormalize boxed particles",default=False)
	parser.add_option("--ptclout",type="string",help="filename to write boxed out particles to",default=None)
	parser.add_option("--parallel",type="int",help="specify more than one processor",default=1)
	parser.add_option("--mergedb",action="store_true",help="A special argument, if true all input arguments are considered to be box files and they are merged into the project database as manually selected particles",default=False)
	
	(options, args) = parser.parse_args()
	if len(args)<1 : parser.error("Input image required")
	
	filenames = []
	for i in range(0,len(args)):
		filenames.append(args[i])

	logid=E2init(sys.argv)
	
	
	if options.mergedb == True:
		mergedb(filenames)
		sys.exit(1)
	
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
				
	if "db" in options.auto:
		print "auto data base boxing"
	
		if len(filenames) == 1:
			projectdb = EMProjectDB()
	
			trimAutoBoxer = projectdb["currentautoboxer"]
			autoBoxer = SwarmAutoBoxer(None)
			autoBoxer.become(trimAutoBoxer)
			autoBoxer.setModeExplicit(SwarmAutoBoxer.COMMANDLINE)
			imagename = filenames[0]
			exists = True
			try:
				oldAutoBoxer = projectdb[imagename+"_autoboxer"]	
			except: exists = False 	
			if exists and autoBoxer.stateTS == oldAutoBoxer.stateTS:
				print "The content in the project data base for",imagename,"is up2date"
				exit(1)
			else:
				print "Auto boxing",imagename
				image = EMData(imagename)
				boxable = Boxable(image,imagename,None,autoBoxer)
				autoBoxer.setBoxable(boxable)
				autoBoxer.autoBox(boxable)
				if options.writedb:
					boxable.writedb(options.force)
				if options.writeboximages:
					boxalbe.writeboximages(options.force)
				exit(1)
		else:
			print "autoboxing using parallel stuff"
			autoboxer = AutoDBBoxer(filenames,options.parallel,options.force)
			try:
				from emimage import EMImage,get_app
			except: 
				print "error"
				exit(1)
				
			a = get_app()
			autoboxer.go(a)
			a.exec_()
			print "done"
				
			exit(1)
			#while autoboxer.working:
				#sleep(1)
				#autoboxer.printcpstatus()
				
				
			#exit(1)

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
		gui=GUIbox(filenames,boxes,boxthr,options.boxsize)
		gui.run()
		
	print "finished running"

	if options.dbout:
		boxes = gui.getBoxes()
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


def mergedb(filenames):
	projectdb = EMProjectDB()
	for filename in filenames:
		f=file(filename,'r')
		lines=f.readlines()
		boxes = []
		for line in lines:
			data = str.split(line)
			b = Box(int(data[0]),int(data[1]),int(data[2]),int(data[3]),False)
			b.ismanual = True
			boxes.append(TrimBox(b))
	
		try:
			manualboxes = projectdb[strip_file_tag(filename)+"_manualboxes"]
		except:
			manualboxes = []
	
		manualboxes.extend(boxes)
		projectdb[strip_file_tag(filename)+"_manualboxes"] = manualboxes
		
	projectdb.close()
		

class AutoDBBoxer(QtCore.QObject):
	def __init__(self,imagenames,nproc,force=False):
		QtCore.QObject.__init__(self)
		self.nproc = nproc
		self.imagenames = imagenames
		self.currentidx = 0	
		self.force = force
		self.working = True
		self.processes = []
		self.jobsdone = 0
		self.cps = []
		self.app = None
		for i in range(0,nproc):
			self.cps.append(None)
		
		for i in range(0,len(self.imagenames)):
			self.processes.append(QtCore.QProcess(self))
			
	def printcpstatus(self):
		for i in self.cps:
			print i.state(),
			print i.pid()
			print kill(i.pid(),signal.SIG_IGN)
			
		print ''
	def go(self,app):
		self.app = app
		
		for i in range(0,self.nproc):
			self.spawn_process()
			self.cps[i] = self.processes[i]
			self.currentidx += 1
			
		
			
	def spawn_process(self):
		if self.currentidx >= len(self.imagenames) :
			return
			
		#process = self.processes[self.currentidx]
			
		program = QtCore.QString("e2boxer.py")
		args = QtCore.QStringList()
		args.append(self.imagenames[self.currentidx])
		args.append("--auto=db")
		
		if self.force:	args.append("-f")
		
		QtCore.QObject.connect(self.processes[self.currentidx], QtCore.SIGNAL("finished(int)"), self.process_finished)
		#QtCore.QObject.connect(self.processes[self.currentidx], QtCore.SIGNAL("finished(int,QProcess.ExitStatus)"), self.process_finished_status)
		#QtCore.QObject.connect(self.processes[self.currentidx], QtCore.SIGNAL("error(QProcess.ProcessError)"), self.process_error)
		#QtCore.QObject.connect(self.processes[self.currentidx], QtCore.SIGNAL("started()"), self.process_start)
		
		#self.processes[self.currentidx].setProcessChannelMode(QtCore.QProcess.ForwardedChannels)
		self.processes[self.currentidx].start(program,args)
		
		
		#self.processes[self.currentidx].waitForStarted()
		print "executing",
		for arg in args: print arg,
		print ''
	
	def process_finished(self,int):
		#print "process finished"
		self.jobsdone += 1
		if self.jobsdone == len(self.imagenames):
			self.app.quit()
		self.spawn_process()
		self.currentidx += 1
		
		
	#def process_start(self):
		#print "process started"
		

class GUIboxMouseEventsObject:
	'''
	A base class for objects that handle mouse events in the GUIbox
	
	Inheriting objects are concerned with supplying their own definitions of 
	the functions mouseup, mousedown, mousedrag, mousemove and mousewheel. 
	They do not have to supply their own definition for all of these functions,
	but it would not make sense to inherit from this class unless the child class
	did not atleast supply one of them.
	
	Also stores a reference to the mediator object which coordinates the requests and signals
	of subclassing objects to the GUIbox class. Modify the mediator class if you ever need
	to extend the scope of the messaging/signaling/requesting system between GUIboxMouseEventsObjects
	and the GUIbox.
	
	The mediator object could be replaced with a direct reference to the GUIbox class, but
	this is not recommended - for documentation and clarity of code purposes.
	'''
	def __init__(self,mediator):
		'''
		Stores only a reference to the mediator
		'''
		if not isinstance(mediator,GUIboxEventsMediator):
			print "error, the mediator should be a GUIboxEventsMediator"
			return
		
		self.mediator = mediator
		
	def get2DGuiImage(self):
		'''
		Ask the mediator to for the EMImage2D object
		'''
		return self.mediator.get2DGuiImage()
	
	def getGuiCtl(self):
		'''
		Ask the mediator to for the main GUI controller
		'''
		return self.mediator.getGuiCtl()
	
	def getMXGuiImage(self):
		'''
		Ask the mediator to for the EMImageMX object
		'''
		return self.mediator.getMXGuiImage()
	
	def mouseup(self,event):
		'''
		Inheriting classes to potentially define this function
		'''
		pass

	def mousedown(self,event):
		'''
		Inheriting classes to potentially define this function
		'''
		pass
		
	def mousedrag(self,event):
		'''
		Inheriting classes to potentially define this function
		'''
		pass
	
	def mousemove(self,event):
		'''
		Inheriting classes to potentially define this function
		'''
		pass

	def mousewheel(self,event):
		'''
		Inheriting classes to potentially define this function
		'''
		pass
	

class GUIboxMouseEraseEvents(GUIboxMouseEventsObject):
	'''
	A class that knows how to handle mouse erase events for a GUIBox
	'''
	def __init__(self,mediator,eraseradius=-1):
		GUIboxMouseEventsObject.__init__(self,mediator)
		
		self.eraseradius=eraseradius	# This is a circular radius 
		self.erasemode = None			# erase mode can be either Boxable.ERASE or Boxable.UNERASE
		
	def setmode(self,mode):
		self.erasemode = mode
		
	def setEraseRadius(self,radius):
		self.eraseradius = radius
	
	def mousemove(self,event):
		m = self.get2DGuiImage().scr2img((event.x(),event.y()))
		self.get2DGuiImage().addShape("eraser",EMShape(["circle",.1,.1,.1,m[0],m[1],self.eraseradius,3]))
		self.mediator.updateImageDisplay()
		
	def mousewheel(self,event):
		if event.modifiers()&Qt.ShiftModifier:
			self.getGuiCtl().adjustEraseRad(event.delta())
			m= self.get2DGuiImage().scr2img((event.x(),event.y()))
			self.get2DGuiImage().addShape("eraser",EMShape(["circle",.1,.1,.1,m[0],m[1],self.eraseradius,3]))
			self.mediator.updateImageDisplay()
	
	def mousedown(self,event) :
		m=self.get2DGuiImage().scr2img((event.x(),event.y()))
		#self.boxable.addExclusionArea("circle",m[0],m[1],self.eraseradius)
		self.get2DGuiImage().addShape("eraser",EMShape(["circle",.9,.9,.9,m[0],m[1],self.eraseradius,3]))
		self.mediator.exclusionAreaAdded("circle",m[0],m[1],self.eraseradius,self.erasemode)	

	def mousedrag(self,event) :
		m=self.get2DGuiImage().scr2img((event.x(),event.y()))
		self.get2DGuiImage().addShape("eraser",EMShape(["circle",.9,.9,.9,m[0],m[1],self.eraseradius,3]))
		self.mediator.exclusionAreaAdded("circle",m[0],m[1],self.eraseradius,self.erasemode)
		# exclusionAreaAdded does the OpenGL update calls, so there is no need to do so here
		
	def mouseup(self,event) :
		# we have finished erasing
		
		# make the eraser shape non visible
		self.get2DGuiImage().addShape("eraser",EMShape(["circle",0,0,0,0,0,0,0.1]))
		self.mediator.erasingDone()
	
class GUIboxParticleManipEvents(GUIboxMouseEventsObject):
	'''
	A class that knows how to add and remove reference and non reference boxes 
	'''
	def __init__(self,mediator):
		GUIboxMouseEventsObject.__init__(self,mediator)
		self.mode =  GUIbox.REFERENCE_ADDING
		self.anchoring = False
		self.moving = None
		self.dynapix = False
		
	def setAnchoring(self,bool):
		self.anchoring = bool
		
	def setmode(self,mode):
		if mode not in [GUIbox.REFERENCE_ADDING,GUIbox.MANUALLY_ADDING]:
			print 'error, that is  an illegal mode'
			return
		
		self.mode = mode
		
	def mousedown(self,event) :
		m= self.get2DGuiImage().scr2img((event.x(),event.y()))
		boxnum = self.mediator.detectBoxCollision(m)
		if boxnum == -1:
			#if we make it here, that means the user has clicked on an area that is not in any box
			
			if event.modifiers()&Qt.ShiftModifier : return # the user tried to delete nothing
			
			# If we get here, we need to add a new reference
			boxsize = self.mediator.getBoxSize()
			
			box = Box(m[0]-boxsize/2,m[1]-boxsize/2,boxsize,boxsize,True)
			box.setImageName(self.mediator.getCurrentImageName())
			
			if self.anchoring:	box.isanchor = False
			else: box.isanchor = True # this is the default behaviour 
			box.changed = True # this is so image2D nows to repaint the shape
			
			if self.mode == GUIbox.REFERENCE_ADDING:
				box.isref = True
				box.ismanual = False
				
			elif self.mode == GUIbox.MANUALLY_ADDING:
				box.isref = False
				box.ismanual = True
			else:
				print 'error, unknown error in mousedown, boxing mode'
			
			x0=box.xcorner+box.xsize/2-1
			y0=box.ycorner+box.ysize/2-1
			self.get2DGuiImage().addShape("cen",EMShape(["rect",.9,.9,.4,x0,y0,x0+2,y0+2,1.0]))
			
			self.mediator.storeBox(box)
			self.mediator.mouseClickUpdatePPC()
		
		elif event.modifiers()&Qt.ShiftModifier :
			# remove the box
			self.mediator.removeBox(boxnum)
			self.mediator.mouseClickUpdatePPC()
			
		else:
			# if we make it here than the we're moving a box
			box = self.mediator.getBox(boxnum)
			self.moving=[box,m,boxnum]
			self.get2DGuiImage().setActive(boxnum,.9,.9,.4)
				
			x0=box.xcorner+box.xsize/2-1
			y0=box.ycorner+box.ysize/2-1
			self.get2DGuiImage().addShape("cen",EMShape(["rect",.9,.9,.4,x0,y0,x0+2,y0+2,1.0]))
			if not self.getMXGuiImage().isVisible(boxnum) : self.getMXGuiImage().scrollTo(boxnum,yonly=1)
			self.getMXGuiImage().setSelected(boxnum)
			self.mediator.updateAllImageDisplay()
		
	def mousedrag(self,event) :
		
		m=self.get2DGuiImage().scr2img((event.x(),event.y()))
		
		if event.modifiers()&Qt.ShiftModifier:
			boxnum = self.mediator.detectBoxCollision(m)
			if ( boxnum != -1):
				self.mediator.removeBox(boxnum)
				self.mediator.mouseClickUpdatePPC()
			
		elif self.moving != None:
			# self.moving[0] is the box, self.moving[1] are the mouse coordinates
			box = self.moving[0]
			# the old m in in self.moving[2]
			oldm = self.moving[1]
			
			self.mediator.moveBox(self.moving[2],m[0]-oldm[0],m[1]-oldm[1])
			self.moving[1] = m
	
	def mouseup(self,event) :
		if self.moving != None:
			box = self.moving[0]
			if box.isref:
				print "calling reference moved"
				self.mediator.referenceMoved(box)
			
		self.moving=None

class GUIboxEventsMediator:
	'''
	This class could just as easily not exist - however it remains in use for documentation
	purposes. If it was removed, the GUIboxMouseEventsObjects could have a reference to the GUIbox
	class instead of this one and the code would still work. But without this class it would be difficult to
	disentangle the relationship between the GUIboxMouseEventsObjects and the GUIbox.
	
	This class coordinates all the requests and 'signals' of the GUIboxMouseEventsObject classes so that they
	are connected to the correct 'slots' and getter functions  in the GUIbox class. This behavior is analogous to the 
	Qt signal/slot mechanism, but no Qt signals/slots are actually used. Instead this class
	supplies a bunch of public functions that accept the requests and signals from the GUIboxMouseEventsObject
	and send them on to the 'slots' and getter functions of the the GUIbox - using basic function interfacing. This
	is also motivated by the Mediator concept in the Gang of Four.
	
	All things considered, the class remains for documentation purposes. It should only be removed if 
	it poses a significant performance hit
	'''
	def __init__(self,parent):
		'''
		Stores only a reference to the parent
		'''
		if not isinstance(parent,GUIbox):
			print "error, the parent of a GUIboxMouseEraseEvents must be a GUIbox"
			return
		
		self.parent = parent	# need a referene to the parent to send it events

	def get2DGuiImage(self):
		'''
		Return the parent's EMImage2D object
		'''
		return self.parent.get2DGuiImage()
	
	def getGuiCtl(self):
		'''
		Return the parent's Controller widgit object
		'''
		return self.parent.getGuiCtl()
	
	def getMXGuiImage(self):
		'''
		Return the parent's EMImageMX object
		'''
		return self.parent.getMXGuiImage()
	
	def updateImageDisplay(self):
		'''
		Send an event to the parent that the EMImage2D should update its display
		'''
		self.parent.updateImageDisplay()
		
	def updateAllImageDisplay(self):
		'''
		Send an event to the parent that the EMImage2D and EMImageMX objects should
		update their displays
		'''
		self.parent.updateAllImageDisplay()
		
	def exclusionAreaAdded(self,typeofexclusion,x,y,radius,mode):
		'''
		Send an event to the parent that an exclusion area was added.
		The parameters define the type of exclusion area. In future the exclusion area
		should probably just be its own class.
		'''
		self.parent.exclusionAreaAdded(typeofexclusion,x,y,radius,mode)

	def erasingDone(self):
		'''
		Send an event to the parent letting it know that the user has stopped adding
		erased area
		'''
		self.parent.erasingDone()
		
	def detectBoxCollision(self,coords):
		'''
		Ask the parent to detect a collision between a resident box and the given coordinates.
		This is in terms of the EMImage2D 
		'''
		return self.parent.detectBoxCollision(coords)
	
	def getCurrentImageName(self):
		'''
		Ask the parent for the current name of the image in the large image view...
		'''
		return self.parent.getCurrentImageName()

	def getBoxSize(self):
		'''
		Ask the parent for the current project boxsize
		'''
		return self.parent.getBoxSize()
	
	def storeBox(self,box):
		'''
		Tell the parent to store a box
		'''
		self.parent.storeBox(box)
	
	def boxDisplayUpdate(self):
		'''
		Tell the parent the a general box display update needs to occur
		'''
		self.parent.boxDisplayUpdate()
		
	def mouseClickUpdatePPC(self):
		'''
		Tell the parent that a mouse click occured and that the PPC metric should be updated
		'''
		self.parent.mouseClickUpdatePPC()
	
	def removeBox(self,boxnum):
		'''
		Tell the parent to remove a box, as given by the boxnum
		'''
		self.parent.removeBox(boxnum)
		
	def getBox(self,boxnum):
		'''
		Ask the parent for the box, as given by the boxnum
		'''
		return self.parent.getBox(boxnum)
	
	def moveBox(self,boxnum,dx,dy):
		'''
		Tell the parent to handle the movement of a box
		'''
		self.parent.moveBox(boxnum,dx,dy)
	
	def referenceMoved(self,box):
		'''
		Tell the parent that a reference was moved - this could trigger automatic boxing
		'''
		self.parent.referenceMoved(box)
	
class GUIbox:
	REFERENCE_ADDING = 0
	ERASING = 1
	MANUALLY_ADDING = 2
	def __init__(self,imagenames,boxes,thr,boxsize=-1):
		"""Implements the 'boxer' GUI."""
		
		try:
			from emimage import EMImage,get_app
		except:
			print "Cannot import EMAN image GUI objects (emimage,etc.)"
			sys.exit(1)

		try:
			os.mkdir(EMProjectDB.outputdir)
		except: pass
		self.dynapix = False
		self.anchoring = False
		
		if len(boxes)>0 and boxsize==-1: self.boxsize=boxes[0][2]
		elif boxsize==-1: self.boxsize=128
		else: self.boxsize=boxsize
		
		self.eraseradius = 2*boxsize
		self.erasemode = None
		self.dynapixp=get_app()
		self.imagenames = imagenames
		self.currentimage = 0
		self.itshrink = -1 # image thumb shrink
		self.imagethumbs = None # image thumbs
		
		try:
			projectdb = EMProjectDB()
			data = projectdb[self.imagenames[0]+"_autoboxer"]
			trimAutoBoxer = projectdb[data["auto_boxer_db_name"]]
			self.autoboxername = data["auto_boxer_unique_id"] 
			self.autoBoxer = SwarmAutoBoxer(self)
			self.autoBoxer.become(trimAutoBoxer)
			self.autoBoxer.setMode(self.dynapix,self.anchoring)
			print 'using cached autoboxer db'
		except:
			self.autoBoxer = SwarmAutoBoxer(self)
			self.autoBoxer.boxsize = boxsize
			self.autoBoxer.setMode(self.dynapix,self.anchoring)
			print self.autoBoxer.boxsize,"in constructor"
			print "loaded a new autoboxer"
		
		self.boxable = Boxable(self.imagenames[0],self,self.autoBoxer)
		self.boxable.addnonrefs(boxes)
		self.boxable.boxsize = boxsize
		self.boxable.setAutoBoxer(self.autoBoxer)
		
		self.eventsmediator = GUIboxEventsMediator(self)
		self.mousehandlers = {}
		self.mousehandlers["boxing"] = GUIboxParticleManipEvents(self.eventsmediator)
		self.mousehandlers["erasing"] = GUIboxMouseEraseEvents(self.eventsmediator,self.eraseradius)
		self.mousehandler = self.mousehandlers["boxing"]
		
		self.ptcl=[]						# list of actual boxed out EMImages
		self.boxm = None
		self.moving=None					# Used during a user box drag
		
		bic = BigImageCache()
		image=bic.getImage(self.imagenames[0])
		self.guiimp=EMImage(image)		# widget for displaying large image
		self.guiim=self.guiimp.child
		self.guiim.setOtherData(self.boxable.getExclusionImage(False),self.autoBoxer.getBestShrink(),True)
		self.guiim.setFrozen(self.boxable.isFrozen())
		
		self.guimxp= None # widget for displaying matrix of smaller imagespaugay
		self.guimx=EMImageMX()	
		
		self.guimxitp = None
		self.genImageThumbnailsWidget()

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
		self.guimx.connect(self.guimx,QtCore.SIGNAL("boxdeleted"),self.boximagedeleted)
		
		self.indisplaylimbo = False	# a flag I am using to solve a problem
		self.mmode = GUIbox.REFERENCE_ADDING
		self.guiim.setmmode(0)
		self.guimx.setmmode("app")
		self.ppc = 1.0
		self.mouseclicks = 0
		
		self.abselmediator = AutoBoxerSelectionsMediator(self)
		self.guictl=GUIboxPanel(self,self.abselmediator)
		
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
			self.guimxit.setSelected(0)

		self.guictl.show()
		self.boxDisplayUpdate()
	
	def changeCurrentAutoBoxer(self, timestamp):
		#print "change current autoboxer"
		projectdb = EMProjectDB()
		trimAutoBoxer = projectdb["autoboxer_"+timestamp]
		#print "changing autoboxer to autoboxer_",timestamp,"and its stamp in the db is",trimAutoBoxer.getCreationTS()
		self.autoBoxer = SwarmAutoBoxer(self)
		self.autoBoxer.become(trimAutoBoxer)
		self.autoBoxer.writeSpecificReferencesToDB(self.boxable.getImageName())
		
		self.ptcl = []
		self.guiim.delShapes()
		
		self.boxable.setAutoBoxer(self.autoBoxer)
		self.boxable.clearAndReloadImages()
		self.indisplaylimbo = True
		self.autoBoxer.regressiveflag = True
		self.autoBoxer.autoBox(self.boxable)
		self.indisplaylimbo = False

		self.boxDisplayUpdate()
		
	def get2DGuiImage(self):
		return self.guiim
	
	def getGuiCtl(self):
		return self.guictl
	
	def getMXGuiImage(self):
		return self.guimx
	
	def getBoxable(self):
		return self.boxable
		
	def exclusionAreaAdded(self,typeofexclusion,x,y,radius,mode):
		self.boxable.addExclusionArea(typeofexclusion,x,y,radius,mode)
		self.guiim.setOtherData(self.boxable.getExclusionImage(False),self.autoBoxer.getBestShrink(),True)
		self.updateImageDisplay()
		
	def erasingDone(self):
		'''
		Call this function after erasing has occured to remove all boxes in the
		erased regions
		'''
		
		# tell the boxable to remove boxes (and return their number)
		[lostboxes,refboxes] = self.boxable.updateExcludedBoxes()
		# after this, make sure the display is correct.
		
		# If any of the boxes are references the autoBoxer needs to know about it...
		# this could potentially trigger auto boxer 
		if len(lostboxes) != 0:
			self.deleteDisplayShapes(lostboxes)
			self.mouseclicks += 1
			self.updateppc()
			self.updateAllImageDisplay()
			
		else: self.updateImageDisplay()
		
		if len(refboxes) != 0: self.autoBoxer.removeReference(refboxes)
		self.boxDisplayUpdate()
	
	def detectBoxCollision(self,coords):
		'''
		Detects a collision of the coordinates with any of the boxes in the current image
		stored in guiim. coords need to be in the image coordinates
		'''
		return  self.collisiondetect(coords,self.getBoxes())
	
	def getCurrentAutoBoxerTS(self):
		return self.autoBoxer.getCreationTS()
	
	def getCurrentImageName(self):
		return self.imagenames[self.currentimage]
	
	def getImageNames(self):
		return self.imagenames
	
	def getBoxSize(self):
		return self.boxsize
	
	def storeBox(self,box):
		self.boxable.addbox(box)
		
		# Should this be here?
		boxnum = len(self.getBoxes())
		if not self.guimx.isVisible(boxnum) : self.guimx.scrollTo(boxnum,yonly=1)
		self.guimx.setSelected(boxnum)
		
		# autoBoxer will autobox depending on the state of its mode
		if box.isref : self.autoBoxer.addReference(box)
		
		self.boxDisplayUpdate()
		
	def getBox(self,boxnum):
		boxes = self.getBoxes()
		return boxes[boxnum]
		
	def removeBox(self,boxnum):
		box = self.delbox(boxnum)
		if (not box.isref or box.ismanual):
			self.boxable.addExclusionParticle(box)
		self.guiim.setOtherData(self.boxable.getExclusionImage(False),self.autoBoxer.getBestShrink(),True)
		self.boxDisplayUpdate()
		
	def referenceMoved(self,box):
		if (self.autoBoxer.referenceMoved(box)):
			self.boxDisplayUpdate()
		
	def mouseClickUpdatePPC(self):
		self.mouseclicks += 1
		self.updateppc()
	
	def setPtclMxData(self,data=None):
		'''
		Call this to set the Ptcl Mx data 
		'''
		if data != None and len(data) != 0:
			self.guimx.setData(data)
			if self.guimxp == None:
				self.guimxp = EMParentWin(self.guimx)
				self.guimxp.show()
				self.guimx.setSelected(0)
	
	def imagesel(self,event,lc):
		#print 'in image select'
		im=lc[0]
		if im != self.currentimage:
			#print 'changing images'
			self.guimxit.setSelected(im)
			
			bic = BigImageCache()
			image=bic.getImage(self.imagenames[im])
			self.guiim.setData(image)
			
			self.boxable.cacheExcToDisk()
			self.boxable = Boxable(self.imagenames[im],self,self.autoBoxer)
			
			self.ptcl = []
			self.guiim.delShapes()
			self.indisplaylimbo = True
			
			try:
				projectdb = EMProjectDB()
				data = projectdb[self.imagenames[im]+"_autoboxer"]
				trimAutoBoxer = projectdb[data["auto_boxer_db_name"]]
				self.autoboxername = data["auto_boxer_unique_id"] 
				self.autoBoxer = SwarmAutoBoxer(self)
				self.autoBoxer.become(trimAutoBoxer)
			except: pass
			
			self.autoBoxer.regressiveflag = True
			self.autoBoxer.autoBox(self.boxable)
			self.autoBoxerDBChanged()
			
			self.indisplaylimbo = False
			
			for box in self.boxable.boxes: box.changed = True
			
			self.currentimage = im
			
			self.guiim.setOtherData(self.boxable.getExclusionImage(False),self.autoBoxer.getBestShrink(),True)
			self.guiim.setFrozen(self.boxable.isFrozen())
			self.boxDisplayUpdate()
			

	def genImageThumbnailsWidget(self):
		'''
		Generates image thumbnails for a single image name
		if there is only one image in the image file on disk
		this function returns 0 and does nothing.
		Else thumbnails are generated, and image matrix is initialized (but not shown),
		and 1 is returned
		'''
		# warnilg here
		try: nim = len(self.imagenames)
		except: 
			# warning - bad hacking going on
			print "the image name ", imagename, "probably doesn't exist"
			raise Exception
		
		if (nim == 1): return 0
		
		n = self.getImageThumbShrink()
		self.imagethumbs = []
		
		for i in range(0,nim):
			thumb = self.getImageThumb(i)
			#print "got thumb",i
			self.imagethumbs.append(thumb)
		
		
		self.guimxit=EMImageMX()		# widget for displaying image thumbs
		self.guimxit.setData(self.imagethumbs)
		self.guimxit.setmmode("app")
		self.guimxitp = EMParentWin(self.guimxit)
		
		return 1
		
	def getImageThumb(self,i):
		n = self.getImageThumbShrink()
		
		bic = BigImageCache()
		image=bic.getImage(self.imagenames[i])
		
		thumb = image.process("math.meanshrink",{"n":n})
		thumb.process_inplace("normalize.edgemean")
		return thumb
		
	def getImageThumbShrink(self):
		if self.itshrink == -1:
			bic = BigImageCache()
			image=bic.getImage(self.imagenames[self.currentimage])
			if image == None:
				print "error - the image is not set, I need it to calculate the image thumb shrink"
				exit(1)
			shrink = 1
			inx =  image.get_xsize()/2
			iny =  image.get_ysize()/2
			while ( inx >= 128 and iny >= 128):
				inx /= 2
				iny /= 2
				shrink *= 2
		
			self.itshrink=shrink
		
		return self.itshrink
		
	def writeBoxesTo(self,imagename,norm=True):
		boxes = self.getBoxes()
		
		# Write EMAN1 style box database
		n = 0
		for box in boxes:
			image = box.getBoxImage()
#			print n,i
#			print i[4]
			image.write_image(imagename,n)
			n += 1

	def boxmove(self,event,scale):
		dx = (self.boxm[0] - event.x())/scale
		dy = (event.y() - self.boxm[1])/scale
		self.moveBox(self.boxm[2],dx,dy)
		self.boxm[0] = event.x()
		self.boxm[1] = event.y()
		self.updateAllImageDisplay()
		
	def boxrelease(self,event):
		boxes = self.getBoxes()
		#print boxes
		box =  boxes[self.boxm[2]]
		if box.isref :
			self.autoBoxer.referenceMoved(box)
		self.boxm = None
		
	def boxsel(self,event,lc):
		im=lc[0]
		self.boxm = [event.x(),event.y(),im]
		self.guiim.setActive(im,.9,.9,.4)
		self.guimx.setSelected(im)
		boxes = self.getBoxes()
		self.guiim.registerScrollMotion(boxes[im].xcorner+boxes[im].xsize/2,boxes[im].ycorner+boxes[im].ysize/2)
		try:
			#self.guiim.scrollTo(boxes[im].xcorner+boxes[im].xsize/2,boxes[im].ycorner+boxes[im].ysize/2)
			pass
			
		except: print "boxsel() scrolling error"

	def getBoxes(self):
		return self.boxable.boxes
	
	def mousemove(self,event):
		self.mousehandler.mousemove(event)

	def mousewheel(self,event):
		self.mousehandler.mousewheel(event)
		
	def mousedown(self,event) :
		self.mousehandler.mousedown(event)
		
	def mouseup(self,event) :
		self.mousehandler.mouseup(event)
	
	def mousedrag(self,event) :
		self.mousehandler.mousedrag(event)
	
	def updateppc(self):
		if self.mouseclicks > 0:
			self.guictl.ppcChanged(len(self.getBoxes())/float(self.mouseclicks))
		else:
			self.guictl.ppcChanged(0)
	
	def collisiondetect(self,m,boxes):
			
		for boxnum,box in enumerate(boxes):
			if m[0]<box.xcorner or m[0]>(box.xcorner +box.xsize) or m[1]<box.ycorner or m[1]>(box.ycorner +box.ysize) :
				# no collision
				continue
			# if we make it here there has been a collision, the box already exists
			return boxnum
		
		return -1

	def moveBox(self,boxnum,dx,dy):
		box = self.getBoxes()[boxnum]
		self.boxable.moveBox(box,dx,dy,boxnum)
			# we have to update the reference also
		self.ptcl[boxnum] = box.getBoxImage()
			
		x0=box.xcorner+box.xsize/2-1
		y0=box.ycorner+box.ysize/2-1
		self.guiim.addShape("cen",EMShape(["rect",.9,.9,.4,x0,y0,x0+2,y0+2,1.0]))
		box.shape = EMShape(["rect",box.r,box.g,box.b,box.xcorner,box.ycorner,box.xcorner+box.xsize,box.ycorner+box.ysize,2.0])
		self.guiim.addShape(boxnum,box.shape)
		self.boxDisplayUpdate()
		#self.updateAllImageDisplay()
		
	def updateboxcolors(self,classify):
		sh=self.guiim.getShapes()
		for i in classify.items():
			color = BoxingTools.get_color(i[1])
			
			sh[int(i[0])].shape[1] = color[0]
			sh[int(i[0])].shape[2] = color[1]
			sh[int(i[0])].shape[3] = color[2]
			sh[int(i[0])].changed=True
			
		self.boxDisplayUpdate()

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
	
	def deleteDisplayShapes(self,numbers):
		'''
		Warning - this won't work unless the numbers go from greatest to smallest- i.e. they are in reverse order
		Deletes shapes displayed by the 2D image viewer
		Pops boxed particles from the list used by the matrix image viewer (for boxes)
		'''
		#print "called delete display shapesS"
		if self.indisplaylimbo: return
		
		sh=self.guiim.getShapes()
		
		for num in numbers:
			sh=self.guiim.getShapes()
			k=sh.keys()
			k.sort()
			del sh[int(num)]
			for j in k:
				if isinstance(j,int):
					if j>num :
						sh[j-1]=sh[j]
						del sh[j]
			self.ptcl.pop(num)
						
		self.guiim.delShapes()
		self.guiim.addShapes(sh)

			
		#self.guiim.delShapes()
		#self.guiim.addShapes(sh)
		#print "now there are",len(sh),"shapes"
		
	def boximagedeleted(self,event,lc,forceimagemxremove=True):
		box = self.delbox(lc[0],forceimagemxremove)
		self.boxable.addExclusionParticle(box)
		self.guiim.setOtherData(self.boxable.getExclusionImage(False),self.autoBoxer.getBestShrink(),True)
		self.updateAllImageDisplay()
		self.updateppc()
		
	def delbox(self,boxnum,forceimagemxremove=True):
		"""
		Deletes the numbered box completely
		Should only be called in the instance where a single box is being deleting - NOT when
		you are deleting a list of boxes sequentially (for that you should use deleteDisplayShapes
		and something to pop the box from the Boxable. See examples in this code)
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
		if forceimagemxremove: self.ptcl.pop(boxnum)
		
		box = self.boxable.boxes[boxnum]
		self.boxable.delbox(boxnum)
		# if the boxable was a reference then the autoboxer needs to be told. It will remove
		# it from its own list and potentially do autoboxing
		if box.isref:
			self.autoBoxer.removeReference(box)
		
		return box
	
	def updateBoxSize(self,boxsize):
		if boxsize != self.boxsize:
			self.boxsize = boxsize
			self.boxable.updateBoxSize(boxsize)
			self.autoBoxer.setBoxSize(boxsize)
			
	def boxDisplayUpdate(self,force=False):
		
		ns = {}
		idx = 0
		#self.ptcl = []
		# get the boxes
		boxes =self.getBoxes()
		for j,box in enumerate(boxes):
	
			if not box.changed and not force:
				idx += 1
				continue
			
			box.changed=False
		
			im=box.getBoxImage()
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

		self.updateAllImageDisplay()

	def run(self):
		"""If you make your own application outside of this object, you are free to use
		your own local app.exec_(). This is a convenience for boxer-only programs."""
		self.dynapixp.exec_()
		
		self.boxable.cacheExcToDisk()
		projectdb = EMProjectDB()
		projectdb.close()
		

		E2saveappwin("boxer","imagegeom",self.guiim)
		E2saveappwin("boxer","matrixgeom",self.guimx)
		E2saveappwin("boxer","controlgeom",self.guictl)
		#E2setappval("boxer","matrixnperrow",self.guimx.nperrow)
		try:
			E2setappval("boxer","imcontrol",self.guiim.inspector.isVisible())
			if self.guiim.inspector.isVisible() : E2saveappwin("boxer","imcontrolgeom",self.guiim.inspector)
		except : E2setappval("boxer","imcontrol",False)
		try:
			E2setappval("boxer","mxcontrol",self.guimx.inspector.isVisible())
			if self.guimx.inspector.isVisible() : E2saveappwin("boxer","mxcontrolgeom",self.guimx.inspector)
		except : E2setappval("boxer","mxcontrol",False)
		
		return (self.getBoxes())

	def autoboxbutton(self,bool):
		if self.autoBoxer.autoBox(self.boxable):
			self.boxDisplayUpdate()

	def toggleDynapix(self,bool):
		self.dynapix = bool
		self.autoBoxer.setMode(self.dynapix,self.anchoring)
		
	def done(self):
		self.boxable.cacheExcToDisk()
		self.dynapixp.quit
		
	def trydata(self,data,thr):
		print 'trydata was pressed, this feature is currently disabled'
		
	def optparamsupdate(self,thresh,profile,radius):
		self.guictl.updatedata(thresh,profile,radius)
		
	def setautobox(self,selmode):
		if self.autoBoxer.setSelectionMode(selmode):
			self.boxDisplayUpdate()
		
	def setprofilecmp(self,cmpmode):
		if self.autoBoxer.setCmpMode(cmpmode):
			self.boxDisplayUpdate()
			
		
	def nocupdate(self,bool):
		self.anchoring = bool
		self.mousehandlers["boxing"].setAnchoring(bool)
		self.autoBoxer.setMode(self.dynapix,self.anchoring)
		
	def classify(self,bool):
		self.boxable.classify()
		
	def erasetoggled(self,bool):
		# for the time being there are only two mouse modes
		
		if bool == True:
			self.mousehandlers["erasing"].setmode(Boxable.ERASE)
			self.mousehandler = self.mousehandlers["erasing"]
		else:
			self.guiim.addShape("eraser",EMShape(["circle",0,0,0,0,0,0,0.1]))
			self.updateImageDisplay()
			self.mousehandler = self.mousehandlers["boxing"]
			
	
	def unerasetoggled(self,bool):
		if bool == True:
			self.mousehandlers["erasing"].setmode(Boxable.UNERASE)
			self.mousehandler = self.mousehandlers["erasing"]
		else:
			self.guiim.addShape("eraser",EMShape(["circle",0,0,0,0,0,0,0.1]))
			self.updateImageDisplay()
			self.mousehandler = self.mousehandlers["boxing"]
	
	def updateEraseRad(self,rad):
		self.mousehandlers["erasing"].setEraseRadius(rad)

	def quit(self):
		self.dynapixp.quit()
	
	def setdummybox(self,box):
		self.autoBoxer.setDummyBox(box)
		self.boxDisplayUpdate()
		
	def setnonedummy(self):
		self.autoBoxer.setDummyBox(None)
		self.boxDisplayUpdate()
		
	def writeboxesimages(self):
		self.boxable.cacheExcToDisk()
		for imagename in self.imagenames:
			
			
			boxable = Boxable(imagename,self,self.autoBoxer)
		
			mode = self.autoBoxer.getMode()
			self.autoBoxer.setModeExplicit(SwarmAutoBoxer.COMMANDLINE)
			self.autoBoxer.autoBox(boxable)
			self.autoBoxer.setModeExplicit(mode)
			
			boxable.writeboximages()
	
	def writeboxesdbs(self):
		self.boxable.cacheExcToDisk()
		for imagename in self.imagenames:
			
			try:
				projectdb = EMProjectDB()
				data = projectdb[imagename+"_autoboxer"]
				trimAutoBoxer = projectdb[data["auto_boxer_db_name"]]
				autoBoxer = SwarmAutoBoxer(self)
				autoBoxer.become(trimAutoBoxer)
				print "writing box coordinates for",imagename,"using",data["auto_boxer_db_name"]
			except:
				autoBoxer = self.autoBoxer
				print "writing box coordinates for",imagename,"using currently stored autoboxer"
			boxable = Boxable(imagename,self,autoBoxer)
			
			mode = autoBoxer.getMode()
			autoBoxer.setModeExplicit(SwarmAutoBoxer.COMMANDLINE)
			autoBoxer.autoBox(boxable)
			autoBoxer.setModeExplicit(mode)
			
			boxable.writedb()
	
	def center(self,technique):
		
		if self.boxable.center(technique):
			self.boxDisplayUpdate()
		else:
			print 'technique',technique,'is unsupported - check back tomorrow'
			
	def toggleFrozen(self,unusedbool):
		self.boxable.toggleFrozen()
		self.guiim.setFrozen(self.boxable.isFrozen())
		self.updateImageDisplay()

	def setBoxingMethod(self,ref,manual):
		'''
		Okay could do it with one argument but leaving it this way makes it obvious
		'''
		
		if ref and manual:
			print 'error, you cant set both ref and manual'
			
		if ref:
			self.mousehandlers["boxing"].setmode(GUIbox.REFERENCE_ADDING)
		elif manual:
			self.mousehandlers["boxing"].setmode(GUIbox.MANUALLY_ADDING)
	
	def autoBoxerDBChanged(self):
		self.guictl.updateABTable()

	def addNewAutoBoxerDB(self):
		autoBoxer = SwarmAutoBoxer(self)
		autoBoxer.boxsize = self.boxsize
		autoBoxer.setMode(self.dynapix,self.anchoring)
		autoboxerdbstring = "autoboxer_"+autoBoxer.getCreationTS()
		trimAutoBoxer = TrimSwarmAutoBoxer(autoBoxer)
		projectdb = EMProjectDB()
		projectdb[autoboxerdbstring] = trimAutoBoxer
	
	def addCopyAutoBoxerDB(self,timestamp):
		#print "adding a copy of the ab with ts",timestamp
		projectdb = EMProjectDB()
		trimAutoBoxer = copy(projectdb["autoboxer_"+timestamp])
		trimAutoBoxer.setCreationTS(gm_time_string())
		autoboxerdbstring = "autoboxer_"+trimAutoBoxer.getCreationTS()
		#print "adding string to db",autoboxerdbstring
		projectdb[autoboxerdbstring] = trimAutoBoxer

class AutoBoxerSelectionsMediator:
	'''
	A class for coordinating the GUIboxPanel and the the GUIbox in relation
	to adding and removing AutoBoxers, and changing which Boxables use which
	AutoBoxer etc
	'''
	def __init__(self,parent):
		if not isinstance(parent,GUIbox):
			print "error, the AutoBoxerSelectionsMediator must be initialized with a GUIbox type as its first constructor argument"
			return
		self.parent=parent
		self.currentimagename = None
		self.imagenames = []
		
		self.setCurrentImageName(parent.getCurrentImageName())
		self.setImageNames(parent.getImageNames())
		
	def setCurrentImageName(self,imagename):
		self.currentimagename = imagename
		
	def setImageNames(self,imagenames):
		self.imagenames = imagenames
		
	def getAutoBoxerData(self):
		projectdb = EMProjectDB()
		dictdata = {}
		
		for i in projectdb.items():
			if i[0][0:10] == "autoboxer_":
				dictdata[i[0]] = []
		
		for imagename in self.imagenames:
			found = False
			try:
				data = projectdb[imagename+"_autoboxer"]
				#trimAutoBoxer = projectdb[data["auto_boxer_db_name"]]
				found = True
			except: pass
			
			if found:
				dictdata[data["auto_boxer_db_name"]].append(imagename)
		
		return dictdata
	
	def getCurrentAutoBoxerTS(self):
		return self.parent.getCurrentAutoBoxerTS()

	def addNewAutoBoxer(self):
		self.parent.addNewAutoBoxerDB()

	def addCopyAutoBoxer(self,timestamp):
		self.parent.addCopyAutoBoxerDB(timestamp)
		
		
class GUIboxPanel(QtGui.QWidget):
	def __init__(self,target,abselmediator) :
		
		QtGui.QWidget.__init__(self,None)
		self.target=target
		self.abselmediator = abselmediator
		
		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		
		self.tabwidget = QtGui.QTabWidget(self)
		self.insertMainTab()
		self.insertAdvancedTab()
		self.vbl.addWidget(self.tabwidget)
		
		self.dummybox = Box()
		self.dummybox.isanchor = False
		self.dummybox.isdummy = True
		self.currentlyselected = -1 # used in the abtable
		
		self.lock = False
		self.connect(self.bs,QtCore.SIGNAL("editingFinished()"),self.newBoxSize)
	
		self.connect(self.thr,QtCore.SIGNAL("valueChanged"),self.newThresh)
		self.connect(self.done,QtCore.SIGNAL("clicked(bool)"),self.target.quit)
		self.connect(self.writeboxesimages,QtCore.SIGNAL("clicked(bool)"),self.target.writeboxesimages)
		self.connect(self.writeboxesdbs,QtCore.SIGNAL("clicked(bool)"),self.target.writeboxesdbs)
		self.connect(self.classifybut,QtCore.SIGNAL("clicked(bool)"),self.target.classify)
		self.connect(self.dynapick,QtCore.SIGNAL("clicked(bool)"),self.dynapickd)
		self.connect(self.trythat,QtCore.SIGNAL("clicked(bool)"),self.trythatd)
		self.connect(self.reset,QtCore.SIGNAL("clicked(bool)"),self.target.setnonedummy)
		self.connect(self.thrbut, QtCore.SIGNAL("clicked(bool)"), self.gboxclick)
		self.connect(self.selbut, QtCore.SIGNAL("clicked(bool)"), self.gboxclick)
		self.connect(self.morselbut, QtCore.SIGNAL("clicked(bool)"), self.gboxclick)
		self.connect(self.ratiobut, QtCore.SIGNAL("clicked(bool)"), self.cmpboxclick)
		#self.connect(self.centerbutton,QtCore.SIGNAL("clicked(bool)"),self.centerpushed)
		self.connect(self.difbut, QtCore.SIGNAL("clicked(bool)"), self.cmpboxclick)
		self.connect(self.nocpick, QtCore.SIGNAL("clicked(bool)"), self.target.nocupdate)
		
#		self.target.connect(self.target,QtCore.SIGNAL("nboxes"),self.nboxesChanged)
	
	#def centerpushed(self,unused):
		#self.target.center(str(self.centerooptions.currentText()))
	
	def insertMainTab(self):
		# this is the box layout that will store everything
		self.main_inspector = QtGui.QWidget()
		self.main_vbl =  QtGui.QVBoxLayout(self.main_inspector)
		
		self.infohbl = QtGui.QHBoxLayout()
		self.info = QtGui.QLabel("%d Boxes"%len(self.target.getBoxes()),self)
		self.ppc = QtGui.QLabel("%f particles per click"%0,self)
		self.infohbl.addWidget(self.info)
		self.infohbl.addWidget(self.ppc)
		
		
		self.statsbox = QtGui.QGroupBox("Stats")
		self.statsbox.setLayout(self.infohbl)
		self.main_vbl.addWidget(self.statsbox)
		
		self.boxingvbl = QtGui.QVBoxLayout()
		
		self.boxinghbl1=QtGui.QHBoxLayout()
		self.boxinghbl1.setMargin(0)
		self.boxinghbl1.setSpacing(2)
		
		self.refbutton=QtGui.QPushButton("Reference")
		self.refbutton.setCheckable(1)
		self.refbutton.setChecked(True)
		self.boxinghbl1.addWidget(self.refbutton)
		
		self.manualbutton=QtGui.QPushButton("Manual")
		self.manualbutton.setCheckable(1)
		self.manualbutton.setChecked(False)
		self.boxinghbl1.addWidget(self.manualbutton)
		
		self.lblbs=QtGui.QLabel("Box Size:",self)
		self.boxinghbl1.addWidget(self.lblbs)
		self.bs = QtGui.QLineEdit(str(self.target.boxsize),self)
		self.boxinghbl1.addWidget(self.bs)
		
		self.boxingvbl.addLayout(self.boxinghbl1)
	
		self.boxinghbl3=QtGui.QHBoxLayout()
		self.dynapick = QtGui.QCheckBox("Dynapix")
		self.dynapick.setChecked(self.target.dynapix)
		self.boxinghbl3.addWidget(self.dynapick)
		self.nocpick = QtGui.QCheckBox("Anchor")
		self.nocpick.setChecked(self.target.anchoring)
		self.boxinghbl3.addWidget(self.nocpick)
		self.autobox=QtGui.QPushButton("Auto Box")
		self.boxinghbl3.addWidget(self.autobox)
		self.togfreeze=QtGui.QPushButton("Toggle Freeze")
		self.boxinghbl3.addWidget(self.togfreeze)
		self.boxingvbl.addLayout(self.boxinghbl3)
	
		self.boxinghbl2=QtGui.QHBoxLayout()
		self.boxinghbl2.setMargin(2)
		self.boxinghbl2.setSpacing(6)
		#self.vbl.addLayout(self.hbl1)
		
		self.erasepic = QtGui.QIcon("/home/d.woolford/erase.png");
		self.erase=QtGui.QPushButton(self.erasepic,Boxable.ERASE)
		self.erase.setCheckable(1)
		self.boxinghbl2.addWidget(self.erase)
		
		self.unerasepic = QtGui.QIcon("/home/d.woolford/erase.png");
		self.unerase=QtGui.QPushButton(self.unerasepic,Boxable.UNERASE)
		self.unerase.setCheckable(1)
		self.boxinghbl2.addWidget(self.unerase)
		
		self.eraseradtext=QtGui.QLabel("Erase Radius",self)
		self.boxinghbl2.addWidget(self.eraseradtext)
		
		self.eraserad = QtGui.QLineEdit(str(self.target.eraseradius),self)
		self.boxinghbl2.addWidget(self.eraserad)
		self.eraserad.setEnabled(False)
		
		self.boxingvbl.addLayout(self.boxinghbl2)
		
		self.interactiveboxing = QtGui.QGroupBox("Interactive Boxing")
		self.interactiveboxing.setLayout(self.boxingvbl)
		self.main_vbl.addWidget(self.interactiveboxing)
		
		self.lock = True
		self.autoboxerhdbl = QtGui.QHBoxLayout()
		# ab means autoboxer
		self.abtable = QtGui.QTableWidget(1,2,self)
		self.abtable.setColumnWidth(1,150)
		self.abcol0title = QtGui.QTableWidgetItem("Autoboxer ID")
		self.abcol1title = QtGui.QTableWidgetItem("Associated Images")
		self.updateABTable()
		self.lock = True
		self.abtable.setHorizontalHeaderItem(0,self.abcol0title)
		self.abtable.setHorizontalHeaderItem(1,self.abcol1title)
		self.autoboxerhdbl.addWidget(self.abtable)
		self.lock = False
		
		self.autoboxervbl1 = QtGui.QVBoxLayout()
		self.abcopy = QtGui.QPushButton("Copy")
		self.autoboxervbl1.addWidget(self.abcopy)
		self.abnew = QtGui.QPushButton("New")
		self.autoboxervbl1.addWidget(self.abnew)
		self.abdelete = QtGui.QPushButton("Delete")
		self.autoboxervbl1.addWidget(self.abdelete)
		self.autoboxerhdbl.addLayout(self.autoboxervbl1)
		
		self.abmanagement = QtGui.QGroupBox("Auto Boxer Management")
		self.abmanagement.setLayout(self.autoboxerhdbl)
		self.main_vbl.addWidget(self.abmanagement)
		
		# output
		self.outputvbl = QtGui.QVBoxLayout()
		self.outputhbl1=QtGui.QHBoxLayout()
		self.writeboxesimages = QtGui.QPushButton("Write Box Images")
		self.outputhbl1.addWidget(self.writeboxesimages)
		self.writeboxesdbs = QtGui.QPushButton("Write Coord Files")
		self.outputhbl1.addWidget(self.writeboxesdbs)
		self.outputvbl.addLayout(self.outputhbl1)
		
		self.outputhbl2=QtGui.QHBoxLayout()
		self.usingboxsizetext=QtGui.QLabel("Using Box Size:",self)
		self.outputhbl2.addWidget(self.usingboxsizetext)
		self.usingboxsize = QtGui.QLineEdit(str(self.target.boxsize),self)
		self.outputhbl2.addWidget(self.usingboxsize)
		self.outputvbl.addLayout(self.outputhbl2)
		
		self.outputbox = QtGui.QGroupBox("Output")
		self.outputbox.setLayout(self.outputvbl)
		self.main_vbl.addWidget(self.outputbox)

		self.classifybut=QtGui.QPushButton("Classify")
		self.main_vbl.addWidget(self.classifybut)
		
		self.done=QtGui.QPushButton("Done")
		self.main_vbl.addWidget(self.done)
		
		self.tabwidget.addTab(self.main_inspector,"Main")
		
		self.connect(self.eraserad,QtCore.SIGNAL("editingFinished()"),self.updateEraseRad)
		self.connect(self.erase, QtCore.SIGNAL("clicked(bool)"), self.erasetoggled)
		self.connect(self.unerase, QtCore.SIGNAL("clicked(bool)"), self.unerasetoggled)
		self.connect(self.autobox,QtCore.SIGNAL("clicked(bool)"),self.target.autoboxbutton)
		self.connect(self.togfreeze,QtCore.SIGNAL("clicked(bool)"),self.target.toggleFrozen)
		
		self.connect(self.refbutton, QtCore.SIGNAL("clicked(bool)"), self.refbuttontoggled)
		self.connect(self.manualbutton, QtCore.SIGNAL("clicked(bool)"), self.manualbuttontoggled)
		self.connect(self.abnew, QtCore.SIGNAL("clicked(bool)"), self.addNewAutoBoxer)
		self.connect(self.abcopy, QtCore.SIGNAL("clicked(bool)"), self.addCopyAutoBoxer)
		self.connect(self.abtable, QtCore.SIGNAL("itemChanged(QtGui.QTableWidgetItem)"), self.abtableItemChanged)
		self.connect(self.abtable, QtCore.SIGNAL("cellChanged(int,int)"), self.abtableCellChanged)
	
	def abtableCellChanged(self,i,j):
		#print i,j
		if self.lock:
			#print "locked,returning"
			return
		if j == 0: # we're in the first row, something happened, maybe a check change
			self.lock = True
			#try:
			if i == self.currentlychecked:
				# just make sure the check stays on
				self.col1[self.currentlychecked].setCheckState(Qt.Checked)
			elif self.col1[i].checkState() == Qt.Checked:
				# uncheck the previously selected one
				self.col1[self.currentlychecked].setCheckState(Qt.Unchecked)
				self.col1[i].setCheckState(Qt.Checked)
				self.currentlychecked = i
				self.target.changeCurrentAutoBoxer(str(self.col1[i].text()))
			else:
				print "error, I don't know what's happening"
			#except: pass
			self.lock = False
				
	def abtableItemChanged(self,item):
		print "item changed"
	
	def updateABTable(self):
		self.lock = True
		data = self.abselmediator.getAutoBoxerData()
		self.col1 = []
		self.col2 = []
		idx = 0
		for d in data.items():
			if idx >= self.abtable.rowCount():
				self.abtable.insertRow(idx)
			col1 = QtGui.QTableWidgetItem(d[0][10:])
			qstr =''
			for i,s in enumerate(d[1]):
				qstr += s
				if i != len(d[1])-1:
					qstr += ', '
			col2 = QtGui.QTableWidgetItem(qstr)
			
			flag1 = Qt.ItemFlags(Qt.ItemIsUserCheckable)
			flag2 = Qt.ItemFlags(Qt.ItemIsSelectable)
			flag3 = Qt.ItemFlags(Qt.ItemIsEnabled)
			flag4 = Qt.ItemFlags(Qt.ItemIsEditable)
			#flags = flags.
			col1.setFlags(flag1|flag2|flag3|flag4) #Qt.ItemIsEnabled+Qt.ItemIsUserCheckable)) #&Qt.ItemIsSelectable))
			col1.setCheckState( Qt.Unchecked)
			col2.setFlags(flag3)
			self.abtable.setItem(idx,0,col1)
			self.abtable.setItem(idx,1,col2)
			self.col1.append(col1)
			self.col2.append(col2)
			idx += 1
		
		# remove any overhanging columns if they already existed
		while (len(data) < self.abtable.rowCount()):
			self.abtable.removeRow(self.abtable.rowCount()-1)
	
		currentsel = self.abselmediator.getCurrentAutoBoxerTS()
		self.currentlychecked = -1
		for i,col in enumerate(self.col1):
			if str(col.text()) == currentsel:
				col.setCheckState( Qt.Checked)
				self.currentlychecked = i
				break
		
		self.lock = False
	
	def addNewAutoBoxer(self,bool):
		self.abselmediator.addNewAutoBoxer()
		self.updateABTable()
		
	def addCopyAutoBoxer(self,bool):
		numsel = 0
		selected = ""
		for col in self.col1:
			if self.abtable.isItemSelected(col):
				numsel += 1
				selected = str(col.text())
				if numsel > 1:
					print "error, more than one autoboxer is selected. Please choose only one"
					return
		if numsel == 0:
			print "no autoboxers were selected, doing nothing"
			return
		else:
			self.abselmediator.addCopyAutoBoxer(selected)
			self.updateABTable()
	
	def insertAdvancedTab(self):
		# this is the box layout that will store everything
		self.adv_inspector = QtGui.QWidget()
		self.advanced_vbl =  QtGui.QVBoxLayout(self.adv_inspector)
		
		#  Insert the plot widget
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
		
		self.advanced_vbl2 = QtGui.QVBoxLayout()
		
		self.advanced_vbl2.addLayout(self.plothbl)
		
		self.thr = ValSlider(self,(0.0,3.0),"Threshold:")
		self.thr.setValue(1.0)
		self.advanced_vbl2.addWidget(self.thr)
		
		
		self.interbox = QtGui.QGroupBox("Interactive Parameters")
		self.interbox.setLayout(self.advanced_vbl2)
		self.advanced_vbl.addWidget(self.interbox)
		
		self.thrbut = QtGui.QRadioButton(SwarmAutoBoxer.THRESHOLD)
		self.selbut = QtGui.QRadioButton(SwarmAutoBoxer.SELECTIVE)
		self.selbut.setChecked(True)
		self.morselbut = QtGui.QRadioButton(SwarmAutoBoxer.MORESELECTIVE)
		
		self.methodhbox = QtGui.QHBoxLayout()
		self.methodhbox.addWidget(self.thrbut)
		self.methodhbox.addWidget(self.selbut)
		self.methodhbox.addWidget(self.morselbut)
		
		self.groupbox = QtGui.QGroupBox("Auto Box Method")
		self.groupbox.setLayout(self.methodhbox)
		
		self.advanced_vbl.addWidget(self.groupbox)

		self.ratiobut = QtGui.QRadioButton("Ratio")
		self.ratiobut.setChecked(True)
		self.difbut = QtGui.QRadioButton("Difference")
		
		self.cmpmethodhbox = QtGui.QHBoxLayout()
		self.cmpmethodhbox.addWidget(self.ratiobut)
		self.cmpmethodhbox.addWidget(self.difbut)
		
		self.cmpgroupbox = QtGui.QGroupBox("Peak Profile Comparitor")
		self.cmpgroupbox.setLayout(self.cmpmethodhbox)
		
		self.advanced_vbl.addWidget(self.cmpgroupbox)

		
		self.tabwidget.addTab(self.adv_inspector,"Advanced")

	def refbuttontoggled(self,bool):
		
		if self.refbutton.isChecked():
			self.manualbutton.setChecked(False)
		
		if not self.refbutton.isChecked():
			self.manualbutton.setChecked(True)

		self.target.setBoxingMethod(self.refbutton.isChecked(),self.manualbutton.isChecked())
		
	def manualbuttontoggled(self,bool):
		if self.manualbutton.isChecked():
			self.refbutton.setChecked(False)
		
		if not self.manualbutton.isChecked():
			self.refbutton.setChecked(True)
			
		self.target.setBoxingMethod(self.refbutton.isChecked(),self.manualbutton.isChecked())

	def erasetoggled(self,bool):
		self.unerase.setChecked(False)
		self.eraserad.setEnabled(bool)
		self.target.guiim.setMouseTracking(bool)
		self.target.erasetoggled(bool)
		
	def unerasetoggled(self,bool):
		self.erase.setChecked(False)
		self.eraserad.setEnabled(bool)
		self.target.guiim.setMouseTracking(bool)
		self.target.unerasetoggled(bool)
		
	def dynapickd(self,bool):
		self.target.toggleDynapix(bool)
	
	def cmpboxclick(self,unusedbool):
		if self.ratiobut.isChecked():
			s = BoxingTools.CmpMode.SWARM_RATIO
		elif self.difbut.isChecked():
			s = BoxingTools.CmpMode.SWARM_DIFFERENCE
		else:
			print "Bug intercepted in e2boxer.py. Please email the development team."
			
		self.target.setprofilecmp(s)
	
	def gboxclick(self,unusedbool):
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
		self.dummybox.optprofile = self.window.getData()
		self.dummybox.correlationscore = float(self.thr.getValue())
		self.target.setdummybox(self.dummybox)
	
	def updatedata(self,thresh,data,datar):
		#print data
		self.window.setData(data,datar)
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
