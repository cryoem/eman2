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
	parser.add_option("--auto","-A",type="string",action="append",help="Autobox using specified method: ref, grid, db",default=[])
	parser.add_option("--writecoords",action="store_true",help="Write data box files",default=False)
	parser.add_option("--writeboximages",action="store_true",help="Write data box files",default=False)
	parser.add_option("--force","-f",action="store_true",help="Force overwrites old files",default=False)
	parser.add_option("--overlap",type="int",help="(auto:grid) number of pixels of overlap between boxes. May be negative.")
	parser.add_option("--nretest",type="int",help="(auto:ref) Number of reference images (starting with the first) to use in the final test for particle quality.",default=-1)
	#parser.add_option("--retestlist",type="string",help="(auto:ref) Comma separated list of image numbers for retest cycle",default="")
	#parser.add_option("--farfocus",type="string",help="filename or 'next', name of an aligned far from focus image for preliminary boxing",default=None)
	parser.add_option("--parallel",type="int",help="specify more than one processor",default=1)
	parser.add_option("--mergedb",action="store_true",help="A special argument, if true all input arguments are considered to be box files and they are merged into the project database as manually selected particles",default=False)
	
	(options, args) = parser.parse_args()
	if len(args)<1 : parser.error("Input image required")
	
	if not options.gui and not options.auto:
		parser.error("Atleast one of the --gui or --auto arguments are required.")
		exit(1)
	
	filenames = []
	for i in range(0,len(args)):
		filenames.append(args[i])

	logid=E2init(sys.argv)
	
	
	if options.mergedb == True:
		#The user wants to add some boxes to the database
		mergedb(filenames)
		sys.exit(1)
	
	# we need to know how big to make the boxes. If nothing is specified, but
	# reference particles are, then we use the reference particle size
	if options.boxsize<5 :
		if not options.boxsize in good_box_sizes:
			print "Note: EMAN2 processing would be more efficient with a boxsize of %d"%good_boxsize(options.boxsize)
	
	boxes=[]
	if len(options.auto)>0 :
		print "Autobox mode ",options.auto[0]
	
		if "db" in options.auto:
			print "auto data base boxing"
		
			autoboxmulti(filenames,options)
			#if len(filenames) == 1:
				#autoboxsingle(filenames[0],options)
				#exit(1)
			#else:
				#print "autoboxing using parallelism - you specified",options.parallel,"processors"
				#autoboxer = AutoDBBoxer(filenames,options.parallel,options,options.force)
				#try:
					#from emimage import get_app
				#except: 
					#print "error, can't import get_app"
					#exit(1)
					
				#a = get_app()
				#autoboxer.go(a)
				#a.exec_()
				#print "done"
					
			exit(1)
		elif "grid" in options.auto:
			image_size=gimme_image_dimensions2D(filenames[0])
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
		else:
			print "unknown autoboxing method:",options.auto
			exit(1)

	E2end(logid)

	# invoke the GUI if requested
	if options.gui:
		gui=GUIbox(filenames,boxes,options.boxsize)
		gui.run()
		
	print "Exiting e2boxer"
	
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
	'''
	Merges a set of .box files into the local database - stores them as manual boxes
	'''
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
			manualboxes = getKeyEntryIDD(filename,"manualboxes")
		except:
			manualboxes = []
	
		manualboxes.extend(boxes)
		setKeyEntryIDD(filename,"manualboxes",manualboxes)


def autoboxmulti(imagenames,options):
	projectdb = EMProjectDB()
	for imagename in imagenames:
		print "autoboxing",imagename
		
		try:
			data = projectdb[getIDDKey(imagename)]
			trimAutoBoxer = projectdb[data["auto_boxer_unique_id"]]["autoboxer"]
			autoBoxer = SwarmAutoBoxer(None)
			autoBoxer.become(trimAutoBoxer)
			print 'using cached autoboxer db'
		except:
			try:
				print "using most recent autoboxer"
				trimAutoBoxer = projectdb["currentautoboxer"]
				autoBoxer = SwarmAutoBoxer(None)
				autoBoxer.become(trimAutoBoxer)
			except:
				print "Error - there seems to be no autoboxing information in the database - autobox interactively first - bailing"
				continue
		
		boxable = Boxable(imagename,None,autoBoxer)
		
		if boxable.isExcluded():
			print "Image",imagename,"is excluded and being ignored"
			continue
		
		autoBoxer.setModeExplicit(SwarmAutoBoxer.COMMANDLINE)
		# Tell the boxer to delete non refs - FIXME - the uniform appraoch needs to occur - see SwarmAutoBoxer.autoBox
		autoBoxer.autoBox(boxable,False)
		if options.writecoords:
			boxable.writecoords(-1,options.force)
		if options.writeboximages:
			boxable.writeboximages(-1,options.force)
	
	
	projectdb.close()

def autoboxsingle(imagename,options):
	
	projectdb = EMProjectDB()
	try:
		data = projectdb[getIDDKey(imagename)]
		trimAutoBoxer = projectdb[data["auto_boxer_unique_id"]]["autoboxer"]
		autoBoxer = SwarmAutoBoxer(None)
		autoBoxer.become(trimAutoBoxer)
		print 'using cached autoboxer db'
	except:
		try:
			trimAutoBoxer = projectdb["currentautoboxer"]
			autoBoxer = SwarmAutoBoxer(None)
			autoBoxer.become(trimAutoBoxer)
		except:
			print "Error - there seems to be no autoboxing information in the database - autobox interactively first - bailing"
			projectdb.close()
			return 0
	
	boxable = Boxable(imagename,None,autoBoxer)
	if boxable.isExcluded():
		print "Image",imagename,"is excluded and being ignored"
		return
	
	autoBoxer.setModeExplicit(SwarmAutoBoxer.COMMANDLINE)
	# Tell the boxer to delete non refs - FIXME - the uniform appraoch needs to occur - see SwarmAutoBoxer.autoBox
	autoBoxer.autoBox(boxable,False)
	if options.writecoords:
		print "writing box coordinates"
		boxable.writecoords(-1,options.force)
	if options.writeboximages:
		print "writing boxed images"
		boxable.writeboximages(-1,options.force)
	
	projectdb.close()
	return 1
	
class AutoDBBoxer(QtCore.QObject):
	'''
	A class for managing the process of spawning many instances of e2boxer singlefile.mrc --auto=db, using parallelism
	If one CPU is specified then it still works. Basically the approach is to spawn the number of the processors, then once
	the process is finished the signal is intercepted and a new process is executed etc.
	'''
	def __init__(self,imagenames,nproc,options,force=False):
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
		self.options = options
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
		if self.options.writecoords != False:
			args.append("--writecoords")
		if self.options.writeboximages != False:
			args.append("--writeboximages")
		if self.options.force != False:
			args.append("--force")
			
		
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
	A class that knows how to add, move and remove reference and non reference boxes 
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
	'''
	This class needs a clean up and needs some comments to be inserted
	'''
	REFERENCE_ADDING = 0
	ERASING = 1
	MANUALLY_ADDING = 2
	def __init__(self,imagenames,boxes,boxsize=-1):
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
			data = projectdb[getIDDKey(self.imagenames[0])]
			autoBoxerID = data["auto_boxer_unique_id"]
			trimAutoBoxer = projectdb[autoBoxerID]["autoboxer"]
			self.autoboxername = autoBoxerID
			self.autoBoxer = SwarmAutoBoxer(self)
			self.autoBoxer.become(trimAutoBoxer)
			self.autoBoxer.setMode(self.dynapix,self.anchoring)
			print 'using cached autoboxer db'
		except:
			self.autoBoxer = SwarmAutoBoxer(self)
			self.autoBoxer.boxsize = boxsize
			self.autoBoxer.setMode(self.dynapix,self.anchoring)
			print "loaded a new autoboxer"
		
		self.boxable = Boxable(self.imagenames[0],self,self.autoBoxer)
		self.boxable.addnonrefs(boxes)
		self.boxable.boxsize = boxsize
		
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
		self.guiim.setExcluded(self.boxable.isExcluded())
		
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
		self.guictl.setImageQuality(self.boxable.getQuality())
		
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
		self.autoBoxer.autoBox(self.boxable,False)
		self.boxDisplayUpdate()
		
	
	def changeCurrentAutoBoxer(self, autoboxerid,autobox=True):
		#print "change current autoboxer"
		projectdb = EMProjectDB()
		trimAutoBoxer = projectdb[autoboxerid]["autoboxer"]
		#print "changing autoboxer to autoboxer_",timestamp,"and its stamp in the db is",trimAutoBoxer.getCreationTS()
		self.autoBoxer = SwarmAutoBoxer(self)
		self.autoBoxer.become(trimAutoBoxer)
		self.autoBoxer.writeSpecificReferencesToDB(self.boxable.getImageName())
		self.autoBoxer.setMode(self.dynapix,self.anchoring)

		self.boxable.setAutoBoxer(self.autoBoxer)
		if autobox:
			print "change current plus autobox"
			self.ptcl = []
			self.guiim.delShapes()
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
		
		if len(refboxes) != 0: 
			val = self.autoBoxer.removeReference(refboxes)
			if val == 2:
				self.boxable.clearAndCache(True)
				self.clearDisplays()
				return # avoid unecessary display update below
			
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
		if not self.boxable.isInteractive():
			return
		
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
		if not self.boxable.isInteractive(): return
		
		box = self.delbox(boxnum)
		if not (box.isref or box.ismanual):
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
		if data != None:
			self.guimx.setData(data)
			if len(data) != 0:
				if self.guimxp == None:
					self.guimxp = EMParentWin(self.guimx)
					self.guimxp.show()
					self.guimx.setSelected(0)
	
	def clearDisplays(self):
		self.ptcl = []
		self.guiim.delShapes()
		self.guimx.setData([])
		self.boxDisplayUpdate() # - the user may still have some manual boxes...
	
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
				data = projectdb[getIDDKey(self.imagenames[im])]
				autoBoxerID = data["auto_boxer_unique_id"]
				trimAutoBoxer = projectdb[autoBoxerID]["autoboxer"]
				self.autoboxername = autoBoxerID
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
			self.guiim.setExcluded(self.boxable.isExcluded())
			self.guictl.setImageQuality(self.boxable.getQuality())
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
		if not self.boxable.isInteractive():
			return
		
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
			val = self.autoBoxer.removeReference(box)
			if val == 2:
				self.boxable.clearAndCache(True)
				self.clearDisplays()
				
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
		'''
		Hit the autobox is like hitting refresh - everything is updated
		'''
		self.boxable.clearAndCache(True)
		self.autoBoxer.writeSpecificReferencesToDB(self.boxable.getImageName())
		self.boxable.getReferencesFromDB()
		self.autoBoxer.autoBox(self.boxable, True,True)
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
		
	def writeboxesimages(self,boxsize,forceoverwrite=False,imageformat="hdf"):
		self.boxable.cacheExcToDisk()
		for imagename in self.imagenames:
			
			
			try:
				projectdb = EMProjectDB()
				data = projectdb[getIDDKey(imagename)]
				trimAutoBoxer = projectdb[data["auto_boxer_unique_id"]]["autoboxer"]
				autoBoxer = SwarmAutoBoxer(self)
				autoBoxer.become(trimAutoBoxer)
				#print "writing box images for",imagename,"using",data["auto_boxer_unique_id"]
			except:
				autoBoxer = self.autoBoxer
				#print "writing box images or",imagename,"using currently stored autoboxer"
				
			boxable = Boxable(imagename,self,autoBoxer)
			if boxable.isExcluded():
				print "Image",imagename,"is excluded and being ignored"
				continue
			
			mode = self.autoBoxer.getMode()
			autoBoxer.setModeExplicit(SwarmAutoBoxer.COMMANDLINE)
			autoBoxer.autoBox(boxable)
			autoBoxer.setModeExplicit(mode)
			
			boxable.writeboximages(boxsize,forceoverwrite,imageformat)
	
	def writeboxesdbs(self,boxsize,forceoverwrite=False):
		self.boxable.cacheExcToDisk()
		for imagename in self.imagenames:
			
			try:
				projectdb = EMProjectDB()
				data = projectdb[getIDDKey(imagename)]
				trimAutoBoxer = projectdb[data["auto_boxer_unique_id"]]["autoboxer"]
				autoBoxer = SwarmAutoBoxer(self)
				autoBoxer.become(trimAutoBoxer)
				#print "writing box coordinates for",imagename,"using",data["auto_boxer_unique_id"]
			except:
				autoBoxer = self.autoBoxer
				#print "writing box coordinates for",imagename,"using currently stored autoboxer"
			boxable = Boxable(imagename,self,autoBoxer)
			
			if boxable.isExcluded():
				print "Image",imagename,"is excluded and being ignored"
				continue
			
			mode = autoBoxer.getMode()
			autoBoxer.setModeExplicit(SwarmAutoBoxer.COMMANDLINE)
			autoBoxer.autoBox(boxable)
			autoBoxer.setModeExplicit(mode)
			
			boxable.writecoords(boxsize,forceoverwrite)
	
	def center(self,technique):
		
		if self.boxable.center(technique):
			self.boxDisplayUpdate()
		else:
			print 'technique',technique,'is unsupported - check back tomorrow'
			
	def toggleFrozen(self):
		if self.boxable.isExcluded() : return
		self.boxable.toggleFrozen()
		self.boxable.writeToDB()
		self.guiim.setFrozen(self.boxable.isFrozen())
		if not self.boxable.isFrozen():
			self.changeCurrentAutoBoxer(self.boxable.getAutoBoxerID(),False)
		else:	
			self.updateImageDisplay()
		
	def changeImageQuality(self,val):
		self.boxable.setQuality(val)
		if val == Boxable.EXCLUDE:
			self.boxable.setFrozen(False) # If it was already frozen then setting to excluded overrides this
			self.guiim.setExcluded(True)
			self.guiim.setFrozen(False)
			self.boxable.clearAndCache(True) # tell boxable to clear its autoboxes and references -
			self.clearDisplays() # tell the display to clear itself
		elif self.guiim.setExcluded(False):
			self.updateImageDisplay()
			
		self.boxable.writeToDB() # make sure the infromation changes that just occured are written to the DB
		

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

	def addNewAutoBoxerDB(self, n):
		if not self.boxable.isInteractive():
			return None
		autoBoxer = SwarmAutoBoxer(self)
		autoBoxer.boxsize = self.boxsize
		autoBoxer.setMode(self.dynapix,self.anchoring)
		autoboxerdbstring = "autoboxer_"+autoBoxer.getCreationTS()
		trimAutoBoxer = TrimSwarmAutoBoxer(autoBoxer)
		convenience_name = "New " + str(n)
		trimAutoBoxer.setConvenienceName(convenience_name)
		trimAutoBoxer.writeToDB()
		return convenience_name
	
	def addCopyAutoBoxerDB(self,autoboxerid,n):
		#print "adding a copy of the ab with id is",autoboxerid
		projectdb = EMProjectDB()
		trimAutoBoxer = copy(projectdb[autoboxerid]["autoboxer"])
		trimAutoBoxer.setCreationTS(gm_time_string())
		convenience_name = "Copy " + str(n)
		trimAutoBoxer.setConvenienceName(convenience_name)
		trimAutoBoxer.writeToDB()
		return convenience_name
		
		#print "done"

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
		self.namemap = {}
		self.dictdata = {}
		self.setCurrentImageName(parent.getCurrentImageName())
		self.setImageNames(parent.getImageNames())
		
	def setCurrentImageName(self,imagename):
		self.currentimagename = imagename
		
	def setImageNames(self,imagenames):
		self.imagenames = imagenames
		
	def getAutoBoxerData(self):
		projectdb = EMProjectDB()
		self.dictdata = {}
		self.namemap = {}
		for i in projectdb:
			print i
		print 'DONE'
		for i in projectdb:
			try:
				if i[0:10] == "autoboxer_":
					tag = projectdb[i]["convenience_name"]
					self.dictdata[tag] = []
					self.namemap[i] = tag
			except:
				print "couldn't handle",i
		print self.dictdata
		print self.namemap
		
		for imagename in self.imagenames:
			found = False
			try:
				data = projectdb[getIDDKey(imagename)]
				#trimAutoBoxer = projectdb[data["auto_boxer_unique_id"]]
				found = True
			except: pass
			
			if found:
				try:
					self.dictdata[self.namemap[data["auto_boxer_unique_id"]]].append(strip_after_dot(imagename))
				except:
					print "error, an autoboxer has been lost, its stamp is",data["auto_boxer_unique_id"]
		return self.dictdata
	
	def getCurrentAutoBoxerTS(self):
		try:
			return self.namemap["autoboxer_"+self.parent.getCurrentAutoBoxerTS()]
		except:
			return None

	def addNewAutoBoxer(self):
		return self.parent.addNewAutoBoxerDB(self.getTotalAutoBoxers())

	def addCopyAutoBoxer(self,tag):
		autoboxerid = self.__getAutoBoxerIDFromTag(tag)
		if autoboxerid != None:
			return self.parent.addCopyAutoBoxerDB(autoboxerid,self.getTotalAutoBoxers())
		else:
			print "error, couldn't find autoboxer from tag",tag
			return None
	
	def changeCurrentAutoBoxer(self,tag):
		# FIXME - is there any way to use a bidirectional map?pyt
		autoboxerid = self.__getAutoBoxerIDFromTag(tag)
		if autoboxerid != None:
			return self.parent.changeCurrentAutoBoxer(autoboxerid)
		else:
			print "error, couldn't get autoboxerid"
			return None
				
	def updateDBConvenienceName(self,newname,oldname):
		'''
		Updates the unique name of the autoboxer in the DB
		if the name is already used then False is returned
		and the calling function should act on this to stop the
		name change, for example within a widget
		'''
		projectdb = EMProjectDB()
		autoboxerid = self.__getAutoBoxerIDFromTag(oldname)
		if autoboxerid != None:
			self.namemap.pop(autoboxerid)
			if self.__nameNotAlreadyPresent(newname):
				self.namemap[autoboxerid] = newname
				autoBoxer = projectdb[autoboxerid]["autoboxer"]
				autoBoxer.setConvenienceName(newname)
				autoBoxer.writeToDB()
				return True
			else:
				self.namemap[autoboxerid]["convenience_name"] = oldname
				return False
					
	def getTotalAutoBoxers(self):
		return len(self.namemap)
	
	def associatedImages(self,tag):
		return self.dictdata[tag]
	
	def remove(self,tag):
		autoboxerid = self.__getAutoBoxerIDFromTag(tag)
		projectdb = EMProjectDB()
		projectdb.pop(autoboxerid)
		self.dictdata.pop(tag)
	
	def __nameNotAlreadyPresent(self,name):
		for names in self.namemap.items():
			if names[1] == name:
			 	return False
		
		return True
	
	def __getAutoBoxerIDFromTag(self,tag):
		for names in self.namemap.items():
			if names[1] == tag:
				return names[0]
			
		print "error, couldn't find",tag,"in the namemap"
		return None
	
	def toggleFrozen(self,tag,bool):
		#frozen = boxable.isFrozen()
		#if frozen:
			#new_name = self.addCopyAutoBoxer(tag)
			#self.getAutoBoxerData() # FIXME this is inefficient, could just add to self.dictdata etc
			#autoboxerid = self.__getAutoBoxerIDFromTag(new_name)
			
			#boxable.setAutoBoxerID(autoboxerid)
			#boxable.writeToDB()
		#boxable.writeToDB()
		self.parent.toggleFrozen()
		
	def clearCurrent(self):
		new_name = self.addNewAutoBoxer()
		
		# if the new_name is none then the current Boxable is frozen!
		if new_name == None: return new_name
		
		self.getAutoBoxerData() # FIXME this is inefficient, could just add to self.dictdata etc
		autoboxerid = self.__getAutoBoxerIDFromTag(new_name)
		boxable = self.parent.getBoxable()
		boxable.setAutoBoxerID(autoboxerid)
		boxable.writeToDB()
		boxable.clearAndCache(True)
		self.parent.clearDisplays()
		return new_name
		
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
		self.pos_int_validator = QtGui.QIntValidator(self)
		self.pos_int_validator.setBottom(1)
		self.bs = QtGui.QLineEdit(str(self.target.boxsize),self)
		self.bs.setValidator(self.pos_int_validator)
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
		
		self.boxinghbl4=QtGui.QHBoxLayout()
		self.togfreeze=QtGui.QPushButton("Toggle Freeze")
		self.boxinghbl4.addWidget(self.togfreeze)
		self.clear=QtGui.QPushButton("Clear")
		self.boxinghbl4.addWidget(self.clear)
		
		self.imagequality=QtGui.QLabel("Image Quality:",self)
		self.boxinghbl4.addWidget(self.imagequality)
		
		self.imagequalities = QtGui.QComboBox(self)
		for metadata in Boxable.QUALITY_META_DATA:
			self.imagequalities.addItem(metadata)
		#self.imagequalities.setSelection(Boxable.AVERAGE)
		self.imagequalities.setCurrentIndex(Boxable.QUALITY_META_DATA_MAP[Boxable.AVERAGE])
		self.boxinghbl4.addWidget(self.imagequalities)
		
		self.boxingvbl.addLayout(self.boxinghbl4)
		
		self.interactiveboxing = QtGui.QGroupBox("Interactive Boxing")
		self.interactiveboxing.setLayout(self.boxingvbl)
		self.main_vbl.addWidget(self.interactiveboxing)
		
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
		self.usingboxsize.setValidator(self.pos_int_validator)
		self.outputhbl2.addWidget(self.usingboxsize)
		
		self.outputvbl.addLayout(self.outputhbl2)
		
		self.outputhbl3=QtGui.QHBoxLayout()
		
		self.outputformat=QtGui.QLabel("Image Format:",self)
		self.outputhbl3.addWidget(self.outputformat)
		
		self.outputformats = QtGui.QComboBox(self)
		self.outputformats.addItem("hdf")
		self.outputformats.addItem("img")
		self.outputhbl3.addWidget(self.outputformats)
		
		
		self.outputforceoverwrite = QtGui.QCheckBox("Force Overwrite")
		self.outputforceoverwrite.setChecked(False)
		self.outputhbl3.addWidget(self.outputforceoverwrite)
		
		self.outputvbl.addLayout(self.outputhbl3)
		
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
		self.connect(self.togfreeze,QtCore.SIGNAL("clicked(bool)"),self.toggleFrozen)
		self.connect(self.clear,QtCore.SIGNAL("clicked(bool)"),self.clearCurrent)
		
		self.connect(self.refbutton, QtCore.SIGNAL("clicked(bool)"), self.refbuttontoggled)
		self.connect(self.manualbutton, QtCore.SIGNAL("clicked(bool)"), self.manualbuttontoggled)

		self.connect(self.writeboxesimages,QtCore.SIGNAL("clicked(bool)"),self.writebimages)
		self.connect(self.writeboxesdbs,QtCore.SIGNAL("clicked(bool)"),self.writebcoords)
		
		QtCore.QObject.connect(self.imagequalities, QtCore.SIGNAL("currentIndexChanged(QString)"), self.imageQualityChanged)

	def setImageQuality(self,integer):
		self.lock = True
		self.imagequalities.setCurrentIndex(integer)
		self.lock = False
	def imageQualityChanged(self,val):
		if self.lock == False:
			self.target.changeImageQuality(str(val))

	def clearCurrent(self,unused):
		self.lock = True
		new_name = self.abselmediator.clearCurrent()
		self.lock = False
		if new_name != None: # if the boxable wasn't frozen...
			self.updateABTable()
			self.setChecked(new_name)
		
	def toggleFrozen(self,bool):
		self.lock = True
		self.abselmediator.toggleFrozen(self.col1[self.currentlyselected].text(),bool)
		self.lock = False
	
	def writebimages(self,unused):
		boxsize = int(str(self.usingboxsize.text()))
		realboxsize = int(str(self.bs.text()))
		if realboxsize == boxsize:
			boxsize = -1 # negative one is a flag that tells the boxes they don't need to be resized... all the way in the Box Class
		self.target.writeboxesimages(boxsize,self.outputforceoverwrite.isChecked(),str(self.outputformats.currentText()))
		
	def writebcoords(self,unused):
		boxsize = int(str(self.usingboxsize.text()))
		realboxsize = int(str(self.bs.text()))
		if realboxsize == boxsize:
			boxsize = -1 # negative one is a flag that tells the boxes they don't need to be resized... all the way in the Box Class
		self.target.writeboxesdbs(boxsize,self.outputforceoverwrite.isChecked())
	
	def setChecked(self,tag):
		
		for i,col in enumerate(self.col1):
			if str(col.text()) == tag:
				col.setCheckState(Qt.Checked)
				self.currentlychecked = i
				break
	
	def abtableCellChanged(self,i,j):
		if i >= len(self.col1): return
		if self.lock:
			return
		#data = self.abselmediator.getAutoBoxerData()
		#data = data.items()
		if str(self.col1[i].text()) != self.colnames[i]:
			if not self.abselmediator.updateDBConvenienceName(str(self.col1[i].text()),self.colnames[i]):
				self.col1[i].setText(self.colnames[i])
			self.lock = False
			return
		
		self.lock = False
		
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
				self.abselmediator.changeCurrentAutoBoxer(str(self.col1[i].text()))
				self.updateABTable()
			else:
				print "error, unforeseen checkstate circumstance. Nothing done"
			#except: pass
		
		self.lock = False
				
	def abtableItemChanged(self,item):
		print "item changed"
	
	def updateABTable(self):
		
		data = self.abselmediator.getAutoBoxerData()
		self.col1 = []
		self.col2 = []
		self.colnames = []
		if len(data) == 0:
			return
			
		self.lock = True
		idx = 0
		for d in data.items():
			if idx >= self.abtable.rowCount():
				self.abtable.insertRow(idx)
			col1 = QtGui.QTableWidgetItem(d[0])
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
			self.colnames.append(col1.text())
			idx += 1
		
		# remove any overhanging columns if they already existed
		while (len(data) < self.abtable.rowCount()):
			self.abtable.removeRow(self.abtable.rowCount()-1)
	
		currentsel = self.abselmediator.getCurrentAutoBoxerTS()
		self.currentlychecked = -1
		
		self.setChecked(currentsel)
		self.lock = False
	
	def deleteAutoBoxer(self,unused):
		items = self.abtable.selectedItems()
		data = self.abselmediator.getAutoBoxerData()
		update = False
		for item in items:
			if item.column() == 0:
				if len(self.abselmediator.associatedImages(str(item.text()))) != 0:
					print "can't delete an autoboxer unless it has no images associated with it"
				else:
					self.abselmediator.remove(str(item.text()))
					update = True
					
		if update:
			self.updateABTable()
	
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
			self.lock=False
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
		self.advanced_vbl.addWidget(self.abmanagement)
		
		
		
		self.tabwidget.addTab(self.adv_inspector,"Advanced")
		
		
		self.connect(self.abnew, QtCore.SIGNAL("clicked(bool)"), self.addNewAutoBoxer)
		self.connect(self.abcopy, QtCore.SIGNAL("clicked(bool)"), self.addCopyAutoBoxer)
		self.connect(self.abdelete, QtCore.SIGNAL("clicked(bool)"), self.deleteAutoBoxer)
		self.connect(self.abtable, QtCore.SIGNAL("itemChanged(QtGui.QTableWidgetItem)"), self.abtableItemChanged)
		self.connect(self.abtable, QtCore.SIGNAL("cellChanged(int,int)"), self.abtableCellChanged)
		

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
