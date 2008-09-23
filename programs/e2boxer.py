#!/usr/bin/env python
#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu) and David Woolford (woolford@bcm.edu)
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
from emimagerotor import EMImageRotor
from emimagemxrotor import EMImageMXRotor
from emrotor import EMRotor
from emimage import EMImage,get_app
from math import *
from time import *
import os
import sys
import signal
from copy import *
from OpenGL import contextdata

# import SPARX definitions
import global_def
from global_def import *
# end import SPARX definitions


from emglplot import *

from time import time,sleep

from sys import getrefcount

from emglobjects import EMOpenGLFlagsAndTools

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
	parser.add_option("--write_coord_file",action="store_true",help="Write data box files",default=False)
	parser.add_option("--write_box_images",action="store_true",help="Write data box files",default=False)
	parser.add_option("--force","-f",action="store_true",help="Force overwrites old files",default=False)
	parser.add_option("--overlap",type="int",help="(auto:grid) number of pixels of overlap between boxes. May be negative.")
	#parser.add_option("--nretest",type="int",help="(auto:ref) Number of reference images (starting with the first) to use in the final test for particle quality.",default=-1)
	#parser.add_option("--retestlist",type="string",help="(auto:ref) Comma separated list of image numbers for retest cycle",default="")
	#parser.add_option("--farfocus",type="string",help="filename or 'next', name of an aligned far from focus image for preliminary boxing",default=None)
	parser.add_option("--parallel",type="int",help="specify more than one processor",default=1)
	parser.add_option("--merge_boxes_to_db",action="store_true",help="A special argument, if true all input arguments are considered to be box files and they are merged into the project database as manually selected particles",default=False)
	parser.add_option("--subsample_method",help="The method used to subsample images prior to generation of the correlation image. Available methods are standard,careful",default="standard")	
	
	(options, args) = parser.parse_args()
	if len(args)<1 : parser.error("Input image required")
	
	filenames = []
	for i in range(0,len(args)):
		filenames.append(args[i])

	logid=E2init(sys.argv)
	
	if options.merge_boxes_to_db == True:
		#The user wants to add some boxes to the database
		merge_boxes_as_manual_to_db(filenames)
		sys.exit(1)
		
	if not options.gui and not options.auto:
		parser.error("Atleast one of the --gui or --auto arguments are required.")
		exit(1)
	
	# we need to know how big to make the boxes. If nothing is specified, but
	# reference particles are, then we use the reference particle size
	#if options.boxsize<5 :
		#if not options.boxsize in good_box_sizes:
			#print "Note: EMAN2 processing would be more efficient with a box_size of %d"%good_box_size(options.boxsize)
	
	boxes=[]
	if len(options.auto)>0 :
		print "Autobox mode ",options.auto[0]
	
		if "db" in options.auto:
			print "auto data base boxing"
		
			autobox_multi(filenames,options)
			#if len(filenames) == 1:
				#autobox_single(filenames[0],options)
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


def merge_boxes_as_manual_to_db(filenames):
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
			manualboxes = get_idd_key_entry(filename,"manual_boxes")
		except:
			manualboxes = []
	
		if manualboxes == None: manualboxes = []
		manualboxes.extend(boxes)
		set_idd_key_entry(filename,"manual_boxes",manualboxes)


def autobox_multi(image_names,options):
	project_db = EMProjectDB()
	for image_name in image_names:
		print "autoboxing",image_name
		
		try:
			data = project_db[get_idd_key(image_name)]
			trim_autoboxer = project_db[data["autoboxer_unique_id"]]["autoboxer"]
			autoboxer = SwarmAutoBoxer(None)
			autoboxer.become(trim_autoboxer)
			print 'using cached autoboxer db'
		except:
			try:
				print "using most recent autoboxer"
				trim_autoboxer = project_db["current_autoboxer"]
				autoboxer = SwarmAutoBoxer(None)
				autoboxer.become(trim_autoboxer)
			except:
				print "Error - there seems to be no autoboxing information in the database - autobox interactively first - bailing"
				continue
		
		boxable = Boxable(image_name,None,autoboxer)
		
		if boxable.is_excluded():
			print "Image",image_name,"is excluded and being ignored"
			continue
		
		autoboxer.set_mode_explicit(SwarmAutoBoxer.COMMANDLINE)
		# Tell the boxer to delete non refs - FIXME - the uniform appraoch needs to occur - see SwarmAutoBoxer.auto_box
		autoboxer.auto_box(boxable,False)
		if options.write_coord_file:
			boxable.write_coord_file(-1,options.force)
		if options.write_box_images:
			boxable.write_box_images(-1,options.force)
	
	
	project_db.close()

def autobox_single(image_name,options):
	
	project_db = EMProjectDB()
	try:
		data = project_db[get_idd_key(image_name)]
		trim_autoboxer = project_db[data["autoboxer_unique_id"]]["autoboxer"]
		autoboxer = SwarmAutoBoxer(None)
		autoboxer.become(trim_autoboxer)
		print 'using cached autoboxer db'
	except:
		try:
			trim_autoboxer = project_db["current_autoboxer"]
			autoboxer = SwarmAutoBoxer(None)
			autoboxer.become(trim_autoboxer)
		except:
			print "Error - there seems to be no autoboxing information in the database - autobox interactively first - bailing"
			project_db.close()
			return 0
	
	boxable = Boxable(image_name,None,autoboxer)
	if boxable.is_excluded():
		print "Image",image_name,"is excluded and being ignored"
		return
	
	autoboxer.set_mode_explicit(SwarmAutoBoxer.COMMANDLINE)
	# Tell the boxer to delete non refs - FIXME - the uniform appraoch needs to occur - see SwarmAutoBoxer.auto_box
	autoboxer.auto_box(boxable,False)
	if options.write_coord_file:
		print "writing box coordinates"
		boxable.write_coord_file(-1,options.force)
	if options.write_box_images:
		print "writing boxed images"
		boxable.write_box_images(-1,options.force)
	
	project_db.close()
	return 1
	
class AutoDBBoxer(QtCore.QObject):
	'''
	A class for managing the process of spawning many instances of e2boxer singlefile.mrc --auto=db, using parallelism
	If one CPU is specified then it still works. Basically the approach is to spawn the number of the processors, then once
	the process is finished the signal is intercepted and a new process is executed etc.
	'''
	def __init__(self,image_names,nproc,options,force=False):
		QtCore.QObject.__init__(self)
		self.nproc = nproc
		self.image_names = image_names
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
		
		for i in range(0,len(self.image_names)):
			self.processes.append(QtCore.QProcess(self))
			
	def print_cp_status(self):
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
		if self.currentidx >= len(self.image_names) :
			return
			
		#process = self.processes[self.currentidx]
			
		program = QtCore.QString("e2boxer.py")
		args = QtCore.QStringList()
		args.append(self.image_names[self.currentidx])
		args.append("--auto=db")
		if self.options.write_coord_file != False:
			args.append("--write_coord_file")
		if self.options.write_box_images != False:
			args.append("--write_box_images")
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
		if self.jobsdone == len(self.image_names):
			self.app.quit()
		self.spawn_process()
		self.currentidx += 1
		
		
	#def process_start(self):
		#print "process started"
		

class GUIboxMouseEventsObject:
	'''
	A base class for objects that handle mouse events in the GUIbox
	
	Inheriting objects are concerned with supplying their own definitions of 
	the functions mouse_up, mouse_down, mouse_drag, mouse_move and mouse_wheel. 
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
		
	def get_2d_gui_image(self):
		'''
		Ask the mediator to for the EMImage2D object
		'''
		return self.mediator.get_2d_gui_image()
	
	def get_gui_ctl(self):
		'''
		Ask the mediator to for the main GUI controller
		'''
		return self.mediator.get_gui_ctl()
	
	def get_mx_gui_image(self):
		'''
		Ask the mediator to for the EMImageMX object
		'''
		return self.mediator.get_mx_gui_image()
	
	def mouse_up(self,event):
		'''
		Inheriting classes to potentially define this function
		'''
		pass

	def mouse_down(self,event):
		'''
		Inheriting classes to potentially define this function
		'''
		pass
		
	def mouse_drag(self,event):
		'''
		Inheriting classes to potentially define this function
		'''
		pass
	
	def mouse_move(self,event):
		'''
		Inheriting classes to potentially define this function
		'''
		pass

	def mouse_wheel(self,event):
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
		
	def set_mode(self,mode):
		self.erasemode = mode
		
	def set_erase_radius(self,radius):
		self.eraseradius = radius
	
	def mouse_move(self,event):
		m = self.get_2d_gui_image().scr2img((event.x(),event.y()))
		self.get_2d_gui_image().addShape("eraser",EMShape(["circle",.1,.1,.1,m[0],m[1],self.eraseradius,3]))
		self.mediator.update_image_display()
		
	def mouse_wheel(self,event):
		if event.modifiers()&Qt.ShiftModifier:
			self.get_gui_ctl().adjust_erase_rad(event.delta())
			m= self.get_2d_gui_image().scr2img((event.x(),event.y()))
			self.get_2d_gui_image().addShape("eraser",EMShape(["circle",.1,.1,.1,m[0],m[1],self.eraseradius,3]))
			self.mediator.update_image_display()
	
	def mouse_down(self,event) :
		m=self.get_2d_gui_image().scr2img((event.x(),event.y()))
		#self.boxable.add_exclusion_area("circle",m[0],m[1],self.eraseradius)
		self.get_2d_gui_image().addShape("eraser",EMShape(["circle",.9,.9,.9,m[0],m[1],self.eraseradius,3]))
		self.mediator.exclusion_area_added("circle",m[0],m[1],self.eraseradius,self.erasemode)	

	def mouse_drag(self,event) :
		m=self.get_2d_gui_image().scr2img((event.x(),event.y()))
		self.get_2d_gui_image().addShape("eraser",EMShape(["circle",.9,.9,.9,m[0],m[1],self.eraseradius,3]))
		self.mediator.exclusion_area_added("circle",m[0],m[1],self.eraseradius,self.erasemode)
		# exclusion_area_added does the OpenGL update calls, so there is no need to do so here
		
	def mouse_up(self,event) :
		# we have finished erasing
		
		# make the eraser shape non visible
		self.get_2d_gui_image().addShape("eraser",EMShape(["circle",0,0,0,0,0,0,0.1]))
		self.mediator.erasing_done(self.erasemode)
	
class GUIboxParticleManipEvents(GUIboxMouseEventsObject):
	'''
	A class that knows how to add, move and remove reference and non reference boxes 
	'''
	def __init__(self,mediator):
		GUIboxMouseEventsObject.__init__(self,mediator)
		self.mode =  GUIbox.REFERENCE_ADDING
		self.moving = None
		self.dynapix = False

	def set_mode(self,mode):
		if mode not in [GUIbox.REFERENCE_ADDING,GUIbox.MANUALLY_ADDING]:
			print 'error, that is  an illegal mode'
			return
		
		self.mode = mode
		
	def mouse_down(self,event) :
		m = self.get_2d_gui_image().scr2img((event.x(),event.y()))
		box_num = self.mediator.detect_box_collision(m)
		if box_num == -1:
			if not self.mediator.within_main_image_bounds(m):	return
			#if we make it here, that means the user has clicked on an area that is not in any box
			
			if event.modifiers()&Qt.ShiftModifier : return # the user tried to delete nothing
			
			# If we get here, we need to add a new reference
			box_size = self.mediator.get_box_size()
			
			box = Box(m[0]-box_size/2,m[1]-box_size/2,box_size,box_size,True)
			box.set_image_name(self.mediator.get_current_image_name())

			box.changed = True # this is so image2D nows to repaint the shape
			
			if self.mode == GUIbox.REFERENCE_ADDING:
				box.isref = True
				box.ismanual = False
				
			elif self.mode == GUIbox.MANUALLY_ADDING:
				box.isref = False
				box.ismanual = True
			else:
				print 'error, unknown error in mouse_down, boxing mode'
			
			x0=box.xcorner+box.xsize/2-1
			y0=box.ycorner+box.ysize/2-1
			self.get_2d_gui_image().addShape("cen",EMShape([self.mediator.get_shape_string(),.9,.9,.4,x0,y0,x0+2,y0+2,1.0]))
			
			self.mediator.add_box(box)
			#self.mediator.mouse_click_update_ppc()
		
		elif event.modifiers()&Qt.ShiftModifier :
			# remove the box
			self.mediator.remove_box(box_num)
			#self.mediator.mouse_click_update_ppc()
			
		else:
			# if we make it here than the we're moving a box
			box = self.mediator.get_box(box_num)
			self.moving=[box,m,box_num]
			self.get_2d_gui_image().setActive(box_num,.9,.9,.4)
				
			x0=box.xcorner+box.xsize/2-1
			y0=box.ycorner+box.ysize/2-1
			self.get_2d_gui_image().addShape("cen",EMShape([self.mediator.get_shape_string(),.9,.9,.4,x0,y0,x0+2,y0+2,1.0]))
			object = self.get_mx_gui_image().get_core_object()
			if object.is_visible(box_num) or True : self.get_mx_gui_image().get_core_object().set_selected([box_num],True)
			self.mediator.update_all_image_displays()
			

	def mouse_drag(self,event) :
		
		m=self.get_2d_gui_image().scr2img((event.x(),event.y()))
		
		if event.modifiers()&Qt.ShiftModifier:
			box_num = self.mediator.detect_box_collision(m)
			if ( box_num != -1):
				self.mediator.remove_box(box_num)
				#self.mediator.mouse_click_update_ppc()
			
		elif self.moving != None:
			# self.moving[0] is the box, self.moving[1] are the mouse coordinates
			box = self.moving[0]
			# the old m in in self.moving[2]
			oldm = self.moving[1]
			
			self.mediator.move_box(self.moving[2],m[0]-oldm[0],m[1]-oldm[1])
			self.moving[1] = m
	
	def mouse_up(self,event) :
		if self.moving != None:
			box = self.moving[0]
			if box.isref:
				self.mediator.reference_moved(box)
			self.moving=None
		
		m = self.get_2d_gui_image().scr2img((event.x(),event.y()))
		box_num = self.mediator.detect_box_collision(m)
		if box_num != -1 and not event.modifiers()&Qt.ShiftModifier:
			object = self.get_mx_gui_image().get_core_object()
			if not object.is_visible(box_num) : object.scroll_to(box_num,True)
			self.get_mx_gui_image().get_core_object().set_selected([box_num],True)
			self.mediator.update_all_image_displays()

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

	def get_2d_gui_image(self):
		'''
		Return the parent's EMImage2D object
		'''
		return self.parent.get_2d_gui_image()
	
	def get_gui_ctl(self):
		'''
		Return the parent's Controller widgit object
		'''
		return self.parent.get_gui_ctl()
	
	def get_mx_gui_image(self):
		'''
		Return the parent's EMImageMX object
		'''
		return self.parent.get_mx_gui_image()
	
	def update_image_display(self):
		'''
		Send an event to the parent that the EMImage2D should update its display
		'''
		self.parent.update_image_display()
		
	def update_all_image_displays(self):
		'''
		Send an event to the parent that the EMImage2D and EMImageMX objects should
		update their displays
		'''
		self.parent.update_all_image_displays()
		
	def exclusion_area_added(self,typeofexclusion,x,y,radius,mode):
		'''
		Send an event to the parent that an exclusion area was added.
		The parameters define the type of exclusion area. In future the exclusion area
		should probably just be its own class.
		'''
		self.parent.exclusion_area_added(typeofexclusion,x,y,radius,mode)

	def erasing_done(self,erase_mode):
		'''
		Send an event to the parent letting it know that the user has stopped adding
		erased area
		'''
		self.parent.erasing_done(erase_mode)
		
	def detect_box_collision(self,coords):
		'''
		Ask the parent to detect a collision between a resident box and the given coordinates.
		This is in terms of the EMImage2D 
		'''
		return self.parent.detect_box_collision(coords)
	
	def get_current_image_name(self):
		'''
		Ask the parent for the current name of the image in the large image view...
		'''
		return self.parent.get_current_image_name()

	def get_box_size(self):
		'''
		Ask the parent for the current project box_size
		'''
		return self.parent.get_box_size()
	
	def add_box(self,box):
		'''
		Tell the parent to store a box
		'''
		self.parent.add_box(box)
	
	def box_display_update(self):
		'''
		Tell the parent the a general box display update needs to occur
		'''
		self.parent.box_display_update()
		
	#def mouse_click_update_ppc(self):
		#'''
		#Tell the parent that a mouse click occured and that the PPC metric should be updated
		#'''
		#self.parent.mouse_click_update_ppc()
	
	def remove_box(self,box_num):
		'''
		Tell the parent to remove a box, as given by the box_num
		'''
		self.parent.remove_box(box_num)
		
	def get_box(self,box_num):
		'''
		Ask the parent for the box, as given by the box_num
		'''
		return self.parent.get_box(box_num)
	
	def move_box(self,box_num,dx,dy):
		'''
		Tell the parent to handle the movement of a box
		'''
		self.parent.move_box(box_num,dx,dy)
	
	def reference_moved(self,box):
		'''
		Tell the parent that a reference was moved - this could trigger automatic boxing
		'''
		self.parent.reference_moved(box)
		
	def within_main_image_bounds(self,coords):
		'''
		Ask the parent to determine if the coords are within the currently display image boundaries
		'''
		return self.parent.within_main_image_bounds(coords)
	
	def get_shape_string(self):
		'''
		Gets the shape string currently used for creating shapes for the 2D image
		'''
		return self.parent.get_shape_string()
	
	
class GUIbox:
	'''
	Cleaning in progress
	'''
	REFERENCE_ADDING = 0
	ERASING = 1
	MANUALLY_ADDING = 2
	FANCY_MODE = 'fancy'
	PLAIN_MODE = 'plain'
	def __init__(self,image_names,boxes,box_size=-1):
		"""Implements the 'boxer' GUI."""
		
		self.dynapixp=get_app() # get the app
		
		# initialize important autoboxer related variables
		self.dynapix = False
		self.image_names = image_names
		self.current_image_idx = 0
		self.box_size = box_size
		# set self.autoboxer 
		self.set_autoboxer(self.image_names[self.current_image_idx])
		
		self.eraseradius = 2*self.box_size # this happens after the autoboxer has been loaded, because the boxsize can change
		self.erasemode = None #stores the erase mode
		self.shape_string = "rectpoint" # the shape of the picked particles
		self.in_display_limbo = False	# a flag I am using to solve a problem
		
		# A boxable is just a class that manages boxes in terms of images
		self.boxable = Boxable(self.image_names[0],self,self.autoboxer)
		self.boxable.add_non_refs(boxes)
		self.boxable.set_box_size(self.box_size)
		
		self.initialize_mouse_event_handlers() # initialize the mouse event handlers

		self.ptcl=[] # list of actual boxed out EMImages. This may be redundant I am working on a better solution.
		self.moving_box_data = None # a vector storing [mouse x, mouse y, box idx]
		self.moving=None # Used during a user box drag. May be redudant could potentially just use self.moving_box_data. FIXME
		
		self.init_guiim() # initialise the 2D image display
		self.__init_guimx() # intialize the matrix display
		self.__init_guimx_thumbs() # initialize the thumbnail diplsy
		self.ab_sel_mediator = AutoBoxerSelectionsMediator(self)
		self.__init_guictl()
		if self.fancy_mode == GUIbox.FANCY_MODE:
			self.__init_ctl_rotor()
		
		self.autoboxer.auto_box(self.boxable,False) # Do the automatic autoboxing - this makes the user see results immediately
		self.box_display_update() # update displays to show boxes etc
		
	def __init_ctl_rotor(self):
		self.ctl_rotor = EMRotor()
		self.ctl_rotor.get_core_object().add_qt_widget(self.guictl)
		self.guictlrotor = EMParentWin(self.ctl_rotor)
		self.guictlrotor.setWindowTitle("e2boxer Controllers")
		self.guictlrotor.show()

	def __init_guictl(self):
		self.guictl=GUIboxPanel(self,self.ab_sel_mediator)
		self.guictl.set_image_quality(self.boxable.get_quality())
		self.guictl.setWindowTitle("e2boxer Controller")
		self.guictl.set_dynapix(self.dynapix)
		self.guictl.show()
		if self.fancy_mode == GUIbox.FANCY_MODE: self.guictl.hide()
	
	def __init_guimx_thumbs(self):
		self.itshrink = -1 # image thumb shrink. Default value of -1 means it has to be calculated when it's first needed
		self.imagethumbs = None # image thumbs - will be a list of tiny images
		self.guimxitp = None
		self.guimxit = None
		self.__gen_image_thumbnails_widget()
		
		if self.guimxitp != None:
			self.guimxitp.show()
			if isinstance(self.guimxit,EMImageRotor):
				self.guimxitp.resize(*self.guimxit.get_optimal_size())
			self.guimxit.connect(self.guimxit,QtCore.SIGNAL("mousedown"),self.image_selected)
			if isinstance(self.guimxit,EMImageRotor):
				self.guimxit.set_frozen(self.boxable.is_frozen(),self.current_image_idx)
	
	def __init_guimx(self):
		self.guimxp= None # widget for displaying matrix of smaller imagespaugay
		glflags = EMOpenGLFlagsAndTools()
		if True or not glflags.npt_textures_unsupported():
			self.guimx=EMImageMXRotor()		# widget for displaying image thumbs
			self.guimx.get_core_object().disable_mx_zoom()
			self.guimx.get_core_object().allow_camera_rotations(False)
			self.guimx.get_core_object().disable_mx_translate()
			#self.fancy_mode = GUIbox.FANCY_MODE
			self.fancy_mode = GUIbox.PLAIN_MODE # for now just make sure the fancy stuff isn't used
		else:
			self.guimx=EMImageMX()	
			self.fancy_mode = GUIbox.PLAIN_MODE
		
		self.guimx.set_mmode("app")
		
		if self.fancy_mode == GUIbox.FANCY_MODE:
			 self.guimx.connect(self.guiim,QtCore.SIGNAL("inspector_shown"),self.guiim_inspector_requested)
		#self.guimx.connect(self.guimx,QtCore.SIGNAL("removeshape"),self.removeshape)
		self.guimx.connect(self.guimx,QtCore.SIGNAL("mousedown"),self.box_selected)
		self.guimx.connect(self.guimx,QtCore.SIGNAL("mousedrag"),self.box_moved)
		self.guimx.connect(self.guimx,QtCore.SIGNAL("mouseup"),self.box_released)
		self.guimx.connect(self.guimx,QtCore.SIGNAL("boxdeleted"),self.box_image_deleted)
		if self.fancy_mode == GUIbox.FANCY_MODE:
			self.guimx.connect(self.guimx,QtCore.SIGNAL("inspector_shown"),self.guimx_inspector_requested)
		
	def init_guiim(self, image=None, imagename=None):
		if image is None:
			imagename = self.image_names[self.current_image_idx]
			image=BigImageCache.get_image_directly(imagename)
		
		self.guiimp=EMImage(image)		# widget for displaying large image
		self.guiimp.setWindowTitle(imagename)
		self.guiim=self.guiimp.child
		self.__update_guiim_states()
		
		self.guiim.set_mmode(0)
		self.guiimp.show()
		
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mousedown"),self.mouse_down)
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mousedrag"),self.mouse_drag)
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mouseup")  ,self.mouse_up  )
		self.guiim.connect(self.guiim,QtCore.SIGNAL("keypress"),self.keypress)
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mousewheel"),self.mouse_wheel)
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mousemove"),self.mouse_move)
	
	def __update_guiim_states(self):
		guiim_core = self.guiim.get_core_object()
		
		guiim_core.set_other_data(self.boxable.get_exclusion_image(False),self.autoboxer.get_subsample_rate(),True)
		
		guiim_core.set_frozen(self.boxable.is_frozen())
		guiim_core.set_excluded(self.boxable.is_excluded())
		guiim_core.set_file_name(self.image_names[self.current_image_idx])
	
	def set_autoboxer(self,imagename):
		'''
		This function sets the self.autoboxer variable. There are 3 cases:
		
		1. The image already has an autoboxer in the database which becomes self.autoboxer
		2. There is a current autoboxer in the database in which case this becomes self.autoboxer
		3. Neither 1 nor 2 - then an empty autoboxer is created
		
		'''
		try:
			project_db = EMProjectDB()
			data = project_db[get_idd_key(imagename)]
			autoboxer_id = data["autoboxer_unique_id"]
			trim_autoboxer = project_db[autoboxer_id]["autoboxer"]
			self.autoboxer_name = autoboxer_id
			self.autoboxer = SwarmAutoBoxer(self)
			self.autoboxer.become(trim_autoboxer)
			self.dynapix = self.autoboxer.dynapix_on()
			if self.box_size==-1: self.box_size = self.autoboxer.get_box_size()
		except:
			try:
				trim_autoboxer = project_db["current_autoboxer"]
				self.autoboxer = SwarmAutoBoxer(self)
				self.autoboxer.become(trim_autoboxer)
				self.autoboxer_name = self.autoboxer.get_unique_stamp()
				self.dynapix = self.autoboxer.dynapix_on()
				if self.box_size==-1: self.box_size = self.autoboxer.get_box_size()
			except:	
				if self.box_size == -1: self.box_size = 128
				self.autoboxer = SwarmAutoBoxer(self)
				self.autoboxer.set_box_size_explicit(self.box_size)
				self.autoboxer.set_interactive_mode(self.dynapix)

	def guiim_inspector_requested(self,event):
		self.ctl_rotor.get_core_object().add_qt_widget(self.guiim.get_core_object().get_inspector())
	
	def guimx_inspector_requested(self,event):
		self.ctl_rotor.get_core_object().add_qt_widget(self.guimx.get_core_object().get_inspector())
	
	def initialize_mouse_event_handlers(self):
		self.eventsmediator = GUIboxEventsMediator(self)
		self.mouse_handlers = {}
		self.mouse_handlers["boxing"] = GUIboxParticleManipEvents(self.eventsmediator)
		self.mouse_handlers["erasing"] = GUIboxMouseEraseEvents(self.eventsmediator,self.eraseradius)
		self.mousehandler = self.mouse_handlers["boxing"]
	
	def get_shape_string(self):
		return self.shape_string
	
	def has_thumbnails(self):
		return self.guimxitp != None
	
	def view_boxes_clicked(self,bool):
		if bool:
			self.guimxp.show()
		else:
			self.guimxp.hide()
			
	def view_image_clicked(self,bool):
		print self.guiimp
		if bool:
			self.guiimp.show()
		else:
			self.guiimp.hide()
			
	def view_thumbs_clicked(self,bool):
		if self.guimxitp == None: return
		if bool:
			self.guimxitp.show()
		else:
			self.guimxitp.hide()
			
	def within_main_image_bounds(self,coords):
		x = coords[0]
		y = coords[1]
		if x < 0 or y < 0:
			return False
		
		main_image = self.get_current_image()
		
		if x >= main_image.get_xsize() or y >= main_image.get_ysize():
			return False
		
		return True
		
	def change_current_autoboxer(self, autoboxer_id,autobox=True):
		#print "change current autoboxer"
		self.autoboxer.write_to_db()
		project_db = EMProjectDB()
		if project_db[autoboxer_id] == None:
			project_db[autoboxer_id] = {}
		trim_autoboxer = project_db[autoboxer_id]["autoboxer"]
		#print "changing autoboxer to autoboxer_",timestamp,"and its stamp in the db is",trim_autoboxer.get_creation_ts()
		self.autoboxer = SwarmAutoBoxer(self)
		self.autoboxer.become(trim_autoboxer)
		self.autoboxer.write_image_specific_references_to_db(self.boxable.get_image_name())
		self.dynapix = self.autoboxer.dynapix_on()
		self.guictl.set_dynapix(self.dynapix)
		
		self.boxable.set_autoboxer(self.autoboxer)
		
		if self.box_size != self.autoboxer.get_box_size():
			self.update_box_size(self.autoboxer.get_box_size())
			#print "box display update"
		
		if autobox:
			self.ptcl = []
			self.guiim.delShapes()
			self.boxable.clear_and_reload_images()
			self.in_display_limbo = True
			#self.autoboxer.regressiveflag = True
			self.autoboxer.auto_box(self.boxable)
			self.in_display_limbo = False
		
		self.box_display_update()	
		
	def update_box_size(self,box_size,mode=0):
		if box_size != self.box_size:
			if mode == 0:
				self.box_size = box_size
				self.boxable.change_box_size(self.box_size)
				self.guictl.set_box_size(self.box_size)
				self.boxable.get_exclusion_image(True)
				self.boxable.reload_boxes() # this may be inefficient
				#print "clearing displays"
				guiim_core = self.guiim.get_core_object()
				guiim_core.set_other_data(self.boxable.get_exclusion_image(False),self.autoboxer.get_subsample_rate(),True)
				self.clear_displays()
			elif mode == 1:
				self.box_size = box_size
				self.autoboxer.set_box_size(self.box_size,self.image_names)
				self.boxable.reload_boxes() # this may be inefficient
				#print "clearing displays"
				guiim_core = self.guiim.get_core_object()
				guiim_core.set_other_data(self.boxable.get_exclusion_image(False),self.autoboxer.get_subsample_rate(),True)
				self.clear_displays()	
			else:
				print "error, unknown mode in update_box_size"
		
	def get_2d_gui_image(self):
		return self.guiim
	
	def get_gui_ctl(self):
		return self.guictl
	
	def get_mx_gui_image(self):
		return self.guimx
	
	def get_boxable(self):
		return self.boxable
		
	def exclusion_area_added(self,typeofexclusion,x,y,radius,mode):
		self.boxable.add_exclusion_area(typeofexclusion,x,y,radius,mode)
		guiim_core = self.guiim.get_core_object()
		guiim_core.set_other_data(self.boxable.get_exclusion_image(False),self.autoboxer.get_subsample_rate(),True)
		self.update_image_display()
		
	def erasing_done(self,erase_mode):
		'''
		Call this function after erasing has occured to remove all boxes in the
		erased regions
		'''
		
		if erase_mode == Boxable.ERASE:
			# tell the boxable to remove boxes (and return their number)
			[lostboxes,refboxes] = self.boxable.update_excluded_boxes()
			# after this, make sure the display is correct.
			
			# If any of the boxes are references the autoboxer needs to know about it...
			# this could potentially trigger auto boxer 
			if len(lostboxes) != 0:
				self.delete_display_boxes(lostboxes)
				#self.mouseclicks += 1
				#self.update_ppc()
				self.update_all_image_displays()
				
			else: self.update_image_display()
			
			if len(refboxes) != 0: 
				val = self.autoboxer.remove_reference(refboxes)
				if val == 2:
					self.boxable.clear_and_cache(True)
					self.clear_displays()
					return # avoid unecessary display update below
			
			self.box_display_update()
			
		elif erase_mode == Boxable.UNERASE:
			[added_boxes,added_refboxes] = self.boxable.update_included_boxes()
			
			# len(added_boxes) is always <= len(added_refboxes)
			if len(added_boxes) != 0:
				self.box_display_update()
				
				if len(added_refboxes) != 0:
					val = self.autoboxer.add_reference(add_refboxes)

				
	def detect_box_collision(self,coords):
		'''
		Detects a collision of the coordinates with any of the boxes in the current image
		stored in guiim. coords need to be in the image coordinates
		'''
		return self.collision_detect(coords,self.get_boxes())
	
	def get_current_autoboxer_ts(self):
		return self.autoboxer.get_creation_ts()
	
	def get_current_image_name(self):
		return self.image_names[self.current_image_idx]
	
	def get_current_image(self):
		return BigImageCache.get_image_directly(self.get_current_image_name())
	
	def get_image_names(self):
		return self.image_names
	
	def get_box_size(self):
		return self.box_size
	
	def add_box(self,box):
		if not self.boxable.is_interactive():
			return
		
		self.boxable.add_box(box)
		
		# Should this be here?
		box_num = len(self.get_boxes())
		#self.guimx.get_core_object().set_selected([box_num])
		
		# autoboxer will autobox depending on the state of its mode
		if box.isref : self.autoboxer.add_reference(box)
		
		self.box_display_update()
		
	def get_box(self,box_num):
		boxes = self.get_boxes()
		return boxes[box_num]
		
	def remove_box(self,box_num):
		if not self.boxable.is_interactive(): return
		
		box = self.delete_box(box_num)
		if not (box.isref or box.ismanual):
			self.boxable.add_exclusion_particle(box)
		guiim_core = self.guiim.get_core_object()
		guiim_core.set_other_data(self.boxable.get_exclusion_image(False),self.autoboxer.get_subsample_rate(),True)
		self.box_display_update()
		
	def reference_moved(self,box):
		if (self.autoboxer.reference_moved(box)):
			self.box_display_update()
		
	#def mouse_click_update_ppc(self):
		#self.mouseclicks += 1
		#self.update_ppc()
	
	def set_ptcl_mx_data(self,data=None):
		'''
		Call this to set the Ptcl Mx data 
		'''
		if data != None:
			self.guimx.setData(data)
			if len(data) != 0:
				if self.guimxp == None:
					self.guimxp = EMParentWin(self.guimx)
					self.guimxp.setWindowTitle("Particles")
					self.guimxp.show()
					#self.guimx.get_core_object().set_selected([0])
	
	def clear_displays(self):
		self.ptcl = []
		self.guiim.delShapes()
		self.guimx.setData([])
		self.box_display_update() # - the user may still have some manual boxes...
	
	def image_selected(self,event,lc):
		app = QtGui.QApplication.instance()
		app.setOverrideCursor(Qt.BusyCursor)
		#try:
		im=lc[0]
		#try:
		debug = False
		if im != self.current_image_idx:
				#print 'changing images'		
			
			if debug: tt = time()
			image=BigImageCache.get_image_directly(self.image_names[im])
			
			
			if debug: tt = time()
			self.guiimp.setWindowTitle(self.image_names[im])
			if debug: print "it took", time() - tt, "to read the prior stuff A1 "
			if debug: tt = time()
			self.guiim.setData(image)
			if debug: print "it took", time() - tt, "to read the prior stuff A2 "
			if debug: tt = time()
			self.boxable.cache_exc_to_db()
			if debug: print "it took", time() - tt, "to read the prior stuff A3 "
			if debug: tt = time()
			self.boxable = Boxable(self.image_names[im],self,self.autoboxer)
			if debug: print "it took", time() - tt, "to read the prior stuff A4 "
			if debug: tt = time()
			self.ptcl = []
			self.guiim.delShapes()
			self.guiim.get_core_object().force_display_update()
			self.in_display_limbo = True
			if debug: print "it took", time() - tt, "to read the prior stuff B "
			if debug: tt = time()
			project_db = EMProjectDB()
			data = project_db[get_idd_key(self.image_names[im])]
			ab_failure = self.autoboxer
			if debug: print "it took", time() - tt, "to read the prior stuff C "
			
			if debug: tt = time()
			if data != None:
				try:
					autoboxer_id = data["autoboxer_unique_id"]
					trim_autoboxer = project_db[autoboxer_id]["autoboxer"]
					self.autoboxer_name = autoboxer_id
					self.autoboxer = SwarmAutoBoxer(self)
					self.autoboxer.become(trim_autoboxer)
					self.dynapix = self.autoboxer.dynapix_on()
					self.guictl.set_dynapix(self.dynapix)
				except:
					try:
						trim_autoboxer = project_db["current_autoboxer"]
						self.autoboxer = SwarmAutoBoxer(self)
						self.autoboxer.become(trim_autoboxer)
						self.autoboxer_name = self.autoboxer.get_unique_stamp()
						self.dynapix = self.autoboxer.dynapix_on()
						self.guictl.set_dynapix(self.dynapix)
					except:
						self.autoboxer = ab_failure
			if debug: print "it took", time() - tt, "to the autobox database stuff"
			
			if self.dynapix:
				#self.autoboxer.regressiveflag = True
				if debug: tt = time()
				self.autoboxer.auto_box(self.boxable,False)
				if debug: print "it took", time() - tt, "to autobox"
				
			self.boxable.set_autoboxer(self.autoboxer)
			
			self.autoboxer_db_changed()
			
		
			if self.box_size != self.autoboxer.get_box_size():
				self.update_box_size(self.autoboxer.get_box_size())

			self.in_display_limbo = False
			
			for box in self.boxable.boxes: box.changed = True
			
			self.current_image_idx = im
			
			self.__update_guiim_states()
			if self.guimxitp != None and isinstance(self.guimxit,EMImageRotor):
				self.guimxit.set_frozen(self.boxable.is_frozen(),self.current_image_idx)
			self.guictl.set_image_quality(self.boxable.get_quality())
			
			if debug: tt = time()
			self.box_display_update()
			if debug: print "it took", time() - tt, "to do the box display update"
		#except: pass
			
		app.setOverrideCursor(Qt.ArrowCursor)
			

	def __gen_image_thumbnails_widget(self):
		'''
		Generates image thumbnails for a single image name
		if there is only one image in the image file on disk
		this function returns 0 and does nothing.
		Else thumbnails are generated, and image matrix is initialized (but not shown),
		and 1 is returned
		'''
		# warnilg here
		try: nim = len(self.image_names)
		except: 
			# warning - bad hacking going on
			print "the image name ", image_name, "probably doesn't exist"
			raise Exception
	
		if (nim == 1): return 0
		
		app = QtGui.QApplication.instance()
		app.setOverrideCursor(Qt.BusyCursor)
	
		try:
			n = self.get_image_thumb_shrink()
			self.imagethumbs = []
			
			a = time()
			self.imagethumbs = []
			
			for i in range(0,nim):
				self.imagethumbs.append(None)
			
			for i in range(nim-1,-1,-1):
				
				thumb = self.get_image_thumb(i)
				#print "got thumb",i
				self.imagethumbs[i] = thumb
			glflags = EMOpenGLFlagsAndTools()
			if True or not glflags.npt_textures_unsupported():
				self.guimxit=EMImageRotor()		# widget for displaying image thumbs
			else:
				self.guimxit=EMImageMX()
				
			self.guimxit.setData(self.imagethumbs)
			
			try:
				for i in range(0,nim):
					print i
					frozen = get_idd_key_entry(self.image_names[i],"frozen_state")
					if frozen != None:
						self.guimxit.set_frozen(frozen,i)
			except: print "setting frozen failed"
					
			self.guimxitp = EMParentWin(self.guimxit)
			self.guimxitp.setWindowTitle("Image Thumbs")
			
			self.guimxit.set_mmode("app")
			app = QtGui.QApplication.instance()
			app.setOverrideCursor(Qt.BusyCursor)
		except: 
			app.setOverrideCursor(Qt.ArrowCursor)
			return 0
			pass
		
		app.setOverrideCursor(Qt.ArrowCursor)
		return 1
		
	def get_image_thumb(self,i):
		
		image = get_idd_key_entry(self.image_names[i],"e2boxer_image_thumb")
		if image == None:
			n = self.get_image_thumb_shrink()

			image=BigImageCache.get_image_directly(self.image_names[i])
			
			#while n > 1:
				#image = image.process("math.meanshrink",{"n":2})
				#n /= 2
			image = image.process("math.meanshrink",{"n":n})
			set_idd_key_entry(self.image_names[i],"e2boxer_image_thumb",image)
		return image
		
	def get_image_thumb_shrink(self):
		if self.itshrink == -1:
			image=BigImageCache.get_image_directly(self.image_names[self.current_image_idx])
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
		
	def box_moved(self,event,scale):
		try:
			dx = (self.moving_box_data[0] - event.x())/scale
			dy = (event.y() - self.moving_box_data[1])/scale
			self.move_box(self.moving_box_data[2],dx,dy)
			self.moving_box_data[0] = event.x()
			self.moving_box_data[1] = event.y()
			self.update_all_image_displays()
		except: pass
		
	def box_released(self,event,lc):
		boxes = self.get_boxes()
		#print boxes
		try: box =  boxes[self.moving_box_data[2]]
		except: return # The user clicked on black
		
		if box.isref :
			self.reference_moved(box)

		im=lc[0]
		self.moving_box_data = [event.x(),event.y(),im]
		self.guiim.setActive(im,.9,.9,.4)
		#self.guimx.get_core_object.set_selected(im)
		boxes = self.get_boxes()
		self.guiim.registerScrollMotion(boxes[im].xcorner+boxes[im].xsize/2,boxes[im].ycorner+boxes[im].ysize/2)
		
		self.moving_box_data = None
		
	def box_selected(self,event,lc):
		im=lc[0]
		self.moving_box_data = [event.x(),event.y(),im]
		self.guiim.setActive(im,.9,.9,.4)
		#object = self.guimx.get_core_object()
		#if not object.is_visible(lc[0]) : object.scroll_to(lc[0],True)
		#self.get_mx_gui_image().get_core_object().set_selected([lc[0]],True)
		boxes = self.get_boxes()
		#self.guiim.registerScrollMotion(boxes[im].xcorner+boxes[im].xsize/2,boxes[im].ycorner+boxes[im].ysize/2)
		#try:
			##self.guiim.scrollTo(boxes[im].xcorner+boxes[im].xsize/2,boxes[im].ycorner+boxes[im].ysize/2)
			#pass
			
		#except: print "box_selected() scrolling error"

	def get_boxes(self):
		return self.boxable.boxes
	
	def mouse_move(self,event):
		self.mousehandler.mouse_move(event)

	def mouse_wheel(self,event):
		self.mousehandler.mouse_wheel(event)
		
	def mouse_down(self,event) :
		self.mousehandler.mouse_down(event)
		
	def mouse_up(self,event) :
		self.mousehandler.mouse_up(event)
	
	def mouse_drag(self,event) :
		self.mousehandler.mouse_drag(event)
	
	#def update_ppc(self):
		#if self.mouseclicks > 0:
			#self.guictl.ppc_changed(len(self.get_boxes())/float(self.mouseclicks))
		#else:
			#self.guictl.ppc_changed(0)
	
	def collision_detect(self,m,boxes):
			
		for box_num,box in enumerate(boxes):
			if m[0]<box.xcorner or m[0]>(box.xcorner +box.xsize) or m[1]<box.ycorner or m[1]>(box.ycorner +box.ysize) :
				# no collision
				continue
			# if we make it here there has been a collision, the box already exists
			return box_num
		
		return -1

	def move_box(self,box_num,dx,dy):
		box = self.get_boxes()[box_num]
		if not self.boxable.is_interactive():
			return
		
		self.boxable.move_box(box,dx,dy,box_num)
			# we have to update the reference also
		self.ptcl[box_num] = box.get_box_image()
			
		x0=box.xcorner+box.xsize/2-1
		y0=box.ycorner+box.ysize/2-1
		self.guiim.addShape("cen",EMShape([self.shape_string,.9,.9,.4,x0,y0,x0+2,y0+2,1.0]))
		box.shape = EMShape([self.shape_string,box.r,box.g,box.b,box.xcorner,box.ycorner,box.xcorner+box.xsize,box.ycorner+box.ysize,2.0])
		self.guiim.addShape(box_num,box.shape)
		self.box_display_update_specific(box_num)

	def change_shapes(self,shape_string):
		if shape_string in ["rectpoint","rect","rcircle","rcirclepoint"]:
			self.shape_string = shape_string
			self.box_display_update(True)
		else:
			print "unknown shape string", shapestring
	
	def update_box_colors(self,classify):
		sh=self.guiim.getShapes()
		for i in classify.items():
			color = BoxingTools.get_color(i[1])
			
			sh[int(i[0])].shape[1] = color[0]
			sh[int(i[0])].shape[2] = color[1]
			sh[int(i[0])].shape[3] = color[2]
			sh[int(i[0])].changed=True
			
		self.box_display_update()

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
	
	def update_all_image_displays(self):
		self.update_image_display()
		self.update_mx_display()
		if self.guimxit != None: self.guimxit.updateGL()
	
		#context = contextdata.getContext(None)
		#print context
		
	def update_image_display(self):
		self.guiim.updateGL()
		if self.guimxit != None: self.guimxit.updateGL()
		
	def update_mx_display(self):
		self.guimx.updateGL()

	def delete_display_boxes(self,numbers):
		'''
		Warning - this won't work unless the numbers go from greatest to smallest- i.e. they are in reverse order
		Deletes shapes displayed by the 2D image viewer
		Pops boxed particles from the list used by the matrix image viewer (for boxes)
		'''
		#print "called delete display shapesS"
		if self.in_display_limbo: return
		
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
		
	def box_image_deleted(self,event,lc,force_image_mx_remove=True):
		box = self.delete_box(lc[0],force_image_mx_remove)
		self.boxable.add_exclusion_particle(box)
		guiim_core = self.guiim.get_core_object()
		guiim_core.set_other_data(self.boxable.get_exclusion_image(False),self.autoboxer.get_subsample_rate(),True)
		self.guictl.num_boxes_changed(len(self.ptcl))
		self.update_all_image_displays()
		#self.update_ppc()
		
	def delete_box(self,box_num,force_image_mx_remove=True):
		"""
		Deletes the numbered box completely
		Should only be called in the instance where a single box is being deleting - NOT when
		you are deleting a list of boxes sequentially (for that you should use delete_display_boxes
		and something to pop the box from the Boxable. See examples in this code)
		"""
		sh=self.guiim.getShapes()
		k=sh.keys()
		k.sort()
		del sh[int(box_num)]
		for j in k:
			if isinstance(j,int):
				if j>box_num :
					sh[j-1]=sh[j]
					del sh[j]
		self.guiim.delShapes()
		self.guiim.addShapes(sh)
		self.guiim.setActive(None,.9,.9,.4)
		
		if force_image_mx_remove: 
			#self.ptcl.pop(box_num)
			#self.guimx.setData(self.ptcl)
			self.guimx.get_core_object().pop_box_image(box_num)

		box = self.boxable.boxes[box_num]
		
		self.boxable.delete_box(box_num)
		# if the boxable was a reference then the autoboxer needs to be told. It will remove
		# it from its own list and potentially do autoboxing
		if box.isref:
			val = self.autoboxer.remove_reference(box)
			if val == 2:
				self.boxable.clear_and_cache(True)
				self.clear_displays()
				
		return box
	
	
	def box_display_update_specific(self, box_num,force=False):
		ns = {}
		box = self.get_boxes()[box_num]
		
		im=box.get_box_image()
		box.shape = EMShape([self.shape_string,box.r,box.g,box.b,box.xcorner,box.ycorner,box.xcorner+box.xsize,box.ycorner+box.ysize,2.0])
		
		#if not force and not box.isref and not box.ismanual:
			#box.shape.isanimated = True
			#box.shape.blend = 0
		
		ns[box_num]=box.shape
		if box_num>=len(self.ptcl) : self.ptcl.append(im)
		else : self.ptcl[box_num]=im
			
		self.guiim.addShapes(ns)
		self.set_ptcl_mx_data(self.ptcl)
		
		self.guictl.num_boxes_changed(len(self.ptcl))
		if self.guimxit != None:
			self.guimxit.get_target().set_extra_hud_data([len(self.ptcl)])
		
		
		if self.guimxitp !=None:
			othershapes = {}
			for shape in ns:
				s = ns[shape].getShape()
				othershapes[shape] = EMShape(["point",s[1],s[2],s[3],(s[4]+s[6])/2,(s[5]+s[7])/2,2])
	
			if isinstance(self.guimxit,EMImageRotor): self.guimxit.set_shapes(othershapes,self.get_image_thumb_shrink(),self.current_image_idx)

		self.update_all_image_displays()
		
		
	def box_display_update(self,force=False):
		
		ns = {}
		idx = 0
		#self.ptcl = []
		# get the boxes
		boxes =self.get_boxes()
		for j,box in enumerate(boxes):
			if not box.changed and not force:
				idx += 1
				continue
			
			box.changed=False
		
			im=box.get_box_image()
			box.shape = EMShape([self.shape_string,box.r,box.g,box.b,box.xcorner,box.ycorner,box.xcorner+box.xsize,box.ycorner+box.ysize,2.0])
			
			if not force and not box.isref and not box.ismanual:
				box.shape.isanimated = True
				box.shape.blend = 0
			ns[idx]=box.shape
			if idx>=len(self.ptcl) : self.ptcl.append(im)
			else : self.ptcl[idx]=im
			idx += 1
		
		self.guiim.addShapes(ns)
		self.set_ptcl_mx_data(self.ptcl)
		
		self.guictl.num_boxes_changed(len(self.ptcl))
		if self.guimxit != None:
			self.guimxit.get_target().set_extra_hud_data([len(self.ptcl)])

		if self.guimxitp !=None:
			othershapes = {}
			for shape in ns:
				s = ns[shape].getShape()
				othershapes[shape] = EMShape(["point",s[1],s[2],s[3],(s[4]+s[6])/2,(s[5]+s[7])/2,2])
	
			if isinstance(self.guimxit,EMImageRotor): self.guimxit.set_shapes(othershapes,self.get_image_thumb_shrink(),self.current_image_idx)
		self.update_all_image_displays()
		
	def run(self):
		"""If you make your own application outside of this object, you are free to use
		your own local app.exec_(). This is a convenience for boxer-only programs."""
		self.dynapixp.exec_()
		
		self.boxable.cache_exc_to_db()
		project_db = EMProjectDB()
		project_db.close()
		

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
		
		return (self.get_boxes())

	def force_autobox(self,bool):
		'''
		like hitting refresh - everything is updated
		'''
		self.guictl.pawel_histogram.clear()
		self.boxable.clear_and_cache(True)
		self.autoboxer.write_image_specific_references_to_db(self.boxable.get_image_name())
		self.boxable.get_references_from_db()
		self.autoboxer.auto_box(self.boxable, True,True)
		self.box_display_update()

	def toggle_dynapix(self,bool):
		self.dynapix = bool
		self.autoboxer.set_interactive_mode(self.dynapix)
		
	def done(self):
		self.boxable.cache_exc_to_db()
		self.dynapixp.quit
		
	def try_data(self,data,thr):
		print 'try that was pressed, this feature is currently disabled'
		
	def opt_params_updated(self,thresh,profile,radius):
		self.guictl.update_data(thresh,profile,radius)
		
	def set_selection_mode(self,selmode):
		if self.autoboxer.set_selection_mode(selmode):
			self.box_display_update()
		
	def set_profile_comparitor(self,cmp_mode):
		if self.autoboxer.set_cmp_mode(cmp_mode):
			self.box_display_update()
			
	def classify(self,bool):
		self.boxable.classify()
		
	def erase_toggled(self,bool):
		# for the time being there are only two mouse modes
		
		if bool == True:
			self.mouse_handlers["erasing"].set_mode(Boxable.ERASE)
			self.mousehandler = self.mouse_handlers["erasing"]
		else:
			self.guiim.addShape("eraser",EMShape(["circle",0,0,0,0,0,0,0.1]))
			self.update_image_display()
			self.mousehandler = self.mouse_handlers["boxing"]
			
	
	def unerase_toggled(self,bool):
		if bool == True:
			self.mouse_handlers["erasing"].set_mode(Boxable.UNERASE)
			self.mousehandler = self.mouse_handlers["erasing"]
		else:
			self.guiim.addShape("eraser",EMShape(["circle",0,0,0,0,0,0,0.1]))
			self.update_image_display()
			self.mousehandler = self.mouse_handlers["boxing"]
	
	def update_erase_rad(self,rad):
		self.mouse_handlers["erasing"].set_erase_radius(rad)

	def quit(self):
		self.dynapixp.quit()
	
	def set_dummy_box(self,box):
		self.autoboxer.set_dummy_box(box)
		self.box_display_update()
		
	def remove_dummy(self):
		self.autoboxer.set_dummy_box(None)
		self.box_display_update()
		
	def write_all_box_image_files(self,box_size,forceoverwrite=False,imageformat="hdf",normalize=True,norm_method="normalize.edgemean"):
		self.boxable.cache_exc_to_db()
		for image_name in self.image_names:
			
			
			try:
				project_db = EMProjectDB()
				data = project_db[get_idd_key(image_name)]
				trim_autoboxer = project_db[data["autoboxer_unique_id"]]["autoboxer"]
				autoboxer = SwarmAutoBoxer(self)
				autoboxer.become(trim_autoboxer)
				#print "writing box images for",image_name,"using",data["autoboxer_unique_id"]
			except:
				autoboxer = self.autoboxer
				#print "writing box images or",image_name,"using currently stored autoboxer"
				
			boxable = Boxable(image_name,self,autoboxer)
			boxable.clear_and_cache()
			if boxable.is_excluded():
				print "Image",image_name,"is excluded and being ignored"
				continue
			
			
			mode = self.autoboxer.get_mode()
			self.autoboxer.set_mode_explicit(SwarmAutoBoxer.COMMANDLINE)
			self.autoboxer.auto_box(boxable,False)
			self.autoboxer.set_mode_explicit(mode)
			
			boxable.write_box_images(box_size,forceoverwrite,imageformat,normalize,norm_method)
	
	def write_all_coord_files(self,box_size,forceoverwrite=False):
		self.boxable.cache_exc_to_db()
		for image_name in self.image_names:
			
			try:
				project_db = EMProjectDB()
				data = project_db[get_idd_key(image_name)]
				trim_autoboxer = project_db[data["autoboxer_unique_id"]]["autoboxer"]
				autoboxer = SwarmAutoBoxer(self)
				autoboxer.become(trim_autoboxer)
				#print "writing box coordinates for",image_name,"using",data["autoboxer_unique_id"]
			except:
				autoboxer = self.autoboxer
				#print "writing box coordinates for",image_name,"using currently stored autoboxer"
			boxable = Boxable(image_name,self,autoboxer)
			
			if boxable.is_excluded():
				print "Image",image_name,"is excluded and being ignored"
				continue
			
			mode = autoboxer.get_mode()
			autoboxer.set_mode_explicit(SwarmAutoBoxer.COMMANDLINE)
			autoboxer.auto_box(boxable,False)
			autoboxer.set_mode_explicit(mode)
			
			boxable.write_coord_file(box_size,forceoverwrite)
	
	def center(self,technique):
		
		if self.boxable.center(technique):
			self.box_display_update()
		else:
			print 'technique',technique,'is unsupported - check back tomorrow'
			
	def toggle_frozen(self):
		if self.boxable.is_excluded() : return
		self.boxable.toggle_frozen()
		self.boxable.write_to_db()
	
		self.guiim.get_core_object().set_frozen(self.boxable.is_frozen())
		
		if self.guimxitp != None  and isinstance(self.guimxit,EMImageRotor): self.guimxit.set_frozen(self.boxable.is_frozen(),self.current_image_idx)
		if not self.boxable.is_frozen():
			self.change_current_autoboxer(self.boxable.get_autoboxer_id(),False)
		else:	
			self.update_image_display()
		
	def change_image_quality(self,val):
		self.boxable.set_quality(val)
		if val == Boxable.EXCLUDE:
			self.boxable.set_frozen(False) # If it was already frozen then setting to excluded overrides this
			self.guiim.set_excluded(True)
			self.guiim.set_frozen(False)
			if self.guimxitp != None and isinstance(self.guimxit,EMImageRotor): self.guimxit.set_frozen(False,self.current_image_idx)
			self.boxable.clear_and_cache(True) # tell boxable to clear its autoboxes and references -
			self.clear_displays() # tell the display to clear itself
		elif self.guiim.set_excluded(False):
			self.update_image_display()
			
		self.boxable.write_to_db() # make sure the infromation changes that just occured are written to the DB
		

	def set_boxing_method(self,ref,manual):
		'''
		Okay could do it with one argument but leaving it this way makes it obvious
		'''
		
		if ref and manual:
			print 'error, you cant set both ref and manual'
			
		if ref:
			self.mouse_handlers["boxing"].set_mode(GUIbox.REFERENCE_ADDING)
		elif manual:
			self.mouse_handlers["boxing"].set_mode(GUIbox.MANUALLY_ADDING)
	
	def autoboxer_db_changed(self):
		self.guictl.update_ab_table()

	def add_new_autoboxer_db(self, n):
		if not self.boxable.is_interactive():
			return None
		autoboxer = SwarmAutoBoxer(self)
		autoboxer.set_box_size_explicit(self.box_size)
		autoboxer.set_interactive_mode(self.dynapix)
		autoboxer_db_string = "autoboxer_"+autoboxer.get_creation_ts()
		trim_autoboxer = TrimSwarmAutoBoxer(autoboxer)
		convenience_name = "New " + str(n)
		trim_autoboxer.set_convenience_name(convenience_name)
		trim_autoboxer.write_to_db()
		return convenience_name
	
	def add_copy_autoboxer_db(self,autoboxer_id,n):
		#print "adding a copy of the ab with id is",autoboxer_id
		project_db = EMProjectDB()
		trim_autoboxer = copy(project_db[autoboxer_id]["autoboxer"])
		trim_autoboxer.set_creation_ts(gm_time_string())
		convenience_name = "Copy " + str(n)
		trim_autoboxer.set_convenience_name(convenience_name)
		trim_autoboxer.write_to_db()
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
		self.current_image_idx_name = None
		self.image_names = []
		self.name_map = {}
		self.dict_data = {}
		self.set_current_image_name(parent.get_current_image_name())
		self.set_image_names(parent.get_image_names())
		
	def set_current_image_name(self,image_name):
		self.current_image_idx_name = image_name
		
	def set_image_names(self,image_names):
		self.image_names = image_names
		
	def get_autoboxer_data(self):
		project_db = EMProjectDB()
		self.dict_data = {}
		self.name_map = {}
		for i in project_db.keys():
			try:
				if i[0:10] == "autoboxer_":
					tag = project_db[i]["convenience_name"]
					self.dict_data[tag] = []
					self.name_map[i] = tag
			except:
				print "error couldn't handle",i
		
		for image_name in self.image_names:
			found = False
			try:
				data = project_db[get_idd_key(image_name)]
				#trim_autoboxer = project_db[data["autoboxer_unique_id"]]
				if data != None: found = True
			except: pass
			
			if found:
				try:
					self.dict_data[self.name_map[data["autoboxer_unique_id"]]].append(strip_after_dot(image_name))
				except: pass
					#print data
					#try:
						#print "error, an autoboxer has been lost, its stamp is",data["autoboxer_unique_id"]
					#except: pass
		return self.dict_data
	
	def get_current_autoboxer_ts(self):
		try:
			return self.name_map["autoboxer_"+self.parent.get_current_autoboxer_ts()]
		except:
			return None

	def add_new_autoboxer(self):
		return self.parent.add_new_autoboxer_db(self.get_total_autoboxers())

	def add_copy_autoboxer(self,tag):
		autoboxer_id = self.__get_autoboxer_id_from_tag(tag)
		if autoboxer_id != None:
			return self.parent.add_copy_autoboxer_db(autoboxer_id,self.get_total_autoboxers())
		else:
			print "error, couldn't find autoboxer from tag",tag
			return None
	
	def change_current_autoboxer(self,tag):
		# FIXME - is there any way to use a bidirectional map?pyt
		autoboxer_id = self.__get_autoboxer_id_from_tag(tag)
		if autoboxer_id != None:
			return self.parent.change_current_autoboxer(autoboxer_id)
		else:
			print "error, couldn't get autoboxer_id"
			return None
				
	def update_db_convenience_name(self,newname,oldname):
		'''
		Updates the unique name of the autoboxer in the DB
		if the name is already used then False is returned
		and the calling function should act on this to stop the
		name change, for example within a widget
		'''
		project_db = EMProjectDB()
		autoboxer_id = self.__get_autoboxer_id_from_tag(oldname)
		if autoboxer_id != None:
			self.name_map.pop(autoboxer_id)
			if self.__name_not_already_present(newname):
				self.name_map[autoboxer_id] = newname
				autoboxer = project_db[autoboxer_id]["autoboxer"]
				autoboxer.set_convenience_name(newname)
				autoboxer.write_to_db()
				return True
			else:
				self.name_map[autoboxer_id]["convenience_name"] = oldname
				return False
					
	def get_total_autoboxers(self):
		return len(self.name_map)
	
	def associated_images(self,tag):
		return self.dict_data[tag]
	
	def remove(self,tag):
		autoboxer_id = self.__get_autoboxer_id_from_tag(tag)
		project_db = EMProjectDB()
		project_db.pop(autoboxer_id)
		self.dict_data.pop(tag)
	
	def __name_not_already_present(self,name):
		for names in self.name_map.items():
			if names[1] == name:
			 	return False
		
		return True
	
	def __get_autoboxer_id_from_tag(self,tag):
		for names in self.name_map.items():
			if names[1] == tag:
				return names[0]
			
		print "error, couldn't find",tag,"in the namemap"
		return None
	
	def toggle_frozen(self,tag,bool):
		#frozen = boxable.is_frozen()
		#if frozen:
			#new_name = self.add_copy_autoboxer(tag)
			#self.get_autoboxer_data() # FIXME this is inefficient, could just add to self.dict_data etc
			#autoboxer_id = self.__get_autoboxer_id_from_tag(new_name)
			
			#boxable.set_autoboxer_id(autoboxer_id)
			#boxable.write_to_db()
		#boxable.write_to_db()
		self.parent.toggle_frozen()
		
	def clear_current(self):
		new_name = self.add_new_autoboxer()
		
		# if the new_name is none then the current Boxable is frozen!
		if new_name == None: return new_name
		
		self.get_autoboxer_data() # FIXME this is inefficient, could just add to self.dict_data etc
		autoboxer_id = self.__get_autoboxer_id_from_tag(new_name)
		boxable = self.parent.get_boxable()
		boxable.set_autoboxer_id(autoboxer_id)
		boxable.write_to_db()
		boxable.clear_and_cache(True)
		self.parent.clear_displays()
		return new_name


		
def histogram1d( data, nbin, presize=0 ) :
	fmax = max( data )
	fmin = min( data )
	binsize = (fmax - fmin)/(nbin-2*presize)
	start = fmin - binsize*presize
	region = [None]*nbin
	hist = [None]*nbin
	for i in xrange(nbin):
		region[i] = start + (i+0.5)*binsize
		hist[i] = 0

	for d in data:
		id = int( (d-start)/binsize )
		hist[id]+=1

	return region,hist



class CcfHistogram(QtGui.QWidget):

	def __init__(self, parent):	
		QtGui.QWidget.__init__(self,parent)
		self.parent=parent
		self.data=None
                self.setMinimumSize(QtCore.QSize(256,256))
		self.PRESIZE = 28

	def clear( self ):
		self.ccfs = None
		self.data = None
		self.parent.ccf_range.setText( "Range: (N/A, N/A)" )
		self.parent.threshold_low.setText( "N/A" )
		self.parent.threshold_hgh.setText( "N/A" )

	def setData( self, data ):
		self.ccfs = data
		self.nbin = self.width()
                self.data = histogram1d( data, self.nbin, self.PRESIZE )

		hmin = self.data[0][0]
		hmax = self.data[0][-1]

		info = "Range: (%8.4f, %8.4f)" % (hmin, hmax)

		self.parent.ccf_range.setText( info )
		self.parent.threshold_low.setText( str(hmin) )
		self.parent.threshold_hgh.setText( str(hmax) )
		
		self.tickers =[1, self.width()-2]

	def paintEvent( self, event ):
		p=QtGui.QPainter()
		p.begin(self)
		p.setBackground(QtGui.QColor(16,16,16))
		p.eraseRect(0,0,self.width(),self.height())
		p.setPen(Qt.darkGray)

		if self.data is None:
			return

		hmax = max( self.data[1] )
                for i in xrange( len(self.data[1]) ):
			h = self.data[1][i]
                        p.drawLine(i, self.height(), i, int(self.height()*(1-0.9*h/hmax)) )

		self.drawTicker( self.tickers[0] )
		self.drawTicker( self.tickers[1] )

	def mousePressEvent(self, event):
		if event.button()==Qt.LeftButton:
			x = event.x()
			dis1 = abs( x - self.tickers[0] )
			dis2 = abs( x - self.tickers[1] )

			if dis1 < dis2:
				self.cur_ticker = 0
			else:
				self.cur_ticker = 1

			if not hasattr( self, "shapes" ):
				self.shapes = self.parent.target.guiim.getShapes().copy()

	def mouseMoveEvent(self, event):
		if event.buttons()&Qt.LeftButton and event.x() > 0 :
			self.tickers[self.cur_ticker] = event.x()
			self.repaint()

			x = event.x()
			if x < 0 :
				thr = self.data[0][0]
			elif x > len(self.data[0]) : 
				thr = self.data[0][-1]
			else :
				thr = self.data[0][x]

			if self.cur_ticker==0:
				self.parent.threshold_low.setText( str(thr) )
			else:
				self.parent.threshold_hgh.setText( str(thr) )


			thr_low = float( self.parent.threshold_low.text() )
			thr_hgh = float( self.parent.threshold_hgh.text() )
	
			guiim = self.parent.target.guiim

			curt_shapes = guiim.getShapes()


			print "# of all shapes: ", len( self.shapes )
			print "# of cur shapes: ", len( curt_shapes )
			print "thr_low: ", thr_low
			print "thr_hgh: ", thr_hgh

			ndelete = 0
			for i in xrange( len(self.ccfs) ):
				score = self.ccfs[i] 

				if (score < thr_low or score > thr_hgh) and curt_shapes.has_key(i):
					ndelete += 1
					guiim.delShape( i )

				if score >= thr_low and score <= thr_hgh and not(curt_shapes.has_key(i)):
					guiim.addShape( i, self.shapes[i] )
			
			guiim.updateGL()


	def drawTicker( self, newpos ) :
		p=QtGui.QPainter()
		p.begin(self)
		p.setPen(Qt.yellow)
		for i in xrange( newpos, newpos+2):
			p.drawLine( i, self.height(), i, int(0.2*self.height()) )


		
class GUIboxPanel(QtGui.QWidget):
	def __init__(self,target,ab_sel_mediator) :
		
		QtGui.QWidget.__init__(self,None)
		self.target=target
		self.ab_sel_mediator = ab_sel_mediator
		
		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		
		self.tabwidget = QtGui.QTabWidget(self)
		self.insert_main_tab()
		self.insert_advanced_tab()
		self.insert_view_tab()
		self.vbl.addWidget(self.tabwidget)
		
		self.dummybox = Box()
		self.dummybox.isdummy = True
		self.currentlyselected = -1 # used in the ab_table
		
		self.lock = False
		self.connect(self.bs,QtCore.SIGNAL("editingFinished()"),self.new_box_size)
	
		self.connect(self.thr,QtCore.SIGNAL("valueChanged"),self.new_threshold)
		self.connect(self.done,QtCore.SIGNAL("clicked(bool)"),self.target.quit)
		self.connect(self.classifybut,QtCore.SIGNAL("clicked(bool)"),self.target.classify)
		self.connect(self.trythat,QtCore.SIGNAL("clicked(bool)"),self.try_dummy_parameters)
		self.connect(self.reset,QtCore.SIGNAL("clicked(bool)"),self.target.remove_dummy)
		self.connect(self.thrbut, QtCore.SIGNAL("clicked(bool)"), self.selection_mode_changed)
		
#		self.target.connect(self.target,QtCore.SIGNAL("nboxes"),self.num_boxes_changed)
	
	#def centerpushed(self,unused):
		#self.target.center(str(self.centerooptions.currentText()))
	
	def insert_main_tab(self):
		# this is the box layout that will store everything
		self.main_inspector = QtGui.QWidget()
		self.main_vbl =  QtGui.QVBoxLayout(self.main_inspector)
		
		self.infohbl = QtGui.QHBoxLayout()
		self.info = QtGui.QLabel("%d Boxes"%len(self.target.get_boxes()),self)
		#self.ppc = QtGui.QLabel("%f particles per click"%0,self)
		self.infohbl.addWidget(self.info)
		#self.infohbl.addWidget(self.ppc)
		
		
		self.statsbox = QtGui.QGroupBox("Stats")
		self.statsbox.setLayout(self.infohbl)
		self.main_vbl.addWidget(self.statsbox)
		
		self.boxingvbl = QtGui.QVBoxLayout()
		
		self.boxinghbl1=QtGui.QHBoxLayout()
		self.boxinghbl1.setMargin(0)
		self.boxinghbl1.setSpacing(2)
		
		self.refbutton=QtGui.QPushButton( QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/black_box.png"), "Reference")
		self.refbutton.setCheckable(1)
		self.refbutton.setChecked(True)
		self.boxinghbl1.addWidget(self.refbutton)
		
		self.manualbutton=QtGui.QPushButton(QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/white_box.png"), "Manual")
		self.manualbutton.setCheckable(1)
		self.manualbutton.setChecked(False)
		self.boxinghbl1.addWidget(self.manualbutton)
		
		self.lblbs=QtGui.QLabel("Box Size:",self)
		self.boxinghbl1.addWidget(self.lblbs)
		self.pos_int_validator = QtGui.QIntValidator(self)
		self.pos_int_validator.setBottom(1)
		self.bs = QtGui.QLineEdit(str(self.target.box_size),self)
		self.bs.setValidator(self.pos_int_validator)
		self.boxinghbl1.addWidget(self.bs)
		
		self.boxingvbl.addLayout(self.boxinghbl1)
	
		self.boxinghbl3=QtGui.QHBoxLayout()
		self.dynapix = QtGui.QCheckBox("Dynapix")
		self.dynapix.setChecked(self.target.dynapix)
		self.boxinghbl3.addWidget(self.dynapix)

		self.method=QtGui.QComboBox()
		self.swarm_icon = QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/swarm_icon.png")
		self.method.addItem( self.swarm_icon, "Swarm" )
		self.setWindowIcon( self.swarm_icon )
		self.pp_icon = QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/pp_boxer_icon.png");
		self.method.addItem( self.pp_icon,"Gauss Conv" )
		self.boxinghbl3.addWidget( self.method )

		self.autobox=QtGui.QPushButton(QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/green_boxes.png"), "Autobox")
		self.boxinghbl3.addWidget(self.autobox)
		self.boxingvbl.addLayout(self.boxinghbl3)
	
		self.boxinghbl2=QtGui.QHBoxLayout()
		self.boxinghbl2.setMargin(2)
		self.boxinghbl2.setSpacing(6)
		#self.vbl.addLayout(self.hbl1)
		
		self.erasepic = QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/boxer_erase.png");
		self.erase=QtGui.QPushButton(self.erasepic,Boxable.ERASE)
		self.erase.setCheckable(1)
		self.boxinghbl2.addWidget(self.erase)
		
		self.unerasepic = QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/boxer_unerase.png");
		self.unerase=QtGui.QPushButton(self.unerasepic,Boxable.UNERASE)
		self.unerase.setCheckable(1)
		self.boxinghbl2.addWidget(self.unerase)
		
		self.eraseradtext=QtGui.QLabel("Erase Radius:",self)
		self.boxinghbl2.addWidget(self.eraseradtext)
		
		self.eraserad = QtGui.QLineEdit(str(self.target.eraseradius),self)
		self.boxinghbl2.addWidget(self.eraserad)
		self.eraserad.setEnabled(False)
		
		self.boxingvbl.addLayout(self.boxinghbl2)

		self.invert_contrast_mic = QtGui.QCheckBox("Invert Contrast")
		self.invert_contrast_mic.setChecked(True)
		self.boxingvbl.addWidget(self.invert_contrast_mic,0, Qt.AlignLeft)
		self.connect(self.invert_contrast_mic,QtCore.SIGNAL("clicked(bool)"),self.invert_contrast_mic_toggled)

		self.boxinghbl4=QtGui.QHBoxLayout()
		self.togfreeze=QtGui.QPushButton(QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/freeze_swirl.png"),"Toggle Freeze")
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
		self.write_all_box_image_files = QtGui.QPushButton("Write Box Images")
		self.outputhbl1.addWidget(self.write_all_box_image_files,0)
		self.normalize_box_images = QtGui.QCheckBox("Normalize Box Images")
		self.normalize_box_images.setChecked(True)
		self.outputhbl1.addWidget(self.normalize_box_images,0, Qt.AlignLeft)
		self.outputvbl.addLayout(self.outputhbl1)
		
		self.outputhbl2=QtGui.QHBoxLayout()
		self.write_all_coord_files = QtGui.QPushButton("Write Coord Files")
		self.outputhbl2.addWidget(self.write_all_coord_files,0)
		self.outputforceoverwrite = QtGui.QCheckBox("Force Overwrite")
		self.outputforceoverwrite.setChecked(False)
		self.outputhbl2.addWidget(self.outputforceoverwrite,0, Qt.AlignRight)
		self.outputvbl.addLayout(self.outputhbl2)
		
		self.output_gridl = QtGui.QGridLayout()
		
		self.usingbox_sizetext=QtGui.QLabel("Using Box Size:",self)
		self.output_gridl.addWidget(self.usingbox_sizetext,0,0,Qt.AlignRight)
		
		self.usingbox_size = QtGui.QLineEdit(str(self.target.box_size),self)
		self.usingbox_size.setValidator(self.pos_int_validator)
		self.output_gridl.addWidget(self.usingbox_size,0,1,Qt.AlignLeft)
		
		self.outputformat=QtGui.QLabel("Image Format:",self)
		self.output_gridl.addWidget(self.outputformat,1,0, Qt.AlignRight)
		
		self.outputformats = QtGui.QComboBox(self)
		self.outputformats.addItem("hdf")
		self.outputformats.addItem("img")
		self.output_gridl.addWidget(self.outputformats,1,1,Qt.AlignLeft)
		
		self.normalization_method=QtGui.QLabel("Normalization Method:",self)
		self.output_gridl.addWidget(self.normalization_method,2,0,Qt.AlignRight)
		
		self.normalization_options = QtGui.QComboBox(self)
		self.normalization_options.addItem("normalize.edgemean")
		self.normalization_options.addItem("normalize.ramp.normvar")
		self.output_gridl.addWidget(self.normalization_options,2,1, Qt.AlignLeft)
		
		self.outputvbl.addLayout(self.output_gridl)
		
		self.outputbox = QtGui.QGroupBox("Output")
		self.outputbox.setLayout(self.outputvbl)
		self.main_vbl.addWidget(self.outputbox)

		self.classifybut=QtGui.QPushButton("Classify")
		self.main_vbl.addWidget(self.classifybut)
		
		self.done=QtGui.QPushButton("Done")
		self.main_vbl.addWidget(self.done)
		
		self.tabwidget.addTab(self.main_inspector,"Main")
		
		self.connect(self.eraserad,QtCore.SIGNAL("editingFinished()"),self.update_erase_rad)
		self.connect(self.erase, QtCore.SIGNAL("clicked(bool)"), self.erase_toggled)
		self.connect(self.unerase, QtCore.SIGNAL("clicked(bool)"), self.unerase_toggled)
		self.connect(self.autobox,QtCore.SIGNAL("clicked(bool)"),self.target.force_autobox)
		self.connect(self.togfreeze,QtCore.SIGNAL("clicked(bool)"),self.toggle_frozen)
		self.connect(self.clear,QtCore.SIGNAL("clicked(bool)"),self.clear_current)
		self.connect(self.dynapix,QtCore.SIGNAL("clicked(bool)"),self.dynapix_toggled)
		
		self.connect(self.refbutton, QtCore.SIGNAL("clicked(bool)"), self.ref_button_toggled)
		self.connect(self.manualbutton, QtCore.SIGNAL("clicked(bool)"), self.manual_button_toggled)

		self.connect(self.write_all_box_image_files,QtCore.SIGNAL("clicked(bool)"),self.write_box_images)
		self.connect(self.write_all_coord_files,QtCore.SIGNAL("clicked(bool)"),self.write_box_coords)
		
		self.connect(self.method, QtCore.SIGNAL("activated(int)"), self.method_changed)

		self.connect(self.normalize_box_images,QtCore.SIGNAL("clicked(bool)"),self.normalize_box_images_toggled)

		QtCore.QObject.connect(self.imagequalities, QtCore.SIGNAL("currentIndexChanged(QString)"), self.image_quality_changed)

	def normalize_box_images_toggled(self):
		val = self.normalize_box_images.isChecked()
		self.normalization_options.setEnabled(val)
		self.normalization_method.setEnabled(val)
	
	def invert_contrast_mic_toggled(self):
		img = self.target.guiim.image2d.data
		avg = img.get_attr("mean")
		invimg = img*(-1.0) #+ 2.0*avg
		self.target.guiim.setData(invimg)
		self.target.guiim.updateGL()
	
	def method_changed(self, methodid):

		name = self.method.itemText( methodid )

		if name[0:5] == "Swarm":
			tabid = self.tabwidget.indexOf( self.pawel_option )
			if tabid != -1:
				self.tabwidget.removeTab( tabid )
				self.tabwidget.insertTab( tabid, self.david_option, "Swarm Advanced" )
				self.autobox.setText("Autobox")
				self.setWindowIcon( self.swarm_icon )

			#self.target.autoboxer = SwarmAutoBoxer(self.target)
			#self.target.autoboxer.set_box_size_explicit(self.target.box_size)
			#self.target.autoboxer.set_interactive_mode(self.target.dynapix)
			self.target.set_autoboxer(self.target.image_names[0])
		else:
			assert name[0:5]=="Gauss"
			tabid = self.tabwidget.indexOf( self.david_option )
			if tabid != -1:
				self.tabwidget.removeTab( tabid )
				self.tabwidget.insertTab( tabid, self.pawel_option, "Gauss Advanced" )
				self.autobox.setText("Run")
				self.setWindowIcon( self.pp_icon )
				
			self.target.autoboxer = PawelAutoBoxer(self.target)


	def set_image_quality(self,integer):
		self.lock = True
		self.imagequalities.setCurrentIndex(int(integer))
		self.lock = False
		
	def image_quality_changed(self,val):
		if self.lock == False:
			self.target.change_image_quality(str(val))

	def clear_current(self,unused):
		self.lock = True
		new_name = self.ab_sel_mediator.clear_current()
		self.lock = False
		if new_name != None: # if the boxable wasn't frozen...
			self.update_ab_table()
			self.setChecked(new_name)
		
	def toggle_frozen(self,bool):
		self.lock = True
		self.ab_sel_mediator.toggle_frozen(self.col1[self.currentlyselected].text(),bool)
		self.lock = False
	
	def write_box_images(self,unused):
		box_size = int(str(self.usingbox_size.text()))
		realbox_size = int(str(self.bs.text()))
		if realbox_size == box_size:
			box_size = -1 # negative one is a flag that tells the boxes they don't need to be resized... all the way in the Box Class
		self.target.write_all_box_image_files(box_size,self.outputforceoverwrite.isChecked(),str(self.outputformats.currentText()),self.normalize_box_images.isChecked(),str(self.normalization_options.currentText()))
		
	def write_box_coords(self,unused):
		box_size = int(str(self.usingbox_size.text()))
		realbox_size = int(str(self.bs.text()))
		if realbox_size == box_size:
			box_size = -1 # negative one is a flag that tells the boxes they don't need to be resized... all the way in the Box Class
		self.target.write_all_coord_files(box_size,self.outputforceoverwrite.isChecked(),)
	
	def setChecked(self,tag):
		
		for i,col in enumerate(self.col1):
			if str(col.text()) == tag:
				col.setCheckState(Qt.Checked)
				self.currentlychecked = i
				break
	
	def ab_table_cell_changed(self,i,j):
		if i >= len(self.col1): return
		if self.lock:
			return
		#data = self.ab_sel_mediator.get_autoboxer_data()
		#data = data.items()
		if str(self.col1[i].text()) != self.colnames[i]:
			if not self.ab_sel_mediator.update_db_convenience_name(str(self.col1[i].text()),self.colnames[i]):
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
				self.ab_sel_mediator.change_current_autoboxer(str(self.col1[i].text()))
				self.update_ab_table()
			else:
				print "error, unforeseen checkstate circumstance. Nothing done"
			#except: pass
		
		self.lock = False
				
	def ab_table_item_changed(self,item):
		print "item changed"
	
	def update_ab_table(self):
		
		data = self.ab_sel_mediator.get_autoboxer_data()
		self.col1 = []
		self.col2 = []
		self.colnames = []
		if len(data) == 0:
			return
			
		self.lock = True
		idx = 0
		for d in data.items():
			if idx >= self.ab_table.rowCount():
				self.ab_table.insertRow(idx)
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
			self.ab_table.setItem(idx,0,col1)
			self.ab_table.setItem(idx,1,col2)
			self.col1.append(col1)
			self.col2.append(col2)
			self.colnames.append(col1.text())
			idx += 1
		
		# remove any overhanging columns if they already existed
		while (len(data) < self.ab_table.rowCount()):
			self.ab_table.removeRow(self.ab_table.rowCount()-1)
	
		currentsel = self.ab_sel_mediator.get_current_autoboxer_ts()
		self.currentlychecked = -1
		
		self.setChecked(currentsel)
		self.lock = False
	
	def delete_autoboxer(self,unused):
		items = self.ab_table.selectedItems()
		data = self.ab_sel_mediator.get_autoboxer_data()
		update = False
		for item in items:
			if item.column() == 0:
				if len(self.ab_sel_mediator.associated_images(str(item.text()))) != 0:
					print "can't delete an autoboxer unless it has no images associated with it"
				else:
					self.ab_sel_mediator.remove(str(item.text()))
					update = True
					
		if update:
			self.update_ab_table()
	
	def add_new_autoboxer(self,bool):
		self.ab_sel_mediator.add_new_autoboxer()
		self.update_ab_table()
		
	def add_copy_autoboxer(self,bool):
		numsel = 0
		selected = ""
		for col in self.col1:
			if self.ab_table.isItemSelected(col):
				numsel += 1
				selected = str(col.text())
				if numsel > 1:
					print "error, more than one autoboxer is selected. Please choose only one"
					return
		if numsel == 0:
			print "no autoboxers were selected, doing nothing"
			return
		else:
			self.ab_sel_mediator.add_copy_autoboxer(selected)
			self.lock=False
			self.update_ab_table()
			
	def insert_view_tab(self):
		# this is the box layout that will store everything
		self.view_inspector = QtGui.QWidget()
		self.view_vbl =  QtGui.QVBoxLayout(self.view_inspector)
		
		#  Insert the plot widget
		self.viewhbl = QtGui.QHBoxLayout()
		
		self.viewboxes = QtGui.QCheckBox("Boxed Particle Window")
		self.viewboxes.setChecked(True)
		self.viewimage = QtGui.QCheckBox("Main Image Window")
		self.viewimage.setChecked(True)
		
		self.viewhbl.addWidget(self.viewboxes)
		self.viewhbl.addWidget(self.viewimage)
		
		if self.target.has_thumbnails():
			self.viewthumbs = QtGui.QCheckBox("Image Thumbnails Window")
			self.viewthumbs.setChecked(True)
			self.viewhbl.addWidget(self.viewthumbs)
		
		self.viewmanagement = QtGui.QGroupBox("Displayed Windows")
		self.viewmanagement.setLayout(self.viewhbl)
		self.view_vbl.addWidget(self.viewmanagement)
		
			
		self.viewhbl2 = QtGui.QHBoxLayout()
		self.boxdisplay =QtGui.QLabel("Box Display Object:",self)
		self.viewhbl2.addWidget(self.boxdisplay)
		
		self.boxformats = QtGui.QComboBox(self)
		self.boxformats.addItem("square with central dot")
		self.boxformats.addItem("square")
		self.boxformats.addItem("circle with central dot")
		self.boxformats.addItem("circle")
		self.viewhbl2.addWidget(self.boxformats)
		
		self.displayboxes = QtGui.QGroupBox("Displayed Boxes")
		self.displayboxes.setLayout(self.viewhbl2)
		self.view_vbl.addWidget(self.displayboxes)
		
		self.tabwidget.addTab(self.view_inspector,"Display Options")
		
		self.connect(self.viewboxes,QtCore.SIGNAL("clicked(bool)"),self.target.view_boxes_clicked)
		self.connect(self.viewimage,QtCore.SIGNAL("clicked(bool)"),self.target.view_image_clicked)
		if self.target.has_thumbnails():
			self.connect(self.viewthumbs,QtCore.SIGNAL("clicked(bool)"),self.target.view_thumbs_clicked)
			
		QtCore.QObject.connect(self.boxformats, QtCore.SIGNAL("currentIndexChanged(QString)"), self.box_format_changed)
	
	def box_format_changed(self,new_format):
		format = str(new_format)
		if format == "square with central dot": format = "rectpoint"
		elif format == "square": format = "rect"
		elif format == "circle with central dot": format = "rcirclepoint"
		elif format == "circle": format = "rcircle"
		else: 
			print "errror, unknown format" 
			return
		print "format changed to ", format
		self.target.change_shapes(format)

	
	def insert_advanced_tab(self):
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
		self.ratio_average_but = QtGui.QRadioButton("Average Ratio")

		self.cmpmethodhbox = QtGui.QHBoxLayout()
		self.cmpmethodhbox.addWidget(self.ratiobut)
		self.cmpmethodhbox.addWidget(self.difbut)
		self.cmpmethodhbox.addWidget(self.ratio_average_but)
		
		self.cmpgroupbox = QtGui.QGroupBox("Peak Profile Comparitor")
		self.cmpgroupbox.setLayout(self.cmpmethodhbox)
		

		self.advanced_vbl.addWidget(self.cmpgroupbox)

		self.lock = True
		self.autoboxerhdbl = QtGui.QHBoxLayout()
		# ab means autoboxer
		self.ab_table = QtGui.QTableWidget(1,2,self)
		self.ab_table.setColumnWidth(1,150)
		self.abcol0title = QtGui.QTableWidgetItem("Autoboxer ID")
		self.abcol1title = QtGui.QTableWidgetItem("Associated Images")
		self.update_ab_table()
		self.lock = True
		self.ab_table.setHorizontalHeaderItem(0,self.abcol0title)
		self.ab_table.setHorizontalHeaderItem(1,self.abcol1title)
		self.autoboxerhdbl.addWidget(self.ab_table)
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
			

		self.david_option = self.adv_inspector
		
		self.pawel_option = QtGui.QWidget()
		self.pawel_option_vbox = QtGui.QVBoxLayout(self.pawel_option)
		self.pawel_option_vbox.addWidget(QtGui.QLabel("Gauss Conv's Parameters") )
		pawel_grid1 = QtGui.QGridLayout( )
		self.pawel_option_vbox.addLayout(pawel_grid1)
		
		pawel_grid1.addWidget( QtGui.QLabel("Input Pixel Size:") , 0, 0 )
		pawel_grid1.addWidget( QtGui.QLabel("Output Pixel Size:"), 1, 0 )
		pawel_grid1.addWidget( QtGui.QLabel("Angstrom"), 0, 2 )
		pawel_grid1.addWidget( QtGui.QLabel("Angstrom"), 1, 2 )		

		self.input_pixel_size = QtGui.QLineEdit("1.0", self)
		self.output_pixel_size = QtGui.QLineEdit("1.0", self)
		pawel_grid1.addWidget( self.input_pixel_size, 0, 1 )
		pawel_grid1.addWidget( self.output_pixel_size, 1, 1 )
	

		self.pawel_option_vbox.addWidget( QtGui.QLabel("CCF Histogram") )
		self.pawel_histogram = CcfHistogram( self )
		self.pawel_option_vbox.addWidget( self.pawel_histogram )

		self.ccf_range = QtGui.QLabel("Range: (N/A, N/A)")
		self.pawel_option_vbox.addWidget( self.ccf_range )

		pawel_grid2 = QtGui.QGridLayout()
		self.pawel_option_vbox.addLayout( pawel_grid2 )

		pawel_grid2.addWidget( QtGui.QLabel("Threshold Low:"),  0, 0 )
		pawel_grid2.addWidget( QtGui.QLabel("Threshold High:"), 1, 0 )
		self.threshold_low = QtGui.QLineEdit()
		self.threshold_hgh = QtGui.QLineEdit()
		pawel_grid2.addWidget( self.threshold_low, 0, 1 )
		pawel_grid2.addWidget( self.threshold_hgh, 1, 1 )

		#self.pawel_table.setVerticalHeaderLabels( ["Input Pixel Size  (Angstrom): ", "Output Pixel size (Angstrom): "] )
		#self.pawel_table.horizontalHeader().hide()
		#self.pawel_table.setItem(0, 0, QtGui.QTableWidgetItem("1.0") )
		#self.pawel_table.setItem(1, 0, QtGui.QTableWidgetItem("1.0") )
		#self.connect(self.pawel_table, QtCore.SIGNAL("cellChanged(int,int)"), self.pawel_parm_changed)


		self.tabwidget.addTab(self.david_option,"Swarm Advanced")
		self.connect(self.abnew, QtCore.SIGNAL("clicked(bool)"), self.add_new_autoboxer)
		self.connect(self.abcopy, QtCore.SIGNAL("clicked(bool)"), self.add_copy_autoboxer)
		self.connect(self.abdelete, QtCore.SIGNAL("clicked(bool)"), self.delete_autoboxer)
		self.connect(self.ab_table, QtCore.SIGNAL("itemChanged(QtGui.QTableWidgetItem)"), self.ab_table_item_changed)
		self.connect(self.ab_table, QtCore.SIGNAL("cellChanged(int,int)"), self.ab_table_cell_changed)
		self.connect(self.selbut, QtCore.SIGNAL("clicked(bool)"), self.selection_mode_changed)
		self.connect(self.morselbut, QtCore.SIGNAL("clicked(bool)"), self.selection_mode_changed)
		self.connect(self.ratiobut, QtCore.SIGNAL("clicked(bool)"), self.cmp_box_changed)
		self.connect(self.ratio_average_but, QtCore.SIGNAL("clicked(bool)"), self.cmp_box_changed)
		#self.connect(self.centerbutton,QtCore.SIGNAL("clicked(bool)"),self.centerpushed)
		self.connect(self.difbut, QtCore.SIGNAL("clicked(bool)"), self.cmp_box_changed)
	
	def pawel_parm_changed(self, row, col ):
		from string import atof

		assert col==0

		t = self.pawel_table.item(row, col).text()

		print 'row, col, text: ', row, col, t

		if row==0:
			self.target.autoboxer.pixel_input = atof( t )
		else:
			assert row==1
			self.target.autoboxer.pixel_output = atof( t )		

	def set_dynapix(self,bool):
		self.dynapix.setChecked(bool)

	def ref_button_toggled(self,bool):
		
		if self.refbutton.isChecked():
			self.manualbutton.setChecked(False)
		
		if not self.refbutton.isChecked():
			self.manualbutton.setChecked(True)

		self.target.set_boxing_method(self.refbutton.isChecked(),self.manualbutton.isChecked())
		
	def manual_button_toggled(self,bool):
		if self.manualbutton.isChecked():
			self.refbutton.setChecked(False)
		
		if not self.manualbutton.isChecked():
			self.refbutton.setChecked(True)
			
		self.target.set_boxing_method(self.refbutton.isChecked(),self.manualbutton.isChecked())

	def erase_toggled(self,bool):
		self.unerase.setChecked(False)
		self.eraserad.setEnabled(bool)
		self.target.guiim.setMouseTracking(bool)
		self.target.erase_toggled(bool)
		
	def unerase_toggled(self,bool):
		self.erase.setChecked(False)
		self.eraserad.setEnabled(bool)
		self.target.guiim.setMouseTracking(bool)
		self.target.unerase_toggled(bool)
		
	def dynapix_toggled(self,bool):
		self.target.toggle_dynapix(bool)
	
	def cmp_box_changed(self,unusedbool):
		if self.ratiobut.isChecked():
			s = BoxingTools.CmpMode.SWARM_RATIO
		elif self.difbut.isChecked():
			s = BoxingTools.CmpMode.SWARM_DIFFERENCE
		elif self.ratio_average_but.isChecked():
			s = BoxingTools.CmpMode.SWARM_AVERAGE_RATIO
		else:
			print "Bug intercepted in e2boxer.py. Please email the development team."
			
		self.target.set_profile_comparitor(s)
	
	def selection_mode_changed(self,unusedbool):
		if self.thrbut.isChecked():
			s = self.thrbut.text()
		elif self.selbut.isChecked():
			s = self.selbut.text()
		elif self.morselbut.isChecked():
			s = self.morselbut.text()
		else:
			print "Bug intercepted in e2boxer.py. Please email the development team."
			
		self.target.set_selection_mode(str(s))
	
	def try_dummy_parameters(self):
		print "option currently disabled"
		return
		self.dummybox.set_opt_profile(self.window.getData())
		self.dummybox.set_correlation_score(float(self.thr.getValue()))
		self.target.set_dummy_box(self.dummybox)
	
	def update_data(self,thresh,data,datar):
		#print data
		self.window.setData(data,datar)
		self.thr.setValue(thresh,True)
		self.resize(self.width(),self.height())
		#self.window.resizeGL(self.window.width(),self.window.height())
		#self.window.updateGL()
	def num_boxes_changed(self,n):
		self.info.setText("%d Boxes"%n)
		
	#def ppc_changed(self,f):
		#self.ppc.setText("%f ppc"%f)
	
	def adjust_erase_rad(self,delta):
		v = float(self.eraserad.text())
		if delta > 0:
			v = 1.1*v
		if delta < 0:
			v = 0.9*v
			
		self.eraserad.setText(str(int(v)))
		# this makes sure the target updates itself 
		# there may be a better approach, seeing as
		# the target called this function
		self.update_erase_rad()
		
	def update_erase_rad(self):
		v = int(self.eraserad.text())
		if ( v < 1 ): raise Exception
		self.target.update_erase_rad(v)
	
	def new_box_size(self):
		try:
			v=int(self.bs.text())
			if v<12 : raise Exception
		except:
			self.bs.setText(str(self.target.box_size))
			return
		
		
		self.usingbox_size.setText(self.bs.text())
		app = QtGui.QApplication.instance()
		app.setOverrideCursor(Qt.BusyCursor)
		self.target.update_box_size(v,1)
		app.setOverrideCursor(Qt.ArrowCursor)
	
	def set_box_size(self,box_size):
		self.bs.setText(str(box_size))
		self.usingbox_size.setText(str(box_size))
	
	def new_threshold(self,val):
		#print "new threshold"
		self.try_dummy_parameters()

if __name__ == "__main__":
	main()
