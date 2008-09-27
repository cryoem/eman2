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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA
#

import PyQt4
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
#from OpenGL import GL,GLU,GLUT
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
from valslider import ValSlider
from math import *
from EMAN2 import *
import EMAN2
import sys
import numpy
from emimageutil import ImgHistogram
from weakref import WeakKeyDictionary
from pickle import dumps,loads
from PyQt4.QtCore import QTimer
from time import time
import math

from emfloatingwidgets import *
from emglobjects import Camera, viewport_width, viewport_height,resize_gl
from e2boxer import *



class EMWindowNode:
	def __init__(self,parent):
		self.parent = parent
		self.children = []
		
	def attach_child(self,new_child):
		for child in self.children:
			if (child == new_child):
				print "error, can not attach the same child to the same parent more than once"
				return
		
		self.children.append(new_child)
	
	def set_parent(self,parent):
		self.parent = parent
		
	def parent_width(self):
		return parent.width()
	
	def parent_height(self):
		return parent.height()

class EMRegion:
	def __init__(self,geometry=Region(0,0,0,0,0,0)):
		self.geometry = geometry
		
	def set_geometry(self,geometry):
		self.geometry = geometry
		
	def get_geomoetry(self): return self.geometry
	
	def width(self):
		return self.geometry.get_width()
	
	def height(self):
		return self.geometry.get_height()
		
	def depth(self):
		return self.geometry.get_depth()
	
	def set_width(self,v):
		self.geometry.set_width(v)
	
	def set_height(self,v):
		self.geometry.set_height(v)
		
	def set_depth(self,v):
		self.geometry.set_depth(v)
	
	def get_size(self):
		return self.geometry.get_size()
	
	def get_origin(self):
		return self.geometry.get_origin()
	
	def set_origin(self,v):
		return self.geometry.set_origin(v)

class Animatable:
	cache_dts = None
	def __init__(self):
		self.time = 0		# time indicates the current time used for the basis of animation.
		self.time_interval = 0.2 # 0.5 seconds for the animation to complete
		self.inverse_time_inverval = 1.0/self.time_interval
		self.time_begin = 0 # records the time at which the animation was begun
		self.animated = True
		self.n = 20
		if Animatable.cache_dts == None:
			self.init_cache_dts()
		
	def init_cache_dts(self):
		Animatable.cache_dts = []
		for i in range(self.n):
			Animatable.cache_dts.append(sin(float(i)/(self.n-1)*math.pi/2))
		
	def set_animated(self,val=True):
		self.animated = val
		
	def is_animated(self): return self.animated
	
	def animate(self,t):
		if not self.animated: return False
		
		if self.time_begin == 0:
			self.time_begin = t
			
		self.time = t -self.time_begin

		if self.time > self.time_interval:
			self.time_begin = 0
			self.time = 0
			self.animated = False
			return 0
		else:
			dt = self.__get_dt()
			self.calculate_animation(dt)
			return 1
	
	def calculate_animation(self,dt): raise
	
	def __get_dt(self):
		idx = int( self.n*self.time*self.inverse_time_inverval)
		#print Animatable.cache_dts
		return Animatable.cache_dts[idx]
		
		
	def transform(self): raise
	
class SingleAxisRotationAnimation(Animatable):
	def __init__(self,target,start,end,axis=[0,1,0]):
		Animatable.__init__(self)
		self.target = target
		self.start = start
		self.current = start
		self.end = end
		self.axis = axis
		
	def get_start(self):
		return self.start
	
	def get_end(self):
		return self.end
	
	def get_current(self):
		return self.current
	
	def calculate_animation(self,dt):
		'''
		based on the assumption that dt goes from 0 to 1
		'''
		if dt > 1: raise
		self.current = (1-dt)*self.start + dt*self.end
		self.target.set_rotation(self.current)
	#def transform(self):
		#glRotate(self.current,*self.axis)
		
class XYScaleAnimation(Animatable):
	def __init__(self,target,start,end):
		Animatable.__init__(self)
		self.start = start
		self.current = start
		self.end = end
		self.target = target
		
	def get_start(self):
		return self.start
	
	def get_end(self):
		return self.end
	
	def get_current(self):
		return self.current
	
	def calculate_animation(self,dt):
		'''
		based on the assumption that dt goes from 0 to 1
		'''
		if dt > 1: raise
		self.current = (1-dt)*self.start + dt*self.end
		self.target.set_xy_scale(self.current)
	def transform(self):
		glScale(self.current,self.current,1.0)

class EMGLViewContainer(EMWindowNode,EMRegion):
	def __init__(self,parent,geometry=Region(0,0,0,0,0,0)):
		EMWindowNode.__init__(self,parent)
		EMRegion.__init__(self,geometry)
		self.current = None
		self.previous = None
		
	def draw(self):
		for child in self.children:
			glPushMatrix()
			child.draw()
			glPopMatrix()
	
	def updateGL(self):
		self.parent.updateGL()
		
	def resizeEvent(self, width, height):
		for child in self.children:
			child.set_update_P_inv()
	
	def mousePressEvent(self, event):
		for child in self.children:
			if ( child.isinwin(event.x(),viewport_height()-event.y()) ):
				child.mousePressEvent(event)
				self.updateGL()
				return
	
	def mouseMoveEvent(self, event):
		for child in self.children:
			if ( child.isinwin(event.x(),viewport_height()-event.y()) ):
				self.current = child
				if (self.current != self.previous ):
					if ( self.previous != None ):
						try: self.previous.leaveEvent()
						except: pass
				child.mouseMoveEvent(event)
				self.previous = child
				self.updateGL()
				return
		
	def mouseReleaseEvent(self, event):
		for child in self.children:
			if ( child.isinwin(event.x(),viewport_height()-event.y()) ):
				child.mouseReleaseEvent(event)
				self.updateGL()
				return
					
		
	def mouseDoubleClickEvent(self, event):
		for child in self.children:
			if ( child.isinwin(event.x(),viewport_height()-event.y()) ):
				child.mouseDoubleClickEvent(event)
				self.updateGL()
				return
		
		
	def wheelEvent(self, event):
		for child in self.children:
			if ( child.isinwin(event.x(),viewport_height()-event.y()) ):
				child.wheelEvent(event)
				self.updateGL()
				return

	def toolTipEvent(self, event):
		for child in self.children:
			if ( child.isinwin(event.x(),viewport_height()-event.y()) ):
				child.toolTipEvent(event)
				self.updateGL()
				return
		
		QtGui.QToolTip.hideText()

	def dragMoveEvent(self,event):
		print "received drag move event"
		
	def event(self,event):
		#print "event"
		#QtGui.QToolTip.hideText()
		if event.type() == QtCore.QEvent.MouseButtonPress: 
			self.mousePressEvent(event)
			return True
		elif event.type() == QtCore.QEvent.MouseButtonRelease:
			self.mouseReleaseEvent(event)
			return True
		elif event.type() == QtCore.QEvent.MouseMove: 
			self.mouseMoveEvent(event)
			return True
		elif event.type() == QtCore.QEvent.MouseButtonDblClick: 
			self.mouseDoubleClickEvent(event)
			return True
		elif event.type() == QtCore.QEvent.Wheel: 
			self.wheelEvent(event)
			return True
		elif event.type() == QtCore.QEvent.ToolTip: 
			self.toolTipEvent(event)
			return True
		else: 
			return QtOpenGL.QGLWidget.event(self,event)

	def hoverEvent(self,event):
		#print "hoverEvent
		for child in self.children:
			if ( child.isinwin(event.x(),self.height()-event.y()) ):
				child.hoverEvent(event)
				break
		self.updateGL()

	def isinwin(self,x,y):
		for child in self.children:
			if child.isinwin(x,y) : return True
			
		return False
			
class EMFrame(EMWindowNode,EMRegion):
	'''
	EMFrame is a base class for windows that have a frame. The frame defines a 3D bounding box, and is defined
	in terms of its origin and its size in each dimension
	'''
	def __init__(self,parent,geometry=Region(0,0,0,0,0,0)):
		EMWindowNode.__init__(self,parent)
		EMRegion.__init__(self,geometry)
		self.children = []
	
	def draw(self):
		for child in self.children:
			glPushMatrix()
			child.draw()
			glPopMatrix()
	
	def mousePressEvent(self, event):
		#YUCK fixme soon this is terribly inefficient
		for i in self.children:
			i.mousePressEvent(event)
	
	def mouseMoveEvent(self, event):
		#YUCK fixme soon this is terribly inefficient
		for i in self.children:
			i.mouseMoveEvent(event)
		
	def mouseReleaseEvent(self, event):
		#YUCK fixme soon this is terribly inefficient
		for i in self.children:
			i.mouseReleaseEvent(event)

	def mouseDoubleClickEvent(self, event):
		#YUCK fixme soon this is terribly inefficient
		for i in self.children:
			i.mouseDoubleClickEvent(event)

	def wheelEvent(self, event):
		#YUCK fixme soon this is terribly inefficient
		for i in self.children:
			i.wheelEvent(event)
	
	def toolTipEvent(self, event):
		#YUCK fixme soon this is terribly inefficient
		for i in self.children:
			i.toolTipEvent(event)

class EMBoxerFrame(EMFrame,GUIbox):
	def __init__(self,parent,geometry=Region(0,0,0,0,0,0),image_names=[],):
		EMFrame.__init__(self,parent,geometry)
		GUIbox.__alt_init__(self,[],[])
		
	def draw(self):
		pass

class LeftSideWidgetBar(EMGLViewContainer):
	def __init__(self,parent):
		EMGLViewContainer.__init__(self,parent)
		self.mouse_on = None
		self.previous_mouse_on = None
		self.active = None
		self.transformers = []
		
		EMDesktop.main_widget.register_resize_aware(self)
		
	def __del__(self):
		try: EMDesktop.main_widget.deregister_resize_aware(self)
		except: pass # this might happen at program death
	def width(self):
		width = 0
		for child in self.children:
			if child.width() > width:
				width = width
		
		return width
		
	def height(self):
		return viewport_height()

	def draw(self):
		glPushMatrix()
		glTranslate(-self.parent.width()/2.0,self.parent.height()/2.0,0)
		for i,child in enumerate(self.children):
			glPushMatrix()
			self.transformers[i].transform()
			child.draw()
			glPopMatrix()
			glTranslate(0,-self.transformers[i].get_xy_scale()*child.height(),0)

		glPopMatrix()

	def i_initialized(self,child):
		if isinstance(child,EMDesktopTaskWidget):
			pass
			#child.set_cam_pos(-self.parent.width()/2.0+child.width()/2.0,self.parent.height()/2.0-child.height()/2.0,0)

	def add_browser_frame(self):
		browser_frame = EMBrowserFrame(self)
		browser_frame.load_browser()
		self.parent.attach_child(browser_frame)
		#print "done"

	def seed_scale_animation(self,i):
		t = self.transformers[i]
		if t.get_xy_scale() != 1.0:
			#if i == 0:
			seed_height = self.children[i].height()
			below_height = 0
			for j in range(0,len(self.children)):
				if j != i: below_height += self.children[j].height()
			
			
			to_height = viewport_height()-seed_height
			below_scale = to_height/float(below_height)
			
			#print "seed a scale event for ", i
			t.seed_scale_animation_event(1.0)
			for j in range(0,len(self.transformers)):
				if j != i: 
					#print "seed a scale event for ", j
					self.transformers[j].seed_scale_animation_event(below_scale)
			#elif i == (len(self.transformers)-1):
				
		
	def resize_gl(self):
		self.reset_scale_animation()	

	def reset_scale_animation(self):
		children_height = 0
		for child in self.children:
			children_height += child.height()
			
		if children_height > viewport_height():
			scale = viewport_height()/float(children_height)
			#i = 0
			for t in self.transformers: 
				#print "resetting scale event for ", i
				#i += 1
				t.seed_scale_animation_event(scale)
	
	def mouseMoveEvent(self, event):
		intercept = False
		for i,child in enumerate(self.children):
			if ( child.isinwin(event.x(),viewport_height()-event.y()) ):
				#if self.animations[i] == None:
				intercept = True
				self.previous_mouse_on = self.mouse_on
				self.mouse_on = i
				self.active = i
				if self.transformers[i].is_animatable():
					t = self.transformers[i]
					#print "mouse entered a window"
					#print "seed a rotation event for ", i
					t.seed_rotation_animation_event(force_active=True)
					self.seed_scale_animation(i)
				
				if self.previous_mouse_on != None and self.previous_mouse_on != self.mouse_on:
					#print "mouse left a window and left another"
					t = self.transformers[self.previous_mouse_on]
					t.seed_rotation_animation_event(force_inactive=True)
					self.previous_mouse_on = None
				
				break
		
		
		
		if not intercept:
			if self.mouse_on != None:
				#print "moust left a window"
				t = self.transformers[self.mouse_on]
				t.seed_rotation_animation_event(force_inactive=True)
				self.mouse_on = None
				self.previous_mouse_on = None
				self.active = None
				self.reset_scale_animation()
				
		EMGLViewContainer.mouseMoveEvent(self,event)
		
	def attach_child(self,new_child):
		self.transformers.append(LeftSideWidgetBar.LeftSideTransform(new_child))
		EMWindowNode.attach_child(self,new_child)
		self.reset_scale_animation()

	class LeftSideTransform:
		ACTIVE = 0
		INACTIVE = 1
		ANIMATED = 2
		def __init__(self,child):
			self.child = child
			self.rotation_animation = None
			self.scale_animation = None
			self.state = LeftSideWidgetBar.LeftSideTransform.INACTIVE
			self.rotation = 90
			self.xy_scale = 1.0
		
		def get_xy_scale(self):
			return self.xy_scale
		
		def set_xy_scale(self,xy_scale):
			self.xy_scale = xy_scale
		
		def set_rotation(self,rotation):
			self.rotation = rotation
		
		def seed_scale_animation_event(self,scale):
			if self.xy_scale == scale: return
			
			if self.scale_animation != None:
				self.scale_animation.set_animated(False) # this will cause the EMDesktop to stop animating
				self.scale_animation = None
			
			animation = XYScaleAnimation(self,self.xy_scale,scale)
			self.scale_animation = animation
			EMDesktop.main_widget.register_animatable(animation)
			
		def seed_rotation_animation_event(self,force_inactive=False,force_active=False):
			
			if self.state == LeftSideWidgetBar.LeftSideTransform.ACTIVE:
				rot = [0,90]
			elif self.state == LeftSideWidgetBar.LeftSideTransform.ANIMATED:
				c = self.rotation
				s = self.rotation_animation.get_start()
				if c < s: s = 90
				elif c > s: s = 0
				else: print "I'm a bad programmer",c,s
				if force_inactive: s = 90
				if force_active: s = 0
				rot =  [c,s]
			elif self.state == LeftSideWidgetBar.LeftSideTransform.INACTIVE:
				rot = [90,0]
				
			if self.rotation_animation != None:
				self.rotation_animation.set_animated(False) # this will cause the EMDesktop to stop animating
				self.rotation_animation = None

			
			#print "adding animation",rot
			animation = SingleAxisRotationAnimation(self,rot[0],rot[1],[0,1,0])
			self.rotation_animation = animation
			self.state = LeftSideWidgetBar.LeftSideTransform.ANIMATED
			EMDesktop.main_widget.register_animatable(animation)
			
		def is_animatable(self):
			if self.rotation_animation != None: return not self.rotation_animation.is_animated()
			elif self.state == LeftSideWidgetBar.LeftSideTransform.ACTIVE:
				return False
			else: return True

		def transform(self):
			if self.rotation_animation != None and not self.rotation_animation.is_animated():
				end = self.rotation_animation.get_end()
				self.rotation_animation = None
				if end == 0.0:
					self.state = LeftSideWidgetBar.LeftSideTransform.ACTIVE
					self.rotation = 0
				elif end == 90.0:
					self.state = LeftSideWidgetBar.LeftSideTransform.INACTIVE
			
			if self.scale_animation != None and not self.scale_animation.is_animated():
				self.scale_animation.set_animated(False) # this will cause the EMDesktop to stop animating
				self.scale_animation = None
			
			if self.rotation_animation == None:
				if self.rotation != 90 and self.rotation != 0:
					self.rotation = 90
			
			glTranslate(0,-self.xy_scale*self.child.height()/2,0)
			glRotate(self.rotation,0,1,0)
			glTranslate(self.xy_scale*self.child.width()/2.0,0,0)
			glScale(self.xy_scale,self.xy_scale,1.0)
		
		
		def draw(self):
			
			glPushMatrix()
			self.transform()
			self.child.draw()
			glPopMatrix()
			glTranslate(0,-self.xy_scale*self.child.height(),0)
		
class EMBrowserFrame(EMGLViewContainer):
	def __init__(self,parent,geometry=Region(0,0,0,0,0,0)):
		EMGLViewContainer.__init__(self,parent,geometry)
		self.file_dialog = None
		
	def get_file_dialog(self):
		if self.file_dialog == None:
			self.file_dialog = QtGui.QFileDialog(EMDesktop.main_widget,"Open File",QtCore.QDir.currentPath(),QtCore.QString("Image files (*.img *.hed *.mrc *.hdf)"))
			#self.fd = EMDesktopFileDialog(self.parent,"Open File",QtCore.QDir.currentPath(),QtCore.QString("Image files (*.img *.hed *.mrc)"))
			QtCore.QObject.connect(self.file_dialog, QtCore.SIGNAL("finished(int)"), self.finished)
			QtCore.QObject.connect(self.file_dialog, QtCore.SIGNAL("currentChanged(QString)"), self.changed)
			self.file_dialog.show()
			self.file_dialog.hide()
			
		return self.file_dialog

	def changed(self,file):
		pass
	
		#try:
			#a=EMData.read_images(str(file))
			#if len(a) == 1:
				#a = a[0]
				#if a.get_zsize() != 1: w = EMGLView3D(self,a)
				#else: w = EMGLView2D(self,a)
			#else: w = EMGLView2D(self,a)
			
			#try: self.cols[1] = []
			#except: self.cols.append([])
			#self.cols[1].append(w)
			#scalex = self.get_col_width()/float(w.width())
			#scaley = self.get_col_height()/float(w.height())
			##print scalex,scaley,yheight,w.height()
			#if scalex > scaley: scalex = scaley
			#try: w.d = scalex*w.d # 3D
			#except: pass
			#w.h = scalex*w.h
			#w.set_width(scalex*w.w)
			
			#try: w.setOptScale()
			#except: pass
			
			#insp = w.getInspector()
			#d = EMGLViewQtWidget(self.parent)
			#d.setQtWidget(insp)
			#try:
				#self.cols[0][2] = d
			#except: self.cols[0].append(d)
			
			#self.layout()
	def emit(self,*args, **kargs):
		EMDesktop.main_widget.emit(*args,**kargs)
	
	def get_near_plane_dims(self):
		return EMDesktop.main_widget.get_near_plane_dims()
	
	def getStartZ(self):
		return EMDesktop.main_widget.getStartZ()
	
	def finished(self,val):
		if ( val == 1 ):
			for i in self.file_dialog.selectedFiles():
				a=EMData.read_images(str(i))
				if len(a) == 1:
					a = a[0]
					if a.get_zsize() != 1:
						w = EMGLView3D(self,a)
						self.attach_child(w)
					else:
						w = EMGLView2D(self,a)
						self.attach_child(w)
				else:
					w = EMGLView2D(self,a)
					self.attach_child(w)
					

	def load_browser(self):
		#self.numcols = 2
		file_dialog = self.get_file_dialog()
		self.fd_widget = EMGLViewQtWidget(EMDesktop.main_widget)
		self.fd_widget.setQtWidget(file_dialog)
		self.parent.attach_child(self.fd_widget)
	
class EMDesktopFrame(EMFrame):
	def __init__(self,parent,geometry=Region(0,0,0,0,0,0)):
		EMFrame.__init__(self,parent,geometry)
		
	def set_geometry(self,geometry):
		EMFrame.set_geometry(self,geometry)
		try:
			for child in self.children:
				if isinstance(child,EMDesktopTaskWidget):
					child.set_cam_pos(-self.parent.width()/2.0+child.width()/2.0,self.parent.height()/2.0-child.height()/2.0,0)
					
		except: pass
		

	def i_initialized(self,child):
		if isinstance(child,EMDesktopTaskWidget):
			child.set_cam_pos(-self.parent.width()/2.0+child.width()/2.0,self.parent.height()/2.0-child.height()/2.0,0)

	def add_browser_frame(self):
		browser_frame = EMBrowserFrame(self,self.geometry,self.parent)
		browser_frame.load_browser()
		self.attach_child(browser_frame)
		print "done"
	
	def updateGL(self):
		self.parent.updateGL()
	
class EMDesktopScreenInfo:
	"""
	A class the figures out how many screen the user has, whether or they are running a virtual desktop etc.
	 
	"""
	
	def __init__(self):
		app=QtGui.QApplication.instance()
		sys_desktop = app.desktop()
		self.__num_screens = sys_desktop.numScreens()
		print "there are this many screens", self.__num_screens
		self.__screens = [] # This will be a list of QtCore.QGeometry objects
		
		print "there are",self.__num_screens,"screen is and the primary screen is", sys_desktop.primaryScreen()
		#if sys_desktop.isVirtualDesktop() or True:
			#print "user is running a virtual desktop"
			##print "now trying to figure out the dimensions of each desktop"
		if self.__num_screens == 1:
			app_screen = sys_desktop.screen(0)
			self.__screens.append(sys_desktop.screenGeometry(app_screen))
			print self.__screens
			print "there is only one screen its dimensions are",app_screen.width(),app_screen.height()
		else:
			x=0
			y=0
			for i in range( self.__num_screens):				
				geom = self.__screens.append(sys_desktop.availableGeometry(QtCore.QPoint(x,y)))
				print "\t  printing available geometry information it is",geom.left(),geom.right(),geom.top(),geom.bottom()
				print "\t geometry starts at",geom.left(),geom.top()," and has dimensions", geom.right()-geom.left()+1,geom.bottom()-geom.top()+1
				x = geom.right()+1
				y = geom.top()
		#else:
			#print "non virtual desktops are not yet supported"

	def get_num_screens():
		return self.__num_screens
	
	def get_screens():
		return self.__screens

class EMDesktop(QtOpenGL.QGLWidget):
	main_widget = None
	"""An OpenGL windowing system, which can contain other EMAN2 widgets and 3-D objects.
	"""
	def __init__(self):
		EMDesktop.main_widget = self
		frame = EMBoxerFrame(self)
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True)
		fmt.setSampleBuffers(True)
		QtOpenGL.QGLWidget.__init__(self,fmt,None)
		
		self.gq=0			# quadric object for cylinders, etc	
		self.app=QtGui.QApplication.instance()
		self.sysdesktop=self.app.desktop()
		self.appscreen=self.sysdesktop.screen(self.sysdesktop.primaryScreen())
		self.frame_dl = 0 # display list of the desktop frame
		self.fov = 40
		self.resize_aware_objects = []
		self.animatables = []
		
		# what is this?
		self.bgob2=ob2dimage(self,self.read_EMAN2_image())
		
		
		self.setMouseTracking(True)
		
		# this float widget has half of the screen (the left)
		self.desktop_frame = EMDesktopFrame(self)
		self.left_side_bar = LeftSideWidgetBar(self.desktop_frame)
		fw1 = EMDesktopTaskWidget(self.left_side_bar)
		self.left_side_bar.attach_child(fw1)
		self.desktop_frame.attach_child(self.left_side_bar)
		
		#print fw1.width(),fw1.height()
		
		fw1.set_gl_widget(self)
		
		
		
		
	
		self.glbasicobjects = EMBasicOpenGLObjects()
		self.borderwidth=10.0
		self.cam = Camera()
		
		self.frame_dl = 0
		
		self.timer = QTimer()
		QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.time_out)
		self.timer.start(10)
		self.time=0
		self.begin_time = time()
		
		# resize finally so that the full screen is used
		self.show()
		self.move(0,0)
		self.resize(self.appscreen.size())
		
	#def get_app_width(self):
		#return self.appwidth
	
	#def get_app_height(self):
		#return self.appheight
		
	def register_resize_aware(self,resize_aware_object):
		self.resize_aware_objects.append(resize_aware_object)
		
	def deregister_resize_aware(self,resize_aware_object):
		for i,obj in enumerate(self.resize_aware_objects):
			if obj == resize_aware_object:
				self.resize_aware_objects.pop(i)
				return
			
		print "warning, can't deregestire resize aware object",resize_aware_object
	
	def register_animatable(self,animatable):
		self.animatables.append(animatable)
	def get_aspect(self):
		return float(self.width())/float(self.height())
	
	def get_z_opt(self):
		return (1.0/tan(self.get_fov()/2.0*pi/180.0))*self.height()/2
	
	def get_fov(self):
		return self.fov
	
	def get_depth_for_height(self, height):
		return 0
		# This function returns the width and height of the renderable 
		# area at the origin of the data volume
		depth = height/(2.0*tan(self.fov/2.0*pi/180.0))
		return depth
	
	def get_render_dims_at_depth(self, depth):
		return 0
		# This function returns the width and height of the renderable 
		# area at the origin of the data volume
		height = -2*tan(self.fov/2.0*pi/180.0)*(depth)
		width = self.aspect*height
		return [width,height]
	
	def draw_frame(self):
		if self.frame_dl == 0:
			#print self.appwidth/2.0,self.appheight/2.0,self.zopt
			glCallList(self.glbasicobjects.getCylinderDL())
			length = self.get_z_opt()
			self.frame_dl=glGenLists(1)
			glNewList(self.frame_dl,GL_COMPILE)
			glPushMatrix()
			glTranslatef(-self.width()/2.0,-self.height()/2.0,0.0)
			glScaled(self.borderwidth,self.borderwidth,length)
			glCallList(self.glbasicobjects.getCylinderDL())
			glPopMatrix()
			glPushMatrix()
			glTranslatef( self.width()/2.0,-self.height()/2.0,0.0)
			glScaled(self.borderwidth,self.borderwidth,length)
			glCallList(self.glbasicobjects.getCylinderDL())
			glPopMatrix()
			
			glPushMatrix()
			glTranslatef( self.width()/2.0, self.height()/2.0,0.0)
			glScaled(self.borderwidth,self.borderwidth,length)
			glCallList(self.glbasicobjects.getCylinderDL())
			glPopMatrix()
			
			glPushMatrix()
			glTranslatef(-self.width()/2.0, self.height()/2.0,0.0)
			glScaled(self.borderwidth,self.borderwidth,length)
			glCallList(self.glbasicobjects.getCylinderDL())
			glPopMatrix()
			
			glEndList()
			
		if self.frame_dl == 0:
			print "error, frame display list failed to compile"
			exit(1)
		glColor(.9,.2,.8)
		## this is a nice light blue color (when lighting is on)
		## and is the default color of the frame
		glMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,(.2,.2,.8,1.0))
		glMaterial(GL_FRONT,GL_SPECULAR,(.8,.8,.8,1.0))
		glMaterial(GL_FRONT,GL_SHININESS,128.0)
		glCallList(self.frame_dl)
		
	def read_EMAN2_image(self):
		self.p = QtGui.QPixmap("EMAN2.0.big.jpg")
		return self.p

	
	def get_time(self):
		return self.time
	
	def time_out(self):
		self.time  = time()-self.begin_time
		rm = []
		if len(self.animatables) != 0:
			#print len(self.animatables)
			
			for a,i in enumerate(self.animatables):
				if not i.animate(self.time): rm.append(a)
				
			rm.reverse()
			for a in rm:
				self.animatables.pop(a)
			
		
			self.updateGL()
		
		
	def initializeGL(self):
		glClearColor(0,0,0,0)
		glEnable(GL_NORMALIZE)
				
		glEnable(GL_LIGHTING)
		glEnable(GL_LIGHT0)
		glEnable(GL_DEPTH_TEST)
		glLightfv(GL_LIGHT0, GL_AMBIENT, [0.9, 0.9, 0.9, 1.0])
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_POSITION, [0.,0,1,0.])

		# get a new Quadric object for drawing cylinders, spheres, etc
		if not self.gq:
			self.gq=gluNewQuadric()
			gluQuadricDrawStyle(self.gq,GLU_FILL)
			gluQuadricNormals(self.gq,GLU_SMOOTH)
			gluQuadricOrientation(self.gq,GLU_OUTSIDE)
			gluQuadricTexture(self.gq,GL_FALSE)
	
		glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST)
		glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST)
		glHint(GL_TEXTURE_COMPRESSION_HINT, GL_NICEST)
	
	
	#def gl_error(self):
		#OpenGL.error.ErrorChecker.glCheckError(
		
	
	def paintGL(self):
		#OpenGL.error.ErrorChecker.registerChecker( self.gl_error )
		#print self.glCheckError()
		
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		#print "dims are ", self.appwidth,self.appheight,self.width(),self.height()

		
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		#self.bgob.render()
		glPushMatrix()
		if (self.get_time() < 0):
			z = self.get_z_opt() + float(self.get_time())/2.0*self.get_z_opt()
			#print z
			glTranslatef(0.,0.,-z)
		else:
			#print -2*self.zopt+0.1
			glTranslatef(0.,0.,-2*self.get_z_opt()+0.1)
			
		glPushMatrix()
		self.draw_frame()
		glPopMatrix()
		

		glPushMatrix()
		glScalef(self.height()/2.0,self.height()/2.0,1.0)
		self.bgob2.render()
		glPopMatrix()
		glPopMatrix()
		
		if self.get_time()>2 or True :
			dx = (self.width() - self.p.width())/2.0
			dy = (self.height() - self.p.height())/2.0
			
			
			#glPushMatrix()
			#glTranslatef(self.width()/2.0+self.p.width(),-dy, -1.8*self.get_z_opt())
			#glScalef(.25,.25,.25)
			#glRotate(self.get_time()*100.0,1.0,0.,0.0)
			#self.bgob2.render2()
			#glPopMatrix()

			#glPushMatrix()
			#dx = (self.width() - self.p.width())/2.0
			#dy = (self.height() - self.p.height())/2.0
			#dz = (self.get_time()%20)/10.0
			#if ( dz > 1): dz = 2-dz
			#glTranslatef(self.width()/2.0,-dy, -self.get_z_opt()-self.p.width()/2.0-dz*(self.get_z_opt()-self.p.width()))
			#glRotatef(-90,0,1,0)
			#self.bgob2.render2()
			#glPopMatrix()

			glPushMatrix()
			glTranslatef(0.,0.,-self.get_z_opt())
			self.desktop_frame.draw()
			glPopMatrix()

	def resizeGL(self, width, height):
		side = min(width, height)
		glViewport(0,0,self.width(),self.height())
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		
		self.zNear = self.get_z_opt()
		self.zFar = 2*self.get_z_opt()
		gluPerspective(self.fov,self.get_aspect(),self.zNear-500,self.zFar)
		
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		if self.frame_dl != 0:
			glDeleteLists(self.frame_dl,1)
			self.frame_dl = 0
		
		self.desktop_frame.set_geometry(Region(0,0,self.width(),self.height()))
		resize_gl()
		
		for obj in self.resize_aware_objects:
			obj.resize_gl()
		
	def mouseReleaseEvent(self, event):
		if event.button()==Qt.LeftButton:
			pass
		
	def mousePressEvent(self, event):
		self.desktop_frame.mousePressEvent(event)
	
	def mouseMoveEvent(self, event):
		#YUCK fixme soon
		self.desktop_frame.mouseMoveEvent(event)
		
	def mouseReleaseEvent(self, event):
		#YUCK fixme soon
		self.desktop_frame.mouseReleaseEvent(event)

	def mouseDoubleClickEvent(self, event):
		#YUCK fixme soon
		self.desktop_frame.mouseDoubleClickEvent(event)

	def wheelEvent(self, event):
		#YUCK fixme soon
		self.desktop_frame.wheelEvent(event)

	def toolTipEvent(self, event):
		#YUCK fixme soon
		self.desktop_frame.toolTipEvent(event)
		QtGui.QToolTip.hideText()
		
	def dragMoveEvent(self,event):
		print "received drag move event, but I don't do anything about it :("
		
	def event(self,event):
		if event.type() == QtCore.QEvent.MouseButtonPress: 
			self.mousePressEvent(event)
			return True
		elif event.type() == QtCore.QEvent.MouseButtonRelease:
			self.mouseReleaseEvent(event)
			return True
		elif event.type() == QtCore.QEvent.MouseMove: 
			self.mouseMoveEvent(event)
			return True
		elif event.type() == QtCore.QEvent.MouseButtonDblClick: 
			self.mouseDoubleClickEvent(event)
			return True
		elif event.type() == QtCore.QEvent.Wheel: 
			self.wheelEvent(event)
			return True
		elif event.type() == QtCore.QEvent.ToolTip: 
			self.toolTipEvent(event)
			return True
		else: 
			return QtOpenGL.QGLWidget.event(self,event)

	def hoverEvent(self,event):
		#YUCK fixme soon
		for i in self.floatwidgets:
			i.hoverEvent(event)

	def get_near_plane_dims(self):
		height = 2.0*self.zNear * tan(self.fov/2.0*pi/180.0)
		width = self.get_aspect() * height
		return [width,height]
		
	def getStartZ(self):
		return self.zNear

class EMDesktopTaskWidget(EMGLViewContainer):
	def __init__(self, parent):
		#print "init"
		EMGLViewContainer.__init__(self,parent)
		self.parent = parent
	
		self.init_flag = True
		
		self.desktop_task_widget = None

		self.cam = Camera()
		
		self.glwidget = None
		self.widget = None
		
		self.animation = None
	
	def register_animation(self,animation):
		self.animation = animation
	
	def set_gl_widget(self,widget):
		self.gl_widget = widget

	def get_depth_for_height(self, height):
		try: 
			return self.gl_widget.get_depth_for_height(height)
		except:
			print "parent can't get height for depth"
			exit(1)
			#return 0
			
	def set_cam_pos(self, x,y,z=0):
		self.cam.cam_x = x
		self.cam.cam_y = y
		self.cam.cam_z = z

	def height(self):
		if self.widget != None: return self.widget.height() + 2*self.desktop_task_widget.get_decoration().get_border_width()
		return 0
		
	def width(self):
		if self.widget != None: return self.widget.width() + 2*self.desktop_task_widget.get_decoration().get_border_height()
		return 0
		
	def updateGL(self):
		try: self.parent.updateGL()
		except: pass
			
	def draw(self):
		if ( self.init_flag == True ):
			self.desktop_task_widget = EMGLViewQtWidget(self.gl_widget)
			self.widget = EMDesktopTaskWidget.EMDesktopTaskInspector(self)
			self.widget.show()
			self.widget.hide()
			self.widget.resize(150,150)
			self.desktop_task_widget.setQtWidget(self.widget)
			self.init_flag = False
			self.attach_child(self.desktop_task_widget)
			#self.parent.i_initialized(self)
		
		if self.animation != None:
			self.animation.transform()
		
		self.cam.position()
		
		for child in self.children:
			glPushMatrix()
			child.draw()
			glPopMatrix()
			
	def bindTexture(self,pixmap):
		return self.gl_widget.bindTexture(pixmap)
	
	def deleteTexture(self,val):
		return self.gl_widget.deleteTexture(val)
	
	def get_render_dims_at_depth(self, depth):
		try: return self.gl_widget.get_render_dims_at_depth(depth)
		except:
			print "parent can't get render dims at for depth"
			return

	def close(self):
		pass

	def add_browser(self):
		self.parent.add_browser_frame()
		print "adding browser"
		
	def add_boxer(self):
		print "adding boxer"

	class EMDesktopTaskInspector(QtGui.QWidget):
		def __init__(self,target) :
			QtGui.QWidget.__init__(self,None)
			self.target=target
			
			
			self.vbl = QtGui.QVBoxLayout(self)
			self.vbl.setMargin(0)
			self.vbl.setSpacing(6)
			self.vbl.setObjectName("vbl")
			
			self.hbl_buttons2 = QtGui.QHBoxLayout()
			
			self.tree_widget = QtGui.QTreeWidget(self)
			self.tree_widget_entries = []
			self.tree_widget_entries.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Browse")))
			self.tree_widget_entries.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Box")))
			self.tree_widget.insertTopLevelItems(0,self.tree_widget_entries)
			self.tree_widget.setHeaderLabel("Choose a task")
			
			self.hbl_buttons2.addWidget(self.tree_widget)
			
			self.close = QtGui.QPushButton("Close")
			
			self.vbl.addLayout(self.hbl_buttons2)
			self.vbl.addWidget(self.close)
			
			QtCore.QObject.connect(self.tree_widget, QtCore.SIGNAL("itemDoubleClicked(QTreeWidgetItem*,int)"), self.tree_widget_double_click)
			QtCore.QObject.connect(self.close, QtCore.SIGNAL("clicked()"), self.target.close)
			
		def tree_widget_double_click(self,tree_item,i):
			task = tree_item.text(0)
			if task == "Browse":
				self.target.add_browser()
			elif task == "Box":
				self.target.add_boxer()
		
class EMDesktopTaskWidget3:
	""" Something that organizes and display EMGLWidgets
	"""
	def __init__(self, parent=None):
		#print "init"
		self.parent = parent
	
		self.current = None
		self.previous = None
	
		self.fd = None
	
		self.init_flag = True
		self.qwidgets = []
		self.imagewidgets = []
		self.inspectorwidgets = []
		self.cols = []
		
		self.desktopwidget = None
		self.suppressUpdateGL = False
		
		self.fdxpos = 0
		self.fdypos = 0

		self.cam = Camera()
		
		self.glwidget = None
		
	def set_gl_widget(self,widget):
		self.gl_widget = widget

	def get_depth_for_height(self, height):
		try: 
			return self.gl_widget.get_depth_for_height(height)
		except:
			print "parent can't get height for depth"
			exit(1)
			#return 0

		
	def set_cam_pos(self, x,y,z=0):
		self.cam.cam_x = x
		self.cam.cam_y = y
		self.cam.cam_z = z

	def height(self):
		return self.parent.height()
	
	def width(self):
		return self.parent.width()

	def updateGL(self):
		if not self.suppressUpdateGL:
			try: self.parent.updateGL()
			except: pass

	def get_col_width(self):
		numcols = len(self.cols)
		if numcols == 0: return 0
		cw = float(self.parent.width())/float(numcols)
		return cw
		
	def get_col_height(self):
		return float(self.parent.height())
	
	def layout(self):
		print "in layout"
		numcols = len(self.cols)
		if numcols == 0: return
		
		cw = self.get_col_width()
		
		colidx = 0
		for i in self.cols:
			ch = 0
			n = len(i)
			if n == 0:
				continue
			for j in i:
				ch += j.height()
				
				if j.width() > cw:
					print "error, can handle wide widgets",j.width(),cw
			
			dif = float(self.parent.height() - ch)
			if dif < 0:
				print "error, can't handle negative diffs atm", dif, colidx, n
				dif = 0
				##exit(1)
			
			
			dif /= float(n)
			dy = dif/2.0
			
			dx = cw/2.0 + colidx*cw
			for j in i:
				dy += j.height()/2.0
				j.cam.cam_y = dy
				j.cam.cam_x = dx
				#print j.cam.cam_y,j.cam.cam_x
				dy += j.height()/2.0
				dy += dif
			
			
			colidx += 1
		
	def draw(self):
		self.cam.position()
		if ( self.init_flag == True ):
			self.desktopwidget = EMGLViewQtWidget(self.gl_widget)
			self.s = EMDesktopInspector(self)
			self.s.show()
			self.s.hide()
			#self.desktopwidget.cam.cam_x = s.width()/2.0
			#self.desktopwidget.cam.cam_y = self.h - s.height()/2.0
			self.desktopwidget.setQtWidget(self.s)
			self.qwidgets.append(self.desktopwidget)
			self.numcols = 1
			self.cols = [[self.desktopwidget]]
			#self.layout()
			self.init_flag = False
		
		#glPushMatrix()
		#self.desktopwidget.paintGL()
		#glPopMatrix()
		for i in self.cols:
			for j in i:
				glPushMatrix()
				j.paintGL()
				glPopMatrix()
		
	def changed(self,file):
	
		#try:
			a=EMData.read_images(str(file))
			if len(a) == 1:
				a = a[0]
				if a.get_zsize() != 1: w = EMGLView3D(self,a)
				else: w = EMGLView2D(self,a)
			else: w = EMGLView2D(self,a)
			
			try: self.cols[1] = []
			except: self.cols.append([])
			self.cols[1].append(w)
			scalex = self.get_col_width()/float(w.width())
			scaley = self.get_col_height()/float(w.height())
			#print scalex,scaley,yheight,w.height()
			if scalex > scaley: scalex = scaley
			try: w.d = scalex*w.d # 3D
			except: pass
			w.h = scalex*w.h
			w.set_width(scalex*w.w)
			
			try: w.setOptScale()
			except: pass
			
			insp = w.getInspector()
			d = EMGLViewQtWidget(self.parent)
			d.setQtWidget(insp)
			try:
				self.cols[0][2] = d
			except: self.cols[0].append(d)
			
			self.layout()
			
		#except:
			#print "could not open"
			#return
		
		
		#print "paintGL done"
	def finished(self,val):
		pass
		if ( val == 1 ):
			for i in self.fd.selectedFiles():
				a=EMData.read_images(str(i))
				if len(a) == 1:
					a = a[0]
					if a.get_zsize() != 1:
						w = EMGLView3D(self,a)
						self.qwidgets.append(w)
					else:
						w = EMGLView2D(self,a)
						self.qwidgets.append(w)
				else:
					w = EMGLView2D(self,a)
					self.qwidgets.append(w)
					
	def timer(self):
		pass
		#self.updateGL()
		
	def bindTexture(self,pixmap):
		return self.gl_widget.bindTexture(pixmap)
	
	def deleteTexture(self,val):
		return self.gl_widget.deleteTexture(val)
	
	def get_render_dims_at_depth(self, depth):
		try: return self.gl_widget.get_render_dims_at_depth(depth)
		except:
			print "parent can't get render dims at for depth"
			return

	def resizeEvent(self, width, height):
		widgets = self.qwidgets + self.imagewidgets
		for i in widgets:
			i.set_update_P_inv()
			
	def getWidgets(self):
		widgets = []
		for i in self.cols:
			for j in i:
				widgets.append(j)
		
		return widgets
	
	def mousePressEvent(self, event):
		widgets = self.getWidgets()
		for i in widgets:
			if ( i.isinwin(event.x(),self.height()-event.y()) ):
				i.mousePressEvent(event)
				intercepted = True
				self.updateGL()
				return
	
	def mouseMoveEvent(self, event):
		widgets = self.getWidgets()
		for i in widgets:
			if ( i.isinwin(event.x(),self.height()-event.y()) ):
				self.current = i
				if (self.current != self.previous ):
					if ( self.previous != None ):
						self.previous.leaveEvent()
				i.mouseMoveEvent(event)
				self.previous = i
				self.updateGL()
				return
		
	def mouseReleaseEvent(self, event):
		widgets = self.getWidgets()
		for i in widgets:
			if ( i.isinwin(event.x(),self.height()-event.y()) ):
				i.mouseReleaseEvent(event)
				self.updateGL()
				return
					
		
	def mouseDoubleClickEvent(self, event):
		widgets = self.getWidgets()
		for i in widgets:
			if ( i.isinwin(event.x(),self.height()-event.y()) ):
				i.mouseDoubleClickEvent(event)
				self.updateGL()
				return
		
		
	def wheelEvent(self, event):
		widgets = self.getWidgets()
		for i in widgets:
			if ( i.isinwin(event.x(),self.height()-event.y()) ):
				i.wheelEvent(event)
				self.updateGL()
				return

	def toolTipEvent(self, event):
		widgets = self.getWidgets()
		for i in widgets:
			if ( i.isinwin(event.x(),self.height()-event.y()) ):
				i.toolTipEvent(event)
				self.updateGL()
				return
		
		QtGui.QToolTip.hideText()

	def dragMoveEvent(self,event):
		print "received drag move event"
		
	def event(self,event):
		#print "event"
		#QtGui.QToolTip.hideText()
		if event.type() == QtCore.QEvent.MouseButtonPress: 
			self.mousePressEvent(event)
			return True
		elif event.type() == QtCore.QEvent.MouseButtonRelease:
			self.mouseReleaseEvent(event)
			return True
		elif event.type() == QtCore.QEvent.MouseMove: 
			self.mouseMoveEvent(event)
			return True
		elif event.type() == QtCore.QEvent.MouseButtonDblClick: 
			self.mouseDoubleClickEvent(event)
			return True
		elif event.type() == QtCore.QEvent.Wheel: 
			self.wheelEvent(event)
			return True
		elif event.type() == QtCore.QEvent.ToolTip: 
			self.toolTipEvent(event)
			return True
		else: 
			return QtOpenGL.QGLWidget.event(self,event)

	def hoverEvent(self,event):
		#print "hoverEvent"
		if self.inspector :
			widgets = self.qwidgets + self.imagewidgets
			for i in widgets:
				if ( i.isinwin(event.x(),self.height()-event.y()) ):
					i.hoverEvent(event)
					break
		self.updateGL()

	def add_boxer(self):
		pass

	def getfd(self):
		if self.fd == None:
			self.fd = QtGui.QFileDialog(self.gl_widget,"Open File",QtCore.QDir.currentPath(),QtCore.QString("Image files (*.img *.hed *.mrc *.hdf)"))
			#self.fd = EMDesktopFileDialog(self.parent,"Open File",QtCore.QDir.currentPath(),QtCore.QString("Image files (*.img *.hed *.mrc)"))
			QtCore.QObject.connect(self.fd, QtCore.SIGNAL("finished(int)"), self.finished)
			QtCore.QObject.connect(self.fd, QtCore.SIGNAL("currentChanged(QString)"), self.changed)
			self.fd.show()
			self.fd.hide()
			
		return self.fd

	def add_browser(self):
		#self.numcols = 2
		fd = self.getfd()
		self.fdwidget = EMGLViewQtWidget(self.gl_widget)
		#self.fdwidget.cam.cam_x = -(self.parent.get_app_width() - fd.width())/2.0
		#self.fdwidget.cam.cam_y = (self.parent.get_app_height() - fd.height())/2.0
		self.fdwidget.setQtWidget(fd)
		self.cols = [[self.desktopwidget, self.fdwidget ]]
		self.layout()
		
	
	def addCompare(self):
		pass
	
	def addFSC(self):
		pass
	
	def close(self):
		self.cols = [[self.desktopwidget]]
		self.layout()
	
	def getNearPlaneDims(self):
		return self.parent.getNearPlaneDims()
		
	def get_start_z(self):
		return self.parent.get_start_z()
		
class ob2dimage:
	def __init__(self,target,pixmap):
		self.pixmap=pixmap
		self.target=target
		self.target.makeCurrent()
		self.itex=self.target.bindTexture(self.pixmap)

	def __del__(self):
		target.deleteTexture(self.itex)

	def setWidget(self,widget,region=None):
		return

	def update(self):
		return
	
	def width(self):
		return self.pixmap.width()
	
	def height(self):
		return self.pixmap.height()
	
	def asp(self):
		return (1.0*self.width())/self.height()
	
	def render2(self):
		if not self.pixmap : return
		glPushMatrix()
		glScalef(self.pixmap.width()/2.0,self.pixmap.height()/2.0,1.0)
		glColor(1.0,1.0,1.0)
		glEnable(GL_TEXTURE_2D)
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE)
		glBindTexture(GL_TEXTURE_2D,self.itex)
		glBegin(GL_QUADS)
		glTexCoord2f(0.,0.)
		glVertex(-1.0,-1.0)
		glTexCoord2f(.999,0.)
		glVertex( 1.0,-1.0)
		glTexCoord2f(.999,0.999)
		glVertex( 1.0, 1.0)
		glTexCoord2f(0.,.999)
		glVertex(-1.0, 1.0)
		glEnd()
		glPopMatrix()
		glDisable(GL_TEXTURE_2D)
	
	def render(self):
		if not self.pixmap : return
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP)
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP)
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE)
		glColor(1.0,1.0,1.0)
		glEnable(GL_TEXTURE_2D)
		glBindTexture(GL_TEXTURE_2D,self.itex)
		glBegin(GL_QUADS)
		glTexCoord2f(0.,0.)
		glVertex(-self.target.get_aspect(),-1.0)
		glTexCoord2f(.999,0.)
		glVertex( self.target.get_aspect(),-1.0)
		glTexCoord2f(.999,0.999)
		glVertex( self.target.get_aspect(), 1.0)
		glTexCoord2f(0.,.999)
		glVertex(-self.target.get_aspect(), 1.0)
		glEnd()
		glDisable(GL_TEXTURE_2D)
		#glPopMatrix()


			
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = EMDesktop()
#	window.showFullScreen()
	window.app.exec_()
