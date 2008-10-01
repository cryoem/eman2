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
		
	def detach_child(self,new_child):
		for i,child in enumerate(self.children):
			if (child == new_child):
				self.children.pop(i)
				return
			
		print "error, can not detach the child from a parent it doesn't belong to"
	
	def set_parent(self,parent):
		self.parent = parent
		
	def parent_width(self):
		return parent.width()
	
	def parent_height(self):
		return parent.height()
	
	def get_children(self):
		return self.children

	def emit(self,*args, **kargs):
		EMDesktop.main_widget.emit(*args,**kargs)
		
	def get_near_plane_dims(self):
		return EMDesktop.main_widget.get_near_plane_dims()
	
	def getStartZ(self):
		return EMDesktop.main_widget.getStartZ()
	
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
		self.time_interval = 0.4 # 0.5 seconds for the animation to complete
		self.inverse_time_inverval = 1.0/self.time_interval
		self.time_begin = 0 # records the time at which the animation was begun
		self.animated = True
		self.n = 100
		if Animatable.cache_dts == None:
			self.init_cache_dts()
		
	def init_cache_dts(self):
		from math import tanh
		Animatable.cache_dts = []
		for i in range(self.n):
			tanh_approach = True
			if tanh_approach:
				val = (1+ (tanh(-4+float(i)/(self.n-1)*8)))/2.0
				Animatable.cache_dts.append(val)
			else:
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
			glTranslate(*self.get_origin())
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
				return True
		
		False
	
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
				return True
		
		return False
		
	def mouseReleaseEvent(self, event):
		for child in self.children:
			if ( child.isinwin(event.x(),viewport_height()-event.y()) ):
				child.mouseReleaseEvent(event)
				self.updateGL()
				return True
			
		return False
					
		
	def mouseDoubleClickEvent(self, event):
		for child in self.children:
			if ( child.isinwin(event.x(),viewport_height()-event.y()) ):
				child.mouseDoubleClickEvent(event)
				self.updateGL()
				return True
		return False
		
	def wheelEvent(self, event):
		for child in self.children:
			if ( child.isinwin(event.x(),viewport_height()-event.y()) ):
				child.wheelEvent(event)
				self.updateGL()
				return True
		
		return False

	def toolTipEvent(self, event):
		for child in self.children:
			if ( child.isinwin(event.x(),viewport_height()-event.y()) ):
				child.toolTipEvent(event)
				self.updateGL()
				QtGui.QToolTip.hideText()
				return True
		
		return False

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
				return True
		
		return False
				#break
		#self.updateGL()

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
		
		self.browser_module = None
		
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
		if self.browser_module == None:
			self.browser_module = EMBrowserModule(self,self.parent.get_display_child(0))
			EMDesktop.main_widget.attach_module(self.browser_module)
		
		self.browser_module.load_browser()

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
			
			if below_scale > 1.0: below_scale = 1.0
			
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
		else: scale = 1.0
		
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
		#print_node_hierarchy(self.parent)
		
	def detach_child(self,new_child):
		for i,child in enumerate(self.children):
			if (child == new_child):
				self.children.pop(i)
				self.transformers.pop(i)
				#self.reset_scale_animation()
				return
			
		print "error, attempt to detach a child that didn't belong to this parent"

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
			
		def __del__(self):
			if self.scale_animation != None:
				self.scale_animation.set_animated(False) # this will cause the EMDesktop to stop animating
				self.scale_animation = None
			
			if self.rotation_animation != None:
				self.rotation_animation.set_animated(False) # this will cause the EMDesktop to stop animating
				self.rotation_animation = None
		
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
		
		
class EMModule:
	def __init__(self):
		pass
	
class EMBrowserModule(EMModule):
	def __init__(self,inspector_target,display_target):
		EMModule.__init__(self)
		
		self.inspector_target = inspector_target
		self.display_target = display_target
		
		self.file_dialogs = []
		self.inspector = None
		self.image_display = None
		self.main_tab = None
		
	def changed(self,val):
		print "changed"
		
	def finished(self,val):
		#print "finished"
		act = False
		if ( val == 1 ):
			if self.image_display != None:
				self.display_target.detach_child(self.image_display)
				self.image_display = None
			if self.inspector != None:
				self.inspector_target.detach_child(self.inspector)
				self.inspector = None
			
			if len(self.file_dialogs) == 1:
				file_dialog = self.file_dialogs[0]
			else:
				file_dialog = self.file_dialogs[self.main_tab.currentIndex()]
			
			for i in file_dialog.selectedFiles():
				a=EMData.read_images(str(i))
				if len(a) == 1:
					a = a[0]
					if a.get_zsize() != 1:
						b = EMGLView3D(self.display_target,a)
						
						self.image_display = EM3DWidget(self,b)
						self.image_display.target_translations_allowed(True)
						self.image_display.allow_camera_rotations(True)
						#self.image_display = EMGLView3D(self.display_target,a)
						act = True
					else:
						self.image_display = EMGLView2D(self.display_target,a)
						act = True
				else:
					self.image_display = EMGLView2D(self.display_target,a)
					act = True
					
				if act: break # only read one image at a time
					
		if act:
			self.display_target.attach_child(self.image_display)
			insp = self.image_display.get_inspector()
			self.inspector = EMGLViewQtWidget(EMDesktop.main_widget)
			self.inspector.setQtWidget(insp)
			self.inspector_target.attach_child(self.inspector)
			if len(self.file_dialogs) != 1:
				self.main_tab.setTabText(self.main_tab.currentIndex(),QtCore.QString(i))
			#file_dialog.show()
			#file_dialog.hide()
			
		if self.main_tab != None: self.main_tab.resize(600,400)
		
	def load_browser(self):
		act = False
		
		file_dialog = self.__get_file_dialog()
		
		if len(self.file_dialogs) == 1:
			self.fd_widget = EMGLViewQtWidget(EMDesktop.main_widget)
			self.fd_widget.setQtWidget(file_dialog)
			file_dialog.resize(600,400)
			self.inspector_target.attach_child(self.fd_widget)
		else:
			if self.main_tab == None:
				self.inspector_target.detach_child(self.fd_widget)
				self.main_tab = QtGui.QTabWidget(EMDesktop.main_widget)
				self.main_tab.addTab(self.file_dialogs[0],"Browse")
				self.fd_widget = EMGLViewQtWidget(EMDesktop.main_widget)
				self.fd_widget.setQtWidget(self.main_tab)
				self.inspector_target.attach_child(self.fd_widget)
				
			#for i in self.file_dialogs:
				#i.setParent(EMDesktop.Widget)
			
			self.main_tab.addTab(file_dialog,"Browse")
			#self.main_tab.show()
			#self.main_tab.hide()
			self.main_tab.resize(600,400)
			self.fd_widget.updateTexture(force=True)
			
			#for i in self.file_dialogs:
				#i.setParent(EMDesktop.main_widget)
			
		EMDesktop.main_widget.updateGL()
	
	def accepted(self):
		print "accepted"
		
	def rejected(self):
		print "rejected"
	
	def files_selected(self,strin):
		print strin
	
	def __get_file_dialog(self):
		file_dialog = QtGui.QFileDialog(EMDesktop.main_widget,"Open File",QtCore.QDir.currentPath(),QtCore.QString("Image files (*.img *.hed *.mrc *.hdf)"))
		#self.fd = EMDesktopFileDialog(self.parent,"Open File",QtCore.QDir.currentPath(),QtCore.QString("Image files (*.img *.hed *.mrc)"))
		QtCore.QObject.connect(file_dialog, QtCore.SIGNAL("finished(int)"), self.finished)
		QtCore.QObject.connect(file_dialog, QtCore.SIGNAL("filesSelected(QString)"), self.files_selected)
		QtCore.QObject.connect(file_dialog, QtCore.SIGNAL("accepted"), self.accepted)
		QtCore.QObject.connect(file_dialog, QtCore.SIGNAL("rejected"), self.rejected)
		QtCore.QObject.connect(file_dialog, QtCore.SIGNAL("currentChanged(QString)"), self.changed)
		#file_dialog.show()
		#file_dialog.hide()
		
		self.file_dialogs.append(file_dialog)
			
		return file_dialog
	
class EMDesktopFrame(EMFrame):
	def __init__(self,parent,geometry=Region(0,0,0,0)):
		EMFrame.__init__(self,parent,geometry)
		self.display_frames = []
		
		EMDesktop.main_widget.register_resize_aware(self)
		
	def set_geometry(self,geometry):
		EMFrame.set_geometry(self,geometry)
		#try:
			#for child in self.children:
				#if isinstance(child,EMDesktopTaskWidget):
					#child.set_cam_pos(-self.parent.width()/2.0+child.width()/2.0,self.parent.height()/2.0-child.height()/2.0,0)
					
		#except: pass

	def updateGL(self):
		self.parent.updateGL()
	
	def attach_display_child(self,child):
		self.display_frames.append(child)
		EMWindowNode.attach_child(self,child)
	
	def get_display_child(self,idx=0):
		try: return self.display_frames[idx]
		except:
			print "warning, attempted to ask for a display child that did not exist"
			print "asked for child number ",idx,"but I have only",len(self.display_frames),"display children"
	
	def resize_gl(self):
	
		self.set_geometry(Region(0,0,int(viewport_width()),int(viewport_height())))
		if len(self.display_frames) != 0:
			self.display_frames[0].set_geometry(Region(200,0,-100,int(viewport_width()-400),int(viewport_height()),100))
	
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

def print_node_hierarchy(node):
	try: children = node.get_children()
	except : children = []
	if len(children) == 0:
		print ""
	else:
		for child in children:
			print child,
			print_node_hierarchy(child)

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
		
		
		self.modules = [] # a list of all the modules that currently exist
		self.gq=0			# quadric object for cylinders, etc	
		self.app=QtGui.QApplication.instance()
		self.sysdesktop=self.app.desktop()
		self.appscreen=self.sysdesktop.screen(self.sysdesktop.primaryScreen())
		self.frame_dl = 0 # display list of the desktop frame
		self.fov = 35
		self.resize_aware_objects = []
		self.animatables = []
		
		# what is this?
		self.bgob2=ob2dimage(self,self.read_EMAN2_image())
		
		self.setMouseTracking(True)
		
		# this float widget has half of the screen (the left)
		self.desktop_frame = EMDesktopFrame(self)
		self.left_side_bar = LeftSideWidgetBar(self.desktop_frame)
		self.display_frame = EMGLViewContainer(self)
		fw1 = EMDesktopTaskWidget(self.left_side_bar)
		self.left_side_bar.attach_child(fw1)
		self.desktop_frame.attach_child(self.left_side_bar)
		self.desktop_frame.attach_display_child(self.display_frame)
		
		print_node_hierarchy(self.desktop_frame)
		#print fw1.width(),fw1.height()
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
	
	def attach_module(self,module):
		self.modules.append(module)
	
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
		glMaterial(GL_FRONT,GL_SHININESS,1.0)
		glDisable(GL_TEXTURE_2D)
		glEnable(GL_LIGHTING)
		glCallList(self.frame_dl)
		
	def read_EMAN2_image(self):
		#self.p = QtGui.QPixmap("EMAN2.0.big2.jpg")
		self.p = QtGui.QPixmap.grabWindow(self.appscreen.winId(),0.0,0.0,self.sysdesktop.width(),self.sysdesktop.height()-30)
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
		glLightfv(GL_LIGHT0, GL_AMBIENT, [0.1, 0.1, 0.1, 1.0])
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_POSITION, [0.1,.1,1.,1.])
		#glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 1.5)
		#glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0.5)
		#glLightf(GL_LIGHT0, GL_QUADRATIC_ATTENUATION, .2)
		#glEnable(GL_LIGHT1)
		#glLightfv(GL_LIGHT1, GL_AMBIENT, [0.1, 0.1, 0.1, 1.0])
		#glLightfv(GL_LIGHT1, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
		#glLightfv(GL_LIGHT1, GL_SPECULAR, [0.1, .1, .1, 1.0])
		#glLightfv(GL_LIGHT1, GL_POSITION, [-0.1,.1,1.,1.])
		#glLightf(GL_LIGHT1, GL_SPOT_DIRECTION, 45)
		#glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, [-0.1,.1,1.,0.])
			
		glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,GL_TRUE)

		glEnable(GL_NORMALIZE)
		#glEnable(GL_RESCALE_NORMAL)
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
		
		#glPushMatrix()
		#glTranslate(self.width()/2.0,0,self.get_z_opt())
		##glLightfv(GL_LIGHT0, GL_POSITION, [0,0,1,1.])
		#glPopMatrix()
		
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
			#glEnable(GL_NORMALIZE)
			#glDisable(GL_RESCALE_NORMAL)
			self.desktop_frame.draw()
			#glDisable(GL_NORMALIZE)
			#glEnable(GL_RESCALE_NORMAL)
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

	def get_depth_for_height(self, height):
		try: 
			return EMDesktop.main_widget.get_depth_for_height(height)
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
			self.desktop_task_widget = EMGLViewQtWidget(EMDesktop.main_widget)
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
		return EMDesktop.main_widget.bindTexture(pixmap)
	
	def deleteTexture(self,val):
		return EMDesktop.main_widget.deleteTexture(val)
	
	def get_render_dims_at_depth(self, depth):
		try: return EMDesktop.main_widget.get_render_dims_at_depth(depth)
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
		
class ob2dimage:
	def __init__(self,target,pixmap):
		self.pixmap=pixmap
		self.target=target
		self.target.makeCurrent()
		self.texture_dl = 0
		self.itex=self.target.bindTexture(self.pixmap)
		
	def __del__(self):
		if self.texture_dl != 0: 
			glDeleteLists(self.texture_dl,1)
			self.texture_dl = 0

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
		if self.texture_dl == 0:
			self.texture_dl=glGenLists(1)
			
			if self.texture_dl == 0:
				return
			
			glNewList(self.texture_dl,GL_COMPILE)
	
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
			glEndList()
		
		if self.texture_dl != 0: glCallList(self.texture_dl)
		#glPopMatrix()


			
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = EMDesktop()
#	window.showFullScreen()
	window.app.exec_()
