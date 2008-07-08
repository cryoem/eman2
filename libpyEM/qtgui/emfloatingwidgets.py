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
#

# EMFloatingWidgets.py  Steve Ludtke  08/06/2006
# An experimental version of emimage.py that displays an image in a 3D context
# using texture mapping. Not a fully fleshed out class at this point

from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from OpenGL.GL import *
from OpenGL.GLU import *
from valslider import ValSlider

from EMAN2 import *
from emimageutil import *
from emimage2d import *
from emimage3d import *
from emimagemx import *
from math import sqrt

from emimage import EMImage

from emglobjects import EMViewportDepthTools, Camera2, EMBasicOpenGLObjects, Camera
from emimage2dtex import *


height_plane = 500

class EM3DWidgetVolume:
	'''
	a EM3DWidgetVolume has width(), height(), and depth() functions, and the associated private variables
	Inheriting functions should define the __determine_dimensions function which will be called implicitly if the update_dims flag is true
	In it they should define the member variables left,right,bottom,top,near and far.
	'''
	def __init__(self):
		self.left = 0
		self.right = 0
		self.bottom = 0
		self.top = 0
		self.near = 0
		self.far = 0

		self.update_dims = True # a flag, when set, cause this widget to re-evaluate its maximum dimensions
	
	def update(self):
		self.update_dims = True
		
	def width(self):
		if self.update_dims:
			self.determine_dimensions()
			
		return int(self.right -self.left)
		
	def height(self):
		if self.update_dims:
			self.determine_dimensions()
			
		return int(self.top - self.bottom)
		
	def depth(self):
		if self.update_dims:
			self.determine_dimensions()
			
		return int(self.near - self.far)
	
	def get_lr_bt_nf(self):
		'''
		get left-right bottom-top near-far
		'''
		if self.update_dims:
			self.determine_dimensions()
			
		return [self.left,self.right,self.bottom,self.top,self.near,self.far]
		
class EM3DWidget:
	'''
	A class for managing a 3D object as an interactive widget
	'''
	def __init__(self,parent,target):
		self.parent = parent
		self.target = target
		
		self.cam = Camera2(self) # a camera/orientation/postion object
		self.vdtools = EMViewportDepthTools(self) # viewport depth tools - essential for rerouting mouse events correctly
		
		self.corner_sets = [] # corner index sets of the 3D volume (?)
		self.planes = []	# string names for visible planes, used in paintGL
		self.model_matrices = [] # stores up to 3 OpenGL model view matrices (4x4 lists or tuples)
		
		self.draw_frame = True
		
	
	def viewportHeight(self):
		return self.parent.height()
	
	def viewportWidth(self):
		return self.parent.width()
	
	def width(self):
		return self.target.width()
		
	def height(self):
		return self.target.height()
	
	def depth(self):
		return self.target.depth()
	
	def get_lr_bt_nf(self):
		return self.target.get_lr_bt_nf()
	
	def __atomic_draw_frame(self,plane_string):
		
		if self.vdtools.drawFrame(True):
			self.corner_sets.append(self.vdtools.getCorners())
			self.planes.append((plane_string))
			self.model_matrices.append(self.vdtools.getModelMatrix())
	
	def paintGL(self):
		#clear everything
		self.corner_sets = []
		self.planes = []
		self.model_matrices = []
		
		glTranslatef(0,0,-self.depth()/2.0) # is this necessary?
		self.cam.position()
		
		lighting = glIsEnabled(GL_LIGHTING)
		glEnable(GL_LIGHTING) # lighting is on to make the borders look nice
		
		glPushMatrix()
		self.target.render()
		glPopMatrix()
		
		if self.draw_frame:
			p = self.get_lr_bt_nf()
			points = []
			points.append((p[0],p[2],p[5]))
			points.append((p[0],p[2],p[4]))
			points.append((p[0],p[3],p[4]))
			points.append((p[0],p[3],p[5]))
			points.append((p[1],p[2],p[5]))
			points.append((p[1],p[2],p[4]))
			points.append((p[1],p[3],p[4]))
			points.append((p[1],p[3],p[5]))
			unprojected = self.vdtools.unproject_points(points)
			
			# left zy plane
			glPushMatrix()
			self.vdtools.set_mouse_coords(unprojected[0],unprojected[1],unprojected[2],unprojected[3])
			self.__atomic_draw_frame('zy')
			glPopMatrix()
			
			glPushMatrix()
			self.vdtools.set_mouse_coords(unprojected[5],unprojected[4],unprojected[7],unprojected[6])
			self.__atomic_draw_frame('yz')
			glPopMatrix()
			
			glPushMatrix()
			self.vdtools.set_mouse_coords(unprojected[0],unprojected[4],unprojected[5],unprojected[1])
			self.__atomic_draw_frame('xz')
			glPopMatrix()
			
			glPushMatrix()
			self.vdtools.set_mouse_coords(unprojected[3],unprojected[2],unprojected[6],unprojected[7])
			self.__atomic_draw_frame('zx')
			glPopMatrix()
			
			glPushMatrix()
			self.vdtools.set_mouse_coords(unprojected[0],unprojected[3],unprojected[7],unprojected[4])
			self.__atomic_draw_frame('yx')
			glPopMatrix()
			
			glPushMatrix()
			self.vdtools.set_mouse_coords(unprojected[1],unprojected[5],unprojected[6],unprojected[2])
			self.__atomic_draw_frame('xy')
			glPopMatrix()
	
	def mousePressEvent(self, event):
		#if event.button()==Qt.MidButton or (event.button()==Qt.LeftButton and event.modifiers()&Qt.ControlModifier and self.inspector == None):	
			#self.drawable.initInspector()
			#self.drawable.inspector.show()
			#self.drawable.inspector.hide()
			#self.parent.addQtWidgetDrawer(self.getInspector())
			
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.mousePressEvent(event)
		else:
			pass
			#l=self.vdtools.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
			#qme=QtGui.QMouseEvent(event.type(),QtCore.QPoint(l[0],l[1]),event.button(),event.buttons(),event.modifiers())
			#self.drawable.mousePressEvent(qme)

	
	def wheelEvent(self,event):
		try: e = event.modifiers()
		except: return
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.wheelEvent(event)
		else:
			self.target.wheelEvent(event)
			
		#self.updateGL()
	
	def mouseMoveEvent(self,event):
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.mouseMoveEvent(event)
		else:
			pass
			#l=self.vdtools.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
			#qme=QtGui.QMouseEvent(event.type(),QtCore.QPoint(l[0],l[1]),event.button(),event.buttons(),event.modifiers())
			#self.drawable.mouseMoveEvent(qme)
			#self.drawable.mouseMoveEvent(event)
		
		#self.updateGL()

	def mouseReleaseEvent(self,event):
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.mouseReleaseEvent(event)
		else:
			pass
			#l=self.vdtools.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
			#qme=QtGui.QMouseEvent(event.type(),QtCore.QPoint(l[0],l[1]),event.button(),event.buttons(),event.modifiers())
			#self.drawable.mouseReleaseEvent(qme)
			#self.drawable.mouseReleaseEvent(event)
	
	def isinwin(self,x,y):
		val = False
		for i,p in enumerate(self.corner_sets):
			if self.vdtools.isinwinpoints(x,y,p):
				val = True
				#self.target.cam.plane = self.planes[i]
				self.vdtools.setModelMatrix(self.model_matrices[i])
				break
		return val
	
	def eyeCoordsDif(self,x1,y1,x2,y2,mdepth=True):
		return self.vdtools.eyeCoordsDif(x1,y1,x2,y2,mdepth)
	
class EMGLRotaryWidget(EM3DWidgetVolume):
	'''
	A display rotary widget - consists of an ellipse  with widgets 'attached' to it.
	Visually, the ellipse would lay in the plane perpendicular to the screen, and the widgets would be display in the plane
	of the screen and be attached to the ellipse and rotate around it interactively
	
	'''
	def __init__(self,parent):
		EM3DWidgetVolume.__init__(self)
		self.parent = parent
		self.widgets = []	# the list of widgets to display
		self.displayed_widget = -1 # the index of the currently displayed widget in self.widgets
		self.target_displayed_widget = -1 # the index of a target displayed widget, used for animation purposes
		self.rotations = -1
		
		# elliptical parameters
		self.ellipse_a = 40 # x direction - widgets displayed in position 0 will be a distance of 20 away from the center
		self.ellipse_b = 400 # y direction - widgets at 90 degrees will be a distance of 100 into the screen, away from the elliptical center
		self.rot = 20*pi/180.0 # this is how much the ellipse should be rotated about the y axis - makes it so that the ellipse isn't 'head on' to the viewer, for visual effects. The orientation of the widgets is not affected (they are always oriented as though the ellipse is not rotated in plane.
		self.cos_rot = cos(self.rot)
		self.sin_rot = sin(self.rot)
		
		self.time = 0		# time indicates the current time used for the basis of animation.
		self.time_interval = .5 # 0.5 seconds for the animation to complete
		self.time_begin = 0 # records the time at which the animation was begun
		self.is_animated = False
		self.angle_information = None # will eventually be list used for both static display and animation
		
		self.dtheta_animations = None  # an array for dtheta animations
		
		self.animation_queue = []

	def add_widget(self,widget,set_current=False):
		'''
		adds a widget to those that already in the rotary
		if set_current is true will also make it so that the newly added widget has the focus of the rotary
		
		return 1 if a redraw is necessary
		return 0 if a redraw is not necessary 
		'''
		if not isinstance(widget,EMGLViewQtWidget):
			print "error, can only add instances of EMGLViewQtWidget to the EMGLRotaryWidget"
			return 0
		else:
			self.widgets.append(widget)
			self.update_dims = True
			if set_current == True and self.is_animated == False:
				self.target_displayed_widget =  len(self.widgets)-1
				self.__start_animation()

			return 1
		
	def animate(self,time):
		if self.time_begin == 0:
			self.time_begin = time
			
		self.time = time -self.time_begin

		if self.time > self.time_interval:
			self.time_begin = 0
			self.time = 0
			self.is_animated = False
			self.dtheta_animations = None
			self.rotations = 0
			return 0
		else:
			for i in range(len(self.widgets)):
				dt = self.__get_dt()
				dtheta = self.angle_information[i][2]
				self.angle_information[i][0] = self.angle_information[i][1]+dt*dtheta
		return 1
	
	def render(self):
		dtheta = self.__get_current_dtheta()
		
		if self.angle_information == None:
			self.angle_information = []
			for i in range(len(self.widgets)):
				angle = i*dtheta
				self.angle_information.append([angle,angle,0])
		
		if self.is_animated: dt = self.__get_dt()
		
		for i,widget in enumerate(self.widgets):
			n = self.angle_information[i][0]
			n_rad = n*pi/180.0

			dx = self.ellipse_a*cos(n_rad)*self.cos_rot-self.ellipse_b*sin(n_rad)*self.sin_rot
			dz = -(self.ellipse_b*sin(n_rad)*self.cos_rot+self.ellipse_a*cos(n_rad)*self.sin_rot)
			h_width = widget.width()/2
			#print "positioning at",dx,dz," angle is", n
			glPushMatrix()
			glTranslate(dx,0,dz)
			glRotate(n,0,1,0)
			glTranslate(h_width,0,0)
			widget.paintGL()
			glPopMatrix()
			
	def __start_animation(self,counter_clockwise=True):
		self.is_animated = True
		self.__gen_animation_dthetas(counter_clockwise)
		self.parent.register_animatable(self)
		
	def wheelEvent(self,event):
		if not self.is_animated: 
			if event.delta() > 0:
				self.target_displayed_widget =  len(self.widgets)-1
			else:
				self.target_displayed_widget = 0
		else:
			if event.delta() > 0:
				if self.target_displayed_widget == 0:
					self.target_displayed_widget =  len(self.widgets)-1
				else: 
					self.target_displayed_widget -= 1
			else:
				if self.target_displayed_widget == len(self.widgets)-1:
					self.target_displayed_widget = 0
				else:
					self.target_displayed_widget += 1
		#print self.target_displayed_widget
		if not self.is_animated: 
			self.__start_animation(event.delta()> 0)
		else: self.__update_animation_dthetas(event.delta() > 0)
	
	def determine_dimensions(self):
		# first determine the ellipse offset - this is related to the inplane rotation of the ellipse. We want the back point corresponding to
		# a rotation of 90 degrees around the boundary of the ellipse (starting at [ellipse_a/2,0,0]) to be visible
		
		xoffset = -self.sin_rot*self.ellipse_b
		yoffset = self.cos_rot*self.ellipse_b
		
		self.left = xoffset
		
		height = -1
		width = -1
		for i in self.widgets:
			if i.width() > width: width = i.width()
			if i.height() > height: height = i.height()
		
		self.right = self.ellipse_a*self.cos_rot+width
		self.bottom = -height/2
		self.top = height/2
		
		self.far = -yoffset -width
		self.near = yoffset + width
		
		self.update_dims = False
		
	def set_update_P_inv(self,val=True):
		self.vdtools.set_update_P_inv(val)
	
	# PRIVATE 
	
	def __get_current_dtheta(self):
		if len(self.widgets) == 1: return 0
		else: return 90/(len(self.widgets)-1)

	def __get_dt(self):
		return sin(self.time/self.time_interval*pi/2)

	def __update_animation_dthetas(self,counter_clockwise=True):
		if not counter_clockwise: #clockwise OpenGL rotations are negative
			self.rotations -= 1
			rotations = -1
		else: # counter clockwise OpenGL rotations are positive
			self.rotations += 1
			rotations = 1
		
		
		n = len(self.widgets)
		dtheta = self.__get_current_dtheta()
		dt = self.__get_dt()
		current_thetas = []
		for i in range(n): 
			current_thetas.append(i*dtheta)
		
		if rotations == -1: # clockwise
			for i in range(n):
				for rot in range(self.rotations,self.rotations-1,-1):
					idx1 = (i+rot+1)%n
					idx2 = (i+rot)%n
					self.angle_information[i][1] =  self.angle_information[i][0]
					if idx2 != (n-1):
						self.angle_information[i][2] +=  current_thetas[idx2] - current_thetas[idx1]
					else:
						self.angle_information[i][2] +=  -( 360-(current_thetas[idx2] - current_thetas[idx1]))
		elif rotations >= 1: # counterclockwise
			for i in range(len(self.widgets)):
				for rot in range(self.rotations,self.rotations+1):
					idx1 = (i+rot-1)%n
					idx2 = (i+rot)%n
					self.angle_information[i][1] =  self.angle_information[i][0]
					if idx1 != (n-1):
						self.angle_information[i][2] +=  current_thetas[idx2] - current_thetas[idx1]
					else:
						self.angle_information[i][2] +=  360-(current_thetas[idx1] - current_thetas[idx2])

	def __gen_animation_dthetas(self,counter_clockwise=True):
		# find the shortest path from the currently displayed to the target displayed widget
				
		#c_distance = self.target_displayed_widget # clockwise distance, looking down -y
		#cc_distance = len(self.widgets) - self.target_displayed_widget #conterclockwise distance, looking down -y
		
		n = len(self.widgets)
		dtheta = self.__get_current_dtheta()
		self.angle_information = []
		current_thetas = []
		for i in range(n): 
			# information stored in current, start, end format
			current_thetas.append(i*dtheta)
			self.angle_information.append([i*dtheta,i*dtheta,0])
	
		if not counter_clockwise: #clockwise OpenGL rotations are negative
			self.rotations = -1
		else: # counter clockwise OpenGL rotations are positive
			self.rotations = 1
		
		
		if self.rotations <= -1: # clockwise
			for i in range(n):
				for rot in range(-1,self.rotations-1,-1):
					idx1 = (i+rot+1)%n
					idx2 = (i+rot)%n
					if idx2 != (n-1):
						self.angle_information[i][2] +=  current_thetas[idx2] - current_thetas[idx1]
					else:
						self.angle_information[i][2] +=  -( 360-(current_thetas[idx2] - current_thetas[idx1]))
		elif self.rotations >= 1: # counterclockwise
			for i in range(len(self.widgets)):
				for rot in range(1,self.rotations+1):
					idx1 = (i+rot-1)%n
					idx2 = (i+rot)%n
					if idx1 != (n-1):
						self.angle_information[i][2] +=  current_thetas[idx2] - current_thetas[idx1]
					else:
						self.angle_information[i][2] +=  360-(current_thetas[idx1] - current_thetas[idx2])
		else:
			print "error - can't rotate when rotations are set to 0"

class EMGLView3D:
	"""
	A view of an EMAN2 3D type, such as an isosurface or a 
	volume rendition, etc.
	"""
	def __init__(self, parent,image=None):
		self.parent = parent
		self.cam = Camera2(self)
		self.cam.motiondull = 3.0
		self.cam.setCamTrans('default_z',-parent.get_depth_for_height(height_plane))
		
		self.w = image.get_xsize()	# width of window
		self.h = image.get_ysize()	# height of window
		self.d = image.get_zsize()	# depth of the window
		self.sizescale = 1.0		# scale/zoom factor
		self.changefactor = 1.1		# used to zoom
		self.invchangefactor = 1.0/self.changefactor # used to invert zoom
		
		self.drawable = EMImage3DCore(image,self)		# the object that is drawable (has a draw function)
		self.drawable.cam.basicmapping = True
		self.drawable.cam.motiondull = 3.0
		#self.drawable.suppressInspector = True
		self.vdtools = EMViewportDepthTools(self)
		
		self.updateFlag = True
		
		self.drawFrame = True
		
		self.inspector = None
		
		self.psets = []
		self.modelmatrices = []
	
	def setOptScale(self):
		dims = self.drawable.getDataDims()
		
		xscale = self.w/dims[0]
		yscale = self.h/dims[1]
		
		if yscale < xscale: xscale = yscale
		
		print xscale
		self.drawable.cam.scale = xscale
	
	def setWidth(self,w):
		self.w = w
		self.drawable.resizeEvent(self.width(),self.height())
		
	def setDepth(self,d):
		self.d = d
		self.drawable.resizeEvent(self.width(),self.height())
		
	def setHeight(self,h):
		self.h = h
		self.drawable.resizeEvent(self.width(),self.height())
		
	def width(self):
		try:
			return int(self.w)
		except:
			return 0
	
	def height(self):
		try:
			return int(self.h)
		except:
			return 0
	
	def depth(self):
		try:
			return int(self.d)
		except:
			return 0
	
	def setData(self,data):
		try: self.drawable.setData(data)
		except: pass
		
	def paintGL(self):
		self.psets = []
		self.planes = []
		self.modelmatrices = []
		glTranslatef(0,0,-self.depth()/2.0)
		self.cam.position()
		lighting = glIsEnabled(GL_LIGHTING)
		glEnable(GL_LIGHTING)
		
		glPushMatrix()
		
		self.drawable.render()
		glPopMatrix()
		
		glPushMatrix()
		glTranslatef(-self.width()/2.0,0,0)
		glRotatef(-90,0,1,0)
		self.vdtools.update(self.depth()/2.0,self.height()/2.0)
		if self.drawFrame: 
			if self.vdtools.drawFrame(True): 
				self.psets.append(self.vdtools.getCorners())
				self.planes.append(('zy'))
				self.modelmatrices.append(self.vdtools.getModelMatrix())
		glPopMatrix()
		
		glPushMatrix()
		glTranslatef(self.width()/2.0,0,0)
		glRotatef(90,0,1,0)
		self.vdtools.update(self.depth()/2.0,self.height()/2.0)
		if self.drawFrame: 
			if self.vdtools.drawFrame(True): 
				self.planes.append(('yz'))
				self.psets.append(self.vdtools.getCorners())
				self.modelmatrices.append(self.vdtools.getModelMatrix())
				
		glPopMatrix()
		
		glPushMatrix()
		glTranslatef(0,self.height()/2.0,0)
		glRotatef(-90,1,0,0)
		self.vdtools.update(self.width()/2.0,self.depth()/2.0)
		if self.drawFrame: 
			if self.vdtools.drawFrame(True): 
				self.planes.append(('xz'))
				self.psets.append(self.vdtools.getCorners())
				self.modelmatrices.append(self.vdtools.getModelMatrix())
		glPopMatrix()
		
		glPushMatrix()
		glTranslatef(0,-self.height()/2.0,0)
		glRotatef(90,1,0,0)
		self.vdtools.update(self.width()/2.0,self.depth()/2.0)
		if self.drawFrame: 
			if self.vdtools.drawFrame(True): 
				self.planes.append(('zx'))
				self.psets.append(self.vdtools.getCorners())
				self.modelmatrices.append(self.vdtools.getModelMatrix())
		glPopMatrix()
		
		glPushMatrix()
		glTranslatef(0,0,-self.depth()/2.0)
		glRotatef(180,0,1,0)
		self.vdtools.update(self.depth()/2.0,self.height()/2.0)
		if self.drawFrame: 
			if self.vdtools.drawFrame(True):
				self.planes.append(('yx'))
				self.psets.append(self.vdtools.getCorners())
				self.modelmatrices.append(self.vdtools.getModelMatrix())
		glPopMatrix()
		
		glPushMatrix()
		glTranslatef(0,0,self.depth()/2.0)
		self.vdtools.update(self.width()/2.0,self.height()/2.0)
		if self.drawFrame: 
			if self.vdtools.drawFrame(True): 
				self.planes.append(('xy'))
				self.psets.append(self.vdtools.getCorners())
				self.modelmatrices.append(self.vdtools.getModelMatrix())
		glPopMatrix()
		
		if not lighting: glDisable(GL_LIGHTING)
		
	def viewportHeight(self):
		return self.parent.height()
	
	def viewportWidth(self):
		return self.parent.width()
	
	def getInspector(self):
		if (self.inspector == None):
			if self.drawable == None:
				return None
			self.drawable.initInspector()
			self.drawable.inspector.show()
			self.drawable.inspector.hide()
			self.inspector = self.drawable.inspector
			
		return self.inspector
	
	def mousePressEvent(self, event):
		if event.button()==Qt.MidButton or (event.button()==Qt.LeftButton and event.modifiers()&Qt.ControlModifier and self.inspector == None):	
			self.drawable.initInspector()
			self.drawable.inspector.show()
			self.drawable.inspector.hide()
			self.parent.addQtWidgetDrawer(self.getInspector())
			
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.mousePressEvent(event)
		else:
			l=self.vdtools.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
			qme=QtGui.QMouseEvent(event.type(),QtCore.QPoint(l[0],l[1]),event.button(),event.buttons(),event.modifiers())
			self.drawable.mousePressEvent(qme)

	
	def wheelEvent(self,event):
		try: e = event.modifiers()
		except: return
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.wheelEvent(event)
		else:
			self.drawable.wheelEvent(event)
			
		#self.updateGL()
	
	def mouseMoveEvent(self,event):
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.mouseMoveEvent(event)
		else:
			l=self.vdtools.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
			qme=QtGui.QMouseEvent(event.type(),QtCore.QPoint(l[0],l[1]),event.button(),event.buttons(),event.modifiers())
			self.drawable.mouseMoveEvent(qme)
			#self.drawable.mouseMoveEvent(event)
		
		#self.updateGL()

	def mouseReleaseEvent(self,event):
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.mouseReleaseEvent(event)
		else:
			l=self.vdtools.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
			qme=QtGui.QMouseEvent(event.type(),QtCore.QPoint(l[0],l[1]),event.button(),event.buttons(),event.modifiers())
			self.drawable.mouseReleaseEvent(qme)
			#self.drawable.mouseReleaseEvent(event)
	
	def update(self):
		pass
		self.parent.updateGL()
	
	def updateGL(self):
		self.parent.updateGL()
	
	def isinwin(self,x,y):
		val = False
		for i,p in enumerate(self.psets):
			if self.vdtools.isinwinpoints(x,y,p):
				val = True
				self.drawable.cam.plane = self.planes[i]
				self.vdtools.setModelMatrix(self.modelmatrices[i])
				break
		return val
	
	def eyeCoordsDif(self,x1,y1,x2,y2,mdepth=True):
		return self.vdtools.eyeCoordsDif(x1,y1,x2,y2,mdepth)
	
	def leaveEvent(self):
		pass
	
	def toolTipEvent(self,event):
		pass
	
	def set_update_P_inv(self,val=True):
		self.vdtools.set_update_P_inv(val)
		
	def get_render_dims_at_depth(self, depth):
		return self.parent.get_render_dims_at_depth(depth)
	
	def getNearPlaneDims(self):
		return self.parent.getNearPlaneDims()
		
	def getStartZ(self):
		return self.parent.getStartZ()
	
class EMGLView2D:
	"""
	A view of a 2D drawable type, such as a single 2D image or a matrix of 2D images
	
	"""
	def __init__(self, parent,image=None):
		self.parent = parent
		self.cam = Camera2(self)
		self.cam.setCamTrans('default_z',-parent.get_depth_for_height(height_plane))
		
		
		if isinstance(image,list):
			if len(image) == 1:
				self.become2DImage(image[0])
			else:
				self.drawable = EMImageMXCore(image,self)
				self.w = image[0].get_xsize()
				self.h = image[0].get_ysize()
		elif isinstance(image,EMData):
			self.become2DImage(image)
		
		self.drawable.suppressInspector = True
		self.initflag = True
		self.vdtools = EMViewportDepthTools(self)
		
		self.updateFlag = True
		
		self.drawFrame = True
		
		self.sizescale = 1.0
		self.changefactor = 1.1
		self.invchangefactor = 1.0/self.changefactor
		
		self.inspector = None
		
	def become2DImage(self,a):
		self.drawable = EMImage2DCore(a,self)
		#self.drawable.originshift = False
		self.w = a.get_xsize()
		self.h = a.get_ysize()
		
	def eyeCoordsDif(self,x1,y1,x2,y2,mdepth=True):
		return self.vdtools.eyeCoordsDif(x1,y1,x2,y2,mdepth)
	
	def set_update_P_inv(self,val=True):
		self.vdtools.set_update_P_inv(val)
	
	def setWidth(self,w):
		self.w = w
		self.drawable.resizeEvent(self.width(),self.height())
		
	def setHeight(self,h):
		self.h = h
		self.drawable.resizeEvent(self.width(),self.height())

	
	def width(self):
		try:
			return int(self.sizescale*self.w)
		except:
			return 0
	
	def height(self):
		try:
			return int(self.sizescale*self.h)
		except:
			return 0
	
	def setData(self,data):
		self.drawable.setData(data)
		
	def initializeGL(self):
		self.drawable.initializeGL()
	
	def viewportHeight(self):
		return self.parent.height()	
	
	def viewportWidth(self):
		return self.parent.width()
	
	def testBoundaries(self):
		'''
		Called when the image is first drawn, this resets the dimensions of this object
		if it is larger than the current size of the viewport. It's somewhat of a hack,
		but it's early stages in the design
		'''
		h = self.vdtools.getMappedHeight()
		w = self.vdtools.getMappedWidth()
		
		if ( w > self.viewportWidth() ):
			self.w = self.viewportWidth()/self.sizescale
		if ( h > self.viewportHeight() ):
			self.h = self.viewportHeight()/self.sizescale
	
	def paintGL(self):
		self.cam.position()
		self.vdtools.update(self.width()/2.0,self.height()/2.0)
		if (self.initflag == True):
			self.testBoundaries()
			self.initflag = False

		
		#self.mediator.checkBoundaryIssues()
		if (self.updateFlag):
			self.drawable.resizeEvent(self.width(),self.height())
			self.updateFlag = False
		glPushMatrix()
		glTranslatef(-self.width()/2.0,-self.height()/2.0,0)
		try: self.drawable.render()
		except Exception, inst:
			print type(inst)     # the exception instance
			print inst.args      # arguments stored in .args
			print int
		glPopMatrix()
		lighting = glIsEnabled(GL_LIGHTING)
		glEnable(GL_LIGHTING)
		if self.drawFrame: self.vdtools.drawFrame()
		if not lighting: glDisable(GL_LIGHTING)
		
	def update(self):
		self.parent.updateGL()
	
	def updateGL(self):
		self.parent.updateGL()
	
	def getInspector(self):
		if (self.inspector == None):
			if self.drawable == None:
				return None
			self.drawable.initInspector()
			self.drawable.inspector.show()
			self.drawable.inspector.hide()
			self.inspector = self.drawable.inspector
			
		return self.inspector
	
	def mousePressEvent(self, event):
		if event.button()==Qt.MidButton or (event.button()==Qt.LeftButton and event.modifiers()&Qt.ControlModifier and self.inspector == None):	
			self.drawable.initInspector()
			self.drawable.inspector.show()
			self.drawable.inspector.hide()
			self.parent.addQtWidgetDrawer(self.getInspector())
			
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.mousePressEvent(event)
		else:
			l=self.vdtools.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
			qme=QtGui.QMouseEvent(event.type(),QtCore.QPoint(l[0],l[1]),event.button(),event.buttons(),event.modifiers())
			self.drawable.mousePressEvent(qme)
			#self.drawable.mousePressEvent(event)
		
		#self.updateGL()
	
	def scaleEvent(self,delta):
		if ( delta > 0 ):
			self.sizescale *= self.changefactor
		elif ( delta < 0 ):
			self.sizescale *= self.invchangefactor

		self.drawable.resizeEvent(self.width(),self.height())
	
	def wheelEvent(self,event):
		if event.modifiers() == Qt.ShiftModifier:
			self.scaleEvent(event.delta())
			#print "updating",self.drawWidth(),self.drawHeight()
		else:
			self.drawable.wheelEvent(event)
			
		#self.updateGL()
	
	def mouseMoveEvent(self,event):
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.mouseMoveEvent(event)
		else:
			l=self.vdtools.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
			qme=QtGui.QMouseEvent(event.type(),QtCore.QPoint(l[0],l[1]),event.button(),event.buttons(),event.modifiers())
			self.drawable.mouseMoveEvent(qme)
			#self.drawable.mouseMoveEvent(event)
		
		#self.updateGL()

	def mouseReleaseEvent(self,event):
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.mouseReleaseEvent(event)
		else:
			l=self.vdtools.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
			qme=QtGui.QMouseEvent(event.type(),QtCore.QPoint(l[0],l[1]),event.button(),event.buttons(),event.modifiers())
			self.drawable.mouseReleaseEvent(qme)
			#self.drawable.mouseReleaseEvent(event)

		#self.updateGL()
	def emit(self, signal, event):
		try:
			QtCore.QObject.emit(signal,event)
		except:
			print "unknown signal", signal, "or unknown event",event
	
	def leaveEvent(self):
		self.drawable.leaveEvent()
	
	def toolTipEvent(self,event):
		pass
	
	def isinwin(self,x,y):
		return self.vdtools.isinwin(x,y)

class EMGLViewQtWidget:
	def __init__(self, parent=None, qwidget=None, widget_parent=None):
		self.parent = parent
		self.qwidget = qwidget
		self.drawFrame = True
		self.mapcoords = True
		self.itex = 0
		self.genTexture = True
		self.click_debug = False
		self.cam = Camera2(self)
		self.cam.setCamTrans('default_z',-parent.get_depth_for_height(height_plane))
		self.cam.motionRotate(0,0)
		self.borderwidth = 3.0
		self.glbasicobjects = EMBasicOpenGLObjects()
		self.setQtWidget(qwidget)
		self.P_inv = None
		self.childreceiver = None
		self.widget_parent = widget_parent
		
		self.current = None
		self.previous = None
		
		self.e2children = []
		self.is_child = False
		
		self.vdtools = EMViewportDepthTools(self)
	
	def __del__(self):
		if (self.itex != 0 ):
			self.parent.deleteTexture(self.itex)
		
	def set_update_P_inv(self,val=True):
		self.vdtools.set_update_P_inv(val)
	
	def width(self):
		return self.qwidget.width()
	
	def height(self):
		return self.qwidget.height()
	
	def viewportHeight(self):
		return self.parent.height()
	
	def viewportWidth(self):
		return self.parent.width()
	
	def setQtWidget(self, widget, delete_current = False):
		if ( delete_current and self.qwidget != None ):
			self.qwidget.deleteLater()
		
		self.qwidget = widget
		
		if ( widget != None ):
			#self.qwidget.setVisible(True)
			self.qwidget.setEnabled(True)
			self.genTexture = True
			self.updateTexture()
			
	def updateTexture(self):
		if ( self.itex == 0 or self.genTexture == True ) : 
			if (self.itex != 0 ):
				#passpyth
				self.parent.deleteTexture(self.itex)
			self.genTexture = False
			##print "binding texture"
			#self.qwidget.setVisible(True)
			#self.qwidget.repaint()
			if ( isinstance(self.qwidget,QtGui.QWidget) ):
				pixmap = QtGui.QPixmap.grabWidget(self.qwidget)
			else:
				pixmap = QtGui.QPixmap.grabWidget(self.qwidget.widget())
			#self.qwidget.setVisible(False)
			if (pixmap.isNull() == True ): print 'error, the pixmap was null'
			self.itex = self.parent.bindTexture(pixmap)
			if ( self.itex == 0 ): print 'Error - I could not generate the texture'
		
	def paintGL(self):
		#print "paintGL children"
		if (self.qwidget == None or self.itex == 0) :
			#print "no widget - paintGL children return" 
			return
		
		self.cam.debug = True
		#self.cam.position()
		
		# make sure the vdtools store the current matrices
		self.vdtools.update(self.width()/2.0,self.height()/2.0)
		
		glPushMatrix()
		glEnable(GL_TEXTURE_2D)
		glBindTexture(GL_TEXTURE_2D,self.itex)
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
		glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE)
		glBegin(GL_QUADS)
		glTexCoord2f(0.,0.)
		glVertex(-self.qwidget.width()/2.0,-self.qwidget.height()/2.0)
		glTexCoord2f(1.,0.)
		glVertex( self.qwidget.width()/2.0,-self.qwidget.height()/2.0)
		glTexCoord2f(1.,1.)
		glVertex( self.qwidget.width()/2.0, self.qwidget.height()/2.0)
		glTexCoord2f(0.,1.)
		glVertex( -self.qwidget.width()/2.0,self.qwidget.height()/2.0)
		glEnd()
		glDisable(GL_TEXTURE_2D)
		glPopMatrix()
	
		lighting = glIsEnabled(GL_LIGHTING)
		glEnable(GL_LIGHTING)
		if self.drawFrame:
			try: self.vdtools.drawFrame()
			except Exception, inst:
				print type(inst)     # the exception instance
				print inst.args      # arguments stored in .args
				print int
		if (not lighting): glDisable(GL_LIGHTING)
		
		# now draw children if necessary - such as a qcombobox list view that has poppud up
		for i in self.e2children:
			glPushMatrix()
			try:
				i.paintGL()
			except Exception, inst:
				print type(inst)     # the exception instance
				print inst.args      # arguments stored in .args
				print int
			glPopMatrix()
	
	def isinwin(self,x,y):
		for i in self.e2children:
			if i.isinwin(x,y):
				self.childreceiver = i
				return True
		
		return self.vdtools.isinwin(x,y)
	
	def eyeCoordsDif(self,x1,y1,x2,y2):
		return self.vdtools.eyeCoordsDif(x1,y1,x2,y2)
			
	def mouseinwin(self,x,y,width,height):
		return self.vdtools.mouseinwin(x,y,width,height)

	def toolTipEvent(self,event):
		if ( self.childreceiver != None ):
			# this means this class already knows that the mouse event is in the child
			# that is being displayed
			self.childreceiver.toolTip(event)
			self.childreceiver = None
			return
		l=self.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
		cw=self.qwidget.childAt(l[0],l[1])
		if cw == None: 
			QtGui.QToolTip.hideText()
			self.genTexture = True
			self.updateTexture()
			return
	
		p1 = QtCore.QPoint(event.x(),event.y())
		p2 = self.parent.mapToGlobal(p1)
		QtGui.QToolTip.showText(p2,cw.toolTip())
	
	def wheelEvent(self,event):
		doElse = True
		try:
			if event.modifiers() == Qt.ShiftModifier:
				self.cam.wheelEvent(event)
				doElse = False
		except: pass
		if doElse:
			if ( self.childreceiver != None ):
				# this means this class already knows that the mouse event is in the child
				# that is being displayed
				self.childreceiver.wheelEvent(event)
				self.childreceiver = None
				return
			else:
				# if we have any children (i.e. a drop down combo box) it should now disappear
				if len(self.e2children) > 0:
					self.e2children.pop()
					return
			l=self.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
			cw=self.qwidget.childAt(l[0],l[1])
			if cw == None: return
			gp=self.qwidget.mapToGlobal(QtCore.QPoint(l[0],l[1]))
			lp=cw.mapFromGlobal(gp)
			qme=QtGui.QWheelEvent(lp,event.delta(),event.buttons(),event.modifiers(),event.orientation())
			QtCore.QCoreApplication.sendEvent(cw,qme)
			self.genTexture = True
			self.updateTexture()
	
	def mouseDoubleClickEvent(self, event):
		if ( self.childreceiver != None ):
			# this means this class already knows that the mouse event is in the child
			# that is being displayed
			self.childreceiver.mouseDoubleClickEvent(event)
			self.childreceiver = None
			return
		l=self.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
		cw=self.qwidget.childAt(l[0],l[1])
		if cw == None: return
		gp=self.qwidget.mapToGlobal(QtCore.QPoint(l[0],l[1]))
		lp=cw.mapFromGlobal(gp)
		if (isinstance(cw,QtGui.QComboBox)):
			print "it's a combo"
		else:
			qme=QtGui.mouseDoubleClickEvent(event.type(),lp,event.button(),event.buttons(),event.modifiers())
			#self.qwidget.setVisible(True)
			QtCore.QCoreApplication.sendEvent(cw,qme)
			#self.qwidget.setVisible(False)
		self.genTexture = True
		self.updateTexture()
		
	def get_depth_for_height(self,height_plane):
		return self.parent.get_depth_for_height(height_plane)
	
	def mousePressEvent(self, event):
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.mousePressEvent(event)
		else:
			if ( self.childreceiver != None ):
				# this means this class already knows that the mouse event is in the child
				# that is being displayed
				self.childreceiver.mousePressEvent(event)
				self.childreceiver = None
				return
			else:
				# if we have any children (i.e. a drop down combo box) it should now disappear
				if len(self.e2children) > 0:
					self.e2children.pop()
					return
				
			l=self.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
			cw=self.qwidget.childAt(l[0],l[1])
			if cw == None: return
			##print cw.objectName()
			gp=self.qwidget.mapToGlobal(QtCore.QPoint(l[0],l[1]))
			lp=cw.mapFromGlobal(gp)
			if (isinstance(cw,QtGui.QComboBox)):
				cw.showPopup()
				cw.hidePopup()
				widget = EMGLViewQtWidget(self.parent,None,cw);
				widget.setQtWidget(cw.view())
				widget.cam.loadIdentity()	
				widget.cam.setCamTrans("x",cw.geometry().x()-self.width()/2.0+cw.view().width()/2.0)
				widget.cam.setCamTrans("y",((self.height()/2.0-cw.geometry().y())-cw.view().height()/2.0))
				widget.cam.setCamTrans("z",0.1)
				widget.drawFrame = False
				self.e2children.append(widget)
				self.e2children[0].is_child = True
			else:
				qme=QtGui.QMouseEvent( event.type(),lp,event.button(),event.buttons(),event.modifiers())
				if (self.is_child): QtCore.QCoreApplication.sendEvent(self.qwidget,qme)
				else: QtCore.QCoreApplication.sendEvent(cw,qme)
				
			self.genTexture = True
			self.updateTexture()
		
	def mouseMoveEvent(self,event):
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.mouseMoveEvent(event)
		else:
			if ( self.childreceiver != None ):
				# this means this class already knows that the mouse event is in the child
				# that is being displayed
				self.childreceiver.mouseMoveEvent(event)
				self.childreceiver = None
				return
			else:
				l=self.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
				cw=self.qwidget.childAt(l[0],l[1])
				self.current = cw
				if ( self.current != self.previous ):
					QtGui.QToolTip.hideText()
					if ( self.current != None ):
						qme=QtCore.QEvent(QtCore.QEvent.Enter)
						QtCore.QCoreApplication.sendEvent(self.current,qme)
						
					if ( self.previous != None ):
						qme=QtCore.QEvent(QtCore.QEvent.Leave)
						QtCore.QCoreApplication.sendEvent(self.previous,qme)
				
				self.previous = self.current
				if cw == None:
					QtGui.QToolTip.hideText()
					if ( self.previous != None ):
						qme=QtCore.QEvent(QtCore.QEvent.Leave)
						QtCore.QCoreApplication.sendEvent(self.previous,qme)
						self.genTexture = True
						self.updateTexture()
					return
				gp=self.qwidget.mapToGlobal(QtCore.QPoint(l[0],l[1]))
				lp=cw.mapFromGlobal(gp)
				qme=QtGui.QMouseEvent(event.type(),lp,event.button(),event.buttons(),event.modifiers())
				QtCore.QCoreApplication.sendEvent(cw,qme)
			# FIXME
			# setting the genTexture flag true here causes the texture to be regenerated
			# when the mouse moves over it, which is inefficient.
			# The fix is to only set the genTexture flag when mouse movement
			# actually causes a change in the appearance of the widget (for instance, list boxes from comboboxes)
			self.genTexture = True
			self.updateTexture()

	def mouseReleaseEvent(self,event):
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.mouseReleaseEvent(event)
		else:
			if ( self.childreceiver != None ):
				# this means this class already knows that the mouse event is in the child
				# that is being displayed
				#try:
				self.childreceiver.mouseReleaseEvent(event)
				self.childreceiver = None
				self.e2children.pop()
				return
			else:
				# if we have any children (i.e. a drop down combo box) it should now disappear
				if len(self.e2children) > 0:
					self.e2children.pop()
					return
			
			l=self.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
			cw=self.qwidget.childAt(l[0],l[1])
			if cw == None: return
			gp=self.qwidget.mapToGlobal(QtCore.QPoint(l[0],l[1]))
			lp=cw.mapFromGlobal(gp)
			if (isinstance(cw,QtGui.QComboBox)):
				print "it's a combo"
				#cw.showPopup()
			else:
				qme=QtGui.QMouseEvent(event.type(),lp,event.button(),event.buttons(),event.modifiers())
				if (self.is_child):
					##print self.qwidget
					##print self.qwidget.currentIndex().row()
					##print self.widget_parent
					##print self.qwidget.rect().left(),self.qwidget.rect().right(),self.qwidget.rect().top(),self.qwidget.rect().bottom()
					##print lp.x(),lp.y()
					self.widget_parent.setCurrentIndex(self.qwidget.currentIndex().row())
					#self.widget_parent.changeEvent(QtCore.QEvent())
					#self.widget_parent.highlighted(self.qwidget.currentIndex().row())
					#self.qwidget.commitData(self.qwidget.parent())
					##print self.qwidget.currentText()
					#self.widget_parent.setVisible(True)
					#self.widget_parent.setEnabled(True)
					#self.qwidget.setVisible(True)
					#QtCore.QCoreApplication.sendEvent(self.widget_parent,qme)
					#self.qwidget.setVisible(False)
					self.widget_parent.emit(QtCore.SIGNAL("activated(QString)"),self.widget_parent.itemText(self.qwidget.currentIndex().row()))
				else:
					#self.qwidget.setVisible(True)
					QtCore.QCoreApplication.sendEvent(cw,qme)
					#self.qwidget.setVisible(False)
			
			self.genTexture = True
			self.updateTexture()
		
	def leaveEvent(self):
		if (self.current != None) : 
			qme = QtCore.QEvent(QtCore.QEvent.Leave)
			QtCore.QCoreApplication.sendEvent(self.current,qme)
			self.current = None
			self.previouse = None
			self.genTexture = True
			self.updateTexture()
			
	def enterEvent():
		pass
	def timerEvent(self,event=None):
		pass
		#self.cam.motionRotate(.2,.2)


class EMFloatingWidgets(QtOpenGL.QGLWidget):
	def __init__(self,parent=None):
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True)
		# enable multisampling to combat aliasing
		fmt.setSampleBuffers(True)
		# stenciling is for object dependent shading
		fmt.setStencil(True)
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		
		self.setMouseTracking(True)
		self.fov = 2*180*atan2(1,5)/pi
		
		self.floatwidget = EMFloatingWidgetsCore(self)
		
		self.cam = Camera()
		
		self.animatables = [] # an array of animatable objects - must have the animate(time) function
		
		self.timer = QTimer()
		QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.timeout)
		self.timer.start(10)
		
	def get_depth_for_height(self, height):
		# This function returns the width and height of the renderable 
		# area at the origin of the data volume
		depth = height/(2.0*tan(self.fov/2.0*pi/180.0))
		return depth
	
	def timeout(self):
		
		if len(self.animatables) == 0: return
		
		for i,animatable in enumerate(self.animatables):
			if not animatable.animate(time.time()):
				# this could be dangerous
				self.animatables.pop(i)
		
		self.updateGL()
	
	def register_animatable(self,animatable):
		self.animatables.append(animatable)
	
	def get_render_dims_at_depth(self, depth):
		# This function returns the width and height of the renderable 
		# area at the origin of the data volume
		height = -2*tan(self.fov/2.0*pi/180.0)*(depth)
		width = self.aspect*height
		return [width,height]

	
	def initializeGL(self):
		#print "initializeGL"
		glClearColor(0,0,0,0)
		
		glLightfv(GL_LIGHT0, GL_AMBIENT, [0.1, 0.1, 0.1, 1.0])
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_POSITION, [0.1,.1,1.,0.])
	
		glEnable(GL_LIGHTING)
		glEnable(GL_LIGHT0)
		glEnable(GL_DEPTH_TEST)
		
		glEnable(GL_NORMALIZE)
		
		glClearStencil(0)
		glEnable(GL_STENCIL_TEST)
		
		glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST)
		glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST)
		glHint(GL_TEXTURE_COMPRESSION_HINT, GL_NICEST)
	
		# enable multisampling to combat aliasing
		if ( "GL_ARB_multisample" in glGetString(GL_EXTENSIONS) ): glEnable(GL_MULTISAMPLE)
		else: glDisable(GL_MULTISAMPLE)
		
	def paintGL(self):
		#print "paintGL"
		glClear(GL_COLOR_BUFFER_BIT)
		if glIsEnabled(GL_DEPTH_TEST):
			glClear(GL_DEPTH_BUFFER_BIT)
		if glIsEnabled(GL_STENCIL_TEST):
			glClear(GL_STENCIL_BUFFER_BIT)
			
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
	
		self.floatwidget.render()
		
	def resizeGL(self, width, height):
		#print "resizeGL"
		glViewport(0,0,self.width(),self.height())
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		#glFrustum(-1.*width/height,1.*width/height, -1.,1., 5.,15.)
		
		# fov angle is the given by
		#self.fov = 2*180*atan2(1,5)/pi
		# aspect ratio is given by
		self.aspect = float(self.width())/float(self.height())
		# this is the same as the glFrustum call above
		depth = self.get_depth_for_height(height_plane)
		#gluPerspective(self.fov,self.aspect,depth-depth/4,depth+depth/4)
		gluPerspective(self.fov,self.aspect,1,10000)
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		try: self.floatwidget.resizeEvent(width,height)
		except: print "couldn't resize floatwidget"
		
	def mousePressEvent(self, event):
		self.floatwidget.mousePressEvent(event)
	
	def mouseMoveEvent(self, event):
		self.floatwidget.mouseMoveEvent(event)
		
	def mouseReleaseEvent(self, event):
		self.floatwidget.mouseReleaseEvent(event)

	def mouseDoubleClickEvent(self, event):
		self.floatwidget.mouseDoubleClickEvent(event)

	def wheelEvent(self, event):
		self.floatwidget.wheelEvent(event)

	def toolTipEvent(self, event):
		self.floatwidget.wheelEvent(event)
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
		self.floatwidget.hoverEvent(event)

class EMFloatingWidgetsCore:
	"""A QT widget for rendering EMData objects. It can display single 2D or 3D images 
	or sets of 2D images.
	"""
	def __init__(self, parent=None):
		#print "init"
		self.parent = parent
	
		self.imtex=0
		self.current = None
		self.previous = None
	
		self.initFlag = True
		self.qwidgets = []
		
		self.suppressUpdateGL = False
		self.rotary = None
		#print "init done"
	
	def register_animatable(self,animatable):
		self.parent.register_animatable(animatable)
		
	def get_depth_for_height(self, height):
		try: 
			return self.parent.get_depth_for_height(height)
		except:
			print "parent can't get height for depth"
			return 0

	def height(self):
		return self.parent.height()
	
	def width(self):
		return self.parent.width()

	def updateGL(self):
		if not self.suppressUpdateGL:
			try: self.parent.updateGL()
			except: pass

	def addQtWidgetDrawer(self,widget):
		w = EMGLViewQtWidget(self)
		w.setQtWidget(widget)
		self.qwidgets.append(w)
		
		#print "initializeGL done"
	def render(self):
		
		if ( self.initFlag == True ):
			self.initFlag = False
			self.fd = QtGui.QFileDialog(self.parent,"Open File",QtCore.QDir.currentPath(),QtCore.QString("Image files (*.img *.hed *.mrc)"))
			QtCore.QObject.connect(self.fd, QtCore.SIGNAL("finished(int)"), self.finished)
			self.fd.show()
			self.fd.hide()
			#self.qwidgets.append(EMGLViewQtWidget(self.parent))
			#self.qwidgets[0].setQtWidget(self.fd)
			##self.qwidgets[0].cam.setCamX(-100)
			
			
			
			self.fd2 = QtGui.QFileDialog(self.parent,"Open File",QtCore.QDir.currentPath(),QtCore.QString("Image files (*.img *.hed *.mrc)"))
			QtCore.QObject.connect(self.fd2, QtCore.SIGNAL("finished(int)"), self.finished)
			self.fd2.show()
			self.fd2.hide()
			#self.qwidgets.append(EMGLViewQtWidget(self.parent))
			#self.qwidgets[1].setQtWidget(self.fd2)
			#self.qwidgets[1].cam.setCamX(-200)
			
			a = EMGLViewQtWidget(self.parent)
			a.setQtWidget(self.fd2)
			b = EMGLViewQtWidget(self.parent)
			b.setQtWidget(self.fd)
			rotary = EMGLRotaryWidget(self)
			rotary.add_widget(b)
			rotary.add_widget(b)
			rotary.add_widget(b)
			rotary.add_widget(b)
			rotary.add_widget(b)
			rotary.add_widget(b)
			rotary.add_widget(a)
			rotary.add_widget(a)
			
			self.qwidgets.append(EM3DWidget(self,rotary))
			#self.rotary.add_widget(self.qwidgets[1])
			#self.rotary.add_widget(self.qwidgets[1])
		glPushMatrix()
		glTranslate(-100,0,-1250)
		for i in self.qwidgets:
			#print "getting opengl matrices"
			glPushMatrix()
			i.paintGL()
			glPopMatrix()
			#print "paint child done"
		
		
		glPopMatrix()
		#print "paintGL done"
	def finished(self,val):
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
		return self.parent.bindTexture(pixmap)
	
	def deleteTexture(self,val):
		return self.parent.deleteTexture(val)
	
	def get_render_dims_at_depth(self, depth):
		try: return self.parent.get_render_dims_at_depth(depth)
		except:
			print "parent can't get render dims at for depth"
			return

	def resizeEvent(self, width, height):
		for i in self.qwidgets:
			i.set_update_P_inv()
	
	def mousePressEvent(self, event):
		for i in self.qwidgets:
			if ( i.isinwin(event.x(),self.height()-event.y()) ):
				i.mousePressEvent(event)
				intercepted = True
				self.updateGL()
				return
	
	def mouseMoveEvent(self, event):
		for i in self.qwidgets:
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
		for i in self.qwidgets:
			if ( i.isinwin(event.x(),self.height()-event.y()) ):
				i.mouseReleaseEvent(event)
				self.updateGL()
				return
					
		
	def mouseDoubleClickEvent(self, event):
		for i in self.qwidgets:
			if ( i.isinwin(event.x(),self.height()-event.y()) ):
				i.mouseReleaseEvent(event)
				self.updateGL()
				return
		
		
	def wheelEvent(self, event):
		for i in self.qwidgets:
				if ( i.isinwin(event.x(),self.height()-event.y()) ):
					i.wheelEvent(event)
					self.updateGL()
					return

	def toolTipEvent(self, event):
		for i in self.qwidgets:
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
			for i in self.qwidgets:
				if ( i.isinwin(event.x(),self.height()-event.y()) ):
					i.hoverEvent(event)
					break
		self.updateGL()

# This is just for testing, of course
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = EMFloatingWidgets()
	window2 = EMParentWin(window)
	window2.show()
	
	sys.exit(app.exec_())
