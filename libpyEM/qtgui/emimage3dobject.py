#!/usr/bin/env python

# Author:  David Woolford 10/26/2007 (woolford@bcm.edu)
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

from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from OpenGL import GL,GLU,GLUT
from OpenGL.GL import *
from OpenGL.GLU import *
from valslider import ValSlider
from math import *
from EMAN2 import *

class Camera2:
	"""\brief A camera object encapsulates 6 degrees of freedom, and a scale factor
	
	The camera object stores x,y,z coordinates and a single transform object.
	For instance, it can be used for positioning and moving the 'camera' in OpenGL,
	however, instead of literally moving a camera, in OpenGL the scene itself is moved
	and the camera is generally thought of as staying still.
	
	Use the public interface of setCamTrans and motionRotate (which is based on mouse movement)_
	to move the camera position
	
	Then call 'position' in your main OpenGL draw function before drawing anything.
	
	"""
	def __init__(self,parent=None):
		# The magnification factor influences how the scale (zoom) is altered when a zoom event is received.
		# The equation is scale *= mag_factor or scale /= mag_factor, depending on the event.
		self.parent=parent
		self.mag_factor = 1.1
		
		self.t3d_stack = []
		self.loadIdentity()
		
		self.mmode = 0
		self.debug = False
		
	def loadIdentity(self):
		self.scale = 1.0
		
		# The camera coordinates
		self.cam_x = 0
		self.cam_y = 0
		self.cam_z = 0
		
		# Camera offsets - generally you want to set default_z to some negative
		# value the OpenGL scene is viewable.
		self.default_x = 0
		self.default_y = 0
		self.default_z = 0
		
		t3d = Transform3D()
		t3d.to_identity()
		self.t3d_stack.append(t3d)
		
	def position(self):
		# position the camera, regualar OpenGL movement.
		if (self.debug):
			print "Camera translational position",self.cam_x,self.cam_y,self.cam_z
		glTranslated(self.cam_x, self.cam_y, self.cam_z)
		
		
		rot = self.t3d_stack[len(self.t3d_stack)-1].get_rotation()
		if (self.debug):
			print "Camera rotation ",float(rot["phi"]),float(rot["alt"]),float(rot["az"])
		glRotate(float(rot["phi"]),0,0,1)
		glRotate(float(rot["alt"]),1,0,0)
		glRotate(float(rot["az"]),0,0,1)
		
		if (self.debug):
			print "Camera scale ",self.scale
		# here is where zoom is applied
		glScalef(self.scale,self.scale,self.scale)
		
	def scale_event(self,delta):
		if delta > 0:
			self.scale *= self.mag_factor
		elif delta < 0:
			self.scale *= 1.0/self.mag_factor
	
	def setCamTrans(self,axis,value):
		if ( axis == 'x'):
			self.setCamX(value)
		elif ( axis == 'y'):
			self.setCamY(value)
		elif ( axis == 'z'):
			self.setCamZ(value)
		elif ( axis == 'default_x'):
			self.default_x = value
		elif ( axis == 'default_y'):
			self.default_y = value
		elif ( axis == 'default_z'):
			self.default_z = value
			self.setCamZ(0)
		else:
			print 'Error, the axis (%s) specified is unknown. No action was taken' %axis
	
	def setCamZ(self,z):
		self.cam_z = self.default_z + z
		
	def setCamY(self,y):
		self.cam_y = self.default_y + y
		
	def setCamX(self,x):
		self.cam_x = self.default_x + x

	def motionRotate(self,x,y):
		# this function implements mouse interactive rotation
		# [x,y] is the vector generating by the mouse movement (in the plane of the screen)
		# Rotation occurs about the vector 90 degrees to [x,y,0]
		# The amount of rotation is linealy proportional to the length of [x,y]
		
		if ( x == 0 and y == 0): return
		
		theta = atan2(-y,x)

		rotaxis_x = sin(theta)
		rotaxis_y = cos(theta)
		rotaxis_z = 0
		
		length = sqrt(x*x + y*y)
		# 8.0 is a magic number - things rotate more if they are closer and slower if they are far away in this appproach
		# Or does it?
		# This magic number could be overcome using a strategy based on the results of get_render_dims_at_depth
		angle = length/8.0*pi
		
		t3d = Transform3D()
		quaternion = {}
		quaternion["Omega"] = angle
		quaternion["n1"] = rotaxis_x
		quaternion["n2"] = rotaxis_y
		quaternion["n3"] = rotaxis_z
		
		t3d.set_rotation( EULER_SPIN, quaternion )
		
		size = len(self.t3d_stack)
		self.t3d_stack[size-1] = t3d*self.t3d_stack[size-1]
		
	def setScale(self,val):
		self.scale = val
	
	def loadRotation(self,t3d):
		self.t3d_stack.append(t3d)

	def getThinCopy(self):
		# this is called a thin copy because it does not copy the entire t3d stack, just the last t3d
		cam = Camera()
		size = len(self.t3d_stack)
		cam.loadRotation(self.t3d_stack[size-1])
		
		cam.scale =	self.scale
		cam.cam_x = self.cam_x
		cam.cam_y = self.cam_y
		cam.cam_z = self.cam_z
		
		cam.default_x = self.default_x
		cam.default_y = self.default_y
		cam.default_z = self.default_z
		
		return cam
	
	def mousePressEvent(self, event):
		self.mpressx = event.x()
		self.mpressy = event.y()
		if event.button()==Qt.LeftButton:
			if self.mmode==0:
				# this is just a way of duplicating the last copy
				tmp =self.t3d_stack.pop()
				t3d = Transform3D(tmp)
				self.t3d_stack.append(tmp)
				self.t3d_stack.append(t3d)
		
	def mouseMoveEvent(self, event):
		if event.buttons()&Qt.LeftButton:
			if self.mmode==0:
				if event.modifiers() == Qt.ControlModifier:
					self.motionTranslate(event.x()-self.mpressx, self.mpressy - event.y())
				else:
					self.motionRotate(self.mpressx - event.x(), self.mpressy - event.y())
				self.mpressx = event.x()
				self.mpressy = event.y()
		elif event.buttons()&Qt.RightButton:
			if self.mmode==0:
				if event.modifiers() == Qt.ControlModifier:
					self.scale_event(event.y()-self.mpressy)	
				else:
					self.motionTranslateLA(self.mpressx, self.mpressy,event)
					
				self.mpressx = event.x()
				self.mpressy = event.y()
	
	def mouseReleaseEvent(self, event):
		if event.button()==Qt.LeftButton:
			if self.mmode==0:
				return
		elif event.button()==Qt.RightButton:
			if self.mmode==0:
				return
			
	def wheelEvent(self, event):
		self.scale_event(event.delta())
	
	def motionTranslateLA(self,prev_x,prev_y,event):
		#print "motion translate"
		[dx,dy] = self.parent.eyeCoordsDif(prev_x,self.parent.parentHeight()-prev_y,event.x(),self.parent.parentHeight()-event.y())
		#[wx2,wy2,wz2] = self.parent.eyeCoords(event.x(),self.parent.parentHeight()-event.y())
		#[wx2,wy2,wz2] =  self.parent.mouseViewportMovement(event.x(),self.parent.parentHeight()-event.y(),wx1,wy1,wz1,zprime)
		#self.parent.mouseViewportMovement(1,2,3,4)
		#[wx1,wy1] = self.parent.mouseinwin(prev_x,self.parent.parentHeight()-prev_y)
		#[wx2,wy2] = self.parent.mouseinwin(event.x(),self.parent.parentHeight()-event.y())
		self.cam_x += dx
		self.cam_y += dy
		#self.cam_z += wz2

class Camera:
	"""\brief A camera object encapsulates 6 degrees of freedom, and a scale factor
	
	The camera object stores x,y,z coordinates and a single transform object.
	For instance, it can be used for positioning and moving the 'camera' in OpenGL,
	however, instead of literally moving a camera, in OpenGL the scene itself is moved
	and the camera is generally thought of as staying still.
	
	Use the public interface of setCamTrans and motionRotate (which is based on mouse movement)_
	to move the camera position
	
	Then call 'position' in your main OpenGL draw function before drawing anything.
	
	"""
	def __init__(self):
		# The magnification factor influences how the scale (zoom) is altered when a zoom event is received.
		# The equation is scale *= mag_factor or scale /= mag_factor, depending on the event.
		self.mag_factor = 1.1
		self.scale = 1.0
		
		# The camera coordinates
		self.cam_x = 0
		self.cam_y = 0
		self.cam_z = 0
		
		# Camera offsets - generally you want to set default_z to some negative
		# value the OpenGL scene is viewable.
		self.default_x = 0
		self.default_y = 0
		self.default_z = 0
		
		# The Transform3D object stores the rotation
		t3d = Transform3D()
		t3d.to_identity()
	
		# At the moment there is a stack of Transform3D objects, this is for the purposes
		# of undoing actions. If the undo functionality was not required, the stack could be
		# removed.
		self.t3d_stack = []
		self.t3d_stack.append(t3d)
		
	def position(self):
		# position the camera, regualar OpenGL movement.
		glTranslated(self.cam_x, self.cam_y, self.cam_z)
		
		rot = self.t3d_stack[len(self.t3d_stack)-1].get_rotation()
		glRotate(float(rot["phi"]),0,0,1)
		glRotate(float(rot["alt"]),1,0,0)
		glRotate(float(rot["az"]),0,0,1)
		
		# here is where zoom is applied
		glScalef(self.scale,self.scale,self.scale)
		
	def scale_event(self,delta):
		if delta > 0:
			self.scale *= self.mag_factor
		elif delta < 0:
			self.scale *= 1.0/self.mag_factor
	
	def setCamTrans(self,axis,value):
		if ( axis == 'x'):
			setCamX(value)
		elif ( axis == 'y'):
			setCamY(value)
		elif ( axis == 'z'):
			setCamZ(value)
		elif ( axis == 'default_x'):
			self.default_x = value
		elif ( axis == 'default_y'):
			self.default_y = value
		elif ( axis == 'default_z'):
			self.default_z = value
		else:
			print 'Error, the axis (%s) specified is unknown. No action was taken' %axis
	
	def setCamZ(self,z):
		self.cam_z = self.default_z + z
		
	def setCamY(self,y):
		self.cam_y = self.default_y + y
		
	def setCamX(self,x):
		self.cam_x = self.default_x + x

	def motionRotate(self,x,y):
		# this function implements mouse interactive rotation
		# [x,y] is the vector generating by the mouse movement (in the plane of the screen)
		# Rotation occurs about the vector 90 degrees to [x,y,0]
		# The amount of rotation is linealy proportional to the length of [x,y]
		
		if ( x == 0 and y == 0): return
		
		theta = atan2(-y,x)

		rotaxis_x = sin(theta)
		rotaxis_y = cos(theta)
		rotaxis_z = 0
		
		length = sqrt(x*x + y*y)
		# 8.0 is a magic number - things rotate more if they are closer and slower if they are far away in this appproach
		# Or does it?
		# This magic number could be overcome using a strategy based on the results of get_render_dims_at_depth
		angle = length/8.0*pi
		
		t3d = Transform3D()
		quaternion = {}
		quaternion["Omega"] = angle
		quaternion["n1"] = rotaxis_x
		quaternion["n2"] = rotaxis_y
		quaternion["n3"] = rotaxis_z
		
		t3d.set_rotation( EULER_SPIN, quaternion )
		
		size = len(self.t3d_stack)
		self.t3d_stack[size-1] = t3d*self.t3d_stack[size-1]
		
	def setScale(self,val):
		self.scale = val
	
	def loadRotation(self,t3d):
		self.t3d_stack.append(t3d)

	def getThinCopy(self):
		# this is called a thin copy because it does not copy the entire t3d stack, just the last t3d
		cam = Camera()
		size = len(self.t3d_stack)
		cam.loadRotation(self.t3d_stack[size-1])
		
		cam.scale =	self.scale
		cam.cam_x = self.cam_x
		cam.cam_y = self.cam_y
		cam.cam_z = self.cam_z
		
		cam.default_x = self.default_x
		cam.default_y = self.default_y
		cam.default_z = self.default_z
		
		return cam

class EMImage3DObject:
	def __init__(self):
		self.rank = 0
		
	def render(self):
		pass

	def closeEvent(self,event) :
		pass
		
	def mousePressEvent(self, event):
		pass
		
	def mouseMoveEvent(self, event):
		pass
	
	def mouseReleaseEvent(self, event):
		pass
	
	def getType(self):
		pass
	
	# this function will be called when OpenGL recieves a resize event
	def resizeEvent(self):
		pass
	
	def setData(self):
		pass
	
	def showInspector(self):
		pass
	
	def getInspector(self):
		pass
	
	def setRank(self,rank):
		self.rank = rank
	
	def setName(self, name):
		self.name = name

	def getName(self):
		return self.name
	
	def getCurrentCamera(self):
		return self.cam.getThinCopy()
	
	def setCamera(self,camera):
		self.cam = camera
	
	def mousePressEvent(self, event):
#		lc=self.scrtoimg((event.x(),event.y()))
		if event.button()==Qt.MidButton:
			if not self.inspector or self.inspector ==None:
				return
			self.inspector.updateRotations(self.cam.t3d_stack[len(self.cam.t3d_stack)-1])
			self.resizeEvent()
			self.showInspector(1)
		
		elif event.button()==Qt.LeftButton:
			if self.mmode==0:
				self.mpressx = event.x()
				self.mpressy = event.y()
				
				# this is just a way of duplicating the last copy
				tmp = self.cam.t3d_stack.pop()
				t3d = Transform3D(tmp)
				self.cam.t3d_stack.append(tmp)
				self.cam.t3d_stack.append(t3d)
				self.updateInspector(t3d)
				
				
				return
		elif event.button()==Qt.RightButton:
			if self.mmode==0:
				self.mpressx = event.x()
				self.mpressy = event.y()
				return
		
	def mouseMoveEvent(self, event):
		if event.buttons()&Qt.LeftButton:
			if self.mmode==0:
				if event.modifiers() == Qt.ControlModifier:
					self.motionTranslate(event.x()-self.mpressx, self.mpressy - event.y())
				else:
					self.motionRotate(self.mpressx - event.x(), self.mpressy - event.y())
				self.parent.updateGL()
				self.mpressx = event.x()
				self.mpressy = event.y()
				return
		if event.buttons()&Qt.RightButton:
			if self.mmode==0:
				if event.modifiers() == Qt.ControlModifier:
					self.scale_event(event.y()-self.mpressy)	
				else:
					self.motionTranslate(event.x()-self.mpressx, self.mpressy - event.y())
					
				self.mpressx = event.x()
				self.mpressy = event.y()
				self.parent.updateGL()
				return
	
	def mouseReleaseEvent(self, event):
		if event.button()==Qt.LeftButton:
			if self.mmode==0:
				return
		elif event.button()==Qt.RightButton:
			if self.mmode==0:
				return
			
	def wheelEvent(self, event):
		self.scale_event(event.delta())
		self.resizeEvent()
		self.parent.updateGL()
		
	def scale_event(self,delta):
		self.cam.scale_event(delta)
		if self.inspector: self.inspector.setScale(self.cam.scale)

	def getTranslateScale(self):
	
		[rx,ry] = self.parent.get_render_dims_at_depth(self.cam.cam_z)
		
		#print "render area is %f %f " %(xx,yy)
		xscale = rx/float(self.parent.width())
		yscale = ry/float(self.parent.height())
		
		return [xscale,yscale]
	
	def motionTranslate(self,x,y):
		[xscale,yscale] = self.getTranslateScale()
		self.cam.cam_x += x*xscale
		self.cam.cam_y += y*yscale
		self.inspector.setXYTrans(self.cam.cam_x, self.cam.cam_y)
		
	def setCamZ(self,z):
		self.cam.setCamZ( z )
		self.parent.updateGL()
		
	def setCamY(self,y):
		self.cam.setCamY( y )
		self.parent.updateGL()
		
	def setCamX(self,x):
		self.cam.setCamX( x )
		self.parent.updateGL()
		
	def motionRotate(self,x,y):
		self.cam.motionRotate(x,y)
		size = len(self.cam.t3d_stack)
		self.updateInspector(self.cam.t3d_stack[size-1])
		
	def setScale(self,val):
		self.cam.scale = val
		self.parent.updateGL()
	
	def loadRotation(self,t3d):
		self.cam.t3d_stack.append(t3d)
		self.parent.updateGL()
		
	def resizeEvent(self):
		if self.inspector == None: return
		[xscale,yscale] = self.getTranslateScale()
		if ( xscale > yscale ): self.inspector.setTranslateScale(xscale,yscale,yscale)
		else: self.inspector.setTranslateScale(xscale,yscale,xscale)
		
	def draw_bc_screen(self):
		if (self.glcontrast == 1 and self.glbrightness == 0 ): return
		
		lighting = glIsEnabled(GL_LIGHTING)
		cull = glIsEnabled(GL_CULL_FACE)
		depth = glIsEnabled(GL_DEPTH_TEST)
		blend = glIsEnabled(GL_BLEND)
		
		polygonmode = glGetIntegerv(GL_POLYGON_MODE)

		glDisable(GL_LIGHTING)
		glDisable(GL_CULL_FACE)
		glDisable(GL_DEPTH_TEST)
		

		glEnable(GL_BLEND)
		glDepthMask(GL_FALSE)
		if ( self.glcontrast > 1 ):
			glBlendFunc(GL_ONE, GL_ONE)
			if self.glbrightness > 0 :
				glBlendEquation(GL_FUNC_ADD);
				glColor4f(self.glbrightness,self.glbrightness,self.glbrightness,1.0)
			else:
				glBlendEquation(GL_FUNC_REVERSE_SUBTRACT);
				glColor4f(-self.glbrightness,-self.glbrightness,-self.glbrightness, 1.0)
			
			glBegin( GL_QUADS )
			glVertex(0, 0)
			glVertex(1, 0)
			glVertex2f(1, 1)
			glVertex2f(0, 1)
			glEnd()
		
			glBlendFunc(GL_DST_COLOR, GL_ONE)
			glBlendEquation(GL_FUNC_ADD)
			
			tmpContrast = self.glcontrast
	
			while ( tmpContrast > 2 ):
				glColor4f(1.0,1.0,1.0,1.0)
				glBegin( GL_QUADS );
				glVertex2f(0, 0)
				glVertex2f(1, 0)
				glVertex2f(1, 1)
				glVertex2f(0, 1)
				glEnd()
				tmpContrast /= 2;
			
	
			glBlendFunc(GL_DST_COLOR, GL_ONE)
			glBlendEquation(GL_FUNC_ADD)
			glColor4f(tmpContrast-1.0,tmpContrast-1.0,tmpContrast-1.0,1.0)
			glBegin( GL_QUADS )
			glVertex2f(0, 0)
			glVertex2f(1, 0)
			glVertex2f(1, 1)
			glVertex2f(0, 1)
			glEnd()
		else:
			if self.glbrightness > 0:
				glBlendEquation(GL_FUNC_ADD)
				glColor4f(self.glbrightness,self.glbrightness,self.glbrightness,self.glcontrast)
			else:
				glBlendEquation(GL_FUNC_REVERSE_SUBTRACT);
				glColor4f(-self.glbrightness,-self.glbrightness,-self.glbrightness,self.glcontrast)
				
			glBlendFunc(GL_ONE, GL_SRC_ALPHA)

			glBegin( GL_QUADS )
			glVertex2f(0, 0)
			glVertex2f(1, 0)
			glVertex2f(1, 1)
			glVertex2f(0, 1)
			glEnd()
		
		glDepthMask(GL_TRUE)
	
		if ( lighting ): glEnable(GL_LIGHTING)
		if ( cull ): glEnable(GL_CULL_FACE)
		if ( depth ): glEnable(GL_DEPTH_TEST)
		if ( not blend ): glDisable(GL_BLEND)
		
		if ( polygonmode[0] == GL_LINE ): glPolygonMode(GL_FRONT, GL_LINE)
		if ( polygonmode[1] == GL_LINE ): glPolygonMode(GL_BACK, GL_LINE)
	
	def setGLBrightness(self,val):
		self.glbrightness = val
		self.parent.updateGL()
		
	def setGLContrast(self,val):
		self.glcontrast = val
		self.parent.updateGL()
		
	def draw_volume_bounds(self):
		# FIXME - should be a display list
		width = self.data.get_xsize()
		height = self.data.get_ysize()
		depth = self.data.get_zsize()
		glTranslate(-width/2.0,-height/2.0,-depth/2.0)
		glLineWidth(0.2)
		glNormal(0,1,0)
		glColor(.2,.1,0.4,1.0)
		glColor(1,1,1,1.0)
		glMaterial(GL_FRONT, GL_AMBIENT, [1, 1, 1,1.0])
		glMaterial(GL_FRONT, GL_DIFFUSE, [1, 1, 1,1.0])
		glMaterial(GL_FRONT, GL_SPECULAR, [0.774597, 0.774597, 0.774597,1.0])
		glMaterial(GL_FRONT, GL_SHININESS, 128.0)

		glBegin(GL_LINE_STRIP)
		glVertex(0,0,0)
		glVertex(width,0,0)
		glVertex(width,0,depth)
		glVertex(0,0,depth)
		glVertex(0,0,0)
		glVertex(0,height,0)
		glVertex(width,height,0)
		glVertex(width,height,depth)
		glVertex(0,height,depth)
		glVertex(0,height,0)
		glEnd()
		
		glBegin(GL_LINES)
		glVertex(width,height,depth)
		glVertex(width,0,depth)
		
		glVertex(width,height,0)
		glVertex(width,0,0)
		
		glVertex(0,height,depth)
		glVertex(0,0,depth)
		glEnd()
	
	def toggleCube(self):
		self.cube = not self.cube
		self.parent.updateGL()
		
	def showInspector(self,force=0):
		if not force and self.inspector==None : return
		
		if not self.inspector : self.inspector=EMVolumeInspector(self)
		self.inspector.show()
			
	
	def closeEvent(self,event) :
		if self.inspector: self.inspector.close()