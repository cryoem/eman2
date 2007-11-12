#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
# and David Woolford 10/26/2007 (woolford@bcm.edu)
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



from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from OpenGL import GL,GLU,GLUT
from OpenGL.GL import *
from OpenGL.GLU import *
from valslider import ValSlider
from math import *
from EMAN2 import *
import sys
import numpy
from emimageutil import ImgHistogram
from weakref import WeakKeyDictionary
from time import time
from PyQt4.QtCore import QTimer

from time import *

from emimage3dobject import EMImage3DObject

MAG_INCREMENT_FACTOR = 1.1


class EMVolume(EMImage3DObject):
	def __init__(self,image=None, parent=None):
		self.parent = parent
		
		self.init()
		self.initialized = True
		
		self.initializedGL= False
		
		self.inspector=None
		
		if image :
			self.setData(image)
	
	def getType(self):
		return "volume"

	def init(self):
		self.data=None
		
		self.aspect=1.0
		self.mmode=0
		
		self.scale = 1.0
		self.cam_x = 0
		self.cam_y = 0
		self.cam_z = 0
		self.cube = False
		
		t3d = Transform3D()
		rot = {}
		rot["az"] = 0.0
		rot["alt"] = 0.0
		rot["phi"] = 0.0
		
		t3d.set_rotation( EULER_EMAN, rot )
		self.t3d_stack = []
		self.t3d_stack.append(t3d)
		
		self.contrast = 1.0
		self.brightness = 0.0
		self.glcontrast = 1.0
		self.glbrightness = 0.0
		self.texsample = 1.0
		
		self.loadColors()
		
		self.tex_name = 0
		
		self.tex_dl = 0
		self.tex_dl_x = 0
		self.tex_dl_y = 0
		self.tex_dl_z = 0
		self.inspector=None
		
	def loadColors(self):
		# There is a bug, if the accompanying functionality to this function is removed we get a seg fault
		ruby = {}
		ruby["ambient"] = [0.1745, 0.01175, 0.01175,1.0]
		ruby["diffuse"] = [0.61424, 0.04136, 0.04136,1.0]
		ruby["specular"] = [0.927811, 0.826959, 0.826959,1.0]
		ruby["shininess"] = 128.0
		self.colors = {}
		self.colors["ruby"] = ruby
		self.isocolor = "ruby"
		
	def setData(self,data):
		"""Pass in a 3D EMData object"""
		
		self.data=data
		if data==None:
			print "Error, the data is empty"
			return
		
	 	if (isinstance(data,EMData) and data.get_zsize()<=1) :
			print "Error, the data is not 3D"
			return
		
		self.default_z = -1.25*data.get_zsize()
		self.cam_z = self.default_z
		
		if not self.inspector or self.inspector ==None:
			self.inspector=EMVolumeInspector(self)
		
		self.inspector.setColors(self.colors,self.isocolor)
		
		self.updateDataAndTexture()
		self.tex_dl = self.tex_dl_z
		
	def render(self):
		lighting = glIsEnabled(GL_LIGHTING)
		cull = glIsEnabled(GL_CULL_FACE)
		depth = glIsEnabled(GL_DEPTH_TEST)
		
		polygonmode = glGetIntegerv(GL_POLYGON_MODE)

		glDisable(GL_LIGHTING)
		glDisable(GL_CULL_FACE)
		glDisable(GL_DEPTH_TEST)
		
		glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
		
		glTranslated(self.cam_x, self.cam_y, self.cam_z)
		
		rot = self.t3d_stack[len(self.t3d_stack)-1].get_rotation()
		glRotate(float(rot["phi"]),0,0,1)
		glRotate(float(rot["alt"]),1,0,0)
		glRotate(float(rot["az"]),0,0,1)
		
		# here is where zoom is applied
		glScalef(self.scale,self.scale,self.scale)

		if ( self.tex_dl == 0 ):
			self.updateDataAndTexture()
		
		# here is where the correct display list (x,y or z direction) is determined
		self.determineTextureView()

		glPushMatrix()
		glTranslate(-self.data.get_xsize()/2.0,-self.data.get_ysize()/2.0,-self.data.get_zsize()/2.0)
		glScalef(self.data.get_xsize(),self.data.get_ysize(),self.data.get_zsize())
		glEnable(GL_BLEND)
		glBlendEquation(GL_FUNC_ADD)
		glDepthMask(GL_FALSE)
		glBlendFunc(GL_ONE, GL_ONE)
		glCallList(self.tex_dl)
		glDepthMask(GL_TRUE)
		glDisable(GL_BLEND)
		glPopMatrix()
	
		glPushMatrix()
		glLoadIdentity()
		glTranslate(-self.data.get_xsize()/2.0,-self.data.get_ysize()/2.0,-10)
		glScalef(self.data.get_xsize(),self.data.get_ysize(),1)
		self.draw_bc_screen()
		glPopMatrix()
		
		if self.cube:
			glPushMatrix()
			self.draw_volume_bounds()
			glPopMatrix()
			
		if ( lighting ): glEnable(GL_LIGHTING)
		if ( cull ): glEnable(GL_CULL_FACE)
		if ( depth ): glEnable(GL_DEPTH_TEST)
		
		if ( polygonmode[0] == GL_LINE ): glPolygonMode(GL_FRONT, GL_LINE)
		if ( polygonmode[1] == GL_LINE ): glPolygonMode(GL_BACK, GL_LINE)
		
	
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
	
	def determineTextureView(self):
		t3d = self.t3d_stack[len(self.t3d_stack)-1]
		
		point = Vec3f(0,0,1)
		
		point = point*t3d
		
		point = [abs(point[0]), abs(point[1]),abs(point[2])]
		
		xyangle = 180*atan2(point[1],point[0])/pi
		xzangle = 180*atan2(point[2],point[0])/pi
		yzangle = 180*atan2(point[2],point[1])/pi
		
		if (xzangle > 45 ):
			
			if ( yzangle > 45 ):
				self.tex_dl = self.tex_dl_z
			else:
				self.tex_dl = self.tex_dl_y
		else:
			if ( xyangle < 45 ):
				self.tex_dl = self.tex_dl_x
			else:
				self.tex_dl = self.tex_dl_y
	
	def genTexture(self):
	
		if ( self.tex_name != 0 ): glDeleteTextures(self.tex_name)
		
		# Get the texture here
		self.tex_name = self.data_copy.gen_gl_texture()
		
		if ( self.tex_dl_z != 0 ): glDeleteLists( self.tex_dl_z, 1)
		
		self.tex_dl_z = glGenLists(1)
		
		glNewList(self.tex_dl_z,GL_COMPILE)
		glEnable(GL_TEXTURE_3D)
		glBindTexture(GL_TEXTURE_3D, self.tex_name)
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1)
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP)
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP)
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP)
		glTexParameter(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
		glTexParameter(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE)
		glBegin(GL_QUADS)
		for z in range(0,int(self.texsample*(self.data.get_zsize()))):
			zz = float(z)/float(self.data.get_zsize()-1)/self.texsample
			glTexCoord3f(0,0,zz)
			glVertex(0,0,zz)
			
			glTexCoord3f(1,0,zz)
			glVertex(1,0,zz)
			
			glTexCoord3f(1,1,zz)
			glVertex(1,1,zz)
			
			glTexCoord3f(0,1,zz)
			glVertex(0,1,zz)
		
		glEnd()
		
		glDisable(GL_TEXTURE_3D)
		glEndList()
		
		if ( self.tex_dl_y != 0 ): glDeleteLists( self.tex_dl_y, 1)
		
		self.tex_dl_y = glGenLists(1)
		
		glNewList(self.tex_dl_y,GL_COMPILE)
		glEnable(GL_TEXTURE_3D)
		glBindTexture(GL_TEXTURE_3D, self.tex_name)
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1)
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP)
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP)
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP)
		glTexParameter(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
		glTexParameter(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE)
		glBegin(GL_QUADS)
		for y in range(0,int(self.texsample*(self.data.get_ysize()))):
			yy = float(y)/float(self.data.get_ysize()-1)/self.texsample
			glTexCoord3f(0,yy,0)
			glVertex(0,yy,0)
			
			glTexCoord3f(1,yy,0)
			glVertex(1,yy,0)
			
			glTexCoord3f(1,yy,1)
			glVertex(1,yy,1)
			
			glTexCoord3f(0,yy,1)
			glVertex(0,yy,1)
		
		glEnd()
		
		glDisable(GL_TEXTURE_3D)
		glEndList()
		
		if ( self.tex_dl_x != 0 ): glDeleteLists( self.tex_dl_x, 1)
		
		self.tex_dl_x = glGenLists(1)
		
		glNewList(self.tex_dl_x,GL_COMPILE)
		glEnable(GL_TEXTURE_3D)
		glBindTexture(GL_TEXTURE_3D, self.tex_name)
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1)
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP)
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP)
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP)
		glTexParameter(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
		glTexParameter(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE)
		glBegin(GL_QUADS)
		for x in range(0,int(self.texsample*(self.data.get_xsize()))):
			xx = float(x)/float(self.data.get_xsize()-1)/self.texsample
			glTexCoord3f(xx,0,0)
			glVertex(xx,0,0)
			
			glTexCoord3f(xx,1,0)
			glVertex(xx,1,0)
			
			glTexCoord3f(xx,1,1)
			glVertex(xx,1,1)
			
			glTexCoord3f(xx,0,1)
			glVertex(xx,0,1)
		
		glEnd()
		
		glDisable(GL_TEXTURE_3D)
		glEndList()
	
	def updateDataAndTexture(self):
	
		if ( not isinstance(self.data,EMData) ): return
		
		self.data_copy = self.data.copy()
		self.data_copy.mult(self.contrast*1.0/self.data.get_zsize())
		self.data_copy.add(self.brightness)

		hist = self.data_copy.calc_hist(256,0,1.0)
		self.inspector.setHist(hist,0,1.0) 

		self.genTexture()
	
	def draw_bc_screen(self):
		if (self.glcontrast == 1 and self.glbrightness == 0 ): return
		
		depth_testing_was_on = glIsEnabled(GL_DEPTH_TEST)

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
		
		glDisable(GL_BLEND)
		if depth_testing_was_on:
			glDepthMask(GL_TRUE)
	
	def showInspector(self,force=0):
		if not force and self.inspector==None : return
		
		if not self.inspector : self.inspector=EMVolumeInspector(self)
		self.inspector.show()
			
	def updateInspector(self,t3d):
		if not self.inspector or self.inspector ==None:
			self.inspector=EMVolumeInspector(self)
		self.inspector.updateRotations(t3d)
	
	def getInspector(self):
		if not self.inspector : self.inspector=EMIsoInspector(self)
		return self.inspector
	
	def closeEvent(self,event) :
		if self.inspector: self.inspector.close()
		
	def mousePressEvent(self, event):
#		lc=self.scrtoimg((event.x(),event.y()))
		if event.button()==Qt.MidButton:
			if not self.inspector or self.inspector ==None:
				self.inspector=EMImageInspector3D(self)
			self.inspector.updateRotations(self.t3d_stack[len(self.t3d_stack)-1])
			self.resizeEvent()
			self.showInspector(1)
		
		elif event.button()==Qt.LeftButton:
			if self.mmode==0:
				self.mpressx = event.x()
				self.mpressy = event.y()
				
				# this is just a way of duplicating the last copy
				tmp = self.t3d_stack.pop()
				t3d = Transform3D(tmp)
				self.t3d_stack.append(tmp)
				self.t3d_stack.append(t3d)
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

	def resizeEvent(self):
		[xscale,yscale] = self.getTranslateScale()
		if ( xscale > yscale ): self.inspector.setTranslateScale(xscale,yscale,yscale)
		else: self.inspector.setTranslateScale(xscale,yscale,xscale)

	def scale_event(self,delta):
		if delta > 0:
			self.scale *= MAG_INCREMENT_FACTOR
		elif delta < 0:
			self.scale *= 1.0/MAG_INCREMENT_FACTOR
		# The self.scale variable is updated now, so just update with that
		if self.inspector: self.inspector.setScale(self.scale)

	def getTranslateScale(self):
	
		[rx,ry] = self.parent.get_render_dims_at_depth(self.cam_z)
		
		#print "render area is %f %f " %(xx,yy)
		xscale = rx/float(self.parent.width())
		yscale = ry/float(self.parent.height())
		
		return [xscale,yscale]

	def motionTranslate(self,x,y):
		[xscale,yscale] = self.getTranslateScale()
		
		self.cam_x += x*xscale
		self.cam_y += y*yscale
		
		self.inspector.setXYTrans(self.cam_x, self.cam_y)
		
	def setCamZ(self,z):
		self.cam_z = self.default_z + z
		self.parent.updateGL()
		
	def setCamY(self,y):
		self.cam_y = y
		self.parent.updateGL()
		
	def setCamX(self,x):
		self.cam_x = x
		self.parent.updateGL()
		
	def motionRotate(self,x,y):
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
		self.updateInspector(self.t3d_stack[size-1])
		
	def setScale(self,val):
		self.scale = val
		self.parent.updateGL()
	
	def loadRotation(self,t3d):
		self.t3d_stack.append(t3d)
		self.parent.updateGL()
		
	def setColor(self,val):
		#print val
		#self.isocolor = str(val)
		self.parent.updateGL()
		
	def toggleCube(self):
		self.cube = not self.cube
		self.parent.updateGL()
		
	def setContrast(self,val):
		self.contrast = val
		self.updateDataAndTexture()
		self.parent.updateGL()
		
	def setBrightness(self,val):
		self.brightness = val/self.data.get_zsize()
		self.updateDataAndTexture()
		self.parent.updateGL()
		
	def setGLBrightness(self,val):
		self.glbrightness = val
		self.parent.updateGL()
		
	def setGLContrast(self,val):
		self.glcontrast = val
		self.parent.updateGL()
		
	def setTextureSample(self,val):
		if ( val < 0 ) :
			print "Error, cannot handle texture sample less than 0"
			return
		
		self.texsample = val
		self.genTexture()
		self.parent.updateGL()

class EMVolumeWidget(QtOpenGL.QGLWidget):
	
	allim=WeakKeyDictionary()
	def __init__(self, image=None, parent=None):
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True)
		fmt.setDepth(1)
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		EMVolumeWidget.allim[self]=0
		
		self.fov = 50 # field of view angle used by gluPerspective
		
		self.volume = EMVolume(image,self)
	def setData(self,data):
		self.volume.setData(data)
	def initializeGL(self):
		
		glEnable(GL_NORMALIZE)
		glEnable(GL_LIGHT0)
		glEnable(GL_DEPTH_TEST)
		#print "Initializing"
		glLightfv(GL_LIGHT0, GL_AMBIENT, [0.9, 0.9, 0.9, 1.0])
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_POSITION, [0.5,0.7,11.,0.])
		GL.glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST)
		
		GL.glClearColor(0,0,0,0)
	
		glShadeModel(GL_SMOOTH)
		
	def paintGL(self):
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		glPushMatrix()
		self.volume.render()
		glPopMatrix()
	
		
	def resizeGL(self, width, height):
		# just use the whole window for rendering
		glViewport(0,0,self.width(),self.height())
		
		# maintain the aspect ratio of the window we have
		self.aspect = float(width)/float(height)
		
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		# using gluPerspective for simplicity
		gluPerspective(self.fov,self.aspect,0.001,1000000)
		
		# switch back to model view mode
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		self.volume.resizeEvent()

	def showInspector(self,force=0):
		self.volume.showInspector(self,force)

	def closeEvent(self,event) :
		self.volume.closeEvent(event)
		
	def mousePressEvent(self, event):
		self.volume.mousePressEvent(event)
		self.emit(QtCore.SIGNAL("mousedown"), event)
		
	def mouseMoveEvent(self, event):
		self.volume.mouseMoveEvent(event)
		self.emit(QtCore.SIGNAL("mousedrag"), event)
	
	def mouseReleaseEvent(self, event):
		self.volume.mouseReleaseEvent(event)
		self.emit(QtCore.SIGNAL("mouseup"), event)
			
	def wheelEvent(self, event):
		self.volume.wheelEvent(event)

	def get_render_dims_at_depth(self,depth):
		# This function returns the width and height of the renderable 
		# area at the origin of the data volume
		height = -2*tan(self.fov/2.0*pi/180.0)*(depth)
		width = self.aspect*height
		
		return [width,height]
		

class EMVolumeInspector(QtGui.QWidget):
	def __init__(self,target) :
		QtGui.QWidget.__init__(self,None)
		self.target=target
		
		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		
		self.hbl = QtGui.QHBoxLayout()
		self.hbl.setMargin(0)
		self.hbl.setSpacing(6)
		self.hbl.setObjectName("hbl")
		self.vbl.addLayout(self.hbl)
		
		self.hist = ImgHistogram(self)
		self.hist.setObjectName("hist")
		self.hbl.addWidget(self.hist)
		
		self.vbl2 = QtGui.QVBoxLayout()
		self.vbl2.setMargin(0)
		self.vbl2.setSpacing(6)
		self.vbl2.setObjectName("vbl2")
		self.hbl.addLayout(self.vbl2)
	
		self.cubetog = QtGui.QPushButton("Cube")
		self.cubetog.setCheckable(1)
		self.vbl2.addWidget(self.cubetog)
		
		self.defaults = QtGui.QPushButton("Defaults")
		self.vbl2.addWidget(self.defaults)
		
		self.tabwidget = QtGui.QTabWidget()
		
		self.tabwidget.addTab(self.getMainTab(), "Main")
		self.tabwidget.addTab(self.getGLTab(),"GL")
		
		self.vbl.addWidget(self.tabwidget)
		
		self.n3_showing = False
		
		self.current_src = EULER_EMAN
		
		QtCore.QObject.connect(self.scale, QtCore.SIGNAL("valueChanged"), target.setScale)
		QtCore.QObject.connect(self.contrast, QtCore.SIGNAL("valueChanged"), target.setContrast)
		QtCore.QObject.connect(self.glcontrast, QtCore.SIGNAL("valueChanged"), target.setGLContrast)
		QtCore.QObject.connect(self.glbrightness, QtCore.SIGNAL("valueChanged"), target.setGLBrightness)
		QtCore.QObject.connect(self.bright, QtCore.SIGNAL("valueChanged"), target.setBrightness)
		QtCore.QObject.connect(self.az, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
		QtCore.QObject.connect(self.alt, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
		QtCore.QObject.connect(self.phi, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
		QtCore.QObject.connect(self.cbb, QtCore.SIGNAL("currentIndexChanged(QString)"), target.setColor)
		QtCore.QObject.connect(self.src, QtCore.SIGNAL("currentIndexChanged(QString)"), self.set_src)
		QtCore.QObject.connect(self.x_trans, QtCore.SIGNAL("valueChanged(double)"), target.setCamX)
		QtCore.QObject.connect(self.y_trans, QtCore.SIGNAL("valueChanged(double)"), target.setCamY)
		QtCore.QObject.connect(self.z_trans, QtCore.SIGNAL("valueChanged(double)"), target.setCamZ)
		QtCore.QObject.connect(self.cubetog, QtCore.SIGNAL("toggled(bool)"), target.toggleCube)
		QtCore.QObject.connect(self.defaults, QtCore.SIGNAL("clicked(bool)"), self.setDefaults)
		QtCore.QObject.connect(self.smp, QtCore.SIGNAL("valueChanged(int)"), target.setTextureSample)
	
	def getGLTab(self):
		self.gltab = QtGui.QWidget()
		gltab = self.gltab
		
		gltab.vbl = QtGui.QVBoxLayout(self.gltab )
		gltab.vbl.setMargin(0)
		gltab.vbl.setSpacing(6)
		gltab.vbl.setObjectName("Main")
		
		self.glcontrast = ValSlider(self,(0.0,20.0),"GLShd:")
		self.glcontrast.setObjectName("GLShade")
		self.glcontrast.setValue(1.0)
		gltab.vbl.addWidget(self.glcontrast)
		
		self.glbrightness = ValSlider(self,(-1.0,1.0),"GLBst:")
		self.glbrightness.setObjectName("GLBoost")
		self.glbrightness.setValue(0.1)
		self.glbrightness.setValue(0.0)
		gltab.vbl.addWidget(self.glbrightness)
	
		return gltab
	
	def getMainTab(self):
	
		self.maintab = QtGui.QWidget()
		maintab = self.maintab
		maintab.vbl = QtGui.QVBoxLayout(self.maintab)
		maintab.vbl.setMargin(0)
		maintab.vbl.setSpacing(6)
		maintab.vbl.setObjectName("Main")
		
		
		self.scale = ValSlider(self,(0.01,30.0),"Zoom:")
		self.scale.setObjectName("scale")
		self.scale.setValue(1.0)
		maintab.vbl.addWidget(self.scale)
		
		self.contrast = ValSlider(maintab,(0.0,20.0),"Cont:")
		self.contrast.setObjectName("contrast")
		self.contrast.setValue(1.0)
		maintab.vbl.addWidget(self.contrast)

		self.bright = ValSlider(maintab,(-5.0,5.0),"Brt:")
		self.bright.setObjectName("bright")
		self.bright.setValue(0.1)
		self.bright.setValue(0.0)
		maintab.vbl.addWidget(self.bright)

		self.hbl_smp = QtGui.QHBoxLayout()
		self.hbl_smp.setMargin(0)
		self.hbl_smp.setSpacing(6)
		self.hbl_smp.setObjectName("Texture Oversampling")
		maintab.vbl.addLayout(self.hbl_smp)
		
		self.smp_label = QtGui.QLabel()
		self.smp_label.setText('Texture Oversampling')
		self.hbl_smp.addWidget(self.smp_label)
		
		self.smp = QtGui.QSpinBox(maintab)
		self.smp.setMaximum(10)
		self.smp.setMinimum(1)
		self.smp.setValue(1)
		self.hbl_smp.addWidget(self.smp)

		self.lowlim=0
		self.highlim=1.0
		self.busy=0
		
		self.cbb = QtGui.QComboBox(maintab)
		self.vbl.addWidget(self.cbb)
		self.cbb.deleteLater()

		self.hbl_trans = QtGui.QHBoxLayout()
		self.hbl_trans.setMargin(0)
		self.hbl_trans.setSpacing(6)
		self.hbl_trans.setObjectName("Trans")
		maintab.vbl.addLayout(self.hbl_trans)
		
		self.x_label = QtGui.QLabel()
		self.x_label.setText('x')
		self.hbl_trans.addWidget(self.x_label)
		
		self.x_trans = QtGui.QDoubleSpinBox(maintab)
		self.x_trans.setMinimum(-10000)
		self.x_trans.setMaximum(10000)
		self.x_trans.setValue(0.0)
		self.hbl_trans.addWidget(self.x_trans)
		
		self.y_label = QtGui.QLabel()
		self.y_label.setText('y')
		self.hbl_trans.addWidget(self.y_label)
		
		self.y_trans = QtGui.QDoubleSpinBox(maintab)
		self.y_trans.setMinimum(-10000)
		self.y_trans.setMaximum(10000)
		self.y_trans.setValue(0.0)
		self.hbl_trans.addWidget(self.y_trans)
		
		self.z_label = QtGui.QLabel()
		self.z_label.setText('z')
		self.hbl_trans.addWidget(self.z_label)
		
		self.z_trans = QtGui.QDoubleSpinBox(maintab)
		self.z_trans.setMinimum(-10000)
		self.z_trans.setMaximum(10000)
		self.z_trans.setValue(0.0)
		self.hbl_trans.addWidget(self.z_trans)
		
		self.hbl_src = QtGui.QHBoxLayout()
		self.hbl_src.setMargin(0)
		self.hbl_src.setSpacing(6)
		self.hbl_src.setObjectName("hbl")
		maintab.vbl.addLayout(self.hbl_src)
		
		self.label_src = QtGui.QLabel()
		self.label_src.setText('Rotation Convention')
		self.hbl_src.addWidget(self.label_src)
		
		self.src = QtGui.QComboBox(maintab)
		self.load_src_options(self.src)
		self.hbl_src.addWidget(self.src)
		
		# set default value -1 ensures that the val slider is updated the first time it is created
		self.az = ValSlider(maintab,(-360.0,360.0),"az",-1)
		self.az.setObjectName("az")
		maintab.vbl.addWidget(self.az)
		
		self.alt = ValSlider(maintab,(-180.0,180.0),"alt",-1)
		self.alt.setObjectName("alt")
		maintab.vbl.addWidget(self.alt)
		
		self.phi = ValSlider(maintab,(-360.0,360.0),"phi",-1)
		self.phi.setObjectName("phi")
		maintab.vbl.addWidget(self.phi)
		
		return maintab
	
	def setDefaults(self):
		self.x_trans.setValue(0.0)
		self.y_trans.setValue(0.0)
		self.z_trans.setValue(0.0)
		self.scale.setValue(1.0)
		self.contrast.setValue(1.0)
		self.bright.setValue(0.0)
		self.glcontrast.setValue(1.0)
		self.glbrightness.setValue(0.0)
		
		self.az.setValue(0.0)
		self.alt.setValue(0.0)
		self.phi.setValue(0.0)

	def setXYTrans(self, x, y):
		self.x_trans.setValue(x)
		self.y_trans.setValue(y)
	
	def setTranslateScale(self, xscale,yscale,zscale):
		self.x_trans.setSingleStep(xscale)
		self.y_trans.setSingleStep(yscale)
		self.z_trans.setSingleStep(zscale)

	def updateRotations(self,t3d):
		rot = t3d.get_rotation(self.src_map[str(self.src.itemText(self.src.currentIndex()))])
		
		convention = self.src.currentText()
		if ( self.src_map[str(convention)] == EULER_SPIN ):
			self.n3.setValue(rot[self.n3.getLabel()],True)
		
		self.az.setValue(rot[self.az.getLabel()],True)
		self.alt.setValue(rot[self.alt.getLabel()],True)
		self.phi.setValue(rot[self.phi.getLabel()],True)
	
	def sliderRotate(self):
		self.target.loadRotation(self.getCurrentRotation())
	
	def getCurrentRotation(self):
		convention = self.src.currentText()
		rot = {}
		if ( self.current_src == EULER_SPIN ):
			rot[self.az.getLabel()] = self.az.getValue()
			
			n1 = self.alt.getValue()
			n2 = self.phi.getValue()
			n3 = self.n3.getValue()
			
			norm = sqrt(n1*n1 + n2*n2 + n3*n3)
			
			n1 /= norm
			n2 /= norm
			n3 /= norm
			
			rot[self.alt.getLabel()] = n1
			rot[self.phi.getLabel()] = n2
			rot[self.n3.getLabel()] = n3
			
		else:
			rot[self.az.getLabel()] = self.az.getValue()
			rot[self.alt.getLabel()] = self.alt.getValue()
			rot[self.phi.getLabel()] = self.phi.getValue()
		
		return Transform3D(self.current_src, rot)
	
	def set_src(self, val):
		t3d = self.getCurrentRotation()
		
		if (self.n3_showing) :
			self.vbl.removeWidget(self.n3)
			self.n3.deleteLater()
			self.n3_showing = False
			self.az.setRange(-360,360)
			self.alt.setRange(-180,180)
			self.phi.setRange(-360,660)
		
		if ( self.src_map[str(val)] == EULER_SPIDER ):
			self.az.setLabel('phi')
			self.alt.setLabel('theta')
			self.phi.setLabel('psi')
		elif ( self.src_map[str(val)] == EULER_EMAN ):
			self.az.setLabel('az')
			self.alt.setLabel('alt')
			self.phi.setLabel('phi')
		elif ( self.src_map[str(val)] == EULER_IMAGIC ):
			self.az.setLabel('alpha')
			self.alt.setLabel('beta')
			self.phi.setLabel('gamma')
		elif ( self.src_map[str(val)] == EULER_XYZ ):
			self.az.setLabel('xtilt')
			self.alt.setLabel('ytilt')
			self.phi.setLabel('ztilt')
		elif ( self.src_map[str(val)] == EULER_MRC ):
			self.az.setLabel('phi')
			self.alt.setLabel('theta')
			self.phi.setLabel('omega')
		elif ( self.src_map[str(val)] == EULER_SPIN ):
			self.az.setLabel('Omega')
			self.alt.setRange(-1,1)
			self.phi.setRange(-1,1)
			
			self.alt.setLabel('n1')
			self.phi.setLabel('n2')
			
			self.n3 = ValSlider(self,(-360.0,360.0),"n3",-1)
			self.n3.setRange(-1,1)
			self.n3.setObjectName("n3")
			self.vbl.addWidget(self.n3)
			QtCore.QObject.connect(self.n3, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
			self.n3_showing = True
		
		self.current_src = self.src_map[str(val)]
		self.updateRotations(t3d)
	
	def load_src_options(self,widgit):
		self.load_src()
		for i in self.src_strings:
			widgit.addItem(i)
	
	# read src as 'supported rotation conventions'
	def load_src(self):
		# supported_rot_conventions
		src_flags = []
		src_flags.append(EULER_EMAN)
		src_flags.append(EULER_SPIDER)
		src_flags.append(EULER_IMAGIC)
		src_flags.append(EULER_MRC)
		src_flags.append(EULER_SPIN)
		src_flags.append(EULER_XYZ)
		
		self.src_strings = []
		self.src_map = {}
		for i in src_flags:
			self.src_strings.append(str(i))
			self.src_map[str(i)] = i
		
	
	def setColors(self,colors,current_color):
		for i in colors:
			self.cbb.addItem(i)

	def setHist(self,hist,minden,maxden):
		self.hist.setData(hist,minden,maxden)

	def setScale(self,newscale):
		self.scale.setValue(newscale)
		
# This is just for testing, of course
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = EMVolumeWidget()
 	if len(sys.argv)==1 : 
		e = EMData()
		e.set_size(128,128,128)
		e.process_inplace('testimage.x')
		window.setData(e)

		# these lines are for testing shape rendering
# 		window.addShape("a",["rect",.2,.8,.2,20,20,80,80,2])
# 		window.addShape("b",["circle",.5,.8,.2,120,50,30.0,2])
# 		window.addShape("c",["line",.2,.8,.5,20,120,100,200,2])
# 		window.addShape("d",["label",.2,.8,.5,220,220,"Testing",14,1])
	else :
		if not os.path.exists(sys.argv[1]):
			print "Error, input file %s does not exist" %sys.argv[1]
			exit(1)
		a=EMData.read_images(sys.argv[1],[0])
		window.setData(a[0])
	window.show()
	
#	w2=QtGui.QWidget()
#	w2.resize(256,128)
	
#	w3=ValSlider(w2)
#	w3.resize(256,24)
#	w2.show()
	
	sys.exit(app.exec_())
