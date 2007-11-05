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

t3d_stack = []
MAG_INCREMENT_FACTOR = 1.1

class EMImage3D(QtOpenGL.QGLWidget):
	""" This class is not yet complete.
	A QT widget for rendering 3D EMData objects.
	"""
	allim=WeakKeyDictionary()
	def __init__(self, image=None, parent=None):
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True)
		fmt.setDepth(1)
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		EMImage3D.allim[self]=0
# 		try: 
# 			if EMImage2D.gq : pass
# 		except:
# 			EMImage2D.gq=GLU.gluNewQuadric()
# 			GLU.gluQuadricDrawStyle(EMImage2D.gq,GLU.GLU_FILL)
# 			GLU.gluQuadricNormals(EMImage2D.gq,GLU.GLU_SMOOTH)
# 			GLU.gluQuadricNormals(EMImage2D.gq,GLU.GLU_NONE)
# 			GLU.gluQuadricOrientation(EMImage2D.gq,GLU.GLU_OUTSIDE)
# 			GLU.gluQuadricOrientation(EMImage2D.gq,GLU.GLU_INSIDE)
# 			GLU.gluQuadricTexture(EMImage2D.gq,GL.GL_FALSE)
		
		self.timer = QTimer()
		QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.timeout)

		
		self.data=None
		
		self.aspect=1.0
		self.gq=0
		self.mmode=0
		self.isothr=0.5
		self.isorender=None
		self.isodl = 0
		self.smpval=-1
		self.griddl = 0
		self.scale = 1.0
		self.cam_x = 0
		self.cam_y = 0
		self.cam_z = 0
		self.cube = False
		
		self.wire = False
		
		self.fov = 50 # field of view angle used by gluPerspective
		self.generate_iso = True
		self.inspector=None
		#self.ctrl_down = False
		
		ruby = {}
		ruby["ambient"] = [0.1745, 0.01175, 0.01175,1.0]
		ruby["diffuse"] = [0.61424, 0.04136, 0.04136,1.0]
		ruby["specular"] = [0.927811, 0.826959, 0.826959,1.0]
		ruby["shininess"] = 128.0
		
		emerald = {}
		emerald["ambient"] = [0.0215, 0.1745, 0.0215,1.0]
		emerald["diffuse"] = [0.07568, 0.61424,  0.07568,1.0]
		emerald["specular"] = [0.833, 0.927811, 0.833,1.0]
		emerald["shininess"] = 128.0
		
		pearl = {}
		pearl["ambient"] = [0.25, 0.20725, 0.20725,1.0]
		pearl["diffuse"] = [1.0, 0.829, 0.829,1.0]
		pearl["specular"] = [0.296648, 0.296648, 0.296648,1.0]
		pearl["shininess"] = 128.0
		
		silver = {}
		silver["ambient"] = [0.25, 0.25, 0.25,1.0]
		silver["diffuse"] = [0.4, 0.4, 0.4,1.0]
		silver["specular"] = [0.974597, 0.974597, 0.974597,1.0]
		silver["shininess"] = 128.0
		
		gold = {}
		gold["ambient"] = [0.24725, 0.2245, 0.0645,1.0]
		gold["diffuse"] = [0.34615, 0.3143, 0.0903,1.0]
		gold["specular"] = [1.000, 0.9079885, 0.26086934,1.0]
		gold["shininess"] = 128.0
		
		copper = {}
		copper["ambient"] = [0.2295, 0.08825, 0.0275,1.0]
		copper["diffuse"] = [0.5508, 0.2118, 0.066,1.0]
		copper["specular"] = [0.9, 0.5, 0.2,1.0]
		copper["shininess"] = 128.0
		
		obsidian = {}
		obsidian["ambient"] = [0.05375,  0.05,     0.06625 ,1.0]
		obsidian["diffuse"] = [0.18275,  0.17,     0.22525,1.0]
		obsidian["specular"] = [0.66, 0.65, 0.69]
		obsidian["shininess"] = 128.0
		
		turquoise = {}
		turquoise["ambient"] = [0.1, 0.18725, 0.1745 ,1.0]
		turquoise["diffuse"] = [0.396, 0.74151, 0.69102,1.0]
		turquoise["specular"] = [0.297254, 0.30829, 0.306678]
		turquoise["shininess"] = 128.0
		
		yellow = {}
		yellow["ambient"] = [0.3, 0.3, 0.0,1]
		yellow["diffuse"] = [0.5, 0.5, 0.0,1]
		yellow["specular"] = [0.7, 0.7, 0.0,1]
		yellow["shininess"] =  60
		
		self.colors = {}
		self.colors["ruby"] = ruby
		self.colors["emerald"] = emerald
		self.colors["pearl"] = pearl
		self.colors["silver"] = silver
		self.colors["gold"] = gold
		self.colors["copper"] = copper
		self.colors["obsidian"] = obsidian
		self.colors["turquoise"] = turquoise
		self.colors["yellow"] = yellow
		
		self.isocolor = "ruby"
		
		if image :
			self.setData(image)
			self.show()
		
	def timeout(self):
		self.updateGL()
		
	
	def setData(self,data):
		"""Pass in a 3D EMData object"""
#		if not self.data and data: self.resize(data.get_xsize(),data.get_ysize())
		
		self.data=data
		if data==None or (isinstance(data,EMData) and data.get_zsize()<=1) :
			self.updateGL()
			return
		
		self.minden=data.get_attr("minimum")
		self.maxden=data.get_attr("maximum")
		mean=data.get_attr("mean")
		sigma=data.get_attr("sigma")

		self.default_z = -1.25*data.get_zsize()
		self.cam_z = self.default_z
		
		if not self.inspector or self.inspector ==None:
			self.inspector=EMImageInspector3D(self)
		
		hist = data.calc_hist(256,self.minden,self.maxden)
		self.inspector.setHist(hist,self.minden,self.maxden) 
	
		self.inspector.setThrs(self.minden,self.maxden,mean+3.0*sigma)
		self.isothr = mean+3.0*sigma
		
		self.isorender=MarchingCubes(data)
		self.inspector.setSamplingRange(self.isorender.get_sampling_range())
		
		self.inspector.setColors(self.colors,self.isocolor)
		
	def myAlternateFunction(a):
		return 0
		#self.timer.start(25)
		
	def initializeGL(self):
		
		glEnable(GL_NORMALIZE)
		
		glEnable(GL_LIGHTING)
		glEnable(GL_LIGHT0)
		glEnable(GL_DEPTH_TEST)
		#print "Initializing"
		glLightfv(GL_LIGHT0, GL_AMBIENT, [0.9, 0.9, 0.9, 1.0])
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_POSITION, [0.5,0.7,11.,0.])

		GL.glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST)
		
		GL.glClearColor(0,0,0,0)
		
		if not self.gq:
			self.gq=gluNewQuadric()
			gluQuadricDrawStyle(self.gq,GLU_FILL)
			gluQuadricNormals(self.gq,GLU_SMOOTH)
			gluQuadricOrientation(self.gq,GLU_OUTSIDE)
			gluQuadricTexture(self.gq,GL_FALSE)
		
		# Precompile a displaylist for the display volume border
		#self.volcubedl=glGenLists(1)
		#glNewList(self.volcubedl,GL_COMPILE)
		#glPushMatrix()
		#glColor(.7,.7,1.0)
##		glRotate(90.,1.,0.,0.)
		#glTranslate(-self.aspect-.01,1.01,-4.0)
		#gluCylinder(self.gq,.01,.01,15.0,12,2)
		#glTranslate(self.aspect*2.0+.02,0.0,0.0)
		#gluCylinder(self.gq,.01,.01,15.0,12,2)
		#glTranslate(0.0,-2.02,0.0)
		#gluCylinder(self.gq,.01,.01,15.0,12,2)
		#glTranslate(-self.aspect*2.0-.02,0.0,0.0)
		#gluCylinder(self.gq,.01,.01,15.0,12,2)		
		#glPopMatrix()
		#glEndList()
		
		#cl = 50
		#self.cylinderdl = glGenLists(1)
		#glNewList(self.cylinderdl,GL_COMPILE)
		#glColor(.7,.7,1.0)
		#glMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, [.2,.1,0.4,1.0])
		#glMaterial(GL_FRONT, GL_SPECULAR, [.2,.2,0.1,1.0])
		#glMaterial(GL_FRONT, GL_SHININESS, 32)
		#gluCylinder(self.gq,5,5,cl,16,16)
		#glEndList()
		
		#self.xshapedl = glGenLists(1)
		#glNewList(self.xshapedl,GL_COMPILE)
		#glPushMatrix()
		#glTranslate(0,0,-cl/2)
		#glCallList(self.cylinderdl)
		#glPopMatrix()
		#glPushMatrix()
		#glRotate(90,0,1,0)
		#glTranslate(0,0,-cl/2)
		#glCallList(self.cylinderdl)
		#glPopMatrix()
		#glPushMatrix()
		#glRotate(90,1,0,0)
		#glTranslate(0,0,-cl/2)
		#glCallList(self.cylinderdl)
		#glPopMatrix()
		#glEndList()
		
		t3d = Transform3D()
		rot = {}
		rot["az"] = 0.01
		rot["alt"] = 0.01
		rot["phi"] = 0.01
		
		t3d.set_rotation( EULER_EMAN, rot )
		self.t3d_stack = []
		self.t3d_stack.append(t3d)
		
		# For the time being
		glEnable(GL_CULL_FACE);
		
		glPolygonMode(GL_FRONT,GL_FILL);
	
		glEnableClientState(GL_VERTEX_ARRAY)
		glEnableClientState(GL_NORMAL_ARRAY)
		#glEnableClientState(GL_COLOR_ARRAY)
		
	def paintGL(self):
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		glLightfv(GL_LIGHT0, GL_POSITION, [0.5,0.7,11.,0.])
		# because the viewing volume extends from -0.001 to -10000 in the z direction,
		# translate the scene back into it so it is visible...
		glTranslated(self.cam_x, self.cam_y, self.cam_z)
		#print "Scene is positioned at %f %f %f" %(self.cam_x, self.cam_y, self.cam_z)
		# get the current rotation from the rotation stack and apply
		rot = self.t3d_stack[len(self.t3d_stack)-1].get_rotation()
				
		glRotate(float(rot["phi"]),0,0,1)
		glRotate(float(rot["alt"]),1,0,0)
		glRotate(float(rot["az"]),0,0,1)
		
		# here is where zoom is applied
		glScalef(self.scale,self.scale,self.scale)

		if ( self.isodl == 0 ):
			self.getIsoDL()

		#if self.griddl == 0:
			#self.createGridDL()
			
		glMaterial(GL_FRONT, GL_AMBIENT, self.colors[self.isocolor]["ambient"])
		glMaterial(GL_FRONT, GL_DIFFUSE, self.colors[self.isocolor]["diffuse"])
		glMaterial(GL_FRONT, GL_SPECULAR, self.colors[self.isocolor]["specular"])
		glMaterial(GL_FRONT, GL_SHININESS, self.colors[self.isocolor]["shininess"])
		glColor(self.colors[self.isocolor]["ambient"])
		glPushMatrix()
		glTranslate(-self.data.get_xsize()/2.0,-self.data.get_ysize()/2.0,-self.data.get_zsize()/2.0)
		glCallList(self.isodl)
		glPopMatrix()
	
		if self.cube:
			glPushMatrix()
			self.draw_volume_bounds()
			glPopMatrix()
		
	def draw_volume_bounds(self):
		
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
		
	def getIsoDL(self):
		# create the isosurface display list
		self.isorender.set_surface_value(self.isothr)
		self.isorender.set_sampling(self.smpval)
		
		#time1 = clock()
		self.isodl = self.isorender.get_isosurface_dl()
		#time2 = clock()
		#dt1 = time2 - time1
		#print "It took %f to render the isosurface" %dt1
		
	def createGridDL(self):
		# create the isosurface display list
		
		if ( self.griddl != 0 ): glDeleteLists(self.griddl,1)
		
		width = self.data.get_xsize()
		height = self.data.get_ysize()
		depth = self.data.get_zsize()
		
		self.griddl = glGenLists(1)
		
		glLineWidth(0.5)
		glNewList(self.griddl,GL_COMPILE)
		
		glColor(1,1,1,1.0)
		glMaterial(GL_FRONT, GL_AMBIENT, [1, 1, 1,1.0])
		glMaterial(GL_FRONT, GL_DIFFUSE, [1, 1, 1,1.0])
		glMaterial(GL_FRONT, GL_SPECULAR, [0.774597, 0.774597, 0.774597,1.0])
		glMaterial(GL_FRONT, GL_SHININESS, 76.8)
		
		glPushMatrix()
		glTranslate(-width/2.0,-height/2.0,-depth/2.0)
		glBegin(GL_LINES)
		
		glNormal(0,1,0)
		for i in range(0,height):
			for j in range(0,width):
				glVertex(j,i,0)
				glVertex(j,i,depth)
			
		glNormal(0,0,-1)
		for i in range(0,depth):
			for j in range(0,width):
				glVertex(j,0,i)
				glVertex(j,height,i)
		glNormal(0,1,0)
		for i in range(0,depth):
			for j in range(0,height):
				glVertex(0,j,i)
				glVertex(width,j,i)
		glEnd()
		
		
		glPopMatrix()
		glEndList()
		
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
		
		self.updateInspectorTranslateScale()
		
	def setupShapes(self):
		# make our own cirle rather than use gluDisk or somesuch
		pass
	
	def showInspector(self,force=0):
		if not force and self.inspector==None : return
		
		if not self.inspector : self.inspector=EMImageInspector3D(self)
		self.inspector.show()
	
	def updateInspector(self,t3d):
		if not self.inspector or self.inspector ==None:
			self.inspector=EMImageInspector3D(self)
		self.inspector.updateRotations(t3d)
	
	def closeEvent(self,event) :
		if self.inspector: self.inspector.close()
		
	def mousePressEvent(self, event):
#		lc=self.scrtoimg((event.x(),event.y()))
		if event.button()==Qt.MidButton:
			if not self.inspector or self.inspector ==None:
				self.inspector=EMImageInspector3D(self)
			self.inspector.updateRotations(self.t3d_stack[len(self.t3d_stack)-1])
			self.updateInspectorTranslateScale()
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
				
				self.emit(QtCore.SIGNAL("mousedown"), event)
				return
		elif event.button()==Qt.RightButton:
			if self.mmode==0:
				self.mpressx = event.x()
				self.mpressy = event.y()
				self.emit(QtCore.SIGNAL("mousedown"), event)
				return
		
	def mouseMoveEvent(self, event):
#		lc=self.scrtoimg((event.x(),event.y()))
# 		if self.rmousedrag and event.buttons()&Qt.RightButton:
# 			self.origin=(self.origin[0]+self.rmousedrag[0]-event.x(),self.origin[1]-self.rmousedrag[1]+event.y())
# 			self.rmousedrag=(event.x(),event.y())
# 			self.update()
		if event.buttons()&Qt.LeftButton:
			if self.mmode==0:
				if event.modifiers() == Qt.ControlModifier:
					self.motionTranslate(event.x()-self.mpressx, self.mpressy - event.y())
				else:
					self.motionRotate(self.mpressx - event.x(), self.mpressy - event.y())
				self.updateGL()
				self.mpressx = event.x()
				self.mpressy = event.y()
				self.emit(QtCore.SIGNAL("mousedrag"), event)
				return
		if event.buttons()&Qt.RightButton:
			if self.mmode==0:
				if event.modifiers() == Qt.ControlModifier:
					self.scale_event(event.y()-self.mpressy)	
				else:
					self.motionTranslate(event.x()-self.mpressx, self.mpressy - event.y())
					
				self.mpressx = event.x()
				self.mpressy = event.y()
				self.updateGL()
				
				self.emit(QtCore.SIGNAL("mousedrag"), event)
				return
	
	def mouseReleaseEvent(self, event):
#		lc=self.scrtoimg((event.x(),event.y()))
# 		if event.button()==Qt.RightButton:
# 			self.rmousedrag=None
		if event.button()==Qt.LeftButton:
			if self.mmode==0:
				return
		elif event.button()==Qt.RightButton:
			if self.mmode==0:
				self.emit(QtCore.SIGNAL("mouseup"), event)
				return
			
	def wheelEvent(self, event):
		self.scale_event(event.delta())
		self.updateInspectorTranslateScale()

	def updateInspectorTranslateScale(self):
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


	def get_render_dims_at_depth(self):
		# This function returns the width and height of the renderable 
		# area at the origin of the data volume
		height = -2*tan(self.fov/2.0*pi/180.0)*(self.cam_z)
		width = self.aspect*height
		
		return [width,height]

	def getTranslateScale(self):
	
		[rx,ry] = self.get_render_dims_at_depth()
		
		#print "render area is %f %f " %(xx,yy)
		xscale = rx/float(self.width())
		yscale = ry/float(self.height())
		
		return [xscale,yscale]

	def motionTranslate(self,x,y):
		[xscale,yscale] = self.getTranslateScale()
		
		self.cam_x += x*xscale
		self.cam_y += y*yscale
		
		self.inspector.setXYTrans(self.cam_x, self.cam_y)
		
	def setCamZ(self,z):
		self.cam_z = self.default_z + z
		self.updateGL()
		
	def setCamY(self,y):
		self.cam_y = y
		self.updateGL()
		
	def setCamX(self,x):
		self.cam_x = x
		self.updateGL()
		
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
		self.updateGL()
	
	def loadRotation(self,t3d):
		self.t3d_stack.append(t3d)
		self.updateGL()
	
	def setThr(self,val):
		if (self.isothr != val):
			self.isothr = val
			self.getIsoDL()
			self.updateGL()
	
	def setSample(self,val):
		if ( self.smpval != int(val)):
			# the minus two is here because the marching cubes thinks -1 is the high level of detail, 0 is the next best and  so forth
			# However the user wants the highest level of detail to be 1, and the next best to be 2 and then 3 etc
			self.smpval = int(val)-2
			self.getIsoDL()
			self.updateGL()
	
	def setColor(self,val):
		#print val
		self.isocolor = str(val)
		self.updateGL()
		
	def toggleCube(self):
		self.cube = not self.cube
		self.updateGL()
	
	def toggleWire(self,val):
		self.wire = not self.wire
		
		if ( self.wire ):
			glPolygonMode(GL_FRONT,GL_LINE);
		else:
			glPolygonMode(GL_FRONT,GL_FILL);
		self.updateGL()
		
	def toggleLight(self,val):
		enabled = glIsEnabled(GL_LIGHTING)
		
		if enabled:
			glDisable(GL_LIGHTING)
		else:
			glEnable(GL_LIGHTING)
			
		self.updateGL()

class EMImageInspector3D(QtGui.QWidget):
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
		
		self.wiretog = QtGui.QPushButton("Wire")
		self.wiretog.setCheckable(1)
		self.vbl2.addWidget(self.wiretog)
		
		self.lighttog = QtGui.QPushButton("Light")
		self.lighttog.setCheckable(1)
		self.vbl2.addWidget(self.lighttog)
		
		self.cubetog = QtGui.QPushButton("Cube")
		self.cubetog.setCheckable(1)
		self.vbl2.addWidget(self.cubetog)
		
		self.scale = ValSlider(self,(0.01,30.0),"Mag:")
		self.scale.setObjectName("scale")
		self.scale.setValue(1.0)
		self.vbl.addWidget(self.scale)
		
		self.lowlim=0
		self.highlim=1.0
		self.busy=0

		self.thr = ValSlider(self,(0.0,2.0),"Thr:")
		self.thr.setObjectName("thr")
		self.thr.setValue(0.5)
		self.vbl.addWidget(self.thr)
		
		self.hbl_smp = QtGui.QHBoxLayout()
		self.hbl_smp.setMargin(0)
		self.hbl_smp.setSpacing(6)
		self.hbl_smp.setObjectName("Sample")
		self.vbl.addLayout(self.hbl_smp)
		
		self.smp_label = QtGui.QLabel()
		self.smp_label.setText('Sample Level')
		self.hbl_smp.addWidget(self.smp_label)
		
		self.smp = QtGui.QSpinBox(self)
		self.smp.setValue(1)
		self.hbl_smp.addWidget(self.smp)
		
		self.hbl_color = QtGui.QHBoxLayout()
		self.hbl_color.setMargin(0)
		self.hbl_color.setSpacing(6)
		self.hbl_color.setObjectName("Material")
		self.vbl.addLayout(self.hbl_color)
		
		self.color_label = QtGui.QLabel()
		self.color_label.setText('Material')
		self.hbl_color.addWidget(self.color_label)
		
		self.cbb = QtGui.QComboBox(self)
		self.hbl_color.addWidget(self.cbb)

		self.hbl_trans = QtGui.QHBoxLayout()
		self.hbl_trans.setMargin(0)
		self.hbl_trans.setSpacing(6)
		self.hbl_trans.setObjectName("Trans")
		self.vbl.addLayout(self.hbl_trans)
		
		self.x_label = QtGui.QLabel()
		self.x_label.setText('x')
		self.hbl_trans.addWidget(self.x_label)
		
		self.x_trans = QtGui.QDoubleSpinBox(self)
		self.x_trans.setMinimum(-10000)
		self.x_trans.setMaximum(10000)
		self.x_trans.setValue(0.0)
		self.hbl_trans.addWidget(self.x_trans)
		
		self.y_label = QtGui.QLabel()
		self.y_label.setText('y')
		self.hbl_trans.addWidget(self.y_label)
		
		self.y_trans = QtGui.QDoubleSpinBox(self)
		self.y_trans.setMinimum(-10000)
		self.y_trans.setMaximum(10000)
		self.y_trans.setValue(0.0)
		self.hbl_trans.addWidget(self.y_trans)
		
		
		self.z_label = QtGui.QLabel()
		self.z_label.setText('z')
		self.hbl_trans.addWidget(self.z_label)
		
		self.z_trans = QtGui.QDoubleSpinBox(self)
		self.z_trans.setMinimum(-10000)
		self.z_trans.setMaximum(10000)
		self.z_trans.setValue(0.0)
		self.hbl_trans.addWidget(self.z_trans)
		
		self.hbl_src = QtGui.QHBoxLayout()
		self.hbl_src.setMargin(0)
		self.hbl_src.setSpacing(6)
		self.hbl_src.setObjectName("hbl")
		self.vbl.addLayout(self.hbl_src)
		
		self.label_src = QtGui.QLabel()
		self.label_src.setText('Rotation Convention')
		self.hbl_src.addWidget(self.label_src)
		
		self.src = QtGui.QComboBox(self)
		self.load_src_options(self.src)
		self.hbl_src.addWidget(self.src)
		
		# set default value -1 ensures that the val slider is updated the first time it is created
		self.az = ValSlider(self,(-360.0,360.0),"az",-1)
		self.az.setObjectName("az")
		self.vbl.addWidget(self.az)
		
		self.alt = ValSlider(self,(-180.0,180.0),"alt",-1)
		self.alt.setObjectName("alt")
		self.vbl.addWidget(self.alt)
		
		self.phi = ValSlider(self,(-360.0,360.0),"phi",-1)
		self.phi.setObjectName("phi")
		self.vbl.addWidget(self.phi)
		
		self.n3_showing = False
		
		self.current_src = EULER_EMAN
		
		QtCore.QObject.connect(self.scale, QtCore.SIGNAL("valueChanged"), target.setScale)
		QtCore.QObject.connect(self.az, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
		QtCore.QObject.connect(self.alt, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
		QtCore.QObject.connect(self.phi, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
		QtCore.QObject.connect(self.thr, QtCore.SIGNAL("valueChanged"), target.setThr)
		QtCore.QObject.connect(self.cbb, QtCore.SIGNAL("currentIndexChanged(QString)"), target.setColor)
		QtCore.QObject.connect(self.src, QtCore.SIGNAL("currentIndexChanged(QString)"), self.set_src)
		QtCore.QObject.connect(self.smp, QtCore.SIGNAL("valueChanged(int)"), target.setSample)
		QtCore.QObject.connect(self.x_trans, QtCore.SIGNAL("valueChanged(double)"), target.setCamX)
		QtCore.QObject.connect(self.y_trans, QtCore.SIGNAL("valueChanged(double)"), target.setCamY)
		QtCore.QObject.connect(self.z_trans, QtCore.SIGNAL("valueChanged(double)"), target.setCamZ)
		QtCore.QObject.connect(self.wiretog, QtCore.SIGNAL("toggled(bool)"), target.toggleWire)
		QtCore.QObject.connect(self.lighttog, QtCore.SIGNAL("toggled(bool)"), target.toggleLight)
		QtCore.QObject.connect(self.cubetog, QtCore.SIGNAL("toggled(bool)"), target.toggleCube)
	
	def setSamplingRange(self,range):
		self.smp.setMinimum(1)
		self.smp.setMaximum(1+range-1)

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
		a = 0
		for i in colors:
			self.cbb.addItem(i)
			if ( i == current_color):
				self.cbb.setCurrentIndex(a)
			a += 1

	def setThrs(self,low,high,val):
		self.thr.setRange(low,high)
		self.thr.setValue(val, True)
	
	def setSamp(self,low,high,val):
		self.smp.setRange(int(low),int(high))
		self.smp.setValue(val, True)
		
	def setHist(self,hist,minden,maxden):
		self.hist.setData(hist,minden,maxden)

	def setScale(self,newscale):
		self.scale.setValue(newscale)
		
# This is just for testing, of course
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = EMImage3D()
 	if len(sys.argv)==1 : 
		e = EMData()
		e.set_size(40,35,30)
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
