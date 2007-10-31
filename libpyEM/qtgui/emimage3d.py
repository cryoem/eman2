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
		self.generate_iso = True
		self.inspector=None
		
		if image :
			print "We have an image"
			self.setData(image)
			self.show()
		
	
		ruby = {}
		ruby["ambient"] = [0.1745, 0.01175, 0.01175,1.0]
		ruby["diffuse"] = [0.61424, 0.04136, 0.04136,1.0]
		ruby["specular"] = [0.727811, 0.626959, 0.626959,1.0]
		ruby["shininess"] = 76.8
		
		emerald = {}
		emerald["ambient"] = [0.0215, 0.1745, 0.0215,1.0]
		emerald["diffuse"] = [0.07568, 0.61424,  0.07568,1.0]
		emerald["specular"] = [0.633, 0.727811, 0.633,1.0]
		emerald["shininess"] = 76.8
		
		pearl = {}
		pearl["ambient"] = [0.25, 0.20725, 0.20725,1.0]
		pearl["diffuse"] = [1.0, 0.829, 0.829,1.0]
		pearl["specular"] = [0.296648, 0.296648, 0.296648,1.0]
		pearl["shininess"] = 11.264
		
		silver = {}
		silver["ambient"] = [0.25, 0.25, 0.25,1.0]
		silver["diffuse"] = [0.4, 0.4, 0.4,1.0]
		silver["specular"] = [0.774597, 0.774597, 0.774597,1.0]
		silver["shininess"] = 76.8
		
		gold = {}
		gold["ambient"] = [0.24725, 0.2245, 0.0645,1.0]
		gold["diffuse"] = [0.34615, 0.3143, 0.0903,1.0]
		gold["specular"] = [0.797357, 0.723991, 0.208006,1.0]
		gold["shininess"] = 83.2
		
		copper = {}
		copper["ambient"] = [0.2295, 0.08825, 0.0275,1.0]
		copper["diffuse"] = [0.5508, 0.2118, 0.066,1.0]
		copper["specular"] = [0.580594, 0.223257, 0.0695701,1.0]
		copper["shininess"] = 51.2
		
		obsidian = {}
		obsidian["ambient"] = [0.05375,  0.05,     0.06625 ,1.0]
		obsidian["diffuse"] = [0.18275,  0.17,     0.22525,1.0]
		obsidian["specular"] = [0.332741, 0.328634, 0.346435]
		obsidian["shininess"] = 38.4
		
		turquoise = {}
		turquoise["ambient"] = [0.1, 0.18725, 0.1745 ,1.0]
		turquoise["diffuse"] = [0.396, 0.74151, 0.69102,1.0]
		turquoise["specular"] = [0.297254, 0.30829, 0.306678]
		turquoise["shininess"] = 12.8
		
		self.colors = {}
		self.colors["ruby"] = ruby
		self.colors["emerald"] = emerald
		self.colors["pearl"] = pearl
		self.colors["silver"] = silver
		self.colors["gold"] = gold
		self.colors["copper"] = copper
		self.colors["obsidian"] = obsidian
		self.colors["turquoise"] = turquoise
		
		
		self.colormap = {}
		self.colormap[0] = "ruby"
		self.colormap[1] = "emerald"
		self.colormap[2] = "pearl"
		self.colormap[3] = "silver"
		self.colormap[4] = "gold"
		self.colormap[5] = "copper"
		self.colormap[6] = "obsidian"
		self.colormap[7] = "turquoise"
		
		self.isocolor = "ruby"
		
	def timeout(self):
		self.updateGL()
		
	
	def setData(self,data):
		"""Pass in a 3D EMData object"""
#		if not self.data and data: self.resize(data.get_xsize(),data.get_ysize())
		
		self.data=data
		if data==None or (isinstance(data,EMData) and data.get_zsize()<=1) :
			self.updateGL()
			return
		
		m0=data.get_attr("minimum")
		m1=data.get_attr("maximum")
		mean=data.get_attr("mean")
		sigma = data.get_attr("sigma")
		
		
		if not self.inspector or self.inspector ==None:
			self.inspector=EMImageInspector3D(self)
		
		
		self.start_z = -1.25*data.get_zsize()
		
		#a=self.data.render_amp8()
		#hist=numpy.fromstring(a[-1024:],'i')
		
		#self.inspector.setHist(hist,self.minden,self.maxden) 
		
		self.inspector.setThrs(m0,m1,mean+3.0*sigma)
		self.isothr = mean+3.0*sigma
		
		self.isorender=MarchingCubes(data,1)
		
		#self.inspector.setSamp(self.isorender.get_root_level(), self.isorender.get_leaf_level(), self.smpval)
		
		self.inspector.setColors(self.colors)
		
		
		#error.ErrorChecker.registerChecker(self.myAlternateFunction )
		
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

		GL.glClearColor(0,0,0,0)
		
		if not self.gq:
			self.gq=gluNewQuadric()
			gluQuadricDrawStyle(self.gq,GLU_FILL)
			gluQuadricNormals(self.gq,GLU_SMOOTH)
			gluQuadricOrientation(self.gq,GLU_OUTSIDE)
			gluQuadricTexture(self.gq,GL_FALSE)
		
		# Precompile a displaylist for the display volume border
		self.volcubedl=glGenLists(1)
		glNewList(self.volcubedl,GL_COMPILE)
		glPushMatrix()
		glColor(.7,.7,1.0)
#		glRotate(90.,1.,0.,0.)
		glTranslate(-self.aspect-.01,1.01,-4.0)
		gluCylinder(self.gq,.01,.01,15.0,12,2)
		glTranslate(self.aspect*2.0+.02,0.0,0.0)
		gluCylinder(self.gq,.01,.01,15.0,12,2)
		glTranslate(0.0,-2.02,0.0)
		gluCylinder(self.gq,.01,.01,15.0,12,2)
		glTranslate(-self.aspect*2.0-.02,0.0,0.0)
		gluCylinder(self.gq,.01,.01,15.0,12,2)		
		glPopMatrix()
		glEndList()
		
		cl = 50
		self.cylinderdl = glGenLists(1)
		glNewList(self.cylinderdl,GL_COMPILE)
		glColor(.7,.7,1.0)
		glMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, [.2,.1,0.4,1.0])
		glMaterial(GL_FRONT, GL_SPECULAR, [.2,.2,0.1,1.0])
		glMaterial(GL_FRONT, GL_SHININESS, 32)
		gluCylinder(self.gq,5,5,cl,16,16)
		glEndList()
		
		self.xshapedl = glGenLists(1)
		glNewList(self.xshapedl,GL_COMPILE)
		glPushMatrix()
		glTranslate(0,0,-cl/2)
		glCallList(self.cylinderdl)
		glPopMatrix()
		glPushMatrix()
		glRotate(90,0,1,0)
		glTranslate(0,0,-cl/2)
		glCallList(self.cylinderdl)
		glPopMatrix()
		glPushMatrix()
		glRotate(90,1,0,0)
		glTranslate(0,0,-cl/2)
		glCallList(self.cylinderdl)
		glPopMatrix()
		glEndList()
		
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
		
		# because the viewing volume extends from -0.001 to -10000 in the z direction,
		# translate the scene back into it so it is visible...
		glTranslated(0, 0, self.start_z)

		# get the current rotation from the rotation stack and apply
		rot = self.t3d_stack[len(self.t3d_stack)-1].get_rotation()
		self.updateInspector(rot)
				
		glRotate(float(rot["phi"]),0,0,1)
		glRotate(float(rot["alt"]),1,0,0)
		glRotate(float(rot["az"]),0,0,1)
		
		# here is where zoom is effectively applied
		glScalef(self.scale,self.scale,self.scale)
		
		if self.isodl == 0:
			#self.createIsoDL()
			self.getIsoDL()

		#if self.griddl == 0:
			#self.createGridDL()
		glMaterial(GL_FRONT, GL_AMBIENT, self.colors[self.isocolor]["ambient"])
		glMaterial(GL_FRONT, GL_DIFFUSE, self.colors[self.isocolor]["diffuse"])
		glMaterial(GL_FRONT, GL_SPECULAR, self.colors[self.isocolor]["specular"])
		glMaterial(GL_FRONT, GL_SHININESS, self.colors[self.isocolor]["shininess"])
		#self.draw_iso()
		#glEnable(GL_COLOR_MATERIAL)
		#glColorMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE)
		glPushMatrix()
		#glScalef(self.data.get_xsize(), self.data.get_ysize(), self.data.get_zsize())
		glTranslate(-self.data.get_xsize()/2.0,-self.data.get_ysize()/2.0,-self.data.get_zsize()/2.0)
		glCallList(self.isodl)
		glPopMatrix()
		#glCallList(self.griddl)
	
	def getIsoDL(self):
		# create the isosurface display list
		
		print "!!!! %f %f" %(self.isothr,self.smpval)
		#if ( self.isodl != 0 ): glDeleteLists(self.isodl,1)
		#time1 = clock()
		time1 = clock()
		self.isorender.set_surface_value(self.isothr)
		self.isorender.set_sample_density(self.smpval)
		time2 = clock()
		dt1 = time2 - time1
		#print "It took %f to traverse the cubes" %dt1
		
		time1 = clock()
		self.isodl = self.isorender.get_isosurface_dl(True)
		time2 = clock()
		dt1 = time2 - time1
		print "It took %f to render the isosurface" %dt1
	
	def createIsoDL(self):
		# create the isosurface display list
		
		print "!!!!"
		#if ( self.isodl != 0 ): glDeleteLists(self.isodl,1)
		#time1 = clock()
		time1 = clock()
		self.isorender.set_surface_value(self.isothr)
		time2 = clock()
		dt1 = time2 - time1
		print "It took %f to traverse the cubes" %dt1
		#time2 = clock()
		#dt1 = time2 - time1
		#print "It took %f to traverse the tree" %dt1
		##self.isorender.set_sample_density(self.smpval)
		#time1 = clock()
		#self.isodl = self.isorender.get_isosurface_dl(True)
		#time2 = clock()
		#dt1 = time2 - time1
		#print "It took %f to get the isosurface" %dt1
		
		
		time1 = clock()
		a=self.isorender.get_isosurface(True)
		time2 = clock()
		dt1 = time2 - time1
		print "It took %f to get the data" %dt1
		
		self.f=a["faces"]
		self.n=a["normals"]
		self.p=a["points"]
	
		self.colors = []
		
		center = [0.5,0.5,0.5]
		gr = 0.5
		
		#for i in range(0,len(self.p),3):
			#print "%f %f %f" %(self.p[i],self.p[i+1],self.p[i+2])
		
		for i in range(0,len(self.p),3):
			v = [self.p[i],self.p[i+1],self.p[i+2]]
			c = center
			distance = sqrt((c[0]-v[0])*(c[0]-v[0]) + (c[1]-v[1])*(c[1]-v[1]) + (c[2]-v[2])*(c[2]-v[2]))
			
			if ( distance < 0.5 ):
				self.colors.append((1.0-distance*2.0))
				self.colors.append(2.0*distance)
				self.colors.append(0)
			if ( distance >= 0.5 ):
				d = distance - 0.5
				self.colors.append(0)
				self.colors.append((1.0-d*2.0))
				self.colors.append(2.0*d)
		
		self.f=[i/3 for i in self.f]
		
		glNormalPointer(GL_FLOAT,0,self.n)
		glVertexPointer(3,GL_FLOAT,0,self.p)
		glColorPointer(3,GL_FLOAT,0,self.colors)
		
		print "Num triangles is %d" %(len(self.f)/3)
		print "Num vertices is %d" %(len(self.p)/3)
		
		self.isodl = glGenLists(1)
		glNewList(self.isodl,GL_COMPILE)
		
		time1 = clock()
		#glBegin(GL_TRIANGLES)
		#for i in self.f:
			#glArrayElement(i)
		#glEnd()
		
		#glDrawArrays(GL_TRIANGLES,0,len(self.f)/3)
		glDrawElements(GL_TRIANGLES,len(self.f),GL_UNSIGNED_INT,self.f)
		
		time2 = clock()
		dt1 = time2 - time1
		print "It took %f to render the isosurface" %dt1
	
		glEndList()
	
	def sq_dist(c,v):
		return (c[0]-v[0])*(c[0]-v[0]) + (c[1]-v[1])*(c[1]-v[1]) + (c[2]-v[2])*(c[2]-v[2])
	
	def createGridDL(self):
		# create the isosurface display list
		
		if ( self.griddl != 0 ): glDeleteLists(self.griddl,1)
		
		width = self.data.get_xsize()
		height = self.data.get_ysize()
		depth = self.data.get_zsize()
		
		self.griddl = glGenLists(1)
		
		glLineWidth(0.5)
		glNewList(self.griddl,GL_COMPILE)
		
		glColor(.2,.1,0.4,1.0)
		glMaterial(GL_FRONT, GL_AMBIENT, [0.25, 0.25, 0.25,1.0])
		glMaterial(GL_FRONT, GL_DIFFUSE, [0.4, 0.4, 0.4,1.0])
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
		aspect = float(width)/float(height)
		
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		# using gluPerspective for simplicity
		gluPerspective(50,aspect,0.001,1000)
		
		# switch back to model view mode
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
	def setupShapes(self):
		# make our own cirle rather than use gluDisk or somesuch
		pass
	
	def showInspector(self,force=0):
		if not force and self.inspector==None : return
		
		if not self.inspector : self.inspector=EMImageInspector3D(self)
		self.inspector.show()
	
	def updateInspector(self,rot):
		if not self.inspector or self.inspector ==None:
			self.inspector=EMImageInspector3D(self)
		self.inspector.newAz(float(rot["az"])+180)
		self.inspector.newAlt(rot["alt"])
		self.inspector.newPhi(float(rot["phi"])+180)
	
	def closeEvent(self,event) :
		if self.inspector: self.inspector.close()
			
	def mousePressEvent(self, event):
#		lc=self.scrtoimg((event.x(),event.y()))
		if event.button()==Qt.MidButton:
			self.showInspector(1)
		elif event.button()==Qt.LeftButton:
			if self.mmode==0:
				self.emit(QtCore.SIGNAL("mousedown"), event)
				return
		elif event.button()==Qt.RightButton:
			if self.mmode==0:
				self.rpressx = event.x()
				self.rpressy = event.y()
				
				tmp = self.t3d_stack.pop()
				t3d = Transform3D(tmp)
				self.t3d_stack.append(tmp)
				self.t3d_stack.append(t3d)
				
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
				self.emit(QtCore.SIGNAL("mousedrag"), event)
				return
		if event.buttons()&Qt.RightButton:
			if self.mmode==0:
				self.motionRotate(self.rpressx - event.x(), self.rpressy - event.y())
				self.updateGL()
				self.rpressx = event.x()
				self.rpressy = event.y()
				self.emit(QtCore.SIGNAL("mousedrag"), event)
				return
	
	def mouseReleaseEvent(self, event):
#		lc=self.scrtoimg((event.x(),event.y()))
# 		if event.button()==Qt.RightButton:
# 			self.rmousedrag=None
		if event.button()==Qt.LeftButton:
			if self.mmode==0:
				self.emit(QtCore.SIGNAL("mouseup"), event)
				return
		elif event.button()==Qt.RightButton:
			if self.mmode==0:
				self.rreleasex = event.x()
				self.rreleasey = event.y()
				self.emit(QtCore.SIGNAL("mousedown"), event)
				return

	def motionRotate(self,x,y):
		if ( x == 0 and y == 0): return
		
		theta = atan2(-y,x)

		rotaxis_x = sin(theta)
		rotaxis_y = cos(theta)
		rotaxis_z = 0
		
		length = sqrt(x*x + y*y)
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
		self.updateGL()
		
	def setAz(self,val):
		rot = self.t3d_stack[len(self.t3d_stack)-1].get_rotation()
		rot["az"] = val - 180.0
		self.t3d_stack.append(Transform3D( EULER_EMAN, rot))
		self.updateGL()
		
	def setAlt(self,val):
		rot = self.t3d_stack[len(self.t3d_stack)-1].get_rotation()
		rot["alt"] = val
		self.t3d_stack.append(Transform3D( EULER_EMAN, rot))
		self.updateGL()
		
	def setPhi(self,val):
		rot = self.t3d_stack[len(self.t3d_stack)-1].get_rotation()
		rot["phi"] = val - 180.0
		self.t3d_stack.append(Transform3D( EULER_EMAN, rot))
		self.updateGL()
	
	def setThr(self,val):
		if (self.isothr != val):
			self.isothr = val
			#self.createIsoDL()
			self.getIsoDL()
			self.updateGL()
	
	def setSample(self,val):
		if ( self.smpval != int(val)):
			self.smpval = int(val)
			self.getIsoDL()
			self.updateGL()
	
	
	def setColor(self,val):
		self.isocolor = self.colormap[val]
		print val
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
		
		self.invtog = QtGui.QPushButton("Invert")
		self.invtog.setCheckable(1)
		self.vbl2.addWidget(self.invtog)
		
		self.ffttog = QtGui.QPushButton("FFT")
		self.ffttog.setCheckable(1)
		self.vbl2.addWidget(self.ffttog)

		self.hbl2 = QtGui.QHBoxLayout()
		self.hbl2.setMargin(0)
		self.hbl2.setSpacing(6)
		self.hbl2.setObjectName("hboxlayout")
		self.vbl.addLayout(self.hbl2)
		
		self.mapp = QtGui.QPushButton("App")
		self.mapp.setCheckable(1)
		self.hbl2.addWidget(self.mapp)

		self.mmeas = QtGui.QPushButton("Meas")
		self.mmeas.setCheckable(1)
		self.hbl2.addWidget(self.mmeas)

		self.mmode=QtGui.QButtonGroup()
		self.mmode.setExclusive(1)
		self.mmode.addButton(self.mapp,0)
		self.mmode.addButton(self.mmeas,1)
		
		self.scale = ValSlider(self,(0.01,30.0),"Mag:")
		self.scale.setObjectName("scale")
		self.scale.setValue(1.0)
		self.vbl.addWidget(self.scale)
		
		self.lowlim=0
		self.highlim=1.0
		self.busy=0

		self.thr = ValSlider(self,(0.0,2.0),"Thr:")
		self.thr.setObjectName("az")
		self.thr.setValue(0.5)
		self.vbl.addWidget(self.thr)
		
		self.smp = QtGui.QSpinBox(self)
		self.smp.setMinimum(-1)
		self.smp.setMaximum(4)
		self.smp.setValue(-1)
		self.vbl.addWidget(self.smp)
			
		self.az = ValSlider(self,(0.0,360.0),"Az:")
		self.az.setObjectName("az")
		self.az.setValue(0.0)
		self.vbl.addWidget(self.az)
		
		self.alt = ValSlider(self,(0.01,180.0),"Alt:")
		self.alt.setObjectName("alt")
		self.alt.setValue(0.0)
		self.vbl.addWidget(self.alt)
		
		self.phi = ValSlider(self,(0.0,360.0),"Phi:")
		self.phi.setObjectName("phi")
		self.phi.setValue(0.0)
		self.vbl.addWidget(self.phi)

		self.cbb = QtGui.QComboBox(self)
		
		self.vbl.addWidget(self.cbb)

		QtCore.QObject.connect(self.scale, QtCore.SIGNAL("valueChanged"), target.setScale)
		QtCore.QObject.connect(self.az, QtCore.SIGNAL("valueChanged"), target.setAz)
		QtCore.QObject.connect(self.alt, QtCore.SIGNAL("valueChanged"), target.setAlt)
		QtCore.QObject.connect(self.phi, QtCore.SIGNAL("valueChanged"), target.setPhi)
		QtCore.QObject.connect(self.thr, QtCore.SIGNAL("valueChanged"), target.setThr)
		QtCore.QObject.connect(self.cbb, QtCore.SIGNAL("currentIndexChanged(int)"), target.setColor)
		QtCore.QObject.connect(self.smp, QtCore.SIGNAL("valueChanged(int)"), target.setSample)
		#QtCore.QObject.connect(self.cbb, QtCore.SIGNAL("textChanged(int)"), target.setColor)
	
	def setColors(self,colors):
		for i in colors:
			self.cbb.addItem(i)
	
	def newAz(self,val):
		if self.busy : return
		self.busy=1
		self.az.setValue(val, True)
		self.busy=0

	def newAlt(self,val):
		if self.busy : return
		self.busy=1
		self.alt.setValue(val, True)
		self.busy=0
	
	def newPhi(self,val):
		if self.busy : return
		self.busy=1
		self.phi.setValue(val, True)
		self.busy=0

	def setThrs(self,low,high,val):
		self.thr.setRange(low,high)
		self.thr.setValue(val, True)
	
	def setSamp(self,low,high,val):
		self.smp.setRange(int(low),int(high))
		self.smp.setValue(val, True)
		
	def setHist(self,hist,minden,maxden):
		self.hist.setData(hist,minden,maxden)

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
