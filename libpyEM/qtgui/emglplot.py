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
from emimageutil import ImgHistogram,EMParentWin
from weakref import WeakKeyDictionary
from time import time
from PyQt4.QtCore import QTimer

from time import *

from emglobjects import EMImage3DObject, Camera2

MAG_INCREMENT_FACTOR = 1.1

class EMGLPlot(EMImage3DObject):
	def __init__(self, parent=None):
		EMImage3DObject.__init__(self)
		self.parent = parent
		
		self.init()
		self.initialized = True
		
		self.cam=Camera2(self)
		
		self.brightness = 0
		self.contrast = 10
		self.glcontrast = 1.0
		self.glbrightness = 0.0
		self.rank = 1
		self.inspector=None
		
		self.data = [1.0,0.93,0.64,0.94,0.30,0.2,0.12,-0.2]
		
		self.cylinderdl = 0
		
		self.gq=gluNewQuadric()
		gluQuadricDrawStyle(self.gq,GLU_FILL)
		gluQuadricNormals(self.gq,GLU_SMOOTH)
		gluQuadricOrientation(self.gq,GLU_OUTSIDE)
		gluQuadricTexture(self.gq,GL_FALSE)
		
		self.marker = None
		
		self.moving = False
		

	def cylinderToFrom(self,next,prev,val=0.5):
		dx = next[0] - prev[0]
		dy = next[1] - prev[1]
		dz = next[2] - prev[2]
		
		length = sqrt(dx**2 + dy**2 + dz**2)
		
		alt = acos(dz/length)*180.0/pi
		phi = atan2(dy,dx)*180.0/pi
		
		glPushMatrix()
		glTranslatef(prev[0],prev[1],prev[2] )
		#print "positioned at", prev[0],prev[1],prev[2]
		glRotatef(90+phi,0,0,1)
		glRotatef(alt,1,0,0)
		
		glScalef(val,val,length)
		glCallList(self.cylinderdl)
		glPopMatrix()

	def getType(self):
		return "emglplot"
	
	def updateGL(self):
		self.parent.updateGL()
		
	def unitcircle(self):
		n = 12.0
		for i in range(0,int(n)):
			i = float(i)
			a = sin(i*2.0*pi/n)
			b = cos(i*2.0*pi/n)
			c = sin((i+1)*2.0*pi/n)
			d = cos((i+1)*2.0*pi/n)
			
			p1 = [a,b,0.]
			p2 = [c,d,0.]
			self.cylinderToFrom(p2,p1,0.3)
	def render(self):
		if ( self.cylinderdl == 0 ):
			self.cylinderdl=glGenLists(1)
				
			glNewList(self.cylinderdl,GL_COMPILE)
			glPushMatrix()
			gluCylinder(self.gq,1.0,1.0,1.0,12,2)
			glPopMatrix()
				
			glEndList()
		
		#if (not isinstance(self.data,EMData)): return
		lighting = glIsEnabled(GL_LIGHTING)
		cull = glIsEnabled(GL_CULL_FACE)
		depth = glIsEnabled(GL_DEPTH_TEST)
		polygonmode = glGetIntegerv(GL_POLYGON_MODE)

		glDisable(GL_CULL_FACE)
		glEnable(GL_DEPTH_TEST)
		
		if ( self.wire ):
			glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
		else:
			glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
		
		if self.light:
			glEnable(GL_LIGHTING)
		else:
			glDisable(GL_LIGHTING)

			
		glShadeModel(GL_SMOOTH)

		glStencilFunc(GL_EQUAL,self.rank,0)
		glStencilOp(GL_KEEP,GL_KEEP,GL_REPLACE)
		glMaterial(GL_FRONT, GL_AMBIENT, self.colors[self.currentcolor]["ambient"])
		glMaterial(GL_FRONT, GL_DIFFUSE, self.colors[self.currentcolor]["diffuse"])
		glMaterial(GL_FRONT, GL_SPECULAR, self.colors[self.currentcolor]["specular"])
		glMaterial(GL_FRONT, GL_SHININESS, self.colors[self.currentcolor]["shininess"])
		glColor(self.colors[self.currentcolor]["ambient"])
		
		glEnable(GL_NORMALIZE)
		#HERE
		glPushMatrix()
		
		n = len(self.data)
		for i in range(0,n-1):
			p1 = [i,self.data[i],0]
			p2 = [i+1,self.data[i+1],0]
			self.cylinderToFrom(p2,p1,0.01)

		glPopMatrix()
		
		if self.marker != None:
			glPushMatrix()
			glTranslate(self.marker[0],self.marker[1],0)
			glScale(self.marker[2],self.marker[3],1.0)
			self.unitcircle()
			glPopMatrix()
		
		glStencilFunc(GL_EQUAL,self.rank,self.rank)
		glStencilOp(GL_KEEP,GL_KEEP,GL_KEEP)
		glPushMatrix()
		glLoadIdentity()
		glScalef(10,10,1)
		glTranslate(-0.5,-0.5,-1)
		self.draw_bc_screen()
		glPopMatrix()
		
		glStencilFunc(GL_ALWAYS,1,1)

			
		if ( lighting ): glEnable(GL_LIGHTING)
		else: glDisable(GL_LIGHTING)	
		if ( cull ): glEnable(GL_CULL_FACE)
		else: glDisable(GL_CULL_FACE)
		if ( depth ): glEnable(GL_DEPTH_TEST)
		else : glDisable(GL_DEPTH_TEST)
		
		if ( polygonmode[0] == GL_LINE ): glPolygonMode(GL_FRONT, GL_LINE)
		else: glPolygonMode(GL_FRONT, GL_FILL)
		if ( polygonmode[1] == GL_LINE ): glPolygonMode(GL_BACK, GL_LINE)
		else: glPolygonMode(GL_BACK, GL_FILL)
		
	def init(self):
		self.mmode = 0
		self.wire = False
		self.light = True
	
	def setInit(self):

		self.cam.default_z = -1.25*32
		self.cam.cam_z = -1.25*32
		
		if not self.inspector or self.inspector ==None:
			self.inspector=EMGLPlotInspector(self)
		
		self.loadColors()
		self.inspector.setColors(self.colors,self.currentcolor)
	def loadColors(self):
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
		
		self.currentcolor = "obsidian"
	
	def setColor(self,val):
		#print val
		self.currentcolor = str(val)
		self.parent.updateGL()
	
	def toggleWire(self,val):
		self.wire = not self.wire
		self.parent.updateGL()
		
	def toggleLight(self,val):
		self.light = not self.light
		self.parent.updateGL()
	
	def updateInspector(self,t3d):
		if not self.inspector or self.inspector ==None:
			self.inspector=EMGLPlotInspector(self)
		self.inspector.updateRotations(t3d)
	
	def getInspector(self):
		if not self.inspector : self.inspector=EMGLPlotInspector(self)
		return self.inspector
		
	def mousePressEvent(self, event):
#		lc=self.scrtoimg((event.x(),event.y()))
		if event.button()==Qt.MidButton:
			if not self.inspector or self.inspector ==None:
				return
			self.inspector.updateRotations(self.cam.t3d_stack[len(self.cam.t3d_stack)-1])
			self.resizeEvent()
			self.showInspector(1)
		else:
			if (self.marker != None ):
				self.moving = True
				self.mc = [event.x(),event.y()]
			else :
				self.moving = False
		
		self.updateGL()
		
	def mouseMoveEvent(self, event):
		wr = float(self.width())/self.parent.width()
		hr = float(self.height())/self.parent.height()
		wx = wr*event.x()
		wy = hr*(self.parent.height()-event.y())+min(self.data)
		
		ratio = hr/wr
		pixelnearness = 200
		
		realx = sqrt((pixelnearness*wr*wr)/(ratio**2+1))
		realy = ratio*realx
		
		if self.moving:
			min1 = min(self.data)
			max1 = max(self.data)
			oldwy = hr*(self.parent.height()-self.mc[1])+min(self.data)
			dy = wy-oldwy
			self.mc[1] = event.y()
			self.data[self.marker[4]] += dy
			self.marker[1] =  self.data[self.marker[4]]
			self.updateGL()
			
			min2 = min(self.data)
			max2 = max(self.data)
			if min2 < min1 or max2 > max1:
				self.parent.resizeGL(self.parent.width(),self.parent.height())
				
			return

		else:
		
			for i in range(0,len(self.data)):
				x = i
				y = self.data[i]
				dx = abs(wx-x)
				dy = abs(wy-y)
				
				dis = (dx/wr)**2+(dy/hr)**2
				
				if dis < pixelnearness:
					self.marker = [x,y,realx,realy,i]
					self.updateGL()
					return
		
		self.marker = None
		#self.cam.mouseMoveEvent(event)
		self.updateGL()
	
	def mouseReleaseEvent(self, event):
		self.moving = False
		self.updateGL()
			
	def wheelEvent(self, event):
		self.cam.wheelEvent(event)
		self.updateGL()
	
	def width(self):
		range = self.xrange()
		return range[1]-range[0]
	
	def height(self):
		range = self.yrange()
		return range[1]-range[0]
	
	def xrange(self):
		return [0,len(self.data)-1]
	
	def yrange(self):
		return [min(self.data),max(self.data)]
		
		
class EMGLPlotWidget(QtOpenGL.QGLWidget):
	""" This class is not yet complete.
	A QT widget for rendering 3D EMData objects.
	"""
	allim=WeakKeyDictionary()
	def __init__(self, parent=None):
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True)
		fmt.setDepth(1)
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		

		EMGLPlotWidget.allim[self]=0
		self.plot = EMGLPlot(self)
		self.timer = QTimer()
		QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.timeout)

		self.setMouseTracking(True)

		self.aspect=1.0
		self.fov = 50 # field of view angle used by gluPerspective
	def timeout(self):
		self.updateGL()
		
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
		
		# For the time being
		
	def paintGL(self):
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		# using gluPerspective for simplicity
		#gluPerspective(self.fov,self.aspect,1,5000)
		
		xr = self.plot.xrange()
		yr = self.plot.yrange()
		glOrtho(xr[0],xr[1],yr[0],yr[1],-1,1)
		
		# switch back to model view mode
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		glPushMatrix()
		self.plot.render()
		glPopMatrix()

	def resizeGL(self, width, height):
		if width<=0 or height<=0 : 
			print "bad size"
			return
		# just use the whole window for rendering
		glViewport(0,0,self.width(),self.height())
		
		# maintain the aspect ratio of the window we have
		self.aspect = float(self.width())/float(self.height())
		
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		# using gluPerspective for simplicity
		#gluPerspective(self.fov,self.aspect,1,5000)
		
		xr = self.plot.xrange()
		yr = self.plot.yrange()
		glOrtho(xr[0],xr[1],yr[0],yr[1],-1,1)
		
		# switch back to model view mode
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		self.plot.resizeEvent()
		
		self.updateGL()

	def setData(self,data):
		self.plot.data = data
		self.updateGL()

	def setInit(self):
		self.plot.setInit()

	def showInspector(self,force=0):
		self.plot.showInspector()
	
	def closeEvent(self,event) :
		self.plot.closeEvent(event)
		
	def mousePressEvent(self, event):
		self.plot.mousePressEvent(event)
		self.updateGL()
		
	def mouseMoveEvent(self, event):
		self.plot.mouseMoveEvent(event)
	
	def mouseReleaseEvent(self, event):
		self.plot.mouseReleaseEvent(event)
		
	def wheelEvent(self, event):
		self.plot.wheelEvent(event)

	def get_render_dims_at_depth(self, depth):
		# This function returns the width and height of the renderable 
		# area at the origin of the data volume
		height = -2*tan(self.fov/2.0*pi/180.0)*(depth)
		width = self.aspect*height
		
		return [width,height]

class EMGLPlotInspector(QtGui.QWidget):
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
		
		self.tabwidget = QtGui.QTabWidget()
		self.maintab = None
		self.tabwidget.addTab(self.getMainTab(), "Main")
		self.tabwidget.addTab(self.getGLTab(),"GL")
		self.vbl.addWidget(self.tabwidget)
		self.n3_showing = False
		
		QtCore.QObject.connect(self.scale, QtCore.SIGNAL("valueChanged"), target.setScale)
		QtCore.QObject.connect(self.az, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
		QtCore.QObject.connect(self.alt, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
		QtCore.QObject.connect(self.phi, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
		QtCore.QObject.connect(self.cbb, QtCore.SIGNAL("currentIndexChanged(QString)"), target.setColor)
		QtCore.QObject.connect(self.src, QtCore.SIGNAL("currentIndexChanged(QString)"), self.set_src)
		QtCore.QObject.connect(self.x_trans, QtCore.SIGNAL("valueChanged(double)"), target.setCamX)
		QtCore.QObject.connect(self.y_trans, QtCore.SIGNAL("valueChanged(double)"), target.setCamY)
		QtCore.QObject.connect(self.z_trans, QtCore.SIGNAL("valueChanged(double)"), target.setCamZ)
		QtCore.QObject.connect(self.wiretog, QtCore.SIGNAL("toggled(bool)"), target.toggleWire)
		QtCore.QObject.connect(self.lighttog, QtCore.SIGNAL("toggled(bool)"), target.toggleLight)
		QtCore.QObject.connect(self.glcontrast, QtCore.SIGNAL("valueChanged"), target.setGLContrast)
		QtCore.QObject.connect(self.glbrightness, QtCore.SIGNAL("valueChanged"), target.setGLBrightness)
	
	def getGLTab(self):
		self.gltab = QtGui.QWidget()
		gltab = self.gltab
		
		gltab.vbl = QtGui.QVBoxLayout(self.gltab )
		gltab.vbl.setMargin(0)
		gltab.vbl.setSpacing(6)
		gltab.vbl.setObjectName("Main")
		
		self.glcontrast = ValSlider(gltab,(1.0,5.0),"GLShd:")
		self.glcontrast.setObjectName("GLShade")
		self.glcontrast.setValue(1.0)
		gltab.vbl.addWidget(self.glcontrast)
		
		self.glbrightness = ValSlider(gltab,(-1.0,0.0),"GLBst:")
		self.glbrightness.setObjectName("GLBoost")
		self.glbrightness.setValue(0.1)
		self.glbrightness.setValue(0.0)
		gltab.vbl.addWidget(self.glbrightness)
	
		return gltab
	
	def getMainTab(self):
		if ( self.maintab == None ):
			self.maintab = QtGui.QWidget()
			maintab = self.maintab
			maintab.vbl = QtGui.QVBoxLayout(self.maintab)
			maintab.vbl.setMargin(0)
			maintab.vbl.setSpacing(6)
			maintab.vbl.setObjectName("Main")
			
			self.scale = ValSlider(maintab,(0.01,30.0),"Zoom:")
			self.scale.setObjectName("scale")
			self.scale.setValue(1.0)
			maintab.vbl.addWidget(self.scale)
			
			self.hbl_color = QtGui.QHBoxLayout()
			self.hbl_color.setMargin(0)
			self.hbl_color.setSpacing(6)
			self.hbl_color.setObjectName("Material")
			maintab.vbl.addLayout(self.hbl_color)
			
			self.color_label = QtGui.QLabel()
			self.color_label.setText('Material')
			self.hbl_color.addWidget(self.color_label)
			
			self.cbb = QtGui.QComboBox(maintab)
			self.hbl_color.addWidget(self.cbb)
	
			self.hbl_trans = QtGui.QHBoxLayout()
			self.hbl_trans.setMargin(0)
			self.hbl_trans.setSpacing(6)
			self.hbl_trans.setObjectName("Trans")
			maintab.vbl.addLayout(self.hbl_trans)
			
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
			self.az = ValSlider(self,(-360.0,360.0),"az",-1)
			self.az.setObjectName("az")
			maintab.vbl.addWidget(self.az)
			
			self.alt = ValSlider(self,(-180.0,180.0),"alt",-1)
			self.alt.setObjectName("alt")
			maintab.vbl.addWidget(self.alt)
			
			self.phi = ValSlider(self,(-360.0,360.0),"phi",-1)
			self.phi.setObjectName("phi")
			maintab.vbl.addWidget(self.phi)
		
			self.current_src = EULER_EMAN
		
		return self.maintab

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

	def setScale(self,newscale):
		self.scale.setValue(newscale)
		
# This is just for testing, of course
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = EMGLPlotWidget()
	window.setInit()
	window2=EMParentWin(window)
	window2.show()
	
	sys.exit(app.exec_())
