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
from emimageutil import ImgHistogram,EMParentWin
from weakref import WeakKeyDictionary
from time import time
from PyQt4.QtCore import QTimer

from time import *

from emglobjects import EMImage3DObject, Camera, EMOpenGLFlagsAndTools

MAG_INCREMENT_FACTOR = 1.1

class EM3DSliceViewer(EMImage3DObject):
	def __init__(self,image=None, parent=None):
		
		EMImage3DObject.__init__(self)
		
		self.parent = parent
		
		self.init()
		self.initialized = True
		
		self.inspector=None
		
		
		self.axes = []
		self.axes.append( Vec3f(1,0,0) )
		self.axes.append( Vec3f(0,1,0) )
		self.axes.append( Vec3f(0,0,1) )
		self.axes_idx = 2
		
		self.track = False
		
		if image :
			self.setData(image)
	
	def getType(self):
		return "Slice Viewer"

	def init(self):
		self.data=None

		self.mmode=0
		self.cam = Camera()
		
		self.cube = False
		self.inspector=None
		
		self.tex_name = 0
		self.tex_dl = 0

		self.glcontrast = 1.0
		self.glbrightness = 0.0
		
		self.rank = 1
		
		self.glflags = EMOpenGLFlagsAndTools()		# OpenGL flags - this is a singleton convenience class for testing texture support
		
	def setData(self,data,fact=1.0):
		"""Pass in a 3D EMData object"""
		
		if data==None:
			print "Error, the data is empty"
			return
		
		if (isinstance(data,EMData) and data.get_zsize()<=1) :
			print "Error, the data is not 3D"
			return
		
		self.data = data.copy()
		
		min = self.data.get_attr("minimum")
		max = self.data.get_attr("maximum")
		
		self.data.add(-min)
		self.data.mult(1/(max-min))
		
		self.cam.default_z = -1.25*data.get_zsize()
		self.cam.cam_z = -1.25*data.get_zsize()
		
		if not self.inspector or self.inspector ==None:
			self.inspector=EMSlice3DInspector(self)

		hist = self.data.calc_hist(256,0,1.0)
		self.inspector.setHist(hist,0,1.0) 

		self.slice = data.get_zsize()/2
		self.zslice = data.get_zsize()/2
		self.yslice = data.get_ysize()/2
		self.xslice = data.get_xsize()/2
		self.trackslice = data.get_xsize()/2
		self.axis = 'z'
		self.inspector.setSliceRange(0,data.get_zsize()-1)
		self.inspector.setSlice(self.zslice)
		self.genCurrentDisplayList()
		
	def getEmanTransform(self,p):
		
		if ( p[2] == 0 ):
			alt = 90
		else :
			alt = acos(p[2])*180.0/pi
		
		phi = atan2(p[0],p[1])
		phi *= 180.0/pi
		
		return [Transform3D(0,alt,phi),alt,phi]
			
	def getDimensionSize(self):
		if ( self.axes_idx == 0 ):
			return self.data.get_zsize()
		elif ( self.axes_idx == 1 ):
			return self.data.get_ysize()
		elif ( self.axes_idx == 2 ):
			return self.data.get_xsize()
		else:
			#print "unsupported axis"
			# this is a hack and needs to be fixed eventually
			return self.data.get_xsize()
			#return 0
	def getCorrectDims2DEMData(self):
		if ( self.axes_idx == 0 ):
			return EMData(self.data.get_xsize(),self.data.get_ysize())
		elif ( self.axes_idx == 1 ):
			return EMData(self.data.get_xsize(),self.data.get_zsize())
		elif ( self.axes_idx == 2 ):
			return EMData(self.data.get_ysize(),self.data.get_zsize())
		else:
			#print "unsupported axis"
			# this is a hack and needs to be fixed eventually
			return EMData(self.data.get_xsize(),self.data.get_zsize())

	def genCurrentDisplayList(self):
		
		if ( self.tex_dl != 0 ): glDeleteLists( self.tex_dl, 1)
		
		self.tex_dl = glGenLists(1)

		if (self.tex_dl == 0): return #OpenGL is initialized yet
		
		if ( self.glflags.threed_texturing_supported()):
			self.gen3DTexture()
		else:
			self.gen2DTexture()
	
	def gen3DTexture(self):

		glNewList(self.tex_dl,GL_COMPILE)
		
		if ( self.tex_name != 0 ): glDeleteTextures(self.tex_name)
		
		self.tex_name = self.glflags.genTextureName(self.data)
		
		glEnable(GL_TEXTURE_3D)
		glBindTexture(GL_TEXTURE_3D, self.tex_name)
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1)
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP)
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP)
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP)
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
		if ( not data_dims_power_of(self.data,2) and self.glflags.npt_textures_unsupported()):
			glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR)
		else:
			glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
			
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE)
		
		glPushMatrix()
		glTranslate(0.5,0.5,0.5)
		
		glBegin(GL_QUADS)
		
		n = self.getDimensionSize()
		v = self.axes[self.axes_idx]
		[t,alt,phi] = self.getEmanTransform(v)
		v1 = t*Vec3f(-0.5,-0.5,0)
		v2 = t*Vec3f(-0.5, 0.5,0)
		v3 = t*Vec3f( 0.5, 0.5,0)
		v4 = t*Vec3f( 0.5,-0.5,0)
		vecs = [v1,v2,v3,v4]
		nn = float(self.slice)/float(n)

		trans = (nn-0.5)*v
		
		for r in vecs:
		
			w = [r[0] + trans[0], r[1] + trans[1], r[2] + trans[2]]
			t = [w[0]+0.5,w[1]+0.5,w[2]+0.5]
			glTexCoord3fv(t)
			glVertex3fv(w)
			
		glEnd()
		glPopMatrix()
		glDisable(GL_TEXTURE_3D)
		
		glEndList()
		
	def gen2DTexture(self):		
		glNewList(self.tex_dl,GL_COMPILE)
		
		n = self.getDimensionSize()
		v = self.axes[self.axes_idx]
		
		[t,alt,phi] = self.getEmanTransform(v)
			
		nn = float(self.slice)/float(n)
		tmp = self.getCorrectDims2DEMData() 

		trans = (nn-0.5)*v
		t.set_posttrans(2.0*int(n/2)*trans)
		tmp.cut_slice(self.data,t,True)
			
		if ( self.tex_name != 0 ): glDeleteTextures(self.tex_name)
		
		self.tex_name = self.glflags.genTextureName(tmp)
		
		glEnable(GL_TEXTURE_2D)
		glBindTexture(GL_TEXTURE_2D, self.tex_name)
			
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1)
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP)
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP)
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
		if ( not data_dims_power_of(self.data,2) and self.glflags.npt_textures_unsupported()):
			glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR)
		else:
			glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE)
		
			
		glPushMatrix()
		glTranslate(trans[0]+0.5,trans[1]+0.5,trans[2]+0.5)
		glRotatef(-phi,0,0,1)
		glRotatef(-alt,1,0,0)
		glBegin(GL_QUADS)
		glTexCoord2f(0,0)
		glVertex2f(-0.5,-0.5)
		
		glTexCoord2f(1,0)
		glVertex2f( 0.5,-0.5)
		
		glTexCoord2f(1,1)
		glVertex2f( 0.5, 0.5)
		
		glTexCoord2f(0,1)
		glVertex2f(-0.5, 0.5)
		glEnd()
		glPopMatrix()
		
		glDisable(GL_TEXTURE_2D)
		glEndList()
		
	def render(self):
		lighting = glIsEnabled(GL_LIGHTING)
		cull = glIsEnabled(GL_CULL_FACE)
		polygonmode = glGetIntegerv(GL_POLYGON_MODE)
		glDisable(GL_LIGHTING)
		glDisable(GL_CULL_FACE)
		
		glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
		
		self.cam.position()
		
		if ( self.tex_dl == 0 ):
			self.genCurrentDisplayList()
		
		glStencilFunc(GL_EQUAL,self.rank,0)
		glStencilOp(GL_KEEP,GL_KEEP,GL_REPLACE)
		glPushMatrix()
		glTranslate(-self.data.get_xsize()/2.0,-self.data.get_ysize()/2.0,-self.data.get_zsize()/2.0)
		glScalef(self.data.get_xsize(),self.data.get_ysize(),self.data.get_zsize())
		glCallList(self.tex_dl)
		glPopMatrix()
		
		glStencilFunc(GL_EQUAL,self.rank,self.rank)
		glStencilOp(GL_KEEP,GL_KEEP,GL_KEEP)
		glPushMatrix()
		glLoadIdentity()
		glTranslate(-self.data.get_xsize()/2.0,-self.data.get_ysize()/2.0,-1)
		glScalef(self.data.get_xsize(),self.data.get_ysize(),1)
		#self.draw_bc_screen()
		glPopMatrix()
		
		glStencilFunc(GL_ALWAYS,1,1)
		glColor3f(1,1,1)
		if self.cube:
			glPushMatrix()
			self.draw_volume_bounds()
			glPopMatrix()
			
		if ( lighting ): glEnable(GL_LIGHTING)
		if ( cull ): glEnable(GL_CULL_FACE)
		
		if ( polygonmode[0] == GL_LINE ): glPolygonMode(GL_FRONT, GL_LINE)
		if ( polygonmode[1] == GL_LINE ): glPolygonMode(GL_BACK, GL_LINE)
	
	def setSlice(self,val):
		self.slice = val
		if self.axis == 'z':
			self.zslice = val
		elif self.axis == 'y':
			self.yslice = val
		elif self.axis == 'x':
			self.xslice = val
		else:
			self.trackslice = val
		
		self.genCurrentDisplayList()
		self.parent.updateGL()
		
	def setAxis(self,val):
		self.axis = str(val).strip()
		
		if (self.inspector != None):
			if self.axis == 'z':
				self.inspector.setSliceRange(0,self.data.get_zsize()-1)
				self.inspector.setSlice(self.zslice)
				self.axes_idx = 2
				self.track = False
			elif self.axis == 'y':
				self.inspector.setSliceRange(0,self.data.get_ysize()-1)
				self.inspector.setSlice(self.yslice)
				self.axes_idx = 1
				self.track = False
			elif self.axis == 'x':
				self.inspector.setSliceRange(0,self.data.get_xsize()-1)
				self.inspector.setSlice(self.xslice)
				self.axes_idx = 0
				self.track = False
			elif self.axis == 'track':
				self.track = True
				self.inspector.setSliceRange(0,self.data.get_xsize()-1)
				self.inspector.setSlice(self.trackslice)
				self.axes_idx = 3
				
				self.loadTrackAxis()
			else:
				print "Error, unknown axis", self.axis, val
		
		self.genCurrentDisplayList()
		self.parent.updateGL()

	def updateInspector(self,t3d):
		if not self.inspector or self.inspector ==None:
			self.inspector=EMSlice3DInspector(self)
		self.inspector.updateRotations(t3d)
	
	def getInspector(self):
		if not self.inspector : self.inspector=EMSlice3DInspector(self)
		return self.inspector

	def loadTrackAxis(self):
		t3d = self.cam.t3d_stack[len(self.cam.t3d_stack)-1]
		
		point = Vec3f(0,0,1)
		
		point = point*t3d
		if ( point[1] != 0 ): point[1] = -point[1]
		
		if len(self.axes) == 3 :
			self.axes.append(point)
		else:
			self.axes[3] = point
	
	def mouseMoveEvent(self, event):
		if ( self.track ):
			self.loadTrackAxis()
		self.genCurrentDisplayList()
		EMImage3DObject.mouseMoveEvent(self,event)
		
class EMSliceViewerWidget(QtOpenGL.QGLWidget):
	
	allim=WeakKeyDictionary()
	def __init__(self, image=None, parent=None):
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True)
		fmt.setDepth(1)
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		EMSliceViewerWidget.allim[self]=0
		
		self.fov = 50 # field of view angle used by gluPerspective
		
		self.sliceviewer = EM3DSliceViewer(image,self)
	def setData(self,data):
		self.sliceviewer.setData(data)
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
		#GL.glClearAccum(0,0,0,0)
	
		glShadeModel(GL_SMOOTH)
		
		glClearStencil(0)
		glEnable(GL_STENCIL_TEST)
		
	def paintGL(self):
		#glClear(GL_ACCUM_BUFFER_BIT)
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT )
		
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		glPushMatrix()
		self.sliceviewer.render()
		glPopMatrix()
		
		#glAccum(GL_ADD, self.sliceviewer.glbrightness)
		#glAccum(GL_ACCUM, self.sliceviewer.glcontrast)
		#glAccum(GL_RETURN, 1.0)
		
	def resizeGL(self, width, height):
		# just use the whole window for rendering
		glViewport(0,0,self.width(),self.height())
		
		# maintain the aspect ratio of the window we have
		self.aspect = float(width)/float(height)
		
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		# using gluPerspective for simplicity
		gluPerspective(self.fov,self.aspect,1,5000)
		
		# switch back to model view mode
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		self.sliceviewer.resizeEvent()

	def showInspector(self,force=0):
		self.sliceviewer.showInspector(self,force)

	def closeEvent(self,event) :
		self.sliceviewer.closeEvent(event)
		
	def mousePressEvent(self, event):
		self.sliceviewer.mousePressEvent(event)
		self.emit(QtCore.SIGNAL("mousedown"), event)
		
	def mouseMoveEvent(self, event):
		self.sliceviewer.mouseMoveEvent(event)
		self.emit(QtCore.SIGNAL("mousedrag"), event)
	
	def mouseReleaseEvent(self, event):
		self.sliceviewer.mouseReleaseEvent(event)
		self.emit(QtCore.SIGNAL("mouseup"), event)
			
	def wheelEvent(self, event):
		self.sliceviewer.wheelEvent(event)

	def get_render_dims_at_depth(self,depth):
		# This function returns the width and height of the renderable 
		# area at the origin of the data volume
		height = -2*tan(self.fov/2.0*pi/180.0)*(depth)
		width = self.aspect*height
		
		return [width,height]
		

class EMSlice3DInspector(QtGui.QWidget):
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
		
		self.vbl.addWidget(self.getMainTab())
		
		self.n3_showing = False
		
		self.current_src = EULER_EMAN
		
		QtCore.QObject.connect(self.scale, QtCore.SIGNAL("valueChanged"), target.setScale)
		QtCore.QObject.connect(self.slice, QtCore.SIGNAL("valueChanged"), target.setSlice)
		QtCore.QObject.connect(self.glcontrast, QtCore.SIGNAL("valueChanged"), target.setGLContrast)
		QtCore.QObject.connect(self.glbrightness, QtCore.SIGNAL("valueChanged"), target.setGLBrightness)
		QtCore.QObject.connect(self.axisCombo, QtCore.SIGNAL("currentIndexChanged(QString)"), target.setAxis)
		QtCore.QObject.connect(self.az, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
		QtCore.QObject.connect(self.alt, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
		QtCore.QObject.connect(self.phi, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
		QtCore.QObject.connect(self.src, QtCore.SIGNAL("currentIndexChanged(QString)"), self.set_src)
		QtCore.QObject.connect(self.x_trans, QtCore.SIGNAL("valueChanged(double)"), target.setCamX)
		QtCore.QObject.connect(self.y_trans, QtCore.SIGNAL("valueChanged(double)"), target.setCamY)
		QtCore.QObject.connect(self.z_trans, QtCore.SIGNAL("valueChanged(double)"), target.setCamZ)
		QtCore.QObject.connect(self.cubetog, QtCore.SIGNAL("toggled(bool)"), target.toggleCube)
		QtCore.QObject.connect(self.defaults, QtCore.SIGNAL("clicked(bool)"), self.setDefaults)
		#QtCore.QObject.connect(self.cbb, QtCore.SIGNAL("currentIndexChanged(QString)"), target.setColor)
		
	def setDefaults(self):
		self.x_trans.setValue(0.0)
		self.y_trans.setValue(0.0)
		self.z_trans.setValue(0.0)
		self.scale.setValue(1.0)
		self.glcontrast.setValue(1.0)
		self.glbrightness.setValue(0.0)
		
		self.az.setValue(0.0)
		self.alt.setValue(0.0)
		self.phi.setValue(0.0)

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
		
		self.hbl_slice = QtGui.QHBoxLayout()
		self.hbl_slice.setMargin(0)
		self.hbl_slice.setSpacing(6)
		self.hbl_slice.setObjectName("Axis")
		maintab.vbl.addLayout(self.hbl_slice)
		
		self.slice = ValSlider(maintab,(0.0,10.0),"Slice:")
		self.slice.setObjectName("slice")
		self.slice.setValue(1.0)
		self.hbl_slice.addWidget(self.slice)
		
		self.axisCombo = QtGui.QComboBox(maintab)
		self.axisCombo.addItem(' z ')
		self.axisCombo.addItem(' y ')
		self.axisCombo.addItem(' x ')
		self.axisCombo.addItem(' track ')
		self.hbl_slice.addWidget(self.axisCombo)
		
		self.glcontrast = ValSlider(maintab,(1.0,5.0),"GLShd:")
		self.glcontrast.setObjectName("GLShade")
		self.glcontrast.setValue(1.0)
		maintab.vbl.addWidget(self.glcontrast)
		
		self.glbrightness = ValSlider(maintab,(-1.0,0.0),"GLBst:")
		self.glbrightness.setObjectName("GLBoost")
		self.glbrightness.setValue(0.1)
		self.glbrightness.setValue(0.0)
		maintab.vbl.addWidget(self.glbrightness)
		
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
			self.maintab.vbl.removeWidget(self.n3)
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
			self.maintab.vbl.addWidget(self.n3)
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
		
	def setSlice(self,val):
		self.slice.setValue(val)
	
	def setSliceRange(self,min,max):
		self.slice.setRange(min,max)
	
		
# This is just for testing, of course
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = EMSliceViewerWidget()
 	if len(sys.argv)==1 : 
		e = EMData()
		e.set_size(7,7,7)
		e.process_inplace('testimage.axes')
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
	window2=EMParentWin(window)
	window2.show()
	
#	w2=QtGui.QWidget()
#	w2.resize(256,128)
	
#	w3=ValSlider(w2)
#	w3.resize(256,24)
#	w2.show()
	
	sys.exit(app.exec_())
