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

from emglobjects import Camera2, EMImage3DObject, EMViewportDepthTools, Camera2, Camera

MAG_INCREMENT_FACTOR = 1.1

class EMIsosurface(EMImage3DObject):
	def eye_coords_dif(self,x1,y1,x2,y2,mdepth=True):
		return self.vdtools.eye_coords_dif(x1,y1,x2,y2,mdepth)

	def viewportHeight(self):
		return self.parent.height()
	
	def viewportWidth(self):
		return self.parent.width()
	
	def __init__(self,image=None, parent=None):
		EMImage3DObject.__init__(self)
		self.parent = parent
		
		self.init()
		self.initialized = True
		
		self.cam=Camera2(self)
		self.tex_name = 0
		self.texture = False

		self.brightness = 0
		self.contrast = 10
		self.glcontrast = 1.0
		self.glbrightness = 0.0
		self.rank = 1
		self.inspector=None
		self.data = None
		self.data_copy = None
		
		self.vdtools = EMViewportDepthTools(self)
		
		if image :
			self.set_data(image)
	
	def getType(self):
		return "Isosurface"
	
	def updateDataAndTexture(self):
		
		self.data_copy = self.data.copy()
		self.data_copy.add(self.brightness)
		self.data_copy.mult(self.contrast)
		
		hist = self.data_copy.calc_hist(256,self.minden,self.maxden)
		self.inspector.set_hist(hist,self.minden,self.maxden) 

		if ( self.texture ): self.genTexture()
	
	def genTexture(self):
		if ( self.texture == False ): return
		if ( self.tex_name != 0 ):
			glDeleteTextures(self.tex_name)
		
		if ( self.data_copy == None ):
			self.tex_name = self.data.gen_gl_texture()
		else:
			self.tex_name = self.data_copy.gen_gl_texture()
			
	def updateGL(self):
		try: self.parent.updateGL()
		except: pass
	
	def render(self):
		if (not isinstance(self.data,EMData)): return
		
		#a = time()
		lighting = glIsEnabled(GL_LIGHTING)
		cull = glIsEnabled(GL_CULL_FACE)
		depth = glIsEnabled(GL_DEPTH_TEST)
		polygonmode = glGetIntegerv(GL_POLYGON_MODE)
		normalize = glIsEnabled(GL_NORMALIZE)
		
		
		glEnable(GL_CULL_FACE)
		glCullFace(GL_BACK)
		glEnable(GL_DEPTH_TEST)
		glEnable(GL_NORMALIZE)
		#glDisable(GL_NORMALIZE)
		if ( self.wire ):
			glPolygonMode(GL_FRONT,GL_LINE);
		else:
			glPolygonMode(GL_FRONT,GL_FILL);
		
		if self.light:
			glEnable(GL_LIGHTING)
		else:
			glDisable(GL_LIGHTING)

		
		glPushMatrix()
		self.cam.position(True)
		# the ones are dummy variables atm... they don't do anything
		self.vdtools.update(1,1)
		glPopMatrix()
		
		self.cam.position()
		glShadeModel(GL_SMOOTH)
		if ( self.isodl == 0 ):
			self.getIsoDL()
		glStencilFunc(GL_EQUAL,self.rank,0)
		glStencilOp(GL_KEEP,GL_KEEP,GL_REPLACE)
		glMaterial(GL_FRONT, GL_AMBIENT, self.colors[self.isocolor]["ambient"])
		glMaterial(GL_FRONT, GL_DIFFUSE, self.colors[self.isocolor]["diffuse"])
		glMaterial(GL_FRONT, GL_SPECULAR, self.colors[self.isocolor]["specular"])
		glMaterial(GL_FRONT, GL_SHININESS, self.colors[self.isocolor]["shininess"])
		glMaterial(GL_FRONT, GL_EMISSION, self.colors[self.isocolor]["emission"])
		glColor(self.colors[self.isocolor]["ambient"])
		glPushMatrix()
		glTranslate(-self.data.get_xsize()/2.0,-self.data.get_ysize()/2.0,-self.data.get_zsize()/2.0)
		if ( self.texture ):
			glScalef(self.data.get_xsize(),self.data.get_ysize(),self.data.get_zsize())
		glCallList(self.isodl)
		glPopMatrix()
		
		glStencilFunc(GL_EQUAL,self.rank,self.rank)
		glStencilOp(GL_KEEP,GL_KEEP,GL_KEEP)
		glPushMatrix()
		glLoadIdentity()
		[width,height] = self.parent.get_near_plane_dims()
		z = self.parent.get_start_z()
		
		glTranslate(-width/2.0,-height/2.0,-z-0.01)
		glScalef(width,height,1.0)
		self.draw_bc_screen()
		glPopMatrix()
		
		glStencilFunc(GL_ALWAYS,1,1)
		if self.cube:
			glPushMatrix()
			self.draw_volume_bounds()
			glPopMatrix()
			
		if ( lighting ): glEnable(GL_LIGHTING)
		else: glDisable(GL_LIGHTING)
		if ( not cull ): glDisable(GL_CULL_FACE)
		else: glDisable(GL_CULL_FACE)
		if ( depth ): glEnable(GL_DEPTH_TEST)
		else : glDisable(GL_DEPTH_TEST)
		
		if ( not normalize ): glDisable(GL_NORMALIZE)
		
		if ( polygonmode[0] == GL_LINE ): glPolygonMode(GL_FRONT, GL_LINE)
		else: glPolygonMode(GL_FRONT, GL_FILL)
		#if ( polygonmode[1] == GL_LINE ): glPolygonMode(GL_BACK, GL_LINE)
		#else: glPolygonMode(GL_BACK, GL_FILL)
		
		#print "total time is", time()-a
			
	def init(self):
		self.mmode = 0
		self.inspector=None
		self.isothr=0.5
		self.isorender=None
		self.isodl = 0
		self.smpval=-1
		self.griddl = 0
		self.scale = 1.0
		self.cube = False
		self.wire = False
		self.light = True
		
		
	def getIsoDL(self):
		# create the isosurface display list
		self.isorender.set_surface_value(self.isothr)
		self.isorender.set_sampling(self.smpval)
		
		if ( self.texture ):
			if ( self.tex_name == 0 ):
				self.updateDataAndTexture()
				
		if ( self.texture  ):
			self.isodl = self.isorender.get_isosurface_dl(self.tex_name)
		else:
			self.isodl = self.isorender.get_isosurface_dl(0)
		#time2 = clock()
		#dt1 = time2 - time1
		#print "It took %f to render the isosurface" %dt1
	
	def updateData(self,data):
		if data==None or (isinstance(data,EMData) and data.get_zsize()<=1) :
			print "Error, tried to set data that is invalid for EMIsosurface"
			return
		self.data=data
		self.isorender=MarchingCubes(data)
		self.getIsoDL()
		self.parent.updateGL()
	
	def set_data(self,data):
		"""Pass in a 3D EMData object"""
		
		if data==None or (isinstance(data,EMData) and data.get_zsize()<=1) :
			print "Error, tried to set data that is invalid for EMIsosurface"
			return
		self.data=data
		
		self.minden=data.get_attr("minimum")
		self.maxden=data.get_attr("maximum")
		mean=data.get_attr("mean")
		sigma=data.get_attr("sigma")
		
		if not self.inspector or self.inspector ==None:
			self.inspector=EMIsoInspector(self)
		
		hist = data.calc_hist(256,self.minden,self.maxden)
		self.inspector.set_hist(hist,self.minden,self.maxden) 
	
		self.inspector.setThrs(self.minden,self.maxden,mean+3.0*sigma)
		self.isothr = mean+3.0*sigma
		self.brightness = -self.isothr
		
		self.isorender=MarchingCubes(data)
		self.inspector.setSamplingRange(self.isorender.get_sampling_range())
		
		self.loadColors()
		self.inspector.setMaterials(self.colors,self.isocolor)
	
	def loadColors(self):
		self.colors = getColors()
		
		self.isocolor = "ruby"
	
	def getMaterial(self):
		return self.colors[self.isocolor]
	
	def setThr(self,val):
		if (self.isothr != val):
			self.isothr = val
			self.brightness = -val
			if ( self.texture ):
				self.updateDataAndTexture()
			self.getIsoDL()
			self.parent.updateGL()
	
	def setSample(self,val):
		if ( self.smpval != int(val)):
			# the minus two is here because the marching cubes thinks -1 is the high level of detail, 0 is the next best and  so forth
			# However the user wants the highest level of detail to be 1, and the next best to be 2 and then 3 etc
			self.smpval = int(val)-2
			self.getIsoDL()
			self.parent.updateGL()
	
	def setMaterial(self,val):
		#print val
		self.isocolor = str(val)
		self.updateGL()
		
	def toggleCube(self):
		self.cube = not self.cube
		self.updateGL()
	
	def toggleWire(self,val):
		self.wire = not self.wire
		self.updateGL()
		
	def toggleLight(self,val):
		self.light = not self.light
		self.updateGL()
	
	def toggleTexture(self):
		self.texture = not self.texture
		if ( self.texture ):
			self.updateDataAndTexture()
		
		self.getIsoDL()
		self.updateGL()
	
	def updateInspector(self,t3d):
		if not self.inspector or self.inspector ==None:
			self.inspector=EMIsoInspector(self)
		self.inspector.updateRotations(t3d)
	
	def get_inspector(self):
		if not self.inspector : self.inspector=EMIsoInspector(self)
		return self.inspector
		
	def setContrast(self,val):
		self.contrast = val
		self.updateDataAndTexture()
		self.updateGL()
		
	def setBrightness(self,val):
		self.brightness = val
		self.updateDataAndTexture()
		self.updateGL()
		
	def mousePressEvent(self, event):
#		lc=self.scrtoimg((event.x(),event.y()))
		if event.button()==Qt.MidButton:
			if not self.inspector or self.inspector ==None:
				return
			self.inspector.updateRotations(self.cam.t3d_stack[len(self.cam.t3d_stack)-1])
			self.resizeEvent()
			self.show_inspector(1)
		else:
			self.cam.mousePressEvent(event)
		
		self.updateGL()
		
	def mouseMoveEvent(self, event):
		self.cam.mouseMoveEvent(event)
		self.updateGL()
	
	def mouseReleaseEvent(self, event):
		self.cam.mouseReleaseEvent(event)
		self.updateGL()
			
	def wheelEvent(self, event):
		self.cam.wheelEvent(event)
		self.updateGL()
		
	def resizeEvent(self,width=0,height=0):
		self.vdtools.set_update_P_inv()

def getColors():
	ruby = {}
	ruby["ambient"] = [0.1745, 0.01175, 0.01175,1.0]
	ruby["diffuse"] = [0.61424, 0.04136, 0.04136,1.0]
	ruby["specular"] = [0.927811, 0.826959, 0.826959,1.0]
	ruby["shininess"] = 32
	ruby["emission"] = [0,0,0]
	
	emerald = {}
	emerald["ambient"] = [0.0215, 0.1745, 0.0215,1.0]
	emerald["diffuse"] = [0.07568, 0.61424,  0.07568,1.0]
	emerald["specular"] = [0.833, 0.927811, 0.833,1.0]
	emerald["shininess"] = 32
	emerald["emission"] = [0,0,0]
	
	pearl = {}
	pearl["ambient"] = [0.25, 0.20725, 0.20725,1.0]
	pearl["diffuse"] = [1.0, 0.829, 0.829,1.0]
	pearl["specular"] = [0.296648, 0.296648, 0.296648,1.0]
	pearl["shininess"] = 128.0
	pearl["emission"] = [0,0,0]
	
	silver = {}
	silver["ambient"] = [0.25, 0.25, 0.25,1.0]
	silver["diffuse"] = [0.4, 0.4, 0.4,1.0]
	silver["specular"] = [0.974597, 0.974597, 0.974597,1.0]
	silver["shininess"] = 4
	silver["emission"] = [0.1,0.1,0.1]
	
	gold = {}
	gold["ambient"] = [0.24725, 0.2245, 0.0645,1.0]
	gold["diffuse"] = [0.34615, 0.3143, 0.0903,1.0]
	gold["specular"] = [1.000, 0.9079885, 0.26086934,1.0]
	gold["shininess"] = 4
	gold["emission"] = [0,0,0]
	
	copper = {}
	copper["ambient"] = [0.2295, 0.08825, 0.0275,1.0]
	copper["diffuse"] = [0.5508, 0.2118, 0.066,1.0]
	copper["specular"] = [0.9, 0.5, 0.2,1.0]
	copper["shininess"] = 20.0
	copper["emission"] = [0,0,0]
	
	obsidian = {}
	obsidian["ambient"] = [0.05375,  0.05,     0.06625 ,1.0]
	obsidian["diffuse"] = [0.18275,  0.17,     0.22525,1.0]
	obsidian["specular"] = [0.66, 0.65, 0.69]
	obsidian["shininess"] = 128.0
	obsidian["emission"] = [0,0,0]
	
	turquoise = {}
	turquoise["ambient"] = [0.1, 0.18725, 0.1745 ,1.0]
	turquoise["diffuse"] = [0.396, 0.74151, 0.69102,1.0]
	turquoise["specular"] = [0.297254, 0.30829, 0.306678]
	turquoise["shininess"] = 128.0
	turquoise["emission"] = [0,0,0]
	
	yellow = {}
	yellow["ambient"] = [0.3, 0.3, 0.0,1]
	yellow["diffuse"] = [0.5, 0.5, 0.0,1]
	yellow["specular"] = [0.7, 0.7, 0.0,1]
	yellow["shininess"] =  60
	yellow["emission"] = [0,0,0]
	
	custom = {}
	custom["custom"] = [0.3, 0.3, 0.0,1]
	custom["custom"] = [0.5, 0.5, 0.0,1]
	custom["custom"] = [0.7, 0.7, 0.0,1]
	custom["custom"] =  60
	custom["emission"] = [0,0,0]
	
	colors = {}
	colors["ruby"] = ruby
	colors["emerald"] = emerald
	colors["pearl"] = pearl
	colors["silver"] = silver
	colors["gold"] = gold
	colors["copper"] = copper
	colors["obsidian"] = obsidian
	colors["turquoise"] = turquoise
	colors["yellow"] = yellow
	colors["custom"] = yellow
	
	return colors
		
class EMIsosurfaceWidget(QtOpenGL.QGLWidget):
	""" This class is not yet complete.
	A QT widget for rendering 3D EMData objects.
	"""
	allim=WeakKeyDictionary()
	def __init__(self, image=None, parent=None):
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True)
		fmt.setDepth(1)
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)

		EMIsosurfaceWidget.allim[self]=0
		self.isosurface = EMIsosurface(image,self)
		self.timer = QTimer()
		QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.timeout)

		self.cam = Camera()
		self.aspect=1.0
		self.fov = 50 # field of view angle used by gluPerspective
		self.startz = 1
		self.endz = 5000
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
		glEnable(GL_CULL_FACE)
		glCullFace(GL_BACK)
		# For the time being
		
	def paintGL(self):
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		self.cam.position()
		
		glPushMatrix()
		self.isosurface.render()
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
		gluPerspective(self.fov,self.aspect,self.startz,self.endz)
		
		# switch back to model view mode
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		self.isosurface.resizeEvent()
	
	def get_start_z(self):
		return self.startz
	
	def get_near_plane_dims(self):
		height = 2.0 * self.startz*tan(self.fov/2.0*pi/180.0)
		width = self.aspect * height
		return [width,height]

	def set_data(self,data):
		self.isosurface.set_data(data)
		self.cam.default_z = -1.25*data.get_zsize()
		self.cam.cam_z = -1.25*data.get_zsize()
	
	def show_inspector(self,force=0):
		self.isosurface.show_inspector()
	
	def closeEvent(self,event) :
		self.isosurface.closeEvent(event)
		
	def mousePressEvent(self, event):
		self.isosurface.mousePressEvent(event)
		
	def mouseMoveEvent(self, event):
		self.isosurface.mouseMoveEvent(event)
	
	def mouseReleaseEvent(self, event):
		self.isosurface.mouseReleaseEvent(event)
			
	def wheelEvent(self, event):
		self.isosurface.wheelEvent(event)

	def get_render_dims_at_depth(self, depth):
		# This function returns the width and height of the renderable 
		# area at the origin of the data volume
		height = -2*tan(self.fov/2.0*pi/180.0)*(depth)
		width = self.aspect*height
		return [width,height]

class EMIsoInspector(QtGui.QWidget):
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
		
		self.texturetog = QtGui.QPushButton("Texture")
		self.texturetog.setCheckable(1)
		self.vbl2.addWidget(self.texturetog)
		self.texture = False
		
		self.tabwidget = QtGui.QTabWidget()
		self.maintab = None
		self.tabwidget.addTab(self.getMainTab(), "Main")
		self.texturetab = None
		self.tabwidget.addTab(self.getGLTab(),"GL")
		self.tabwidget.addTab(self.getTextureTab(),"Texture")
		self.getTextureTab().setEnabled(False)
		self.vbl.addWidget(self.tabwidget)
		self.n3_showing = False
		
		QtCore.QObject.connect(self.scale, QtCore.SIGNAL("valueChanged"), target.set_scale)
		QtCore.QObject.connect(self.az, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
		QtCore.QObject.connect(self.alt, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
		QtCore.QObject.connect(self.phi, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
		QtCore.QObject.connect(self.thr, QtCore.SIGNAL("valueChanged"), self.onThrSlider)
		QtCore.QObject.connect(self.contrast, QtCore.SIGNAL("valueChanged"), target.setContrast)
		QtCore.QObject.connect(self.bright, QtCore.SIGNAL("valueChanged"), target.setBrightness)
		QtCore.QObject.connect(self.cbb, QtCore.SIGNAL("currentIndexChanged(QString)"), self.setMaterial)
		QtCore.QObject.connect(self.src, QtCore.SIGNAL("currentIndexChanged(QString)"), self.set_src)
		QtCore.QObject.connect(self.smp, QtCore.SIGNAL("valueChanged(int)"), target.setSample)
		QtCore.QObject.connect(self.x_trans, QtCore.SIGNAL("valueChanged(double)"), target.set_cam_x)
		QtCore.QObject.connect(self.y_trans, QtCore.SIGNAL("valueChanged(double)"), target.set_cam_y)
		QtCore.QObject.connect(self.z_trans, QtCore.SIGNAL("valueChanged(double)"), target.set_cam_z)
		QtCore.QObject.connect(self.wiretog, QtCore.SIGNAL("toggled(bool)"), target.toggleWire)
		QtCore.QObject.connect(self.lighttog, QtCore.SIGNAL("toggled(bool)"), target.toggleLight)
		QtCore.QObject.connect(self.texturetog, QtCore.SIGNAL("toggled(bool)"), self.toggleTexture)
		QtCore.QObject.connect(self.cubetog, QtCore.SIGNAL("toggled(bool)"), target.toggleCube)
		QtCore.QObject.connect(self.glcontrast, QtCore.SIGNAL("valueChanged"), target.setGLContrast)
		QtCore.QObject.connect(self.glbrightness, QtCore.SIGNAL("valueChanged"), target.setGLBrightness)
		
		QtCore.QObject.connect(self.ambient_tab.r, QtCore.SIGNAL("valueChanged"), self.updateMaterial)
		QtCore.QObject.connect(self.ambient_tab.g, QtCore.SIGNAL("valueChanged"), self.updateMaterial)
		QtCore.QObject.connect(self.ambient_tab.b, QtCore.SIGNAL("valueChanged"), self.updateMaterial)
		QtCore.QObject.connect(self.diffuse_tab.r, QtCore.SIGNAL("valueChanged"), self.updateMaterial)
		QtCore.QObject.connect(self.diffuse_tab.g, QtCore.SIGNAL("valueChanged"), self.updateMaterial)
		QtCore.QObject.connect(self.diffuse_tab.b, QtCore.SIGNAL("valueChanged"), self.updateMaterial)
		QtCore.QObject.connect(self.specular_tab.r, QtCore.SIGNAL("valueChanged"), self.updateMaterial)
		QtCore.QObject.connect(self.specular_tab.g, QtCore.SIGNAL("valueChanged"), self.updateMaterial)
		QtCore.QObject.connect(self.specular_tab.b, QtCore.SIGNAL("valueChanged"), self.updateMaterial)
		QtCore.QObject.connect(self.emission_tab.r, QtCore.SIGNAL("valueChanged"), self.updateMaterial)
		QtCore.QObject.connect(self.emission_tab.g, QtCore.SIGNAL("valueChanged"), self.updateMaterial)
		QtCore.QObject.connect(self.emission_tab.b, QtCore.SIGNAL("valueChanged"), self.updateMaterial)
		QtCore.QObject.connect(self.shininess, QtCore.SIGNAL("valueChanged"), self.updateMaterial)
		
	def updateMaterial(self):
		self.target.isocolor = "custom"
		custom = {}
		
		custom["ambient"] = [self.ambient_tab.r.getValue(), self.ambient_tab.g.getValue(), self.ambient_tab.b.getValue(),1.0]
		custom["diffuse"] = [self.diffuse_tab.r.getValue(), self.diffuse_tab.g.getValue(), self.diffuse_tab.b.getValue(),1.0]
		custom["specular"] = [self.specular_tab.r.getValue(), self.specular_tab.g.getValue(), self.specular_tab.b.getValue(),1.0]
		custom["emission"] = [self.emission_tab.r.getValue(), self.emission_tab.g.getValue(), self.emission_tab.b.getValue(),1.0]
		custom["shininess"] = self.shininess.getValue()
		self.target.colors["custom"] = custom
		
		n = self.cbb.findText(QtCore.QString("custom"))
		if n < 0: return
		self.cbb.setCurrentIndex(n)
		self.target.updateGL()
	
	def setMaterial(self,color):
		self.target.setMaterial(color)
		material = self.target.getMaterial()
		
		self.ambient_tab.r.setValue(material["ambient"][0])
		self.ambient_tab.g.setValue(material["ambient"][1])
		self.ambient_tab.b.setValue(material["ambient"][2])
		
		self.diffuse_tab.r.setValue(material["diffuse"][0])
		self.diffuse_tab.g.setValue(material["diffuse"][1])
		self.diffuse_tab.b.setValue(material["diffuse"][2])
		
		self.specular_tab.r.setValue(material["specular"][0])
		self.specular_tab.g.setValue(material["specular"][1])
		self.specular_tab.b.setValue(material["specular"][2])
		
		self.emission_tab.r.setValue(material["emission"][0])
		self.emission_tab.g.setValue(material["emission"][1])
		self.emission_tab.b.setValue(material["emission"][2])
		
		self.shininess.setValue(material["shininess"])
	
	def getRGBTab(self, name=""):
		rgbtab = QtGui.QWidget(self)
		rgbtab.vbl = QtGui.QVBoxLayout(rgbtab)
		rgbtab.vbl.setMargin(0)
		rgbtab.vbl.setSpacing(6)
		rgbtab.vbl.setObjectName(name)
		
		rgbtab.r = ValSlider(rgbtab,(0.0,1.0),"R:")
		rgbtab.r.setObjectName("R")
		rgbtab.r.setValue(0.5)
		rgbtab.vbl.addWidget(rgbtab.r)
		
		rgbtab.g = ValSlider(rgbtab,(0.0,1.0),"G:")
		rgbtab.g.setObjectName("G")
		rgbtab.g.setValue(0.5)
		rgbtab.vbl.addWidget(rgbtab.g)
		
		rgbtab.b = ValSlider(rgbtab,(0.0,1.0),"B:")
		rgbtab.b.setObjectName("B")
		rgbtab.b.setValue(0.5)
		rgbtab.vbl.addWidget(rgbtab.b)
		
		return rgbtab
	
	def getGLTab(self):
		self.gltab = QtGui.QWidget()
		gltab = self.gltab
		
		gltab.vbl = QtGui.QVBoxLayout(self.gltab )
		gltab.vbl.setMargin(0)
		gltab.vbl.setSpacing(6)
		gltab.vbl.setObjectName("GL")
		
		
		self.glcontrast = ValSlider(gltab,(1.0,5.0),"GLShd:")
		self.glcontrast.setObjectName("GLShade")
		self.glcontrast.setValue(1.0)
		gltab.vbl.addWidget(self.glcontrast)
		
		self.glbrightness = ValSlider(gltab,(-1.0,0.0),"GLBst:")
		self.glbrightness.setObjectName("GLBoost")
		self.glbrightness.setValue(0.1)
		self.glbrightness.setValue(0.0)
		gltab.vbl.addWidget(self.glbrightness)
	
		self.material_tab_widget = QtGui.QTabWidget()
		self.ambient_tab = self.getRGBTab("ambient")
		self.material_tab_widget.addTab(self.ambient_tab, "Ambient")
		
		self.diffuse_tab = self.getRGBTab("diffuse")
		self.material_tab_widget.addTab(self.diffuse_tab, "Diffuse")
		
		self.specular_tab = self.getRGBTab("specular")
		self.material_tab_widget.addTab(self.specular_tab, "Specular")
		
		self.emission_tab = self.getRGBTab("emission")
		self.material_tab_widget.addTab(self.emission_tab, "Emission")
		
		gltab.vbl.addWidget(self.material_tab_widget)

		self.shininess = ValSlider(gltab,(0,128),"Shininess:")
		self.shininess.setObjectName("Shininess")
		self.shininess.setValue(64)
		gltab.vbl.addWidget(self.shininess)

		self.hbl_color = QtGui.QHBoxLayout()
		self.hbl_color.setMargin(0)
		self.hbl_color.setSpacing(6)
		self.hbl_color.setObjectName("Material")
		gltab.vbl.addLayout(self.hbl_color)
		
		self.color_label = QtGui.QLabel()
		self.color_label.setText('Material')
		self.hbl_color.addWidget(self.color_label)
		
		self.cbb = QtGui.QComboBox(gltab)
		self.hbl_color.addWidget(self.cbb)
		
		return gltab
	
	def toggleTexture(self):
		self.texture = not self.texture
		self.target.toggleTexture()
		self.getTextureTab().setEnabled(self.texture)
	
	def getTextureTab(self):
		if ( self.texturetab == None ):
			self.texturetab = QtGui.QWidget()
			texturetab = self.texturetab
			texturetab.vbl = QtGui.QVBoxLayout(self.texturetab)
			texturetab.vbl.setMargin(0)
			texturetab.vbl.setSpacing(6)
			texturetab.vbl.setObjectName("Main")
		
			self.contrast = ValSlider(texturetab,(0.0,20.0),"Cont:")
			self.contrast.setObjectName("contrast")
			self.contrast.setValue(10.0)
			texturetab.vbl.addWidget(self.contrast)
	
			self.bright = ValSlider(texturetab,(-5.0,5.0),"Brt:")
			self.bright.setObjectName("bright")
			self.bright.setValue(0.1)
			self.bright.setValue(0.0)
			texturetab.vbl.addWidget(self.bright)
			
			#self.glcontrast = ValSlider(texturetab,(1.0,5.0),"GLShd:")
			#self.glcontrast.setObjectName("GLShade")
			#self.glcontrast.setValue(1.0)
			#texturetab.vbl.addWidget(self.glcontrast)
			
			#self.glbrightness = ValSlider(texturetab,(-1.0,0.0),"GLBst:")
			#self.glbrightness.setObjectName("GLBoost")
			#self.glbrightness.setValue(0.1)
			#self.glbrightness.setValue(0.0)
			#texturetab.vbl.addWidget(self.glbrightness)
			
		return self.texturetab
	
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
			
			self.thr = ValSlider(maintab,(0.0,4.0),"Thr:")
			self.thr.setObjectName("thr")
			self.thr.setValue(0.5)
			maintab.vbl.addWidget(self.thr)
			
			self.hbl_smp = QtGui.QHBoxLayout()
			self.hbl_smp.setMargin(0)
			self.hbl_smp.setSpacing(6)
			self.hbl_smp.setObjectName("Sample")
			maintab.vbl.addLayout(self.hbl_smp)
			
			self.smp_label = QtGui.QLabel()
			self.smp_label.setText('Sample Level')
			self.hbl_smp.addWidget(self.smp_label)
			
			self.smp = QtGui.QSpinBox(maintab)
			self.smp.setValue(1)
			self.hbl_smp.addWidget(self.smp)
	
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
			self.az.setValue(0.0)
			maintab.vbl.addWidget(self.az)
			
			self.alt = ValSlider(self,(-180.0,180.0),"alt",-1)
			self.alt.setObjectName("alt")
			self.alt.setValue(0.0)
			maintab.vbl.addWidget(self.alt)
			
			self.phi = ValSlider(self,(-360.0,360.0),"phi",-1)
			self.phi.setObjectName("phi")
			self.phi.setValue(0.0)
			maintab.vbl.addWidget(self.phi)
		
			self.current_src = EULER_EMAN
		
		return self.maintab
	
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
		self.target.load_rotation(self.getCurrentRotation())
	
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
		
	
	def setMaterials(self,colors,current_color):
		a = 0
		for i in colors:
			self.cbb.addItem(i)
			if ( i == current_color):
				self.cbb.setCurrentIndex(a)
			a += 1

	def onThrSlider(self,val):
		self.target.setThr(val)
		self.bright.setValue(-val,True)
		
	def setThrs(self,low,high,val):
		self.thr.setRange(low,high)
		self.thr.setValue(val, True)
		self.bright.setValue(-val,True)
	
	def setSamp(self,low,high,val):
		self.smp.setRange(int(low),int(high))
		self.smp.setValue(val, True)
		
	def set_hist(self,hist,minden,maxden):
		self.hist.set_data(hist,minden,maxden)

	def set_scale(self,newscale):
		self.scale.setValue(newscale)
		
# This is just for testing, of course
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = EMIsosurfaceWidget()
 	if len(sys.argv)==1 : 
		e = EMData()
		e.set_size(40,35,30)
		e.process_inplace('testimage.axes')
 		window.set_data(e)

	else :
		if not os.path.exists(sys.argv[1]):
			print "Error, input file %s does not exist" %sys.argv[1]
			exit(1)
		a=EMData.read_images(sys.argv[1],[0])
		window.set_data(a[0])
	window2=EMParentWin(window)
	window2.show()
	
	sys.exit(app.exec_())
