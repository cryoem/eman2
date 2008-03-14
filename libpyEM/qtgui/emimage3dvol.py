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


class EMVolume(EMImage3DObject):
	def __init__(self,image=None, parent=None):
		EMImage3DObject.__init__(self)
		self.parent = parent
		
		self.init()
		self.initialized = True
		
		self.initializedGL= False
		
		self.inspector=None
		
		self.tex_names_list = []		# A storage object, used to remember and later delete texture names
		
		self.axes_idx = -1
		self.axes = []
		self.axes.append( Vec3f(1,0,0) )
		self.axes.append( Vec3f(0,1,0) )
		self.axes.append( Vec3f(0,0,1) )
		#self.axes.append( Vec3f(-1,0,0) )
		#self.axes.append( Vec3f(0,-1,0) )
		#self.addRenderAxis(1,1,1)
		#self.addRenderAxis(-1,1,1)
		#self.addRenderAxis(-1,-1,1)
		#self.addRenderAxis(1,-1,1)
		
		#self.addRenderAxis(1,1,0)
		#self.addRenderAxis(-1,1,0)
		#self.addRenderAxis(-1,-1,0)
		#self.addRenderAxis(1,-1,0)
		
		#self.addRenderAxis(0,1,1)
		#self.addRenderAxis(0,-1,1)
		
		#self.addRenderAxis(1,0,1)
		#self.addRenderAxis(-1,0,1)
	
		
		
		if image :
			self.setData(image)

	def addRenderAxis(self,a,b,c):
		v = Vec3f(a,b,c);
		v.normalize()
		self.axes.append( v )
	
	def getType(self):
		return "Volume"

	def init(self):
		self.data=None

		self.mmode=0
		
		self.cam = Camera()
		self.cube = False
		
		self.contrast = 1.0
		self.brightness = 0.0
		self.texsample = 1.0
		self.glcontrast = 1.0
		self.glbrightness = 0.0
		self.cube = False
		
		self.tex_name = 0
		
		self.rank = 1
		
		self.tex_dl = 0
		self.inspector=None
		
		self.force_texture_update = False

		self.glflags = EMOpenGLFlagsAndTools()		# OpenGL flags - this is a singleton convenience class for testing texture support
		
	def setData(self,data):
		"""Pass in a 3D EMData object"""
		
		self.data=data
		if data==None:
			print "Error, the data is empty"
			return
		
	 	if (isinstance(data,EMData) and data.get_zsize()<=1) :
			print "Error, the data is not 3D"
			return
		
		self.cam.default_z = -1.25*data.get_zsize()
		self.cam.cam_z = -1.25*data.get_zsize()
		
		if not self.inspector or self.inspector ==None:
			self.inspector=EMVolumeInspector(self)
		
		self.updateDataAndTexture()
		
	def test_accum(self):
		# this code will do volume rendering using the accumulation buffer
		# I opted not to go this way because you can't retain depth in the accumulation buffer
		# Note that it only works in the z-direction
		glClear(GL_ACCUM_BUFFER_BIT)
		
		self.accum = True
		self.zsample = self.texsample*(self.data.get_zsize())
		
		if self.tex_name == 0:
			print "Error, can not render 3D texture - texture name is 0"
			return
		
		
		for z in range(0,int(self.texsample*(self.data.get_zsize()))):
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
			
			zz = float(z)/float(self.data.get_zsize()-1)/self.texsample
			glTexCoord3f(0,0,zz)
			glVertex3f(0,0,zz)
			
			glTexCoord3f(1,0,zz)
			glVertex3f(1,0,zz)
			
			glTexCoord3f(1,1,zz)
			glVertex3f(1,1,zz)
			
			glTexCoord3f(0,1,zz)
			glVertex3f(0,1,zz)
		
			glEnd()
			glDisable(GL_TEXTURE_3D)
		
			if ( self.accum ):
				glAccum(GL_ADD, 1.0/self.zsample*self.brightness)
				glAccum(GL_ACCUM, 1.0/self.zsample*self.contrast)
		
		
		glAccum(GL_RETURN, 1.0)
		
	def render(self):
		lighting = glIsEnabled(GL_LIGHTING)
		cull = glIsEnabled(GL_CULL_FACE)
		depth = glIsEnabled(GL_DEPTH_TEST)
		
		polygonmode = glGetIntegerv(GL_POLYGON_MODE)

		glDisable(GL_LIGHTING)
		glDisable(GL_CULL_FACE)
		glDisable(GL_DEPTH_TEST)
		
		glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
		
		self.cam.position()

		# here is where the correct display list (x,y or z direction) is determined
		self.textureUpdateIfNecessary()

		glStencilFunc(GL_EQUAL,self.rank,0)
		glStencilOp(GL_KEEP,GL_KEEP,GL_REPLACE)
		glPushMatrix()
		glTranslate(-self.data.get_xsize()/2.0,-self.data.get_ysize()/2.0,-self.data.get_zsize()/2.0)
		glScalef(self.data.get_xsize(),self.data.get_ysize(),self.data.get_zsize())
		glEnable(GL_BLEND)
		#glBlendEquation(GL_MAX)
		if self.glflags.blend_equation_supported():
			glBlendEquation(GL_FUNC_ADD)
		glDepthMask(GL_FALSE)
		glBlendFunc(GL_ONE, GL_ONE)
		glCallList(self.tex_dl)
		glDepthMask(GL_TRUE)
		glDisable(GL_BLEND)
		glPopMatrix()

		# this is the accumulation buffer version of the volume renderer - it was for testing purposes
		# and is left here commented out incase anyone wants to investigate it in the future
		#glPushMatrix()
		#glTranslate(-self.data.get_xsize()/2.0,-self.data.get_ysize()/2.0,-self.data.get_zsize()/2.0)
		#glScalef(self.data.get_xsize(),self.data.get_ysize(),self.data.get_zsize())
		#self.test_accum()
		#glPopMatrix()
		
		
		glStencilFunc(GL_EQUAL,self.rank,self.rank)
		glStencilOp(GL_KEEP,GL_KEEP,GL_KEEP)
		glPushMatrix()
		glLoadIdentity()
		glTranslate(-self.data.get_xsize()/2.0,-self.data.get_ysize()/2.0,-10)
		glScalef(self.data.get_xsize(),self.data.get_ysize(),1)
		self.draw_bc_screen()
		glPopMatrix()
		
		glStencilFunc(GL_ALWAYS,1,1)
		if self.cube:
			glPushMatrix()
			self.draw_volume_bounds()
			glPopMatrix()
			
		if ( lighting ): glEnable(GL_LIGHTING)
		if ( cull ): glEnable(GL_CULL_FACE)
		if ( depth ): glEnable(GL_DEPTH_TEST)
		
		if ( polygonmode[0] == GL_LINE ): glPolygonMode(GL_FRONT, GL_LINE)
		if ( polygonmode[1] == GL_LINE ): glPolygonMode(GL_BACK, GL_LINE)
	
	def textureUpdateIfNecessary(self):
		
		t3d = self.cam.t3d_stack[len(self.cam.t3d_stack)-1]
		
		point = Vec3f(0,0,1)
		
		point = point*t3d
		
		point[0] = abs(point[0])
		point[1] = abs(point[1])
		point[2] = abs(point[2])
		#point[1] = -point[1]
		#if ( point[2] < 0 ):
			#point[2] = -point[2]
			#point[1] = -point[1]
			#point[0] = -point[0]
	
		currentaxis = self.axes_idx
		
		closest = 2*pi
		lp = point.length()
		idx = 0
		for i in self.axes:
			angle = abs(acos(point.dot(i)))
			if (angle < closest):
				closest = angle
				self.axes_idx = idx
			
			idx += 1

		if (currentaxis != self.axes_idx or self.force_texture_update):
			#print self.axes[self.axes_idx]
			self.genTexture()
			
	def genTexture(self):
		if ( self.glflags.threed_texturing_supported() ):
			self.gen3DTexture()
		else:
			self.gen2DTexture()
			
	def gen3DTexture(self):
	
		if ( self.tex_dl != 0 ): glDeleteLists( self.tex_dl, 1)
		
		self.tex_dl = glGenLists(1)
		
		if self.tex_dl == 0:
			print "Error, failed to generate display list"
			return
		
		if ( self.force_texture_update ):
			if self.tex_name != 0:
				glDeleteTextures(self.tex_name)
			
			self.tex_name = self.glflags.genTextureName(self.data_copy)
			
			self.force_texture_update = False
			
		
		if self.tex_name == 0:
			print "Error, can not render 3D texture - texture name is 0"
			return
		
		
		glNewList(self.tex_dl,GL_COMPILE)
		glEnable(GL_TEXTURE_3D)
		glBindTexture(GL_TEXTURE_3D, self.tex_name)
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1)
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE)
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP)
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP)
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP)
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
		if ( not data_dims_power_of(self.data_copy,2) and self.glflags.npt_textures_unsupported()):
			glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR)
		else:
			glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)

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
		for i in range(0,int(self.texsample*n)):
			nn = float(i)/float(n)/self.texsample

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
			return self.data_copy.get_zsize()
		elif ( self.axes_idx == 1 ):
			return self.data_copy.get_ysize()
		elif ( self.axes_idx == 2 ):
			return self.data_copy.get_xsize()
		else:
			#print "unsupported axis"
			# this is a hack and needs to be fixed eventually
			return self.data_copy.get_xsize()
			#return 0
		
	def getCorrectDims2DEMData(self):
		if ( self.axes_idx == 0 ):
			return EMData(self.data_copy.get_xsize(),self.data_copy.get_ysize())
		elif ( self.axes_idx == 1 ):
			return EMData(self.data_copy.get_xsize(),self.data_copy.get_zsize())
		elif ( self.axes_idx == 2 ):
			return EMData(self.data_copy.get_ysize(),self.data_copy.get_zsize())
		else:
			#print "unsupported axis"
			# this is a hack and needs to be fixed eventually
			return EMData(self.data_copy.get_xsize(),self.data_copy.get_zsize())

	
	def gen2DTexture(self):
			
		if ( self.tex_dl != 0 ): 
			glDeleteLists( self.tex_dl, 1)
		
		for i in self.tex_names_list:
			glDeleteTextures(i)
			
		self.tex_dl = glGenLists(1)
		if (self.tex_dl == 0 ):
			print "error, could not generate list"
			return

		glNewList(self.tex_dl,GL_COMPILE)
		glEnable(GL_TEXTURE_2D)

		n = self.getDimensionSize()
		v = self.axes[self.axes_idx]
		[t,alt,phi] = self.getEmanTransform(v)
		for i in range(0,int(self.texsample*n)):
			nn = float(i)/float(n)/self.texsample
			tmp = self.getCorrectDims2DEMData() 
			
			trans = (nn-0.5)*v
			t.set_posttrans(2.0*int(n/2)*trans)
			tmp.cut_slice(self.data_copy,t,True)
			#tmp.write_image("tmp.img",-1)
			
			# get the texture name, store it, and bind it in OpenGL
			tex_name = self.glflags.genTextureName(tmp)
			self.tex_names_list.append(tex_name)
			glBindTexture(GL_TEXTURE_2D, tex_name)
			
			self.loadDefault2DTextureParms()
			
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
		
		# this may have been toggled (i.e. if the image contrast or brightness changed)
		if self.force_texture_update == True:
			self.force_texture_update = False
	
	def loadDefault2DTextureParms(self):
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1)
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP)
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP)
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE)
		if ( not data_dims_power_of(self.data_copy,2) and self.glflags.npt_textures_unsupported()):
			glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR)
		else:
			glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR) 
		
	def updateDataAndTexture(self):
	
		if ( not isinstance(self.data,EMData) ): return
		
		self.data_copy = self.data.copy()
		self.data_copy.add(self.brightness)
		self.data_copy.mult(self.contrast*1.0/self.data.get_zsize())
		
		hist = self.data_copy.calc_hist(256,0,1.0)
		self.inspector.setHist(hist,0,1.0)

		self.force_texture_update = True

	def setContrast(self,val):
		self.contrast = val
		self.updateDataAndTexture()
		self.parent.updateGL()
		
	def setBrightness(self,val):
		self.brightness = val
		self.updateDataAndTexture()
		self.parent.updateGL()
		
	def setTextureSample(self,val):
		if ( val < 0 ) :
			print "Error, cannot handle texture sample less than 0"
			return
		
		self.texsample = val
		self.force_texture_update = True
		self.parent.updateGL()

	def updateInspector(self,t3d):
		if not self.inspector or self.inspector ==None:
			self.inspector=EMVolumeInspector(self)
		self.inspector.updateRotations(t3d)
	
	def getInspector(self):
		if not self.inspector : self.inspector=EMVolumeInspector(self)
		return self.inspector

class EMVolumeWidget(QtOpenGL.QGLWidget):
	
	allim=WeakKeyDictionary()
	def __init__(self, image=None, parent=None):
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True)
		fmt.setDepth(1)
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		self.initializeGL()
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
		glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST)
		
		glClearColor(0,0,0,0)
	
		glShadeModel(GL_SMOOTH)
	def paintGL(self):
		glClear(GL_ACCUM_BUFFER_BIT)
		glClear(GL_STENCIL_BUFFER_BIT)
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		glPushMatrix()
		self.volume.render()
		glPopMatrix()
	
		
	def resizeGL(self, width, height):
		if width<=0 or height<=0 : return
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
		e.process_inplace('testimage.axes')
		
		window2=EMParentWin(window)
		window2.show()
		window.setData(e)
	else :
		if not os.path.exists(sys.argv[1]):
			print "Error, input file %s does not exist" %sys.argv[1]
			exit(1)
		a=EMData.read_images(sys.argv[1],[0])
		window2=EMParentWin(window)
		window2.show()
		window.setData(a[0])
	
	
#	w2=QtGui.QWidget()
#	w2.resize(256,128)
	
#	w3=ValSlider(w2)
#	w3.resize(256,24)
#	w2.show()
	
	sys.exit(app.exec_())
