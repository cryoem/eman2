#!/usr/bin/env python

#
# Author: David Woolford 10/26/2007 (woolford@bcm.edu)
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
from weakref import WeakKeyDictionary
from emimageutil import ImgHistogram,EMParentWin
import copy

from emimage3dobject import Camera2

# this class tracks depth in opengl
class DepthTracker:
	def __init__(self,parent):
		self.parent = parent
		self.update_P_inv = True
		
	def set_update_P_inv(self,val=True):
		self.update_P_inv = val
		
	def storeMatrices(self):
		self.wmodel = glGetDoublev(GL_MODELVIEW_MATRIX)
		if self.update_P_inv == True:
			self.wproj = glGetDoublev(GL_PROJECTION_MATRIX)
		self.wview = glGetIntegerv(GL_VIEWPORT)

	def eyeCoordsDif(self,x1,y1,x2,y2,maintaindepth=True):
		# get x and y normalized device coordinates
		xNDC1 = 2.0*(x1-self.wview[0])/self.wview[2] - 1
		yNDC1 = 2.0*(y1-self.wview[1])/self.wview[3] - 1
		
		xNDC2 = 2.0*(x2-self.wview[0])/self.wview[2] - 1
		yNDC2 = 2.0*(y2-self.wview[1])/self.wview[3] - 1
		
		## invert the projection and model view matrices, they will be used shortly
		## note the OpenGL returns matrices are in column major format -  the calculations below 
		## are done with this in  mind - this saves the need to transpose the matrices
		
		if ( self.update_P_inv == True):
			P = numpy.matrix(self.wproj)
			self.update_P_inv = False
			self.P_inv = P.I
		
		M = numpy.matrix(self.wmodel)
		M_inv = M.I
		
		#PM_inv = numpy.matrixmultiply(P_inv,M_inv)
		PM_inv = self.P_inv*M_inv
		# If the widget is planar (which obviosuly holds), and along z=0, then the following holds
		zNDC1 = (PM_inv[0,2]*xNDC1 + PM_inv[1,2]*yNDC1 + PM_inv[3,2])/(-PM_inv[2,2])
		if ( maintaindepth == False):
			zNDC2 = (PM_inv[0,2]*xNDC2 + PM_inv[1,2]*yNDC2 + PM_inv[3,2])/(-PM_inv[2,2])
		else:
			zNDC2 = zNDC1
	
		# We need zprime, which is really 'eye_z' in OpenGL lingo
		zprime1 = 1.0/(xNDC1*self.P_inv[0,3]+yNDC1*self.P_inv[1,3]+zNDC1*self.P_inv[2,3]+self.P_inv[3,3])
		zprime2 = 1.0/(xNDC2*self.P_inv[0,3]+yNDC2*self.P_inv[1,3]+zNDC2*self.P_inv[2,3]+self.P_inv[3,3])

		ex1 = (self.P_inv[0,0]*xNDC1 + self.P_inv[1,0]*yNDC1 + self.P_inv[2,0]*zNDC1+self.P_inv[3,0])*zprime1;
		ey1 = (self.P_inv[0,1]*xNDC1 + self.P_inv[1,1]*yNDC1 + self.P_inv[2,1]*zNDC1+self.P_inv[3,1])*zprime1;
		#ez1 = (self.P_inv[0,2]*xNDC1 + self.P_inv[1,2]*yNDC1 + self.P_inv[2,2]*zNDC1+self.P_inv[3,2])*zprime1;
		
		ex2 = (self.P_inv[0,0]*xNDC2 + self.P_inv[1,0]*yNDC2 + self.P_inv[2,0]*zNDC2+self.P_inv[3,0])*zprime2;
		ey2 = (self.P_inv[0,1]*xNDC2 + self.P_inv[1,1]*yNDC2 + self.P_inv[2,1]*zNDC2+self.P_inv[3,1])*zprime2;
		
		return [ex2-ex1,ey2-ey1]

class BrightContrastScreen:
	def __init__(self):
		# this class draws a brightness/contrast screen on the zplane,
		# on a square polygon from [0,0] to [1,1]
		self.glcontrast = 1.0
		self.glbrightness = 0.0

	def setGLBrightness(self,val):
		self.glbrightness = val
		
	def setGLContrast(self,val):
		self.glcontrast = val

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
		
		glShadeModel(GL_SMOOTH)

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

class EMImage2DTexInspector(QtGui.QWidget):
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
		
		self.defaults = QtGui.QPushButton("A Button")
		self.vbl2.addWidget(self.defaults)
	
		self.glcontrast = ValSlider(self,(1.0,5.0),"GLShd:")
		self.glcontrast.setObjectName("GLShade")
		self.glcontrast.setValue(1.0)
		self.vbl.addWidget(self.glcontrast)
		
		self.glbrightness = ValSlider(self,(-1.0,0.0),"GLBst:")
		self.glbrightness.setObjectName("GLBoost")
		self.glbrightness.setValue(0.1)
		self.glbrightness.setValue(0.0)
		self.vbl.addWidget(self.glbrightness)
	
		QtCore.QObject.connect(self.glcontrast, QtCore.SIGNAL("valueChanged"), target.setGLContrast)
		QtCore.QObject.connect(self.glbrightness, QtCore.SIGNAL("valueChanged"), target.setGLBrightness)
		
	def setHist(self,hist,minden,maxden):
		self.hist.setData(hist,minden,maxden)


class EMImage2DTex:
	def __init__(self,image=None, parent=None):
		
		self.parent = parent
		
		# the OpenGL texture id
		self.tex_id = 0
		# the OpenGL display list number
		self.tex_dl = 0
		
		# this rank is for stencil operations
		self.rank = 1
		
		# we have a brightness contrast screen - it's used in render
		self.bcscreen = BrightContrastScreen()
		
		# need a corner mapper for accurate mouse movement
		self.depthtracker = DepthTracker(self)
		
		# a qt widget inspector
		self.inspector = None
		
		# a camera (encapsulates rotation, movement, and scale)
		self.cam = Camera2(self)
		
		#this disables interactive rotation of the 2D texture
		self.cam.enablerotation = False
		
		#depth tracking should be enabled if you're running python emimage2dtex.py, but should be
		#disabled if you're using this object in emfloatingwidgets
		#similarly - if you disable depthtracking, then you should probably set the camera object
		#variable 'basicmapping' to True
		self.depthtracking = True
		
		# the histogram
		self.hist = None
		self.histbins = 256
		self.histmin = 0.0
		self.histmax = 1.0
		
		if image != None :
			self.setData(image)
	
	
	def isinwin(self,x,y):
		return self.depthtracker.isinwin(x,y)
	
	def set_update_P_inv(self,val=True):
		self.depthtracker.set_update_P_inv(val)
	
	def viewportWidth(self):
		try:
			return self.parent.width()
		except:
			return 0
	
	def viewportHeight(self):
		try:
			return self.parent.height()
		except:
			return 0
	
	def width(self):
		try:
			return self.data.get_xsize()
		except:
			return 0
	
	def height(self):
		try:
			return self.data.get_ysize()
		except:
			return 0
	
	def eyeCoordsDif(self,x1,y1,x2,y2):
		return self.depthtracker.eyeCoordsDif(x1,y1,x2,y2)
	
	def setData(self,data):
		""" Perform the initial data set up"""
		
		if ( not isinstance(data,EMData) ): 
			print "Error, tried to set data using something that was not an EMData"
			return
		
		if ( data.get_zsize() != 1) :
			print "Error, the data is not 2D - the number of z slices is",data.get_zsize()
			return
		
		self.data = data
		
		#normalize the image so that it's on the interval 0,1
		#this will ensure the best starting conditions for
		#texturing
		min = self.data.get_attr("minimum")
		max = self.data.get_attr("maximum")
		self.data.add(-min)
		self.data.mult(1/(max-min))
	
	def setGLContrast(self,val):
		self.bcscreen.setGLContrast(val)
		self.updateHist()
		try:
			self.parent.updateGL()
		except:
			# the parent may not have been set
			pass
	
	def setGLBrightness(self,val):
		self.bcscreen.setGLBrightness(val)
		self.updateHist()
		try:
			self.parent.updateGL()
		except:
			# the parent may not have been set
			pass
	
	def initializeGL(self):
		# This class requires the stencil test
		glClearStencil(0)
		glEnable(GL_STENCIL_TEST)
	
	def initInspector(self):
		if not self.inspector or self.inspector ==None:
			self.inspector=EMImage2DTexInspector(self)

		try:
			self.hist = self.data.calc_hist(self.histbins ,self.histmin,self.histmax)
			self.inspector.setHist(self.hist,self.histmin,self.histmax)
		except:
			print "Error getting histogram from data"
	
	def updateHist(self):
		
		dhist = 1.0/(self.histbins-1.0)
		
		if self.bcscreen.glbrightness == 0: new_zero = 0
		else: new_zero = -self.bcscreen.glbrightness/dhist
		new_zero_val = 0.0
		if new_zero < 0:
			new_zero=0
			new_zero_val = self.bcscreen.glbrightness*self.bcscreen.glcontrast
			
		new_one = ((1.0/self.bcscreen.glcontrast)-self.bcscreen.glbrightness)/dhist
		
		self.inspector.setHist(self.hist[int(new_zero):int(new_one+1)],new_zero_val,1.0) 
	
	def showInspector(self,force):
		try:
			self.inspector.show()
		except:
			if force:
				self.initInspector()
				try:
					self.inspector.show()
				except:
					print "Error - could not force-show the inspector"
			else:
				# no warning
				pass
			
	def genTexture(self):
		if ( self.tex_id != 0 ):
			glDeleteTextures(self.tex_id)
		
		self.tex_id = self.data.gen_gl_texture()
		if (self.tex_id == 0):
			 #OpenGL is not initialized yet?
			print "warning, could not get texture id"
			
	def genCurrentDisplayList(self):
		
		if (self.tex_id == 0): self.genTexture()
		
		if ( self.tex_dl != 0 ): glDeleteLists( self.tex_dl, 1)
		
		self.tex_dl = glGenLists(1)
	
		if (self.tex_dl == 0): 
			 #OpenGL is not initialized yet?
			print "warning, could not get display list"
			return
		
		glNewList(self.tex_dl,GL_COMPILE)
		glEnable(GL_TEXTURE_2D)
		glBindTexture(GL_TEXTURE_2D, self.tex_id)
		#glPixelStorei(GL_UNPACK_ALIGNMENT, 1)
		#glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP)
		#glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP)
		# using GL_NEAREST ensures pixel granularity
		# using GL_LINEAR blurs textures and makes them more difficult
		# to interpret (in cryo-em)
		glTexParameter(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)
		glTexParameter(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)
		# this makes it so that the texture is impervious to lighting
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE)
		
		
		# POSITIONING POLICY - the texture is centered on the origin,
		# and has the height and width of the data
		glBegin(GL_QUADS)
		
		width = self.data.get_xsize()/2.0
		height = self.data.get_ysize()/2.0
		
		glTexCoord2f(0,0)
		glVertex(-width,-height)
		
		glTexCoord2f(1,0)
		glVertex( width,-height)
			
		glTexCoord2f(1,1)
		glVertex( width, height)
		
		glTexCoord2f(0,1)
		glVertex(-width, height)
			
		glEnd()
		
		glDisable(GL_TEXTURE_2D)
		glEndList()
		
	def render(self):
		cull = glIsEnabled(GL_CULL_FACE)
		polygonmode = glGetIntegerv(GL_POLYGON_MODE)
		
		glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
		glDisable(GL_CULL_FACE)

		if ( self.tex_dl == 0 ):
			self.genCurrentDisplayList()
			
		self.cam.position()
		
		
		if ( self.depthtracking == True ): self.depthtracker.storeMatrices()


		glStencilFunc(GL_EQUAL,self.rank,0)
		glStencilOp(GL_KEEP,GL_KEEP,GL_REPLACE)
		glPushMatrix()
		glCallList(self.tex_dl)
		glPopMatrix()
		
		# now do contrast stretching 
		glStencilFunc(GL_EQUAL,self.rank,self.rank)
		glStencilOp(GL_KEEP,GL_KEEP,GL_KEEP)
		glPushMatrix()
		glTranslate(-self.data.get_xsize()/2.0,-self.data.get_ysize()/2.0,.01)
		glScalef(self.data.get_xsize(),self.data.get_ysize(),1)
		self.bcscreen.draw_bc_screen()
		glPopMatrix()
		
		glStencilFunc(GL_ALWAYS,1,1)

		if ( cull ): glEnable(GL_CULL_FACE)
		
		if ( polygonmode[0] == GL_LINE ): glPolygonMode(GL_FRONT, GL_LINE)
		if ( polygonmode[1] == GL_LINE ): glPolygonMode(GL_BACK, GL_LINE)
		
	def wheelEvent(self,event):
		self.cam.wheelEvent(event)

	def mousePressEvent(self, event):
		self.cam.mousePressEvent(event)
		
	def mouseMoveEvent(self,event):
		self.cam.mouseMoveEvent(event)

	def mouseReleaseEvent(self,event):
		self.cam.mouseReleaseEvent(event)
		
class EMImage2DTexGLWidget(QtOpenGL.QGLWidget):
	allim=WeakKeyDictionary()
	def __init__(self, image=None, parent=None):
		
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True)
		fmt.setDepth(1)
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		EMImage2DTexGLWidget.allim[self]=0
		
		self.fov = 50 # field of view angle used by gluPerspective
		
		# a camera or positioning object
		self.cam = Camera2(self)
		
		self.image2dtex = EMImage2DTex(image,self)
	def setData(self,data):
		self.image2dtex.setData(data)
		self.cam.setCamTrans("default_z",-data.get_ysize())
		
	def initializeGL(self):
		GL.glClearColor(0,0,0,0)
		self.image2dtex.initializeGL()
	
	def paintGL(self):
		glClear(GL_COLOR_BUFFER_BIT)
		if glIsEnabled(GL_DEPTH_TEST):
			glClear(GL_DEPTH_BUFFER_BIT)
		if glIsEnabled(GL_STENCIL_TEST):
			glClear(GL_STENCIL_BUFFER_BIT)
		
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		self.cam.position()
		
		glPushMatrix()
		self.image2dtex.render()
		glPopMatrix()
		
	def resizeGL(self, width, height):
		# just use the whole window for rendering
		glViewport(0,0,self.width(),self.height())
		
		# maintain the aspect ratio of the window we have
		self.aspect = float(width)/float(height)
		
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		# using gluPerspective for simplicity
		gluPerspective(self.fov,self.aspect,1,1000)
		self.image2dtex.set_update_P_inv()
		
		# switch back to model view mode
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
	
	def mousePressEvent(self, event):
		if event.button()==Qt.MidButton or (event.button()==Qt.LeftButton and event.modifiers()&Qt.ControlModifier):
			self.image2dtex.showInspector(1)
		else:
			self.image2dtex.mousePressEvent(event)
		
		self.updateGL()
			
	def wheelEvent(self,event):
		self.image2dtex.wheelEvent(event)
		self.updateGL()
	
	def mouseMoveEvent(self,event):
		self.image2dtex.mouseMoveEvent(event)
		self.updateGL()

	def mouseReleaseEvent(self,event):
		self.image2dtex.mouseReleaseEvent(event)
		self.updateGL()

if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = EMImage2DTexGLWidget()
 	if len(sys.argv)==1 : 
		e = EMData()
		e.set_size(512,512,1)
		e.process_inplace('testimage.scurve')
		window.setData(e)
	else :
		if not os.path.exists(sys.argv[1]):
			print "Error, input file %s does not exist" %sys.argv[1]
			exit(1)
		a=EMData.read_images(sys.argv[1],[0])
		window.setData(a[0])
	window2=EMParentWin(window)
	window2.show()
	
	sys.exit(app.exec_())

		