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

# EMFXImage.py  Steve Ludtke  08/06/2006
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
from emimagemx import *
from math import sqrt

from emimage import EMImage

from emimage3dobject import Camera2
from emimage2dtex import *

import numpy
import sys
import array

height_plane = 500

class EMBasicOpenGLObjects:
	"""
	This class is supposed to encapsulate basic and common objects
	used by various other OpenGL classes in EMAN2 pyqt interfaces.
	
	It currently supplies display list ids - one for a sphere, one for
	a cylinder
	
	It's implemented as a singleton
	"""
	class __impl:
		""" Implementation of the singleton interface """

		def __init__(self):
			# this display list ids
			self.cylinderdl = 0
			self.spheredl = 0
		
			# cylinder parameters
			self.cylinder_around_z = 12
			self.cylinder_along_z = 2
			
			# sphere parameters
			self.sphere_around_z = 6
			self.sphere_along_z = 6
			
			self.gq=gluNewQuadric()
			gluQuadricDrawStyle(self.gq,GLU_FILL)
			gluQuadricNormals(self.gq,GLU_SMOOTH)
			gluQuadricOrientation(self.gq,GLU_OUTSIDE)
			gluQuadricTexture(self.gq,GL_FALSE)
		
		def __del__(self):
			if self.cylinderdl != 0:
				glDeleteLists(self.cylinderdl,1)
			if self.spheredl != 0:
				glDeleteLists(self.spheredl,1)
			
		def getSphereDL(self):
			if ( self.spheredl == 0 ):
				self.spheredl=glGenLists(1)
				
				glNewList(self.spheredl,GL_COMPILE)
				glPushMatrix()
				gluSphere(self.gq,.5,self.sphere_along_z,self.sphere_around_z)
				glPopMatrix()
				
				glEndList()
				
			return self.spheredl
		
		def getCylinderDL(self):
			if ( self.cylinderdl == 0 ):
				self.cylinderdl=glGenLists(1)
				
				glNewList(self.cylinderdl,GL_COMPILE)
				glPushMatrix()
				gluCylinder(self.gq,1.0,1.0,1.0,self.cylinder_around_z,self.cylinder_along_z)
				glPopMatrix()
				
				glEndList()
				
			return self.cylinderdl

	# storage for the instance reference
	__instance = None


	def __init__(self):
		""" Create singleton instance """
		# Check whether we already have an instance
		if EMBasicOpenGLObjects.__instance is None:
			# Create and remember instance
			EMBasicOpenGLObjects.__instance = EMBasicOpenGLObjects.__impl()
	
	def __getattr__(self, attr):
		""" Delegate access to implementation """
		return getattr(self.__instance, attr)

	def __setattr__(self, attr, value):
		""" Delegate access to implementation """
		return setattr(self.__instance, attr, value)

class EMGLProjViewMatrices:
	"""
	Controls a static single instance of OpenGL projection and view(port) matrices
	Also stores a static single instance of the inverse of the projection matrix, which
	is calculates using numpy/scipy. Other classes interface with this by calling
	setUpdate when the viewport or projection matrix has changed (i.e., a resize event)
	"""
	class __impl:
		""" Implementation of the singleton interface """

		def __init__(self):
			self.update = True
			self.proj_matrix = None
			self.inv_proj_matrix = None
			self.view_matrix = None
		
		def updateOpenGLMatrices(self, val=True):
			self.updateViewMatrix()
			self.updateProjectionMatrix()
			self.updateProjectionInvMatrix()
			
		def setUpdate(self,val=True):
			self.update = val
		
		def updateProjectionInvMatrix(self):
			P = numpy.matrix(self.proj_matrix)
			self.inv_proj_matrix = P.I
			
		def updateProjectionMatrix(self):
			self.proj_matrix = glGetDoublev(GL_PROJECTION_MATRIX)
		
		def updateViewMatrix(self):
			 self.view_matrix = glGetIntegerv(GL_VIEWPORT)
		
		def checkUpdate(self):
			if (self.update):
				self.updateOpenGLMatrices()
				self.update = False
		
		def getViewMatrix(self):
			self.checkUpdate()
			return self.view_matrix
		
		def getProjMatrix(self):
			self.checkUpdate()
			return self.proj_matrix
		
		def getProjInvMatrix(self):
			self.checkUpdate()
			return self.inv_proj_matrix
			
		
			

	# storage for the instance reference
	__instance = None


	def __init__(self):
		""" Create singleton instance """
		# Check whether we already have an instance
		if EMGLProjViewMatrices.__instance is None:
			# Create and remember instance
			EMGLProjViewMatrices.__instance = EMGLProjViewMatrices.__impl()
	
	def __getattr__(self, attr):
		""" Delegate access to implementation """
		return getattr(self.__instance, attr)

	def __setattr__(self, attr, value):
		""" Delegate access to implementation """
		return setattr(self.__instance, attr, value)

class EMViewportDepthTools:
	"""
	This class provides important tools for EMAN2 floating widgets -
	these are either qt widgets that get mapped as textures to
	polygons situated in OpenGL volumes, or are OpenGL objects themselves
	(such as a 2D image texture, or a 3D isosurface), that get drawn
	within a frame somewhere in the 3D world.
	
	These objects need to have mouse events (etc) correctly rerouted to them,
	and this is not trivial (but not difficult) to do, considering that the
	mouse events are always in terms of the viewport, but the texture mapped
	widget is somewhere in 3D space. The function eyeCoordsDif is primarily
	for positioning the widgits in 3D space (translation), whereas the 
	mouseinwin function is specifically for mapping the mouse event coordinates
	into the widgit's transformed coordinate system.

	This class also provides important collision detection functionality - it
	does this by mapping the corners of the (widgit mapped) polygon to the viewport,
	and then determing if the (mouse) coordinate is within this area. Mapping 
	polygon vertices to the veiwport is done using gluUnproject, whereas converting
	viewport coordinates into polygon coordinates is done doing something similar to
	gluProject (the opposite operation).

	This class also provides widget frame drawing functionality (which may become redundant)
	It does this by drawing cylinders in terms of how the polygon corners mapped
	to the viewport - this causes the loss of depth information
	
	The only important behaviour expected of something that uses this class is
	1 - you must call update() just before you draw the textured widgit polygon
	other object (i.e. when the contents of the OpenGL modelview matrix reflect 
	all of the operations that are applied before rendering)
	2 - you should call set_update_P_inv() if the OpenGL projection matrix is altered,
	this typically happens when resizeGL is called in the root widgit.
	"""
	def __init__(self, parent):
		self.parent = parent
		
		self.glbasicobjects = EMBasicOpenGLObjects()		# need basic objects for drawing the frame
		self.matrices = EMGLProjViewMatrices()				# an object that stores static instances of the important viewing matrices
		self.borderwidth = 3.0								# effects the border width of the frame decoration
		
	def set_update_P_inv(self,val=True):
		self.matrices.setUpdate(val)
	
	def drawFrame(self, ftest=False):
		
		if (ftest):
			a = Vec3f(self.mc10[0]-self.mc00[0],self.mc10[1]-self.mc00[1],self.mc10[2]-self.mc00[2])
			b = Vec3f(self.mc01[0]-self.mc00[0],self.mc01[1]-self.mc00[1],self.mc01[2]-self.mc00[2])
			c = a.cross(b)
			if ( c[2] < 0 ):
				#print "facing backward"
				return False
		glMatrixMode(GL_PROJECTION)
		glPushMatrix()
		glLoadIdentity()
		glOrtho(0.0,self.parent.viewportWidth(),0.0,self.parent.viewportHeight(),-5,5)

		glMatrixMode(GL_MODELVIEW)
		glPushMatrix()
		glLoadIdentity()
		
		glColor(.9,.2,.8)
		# this is a nice light blue color (when lighting is on)
		# and is the default color of the frame
		glMaterial(GL_FRONT,GL_AMBIENT,(.2,.2,.8,1.0))
		glMaterial(GL_FRONT,GL_SPECULAR,(.8,.8,.8,1.0))
		glMaterial(GL_FRONT,GL_SHININESS,50.0)
		
		#draw the cylinders around the edges of the frame
		glPushMatrix()
		self.cylinderToFrom(self.mc00,self.mc10)
		glPopMatrix()
		glPushMatrix()
		self.cylinderToFrom(self.mc10,self.mc11)
		glPopMatrix()
		glPushMatrix()
		self.cylinderToFrom(self.mc11,self.mc01)
		glPopMatrix()
		glPushMatrix()
		self.cylinderToFrom(self.mc01,self.mc00)
		glPopMatrix()
		
		# draw spheres as nodes at the corners of the frame
		glPushMatrix()
		self.sphereAt(self.mc00)
		glPopMatrix()
		glPushMatrix()
		self.sphereAt(self.mc10)
		glPopMatrix()
		glPushMatrix()
		self.sphereAt(self.mc11)
		glPopMatrix()
		glPushMatrix()
		self.sphereAt(self.mc01)
		glPopMatrix()
		
		glPopMatrix()
		
		glMatrixMode(GL_PROJECTION)
		# pop the temporary orthographic matrix from the GL_PROJECTION stack
		glPopMatrix()
		glMatrixMode(GL_MODELVIEW)
		
		return True
	def cylinderToFrom(self,To,From):
		#print "To",To[0],To[1]
		#print "From",From[0],From[1]
		# draw a cylinder to To from From
		dx = To[0] - From[0]
		dy = To[1] - From[1]
		
		length = sqrt(dx*dx+dy*dy)
		
		angle = 180.0*atan2(dy,dx)/pi
		
		glTranslated(From[0],From[1],0)
		glRotated(90.0+angle,0.,0.,1.0)
		glRotated(90.,1.,0.,0.)
		glTranslated(0.0,0.0,length/2.0)
		glScaled(self.borderwidth,self.borderwidth,length)
		glTranslated(0.0,0.0,-0.5)
		glCallList(self.glbasicobjects.getCylinderDL())
		
	def sphereAt(self,at):
		glTranslate(at[0],at[1],0)
		glScalef(3.0*self.borderwidth,3.0*self.borderwidth,3.0*self.borderwidth)
		glCallList(self.glbasicobjects.getSphereDL())
	
	def update(self,width,height):
		
		self.wmodel= glGetDoublev(GL_MODELVIEW_MATRIX)
		self.wproj = self.matrices.getProjMatrix()
		self.wview = self.matrices.getViewMatrix()

		try:
			self.mc00=gluProject(-width,-height,0.,self.wmodel,self.wproj,self.wview)
			self.mc10=gluProject( width,-height,0.,self.wmodel,self.wproj,self.wview)
			self.mc11=gluProject( width, height,0.,self.wmodel,self.wproj,self.wview)
			self.mc01=gluProject(-width, height,0.,self.wmodel,self.wproj,self.wview)
			
		except:
			self.mc00 = [0,0,0]
			self.mc10 = [0,0,0]
			self.mc11 = [0,0,0]
			self.mc01 = [0,0,0]
			
	def getCorners(self):
		return [self.mc00,self.mc01,self.mc11,self.mc10]
	
	def getModelMatrix(self):
		return self.wmodel
	
	def setModelMatrix(self, model):
		self.wmodel = model
	
	def setCorners(self,points):
		self.mc00 = points[0]
		self.mc01 = points[1]
		self.mc11 = points[2]
		self.mc10 = points[3]
		
	def storeModel(self):
		self.wmodel= glGetDoublev(GL_MODELVIEW_MATRIX)
	
	def printUnproj(self,x,y):
		self.wmodel= glGetDoublev(GL_MODELVIEW_MATRIX)
		self.wproj = self.matrices.getProjMatrix()
		self.wview = self.matrices.getViewMatrix()
	
		t = gluProject(x,y,0.,self.wmodel,self.wproj,self.wview)
		print t[0],t[1],t[2]
		
		
	def getMappedWidth(self):
		return int(self.mc11[0]-self.mc00[0])
		
	def getMappedHeight(self):
		return int(self.mc11[1]-self.mc00[1])
	
	def isinwinpoints(self,x,y,points):
		try:
			a = [points[0][0]-x, points[0][1]-y]
			b = [points[1][0]-x, points[1][1]-y]
			c = [points[2][0]-x, points[2][1]-y]
			d = [points[3][0]-x, points[3][1]-y]
			
			aeb = self.getsubtendingangle(a,b)
			bec = self.getsubtendingangle(b,c)
			ced = self.getsubtendingangle(c,d)
			dea = self.getsubtendingangle(d,a)
			if abs(aeb + bec + ced + dea) > 0.1:
				return True 
			else:
				return False
		except:
			return False
	
	def isinwin(self,x,y):
		# this function can be called to determine
		# if the event at x,y (in terms of the viewport)
		# was within the frame of the object being drawn
		try:
			a = [self.mc00[0]-x, self.mc00[1]-y]
			b = [self.mc01[0]-x, self.mc01[1]-y]
			c = [self.mc11[0]-x, self.mc11[1]-y]
			d = [self.mc10[0]-x, self.mc10[1]-y]
			
			aeb = self.getsubtendingangle(a,b)
			bec = self.getsubtendingangle(b,c)
			ced = self.getsubtendingangle(c,d)
			dea = self.getsubtendingangle(d,a)
			if abs(aeb + bec + ced + dea) > 0.1:
				return True 
			else:
				return False
		except:
			return False
	
	def getsubtendingangle(self,a,b):
		sinaeb = a[0]*b[1]-a[1]*b[0]
		cosaeb = a[0]*b[0]+a[1]*b[1]
		
		return atan2(sinaeb,cosaeb)
	
	def eyeCoordsDif(self,x1,y1,x2,y2,maintaindepth=True):
		# get x and y normalized device coordinates
		xNDC1 = 2.0*(x1-self.wview[0])/self.wview[2] - 1
		yNDC1 = 2.0*(y1-self.wview[1])/self.wview[3] - 1
		
		xNDC2 = 2.0*(x2-self.wview[0])/self.wview[2] - 1
		yNDC2 = 2.0*(y2-self.wview[1])/self.wview[3] - 1
		

		self.P_inv = self.matrices.getProjInvMatrix()
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

	def mouseinwin(self,x,y,width,height):
		# to determine the mouse coordinates in the window we carefully perform
		# linear algebra similar to what's done in gluUnProject

		# the problem is determining what the z coordinate of the mouse event should have
		# been, given that we know that the widget itself is located in the x,y plane, along z=0.
		
		# get x and y normalized device coordinates
		xNDC = 2.0*(x-self.wview[0])/self.wview[2] - 1
		yNDC = 2.0*(y-self.wview[1])/self.wview[3] - 1
		
		# invert the projection and model view matrices, they will be used shortly
		# note the OpenGL returns matrices are in column major format -  the calculations below 
		# are done with this in  mind - this saves the need to transpose the matrices
		self.P_inv = self.matrices.getProjInvMatrix()
		M = numpy.matrix(self.wmodel)
		
		M_inv = M.I
		#PM_inv = numpy.matrixmultiply(P_inv,M_inv)
		PM_inv = self.P_inv*M_inv
		
		# If the widget is planar (which obviosuly holds), and along z=0, then the following holds
		zNDC = (PM_inv[0,2]*xNDC + PM_inv[1,2]*yNDC + PM_inv[3,2])/(-PM_inv[2,2])
	
		# We need zprime, which is really 'eye_z' in OpenGL lingo
		zprime = 1.0/(xNDC*self.P_inv[0,3]+yNDC*self.P_inv[1,3]+zNDC*self.P_inv[2,3]+self.P_inv[3,3])
		
		# Now we compute the x and y coordinates - these are precisely what we're after
		xcoord = zprime*(xNDC*PM_inv[0,0]+yNDC*PM_inv[1,0]+zNDC*PM_inv[2,0]+PM_inv[3,0])
		ycoord = zprime*(xNDC*PM_inv[0,1]+yNDC*PM_inv[1,1]+zNDC*PM_inv[2,1]+PM_inv[3,1])

		return (xcoord + width*0.5, 0.5*height-ycoord)

class BoundaryMediator:
	def __init__(self,drawer,depthtools):
		self.drawer = drawer
		self.depthtools = depthtools
		
	def checkBoundaryIssues(self):
		corners = self.depthtools.getCorners()
		for i in corners:
			if i[0] < 0:
				print "out to the left"
			elif i[0] > self.drawer.viewportWidth():
				print "out to the right"
				
			if i[1] < 0:
				print "out below"
			elif i[1] > self.drawer.viewportHeight():
				print "out above"

class EMGLView3D:
	"""
	A view of an EMAN2 3D type, such as an isosurface or a 
	volume rendition, etc.
	"""
	def __init__(self, parent,image=None):
		self.parent = parent
		self.cam = Camera2(self)
		self.cam.setCamTrans('default_z',-parent.get_depth_for_height(height_plane))
		
		self.w = image.get_xsize()	# width of window
		self.h = image.get_ysize()	# height of window
		self.d = image.get_zsize()	# depth of the window
		self.sizescale = 1.0		# scale/zoom factor
		self.changefactor = 1.1		# used to zoom
		self.invchangefactor = 1.0/self.changefactor # used to invert zoom
		
		self.drawable = None		# the object that is drawable (has a draw function)
		
		self.vdtools = EMViewportDepthTools(self)
		
		self.updateFlag = True
		
		self.drawFrame = True
		
		self.vdtools
		
		self.psets = []
		self.modelmatrices = []
		
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
	
	def depth(self):
		try:
			return int(self.sizescale*self.d)
		except:
			return 0
	
	def setData(self,data):
		try: self.drawable.setData(data)
		except: pass
		
		
	def paintGL(self):
		self.psets = []
		self.modelmatrices = []
		self.cam.position()
		lighting = glIsEnabled(GL_LIGHTING)
		glEnable(GL_LIGHTING)
		
		glPushMatrix()
		glTranslatef(-self.width()/2.0,0,0)
		glRotatef(-90,0,1,0)
		self.vdtools.update(self.depth()/2.0,self.height()/2.0)
		if self.drawFrame: 
			if self.vdtools.drawFrame(True): 
				self.psets.append(self.vdtools.getCorners())
				self.modelmatrices.append(self.vdtools.getModelMatrix())
		glPopMatrix()
		
		glPushMatrix()
		glTranslatef(self.width()/2.0,0,0)
		glRotatef(90,0,1,0)
		self.vdtools.update(self.depth()/2.0,self.height()/2.0)
		if self.drawFrame: 
			if self.vdtools.drawFrame(True): 
				self.psets.append(self.vdtools.getCorners())
				self.modelmatrices.append(self.vdtools.getModelMatrix())
				
		glPopMatrix()
		
		glPushMatrix()
		glTranslatef(0,self.height()/2.0,0)
		glRotatef(-90,1,0,0)
		self.vdtools.update(self.width()/2.0,self.depth()/2.0)
		if self.drawFrame: 
			if self.vdtools.drawFrame(True): 
				self.psets.append(self.vdtools.getCorners())
				self.modelmatrices.append(self.vdtools.getModelMatrix())
		glPopMatrix()
		
		glPushMatrix()
		glTranslatef(0,-self.height()/2.0,0)
		glRotatef(90,1,0,0)
		self.vdtools.update(self.width()/2.0,self.depth()/2.0)
		if self.drawFrame: 
			if self.vdtools.drawFrame(True): 
				self.psets.append(self.vdtools.getCorners())
				self.modelmatrices.append(self.vdtools.getModelMatrix())
		glPopMatrix()
		
		glPushMatrix()
		glTranslatef(0,0,-self.depth()/2.0)
		glRotatef(180,0,1,0)
		self.vdtools.update(self.depth()/2.0,self.height()/2.0)
		if self.drawFrame: 
			if self.vdtools.drawFrame(True): 
				self.psets.append(self.vdtools.getCorners())
				self.modelmatrices.append(self.vdtools.getModelMatrix())
		glPopMatrix()
		
		glPushMatrix()
		glTranslatef(0,0,self.depth()/2.0)
		self.vdtools.update(self.width()/2.0,self.height()/2.0)
		if self.drawFrame: 
			if self.vdtools.drawFrame(True): 
				self.psets.append(self.vdtools.getCorners())
				self.modelmatrices.append(self.vdtools.getModelMatrix())
		glPopMatrix()
		
		if not lighting: glDisable(GL_LIGHTING)
		
	def viewportHeight(self):
		return self.parent.height()
	
	def viewportWidth(self):
		return self.parent.width()
	
	def mousePressEvent(self, event):
		if event.button()==Qt.MidButton or (event.button()==Qt.LeftButton and event.modifiers()&Qt.ControlModifier):	
			return
			self.drawable.initInspector()
			self.drawable.inspector.show()
			self.drawable.inspector.hide()
			self.parent.addQtWidgetDrawer(self.drawable.inspector)
			
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.mousePressEvent(event)
		else:
			return
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

		#self.drawable.resizeEvent(self.width(),self.height())
	
	def wheelEvent(self,event):
		if event.modifiers() == Qt.ShiftModifier:
			self.scaleEvent(event.delta())
			#print "updating",self.drawWidth(),self.drawHeight()
			self.updateGL()
		else:
			return
			self.drawable.wheelEvent(event)
			
		#self.updateGL()
	
	def mouseMoveEvent(self,event):
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.mouseMoveEvent(event)
		else:
			return
			#l=self.vdtools.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
			#qme=QtGui.QMouseEvent(event.type(),QtCore.QPoint(l[0],l[1]),event.button(),event.buttons(),event.modifiers())
			#self.drawable.mouseMoveEvent(qme)
			#self.drawable.mouseMoveEvent(event)
		
		#self.updateGL()

	def mouseReleaseEvent(self,event):
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.mouseReleaseEvent(event)
		else:
			return
			l=self.vdtools.mouseinwin(event.x(),self.parent.height()-event.y(),self.width(),self.height())
			qme=QtGui.QMouseEvent(event.type(),QtCore.QPoint(l[0],l[1]),event.button(),event.buttons(),event.modifiers())
			self.drawable.mouseReleaseEvent(qme)
			#self.drawable.mouseReleaseEvent(event)
	
	def update(self):
		self.parent.updateGL()
	
	def updateGL(self):
		self.parent.updateGL()
	
	def isinwin(self,x,y):
		val = False
		for i,p in enumerate(self.psets):
			if self.vdtools.isinwinpoints(x,y,p):
				val = True
				self.vdtools.setModelMatrix(self.modelmatrices[i])
				break
		return val
	
	def eyeCoordsDif(self,x1,y1,x2,y2,mdepth=True):
		return self.vdtools.eyeCoordsDif(x1,y1,x2,y2,mdepth)
	
	def leaveEvent(self):
		return
		self.drawable.leaveEvent()
	
	def toolTipEvent(self,event):
		pass
	
	def set_update_P_inv(self,val=True):
		self.vdtools.set_update_P_inv(val)
	
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
		
		self.drawable.supressInspector = True
		self.initflag = True
		self.vdtools = EMViewportDepthTools(self)
		
		self.updateFlag = True
		
		self.drawFrame = True
		
		self.sizescale = 1.0
		self.changefactor = 1.1
		self.invchangefactor = 1.0/self.changefactor
		
		self.mediator = BoundaryMediator(self,self.vdtools)
		
	def become2DImage(self,a):
		self.drawable = EMImage2DCore(a,self)
		#self.drawable.originshift = False
		self.w = a.get_xsize()
		self.h = a.get_ysize()
		
	def eyeCoordsDif(self,x1,y1,x2,y2,mdepth=True):
		return self.vdtools.eyeCoordsDif(x1,y1,x2,y2,mdepth)
	
	def set_update_P_inv(self,val=True):
		self.vdtools.set_update_P_inv(val)
	
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
	
	def mousePressEvent(self, event):
		if event.button()==Qt.MidButton or (event.button()==Qt.LeftButton and event.modifiers()&Qt.ControlModifier):	
			self.drawable.initInspector()
			self.drawable.inspector.show()
			self.drawable.inspector.hide()
			self.parent.addQtWidgetDrawer(self.drawable.inspector)
			
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
			self.updateGL()
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
		self.parent.emit(signal,event)
	
	def leaveEvent(self):
		self.drawable.leaveEvent()
	
	def toolTipEvent(self,event):
		pass
	
	def isinwin(self,x,y):
		return self.vdtools.isinwin(x,y)

class EMQtWidgetDrawer:
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
			pixmap = QtGui.QPixmap.grabWidget(self.qwidget)
			if (pixmap.isNull() == True ): print 'error, the pixmap was null'
			self.itex = self.parent.bindTexture(pixmap)
			if ( self.itex == 0 ): print 'Error - I could not generate the texture'
		
	def paintGL(self):
		#print "paintGL children"
		if (self.qwidget == None or self.itex == 0) :
			#print "no widget - paintGL children return" 
			return
		
		self.cam.position()
		
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
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.wheelEvent(event)
		else:
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
				widget = EMQtWidgetDrawer(self.parent,None,cw);
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


class EMFXImage(QtOpenGL.QGLWidget):
	"""A QT widget for rendering EMData objects. It can display single 2D or 3D images 
	or sets of 2D images.
	"""
	def __init__(self, parent=None):
		#print "init"
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True)
		# enable multisampling to combat aliasing
		fmt.setSampleBuffers(True)
		# stenciling is for object dependent shading
		fmt.setStencil(True)
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		
		self.fov = 2*180*atan2(1,5)/pi
		self.imtex=0
		
		self.setMouseTracking(True)
		self.current = None
		self.previous = None
	
		self.initFlag = True
		self.qwidgets = []
		
		#print "init done"
	
		
	def get_depth_for_height(self, height):
		# This function returns the width and height of the renderable 
		# area at the origin of the data volume
		depth = height/(2.0*tan(self.fov/2.0*pi/180.0))
		return depth
		
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
	def addQtWidgetDrawer(self,widget):
		w = EMQtWidgetDrawer(self)
		w.setQtWidget(widget)
		self.qwidgets.append(w)
		
		#print "initializeGL done"
	def paintGL(self):
		if ( self.initFlag == True ):
			self.fd = QtGui.QFileDialog(self,"Open File",QtCore.QDir.currentPath(),QtCore.QString("Image files (*.img *.hed *.mrc)"))
			QtCore.QObject.connect(self.fd, QtCore.SIGNAL("finished(int)"), self.finished)
			self.fd.show()
			self.fd.hide()
			self.qwidgets.append(EMQtWidgetDrawer(self))
			self.qwidgets[0].setQtWidget(self.fd)
			self.qwidgets[0].cam.setCamX(-100)
			self.initFlag = False
			
			e = EMData(256,256,256)
			e.process_inplace("testimage.axes")
			w = EMGLView3D(self,e)
			self.qwidgets.append(w)
		#print "paintGL"
		glClear(GL_COLOR_BUFFER_BIT)
		if glIsEnabled(GL_DEPTH_TEST):
			glClear(GL_DEPTH_BUFFER_BIT)
		if glIsEnabled(GL_STENCIL_TEST):
			glClear(GL_STENCIL_BUFFER_BIT)
			
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
	
		for i in self.qwidgets:
			#print "getting opengl matrices"
			glPushMatrix()
			i.paintGL()
			glPopMatrix()
			#print "paint child done"
		
		#print "paintGL done"
	def finished(self,val):
		#print "file dialog finished with code",val
		if ( val == 1 ):
			for i in self.fd.selectedFiles():
				#try:
					#print i,str(i)
					a=EMData.read_images(str(i))
					#if len(a)==1 : a=a[0]
					#w.setWindowTitle("EMImage (%s)"%f)
					#w.show()
					w = EMGLView2D(self,a)
					#w.setData(a[0])
					self.qwidgets.append(w)
					
					#self.qwidgets[0].cam.setCamX(100)
					#self.initFlag = False
				#except:
					#print "error, could not open",i
			
	def timer(self):
		pass
		#self.updateGL()
	
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
		gluPerspective(self.fov,self.aspect,depth-depth/4,depth+depth/4)
		for i in self.qwidgets:
			i.set_update_P_inv()
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		#self.updateGL()
		#print "resizeGL done"
	
	
	def mousePressEvent(self, event):
		for i in self.qwidgets:
			if ( i.isinwin(event.x(),self.height()-event.y()) ):
				i.mousePressEvent(event)
				intercepted = True
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
	window = EMFXImage()
	window2 = EMParentWin(window)
	window2.show()
	
	sys.exit(app.exec_())
