#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
# Author: David Woolford Early 2008 woolford@bcm.edu
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
from math import sqrt

from emimage3dobject import Camera2

import numpy
import sys
import array

height_plane = 500

class EMBasicObjects:
	def __init__(self):
		self.framedl = 0
		self.cylinderdl = 0
		self.spheredl = 0
	
		self.cylinder_around_z = 12
		self.cylinder_along_z = 2
		
		self.sphere_around_z = 6
		self.sphere_along_z = 6
		
		self.gq=gluNewQuadric()
		gluQuadricDrawStyle(self.gq,GLU_FILL)
		gluQuadricNormals(self.gq,GLU_SMOOTH)
		gluQuadricOrientation(self.gq,GLU_OUTSIDE)
		gluQuadricTexture(self.gq,GL_FALSE)
		
		# call setter methods to changes these
		self.width = 1.0
		self.height = 1.0
	
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
	
	def getFrameDL(self):
		# draws a frame using Quadrics
		if ( self.framedl == 0 ):
			self.framedl=glGenLists(1)
			
			glNewList(self.framedl,GL_COMPILE)
			glPushMatrix()
			
			glRotate(90.,1.,0.,0.)
			glTranslate(1.0,0.,-1.0)
			gluCylinder(self.gq,.02,.02,2.04,12,2)
			glTranslate(-2.04,0.,0.)
			gluCylinder(self.gq,.02,.02,2.04,12,2)
			
			glPopMatrix()
			
			glPushMatrix()
				
			glRotate(90.,0.,1.,0.)
			glTranslate(0.,1.02,-1.02)
			gluCylinder(self.gq,.02,.02,2.04,12,2)
			glTranslate(0.,-2.04,0.)
			gluCylinder(self.gq,.02,.02,2.04,12,2)
			
			glPopMatrix()
			
			glPushMatrix()
			
			glTranslate(1.02,1.02,0.)
			gluSphere(self.gq,.02,6,6)
			glTranslate(-2.04,0.,0.)
			gluSphere(self.gq,.02,6,6)
			glTranslate(0.,-2.04,0.)
			gluSphere(self.gq,.02,6,6)
			glTranslate(2.04,0.,0.)
			gluSphere(self.gq,.02,6,6)
			
			glPopMatrix()
			
			glEndList()
			
		return self.framedl
		
	def setWidth(self,width):
		self.width = width
	def setHeight(self,height):
		self.height = height
	

class EMQtWidgetDrawer:
	def __init__(self, parent=None, qwidget=None, widget_parent=None):
		self.parent = parent
		self.qwidget = qwidget
		self.drawframe = True
		self.mapcoords = True
		self.itex = 0
		self.genTexture = True
		self.click_debug = False
		self.cam = Camera2(self)
		self.cam.setCamTrans('default_z',-1250)
		self.cam.motionRotate(0,0)
		self.borderwidth = 3.0
		self.glbasicobjects = EMBasicObjects()
		self.setQtWidget(qwidget)
		self.P_inv = None
		self.update_P_inv = True
		self.childreceiver = None
		self.widget_parent = widget_parent
		
		self.current = None
		self.previous = None
		
		self.e2children = []
		self.is_child = False
	
	def set_update_P_inv(self,val=True):
		self.update_P_inv = val
		for i in self.e2children:
			i.set_update_P_inv(val)
	
	def width(self):
		return self.qwidget.width()
	
	def height(self):
		return self.qwidget.height()
	
	def parentHeight(self):
		return self.parent.height()
	def parentWidth(self):
		return self.parent.width()
	
	def setQtWidget(self, widget, delete_current = False):
		if ( delete_current and self.qwidget != None ):
			self.qwidget.deleteLater()
		
		self.qwidget = widget
		
		if ( widget != None ):
			self.glbasicobjects.setWidth(self.qwidget.width())
			self.glbasicobjects.setHeight(self.qwidget.height())
			self.genTexture = True
			self.updateTexture()
			
	def updateTexture(self):
		if ( self.itex == 0 or self.genTexture == True ) : 
			if (self.itex != 0 ):
				#passpyth
				self.parent.deleteTexture(self.itex)
			self.genTexture = False
			#print "binding texture"
			self.itex = self.parent.bindTexture(QtGui.QPixmap.grabWidget(self.qwidget))
			if ( self.itex == 0 ): print 'NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO'
			#print "I have", len(self.qwidget.children())
			#for i in self.qwidget.children():
				#if len(i.children()) > 0:
					#for j in i.children():
						#print j
			#print self.itex
			#print 'done'
		
	def print_angles(self):
		self.print_angles_d(self.mc00,self.mc10)
		self.print_angles_d(self.mc10,self.mc11)
		self.print_angles_d(self.mc11,self.mc01)
		self.print_angles_d(self.mc01,self.mc00)
	
	def print_angles_d(self,p1,p2,u=1,v=0):
		x = p2[0]-p1[0]
		y = p2[1]-p1[1]
		
		l1 = sqrt(x*x+y*y)
		l2 = sqrt(u*u+v*v)
		
		theta = acos((x*u+v*y)/(l1*l2))
		print "the angle between [%f,%f] and [%f,%f] is %f" %(p1[0],p1[1],p2[0],p2[1],180.0*theta/pi)
		
	def paintGL(self):
		#print "in paint", self.qwidget
		if (self.qwidget == None ) : return
		
		
		self.cam.position()
		
		self.wmodel=glGetDoublev(GL_MODELVIEW_MATRIX)
		self.wproj=glGetDoublev(GL_PROJECTION_MATRIX)
		self.wview=glGetIntegerv(GL_VIEWPORT)
		
		glPushMatrix()
		glEnable(GL_TEXTURE_2D)
		glBindTexture(GL_TEXTURE_2D,self.itex)
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
		#self.parent.deleteTexture(self.itex)
	
		if ( self.mapcoords ):
			self.mc00=gluProject(-self.qwidget.width()/2.0,-self.qwidget.height()/2.0,0.,self.wmodel,self.wproj,self.wview)
			self.mc10=gluProject( self.qwidget.width()/2.0,-self.qwidget.height()/2.0,0.,self.wmodel,self.wproj,self.wview)
			self.mc11=gluProject( self.qwidget.width()/2.0, self.qwidget.height()/2.0,0.,self.wmodel,self.wproj,self.wview)
			self.mc01=gluProject(-self.qwidget.width()/2.0, self.qwidget.height()/2.0,0.,self.wmodel,self.wproj,self.wview)
	
			if ( self.drawframe ):
				glMatrixMode(GL_PROJECTION)
				glPushMatrix()
				glLoadIdentity()
				glOrtho(0.0,self.parent.width(),0.0,self.parent.height(),-5,5)
	
				glMatrixMode(GL_MODELVIEW)
				glPushMatrix()
				glLoadIdentity()
				
				glColor(.9,.2,.8)
				glMaterial(GL_FRONT,GL_AMBIENT,(.2,.2,.8,1.0))
				glMaterial(GL_FRONT,GL_SPECULAR,(.8,.8,.8,1.0))
				glMaterial(GL_FRONT,GL_SHININESS,50.0)
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
				
				#glBegin(GL_LINES)
				#glVertex(self.mc00[0],self.mc00[1])
				#glVertex(self.mc11[0],self.mc11[1])
				#glVertex(self.mc01[0],self.mc01[1])
				#glVertex(self.mc10[0],self.mc10[1])
				#glEnd()
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
				
		#if ( self.drawframe ):
			#glPushMatrix()
			#glCallList(self.glbasicobjects.getFrameDL())
			#glPopMatrix()
			
		for i in self.e2children:
			glPushMatrix()
			i.paintGL()
			glPopMatrix()
			
		#print "leave paint"
	
	def sphereAt(self,at):
		glTranslate(at[0],at[1],0)
		glScalef(3.0*self.borderwidth,3.0*self.borderwidth,3.0*self.borderwidth)
		glCallList(self.glbasicobjects.getSphereDL())
		
	def cylinderToFrom(self,To,From):
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

	def isinwin(self,x,y):
		for i in self.e2children:
			if i.isinwin(x,y):
				self.childreceiver = i
				return True
		# this works by simple geometry - if the mouse point e is within the four points (a,b,c,d)
		# at the extremities of  the qtwidget, then aed + bec + ced + dea is +/- 360 degrees. 
		# If e is outside the four points then the sum is zero...
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
		
	def check(self):
		[A,B,C,D] = self.equation_of_plane(self.mc00,self.mc01,self.mc11)
		[A,B,C,D] = self.equation_of_plane(self.mc11,self.mc01,self.mc00)
		[A,B,C,D] = self.equation_of_plane(self.mc10,self.mc01,self.mc11)
		[A,B,C,D] = self.equation_of_plane(self.mc10,self.mc01,self.mc00)
		
	def equation_of_plane(self,a,b,c):
		x1,y1,z1=a[0],a[1],a[2]
		x2,y2,z2=b[0],b[1],b[2]
		x3,y3,z3=c[0],c[1],c[2]
		A = y1*(z2-z3)+y2*(z3-z1)+y3*(z1-z2)
		B = z1*(x2-x3)+z2*(x3-x1)+z3*(x1-x2)
		C = x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2)
		D = -(x1*(y2*z3-y3*z2)+x2*(y3*z1-y1*z3)+x3*(y1*z2-y2*z1))
		
		print A,B,C,D
		return [A,B,C,D]
	
	def mouseViewportMovement(self,x,y,ex,ey,ez,zprime):
		xNDC = 2.0*(x-self.wview[0])/self.wview[2] - 1
		yNDC = 2.0*(y-self.wview[1])/self.wview[3] - 1
		
		P = numpy.matrix(self.wproj)
		M = numpy.matrix(self.wmodel)
		PM = P*M
		## invert the projection and model view matrices, they will be used shortly
		## note the OpenGL returns matrices are in column major format -  the calculations below 
		## are done with this in  mind - this saves the need to transpose the matrices
		if ( self.update_P_inv == True ):
			P = numpy.matrix(self.wproj)
			self.update_P_inv = False
			self.P_inv = P.I
		M = numpy.matrix(self.wmodel)
		M_inv = M.I
		
		#PM_inv = numpy.matrixmultiply(P_inv,M_inv)
		PM_inv = self.P_inv*M_inv
		
		# If the widget is planar (which obviosuly holds), and along z=0, then the following holds
		zNDC = (PM_inv[0,2]*xNDC + PM_inv[1,2]*yNDC + PM_inv[3,2])/(-PM_inv[2,2])
	
		# We need zprime, which is really 'eye_z' in OpenGL lingo
		print "zprimes"
		print 1.0/(xNDC*self.P_inv[0,3]+yNDC*self.P_inv[1,3]+zNDC*self.P_inv[2,3]+self.P_inv[3,3])
		print zprime
		zp = M[0,3]*ex + M[1,3]*ey + M[2,3]*ez + M[3,3]
		print zp
		
		xtrans = (self.P_inv[0,0]*xNDC + self.P_inv[1,0]*yNDC + self.P_inv[2,0]*zNDC + self.P_inv[3,0])*zprime
		xtrans = xtrans - M[0,0]*ex - M[1,0]*ey - M[2,0]*ez
		
		ytrans = (self.P_inv[0,1]*xNDC + self.P_inv[1,1]*yNDC + self.P_inv[2,1]*zNDC + self.P_inv[3,1])*zprime
		ytrans = ytrans - M[0,1]*ex - M[1,1]*ey - M[2,1]*ez
		
		ztrans = (self.P_inv[0,2]*xNDC + self.P_inv[1,2]*yNDC + self.P_inv[2,2]*zNDC + self.P_inv[3,2])*zprime
		ztrans = ytrans - M[0,2]*ex - M[1,2]*ey - M[2,2]*ez
		
		#print xNDC,yNDC,zNDC
		#print xtrans,ytrans,ztrans
		return [xtrans,ytrans,ztrans]
	
	def eyeCoordsDif(self,x1,y1,x2,y2):
		# get x and y normalized device coordinates
		xNDC1 = 2.0*(x1-self.wview[0])/self.wview[2] - 1
		yNDC1 = 2.0*(y1-self.wview[1])/self.wview[3] - 1
		
		xNDC2 = 2.0*(x2-self.wview[0])/self.wview[2] - 1
		yNDC2 = 2.0*(y2-self.wview[1])/self.wview[3] - 1
		
		## invert the projection and model view matrices, they will be used shortly
		## note the OpenGL returns matrices are in column major format -  the calculations below 
		## are done with this in  mind - this saves the need to transpose the matrices
		if ( self.update_P_inv == True ):
			P = numpy.matrix(self.wproj)
			self.update_P_inv = False
			self.P_inv = P.I
		M = numpy.matrix(self.wmodel)
		M_inv = M.I
		
		#PM_inv = numpy.matrixmultiply(P_inv,M_inv)
		PM_inv = self.P_inv*M_inv
		
		# If the widget is planar (which obviosuly holds), and along z=0, then the following holds
		zNDC1 = (PM_inv[0,2]*xNDC1 + PM_inv[1,2]*yNDC1 + PM_inv[3,2])/(-PM_inv[2,2])
		#zNDC2 = (PM_inv[0,2]*xNDC2 + PM_inv[1,2]*yNDC2 + PM_inv[3,2])/(-PM_inv[2,2])
	
		# We need zprime, which is really 'eye_z' in OpenGL lingo
		zprime1 = 1.0/(xNDC1*self.P_inv[0,3]+yNDC1*self.P_inv[1,3]+zNDC1*self.P_inv[2,3]+self.P_inv[3,3])
		zprime2 = 1.0/(xNDC2*self.P_inv[0,3]+yNDC2*self.P_inv[1,3]+zNDC1*self.P_inv[2,3]+self.P_inv[3,3])

		ex1 = (self.P_inv[0,0]*xNDC1 + self.P_inv[1,0]*yNDC1 + self.P_inv[2,0]*zNDC1+self.P_inv[3,0])*zprime1;
		ey1 = (self.P_inv[0,1]*xNDC1 + self.P_inv[1,1]*yNDC1 + self.P_inv[2,1]*zNDC1+self.P_inv[3,1])*zprime1;
		#ez1 = (self.P_inv[0,2]*xNDC1 + self.P_inv[1,2]*yNDC1 + self.P_inv[2,2]*zNDC1+self.P_inv[3,2])*zprime1;
		
		ex2 = (self.P_inv[0,0]*xNDC2 + self.P_inv[1,0]*yNDC2 + self.P_inv[2,0]*zNDC1+self.P_inv[3,0])*zprime2;
		ey2 = (self.P_inv[0,1]*xNDC2 + self.P_inv[1,1]*yNDC2 + self.P_inv[2,1]*zNDC1+self.P_inv[3,1])*zprime2;
		#ez2 = (self.P_inv[0,2]*xNDC2 + self.P_inv[1,2]*yNDC2 + self.P_inv[2,2]*zNDC1+self.P_inv[3,2])*zprime2;
		
		#return [ex2-ex1,ey2-ey1,ez2-ez1] # ez2-ez1 should be zero
		return [ex2-ex1,ey2-ey1]
	def eyeCoords(self,x,y):
		
		
		# get x and y normalized device coordinates
		xNDC = 2.0*(x-self.wview[0])/self.wview[2] - 1
		yNDC = 2.0*(y-self.wview[1])/self.wview[3] - 1
		
		P = numpy.matrix(self.wproj)
		M = numpy.matrix(self.wmodel)
		PM = P*M
		## invert the projection and model view matrices, they will be used shortly
		## note the OpenGL returns matrices are in column major format -  the calculations below 
		## are done with this in  mind - this saves the need to transpose the matrices
		if ( self.update_P_inv == True ):
			P = numpy.matrix(self.wproj)
			self.update_P_inv = False
			self.P_inv = P.I
		M = numpy.matrix(self.wmodel)
		M_inv = M.I
		
		#PM_inv = numpy.matrixmultiply(P_inv,M_inv)
		PM_inv = self.P_inv*M_inv
		
		# If the widget is planar (which obviosuly holds), and along z=0, then the following holds
		zNDC = (PM_inv[0,2]*xNDC + PM_inv[1,2]*yNDC + PM_inv[3,2])/(-PM_inv[2,2])
	
		# We need zprime, which is really 'eye_z' in OpenGL lingo
		zprime = 1.0/(xNDC*self.P_inv[0,3]+yNDC*self.P_inv[1,3]+zNDC*self.P_inv[2,3]+self.P_inv[3,3])
		
		
		ex = (self.P_inv[0,0]*xNDC + self.P_inv[1,0]*yNDC + self.P_inv[2,0]*zNDC+self.P_inv[3,0])*zprime;
		ey = (self.P_inv[0,1]*xNDC + self.P_inv[1,1]*yNDC + self.P_inv[2,1]*zNDC+self.P_inv[3,1])*zprime;
		ez = (self.P_inv[0,2]*xNDC + self.P_inv[1,2]*yNDC + self.P_inv[2,2]*zNDC+self.P_inv[3,2])*zprime;
		
		#print xNDC,yNDC,zNDC
		#print ex,ey,ez
		return [ex,ey,ez]
		
	def mouseinwin(self,x,y):
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
		if ( self.update_P_inv == True ):
			P = numpy.matrix(self.wproj)
			self.update_P_inv = False
			self.P_inv = P.I
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

		return (xcoord + self.qwidget.width()/2.0, 0.5*self.qwidget.height()-ycoord)
	
	def toolTipEvent(self,event):
		if ( self.childreceiver != None ):
			# this means this class already knows that the mouse event is in the child
			# that is being displayed
			self.childreceiver.toolTip(event)
			self.childreceiver = None
			return
		l=self.mouseinwin(event.x(),self.parent.height()-event.y())
		cw=self.qwidget.childAt(l[0],l[1])
		if cw == None: 
			QtGui.QToolTip.hideText()
			self.genTexture = True
			self.updateTexture()
			return
	
		p1 = QtCore.QPoint(event.x(),event.y())
		p2 = self.parent.mapToGlobal(p1)
		QtGui.QToolTip.showText(p2,cw.toolTip())
		#QtGui.QToolTip.showText(p1,cw.toolTip())
	
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
			l=self.mouseinwin(event.x(),self.parent.height()-event.y())
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
		l=self.mouseinwin(event.x(),self.parent.height()-event.y())
		cw=self.qwidget.childAt(l[0],l[1])
		if cw == None: return
		gp=self.qwidget.mapToGlobal(QtCore.QPoint(l[0],l[1]))
		lp=cw.mapFromGlobal(gp)
		if (isinstance(cw,QtGui.QComboBox)):
			print "it's a combo"
			#cw.showPopup()
		else:
			qme=QtGui.mouseDoubleClickEvent(event.type(),lp,event.button(),event.buttons(),event.modifiers())
			QtCore.QCoreApplication.sendEvent(cw,qme)
		self.genTexture = True
		self.updateTexture()
		
	#def get_depth_for_height(self,height_plane):
		#return self.parent.get_depth_for_height(height_plane)
	
	def mousePressEvent(self, event):
		if event.modifiers() == Qt.ShiftModifier:
			#qme=QtGui.QMouseEvent(event.Type(),QtCore.QPoint(l[0],l[1]),event.button(),event.buttons(),event.modifiers())
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
				
			l=self.mouseinwin(event.x(),self.parent.height()-event.y())
			cw=self.qwidget.childAt(l[0],l[1])
			#print cw
			if cw == None: return
			print cw.objectName()
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
				widget.drawframe = False
				self.e2children.append(widget)
				self.e2children[0].is_child = True
			else:
				qme=QtGui.QMouseEvent( event.type(),lp,event.button(),event.buttons(),event.modifiers())
				if (self.is_child):
					print self.qwidget
					print self.qwidget.currentIndex()
					#self.qwidget.commitData(self.qwidget.parent())
					#print self.qwidget.currentText()
					QtCore.QCoreApplication.sendEvent(self.qwidget,qme)
				else:
					QtCore.QCoreApplication.sendEvent(cw,qme)
				
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
			l=self.mouseinwin(event.x(),self.parent.height()-event.y())
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
			# The fix is to only set the genTexture flag when mouse movement without modifiers
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
			
			l=self.mouseinwin(event.x(),self.parent.height()-event.y())
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
					print self.qwidget
					print self.qwidget.currentIndex().row()
					print self.widget_parent
					print self.qwidget.rect().left(),self.qwidget.rect().right(),self.qwidget.rect().top(),self.qwidget.rect().bottom()
					print lp.x(),lp.y()
					self.widget_parent.setCurrentIndex(self.qwidget.currentIndex().row())
					#self.widget_parent.changeEvent(QtCore.QEvent())
					#self.widget_parent.highlighted(self.qwidget.currentIndex().row())
					#self.qwidget.commitData(self.qwidget.parent())
					#print self.qwidget.currentText()
					self.widget_parent.setVisible(True)
					self.widget_parent.setEnabled(True)
					QtCore.QCoreApplication.sendEvent(self.widget_parent,qme)
				else:
					QtCore.QCoreApplication.sendEvent(cw,qme)
			
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
		
	def getsubtendingangle(self,a,b):
		sinaeb = a[0]*b[1]-a[1]*b[0]
		cosaeb = a[0]*b[0]+a[1]*b[1]
		
		return atan2(sinaeb,cosaeb)

class EMFXImage(QtOpenGL.QGLWidget):
	"""A QT widget for rendering EMData objects. It can display single 2D or 3D images 
	or sets of 2D images.
	"""
	def __init__(self, parent=None):
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True);
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		
		self.fov = 2*180*atan2(1,5)/pi
		self.imtex=0
		self.spinang=0.0
		self.insang=0.0
		self.gq=0			# quadric object for cylinders, etc
		self.framedl=0		# display list for an image frame
		
		self.data=None
		self.datasize=(1,1)
		self.scale=1.0
		self.minden=0
		self.maxden=1.0
		self.mindeng=0
		self.maxdeng=1.0
		self.origin=(0,0)
		self.nperrow=6
		self.nshow=0
		self.mousedrag=None
		
		self.setMouseTracking(True)
		
		self.inspector=None
		self.inspectorl=None
		self.inspector3=None
		
		self.current = None
		self.previous = None
	
		self.qwidgets = []
		self.qwidgets.append(EMQtWidgetDrawer(self))
		self.qwidgets.append(EMQtWidgetDrawer(self))
	
	def setData(self,data):
		"""You may pass a single 2D image, a list of 2D images or a single 3D image"""
		self.data=data
		if data==None:
			self.updateGL()
			return
		
		# If we have a list of 2D images
		if isinstance(data,list) :
			self.minden=data[0].get_attr("mean")
			self.maxden=self.minden
			self.mindeng=self.minden
			self.maxdeng=self.minden
			for i in data:
				if i.get_zsize()!=1 :
					self.data=None
					self.updateGL()
					return
				mean=i.get_attr("mean")
				sigma=i.get_attr("sigma")
				m0=i.get_attr("minimum")
				m1=i.get_attr("maximum")
			
				self.minden=min(self.minden,max(m0,mean-3.0*sigma))
				self.maxden=max(self.maxden,min(m1,mean+3.0*sigma))
				self.mindeng=min(self.mindeng,max(m0,mean-5.0*sigma))
				self.maxdeng=max(self.maxdeng,min(m1,mean+5.0*sigma))
		# If we have a single 2D image
		elif data.get_zsize()==1:
			mean=data.get_attr("mean")
			sigma=data.get_attr("sigma")
			m0=data.get_attr("minimum")
			m1=data.get_attr("maximum")
			
			self.minden=max(m0,mean-3.0*sigma)
			self.maxden=min(m1,mean+3.0*sigma)
			self.mindeng=max(m0,mean-5.0*sigma)
			self.maxdeng=min(m1,mean+5.0*sigma)

			self.datasize=(data.get_xsize(),data.get_ysize())
		# if we have a single 3D image
		elif data.get_zsize()>1 :
			pass
		# Someone passed something wierd
		else :
			self.data=None
			self.updateGL()
			return
		
		self.showInspector()		# shows the correct inspector if already open
		self.updateGL()
		
	def setDenRange(self,x0,x1):
		"""Set the range of densities to be mapped to the 0-255 pixel value range"""
		self.minden=x0
		self.maxden=x1
		self.updateGL()
		
	def get_depth_for_height(self, height):
		# This function returns the width and height of the renderable 
		# area at the origin of the data volume
		depth = height/(2.0*tan(self.fov/2.0*pi/180.0))
		return depth
	
	def setOrigin(self,x,y):
		"""Set the display origin within the image"""
		self.origin=(x,y)
		self.updateGL()
		
	def setScale(self,newscale):
		"""Adjusts the scale of the display. Tries to maintain the center of the image at the center"""
		if isinstance(self.data,list) :
			yo=self.height()-self.origin[1]-1
			self.origin=(newscale/self.scale*(self.width()/2+self.origin[0])-self.width()/2,newscale/self.scale*(self.height()/2+yo)-self.height()/2)
		else : self.origin=(newscale/self.scale*(self.width()/2+self.origin[0])-self.width()/2,newscale/self.scale*(self.height()/2+self.origin[1])-self.height()/2)
		self.scale=newscale
		self.updateGL()
		
	def setDenMin(self,val):
		self.minden=val
		self.updateGL()
		
	def setDenMax(self,val):
		self.maxden=val
		self.updateGL()

	def setNPerRow(self,val):
		self.nperrow=val
		self.updateGL()
		
	def setNShow(self,val):
		self.nshow=val
		self.updateGL()

	def initializeGL(self):
		glClearColor(0,0,0,0)
		
		
		glLightfv(GL_LIGHT0, GL_AMBIENT, [0.1, 0.1, 0.1, 1.0])
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_POSITION, [0.1,.1,1.,0.])
	
		glEnable(GL_LIGHTING)
		glEnable(GL_LIGHT0)
		glEnable(GL_DEPTH_TEST)
		
		glEnable(GL_NORMALIZE)
	
	def paintGL(self):
		#print "in main paint"
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		self.makeCurrent()
		
		
		# define the texture used to render the image on the screen
		if self.inspector:
			if ( self.qwidgets[0].qwidget == None ):
				print "setting Q widget"
				self.qwidgets[0].setQtWidget(self.inspector)
				self.qtc = QtCore.QCoreApplication
				self.fd = QtGui.QFileDialog(self,"Open File")
				self.fd.show()
				#self.fd.hide()
				self.qwidgets[1].setQtWidget(self.fd)
				self.qwidgets[1].cam.setCamX(-100)
			
			glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
			glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
			glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE)
			for i in self.qwidgets:
				glPushMatrix()
				i.paintGL()
				glPopMatrix()
		
		#print "exiting main paint"
		
	def timer(self):
		pass
		#self.updateGL()
	
	def resizeGL(self, width, height):
		
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
		
		self.updateGL()
	
	def showInspector(self,force=0):
		if not force and self.inspector==None and self.inspectorl==None and self.inspector3==None : return
		if isinstance(self.data,list) or isinstance(self.data,tuple) :
			if self.inspector : self.inspector.hide()
			if self.inspector3 : self.inspector3.hide()
			if not self.inspectorl : self.inspectorl=EMImageMxInspector2D(self)
			self.inspectorl.setLimits(self.mindeng,self.maxdeng,self.minden,self.maxden)
#			self.inspectorl.show()
		elif self.data.get_zsize()==1 :
			if self.inspectorl : self.inspectorl.hide()
			if self.inspector3 : self.inspector3.hide()
			if not self.inspector : self.inspector=EMImageInspector2D(self)
			self.inspector.setLimits(self.mindeng,self.maxdeng,self.minden,self.maxden)
			self.inspector.show()
			self.inspector.hide()
			print "told gen texture"
			self.qwidgets[0].genTexture = True
			self.qwidgets[0].updateTexture()
		else:
			pass	# 3d not done yet
	
	def mousePressEvent(self, event):
		if event.button()==Qt.MidButton:
			self.showInspector(1)
		elif event.button()==Qt.RightButton:
			intercepted = False
			for i in self.qwidgets:
				if ( i.isinwin(event.x(),self.height()-event.y()) ):
					i.mousePressEvent(event)
					intercepted = True
					return
			if intercepted == False:
				self.mousedrag=(event.x(),event.y())
		elif event.button()==Qt.LeftButton:
#			self.mousedrag=(event.x(),event.y())
			app=QtGui.QApplication.instance()
			if self.inspector :
				for i in self.qwidgets:
					if ( i.isinwin(event.x(),self.height()-event.y()) ):
						i.mousePressEvent(event)
						self.updateGL()
						return
#				print app.sendEvent(self.inspector.childAt(l[0],l[1]),qme)
	
	def mouseMoveEvent(self, event):
		if self.mousedrag:
			self.origin=(self.origin[0]+self.mousedrag[0]-event.x(),self.origin[1]-self.mousedrag[1]+event.y())
			self.mousedrag=(event.x(),event.y())
			self.update()
			return
		else :
#			self.mousedrag=(event.x(),event.y())
			if self.inspector :
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
			
#				print app.sendEvent(self.inspector.childAt(l[0],l[1]),qme)
				#print qme.x(),qme.y(),l,gp.x(),gp.y()
		
	def mouseReleaseEvent(self, event):
		if event.button()==Qt.RightButton:
			self.mousedrag=None
			return
		elif event.button()==Qt.LeftButton:
#			self.mousedrag=(event.x(),event.y())
			if self.inspector :
				for i in self.qwidgets:
					if ( i.isinwin(event.x(),self.height()-event.y()) ):
						i.mouseReleaseEvent(event)
						self.updateGL()
						return
					
		
	def mouseDoubleClickEvent(self, event):
		if self.inspector :
			for i in self.qwidgets:
				if ( i.isinwin(event.x(),self.height()-event.y()) ):
					i.mouseReleaseEvent(event)
					self.updateGL()
					return
		
		
	def wheelEvent(self, event):
		if self.inspector :
			for i in self.qwidgets:
					if ( i.isinwin(event.x(),self.height()-event.y()) ):
						i.wheelEvent(event)
						self.updateGL()
						return

	def toolTipEvent(self, event):
		if self.inspector :
			for i in self.qwidgets:
					if ( i.isinwin(event.x(),self.height()-event.y()) ):
						i.toolTipEvent(event)
						self.updateGL()
						return
		
		QtGui.QToolTip.hideText()
		

	def dragMoveEvent(self,event):
		print "received drag move event"
		
	def event(self,event):
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
		print "hoverEvent"
		if self.inspector :
			for i in self.qwidgets:
					if ( i.isinwin(event.x(),self.height()-event.y()) ):
						i.hoverEvent(event)
						break
		self.updateGL()

class EMImageMxInspector2D(QtGui.QWidget):
	def __init__(self,target) :
		QtGui.QWidget.__init__(self,None)
		self.target=target
		
		self.vboxlayout = QtGui.QVBoxLayout(self)
		self.vboxlayout.setMargin(0)
		self.vboxlayout.setSpacing(6)
		self.vboxlayout.setObjectName("vboxlayout")
		
		self.hist = ImgHistogram(self)
		self.hist.setObjectName("hist")
		self.vboxlayout.addWidget(self.hist)
		
		self.hboxlayout = QtGui.QHBoxLayout()
		self.hboxlayout.setMargin(0)
		self.hboxlayout.setSpacing(6)
		self.hboxlayout.setObjectName("hboxlayout")
		self.vboxlayout.addLayout(self.hboxlayout)
		
		self.lbl = QtGui.QLabel("#/row:")
		self.lbl.setAlignment(Qt.AlignRight|Qt.AlignVCenter)
		self.hboxlayout.addWidget(self.lbl)
		
		self.nrow = QtGui.QSpinBox(self)
		self.nrow.setObjectName("nrow")
		self.nrow.setRange(1,50)
		self.nrow.setValue(6)
		self.hboxlayout.addWidget(self.nrow)
		
		self.lbl = QtGui.QLabel("N:")
		self.lbl.setAlignment(Qt.AlignRight|Qt.AlignVCenter)
		self.hboxlayout.addWidget(self.lbl)
		
		self.imgn = QtGui.QSpinBox(self)
		self.imgn.setObjectName("imgn")
		self.imgn.setRange(-1,50)
		self.imgn.setValue(0)
		self.imgn.setSpecialValueText("All")
		self.hboxlayout.addWidget(self.imgn)
		
		self.scale = ValSlider(self,(0.1,5.0),"Mag:")
		self.scale.setObjectName("scale")
		self.scale.setValue(1.0)
		self.vboxlayout.addWidget(self.scale)
		
		self.mins = ValSlider(self,label="Min:")
		self.mins.setObjectName("mins")
		self.vboxlayout.addWidget(self.mins)
		
		self.maxs = ValSlider(self,label="Max:")
		self.maxs.setObjectName("maxs")
		self.vboxlayout.addWidget(self.maxs)
		
		self.brts = ValSlider(self,(-1.0,1.0),"Brt:")
		self.brts.setObjectName("brts")
		self.vboxlayout.addWidget(self.brts)
		
		self.conts = ValSlider(self,(0.0,1.0),"Cont:")
		self.conts.setObjectName("conts")
		self.vboxlayout.addWidget(self.conts)
		
		self.lowlim=0
		self.highlim=1.0
		self.busy=0
		
		QtCore.QObject.connect(self.nrow, QtCore.SIGNAL("valueChanged(int)"), target.setNPerRow)
		QtCore.QObject.connect(self.imgn, QtCore.SIGNAL("valueChanged(int)"), target.setNShow)
		QtCore.QObject.connect(self.scale, QtCore.SIGNAL("valueChanged"), target.setScale)
		QtCore.QObject.connect(self.mins, QtCore.SIGNAL("valueChanged"), self.newMin)
		QtCore.QObject.connect(self.maxs, QtCore.SIGNAL("valueChanged"), self.newMax)
		QtCore.QObject.connect(self.brts, QtCore.SIGNAL("valueChanged"), self.newBrt)
		QtCore.QObject.connect(self.conts, QtCore.SIGNAL("valueChanged"), self.newCont)


	def newMin(self,val):
		if self.busy : return
		self.busy=1
		self.target.setDenMin(val)

		self.updBC()
		self.busy=0
		
	def newMax(self,val):
		if self.busy : return
		self.busy=1
		self.target.setDenMax(val)
		self.updBC()
		self.busy=0
	
	def newBrt(self,val):
		if self.busy : return
		self.busy=1
		self.updMM()
		self.busy=0
		
	def newCont(self,val):
		if self.busy : return
		self.busy=1
		self.updMM()
		self.busy=0

	def updBC(self):
		b=0.5*(self.mins.value+self.maxs.value-(self.lowlim+self.highlim))
		c=(self.mins.value-self.maxs.value)/(2.0*(self.lowlim-self.highlim))
		self.brts.setValue(-b)
		self.conts.setValue(1.0-c)
		
	def updMM(self):
		x0=((self.lowlim+self.highlim)/2.0-(self.highlim-self.lowlim)*(1.0-self.conts.value+self.brts.value))
		x1=((self.lowlim+self.highlim)/2.0+(self.highlim-self.lowlim)*(1.0-self.conts.value-self.brts.value))
		self.mins.setValue(x0)
		self.maxs.setValue(x1)
		self.target.setDenRange(x0,x1)
		
	def setHist(self,hist,minden,maxden):
		self.hist.setData(hist,minden,maxden)

	def setLimits(self,lowlim,highlim,curmin,curmax):
		self.lowlim=lowlim
		self.highlim=highlim
		self.mins.setRange(lowlim,highlim)
		self.maxs.setRange(lowlim,highlim)
		self.mins.setValue(curmin)
		self.maxs.setValue(curmax)


class EMImageInspector2D(QtGui.QWidget):
	def __init__(self,target) :
		QtGui.QWidget.__init__(self,None)
		self.target=target
		
		self.vboxlayout = QtGui.QVBoxLayout(self)
		self.vboxlayout.setMargin(0)
		self.vboxlayout.setSpacing(6)
		self.vboxlayout.setObjectName("vboxlayout")
		
		self.hist = ImgHistogram(self)
		self.hist.setObjectName("hist")
		self.vboxlayout.addWidget(self.hist)
		
		self.scale = ValSlider(self,(0.1,5.0),"Mag:")
		self.scale.setObjectName("scale")
		self.scale.setValue(1.0)
		self.vboxlayout.addWidget(self.scale)
		
		self.mins = ValSlider(self,label="Min:")
		self.mins.setObjectName("mins")
		self.vboxlayout.addWidget(self.mins)
		
		self.combo = QtGui.QComboBox(self)
		
		for i in range(0,10):
			self.combo.addItem(str(i))
		self.vboxlayout.addWidget(self.combo)
		
		self.maxs = ValSlider(self,label="Max:")
		self.maxs.setObjectName("maxs")
		self.vboxlayout.addWidget(self.maxs)
		
		self.brts = ValSlider(self,(-1.0,1.0),"Brt:")
		self.brts.setObjectName("brts")
		self.vboxlayout.addWidget(self.brts)
		
		self.conts = ValSlider(self,(0.0,1.0),"Cont:")
		self.conts.setObjectName("conts")
		self.vboxlayout.addWidget(self.conts)
		
		self.lowlim=0
		self.highlim=1.0
		self.busy=0
		
		QtCore.QObject.connect(self.scale, QtCore.SIGNAL("valueChanged"), target.setScale)
		QtCore.QObject.connect(self.mins, QtCore.SIGNAL("valueChanged"), self.newMin)
		QtCore.QObject.connect(self.maxs, QtCore.SIGNAL("valueChanged"), self.newMax)
		QtCore.QObject.connect(self.brts, QtCore.SIGNAL("valueChanged"), self.newBrt)
		QtCore.QObject.connect(self.conts, QtCore.SIGNAL("valueChanged"), self.newCont)
		QtCore.QObject.connect(self.combo, QtCore.SIGNAL("currentIndexChanged(QString)"), self.setCombo)
		
	def setCombo(self,val):
		print val
		print "yeealllow"
	
	def newMin(self,val):
		if self.busy : return
		self.busy=1
		self.target.setDenMin(val)

		self.updBC()
		self.busy=0
		
	def newMax(self,val):
		if self.busy : return
		self.busy=1
		self.target.setDenMax(val)
		self.updBC()
		self.busy=0
	
	def newBrt(self,val):
		if self.busy : return
		self.busy=1
		self.updMM()
		self.busy=0
		
	def newCont(self,val):
		if self.busy : return
		self.busy=1
		self.updMM()
		self.busy=0

	def updBC(self):
		b=0.5*(self.mins.value+self.maxs.value-(self.lowlim+self.highlim))
		c=(self.mins.value-self.maxs.value)/(2.0*(self.lowlim-self.highlim))
		self.brts.setValue(-b)
		self.conts.setValue(1.0-c)
		
	def updMM(self):
		x0=((self.lowlim+self.highlim)/2.0-(self.highlim-self.lowlim)*(1.0-self.conts.value+self.brts.value))
		x1=((self.lowlim+self.highlim)/2.0+(self.highlim-self.lowlim)*(1.0-self.conts.value-self.brts.value))
		self.mins.setValue(x0)
		self.maxs.setValue(x1)
		self.target.setDenRange(x0,x1)
		
	def setHist(self,hist,minden,maxden):
		self.hist.setData(hist,minden,maxden)

	def setLimits(self,lowlim,highlim,curmin,curmax):
		self.lowlim=lowlim
		self.highlim=highlim
		self.mins.setRange(lowlim,highlim)
		self.maxs.setRange(lowlim,highlim)
		self.mins.setValue(curmin)
		self.maxs.setValue(curmax)
	
	#def event(self,event):
		##if event.spontaneous(): return False
		#print "hi from here"
		#print event.type()
		#if event.type() in mouseEvents:
			#cw=self.childAt(event.x,event.y)
			#if cw == None:
				#print "no child"
				#return False
			#gp=self.mapToGlobal(QtCore.QPoint(event.x,event.y))
			#lp=cw.mapFromGlobal(gp)
			#qme=QtGui.QMouseEvent(event.type(),lp,event.button,event.buttons,event.modifiers)
			#print "returning event"
			#return cw.event(qme)
		#elif event.Type() == QtCore.QEvent.Wheel:
			#cw=self.childAt(event.x,event.y)
			#if cw == None:
				#print "no child"
				#return False
			#gp=self.mapToGlobal(QtCore.QPoint(event.x,event.y))
			#lp=cw.mapFromGlobal(gp)
			#qme=QtGui.QWheelEvent(lp,event.delta,event.buttons,event.modifiers,event.orientation)
			#return cw.event(qme)
		
	
		##print event.t
		#print event.Type()
		#print mouseEvents
		#print "no event"
		#return False
	
	#def mouseMoveEvent(self,event):
		#print str(event.__dict__)
		#print "received event"
		
	#def mousePressEvent(self, event):
		#print "RECIEVED EVENT"
		
# 	def mousePressEvent(self, event):
# 		print str(event.__dict__)
# 		QtGui.QWidget.mousePressEvent(self,event)
		

# This is just for testing, of course
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = EMFXImage()
	if len(sys.argv)==1 : window.setData(test_image(size=(512,512)))
	else :
		a=EMData.read_images(sys.argv[1])
		if len(a)==1 : window.setData(a[0])
		else : window.setData(a)
	window2 = EMParentWin(window)
	window2.show()
	
	ti=QtCore.QTimer()
	ti.setInterval(50.)
	QtCore.QObject.connect(ti, QtCore.SIGNAL("timeout()"), window.timer)
	ti.start()

	
#	w2=QtGui.QWidget()
#	w2.resize(256,128)
	
#	w3=ValSlider(w2)
#	w3.resize(256,24)
#	w2.show()
	
	sys.exit(app.exec_())
