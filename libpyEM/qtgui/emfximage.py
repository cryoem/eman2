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
	def __init__(self, parent=None, qwidget=None):
		self.parent = parent
		self.drawframe = False
		self.mapcoords = True
		self.debugcoords = True
		self.itex = 0
		self.genTexture = True
		self.click_debug = False
		self.cam = Camera2(self)
		self.cam.setCamTrans('default_z',-self.parent.get_depth_for_height(height_plane))
		self.cam.motionRotate(0,0)
		self.borderwidth = 3.0
		self.glbasicobjects = EMBasicObjects()
		self.setQtWidget(qwidget)
		self.P_inv = None
		self.update_P_inv = True
	
	def set_update_P_inv(self,val=True):
		self.update_P_inv = val
	
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
		
		if (self.qwidget == None ) : return
		
		self.cam.position()
		if ( self.itex == 0 or self.genTexture == True ) : 
			self.parent.deleteTexture(self.itex)
			self.genTexture = False
			self.itex = self.parent.bindTexture(QtGui.QPixmap.grabWidget(self.qwidget))
		
		self.wmodel=glGetDoublev(GL_MODELVIEW_MATRIX)
		self.wproj=glGetDoublev(GL_PROJECTION_MATRIX)
		self.wview=glGetIntegerv(GL_VIEWPORT)
		
		glPushMatrix()
		glTranslate(-self.qwidget.width()/2.0,-self.qwidget.height()/2.0,0.0)
		glEnable(GL_TEXTURE_2D)
		glBindTexture(GL_TEXTURE_2D,self.itex)
		glBegin(GL_QUADS)
		glTexCoord2f(0.,0.)
		glVertex(0,0)
		glTexCoord2f(1.,0.)
		glVertex( self.qwidget.width(),0)
		glTexCoord2f(1.,1.)
		glVertex( self.qwidget.width(), self.qwidget.height())
		glTexCoord2f(0.,1.)
		glVertex(0, self.qwidget.height())
		glEnd()
		glDisable(GL_TEXTURE_2D)
		glPopMatrix()
		#self.parent.deleteTexture(self.itex)
	
		if ( self.mapcoords ):
			
	
			self.mc00=gluProject(-self.qwidget.width()/2.0,-self.qwidget.height()/2.0,0.,self.wmodel,self.wproj,self.wview)
			self.mc10=gluProject( self.qwidget.width()/2.0,-self.qwidget.height()/2.0,0.,self.wmodel,self.wproj,self.wview)
			self.mc11=gluProject( self.qwidget.width()/2.0, self.qwidget.height()/2.0,0.,self.wmodel,self.wproj,self.wview)
			self.mc01=gluProject(-self.qwidget.width()/2.0, self.qwidget.height()/2.0,0.,self.wmodel,self.wproj,self.wview)
	
			if ( self.debugcoords ):
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
				
		if ( self.drawframe ):
			glPushMatrix()
			glCallList(self.glbasicobjects.getFrameDL())
			glPopMatrix()
	
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
	
	def wheelEvent(self,event):
		
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.wheelEvent(event)
		else:
			l=self.mouseinwin(event.x(),self.parent.height()-event.y())
			cw=self.qwidget.childAt(l[0],l[1])
			if cw == None: return
			gp=self.qwidget.mapToGlobal(QtCore.QPoint(l[0],l[1]))
			lp=cw.mapFromGlobal(gp)
			qme=QtGui.QWheelEvent(lp,event.delta(),event.buttons(),event.modifiers(),event.orientation())
			print  QtCore.QCoreApplication.sendEvent(cw,qme)
			self.genTexture = True
	
	def mousePressEvent(self, event):
		if event.modifiers() == Qt.ShiftModifier:
			#qme=QtGui.QMouseEvent(event.Type(),QtCore.QPoint(l[0],l[1]),event.button(),event.buttons(),event.modifiers())
			self.cam.mousePressEvent(event)
		else:
			l=self.mouseinwin(event.x(),self.parent.height()-event.y())
			cw=self.qwidget.childAt(l[0],l[1])
			#print cw
			if cw == None: return
			print cw.objectName()
			gp=self.qwidget.mapToGlobal(QtCore.QPoint(l[0],l[1]))
			lp=cw.mapFromGlobal(gp)
			if (isinstance(cw,QtGui.QComboBox)):
				print "it's a combo"
				cw.showPopup()
				#QtGui.QComboBox(cw).showPopup()
			else:
				qme=QtGui.QMouseEvent( event.type(),lp,event.button(),event.buttons(),event.modifiers())
				print  QtCore.QCoreApplication.sendEvent(cw,qme)
				
			self.genTexture = True
		
	def mouseMoveEvent(self,event):
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.mouseMoveEvent(event)
		else:
			l=self.mouseinwin(event.x(),self.parent.height()-event.y())
			cw=self.qwidget.childAt(l[0],l[1])
			if cw == None: return
			gp=self.qwidget.mapToGlobal(QtCore.QPoint(l[0],l[1]))
			lp=cw.mapFromGlobal(gp)
			qme=QtGui.QMouseEvent(event.type(),lp,event.button(),event.buttons(),event.modifiers())
			print  QtCore.QCoreApplication.sendEvent(cw,qme)
			self.genTexture = True

	def mouseReleaseEvent(self,event):
		if event.modifiers() == Qt.ShiftModifier:
			self.cam.mouseReleaseEvent(event)
		else:
			l=self.mouseinwin(event.x(),self.parent.height()-event.y())
			cw=self.qwidget.childAt(l[0],l[1])
			if cw == None: return
			gp=self.qwidget.mapToGlobal(QtCore.QPoint(l[0],l[1]))
			lp=cw.mapFromGlobal(gp)
			if (isinstance(cw,QtGui.QComboBox)):
				print "it's a combo"
				cw.showPopup()
			else:
				qme=QtGui.QMouseEvent(event.type(),lp,event.button(),event.buttons(),event.modifiers())
				print  QtCore.QCoreApplication.sendEvent(cw,qme)
			self.genTexture = True
		
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
		
		self.inspector=None
		self.inspectorl=None
		self.inspector3=None
	
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
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		#glTranslated(0.0, 0.0, -10.0)
		
		if not self.data : return
		
		# define the texture used to render the image on the screen
		if not self.imtex :
#			 glDeleteTextures([self.imtex])
			self.imtex=glGenTextures(1)				# 'name' of the image display texture
#			print "tex=",self.imtex
			glBindTexture(GL_TEXTURE_2D,self.imtex)
			glTexImage2D(GL_TEXTURE_2D,0,GL_LUMINANCE8,512,512,0,GL_LUMINANCE,GL_UNSIGNED_BYTE,"\000"*(512*512))	# start with an empty texture
			glBindTexture(GL_TEXTURE_2D,0)

		# get a new Quadric object for drawing cylinders, spheres, etc
		if not self.gq:
			self.gq=gluNewQuadric()
			gluQuadricDrawStyle(self.gq,GLU_FILL)
			gluQuadricNormals(self.gq,GLU_SMOOTH)
			gluQuadricOrientation(self.gq,GLU_OUTSIDE)
			gluQuadricTexture(self.gq,GL_FALSE)

		glPushMatrix()
		# define a square frame of rounded components bordering -1,-1 to 1,1
		if not self.framedl :
			self.framedl=glGenLists(1)
			glNewList(self.framedl,GL_COMPILE)
			glPushMatrix()
			glRotate(90.,1.,0.,0.)
			glTranslate(1.02,0.,-1.02)
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
			
			glTranslate(1.02,1.02,0.)
			gluSphere(self.gq,.02,6,6)
			glTranslate(-2.04,0.,0.)
			gluSphere(self.gq,.02,6,6)
			glTranslate(0.,-2.04,0.)
			gluSphere(self.gq,.02,6,6)
			glTranslate(2.04,0.,0.)
			gluSphere(self.gq,.02,6,6)
			glEndList()

		if isinstance(self.data,list) :
			if self.nshow==-1 :
				glPixelZoom(1.0,-1.0)
				n=len(self.data)
				x,y=-self.origin[0],self.height()-self.origin[1]-1
				for i in range(n):
					w=int(min(self.data[i].get_xsize()*self.scale,self.width()))
					h=int(min(self.data[i].get_ysize()*self.scale,self.height()))
#					print i,x,y,w,h
					if x>0 and x<self.width() and y>0 and y<self.height() :
						a=self.data[i].render_amp8(0,0,w,h,(w-1)/4*4+4,self.scale,0,255,self.minden,self.maxden,1.0,2)
						glRasterPos(x,y)
						glDrawPixels(w,h,GL_LUMINANCE,GL_UNSIGNED_BYTE,a)
					elif x+w>0 and y-h<self.height() and x<self.width() and y>0:
						if x<0 : 
							x0=-x/self.scale
							x1=w+x
						else : 
							x0=0
							x1=w
						if y>self.height()-1 : y1=h-y+self.height()-1
						else : y1=h
						x0,x1,y1=int(x0),int(x1),int(y1)
						a=self.data[i].render_amp8(x0,0,x1,y1,(x1-1)/4*4+4,self.scale,0,255,self.minden,self.maxden,1.0,2)
#						a=self.data[i].render_amp8(int(-x),int(y-self.height()),min(w,int(x+w)),int(h-y+self.height()),min(w-1,int(x+w-1))/4*4+4,self.scale,0,255,self.minden,self.maxden,2)
						if x<0 : xx=0
						else : xx=x
						if y>=self.height() : yy=self.height()-1
						else : yy=y
						glRasterPos(xx,yy)
						glDrawPixels(x1,y1,GL_LUMINANCE,GL_UNSIGNED_BYTE,a)
						
					
					if (i+1)%self.nperrow==0 : 
						y-=h+2.0
						x=-self.origin[0]
					else: x+=w+2.0
			else:
				a=self.data[self.nshow].render_amp8(int(self.origin[0]/self.scale),int(self.origin[1]/self.scale),self.width(),self.height(),(self.width()-1)/4*4+4,self.scale,0,255,self.minden,self.maxden,1.0,2)
				glRasterPos(0,self.height()-1)
				glPixelZoom(1.0,-1.0)
				glDrawPixels(self.width(),self.height(),GL_LUMINANCE,GL_UNSIGNED_BYTE,a)
				hist=array.array("I")
				hist.fromstring(a[-1024:])
				if self.inspectorl : self.inspectorl.setHist(hist,self.minden,self.maxden)
		else :
			a=self.data.render_amp8(int(self.origin[0]/self.scale),int(self.origin[1]/self.scale),512,512,512,self.scale,0,255,self.minden,self.maxden,1.0,2)
			hist=array.array("I")
			hist.fromstring(a[-1024:])
			if self.inspector : self.inspector.setHist(hist,self.minden,self.maxden)

			glEnable(GL_TEXTURE_2D)
			glBindTexture(GL_TEXTURE_2D,self.imtex)
#			glTexImage2D(GL_TEXTURE_2D,0,GL_LUMINANCE8,512,512,0,GL_LUMINANCE,GL_UNSIGNED_BYTE,a)	# start with an empty texture
			glTexSubImage2D(GL_TEXTURE_2D,0,0,0,512,512,GL_LUMINANCE,GL_UNSIGNED_BYTE,a)
			
			glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP)
			glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP)
			glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)
			glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)
			glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE)

			# This is the textured square showing the actual image
			glPushMatrix()
			glRotate(self.spinang,1.,.2,0.)
			glBegin(GL_QUADS)
			glTexCoord2f(0.,.999)
			glVertex(-1.,-1.)
			glTexCoord2f(.999,.999)
			glVertex( 1.,-1.)
			glTexCoord2f(.999,0.)
			glVertex( 1., 1.)
			glTexCoord2f(0.,0.)
			glVertex(-1., 1.)
			glEnd()

			glBindTexture(GL_TEXTURE_2D,0)
			
			glDisable(GL_TEXTURE_2D)

			glColor(.2,.2,.8)
			glMaterial(GL_FRONT,GL_AMBIENT,(.2,.2,.8,1.0))
			glMaterial(GL_FRONT,GL_SPECULAR,(.8,.8,.8,1.0))
			glMaterial(GL_FRONT,GL_SHININESS,50.0)
			
			glPushMatrix()
			glCallList(self.framedl)
			glPopMatrix()
			
			glPopMatrix()
			
		glPopMatrix()
		if self.inspector:
			if ( self.qwidgets[0].qwidget == None ):
				print "setting Q widget"
				self.qwidgets[0].setQtWidget(self.inspector)
				self.qtc = QtCore.QCoreApplication
				self.fd = QtGui.QFileDialog(self,"Open File")
				#self.fd.show()
				self.qwidgets[1].setQtWidget(self.fd)
				self.qwidgets[1].cam.setCamX(-100)
			glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)
			glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)
			glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE)
			for i in self.qwidgets:
				glPushMatrix()
				i.paintGL()
				glPopMatrix()
				
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
		print depth
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
			print "told gen texture"
			self.qwidgets[0].genTexture = True
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
					break
			if intercepted == False:
				self.mousedrag=(event.x(),event.y())
		elif event.button()==Qt.LeftButton:
#			self.mousedrag=(event.x(),event.y())
			app=QtGui.QApplication.instance()
			if self.inspector :
				for i in self.qwidgets:
					if ( i.isinwin(event.x(),self.height()-event.y()) ):
						i.mousePressEvent(event)
						break
#				print app.sendEvent(self.inspector.childAt(l[0],l[1]),qme)
		self.updateGL()
	
	def mouseMoveEvent(self, event):
		if self.mousedrag:
			self.origin=(self.origin[0]+self.mousedrag[0]-event.x(),self.origin[1]-self.mousedrag[1]+event.y())
			self.mousedrag=(event.x(),event.y())
			self.update()
		else :
#			self.mousedrag=(event.x(),event.y())
			if self.inspector :
				for i in self.qwidgets:
					if ( i.isinwin(event.x(),self.height()-event.y()) ):
						i.mouseMoveEvent(event)
						break
		
		self.updateGL()
#				print app.sendEvent(self.inspector.childAt(l[0],l[1]),qme)
				#print qme.x(),qme.y(),l,gp.x(),gp.y()
		
	def mouseReleaseEvent(self, event):
		if event.button()==Qt.RightButton:
			self.mousedrag=None
		elif event.button()==Qt.LeftButton:
#			self.mousedrag=(event.x(),event.y())
			if self.inspector :
				for i in self.qwidgets:
					if ( i.isinwin(event.x(),self.height()-event.y()) ):
						i.mouseReleaseEvent(event)
						break
					
		self.updateGL()
					
	def wheelEvent(self, event):
		if self.inspector :
			for i in self.qwidgets:
					if ( i.isinwin(event.x(),self.height()-event.y()) ):
						i.wheelEvent(event)
						break
		self.updateGL()

class EMFxTexture:
	def __init__(self,parent):
		self.parent = parent

class ImgHistogram(QtGui.QWidget):
	def __init__(self,parent):
		QtGui.QWidget.__init__(self,parent)
		self.brush=QtGui.QBrush(Qt.black)
		
		self.font=QtGui.QFont("Helvetica", 10);
		self.probe=None
		self.histdata=None
		self.setMinimumSize(QtCore.QSize(258,128))
	
	def setData(self,data,minden,maxden):
		self.histdata=data
#		self.norm=max(self.histdata)
		self.norm=0
		self.minden=minden
		self.maxden=maxden
		for i in self.histdata: self.norm+=i*i
		self.norm-=max(self.histdata)**2
		self.norm=sqrt(self.norm/255)*3.0
		self.total=sum(self.histdata)
		if self.norm==0 : self.norm=1.0
		if self.total==0 : self.histdata=None
		self.update()
	
	def paintEvent (self, event):
		if not self.histdata : return
		p=QtGui.QPainter()
		p.begin(self)
		p.setPen(Qt.darkGray)
		for i,j in enumerate(self.histdata):
			p.drawLine(i,127,i,127-j*126/self.norm)
		
		# If the user has dragged, we need to show a value
		if self.probe :
			p.setPen(Qt.blue)
			p.drawLine(self.probe[0]+1,0,self.probe[0]+1,127-self.probe[1]*126/self.norm)
			p.setPen(Qt.darkRed)
			p.drawLine(self.probe[0]+1,127,self.probe[0]+1,127-self.probe[1]*126/self.norm)
			p.setFont(self.font)
			p.drawText(200,20,"x=%d"%(self.probe[0]))
			p.drawText(200,34,"%1.2f"%(self.probe[0]/255.0*(self.maxden-self.minden)+self.minden))
			p.drawText(200,48,"y=%d"%(self.probe[1]))
			p.drawText(200,62,"%1.2f%%"%(100.0*float(self.probe[1])/self.total))
		
		p.setPen(Qt.black)
		p.drawRect(0,0,257,128)
		p.end()

	def mousePressEvent(self, event):
		if event.button()==Qt.LeftButton:
			x=max(min(event.x()-1,255),0)
			self.probe=(x,self.histdata[x])
			self.update()
			
	def mouseMoveEvent(self, event):
		if event.buttons()&Qt.LeftButton:
			x=max(min(event.x()-1,255),0)
			self.probe=(x,self.histdata[x])
			self.update()
	
	def mouseReleaseEvent(self, event):
		self.probe=None

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
