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
from math import sqrt

from emimage3dobject import Camera

import Numeric
import LinearAlgebra
import sys
import array

class EMBasicObjects:
	def __init__(self):
		self.framedl = 0
		
		self.gq=gluNewQuadric()
		gluQuadricDrawStyle(self.gq,GLU_FILL)
		gluQuadricNormals(self.gq,GLU_SMOOTH)
		gluQuadricOrientation(self.gq,GLU_OUTSIDE)
		gluQuadricTexture(self.gq,GL_FALSE)

	def getFrameDL(self):
		# draws a frame using Quadrics
		if ( self.framedl == 0 ):
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


class EMQtWidgetDrawer:
	def __init__(self, parent=None, qwidget=None):
		self.glparent = parent
		self.qwidget = qwidget
		self.drawframe = True
		self.mapcoords = True
		self.debugcoords = True
		self.itex = 0
		self.click_debug = False
		self.cam = Camera()
		self.cam.motionRotate(25,25)
		
		self.glbasicobjects = EMBasicObjects()
		
	def setQtWidget(self, widget, delete_current = False):
		
		if ( delete_current ):
			self.qwidget.deleteLater()
		
		self.qwidget = widget
	
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
		
		self.cam.position()
		self.wmodel=glGetDoublev(GL_MODELVIEW_MATRIX)
		self.wproj=glGetDoublev(GL_PROJECTION_MATRIX)
		self.wview=glGetIntegerv(GL_VIEWPORT)
		self.itex = self.glparent.bindTexture(QtGui.QPixmap.grabWidget(self.qwidget))
		glEnable(GL_TEXTURE_2D)
		glBindTexture(GL_TEXTURE_2D,self.itex)
		#print glGetTexLevelParameteriv(GL_TEXTURE_2D,0,GL_TEXTURE_WIDTH),glGetTexLevelParameteriv(GL_TEXTURE_2D,0,GL_TEXTURE_HEIGHT)
		glBegin(GL_QUADS)
		glTexCoord2f(0.,0.)
		glVertex(-1.,-1.)
		glTexCoord2f(1.,0.)
		glVertex( 1.,-1.)
		glTexCoord2f(1.,1.)
		glVertex( 1., 1.)
		glTexCoord2f(0.,1.)
		glVertex(-1., 1.)
		glEnd()
		glDisable(GL_TEXTURE_2D)
		self.glparent.deleteTexture(self.itex)
	
		if ( self.mapcoords ):
			
			
	
			self.mc00=gluProject(-1.,-1.,0.,self.wmodel,self.wproj,self.wview)
			self.mc10=gluProject( 1.,-1.,0.,self.wmodel,self.wproj,self.wview)
			self.mc11=gluProject( 1., 1.,0.,self.wmodel,self.wproj,self.wview)
			self.mc01=gluProject(-1., 1.,0.,self.wmodel,self.wproj,self.wview)
	
			#print "CHECKING",self.mc00[0],self.mc00[1],self.mc00[2]
			#print gluUnProject(self.mc00[0],self.mc00[1],self.mc00[2],self.wmodel,self.wproj,self.wview)
			#self.mouseinwin(self.mc00[0],self.mc00[1])
	
			if ( self.debugcoords ):
				glMatrixMode(GL_PROJECTION)
				glPushMatrix()
				glLoadIdentity()
				glOrtho(0.0,self.glparent.width(),0.0,self.glparent.height(),-1,1)
	
				glMatrixMode(GL_MODELVIEW)
				glPushMatrix()
				
				glLoadIdentity()
				glPointSize(10)
				glColor(1.0,1.0,1.0,1.0)
				glBegin(GL_POINTS)
				glVertex(self.mc00[0], self.mc00[1])
				glVertex(self.mc10[0], self.mc10[1])
				glVertex(self.mc11[0], self.mc11[1])
				glVertex(self.mc01[0], self.mc01[1])
				if (self.click_debug == True):
					glVertex(self.mcd00[0], self.mcd00[1])
					glVertex(self.mcd10[0], self.mcd10[1])
					glVertex(self.mcd11[0], self.mcd11[1])
					glVertex(self.mcd01[0], self.mcd01[1])
					
				glEnd()
				
				glBegin(GL_LINES)
				glVertex(self.mc00[0], self.mc00[1])
				glVertex(self.mc11[0], self.mc11[1])
				glVertex(self.mc10[0], self.mc10[1])
				glVertex(self.mc01[0], self.mc01[1])
				if (self.click_debug == True):
					glVertex(self.mcd00[0], self.mcd00[1])
					glVertex(self.mcd10[0], self.mcd10[1])
					glVertex(self.mcd11[0], self.mcd11[1])
					glVertex(self.mcd01[0], self.mcd01[1])
					
				glEnd()
				
				glPopMatrix()
				
				glMatrixMode(GL_PROJECTION)
				# pop the temporary orthographic matrix from the GL_PROJECTION stack
				glPopMatrix()
				glMatrixMode(GL_MODELVIEW)
				
		if ( self.drawframe ):
			glPushMatrix()
			glCallList(self.glbasicobjects.getFrameDL())
			glPopMatrix()

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
	
	def mouseinwin(self,x,y):
		# to determine the mouse coordinates in the window we carefully perform
		# linear algebra similar to what's done in gluUnProject

		# the problem is determining what the z coordinate of the mouse event should have
		# been, given that we know that  thewidget itself is located in the x,y plane, along z=0.
		
		# get x and y normalized device coordinates
		xNDC = 2.0*(x-self.wview[0])/self.wview[2] - 1
		yNDC = 2.0*(y-self.wview[1])/self.wview[3] - 1
		
		# invert the projection and model view matrices, they will be used shortly
		# note the OpenGL returns matrices in column major format -  the calculations below 
		# are done with this in  mind - this saves the need to transpose the matrices
		P = Numeric.array(self.wproj)
		M = Numeric.array(self.wmodel)
		P_inv = LinearAlgebra.inverse(P)
		M_inv = LinearAlgebra.inverse(M)
		PM_inv = Numeric.matrixmultiply(P_inv,M_inv)
		
		# If the widget is planar (which obviosuly holds), and along z=0, then the following holds
		zNDC = (PM_inv[0,2]*xNDC + PM_inv[1,2]*yNDC + PM_inv[3,2])/(-PM_inv[2,2])
	
		# We need zprime, which is really 'eye_z' in OpenGL lingo
		zprime = 1.0/(xNDC*P_inv[0,3]+yNDC*P_inv[1,3]+zNDC*P_inv[2,3]+P_inv[3,3])
		
		# Now we compute the x and y coordinates - these are precisely what we're after
		xcoord = zprime*(xNDC*PM_inv[0,0]+yNDC*PM_inv[1,0]+zNDC*PM_inv[2,0]+PM_inv[3,0])
		ycoord = zprime*(xNDC*PM_inv[0,1]+yNDC*PM_inv[1,1]+zNDC*PM_inv[2,1]+PM_inv[3,1])

		return ((xcoord+1)/2.0*self.qwidget.width(),(1-(ycoord+1)/2.0)*self.qwidget.height())
		
	def mousePressEvent(self, event):
		l=self.mouseinwin(event.x(),self.glparent.height()-event.y())
		cw=self.qwidget.childAt(l[0],l[1])
		gp=self.qwidget.mapToGlobal(QtCore.QPoint(l[0],l[1]))
		qme=QtGui.QMouseEvent(event.Type(),QtCore.QPoint(l[0]-cw.x(),l[1]-cw.y()),gp,event.button(),event.buttons(),event.modifiers())
		cw.mousePressEvent(qme)
		
	def mouseMoveEvent(self,event):
		l=self.mouseinwin(event.x(),self.glparent.height()-event.y())
		cw=self.qwidget.childAt(l[0],l[1])
		gp=self.qwidget.mapToGlobal(QtCore.QPoint(l[0],l[1]))
		qme=QtGui.QMouseEvent(event.Type(),QtCore.QPoint(l[0]-cw.x(),l[1]-cw.y()),gp,event.button(),event.buttons(),event.modifiers())
		#self.inspector.mousePressEvent(qme)
		cw.mouseMoveEvent(qme)
		
	def mouseReleaseEvent(self,event):
		l=self.mouseinwin(event.x(),self.glparent.height()-event.y())
		cw=self.qwidget.childAt(l[0],l[1])
		gp=self.qwidget.mapToGlobal(QtCore.QPoint(l[0],l[1]))
		qme=QtGui.QMouseEvent(event.Type(),QtCore.QPoint(l[0]-cw.x(),l[1]-cw.y()),gp,event.button(),event.buttons(),event.modifiers())
		#self.inspector.mousePressEvent(qme)
		cw.mouseReleaseEvent(qme)
		#print app.sendEvent(self.inspector.childAt(l[0],l[1]),qme)
		#print qme.x(),qme.y(),l,gp.x(),gp.y()
		
	def timerEvent(self,event=None):
		# event = None is just here incase anyone ever actually wants to pass the event
		self.cam.motionRotate(.2,.2)
		pass
		
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
	
		self.qwidgetdrawer = EMQtWidgetDrawer(self)
	
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
	
	def paintGL(self):
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		glTranslated(0.0, 0.0, -10.0)
		
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
				
		if self.inspector:
			if ( self.qwidgetdrawer.qwidget == None ):
					self.qwidgetdrawer.setQtWidget(self.inspector)
			
			glPushMatrix()
			glTranslate(-2.5,0.0,0.0)
			self.qwidgetdrawer.paintGL()
			glPopMatrix()

	def isinwin(self,x,y):
		# this works by simple geometry - if the mouse point e is within the four points (a,b,c,d)
		# at the extremities of  the qtwidget, then aed + bec + ced + dea is +/- 360 degrees. 
		# If e is outside the four points then the sum is zero...
		
		ly = self.height() - y
		a = [self.mc00[0]-x, self.mc00[1]-ly]
		b = [self.mc01[0]-x, self.mc01[1]-ly]
		c = [self.mc11[0]-x, self.mc11[1]-ly]
		d = [self.mc10[0]-x, self.mc10[1]-ly]
		
		
		aeb = self.getsubtendingangle(a,b)
		bec = self.getsubtendingangle(b,c)
		ced = self.getsubtendingangle(c,d)
		dea = self.getsubtendingangle(d,a)
		if abs(aeb + bec + ced + dea) > 0.1:
			return True 
		else:
			return False
		
	def getsubtendingangle(self,a,b):
		sinaeb = a[0]*b[1]-a[1]*b[0]
		cosaeb = a[0]*b[0]+a[1]*b[1]
		
		return atan2(sinaeb,cosaeb)

	def mouseinwin(self,x,y):

		x00=self.mc00[0]
		x01=self.mc01[0]-x00
		x10=self.mc10[0]-x00
		x11=self.mc11[0]-x00
		y00=self.mc00[1]
		y01=self.mc01[1]-y00
		y10=self.mc10[1]-y00
		y11=self.mc11[1]-y00
		x-=x00
		y-=y00
		
# 		print "%f,%f  %f,%f  %f,%f  %f,%f"%(x00,y00,x01,y01,x11,y11,x10,y10)
		
		try: xx=(x01*y + x10*y - x11*y - x*y01 - x10*y01 - x*y10 + x01*y10 +
		x*y11 + sqrt(pow(x11*y + x*y01 - x10*(y + y01) + x*y10 + x01*(-y + y10) - x*y11,2) -
		4*(x10*y - x*y10)*(x10*y01 - x11*y01 + x01*(-y10 + y11))))/(2.*((x01 - x11)*y10 + x10*(-y01 + y11)))
		except: xx=x/x10
			
		try: yy=(x01*y + x10*y - x11*y - x*y01 + x10*y01 - x*y10 - x01*y10 +
		x*y11 - sqrt(pow(x11*y + x*y01 - x10*(y + y01) + x*y10 + x01*(-y + y10) - x*y11,2) -
		4*(x10*y - x*y10)*(x10*y01 - x11*y01 + x01*(-y10 + y11))))/(2.*(x10*y01 - x11*y01 + x01*(-y10 + y11)))
		except: yy=y/y01
			
		return (xx*self.inspector.width(),yy*self.inspector.height())
		
		
# 		return (x01*y + x10*y - x11*y - x*y01 - x10*y01 - x*y10 + x01*y10 +
# 		x*y11 + sqrt(pow(x11*y + x*y01 - x10*(y + y01) + x*y10 + x01*(-y + y10) - x*y11,2) -
# 		4*(x10*y - x*y10)*(x10*y01 - x11*y01 + x01*(-y10 + y11))))/(2.*((x01 - x11)*y10 + x10*(-y01 + y11)))


	def timer(self):
		self.spinang+=0.5
		self.insang+=0.2
		if ( self.qwidgetdrawer.qwidget != None ):
			self.qwidgetdrawer.timerEvent()
		self.updateGL()
	
	def resizeGL(self, width, height):
		glLightfv(GL_LIGHT0, GL_AMBIENT, [0.1, 0.1, 0.1, 1.0])
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_POSITION, [0.1,.1,1.,0.])
	
		glEnable(GL_LIGHTING)
		glEnable(GL_LIGHT0)
		glEnable(GL_DEPTH_TEST)


		side = min(width, height)
#		glViewport((width - side) / 2, (height - side) / 2, side, side)
		glViewport(0,0,self.width(),self.height())
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		#glFrustum(-1.*width/height,1.*width/height, -1.,1., 5.,15.)
		
		
		# fov angle is the given by
		self.fov = 2*180*atan2(1,5)/pi
		# aspect ratio is given by
		self.aspect = float(self.width())/float(self.height())
		# this is the same as the glFrustum call above
		gluPerspective(self.fov,self.aspect,5,15)
		
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()


	
# 		glMatrixMode(GL_PROJECTION)
# 		glLoadIdentity()
# 		GLU.gluOrtho2D(0.0,self.width(),0.0,self.height())
# 		glMatrixMode(GL_MODELVIEW)
# 		glLoadIdentity()
		
#		glMatrixMode(GL_PROJECTION)
#		glLoadIdentity()
#		glOrtho(-0.5, +0.5, +0.5, -0.5, 4.0, 15.0)
#		glMatrixMode(GL_MODELVIEW)
	
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
#			self.inspector.hide()
		else:
			pass	# 3d not done yet
	
	def mousePressEvent(self, event):
		if event.button()==Qt.MidButton:
			self.showInspector(1)
		elif event.button()==Qt.RightButton:
			self.mousedrag=(event.x(),event.y())
		elif event.button()==Qt.LeftButton:
#			self.mousedrag=(event.x(),event.y())
			app=QtGui.QApplication.instance()
			if self.inspector :
				if ( self.qwidgetdrawer.isinwin(event.x(),self.height()-event.y()) ):
					self.qwidgetdrawer.mousePressEvent(event)
#				print app.sendEvent(self.inspector.childAt(l[0],l[1]),qme)
	
	def mouseMoveEvent(self, event):
		if self.mousedrag:
			self.origin=(self.origin[0]+self.mousedrag[0]-event.x(),self.origin[1]-self.mousedrag[1]+event.y())
			self.mousedrag=(event.x(),event.y())
			self.update()
		elif event.buttons()==Qt.LeftButton:
#			self.mousedrag=(event.x(),event.y())
			if self.inspector :
				if ( self.qwidgetdrawer.isinwin(event.x(),self.height()-event.y()) ):
					self.qwidgetdrawer.mouseMoveEvent(event)
				
#				print app.sendEvent(self.inspector.childAt(l[0],l[1]),qme)
				#print qme.x(),qme.y(),l,gp.x(),gp.y()
		
	def mouseReleaseEvent(self, event):
		if event.button()==Qt.RightButton:
			self.mousedrag=None
		elif event.button()==Qt.LeftButton:
#			self.mousedrag=(event.x(),event.y())
			if self.inspector :
				if ( self.qwidgetdrawer.isinwin(event.x(),self.height()-event.y()) ):
					self.qwidgetdrawer.mouseReleaseEvent(event)

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
