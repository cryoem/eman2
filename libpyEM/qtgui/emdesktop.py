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

import PyQt4
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
#from OpenGL import GL,GLU,GLUT
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
from valslider import ValSlider
from math import *
from EMAN2 import *
import EMAN2
import sys
import numpy
from emimageutil import ImgHistogram
from weakref import WeakKeyDictionary
from pickle import dumps,loads
from PyQt4.QtCore import QTimer

from emfloatingwidgets import EMFloatingWidgetsCore, EMBasicOpenGLObjects
from emglobjects import Camera


class EMDesktop(QtOpenGL.QGLWidget):
	"""An opengl windowing system, which can contain other EMAN2 widgets and 3-D objects.
	"""
	def __init__(self):
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True);
		QtOpenGL.QGLWidget.__init__(self,fmt,None)
		
		self.imtex=0
		self.spinang=0.0
		self.insang=0.0
		self.gq=0			# quadric object for cylinders, etc
		self.sqframedl=1	# display list for an image frame
		self.bigcubedl=0	# display list for frame around volume
	
		self.obs2d=[]
		self.obs3d=[]
		
		self.app=QtGui.QApplication.instance()
		self.sysdesktop=self.app.desktop()
		self.appscreen=self.sysdesktop.screen(self.sysdesktop.primaryScreen())
		self.aspect=float(self.appscreen.width())/self.appscreen.height()
		self.appwidth = self.appscreen.width()
		self.appheight = self.appscreen.height()
		self.fov = 30
		self.zopt = (1.0/tan(self.fov/2.0*pi/180.0))*self.appheight/2.0
		self.resizeGL(self.appwidth,self.appheight)
		#self.bgob=ob2dimage(self,QtGui.QPixmap.grabWindow(self.appscreen.winId(),0.0,0.0,self.sysdesktop.width(),self.sysdesktop.height()-30),self.aspect)
		self.bgob2=ob2dimage(self,self.readEMAN2Image(),self.aspect)
#		self.setWindowFlags(Qt.FramelessWindowHint)
#		print self.sysdesktop.primaryScreen(),self.sysdesktop.numScreens(),self.appscreen.size().width(),self.sysdesktop.screenGeometry(0).width(),self.sysdesktop.screenGeometry(1).width()
		self.norender=1
		self.show()
		self.move(0,0)
		self.resize(self.appscreen.size())
		self.norender=0
		
		self.timer = QTimer()
		QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.timeout)
		self.timer.start(10)
		self.time=0
	
		self.setMouseTracking(True)
		
		self.floatwidget = EMFloatingWidgetsCore(self)
		
		self.glbasicobjects = EMBasicOpenGLObjects()
		self.borderwidth=100.0
		self.cam = Camera()
		
	def get_depth_for_height(self, height):
		return 0
		# This function returns the width and height of the renderable 
		# area at the origin of the data volume
		depth = height/(2.0*tan(self.fov/2.0*pi/180.0))
		return depth
	
	def get_render_dims_at_depth(self, depth):
		return 0
		# This function returns the width and height of the renderable 
		# area at the origin of the data volume
		height = -2*tan(self.fov/2.0*pi/180.0)*(depth)
		width = self.aspect*height
		return [width,height]
	
	def cylinderToFrom(self,To,From):
		#print "To",To[0],To[1]
		#print "From",From[0],From[1]
		# draw a cylinder to To from From
		dx = To[0] - From[0]
		dy = To[1] - From[1]
		dz = To[2] - From[2]
		
		length = sqrt(dx*dx+dy*dy+dz*dz)
		
		#angle = 180.0*atan2(dy,dx)/pi
		print From[0],From[1],From[2]
		glTranslatef(From[0],From[1],From[2])
		#glRotated(90.0+angle,0.,0.,1.0)
		#glRotated(90.,1.,0.,0.)
		glTranslated(0.0,0.0,length/2.0)
		glScaled(self.borderwidth,self.borderwidth,length)
		glTranslated(0.0,0.0,-0.5)
		glCallList(self.glbasicobjects.getCylinderDL())

	
	def readEMAN2Image(self):
		self.p = QtGui.QPixmap("EMAN2.0.big.jpg")
		return self.p
	
	def timeout(self):
		self.time+=1
		if (self.time<150) :
			strat = 2
			if ( strat == 1  ):
				self.zNear = 15.-self.time/20.0
				self.zFar = 15.0
				self.fov = atan(2.0/self.zNear)*180.0/pi
				self.makeCurrent()
				glMatrixMode(GL_PROJECTION)
				glLoadIdentity()
				#glFrustum(-self.aspect,self.aspect, -1.,1., 15.-self.time/20.0,15.)
				gluPerspective(self.fov,self.aspect,self.zNear,self.zFar)
				self.updateGL()
			else:
				#self.zNear = 2*self.zopt - self.time/150.0*self.zopt
				#self.zFar = 2*self.zopt
				#glMatrixMode(GL_PROJECTION)
				#glLoadIdentity()
				#print self.zNear, self.zFar
				#gluPerspective(self.fov,self.aspect,self.zNear,self.zFar)
				#glMatrixMode(GL_MODELVIEW)
				pass
				self.updateGL()
				
				
				
		#elif (self.time<300) :
			#self.fov =  atan(2.0/self.zNear)*180.0/pi + 100.0*(self.time-150.0)/150.
			#gluPerspective(self.fov,self.aspect,self.zNear,self.zFar)
			#self.updateGL()
		else:
			self.updateGL()
		
	
	def initializeGL(self):
		glClearColor(0,0,0,0)
		glEnable(GL_NORMALIZE)
		
		# get a new Quadric object for drawing cylinders, spheres, etc
		if not self.gq:
			self.gq=gluNewQuadric()
			gluQuadricDrawStyle(self.gq,GLU_FILL)
			gluQuadricNormals(self.gq,GLU_SMOOTH)
			gluQuadricOrientation(self.gq,GLU_OUTSIDE)
			gluQuadricTexture(self.gq,GL_FALSE)
		
		# Precompile a displaylist for a frame
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

	def paintGL(self):
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
#		glLoadIdentity()
#		glTranslated(0.0, 0.0, -10.0)
		if self.norender : return
		
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
	
		glDisable(GL_LIGHTING)
		
		
		#self.bgob.render()
		glPushMatrix()
		if (self.time < 150):
			z = self.zopt + float(self.time)/150.0*self.zopt
			#print z
			glTranslatef(0.,0.,-z)
		else:
			#print -2*self.zopt+0.1
			glTranslatef(0.,0.,-2*self.zopt+0.1)
		#print "translating back", -2*self.zopt+1
		#glTranslatef(0,0,-16)
		glScalef(self.appheight/2.0,self.appheight/2.0,1.0)
		#print self.appheight
		self.bgob2.render()
		glPopMatrix()
		glEnable(GL_LIGHTING)
		glEnable(GL_LIGHT0)
		glPushMatrix()
		glScalef(self.appheight/2.0,self.appheight/2.0,1.0)
		glCallList(self.volcubedl)
		glPopMatrix()
		
		
		glColor(.9,.2,.8)
		# this is a nice light blue color (when lighting is on)
		# and is the default color of the frame
		glMaterial(GL_FRONT,GL_AMBIENT,(.2,.2,.8,1.0))
		glMaterial(GL_FRONT,GL_SPECULAR,(.8,.8,.8,1.0))
		glMaterial(GL_FRONT,GL_SHININESS,50.0)
		if self.time>150 :
			dx = (self.appwidth - self.p.width())/2.0
			dy = (self.appheight - self.p.height())/2.0
			
			
			glPushMatrix()
			glTranslatef(0.,0,-self.p.height())
			glTranslatef(-dx,-dy, -self.zopt)
			#glScalef(.25,.25,.25)
			#glTranslate(5.5,-2.5,2.0)
			glRotate(self.time,1.0,0.,0.0)
			#glScalef(self.p.width(),self.p.height(),1.0)
			#glScalef(self.bgob2.asp(),1.0,1.0)
			self.bgob2.render2()
			glPopMatrix()
			
			glPushMatrix()
			glTranslatef(-dx, dy, -self.zopt)
			self.bgob2.render2()
			glPopMatrix()

			glPushMatrix()
			dx = (self.appwidth - self.p.width())/2.0
			dy = (self.appheight - self.p.height())/2.0
			glTranslatef(dx,-dy, -3.0/2.0*self.zopt)
			self.bgob2.render2()
			glPopMatrix()

			glPushMatrix()
			glTranslatef(0.,0.,-self.zopt)
			self.floatwidget.render()
			glPopMatrix()
			
			
			glColor(.9,.2,.8)
			# this is a nice light blue color (when lighting is on)
			# and is the default color of the frame
			glMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,(.2,.2,.8,1.0))
			glMaterial(GL_FRONT,GL_SPECULAR,(.8,.8,.8,1.0))
			glMaterial(GL_FRONT,GL_SHININESS,50.0)
			glPushMatrix()
			to = [-dx,dy,-2*self.zopt]
			f = [-dx,dy,-self.zopt]
			self.cylinderToFrom(to,f)
			glPopMatrix()
	def resizeGL(self, width, height):
	
		glEnable(GL_LIGHTING)
		glEnable(GL_LIGHT0)
		glEnable(GL_DEPTH_TEST)
		glLightfv(GL_LIGHT0, GL_AMBIENT, [0.9, 0.9, 0.9, 1.0])
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_POSITION, [0.5,0.7,11.,0.])


		side = min(width, height)
#		glViewport((width - side) / 2, (height - side) / 2, side, side)
		glViewport(0,0,self.width(),self.height())
		glMatrixMode(GL_PROJECTION)
		#self.zNear = 7.5
		#self.zFar = 15.0
		#self.fov = atan(2.0/self.zNear)*180.0/pi
		glLoadIdentity()
		
		self.zNear = self.zopt
		self.zFar = 2*self.zopt
				#glMatrixMode(GL_PROJECTION)
				#glLoadIdentity()
				#print self.zNear, self.zFar
		gluPerspective(self.fov,self.aspect,self.zNear,self.zFar)
				#glMatrixMode(GL_MODELVIEW)
				#self.updateGL()
		
		#gluPerspective(self.fov,self.aspect,self.zNear,self.zFar)
		
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
	
	def mouseReleaseEvent(self, event):
		if event.button()==Qt.LeftButton:
			pass
		
	def mousePressEvent(self, event):
		self.floatwidget.mousePressEvent(event)
	
	def mouseMoveEvent(self, event):
		self.floatwidget.mouseMoveEvent(event)
		
	def mouseReleaseEvent(self, event):
		self.floatwidget.mouseReleaseEvent(event)

	def mouseDoubleClickEvent(self, event):
		self.floatwidget.mouseDoubleClickEvent(event)

	def wheelEvent(self, event):
		self.floatwidget.wheelEvent(event)

	def toolTipEvent(self, event):
		self.floatwidget.wheelEvent(event)
		QtGui.QToolTip.hideText()
		

	def dragMoveEvent(self,event):
		print "received drag move event, but I don't do anything about it :("
		
	def event(self,event):
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
		self.floatwidget.hoverEvent(event)

class ob2dimage:
	def __init__(self,target,pixmap,aspect):
		self.pixmap=pixmap
		self.target=target
		self.aspect=aspect
		self.target.makeCurrent()
		self.itex=self.target.bindTexture(self.pixmap)

	def __del__(self):
		target.deleteTexture(self.itex)

	def setWidget(self,widget,region=None):
		return

	def update(self):
		return
	
	def width(self):
		return self.pixmap.width()
	
	def height(self):
		return self.pixmap.height()
	
	def asp(self):
		return (1.0*self.width())/self.height()
	
	def render2(self):
		if not self.pixmap : return
		glPushMatrix()
		glScalef(self.pixmap.width()/2.0,self.pixmap.height()/2.0,1.0)
		glColor(1.0,1.0,1.0)
		glEnable(GL_TEXTURE_2D)
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE)
		glBindTexture(GL_TEXTURE_2D,self.itex)
		glBegin(GL_QUADS)
		glTexCoord2f(0.,0.)
		glVertex(-1.0,-1.0)
		glTexCoord2f(.999,0.)
		glVertex( 1.0,-1.0)
		glTexCoord2f(.999,0.999)
		glVertex( 1.0, 1.0)
		glTexCoord2f(0.,.999)
		glVertex(-1.0, 1.0)
		glEnd()
		glPopMatrix()
		glDisable(GL_TEXTURE_2D)
	
	def render(self):
		if not self.pixmap : return
		#glPushMatrix()
		#glTranslate(-2.5,0,0)
		#glRotate(self.insang,0.1,1.0,0.0)
		#glTranslate(2.5,0,0)
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP)
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP)
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE)
		glColor(1.0,1.0,1.0)
		glEnable(GL_TEXTURE_2D)
		glBindTexture(GL_TEXTURE_2D,self.itex)
		glBegin(GL_QUADS)
		glTexCoord2f(0.,0.)
		glVertex(-self.aspect,-1.0)
		glTexCoord2f(.999,0.)
		glVertex( self.aspect,-1.0)
		glTexCoord2f(.999,0.999)
		glVertex( self.aspect, 1.0)
		glTexCoord2f(0.,.999)
		glVertex(-self.aspect, 1.0)
		glEnd()
		glDisable(GL_TEXTURE_2D)
		#glPopMatrix()

if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = EMDesktop()
#	window.showFullScreen()
	window.app.exec_()
