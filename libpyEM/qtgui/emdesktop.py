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
		self.bgob=ob2dimage(self,QtGui.QPixmap.grabWindow(self.appscreen.winId(),0.0,0.0,self.sysdesktop.width(),self.sysdesktop.height()-30),self.aspect)
		
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
	
	def timeout(self):
		self.time+=1
		if (self.time<150) :
			self.makeCurrent()
			glMatrixMode(GL_PROJECTION)
			glLoadIdentity()
			glFrustum(-self.aspect,self.aspect, -1.,1., 15.-self.time/20.0,15.)
			glTranslatef(0.,0.,-14.99)
			self.updateGL()
		else:
			self.updateGL()
		
	
	def initializeGL(self):
		glClearColor(0,0,0,0)
		
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
		
		glDisable(GL_LIGHTING)
		self.bgob.render()
		glEnable(GL_LIGHTING)
		glEnable(GL_LIGHT0)
		glCallList(self.volcubedl)
		
		if self.time>150 :
			glPushMatrix()
			glScalef(.25,.25,.25)
			glTranslate(-2.5,0.0,2.0)
			glRotate(self.time,1.0,-.2,0.0)
			glCallList(self.framedl)
			glPopMatrix()

	def mouseinwin(self,x,y):
		x00=self.mc00[0]
		x01=self.mc01[0]-x00
		x10=self.mc10[0]-x00
		x11=self.mc11[0]-x00
		y00=self.mc00[1]
		y01=self.mc01[1]-y00
		y10=self.mc10[1]-y00
		y11=self.mc11[1]-y00
		x-=x0015
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
		glLoadIdentity()
		glFrustum(-self.aspect,self.aspect, -1.,1., 5.,15.)
		glTranslatef(0.,0.,-14.9)
		
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
		
	def mousePressEvent(self, event):
		if event.button()==Qt.MidButton:
			self.showInspector(1)
		elif event.button()==Qt.LeftButton:
#			self.mousedrag=(event.x(),event.y())
			app=QtGui.QApplication.instance()
			#if self.inspector :
				#l=self.mouseinwin(event.x(),event.y())
				#cw=self.inspector.childAt(l[0],l[1])
				#gp=self.inspector.mapToGlobal(QtCore.QPoint(l[0],l[1]))
				#qme=QtGui.QMouseEvent(event.Type(),QtCore.QPoint(l[0]-cw.x(),l[1]-cw.y()),gp,event.button(),event.buttons(),event.modifiers())
##				self.inspector.mousePressEvent(qme)
				#cw.mousePressEvent(qme)
##				print app.sendEvent(self.inspector.childAt(l[0],l[1]),qme)
				#print qme.x(),qme.y(),l,gp.x(),gp.y()
	
	def mouseMoveEvent(self, event):
		if event.buttons()==Qt.LeftButton:
			pass
#			self.mousedrag=(event.x(),event.y())
			#if self.inspector :
				#l=self.mouseinwin(event.x(),event.y())
				#cw=self.inspector.childAt(l[0],l[1])
				#gp=self.inspector.mapToGlobal(QtCore.QPoint(l[0],l[1]))
				#qme=QtGui.QMouseEvent(event.Type(),QtCore.QPoint(l[0]-cw.x(),l[1]-cw.y()),gp,event.button(),event.buttons(),event.modifiers())
				#cw.mouseMoveEvent(qme)
				#print qme.x(),qme.y(),l,gp.x(),gp.y()
		
	def mouseReleaseEvent(self, event):
		if event.button()==Qt.LeftButton:
			pass
			#if self.inspector :
				#l=self.mouseinwin(event.x(),event.y())
				#cw=self.inspector.childAt(l[0],l[1])
				#gp=self.inspector.mapToGlobal(QtCore.QPoint(l[0],l[1]))
				#qme=QtGui.QMouseEvent(event.Type(),QtCore.QPoint(l[0]-cw.x(),l[1]-cw.y()),gp,event.button(),event.buttons(),event.modifiers())
##				self.inspector.mousePressEvent(qme)
				#cw.mouseReleaseEvent(qme)
##				print app.sendEvent(self.inspector.childAt(l[0],l[1]),qme)
				#print qme.x(),qme.y(),l,gp.x(),gp.y()

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
	
	def render(self):
		if not self.pixmap : return
		#glPushMatrix()
		#glTranslate(-2.5,0,0)
		#glRotate(self.insang,0.1,1.0,0.0)
		#glTranslate(2.5,0,0)
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
		
		#glPopMatrix()


class ob2d:
	def __init__(self,target):
		self.mc00=(0,0,0)
		self.mc10=(0,0,0)
		self.mc11=(0,0,0)
		self.mc01=(0,0,0)
		self.pixmap=None
		self.widget=None
		self.target=target
		self.itex=None

	def __del__(self):
		target.deleteTexture(self.itex)


	def setWidget(self,widget,region=None):
		self.widget=widget
		self.region=region
		self.update()

	def update(self):
#		self.pixmap=QtGui.QPixmap.grabWidget(self.widget)
		if self.region : self.pixmap=QtGui.QPixmap.grabWindow(self.widget.winId(),self.region.x(),self.region.y(),self.region.width(),self.region.height())
		else : self.pixmap=QtGui.QPixmap.grabWindow(self.widget.winId())
		self.aspect=float(self.pixmap.width())/self.pixmap.height()
#		self.pixmap.save("test.jpg")
		self.target.makeCurrent()
		self.itex=self.target.bindTexture(self.pixmap)
		print "---> ",self.itex

	def render(self):
		if not self.pixmap : return
		glPushMatrix()
		#glTranslate(-2.5,0,0)
		#glRotate(self.insang,0.1,1.0,0.0)
		#glTranslate(2.5,0,0)
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
		
		wmodel=glGetDoublev(GL_MODELVIEW_MATRIX)
		wproj=glGetDoublev(GL_PROJECTION_MATRIX)
		wview=glGetIntegerv(GL_VIEWPORT)

		self.mc00=gluProject(-self.aspect,-1.0,0.0,wmodel,wproj,wview)
		self.mc10=gluProject(self.aspect,-1.0,0.0,wmodel,wproj,wview)
		self.mc11=gluProject(self.aspect, 1.0,0.0,wmodel,wproj,wview)
		self.mc01=gluProject(-self.aspect, 1.0,0.0,wmodel,wproj,wview)
		glPopMatrix()


if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = EMDesktop()
#	window.showFullScreen()
	window.app.exec_()
