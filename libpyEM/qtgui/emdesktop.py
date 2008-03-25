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

from emfloatingwidgets import *
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
		
		self.floatwidget = EMDesktopWidgets(self)
		self.floatwidget.suppressUpdateGL = True
		
		self.glbasicobjects = EMBasicOpenGLObjects()
		self.borderwidth=10.0
		self.cam = Camera()
		
		#self.fd = EMDesktopFileDialog(self,"Open File",QtCore.QDir.currentPath(),QtCore.QString("Image files (*.img *.hed *.mrc)"))
		#exit(1)
		self.framedl = 0

	def get_app_width(self):
		return self.appwidth
	
	def get_app_height(self):
		return self.appheight
	
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
	
	def drawFrame(self):
		if self.framedl == 0:
			#print self.appwidth/2.0,self.appheight/2.0,self.zopt
			glCallList(self.glbasicobjects.getCylinderDL())
			length = self.zopt
			self.framedl=glGenLists(1)
			glNewList(self.framedl,GL_COMPILE)
			glPushMatrix()
			glTranslatef(-self.appwidth/2.0,-self.appheight/2.0,0.0)
			glScaled(self.borderwidth,self.borderwidth,length)
			glCallList(self.glbasicobjects.getCylinderDL())
			glPopMatrix()
			glPushMatrix()
			glTranslatef( self.appwidth/2.0,-self.appheight/2.0,0.0)
			glScaled(self.borderwidth,self.borderwidth,length)
			glCallList(self.glbasicobjects.getCylinderDL())
			glPopMatrix()
			
			glPushMatrix()
			glTranslatef( self.appwidth/2.0, self.appheight/2.0,0.0)
			glScaled(self.borderwidth,self.borderwidth,length)
			glCallList(self.glbasicobjects.getCylinderDL())
			glPopMatrix()
			
			glPushMatrix()
			glTranslatef(-self.appwidth/2.0, self.appheight/2.0,0.0)
			glScaled(self.borderwidth,self.borderwidth,length)
			glCallList(self.glbasicobjects.getCylinderDL())
			glPopMatrix()
			
			glEndList()
			
		if self.framedl == 0:
			print "error, frame display list failed to compile"
			exit(1)
		glColor(.9,.2,.8)
		## this is a nice light blue color (when lighting is on)
		## and is the default color of the frame
		glMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,(.2,.2,.8,1.0))
		glMaterial(GL_FRONT,GL_SPECULAR,(.8,.8,.8,1.0))
		glMaterial(GL_FRONT,GL_SHININESS,128.0)
		glCallList(self.framedl)
		
	def readEMAN2Image(self):
		self.p = QtGui.QPixmap("EMAN2.0.big.jpg")
		return self.p
	
	def timeout(self):
		self.time+=1
		self.updateGL()
		
	def initializeGL(self):
		glClearColor(0,0,0,0)
		glEnable(GL_NORMALIZE)
				
		glEnable(GL_LIGHTING)
		glEnable(GL_LIGHT0)
		glEnable(GL_DEPTH_TEST)
		glLightfv(GL_LIGHT0, GL_AMBIENT, [0.9, 0.9, 0.9, 1.0])
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_POSITION, [0.,1,0.,0.])

		# get a new Quadric object for drawing cylinders, spheres, etc
		if not self.gq:
			self.gq=gluNewQuadric()
			gluQuadricDrawStyle(self.gq,GLU_FILL)
			gluQuadricNormals(self.gq,GLU_SMOOTH)
			gluQuadricOrientation(self.gq,GLU_OUTSIDE)
			gluQuadricTexture(self.gq,GL_FALSE)
	
	def paintGL(self):
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
#		glLoadIdentity()
#		glTranslatef(0.0, 0.0, -10.0)
		if self.norender : return
		
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		#self.bgob.render()
		glPushMatrix()
		if (self.time < 150):
			z = self.zopt + float(self.time)/150.0*self.zopt
			#print z
			glTranslatef(0.,0.,-z)
		else:
			#print -2*self.zopt+0.1
			glTranslatef(0.,0.,-2*self.zopt+0.1)
			
		glPushMatrix()
		self.drawFrame()
		glPopMatrix()
		

		glPushMatrix()
		glScalef(self.appheight/2.0,self.appheight/2.0,1.0)
		self.bgob2.render()
		glPopMatrix()
		glPopMatrix()
		
		if self.time>150 :
			dx = (self.appwidth - self.p.width())/2.0
			dy = (self.appheight - self.p.height())/2.0
			
			
			glPushMatrix()
			glTranslatef(self.appwidth/2.0+self.p.width(),-dy, -1.8*self.zopt)
			glScalef(.25,.25,.25)
			glRotate(self.time,1.0,0.,0.0)
			self.bgob2.render2()
			glPopMatrix()
			
			#glPushMatrix()
			#glTranslatef(-dx, dy, -self.zopt)
			#self.bgob2.render2()
			#glPopMatrix()

			glPushMatrix()
			dx = (self.appwidth - self.p.width())/2.0
			dy = (self.appheight - self.p.height())/2.0
			
			dz = (self.time%10000)/5000.0
			if ( dz > 1): dz = 2-dz
			glTranslatef(self.appwidth/2.0,-dy, -self.zopt-self.p.width()/2.0-dz*(self.zopt-self.p.width()))
			glRotatef(-90,0,1,0)
			self.bgob2.render2()
			glPopMatrix()

			glPushMatrix()
			glTranslatef(0.,0.,-self.zopt)
			self.floatwidget.render()
			glPopMatrix()

	def resizeGL(self, width, height):
		side = min(width, height)
		glViewport(0,0,self.width(),self.height())
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		
		self.zNear = self.zopt
		self.zFar = 2*self.zopt
		gluPerspective(self.fov,self.aspect,self.zNear-500,self.zFar)
		
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

class EMDesktopFileDialog:
	'''
	This is just a wrapper for a QtGui Filedialog, that does some extra stuff...
	'''
	def __init__(self, parent,title,directory,wildcard):
		#QtGui.QFileDialog.__init__(parent,message,directory,wildcard)
		self.fd = QtGui.QFileDialog(parent,title,directory,wildcard)
		self.fname = ''
	def show(self):
		self.fd.show()
		
	def hide(self):
		self.fd.hide()
		
	def getfd(self):
		return self.fd
	
	def width(self):
		return self.fd.width()
	
	def height(self):
		return self.fd.height()
	
	def selectedFiles(self):
		return self.fd.selectedFiles()
	
	def deleteLater(self):
		self.fd.deleteLater()
		
	def setEnabled(self,val=True):
		self.fd.setEnabled(val)
	
	def widget(self):
		return self.fd
	
	def childAt(self,x,y):
		return self.fd.childAt(x,y)
	
	def mapToGlobal(self,p):
		return self.fd.mapToGlobal(p)
	
	def currentIndex(self):
		return self.fd.currentIndex()
	
	def mousePressEvent(self, event):
		print "press event"
		self.fd.mousePressEvent(event)
		if event.button()==Qt.LeftButton:
			self.fnames = self.fd.selectedFiles()
			for i in self.fnames: print i
			
	def mouseReleaseEvent(self,event):
		print "release"
		self.fd.mouseReleaseEvent(event)
		if event.button()==Qt.LeftButton:
			self.fnames = self.fd.selectedFiles()
			print "release"
			for i in self.fnames: print i
		
	def mouseMoveEvent(self, event):
		self.fd.mouseMoveEvent(event)
		
	def mouseReleaseEvent(self, event):
		self.fd.mouseReleaseEvent(event)

	def mouseDoubleClickEvent(self, event):
		self.fd.mouseDoubleClickEvent(event)

	def wheelEvent(self, event):
		self.fd.wheelEvent(event)

	def toolTipEvent(self, event):
		self.fd.wheelEvent(event)
		QtGui.QToolTip.hideText()
		
	def dragMoveEvent(self,event):
		print "received drag move event, but I don't do anything about it :("
		
	#def __getattr__(self, attr):
		#""" Delegate access to implementation """
		#return getattr(self.fd, attr)

	#def __setattr__(self, attr, value):
		#""" Delegate access to implementation """
		#return setattr(self.fd, attr, value)
		
	

class EMDesktopWidgets:
	"""A QT widget for rendering EMData objects. It can display single 2D or 3D images 
	or sets of 2D images.
	"""
	def __init__(self, parent=None):
		#print "init"
		self.parent = parent
	
		self.imtex=0
		self.current = None
		self.previous = None
	
		self.initFlag = True
		self.qwidgets = []
		self.imagewidgets = []
		
		self.suppressUpdateGL = False
		
		self.fdxpos = 0
		self.fdypos = 0

	def get_depth_for_height(self, height):
		try: 
			return self.parent.get_depth_for_height(height)
		except:
			print "parent can't get height for depth"
			return 0

	def height(self):
		return self.parent.height()
	
	def width(self):
		return self.parent.width()

	def updateGL(self):
		if not self.suppressUpdateGL:
			try: self.parent.updateGL()
			except: pass

	def addQtWidgetDrawer(self,widget):
		w = EMGLViewQtWidget(self)
		w.setQtWidget(widget)
		self.qwidgets.append(w)
		
		#print "initializeGL done"
	def render(self):
		
		if ( self.initFlag == True ):
			self.initFlag = False
			self.fd = EMDesktopFileDialog(self.parent,"Open File",QtCore.QDir.currentPath(),QtCore.QString("Image files (*.img *.hed *.mrc)"))
			QtCore.QObject.connect(self.fd.getfd(), QtCore.SIGNAL("finished(int)"), self.finished)
			QtCore.QObject.connect(self.fd.getfd(), QtCore.SIGNAL("currentChanged(QString)"), self.changed)
			self.fd.show()
			self.fd.hide()
			self.fdwidget = EMGLViewQtWidget(self.parent)
			self.fdwidget.cam.cam_x = -(self.parent.get_app_width() - self.fd.width())/2.0
			self.fdwidget.cam.cam_y = (self.parent.get_app_height() - self.fd.height())/2.0
			self.qwidgets.append(self.fdwidget)
			self.qwidgets[0].setQtWidget(self.fd)
			self.initFlag = False

		for i in self.qwidgets:
			glPushMatrix()
			i.paintGL()
			glPopMatrix()
			
		for i in self.imagewidgets:
			glPushMatrix()
			i.paintGL()
			glPopMatrix()
		
	def changed(self,file):
		print "changing to", file
		# get the center coord for the next image to the right
		x1 = -self.parent.get_app_width()/2.0 + self.fd.width()
		xwidth = self.parent.get_app_width()/2.0  - self.fd.width()
		
		h = self.parent.get_app_height()
		newposx = x1+xwidth/2.0
		newposy = h/4.0
		yheight = h/2.0
		
		try:
			a=EMData.read_images(str(file))
			if len(a) == 1:
				a = a[0]
				if a.get_zsize() != 1: w = EMGLView3D(self,a)
				else: w = EMGLView2D(self,a)
			else: w = EMGLView2D(self,a)
			
			self.imagewidgets = []
			w.cam.cam_x = newposx
			w.cam.cam_y = newposy
			scalex = xwidth/w.width()
			scaley = yheight/w.height()
			if scalex > scaley: scalex = scaley
			w.cam.scale = scalex
			self.imagewidgets.append(w)
		
			newposy = -h/4.0
		except:
			print "could not open"
			return
		
		
		#print "paintGL done"
	def finished(self,val):
		if ( val == 1 ):
			for i in self.fd.selectedFiles():
				a=EMData.read_images(str(i))
				if len(a) == 1:
					a = a[0]
					if a.get_zsize() != 1:
						w = EMGLView3D(self,a)
						self.qwidgets.append(w)
					else:
						w = EMGLView2D(self,a)
						self.qwidgets.append(w)
				else:
					w = EMGLView2D(self,a)
					self.qwidgets.append(w)
					
	def timer(self):
		pass
		#self.updateGL()
		
	def bindTexture(self,pixmap):
		return self.parent.bindTexture(pixmap)
	
	def deleteTexture(self,val):
		return self.parent.deleteTexture(val)
	
	def get_render_dims_at_depth(self, depth):
		try: return self.parent.get_render_dims_at_depth(depth)
		except:
			print "parent can't get render dims at for depth"
			return

	def resizeEvent(self, width, height):
		for i in self.qwidgets:
			i.set_update_P_inv()
	
	def mousePressEvent(self, event):
		for i in self.qwidgets:
			if ( i.isinwin(event.x(),self.height()-event.y()) ):
				i.mousePressEvent(event)
				intercepted = True
				self.updateGL()
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
