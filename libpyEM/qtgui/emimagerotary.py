#!/usr/bin/env python

#
# Author: Steven Ludtke (sludtke@bcm.edu)
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
import EMAN2
import sys
import numpy
from emimageutil import ImgHistogram,EMParentWin
from weakref import WeakKeyDictionary
from pickle import dumps,loads
from PyQt4.QtGui import QImage
from PyQt4.QtCore import QTimer

from emglobjects import EMOpenGLFlagsAndTools

from emfloatingwidgets import EMGLRotaryWidget, EMGLView2D,EM3DWidget


class EMImageRotary(QtOpenGL.QGLWidget):
	"""A QT widget for rendering EMData objects. It can display stacks of 2D images
	in 'matrix' form on the display. The middle mouse button will bring up a
	control-panel. The QT event loop must be running for this object to function
	properly.
	"""
	allim=WeakKeyDictionary()
	def __init__(self, data=None,parent=None):
		self.image_rotary = None
		#self.initflag = True
		self.mmode = "drag"

		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True);
		#fmt.setDepthBuffer(True)
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		EMImageRotary.allim[self]=0
		
		
		self.image_rotary = EMImageRotaryCore(data,self)
		
		self.imagefilename = None
		
		self.fov = 20
		self.aspect = 1.0
		self.zNear = 1
		self.zFar = 5000
		
		self.animatables = []
		
		self.timer = QTimer()
		QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.timeout)
		self.timer.start(10)
		
	def setData(self,data):
		self.image_rotary.setData(data)
	
	def get_optimal_size(self):
		lr = self.image_rotary.rotary.get_suggested_lr_bt_nf()
		width = lr[1] - lr[0]
		height = lr[3] - lr[2]
		return [width+80,height+20]
	
	def timeout(self):
		
		if len(self.animatables) == 0: return
		
		for i,animatable in enumerate(self.animatables):
			if not animatable.animate(time.time()):
				# this could be dangerous
				self.animatables.pop(i)
		
		self.updateGL()
		
	def register_animatable(self,animatable):
		self.animatables.append(animatable)
	
	def setImageFileName(self,name):
		#print "set image file name",name
		self.imagefilename = name
		
	def getImageFileName(self):
		return self.imagefilename
	
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
		
		glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST)
		glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST)
		glHint(GL_TEXTURE_COMPRESSION_HINT, GL_NICEST)
		
	def paintGL(self):
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		
		if ( self.image_rotary == None ): return
		self.image_rotary.render()

	
	def resizeGL(self, width, height):
		if width <= 0 or height <= 0: return None
		GL.glViewport(0,0,width,height)
	
		GL.glMatrixMode(GL.GL_PROJECTION)
		GL.glLoadIdentity()
		self.aspect = float(width)/float(height)
		GLU.gluPerspective(self.fov,self.aspect,1,2000)
		#GL.glOrtho(0.0,width,0.0,height,-width,width)
		GL.glMatrixMode(GL.GL_MODELVIEW)
		GL.glLoadIdentity()
		
	def get_depth_for_height(self, height):
		# This function returns the width and height of the renderable 
		# area at the origin of the data volume
		depth = height/(2.0*tan(self.fov/2.0*pi/180.0))
	
		return depth
	def set_mmode(self,mode):
		self.mmode = mode
		self.image_rotary.set_mmode(mode)
	
	def mousePressEvent(self, event):
		self.image_rotary.mousePressEvent(event)
			
	def wheelEvent(self,event):
		self.image_rotary.wheelEvent(event)
	
	def mouseMoveEvent(self,event):
		self.image_rotary.mouseMoveEvent(event)

	def mouseReleaseEvent(self,event):
		self.image_rotary.mouseReleaseEvent(event)
	
	def keyPressEvent(self,event):
		if self.mmode == "app":
			self.emit(QtCore.SIGNAL("keypress"),event)

	def dropEvent(self,event):
		self.image_rotary.dropEvent(event)
		
	def closeEvent(self,event) :
		self.image_rotary.closeEvent(event)
		
	def dragEnterEvent(self,event):
		self.image_rotary.dragEnterEvent(event)

	def dropEvent(self,event):
		self.image_rotary.dropEvent(event)
	
	
	def set_shapes(self,shapes,shrink):
		self.image_rotary.set_shapes(shapes,shrink)
	
	def set_frozen(self,frozen):
		self.image_rotary.set_frozen(frozen)
	
class EMImageRotaryCore:

	#allim=WeakKeyDictionary()
	def __init__(self, data=None,parent=None):
		self.parent = parent
		self.data=None
		try: self.parent.setAcceptDrops(True)
		except:	pass

		self.initsizeflag = True
		if data:
			self.setData(data)
		
		self.rotary = EMGLRotaryWidget(self,-25,10,40,EMGLRotaryWidget.LEFT_ROTARY)
		self.widget = EM3DWidget(self,self.rotary)
		self.widget.set_draw_frame(False)
		
		#self.rotary.set_shapes([],1.01)
	def context(self):
		# asking for the OpenGL context from the parent
		return self.parent.context()
	
	def emit(self,signal,event,integer=None):
		if integer != None:
			self.parent.emit(signal,event,integer)
		else:
			self.parent.emit(signal,event)	
	def set_mmode(self,mode):
		self.mmode = mode
		self.rotary.set_mmode(mode)

	def set_frozen(self,frozen):
		self.rotary.set_frozen(frozen)

	def set_shapes(self,shapes,shrink):
		self.rotary.set_shapes(shapes,shrink)

	def register_animatable(self,animatable):
		self.parent.register_animatable(animatable)
		
	def width(self):
		return self.parent.width()
	
	def height(self):
		return self.parent.height()

	def getImageFileName(self):
		''' warning - could return none in some circumstances'''
		try: return self.parent.getImageFileName()
		except: return None
		
	def setData(self,data):
		if data == None or not isinstance(data,list) or len(data)==0:
			self.data = [] 
			return
		
		self.data = data
		
		self.rotary.clear_widgets()
		for d in self.data:
			w = EMGLView2D(self,d)
			self.rotary.add_widget(w)
			
		#self.showInspector()		# shows the correct inspector if already open
		#self.timer.start(25)
		
		# experimental for lst file writing
		for i,d in enumerate(data):
			d.set_attr("original_number",i)

	def updateGL(self):
		try: self.parent.updateGL()
		except: pass


	def render(self):
		if not self.data : return
		
		glLoadIdentity()
		
		lr = self.rotary.get_suggested_lr_bt_nf()
		#print lr
		GL.glEnable(GL.GL_DEPTH_TEST)
		GL.glEnable(GL.GL_LIGHTING)
		#print self.parent.get_depth_for_height(self.height())
		#lr = self.rotary.get_lr_bt_nf()

		GL.glPushMatrix()
		glTranslate(-(lr[1]+lr[0])/2.0,-(lr[3]+lr[2])/2.0,-self.parent.get_depth_for_height(abs(lr[3]-lr[2])))
		self.widget.paintGL()
		GL.glPopMatrix()
	
	def dragEnterEvent(self,event):
		pass

	
	def dropEvent(self,event):
		pass


	def mousePressEvent(self, event):
		self.widget.mousePressEvent(event)
		self.updateGL()
	
	def mouseMoveEvent(self, event):
		self.widget.mouseMoveEvent(event)
		self.updateGL()
		
	def mouseReleaseEvent(self, event):
		self.widget.mouseReleaseEvent(event)
		self.updateGL()
		
	def wheelEvent(self, event):
#		if event.delta() > 0:
#			self.setScale( self.scale * self.mag )
#		elif event.delta() < 0:
#			self.setScale(self.scale * self.invmag )
#		self.resizeEvent(self.parent.width(),self.parent.height())
#		# The self.scale variable is updated now, so just update with that
#		if self.inspector: self.inspector.setScale(self.scale)
		self.widget.wheelEvent(event)
		self.updateGL()
	def leaveEvent(self):
		pass


# This is just for testing, of course
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	GLUT.glutInit("")
	print "a"
	window = EMImageRotary()
	if len(sys.argv)==1 : window.setData([test_image(),test_image(1),test_image(2),test_image(3)]*4)
	else :
		a=EMData.read_images(sys.argv[1])
		window.setImageFileName(sys.argv[1])
		window.setData(a)
	window2=EMParentWin(window)
	window2.show()
	window2.resize(*window.get_optimal_size())
#	w2=QtGui.QWidget()
#	w2.resize(256,128)
	
#	w3=ValSlider(w2)
#	w3.resize(256,24)
#	w2.show()
	
	sys.exit(app.exec_())
