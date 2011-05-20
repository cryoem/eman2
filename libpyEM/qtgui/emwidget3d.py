#!/usr/bin/env python
#
# Author: John Flanagan (jfflanag@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine


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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

from EMAN2 import *
from OpenGL.GL import *
from OpenGL import GLU
from PyQt4 import QtCore, QtGui, QtOpenGL 
from PyQt4.QtCore import Qt
from emapplication import EMGLWidget
from emitem3d import EMItem3D
from libpyGLUtils2 import GLUtil
import math

# XPM format Cursors
zrotatecursor = [
    '15 14 2 1',
    'b c #00ff00',
    'c c None',
    'ccccccccccccccc',
    'ccccbbbbbbccbcc',
    'ccbbbbbbbbbbbbc',
    'cbbbcccccbbbbbc',
    'bbbcccccbbbbbbb',
    'bbcccccbbbbbbbc',
    'ccccccccccccccc',
    'ccccccccccccccc',
    'cbbbbbbbcccccbb',
    'bbbbbbbcccccbbb',
    'cbbbbbcccccbbbc',
    'cbbbbbbbbbbbbcc',
    'ccbccbbbbbbcccc',
    'ccccccccccccccc'
]

xyrotatecursor = [
    '14 13 2 1',
    'b c #00ff00',
    'c c None',
    'cccccccccccccc',
    'ccccbbbbbbcccc',
    'ccbbbbbbbbbbcc',
    'cbbbccccccbbbc',
    'bbbccccccccbbb',
    'bbbccccccccbbb',
    'bbbccccccccbbb',
    'bbbccccccccbbb',
    'bbbccccccccbbb',
    'cbbbccccccbbbc',
    'ccbbbbbbbbbbcc',
    'ccccbbbbbbcccc',
    'cccccccccccccc'
]

crosshairscursor = [
    '16 16 2 1',
    'b c #00ff00',
    'c c None',
    'cccccccbbcccccccc',
    'ccccccbbbbccccccc',
    'cccccbbbbbbcccccc',
    'cccccccbbcccccccc',
    'cccccccbbcccccccc',
    'ccbccccbbccccbccc',
    'cbbccccbbccccbbcc',
    'bbbbbbbbbbbbbbbbb',
    'bbbbbbbbbbbbbbbbb',
    'cbbccccbbccccbbcc',
    'ccbccccbbccccbccc',
    'cccccccbbcccccccc',
    'cccccccbbcccccccc',
    'cccccbbbbbbcccccc',
    'ccccccbbbbccccccc',
    'cccccccbbcccccccc'
]   
 
zhaircursor = [
    '16 16 2 1',
    'b c #00ff00',
    'c c None',
    'cccccccbbcccccccc',
    'ccccccbbbbccccccc',
    'cccccbbbbbbcccccc',
    'ccccbbcbbcbbccccc',
    'cccbbccbbccbbcccc',
    'cccccccbbcccccccc',
    'cccccccbbcccccccc',
    'cccccccbbcccccccc',
    'cccccccbbcccccccc',
    'cccccccbbcccccccc',
    'cccccccbbcccccccc',
    'cccbbccbbccbbcccc',
    'ccccbbcbbcbbccccc',
    'cccccbbbbbbcccccc',
    'ccccccbbbbccccccc',
    'cccccccbbcccccccc'
]   
scalecursor = [
    '16 16 2 1',
    'b c #00ff00',
    'c c None',
    'bbbbbbbbcccccccc',
    'bccccccccccccccc',
    'bccccccccccccccc',
    'bccccccccccccccc',
    'bccccccccccccccc',
    'bccccbbbbbbccccc',
    'bccccbccccbccccc',
    'bccccbccccbccccc',
    'cccccbccccbccccb',
    'cccccbccccbccccb',
    'cccccbbbbbbccccb',
    'cccccccccccbcccb',
    'ccccccccccccbccb',
    'cccccccccccccbcb',
    'ccccccccccccccbb',
    'ccccccccbbbbbbbb'
]   

class EMScene3DWidget(EMItem3D, EMGLWidget):
	"""
	Widget for rendering 3D objects. Uses a scne graph for rendering
	"""
	def __init__(self, parentwidget=None, SGactivenodeset=set(), scalestep=0.5):
		"""
		@param parent: The parent of the widget
		@param SGnodelist: a list enumerating all the SGnodes
		@param SGactivenodeset: a set enumerating the list of active nodes
		@param scalestep: The step to increment the object scaling
		"""
		EMItem3D.__init__(self, parent=None, transform=Transform())
		EMGLWidget.__init__(self,parentwidget)
		QtOpenGL.QGLFormat().setDoubleBuffer(True)
		self.camera = EMCamera(1.0, 10000.0)	# Default near,far, and zclip values
		#self.SGactivenodeset = SGactivenodeset			# A set of all active nodes (currently not used)
		self.scalestep = scalestep				# The scale factor stepsize
		self.zrotatecursor = QtGui.QCursor(QtGui.QPixmap(zrotatecursor),-1,-1)
		self.xyrotatecursor = QtGui.QCursor(QtGui.QPixmap(xyrotatecursor),-1,-1)
		self.crosshaircursor = QtGui.QCursor(QtGui.QPixmap(crosshairscursor),-1,-1)
		self.scalecursor = QtGui.QCursor(QtGui.QPixmap(scalecursor),-1,-1)
		self.zhaircursor = QtGui.QCursor(QtGui.QPixmap(zhaircursor),-1,-1)

	def initializeGL(self):
		glClearColor(0.0, 0.0, 0.0, 0.0)		# Default clear color is black
		glShadeModel(GL_SMOOTH)
		glEnable(GL_DEPTH_TEST)
		self.firstlight = EMLight(GL_LIGHT0)
		self.firstlight.enablelighting()

	def paintGL(self):
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		glColor3f(1.0, 1.0, 1.0)	# Default color is white
		#Call rendering
		self.render()		# SG nodes must have a render method
		glFlush()			# Finish rendering

	def resizeGL(self, width, height):
		self.camera.update(width, height)
	
	# Functions to control the list of nodes whose transforms will be updated
	#def activatenode(self, node):
	#	"""
	#	@param node: A SG node to activate
	#	Adds a node to the set, SGactivenodeset, enumerating the active nodes
	#	"""
	#	self.SGactivenodeset.add(node)
	#	node.isactive = True
	#	# When a node is actived we need to do some error checking to make sure that both parent and child node are in the set (otherwise thing will be moved multiple times)
	#	parentnode = node.parent
	#	while(parentnode):
	#		if parentnode.isactive:
	#			self.deactivatenode(parentnode)
	#		parentnode = parentnode.parent
	#	# Remove any children of the active nodelist.......
	#	
	#def deactivatenode(self, node):
	#	"""
	#	@param node: A SG node to deactivatenodeDeactivate the nodelist"
	#	"""
	#	if node.isactive:
	#		self.SGactivenodeset.remove(node)
	#		node.isactive = False
	#
			
	# Event subclassing
	def mousePressEvent(self, event):
		"""
		QT event handler. Records the coords when a mouse button is pressed and sets the cursor depending on what button(s) are pressed
		"""
		self.previous_x = event.x()
		self.previous_y = event.y()
		if event.buttons()&Qt.LeftButton:
			if  event.y() > 0.95*self.size().height():
				self.setCursor(self.zrotatecursor)
			else:
				self.setCursor(self.xyrotatecursor)
		if event.buttons()&Qt.MidButton:
			if event.modifiers()&Qt.ControlModifier:
				self.setCursor(self.zhaircursor)
			else:
				self.setCursor(self.crosshaircursor)
		if event.buttons()&Qt.RightButton:
			self.setCursor(self.scalecursor)		
			
	def mouseMoveEvent(self, event):
		"""
		Qt event handler. Scales the SG depending on what mouse button(s) are pressed when dragged
		"""
		dx = event.x() - self.previous_x
		dy = event.y() - self.previous_y
		if event.buttons()&Qt.LeftButton:
			magnitude = math.sqrt(dx*dx + dy*dy)
			#Check to see if the cursor is in the 'virtual slider pannel'
			if  event.y() > 0.95*self.size().height(): # The lowest 5% of the screen is reserved from the Z spin virtual slider
				self.setCursor(self.zrotatecursor)
				self.update_matrices([magnitude,0,0,-dx/magnitude], "rotate")
			else:
				self.setCursor(self.xyrotatecursor) 
				self.update_matrices([magnitude,-dy/magnitude,-dx/magnitude,0], "rotate")
			self.updateSG()
		if event.buttons()&Qt.MidButton:
			if event.modifiers()&Qt.ControlModifier:
				self.update_matrices([0,0,(dx+dy)], "translate")
			else:
				self.update_matrices([dx,-dy,0], "translate")
			self.updateSG()	
		if event.buttons()&Qt.RightButton:
			self.update_matrices([self.scalestep*0.1*(dx+dy)], "scale")
			self.setCursor(self.scalecursor)
			self.updateSG()
		self.previous_x =  event.x()
		self.previous_y =  event.y()
			
			
	def mouseReleaseEvent(self, event):
		"""
		Qt event handler. Returns the cursor to arrow unpn mouse button release
		"""
		self.setCursor(Qt.ArrowCursor)
	
	def wheelEvent(self, event):
		"""
		QT event handler. Scales the SG unpon wheel movement
		"""
		if event.orientation() & Qt.Vertical:
			if event.delta() > 0:
				self.update_matrices([self.scalestep], "scale")
			else:
				self.update_matrices([-self.scalestep], "scale")
			self.updateSG()
			
	def mouseDoubleClickEvent(self,event):
		print "Mouse Double Click Event"
	
	def keyPressEvent(self,event):
		print "Mouse KeyPress Event"
	
	def updateSG(self):
		"""
		Update the SG
		"""
		QtOpenGL.QGLWidget.updateGL(self)
		
	# Maybe add methods to control the lights

class EMLight:
	def __init__(self, light):
		"""
		@type light: GL_LIGHTX, where 0 =< X <= 8
		@param light: an OpenGL light
		The light properties are set to reasnonale defaults.
		"""
		self.light = light
		self.setambient(0.1, 0.1, 0.1, 1.0)		# Default ambient color is light grey
		self.setdiffuse(1.0, 1.0, 1.0, 1.0)		# Default diffuse color is white
		self.setspecualar(1.0, 1.0, 1.0, 1.0)		# Default specular color is white
		self.setposition(0.0, 0.0, 1.0, 0.0)		# Defulat position is 0, 0, 1.0 and light is directional (w=0)
		if not glIsEnabled(GL_LIGHTING):
			glEnable(GL_LIGHTING)

	def setambient(self, r, g, b, a):
		"""
		@param r: the red component of the ambient light
		@param g: the green component of the ambient light
		@param b: the blue component of the ambient light
		@param a: the alpha component of the ambient light
		Set the ambient light color
		"""
		self.colorambient = [r, g, b, a]
		glLightfv(self.light, GL_AMBIENT, self.colorambient)

	def setdiffuse(self, r, g, b, a):
		"""
		@param r: the red component of the diffuse and specular light
		@param g: the green component of the diffuse and specular light
		@param b: the blue component of the diffuse and specular light
		@param a: the alpha component of the diffuse and specular light
		Set the diffuse light color
		"""
		self.colordiffuse = [r, g, b, a]
		glLightfv(self.light, GL_DIFFUSE, self.colordiffuse)
		
	def setspecualar(self, r, g, b, a):
		"""
		@param r: the red component of the diffuse and specular light
		@param g: the green component of the diffuse and specular light
		@param b: the blue component of the diffuse and specular light
		@param a: the alpha component of the diffuse and specular light
		Set the specualr light color
		"""
		self.colorspecular = [r, g, b, a]
		glLightfv(self.light, GL_SPECULAR, self.colorspecular)

	def setposition(self, x, y, z, w):
		"""
		@param x: The x component of the light position
		@param y: The y component of the light position
		@param z: The z component of the light position
		@param w: The w component of the light position
		Set the light position, in gomogenious corrds
		"""
		self.position = [x, y, z, w]
		glLightfv(self.light, GL_POSITION, self.position)

	def enablelighting(self):
		"""
		Enables this light
		"""
		if not glIsEnabled(self.light):
			glEnable(self.light)

	def disablelighting(self):
		"""
		Disables this light
		"""
		if glIsEnabled(self.light):
			glDisable(self.light)
class EMCamera:
	"""Implmentation of the camera"""
	def __init__(self, near, far, zclip=-1000.0, useothro=True, fovy=60.0, boundingbox=50.0, screenfraction=0.5):
		"""
		@param fovy: The field of view angle
		@param near: The volume view near position
		@param far: The volume view far position
		@param zclip: The zclipping plane (basicaly how far back the camera is)
		@param useothro: Use orthographic projection
		@param boundingbox: The dimension of the bounding for the object to be rendered
		@param screenfraction: The fraction of the screen height to occupy
		"""
		self.far = far
		self.near = near
		if useothro:
			self.useortho(zclip)
		else:
			self.useprespective(boundingbox, screenfraction, fovy)

	def update(self, width, height):
		"""
		@param width: The width of the window in pixels
		@param height: The height of the window in pixels
		updates the camera and viewport after windowresize
		"""
		if self.useothro:
			glViewport(0,0,width,height)
			glMatrixMode(GL_PROJECTION)
			glLoadIdentity()
			glOrtho(-width/2, width/2, -height/2, height/2, self.near, self.far);
			glMatrixMode(GL_MODELVIEW)
			glLoadIdentity()
			glTranslate(0,0,self.zclip)
		else:
			# This may need some work to get it to behave
			glViewport(0,0,width,height)
			glMatrixMode(GL_PROJECTION)
			glLoadIdentity()
			GLU.gluPerspective(self.fovy, (float(width)/float(height)), self.near, self.far);
			glMatrixMode(GL_MODELVIEW)
			glLoadIdentity()
			glTranslate(0,0,-self.perspective_z) #How much to set the camera back depends on how big the object is
			
	def useprespective(self, boundingbox, screenfraction, fovy=60.0):
		""" 
		@param boundingbox: The dimension of the bounding for the object to be rendered
		@param screenfraction: The fraction of the screen height to occupy
		Changes projection matrix to perspective
		"""
		self.fovy = fovy
		self.perspective_z = (boundingbox/screenfraction)/(2*math.tan(math.radians(self.fovy/2)))  + boundingbox/2
		self.useothro = False
		

	def useortho(self, zclip):
		"""
		Changes projection matrix to orthographic
		"""
		self.useothro = True
		self.zclip = zclip

	def setclipfar(self, far):
		"""
		Set the far aspect of the viewing volume
		"""
		self.far = far

	def setclipnear(self, near):
		"""
		Set the near aspect of the viewing volume
		"""
		self.near = near

	def setfovy(self, fovy):
		"""
		Set the field of view angle aspect of the viewing volume
		"""
		self.fovy = fovy

	# Maybe other methods to control the camera
			

###################################### TEST CODE, THIS WILL NOT APPEAR IN THE WIDGET3D MODULE ##################################################
		
# All object that are rendered inherit from abstractSGnode and implement the render method
# In this example I use a cube, but any object can be drawn and so long as the object class inherits from abstractSGnode
class glCube(EMItem3D):
	def __init__(self, x, y, z, size):
		EMItem3D.__init__(self, parent=None, transform=Transform())
		self.xi = -size/2
		self.yi = -size/2
		self.zi = -size/2
		self.xf = size/2
		self.yf = size/2
		self.zf = size/2
		
		self.isactive = False
				
		# Matrix initialization
		self.nodematrix = Transform()
		self.parent = None
		
	def render_node(self):
		# So I can see the box
		
		glPushMatrix()
		GLUtil.glMultMatrix(self.nodematrix)
		
		# Material properties of the box
		glMaterialfv(GL_FRONT, GL_DIFFUSE, [0.5,0.5,0.5,1.0])
		glMaterialfv(GL_FRONT, GL_SPECULAR, [1.0,1.0,1.0,1.0])
		glMaterialf(GL_FRONT, GL_SHININESS, 25.0)
		glMaterialfv(GL_FRONT, GL_AMBIENT, [1.0, 1.0, 1.0, 1.0])
		# The box itself anlong with normal vectors
		glBegin(GL_QUADS)

		glNormal3f(self.xi, self.yi, self.zi + 1)
		glVertex3f(self.xi, self.yi, self.zi)
		glNormal3f(self.xf, self.yi, self.zi + 1)
		glVertex3f(self.xf, self.yi, self.zi)
		glNormal3f(self.xf, self.yf, self.zi + 1)
		glVertex3f(self.xf, self.yf, self.zi)
		glNormal3f(self.xi, self.yf, self.zi + 1)
		glVertex3f(self.xi, self.yf, self.zi)

		glNormal3f(self.xi - 1, self.yi, self.zi)
		glVertex3f(self.xi, self.yi, self.zi)
		glNormal3f(self.xi - 1, self.yi, self.zf)
		glVertex3f(self.xi, self.yi, self.zf)
		glNormal3f(self.xi - 1, self.yf, self.zf)
		glVertex3f(self.xi, self.yf, self.zf)
		glNormal3f(self.xi - 1, self.yf, self.zi)
		glVertex3f(self.xi, self.yf, self.zi)

		glNormal3f(self.xi, self.yi, self.zf - 1)
		glVertex3f(self.xi, self.yi, self.zf)
		glNormal3f(self.xi, self.yf, self.zf - 1)
		glVertex3f(self.xi, self.yf, self.zf)
		glNormal3f(self.xf, self.yf, self.zf - 1)
		glVertex3f(self.xf, self.yf, self.zf)
		glNormal3f(self.xf, self.yi, self.zf - 1)
		glVertex3f(self.xf, self.yi, self.zf)

		glNormal3f(self.xf + 1, self.yf, self.zf)
		glVertex3f(self.xf, self.yf, self.zf)
		glNormal3f(self.xf + 1, self.yf, self.zi)
		glVertex3f(self.xf, self.yf, self.zi)
		glNormal3f(self.xf + 1, self.yi, self.zi)
		glVertex3f(self.xf, self.yi, self.zi)
		glNormal3f(self.xf + 1, self.yi, self.zf)
		glVertex3f(self.xf, self.yi, self.zf)

		glNormal3f(self.xi, self.yf + 1, self.zi)
		glVertex3f(self.xi, self.yf, self.zi)
		glNormal3f(self.xi, self.yf + 1, self.zf)
		glVertex3f(self.xi, self.yf, self.zf)
		glNormal3f(self.xf, self.yf + 1, self.zf)
		glVertex3f(self.xf, self.yf, self.zf)
		glNormal3f(self.xf, self.yf + 1, self.zi)
		glVertex3f(self.xf, self.yf, self.zi)

		glNormal3f(self.xi, self.yi - 1, self.zi)
		glVertex3f(self.xi, self.yi, self.zi)
		glNormal3f(self.xi, self.yi - 1, self.zf)
		glVertex3f(self.xi, self.yi, self.zf)
		glNormal3f(self.xf, self.yi - 1, self.zf)
		glVertex3f(self.xf, self.yi, self.zf)
		glNormal3f(self.xf, self.yi - 1, self.zi)
		glVertex3f(self.xf, self.yi, self.zi)

		glEnd()
		glPopMatrix()
		
class GLdemo(QtGui.QWidget):
	def __init__(self):
		QtGui.QWidget.__init__(self)
		self.widget = EMScene3DWidget()
		self.widget.camera.useprespective(50, 0.5)
		self.cube1 = glCube(100.0, 100.0, -1000.0, 50.0)
		self.widget.add_child(self.cube1)    # Something to Render something..... (this could just as well be one of Ross's SGnodes)
		#self.widget.activatenode(cube1)
		self.cube2 = glCube(200.0, 200.0, -1000.0, 50.0)
		self.widget.add_child(self.cube2)
		#self.widget.activatenode(cube2)

		# QT stuff to display the widget
		self.cube1box = QtGui.QCheckBox("Cube1")
		self.cube2box = QtGui.QCheckBox("Cube2")
		self.rootbox = QtGui.QCheckBox("root")
		vbox = QtGui.QVBoxLayout()
		vbox.addWidget(self.widget)
		hbox = QtGui.QHBoxLayout()
		hbox.addWidget(self.cube1box)
		hbox.addWidget(self.cube2box)
		hbox.addWidget(self.rootbox)
		vbox.addLayout(hbox)
		self.setLayout(vbox)
		self.setGeometry(300, 300, 600, 600)
		self.setWindowTitle('BCM Chimera')
		
		self.connect(self.cube1box,QtCore.SIGNAL("stateChanged(int)"),self.on_updatecb1)
		self.connect(self.cube2box,QtCore.SIGNAL("stateChanged(int)"),self.on_updatecb2)
		self.connect(self.rootbox,QtCore.SIGNAL("stateChanged(int)"),self.on_updatecb3)
		
	def on_updatecb1(self):
		if self.cube1box.isChecked():
			self.cube1.is_selected = True
			self.rootbox.setCheckState(0)
			self.widget.is_selected = False
			
		else:
			self.cube1.is_selected = False
	
	def on_updatecb2(self):
		if self.cube2box.isChecked():
			self.cube2.is_selected = True
			self.rootbox.setCheckState(0)
			self.widget.is_selected = False
			
		else:
			self.cube2.is_selected = False
			
	def on_updatecb3(self):
		if self.rootbox.isChecked():
			self.widget.is_selected = True
			self.cube1box.setCheckState(0)
			self.cube2box.setCheckState(0)
			self.cube1.is_selected = False
			self.cube2.is_selected = False
		else:
			self.widget.is_selected = False
			
if __name__ == "__main__":
	import sys
	app = QtGui.QApplication(sys.argv)
	window = GLdemo()
	window.show()
	app.exec_()
