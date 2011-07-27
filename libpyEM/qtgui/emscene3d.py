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
from emitem3d import EMItem3D, EMItem3DInspector
from libpyGLUtils2 import GLUtil
from valslider import ValSlider, EMLightControls, CameraControls, EMSpinWidget, EMQTColorWidget
import math, weakref, os, pickle

# XPM format Cursors

visibleicon = [
    '16 12 3 1',
    'a c #0000ff',
    'b c #000000',
    'c c None',
    'cccccccccccccccc',
    'ccccccbbbbcccccc',
    'ccccbbbbbbbbcccc',
    'ccbbbccccccbbbcc',
    'cbbccccaaccccbbc',
    'cbccccaaaaccccbc',
    'cbccccaaaaccccbc',
    'cbbccccaaccccbbc',
    'ccbbbccccccbbbcc',
    'ccccbbbbbbbbcccc',
    'ccccccbbbbcccccc',
    'cccccccccccccccc'
]
    
invisibleicon = [
    '16 12 2 1',
    'b c #000000',
    'c c None',
    'cbbcccccccccbbcc',
    'ccbbcccccccbbccc',
    'cccbbcccccbbcccc',
    'ccccbbcccbbccccc',
    'cccccbbcbbcccccc',
    'ccccccbbbccccccc',
    'ccccccbbbccccccc',
    'cccccbbcbbcccccc',
    'ccccbbcccbbccccc',
    'cccbbcccccbbcccc',
    'ccbbcccccccbbccc',
    'cbbcccccccccbbcc'
]

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

selectorcursor = [
    '16 16 2 1',
    'b c #00ff00',
    'c c None',
    'cbbbbbbbbbcccccc',
    'bcccccccccbccccc',
    'cbbbbbbbcccbcccc',
    'cccbccccccccbccc',
    'ccccbbbbccccbccc',
    'cccbccccccccbccc',
    'ccccbbbbcccbcbcc',
    'cccbccccccbcccbc',
    'ccccbbbbbbcccccb',
    'cccccccbccbcccbc',
    'ccccccccbccccbcc',
    'cccccccccbccbccc',
    'ccccccccccbbcccc',
    'cccccccccccccccc',
    'cccccccccccccccc',
    'cccccccccccccccc'
]

class EMScene3D(EMItem3D, EMGLWidget):
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
		self.camera = EMCamera(1.0, 500.0)	# Default near,far, and zclip values
		self.clearcolor = [0.0, 0.0, 0.0, 0.0]	# Back ground color
		self.main_3d_inspector = None
		self.item_inspector = None				# Get the inspector GUI
		self.reset_camera = False			# Toogle flag to deterine if the clipping plane has changed and needs redrawing
		self.fuzzyselectionfactor = 2.0			# A fudge factor to determine the selection box 'thickness'
		#self.SGactivenodeset = SGactivenodeset			# A set of all active nodes (currently not used)
		self.scalestep = scalestep				# The scale factor stepsize
		self.toggle_render_selectedarea = False			# Don't render the selection box by default
		self.zrotatecursor = QtGui.QCursor(QtGui.QPixmap(zrotatecursor),-1,-1)
		self.xyrotatecursor = QtGui.QCursor(QtGui.QPixmap(xyrotatecursor),-1,-1)
		self.crosshaircursor = QtGui.QCursor(QtGui.QPixmap(crosshairscursor),-1,-1)
		self.scalecursor = QtGui.QCursor(QtGui.QPixmap(scalecursor),-1,-1)
		self.zhaircursor = QtGui.QCursor(QtGui.QPixmap(zhaircursor),-1,-1)
		self.selectorcursor = QtGui.QCursor(QtGui.QPixmap(selectorcursor),-1,-1)

	def getEvalString(self):
		"""
		Retrun a string that after eval can reinstatiate the object
		"""
		return "SG"
		
	def initializeGL(self):
		glClearColor(self.clearcolor[0], self.clearcolor[1], self.clearcolor[2], self.clearcolor[3])		# Default clear color is black
		glShadeModel(GL_SMOOTH)
		glEnable(GL_DEPTH_TEST)
		glEnable(GL_NORMALIZE)
		self.firstlight = EMLight(GL_LIGHT0)
		self.firstlight.enableLighting()
		if self.main_3d_inspector: self.main_3d_inspector.postGLInitialization()
        
	def paintGL(self):
		if self.reset_camera: self.camera.update()
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)		
		glColor3f(1.0, 1.0, 1.0)	# Default color is white
		#Call rendering
		self.renderSelectedArea() 	# Draw the selection box if needed
		self.render()			# SG nodes must have a render method
		glFlush()			# Finish rendering
		self.reset_camera = False
		
	def resizeGL(self, width, height):
		self.camera.update(width, height)
		if self.main_3d_inspector: self.main_3d_inspector.updateInspector()
	
	def getItemInspector(self):
		"""
		Return a Qt widget that controls the scene item
		"""	
		if not self.item_inspector: self.item_inspector = EMItem3DInspector("SG", self)
		return self.item_inspector
		
	def renderNode(self):
		pass
	
	def setInspector(self, inspector):
		"""
		Set the main 3d inspector
		"""
		self.main_3d_inspector = inspector
		
	def pickItem(self):
		"""
		Pick an item on the screen using openGL's selection mechanism
		"""
		viewport = glGetIntegerv(GL_VIEWPORT)
		glSelectBuffer(1024)	# The buffer size for selection
		glRenderMode(GL_SELECT)
		glInitNames()
		glMatrixMode(GL_PROJECTION)
		glPushMatrix()
		glLoadIdentity()
		
		# Find the selection box. Go from Volume view coords to viewport coords. sa = selection area
		x = self.sa_xi + self.camera.getWidth()/2
		y = self.camera.getHeight()/2 - self.sa_yi
		dx = self.fuzzyselectionfactor*math.fabs(self.sa_xi - self.sa_xf) # The 2x is a hack.....
		dy = self.fuzzyselectionfactor*math.fabs(self.sa_yi - self.sa_yf) # The 2x is a hack.....

		# Apply selection box
		GLU.gluPickMatrix(x, viewport[3] - y, dx, dy, viewport)
		self.camera.setProjectionMatrix()
		
		#drawstuff, but first we need to remove the influence of any previous xforms which ^$#*$ the selection
		glMatrixMode(GL_MODELVIEW)
		glPushMatrix()
		glLoadIdentity()
		self.camera.setCameraPosition(sfactor=1) # Factor of two to compensate for the samera already being set
		self.render()
		glPopMatrix()
		
		# Return to default state
		glMatrixMode(GL_PROJECTION)
		glPopMatrix()
		glMatrixMode(GL_MODELVIEW)
		records = glRenderMode(GL_RENDER)
		
		# process records
		self.processSelection(records)
	
	def selectArea(self, xi, xf, yi, yf):
		"""
		Set an area for selection. Need to switch bewteen viewport coords, where (0,0 is bottom left) to
		volume view coords where 0,0) is center of the screen.
		"""
		self.sa_xi = xi - self.camera.getWidth()/2
		self.sa_xf = xf - self.camera.getWidth()/2
		self.sa_yi = -yi + self.camera.getHeight()/2
		self.sa_yf = -yf + self.camera.getHeight()/2
		self.toggle_render_selectedarea = True
		
	def deselectArea(self):
		"""
		Turn off selectin box
		"""
		self.sa_xi = 0.0
		self.sa_xf = 0.0
		self.sa_yi = 0.0
		self.sa_yf = 0.0
		self.toggle_render_selectedarea = False
		
	def renderSelectedArea(self):
		"""
		Draw the selection box, box is always drawn orthographically
		"""
		if self.toggle_render_selectedarea: 
			glMatrixMode(GL_PROJECTION)
			glPushMatrix()
			glLoadIdentity()
			self.camera.setOrthoProjectionMatrix()
			glColor3f(0.0,1.0,0.0)
			glMaterialfv(GL_FRONT, GL_AMBIENT, [0.0,0.0,0.0,1.0])
			glMaterialfv(GL_FRONT, GL_DIFFUSE, [0.0,0.0,0.0,1.0])
			glMaterialfv(GL_FRONT, GL_SPECULAR, [0.0,0.0,0.0,1.0])
			glMaterialfv(GL_FRONT, GL_EMISSION, [0.0,1.0,0.0,1.0])
			glBegin(GL_LINE_LOOP)
			# set the box just in front of the cliping plane
			z = -self.camera.getClipNear() - self.camera.getZclip()
			glVertex3f(self.sa_xi, self.sa_yi, z)
			glVertex3f(self.sa_xi, self.sa_yf, z)
			glVertex3f(self.sa_xf, self.sa_yf, z)
			glVertex3f(self.sa_xf, self.sa_yi, z)
			glEnd()
			glPopMatrix()
			glMatrixMode(GL_MODELVIEW)
			glMaterialfv(GL_FRONT, GL_EMISSION, [0.0,0.0,0.0,1.0])
		
	def processSelection(self, records):
		"""
		Process the selection records
		"""
		# Remove old selection if not in append mode
		if not self.appendselection:
			for selected in self.getAllSelectedNodes():
				selected.setSelectedItem(False)
				# Inspector tree management
				if EMQTreeWidgetItem:
					selected.EMQTreeWidgetItem.setSelectionStateBox()
					self.main_3d_inspector.tree_widget.setCurrentItem(selected.EMQTreeWidgetItem)
		# Select the desired items	
		closestitem = None
		bestdistance = 1.0
		for record in records:
			selecteditem = EMItem3D.selection_idx_dict[record.names[len(record.names)-1]]()
			selecteditem.setSelectedItem(True)
			# Inspector tree management
			if self.EMQTreeWidgetItem:
				selecteditem.EMQTreeWidgetItem.setSelectionStateBox()
				self.main_3d_inspector.tree_widget.setCurrentItem(selecteditem.EMQTreeWidgetItem)
			try:
				self.main_3d_inspector.stacked_widget.setCurrentWidget(selecteditem.getItemInspector())
				self.main_3d_inspector.tree_widget.setCurrentItem(selecteditem.getItemInspector().treeitem)	# Hmmm... tak about tight coupling! It would be better to use getters
			except:
				pass
			#if record.near < bestdistance:
				#bestdistance = record.near
				#closestitem = record.names
			
	# Event subclassing
	def mousePressEvent(self, event):
		"""
		QT event handler. Records the coords when a mouse button is pressed and sets the cursor depending on what button(s) are pressed
		"""
		# The previous x,y records where the mouse was prior to mouse move, rest upon mouse move
		self.previous_x = event.x()
		self.previous_y = event.y()
		# The first x,y records where the mouse was first pressed
		self.first_x = self.previous_x
		self.first_y = self.previous_y
		if event.buttons()&Qt.LeftButton:
			if event.modifiers()&Qt.ControlModifier:
				self.setCursor(self.selectorcursor)
				self.appendselection = False
				if event.modifiers()&Qt.ShiftModifier:
					self.appendselection = True
			else:
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
			if event.modifiers()&Qt.ControlModifier:
				self.setCursor(self.selectorcursor)
				self.selectArea(self.first_x, event.x(), self.first_y, event.y())
			else:
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
		Qt event handler. Returns the cursor to arrow unpon mouse button release
		"""
		self.setCursor(Qt.ArrowCursor)
		# Select using the selection box
		if self.toggle_render_selectedarea:
			self.pickItem()
			self.deselectArea()
			self.updateSG()
		else:	# Pick without selection box, juwst click
			if event.modifiers()&Qt.ControlModifier:
				#self.selectArea(event.x(), event.x()+1, event.y(), event.y()+1)
				self.sa_xi = event.x() - self.camera.getWidth()/2
				self.sa_xf = event.x() + 1 - self.camera.getWidth()/2
				self.sa_yi = -event.y() + self.camera.getHeight()/2
				self.sa_yf = -event.y() + 1 + self.camera.getHeight()/2
				self.pickItem()
				self.updateSG()
				
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
	
	def saveSnapShot(self, filename, format="tiff"):
		image = self.grabFrameBuffer()
		image.save(filename, format)
		
	def getFuzzyFactor(self):
		"""Get the fuzzy selection factor """
		return self.fuzzyselectionfactor
		
	def setFuzzyFactor(self, factor):
		"""Set the fuzzy selection factor"""
		self.fuzzyselectionfactor = factor
		
	def setZclip(self, zclip):
		""" Set the Z clipping plane """
		self.camera.setZclip(zclip)
		self.reset_camera = True
		
	def setZslice(self):
		"""
		Get a Z slice to display in the camera widget, only works for orthoganal mode
		"""
		# Getting the Z slice will have problems when using perspective viewing
		self.setAutoBufferSwap(False)
		oldnear = self.camera.getClipNear() 
		oldfar = self.camera.getClipFar()
		# We want to see the full volume data rather than a clip, so move clipping planes to BIG
		# BIG is a bit differnt for perspective and orthgraphic volumes
		if self.camera.usingortho:
			self.camera.setClipNear(-1000)
			self.camera.setClipFar(1000)
		else:
			self.camera.setClipNear(1)
			self.camera.setClipFar(1000)
		self.reset_camera = True
		self.getTransform().rotate_origin(Transform({"type":"spin","Omega":90,"n1":0,"n2":1,"n3":0}))
		QtOpenGL.QGLWidget.updateGL(self)
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1)
		pixeldata = glReadPixels(1,1,self.camera.width,self.camera.height,GL_RGB,GL_UNSIGNED_BYTE)
		self.getTransform().rotate_origin(Transform({"type":"spin","Omega":-90,"n1":0,"n2":1,"n3":0}))
		# Then move back
		self.camera.setClipNear(oldnear)
		self.camera.setClipFar(oldfar)
		self.reset_camera = True	# Reset the camera when redering the real scene
		self.setAutoBufferSwap(True)
		self.pixels = []
		self.pixels.append(0)
		self.pixels.append(0)
		self.pixels.append(self.camera.width)
		self.pixels.append(self.camera.height)
		self.pixels.append(pixeldata)
	
	def cameraNeedsanUpdate(self):
		"""
		Tell the SG to restet the camera
		"""
		self.reset_camera = True
		
	def setClearColor(self, r, g, b):
		"""
		Set the background colorambient
		"""
		self.clearcolor = [r, g, b, 0.0]
		glClearColor(r, g, b, 0.0)
	
	def getClearColor(self):
		"""
		Return the clear color
		"""
		return self.clearcolor
		
	def updateSG(self):
		"""
		Update the SG
		"""
		QtOpenGL.QGLWidget.updateGL(self)
		if self.main_3d_inspector: self.main_3d_inspector.updateInspector()
	# Maybe add methods to control the lights

class EMLight:
	def __init__(self, light):
		"""
		@type light: GL_LIGHTX, where 0 =< X <= 8
		@param light: an OpenGL light
		The light properties are set to reasnonale defaults.
		"""
		self.light = light
		self.setAmbient(0.3, 0.3, 0.3, 1.0)		# Default ambient color is light grey
		self.setGlobalAmbient(0.0, 0.0, 0.0, 1.0)		# Default global ambient
		self.setDiffuse(1.0, 1.0, 1.0, 1.0)		# Default diffuse color is white
		self.setSpecualar(1.0, 1.0, 1.0, 1.0)		# Default specular color is white
		self.setPosition(0.0, 0.0, 1.0, 0.0)		# Defulat position is 0, 0, 1.0 and light is directional (w=0)
		if not glIsEnabled(GL_LIGHTING):
			glEnable(GL_LIGHTING)

	def setAmbient(self, r, g, b, a):
		"""
		@param r: the red component of the ambient light
		@param g: the green component of the ambient light
		@param b: the blue component of the ambient light
		@param a: the alpha component of the ambient light
		Set the ambient light color
		"""
		self.colorambient = [r, g, b, a]
		glLightfv(self.light, GL_AMBIENT, self.colorambient)
		

	def setGlobalAmbient(self, r, g, b, a):
		"""
		@param r: the red component of the ambient light
		@param g: the green component of the ambient light
		@param b: the blue component of the ambient light
		@param a: the alpha component of the ambient light
		Set the global ambient light color
		"""
		self.colorglobalambient = [r, g, b, a]
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT, self.colorglobalambient)	# We want no global ambient light
		
	def setDiffuse(self, r, g, b, a):
		"""
		@param r: the red component of the diffuse and specular light
		@param g: the green component of the diffuse and specular light
		@param b: the blue component of the diffuse and specular light
		@param a: the alpha component of the diffuse and specular light
		Set the diffuse light color
		"""
		self.colordiffuse = [r, g, b, a]
		glLightfv(self.light, GL_DIFFUSE, self.colordiffuse)
		
	def setSpecualar(self, r, g, b, a):
		"""
		@param r: the red component of the diffuse and specular light
		@param g: the green component of the diffuse and specular light
		@param b: the blue component of the diffuse and specular light
		@param a: the alpha component of the diffuse and specular light
		Set the specualr light color
		"""
		self.colorspecular = [r, g, b, a]
		glLightfv(self.light, GL_SPECULAR, self.colorspecular)

	def setPosition(self, x, y, z, w):
		"""
		@param x: The x component of the light position
		@param y: The y component of the light position
		@param z: The z component of the light position
		@param w: The w component of the light position
		Set the light position, in gomogenious corrds
		"""
		self.position = [x, y, z, w]
		glLightfv(self.light, GL_POSITION, self.position)

	def getAmbient(self):
		"""
		Return the ambient lighting
		"""
		return self.colorambient
	
	def getGlobalAmbient(self):
		"""
		Return the global ambient color
		"""
		return self.colorglobalambient
		
	def getDiffuse(self):
		"""
		Return the diffuse lighting
		"""
		return self.colordiffuse
		
	def getSpecular(self):
		"""
		Return the specular lighting
		"""
		returnself.colorspecular
	
	def getPosition(self):
		"""
		Return the light position
		"""
		return self.position
		
	def enableLighting(self):
		"""
		Enables this light
		"""
		if not glIsEnabled(self.light):
			glEnable(self.light)

	def disableLighting(self):
		"""
		Disables this light
		"""
		if glIsEnabled(self.light):
			glDisable(self.light)
class EMCamera:
	"""Implmentation of the camera"""
	def __init__(self, near, far, usingortho=True, fovy=60.0, boundingbox=50.0, screenfraction=0.5):
		"""
		@param fovy: The field of view angle
		@param near: The volume view near position
		@param far: The volume view far position
		@param usingortho: Use orthographic projection
		@param boundingbox: The dimension of the bounding for the object to be rendered
		@param screenfraction: The fraction of the screen height to occupy
		"""
		self.far = far
		self.near = near
		self.fovy = fovy
		zclip = (self.near-self.far)/2.0	# Puts things in the center of the viewing volume
		if usingortho:
			self.useOrtho(zclip)
		else:
			self.usePrespective(boundingbox, screenfraction, fovy)

	def update(self, width=None, height=None):
		"""
		@param width: The width of the window in pixels
		@param height: The height of the window in pixels
		updates the camera and viewport after windowresize
		"""
		if width: self.width = width
		if height: self.height = height
		if self.usingortho:
			glViewport(0,0,self.width,self.height)
			glMatrixMode(GL_PROJECTION)
			glLoadIdentity()
			self.setOrthoProjectionMatrix()
			glMatrixMode(GL_MODELVIEW)
			glLoadIdentity()
			self.setCameraPosition()
		else:
			# This may need some work to get it to behave
			glViewport(0,0,self.width,self.height)
			glMatrixMode(GL_PROJECTION)
			glLoadIdentity()
			self.setPerspectiveProjectionMatrix()
			glMatrixMode(GL_MODELVIEW)
			glLoadIdentity()
			self.setCameraPosition()
	
	def setCameraPosition(self, sfactor=1):
		"""
		Set the default camera position
		"""
		glTranslate(0,0,sfactor*self.getZclip())
		
	def setProjectionMatrix(self):
		"""
		Set the projection matrix
		"""
		if self.usingortho:
			self.setOrthoProjectionMatrix()
		else:
			self.setPerspectiveProjectionMatrix()
			
	def setOrthoProjectionMatrix(self):
		"""
		Set the orthographic projection matrix. Volume view origin (0,0) is center of screen
		"""
		glOrtho(-self.width/2, self.width/2, -self.height/2, self.height/2, self.near, self.far)
		
	def setPerspectiveProjectionMatrix(self):
		"""
		Set the perspective projection matrix. Volume view origin (0,0) is center of screen
		"""
		GLU.gluPerspective(self.fovy, (float(self.width)/float(self.height)), self.near, self.far)
			
	def usePrespective(self, boundingbox, screenfraction, fovy=60.0):
		""" 
		@param boundingbox: The dimension of the bounding for the object to be rendered
		@param screenfraction: The fraction of the screen height to occupy
		Changes projection matrix to perspective
		"""
		self.fovy = fovy
		self.perspective_z = -(boundingbox*2/screenfraction)/(2*math.tan(math.radians(self.fovy/2)))  + boundingbox
		self.usingortho = False
		
	def useOrtho(self, zclip):
		"""
		Changes projection matrix to orthographic
		"""
		self.usingortho = True
		self.zclip = zclip
		
	def getUseOrtho(self):
		"""
		Returns the projectrion state
		"""
		return self.usingortho
		
	def setClipFar(self, far):
		"""
		Set the far aspect of the viewing volume
		"""
		self.far = far

	def setClipNear(self, near):
		"""
		Set the near aspect of the viewing volume
		"""
		self.near = near

	def getClipNear(self):
		"""
		Get the near clipping plane
		"""
		return self.near
		
	def getClipFar(self):
		"""
		Get the far clipping plane
		"""
		return self.far
		
	def setFovy(self, fovy):
		"""
		Set the field of view angle aspect of the viewing volume
		"""
		self.fovy = fovy
	
	def getFovy(self):
		""" Return FOVY """
		return self.fovy
		
	def getHeight(self):
		""" Get the viewport height """
		return self.height
	
	def getWidth(self):
		""" Get the viewport width """
		return self.width
		
	def getZclip(self, disregardvv=False):
		""" Get the zclip """
		if self.usingortho or disregardvv:
			return self.zclip
		else:
			return self.perspective_z
			
	def setZclip(self, clip):
		""" Set the Z clip """
		if self.usingortho:
			self.zclip = clip
		else:
			self.perspective_z = clip
	# Maybe other methods to control the camera

###################################### Inspector Code #########################################################################################

class EMInspector3D(QtGui.QWidget):
	def __init__(self, scenegraph):
		"""
		The inspector for the 3D widget
		"""
		QtGui.QWidget.__init__(self)
		self.scenegraph = scenegraph
		self.mintreewidth = 200		# minimum width of the tree
		self.mincontrolwidth = 250
		
		vbox = QtGui.QVBoxLayout(self)
		
		self.inspectortab = QtGui.QTabWidget()
		self.inspectortab.addTab(self.getTreeWidget(), "Tree View")
		self.inspectortab.addTab(self.getLightsWidget(), "Lights")
		self.inspectortab.addTab(self.getCameraWidget(), "Camera")
		self.inspectortab.addTab(self.getUtilsWidget(), "Utils")

		vbox.addWidget(self.inspectortab)
		QtCore.QObject.connect(self.inspectortab, QtCore.SIGNAL("currentChanged(int)"), self._on_load_camera)
		
		self.setLayout(vbox)
		self.updateGeometry()
		
	def getTreeWidget(self):
		"""
		This returns the treeview-control panel widget
		"""
		widget = QtGui.QWidget()
		hbox = QtGui.QHBoxLayout(widget)
		treesplitter = QtGui.QFrame()
		treesplitter.setFrameShape(QtGui.QFrame.StyledPanel)
		treesplitter.setLayout(self._get_tree_layout(widget))
		treesplitter.setMinimumWidth(self.mintreewidth)
		hbox.addWidget(treesplitter)
		controlsplitter = QtGui.QFrame()
		controlsplitter.setFrameShape(QtGui.QFrame.StyledPanel)
		controlsplitter.setLayout(self._get_controler_layout(widget))
		controlsplitter.setMinimumWidth(self.mincontrolwidth)
		hbox.addWidget(controlsplitter)
		widget.setLayout(hbox)
		
		return widget
		
	def _get_tree_layout(self, parent):
		"""
		Returns the tree layout
		"""
		tvbox = QtGui.QVBoxLayout()
		self.tree_widget = EMQTreeWidget(parent)
		self.tree_widget.setHeaderLabel("Choose a item")
		tvbox.addWidget(self.tree_widget)
		self.tree_node_button_add = QtGui.QPushButton("Add Node")
		self.tree_node_button_remove = QtGui.QPushButton("Remove Node")
		tvbox.addWidget(self.tree_node_button_add)
		tvbox.addWidget(self.tree_node_button_remove)
		
		QtCore.QObject.connect(self.tree_widget, QtCore.SIGNAL("itemClicked(QTreeWidgetItem*,int)"), self._tree_widget_click)
		QtCore.QObject.connect(self.tree_widget, QtCore.SIGNAL("visibleItem(QTreeWidgetItem*)"), self._tree_widget_visible)
		QtCore.QObject.connect(self.tree_node_button_add, QtCore.SIGNAL("clicked()"), self._tree_widget_add)
		QtCore.QObject.connect(self.tree_node_button_remove, QtCore.SIGNAL("clicked()"), self._tree_widget_remove)
		
		return tvbox
	
	def _get_controler_layout(self, parent):
		"""
		Returns the control layout
		"""
		cvbox = QtGui.QVBoxLayout()
		self.stacked_widget = QtGui.QStackedWidget()
		cvbox.addWidget(self.stacked_widget)
		
		return cvbox
		
	def addTreeNode(self, name, item3d, parentnode=None, insertionindex=-1):
		"""
		Add a node (item3d) to the TreeWidget if not parent node, otherwise add a child to parent node
		We need to get a GUI for the treeitem. The treeitem and the GUI need know each other so they can talk
		The Treeitem also needs to know the node, so it can talk to the node.
		You can think of this as a three way conversation (the alterative it to use a mediator, but that is not worth it w/ only three players)
		"""
		tree_item = EMQTreeWidgetItem(QtCore.QStringList(name), item3d, parentnode)	# Make a QTreeItem widget, and let the TreeItem talk to the scenegraph node and its GUI
		item3d.setEMQTreeWidgetItem(tree_item)				# Reference to the EMQTreeWidgetItem
		item_inspector = item3d.getItemInspector()				# Get the node GUI controls 
		item_inspector.setInspector(self)					# Associate the item GUI with the inspector
		self.stacked_widget.addWidget(item_inspector)			# Add a widget to the stack
		# Set icon status
		tree_item.setSelectionStateBox()
		# Set parent if one exists	
		if not parentnode:
			self.tree_widget.insertTopLevelItem(0, tree_item)
		else:
			if insertionindex >= 0:
				parentnode.insertChild(insertionindex, tree_item)
			else:
				parentnode.addChild(tree_item)
		return tree_item
	
	def removeTreeNode(self, parentitem, childindex):
		# I am using the parent item rather than the item itself b/c the stupid widget has no , remove self function...
		# Remove both the QTreeWidgetItem and the widget from the WidgetStack, otherwise we'll get memory leaks 
		if parentitem.child(childindex).item3d(): self.stacked_widget.removeWidget(parentitem.child(childindex).item3d().getItemInspector())
		parentitem.takeChild(childindex)
	
	def getTreeWidgetChildrenData(self, item):
		"""
		Return a list of all the child items (actually a tree of sorts) along with all metadata. Used for saving a session
		item is the item you want to start at
		"""
		childlist = []
		dictionary = {"CONSTRUCTOR":item.item3d().getEvalString(),"NAME":str(item.name),"VISIBLE":item.item3d().isVisibleItem(),"SELECTED":item.item3d().isSelectedItem(),"NODETYPE":item.item3d().nodetype}
		if item.item3d().getTransform(): dictionary["TRANSFORMATION"] = item.item3d().getTransform().get_params("eman")
		# Process specific node types
		if item.item3d().nodetype == "ShapeNode" or item.item3d().nodetype == "DataChild":
			dictionary["COLOR"] = [item.item3d().ambient, item.item3d().diffuse, item.item3d().specular, item.item3d().shininess]
		if item.item3d().name == "Isosurface": dictionary["ISOPARS"] = [item.item3d().wire, item.item3d().cullbackfaces, item.item3d().isothr]
		
		# Process a SceneGraph if needed
		if item.item3d().getEvalString() == "SG":	
			dictionary["AMBIENTLIGHT"] = self.scenegraph.firstlight.getAmbient()
			dictionary["LIGHTPOSITION"] = self.lightwidget.getAngularPosition()
			dictionary["CAMERACLIP"] = [self.scenegraph.camera.getClipNear(), self.scenegraph.camera.getClipFar()]
			dictionary["CAMERAPM"] = self.scenegraph.camera.getUseOrtho()
			dictionary["CLEARCOLOR"] = self.scenegraph.getClearColor()
			dictionary["FUZZYFACTOR"] = self.scenegraph.getFuzzyFactor()
		# Now do the recursion to build tree
		childlist.append(dictionary)
		for i in xrange(item.childCount()):
			ichild = self.getTreeWidgetChildrenData(item.child(i))
			childlist.append(ichild)
		
		return childlist
		
	def _tree_widget_click(self, item, col):
		"""
		When a user clicks on the selection tree check box
		"""
		self.stacked_widget.setCurrentWidget(item.item3d().getItemInspector())
		item.setSelectionState(item.checkState(0))
		# This code is to prevent both decendents and childer from being selected....
		for ancestor in item.item3d().getSelectedAncestorNodes():
			if ancestor.EMQTreeWidgetItem:			# Not al ancestors are listed on the inspector tree (such as a data node)
				ancestor.EMQTreeWidgetItem.setSelectionState(False)
		for child in item.item3d().getAllSelectedNodes()[1:]: 	# Lop the node itself off
			child.EMQTreeWidgetItem.setSelectionState(False)
		self.updateSceneGraph()
		
	def _tree_widget_visible(self, item):
		"""
		When a user clicks on the visible icon
		"""
		item.toggleVisibleState()
		self.updateSceneGraph()
	
	def _tree_widget_add(self):
		"""
		when a user wants to add an item3d
		"""
		nodedialog = NodeDialog(self, self.tree_widget.currentItem())
		nodedialog.exec_()
		self.activateWindow()
		
	def _tree_widget_remove(self):
		"""
		When a use wnats to remove a node_name
		"""
		item = self.tree_widget.currentItem()
		if item.parent:
			self.removeTreeNode(item.parent(), item.parent().indexOfChild(item)) 
			item.parent().item3d().removeChild(item.item3d())
			self.updateSceneGraph()
		else:
			print "Error cannot remove root node!!"
		
	def getLightsWidget(self):
		"""
		Returns the lights control widget
		"""
		lwidget = QtGui.QWidget()
		lvbox = QtGui.QVBoxLayout()
		lightslabel = QtGui.QLabel("Lights", lwidget)
		lightslabel.setAlignment(QtCore.Qt.AlignCenter)
		lightslabel.setMaximumHeight(30.0)
		font = QtGui.QFont()
		font.setBold(True)
		lightslabel.setFont(font)
		lvbox.addWidget(lightslabel)
		self.lightwidget = EMLightControls(GL_LIGHT1)
		positionlabel = QtGui.QLabel("Position", lwidget)
		positionlabel.setMaximumHeight(20.0)
		positionlabel.setAlignment(QtCore.Qt.AlignCenter)
		valslidersplitter = QtGui.QFrame()
		valslidersplitter.setFrameShape(QtGui.QFrame.StyledPanel)
		valslidersplitter.setMaximumHeight(80)
		valvbox = QtGui.QVBoxLayout()
		self.hvalslider = ValSlider(lwidget,(0.0,360.0),"Horizontal")
		self.vvalslider = ValSlider(lwidget,(0.0,360.0),"Vertical")
		valvbox.addWidget(self.hvalslider)
		valvbox.addWidget(self.vvalslider)
		valslidersplitter.setLayout(valvbox)
		lvbox.addWidget(self.lightwidget)
		lvbox.addWidget(positionlabel)
		lvbox.addWidget(valslidersplitter)
		self.ambientlighting = ValSlider(lwidget,(0.0,1.0),"Ambient Lighting")
		self.ambientlighting.setMaximumHeight(30)
		lvbox.addWidget(self.ambientlighting)
		lwidget.setLayout(lvbox)
		
		QtCore.QObject.connect(self.lightwidget, QtCore.SIGNAL("lightPositionMoved"), self._light_position_moved)
		QtCore.QObject.connect(self.hvalslider,QtCore.SIGNAL("valueChanged"),self._on_light_slider)
		QtCore.QObject.connect(self.vvalslider,QtCore.SIGNAL("valueChanged"),self._on_light_slider)
		QtCore.QObject.connect(self.ambientlighting,QtCore.SIGNAL("valueChanged"),self._on_light_ambient)

		return lwidget
	
	def _light_position_moved(self, position): 
		angularposition = self.lightwidget.getAngularPosition()
		self.hvalslider.setValue(angularposition[0], quiet=1)
		self.vvalslider.setValue(angularposition[1], quiet=1)
		self.scenegraph.firstlight.setPosition(position[0], position[1], position[2], position[3])
		self.scenegraph.update()
	
	def _on_light_slider(self, value):
		self.lightwidget.setAngularPosition(self.hvalslider.getValue(), self.vvalslider.getValue())
		position = self.lightwidget.getPosition()
		self.scenegraph.firstlight.setPosition(position[0], position[1], position[2], position[3])
		self.scenegraph.update()
	
	def _on_light_ambient(self):
		ambient = self.ambientlighting.getValue()
		self.scenegraph.firstlight.setAmbient(ambient, ambient, ambient, 1.0)
		self.scenegraph.update()
	
	def getCameraWidget(self):
		"""
		Returns the camera control widget
		"""
		self.cameratab_open = False
		cwidget = QtGui.QWidget()
		grid = QtGui.QGridLayout()
		self.camerawidget = CameraControls(scenegraph=self.scenegraph)
		grid.addWidget(self.camerawidget, 0, 0, 1, 2)
		nlabel = QtGui.QLabel("Near clipping plane", cwidget)
		nlabel.setMaximumHeight(30.0)
		nlabel.setAlignment(QtCore.Qt.AlignCenter)
		self.near = EMSpinWidget(self.scenegraph.camera.getClipNear(), 1.0)
		self.near.setToolTip("In the Window above:\nClick 'n' drag, near the near clipping, to move near clipping plane")
		self.near.setMaximumHeight(40.0)
		grid.addWidget(nlabel, 1, 0)
		grid.addWidget(self.near ,1, 1)
		flabel = QtGui.QLabel("Far clipping plane", cwidget)
		flabel.setMaximumHeight(30.0)
		flabel.setAlignment(QtCore.Qt.AlignCenter)
		self.far = EMSpinWidget(self.scenegraph.camera.getClipFar(), 1.0)
		self.far.setToolTip("In the Window above:\nClick 'n' drag, near the far clipping, to move far clipping plane")
		self.far.setMaximumHeight(40.0)
		grid.addWidget(flabel, 2, 0)
		grid.addWidget(self.far, 2, 1)
		frame = QtGui.QFrame()
		frame.setMaximumHeight(40.0)
		frame.setFrameShape(QtGui.QFrame.StyledPanel)
		hbox = QtGui.QHBoxLayout()
		vvlabel = QtGui.QLabel("Viewing Volume")
		self.orthoradio = QtGui.QRadioButton("Orthographic")
		self.perspectiveradio = QtGui.QRadioButton("Perspective")
		hbox.addWidget(vvlabel)
		hbox.addWidget(self.orthoradio)
		hbox.addWidget(self.perspectiveradio)
		frame.setLayout(hbox)
		grid.addWidget(frame, 3, 0,1,2)
		cwidget.setLayout(grid)

		QtCore.QObject.connect(self.near,QtCore.SIGNAL("valueChanged(int)"),self._on_near)
		QtCore.QObject.connect(self.far,QtCore.SIGNAL("valueChanged(int)"),self._on_far)
		QtCore.QObject.connect(self.camerawidget,QtCore.SIGNAL("nearMoved(int)"),self._on_near_move)
		QtCore.QObject.connect(self.camerawidget,QtCore.SIGNAL("farMoved(int)"),self._on_far_move)
		QtCore.QObject.connect(self.orthoradio,QtCore.SIGNAL("clicked()"),self._on_radio_click)
		QtCore.QObject.connect(self.perspectiveradio,QtCore.SIGNAL("clicked()"),self._on_radio_click)
		
		return cwidget
	
	def _on_near(self, value):
		if not self.scenegraph.camera.usingortho and value <= 0:
			return
		if self.scenegraph.camera.getClipFar() > value:
			self.camerawidget.updateWidget()
			self.scenegraph.camera.setClipNear(value)
			self.scenegraph.updateSG()
	
	def _on_far(self, value):
		if value > self.scenegraph.camera.getClipNear():
			self.camerawidget.updateWidget()
			self.scenegraph.camera.setClipFar(value)
			self.scenegraph.updateSG()
		
	def _on_near_move(self, movement):
		value = self.scenegraph.camera.getClipNear() + movement
		if not self.scenegraph.camera.usingortho and value <= 0:
			return
		if self.scenegraph.camera.getClipFar() > value:
			self.near.setValue(value, quiet=1)
			self.scenegraph.camera.setClipNear(value)
			self.scenegraph.updateSG()
			self.camerawidget.updateWidget()
		
	def _on_far_move(self, movement):
		value = self.scenegraph.camera.getClipFar() + movement
		if value > self.scenegraph.camera.getClipNear():
			self.far.setValue(value, quiet=1)
			self.scenegraph.camera.setClipFar(value)
			self.scenegraph.updateSG()
			self.camerawidget.updateWidget()
		
	def _on_load_camera(self, idx):
		"""
		Load the SG Z slice
		"""
		if idx == 2: # '2' is the camera tab
			self.scenegraph.setZslice()
			self.cameratab_open = True
		else:
			self.cameratab_open = False
			
	def _get_vv_state(self):
		"""
		Get the viewing volume of the camera
		"""
		if self.scenegraph.camera.usingortho: 
			self.orthoradio.setChecked(True)
		else:
			self.perspectiveradio.setChecked(True)
	
	def _on_radio_click(self):
		"""
		set the viewing volume. objectsize is the size of the bounding box of the largest displayed object
		"""
		objectsize = 50
		if self.orthoradio.isChecked(): self.scenegraph.camera.useOrtho(self.scenegraph.camera.getZclip(disregardvv=True))
		if self.perspectiveradio.isChecked(): self.scenegraph.camera.usePrespective(objectsize, 0.25, 60.0)
		self.scenegraph.updateSG()
		
	def getUtilsWidget(self):
		"""
		Return the utilites widget
		"""
		uwidget = QtGui.QWidget()
		uvbox = QtGui.QVBoxLayout()
		frame = QtGui.QFrame()
		frame.setFrameShape(QtGui.QFrame.StyledPanel)
		fvbox = QtGui.QHBoxLayout()
		backgroundcolor_label = QtGui.QLabel("Background Color", frame)
		self.backgroundcolor = EMQTColorWidget(parent=frame)
		self.backgroundcolor.setColor(QtGui.QColor(255*self.scenegraph.clearcolor[0],255*self.scenegraph.clearcolor[1],255*self.scenegraph.clearcolor[2]))
		fvbox.addWidget(backgroundcolor_label)
		fvbox.addWidget(self.backgroundcolor)
		fvbox.setAlignment(QtCore.Qt.AlignCenter)
		frame.setLayout(fvbox)
		uvbox.addWidget(frame)
		self.fuzzy_slider = ValSlider(uwidget, (0.2, 5.0), "Fuzzy Selection factor")
		self.fuzzy_slider.setValue(self.scenegraph.getFuzzyFactor(), quiet=1)
		self.opensession_button = QtGui.QPushButton("Open Session")
		self.savesession_button = QtGui.QPushButton("Save Session")
		self.savebutton = QtGui.QPushButton("Save Image Snapshot")
		uvbox.addWidget(self.fuzzy_slider)
		uvbox.addWidget(self.opensession_button)
		uvbox.addWidget(self.savesession_button)
		uvbox.addWidget(self.savebutton)
		uwidget.setLayout(uvbox)
		
		QtCore.QObject.connect(self.fuzzy_slider, QtCore.SIGNAL("valueChanged"),self._on_fuzzy_slider)
		QtCore.QObject.connect(self.backgroundcolor,QtCore.SIGNAL("newcolor(QColor)"),self._on_bg_color)
		QtCore.QObject.connect(self.savebutton, QtCore.SIGNAL("clicked()"),self._on_save)
		QtCore.QObject.connect(self.savesession_button, QtCore.SIGNAL("clicked()"),self._on_save_session)
		QtCore.QObject.connect(self.opensession_button, QtCore.SIGNAL("clicked()"),self._on_open_session)
		
		return uwidget
	
	def _process_session_load(self, line):
		if type([]) == type(line):
			if len(line) > 1:
				for cline in line:
					self._process_session_load(cline)
				self.parentnodestack.pop()	
			else:
				# These nodes are leaves
				item3dobject = eval(line[0]["CONSTRUCTOR"])
				if line[0]["NODETYPE"] == "DataChild": # Data child need to have a parent set
					item3dobject.setParent(self.parentnodestack[-1:][0].item3d())
				self._constructitem3d(line[0], item3dobject)

		else:
			# These nodes have children to process
			if line["CONSTRUCTOR"] != "SG": 
				print line["CONSTRUCTOR"]
				item3dobject = eval(line["CONSTRUCTOR"])
				widgetitem = self._constructitem3d(line, item3dobject)
				self.parentnodestack.append(widgetitem)
			else:
				# This is the SG, so treat special
				self.scenegraph.setTransform(Transform(line["TRANSFORMATION"]))
				self.scenegraph.setVisibleItem(line["VISIBLE"])
				self.scenegraph.setSelectedItem(line["SELECTED"])
				# Set up the light
				self.scenegraph.firstlight.setAmbient(line["AMBIENTLIGHT"][0],line["AMBIENTLIGHT"][1],line["AMBIENTLIGHT"][2],line["AMBIENTLIGHT"][3])
				al = (line["AMBIENTLIGHT"][0] + line["AMBIENTLIGHT"][1] + line["AMBIENTLIGHT"][2])/3
				self.ambientlighting.setValue(al,quiet=1)
				self.lightwidget.setAngularPosition(line["LIGHTPOSITION"][0],line["LIGHTPOSITION"][1], quiet=False)
				p = self.lightwidget.getPosition()
				self.scenegraph.firstlight.setPosition(p[0],p[1],p[2],p[3])
				# Set up the Camera
				self.scenegraph.camera.setClipNear(line["CAMERACLIP"][0])
				self.scenegraph.camera.setClipFar(line["CAMERACLIP"][1])
				self.near.setValue(line["CAMERACLIP"][0], quiet=1)
				self.far.setValue(line["CAMERACLIP"][1], quiet=1)
				if line["CAMERAPM"]: 
					self.orthoradio.setChecked(True)
					self.scenegraph.camera.useOrtho(self.scenegraph.camera.getZclip(disregardvv=True))
				else:
					objectsize = 50
					self.perspectiveradio.setChecked(True)
					self.scenegraph.camera.usePrespective(objectsize, 0.25, 60.0)
				self.scenegraph.cameraNeedsanUpdate()
				# Set up the utils
				self.scenegraph.setClearColor(line["CLEARCOLOR"][0],line["CLEARCOLOR"][1],line["CLEARCOLOR"][2])
				self.backgroundcolor.setColor(QtGui.QColor(255*line["CLEARCOLOR"][0],255*line["CLEARCOLOR"][1],255*line["CLEARCOLOR"][2]))
				self.scenegraph.setFuzzyFactor(line["FUZZYFACTOR"])
				self.fuzzy_slider.setValue(line["FUZZYFACTOR"])
			
	def _constructitem3d(self, datadict, item3dobject):
		"""
		Helper function for _process_session_load
		"""
		if datadict.has_key("TRANSFORMATION"):
			item3dobject.setTransform(Transform(datadict["TRANSFORMATION"]))
		if datadict.has_key("COLOR"):
			item3dobject.setAmbientColor(datadict["COLOR"][0][0], datadict["COLOR"][0][1], datadict["COLOR"][0][2], datadict["COLOR"][0][3])
			item3dobject.setDiffuseColor(datadict["COLOR"][1][0], datadict["COLOR"][1][1], datadict["COLOR"][1][2], datadict["COLOR"][1][3])
			item3dobject.setSpecularColor(datadict["COLOR"][2][0], datadict["COLOR"][2][1], datadict["COLOR"][2][2], datadict["COLOR"][2][3])
			item3dobject.setShininess(datadict["COLOR"][3])
		if datadict.has_key("ISOPARS"):
			item3dobject.wire = datadict["ISOPARS"][0]
			item3dobject.cullbackfaces = datadict["ISOPARS"][0]
			item3dobject.cullbackfaces = datadict["ISOPARS"][1]
			item3dobject.setThreshold(datadict["ISOPARS"][2])
			
		item3dobject.setVisibleItem(datadict["VISIBLE"])
		item3dobject.setSelectedItem(datadict["SELECTED"])
		self.parentnodestack[-1:][0].item3d().addChild(item3dobject)
		return self.addTreeNode(datadict["NAME"], item3dobject, self.parentnodestack[-1:][0])
				
	def _on_open_session(self):
		"""
		Open a session
		"""
		# Open the file
		filename = QtGui.QFileDialog.getOpenFileName(self, 'Open Session', os.getcwd(), "*.eman")
		self.loadSession(filename)
	
	def loadSession(self, filename):
		"""
		Load a session
		"""
		rfile = open(filename, 'rb')
		try:
			tree = pickle.load(rfile)
		except:
			print "ERROR!!! Couldn't load the session file"
			rfile.close()
			return	
		rfile.close()
		# Clear the current tree
		self.tree_widget.topLevelItem(0).removeAllChildren(self)
		# Load the new data
		self.parentnodestack = [self.tree_widget.topLevelItem(0)]
		self._process_session_load(tree)
		self.updateSceneGraph()
		
	def _on_save_session(self):
		"""
		Return a list of all the child items (actually a tree of sorts)
		"""
		filename = QtGui.QFileDialog.getSaveFileName(self, 'Save Session', os.getcwd(), "*.eman")
		if filename: # if we cancel
			wfile = open(filename, 'wb')
			childrentree = self.getTreeWidgetChildrenData(self.tree_widget.topLevelItem(0))
			pickle.dump(childrentree, wfile)
			wfile.close()
		
	def _on_save(self):
		"""
		Save a snapshot of the scene
		"""
		filename = QtGui.QFileDialog.getSaveFileName(self, 'Save Image', os.getcwd(), "*.tiff")
		if filename: # if we cancel
			self.scenegraph.saveSnapShot(filename)
			print "Saved %s to disk"%os.path.basename(str(filename))
		
	def _on_fuzzy_slider(self, value):
		"""
		Set the fuzzy selection coeff
		"""
		self.scenegraph.setFuzzyFactor(value)
	
	def _on_bg_color(self, color):
		rgb = color.getRgb()
		self.scenegraph.setClearColor(float(rgb[0])/255.0, float(rgb[1])/255.0, float(rgb[2])/255.0)
		self.updateSceneGraph()
		# No, this is not  a bug.... It is a hack b/c sometimes openGL does not change its clear color the first time around!!!!!!
		self.scenegraph.setClearColor(float(rgb[0])/255.0, float(rgb[1])/255.0, float(rgb[2])/255.0)
		self.updateSceneGraph()
	
	def postGLInitialization(self):
		"""
		We run this method after openGL is intialized so that we can get GL params.... Otherwise we get CRAP!!!
		"""
		ambient = self.scenegraph.firstlight.getAmbient()
		avgambient = (ambient[0] + ambient[1] + ambient[2])/3
		self.ambientlighting.setValue(avgambient)
		self._get_vv_state()
		
	def updateInspector(self):
		"""
		Update Inspector, for now this only pertains to the camera widget
		"""
		if self.cameratab_open:
			self.scenegraph.setZslice()
			self.camerawidget.updateWidget()
			
	def updateSceneGraph(self):
		""" 
		Updates SG, in the near future this will be improved to allow for slow operations
		"""
		self.scenegraph.updateSG()

class EMQTreeWidget(QtGui.QTreeWidget):
	"""
	Subclassing the QTreeWidget to enable is_visible toggling
	"""
	def __init__(self, parent=None):
		QtGui.QTreeWidget.__init__(self, parent)
			
	def mousePressEvent(self, e):
		QtGui.QTreeWidget.mousePressEvent(self, e)
		if e.button()==Qt.RightButton:
			self.emit(QtCore.SIGNAL("visibleItem(QTreeWidgetItem*)"), self.currentItem())
			
class EMQTreeWidgetItem(QtGui.QTreeWidgetItem):
	"""
	Subclass of QTreeWidgetItem
	adds functionality
	"""
	def __init__(self, qstring, item3d, parentnode):
		QtGui.QTreeWidgetItem.__init__(self, qstring)
		self.name = qstring.join('')
		self.item3d = weakref.ref(item3d)
		if parentnode: self.parent = weakref.ref(parentnode)
		else: self.parent = None
		self.setCheckState(0, QtCore.Qt.Unchecked)
		self.visible = QtGui.QIcon(QtGui.QPixmap(visibleicon))
		self.invisible = QtGui.QIcon(QtGui.QPixmap(invisibleicon))
		self.getVisibleState()
		self.setToolTip(0, 'Click on the checkbox to select\nRight click to toogle visible')
	
	def setSelectionState(self, state):
		""" 
		Toogle selection state on and off
		"""
		if state == QtCore.Qt.Checked:
			self.item3d().setSelectedItem(True)
		else:
			self.item3d().setSelectedItem(False)
		self.setSelectionStateBox() # set state of TreeItemwidget
		
	def toggleVisibleState(self):
		self.item3d().setVisibleItem(not self.item3d().isVisibleItem())
		self.getVisibleState()
		
	def getVisibleState(self):
		"""
		Toogle the visble state
		"""
		if self.item3d().isVisibleItem():
			self.setIcon(0, self.visible)
		else:
			self.setIcon(0, self.invisible)
	
	def setSelectionStateBox(self):
		"""
		Set the selection state icon
		"""
		if self.item3d().isSelectedItem():
			self.setCheckState(0, QtCore.Qt.Checked)
		else:
			self.setCheckState(0, QtCore.Qt.Unchecked)
	
	def removeAllChildren(self, inspector):
		"""
		Remove all children from the SG
		"""
		for i in xrange(self.childCount()):
			self.child(0).removeAllChildren(inspector)
			self.item3d().removeChild(self.child(0).item3d())
			inspector.removeTreeNode(self, 0) 
			
class NodeDialog(QtGui.QDialog):
	"""
	Generate a dialog to add or remove node. If reome is chosen 'item' node is removed
	If add node is chosen, a node is inserted just below this node.
	"""
	def __init__(self, inspector, item):
		QtGui.QDialog.__init__(self)
		self.item = item
		self.inspector = weakref.ref(inspector)
		self.setWindowTitle('Node Controler')
		#grid = QtGui.QGridLayout()
		vbox = QtGui.QVBoxLayout(self)
		# Stuff within the frame
		frame = QtGui.QFrame()
		frame.setFrameStyle(QtGui.QFrame.StyledPanel)
		fvbox = QtGui.QVBoxLayout(frame)
		label = QtGui.QLabel("Node Type to add")
		self.node_type_combo = QtGui.QComboBox() 
		self.node_stacked_widget = QtGui.QStackedWidget()
		self.node_stacked_widget.setFrameStyle(QtGui.QFrame.StyledPanel)
		self.addnode_button = QtGui.QPushButton("Add Node")
		fvbox.addWidget(label)
		fvbox.addWidget(self.node_type_combo)
		fvbox.addWidget(self.node_stacked_widget)
		fvbox.addWidget(self.addnode_button)
		frame.setLayout(fvbox)
		# vbox widgets
		vbox.addWidget(frame)
		self.cancel_button = QtGui.QPushButton("Cancel")
		vbox.addWidget(self.cancel_button)
		self.setLayout(vbox)
		
		# Populate node types		
		self.node_type_combo.addItem("Cube")
		self.node_stacked_widget.addWidget(self.getCubeWidget())
		self.node_type_combo.addItem("Sphere")
		self.node_stacked_widget.addWidget(self.getSphereWidget())
		self.node_type_combo.addItem("Cylinder")
		self.node_stacked_widget.addWidget(self.getCylinderWidget())
		self.node_type_combo.addItem("Data")
		self.node_stacked_widget.addWidget(self.getDataWidget())
		if self.item.item3d().name == "Data":
			self.node_type_combo.addItem("Isosurface")
			self.node_stacked_widget.addWidget(self.getIsosurfaceWidget())
		
		self.connect(self.addnode_button, QtCore.SIGNAL('clicked()'), self._on_add_node)
		self.connect(self.cancel_button, QtCore.SIGNAL('clicked()'), self._on_cancel)
		self.connect(self.node_type_combo, QtCore.SIGNAL("activated(int)"), self._node_combobox_changed)
	
	def getCubeWidget(self):
		"""
		Return a cube control widget for the stacked_widget
		"""
		cubewidget = QtGui.QWidget()
		grid = QtGui.QGridLayout()
		cube_dim_label = QtGui.QLabel("Cube Dimension")
		self.cube_dim = QtGui.QLineEdit()
		node_name_label = QtGui.QLabel("Cube Name")
		self.node_name_cube = QtGui.QLineEdit()
		grid.addWidget(cube_dim_label, 0, 0)
		grid.addWidget(self.cube_dim, 0, 1)
		grid.addWidget(node_name_label , 1, 0)
		grid.addWidget(self.node_name_cube, 1, 1)
		cubewidget.setLayout(grid)
		return cubewidget
		
	def getSphereWidget(self):
		"""
		Return a sphere control widget for the stacked_widget
		"""
		spherewidget = QtGui.QWidget()
		grid = QtGui.QGridLayout()
		sphere_dim_label = QtGui.QLabel("Sphere Dimension")
		self.sphere_dim = QtGui.QLineEdit()
		node_name_label = QtGui.QLabel("Sphere Name")
		self.node_name_sphere = QtGui.QLineEdit()
		grid.addWidget(sphere_dim_label, 0, 0)
		grid.addWidget(self.sphere_dim, 0, 1)
		grid.addWidget(node_name_label , 1, 0)
		grid.addWidget(self.node_name_sphere, 1, 1)
		spherewidget.setLayout(grid)
		return spherewidget
		
	def getCylinderWidget(self):
		"""
		Return a cylider control widget for the stacked_widget
		"""
		cyliderwidget = QtGui.QWidget()
		grid = QtGui.QGridLayout()
		cylider_radius_label = QtGui.QLabel("Cylider Radius")
		self.cylider_radius = QtGui.QLineEdit()
		grid.addWidget(cylider_radius_label, 0, 0)
		grid.addWidget(self.cylider_radius, 0, 1)
		cylider_height_label = QtGui.QLabel("Cylider Height")
		self.cylider_height = QtGui.QLineEdit()
		node_name_label = QtGui.QLabel("Cylider Name")
		self.node_name_cylinder = QtGui.QLineEdit()
		grid.addWidget(cylider_height_label, 1, 0)
		grid.addWidget(self.cylider_height, 1, 1)
		grid.addWidget(node_name_label , 2, 0)
		grid.addWidget(self.node_name_cylinder, 2, 1)
		cyliderwidget.setLayout(grid)
		return cyliderwidget
	
	def getDataWidget(self):
		"""
		Get Data Widget
		"""
		datawidget = QtGui.QWidget()
		grid = QtGui.QGridLayout()
		node_name_data_label = QtGui.QLabel("Data Label")
		self.node_name_data = QtGui.QLineEdit()
		data_path_label = QtGui.QLabel("Data Path")
		self.data_path = QtGui.QLineEdit()
		self.browse_button = QtGui.QPushButton("Browse")
		grid.addWidget(node_name_data_label, 0, 0)
		grid.addWidget(self.node_name_data, 0, 1)
		grid.addWidget(data_path_label, 1, 0)
		grid.addWidget(self.data_path, 1, 1)
		grid.addWidget(self.browse_button, 2, 0, 1, 2)
		datawidget.setLayout(grid)
		
		self.connect(self.browse_button, QtCore.SIGNAL('clicked()'), self._on_browse)
		
		return datawidget
	
	def getIsosurfaceWidget(self):
		"""
		Get Isosurface Widget
		"""
		isosurfacewidget = QtGui.QWidget()
		grid = QtGui.QGridLayout()
		node_name_data_label = QtGui.QLabel("Isosurface Name")
		self.node_name_data = QtGui.QLineEdit()
		grid.addWidget(node_name_data_label, 0, 0)
		grid.addWidget(self.node_name_data, 0, 1)
		isosurfacewidget.setLayout(grid)
		
		return isosurfacewidget
		
	def _on_browse(self):
		filename = QtGui.QFileDialog.getOpenFileName(self, 'Get file', os.getcwd())
		self.data_path.setText(filename)
		name = os.path.basename(str(filename))
		self.node_name_data.setText(name)
	
	def _on_add_node(self):
		if not self.item.parent:
			print "Cannot add another top level node!!!"
			self.done(0)
			return
		insertion_node = None
		node_name = "default"
		itemparentnode = self.item.parent()
		# Cube
		if self.node_type_combo.currentText() == "Cube":
			insertion_node = EMCube(float(self.cube_dim.text()))
			if self.node_name_cube.text() != "": node_name = self.node_name_cube.text()
		# Sphere
		if self.node_type_combo.currentText() == "Sphere":
			insertion_node = EMSphere(float(self.sphere_dim.text()))
			if self.node_name_sphere.text() != "": node_name = self.node_name_sphere.text()
		# Cylinder
		if self.node_type_combo.currentText() == "Cylinder":
			insertion_node = EMCylinder(float(self.cylider_radius.text()), float(self.cylider_height.text()))
			if self.node_name_cylinder.text() != "": node_name = self.node_name_cylinder.text()
		# Data
		if self.node_type_combo.currentText() == "Data": 
			if str(self.data_path.text()) == "": self._on_browse()
			insertion_node =  EMDataItem3D(str(self.data_path.text()), transform=Transform())
			if self.node_name_data.text() != "": node_name = self.node_name_data.text()
		# Isosurface
		if self.node_type_combo.currentText() == "Isosurface": 
			insertion_node =  EMIsosurface(self.item.item3d(), transform=Transform())
			if self.node_name_data.text() != "": node_name = self.node_name_data.text()
			itemparentnode=self.item
			
		thisnode_index = itemparentnode.indexOfChild(self.item)
		self.inspector().addTreeNode(node_name, insertion_node, parentnode=itemparentnode, insertionindex=(thisnode_index+1))
		itemparentnode.item3d().addChild(insertion_node)
		self.inspector().updateSceneGraph()
		self.done(0)
	
	def _on_cancel(self):
		self.done(1)
	
	def _node_combobox_changed(self, nodetype):
		self.node_stacked_widget.setCurrentIndex(nodetype)
	
class EMInspectorControlShape(EMItem3DInspector):
	"""
	Class to make EMItem GUI SHAPE Inspector
	"""
	def __init__(self, name, item3d):
		EMItem3DInspector.__init__(self, name, item3d)
		
	def addControls(self, igvbox):
		pass
		
	def addColorControls(self, box):
		colorframe = QtGui.QFrame()
		colorframe.setFrameShape(QtGui.QFrame.StyledPanel)
		colorvbox = QtGui.QVBoxLayout()
		lfont = QtGui.QFont()
		lfont.setBold(True)
		colorlabel = QtGui.QLabel("Color",colorframe)
		colorlabel.setFont(lfont)
		colorlabel.setAlignment(QtCore.Qt.AlignCenter)

		# These boxes are a pain maybe I should use a Grid?
		cdialoghbox = QtGui.QHBoxLayout()
		cabox = QtGui.QHBoxLayout()
		self.ambcolorbox = EMQTColorWidget(parent=colorframe)
		self.ambcolorbox.setColor(QtGui.QColor(255*self.item3d().ambient[0],255*self.item3d().ambient[1],255*self.item3d().ambient[2]))
		cabox.addWidget(self.ambcolorbox)
		cabox.setAlignment(QtCore.Qt.AlignCenter)
		cdbox = QtGui.QHBoxLayout()
		self.diffusecolorbox = EMQTColorWidget(parent=colorframe)
		self.diffusecolorbox.setColor(QtGui.QColor(255*self.item3d().diffuse[0],255*self.item3d().diffuse[1],255*self.item3d().diffuse[2]))
		cdbox.addWidget(self.diffusecolorbox)
		cdbox.setAlignment(QtCore.Qt.AlignCenter)
		csbox = QtGui.QHBoxLayout()
		self.specularcolorbox = EMQTColorWidget(parent=colorframe)
		self.specularcolorbox.setColor(QtGui.QColor(255*self.item3d().specular[0],255*self.item3d().specular[1],255*self.item3d().specular[2]))
		csbox.addWidget(self.specularcolorbox)
		csbox.setAlignment(QtCore.Qt.AlignCenter)
		cdialoghbox.addLayout(cabox)
		cdialoghbox.addLayout(cdbox)
		cdialoghbox.addLayout(csbox)
		
		colorhbox = QtGui.QHBoxLayout()
		self.ambient = QtGui.QLabel("Ambient", colorframe)
		self.ambient.setAlignment(QtCore.Qt.AlignCenter)
		self.diffuse = QtGui.QLabel("Diffuse", colorframe)
		self.diffuse.setAlignment(QtCore.Qt.AlignCenter)
		self.specular = QtGui.QLabel("Specular", colorframe)
		self.specular.setAlignment(QtCore.Qt.AlignCenter)
		colorhbox.addWidget(self.ambient)
		colorhbox.addWidget(self.diffuse)
		colorhbox.addWidget(self.specular)
		
		self.shininess = ValSlider(colorframe, (0.0, 50.0), "Shininess")
		self.shininess.setValue(self.item3d().shininess)
		
		colorvbox.addWidget(colorlabel)
		colorvbox.addLayout(cdialoghbox)
		colorvbox.addLayout(colorhbox)
		colorvbox.addWidget(self.shininess)
		colorframe.setLayout(colorvbox)
		box.addWidget(colorframe)
		
		QtCore.QObject.connect(self.ambcolorbox,QtCore.SIGNAL("newcolor(QColor)"),self._on_ambient_color)
		QtCore.QObject.connect(self.diffusecolorbox,QtCore.SIGNAL("newcolor(QColor)"),self._on_diffuse_color)
		QtCore.QObject.connect(self.specularcolorbox,QtCore.SIGNAL("newcolor(QColor)"),self._on_specular_color)
		QtCore.QObject.connect(self.shininess,QtCore.SIGNAL("valueChanged"),self._on_shininess)
		
	def _on_ambient_color(self, color):
		rgb = color.getRgb()
		self.item3d().setAmbientColor((float(rgb[0])/255.0),(float(rgb[1])/255.0),(float(rgb[2])/255.0))
		self.inspector.updateSceneGraph()
		
	def _on_diffuse_color(self, color):
		rgb = color.getRgb()
		self.item3d().setDiffuseColor((float(rgb[0])/255.0),(float(rgb[1])/255.0),(float(rgb[2])/255.0))
		self.inspector.updateSceneGraph()
		
	def _on_specular_color(self, color):
		rgb = color.getRgb()
		self.item3d().setSpecularColor((float(rgb[0])/255.0),(float(rgb[1])/255.0),(float(rgb[2])/255.0))
		self.inspector.updateSceneGraph()
		
	def _on_shininess(self, shininess):
		self.item3d().setShininess(shininess)
		self.inspector.updateSceneGraph()
		
###################################### TEST CODE, THIS WILL NOT APPEAR IN THE WIDGET3D MODULE ##################################################
		
# All object that are rendered inherit from abstractSGnode and implement the render method
from emshapeitem3d import *
from emdataitem3d import EMDataItem3D, EMIsosurface

class GLdemo(QtGui.QWidget):
	def __init__(self):
		QtGui.QWidget.__init__(self)
		self.widget = EMScene3D()
		#self.widget.camera.usePrespective(50, 0.5)
		self.cube = EMCube(50.0)
		self.widget.addChild(self.cube)    # Something to Render something..... (this could just as well be one of Ross's SGnodes)
		self.sphere = EMSphere(50.0)
		self.widget.addChild(self.sphere)
		self.cylider = EMCylinder(50.0, 50.0)
		self.widget.addChild(self.cylider)
		#data = EMData("/home/john/Bo_data/simulated_data/3DmapIP3R1_small.mrc")
		#self.emdata = EMDataItem3D("/home/john/Bo_data/simulated_data/3DmapIP3R1_small.mrc", transform=Transform())
		#self.widget.addChild(self.emdata)
		#self.isosurface = EMIsosurface(self.emdata, transform=Transform())
		#self.emdata.addChild(self.isosurface)
		
		self.inspector = EMInspector3D(self.widget)
		self.widget.setInspector(self.inspector)
		
		rootnode = self.inspector.addTreeNode("root node", self.widget)
		self.inspector.addTreeNode("cube", self.cube, rootnode)
		self.inspector.addTreeNode("sphere", self.sphere, rootnode)
		self.inspector.addTreeNode("cylider", self.cylider, rootnode)
		#datanode = self.inspector.addTreeNode("data", self.emdata, rootnode)
		#self.inspector.addTreeNode("isosurface", self.isosurface, datanode)
		
		# QT stuff to display the widget
		vbox = QtGui.QVBoxLayout()
		vbox.addWidget(self.widget)
		self.setLayout(vbox)
		self.setGeometry(300, 300, 600, 600)
		self.setWindowTitle('BCM EM Viewer')
	
	def show_inspector(self):
		self.inspector.show()
		
if __name__ == "__main__":
	import sys
	app = QtGui.QApplication(sys.argv)
	window = GLdemo()
	window.show()
	window.show_inspector()
	app.exec_()
