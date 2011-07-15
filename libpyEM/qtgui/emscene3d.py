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
from valslider import ValSlider, EMSpinWidget, EMQTColorWidget, EMLightControls, CameraControls
import math
import weakref

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
		self.camera = EMCamera(1.0, 10000.0)	# Default near,far, and zclip values
		self.main_3d_inspector = None
		self.widget = None				# Get the inspector GUI
		self.z_clipped_changed = False			# Toogle flag to deterine if the clipping plane has changed and needs redrawing
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

	def initializeGL(self):
		glClearColor(0.0, 0.0, 0.0, 0.0)		# Default clear color is black
		glShadeModel(GL_SMOOTH)
		glEnable(GL_DEPTH_TEST)
		self.firstlight = EMLight(GL_LIGHT0)
		self.firstlight.enableLighting()
        
	def paintGL(self):
		if self.z_clipped_changed: self.camera.update()
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)		
		glColor3f(1.0, 1.0, 1.0)	# Default color is white
		#Call rendering
		self.renderSelectedArea() 	# Draw the selection box if needed
		self.render()			# SG nodes must have a render method
		glFlush()			# Finish rendering
		self.z_clipped_changed = False
		
	def resizeGL(self, width, height):
		self.camera.update(width, height)
	
	def getSceneGui(self):
		"""
		Return a Qt widget that controls the scene item
		"""	
		if not self.widget: self.widget = EMItem3DInspector("SG", self)
		return self.widget
		
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
		
		# Find the selection box. Go from Volume view coords to viewport coords
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
		self.camera.setCameraPosition(sfactor=2) # Factor of two to compensate for the samera already being set
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
			z = -self.camera.getZclip() - 1
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
				self.main_3d_inspector.stacked_widget.setCurrentWidget(selecteditem.getSceneGui())
				self.main_3d_inspector.tree_widget.setCurrentItem(selecteditem.getSceneGui().treeitem)	# Hmmm... tak about tight coupling! It would be better to use getters
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
		Qt event handler. Returns the cursor to arrow unpn mouse button release
		"""
		self.setCursor(Qt.ArrowCursor)
		if self.toggle_render_selectedarea:
			self.pickItem()
			self.deselectArea()
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
	
	def getFuzzyFactor(self):
		"""Get the fuzzy selection factor """
		return self.fuzzyselectionfactor
		
	def setFuzzyFactor(self, factor):
		"""Set the fuzzy selection factor"""
		self.fuzzyselectionfactor = factor
		
	def setZclip(self, zclip):
		""" Set the Z clipping plane """
		self.camera.setZclip(zclip)
		self.z_clipped_changed = True
		
	def setZslice(self):
		"""
		Get a Z slice to display in the camera widget, only works for orthoganal mode
		"""
		# Getting the Z slice will have problems when using perspective viewing
		self.setAutoBufferSwap(False)
		self.getTransform().rotate_origin(Transform({"type":"spin","Omega":90,"n1":0,"n2":1,"n3":0}))
		oldzclip = self.camera.getZclip()
		self.camera.setZclip(-1000) 		# We must move the camera back when we take a side slice
		self.z_clipped_changed = True
		QtOpenGL.QGLWidget.updateGL(self)
		pixeldata = glReadPixels(10,10,self.camera.width,self.camera.height,GL_RGB,GL_UNSIGNED_BYTE)
		self.getTransform().rotate_origin(Transform({"type":"spin","Omega":-90,"n1":0,"n2":1,"n3":0}))
		self.camera.setZclip(oldzclip)
		self.z_clipped_changed = True
		self.setAutoBufferSwap(True)
		self.pixels = []
		self.pixels.append(0)
		self.pixels.append(0)
		self.pixels.append(self.camera.width)
		self.pixels.append(self.camera.height)
		self.pixels.append(pixeldata)
		
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
		self.setAmbient(0.1, 0.1, 0.1, 1.0)		# Default ambient color is light grey
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
	def __init__(self, near, far, zclip=-500.0, usingortho=True, fovy=60.0, boundingbox=50.0, screenfraction=0.5):
		"""
		@param fovy: The field of view angle
		@param near: The volume view near position
		@param far: The volume view far position
		@param zclip: The zclipping plane (basicaly how far back the camera is)
		@param usingortho: Use orthographic projection
		@param boundingbox: The dimension of the bounding for the object to be rendered
		@param screenfraction: The fraction of the screen height to occupy
		"""
		self.far = far
		self.near = near
		self.fovy = fovy
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
		self.perspective_z = -(boundingbox/screenfraction)/(2*math.tan(math.radians(self.fovy/2)))  + boundingbox/2
		self.usingortho = False
		

	def useOrtho(self, zclip):
		"""
		Changes projection matrix to orthographic
		"""
		self.usingortho = True
		self.zclip = zclip

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
		
	def getZclip(self):
		""" Get the zclip """
		if self.usingortho:
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
		treesplitter = QtGui.QSplitter(QtCore.Qt.Horizontal)
		treesplitter.setFrameShape(QtGui.QFrame.StyledPanel)
		treesplitter.setLayout(self._get_tree_layout(widget))
		treesplitter.setMinimumWidth(self.mintreewidth)
		hbox.addWidget(treesplitter)
		controlsplitter = QtGui.QSplitter(QtCore.Qt.Horizontal)
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
		
		QtCore.QObject.connect(self.tree_widget, QtCore.SIGNAL("itemClicked(QTreeWidgetItem*,int)"), self._tree_widget_click)
		QtCore.QObject.connect(self.tree_widget, QtCore.SIGNAL("visibleItem(QTreeWidgetItem*)"), self._tree_widget_visible)
		
		return tvbox
	
	def addTreeNode(self, name, item3d, parentnode=None):
		"""
		Add a node (item3d) to the TreeWidget if not parent node, otherwise add a child to parent node
		We need to get a GUI for the treeitem. The treeitem and the GUI need know each other so they can talk
		The Treeitem also needs to know the node, so it can talk to the node.
		You can think of this as a three way conversation (the alterative it to use a mediator, but that is not worth it w/ only three players)
		"""
		tree_item = EMQTreeWidgetItem(QtCore.QStringList(name), item3d)	# Make a QTreeItem widget, and let the TreeItem talk to the scenegraph node and its GUI
		item3d.setEMQTreeWidgetItem(tree_item)				# Reference to the EMQTreeWidgetItem
		item3dgui = item3d.getSceneGui()				# Get the node GUI controls 
		item3dgui.setInspector(self)					# Associate the item GUI with the inspector
		self.stacked_widget.addWidget(item3dgui)			# Add a widget to the stack
		# Set icon status
		tree_item.setSelectionStateBox()
		# Set parent if one exists	
		if not parentnode:
			self.tree_widget.insertTopLevelItem(0, tree_item)
		else:
			parentnode.addChild(tree_item)
		return tree_item
			
	def _tree_widget_click(self, item, col):
		"""
		When a user clicks on the selection tree check box
		"""
		self.stacked_widget.setCurrentWidget(item.item3d.getSceneGui())
		item.setSelectionState(item.checkState(0))
		for ancestor in item.item3d.getSelectedAncestorNodes():
			ancestor.EMQTreeWidgetItem.setSelectionState(False)
		for child in item.item3d.getAllSelectedNodes()[1:]: # Lop the node itself off
			child.EMQTreeWidgetItem.setSelectionState(False)
		
	def _tree_widget_visible(self, item):
		"""
		When a user clicks on the visible icon
		"""
		item.toggleVisibleState()
		self.updateSceneGraph()
		
	def _get_controler_layout(self, parent):
		"""
		Returns the control layout
		"""
		cvbox = QtGui.QVBoxLayout()
		self.stacked_widget = QtGui.QStackedWidget()
		cvbox.addWidget(self.stacked_widget)
		
		return cvbox
		
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
		valslidersplitter = QtGui.QSplitter(QtCore.Qt.Vertical)
		valslidersplitter.setFrameShape(QtGui.QFrame.StyledPanel)
		valvbox = QtGui.QVBoxLayout()
		self.hvalslider = ValSlider(lwidget,(0.0,360.0),"Horizontal")
		self.hvalslider.setMaximumHeight(40.0)
		self.vvalslider = ValSlider(lwidget,(0.0,360.0),"Vertical")
		self.vvalslider.setMaximumHeight(40.0)
		valvbox.addWidget(self.hvalslider)
		valvbox.addWidget(self.vvalslider)
		valslidersplitter.setLayout(valvbox)
		lvbox.addWidget(self.lightwidget)
		lvbox.addWidget(positionlabel)
		lvbox.addWidget(valslidersplitter)
		lwidget.setLayout(lvbox)
		
		QtCore.QObject.connect(self.lightwidget, QtCore.SIGNAL("lightPositionMoved"), self._light_position_moved)
		QtCore.QObject.connect(self.hvalslider,QtCore.SIGNAL("valueChanged"),self._on_light_slider)
		QtCore.QObject.connect(self.vvalslider,QtCore.SIGNAL("valueChanged"),self._on_light_slider)

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
		
	def getCameraWidget(self):
		"""
		Returns the camera control widget
		"""
		self.cameratab_open = False
		cwidget = QtGui.QWidget()
		cvbox = QtGui.QVBoxLayout()
		self.camerawidget = CameraControls(scenegraph=self.scenegraph)
		
		cvbox.addWidget(self.camerawidget)
		clabel = QtGui.QLabel("Camera position", cwidget)
		clabel.setMaximumHeight(30.0)
		clabel.setAlignment(QtCore.Qt.AlignCenter)
		self.camera_z = EMSpinWidget(self.scenegraph.camera.getZclip(), 1.0)
		self.camera_z.setMaximumHeight(40.0)
		cvbox.addWidget(clabel)
		cvbox.addWidget(self.camera_z)
		cwidget.setLayout(cvbox)

		QtCore.QObject.connect(self.camera_z,QtCore.SIGNAL("valueChanged(int)"),self._on_camera_z)
		QtCore.QObject.connect(self.camerawidget,QtCore.SIGNAL("cameraMoved(int)"),self._on_camera_move)
		
		return cwidget
	
	def _on_camera_z(self, value):
		self.camerawidget.updateWidget()
		self.scenegraph.setZclip(value)
		self.scenegraph.updateSG()
		
	def _on_camera_move(self, movement):
		value = self.scenegraph.camera.getZclip() + movement
		self.camera_z.setValue(value, quiet=1)
		self.scenegraph.setZclip(value)
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
			
	def getUtilsWidget(self):
		"""
		Return the utilites widget
		"""
		uwidget = QtGui.QWidget()
		uvbox = QtGui.QVBoxLayout()
		self.fuzzy_slider = ValSlider(uwidget, (0.2, 5.0), "Fuzzy Selection factor")
		self.fuzzy_slider.setValue(self.scenegraph.getFuzzyFactor(), quiet=1)
		uvbox.addWidget(self.fuzzy_slider)
		uwidget.setLayout(uvbox)
		
		QtCore.QObject.connect(self.fuzzy_slider, QtCore.SIGNAL("valueChanged"),self._on_fuzzy_slider)
		
		return uwidget
	
	def _on_fuzzy_slider(self, value):
		"""
		Set the fuzzy selection coeff
		"""
		self.scenegraph.setFuzzyFactor(value)
		
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
	def __init__(self, qstring, item3d):
		QtGui.QTreeWidgetItem.__init__(self, qstring)
		self.item3d = item3d
		self.setCheckState(0, QtCore.Qt.Unchecked)
		self.visible = QtGui.QIcon(QtGui.QPixmap(visibleicon))
		self.invisible = QtGui.QIcon(QtGui.QPixmap(invisibleicon))
		self.getVisibleState()
	
	def setSelectionState(self, state):
		""" 
		Toogle selection state on and off
		"""
		if state == QtCore.Qt.Checked:
			self.item3d.setSelectedItem(True)
		else:
			self.item3d.setSelectedItem(False)
		self.setSelectionStateBox() # set state of TreeItemwidget
		
	def toggleVisibleState(self):
		self.item3d.setVisibleItem(not self.item3d.isVisibleItem())
		self.getVisibleState()
		
	def getVisibleState(self):
		"""
		Toogle the visble state
		"""
		if self.item3d.isVisibleItem():
			self.setIcon(0, self.visible)
		else:
			self.setIcon(0, self.invisible)
	
	def setSelectionStateBox(self):
		"""
		Set the selection state icon
		"""
		if self.item3d.isSelectedItem():
			self.setCheckState(0, QtCore.Qt.Checked)
		else:
			self.setCheckState(0, QtCore.Qt.Unchecked)

class EMInspectorControlShape(EMItem3DInspector):
	"""
	Class to make EMItem GUI SHAPE
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
		self.ambcolorbox.setColor(QtGui.QColor(255*self.item3d.ambient[0],255*self.item3d.ambient[1],255*self.item3d.ambient[2]))
		cabox.addWidget(self.ambcolorbox)
		cabox.setAlignment(QtCore.Qt.AlignCenter)
		cdbox = QtGui.QHBoxLayout()
		self.diffusecolorbox = EMQTColorWidget(parent=colorframe)
		self.diffusecolorbox.setColor(QtGui.QColor(255*self.item3d.diffuse[0],255*self.item3d.diffuse[1],255*self.item3d.diffuse[2]))
		cdbox.addWidget(self.diffusecolorbox)
		cdbox.setAlignment(QtCore.Qt.AlignCenter)
		csbox = QtGui.QHBoxLayout()
		self.specularcolorbox = EMQTColorWidget(parent=colorframe)
		self.specularcolorbox.setColor(QtGui.QColor(255*self.item3d.specular[0],255*self.item3d.specular[1],255*self.item3d.specular[2]))
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
		self.shininess.setValue(self.item3d.shininess)
		
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
		self.item3d.setAmbientColor((float(rgb[0])/255.0),(float(rgb[1])/255.0),(float(rgb[2])/255.0))
		self.inspector.updateSceneGraph()
		
	def _on_diffuse_color(self, color):
		rgb = color.getRgb()
		self.item3d.setDiffuseColor((float(rgb[0])/255.0),(float(rgb[1])/255.0),(float(rgb[2])/255.0))
		self.inspector.updateSceneGraph()
		
	def _on_specular_color(self, color):
		rgb = color.getRgb()
		self.item3d.setSpecularColor((float(rgb[0])/255.0),(float(rgb[1])/255.0),(float(rgb[2])/255.0))
		self.inspector.updateSceneGraph()
		
	def _on_shininess(self, shininess):
		self.item3d.setShininess(shininess)
		self.inspector.updateSceneGraph()
		
###################################### TEST CODE, THIS WILL NOT APPEAR IN THE WIDGET3D MODULE ##################################################
		
# All object that are rendered inherit from abstractSGnode and implement the render method
from emshapeitem3d import *

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

		self.inspector = EMInspector3D(self.widget)
		self.widget.setInspector(self.inspector)
		
		rootnode = self.inspector.addTreeNode("root node", self.widget)
		self.inspector.addTreeNode("cube", self.cube, rootnode)
		self.inspector.addTreeNode("sphere", self.sphere, rootnode)
		self.inspector.addTreeNode("cylider", self.cylider, rootnode)
		
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
