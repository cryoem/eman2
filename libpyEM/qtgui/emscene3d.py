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
import weakref

# XPM format Cursors

selectedicon = [
    '13 12 2 1',
    'a c #00ff00',
    'c c None',
    'cccccccccccaa',
    'cccccccccccaa',
    'ccccccccccaac',
    'ccccccccccaac',
    'cccccccccaacc',
    'cccccccccaacc',
    'ccccccccaaccc',
    'ccccccccaaccc',
    'aacccccaacccc',
    'caaccccaacccc',
    'ccaaaaaaccccc',
    'ccccaaacccccc'
]
    
deselectedicon = [
    '12 11 2 1',
    'a c #ff0000',
    'c c None',
    'aaccccccccaa',
    'caaccccccaac',
    'ccaaccccaacc',
    'cccaaccaaccc',
    'ccccaaaacccc',
    'cccccaaccccc',
    'ccccaaaacccc',
    'cccaaccaaccc',
    'ccaaccccaacc',
    'caaccccccaac',
    'aaccccccccaa'
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
		self.widget = EMInspectorControl("SG", self)	# Get the inspector GUI
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
		self.firstlight.enablelighting()
        
	def paintGL(self):
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)		
		glColor3f(1.0, 1.0, 1.0)	# Default color is white
		#Call rendering
		self.render_selectedarea() 	# Draw the selection box if needed
		self.render()			# SG nodes must have a render method
		glFlush()			# Finish rendering

	def resizeGL(self, width, height):
		self.camera.update(width, height)
	
	def get_scene_gui(self):
		"""
		Return a Qt widget that controls the scene item
		"""	
		return self.widget
		
	def render_node(self):
		pass
	
	def set_inspector(self, inspector):
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
		x = self.sa_xi + self.camera.get_width()/2
		y = self.camera.get_height()/2 - self.sa_yi
		dx = 2*math.fabs(self.sa_xi - self.sa_xf) # The 2x is a hack.....
		dy = 2*math.fabs(self.sa_yi - self.sa_yf) # The 2x is a hack.....
		
		# Apply selection box
		GLU.gluPickMatrix(x, viewport[3] - y, dx, dy, viewport)
		self.camera.setprojectionmatrix()
		
		#drawstuff, but first we need to remove the influence of any previous xforms which ^$#*$ the selection
		glMatrixMode(GL_MODELVIEW)
		glPushMatrix()
		glLoadIdentity()
		self.camera.set_camera_position(sfactor=2) # Factor of two to compensate for the samera already being set
		self.render()
		glPopMatrix()
		
		# Return to default state
		glMatrixMode(GL_PROJECTION)
		glPopMatrix()
		glMatrixMode(GL_MODELVIEW)
		records = glRenderMode(GL_RENDER)
		
		# process records
		self.processselection(records)
	
	def selectarea(self, xi, xf, yi, yf):
		"""
		Set an area for selection. Need to switch bewteen viewport coords, where (0,0 is bottom left) to
		volume view coords where 0,0) is center of the screen.
		"""
		self.sa_xi = xi - self.camera.get_width()/2
		self.sa_xf = xf - self.camera.get_width()/2
		self.sa_yi = -yi + self.camera.get_height()/2
		self.sa_yf = -yf + self.camera.get_height()/2
		self.toggle_render_selectedarea = True
		
	def deselectarea(self):
		"""
		Turn off selectin box
		"""
		self.sa_xi = 0.0
		self.sa_xf = 0.0
		self.sa_yi = 0.0
		self.sa_yf = 0.0
		self.toggle_render_selectedarea = False
		
	def render_selectedarea(self):
		"""
		Draw the selection box, box is always drawn orthographically
		"""
		if self.toggle_render_selectedarea: 
			glMatrixMode(GL_PROJECTION)
			glPushMatrix()
			glLoadIdentity()
			self.camera.set_ortho_projectionmatrix()
			glColor3f(0.0,1.0,0.0)
			glMaterialfv(GL_FRONT, GL_EMISSION, [0.0,1.0,0.0,1.0])
			glBegin(GL_LINE_LOOP)
			z = -self.camera.get_zclip() - 1
			glVertex3f(self.sa_xi, self.sa_yi, z)
			glVertex3f(self.sa_xi, self.sa_yf, z)
			glVertex3f(self.sa_xf, self.sa_yf, z)
			glVertex3f(self.sa_xf, self.sa_yi, z)
			glEnd()
			glPopMatrix()
			glMatrixMode(GL_MODELVIEW)
			glMaterialfv(GL_FRONT, GL_EMISSION, [0.0,0.0,0.0,1.0])
		
	def processselection(self, records):
		"""
		Process the selection records
		"""
		# Remove old selection if not in append mode
		if not self.appendselection:
			for selected in self.get_all_selected_nodes():
				selected.is_selected = False
				selected.widget.update_gui()
		# Select the desired items	
		closestitem = None
		bestdistance = 1.0
		for record in records:
			selecteditem = EMItem3D.selection_idx_dict[record.names[len(record.names)-1]]()
			selecteditem.is_selected = True
			selecteditem.widget.update_gui()
			try:
				self.main_3d_inspector.stacked_widget.setCurrentWidget(selecteditem.widget)
				self.main_3d_inspector.tree_widget.setCurrentItem(selecteditem.widget.treeitem)	# Hmmm... tak about tight coupling! It would be better to use getters
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
				self.selectarea(self.first_x, event.x(), self.first_y, event.y())
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
			self.deselectarea()
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
	def __init__(self, near, far, zclip=-1000.0, usingortho=True, fovy=60.0, boundingbox=50.0, screenfraction=0.5):
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
		if usingortho:
			self.useortho(zclip)
		else:
			self.useprespective(boundingbox, screenfraction, fovy)

	def update(self, width, height):
		"""
		@param width: The width of the window in pixels
		@param height: The height of the window in pixels
		updates the camera and viewport after windowresize
		"""
		self.width = width
		self.height = height
		if self.usingortho:
			glViewport(0,0,width,height)
			glMatrixMode(GL_PROJECTION)
			glLoadIdentity()
			self.set_ortho_projectionmatrix()
			glMatrixMode(GL_MODELVIEW)
			glLoadIdentity()
			glTranslate(0,0,self.zclip)
			self.set_camera_position()
		else:
			# This may need some work to get it to behave
			glViewport(0,0,width,height)
			glMatrixMode(GL_PROJECTION)
			glLoadIdentity()
			self.set_perspective_projectionmatrix()
			glMatrixMode(GL_MODELVIEW)
			glLoadIdentity()
			glTranslate(0,0,self.perspective_z) #How much to set the camera back depends on how big the object is
			self.set_camera_position()
	
	def set_camera_position(self, sfactor=1):
		"""
		Set the default camera position
		"""
		glTranslate(0,0,sfactor*self.get_zclip())
		
	def setprojectionmatrix(self):
		"""
		Set the projection matrix
		"""
		if self.usingortho:
			self.set_ortho_projectionmatrix()
		else:
			self.set_perspective_projectionmatrix()
			
	def set_ortho_projectionmatrix(self):
		"""
		Set the orthographic projection matrix. Volume view origin (0,0) is center of screen
		"""
		glOrtho(-self.width/2, self.width/2, -self.height/2, self.height/2, self.near, self.far)
		
	def set_perspective_projectionmatrix(self):
		"""
		Set the perspective projection matrix. Volume view origin (0,0) is center of screen
		"""
		GLU.gluPerspective(self.fovy, (float(self.width)/float(self.height)), self.near, self.far)
			
	def useprespective(self, boundingbox, screenfraction, fovy=60.0):
		""" 
		@param boundingbox: The dimension of the bounding for the object to be rendered
		@param screenfraction: The fraction of the screen height to occupy
		Changes projection matrix to perspective
		"""
		self.fovy = fovy
		self.perspective_z = -(boundingbox/screenfraction)/(2*math.tan(math.radians(self.fovy/2)))  + boundingbox/2
		self.usingortho = False
		

	def useortho(self, zclip):
		"""
		Changes projection matrix to orthographic
		"""
		self.usingortho = True
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
		
	def get_height(self):
		""" Get the viewport height """
		return self.height
	
	def get_width(self):
		""" Get the viewport width"""
		return self.width
		
	def get_zclip(self):
		""" Get the zclip """
		if self.usingortho:
			return self.zclip
		else:
			return self.perspective_z
	# Maybe other methods to control the camera

###################################### Inspector Code #########################################################################################

class EMInspector3D(QtGui.QWidget):
	def __init__(self):
		"""
		"""
		QtGui.QWidget.__init__(self)
		
		vbox = QtGui.QVBoxLayout(self)
		
		self.inspectortab = QtGui.QTabWidget()
		self.inspectortab.addTab(self.get_tree_widget(), "Tree View")
		self.inspectortab.addTab(self.get_lights_widget(), "Lights")
		self.inspectortab.addTab(self.get_camera_widget(), "Camera")
		self.inspectortab.addTab(self.get_utils_widget(), "Utils")

		vbox.addWidget(self.inspectortab)
		
		self.setLayout(vbox)
		self.updateGeometry()

	def get_tree_widget(self):
		"""
		This returns the treeview-control panel widget
		"""
		widget = QtGui.QWidget()
		hbox = QtGui.QHBoxLayout(widget)
		treesplitter = QtGui.QSplitter(QtCore.Qt.Vertical)
		treesplitter.setFrameShape(QtGui.QFrame.StyledPanel)
		treesplitter.setLayout(self.get_tree_layout(widget))
		hbox.addWidget(treesplitter)
		controlsplitter = QtGui.QSplitter(QtCore.Qt.Vertical)
		controlsplitter.setFrameShape(QtGui.QFrame.StyledPanel)
		controlsplitter.setLayout(self.get_controler_layout(widget))
		hbox.addWidget(controlsplitter)
		widget.setLayout(hbox)
		
		return widget
		
	def get_tree_layout(self, parent):
		"""
		Returns the tree layout
		"""
		tvbox = QtGui.QVBoxLayout()
		self.tree_widget = QtGui.QTreeWidget(parent)
		self.tree_widget.setHeaderLabel("Choose a item")
		tvbox.addWidget(self.tree_widget)
		
		QtCore.QObject.connect(self.tree_widget, QtCore.SIGNAL("itemClicked(QTreeWidgetItem*,int)"), self.tree_widget_click)
		
		return tvbox
	
	def add_tree_node(self, name, sgnode, parentnode=None):
		"""
		Add a node to the TreeWidget if not parent node, otherwise add a child to parent node
		We need to get a GUI for the treeitem. The treeitem and the GUI need know each other so they can talk
		The Treeitem also needs to know the SGnode, so it can talk to the SGnode.
		You can think of this as a three way conversation(the alterative it to use a meniator, but that is not worth it w/ only three players
		"""
		sgnodegui = sgnode.get_scene_gui()				# Get the SG node GUI controls 
		node = EMQTreeWidgetItem(QtCore.QStringList(name), sgnode, sgnodegui)	# Make a QTreeItem widget, and let the TreeItem talk to the scenegraph node and its GUI 
		sgnodegui.set_treeitem(node)					# Relate the SG node GUI controls to the QTreeItem
		self.stacked_widget.addWidget(sgnodegui)			# Add a widget to the stack
		# Set icon status
		node.set_selection_state_box()
		# Set parent if one exists	
		if not parentnode:
			self.tree_widget.insertTopLevelItem(0, node)
		else:
			parentnode.addChild(node)
		return node
			
	def tree_widget_click(self, item, col):
		self.stacked_widget.setCurrentWidget(item.sgnode.get_scene_gui())
		item.set_selection_state(item.checkState(0))
		
	def get_controler_layout(self, parent):
		"""
		Returns the control layout
		"""
		cvbox = QtGui.QVBoxLayout()
		self.stacked_widget = QtGui.QStackedWidget()
		cvbox.addWidget(self.stacked_widget)
		return cvbox
		
	def get_lights_widget(self):
		"""
		Returns the lights control widget
		"""
		lwidget = QtGui.QWidget()
		
		return lwidget
		
	def get_camera_widget(self):
		"""
		Returns the camera control widget
		"""
		cwidget = QtGui.QWidget()
		
		return cwidget
	
	def get_utils_widget(self):
		"""
		Retrusn the utilites widget
		"""
		uwidget = QtGui.QWidget()
		
		return uwidget

		
class EMQTreeWidgetItem(QtGui.QTreeWidgetItem):
	"""
	Subclass of QTreeWidgetItem
	adds functionality
	"""
	def __init__(self, qstring, sgnode, sgnodegui):
		QtGui.QTreeWidgetItem.__init__(self, qstring)
		self.sgnode = sgnode
		self.sgnodegui = sgnodegui
		self.setCheckState(0, QtCore.Qt.Unchecked)
	
	def set_selection_state(self, state):
		""" 
		Toogle selection state on and off
		"""
		if state == QtCore.Qt.Checked:
			self.sgnode.is_selected = True
		else:
			self.sgnode.is_selected = False
		self.set_selection_state_box() # set state of TreeItemwidget
		self.sgnodegui.set_selection_checkbox()	# Set the GUI checked box
		
	def set_selection_state_box(self):
		"""
		Set the selection state icon
		"""
		if self.sgnode.is_selected:
			self.setCheckState(0, QtCore.Qt.Checked)
		else:
			self.setCheckState(0, QtCore.Qt.Unchecked)
		
class EMInspectorControl(QtGui.QWidget):
	"""
	Class to make the EMItem GUI controls
	"""
	def __init__(self, name, sgnode):
		QtGui.QWidget.__init__(self)
		self.sgnode = sgnode
		
		igvbox = QtGui.QVBoxLayout()
		igvbox.addWidget(QtGui.QLabel(name,self))
		self.isselectedbox = QtGui.QCheckBox("Is Selected", self)
		self.set_selection_checkbox()
		igvbox.addWidget(self.isselectedbox)
		test =  EMColorWidget(150, 150, 0)
		igvbox.addWidget(test)
		self.setLayout(igvbox)
		
		self.connect(self.isselectedbox,QtCore.SIGNAL("stateChanged(int)"),self.on_isselected)
		
	def set_treeitem(self, treeitem):
		"""
		Relate the GUI to the treeitem so it can talk directy to the tree item
		"""
		self.treeitem = treeitem
		
	def set_selection_checkbox(self):
		"""
		Set the selection state
		"""
		if self.sgnode.is_selected:	
			self.isselectedbox.setCheckState(QtCore.Qt.Checked)
		else:
			self.isselectedbox.setCheckState(QtCore.Qt.Unchecked)
		
	def on_isselected(self):
		"""
		Set the selection state when the check boix is clicked
		"""
		if self.isselectedbox.isChecked():
			self.sgnode.is_selected = True
		else:
			self.sgnode.is_selected = False
		self.treeitem.set_selection_state_box()
		
	def update_gui(self):
		self.set_selection_checkbox()
		self.treeitem.set_selection_state_box()

class EMColorWidget(QtGui.QWidget):
	def __init__(self, redcolor, greencolor, bluecolor):
		QtGui.QWidget.__init__(self)
		# Make redbar
		cvbox = QtGui.QVBoxLayout()
		# Red
		rhbox = QtGui.QHBoxLayout()
		self.redbar = EMColorBar("red", redcolor)
		self.redspin = QtGui.QSpinBox(self)
		self.redspin.setRange(0,255)
		self.redspin.setMaximumWidth(60)
		self.redspin.setValue(redcolor)
		self.redbar.setHeight(self.redspin.height())
		rhbox.addWidget(self.redbar)
		rhbox.addWidget(self.redspin)
		cvbox.addLayout(rhbox)
		# Green
		ghbox = QtGui.QHBoxLayout()
		self.greenbar = EMColorBar("green", greencolor)
		self.greenspin = QtGui.QSpinBox(self)
		self.greenspin.setRange(0,255)
		self.greenspin.setMaximumWidth(60)
		self.greenspin.setValue(greencolor)
		self.greenbar.setHeight(self.greenspin.height())
		ghbox.addWidget(self.greenbar)
		ghbox.addWidget(self.greenspin)
		cvbox.addLayout(ghbox)
		# Blue
		bhbox = QtGui.QHBoxLayout()
		self.bluebar = EMColorBar("blue", bluecolor)
		self.bluespin = QtGui.QSpinBox(self)
		self.bluespin.setRange(0,255)
		self.bluespin.setMaximumWidth(60)
		self.bluespin.setValue(bluecolor)
		self.bluebar.setHeight(self.bluespin.height())
		bhbox.addWidget(self.bluebar)
		bhbox.addWidget(self.bluespin)
		cvbox.addLayout(bhbox)
		
		self.setLayout(cvbox)
		
		QtCore.QObject.connect(self.redbar,QtCore.SIGNAL("colorChanged(int)"),self.on_redbar)
		QtCore.QObject.connect(self.redspin,QtCore.SIGNAL("valueChanged(int)"),self.on_redspin)
		QtCore.QObject.connect(self.greenbar,QtCore.SIGNAL("colorChanged(int)"),self.on_greenbar)
		QtCore.QObject.connect(self.greenspin,QtCore.SIGNAL("valueChanged(int)"),self.on_greenspin)
		QtCore.QObject.connect(self.bluebar,QtCore.SIGNAL("colorChanged(int)"),self.on_bluebar)
		QtCore.QObject.connect(self.bluespin,QtCore.SIGNAL("valueChanged(int)"),self.on_bluespin)
		
	def on_redbar(self, color):
		self.redspin.setValue(color)
		self.redbar.setMarker(color)
		self.emit(QtCore.SIGNAL("redcolorChanged(int)"),color)
		self.redbar.update()
		
	def on_redspin(self, value):
		self.redbar.setMarker(value)
		self.emit(QtCore.SIGNAL("redcolorChanged(int)"),value)
		self.redbar.update()
		
	def on_greenbar(self, color):
		self.greenspin.setValue(color)
		self.greenbar.setMarker(color)
		self.emit(QtCore.SIGNAL("greencolorChanged(int)"),color)
		self.greenbar.update()
		
	def on_greenspin(self, value):
		self.greenbar.setMarker(value)
		self.emit(QtCore.SIGNAL("greencolorChanged(int)"),value)
		self.greenbar.update()
		
	def on_bluebar(self, color):
		self.bluespin.setValue(color)
		self.bluebar.setMarker(color)
		self.emit(QtCore.SIGNAL("greencolorChanged(int)"),color)
		self.bluebar.update()
		
	def on_bluespin(self, value):
		self.bluebar.setMarker(value)
		self.emit(QtCore.SIGNAL("greencolorChanged(int)"),value)
		self.bluebar.update()
	
class EMColorBar(QtGui.QWidget):
	def __init__(self, color, brightness):
		QtGui.QWidget.__init__(self)
		minwidth = 50
		self.setMinimumSize(minwidth, 2)
		self.color = color
		self.radius = 6
		self.w = minwidth
		self.setMarker(brightness)
		
	def setHeight(self, height):
		self.setMaximumHeight(height+self.radius/2)
		
	def setMarker(self, position):
		self.absposition = position
		self.position = int((float(self.w-self.radius)/255.0)*position)
		if self.position < 0: self.position = 0
		if self.position > (self.w-self.radius): self.position = (self.w-self.radius)
		
	def paintEvent(self, e):
		qp = QtGui.QPainter()
		qp.begin(self)
		self.drawWidget(qp)
		qp.end()
		
	def drawWidget(self, qp):
		size = self.size()
		self.w = size.width()
		self.h = size.height()
		self.setMarker(self.absposition)
		

		qp.setPen(QtGui.QColor(255, 255, 255))
		qp.setBrush(QtGui.QColor(0, 0, 0))
		qp.drawEllipse(self.position,0,self.radius,self.radius)
		
		if self.color == "red":
			for i in xrange(self.w-self.radius):
				color = int((255.0/float(self.w))*i)
				qcolor = QtGui.QColor(color, 0, 0)
				qp.setPen(qcolor)
				qp.setBrush(qcolor)
				qp.drawRect(i+self.radius/2, self.radius, 1, self.h-self.radius)
		if self.color == "green":
			for i in xrange(self.w-self.radius):
				color = int((255.0/float(self.w))*i)
				qcolor = QtGui.QColor(0, color, 0)
				qp.setPen(qcolor)
				qp.setBrush(qcolor)
				qp.drawRect(i+self.radius/2, self.radius, 1, self.h-self.radius)
		if self.color == "blue":
			for i in xrange(self.w-self.radius):
				color = int((255.0/float(self.w))*i)
				qcolor = QtGui.QColor(0, 0, color)
				qp.setPen(qcolor)
				qp.setBrush(qcolor)
				qp.drawRect(i+self.radius/2, self.radius, 1, self.h-self.radius)	
	
	def mousePressEvent(self, event):
		color = int((255.0/float(self.w-self.radius))*event.x())
		self.emit(QtCore.SIGNAL("colorChanged(int)"),color)
		
	def mouseMoveEvent(self, event):
		color = int((255.0/float(self.w-self.radius))*event.x())
		self.emit(QtCore.SIGNAL("colorChanged(int)"),color)

###################################### TEST CODE, THIS WILL NOT APPEAR IN THE WIDGET3D MODULE ##################################################
		
# All object that are rendered inherit from abstractSGnode and implement the render method
# In this example I use a cube, but any object can be drawn and so long as the object class inherits from abstractSGnode
class glCube(EMItem3D):
	def __init__(self, size):
		EMItem3D.__init__(self, parent=None, transform=Transform())
		# size
		self.xi = -size/2
		self.yi = -size/2
		self.zi = -size/2
		self.xf = size/2
		self.yf = size/2
		self.zf = size/2
		
		# color
		self.diffuse = [0.5,0.5,0.5,1.0]
		self.specular = [1.0,1.0,1.0,1.0]
		self.ambient = [1.0, 1.0, 1.0, 1.0]
		
		# GUI contols
		self.widget = EMInspectorControl("CUBE", self)
	
	#def on_selection(self):
		#self.diffuse = [0.0,0.5,0.0,1.0]
		#self.specular = [0.0,1.0,0.0,1.0]
		#self.ambient = [0.0, 1.0, 0.0, 1.0]
		
	def get_scene_gui(self):
		"""
		Return a Qt widget that controls the scene item
		"""
		return self.widget	
		
	def render_node(self):
			
		# Material properties of the box
		glMaterialfv(GL_FRONT, GL_DIFFUSE, self.diffuse)
		glMaterialfv(GL_FRONT, GL_SPECULAR, self.specular)
		glMaterialf(GL_FRONT, GL_SHININESS, 25.0)
		glMaterialfv(GL_FRONT, GL_AMBIENT, self.ambient)
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

		
class GLdemo(QtGui.QWidget):
	def __init__(self):
		QtGui.QWidget.__init__(self)
		self.widget = EMScene3D()
		#self.widget.camera.useprespective(50, 0.5)
		self.cube1 = glCube(50.0)
		self.widget.add_child(self.cube1)    # Something to Render something..... (this could just as well be one of Ross's SGnodes)
		#self.widget.activatenode(cube1)
		self.cube2 = glCube(50.0)
		self.widget.add_child(self.cube2)
		#self.widget.activatenode(cube2)

		self.inspector = EMInspector3D()
		self.inspector.setGeometry(300, 300, 400, 300)
		self.widget.set_inspector(self.inspector)
		
		rootnode = self.inspector.add_tree_node("root node", self.widget)
		self.inspector.add_tree_node("cube1", self.cube1, rootnode)
		self.inspector.add_tree_node("cube2", self.cube2, rootnode)
		
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
