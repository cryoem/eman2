#!/usr/bin/env python
#
# Author: Grant Tang (gtang@bcm.edu)
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

from PyQt4 import QtCore, QtGui, QtOpenGL 
from OpenGL.GLU import *
from OpenGL.GL import *
import math 
from emglobjects import init_glut
from emitem3d import EMItem3D, EMItem3DInspector
from EMAN2 import Transform
from valslider import EMQTColorWidget, ValSlider

class EMCube(EMItem3D):
	name = "Cube"
	nodetype = "ShapeNode"
	def __init__(self, size, transform=Transform()):
		EMItem3D.__init__(self, parent=None, children=set(), transform=transform)
		# size
		self.setSize(size)
		
		# color
		self.diffuse = [0.5,0.5,0.5,1.0]
		self.specular = [1.0,1.0,1.0,1.0]
		self.ambient = [1.0, 1.0, 1.0, 1.0]
		self.shininess = 25.0		
	
	def setSize(self, size):
		self.size = size
		self.boundingboxsize = size
		self.xi = -size/2
		self.yi = -size/2
		self.zi = -size/2
		self.xf = size/2
		self.yf = size/2
		self.zf = size/2
		if self.item_inspector: self.item_inspector.updateMetaData()
		
	def getEvalString(self):
		return "EMCube(%s)"%self.size
		
	def setAmbientColor(self, red, green, blue, alpha=1.0):
		self.ambient = [red, green, blue, alpha]

	def setDiffuseColor(self, red, green, blue, alpha=1.0):
		self.diffuse = [red, green, blue, alpha]
		
	def setSpecularColor(self, red, green, blue, alpha=1.0):
		self.specular = [red, green, blue, alpha]
	
	def setShininess(self, shininess):
		self.shininess = shininess
		
	def getItemInspector(self):
		"""
		Return a Qt widget that controls the scene item
		"""
		if not self.item_inspector: self.item_inspector = EMInspectorControlShape("CUBE", self)
		return self.item_inspector
	
	def renderNode(self):
		if self.is_selected and glGetIntegerv(GL_RENDER_MODE) == GL_RENDER: # No need to draw outline in selection mode
			#if glGetIntegerv(GL_RENDER_MODE) == GL_RENDER: print "X"
			glPushAttrib( GL_ALL_ATTRIB_BITS )
		
			# First render the cylinder, writing the outline to the stencil buffer
			glClearStencil(0)
			glClear( GL_STENCIL_BUFFER_BIT )
			glEnable( GL_STENCIL_TEST )
			glStencilFunc( GL_ALWAYS, 1, 0xFFFF )		# Write to stencil buffer
			glStencilOp( GL_KEEP, GL_KEEP, GL_REPLACE )	# Only pixels that pass the depth test are written to the stencil buffer
			glPolygonMode( GL_FRONT_AND_BACK, GL_FILL )	
			self.renderCube()
		
			# Then render the outline
			glStencilFunc( GL_NOTEQUAL, 1, 0xFFFF )		# The object itself is stenciled out
			glStencilOp( GL_KEEP, GL_KEEP, GL_REPLACE )
			glLineWidth( 4.0 )				# By increasing the line width only the outline is drawn
			glPolygonMode( GL_FRONT_AND_BACK, GL_LINE )
			glMaterialfv(GL_FRONT, GL_EMISSION, [0.0, 1.0, 0.0, 1.0])
			self.renderCube()
	
			glPopAttrib()
		else:
			self.renderCube()
			
	def renderCube(self):        
		# Material properties of the box
		glMaterialfv(GL_FRONT, GL_DIFFUSE, self.diffuse)
		glMaterialfv(GL_FRONT, GL_SPECULAR, self.specular)
		glMaterialf(GL_FRONT, GL_SHININESS, self.shininess)
		glMaterialfv(GL_FRONT, GL_AMBIENT, self.ambient)
		
		# The box itself along with normal vectors
		glBegin(GL_QUADS)	# Start Drawing The Cube

		# Front Face (note that the texture's corners have to match the quad's corners)
		glNormal3f(self.xi, self.yi, self.zi + 1); glVertex3f(self.xi, self.yi, self.zi)
		glNormal3f(self.xf, self.yi, self.zi + 1); glVertex3f(self.xf, self.yi, self.zi)
		glNormal3f(self.xf, self.yf, self.zi + 1); glVertex3f(self.xf, self.yf, self.zi)
		glNormal3f(self.xi, self.yf, self.zi + 1); glVertex3f(self.xi, self.yf, self.zi)

		# Back Face
		glNormal3f(self.xi - 1, self.yi, self.zi); glVertex3f(self.xi, self.yi, self.zi)
		glNormal3f(self.xi - 1, self.yi, self.zf); glVertex3f(self.xi, self.yi, self.zf)
		glNormal3f(self.xi - 1, self.yf, self.zf); glVertex3f(self.xi, self.yf, self.zf)
		glNormal3f(self.xi - 1, self.yf, self.zi); glVertex3f(self.xi, self.yf, self.zi)

		# Top Face
		glNormal3f(self.xi, self.yi, self.zf - 1); glVertex3f(self.xi, self.yi, self.zf)
		glNormal3f(self.xi, self.yf, self.zf - 1); glVertex3f(self.xi, self.yf, self.zf)
		glNormal3f(self.xf, self.yf, self.zf - 1); glVertex3f(self.xf, self.yf, self.zf)
		glNormal3f(self.xf, self.yi, self.zf - 1); glVertex3f(self.xf, self.yi, self.zf)

		# Bottom Face
		glNormal3f(self.xf + 1, self.yf, self.zf); glVertex3f(self.xf, self.yf, self.zf)
		glNormal3f(self.xf + 1, self.yf, self.zi); glVertex3f(self.xf, self.yf, self.zi)
		glNormal3f(self.xf + 1, self.yi, self.zi); glVertex3f(self.xf, self.yi, self.zi)
		glNormal3f(self.xf + 1, self.yi, self.zf); glVertex3f(self.xf, self.yi, self.zf)

		# Right face
		glNormal3f(self.xi, self.yf + 1, self.zi); glVertex3f(self.xi, self.yf, self.zi)
		glNormal3f(self.xi, self.yf + 1, self.zf); glVertex3f(self.xi, self.yf, self.zf)
		glNormal3f(self.xf, self.yf + 1, self.zf); glVertex3f(self.xf, self.yf, self.zf)
		glNormal3f(self.xf, self.yf + 1, self.zi); glVertex3f(self.xf, self.yf, self.zi)

		# Left Face
		glNormal3f(self.xi, self.yi - 1, self.zi); glVertex3f(self.xi, self.yi, self.zi)
		glNormal3f(self.xi, self.yi - 1, self.zf); glVertex3f(self.xi, self.yi, self.zf)
		glNormal3f(self.xf, self.yi - 1, self.zf); glVertex3f(self.xf, self.yi, self.zf)
		glNormal3f(self.xf, self.yi - 1, self.zi); glVertex3f(self.xf, self.yi, self.zi)

		glEnd()	# Done Drawing The Cube


class EMSphere(EMItem3D):
	name = "Sphere"
	nodetype = "ShapeNode"
	def __init__(self, radius, transform=Transform()):
		EMItem3D.__init__(self, parent=None, children=set(), transform=transform)
		# size
		self.setRadius(radius)
		
		# color
		self.diffuse = [0.5,0.5,0.5,1.0]
		self.specular = [1.0,1.0,1.0,1.0]
		self.ambient = [1.0, 1.0, 1.0, 1.0]
		self.shininess = 25.0
	
	def setRadius(self, radius):
		self.radius = radius
		self.boundingboxsize = radius
		self.slices = int(radius)
		self.stacks = int(radius)
		if self.item_inspector: self.item_inspector.updateMetaData()
		
	def getEvalString(self):
		return "EMSphere(%s)"%self.radius
		
	def setAmbientColor(self, red, green, blue, alpha=1.0):
		self.ambient = [red, green, blue, alpha]

	def setDiffuseColor(self, red, green, blue, alpha=1.0):
		self.diffuse = [red, green, blue, alpha]
		
	def setSpecularColor(self, red, green, blue, alpha=1.0):
		self.specular = [red, green, blue, alpha]
	
	def setShininess(self, shininess):
		self.shininess = shininess
		
	def getItemInspector(self):
		"""
		Return a Qt widget that controls the scene item
		"""
		if not self.item_inspector: self.item_inspector = EMInspectorControlShape("SPHERE", self)
		return self.item_inspector
	
	def renderNode(self):
		if self.is_selected and glGetIntegerv(GL_RENDER_MODE) == GL_RENDER:
			glPushAttrib( GL_ALL_ATTRIB_BITS )
		
			# First render the cylinder, writing the outline to the stencil buffer
			glClearStencil(0)
			glClear( GL_STENCIL_BUFFER_BIT )
			glEnable( GL_STENCIL_TEST )
			glStencilFunc( GL_ALWAYS, 1, 0xFFFF )		# Write to stencil buffer
			glStencilOp( GL_KEEP, GL_KEEP, GL_REPLACE )	# Only pixels that pass the depth test are written to the stencil buffer
			glPolygonMode( GL_FRONT_AND_BACK, GL_FILL )	
			self.renderSphere()
		
			# Then render the outline
			glStencilFunc( GL_NOTEQUAL, 1, 0xFFFF )		# The object itself is stenciled out
			glStencilOp( GL_KEEP, GL_KEEP, GL_REPLACE )
			glLineWidth( 4.0 )				# By increasing the line width only the outline is drawn
			glPolygonMode( GL_FRONT_AND_BACK, GL_LINE )
			glMaterialfv(GL_FRONT, GL_EMISSION, [0.0, 1.0, 0.0, 1.0])
			self.renderSphere()
	
			glPopAttrib()
		else:
			self.renderSphere()
			
	def renderSphere(self):
		# Material properties of the sphere
		glMaterialfv(GL_FRONT, GL_DIFFUSE, self.diffuse)
		glMaterialfv(GL_FRONT, GL_SPECULAR, self.specular)
		glMaterialf(GL_FRONT, GL_SHININESS, self.shininess)
		glMaterialfv(GL_FRONT, GL_AMBIENT, self.ambient)
		
		quadratic = gluNewQuadric()
		gluQuadricNormals(quadratic, GLU_SMOOTH)    # Create Smooth Normals (NEW) 
		gluQuadricTexture(quadratic, GL_TRUE)      # Create Texture Coords (NEW) 

		gluSphere(quadratic,self.radius,self.slices,self.stacks)
		
class EMCylinder(EMItem3D):
	name = "Cylinder"
	nodetype = "ShapeNode"
	def __init__(self, radius, height, transform=Transform()):
		EMItem3D.__init__(self, parent=None, children=set(), transform=transform)
		#size
		self.setRadiusAndHeight(radius, height)
		
		# color
		self.diffuse = [0.5,0.5,0.5,1.0]
		self.specular = [1.0,1.0,1.0,1.0]
		self.ambient = [1.0, 1.0, 1.0, 1.0]		
		self.shininess = 25.0
	
	def setRadiusAndHeight(self, radius, height):
		self.radius = radius
		self.height = height
		self.slices = int(radius)
		self.stacks = int(radius)
		self.boundingboxsize = max(self.radius, self.height)
		if self.item_inspector: self.item_inspector.updateMetaData()
		
	def getEvalString(self):
		return "EMCylinder(%s, %s)"%(self.radius, self.height)
		
	def setAmbientColor(self, red, green, blue, alpha=1.0):
		self.ambient = [red, green, blue, alpha]

	def setDiffuseColor(self, red, green, blue, alpha=1.0):
		self.diffuse = [red, green, blue, alpha]
		
	def setSpecularColor(self, red, green, blue, alpha=1.0):
		self.specular = [red, green, blue, alpha]
	
	def setShininess(self, shininess):
		self.shininess = shininess
		
	def getItemInspector(self):
		"""
		Return a Qt widget that controls the scene item
		"""
		if not self.item_inspector: self.item_inspector = EMInspectorControlShape("CYLINDER", self)
		return self.item_inspector
	
	def renderNode(self):
		if self.is_selected and glGetIntegerv(GL_RENDER_MODE) == GL_RENDER:
			glPushAttrib( GL_ALL_ATTRIB_BITS )
		
			# First render the cylinder, writing the outline to the stencil buffer
			glClearStencil(0)
			glClear( GL_STENCIL_BUFFER_BIT )
			glEnable( GL_STENCIL_TEST )
			glStencilFunc( GL_ALWAYS, 1, 0xFFFF )		# Write to stencil buffer
			glStencilOp( GL_KEEP, GL_KEEP, GL_REPLACE )	# Only pixels that pass the depth test are written to the stencil buffer
			glPolygonMode( GL_FRONT_AND_BACK, GL_FILL )	
			self.renderCylinder()
		
			# Then render the outline
			glStencilFunc( GL_NOTEQUAL, 1, 0xFFFF )		# The object itself is stenciled out
			glStencilOp( GL_KEEP, GL_KEEP, GL_REPLACE )
			glLineWidth( 4.0 )				# By increasing the line width only the outline is drawn
			glPolygonMode( GL_FRONT_AND_BACK, GL_LINE )
			glMaterialfv(GL_FRONT, GL_EMISSION, [0.0, 1.0, 0.0, 1.0])
			self.renderCylinder()
	
			glPopAttrib()
		else:
			self.renderCylinder()
	
	def renderCylinder(self):
	#def renderNode(self):	
		# Material properties of the cylinder
		glMaterialfv(GL_FRONT, GL_DIFFUSE, self.diffuse)
		glMaterialfv(GL_FRONT, GL_SPECULAR, self.specular)
		glMaterialf(GL_FRONT, GL_SHININESS, self.shininess)
		glMaterialfv(GL_FRONT, GL_AMBIENT, self.ambient)
		
		quadratic = gluNewQuadric()
		gluQuadricDrawStyle(quadratic, GLU_FILL)
		gluQuadricNormals(quadratic, GLU_SMOOTH)    # Create Smooth Normals (NEW) 
		gluQuadricTexture(quadratic, GL_TRUE)      # Create Texture Coords (NEW)
		gluCylinder(quadratic,self.radius,self.radius,self.height,self.slices,self.stacks)
		gluQuadricOrientation(quadratic,GLU_INSIDE)
		gluDisk( quadratic, 0.0, self.radius, self.slices, 1)
		glPushMatrix()
		glTranslatef( 0,0,self.height)
		gluQuadricOrientation(quadratic,GLU_OUTSIDE)
		gluDisk( quadratic, 0.0, self.radius, self.slices, 1)
		glPopMatrix()

class EMLine(EMItem3D):
	name = "Line"
	nodetype = "ShapeNode"
	pass

class EMCone(EMItem3D):
	name = "Cone"
	nodetype = "ShapeNode"
	def __init__(self, radius, height, transform=Transform()):
		EMItem3D.__init__(self, parent=None, children=set(), transform=transform)
		#size
		self.setRadiusAndHeight(radius, height)
		
		# color
		self.diffuse = [0.5,0.5,0.5,1.0]
		self.specular = [1.0,1.0,1.0,1.0]
		self.ambient = [1.0, 1.0, 1.0, 1.0]		
		self.shininess = 25.0
	
	def setRadiusAndHeight(self, radius, height):
		self.radius = radius
		self.height = height
		self.slices = int(radius)
		self.stacks = int(radius)
		self.boundingboxsize = max(self.radius, self.height)
		if self.item_inspector: self.item_inspector.updateMetaData()
		
	def getEvalString(self):
		return "EMCone(%s, %s)"%(self.radius, self.height)
		
	def setAmbientColor(self, red, green, blue, alpha=1.0):
		self.ambient = [red, green, blue, alpha]
	
	def setDiffuseColor(self, red, green, blue, alpha=1.0):
		self.diffuse = [red, green, blue, alpha]
		
	def setSpecularColor(self, red, green, blue, alpha=1.0):
		self.specular = [red, green, blue, alpha]
	
	def setShininess(self, shininess):
		self.shininess = shininess
		
	def getItemInspector(self):
		"""
		Return a Qt widget that controls the scene item
		"""
		if not self.item_inspector: self.item_inspector = EMInspectorControlShape("CONE", self)
		return self.item_inspector
	
	def renderNode(self):
		if self.is_selected:
			glPushAttrib( GL_ALL_ATTRIB_BITS )
		
			# First render the cone, writing the outline to the stencil buffer
			glClearStencil(0)
			glClear( GL_STENCIL_BUFFER_BIT )
			glEnable( GL_STENCIL_TEST )
			glStencilFunc( GL_ALWAYS, 1, 0xFFFF )		# Write to stencil buffer
			glStencilOp( GL_KEEP, GL_KEEP, GL_REPLACE )	# Only pixels that pass the depth test are written to the stencil buffer
			glPolygonMode( GL_FRONT_AND_BACK, GL_FILL )	
			self.renderCone()
		
			# Then render the outline
			glStencilFunc( GL_NOTEQUAL, 1, 0xFFFF )		# The object itself is stenciled out
			glStencilOp( GL_KEEP, GL_KEEP, GL_REPLACE )
			glLineWidth( 4.0 )				# By increasing the line width only the outline is drawn
			glPolygonMode( GL_FRONT_AND_BACK, GL_LINE )
			glMaterialfv(GL_FRONT, GL_EMISSION, [0.0, 1.0, 0.0, 1.0])
			self.renderCone()
	
			glPopAttrib()
		else:
			self.renderCone()
	
	def renderCone(self):
	#def renderNode(self):	
		# Material properties of the cone
		glMaterialfv(GL_FRONT, GL_DIFFUSE, self.diffuse)
		glMaterialfv(GL_FRONT, GL_SPECULAR, self.specular)
		glMaterialf(GL_FRONT, GL_SHININESS, self.shininess)
		glMaterialfv(GL_FRONT, GL_AMBIENT, self.ambient)
		
		quadratic = gluNewQuadric()
		gluQuadricDrawStyle(quadratic, GLU_FILL)
		gluQuadricNormals(quadratic, GLU_SMOOTH)    # Create Smooth Normals (NEW) 
		gluQuadricTexture(quadratic, GL_TRUE)      # Create Texture Coords (NEW)
		gluCylinder(quadratic,0,self.radius,self.height,self.slices,self.stacks)
		glPushMatrix()
		glTranslatef( 0,0,self.height)
		gluQuadricOrientation(quadratic,GLU_OUTSIDE)
		gluDisk( quadratic, 0.0, self.radius, self.slices, 1)
		glPopMatrix()

class EM3DText(EMItem3D):
	name = "3DText"
	nodetype = "ShapeNode"
	pass

class EMInspectorControlShape(EMItem3DInspector):
	"""
	Class to make EMItem GUI SHAPE Inspector
	"""
	def __init__(self, name, item3d):
		EMItem3DInspector.__init__(self, name, item3d)
	
	def updateItemControls(self):
		super(EMInspectorControlShape, self).updateItemControls()
		self.ambcolorbox.setColor(QtGui.QColor(255*self.item3d().ambient[0],255*self.item3d().ambient[1],255*self.item3d().ambient[2]))
		self.diffusecolorbox.setColor(QtGui.QColor(255*self.item3d().diffuse[0],255*self.item3d().diffuse[1],255*self.item3d().diffuse[2]))
		self.specularcolorbox.setColor(QtGui.QColor(255*self.item3d().specular[0],255*self.item3d().specular[1],255*self.item3d().specular[2]))
		
	def addControls(self, gridbox):
		super(EMInspectorControlShape, self).addControls(gridbox)
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
		cabox.addWidget(self.ambcolorbox)
		cabox.setAlignment(QtCore.Qt.AlignCenter)
		cdbox = QtGui.QHBoxLayout()
		self.diffusecolorbox = EMQTColorWidget(parent=colorframe)
		cdbox.addWidget(self.diffusecolorbox)
		cdbox.setAlignment(QtCore.Qt.AlignCenter)
		csbox = QtGui.QHBoxLayout()
		self.specularcolorbox = EMQTColorWidget(parent=colorframe)
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
		
		self.shininess = ValSlider(colorframe, (0.0, 50.0), "Shine")
		self.shininess.setValue(self.item3d().shininess)
		
		colorvbox.addWidget(colorlabel)
		colorvbox.addLayout(cdialoghbox)
		colorvbox.addLayout(colorhbox)
		colorvbox.addWidget(self.shininess)
		colorframe.setLayout(colorvbox)
		colorframe.setMaximumWidth(350)
		gridbox.addWidget(colorframe, 3, 0, 1, 1)
		
		# Set to default, but do not run if being inherited
		if type(self) == EMInspectorControlShape: self.updateItemControls()
		
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