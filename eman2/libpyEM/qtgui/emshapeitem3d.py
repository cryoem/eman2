#!/usr/bin/env python
#
# Author: Grant Tang (gtang@bcm.edu)
# Author: John Flanagan (jfflanag@bcm.edu)
# Copyright (c) 2011- Baylor College of Medicine


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
from emglobjects import init_glut, get_default_gl_colors
from emitem3d import EMItem3D, EMItem3DInspector, drawBoundingBox
from libpyGLUtils2 import FTGLFontMode
import math
from valslider import EMQTColorWidget, ValSlider, EMSpinWidget

from OpenGL.GL import *
from OpenGL.GLU import *
from PyQt4 import QtCore, QtGui, QtOpenGL

import numpy as np

class EMShapeBase(EMItem3D):
	""" Base class that all shape object need to inherit """
	name = "shapebase"
	nodetype = "ShapeNode"
	
	def __init__(self, parent=None, children=None, transform=None):
		if not transform: transform = Transform()
		EMItem3D.__init__(self, parent=parent, children=children, transform=transform)
	
		# initial color
		self.diffuse = [0.5,0.5,0.5,1.0]
		self.specular = [1.0,1.0,1.0,1.0]
		self.ambient = [1.0, 1.0, 1.0, 1.0]
		self.shininess = 25.0		
		
	def getEvalString(self):
		raise NotImplementedError("Not to reimplemnt this function")
	
	def getItemInspector(self):
		raise NotImplementedError("Not to reimplemnt this function")
	
	def setAmbientColor(self, red, green, blue, alpha=1.0):
		self.ambient = [red, green, blue, alpha]

	def setDiffuseColor(self, red, green, blue, alpha=1.0):
		self.diffuse = [red, green, blue, alpha]
		
	def setSpecularColor(self, red, green, blue, alpha=1.0):
		self.specular = [red, green, blue, alpha]
	
	def setShininess(self, shininess):
		self.shininess = shininess
	
	def getItemDictionary(self):
		"""
		Return a dictionary of item parameters (used for restoring sessions
		"""
		dictionary = super(EMShapeBase, self).getItemDictionary()
		dictionary.update({"COLOR":[self.ambient, self.diffuse, self.specular, self.shininess]})
		return dictionary
	
	def setUsingDictionary(self, dictionary):
		"""
		Set item attributes using a dictionary, used in session restoration
		"""
		super(EMShapeBase, self).setUsingDictionary(dictionary)
		self.setAmbientColor(*dictionary["COLOR"][0])
		self.setDiffuseColor(*dictionary["COLOR"][1])
		self.setSpecularColor(*dictionary["COLOR"][2])
		self.setShininess(dictionary["COLOR"][3])
		
	def renderNode(self):
		if self.is_selected and glGetIntegerv(GL_RENDER_MODE) == GL_RENDER and not self.isSelectionHidded(): # No need to draw outline in selection mode
			#if glGetIntegerv(GL_RENDER_MODE) == GL_RENDER: print "X"
			glPushAttrib( GL_ALL_ATTRIB_BITS )
		
			# First render the cylinder, writing the outline to the stencil buffer
			glClearStencil(0)
			glClear( GL_STENCIL_BUFFER_BIT )
			glEnable( GL_STENCIL_TEST )
			glStencilFunc( GL_ALWAYS, 1, 0xFFFF )		# Write to stencil buffer
			glStencilOp( GL_KEEP, GL_KEEP, GL_REPLACE )	# Only pixels that pass the depth test are written to the stencil buffer
			glPolygonMode( GL_FRONT_AND_BACK, GL_FILL )	
			self.renderShape()
		
			# Then render the outline
			glStencilFunc( GL_NOTEQUAL, 1, 0xFFFF )		# The object itself is stenciled out
			glStencilOp( GL_KEEP, GL_KEEP, GL_REPLACE )
			glLineWidth( 4.0 )				# By increasing the line width only the outline is drawn
			glPolygonMode( GL_FRONT_AND_BACK, GL_LINE )
			glMaterialfv(GL_FRONT, GL_EMISSION, [0.0, 1.0, 0.0, 1.0])
			self.renderShape()
	
			glPopAttrib()	
		else:
			glPushAttrib( GL_ALL_ATTRIB_BITS )
			
			self.renderShape()
			
			glPopAttrib()
			
	def renderShape(self):
		raise NotImplementedError("Not to reimplemnt this function")	
	
class EMRuler(EMShapeBase):
	"""
	Widget to make a ruler for measuring distances and drawing a little ruler on the screen
	"""
	name = "Ruler"
	nodetype = "ShapeNode"
	
	@staticmethod
	def getNodeDialogWidget(attribdict):
		"""
		Return a cube control widget for the stacked_widget
		"""
		raise NotImplementedError("Not yet implmented")
	
	@staticmethod
	def getNodeForDialog(attribdict):
		"""
		Create a new node using a attribdict
		"""
		raise NotImplementedError("Not yet implmented")
			
	def __init__(self, x1, y1, z1, x2, y2, z2, apix, scaling, transform=None):
		EMShapeBase.__init__(self, parent=None, children=set(), transform=transform)
		self.setRulerAPix(apix)
		self.setRulerScaling(scaling)
		self.barwidth = 10.0
		self.setRuler(x1, y1, z1, x2, y2, z2)
		
	def getEvalString(self):
		return "EMRuler(%s,%s,%s,%s,%s,%s,%s,%s)"%(self.xi,self.yi,self.zi,self.xf,self.yf,self.zf,self.getRulerAPix(),self.getRulerScaling())	
		
	def updateInfo(self):
		""" 
		Sets the ruler length 
		"""
		scaledapix = self.getRulerAPix()*self.getRulerScaling()
		length = self.getLength()*scaledapix
		smallbar = 2*scaledapix*self.barwidth
		self.boundingboxsize = 'length='+str(round(length, 2))+u'\u212B'+';  smallbar='+str(round(smallbar, 2))+u'\u212B'+';  apix='+str(round(self.getRulerAPix(), 2))+u'\u212B'
		if self.item_inspector: self.item_inspector.updateMetaData()
		
	def getRulerAPix(self):
		"""Return the apix for this ruler"""
		return self.apix

	def setRulerAPix(self, apix):
		"""Set the apix for this ruler, overrides SG value"""
		self.apix = apix
	
	def getRulerScaling(self):
		""" Set ruler scaling, overrides SG value """
		return self.scaling
		
	def setRulerScaling(self, scaling):
		""" Get ruler scaling """
		self.scaling = scaling
		
	def getLength(self):
		"""Return the length"""
		return self.pixlen
		
	def setLength(self, length):
		"""Set the length in pixels"""
		self.pixlen = length
		
	def getItemInspector(self):
		"""Return a Qt widget that controls the scene item"""
		if not self.item_inspector: self.item_inspector = EMInspectorControlShape("RULER", self)
		return self.item_inspector
	
	def setRuler(self, x1, y1, z1, x2, y2, z2):
		"""Set the ruler based on end and starting points"""
		self.xi = x1
		self.yi = y1
		self.zi = z1
		self.xf = x2
		self.yf = y2
		self.zf = z2
		self.setLength(math.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2))
		self.updateInfo()
		# Compute bars using parametric equations, and vector
		angle = -math.atan2((y2-y1),(x2-x1))
		if self.getLength() > 0:
			self.direction = [(x2-x1)/self.getLength(), (y2-y1)/self.getLength(), (z2-z1)/self.getLength()]
		else:
			self.direction = [0.0,0.0,0.0]
		self.rsinO = self.barwidth*math.sin(angle)
		self.rcosO = self.barwidth*math.cos(angle)
		self.smallbars = [[i*self.direction[0],i*self.direction[1]] for i in xrange(0,int(self.pixlen),int(2*self.barwidth))]
		
	def renderShape(self):        
		# Material properties of the box
		glDisable(GL_COLOR_MATERIAL)
		glMaterialfv(GL_FRONT, GL_DIFFUSE, self.diffuse)
		glMaterialfv(GL_FRONT, GL_SPECULAR, self.specular)
		glMaterialf(GL_FRONT, GL_SHININESS, self.shininess)
		glMaterialfv(GL_FRONT, GL_AMBIENT, self.ambient)
		
		glEnable(GL_BLEND)
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
		#glLineWidth(1.5)
		glEnable(GL_LINE_SMOOTH)
	
		glBegin(GL_LINES)
		glVertex3f(self.xi, self.yi, self.zi)
		glVertex3f(self.xf, self.yf, self.zf)
		glVertex3f(self.xi+self.rsinO, self.yi+self.rcosO, self.zf)
		glVertex3f(self.xi-self.rsinO, self.yi-self.rcosO, self.zf)
		glVertex3f(self.xf+self.rsinO, self.yf+self.rcosO, self.zf)
		glVertex3f(self.xf-self.rsinO, self.yf-self.rcosO, self.zf)
		for i in self.smallbars:
			glVertex3f(self.xi+i[0]+self.rsinO/2.0, self.yi+i[1]+self.rcosO/2.0, self.zf)
			glVertex3f(self.xi+i[0]-self.rsinO/2.0, self.yi+i[1]-self.rcosO/2.0, self.zf)
		glEnd()	# Done Drawing The Cube
		
class EMCube(EMShapeBase):
	name = "Cube"
	nodetype = "ShapeNode"
	
	@staticmethod
	def getNodeDialogWidget(attribdict):
		"""
		Return a cube control widget for the stacked_widget
		"""
		cubewidget = QtGui.QWidget()
		grid = QtGui.QGridLayout()
		cube_dim_label = QtGui.QLabel("Cube Dimension")
		attribdict["cube_dim"] = QtGui.QLineEdit("50")
		node_name_label = QtGui.QLabel("Cube Name")
		attribdict["node_name"] = QtGui.QLineEdit(str(EMCube.name))
		grid.addWidget(cube_dim_label, 0, 0, 1, 2)
		grid.addWidget(attribdict["cube_dim"], 0, 2, 1, 2)
		grid.addWidget(node_name_label , 1, 0, 1, 2)
		grid.addWidget(attribdict["node_name"], 1, 2, 1, 2)
		EMItem3D.get_transformlayout(grid, 2, attribdict)
		cubewidget.setLayout(grid)
		return cubewidget
	
	@staticmethod
	def getNodeForDialog(attribdict):
		"""
		Create a new node using a attribdict
		"""
		return EMCube(float(attribdict["cube_dim"].text()), transform=EMItem3D.getTransformFromDict(attribdict))
		
	def __init__(self, size, transform=None):
		EMShapeBase.__init__(self, parent=None, children=set(), transform=transform)
		# size
		self.setSize(size)
		
	def setSize(self, size):
		self.size = size
		self.boundingboxsize = str(round(size, 2))+u'\u00B3'
		self.xi = -size/2
		self.yi = -size/2
		self.zi = -size/2
		self.xf = size/2
		self.yf = size/2
		self.zf = size/2
		if self.item_inspector: self.item_inspector.updateMetaData()
		
	def getEvalString(self):
		return "EMCube(%s)"%self.size
		
	def getItemInspector(self):
		"""
		Return a Qt widget that controls the scene item
		"""
		if not self.item_inspector: self.item_inspector = EMInspectorControlShape("CUBE", self)
		return self.item_inspector
			
	def renderShape(self):        
		# Material properties of the box
		
		glDisable(GL_COLOR_MATERIAL)
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

class EMSphere(EMShapeBase):
	name = "Sphere"
	nodetype = "ShapeNode"
	
	@staticmethod
	def getNodeDialogWidget(attribdict):
		"""
		Return a sphere control widget for the stacked_widget
		"""
		spherewidget = QtGui.QWidget()
		grid = QtGui.QGridLayout()
		sphere_dim_label = QtGui.QLabel("Sphere Dimension")
		attribdict["sphere_dim"] = QtGui.QLineEdit("50")
		node_name_label = QtGui.QLabel("Sphere Name")
		attribdict["node_name"] = QtGui.QLineEdit(str(EMSphere.name))
		grid.addWidget(sphere_dim_label, 0, 0, 1, 2)
		grid.addWidget(attribdict["sphere_dim"], 0, 2, 1, 2)
		grid.addWidget(node_name_label , 1, 0, 1, 2)
		grid.addWidget(attribdict["node_name"], 1, 2, 1, 2)
		EMItem3D.get_transformlayout(grid, 2, attribdict)
		spherewidget.setLayout(grid)
		return spherewidget
	
	@staticmethod
	def getNodeForDialog(attribdict):
		"""
		Create a new node using a attribdict
		"""
		return EMSphere(float(attribdict["sphere_dim"].text()), transform=EMItem3D.getTransformFromDict(attribdict))
		
	def __init__(self, radius, transform=None):
		EMShapeBase.__init__(self, parent=None, children=set(), transform=transform)
		# size
		self.setRadius(radius)
	
	def setRadius(self, radius):
		self.radius = radius
		self.boundingboxsize = str(round(radius, 2))+u'\u00B3'
		self.slices = int(radius)
		self.stacks = int(radius)
		if self.item_inspector: self.item_inspector.updateMetaData()
		
	def getEvalString(self):
		return "EMSphere(%s)"%self.radius		
		
	def getItemInspector(self):
		"""
		Return a Qt widget that controls the scene item
		"""
		if not self.item_inspector: self.item_inspector = EMInspectorControlShape("SPHERE", self)
		return self.item_inspector
			
	def renderShape(self):
		# Material properties of the sphere
		glDisable(GL_COLOR_MATERIAL)
		glMaterialfv(GL_FRONT, GL_DIFFUSE, self.diffuse)
		glMaterialfv(GL_FRONT, GL_SPECULAR, self.specular)
		glMaterialf(GL_FRONT, GL_SHININESS, self.shininess)
		glMaterialfv(GL_FRONT, GL_AMBIENT, self.ambient)
		
		quadratic = gluNewQuadric()
		gluQuadricNormals(quadratic, GLU_SMOOTH)    # Create Smooth Normals (NEW) 
		gluQuadricTexture(quadratic, GL_TRUE)      # Create Texture Coords (NEW) 

		gluSphere(quadratic,self.radius,self.slices,self.stacks)

class EMScatterPlot3D(EMShapeBase):
	name = "Plot3D"
	nodetype = "ShapeNode"
	
	def __init__(self, transform=None):
		EMShapeBase.__init__(self, parent=None, children=set(), transform=transform)
		
	def setData(self, data, pointsize):
		""" Set the dat to plot. Format is a [X, Y, Z] whereX Y and Z are lists of the smae length """
		if len(data) != 3:
			raise ValueError("Data must have 3 dimensions")
		if len(data[0]) != len(data[1]) or len(data[1]) != len(data[2]):
			raise ValueError("Dimensions must be of the same length")
		self.data = data
		self.setPointSize(pointsize)

	def getPointSize(self):
		return self.pointsize
		
	def setPointSize(self, pointsize):
		self.pointsize = pointsize
		self.slices = int(pointsize)
		self.stacks = int(pointsize)
		
	def getEvalString(self):
		"""Implement me"""
		raise NotImplementedError("Implement me")
	
	def getItemInspector(self):
		""" Implement me """
		if not self.item_inspector: self.item_inspector = EMInspectorControlScatterPlot("Scatter Plot", self)
		return self.item_inspector
	
	def renderShape(self):
		# Material properties of the sphere
		glDisable(GL_COLOR_MATERIAL)
		glMaterialfv(GL_FRONT, GL_DIFFUSE, self.diffuse)
		glMaterialfv(GL_FRONT, GL_SPECULAR, self.specular)
		glMaterialf(GL_FRONT, GL_SHININESS, self.shininess)
		glMaterialfv(GL_FRONT, GL_AMBIENT, self.ambient)
		
		quadratic = gluNewQuadric()
		gluQuadricNormals(quadratic, GLU_SMOOTH)    # Create Smooth Normals (NEW) 
		gluQuadricTexture(quadratic, GL_TRUE)      # Create Texture Coords (NEW)
		
		for i in xrange(len(self.data[0])):
			glPushMatrix()
			
			glTranslatef(self.data[0][i],self.data[1][i],self.data[2][i])
			gluSphere(quadratic,self.pointsize,self.slices,self.stacks)
			
			glPopMatrix()

			
class EMCylinder(EMShapeBase):
	name = "Cylinder"
	nodetype = "ShapeNode"
	
	@staticmethod
	def getNodeDialogWidget(attribdict):
		"""
		Return a cylinder control widget for the stacked_widget
		"""
		cyliderwidget = QtGui.QWidget()
		grid = QtGui.QGridLayout()
		cylider_radius_label = QtGui.QLabel("Cylider Radius")
		attribdict["cylider_radius"] = QtGui.QLineEdit("50")
		grid.addWidget(cylider_radius_label, 0, 0, 1, 2)
		grid.addWidget(attribdict["cylider_radius"], 0, 2, 1, 2)
		cylider_height_label = QtGui.QLabel("Cylider Height")
		attribdict["cylider_height"] = QtGui.QLineEdit("50")
		node_name_label = QtGui.QLabel("Cylider Name")
		attribdict["node_name"] = QtGui.QLineEdit(str(EMCylinder.name))
		grid.addWidget(cylider_height_label, 1, 0, 1, 2)
		grid.addWidget(attribdict["cylider_height"], 1, 2, 1, 2)
		grid.addWidget(node_name_label , 2, 0, 1, 2)
		grid.addWidget(attribdict["node_name"], 2, 2, 1, 2)
		EMItem3D.get_transformlayout(grid, 4, attribdict)
		cyliderwidget.setLayout(grid)
		return cyliderwidget
	
	@staticmethod
	def getNodeForDialog(attribdict):
		"""
		Create a new node using a attribdict
		"""
		return EMCylinder(float(attribdict["cylider_radius"].text()), float(attribdict["cylider_height"].text()), transform=EMItem3D.getTransformFromDict(attribdict))
		
	def __init__(self, radius, height, transform=None):
		EMShapeBase.__init__(self, parent=None, children=set(), transform=transform)
		#size
		self.setRadiusAndHeight(radius, height)
		
	def setRadiusAndHeight(self, radius, height):
		self.radius = radius
		self.height = height
		self.slices = int(radius)
		self.stacks = int(radius)
		self.boundingboxsize = str(round(radius, 2))+'x'+str(round(height, 2))
		if self.item_inspector: self.item_inspector.updateMetaData()
		
	def getEvalString(self):
		return "EMCylinder(%s, %s)"%(self.radius, self.height)
		
	def getItemInspector(self):
		"""
		Return a Qt widget that controls the scene item
		"""
		if not self.item_inspector: self.item_inspector = EMInspectorControlShape("CYLINDER", self)
		return self.item_inspector
	
	def renderShape(self):
		# Material properties of the cylinder
		glDisable(GL_COLOR_MATERIAL)
		glMaterialfv(GL_FRONT, GL_DIFFUSE, self.diffuse)
		glMaterialfv(GL_FRONT, GL_SPECULAR, self.specular)
		glMaterialf(GL_FRONT, GL_SHININESS, self.shininess)
		glMaterialfv(GL_FRONT, GL_AMBIENT, self.ambient)
		
		glPushMatrix()
		quadratic = gluNewQuadric()
		gluQuadricDrawStyle(quadratic, GLU_FILL)
		gluQuadricNormals(quadratic, GLU_SMOOTH)    # Create Smooth Normals (NEW) 
		gluQuadricTexture(quadratic, GL_TRUE)      # Create Texture Coords (NEW)
		glTranslatef( 0,0,-self.height/2)
		gluCylinder(quadratic,self.radius,self.radius,self.height,self.slices,self.stacks)
		gluQuadricOrientation(quadratic,GLU_INSIDE)
		gluDisk( quadratic, 0.0, self.radius, self.slices, 1)
		glPushMatrix()
		glTranslatef( 0,0,self.height)
		gluQuadricOrientation(quadratic,GLU_OUTSIDE)
		gluDisk( quadratic, 0.0, self.radius, self.slices, 1)
		glPopMatrix()
		glPopMatrix()

class EMLine(EMShapeBase):
	name = "Line"
	nodetype = "ShapeNode"
	
	@staticmethod
	def getNodeDialogWidget(attribdict):
		"""
		Return a line control widget for the stacked_widget
		"""
		linewidget = QtGui.QWidget()
		grid = QtGui.QGridLayout()
		line_xyzi_label = QtGui.QLabel("Line start, X, Y, Z")
		attribdict["linexi"] = QtGui.QLineEdit("0.0")
		attribdict["lineyi"] = QtGui.QLineEdit("0.0")
		attribdict["linezi"] = QtGui.QLineEdit("0.0")
		grid.addWidget(line_xyzi_label, 0, 0, 1, 3)
		grid.addWidget(attribdict["linexi"], 1, 0, 1, 1)
		grid.addWidget(attribdict["lineyi"], 1, 1, 1, 1)
		grid.addWidget(attribdict["linezi"], 1, 2, 1, 1)
		line_xyzf_label = QtGui.QLabel("Line end, X, Y, Z")
		attribdict["linexf"] = QtGui.QLineEdit("0.0")
		attribdict["lineyf"] = QtGui.QLineEdit("0.0")
		attribdict["linezf"] = QtGui.QLineEdit("0.0")
		grid.addWidget(line_xyzf_label, 2, 0, 1, 3)
		grid.addWidget(attribdict["linexf"], 3, 0, 1, 1)
		grid.addWidget(attribdict["lineyf"], 3, 1, 1, 1)
		grid.addWidget(attribdict["linezf"], 3, 2, 1, 1)
		line_width = QtGui.QLabel("Line Width")
		line_width.setAlignment(QtCore.Qt.AlignCenter)
		attribdict["linewidth"] = QtGui.QLineEdit("10.0")
		grid.addWidget(line_width, 4, 0, 1, 2)
		grid.addWidget(attribdict["linewidth"], 4, 2, 1, 1)
		node_name_label = QtGui.QLabel("Line Name")
		attribdict["node_name"] = QtGui.QLineEdit(str(EMLine.name))
		grid.addWidget(node_name_label , 5, 0, 1, 3)
		grid.addWidget(attribdict["node_name"], 6, 0, 1, 3)
		linewidget.setLayout(grid)
		return linewidget
	
	@staticmethod
	def getNodeForDialog(attribdict):
		"""
		Create a new node using a attribdict
		"""
		transform=Transform()
		return EMLine(float(attribdict["linexi"].text()), float(attribdict["lineyi"].text()), float(attribdict["linezi"].text()), float(attribdict["linexf"].text()), float(attribdict["lineyf"].text()), float(attribdict["linezf"].text()), float(attribdict["linewidth"].text()), transform=transform)
		
	def __init__(self, x1, y1, z1, x2, y2, z2, linewidth, transform=None):
		EMShapeBase.__init__(self, parent=None, children=set(), transform=transform)
		# size
		self.setEndAndWidth(x1, y1, z1, x2, y2, z2, linewidth)
		
		self.setShowLeftArrow(True)
		self.setShowRightArrow(True)

	def setEndAndWidth(self, x1, y1, z1, x2, y2, z2, width, comparrows=True):
		self.x1 = x1
		self.y1 = y1
		self.z1 = z1
		self.x2 = x2
		self.y2 = y2
		self.z2 = z2
		self.width = width
		self.slices = int(width/2)
		self.stacks = int(width/2)
		
		#arrow size
		dx = self.x1 - self.x2
		dy = self.y1 - self.y2
		dz = self.z1 - self.z2
		self.length = math.sqrt(dx*dx + dy*dy + dz*dz)	#cylinder length
		if comparrows: 
			self.leftArrowSize = self.width
			self.leftArrowLength = self.length/10.0
			self.rightArrowSize = self.width
			self.rightArrowLength = self.length/10.0
		self.boundingboxsize = 'length='+str(round(self.length, 2))+', width='+str(round(self.width, 2))
		
		if self.item_inspector: self.item_inspector.updateMetaData()	
	
	def setLength(self, length):
		v = Vec3f((self.x2-self.x1), (self.y2-self.y1), (self.z2-self.z1))
		v.normalize()
		self.setEndAndWidth(self.x1,self.y1,self.z1,v[0]*length,v[1]*length,v[2]*length,self.width,comparrows=False)
		
	def setSlices(self, slices):
		if slices>0:
			self.slices = int(slices)
		else:
			self.slices = int(self.width/2)
	
	def setWidth(self, width):
		self.width = width
		
	def setStacks(self, stacks):
		if stacks>0:
			self.stacks = int(stacks)
		else:
			self.stacks = int(self.width/2)
	
	def getEvalString(self):
		return "EMLine(%s, %s, %s, %s, %s, %s, %s)"%(self.x1, self.y1, self.z1, self.x2, self.y2,self.z2, self.width)

	def setShowLeftArrow(self, state):
		self.showLeftArrow = state

	def setShowRightArrow(self, state):
		self.showRightArrow = state
		
	def setLeftArrowSize(self, size):
		self.leftArrowSize = size
		
	def setLeftArrowLength(self, length):
		self.leftArrowLength = length

	def setRightArrowSize(self, size):
		self.rightArrowSize = size
		
	def setRightArrowLength(self, length):
		self.rightArrowLength = length
		
	def getItemInspector(self):
		"""
		Return a Qt widget that controls the scene item
		"""
		if not self.item_inspector: self.item_inspector = EMInspectorControlLine("LINE", self)
		return self.item_inspector
		
	def getItemDictionary(self):
		"""
		Return a dictionary of item parameters (used for restoring sessions
		"""
		dictionary = super(EMLine, self).getItemDictionary()
		dictionary.update({"LINEPARS":[self.leftArrowSize, self.leftArrowLength, self.showLeftArrow, self.rightArrowSize, self.rightArrowLength, self.showRightArrow, self.slices, self.stacks, self.width]})
		return dictionary
		
	def setUsingDictionary(self, dictionary):
		"""
		Set item attributes using a dictionary, used in session restoration
		"""
		super(EMLine, self).setUsingDictionary(dictionary)
		self.leftArrowSize = dictionary["LINEPARS"][0]
		self.leftArrowLength = dictionary["LINEPARS"][1]
		self.showLeftArrow = dictionary["LINEPARS"][2]
		self.rightArrowSize = dictionary["LINEPARS"][3]
		self.rightArrowLength = dictionary["LINEPARS"][4]
		self.showRightArrow = dictionary["LINEPARS"][5]
		self.setSlices(dictionary["LINEPARS"][6])
		self.setStacks(dictionary["LINEPARS"][7])
		self.setWidth(dictionary["LINEPARS"][8])	
			
	def renderShape(self):
		r2d = 180.0/math.pi
				
		glPushMatrix()
		glTranslatef(self.x1, self.y1, self.z1)
		
		#orentation vector
		vx = self.x2 - self.x1
		vy = self.y2 - self.y1
		vz = self.z2 - self.z1
		
		#rotation vector, z x r
		rx = -vy*vz
		ry = vx*vz
		ax = 0.0
		
		if vz == 0:
			ax = r2d*math.atan2(vy, vx) # Find trig, John F
		else:
			ax = r2d*math.acos(vz/self.length)
			if vz<=0: ax = -ax

		if vz==0:
			glRotated(90.0, 0, 1, 0.0)	#Rotate & align with x axis
			glRotated(ax, -1.0, 0.0, 0.0)	#Rotate to point 2 in x-y plane
		else:
			glRotated(ax, rx, ry, 0)
			
		glDisable(GL_COLOR_MATERIAL)
		glMaterialfv(GL_FRONT, GL_DIFFUSE, self.diffuse)
		glMaterialfv(GL_FRONT, GL_SPECULAR, self.specular)
		glMaterialf(GL_FRONT, GL_SHININESS, self.shininess)
		glMaterialfv(GL_FRONT, GL_AMBIENT, self.ambient)
		
		glPushMatrix()
		quadratic = gluNewQuadric()
		gluQuadricDrawStyle(quadratic, GLU_FILL)
		gluQuadricNormals(quadratic, GLU_SMOOTH)    # Create Smooth Normals (NEW) 
		gluQuadricTexture(quadratic, GL_TRUE)      # Create Texture Coords (NEW)
		#glTranslatef( 0,0,-self.length/2)			# Do we want the line to be centered on the origin?
		gluCylinder(quadratic,self.width/2,self.width/2,self.length,self.slices,self.stacks)
		gluQuadricOrientation(quadratic,GLU_OUTSIDE)
		glTranslatef( 0,0,-self.leftArrowLength)
		if self.showLeftArrow:
			gluCylinder(quadratic,0,self.leftArrowSize,self.leftArrowLength,self.slices,self.stacks)
		glTranslatef( 0,0,self.length+self.leftArrowLength)
		gluQuadricOrientation(quadratic,GLU_OUTSIDE)
		if self.showRightArrow:
			gluCylinder(quadratic,self.rightArrowSize,0,self.rightArrowLength,self.slices,self.stacks)
		glPopMatrix()
		glPopMatrix()

class EMCone(EMShapeBase):
	name = "Cone"
	nodetype = "ShapeNode"
	
	@staticmethod
	def getNodeDialogWidget(attribdict):
		"""
		Return a cone control widget for the stacked_widget
		"""
		conewidget = QtGui.QWidget()
		grid = QtGui.QGridLayout()
		cone_radius_label = QtGui.QLabel("Cone Radius")
		attribdict["cone_radius"] = QtGui.QLineEdit("50")
		grid.addWidget(cone_radius_label, 0, 0, 1, 2)
		grid.addWidget(attribdict["cone_radius"], 0, 2, 1, 2)
		cone_height_label = QtGui.QLabel("Cone Height")
		attribdict["cone_height"] = QtGui.QLineEdit("50")
		node_name_label = QtGui.QLabel("Cone Name")
		attribdict["node_name"] = QtGui.QLineEdit(str(EMCone.name))
		grid.addWidget(cone_height_label, 1, 0, 1, 2)
		grid.addWidget(attribdict["cone_height"], 1, 2, 1, 2)
		grid.addWidget(node_name_label , 2, 0, 1, 2)
		grid.addWidget(attribdict["node_name"], 2, 2, 1, 2)
		EMItem3D.get_transformlayout(grid, 4, attribdict)
		conewidget.setLayout(grid)
		return conewidget
	
	@staticmethod
	def getNodeForDialog(attribdict):
		"""
		Create a new node using a attribdict
		"""
		return EMCone(float(attribdict["cone_radius"].text()), float(attribdict["cone_height"].text()), transform=EMItem3D.getTransformFromDict(attribdict))
		
	def __init__(self, radius, height, transform=None):
		EMShapeBase.__init__(self, parent=None, children=set(), transform=transform)
		#size
		self.setRadiusAndHeight(radius, height)
		
	def setRadiusAndHeight(self, radius, height):
		self.radius = radius
		self.height = height
		self.slices = int(radius)
		self.stacks = int(radius)
		self.boundingboxsize = str(round(radius, 2))+'x'+str(round(height, 2))
		if self.item_inspector: self.item_inspector.updateMetaData()
		
	def getEvalString(self):
		return "EMCone(%s, %s)"%(self.radius, self.height)
		
	def getItemInspector(self):
		"""
		Return a Qt widget that controls the scene item
		"""
		if not self.item_inspector: self.item_inspector = EMInspectorControlShape("CONE", self)
		return self.item_inspector
	
	def renderShape(self):
		# Material properties of the cone
		glDisable(GL_COLOR_MATERIAL)
		glMaterialfv(GL_FRONT, GL_DIFFUSE, self.diffuse)
		glMaterialfv(GL_FRONT, GL_SPECULAR, self.specular)
		glMaterialf(GL_FRONT, GL_SHININESS, self.shininess)
		glMaterialfv(GL_FRONT, GL_AMBIENT, self.ambient)
		
		glPushMatrix()
		quadratic = gluNewQuadric()
		gluQuadricDrawStyle(quadratic, GLU_FILL)
		gluQuadricNormals(quadratic, GLU_SMOOTH)    # Create Smooth Normals (NEW) 
		gluQuadricTexture(quadratic, GL_TRUE)      # Create Texture Coords (NEW)
		glTranslatef( 0,0,-self.height/2)
		gluCylinder(quadratic,0,self.radius,self.height,self.slices,self.stacks)
		glPushMatrix()
		glTranslatef( 0,0,self.height)
		gluQuadricOrientation(quadratic,GLU_OUTSIDE)
		gluDisk( quadratic, 0.0, self.radius, self.slices, 1)
		glPopMatrix()
		glPopMatrix()

class EM3DText(EMShapeBase):
	name = "3DText"
	nodetype = "ShapeNode"
	
	@staticmethod
	def getNodeDialogWidget(attribdict):
		"""
		Return a text control widget for the stacked_widget
		"""
		textwidget = QtGui.QWidget()
		grid = QtGui.QGridLayout()
		text_label = QtGui.QLabel("Text")
		attribdict["text_content"] = QtGui.QLineEdit()
		grid.addWidget(text_label, 0, 0, 1, 2)
		grid.addWidget(attribdict["text_content"], 0, 2, 1, 2)
		fontsize_label = QtGui.QLabel("Font Size")
		attribdict["fontsize"] = QtGui.QLineEdit("32.0")
		grid.addWidget(fontsize_label , 1, 0, 1, 2)
		grid.addWidget(attribdict["fontsize"], 1, 2, 1, 2)
		node_name_label = QtGui.QLabel("Text Name")
		attribdict["node_name"] = QtGui.QLineEdit(str(EM3DText.name))
		grid.addWidget(node_name_label , 2, 0, 1, 2)
		grid.addWidget(attribdict["node_name"], 2, 2, 1, 2)
		EMItem3D.get_transformlayout(grid, 4, attribdict)
		textwidget.setLayout(grid)
		return textwidget
	
	@staticmethod
	def getNodeForDialog(attribdict):
		"""
		Create a new node using a attribdict
		"""
		return EM3DText(str(attribdict["text_content"].text()), float(attribdict["fontsize"].text()), transform=EMItem3D.getTransformFromDict(attribdict))
		
	def __init__(self, string, fontSize, fontMode=FTGLFontMode.TEXTURE, depth=10, transform=None):
		EMShapeBase.__init__(self, parent=None, children=set(), transform=transform)
		#size
		self.setRenderString(string, fontSize)
		self.setFontMode(fontMode)
		self.setFontDepth(depth)
		
		self.font_renderer = get_3d_font_renderer()
		
	def setFontMode(self, fontMode):
		self.fontMode = FTGLFontMode(fontMode)	# Cast to enum
		if self.item_inspector: self.item_inspector.updateMetaData()
		
	def setFontDepth(self, fontDepth):
		self.fontDepth = fontDepth
		if self.item_inspector: self.item_inspector.updateMetaData()
	
	def setRenderString(self, string, fontSize):
		self.renderString = string
		self.fontSize = fontSize
		if self.item_inspector: self.item_inspector.updateMetaData()
		
	def getEvalString(self):
		return "EM3DText('%s', %s, fontMode=%d, depth=%d)"%(self.renderString, self.fontSize, self.fontMode, self.fontDepth)
	
	def getFontMode(self):
		return self.fontMode
		
	def getFontDepth(self):
		return self.fontDepth
	
	def getFontSize(self):
		return self.fontSize
	
	def getRenderString(self):
		return self.renderString
		
		
	def getItemInspector(self):
		"""
		Return a Qt widget that controls the scene item
		"""
		if not self.item_inspector: self.item_inspector = EMInspectorControl3DText("3DText", self)
		return self.item_inspector
		
	def getItemDictionary(self):
		"""
		Return a dictionary of item parameters (used for restoring sessions
		"""
		dictionary = super(EM3DText, self).getItemDictionary()
		dictionary.update({"TEXT":[self.getFontMode(), self.getFontDepth(), self.getFontSize(), self.getRenderString()]})
		return dictionary
	
	def setUsingDictionary(self, dictionary):
		"""
		Set item attributes using a dictionary, used in session restoration
		"""
		super(EM3DText, self).setUsingDictionary(dictionary)
		self.setFontMode(dictionary["TEXT"][0])
		self.setFontDepth(dictionary["TEXT"][1])
		self.setRenderString(dictionary["TEXT"][3],dictionary["TEXT"][2])	
	
	def renderShape(self):
		# Material properties of the 3D text
		glDisable(GL_COLOR_MATERIAL)
		glMaterialfv(GL_FRONT, GL_DIFFUSE, self.diffuse)
		glMaterialfv(GL_FRONT, GL_SPECULAR, self.specular)
		glMaterialf(GL_FRONT, GL_SHININESS, self.shininess)
		glMaterialfv(GL_FRONT, GL_AMBIENT, self.ambient)
		
		glEnable(GL_NORMALIZE)
		#HERE
		glPushMatrix()
		glNormal(0,0,1)
		glEnable(GL_TEXTURE_2D)
		
		self.font_renderer.set_font_mode(self.fontMode)
		self.font_renderer.set_depth(self.fontDepth)
		self.font_renderer.set_face_size(int(self.fontSize))
		
		#make 3D text rotate at the center
		tvar = self.font_renderer.bounding_box(self.renderString)
		glTranslate((tvar[0]-tvar[3])/2,(tvar[1]-tvar[4])/2,-(tvar[2]-tvar[5])/2)
		
		self.font_renderer.render_string(self.renderString)
		glPopMatrix()

class EMInspectorControlShape(EMItem3DInspector):
	"""
	Class to make EMItem GUI SHAPE Inspector
	"""
	def __init__(self, name, item3d):
		EMItem3DInspector.__init__(self, name, item3d)
		
	def updateItemControls(self):
		""" Updates this item inspector. Function is called by the item it observes"""
		super(EMInspectorControlShape, self).updateItemControls()
		self.ambcolorbox.setColor(QtGui.QColor(255*self.item3d().ambient[0],255*self.item3d().ambient[1],255*self.item3d().ambient[2]))
		self.diffusecolorbox.setColor(QtGui.QColor(255*self.item3d().diffuse[0],255*self.item3d().diffuse[1],255*self.item3d().diffuse[2]))
		self.specularcolorbox.setColor(QtGui.QColor(255*self.item3d().specular[0],255*self.item3d().specular[1],255*self.item3d().specular[2]))
	
	def addTabs(self):
		""" Add a tab for each 'column' """
		tabwidget = QtGui.QWidget()
		gridbox = QtGui.QGridLayout()
		
		EMInspectorControlShape.addControls(self, gridbox)
		
		tabwidget.setLayout(gridbox)
		self.addTab(tabwidget, "basic")
			
	def addControls(self, gridbox):
		""" Construct all the widgets in this Item Inspector """
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
		if self.item3d():
			self.item3d().setAmbientColor((float(rgb[0])/255.0),(float(rgb[1])/255.0),(float(rgb[2])/255.0))
			self.inspector().updateSceneGraph()
		
	def _on_diffuse_color(self, color):
		rgb = color.getRgb()
		if self.item3d():
			self.item3d().setDiffuseColor((float(rgb[0])/255.0),(float(rgb[1])/255.0),(float(rgb[2])/255.0))
			self.inspector().updateSceneGraph()
		
	def _on_specular_color(self, color):
		rgb = color.getRgb()
		if self.item3d():	
			self.item3d().setSpecularColor((float(rgb[0])/255.0),(float(rgb[1])/255.0),(float(rgb[2])/255.0))
			self.inspector().updateSceneGraph()
		
	def _on_shininess(self, shininess):
		if self.item3d():
			self.item3d().setShininess(shininess)
			self.inspector().updateSceneGraph()

class EMInspectorControlScatterPlot(EMInspectorControlShape):
	"""
	Class to make scatter plot 3D inspector 
	"""
	def __init__(self, name, item3d):
		EMInspectorControlShape.__init__(self, name, item3d)
		
	def updateItemControls(self):
		""" Updates this item inspector. Function is called by the item it observes"""
		super(EMInspectorControlScatterPlot, self).updateItemControls()
		self.pointsize.setValue(self.item3d().getPointSize())
		
	def addTabs(self):
		""" Add a tab for each 'column' """
		super(EMInspectorControlScatterPlot, self).addTabs()
		tabwidget = QtGui.QWidget()
		gridbox = QtGui.QGridLayout()
		
		EMInspectorControlScatterPlot.addControls(self, gridbox)
		
		tabwidget.setLayout(gridbox)
		self.addTab(tabwidget, "scatterplot")
		
	def addControls(self, gridbox):
		""" Construct all the widgets in this Item Inspector """
		
		scatterframe = QtGui.QFrame()
		scatterframe.setFrameShape(QtGui.QFrame.StyledPanel)
		lfont = QtGui.QFont()
		lfont.setBold(True)
		scattergridbox = QtGui.QGridLayout()
		scattergridbox.setAlignment(QtCore.Qt.AlignTop)
		
		# Add widgets to frame
		pslabel = QtGui.QLabel("Point Size")
		pslabel.setFont(lfont)
		pslabel.setAlignment(QtCore.Qt.AlignCenter)
		scattergridbox.addWidget(pslabel, 0, 0, 1, 1)
		
		self.pointsize = EMSpinWidget(int(self.item3d().getPointSize()), 1.0, rounding=0)
		self.pointsize.setMinimumWidth(120)
		scattergridbox.addWidget(self.pointsize, 1, 0, 1, 1)
		
		scatterframe.setLayout(scattergridbox)
		gridbox.addWidget(scatterframe, 3, 0)
		
		# set to default, but run only as a base class
		if type(self) == EMInspectorControlScatterPlot: 
			self.updateItemControls()
		
		QtCore.QObject.connect(self.pointsize,QtCore.SIGNAL("valueChanged(int)"),self.onPointSizeChanged)
	
	def onPointSizeChanged(self):
		self.item3d().setPointSize(self.pointsize.getValue())
		if self.inspector: self.inspector().updateSceneGraph()
		
class EMInspectorControl3DText(EMInspectorControlShape):
	"""
	Class to make EMItem GUI SHAPE 3DText Inspector
	"""
	def __init__(self, name, item3d):
		EMInspectorControlShape.__init__(self, name, item3d)

	def updateItemControls(self):
		""" Updates this item inspector. Function is called by the item it observes"""
		super(EMInspectorControl3DText, self).updateItemControls()
		self.fontDepth.setValue(int(self.item3d().getFontDepth()))
		self.fontSize.setValue(int(self.item3d().getFontSize()))
	
	def updateMetaData(self):
		""" Updates the items metadata, such as line length, width. Function is called by the item it observes when the items meta data changes """
		super(EMInspectorControlShape, self).updateMetaData()
		if self.item3d().getFontMode() == FTGLFontMode.EXTRUDE:
			self.textModeBox.setCurrentIndex(0)
		if self.item3d().getFontMode() == FTGLFontMode.TEXTURE:
			self.textModeBox.setCurrentIndex(1)
		if self.item3d().getFontMode() == FTGLFontMode.POLYGON:
			self.textModeBox.setCurrentIndex(2)
		if self.item3d().getFontMode() == FTGLFontMode.OUTLINE:
			self.textModeBox.setCurrentIndex(3)
	
	def addTabs(self):
		""" Add a tab for each 'column' """
		super(EMInspectorControl3DText, self).addTabs()
		tabwidget = QtGui.QWidget()
		gridbox = QtGui.QGridLayout()
		
		EMInspectorControl3DText.addControls(self, gridbox)
		
		tabwidget.setLayout(gridbox)
		self.addTab(tabwidget, "text")
		
	def addControls(self, gridbox):
		""" Construct all the widgets in this Item Inspector """
	
		textframe = QtGui.QFrame()
		textframe.setFrameShape(QtGui.QFrame.StyledPanel)
		lfont = QtGui.QFont()
		lfont.setBold(True)
		textgridbox = QtGui.QGridLayout()
		
		# Add widgets to textframe
		textlabel = QtGui.QLabel("3D Font Mode")
		textlabel.setFont(lfont)
		textlabel.setAlignment(QtCore.Qt.AlignCenter)
		textgridbox.addWidget(textlabel, 0, 0, 1, 1)
		
		self.textModeBox = QtGui.QComboBox()
		self.textModeBox.addItems(["EXTRUDE", "TEXTURE", "POLYGON", "OUTLINE"])
		textgridbox.addWidget(self.textModeBox, 0, 1, 1, 1)
			
		textlabel2 = QtGui.QLabel("3D Font Depth")
		textlabel2.setFont(lfont)
		textlabel2.setAlignment(QtCore.Qt.AlignCenter)
		textgridbox.addWidget(textlabel2, 1, 0, 1, 1)
		
		self.fontDepth = EMSpinWidget(int(self.item3d().getFontDepth()), 1.0, rounding=0)
		self.fontDepth.setMinimumWidth(120)
		textgridbox.addWidget(self.fontDepth, 1, 1, 1, 1)
		
		textlabel2 = QtGui.QLabel("3D Font Size")
		textlabel2.setFont(lfont)
		textlabel2.setAlignment(QtCore.Qt.AlignCenter)
		textgridbox.addWidget(textlabel2, 2, 0, 1, 1)
		
		self.fontSize = EMSpinWidget(int(self.item3d().getFontSize()), 1.0, rounding=0)
		self.fontSize.setMinimumWidth(120)
		textgridbox.addWidget(self.fontSize, 2, 1, 1, 1)
		
		textframe.setLayout(textgridbox)	
		gridbox.addWidget(textframe, 2, 0)
		
		# Add text
		text3dframe = QtGui.QFrame()
		text3dframe.setFrameShape(QtGui.QFrame.StyledPanel)
		text3dgridbox = QtGui.QGridLayout()
		
		textlabel3 = QtGui.QLabel("3D Text")
		textlabel3.setFont(lfont)
		text3dgridbox.addWidget(textlabel3, 3, 0, 2, 1)
		
		self.text3d = QtGui.QLineEdit(self.item3d().getRenderString())
		text3dgridbox.addWidget(self.text3d, 3, 1, 2, 1)
		
		text3dframe.setLayout(text3dgridbox)
		gridbox.addWidget(text3dframe, 3, 0)
		
		# set to default, but run only as a base class
		if type(self) == EMInspectorControl3DText: 
			self.updateItemControls()
			self.updateMetaData()
			
		self.textModeBox.currentIndexChanged.connect(self.on3DTextModeChanged)
		QtCore.QObject.connect(self.fontDepth,QtCore.SIGNAL("valueChanged(int)"),self.on3DTextDepthChanged)
		QtCore.QObject.connect(self.fontSize,QtCore.SIGNAL("valueChanged(int)"),self.on3DTextFontChanged)
		QtCore.QObject.connect(self.text3d,QtCore.SIGNAL("textChanged(const QString&)"),self.on3DTextChanged)
		
	def on3DTextModeChanged(self):
		textMode = str(self.textModeBox.currentText())
		if textMode == "EXTRUDE":
			self.item3d().setFontMode(FTGLFontMode.EXTRUDE)
		elif textMode == "TEXTURE":
			self.item3d().setFontMode(FTGLFontMode.TEXTURE)
		elif textMode == "POLYGON":
			self.item3d().setFontMode(FTGLFontMode.POLYGON)
		elif textMode == "OUTLINE":
			self.item3d().setFontMode(FTGLFontMode.OUTLINE)	
		if self.inspector: self.inspector().updateSceneGraph()
		
	def on3DTextDepthChanged(self):
		self.item3d().setFontDepth(int(self.fontDepth.getValue()))
		if self.inspector: self.inspector().updateSceneGraph()
	
	def on3DTextFontChanged(self):
		self.item3d().setRenderString(self.item3d().getRenderString(), int(self.fontSize.getValue()))
		if self.inspector: self.inspector().updateSceneGraph()
		
	def on3DTextChanged(self, string):
		self.item3d().setRenderString(str(string), self.item3d().getFontSize())
		if self.inspector: self.inspector().updateSceneGraph()
		
class EMInspectorControlLine(EMInspectorControlShape):
	"""
	Class to make EMItem GUI SHAPE Line Inspector
	"""
	def __init__(self, name, item3d):
		EMInspectorControlShape.__init__(self, name, item3d)
		
	def updateItemControls(self):
		"""Updates this item inspector. Function is called by the item it observes"""
		super(EMInspectorControlLine, self).updateItemControls()
	
	def updateMetaData(self):
		"""Updates the items metadata, such as line length, width. Function is called by the item it observes when the items meta data changes"""
		super(EMInspectorControlLine, self).updateMetaData()
		self.leftArrowSize.setValue(self.item3d().leftArrowSize, quiet=1)
		self.leftArrowLength.setValue(self.item3d().leftArrowLength, quiet=1)
		self.rightArrowSize.setValue(self.item3d().rightArrowSize, quiet=1)
		self.rightArrowLength.setValue(self.item3d().rightArrowLength, quiet=1)
		self.slice.setValue(self.item3d().slices, quiet=1)
		self.stack.setValue(self.item3d().stacks, quiet=1)
		self.linelength.setValue(int(self.item3d().length), quiet=1)
	
	def addTabs(self):
		""" Add a tab for each 'column' """
		super(EMInspectorControlLine, self).addTabs()
		tabwidget = QtGui.QWidget()
		gridbox = QtGui.QGridLayout()
		
		EMInspectorControlLine.addControls(self, gridbox)
		
		tabwidget.setLayout(gridbox)
		self.addTab(tabwidget, "line")
		
	def addControls(self, gridbox):
		""" Construct all the widgets in this Item Inspector """
		#frame to control properties of left/right arrows
		lineframe = QtGui.QFrame()
		lineframe.setFrameShape(QtGui.QFrame.StyledPanel)
		lfont = QtGui.QFont()
		lfont.setBold(True)
		linegridbox = QtGui.QGridLayout()
		
		leftlabel = QtGui.QLabel("Left arrow")
		leftlabel.setFont(lfont)
		leftlabel.setAlignment(QtCore.Qt.AlignCenter)
		linegridbox.addWidget(leftlabel, 0, 1, 1, 1)
		
		sidelabel1 = QtGui.QLabel("Size")
		sidelabel1.setFont(lfont)
		sidelabel1.setAlignment(QtCore.Qt.AlignVCenter)
		linegridbox.addWidget(sidelabel1, 2, 0, 1, 1)
		
		sidelabel2 = QtGui.QLabel("Length")
		sidelabel2.setFont(lfont)
		sidelabel2.setAlignment(QtCore.Qt.AlignVCenter)
		linegridbox.addWidget(sidelabel2, 3, 0, 1, 1)
		
		self.leftShowArrow = QtGui.QCheckBox("Show")
		self.leftShowArrow.setChecked(self.item3d().showLeftArrow)
		linegridbox.addWidget(self.leftShowArrow, 1, 1, 1, 1)
		
		self.leftArrowSize = EMSpinWidget(int(self.item3d().leftArrowSize), 1.0, rounding=0)
		self.leftArrowSize.setMinimumWidth(120)
		linegridbox.addWidget(self.leftArrowSize, 2, 1, 1, 1)
		
		self.leftArrowLength = EMSpinWidget(int(self.item3d().leftArrowLength), 1.0, rounding=0)
		self.leftArrowLength.setMinimumWidth(120)
		linegridbox.addWidget(self.leftArrowLength, 3, 1, 1, 1)
		
		rightlabel = QtGui.QLabel("Right arrow")
		rightlabel.setFont(lfont)
		rightlabel.setAlignment(QtCore.Qt.AlignCenter)
		linegridbox.addWidget(rightlabel, 0, 2, 1, 1)
		
		self.rightShowArrow = QtGui.QCheckBox("Show")
		self.rightShowArrow.setChecked(self.item3d().showRightArrow)
		linegridbox.addWidget(self.rightShowArrow, 1, 2, 1, 1)
		
		self.rightArrowSize = EMSpinWidget(int(self.item3d().rightArrowSize), 1.0, rounding=0)
		self.rightArrowSize.setMinimumWidth(120)
		linegridbox.addWidget(self.rightArrowSize, 2, 2, 1, 1)
		
		self.rightArrowLength = EMSpinWidget(int(self.item3d().rightArrowLength), 1.0, rounding=0)
		self.rightArrowLength.setMinimumWidth(120)
		linegridbox.addWidget(self.rightArrowLength, 3, 2, 1, 1)
		
		linelengthlabel = QtGui.QLabel("Line Length")
		linelengthlabel.setFont(lfont)
		linelengthlabel.setAlignment(QtCore.Qt.AlignCenter)
		linegridbox.addWidget(linelengthlabel, 4, 0, 2, 2)
		
		self.linelength = EMSpinWidget(int(self.item3d().length), 1.0, rounding=0)
		linegridbox.addWidget(self.linelength, 4, 2, 2, 2)
		
		linewidthlabel = QtGui.QLabel("Line Width")
		linewidthlabel.setFont(lfont)
		linewidthlabel.setAlignment(QtCore.Qt.AlignCenter)
		linegridbox.addWidget(linewidthlabel, 5, 0, 1, 2)
		
		self.linewidth = EMSpinWidget(int(self.item3d().width), 1.0, rounding=0)
		linegridbox.addWidget(self.linewidth, 5, 2, 1, 2)
		
		lineframe.setLayout(linegridbox)	
		gridbox.addWidget(lineframe, 2, 0)
		
		#frame to control slice/stack of the line
		lineframe2 = QtGui.QFrame()
		lineframe2.setFrameShape(QtGui.QFrame.StyledPanel)
		linehbox = QtGui.QVBoxLayout()
				
		self.slice = ValSlider(lineframe2, (1, 100), "Slice", rounding=0)
		self.slice.setValue(self.item3d().slices)
		
		self.stack = ValSlider(lineframe2, (1, 100), "Stack", rounding=0)
		self.slice.setValue(self.item3d().stacks)
		
		linehbox.addWidget(self.slice)
		linehbox.addWidget(self.stack)
		
		lineframe2.setLayout(linehbox)
		gridbox.addWidget(lineframe2, 3, 0)
		
		# set to default, but run only as a base class
		if type(self) == EMInspectorControl3DText: 
			self.updateItemControls()
			self.updateMetaData()
		
		QtCore.QObject.connect(self.leftShowArrow, QtCore.SIGNAL("stateChanged(int)"), self.redraw)
		QtCore.QObject.connect(self.rightShowArrow, QtCore.SIGNAL("stateChanged(int)"), self.redraw)
		QtCore.QObject.connect(self.leftArrowSize,QtCore.SIGNAL("valueChanged(int)"),self.redraw)
		QtCore.QObject.connect(self.leftArrowLength,QtCore.SIGNAL("valueChanged(int)"),self.redraw)
		QtCore.QObject.connect(self.rightArrowSize,QtCore.SIGNAL("valueChanged(int)"),self.redraw)
		QtCore.QObject.connect(self.rightArrowLength,QtCore.SIGNAL("valueChanged(int)"),self.redraw)
		QtCore.QObject.connect(self.linelength,QtCore.SIGNAL("valueChanged(int)"),self.redraw)
		QtCore.QObject.connect(self.linewidth,QtCore.SIGNAL("valueChanged(int)"),self.redraw)
		
		QtCore.QObject.connect(self.slice,QtCore.SIGNAL("valueChanged"),self.redraw)
		QtCore.QObject.connect(self.stack,QtCore.SIGNAL("valueChanged"),self.redraw)
	
	def redraw(self):
		self.item3d().setShowLeftArrow(self.leftShowArrow.isChecked())
		self.item3d().setShowRightArrow(self.rightShowArrow.isChecked())
		self.item3d().leftArrowSize = self.leftArrowSize.getValue()
		self.item3d().leftArrowLength = self.leftArrowLength.getValue()
		self.item3d().rightArrowSize = self.rightArrowSize.getValue()
		self.item3d().rightArrowLength = self.rightArrowLength.getValue()
		self.item3d().setLength(self.linelength.getValue())
		self.item3d().setWidth(self.linewidth.getValue())
		
		self.item3d().setSlices(self.slice.getValue())
		self.item3d().setStacks(self.stack.getValue())
		
		if self.inspector:
			self.inspector().updateSceneGraph()