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
from OpenGL import GLUT

from emglobjects import init_glut
from emitem3d import EMItem3D
from emscene3d import EMScene3D, EMInspector3D, EMInspectorControlShape
from EMAN2 import Transform

class EMCube(EMItem3D):
	def __init__(self, size):
		EMItem3D.__init__(self, parent=None, children=set(), transform=Transform())
		# size
		self.size = size  
		self.boundingboxsize = size
		# color
		self.diffuse = [0.5,0.5,0.5,1.0]
		self.specular = [1.0,1.0,1.0,1.0]
		self.ambient = [1.0, 1.0, 1.0, 1.0]
		self.shininess = 25.0
		
	def setAmbientColor(self, red, green, blue, alpha=1.0):
		self.ambient = [red, green, blue, alpha]

	def setDiffuseColor(self, red, green, blue, alpha=1.0):
		self.diffuse = [red, green, blue, alpha]
		
	def setSpecularColor(self, red, green, blue, alpha=1.0):
		self.specular = [red, green, blue, alpha]
	
	def setShininess(self, shininess):
		self.shininess = shininess
		
	def getSceneGui(self):
		"""
		Return a Qt widget that controls the scene item
		"""
		if not self.widget: self.widget = EMInspectorControlShape("CUBE", self)
		return self.widget
		
	def renderNode(self):        
		# Material properties of the box
		glMaterialfv(GL_FRONT, GL_DIFFUSE, self.diffuse)
		glMaterialfv(GL_FRONT, GL_SPECULAR, self.specular)
		glMaterialf(GL_FRONT, GL_SHININESS, 25.0)
		glMaterialfv(GL_FRONT, GL_AMBIENT, self.ambient)
		
		init_glut()
		GLUT.glutSolidCube(self.size)


class EMSphere(EMItem3D):
	def __init__(self, radius):
		EMItem3D.__init__(self, parent=None, children=set(), transform=Transform())
		# size
		self.radius = radius
		self.slices = int(radius)
		self.stacks = int(radius)
		self.boundingboxsize = radius
		
		# color
		self.diffuse = [0.5,0.5,0.5,1.0]
		self.specular = [1.0,1.0,1.0,1.0]
		self.ambient = [1.0, 1.0, 1.0, 1.0]
		self.shininess = 25.0
		
	def setAmbientColor(self, red, green, blue, alpha=1.0):
		self.ambient = [red, green, blue, alpha]

	def setDiffuseColor(self, red, green, blue, alpha=1.0):
		self.diffuse = [red, green, blue, alpha]
		
	def setSpecularColor(self, red, green, blue, alpha=1.0):
		self.specular = [red, green, blue, alpha]
	
	def setShininess(self, shininess):
		self.shininess = shininess
		
	def getSceneGui(self):
		"""
		Return a Qt widget that controls the scene item
		"""
		if not self.widget: self.widget = EMInspectorControlShape("SPHERE", self)
		return self.widget
				
	def renderNode(self):
		# Material properties of the sphere
		glMaterialfv(GL_FRONT, GL_DIFFUSE, self.diffuse)
		glMaterialfv(GL_FRONT, GL_SPECULAR, self.specular)
		glMaterialf(GL_FRONT, GL_SHININESS, 25.0)
		glMaterialfv(GL_FRONT, GL_AMBIENT, self.ambient)
		
		init_glut()
		GLUT.glutSolidSphere(self.radius, self.slices, self.stacks)

class EMCylinder(EMItem3D):
	def __init__(self, radius, height):
		EMItem3D.__init__(self, parent=None, children=set(), transform=Transform())
		#size
		self.radius = radius
		self.height = height
		self.slices = int(radius)
		self.stacks = int(radius)
		self.boundingboxsize = max(self.radius, self.height)

		# color
		self.diffuse = [0.5,0.5,0.5,1.0]
		self.specular = [1.0,1.0,1.0,1.0]
		self.ambient = [1.0, 1.0, 1.0, 1.0]		
		self.shininess = 25.0
		
	def setAmbientColor(self, red, green, blue, alpha=1.0):
		self.ambient = [red, green, blue, alpha]

	def setDiffuseColor(self, red, green, blue, alpha=1.0):
		self.diffuse = [red, green, blue, alpha]
		
	def setSpecularColor(self, red, green, blue, alpha=1.0):
		self.specular = [red, green, blue, alpha]
	
	def setShininess(self, shininess):
		self.shininess = shininess
		
	def getSceneGui(self):
		"""
		Return a Qt widget that controls the scene item
		"""
		if not self.widget: self.widget = EMInspectorControlShape("CYLINDER", self)
		return self.widget
		
	def renderNode(self):
		# Material properties of the cylinder
		glMaterialfv(GL_FRONT, GL_DIFFUSE, self.diffuse)
		glMaterialfv(GL_FRONT, GL_SPECULAR, self.specular)
		glMaterialf(GL_FRONT, GL_SHININESS, 25.0)
		glMaterialfv(GL_FRONT, GL_AMBIENT, self.ambient)
		
		init_glut()
		GLUT.glutSolidCylinder(self.radius, self.height, self.slices, self.stacks)

class EMLine(EMItem3D):
	pass
