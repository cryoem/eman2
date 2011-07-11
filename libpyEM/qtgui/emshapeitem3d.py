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
		
		quadratic = gluNewQuadric()
		gluQuadricNormals(quadratic, GLU_SMOOTH)    # Create Smooth Normals (NEW) 
		gluQuadricTexture(quadratic, GL_TRUE)      # Create Texture Coords (NEW) 

		gluSphere(quadratic,self.radius,self.slices,self.stacks)
		
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
		
		quadratic = gluNewQuadric()
		gluQuadricNormals(quadratic, GLU_SMOOTH)    # Create Smooth Normals (NEW) 
		gluQuadricTexture(quadratic, GL_TRUE)      # Create Texture Coords (NEW) 
		gluCylinder(quadratic,self.radius,self.radius,self.height,self.slices,self.stacks)

class EMLine(EMItem3D):
	pass
