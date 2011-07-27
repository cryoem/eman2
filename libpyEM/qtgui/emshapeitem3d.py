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
from emscene3d import EMScene3D, EMInspectorControlShape
from EMAN2 import Transform

class EMCube(EMItem3D):
	name = "Cube"
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
		if self.is_selected:
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
	pass
