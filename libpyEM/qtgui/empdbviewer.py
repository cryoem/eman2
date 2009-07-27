#!/usr/bin/env python

# Author: Muthu Alagappan, 07/22/09
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

from emglobjects import EMImage3DGUIModule,Camera2,EMViewportDepthTools,get_default_gl_colors,get_3d_font_renderer
from emlights import EMLightsDrawer

from EMAN2 import get_image_directory
from EMAN2 import *
from libpyGLUtils2 import *
from math import *

import weakref

class EM3DModule(EMLightsDrawer,EMImage3DGUIModule):
	
	def get_qt_widget(self):
		if self.qt_context_parent == None:	
			from emimageutil import EMParentWin
			from emimage3d import EMImage3DWidget
			self.under_qt_control = True
			self.gl_context_parent = EMImage3DWidget(self)
			self.qt_context_parent = EMParentWin(self.gl_context_parent)
			self.gl_widget = self.gl_context_parent			
			self.qt_context_parent.setWindowIcon(QtGui.QIcon(get_image_directory() +"single_image_3d.png"))
		
		return self.qt_context_parent
	
	def get_gl_widget(self,qt_context_parent,gl_context_parent):
		self.under_qt_control = False
		ret = EMImage3DGUIModule.get_gl_widget(self,qt_context_parent,gl_context_parent)
		self.gl_widget.setWindowTitle(remove_directories_from_name(self.file_name))
		self.__set_module_contexts()
		return ret
	
	def get_desktop_hint(self):
		return "image"

	def __init__(self,application=None):
		EMImage3DGUIModule.__init__(self,application,ensure_gl_context=True)
		EMLightsDrawer.__init__(self)
		self.cam = Camera2(self)
		self.vdtools = EMViewportDepthTools(self)
		self.perspective = False
		self.colors = get_default_gl_colors()
		self.perspective = True
		
		# basic shapes will be stored in these lists
		self.gq = None # will be a glu quadric
		self.cylinderdl = 0 # will be a cylinder with no caps
		self.diskdl = 0 # this will be a flat disk
		self.spheredl = 0 # this will be a low resolution sphere
		self.highresspheredl = 0 # high resolution sphere
		self.cappedcylinderdl = 0 # a capped cylinder
		self.first_render_flag = True # this is used to catch the first call to the render function - so you can do an GL context sensitive initialization when you know there is a valid context
	
		self.inspector = None # will be the inspector, i.e. an instance of an EM3DInspector
		self.radius = 100
		
		self.font_renderer = None # will be a 3D fonts renderer
		
		
	def __del__(self):
		if self.under_qt_control and not self.dont_delete_parent:
			self.qt_context_parent.deleteLater()
		self.core_object.deleteLater()

	def width(self):
		try: return self.gl_widget.width()
		except: return 0
		
	def height(self):
		try: return self.gl_widget.height()
		except: return 0
	
	def initializeGL(self):
		# put your own initialization things in here
		glEnable(GL_LIGHTING)
		glEnable(GL_NORMALIZE)
		
	def load_gl_color(self,name):
		color = self.colors[name]
		glColor(color["ambient"])
		glMaterial(GL_FRONT,GL_AMBIENT,color["ambient"])
		glMaterial(GL_FRONT,GL_DIFFUSE,color["diffuse"])
		glMaterial(GL_FRONT,GL_SPECULAR,color["specular"])
		glMaterial(GL_FRONT,GL_EMISSION,color["emission"])
		glMaterial(GL_FRONT,GL_SHININESS,color["shininess"])
	
	
	def render(self):
		if self.first_render_flag:
			self.initializeGL()
			self.init_basic_shapes() # only does something the first time you call it
			self.init_font_renderer()
			if not self.perspective:self.gl_context_parent.load_orthographic()
			else: self.gl_context_parent.load_perspective()
					
		#self.vdtools.set_update_P_inv()
		glPushMatrix()
		self.cam.position(True)
		# the ones are dummy variables atm... they don't do anything
		self.vdtools.update(1,1)
		glPopMatrix()
		
		glPushMatrix()
		self.cam.position()
	
		self.draw_objects()
		
		glPopMatrix()
		
		glPushMatrix()
		self.cam.translate_only()
		EMLightsDrawer.draw(self)
		glPopMatrix()
		
	
	def draw_objects(self):
#		glPushMatrix()
#		glScale(50,50,50)
#		self.load_gl_color("red")
#		glCallList(self.highresspheredl)
#		glPopMatrix()
		
		glPushMatrix()
		glTranslate(100,0,0)
		glScale(50,50,50)
		self.load_gl_color("blue")
		glCallList(self.spheredl)
		glPopMatrix()
#		
#		glPushMatrix()
#		glTranslate(0,100,0)
#		glScale(50,10,10)
#		glTranslate(-0.5,0,0)
#		glRotate(90,0,1,0)
#		self.load_gl_color("emerald")
#		glCallList(self.cylinderdl)
#		glPopMatrix()
#		
#		glPushMatrix()
#		glTranslate(-100,25,0)
#		glScale(10,50,10)
#		glTranslate(-0.5,0,0)
#		glRotate(90,1,0,0)
#		self.load_gl_color("gold")
#		glCallList(self.cappedcylinderdl)
#		glPopMatrix()
#		
#		glPushMatrix()
#		glTranslate(0,-100,0)
#		glScale(15,15,15)
#		self.load_gl_color("copper")
#		glCallList(self.diskdl)
#		glPopMatrix()
#		
#		glPushMatrix()
#		glTranslate(0,-100,-10)
#		glScale(15,15,15)
#		glRotate(180,0,1,0)
#		self.load_gl_color("silver")
#		glCallList(self.diskdl)
#		glPopMatrix()
		
		glPushMatrix()
		glTranslate(0,0,50)
		s = "bdb:EMAN2"
		
		bbox = self.font_renderer.bounding_box(s)
		glTranslate(-(bbox[3]-bbox[0])/2, -(bbox[4]-bbox[1])/2,-(bbox[5]-bbox[02])/2)
		self.font_renderer.render_string(s)
		glPopMatrix()
		
		
#		glPushMatrix()
#		glTranslate(0,-50,-50)
#		s = "=^_^="
#		bbox = self.font_renderer.bounding_box(s)
#		glTranslate(-(bbox[3]-bbox[0])/2, -(bbox[4]-bbox[1])/2,-(bbox[5]-bbox[02])/2)
#		self.font_renderer.render_string(s)
#		glPopMatrix()
	
	
	def init_font_renderer(self):
		if self.font_renderer == None:
			self.font_renderer = get_3d_font_renderer()
			self.font_renderer.set_face_size(20)
			self.font_renderer.set_depth(12)
			self.font_renderer.set_font_mode(FTGLFontMode.EXTRUDE)
		
		
	def init_basic_shapes(self):
		self.gl_context_parent.makeCurrent()
		if self.gq == None:
			
			self.gq=gluNewQuadric() # a quadric for general use
			gluQuadricDrawStyle(self.gq,GLU_FILL)
			gluQuadricNormals(self.gq,GLU_SMOOTH)
			gluQuadricOrientation(self.gq,GLU_OUTSIDE)
			gluQuadricTexture(self.gq,GL_FALSE)
		
		if ( self.cylinderdl == 0 ):
			self.cylinderdl=glGenLists(1)
				
			glNewList(self.cylinderdl,GL_COMPILE)
			glPushMatrix()
			gluCylinder(self.gq,1.0,1.0,1.0,12,2)
			glPopMatrix()
				
			glEndList()
		
		if self.diskdl == 0:
			self.diskdl=glGenLists(1)
				
			glNewList(self.diskdl,GL_COMPILE)
			gluDisk(self.gq,0,1,12,2)
			glEndList()
		
		if self.spheredl == 0:
			self.spheredl=glGenLists(1)
				
			glNewList(self.spheredl,GL_COMPILE)
			gluSphere(self.gq,.5,4,2)
			glEndList()

		
		if self.highresspheredl == 0:
			self.highresspheredl=glGenLists(1)
				
			glNewList(self.highresspheredl,GL_COMPILE)
			gluSphere(self.gq,.5,16,16)
			glEndList()
			
		if ( self.cappedcylinderdl == 0 ):
			self.cappedcylinderdl=glGenLists(1)
			glNewList(self.cappedcylinderdl,GL_COMPILE)
			glCallList(self.cylinderdl)
			glPushMatrix()
			glTranslate(0,0,1)
			glCallList(self.diskdl)
			glPopMatrix()
			glPushMatrix()
			glRotate(180,0,1,0)
			glCallList(self.diskdl)
			glPopMatrix()
			glEndList()
			
	def eye_coords_dif(self,x1,y1,x2,y2,mdepth=True):
		return self.vdtools.eye_coords_dif(x1,y1,x2,y2,mdepth)
	
	def resizeEvent(self, width, height):
		for i in self.viewables:
			i.resizeEvent()
	
	def get_inspector(self):
		if self.inspector == None:
			self.inspector = EM3DInspector(self)
		return self.inspector

	def set_cam_z(self,z):
		self.cam.set_cam_z( z )
		self.updateGL()
		
	def set_cam_y(self,y):
		self.cam.set_cam_y( y )
		self.updateGL()
		
	def set_cam_x(self,x):
		self.cam.set_cam_x( x )
		self.updateGL()
	
	def set_scale(self,val):
		self.cam.scale = val
		self.updateGL()	

	def resizeEvent(self,width=0,height=0):
		self.vdtools.set_update_P_inv()
	
	def load_rotation(self,t3d):
		self.cam.load_rotation(t3d)
		self.updateGL()

	def get_start_z(self):
		return self.gl_context_parent.get_start_z()
	
	def get_near_plane_dims(self):
		return self.gl_context_parent.get_near_plane_dims()
	
	def set_perspective(self,bool):
		self.perspective = bool
		
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		if not self.perspective:self.gl_context_parent.load_orthographic()
		else: self.gl_context_parent.load_perspective()
		glMatrixMode(GL_MODELVIEW)
		
		self.updateGL()


class EM3DInspector(QtGui.QWidget):
	def __init__(self,target):
		QtGui.QWidget.__init__(self,None)
		self.target = weakref.ref(target) # prevent a strong cycle - this target object should be an EM3DModule, but that could change depending on who builds on this object
		
		self.vbl = QtGui.QVBoxLayout(self) # this is the main vbl
		
		self.advanced_tab = None
		
		self.tabwidget = QtGui.QTabWidget()
		
		#self.vbl.addLayout(self.hbl_check)
		self.vbl.addWidget(self.tabwidget)
		
		self.insert_advanced_tab()
		
	def insert_advanced_tab(self):
		if self.advanced_tab == None:
			from emimage3d import EM3DAdvancedInspector
			self.advanced_tab = EM3DAdvancedInspector(self.target(), self)
			
		self.advanced_tab.update_rotations(self.target().get_current_transform())
		self.advanced_tab.set_scale(self.target().cam.scale)
		self.tabwidget.addTab(self.advanced_tab,"Advanced")
		self.settingsrow = self.tabwidget.count()-1
		self.tabwidget.setCurrentIndex(self.settingsrow)

	
	def update_rotations(self,t3d):
		self.advanced_tab.update_rotations(t3d)
	
	def set_scale(self,val):
		self.advanced_tab.set_scale(val)
	
	def set_xy_trans(self, x, y):
		self.advanced_tab.set_xy_trans(x,y)
	
	def set_xyz_trans(self,x,y,z):
		self.advanced_tab.set_xyz_trans(x,y,z)
		
	def set_directional_light_dir(self,d):
		self.advanced_tab.set_directional_light_dir(d)
	
	def set_positional_light_pos(self,d):
		self.advanced_tab.set_positional_light_pos(d)
		
	def set_positional_light_dir(self,d):
		self.advanced_tab.set_positional_light_dir(d)
	


















class EMPDBViewer(EM3DModule):
	def __init__(self, application=None):
		EM3DModule.__init__(self,application)
		#self.fName = raw_input ("Enter the file name of a pdb file: ")
		self.fName = ""
		self.text = self.fName
		self.dl = None
	
	def current_text(self): return self.text
	
	def set_current_text(self,text):
		self.text = text
	
	def get_inspector(self):
		if self.inspector == None:
			self.inspector = EMPDBInspector(self)
		return self.inspector
		
	def draw_objects(self):

		if (self.text == ""): return
		
		if (self.text != self.fName): 
			self.dl=None
			self.fName = self.text

		if (self.dl == None):
			self.dl=glGenLists(1)
			glNewList(self.dl,GL_COMPILE)
			self.buildResList()
			
			for res in self.allResidues:
				for i in range (0, len(res[0])):
					glPushMatrix()
					glTranslate(res[0][i], res[1][i], res[2][i])
					glScale(1,1,1)
					if (str(res[3][i])[0] == 'C'): self.load_gl_color("white")
					elif (str(res[3][i])[0] == 'N'): self.load_gl_color("green")
					elif (str(res[3][i])[0] == 'O'): self.load_gl_color("blue")
					elif (str(res[3][i])[0] == 'S'): self.load_gl_color("red")
					else: self.load_gl_color("silver")
					glCallList(self.spheredl)
					glPopMatrix()
		
			for k in range (0, len(self.allResidues)):
				res = self.allResidues[k]
				if res[4][0] == "ALA":
					#t1 = res[3].index('N')
					#t2 = res[3].index('CA')
					#t3 = res[3].index('C')
					#t4 = res[3].index('O')
					t1 = res[3].index('CB')

					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
				
			
				elif res[4][0] == "ARG":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG')
					t3 = res[3].index('CD')
					t4 = res[3].index('NE')
					t5 = res[3].index('CZ')
					t6 = res[3].index('NH1')
					t7 = res[3].index('NH2')


					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t2, t3)
					self.makeStick(res, t3, t4)
					self.makeStick(res, t4, t5)
					self.makeStick(res, t5, t6)
					self.makeStick(res, t5, t7)

			
				elif res[4][0] == "ASP":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG')
					t3 = res[3].index('OD1')
					t4 = res[3].index('OD2')


					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t2, t3)
					self.makeStick(res, t2, t4)
				
			
				elif res[4][0] == "ASN":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG')
					t3 = res[3].index('OD1')
					t4 = res[3].index('ND2')


					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t2, t3)
					self.makeStick(res, t2, t4)

				elif res[4][0] == "CYS":
					t1 = res[3].index('CB')
					t2 = res[3].index('SG')

					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)


				elif res[4][0] == "GLY":
					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

				elif res[4][0] == "GLN":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG')
					t3 = res[3].index('CD')
					t4 = res[3].index('OE1')
					t5 = res[3].index('NE2')


					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t2, t3)
					self.makeStick(res, t3, t4)
					self.makeStick(res, t3, t5)

				elif res[4][0] == "GLU":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG')
					t3 = res[3].index('CD')
					t4 = res[3].index('OE1')
					t5 = res[3].index('OE2')


					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t2, t3)
					self.makeStick(res, t3, t4)
					self.makeStick(res, t3, t5)

			
				elif res[4][0] == "HIS":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG')
					t3 = res[3].index('CD2')
					t4 = res[3].index('ND1')
					t5 = res[3].index('NE2')
					t6 = res[3].index('CE1')


					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t2, t3)
					self.makeStick(res, t2, t4)
					self.makeStick(res, t3, t5)
					self.makeStick(res, t5, t6)
					self.makeStick(res, t4, t6)


				elif res[4][0] == "ILE":	
					t1 = res[3].index('CB')
					t2 = res[3].index('CG1')
					t3 = res[3].index('CG2')
					t4 = res[3].index('CD1')


					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t1, t3)
					self.makeStick(res, t2, t4)

			
				elif res[4][0] == "LEU":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG')
					t3 = res[3].index('CD1')
					t4 = res[3].index('CD2')


					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t2, t3)
					self.makeStick(res, t2, t4)

				elif res[4][0] == "LYS":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG')
					t3 = res[3].index('CD')
					t4 = res[3].index('CE')
					t5 = res[3].index('NZ')


					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t2, t3)
					self.makeStick(res, t3, t4)
					self.makeStick(res, t4, t5)

			
				elif res[4][0] == "MET":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG')
					t3 = res[3].index('SD')
					t4 = res[3].index('CE')


					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t2, t3)
					self.makeStick(res, t3, t4)

				elif res[4][0] == "PHE":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG')
					t3 = res[3].index('CD1')
					t4 = res[3].index('CD2')
					t5 = res[3].index('CE1')
					t6 = res[3].index('CE2')
					t7 = res[3].index('CZ')


					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t2, t3)
					self.makeStick(res, t2, t4)
					self.makeStick(res, t3, t5)
					self.makeStick(res, t4, t6)
					self.makeStick(res, t5, t7)
					self.makeStick(res, t6, t7)

				elif res[4][0] == "PRO":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG')
					t3 = res[3].index('CD')
					t4 = res[3].index('N')


					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t2, t3)
					self.makeStick(res, t3, t4)

				elif res[4][0] == "SER":
					t1 = res[3].index('CB')
					t2 = res[3].index('OG')

					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)

			
				elif res[4][0] == "THR":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG2')
					t3 = res[3].index('OG1')

					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t1, t3)

				elif res[4][0] == "TRP":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG')
					t3 = res[3].index('CD1')
					t4 = res[3].index('CD2')
					t5 = res[3].index('NE1')
					t6 = res[3].index('CE2')
					t7 = res[3].index('CE3')
					t8 = res[3].index('CZ3')
					t9 = res[3].index('CH2')
					t10 = res[3].index('CZ2')

					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t2, t3)
					self.makeStick(res, t2, t4)
					self.makeStick(res, t3, t5)
					self.makeStick(res, t5, t6)
					self.makeStick(res, t4, t6)
					self.makeStick(res, t4, t7)
					self.makeStick(res, t7, t8)
					self.makeStick(res, t8, t9)
					self.makeStick(res, t10, t9)

				elif res[4][0] == "TYR":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG')
					t3 = res[3].index('CD1')
					t4 = res[3].index('CD2')
					t5 = res[3].index('CE1')
					t6 = res[3].index('CE2')
					t7 = res[3].index('CZ')
					t8 = res[3].index('OH')


					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t2, t3)
					self.makeStick(res, t2, t4)
					self.makeStick(res, t3, t5)
					self.makeStick(res, t4, t6)
					self.makeStick(res, t5, t7)
					self.makeStick(res, t6, t7)
					self.makeStick(res, t7, t8)

				elif res[4][0] == "VAL":
					t1 = res[3].index('CB')
					t2 = res[3].index('CG2')
					t3 = res[3].index('CG1')

					self.makeStick(res, 0, 1)
					self.makeStick(res, 1, 2)
					self.makeStick(res, 2, 3)

					self.makeStick(res, 1, t1)
					self.makeStick(res, t1, t2)
					self.makeStick(res, t1, t3)

				if (k!=0):
				
					nt = [0,0,0]
					pt = [0,0,0]
					nt[0] = res[0][0]
					nt[1] = res[1][0]
					nt[2] = res[2][0]

					pt[0] = self.allResidues[(k-1)][0][2]
					pt[1] = self.allResidues[(k-1)][1][2]
					pt[2] = self.allResidues[(k-1)][2][2]
					self.cylinder_to_from(nt, pt, 0.2)
			glEndList()

		glCallList(self.dl)


	def makeStick (self, res, index1, index2):
		n = [0,0,0]
		p = [0,0,0]
		p[0] = res[0][index1]
		p[1] = res[1][index1]
		p[2] = res[2][index1]

		n[0] = res[0][index2]
		n[1] = res[1][index2]
		n[2] = res[2][index2]
		self.cylinder_to_from(n, p, 0.2)	

	def buildResList (self):

		self.allResidues = []
		
		try:
			f = open(self.fName)
			f.close()
		except IOError:	
			print "Sorry, the file name \"" + str(self.fName) + "\" does not exist"
			sys.exit()
		
   		self.a = PDBReader()
    		self.a.read_from_pdb(self.fName)
    		point_x = self.a.get_x()
   		point_y = self.a.get_y()
	        point_z = self.a.get_z()
		point_atomName = self.a.get_atomName()
		point_resName = self.a.get_resName()
		point_resNum = self.a.get_resNum()
		x =[]
		y =[]
		z =[]
		atomName =[]
		resName = []
		amino = []
		currentRes = 1

    		for i in range(0, len(point_x)):
        		if (point_resNum[i]==currentRes):
           			x.append(point_x[i])
            			y.append(point_y[i])
            			z.append(point_z[i])
				temp = point_atomName[i]
				temp2 = temp.strip()
				atomName.append(temp2)
            			resName.append(point_resName[i])
       			else:
            			currentRes = point_resNum[i]
				amino.append(x[:])
				amino.append(y[:])
				amino.append(z[:])
				amino.append(atomName[:])
				amino.append(resName[:])
				self.allResidues.append(amino[:])
				del amino[:]
            			del x[:]
            			del y[:]
            			del z[:]
            			del atomName[:]
            			del resName[:]
           			x.append(point_x[i])
            			y.append(point_y[i])
            			z.append(point_z[i])
				temp = point_atomName[i]
				temp2 = temp.strip()
				atomName.append(temp2)
            			resName.append(point_resName[i])
			if (i == (len(point_x)-1)): 
				amino.append(x[:])
				amino.append(y[:])
				amino.append(z[:])
				amino.append(atomName[:])
				amino.append(resName[:])
				self.allResidues.append(amino[:])
				break


	def cylinder_to_from(self,next,prev,scale=0.5):

		dx = next[0] - prev[0]
		dy = next[1] - prev[1]
		dz = next[2] - prev[2]
		try:
			length = sqrt(dx**2 + dy**2 + dz**2)
		except: return
		if length == 0: return
		

		alt = acos(dz/length)*180.0/pi
		phi = atan2(dy,dx)*180.0/pi
		
		glPushMatrix()
		glTranslatef(prev[0],prev[1],prev[2] )
		glRotatef(90+phi,0,0,1)
		glRotatef(alt,1,0,0)
		glScalef(scale,scale,length)
		self.load_gl_color("silver")
		glCallList(self.cylinderdl)

		glPopMatrix()		
		
class EMPDBInspector(EM3DInspector):
	def __init__(self,target):
		EM3DInspector.__init__(self,target)
		self.tabwidget.insertTab(0,self.get_example_tab(),"Main")
		self.tabwidget.setCurrentIndex(0)
			
	def get_example_tab(self):
		'''
		@return an QtGui.QWidget, i.e. for insertion into a tab widget, or layour, etc
		'''
		widget = QtGui.QWidget()
		vbl = QtGui.QVBoxLayout(widget)
		vbl.setMargin(0)
		vbl.setSpacing(6)

		hbl1 = QtGui.QHBoxLayout()
		self.text = QtGui.QLineEdit()
		self.text.setText(self.target().current_text())
		#text_label = QtGui.QLabel("Enter Text:",self)
		#hbl1.addWidget(text_label)
		hbl1.addWidget(self.text)
		self.browse = QtGui.QPushButton("Browse")
		hbl1.addWidget(self.browse)
		vbl.addLayout(hbl1)
		
		QtCore.QObject.connect(self.text, QtCore.SIGNAL("textChanged(const QString&)"), self.on_text_change)
		QtCore.QObject.connect(self.browse, QtCore.SIGNAL("clicked(bool)"), self.on_browse)
		
		return widget
	
	def on_text_change(self,text):
		self.target().set_current_text(str(text))
		self.target().updateGL()

	def on_browse(self):
		self.fileName = QtGui.QFileDialog.getOpenFileName(self, "open file", "/home", "Text files (*.pdb)")
		if (self.fileName == ""): return
		self.target().set_current_text(str(self.fileName)) #self.target().text and self.text are what the user sees. 
		self.text.setText(self.fileName) #if self.text changes, then self.fName becomes self.text and the image regenerates
		self.target().updateGL()
	
if __name__ == '__main__':
	from emapplication import EMStandAloneApplication
	em_app = EMStandAloneApplication()
	window = EMPDBViewer()
	em_app.show()
	em_app.execute()

