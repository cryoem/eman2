#!/usr/bin/env python

# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
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
from libpyGLUtils2 import *

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
	def __init__(self,target,enable_advanced=True):
		QtGui.QWidget.__init__(self,None)
		self.target = weakref.ref(target) # prevent a strong cycle - this target object should be an EM3DModule, but that could change depending on who builds on this object
		
		self.vbl = QtGui.QVBoxLayout(self) # this is the main vbl
		
		self.advanced_tab = None
		
		self.tabwidget = QtGui.QTabWidget()
		
		#self.vbl.addLayout(self.hbl_check)
		self.vbl.addWidget(self.tabwidget)
		
		if enable_advanced:	self.insert_advanced_tab()
		
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
		if self.advanced_tab:
			self.advanced_tab.update_rotations(t3d)
	
	def set_scale(self,val):
		if self.advanced_tab: self.advanced_tab.set_scale(val)
	
	def set_xy_trans(self, x, y):
		if self.advanced_tab: self.advanced_tab.set_xy_trans(x,y)
	
	def set_xyz_trans(self,x,y,z):
		if self.advanced_tab: self.advanced_tab.set_xyz_trans(x,y,z)
		
	def set_directional_light_dir(self,d):
		if self.advanced_tab: self.advanced_tab.set_directional_light_dir(d)
	
	def set_positional_light_pos(self,d):
		if self.advanced_tab: self.advanced_tab.set_positional_light_pos(d)
		
	def set_positional_light_dir(self,d):
		if self.advanced_tab: self.advanced_tab.set_positional_light_dir(d)
	
class EM3DExampleModule(EM3DModule):
	def __init__(self,application=None):
		EM3DModule.__init__(self,application)
		self.text = "bdb:EMAN2"
		
		self.dl = None
	
	def current_text(self): return self.text
	
	def set_current_text(self,text):
		self.text = text
	
	def get_inspector(self):
		if self.inspector == None:
			self.inspector = EM3DExampleInspector(self)
		return self.inspector
		
	def draw_objects(self):
		
		glPushMatrix()
		glTranslate(100,0,0)
		glScale(50,50,50)
		self.load_gl_color("blue")
		glCallList(self.spheredl)
		glPopMatrix()
		
		glPushMatrix()
		glTranslate(0,100,0)
		glScale(50,10,10)
		glTranslate(-0.5,0,0)
		glRotate(90,0,1,0)
		self.load_gl_color("emerald")
		glCallList(self.cylinderdl)
		glPopMatrix()
		
		glPushMatrix()
		glTranslate(-100,25,0)
		glScale(10,50,10)
		glTranslate(-0.5,0,0)
		glRotate(90,1,0,0)
		self.load_gl_color("gold")
		glCallList(self.cappedcylinderdl)
		glPopMatrix()
		
		glPushMatrix()
		glTranslate(0,-100,0)
		glScale(15,15,15)
		self.load_gl_color("copper")
		glCallList(self.diskdl)
		glPopMatrix()
		
		glPushMatrix()
		glTranslate(0,-100,-10)
		glScale(15,15,15)
		glRotate(180,0,1,0)
		self.load_gl_color("silver")
		glCallList(self.diskdl)
		glPopMatrix()
		
		glPushMatrix()
		glTranslate(0,0,50)	
		bbox = self.font_renderer.bounding_box(self.text)
		glTranslate(-(bbox[3]-bbox[0])/2, -(bbox[4]-bbox[1])/2,-(bbox[5]-bbox[02])/2)
		self.font_renderer.render_string(self.text)
		glPopMatrix()
		
		
		glPushMatrix()
		glTranslate(0,-50,-50)
		s = "=^_^="
		bbox = self.font_renderer.bounding_box(s)
		glTranslate(-(bbox[3]-bbox[0])/2, -(bbox[4]-bbox[1])/2,-(bbox[5]-bbox[02])/2)
		self.font_renderer.render_string(s)
		glPopMatrix()
		
class EM3DExampleInspector(EM3DInspector):
	def __init__(self,target):
		EM3DInspector.__init__(self,target)
		self.tabwidget.insertTab(0,self.get_example_tab(),"Example")
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
		text_label = QtGui.QLabel("Enter Text:",self)
		hbl1.addWidget(text_label)
		hbl1.addWidget(self.text)
		vbl.addLayout(hbl1)
		
		QtCore.QObject.connect(self.text, QtCore.SIGNAL("textChanged(const QString&)"), self.on_text_change)
		
		return widget
	
	def on_text_change(self,text):
		self.target().set_current_text(str(text))
		self.target().updateGL()
	
if __name__ == '__main__':
	from emapplication import EMStandAloneApplication
	em_app = EMStandAloneApplication()
	window = EM3DExampleModule()
	em_app.show()
	em_app.execute()