#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
# and David Woolford 10/26/2007 (woolford@bcm.edu)
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
from valslider import ValSlider
from math import *
from EMAN2 import *
import sys
import numpy
from emimageutil import ImgHistogram,EMParentWin,EMTransformPanel
from weakref import WeakKeyDictionary
from time import time
from PyQt4.QtCore import QTimer
import weakref
from time import *

from emglobjects import EMImage3DGUIModule, Camera2,EMViewportDepthTools, EMGLProjectionViewMatrices, get_default_gl_colors

MAG_INCREMENT_FACTOR = 1.1

class EM3DFontWidget(EMImage3DGUIModule):
	def __init__(self):
		EMImage3DGUIModule.__init__(self)
		#self.parent = parent
		
		self.init()
		self.initialized = True
		
		self.cam=Camera2(self)
		
		self.brightness = 0
		self.contrast = 10
		self.glcontrast = 1.0
		self.glbrightness = 0.0
		self.rank = 1
		self.inspector=None
		
		self.vdtools = EMViewportDepthTools(self)
		self.font_renderer = get_3d_font_renderer()
		#self.font_renderer.set_font_mode(FTGLFontMode.POLYGON) # or EXTRUDE, PIXMAP, BITMAP, POLYGON or OUTLINE
		self.font_renderer.set_font_mode(FTGLFontMode.EXTRUDE) # or EXTRUDE, PIXMAP, BITMAP, POLYGON or OUTLINE
		self.font_renderer.set_depth(100)
		self.setInit()
		#self.font_renderer.set_font_mode(FTGLFontMode.OUTLINE)
		
	def get_type(self):
		return "Font"

	def render(self):
		#if (not isinstance(self.data,EMData)): return
		
		
		glEnable(GL_NORMALIZE)
		lighting = glIsEnabled(GL_LIGHTING)
		cull = glIsEnabled(GL_CULL_FACE)
		depth = glIsEnabled(GL_DEPTH_TEST)
		polygonmode = glGetIntegerv(GL_POLYGON_MODE)

		glDisable(GL_CULL_FACE)
		glEnable(GL_DEPTH_TEST)
		
		if ( self.wire ):
			glPolygonMode(GL_FRONT_AND_BACK,GL_LINE)
		else:
			glPolygonMode(GL_FRONT_AND_BACK,GL_FILL)
		
		if self.light:
			glEnable(GL_LIGHTING)
		else:
			glDisable(GL_LIGHTING)

		glPushMatrix()
		self.cam.position(True)
		# the ones are dummy variables atm... they don't do anything
		self.vdtools.update(1,1)
		glPopMatrix()

		self.cam.position()
			
		glShadeModel(GL_SMOOTH)

		glStencilFunc(GL_EQUAL,self.rank,0)
		glStencilOp(GL_KEEP,GL_KEEP,GL_REPLACE)
		glMaterial(GL_FRONT, GL_AMBIENT, self.colors[self.currentcolor]["ambient"])
		glMaterial(GL_FRONT, GL_DIFFUSE, self.colors[self.currentcolor]["diffuse"])
		glMaterial(GL_FRONT, GL_SPECULAR, self.colors[self.currentcolor]["specular"])
		glMaterial(GL_FRONT, GL_SHININESS, self.colors[self.currentcolor]["shininess"])
		glColor(self.colors[self.currentcolor]["ambient"])
		
		glEnable(GL_NORMALIZE)
		#HERE
		glPushMatrix()
		glNormal(0,0,1)
		glEnable(GL_TEXTURE_2D)
		bbox = self.font_renderer.bounding_box("hello world")
		glTranslate((bbox[0]-bbox[3])/2,(bbox[1]-bbox[4])/2,-(bbox[2]-bbox[5])/2)
		self.font_renderer.render_string("hello world");

		glPopMatrix()
		
		glStencilFunc(GL_EQUAL,self.rank,self.rank)
		glStencilOp(GL_KEEP,GL_KEEP,GL_KEEP)
		glPushMatrix()
		glLoadIdentity()
		glScalef(10,10,1)
		glTranslate(-0.5,-0.5,-1)
		self.draw_bc_screen()
		glPopMatrix()
		
		glStencilFunc(GL_ALWAYS,1,1)

			
		if ( lighting ): glEnable(GL_LIGHTING)
		else: glDisable(GL_LIGHTING)
		if ( cull ): glEnable(GL_CULL_FACE)
		else: glDisable(GL_CULL_FACE)
		if ( depth ): glEnable(GL_DEPTH_TEST)
		else : glDisable(GL_DEPTH_TEST)
		
		if ( polygonmode[0] == GL_LINE ): glPolygonMode(GL_FRONT, GL_LINE)
		else: glPolygonMode(GL_FRONT, GL_FILL)
		if ( polygonmode[1] == GL_LINE ): glPolygonMode(GL_BACK, GL_LINE)
		else: glPolygonMode(GL_BACK, GL_FILL)
			
	def init(self):
		self.mmode = 0
		self.wire = False
		self.light = True
	
	def setInit(self):

		self.cam.default_z = -10*32
		self.cam.cam_z = -10*32
		
		if not self.inspector or self.inspector ==None:
			self.inspector=EMFontInspector(self)
		
		self.load_colors()
		self.inspector.setColors(self.colors,self.currentcolor)
	def load_colors(self):
		self.colors = get_default_gl_colors()
#		ruby = {}
#		ruby["ambient"] = [0.1745, 0.01175, 0.01175,1.0]
#		ruby["diffuse"] = [0.61424, 0.04136, 0.04136,1.0]
#		ruby["specular"] = [0.927811, 0.826959, 0.826959,1.0]
#		ruby["shininess"] = 128.0
#		
#		emerald = {}
#		emerald["ambient"] = [0.0215, 0.1745, 0.0215,1.0]
#		emerald["diffuse"] = [0.07568, 0.61424,  0.07568,1.0]
#		emerald["specular"] = [0.833, 0.927811, 0.833,1.0]
#		emerald["shininess"] = 128.0
#		
#		pearl = {}
#		pearl["ambient"] = [0.25, 0.20725, 0.20725,1.0]
#		pearl["diffuse"] = [1.0, 0.829, 0.829,1.0]
#		pearl["specular"] = [0.296648, 0.296648, 0.296648,1.0]
#		pearl["shininess"] = 128.0
#		
#		silver = {}
#		silver["ambient"] = [0.25, 0.25, 0.25,1.0]
#		silver["diffuse"] = [0.4, 0.4, 0.4,1.0]
#		silver["specular"] = [0.974597, 0.974597, 0.974597,1.0]
#		silver["shininess"] = 128.0
#		
#		gold = {}
#		gold["ambient"] = [0.24725, 0.2245, 0.0645,1.0]
#		gold["diffuse"] = [0.34615, 0.3143, 0.0903,1.0]
#		gold["specular"] = [1.000, 0.9079885, 0.26086934,1.0]
#		gold["shininess"] = 128.0
#		
#		copper = {}
#		copper["ambient"] = [0.2295, 0.08825, 0.0275,1.0]
#		copper["diffuse"] = [0.5508, 0.2118, 0.066,1.0]
#		copper["specular"] = [0.9, 0.5, 0.2,1.0]
#		copper["shininess"] = 128.0
#		
#		obsidian = {}
#		obsidian["ambient"] = [0.05375,  0.05,     0.06625 ,1.0]
#		obsidian["diffuse"] = [0.18275,  0.17,     0.22525,1.0]
#		obsidian["specular"] = [0.66, 0.65, 0.69]
#		obsidian["shininess"] = 128.0
#		
#		turquoise = {}
#		turquoise["ambient"] = [0.1, 0.18725, 0.1745 ,1.0]
#		turquoise["diffuse"] = [0.396, 0.74151, 0.69102,1.0]
#		turquoise["specular"] = [0.297254, 0.30829, 0.306678]
#		turquoise["shininess"] = 128.0
#		
#		yellow = {}
#		yellow["ambient"] = [0.3, 0.3, 0.0,1]
#		yellow["diffuse"] = [0.5, 0.5, 0.0,1]
#		yellow["specular"] = [0.7, 0.7, 0.0,1]
#		yellow["shininess"] =  60
#		
#		self.colors = {}
#		self.colors["ruby"] = ruby
#		self.colors["emerald"] = emerald
#		self.colors["pearl"] = pearl
#		self.colors["silver"] = silver
#		self.colors["gold"] = gold
#		self.colors["copper"] = copper
#		self.colors["obsidian"] = obsidian
#		self.colors["turquoise"] = turquoise
#		self.colors["yellow"] = yellow
		
		self.currentcolor = "turquoise"
	
	def setColor(self,val):
		#print val
		self.currentcolor = str(val)
		self.updateGL()
	
	def toggle_wire(self,val):
		self.wire = not self.wire
		self.updateGL()
		
	def toggle_light(self,val):
		self.light = not self.light
		self.updateGL()
	
	def update_inspector(self,t3d):
		if not self.inspector or self.inspector ==None:
			self.inspector=EMFontInspector(self)
		self.inspector.update_rotations(t3d)
	
	def get_inspector(self):
		if not self.inspector : self.inspector=EMFontInspector(self)
		return self.inspector
		
#	def mousePressEvent(self, event):
##		lc=self.scrtoimg((event.x(),event.y()))
#		if event.button()==Qt.MidButton:
#			if not self.inspector or self.inspector ==None:
#				return
#			self.inspector.update_rotations(self.cam.t3d_stack[len(self.cam.t3d_stack)-1])
#			self.resizeEvent()
#			self.show_inspector(1)
#		else:
#			self.cam.mousePressEvent(event)
#		
#		self.updateGL()
		
#	def mouseMoveEvent(self, event):
#		self.cam.mouseMoveEvent(event)
#		self.updateGL()
#	
#	def mouseReleaseEvent(self, event):
#		self.cam.mouseReleaseEvent(event)
#		self.updateGL()
#			
#	def wheelEvent(self, event):
#		self.cam.wheelEvent(event)
#		self.updateGL()
#		
	def eye_coords_dif(self,x1,y1,x2,y2,mdepth=True):
		return self.vdtools.eye_coords_dif(x1,y1,x2,y2,mdepth)
	
	def resizeEvent(self,width=0,height=0):
		self.vdtools.set_update_P_inv()

class EMFontInspector(QtGui.QWidget):
	def __init__(self,target) :
		QtGui.QWidget.__init__(self,None)
		self.target=target
		self.transform_panel = EMTransformPanel(target,self)
		
		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		
		self.hbl = QtGui.QHBoxLayout()
		self.hbl.setMargin(0)
		self.hbl.setSpacing(6)
		self.hbl.setObjectName("hbl")
		self.vbl.addLayout(self.hbl)

		self.vbl2 = QtGui.QVBoxLayout()
		self.vbl2.setMargin(0)
		self.vbl2.setSpacing(6)
		self.vbl2.setObjectName("vbl2")
		self.hbl.addLayout(self.vbl2)
		
		self.wiretog = QtGui.QPushButton("Wire")
		self.wiretog.setCheckable(1)
		self.vbl2.addWidget(self.wiretog)
		
		self.lighttog = QtGui.QPushButton("Light")
		self.lighttog.setCheckable(1)
		self.vbl2.addWidget(self.lighttog)
		
		self.tabwidget = QtGui.QTabWidget()
		self.maintab = None
		self.tabwidget.addTab(self.get_main_tab(), "Main")
		self.tabwidget.addTab(self.get_GL_tab(),"GL")
		self.vbl.addWidget(self.tabwidget)
		self.n3_showing = False
		
		
		QtCore.QObject.connect(self.cbb, QtCore.SIGNAL("currentIndexChanged(QString)"), target.setColor)
		QtCore.QObject.connect(self.wiretog, QtCore.SIGNAL("toggled(bool)"), target.toggle_wire)
		QtCore.QObject.connect(self.lighttog, QtCore.SIGNAL("toggled(bool)"), target.toggle_light)
		QtCore.QObject.connect(self.glcontrast, QtCore.SIGNAL("valueChanged"), target.set_GL_contrast)
		QtCore.QObject.connect(self.glbrightness, QtCore.SIGNAL("valueChanged"), target.set_GL_brightness)
	
	def update_rotations(self,t3d):
		self.transform_panel.update_rotations(t3d)
	
	def set_scale(self,val):
		self.transform_panel.set_scale(val)
	
	def set_xy_trans(self, x, y):
		self.transform_panel.set_xy_trans(x,y)
		
	def set_xyz_trans(self,x,y,z):
		self.transform_panel.set_xyz_trans(x,y,z)	
	
	def get_GL_tab(self):
		self.gltab = QtGui.QWidget()
		gltab = self.gltab
		
		gltab.vbl = QtGui.QVBoxLayout(self.gltab )
		gltab.vbl.setMargin(0)
		gltab.vbl.setSpacing(6)
		gltab.vbl.setObjectName("Main")
		
		self.glcontrast = ValSlider(gltab,(1.0,5.0),"GLShd:")
		self.glcontrast.setObjectName("GLShade")
		self.glcontrast.setValue(1.0)
		gltab.vbl.addWidget(self.glcontrast)
		
		self.glbrightness = ValSlider(gltab,(-1.0,0.0),"GLBst:")
		self.glbrightness.setObjectName("GLBoost")
		self.glbrightness.setValue(0.1)
		self.glbrightness.setValue(0.0)
		gltab.vbl.addWidget(self.glbrightness)
	
		return gltab
	
	def get_main_tab(self):
		if ( self.maintab == None ):
			self.maintab = QtGui.QWidget()
			maintab = self.maintab
			maintab.vbl = QtGui.QVBoxLayout(self.maintab)
			maintab.vbl.setMargin(0)
			maintab.vbl.setSpacing(6)
			maintab.vbl.setObjectName("Main")
			
			self.hbl_color = QtGui.QHBoxLayout()
			self.hbl_color.setMargin(0)
			self.hbl_color.setSpacing(6)
			self.hbl_color.setObjectName("Material")
			maintab.vbl.addLayout(self.hbl_color)
			
			self.color_label = QtGui.QLabel()
			self.color_label.setText('Material')
			self.hbl_color.addWidget(self.color_label)
			
			self.cbb = QtGui.QComboBox(maintab)
			self.hbl_color.addWidget(self.cbb)
	
			self.transform_panel.addWidgets(maintab.vbl)
		
		return maintab
	
#	def slider_rotate(self):
#		self.target.load_rotation(self.get_current_rotation())

#	def set_xy_trans(self, x, y):
#		self.x_trans.setValue(x)
#		self.y_trans.setValue(y)
	
#	def set_translate_scale(self, xscale,yscale,zscale):
#		self.x_trans.setSingleStep(xscale)
#		self.y_trans.setSingleStep(yscale)
#		self.z_trans.setSingleStep(zscale)

#	def update_rotations(self,t3d):
#		rot = t3d.get_rotation(self.src_map[str(self.src.itemText(self.src.currentIndex()))])
#		
#		convention = self.src.currentText()
#		if ( self.src_map[str(convention)] == EULER_SPIN ):
#			self.n3.setValue(rot[self.n3.getLabel()],True)
#		
#		self.az.setValue(rot[self.az.getLabel()],True)
#		self.alt.setValue(rot[self.alt.getLabel()],True)
#		self.phi.setValue(rot[self.phi.getLabel()],True)
	
#	def slider_rotate(self):
#		self.target.load_rotation(self.get_current_rotation())
#	
#	def get_current_rotation(self):
#		convention = self.src.currentText()
#		rot = {}
#		if ( self.current_src == EULER_SPIN ):
#			rot[self.az.getLabel()] = self.az.getValue()
#			
#			n1 = self.alt.getValue()
#			n2 = self.phi.getValue()
#			n3 = self.n3.getValue()
#			
#			norm = sqrt(n1*n1 + n2*n2 + n3*n3)
#			
#			n1 /= norm
#			n2 /= norm
#			n3 /= norm
#			
#			rot[self.alt.getLabel()] = n1
#			rot[self.phi.getLabel()] = n2
#			rot[self.n3.getLabel()] = n3
#			
#		else:
#			rot[self.az.getLabel()] = self.az.getValue()
#			rot[self.alt.getLabel()] = self.alt.getValue()
#			rot[self.phi.getLabel()] = self.phi.getValue()
#		
#		return Transform3D(self.current_src, rot)
	
	def setColors(self,colors,current_color):
		a = 0
		for i in colors:
			self.cbb.addItem(i)
			if ( i == current_color):
				self.cbb.setCurrentIndex(a)
			a += 1

#	def set_scale(self,newscale):
#		self.scale.setValue(newscale)
#		
# This is just for testing, of course
if __name__ == '__main__':
	from emapplication import EMStandAloneApplication
	em_app = EMStandAloneApplication()
	window = EM3DFontWidget()
#	window.setInit()
	em_app.show()
	em_app.execute()
