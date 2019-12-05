#!/usr/bin/env python
from __future__ import absolute_import
#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
# and David Woolford 10/26/2007 (woolford@bcm.edu)
# and Ahmad Qamar July 2009 (qamar@uchicago.edu)
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

from past.utils import old_div
from builtins import object
from PyQt5 import QtCore, QtGui, QtWidgets, QtOpenGL
from PyQt5.QtCore import Qt
import OpenGL
OpenGL.ERROR_CHECKING = False
from OpenGL import GL,GLU,GLUT
from OpenGL.GL import *
from OpenGL.GLU import *
from .valslider import ValSlider
from math import *
from EMAN2 import *
from .emimageutil import EMTransformPanel
import weakref
from time import *
from libpyGLUtils2 import *

from .emglobjects import EM3DModel, Camera2,EMViewportDepthTools, get_default_gl_colors
from .emlights import *


class DynamicFonts(object):
	def __init__(self):
		self.font_renderer = get_3d_font_renderer()

	def set_depth(self,length):
		self.font_renderer.set_depth(length)

class EM3DFontModel(EMLightsDrawer,EM3DModel,DynamicFonts):
	def __init__(self, gl_widget):
		EM3DModel.__init__(self, gl_widget)
		DynamicFonts.__init__(self)

		self.init()
		self.initialized = True
		self.load_colors()
		self.cam=Camera2(self)
		self.cam.basicmapping = True #Ross's experiment... fixes translation

		self.brightness = 0
		self.contrast = 10
		self.glcontrast = 1.0
		self.glbrightness = 0.0
		self.rank = 1
		self.inspector=None
		self.bgR = .85
		self.bgG = .85
		self.bgB = 1.0
		self.bg_a = 1
		self.lspacing = 75
#		self.get_gl_widget().cam.default_z = -25	# this is me hacking
#		self.get_gl_widget().cam.cam_z = -25 		# this is me hacking
		self.vdtools = EMViewportDepthTools(self)
		self.font_renderer = get_3d_font_renderer()
		self.font_renderer.set_font_mode(FTGLFontMode.EXTRUDE)
		self.font_renderer.set_depth(75)
		self.render_string = "hello world"

		EMLightsDrawer.__init__(self)
		self.setInit()

	def set_render_string(self,string):
		self.render_string = string
	def set_bg_r(self,bgR):
		self.bgR = bgR
	def set_bg_g(self,bgG):
		self.bgG = bgG
	def set_bg_b(self,bgB):
		self.bgB = bgB
	def set_bg_a(self,bg_a):
		self.bg_a = bg_a
	def set_lspacing(self,lspacing):
		self.lspacing = lspacing

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

		glPushMatrix()
		self.cam.position()
	
		glPushMatrix()
		self.cam.position(True)
		# the ones are dummy variables atm... they don't do anything
		self.vdtools.update(1,1)
		glPopMatrix()

		glShadeModel(GL_SMOOTH)

		glStencilFunc(GL_EQUAL,self.rank,0)
		glStencilOp(GL_KEEP,GL_KEEP,GL_REPLACE)
		glMaterial(GL_FRONT, GL_AMBIENT, self.colors[self.currentcolor]["ambient"])
		glMaterial(GL_FRONT, GL_DIFFUSE, self.colors[self.currentcolor]["diffuse"])
		glMaterial(GL_FRONT, GL_SPECULAR, self.colors[self.currentcolor]["specular"])
		glMaterial(GL_FRONT, GL_SHININESS, self.colors[self.currentcolor]["shininess"])
		glColor(self.colors[self.currentcolor]["diffuse"])

		glClearColor(self.bgR,self.bgG,self.bgB,self.bg_a)

		glEnable(GL_NORMALIZE)
		#HERE
		glPushMatrix()
		glNormal(0,0,1)
		glEnable(GL_TEXTURE_2D)

		i = 0
		self.ifevalstr = self.render_string.split("\n")
		ifevalref = len(self.ifevalstr)-1
		spfac_i = -1*int(old_div((ifevalref+1),2))
		spfac_f = int(old_div((ifevalref+1),2))
		spfac = [-1*int(old_div((ifevalref+1),2))]
		while spfac_i<spfac_f:
			spfac_i = spfac_i+1
			spfac.append(spfac_i)
		if ifevalref%2!=0:
			spfac.remove(0)
			while i<old_div(len(spfac),2):
				i = i+1
				spfac[i-1]=spfac[i-1]+0.5
			while (i<len(spfac)):
				i = i+1
				spfac[i-1]=spfac[i-1]-0.5
		i = 0
		while i<=ifevalref:
			i = i+1
			tvar = str("bbox"+str(i))
			tvar = self.font_renderer.bounding_box(self.ifevalstr[i-1])
			glPushMatrix()
			glTranslate(old_div((tvar[0]-tvar[3]),2),old_div((tvar[1]-tvar[4]-(((spfac[i-1])*self.lspacing)-0)),2),old_div(-(tvar[2]-tvar[5]),2))
			self.font_renderer.render_string(self.ifevalstr[i-1]);
			glPopMatrix()	

		glPopMatrix()

		glPopMatrix()
		glStencilFunc(GL_EQUAL,self.rank,self.rank)
		glStencilOp(GL_KEEP,GL_KEEP,GL_KEEP)
		glPushMatrix()
		glLoadIdentity()
		glScalef(10,10,1)
		glTranslate(-90.5,-90.5,-91)
		self.draw_bc_screen()
		glPopMatrix()

		EMLightsDrawer.draw(self)		
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
		self.cam.default_z = 0
		self.cam.cam_z = -10*32
		if not self.inspector or self.inspector ==None:
			self.inspector=EMFontInspector(self)
		self.inspector.setColors(self.colors,self.currentcolor)

	def load_colors(self):
		self.colors = get_default_gl_colors()
		self.currentcolor = "ruby"
	def mouseDoubleClickEvent(self, event):
		if self.current_mouse_mode:
			EMLightsDrawer.mouseDoubleClickEvent(self, event)
		else:
			EM3DModel.mouseDoubleClickEvent(self, event)
	def mouseMoveEvent(self, event):
		if self.current_mouse_mode:
			EMLightsDrawer.mouseMoveEvent(self, event)
		else:
			EM3DModel.mouseMoveEvent(self, event)
	def mousePressEvent(self, event):
		if self.current_mouse_mode:
			EMLightsDrawer.mousePressEvent(self, event)
		else:
			EM3DModel.mousePressEvent(self, event)
	def mouseReleaseEvent(self, event):
		if self.current_mouse_mode:
			EMLightsDrawer.mouseReleaseEvent(self, event)
		else:
			EM3DModel.mouseReleaseEvent(self, event)
	def setColor(self,val):
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
	
	def eye_coords_dif(self,x1,y1,x2,y2,mdepth=True):
		return self.vdtools.eye_coords_dif(x1,y1,x2,y2,mdepth)

#	def resize(self):
#		self.vdtools.set_update_P_inv()
	def get_type(self):
		return "EM3DFontModel"


class EMFontInspector(QtWidgets.QWidget, EMLightsInspectorBase):
	def __init__(self,target) :
		QtWidgets.QWidget.__init__(self,None)
		EMLightsInspectorBase.__init__(self)
		self.target=weakref.ref(target)
		self.transform_panel = EMTransformPanel(target,self)
		self.transform_vbl = None # This will eventually be a vertical box layout for the transform panel
		self.init_fonts()

		self.vbl = QtWidgets.QVBoxLayout(self)
		self.vbl.setContentsMargins(0, 0, 0, 0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")

		self.hbl = QtWidgets.QHBoxLayout()
		self.hbl.setContentsMargins(0, 0, 0, 0)
		self.hbl.setSpacing(6)
		self.hbl.setObjectName("hbl")
		self.vbl.addLayout(self.hbl)

		self.vbl2 = QtWidgets.QVBoxLayout()
		self.vbl2.setContentsMargins(0, 0, 0, 0)
		self.vbl2.setSpacing(6)
		self.vbl2.setObjectName("vbl2")
		self.hbl.addLayout(self.vbl2)

		self.wiretog = QtWidgets.QPushButton("Wire")
		self.wiretog.setCheckable(1)
		self.vbl2.addWidget(self.wiretog)

		self.lighttog = QtWidgets.QPushButton("Light")
		self.lighttog.setCheckable(1)
		self.vbl2.addWidget(self.lighttog)

		self.tabwidget2 = QtWidgets.QTabWidget()
		self.maintab = None
		self.tabwidget2.addTab(self.get_main_tab(), "Main")
		#self.tabwidget2.addTab(self.get_GL_tab(),"GL")
		self.tabwidget2.addTab(self.get_format_tab(),"Formatting")
		self.tabwidget2.addTab(self.get_light_tab(), "Lights")
		self.vbl.addWidget(self.tabwidget2)
		self.n3_showing = False

		self.cbb.currentIndexChanged[str].connect(target.setColor)
		self.wiretog.toggled[bool].connect(target.toggle_wire)
		self.lighttog.toggled[bool].connect(target.toggle_light)
		self.glcontrast.valueChanged.connect(target.set_GL_contrast)
		self.glbrightness.valueChanged.connect(target.set_GL_brightness)
		self.combo.currentIndexChanged [str].connect(self.on_combo_change)
		self.text.textChanged[str].connect(self.on_text_change)
		self.lspacing.valueChanged.connect(self.set_GL_lspacing)
		self.length.valueChanged.connect(self.set_GL_length)
		self.tsize.valueChanged[int].connect(self.set_GL_tsize)
		self.Dfont.currentIndexChanged [str].connect(self.on_Dfont_change)
		self.bgR.valueChanged.connect(self.set_GL_bgR)
		self.bgG.valueChanged.connect(self.set_GL_bgG)
		self.bgB.valueChanged.connect(self.set_GL_bgB)
		self.bg_a.valueChanged.connect(self.set_GL_bg_a)
	
	def get_transform_layout(self):
		return self.transform_vbl
		
	def set_GL_bgR(self,bgR):
		self.target().set_bg_r(bgR)
		self.target().updateGL()

	def set_GL_bgG(self,bgG):
		self.target().set_bg_g(bgG)
		self.target().updateGL()

	def set_GL_bgB(self,bgB):
		self.target().set_bg_b(bgB)
		self.target().updateGL()

	def set_GL_bg_a(self,bg_a):
		self.target().set_bg_a(bg_a)
		self.target().updateGL()

	def init_fonts(self):
		self.d = {}
		self.l = []
		platform = get_platform()
		if platform == "Linux":
			f_dir = "/usr/share/fonts/"
		elif platform == "Windows" or platform == "win32":
			f_dir = ":/windows/fonts/"
		elif platform in ["Apple", "Darwin"]:
			f_dir = "/Library/Fonts/"
		else:
			raise RuntimeError("Platform %s is not supported" %platform )
		
		for root, dirs, files in os.walk(f_dir):
			for name in files:
				if name.find("ttf")!=-1:
					filename = os.path.join(root, name)
					self.d[name] = filename
					self.l.extend([name])
		return self.l, self.d

	def on_Dfont_change(self,Dfont):
		self.target().font_renderer.set_font_file_name(self.d[str(Dfont)])
		self.target().updateGL()

	def set_GL_lspacing(self,lspacing):
		self.target().set_lspacing(lspacing)
		#THE FOLLOWING IF STATEMENT DOES IS NOT EFFECTIVE
		if len(self.target().render_string.split("\n")) != 1:
			self.lspacing.setEnabled(True)
		else:
			self.lspacing.setEnabled(False)
		self.target().updateGL()

	def set_GL_length(self,length):
		self.target().font_renderer.set_depth(int(length))	
		self.target().updateGL()

	def set_GL_tsize(self,tsize):
		self.target().font_renderer.set_face_size(tsize)
		self.target().updateGL()

	def on_text_change(self,text):
		try:
			evalt=str(eval(str(text)))
			self.target().set_render_string(evalt)
		except:
			self.target().set_render_string(str(text))
			
		if len(self.target().render_string.split("\n")) != 1:
			self.lspacing.setEnabled(True)
		else:
			self.lspacing.setEnabled(False)
		self.target().updateGL()

	def on_combo_change(self,mode):
		d = {}
		d["Extrude"] = FTGLFontMode.EXTRUDE
		d["Pixmap"] = FTGLFontMode.PIXMAP
		d["Bitmap"] = FTGLFontMode.BITMAP
		d["Polygon"] = FTGLFontMode.POLYGON
		d["Outline"] = FTGLFontMode.OUTLINE
		d["Texture"] = FTGLFontMode.TEXTURE
		self.target().font_renderer.set_font_mode(d[str(mode)])
		if mode == "Extrude":
			self.length.setEnabled(True)
		else:
			self.length.setEnabled(False)
		self.target().updateGL()

	def update_rotations(self,t3d):
		self.transform_panel.update_rotations(t3d)

	def set_scale(self,val):
		self.transform_panel.set_scale(val)

	def set_xy_trans(self, x, y):
		self.transform_panel.set_xy_trans(x,y)

	def set_xyz_trans(self,x,y,z):
		self.transform_panel.set_xyz_trans(x,y,z)	

	def get_GL_tab(self):
		self.gltab = QtWidgets.QWidget()
		gltab = self.gltab

		gltab.vbl = QtWidgets.QVBoxLayout(self.gltab)
		gltab.vbl.setContentsMargins(0, 0, 0, 0)
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
			self.maintab = QtWidgets.QWidget()
			maintab = self.maintab
			maintab.vbl = QtWidgets.QVBoxLayout(self.maintab)
			maintab.vbl.setContentsMargins(0, 0, 0, 0)
			maintab.vbl.setSpacing(6)
			maintab.vbl.setObjectName("Main")
			
			self.transform_vbl = QtWidgets.QVBoxLayout()
			self.transform_panel.addWidgets(self.transform_vbl)
			maintab.vbl.addLayout(self.transform_vbl)
			self.glwidget = QtWidgets.QTabWidget()
			self.glwidget.addTab(self.get_GL_tab(),"GL")
			maintab.vbl.addWidget(self.glwidget)

		return maintab

	def get_format_tab(self):
		self.formattab = QtWidgets.QWidget()
		formattab = self.formattab
		formattab.vbl = QtWidgets.QVBoxLayout(self.formattab)
		formattab.vbl.setContentsMargins(0, 0, 0, 0)
		formattab.vbl.setSpacing(6)
		formattab.vbl.setObjectName("Format")

		self.hbl1 = QtWidgets.QHBoxLayout()
		self.text = QtWidgets.QLineEdit()
		self.text.setText("hello world")
		text_label = QtWidgets.QLabel("Enter Text:",self)
		text_label.setToolTip("Enters quotes to evaluate new line e.g. \"hello\\nworld\". Evaluates numerical expressions e.g. 9*9 (with out quotes)")
		self.hbl1.addWidget(text_label)
		self.hbl1.addWidget(self.text)
		formattab.vbl.addLayout(self.hbl1)

		self.hbl1 = QtWidgets.QHBoxLayout()
		self.Dfont = QtWidgets.QComboBox()
		for k in self.l: self.Dfont.addItem(k)
		self.hbl1.addWidget(QtWidgets.QLabel("Fonts:",self))
		self.hbl1.addWidget(self.Dfont)
		formattab.vbl.addLayout(self.hbl1)

		self.hbl1 = QtWidgets.QHBoxLayout()
		self.tsize = QtWidgets.QSpinBox()
		self.tsize.setRange(0,500)
		self.tsize.setValue(32)
		self.hbl1.addWidget(QtWidgets.QLabel("Size:",self),Qt.AlignLeft)
		self.hbl1.addWidget(self.tsize,Qt.AlignRight)
		self.combo = QtWidgets.QComboBox()
		self.items = ["Extrude","Pixmap","Bitmap","Polygon","Outline","Texture"]
		for k in self.items: self.combo.addItem(k)
		self.hbl1.addWidget(QtWidgets.QLabel("Style:",self),Qt.AlignLeft)
		self.hbl1.addWidget(self.combo,Qt.AlignRight)
		formattab.vbl.addLayout(self.hbl1)

		self.hbl1 = QtWidgets.QHBoxLayout()
		self.lspacing = ValSlider(self,(-100.0,100.0),"Line Spacing:")
		self.lspacing.setObjectName("Length")
		self.lspacing.setValue(75.0)
		self.lspacing.setEnabled(False)
		self.hbl1.addWidget(self.lspacing)
		formattab.vbl.addLayout(self.hbl1)

		self.hbl1 = QtWidgets.QHBoxLayout()
		self.length = ValSlider(self,(0.0,500.0),"Length:")
		self.length.setObjectName("Length")
		self.length.setValue(75.0)
		self.hbl1.addWidget(self.length)
		formattab.vbl.addLayout(self.hbl1)

		self.hbl1 = QtWidgets.QHBoxLayout()
		self.cbb = QtWidgets.QComboBox()
		self.hbl1.addWidget(QtWidgets.QLabel("Material:",self))
		self.hbl1.addWidget(self.cbb)
		formattab.vbl.addLayout(self.hbl1)

		self.hbl1 = QtWidgets.QHBoxLayout()
		self.bgtabwidget = QtWidgets.QTabWidget()
		self.maintab = None
		self.bgtabwidget.addTab(self.get_bgRGB_tab(), "BG RGB")
		self.hbl1.addWidget(self.bgtabwidget)
		self.n3_showing = False
		formattab.vbl.addLayout(self.hbl1)

		return formattab

	def get_bgRGB_tab(self):
		self.bgRGBtab = QtWidgets.QWidget()
		bgRGBtab = self.bgRGBtab
		bgRGBtab.vbl2 = QtWidgets.QVBoxLayout(self.bgRGBtab)
		bgRGBtab.vbl2.setContentsMargins(0, 0, 0, 0)
		bgRGBtab.vbl2.setSpacing(6)
		bgRGBtab.vbl2.setObjectName("BG RGB")

		self.hbl2 = QtWidgets.QHBoxLayout()
		self.bgR = ValSlider(self,(0,1),"R:")
		self.bgR.setObjectName("R")
		self.bgR.setValue(0.5)
		self.hbl2.addWidget(self.bgR)
		bgRGBtab.vbl2.addLayout(self.hbl2)

		self.hbl2 = QtWidgets.QHBoxLayout()
		self.bgG = ValSlider(self,(0,1),"G:")
		self.bgG.setObjectName("G")
		self.bgG.setValue(0.5)
		self.hbl2.addWidget(self.bgG)
		bgRGBtab.vbl2.addLayout(self.hbl2)

		self.hbl2 = QtWidgets.QHBoxLayout()
		self.bgB = ValSlider(self,(0,1),"B:")
		self.bgB.setObjectName("B")
		self.bgB.setValue(0.5)
		self.hbl2.addWidget(self.bgB)
		bgRGBtab.vbl2.addLayout(self.hbl2)		

		self.hbl2 = QtWidgets.QHBoxLayout()
		self.bg_a = ValSlider(self,(0,1),"Alpha:")
		self.bg_a.setObjectName("Alpha")
		self.bg_a.setValue(1.0)
		self.hbl2.addWidget(self.bg_a)
		bgRGBtab.vbl2.addLayout(self.hbl2)

		return bgRGBtab

	def setColors(self,colors,current_color):
		a = 0
		for i in colors:
			self.cbb.addItem(i)
			if ( i == current_color):
				self.cbb.setCurrentIndex(a)
			a += 1

# This is just for testing, of course
def main():
	from .emapplication import EMApp
	from .emimage3d import EMImage3DWidget
	em_app = EMApp()
	window = EMImage3DWidget()
	font_model = EM3DFontModel(window)
	window.add_model(font_model)
	window.cam.default_z = -25	# From David's "this is me hacking"
	window.cam.cam_z = -25 		# From David's "this is me hacking"
	em_app.show()
	em_app.execute()


if __name__ == '__main__':
	main()
