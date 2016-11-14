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
from time import *
from emglobjects import EM3DModel, Camera2,get_default_gl_colors,EMViewportDepthTools2

MAG_INCREMENT_FACTOR = 1.1

class EM3DHelloWorld(EM3DModel):
	def eye_coords_dif(self,x1,y1,x2,y2,mdepth=True):
		return self.vdtools.eye_coords_dif(x1,y1,x2,y2,mdepth)
	
	def __init__(self, gl_widget):
		EM3DModel.__init__(self, gl_widget)
		
		self.init()
		self.initialized = True
		
		self.cam=Camera2(self)
		
		self.brightness = 0
		self.contrast = 10
		self.glcontrast = 1.0
		self.glbrightness = 0.0
		self.rank = 1
		self.inspector=None
		self.colors = get_default_gl_colors()
		self.currentcolor = "bluewhite"
		
		self.vdtools = EMViewportDepthTools2(self.get_gl_widget())
		
		self.get_gl_widget().cam.default_z = -1.25*50 # this is me hacking
		self.get_gl_widget().cam.cam_z = -1.25*50 # this is me hacking
	
		self.dl = None # display list id
	
	def get_type(self):
		return "helloworld"
	
	def render(self):
		#if (not isinstance(self.data,EMData)): return
		lighting = glIsEnabled(GL_LIGHTING)
		cull = glIsEnabled(GL_CULL_FACE)
		depth = glIsEnabled(GL_DEPTH_TEST)
		polygonmode = glGetIntegerv(GL_POLYGON_MODE)

		glDisable(GL_CULL_FACE)
		glEnable(GL_DEPTH_TEST)
		
		if ( self.wire ):
			glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
		else:
			glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
		
		if self.light:
			glEnable(GL_LIGHTING)
		else:
			glDisable(GL_LIGHTING)

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
		
		if self.dl == None:
			self.dl = glGenLists(1)
			glNewList(self.dl,GL_COMPILE)
		
			glPushMatrix()
			glScalef(10,10,1)
			
			glBegin(GL_QUADS)
			
			n = 15.0
			d = 2.0/(n+1)
			for i in range(0,int(n)):
				for j in range(0,int(n)):
					
					x1 = -1.0 + j*d
					x2 = -1.0 + (j+1)*d
					y1 = -1.0 + i*d
					y2 = -1.0 + (i+1)*d
					d1 = d*(j**2)
					d2 = d*((j+1)**2)
					glNormal(j*d,0,1)
					glVertex(x1,y1,d1)
					glNormal((j+1)*d,.01,1)
					glVertex(x2,y1,d2)
					glNormal((j+1)*d,.01,1)
					glVertex(x2,y2,d2)
					glNormal(j*d,-.01,1)
					glVertex(x1,y2,d1)
				
			glEnd()
			glPopMatrix()
			
			glEndList()
			
		glCallList(self.dl)
			
		self.draw_bc_screen()


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

		self.cam.default_z = -1.25*32
		self.cam.cam_z = -1.25*32
		
		if not self.inspector or self.inspector ==None:
			self.inspector=EMHelloWorldInspector(self)
		
		self.load_colors()
		self.inspector.setColors(self.colors,self.currentcolor)
	
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
			self.inspector=EMHelloWorldInspector(self)
		self.inspector.update_rotations(t3d)
	
	def get_inspector(self):
		if not self.inspector : 
			self.inspector=EMHelloWorldInspector(self)
			self.inspector.setColors(self.colors,self.currentcolor)
		return self.inspector
		

class EMHelloWorldInspector(QtGui.QWidget):
	def __init__(self,target) :
		QtGui.QWidget.__init__(self,None)
		self.target=target
		
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
		
		QtCore.QObject.connect(self.scale, QtCore.SIGNAL("valueChanged"), target.set_scale)
		QtCore.QObject.connect(self.az, QtCore.SIGNAL("valueChanged"), self.slider_rotate)
		QtCore.QObject.connect(self.alt, QtCore.SIGNAL("valueChanged"), self.slider_rotate)
		QtCore.QObject.connect(self.phi, QtCore.SIGNAL("valueChanged"), self.slider_rotate)
		QtCore.QObject.connect(self.cbb, QtCore.SIGNAL("currentIndexChanged(QString)"), target.setColor)
		QtCore.QObject.connect(self.src, QtCore.SIGNAL("currentIndexChanged(QString)"), self.set_src)
		QtCore.QObject.connect(self.x_trans, QtCore.SIGNAL("valueChanged(double)"), target.set_cam_x)
		QtCore.QObject.connect(self.y_trans, QtCore.SIGNAL("valueChanged(double)"), target.set_cam_y)
		QtCore.QObject.connect(self.z_trans, QtCore.SIGNAL("valueChanged(double)"), target.set_cam_z)
		QtCore.QObject.connect(self.wiretog, QtCore.SIGNAL("toggled(bool)"), target.toggle_wire)
		QtCore.QObject.connect(self.lighttog, QtCore.SIGNAL("toggled(bool)"), target.toggle_light)
		QtCore.QObject.connect(self.glcontrast, QtCore.SIGNAL("valueChanged"), target.set_GL_contrast)
		QtCore.QObject.connect(self.glbrightness, QtCore.SIGNAL("valueChanged"), target.set_GL_brightness)
	
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
			
			self.scale = ValSlider(maintab,(0.01,30.0),"Zoom:")
			self.scale.setObjectName("scale")
			self.scale.setValue(1.0)
			maintab.vbl.addWidget(self.scale)
			
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
	
			self.hbl_trans = QtGui.QHBoxLayout()
			self.hbl_trans.setMargin(0)
			self.hbl_trans.setSpacing(6)
			self.hbl_trans.setObjectName("Trans")
			maintab.vbl.addLayout(self.hbl_trans)
			
			self.x_label = QtGui.QLabel()
			self.x_label.setText('x')
			self.hbl_trans.addWidget(self.x_label)
			
			self.x_trans = QtGui.QDoubleSpinBox(self)
			self.x_trans.setMinimum(-10000)
			self.x_trans.setMaximum(10000)
			self.x_trans.setValue(0.0)
			self.hbl_trans.addWidget(self.x_trans)
			
			self.y_label = QtGui.QLabel()
			self.y_label.setText('y')
			self.hbl_trans.addWidget(self.y_label)
			
			self.y_trans = QtGui.QDoubleSpinBox(maintab)
			self.y_trans.setMinimum(-10000)
			self.y_trans.setMaximum(10000)
			self.y_trans.setValue(0.0)
			self.hbl_trans.addWidget(self.y_trans)
			
			
			self.z_label = QtGui.QLabel()
			self.z_label.setText('z')
			self.hbl_trans.addWidget(self.z_label)
			
			self.z_trans = QtGui.QDoubleSpinBox(maintab)
			self.z_trans.setMinimum(-10000)
			self.z_trans.setMaximum(10000)
			self.z_trans.setValue(0.0)
			self.hbl_trans.addWidget(self.z_trans)
			
			self.hbl_src = QtGui.QHBoxLayout()
			self.hbl_src.setMargin(0)
			self.hbl_src.setSpacing(6)
			self.hbl_src.setObjectName("hbl")
			maintab.vbl.addLayout(self.hbl_src)
			
			self.label_src = QtGui.QLabel()
			self.label_src.setText('Rotation Convention')
			self.hbl_src.addWidget(self.label_src)
			
			self.src = QtGui.QComboBox(maintab)
			self.load_src_options(self.src)
			self.hbl_src.addWidget(self.src)
			
			# set default value -1 ensures that the val slider is updated the first time it is created
			self.az = ValSlider(self,(-360.0,360.0),"az",-1)
			self.az.setObjectName("az")
			maintab.vbl.addWidget(self.az)
			
			self.alt = ValSlider(self,(-180.0,180.0),"alt",-1)
			self.alt.setObjectName("alt")
			maintab.vbl.addWidget(self.alt)
			
			self.phi = ValSlider(self,(-360.0,360.0),"phi",-1)
			self.phi.setObjectName("phi")
			maintab.vbl.addWidget(self.phi)
		
			self.current_src = EULER_EMAN
		
		return self.maintab

	def set_xy_trans(self, x, y):
		self.x_trans.setValue(x)
		self.y_trans.setValue(y)
	
	def set_xyz_trans(self,x,y,z):
		self.x_trans.setValue(x)
		self.y_trans.setValue(y)
		self.z_trans.setValue(z)
		
	def set_translate_scale(self, xscale,yscale,zscale):
		self.x_trans.setSingleStep(xscale)
		self.y_trans.setSingleStep(yscale)
		self.z_trans.setSingleStep(zscale)

	def update_rotations(self,t3d):
		convention = str( self.src.currentText() )
		#FIXME: Transform.get_rotation() wants a string sometimes and a EulerType other times
		try:
			rot = t3d.get_rotation(str(self.src_map[convention]))
		except StandardError, e: #doing a quick fix
			print e
			print "Developers: This catches a large range of exceptions... a better way surely exists"
			rot = t3d.get_rotation(self.src_map[convention])
		
		if ( self.src_map[convention] == EULER_SPIN ):
			self.n3.setValue(rot[self.n3.getLabel()],True)
		
		self.az.setValue(rot[self.az.getLabel()],True)
		self.alt.setValue(rot[self.alt.getLabel()],True)
		self.phi.setValue(rot[self.phi.getLabel()],True)
	
	def slider_rotate(self):
		self.target.load_rotation(self.get_current_rotation())
	
	def get_current_rotation(self):
		convention = self.src.currentText()
		rot = {}
		if ( self.current_src == EULER_SPIN ):
			rot[self.az.getLabel()] = self.az.getValue()
			
			n1 = self.alt.getValue()
			n2 = self.phi.getValue()
			n3 = self.n3.getValue()
			
			norm = sqrt(n1*n1 + n2*n2 + n3*n3)
			
			n1 /= norm
			n2 /= norm
			n3 /= norm
			
			rot[self.alt.getLabel()] = n1
			rot[self.phi.getLabel()] = n2
			rot[self.n3.getLabel()] = n3
			
		else:
			rot[self.az.getLabel()] = self.az.getValue()
			rot[self.alt.getLabel()] = self.alt.getValue()
			rot[self.phi.getLabel()] = self.phi.getValue()
		
		return Transform(self.current_src, rot)
	
	def set_src(self, val):
		t3d = self.get_current_rotation()
		
		if (self.n3_showing) :
			self.vbl.removeWidget(self.n3)
			self.n3.deleteLater()
			self.n3_showing = False
			self.az.setRange(-360,360)
			self.alt.setRange(-180,180)
			self.phi.setRange(-360,660)
		
		if ( self.src_map[str(val)] == EULER_SPIDER ):
			self.az.setLabel('phi')
			self.alt.setLabel('theta')
			self.phi.setLabel('psi')
		elif ( self.src_map[str(val)] == EULER_EMAN ):
			self.az.setLabel('az')
			self.alt.setLabel('alt')
			self.phi.setLabel('phi')
		elif ( self.src_map[str(val)] == EULER_IMAGIC ):
			self.az.setLabel('alpha')
			self.alt.setLabel('beta')
			self.phi.setLabel('gamma')
		elif ( self.src_map[str(val)] == EULER_XYZ ):
			self.az.setLabel('xtilt')
			self.alt.setLabel('ytilt')
			self.phi.setLabel('ztilt')
		elif ( self.src_map[str(val)] == EULER_MRC ):
			self.az.setLabel('phi')
			self.alt.setLabel('theta')
			self.phi.setLabel('omega')
		elif ( self.src_map[str(val)] == EULER_SPIN ):
			self.az.setLabel('omega')
			self.alt.setRange(-1,1)
			self.phi.setRange(-1,1)
			
			self.alt.setLabel('n1')
			self.phi.setLabel('n2')
			
			self.n3 = ValSlider(self,(-360.0,360.0),"n3",-1)
			self.n3.setRange(-1,1)
			self.n3.setObjectName("n3")
			self.vbl.addWidget(self.n3)
			QtCore.QObject.connect(self.n3, QtCore.SIGNAL("valueChanged"), self.slider_rotate)
			self.n3_showing = True
		
		self.current_src = self.src_map[str(val)]
		self.update_rotations(t3d)
	
	def load_src_options(self,widgit):
		self.load_src()
		for i in self.src_strings:
			widgit.addItem(i)
	
	# read src as 'supported rotation conventions'
	def load_src(self):
		# supported_rot_conventions
		src_flags = []
		src_flags.append(EULER_EMAN)
		src_flags.append(EULER_SPIDER)
		src_flags.append(EULER_IMAGIC)
		src_flags.append(EULER_MRC)
		src_flags.append(EULER_SPIN)
		src_flags.append(EULER_XYZ)
		
		self.src_strings = []
		self.src_map = {}
		for i in src_flags:
			self.src_strings.append(str(i))
			self.src_map[str(i)] = i
		
	
	def setColors(self,colors,current_color):
		a = 0
		for i in colors:
			self.cbb.addItem(i)
			if ( i == current_color):
				self.cbb.setCurrentIndex(a)
			a += 1

	def set_scale(self,newscale):
		self.scale.setValue(newscale)
		
# This is just for testing, of course
if __name__ == '__main__':
	from emapplication import EMApp
	from emimage3d import EMImage3DWidget
	em_app = EMApp()
	window = EMImage3DWidget()
	hello_world = EM3DHelloWorld(window)
	window.add_model(hello_world)
	em_app.show()
	em_app.execute()
	
