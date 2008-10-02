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

from weakref import WeakKeyDictionary
from time import time
from PyQt4.QtCore import QTimer

from emimage3diso import EMIsosurface
from emimage3dvol import EMVolume
from emimage3dslice import EM3DSliceViewer
from emimage3dsym import EM3DSymViewer

from emglobjects import Camera2, EMViewportDepthTools, Camera
from emimageutil import EventRerouter
from emapplication import EMStandAloneApplication, EMQtWidgetModule, EMModule

MAG_INCREMENT_FACTOR = 1.1

class EMImage3DWidget(QtOpenGL.QGLWidget,EventRerouter):
	""" 
	A QT widget for rendering 3D EMData objects
	"""
	allim=WeakKeyDictionary()
	def __init__(self, image_3d_module, parent=None):
		self.target = None
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True)
		fmt.setDepth(True)
		fmt.setStencil(True)
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		EMImage3DWidget.allim[self]=0
		EventRerouter.__init__(self,image_3d_module)
		
		self.aspect=1.0
		self.fov = 10 # field of view angle used by gluPerspective
		self.d = 0
		self.zwidth = 0
		self.perspective = True
		
		self.target = image_3d_module
		self.initGL = True
		self.cam = Camera()

		self.resize(640,640)
		self.startz = 0
		self.endz = 0
	
	def get_fov(self):
		return self.fov

	def set_cam_z(self,fov,image):
		self.d = (image.get_ysize()/2.0)/tan(fov/2.0*pi/180.0)
		self.zwidth = image.get_zsize()
		self.yheight = image.get_ysize()
		self.xwidth = image.get_xsize()
		self.cam.default_z = -self.d
		self.cam.cam_z = -self.d
	
	def set_data(self,data):
		self.target.set_data(data)
		if ( data != None and isinstance(data,EMData)):
			self.set_cam_z(self.fov,data)
			
		self.resize(640,640)
		
	def initializeGL(self):
		glEnable(GL_LIGHTING)
		glEnable(GL_LIGHT0)
		glEnable(GL_DEPTH_TEST)
		glLightfv(GL_LIGHT0, GL_AMBIENT, [0.9, 0.9, 0.9, 1.0])
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_POSITION, [0.5,0.7,11.,0.])
		glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST)
		glClearStencil(0)
		glEnable(GL_STENCIL_TEST)
		glClearColor(0,0,0,0)
		try:
			self.target.initializeGL()
			self.initGL = False
		except:
			pass
		
		
	def paintGL(self):
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT )
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		try:
			self.cam.position()
		except:
			return
		
		
		if ( self.initGL ):
			self.target.initializeGL()
			self.initGL = False

		if ( self.target != None ):
			self.target.render()


	def resizeGL(self, width, height):
		# just use the whole window for rendering
		glViewport(0,0,self.width(),self.height())
		
		# maintain the aspect ratio of the window we have
		self.aspect = float(self.width())/float(self.height())
		
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		
		if (self.zwidth == 0):
			# We've received  a resize event but no data has been set
			# in which case nothing is being rendered.
			# Therefore just leave the identity as the projection matrix.
			# This is an exceptional circumstance which probably 
			# highlights the need for some redesigning (d.woolford)
			glMatrixMode(GL_MODELVIEW)
			#glLoadIdentity()
			return
		
		self.startz = self.d - 2.0*self.zwidth
		self.endz = self.d + 2.0*self.zwidth
		if self.perspective:
			# using gluPerspective for simplicity
			
			if self.startz < 0: self.startz = 1
			gluPerspective(self.fov,self.aspect,self.startz,self.endz)
		else:
			self.xwidth = self.aspect*self.yheight
			glOrtho(-self.xwidth/2.0,self.xwidth/2.0,-self.yheight/2.0,self.yheight/2.0,self.startz,self.endz)
			
		# switch back to model view mode
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		if (self.target != None):
			try: self.target.resizeEvent(width,height)
			except: pass
		
		self.updateGL()
		
	#def get_render_dims_at_depth(self, depth):
		## This function returns the width and height of the renderable 
		## area at the origin of the data volume
		#height = -2*tan(self.fov/2.0*pi/180.0)*(depth)
		#width = self.aspect*height
		#return [width,height]
			
	
	def set_perspective(self,bool):
		self.perspective = bool
		self.resizeGL(self.width(),self.height())
		
		
	def get_start_z(self):
		return self.startz
	
	def get_near_plane_dims(self):
		if self.perspective:
			height = 2.0*self.startz * tan(self.fov/2.0*pi/180.0)
			width = self.aspect * height
			return [width,height]
		else:
			return [self.xwidth,self.yheight]

class EMImage3DModule(EMModule):
	
	def get_qt_widget(self):
		if self.parent == None:	
			self.parent = EMImage3DWidget(self)
			if isinstance(self.image,EMData):
				self.parent.set_cam_z(self.parent.get_fov(),self.image)
		return self.parent
	
	def __init__(self, image=None,application=None):
		self.parent = None
		self.image = None
		self.currentselection = -1
		self.inspector = None
		#self.isosurface = EMIsosurface(image,self)
		#self.volume = EMVolume(image,self)
		self.viewables = []
		self.num_iso = 0
		self.num_vol = 0
		self.num_sli = 0
		self.num_sym = 0
		self.suppress_inspector = False 	
		self.cam = Camera2(self)
		self.vdtools = EMViewportDepthTools(self)
		
		self.rottarget = None
		self.set_data(image)
			
		if application != None:
			application.attach_child(self)
			EMModule.set_app(self,application)
		self.em_qt_inspector_widget = None
			
	def width(self):
		try: return self.parent.width()
		except: return 0
		
	def height(self):
		try: return self.parent.height()
		except: return 0
	
	def updateGL(self):
		try: self.parent.updateGL()
		except: pass
	
	def eye_coords_dif(self,x1,y1,x2,y2,mdepth=True):
		return self.vdtools.eye_coords_dif(x1,y1,x2,y2,mdepth)
	
	def initializeGL(self):
		glEnable(GL_NORMALIZE)
	
	def render(self):
		glPushMatrix()
		self.cam.position(True)
		# the ones are dummy variables atm... they don't do anything
		self.vdtools.update(1,1)
		glPopMatrix()
		
		self.cam.position()
		
		for i in self.viewables:
			glPushMatrix()
			i.render()
			glPopMatrix()

	def resizeEvent(self, width, height):
		for i in self.viewables:
			i.resizeEvent()
	
	def get_data_dims(self):
		return [self.image.get_xsize(),self.image.get_ysize(),self.image.get_zsize()]

	def set_data(self,data):
		if data == None: return
		self.image = data
		for i in self.viewables:
			i.set_data(data)
		
		if self.parent != None: 
			self.resizeEvent(self.parent.width(),self.parent.height())
			self.parent.set_cam_z(self.parent.get_fov(),self.data)
		
		if self.inspector == None:
			self.inspector=EMImageInspector3D(self)
		self.inspector.add_isosurface()
	
	def show_inspector(self,force=0):
		if self.suppress_inspector: return
		if not force and self.inspector==None : return
		self.init_inspector()
		self.application.show_specific(self.em_qt_inspector_widget)
		
	def init_inspector(self):
		if not self.inspector : 
			self.inspector=EMImageInspector3D(self)
		if not self.em_qt_inspector_widget:
			self.em_qt_inspector_widget = EMQtWidgetModule(self.inspector,self.application)
			
	def closeEvent(self,event) :
		#for i in self.viewables:
			#i.closeEvent(event)
		if self.inspector: self.inspector.close()
		
	def mouseMoveEvent(self, event):
		self.cam.mouseMoveEvent(event)
		if self.rottarget != None :
			if event.buttons()&Qt.LeftButton:
				self.rottarget.updateRotations(self.getCurrentT3d())
			elif event.buttons()&Qt.RightButton:
				self.rottarget.setXYTrans(self.cam.cam_x, self.cam.cam_y)
		self.updateGL()
	
			
	def set_cam_z(self,z):
		self.cam.set_cam_z( z )
		self.updateGL()
		
	def setCamY(self,y):
		self.cam.setCamY( y )
		self.updateGL()
		
	def setCamX(self,x):
		self.cam.setCamX( x )
		self.updateGL()
	
	def mouseReleaseEvent(self, event):
		self.cam.mouseReleaseEvent(event)
		self.updateGL()
			
	def wheelEvent(self, event):
		self.cam.wheelEvent(event)
		if self.rottarget != None :
				self.rottarget.set_scale(self.cam.scale)
		self.updateGL()
	
	def set_scale(self,val):
		self.cam.scale = val
		self.updateGL()
	
	def mousePressEvent(self, event):
		if event.button()==Qt.MidButton:
			self.show_inspector(1)
		else:
			self.cam.mousePressEvent(event)
				
		self.updateGL()
	
	def get_render_dims_at_depth(self, depth):
		return self.parent.get_render_dims_at_depth(depth)

	def getSundryInspector(self):
		return self.viewables[self.currentselection].get_inspector()
	
	def add_sym(self):
		sym = EM3DSymViewer(self)
		self.viewables.append(sym)
		#if (len(self.viewables)==1):
			#sym.cam.default_z = -1.25*self.image.get_zsize()
			#sym.cam.cam_z = -1.25*self.image.get_zsize()
		#else:
			#pass
			##self.load_last_viewable_camera()
		sym.setRadius(self.image.get_zsize()/2.0)
		self.num_sym += 1
		name = "Sym " + str(self.num_sym)
		self.viewables[len(self.viewables)-1].setName(name)
		self.viewables[len(self.viewables)-1].setRank(len(self.viewables))
		self.currentselection = len(self.viewables)-1
		self.updateGL()
	
	def add_isosurface(self):
		self.viewables.append(EMIsosurface(self.image,self))
		#self.load_last_viewable_camera()
		self.num_iso += 1
		name = "Isosurface " + str(self.num_iso)
		self.viewables[len(self.viewables)-1].setName(name)
		self.viewables[len(self.viewables)-1].setRank(len(self.viewables))
		self.currentselection = len(self.viewables)-1
		self.updateGL()
		
	def add_volume(self):
		self.viewables.append(EMVolume(self.image,self))
		#self.load_last_viewable_camera()
		self.num_vol += 1
		name = "Volume " + str(self.num_vol)
		self.viewables[len(self.viewables)-1].setName(name)
		self.viewables[len(self.viewables)-1].setRank(len(self.viewables))
		self.currentselection = len(self.viewables)-1
		self.updateGL()
		
	def add_slice_viewer(self):
		self.viewables.append(EM3DSliceViewer(self.image,self))
		#self.load_last_viewable_camera()
		self.num_sli += 1
		name = "Slices " + str(self.num_sli)
		self.viewables[len(self.viewables)-1].setName(name)
		self.viewables[len(self.viewables)-1].setRank(len(self.viewables))
		self.currentselection = len(self.viewables)-1
		self.updateGL()
		
	def load_last_viewable_camera(self):
		return
		size = len(self.viewables)
		if ( size <= 1 ): return
		self.viewables[size-1].setCamera(self.viewables[0].getCurrentCamera())

	def rowChanged(self,row):
		if ( row == self.currentselection ): return
		self.currentselection=row
		self.updateGL()
		
	def get_current_idx(self):
		return self.currentselection
		
	def get_current_name(self):
		if self.currentselection == -1 : return ""
		elif self.currentselection >= len(self.viewables):
			print "error, current seletion too large", self.currentselection,len(self.viewables)
			return ""
		return self.viewables[self.currentselection].getName()
	
	def get_current_inspector(self):
		if self.currentselection == -1 : return None
		elif self.currentselection >= len(self.viewables):
			print "error, current seletion too large", self.currentselection,len(self.viewables)
			return None
		return self.viewables[self.currentselection].get_inspector()
	
	def delete_current(self, val):
		if ( len(self.viewables) == 0 ): return
		
		self.viewables.pop(val)
		if (len(self.viewables) == 0 ) : 
			self.currentselection = -1
		elif ( len(self.viewables) == 1):
			self.currentselection = 0
		elif ( val == 0):
			pass
		else:
			self.currentselection = val - 1
		
		
		# Need to set the rank appropriately
		for i in range(0,len(self.viewables)):
			self.viewables[i].setRank(i+1)
		
		self.updateGL()
	
	def resizeEvent(self,width=0,height=0):
		self.vdtools.set_update_P_inv()
		
	def set_perspective(self,bool):
		self.parent.set_perspective(bool)
		
	def load_rotation(self,t3d):
		self.cam.load_rotation(t3d)
		self.updateGL()

	def getCurrentT3d(self):
		size = len(self.cam.t3d_stack)
		return self.cam.t3d_stack[size-1]
	
	def registerRotTarget(self, targ):
		self.rottarget = targ
	
	def get_start_z(self):
		return self.parent.get_start_z()
	
	def get_near_plane_dims(self):
		return self.parent.get_near_plane_dims()
	
class EMImageInspector3D(QtGui.QWidget):
	def __init__(self,target) :
		QtGui.QWidget.__init__(self,None)
		self.target=target
		
		
		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		
		self.hbl = QtGui.QHBoxLayout()
		self.hbl.setMargin(2)
		self.hbl.setSpacing(6)
		self.hbl.setObjectName("hbl")
		
		#self.listwidget = QtGui.QListWidget(self)
		#self.vbl.addWidget(self.listwidget)
		
		self.tabwidget = QtGui.QTabWidget(self)
		
		self.hbl_check = QtGui.QHBoxLayout()
		self.hbl_check.setMargin(0)
		self.hbl_check.setSpacing(6)
		self.hbl_check.setObjectName("hbl_check")
		
		#self.advancedcheck = QtGui.QCheckBox("Advanced",self)
		#self.hbl_check.addWidget(self.advancedcheck)
		
		self.hbl_buttons = QtGui.QHBoxLayout()
		self.hbl_buttons.setMargin(0)
		self.hbl_buttons.setSpacing(6)
		self.hbl_buttons.setObjectName("hbl_buttons")
		
		self.hbl_buttons2 = QtGui.QHBoxLayout()
		self.hbl_buttons2.setMargin(0)
		self.hbl_buttons2.setSpacing(6)
		self.hbl_buttons2.setObjectName("hbl_buttons2")
		
		self.addIso = QtGui.QPushButton("Isosurface")
		self.hbl_buttons.addWidget(self.addIso)
		
		self.addVol = QtGui.QPushButton("Volume")
		self.hbl_buttons.addWidget(self.addVol)
		
		self.addSli = QtGui.QPushButton("Slices")
		self.hbl_buttons2.addWidget(self.addSli)
		
		self.add_sym = QtGui.QPushButton("Sym")
		self.hbl_buttons2.addWidget(self.add_sym)

		self.vbl.addLayout(self.hbl_buttons)
		self.vbl.addLayout(self.hbl_buttons2)
		
		self.hbl_buttons3 = QtGui.QHBoxLayout()
		self.delete = QtGui.QPushButton("Delete")
		self.hbl_buttons3.addWidget(self.delete)
		self.vbl.addLayout(self.hbl_buttons3)
		
		self.vbl.addLayout(self.hbl_check)
		self.vbl.addWidget(self.tabwidget)
		
		self.setinspector = None
		
		self.currentselection = -1
		self.settingsrow = -2
		self.targetidxmap = {}
		
		#self.advancedcheck.click()
		self.insertAdvancedTab()
		
		QtCore.QObject.connect(self.addIso, QtCore.SIGNAL("clicked()"), self.add_isosurface)
		QtCore.QObject.connect(self.addVol, QtCore.SIGNAL("clicked()"), self.add_volume)
		QtCore.QObject.connect(self.addSli, QtCore.SIGNAL("clicked()"), self.addSlices)
		QtCore.QObject.connect(self.add_sym, QtCore.SIGNAL("clicked()"), self.add_symmetry)
		QtCore.QObject.connect(self.delete, QtCore.SIGNAL("clicked()"), self.deleteSelection)
		#QtCore.QObject.connect(self.advancedcheck, QtCore.SIGNAL("stateChanged(int)"), self.advancedClicked)
		
		#QtCore.QObject.connect(self.listwidget, QtCore.SIGNAL("currentRowChanged(int)"), self.rowChanged)
		#QtCore.QObject.connect(self.tabwidget, QtCore.SIGNAL("currentChanged(int)"), self.tabChanged)
		
		
	def insertAdvancedTab(self):
		if self.setinspector == None:
			self.setinspector = EM3DAdvancedInspector(self.target, self)
			
		self.target.registerRotTarget(self.setinspector)
		self.setinspector.updateRotations(self.target.getCurrentT3d())
		self.setinspector.set_scale(self.target.cam.scale)
		self.tabwidget.addTab(self.setinspector,"Advanced")
		self.settingsrow = self.tabwidget.count()-1
		self.targetidxmap[self.settingsrow] = -1
		self.tabwidget.setCurrentIndex(self.settingsrow)
	

	def add_isosurface(self):
		self.target.add_isosurface()
		self.updateSelection()
	
	def add_symmetry(self):
		self.target.add_sym()
		self.updateSelection()
	
	def add_volume(self):
		self.target.add_volume()
		self.updateSelection()
	
	def updateSelection(self):
		n = self.tabwidget.count()
		if n > 0: n = n - 1
		self.tabwidget.insertTab(n, self.target.get_current_inspector(), self.target.get_current_name())
		self.targetidxmap[n] = self.target.currentselection
		self.tabwidget.setCurrentIndex(n)

	def addSlices(self):
		self.target.add_slice_viewer()
		self.updateSelection()
	
	def deleteSelection(self):
		idx = self.tabwidget.currentIndex()
		n = self.tabwidget.count()
		if n <= 1: return
		if idx == n-1: return
		
		self.tabwidget.removeTab(idx)
		self.target.delete_current(self.targetidxmap[idx])

class EMRotateSliders:
	def __init__(self,target,parent):
		self.target = target
		self.parent = parent
		
		self.label_src = QtGui.QLabel(parent)
		self.label_src.setText('Rotation Convention')
		
		self.src = QtGui.QComboBox(parent)
		self.load_src_options(self.src)
		
		self.x_label = QtGui.QLabel()
		self.x_label.setText('x')
		
		self.x_trans = QtGui.QDoubleSpinBox(parent)
		self.x_trans.setMinimum(-10000)
		self.x_trans.setMaximum(10000)
		self.x_trans.setValue(0.0)
	
		self.y_label = QtGui.QLabel()
		self.y_label.setText('y')
		
		self.y_trans = QtGui.QDoubleSpinBox(parent)
		self.y_trans.setMinimum(-10000)
		self.y_trans.setMaximum(10000)
		self.y_trans.setValue(0.0)
		
		self.z_label = QtGui.QLabel()
		self.z_label.setText('z')
		
		self.z_trans = QtGui.QDoubleSpinBox(parent)
		self.z_trans.setMinimum(-10000)
		self.z_trans.setMaximum(10000)
		self.z_trans.setValue(0.0)
		
		self.az = ValSlider(parent,(-360.0,360.0),"az",-1)
		self.az.setObjectName("az")
		self.az.setValue(0.0)
		
		self.alt = ValSlider(parent,(-180.0,180.0),"alt",-1)
		self.alt.setObjectName("alt")
		self.alt.setValue(0.0)
		
		self.phi = ValSlider(parent,(-360.0,360.0),"phi",-1)
		self.phi.setObjectName("phi")
		self.phi.setValue(0.0)
		
		self.scale = ValSlider(parent,(0.01,30.0),"Zoom:")
		self.scale.setObjectName("scale")
		self.scale.setValue(1.0)
		
		self.n3_showing = False
		
		self.current_src = EULER_EMAN
		
		QtCore.QObject.connect(self.az, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
		QtCore.QObject.connect(self.alt, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
		QtCore.QObject.connect(self.phi, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
		QtCore.QObject.connect(self.src, QtCore.SIGNAL("currentIndexChanged(QString)"), self.set_src)
		QtCore.QObject.connect(self.scale, QtCore.SIGNAL("valueChanged"), target.set_scale)
		QtCore.QObject.connect(self.x_trans, QtCore.SIGNAL("valueChanged(double)"), target.setCamX)
		QtCore.QObject.connect(self.y_trans, QtCore.SIGNAL("valueChanged(double)"), target.setCamY)
		QtCore.QObject.connect(self.z_trans, QtCore.SIGNAL("valueChanged(double)"), target.set_cam_z)
		
	def sliderRotate(self):
		self.target.load_rotation(self.getCurrentRotation())
		
	def getCurrentRotation(self):
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
		
		return Transform3D(self.current_src, rot)
	
	def addWidgets(self,target):
		self.hbl_trans = QtGui.QHBoxLayout()
		self.hbl_trans.setMargin(0)
		self.hbl_trans.setSpacing(6)
		self.hbl_trans.setObjectName("Trans")
		self.hbl_trans.addWidget(self.x_label)
		self.hbl_trans.addWidget(self.x_trans)
		self.hbl_trans.addWidget(self.y_label)
		self.hbl_trans.addWidget(self.y_trans)
		self.hbl_trans.addWidget(self.z_label)
		self.hbl_trans.addWidget(self.z_trans)
		
		target.addLayout(self.hbl_trans)
		
		self.hbl_src = QtGui.QHBoxLayout()
		self.hbl_src.setMargin(0)
		self.hbl_src.setSpacing(6)
		self.hbl_src.setObjectName("hbl")
		self.hbl_src.addWidget(self.label_src)
		self.hbl_src.addWidget(self.src)
		
		target.addWidget(self.scale)
		target.addLayout(self.hbl_src)
		target.addWidget(self.az)
		target.addWidget(self.alt)
		target.addWidget(self.phi)
	
	def set_src(self, val):
		t3d = self.getCurrentRotation()
		
		if (self.n3_showing) :
			self.maintab.vbl.removeWidget(self.n3)
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
			self.az.setLabel('Omega')
			self.alt.setRange(-1,1)
			self.phi.setRange(-1,1)
			
			self.alt.setLabel('n1')
			self.phi.setLabel('n2')
			
			self.n3 = ValSlider(self,(-360.0,360.0),"n3",-1)
			self.n3.setRange(-1,1)
			self.n3.setObjectName("n3")
			self.maintab.vbl.addWidget(self.n3)
			QtCore.QObject.connect(self.n3, QtCore.SIGNAL("valueChanged"), self.sliderRotate)
			self.n3_showing = True
		
		self.current_src = self.src_map[str(val)]
		self.updateRotations(t3d)
	
	def load_src_options(self,widgit):
		self.load_src()
		for i in self.src_strings:
			widgit.addItem(i)
			
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
			
	def updateRotations(self,t3d):
		rot = t3d.get_rotation(self.src_map[str(self.src.itemText(self.src.currentIndex()))])
		
		convention = self.src.currentText()
		if ( self.src_map[str(convention)] == EULER_SPIN ):
			self.n3.setValue(rot[self.n3.getLabel()],True)
		
		self.az.setValue(rot[self.az.getLabel()],True)
		self.alt.setValue(rot[self.alt.getLabel()],True)
		self.phi.setValue(rot[self.phi.getLabel()],True)
		
	def set_scale(self,newscale):
		self.scale.setValue(newscale)
		
	def setXYTrans(self, x, y):
		self.x_trans.setValue(x)
		self.y_trans.setValue(y)
	

class EM3DAdvancedInspector(QtGui.QWidget):
	def __init__(self,target,parent=None):
		QtGui.QWidget.__init__(self,None)
		self.target=target
		self.parent=parent
		
		self.rotsliders = EMRotateSliders(target,self)
		
		self.hbl = QtGui.QHBoxLayout()
		self.hbl.setMargin(2)
		self.hbl.setSpacing(6)
		self.hbl.setObjectName("hbl")
		
		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		
		self.persbut = QtGui.QRadioButton("Perspective")
		self.persbut.setChecked(True)
		self.orthbut = QtGui.QRadioButton("Orthographic")
		
		self.groupbox = QtGui.QVBoxLayout()
		self.groupbox.addWidget(self.persbut)
		self.groupbox.addWidget(self.orthbut)
		
		self.viewingvol = QtGui.QGroupBox("Viewing Volume")
		self.viewingvol.setLayout(self.groupbox)
		
		self.hbl.addWidget(self.viewingvol)
		
		self.vbl.addLayout(self.hbl)
		
		self.rotsliders.addWidgets(self.vbl)
		
		QtCore.QObject.connect(self.persbut, QtCore.SIGNAL("pressed()"), self.persClicked)
		QtCore.QObject.connect(self.orthbut, QtCore.SIGNAL("pressed()"), self.orthClicked)
		
	def updateRotations(self,t3d):
		self.rotsliders.updateRotations(t3d)
		
	def set_scale(self,val):
		self.rotsliders.set_scale(val)
	
	def setXYTrans(self, x, y):
		self.rotsliders.setXYTrans(x,y)
		
	def persClicked(self):
		self.target.set_perspective(True)
		
	def orthClicked(self):
		self.target.set_perspective(False)
		
	
	
if __name__ == '__main__':
	em_app = EMStandAloneApplication()
	window = EMImage3DModule(application=em_app)
	
	if len(sys.argv)==1 : 
		data = []
		#for i in range(0,200):
		e = EMData(64,64,64)
		e.process_inplace('testimage.axes')
		window.set_data(e)
	else :
		a=EMData.read_images(sys.argv[1])
		window.set_file_name(sys.argv[1])
		window.set_data(a)
		
	em_app.show()
	em_app.execute()

# This is just for testing, of course
#if __name__ == '__main__':
	#app = QtGui.QApplication(sys.argv)
	#window = EMImage3D()
 	#if len(sys.argv)==1 : 
		#pass
		#e = EMData()
		#e.set_size(64,64,64)
		#e.process_inplace('testimage.axes')
 		#window.set_data(e)
	#else :
		#if not os.path.exists(sys.argv[1]):
			#print "Error, input file %s does not exist" %sys.argv[1]
			#exit(1)
		#a=EMData.read_images(sys.argv[1],[0])
		#window.set_data(a[0])
	#window2=EMParentWin(window)
	#window2.show()

	#sys.exit(app.exec_())
