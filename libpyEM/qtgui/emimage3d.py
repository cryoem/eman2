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
from emimageutil import ImgHistogram,EMParentWin
from weakref import WeakKeyDictionary
from time import time
from PyQt4.QtCore import QTimer

from emimage3diso import EMIsosurface
from emimage3dvol import EMVolume
from emimage3dslice import EM3DSliceViewer
from emimage3dsym import EM3DSymViewer

from emglobjects import Camera2, EMViewportDepthTools, Camera

MAG_INCREMENT_FACTOR = 1.1

class EMImage3D(QtOpenGL.QGLWidget):
	""" 
	A QT widget for rendering 3D EMData objects
	"""
	allim=WeakKeyDictionary()
	def __init__(self, image=None, parent=None):
		self.image3d = None
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True)
		fmt.setDepth(True)
		fmt.setStencil(True)
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		EMImage3D.allim[self]=0
		
		self.aspect=1.0
		self.fov = 10 # field of view angle used by gluPerspective
		self.d = 0
		self.zwidth = 0
		self.perspective = True
		
		self.image3d = EMImage3DCore(image,self)
		self.initGL = True
		self.cam = Camera()
		
		if ( image != None and isinstance(image,EMData)):
			self.setCamZ(self.fov,image)
		
		self.resize(640,640)
		
		
		
	def setCamZ(self,fov,image):
		self.d = (image.get_ysize()/2.0)/tan(fov/2.0*pi/180.0)
		self.zwidth = image.get_zsize()
		self.ywidth = image.get_ysize()
		self.xwidth = image.get_xsize()
		self.cam.default_z = -self.d
		self.cam.cam_z = -self.d
	
	
	def setData(self,data):
		self.image3d.setData(data)
		if ( data != None and isinstance(data,EMData)):
			self.setCamZ(self.fov,data)
			
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
			self.image3d.initializeGL()
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
			self.image3d.initializeGL()
			self.initGL = False

		if ( self.image3d != None ):
			self.image3d.render()


	def resizeGL(self, width, height):
		# just use the whole window for rendering
		glViewport(0,0,self.width(),self.height())
		
		# maintain the aspect ratio of the window we have
		self.aspect = float(self.width())/float(self.height())
		
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		startz = self.d - 2.0*self.zwidth
		endz = self.d + 2.0*self.zwidth
		if self.perspective:
			# using gluPerspective for simplicity
			
			if startz < 0: startz = 1
			gluPerspective(self.fov,self.aspect,startz,endz)
		else:
			width = self.aspect*self.ywidth
			glOrtho(-width/2.0,width/2.0,-self.ywidth/2.0,self.ywidth/2.0,startz,endz)
			
		# switch back to model view mode
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		if (self.image3d != None):
			try: self.image3d.resizeEvent(width,height)
			except: pass
		
		self.updateGL()
		
	#def get_render_dims_at_depth(self, depth):
		## This function returns the width and height of the renderable 
		## area at the origin of the data volume
		#height = -2*tan(self.fov/2.0*pi/180.0)*(depth)
		#width = self.aspect*height
		#return [width,height]
			
	def mousePressEvent(self, event):
		self.image3d.mousePressEvent(event)
			
	def wheelEvent(self,event):
		self.image3d.wheelEvent(event)
	
	def mouseMoveEvent(self,event):
		self.image3d.mouseMoveEvent(event)

		
	def mouseReleaseEvent(self,event):
		self.image3d.mouseReleaseEvent(event)
		
	#def dropEvent(self,event):
		#self.image3d.dropEvent(event)
		
	def closeEvent(self,event) :
		self.image3d.closeEvent(event)
	
	def setPerspective(self,bool):
		self.perspective = bool
		self.resizeGL(self.width(),self.height())
	#def dragEnterEvent(self,event):
		#self.image3d.dragEnterEvent(event)

class EMImage3DCore:

	def __init__(self, image=None, parent=None):
		self.parent = parent
		
		self.currentselection = -1
		self.inspector = None
		#self.isosurface = EMIsosurface(image,self)
		#self.volume = EMVolume(image,self)
		self.viewables = []
		self.num_iso = 0
		self.num_vol = 0
		self.num_sli = 0
		self.num_sym = 0
		self.supressInspector = False 	# Suppresses showing the inspector - switched on in emfloatingwidgets
		
		#self.timer = QTimer()
		#QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.timeout)

		self.cam = Camera2(self)
		self.vdtools = EMViewportDepthTools(self)
		
		self.rottarget = None
		self.setData(image)
		#self.inspector.addIso()
	#def timeout(self):
		#self.updateGL()
	def width(self):
		try: return self.parent.width()
		except: return 0
		
	def height(self):
		try: return self.parent.height()
		except: return 0
	
	def updateGL(self):
		try: self.parent.updateGL()
		except: pass
	
	def eyeCoordsDif(self,x1,y1,x2,y2,mdepth=True):
		return self.vdtools.eyeCoordsDif(x1,y1,x2,y2,mdepth)

	def viewportHeight(self):
		return self.parent.height()
	
	def viewportWidth(self):
		return self.parent.width()
	
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
			
	def setData(self,data):
		if data == None: return
		self.image = data
		for i in self.viewables:
			i.setData(data)
			
		self.resizeEvent(self.parent.width(),self.parent.height())
		#self.volume.setData(data)
		
		if self.inspector == None:
			self.inspector=EMImageInspector3D(self)
		self.inspector.addIsosurface()
	
	def showInspector(self,force=0):
		if self.supressInspector: return
		if not force and self.inspector==None : return
		self.initInspector()
		self.inspector.show()
		
	def initInspector(self):
		if not self.inspector : self.inspector=EMImageInspector3D(self)
	
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
	
			
	def setCamZ(self,z):
		self.cam.setCamZ( z )
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
				self.rottarget.setScale(self.cam.scale)
		self.updateGL()
	
	def setScale(self,val):
		self.cam.scale = val
		self.updateGL()
	
	def mousePressEvent(self, event):
		if event.button()==Qt.MidButton:
			self.showInspector(1)
		else:
			self.cam.mousePressEvent(event)
				
		self.updateGL()
	
	def get_render_dims_at_depth(self, depth):
		return self.parent.get_render_dims_at_depth(depth)

	def getSundryInspector(self):
		return self.viewables[self.currentselection].getInspector()
	
	def addSym(self):
		sym = EM3DSymViewer(self)
		self.viewables.append(sym)
		#if (len(self.viewables)==1):
			#sym.cam.default_z = -1.25*self.image.get_zsize()
			#sym.cam.cam_z = -1.25*self.image.get_zsize()
		#else:
			#pass
			##self.loadLastViewableCamera()
		sym.setRadius(self.image.get_zsize()/2.0)
		self.num_sym += 1
		name = "Sym " + str(self.num_sym)
		self.viewables[len(self.viewables)-1].setName(name)
		self.viewables[len(self.viewables)-1].setRank(len(self.viewables))
		self.currentselection = len(self.viewables)-1
		self.updateGL()
	
	def addIsosurface(self):
		self.viewables.append(EMIsosurface(self.image,self))
		#self.loadLastViewableCamera()
		self.num_iso += 1
		name = "Isosurface " + str(self.num_iso)
		self.viewables[len(self.viewables)-1].setName(name)
		self.viewables[len(self.viewables)-1].setRank(len(self.viewables))
		self.currentselection = len(self.viewables)-1
		self.updateGL()
		
	def addVolume(self):
		self.viewables.append(EMVolume(self.image,self))
		#self.loadLastViewableCamera()
		self.num_vol += 1
		name = "Volume " + str(self.num_vol)
		self.viewables[len(self.viewables)-1].setName(name)
		self.viewables[len(self.viewables)-1].setRank(len(self.viewables))
		self.currentselection = len(self.viewables)-1
		self.updateGL()
		
	def addSliceViewer(self):
		self.viewables.append(EM3DSliceViewer(self.image,self))
		#self.loadLastViewableCamera()
		self.num_sli += 1
		name = "Slices " + str(self.num_sli)
		self.viewables[len(self.viewables)-1].setName(name)
		self.viewables[len(self.viewables)-1].setRank(len(self.viewables))
		self.currentselection = len(self.viewables)-1
		self.updateGL()
		
	def loadLastViewableCamera(self):
		return
		size = len(self.viewables)
		if ( size <= 1 ): return
		self.viewables[size-1].setCamera(self.viewables[0].getCurrentCamera())

	def rowChanged(self,row):
		if ( row == self.currentselection ): return
		self.currentselection=row
		self.updateGL()
		
	def getcuridx(self):
		return self.currentselection
		
	def getCurrentName(self):
		if self.currentselection == -1 : return ""
		elif self.currentselection >= len(self.viewables):
			print "error, current seletion too large", self.currentselection,len(self.viewables)
			return ""
		return self.viewables[self.currentselection].getName()
	
	def getCurrentInspector(self):
		if self.currentselection == -1 : return None
		elif self.currentselection >= len(self.viewables):
			print "error, current seletion too large", self.currentselection,len(self.viewables)
			return None
		return self.viewables[self.currentselection].getInspector()
	
	def deleteCurrent(self, val):
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
		
	def setPerspective(self,bool):
		self.parent.setPerspective(bool)
		
	def loadRotation(self,t3d):
		self.cam.loadRotation(t3d)
		self.updateGL()

	def getCurrentT3d(self):
		size = len(self.cam.t3d_stack)
		return self.cam.t3d_stack[size-1]
	
	def registerRotTarget(self, targ):
		self.rottarget = targ
	
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
		self.vbl.addWidget(self.tabwidget)
		
		self.hbl_check = QtGui.QHBoxLayout()
		self.hbl_check.setMargin(0)
		self.hbl_check.setSpacing(6)
		self.hbl_check.setObjectName("hbl_check")
		
		self.setcheck = QtGui.QCheckBox("Properties",self)
		self.hbl_check.addWidget(self.setcheck)
		
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
		
		self.addSym = QtGui.QPushButton("Sym")
		self.hbl_buttons2.addWidget(self.addSym)

		self.vbl.addLayout(self.hbl_buttons)
		self.vbl.addLayout(self.hbl_buttons2)
		
		self.hbl_buttons3 = QtGui.QHBoxLayout()
		self.delete = QtGui.QPushButton("Delete")
		self.hbl_buttons3.addWidget(self.delete)
		self.vbl.addLayout(self.hbl_buttons3)
		
		self.vbl.addLayout(self.hbl_check)
		
		self.setinspector = None
		
		self.currentselection = -1
		self.settingsrow = -2
		self.targetidxmap = {}
		
		QtCore.QObject.connect(self.addIso, QtCore.SIGNAL("clicked()"), self.addIsosurface)
		QtCore.QObject.connect(self.addVol, QtCore.SIGNAL("clicked()"), self.addVolume)
		QtCore.QObject.connect(self.addSli, QtCore.SIGNAL("clicked()"), self.addSlices)
		QtCore.QObject.connect(self.addSym, QtCore.SIGNAL("clicked()"), self.addSymmetry)
		QtCore.QObject.connect(self.delete, QtCore.SIGNAL("clicked()"), self.deleteSelection)
		QtCore.QObject.connect(self.setcheck, QtCore.SIGNAL("stateChanged(int)"), self.setChanged)
		
		#QtCore.QObject.connect(self.listwidget, QtCore.SIGNAL("currentRowChanged(int)"), self.rowChanged)
		#QtCore.QObject.connect(self.tabwidget, QtCore.SIGNAL("currentChanged(int)"), self.tabChanged)
		
		
	def setChanged(self, val):
		if val > 0:
			if self.setinspector == None:
				self.setinspector = EM3DSettingsInspector(self.target, self)
			
			self.target.registerRotTarget(self.setinspector)
			self.setinspector.updateRotations(self.target.getCurrentT3d())
			self.setinspector.setScale(self.target.cam.scale)
			self.tabwidget.addTab(self.setinspector,"Properties")
			self.settingsrow = self.tabwidget.count()-1
			self.targetidxmap[self.settingsrow] = -1
			self.tabwidget.setCurrentIndex(self.settingsrow)
		else:
			self.tabwidget.removeTab(self.settingsrow)
			self.target.registerRotTarget(None)
			# this isn't -1, because self.tabwidget.currentIndex() will possibly return -1 if
			# it has no tabs
			self.settingsrow = -2

	def addIsosurface(self):
		self.target.addIsosurface()
		self.updateSelection()
	
	def addSymmetry(self):
		self.target.addSym()
		self.updateSelection()
	
	def addVolume(self):
		self.target.addVolume()
		self.updateSelection()
	
	def updateSelection(self):
		self.tabwidget.addTab(self.target.getCurrentInspector(), self.target.getCurrentName())
		row = self.tabwidget.count()-1
		self.targetidxmap[row] = self.target.currentselection
		self.tabwidget.setCurrentIndex(row)

	def addSlices(self):
		self.target.addSliceViewer()
		self.updateSelection()
	
	def deleteSelection(self):
		idx = self.tabwidget.currentIndex()
		if ( idx == self.settingsrow ):
			self.setcheck.setCheckState(Qt.Unchecked)
			return
		self.tabwidget.removeTab(idx)
		self.target.deleteCurrent(self.targetidxmap[idx])

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
		QtCore.QObject.connect(self.scale, QtCore.SIGNAL("valueChanged"), target.setScale)
		QtCore.QObject.connect(self.x_trans, QtCore.SIGNAL("valueChanged(double)"), target.setCamX)
		QtCore.QObject.connect(self.y_trans, QtCore.SIGNAL("valueChanged(double)"), target.setCamY)
		QtCore.QObject.connect(self.z_trans, QtCore.SIGNAL("valueChanged(double)"), target.setCamZ)
		
	def sliderRotate(self):
		self.target.loadRotation(self.getCurrentRotation())
		
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
		
	def setScale(self,newscale):
		self.scale.setValue(newscale)
		
	def setXYTrans(self, x, y):
		self.x_trans.setValue(x)
		self.y_trans.setValue(y)
	

class EM3DSettingsInspector(QtGui.QWidget):
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
		self.hbl.addLayout(self.groupbox)
		
		self.vbl.addLayout(self.hbl)
		
		self.rotsliders.addWidgets(self.vbl)
		
		QtCore.QObject.connect(self.persbut, QtCore.SIGNAL("pressed()"), self.persClicked)
		QtCore.QObject.connect(self.orthbut, QtCore.SIGNAL("pressed()"), self.orthClicked)
		
	def updateRotations(self,t3d):
		self.rotsliders.updateRotations(t3d)
		
	def setScale(self,val):
		self.rotsliders.setScale(val)
	
	def setXYTrans(self, x, y):
		self.rotsliders.setXYTrans(x,y)
		
	def persClicked(self):
		self.target.setPerspective(True)
		
	def orthClicked(self):
		self.target.setPerspective(False)
		
	
	
# This is just for testing, of course
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = EMImage3D()
 	if len(sys.argv)==1 : 
		pass
		e = EMData()
		e.set_size(64,64,64)
		e.process_inplace('testimage.axes')
 		window.setData(e)

		# these lines are for testing shape rendering
# 		window.addShape("a",["rect",.2,.8,.2,20,20,80,80,2])
# 		window.addShape("b",["circle",.5,.8,.2,120,50,30.0,2])
# 		window.addShape("c",["line",.2,.8,.5,20,120,100,200,2])
# 		window.addShape("d",["label",.2,.8,.5,220,220,"Testing",14,1])
	else :
		if not os.path.exists(sys.argv[1]):
			print "Error, input file %s does not exist" %sys.argv[1]
			exit(1)
		a=EMData.read_images(sys.argv[1],[0])
		window.setData(a[0])
	window2=EMParentWin(window)
	window2.show()
	
#	w2=QtGui.QWidget()
#	w2.resize(256,128)
	
#	w3=ValSlider(w2)
#	w3.resize(256,24)
#	w2.show()
	
	sys.exit(app.exec_())
