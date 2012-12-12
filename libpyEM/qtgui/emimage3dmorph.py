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

#from PyQt4 import QtCore, QtGui, QtOpenGL
#from PyQt4.QtCore import Qt
#from OpenGL import GL,GLU,GLUT
#from OpenGL.GL import *
#from OpenGL.GLU import *
#from valslider import ValSlider
#from math import *
#from EMAN2 import *
#import sys
#import numpy
#from emimageutil import ImgHistogram,EMParentWin
#from weakref import WeakKeyDictionary
#from time import time
#from PyQt4.QtCore import QTimer

#from emimage3diso import EMIsosurface
#from emimage3dvol import EMVolume
#from emimage3dslice import EM3DSliceViewerModule
#from emimage3dsym import EM3DSymViewer

#from emglobjects import Camera2, EMViewportDepthTools, Camera

#MAG_INCREMENT_FACTOR = 1.1

#class EMImageMorph3D(QtOpenGL.QGLWidget):
	#""" 
	#A QT widget for rendering 3D EMData objects
	#"""
	#allim=WeakKeyDictionary()
	#def __init__(self, image1=None,image2=None, parent=None):
		#self.image3d = None
		#fmt=QtOpenGL.QGLFormat()
		#fmt.setDoubleBuffer(True)
		#fmt.setDepth(True)
		#fmt.setStencil(True)
		#QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		#EMImageMorph3D.allim[self]=0
		
		#self.aspect=1.0
		#self.fov = 10 # field of view angle used by gluPerspective
		#self.d = 0
		#self.zwidth = 0
		#self.perspective = True
		
		#self.image3d = EMImageMorph3DCore(image1,image2,self)
		#self.initGL = True
		#self.cam = Camera()
		
		#if ( image1 != None and isinstance(image,EMData)):
			#self.set_cam_z(self.fov,image1)
		
		#self.resize(640,640)
		#self.startz = 0
		#self.endz = 0
		
		
	#def set_cam_z(self,fov,image):
		#self.d = (image.get_ysize()/2.0)/tan(fov/2.0*pi/180.0)
		#self.zwidth = image.get_zsize()
		#self.yheight = image.get_ysize()
		#self.xwidth = image.get_xsize()
		#self.cam.default_z = -self.d
		#self.cam.cam_z = -self.d
	
	#def set_data(self,dataA,dataB):
		#self.image3d.set_data(dataA,dataB)
		#if ( dataA != None and dataB != None and isinstance(dataA,EMData) and isinstance(dataB,EMData)):
			#self.set_cam_z(self.fov,dataA)
			
		#self.resize(640,640)
		
	#def initializeGL(self):
		#glEnable(GL_LIGHTING)
		#glEnable(GL_LIGHT0)
		#glEnable(GL_DEPTH_TEST)
		#glLightfv(GL_LIGHT0, GL_AMBIENT, [0.9, 0.9, 0.9, 1.0])
		#glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
		#glLightfv(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
		#glLightfv(GL_LIGHT0, GL_POSITION, [0.5,0.7,11.,0.])
		#glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST)
		#glClearStencil(0)
		#glEnable(GL_STENCIL_TEST)
		#glClearColor(0,0,0,0)
		#try:
			#self.image3d.initializeGL()
			#self.initGL = False
		#except:
			#pass
		
		
	#def paintGL(self):
		#glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT )
		#glMatrixMode(GL_MODELVIEW)
		#glLoadIdentity()
		#try:
			#self.cam.position()
		#except:
			#return
		
		
		#if ( self.initGL ):
			#self.image3d.initializeGL()
			#self.initGL = False

		#if ( self.image3d != None ):
			#self.image3d.render()


	#def resizeGL(self, width, height):
		## just use the whole window for rendering
		#glViewport(0,0,self.width(),self.height())
		
		## maintain the aspect ratio of the window we have
		#self.aspect = float(self.width())/float(self.height())
		
		#glMatrixMode(GL_PROJECTION)
		#glLoadIdentity()
		
		#if (self.zwidth == 0):
			## We've received  a resize event but no data has been set
			## in which case nothing is being rendered.
			## Therefore just leave the identity as the projection matrix.
			## This is an exceptional circumstance which probably 
			## highlights the need for some redesigning (d.woolford)
			#glMatrixMode(GL_MODELVIEW)
			##glLoadIdentity()
			#return
		
		#self.startz = self.d - 2.0*self.zwidth
		#self.endz = self.d + 2.0*self.zwidth
		#if self.perspective:
			## using gluPerspective for simplicity
			
			#if self.startz < 0: self.startz = 1
			#gluPerspective(self.fov,self.aspect,self.startz,self.endz)
		#else:
			#self.xwidth = self.aspect*self.yheight
			#glOrtho(-self.xwidth/2.0,self.xwidth/2.0,-self.yheight/2.0,self.yheight/2.0,self.startz,self.endz)
			
		## switch back to model view mode
		#glMatrixMode(GL_MODELVIEW)
		#glLoadIdentity()
		
		#if (self.image3d != None):
			#try: self.image3d.resizeEvent(width,height)
			#except: pass
		
		#self.updateGL()
		
	##def get_render_dims_at_depth(self, depth):
		### This function returns the width and height of the renderable 
		### area at the origin of the data volume
		##height = -2*tan(self.fov/2.0*pi/180.0)*(depth)
		##width = self.aspect*height
		##return [width,height]
			
	#def mousePressEvent(self, event):
		#self.image3d.mousePressEvent(event)
			
	#def wheelEvent(self,event):
		#self.image3d.wheelEvent(event)
	
	#def mouseMoveEvent(self,event):
		#self.image3d.mouseMoveEvent(event)

		
	#def mouseReleaseEvent(self,event):
		#self.image3d.mouseReleaseEvent(event)
		
	##def dropEvent(self,event):
		##self.image3d.dropEvent(event)
		
	#def closeEvent(self,event) :
		#self.image3d.closeEvent(event)
	
	#def set_perspective(self,bool):
		#self.perspective = bool
		#self.resizeGL(self.width(),self.height())
		
		
	#def get_start_z(self):
		#return self.startz
	
	#def get_near_plane_dims(self):
		#if self.perspective:
			#height = 2.0*self.startz * tan(self.fov/2.0*pi/180.0)
			#width = self.aspect * height
			#return [width,height]
		#else:
			#return [self.xwidth,self.yheight]
		
	##def dragEnterEvent(self,event):
		##self.image3d.dragEnterEvent(event)

#class EMImageMorph3DCore:

	#def __init__(self, image1=None,image2=None, parent=None):
		#self.parent = parent
		
		#self.currentselection = -1
		#self.inspector = None
		##self.isosurface = EMIsosurface(image,self)
		##self.volume = EMVolumeModule(image,self)
		#self.viewables = []
		#self.num_iso = 0
		#self.num_vol = 0
		#self.num_sli = 0
		#self.num_sym = 0
		#self.suppress_inspector = False 	# Suppresses showing the inspector - switched on in emfloatingwidgets
		
		##self.timer = QTimer()
		##QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.timeout)

		#self.cam = Camera2(self)
		#self.vdtools = EMViewportDepthTools(self)
		
		#self.rottarget = None
		#self.set_data(image1,image2)
		#self.ratio = 0.5
		##self.inspector.addIso()
	##def timeout(self):
		##self.updateGL()
	#def width(self):
		#try: return self.parent.width()
		#except: return 0
		
	#def height(self):
		#try: return self.parent.height()
		#except: return 0
	
	#def updateGL(self):
		#try: self.parent.updateGL()
		#except: pass
	
	#def eye_coords_dif(self,x1,y1,x2,y2,mdepth=True):
		#return self.vdtools.eye_coords_dif(x1,y1,x2,y2,mdepth)

	#def initializeGL(self):
		#glEnable(GL_NORMALIZE)
	
	#def render(self):
		#glPushMatrix()
		#self.cam.position(True)
		## the ones are dummy variables atm... they don't do anything
		#self.vdtools.update(1,1)
		#glPopMatrix()
		
		#self.cam.position()
		
		#for i in self.viewables:
			#glPushMatrix()
			#i.render()
			#glPopMatrix()

	#def resizeEvent(self, width, height):
		#for i in self.viewables:
			#i.resizeEvent()
	
	#def get_data_dims(self):
		#return [self.image.get_xsize(),self.image.get_ysize(),self.image.get_zsize()]

	#def updateRatio(self,ratio):
		#if ratio > 1:
			#self.ratio = 1
		#elif ratio < 0:
			#self.ratio = 0
		#else:
			#self.ratio = ratio
		
		#self.image = self.ratio*self.data1 + (1.0 - self.ratio)*self.data2
		#for i in self.viewables:
			#i.update_data(self.image)
		

	#def set_data(self,dataA,dataB):
		#if dataA == None or dataB == None: return
		#self.data1 = dataA
		#self.data2 = dataB
		#self.image = self.ratio*self.data1 + (1.0 - self.ratio)*self.data2
		#for i in self.viewables:
			#i.set_data(data)
			
		#self.resizeEvent(self.parent.width(),self.parent.height())
		##self.volume.set_data(data)
		
		#if self.inspector == None:
			#self.inspector=EMImageMorphInspector3D(self)
		#self.inspector.add_isosurface()
	
	#def show_inspector(self,force=0):
		#if self.suppress_inspector: return
		#if not force and self.inspector==None : return
		#self.init_inspector()
		#self.inspector.show()
		
	#def init_inspector(self):
		#if not self.inspector : self.inspector=EMImageMorphInspector3D(self)
	
	#def closeEvent(self,event) :
		##for i in self.viewables:
			##i.closeEvent(event)
		#if self.inspector: self.inspector.close()
		
	#def mouseMoveEvent(self, event):
		#self.cam.mouseMoveEvent(event)
		#if self.rottarget != None :
			#if event.buttons()&Qt.LeftButton:
				#self.rottarget.update_rotations(self.get_current_transform())
			#elif event.buttons()&Qt.RightButton:
				#self.rottarget.set_xy_trans(self.cam.cam_x, self.cam.cam_y)
		#self.updateGL()
	
			
	#def set_cam_z(self,z):
		#self.cam.set_cam_z( z )
		#self.updateGL()
		
	#def set_cam_y(self,y):
		#self.cam.set_cam_y( y )
		#self.updateGL()
		
	#def set_cam_x(self,x):
		#self.cam.set_cam_x( x )
		#self.updateGL()
	
	#def mouseReleaseEvent(self, event):
		#self.cam.mouseReleaseEvent(event)
		#self.updateGL()
			
	#def wheelEvent(self, event):
		#self.cam.wheelEvent(event)
		#if self.rottarget != None :
				#self.rottarget.set_scale(self.cam.scale)
		#self.updateGL()
	
	#def set_scale(self,val):
		#self.cam.scale = val
		#self.updateGL()
	
	#def mousePressEvent(self, event):
		#if event.button()==Qt.MidButton:
			#self.show_inspector(1)
		#else:
			#self.cam.mousePressEvent(event)
				
		#self.updateGL()
	
	#def get_render_dims_at_depth(self, depth):
		#return self.parent.get_render_dims_at_depth(depth)

	#def get_sundry_inspector(self):
		#return self.viewables[self.currentselection].get_inspector()
	
	#def add_sym(self):
		#sym = EM3DSymViewer(self)
		#self.viewables.append(sym)
		##if (len(self.viewables)==1):
			##sym.cam.default_z = -1.25*self.image.get_zsize()
			##sym.cam.cam_z = -1.25*self.image.get_zsize()
		##else:
			##pass
			###self.load_last_viewable_camera()
		#sym.set_radius(self.image.get_zsize()/2.0)
		#self.num_sym += 1
		#name = "Sym " + str(self.num_sym)
		#self.viewables[len(self.viewables)-1].set_name(name)
		#self.viewables[len(self.viewables)-1].set_rank(len(self.viewables))
		#self.currentselection = len(self.viewables)-1
		#self.updateGL()
	
	#def add_isosurface(self):
		#self.viewables.append(EMIsosurface(self.image,self))
		##self.load_last_viewable_camera()
		#self.num_iso += 1
		#name = "Isosurface " + str(self.num_iso)
		#self.viewables[len(self.viewables)-1].set_name(name)
		#self.viewables[len(self.viewables)-1].set_rank(len(self.viewables))
		#self.currentselection = len(self.viewables)-1
		#self.updateGL()
		
	#def add_volume(self):
		#self.viewables.append(EMVolumeModule(self.image,self))
		##self.load_last_viewable_camera()
		#self.num_vol += 1
		#name = "Volume " + str(self.num_vol)
		#self.viewables[len(self.viewables)-1].set_name(name)
		#self.viewables[len(self.viewables)-1].set_rank(len(self.viewables))
		#self.currentselection = len(self.viewables)-1
		#self.updateGL()
		
	#def add_slice_viewer(self):
		#self.viewables.append(EM3DSliceViewerModule(self.image,self))
		##self.load_last_viewable_camera()
		#self.num_sli += 1
		#name = "Slices " + str(self.num_sli)
		#self.viewables[len(self.viewables)-1].set_name(name)
		#self.viewables[len(self.viewables)-1].set_rank(len(self.viewables))
		#self.currentselection = len(self.viewables)-1
		#self.updateGL()
		
	#def load_last_viewable_camera(self):
		#return
		#size = len(self.viewables)
		#if ( size <= 1 ): return
		#self.viewables[size-1].set_camera(self.viewables[0].get_current_camera())

	#def rowChanged(self,row):
		#if ( row == self.currentselection ): return
		#self.currentselection=row
		#self.updateGL()
		
	#def get_current_idx(self):
		#return self.currentselection
		
	#def get_current_name(self):
		#if self.currentselection == -1 : return ""
		#elif self.currentselection >= len(self.viewables):
			#print "error, current seletion too large", self.currentselection,len(self.viewables)
			#return ""
		#return self.viewables[self.currentselection].get_name()
	
	#def getCurrentInspector(self):
		#if self.currentselection == -1 : return None
		#elif self.currentselection >= len(self.viewables):
			#print "error, current seletion too large", self.currentselection,len(self.viewables)
			#return None
		#return self.viewables[self.currentselection].get_inspector()
	
	#def delete_current(self, val):
		#if ( len(self.viewables) == 0 ): return
		
		#self.viewables.pop(val)
		#if (len(self.viewables) == 0 ) : 
			#self.currentselection = -1
		#elif ( len(self.viewables) == 1):
			#self.currentselection = 0
		#elif ( val == 0):
			#pass
		#else:
			#self.currentselection = val - 1
		
		
		## Need to set the rank appropriately
		#for i in range(0,len(self.viewables)):
			#self.viewables[i].set_rank(i+1)
		
		#self.updateGL()
	
	#def resizeEvent(self,width=0,height=0):
		#self.vdtools.set_update_P_inv()
		
	#def set_perspective(self,bool):
		#self.parent.set_perspective(bool)
		
	#def load_rotation(self,t3d):
		#self.cam.load_rotation(t3d)
		#self.updateGL()

	#def get_current_transform(self):
		#size = len(self.cam.t3d_stack)
		#return self.cam.t3d_stack[size-1]
	
	#def register_rotation_target(self, targ):
		#self.rottarget = targ
	
	#def get_start_z(self):
		#return self.parent.get_start_z()
	
	#def get_near_plane_dims(self):
		#return self.parent.get_near_plane_dims()
	
#class EMImageMorphInspector3D(QtGui.QWidget):
	#def __init__(self,target) :
		#QtGui.QWidget.__init__(self,None)
		#self.target=target
		
		
		#self.vbl = QtGui.QVBoxLayout(self)
		#self.vbl.setMargin(0)
		#self.vbl.setSpacing(6)
		#self.vbl.setObjectName("vbl")
		
		#self.hbl = QtGui.QHBoxLayout()
		#self.hbl.setMargin(2)
		#self.hbl.setSpacing(6)
		#self.hbl.setObjectName("hbl")
		
		##self.listwidget = QtGui.QListWidget(self)
		##self.vbl.addWidget(self.listwidget)
		
		#self.tabwidget = QtGui.QTabWidget(self)
		#self.vbl.addWidget(self.tabwidget)
		
		#self.ratio = ValSlider(self,(0.0,1.0),"Ratio:")
		#self.ratio.setObjectName("ratio")
		#self.ratio.setValue(0.5)
		#self.vbl.addWidget(self.ratio)
		
		#self.hbl_check = QtGui.QHBoxLayout()
		#self.hbl_check.setMargin(0)
		#self.hbl_check.setSpacing(6)
		#self.hbl_check.setObjectName("hbl_check")
		
		##self.advancedcheck = QtGui.QCheckBox("Advanced",self)
		##self.hbl_check.addWidget(self.advancedcheck)
		
		#self.hbl_buttons = QtGui.QHBoxLayout()
		#self.hbl_buttons.setMargin(0)
		#self.hbl_buttons.setSpacing(6)
		#self.hbl_buttons.setObjectName("hbl_buttons")
		
		#self.hbl_buttons2 = QtGui.QHBoxLayout()
		#self.hbl_buttons2.setMargin(0)
		#self.hbl_buttons2.setSpacing(6)
		#self.hbl_buttons2.setObjectName("hbl_buttons2")
		
		#self.addIso = QtGui.QPushButton("Isosurface")
		#self.hbl_buttons.addWidget(self.addIso)
		
		#self.addVol = QtGui.QPushButton("Volume")
		#self.hbl_buttons.addWidget(self.addVol)
		
		#self.addSli = QtGui.QPushButton("Slices")
		#self.hbl_buttons2.addWidget(self.addSli)
		
		#self.add_sym = QtGui.QPushButton("Sym")
		#self.hbl_buttons2.addWidget(self.add_sym)

		#self.vbl.addLayout(self.hbl_buttons)
		#self.vbl.addLayout(self.hbl_buttons2)
		
		#self.hbl_buttons3 = QtGui.QHBoxLayout()
		#self.delete = QtGui.QPushButton("Delete")
		#self.hbl_buttons3.addWidget(self.delete)
		#self.vbl.addLayout(self.hbl_buttons3)
		
		#self.vbl.addLayout(self.hbl_check)
		
		#self.setinspector = None
		
		#self.currentselection = -1
		#self.settingsrow = -2
		#self.targetidxmap = {}
		
		##self.advancedcheck.click()
		#self.insert_advance_tab()
		
		#QtCore.QObject.connect(self.addIso, QtCore.SIGNAL("clicked()"), self.add_isosurface)
		#QtCore.QObject.connect(self.addVol, QtCore.SIGNAL("clicked()"), self.add_volume)
		#QtCore.QObject.connect(self.addSli, QtCore.SIGNAL("clicked()"), self.add_slices)
		#QtCore.QObject.connect(self.add_sym, QtCore.SIGNAL("clicked()"), self.add_symmetry)
		#QtCore.QObject.connect(self.delete, QtCore.SIGNAL("clicked()"), self.delete_selection)
		#QtCore.QObject.connect(self.ratio, QtCore.SIGNAL("valueChanged"), self.target.updateRatio)
		##QtCore.QObject.connect(self.advancedcheck, QtCore.SIGNAL("stateChanged(int)"), self.advancedClicked)
		
		##QtCore.QObject.connect(self.listwidget, QtCore.SIGNAL("currentRowChanged(int)"), self.rowChanged)
		##QtCore.QObject.connect(self.tabwidget, QtCore.SIGNAL("currentChanged(int)"), self.tabChanged)
		
		
	#def insert_advance_tab(self):
		#if self.setinspector == None:
			#self.setinspector = EM3DAdvancedInspector(self.target, self)
			
		#self.target.register_rotation_target(self.setinspector)
		#self.setinspector.update_rotations(self.target.get_current_transform())
		#self.setinspector.set_scale(self.target.cam.scale)
		#self.tabwidget.addTab(self.setinspector,"Advanced")
		#self.settingsrow = self.tabwidget.count()-1
		#self.targetidxmap[self.settingsrow] = -1
		#self.tabwidget.setCurrentIndex(self.settingsrow)
	

	#def add_isosurface(self):
		#self.target.add_isosurface()
		#self.update_selection()
	
	#def add_symmetry(self):
		#self.target.add_sym()
		#self.update_selection()
	
	#def add_volume(self):
		#self.target.add_volume()
		#self.update_selection()
	
	#def update_selection(self):
		#n = self.tabwidget.count()
		#if n > 0: n = n - 1
		#self.tabwidget.insertTab(n, self.target.getCurrentInspector(), self.target.get_current_name())
		#self.targetidxmap[n] = self.target.currentselection
		#self.tabwidget.setCurrentIndex(n)

	#def add_slices(self):
		#self.target.add_slice_viewer()
		#self.update_selection()
	
	#def delete_selection(self):
		#idx = self.tabwidget.currentIndex()
		#n = self.tabwidget.count()
		#if n <= 1: return
		#if idx == n-1: return
		
		#self.tabwidget.removeTab(idx)
		#self.target.delete_current(self.targetidxmap[idx])

#class EMTransformPanel:
	#def __init__(self,target,parent):
		#self.target = target
		#self.parent = parent
		
		#self.label_src = QtGui.QLabel(parent)
		#self.label_src.setText('Rotation Convention')
		
		#self.src = QtGui.QComboBox(parent)
		#self.load_src_options(self.src)
		
		#self.x_label = QtGui.QLabel()
		#self.x_label.setText('x')
		
		#self.x_trans = QtGui.QDoubleSpinBox(parent)
		#self.x_trans.setMinimum(-10000)
		#self.x_trans.setMaximum(10000)
		#self.x_trans.setValue(0.0)
	
		#self.y_label = QtGui.QLabel()
		#self.y_label.setText('y')
		
		#self.y_trans = QtGui.QDoubleSpinBox(parent)
		#self.y_trans.setMinimum(-10000)
		#self.y_trans.setMaximum(10000)
		#self.y_trans.setValue(0.0)
		
		#self.z_label = QtGui.QLabel()
		#self.z_label.setText('z')
		
		#self.z_trans = QtGui.QDoubleSpinBox(parent)
		#self.z_trans.setMinimum(-10000)
		#self.z_trans.setMaximum(10000)
		#self.z_trans.setValue(0.0)
		
		#self.az = ValSlider(parent,(-360.0,360.0),"az",-1)
		#self.az.setObjectName("az")
		#self.az.setValue(0.0)
		
		#self.alt = ValSlider(parent,(-180.0,180.0),"alt",-1)
		#self.alt.setObjectName("alt")
		#self.alt.setValue(0.0)
		
		#self.phi = ValSlider(parent,(-360.0,360.0),"phi",-1)
		#self.phi.setObjectName("phi")
		#self.phi.setValue(0.0)
		
		#self.scale = ValSlider(parent,(0.01,30.0),"Zoom:")
		#self.scale.setObjectName("scale")
		#self.scale.setValue(1.0)
		
		#self.n3_showing = False
		
		#self.current_src = EULER_EMAN
		
		#QtCore.QObject.connect(self.az, QtCore.SIGNAL("valueChanged"), self.slider_rotate)
		#QtCore.QObject.connect(self.alt, QtCore.SIGNAL("valueChanged"), self.slider_rotate)
		#QtCore.QObject.connect(self.phi, QtCore.SIGNAL("valueChanged"), self.slider_rotate)
		#QtCore.QObject.connect(self.src, QtCore.SIGNAL("currentIndexChanged(QString)"), self.set_src)
		#QtCore.QObject.connect(self.scale, QtCore.SIGNAL("valueChanged"), target.set_scale)
		#QtCore.QObject.connect(self.x_trans, QtCore.SIGNAL("valueChanged(double)"), target.set_cam_x)
		#QtCore.QObject.connect(self.y_trans, QtCore.SIGNAL("valueChanged(double)"), target.set_cam_y)
		#QtCore.QObject.connect(self.z_trans, QtCore.SIGNAL("valueChanged(double)"), target.set_cam_z)
		
	#def slider_rotate(self):
		#self.target.load_rotation(self.get_current_rotation())
		
	#def get_current_rotation(self):
		#convention = self.src.currentText()
		#rot = {}
		#if ( self.current_src == EULER_SPIN ):
			#rot[self.az.getLabel()] = self.az.getValue()
			
			#n1 = self.alt.getValue()
			#n2 = self.phi.getValue()
			#n3 = self.n3.getValue()
			
			#norm = sqrt(n1*n1 + n2*n2 + n3*n3)
			
			#n1 /= norm
			#n2 /= norm
			#n3 /= norm
			
			#rot[self.alt.getLabel()] = n1
			#rot[self.phi.getLabel()] = n2
			#rot[self.n3.getLabel()] = n3
			
		#else:
			#rot[self.az.getLabel()] = self.az.getValue()
			#rot[self.alt.getLabel()] = self.alt.getValue()
			#rot[self.phi.getLabel()] = self.phi.getValue()
		
		#return Transform3D(self.current_src, rot)
	
	#def addWidgets(self,target):
		#self.hbl_trans = QtGui.QHBoxLayout()
		#self.hbl_trans.setMargin(0)
		#self.hbl_trans.setSpacing(6)
		#self.hbl_trans.setObjectName("Trans")
		#self.hbl_trans.addWidget(self.x_label)
		#self.hbl_trans.addWidget(self.x_trans)
		#self.hbl_trans.addWidget(self.y_label)
		#self.hbl_trans.addWidget(self.y_trans)
		#self.hbl_trans.addWidget(self.z_label)
		#self.hbl_trans.addWidget(self.z_trans)
		
		#target.addLayout(self.hbl_trans)
		
		#self.hbl_src = QtGui.QHBoxLayout()
		#self.hbl_src.setMargin(0)
		#self.hbl_src.setSpacing(6)
		#self.hbl_src.setObjectName("hbl")
		#self.hbl_src.addWidget(self.label_src)
		#self.hbl_src.addWidget(self.src)
		
		#target.addWidget(self.scale)
		#target.addLayout(self.hbl_src)
		#target.addWidget(self.az)
		#target.addWidget(self.alt)
		#target.addWidget(self.phi)
	
	#def set_src(self, val):
		#t3d = self.get_current_rotation()
		
		#if (self.n3_showing) :
			#self.maintab.vbl.removeWidget(self.n3)
			#self.n3.deleteLater()
			#self.n3_showing = False
			#self.az.setRange(-360,360)
			#self.alt.setRange(-180,180)
			#self.phi.setRange(-360,660)
		
		#if ( self.src_map[str(val)] == EULER_SPIDER ):
			#self.az.setLabel('phi')
			#self.alt.setLabel('theta')
			#self.phi.setLabel('psi')
		#elif ( self.src_map[str(val)] == EULER_EMAN ):
			#self.az.setLabel('az')
			#self.alt.setLabel('alt')
			#self.phi.setLabel('phi')
		#elif ( self.src_map[str(val)] == EULER_IMAGIC ):
			#self.az.setLabel('alpha')
			#self.alt.setLabel('beta')
			#self.phi.setLabel('gamma')
		#elif ( self.src_map[str(val)] == EULER_XYZ ):
			#self.az.setLabel('xtilt')
			#self.alt.setLabel('ytilt')
			#self.phi.setLabel('ztilt')
		#elif ( self.src_map[str(val)] == EULER_MRC ):
			#self.az.setLabel('phi')
			#self.alt.setLabel('theta')
			#self.phi.setLabel('omega')
		#elif ( self.src_map[str(val)] == EULER_SPIN ):
			#self.az.setLabel('omega')
			#self.alt.setRange(-1,1)
			#self.phi.setRange(-1,1)
			
			#self.alt.setLabel('n1')
			#self.phi.setLabel('n2')
			
			#self.n3 = ValSlider(self,(-360.0,360.0),"n3",-1)
			#self.n3.setRange(-1,1)
			#self.n3.setObjectName("n3")
			#self.maintab.vbl.addWidget(self.n3)
			#QtCore.QObject.connect(self.n3, QtCore.SIGNAL("valueChanged"), self.slider_rotate)
			#self.n3_showing = True
		
		#self.current_src = self.src_map[str(val)]
		#self.update_rotations(t3d)
	
	#def load_src_options(self,widgit):
		#self.load_src()
		#for i in self.src_strings:
			#widgit.addItem(i)
			
	#def load_src(self):
		## supported_rot_conventions
		#src_flags = []
		#src_flags.append(EULER_EMAN)
		#src_flags.append(EULER_SPIDER)
		#src_flags.append(EULER_IMAGIC)
		#src_flags.append(EULER_MRC)
		#src_flags.append(EULER_SPIN)
		#src_flags.append(EULER_XYZ)
		
		#self.src_strings = []
		#self.src_map = {}
		#for i in src_flags:
			#self.src_strings.append(str(i))
			#self.src_map[str(i)] = i
			
	#def update_rotations(self,t3d):
		#rot = t3d.get_rotation(self.src_map[str(self.src.itemText(self.src.currentIndex()))])
		
		#convention = self.src.currentText()
		#if ( self.src_map[str(convention)] == EULER_SPIN ):
			#self.n3.setValue(rot[self.n3.getLabel()],True)
		
		#self.az.setValue(rot[self.az.getLabel()],True)
		#self.alt.setValue(rot[self.alt.getLabel()],True)
		#self.phi.setValue(rot[self.phi.getLabel()],True)
		
	#def set_scale(self,newscale):
		#self.scale.setValue(newscale)
		
	#def set_xy_trans(self, x, y):
		#self.x_trans.setValue(x)
		#self.y_trans.setValue(y)
	

#class EM3DAdvancedInspector(QtGui.QWidget):
	#def __init__(self,target,parent=None):
		#QtGui.QWidget.__init__(self,None)
		#self.target=target
		#self.parent=parent
		
		#self.rotsliders = EMTransformPanel(target,self)
		
		#self.hbl = QtGui.QHBoxLayout()
		#self.hbl.setMargin(2)
		#self.hbl.setSpacing(6)
		#self.hbl.setObjectName("hbl")
		
		#self.vbl = QtGui.QVBoxLayout(self)
		#self.vbl.setMargin(0)
		#self.vbl.setSpacing(6)
		#self.vbl.setObjectName("vbl")
		
		#self.persbut = QtGui.QRadioButton("Perspective")
		#self.persbut.setChecked(True)
		#self.orthbut = QtGui.QRadioButton("Orthographic")
		
		#self.groupbox = QtGui.QVBoxLayout()
		#self.groupbox.addWidget(self.persbut)
		#self.groupbox.addWidget(self.orthbut)
		#self.hbl.addLayout(self.groupbox)
		
		#self.vbl.addLayout(self.hbl)
		
		#self.rotsliders.addWidgets(self.vbl)
		
		#QtCore.QObject.connect(self.persbut, QtCore.SIGNAL("pressed()"), self.perspective_clicked)
		#QtCore.QObject.connect(self.orthbut, QtCore.SIGNAL("pressed()"), self.ortho_clicked)
		
	#def update_rotations(self,t3d):
		#self.rotsliders.update_rotations(t3d)
		
	#def set_scale(self,val):
		#self.rotsliders.set_scale(val)
	
	#def set_xy_trans(self, x, y):
		#self.rotsliders.set_xy_trans(x,y)
		
	#def perspective_clicked(self):
		#self.target.set_perspective(True)
		
	#def ortho_clicked(self):
		#self.target.set_perspective(False)
		
	
	
# This is just for testing, of course
if __name__ == '__main__':
	print "emimage3dmorph.py is not currently maintained. We want its functionality to be avaialable, however the code needs a complete rewrite"
	#app = QtGui.QApplication(sys.argv)
	#window = EMImageMorph3D()
	#if len(sys.argv) != 3:
		#print "Error, you must specify two images"
		#exit(1)
 	#if len(sys.argv)==1 : 
		#pass
		#e = EMData()
		#e.set_size(64,64,64)
		#e.process_inplace('testimage.axes')
 		#window.set_data(e)

		## these lines are for testing shape rendering
## 		window.add_shape("a",["rect",.2,.8,.2,20,20,80,80,2])
## 		window.add_shape("b",["circle",.5,.8,.2,120,50,30.0,2])
## 		window.add_shape("c",["line",.2,.8,.5,20,120,100,200,2])
## 		window.add_shape("d",["label",.2,.8,.5,220,220,"Testing",14,1])
	#else :
		#if not os.path.exists(sys.argv[1]):
			#print "Error, input file %s does not exist" %sys.argv[1]
			#exit(1)
		#if not os.path.exists(sys.argv[2]):
			#print "Error, input file %s does not exist" %sys.argv[2]
			#exit(1)
		#a=EMData.read_images(sys.argv[1],[0])
		#b=EMData.read_images(sys.argv[2],[0])
		#window.set_data(a[0],b[0])
	#window2=EMParentWin(window)
	#window2.show()
	
#	w2=QtGui.QWidget()
#	w2.resize(256,128)
	
#	w3=ValSlider(w2)
#	w3.resize(256,24)
#	w2.show()
	
	#sys.exit(app.exec_())
