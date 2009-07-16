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
from math import *
from EMAN2 import *
import sys
import numpy
import weakref

from weakref import WeakKeyDictionary
from time import time
from PyQt4.QtCore import QTimer

from emimage3diso import EMIsosurfaceModule
from emimage3dvol import EMVolumeModule
from emimage3dslice import EM3DSliceViewerModule
from emimage3dsym import EM3DSymViewerModule
from e2eulerxplor import EMAsymmetricUnitViewer
from emlights import EMLightsInspectorBase,EMLightsDrawer

from emglobjects import Camera2, EMViewportDepthTools, Camera, EMImage3DGUIModule,EMGLProjectionViewMatrices,EMOpenGLFlagsAndTools
from emimageutil import EMEventRerouter, EMTransformPanel, EMParentWin
from emapplication import EMStandAloneApplication, EMGUIModule,get_application


MAG_INCREMENT_FACTOR = 1.1


class EMImage3DGeneralWidget(QtOpenGL.QGLWidget,EMEventRerouter,EMGLProjectionViewMatrices):
	def __init__(self, em_3d_module):
		assert(isinstance(em_3d_module,EMImage3DGUIModule))
		
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True)
		fmt.setDepth(1)
		fmt.setSampleBuffers(True)
		QtOpenGL.QGLWidget.__init__(self,fmt)
		EMEventRerouter.__init__(self,em_3d_module)
		EMGLProjectionViewMatrices.__init__(self)
		
		self.fov = 50 # field of view angle used by gluPerspective
		self.startz = 1
		self.endz = 5000
		self.cam = Camera()
		self.resize(480,480)
		
	def set_camera_defaults(self,data):
		if isinstance(data,EMData):
			self.cam.default_z = -1.25*data.get_xsize()
			self.cam.cam_z = -1.25*data.get_xsize()
		elif isinstance(data,float):
			self.cam.default_z = -1.25*data
			self.cam.cam_z = -1.25*data
		
	def set_data(self,data):
		self.target().set_data(data)
		self.set_camera_defaults(data)
		
		
	def initializeGL(self):
		glEnable(GL_NORMALIZE)
		glEnable(GL_LIGHT0)
		glEnable(GL_DEPTH_TEST)
		#print "Initializing"
		glLightfv(GL_LIGHT0, GL_AMBIENT, [0.3, 0.3, 0.3, 1.0])
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_POSITION, [0.1,.1,1.,0.])
		GL.glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST)
		glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,GL_TRUE)
		GL.glClearColor(0,0,0,0)
		#GL.glClearAccum(0,0,0,0)
	
		glShadeModel(GL_SMOOTH)
		
		glClearStencil(0)
		glEnable(GL_STENCIL_TEST)
		
	def paintGL(self):
		#glClear(GL_ACCUM_BUFFER_BIT)
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT )
		
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		self.cam.position()
		
		glPushMatrix()
		self.target().render()
		glPopMatrix()
		
		#glAccum(GL_ADD, self.target.glbrightness)
		#glAccum(GL_ACCUM, self.target.glcontrast)
		#glAccum(GL_RETURN, 1.0)
		
	def resizeGL(self, width, height):
		if width<=0 or height<=0 : return # this is fine, the window has be size interactively to zero, for example
		# just use the whole window for rendering
		glViewport(0,0,self.width(),self.height())
		
		# maintain the aspect ratio of the window we have
		self.aspect = float(self.width())/float(self.height())
		
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		# using gluPerspective for simplicity
		gluPerspective(self.fov,self.aspect,self.startz,self.endz)
		
		# switch back to model view mode
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		self.set_projection_view_update()
		self.target().resizeEvent()

	def get_start_z(self):
		return self.startz
	
	def get_near_plane_dims(self):
		height = 2.0 * self.startz*tan(self.fov/2.0*pi/180.0)
		width = self.aspect * height
		return [width,height]

	def show_inspector(self,force=0):
		self.target().show_inspector(self,force)

	def get_render_dims_at_depth(self,depth):
		# This function returns the width and height of the renderable 
		# area at the origin of the data volume
		height = -2*tan(self.fov/2.0*pi/180.0)*(depth)
		width = self.aspect*height
		
		return [width,height]

class EMImage3DWidget(QtOpenGL.QGLWidget,EMEventRerouter,EMGLProjectionViewMatrices):
	""" 
	A QT widget for rendering 3D EMData objects
	"""
	def __init__(self, image_3d_module, parent=None):
		self.target = None
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True)
		fmt.setDepth(True)
		fmt.setStencil(True)
		fmt.setSampleBuffers(True)
		QtOpenGL.QGLWidget.__init__(self,fmt)
		EMEventRerouter.__init__(self,image_3d_module)
		EMGLProjectionViewMatrices.__init__(self)
		self.aspect=1.0
		self.fov = 50 # field of view angle used by gluPerspective
		self.d = 0
		self.zwidth = 0
		self.perspective = True

		self.cam = Camera()

		self.resize(480,480)
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
		
		max = self.zwidth
		if self.yheight > max: max = self.yheight
		if self.xwidth > max: mas = self.xwidth
		self.startz = self.d - 2.0*max
		self.endz = self.d + 2.0*max
	def set_data(self,data):
		self.target().set_data(data)
		if ( data != None and isinstance(data,EMData)):
			self.set_cam_z(self.fov,data)
			
		self.resize(640,640)
		
	def initializeGL(self):
		glEnable(GL_LIGHTING)
		glEnable(GL_LIGHT0)
		#glEnable(GL_LIGHT1)
		glEnable(GL_DEPTH_TEST)
		glLightfv(GL_LIGHT0, GL_AMBIENT, [0.0, 0.0, 0.0, 1.0])
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0,1.0,1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_POSITION, [0.1,.1,1.,0.])
		
#		glLightfv(GL_LIGHT1, GL_AMBIENT, [0.2, 0.2, 0.0, 1.0])
#		glLightfv(GL_LIGHT1, GL_DIFFUSE, [0.8,0.8,0.8, 1.0])
#		glLightfv(GL_LIGHT1, GL_SPECULAR, [0.0, 0.0, 0.0, 1.0])
#		glLightfv(GL_LIGHT1, GL_POSITION, [0,0,1,1]) # set the is self.radius when it's known
#		glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, [0,0,-1])
#		glLightfv(GL_LIGHT1, GL_QUADRATIC_ATTENUATION,0.0037)
		#GL_SPOT_DIRECTION,GL_SPOT_CUTOFF,GL_QUADRATIC_ATTENUATION
		
		
		GL.glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST)
		glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,GL_TRUE)
		glShadeModel(GL_SMOOTH)
		#glLightModelfv(GL_LIGHT_MODEL_AMBIENT, [0.1,0.1,0.1,1.0]);
		
		glClearStencil(0)
		glEnable(GL_STENCIL_TEST)
		GL.glClearColor(0,0,0,0)
		try:
			self.target().initializeGL()
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

		if ( self.target != None ):
			self.target().render()


	def resizeGL(self, width, height):
		# just use the whole window for rendering
		
		if width == 0 or height == 0: return # this is fine
		
		glViewport(0,0,self.width(),self.height())
		
		# maintain the aspect ratio of the window we have
		#self.aspect = float(self.width())/float(self.height())
		
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		
		#self.startz = self.d - 2.0*self.zwidth
		#self.endz = self.d + 2.0*self.zwidth
		
		if (self.zwidth == 0):
			# We've received  a resize event but no data has been set
			# in which case nothing is being rendered.
			# Therefore just leave the identity as the projection matrix.
			# This is an exceptional circumstance which probably 
			# highlights the need for some redesigning (d.woolford)
			glMatrixMode(GL_MODELVIEW)
			#glLoadIdentity()
			return
		
		#if self.perspective:
			# using gluPerspective for simplicity
			
		self.load_perspective()
		#else:
			#self.load_orthographic()
			#self.xwidth = self.aspect*self.yheight
			#glOrtho(-self.xwidth/2.0,self.xwidth/2.0,-self.yheight/2.0,self.yheight/2.0,self.startz,self.endz)
			
		# switch back to model view mode
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		if (self.target != None):
			try: self.target().resizeEvent(width,height)
			except: pass
		self.set_projection_view_update()
		self.updateGL()

	def load_orthographic(self):
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		self.aspect = float(self.width())/float(self.height())
		self.xwidth = self.aspect*self.yheight
		if self.xwidth == 0 or self.yheight == 0: return # probably startup
		glOrtho(-self.xwidth/2.0,self.xwidth/2.0,-self.yheight/2.0,self.yheight/2.0,self.startz,self.endz)
		glMatrixMode(GL_MODELVIEW)
		
	def load_perspective(self):
		self.aspect = float(self.width())/float(self.height())
		
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		if self.startz < 0: self.startz = 1
		gluPerspective(self.fov,self.aspect,self.startz,self.endz)
		glMatrixMode(GL_MODELVIEW)
		
		
	#def set_perspective(self,bool):
		#self.perspective = bool
		#self.resizeGL(self.width(),self.height())
		
	def get_start_z(self):
		return self.startz
	
	def get_near_plane_dims(self):
		if self.perspective:
			height = 2.0*self.startz * tan(self.fov/2.0*pi/180.0)
			width = self.aspect * height
			return [width,height]
		else:
			return [self.xwidth,self.yheight]
		
	def set_camera_defaults(self,data):
		self.cam.default_z = -1.25*data.get_xsize()
		self.cam.cam_z = -1.25*data.get_xsize()
	
	def set_data(self,data):
		self.target().set_data(data)
		self.set_camera_defaults()
	
	def show_inspector(self,force=0):
		self.target().show_inspector()

	def get_render_dims_at_depth(self, depth):
		# This function returns the width and height of the renderable 
		# area at the origin of the data volume
		height = -2*tan(self.fov/2.0*pi/180.0)*(depth)
		width = self.aspect*height
		return [width,height]		
		
class EMImage3DModule(EMLightsDrawer,EMImage3DGUIModule):
	
	def get_qt_widget(self):
		if self.qt_context_parent == None:	
			from emimageutil import EMParentWin
			self.under_qt_control = True
			self.gl_context_parent = EMImage3DWidget(self)
			self.qt_context_parent = EMParentWin(self.gl_context_parent)
			self.gl_widget = self.gl_context_parent
			
			for i in self.viewables:
				i.set_qt_context_parent(self.qt_context_parent)
				i.set_gl_context_parent(self.gl_context_parent)
				i.set_gl_widget(self.gl_context_parent)
				i.under_qt_control = True
				
			if isinstance(self.data,EMData):
				self.gl_context_parent.set_cam_z(self.gl_context_parent.get_fov(),self.data)
			
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
	
	allim=WeakKeyDictionary()
	def __init__(self, image=None,application=None):
		self.viewables = []
		EMImage3DGUIModule.__init__(self,application,ensure_gl_context=True)
		EMLightsDrawer.__init__(self)
		EMImage3DModule.allim[self] = 0
		self.currentselection = -1
		self.inspector = None
		#self.isosurface = EMIsosurfaceModule(image,self)
		#self.volume = EMVolumeModule(image,self)
		self.viewables = []
		self.num_iso = 0
		self.num_vol = 0
		self.num_sli = 0
		self.num_sym = 0
		self.suppress_inspector = False 	
		self.cam = Camera2(self)
		self.vdtools = EMViewportDepthTools(self)
		
		if image != None: self.set_data(image)
			
		self.em_qt_inspector_widget = None
		
		self.last_window_width = -1 # used for automatic resizing from the desktop
		self.last_window_height = -1 # used for automatic resizing from the desktop
		
		self.file_name = None
		
		self.emit_events = False
		
		self.perspective = False
		
	def __del__(self):
		if self.under_qt_control and not self.dont_delete_parent:
			self.qt_context_parent.deleteLater()
		self.core_object.deleteLater()
	
	
	def enable_emit_events(self,val=True):
		for v in self.viewables: v.enable_emit_events(val)
		self.emit_events = val
		self.cam.enable_emit_events(val)
		
	def is_emitting(self): return self.emit_events
	
	def get_emit_signals_and_connections(self):
		ret = {}
		for v in self.viewables: ret.update(v.get_emit_signals_and_connections())
		ret.update(self.cam.get_emit_signals_and_connections())
		ret.update({"set_perspective":self.set_perspective})
		
		return ret
		
	def set_file_name(self,name):
		self.file_name = name
		if self.qt_context_parent != None:
			self.qt_context_parent.setWindowTitle(remove_directories_from_name(self.file_name))
	
	def width(self):
		try: return self.gl_widget.width()
		except: return 0
		
	def height(self):
		try: return self.gl_widget.height()
		except: return 0
	
	#def updateGL(self):
		#from emfloatingwidgets import EM3DGLWindow
		#if self.gl_widget != None and not isinstance(self.gl_widget,EM3DGLWindow):
			#self.gl_widget.updateGL()
	
	def eye_coords_dif(self,x1,y1,x2,y2,mdepth=True):
		return self.vdtools.eye_coords_dif(x1,y1,x2,y2,mdepth)
	
	def initializeGL(self):
		glEnable(GL_NORMALIZE)
	
	def render(self):
		self.image_change_count = self.data.get_changecount() # this is important when the user has more than one display instance of the same image, for instance in e2.py if 
		glPushMatrix()
		self.cam.position(True)
		# the ones are dummy variables atm... they don't do anything
		self.vdtools.update(1,1)
		glPopMatrix()
		
		dz = None
		if not self.perspective:
			glMatrixMode(GL_PROJECTION)
			glPushMatrix() 
			self.gl_context_parent.load_orthographic()
			glMatrixMode(GL_MODELVIEW)
		
		glPushMatrix()
		self.cam.position()
		
		for i in self.viewables:
			glPushMatrix()
			i.render()
			glPopMatrix()
		glPopMatrix()
		
		
		glPushMatrix()
		self.cam.translate_only()
		EMLightsDrawer.draw(self)
		glPopMatrix()
		
		if not self.perspective:
			glMatrixMode(GL_PROJECTION)
			glPopMatrix()
			glMatrixMode(GL_MODELVIEW)
		
	def resizeEvent(self, width, height):
		for i in self.viewables:
			i.resizeEvent()
	
	def get_data_dims(self):
		if self.data != None:
			return [self.data.get_xsize(),self.data.get_ysize(),self.data.get_zsize()]
		else: return [0,0,0]

	def set_data(self,data,file_name="",replace=True):
		self.file_name = file_name # fixme fix this later
		if self.qt_context_parent != None:
			self.qt_context_parent.setWindowTitle(remove_directories_from_name(self.file_name))

		if data == None: return
		self.data = data
		self.data.process_inplace("normalize.edgemean")
		
		nx,ny,nz = self.data.get_xsize(),self.data.get_ysize(),self.data.get_zsize()
		
		self.radius = nz
		if ny > nz: self.radius = ny
		if nx > ny: self.radius = nx
		self.radius /= 2
		#for i in self.viewables:
			#i.set_data(data)
		
		if isinstance(self.gl_context_parent,EMImage3DWidget):
			self.resizeEvent(self.gl_context_parent.width(),self.gl_context_parent.height())
			self.gl_context_parent.set_cam_z(self.gl_context_parent.get_fov(),self.data)
		
		if self.inspector == None:
			self.inspector=EMImageInspector3D(self)
		
		if replace:
			self.inspector.delete_all()
		self.inspector.add_isosurface()
	
	def get_inspector(self):
		if not self.inspector :  self.inspector=EMImageInspector3D(self)
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

	def get_render_dims_at_depth(self, depth):
		return self.gl_context_parent.get_render_dims_at_depth(depth)

	def get_sundry_inspector(self):
		return self.viewables[self.currentselection].get_inspector()
	
	def add_sym(self):
		# the difference between the EMAsymmetricUnitViewer and the EM3DSymViewerModule
		# is only that the EMAsymmetricUnitViewer will look in the current directory for refinement directories and
		# display related information. Simply change from one to the other if you don't like it
		module = EMAsymmetricUnitViewer(None,True,False)
		#module = EM3DSymViewerModule(None,True,False)
		module.set_radius(self.radius)
		self.num_sym += 1
		self.__add_module(module,self.num_sym)
	
	def add_isosurface(self):
		module = EMIsosurfaceModule(self.data,None,False,False)
		self.num_iso += 1
		self.__add_module(module,self.num_iso)
		
	def add_volume(self):
		module = EMVolumeModule(self.data,None,False,False)
		self.num_vol += 1
		self.__add_module(module,self.num_vol)
	
	def add_slice_viewer(self):
		module = EM3DSliceViewerModule(self.data,None,False,False)
		self.num_sli += 1
		self.__add_module(module,self.num_sli)
	
	def __add_module(self,module,num=0):
		
		module.set_qt_context_parent(self.qt_context_parent)
		module.set_gl_context_parent(self.gl_context_parent)
		
		module.set_gl_widget(self.gl_context_parent)
		module.set_dont_delete_parent() # stops a RunTimeError
		module.under_qt_control = self.under_qt_control
		
		self.viewables.append(module)
		module.set_gl_context_parent(self.gl_context_parent)
		name = module.get_type()+" " + str(num)
		self.viewables[len(self.viewables)-1].set_name(name)
		self.viewables[len(self.viewables)-1].set_rank(len(self.viewables))
		self.currentselection = len(self.viewables)-1
		self.updateGL()
	
	
	def __set_module_contexts(self):
		for v in self.viewables:
			v.set_qt_context_parent(self.qt_context_parent)
			v.set_gl_context_parent(self.gl_context_parent)
			v.set_gl_widget(self.gl_context_parent)
			#self.application.register_qt_emitter(v,self.application.get_qt_emitter(self))
	
	def load_last_viewable_camera(self):
		return
		size = len(self.viewables)
		if ( size <= 1 ): return
		self.viewables[size-1].set_camera(self.viewables[0].get_current_camera())

	def rowChanged(self,row):
		if ( row == self.currentselection ): return
		self.currentselection=row
		self.updateGL()
		
	def get_current_idx(self):
		return self.currentselection
		
	def get_current_name(self):
		if self.currentselection == -1 : return ""
		elif self.currentselection >= len(self.viewables):
			print "error, current selection too large", self.currentselection,len(self.viewables)
			return ""
		return self.viewables[self.currentselection].get_name()
	
	def get_current_inspector(self):
		if self.currentselection == -1 : return None
		elif self.currentselection >= len(self.viewables):
			print "error, current selection too large", self.currentselection,len(self.viewables)
			return None
		return self.viewables[self.currentselection].get_inspector()
	
	def delete_current(self, val):
		if ( len(self.viewables) == 0 ): return
		
		v = self.viewables.pop(val)
		
		
		#self.application.deregister_qt_emitter(v)
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
			self.viewables[i].set_rank(i+1)
			
	
	def resizeEvent(self,width=0,height=0):
		self.vdtools.set_update_P_inv()
	
	def resize_event(self,width,height):
		if self.last_window_width == -1:
			self.last_window_width = width
			self.last_window_height = height
		else:
			height_scale = height/float(self.last_window_height)
			width_scale = width/float(self.last_window_width)
			
			if height_scale < width_scale: width_scale = height_scale
			#print width_scale, "is the factor"
			self.cam.scale *= width_scale
			self.last_window_width = width
			self.lsat_window_height = height
	
	def set_perspective(self,bool):
		self.perspective = bool
		if self.emit_events: self.emit(QtCore.SIGNAL("set_perspective"),bool)
		self.updateGL()
		#self.gl_context_parent.set_perspective(bool)
		
	def load_rotation(self,t3d):
		self.cam.load_rotation(t3d)
		self.updateGL()

	def get_current_transform(self):
		size = len(self.cam.t3d_stack)
		return self.cam.t3d_stack[size-1]
	
	def get_start_z(self):
		return self.gl_context_parent.get_start_z()
	
	def get_near_plane_dims(self):
		return self.gl_context_parent.get_near_plane_dims()
	
class EMImageInspector3D(QtGui.QWidget):
	def get_desktop_hint(self):
		return "inspector"
	
	def set_directional_light_dir(self,d):
		self.advanced_tab.set_directional_light_dir(d)
	
	def set_positional_light_pos(self,d):
		self.advanced_tab.set_positional_light_pos(d)
		
	def set_positional_light_dir(self,d):
		self.advanced_tab.set_positional_light_dir(d)
	
	def __init__(self,target) :
		QtGui.QWidget.__init__(self,None)
		self.target=weakref.ref(target)
		self.setWindowIcon(QtGui.QIcon(get_image_directory() +"desktop.png"))
		
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
		
		glflags = EMOpenGLFlagsAndTools()
		if glflags.npt_textures_unsupported(): self.addVol.setEnabled(False)
		
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
		
		self.advanced_tab = None
		
		self.currentselection = -1
		self.settingsrow = -2
		self.targetidxmap = {}

		self.insert_advance_tab()
		
		QtCore.QObject.connect(self.addIso, QtCore.SIGNAL("clicked()"), self.add_isosurface)
		QtCore.QObject.connect(self.addVol, QtCore.SIGNAL("clicked()"), self.add_volume)
		QtCore.QObject.connect(self.addSli, QtCore.SIGNAL("clicked()"), self.add_slices)
		QtCore.QObject.connect(self.add_sym, QtCore.SIGNAL("clicked()"), self.add_symmetry)
		QtCore.QObject.connect(self.delete, QtCore.SIGNAL("clicked()"), self.delete_selection)
		
	def update_rotations(self,t3d):
		self.advanced_tab.update_rotations(t3d)
	
	def set_scale(self,val):
		self.advanced_tab.set_scale(val)
	
	def set_xy_trans(self, x, y):
		self.advanced_tab.set_xy_trans(x,y)
	
	def set_xyz_trans(self,x,y,z):
		self.advanced_tab.set_xyz_trans(x,y,z)
	
	def insert_advance_tab(self):
		if self.advanced_tab == None:
			self.advanced_tab = EM3DAdvancedInspector(self.target(), self)
			
		self.advanced_tab.update_rotations(self.target().get_current_transform())
		self.advanced_tab.set_scale(self.target().cam.scale)
		self.tabwidget.addTab(self.advanced_tab,"Advanced")
		self.settingsrow = self.tabwidget.count()-1
		self.targetidxmap[self.settingsrow] = -1
		self.tabwidget.setCurrentIndex(self.settingsrow)
	
	def add_isosurface(self):
		self.target().add_isosurface()
		self.update_selection()
	
	def add_symmetry(self):
		self.target().add_sym()
		self.update_selection()
	
	def add_volume(self):
		self.target().add_volume()
		self.update_selection()
	
	def update_selection(self):
		n = self.tabwidget.count()
		if n > 0: n = n - 1
		self.tabwidget.insertTab(n, self.target().get_current_inspector(), self.target().get_current_name())
		self.targetidxmap[n] = self.target().currentselection
		self.tabwidget.setCurrentIndex(n)

	def add_slices(self):
		self.target().add_slice_viewer()
		self.update_selection()
	
	def delete_selection(self):
		idx = self.tabwidget.currentIndex()
		n = self.tabwidget.count()
		if n <= 1: return
		if idx == n-1: return
		
		self.tabwidget.removeTab(idx)
		self.target().delete_current(self.targetidxmap[idx])
		
		self.target().updateGL()
	
	def delete_all(self):
		n = self.tabwidget.count()
		if n <= 1: return
		
		for idx in range(n-2,-1,-1):
			self.tabwidget.removeTab(idx)
			self.target().delete_current(self.targetidxmap[idx])
	


class EM3DAdvancedInspector(QtGui.QWidget,EMLightsInspectorBase):
	
	
	def __init__(self,target,parent=None):
		QtGui.QWidget.__init__(self,None)
		EMLightsInspectorBase.__init__(self)
		self.target=weakref.ref(target)
		self.parent=weakref.ref(parent)

		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")

		self.tabwidget = QtGui.QTabWidget()
		self.tabwidget.addTab(self.get_main_tab(), "Transform")
		self.tabwidget.addTab(self.get_light_tab(), "Lights")
		
		#self.tabwidget.addTab(self.get_GL_tab(),"GL")
		
		self.vbl.addWidget(self.tabwidget)
		

		QtCore.QObject.connect(self.persbut, QtCore.SIGNAL("pressed()"), self.perspective_clicked)
		QtCore.QObject.connect(self.orthbut, QtCore.SIGNAL("pressed()"), self.ortho_clicked)
	
	
	def get_main_tab(self):
		self.maintab = QtGui.QWidget()
		maintab = self.maintab
		maintab.vbl = QtGui.QVBoxLayout(self.maintab)
		
		self.hbl = QtGui.QHBoxLayout()
		self.hbl.setMargin(2)
		self.hbl.setSpacing(6)
		self.hbl.setObjectName("hbl")
		
		self.persbut = QtGui.QRadioButton("Perspective")
		
		self.orthbut = QtGui.QRadioButton("Orthographic")
		self.orthbut.setChecked(True)
		
		self.groupbox = QtGui.QVBoxLayout()
		self.groupbox.addWidget(self.persbut)
		self.groupbox.addWidget(self.orthbut)
		
		self.viewingvol = QtGui.QGroupBox("Viewing Volume")
		self.viewingvol.setLayout(self.groupbox)
		
		self.hbl.addWidget(self.viewingvol)
		
		maintab.vbl.addLayout(self.hbl)
		
		self.rotation_sliders = EMTransformPanel(self.target(),self)
		self.rotation_sliders.addWidgets(maintab.vbl)
		
		return self.maintab
		
	def get_transform_layout(self):
		return self.vbl
		
	def update_rotations(self,t3d):
		self.rotation_sliders.update_rotations(t3d)
		
	def set_scale(self,val):
		self.rotation_sliders.set_scale(val)
	
	def set_xy_trans(self, x, y):
		self.rotation_sliders.set_xy_trans(x,y)
	
	def set_xyz_trans(self, x, y,z):
		self.rotation_sliders.set_xyz_trans(x,y,z)
	
	def perspective_clicked(self):
		self.target().set_perspective(True)
		
	def ortho_clicked(self):
		self.target().set_perspective(False)
	
if __name__ == '__main__':
	em_app = EMStandAloneApplication()
	window = EMImage3DModule(application=em_app)
	
	if len(sys.argv)==1 : 
		data = []
		#for i in range(0,200):
		e = test_image_3d(1,size=(64,64,64))
		window.set_data(e)
	else :
		a=EMData(sys.argv[1])
		window.set_file_name(sys.argv[1])
		window.set_data(a)
		
	em_app.show()
	em_app.execute()

