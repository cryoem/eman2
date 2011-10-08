#!/usr/bin/env python

#
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

from EMAN2 import *
from OpenGL import GL, GLU, GLUT
from OpenGL.GL import *
from OpenGL.GLU import *
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import QTimer, Qt
from e2eulerxplor import EMEulerExplorer
from emglobjects import Camera, Camera2, EMGLWidget, EMViewportDepthTools, EMGLProjectionViewMatrices, EMOpenGLFlagsAndTools
from emimage3diso import EMIsosurfaceModel
from emimage3dslice import EM3DSliceModel
from emimage3dsym import EM3DSymModel
from emimage3dvol import EMVolumeModel
from emimageutil import EMTransformPanel
from emlights import EMLightsInspectorBase, EMLightsDrawer
from math import *
import weakref





MAG_INCREMENT_FACTOR = 1.1

class EMImage3DWidget(EMGLWidget, EMLightsDrawer, EMGLProjectionViewMatrices):
	""" 
	A QT widget for rendering 3D EMData objects
	"""
	allim=weakref.WeakKeyDictionary()
	def add_model(self,model,num=0):
		model.set_gl_widget(self)
		model.set_dont_delete_parent() # stops a RunTimeError
		
		if self.viewables: 
			#TODO: find a better way to work with cameras
			model.cam.scale = self.viewables[0].cam.scale #Make the new model have the same scale as the first "viewable"
			model.cam.t3d_stack = self.viewables[0].cam.t3d_stack[:] #Make the new model have the same orientation and position as the first "viewable"
			model.get_inspector().update_rotations(model.cam.t3d_stack[-1])
			model.get_inspector().set_xyz_trans(model.cam.cam_x,model.cam.cam_y,model.cam.cam_z)
			model.get_inspector().set_scale(model.cam.scale)
		
		self.viewables.append(model)
		name = model.get_type()+" " + str(num)
		self.viewables[-1].set_name(name)
		self.viewables[-1].set_rank(len(self.viewables))
		self.currentselection = len(self.viewables)-1
		self.updateGL()
		
	def __init__(self, parent=None, image=None,application=None,winid=None):
		EMImage3DWidget.allim[self] = 0
		EMGLWidget.__init__(self,parent)
		EMLightsDrawer.__init__(self)
		EMGLProjectionViewMatrices.__init__(self)
		
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True)
		fmt.setDepth(True)
		fmt.setStencil(True)
		fmt.setSampleBuffers(True)
		self.setFormat(fmt)
		
		self.aspect=1.0
		self.fov = 50 # field of view angle used by gluPerspective
		self.d = 0
		self.zwidth = 0
		self.yheight = None
		
		self.data = None # should eventually be an EMData object

#		self.cam = Camera()
		self.cam = Camera2(self)
		self.cam.cam_z = -250
		
		self.resize(480,480)
		self.startz = 1
		self.endz = 500
		
		self.currentselection = -1
		self.inspector = None
		
		self.viewables = []
		self.num_iso = 0
		self.num_vol = 0
		self.num_sli = 0
		self.num_sym = 0

		self.vdtools = EMViewportDepthTools(self)
		
		self.last_window_width = -1 # used for automatic resizing from the desktop
		self.last_window_height = -1 # used for automatic resizing from the desktop
		
		self.file_name = None		
		self.emit_events = False		
		self.perspective = False

		if image != None: 
			self.set_data(image)
		
		#From get_qt_widget...
		if isinstance(self.data,EMData):
			self.set_cam_z_from_fov_image(self.get_fov(),self.data)
		
		self.qt_parent.setWindowIcon(QtGui.QIcon(get_image_directory() +"single_image_3d.png"))
		#End from get_qt_widget
		
		self.updateGL() #Solves "error, OpenGL seems not to be initialized" message
		
	def __set_model_contexts(self):
		for v in self.viewables:
			v.set_gl_widget(self)
	def add_isosurface(self):
		model = EMIsosurfaceModel(self, self.data, False)
		self.num_iso += 1
		self.add_model(model,self.num_iso)
	def add_slice_viewer(self):
		model = EM3DSliceModel(self, self.data)
		self.num_sli += 1
		self.add_model(model,self.num_sli)
	def add_sym(self):
		# the difference between the EMEulerExplorer and the EM3DSymModel
		# is only that the EMEulerExplorer will look in the current directory for refinement directories and
		# display related information. Simply change from one to the other if you don't like it
		model = EMEulerExplorer(self,True,False)
		#model = EM3DSymModel(self)
		model.set_radius(self.radius)
		self.num_sym += 1
		self.add_model(model,self.num_sym)
	def add_volume(self):
		model = EMVolumeModel(self, self.data)
		self.num_vol += 1
		self.add_model(model,self.num_vol)
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
	def enable_emit_events(self,val=True):
		for v in self.viewables: v.enable_emit_events(val)
		self.emit_events = val
		self.cam.enable_emit_events(val)
	def eye_coords_dif(self,x1,y1,x2,y2,mdepth=True):
		return self.vdtools.eye_coords_dif(x1,y1,x2,y2,mdepth)
	def get_current_idx(self):
		return self.currentselection
	def get_current_inspector(self):
		if self.currentselection == -1 : return None
		elif self.currentselection >= len(self.viewables):
			print "error, current selection too large", self.currentselection,len(self.viewables)
			return None
		return self.viewables[self.currentselection].get_inspector()		
	def get_current_name(self):
		if self.currentselection == -1 : return ""
		elif self.currentselection >= len(self.viewables):
			print "error, current selection too large", self.currentselection,len(self.viewables)
			return ""
		return self.viewables[self.currentselection].get_name()
	def get_current_transform(self):
		size = len(self.cam.t3d_stack)
		return self.cam.t3d_stack[size-1]
	def get_data_dims(self):
		if self.data != None:
			return [self.data.get_xsize(),self.data.get_ysize(),self.data.get_zsize()]
		else: return [0,0,0]
	def get_emit_signals_and_connections(self):
		ret = {}
		for v in self.viewables: ret.update(v.get_emit_signals_and_connections())
		ret.update(self.cam.get_emit_signals_and_connections())
		ret.update({"set_perspective":self.set_perspective})
		
		return ret
	def get_fov(self):
		return self.fov
	def get_inspector(self):
		if not self.inspector :  self.inspector=EMImageInspector3D(self)
		return self.inspector
	def get_near_plane_dims(self):
		if self.perspective:
			height = 2.0*self.startz * tan(self.fov/2.0*pi/180.0)
			width = self.aspect * height
			return [width,height]
		else:
			return [self.xwidth,self.yheight]
	def get_render_dims_at_depth(self, depth):
		# This function returns the width and height of the renderable 
		# area at the origin of the data volume
		height = -2*tan(self.fov/2.0*pi/180.0)*(depth)
		width = self.aspect*height
		return [width,height]
	def get_start_z(self):
		return self.startz
	def get_sundry_inspector(self):
		return self.viewables[self.currentselection].get_inspector()
	def initializeGL(self):
		glEnable(GL_LIGHTING)
		glEnable(GL_LIGHT0)
		#glEnable(GL_LIGHT1)
		glEnable(GL_DEPTH_TEST)
		glLightfv(GL_LIGHT0, GL_AMBIENT, [0.0, 0.0, 0.0, 1.0])
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0,1.0,1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_POSITION, [0.1,.1,1.,0.])
		
		glLightfv(GL_LIGHT1, GL_AMBIENT, [0.0, 0.0, 0.0, 1.0])
		glLightfv(GL_LIGHT1, GL_DIFFUSE, [0.5,0.5,0.5, 1.0])
		glLightfv(GL_LIGHT1, GL_SPECULAR, [0.0, 0.0, 0.0, 1.0])
		glLightfv(GL_LIGHT1, GL_POSITION, [0,0,1,1]) # set the is self.radius when it's known
		glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, [0,0,-1])
		glLightfv(GL_LIGHT1, GL_QUADRATIC_ATTENUATION,0.0037)
		
		glLightfv(GL_LIGHT2, GL_AMBIENT, [0.0, 0.0, 0.0, 1.0])
		glLightfv(GL_LIGHT2, GL_DIFFUSE, [0.5,0.5,0.5, 1.0])
		glLightfv(GL_LIGHT2, GL_SPECULAR, [0.0, 1.0, 0.0, 1.0])
		glLightfv(GL_LIGHT2, GL_POSITION, [0.1,.1,1.,0.])
		
		glLightfv(GL_LIGHT3, GL_AMBIENT, [0.0, 0.0, 0.0, 1.0])
		glLightfv(GL_LIGHT3, GL_DIFFUSE, [0.5,0.5,0.5, 1.0])
		glLightfv(GL_LIGHT3, GL_SPECULAR, [0.0, 1.0, 0.0, 1.0])
		glLightfv(GL_LIGHT3, GL_POSITION, [0.1,.1,1.,0.])
		#GL_SPOT_DIRECTION,GL_SPOT_CUTOFF,GL_QUADRATIC_ATTENUATION
		
		
		GL.glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST)
		glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,GL_TRUE)
		glShadeModel(GL_SMOOTH)
		#glLightModelfv(GL_LIGHT_MODEL_AMBIENT, [0.1,0.1,0.1,1.0]);
		
		glClearStencil(0)
		glEnable(GL_STENCIL_TEST)
		GL.glClearColor(0,0,0,0)
		
		
		glEnable(GL_NORMALIZE)
	def is_emitting(self): return self.emit_events
	def load_last_viewable_camera(self):
		return
		size = len(self.viewables)
		if ( size <= 1 ): return
		self.viewables[size-1].set_camera(self.viewables[0].get_current_camera())
	def load_orthographic(self):
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		
		if self.yheight == None: self.yheight = self.height()
		
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
	def load_rotation(self,t3d):
		self.cam.load_rotation(t3d)
		self.updateGL()
	def mouseMoveEvent(self, event):
		if self.current_mouse_mode:
			EMLightsDrawer.mouseMoveEvent(self, event)
		else:
			for model in self.viewables:
				try:
					model.mouseMoveEvent(event)
				except AttributeError, e:
					pass
		self.updateGL()
	def mousePressEvent(self, event):
		if event.button()==Qt.MidButton or (event.button()==Qt.LeftButton and event.modifiers()&Qt.AltModifier):
			self.show_inspector()
		if self.current_mouse_mode:
			EMLightsDrawer.mousePressEvent(self, event)
		else:
			for model in self.viewables:
				try:
					model.mousePressEvent(event)
				except AttributeError, e:
					pass
		self.updateGL()
	def mouseReleaseEvent(self, event):
		if self.current_mouse_mode:
			EMLightsDrawer.mouseReleaseEvent(self, event)
		else:
			for model in self.viewables:
				try:
					model.mouseReleaseEvent(event)
				except AttributeError, e:
					pass
		self.updateGL()
	def paintGL(self):
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT )
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		try:
			self.cam.position()
		except:
			return

		self.render()
	def render(self):
		try: 
			self.image_change_count = self.data["changecount"] # this is important when the user has more than one display instance of the same image, for instance in e2.py if 
		except:
			try: self.image_change_count = self.data[0]["changecount"]
			except: pass
			
		glPushMatrix()
		self.cam.position(True)
		# the ones are dummy variables atm... they don't do anything
		self.vdtools.update(1,1)
		glPopMatrix()
		
		dz = None
		if not self.perspective:
			glMatrixMode(GL_PROJECTION)
			glPushMatrix() 
			self.load_orthographic()
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
	def resizeEvent(self, event):
		self.vdtools.set_update_P_inv()
		EMGLWidget.resizeEvent(self, event)
	def resizeGL(self, width, height):
		# just use the whole window for rendering
		
		if width == 0 or height == 0: return # this is fine
		
		glViewport(0,0,width,height)
		
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
		
		self.set_projection_view_update()
		self.updateGL()
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
			self.last_window_height = height
	def rowChanged(self,row):
		if ( row == self.currentselection ): return
		self.currentselection=row
		self.updateGL()
	def set_cam_x(self,x):
		self.cam.set_cam_x( x )
		self.updateGL()	
	def set_cam_y(self,y):
		self.cam.set_cam_y( y )
		self.updateGL()		
	def set_cam_z(self,z):
		self.cam.set_cam_z( z )
		self.updateGL()
	def set_cam_z_from_fov_image(self,fov,image):
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
	def set_camera_defaults(self,data):
		if data:
			self.cam.default_z = -1.25*data.get_xsize()
			self.cam.cam_z = -1.25*data.get_xsize()
	def set_data(self,data,file_name="",replace=False):
		self.busy = True #Prevent updateGL() from doing anything until the end of this method
		was_previous_data = bool(self.data)
		
		###########TODO: update this code once we start doing mouse rotations and translations on the Widget camera
		previous_size = None
		previous_cam_data = {}
		previous_normalized_threshold = None #will be (threshold - mean)/(standard deviation)
		if was_previous_data:
			previous_size = (self.data["nx"], self.data["ny"], self.data["nz"])
			for model in self.viewables:
				if model.get_type() == "Isosurface":
					try: previous_normalized_threshold = (model.isothr - self.data["mean"])/self.data["sigma"]
					except: previous_normalized_threshold =1.0
					previous_cam_data = {"rot": model.cam.t3d_stack[-1], "pos":(model.cam.cam_x, model.cam.cam_y, model.cam.cam_z), "scale":model.cam.scale}
					break
		#############
		
		self.file_name = file_name # fixme fix this later
		if self.qt_parent != None:
			self.qt_parent.setWindowTitle(remove_directories_from_name(self.file_name))

		if data == None: 
			self.set_camera_defaults(data)
			return
		
		if data.get_zsize() == 1:
			new_data = EMData(data.get_xsize(),data.get_ysize(),2)
			new_data.insert_clip(data,[0,0,0])
			new_data.insert_clip(data,[0,0,1])
			self.data = new_data
		else: self.data = data
		
		#self.data.process_inplace("normalize.edgemean")
		
		nx,ny,nz = self.data.get_xsize(),self.data.get_ysize(),self.data.get_zsize()
		
		self.radius = nz
		if ny > nz: self.radius = ny
		if nx > ny: self.radius = nx
		self.radius /= 2
		#for i in self.viewables:
			#i.set_data(data)
		
		if self.inspector == None:
			self.inspector=EMImageInspector3D(self)
		
		if replace:
			self.inspector.delete_all()
			self.inspector.add_isosurface()
			###########TODO: update this code when mouse rotations and translations are performed on the Widget camera
			if previous_cam_data:
				if (nx, ny, nz) != previous_size:
					self.set_cam_z_from_fov_image(self.get_fov(), self.data)
				for model in self.viewables:
					model.busy = True
					inspector = model.get_inspector()
					inspector.update_rotations( previous_cam_data["rot"] )
#					inspector.rotation_sliders.slider_rotate()
					inspector.set_xyz_trans(*previous_cam_data["pos"])				
					if (nx,ny,nz) == previous_size:
						inspector.rotation_sliders.set_scale(previous_cam_data["scale"])
					if model.get_type() == "Isosurface":
						new_threshold = previous_normalized_threshold*data["sigma"] + data["mean"]
						inspector.thr.setValue(new_threshold) # Causes a GL update
						model.set_force_update(True)

					inspector.rotation_sliders.slider_rotate() #Causes a GL update
			#############
		else:
			if not self.viewables:
				self.inspector.add_isosurface()
			else:
				for model in self.viewables:
					model.set_data(self.data)
					inspector = model.get_inspector()
					if previous_cam_data:
						if model.get_type() == "Isosurface":
							new_threshold = previous_normalized_threshold*data["sigma"] + data["mean"]
							inspector.thr.setValue(new_threshold) # Causes a GL update
		
		if isinstance(self,EMImage3DWidget) and not was_previous_data:
			#self.resize(self.width(),self.height())
			self.set_cam_z_from_fov_image(self.get_fov(),self.data) #perhaps this should be done every time

		self.busy = False
		self.updateGL()
		for model in self.viewables:
			model.busy = False
			model.updateGL()
	def set_file_name(self,name):
		self.file_name = name
		if self.qt_parent != None:
			self.qt_parent.setWindowTitle(remove_directories_from_name(self.file_name))
	def set_perspective(self,bool):
		self.perspective = bool
		if self.emit_events: self.emit(QtCore.SIGNAL("set_perspective"),bool)
		self.updateGL()
		#self.set_perspective(bool)
	def set_scale(self,val):
		self.cam.scale = val
		self.updateGL()
		

	def viewing_volume(self):
		'''
		A useful function for getting the viewing volume as a list - so
		the calling function can figure out what to draw
		'''
		
		return [-self.xwidth/2.0,-self.yheight/2.0,self.startz, self.xwidth/2.0,self.yheight/2.0,self.endz]
	def wheelEvent(self, event):
		for model in self.viewables:
			try:
				model.wheelEvent(event)
			except AttributeError, e:
				pass
		self.updateGL()

class EMImageInspector3D(QtGui.QWidget):
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
		
		self.tabwidget = QtGui.QTabWidget()
		
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
		
		widget = self.tabwidget.widget(idx) #Ross: trying this to fix memory leak
		self.tabwidget.removeTab(idx)
		widget.deleteLater() #Ross: trying this to fix memory leak
		self.target().delete_current(self.targetidxmap[idx])

		
		self.target().updateGL()
	
	def delete_all(self):
		n = self.tabwidget.count()
		if n <= 1: return
		
		for idx in range(n-2,-1,-1):	
			widget = self.tabwidget.widget(idx) #Ross: trying this to fix memory leak
			self.tabwidget.removeTab(idx)
			widget.deleteLater() #Ross: trying this to fix memory leak
		
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
		self.persbut.setChecked(self.target().perspective==True)
		
		self.orthbut = QtGui.QRadioButton("Orthographic")
		self.orthbut.setChecked(self.target().perspective==False)
		
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


class EMImage3DModule(EMImage3DWidget):
	def __init__(self, parent=None, image=None,application=None,winid=None):
		import warnings	
		warnings.warn("convert EMImage3DModule to EMImage3DWidget", DeprecationWarning)
		EMImage3DWidget.__init__(self, parent, image, application, winid)
	
if __name__ == '__main__':
	from emapplication import EMApp
	import sys
	em_app = EMApp()
	window = EMImage3DWidget(application=em_app)

	if len(sys.argv)==1 : 
		data = []
		#for i in range(0,200):
		e = test_image_3d(1,size=(64,64,64))
		window.set_data(e)
	else :
		a=EMData(sys.argv[1])
		window.set_data(a,sys.argv[1])
	em_app.show()
	em_app.execute()
