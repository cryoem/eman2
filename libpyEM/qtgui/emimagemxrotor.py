#!/usr/bin/env python

#
# Author: David Woolford (sludtke@bcm.edu)
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

from math import ceil,tan,pi
import sys
import time
from weakref import WeakKeyDictionary

from EMAN2 import *
from emglobjects import *
from EMAN2db import EMAN2DB

from emfloatingwidgets import EMGLRotorWidget, EM2DGLView, EM3DGLWindowOverride, EM2DGLWindow
from emimagemx import EMImageInspectorMX, EMImageMXModule,EMDataListCache
from emimageutil import  EMEventRerouter
from emglobjects import EMOpenGLFlagsAndTools, EMGLProjectionViewMatrices
from emapplication import EMStandAloneApplication, EMGUIModule,get_application

import warnings
warnings.warn("emimagemxrotor.py is deprecated.", DeprecationWarning)

class EMImageMXRotorWidget(EMEventRerouter,QtOpenGL.QGLWidget,EMGLProjectionViewMatrices):
	"""
	"""
	allim=WeakKeyDictionary()
	def __init__(self, em_mx_rotor_module):
		
		assert(isinstance(em_mx_rotor_module,EMImageMXRotorModule))
		self.mmode = "drag"

		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True);
		fmt.setSampleBuffers(True)
		#fmt.setDepthBuffer(True)
		QtOpenGL.QGLWidget.__init__(self,fmt)
		EMEventRerouter.__init__(self,em_mx_rotor_module)
		EMGLProjectionViewMatrices.__init__(self)
		EMImageMXRotorWidget.allim[self]=0
		
		self.setFocusPolicy(Qt.StrongFocus)
		
		self.imagefilename = None
		
		self.fov = 30
		self.aspect = 1.0
		self.z_near = 6000
		self.z_far = 13000

		self.light_0_pos = [0.1,.1,1.,0.]
		
		self.resize(480,480)
		
	def set_data(self,data):
		self.target().set_data(data)
	
	def get_optimal_size(self):
		lr = self.target().rotor.get_suggested_lr_bt_nf()
		width = lr[1] - lr[0]
		height = lr[3] - lr[2]
		return [width+20,height+20]

	
	def set_image_file_name(self,name):
		#print "set image file name",name
		self.imagefilename = name
		self.target().set_image_file_name(name)
		
	def get_image_file_name(self):
		return self.imagefilename
	
	def initializeGL(self):
		glClearColor(0,0,0,0)
		
		glEnable(GL_LIGHTING)
		glEnable(GL_LIGHT0)
		glEnable(GL_DEPTH_TEST)
		glLightfv(GL_LIGHT0, GL_AMBIENT, [0.1, 0.1, 0.1, 1.0])
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_POSITION, [.1,.1,1.,1.])
		glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,GL_TRUE)

		
		glEnable(GL_DEPTH_TEST)
		
		glEnable(GL_NORMALIZE)
		
		glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST)
		glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST)
		glHint(GL_TEXTURE_COMPRESSION_HINT, GL_NICEST)
		
	def paintGL(self):
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		#glLoadIdentity()
		if ( self.target() == None ): return
		self.target().render()

	
	def resizeGL(self, width, height):
		if width <= 0 or height <= 0: return None
		GL.glViewport(0,0,width,height)
	
		GL.glMatrixMode(GL.GL_PROJECTION)
		GL.glLoadIdentity()
		self.aspect = float(width)/float(height)
		GLU.gluPerspective(self.fov,self.aspect,self.z_near,self.z_far)
		#GL.glOrtho(0.0,width,0.0,height,-width,width)
		GL.glMatrixMode(GL.GL_MODELVIEW)
		GL.glLoadIdentity()
		
		self.set_projection_view_update()
		self.target().resize_event(width,height)
	
	def set_near_far(self,near,far):
		self.z_near = near
		self.z_far = far
		
		GL.glMatrixMode(GL.GL_PROJECTION)
		GL.glLoadIdentity()
		GLU.gluPerspective(self.fov,self.aspect,self.z_near,self.z_far)
		#GL.glOrtho(0.0,width,0.0,height,-width,width)
		GL.glMatrixMode(GL.GL_MODELVIEW)
		GL.glLoadIdentity()
		
		self.target().projection_or_viewport_changed()
		
	def get_depth_for_height(self, height):
		depth = height/(2.0*tan(self.fov/2.0*pi/180.0))
	
		return depth
	
	def get_depth_for_width(self, width):
		equiv_height = width/self.aspect
		return self.get_depth_for_height(equiv_height)

	def set_mouse_mode(self,mode):
		self.mmode = mode
		self.target().set_mouse_mode(mode)
	
	def dropEvent(self,event):
		self.target().dropEvent(event)
	
	def set_shapes(self,shapes,shrink):
		self.target().set_shapes(shapes,shrink)
	
	def set_frozen(self,frozen):
		self.target().set_frozen(frozen)
	
	def get_frame_buffer(self):
		# THIS WILL FAIL ON WINDOWS APPARENTLY, because Windows requires a temporary context to be created and this is what the True flag
		# trying to stop.
		return None
	
		return self.renderPixmap(0,0,True)
		# to get around it we would have to render everything without display lists (a supreme pain).
		
	def set_selected(self,n):
		return self.target().set_selected(n)
	
	def get_core_object(self):
		return self.target()
	
	def keyPressEvent(self,event):
		self.target().keyPressEvent(event)
	
class EMImageMXRotorModule(EMGUIModule):
	
	def get_desktop_hint(self):
		return "rotor"
	
	
	def get_gl_widget(self,qt_context_parent,gl_context_parent):
		from emfloatingwidgets import EM3DGLWindowOverride
		if self.gl_widget == None:
			
			self.gl_context_parent = gl_context_parent
			self.qt_context_parent = qt_context_parent
			
			self.gl_widget = EM3DGLWindowOverride(self,self.rotor)
			self.gl_widget.set_enable_clip(False)
			self.widget = self.gl_widget
			self.gl_widget.resize(640,640)
			self.disable_mx_zoom()
			self.disable_mx_translate()
			
		return self.gl_widget
		
	def get_qt_widget(self):
		if self.qt_context_parent == None:	
			from emimageutil import EMParentWin
			self.gl_context_parent = EMImageMXRotorWidget(self)
			self.qt_context_parent = EMParentWin(self.gl_context_parent)
			self.gl_widget = self.gl_context_parent
		
			self.qt_context_parent.setAcceptDrops(True)

		return self.qt_context_parent
	
		
	def using_ftgl(self): return False
	
	def __init__(self, data=None,application=None):
		self.widget = None
		self.data=None
		self.rotor = EMGLRotorWidget(self,-15,70,-15,EMGLRotorWidget.BOTTOM_ROTARY,200)
		#self.rotor.set_angle_range(110.0)
		#self.rotor.set_child_mouse_events(False)
		self.rotor.set_mouse_mode("mxrotor")
		
		
		self.image_file_name = None	# keeps track of the image file name (if any) - book keeping purposes only
		self.emdata_list_cache = None # all import emdata list cache, the object that stores emdata objects efficiently. Must be initialized via set_data or set_image_file_name
		
		self.default_mxs = 3
		self.visible_mxs = self.default_mxs	# the number of visible imagemxs in the rotor
		self.mx_rows = 4 # the number of rows in any given imagemx
		self.mx_cols = 4 # the number of columns in any given imagemx
		self.start_mx = 0 # the starting index for the currently visible set of imagemxs
		
		EMGUIModule.__init__(self,ensure_gl_context=True)
		
		self.inspector = None
		self.minden=0
		self.maxden=1.0
		self.mindeng=0
		self.maxdeng=1.0
		self.hist = None
		self.gamma=1.0
		self.mmode = 'app'
		self.rot_mode = 'app'
		
		self.vals_to_display = ["Img #"]
		
		self.z_near = 0
		self.z_far = 0

		self.init_flag = True
		
		self.__init_font_renderer()
		
		if data:
			self.set_data(data)

	def __del__(self):
		for widget in self.rotor.get_widgets():
			get_application().deregister_qt_emitter(widget.get_drawable().get_drawable())
		
	def __init_gl_widget(self):
		self.widget = EM3DGLWindowOverride(self,self.rotor)
		self.widget.set_enable_clip(False)
		self.widget.set_draw_frame(False)
		self.disable_mx_zoom()
		self.disable_mx_translate()
			
	def __init_font_renderer(self):
		try:
			self.font_renderer = get_3d_font_renderer()
			self.font_renderer.set_face_size(32)
			self.font_renderer.set_depth(8)
#			self.font_renderer.set_font_mode(FTGLFontMode.EXTRUDE)
			self.font_render_mode = EMGUIModule.FTGL
		except:
			self.font_render_mode = EMGUIModule.GLUT
		
	def get_inspector(self):
		return self.inspector
	
	#def emit(self,*args,**kargs):
		#self.application.get_qt_emitter(self).emit(*args,**kargs)
		
	def set_selected(self,n,update_gl=True):
		widget = self.rotor[0]
		if widget != None:
			widget.get_drawable().get_drawable().set_selected(n,update_gl)
			
	def get_rows(self):
		return self.mx_rows
	
	def get_cols(self):
		return self.mx_cols
	
	def get_mxs(self):
		return self.visible_mxs

	def context(self):
		# asking for the OpenGL context from the parent
		return self.gl_context_parent.context()
	
	def is_visible(self,n):
		widget = self.rotor[0]
		if widget != None:
			img_offset = widget.get_drawable().get_drawable().get_img_num_offset()
			if n >= img_offset and n < (img_offset+self.mx_rows*self.mx_cols):
				return True
			else: return False
		else: return False
	
	def scroll_to(self,n,unused):
		widget = self.rotor[0]
		if widget != None:
			img_offset = widget.get_drawable().get_drawable().get_img_num_offset()
			scroll = (n-img_offset)/(self.mx_rows*self.mx_cols)
			self.rotor.explicit_animation(-scroll)
	
	def set_mx_cols(self,cols):
		self.mx_cols = cols
		self.emdata_list_cache.set_cache_size(8*self.visible_mxs*self.mx_rows*self.mx_cols)
		self.__refresh_rotor(True)
		self.updateGL()
		
	def set_mx_rows(self,rows):
		self.mx_rows = rows
		self.emdata_list_cache.set_cache_size(8*self.visible_mxs*self.mx_rows*self.mx_cols)
		self.__refresh_rotor(True)
		self.updateGL()
	
	def set_display_values(self,v2d):
		self.vals_to_display = v2d
		self.__refresh_rotor(False)
		self.updateGL()

	def set_mxs(self,mxs):
		self.visible_mxs = mxs
		self.emdata_list_cache.set_cache_size(8*self.visible_mxs*self.mx_rows*self.mx_cols)
		self.__regenerate_rotor() # this could be done more efficiently
		self.__refresh_rotor_size()
		self.updateGL()
	
	def set_mouse_mode(self,mode):
		self.mmode = mode
		try:
			for i in range(self.visible_mxs):
				w = self.rotor[i].get_drawable()
				w.set_mouse_mode(self.mmode)
		except: pass
	def set_scale(self,scale):
		pass
	
	def get_density_min(self):
		return self.rotor[0].get_drawable().get_drawable().get_density_min()
	
	def get_density_max(self):
		return self.rotor[0].get_drawable().get_drawable().get_density_max()
	
	def get_hist(self):
		return self.rotor[0].get_drawable().get_drawable().get_hist()
	
	def set_density_max(self,val):
		self.maxden=val
		self.update_min_max_gamma()
	
	def set_density_min(self,val):
		self.minden=val
		self.update_min_max_gamma()
		
	def set_gamma(self,val):
		self.gamma=val
		self.update_min_max_gamma()
	
	def set_plane(self,plane):
		self.plane = plane
	
	def disable_mx_zoom(self):
		'''
		Disable mx zoom.
		'''
		self.rotor.target_zoom_events_allowed(False)
		
	def disable_mx_translate(self):
		self.widget.target_translations_allowed(False)
	
	def allow_camera_rotations(self,bool=False):
		self.widget.allow_camera_rotations(bool)
	
	def pop_box_image(self,idx):
		val = self.emdata_list_cache.delete_box(idx)
		if val == 1:
			self.max_idx = self.emdata_list_cache.get_max_idx()
			self.__refresh_rotor()
		elif val == 2:
			w = self.rotor[0].get_drawable()
			w.force_display_update()
		else:
			print 'failed to delete box image'
		
	def update_min_max_gamma(self,update_gl=True):
		for i in range(self.visible_mxs):
			w = self.rotor[i].get_drawable().get_drawable()
			w.set_min_max_gamma(self.minden,self.maxden,self.gamma,False)
			if  i == 0 and  self.inspector != None:	
				try:
					self.inspector.set_hist(w.get_hist(),self.minden,self.maxden)
				except:
					# the histogram isn't created yet - this is a FIXME
					pass
		
		if update_gl: self.updateGL()
	
	def set_den_range(self,minden,maxden):
		self.minden=minden
		self.maxden=maxden
		self.update_min_max_gamma()
				
	def update_rotor_position(self,inc):
		self.start_mx += inc
	
		if inc > 0:
			end_changed = self.visible_mxs
			start_changed = self.visible_mxs-inc
			
			#print start_changed,end_changed
			#for idx in range(start_changed,end_changed):
			
		elif inc < 0:
			
			end_changed = -inc
			start_changed = 0
			
		else: return
		
		#print "iterating through",start_changed,end_changed,"start_mx is",self.start_mx
		
		self.__update_rotor_range(start_changed ,end_changed)
		
		try:
			w = self.rotor[0].get_drawable().get_drawable()
			self.inspector.set_hist(w.get_hist(),self.minden,self.maxden)
		except: pass
	
	def __update_rotor_range(self,start_changed,end_changed):
		
		num_per_view = self.mx_rows*self.mx_cols
		
		for idx in range(start_changed,end_changed):
			n = idx+ self.start_mx
			
			panel = n
			if panel != 0: panel %= self.get_num_panels()
			
			start_idx = panel*num_per_view
			image_offset = start_idx 
			if image_offset != 0: image_offset %= self.emdata_list_cache.get_max_idx()
			
			d = []
			if image_offset > (self.emdata_list_cache.get_max_idx()-num_per_view):
				num_visible = self.emdata_list_cache.get_max_idx() - image_offset
				num_none = num_per_view - num_visible
				for i in range(start_idx,start_idx+num_visible): d.append(self.emdata_list_cache[i])
				for i in range(num_none): d.append(None)
			else:
				for i in range(start_idx,start_idx+num_per_view): d.append(self.emdata_list_cache[i])
			

			w = self.rotor[idx].get_drawable().get_drawable()
			w.set_data(d,"",False)

			w.set_img_num_offset(image_offset)
			w.set_max_idx(self.emdata_list_cache.get_max_idx())
			w.set_min_max_gamma(self.minden,self.maxden,self.gamma,False)
	
	def get_num_panels(self):
		return int(ceil(float(self.emdata_list_cache.get_max_idx())/ (self.mx_rows*self.mx_cols)))

	def set_rotor_mode(self,mode):
		self.rot_mode = mode
		self.rotor.set_mouse_mode(mode)

	def optimize_fit(self,update_gl=True):
		render_width = self.gl_widget.width()
		render_height = self.gl_widget.height()
		try:
			self.mx_rows = render_width/self.emdata_list_cache.get_xsize()
			self.mx_cols = render_height/self.emdata_list_cache.get_ysize()
			if self.mx_rows == 0: self.mx_rows = 1
			if self.mx_cols == 0: self.mx_cols = 1
		except: return
		
		try:
			self.inspector.set_n_cols(self.mx_cols)
			self.inspector.set_n_rows(self.mx_rows)
		except: pass

		#try:
			#self.widget.set_width(render_width)
			#self.widget.set_height(render_height)
		#except: pass
		
		self.__refresh_rotor(True)
		if update_gl: self.updateGL()
		
	def set_frozen(self,frozen):
		self.rotor.set_frozen(frozen)

	def set_shapes(self,shapes,shrink):
		self.rotor.set_shapes(shapes,shrink)

	def register_animatable(self,animatable):
		self.qt_context_parent.register_animatable(animatable)
		
	def width(self):
		try: return self.gl_widget.width()
		except: return 0
	
	def height(self):
		try: return self.gl_widget.height()
		except: return 0
	
	def get_image_file_name(self):
		''' warning - could return none in some circumstances'''
		try: return self.gl_parent.get_image_file_name()
		except: return None
	
	def get_image(self,idx):
		return self.emdata_list_cache[idx]
	
	def set_data(self,data):
		if data == None or not isinstance(data,list) or len(data)==0:
			self.data = [] 
			return
		fac = int(ceil(float(len(data))/(self.mx_rows*self.mx_cols)))
		if fac == 0: fac = 1
		self.emdata_list_cache = EMDataListCache(data,8*self.visible_mxs*self.mx_rows*self.mx_cols)
		if self.init_flag:
			if fac < self.visible_mxs:
				self.visible_mxs = fac
			self.__regenerate_rotor()
			self.init_flag = False
		else:
#			print "here we are", fac,self.visible_mxs,self.default_mxs
			if fac < self.visible_mxs:
				self.visible_mxs = fac
				self.__regenerate_rotor()
				self.__refresh_rotor_size()
			elif self.visible_mxs < self.default_mxs and fac >= self.default_mxs:
#				print "setting visible mxs"
				self.visible_mxs = self.default_mxs
				self.__regenerate_rotor()
				self.__refresh_rotor_size()
			elif fac > self.visible_mxs and self.visible_mxs < self.default_mxs:
				self.visible_mxs = fac
				self.__regenerate_rotor()
				self.__refresh_rotor_size()
			else:
				self.__refresh_rotor(True)
			
	
	def set_image_file_name(self,name):
		#print "set image file name",name
		self.image_file_name = name
		self.emdata_list_cache = EMDataListCache(name)
		self.__regenerate_rotor()
	
	def __refresh_rotor(self,inc_size=False):
			
		self.__update_rotor_range(0,self.visible_mxs)
		
		if inc_size:
			self.__refresh_rotor_size()

	def __refresh_rotor_size(self):
		if len(self.rotor) == 0: return
		
		
		
		self.render_width = self.gl_widget.width()
		self.render_height = self.gl_widget.height()
		width = self.mx_rows*(self.emdata_list_cache.get_xsize()+2)-2
		height = self.mx_cols*(self.emdata_list_cache.get_ysize()+2)-2
		scale1 = self.render_height/float(height)
		scale2 = self.render_width/float(width)
		
		if self.render_width > self.render_height:
			self.render_height = float(height)/width*self.render_width
			scale = scale2
		else:
			self.render_width = float(width)/height*self.render_height
			scale = scale1
		
		for idx in range(0,self.visible_mxs):
			e = self.rotor[idx]
			w = e.get_drawable().get_drawable()

			w.set_scale(scale,False)
			e.set_width(self.render_width)
			e.set_height(self.render_height)
			e.set_update_frame(True)
			##w.set_mx_cols(self.mx_rows,False)
			
		self.rotor.update()

	def __regenerate_rotor(self):
		for widget in self.rotor.get_widgets():
			get_application().deregister_qt_emitter(widget.get_drawable().get_drawable())
			self.rotor.clear_widgets()
#		self.parent.updateGL() # i can't figure out why I have to do this (when this function is called from set_mxs
		num_per_view = self.mx_rows*self.mx_cols
		for idx in range(self.start_mx,self.start_mx+self.visible_mxs):
			n = idx
			
			panel = n
			if panel != 0: panel %= self.get_num_panels()
			
			start_idx = panel*num_per_view
		
			image_offset = start_idx 
			if image_offset != 0: image_offset %= self.emdata_list_cache.get_max_idx()
			#print idx,panel,image_offset
			
			d = []
			if image_offset > (self.emdata_list_cache.get_max_idx()-num_per_view):
				num_visible = self.emdata_list_cache.get_max_idx() - image_offset
				num_none = num_per_view - num_visible
				#print num_visible,num_none
				for i in range(start_idx,start_idx+num_visible): d.append(self.emdata_list_cache[i])
				for i in range(num_none): d.append(None)
			else:
				for i in range(start_idx,start_idx+num_per_view): d.append(self.emdata_list_cache[i])

			e = EM2DGLView(self,d)
			e.get_drawable().set_app(get_application())
			get_application().register_qt_emitter(e.get_drawable(),get_application().get_qt_emitter(self))
			x = EM2DGLWindow(self,e)
			x.decoration.draw_x_enabled = False
			x.decoration.color_flag = "black"
			self.rotor.add_widget(x)
			
			w = e.get_drawable()
	
			w.set_use_display_list(True) # saves HEAPS of time, makes interaction much smoother
			w.set_reroute_delete_target(self)
			w.set_draw_background(True)
			w.set_mouse_mode(self.mmode)
			w.set_display_values(self.vals_to_display,False)	

			w.set_img_num_offset(image_offset)
			w.set_max_idx(self.emdata_list_cache.get_max_idx())

		e = self.rotor[0]
		w = e.get_drawable().get_drawable()
		self.minden = w.get_density_min()
		self.maxden = w.get_density_max()
		self.hist = w.get_hist()
		self.mindeng = self.minden
		self.maxdeng = self.maxden
		self.gamma = w.get_gamma()
		self.update_min_max_gamma(False)
		
		self.__refresh_rotor_size()

	def updateGL(self):
		try: self.gl_widget.updateGL()
		except: pass


	def render(self):
		if self.widget == None: self.__init_gl_widget()
		
		if not self.get_qt_context_parent().isVisible(): return
		if self.emdata_list_cache == None: return
		
		
		lrt = self.widget.get_lr_bt_nf()

		
		GL.glEnable(GL.GL_DEPTH_TEST)
		GL.glEnable(GL.GL_LIGHTING)
		z = self.gl_context_parent.get_depth_for_height(abs(lrt[3]-lrt[2]))
		z2 = self.gl_context_parent.get_depth_for_width(abs(lrt[1]-lrt[0]))
		
		if z2 > z: z = z2
		
		
		z_near = z-lrt[4]-100
		z_trans = 0
		z_far = z-lrt[5]
		if z_near < 0:
			z_trans = z_near
			z_near = 1
			z_far -= z_trans
		if z_far < 0: z_far = 0.1 # hacking alert
		z_far += abs(lrt[3]-lrt[2]) # hacking alert
		if self.z_near != z_near or self.z_far != z_far:
			self.z_near = z_near
			self.z_far = z_far
			if isinstance(self.gl_context_parent,EMImageMXRotorWidget):
				self.gl_context_parent.set_near_far(self.z_near,self.z_far)
			else: print "bug 2"

		GL.glPushMatrix()
		glTranslate(-(lrt[1]+lrt[0])/2.0,-(lrt[3]+lrt[2])/2.0,-z)
		self.widget.draw()
		GL.glPopMatrix()
	
		self.draw_hud()
	
	def dragEnterEvent(self,event):
		pass

	def dropEvent(self,event):
		pass

	def get_inspector(self):
		if not self.inspector: 
			self.inspector=EMImageInspectorMX(self,allow_col_variation=True,allow_window_variation=True,allow_opt_button=True)
			self.inspector.set_limits(self.mindeng,self.maxdeng,self.minden,self.maxden)
			
		return self.inspector
	def mousePressEvent(self, event):
		if event.button()==Qt.MidButton or (event.button()==Qt.LeftButton and event.modifiers()&Qt.AltModifier):
			self.show_inspector(True)
		else:
			self.widget.mousePressEvent(event)
			
		self.updateGL()
		
	def mouseMoveEvent(self, event):
		if self.widget.isinwin(event.x(),self.gl_context_parent.viewport_height()-event.y()):
			self.widget.mouseMoveEvent(event)
			self.updateGL()
		
	def mouseReleaseEvent(self, event):
		self.widget.mouseReleaseEvent(event)
		self.updateGL()
		
	def wheelEvent(self, event):
		if self.widget.isinwin(event.x(),self.gl_context_parent.viewport_height()-event.y()):
			self.widget.wheelEvent(event)
			self.updateGL()
		
	def leaveEvent(self,event):
		pass
	
	def keyPressEvent(self,event):
		self.widget.keyPressEvent(event)
		self.updateGL()
		
	def resize_event(self, width, height):
		if self.widget != None: 
			self.widget.update()
		self.rotor.resize_event(width,height)
		self.optimize_fit(False)
	
	def projection_or_viewport_changed(self):
		self.rotor.resize_event(-1,1)
	
	def save_lst(self):
		self.emdata_list_cache.save_lst()
		
	def save_data(self):
		self.emdata_list_cache.save_data()
	
	def get_frame_buffer(self):
		return self.gl_widget.get_frame_buffer()
	
	
	ligh_yellow_diffuse = (.84,.82,.38,1.0)
	ligh_yellow_ambient = (.83,.83,.38,1.0)
	ligh_yellow_specular = (.76,.75,.39,1.0)
	def draw_hud(self):
		width = self.gl_widget.viewport_width()
		height =self.gl_widget.viewport_height()
		glMatrixMode(GL_PROJECTION)
		glPushMatrix()
		glLoadIdentity()
		glOrtho(0,width,0,height,-200,200)
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		glDisable(GL_LIGHTING)
		glEnable(GL_NORMALIZE)

		glMaterial(GL_FRONT,GL_AMBIENT,EMImageMXRotorModule.ligh_yellow_ambient)
		glMaterial(GL_FRONT,GL_DIFFUSE,EMImageMXRotorModule.ligh_yellow_diffuse)
		glMaterial(GL_FRONT,GL_SPECULAR,EMImageMXRotorModule.ligh_yellow_specular)
		glMaterial(GL_FRONT,GL_SHININESS,30.0)
		
		enable_depth = glIsEnabled(GL_DEPTH_TEST)
		glDisable(GL_DEPTH_TEST)
		glEnable(GL_TEXTURE_2D)
		glColor(1.0,1.0,1.0)
		glColor(*EMImageMXRotorModule.ligh_yellow_specular)
		
		if self.font_render_mode == EMImageMXModule.FTGL:
			panels = self.get_num_panels()
			idx = self.start_mx
			if idx != 0: idx %= panels
			string = str(idx+1) + ' / ' + str(panels)
			glPushMatrix()
			glTranslate(10,height-1.2*self.font_renderer.get_face_size(),0)
			glRotate(20,0,1,0)
			self.font_renderer.render_string(string)
			glPopMatrix()
		else:
			pass
		
		if enable_depth: glEnable(GL_DEPTH_TEST)
		
		glMatrixMode(GL_PROJECTION)
		glPopMatrix()
		glMatrixMode(GL_MODELVIEW)
		
	
	
# This is just for testing, of course
if __name__ == '__main__':
	em_app = EMStandAloneApplication()
	window = EMImageMXRotorModule(application=em_app)
	if len(sys.argv)==1 : 
		data = []
		for i in range(0,200):
			e = test_image(Util.get_irand(0,9))
			if ( Util.get_irand(0,4) == 0):	e.set_attr("excluded",True)
			data.append(e)
			
		window.set_data(data) 
	else :
#		a=EMData.read_images(sys.argv[1])
		window.set_image_file_name(sys.argv[1])
#		window.set_file_name(sys.argv[1])
#		window.set_data(a)
	
	em_app.show()
	#window.get_qt_widget().resize(640,640)
	#window.optimize_fit()
	em_app.execute()
