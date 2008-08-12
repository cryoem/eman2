#!/usr/bin/env python

#
# Author: Steven Ludtke (sludtke@bcm.edu)
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
from emimageutil import EMParentWin
from emglobjects import EMOpenGLFlagsAndTools
from emfloatingwidgets import EMGLRotorWidget, EMGLView2D, EM3DWidget
from emimagemx import EMImageMxInspector2D, EMImageMXCore
from EMAN2db import EMAN2DB

class EMImageMXRotor(QtOpenGL.QGLWidget):
	"""
	"""
	allim=WeakKeyDictionary()
	def __init__(self, data=None,parent=None):
		self.image_rotor = None
		#self.initflag = True
		self.mmode = "drag"

		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True);
		fmt.setSampleBuffers(True)
		#fmt.setDepthBuffer(True)
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		EMImageMXRotor.allim[self]=0
		
		self.setFocusPolicy(Qt.StrongFocus)
		
		self.image_rotor = EMImageMXRotorCore(data,self)
		
		self.imagefilename = None
		
		self.fov = 30
		self.aspect = 1.0
		self.z_near = 6000
		self.z_far = 13000
		
		self.animatables = []
		
		self.timer = QtCore.QTimer()
		QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.timeout)
		self.timer.start(10)
		
		self.light_0_pos = [0.1,.1,1.,0.]
		
	def setData(self,data):
		self.image_rotor.setData(data)
	
	def get_optimal_size(self):
		lr = self.image_rotor.rotor.get_suggested_lr_bt_nf()
		width = lr[1] - lr[0]
		height = lr[3] - lr[2]
		return [width+80,height+20]
	
	def timeout(self):
		
		if len(self.animatables) == 0: 
			#self.light_0_pos[0] +=  sin(Util.get_frand(-0.2,0.2))
			#self.light_0_pos[1] +=  sin(Util.get_frand(-0.2,0.2))
			#self.light_0_pos[2] +=  sin(Util.get_frand(-0.2,0.2))
			#glLightfv(GL_LIGHT0, GL_POSITION, self.light_0_pos)
			#self.updateGL()
			return
		for i,animatable in enumerate(self.animatables):
			try:
				if not animatable.animate(time.time()):
					# this could be dangerous
					self.animatables.pop(i)
			
			except:
				self.animatables.pop(i)
		self.updateGL()
		
	def register_animatable(self,animatable):
		self.animatables.append(animatable)
	
	def set_image_file_name(self,name):
		#print "set image file name",name
		self.imagefilename = name
		self.image_rotor.set_image_file_name(name)
		
	def get_image_file_name(self):
		return self.imagefilename
	
	def initializeGL(self):
		glClearColor(0,0,0,0)
		
		glLightfv(GL_LIGHT0, GL_AMBIENT, [0.1, 0.1, 0.1, 1.0])
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_POSITION, self.light_0_pos)
	
		glEnable(GL_LIGHTING)
		glEnable(GL_LIGHT0)
		
		glEnable(GL_DEPTH_TEST)
		
		glEnable(GL_NORMALIZE)
		
		glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST)
		glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST)
		glHint(GL_TEXTURE_COMPRESSION_HINT, GL_NICEST)
		
	def paintGL(self):
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		#glLoadIdentity()
		if ( self.image_rotor == None ): return
		self.image_rotor.render()

	
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
		
		self.image_rotor.resize_event(width,height)
	
	def set_near_far(self,near,far):
		self.z_near = near
		self.z_far = far
		
		GL.glMatrixMode(GL.GL_PROJECTION)
		GL.glLoadIdentity()
		GLU.gluPerspective(self.fov,self.aspect,self.z_near,self.z_far)
		#GL.glOrtho(0.0,width,0.0,height,-width,width)
		GL.glMatrixMode(GL.GL_MODELVIEW)
		GL.glLoadIdentity()
		
	
		self.image_rotor.projection_or_viewport_changed()
		
	def get_depth_for_height(self, height):
		# This function returns the width and height of the renderable 
		# area at the origin of the data volume
		depth = height/(2.0*tan(self.fov/2.0*pi/180.0))
	
		return depth
	
	def set_mmode(self,mode):
		self.mmode = mode
		self.image_rotor.set_mmode(mode)
	
	def mousePressEvent(self, event):
		self.image_rotor.mousePressEvent(event)
		self.updateGL()
		
	def wheelEvent(self,event):
		self.image_rotor.wheelEvent(event)
		self.updateGL()
		
	def mouseMoveEvent(self,event):
		self.image_rotor.mouseMoveEvent(event)
		self.updateGL()
		
	def mouseReleaseEvent(self,event):
		self.image_rotor.mouseReleaseEvent(event)
		self.updateGL()
		
	def closeEvent(self,event) :
		self.image_rotor.closeEvent(event)
		
	def dragEnterEvent(self,event):
		self.image_rotor.dragEnterEvent(event)

	def dropEvent(self,event):
		self.image_rotor.dropEvent(event)
	
	def set_shapes(self,shapes,shrink):
		self.image_rotor.set_shapes(shapes,shrink)
	
	def set_frozen(self,frozen):
		self.image_rotor.set_frozen(frozen)
	
	def get_frame_buffer(self):
		# THIS WILL FAIL ON WINDOWS APPARENTLY, because Windows requires a temporary context to be created and this is what the True flag
		# trying to stop.
		return None
	
		return self.renderPixmap(0,0,True)
		# to get around it we would have to render everything without display lists (a supreme pain).
		
	def set_selected(self,n):
		return self.image_rotor.set_selected(n)
	
	def get_core_object(self):
		return self.image_rotor
	
	def keyPressEvent(self,event):
		self.image_rotor.keyPressEvent(event)
	
class EMImageMXRotorCore:
	def __init__(self, data=None,parent=None):
		self.parent = parent
		self.data=None
		try: self.parent.setAcceptDrops(True)
		except:	pass

		self.initsizeflag = True
	
		
		self.rotor = EMGLRotorWidget(self,0,-70,-15,EMGLRotorWidget.TOP_ROTARY,100)
		self.rotor.set_angle_range(110.0)
		#self.rotor.set_child_mouse_events(False)
		self.rotor.set_mmode("mxrotor")
		self.widget = EM3DWidget(self,self.rotor)
		self.widget.set_draw_frame(False)
		
		self.image_file_name = None	# keeps track of the image file name (if any) - book keeping purposes only
		self.emdata_list_cache = None # all import emdata list cache, the object that stores emdata objects efficiently. Must be initialized via setData or set_image_file_name
		
		self.default_mxs = 3
		self.visible_mxs = self.default_mxs	# the number of visible imagemxs in the rotor
		self.mx_rows = 4 # the number of rows in any given imagemx
		self.mx_cols = 4 # the number of columns in any given imagemx
		self.start_mx = 0 # the starting index for the currently visible set of imagemxs
		
		self.inspector = None
		self.minden=0
		self.maxden=1.0
		self.mindeng=0
		self.maxdeng=1.0
		self.gamma=1.0
		self.mmode = 'app'
		self.rot_mode = 'app'
		
		self.display_help_hud = False
		self.display_help = ["Wheel - Rotor Rotate","Ctrl+Wheel - Zoom", "Ctrl+Right Click - Move Rotor","Shift+Left Click - Delete Box"] 
		
		self.vals_to_display = ["Img #"]
		
		self.z_near = 0
		self.z_far = 0

		self.init_flag = True
		try:
			self.font_renderer = get_3d_font_renderer()
			self.font_renderer.set_face_size(32)
			self.font_renderer.set_depth(8)
			self.font_renderer.set_font_mode(FTGLFontMode.EXTRUDE)
			self.font_render_mode = EMImageMXCore.FTGL
		except:
			self.font_render_mode = EMImageMXCore.GLUT
	
		self.disable_mx_zoom()
		self.disable_mx_translate()
		if data:
			self.setData(data)
	def get_inspector(self):
		return self.inspector
	
	def emit(self,signal,event,data=None,bool=None):
		if bool==None:
			self.parent.emit(signal,event,data)
		elif data == None:
			self.parent.emit(signal,event)
		else:
			self.parent.emit(signal,event,data,bool)
	
	def set_selected(self,n,update_gl=True):
		self.rotor[0].get_drawable().set_selected(n,update_gl)
	
	
	def get_qt_parent(self):
		return self.parent
	
	def get_rows(self):
		return self.mx_rows
	
	def get_cols(self):
		return self.mx_cols
	
	def get_mxs(self):
		return self.visible_mxs

	def context(self):
		# asking for the OpenGL context from the parent
		return self.parent.context()
	
	def is_visible(self,n):
		img_offset = self.rotor[0].get_drawable().get_img_num_offset()
		if n >= img_offset and n < (img_offset+self.mx_rows*self.mx_cols):
			return True
		else: return False
	
	def scroll_to(self,n,unused):
		img_offset = self.rotor[0].get_drawable().get_img_num_offset()
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
	
	def set_mmode(self,mode):
		self.mmode = mode
		try:
			for i in range(self.visible_mxs):
				w = self.rotor[i].get_drawable()
				w.set_mmode(self.mmode)
		except: pass
	def set_scale(self,scale):
		pass
	
	def set_density_max(self,val):
		self.maxden=val
		self.update_min_max_gamma()
	
	def set_density_min(self,val):
		self.minden=val
		self.update_min_max_gamma()
		
	def set_gamma(self,val):
		self.gamma=val
		self.update_min_max_gamma()
	
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
			w.force_dl_update()
		else:
			print 'failed to delete box image'
		
	def update_min_max_gamma(self,update_gl=True):
		for i in range(self.visible_mxs):
			w = self.rotor[i].get_drawable()
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
		
		w = self.rotor[idx].get_drawable()
		self.inspector.set_hist(w.get_hist(),self.minden,self.maxden)
	
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
			

			w = self.rotor[idx].get_drawable()
			w.setData(d,False)

			w.set_img_num_offset(image_offset)
			w.set_max_idx(self.emdata_list_cache.get_max_idx())
			w.set_min_max_gamma(self.minden,self.maxden,self.gamma,False)
	
	def get_num_panels(self):
		return int(ceil(float(self.emdata_list_cache.get_max_idx())/ (self.mx_rows*self.mx_cols)))

	def set_rotor_mode(self,mode):
		self.rot_mode = mode
		self.rotor.set_mmode(mode)

	def optimize_fit(self,update_gl=True):
		render_width = self.parent.width()
		render_height = self.parent.height()
		try:
			self.mx_rows = render_width/self.emdata_list_cache.get_image_width()
			self.mx_cols = render_height/self.emdata_list_cache.get_image_height()
			if self.mx_rows == 0: self.mx_rows = 1
			if self.mx_cols == 0: self.mx_cols = 1
		except: return
		
		try:
			self.inspector.set_n_cols(self.mx_cols)
			self.inspector.set_n_rows(self.mx_rows)
		except: pass

		self.__refresh_rotor(True)
		if update_gl: self.updateGL()
		
	def set_frozen(self,frozen):
		self.rotor.set_frozen(frozen)

	def set_shapes(self,shapes,shrink):
		self.rotor.set_shapes(shapes,shrink)

	def register_animatable(self,animatable):
		self.parent.register_animatable(animatable)
		
	def width(self):
		return self.parent.width()
	
	def height(self):
		return self.parent.height()

	def viewport_width(self):
		return self.parent.width()
	
	def viewport_height(self):
		return self.parent.height()
	
	def get_image_file_name(self):
		''' warning - could return none in some circumstances'''
		try: return self.parent.get_image_file_name()
		except: return None
	
	def get_image(self,idx):
		return self.emdata_list_cache[idx]
	
	def setData(self,data):
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
		
		self.render_width = self.parent.width()
		self.render_height = self.parent.height()
		width = self.mx_rows*(self.emdata_list_cache.get_image_width()+2)-2
		height = self.mx_cols*(self.emdata_list_cache.get_image_height()+2)-2
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
			w = e.get_drawable()

			w.set_scale(scale,False,False)
			e.set_width(self.render_width,False)
			e.set_height(self.render_height,False)
			w.set_mx_cols(self.mx_rows,False)
			
		self.rotor.update()

	def __regenerate_rotor(self):
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

			e = EMGLView2D(self,d)
			self.rotor.add_widget(e)
			w = e.get_drawable()
	
			w.set_use_display_list(True) # saves HEAPS of time, makes interaction much smoother
			w.set_reroute_delete_target(self)
			w.set_draw_background(True)
			w.set_mmode(self.mmode)
			w.set_display_values(self.vals_to_display,False)	

			w.set_img_num_offset(image_offset)
			w.set_max_idx(self.emdata_list_cache.get_max_idx())

		e = self.rotor[0]
		w = e.get_drawable()
		self.minden = w.get_density_min()
		self.maxden = w.get_density_max()
		self.mindeng = self.minden
		self.maxdeng = self.maxden
		self.gamma = w.get_gamma()
		self.update_min_max_gamma(False)
		
		self.__refresh_rotor_size()

	def updateGL(self):
		try: self.parent.updateGL()
		except: pass


	def render(self):
		if not self.parent.isVisible(): return
		if self.emdata_list_cache == None: return
		
		lr = self.rotor.get_suggested_lr_bt_nf()
		
		GL.glEnable(GL.GL_DEPTH_TEST)
		GL.glEnable(GL.GL_LIGHTING)
		z = self.parent.get_depth_for_height(abs(lr[3]-lr[2]))
		lrt = self.widget.get_lr_bt_nf()
		
		z_near = z-lrt[4]
		z_trans = 0
		z_far = z-lrt[5]
		if z_near < 0:
			z_trans = z_near
			z_near = 1
			z_far -= z_trans
		if z_far < 0: z_far = 0.1 # hacking alert
		z_far += abs(lr[3]-lr[2]) # hacking alert
		if self.z_near != z_near or self.z_far != z_far:
			self.z_near = z_near
			self.z_far = z_far
			self.parent.set_near_far(self.z_near,self.z_far)

		GL.glPushMatrix()
		#print -self.parent.get_depth_for_height(abs(lr[3]-lr[2])),self.z_near,self.z_far,abs(lr[3]-lr[2])
		glTranslate(-(lr[1]+lr[0])/2.0,-(lr[3]+lr[2])/2.0,-self.parent.get_depth_for_height(abs(lr[3]-lr[2]))+z_trans+abs(lr[3]-lr[2]))
		self.widget.paintGL()
		GL.glPopMatrix()
		
		self.draw_hud()
	
	def dragEnterEvent(self,event):
		pass

	def dropEvent(self,event):
		pass

	def show_inspector(self,force=False):
		if not force and self.inspector==None : return
		self.init_inspector()
		self.inspector.show()

	def init_inspector(self):
		allow_col_variation = True
		allow_mx_variation = True
		opt_fit_button_on = True
		if not self.inspector : self.inspector=EMImageMxInspector2D(self,allow_col_variation,allow_mx_variation,opt_fit_button_on)
		self.inspector.set_limits(self.mindeng,self.maxdeng,self.minden,self.maxden)

	def mousePressEvent(self, event):
		if event.button()==Qt.MidButton:
			self.show_inspector(True)
			self.emit(QtCore.SIGNAL("inspector_shown"),event)
		else:
			self.widget.mousePressEvent(event)
			
	
	def mouseMoveEvent(self, event):
		self.widget.mouseMoveEvent(event)
		
	def mouseReleaseEvent(self, event):
		self.widget.mouseReleaseEvent(event)
		
	def wheelEvent(self, event):
		self.widget.wheelEvent(event)
		
	def leaveEvent(self):
		pass
	
	def keyPressEvent(self,event):
		if event.key() == Qt.Key_F1:
			self.display_help_hud = not self.display_help_hud
			self.updateGL()

	def resize_event(self, width, height):
		self.rotor.resize_event(width,height)
		self.optimize_fit(False)
	
	def projection_or_viewport_changed(self):
		self.rotor.resize_event(-1,1)
	
	def save_lst(self):
		self.emdata_list_cache.save_lst()
		
	def save_data(self):
		self.emdata_list_cache.save_data()
	
	def get_frame_buffer(self):
		return self.parent.get_frame_buffer()
	
	def draw_hud(self):
		width = self.parent.width()
		height = self.parent.height()
		glMatrixMode(GL_PROJECTION)
		glPushMatrix()
		glLoadIdentity()
		glOrtho(0,width,0,height,-200,200)
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		glEnable(GL_LIGHTING)
		glEnable(GL_NORMALIZE)
		glMaterial(GL_FRONT,GL_AMBIENT,(0.2, 0.9, 0.2,1.0))
		glMaterial(GL_FRONT,GL_DIFFUSE,(0.2, 0.9, 0.9,1.0))
		glMaterial(GL_FRONT,GL_SPECULAR,(0.9, 0.5, 0.2,1.0))
		glMaterial(GL_FRONT,GL_SHININESS,20.0)
		
		glDisable(GL_DEPTH_TEST)
		glColor(1.0,1.0,1.0)
		
		
		if self.font_render_mode == EMImageMXCore.FTGL:
			panels = self.get_num_panels()
			idx = self.start_mx
			if idx != 0: idx %= panels
			string = str(idx+1) + ' / ' + str(panels)
			glPushMatrix()
			glTranslate(10,height-1.2*self.font_renderer.get_face_size(),0)
			glRotate(20,0,1,0)
			self.font_renderer.render_string(string)
			glPopMatrix()
			if self.display_help_hud:
				for i,s in enumerate(self.display_help):
					glPushMatrix()
					glTranslate(10,height-(i+2)*1.2*self.font_renderer.get_face_size(),0)
					glRotate(20,0,1,0)
					self.font_renderer.render_string(s)
					glPopMatrix()
					#string = str(s)
					#bbox = self.font_renderer.bounding_box(string)
					#x_offset = width-(bbox[3]-bbox[0]) - 10
					#y_offset += 10
					#glPushMatrix()
					#glTranslate(x_offset,y_offset,0)
					#glRotate(20,0,1,0)
					#self.font_renderer.render_string(string)
					#glPopMatrix()
					#y_offset += bbox[4]-bbox[1]
		else:
			pass
		
		glMatrixMode(GL_PROJECTION)
		glPopMatrix()
		glMatrixMode(GL_MODELVIEW)
		
			
class EMDataListCache:
	'''
	This class designed primarily for memory management in the context of large lists of EMData objects.
	
	The main public interface is to acquired a list of EMData objects in a specific range
	'''
	LIST_MODE = 'list_mode'
	FILE_MODE = 'file_mode'
	def __init__(self,object,cache_size=256,start_idx=0):
		DB = EMAN2db.EMAN2DB.open_db(".")
		if isinstance(object,list):
			# in list mode there is no real caching
			self.mode = EMDataListCache.LIST_MODE
			self.max_idx = len(object)
			self.cache_size = self.max_idx
			self.images = object
			self.start_idx = 0
			DB.open_dict("emimage_mx_rotor_cache")
			self.db = DB.emimage_mx_rotor_cache
			self.exclusions_key = "interactive_exclusions"
			self.db[self.exclusions_key] = []
			self.exclusions = self.db["interactive_exclusions"]
			
			for i,d in enumerate(self.images):	d.set_attr("original_number",i)

		elif isinstance(object,str):
			#print "file mode"
			self.mode = EMDataListCache.FILE_MODE
			if not os.path.exists(object):
				print "error, the file you specified does not exist:",object
				return
			self.file_name = object
			self.max_idx = EMUtil.get_image_count(self.file_name)
			self.images = {}
			if self.max_idx < cache_size:
				self.cache_size = self.max_idx
			else:
				self.cache_size = cache_size
			self.start_idx = start_idx - self.cache_size/2
			
			DB.open_dict("emimage_mx_rotor_cache")
			self.db = DB.emimage_mx_rotor_cache
			self.exclusions_key = self.file_name+"_interactive_exclusions"
			self.exclusions  = self.db[self.exclusions_key]
			if self.exclusions == None: self.exclusions = []
			self.__refresh_cache()
		else:
			print "the object used to construct the EMDataListCache is not a string (filename) or a list (of EMData objects). Can't proceed"
			return
		
		self.soft_delete = False # toggle to prevent permanent deletion of particles
		self.image_width = -1
		self.image_height = -1
	
	def __del__(self):
		DB = EMAN2db.EMAN2DB.open_db(".")
		try:
			DB.close_dict("emimage_mx_rotor_cache")
		except: pass
	
	def get_image_width(self):
		if self.mode == EMDataListCache.FILE_MODE:
			for i in self.images:
				try:
					if self.images[i] != None:
						return self.images[i].get_xsize()
				except: pass
				
			return 0
		elif self.mode == EMDataListCache.LIST_MODE:
			for i in self.images:
				try: return i.get_xsize()
				except: pass
				
				
			return 0
		
	def get_image_height(self):
		if self.mode == EMDataListCache.FILE_MODE:
			for i in self.images:
				try:
					if self.images[i] != None:
						return self.images[i].get_ysize()
				except: pass
				
			return 0
		elif self.mode == EMDataListCache.LIST_MODE:
			for i in self.images:
				try: return i.get_ysize()
				except: pass
				
				
			return 0
		
	def delete_box(self,idx):
		if self.mode == EMDataListCache.LIST_MODE and not self.soft_delete:
			# we can actually delete the emdata object
			image = self.images.pop(idx)
			self.max_idx = len(self.images)
			self.cache_size = self.max_idx
			return 1
		elif self.mode == EMDataListCache.FILE_MODE or self.soft_delete:
			im = self.images[idx]
			try:
				val = im.get_attr("excluded")
				if val == True:
					im.del_attr("excluded")
					self.exclusions.remove(idx)
				else: raise
			except: 
					im.set_attr("excluded",True)
					self.exclusions.append(idx)
			
			if self.mode == EMDataListCache.FILE_MODE: self.db[self.exclusions_key] = self.exclusions
			
			if self.mode == EMDataListCache.FILE_MODE: im.write_image(self.file_name,idx)
			return 2
		
		return 0
				
	def save_lst(self):
		# Get the output filespec
		fsp=QtGui.QFileDialog.getSaveFileName(None, "Specify lst file to save","","","")
		fsp=str(fsp)
		
		if fsp != '':
			f = file(fsp,'w')
			f.write('#LST\n')
			
			if self.mode == EMDataListCache.LIST_MODE and not self.soft_delete:
				for d in self.data:
					f.write(str(d.get_attr('original_number'))+'\n')
			elif self.mode == EMDataListCache.FILE_MODE or self.soft_delete:
				indices = [i for i in range(self.max_idx)]
				for exc in self.exclusions: indices.remove(exc)
				
				if self.mode ==  EMDataListCache.FILE_MODE:
					for idx in indices:	f.write(str(idx)+'\t'+self.file_name+'\n')
				elif self.soft_delete:
					for idx in indices:	f.write(str(idx)+'\n')
		
			f.close()
	def save_data(self):
		
		fsp=QtGui.QFileDialog.getSaveFileName(None, "Specify name of file to save","","","")
		fsp=str(fsp)
		
		try:
			if fsp == self.file_name:
				print "writing over the same file is currently not supported"
				return
		except: pass
		
		if fsp != '':
			for idx in range(self.max_idx):
				d = self.__getitem__(idx)
				try:
					d.get_attr("excluded")
					continue
				except: pass
				d.write_image(fsp,-1)

	def get_max_idx(self):
		return self.max_idx
	
	def get_num_images(self):
		return len(self.images)
	
	def set_cache_size(self,cache_size,refresh=False):
		if self.mode != EMDataListCache.LIST_MODE:
			if cache_size > self.max_idx: self.cache_size = self.max_idx
			else: self.cache_size = cache_size
			self.start_idx = self.start_idx - self.cache_size/2
			if refresh: self.__refresh_cache()
		else:
			if self.cache_size != self.max_idx:
				print "error, in list mode the cache size is always equal to the max idx"
				return
	def set_start_idx(self,start_idx,refresh=True):
		self.start_idx = start_idx
		if refresh: self.__refresh_cache()
	
	def __refresh_cache(self):
		#print "regenerating CACHE"
		app = QtGui.QApplication.instance()
		app.setOverrideCursor(Qt.BusyCursor)
		
		try:
			cache = {}
			for i in range(self.start_idx,self.start_idx+self.cache_size,1):
				if i != 0:
					idx = i % self.max_idx
				else: idx = 0
				try: 
					cache[idx] = self.images[idx]
				except:
					try:
						if self.mode ==  EMDataListCache.FILE_MODE:
							cache[idx] = EMData(self.file_name,idx)
							if self.exclusions.count(idx) != 0:
								cache[idx].set_attr("excluded",True)
						else:
							print "data has been lost"
							raise
					except: print "couldn't access",idx,"the max idx was",self.max_idx,"i was",i,"start idx",self.start_idx,"cache size",self.cache_size,len(self.images)
				#print i,idx
					
			self.images = cache
		except:
			print "there was an error in cache regeneration. Suggest restarting"
			
		app.setOverrideCursor(Qt.ArrowCursor)
	def __getitem__(self,idx):
		
		i = 0
		if idx != 0: i = idx%self.max_idx
		try:
			return self.images[i]
		except:
			self.start_idx = i - self.cache_size/2
			#if self.start_idx < 0: 
				#self.start_idx = self.start_idx % self.max_idx
			#elif self.start_idx+self.cache_size >= self.max_idx:
				#self.start_idx =  self.max_idx - self.cache_size/2 -1
			self.__refresh_cache()
			try:
				return self.images[i]
			except:
				print "error, couldn't get image",i,self.start_idx,self.max_idx,self.cache_size
				#for i in self.images:
					#print i,
				#print ''
			
# This is just for testing, of course
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	GLUT.glutInit("")
	window = EMImageMXRotor()
	if len(sys.argv)==1 : 
		data = []
		for i in range(500):
			#if i == 0: idx = 0
			#else: idx = i%64
			#e = EMData(64,64)
			#e.set_size(64,64,1)
			#e.to_zero()
			#e.add(sin( (i/10.0) % (pi/2)))
			data.append(test_image(Util.get_irand(0,3)))
			#data.append(e)
		
		
		window.setData(data)
	else :
		print "reading image"
		#a=EMData.read_images(sys.argv[1])
		window.set_image_file_name(sys.argv[1])
	window2=EMParentWin(window)
	window2.resize(*window.get_optimal_size())
	window2.show()

	
	sys.exit(app.exec_())
