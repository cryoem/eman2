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
from EMAN2 import *
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from OpenGL import GL,GLU,GLUT
from OpenGL.GL import *
from OpenGL.GLU import *
from valslider import ValSlider
from math import *
import sys
import numpy
from emimageutil import ImgHistogram,EMParentWin
from weakref import WeakKeyDictionary
from pickle import dumps,loads
from PyQt4.QtGui import QImage
from PyQt4.QtCore import QTimer
from emglobjects import EMOpenGLFlagsAndTools
from emfloatingwidgets import EMGLRotaryWidget, EMGLView2D, EM3DWidget
from emimagemx import EMImageMxInspector2D,EMImageMXCore
from EMAN2db import EMAN2DB

class EMImageMXRotary(QtOpenGL.QGLWidget):
	"""A QT widget for rendering EMData objects. It can display stacks of 2D images
	in 'matrix' form on the display. The middle mouse button will bring up a
	control-panel. The QT event loop must be running for this object to function
	properly.
	"""
	allim=WeakKeyDictionary()
	def __init__(self, data=None,parent=None):
		self.image_rotary = None
		#self.initflag = True
		self.mmode = "drag"

		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True);
		#fmt.setDepthBuffer(True)
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		EMImageMXRotary.allim[self]=0
		
		self.image_rotary = EMImageMXRotaryCore(data,self)
		
		self.imagefilename = None
		
		self.fov = 30
		self.aspect = 1.0
		self.z_near = 6000
		self.z_far = 13000
		
		self.animatables = []
		
		self.timer = QTimer()
		QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.timeout)
		self.timer.start(10)
		
		self.light_0_pos = [0.1,.1,1.,0.]
		
	def setData(self,data):
		self.image_rotary.setData(data)
	
	def get_optimal_size(self):
		lr = self.image_rotary.rotary.get_suggested_lr_bt_nf()
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
		self.image_rotary.set_image_file_name(name)
		
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
		
		if ( self.image_rotary == None ): return
		self.image_rotary.render()

	
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
		
		self.image_rotary.resize_event(width,height)
	
	def set_near_far(self,near,far):
		self.z_near = near
		self.z_far = far
		
		GL.glMatrixMode(GL.GL_PROJECTION)
		GL.glLoadIdentity()
		GLU.gluPerspective(self.fov,self.aspect,self.z_near,self.z_far)
		#GL.glOrtho(0.0,width,0.0,height,-width,width)
		GL.glMatrixMode(GL.GL_MODELVIEW)
		GL.glLoadIdentity()
		
	
		self.image_rotary.resize_event(-1,-1)
	def get_depth_for_height(self, height):
		# This function returns the width and height of the renderable 
		# area at the origin of the data volume
		depth = height/(2.0*tan(self.fov/2.0*pi/180.0))
	
		return depth
	
	def set_mmode(self,mode):
		self.mmode = mode
		self.image_rotary.set_mmode(mode)
	
	def mousePressEvent(self, event):
		self.image_rotary.mousePressEvent(event)
			
	def wheelEvent(self,event):
		self.image_rotary.wheelEvent(event)
	
	def mouseMoveEvent(self,event):
		self.image_rotary.mouseMoveEvent(event)

	def mouseReleaseEvent(self,event):
		self.image_rotary.mouseReleaseEvent(event)
	
	def keyPressEvent(self,event):
		if self.mmode == "app":
			self.emit(QtCore.SIGNAL("keypress"),event)

	def dropEvent(self,event):
		self.image_rotary.dropEvent(event)
		
	def closeEvent(self,event) :
		self.image_rotary.closeEvent(event)
		
	def dragEnterEvent(self,event):
		self.image_rotary.dragEnterEvent(event)

	def dropEvent(self,event):
		self.image_rotary.dropEvent(event)
	
	
	def set_shapes(self,shapes,shrink):
		self.image_rotary.set_shapes(shapes,shrink)
	
	def set_frozen(self,frozen):
		self.image_rotary.set_frozen(frozen)
	
	def get_frame_buffer(self):
		# THIS WILL FAIL ON WINDOWS APPARENTLY, because Windows requires a temporary context to be created and this is what the True flag
		# trying to stop.
		return self.renderPixmap(0,0,True)
		# to get around it we would have to render everything without display lists (a supreme pain).
	
class EMImageMXRotaryCore:

	#allim=WeakKeyDictionary()
	def __init__(self, data=None,parent=None):
		self.parent = parent
		self.data=None
		try: self.parent.setAcceptDrops(True)
		except:	pass

		self.initsizeflag = True
		if data:
			self.setData(data)
		
		self.rotary = EMGLRotaryWidget(self,-15,50,-15,EMGLRotaryWidget.TOP_ROTARY,60)
		#self.rotary.set_child_mouse_events(False)
		self.rotary.set_mmode("mxrotary")
		self.widget = EM3DWidget(self,self.rotary)
		self.widget.set_draw_frame(False)
		
		self.image_file_name = None	# keeps track of the image file name (if any) - book keeping purposes only
		self.emdata_list_cache = None # all import emdata list cache, the object that stores emdata objects efficiently. Must be initialized via setData or set_image_file_name
		
		self.visible_mxs = 3	# the number of visible imagemxs in the rotary
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
		
		self.vals_to_display = ["Img #"]
		
		self.z_near = 0
		self.z_far = 0
		
		self.s_flag = False
	
		try:
			self.font_renderer = EMFTGL()
			self.font_renderer.set_face_size(32)
			self.font_renderer.set_depth(8)
			self.font_renderer.set_using_display_lists(True)
			self.font_renderer.set_font_mode(FTGLFontMode.EXTRUDE)
			
			self.font_renderer.set_font_file_name("/usr/share/fonts/dejavu/DejaVuSerif.ttf")
			self.font_render_mode = EMImageMXCore.FTGL
		except:
			self.font_render_mode = EMImageMXCore.GLUT
	
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
	
	def set_mx_cols(self,cols):
		self.mx_cols = cols
		self.emdata_list_cache.set_cache_size(8*self.visible_mxs*self.mx_rows*self.mx_cols)
		self.__refresh_rotary(True)
		self.updateGL()
		
	def set_mx_rows(self,rows):
		self.mx_rows = rows
		self.emdata_list_cache.set_cache_size(8*self.visible_mxs*self.mx_rows*self.mx_cols)
		self.__refresh_rotary(True)
		self.updateGL()
	
	def set_display_values(self,v2d):
		self.vals_to_display = v2d
		self.__refresh_rotary(False)
		self.updateGL()

	def set_mxs(self,mxs):
		self.visible_mxs = mxs
		self.emdata_list_cache.set_cache_size(8*self.visible_mxs*self.mx_rows*self.mx_cols)
		self.__regenerate_rotary() # this could be done more efficiently
		self.updateGL()
	
	def set_mmode(self,mode):
		self.mmode = mode
		for i in range(self.visible_mxs):
			w = self.rotary[i].get_drawable()
			w.set_mmode(self.mmode)
	
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
	
	def pop_box_image(self,idx):
		val = self.emdata_list_cache.delete_box(idx)
		if val == 1:
			self.max_idx = self.emdata_list_cache.get_max_idx()
			
			self.__refresh_rotary()
		elif val == 2:
			w = self.rotary[0].get_drawable()
			w.force_dl_update()
		else:
			print 'failed to delete box image'
		
		
	def update_min_max_gamma(self):
		for i in range(self.visible_mxs):
			w = self.rotary[i].get_drawable()
			w.set_min_max_gamma(self.minden,self.maxden,self.gamma)
			if  i == 0 and  self.inspector != None:	
				self.inspector.set_hist(w.get_hist(),self.minden,self.maxden)
		
		self.updateGL()
	
	def set_den_range(self,minden,maxden):
		self.minden=minden
		self.maxden=maxden
		self.update_min_max_gamma()
		
	def emit(self,signal,event,integer=None):
		if integer != None:
			self.parent.emit(signal,event,integer)
		else:pass
				
	def update_rotary_position(self,inc):
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
		
		self.__update_rotary_range(start_changed,end_changed)
		
		w = self.rotary[idx].get_drawable()
		self.inspector.set_hist(w.get_hist(),self.minden,self.maxden)
	
	def __update_rotary_range(self,start_changed,end_changed):
		num_per_view = self.mx_rows*self.mx_cols
		for idx in range(start_changed,end_changed):
			n = idx + self.start_mx
			
			panel = n
			if panel != 0: panel %= self.get_num_panels()
			
			start_idx = panel*num_per_view
			w = self.rotary[idx].get_drawable()
			
			image_offset = start_idx 
			if image_offset != 0: image_offset %= self.emdata_list_cache.get_max_idx()
			#d = []
			#for i in range(start_idx,start_idx+num_per_view): d.append(self.emdata_list_cache[i])
			
			d = []
			if image_offset > (self.emdata_list_cache.get_max_idx()-num_per_view):
				num_visible = self.emdata_list_cache.get_max_idx() - image_offset
				num_none = num_per_view - num_visible
				print num_visible,num_none
				for i in range(start_idx,start_idx+num_visible): d.append(self.emdata_list_cache[i])
				for i in range(num_none): d.append(None)
			else:
				for i in range(start_idx,start_idx+num_per_view): d.append(self.emdata_list_cache[i])
				
			w.set_img_num_offset(image_offset)	
			w.setData(d,False)
			w.set_min_max_gamma(self.minden,self.maxden,self.gamma)
			w.set_max_idx(self.emdata_list_cache.get_max_idx())
			
	
	def get_num_panels(self):
		return int(ceil(float(self.emdata_list_cache.get_max_idx())/ (self.mx_rows*self.mx_cols)))

	def set_rotary_mode(self,mode):
		self.rot_mode = mode
		self.rotary.set_mmode(mode)

	def optimize_fit(self):
		render_width = self.parent.width()
		render_height = self.parent.height()
		
		self.mx_rows = render_width/self.emdata_list_cache.get_image_width()
		self.mx_cols = render_height/self.emdata_list_cache.get_image_height()
		
		self.inspector.set_n_cols(self.mx_cols)
		self.inspector.set_n_rows(self.mx_rows)

		self.__refresh_rotary(True)
		self.updateGL()
		
	def set_frozen(self,frozen):
		self.rotary.set_frozen(frozen)

	def set_shapes(self,shapes,shrink):
		self.rotary.set_shapes(shapes,shrink)

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
		
		self.emdata_list_cache = EMDataListCache(data,8*self.visible_mxs*self.mx_rows*self.mx_cols)
		self.__regenerate_rotary()
	
	def set_image_file_name(self,name):
		#print "set image file name",name
		self.image_file_name = name
		self.emdata_list_cache = EMDataListCache(name)
		self.__regenerate_rotary()
	
	def __refresh_rotary(self,inc_size=False):
		
		
		self.__update_rotary_range(0,self.visible_mxs)
		#for idx in range(0,self.visible_mxs):
			
			#panel = idx
			#if panel != 0: panel %= self.get_num_panels()		
			#start_idx = panel*num_per_view

			#e = self.rotary[idx]
			#w = e.get_drawable()

			#image_offset = start_idx
			#if image_offset != 0: image_offset %= self.emdata_list_cache.get_max_idx()
			#w.set_img_num_offset(image_offset)
			#d = []
			#if image_offset > (self.emdata_list_cache.get_max_idx()-num_per_view):
				#num_visible = self.emdata_list_cache.get_max_idx() - image_offset
				#num_none = num_per_view - num_visible
				#print num_visible,num_none
				#for i in range(start_idx,start_idx+num_visible): d.append(self.emdata_list_cache[i])
				#for i in range(num_none): d.append(None)
			#else:
				#for i in range(start_idx,start_idx+num_per_view): d.append(self.emdata_list_cache[i])

			#w.setData(d,False)

			#w.set_min_max_gamma(self.minden,self.maxden,self.gamma,False)
			#w.set_max_idx(self.emdata_list_cache.get_max_idx())
			#w.set_display_values(self.vals_to_display,False)
			
		if inc_size:
			self.__refresh_rotary_size()

	def __refresh_rotary_size(self):
		if len(self.rotary) == 0: return
		
		self.render_width = self.parent.width()
		self.render_height = self.parent.height()
		width = self.mx_rows*self.emdata_list_cache.get_image_width()
		height = self.mx_cols*self.emdata_list_cache.get_image_height()
		scale1 = self.render_height/float(height)
		scale2 = self.render_width/float(width)
		
		if self.render_width > self.render_height:
			self.render_height = float(height)/width*self.render_width
			scale = scale2
		else:
			self.render_width = float(width)/height*self.render_height
			scale = scale1
		
		for idx in range(0,self.visible_mxs):
			e = self.rotary[idx]
			w = e.get_drawable()

			w.set_scale(scale,False,False)
			e.set_width(self.render_width,False)
			e.set_height(self.render_height,False)
			w.set_mx_cols(self.mx_rows,False)
			
		self.rotary.update()

	def __regenerate_rotary(self):
		self.rotary.clear_widgets()
		num_per_view = self.mx_rows*self.mx_cols
		
		for idx in range(self.start_mx,self.visible_mxs+self.start_mx):
			panel = idx
			if panel != 0: panel %= self.get_num_panels()		
			start_idx = panel*num_per_view


			image_offset = start_idx
			if image_offset != 0: image_offset %= self.emdata_list_cache.get_max_idx()
			
			d = []
			if image_offset > (self.emdata_list_cache.get_max_idx()-num_per_view):
				num_visible = self.emdata_list_cache.get_max_idx() - image_offset
				num_none = num_per_view - num_visible
				print num_visible,num_none
				for i in range(start_idx,start_idx+num_visible): d.append(self.emdata_list_cache[i])
				for i in range(num_none): d.append(None)
			else:
				for i in range(start_idx,start_idx+num_per_view): d.append(self.emdata_list_cache[i])
			e = EMGLView2D(self,d)
			w = e.get_drawable()
			w.set_img_num_offset(image_offset)
			w.set_max_idx(self.emdata_list_cache.get_max_idx())
			w.set_use_display_list(True) # saves HEAPS of time, makes interaction much smoother
			w.set_draw_background(True)
			w.set_reroute_delete_target(self)
			w.set_mmode(self.mmode)
			w.set_display_values(self.vals_to_display,False)
			self.rotary.add_widget(e)

		self.__refresh_rotary_size()
		
		w = self.rotary[0].get_drawable()
		self.minden = w.get_density_min()
		self.maxden = w.get_density_max()
		self.mindeng = self.minden
		self.maxdeng = self.maxden
		self.gamma = w.get_gamma()
		
	def updateGL(self):
		try: self.parent.updateGL()
		except: pass


	def render(self):
		if not self.parent.isVisible(): return
		if self.emdata_list_cache == None: return
	
		glLoadIdentity()
		
		lr = self.rotary.get_suggested_lr_bt_nf()
		suppress = False
		GL.glEnable(GL.GL_DEPTH_TEST)
		GL.glEnable(GL.GL_LIGHTING)
		z = self.parent.get_depth_for_height(abs(lr[3]-lr[2]))
		lrt = self.widget.get_lr_bt_nf()
		z_near = z-lrt[4]
		z_trans = 0
		z_far = z-lrt[5]
		if z_near < 0:
			z_trans = z_near
			z_near = 0.1
			z_far -= z_trans
		if z_far < 0: z_far = 0.1 # hacking alert
		if self.z_near != z_near or self.z_far != z_far:
			self.z_near = z_near
			self.z_far = z_far
			self.parent.set_near_far(self.z_near,self.z_far)
			suppress = True
		
		
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
		if event.button()==Qt.MidButton or (event.button()==Qt.LeftButton and event.modifiers()&Qt.ControlModifier):
			self.show_inspector(True)
		else:
			self.widget.mousePressEvent(event)
			self.updateGL()
	
	def mouseMoveEvent(self, event):
		self.widget.mouseMoveEvent(event)
		self.updateGL()
		
	def mouseReleaseEvent(self, event):
		self.widget.mouseReleaseEvent(event)
		self.updateGL()
		
	def wheelEvent(self, event):
#		if event.delta() > 0:
#			self.setScale( self.scale * self.mag )
#		elif event.delta() < 0:
#			self.setScale(self.scale * self.invmag )
#		self.resizeEvent(self.parent.width(),self.parent.height())
#		# The self.scale variable is updated now, so just update with that
#		if self.inspector: self.inspector.setScale(self.scale)
		self.widget.wheelEvent(event)
		self.updateGL()
	def leaveEvent(self):
		pass

	def resize_event(self, width, height):
		self.rotary.resize_event(width,height)
		self.__refresh_rotary_size()
	
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
		glOrtho(0,width,0,height,-100,100)
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		glEnable(GL_LIGHTING)
		glEnable(GL_NORMALIZE)
		glMaterial(GL_FRONT,GL_AMBIENT,(0.9, 0.5, 0.2,1.0))
		glMaterial(GL_FRONT,GL_DIFFUSE,(0.9, 0.5, 0.2,1.0))
		glMaterial(GL_FRONT,GL_SPECULAR,(0.9, 0.5, 0.2,1.0))
		glMaterial(GL_FRONT,GL_SHININESS,20.0)
		
		glDisable(GL_DEPTH_TEST)
		glColor(1.0,1.0,1.0)
		
		glTranslate(10,height-1.2*self.font_renderer.get_face_size(),0)
		glRotate(20,0,1,0)
		if self.font_render_mode == EMImageMXCore.FTGL:
			panels = self.get_num_panels()
			idx = self.start_mx
			if idx != 0: idx %= panels
			string = str(idx+1) + ' / ' + str(panels)
			self.font_renderer.render_string(string)
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
		global DB
		if isinstance(object,list):
			# in list mode there is no real caching
			self.mode = EMDataListCache.LIST_MODE
			self.max_idx = len(object)
			self.cache_size = self.max_idx
			self.images = object
			self.start_idx = 0
			DB.open_dict("emimage_mx_rotary_cache")
			self.db = DB.emimage_mx_rotary_cache
			self.exclusions_key = "interactive_exclusions"
			self.db[self.exclusions_key] = []
			self.exclusions = self.db["interactive_exclusions"]
			
			for i,d in enumerate(data):
				d.set_attr("original_number",i)
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
			
			DB.open_dict("emimage_mx_rotary_cache")
			self.db = DB.emimage_mx_rotary_cache
			self.exclusions_key = self.file_name+"_interactive_exclusions"
			self.exclusions  = self.db[self.exclusions_key]
			if self.exclusions == None: self.exclusions = []
			self.__refresh_cache()
		else:
			print "the object used to construct the EMDataListCache is not a string (filename) or a list (of EMData objects). Can't proceed"
			return
		
		self.soft_delete = True # toggle to prevent permanent deletion of particles
		self.image_width = -1
		self.image_height = -1
	
	def __del__(self):
		global DB
		try:
			print "closign dict"
			DB.close_dict("emimage_mx_rotary_cache")
			print "dict closed"
		except: pass
	
	def get_image_width(self):
		idx = self.start_idx
		if idx != 0: idx = idx % self.max_idx
		return self.images[idx].get_xsize()
		
	def get_image_height(self):
		idx = self.start_idx
		if idx != 0: idx = idx % self.max_idx
		return self.images[idx].get_ysize()
		
	def delete_box(self,idx):
		if self.mode == EMDataListCache.LIST_MODE and not self.soft_delete:
			# we can actually delete the emdata object
			image = self.images.pop(idx)
			self.max_idx = len(self.images)
			self.cache_size = self.max_idx
			#self.__refresh_cache()
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
	
	def set_cache_size(self,cache_size,refresh=True):
		if cache_size > self.max_idx: self.cache_size = self.max_idx
		else: self.cache_size = cache_size
		if refresh: self.__refresh_cache()
		
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
					except: print "couldn't access",idx	
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
			self.start_idx = idx - self.cache_size/2
			#if self.start_idx < 0: 
				#self.start_idx = self.start_idx % self.max_idx
			#elif self.start_idx+self.cache_size >= self.max_idx:
				#self.start_idx =  self.max_idx - self.cache_size/2 -1
			self.__refresh_cache()
			try:
				return self.images[i]
			except:
				print "error, couldn't get image",i,self.start_idx,self.max_idx,self.cache_size
				for i in self.images:
					print i,
				print ''
			
# This is just for testing, of course
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	GLUT.glutInit("")
	window = EMImageMXRotary()
	if len(sys.argv)==1 : 
		data = []
		for i in range(500):
			data.append(test_image(Util.get_irand(0,3)))
		
		
		window.setData(data)
	else :
		print "reading image"
		#a=EMData.read_images(sys.argv[1])
		window.set_image_file_name(sys.argv[1])
	window2=EMParentWin(window)
	window2.resize(*window.get_optimal_size())
	window2.show()
#	w2=QtGui.QWidget()
#	w2.resize(256,128)
	
#	w3=ValSlider(w2)
#	w3.resize(256,24)
#	w2.show()
	
	sys.exit(app.exec_())
