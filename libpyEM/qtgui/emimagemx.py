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
from valslider import ValSlider
from math import *
from EMAN2 import *
import EMAN2
import copy
import sys
import numpy
from emimageutil import ImgHistogram,EMEventRerouter,EMParentWin
from weakref import WeakKeyDictionary
from pickle import dumps,loads
from PyQt4.QtGui import QImage
from PyQt4.QtCore import QTimer
from libpyGLUtils2 import *

from emglobjects import EMOpenGLFlagsAndTools,EMGLProjectionViewMatrices,EMBasicOpenGLObjects
from emapplication import EMStandAloneApplication, EMGUIModule,get_application
from emanimationutil import LineAnimation
import weakref

from emapplication import EMProgressDialogModule

class EMImageMXWidget(QtOpenGL.QGLWidget,EMEventRerouter,EMGLProjectionViewMatrices):
	"""
	"""
	def __init__(self, em_mx_module,parent=None):
		#self.initflag = True
		self.mmode = "drag"

		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True)
		#fmt.setSampleBuffers(True)
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		EMEventRerouter.__init__(self,em_mx_module)
		EMGLProjectionViewMatrices.__init__(self)
		
		self.imagefilename = None
		
		#self.resize(480,480)
	
	
	def get_target(self):
		return self.target()
	
	def set_data(self,data):
		self.target().set_data(data)
	
	def set_file_name(self,name):
		#print "set image file name",name
		self.imagefilename = name
		
	def get_image_file_name(self):
		return self.imagefilename
	
	def initializeGL(self):
		glClearColor(0,0,0,0)
		
		
		glEnable(GL_LIGHTING)
		glEnable(GL_LIGHT0)
		#glEnable(GL_DEPTH_TEST)
		#print "Initializing"
		glLightfv(GL_LIGHT0, GL_AMBIENT, [0.0, 0.0, 0.0, 1.0])
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0,1.0,1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_POSITION, [1,1,1.,0.])
		glLightfv(GL_LIGHT0, GL_AMBIENT, [0.4, 0.4, 0.4, 1.0])
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, [1,0,-1.,1.])
		glLightfv(GL_LIGHT0, GL_POSITION, [40,40,100.,1.])
		#glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,GL_TRUE) # this is intenionally turned off

		glEnable(GL_CULL_FACE)
		glCullFace(GL_BACK)
	def paintGL(self):
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		
		if ( self.target == None ): return
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		#context = OpenGL.contextdata.getContext(None)
		#print "Matrix context is", context
		self.target().render()

	
	def resizeGL(self, width, height):
		if width <= 0 or height <= 0: return
		GL.glViewport(0,0,width,height)
	
		GL.glMatrixMode(GL.GL_PROJECTION)
		GL.glLoadIdentity()
		GL.glOrtho(0.0,width,0.0,height,-50,50)
		GL.glMatrixMode(GL.GL_MODELVIEW)
		GL.glLoadIdentity()
		
		glLightfv(GL_LIGHT0, GL_POSITION, [width/2,height/2,100.,1.])
		
		self.set_projection_view_update()
		try: self.target().resize_event(width,height)
		except: pass
	def set_mouse_mode(self,mode):
		self.mmode = mode
		self.target().set_mouse_mode(mode)
	
	def get_frame_buffer(self):
		# THIS WILL FAIL ON WINDOWS APPARENTLY, because Windows requires a temporary context - but the True flag is stopping the creation of a temporary context
		# (because display lists are involved)
		# December 2008 - tested in MAC it doesn't work,only returns a blank image
		return self.renderPixmap(0,0,True)

	def closeEvent(self,event):
		self.target().clear_gl_memory()
		self.target().closeEvent(None)
		
		

class EMMXCoreMouseEvents:
	'''
	A base class for objects that handle mouse events in the EMImageMXModule
	'''
	def __init__(self,mediator):
		'''
		Stores only a reference to the mediator
		'''
		if not isinstance(mediator,EMMXCoreMouseEventsMediator):
			print "error, the mediator should be a EMMXCoreMouseEventsMediator"
			return
		self.mediator = mediator
		
	def mouse_up(self,event):
		'''
		Inheriting classes to potentially define this function
		'''
		pass

	def mouse_down(self,event):
		'''
		Inheriting classes to potentially define this function
		'''
		pass
		
	def mouse_drag(self,event):
		'''
		Inheriting classes to potentially define this function
		'''
		pass
	
	def mouse_move(self,event):
		'''
		Inheriting classes to potentially define this function
		'''
		pass

	def mouse_wheel(self,event):
		'''
		Inheriting classes to potentially define this function
		'''
		pass
	
	def mouse_double_click(self,event):
		'''
		Inheriting classes to potentially define this function
		'''
		pass

class EMMXCoreMouseEventsMediator:
	def __init__(self,target):
		if not isinstance(target,EMImageMXModule):
			print "error, the target should be a EMImageMXModule"
			return
		
		self.target = weakref.ref(target)
		
	def get_image_file_name(self):
		return self.target().get_image_file_name()
	
	def scr_to_img(self,vec):
		return self.target().scr_to_img(vec)
	
	def get_parent(self):
		return self.target().get_parent()

	def get_box_image(self,idx):
		return self.target().get_box_image(idx)
	
	def pop_box_image(self,idx,event=None,redraw=False):
		self.target().pop_box_image(idx,event,redraw)
		
	def image_set_associate(self,idx,event=None,redraw=False):
		self.target().image_set_associate(idx,event,redraw)
	
	def get_density_max(self):
		return self.target().get_density_max()
	
	def get_density_min(self):
		return self.target().get_density_min()
	
	def emit(self,*args,**kargs):
		self.target().emit(*args,**kargs)

	def set_selected(self,selected,update_gl=True):
		self.target().set_selected(selected,update_gl)
		
	def force_display_update(self):
		self.target().force_display_update() 
		
	def get_scale(self):
		return self.target().get_scale()
	
	def update_inspector_texture(self):
		self.target().update_inspector_texture()

class EMMXDelMouseEvents(EMMXCoreMouseEvents):
	def __init__(self,mediator):
		EMMXCoreMouseEvents.__init__(self,mediator)
		self.lc = None
		
	def mouse_down(self,event):
		if event.button()==Qt.LeftButton:
			self.lc=self.mediator.scr_to_img((event.x(),event.y()))
			
	def mouse_up(self,event):
		if event.button()==Qt.LeftButton:
			lc=self.mediator.scr_to_img((event.x(),event.y()))
			if lc != None and self.lc != None and lc[0] == self.lc[0]:
				self.mediator.pop_box_image(lc[0],event,True)
				#self.mediator.force_display_update()
				
			self.lc = None
				
class EMMXSetMouseEvents(EMMXCoreMouseEvents):
	def __init__(self,mediator):
		EMMXCoreMouseEvents.__init__(self,mediator)
		self.lc = None
	def mouse_down(self,event):
		if event.button()==Qt.LeftButton:
			self.lc= self.mediator.scr_to_img((event.x(),event.y()))
			
	def mouse_up(self,event):
		if event.button()==Qt.LeftButton:
			lc=self.mediator.scr_to_img((event.x(),event.y()))
			if lc != None and self.lc != None and lc[0] == self.lc[0]:
				self.mediator.image_set_associate(lc[0],event,True)
				#self.mediator.force_display_update()
				
			self.lc = None
				
				
class EMMXDragMouseEvents(EMMXCoreMouseEvents):
	def __init__(self,mediator):
		EMMXCoreMouseEvents.__init__(self,mediator)
		self.class_window = None # used if people are looking at class averages and they double click - which case a second window is opend showing the particles in the class
		
	def mouse_down(self,event):
		return
	   	# this is currently disabled because it causes seg faults on MAC. FIXME investigate and establish the functionality that we want for mouse dragging and dropping
		if event.button()==Qt.LeftButton:
			lc= self.mediator.scr_to_img((event.x(),event.y()))
			if lc == None: 
				print "strange lc error"
				return
			box_image = self.mediator.get_box_image(lc[0])
			xs=int(box_image.get_xsize())
			ys=int(box_image.get_ysize())
			drag = QtGui.QDrag(self.mediator.get_parent())
			mime_data = QtCore.QMimeData()
			
			mime_data.setData("application/x-eman", dumps(box_image))
			
			EMAN2.GUIbeingdragged= box_image	# This deals with within-application dragging between windows
			mime_data.setText( str(lc[0])+"\n")
			di=QImage(GLUtil.render_amp8(box_image, 0,0,xs,ys,xs*4,1.0,0,255,self.mediator.get_density_min(),self.mediator.get_density_max(),1.0,14),xs,ys,QImage.Format_RGB32)
			mime_data.setImageData(QtCore.QVariant(di))
			drag.setMimeData(mime_data)
	
			# This (mini image drag) looks cool, but seems to cause crashing sometimes in the pixmap creation process  :^(
			#di=QImage(GLUtil.render_amp8(elf.data[lc[0]],0,0,xs,ys,xs*4,1.0,0,255,self.minden,self.maxden,14),xs,ys,QImage.Format_RGB32)
			#if xs>64 : pm=QtGui.QPixmap.fromImage(di).scaledToWidth(64)
			#else: pm=QtGui.QPixmap.fromImage(di)
			#drag.setPixmap(pm)
			#drag.setHotSpot(QtCore.QPoint(12,12))
					
			dropAction = drag.start()
	
	def mouse_double_click(self,event):
		'''
		Inheriting classes to potentially define this function
		'''
		lc=self.mediator.scr_to_img((event.x(),event.y()))
		if lc != None:
			a = self.mediator.get_box_image(lc[0])
			d = a.get_attr_dict()
			if d.has_key("class_ptcl_src") and d.has_key("class_ptcl_idxs"):
				data = []
				idxs = d["class_ptcl_idxs"]
				name = d["class_ptcl_src"]
				progress = QtGui.QProgressDialog("Reading images from %s" %get_file_tag(name), "Cancel", 0, len(idxs),None)
				progress.show()
				get_application().setOverrideCursor(Qt.BusyCursor)
				
				i = 0
				for idx in idxs:
					data.append(EMData(name,idx))
					i+=1
					progress.setValue(i)
 					get_application().processEvents()
 				
		 			if progress.wasCanceled():
		 				progress.close()
		 				get_application().setOverrideCursor(Qt.ArrowCursor)
			 			return
		
				progress.close()
				get_application().setOverrideCursor(Qt.ArrowCursor)
	
				resize_necessary = False
				if self.class_window == None:
					self.class_window = EMImageMXModule()
					QtCore.QObject.connect(self.class_window,QtCore.SIGNAL("module_closed"),self.on_class_window_closed)
					resize_necessary = True
				
				self.class_window.set_data(data,"Class Particles")
				
				if resize_necessary:
					get_application().show_specific(self.class_window)
					self.class_window.optimally_resize()
				else:
					self.class_window.updateGL()
	def on_class_window_closed(self):
		self.class_window = None
	
class EMMAppMouseEvents(EMMXCoreMouseEvents):
	def __init__(self,mediator):
		EMMXCoreMouseEvents.__init__(self,mediator)
	
	def mouse_down(self,event):
		if event.button()==Qt.LeftButton:
			lc=self.mediator.scr_to_img((event.x(),event.y()))
			if lc:
				self.mediator.emit(QtCore.SIGNAL("mx_image_selected"),event,lc)
				#print "setting selected"
				self.mediator.set_selected([lc[0]],True)
			xians_stuff = False
			if xians_stuff:
				if lc[0] != None:
					image = self.mediator.get_box_image(lc[0])
    				cx = image.get_xsize()/2
    				cy = image.get_ysize()/2
    				x = lc[1]-cx
    				y = lc[2]-cy
    				angle = atan2(y,x)*180.0/pi
    				# convert from clockwise to anti clockwise and convert to Xian's axis definition
    				angle = 270-angle
    				angle %= 360
    				print "Filename: ", self.mediator.get_image_file_name()
    				print "Sequence#: ",lc[0]
    				print "Angle: ", angle
    			 	
			
	def mouse_move(self,event):
		if event.buttons()&Qt.LeftButton:
			self.mediator.emit(QtCore.SIGNAL("mx_mousedrag"),event,self.mediator.get_scale())
	
	def mouse_up(self,event):
		if event.button()==Qt.LeftButton:
			lc=self.mediator.scr_to_img((event.x(),event.y()))
		
			if  not event.modifiers()&Qt.ShiftModifier:
				self.mediator.emit(QtCore.SIGNAL("mx_mouseup"),event,lc)
			else:
				if lc != None:
					self.mediator.pop_box_image(lc[0],event,True)
					self.mediator.force_display_update()
			
class EMMatrixPanel:
	'''
	A class for managing the parameters of displaying a matrix panel
	'''
	def __init__(self):
		self.min_sep = 2 # minimum separation
		self.max_y = 0
		self.height = 0 # height of the panel
		self.visiblecols = 0 # the number of visible images in the x direction
		self.xoffset = 0 # Offset for centering the display images in x
		self.visiblerows = 0 # the number of visible rows
		self.rowstart = 0
		
		self.scale_cache = {}
	
	def update_panel_params(self,view_width,view_height,view_scale,view_data,y,target):
		rendered_image_width = (view_data.get_xsize())*view_scale
		rendered_image_height = view_data.get_ysize()*view_scale
		
		[self.ystart,self.visiblerows,self.visiblecols] = self.visible_row_col(view_width,view_height,view_scale,view_data,y)
		if self.ystart == None: 
			scale = self.get_min_scale(view_width,view_height,view_scale,view_data)
			target.scale = scale
			view_scale = scale
			[self.ystart,self.visiblerows,self.visiblecols] = self.visible_row_col(view_width,view_height,view_scale,view_data,y)
			if self.ystart == None: 
				print "kill me now....please"
		xsep = view_width - self.visiblecols*(rendered_image_width+self.min_sep)
		self.xoffset = xsep/2		
		
		self.height = ceil(len(view_data)/float(self.visiblecols))*(rendered_image_height+self.min_sep) + self.min_sep
		self.max_y = self.height - view_height # adjusted height is the maximum value for current y!
	
	def visible_row_col(self,view_width,view_height,view_scale,view_data,y):
		rendered_image_width = (view_data.get_xsize())*view_scale
		rendered_image_height = view_data.get_ysize()*view_scale
		
		visiblecols = int(floor(view_width/(rendered_image_width+self.min_sep)))
		
		if visiblecols == 0: 
			#print "scale is too large - the panel can not be rendered"
			# The calling function should be sophisticated enough to handle this
			return [None,None,None]
		
		
		w=int(min(rendered_image_width,view_width))
		h=int(min(rendered_image_height,view_height))
		
		yoff = 0
		if y < 0:
			ybelow = floor(-y/(h+2))
			yoff = ybelow*(h+2)+y
			visiblerows = int(ceil(float(view_height-yoff)/(h+2)))
		else: visiblerows = int(ceil(float(view_height-y)/(h+2)))
				
		maxrow = int(ceil(float(len(view_data))/visiblecols))
		rowstart =-y/(h+2)
		if rowstart < 0: rowstart = 0
		elif rowstart > 0:
			rowstart = int(rowstart)
			visiblerows = visiblerows + rowstart
		if visiblerows > maxrow: visiblerows = maxrow
		
		visiblerows = int(visiblerows)
		rowstart = int(rowstart)
		
		return [rowstart,visiblerows,visiblecols]
	
	
	def basic_visible_row_col(self,view_width,view_height,view_scale,view_data):
		rendered_image_width = (view_data.get_xsize())*view_scale
		rendered_image_height = view_data.get_ysize()*view_scale
		
		visiblecols = int(floor(view_width/(rendered_image_width+self.min_sep)))
		visiblerows = int(floor(view_height/(rendered_image_height+self.min_sep)))
		
		return [visiblerows,visiblecols]
		
	def get_min_scale(self,view_width,view_height,view_scale,view_data):
		'''
		Gets min scale using something like a bifurcation algorithm
		
		'''
		
		s = str(view_width) + str(view_height)
		if self.scale_cache.has_key(s):
			return self.scale_cache[s]
		
		n = len(view_data)
		
		[visiblerows,visiblecols] = self.basic_visible_row_col(view_width,view_height,view_scale,view_data)
		
		if self.auto_scale_logic(visiblerows,visiblecols,n) == 1:
			min_scale = view_scale
			max_scale = min_scale
			while (True):
				max_scale *= 2
				[visiblerows,visiblecols] = self.basic_visible_row_col(view_width,view_height,max_scale,view_data)
				if self.auto_scale_logic(visiblerows,visiblecols,n) == -1:
					break
		else:
			max_scale = view_scale
			min_scale = max_scale
			while (True):
				min_scale /= 2
				[visiblerows,visiblecols] = self.basic_visible_row_col(view_width,view_height,min_scale,view_data)
				if self.auto_scale_logic(visiblerows,visiblecols,n) == 1:
					break
				
		
		while (True):
			estimate_scale = (min_scale + max_scale)*0.5
			[visiblerows,visiblecols] = self.basic_visible_row_col(view_width,view_height,estimate_scale,view_data)
			#[rowstart,visiblerows,visiblecols] = self.visible_row_col(view_width,view_height,estimate_scale,view_data,-self.min_sep)
			if self.auto_scale_logic(visiblerows,visiblecols,n) >= 0:
				min_scale = estimate_scale
			elif self.auto_scale_logic(visiblerows,visiblecols,n) == -1:
				max_scale = estimate_scale
		
			if abs(min_scale-max_scale) < 0.01:
				self.scale_cache[s] = max_scale
				return max_scale
				
	
		
	def auto_scale_logic(self,visiblerows,visiblecols,n):
		if (visiblerows*visiblecols) >= n :
			if ((visiblerows-1)*visiblecols) < n:

				return 0 # scale is perfect
			else:
				return 1 # scale is too small
		else: 
			return -1 # scale is too large
		
		
#		print self.height,self.xsep,self.visiblecols
	
		

class EMImageMXModule(EMGUIModule):
	
	def enable_set(self,idx,name,display=True,lst=[]):
		'''
		Called from e2eulerxplor
		'''
		self.get_inspector()
		self.inspector.add_set(idx,name,display)
		self.data.associate_set(idx,lst)
	
	def clear_sets(self):
		self.inspector.clear_sets()
		self.data.clear_sets()
		self.updateGL()
	
	def set_single_active_set(self,db_name,idx=0,exclusive=0):
		'''
		Called from emform
		'''
		db = db_open_dict("bdb:select")
		set_list = db[db_name]
		if set_list == None: set_list = []
		
		
		self.get_inspector()
		self.inspector.clear_sets()
		self.inspector.add_set(idx,db_name,True)
		self.data.clear_sets()
		self.data.associate_set(idx,set_list,True)
		self.inspector.set_set_mode()
		self.updateGL()
		
	def load_font_renderer(self):
		try:
			self.font_render_mode = EMGUIModule.FTGL
			self.font_renderer = get_3d_font_renderer()
			self.font_renderer.set_face_size(self.font_size)
			self.font_renderer.set_font_mode(FTGLFontMode.TEXTURE)
		except:
			self.font_render_mode = EMGUIModule.GLUT
	
	def get_qt_widget(self):
		if self.qt_context_parent == None:	
			self.under_qt_control = True
			self.gl_context_parent = EMImageMXWidget(self)
			self.qt_context_parent = EMParentWin(self.gl_context_parent)
			self.gl_widget = self.gl_context_parent
			self.qt_context_parent.setWindowIcon(QtGui.QIcon(get_image_directory() +"multiple_images.png"))
			#self.optimally_resize()
			#self.qt_context_parent.setAcceptDrops(True)
		
		return self.qt_context_parent
	
	def get_gl_widget(self,qt_context_parent,gl_context_parent):
		from emfloatingwidgets import EM2DGLView, EM2DGLWindow
		self.init_size_flag = False
		if self.gl_widget == None:
			self.under_qt_control = False
			self.gl_context_parent = gl_context_parent
			self.qt_context_parent = qt_context_parent
			self.draw_background = True
			
			gl_view = EM2DGLView(self,image=None)
			self.gl_widget = EM2DGLWindow(self,gl_view)
			self.gl_widget.target_translations_allowed(True)
			self.update_window_title(self.file_name)
		return self.gl_widget
	
	def get_desktop_hint(self):
		return self.desktop_hint
	allim=WeakKeyDictionary()
	def __init__(self, data=None,application=None):
		self.desktop_hint = "image"
		self.init_size_flag = True
		self.data=None
		EMGUIModule.__init__(self,ensure_gl_context=True)
		EMImageMXModule.allim[self] = 0
		self.file_name = ''
		self.datasize=(1,1)
		self.scale=1.0
		self.minden=0
		self.maxden=1.0
		self.invert=0
		self.fft=None
		self.mindeng=0
		self.maxdeng=1.0
		self.gamma=1.0
		self.nshow=-1
		self.mousedrag=None
		self.nimg=0
		self.changec={}
		self.mmode="drag"
		self.selected=[]
		self.hist = []
		self.targetorigin=None
		self.targetspeed=20.0
		self.mag = 1.1				# magnification factor
		self.invmag = 1.0/self.mag	# inverse magnification factor
		self.glflags = EMOpenGLFlagsAndTools() 	# supplies power of two texturing flags
		self.tex_names = [] 		# tex_names stores texture handles which are no longer used, and must be deleted
		self.suppress_inspector = False 	# Suppresses showing the inspector - switched on in emfloatingwidgets
		self.image_file_name = None
		self.first_render = True # a hack, something is slowing us down in FTGL
		self.scroll_bar = EMGLScrollBar(self)
		self.draw_scroll = True
		self.scroll_bar_has_mouse = False
		self.matrix_panel = EMMatrixPanel()
		self.origin=(self.matrix_panel.min_sep,self.matrix_panel.min_sep)
		self.emdata_list_cache = None # all import emdata list cache, the object that stores emdata objects efficiently. Must be initialized via set_data or set_image_file_name
		self.animation_enabled = False
		self.line_animation = None
		
		self.coords={}
		self.nshown=0
		self.valstodisp=["Img #"]
		
		self.inspector=None
		
		self.font_size = 11
		self.load_font_renderer()
		if data:
			self.set_data(data,False)
			
		self.text_bbs = {} # bounding box cache - key is a string, entry is a list of 6 values defining a 
		
		self.use_display_list = True # whether or not a display list should be used to render the main view - if on, this will save on time if the view is unchanged
		self.main_display_list = 0	# if using display lists, the stores the display list
		self.display_states = [] # if using display lists, this stores the states that are checked, and if different, will cause regeneration of the display list
		self.draw_background = False # if true  will paint the background behind the images black using a polygon - useful in 3D contexts, ie i the emimagemxrotary
		
		
		self.img_num_offset = 0		# used by emimagemxrotary for display correct image numbers
		self.max_idx = 99999999		# used by emimagemxrotary for display correct image numbers
	
		self.__init_mouse_handlers()
		
		self.reroute_delete_target = None
	
	def __del__(self):
		#handler = QtCore.QObjectCleanupHandler()
		#handler.add(self)
		#handler.remove(self)
		#self.setObjectName("me")
		self.clear_gl_memory()
	
	def get_emit_signals_and_connections(self):
		return {"set_origin":self.set_origin,"set_scale":self.set_scale,"origin_update":self.origin_update}
	
	def width(self):
		if self.gl_widget != None:
			return self.view_width()
		else:
			return 0
	
	def using_ftgl(self):
		return self.font_render_mode == EMGUIModule.FTGL
	
	def get_font_size(self):
		return self.font_renderer.get_face_size()
		
	def set_font_size(self,value):
		self.font_renderer.set_face_size(value)
		self.force_display_update() # only for redoing the fonts, this could be made more efficient :(
		self.updateGL()
		
	def height(self):
		if self.gl_widget != None:
			return self.gl_widget.height()
		else:
			return 0
	
	def __init_mouse_handlers(self):
		self.mouse_events_mediator = EMMXCoreMouseEventsMediator(self)
		self.mouse_event_handlers = {}
		self.mouse_event_handlers["app"] = EMMAppMouseEvents(self.mouse_events_mediator)
		self.mouse_event_handlers["del"] = EMMXDelMouseEvents(self.mouse_events_mediator)
		self.mouse_event_handlers["drag"] = EMMXDragMouseEvents(self.mouse_events_mediator)
		self.mouse_event_handlers["set"] = EMMXSetMouseEvents(self.mouse_events_mediator)
		self.mouse_event_handler = self.mouse_event_handlers[self.mmode]
		
	def set_mouse_mode(self,mode):
		self.mmode = mode
		meh  = self.mouse_event_handler
		try:
			self.mouse_event_handler = self.mouse_event_handlers[self.mmode]
		except:
			print "unknown mode:",mode
			self.mouse_event_handler = meh # just keep the old one
	
	def set_file_name(self,name):
		#print "set image file name",name
		self.image_file_name = name
	
	def get_inspector(self):
		if not self.inspector : 
			self.inspector=EMImageInspectorMX(self)
			self.inspector_update()
		return self.inspector
	
	def inspector_update(self):
		#FIXME
		pass
		#print "do inspector update"
		
	
	def set_reroute_delete_target(self,target):
		self.reroute_delete_target = target
	
	def get_scale(self):
		return self.scale
	
	def pop_box_image(self,idx,event=None,update_gl=False):
		val = self.data.delete_box(idx)
		if val > 0 and update_gl:
			self.force_display_update()
			self.updateGL()
			if event != None: self.emit(QtCore.SIGNAL("mx_boxdeleted"),event,[idx],False) 
	
	
	def image_set_associate(self,idx,event=None,update_gl=False):
		val = self.data.image_set_associate(idx)
		if val > 0 and update_gl:
			self.force_display_update()
			self.updateGL()
		self.get_inspector()
		self.inspector.save_set()
		
#	def pop_box_image(self,idx,event=None,update_gl=False):
#		if self.reroute_delete_target  == None:
#			d = self.data.pop(idx)
#			self.display_states = []
#			if event != None: self.emit(QtCore.SIGNAL("mx_boxdeleted"),event,[idx],False)
#			if update_gl:
#				self.display_states = [] 
#				self.updateGL()
#			return d
#		else:
#			self.reroute_delete_target.pop_box_image(idx)
#			if event != None: self.emit(QtCore.SIGNAL("mx_boxdeleted"),event,[idx],False)

	def get_box_image(self,idx):
		return self.data[idx]

	#def emit(self,*args,**kargs):
		#qt_widget = self.application.get_qt_emitter(self)
		#qt_widget.emit(*args,**kargs)
	
	def clear_gl_memory(self):
		if self.main_display_list != 0:
			glDeleteLists(self.main_display_list,1)
			self.main_display_list = 0
		if self.tex_names != None and ( len(self.tex_names) > 0 ):	
			glDeleteTextures(self.tex_names)
			self.tex_names = []	
	
	def force_display_update(self):
		''' If display lists are being used this will force a regeneration'''
		self.display_states = []
		self.update_inspector_texture()
	
	def set_img_num_offset(self,n):
		self.img_num_offset = n
		self.force_display_update() # empty display lists causes an automatic regeneration of the display list
	
	def get_img_num_offset(self):
		return self.img_num_offset
	
	def set_draw_background(self,bool):
		self.draw_background = bool
		self.force_display_update()# empty display lists causes an automatic regeneration of the display list
		
	def set_use_display_list(self,bool):
		self.use_display_list = bool
		
	def set_max_idx(self,n):
		self.force_display_update()# empty display lists causes an automatic regeneration of the display list
		self.max_idx = n
		
	def set_min_max_gamma(self,minden,maxden,gamma,update_gl=True):
		self.minden= minden
		self.maxden= maxden
		self.gamma = gamma
		if update_gl: self.updateGL()
	def get_hist(self): return self.hist
	def get_image(self,idx): return self.data[idx]
	def get_image_file_name(self):
		''' warning - could return none in some circumstances'''
		return self.file_name
	
	def optimally_resize(self):
		
	
		if isinstance(self.gl_context_parent,EMImageMXWidget):
			self.qt_context_parent.resize(*self.get_parent_suggested_size())
		else:
			self.gl_widget.resize(*self.get_parent_suggested_size())
			
		# I disabled this because it was giving me problems
#		self.scale =  self.matrix_panel.get_min_scale(self.view_width(),self.gl_widget.height(),self.scale,self.data) # this is to prevent locking
#		print self.scale
		
	def get_parent_suggested_size(self):
		if self.data != None and isinstance(self.data[0],EMData): 
			if len(self.data)<self.matrix_panel.visiblecols :
				w=len(self.data)*(self.data.get_xsize()+2)
				hfac = 1
			else : 
				w=self.matrix_panel.visiblecols*(self.data.get_xsize()+2)
				hfac = len(self.data)/self.matrix_panel.visiblecols+1
			hfac *= self.data.get_ysize()
			if hfac > 512: hfac = 512
			if w > 512: w = 512
			return (int(w)+12,int(hfac)+12) # the 12 is related to the EMParentWin... hack...
		else: return (512+12,512+12)
	
	def update_window_title(self,filename):
		if isinstance(self.gl_context_parent,EMImageMXWidget):
			self.qt_context_parent.setWindowTitle(remove_directories_from_name(filename))
		else:
			if self.gl_widget != None:
				self.gl_widget.setWindowTitle(remove_directories_from_name(filename))
	
	def set_data(self,data_or_filename,filename='',update_gl=True):
		'''
		This function will work if you give it a list of EMData objects, or if you give it the file name (as the first argument)
		If this solution is undesirable one could easily split this function into two equivalents.
		'''
		cache_size = -1
		if isinstance(data_or_filename,str):
			nx,ny,nz = gimme_image_dimensions3D(data_or_filename)
			bytes_per_image = nx*ny*4.0
			cache_size = int(75000000/bytes_per_image) # 75 MB
			
		
		self.data = EMDataListCache(data_or_filename,cache_size)

		self.file_name = filename
		self.update_window_title(self.file_name)
		
		self.force_display_update()
		self.nimg=len(self.data)
		self.max_idx = len(self.data)
		if self.nimg == 0: return # the list is empty
			
		global HOMEDB
		HOMEDB=EMAN2db.EMAN2DB.open_db()
		HOMEDB.open_dict("display_preferences")
		db = HOMEDB.display_preferences
		auto_contrast = db.get("display_stack_auto_contrast",dfl=True)
		start_guess = db.get("display_stack_np_for_auto",dfl=5)
		
		d = self.data[0]

	   	mean = d.get_attr("mean")
	   	sigma = d.get_attr("sigma")
		m0=d.get_attr("minimum")
		m1=d.get_attr("maximum")
		
		if auto_contrast:
			mn=max(m0,mean-3.0*sigma)
			mx=min(m1,mean+3.0*sigma)
		else:
			mn=m0
			mx=m1
			
		self.minden=mn
		self.maxden=mx
		self.mindeng=m0
		self.maxdeng=m1

		if start_guess > len(self.data) or start_guess < 0 :start_guess = len(self.data)
		for j in range(1,start_guess): #
			d = self.data.get_image_header(j)
			if d == None: continue
			if d["nz"] !=1 :
				self.data=None
				if update_gl: self.updateGL()
				return
			try:
				mean=d["mean"]
				sigma=d["sigma"]
				m0=d["minimum"]
				m1=d["maximum"]
			except:
				print d
				print "failed key"
			
			if sigma == 0: continue


			if auto_contrast:
				mn=max(m0,mean-3.0*sigma)
				mx=min(m1,mean+3.0*sigma)
			else:
				mn=m0
				mx=m1
		
			
			self.minden=min(self.minden,mn)
			self.maxden=max(self.maxden,mx)
			self.mindeng=min(self.mindeng,m0)
			self.maxdeng=max(self.maxdeng,m1)

		if self.inspector: self.inspector.set_limits(self.mindeng,self.maxdeng,self.minden,self.maxden)
		
		#if update_gl: self.updateGL()

	def updateGL(self):
		if self.gl_widget != None and self.under_qt_control:
			self.gl_widget.updateGL()
			
	def set_den_range(self,x0,x1,update_gl=True):
		"""Set the range of densities to be mapped to the 0-255 pixel value range"""
		self.minden=x0
		self.maxden=x1
		if update_gl: self.updateGL()
	
	def get_density_min(self):
		return self.minden
	
	def get_density_max(self):
		return self.maxden
	
	def get_gamma(self):
		return self.gamma
	
	def set_origin(self,x,y,update_gl=True):
		"""Set the display origin within the image"""
		if self.animation_enabled:
			if self.line_animation != None and self.line_animation.animated: return # this is so the current animation has to end before starting another one. It could be the other way but I like it this way
			self.line_animation = LineAnimation(self,self.origin,(x,y))
			self.get_qt_context_parent().register_animatable(self.line_animation)
			return True
		else:
			self.origin=(x,y)
			self.targetorigin=None
		if update_gl: self.updateGL()
		
	def set_line_animation(self,x,y):
		self.origin=(x,y)
		self.updateGL()
	
	def animation_done_event(self,animation):
		pass
		
	def set_scale(self,newscale,update_gl=True):
		"""Adjusts the scale of the display. Tries to maintain the center of the image at the center"""
			
		if self.data and len(self.data)>0 and (self.data.get_ysize()*newscale>self.gl_widget.height() or self.data.get_xsize()*newscale>self.view_width()):
			newscale=min(float(self.gl_widget.height())/self.data.get_ysize(),float(self.view_width())/self.data.get_xsize())
			if self.inspector: self.inspector.scale.setValue(newscale)
		
		oldscale = self.scale
		self.scale=newscale
		self.matrix_panel.update_panel_params(self.view_width(),self.gl_widget.height(),self.scale,self.data,self.origin[1],self)
		
		
		view_height = self.gl_widget.height()
		panel_height = self.matrix_panel.height
		if panel_height < view_height :
			# This auto rescaling stuff was disabled at the requeset of Steve Ludtke. Uncomment to see how it works
			#if oldscale > newscale: self.scale =  self.matrix_panel.get_min_scale(self.view_width(),self.gl_widget.height(),self.scale,self.data) # this is to prevent locking
			self.draw_scroll = False
			self.origin=(self.matrix_panel.min_sep,self.matrix_panel.min_sep)
			self.matrix_panel.update_panel_params(self.view_width(),self.gl_widget.height(),self.scale,self.data,self.origin[1],self)
		else:
			self.draw_scroll = True
			self.scroll_bar.update_target_ypos()	
		
		if self.emit_events: self.emit(QtCore.SIGNAL("set_scale"),self.scale,adjust,update_gl)
		if update_gl: self.updateGL()
	
	def resize_event(self, width, height):
		self.scroll_bar.height = height
		
		self.matrix_panel.update_panel_params(self.view_width(),height,self.scale,self.data,self.origin[1],self)
		view_height = self.gl_widget.height()
		panel_height = self.matrix_panel.height
		if panel_height < view_heights :
			# This auto rescaling stuff was disabled at the requeset of Steve Ludtke. Uncomment to see how it works
			#self.scale =  self.matrix_panel.get_min_scale(self.view_width(),self.gl_widget.height(),self.scale,self.data) # this is to prevent locking
			self.draw_scroll = False
			self.origin=(self.matrix_panel.min_sep,self.matrix_panel.min_sep)
			self.matrix_panel.update_panel_params(self.view_width(),self.gl_widget.height(),self.scale,self.data,self.origin[1],self)
		else:
			self.draw_scroll = True
			self.scroll_bar.update_target_ypos()	
			
		
	def set_density_min(self,val,update_gl=True):
		self.minden=val
		if update_gl: self.updateGL()
		
	def set_density_max(self,val,update_gl=True):
		self.maxden=val
		if update_gl: self.updateGL()

	def get_mmode(self):
		return self.mmode
		
	def set_gamma(self,val,update_gl=True):
		self.gamma=val
		if update_gl:self.updateGL()
	
	def set_n_show(self,val,update_gl=True):
		self.nshow=val
		if update_gl: self.updateGL()

	def set_invert(self,val,update_gl=True):
		if val: self.invert=1
		else : self.invert=0
		if update_gl: self.updateGL()

	def timeout(self):
		"""Called a few times each second when idle for things like automatic scrolling"""
		if self.targetorigin :
			vec=(self.targetorigin[0]-self.origin[0],self.targetorigin[1]-self.origin[1])
			h=hypot(vec[0],vec[1])
			if h<=self.targetspeed :
				self.origin=self.targetorigin
				self.targetorigin=None
			else :
				vec=(vec[0]/h,vec[1]/h)
				self.origin=(self.origin[0]+vec[0]*self.targetspeed,self.origin[1]+vec[1]*self.targetspeed)
			#self.updateGL()
	
	def display_state_changed(self):
		display_states = []
		display_states.append(self.gl_widget.width())
		display_states.append(self.gl_widget.height())
		display_states.append(self.origin[0])
		display_states.append(self.origin[1])
		display_states.append(self.scale)
		display_states.append(self.invert)
		display_states.append(self.minden)
		display_states.append(self.maxden)
		display_states.append(self.gamma)
		display_states.append(self.draw_background)
		display_states.append(self.img_num_offset)
		if len(self.display_states) == 0:
			self.display_states = display_states
			return True
		else:
			for i in range(len(display_states)):
				
				if display_states[i] != self.display_states[i]:
					self.display_states = display_states
					return True
		
		return False
			
	
	def set_font_render_resolution(self):
		if self.font_render_mode != EMGUIModule.FTGL:
			print "error, can't call set_font_render_resolution if the mode isn't FTGL"
			
		#self.font_renderer.set_face_size(int(self.gl_widget.height()*0.015))
		#print "scale is",self.scale
	
	def __draw_backdrop(self):
		light = glIsEnabled(GL_LIGHTING)
		glDisable(GL_LIGHTING)
	
		glColor(.9,.9,.9)
		glBegin(GL_QUADS)
		glVertex(0,0,-1)
		glColor(.9,.9,.9)
		glVertex(self.gl_widget.width(),0,-1)
		glColor(.9,.9,.9)
		glVertex(self.gl_widget.width(),self.gl_widget.height(),-1)
		glColor(.9,.9,.9)
		glVertex(0,self.gl_widget.height(),-1)
		glEnd()
		if light: glEnable(GL_LIGHTING)
	
	
	def view_width(self):
		return self.gl_widget.width() - self.draw_scroll*self.scroll_bar.width
	
	def render(self):
		if not self.data : return
		if self.font_render_mode == EMGUIModule.FTGL: self.set_font_render_resolution()
#		self.image_change_count = self.data.get_changecount() # this is important when the user has more than one display instance of the same image, for instance in e2.py if 
		render = False
		
		update = False
		if self.display_state_changed():
			update = True
			
		if update:
			self.update_inspector_texture() # important for this to occur in term of the e2desktop only
			
		if self.use_display_list:
			
			if update:
				if self.main_display_list != 0:
					glDeleteLists(self.main_display_list,1)
					self.main_display_list = 0

			if self.main_display_list == 0:
				self.main_display_list = glGenLists(1)
				glNewList(self.main_display_list,GL_COMPILE)
				render = True
		else: render = True
		
		self.matrix_panel.update_panel_params(self.view_width(),self.gl_widget.height(),self.scale,self.data,self.origin[1],self)
		
		if render: 
			if self.draw_background:
				self.__draw_backdrop()
		
			#for i in self.data:
				#self.changec[i]=i.get_attr("changecount")
			
			if not self.invert : pixden=(0,255)
			else: pixden=(255,0)
			
			n=len(self.data)
			self.hist=numpy.zeros(256)
			#if len(self.coords)>n : self.coords=self.coords[:n] # dont know what this does? Had to comment out, changing from a list to a dictionary
			glColor(0.5,1.0,0.5)
			glLineWidth(2)

			if self.font_render_mode == EMGUIModule.GLUT:
				init_glut()
				try:
					# we render the 16x16 corner of the image and decide if it's light or dark to decide the best way to 
					# contrast the text labels...
					a=GLUtil.render_amp8(self.data[0],0,0,16,16,16,self.scale,pixden[0],pixden[1],self.minden,self.maxden,self.gamma,4)
					ims=[ord(pv) for pv in a]
					if sum(ims)>32768 : txtcol=(0.0,0.0,0.0,1.0)
					else : txtcol=(1.0,1.0,1.0,1.0)
				except: txtcol=(1.0,1.0,1.0,1.0)
			else: txtcol=(1,1,1,1)
			
			if ( len(self.tex_names) > 0 ):	glDeleteTextures(self.tex_names)
			self.tex_names = []
	
			self.nshown=0
			
			#x,y=-self.origin[0],-self.origin[1]
			x,y = self.matrix_panel.xoffset,int(self.origin[1])
			w=int(min(self.data.get_xsize()*self.scale,self.view_width()))
			h=int(min(self.data.get_ysize()*self.scale,self.gl_widget.height()))
				
			invscale=1.0/self.scale
			self.set_label_ratio = 0.1
			self.coords = {}
			nsets = len(self.data.visible_sets)
			current_sets = copy.copy(self.data.visible_sets)
			current_sets.sort()
			
			for row in range(self.matrix_panel.ystart,self.matrix_panel.visiblerows):
				for col in range(0,self.matrix_panel.visiblecols):
					i = int((row)*self.matrix_panel.visiblecols+col)
					
					#print i,n
					if i >= n : break
					tx = int((w+self.matrix_panel.min_sep)*(col) + x)
					ty = int((h+self.matrix_panel.min_sep)*(row) + y)
					real_y = ty # need this for set display
					tw = w
					th = h
					rx = 0	#render x
					ry = 0	#render y
					#print "Prior",i,':',row,col,tx,ty,tw,th,y,x
					drawlabel = True
					if (tx+tw) > self.view_width():
						tw = int(self.view_width()-tx)
					elif tx<0:
						drawlabel=False
						rx = int(ceil(-tx*invscale))
						tw=int(w+tx)
						tx = 0
						
	
					#print h,row,y
					#print "Prior",i,':',row,col,tx,ty,tw,th,'offsets',yoffset,xoffset
					if (ty+th) > self.gl_widget.height():
						#print "ty + th was greater than",self.gl_widget.height()
						th = int(self.gl_widget.height()-ty)
					elif ty<0:
						drawlabel = False
						ry = int(ceil(-ty*invscale))
						th=int(h+ty)
						ty = 0
						
					#print i,':',row,col,tx,ty,tw,th,'offsets',yoffset,xoffset
					if th < 0 or tw < 0:
						#weird += 1
						#print "weirdness"
						#print col,row,
						continue
					
					self.coords[i]=(tx,ty,tw,th)
					
					draw_tex = True
					d = self.data[i]
					if not d.has_attr("excluded"):
						#print rx,ry,tw,th,self.gl_widget.width(),self.gl_widget.height(),self.origin
						if not self.glflags.npt_textures_unsupported():
							a=GLUtil.render_amp8(self.data[i],rx,ry,tw,th,(tw-1)/4*4+4,self.scale,pixden[0],pixden[1],self.minden,self.maxden,self.gamma,2)
							self.texture(a,tx,ty,tw,th)
						else:
							a=GLUtil.render_amp8(self.data[i],rx,ry,tw,th,(tw-1)/4*4+4,self.scale,pixden[0],pixden[1],self.minden,self.maxden,self.gamma,6)
							glRasterPos(tx,ty)
							glDrawPixels(tw,th,GL_LUMINANCE,GL_UNSIGNED_BYTE,a)
							
						hist2=numpy.fromstring(a[-1024:],'i')
						self.hist+=hist2	
					else:
						# d is excluded/deleted
						self.load_set_color(0)
						glPushMatrix()
						glTranslatef(tx+tw/2.0,ty+th/2.0,0)
						glScale(tw/2.0,th/2.0,1.0)
						self.__render_excluded_square()
						glPopMatrix()
						
								
						
						
					if drawlabel:
						self.__draw_mx_text(tx,ty,txtcol,i)
						
					if draw_tex : self.nshown+=1
					
					if nsets>0 and hasattr(d,"mxset") and len(d.mxset) > 0:
						light = glIsEnabled(GL_LIGHTING)
						glEnable(GL_LIGHTING)
						
						iss = 0
						for rot,set in enumerate(current_sets):
							#s
							if set in d.mxset:
								x_pos_ratio = 1-self.set_label_ratio
								y_pos_ratio = x_pos_ratio
								y_pos_ratio -= 2*float(iss)/nsets*self.set_label_ratio
								y_pos_ratio *= 2
								x_pos_ratio *= 2
								self.load_set_color(set)
								width = w/2.0
								height = h/2.0
								
								glPushMatrix()
								glTranslatef(tx+x_pos_ratio*width,real_y+y_pos_ratio*height,0)
								#glScale((1-0.2*rot/nsets)*width,(1-0.2*rot/nsets)height,1.0)
								#self.__render_excluded_hollow_square()
								glScale(self.set_label_ratio*width,self.set_label_ratio*height,1.0)
								self.__render_excluded_square()
								
								glPopMatrix()
								iss  += 1
								
						if not light: glEnable(GL_LIGHTING)
			for i in self.selected:
				try:
					data = self.coords[i]	
					glColor(0.5,0.5,1.0)
					glBegin(GL_LINE_LOOP)
					glVertex(data[0],data[1])
					glVertex(data[0]+data[2],data[1])
					glVertex(data[0]+data[2],data[1]+data[3])
					glVertex(data[0],data[1]+data[3])
					glEnd()
				except:
					# this means the box isn't visible!
					pass
			
			if self.inspector : self.inspector.set_hist(self.hist,self.minden,self.maxden)
		else:
			try:
				glCallList(self.main_display_list)
			except: pass
		
		if self.use_display_list and render :
			glEndList()
			glCallList(self.main_display_list)
			if self.first_render: #A hack, FTGL is slowing us down
				self.display_states = []
				self.first_render = False 
		if self.draw_scroll: self.draw_scroll_bar()
	
	def load_set_color(self,set):
		color = BoxingTools.get_color(set+1)
		c = [color[0],color[1],color[2],1.0]
		glMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,c)
		glMaterial(GL_FRONT,GL_SPECULAR,c)
	  	glMaterial(GL_FRONT,GL_SHININESS,100.0)
#		if set == 0:
#			glMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,(.2,.2,.8,1.0))
#			glMaterial(GL_FRONT,GL_SPECULAR,(.2,.2,.8,1.0))
#		  	glMaterial(GL_FRONT,GL_SHININESS,100.0)
#		else:
#			glMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,(.8,.2,.2,1.0))
#			glMaterial(GL_FRONT,GL_SPECULAR,(.8,.2,.2,1.0))
#		  	glMaterial(GL_FRONT,GL_SHININESS,100.0)
	def __render_excluded_square(self):
		
		glBegin(GL_QUADS)
		#glColor(0,0,0)
		glNormal(-.1,.1,0.1)
		glVertex(-1,1,0.1)
		#glColor(0.15,0.15,0.15)
		glNormal(-1,1,-0.1)
		glVertex(-1,-1,0.1)	
		#glColor(0.3,0.3,0.3)
		glNormal(.1,-.1,1)
		glVertex(1,-1,0.1)
		#glColor(0.22,0.22,0.22)
		glNormal(1,-1,0.1)
		glVertex(1,1,0.1)
		glEnd()
		#glPopMatrix()
	
	def __render_excluded_hollow_square(self):
		'''
		Renders a hollow square that goes from -1 to 1 and has a thickness of 10%
		'''
		d = .1
		glBegin(GL_QUADS)
		#glColor(0,0,0)
		glNormal(-.1,.1,1)
		glVertex(-1,1,0.1)
		#glColor(0.15,0.15,0.15)
		glNormal(-1,1,-0.1)
		glVertex(-1,-1,0.1)	
		#glColor(0.3,0.3,0.3)
		glNormal(.1,-.1,1)
		glVertex(-1+d,-1,0.1)
		#glColor(0.22,0.22,0.22)
		glNormal(1,-1,0.1)
		glVertex(-1+d,1,0.1)
		
		
		
		glNormal(-.1,.1,1)
		glVertex(1,1,0.1)
		#glColor(0.15,0.15,0.15)
		#glColor(0.22,0.22,0.22)
		glNormal(1,-1,0.1)
		glVertex(1-d,1,0.1)
		
		#glColor(0.3,0.3,0.3)
		glNormal(.1,-.1,1)
		glVertex(1-d,-1,0.1)
		
		glNormal(-1,1,-0.1)
		glVertex(1,-1,0.1)	
		
		
		glNormal(1,-1,0.1)
		glVertex(-1+d,1,0.1)
		
		glNormal(1,1,0.1)
		glVertex(-1+d,1-d,0.1)
		
		glNormal(1,1,0.1)
		glVertex(1-d,1-d,0.1)
		
		glNormal(1,-1,0.1)
		glVertex(1-d,1,0.1)
		
		
		glNormal(1,1,0.1)
		glVertex(-1+d,-1+d,0.1)
		
		glNormal(1,-1,0.1)
		glVertex(-1+d,-1,0.1)
		
		glNormal(1,1,0.1)
		glVertex(1-d,-1,0.1)
		
		glNormal(1,-1,0.1)
		glVertex(1-d,-1+d,0.1)
		
		
		glEnd()
	
	def __draw_mx_text(self,tx,ty,txtcol,i):
		bgcol = [1-v for v in txtcol]
		if self.font_render_mode == EMGUIModule.FTGL:
			
			glEnable(GL_TEXTURE_2D)
			glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE)
			lighting = glIsEnabled(GL_LIGHTING)
			glDisable(GL_LIGHTING)
			#glEnable(GL_NORMALIZE)
			tagy = ty
			
			glMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,txtcol)
			glMaterial(GL_FRONT,GL_SPECULAR,txtcol)
			glMaterial(GL_FRONT,GL_SHININESS,1.0)
			for v in self.valstodisp:
				glPushMatrix()
				glTranslate(tx,tagy,0)
				glTranslate(0,1,0.2)
				if v=="Img #" : 
					#print i,self.img_num_offset,self.max_idx,(i+self.img_num_offset)%self.max_idx,
					idx = i+self.img_num_offset
					if idx != 0: idx = idx%self.max_idx
					sidx = str(idx)
					bbox = self.font_renderer.bounding_box(sidx)
					GLUtil.mx_bbox(bbox,txtcol,bgcol)
					self.font_renderer.render_string(sidx)
					
				else : 
					av=self.data[i].get_attr(v)
					try:
						if isinstance(av,float) : avs="%1.4g"%av
						elif isinstance(av,Transform):
							t = av.get_rotation("eman")
							avs = "%1.1f,%1.1f,%1.1f" %(t["az"],t["alt"],t["phi"])
						elif isinstance(av,EMAN2Ctf):
							avs = "%1.3f,%1.1f" %(av.defocus,av.bfactor)
						else: avs=str(av)
					except:avs ="---"
					bbox = self.font_renderer.bounding_box(avs)
					
					GLUtil.mx_bbox(bbox,txtcol,bgcol)
					self.font_renderer.render_string(avs)

				tagy+=self.font_renderer.get_face_size()

				glPopMatrix()
			if lighting:
				glEnable(GL_LIGHTING)
			glDisable(GL_TEXTURE_2D)
		elif self.font_render_mode == EMGUIModule.GLUT:
			tagy = ty
			glColor(*txtcol)
			glDisable(GL_LIGHTING)
			for v in self.valstodisp:
				if v=="Img #" :
					idx = i+self.img_num_offset
					if idx != 0: idx = idx%self.max_idx
					self.render_text(tx,tagy,"%d"%idx)
				else : 
					av=self.data[i].get_attr(v)
					if isinstance(av,float) : avs="%1.4g"%av
					else: avs=str(av)
					try: self.render_text(tx,tagy,str(avs))
					except:	self.render_text(tx,tagy,"------")
				tagy+=16
								
	
	def bounding_box(self,character):
		try: self.text_bbs[character]
		except:
			 self.text_bbs[character] = self.font_renderer.bounding_box(character)
			 
		return self.text_bbs[character]
	
	def texture(self,a,x,y,w,h):
		
		tex_name = glGenTextures(1)
		if ( tex_name <= 0 ):
			raise("failed to generate texture name")
		
		width = w/2.0
		height = h/2.0
		
		glPushMatrix()
		glTranslatef(x+width,y+height,0)
			
		glBindTexture(GL_TEXTURE_2D,tex_name)
		glPixelStorei(GL_UNPACK_ALIGNMENT,4)
		glTexImage2D(GL_TEXTURE_2D,0,GL_LUMINANCE,w,h,0,GL_LUMINANCE,GL_UNSIGNED_BYTE, a)
		
		glEnable(GL_TEXTURE_2D)
		glBindTexture(GL_TEXTURE_2D, tex_name)
		
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP)
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP)
		# using GL_NEAREST ensures pixel granularity
		# using GL_LINEAR blurs textures and makes them more difficult
		# to interpret (in cryo-em)
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)
		# this makes it so that the texture is impervious to lighting
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE)
		
		
		# POSITIONING POLICY - the texture occupies the entire screen area
		glBegin(GL_QUADS)
		
		glTexCoord2f(0,0)
		glVertex2f(-width,height)
				
		glTexCoord2f(0,1)
		glVertex2f(-width,-height)
			
		glTexCoord2f(1,1)
		glVertex2f(width,-height)
		
		glTexCoord2f(1,0)
		glVertex2f(width,height)
		glEnd()
		
		glDisable(GL_TEXTURE_2D)
		
		glPopMatrix()
		self.tex_names.append(tex_name)
	
	def render_text(self,x,y,s):
#	     print 'in render Text'
		glRasterPos(x+2,y+2)
		for c in s:
			GLUT.glutBitmapCharacter(GLUT.GLUT_BITMAP_9_BY_15,ord(c))

				
	def is_visible(self,n):
		try:
			val = self.coords[n][4]
			return True
		except: return False
	
	def scroll_to(self,n,yonly=0):
		"""Moves image 'n' to the center of the display"""
#		print self.origin,self.coords[0],self.coords[1]
#		try: self.origin=(self.coords[n][0]-self.width()/2,self.coords[n][1]+self.height()/2)
#		try: self.origin=(self.coords[8][0]-self.width()/2-self.origin[0],self.coords[8][1]+self.height()/2-self.origin[1])
		if yonly :
			try: 
				self.targetorigin=(0,self.coords[n][1]-self.gl_widget.height()/2+self.data.get_ysize()*self.scale/2)
			except: return
		else:
			try: self.targetorigin=(self.coords[n][0]-self.view_width()/2+self.data.get_xsize()*self.scale/2,self.coords[n][1]-self.gl_widget.height()/2+self.data.get_ysize()*self.scale/2)
			except: return
		self.targetspeed=hypot(self.targetorigin[0]-self.origin[0],self.targetorigin[1]-self.origin[1])/20.0
#		print n,self.origin
#		self.updateGL()
	
	def set_selected(self,numlist,update_gl=True):
		"""pass an integer or a list/tuple of integers which should be marked as 'selected' in the
		display"""
		real_numlist = []
		for i in numlist:
			t = i-self.get_img_num_offset()
			if t != 0:
				t %= len(self.data)
			real_numlist.append(t)

		if isinstance(numlist,int) : numlist=[real_numlist]
		if isinstance(numlist,list) or isinstance(real_numlist,tuple) : self.selected=real_numlist
		else : self.selected=[]
		self.force_display_update()
		if update_gl: self.updateGL()
	
	def set_display_values(self,v2d,update_gl=True):
		"""Pass in a list of strings describing image attributes to overlay on the image, in order of display"""
		v2d.reverse()
		self.valstodisp=v2d
		self.display_states = []
		if update_gl: self.updateGL()

	def scr_to_img(self,vec):
		"""Converts screen location (ie - mouse event) to pixel coordinates within a single
		image from the matrix. Returns (image number,x,y) or None if the location is not within any
		of the contained images. """ 
		if  self.max_idx == 0: return # there is no data
		
		absloc=((vec[0]),(self.gl_widget.height()-(vec[1])))
		for item in self.coords.items():
			index = item[0]+self.img_num_offset
			if index != 0: index %= self.max_idx
			data = item[1]
			if absloc[0]>data[0] and absloc[1]>data[1] and absloc[0]<data[0]+data[2] and absloc[1]<data[1]+data[3] :
				return (index,(absloc[0]-data[0])/self.scale,(absloc[1]-data[1])/self.scale)
		return None
		
	def dragEnterEvent(self,event):
#		f=event.mime_data().formats()
#		for i in f:
#			print str(i)
		
		if event.source()==self:
			event.setDropAction(Qt.MoveAction)
			event.accept()
		elif event.provides("application/x-eman"):
			event.setDropAction(Qt.CopyAction)
			event.accept()

	def save_data(self):
		'''
		Overwrites existing data (the old file is removed)
		'''
		
		msg = QtGui.QMessageBox()
		msg.setWindowTitle("Woops")
		if self.data==None or len(self.data)==0:
			msg.setText("there is no data to save" %fsp)
			msg.exec_()
			return
		
		from emsave import save_data
		file_name = save_data(self.data)
		if file_name == self.file_name and file_exists(file_name): # the file we are working with was overwritten
			self.set_data(file_name)

	def save_lst(self,fsp):
		'''
		If we make it here the dialog has taken care of check whether or not overwrite should occur
		'''
		
		origname = self.get_image_file_name()

		f = file(fsp,'w')
		f.write('#LST\n')
		
		progress = EMProgressDialogModule(get_application(),"Writing files", "abort", 0, len(self.data),None)
		progress.qt_widget.show()
		for i in xrange(0,len(self.data)):
			d = self.data[i]
			if d == None: continue # the image has been excluded
			progress.qt_widget.setValue(i)
			get_application().processEvents()
			if progress.qt_widget.wasCanceled():
				#remove_file(fsp)# we could do this but if they're overwriting the original data then they lose it all
				f.close()
				progress.qt_widget.close()
				return
				#pass
		progress.qt_widget.close()
		f.close()
	
	def dropEvent(self,event):
		lc=self.scr_to_img((event.pos().x(),event.pos().y()))
		if event.source()==self:
#			print lc
			n=int(event.mime_data().text())
			if not lc : lc=[len(self.data)]
			if n>lc[0] : 
				self.data.insert(lc[0],self.data[n])
				del self.data[n+1]
			else : 
				self.data.insert(lc[0]+1,self.data[n])
				del self.data[n]
			event.setDropAction(Qt.MoveAction)
			event.accept()
		elif EMAN2.GUIbeingdragged:
			self.data.append(EMAN2.GUIbeingdragged)
			self.set_data(self.data)
			EMAN2.GUIbeingdragged=None
		elif event.provides("application/x-eman"):
			x=loads(event.mime_data().data("application/x-eman"))
			if not lc : self.data.append(x)
			else : self.data.insert(lc[0],x)
			self.set_data(self.data)
			event.acceptProposedAction()

	def keyPressEvent(self,event):
		if self.data == None: return
	
		if event.key() == Qt.Key_F1:
			self.display_web_help("http://blake.bcm.edu/emanwiki/EMAN2/Programs/emimagemx")
		elif event.key()==Qt.Key_Up :
			self.scroll_bar.up_button_pressed(False)
		elif event.key()==Qt.Key_Down:
			self.scroll_bar.down_button_pressed(False)
		elif event.key()==Qt.Key_W or event.key()==Qt.Key_PageUp:
			self.scroll_bar.scroll_up(False)
		elif event.key()==Qt.Key_S or event.key()==Qt.Key_PageDown:
			self.scroll_bar.scroll_down(False)
		elif event.key()==Qt.Key_Left :
			pass
			#self.origin=(self.origin[0]+xstep,self.origin[1])
			#self.updateGL()
		elif event.key()==Qt.Key_Right:
			pass
			#self.origin=(self.origin[0]-xstep,self.origin[1])
			#self.updateGL()
		else:
			return
		
		if self.emit_events: self.emit(QtCore.SIGNAL("origin_update"),self.origin)
		
	def check_newy(self,y):
		newy = y
		if newy > self.matrix_panel.min_sep: newy = self.matrix_panel.min_sep # this is enforcing a strong boundary at the bottom
		elif newy < -self.matrix_panel.max_y: newy = -self.matrix_panel.max_y
		return newy
	
	def origin_update(self,new_origin):
		self.origin = new_origin
	
	def mouseDoubleClickEvent(self,event):
		self.mouse_event_handler.mouse_double_click(event)

	def mousePressEvent(self, event):
		if (self.gl_widget.width()-event.x() <= self.scroll_bar.width):
			self.scroll_bar_has_mouse = True
			self.scroll_bar.mousePressEvent(event)
			return
		
		if event.button()==Qt.MidButton or (event.button()==Qt.LeftButton and event.modifiers()&Qt.AltModifier):
			self.show_inspector(1)
			self.inspector.set_limits(self.mindeng,self.maxdeng,self.minden,self.maxden)

#			self.emit(QtCore.SIGNAL("inspector_shown"),event)
		elif event.button()==Qt.RightButton or (event.button()==Qt.LeftButton and event.modifiers()&Qt.AltModifier):
			if not self.draw_scroll: return # if the (vertical) scroll bar isn't drawn then mouse movement is disabled (because the images occupy the whole view)
			app =  QtGui.QApplication.instance()
			try:
				get_application().setOverrideCursor(Qt.ClosedHandCursor)
				#app.setOverrideCursor(Qt.ClosedHandCursor)
			except: # if we're using a version of qt older than 4.2 than we have to use this...
				get_application().setOverrideCursor(Qt.SizeAllCursor)
				#app.setOverrideCursor(Qt.SizeAllCursor)
				
			self.mousedrag=(event.x(),event.y())
		else:
			self.mouse_event_handler.mouse_down(event)
		
	def mouseMoveEvent(self, event):
		if self.scroll_bar_has_mouse:
			self.scroll_bar.mouseMoveEvent(event)
			return
		
		if self.mousedrag and self.draw_scroll: # if the (vertical) scroll bar isn't drawn then mouse movement is disabled (because the images occupy the whole view)
			oldy = self.origin[1]
			newy = self.origin[1]+self.mousedrag[1]-event.y()
			newy = self.check_newy(newy)
			if newy == oldy: return # this won't hold if the zoom level is great and an image occupies more than the size of the display
	
			#self.origin=(self.origin[0]+self.mousedrag[0]-event.x(),self.origin[1]-self.mousedrag[1]+event.y())
			self.origin=(self.matrix_panel.xoffset,newy)
			if self.emit_events: self.emit(QtCore.SIGNAL("set_origin"),self.origin[0],self.origin[1],False)
			self.mousedrag=(event.x(),event.y())
			try:self.gl_widget.updateGL()
			except: pass
		else: 
			self.mouse_event_handler.mouse_move(event)
		
	def mouseReleaseEvent(self, event):
		if self.scroll_bar_has_mouse:
			self.scroll_bar.mouseReleaseEvent(event)
			self.scroll_bar_has_mouse = False
			return
		
		get_application().setOverrideCursor(Qt.ArrowCursor)
		lc=self.scr_to_img((event.x(),event.y()))
		if self.mousedrag:
			self.mousedrag=None
		else: self.mouse_event_handler.mouse_up(event)
			
	def wheelEvent(self, event):
		if event.delta() > 0:
			self.set_scale( self.scale * self.mag )
		elif event.delta() < 0:
			self.set_scale(self.scale * self.invmag)
		#self.resize_event(self.gl_widget.width(),self.gl_widget.height())
		# The self.scale variable is updated now, so just update with that
		if self.inspector: self.inspector.set_scale(self.scale)
	
	def leaveEvent(self,event):
		get_application().setOverrideCursor(Qt.ArrowCursor)
		if self.mousedrag:
			self.mousedrag=None
			
	def get_frame_buffer(self):
		return self.gl_widget.get_frame_buffer()
	
	def draw_scroll_bar(self):
		width = self.gl_widget.width()
		glEnable(GL_LIGHTING)
		glEnable(GL_NORMALIZE)
		glDisable(GL_TEXTURE_2D)

		glPushMatrix()
		glTranslate(width-self.scroll_bar.width,0,0)
		self.scroll_bar.draw()
		glPopMatrix()
		
class EMGLScrollBar:
	def __init__(self,target):
		self.min = 0
		self.max = 0
		self.current_pos = 0
		self.width = 12
		self.height = 0
		self.startx = 0
		
		self.cylinder_around_z = 30
		
		self.arrow_button_height = 12
		self.arrow_part_offset = .5
		self.arrow_part_thickness = 1
		self.arrow_width = (self.width/2-self.arrow_part_offset)
		self.arrow_height = (self.arrow_button_height-2*self.arrow_part_offset)
		self.arrow_part_length = sqrt(self.arrow_width**2 +self.arrow_height**2)
		
		self.arrow_theta = 90-atan(self.arrow_height/self.arrow_width)*180.0/pi
		
		self.starty = 2*self.arrow_button_height
		
		self.target = weakref.ref(target)
		self.scroll_bit_position_ratio = 0
		self.scroll_bit_position = 0
		
		self.mouse_scroll_pos = None # used for moving the scroll bar with the mouse
		
		self.min_scroll_bar_size = 30
		
		self.scroll_bar_press_color = (.2,.2,.3,0)
		self.scroll_bar_idle_color = (0,0,0,0)
		self.scroll_bar_color = self.scroll_bar_idle_color
		
		self.scroll_bit_press_color = (0,0,.5,0)
		self.scroll_bit_idle_color = (.5,0,0,0)
		self.scroll_bit_color = self.scroll_bit_idle_color
		
		self.up_arrow_color = self.scroll_bar_idle_color
		self.down_arrow_color = self.scroll_bar_idle_color
		
		glShadeModel(GL_SMOOTH) # because we want smooth shading for the scroll bar components
		
	def update_stuff(self):
	
		view_height = self.target().gl_widget.height()
		current_y = self.target().origin[1]
		panel_height = self.target().matrix_panel.height
		adjusted_height = panel_height - view_height # adjusted height is the maximum value for current y!
		
		if adjusted_height <= 0:
			return False
		
		self.scroll_bar_height = self.height - self.starty
		self.scroll_bit_height = view_height/float(panel_height)*self.scroll_bar_height
		if self.scroll_bit_height < self.min_scroll_bar_size: self.scroll_bit_height = self.min_scroll_bar_size
		
		adjusted_scroll_bar_height = self.scroll_bar_height - self.scroll_bit_height
		
		self.scroll_bit_position_ratio = (-float(current_y)/adjusted_height)
		self.scroll_bit_position = self.scroll_bit_position_ratio *adjusted_scroll_bar_height
		
		return True
	def update_target_ypos(self):
		
		view_height = self.target().gl_widget.height()
		panel_height = self.target().matrix_panel.height
		adjusted_height = panel_height - view_height # adjusted height is the maximum value for current y!
	
		new_y = -self.scroll_bit_position_ratio*adjusted_height
		new_y = self.target().check_newy(new_y)
		
		old_x = self.target().origin[0]
		self.target().set_origin(old_x,new_y,False) # suppress gl updates because the calling functions always do this (at least, for now)
		
	def draw(self):
		
		if not self.update_stuff(): return # it doesn't make sense to draw the bar
		
		sx = self.startx
		sy = self.starty
		ex = self.width
		ey = self.height
		
		glNormal(0,0,1) # this normal is fine for everything that is drawn here
		# this provides some defaults for specular and shininess
		glMaterial(GL_FRONT,GL_SPECULAR,(.8,1,1,1.0))
		glMaterial(GL_FRONT,GL_SHININESS,20.0)

		# The scroll bar - the long bar that defines the extent
		glBegin(GL_QUADS)
		glNormal(0,0,1)
		x1 = self.startx
		x2 = self.width
		y1 = self.starty
		y2 = self.height
		glMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,self.scroll_bar_color )
		glVertex(x1,y1,0)
		glMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,(1,1,1,1.0))
		glVertex(x2,y1,0)
		glMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,(1,1,1,1.0))
		glVertex(x2,y2,0)
		glMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,self.scroll_bar_color )
		glVertex(x1,y2,0)
		glEnd()

		# The up pointing arrow
		glPushMatrix()
		glTranslate(0,self.arrow_button_height,0)
		glBegin(GL_TRIANGLES)
		glMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,self.up_arrow_color)
		glVertex(self.arrow_part_offset,self.arrow_part_offset,0)
		glMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,self.up_arrow_color)
		glVertex(self.width-self.arrow_part_offset,self.arrow_part_offset,0)
		glMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,(1,1,1,1.0))
		glVertex(self.width/2,self.arrow_button_height-self.arrow_part_offset,0)
		glEnd()
		glPopMatrix()
		
		# the down pointing arrow
		glBegin(GL_TRIANGLES) 
		glMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,self.down_arrow_color)
		glVertex(self.arrow_part_offset,self.arrow_button_height-self.arrow_part_offset,0)
		glMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,(1,1,1,1.0))
		glVertex(self.width/2,self.arrow_part_offset,0)
		glMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,self.down_arrow_color)
		glVertex(self.width-self.arrow_part_offset,self.arrow_button_height-self.arrow_part_offset,0)
		glEnd()
		
		# The scroll bit
		glBegin(GL_QUADS)
		x1 = sx
		x2 = ex
		y1 = sy+self.scroll_bit_position
		y2 = y1+self.scroll_bit_height
		glMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,(0,.5,0.5,1.0))
		glVertex(x1,y1,1)
		glMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,(1,1,1,1.0))
		glVertex(x2,y1,1)
		glMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,(1,1,1,1.0))
		glVertex(x2,y2,1)
		glMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,self.scroll_bit_color)
		glVertex(x1,y2,1)		
		glEnd()

	def down_button_pressed(self,color_change=True):
		if color_change: self.down_arrow_color = self.scroll_bar_press_color
		image_height = self.target().data.get_ysize()*self.target().scale+self.target().matrix_panel.min_sep
		newy = self.target().origin[1]+ image_height
		newy = self.target().check_newy(newy)
		oldx = self.target().origin[0]
		self.target().set_origin(oldx,newy)

	def up_button_pressed(self,color_change=True):
		if color_change: self.up_arrow_color = self.scroll_bar_press_color
		image_height = self.target().data.get_ysize()*self.target().scale+self.target().matrix_panel.min_sep
		newy  = self.target().origin[1] - image_height
		newy = self.target().check_newy(newy)
		oldx = self.target().origin[0]
		self.target().set_origin(oldx,newy)
	
	def scroll_down(self,color_change=True):
		if color_change: self.scroll_bar_color = self.scroll_bar_press_color
		newy  = self.target().origin[1]+ self.target().gl_widget.height()
		newy = self.target().check_newy(newy)
		oldx = self.target().origin[0]
		self.target().set_origin(oldx,newy)

	def scroll_up(self,color_change=True):
		if color_change: self.scroll_bar_color = self.scroll_bar_press_color
		newy  = self.target().origin[1] - self.target().gl_widget.height()
		newy = self.target().check_newy(newy)
		oldx = self.target().origin[0]
		self.target().set_origin(oldx,newy)
		
	def scroll_move(self,dy):
		ratio = dy/float(self.scroll_bar_height)
		panel_height = self.target().matrix_panel.height
		newy  = self.target().origin[1] - ratio*panel_height
		newy = self.target().check_newy(newy)
		oldx = self.target().origin[0]
		self.target().origin = (oldx,newy)
		self.target().updateGL()
		
	def mousePressEvent(self,event):
		x = self.target().gl_widget.width()-event.x()
		y = self.target().gl_widget.height()-event.y()
		if x < 0 or x > self.width: return # this shouldn't happen but it's nice to check I guess
		
		if y < self.arrow_button_height:
			self.down_button_pressed()
		elif y < 2*self.arrow_button_height:
			self.up_button_pressed()
		else:
			scroll_bit_ymin = self.starty +self.scroll_bit_position
			scroll_bit_ymax = scroll_bit_ymin + self.scroll_bit_height
			if y >= scroll_bit_ymin and y < scroll_bit_ymax:
				self.scroll_bit_color = self.scroll_bit_press_color
				self.mouse_scroll_pos = y
				self.target().updateGL()
			elif y < scroll_bit_ymin:
				self.scroll_down()
			else:
				self.scroll_up()
	
	def mouseMoveEvent(self,event):
		if self.mouse_scroll_pos != None:
			y = self.target().gl_widget.height()-event.y()
			dy = y - self.mouse_scroll_pos
			self.scroll_move(dy)
			self.mouse_scroll_pos = y
		
	def mouseReleaseEvent(self,event):
		self.mouse_scroll_pos = None
		self.scroll_bar_color = self.scroll_bar_idle_color
		self.scroll_bit_color = self.scroll_bit_idle_color
		self.up_arrow_color = self.scroll_bar_idle_color
		self.down_arrow_color = self.scroll_bar_idle_color
		self.target().updateGL()
		
class EMImageInspectorMX(QtGui.QWidget):
	def __init__(self,target,allow_col_variation=False,allow_window_variation=False,allow_opt_button=False):
		QtGui.QWidget.__init__(self,None)
		self.setWindowIcon(QtGui.QIcon(get_image_directory() +"multiple_images.png"))
		
		self.target=weakref.ref(target)
		self.busy = 1
		self.vals = QtGui.QMenu()
		self.valsbut = QtGui.QPushButton("Values")
		self.valsbut.setMenu(self.vals)
		
		try:
			self.vals.clear()
			vn=self.target().data.get_image_keys()
			vn.sort()
			for i in vn:
				action=self.vals.addAction(i)
				action.setCheckable(1)
				action.setChecked(0)
		except Exception, inst:
			print type(inst)     # the exception instance
			print inst.args      # arguments stored in .args
			print int
		
		action=self.vals.addAction("Img #")
		action.setCheckable(1)
		action.setChecked(1)
		
		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(2)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vboxlayout")
		
			
		self.hbl3 = QtGui.QHBoxLayout()
		self.hbl3.setMargin(0)
		self.hbl3.setSpacing(6)
		self.hbl3.setObjectName("hboxlayout")
		self.vbl.addLayout(self.hbl3)
		
		self.hist = ImgHistogram(self)
		self.hist.setObjectName("hist")
		self.hbl3.addWidget(self.hist)

		self.vbl2 = QtGui.QVBoxLayout()
		self.vbl2.setMargin(0)
		self.vbl2.setSpacing(6)
		self.vbl2.setObjectName("vboxlayout")
		self.hbl3.addLayout(self.vbl2)

		self.bsavedata = QtGui.QPushButton("Save")
		self.vbl2.addWidget(self.bsavedata)

		if allow_opt_button:
			self.opt_fit = QtGui.QPushButton("Opt. Fit")
			self.vbl2.addWidget(self.opt_fit)

		self.bsnapshot = QtGui.QPushButton("Snap")
		self.vbl2.addWidget(self.bsnapshot)
		if get_platform() != "Linux":
			self.bsnapshot.setEnabled(False)
			self.bsnapshot.setToolTip("Snap only works on Linux")

		# This shows the mouse mode buttons
		self.hbl2 = QtGui.QHBoxLayout()
		self.hbl2.setMargin(0)
		self.hbl2.setSpacing(6)
		self.hbl2.setObjectName("hboxlayout")
		self.vbl.addLayout(self.hbl2)

		self.mapp = QtGui.QPushButton("App")
		self.mapp.setCheckable(1)
		self.hbl2.addWidget(self.mapp)
		
		self.mdel = QtGui.QPushButton("Del")
		self.mdel.setCheckable(1)
		self.hbl2.addWidget(self.mdel)

		self.mdrag = QtGui.QPushButton("Drag")
		self.mdrag.setCheckable(1)
		self.hbl2.addWidget(self.mdrag)
		
		self.mset = QtGui.QPushButton("Sets")
		self.mset.setCheckable(1)
		self.hbl2.addWidget(self.mset)

		self.bg=QtGui.QButtonGroup()
		self.bg.setExclusive(1)
		self.bg.addButton(self.mapp)
		self.bg.addButton(self.mdel)
		self.bg.addButton(self.mdrag)
		self.bg.addButton(self.mset)
		
		self.mdrag.setChecked(True)
		
		
		self.hbl = QtGui.QHBoxLayout()
		self.hbl.setMargin(0)
		self.hbl.setSpacing(6)
		self.hbl.setObjectName("hboxlayout")
		self.vbl.addLayout(self.hbl)
		
		self.hbl.addWidget(self.valsbut)
		
		if self.target().using_ftgl():
			self.font_label = QtGui.QLabel("font size:")
			self.font_label.setAlignment(Qt.AlignRight|Qt.AlignVCenter)
			self.hbl.addWidget(self.font_label)
		
			self.font_size = QtGui.QSpinBox(self)
			self.font_size.setObjectName("nrow")
			self.font_size.setRange(1,50)
			self.font_size.setValue(int(self.target().get_font_size()))
			self.hbl.addWidget(self.font_size)
			
			QtCore.QObject.connect(self.font_size, QtCore.SIGNAL("valueChanged(int)"), self.target().set_font_size)
		
		
		self.banim = QtGui.QPushButton("Animate")
		self.banim.setCheckable(True)
		self.banim.setChecked(self.target().animation_enabled)
		self.hbl.addWidget(self.banim)
		
		self.tabwidget = QtGui.QTabWidget(self)
			
		self.tabwidget.addTab(self.get_image_manip_page(),"Main")
		self.tabwidget.addTab(self.get_sets_page(),"Sets")
		self.vbl.addWidget(self.tabwidget)

		
		# this stuff was used by flick... now it's here until I redo flick (for reference)
#		if allow_col_variation:
#			self.lbl2 = QtGui.QLabel("#/col:")
#			self.lbl2.setAlignment(Qt.AlignRight|Qt.AlignVCenter)
#			self.hbl.addWidget(self.lbl2)
#			
#			self.ncol = QtGui.QSpinBox(self)
#			self.ncol.setObjectName("ncol")
#			self.ncol.setRange(1,50)
#			self.ncol.setValue(self.target().get_rows())
#			self.hbl.addWidget(self.ncol)
#			
#		if allow_window_variation:
#			self.lbl3 = QtGui.QLabel("#/mx:")
#			self.lbl3.setAlignment(Qt.AlignRight|Qt.AlignVCenter)
#			self.hbl.addWidget(self.lbl3)
#			
#			self.nmx = QtGui.QSpinBox(self)
#			self.nmx.setObjectName("ncol")
#			self.nmx.setRange(1,50)
#			self.nmx.setValue(self.target().get_mxs())
#			self.hbl.addWidget(self.nmx)
		
		self.lowlim=0
		self.highlim=1.0
		self.del_mode_set = None # Because a selection set can have any name, the deletion mode set can be anyhting
		
		self.update_brightness_contrast()
		minden = self.target().get_density_min()
		maxden = self.target().get_density_max()
		self.hist.set_data(self.target().get_hist(),minden,maxden)
		self.busy=0
		
		QtCore.QObject.connect(self.vals, QtCore.SIGNAL("triggered(QAction*)"), self.newValDisp)
#		if allow_col_variation:
#			QtCore.QObject.connect(self.ncol, QtCore.SIGNAL("valueChanged(int)"), self.target().set_mx_rows)
#		if allow_window_variation:
#			QtCore.QObject.connect(self.nmx, QtCore.SIGNAL("valueChanged(int)"), self.target().set_mxs)
		
		#QtCore.QObject.connect(self.mmeas, QtCore.SIGNAL("clicked(bool)"), self.setMeasMode)
		QtCore.QObject.connect(self.mapp, QtCore.SIGNAL("clicked(bool)"), self.set_app_mode)
		QtCore.QObject.connect(self.mdel, QtCore.SIGNAL("clicked(bool)"), self.set_del_mode)
		QtCore.QObject.connect(self.mdrag, QtCore.SIGNAL("clicked(bool)"), self.set_drag_mode)
		QtCore.QObject.connect(self.mset, QtCore.SIGNAL("clicked(bool)"), self.set_set_mode)

		QtCore.QObject.connect(self.bsavedata, QtCore.SIGNAL("clicked(bool)"), self.save_data)
		if allow_opt_button:
			QtCore.QObject.connect(self.opt_fit, QtCore.SIGNAL("clicked(bool)"), self.target().optimize_fit)
		QtCore.QObject.connect(self.bsnapshot, QtCore.SIGNAL("clicked(bool)"), self.snapShot)
		QtCore.QObject.connect(self.banim, QtCore.SIGNAL("clicked(bool)"), self.animation_clicked)
		
	def animation_clicked(self,bool):
		self.target().animation_enabled = bool
	
	def get_desktop_hint(self):
		return "inspector"
	
	def get_image_manip_page(self):
		self.impage = QtGui.QWidget(self)
		vbl = QtGui.QVBoxLayout(self.impage )
		vbl.setMargin(2)
		vbl.setSpacing(6)
		vbl.setObjectName("get_image_manip_page")
		
		self.scale = ValSlider(self,(0.1,5.0),"Mag:")
		self.scale.setObjectName("scale")
		self.scale.setValue(self.target().scale)
		vbl.addWidget(self.scale)
		
		self.mins = ValSlider(self,label="Min:")
		self.mins.setObjectName("mins")
		vbl.addWidget(self.mins)
		minden = self.target().get_density_min()
		maxden = self.target().get_density_max()
		self.mins.setRange(minden,maxden)
		self.mins.setValue(minden)
		
		self.maxs = ValSlider(self,label="Max:")
		self.maxs.setObjectName("maxs")
		vbl.addWidget(self.maxs)
		self.maxs.setRange(minden,maxden)
		self.maxs.setValue(maxden)
		
		self.brts = ValSlider(self,label="Brt:")
		self.brts.setObjectName("brts")
		vbl.addWidget(self.brts)
		self.brts.setValue(0.0)
		self.brts.setRange(-1.0,1.0)
		
		self.conts = ValSlider(self,label="Cont:")
		self.conts.setObjectName("conts")
		vbl.addWidget(self.conts)
		self.conts.setValue(0.5)
		self.conts.setRange(-1.0,1.0)
		
		self.gammas = ValSlider(self,(.5,2.0),"Gam:")
		self.gammas.setObjectName("gamma")
		self.gammas.setValue(1.0)
		vbl.addWidget(self.gammas)

		QtCore.QObject.connect(self.scale, QtCore.SIGNAL("valueChanged"), self.target().set_scale)
		QtCore.QObject.connect(self.mins, QtCore.SIGNAL("valueChanged"), self.newMin)
		QtCore.QObject.connect(self.maxs, QtCore.SIGNAL("valueChanged"), self.newMax)
		QtCore.QObject.connect(self.brts, QtCore.SIGNAL("valueChanged"), self.newBrt)
		QtCore.QObject.connect(self.conts, QtCore.SIGNAL("valueChanged"), self.newCont)
		QtCore.QObject.connect(self.gammas, QtCore.SIGNAL("valueChanged"), self.newGamma)

		
		return self.impage
		
	def get_sets_page(self):
		self.setspage = QtGui.QWidget(self)
		hbl = QtGui.QHBoxLayout(self.setspage)
		self.setlist=QtGui.QListWidget(self)
		self.setlist.setSizePolicy(QtGui.QSizePolicy.Preferred,QtGui.QSizePolicy.Expanding)
		hbl.addWidget(self.setlist)
		
		available = self.target().data.sets.keys()
		first =  True
		for k in available:
			a = QtGui.QListWidgetItem("set"+str(k))
				
			flag2 = Qt.ItemFlags(Qt.ItemIsSelectable)
			flag3 = Qt.ItemFlags(Qt.ItemIsEnabled)
		  	flag4 = Qt.ItemFlags(Qt.ItemIsUserCheckable)
			
			a.setFlags(flag2|flag3|flag4)
			a.setCheckState(Qt.Unchecked)
			
			self.setlist.addItem(a)
			if first: 
				a.setSelected(True)
				first = False
			
			
		vbl = QtGui.QVBoxLayout()
		
		self.new_set_button = QtGui.QPushButton("New")
		vbl.addWidget(self.new_set_button)
		self.delete_set_button = QtGui.QPushButton("Delete")
		vbl.addWidget(self.delete_set_button)
		self.save_set_button = QtGui.QPushButton("Save")
		vbl.addWidget(self.save_set_button)
		
		hbl.addLayout(vbl)
		
		QtCore.QObject.connect(self.save_set_button, QtCore.SIGNAL("clicked(bool)"), self.save_set)
		QtCore.QObject.connect(self.new_set_button, QtCore.SIGNAL("clicked(bool)"), self.new_set)
		QtCore.QObject.connect(self.delete_set_button, QtCore.SIGNAL("clicked(bool)"), self.delete_set)
		QtCore.QObject.connect(self.setlist,QtCore.SIGNAL("itemChanged(QListWidgetItem*)"),self.set_list_item_changed)
		QtCore.QObject.connect(self.setlist,QtCore.SIGNAL("currentRowChanged(int)"),self.set_list_row_changed)
			
		return self.setspage
	
	def set_list_row_changed(self,i):
		if self.busy: return
		a = self.setlist.item(i)
		idx = a.index_key
		self.target().data.set_current_set(idx)
		
	def set_list_item_changed(self,item):
		checked = False
		if item.checkState() == Qt.Checked: checked = True
		
		i = item.index_key
				
		if checked:
			self.target().data.set_current_set(i)
			self.target().data.make_set_visible(i)
		else:
			self.target().data.make_set_not_visible(i)
			
		self.target().force_display_update()
		self.target().updateGL()
	
	def save_set(self,unused=None):
		selections = self.setlist.selectedItems()
		for i in selections:
			idx = i.index_key
			if self.target().data.sets.has_key(idx):
				s = self.target().data.sets[idx]
				db = db_open_dict("bdb:select")
				db[str(i.text())] = s
			else:
				print "there should be a warning message saying the set is empty"
	def delete_set(self,unused):
		self.busy = 1
		selections = self.setlist.selectedItems()
		
		for i in selections:
			idx = i.index_key
			self.target().data.delete_set(idx)
			db_name = "bdb:select#"+text
			if db_check_dict(db_name): db_remove_dict(db_name)
				
		items = self.get_set_list_items_copy()
		
		stext = [s.text() for s in selections]
		
		new_items = []
		for i in items:
			if i.text() not in stext:
				new_items.append(i)
				
		self.setlist.clear()
		for i in new_items: self.setlist.addItem(i)
#		if len(new_items) > 0:
#			new_items[-1].setSelected(True)
		
		self.busy = 0
		
		self.target().force_display_update()
		self.target().updateGL()
	
	def get_set_list_items_copy(self):
		items = [QtGui.QListWidgetItem(self.setlist.item(i)) for i in range(self.setlist.count())]
		return items

	def clear_sets(self):
		self.setlist.clear()

	def add_set(self,set_idx,set_name,display=True):
		items = self.get_set_list_items_copy()
		for item in items:
			if str(item.text()) == set_name:
				if item.checkState() == Qt.Checked:
					
					if display:
						self.target().data.set_current_set(set_idx)
						self.target().data.make_set_visible(set_idx)
				return
			 
		
		flag1 = Qt.ItemFlags(Qt.ItemIsEditable)
		flag2 = Qt.ItemFlags(Qt.ItemIsSelectable)
		flag3 = Qt.ItemFlags(Qt.ItemIsEnabled)
	  	flag4 = Qt.ItemFlags(Qt.ItemIsUserCheckable)
		a = QtGui.QListWidgetItem(set_name)
		a.setFlags(flag1|flag2|flag3|flag4)
		#a.setTextColor(qt_color_map[colortypes[parms[j][0]]])
		#if visible[j]:
		if display:	a.setCheckState(Qt.Checked)
		else: a.setCheckState(Qt.Unchecked)
		a.index_key = set_idx
	
	#a.setCheckState(Qt.Unchecked)
		
		self.setlist.addItem(a)
	
		a.setSelected(True)
		
		if display:
			self.target().data.set_current_set(set_idx)
			self.target().data.make_set_visible(set_idx)
		
		
	def new_set(self,unused=None):
		set = self.get_available_set_name()
		flag1 = Qt.ItemFlags(Qt.ItemIsEditable)
		flag2 = Qt.ItemFlags(Qt.ItemIsSelectable)
		flag3 = Qt.ItemFlags(Qt.ItemIsEnabled)
	  	flag4 = Qt.ItemFlags(Qt.ItemIsUserCheckable)
		a = QtGui.QListWidgetItem(set)
		a.setFlags(flag1|flag2|flag3|flag4)
		#a.setTextColor(qt_color_map[colortypes[parms[j][0]]])
		#if visible[j]:
		a.setCheckState(Qt.Checked)
		i = int(a.text()[-1])
		a.index_key = i
		
		#a.setCheckState(Qt.Unchecked)
			
		self.setlist.addItem(a) 
		a.setSelected(True)
		self.target().data.set_current_set(i)
		self.target().data.make_set_visible(i)
	
	def get_available_set_name(self):
		
		unavailable = []
		for i in range(self.setlist.count()):
			idx = self.setlist.item(i).index_key
			unavailable.append(idx)
		
		i = 0
		while True:
			if unavailable.count(i) == 0:
				return "set"+str(i)
			else:
				i += 1
			
	def set_del_mode(self,i):
		self.target().set_mouse_mode("del")
		
	def set_app_mode(self,i):
		self.target().set_mouse_mode("app")
	
	#def setMeasMode(self,i):
		#self.target.set_mouse_mode("meas")
	
	def set_drag_mode(self,i):
		self.target().set_mouse_mode("drag")
		
	def set_set_mode(self,unused=None):
		if self.busy: return
		self.busy = 1
		if not self.mset.isChecked():
			self.mset.setChecked(True)
		n = self.setlist.count()
		if n == 0:
			self.new_set() # enables the new set too
		else: # make sure a set is selected
			for i in range(n):
				a = self.setlist.item(i)
				if a.isSelected(): # something is selected, make sure its current in the data object
					a.setCheckState(Qt.Checked)
					i = a.index_key
					self.target().data.set_current_set(i)
					self.target().data.make_set_visible(i)
					break
			else:
				a = self.setlist.item(0)
				a.setCheckState(Qt.Checked)
				a.setSelected(True) # just make the first one the selected set
				self.target().data.set_current_set(0)
				self.target().data.make_set_visible(0)
					
			
		self.target().set_mouse_mode("set")
		self.tabwidget.setCurrentIndex(1) # This is the Sets page, but if someone could change the tab set up and break this behaviou
		self.busy = 0
	def set_scale(self,val):
		if self.busy : return
		self.busy=1
		self.scale.setValue(val)
		self.busy=0
		
	def set_n_cols(self,val):
		self.nrow.setValue(val)
		
	def set_n_rows(self,val):
		self.ncol.setValue(val)
		
	def set_mxs(self,val):
		self.nmx = val
	
	def save_data(self):
		self.target().save_data()
		
	def save_lst(self):
		self.target().save_lst()
			
	def snapShot(self):
		"Save a screenshot of the current image display"
		
		#try:
		qim=self.target().get_frame_buffer()
		#except:
			#QtGui.QMessageBox.warning ( self, "Framebuffer ?", "Could not read framebuffer")
		
		# Get the output filespec
		fsp=QtGui.QFileDialog.getSaveFileName(self, "Select File")
		fsp=str(fsp)
		
		qim.save(fsp,None,90)
		
	def newValDisp(self):
		v2d=[str(i.text()) for i in self.vals.actions() if i.isChecked()]
		self.target().set_display_values(v2d)

	def newMin(self,val):
		if self.busy : return
		self.busy=1
		self.target().set_density_min(val)

		self.update_brightness_contrast()
		self.busy=0
		
	def newMax(self,val):
		if self.busy : return
		self.busy=1
		self.target().set_density_max(val)
		self.update_brightness_contrast()
		self.busy=0
	
	def newBrt(self,val):
		if self.busy : return
		self.busy=1
		self.update_min_max()
		self.busy=0
		
	def newCont(self,val):
		if self.busy : return
		self.busy=1
		self.update_min_max()
		self.busy=0
	
	def newGamma(self,val):
		if self.busy : return
		self.busy=1
		self.target().set_gamma(val)
		self.busy=0

	def update_brightness_contrast(self):
		b=0.5*(self.mins.value+self.maxs.value-(self.lowlim+self.highlim))/((self.highlim-self.lowlim))
		c=(self.mins.value-self.maxs.value)/(2.0*(self.lowlim-self.highlim))
		self.brts.setValue(-b,1)
		self.conts.setValue(1.0-c,1)
		
	def update_min_max(self):
		x0=((self.lowlim+self.highlim)/2.0-(self.highlim-self.lowlim)*(1.0-self.conts.value)-self.brts.value*(self.highlim-self.lowlim))
		x1=((self.lowlim+self.highlim)/2.0+(self.highlim-self.lowlim)*(1.0-self.conts.value)-self.brts.value*(self.highlim-self.lowlim))
		self.mins.setValue(x0,1)
		self.maxs.setValue(x1,1)
		self.target().set_den_range(x0,x1)
		
	def set_hist(self,hist,minden,maxden):
		self.hist.set_data(hist,minden,maxden)

	def set_limits(self,lowlim,highlim,curmin,curmax):
		self.lowlim=lowlim
		self.highlim=highlim
		self.mins.setRange(lowlim,highlim)
		self.maxs.setRange(lowlim,highlim)
		self.mins.setValue(curmin,1)
		self.maxs.setValue(curmax,1)
		self.brts.setRange(-1.0,1.0)
		self.conts.setRange(0,1.0)
		self.update_brightness_contrast()


class EMDataListCache:
	'''
	This class designed primarily for memory management in the context of large lists of EMData objects. It is only efficient
	if accessing contiguous blocks of images - it is not effecient if randomly accessing images.
	
	The basic strategy is this - if asked for an image that is not currently in memory, then the cache is refocused about the newly
	requested image index.
	
	You can initialize this object with a filename of an image matrix. Then you can treat it as though it's
	a list of EMData objects. For example,
	
	data = EMDataListCache("particles.hdf")
	image1 = data[1] # just treat it like a list
	
	for i in data: i.write_image("test.hdf",-1) # iterable and enumerable
	
	'''
	LIST_MODE = 'list_mode'
	FILE_MODE = 'file_mode'
	def __init__(self,object,cache_size=256,start_idx=0):
		#DB = EMAN2db.EMAN2DB.open_db(".")
		self.xsize = -1
		self.ysize = -1
		self.keys = None
		self.current_set = None
		self.sets = {}
		self.set_init_flag = {} # sets stored on disk need an initialization flag, so action is only ever taken if the user chooses the set
		self.visible_sets = [] # stores the current sets, if there are more than one
		self.exclusions = []
		if isinstance(object,list):
			# in list mode there is no real caching
			self.mode = EMDataListCache.LIST_MODE
			self.max_idx = len(object)
			self.cache_size = self.max_idx
			self.images = object
			self.start_idx = 0

		elif isinstance(object,str):
			#print "file mode"
			self.mode = EMDataListCache.FILE_MODE
			if not os.path.exists(object) and not db_check_dict(object):
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
			if self.start_idx < 0: self.start_idx = 0
		
			self.__refresh_cache()
		else:
			print "the object used to construct the EMDataListCache is not a string (filename) or a list (of EMData objects). Can't proceed"
			return
		
		self.__init_sets()
		
		self.current_iter = 0 # For iteration support
		self.soft_delete = False # toggle to prevent permanent deletion of particles
		self.image_width = -1
		self.image_height = -1
	
	def __del__(self):
	   pass
	
	def get_xsize(self):
		if self.xsize == -1:
			
			if self.mode == EMDataListCache.FILE_MODE:
				for i in self.images:
					try:
						if self.images[i] != None:
							self.xsize = self.images[i].get_xsize()
							break
					except: pass
					
			elif self.mode == EMDataListCache.LIST_MODE:
				for i in self.images:
					try:
						self.xsize = i.get_xsize()
						break
					except: pass
					
					
		return self.xsize
	
	def get_ysize(self):
		if self.ysize == -1:
			
			if self.mode == EMDataListCache.FILE_MODE:
				for i in self.images:
					try:
						if self.images[i] != None:
							self.ysize = self.images[i].get_ysize()
							break
					except: pass
					
			elif self.mode == EMDataListCache.LIST_MODE:
				for i in self.images:
					try:
						self.ysize = i.get_ysize()
						break
					except: pass
		return self.ysize
	
	def get_image_header(self,idx):
		if self.mode == EMDataListCache.FILE_MODE:
			if len(self.file_name) > 3 and self.file_name[:4] == "bdb:":
				db = db_open_dict(self.file_name)
				return db.get_header(idx)
			else:
				e = EMData()
				e.read_image(self.file_name,idx,True)
				return e.get_attr_dict()
		elif self.mode == EMDataListCache.LIST_MODE:
			return self.images[idx].get_attr_dict()

	
	def get_image_keys(self):
		if self.keys == None:
			if self.mode == EMDataListCache.FILE_MODE:
				for i in self.images:
					try:
						if self.images[i] != None:
							self.keys = self.images[i].get_attr_dict().keys()
							break
					except: pass
					
			elif self.mode == EMDataListCache.LIST_MODE:
				for i in self.images:
					try:
						 self.keys = i.get_attr_dict().keys()
						 break
					except: pass
				
		return self.keys
	
	def __init_sets(self):
		if db_check_dict("bdb:select"):
			db = db_open_dict("bdb:select")
			keys = db.keys()
			for k in keys:
				try:
					idx = int(k[-1])
					self.sets[idx] = db[k]
					self.set_init_flag[idx] = True
				except: pass
	
	def make_set_visible(self,i):
		if i not in self.visible_sets:
			self.visible_sets.append(i)
		
	def set_current_set(self,i):
		if i == self.current_set: return # already this set
		self.current_set = i
		
		if not self.sets.has_key(i):
			self.sets[i] = []
			self.set_init_flag[i] = False
		elif self.set_init_flag[i]:
			for k in self.images.keys():
				if self.sets[i].count(k) != 0:
					im = self.images[k]
					if hasattr(im,"mxset"):
						im.mxset.append(i)
					else:
						im.mxset = [i]
			self.set_init_flag[i] = False
			
	def make_set_not_visible(self,i):
		if i in self.visible_sets: self.visible_sets.remove(i)
		if i == self.current_set:
			self.current_set = None
	
	def delete_set(self,idx):
		if self.sets.has_key(idx):
			self.sets.pop(idx)
		self.make_set_not_visible(idx)
		
		for im in self.images:
			if hasattr(im,"mxset") and idx in im.mxset: 
				print im.mxset
				im.mxset.remove(idx)
				if len(im.mxset) == 0: delattr(im,"mxset")

	def delete_box(self,idx):
		if self.mode == EMDataListCache.LIST_MODE and not self.soft_delete:
			# we can actually delete the emdata object
			self.exclusions.append(idx)
			image = self.images.pop(idx)
			self.max_idx = len(self.images)
			self.cache_size = self.max_idx
			return 1
		elif self.mode == EMDataListCache.FILE_MODE or self.soft_delete:
			if self.images[idx].has_attr("excluded"):
				self.images[idx].del_attr("excluded")
				self.exclusions.remove(idx)
			else:
				self.images[idx]["excluded"] = True
				self.exclusions.append(idx)
			return 2
		return 0
	
	def clear_sets(self):
		self.sets = {}
	
	def associate_set(self,idx,lst,mxset=False):
		self.sets[idx] = lst
		if mxset:
			for key,im in self.images.items():
				if key in lst:
					if hasattr(im,"mxset"):
						im.mxset.append(idx)
					else:
						im.mxset = [idx]
			
	
	def image_set_associate(self,idx):
		'''
		Associate the image with the set (or disassociate)
		'''
		if self.current_set == None: return 0 # there is no set, nothing happens
		
		im = self.images[idx]
		
		if hasattr(im,"mxset") and self.current_set in im.mxset:
			#print self.current_set,self.sets,im.mxset,idx 
			im.mxset.remove(self.current_set)
			if len(im.mxset) == 0: delattr(im,"mxset")
			self.sets[self.current_set].remove(idx)
		else:
			if hasattr(im,"mxset"):
				im.mxset.append(self.current_set)
			else:
				im.mxset = [self.current_set]
#				im.mxset = self.current_set
			self.sets[self.current_set].append(idx)
			
		return 1
				
	def save_lst(self,fsp):
		# Get the output filespec
		return # doesn't work yet
		f = file(fsp,'w')
		f.write('#LST\n')
		
		if self.mode == EMDataListCache.LIST_MODE and not self.soft_delete:
			raise # fixme - saving an lst in list mode is something that needs thought. EMAN2 only supports sets, anyway.
		elif self.mode == EMDataListCache.FILE_MODE or self.soft_delete:
			indices = [i for i in range(self.max_idx)]
			if self.sets.has_key(0):
				for exc in self.sets[0]: indices.remove(exc)
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
		''' Get the maximum image index  '''
		return self.max_idx
	
	def get_num_images(self):
		''' Get the number of images currently cached '''
		return len(self.images)
	
	def set_cache_size(self,cache_size,refresh=False):
		''' Set the cache size. May cause the cache to be refreshed, which could take a few moments '''
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
		''' Set the starting index of the cache, '''
		self.start_idx = start_idx
		if refresh: self.__refresh_cache()
	
	def __refresh_cache(self):
		app = QtGui.QApplication.instance()
		app.setOverrideCursor(Qt.BusyCursor)
		
		try:
			cache = {}
			i = self.start_idx
			for i in range(self.start_idx,self.start_idx+self.cache_size,1):
				if i != 0:
					idx = i % self.max_idx
				else: idx = 0
				try: 
					cache[idx] = self.images[idx]
				except:
					try:
						if self.mode ==  EMDataListCache.FILE_MODE:
							a = EMData()
							a.read_image(self.file_name,idx)
							if idx in self.exclusions: a["excluded"] = True
							cache[idx] = a
							if self.current_set != None:
								sets = []
								for set in self.current_set:
									
									if not hasattr(cache[idx],"mxset") and self.sets[set].count(idx) != 0:
										sets.append(set)
								if len(sets) != 0: cache[idx].mxset = sets
						else:
							print "data has been lost"
							raise
					except: print "couldn't access",idx,"the max idx was",self.max_idx,"i was",i,"start idx",self.start_idx,"cache size",self.cache_size,len(self.images)
				i += 1
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
			if self.start_idx < 0: self.start_idx = 0 
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
	def __len__(self):
		return self.max_idx
	
	def __iter__(self):
		''' Iteration support '''
		self.current_iter = 0
		return self
	
	def next(self):
		''' Iteration support '''
		if self.current_iter > self.max_idx:
			raise StopIteration
		else:
			self.current_iter += 1
			return self[self.current_iter-1]


if __name__ == '__main__':
	em_app = EMStandAloneApplication()
	window = EMImageMXModule(application=em_app)
	
	if len(sys.argv)==1 : 
		data = []
		for i in range(0,200):
			e = test_image(Util.get_irand(0,9))
			e.set_attr("excluded",Util.get_irand(0,1))
			data.append(e)
			
		window.set_data(data) 
	else :
		a=EMData.read_images(sys.argv[1])
		window.set_file_name(sys.argv[1])
		window.set_data(a)
		
	em_app.show()
	window.optimally_resize()
	em_app.execute()

