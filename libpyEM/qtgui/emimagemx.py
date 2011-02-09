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
from emimageutil import ImgHistogram, EMParentWin
from weakref import WeakKeyDictionary
from pickle import dumps,loads
from PyQt4.QtGui import QImage
from PyQt4.QtCore import QTimer
from libpyGLUtils2 import *

from emglobjects import EMOpenGLFlagsAndTools,EMGLProjectionViewMatrices,EMBasicOpenGLObjects,init_glut
from emapplication import EMGLWidget, get_application, EMApp
from emanimationutil import LineAnimation
import weakref

from emapplication import EMProgressDialog

class EMMXCoreMouseEvents:
	'''
	A base class for objects that handle mouse events in the EMImageMXWidget
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
	'''
	This class is like documentation. It doesn't really need to exist. MouseEventHandlers talk to this
	and not to the EMImageMXWidget. This talks to the EMImageMXWidget. So instead the MouseEventHandlers could
	just talk directly to EMImageMXWidget - I like it here because I can see exactly what's required of the
	MouseEventHandlers - however the problem could be solved using subclassing as a method of documentation
	It's much of a muchness.
	'''
	def __init__(self,target):
		if not isinstance(target,EMImageMXWidget):
			print "error, the target should be a EMImageMXWidget"
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
	
	def remove_particle_image(self,idx,event=None,redraw=False):
		self.target().remove_particle_image(idx,event,redraw)
		
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
				self.mediator.remove_particle_image(lc[0],event,True)
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
				
				# With shift-click we try to show rotated/translated images
				if  event.modifiers()&Qt.ShiftModifier:			# If shift is pressed, transform the particle orientations
					curim=[a["source_path"],a["source_n"]]		# this is the class-average
					if "EMAN2DB" in curim[0] : curim[0]="bdb:"+curim[0].replace("/EMAN2DB","")
					if curim[0][-10:-3]=="classes" : 
						mxpath=curim[0][:-11]+"#classmx_"+curim[0][-2:]
						try: mx=(EMData(mxpath,2),EMData(mxpath,3),EMData(mxpath,4),EMData(mxpath,5))
						except:
							mxpath=curim[0][:-11]+"#classify_"+curim[0][-2:]
							mx=(EMData(mxpath,2),EMData(mxpath,3),EMData(mxpath,4),EMData(mxpath,5))
					else: mx=None
				else: mx=None
				
				data = []
				idxs = d["class_ptcl_idxs"]
				name = d["class_ptcl_src"]
				progress = QtGui.QProgressDialog("Reading images from %s" %get_file_tag(name), "Cancel", 0, len(idxs),None)
				progress.show()
				get_application().setOverrideCursor(Qt.BusyCursor)
				
				i = 0
				for idx in idxs:
					data.append(EMData(name,idx))
					if mx!=None:
						xfm=Transform({"type":"2d","tx":mx[0][idx],"ty":mx[1][idx],"alpha":mx[2][idx],"mirror":bool(mx[3][idx])})
						data[-1].transform(xfm)
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
					self.class_window = EMImageMXWidget()
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
					self.mediator.remove_particle_image(lc[0],event,True)
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
			return False 
		  	# if you uncomment this code it will automatically set the scale in the main window so that the mxs stay visible
		  	# it's not what we wanted but it's left here in case anyone wants to experiment
#			scale = self.get_min_scale(view_width,view_height,view_scale,view_data)
#			target.scale = scale
#			view_scale = taget.scale
#			[self.ystart,self.visiblerows,self.visiblecols] = self.visible_row_col(view_width,view_height,view_scale,view_data,y)
		xsep = view_width - self.visiblecols*(rendered_image_width+self.min_sep)
		self.xoffset = xsep/2		
		
		self.height = ceil(len(view_data)/float(self.visiblecols))*(rendered_image_height+self.min_sep) + self.min_sep
		self.max_y = self.height - view_height # adjusted height is the maximum value for current y!
		return True
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
	
		

class EMImageMXWidget(EMGLWidget, EMGLProjectionViewMatrices):
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
		
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		#context = OpenGL.contextdata.getContext(None)
		#print "Matrix context is", context
		self.render()

	
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
		try: self.resize_event(width,height)
		except: pass
	
	def get_frame_buffer(self):
		# THIS WILL FAIL ON WINDOWS APPARENTLY, because Windows requires a temporary context - but the True flag is stopping the creation of a temporary context
		# (because display lists are involved)
		# December 2008 - tested in MAC it doesn't work,only returns a blank image
		return self.renderPixmap(0,0,True)

	def closeEvent(self,event):
		self.clear_gl_memory()
		EMGLWidget.closeEvent(self, event)		
		
	
	
	def enable_set(self,name,lst=[],display=True,update=True):
		'''
		Called from e2eulerxplor
		'''
		self.get_inspector()
		self.sets_manager.associate_set(name,lst,display)
		self.force_display_update()
		if update: self.updateGL()
		
	def clear_sets(self,update=True):
		self.sets_manager.clear_sets()
		self.force_display_update()
		if update: self.updateGL()
	
	def set_single_active_set(self,db_name):
		'''
		Called from emform
		'''
		self.get_inspector()
		self.sets_manager.clear_sets()
		self.sets_manager.enable_set(db_name)
		
		self.sets_manager.enable_sets_mode()
		self.force_display_update()
		self.updateGL()
		
	def load_font_renderer(self):
		try:
			self.font_render_mode = EMGLWidget.FTGL
			self.font_renderer = get_3d_font_renderer()
			self.font_renderer.set_face_size(self.font_size)
			self.font_renderer.set_font_mode(FTGLFontMode.TEXTURE)
		except:
			self.font_render_mode = EMGLWidget.GLUT
		
	def get_desktop_hint(self):
		return self.desktop_hint
	allim=WeakKeyDictionary()
	def __init__(self, data=None,application=None,winid=None, parent=None):
		
		self.emit_events = False

		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True)
		#fmt.setSampleBuffers(True)
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		EMGLProjectionViewMatrices.__init__(self)
		
		
		self.desktop_hint = "image"
		self.init_size_flag = True
		self.data=None
#		EMGLWidget.__init__(self,ensure_gl_context=True,winid=winid)
		EMGLWidget.__init__(self,winid=winid)
		EMImageMXWidget.allim[self] = 0
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
		self.mmode="Drag"
		self.selected=[]
		self.hist = []
		self.targetorigin=None
		self.targetspeed=20.0
		self.mag = 1.1				# magnification factor
		self.invmag = 1.0/self.mag	# inverse magnification factor
		self.glflags = EMOpenGLFlagsAndTools() 	# supplies power of two texturing flags
		self.tex_names = [] 		# tex_names stores texture handles which are no longer used, and must be deleted
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
		self.sets_manager = EMMXSetsManager(self) # we need this for managing sets
		self.deletion_manager = EMMXDeletionManager(self) # we need this for managing deleted particles
		
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
		
		self.reroute_delete = False
	
	def qt_parent_destroyed(self,object):
		self.under_qt_control = False
	
	def __del__(self):
		#self.clear_gl_memory() # this is intentionally commented out, it makes sense to clear the memory but not here
		if self.under_qt_control:
			self.deleteLater()
		self.core_object.deleteLater()
	def get_emit_signals_and_connections(self):
		return {"set_origin":self.set_origin,"set_scale":self.set_scale,"origin_update":self.origin_update}
	
	def get_data(self):
		'''
		Gets the current data object, this is either an EMDataListCache, an EM3DDataListCache, an EMLightWeightParticleCache, or None
		'''	
		return self.data
	
#	def width(self):
#		return self.view_width()
	
	def using_ftgl(self):
		return self.font_render_mode == EMGLWidget.FTGL
	
	def get_font_size(self):
		return self.font_renderer.get_face_size()
		
	def set_font_size(self,value):
		self.font_renderer.set_face_size(value)
		self.force_display_update() # only for redoing the fonts, this could be made more efficient :(
		self.updateGL()
			
	def __init_mouse_handlers(self):
		self.mouse_events_mediator = EMMXCoreMouseEventsMediator(self)
		self.mouse_event_handlers = {}
		self.mouse_event_handlers["App"] = EMMAppMouseEvents(self.mouse_events_mediator)
		self.mouse_event_handlers["Del"] = EMMXDelMouseEvents(self.mouse_events_mediator)
		self.mouse_event_handlers["Drag"] = EMMXDragMouseEvents(self.mouse_events_mediator)
		self.mouse_event_handlers[self.sets_manager.unique_name()] = EMMXSetMouseEvents(self.mouse_events_mediator)
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
			button = self.inspector.add_panel(self.sets_manager.get_panel_widget(),self.sets_manager.unique_name())
			self.sets_manager.set_button(button)
			self.inspector_update()
		return self.inspector
	
	def inspector_update(self):
		#FIXME
		pass
		#print "do inspector update"
		
	def set_reroute_delete(self,val=True):
		self.reroute_delete = val
	
	def get_scale(self):
		return self.scale
	
	def remove_particle_image(self,idx,event=None,update_gl=False):
		if self.reroute_delete == False:
			self.deletion_manager.delete_box(idx)
			if update_gl:
				self.force_display_update()
				self.updateGL()
				if event != None: self.emit(QtCore.SIGNAL("mx_boxdeleted"),event,[idx],False) 
		else:
			self.emit(QtCore.SIGNAL("mx_boxdeleted"),event,[idx],False)
	
	
	def image_set_associate(self,idx,event=None,update_gl=False):
		self.sets_manager.image_set_associate(idx)
		if update_gl:
			self.force_display_update()
			self.updateGL()

	def get_box_image(self,idx):
		return self.data[idx]

	def clear_gl_memory(self):
		self.makeCurrent() # this is important  when you have more than one OpenGL context operating at the same time
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
		self.resize(*self.get_parent_suggested_size())
			
		# I disabled this because it was giving me problems
#		self.scale =  self.matrix_panel.get_min_scale(self.view_width(),.height(),self.scale,self.data) # this is to prevent locking
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
	
	def setWindowTitle(self,filename):
		EMGLWidget.setWindowTitle(self, remove_directories_from_name(filename))
	
	def xyz_changed(self,xyz):
		if self.data.set_xyz(str(xyz)):
			self.display_states = []
			self.updateGL()
			
	def __get_cache(self,obj,soft_delete=False):
		'''
		Get the correct cache for the given obj
		@param obj a string or a list of EMData objects
		@param soft_delete only applicable if obj is a list
		'''
		if isinstance(obj,str):
			nx,ny,nz = gimme_image_dimensions3D(obj)
			if nz == 1:
				return EMLightWeightParticleCache.from_file(obj)
			else:
				return EM3DDataListCache(obj)
		else:
			if isinstance(obj,list):
				return EMDataListCache(obj,-1,soft_delete=soft_delete)
			else:
				return None
		
	def set_data(self,obj,filename='',update_gl=True,soft_delete=False):
		'''
		This function will work if you give it a list of EMData objects, or if you give it the file name (as the first argument)
		If this solution is undesirable one could easily split this function into two equivalents.
		'''
		cache_size = -1
		if isinstance(obj,EMMXDataCache): self.data = obj
		else: self.data = self.__get_cache(obj,soft_delete)
		
		if self.data == None: return
		
		if self.data.is_3d():
			self.get_inspector()
			self.inspector.enable_xyz()
			self.data.set_xyz(str(self.inspector.xyz.currentText()))
		else:
			self.get_inspector()
			self.inspector.disable_xyz()

		self.file_name = filename
		if self.file_name != None and len(self.file_name) > 0:self.setWindowTitle(self.file_name)
		
		self.force_display_update()
		self.nimg=len(self.data)
		self.max_idx = len(self.data)
		if self.nimg == 0: return # the list is empty
			
		global HOMEDB
		HOMEDB=EMAN2db.EMAN2DB.open_db()
		HOMEDB.open_dict("display_preferences")
		db = HOMEDB.display_preferences
		auto_contrast = db.get("display_stack_auto_contrast",dfl=True)
		start_guess = db.get("display_stack_np_for_auto",dfl=20)
		
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
			self.qt_parent.register_animatable(self.line_animation)
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
			
		if self.data and len(self.data)>0 and (self.data.get_ysize()*newscale>self.height() or self.data.get_xsize()*newscale>self.view_width()):
			newscale=min(float(self.height())/self.data.get_ysize(),float(self.view_width())/self.data.get_xsize())
			if self.inspector: self.inspector.scale.setValue(newscale)
		
		oldscale = self.scale
		self.scale=newscale
		self.matrix_panel.update_panel_params(self.view_width(),self.height(),self.scale,self.data,self.origin[1],self)
		
		
		view_height = self.height()
		panel_height = self.matrix_panel.height
		if panel_height < view_height :
			# This auto rescaling stuff was disabled at the requeset of Steve Ludtke. Uncomment to see how it works
			#if oldscale > newscale: self.scale =  self.matrix_panel.get_min_scale(self.view_width(),self.height(),self.scale,self.data) # this is to prevent locking
			self.draw_scroll = False
			self.origin=(self.matrix_panel.min_sep,self.matrix_panel.min_sep)
			self.matrix_panel.update_panel_params(self.view_width(),self.height(),self.scale,self.data,self.origin[1],self)
		else:
			self.draw_scroll = True
			self.scroll_bar.update_target_ypos()	
		
		if self.emit_events: self.emit(QtCore.SIGNAL("set_scale"),self.scale,adjust,update_gl)
		if update_gl: self.updateGL()
	
	def resize_event(self, width, height):
		self.scroll_bar.height = height
		
		self.matrix_panel.update_panel_params(self.view_width(),height,self.scale,self.data,self.origin[1],self)
		view_height = self.height()
		panel_height = self.matrix_panel.height
		if panel_height < view_heights :
			# This auto rescaling stuff was disabled at the requeset of Steve Ludtke. Uncomment to see how it works
			#self.scale =  self.matrix_panel.get_min_scale(self.view_width(),self.height(),self.scale,self.data) # this is to prevent locking
			self.draw_scroll = False
			self.origin=(self.matrix_panel.min_sep,self.matrix_panel.min_sep)
			self.matrix_panel.update_panel_params(self.view_width(),self.height(),self.scale,self.data,self.origin[1],self)
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
		# need to think this through yet...	
		#self.data.on_idle()
	
	def display_state_changed(self):
		display_states = []
		display_states.append(self.width())
		display_states.append(self.height())
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
		if self.font_render_mode != EMGLWidget.FTGL:
			print "error, can't call set_font_render_resolution if the mode isn't FTGL"
			
		#self.font_renderer.set_face_size(int(self.height()*0.015))
		#print "scale is",self.scale
	
	def __draw_backdrop(self):
		light = glIsEnabled(GL_LIGHTING)
		glDisable(GL_LIGHTING)
	
		glColor(.9,.9,.9)
		glBegin(GL_QUADS)
		glVertex(0,0,-1)
		glColor(.9,.9,.9)
		glVertex(self.width(),0,-1)
		glColor(.9,.9,.9)
		glVertex(self.width(),self.height(),-1)
		glColor(.9,.9,.9)
		glVertex(0,self.height(),-1)
		glEnd()
		if light: glEnable(GL_LIGHTING)
	
	
	def view_width(self):
		return EMGLWidget.width(self) - self.draw_scroll*self.scroll_bar.width
	
	def render(self):
		if not self.data : return
		if self.font_render_mode == EMGLWidget.FTGL: self.set_font_render_resolution()
		try: 
			self.image_change_count = self.data[0]["changecount"] 		# this is important when the user has more than one display instance of the same image, for instance in e2.py if 
		except:
			try: self.image_change_count = self.data["changecount"]
			except: pass

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
		
		self.matrix_panel.update_panel_params(self.view_width(),self.height(),self.scale,self.data,self.origin[1],self)
		
		if render: 
			deleted_idxs = self.deletion_manager.deleted_ptcls()
			visible_set_data = self.sets_manager.visible_set_data()
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

			if self.font_render_mode == EMGLWidget.GLUT:
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
			h=int(min(self.data.get_ysize()*self.scale,self.height()))
				
			invscale=1.0/self.scale
			self.set_label_ratio = 0.1
			self.coords = {}
			
			if self.matrix_panel.visiblerows:
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
						if (ty+th) > self.height():
							#print "ty + th was greater than",self.height()
							th = int(self.height()-ty)
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
						excluded = False
						try:
							(ii for ii in deleted_idxs if ii == i ).next()
							excluded = True
						except: pass
						if not excluded:
							#print rx,ry,tw,th,self.width(),self.height(),self.origin
							if not self.glflags.npt_textures_unsupported():
								if self.data[i]==None : print "bad image (tex) ",i
								a=GLUtil.render_amp8(self.data[i],rx,ry,tw,th,(tw-1)/4*4+4,self.scale,pixden[0],pixden[1],self.minden,self.maxden,self.gamma,2)
								self.texture(a,tx,ty,tw,th)
							else:
								if self.data[i]==None : print "bad image ",i								
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
						
						light = glIsEnabled(GL_LIGHTING)
						glEnable(GL_LIGHTING)
						iss = 0
						for set_name,set_list in visible_set_data.items():
							in_set = False
							try:
								(ii for ii in set_list  if ii == i ).next()
								in_set = True
							except: pass
							
							if in_set:
								x_pos_ratio = 1-self.set_label_ratio
								y_pos_ratio = x_pos_ratio
								y_pos_ratio -= 2*float(iss)/len(visible_set_data)*self.set_label_ratio
								y_pos_ratio *= 2
								x_pos_ratio *= 2
								self.sets_manager.load_set_color(set_name)
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
#									
							
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
		if self.font_render_mode == EMGLWidget.FTGL:
			
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
				if v=="Img #" and not self.data[i].has_attr("Img #") : 
					idx = i+self.img_num_offset
					if idx != 0: idx = idx%self.max_idx
					sidx = str(idx)
					bbox = self.font_renderer.bounding_box(sidx)
					GLUtil.mx_bbox(bbox,txtcol,bgcol)
					self.font_renderer.render_string(sidx)
					
				else : 
					try:
						av=self.data[i].get_attr(v)
						try:
							if isinstance(av,float) : avs="%1.4g"%av
							elif isinstance(av,Transform):
								t = av.get_rotation("eman")
								avs = "%1.1f,%1.1f,%1.1f" %(t["az"],t["alt"],t["phi"])
							elif isinstance(av,EMAN2Ctf):
								avs = "%1.3f,%1.1f" %(av.defocus,av.bfactor)
							else: avs=str(av)
						except:avs ="???"
					except: avs = ""
					bbox = self.font_renderer.bounding_box(avs)
					
					GLUtil.mx_bbox(bbox,txtcol,bgcol)
					self.font_renderer.render_string(avs)

				tagy+=self.font_renderer.get_face_size()

				glPopMatrix()
			if lighting:
				glEnable(GL_LIGHTING)
			glDisable(GL_TEXTURE_2D)
		elif self.font_render_mode == EMGLWidget.GLUT:
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
#		 print 'in render Text'
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
				self.targetorigin=(0,self.coords[n][1]-self.height()/2+self.data.get_ysize()*self.scale/2)
			except: return
		else:
			try: self.targetorigin=(self.coords[n][0]-self.view_width()/2+self.data.get_xsize()*self.scale/2,self.coords[n][1]-self.height()/2+self.data.get_ysize()*self.scale/2)
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
		
		absloc=((vec[0]),(self.height()-(vec[1])))
		for item in self.coords.items():
			index = item[0]+self.img_num_offset
			if index != 0: index %= self.max_idx
			data = item[1]
			if absloc[0]>data[0] and absloc[1]>data[1] and absloc[0]<data[0]+data[2] and absloc[1]<data[1]+data[3] :
				return (index,(absloc[0]-data[0])/self.scale,(absloc[1]-data[1])/self.scale, self.data[index].get_attr_dict())
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
		
		self.data.set_excluded_ptcls(self.deletion_manager.deleted_ptcls())
		from emsave import save_data
		file_name = save_data(self.data)
		if file_name == self.file_name and file_exists(file_name): # the file we are working with was overwritten
			self.set_data(file_name)
		self.data.set_excluded_ptcls(None)

	def save_lst(self,fsp):
		'''
		If we make it here the dialog has taken care of check whether or not overwrite should occur
		'''
		
		origname = self.get_image_file_name()

		f = file(fsp,'w')
		f.write('#LST\n')
		
		progress = EMProgressDialog("Writing files", "abort", 0, len(self.data),None)
		progress.show()
		for i in xrange(0,len(self.data)):
			d = self.data[i]
			if d == None: continue # the image has been excluded
			progress.setValue(i)
			get_application().processEvents()
			if progress.wasCanceled():
				#remove_file(fsp)# we could do this but if they're overwriting the original data then they lose it all
				f.close()
				progress.close()
				return
				#pass
		progress.close()
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
		if not self.data: return
		self.mouse_event_handler.mouse_double_click(event)

	def mousePressEvent(self, event):
		if not self.data: return
		if (self.width()-event.x() <= self.scroll_bar.width):
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
		if not self.data: return
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
			try:self.updateGL()
			except: pass
		else: 
			self.mouse_event_handler.mouse_move(event)
		
	def mouseReleaseEvent(self, event):
		if not self.data: return
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
		if not self.data: return
		if event.delta() > 0:
			self.set_scale( self.scale * self.mag )
		elif event.delta() < 0:
			self.set_scale(self.scale * self.invmag)
		#self.resize_event(self.width(),self.height())
		# The self.scale variable is updated now, so just update with that
		if self.inspector: self.inspector.set_scale(self.scale)
	
	def leaveEvent(self,event):
		get_application().setOverrideCursor(Qt.ArrowCursor)
		if self.mousedrag:
			self.mousedrag=None
			
	def get_frame_buffer(self):
		return self.get_frame_buffer()
	
	def draw_scroll_bar(self):
		width = self.width()
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
		
		
	def update_stuff(self):
	
		view_height = self.target().height()
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
		
		view_height = self.target().height()
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
		
		glShadeModel(GL_SMOOTH) # because we want smooth shading for the scroll bar components
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
		newy  = self.target().origin[1]+ self.target().height()
		newy = self.target().check_newy(newy)
		oldx = self.target().origin[0]
		self.target().set_origin(oldx,newy)

	def scroll_up(self,color_change=True):
		if color_change: self.scroll_bar_color = self.scroll_bar_press_color
		newy  = self.target().origin[1] - self.target().height()
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
		x = self.target().width()-event.x()
		y = self.target().height()-event.y()
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
			y = self.target().height()-event.y()
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
			vn=self.target().data.get_image_header_keys()
			vn.sort()
			for i in vn:
				action=self.vals.addAction(i)
				action.setCheckable(1)
				action.setChecked(0)
		except Exception, inst:
			print type(inst)	 # the exception instance
			print inst.args	  # arguments stored in .args
		
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

		#self.bnorm = QtGui.QPushButton("Norm")
		#self.bnorm.setCheckable(1)
		#self.vbl2.addWidget(self.bsnapshot)
		#self.bnorm.setToolTip("Normalize particles (in memory)")


		# This shows the mouse mode buttons
		self.hbl2 = QtGui.QHBoxLayout()
		self.hbl2.setMargin(0)
		self.hbl2.setSpacing(6)
		self.hbl2.setObjectName("hboxlayout")
		self.vbl.addLayout(self.hbl2)

		self.mapp = QtGui.QPushButton("App")
		self.mapp.setCheckable(1)
		self.hbl2.addWidget(self.mapp)
		
		self.mDel = QtGui.QPushButton("Del")
		self.mDel.setCheckable(1)
		self.hbl2.addWidget(self.mDel)

		self.mdrag = QtGui.QPushButton("Drag")
		self.mdrag.setCheckable(1)
		self.hbl2.addWidget(self.mdrag)
		
		
		self.mouse_mode_but_grp=QtGui.QButtonGroup()
		self.mouse_mode_but_grp.setExclusive(1)
		self.mouse_mode_but_grp.addButton(self.mapp)
		self.mouse_mode_but_grp.addButton(self.mDel)
		self.mouse_mode_but_grp.addButton(self.mdrag)
#		self.mouse_mode_but_grp.addButton(self.mset)
		
		self.mdrag.setChecked(True)
		
		
		self.hbl = QtGui.QHBoxLayout()
		self.hbl.setMargin(0)
		self.hbl.setSpacing(6)
		self.hbl.setObjectName("hboxlayout")
		self.vbl.addLayout(self.hbl)
		
		self.hbl.addWidget(self.valsbut)
		
		self.xyz = None
		if self.target().using_ftgl():
			self.font_label = QtGui.QLabel("font size:")
			self.font_label.setAlignment(Qt.AlignRight|Qt.AlignVCenter)
			self.hbl.addWidget(self.font_label)
		
			self.font_size = QtGui.QSpinBox()
			self.font_size.setObjectName("nrow")
			self.font_size.setRange(1,50)
			self.font_size.setValue(int(self.target().get_font_size()))
			self.hbl.addWidget(self.font_size)
			
			QtCore.QObject.connect(self.font_size, QtCore.SIGNAL("valueChanged(int)"), self.target().set_font_size)
		
		
		self.banim = QtGui.QPushButton("Animate")
		self.banim.setCheckable(True)
		self.banim.setChecked(self.target().animation_enabled)
		self.hbl.addWidget(self.banim)
		
		self.tabwidget = QtGui.QTabWidget()
			
		self.tabwidget.addTab(self.get_image_manip_page(),"Main")
		
		self.vbl.addWidget(self.tabwidget)
		
		self.lowlim=0
		self.highlim=1.0
		self.Del_mode_set = None # Because a selection set can have any name, the Deletion mode set can be anyhting
		
		self.update_brightness_contrast()
		minden = self.target().get_density_min()
		maxden = self.target().get_density_max()

		self.busy=0
		
		QtCore.QObject.connect(self.vals, QtCore.SIGNAL("triggered(QAction*)"), self.newValDisp)
#		QtCore.QObject.connect(self.mapp, QtCore.SIGNAL("clicked(bool)"), self.set_app_mode)
#		QtCore.QObject.connect(self.mDel, QtCore.SIGNAL("clicked(bool)"), self.set_Del_mode)
#		QtCore.QObject.connect(self.mdrag, QtCore.SIGNAL("clicked(bool)"), self.set_drag_mode)
#		QtCore.QObject.connect(self.mset, QtCore.SIGNAL("clicked(bool)"), self.set_set_mode)
		QtCore.QObject.connect(self.mouse_mode_but_grp,QtCore.SIGNAL("buttonClicked(QAbstractButton *)"),self.mouse_mode_button_clicked)

		QtCore.QObject.connect(self.bsavedata, QtCore.SIGNAL("clicked(bool)"), self.save_data)
		if allow_opt_button:
			QtCore.QObject.connect(self.opt_fit, QtCore.SIGNAL("clicked(bool)"), self.target().optimize_fit)
		QtCore.QObject.connect(self.bsnapshot, QtCore.SIGNAL("clicked(bool)"), self.snapShot)
		#QtCore.QObject.connect(self.bnorm, QtCore.SIGNAL("clicked(bool)"), self.setNorm)
		QtCore.QObject.connect(self.banim, QtCore.SIGNAL("clicked(bool)"), self.animation_clicked)
	
	def add_panel(self,widget,name):
		self.tabwidget.addTab(widget,name)
		
		
		button = QtGui.QPushButton(name)
		button.setCheckable(1)
		self.hbl2.addWidget(button)
		
		self.mouse_mode_but_grp.addButton(button)
		return button
	
	def set_mouse_mode(self,mode):
		b = (button for button in self.mouse_mode_but_grp.buttons() if str(button.text()) == mode ).next() # raises if it's not there, as it should
		b.setChecked(True) # triggers an event telling the EMImageMXWidget to changes its mouse event handler

	def set_current_tab(self,widget):
		self.tabwidget.setCurrentWidget(widget)

	def enable_xyz(self):
			
		if self.xyz == None:
			self.xyz = QtGui.QComboBox()
			self.xyz.addItems(["x","y","z"])
			self.hbl.addWidget(self.xyz)
			self.xyz.setCurrentIndex(2)
			QtCore.QObject.connect(self.xyz, QtCore.SIGNAL("currentIndexChanged(const QString&)"), self.target().xyz_changed)

	def disable_xyz(self):
		if self.xyz != None:
			self.hbl.removeWidget(self.xyz)
			self.xyz.deleteLater()
			self.xyz = None
	
	def animation_clicked(self,bool):
		self.target().animation_enabled = bool
	
	def get_desktop_hint(self):
		return "inspector"
	
	def get_image_manip_page(self):
		self.impage = QtGui.QWidget()
		vbl = QtGui.QVBoxLayout(self.impage )
		vbl.setMargin(2)
		vbl.setSpacing(6)
		vbl.setObjectName("get_image_manip_page")
		
		self.scale = ValSlider(None,(0.1,5.0),"Mag:")
		self.scale.setObjectName("scale")
		self.scale.setValue(self.target().scale)
		vbl.addWidget(self.scale)
		
		self.mins = ValSlider(None,label="Min:")
		self.mins.setObjectName("mins")
		vbl.addWidget(self.mins)
		minden = self.target().get_density_min()
		maxden = self.target().get_density_max()
		self.mins.setRange(minden,maxden)
		self.mins.setValue(minden)
		
		self.maxs = ValSlider(None,label="Max:")
		self.maxs.setObjectName("maxs")
		vbl.addWidget(self.maxs)
		self.maxs.setRange(minden,maxden)
		self.maxs.setValue(maxden)
		
		self.brts = ValSlider(None,label="Brt:")
		self.brts.setObjectName("brts")
		vbl.addWidget(self.brts)
		self.brts.setValue(0.0)
		self.brts.setRange(-1.0,1.0)
		
		self.conts = ValSlider(None,label="Cont:")
		self.conts.setObjectName("conts")
		vbl.addWidget(self.conts)
		self.conts.setValue(0.5)
		self.conts.setRange(-1.0,1.0)
		
		self.gammas = ValSlider(None,(.5,2.0),"Gam:")
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
	
	def mouse_mode_button_clicked(self,button):
		s = str(button.text())
		self.target().set_mouse_mode(s)
	
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
			
	#def setNorm(self,state):
		#"Set normalization mode"
		#self.target().set_norm(state)

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
		try: 
			b=0.5*(self.mins.value+self.maxs.value-(self.lowlim+self.highlim))/((self.highlim-self.lowlim))
			c=(self.mins.value-self.maxs.value)/(2.0*(self.lowlim-self.highlim))
		except:
			b=0
			c=0.5
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
		

class EMMXDeletionManager:
	'''
	This class handles everything to do with Deleting particles
	'''
	def __init__(self,target):
		self.target = weakref.ref(target)
		self.deleted_idxs = []
	def delete_box(self,idx):
		val = self.target().get_data().delete_box(idx)
		if val == 0: # value == 0 means responsibility is deferred to the EMMXDeletionManager
			try:
				#http://dev.ionous.net/2009/01/python-find-item-in-list.html
				self.deleted_idxs.remove(idx)
			except:	self.deleted_idxs.append(idx)
			
	def deleted_ptcls(self):
		return self.deleted_idxs

class EMMXSetsListWidgetItem(QtGui.QListWidgetItem):
	def __init__(self,*args):
		self.set_index = -1
		if len(args) > 0:
			if isinstance(args[0],EMMXSetsListWidgetItem):
				self.set_index = args[0].set_idx()
			
		QtGui.QListWidgetItem.__init__(self,*args)
		
	def set_idx(self): return self.set_index
	
	def assign_set_idx(self,val): self.set_index = val

class EMMXSetsPanel:
	'''
	Encapsulates the widgets that are used by the sets management system
	Talks directly to the EMMXSetsManager
	'''
	def __init__(self,target):
		self.target = weakref.ref(target) # this should be an EMMXSetsManager
		self.widget = None # this will be a qt widget
		self.busy = False
		
	def get_widget(self):
		self.busy = True
		if self.widget == None:
			
			self.widget = QtGui.QWidget()
			hbl = QtGui.QHBoxLayout(self.widget)
			self.setlist=QtGui.QListWidget()
			self.setlist.setSizePolicy(QtGui.QSizePolicy.Preferred,QtGui.QSizePolicy.Expanding)
			hbl.addWidget(self.setlist)
				
			vbl = QtGui.QVBoxLayout()
			
			self.new_set_button = QtGui.QPushButton("New")
			vbl.addWidget(self.new_set_button)
			self.delete_set_button = QtGui.QPushButton("delete")
			vbl.addWidget(self.delete_set_button)
			self.save_set_button = QtGui.QPushButton("Save")
			vbl.addWidget(self.save_set_button)
			
			hbl.addLayout(vbl)
			
			QtCore.QObject.connect(self.save_set_button, QtCore.SIGNAL("clicked(bool)"), self.save_set)
			QtCore.QObject.connect(self.new_set_button, QtCore.SIGNAL("clicked(bool)"), self.new_set)
			QtCore.QObject.connect(self.delete_set_button, QtCore.SIGNAL("clicked(bool)"), self.delete_set)
			QtCore.QObject.connect(self.setlist,QtCore.SIGNAL("itemChanged(QListWidgetItem*)"),self.set_list_item_changed)
			QtCore.QObject.connect(self.setlist,QtCore.SIGNAL("currentRowChanged(int)"),self.set_list_row_changed)
		
		self.busy = False
		return self.widget
	
	def set_list_row_changed(self,i):
		if self.busy: return
		a = self.setlist.item(i)
		#idx = a.set_idx()
		idx = str(a.text())
		self.target().set_current_set(idx)
		
	def set_list_item_changed(self,item):
		checked = False
		if item.checkState() == Qt.Checked: checked = True
		
		#i = item.set_idx()
		
		i = str(item.text())
				
		if checked:
			self.target().set_current_set(i)
			self.target().make_set_visible(i)
		else:
			self.target().make_set_invisible(i)
			
		self.target().request_display_update()
				
	def delete_set(self,unused):
		self.busy = True
		selections = self.setlist.selectedItems()
		
		for i in selections:
#			idx = i.set_idx()
			idx = str(i.text())
			self.target().delete_set(idx)
				
		items = self.get_set_list_items_copy()
		
		stext = [s.text() for s in selections]
		
		new_items = []
		for i in items:
			if i.text() not in stext:
				new_items.append(i)
				
		self.setlist.clear()
		for i in new_items: self.setlist.addItem(i)
		
		self.busy = False
		
		self.target().request_display_update()
		
		
	def new_set(self,unused=None):
		set = self.get_available_set_name()
		flag1 = Qt.ItemFlags(Qt.ItemIsEditable)
		flag2 = Qt.ItemFlags(Qt.ItemIsSelectable)
		flag3 = Qt.ItemFlags(Qt.ItemIsEnabled)
	  	flag4 = Qt.ItemFlags(Qt.ItemIsUserCheckable)
		a = EMMXSetsListWidgetItem(set)
		a.setFlags(flag1|flag2|flag3|flag4)
		color = self.target().get_set_color(set)
		a.setTextColor(self.to_qt_color(color))
		#a.setTextColor(qt_color_map[colortypes[parms[j][0]]])
		#if visible[j]:
		a.setCheckState(Qt.Checked)
#		i = int(a.text()[-1])
#		a.index_key = i
#		a.setProperty("index_key",i)
#		a.assign_set_idx(i)
		
		#a.setCheckState(Qt.Unchecked)
			
		self.setlist.addItem(a) 
		a.setSelected(True)
		self.target().set_current_set(str(a.text()))
		self.target().make_set_visible(str(a.text()))
	
		
	def get_set_list_items_copy(self):
		items = [EMMXSetsListWidgetItem(self.setlist.item(i)) for i in range(self.setlist.count())]
		
		return items
	
	def get_available_set_name(self):
		unavailable = []
		for i in range(self.setlist.count()):
		#	idx = self.setlist.item(i).set_idx()
			idx = str(self.setlist.item(i).text())
			unavailable.append(idx)
			
		
		i = 0
		while True:
			new_name = "Set"+str(i)
			if unavailable.count(new_name) == 0:
				return new_name
			else:
				i += 1

	def add_set(self,set_name,display=True):
		
		items = self.get_set_list_items_copy()
		for item in items:
			if str(item.text()) == set_name:
				if item.checkState() == Qt.Checked:
					if display:
						self.target().set_current_set(set_name)
						self.target().make_set_visible(set_name)
				return
			 
		flag1 = Qt.ItemFlags(Qt.ItemIsEditable)
		flag2 = Qt.ItemFlags(Qt.ItemIsSelectable)
		flag3 = Qt.ItemFlags(Qt.ItemIsEnabled)
	  	flag4 = Qt.ItemFlags(Qt.ItemIsUserCheckable)
		a = EMMXSetsListWidgetItem(set_name)
		a.setFlags(flag1|flag2|flag3|flag4)
		color = self.target().get_set_color(set_name)
		a.setTextColor(self.to_qt_color(color))
		if display:	a.setCheckState(Qt.Checked)
		else: a.setCheckState(Qt.Unchecked)
		
#		a.index_key = set_idx
#		a.setProperty("index_key", set_idx)
		#a.assign_set_idx(set_idx)
		
		
		self.setlist.addItem(a)
	
		a.setSelected(True)
		
		if display:
			self.target().set_current_set(set_name)
			self.target().make_set_visible(set_name)
	
	def to_qt_color(self,color):
		return QtGui.QColor(int(255*color[0]),int(255*color[1]),int(255*color[2]))
	
	def add_init_sets(self,set_names):
		
		for set_name in set_names:
			
			flag1 = Qt.ItemFlags(Qt.ItemIsEditable)
			flag2 = Qt.ItemFlags(Qt.ItemIsSelectable)
			flag3 = Qt.ItemFlags(Qt.ItemIsEnabled)
		  	flag4 = Qt.ItemFlags(Qt.ItemIsUserCheckable)
			a = EMMXSetsListWidgetItem(set_name)
			color = self.target().get_set_color(set_name)
			a.setTextColor(self.to_qt_color(color))
			a.setFlags(flag1|flag2|flag3|flag4)
			a.setCheckState(Qt.Unchecked)
#			a.index_key = set_name
#			a.setProperty("index_key", set_names)
			#a.assign_set_idx(set_name)
			
			self.setlist.addItem(a)
		
			a.setSelected(True)
			
	def save_set(self):
		selections = self.setlist.selectedItems()
		for item in selections:
			self.target().save_set(str(item.text()))

	def clear_sets(self):
		self.setlist.clear()
		
	def set_mode_enabled(self):
		n = self.setlist.count()
		if n == 0:
			self.new_set() # enables the new set too
		else: # make sure a set is selected
			for i in range(n):
				a = self.setlist.item(i)
				if a.isSelected(): # something is selected, make sure its current in the data object
					a.setCheckState(Qt.Checked)
#					i = a.set_idx()
					i = str(a.text())
					self.target().set_current_set(i)
					self.target().make_set_visible(i)
					break
			else:
				a = self.setlist.item(0)
				a.setCheckState(Qt.Checked)
				a.setSelected(True) # just make the first one the selected set
				self.target().set_current_set(str(a.text()))
				self.target().make_set_visible(str(a.text()))
					
class EMMXSetsManager:
	'''
	This is a class that separates the management of sets from the management of data
	'''
	
	def __init__(self,target):
		self.target = weakref.ref(target) # target should be an EMImageMXWidget
		# set stuff
		self.visible_sets = [] # this is a list of visible sets, ints
		self.sets = {} # this is a dictionary, keys are set numbers, values are list of indices of particles in the set 
		self.sets_color_map = {} # this is a dictionary mapping a set number (int) to a an int which is used for a color 
		self.keys = None
		self.current_set = None # an int represent the curret set
		self.panel = None # this will be the associated widget
		self.button = None # this will be a button in the inspector
		self.__init_sets()
	
	def get_set_color(self,set_name):
		if not self.sets_color_map.has_key(set_name):
			self.sets_color_map[set_name] = len(self.sets_color_map)
	
		return BoxingTools.get_color(self.sets_color_map[set_name]+1)
	
	def load_set_color(self,set_name):
		
		color = self.get_set_color(set_name)
		c = [color[0],color[1],color[2],1.0]
		glMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,c)
		glMaterial(GL_FRONT,GL_SPECULAR,c)
	  	glMaterial(GL_FRONT,GL_SHININESS,100.0)
	
	def set_button(self,button):
		self.button = button
	
	def add_set(self,set_name,display=True):
		if self.panel == None: self.__init_panel()
		self.panel.add_set(set_name,display)
	
	def unique_name(self): 
		'''
		A way to get a unique name for this object. If we end up adding more modules in the same
		fashion then they should have the same function, etc. For the time being this is the only
		object that is added in a modular fashion to the EMImageMXWidget
		'''
		return "Sets"
	
	def __init_panel(self):
		if self.panel == None:
			self.panel = EMMXSetsPanel(self)
			self.panel.get_widget()
			self.panel.add_init_sets(self.sets.keys())
		
	def get_panel_widget(self):
		'''
		Gets the panel widget for addition to the EMImageMXInspector
		'''
		if self.panel == None: self.__init_panel()
			
		return self.panel.get_widget()
	
	def __init_sets(self):
		'''
		Tries to read an previously establish sets from the local database
		'''
		if db_check_dict("bdb:select"):
			db = db_open_dict("bdb:select")
			for k,value in db.items():
				if isinstance(value,list):
					if len(value) > 0 and isinstance(value[0],int):
						self.sets[k] = value
						
	def request_display_update(self):
		'''
		A way for children objects to get the display to update
		'''
		self.target().force_display_update()
		self.target().updateGL()
		
	def get_available_sets(self,data_too=False):
		'''
		Gets the currently known sets
		'''
		return self.sets.keys()
	
	def has_set(self,set_name):
		'''
		Ask whether we have the given set in memory
		'''
		return self.sets.has_key(set_name)
	
	def get_set(self,set_name):
		'''
		Get the given se
		'''
		return self.sets[set_name]
	
	def visible_set_data(self):
		'''
		Get the visible sets as a dict
		'''
		d = {}
		for s in self.visible_sets:
			d[s] = self.sets[s]
		return d
	
	def make_set_visible(self,set_name):
		'''
		Make a particular set visible
		'''
		if set_name not in self.visible_sets:
			self.visible_sets.append(set_name)
#		
	def set_current_set(self,set_name):
		if set_name == self.current_set: return # already this set
		
		self.current_set = set_name
		
		if not self.sets.has_key(set_name):
			self.sets[set_name] = []
			
		if set_name not in self.visible_sets: self.visible_sets.append(set_name)

	def enable_set(self,set_name,display=True):
		act = True
		if db_check_dict("bdb:select"):
			db = db_open_dict("bdb:select")
			if db.has_key(set_name):
				self.sets[set_name] = db[set_name]
				act = False
		if act:
			self.sets[set_name] = []
			
		if self.panel == None: self.__init_panel()
		self.panel.add_set(set_name,display)
	
	def make_set_invisible(self,set_name):
		'''
		Remove a set from the list of visible sets
		'''
		if set_name in self.visible_sets: self.visible_sets.remove(set_name)
		if set_name == self.current_set: self.current_set = None
	
	def delete_set(self,set_name):
		if self.sets.has_key(set_name):
			self.sets.pop(set_name)
		else:
			raise RuntimeError("Error, attempt to delete a set that didn't exist")
		
		self.make_set_invisible(set_name)

	def clear_sets(self):
		'''
		A way to remove all information stored by this class - re-establish a blank slate
		'''
		self.sets = {}
		self.visible_sets = []
		if self.panel != None:
			self.panel.clear_sets()
			
		self.sets_color_map = {} 
	
	def associate_set(self,set_name,set_lst,display=True,force=False):
		'''
		@param set_name the name of the set, usually an int starting at 0
		@param set_list a list of particle indices indicating particles that are in the set
		@param force - needs to be specified if a set with the given set_name already exists
		'''
		if self.sets.has_key(set_name) and not force:
			raise RuntimeError("Cannot replace a set without supplying the force flag")
		
		self.sets[set_name] = set_lst
#			
		self.panel.add_set(set_name,display)
	
	def image_set_associate(self,idx):
		'''
		Associate/disassociate the image with self.current_set 
		@param idx particle index
		'''
		if self.current_set == None: return 0
		
		val = -1
		l = self.sets[self.current_set]
		try:
			#http://dev.ionous.net/2009/01/python-find-item-in-list.html
			val = (i for i in l if i == idx ).next()
		except:
			pass
		
		if val == -1: # it wasn't in the list - we add the particle to the list
			l.append(idx)
		else: # it is in the list so we remove it
			l.remove(idx)
			
		self.save_current_set()
		
	def enable_sets_mode(self,unused=None):
		inspector = self.target().get_inspector()
		inspector.set_mouse_mode(self.unique_name())
		self.target().set_mouse_mode(self.unique_name())
		inspector.set_current_tab(self.get_panel_widget())

		if self.panel == None: self.__init_panel()
		self.panel.set_mode_enabled()
	

	def save_current_set(self):
		if self.current_set == None: raise RuntimeError("Can't save current set, it doesn't exist")
		db = db_open_dict("bdb:select")
		db[self.current_set] = self.sets[self.current_set]
	
	def save_set(self,set_name):
		if not self.sets.has_key(set_name): raise RuntimeError("Unknown set %s" %set_name)
		db = db_open_dict("bdb:select")
		db[set_name] =self.sets[set_name]
		
	def save_sets(self):
		db = db_open_dict("bdb:select")
		for key,value in self.sets.items():
			db[key] = value
			
class EMMXDataCache:
	'''
	Base class for EMMXDataCaches
	'''
	def __init__(self):
		self.excluded_list = [] # a list of excluded idxs, used when saving the data to disk
		
	
	def delete_box(self,idx):
		'''
		@ must return a value = 1 indicates the box is permanently gone, 0 indicates the class is happy to do nothing
		and let the calling program display the deleted box differently
		'''
		raise NotImplementedException
	def __getitem__(self,idx): 
		'''
		operator[] must be supported
		'''
		raise NotImplementedException
	
	def __len__(self,idx): 
		'''
		must be able to get the length 
		'''
		raise NotImplementedException
	
	def get_xsize(self):
		'''
		must be able to get the xsize of the data
		'''
		raise NotImplementedException
	
	def get_ysize(self):
		'''
		must be able to get the zsize of the data
		'''
		raise NotImplementedException
	
	def get_zsize(self):
		'''
		must be able to get the zsize of the data
		'''
		raise NotImplementedException
	
	def is_complex(self):
		raise NotImplementedException
	
	def get_attr(self,attr):
		return self.get_image_header(0)[attr]
	
	def get_image_header_keys(self):
		'''
		Must be able to get the keys of in a typical EMData header, usually you just get the
		header of first image, cache it, and return it whenever it is asked for
		'''
		raise NotImplementedException
	
	def get_image_header(self,idx):
		'''
		Must be able to get the header of the image at the given index. Suggest reading only header
		from disk if the image is not already in memory
		'''
		raise NotImplementedException
	def on_idle(self):
		'''
		the EMImageMXWidget is at liberty to call this function when it becomes idle, etc.
		This function is useful for miserly caching-strategies, i.e. you can load images that
		are not already in memory
		'''
		raise NotImplementedException
	
	# CONCRETE FUNCTIONALITY	
	def set_excluded_ptcls(self,excluded_list):
		self.excluded_list = excluded_list
	
	
	def get_item_from_emsave(self,idx):
		try:
			(i for i in self.excluded_list if i == idx).next()
			return None
		except:
			return self[idx]
		
	def is_3d(self):
		'''
		Asks whether the cache is of type 3D - so the inspector can have an x/y/z combo box
		'''
		raise NotImplementedException
	
	def set_xyz(self,x_y_or_z):
		'''
		Must be supplied by 3d type caches
		@param x_y_or_z a string
		'''
		raise NotImplementedException

class ApplyTransform:
	def __init__(self,transform):
		self.transform = transform
	
	def __call__(self,emdata):
		emdata.transform(self.transform)
		
class ApplyAttribute:
	def __init__(self,attribute,value):
		self.attribute = attribute
		self.value = value
	
	def __call__(self,emdata):
		emdata.set_attr(self.attribute,self.value)

class ApplyProcessor:
	def __init__(self,processor="",processor_args={}):
		self.processor = processor
		self.processor_args = processor_args
		
	def __call__(self,data):
		data.process_inplace(self.processor,self.processor_args)

class EMLightWeightParticleCache(EMMXDataCache):
	'''
	A light weight particle cache is exactly that and more. Initialize it with a list of filenames and particle indices 
	corresponding to the particles that you want to view. So the calling function doesn't do any image reading, you just tell this
	thing what will need to be (or not) read.
	
	Primary data is basically a list like this: [[filename, idx, [func1,func2,...]], [filename, idx2, [func1,func2,...]],...]
	the filename and idx variables should be obvious, however the extra list contains functions that take an EMData as the argument -
	I used this, for example, for assignint attributes to images once they are in memory, and for transforming them, etc.
	
	A big advantage of this cache is that it only displays the images that are asked for. Additionally, it has a maximum cache size,
	and refocuses the cache when asked for an image outside its current index bounds. This makes this cache only suitable for linear access
	schemes, not random.
	
	'''
	def from_file(file_name):
		'''
		If this was C++ this would be the constructor for this class that took a single file name
		@param file_name the name of a particle stack file
		'''
		
		n = EMUtil.get_image_count(file_name)
		data = [[file_name,i,[]] for i in xrange(n)]
		
		return EMLightWeightParticleCache(data,len(data))
		
	from_file = staticmethod(from_file)
	
	def __init__(self,data,cache_max=2048):
		'''
		@param data list of lists - lists in in the list are of the form [image_name, idx, [list of functions that take an EMData as the first argument]]
		@param cache_max the maximum number of stored images - you might have to shorten this if you have very large images
		'''
		EMMXDataCache.__init__(self)
		self.data = data
		self.cache_start = 0
		self.cache_max = cache_max
		self.cache = [None for i in range(self.cache_max)]
		self.xsize = None
		self.ysize = None
		self.zsize = None
		self.header_keys = None
		# set stuff
		self.visible_sets = [] 
		self.sets = {}

	def __len__(self):
		'''
		support for len
		'''
		return len(self.data)
	
	def is_complex(self): return False
	
	def delete_box(self,idx):
		'''
		@ must return a value = 1 indicates the box is permanently gone, 0 indicates the class is happy to do nothing
		and let the calling program display the deleted box differently
		'''
		return 0
	
	def get_xsize(self):
		'''
		Get the get_xsize of the particles. Assumes all particle have the same size, which is potentially flawed
		'''
		if self.xsize == None:
			image = self[self.cache_start]
			self.xsize = image.get_xsize()
		
		return self.xsize
	
	def get_ysize(self):
		'''
		Get the get_ysize of the particles. Assumes all particle have the same size, which is potentially flawed
		'''
		if self.ysize == None:
			image = self[self.cache_start]
			self.ysize = image.get_ysize()
		
		return self.ysize
	
	def get_zsize(self):
		'''
		Get the get_ysize of the particles. Assumes all particle have the same size, which is potentially flawed
		'''
		if self.zsize == None:
			image = self[self.cache_start]
			self.zsize = image.get_zsize()
		
		return self.zsize
	
	def get_image_header(self,idx):
		'''
		Gets the header of the ith particle. Does not read the full image into memory if it's not stored, instead
		just reading the header and returning it. This can give significant speeds ups where only headers are needed,
		i.e. e2eulerxplor
		'''
#		return self[idx].get_attr_dict()
		adj_idx = idx-self.cache_start
		image = self.cache[adj_idx]
		if image == None:
			data = self.data[idx]
			h = get_header(data[0],data[1])
			return h
			#e.read_image(data[0],data[1],True)
			#return e.get_attr_dict()
		else:return image.get_attr_dict()
	
	def get_image_header_keys(self):
		'''
		Gets the keys in the header of the first image
		'''
		if self.header_keys == None:
			self.header_keys = self.get_image_header(self.cache_start).keys()
		return self.header_keys
	
	def refocus_cache(self,new_focus):
		'''
		Called internally to refocus the cache on a new focal point.
		@param new_focus the value at which the current cache failed - i.e. the first value that was beyond the current cache limits
		'''
		new_cache_start = new_focus-self.cache_max/2
		if new_cache_start < self.cache_start and (new_cache_start + self.cache_max) > self.cache_start:
			overlap = new_cache_start + self.cache_max - self.cache_start
			cache = [None for i in range(0,self.cache_max-overlap)]
			cache.extend(self.cache[0:overlap])
			self.cache = cache
		elif new_cache_start > self.cache_start and  new_cache_start < (self.cache_start + self.cache_max):
			cache = self.cache[new_cache_start:self.cache_start+self.cache_max]
			overlap =self.cache_max - ( new_cache_start - self.cache_start)
			cache = [None for i in range(0,self.cache_max-overlap)]
			cache.extend(self.cache[0:overlap])
			self.cache = cache
		else:
			self.cache = [None for i in range(self.cache_max)]
			
		self.cache_start = new_cache_start
			
	
	def __getitem__(self,idx):
		'''
		operator[] support - the main interface
		'''
		if idx < self.cache_start or idx > self.cache_start+self.cache_max:
			self.refocus_cache(idx)
			
		adj_idx = idx-self.cache_start
		image = self.cache[adj_idx]
		if image == None:
			a = self.__load_item(idx,adj_idx)
			return a
		else: return image
	
	def __load_item(self,idx,adj_idx):
		'''
		Work horse function for reading an image and applying any of the supplied functions 
		'''
		data = self.data[idx]
		
		try: 
			a = EMData(data[0],data[1])
			if a==None : raise Exception
		except :
			for i in range(10): 
				try:
					a=EMData(data[0],i)
					if a==None: raise Exception
				except: continue
				break
			a.to_zero()
		
		for func in data[2]: func(a) 
		self.cache[adj_idx] = a
		return a
	
	def on_idle(self):
		'''
		call this to load unloaded images in the cache
		This needs to be rethought, maybe use a thread?
		'''
		for idx,i in enumerate(self.cache):
			if i == None:
				# only does one at a time
				self.__load_item(idx,idx+self.cache_start)
				return

	def is_3d(self): return False
	
class EMDataListCache(EMMXDataCache):
	'''
	This class became semi-redundant after the introduction of the EMLightWeightParticleCache, however it is still used
	as the primary container for lists of the form [EMData,EMData,EMData,EMData,...]. 
	
	You can also initialize this object with a filename of an image matrix. Then you can treat it as though it's
	a list of EMData objects. For example,
	
	data = EMDataListCache("particles.hdf")
	image1 = data[1] # just treat it like a list
	
	for i in data: i.write_image("test.hdf",-1) # iterable and enumerable
	
	'''
	LIST_MODE = 'list_mode'
	FILE_MODE = 'file_mode'
	def __init__(self,object,cache_size=256,start_idx=0,soft_delete=False):
		EMMXDataCache.__init__(self)
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
		
#		self.__init_sets()
		
		self.current_iter = 0 # For iteration support
		self.soft_delete = soft_delete # toggle to prevent permanent Deletion of particles
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
	
	def get_zsize(self): return 1
	
	def is_complex(self): return False
	
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
	
	def get_image_header_keys(self):
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
	
	def on_idle(self):
		'''
		call this to load unloaded images in the cache
		'''
		for idx,i in enumerate(self.cache):
			if i == None:
				# only does one at a time
				self.__load_item(idx,idx+self.cache_start)
				return
	
	def delete_box(self,idx):
		'''
		@ must return a value = 1 indicates the box is permanently gone, 0 indicates the class is happy to do nothing
		and let the calling program display the deleted box differently
		'''
		if self.mode == EMDataListCache.LIST_MODE and not self.soft_delete:
			# we can actually delete the emdata object
			image = self.images.pop(idx)
			self.max_idx = len(self.images)
			self.cache_size = self.max_idx
			return 1
		elif self.mode == EMDataListCache.FILE_MODE or self.soft_delete:
			return 0
		
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
#							if idx in self.exclusions: a["excluded"] = True
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
	def on_idle(self):
		'''
		call this to load unloaded images in the cache, for example
		'''
		pass
	
	def is_3d(self): return False
	
class EM3DDataListCache(EMMXDataCache):
	'''
	A class that looks like a list to the outside word
	automated way of handling 3d images for the EMImageMXWidget
	
	'''
	def __init__(self,filename):
		EMMXDataCache.__init__(self)
		self.filename = filename
		self.nx,self.ny,self.nz = gimme_image_dimensions3D(filename)
		if self.nz == 1: raise RuntimeError("EM3DDataForMx class is only meant to be used with 3D images")
		self.keys = None
		self.header = None
		self.exclusions = []
		self.nz
		self.images = {}
		self.major_axis = "z"
		self.max_idx = self.nz
		
		
	def delete_box(self,idx):
		'''
		@ must return a value = 1 indicates the box is permanently gone, 0 indicates the class is happy to do nothing
		and let the calling program display the deleted box differently
		'''
		return 0
	
	def is_complex(self): return False
	
	def set_xyz(self,xyz):
		if xyz != self.major_axis:
			self.major_axis = xyz
			self.images = {}
			return True
		return False
	def get_xsize(self):
		return self.nx

	def get_ysize(self):
		return self.ny
	
	def get_zsize(self):
		return self.nz
	
	def get_image_header(self,idx):
		if self.header ==None:
			image = self[self.nz/2]
			self.header = image.get_attr_dict()
			
		return self.header

	
	def get_image_header_keys(self):
		if self.keys == None:
			self.keys = self[0].get_attr_dict().keys()
				
		return self.keys
	
	def get_max_idx(self):
		''' Get the maximum image index  '''
		return self.max_idx
	
	def get_num_images(self):
		''' Get the number of images currently cached '''
		return self.max_idx
	
	def set_cache_size(self,cache_size,refresh=False):
		''' Set the cache size. May cause the cache to be refreshed, which could take a few moments '''
		pass
#		if self.mode != EMDataListCache.LIST_MODE:
#			if cache_size > self.max_idx: self.cache_size = self.max_idx
#			else: self.cache_size = cache_size
#			self.start_idx = self.start_idx - self.cache_size/2
#			if refresh: self.__refresh_cache()
#		else:
#			if self.cache_size != self.max_idx:
#				print "error, in list mode the cache size is always equal to the max idx"
#				return
	def set_start_idx(self,start_idx,refresh=True):
		''' Set the starting index of the cache, '''
		pass
#		self.start_idx = start_idx
#		if refresh: self.__refresh_cache()
	
	def __refresh_cache(self):
		pass
	
	def __getitem__(self,idx):
		
		if not self.images.has_key(idx):
			a = EMData()
			if self.major_axis == "z":
				r = Region(0,0,idx,self.nx,self.ny,1)
				a.read_image(self.filename,0,False,r)
			elif self.major_axis == "y":
				r = Region(0,idx,0,self.nx,1,self.nz)
				a.read_image(self.filename,0,False,r)
				a.set_size(self.nx,self.nz)
			elif self.major_axis == "x":
				r = Region(idx,0,0,1,self.ny,self.nz)
				a.read_image(self.filename,0,False,r)
				a.set_size(self.ny,self.nz)
			
			self.images[idx] = a
			
		return self.images[idx]
			
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
		
	def on_idle(self):
		'''
		call this to load unloaded images in the cache, for example
		'''
		pass
	
	def is_3d(self): return True

class EMImageMXModule(EMImageMXWidget):
	def __init__(self, data=None,application=None,winid=None, parent=None):
		EMImageMXWidget.__init__(self, data, application, winid, parent)
		import warnings
		warnings.warn("convert EMImageMXModule to EMImageMXWidget", DeprecationWarning)

if __name__ == '__main__':
	em_app = EMApp()
	window = EMImageMXWidget(application=em_app)
	
	if len(sys.argv)==1 : 
		data = []
		for i in range(0,200):
			e = test_image(Util.get_irand(0,9))
			data.append(e)
			
		window.set_data(data,soft_delete=True) 
	else :
		a=EMData.read_images(sys.argv[1])
		window.set_file_name(sys.argv[1])
		window.set_data(a)
		
#	widget.show()
	em_app.show()
	window.optimally_resize()
	em_app.execute()
