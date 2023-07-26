#!/usr/bin/env python
#
# Author: Lan Dang, 03/17/2022 (dlan@bcm.edu)

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

from past.utils import old_div
from builtins import range
from PyQt5 import QtCore, QtGui, QtWidgets, QtOpenGL
from PyQt5.QtCore import Qt
import OpenGL
OpenGL.ERROR_CHECKING = False
from OpenGL import GL,GLU,GLUT
from OpenGL.GL import *
from eman2_gui.valslider import ValSlider,ValBox,StringBox,EMSpinWidget
from math import *
from EMAN2 import *

import EMAN2
import sys
import numpy as np
import struct
from eman2_gui.emimageutil import ImgHistogram, EMParentWin
from eman2_gui import emshape
from eman2_gui.emshape import EMShape
from weakref import WeakKeyDictionary
import weakref
from pickle import dumps,loads
from libpyGLUtils2 import *

from eman2_gui.emglobjects import EMOpenGLFlagsAndTools
from eman2_gui.emapplication import get_application, EMGLWidget
from eman2_gui.emimageutil import EMMetaDataTable
from eman2_gui.embrowser import EMBrowserWidget
from eman2_gui.empmwidgets import *
from eman2_gui.emanimationutil import SingleValueIncrementAnimation, LineAnimation
import json

import platform

from eman2_gui.emglobjects import EMOpenGLFlagsAndTools

class EMAnnotate2DWidget(EMGLWidget):
	"""
	"""
	origin_update = QtCore.pyqtSignal(tuple)
	signal_set_scale = QtCore.pyqtSignal(float)
	mousedown = QtCore.pyqtSignal(QtGui.QMouseEvent,tuple)
	mousedrag = QtCore.pyqtSignal(QtGui.QMouseEvent,tuple)
	mousemove = QtCore.pyqtSignal(QtGui.QMouseEvent,tuple)
	mouseup = QtCore.pyqtSignal(QtGui.QMouseEvent,tuple)
	mousewheel = QtCore.pyqtSignal(QtGui.QWheelEvent)
	signal_increment_list_data = QtCore.pyqtSignal(float)
	keypress = QtCore.pyqtSignal(QtGui.QKeyEvent)
	#itemdrop = QtCore.pyqtSignal(QtGui.QDropEvent)


	allim=WeakKeyDictionary()

	def __init__(self, image=None, annotation=None, application=get_application(),winid=None, parent=None,sizehint=(512,512)):

		self.inspector = None # this should be a qt widget, otherwise referred to as an inspector in eman
		print("Using local version pf emannotate widget")
		EMGLWidget.__init__(self,parent)
		self.setFocusPolicy(Qt.StrongFocus)
		self.setMouseTracking(True)
		self.initimageflag = True
		self.initsizehint = (sizehint[0],sizehint[1])	# this is used when no data has been set yet

		self.fftorigincenter = E2getappval("emimage2d","origincenter")
		if self.fftorigincenter == None : self.fftorigincenter=False
		emshape.pixelratio=self.devicePixelRatio()	# not optimal. Setting this factor globally, but should really be per-window

		self.data = None				# The currently displayed image/volume slice
		self.annotation = None			# The currently displayed annotation
		self.list_data = None
		self.list_idx = 0
		self.full_data = None 			# The actual image/volume
		self.full_annotation = None		# This is the actual annotation image/volume
		self.file_name = ""# stores the filename of the image, if None then member functions should be smart enough to handle it
		self.enable_clip = False
		EMAnnotate2DWidget.allim[self] = 0

		self.jsonfile=info_name(self.full_data)
		#info=js_open_dict(self.jsonfile)

		self.init_gl_flag = True
		self.oldsize=(-1,-1)
		self.scale=1.0				# Scale factor for display
		self.origin=(0,0)			# Current display origin
		self.invert=0				# invert image on display
		self.az = 0.0;				# used to transform slice orientation
		self.alt = 0.0;				# used to transform slice orientation
		self.zpos = 0				# z position relative to center
		self.minden=0
		self.maxden=1.0
		self.curmin=0.0
		self.curmax=0.0
		self.disp_proc=[] 			# a list/set of Processor objects to apply before rendering
		#self.tree_sels = []			# a list/item objects selected in the treeset to display with original color in rendering
		self.rmousedrag=None		# coordinates during a right-drag operation
		#self.mouse_mode_dict = {0:"emit", 1:"emit", 2:"emit", 3:"probe", 4:"measure", 5:"draw", 6:"emit", 7:"emit",8:"seg"}
		self.mouse_mode_dict = {0:"emit", 1:"emit", 2:"probe", 3:"measure", 4:"emit", 5:"emit",6:"seg"}

		self.mouse_mode = 6         # current mouse mode as selected by the inspector
		self.mag = 1.1				# magnification factor
		self.invmag = 1.0/self.mag	# inverse magnification factor

		self.shapes={}				# dictionary of shapes to draw, see add_shapes
		self.shapechange=1			# Set to 1 when shapes need to be redrawn
		self.active=(None,0,0,0)	# The active shape and a highlight color (n,r,g,b)

		self.extras = []			# an empty set of extras - other images that can be rendered over this one

		self.startorigin = None
		self.endorigin = None
		self.isanimated = False
		self.time = 1
		self.timeinc = 0.125
		self.key_mvt_animation = None

		self.current_class = 1

		self.init_size = True		# A flag used to set the initial origin offset

		self.shapelist = 0			# a display list identify

		self.glflags = EMOpenGLFlagsAndTools() 	# supplies power of two texturing flags

		self.tex_name = 0			# an OpenGL texture handle

		self.window_width = None # Used for intelligently managing resize events
		self.window_height = None # Used for intelligently managing resize events
		#self.setSizePolicy(QtWidgets.QSizePolicy.Preferred,QtWidgets.QSizePolicy.Expanding)


		self.init_size_flag = True
		self.frozen = False
		self.isexcluded = False
		self.parent_geometry = None

		self.eraser_shape = None # a single circle shape used 0for erasing in e2boxer


		self.use_display_list = True # whether or not a display list should be used to render the image pixelsw - if on, this will save on time if the view of the image is unchanged, which can quite often be the case
		self.main_display_list = -1	# if using display lists, the stores the display list
		self.display_states = [] # if using display lists, this stores the states that are checked, and if different, will cause regeneration of the display list
		self.hist = []

		self.display_shapes = True # A flag that can be used to turn of the display of shapes - useful to e2boxer
		self.display_z_label = True
		self.circle_dl = None # used for a circle list, for displaying circled particles, for example
		#self.xform = Transform({"type":"eman","alt":self.alt,"az":self.az,"tx":self.full_data["nx"]//2,"ty":self.full_data["ny"]//2,"tz":self.full_data["nz"]//2+self.zpos})
		self.xform = Transform()
		self.display_group = False

		s=np.arange(0, 15, 0.5)
		#print(s.type)
		sf=XYData()
		sf.set_xy_list(s.tolist(),s.tolist())

		#self.ctable = s.tolist()
		self.need_new_RGB = 1
		self.colors = self.create_palette(256)
		self.ctable = self.create_RGB_list()


		#self.external_seg_tab = None

		self.setAcceptDrops(True) #TODO: figure out the purpose of this (moved) line of code
		self.setWindowIcon(QtGui.QIcon(get_image_directory() +"single_image.png")) #TODO: figure out why this icon doesn't work

		if image :
			if annotation:
				self.set_data(image, annotation)
			else:
				self.set_data(image, None)
		self.auto_contrast(inspector_update=True,display_update=True)
		#print("Full data nz",self.full_data["nz"])


			#print("data", (self.data))
			#print("ann", (self.annotation))
			#print("full data", (self.full_data))
			#print("full ann", (self.full_annotation))

	def get_annotation(self) :
		if self.annotation is None:
			print("Annotation is", self.annotation)
		return self.annotation
		#return self.full_annotation
	def get_full_annotation(self):
		return self.full_annotation

	def initializeGL(self):
		GL.glClearColor(0,0,0,0)

		glLightfv(GL_LIGHT0, GL_AMBIENT, [0.1, 0.1, 0.1, 1.0])
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_POSITION,  [.1,.1,1,0.])

		glEnable(GL_LIGHTING)
		glEnable(GL_LIGHT0)

	def paintGL(self):
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		# glClear(GL_COLOR_BUFFER_BIT) # throws error.
		glClearColor(0.0, 0.0, 0.0, 0.0)
		if glIsEnabled(GL_DEPTH_TEST):
			glClear(GL_DEPTH_BUFFER_BIT)
		if glIsEnabled(GL_STENCIL_TEST):
			glClear(GL_STENCIL_BUFFER_BIT)
		#self.cam.position()
		#context = OpenGL.contextdata.getContext(None)
		#print "Image2D context is", context
		glPushMatrix()
		self.render()
		glPopMatrix()

	def resizeGL(self, width, height):
		if width == 0 or height == 0: return # this is okay, nothing needs to be drawn
		width = width // self.devicePixelRatio()
		height = height // self.devicePixelRatio()
		side = min(width, height)
		dpr=self.devicePixelRatio()
		GL.glViewport(0,0,self.width()*dpr,self.height()*dpr)

		GL.glMatrixMode(GL.GL_PROJECTION)
		GL.glLoadIdentity()
		GLU.gluOrtho2D(0.0,width,0.0,height)
		GL.glMatrixMode(GL.GL_MODELVIEW)
		GL.glLoadIdentity()

		self.resize_event(width,height)
		#except: pass

	def optimally_resize(self):
		if self.parent_geometry != None:
			#self.load_default_scale_origin()
			self.restoreGeometry(self.parent_geometry)
		else:
			new_size = self.get_parent_suggested_size()
			self.resize(*new_size)
			self.load_default_scale_origin(new_size)

	def get_parent_suggested_size(self):


		if self.data==None : return (self.initsizehint[0]+12,self.initsizehint[1]+12)
		data = self.data

		try: return (data["nx"]+12,data["ny"]+12)
		except : return (self.initsizehint[0]+12,self.initsizehint[1]+12)

	def sizeHint(self):
#		print self.get_parent_suggested_size()
		if self.data==None: return QtCore.QSize(*self.initsizehint)
		return QtCore.QSize(*self.get_parent_suggested_size())

	def set_disp_proc(self,procs):
		print('set_disp_proc')
		self.disp_proc=procs
		self.force_display_update(set_clip=True)
		self.updateGL()




	def set_enable_clip(self,val=True):
		self.enable_clip = val

	def clear_gl_memory(self):
		self.makeCurrent() # this is important  when you have more than one OpenGL context operating at the same time
		if (self.shapelist != 0):
			glDeleteLists(self.shapelist,1)
			self.shapelist = 0
		if self.main_display_list != -1:
			glDeleteLists(self.main_display_list,1)
			self.main_display_list = -1


	def set_mouse_mode(self,mode_num):
		if self.mouse_mode == mode_num:
			return
		if mode_num not in self.mouse_mode_dict:
			print("unknown mouse mode:",mode_num)
			return
		self.mouse_mode = mode_num
		print("Mouse mode is: ",mode_num,"which is", self.mouse_mode_dict[self.mouse_mode])
		self.del_shapes()

	def set_xform(self,tx,ty,tz,alt,az):
		self.xform = Transform({"type":"eman","alt":alt,"az":az,"tx":tx,"ty":ty,"tz":tz})
		#print("New xform", self.xform)
		#self.signal_set_scale.emit(1)
		self.inspector_update()
		self.force_display_update()
		self.updateGL()



	def get_minden(self): return self.minden
	def get_maxden(self): return self.maxden
	def get_alt(self): return self.alt
	def get_az(self): return self.az
	def get_shapes(self): return self.shapes
	def get_xform(self): return self.xform

	def coords_within_image_bounds(self,coords):
		x = coords[0]
		y = coords[1]
		if x < 0 or y < 0:
			return False

		if x >= self.data.get_xsize() or y >= self.data.get_ysize():
			return False

		return True

	def set_density_range(self,x0,x1):
		"""Set the range of densities to be mapped to the 0-255 pixel value range"""
		self.curmin=x0
		self.curmax=x1
		self.force_display_update()
		self.updateGL()

	def set_density_min(self,val):
		self.curmin=val
		self.force_display_update()
		self.updateGL()

	def set_density_max(self,val):
		self.curmax=val
		self.force_display_update()
		self.updateGL()

	def set_alt(self,val):
		#self.set_clip()
		self.alt=val
		#self.set_xform(self.nx//2,self.ny//2,self.nz//2+self.zpos,self.alt,self.az)
		self.force_display_update(set_clip=True)
		self.updateGL()

	def set_az(self,val):
		#self.set_clip()
		self.az=val
		#self.set_xform(self.nx//2,self.ny//2,self.nz//2+self.zpos,self.alt,self.az)
		self.force_display_update(set_clip=True)
		self.updateGL()

	def set_n(self,val):
		#self.set_clip()
		self.zpos = val
		#self.set_xform(self.nx//2,self.ny//2,self.nz//2+self.zpos,self.alt,self.az)
		self.force_display_update(set_clip=True)
		self.updateGL()

	def set_file_name(self,file_name):
		self.file_name = file_name

	def get_file_name(self):
		return self.file_name

	def get_data_dims(self):
		data = None

		if self.data != None: data = self.data
		else: return [0,0,0]

		return [data.get_xsize(),data.get_ysize(),data.get_zsize()]

	# def updateGL(self):
	# 	if self.gl_widget != None:
	# 		self.gl_widget.updateGL()

	def set_frozen(self,frozen):
		self.frozen = frozen

	def set_excluded(self,isexcluded):
		wasexcluded = self.isexcluded

		self.isexcluded = isexcluded

		if wasexcluded or self.isexcluded: return True
		else: return False

	def get_data(self):

		#return self.full_data
		if self.data is None:
			self.render_bitmap()
			print("Data is None. Render")
		return self.data

	def get_full_data(self):
		return self.full_data

	def set_data(self,data,annotation=None,file_name="",retain_current_settings=True, keepcontrast=False):
		"""You may pass a single 2D image or a single 3-D volume
		incoming_annotation must have the same dimensionality as incoming data or be None"""

		self.set_file_name(file_name)
		if self.file_name != "": self.setWindowTitle(remove_directories_from_name(self.file_name))

		if data == None:
			self.data=self.full_data = None
			self.annotation=self.full_annotation=None
			return

		# make an empty annotation if none provided
		if annotation is None :
			annotation=data.copy_head()
			annotation.to_zero()

		needresize=True if self.data==None else False

		apix=data["apix_x"]
		if not isinstance(data,EMData):
			raise Exception("EMAnnotate2D only supports a single 2-D or 3-D image")

		if data.get_sizes()!=annotation.get_sizes() :
			print("EMAnnotate2D requires data and annotation to be the same size")
			raise Exception("EMAnnotate2D requires data and annotation to be the same size")

		self.full_data=data
		self.full_annotation=annotation
		self.nx = self.full_annotation["nx"]
		self.ny = self.full_annotation["ny"]
		self.nz = self.full_annotation["nz"]

		#TOREAD
		#self.xform = Transform({"type":"eman","alt":self.alt,"az":self.az,"tx":self.full_data["nx"]//2,"ty":self.full_data["ny"]//2,"tz":self.full_data["nz"]//2+self.zpos})
		self.set_xform(self.nx//2,self.ny//2,self.nz//2+self.zpos,self.alt,self.az)
		#self.data = self.full_data.copy_head()
		#self.annotation = self.full_annotation.copy_head()
		#self.data=self.full_data.get_clip(Region(0,0,0,self.full_data["nx"],self.full_data["ny"],1))
		#self.annotation=self.full_annotation.get_clip(Region(0,0,0,self.full_data["nx"],self.full_data["ny"],1))
		#print("Set initial data")



		self.image_change_count = 0

		if not keepcontrast:
			self.auto_contrast(inspector_update=False,display_update=False)
		try:
			if needresize:
				x=self.data["nx"]
				y=self.data["ny"]
				xys=QtWidgets.QApplication.desktop().availableGeometry()
				mx=xys.width()*2//3
				my=xys.height()*2//3

				self.resize(min(x,mx),min(y,my))
		except: pass


		self.inspector_update()
		self.force_display_update()
		self.updateGL()

	def load_default_scale_origin(self,size_specified=None):
		'''
		Size specified is just a way of saying that the height is not
		currently correct, for example. It's used in self.optimally_resize
		'''
		self.scale=1.0				# self.origin=(0,0)Scale factor for display
		self.origin = (0,0)
		##try:
		if size_specified:
			w = size_specified[0]
			h = size_specified[1]
		else:
			w = self.width()
			h = self.height()
		data = self.get_data_dims()
		if data[0] == 0 or data[1] == 0:
			self.scale=1.0
			return
		scalew = old_div(float(w),data[0])
		scaleh = old_div(float(h),data[1])
		if scaleh < scalew:
			self.scale = scaleh
		else: self.scale = scalew


		#except: pass

	def full_contrast(self,boolv=False,inspector_update=True,display_update=True):

		if self.data == None: return
		# histogram is impacted by downsampling, so we need to compensate
		if self.scale<=0.5 :
			down=self.data.process("math.meanshrink",{"n":int(floor(1.0/self.scale))})
			mean=down["mean"]
			sigma=down["sigma"]
			m0=down["minimum"]
			m1=down["maximum"]
		else:

			mean=self.data.get_attr("mean")
			sigma=self.data.get_attr("sigma")
			m0=self.data.get_attr("minimum")
			m1=self.data.get_attr("maximum")
		self.minden=m0
		self.maxden=m1
		self.curmin = m0
		self.curmax = m1
		if inspector_update: self.inspector_update()
		if display_update:
			self.force_display_update()
			self.updateGL()

	def auto_contrast(self,boolv=False,inspector_update=True,display_update=True):
		auto_contrast = E2getappval("display2d","autocontrast",True)

		if self.full_data == None:
			#print("data is None")
			return
		# histogram is impacted by downsampling, so we need to compensate
		if self.scale<=0.5 :
			down=self.full_data.process("math.meanshrink",{"n":int(floor(1.0/self.scale))})
			mean=down["mean"]
			sigma=down["sigma"]
			m0=down["minimum"]
			m1=down["maximum"]
		else:
			mean=self.full_data.get_attr("mean")
			sigma=self.full_data.get_attr("sigma")
			m0=self.full_data.get_attr("minimum")
			m1=self.full_data.get_attr("maximum")
		self.minden=m0
		self.maxden=m1
		if auto_contrast:
			self.curmin = max(m0,mean-3.0*sigma)
			self.curmax = min(m1,mean+3.0*sigma)
		else:
			self.curmin = m0
			self.curmax = m1
		if inspector_update: self.inspector_update()
		if display_update:
			self.force_display_update()
			self.updateGL()

	def set_origin(self,x,y,quiet=False):
		"""Set the display origin within the image"""
		if self.origin==(x,y) : return
		self.origin=(x,y)
		if not quiet : self.origin_update.emit((x,y))
		self.updateGL()

	def get_origin(self) : return self.origin



	def scroll_to(self,x=None,y=None,quiet=False):
		"""center the point on the screen"""
		if x==None:
			if y==None: return
			self.set_origin(self.origin[0],y*self.scale-old_div(self.height(),2),quiet)
		elif y==None:
			self.set_origin(x*self.scale-old_div(self.width(),2),self.origin[1],quiet)
		else: self.set_origin(x*self.scale-old_div(self.width(),2),y*self.scale-old_div(self.height(),2),quiet)

	def set_shapes(self,shapes):
		self.shapes = shapes
		self.shapechange=1

	def update_shapes(self,shapes):
		self.shapes.update(shapes)
		self.shapechange=1

	def register_scroll_motion(self,x,y,z):
		print("Self.list_data", self.list_data)
		if self.list_data!=None:

			self.ns_changed(z)
			#self.setup_shapes()
		animation = LineAnimation(self,self.origin,(x*self.scale-old_div(self.width(),2),y*self.scale-old_div(self.height(),2)))
		self.qt_parent.register_animatable(animation)
		return True


	#TODO: Scrolling change scale keep reset to center image.
	def set_scale(self,newscale,quiet=False):
		"""Adjusts the scale of the display. Tries to maintain the center of the image at the center"""
		if self.scale==newscale: return
		try:
			#print("Self.origin", self.origin)
			self.origin=(old_div(newscale,self.scale)*(old_div(self.width(),2.0)+self.origin[0])-old_div(self.width(),2.0),old_div(newscale,self.scale)*(old_div(self.height(),2.0)+self.origin[1])-old_div(self.height(),2.0))
			self.scale=newscale
			if not quiet : self.signal_set_scale.emit(newscale)
			self.inspector_update()
			self.force_display_update(set_clip=True)
			self.updateGL()
		except: pass

	def set_invert(self,val):
		if val: self.invert=1
		else : self.invert=0
		self.updateGL()

	def set_histogram(self,mode):
		self.histogram=mode
		self.updateGL()

	# def __set_display_image(self,val,alt,az):
	#
	# 	if self.full_data["nz"]==1:
	# 		print("Full data nz = 1")
	#
	# 		return False
	#
	#
	# 	# faster, so do it this way when not tilted
	# 	if alt==0 and az==0 :
	# 		self.set_clip()
	# 		print("Setting clip in the alt az 0")
	# 		self.data=self.full_data.get_clip(Region(0,0,val,self.full_data["nx"],self.full_data["ny"],1))
	# 		#print("Set data from set display")
	# 		#self.data=self.full_data
	# 		self.annotation=self.full_annotation.get_clip(Region(0,0,val,self.full_data["nx"],self.full_data["ny"],1))
	# 		#self.annotation=self.full_annotation
	# 		return True
	#
	# 	return False

	def force_display_update(self,set_clip=False):
		if set_clip:
			self.set_clip()
		#self.__set_display_image(zcent+self.zpos,self.az, self.alt)
		self.display_states = []

	def display_state_changed(self):
		display_states = []
		# FIXME - this should really be static
		display_states.append(self.width())
		display_states.append(self.height())
		display_states.append(self.origin[0])
		display_states.append(self.origin[1])
		display_states.append(self.scale)
		display_states.append(self.invert)
		display_states.append(self.minden)
		display_states.append(self.maxden)
		display_states.append(self.curmin)
		display_states.append(self.curmax)
		display_states.append(self.az)
		display_states.append(self.alt)
		if len(self.display_states) == 0:
			self.display_states = display_states
			return True
		else:
			for i in range(len(display_states)):

				if display_states[i] != self.display_states[i]:
					self.display_states = display_states
					return True

		return False

	def create_palette(self, n):
		#Create color palette for table items and indices
		l = [QtGui.QColor.fromHsv(0,0,255)]
		#l = []
		for i in range(1,n):
			h = ((130 + 139*i)%360)
			#s = (255+(1 - i/(n*4))*256)%256
			#v = (255+(1-i/(n*8))*256)%256
			#s = ((1-4*i/(n))*255)%256
			s = ((0.8-3*i/(n))*255)%256
			v = 255
			l.append(QtGui.QColor.fromHsv(h,s,v))
		return l

	def get_color_palette(self):
		return self.colors


	# def create_palette(self,n=256):
	# 	#Create color palette for table items and indices
	# 	l = []
	# 	#ctable = list(to_numpy(self.target.ctable))
	# 	ctable = self.target.ctable
	# 	print("Len ctable", len(ctable))
	# 	for i in range(0,256*3,3):
	# 		r=ctable[255*256*3+i]
	# 		g=ctable[255*256*3+i+1]
	# 		b=ctable[255*256*3+i+2]
	# 		l.append(QtGui.QColor.fromRgb(r,g,b))
	# 		# if i < 10:
	# 		# 	print("rgb2",r,g,b)
	# 	return l

	def create_RGB_list(self):
		if self.need_new_RGB == 0:
			return self.ctable

		else:
			ctable = []
			for v in range(256):
				# for i in (range(n)):
				# 	if i == 0:
				# 		h = 0
				# 		s = 0
				# 	else:
				#
				# 		h = ((130 + 139*i)%360)
				# 		#s = (255+(1 - i/(n*4))*256)%256
				# 		#v = (255+(1-i/(n*8))*256)%256
				# 		#s = ((1-4*i/(n))*255)%256
				# 		s = ((0.8-3*i/(n))*255)%256
				# 		#s = 255
				for fullv_color in self.colors:
					color = QtGui.QColor.fromHsv(fullv_color.hue(),fullv_color.saturation(),v)
					#color = QtGui.QColor.fromHsv(h,s,v)
					ctable.append(color.red())
					ctable.append(color.green())
					ctable.append(color.blue())
					# if i == 0 and v >=250:
					# 	print("I=0",color.red(),color.green(),color.blue())

			return (np.array(ctable)).tolist()




				#l.append(QtGui.QColor.fromHsv(h,s,v))


		#Create color palette for table items and indices


	def render_bitmap (self) :
		"""This will render the current bitmap into a string and
			return a tuple with 1 or 3 (grey vs RGB), width, height,
			and the raw data string, with optional histogram data appended
			(1024 bytes)."""

		if self.invert :
			pixden = (255, 0)
		else :
			pixden = (0, 255)


		# min_val = self.curmin
		# max_val = self.curmax
		# gam_val = self.gamma

		have_textures = (not self.glflags.npt_textures_unsupported ( ))

		value_size = 1  # 1 8-bit byte

		flags = 2  # histogram

		min_val = self.curmin
		max_val = self.curmax



		wid =  ((self.width() * value_size - 1) // 4) * 4+ 4
		wdt =  self.width()
		hgt =  self.height()

		x0  = 1 + int(old_div(self.origin[0], self.scale))
		y0  = 1 + int(old_div(self.origin[1], self.scale))

		#print("Zpos:",self.zpos)


		self.xform=Transform({"type":"eman","alt":self.alt,"az":self.az,"tx":self.full_data["nx"]//2,"ty":self.full_data["ny"]//2,"tz":self.full_data["nz"]//2+self.zpos})
		#print(self.xform)

		#print("Called setting clip")
		#self.set_xform(self.nx//2,self.ny//2,self.nz//2+self.zpos,self.alt,self.az)
		#self.xform=Transform({"type":"eman","alt":self.alt,"az":self.az,"tx":self.full_data["nx"]//2,"ty":self.full_data["ny"]//2,"tz":self.full_data["nz"]//2+self.zpos})

		#self.data=self.full_data.get_rotated_clip(Transform({"type":"eman","alt":self.alt,"az":self.az,"tx":self.full_data["nx"]//2,"ty":self.full_data["ny"]//2,"tz":self.full_data["nz"]//2+self.zpos}),(self.full_data["nx"],self.full_data["ny"],1))
		self.data=self.full_data.get_rotated_clip(self.xform,(self.full_data["nx"],self.full_data["ny"],1))

		#self.xform = Transform({"type":"eman","alt":self.alt,"az":self.az})
		#self.data=self.full_data.process("xform",{"transform":self.xform}).get_clip(Region(0,0,self.zpos,self.full_data["nx"],self.full_data["ny"],1))
		#self.annotation=self.full_annotation.get_rotated_clip(Transform({"type":"eman","alt":self.alt,"az":self.az,"tx":self.full_data["nx"]//2,"ty":self.full_data["ny"]//2,"tz":self.full_data["nz"]//2+self.zpos}),(self.full_data["nx"],self.full_data["ny"],1))
		self.annotation=self.full_annotation.get_rotated_clip(self.xform,(self.full_data["nx"],self.full_data["ny"],1),0)
		#self.annotation=self.full_annotation.process("xform",{"transform":self.xform}).get_clip(Region(0,0,self.zpos,self.full_data["nx"],self.full_data["ny"],1))
		#print("Done getting clip")
		#print(self.full_data,self.full_annotation,self.alt,self.az,self.zpos,x0,y0,wdt,hgt,wid)


		#print((to_numpy(self.data)[50,50]),(to_numpy(self.annotation)[50,50]))
		values  = self.data
		if len(self.disp_proc)>0 :
			#print(self.disp_proc)
			tmp=values.copy()
			for p in self.disp_proc: p.process_inplace(tmp)
			values=tmp

		annotation = self.annotation.copy()

		#RENDER BY GROUP : May not be necessary
		#if self.display_group:
		#if self.get_inspector():
		if self.display_group:

			temp = annotation.copy_head()
			temp.to_zero()

			it = QtWidgets.QTreeWidgetItemIterator(self.get_inspector().seg_tab.tree_set)

			while it.value():
				item = it.value()
				#if item.text(2) == "-1" or item in self.tree_sels:
				if item.text(2) == "-1" or item in self.get_inspector().seg_tab.get_selected_item():
					it +=1
					#pass
				else:
					index = int(item.text(0))
					val = int(item.text(2))
					temp += val*(annotation.process("threshold.binaryrange",{"high":index+0.1,"low":index-0.1}))
					annotation.process_inplace("threshold.rangetozero",{"maxval":(index+0.1),"minval":(index-0.1)})
					it +=1
			annotation += temp

		else:
			#print("No group to concern")
			pass

		# val = int(item.text(0))
		# print("Group_val",group_val)
		# temp += group_val*(ori.process("threshold.binaryrange",{"high":val+0.1,"low":val-0.1}))
		# self.target.get_full_annotation().process_inplace("threshold.rangetozero",{"maxval":(val+0.1),"minval":(val-0.1)})
		# it += 1

			#it +=1

			# for row in range(self.inspector.seg_tab.table_set.rowCount()):
			# 	index = int(self.inspector.seg_tab.table_set.item(row,0).text())
			# 	#print("index",index)
			# 	#if index <100:
			# 	print("Nodes",self.inspector.seg_tab.nodes)
			# 	self.inspector.seg_tab.nodes[index]
			# 	val = self.inspector.seg_tab.nodes[index].get_value()
			# 	temp += annotation.process("threshold.binaryrange",{"high":index+0.1,"low":index-0.1})*(val)
			# 	annotation.process_inplace("threshold.rangetozero",{"maxval":(index+0.1),"minval":(index-0.1)})
			# annotation += temp


		#
		# #it +=1
		# while it.value():
		# 	#print(item.text(0))
		# 	item = it.value()
		# 	#print("CP",item.text(0))
		# 	if new_parent.text(0) == "":
		# 		print("EMPTY")
		# 		item.setText(2,"-1")
		# 	else:
		# 		print("NOT EMPTY")
		# 		item.setText(2,new_parent.text(0))
		# 	it += 1


		# return_data = ( wid*3, hgt,
		# 					GLUtil.render_annotated24 (self.data, self.annotation,
		# 					x0, y0, wdt, hgt, wid*3,
		# 					self.scale, pixden[0], pixden[1],
		# 					min_val, max_val, flags))
		# return_data = ( wid*3, hgt,
		# 					GLUtil.render_annotated24(values, self.annotation,
		# 					x0, y0, wdt, hgt, wid*3,
		# 					self.scale, pixden[0], pixden[1],
		# 					min_val, max_val, flags))

		###CORRECT
		# return_data = ( wid*3, hgt,
		# 					GLUtil.render_annotated24(values, annotation,
		# 					x0, y0, wdt, hgt, wid*3,
		# 					self.scale, pixden[0], pixden[1],
		# 					min_val, max_val, flags))
		self.ctable = self.create_RGB_list()
		return_data = ( wid*3, hgt,
							GLUtil.render_annotated24(values, annotation,
							x0, y0, wdt, hgt, wid*3,
							self.scale, pixden[0], pixden[1],
							min_val, max_val, flags, self.ctable))
		# return_data = ( wid*3, hgt,
		# 					GLUtil.render_annotated24(values, annotation,
		# 					x0, y0, wdt, hgt, wid*3,
		# 					self.scale, pixden[0], pixden[1],
		# 					min_val, max_val, flags))
		return return_data


	# def recolor_by_group(self):
	#
	# 	values  = self.data
	# 	if len(self.disp_proc)>0 :
	# 		print(self.disp_proc)
	# 		tmp=values.copy()
	# 		for p in self.disp_proc: p.process_inplace(tmp)
	# 		values=tmp
	# 	#self.target.full_annotation=EMData('./segs/temp_temp.hdf')
	#
	# 	for row in range(self.inspector.seg_tab.table_set.rowCount()):
	# 		index = int(self.inspector.seg_tab.table_set.item(row,0).text())
	# 		if index <100:
	# 			val = self.nodes[index].get_value()
	# 			temp += self.get_annotation().process("threshold.binaryrange",{"high":index+0.1,"low":index-0.1})*(val)
	# 			self.annotation.process_inplace("threshold.rangetozero",{"maxval":(index+0.1),"minval":(index-0.1)})
	# 	self.annotation += temp




	def set_clip(self):
		#self.xform=Transform({"type":"eman","alt":self.alt,"az":self.az,"tx":self.full_data["nx"]//2,"ty":self.full_data["ny"]//2,"tz":self.full_data["nz"]//2+self.zpos})
		#print("Setting clip")
		if self.annotation:
			self.full_annotation.set_rotated_clip(self.xform,self.annotation)
		if self.data:
			self.full_data.set_rotated_clip(self.xform,self.data)


	def render_bitmap_old(self):   # no longer used - use new render_bitmap
		"""This will render the current bitmap into a string and return a tuple with 1 or 3 (grey vs RGB), width, height, and the raw data string"""
		if not self.invert : pixden=(0,255)
		else: pixden=(255,0)


		if self.histogram==1:
			a=(1,old_div((self.width()-1),4)*4+4,self.height(),GLUtil.render_amp8(self.data, 1+int(old_div(self.origin[0],self.scale)),1+int(old_div(self.origin[1],self.scale)),self.width(),self.height(),old_div((self.width()-1),4)*4+4,self.scale,pixden[0],pixden[1],self.curmin,self.curmax,self.gamma,34))
		elif self.histogram==2:
			a=(1,old_div((self.width()-1),4)*4+4,self.height(),GLUtil.render_amp8(self.data, 1+int(old_div(self.origin[0],self.scale)),1+int(old_div(self.origin[1],self.scale)),self.width(),self.height(),old_div((self.width()-1),4)*4+4,self.scale,pixden[0],pixden[1],self.curmin,self.curmax,self.gamma,98))
		else :
			a=(1,old_div((self.width()-1),4)*4+4,self.height(),GLUtil.render_amp8(self.data, 1+int(old_div(self.origin[0],self.scale)),1+int(old_div(self.origin[1],self.scale)),self.width(),self.height(),old_div((self.width()-1),4)*4+4,self.scale,pixden[0],pixden[1],self.curmin,self.curmax,self.gamma,2))
#			else :
#				a=(1,(self.width()-1)/4*4+4,self.height(),GLUtil.render_amp8(self.data, 1+int(self.origin[0]/self.scale),1+int(self.origin[1]/self.scale),self.width(),self.height(),(self.width()-1)/4*4+4,self.scale,pixden[0],pixden[1],self.curmin,self.curmax,self.gamma,6))

		return a

	def render(self):
		if self.full_data is None : return

		if not self.isVisible():
			return

		try:
			self.image_change_count = self.data["changecount"] # this is important when the user has more than one display instance of the same image, for instance in e2.py if
		except: pass # probably looking at an FFT image


		lighting = glIsEnabled(GL_LIGHTING)
		glDisable(GL_LIGHTING)

		if self.shapechange:
			self.setup_shapes()
			self.shapechange=0

		width = self.width()/2.0
		height = self.height()/2.0

		if not self.invert : pixden=(0,255)
		else: pixden=(255,0)


		update = False
		if self.display_state_changed():
			update = True

		if update:
			self.update_inspector_texture() # important for this to occur in term of the e2desktop only

#		print "render",update,self.image_change_count


		render = False
		if self.use_display_list:
			if update:
				if self.main_display_list != -1:
					glDeleteLists(self.main_display_list,1)
					self.main_display_list = -1

			if self.main_display_list == -1:
				self.main_display_list = glGenLists(1)
				render = True
		else: render = True

		if render:
			(w,h,a)=self.render_bitmap()
			gl_render_type = GL_RGB

			#TODO - not fixed below for RGB
			if not self.glflags.npt_textures_unsupported():

				self.hist=struct.unpack('256i',a[-1024:])

				if self.tex_name != 0: glDeleteTextures(self.tex_name)
				self.tex_name = glGenTextures(1)
				if ( self.tex_name <= 0 ):
					raise("failed to generate texture name")

				GL.glBindTexture(GL.GL_TEXTURE_2D,self.tex_name)
				glPixelStorei(GL_UNPACK_ALIGNMENT,4)
				GL.glTexImage2D(GL.GL_TEXTURE_2D,0,gl_render_type,old_div(w,3),h,0,gl_render_type, GL.GL_UNSIGNED_BYTE, a)

				glNewList(self.main_display_list,GL_COMPILE)
				GL.glBindTexture(GL.GL_TEXTURE_2D,self.tex_name)
				glPushMatrix()
				glTranslatef(width,height,0)
				self.__draw_texture(self.tex_name,-width,-height,width,height)
				glPopMatrix()

			else:
				self.hist=struct.unpack('256i',a[-1024:])
				glNewList(self.main_display_list,GL_COMPILE)
				#GL.glRasterPos(0,self.height()-1)
				#GL.glPixelZoom(1.0,-1.0)
				GL.glRasterPos(0,0)
				GL.glDrawPixels(self.width(),self.height(),gl_render_type,GL.GL_UNSIGNED_BYTE,a)
		else:
			glCallList(self.main_display_list)

		if self.use_display_list and render :
			glEndList()
			glCallList(self.main_display_list)


		if self.frozen or self.isexcluded:
			if self.isexcluded:	glColor(0,0.1,0,1)
			else:glColor(0,0,0.1,1)

			glPushMatrix()
			glTranslatef(width,height,0)
			self.__draw_square_shaded_region(-width,-height,width,height)
			glPopMatrix()

		if self.inspector:
			if self.invert: self.inspector.set_hist(self.hist,self.maxden,self.minden)
			else: self.inspector.set_hist(self.hist,self.minden,self.maxden)

		if self.shapelist != 0 and self.display_shapes:
			emshape.EMShape.font_renderer=self.font_renderer		# Important !  Each window has to have its own font_renderer. Only one context active at a time, so this is ok.
			GL.glColor(.7,1.0,.7,1.0)
			GL.glPushMatrix()
			#TODO: change self.render_bitmap()'s C++ functions so this ugly hack isn't necessary.
			GL.glTranslate(-self.scale*(int(old_div(self.origin[0],self.scale))+0.5),-self.scale*(int(old_div(self.origin[1],self.scale))+0.5),0.1)
			GL.glScalef(self.scale,self.scale,1.0)
			GL.glCallList(self.shapelist)
			GL.glPopMatrix()

		self.__draw_hud()

		if ( lighting ): glEnable(GL_LIGHTING)

	def __draw_texture(self,texture_handle,xmin,ymin,xmax,ymax):

		texture_2d_was_enabled = GL.glIsEnabled(GL_TEXTURE_2D)
		if not texture_2d_was_enabled:glEnable(GL_TEXTURE_2D)

		glBindTexture(GL_TEXTURE_2D, texture_handle)
		#glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP)
		#glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP)
		# using GL_NEAREST ensures pixel granularity
		# using GL_LINEAR blurs textures and makes them more difficult
		# to interpret (in cryo-em)
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)
		# this makes it so that the texture is impervious to lighting
		#glTexEnvfv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_COLOR, [0.2,0.4,0.8,0.5])
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE)

		# POSITIONING POLICY - the texture occupies the entire screen area
		glBegin(GL_QUADS)

		glTexCoord2f(0,0)
		glVertex2f(xmin,ymax)

		glTexCoord2f(1,0)
		glVertex2f(xmax,ymax)

		glTexCoord2f(1,1)
		glVertex2f(xmax,ymin)

		glTexCoord2f(0,1)
		glVertex2f(xmin,ymin)

		glEnd()

		if not texture_2d_was_enabled:glDisable(GL_TEXTURE_2D)


	def __draw_square_shaded_region(self,xmin,ymin,xmax,ymax):

		GL.glEnable(GL.GL_BLEND);
		depth_testing_was_on = GL.glIsEnabled(GL.GL_DEPTH_TEST);
		GL.glDisable(GL.GL_DEPTH_TEST)
		GL.glBlendEquation(GL.GL_FUNC_ADD)
		#GL.glBlendFunc(GL.GL_SRC_ALPHA,GL.GL_ONE_MINUS_SRC_ALPHA);
		GL.glBlendFunc(GL.GL_ONE,GL.GL_ONE)
		glBegin(GL_QUADS)

		glVertex2f(xmin,ymax)

		glVertex2f(xmax,ymax)

		glVertex2f(xmax,ymin)

		glVertex2f(xmin,ymin)

		glEnd()

		glDisable(GL_TEXTURE_2D)

		glDisable( GL.GL_BLEND)
		if (depth_testing_was_on):	GL.glEnable(GL.GL_DEPTH_TEST)

	def resize_event(self,width,height):
		if self.init_size :
			if self.origin == (0,0):
				#self.origin = ((self.data.get_xsize() - self.width())/2.0, (self.data.get_ysize() - self.height())/2.0 )
				self.oldsize=(width,height)
			self.init_size = False
		else:
			if self.origin == (0,0):
				#self.origin=((self.oldsize[0]/2.0+self.origin[0])-self.width()/2.0,(self.oldsize[1]/2.0+self.origin[1])-self.height()/2.0)
				self.oldsize=(width,height)

		display_width = self.width()
		display_height = self.height()

		data = self.data

		#print data
		if data != None:
			try:
				pixel_x = self.scale*data.get_xsize()
				pixel_y = self.scale*data.get_ysize()
			except:
				try:
					pixel_x = self.scale*data[0].get_xsize()
					pixel_y = self.scale*data[0].get_ysize()
				except:
					pixel_x = 256
					pixel_y = 256
		else: return


		ox = self.origin[0]
		oy = self.origin[1]

		# if the image is off the screen automatically center it
		if pixel_x + ox < 0: ox = old_div(- (display_width-pixel_x),2.0)
		if pixel_y + oy < 0: oy = old_div(- (display_height-pixel_y),2.0)


		# this operation keeps the image iff it is completely visible
		if pixel_x < display_width:	ox = old_div(- (display_width-pixel_x),2.0)
		if pixel_y < display_width: oy =  old_div(- (display_height-pixel_y),2.0)

		self.origin = (ox,oy)

	def init_circle_list(self):
		if self.circle_dl == None:
			self.circle_dl = glGenLists(1)
			glNewList(self.circle_dl,GL_COMPILE)
			glBegin(GL_LINE_LOOP)
			d2r=old_div(pi,180.0)
			for i in range(90): glVertex(sin(i*d2r*4.0),cos(i*d2r*4.0))
			glEnd()
			glEndList()

	def setup_shapes(self):
		if self.shapelist != 0: GL.glDeleteLists(self.shapelist,1)
		self.init_circle_list()
		self.shapelist = glGenLists(1)

		#context = OpenGL.contextdata.getContext(None)
		#print "Image2D context is", context,"display list is",self.shapelist

		# make our own circle rather than use gluDisk or somesuch
		emshape.EMShape.font_renderer=self.font_renderer		# Important !  Each window has to have its own font_renderer. Only one context active at a time, so this is ok.
		glNewList(self.shapelist,GL_COMPILE)

		isanimated = False
		alpha = 1.0
		if len(self.shapes) > 0:

			for k in list(self.shapes.keys()):
				shape = self.shapes[k]
				if not isinstance(shape,EMShape) : continue
				glLineWidth(2*self.devicePixelRatio())
				if shape.isanimated:
					isanimated = True
					alpha = shape.blend

					break

		if isanimated:
			GL.glEnable(GL.GL_BLEND);
			depth_testing_was_on = GL.glIsEnabled(GL.GL_DEPTH_TEST);
			GL.glDisable(GL.GL_DEPTH_TEST);
			try:GL.glBlendEquation(GL.GL_FUNC_ADD)
			except:pass
			GL.glBlendFunc(GL.GL_SRC_ALPHA,GL.GL_ONE_MINUS_SRC_ALPHA);

		glPointSize(2*self.devicePixelRatio())
		for k,s in list(self.shapes.items()):
			### handle boxes for 3D images
			if s.shape[0] == "ellipse":
				mxlen=11
			else:
				mxlen=10

			if self.list_data!=None and len(s.shape)==mxlen:
				z_idx=s[mxlen-1]
				if z_idx!=self.zpos:
					continue

			if k == self.active[0]:
				if not isinstance(s,EMShape) :
					print("Invalid shape in EMImage : ",s)
					continue
				vals = s.shape[1:8]
				vals[0:3] = self.active[1:4]
				if s.shape[0] == "rect":
					GLUtil.colored_rectangle(vals,alpha,True)
				elif s.shape[0] == "rectpoint":
					GLUtil.colored_rectangle(vals,alpha,True)
				continue
			try:
				if s.shape[0] == "rectpoint":
					GLUtil.colored_rectangle(s.shape[1:8],alpha,True)
				elif  s.shape[0] == "rcircle":
					glPushMatrix()
					glColor(*s.shape[1:4])
					x1 = s.shape[4]
					x2 = s.shape[6]
					y1 = s.shape[5]
					y2 = s.shape[7]
					glTranslate(old_div((x1+x2),2.0), old_div((y1+y2),2.0),0)
					glScalef(old_div((x1-x2),2.0), old_div((y1-y2),2.0),1.0)
					glCallList(self.circle_dl)
					glPopMatrix()
				elif  s.shape[0] == "rcirclepoint":
					glColor(*s.shape[1:4])
					glPushMatrix()
					x1 = s.shape[4]
					x2 = s.shape[6]
					y1 = s.shape[5]
					y2 = s.shape[7]
					glTranslate(old_div((x1+x2),2.0), old_div((y1+y2),2.0),0)
					glScalef(old_div((x1-x2),2.0), old_div((y1-y2),2.0),1.0)
					glCallList(self.circle_dl)
					glPopMatrix()
					glBegin(GL_POINTS)
					glVertex( old_div((x1+x2),2.0), old_div((y1+y2),2.0),0);
					glEnd()
				elif s.shape[0] == "circle":
					GL.glPushMatrix()
					p = s.shape
					GL.glColor(*p[1:4])
					v= (p[4],p[5])
					GL.glLineWidth(p[7]*self.devicePixelRatio())
					GL.glTranslate(v[0],v[1],0)
					GL.glScalef(p[6],p[6],1.0)
					glCallList(self.circle_dl)
					GL.glPopMatrix()
				elif s.shape[0] == "ellipse":
					GL.glPushMatrix()
					p = s.shape
					GL.glColor(*p[1:4])
					v= (p[4],p[5])
					v2=(p[4]+1,p[5]+1)
					sc=v2[0]-v[0]
					GL.glLineWidth(p[9]*self.devicePixelRatio())
					GL.glTranslate(v[0],v[1],0)
					GL.glScalef((v2[0]-v[0]),(v2[1]-v[1]),1.0)
					GL.glRotatef(p[8],0,0,1.0)
					GL.glScalef(p[6],p[7],1.0)
					glCallList(self.circle_dl)
					GL.glPopMatrix()
				elif s.shape[0][:3]!="scr":
#					print "shape",s.shape
					GL.glPushMatrix()		# The push/pop here breaks the 'scr*' shapes !

					s.draw()		# No need for coordinate transform any more
					GL.glPopMatrix()
#					GLUtil.colored_rectangle(s.shape[1:8],alpha)
			except: pass


		# We do the scr* shapes last since they mess up the matrix
		for k,s in list(self.shapes.items()):
			try:
				if s.shape[0][:3]=="scr":
#					print "shape",s.shape
#					GL.glPushMatrix()		# The push/pop here breaks the 'scr*' shapes !
					s.draw()		# No need for coordinate transform any more
#					GL.glPopMatrix()
#					GLUtil.colored_rectangle(s.shape[1:8],alpha)
			except: pass


		if isanimated:
			GL.glDisable( GL.GL_BLEND);
			if (depth_testing_was_on):
				GL.glEnable(GL.GL_DEPTH_TEST)

		if self.eraser_shape != None:

			GL.glPushMatrix()
			p=self.eraser_shape.shape
			GL.glColor(*p[1:4])
			v= (p[4],p[5])
			GL.glLineWidth(p[7]*self.devicePixelRatio())
			GL.glTranslate(v[0],v[1],0)
			GL.glScalef(p[6],p[6],1.0)
			glCallList(self.circle_dl)
			GL.glPopMatrix()
			#self.eraser_shape.draw()

		glEndList()

	def inspector_update(self):
		if not self.inspector is None and not self.data is None:
			self.inspector.set_limits(self.minden,self.maxden,self.curmin,self.curmax)
			self.inspector.set_scale(self.scale)
			self.inspector.update_brightness_contrast()
			self.inspector.mtapix.setValue(self.data["apix_x"])
			self.inspector.update_zrange()

	def get_inspector(self):
		if not self.inspector:
			self.inspector=EMAnnotateInspector2D(target=self)
			self.inspector_update()
		return self.inspector


	def set_active(self,n,r,g,b):
		self.active=(n,r,g,b)
		self.shapechange=1
		#self.updateGL()

	def update_animation(self):

		if not self.isanimated:
			return False

		self.time += self.timeinc
		if self.time > 1:
			self.time = 1
			self.isanimated = False
			self.set_origin(self.endorigin[0],self.endorigin[1])
			return True

		# get time squared
		tinv = 1-self.time
		t1 = tinv**2
		t2 = 1-t1

		x = t1*self.startorigin[0] + t2*self.endorigin[0]
		y = t1*self.startorigin[1] + t2*self.endorigin[1]
		self.set_origin(x,y)
		return True


	def animation_done_event(self,animation):
		if animation == self.key_mvt_animation:
			self.key_mvt_animation = None

		if isinstance(animation,SingleValueIncrementAnimation):
			self.set_animation_increment(animation.get_end())
		elif isinstance(animation,LineAnimation):
			self.set_line_animation(*animation.get_end())

	def set_animation_increment(self,increment):
		for shape in list(self.shapes.items()):
			shape[1].set_blend(increment)

		self.shapechange = True

	def set_line_animation(self,x,y):
		self.origin=(x,y)
		self.display_states = [] #forces a display list update

	def update_blend(self):
		ret = False
		for shape in list(self.shapes.items()):
			s = shape[1]
			if s.isanimated:
				v = s.incblend()
				ret = True

		return ret

	def add_vector_overlay(self,x0,y0,x1,y1,v=None):
		"""Add an array of vectors to the image as a shape overlay. This can be used to make pixel-centered vector plots
		on top of images. This will be represented as additions to the shape list for the image. Only a single vector
		overlay may be added to an image at a time. If provided, v will map to color (red low, blue high)"""

		if len(x0)!=len(y0) or len(x0)!=len(x1) or len(x0)!=len(y1):
			raise Exception("Error, vector overlays must have 4 equal-length arrays")

		if v is None: newshapes={f"vl({i})":EMShape(["line",0.2,0.8,0.2,x0[i],y0[i],x1[i],y1[i],1.0]) for i in range(len(x0))}
		else:
			vmin=min(v)
			vmax=max(v)
			v=(np.array(v)-vmin)/(vmax-vmin)
			newshapes={f"vl({i})":EMShape(["line",1.0-v[i],0.2,v[i],x0[i],y0[i],x1[i],y1[i],1.5]) for i in range(len(x0))}
			newshapes2={f"vp({i})":EMShape(["point",1.0-v[i],0.2,v[i],x0[i],y0[i],1.5]) for i in range(len(x0))}
		self.del_shapes()
		self.add_shapes(newshapes)
		self.add_shapes(newshapes2)
		self.updateGL()

	def add_shape(self,k,s):
		"""Add an EMShape object to be overlaid on the image. Each shape is
		keyed into a dictionary, so different types of shapes for different
		purposes may be simultaneously displayed.

		"""
		self.shapes[k]=s
		self.shapechange=1

	def add_eraser_shape(self,k,args):
		"""Add an EMShape object to be overlaid on the image. Each shape is
		keyed into a dictionary, so different types of shapes for different
		purposes may be simultaneously displayed.

		"""
		self.makeCurrent() # this is important  when you have more than one OpenGL context operating at the same time

		if k != "None":
			#self.eraser_shape=EMShape(args)
			self.shapes["eraser"] = EMShape(args)
		else:
			self.eraser_shape = None
			try:
				self.shapes.pop("eraser")
			except: pass

		self.shapechange=1
		#self.updateGL()

	def add_shapes(self,d,register_animation=False):
		if register_animation:
			animation = SingleValueIncrementAnimation(self,0,1)

			self.qt_parent.register_animatable(animation)
		self.shapes.update(d)
		self.shapechange=1
		#self.updateGL()


	def del_shape(self,p):
		try:
			self.shapes.pop(p)
			self.shapechange=1
		except:pass

	def del_shapes(self,k=None):
		if k:
			try:
				for i in k:
					del self.shapes[i]
			except: del self.shapes[k]
		else:
			self.shapes={}
		self.shapechange=1
		#self.updateGL()

	def scr_to_img(self,v0,v1=None):
		#TODO: origin_x and origin_y are part of the hack in self.render() and self.render_bitmap()
		origin_x = self.scale*(int(self.origin[0]/self.scale)+0.5)
		origin_y = self.scale*(int(self.origin[1]/self.scale)+0.5)

		try: img_coords = (((v0+origin_x)/self.scale), ((self.height()-(v1-origin_y))/self.scale) )
		except:	img_coords = (((v0[0]+origin_x)/self.scale),((self.height()-(v0[1]-origin_y))/self.scale))

#		print "Screen:", v0, v1
#		print "Img:", img_coords
#		print "Screen:", self.img_to_scr(img_coords)

		return img_coords

	def img_to_scr(self,v0,v1=None):
		#TODO: origin_x and origin_y are part of the hack in self.render() and self.render_bitmap()
		origin_x = self.scale*(int(old_div(self.origin[0],self.scale))+0.5)
		origin_y = self.scale*(int(old_div(self.origin[1],self.scale))+0.5)

		try: return (v0*self.scale-origin_x,self.height()-v1*self.scale+origin_y)
		except: return (v0[0]*self.scale-origin_x,self.height()-v0[1]*self.scale+origin_y)


	def closeEvent(self,event) :

		EMGLWidget.closeEvent(self,event)
		try:
			os.rmdir('./segs/temp')
		except:
			pass

		try:
			for w in self.inspector.pspecwins: w.close()
		except: pass

	def dragEnterEvent(self,event):
#		f=event.mimeData().formats()
#		for i in f:
#			print str(i)

		if event.mimeData().hasFormat("application/x-eman"):
			event.setDropAction(Qt.CopyAction)
			event.accept()

	def dropEvent(self,event):
		if EMAN2.GUIbeingdragged:
			self.set_data(EMAN2.GUIbeingdragged)
			EMAN2.GUIbeingdragged=None
		if event.mimeData().hasFormat("application/x-eman"):
			x=loads(event.mimeData().data("application/x-eman"))
			self.set_data(x)
			event.acceptProposedAction()

	def do_probe(self,x,y):
		"response to a probe mouse click/drag"
		try: sz=int(self.inspector.ptareasize.getValue())
		except: sz=16
		x,y=int(x),int(y)

		self.del_shape("PROBE")
		self.add_shape("PROBE",EMShape(("rectpoint",.5,.5,.1,x-old_div(sz,2),y-old_div(sz,2),x+old_div((sz+1),2),y+old_div((sz+1),2),2)))
		self.updateGL()

		clp=self.get_data().get_clip(Region(x-old_div(sz,2),y-old_div(sz,2),sz,sz))
		self.inspector.ptpointval.setText("Point Value: %1.3f"%(self.get_data())[x,y])
		self.inspector.ptareaavg.setText("Area Avg: %1.3f"%clp["mean"])
		self.inspector.ptareaavgnz.setText("Area Avg (!=0): %1.3f"%clp["mean_nonzero"])
		self.inspector.ptareasig.setText("Area Sig: %1.3f"%clp["sigma"])
		self.inspector.ptareasignz.setText("Area Sig (!=0): %1.3f"%clp["sigma_nonzero"])
		self.inspector.ptareaskew.setText("Area Skewness: %1.3f"%clp["skewness"])
		self.inspector.ptareakurt.setText("Area Kurtosis: %1.3f"%clp["kurtosis"])
		self.inspector.ptcoord.setText("Center Coord: %d, %d"%(x,y))
		self.inspector.ptcoord2.setText("dcen (%d, %d)"%(x-self.data["nx"]//2,y-self.data["ny"]//2))

	def mousePressEvent(self, event):
		lc=self.scr_to_img(event.x(),event.y())
		if event.button()==Qt.MidButton or (event.button()==Qt.LeftButton and event.modifiers()&Qt.AltModifier):
			self.show_inspector(1)
		elif event.button()==Qt.RightButton or (event.button()==Qt.LeftButton and event.modifiers()&Qt.AltModifier):
			try:
				get_application().setOverrideCursor(Qt.ClosedHandCursor)
			except: # if we're using a version of qt older than 4.2 than we have to use this...
				get_application().setOverrideCursor(Qt.SizeAllCursor)
			self.rmousedrag=(event.x(),event.y())
			if event.buttons()&Qt.RightButton:
				self.mousedrag.emit(event, lc)
				#print("Rmousedrag:",self.rmousedrag)
		else:
			if self.mouse_mode_dict[self.mouse_mode] == "emit":
				lc=self.scr_to_img(event.x(),event.y())
				self.mousedown.emit(event, lc)
			elif self.mouse_mode_dict[self.mouse_mode] == "probe":
				if event.buttons()&Qt.LeftButton:
					lc=self.scr_to_img(event.x(),event.y())
					self.do_probe(lc[0],lc[1])
			elif self.mouse_mode_dict[self.mouse_mode] == "measure":
				if event.buttons()&Qt.LeftButton:
					lc=self.scr_to_img(event.x(),event.y())
		#			self.del_shape("MEASL")
					self.del_shape("MEAS")
					self.add_shape("MEAS",EMShape(("line",.5,.1,.5,lc[0],lc[1],lc[0]+1,lc[1],2)))
					self.updateGL()
###TO READ
			elif self.mouse_mode_dict[self.mouse_mode] == "draw":
				if event.buttons()&Qt.LeftButton:
					inspector = self.get_inspector()
					lc=self.scr_to_img(event.x(),event.y())
					if inspector:
						self.drawr1=int(float(inspector.dtpen.text()))
						self.drawv1=float(inspector.dtpenv.text())
						self.drawr2=int(float(inspector.dtpen2.text()))
						self.drawv2=float(inspector.dtpenv2.text())
						self.get_data().process_inplace("mask.paint",{"x":lc[0],"y":lc[1],"z":0,"r1":self.drawr1,"v1":self.drawv1,"r2":self.drawr2,"v2":self.drawv2})
						self.force_display_update()
						self.updateGL()
			elif self.mouse_mode_dict[self.mouse_mode]=="seg":
				if event.buttons()&Qt.LeftButton:
					inspector = self.get_inspector()
					#inspector.show()
					lc=self.scr_to_img(event.x(),event.y())
					#current_class = self.current_class
					if inspector:
						current_class = inspector.seg_tab.get_current_class()
						pen_width = inspector.seg_tab.get_pen_width()


							#self.get_annotation().process_inplace("mask.paint",{"x":lc[0],"y":lc[1],"z":0,"r1":pen_width,"v1":current_class,"r2":pen_width,"v2":0})
						self.get_annotation().process_inplace("mask.paint",{"x":lc[0],"y":lc[1],"z":0,"r1":pen_width,"v1":current_class,"r2":pen_width,"v2":0})
						#print("Turn point", int(lc[0]),int(lc[1]),"to 2")
						self.force_display_update(set_clip=True)
						self.updateGL()


	def mouseMoveEvent(self, event):
		lc=self.scr_to_img(event.x(),event.y())

		if self.rmousedrag:
			self.set_origin(self.origin[0]+self.rmousedrag[0]-event.x(),self.origin[1]-self.rmousedrag[1]+event.y())
			self.rmousedrag=(event.x(),event.y())
			if event.buttons()&Qt.RightButton:
				self.mousedrag.emit(event, lc)
				#print("Rmousedrag:",self.rmousedrag)
			#self.emit(QtCore.SIGNAL("origin_update"),self.origin)
			try: self.updateGL()
			except: pass
		else:
			if self.mouse_mode_dict[self.mouse_mode] == "emit":
				lc=self.scr_to_img(event.x(),event.y())
				if event.buttons()&Qt.LeftButton:
					self.mousedrag.emit(event, lc)
				else:
					self.mousemove.emit(event, lc)
			elif self.mouse_mode_dict[self.mouse_mode] == "probe":
				if event.buttons()&Qt.LeftButton:
					lc=self.scr_to_img(event.x(),event.y())
					self.do_probe(lc[0],lc[1])
			elif self.mouse_mode_dict[self.mouse_mode] == "measure":
				if event.buttons()&Qt.LeftButton:
					lc=self.scr_to_img(event.x(),event.y())
					current_shapes = self.get_shapes()
					self.add_shape("MEAS",EMShape(("line",.5,.1,.5,current_shapes["MEAS"].shape[4],current_shapes["MEAS"].shape[5],lc[0],lc[1],2)))

					dx=lc[0]-current_shapes["MEAS"].shape[4]
					dy=lc[1]-current_shapes["MEAS"].shape[5]
					#self.add_shape("MEASL",EMShape(("label",.1,.1,.1,lc[0]+2,lc[1]+2,"%d,%d - %d,%d"%(current_shapes["MEAS"].shape[4],current_shapes["MEAS"].shape[5],lc[0],lc[1]),9,-1)))

					inspector = self.get_inspector()
					if inspector:
						apix=inspector.mtapix.value

						inspector.mtshoworigin.setText("Start: %d , %d"%(current_shapes["MEAS"].shape[4],current_shapes["MEAS"].shape[5]))
						inspector.mtshowend.setText("  End: %d , %d"%(lc[0],lc[1]))
						inspector.mtshowlen.setText("dx,dy: %1.2f A, %1.2f A"%(dx*apix,dy*apix))
						inspector.mtshowlen2.setText("Len: %1.3f A"%(hypot(dx,dy)*apix))

						try: inspector.mtshowval.setText("Value: %1.4g"%inspector.target().data[int(lc[0]),int(lc[1])])
						except:
							idx=inspector.target().zpos
							inspector.mtshowval.setText("Value: %1.4g"%inspector.target().list_data[idx][int(lc[0]),int(lc[1])])
						inspector.mtshowval2.setText("  ")

						#except: pass

					self.update_inspector_texture()
					self.updateGL()
			elif self.mouse_mode_dict[self.mouse_mode] == "draw":
				if event.buttons()&Qt.LeftButton:
					lc=self.scr_to_img(event.x(),event.y())
					self.get_data().process_inplace("mask.paint",{"x":lc[0],"y":lc[1],"z":0,"r1":self.drawr1,"v1":self.drawv1,"r2":self.drawr2,"v2":self.drawv2})
					self.force_display_update()
					self.updateGL()

			elif self.mouse_mode_dict[self.mouse_mode] =="seg":

				if event.buttons()&Qt.LeftButton:
					inspector = self.get_inspector()
					#inspector.show()
					#lc=self.scr_to_img(event.x(),event.y())
					#current_class = self.current_class
					#if inspector:
					current_class = inspector.seg_tab.get_current_class()
					pen_width = inspector.seg_tab.get_pen_width()
						#print("Seg:", lc)
						#NEED TO DO SOMETHING HERE???
						#print(self.get_annotation()[int(lc[0]),int(lc[1])])
						#print(self.get_annotation())
						#a = to_numpy(self.get_annotation())
						#for i in range(int(lc[0])-20,int(lc[0])+20):
							#for j in range(int(lc[1])-20,int(lc[1])+20):
						#self.get_annotation().process_inplace("mask.sharp",{"dx":lc[0],"dy":lc[1],"dz":0,"inner_radius":0,"outer_radius":10,"value":2})
					self.get_annotation().process_inplace("mask.paint",{"x":lc[0],"y":lc[1],"z":0,"r1":pen_width,"v1":current_class,"r2":pen_width,"v2":0})
						#print("Turn point", int(lc[0]),int(lc[1]),"to 2")
					self.force_display_update(set_clip=1)
					self.updateGL()


	def mouseReleaseEvent(self, event):
		get_application().setOverrideCursor(Qt.ArrowCursor)
		lc=self.scr_to_img(event.x(),event.y())
		current_shapes = self.get_shapes()
		if self.rmousedrag:
			self.rmousedrag=None
		else:
			if self.mouse_mode_dict[self.mouse_mode] == "emit":
				lc=self.scr_to_img(event.x(),event.y())
				self.mouseup.emit(event, lc)
			elif self.mouse_mode_dict[self.mouse_mode] == "measure":
				if event.buttons()&Qt.LeftButton:
					self.add_shape("MEAS",EMShape(("line",.5,.1,.5,current_shapes["MEAS"].shape[4],current_shapes["MEAS"].shape[5],lc[0],lc[1],2)))
			elif self.mouse_mode_dict[self.mouse_mode] == "draw":
				if event.button()==Qt.LeftButton:
					self.force_display_update()
					self.updateGL()
			elif self.mouse_mode_dict[self.mouse_mode] == "seg":
				if event.button()==Qt.LeftButton:
					self.force_display_update(set_clip=1)
					self.updateGL()

	def wheelEvent(self, event):
		if self.mouse_mode==0 and event.modifiers()&Qt.ShiftModifier:
			self.mousewheel.emit(event)
			return
		#print(event.angleDelta().y())
		if event.angleDelta().y() > 0:
			self.set_scale(self.scale * self.mag )
			#print("Scale", self.scale,"Mag", self.mag)
			self.mousewheel.emit(event)
		elif event.angleDelta().y() < 0:
			self.set_scale(self.scale * self.invmag )
			#print("Scale", self.scale,"Mag", self.mag)
			self.mousewheel.emit(event)
		# The self.scale variable is updated now, so just update with that
		if self.inspector: self.inspector.set_scale(self.scale)


	def get_z_cent(self):
		return (self.get_full_data()["nz"])//2

	def mouseDoubleClickEvent(self,event):
		return
#		if platform.system() == "Darwin":
#			self.wheel_navigate = not self.wheel_navigate
#		else:
#			print "double click only performs a function on Mac"

	def __key_mvt_animation(self,dx,dy):
		if self.key_mvt_animation == None:
			new_origin=(self.origin[0]+dx,self.origin[1]+dy)
			self.key_mvt_animation = LineAnimation(self,self.origin,new_origin)
			self.qt_parent.register_animatable(self.key_mvt_animation)
		else:
			new_origin = self.key_mvt_animation.get_end()
			new_origin = (new_origin[0]+dx,new_origin[1]+dy)
			self.key_mvt_animation.set_end(new_origin)

	def keyPressEvent(self,event):
		if event.key() == Qt.Key_F1:
			self.display_web_help("http://blake.bcm.edu/emanwiki/EMAN2/Programs/emimage2d")

		elif event.key() == Qt.Key_Up:
			if self.data != None:
				if (self.zpos + self.get_z_cent()) < self.get_full_data()["nz"]-1:
					self.zpos += 1
					if self.inspector:
						self.inspector.ns.setValue(self.zpos)
				self.force_display_update(set_clip=True)
				self.updateGL()
				self.keypress.emit(event)

#			else:
#				self.__key_mvt_animation(0,self.height()*.1)
			#self.signal_increment_list_data.emit(1)

		elif event.key() == Qt.Key_Down:
			#if self.list_data != None:
				#self.increment_list_data(-1)
			if self.data != None:
				if (self.zpos + self.get_z_cent()) > 0:
					self.zpos -= 1
					if self.inspector:
						self.inspector.ns.setValue(self.zpos)
				self.force_display_update(set_clip=True)
				self.updateGL()
				self.keypress.emit(event)

#			else:
#				self.__key_mvt_animation(0,-self.height()*.1)
			#self.signal_increment_list_data.emit(-1)

		elif event.key() == Qt.Key_Right:
			self.__key_mvt_animation(self.width()*.1,0)
		elif event.key() == Qt.Key_Left:
			self.__key_mvt_animation(-self.width()*.1,0)
		elif event.key()==Qt.Key_W :
			self.__key_mvt_animation(0,self.height())
		elif event.key()==Qt.Key_S :
			self.__key_mvt_animation(0,-self.height())
		elif event.key()==Qt.Key_D  :
			self.__key_mvt_animation(self.width(),0)
		elif event.key()==Qt.Key_A :
			self.__key_mvt_animation(-self.width(),0)
		elif event.key()==Qt.Key_Space:
			self.display_shapes = not self.display_shapes
			self.updateGL()
		elif event.key()==Qt.Key_P :
			try: self.get_inspector().do_pspec_single(0)
			except: pass
		elif event.key()==Qt.Key_C:
			self.auto_contrast()
		elif event.key()==Qt.Key_I:
			self.show_inspector(6)
		else:
			self.keypress.emit(event)



	def increment_list_data(self,delta):
		'''
		delta is negative or positive, indicating forward and backwards movement
		'''
		if self.list_data != None:
			if delta > 0:
				if (self.zpos < (len(self.list_data)-1)):
					self.zpos += 1
					self.get_inspector().set_image_idx(self.zpos+1,1)
					self.setup_shapes()
					self.force_display_update()

			elif delta < 0:
				if (self.zpos > 0):
					self.zpos -= 1
					self.get_inspector().set_image_idx(self.zpos-1,1)
					self.setup_shapes()
					self.force_display_update()

	def leaveEvent(self,event):
		get_application().setOverrideCursor(Qt.ArrowCursor)
		if self.rmousedrag:
			self.rmousedrag=None

	def __draw_hud(self):
		if self.full_data is None : return

		width = self.width()
		height = self.height()
		glMatrixMode(GL_PROJECTION)
		glPushMatrix()
		glLoadIdentity()
		glOrtho(0,width,0,height,-100,100)
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		glEnable(GL_LIGHTING)
		glEnable(GL_NORMALIZE)
		glMaterial(GL_FRONT,GL_AMBIENT,(0.2, 1.0, 0.2,1.0))
		glMaterial(GL_FRONT,GL_DIFFUSE,(0.2, 1.0, 0.9,1.0))
		glMaterial(GL_FRONT,GL_SPECULAR,(1.0, 0.5, 0.2,1.0))
		glMaterial(GL_FRONT,GL_SHININESS,20.0)
		enable_depth = glIsEnabled(GL_DEPTH_TEST)
		glDisable(GL_DEPTH_TEST)
		glColor(1.0,1.0,1.0)
#		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE)
		glNormal(0,0,1)
		glEnable(GL_TEXTURE_2D)
		n = self.full_data["nz"]
		pos = (n)//2 +self.zpos
		#string = f"{self.zpos} ({n})"

		string = f"Z: {round(pos,2)}"
		bbox = self.font_renderer.bounding_box(string)
		x_offset = width-(bbox[3]-bbox[0]) - 10
		y_offset = 10


		glPushMatrix()
		glTranslate(x_offset,y_offset,0)
		glRotate(20,0,1,0)
		self.font_renderer.render_string(string)
		glPopMatrix()
		#y_offset += bbox[4]-bbox[1]

		if enable_depth: glEnable(GL_DEPTH_TEST)

		glMatrixMode(GL_PROJECTION)
		glPopMatrix()
		glMatrixMode(GL_MODELVIEW)
		glDisable(GL_TEXTURE_2D)


class EMAnnotateInspector2D(QtWidgets.QWidget):
	def __init__(self,target) :
		QtWidgets.QWidget.__init__(self,None)
		self.target=weakref.ref(target)
		self.setFixedWidth(390)
		self.vbl = QtWidgets.QVBoxLayout(self)

		self.vbl.setContentsMargins(2, 2, 2, 2)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")

		# This is the tab-bar for mouse mode selection
		self.mmtab = QtWidgets.QTabWidget()

		# App tab
		#self.apptab = QtWidgets.QWidget()
		self.apptablab = QtWidgets.QTextEdit("Application specific mouse functions")
		self.mmtab.addTab(self.apptablab,"App")



		# Save tab
		self.savetab = QtWidgets.QWidget()
		self.stlay = QtWidgets.QGridLayout(self.savetab)

		self.stsnapbut = QtWidgets.QPushButton("Snapshot")
		self.stsnapbut.setToolTip(".pgm, .ppm, .jpeg, .png, or .tiff format only")
		self.stwholebut = QtWidgets.QPushButton("Save Img")
		self.stwholebut.setToolTip("save EMData in any EMAN2 format")
		self.ststackbut = QtWidgets.QPushButton("Save Stack")
		self.ststackbut.setToolTip("save EMData objects as stack in any EMAN2 format")
		self.stfpssb = QtWidgets.QSpinBox()
		self.stfpssb.setRange(1,60)
		self.stfpssb.setValue(8)
		self.stfpssb.setToolTip("Frames per second for movie generation")
		self.stmoviebut = QtWidgets.QPushButton("Movie")
		self.stanimgif = QtWidgets.QPushButton("GIF Anim")

		self.stlay.addWidget(self.stsnapbut,0,0)
		self.stlay.addWidget(self.stwholebut,1,0)
		self.stlay.addWidget(self.ststackbut,1,1)
		self.stlay.addWidget(self.stfpssb,0,1)
		self.stlay.addWidget(self.stmoviebut,0,2)
		self.stlay.addWidget(self.stanimgif,1,2)
		self.stmoviebut.setEnabled(False)
		self.stanimgif.setEnabled(False)

		self.rngbl = QtWidgets.QHBoxLayout()
		self.stlay.addLayout(self.rngbl,2,0,1,3)
		self.stmmlbl = QtWidgets.QLabel("Img Range :")
		self.stminsb = QtWidgets.QSpinBox()
		self.stminsb.setRange(0,0)
		self.stminsb.setValue(0)
		self.stmaxsb = QtWidgets.QSpinBox()
		self.stmaxsb.setRange(0,0)
		self.stmaxsb.setValue(0)

		self.rngbl.addWidget(self.stmmlbl,Qt.AlignRight)
		self.rngbl.addWidget(self.stminsb)
		self.rngbl.addWidget(self.stmaxsb)

		self.mmtab.addTab(self.savetab,"Save")

		self.stsnapbut.clicked[bool].connect(self.do_snapshot)
		self.stwholebut.clicked[bool].connect(self.do_saveimg)
		self.ststackbut.clicked[bool].connect(self.do_savestack)
		self.stmoviebut.clicked[bool].connect(self.do_makemovie)
		self.stanimgif.clicked[bool].connect(self.do_makegifanim)

		# # Filter tab
		# self.filttab = QtWidgets.QWidget()
		# self.ftlay=QtWidgets.QGridLayout(self.filttab)
		#
		# self.procbox1=StringBox(label="Process1:",value="filter.lowpass.gauss:cutoff_abs=0.125",showenable=0)
		# self.ftlay.addWidget(self.procbox1,2,0)
		#
		# self.procbox2=StringBox(label="Process2:",value="filter.highpass.gauss:cutoff_pixels=3",showenable=0)
		# self.ftlay.addWidget(self.procbox2,6,0)
		#
		# self.procbox3=StringBox(label="Process3:",value="math.linear:scale=5:shift=0",showenable=0)
		# self.ftlay.addWidget(self.procbox3,10,0)
		#
		# self.proclbl1=QtWidgets.QLabel("Image unchanged, display only!")
		# self.ftlay.addWidget(self.proclbl1,12,0)
		#
		# self.mmtab.addTab(self.filttab,"Filt")
		#
		# self.procbox1.enableChanged.connect(self.do_filters)
		# self.procbox1.textChanged.connect(self.do_filters)
		# self.procbox2.enableChanged.connect(self.do_filters)
		# self.procbox2.textChanged.connect(self.do_filters)
		# self.procbox3.enableChanged.connect(self.do_filters)
		# self.procbox3.textChanged.connect(self.do_filters)

		# Probe tab

		self.probetab = QtWidgets.QWidget()
		self.ptlay=QtWidgets.QGridLayout(self.probetab)

		self.ptareasize= ValBox(label="Probe Size:",value=32)
		self.ptareasize.setIntonly(True)
		self.ptlay.addWidget(self.ptareasize,0,0,1,2)

		self.ptpointval= QtWidgets.QLabel("Point Value (ctr pix): ")
		self.ptlay.addWidget(self.ptpointval,1,0,1,2,Qt.AlignLeft)

		self.ptareaavg= QtWidgets.QLabel("Area Avg: ")
		self.ptlay.addWidget(self.ptareaavg,2,0,Qt.AlignLeft)

		self.ptareaavgnz= QtWidgets.QLabel("Area Avg (!=0): ")
		self.ptlay.addWidget(self.ptareaavgnz,2,1,Qt.AlignLeft)

		self.ptareasig= QtWidgets.QLabel("Area Sig: ")
		self.ptlay.addWidget(self.ptareasig,3,0,Qt.AlignLeft)

		self.ptareasignz= QtWidgets.QLabel("Area Sig (!=0): ")
		self.ptlay.addWidget(self.ptareasignz,3,1,Qt.AlignLeft)

		self.ptareaskew= QtWidgets.QLabel("Skewness: ")
		self.ptlay.addWidget(self.ptareaskew,4,0,Qt.AlignLeft)

		self.ptcoord= QtWidgets.QLabel("Center Coord: ")
		self.ptlay.addWidget(self.ptcoord,4,1,Qt.AlignLeft)

		self.ptareakurt= QtWidgets.QLabel("Kurtosis: ")
		self.ptlay.addWidget(self.ptareakurt,5,0,Qt.AlignLeft)

		self.ptcoord2= QtWidgets.QLabel("( ) ")
		self.ptlay.addWidget(self.ptcoord2,5,1,Qt.AlignLeft)

		# not really necessary since the pointbox accurately labels the pixel when zoomed in
		#self.ptpixels= QtWidgets.QWidget()
		#self.ptlay.addWidget(self.ptpixels,0,2)

		self.mmtab.addTab(self.probetab,"Probe")

		# Measure tab
		self.meastab = QtWidgets.QWidget()
		self.mtlay = QtWidgets.QGridLayout(self.meastab)

		#self.mtl1= QtWidgets.QLabel("A/Pix")
		#self.mtl1.setAlignment(Qt.AlignRight)
		#self.mtlay.addWidget(self.mtl1,0,0)

		self.mtapix = ValSlider(self,label="A/Pix")
		self.mtapix.setRange(0.5,10.0)
		self.mtapix.setValue(1.0)
		self.mtlay.addWidget(self.mtapix,0,0,1,2)
#		self.mtapix.setSizePolicy(QtWidgets.QSizePolicy.Fixed,QtWidgets.QSizePolicy.Fixed)
#		self.mtapix.resize(60,21)
#		print self.mtapix.sizeHint().width(),self.mtapix.sizeHint().height()


		self.mtshoworigin= QtWidgets.QLabel("Origin: 0,0")
		self.mtlay.addWidget(self.mtshoworigin,1,0,Qt.AlignLeft)

		self.mtshowend= QtWidgets.QLabel("End: 0,0")
		self.mtlay.addWidget(self.mtshowend,1,1,Qt.AlignLeft)

		self.mtshowlen= QtWidgets.QLabel("dx,dy: 0")
		self.mtlay.addWidget(self.mtshowlen,2,0,Qt.AlignLeft)

		self.mtshowlen2= QtWidgets.QLabel("Length: 0")
		self.mtlay.addWidget(self.mtshowlen2,2,1,Qt.AlignLeft)

		self.mtshowval= QtWidgets.QLabel("Value: ?")
		self.mtlay.addWidget(self.mtshowval,3,0,1,2,Qt.AlignLeft)

		self.mtshowval2= QtWidgets.QLabel(" ")
		self.mtlay.addWidget(self.mtshowval2,4,0,1,2,Qt.AlignLeft)


		self.mmtab.addTab(self.meastab,"Meas")

		# Draw tab
		# self.drawtab = QtWidgets.QWidget()
		# self.drawlay = QtWidgets.QGridLayout(self.drawtab)
		#
		# self.dtl1 = QtWidgets.QLabel("Pen Size:")
		# self.dtl1.setAlignment(Qt.AlignRight)
		# self.drawlay.addWidget(self.dtl1,0,0)
		#
		# self.dtpen = QtWidgets.QLineEdit("5")
		# self.drawlay.addWidget(self.dtpen,0,1)
		#
		# self.dtl2 = QtWidgets.QLabel("Pen Val:")
		# self.dtl2.setAlignment(Qt.AlignRight)
		# self.drawlay.addWidget(self.dtl2,1,0)
		#
		# self.dtpenv = QtWidgets.QLineEdit("1.0")
		# self.drawlay.addWidget(self.dtpenv,1,1)
		#
		# self.dtl3 = QtWidgets.QLabel("Pen Size2:")
		# self.dtl3.setAlignment(Qt.AlignRight)
		# self.drawlay.addWidget(self.dtl3,0,2)
		#
		# self.dtpen2 = QtWidgets.QLineEdit("5")
		# self.drawlay.addWidget(self.dtpen2,0,3)
		#
		# self.dtl4 = QtWidgets.QLabel("Pen Val2:")
		# self.dtl4.setAlignment(Qt.AlignRight)
		# self.drawlay.addWidget(self.dtl4,1,2)
		#
		# self.dtpenv2 = QtWidgets.QLineEdit("0")
		# self.drawlay.addWidget(self.dtpenv2,1,3)
		#
		# self.mmtab.addTab(self.drawtab,"Draw")

		# PSpec tab
		self.pstab = QtWidgets.QWidget()
		self.pstlay = QtWidgets.QGridLayout(self.pstab)

		self.psbsing = QtWidgets.QPushButton("Single")
		self.pstlay.addWidget(self.psbsing,0,0)

		self.psbaz = QtWidgets.QPushButton("Azimuthal")
		self.pstlay.addWidget(self.psbaz,1,0)

		self.psbstack = QtWidgets.QPushButton("Stack")
		self.pstlay.addWidget(self.psbstack,0,1)

		#self.mmtab.addTab(self.pstab,"PSpec")
		self.pspecwins=[]

		# Python tab
		self.pytab = QtWidgets.QWidget()
		self.pytlay = QtWidgets.QGridLayout(self.pytab)

		self.pyinp = QtWidgets.QLineEdit()
		self.pytlay.addWidget(self.pyinp,0,0)

		self.pyout = QtWidgets.QTextEdit()
		self.pyout.setText("The displayed image is 'img'.\nEnter an expression above, like img['sigma']")
		self.pyout.setReadOnly(True)
		self.pytlay.addWidget(self.pyout,1,0,4,1)

		self.mmtab.addTab(self.pytab,"Python")

		# if self.target.external_seg_tab:
		# 	self.seg_tab = self.target.external_seg_tab
		# else:
		self.seg_tab = EMSegTab(target=target)
		self.mmtab.addTab(self.seg_tab,"Seg" )
		self.mmtab.setCurrentWidget(self.seg_tab)


		try:
			md = self.target().data.get_attr_dict()
			self.metadatatab = EMMetaDataTable(None,md)
			self.mmtab.addTab(self.metadatatab,"Meta Data")

		except:
			# the target doesn't have data yet?
			self.metadatatab = None
			pass

		self.vbl.addWidget(self.mmtab)

		# histogram level horiz layout
		self.hbl = QtWidgets.QGridLayout()
		self.hbl.setColumnStretch(290,120)
		self.hbl.setContentsMargins(10, 2, 2, 2)
		self.hbl.setSpacing(6)
		self.hbl.setObjectName("hbl")
		self.vbl.addLayout(self.hbl)

		self.hist = ImgHistogram(self)
		#self.hist.setFixedWidth(250)
		self.hist.setObjectName("hist")
		self.hbl.addWidget(self.hist,0,0,1,1)

		# Buttons next to the histogram
		self.vbl2 = QtWidgets.QGridLayout()
		self.vbl2.setColumnStretch(50,50)
		self.vbl2.setContentsMargins(2, 2, 2, 2)
		self.vbl2.setSpacing(6)
		self.vbl2.setObjectName("vbl2")
		self.hbl.addLayout(self.vbl2,0,1,1,1)

		self.invtog = QtWidgets.QPushButton("Invert")
		self.invtog.setCheckable(1)
		self.vbl2.addWidget(self.invtog,0,0,1,1)#0012


		self.histoequal = QtWidgets.QComboBox(self)
		self.histoequal.addItem("Normal")
		self.histoequal.addItem("Hist Flat")
		self.histoequal.addItem("Hist Gauss")
		self.vbl2.addWidget(self.histoequal,1,0,1,1)

		self.auto_contrast_button = QtWidgets.QPushButton("AutoC")
		self.vbl2.addWidget(self.auto_contrast_button,2,0,1,1)

		self.full_contrast_button = QtWidgets.QPushButton("FullC")
		self.vbl2.addWidget(self.full_contrast_button,3,0,1,1)

		## FFT Buttons
		#self.fftg=QtWidgets.QButtonGroup()
		#self.fftg.setExclusive(1)

		#self.ffttog0 = QtWidgets.QPushButton("Real")
		#self.ffttog0.setCheckable(1)
		#self.ffttog0.setChecked(1)
		#self.vbl2.addWidget(self.ffttog0,2,0)
		#self.fftg.addButton(self.ffttog0,0)

		#self.ffttog1 = QtWidgets.QPushButton("FFT")
		#self.ffttog1.setCheckable(1)
		#self.vbl2.addWidget(self.ffttog1,2,1)
		#self.fftg.addButton(self.ffttog1,1)

		#self.ffttog2 = QtWidgets.QPushButton("Amp")
		#self.ffttog2.setCheckable(1)
		#self.vbl2.addWidget(self.ffttog2,3,0)
		#self.fftg.addButton(self.ffttog2,2)

		#self.ffttog3 = QtWidgets.QPushButton("Pha")
		#self.ffttog3.setCheckable(1)
		#self.vbl2.addWidget(self.ffttog3,3,1)
		#self.fftg.addButton(self.ffttog3,3)

		self.scale = ValSlider(self,(0.1,5.0),"Mag:")
		self.scale.setObjectName("scale")
		self.scale.setValue(self.target().scale)
		self.vbl.addWidget(self.scale)

		self.mins = ValSlider(self,label="Min:")
		self.mins.setObjectName("mins")
		self.mins.setValue(self.target().get_minden())
		self.vbl.addWidget(self.mins)

		self.maxs = ValSlider(self,label="Max:")
		self.maxs.setObjectName("maxs")
		self.maxs.setValue(self.target().get_maxden())
		self.vbl.addWidget(self.maxs)

		self.brts = ValSlider(self,(-1.0,1.0),"Brt:")
		self.brts.setObjectName("brts")
		#self.brts.setValue(0.0)
		self.vbl.addWidget(self.brts)

		self.conts = ValSlider(self,(0.0,1.0),"Cont:")
		self.conts.setObjectName("conts")
		self.conts.setValue(1.0)
		self.vbl.addWidget(self.conts)

		self.alts = ValSlider(self,(.1,5.0),"Alt:",rounding=2)
		self.alts.setObjectName("alt")
		self.alts.setValue(self.target().get_alt())
		self.vbl.addWidget(self.alts)

		self.azs = ValSlider(self,(.1,5.0),"Az:",rounding=2)
		self.azs.setObjectName("az")
		self.azs.setValue(self.target().get_az())
		self.vbl.addWidget(self.azs)

		self.setWindowIcon(QtGui.QIcon(get_image_directory() +"eman.png"))

		self.ns = ValSlider(self,label="zpos:",rounding=2)
		nz = self.target().get_full_data()["nz"]
		zr = nz*3//4

		print("Rotate by", self.alts.getValue(),self.azs.getValue())
		self.ns.setRange(-zr,zr)
		#self.ns.setIntonly(True)

		self.ns.setValue(0)


		self.vbl.addWidget(self.ns)




		self.stmoviebut.setEnabled(True)
		self.stanimgif.setEnabled(True)

		self.lowlim=0
		self.highlim=1.0
		self.busy=0

		self.psbsing.clicked[bool].connect(self.do_pspec_single)
		self.psbstack.clicked[bool].connect(self.do_pspec_stack)
		self.psbaz.clicked[bool].connect(self.do_pspec_az)
		self.scale.valueChanged.connect(target.set_scale)
		self.mins.valueChanged.connect(self.new_min)
		self.maxs.valueChanged.connect(self.new_max)
		self.brts.valueChanged.connect(self.new_brt)
		self.conts.valueChanged.connect(self.new_cont)
		self.pyinp.returnPressed.connect(self.do_python)
		self.invtog.toggled[bool].connect(target.set_invert)
		self.histoequal.currentIndexChanged[int].connect(target.set_histogram)
		self.mmtab.currentChanged[int].connect(target.set_mouse_mode)
		self.auto_contrast_button.clicked[bool].connect(target.auto_contrast)
		self.full_contrast_button.clicked[bool].connect(target.full_contrast)
		self.alts.valueChanged.connect(self.target().set_alt)
		self.azs.valueChanged.connect(self.target().set_az)
		self.ns.valueChanged.connect(self.target().set_n)

		#self.resize(400,400) # d.woolford thinks this is a good starting size as of Nov 2008 (especially on MAC)

	# def get_seg_tab(self):
	# 	return self.seg_tab

	def update_zrange(self):
		if self.target().nz > 1:
			zr=self.target().nz*3//4
			#print(zr)		# approximate max z range with high tilt
			self.ns.setRange(-zr,zr)
			#self.ns.setValue(0)
			#self.ns.setRange(0,nz-1)
			self.stminsb.setRange(-zr,zr)
			self.stminsb.setValue(-zr)
			self.stmaxsb.setRange(-zr,zr)
			self.stmaxsb.setValue(zr)

		# self.stminsb.setRange(-zr,zr)
		# self.stminsb.setValue(-nz)
		# self.stmaxsb.setRange(-zr,zr)
		# self.stmaxsb.setValue(nz)

	def do_pspec_single(self,ign):
		"""Compute 1D power spectrum of single image and plot"""
		try:
			data=self.target().list_data[self.target().zpos]
			imgn=self.target().zpos
			lbl=f"img_{imgn}"
		except:
			data=self.target().get_data()
			lbl=time.strftime("%H:%M:%S")
		if data==None: return
#		print data
		fft=data.do_fft()
		pspec=fft.calc_radial_dist(fft["ny"]//2,0.0,1.0,1)
		ds=1.0/(fft["ny"]*data["apix_x"])
		s=[ds*i for i in range(fft["ny"]//2)]

		from .emplot2d import EMDataFnPlotter

		try:
			dfp=self.pspecwins[-1]
			dfp.set_data((s,pspec),lbl)
		except:
			dfp=EMDataFnPlotter(data=(s,pspec),key=lbl)
			self.pspecwins.append(dfp)
		dfp.show()

	def do_pspec_az(self,ign):
		"""Compute azimuthal power spectrum of single image and plot"""
		try:
			data=self.target().list_data[self.target().zpos]
			imgn=self.target().zpos
			lbl=f"img_{imgn}"
		except:
			data=self.target().get_data()
			lbl=time.strftime("%H:%M:%S")
		if data==None: return
#		print data
		fft=data.do_fft()
		pspec=fft.calc_az_dist(180,-90,1.0,4,fft["ny"]//3)
		s=np.arange(-90.0,90.0,1)

		from .emplot2d import EMDataFnPlotter

		try:
			dfp=self.pspecwins[-1]
			dfp.set_data((s,pspec),lbl)
		except:
			dfp=EMDataFnPlotter(data=(s,pspec),key=lbl)
			self.pspecwins.append(dfp)
		dfp.show()

	def do_pspec_stack(self,ign):
		"""compute average 1D power spectrum of all images and plot"""
		from .emplot2d import EMDataFnPlotter

		if len(self.target().list_data)<2 : return
		for im in self.target().list_data:
			fft=im.do_fft()
			fft.ri2inten()
			try: fftsum.add(fft)
			except: fftsum=fft

		fftsum.mult(old_div(1.0,len(self.target().list_data)))
		pspec=fftsum.calc_radial_dist(old_div(fft["ny"],2),0.0,1.0,1)	# note that this method knows about is_inten() image flag
		ds=old_div(1.0,(fft["ny"]*self.target().get_data()["apix_x"]))
		s=[ds*i for i in range(old_div(fft["ny"],2))]

		from .emplot2d import EMDataFnPlotter

		dfp=EMDataFnPlotter(data=(s,pspec))
		dfp.show()
		self.pspecwins.append(dfp)

	def do_filters(self,val=None):
		ret=[]
		if self.procbox1.getEnabled():
			try:
				nm,op=parsemodopt(self.procbox1.getValue())
				ret.append(Processors.get(nm,op))
			except:
				print("Error with processor: ",self.procbox1.getValue())

		if self.procbox2.getEnabled():
			try:
				nm,op=parsemodopt(self.procbox2.getValue())
				ret.append(Processors.get(nm,op))
			except:
				print("Error with processor2: ",self.procbox2.getValue())

		if self.procbox3.getEnabled():
			try:
				nm,op=parsemodopt(self.procbox3.getValue())
				ret.append(Processors.get(nm,op))
			except:
				print("Error with processor: ",self.procbox3.getValue())

		self.target().set_disp_proc(ret)

	def do_snapshot(self,du) :
		if self.target().data==None or self.target() == None: return
		fsp=QtWidgets.QFileDialog.getSaveFileName(self, "Select output file, .pgm, .ppm, .jpeg, .png or .tiff only")[0]
		fsp=str(fsp)
		# Just grab the framebuffer, as a QTImage, and save as tiff
		self.target().update()
		if fsp[-5:]==".jpeg" or fsp[-5:]==".tiff":
			image = self.target().grabFrameBuffer()
			image.save(fsp,str(fsp[-4:]))
			return
		if fsp[-4:]==".png" or fsp[-4:]==".tif" or fsp[-4:]==".jpg":
			image = self.target().grabFrameBuffer()
			image.save(fsp,str(fsp[-3:]))
			return
		if fsp[-4:]==".pgm" or fsp[-4:]==".ppm" : fsp=fsp[:-4]
		(depth,w,h,bmap)=self.target().render_bitmap()
		if depth==3 :
			out=open(fsp+".ppm","w")
			out.write("P6 %d %d 255\n"%(w,h))
			out.write(bmap)
			out.close()
		elif depth==1 :
			out=open(fsp+".pgm","w")
			out.write("P5 %d %d 255\n"%(w,h))
			print(w,h,w*h,len(bmap))
			out.write(bmap)
			out.close()

	def do_saveimg(self,du) :
		if self.target().data==None : return
		fsp=QtWidgets.QFileDialog.getSaveFileName(self, "Select output file, format extrapolated from file extenstion")[0]
		fsp=str(fsp)
		self.target().data.write_image(fsp,-1)

	def do_savestack(self,du) :
		if self.target().list_data==None : return
		fsp=str(QtWidgets.QFileDialog.getSaveFileName(self, "Select root output file, format extrapolated from file extenstion")[0])
		#fsp=str(fsp).split(".")
		#if len(fsp)==1 :
			#fsp1=fsp[0]
			#fsp2="mrc"
		#else :
			#fsp1=".".join(fsp[:-1])
			#fsp2=fsp[-1]

		for i in range(self.stminsb.value()-1,self.stmaxsb.value()):
			im=self.target().list_data[i]
			im.write_image(fsp,-1)

	def do_makemovie(self,du) :
		fsp=QtWidgets.QFileDialog.getSaveFileName(self, "Select output file, format extrapolated from file extenstion. ffmpeg must be installed")[0]
		if self.target().list_data==None : return

		for i in range(self.stminsb.value()-1,self.stmaxsb.value()):
			im=self.target().list_data[i]
			im["render_min"]=im["mean"]-im["sigma"]*2.5
			im["render_max"]=im["mean"]+im["sigma"]*2.5
			im.write_image("tmp.%03d.png"%(i-self.stminsb.value()+1))

		# vcodec and pix_fmt are for quicktime compatibility. -r 2 is 2 FPS
		rate=int(self.stfpssb.value())
		ret= os.system("ffmpeg -r %d -i tmp.%%03d.png -vcodec libx264 -pix_fmt yuv420p -r %d %s"%(rate,rate,fsp))
		if ret!=0 :
			QtWidgets.QMessageBox.warning(None,"Error","Movie conversion (ffmpeg) failed. Please make sure ffmpeg is in your path. Frames not deleted.")
			return

		for i in range(self.stminsb.value()-1,self.stmaxsb.value()):
			os.unlink("tmp.%03d.png"%(i-self.stminsb.value()+1))


	def do_makegifanim(self,du) :
		fsp=QtWidgets.QFileDialog.getSaveFileName(self, "Select output gif file")[0]
		if self.target().list_data==None : return

		for i in range(self.stminsb.value()-1,self.stmaxsb.value()):
			im=self.target().list_data[i]
			im.write_image("tmp.%03d.png"%(i-self.stminsb.value()+1))

		ret= os.system("gm convert tmp.???.png %s"%fsp)
		if ret!=0 :
			QtWidgets.QMessageBox.warning(None,"Error","GIF conversion failed. Please make sure GraphicsMagick ('gm convert' program) is installed and in your path. Frames not deleted.")
			return

		for i in range(self.stminsb.value()-1,self.stmaxsb.value()):
			os.unlink("tmp.%03d.png"%(i-self.stminsb.value()+1))

	def do_python(self):
		img=self.target().data
		try:
			r=str(eval(str(self.pyinp.text())))
			self.pyinp.setText("")
		except:
			import traceback
			traceback.print_exc()
			print("Error executing. Access the image as 'img', for example \nimg['mean'] will yield the mean image value")

		self.pyout.setText(str(r))

	def set_image_idx(self,val,quiet=0):
		self.ns.setValue(val,quiet=quiet)

	def get_contrast(self):
		return float(self.conts.getValue())

	def get_brightness(self):
		return float(self.brts.getValue())

	def set_maxden(self,value,quiet=1):
		self.maxs.setValue(value,quiet)

	def set_minden(self,value,quiet=1):
		self.mins.setValue(value,quiet)


	def set_scale(self,val):
		if self.busy : return

		self.busy=1
		self.scale.setValue(val)
		#print("Scale in set_scale", self.scale.value)
		self.busy=0

	def new_min(self,val):
		if self.busy : return
		self.busy=1
		self.target().set_density_min(val)
		self.update_brightness_contrast()
		self.busy=0

	def new_max(self,val):
		if self.busy : return
		self.busy=1
		self.target().set_density_max(val)
		self.update_brightness_contrast()
		self.busy=0

	def new_brt(self,val):
		if self.busy : return
		self.busy=1
		self.update_min_max()
		self.busy=0

	def new_cont(self,val):
		if self.busy : return
		self.busy=1
		self.update_min_max()
		self.busy=0


	def update_brightness_contrast(self):
		b=0.5*(self.mins.value+self.maxs.value-(self.lowlim+self.highlim))/((self.highlim-self.lowlim))
		c=old_div((self.mins.value-self.maxs.value),(2.0*(self.lowlim-self.highlim)))
		brts = -b
		conts = 1.0-c
		self.brts.setValue(brts,1)
		self.conts.setValue(conts,1)

	def update_min_max(self):
		x0=(old_div((self.lowlim+self.highlim),2.0)-(self.highlim-self.lowlim)*(1.0-self.conts.value)-self.brts.value*(self.highlim-self.lowlim))
		x1=(old_div((self.lowlim+self.highlim),2.0)+(self.highlim-self.lowlim)*(1.0-self.conts.value)-self.brts.value*(self.highlim-self.lowlim))
		self.mins.setValue(x0,1)
		self.maxs.setValue(x1,1)
		self.target().set_density_range(x0,x1)

	def set_hist(self,hist,minden,maxden):
		if hist != None and len(hist) != 0:self.hist.set_data(hist,minden,maxden)

		# We also update the image header info, since it is coordinated
		d=self.target().get_data()
		if d==None: return
		try:
			ctfdef=d["ctf"].defocus
			defocus=f"&Delta;Z={ctfdef:1.5g}"
		except: defocus=""

		try: ptclrepr=f"ptcl_repr={d['ptcl_repr']}"
		except: ptclrepr=""

		header=f'''<html><body>
<p>Mouse buttons -> Application</p>
<p/>
<table style="width: 300">
<tr><td width="80">nx={d["nx"]:d}</td><td width="120">min={d["minimum"]:1.4g}</td><td width="120">apix_x={d["apix_x"]:1.3f}</td></tr>
<tr><td width="80">ny={d["ny"]:d}</td><td width="120">max={d["maximum"]:1.4g}</td><td width="120">apix_y={d["apix_y"]:1.3f}</td></tr>
<tr><td width="80">nz={d["nz"]:d}</td><td width="120">mean={d["mean"]:1.5g}</td><td width="120">apix_z={d["apix_z"]:1.3f}</td></tr>
<tr><td width="80">{defocus}</td><td width="120">sigma={d["sigma"]:1.5g}</td><td>{ptclrepr}</td></tr>
</table></html>
'''
		self.apptablab.setHtml(header)


	def set_bc_range(self,lowlim,highlim):
		#print "set range",lowlim,highlim
		self.lowlim=lowlim
		self.highlim=highlim
		self.mins.setRange(lowlim,highlim)
		self.maxs.setRange(lowlim,highlim)

	def set_limits(self,lowlim,highlim,curmin,curmax):
		if highlim<=lowlim : highlim=lowlim+.001
		#print "in set limits", self.conts.getValue(), self.conts.getValue()
		self.lowlim=lowlim
		self.highlim=highlim
		self.mins.setRange(lowlim,highlim)
		self.maxs.setRange(lowlim,highlim)
		self.mins.setValue(curmin,1)
		self.maxs.setValue(curmax,1)
		#print "leave set limits", self.conts.getValue(), self.conts.getValue()

class CustomTreeSet(QtWidgets.QTreeWidget):
	def __init__(self, target = None):
		super().__init__()
		self.target = target

	def dropEvent(self,event):

		super().dropEvent(event)
		#print("DROPPING")
		self.target.update_sets()
		#self.target.target.tree_sels = self.target.get_whole_branch(item)
		self.target.target.force_display_update(set_clip=False)
		self.target.target.updateGL()

class CustomTreeWidgetItem(QtWidgets.QTreeWidgetItem):
	def __init__(self, color = None):
		super().__init__()
		if color:
			self.color = color




		#


class EMSegTab(QtWidgets.QWidget):
	'''
	This is the set display panel

	'''
	def __init__(self, target = None):
		super().__init__()

		self.eraser=QtWidgets.QRadioButton("Eraser")
		self.classes_cb = QtWidgets.QRadioButton("Classes")
		self.classes_cb.setChecked(True)
		self.cb_group = QtWidgets.QButtonGroup()
		self.cb_group.addButton(self.classes_cb,1)
		self.cb_group.addButton(self.eraser,2)
		self.pen_width=ValSlider(label="Pen Sz",labelwidth=30,value=15,rounding=0,rng=(0,60))
		self.pen_width.setIntonly(1)
		self.pen_width.setEnabled(1)
		self.target = target
		self.og = None #to store annotation image before grouping

		#self.background=False



		# self.set_list=QtWidgets.QListWidget()
		# self.set_list.setSizePolicy(QtWidgets.QSizePolicy.Preferred,QtWidgets.QSizePolicy.Expanding)
		# self.set_list.addItem("Void")
		# self.set_list.item(0).setHidden(True)

		#self.tree_set=QtWidgets.QTreeWidget(self)
		self.tree_set=CustomTreeSet(target=self)
		self.tree_set.setColumnCount(3)
		self.tree_set.setHeaderLabels(['Index','Class Name','Group'])
		self.tree_set.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
		self.tree_set.setDragEnabled(True)
		self.tree_set.setDragDropMode(QtWidgets.QAbstractItemView.InternalMove)
		self.tree_root = self.tree_set.invisibleRootItem()
		# self.tree_set.setColumnWidth(0,60)
		# self.tree_set.setColumnWidth(1,120)
		# self.tree_set.setColumnWidth(2,60)

		#REMOVE TABLE_SET
		# self.table_set=QtWidgets.QTableWidget(0, 3, self)
		# #self.table_set=QtWidgets.QTreeWidget(self)
		# self.table_set.setColumnWidth(0,60)
		# self.table_set.setColumnWidth(1,120)
		# self.table_set.setColumnWidth(2,60)
		# self.table_set.setHorizontalHeaderLabels(['Index','Class Name','Group'])
		# for i in range(3):
		# 	self.table_set.horizontalHeaderItem(i).setTextAlignment(Qt.AlignLeft)


		#self.itemflags = Qt.ItemFlags(Qt.ItemIsEditable)|Qt.ItemFlags(Qt.ItemIsSelectable)|Qt.ItemFlags(Qt.ItemIsEnabled)|Qt.ItemFlags(Qt.ItemIsUserCheckable)
		self.itemflags = Qt.ItemFlags(Qt.ItemIsSelectable)|Qt.ItemFlags(Qt.ItemIsEnabled)|Qt.ItemFlags(Qt.ItemIsUserCheckable)|Qt.ItemFlags(Qt.ItemIsDragEnabled) |Qt.ItemFlags(Qt.ItemIsDropEnabled)
		self.indexflags = Qt.ItemFlags(Qt.ItemIsEnabled)|Qt.ItemFlags(Qt.ItemIsUserCheckable)|Qt.ItemFlags(Qt.ItemIsSelectable)
		#self.used_group_index = [100]
		self.used_index=[0]
		#self.nodes = {}
		#self.colors = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000']
		self.colors = self.target.get_color_palette()

		#self.read_header(self.target.get_full_annotation())



		self.browser_ret = None
		button_vbl = QtWidgets.QVBoxLayout()

		self.new_class_button = QtWidgets.QPushButton("New Class")
		button_vbl.addWidget(self.new_class_button)
		self.group_button = QtWidgets.QPushButton("Group")
		button_vbl.addWidget(self.group_button)
		self.ungroup_button = QtWidgets.QPushButton("Ungroup")
		button_vbl.addWidget(self.ungroup_button)
		self.delete_sel_button = QtWidgets.QPushButton("Delete Selected")
		button_vbl.addWidget(self.delete_sel_button)
		self.load_mask_button = QtWidgets.QPushButton("Load Mask")
		button_vbl.addWidget(self.load_mask_button)
		self.load_class_button = QtWidgets.QPushButton("Load Class")
		button_vbl.addWidget(self.load_class_button)
		self.save_mask_button = QtWidgets.QPushButton("Save Mask")
		button_vbl.addWidget(self.save_mask_button)
		self.append_ann_button = QtWidgets.QPushButton("Append Ann")
		button_vbl.addWidget(self.append_ann_button)
		self.load_all_button = QtWidgets.QPushButton("Load All")
		button_vbl.addWidget(self.load_all_button)
		self.save_all_button = QtWidgets.QPushButton("Save All")
		button_vbl.addWidget(self.save_all_button)

		#For testing only
		self.test_widget_button = QtWidgets.QPushButton("Test Button")
		button_vbl.addWidget(self.test_widget_button)

		hbl = QtWidgets.QHBoxLayout()
		#hbl.addWidget(self.table_set)
		hbl.addWidget(self.tree_set)
		#hbl.addWidget(self.set_list)
		hbl.addLayout(button_vbl)

		self.button_list = [self.load_all_button,self.save_mask_button,self.load_mask_button,self.new_class_button,self.delete_sel_button]


		segtab_vbl = QtWidgets.QVBoxLayout(self)
		eraser_lay = QtWidgets.QHBoxLayout()
		eraser_lay.addWidget(self.classes_cb)
		eraser_lay.addWidget(self.eraser)
		eraser_lay.addWidget(self.pen_width)

		segtab_vbl.addLayout(eraser_lay)
		segtab_vbl.addLayout(hbl)


		self.new_class_button.clicked[bool].connect(self.new_class)
		self.delete_sel_button.clicked[bool].connect(self.delete_sel)
		self.load_mask_button.clicked[bool].connect(self.load_mask)
		self.load_class_button.clicked[bool].connect(self.load_class)
		self.save_mask_button.clicked[bool].connect(self.save_mask)
		self.append_ann_button.clicked[bool].connect(self.append_ann)
		self.load_all_button.clicked[bool].connect(self.load_all)
		self.save_all_button.clicked[bool].connect(self.save_all)
		self.cb_group.buttonClicked[QtWidgets.QAbstractButton].connect(self.on_check_cb_group)
		# self.table_set.cellClicked[int,int].connect(self.on_table_set)

		self.group_button.clicked[bool].connect(self.group_sel)
		self.ungroup_button.clicked[bool].connect(self.ungroup_sel)
		self.test_widget_button.clicked[bool].connect(self.test_widget_button_clicked)

		self.tree_set.itemClicked[QtWidgets.QTreeWidgetItem,int].connect(self.on_tree_item_clicked)
		self.tree_set.itemDoubleClicked[QtWidgets.QTreeWidgetItem,int].connect(self.on_tree_item_double_clicked)

		self.tree_set.currentItemChanged[QtWidgets.QTreeWidgetItem,QtWidgets.QTreeWidgetItem].connect(self.on_tree_item_changed)
		self.tree_set.itemCollapsed[QtWidgets.QTreeWidgetItem].connect(self.on_tree_item_collapsed)
		#self.tree_set.itemChanged[QtWidgets.QTreeWidgetItem,int].connect(self.update_sets)









	def get_current_class(self):
		#Return the index of current class
		if self.cb_group.checkedId() == 1:
			#sels = self.table_set.selectedItems()
			sels = self.tree_set.selectedItems()
			if len(sels) == 0:
				print("No class selected. Use manual annotate panels to create class")
				return 0
			else:
				# row = self.table_set.row(sels[0])
				# if int(self.table_set.item(row,2).text()) == -1:
				# 	return int(self.table_set.item(row,0).text())
				# else:
				# 	return int(self.table_set.item(row,2).text())+100
				try:
					return self.tree_set.currentItem().text(0)
				except:
					print("No class selected. Use manual annotate panels to create class")
					return 0
		else:
			return 0

	def get_pen_width(self):
		#Return the pen width
		return int(self.pen_width.value)

	def on_check_cb_group(self,cb):
		#Return if class or eraser is selected
		print(cb.text()+" is selected")
		if cb.text() == "Eraser":
			# if self.table_set.rowCount() == 0:
			# 	print("There are no annotation to erase")
			if self.tree_set.topLevelItemCount() == 0:
				print("There are no annotation to erase")
			else:
				#self.table_set.currentItem().setSelected(False)
				try:
					self.tree_set.currentItem().setSelected(False)
				except:
					pass
			for button in self.button_list:
				button.setEnabled(False)
		else:
			#self.table_set.setEnabled(True)
			self.tree_set.setEnabled(True)
			# try:
			# 	#self.table_set.item(0,1).setSelected(True)
			# 	#self.table_set.topLevelItem(0).setSelected(True)
			# except:
			# 	pass
			for button in self.button_list:
				button.setEnabled(True)
		self.target.force_display_update()
		self.target.updateGL()

	def on_tree_item_clicked(self,item,col):
		#print("Row", row, "Col", col)
		self.classes_cb.setChecked(True)
		for button in self.button_list:
			button.setEnabled(True)
		#self.target.tree_sels = self.get_whole_branch(item)
		self.target.force_display_update()
		self.target.updateGL()

	def get_selected_item(self):
		sels = self.tree_set.selectedItems()
		sel_l = []
		for sel in sels:
			# print("Sel",sel.text(0))
			# print(self.get_whole_branch(sel))
			for item in self.get_whole_branch(sel):
				sel_l.append(item)
				#print("Get_selected_item",item.text(0))
		return sel_l




	# def on_tree_item_changed(self,item,col):
	# 	#print("Row", row, "Col", col)
	# 	self.update_sets()
	# 	self.target.tree_sels = self.get_whole_branch(item)
	# 	self.target.force_display_update()
	# 	self.target.updateGL()
	#
	# def get_whole_branch(self,item):
	# 	if item.childCount() == 0:
	# 		#print(item.text(0))
	# 		return [item]
	# 	else:
	# 		group_val = int(item.text(0))
	# 		#print("Group Val")
	# 		item_l = [item]
	# 		it = QtWidgets.QTreeWidgetItemIterator(item,flags=QtWidgets.QTreeWidgetItemIterator.All)
	# 		it +=1
	# 		#it.setFlags()
	# 		# Qt.ItemFlags(Qt.ItemIsSelectable)
	# 		# it.QTreeWidgetItemIterator::Selected
	# 		while it.value():
	# 			#if int(it.value().text(2)) == group_val:
	#
	# 			item_l.append(it.value())
	# 			it +=1
	#
	# 			try:
	# 				print("Parent",it.value().parent().text(0))
	# 				if it.value().text(2) == "-1" or it.value().parent().text(0) == item.parent().text(0):
	# 					break
	# 			except:
	# 				continue
	# 		return item_l

	def get_whole_branch(self,tree_widget_item):
		"""Returns all QTreeWidgetItems in the subtree rooted at the given node."""
		nodes = []
		nodes.append(tree_widget_item)
		for i in range(tree_widget_item.childCount()):
			nodes.extend(self.get_whole_branch(tree_widget_item.child(i)))
		return nodes









	def on_tree_item_double_clicked(self,item,col):
		index = int(item.text(0))
		#print("Row", row, "Col", col)
		if col == 1:
			self.tree_set.openPersistentEditor(item,col)
		elif col == 0:
			#name,ok=QtWidgets.QInputDialog.getText( self, "Class Name", "Enter name for the new class:")

			color_dialog = QtWidgets.QColorDialog()
			color = color_dialog.getColor()
			if not color.isValid():
				 return
			#print(color.hue())

			self.target.colors[index] = QtGui.QColor.fromHsv(color.hue(),color.saturation(),255)

			self.update_sets()
			self.target.need_new_RGB = 1
			self.target.force_display_update()
			self.target.updateGL()

			return
			#self.update_sets()
		else:
			return
	def on_tree_item_changed(self,current,prev):
		self.tree_set.closePersistentEditor(prev,1)
		self.update_sets()
		# try:
		# 	self.target.tree_sels = self.get_whole_branch(current)
		# except:
		# 	self.target.tree_sels = []
		self.target.force_display_update()
		self.target.updateGL()
	# def keyPressEvent(self,event):
	# 	if event.key() == Qt.Key_Enter:
	# 		self.tree_set.closePersistentEditor(item,col)
	# 	else:
	# 		self.target.keypress.emit(event)



	def on_tree_item_collapsed(self,item):
		self.tree_set.setCurrentItem(None)
		self.target.force_display_update()
		self.target.updateGL()


		return


	def get_unused_index(self):
		# if not group:
			# n = self.table_set.rowCount()
			# self.used_index = []
			# for i in range(n):
			# 	self.used_index.append(int(self.table_set.item(i,0).text()))
			# try:
			# 	self.unused = (set(range(min(self.used_index)+1, max(self.used_index)+2)) - set(self.used_index))
			# 	new_index = sorted(list(self.unused))[0]
			# except:
			# 	new_index = 1
			# return new_index
			# n = self.tree_set.topLevelItemCount()
			# self.used_index = []
			# for i in range(n):
			# 	self.used_index.append(int(self.tree_set.topLevelItem(i).text(0)))
		try:
			self.unused = (set(range(min(self.used_index)+1, max(self.used_index)+2)) - set(self.used_index))
			new_index = sorted(list(self.unused))[0]
		except:
			new_index = 1
		return new_index
		# else:
		# 	self.used_group_index.append(self.used_group_index[-1]+1)
		# 	return self.used_group_index[-1]


	def new_class(self):
		#Add a new class to the table set and annotation file
		name,ok=QtWidgets.QInputDialog.getText( self, "Class Name", "Enter name for the new class:")
		if not ok : return
		new_index = self.get_unused_index()
		self.add_new_row(new_index,str(name))
		#self.nodes[new_index] = Node(new_index,new_index)
		#self.table_set.setCurrentCell(self.table_set.rowCount()-1,1)

		self.tree_set.setCurrentItem(self.tree_set.topLevelItem(self.tree_set.topLevelItemCount()-1))
		self.new_added = 1
		self.update_sets()


	def test_widget_button_clicked(self):

		#self.tree_set.currentItem().addChild(QtWidgets.QTreeWidgetItem(["17","b","-1"]))
		#self.add_child(child_l=["17","b","-1"])
		self.color_label = QtWidgets.QLabel()
		self.color_label.setGeometry(100, 100, 200, 60)
		#self.color_label.setAutoFillBackground(True)
		self.color_label.show()
		color_dialog = QtWidgets.QColorDialog()
		color_dialog.setOption(QtWidgets.QColorDialog.NoButtons)
		color = color_dialog.getColor()

		# setting graphic effect to the label
		# graphic = QtWidgets.QGraphicsColorizeEffect(self)
		# graphic.setColor(self.color)
		alpha  = 140
		values = "{r}, {g}, {b}".format(r = color.hue(),
												g = color.saturation(), b = 120
												#b = color.value()
												#a = alpha
												)
		self.color_label.setText(str(values))
		self.color_label.setStyleSheet("QLabel { background-color: hsv("+values+");  }")
		# a = self.target.create_RGB_list(3).numpy()
		# print(a, len(a))




	def add_child(self,child_item=None,child_l=[]):
		if child_item:
			self.tree_set.currentItem().addChild(child_item)
		elif len(child_l) > 0:
			child_item = QtWidgets.QTreeWidgetItem(child_l)
			self.tree_set.currentItem().addChild(child_item)
		else:
			return





	def change_parent(self,item, new_parent):
		old_parent = item.parent()
		#print("Old parent",old_parent)
		try:
			ix = old_parent.indexOfChild(item)
			item_without_parent = old_parent.takeChild(ix)
		except:
			self.tree_root.addChild(item)
			ix = self.tree_root.indexOfChild(item)
			item_without_parent = self.tree_root.takeChild(ix)
		new_parent.addChild(item_without_parent)
		# if item_without_parent.childCount() > 0:
		# it = QtWidgets.QTreeWidgetItemIterator(item_without_parent)
		#
		# #it +=1
		# while it.value():
		# 	#print(item.text(0))
		# 	item = it.value()
		# 	#print("CP",item.text(0))
		# 	if new_parent.text(0) == "":
		# 		print("EMPTY")
		# 		item.setText(2,"-1")
		# 	else:
		# 		print("NOT EMPTY")
		# 		item.setText(2,new_parent.text(0))
		# 	it += 1



	def color_by_group(self):
		temp = self.target.get_full_annotation().copy_head()
		ori = EMData('./segs/temp/single_temp.hdf')
		self.target.full_annotation = ori
		it = QtWidgets.QTreeWidgetItemIterator(self.tree_root)
		while it.value():
			item = it.value()
			print("item", item.text(1),"and",item.text(2))
			try:
				group_val = int(item.text(2))
			except:
				group_val = -1
			if group_val == -1:
				# print("No parent")
				# print(item.text(1))
				it += 1
			else:
				self.color_tree_item(item)
				val = int(item.text(0))
				print("Group_val",group_val)
				temp += group_val*(ori.process("threshold.binaryrange",{"high":val+0.1,"low":val-0.1}))
				self.target.get_full_annotation().process_inplace("threshold.rangetozero",{"maxval":(val+0.1),"minval":(val-0.1)})
				it += 1
		self.target.full_annotation += temp
		self.target.force_display_update()
		self.update_sets()
		self.target.updateGL()

	def group_sel(self):

		#self.target.get_full_annotation().write_image('temp_temp.hdf')
		# if len(self.nodes)==0:
		# 	self.make_nodes()


		#group_index=self.get_unused_index(group=True)
		group_index=self.get_unused_index()
		#group_node = Node(group_index,group_index)
		#self.nodes[group_index] = group_node

		temp = self.target.get_full_annotation().copy_head()
		# sels = self.table_set.selectedItems()
		# if len(sels) == 0:
		# 	print("Must select class to group")
		# 	return
		# for sel in sels:
		# 	row = self.table_set.row(sel)
		# 	val = int(self.table_set.item(row,0).text())
		# 	group_node.add_child(self.nodes[val])
			#print("add child at row", row, " to group node")
			#self.table_set.item(row,2).setText(str(group_index%100))

		sels = self.tree_set.selectedItems()
		if len(sels) == 0:
			print("Must select class to group")
			return
		else:
			name,ok=QtWidgets.QInputDialog.getText( self, "Group Name", "Enter name for group:")
			if not ok : return


		if not os.path.exists('./segs/temp') :
			os.mkdir('./segs/temp')

		#self.target.get_full_annotation().write_image('./segs/temp/temp_'+str(group_index)+'.hdf')
		if not os.path.isfile('./segs/temp/single_temp.hdf'):
			self.target.get_full_annotation().write_image('./segs/temp/single_temp.hdf')

		self.add_new_row(group_index,name)
		group_item = self.tree_set.topLevelItem(self.tree_set.topLevelItemCount()-1)
		for sel in sels:
			sel.setText(2,str(group_index))
			self.change_parent(sel,group_item)

			# group_item.addChild(sel)
			# print("parent",sel.parent())
			# row = self.table_set.row(sel)
			# val = int(self.table_set.item(row,0).text())
			# group_node.add_child(self.nodes[val])
		# it = QtWidgets.QTreeWidgetItemIterator(self.tree_set)
		# while it.value():
		# 	item = it.value()
		# 	if int(item.text(2)) == -1:
		# 		print("No parent")
		# 		print(item.text(1))
		# 		it += 1
		# 	else:
		# 		self.color_tree_item(item)
		# 		val = int(item.text(0))
		# 		group_val = int(item.text(2))
		# 		temp += (self.target.get_full_annotation().process("threshold.binaryrange",{"high":val+0.1,"low":val-0.1}))
		# 		self.target.get_full_annotation().process_inplace("threshold.rangetozero",{"maxval":(val+0.1),"minval":(val-0.1)})
		# 		it += 1
		#
		# self.target.full_annotation += temp*group_index
		#self.color_by_group()
		#self.target.force_display_update()
		#self.add_new_row(group_index,str("GROUP ")+str(group_index%100))
		#self.update_group_nums()
		self.update_sets()


		#self.target.display_group.append(group_index)
		#self.recolor_by_group()
		#self.target.tree_sels = []
		self.target.force_display_update(set_clip=False)
		self.target.updateGL()

		return


	def ungroup_sel(self):


		# sels = self.table_set.selectedItems()
		# if len(sels) == 0 or len(sels) >1:
		# 	print("Ungroup once at a time")
		# 	return
		# sel_row = self.table_set.row(sels[0])
		# if not self.table_set.item(sel_row,1).text().startswith('GROUP'):
		# 	print("Select a group item to ungroup")
		# 	return
		sels = self.tree_set.selectedItems()
		if len(sels) == 0 or len(sels) >1:
			print("Ungroup once at a time")
			return
		sel = sels[0]
		if sel.childCount() == 0:
			print("Not a group node")
			return
		else:
			#group_val = int(sel.text(0))
			# ori = EMData('./segs/temp/temp_'+str(group_val)+'.hdf')
			# temp = self.target.get_full_annotation().copy_head()
			for orphan in sel.takeChildren():


				#self.change_parent(orphan,self.tree_root)

				self.change_parent(orphan,sel.parent() or self.tree_root)
					#orphan.setText(2,(orphan.parent().text(0)))

					#orphan.setText(2, "-1")
				#self.color_by_group()

			# 	val = int(orphan.text(0))
			# 	temp += (val-group_val)*(ori.process("threshold.binaryrange",{"high":val+0.1,"low":val-0.1}))
			#
			# #temp  += self.target.get_full_annotation().process("threshold.rangetozero",{"maxval":(group_val+0.1),"minval":(group_val-0.1)})
			# temp  += self.target.get_full_annotation()
			#
			# self.target.full_annotation = temp
			#
			# #self.tree_root.removeChild(sel)
			# del temp,ori
			# 	return


		#group_index = int(self.table_set.item(sel_row,0).text())
		#group_index = int(sel.text(0))




		#self.table_set.removeRow(sel_row)

		# for child in self.nodes[group_index].children:
		# 	child.make_orphan()
		# self.update_group_nums()
		# try:
		# 	self.used_group_index.remove(group_index)
		# 	self.target.display_group.remove(group_index)
		# except:
		# 	pass
		self.update_sets()
		#self.recolor_by_group()
		self.target.force_display_update()
		self.target.updateGL()



	# def recolor_by_group(self):
	# 	self.target.full_annotation=EMData('./segs/temp_temp.hdf')
	# 	temp = self.target.get_full_annotation().copy_head()
	# 	for row in range(self.table_set.rowCount()):
	# 		index = int(self.table_set.item(row,0).text())
	# 		if index <100:
	# 			val = self.nodes[index].get_value()
	# 			temp += self.target.get_full_annotation().process("threshold.binaryrange",{"high":index+0.1,"low":index-0.1})*(val)
	# 			self.target.get_full_annotation().process_inplace("threshold.rangetozero",{"maxval":(index+0.1),"minval":(index-0.1)})
	#
	# 	self.target.full_annotation += temp
	# 	self.target.force_display_update()
	# 	#self.target.force_display_update(set_clip=False)
	# 	self.target.updateGL()
	# 	self.update_sets()
	# 	del temp



	# def make_nodes(self):
	# 	rows = self.table_set.rowCount()
	# 	self.nodes = {}
	# 	for row in range(rows):
	# 		index = int(self.table_set.item(row,0).text())
	# 		node = Node(index,index)
	# 		self.nodes[index] = node
	# 		#node.print_tree()

	# def update_group_nums(self):
	# 	for row in range(self.table_set.rowCount()):
	# 		index = int(self.table_set.item(row,0).text())
	# 		val  = self.nodes[index].get_value()%100
	# 		if self.nodes[index].parent:
	# 			self.table_set.item(row,2).setText(str(val))
	# 		else:
	# 			self.table_set.item(row,2).setText("-1")



	def delete_sel(self):

		#Remove classes to the table set and annotation file. The deleted classes can be saved into separated annotation file.
		#sels = self.table_set.selectedItems()
		sels = []
		for sel in self.tree_set.selectedItems():
			for item in self.get_whole_branch(sel):

				sels.append(item)
		#print("Length", len(sels))


		if len(sels) == 0:
			print("Must select class to delete")
			return
		#sels_name = [sel.text() for sel in sels]
		sels_name = [sel.text(1) for sel in sels]
		msg_box=QtWidgets.QMessageBox()
		msg_box.setText("The class "+",".join(sels_name)+" will be deleted")
		msg_box.setInformativeText("Do you want to save these classes to disk as an annotation file?")
		msg_box.setStandardButtons(QtWidgets.QMessageBox.Save |  QtWidgets.QMessageBox.Cancel | QtWidgets.QMessageBox.Discard )
		msg_box.setDefaultButton(QtWidgets.QMessageBox.Save)
		msg_box.exec_()
		if msg_box.clickedButton().text() == "Save":
			self.save_mask(multiple_class = True)
		elif msg_box.clickedButton().text() == "Cancel":
			return
		else:
			pass
		for sel in sels:

		#for sel in self.get_whole_branch(sel_node):
			val = int(sel.text(0))
			(sel.parent() or self.tree_root).removeChild(sel)
			self.target.get_full_annotation().process_inplace("threshold.rangetozero",{"maxval":(val+0.1),"minval":(val-0.1)})


		self.target.force_display_update(set_clip=False)
		self.target.updateGL()
		self.update_sets()

	def load_mask(self):
		#Initiate loading a binary mask on the current annotation by popping up a EMBrowserWidget
		self.openbrowser = EMBrowserWidget(withmodal=True,multiselect=False)
		self.openbrowser.ok.connect(self.load_mask_browser_ok)
		self.openbrowser.show()
		return

	def load_mask_browser_ok(self):
		#Method to actual load a binary mask selected from the browser on the current annotation
		self.browser_ret = (self.openbrowser.getResult())
		print(self.browser_ret)
		in_f = EMData(self.browser_ret[0]).process("threshold.binary",{"value":0.3})
		self.update_sets()
		self.target.full_annotation *= in_f
		self.target.force_display_update()
		self.target.updateGL()
		return

	def load_class(self):
		#Initiate loading a binary mask on the current annotation by popping up a EMBrowserWidget
		self.openbrowser = EMBrowserWidget(withmodal=True,multiselect=False)
		self.openbrowser.ok.connect(self.load_class_browser_ok)
		self.openbrowser.show()
		return

	def load_class_browser_ok(self):
		#Method to actual load a binary mask selected from the browser as the selected class
		sels = self.tree_set.selectedItems()
		if len(sels) == 0:
			print("Must select class to delete")
			return

		val = int(sels[0].text(0))


		self.browser_ret = (self.openbrowser.getResult())
		print(self.browser_ret)
		in_f = EMData(self.browser_ret[0]).process("threshold.binary",{"value":0.3})
		if in_f.get_sizes() != self.target.full_data.get_sizes():
			print("Annotation file must have the same dimension with the data")
			return
		else:
			self.target.full_annotation += in_f*val
			self.target.force_display_update()
			self.target.updateGL()
		return


	def get_selected_annotate(self):
		sels = self.tree_set.selectedItems()
		if len(sels) == 0:
			print("Must select class to save")
			return
		annotation_out=self.target.get_full_annotation().copy_head()
		annotation_out.to_zero()
		for sel in sels:
			num = int(sel.text(0))
			annotation_out += (self.target.get_full_annotation().process("threshold.binaryrange",{"high":num+0.1,"low":num-0.1}))
		return annotation_out

	def get_group_annotate(self):
		sels = self.tree_set.selectedItems()
		if len(sels) == 0 or len(sels) > 1:
			print("Must select single group")
			return
		if sels[0].childCount() == 0:
			print("Not a group node")
			return self.get_selected_annotate()

		annotation_out=self.target.get_full_annotation().copy_head()
		annotation_out.to_zero()
		for sel in self.get_whole_branch(sels[0]):
			num = int(sel.text(0))
			annotation_out += (self.target.get_full_annotation().process("threshold.binaryrange",{"high":num+0.1,"low":num-0.1}))
		return annotation_out

	def save_mask(self, multiple_class = False, ret = False):
		#Save the selected annotations as a binary mask (default)
		#or a multilabel annotation file (to be called when calling delete sel method)
		#sels = self.table_set.selectedItems()
		#sels = self.tree_set.selectedItems()
		sels = []
		for sel in self.tree_set.selectedItems():

			for item in self.get_whole_branch(sel):

				sels.append(item)
		#print("Length", len(sels))

		if len(sels) == 0:
			print("Must select class to save")
			return
		out_name,ok=QtWidgets.QInputDialog.getText( self, "Save Selected", "Save selected annotation classes")
		if not ok : return

		nums = []
		names = []
		annotation_out=self.target.get_full_annotation().copy_head()
		annotation_out.to_zero()
		for sel in sels:
			num = int(sel.text(0))
			if multiple_class:
				annotation_out += num*(self.target.get_full_annotation().process("threshold.binaryrange",{"high":num+0.1,"low":num-0.1}))
			else:
				annotation_out += (self.target.get_full_annotation().process("threshold.binaryrange",{"high":num+0.1,"low":num-0.1}))

			nums.append(num)
			# names.append(str(self.table_set.item(row,1).text()))
			names.append(sel.text(1))
		if ret:
			return annotation_out
		else:

			name_str=names[0]+""
			#for i in range(1,len(names)):
				#name_str = name_str + "," + names[i]
			serialize_name = json.dumps(names, default=lambda a: "[%s,%s]" % (str(type(a)), a.pk))
			annotation_out["ann_name"] = serialize_name
			annotation_out["ann_num"] = nums
			annotation_out.write_image(out_name)
			print("Annotation is saved to", out_name)
			return

	def append_ann(self):
		#Initiate loading an annotation file to append on the current annotation
		self.openbrowser = EMBrowserWidget(withmodal=True,multiselect=False)
		self.openbrowser.ok.connect(self.append_ann_browser_ok)
		self.openbrowser.show()
		return

	def append_ann_browser_ok(self):
		#Method to actually load an annotation file to append on the current annotation.
		#Helpful when user want to display multiple annotation on tomograms. May looks weird if the annotations files have overlays.
		self.browser_ret = (self.openbrowser.getResult())
		#print(self.browser_ret)
		in_f = EMData(self.browser_ret[0])
		if in_f.get_sizes() != self.target.full_data.get_sizes():
			print("Annotation file must have the same dimension with the data")
			return
		else:
			try:
				append_dict = self.read_header(in_f,ret=True)
				for key, value in append_dict.items():
					self.add_new_row(key,value)
				self.update_sets()
			except:
				pass
			self.target.full_annotation += in_f
			self.target.force_display_update()
			self.target.updateGL()
			return




	def load_all(self):
		#Initiate loading a new annotation file to display.
		self.openbrowser = EMBrowserWidget(withmodal=True,multiselect=False)
		self.openbrowser.ok.connect(self.load_all_browser_ok)
		self.openbrowser.show()
		return

	def load_all_browser_ok(self):
		#Method to actually load a new annotation file to display.
		self.browser_ret = (self.openbrowser.getResult())
		inf = EMData(self.browser_ret[0])
		if inf.get_sizes() != self.target.full_data.get_sizes():
			print("Annotation file must have the same dimension with the data")
			return
		else:
			# row_count = self.table_set.rowCount()
			# for i in range(row_count):
			# 	self.table_set.removeRow(row_count - i - 1)
			self.tree_set.clear()
			self.target.full_annotation=inf
			self.target.force_display_update(set_clip=0)
			#self.read_header(self.target.get_full_annotation())
			self.update_sets()
			self.target.updateGL()
			return



	def save_all(self):
		#Save the current annotation to disk as a multiclass annotation file.
		out_name,ok=QtWidgets.QInputDialog.getText( self, "Save Full Annotation", "Save full annotation to:")
		if not ok : return

		# nums = [int(self.table_set.item(row,0).text()) for row in range(self.table_set.rowCount())]
		# names = [str(self.table_set.item(row,1).text()) for row in range(self.table_set.rowCount())]
		# group_id = [int(self.table_set.item(row,2).text()) for row in range(self.table_set.rowCount())]
		# #name_str=names[0]+""
		# #for i in range(1,len(names)):
		# 	#name_str = name_str + "," + names[i]
		#
		# #TOREAD
		# # self.xform = Transform({"type":"eman","tx":self.full_data["nx"]//2,"ty":self.full_data["ny"]//2,"tz":self.full_data["nz"]//2+self.zpos})
		# # self.full_data.set_rotated_clip(self.xform, self.data)
		# # self.full_annotation.set_rotated_clip(self.xform, self.annotation)
		#
		# #self.target.get_full_annotation()["ann_name"] = name_str
		# serialize_name = json.dumps(names, default=lambda a: "[%s,%s]" % (str(type(a)), a.pk))
		# self.target.get_full_annotation()["ann_name"] = serialize_name
		# self.target.get_full_annotation()["ann_num"] = nums
		# self.target.get_full_annotation()["ann_group"] = group_id
		self.target.get_full_annotation().write_image(out_name)
		print("Annotation is really saved to", out_name)
		return

	def update_sets(self):
		self.target.need_new_RGB = 0
		#print("NEED NEW RGB", self.target.need_new_RGB)
		it = QtWidgets.QTreeWidgetItemIterator(self.tree_set)
		self.used_index = [0]
		while it.value():
			item = it.value()
			self.used_index.append(int(item.text(0)))
			if item.text(2) != "-1":
				self.target.display_group = True
			if item.parent():
				if item.parent().text(2) != "-1":
					item.setText(2,item.parent().text(2))
				else:
					item.setText(2,item.parent().text(0))
			else:
				item.setText(2,"-1")
			self.color_tree_item(item)
			it += 1
			# key = int(self.tree_set.topLevelItem(i).text(0))
			# self.tree_set.topLevelItem(i).setFlags(self.itemflags)
			# self.tree_set.topLevelItem(i).setForeground(0,self.colors[key])
			# self.tree_set.topLevelItem(i).setForeground(1,self.colors[key])



			#group_key = int(self.tree_set.item(i,2).text())
			#if group_key != -1:
			# 		self.table_set.item(i,1).setForeground(self.colors[group_key+100])
			# 		self.table_set.item(i,2).setForeground(self.colors[group_key+100])
			# 	else:
			# 		self.table_set.item(i,1).setForeground(self.colors[key])
			# 		self.table_set.item(i,2).setForeground(QtGui.QColor.fromRgb(120,120,120))

		#Set the colors and flags of table set items
		# for i in range(self.table_set.rowCount()):
		# 	key = int(self.table_set.item(i,0).text())
		#
		# 	self.table_set.item(i,0).setFlags(self.indexflags)
		# 	self.table_set.item(i,1).setFlags(self.itemflags)
		# 	self.table_set.item(i,0).setForeground(self.colors[key])
		# 	#self.table_set.item(i,1).setForeground(self.colors[key])
		# 	group_key = int(self.table_set.item(i,2).text())
		# 	self.table_set.item(i,2).setFlags(self.indexflags)
		# 	if group_key != -1:
		# 		self.table_set.item(i,1).setForeground(self.colors[group_key+100])
		# 		self.table_set.item(i,2).setForeground(self.colors[group_key+100])
		# 	else:
		# 		self.table_set.item(i,1).setForeground(self.colors[key])
		# 		self.table_set.item(i,2).setForeground(QtGui.QColor.fromRgb(120,120,120))

				# self.tree_set.itemAt(i,1).setForeground(self.colors[key])
				# self.tree_set.itemAt(i,2).setForeground(QtGui.QColor.fromRgb(120,120,120))
	def color_tree_item(self,item):
		#item.setFlags(self.itemflags)
		key = int(item.text(0))
		item.setForeground(0,self.colors[key])
		if int(item.text(2)) == -1:
			key = int(item.text(0))
				# print(key)
			#item.setForeground(0,self.colors[key])
			item.setForeground(1,self.colors[key])
			item.setForeground(2,QtGui.QColor.fromRgb(120,120,120))
		else:
			key = int(item.text(2))
			item.setForeground(1,self.colors[key])
			item.setForeground(2,self.colors[key])




	def read_header(self,file_in,ret = False):
		#Read the header information from annotation file
		group_ids = []
		try:
			group_ids = file_in["ann_group"]
		except:
			print("No group information. Pass")
			pass
		try:
			keys = file_in["ann_num"]
			#values = file_in["ann_name"].split(",")
			values = json.loads(file_in["ann_name"])
			#print(values)
			if type(keys) != list: keys = [keys]
			item_dict=(dict(zip(keys, values)))
			print(item_dict)

			if ret:
				return item_dict
			else:
				# for key, value in item_dict.items():
				# 	self.add_new_row(key,value)
				if len(group_ids)>0:
					for i in range(len(keys)):
						self.add_new_row(keys[i],values[i],group_ids[i])
					# if len(self.nodes)==0:
					# 	self.make_nodes()
					# print("nodes",[node.get_value() for node in (self.nodes.values())])
					# print("rowCount",self.table_set.rowCount())
					for row in range(self.table_set.rowCount()):
						print(int(self.table_set.item(row,0).text()))
						if int(self.table_set.item(row,2).text()) != -1:
							#print("Row", row)
							group_id = int(self.table_set.item(row,2).text())
							#self.nodes[100+group_id].add_child(self.nodes[int(self.table_set.item(row,0).text())])
							continue
						if int(self.table_set.item(row,0).text()) >=101:
							#print("Group",int(self.table_set.item(row,0).text()),"detected")
							self.target.display_group.append(int(self.table_set.item(row,0).text()))
							self.used_group_index.append(int(self.table_set.item(row,0).text()))

				else:
					for i in range(len(keys)):
						self.add_new_row(keys[i],values[i])
				self.update_sets()
				self.target.force_display_update()
				self.target.updateGL()
		except:
			#print("Trouble reading headers information. Continue...")
			pass


	def reset(self):
		#reset seg tab when set new data
		return




	def add_new_row(self,num,name,group_num=-1):
		next = self.tree_set.topLevelItemCount()
		#print("Next", next)
		item = QtWidgets.QTreeWidgetItem([str(num),name,str(group_num)])
		self.tree_set.insertTopLevelItem(next, item)


	def write_treeset_json(self,json_file=""):
		if json_file == "" or not (json_file.endswith(".json")):
			print("Must specify info json file")
			return
		else:
			pass
		def tree_to_dict(parent):
			childCount = parent.childCount()
			if not childCount:
				return None
			content = {}
			for row in range(childCount):
				child = parent.child(row)

				text = [child.text(0),child.text(1),child.text(2)]

				ser_text =  json.dumps(text, default=lambda a: "[%s,%s]" % (str(type(a)), a.pk))

				#content[child.text(0)] = tree_to_dict(child)
				content[ser_text] = tree_to_dict(child)
			return content

		json_str = tree_to_dict(self.tree_set.invisibleRootItem())
		js=js_open_dict(json_file)
		js['tree_dict'] = json_str

	def read_json_treeset(self, json_file=""):
		if json_file == "" or not (json_file.endswith(".json")):
			print("Must specify info json file")
			return
		else:
			pass
		def fill_item(item, value):
			def new_item(parent, text, val=None):

				ser_text = list(json.loads(text))
				child = QtWidgets.QTreeWidgetItem(ser_text)
				fill_item(child, val)
				parent.addChild(child)
				#child.setExpanded(True)
			if value is None: return
			elif isinstance(value, dict):
				for key, val in sorted(value.items()):
					new_item(item, str(key), val)
			elif isinstance(value, (list, tuple)):
				for val in value:
					text = (str(val) if not isinstance(val, (dict, list, tuple))
							else '[%s]' % type(val).__name__)
					new_item(item, text, val)
			else:
				new_item(item, str(value))
		#js=js_open_dict("./test_json.json")
		js=js_open_dict(json_file)
		json_str = js['tree_dict']
		fill_item(self.tree_set.invisibleRootItem(),json_str)


def main():
	from eman2_gui.emapplication import EMApp
	em_app = EMApp()
	window = EMAnnotate2DWidget(application=em_app)

	if len(sys.argv)==1 :
		window.set_data(test_image(size=(512,512)))
	elif len(sys.argv)==3:
		a=EMData(sys.argv[1])
		b=EMData(sys.argv[2])
		window.set_data(a,b)
	else:
		a=EMData(sys.argv[1])
		window.set_data(a,None)



	em_app.show()
	window.optimally_resize()
	sys.exit(em_app.exec_())





if __name__ == '__main__':
	main()
