#!/usr/bin/env python
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

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

from past.utils import old_div
from builtins import range
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from OpenGL import GL,GLU,GLUT
from OpenGL.GL import *
from .valslider import ValSlider,ValBox,StringBox
from math import *
import EMAN2db
from EMAN2 import *
import EMAN2
import sys
import numpy
import struct
from .emimageutil import ImgHistogram, EMParentWin
from . import emshape
from .emshape import EMShape
from weakref import WeakKeyDictionary
import weakref
from pickle import dumps,loads
from libpyGLUtils2 import *

from .emglobjects import EMOpenGLFlagsAndTools
from .emapplication import get_application, EMGLWidget
from .emimageutil import EMMetaDataTable

from .emanimationutil import SingleValueIncrementAnimation, LineAnimation

import platform

from .emglobjects import EMOpenGLFlagsAndTools

class EMImage2DWidget(EMGLWidget):
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

	allim=WeakKeyDictionary()

	def __init__(self, image=None, application=get_application(),winid=None, parent=None):

		self.inspector = None # this should be a qt widget, otherwise referred to as an inspector in eman

		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True)
		#fmt.setSampleBuffers(True)
		fmt.setDepth(1)
		EMGLWidget.__init__(self,parent)
		self.setFormat(fmt)
		self.setFocusPolicy(Qt.StrongFocus)
		self.setMouseTracking(True)
		self.initimageflag = True

		self.fftorigincenter = E2getappval("emimage2d","origincenter")
		if self.fftorigincenter == None : self.fftorigincenter=False

		#sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Policy(7),QtGui.QSizePolicy.Policy(7))
		#sizePolicy.setHorizontalStretch(7)
		#sizePolicy.setVerticalStretch(7)
		#sizePolicy.setHeightForWidth(False)
		#self.setSizePolicy(sizePolicy)

#		self.data = image 	   # EMData object to display
		self.data = None		# set to image below
		self.file_name = ""# stores the filename of the image, if None then member functions should be smart enough to handle it
		self.enable_clip = False
		EMImage2DWidget.allim[self] = 0

		self.init_gl_flag = True
		self.oldsize=(-1,-1)
		self.scale=1.0				# Scale factor for display
		self.origin=(0,0)			# Current display origin
		self.invert=0				# invert image on display
		self.histogram=0            # histogram equalization
		self.gamma=1.0				# gamma for display (impact on inverted contrast ?
		self.minden=0
		self.maxden=1.0
		self.curmin=0.0
		self.curmax=0.0
		self.maxden=1.0
		self.fgamma = 1.0
		self.fminden=0
		self.fmaxden=1.0
		self.fcurmin=0.0
		self.fcurmax=0.0
		self.disp_proc=[]			# a list/set of Processor objects to apply before rendering
		self.display_fft = None		# a cached version of the FFT
		self.fft=None				# The FFT of the current target if currently displayed
		self.rmousedrag=None		# coordinates during a right-drag operation
		self.mouse_mode_dict = {0:"emit", 1:"emit", 2:"emit", 3:"probe", 4:"measure", 5:"draw", 6:"emit"}
		self.mouse_mode = 0         # current mouse mode as selected by the inspector
		self.curfft=0				# current FFT mode (when starting with real images only)
		self.mag = 1.1				# magnification factor
		self.invmag = old_div(1.0,self.mag)	# inverse magnification factor

		self.shapes={}				# dictionary of shapes to draw, see add_shapes
		self.shapechange=1			# Set to 1 when shapes need to be redrawn
		self.active=(None,0,0,0)	# The active shape and a hilight color (n,r,g,b)

		self.extras = []			# an empty set of extras - other images that can be rendered over this one

		self.startorigin = None
		self.endorigin = None
		self.isanimated = False
		self.time = 1
		self.timeinc = 0.125
		self.key_mvt_animation = None

		self.init_size = True		# A flag used to set the initial origin offset

		self.shapelist = 0			# a display list identify

		self.glflags = EMOpenGLFlagsAndTools() 	# supplies power of two texturing flags

		self.tex_name = 0			# an OpenGL texture handle

		self.window_width = None # Used for intelligently managing resize events
		self.window_height = None # Used for intelligently managing resize events

		self.otherdata = None
		self.otherdatascale = -1
		self.otherdatablend = False
		self.otherdataupdate = False # Used for forcing a
		self.otherdatadl = 0
		self.other_tex_name = None
		self.init_size_flag = True
		self.frozen = False
		self.isexcluded = False
		self.parent_geometry = None

		self.eraser_shape = None # a single circle shape used 0for erasing in e2boxer

		self.list_data = None 			# this can be used for viewing lists of data
		self.list_fft_data = None		# this is used for doing the ffts of list data
		self.list_idx = 0	# and idx to the list_data

		self.use_display_list = True # whether or not a display list should be used to render the image pixelsw - if on, this will save on time if the view of the image is unchanged, which can quite often be the case
		self.main_display_list = -1	# if using display lists, the stores the display list
		self.display_states = [] # if using display lists, this stores the states that are checked, and if different, will cause regeneration of the display list
		self.hist = []

		self.display_shapes = True # A flag that can be used to turn of the display of shapes - useful to e2boxer

		self.wheel_navigate = False # useful on Mac laptops
		self.circle_dl = None # used for a circle list, for displaying circled particles, for example

		self.setAcceptDrops(True) #TODO: figure out the purpose of this (moved) line of code
		self.setWindowIcon(QtGui.QIcon(get_image_directory() +"single_image.png")) #TODO: figure out why this icon doesn't work

		if image : self.set_data(image)
#		else:self.__load_display_settings_from_db()

#	def __del__(self):
#		#self.clear_gl_memory() # this is intentionally commented out, it makes sense to clear the memory but not here
#		self.qt_parent.deleteLater()

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
		side = min(width, height)
		GL.glViewport(0,0,width,height)

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

		data = self.data
		if data == None: data = self.fft

		if data != None and  data.get_xsize()<640 and data.get_ysize()<640:
			try : return (data.get_xsize()+12,data.get_ysize()+12)
			except : return (640,640)
		else:
			return (640,640)

	def sizeHint(self):
#		print self.get_parent_suggested_size()
		if self.data==None : return QtCore.QSize(512,512)
		return QtCore.QSize(*self.get_parent_suggested_size())

	def set_disp_proc(self,procs):
		self.disp_proc=procs
		self.force_display_update()
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
		self.del_shapes()

	def get_minden(self): return self.minden
	def get_maxden(self): return self.maxden
	def get_gamma(self): return self.gamma
	def get_shapes(self): return self.shapes

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
		if self.curfft == 0:
			self.curmin=x0
			self.curmax=x1
		else:
			self.fcurmin=x0
			self.fcurnax=x1
		self.force_display_update()
		self.updateGL()

	def set_density_min(self,val):
		if self.curfft == 0:
			self.curmin=val
		else:
			self.curmax=val
		self.force_display_update()
		self.updateGL()

	def set_density_max(self,val):
		if self.curfft == 0:
			self.curmax=val
		else:
			self.fcurmax=val
		self.force_display_update()
		self.updateGL()

	def set_gamma(self,val):
		if self.curfft == 0:
			self.gamma=val
		else:
			self.fgamma=val
		self.force_display_update()
		self.updateGL()

	def set_file_name(self,file_name,load_cache_settings=True):
		self.file_name = file_name

		#if load_cache_settings:
			#self.__load_display_settings_from_db()

	def get_file_name(self):
		return self.file_name

	def set_other_data(self,data,scale,blend=False):
		self.otherdata = data
		self.otherdatascale = scale
		self.otherdatablend = blend
		self.otherdataupdate=True

	def get_data_dims(self):
		data = None

		if self.data != None: data = self.data
		elif self.fft != None: data = self.fft
		else: return [0,0,0]

		return [data.get_xsize(),data.get_ysize(),data.get_zsize()]

#	def updateGL(self):
#		if self.gl_widget != None:
#			self.gl_widget.updateGL()

	def set_frozen(self,frozen):
		self.frozen = frozen

	def set_excluded(self,isexcluded):
		wasexcluded = self.isexcluded

		self.isexcluded = isexcluded

		if wasexcluded or self.isexcluded: return True
		else: return False

	def get_data(self):
		return self.data

	def set_data(self,incoming_data,file_name="",retain_current_settings=True, keepcontrast=False):
		"""You may pass a single 2D image or a list of images"""
		from .emimagemx import EMDataListCache,EMLightWeightParticleCache
		#if self.data != None and self.file_name != "":
			#self.__write_display_settings_to_db()

		self.set_file_name(file_name,load_cache_settings=False)
		if self.file_name != "": self.setWindowTitle(remove_directories_from_name(self.file_name))

		data = incoming_data
		if data == None:
			self.data = None
			return

		if self.data==None :
			needresize=True
		else:
			needresize=False

		fourier = False

		# it's a 3D image
		if not isinstance(data,list) and not isinstance(data,EMDataListCache) and not isinstance(data,EMLightWeightParticleCache) and data.get_zsize() != 1:
			data = []
			for z in range(incoming_data.get_zsize()):
				image = incoming_data.get_clip(Region(0,0,z,incoming_data.get_xsize(),incoming_data.get_ysize(),1))
				data.append(image)


		if isinstance(data,list) or isinstance(data,EMDataListCache) or isinstance(data,EMLightWeightParticleCache):
			if self.list_data == None and self.list_idx > len(data): self.list_idx = old_div(len(data),2) #otherwise we use the list idx from the previous list data, as in when being used from the emselector
			d = data[0]
			if d.is_complex():
				self.list_data = []
				self.list_fft_data = data
				for i in range(len(data)):self.list_data.append(None)
				self.curfft = 2
				self.__set_display_image(self.curfft)
				fourier = True
			elif self.curfft in [1,2,3]:
				self.list_data = data
				self.data = self.list_data[self.list_idx]
				self.list_fft_data = [d.do_fft() for d in data]
				if self.fftorigincenter :
					for im in self.list_fft_data : im.process_inplace("xform.phaseorigin.tocorner")
				self.__set_display_image(self.curfft)
			else:
				self.list_data = data
				self.data = self.list_data[self.list_idx]
				self.list_fft_data = []
				for i in range(len(data)):self.list_fft_data.append(None)

			self.get_inspector().enable_image_range(1,len(data),self.list_idx+1)
		else:
			self.list_data = None
			self.list_fft_data = None
			if data.is_complex() or self.curfft in [1,2,3]:
				self.display_fft = None
				self.data = None
				fourier = True
				if data.is_complex():
					self.fft = data.copy()# have to make copies here because we alter it!
					if self.curfft == 0:
						self.curfft = 2 # switch to displaying amplitude automatically
						inspector = self.get_inspector()
						inspector.set_fft_amp_pressed()
				else:
					self.fft = data.do_fft()
					if self.fftorigincenter : self.fft.process_inplace("xform.phaseorigin.tocorner")
				self.fft.set_value_at(0,0,0,0) # get rid of the DC component
				self.fft.set_value_at(1,0,0,0) # this should already by 0... ?

				self.__set_display_image(self.curfft)
				fourier = True
			else:
				self.data = data
				self.display_fft = None
				self.fft = None

		self.image_change_count = 0
		
		if not keepcontrast:
			self.auto_contrast(inspector_update=False,display_update=False)

		#if not retain_current_settings:
			#self.__load_display_settings_from_db(inspector_update=False,display_update=False)

		try:
			if needresize:
				x=self.data["nx"]
				y=self.data["ny"]
				xys=QtWidgets.QApplication.desktop().availableGeometry()
				mx=old_div(xys.width()*2,3)
				my=old_div(xys.height()*2,3)

				self.resize(min(x,mx),min(y,my))
		except: pass

		self.inspector_update(use_fourier=fourier)
		if self.curfft in [1,2,3] and self.data!=None and self.data.is_complex() : self.redo_fft()
		self.force_display_update()
		self.updateGL()

		if not isinstance(data,list) and not isinstance(data,EMDataListCache) and not isinstance(data,EMLightWeightParticleCache):
			self.get_inspector().disable_image_range()

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

		#print "scale is ", self.scale
		#except: pass

	def auto_contrast(self,bool=False,inspector_update=True,display_update=True):
		auto_contrast = E2getappval("display2d","autocontrast",True)
		
		if self.curfft == 0:
			if self.data == None: return
			# histogram is impacted by downsampling, so we need to compensate
			if self.scale<=0.5 :
				down=self.data.process("math.meanshrink",{"n":int(floor(old_div(1.0,self.scale)))})
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
		else:
			if self.display_fft == None: return

			mean=self.display_fft.get_attr("mean")
			sigma=self.display_fft.get_attr("sigma")
			m0=self.display_fft.get_attr("minimum")
			m1=self.display_fft.get_attr("maximum")

			self.fminden=0
			self.fmaxden=min(m1,mean+20.0*sigma)
			self.fcurmin = 0
			if auto_contrast:
				self.fcurmax = min(m1,mean+4.0*sigma)
			else:
				self.fcurmax = m1

			self.force_display_update()

			if inspector_update: self.inspector_update(use_fourier=True)
			if display_update:
				self.force_display_update()
				self.updateGL()

	def __load_display_settings_from_db(self,inspector_update=True,display_update=True):
		if self.file_name == "": return # there is no file name, we have no means to stores information
		try:
			global HOMEDB
			HOMEDB=EMAN2db.EMAN2DB.open_db()
			HOMEDB.open_dict("image_2d_display_settings")
		except:
			# something wrong with the HOMEDB?
			return

		db = HOMEDB.image_2d_display_settings

		data = db[self.file_name]
		if data == None: return
		try:
			self.minden = data["min"]
			self.maxden = data["max"]
			self.minden = data["min"]
			self.fminden = data["fourier_min"]
			self.fmaxden = data["fourier_max"]
			self.curmin = data["curmin"]
			self.curmax = data["curmax"]
			self.fcurmin = data["fcurmin"]
			self.fcurmax = data["fcurmax"]
			self.fgamma = data["fourier_gamma"]
			self.gamma = data["gamma"]
			self.scale = data["scale"]
			self.origin = data["origin"]
		except:
			# perhaps in inconsistency, something wasn't set. This should not happen in general
			pass

		try:
			self.parent_geometry = data["parent_geometry"]
			if self.qt_parent != None:
				try:
					self.qt_parent.restoreGeometry(self.parent_geometry)
				except: pass
		except:pass

		if inspector_update: self.inspector_update()
		if display_update: self.force_display_update()

	def __write_display_settings_to_db(self):
		'''
		writes the min,max, brightness, contrast and gamma values associated with
		the current image to the homedb. The full path of the image is used
		'''

		if self.file_name == None: return # there is no file name, we have no means to stores information

		try:
			DB = HOMEDB
			DB.open_dict("image_2d_display_settings")
		except:
			# Databasing is not supported, in which case we do nothing
			return

		data = {}
		data["min"] = self.minden
		data["max"] = self.maxden
		data["curmin"] = self.curmin
		data["curmax"] = self.curmax
		data["fourier_min"] = self.fminden
		data["fourier_max"] = self.fmaxden
		data["fcurmin"] = self.fcurmin
		data["fcurmax"] = self.fcurmax
		data["fourier_gamma"] = self.fgamma
		data["gamma"] = self.gamma
		data["origin"] = self.origin
		data["scale"] = self.scale

		try:
			data["parent_geometry"] = self.qt_parent.saveGeometry()
		except: pass

		db = DB.image_2d_display_settings
		db[self.file_name] = data

	def set_origin(self,x,y,quiet=False):
		"""Set the display origin within the image"""
		if self.origin==(x,y) : return
		self.origin=(x,y)
		if not quiet : self.origin_update.emit((x,y))
		self.updateGL()

	def get_origin(self) : return self.origin

	def scroll_to(self,x=None,y=None):
		"""center the point on the screen"""
		if x==None:
			if y==None: return
			self.set_origin(self.origin[0],y*self.scale-old_div(self.height(),2))
		elif y==None:
			self.set_origin(x*self.scale-old_div(self.width(),2),self.origin[1])
		else: self.set_origin(x*self.scale-old_div(self.width(),2),y*self.scale-old_div(self.height(),2))

	def set_shapes(self,shapes):
		self.shapes = shapes
		self.shapechange=1

	def update_shapes(self,shapes):
		self.shapes.update(shapes)
		self.shapechange=1

	def register_scroll_motion(self,x,y,z=0):
		if self.list_data!=None:
			self.image_range_changed(z+1)
			#self.setup_shapes()
		animation = LineAnimation(self,self.origin,(x*self.scale-old_div(self.width(),2),y*self.scale-old_div(self.height(),2)))
		self.qt_parent.register_animatable(animation)
		return True

	def set_scale(self,newscale,quiet=False):
		"""Adjusts the scale of the display. Tries to maintain the center of the image at the center"""
		if self.scale==newscale: return
		try:
			self.origin=(old_div(newscale,self.scale)*(old_div(self.width(),2.0)+self.origin[0])-old_div(self.width(),2.0),old_div(newscale,self.scale)*(old_div(self.height(),2.0)+self.origin[1])-old_div(self.height(),2.0))
			self.scale=newscale
			if not quiet : self.signal_set_scale.emit(newscale)
			self.updateGL()
		except: pass

	def set_invert(self,val):
		if val: self.invert=1
		else : self.invert=0
		self.updateGL()

	def set_histogram(self,mode):
		self.histogram=mode
		self.updateGL()

	def set_FFT(self,val):
		if self.data != None and self.data.is_complex():
			print(" I am returning")
			return

		self.curfft=val
		fourier = self.__set_display_image(self.curfft)

		self.inspector_update(use_fourier=fourier)

		self.force_display_update()
		self.updateGL()

	def redo_fft(self):
		if self.list_data == None:
			self.fft = None
		else:
			self.list_fft_data[self.list_idx] = None

		if self.curfft > 0: self.__set_display_image(self.curfft)

	def force_fft_redo(self):
		'''
		Called from image_update in emimage.py, this is useful if you
		have two display windows displaying the same image, but one shows
		the fft and the other the real space image. If you draw on the real
		space image you want the FFT to update in real time.
		It only redoes the FFT if the FFT is currently being displayed. Otherwise
		it sets a flag for the FFT to be recalculated next time it is displayed.
		'''
		if self.curfft > 0:
			self.fft = None
			self.display_fft = None
			self.__set_display_image(self.curfft)
		else:
			self.fft = None
			self.display_fft = None

	def __set_display_image(self,val):
		if self.list_data == None:
			if val > 0 :
				try:
					if self.fft == None:
						self.fft = self.data.do_fft()
						self.fft.set_value_at(0,0,0,0)
						self.fft.set_value_at(1,0,0,0)
						if self.fftorigincenter : self.fft.process_inplace("xform.phaseorigin.tocorner")
					if val==1 :
#						self.display_fft = self.fft.process("xform.phaseorigin.tocorner")
						self.display_fft = self.fft.copy()
						return True
					elif val==2 :
						self.display_fft = self.fft.process("xform.fourierorigin.tocenter")
						self.display_fft = self.display_fft.get_fft_amplitude()
						return True
					elif val==3 :
						self.display_fft = self.fft.process("xform.fourierorigin.tocenter")
						self.display_fft = self.display_fft.get_fft_phase()
						return True
				except:
					self.display_fft=None
			elif val == 0:
				if self.data == None:
					self.data = self.fft.do_ift()
				return False
			else:						# val is negative!?!?
				self.display_fft=None

			return False
		else:
			if self.list_idx>=len(self.list_data): self.list_idx=len(self.list_data)-1
			if val > 0 :
				try:
					if len(self.list_fft_data)!=len(self.list_data) : self.list_fft_data=[None for i in range(len(self.list_data))]
					if self.list_fft_data[self.list_idx] == None:
						self.list_fft_data[self.list_idx] = self.list_data[self.list_idx].do_fft()
						if self.fftorigincenter :
							self.list_fft_data[self.list_idx].process_inplace("xform.phaseorigin.tocorner")
					fft = self.list_fft_data[self.list_idx]
					if val==1 :
#						self.display_fft = fft.process("xform.phaseorigin.tocorner")
						self.display_fft = fft.copy()
						return True
					elif val==2 :
						self.display_fft = fft.process("xform.fourierorigin.tocenter")
						self.display_fft = self.display_fft.get_fft_amplitude()
						return True
					elif val==3 :
						self.display_fft = fft.process("xform.fourierorigin.tocenter")
						self.display_fft = self.display_fft.get_fft_phase()
						return True
				except:
					self.display_fft=None
			elif val == 0:
				if self.list_data[self.list_idx] == None:
					self.list_data[self.list_idx] = self.list_fft_data[self.list_idx].do_ift()

				self.data = self.list_data[self.list_idx]
				return False

			else:
				self.display_fft=None

			return False

	def force_display_update(self):
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
		display_states.append(self.histogram)
		display_states.append(self.minden)
		display_states.append(self.maxden)
		display_states.append(self.gamma)
		display_states.append(self.curfft)
		display_states.append(self.curmin)
		display_states.append(self.fcurmin)
		display_states.append(self.fcurmax)
		display_states.append(self.curmax)
		if len(self.display_states) == 0:
			self.display_states = display_states
			return True
		else:
			for i in range(len(display_states)):

				if display_states[i] != self.display_states[i]:
					self.display_states = display_states
					return True

		return False

	def render_bitmap (self) :
		"""This will render the current bitmap into a string and
			return a tuple with 1 or 3 (grey vs RGB), width, height,
			and the raw data string, with optional histogram data appended
			(1024 bytes)."""

		if self.invert :
			pixden = (255, 0)
		else :
			pixden = (0, 255)

		use_fft = (self.curfft in (1,2,3))

		have_textures = (not self.glflags.npt_textures_unsupported ( ))

		if self.curfft == 1 :
			# doing a 24-bit color FFT:

			if not self.display_fft.is_complex ( ) :
				print("Error, the FFT is not complex; internal error")
				return None

			value_size = 3  # 3 8-bit bytes
			flags = 3       # 2 (histogram) + 1 (an R,G,B flag)
		else :
			value_size = 1  # 1 8-bit byte

			if self.curfft in (2,3) :
				# doing a grey scale FFT (of amplitude or phase):

#				if have_textures :
					flags = 2  # histogram
#				else :
#					flags = 6  # 2 (histogram) + 4 (invert y)
			else :
#				if have_textures :
					# doing a histogram:

					if self.histogram == 1 :
						# flat histogram:

						flags = 34  # 2 (histogram) + 32 (not ordinary histogram)
					elif self.histogram == 2 :
						# Gaussian histogram:

						flags = 98  # 2 + 32 + 64 (Gaussian histogram)
					else :  ## self.histogram == 0
						# ordinary histogram:

						flags = 2   # 2 (histogram)
#				else :
#					flags = 6  # 2 (histogram) + 4 (invert y)

		if use_fft :
			# doing an FFT:

			if value_size != 3 :
				values = self.display_fft
			else:
				values = self.display_fft

			min_val = self.fcurmin
			max_val = self.fcurmax
			gam_val = self.fgamma
		else :
			# using real data

			values  = self.data
			min_val = self.curmin
			max_val = self.curmax
			gam_val = self.gamma

		if len(self.disp_proc)>0 :
			tmp=values.copy()
			for p in self.disp_proc: p.process_inplace(tmp)
			values=tmp

		wid =  old_div((self.width() * value_size - 1) , 4) * 4+ 4
		wdt =  self.width()
		hgt =  self.height()

		x0  = 1 + int(old_div(self.origin[0], self.scale))
		y0  = 1 + int(old_div(self.origin[1], self.scale))

#		print "--------------------------------------------------------------"
#		print "invert, curfft, histogram:", \
#				 self.invert, self.curfft, self.histogram
#		print "size, x0, y0, wdt, hgt, wid, flags:", \
#				 value_size, x0, y0, wdt, hgt, wid, flags
#		print "scl, pix0, pix1, min, max, gam:", \
#				 self.scale, pixden[0], pixden[1], min_val, max_val, gam_val

		if value_size == 3 :
			# get color FFT data:

			return_data = (value_size, wid, hgt,
								self.display_fft.render_ap24 (
								x0, y0, wdt, hgt, wid,
								self.scale, pixden[0], pixden[1],
								min_val, max_val, gam_val, flags))
		else :
			# get grey scale data:

			return_data = (value_size, wid, hgt,
								GLUtil.render_amp8 (values,
								x0, y0, wdt, hgt, wid,
								self.scale, pixden[0], pixden[1],
								min_val, max_val, gam_val, flags))

		return return_data

	def render_bitmap_old(self):   # no longer used - use new render_bitmap
		"""This will render the current bitmap into a string and return a tuple with 1 or 3 (grey vs RGB), width, height, and the raw data string"""
		if not self.invert : pixden=(0,255)
		else: pixden=(255,0)


		if self.curfft==1 :
			if self.display_fft.is_complex() == False:
				print("error, the fft is not complex, internal error")
				return
			a=(3,old_div((self.width()*3-1),4)*4+4,self.height(),self.display_fft.render_ap24(1+int(old_div(self.origin[0],self.scale)),1+int(old_div(self.origin[1],self.scale)),self.width(),self.height(),old_div((self.width()*3-1),4)*4+4,self.scale,pixden[0],pixden[1],self.fcurmin,self.fcurmax,self.fgamma,3))
		elif self.curfft in (2,3) :
#			if not self.glflags.npt_textures_unsupported():
				a=(1,old_div((self.width()-1),4)*4+4,self.height(),GLUtil.render_amp8(self.display_fft, 1+int(old_div(self.origin[0],self.scale)),1+int(old_div(self.origin[1],self.scale)),self.width(),self.height(),old_div((self.width()-1),4)*4+4,self.scale,pixden[0],pixden[1],self.fcurmin,self.fcurmax,self.fgamma,2))
#			else :
#				a=(1,(self.width()-1)/4*4+4,self.height(),GLUtil.render_amp8(self.display_fft, 1+int(self.origin[0]/self.scale),1+int(self.origin[1]/self.scale),self.width(),self.height(),(self.width()-1)/4*4+4,self.scale,pixden[0],pixden[1],self.fcurmin,self.fcurmax,self.fgamma,6))
		else :
#			if not self.glflags.npt_textures_unsupported():
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
		if self.curfft and not self.display_fft : return
		if not self.curfft and not self.data : return

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

		width = old_div(self.width(),2.0)
		height = old_div(self.height(),2.0)

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
			(bpp,w,h,a)=self.render_bitmap()
			if bpp==3 : gl_render_type = GL_RGB
			elif bpp==4 : gl_render_type = GL_RGBA
			else : gl_render_type = GL_LUMINANCE

			if not self.glflags.npt_textures_unsupported():

				self.hist=struct.unpack('256i',a[-1024:])

				if self.tex_name != 0: glDeleteTextures(self.tex_name)
				self.tex_name = glGenTextures(1)
				if ( self.tex_name <= 0 ):
					raise("failed to generate texture name")

				#if self.otherdatablend and self.otherdata != None:

					#glBlendFunc(GL_DST_ALPHA,GL_ONE_MINUS_DST_ALPHA)

				GL.glBindTexture(GL.GL_TEXTURE_2D,self.tex_name)
				glPixelStorei(GL_UNPACK_ALIGNMENT,4)
				GL.glTexImage2D(GL.GL_TEXTURE_2D,0,gl_render_type,old_div(w,bpp),h,0,gl_render_type, GL.GL_UNSIGNED_BYTE, a)

				glNewList(self.main_display_list,GL_COMPILE)
				GL.glBindTexture(GL.GL_TEXTURE_2D,self.tex_name)
				glPushMatrix()
				glTranslatef(width,height,0)
				self.__draw_texture(self.tex_name,-width,-height,width,height)
				glPopMatrix()

				#

				#if self.otherdatablend and self.otherdata != None:
				#	GL.glDisable( GL.GL_BLEND);
				#	if (depth_testing_was_on):	GL.glEnable(GL.GL_DEPTH_TEST)

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


		if self.otherdata != None and isinstance(self.otherdata,EMData) and not self.glflags.npt_textures_unsupported():
			if self.otherdataupdate or self.otherdatadl == 0:


				if self.other_tex_name != 0: glDeleteTextures(self.other_tex_name)

				scale = self.scale*self.otherdatascale
				b=GLUtil.render_amp8(self.otherdata, int(old_div(self.origin[0],scale)),int(old_div(self.origin[1],scale)),self.width(),self.height(),old_div((self.width()-1),4)*4+4,scale,pixden[0],pixden[1],0,1,1,2)
				gl_render_type = GL_LUMINANCE

				if self.other_tex_name != 0: GL.glDeleteTextures(self.other_tex_name)
				self.other_tex_name = GL.glGenTextures(1)
				if ( self.other_tex_name <= 0 ):
					raise("failed to generate texture name")

				glBindTexture(GL.GL_TEXTURE_2D,self.other_tex_name)
				glPixelStorei(GL_UNPACK_ALIGNMENT,4)
				glTexImage2D(GL.GL_TEXTURE_2D,0,gl_render_type,self.width(),self.height(),0,gl_render_type, GL.GL_UNSIGNED_BYTE, b)


				if self.otherdatadl != 0: glDeleteLists(self.otherdatadl,1)

				self.otherdatadl = glGenLists(1)

				glNewList(self.otherdatadl,GL_COMPILE)
				glBindTexture(GL.GL_TEXTURE_2D,self.other_tex_name)
				glDisable(GL_DEPTH_TEST)
				GL.glEnable(GL.GL_BLEND);
						#depth_testing_was_on = GL.glIsEnabled(GL.GL_DEPTH_TEST);
				try:
					GL.glBlendEquation(GL.GL_FUNC_REVERSE_SUBTRACT);
				except: pass
	#					GL.glBlendFunc(GL.GL_SRC_ALPHA,GL.GL_ONE_MINUS_SRC_ALPHA)
				GL.glBlendFunc(GL.GL_ONE,GL.GL_ONE);
				glPushMatrix()
				glTranslatef(width,height,0)
				self.__draw_texture(self.other_tex_name,-width,-height,width,height)
				glPopMatrix()

				GL.glDisable(GL_BLEND)
				glEnable(GL_DEPTH_TEST)
				glEndList()

			glCallList(self.otherdatadl)



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
		if data == None: data = self.fft
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

		# make our own cirle rather than use gluDisk or somesuch
		emshape.EMShape.font_renderer=self.font_renderer		# Important !  Each window has to have its own font_renderer. Only one context active at a time, so this is ok.
		glNewList(self.shapelist,GL_COMPILE)

		isanimated = False
		alpha = 1.0
		if len(self.shapes) > 0:

			for k in list(self.shapes.keys()):
				shape = self.shapes[k]
				if not isinstance(shape,EMShape) : continue
				glLineWidth(2)
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

		glPointSize(2)
		for k,s in list(self.shapes.items()):
			### handle boxes for 3D images
			if s.shape[0] == "ellipse":
				mxlen=11
			else:
				mxlen=10

			if self.list_data!=None and len(s.shape)==mxlen:
				z_idx=s[mxlen-1]
				if z_idx!=self.list_idx:
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
					GL.glLineWidth(p[7])
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
					GL.glLineWidth(p[9])
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
			self.eraser_shape.draw()

		glEndList()

	def inspector_update(self,use_fourier=False):
		if self.inspector:
			if not use_fourier:
				self.inspector.set_limits(self.minden,self.maxden,self.curmin,self.curmax)
				self.inspector.set_gamma(self.gamma)
			else:
				self.inspector.set_limits(self.fminden,self.fmaxden,self.fcurmin,self.fcurmax)

				self.inspector.set_gamma(self.fgamma)

			self.inspector.set_scale(self.scale)
			self.inspector.update_brightness_contrast()

	def get_inspector(self):
		if not self.inspector:
			self.inspector=EMImageInspector2D(self)
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
		origin_x = self.scale*(int(old_div(self.origin[0],self.scale))+0.5)
		origin_y = self.scale*(int(old_div(self.origin[1],self.scale))+0.5)

		try: img_coords = ( old_div((v0+origin_x),self.scale), old_div((self.height()-(v1-origin_y)),self.scale) )
		except:	img_coords = (old_div((v0[0]+origin_x),self.scale),old_div((self.height()-(v0[1]-origin_y)),self.scale))

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
		self.__write_display_settings_to_db()
		EMGLWidget.closeEvent(self,event)

	def dragEnterEvent(self,event):
#		f=event.mimeData().formats()
#		for i in f:
#			print str(i)

		if event.provides("application/x-eman"):
			event.setDropAction(Qt.CopyAction)
			event.accept()

	def dropEvent(self,event):
		if EMAN2.GUIbeingdragged:
			self.set_data(EMAN2.GUIbeingdragged)
			EMAN2.GUIbeingdragged=None
		elif event.provides("application/x-eman"):
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
			self.rmousedrag=(event.x(),event.y() )
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

	def mouseMoveEvent(self, event):
		lc=self.scr_to_img(event.x(),event.y())
		if self.rmousedrag:
			self.set_origin(self.origin[0]+self.rmousedrag[0]-event.x(),self.origin[1]-self.rmousedrag[1]+event.y())
			self.rmousedrag=(event.x(),event.y())
#			self.emit(QtCore.SIGNAL("origin_update"),self.origin)
			#try: self.updateGL()
			#except: pass
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

						# displays the pixel value at the current endpoint
						if inspector.target().curfft :
							fft=inspector.target().fft
							if fft==None :
								fft=inspector.target().list_fft_data[inspector.target().list_idx]
							xs=fft.get_xsize()
							ys=fft.get_ysize()
							x,y=int(lc[0])+1,int(lc[1])
							if x<old_div(xs,2) : x=old_div(xs,2)-x
							else : x-=old_div(xs,2)
							if y<old_div(ys,2) : y+=old_div(ys,2)
							else: y-=old_div(ys,2)
							val=fft[x,y]
							inspector.mtshowval.setText("Value: %1.4g + %1.4g i  @(%d,%d)"%(val.real,val.imag,x,y))
							inspector.mtshowval2.setText("       (%1.4g, %1.4g)"%(abs(val),atan2(val.imag,val.real)*57.295779513))
						else :
							try: inspector.mtshowval.setText("Value: %1.4g"%inspector.target().data[int(lc[0]),int(lc[1])])
							except:
								idx=inspector.target().list_idx
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

	def mouseReleaseEvent(self, event):
		get_application().setOverrideCursor(Qt.ArrowCursor)
		lc=self.scr_to_img(event.x(),event.y())
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
					self.redo_fft()
					self.force_display_update()
					self.updateGL()

	def wheelEvent(self, event):
		if not self.wheel_navigate:
			if event.orientation() & Qt.Vertical:
				if self.mouse_mode==0 and event.modifiers()&Qt.ShiftModifier:
					self.mousewheel.emit(event)
					return
				if event.delta() > 0:
					self.set_scale( self.scale * self.mag )
				elif event.delta() < 0:
					self.set_scale(self.scale * self.invmag )
				# The self.scale variable is updated now, so just update with that
				if self.inspector: self.inspector.set_scale(self.scale)
		else:
			move_fac = old_div(1.0,20.0)
			delta = old_div(event.delta(),120.0)

#			print self.origin, self.data.get_xsize(),self.data.get_ysize(),self.scale,self.width(),self.height()

#			print self.origin
			if event.orientation() & Qt.Vertical:
				visible_vertical_pixels = old_div(self.height(),sqrt(self.scale))
				shift_per_delta = move_fac*visible_vertical_pixels
#				print "there are this many visible vertical pixels",visible_vertical_pixels, "deltas", delta, "shift per delta",shift_per_delta
#				print "shifting vertical",event.delta(),shift_per_delta
				self.origin=(self.origin[0],self.origin[1]-delta*shift_per_delta)
			elif event.orientation() & Qt.Horizontal:
				visible_horizontal_pixels = old_div(self.width(),sqrt(self.scale))
				shift_per_delta = move_fac*visible_horizontal_pixels
#				print "shifting horizontal",event.delta(),shift_per_delta
#	   	   	   	print "there are this many visible horizontal pixels",visible_horizontal_pixels, "deltas", delta, "shift per delta",shift_per_delta
				self.origin=(self.origin[0]+delta*shift_per_delta,self.origin[1])
			try: self.updateGL()
			except: pass
#			print "exit",self.origin


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
			if self.list_data != None:
				self.increment_list_data(1)
				self.updateGL()
#			else:
#				self.__key_mvt_animation(0,self.height()*.1)
			self.signal_increment_list_data.emit(1)

		elif event.key() == Qt.Key_Down:
			if self.list_data != None:
				self.increment_list_data(-1)
				self.updateGL()
#			else:
#				self.__key_mvt_animation(0,-self.height()*.1)
			self.signal_increment_list_data.emit(-1)

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

		else:
			self.keypress.emit(event)



	def increment_list_data(self,delta):
		'''
		delta is negative or positive, indicating forward and backwards movement
		'''
		if self.list_data != None:
			if delta > 0:
				if (self.list_idx < (len(self.list_data)-1)):
					self.list_idx += 1
					self.get_inspector().set_image_idx(self.list_idx+1)
					self.__set_display_image(self.curfft)
					self.setup_shapes()
					self.force_display_update()

			elif delta < 0:
				if (self.list_idx > 0):
					self.list_idx -= 1
					self.get_inspector().set_image_idx(self.list_idx+1)
					self.__set_display_image(self.curfft)
					self.setup_shapes()
					self.force_display_update()


	def image_range_changed(self,val):
		l_val = int(val-1)

		if l_val == self.list_idx: return
		else:
			self.list_idx = l_val
			self.__set_display_image(self.curfft)
			self.force_display_update()
			self.updateGL()
			self.setup_shapes()

	def leaveEvent(self,event):
		get_application().setOverrideCursor(Qt.ArrowCursor)
		if self.rmousedrag:
			self.rmousedrag=None

	def __draw_hud(self):
		if self.list_data == None : return

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
		glMaterial(GL_FRONT,GL_SPECULAR,(1.0	, 0.5, 0.2,1.0))
		glMaterial(GL_FRONT,GL_SHININESS,20.0)
		enable_depth = glIsEnabled(GL_DEPTH_TEST)
		glDisable(GL_DEPTH_TEST)
		glColor(1.0,1.0,1.0)
#		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE)
		glNormal(0,0,1)
		glEnable(GL_TEXTURE_2D)
		n = len(self.list_data)
		string = str(self.list_idx+1) + ' / ' + str(n)
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


class EMImageInspector2D(QtGui.QWidget):
	def __init__(self,target) :
		QtGui.QWidget.__init__(self,None)
		self.target=weakref.ref(target)

		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(2)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")

		# This is the tab-bar for mouse mode selection
		self.mmtab = QtGui.QTabWidget()

		# App tab
		self.apptab = QtGui.QWidget()
		self.apptablab = QtGui.QLabel("Application specific mouse functions",self.apptab)
		self.mmtab.addTab(self.apptab,"App")

		# Save tab
		self.savetab = QtGui.QWidget()
		self.stlay = QtGui.QGridLayout(self.savetab)

		self.stsnapbut = QtGui.QPushButton("Snapshot")
		self.stsnapbut.setToolTip(".pgm, .ppm, .jpeg, .png, or .tiff format only")
		self.stwholebut = QtGui.QPushButton("Save Img")
		self.stwholebut.setToolTip("save EMData in any EMAN2 format")
		self.ststackbut = QtGui.QPushButton("Save Stack")
		self.ststackbut.setToolTip("save EMData objects as stack in any EMAN2 format")
		self.stmoviebut = QtGui.QPushButton("Movie")
		self.stanimgif = QtGui.QPushButton("GIF Anim")

		self.stlay.addWidget(self.stsnapbut,0,0)
		self.stlay.addWidget(self.stwholebut,1,0)
		self.stlay.addWidget(self.ststackbut,1,1)
		self.stlay.addWidget(self.stmoviebut,0,2)
		self.stlay.addWidget(self.stanimgif,1,2)
		self.stmoviebut.setEnabled(False)
		self.stanimgif.setEnabled(False)

		self.rngbl = QtGui.QHBoxLayout()
		self.stlay.addLayout(self.rngbl,2,0,1,3)
		self.stmmlbl = QtGui.QLabel("Img Range :")
		self.stminsb = QtGui.QSpinBox()
		self.stminsb.setRange(0,0)
		self.stminsb.setValue(0)
		self.stmaxsb = QtGui.QSpinBox()
		self.stmaxsb.setRange(0,0)
		self.stmaxsb.setValue(0)

		self.rngbl.addWidget(self.stmmlbl)
		self.rngbl.addWidget(self.stminsb)
		self.rngbl.addWidget(self.stmaxsb)

		self.mmtab.addTab(self.savetab,"Save")

		self.stsnapbut.clicked[bool].connect(self.do_snapshot)
		self.stwholebut.clicked[bool].connect(self.do_saveimg)
		self.ststackbut.clicked[bool].connect(self.do_savestack)
		self.stmoviebut.clicked[bool].connect(self.do_makemovie)
		self.stanimgif.clicked[bool].connect(self.do_makegifanim)

		# Filter tab
		self.filttab = QtGui.QWidget()
		self.ftlay=QtGui.QGridLayout(self.filttab)

		self.procbox1=StringBox(label="Process1:",value="filter.lowpass.gauss:cutoff_abs=0.125",showenable=0)
		self.ftlay.addWidget(self.procbox1,2,0)

		self.procbox2=StringBox(label="Process2:",value="filter.highpass.gauss:cutoff_pixels=3",showenable=0)
		self.ftlay.addWidget(self.procbox2,6,0)

		self.procbox3=StringBox(label="Process3:",value="math.linear:scale=5:shift=0",showenable=0)
		self.ftlay.addWidget(self.procbox3,10,0)

		self.proclbl1=QtGui.QLabel("Image unchanged, display only!")
		self.ftlay.addWidget(self.proclbl1,12,0)

		self.mmtab.addTab(self.filttab,"Filt")

		self.procbox1.enableChanged.connect(self.do_filters)
		self.procbox1.textChanged.connect(self.do_filters)
		self.procbox2.enableChanged.connect(self.do_filters)
		self.procbox2.textChanged.connect(self.do_filters)
		self.procbox3.enableChanged.connect(self.do_filters)
		self.procbox3.textChanged.connect(self.do_filters)

		# Probe tab
		self.probetab = QtGui.QWidget()
		self.ptlay=QtGui.QGridLayout(self.probetab)

		self.ptareasize= ValBox(label="Probe Size:",value=32)
		self.ptareasize.setIntonly(True)
		self.ptlay.addWidget(self.ptareasize,0,0,1,2)

		self.ptpointval= QtGui.QLabel("Point Value (ctr pix): ")
		self.ptlay.addWidget(self.ptpointval,1,0,1,2,Qt.AlignLeft)

		self.ptareaavg= QtGui.QLabel("Area Avg: ")
		self.ptlay.addWidget(self.ptareaavg,2,0,Qt.AlignLeft)

		self.ptareaavgnz= QtGui.QLabel("Area Avg (!=0): ")
		self.ptlay.addWidget(self.ptareaavgnz,2,1,Qt.AlignLeft)

		self.ptareasig= QtGui.QLabel("Area Sig: ")
		self.ptlay.addWidget(self.ptareasig,3,0,Qt.AlignLeft)

		self.ptareasignz= QtGui.QLabel("Area Sig (!=0): ")
		self.ptlay.addWidget(self.ptareasignz,3,1,Qt.AlignLeft)

		self.ptareaskew= QtGui.QLabel("Skewness: ")
		self.ptlay.addWidget(self.ptareaskew,4,0,Qt.AlignLeft)

		self.ptcoord= QtGui.QLabel("Center Coord: ")
		self.ptlay.addWidget(self.ptcoord,4,1,Qt.AlignLeft)

		self.ptareakurt= QtGui.QLabel("Kurtosis: ")
		self.ptlay.addWidget(self.ptareakurt,5,0,Qt.AlignLeft)

		self.ptcoord2= QtGui.QLabel("( ) ")
		self.ptlay.addWidget(self.ptcoord2,5,1,Qt.AlignLeft)

		# not really necessary since the pointbox accurately labels the pixel when zoomed in
		#self.ptpixels= QtGui.QWidget()
		#self.ptlay.addWidget(self.ptpixels,0,2)

		self.mmtab.addTab(self.probetab,"Probe")

		# Measure tab
		self.meastab = QtGui.QWidget()
		self.mtlay = QtGui.QGridLayout(self.meastab)

		#self.mtl1= QtGui.QLabel("A/Pix")
		#self.mtl1.setAlignment(Qt.AlignRight)
		#self.mtlay.addWidget(self.mtl1,0,0)

		self.mtapix = ValSlider(self,label="A/Pix")
		self.mtapix.setRange(0.5,10.0)
		self.mtapix.setValue(1.0)
		self.mtlay.addWidget(self.mtapix,0,0,1,2)
#		self.mtapix.setSizePolicy(QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed)
#		self.mtapix.resize(60,21)
#		print self.mtapix.sizeHint().width(),self.mtapix.sizeHint().height()


		self.mtshoworigin= QtGui.QLabel("Origin: 0,0")
		self.mtlay.addWidget(self.mtshoworigin,1,0,Qt.AlignLeft)

		self.mtshowend= QtGui.QLabel("End: 0,0")
		self.mtlay.addWidget(self.mtshowend,1,1,Qt.AlignLeft)

		self.mtshowlen= QtGui.QLabel("dx,dy: 0")
		self.mtlay.addWidget(self.mtshowlen,2,0,Qt.AlignLeft)

		self.mtshowlen2= QtGui.QLabel("Length: 0")
		self.mtlay.addWidget(self.mtshowlen2,2,1,Qt.AlignLeft)

		self.mtshowval= QtGui.QLabel("Value: ?")
		self.mtlay.addWidget(self.mtshowval,3,0,1,2,Qt.AlignLeft)

		self.mtshowval2= QtGui.QLabel(" ")
		self.mtlay.addWidget(self.mtshowval2,4,0,1,2,Qt.AlignLeft)


		self.mmtab.addTab(self.meastab,"Meas")

		# Draw tab
		self.drawtab = QtGui.QWidget()
		self.drawlay = QtGui.QGridLayout(self.drawtab)

		self.dtl1 = QtGui.QLabel("Pen Size:")
		self.dtl1.setAlignment(Qt.AlignRight)
		self.drawlay.addWidget(self.dtl1,0,0)

		self.dtpen = QtGui.QLineEdit("5")
		self.drawlay.addWidget(self.dtpen,0,1)

		self.dtl2 = QtGui.QLabel("Pen Val:")
		self.dtl2.setAlignment(Qt.AlignRight)
		self.drawlay.addWidget(self.dtl2,1,0)

		self.dtpenv = QtGui.QLineEdit("1.0")
		self.drawlay.addWidget(self.dtpenv,1,1)

		self.dtl3 = QtGui.QLabel("Pen Size2:")
		self.dtl3.setAlignment(Qt.AlignRight)
		self.drawlay.addWidget(self.dtl3,0,2)

		self.dtpen2 = QtGui.QLineEdit("5")
		self.drawlay.addWidget(self.dtpen2,0,3)

		self.dtl4 = QtGui.QLabel("Pen Val2:")
		self.dtl4.setAlignment(Qt.AlignRight)
		self.drawlay.addWidget(self.dtl4,1,2)

		self.dtpenv2 = QtGui.QLineEdit("0")
		self.drawlay.addWidget(self.dtpenv2,1,3)

		self.mmtab.addTab(self.drawtab,"Draw")

		# PSpec tab
		self.pstab = QtGui.QWidget()
		self.pstlay = QtGui.QGridLayout(self.pstab)

		self.psbsing = QtGui.QPushButton("Single")
		self.pstlay.addWidget(self.psbsing,0,0)

		self.psbstack = QtGui.QPushButton("Stack")
		self.pstlay.addWidget(self.psbstack,0,1)

		self.mmtab.addTab(self.pstab,"PSpec")
		self.pspecwins=[]

		# Python tab
		self.pytab = QtGui.QWidget()
		self.pytlay = QtGui.QGridLayout(self.pytab)

		self.pyinp = QtGui.QLineEdit()
		self.pytlay.addWidget(self.pyinp,0,0)

		self.pyout = QtGui.QTextEdit()
		self.pyout.setText("The displayed image is 'img'.\nEnter an expression above, like img['sigma']")
		self.pyout.setReadOnly(True)
		self.pytlay.addWidget(self.pyout,1,0,4,1)

		self.mmtab.addTab(self.pytab,"Python")


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
		self.hbl = QtGui.QHBoxLayout()
		self.hbl.setMargin(0)
		self.hbl.setSpacing(6)
		self.hbl.setObjectName("hbl")
		self.vbl.addLayout(self.hbl)

		self.hist = ImgHistogram(self)
		self.hist.setObjectName("hist")
		self.hbl.addWidget(self.hist)

		# Buttons next to the histogram
		self.vbl2 = QtGui.QGridLayout()
		self.vbl2.setMargin(0)
		self.vbl2.setSpacing(6)
		self.vbl2.setObjectName("vbl2")
		self.hbl.addLayout(self.vbl2)

		self.invtog = QtGui.QPushButton("Invert")
		self.invtog.setCheckable(1)
		self.vbl2.addWidget(self.invtog,0,0,1,1)#0012


		self.histoequal = QtGui.QComboBox(self)
		self.histoequal.addItem("Normal")
		self.histoequal.addItem("Hist Flat")
		self.histoequal.addItem("Hist Gauss")
		self.vbl2.addWidget(self.histoequal,0,1,1,1)

		self.auto_contrast_button = QtGui.QPushButton("Auto contrast")
		self.vbl2.addWidget(self.auto_contrast_button,1,0,1,2)

		# FFT Buttons
		self.fftg=QtWidgets.QButtonGroup()
		self.fftg.setExclusive(1)

		self.ffttog0 = QtGui.QPushButton("Real")
		self.ffttog0.setCheckable(1)
		self.ffttog0.setChecked(1)
		self.vbl2.addWidget(self.ffttog0,2,0)
		self.fftg.addButton(self.ffttog0,0)

		self.ffttog1 = QtGui.QPushButton("FFT")
		self.ffttog1.setCheckable(1)
		self.vbl2.addWidget(self.ffttog1,2,1)
		self.fftg.addButton(self.ffttog1,1)

		self.ffttog2 = QtGui.QPushButton("Amp")
		self.ffttog2.setCheckable(1)
		self.vbl2.addWidget(self.ffttog2,3,0)
		self.fftg.addButton(self.ffttog2,2)

		self.ffttog3 = QtGui.QPushButton("Pha")
		self.ffttog3.setCheckable(1)
		self.vbl2.addWidget(self.ffttog3,3,1)
		self.fftg.addButton(self.ffttog3,3)

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

		self.gammas = ValSlider(self,(.1,5.0),"Gam:")
		self.gammas.setObjectName("gamma")
		self.gammas.setValue(self.target().get_gamma())
		#self.gammas.setValue(1.0)
		self.vbl.addWidget(self.gammas)

		self.setWindowIcon(QtGui.QIcon(get_image_directory() +"eman.png"))

		self.lowlim=0
		self.highlim=1.0
		self.image_range = None
		#self.update_min_max()
		#self.update_brightness_contrast()
		self.busy=0

		self.psbsing.clicked[bool].connect(self.do_pspec_single)
		self.psbstack.clicked[bool].connect(self.do_pspec_stack)
		self.scale.valueChanged.connect(target.set_scale)
		self.mins.valueChanged.connect(self.new_min)
		self.maxs.valueChanged.connect(self.new_max)
		self.brts.valueChanged.connect(self.new_brt)
		self.conts.valueChanged.connect(self.new_cont)
		self.gammas.valueChanged.connect(self.new_gamma)
		self.pyinp.returnPressed.connect(self.do_python)
		self.invtog.toggled[bool].connect(target.set_invert)
		self.histoequal.currentIndexChanged[int].connect(target.set_histogram)
		self.fftg.buttonClicked[int].connect(target.set_FFT)
		self.mmtab.currentChanged[int].connect(target.set_mouse_mode)
		self.auto_contrast_button.clicked[bool].connect(target.auto_contrast)

		self.resize(400,440) # d.woolford thinks this is a good starting size as of Nov 2008 (especially on MAC)

	def do_pspec_single(self,ign):
		"""Compute 1D power spectrum of single image and plot"""
		try: data=self.target().list_data[self.target().list_idx]
		except: data=self.target().get_data()
		if data==None: return
#		print data
		fft=data.do_fft()
		pspec=fft.calc_radial_dist(old_div(fft["ny"],2),0.0,1.0,1)
		ds=old_div(1.0,(fft["ny"]*data["apix_x"]))
		s=[ds*i for i in range(old_div(fft["ny"],2))]

		from .emplot2d import EMDataFnPlotter

		dfp=EMDataFnPlotter(data=(s,pspec))
		dfp.show()
		self.pspecwins.append(dfp)

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
		fsp=QtGui.QFileDialog.getSaveFileName(self, "Select output file, .pgm, .ppm, .jpeg, .png or .tiff only")
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
		fsp=QtGui.QFileDialog.getSaveFileName(self, "Select output file, format extrapolated from file extenstion")
		fsp=str(fsp)
		self.target().data.write_image(fsp,-1)

	def do_savestack(self,du) :
		if self.target().list_data==None : return
		fsp=str(QtGui.QFileDialog.getSaveFileName(self, "Select root output file, format extrapolated from file extenstion"))
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
		fsp=QtGui.QFileDialog.getSaveFileName(self, "Select output file, format extrapolated from file extenstion. ffmpeg must be installed")
		if self.target().list_data==None : return

		for i in range(self.stminsb.value()-1,self.stmaxsb.value()):
			im=self.target().list_data[i]
			im["render_min"]=im["mean"]-im["sigma"]*2.5
			im["render_max"]=im["mean"]+im["sigma"]*2.5
			im.write_image("tmp.%03d.png"%(i-self.stminsb.value()+1))

		# vcodec and pix_fmt are for quicktime compatibility. r 2 is 2 FPS
		ret= os.system("ffmpeg -r 5 -i tmp.%%03d.png -vcodec libx264 -pix_fmt yuv420p  %s"%fsp)
		if ret!=0 :
			QtGui.QMessageBox.warning(None,"Error","Movie conversion (ffmpeg) failed. Please make sure ffmpeg is in your path. Frames not deleted.")
			return

		for i in range(self.stminsb.value()-1,self.stmaxsb.value()):
			os.unlink("tmp.%03d.png"%(i-self.stminsb.value()+1))


	def do_makegifanim(self,du) :
		fsp=QtGui.QFileDialog.getSaveFileName(self, "Select output gif file")
		if self.target().list_data==None : return

		for i in range(self.stminsb.value()-1,self.stmaxsb.value()):
			im=self.target().list_data[i]
			im.write_image("tmp.%03d.png"%(i-self.stminsb.value()+1))

		ret= os.system("convert tmp.???.png %s"%fsp)
		if ret!=0 :
			QtGui.QMessageBox.warning(None,"Error","GIF conversion failed. Please make sure ImageMagick (convert program) is installed and in your path. Frames not deleted.")
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
			r="Error executing. Access the image as 'img', for example \nimg['mean'] will yield the mean image value"
#		print r

		self.pyout.setText(str(r))

	def disable_image_range(self):
		self.stminsb.setRange(0,0)
		self.stmaxsb.setRange(0,0)
		self.stmoviebut.setEnabled(False)
		self.stanimgif.setEnabled(False)

		if self.image_range != None:
			self.vbl.removeWidget(self.image_range)
			self.image_range.deleteLater()
			self.image_range = None
		else:
			# this is fine
			pass
			#print "warning, attempted to disable image range when there was none!"

	def enable_image_range(self,minimum,maximum,current_idx):
		if self.image_range == None:
			self.image_range = ValSlider(self,label="N#:")
			self.image_range.setIntonly(True)
			self.vbl.addWidget(self.image_range)

		self.image_range.setRange(minimum,maximum)
		self.image_range.setValue(current_idx)
		self.stmoviebut.setEnabled(True)
		self.stanimgif.setEnabled(True)

		self.stminsb.setRange(minimum,maximum)
		self.stmaxsb.setRange(minimum,maximum)
		self.stmaxsb.setValue(maximum)

		self.image_range.valueChanged.connect(self.target().image_range_changed)

	def set_image_idx(self,val):
		self.image_range.setValue(val)

	def set_fft_amp_pressed(self):
		self.ffttog2.setChecked(1)

	def get_contrast(self):
		return float(self.conts.getValue())

	def get_brightness(self):
		return float(self.brts.getValue())

	#def set_contrast(self,value,quiet=1):
		#self.conts.setValue(value,quiet)

	#def set_brightness(self,value,quiet=1):
		#self.brts.setValue(value,quiet)

	def set_maxden(self,value,quiet=1):
		self.maxs.setValue(value,quiet)

	def set_minden(self,value,quiet=1):
		self.mins.setValue(value,quiet)

	def set_gamma(self,value,quiet=1):
		self.gammas.setValue(value,quiet)

	def set_scale(self,val):
		if self.busy : return
		self.busy=1
		self.scale.setValue(val)
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

	def new_gamma(self,val):
		if self.busy : return
		self.busy=1
		self.target().set_gamma(val)
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


# This is just for testing, of course
def main():
	from .emapplication import EMApp
	em_app = EMApp()
	window = EMImage2DWidget(application=em_app)

	if len(sys.argv)==1 :
		window.set_data(test_image(size=(128,128)))
	else :
		a=EMData.read_images(sys.argv[1])
		if len(a) == 1:	a = a[0]
		window.set_data(a,sys.argv[1])

	em_app.show()
	window.optimally_resize()
	sys.exit(em_app.exec_())


if __name__ == '__main__':
	main()
