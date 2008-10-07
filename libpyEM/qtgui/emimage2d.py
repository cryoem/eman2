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

import PyQt4
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from OpenGL import GL,GLU,GLUT
from OpenGL.GL import *
from valslider import ValSlider
from math import *
from EMAN2 import *
import EMAN2
import sys
import numpy
import struct
from emimageutil import ImgHistogram,EMEventRerouter
import emshape 
from emshape import EMShape
from weakref import WeakKeyDictionary
from pickle import dumps,loads

from emglobjects import EMOpenGLFlagsAndTools, EMImage2DGUIModule
from emapplication import EMStandAloneApplication, EMQtWidgetModule, EMGUIModule

try: from PyQt4 import QtWebKit
except: pass
import platform

MAG_INC = 1.1

from emglobjects import EMOpenGLFlagsAndTools

GLUT.glutInit(sys.argv )

class EMImage2DWidget(QtOpenGL.QGLWidget,EMEventRerouter):
	allim=WeakKeyDictionary()
	def __init__(self, em_image_2d_module,enable_timer=False):
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True)
		fmt.setSampleBuffers(True)
		fmt.setDepth(1)
		QtOpenGL.QGLWidget.__init__(self,fmt)
		EMEventRerouter.__init__(self,em_image_2d_module)
		self.initimageflag = True
		
		self.setFocusPolicy(Qt.StrongFocus)
		
		if enable_timer:
			self.enable_timer()
		else: self.time_enabled = False
		
	def enable_timer(self):
		if self.time_enabled == False:
			self.timer = QtCore.QTimer()
			QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.timeout)
			self.timeinterval = 50
			self.timer.start(50)
			self.time_enabled = True
		else: print "timer already enabled in emimage2d"
		
	def set_parent(self,parent):
		self.parent = parent
		
	def timeout(self):
		update = False
		if self.target.update_blend() :
			self.target.shapechange = 1
			update = True
		
		if self.target.update_animation():
			update = True
		
		if update:
			self.updateGL()

	def set_data(self,data):
		self.target.set_data(data)
		
	def initializeGL(self):
		GL.glClearColor(0,0,0,0)
		
		glLightfv(GL_LIGHT0, GL_AMBIENT, [0.1, 0.1, 0.1, 1.0])
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_POSITION,  [.1,.1,1,0.])
	
		glEnable(GL_LIGHTING)
		glEnable(GL_LIGHT0)
		try:
			self.target.initializeGL()
			self.initimageflag = False
		except:
			pass
	
	def paintGL(self):
		if not self.target: return
		
		glClear(GL_COLOR_BUFFER_BIT)
		if glIsEnabled(GL_DEPTH_TEST):
			glClear(GL_DEPTH_BUFFER_BIT)
		if glIsEnabled(GL_STENCIL_TEST):
			glClear(GL_STENCIL_BUFFER_BIT)

		
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		#self.cam.position()
		if self.initimageflag == True:
			self.target.initializeGL()
			self.initimageflag = False
		#context = OpenGL.contextdata.getContext(None)
		#print "Image2D context is", context
		glPushMatrix()
		self.target.render()
		glPopMatrix()
		
	def resizeGL(self, width, height):
		side = min(width, height)
		GL.glViewport(0,0,self.width(),self.height())
	
		GL.glMatrixMode(GL.GL_PROJECTION)
		GL.glLoadIdentity()
		GLU.gluOrtho2D(0.0,self.width(),0.0,self.height())
		GL.glMatrixMode(GL.GL_MODELVIEW)
		GL.glLoadIdentity()
		
		try: self.target.resizeEvent(width,height)
		except: pass
		
	def add_shapes(self,s):
		self.target.add_shapes(s)
		
	def add_shape(self,name,shape):
		return self.target.add_shape(name,shape)

	def scr_to_img(self,p):
		return self.target.scr_to_img(p)
	
	def set_active(self,a,b,c,d):
		return self.target.set_active(a,b,c,d)
	
	def get_shapes(self):
		return self.target.get_shapes()
	
	def del_shapes(self):
		return self.target.del_shapes()
	
	def del_shape(self,p):
		return self.target.del_shape(p)

	def scroll_to(self,x,y):
		return self.target.scroll_to(x,y)
	
	def register_scroll_motion(self,x,y):
		return self.target.register_scroll_motion(x,y)
	
	def get_depth_for_height(self, height):
		return 0
	
	def get_shapes(self):
		return self.target.get_shapes()
	
class EMImage2DModule(EMImage2DGUIModule):
	"""
	"""
	allim=WeakKeyDictionary()
	
	def get_qt_widget(self):
		if self.parent == None:	
			self.parent = EMImage2DWidget(self)
			self.parent.setWindowTitle(self.file_name)
			if isinstance(self.data,EMData):
				self.__parent_resize()
		return self.parent
	
	def __parent_resize(self):
		try:
			if self.parent_geometry != None:
				self.parent.restoreGeometry(self.parent_geometry)
			elif self.data.get_xsize()<1024 and self.data.get_ysize()<1024: self.parent.resize(self.data.get_xsize(),self.data.get_ysize())
			else: self.parent.resize(800,800)
			self.init_size_flag = False
		except: pass
			
	
	def __init__(self, image=None,application=None):
		EMImage2DGUIModule.__init__(self,application)
		
		self.parent = None
		self.data=image				# EMData object to display
		self.oldsize=(-1,-1)
		self.datasize=(1,1)			# Dimensions of current image
		self.scale=1.0				# Scale factor for display
		self.origin=(0,0)			# Current display origin
		self.invert=0				# invert image on display
		self.gamma=1.0				# gamma for display (impact on inverted contrast ?
		self.minden=0
		self.maxden=1.0
		self.mindeng=0
		self.maxdeng=1.0
		self.brts=0 # stored purely for persistence reasons - the database needs this object to remember the value
		self.conts=0 # stored purely for persistence reasons - the database needs this object to remember the value
		self.fft=None				# The FFT of the current target if currently displayed
		self.rmousedrag=None		# coordinates during a right-drag operation
		self.mmode=0				# current mouse mode as selected by the inspector
		self.curfft=0				# current FFT mode (when starting with real images only)
		self.mag = 1.1				# magnification factor
		self.invmag = 1.0/self.mag	# inverse magnification factor
		
		self.shapes={}				# dictionary of shapes to draw, see add_shapes
		self.shapechange=1			# Set to 1 when shapes need to be redrawn
		self.active=(None,0,0,0)	# The active shape and a hilight color (n,r,g,b)
		
		self.extras = []			# an empty set of extras - other images that can be rendered over this one
		
		self.startorigin = None
		self.endorigin = None
		self.isanimated = False
		self.time = 1
		self.timeinc = 0.125
		
		self.inspector=None			# set to inspector panel widget when exists
		self.em_qt_inspector_widget = None # eventually a EMQtWidgetModule storing the inspector
		
		self.init_size = True		# A flag used to set the initial origin offset
		
		self.shapelist = 0			# a display list identify
		
		self.glflags = EMOpenGLFlagsAndTools() 	# supplies power of two texturing flags
		
		self.suppress_inspector = False 	# Suppresses showing the inspector - switched on in emfloatingwidgets
		self.tex_name = 0			# an OpenGL texture handle
		
		self.otherdata = None
		self.otherdatascale = -1
		self.otherdatablend = False
		self.other_tex_name = None
		self.init_size_flag = True
		self.frozen = False
		self.isexcluded = False
		self.hack_shrink = 1
		self.parent_geometry = None
		
		self.list_data = None # this can be used for viewing lists of data
		self.list_idx = 0	# and idx to the list_data
		
		self.use_display_list = True # whether or not a display list should be used to render the image pixelsw - if on, this will save on time if the view of the image is unchanged, which can quite often be the case
		self.main_display_list = 0	# if using display lists, the stores the display list
		self.display_states = [] # if using display lists, this stores the states that are checked, and if different, will cause regeneration of the display list
		self.hist = []
		self.file_name = ""# stores the filename of the image, if None then member functions should be smart enough to handle it
		
		self.wheel_navigate = False # useful on Mac laptops
		
		try: self.parent.setAcceptDrops(True)
		except:	pass
		

		if image : self.set_data(image)
		else:self.__load_display_settings_from_db()
	def __del__(self):
		if (self.shapelist != 0):
			glDeleteLists(self.shapelist,1)
			self.shapelist = 0
		if self.main_display_list != 0:
			glDeleteLists(self.main_display_list,1)
			self.main_display_list = 0
	
	def set_mmode(self,mode): self.mmode = mode
	def get_minden(self): return self.minden
	def get_maxden(self): return self.maxden
	def get_gamma(self): return self.gamma
	def get_brts(self): return self.brts
	def get_conts(self): return self.conts
	def get_shapes(self): return self.shapes
	
	def set_brightness(self,brts):
		self.brts = brts
		self.updateGL()
		
	def set_contrast(self,conts):
		self.conts = conts
		self.updateGL()
		
	def set_brightness_contrast(self,brts,conts):
		'''
		Setting the brightness contrast does not require an updateGL but does cause information to
		be stored to the data base
		'''
		self.brts = brts
		self.const = conts
		self.updateGL()
		
	def set_density_range(self,x0,x1):
		"""Set the range of densities to be mapped to the 0-255 pixel value range"""
		self.minden=x0
		self.maxden=x1
		#self.__write_display_settings_to_db()
		self.updateGL()
	
	def set_density_min(self,val):
		self.minden=val
		self.updateGL()
		
	def set_density_max(self,val):
		self.maxden=val
		#self.__write_display_settings_to_db()
		self.updateGL()
	
	def set_gamma(self,val):
		self.gamma=val
		#self.__write_display_settings_to_db()
		self.updateGL()
	
	def set_file_name(self,file_name):
		self.file_name = file_name
		try:self.parent.setWindowTitle(file_name)
		except:pass
		self.__load_display_settings_from_db()
		
	def get_file_name(self):
		return self.file_name
	
	def set_other_data(self,data,scale,blend=False):
		self.otherdata = data
		self.otherdatascale = scale
		self.otherdatablend = blend
	
	def width(self):
		try:
			return self.data.get_xsize()
		except:
			return 0
	
	def height(self):
		try:
			return self.data.get_ysize()
		except:
			return 0
	
	def updateGL(self):
		try: self.parent.updateGL()
		except: pass
		
	def set_frozen(self,frozen):
		self.frozen = frozen
		
	def set_excluded(self,isexcluded):
		wasexcluded = self.isexcluded
		
		self.isexcluded = isexcluded
		
		if wasexcluded or self.isexcluded: return True
		else: return False
	
	def set_data(self,data):
		"""You may pass a single 2D image"""
		if self.data != None:
			self.__write_display_settings_to_db()
	
		self.data=data
		if data==None:
			self.list_data = None
			self.updateGL()
			return
		elif isinstance(data,list):
			self.list_data = data
			self.list_idx = len(self.list_data)/2
			self.data = data[self.list_idx]
		else:
			self.list_data = None
		
		if self.init_size_flag and isinstance(self.parent,QtGui.QWidget):
			self.__parent_resize()
		
		mean=self.data.get_attr("mean")
		sigma=self.data.get_attr("sigma")
		m0=self.data.get_attr("minimum")
		m1=self.data.get_attr("maximum")
		
		self.minden=max(m0,mean-3.0*sigma)
		self.maxden=min(m1,mean+3.0*sigma)
		self.mindeng=max(m0,mean-5.0*sigma)
		self.maxdeng=min(m1,mean+5.0*sigma)

		self.datasize=(self.data.get_xsize(),self.data.get_ysize())
		self.scale=1.0				# Scale factor for display
		self.origin=(0,0)			# Current display origin
		self.__load_display_settings_from_db()
		if self.curfft : 
			self.set_FFT(self.curfft)

	def __load_display_settings_from_db(self):
		if self.file_name == "": return # there is no file name, we have no means to stores information
		
		try:
			DB = EMAN2db.EMAN2DB.open_db(".")
			DB.open_dict("image_2d_display_settings")
		except:
			# Databasing is not supported, in which case w
			return
		
		db = DB.image_2d_display_settings
	
		data = db[self.file_name]
		if data == None:
			data = db["latest_display_settings"] # if there isn't already information try for the latest
			if data == None:
				return # there are no settings we can use
			else:
				self.__write_display_settings_to_db() # we should store the information if we are suddenly using it
		self.minden = data["min"]
		self.maxden = data["max"]
		self.gamma = data["gamma"]
		self.scale = data["scale"] 
		self.origin = data["origin"]
		try:
			self.parent_geometry = data["parent_geometry"]
			if self.parent != None:
				try: self.parent.restoreGeometry(self.parent_geometry)
				except: pass
		except:pass
		
		self.inspector_update()
		self.force_display_update()
	
	def __write_display_settings_to_db(self):
		'''
		writes the min,max, brightness, contrast and gamma values associated with
		the current image to the homedb. The full path of the image is used
		'''
		
		if self.file_name == None: return # there is no file name, we have no means to stores information
		
		try:
			DB = EMAN2db.EMAN2DB.open_db(".")
			DB.open_dict("image_2d_display_settings")
		except:
			# Databasing is not supported, in which case we do nothing
			return
		
		data = {}	
		data["min"] = self.minden
		data["max"] = self.maxden
		data["gamma"] = self.gamma
		data["origin"] = self.origin
		data["scale"] = self.scale
		try: data["parent_geometry"] = self.parent.saveGeometry()
		except: pass
		
		db = DB.image_2d_display_settings
		db[self.file_name] = data
		db["latest_display_settings"] = data #store this to automatically apply previously used settings to other images - this was originally a request of Yao Cong


	def set_origin(self,x,y):
		"""Set the display origin within the image"""
		self.origin=(x,y)
		self.updateGL()
	
	def scroll_to(self,x,y):
		"""center the point on the screen"""
		self.set_origin(x*self.scale-self.parent.width()/2,y*self.scale-self.parent.height()/2)

	def set_shapes(self,shapes,shrink):
		self.hack_shrink = shrink
		self.shapes = shapes
		self.shapechange=1
		
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
		
	def register_scroll_motion(self,x,y):
		self.startorigin = self.origin
		self.endorigin = (x*self.scale-self.parent.width()/2,y*self.scale-self.parent.height()/2)
		self.isanimated = True
		self.time = 0
		return True

	def set_scale(self,newscale):
		"""Adjusts the scale of the display. Tries to maintain the center of the image at the center"""
		self.origin=(newscale/self.scale*(self.parent.width()/2.0+self.origin[0])-self.parent.width()/2.0,newscale/self.scale*(self.parent.height()/2.0+self.origin[1])-self.parent.height()/2.0)
		self.scale=newscale
		self.updateGL()
		
	def set_invert(self,val):
		if val: self.invert=1
		else : self.invert=0
		self.updateGL()
		
	def set_FFT(self,val):
		if self.data.is_complex(): return
		self.curfft=val
		if val>0 :
			try:
				self.fft=self.data.do_fft()
				self.fft.set_value_at(0,0,0,0)
				self.fft.set_value_at(1,0,0,0)
				if val==1 :
					self.fft.process_inplace("xform.phaseorigin.tocorner")
				elif val==2 :
					self.fft.process_inplace("xform.fourierorigin.tocenter")
					self.fft=self.fft.get_fft_amplitude()
				elif val==3 :
					self.fft.process_inplace("xform.fourierorigin.tocenter")
					self.fft=self.fft.get_fft_phase()

			
				mean=self.fft.get_attr("mean")
				sigma=self.fft.get_attr("sigma")
				m0=self.fft.get_attr("minimum")
				m1=self.fft.get_attr("maximum")
			
				self.fminden=0
				self.fmaxden=min(m1,mean+5.0*sigma)
				self.fmindeng=0
				self.fmaxdeng=min(m1,mean+8.0*sigma)
				
				self.ominden=self.minden
				self.omaxden=self.maxden
				
				self.show_inspector()
			except: 
				self.fft=None
		else: 
			self.fft=None
			self.minden=self.ominden
			self.maxden=self.omaxden
			self.show_inspector()
		self.updateGL()

	def initializeGL(self):
		#GL.glClearColor(0,0,0,0)
		emshape.initGL()
		self.shapelist=GL.glGenLists(1)		# displaylist for shapes displayed over the image
		GL.glNewList(self.shapelist,GL.GL_COMPILE)
		GL.glEndList()

	def force_display_update(self):
		self.display_states = []

	def display_state_changed(self):
		
		display_states = []
		display_states.append(self.parent.width())
		display_states.append(self.parent.height())
		display_states.append(self.origin[0])
		display_states.append(self.origin[1])
		display_states.append(self.scale)
		display_states.append(self.invert)
		display_states.append(self.minden)
		display_states.append(self.maxden)
		display_states.append(self.gamma)
		display_states.append(self.curfft)
		if len(self.display_states) == 0:
			self.display_states = display_states
			return True
		else:
			for i in range(len(display_states)):
				
				if display_states[i] != self.display_states[i]:
					self.display_states = display_states
					return True
		
		return False

	def render(self):
		
		if not self.data : return
		
		lighting = glIsEnabled(GL_LIGHTING)
		glDisable(GL_LIGHTING)

		if self.shapechange:
			self.setupShapes()
			self.shapechange=0

		width = self.parent.width()/2.0
		height = self.parent.height()/2.0
		
		if not self.invert : pixden=(0,255)
		else: pixden=(255,0)
		
		if self.otherdata != None and isinstance(self.otherdata,EMData) and not self.glflags.npt_textures_unsupported():
			scale = self.scale*self.otherdatascale
			b=self.otherdata.render_amp8(int(self.origin[0]/scale),int(self.origin[1]/scale),self.parent.width(),self.parent.height(),(self.parent.width()-1)/4*4+4,scale,pixden[0],pixden[1],0,1,1,2)
			gl_render_type = GL_LUMINANCE
			
			if self.other_tex_name != 0: GL.glDeleteTextures(self.other_tex_name)
			self.other_tex_name = GL.glGenTextures(1)
			if ( self.other_tex_name <= 0 ):
				raise("failed to generate texture name")
			
			glBindTexture(GL.GL_TEXTURE_2D,self.other_tex_name)
			glTexImage2D(GL.GL_TEXTURE_2D,0,gl_render_type,self.parent.width(),self.parent.height(),0,gl_render_type, GL.GL_UNSIGNED_BYTE, b)
			
			glPushMatrix()
			glTranslatef(width,height,0)
			self.__draw_texture(self.other_tex_name,-width,-height,width,height)
			glPopMatrix()

		render = False
		if self.use_display_list:
			if self.display_state_changed():
				if self.main_display_list != 0:
					glDeleteLists(self.main_display_list,1)
					self.main_display_list = 0

			if self.main_display_list == 0:
				self.main_display_list = glGenLists(1)
				glNewList(self.main_display_list,GL_COMPILE)
				render = True
		else: render = True
		
		if render:
			if self.curfft==1 :
				if self.fft.is_complex() == False:
					print "error, the fft is not complex, internal error"
					return
				a=self.fft.render_ap24(int(self.origin[0]/self.scale),int(self.origin[1]/self.scale),self.parent.width(),self.parent.height(),(self.parent.width()*3-1)/4*4+4,self.scale,pixden[0],pixden[1],self.minden,self.maxden,self.gamma,3)
				gl_render_type = GL_RGB
				
			elif self.curfft in (2,3) :
				a=self.fft.render_amp8(int(self.origin[0]/self.scale),int(self.origin[1]/self.scale),self.parent.width(),self.parent.height(),(self.parent.width()-1)/4*4+4,self.scale,pixden[0],pixden[1],self.minden,self.maxden,self.gamma,2)
				gl_render_type = GL_LUMINANCE
			else : 
				a=self.data.render_amp8(int(self.origin[0]/self.scale),int(self.origin[1]/self.scale),self.parent.width(),self.parent.height(),(self.parent.width()-1)/4*4+4,self.scale,pixden[0],pixden[1],self.minden,self.maxden,self.gamma,2)
				gl_render_type = GL_LUMINANCE
			if  not self.glflags.npt_textures_unsupported():
				
				self.hist=struct.unpack('256i',a[-1024:])
			
				if self.tex_name != 0: glDeleteTextures(self.tex_name)
				self.tex_name = glGenTextures(1)
				if ( self.tex_name <= 0 ):
					raise("failed to generate texture name")
				
				if self.otherdatablend and self.otherdata != None:
					GL.glEnable(GL.GL_BLEND);
					depth_testing_was_on = GL.glIsEnabled(GL.GL_DEPTH_TEST);
					GL.glDisable(GL.GL_DEPTH_TEST);
					try:
						GL.glBlendEquation(GL.GL_FUNC_SUBTRACT);
					except: pass
					#GL.glBlendFunc(GL.GL_SRC_ALPHA,GL.GL_ONE_MINUS_SRC_ALPHA);
					GL.glBlendFunc(GL.GL_ONE,GL.GL_ONE);
	
				GL.glBindTexture(GL.GL_TEXTURE_2D,self.tex_name)
				GL.glTexImage2D(GL.GL_TEXTURE_2D,0,gl_render_type,self.parent.width(),self.parent.height(),0,gl_render_type, GL.GL_UNSIGNED_BYTE, a)
				
				glPushMatrix()
				glTranslatef(width,height,0)
				self.__draw_texture(self.tex_name,-width,-height,width,height)
				glPopMatrix()
				
				if self.otherdatablend and self.otherdata != None:
					GL.glDisable( GL.GL_BLEND);
					if (depth_testing_was_on):	GL.glEnable(GL.GL_DEPTH_TEST)
			
			else:
				GL.glRasterPos(0,self.parent.height()-1)
				GL.glPixelZoom(1.0,-1.0)
				GL.glDrawPixels(self.parent.width(),self.parent.height(),gl_render_type,GL.GL_UNSIGNED_BYTE,a)
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
	
		if self.shapelist != 0:
			GL.glPushMatrix()
			GL.glTranslate(-self.origin[0],-self.origin[1],0.1)
			GL.glScalef(self.scale,self.scale,1.0)
			if self.hack_shrink != 1:
	#			print len(self.shapes),self.hack_shrink
	#			print self.hack_shrink
				GL.glScale(1.0/self.hack_shrink,1.0/self.hack_shrink,1.0)
				#for k,s in self.shapes.items():
					#if self.active[0]==k: s.draw(None,self.active[1:])
					#else: s.draw()
			#else:
				##print self.shapelist
			GL.glCallList(self.shapelist)
			GL.glPopMatrix()
		self.changec=self.data.get_attr("changecount")
		

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
		glVertex2f(xmin,ymin,)
			
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
		
	def resizeEvent(self,width,height):
		if self.init_size :
			if self.origin == (0,0):
				self.origin = ((self.data.get_xsize() - self.parent.width())/2.0, (self.data.get_ysize() - self.parent.height())/2.0 )
				self.oldsize=(width,height)
			self.init_size = False
		else:
			if self.origin == (0,0):
				self.origin=((self.oldsize[0]/2.0+self.origin[0])-self.parent.width()/2.0,(self.oldsize[1]/2.0+self.origin[1])-self.parent.height()/2.0)
				self.oldsize=(width,height)
				
		self.window_width = width
		self.window_height = height

	def get_shapes(self):
		return self.shapes
	
	def setupShapes(self):
		if self.shapelist != 0: GL.glDeleteLists(self.shapelist,1)
		
		self.shapelist = GL.glGenLists(1)
		
		#context = OpenGL.contextdata.getContext(None)
		#print "Image2D context is", context,"display list is",self.shapelist
		
		# make our own cirle rather than use gluDisk or somesuch
		GL.glNewList(self.shapelist,GL.GL_COMPILE)
		for k,s in self.shapes.items():
			if self.active[0]==k: s.draw(None,self.active[1:])
			else: s.draw()
		GL.glEndList()
	
	def inspector_update(self):
		if self.inspector:
			self.inspector.set_minden(self.minden)
			self.inspector.set_maxden(self.maxden)
			self.inspector.set_gamma(self.gamma)
			self.inspector.set_scale(self.scale)
			self.inspector.update_brightness_contrast()
	
	def get_inspector(self):
		if not self.inspector:
			self.inspector=EMImageInspector2D(self)
			self.inspector_update()
		return self.inspector
	
	def set_mouse_mode(self,m):
		self.mmode=m
	
	def set_active(self,n,r,g,b):
		self.active=(n,r,g,b)
		self.shapechange=1
		#self.updateGL()

	def update_blend(self):
		ret = False
		for shape in self.shapes.items():
			s = shape[1]
			if s.isanimated:
				v = s.incblend()
				if v == 2:
					self.parent.emit(QtCore.SIGNAL("removeshape"), shape[0])
				ret = True
		
		return ret

	def add_shape(self,k,s):
		"""Add an EMShape object to be overlaid on the image. Each shape is
		keyed into a dictionary, so different types of shapes for different
		purposes may be simultaneously displayed.
		
		"""
		self.shapes[k]=s
		self.shapechange=1
		#self.updateGL()
	
	def add_shapes(self,d):
		self.shapes.update(d)
		self.shapechange=1
		#self.updateGL()
	
	
	def del_shape(self,p):
		self.shapes.pop(p)
		self.shapechange=1

	def del_shapes(self,k=None):
		if k:
			try:
				for i in k:
					del self.shapes[k]
			except: del self.shapes[k]
		else:
			self.shapes={}
		self.shapechange=1
		#self.updateGL()
	
	def scr_to_img(self,v0,v1=None):
		try: return ((v0+self.origin[0])/self.scale,(self.parent.height()-(v1-self.origin[1]))/self.scale)
		except:	return ((v0[0]+self.origin[0])/self.scale,(self.parent.height()-(v0[1]-self.origin[1]))/self.scale)

	def closeEvent(self,event) :
		self.__write_display_settings_to_db()
		if self.inspector: 
			self.application.close_child(self.em_qt_inspector_widget)
		
	def dragEnterEvent(self,event):
#		f=event.mimeData().formats()
#		for i in f:
#			print str(i)
		
		if event.provides("application/x-eman"):
			event.setDropAction(Qt.CopyAction)
			event.accept()

	def dropEvent(self,event):
#		lc=self.scr_to_img((event.pos().x(),event.pos().y()))
		if EMAN2.GUIbeingdragged:
			self.set_data(EMAN2.GUIbeingdragged)
			EMAN2.GUIbeingdragged=None
		elif event.provides("application/x-eman"):
			x=loads(event.mimeData().data("application/x-eman"))
			self.set_data(x)
			event.acceptProposedAction()

	def mousePressEvent(self, event):
		lc=self.scr_to_img(event.x(),event.y())
		if event.button()==Qt.MidButton or (event.button()==Qt.LeftButton and event.modifiers()&Qt.ControlModifier):
			self.show_inspector(1)
			self.parent.emit(QtCore.SIGNAL("inspector_shown"),event)
		elif event.button()==Qt.RightButton or (event.button()==Qt.LeftButton and event.modifiers()&Qt.AltModifier):
			app =  QtGui.QApplication.instance()
			try:
				app.setOverrideCursor(Qt.ClosedHandCursor)
			except: # if we're using a version of qt older than 4.2 than we have to use this...
				app.setOverrideCursor(Qt.SizeAllCursor)
			self.rmousedrag=(event.x(),event.y() )
		elif event.button()==Qt.LeftButton:
			if self.mmode==0:
				self.parent.emit(QtCore.SIGNAL("mousedown"), event)
				return
			elif self.mmode==1 :
				try: 
					del self.shapes["MEASL"]
				except: pass
				self.add_shape("MEAS",EMShape(("line",.5,.1,.5,lc[0],lc[1],lc[0]+1,lc[1],2)))
			elif self.mmode==2 and self.inspector:
				#try:
#					print "paint ",lc
					self.drawr1=int(float(self.inspector.dtpen.text()))
					self.drawv1=float(self.inspector.dtpenv.text())
					self.drawr2=int(float(self.inspector.dtpen2.text()))
					self.drawv2=float(self.inspector.dtpenv2.text())
					self.data.process_inplace("mask.paint",{"x":lc[0],"y":lc[1],"z":0,"r1":self.drawr1,"v1":self.drawv1,"r2":self.drawr2,"v2":self.drawv2})
					self.parent.update()
				#except:
					#print "paint error"
					#return
				
	def mouseMoveEvent(self, event):
		lc=self.scr_to_img(event.x(),event.y())
		if self.rmousedrag:
			self.origin=(self.origin[0]+self.rmousedrag[0]-event.x(),self.origin[1]-self.rmousedrag[1]+event.y())
			self.rmousedrag=(event.x(),event.y())
			try: self.parent.update()
			except: pass
		elif event.buttons()&Qt.LeftButton:
			if self.mmode==0:
				self.parent.emit(QtCore.SIGNAL("mousedrag"), event)
				return
			elif self.mmode==1 :
				self.add_shape("MEAS",EMShape(("line",.5,.1,.5,self.shapes["MEAS"].shape[4],self.shapes["MEAS"].shape[5],lc[0],lc[1],2)))
				dx=lc[0]-self.shapes["MEAS"].shape[4]
				dy=lc[1]-self.shapes["MEAS"].shape[5]
				self.add_shape("MEASL",EMShape(("label",.1,.1,.1,lc[0]+2,lc[1]+2,"%d,%d - %d,%d"%(self.shapes["MEAS"].shape[4],self.shapes["MEAS"].shape[5],lc[0],lc[1]),9,-1)))
				if self.inspector:
					apix=self.inspector.mtapix.value
					self.inspector.mtshoworigin.setText("Start: %d , %d"%(self.shapes["MEAS"].shape[4],self.shapes["MEAS"].shape[5]))
					self.inspector.mtshowend.setText("  End: %d , %d"%(lc[0],lc[1]))
					self.inspector.mtshowlen.setText("dx,dy (len): %1.2f , %1.2f (%1.3f)"%(dx*apix,dy*apix,hypot(dx,dy)*apix))
			elif self.mmode==2 and self.inspector:
				self.data.process_inplace("mask.paint",{"x":lc[0],"y":lc[1],"z":0,"r1":self.drawr1,"v1":self.drawv1,"r2":self.drawr2,"v2":self.drawv2})
				self.parent.update()
		elif self.mmode == 0:
			self.parent.emit(QtCore.SIGNAL("mousemove"), event)
	
	def keyPressEvent(self,event):
		if event.key() == Qt.Key_F1:
			try:
				try:
					test = self.browser
				except: 
					self.browser = QtWebKit.QWebView()
					self.browser.load(QtCore.QUrl("http://blake.bcm.edu/emanwiki/e2display"))
					self.browser.resize(800,800)
				
				if not self.browser.isVisible(): self.browser.show()
			except:
				print "in the middle of getting help working"
				self.browser2 = QtGui.QTextBrowser()
				#url = QtCore.QUrl("http://blake.bcm.edu/emanwiki/e2display")
				url = QtCore.QUrl("http://www.google.com")
				url.setPort(80)
				#print url.port()
				self.browser2.setSource(url)
				#print browser2.port()
				self.browser2.show()
				self.browser2.resize(800,800)
				
		elif event.key() == Qt.Key_Up:
			if self.list_data != None:
				if (self.list_idx < (len(self.list_data)-1)):
					self.list_idx += 1
					self.data = self.list_data[self.list_idx]
					self.force_display_update()
					self.updateGL()
		elif event.key() == Qt.Key_Down:
			if self.list_data != None:
				if (self.list_idx > 0):
					self.list_idx -= 1
					self.data = self.list_data[self.list_idx]
					self.force_display_update()
					self.updateGL()
			#window = EMParentWin(browser)
			##browser.setParent(
			#window.show()
			#window.resize(800,800)
			

	def leaveEvent(self):
		if self.rmousedrag:
			self.rmousedrag=None

	def mouseReleaseEvent(self, event):
		app =  QtGui.QApplication.instance()
		app.setOverrideCursor(Qt.ArrowCursor)
		lc=self.scr_to_img(event.x(),event.y())
		if self.rmousedrag:
			self.rmousedrag=None
		elif event.button()==Qt.LeftButton:
			if self.mmode==0:
				self.parent.emit(QtCore.SIGNAL("mouseup"), event)
				return
			elif self.mmode==1 :
				self.add_shape("MEAS",EMShape(("line",.5,.1,.5,self.shapes["MEAS"].shape[4],self.shapes["MEAS"].shape[5],lc[0],lc[1],2)))
			elif self.mmode==2 and self.inspector:
				self.set_data(self.data)

	def wheelEvent(self, event):
		if not self.wheel_navigate:
			if event.orientation() & Qt.Vertical:
				if self.mmode==0 and event.modifiers()&Qt.ShiftModifier:
					self.parent.emit(QtCore.SIGNAL("mousewheel"), event)
					return
				if event.delta() > 0:
					self.set_scale( self.scale * self.mag )
				elif event.delta() < 0:
					self.set_scale(self.scale * self.invmag )
				# The self.scale variable is updated now, so just update with that
				if self.inspector: self.inspector.set_scale(self.scale)
		else:
			move_fac = 1.0/20.0
			delta = event.delta()/120.0
			
#			print self.origin, self.data.get_xsize(),self.data.get_ysize(),self.scale,self.parent.width(),self.parent.height()

#			print self.origin
			if event.orientation() & Qt.Vertical:
				visible_vertical_pixels = self.parent.height()/sqrt(self.scale)
				shift_per_delta = move_fac*visible_vertical_pixels
#				print "there are this many visible vertical pixels",visible_vertical_pixels, "deltas", delta, "shift per delta",shift_per_delta
#				print "shifting vertical",event.delta(),shift_per_delta
				self.origin=(self.origin[0],self.origin[1]-delta*shift_per_delta)
			elif event.orientation() & Qt.Horizontal:
				visible_horizontal_pixels = self.parent.width()/sqrt(self.scale)
				shift_per_delta = move_fac*visible_horizontal_pixels
#				print "shifting horizontal",event.delta(),shift_per_delta
#	   	   	   	print "there are this many visible horizontal pixels",visible_horizontal_pixels, "deltas", delta, "shift per delta",shift_per_delta
				self.origin=(self.origin[0]+delta*shift_per_delta,self.origin[1])
			try: self.parent.update()
			except: pass
#			print "exit",self.origin
			
	
	def mouseDoubleClickEvent(self,event):
		if platform.system() == "Darwin":
			self.wheel_navigate = not self.wheel_navigate
		else:
			print "double click only performs a function on Mac"
		
		
class EMImageInspector2D(QtGui.QWidget):
	def __init__(self,target) :
		QtGui.QWidget.__init__(self,None)
		self.target=target
		
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
		
		self.mtshowlen= QtGui.QLabel("Length:")
		self.mtlay.addWidget(self.mtshowlen,2,0,1,2,Qt.AlignHCenter)
		
		self.mtshoworigin= QtGui.QLabel("Origin:")
		self.mtlay.addWidget(self.mtshoworigin,1,0,Qt.AlignHCenter)
		
		self.mtshowend= QtGui.QLabel("End:")
		self.mtlay.addWidget(self.mtshowend,1,1,Qt.AlignHCenter)
		
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
		self.vbl2.addWidget(self.invtog,0,0,1,2)
		
		# FFT Buttons
		self.fftg=QtGui.QButtonGroup()
		self.fftg.setExclusive(1)
		
		self.ffttog0 = QtGui.QPushButton("Real")
		self.ffttog0.setCheckable(1)
		self.ffttog0.setChecked(1)
		self.vbl2.addWidget(self.ffttog0,1,0)
		self.fftg.addButton(self.ffttog0,0)

		self.ffttog1 = QtGui.QPushButton("FFT")
		self.ffttog1.setCheckable(1)
		self.vbl2.addWidget(self.ffttog1,1,1)
		self.fftg.addButton(self.ffttog1,1)

		self.ffttog2 = QtGui.QPushButton("Amp")
		self.ffttog2.setCheckable(1)
		self.vbl2.addWidget(self.ffttog2,2,0)
		self.fftg.addButton(self.ffttog2,2)
		
		self.ffttog3 = QtGui.QPushButton("Pha")
		self.ffttog3.setCheckable(1)
		self.vbl2.addWidget(self.ffttog3,2,1)
		self.fftg.addButton(self.ffttog3,3)
	
		self.scale = ValSlider(self,(0.1,5.0),"Mag:")
		self.scale.setObjectName("scale")
		self.scale.setValue(1.0)
		self.vbl.addWidget(self.scale)
		
		self.mins = ValSlider(self,label="Min:")
		self.mins.setObjectName("mins")
		self.mins.setValue(self.target.get_minden())
		self.vbl.addWidget(self.mins)
		
		self.maxs = ValSlider(self,label="Max:")
		self.maxs.setObjectName("maxs")
		self.maxs.setValue(self.target.get_maxden())
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
		self.gammas.setValue(self.target.get_gamma())
		#self.gammas.setValue(1.0)
		self.vbl.addWidget(self.gammas)

		self.setWindowIcon(QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/eman.png"))

		self.lowlim=0
		self.highlim=1.0
		#self.update_min_max()
		#self.update_brightness_contrast()
		self.busy=0
		
		QtCore.QObject.connect(self.scale, QtCore.SIGNAL("valueChanged"), target.set_scale)
		QtCore.QObject.connect(self.mins, QtCore.SIGNAL("valueChanged"), self.new_min)
		QtCore.QObject.connect(self.maxs, QtCore.SIGNAL("valueChanged"), self.new_max)
		QtCore.QObject.connect(self.brts, QtCore.SIGNAL("valueChanged"), self.new_brt)
		QtCore.QObject.connect(self.conts, QtCore.SIGNAL("valueChanged"), self.new_cont)
		QtCore.QObject.connect(self.gammas, QtCore.SIGNAL("valueChanged"), self.new_gamma)
		QtCore.QObject.connect(self.invtog, QtCore.SIGNAL("toggled(bool)"), target.set_invert)
		QtCore.QObject.connect(self.fftg, QtCore.SIGNAL("buttonClicked(int)"), target.set_FFT)
		QtCore.QObject.connect(self.mmtab, QtCore.SIGNAL("currentChanged(int)"), target.set_mouse_mode)
#		QtCore.QObject.connect(self.mmode, QtCore.SIGNAL("buttonClicked(int)"), target.set_mouse_mode)

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
		self.target.set_density_min(val)
		self.update_brightness_contrast()
		self.busy=0
		
	def new_max(self,val):
		if self.busy : return
		self.busy=1
		self.target.set_density_max(val)
		self.update_brightness_contrast()
		self.busy=0
	
	def new_brt(self,val):
		if self.busy : return
		self.busy=1
		self.target.set_brightness(val)
		self.update_min_max()
		self.busy=0
		
	def new_cont(self,val):
		if self.busy : return
		self.busy=1
		self.target.set_contrast(val)
		self.update_min_max()
		self.busy=0
		
	def new_gamma(self,val):
		if self.busy : return
		self.busy=1
		self.target.set_gamma(val)
		self.busy=0

	def update_brightness_contrast(self):
		b=0.5*(self.mins.value+self.maxs.value-(self.lowlim+self.highlim))/((self.highlim-self.lowlim))
		c=(self.mins.value-self.maxs.value)/(2.0*(self.lowlim-self.highlim))
		brts = -b
		conts = 1.0-c
		self.brts.setValue(brts)
		self.conts.setValue(conts)
		self.target.set_brightness_contrast(brts,conts)
		
	def update_min_max(self):
		
		x0=((self.lowlim+self.highlim)/2.0-(self.highlim-self.lowlim)*(1.0-self.conts.value)-self.brts.value*(self.highlim-self.lowlim))
		x1=((self.lowlim+self.highlim)/2.0+(self.highlim-self.lowlim)*(1.0-self.conts.value)-self.brts.value*(self.highlim-self.lowlim))
		self.mins.setValue(x0)
		self.maxs.setValue(x1)
		self.target.set_density_range(x0,x1)
		
	def set_hist(self,hist,minden,maxden):
		self.hist.set_data(hist,minden,maxden)

	def set_limits(self,lowlim,highlim,curmin,curmax):
		if highlim<=lowlim : highlim=lowlim+.001
		#print "in set limits", self.conts.getValue(), self.conts.getValue()
		self.lowlim=lowlim
		self.highlim=highlim
		self.mins.setRange(lowlim,highlim)
		self.maxs.setRange(lowlim,highlim)
		self.mins.setValue(curmin)
		self.maxs.setValue(curmax)
		self.target.set_density_range(curmin,curmax)
		#print "leave set limits", self.conts.getValue(), self.conts.getValue()
		

# This is just for testing, of course
if __name__ == '__main__':
	em_app = EMStandAloneApplication()
	window = EMImage2DModule(application=em_app)
	
	if len(sys.argv)==1 : 
		window.set_data(test_image(size=(128,128)))
	else :
		a=EMData.read_images(sys.argv[1],[0])
		if len(a) == 1:
			a = a[0]
			data = []
			if a.get_ndim() == 3:
				for z in range(a.get_zsize()):
					image = a.get_clip(Region(0,0,z,a.get_xsize(),a.get_ysize(),1))
					data.append(image)
				a  = data
		window.set_data(a)
		
		window.set_file_name(sys.argv[1])
		
	em_app.show()
	em_app.execute()
	
