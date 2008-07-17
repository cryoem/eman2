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
import sys
import numpy
from emimageutil import ImgHistogram,EMParentWin
from weakref import WeakKeyDictionary
from pickle import dumps,loads
from PyQt4.QtGui import QImage
from PyQt4.QtCore import QTimer

from emglobjects import EMOpenGLFlagsAndTools

from emfloatingwidgets import EMGLRotaryWidget, EMGLView2D,EM3DWidget


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
		
		self.fov = 20
		self.aspect = 1.0
		self.zNear = 2000
		self.zFar = 10000
		
		self.animatables = []
		
		self.timer = QTimer()
		QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.timeout)
		self.timer.start(10)
		
	def setData(self,data):
		self.image_rotary.setData(data)
	
	def get_optimal_size(self):
		lr = self.image_rotary.rotary.get_suggested_lr_bt_nf()
		width = lr[1] - lr[0]
		height = lr[3] - lr[2]
		return [width+80,height+20]
	
	def timeout(self):
		
		if len(self.animatables) == 0: return
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
		glLightfv(GL_LIGHT0, GL_POSITION, [0.1,.1,1.,0.])
	
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
		GLU.gluPerspective(self.fov,self.aspect,self.zNear,self.zFar)
		#GL.glOrtho(0.0,width,0.0,height,-width,width)
		GL.glMatrixMode(GL.GL_MODELVIEW)
		GL.glLoadIdentity()
		
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
		
		self.rotary = EMGLRotaryWidget(self,-25,10,40,EMGLRotaryWidget.LEFT_ROTARY)
		self.rotary.set_mmode("mxrotary")
		self.widget = EM3DWidget(self,self.rotary)
		self.widget.set_draw_frame(False)
		
		self.image_file_name = None	# keeps track of the image file name (if any) - book keeping purposes only
		self.emdata_list_cache = None # all import emdata list cache, the object that stores emdata objects efficiently. Must be initialized via setData or set_image_file_name
		
		self.visible_mxs = 4	# the number of visible imagemxs in the rotary
		self.mx_rows = 8	# the number of rows in any given imagemx
		self.mx_cols = 8 # the number of columns in any given imagemx
		self.start_mx = 0 # the starting index for the currently visible set of imagemxs

	def context(self):
		# asking for the OpenGL context from the parent
		return self.parent.context()
	
	def emit(self,signal,event,integer=None):
		if integer != None:
			self.parent.emit(signal,event,integer)
		else:pass
				
	def update_rotary_position(self,inc):
		self.start_mx += inc
		num_per_view = self.mx_rows*self.mx_cols
		
		if inc > 0:
			end_changed = self.visible_mxs
			start_changed = self.visible_mxs-inc
			
			#print start_changed,end_changed
			#for idx in range(start_changed,end_changed):
			
		elif inc < 0:
			
			end_changed = -inc
			start_changed = 0
			
		else: return
		
		for idx in range(start_changed,end_changed):
			n = idx + self.start_mx
			start_idx = n*num_per_view
			d = []
			for i in range(start_idx,start_idx+num_per_view): d.append(self.emdata_list_cache[i])
				
			self.rotary[idx].setData(d)
			w = self.rotary[idx].get_drawable()
			#print start_idx%self.emdata_list_cache.get_max_idx()
			
			w.set_img_num_offset(start_idx%self.emdata_list_cache.get_max_idx())
		
	def set_mmode(self,mode):
		self.mmode = mode
		self.rotary.set_mmode(mode)

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
		
	def setData(self,data):
		if data == None or not isinstance(data,list) or len(data)==0:
			self.data = [] 
			return
		
		self.emdata_list_cache = EMDataListCache(data)
		self.__refresh_rotary()
	
	def set_image_file_name(self,name):
		#print "set image file name",name
		self.image_file_name = name
		self.emdata_list_cache = EMDataListCache(name)
		self.__refresh_rotary()
	
	def __refresh_rotary(self):
		self.rotary.clear_widgets()
		num_per_view = self.mx_rows*self.mx_cols
		for idx in range(self.start_mx,self.visible_mxs+self.start_mx):
			start_idx = idx*num_per_view
			d = []
			for i in range(start_idx,start_idx+num_per_view): d.append(self.emdata_list_cache[i])
			e = EMGLView2D(self,d)
			e.setWidth(self.mx_rows*d[0].get_xsize())
			e.setHeight(self.mx_cols*d[0].get_ysize())
			e.get_drawable().set_img_num_offset(start_idx%self.emdata_list_cache.get_max_idx())
			e.get_drawable().set_max_idx(self.emdata_list_cache.get_max_idx())
			self.rotary.add_widget(e)

	def updateGL(self):
		try: self.parent.updateGL()
		except: pass


	def render(self):
		
		glLoadIdentity()
		
		lr = self.rotary.get_suggested_lr_bt_nf()
		#print lr
		GL.glEnable(GL.GL_DEPTH_TEST)
		GL.glEnable(GL.GL_LIGHTING)
		#print self.parent.get_depth_for_height(self.height())
		#lr = self.rotary.get_lr_bt_nf()

		GL.glPushMatrix()
		#print -self.parent.get_depth_for_height(abs(lr[3]-lr[2]))
		glTranslate(-(lr[1]+lr[0])/2.0,-(lr[3]+lr[2])/2.0,-self.parent.get_depth_for_height(abs(lr[3]-lr[2])))
		self.widget.paintGL()
		GL.glPopMatrix()
	
	def dragEnterEvent(self,event):
		pass

	
	def dropEvent(self,event):
		pass


	def mousePressEvent(self, event):
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


class EMDataListCache:
	'''
	This class designed primarily for memory management in the context of large lists of EMData objects.
	
	The main public interface is to acquired a list of EMData objects in a specific range
	'''
	#LIST_MODE = 'list_mode'
	#FILE_MODE = 'file_mode'
	def __init__(self,object,cache_size=100,start_idx=0):
		if isinstance(object,list):
			# in list mode there is no real caching
			#self.mode = EMDataListCache.LIST_MODE
			self.max_idx = len(object)
			self.cache_size = self.max_idx
			self.images = object
			self.start_idx = 0
		elif isinstance(object,str):
			print "file mode"
			#self.mode = EMDataListCache.FILE_MODE
			if not os.path.exists(object):
				print "error, the file you specified does not exist:",object
				return
			self.file_name = object
			self.max_idx = EMUtil.get_image_count(self.file_name)
			self.images = {}
			self.start_idx = start_idx
			if self.max_idx < cache_size:
				self.cache_size = self.max_idx
			else:
				self.cache_size = cache_size
			self.__refresh_cache()
		else:
			print "the object used to construct the EMDataListCache is not a string (filename) or a list (of EMData objects). Can't proceed"
			return
	
	def get_max_idx(self):
		return self.max_idx
	
	def set_cache_size(self,cache_size):
		self.cache_size = cache_size
		if refresh: self.__refresh_cache()
		
	def set_start_idx(self,start_idx,refresh=True):
		self.start_idx = start_idx
		if refresh: self.__refresh_cache()
	
	def __refresh_cache(self):
		cache = {}
		for i in range(self.start_idx,self.start_idx+self.cache_size,1):
			if i != 0:
				idx = i % self.max_idx
			else: idx = 0
			try: 
				cache[idx] = self.images[idx]
			except:
				cache[idx] = EMData(self.file_name,idx)
			#print i,idx
				
		self.images = cache
	
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
				print "error, could get image",i,self.start_idx,self.max_idx
				for i in self.images:
					print i,
				print ''
			
# This is just for testing, of course
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	GLUT.glutInit("")
	window = EMImageMXRotary()
	if len(sys.argv)==1 : window.setData([test_image(),test_image(1),test_image(2),test_image(3),test_image(2),test_image(1),test_image()]*50)
	else :
		print "reading image"
		#a=EMData.read_images(sys.argv[1])
		window.set_image_file_name(sys.argv[1])
	window2=EMParentWin(window)
	window2.show()
	window2.resize(*window.get_optimal_size())
#	w2=QtGui.QWidget()
#	w2.resize(256,128)
	
#	w3=ValSlider(w2)
#	w3.resize(256,24)
#	w2.show()
	
	sys.exit(app.exec_())
