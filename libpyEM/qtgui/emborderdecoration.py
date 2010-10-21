#!/usr/bin/env python

#
# Author: David Woolford 2008 (woolford@bcm.edu)
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
from OpenGL.GL import *
from OpenGL.GLU import *
import os

import warnings
warnings.warn("emborderdecoration.py is deprecated.", DeprecationWarning)


white = (1.0,1.0,1.0,1.0)
grey = (0.5,0.5,0.5,1.0)
yellow = (1.0,1.0,0.0,1.0)
red = (1.0,1.0,1.0,1.0)
blue = (0.0,0.0,1.0,1.0)
black = (0.2,0.2,0.2,1.0)
grey = (0.6,0.6,0.6,0.0)
green = (0.0,1.0,0.0,1.0)
purple = (1.0,0.0,1.0,1.0)
black = (0.0,0.0,0.0,1.0)
dark_purple = (0.4,0.0,0.4,1.0)
light_blue_diffuse = (.4,.53,.86,1.0)
light_blue_ambient = (.4,.7,.89,1.0)

ligh_yellow_diffuse = (.84,.82,.38,1.0)
ligh_yellow_ambient = (.83,.83,.38,1.0)
ligh_yellow_specular = (.76,.75,.39,1.0)

from EMAN2 import Vec3f, get_image_directory
import weakref


class EMBorderDecoration:
	'''
	An class for drawing borders around widgets
	The implementation of border decorations in EMAN2 floating widgets is based on the 
	Decorator pattern in the Gang of Four. Inheriting classes should provide draw(object)
	which draws a border around the given object.
	'''
	FROZEN_COLOR = "frozen"
	DEFAULT_COLOR = "default"
	YELLOW_COLOR = "yellow"
	BLACK_COLOR = "black"
	x_texture_dl = None
	PERMISSABLE_COLOR_FLAGS = [FROZEN_COLOR,DEFAULT_COLOR,YELLOW_COLOR,BLACK_COLOR]
	def __init__(self,object):
		self.object = weakref.ref(object)
		self.color_flag = EMBorderDecoration.DEFAULT_COLOR
		self.display_list = None
		self.unprojected_points = [] # for determing mouse-border collision detection
		self.pressed = -1
		self.current_frame = ""
		self.moving = [-1,-1]
		self.init_frame_unproject_order()
		self.corner_sets = []
		self.bottom_border_height = 7
		self.top_border_height = 14
		self.border_width = 7
		self.border_depth = 6
		
		self.draw_x_enabled = True
		
		try:
			from EMAN2 import get_3d_font_renderer,FTGLFontMode
			self.font_renderer = get_3d_font_renderer()
			self.font_renderer.set_face_size(self.top_border_height-2)
			self.font_renderer.set_font_mode(FTGLFontMode.TEXTURE)
		except:
			self.font_renderer = None
			
		self.window_title = ""
		
		self.init_x_texture()
		
		self.do_clip = False
		
		#gl_context_parent = object.get_connection_object()
		#QtCore.QObject.connect(gl_context_parent, QtCore.SIGNAL("window_selected_added"), gl_context_parent.window_selected_added)
		
	def enable_clip(self,val=True):
		self.do_clip = val
		
	def __del__(self):
		self.delete_list()
	
	def draw(self,object): raise
	
	def total_width(self):
		return 2*self.border_width
	
	def setWindowTitle(self,title):
		self.window_title = title
	
	def total_height(self):
		return self.top_border_height+self.bottom_border_height

	def delete_list(self):
		if self.display_list != None:
			glDeleteLists(self.display_list,1)
			self.display_list = None

	def frame_inner_shell_basic(self,width,height,front,back):
		w = width/2.0
		h = height/2.0
		glBegin(GL_TRIANGLE_STRIP)
		glNormal(-1,0,0)
		glVertex(w,-h,back)
		glVertex(w,-h,front)
		glVertex(w,h,back)
		glVertex(w,h,front)
		glEnd()
		
		glNormal(0,-1,0)
		glBegin(GL_TRIANGLE_STRIP)
		glVertex(w,h,back)
		glVertex(w,h,front)
		glVertex(-w,h,back)
		glVertex(-w,h,front)
		glEnd()
		
		glNormal(1,0,0)
		glBegin(GL_TRIANGLE_STRIP)
		glVertex(-w,h,back)
		glVertex(-w,h,front)
		glVertex(-w,-h,back)
		glVertex(-w,-h,front)
		glEnd()
		
		glNormal(0,1,0)
		glBegin(GL_TRIANGLE_STRIP)
		glVertex(-w,-h,back)
		glVertex(-w,-h,front)
		glVertex(w,-h,back)
		glVertex(w,-h,front)
		glEnd()	
	
	def frame_outer_shell_basic(self,width_plus,height_plus,front,back):
		
		wp = width_plus/2.0
		hp = height_plus/2.0
		
		glNormal(-1,0,0)
		glBegin(GL_TRIANGLE_STRIP)
		glVertex(-wp,-hp,back)
		glVertex(-wp,-hp,front)
		glVertex(-wp,hp,back)
		glVertex(-wp,hp,front)
		glEnd()
		
		glNormal(0,1,0)
		glBegin(GL_TRIANGLE_STRIP)
		glVertex(-wp,hp,back)
		glVertex(-wp,hp,front)
		glVertex(wp,hp,back)
		glVertex(wp,hp,front)
		glEnd()
		
		glNormal(1,0,0)
		glBegin(GL_TRIANGLE_STRIP)
		glVertex(wp,hp,back)
		glVertex(wp,hp,front)
		glVertex(wp,-hp,back)
		glVertex(wp,-hp,front)
		glEnd()
		
		glNormal(0,-1,0)
		glBegin(GL_TRIANGLE_STRIP)
		glVertex(wp,-hp,back)
		glVertex(wp,-hp,front)
		glVertex(-wp,-hp,back)
		glVertex(-wp,-hp,front)
		glEnd()
		
	def frame_face_basic(self,width,width_plus,height,height_plus):
		w = width/2.0
		wp = width_plus/2.0
		h = height/2.0
		hp = height_plus/2.0
		
		glNormal(0,0,1)
		glBegin(GL_TRIANGLE_STRIP)
		glVertex(-wp,-hp,0)
		glVertex(-w,-hp,0)
		glVertex(-wp,hp,0)
		glVertex(-w,hp,0)
		glEnd()
		
		glBegin(GL_TRIANGLE_STRIP)
		glVertex(-wp,hp,0)
		glVertex(-w,h,0)
		glVertex(w,hp,0)
		glVertex(w,h,0)
		glEnd()
		
		glBegin(GL_TRIANGLE_STRIP)
		glVertex(wp,hp,0)
		glVertex(w,hp,0)
		glVertex(wp,-hp,0)
		glVertex(w,-hp,0)
		glEnd()
		
		glBegin(GL_TRIANGLE_STRIP)
		glVertex(w,-h,0)
		glVertex(w,-hp,0)
		glVertex(-w,-h,0)
		glVertex(-w,-hp,0)
		glEnd()
		
	
	def frame_face(self,left,left_plus,right,right_plus,bottom,bottom_plus,top,top_plus):
		glNormal(0,0,1)
		glBegin(GL_TRIANGLE_STRIP)
		glVertex(left_plus,bottom_plus,0)
		glVertex(left,bottom_plus,0)
		glVertex(left_plus,top_plus,0)
		glVertex(left,top_plus,0)
		glEnd()
		
		glBegin(GL_TRIANGLE_STRIP)
		glVertex(left_plus,top_plus,0)
		glVertex(left,top,0)
		glVertex(right,top_plus,0)
		glVertex(right,top,0)
		glEnd()
		
		glBegin(GL_TRIANGLE_STRIP)
		glVertex(right_plus,top_plus,0)
		glVertex(right,top_plus,0)
		glVertex(right_plus,bottom_plus,0)
		glVertex(right,bottom_plus,0)
		glEnd()
		
		glBegin(GL_TRIANGLE_STRIP)
		glVertex(right,bottom,0)
		glVertex(right,bottom_plus,0)
		glVertex(left,bottom,0)
		glVertex(left,bottom_plus,0)
		glEnd()
	
	def frame_outer_shell(self,left_plus,right_plus,bottom_plus,top_plus,front,back):
		
		glNormal(-1,0,0)
		glBegin(GL_TRIANGLE_STRIP)
		glVertex(left_plus,bottom_plus,back)
		glVertex(left_plus,bottom_plus,front)
		glVertex(left_plus,top_plus,back)
		glVertex(left_plus,top_plus,front)
		glEnd()
		
		glNormal(0,1,0)
		glBegin(GL_TRIANGLE_STRIP)
		glVertex(left_plus,top_plus,back)
		glVertex(left_plus,top_plus,front)
		glVertex(right_plus,top_plus,back)
		glVertex(right_plus,top_plus,front)
		glEnd()
		
		glNormal(1,0,0)
		glBegin(GL_TRIANGLE_STRIP)
		glVertex(right_plus,top_plus,back)
		glVertex(right_plus,top_plus,front)
		glVertex(right_plus,bottom_plus,back)
		glVertex(right_plus,bottom_plus,front)
		glEnd()
		
		glNormal(0,-1,0)
		glBegin(GL_TRIANGLE_STRIP)
		glVertex(right_plus,bottom_plus,back)
		glVertex(right_plus,bottom_plus,front)
		glVertex(left_plus,bottom_plus,back)
		glVertex(left_plus,bottom_plus,front)
		glEnd()
		
	def frame_inner_shell(self,left,right,bottom,top,front,back):
		glBegin(GL_TRIANGLE_STRIP)
		glNormal(-1,0,0)
		glVertex(right,bottom,back)
		glVertex(right,bottom,front)
		glVertex(right,top,back)
		glVertex(right,top,front)
		glEnd()
		
		glNormal(0,-1,0)
		glBegin(GL_TRIANGLE_STRIP)
		glVertex(right,top,back)
		glVertex(right,top,front)
		glVertex(left,top,back)
		glVertex(left,top,front)
		glEnd()
		
		glNormal(1,0,0)
		glBegin(GL_TRIANGLE_STRIP)
		glVertex(left,top,back)
		glVertex(left,top,front)
		glVertex(left,bottom,back)
		glVertex(left,bottom,front)
		glEnd()
		
		glNormal(0,1,0)
		glBegin(GL_TRIANGLE_STRIP)
		glVertex(left,bottom,back)
		glVertex(left,bottom,front)
		glVertex(right,bottom,back)
		glVertex(right,bottom,front)
		glEnd()
		
	def set_color_flag(self,flag):
		if flag not in EMBorderDecoration.PERMISSABLE_COLOR_FLAGS:
			print 'unknown color flag'
		else:
			self.color_flag = flag
			
	def load_materials(self):
		glMaterial(GL_FRONT,GL_EMISSION,(0,0,0,1))
		if self.color_flag ==  EMBorderDecoration.DEFAULT_COLOR:
			glMaterial(GL_FRONT,GL_AMBIENT,light_blue_ambient)
			glMaterial(GL_FRONT,GL_DIFFUSE,light_blue_diffuse)
			glMaterial(GL_FRONT,GL_SPECULAR,light_blue_diffuse)
			glMaterial(GL_FRONT,GL_SHININESS,20.0)
			glColor(*white)
		elif self.color_flag ==  EMBorderDecoration.FROZEN_COLOR:
			glMaterial(GL_FRONT,GL_AMBIENT,dark_purple)
			glMaterial(GL_FRONT,GL_DIFFUSE,blue)
			glMaterial(GL_FRONT,GL_SPECULAR,light_blue_ambient)
			glMaterial(GL_FRONT,GL_SHININESS,20.0)
			glColor(*blue)
		elif self.color_flag ==  EMBorderDecoration.YELLOW_COLOR:
			glMaterial(GL_FRONT,GL_AMBIENT,ligh_yellow_ambient)
			glMaterial(GL_FRONT,GL_DIFFUSE,ligh_yellow_diffuse)
			glMaterial(GL_FRONT,GL_SPECULAR,ligh_yellow_specular)
			glMaterial(GL_FRONT,GL_SHININESS,18.0)
			glColor(*blue)
		elif self.color_flag ==  EMBorderDecoration.BLACK_COLOR:
			glMaterial(GL_FRONT,GL_AMBIENT,black)
			glMaterial(GL_FRONT,GL_DIFFUSE,black)
			glMaterial(GL_FRONT,GL_SPECULAR,white)
			glMaterial(GL_FRONT,GL_SHININESS,100.0)
			glColor(*blue)
		else:
			print "warning, unknown color flag, coloring failed"

	
	def is_selected(self):
		return  self.color_flag == EMBorderDecoration.YELLOW_COLOR
	
	def set_selected(self,bool=True):
		if bool: self.color_flag = EMBorderDecoration.YELLOW_COLOR
		else: self.color_flag =  EMBorderDecoration.DEFAULT_COLOR
			
	def init_frame_unproject_order(self):
		self.frame_unproject_order = []
		self.currently_selected_frame = -1
		self.frame_unproject_order.append("bottom_left")
		self.frame_unproject_order.append("bottom")
		self.frame_unproject_order.append("bottom_right")
		self.frame_unproject_order.append("right")
		self.frame_unproject_order.append("top_right")
		self.frame_unproject_order.append("top")
		self.frame_unproject_order.append("top_left")
		self.frame_unproject_order.append("left")
		
	def init_x_texture(self):
		if EMBorderDecoration.x_texture_dl == None or EMBorderDecoration.x_texture_dl < 0:
			#glEnable(GL_TEXTURE_2D)
			pixmap = QtGui.QPixmap(get_image_directory() + "/Close.png")
			#self.qwidget.setVisible(False)
			if (pixmap.isNull() == True ): print 'error, the pixmap was null'
			self.texture = self.object().parent().get_gl_context_parent().bindTexture(pixmap)
			EMBorderDecoration.x_texture_dl=glGenLists(1)
			
			if EMBorderDecoration.x_texture_dl == 0:
				print "BORDER GENERATION FAILED"
				import sys
				sys.exit(1)
				return
			glNewList( EMBorderDecoration.x_texture_dl,GL_COMPILE)
			glEnable(GL_TEXTURE_2D)
			
			texture = glGenTextures(1)
			glPixelStorei(GL_UNPACK_ALIGNMENT,1)
			#print glGetInteger(GL_UNPACK_ALIGNMENT)
			glBindTexture(GL_TEXTURE_2D, self.texture)
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)
			# this makes it so that the texture is impervious to lighting
			#glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE)
			
			xmin = -0.5
			xmax =  0.5
			ymin = -0.5
			ymax =  0.5
			glBegin(GL_QUADS)
		
			glNormal(0,0,1)
			glTexCoord2f(0,0)
			glVertex2f(xmin,ymax)
			
			glTexCoord2f(1,0)
			glVertex2f(xmax,ymax)
				
			glTexCoord2f(1,1)
			glVertex2f(xmax,ymin)
			
			glTexCoord2f(0,1)
			glVertex2f(xmin,ymin)
			
			glEnd()
			glDisable(GL_TEXTURE_2D)
			
			glEndList()
			
			
			#EM3DPlainBorderDecoration.x_texture_dl

	def mousePressEvent(self,event):
		self.pressed = self.currently_selected_frame 
		if self.currently_selected_frame == 5: # top
			self.object().set_mouse_lock(self)
			self.object().cam.mousePressEvent(event)
			self.object().updateGL()
			self.current_frame = self.frame_unproject_order[self.currently_selected_frame]
		elif self.currently_selected_frame in [0,1,2,3,4,6,7]:  # right
			if (event.buttons()&Qt.LeftButton):
				self.object().set_mouse_lock(self)
				self.current_frame = self.frame_unproject_order[self.currently_selected_frame]
				self.moving = [event.x(), event.y()]
			
		
	def mouseMoveEvent(self,event):
		if self.current_frame != "": # establishing and destablished in mousePress and mouseRelease
			if self.currently_selected_frame == 5: # top
				self.object().cam.mouseMoveEvent(event)
				self.object().updateGL()
			if (event.buttons()&Qt.LeftButton):
				movement = [event.x()-self.moving[0],event.y()-self.moving[1]]
				self.moving = [event.x(), event.y()]
				if self.current_frame == "right":
					self.object().add_width_right(movement[0])
					self.object().updateGL()
				if self.current_frame == "left":
					self.object().add_width_left(-movement[0])
					self.object().updateGL()
				if self.current_frame == "bottom":
					self.object().add_height_bottom(movement[1])
					self.object().updateGL()
				if self.current_frame == "bottom_left":
					self.object().add_width_left(-movement[0])
					self.object().add_height_bottom(movement[1])
					self.object().add_depth((movement[1]+movement[0])/2)
					self.object().updateGL()
				if self.current_frame == "bottom_right":
					self.object().add_width_right(movement[0])
					self.object().add_height_bottom(movement[1])
					self.object().add_depth((movement[1]+movement[0])/2)
					self.object().updateGL()
				if self.current_frame == "top_right":
					self.object().add_width_right(movement[0])
					self.object().add_height_top(-movement[1])
					self.object().add_depth((movement[1]+movement[0])/2)
					self.object().updateGL()
				if self.current_frame == "top_left":
					self.object().add_width_left(-movement[0])
					self.object().add_height_top(-movement[1])
					self.object().add_depth((movement[1]+movement[0])/2)
					self.object().updateGL()
			
	def mouseReleaseEvent(self,event):
		if self.currently_selected_frame == 5: # top
			self.object().release_mouse_lock()
			self.object().cam.mouseReleaseEvent(event)
			self.object().updateGL()
		elif self.currently_selected_frame == 8 and self.pressed == 8:
			self.object().closeEvent(event)
		
		if self.current_frame != "":
			self.object().release_mouse_lock()
			self.current_frame = ""
			self.object().updateGL()
	
	def mouseDoubleClickEvent(self,event):
		pass
	
	def wheelEvent(self,event):
		if self.currently_selected_frame == 5: # top
			self.object().top_frame_wheel(event.delta())
			self.object().updateGL()
			
	def keyPressEvent(self,event):
		pass
	

class EM2DPlainBorderDecoration(EMBorderDecoration):
	def __init__(self,object,gl_context_parent): 
		from emglobjects import EMViewportDepthTools2
		from emfloatingwidgets import EM2DGLWindow
		EMBorderDecoration.__init__(self,object)
		self.gl_context_parent = gl_context_parent
		self.vdtools = EMViewportDepthTools2(gl_context_parent)
		
		
		self.force_update = False
		
		self.faulty = False
		if not isinstance(object,EM2DGLWindow) and not isinstance(object,EMGLViewQtWidget): 
			print "error, EM2DPlainBorderDecoration works only for EM2DGLWindow and EMGLViewQtWidget"
			self.faulty = True
			return
		else: self.object = weakref.ref(object)
		
	def draw(self,force_update=False):
		#self.init_x_texture()
		if force_update or self.force_update:
			self.delete_list()
			self.force_update = False
		
		if self.display_list == None:
			self.__gen_2d_object_border_list()
		
		if self.display_list == None:
			print "display list generation failed" 
			return
		
		self.viewport_update()
		
		self.load_materials()
			
		if self.display_list != None and self.display_list != 0:
			glCallList(self.display_list)
		else: print "An error occured"
	
		glPushMatrix()
		glTranslate(self.object().width()/2.0-self.top_border_height/2,self.object().height()/2.0+self.top_border_height/2.0,self.border_depth+.1)
		glScale(self.top_border_height,self.top_border_height,1.0)
		#glEnable(GL_BLEND)
		#glBlendFunc(GL_DST_ALPHA,GL_ONE_MINUS_DST_ALPHA)
		glCallList( EMBorderDecoration.x_texture_dl)
		#glDisable(GL_BLEND)
		glPopMatrix()
		
		if self.font_renderer != None:
			glPushMatrix()
			glDisable(GL_LIGHTING)
			glEnable(GL_TEXTURE_2D)
			bbox = self.font_renderer.bounding_box(self.window_title)
			glTranslate(4,4,0.2)

			glTranslate((bbox[0]-bbox[3])/2,0,0)
			glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE)
			glTranslate(0,self.object().height()/2.0,self.border_depth+.1)
			self.font_renderer.render_string(self.window_title)
			glPopMatrix()
		
	def viewport_update(self):
		self.corner_sets = []
		width = self.object().width()
		height = self.object().height()
			
			# plus, i.e. plus the border
		left = -width/2.0
		left_plus = left-self.border_width
		right = -left
		right_plus = -left_plus
			
		bottom = -height/2.0
		bottom_plus = bottom - self.bottom_border_height
		top = -bottom
		top_plus = top + self.top_border_height
			#top_plus = -bottom_plus
			
		front = self.border_depth/2.0
		back = -self.border_depth/2.0
		
		points = []
		points.append((left_plus,bottom_plus,0))
		points.append((left,bottom_plus,0))
		points.append((left,bottom,0))
		points.append((left_plus,bottom,0))
		
		points.append((right,bottom_plus,0))
		points.append((right_plus,bottom_plus,0))
		points.append((right_plus,bottom,0))
		points.append((right,bottom,0))
		
		points.append((right,top,0))
		points.append((right_plus,top,0))
		points.append((right_plus,top_plus,0))
		points.append((right,top_plus,0))
		
		points.append((left_plus,top,0))
		points.append((left,top,0))
		points.append((left,top_plus,0))
		points.append((left_plus,top_plus,0))	
		
		points.append((right-self.top_border_height,top,0))
		points.append((right-self.top_border_height,top_plus,0))
		
		self.unprojected_points = self.vdtools.unproject_points(points)
		unprojected = self.unprojected_points
		self.vdtools.set_mouse_coords(unprojected[0],unprojected[1],unprojected[2],unprojected[3])
		self.corner_sets.append( self.vdtools.get_corners() )
		
		self.vdtools.set_mouse_coords(unprojected[1],unprojected[4],unprojected[7],unprojected[2])
		self.corner_sets.append( self.vdtools.get_corners() )
		
		self.vdtools.set_mouse_coords(unprojected[4],unprojected[5],unprojected[6],unprojected[7])
		self.corner_sets.append( self.vdtools.get_corners() )
		
		self.vdtools.set_mouse_coords(unprojected[7],unprojected[6],unprojected[9],unprojected[8])
		self.corner_sets.append( self.vdtools.get_corners() )
		
		self.vdtools.set_mouse_coords(unprojected[8],unprojected[9],unprojected[10],unprojected[11])
		self.corner_sets.append( self.vdtools.get_corners() )
		
		self.vdtools.set_mouse_coords(unprojected[13],unprojected[16],unprojected[17],unprojected[14])
		self.corner_sets.append( self.vdtools.get_corners() )
		
		self.vdtools.set_mouse_coords(unprojected[12],unprojected[13],unprojected[14],unprojected[15])
		self.corner_sets.append( self.vdtools.get_corners() )
		
		self.vdtools.set_mouse_coords(unprojected[3],unprojected[2],unprojected[13],unprojected[12])
		self.corner_sets.append( self.vdtools.get_corners() )
		
		self.vdtools.set_mouse_coords(unprojected[16],unprojected[8],unprojected[11],unprojected[17])
		self.corner_sets.append( self.vdtools.get_corners() )
	
	def isinwin(self,x,y):
		interception = False
		#print "I have this many corner sets",len(self.corner_sets)
		for i,p in enumerate(self.corner_sets):
			if self.vdtools.isinwinpoints(x,y,p):
				interception = True
				self.currently_selected_frame = i
				break
		return interception
		
	def get_border_width(self): return self.border_width
	#def get_border_height(self): return self.border_height
	
	def __gen_2d_object_border_list(self):
		
		#context = self.object.context()
		#context.makeCurrent()
		#print "made",context,"current"
		
		#if EM3DPlainBorderDecoration.x_texture_dl == None:
			#self.__init_x_texture()
		
		self.delete_list()
		if self.display_list == None:
			
			width = self.object().width()
			height = self.object().height()
			
			# plus, i.e. plus the border
			left = -width/2.0
			left_plus = left-self.border_width
			right = -left
			right_plus = -left_plus
			
			bottom = -height/2.0
			bottom_plus = bottom - self.bottom_border_height
			top = -bottom
			top_plus = top + self.top_border_height
			#top_plus = -bottom_plus
			
			front = self.border_depth/2.0
			back = -self.border_depth/2.0
			
			self.display_list=glGenLists(1)
			
			if self.display_list == 0:
				return
			
			glNewList(self.display_list,GL_COMPILE)

	
			# All triangles are drawn in counter clockwise direction
			# Do the front facing strip first
			glPushMatrix()
			glTranslate(0,0,front)
			self.frame_face(left,left_plus,right,right_plus,bottom,bottom_plus,top,top_plus)
			glPopMatrix()
			
			# Do the back facing part
			glPushMatrix()
			glTranslate(0,0,back)
			glRotate(180,0,1,0)
			self.frame_face(left,left_plus,right,right_plus,bottom,bottom_plus,top,top_plus)
			glPopMatrix()
			
			# Now do the border around the edges
			glPushMatrix()
			self.frame_outer_shell(left_plus,right_plus,bottom_plus,top_plus,front,back)
			glPopMatrix()
			
			# Now do the border around the inside edges
			glPushMatrix()
			self.frame_inner_shell(left,right,bottom,top,front,back)
			glPopMatrix()
			
			glEndList()
		else: print "error, the delete list operation failed"
	
		
class EM3DPlainBorderDecoration(EMBorderDecoration):
	'''
	A plain border decoration
	'''
	
	
	
	#x_texture_size = (-1,1)
	def __init__(self, object,gl_context_parent):
		from emglobjects import EMViewportDepthTools2
		from emfloatingwidgets import EM3DGLWindow
		EMBorderDecoration.__init__(self,object)
		self.gl_context_parent = gl_context_parent
		self.vdtools = EMViewportDepthTools2(gl_context_parent)
		
		self.force_update = False
		
		self.faulty = False
		#if not isinstance(object,EM3DGLWindow):
			#print "error, border construction"
			#self.faulty = True
			#return
		#else:
		self.object = weakref.ref(object)
		
	
	def get_border_width(self): return self.border_width
	#def get_border_height(self): return self.border_height
	def get_border_depth(self): return self.border_depth
	
	def set_force_update(self,val=True):
		self.force_update = val
	
	
	
	
	def enable_clip_planes(self,two_d_only=False):
		self.calculate_clip_planes(two_d_only)
		glEnable(GL_CLIP_PLANE0)
		glEnable(GL_CLIP_PLANE1)
		glEnable(GL_CLIP_PLANE2)
		glEnable(GL_CLIP_PLANE3)
		if not two_d_only:
			glEnable(GL_CLIP_PLANE4)
			glEnable(GL_CLIP_PLANE5)
	
	def disable_clip_planes(self,two_d_only=False):
		glDisable(GL_CLIP_PLANE0)
		glDisable(GL_CLIP_PLANE1)
		glDisable(GL_CLIP_PLANE2)
		glDisable(GL_CLIP_PLANE3)
		if not two_d_only:
			glDisable(GL_CLIP_PLANE4)
			glDisable(GL_CLIP_PLANE5)

	def draw(self,force_update=False):
		
			
			
		if force_update or self.force_update:
			self.delete_list()
			
			self.force_update = False
		
		if self.display_list == None:
			self.__gen_3d_object_border_list()
			
	
		if self.display_list == None:
			print "display list generation failed" 
			return
		
		
		self.viewport_update()
		
		self.load_materials()
			
		if self.display_list != None and self.display_list != 0:
			glCallList(self.display_list)
		else: print "An error occured"
	
	
		if self.font_renderer != None:
			dims = self.object().get_lr_bt_nf()
			glPushMatrix()
			light_on = glIsEnabled(GL_LIGHTING)
			texture_on = glIsEnabled(GL_TEXTURE_2D)
			glDisable(GL_LIGHTING)
			glEnable(GL_TEXTURE_2D)
			bbox = self.font_renderer.bounding_box(self.window_title)
			glTranslate((dims[0]+dims[1])/2 , (dims[3]+dims[2])/2+ 2, dims[4]+self.border_depth)

			glTranslate((bbox[0]-bbox[3])/2,0,0)
			glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE)
			glTranslate(0,self.object().height()/2.0,self.border_depth+.2)
			self.font_renderer.render_string(self.window_title)
			glPopMatrix()
			
			if not texture_on: glDisable(GL_TEXTURE_2D)
			if light_on: glEnable(GL_LIGHTING)
	
		if self.draw_x_enabled:
			self.draw_x()
		
	def draw_x(self):
		dims = self.object().get_lr_bt_nf()
		# plus, i.e. plus the border
		
		right = dims[1]
		top = dims[3]
		width =  dims[1]- dims[0]
		height = dims[3] -  dims[2]
		front =  dims[4]
		glPushMatrix()
		glTranslate(right-self.top_border_height/2,top+self.top_border_height/2.0,front+self.border_depth+.1)
		glScale(self.top_border_height,self.top_border_height,1.0)
		glCallList( EMBorderDecoration.x_texture_dl)
		glPopMatrix()
		
	def calculate_clip_planes(self,two_d_only=False):
		dims = self.object().get_lr_bt_nf()
		if two_d_only:
			front = 0
			back = -1
		else:
			front =  dims[4]
			back =  dims[5]
		left =  dims[0]
		right =  dims[1]
		
		bottom = dims[2]
		top =  dims[3]
		
		# left plane
		
		p1 = Vec3f(left,bottom,front)
		p2 =  Vec3f(left,bottom,back)
		p3 = Vec3f(left,top,front)
		
		v1 = p3 - p1
		v2 = p2 - p1
		
		cross = v2.cross(v1)
		cross.normalize()
		d = -left
		plane = [cross[0],cross[1],cross[2],d]
		#print plane
		
		glClipPlane(GL_CLIP_PLANE0,plane)
		
		# right plane
		
		p1 = Vec3f(right,bottom,front)
		p2 =  Vec3f(right,bottom,back)
		p3 = Vec3f(right,top,front)
		
		v1 = p3 - p1
		v2 = p2 - p1
		
		cross = v1.cross(v2)
		cross.normalize()
		d = right
		plane2 = [cross[0],cross[1],cross[2],d]
		#print plane
		
		glClipPlane(GL_CLIP_PLANE1,plane2)
		
		
		# top
		p1 = Vec3f(right,top,front)
		p2 = Vec3f(right,top,back)
		p3 = Vec3f(left,top,front)
		
		v1 = p3-p1
		v2 = p2-p1
		
		cross = v1.cross(v2)
		cross.normalize()
		d = top
		plane3 = [cross[0],cross[1],cross[2],d]
		#print plane
		
		glClipPlane(GL_CLIP_PLANE2,plane3)
		
		
		# bottom
		p1 = Vec3f(right,bottom,front)
		p2 = Vec3f(right,bottom,back)
		p3 = Vec3f(left,bottom,front)
		
		v1 = p3-p1
		v2 = p2-p1
		
		cross = v2.cross(v1)
		cross.normalize()
		d = bottom
		plane4 = [cross[0],cross[1],cross[2],-d]
		#print plane
		
		glClipPlane(GL_CLIP_PLANE3,plane4)
		
		if not two_d_only:
			#front
			p1 = Vec3f(right,bottom,front)
			p2 = Vec3f(right,top,front)
			p3 = Vec3f(left,bottom,front)
			
			v1 = p3-p1
			v2 = p2-p1
			
			cross = v1.cross(v2)
			cross.normalize()
			d = front
			plane5 = [cross[0],cross[1],cross[2],d]
			#print plane
			
			glClipPlane(GL_CLIP_PLANE4,plane5)
			
			#back
			p1 = Vec3f(right,bottom,back)
			p2 = Vec3f(right,top,back)
			p3 = Vec3f(left,bottom,back)
			
			v1 = p3-p1
			v2 = p2-p1
			
			cross = v2.cross(v1)
			cross.normalize()
			d = back
			plane6 = [cross[0],cross[1],cross[2],-d]
			#print plane
			
			glClipPlane(GL_CLIP_PLANE5,plane6)
			
		
	def __gen_3d_object_border_list(self):
		self.delete_list()
	
		#if EM3DPlainBorderDecoration.x_texture_dl == None:
			#self.__init_x_texture()
	
		if self.display_list == None:
			thick_front = self.border_depth
			thick_back = -self.border_depth
				
			dims = self.object().get_lr_bt_nf()
			# plus, i.e. plus the border
			left =  dims[0]
			left_plus = left-self.border_width
			right = dims[1]
			right_plus = right + self.border_width
			
			bottom = dims[2]
			bottom_plus = bottom - self.bottom_border_height
			top =  dims[3]
			top_plus = top + self.top_border_height
			
			
			front =  dims[4]
			front_plus = front + self.border_depth
			back = dims[5]
			back_plus = back - self.border_depth


			x_center = (right+left)/2.0
			y_center = (top+bottom)/2.0
			z_center = (front+back)/2.0
			self.display_list=glGenLists(1)
				
			glNewList(self.display_list,GL_COMPILE)
			
			width = right-left
			width_plus = width + self.border_width*2
			height = top - bottom
			#height_plus = height + self.border_height*2
			depth = front-back
			depth_plus = depth + +self.border_depth*2
			
			
			# back
#			glPushMatrix()
#			glTranslate(x_center,y_center,back_plus)
#			glRotate(180,0,1,0)
#			#self.frame_face_basic(width,width_plus,height,height_plus)
#			self.frame_face(-width/2,-width/2-self.border_width,width/2,width/2+self.border_width,-height/2,-height/2-self.bottom_border_height,height/2,height/2+self.top_border_height)
#			self.frame_inner_shell_basic(width,height,thick_front,thick_back)
#			glPopMatrix()
#			
#			
#			#bottom
#			glPushMatrix()
#			glTranslate(x_center,bottom_plus,z_center)	
#			glRotate(90,1,0,0)
#			self.frame_face_basic(width,width_plus,depth,depth_plus)
#			self.frame_inner_shell_basic(width,depth,0,-self.bottom_border_height)
#			glPopMatrix()
#			
#			#top
#			glPushMatrix()
#			glTranslate(x_center,top_plus,z_center)	
#			glRotate(-90,1,0,0)
#			self.frame_face_basic(width,width_plus,depth,depth_plus)
#			self.frame_inner_shell_basic(width,depth,0,-self.top_border_height)
#			glPopMatrix()
#			
#			#right side
#			glPushMatrix()
#			glTranslate(right_plus,y_center,z_center)
#			glRotate(90,0,1,0)
#			self.frame_face(-depth/2,-depth/2-self.border_depth,depth/2,depth/2+self.border_depth,-height/2,-height/2-self.bottom_border_height,height/2,height/2+self.top_border_height)
#			#self.frame_face_basic(depth,depth_plus,height,height_plus)
#			self.frame_inner_shell_basic(depth,height,thick_front,thick_back)
#			glPopMatrix()
#			
#			# left
#			glPushMatrix()
#			glTranslate(left_plus,y_center,z_center)
#			glRotate(-90,0,1,0)
#			self.frame_face(-depth/2,-depth/2-self.border_depth,depth/2,depth/2+self.border_depth,-height/2,-height/2-self.bottom_border_height,height/2,height/2+self.top_border_height)
#			self.frame_inner_shell_basic(depth,height,thick_front,thick_back)
#			glPopMatrix()
			
			# front
			glPushMatrix()
			glTranslate(x_center,y_center,thick_front)
			#self.frame_face_basic(width,width_plus,height,height_plus)
			self.frame_face(-width/2,-width/2-self.border_width,width/2,width/2+self.border_width,-height/2,-height/2-self.bottom_border_height,height/2,height/2+self.top_border_height)
			glPopMatrix()
			
			glPushMatrix()
			glTranslate(x_center,y_center,0)
			self.frame_inner_shell_basic(width,height,thick_front,thick_back)
			glPopMatrix()
			
			
			
			
#			# Do the front facing strip first
#			glPushMatrix()
#			glTranslate(0,0,front)
#			self.frame_face(left,left_plus,right,right_plus,bottom,bottom_plus,top,top_plus)
#			glPopMatrix()
			
			# Do the back facing part
			glPushMatrix()
			glTranslate(x_center,y_center,thick_back)
			glRotate(180,0,1,0)
			self.frame_face(-width/2,-width/2-self.border_width,width/2,width/2+self.border_width,-height/2,-height/2-self.bottom_border_height,height/2,height/2+self.top_border_height)
			glPopMatrix()
			
			# Now do the border around the edges
			glPushMatrix()
			glTranslate(x_center,y_center,0)
			self.frame_outer_shell(-width/2-self.border_width,width/2+self.border_width,-height/2-self.bottom_border_height,height/2+self.top_border_height,thick_front,thick_back)
			glPopMatrix()
			
			# Now do the border around the inside edges
#			glPushMatrix()
#			self.frame_inner_shell(left,right,bottom,top,front,back)
#			glPopMatrix()
#				
			glEndList()
		else: print "error, the delete list operation failed"
		
	def viewport_update(self):
		self.corner_sets = []
		
		dims = self.object().get_lr_bt_nf()
		# plus, i.e. plus the border
		left =  dims[0]
		left_plus = left-self.border_width
		right = dims[1]
		right_plus = right + self.border_width
		
		bottom = dims[2]
		bottom_plus = bottom - self.bottom_border_height
		top =  dims[3]
		top_plus = top + self.top_border_height
		
		
		front =  dims[4]
		front_plus = front + self.border_depth
		back = dims[5]
		back_plus = back - self.border_depth
				
		points = []
		points.append((left_plus,bottom_plus,front_plus))
		points.append((left,bottom_plus,front_plus))
		points.append((left,bottom,front_plus))
		points.append((left_plus,bottom,front_plus))
		
		points.append((right,bottom_plus,front_plus))
		points.append((right_plus,bottom_plus,front_plus))
		points.append((right_plus,bottom,front_plus))
		points.append((right,bottom,front_plus))
		
		points.append((right,top,front_plus))
		points.append((right_plus,top,front_plus))
		points.append((right_plus,top_plus,front_plus))
		points.append((right,top_plus,front_plus))
		
		points.append((left_plus,top,front_plus))
		points.append((left,top,front_plus))
		points.append((left,top_plus,front_plus))
		points.append((left_plus,top_plus,front_plus))
		
		points.append((right-self.top_border_height,top,front_plus))
		points.append((right-self.top_border_height,top_plus,front_plus))
		
		self.unprojected_points = self.vdtools.unproject_points(points)
		unprojected = self.unprojected_points
		self.vdtools.set_mouse_coords(unprojected[0],unprojected[1],unprojected[2],unprojected[3])
		self.corner_sets.append( self.vdtools.get_corners() )
		
		self.vdtools.set_mouse_coords(unprojected[1],unprojected[4],unprojected[7],unprojected[2])
		self.corner_sets.append( self.vdtools.get_corners() )
		
		self.vdtools.set_mouse_coords(unprojected[4],unprojected[5],unprojected[6],unprojected[7])
		self.corner_sets.append( self.vdtools.get_corners() )
		
		self.vdtools.set_mouse_coords(unprojected[7],unprojected[6],unprojected[9],unprojected[8])
		self.corner_sets.append( self.vdtools.get_corners() )
		
		self.vdtools.set_mouse_coords(unprojected[8],unprojected[9],unprojected[10],unprojected[11])
		self.corner_sets.append( self.vdtools.get_corners() )
		
		self.vdtools.set_mouse_coords(unprojected[13],unprojected[16],unprojected[17],unprojected[14])
		self.corner_sets.append( self.vdtools.get_corners() )
		
		self.vdtools.set_mouse_coords(unprojected[12],unprojected[13],unprojected[14],unprojected[15])
		self.corner_sets.append( self.vdtools.get_corners() )
		
		self.vdtools.set_mouse_coords(unprojected[3],unprojected[2],unprojected[13],unprojected[12])
		self.corner_sets.append( self.vdtools.get_corners() )
		
		self.vdtools.set_mouse_coords(unprojected[16],unprojected[8],unprojected[11],unprojected[17])
		self.corner_sets.append( self.vdtools.get_corners() )
		
		# the whole thing this can save time
		self.vdtools.set_mouse_coords(unprojected[0],unprojected[5],unprojected[10],unprojected[15])
		self.corner_sets.append( self.vdtools.get_corners() )
		
	def isinwin(self,x,y):
		interception = False
		#print "I have this many corner sets",len(self.corner_sets)
		for i,p in enumerate(self.corner_sets):
			if self.vdtools.isinwinpoints(x,y,p):
				interception = True
				self.currently_selected_frame = i
				break
		return interception

	def isinwin_broadly(self,x,y):
		p = self.corner_sets[-1]
		if self.vdtools.isinwinpoints(x,y,p): return True
		return False
		
		