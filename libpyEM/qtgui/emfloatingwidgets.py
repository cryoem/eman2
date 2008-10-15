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

# EMFloatingWidgets.py  Steve Ludtke  08/06/2006
# An experimental version of emimage.py that displays an image in a 3D context
# using texture mapping. Not a fully fleshed out class at this point

from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from OpenGL.GL import *
from OpenGL.GLU import *
from valslider import ValSlider

from EMAN2 import *
from emimageutil import *
from emimage2d import *
from emimage3d import *
from emimagemx import *
from math import sqrt

try: from emimage import EMImage
except: pass

from emglobjects import EMViewportDepthTools, Camera2, EMBasicOpenGLObjects, Camera, viewport_width,viewport_height
from emimageutil import EMEventRerouter

height_plane = 500

white = (1.0,1.0,1.0,1.0)
yellow = (1.0,1.0,0.0,1.0)
red = (1.0,1.0,1.0,1.0)
blue = (0.0,0.0,1.0,1.0)
black = (0.0,0.0,0.0,0.0)
grey = (0.6,0.6,0.6,0.0)
green = (0.0,1.0,0.0,1.0)
purple = (1.0,0.0,1.0,1.0)
black = (0.0,0.0,0.0,1.0)
dark_purple = (0.4,0.0,0.4,1.0)
light_blue = (.2,.2,.6,1.0)

class EM3DBorderDecoration:
	'''
	An class for drawing borders around widgets
	The implementation of border decorations in EMAN2 floating widgets is based on the 
	Decorator pattern in the Gang of Four. Inheriting classes should provide draw(object)
	which draws a border around the given object.
	'''
	def __init__(self):
		pass
	
	def draw(self,object):
		pass

class EM3DPlainBorderDecoration(EM3DBorderDecoration):
	'''
	A plain border decoration
	'''
	
	FROZEN_COLOR = "frozen"
	DEFAULT_COLOR = "default"
	PERMISSABLE_COLOR_FLAGS = [FROZEN_COLOR,DEFAULT_COLOR]
	def __init__(self, object):
		EM3DBorderDecoration.__init__(self)
		self.border_width = 6
		self.border_height = 6
		self.border_depth = 6
		
		self.display_list = None
		self.force_update = False
		
		self.faulty = False
		if not isinstance(object,EMGLView3D) and not isinstance(object,EMGLView2D) and not isinstance(object,EMGLViewQtWidget) and not isinstance(object,EM3DGLVolume) and not isinstance(object,EM3DGLWindow) and not isinstance(object,EMGLWindow):
			print "error, border construction works only for EMGLView3D, EMGLView2D, EM3DGLVolume, and EMGLViewQtWidget objects"
			self.faulty = True
			return
		else: self.object = object
		
		self.color_flag = EM3DPlainBorderDecoration.DEFAULT_COLOR
	def __del__(self):
		self.__delete_list()
	
	
	def get_border_width(self): return self.border_width
	def get_border_height(self): return self.border_height
	def get_border_depth(self): return self.border_depth
	
	def set_color_flag(self,flag):
		if flag not in EM3DPlainBorderDecoration.PERMISSABLE_COLOR_FLAGS:
			print 'unknown color flag'
		else:
			self.color_flag = flag
	
	def set_force_update(self,val=True):
		self.force_update = val
	
	def draw(self,force_update=False):
		if force_update or self.force_update:
			self.__delete_list()
			self.force_update = False
		
		if self.display_list == None:
			if isinstance(self.object,EMGLView3D) or isinstance(self.object,EM3DGLVolume) or isinstance(self.object,EM3DGLWindow):
				self.__gen_3d_object_border_list()
			elif isinstance(self.object,EMGLView2D) or isinstance(self.object,EMGLViewQtWidget) or isinstance(self.object,EM2DGLWindow):
				self.__gen_2d_object_border_list()
			else:
				print "error, border decoration works only for EMGLView3D, EMGLView2D and EMGLViewQtWidget objects"
				return
		
		if self.display_list == None: return
		
		if self.color_flag ==  EM3DPlainBorderDecoration.DEFAULT_COLOR:
			glMaterial(GL_FRONT,GL_AMBIENT,green)
			glMaterial(GL_FRONT,GL_DIFFUSE,black)
			glMaterial(GL_FRONT,GL_SPECULAR,grey)
			glMaterial(GL_FRONT,GL_SHININESS,20.0)
			glColor(*white)
		elif self.color_flag ==  EM3DPlainBorderDecoration.FROZEN_COLOR:
			glMaterial(GL_FRONT,GL_AMBIENT,dark_purple)
			glMaterial(GL_FRONT,GL_DIFFUSE,blue)
			glMaterial(GL_FRONT,GL_SPECULAR,light_blue)
			glMaterial(GL_FRONT,GL_SHININESS,20.0)
			glColor(*blue)
		else:
			print "warning, unknown color flag, coloring failed"
			
		if self.display_list != None and self.display_list != 0:
			glCallList(self.display_list)
		else: print "weird"
	
	def __delete_list(self):
		if self.display_list != None:
			glDeleteLists(self.display_list,1)
			self.display_list = None
	
	def __gen_2d_object_border_list(self):
		
		#context = self.object.context()
		#context.makeCurrent()
		#print "made",context,"current"
		
		self.__delete_list()
		if self.display_list == None:
			
			width = self.object.width()
			height = self.object.height()
			
			# plus, i.e. plus the border
			left = -width/2.0
			left_plus = left-self.border_width
			right = -left
			right_plus = -left_plus
			
			bottom = -height/2.0
			bottom_plus = bottom - self.border_height
			top = -bottom
			top_plus = -bottom_plus
			
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
			self.__frame_face(left,left_plus,right,right_plus,bottom,bottom_plus,top,top_plus)
			glPopMatrix()
			
			# Do the back facing part
			glPushMatrix()
			glTranslate(0,0,back)
			glRotate(180,0,1,0)
			self.__frame_face(left,left_plus,right,right_plus,bottom,bottom_plus,top,top_plus)
			glPopMatrix()
			
			# Now do the border around the edges
			glPushMatrix()
			self.__frame_outer_shell(left_plus,right_plus,bottom_plus,top_plus,front,back)
			glPopMatrix()
			
			# Now do the border around the inside edges
			glPushMatrix()
			self.__frame_inner_shell(left,right,bottom,top,front,back)
			glPopMatrix()
			
			glEndList()
			
		else: print "error, the delete list operation failed"
	
	def __gen_3d_object_border_list(self):
		self.__delete_list()
	
		if self.display_list == None:
			thick_front = 0
			thick_back = -self.border_depth
				
			dims = self.object.get_lr_bt_nf()
			# plus, i.e. plus the border
			left =  dims[0]
			left_plus = left-self.border_width
			right = dims[1]
			right_plus = right + self.border_width
			
			bottom = dims[2]
			bottom_plus = bottom - self.border_height
			top =  dims[3]
			top_plus = top + self.border_height
			
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
			width_plus = width + self.border_height*2
			height = top - bottom
			height_plus = height + self.border_height*2
			depth = front-back
			depth_plus = depth + +self.border_height*2
			# front
			glPushMatrix()
			glTranslate(x_center,y_center,front_plus)
			self.__frame_face_basic(width,width_plus,height,height_plus)
			self.__frame_inner_shell_basic(width,height,thick_front,thick_back)
			glPopMatrix()
			
			# back
			glPushMatrix()
			glTranslate(x_center,y_center,back_plus)
			glRotate(180,0,1,0)
			self.__frame_face_basic(width,width_plus,height,height_plus)
			self.__frame_inner_shell_basic(width,height,thick_front,thick_back)
			glPopMatrix()
			
			#right side
			glPushMatrix()
			glTranslate(right_plus,y_center,z_center)
			glRotate(90,0,1,0)
			self.__frame_face_basic(depth,depth_plus,height,height_plus)
			self.__frame_inner_shell_basic(depth,height,thick_front,thick_back)
			glPopMatrix()
			
			# left
			glPushMatrix()
			glTranslate(left_plus,y_center,z_center)
			glRotate(-90,0,1,0)
			self.__frame_face_basic(depth,depth_plus,height,height_plus)
			self.__frame_inner_shell_basic(depth,height,thick_front,thick_back)
			glPopMatrix()
			
			#bottom
			glPushMatrix()
			glTranslate(x_center,bottom_plus,z_center)	
			glRotate(90,1,0,0)
			self.__frame_face_basic(width,width_plus,depth,depth_plus)
			self.__frame_inner_shell_basic(width,depth,thick_front,thick_back)
			glPopMatrix()
			
			#top
			glPushMatrix()
			glTranslate(x_center,top_plus,z_center)	
			glRotate(-90,1,0,0)
			self.__frame_face_basic(width,width_plus,depth,depth_plus)
			self.__frame_inner_shell_basic(width,depth,thick_front,thick_back)
			glPopMatrix()
			
			glEndList()
		else: print "error, the delete list operation failed"
		
	def __frame_inner_shell_basic(self,width,height,front,back):
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
	
	def __frame_outer_shell_basic(self,width_plus,height_plus,front,back):
		
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
		
	def __frame_face_basic(self,width,width_plus,height,height_plus):
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
	
	def __frame_face(self,left,left_plus,right,right_plus,bottom,bottom_plus,top,top_plus):
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
	
	def __frame_outer_shell(self,left_plus,right_plus,bottom_plus,top_plus,front,back):
		
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
		
	def __frame_inner_shell(self,left,right,bottom,top,front,back):
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
	
class EM3DGLVolume:
	'''
	a EM3DGLVolume has width(), height(), and depth() functions, and the associated private variables
	Inheriting functions should define the __determine_dimensions function which will be called implicitly if the update_dims flag is true
	In it they should define the member variables left,right,bottom,top,near and far.
	'''
	def __init__(self):
		self.left = 0
		self.right = 0
		self.bottom = 0
		self.top = 0
		self.near = 0
		self.far = 0

		self.update_dims = True # a flag, when set, cause this widget to re-evaluate its maximum dimensions
	
	def update(self):
		self.update_dims = True
		
	def width(self):
		if self.update_dims:
			self.determine_dimensions()
			
		return int(self.right -self.left)
		
	def height(self):
		if self.update_dims:
			self.determine_dimensions()
			
		return int(self.top - self.bottom)
		
	def depth(self):
		if self.update_dims:
			self.determine_dimensions()
			
		return int(self.near - self.far)
	
	def get_lr_bt_nf(self):
		'''
		get left-right bottom-top near-far
		'''
		if self.update_dims:
			self.determine_dimensions()
			
		return [self.left,self.right,self.bottom,self.top,0,self.far]
		
	def get_suggested_lr_bt_nf(self):
		'''
		A function that can be 
		'''
		# Either the deriving class provides this function (which is meant to provide optimal display dimensions)
		try: return self.nice_lr_bt_nf()
		except: return self.get_lr_bt_nf()
		




class EMGLRotorWidget(EM3DGLVolume):
	'''
	A display rotor widget - consists of an ellipse  with widgets 'attached' to it.
	Visually, the ellipse would lay in the plane perpendicular to the screen, and the widgets would be displayed in the plane
	of the screen and be attached to the ellipse and rotate around it interactively
	
	Actually, this is not true - I eventually made it so the ellipse could be rotated in all three directions, and starts in the 
	xz or yz planes.
	
	Think of the ellipse as being attached to the side of stretched cube. Then rotate the cube so that it points up, or to the left,
	or up and then a little to the left etc - the rotations are parameter that are established in the constructor.
	The widgets are stacked onto the ellipse for theta in [0,90], and travel alongs its boundary when the animate function is called.
	The leading widget is always displayed such that it is in xy plane (parallel to the screen).
		
	At the moment it works but it is incomplete and is a work in progress. Contact David Woolford for details.
	
	'''
	LEFT_ROTARY = "left_rotor"
	RIGHT_ROTARY = "right_rotor"
	TOP_ROTARY = "top_rotor"
	BOTTOM_ROTARY = "bottom_rotor"
	def __init__(self,parent,y_rot=15,z_rot=70,x_rot=-45,rotor_type=RIGHT_ROTARY,ellipse_a=200,ellipse_b=400):
		EM3DGLVolume.__init__(self)
		self.parent = parent
		self.widgets = []	# the list of widgets to display
		self.displayed_widget = -1 # the index of the currently displayed widget in self.widgets
		self.target_displayed_widget = -1 # the index of a target displayed widget, used for animation purposes
		self.rotations = 0 # an index
		self.anim_rotations = 0
		
		# elliptical parameters
		self.ellipse_a = ellipse_a # x direction - widgets displayed in position 0 will be a distance of 20 away from the center
		self.ellipse_b = ellipse_b # y direction - widgets at 90 degrees will be a distance of 100 into the screen, away from the elliptical center
		self.y_rot = y_rot*pi/180.0 # this is how much the ellipse should be rotated about the y axis - makes it so that the ellipse isn't 'head on' to the viewer, for visual effects. The orientation of the widgets is not affected (they are always oriented as though the ellipse is not rotated in plane.
		self.cos_y_rot = cos(self.y_rot)
		self.sin_y_rot = sin(self.y_rot)
		
		self.z_rot = z_rot*pi/180.0 # this is how much the ellipse should be rotated about the Z axis - makes it so that the ellipse isn't 'head on' to the viewer, for visual effects. The orientation of the widgets is not affected (they are always oriented as though the ellipse is not rotated in plane.
		self.cos_z_rot = cos(self.z_rot)
		self.sin_z_rot = sin(self.z_rot)
		
		self.x_rot = x_rot*pi/180.0 # this is how much the ellipse should be rotated about the Z axis - makes it so that the ellipse isn't 'head on' to the viewer, for visual effects. The orientation of the widgets is not affected (they are always oriented as though the ellipse is not rotated in plane.
		self.cos_x_rot = cos(self.x_rot)
		self.sin_x_rot = sin(self.x_rot)
		
		self.time = 0		# time indicates the current time used for the basis of animation.
		self.time_interval = .4 # 0.5 seconds for the animation to complete
		self.time_begin = 0 # records the time at which the animation was begun
		self.is_animated = False
		self.angle_information = None # will eventually be list used for both static display and animation
		
		self.dtheta_animations = None  # an array for dtheta animations
		
		self.animation_queue = []
		
		self.rotor_type = rotor_type # the rotor type, LEFT_ROTARY etc
	
		self.mmode = None	# mouse mode - used for potentially emitting signal
		
		self.limits = None # Potentially a tuple of 2 values denoting rotation limits
	
		self.allow_child_mouse_events = True # turn off if you dont want the widgets in the rotor to receive event
		
		self.angle_range = 90.0	# the angular amount of the ellipse that is occupied by resident widgets
		self.allow_target_wheel_events = True
		self.allow_target_translation = True
		
		
		t1 = Transform({"type":"xyz","ztilt":z_rot})
		t2 = Transform({"type":"xyz","ytilt":y_rot})
		t3 = Transform({"type":"xyz","xtilt":x_rot})
		self.rotation = t3*t2*t1

		#self.cam = Camera2(self) # currently unused
	def target_zoom_events_allowed(self,bool):
		self.allow_target_wheel_events = bool
	
	def target_translations_allowed(self,bool):
		self.allow_target_translation = bool
	
	def set_angle_range(self,angle_range):
		self.angle_range = angle_range
	
	def set_plane(self,plane):
		pass
	
	def set_model_matrix(self,matrix):
		pass
	
	def set_rotations(self,r):
		self.rotations = r	
	
	def __getitem__(self,idx):
		if len(self.widgets) == 0: return None
		i = idx-self.rotations
		if i != 0:
			i = i % len(self.widgets)
		#print len(self.widgets), "asking for i", i
		return self.widgets[i]
	
	def current_index(self):
		i = -self.rotations
		if i != 0:
			i = i % len(self.widgets)
		return i
	
	def __len__(self):
		return len(self.widgets)
	
	def set_child_mouse_events(self,bool):
		self.allow_child_mouse_events = bool
	
	def set_limits(self,limits):
		
		if len(limits) != 2:
			print "error, the limits are supposed to consist of 2 values"
			self.limits = None
		else:
			self.limits = limits
		
	def context(self):
		return self.parent.context()
	
	def set_frozen(self,frozen,idx=0):
		'''
		default is for the leading image to get the signal (idx=0)
		'''
		self.widgets[idx].set_frozen(frozen)
	
	def set_shapes(self,shapes,shrink,idx=0):
		self.widgets[idx].get_drawable().set_shapes(shapes,shrink)

	def set_mouse_mode(self,mode):
		self.mmode = mode
	
	def clear_widgets(self):
		self.widgets = []
		self.time = 0		# time indicates the current time used for the basis of animation.
		self.time_begin = 0 # records the time at which the animation was begun
		self.is_animated = False
		self.angle_information = None # will eventually be list used for both static display and animation
		self.dtheta_animations = None  # an array for dtheta animations
		self.animation_queue = []
		self.rotations = 0
	
	def add_widget(self,widget,set_current=False):
		'''
		adds a widget to those that already in the rotor
		if set_current is true will also make it so that the newly added widget has the focus of the rotor
		
		return 1 if a redraw is necessary
		return 0 if a redraw is not necessary 
		'''
		if not isinstance(widget,EMGLViewQtWidget) and not isinstance(widget,EMGLView2D) and not isinstance(widget,EMGLView3D) and not isinstance(widget,EM3DGLWindow)  and not isinstance(widget,EM2DGLWindow):
			print "error, can only add instances of EMGLViewQtWidget to the EMGLRotorWidget"
			return 0
		else:
			self.widgets.append(widget)
			self.update_dims = True
			if set_current == True and self.is_animated == False:
				self.target_displayed_widget =  len(self.widgets)-1
				self.__start_animation()

		return 1
		
	def is_flat(self):
		'''
		A hack - used by the calling function to determine if the first object in the rotor is flat against the screen.
		At the moment all that is returned is whether or not this object is animated, and this indeed indicates the we have
		'flatness', but only if the rotor is not rotated out of its original position by the user
		'''
		return not self.is_animated
	
	def animate(self,time):
		if self.time_begin == 0:
			self.time_begin = time
			
		self.time = time -self.time_begin

		if self.time > self.time_interval:
			self.time_begin = 0
			self.time = 0
			self.is_animated = False
			self.angle_information = None
			self.rotations += self.anim_rotations
			if self.mmode=="app":
				self.parent.emit(QtCore.SIGNAL("image_selected"),None,[(-self.rotations)%len(self.widgets)])
			elif self.mmode == "mxrotor":
				self.parent.update_rotor_position(-self.anim_rotations)
			return 0
			
		else:
			for i in range(len(self.widgets)):
				dt = self.__get_dt()
				dtheta = self.angle_information[i][2]-self.angle_information[i][1]
				self.angle_information[i][0] = self.angle_information[i][1]+dt*dtheta
				
			
		return 1
	
	def render(self):
		size = len(self.widgets)
		if size == 0:
			return
		dtheta = self.__get_current_dtheta()
		if self.angle_information == None or  size > len(self.angle_information):
			self.angle_information = []
			for i in range(len(self.widgets)):
				angle = i*dtheta
				self.angle_information.append([angle,angle,0])
		
		if self.is_animated: dt = self.__get_dt()
		
		points = []
		for i in range(size):
			n = self.angle_information[i][0]
			n_rad = n*pi/180.0
			if self.rotor_type == EMGLRotorWidget.LEFT_ROTARY:
				dx = self.ellipse_a*cos(n_rad)
				dy = 0
				dz = -self.ellipse_b*sin(n_rad)
				rot_v = [0,1,0]
			elif self.rotor_type == EMGLRotorWidget.RIGHT_ROTARY:
				dx = -self.ellipse_a*cos(n_rad)
				dy = 0
				dz = self.ellipse_b*sin(n_rad)
				rot_v = [0,1,0]
			elif self.rotor_type == EMGLRotorWidget.BOTTOM_ROTARY:
				dx = 0
				dy = self.ellipse_a*cos(n_rad)
				dz = self.ellipse_b*sin(n_rad)
				rot_v = [1,0,0]
			elif self.rotor_type == EMGLRotorWidget.TOP_ROTARY:
				dx = 0
				dy = -self.ellipse_a*cos(n_rad)
				dz = -self.ellipse_b*sin(n_rad)
				rot_v = [1,0,0]
			else:
				print "unsupported"
				return
			points.append([dx,dy,dz])
		
		self.__rotate_elliptical_points(points)
		
		for i in range(size):
			n = self.angle_information[i][0]

			idx = i-self.rotations
			if idx != 0:
				idx = idx % size
			#print idx,
			widget = self.widgets[idx]
			
			if self.rotor_type == EMGLRotorWidget.LEFT_ROTARY:
				h_width = widget.width()*0.5
				h_height = 0
			elif self.rotor_type == EMGLRotorWidget.RIGHT_ROTARY:
				h_width = -widget.width()*0.5
				h_height = 0
			elif self.rotor_type == EMGLRotorWidget.BOTTOM_ROTARY:
				h_width = 0
				h_height =  widget.height()*0.5
			elif self.rotor_type == EMGLRotorWidget.TOP_ROTARY:
				h_width = 0
				h_height =  -widget.height()*0.5
			else:
				print "unsupported"
				return
		
			glPushMatrix()
			glTranslate(points[i][0],points[i][1],points[i][2])
			glRotate(n,*rot_v)
			glTranslate(h_width,h_height,0)
			try: widget.draw()
			except: pass
			glPopMatrix()
			
	
	def __get_visible(self):
		'''
		returns the visible widget
		'''
		idx =  (-self.rotations) % (len(self.widgets))
		#print idx,self.rotations
		return self.widgets[idx]
	
	def __start_animation(self,counter_clockwise=True):
		self.is_animated = True
		self.__gen_animation_dthetas(counter_clockwise)
		self.parent.register_animatable(self)
		
	def wheelEvent(self,event):
		
		if not self.is_animated and self.allow_child_mouse_events and self.allow_target_wheel_events: 
			if self.__get_visible().isinwin(event.x(),viewport_height()-event.y()) :
				self.__get_visible().wheelEvent(event)
				return

		if event.delta() > 0:
			if self.target_displayed_widget == 0:
				self.target_displayed_widget =  len(self.widgets)-1
			else: 
				self.target_displayed_widget -= 1
		else:
			if self.target_displayed_widget == len(self.widgets)-1:
				self.target_displayed_widget = 0
			else:
				self.target_displayed_widget += 1
		#print self.target_displayed_widget
		if not self.is_animated: 
			self.__start_animation(event.delta()> 0)
		else: self.__update_animation_dthetas(event.delta() > 0)
		
	
	def mousePressEvent(self,event):
		if self.is_animated: return
		if self.allow_child_mouse_events and self.__get_visible().isinwin(event.x(),viewport_height()-event.y()) :
			self.__get_visible().mousePressEvent(event)
	
	def mouseReleaseEvent(self,event):
		if self.is_animated: return
		if self.allow_child_mouse_events and self.__get_visible().isinwin(event.x(),viewport_height()-event.y()) :
			self.__get_visible().mouseReleaseEvent(event)
	
	def mouseMoveEvent(self,event):
		if self.is_animated: return
		if self.allow_child_mouse_events and self.__get_visible().isinwin(event.x(),viewport_height()-event.y()) :
			self.__get_visible().mouseMoveEvent(event)
	
	def toolTipEvent(self,event):
		if self.is_animated: return
		if self.allow_child_mouse_events and self.__get_visible().isinwin(event.x(),viewport_height()-event.y()) :
			self.__get_visible().toolTipEvent(event)
	
	
	def keyPressEvent(self,event):
		print "should act on up and down"
	
	def determine_dimensions(self):
		print "in determine dimensions"
		# first determine the ellipse offset - this is related to the inplane rotation of the ellipse. We want the back point corresponding to
		# a rotation of 90 degrees around the boundary of the ellipse (starting at [ellipse_a/2,0,0]) to be visible
		
		if self.rotor_type == EMGLRotorWidget.LEFT_ROTARY or self.rotor_type == EMGLRotorWidget.RIGHT_ROTARY :
			interesting_points = [[-self.ellipse_a,0,0],[self.ellipse_a,0,0],[0,0,-self.ellipse_b],[0,0,self.ellipse_b]]
		elif self.rotor_type == EMGLRotorWidget.BOTTOM_ROTARY or self.rotor_type == EMGLRotorWidget.TOP_ROTARY :	
			interesting_points = [[0,-self.ellipse_a,0],[0,self.ellipse_a,0],[0,0,-self.ellipse_b],[0,0,self.ellipse_b]]
		else:
			print "unsupported"
		
		self.__rotate_elliptical_points(interesting_points)
			
		
		max_z = interesting_points[0][2]
		min_z = interesting_points[0][2]
		min_y = interesting_points[0][1]
		max_y = interesting_points[0][1]
		min_x = interesting_points[0][0]
		max_x = interesting_points[0][0]
		
		for n in range(1,len(interesting_points)):
			val = interesting_points[n][2]
			if val > max_z: max_z = val
			if val < min_z: min_z = val
		
		height = -1
		width = -1
		for i in self.widgets:
			if i.width() > width: width = i.width()
			if i.height() > height: height = i.height()
		
		if self.rotor_type == EMGLRotorWidget.LEFT_ROTARY:
			
			self.bottom = -height/2.0 + min_z
			self.top = height/2.0 + max_z

			if self.y_rot < 0:
				self.left  = interesting_points[1][0]
			elif self.y_rot >= 0:
				self.left = interesting_points[2][0]
			self.right = interesting_points[1][0]+width
	
			self.far = -interesting_points[3][1] -width
			self.near = interesting_points[3][1] + width
			
		elif self.rotor_type == EMGLRotorWidget.RIGHT_ROTARY:
			
			self.bottom = -height/2.0 + min_z
			self.top = height/2.0 + max_z
			
			if self.y_rot < 0:
				self.right = interesting_points[2][0]
			elif self.y_rot >= 0:
				self.right = interesting_points[0][0]
			self.left = interesting_points[0][0] - width

			self.far = -interesting_points[3][1] - width
			self.near = interesting_points[3][1] + width
			
		elif self.rotor_type == EMGLRotorWidget.BOTTOM_ROTARY:
			
			self.left = -width/2 + min_x
			self.right = width/2 + max_x
			
			if self.x_rot < 0:
				self.bottom = -interesting_points[3][1]
			elif self.x_rot >= 0:
				self.bottom = interesting_points[2][1]
				
			self.top = height + interesting_points[1][1]
			
			self.near = interesting_points[3][2]+height
			self.far =  -interesting_points[3][2] -height
			
		elif self.rotor_type == EMGLRotorWidget.TOP_ROTARY:
			
			self.left = -width/2 + min_x
			self.right = width/2 + max_x
			
			if self.x_rot < 0:
				self.top =  interesting_points[2][1]
			elif self.x_rot >= 0:
				self.top = -interesting_points[3][1]
				
			self.bottom = -height + interesting_points[0][1]
			
			self.near = interesting_points[3][2] + height
			self.far =  -interesting_points[3][2] - height
		
		else:
			print "unsupported"
			return
		
		self.update_dims = False

	def get_lr_bt_nf_2(self):
		from math import sqrt
		zmax = (self.rotation.at(2,0)*self.ellipse_a)**2
		zmax += (self.rotation.at(2,1)*self.ellipse_b)**2
		zmax = sqrt(zmax)
		zmin = -zmax
		
		ymax = (self.rotation.at(1,0)*self.ellipse_a)**2
		ymax += (self.rotation.at(1,1)*self.ellipse_b)**2
		ymax = sqrt(ymax)
		ymin = -ymax
		
		xmax = (self.rotation.at(0,0)*self.ellipse_a)**2
		xmax += (self.rotation.at(0,1)*self.ellipse_b)**2
		xmax = sqrt(xmax)
		xmin = -xmax
		
		return [xmin,xmax,ymin,ymax,zmin,zmax]
	
	def get_lr_bt_nf(self):
		size = len(self.widgets)
		
		#if size == 0:
		#return
		#dtheta = self.__get_current_dtheta()
		#if self.angle_information == None or  size > len(self.angle_information):
			#self.angle_information = []
			#for i in range(len(self.widgets)):
				#angle = i*dtheta
				#self.angle_information.append([angle,angle,0])
		points = []
		for n in [0,90]:
			#n = self.angle_information[i][0]
			n_rad = n*pi/180.0
			if self.rotor_type == EMGLRotorWidget.LEFT_ROTARY:
				dx = self.ellipse_a*cos(n_rad)
				dy = 0
				dz = -self.ellipse_b*sin(n_rad)
				rot_v = [0,1,0]
			elif self.rotor_type == EMGLRotorWidget.RIGHT_ROTARY:
				dx = -self.ellipse_a*cos(n_rad)
				dy = 0
				dz = self.ellipse_b*sin(n_rad)
				rot_v = [0,1,0]
			elif self.rotor_type == EMGLRotorWidget.BOTTOM_ROTARY:
				dx = 0
				dy = self.ellipse_a*cos(n_rad)
				dz = self.ellipse_b*sin(n_rad)
				rot_v = [1,0,0]
			elif self.rotor_type == EMGLRotorWidget.TOP_ROTARY:
				dx = 0
				dy = -self.ellipse_a*cos(n_rad)
				dz = -self.ellipse_b*sin(n_rad)
				rot_v = [1,0,0]
			else:
				print "unsupported"
				return
			points.append([dx,dy,dz])
		
		self.__rotate_elliptical_points(points)
		
		p1 = points[0]
		pn = points[len(points)-1]
		#print "points",p1,pn
		if self.rotor_type == EMGLRotorWidget.LEFT_ROTARY:
			left = pn[0]
			right = p1[0]+self.widgets[0].width()+ 10
			top = p1[1]+ self.widgets[0].height()/2 + 10
			bottom = pn[1]-self.widgets[0].height()/2 - 10
			near = p1[2]
			far = pn[2]-self.widgets[0].width()
			
			#print "return",[left,right,bottom,top,near,far],
			return [left,right,bottom,top,near,far]
		elif self.rotor_type == EMGLRotorWidget.TOP_ROTARY:
			if p1[0] < pn[0] :
				left = p1[0] -self.widgets[0].width()/2 - 10
				right = pn[0]+ self.widgets[0].widths()/2 + 10
			else:
				left = pn[0] -self.widgets[0].width()/2 - 10
				right = p1[0]+ self.widgets[0].width()/2 + 10
				
			top = pn[1]
			bottom = p1[1]-self.widgets[0].height() - 10
			
			near = p1[2]
			far = pn[2]-self.widgets[0].width()
			
			#print "return",[left,right,bottom,top,near,far],self.widgets[0].height()/2
			return [left,right,bottom,top,near,far]
		
		
	def nice_lr_bt_nf(self):
	
		'''
		This function returns what it thinks are optimal values for left, right, top, bottom, near and far,
		for the purposes of the display
		'''
		
		if self.rotor_type == EMGLRotorWidget.LEFT_ROTARY or self.rotor_type == EMGLRotorWidget.RIGHT_ROTARY :
			interesting_points = [[-self.ellipse_a,0,0],[self.ellipse_a,0,0],[0,0,-self.ellipse_b],[0,0,self.ellipse_b]]
		elif self.rotor_type == EMGLRotorWidget.BOTTOM_ROTARY or self.rotor_type == EMGLRotorWidget.TOP_ROTARY :	
			interesting_points = [[0,-self.ellipse_a,0],[0,self.ellipse_a,0],[0,0,-self.ellipse_b],[0,0,self.ellipse_b]]
		else:
			print "unsupported"
		
		self.__rotate_elliptical_points(interesting_points)
			
		
		max_z = interesting_points[0][2]
		min_z = interesting_points[0][2]
		#min_y = interesting_points[0][1]
		#max_y = interesting_points[0][1]
		min_x = interesting_points[0][0]
		max_x = interesting_points[0][0]
		
		for n in range(1,len(interesting_points)):
			val = interesting_points[n][2]
			if val > max_z: max_z = val
			if val < min_z: min_z = val
		
		height = -1
		width = -1
		for i in self.widgets:
			if i.width() > width: width = i.width()
			if i.height() > height: height = i.height()
		
		if self.rotor_type == EMGLRotorWidget.LEFT_ROTARY:
			if interesting_points[3][2] > 0:
				bottom = -height/2.0 + min_z
				top = height/2.0
			else:
				bottom = -height/2.0 
				top = height/2.0 + max_z

			if self.y_rot < 0:
				left  = interesting_points[1][0]
			elif self.y_rot >= 0:
				left = interesting_points[2][0]
			right = interesting_points[1][0]+width
	
			far = -interesting_points[3][1] -width
			near = interesting_points[3][1] + width
			
			return [left,right,bottom,top,near,far]
			
		elif self.rotor_type == EMGLRotorWidget.TOP_ROTARY:
			left = -width/2 + min_x
			right = width/2 + max_x
			
			if self.x_rot < 0:
				top =  interesting_points[2][1]
			elif self.x_rot >= 0:
				top = -interesting_points[3][1]
				
			bottom = -height + interesting_points[0][1]
			
			near = interesting_points[3][2] + height
			far =  -interesting_points[3][2] - height
			print get_lr_bt_nf()
			print "and",[left,right,bottom,top,near,far]
			return [left,right,bottom,top,near,far]
		else:
			error_string = "error"
			raise error_string
		
	def resize_event(self,width,height):
		self.set_update_P_inv(True)
		
	def set_update_P_inv(self,val=True):
		for widget in self.widgets:
			widget.set_update_P_inv(val)
	
	# PRIVATE 
	
	def __rotate_elliptical_points(self,points):
		
		#if self.z_rot != 0:
			#for p in points:
				#x = p[0]
				#y = p[1]
				#p[0] = self.cos_z_rot*x - self.sin_z_rot*y
				#p[1] = -self.sin_z_rot*x + self.cos_z_rot*y
		
		#if self.y_rot != 0:
			#for p in points:
				#x = p[0]
				#z = p[2]
				#p[0] = self.cos_y_rot*x + self.sin_y_rot*z
				#p[2] = -self.sin_y_rot*x +self.cos_y_rot*z
		
		#if self.x_rot != 0:
			#for p in points:
				#y = p[1]
				#z = p[2]
				#p[1] = self.cos_x_rot*y + self.sin_x_rot*z
				#p[2] = -self.sin_x_rot*y +self.cos_x_rot*z
				
		
		for p in points:
			v = Vec3f(p)
			v_d = self.rotation*v
			p[0] = v_d[0]
			p[1] = v_d[1]
			p[2] = v_d[2]
			
	def __get_current_dtheta_range(self):
		if len(self.widgets) == 1: return [0,0]
		elif self.rotor_type == EMGLRotorWidget.LEFT_ROTARY or self.rotor_type == EMGLRotorWidget.TOP_ROTARY: 
			return [self.angle_range/(len(self.widgets)-1),-self.angle_range/(len(self.widgets)-1)]
		elif self.rotor_type == EMGLRotorWidget.RIGHT_ROTARY or self.rotor_type == EMGLRotorWidget.BOTTOM_ROTARY:
			 return [-self.angle_range/(len(self.widgets)-1),self.angle_range/(len(self.widgets)-1)]
		else: 
			print "unsupported"
			return 0
		
	def __get_current_dtheta(self):
		if len(self.widgets) == 1: return 0
		elif self.rotor_type == EMGLRotorWidget.LEFT_ROTARY or self.rotor_type == EMGLRotorWidget.TOP_ROTARY: 
			return self.angle_range/(len(self.widgets)-1)
		elif self.rotor_type == EMGLRotorWidget.RIGHT_ROTARY or self.rotor_type == EMGLRotorWidget.BOTTOM_ROTARY:
			return -self.angle_range/(len(self.widgets)-1)
		else: 
			print "unsupported"
			return 0
	def __get_dt(self):
		return sin(self.time/self.time_interval*pi/2)
	
	def __get_dt_2(self):
		return (cos(2*pi*self.time/self.time_interval) + 1)/32.0 + 15.0/16.0
	
	def explicit_animation(self,rotations):
		if self.is_animated == False: 
			self.is_animated = True
			self.__gen_animation_dthetas(True,rotations)
			self.parent.register_animatable(self)
		else:
			self.__update_animation_dthetas(True,rotations)

	def __update_animation_dthetas(self,counter_clockwise=True,rotations=1):
		if not counter_clockwise: #clockwise OpenGL rotations are negative
			#self.rotations -= 1
			self.anim_rotations -= rotations
			rotations = -rotations
		else: # counter clockwise OpenGL rotations are positive
			#self.rotations += 1
			self.anim_rotations += rotations
			rotations = rotations
		
		self.time_begin = time.time()
		n = len(self.widgets)
		dtheta = self.__get_current_dtheta()
		dt = self.__get_dt()
		current_thetas = []
		for i in range(n): 
			current_thetas.append(i*dtheta)
		
		if rotations == -1: # clockwise
			for i in range(n):
				for rot in range(self.anim_rotations,self.anim_rotations-1,-1):
					idx1 = (i+rot+1)%n
					idx2 = (i+rot)%n
					self.angle_information[i][1] =  self.angle_information[i][0]
					if idx2 != (n-1):
						self.angle_information[i][2] +=  current_thetas[idx2] - current_thetas[idx1]
					else:
						if self.rotor_type == EMGLRotorWidget.LEFT_ROTARY or self.rotor_type == EMGLRotorWidget.TOP_ROTARY:
							self.angle_information[i][2] +=  -( 360-(current_thetas[idx2] - current_thetas[idx1]))
						elif  self.rotor_type == EMGLRotorWidget.RIGHT_ROTARY  or self.rotor_type == EMGLRotorWidget.BOTTOM_ROTARY:
							self.angle_information[i][2] +=  360-(current_thetas[idx1] - current_thetas[idx2])
		elif rotations >= 1: # counterclockwise
			for i in range(len(self.widgets)):
				for rot in range(self.anim_rotations,self.anim_rotations+1):
					idx1 = (i+rot-1)%n
					idx2 = (i+rot)%n
					self.angle_information[i][1] =  self.angle_information[i][0]
					if idx1 != (n-1):
						self.angle_information[i][2] +=  current_thetas[idx2] - current_thetas[idx1]
					else:
						if self.rotor_type == EMGLRotorWidget.LEFT_ROTARY or self.rotor_type == EMGLRotorWidget.TOP_ROTARY:
							self.angle_information[i][2] +=  360-(current_thetas[idx1] - current_thetas[idx2])
						elif  self.rotor_type == EMGLRotorWidget.RIGHT_ROTARY or self.rotor_type == EMGLRotorWidget.BOTTOM_ROTARY:
							self.angle_information[i][2] +=  -( 360-(current_thetas[idx2] - current_thetas[idx1]))
						

	def __gen_animation_dthetas(self,counter_clockwise=True,rotations=1):
		# find the shortest path from the currently displayed to the target displayed widget
				
		#c_distance = self.target_displayed_widget # clockwise distance, looking down -y
		#cc_distance = len(self.widgets) - self.target_displayed_widget #conterclockwise distance, looking down -y
		
		n = len(self.widgets)
		dtheta = self.__get_current_dtheta()
		self.angle_information = []
		current_thetas = []
		for i in range(n): 
			# information stored in current, start, end format
			current_thetas.append(i*dtheta)
			self.angle_information.append([i*dtheta,i*dtheta,i*dtheta])
	
		if not counter_clockwise: #clockwise OpenGL rotations are negative
			#self.rotations -= 1
			self.anim_rotations = -rotations
			rotations = -rotations
		else: # counter clockwise OpenGL rotations are positive
			#self.rotations += 1
			self.anim_rotations = rotations
			rotations = rotations
		
		
		if rotations <= -1: # clockwise
			for i in range(n):
				for rot in range(-1,rotations-1,-1):
					idx1 = (i+rot+1)%n
					idx2 = (i+rot)%n
					if idx2 != (n-1):
						self.angle_information[i][2] +=  current_thetas[idx2] - current_thetas[idx1]
					else:
						if self.rotor_type == EMGLRotorWidget.LEFT_ROTARY or self.rotor_type == EMGLRotorWidget.TOP_ROTARY:
							self.angle_information[i][2] +=  -( 360-(current_thetas[idx2] - current_thetas[idx1]))
						elif  self.rotor_type == EMGLRotorWidget.RIGHT_ROTARY or self.rotor_type == EMGLRotorWidget.BOTTOM_ROTARY:
							self.angle_information[i][2] +=   360-(current_thetas[idx1] - current_thetas[idx2])
		elif rotations >= 1: # counterclockwise
			for i in range(len(self.widgets)):
				for rot in range(1,rotations+1):
					idx1 = (i+rot-1)%n
					idx2 = (i+rot)%n
					if idx1 != (n-1):
						self.angle_information[i][2] +=  current_thetas[idx2] - current_thetas[idx1]
					else:
						if self.rotor_type == EMGLRotorWidget.LEFT_ROTARY or self.rotor_type == EMGLRotorWidget.TOP_ROTARY:
							self.angle_information[i][2] +=  360-(current_thetas[idx1] - current_thetas[idx2])
						elif  self.rotor_type == EMGLRotorWidget.RIGHT_ROTARY or self.rotor_type == EMGLRotorWidget.BOTTOM_ROTARY:
							self.angle_information[i][2] +=   -( 360-(current_thetas[idx2] - current_thetas[idx1]))
						
		else:
			print "error - can't rotate when rotations are set to 0"

class EMGLWindow:
	def __init__(self,parent,drawable):
		self.parent = parent
		self.drawable = drawable
		self.cam = Camera2(self) # a camera/orientation/postion object
		self.vdtools = EMViewportDepthTools(self) # viewport depth tools - essential for rerouting mouse events
		
		self.allow_target_wheel_events = True
		
		self.allow_target_translations = True
		
		self.texture_lock = 0
	
	def lock_texture(self):
		self.texture_lock += 1
	
	def unlock_texture(self):
		self.texture_lock -= 1
	
	def get_drawable(self): return self.drawable
	
	def isinwin(self,x,y): raise
	
	def context(self):
		return self.parent.context() # Fixme, this will raise if the parent isn't a QtOpenGL.QGLWidget, which is a more recent development. However, this function may become redundant, it is not currently used but it potentially could be
	
	def target_wheel_events_allowed(self,bool):
		self.allow_target_wheel_events = bool
	
	def target_translations_allowed(self,bool):
		self.allow_target_translations = bool
	
	def allow_camera_rotations(self,bool):
		self.cam.allow_camera_rotations(bool)
	
	
	def get_inspector(self):
		return self.drawable.get_inspector()
	
	def mousePressEvent(self, event):
		if self.texture_lock > 0: return
		
		if event.modifiers() == Qt.ControlModifier or (event.button()==Qt.RightButton and not self.allow_target_translations):
			self.cam.set_plane('xy')
			self.cam.mousePressEvent(event)
		else:
			l=self.vdtools.mouseinwin(event.x(),viewport_height()-event.y(),self.width(),self.height())
			qme=QtGui.QMouseEvent(event.type(),QtCore.QPoint(l[0],l[1]),event.button(),event.buttons(),event.modifiers())
			self.drawable.mousePressEvent(qme)

	def wheelEvent(self,event):
		
		if self.texture_lock > 0: return
		
		try: e = event.modifiers()
		except: return
		if event.modifiers() == Qt.ControlModifier:
			self.cam.set_plane('xy')
			self.cam.wheelEvent(event)
		else:
			if self.allow_target_wheel_events:
				self.drawable.wheelEvent(event)
			
		#self.updateGL()
	
	def mouseMoveEvent(self,event):
		if self.texture_lock > 0: return
		
		if event.modifiers() == Qt.ControlModifier:
			self.cam.set_plane('xy')
			self.cam.mouseMoveEvent(event)
		else:
			l=self.vdtools.mouseinwin(event.x(),viewport_height()-event.y(),self.width(),self.height())
			qme=QtGui.QMouseEvent(event.type(),QtCore.QPoint(l[0],l[1]),event.button(),event.buttons(),event.modifiers())
			self.drawable.mouseMoveEvent(qme)
	
	def mouseDoubleClickEvent(self, event):
		if self.texture_lock > 0: return
		
		if event.modifiers() == Qt.ControlModifier:
			self.cam.set_plane('xy')
			self.cam.mouseMoveEvent(event)
		else:
			l=self.vdtools.mouseinwin(event.x(),viewport_height()-event.y(),self.width(),self.height())
			qme=QtGui.QMouseEvent(event.type(),QtCore.QPoint(l[0],l[1]),event.button(),event.buttons(),event.modifiers())
			self.drawable.mouseMoveEvent(qme)
	
		#self.updateGL()

	def mouseReleaseEvent(self,event):
		if self.texture_lock > 0: return
		
		if event.modifiers() == Qt.ControlModifier:
			self.cam.set_plane('xy')
			self.cam.mouseReleaseEvent(event)
		else:
			l=self.vdtools.mouseinwin(event.x(),viewport_height()-event.y(),self.width(),self.height())
			qme=QtGui.QMouseEvent(event.type(),QtCore.QPoint(l[0],l[1]),event.button(),event.buttons(),event.modifiers())
			self.drawable.mouseReleaseEvent(qme)
	
	def toolTipEvent(self,event):
		if self.texture_lock > 0: return
		self.drawable.toolTipEvent(event)
	
	def keyPressEvent(self,event):
		if self.texture_lock > 0: return
		self.drawable.keyPressEvent(event)
	
	def eye_coords_dif(self,x1,y1,x2,y2,mdepth=True):
		return self.vdtools.eye_coords_dif(x1,y1,x2,y2,mdepth)
	
	def set_update_P_inv(self,val=True):
		self.vdtools.set_update_P_inv(val)
	
	def leaveEvent(self):
		pass
	
	def say_something(self):
		print "hello from EMGLVolume"

class EM3DGLWindow(EMGLWindow):
	'''
	A class for managing a 3D object as an interactive widget
	'''
	def __init__(self,parent,gl_view):
		EMGLWindow.__init__(self,parent,gl_view)

		self.corner_sets = [] # corner index sets of the 3D volume (?)
		self.planes = []	# string names for visible planes, used in draw
		self.model_matrices = [] # stores up to 3 OpenGL model view matrices (4x4 lists or tuples)
		
		self.draw_frame = True
		
		self.decoration = EM3DPlainBorderDecoration(self)
		
		self.w = self.parent.width()
		self.h = self.parent.height()
		self.d = self.parent.height()
		self.init_flag = True
		self.texure_lock = 0
	
	def lock_texture(self):
		self.texture_lock += 1
	
	def unlock_texture(self):
		self.texture_lock -= 1
	
	def context(self):
		return self.parent.context()
	
	def set_draw_frame(self,bool):
		self.draw_frame = bool

	def width(self):
		return self.w
		
	def height(self):
		return self.h
	
	def depth(self):
		return self.d
	
	def resize(self,width,height):
		self.w = width
		self.h = height
	
	def get_lr_bt_nf(self):
		
		#return [-self.w/2,self.w/2,-self.h/2,self.h/2,self.d/2,-self.d/2]
		return self.drawable.get_lr_bt_nf()
	def get_my_dims(self):
		
		return [-self.w/2,self.w/2,-self.h/2,self.h/2,self.d/2,-self.d/2]
	def __atomic_draw_frame(self,plane_string):
		
		[mc00,mc01,mc11,mc10] = self.vdtools.get_corners()
		a = Vec3f(mc10[0]-mc00[0],mc10[1]-mc00[1],mc10[2]-mc00[2])
		b = Vec3f(mc01[0]-mc00[0],mc01[1]-mc00[1],mc01[2]-mc00[2])
		c = a.cross(b)
		if ( c[2] < 0 ):
			#print "facing backward"
			return False
		self.corner_sets.append([mc00,mc01,mc11,mc10] )
		
		self.planes.append((plane_string))
		self.model_matrices.append(self.vdtools.getModelMatrix())
		return True
	
	#HACK ALERT
	def render(self):
		self.draw()
	
	def draw(self):
		#clear everything
		
		if self.init_flag == True:
			lrt = self.drawable.get_lr_bt_nf()
			self.w = lrt[1]-lrt[0]
			self.h = lrt[3]-lrt[2]
			self.d = lrt[4]-lrt[5]
			self.init_flag = False
	
		self.corner_sets = []
		self.planes = []
		self.model_matrices = []
		
		#glTranslatef(0,0,-self.depth()/2.0) # is this necessary?
		self.cam.position()
		
		lighting = glIsEnabled(GL_LIGHTING)
		glEnable(GL_LIGHTING) # lighting is on to make the borders look nice
			
		p = self.get_my_dims()
		points = []
		points.append((p[0],p[2],p[5]))
		points.append((p[0],p[2],p[4]))
		points.append((p[0],p[3],p[4]))
		points.append((p[0],p[3],p[5]))
		points.append((p[1],p[2],p[5]))
		points.append((p[1],p[2],p[4]))
		points.append((p[1],p[3],p[4]))
		points.append((p[1],p[3],p[5]))
		unprojected = self.vdtools.unproject_points(points)
		
		#print "the unprojected points are"
		
		# left zy plane
		glPushMatrix()
		self.vdtools.set_mouse_coords(unprojected[0],unprojected[1],unprojected[2],unprojected[3])
		self.__atomic_draw_frame('zy')
		glPopMatrix()
		
		glPushMatrix()
		self.vdtools.set_mouse_coords(unprojected[5],unprojected[4],unprojected[7],unprojected[6])
		self.__atomic_draw_frame('yz')
		glPopMatrix()
		
		glPushMatrix()
		self.vdtools.set_mouse_coords(unprojected[0],unprojected[4],unprojected[5],unprojected[1])
		self.__atomic_draw_frame('xz')
		glPopMatrix()
		
		glPushMatrix()
		self.vdtools.set_mouse_coords(unprojected[3],unprojected[2],unprojected[6],unprojected[7])
		self.__atomic_draw_frame('zx')
		glPopMatrix()
		
		glPushMatrix()
		self.vdtools.set_mouse_coords(unprojected[0],unprojected[3],unprojected[7],unprojected[4])
		self.__atomic_draw_frame('yx')
		glPopMatrix()
		
		glPushMatrix()
		self.vdtools.set_mouse_coords(unprojected[1],unprojected[5],unprojected[6],unprojected[2])
		self.__atomic_draw_frame('xy')
		glPopMatrix()
	
		glPushMatrix()
		lrt = self.drawable.get_lr_bt_nf()
		glTranslate(-(lrt[1]+lrt[0])/2.0,-(lrt[3]+lrt[2])/2.0,-lrt[4])
		
		glPushMatrix()
		self.drawable.render()
		glPopMatrix()

		glEnable(GL_LIGHTING) # lighting is on to make the borders look nice
		glEnable(GL_DEPTH_TEST) # lighting is on to make the borders look nice
		if self.draw_frame:
			self.decoration.draw()
			
		glPopMatrix()

	def isinwin(self,x,y):
		interception = False
		#print "I have this many corner sets",len(self.corner_sets)
		for i,p in enumerate(self.corner_sets):
			if self.vdtools.isinwinpoints(x,y,p):
				interception = True
				#self.cam.set_plane(self.planes[i])
				#print "plane is",self.planes[i]
				self.vdtools.setModelMatrix(self.model_matrices[i])
				break
		return interception
		
	def updateGL(self):
		self.parent.updateGL()

	def emit(self,*args,**kargs):
		#print "emit me"
		self.parent.emit(*args,**kargs)
		
class EM3DGLWindowOverride(EM3DGLWindow):
	'''
	A class for managing a 3D object as an interactive widget
	'''
	def __init__(self,parent,gl_view):
		EM3DGLWindow.__init__(self,parent,gl_view)
		
	def mousePressEvent(self, event):
		if self.texture_lock > 0: return
		if event.modifiers() == Qt.ControlModifier or (event.button()==Qt.RightButton and not self.allow_target_translations):
			self.cam.set_plane('xy')
			self.cam.mousePressEvent(event)
		else: self.drawable.mousePressEvent(event)

	
	def mouseMoveEvent(self,event):
		if self.texture_lock > 0: return
		if event.modifiers() == Qt.ControlModifier:
			self.cam.set_plane('xy')
			self.cam.mouseMoveEvent(event)
		else:
			self.drawable.mouseMoveEvent(event)
	
	def mouseDoubleClickEvent(self, event):
		if self.texture_lock > 0: return
		if event.modifiers() == Qt.ControlModifier:
			self.cam.set_plane('xy')
			self.cam.mouseMoveEvent(event)
		else: self.drawable.mouseMoveEvent(event)
	
		#self.updateGL()

	def mouseReleaseEvent(self,event):
		if self.texture_lock > 0: return
		if event.modifiers() == Qt.ControlModifier:
			self.cam.set_plane('xy')
			self.cam.mouseReleaseEvent(event)
		else: self.drawable.mouseReleaseEvent(event)
			
class EMGLView3D(EM3DGLVolume,EMEventRerouter):
	"""
	A view of an EMAN2 3D type, such as an isosurface or a 
	volume rendition, etc.
	"""
	def __init__(self, parent,image=None):
		EM3DGLVolume.__init__(self)
	
		self.parent = parent
		self.cam = Camera2(self)
		self.cam.motiondull = 3.0
		
		if isinstance(parent,EMImage3DGUIModule):
			self.drawable = parent
			self.w = self.drawable.width()
			self.h = self.drawable.height()
			self.d = self.drawable.height() # height of window
		else:
		#self.cam.setCamTrans('default_z',-parent.get_depth_for_height(height_plane))
		
			#self.w = image.get_xsize()	# width of window
			#self.h = image.get_ysize()	# height of window
			#self.d = image.get_zsize()	# depth of the window
					
			self.w = self.parent.width()	# width of window
			self.h = self.parent.height() # height of window
			self.d = self.parent.height() # height of window
			
	
			self.image = image # FIXME is this allright?
			self.drawable = EMImage3DModule(image,self)		# the object that is drawable (has a draw function)
			
			self.drawable.cam.basicmapping = True
			self.drawable.cam.motiondull = 3.0
		
		EMEventRerouter.__init__(self,self.drawable)
		self.sizescale = 1.0		# scale/zoom factor
		self.changefactor = 1.1		# used to zoom
		self.invchangefactor = 1.0/self.changefactor # used to invert zoom
		
		self.update_flag = True
		
		self.inspector = None
		
		self.model_matrix = None
		
		self.init_scale_flag = True
		
	def set_opt_scale(self):
		dims = self.drawable.get_data_dims()
		
		xscale = self.parent.width()/dims[0]
		yscale = self.parent.height()/dims[1]
		
		if yscale < xscale: xscale = yscale
		
		self.drawable.cam.scale = xscale
		
	def get_data_dims(self):
		return self.drawable.get_data_dims()
	
	def determine_dimensions(self):
		width = self.width()
		height= self.height()
		depth= self.depth()
		
		self.left = -width/2.0
		self.right = width/2.0
		self.bottom = -height/2.0
		self.top =  height/2.0
		self.near = depth/2.0
		self.far =  -depth/2.0
	
		self.update_dims = False

	def set_width(self,w):
		self.w = w
		self.drawable.resizeEvent(self.width(),self.height())
		
	def set_depth(self,d):
		self.d = d
		self.drawable.resizeEvent(self.width(),self.height())
		
	def set_height(self,h):
		self.h = h
		self.drawable.resizeEvent(self.width(),self.height())
		
	def width(self):
		try:
			return int(self.w)
		except:
			return 0
	
	def height(self):
		try:
			return int(self.h)
		except:
			return 0
	
	def depth(self):
		try:
			return int(self.d)
		except:
			return 0
	
	def set_data(self,data):
		try:
			self.drawable.set_data(data)
			self.image = image
			self.update()
		except: pass
	
	def render(self):
		self.draw()
	
	def draw(self):
		
		if self.init_scale_flag:
			self.set_opt_scale()
			self.init_scale_flag = False
		
		lighting = glIsEnabled(GL_LIGHTING)
		glEnable(GL_LIGHTING)
		
		glPushMatrix()
		self.cam.position()
		self.drawable.render()
		glPopMatrix()
		
		if not lighting: glDisable(GL_LIGHTING)
		
		
	def set_plane(self,plane):
		self.cam.set_plane(plane)
	
	def updateGL(self):
		self.parent.updateGL()
		
	def get_render_dims_at_depth(self, depth):
		return self.parent.get_render_dims_at_depth(depth)
	
	def get_near_plane_dims(self):
		return self.parent.get_near_plane_dims()
		
	def get_start_z(self):
		return self.parent.get_start_z()
	

class EM2DGLWindow(EMGLWindow):
	'''
	A class for managing a 3D object as an interactive widget
	'''
	def __init__(self,parent,gl_view):
		EMGLWindow.__init__(self,parent,gl_view)

		self.draw_frame = True
		
		self.decoration = EM3DPlainBorderDecoration(self)
		
		self.w = self.parent.width()
		self.h = self.parent.height()
		self.draw_vd_frame = False
		self.update_border_flag = False
		
	def set_draw_frame(self,bool):
		self.draw_frame = bool

	def width(self):
		return self.w
		
	def height(self):
		return self.h
	
	def set_width(self,w): self.w = w
	def set_height(self,h): self.h = h
	
	def set_frozen(self,frozen):
		#self.drawable.set_frozen(frozen)
		if frozen: self.decoration.set_color_flag(EM3DPlainBorderDecoration.FROZEN_COLOR)
		else : self.decoration.set_color_flag(EM3DPlainBorderDecoration.DEFAULT_COLOR)
	
	def draw(self):
		
		self.cam.position()
		
		self.vdtools.update(self.width()/2.0,self.height()/2.0)
		
		glEnable(GL_DEPTH_TEST)
		glPushMatrix()
		glTranslatef(-self.width()/2.0,-self.height()/2.0,0)
		self.drawable.draw()
		glPopMatrix()
		
		lighting = glIsEnabled(GL_LIGHTING)
		glEnable(GL_LIGHTING)
		
		self.decoration.draw(self.update_border_flag)
		self.update_border_flag = False
		
		if self.draw_vd_frame: self.vdtools.draw_frame()
		if not lighting: glDisable(GL_LIGHTING)

	def isinwin(self,x,y):
		interception = False
		#print "I have this many corner sets",len(self.corner_sets)
		for i,p in enumerate(self.corner_sets):
			if self.vdtools.isinwinpoints(x,y,p):
				interception = True
				#self.cam.set_plane(self.planes[i])
				#print "plane is",self.planes[i]
				self.vdtools.setModelMatrix(self.model_matrices[i])
				break
		return interception
		
	def isinwin(self,x,y):
		return self.vdtools.isinwin(x,y)
	
class EMGLView2D_v2(EMEventRerouter):
	"""
	A view of a 2D drawable type, such as a single 2D image or a matrix of 2D images
	
	"""
	def __init__(self, parent,image):
		self.parent = parent
		self.cam = Camera2(self)
		#self.cam.setCamTrans('default_z',-parent.get_depth_for_height(height_plane))
		
		if isinstance(parent,EMImageMXModule) or isinstance(parent,EMImage2DModule):
			self.drawable = parent
			self.w = self.parent.width()
			self.h = self.parent.height()
		else:
			if isinstance(image,list):
				if len(image) == 1:
					self.become_2d_image(image[0])
				else:
					self.drawable = EMImageMXModule(image)
					self.drawable.set_parent(self)
					self.w = self.parent.width()
					self.h = self.parent.height()
			elif isinstance(image,EMData):
				self.become_2d_image(image)
			else:
				print "error, the EMGLView2D class must be initialized with data"
				return
			
		self.drawable.suppressInspector = True
		
		EMEventRerouter.__init__(self,self.drawable)
		self.initflag = True
		
		self.update_flag = True
		
		
		self.sizescale = 1.0
		self.changefactor = 1.1
		self.invchangefactor = 1.0/self.changefactor

	def get_drawable(self):
		return self.drawable

	def set_sizescale(self,scale):
		self.sizescale = scale
	
	def get_drawable(self):
		return self.drawable
	
	def set_shapes(self,shapes,shrink):
		self.drawable.set_shapes(shapes,shrink)
		
	def become_2d_image(self,a):
		self.drawable = EMImage2DModule(a)
		self.drawable.set_parent(self)
		#self.drawable.originshift = False
		self.w = a.get_xsize()
		#if self.w > self.parent.width(): self.w = self.parent.width()
		self.h = a.get_ysize()
		#if self.h > self.parent.height(): self.h = self.parent.height()
		#print self.w,self.h

	def set_width(self,w,resize_event=True):
		self.w = w
		if resize_event: self.drawable.resizeEvent(self.width(),self.height())
		self.update_border_flag = True
		
	def set_height(self,h,resize_event=True):
		self.h = h
		if resize_event: self.drawable.resizeEvent(self.width(),self.height())
		self.update_border_flag = True

	def width(self):
		try:
			return int(self.sizescale*self.w)
		except:
			return 0
		#return self.drawable.width()
	
	def height(self):
		try:
			return int(self.sizescale*self.h)
		except:
			return 0

	def set_data(self,data):
		self.drawable.set_data(data)
		
	def initializeGL(self):
		self.drawable.initializeGL()
	
	def draw(self):
		
		self.drawable.render()
		


class EMGLView2D:
	"""
	A view of a 2D drawable type, such as a single 2D image or a matrix of 2D images
	
	"""
	def __init__(self, parent,image):
		self.parent = parent
		self.cam = Camera2(self)
		#self.cam.setCamTrans('default_z',-parent.get_depth_for_height(height_plane))
		
		if isinstance(parent,EMImageMXModule) or isinstance(parent,EMImage2DModule):
			self.drawable = parent
			self.w = self.parent.width()
			self.h = self.parent.height()
		else:
			if isinstance(image,list):
				if len(image) == 1:
					self.become_2d_image(image[0])
				else:
					self.drawable = EMImageMXModule(image)
					self.drawable.set_parent(self)
					self.w = self.parent.width()
					self.h = self.parent.height()
			elif isinstance(image,EMData):
				self.become_2d_image(image)
			else:
				print "error, the EMGLView2D class must be initialized with data"
				return
			
		self.drawable.suppressInspector = True
		self.initflag = True
		self.vdtools = EMViewportDepthTools(self)
		
		self.update_flag = True
		
		self.draw_frame = False
		
		self.sizescale = 1.0
		self.changefactor = 1.1
		self.invchangefactor = 1.0/self.changefactor
		
		self.inspector = None
		
		self.update_border_flag = False
		
		self.decoration = EM3DPlainBorderDecoration(self)
	
	def set_sizescale(self,scale):
		self.sizescale = scale
	
	def get_drawable(self):
		return self.drawable
	
	def set_frozen(self,frozen):
		#self.drawable.set_frozen(frozen)
		if frozen: self.decoration.set_color_flag(EM3DPlainBorderDecoration.FROZEN_COLOR)
		else : self.decoration.set_color_flag(EM3DPlainBorderDecoration.DEFAULT_COLOR)
	
	def context(self):
		# asking for the OpenGL context from the parent
		return self.parent.context()
	
	def set_shapes(self,shapes,shrink):
		self.drawable.set_shapes(shapes,shrink)
		
	def become_2d_image(self,a):
		self.drawable = EMImage2DModule(a)
		self.drawable.set_parent(self)
		#self.drawable.originshift = False
		self.w = a.get_xsize()
		if self.w > self.parent.width(): self.w = self.parent.width()
		self.h = a.get_ysize()
		if self.h > self.parent.height(): self.h = self.parent.height()
		
	def eye_coords_dif(self,x1,y1,x2,y2,mdepth=True):
		return self.vdtools.eye_coords_dif(x1,y1,x2,y2,mdepth)
	
	def set_update_P_inv(self,val=True):
		self.vdtools.set_update_P_inv(val)
	
	def set_width(self,w,resize_event=True):
		self.w = w
		if resize_event: self.drawable.resizeEvent(self.width(),self.height())
		self.update_border_flag = True
		
	def set_height(self,h,resize_event=True):
		self.h = h
		if resize_event: self.drawable.resizeEvent(self.width(),self.height())
		self.update_border_flag = True

	def width(self):
		try:
			return int(self.sizescale*self.w)
		except:
			return 0
		#return self.drawable.width()
	
	def height(self):
		try:
			return int(self.sizescale*self.h)
		except:
			return 0
		#return self.drawable.height()
	
	def get_render_area_coords(self):
		return self.vdtools.get_corners()
	
	def set_data(self,data):
		self.drawable.set_data(data)
		
	def initializeGL(self):
		self.drawable.initializeGL()
	
	def testBoundaries(self):
		'''
		Called when the image is first drawn, this resets the dimensions of this object
		if it is larger than the current size of the viewport. It's somewhat of a hack,
		but it's early stages in the design
		'''
		h = self.vdtools.getMappedHeight()
		w = self.vdtools.getMappedWidth()
		
		if ( w > viewport_width() ):
			self.w = viewport_width()/self.sizescale
		if ( h > viewport_height() ):
			self.h = viewport_height()/self.sizescale
	
	def draw(self):
		
		self.cam.position()
		
		self.vdtools.update(self.width()/2.0,self.height()/2.0)
		
		if (self.initflag == True):
			self.testBoundaries()
			self.initflag = False

		
		#self.mediator.checkBoundaryIssues()
		if (self.update_flag):
			self.drawable.resizeEvent(self.width(),self.height())
			self.update_flag = False
			
		glPushMatrix()
		glTranslatef(-self.width()/2.0,-self.height()/2.0,0)
		#try: 
		self.drawable.render()
		#except Exception, inst:
			#print type(inst)     # the exception instance
			#print inst.args      # arguments stored in .args
			#print int
		glPopMatrix()
		
		lighting = glIsEnabled(GL_LIGHTING)
		glEnable(GL_LIGHTING)
		self.decoration.draw(self.update_border_flag)
		self.update_border_flag = False
		
		if self.draw_frame: self.vdtools.draw_frame()
		if not lighting: glDisable(GL_LIGHTING)
		
	#def update(self):
		#self.parent.updateGL()
	
	#def updateGL(self):
		#self.parent.updateGL()
	
	def get_inspector(self):
		if (self.inspector == None):
			if self.drawable == None:
				return None
			self.drawable.init_inspector()
			self.drawable.inspector.show()
			self.drawable.inspector.hide()
			self.inspector = self.drawable.inspector
			
		return self.inspector
	
	def mousePressEvent(self, event):
		if event.button()==Qt.MidButton or (event.button()==Qt.LeftButton and self.inspector == None):
			try:	
				self.drawable.init_inspector()
				self.drawable.inspector.show()
				self.drawable.inspector.hide()
				self.parent.addQtWidgetDrawer(self.get_inspector())
			except:
				pass
				#print "the mouse event had no effect"
			
		if event.modifiers() == Qt.ControlModifier:
			self.cam.mousePressEvent(event)
		else:
			l=self.vdtools.mouseinwin(event.x(),viewport_height()-event.y(),self.width(),self.height())
			qme=QtGui.QMouseEvent(event.type(),QtCore.QPoint(l[0],l[1]),event.button(),event.buttons(),event.modifiers())
			self.drawable.mousePressEvent(qme)
			#self.drawable.mousePressEvent(event)
		
		#self.updateGL()
	
	def scaleEvent(self,delta):
		if ( delta > 0 ):
			self.sizescale *= self.changefactor
		elif ( delta < 0 ):
			self.sizescale *= self.invchangefactor

		self.drawable.resizeEvent(self.width(),self.height())
	
	def wheelEvent(self,event):
		if event.modifiers() == Qt.ControlModifier:
			self.scaleEvent(event.delta())
			self.decoration.set_force_update()
			#print "updating",self.drawWidth(),self.drawHeight()
		else:
			self.drawable.wheelEvent(event)
			
		#self.updateGL()
	
	def mouseMoveEvent(self,event):
		if event.modifiers() == Qt.ControlModifier:
			self.cam.mouseMoveEvent(event)
		else:
			l=self.vdtools.mouseinwin(event.x(),viewport_height()-event.y(),self.width(),self.height())
			qme=QtGui.QMouseEvent(event.type(),QtCore.QPoint(l[0],l[1]),event.button(),event.buttons(),event.modifiers())
			self.drawable.mouseMoveEvent(qme)
			#self.drawable.mouseMoveEvent(event)
		
		#self.updateGL()

	def mouseReleaseEvent(self,event):
		if event.modifiers() == Qt.ControlModifier:
			self.cam.mouseReleaseEvent(event)
		else:
			l=self.vdtools.mouseinwin(event.x(),viewport_height()-event.y(),self.width(),self.height())
			qme=QtGui.QMouseEvent(event.type(),QtCore.QPoint(l[0],l[1]),event.button(),event.buttons(),event.modifiers())
			self.drawable.mouseReleaseEvent(qme)
			#self.drawable.mouseReleaseEvent(event)

	def mouseDoubleClickEvent(self, event):
		if event.modifiers() == Qt.ControlModifier:
			self.cam.mouseMoveEvent(event)
		else:
			l=self.vdtools.mouseinwin(event.x(),viewport_height()-event.y(),self.width(),self.height())
			qme=QtGui.QMouseEvent(event.type(),QtCore.QPoint(l[0],l[1]),event.button(),event.buttons(),event.modifiers())
			self.drawable.mouseDoubleClickEvent(qme)
			
	def keyPressEvent(self,event):
		self.drawable.keyPressEvent(event)

		#self.updateGL()
	def emit(self, *args,**kargs):
		self.parent.emit(*args,**kargs)

	def leaveEvent(self):
		self.drawable.leaveEvent()
	
	def toolTipEvent(self,event):
		pass
	
	def isinwin(self,x,y):
		return self.vdtools.isinwin(x,y)


class EMGLViewQtWidget:
	def __init__(self, parent=None, qwidget=None, widget_parent=None):
		self.parent = parent
		self.qwidget = qwidget
		self.draw_frame = False
		self.mapcoords = True
		self.itex = 0
		self.gen_texture = True
		self.click_debug = False
		self.cam = Camera2(self)
		#self.cam.setCamTrans('default_z',-parent.get_depth_for_height(height_plane))
		self.cam.motion_rotate(0,0)
		self.borderwidth = 3.0
		self.glbasicobjects = EMBasicOpenGLObjects()
		self.setQtWidget(qwidget)
		self.P_inv = None
		self.childreceiver = None
		self.widget_parent = widget_parent
		
		self.current = None
		self.previous = None
		
		self.e2children = []
		self.is_child = False
		
		self.vdtools = EMViewportDepthTools(self)
	
		self.refresh_dl = True
		self.texture_dl = 0
		self.texture_lock = 0
		self.decoration = EM3DPlainBorderDecoration(self)
	def get_decoration(self):
		return self.decoration
	
	def lock_texture(self):
		self.texture_lock += 1
	
	def unlock_texture(self):
		self.texture_lock -= 1
	
	def __del__(self):
		if (self.itex != 0 ):
			self.parent.deleteTexture(self.itex)
		
	def set_update_P_inv(self,val=True):
		self.vdtools.set_update_P_inv(val)
	
	def set_refresh_dl(self,val=True):
		self.refresh_dl = val
	
	def width(self):
		return self.qwidget.width()
	
	def height(self):
		return self.qwidget.height()

	def setQtWidget(self, widget, delete_current = False):
		if ( delete_current and self.qwidget != None ):
			self.qwidget.deleteLater()
		
		self.qwidget = widget
		
		if ( widget != None ):
			#self.qwidget.setVisible(True)
			self.qwidget.setEnabled(True)
			self.gen_texture = True
			self.updateTexture()
			
	def updateTexture(self,force=False):
		if self.texture_lock > 0:
			return
		
		if ( self.itex == 0 or self.gen_texture == True or force) : 
			self.refresh_dl = True
			if (self.itex != 0 ):
				#passpyth
				self.parent.deleteTexture(self.itex)
			self.gen_texture = False
			##print "binding texture"
			#self.qwidget.setVisible(True)
			#self.qwidget.repaint()
			if ( isinstance(self.qwidget,QtGui.QWidget) ):
				pixmap = QtGui.QPixmap.grabWidget(self.qwidget)
			else:
				pixmap = QtGui.QPixmap.grabWidget(self.qwidget.widget())
			#self.qwidget.setVisible(False)
			if (pixmap.isNull() == True ): print 'error, the pixmap was null'
			self.itex = self.parent.bindTexture(pixmap)
			if ( self.itex == 0 ): print 'Error - I could not generate the texture'
		
			#self.decoration.set_force_update()
	
	def draw(self):
		#print "draw children"
		if (self.qwidget == None or self.itex == 0) :
			#print "no widget - draw children return" 
			return
		
		#self.cam.debug = True
		self.cam.position()
		
		# make sure the vdtools store the current matrices
		self.vdtools.update(self.width()/2.0,self.height()/2.0)
		
		if self.refresh_dl == True:
			if self.texture_dl != 0:
				glDeleteLists(self.texture_dl,1)
			
			self.texture_dl=glGenLists(1)
			glNewList(self.texture_dl,GL_COMPILE)
		
			glPushMatrix()
			glEnable(GL_TEXTURE_2D)
			glBindTexture(GL_TEXTURE_2D,self.itex)
			glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
			glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
			glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE)
			glBegin(GL_QUADS)
			glTexCoord2f(0.,0.)
			glVertex(-self.qwidget.width()/2.0,-self.qwidget.height()/2.0)
			glTexCoord2f(1.,0.)
			glVertex( self.qwidget.width()/2.0,-self.qwidget.height()/2.0)
			glTexCoord2f(1.,1.)
			glVertex( self.qwidget.width()/2.0, self.qwidget.height()/2.0)
			glTexCoord2f(0.,1.)
			glVertex( -self.qwidget.width()/2.0,self.qwidget.height()/2.0)
			glEnd()
			glDisable(GL_TEXTURE_2D)
			glPopMatrix()
			
			glEndList()
		
		if self.texture_dl == 0: return
		glCallList(self.texture_dl)
		self.refresh_dl = False
	
		lighting = glIsEnabled(GL_LIGHTING)
		glEnable(GL_LIGHTING)
		self.decoration.draw()
		if self.draw_frame:
			self.vdtools.draw_frame()
			#except Exception, inst:
				#print type(inst)     # the exception instance
				#print inst.args      # arguments stored in .args
				#print int
		if (not lighting): glDisable(GL_LIGHTING)
		
		# now draw children if necessary - such as a qcombobox list view that has poppud up
		for i in self.e2children:
			glPushMatrix()
			#try:
			i.draw()
			#except Exception, inst:
				#print type(inst)     # the exception instance
				#print inst.args      # arguments stored in .args
				#print int
			glPopMatrix()

	def isinwin(self,x,y):
		for i in self.e2children:
			if i.isinwin(x,y):
				self.childreceiver = i
				return True
		
		return self.vdtools.isinwin(x,y)
	
	def eye_coords_dif(self,x1,y1,x2,y2):
		return self.vdtools.eye_coords_dif(x1,y1,x2,y2)
			
	def mouseinwin(self,x,y,width,height):
		return self.vdtools.mouseinwin(x,y,width,height)

	def toolTipEvent(self,event):
		if ( self.childreceiver != None ):
			# this means this class already knows that the mouse event is in the child
			# that is being displayed
			self.childreceiver.toolTip(event)
			self.childreceiver = None
			return
		l=self.mouseinwin(event.x(),viewport_height()-event.y(),self.width(),self.height())
		cw=self.qwidget.childAt(l[0],l[1])
		if cw == None: 
			QtGui.QToolTip.hideText()
			self.gen_texture = True
			self.updateTexture()
			return
	
		p1 = QtCore.QPoint(event.x(),event.y())
		p2 = self.parent.mapToGlobal(p1)
		QtGui.QToolTip.showText(p2,cw.toolTip())
	
	def wheelEvent(self,event):
		doElse = True
		try:
			if event.modifiers() == Qt.ControlModifier:
				self.cam.wheelEvent(event)
				doElse = False
		except: pass
		if doElse:
			if ( self.childreceiver != None ):
				# this means this class already knows that the mouse event is in the child
				# that is being displayed
				self.childreceiver.wheelEvent(event)
				self.childreceiver = None
				return
			else:
				# if we have any children (i.e. a drop down combo box) it should now disappear
				if len(self.e2children) > 0:
					self.e2children.pop()
					return
			l=self.mouseinwin(event.x(),viewport_height()-event.y(),self.width(),self.height())
			cw=self.qwidget.childAt(l[0],l[1])
			if cw == None: return
			gp=self.qwidget.mapToGlobal(QtCore.QPoint(l[0],l[1]))
			lp=cw.mapFromGlobal(gp)
			qme=QtGui.QWheelEvent(lp,event.delta(),event.buttons(),event.modifiers(),event.orientation())
			QtCore.QCoreApplication.sendEvent(cw,qme)
			self.gen_texture = True
			self.updateTexture()
	
	def mouseDoubleClickEvent(self, event):
		if ( self.childreceiver != None ):
			# this means this class already knows that the mouse event is in the child
			# that is being displayed
			self.childreceiver.mouseDoubleClickEvent(event)
			self.childreceiver = None
			return
		l=self.mouseinwin(event.x(),viewport_height()-event.y(),self.width(),self.height())
		cw=self.qwidget.childAt(l[0],l[1])
		if cw == None: return
		gp=self.qwidget.mapToGlobal(QtCore.QPoint(l[0],l[1]))
		lp=cw.mapFromGlobal(gp)
		if (isinstance(cw,QtGui.QComboBox)):
			print "it's a combo"
		else:
			qme=QtGui.QMouseEvent(event.type(),lp,event.button(),event.buttons(),event.modifiers())
			#self.qwidget.setVisible(True)
			QtCore.QCoreApplication.sendEvent(cw,qme)
			#self.qwidget.setVisible(False)
		self.gen_texture = True
		self.updateTexture()
		
	def get_depth_for_height(self,height_plane):
		return self.parent.get_depth_for_height(height_plane)
	
	def mousePressEvent(self, event):
		if event.modifiers() == Qt.ControlModifier:
			self.cam.mousePressEvent(event)
		else:
			if ( self.childreceiver != None ):
				# this means this class already knows that the mouse event is in the child
				# that is being displayed
				self.childreceiver.mousePressEvent(event)
				self.childreceiver = None
				return
			else:
				# if we have any children (i.e. a drop down combo box) it should now disappear
				if len(self.e2children) > 0:
					self.e2children.pop()
					return
				
			l=self.mouseinwin(event.x(),viewport_height()-event.y(),self.width(),self.height())
			cw=self.qwidget.childAt(l[0],l[1])
			if cw == None: return
			##print cw.objectName()
			gp=self.qwidget.mapToGlobal(QtCore.QPoint(l[0],l[1]))
			lp=cw.mapFromGlobal(gp)
			if (isinstance(cw,QtGui.QComboBox)):
				cw.showPopup()
				cw.hidePopup()
				widget = EMGLViewQtWidget(self.parent,None,cw);
				widget.setQtWidget(cw.view())
				widget.cam.loadIdentity()	
				widget.cam.setCamTrans("x",cw.geometry().x()-self.width()/2.0+cw.view().width()/2.0)
				widget.cam.setCamTrans("y",((self.height()/2.0-cw.geometry().y())-cw.view().height()/2.0))
				widget.cam.setCamTrans("z",0.1)
				widget.draw_frame = False
				self.e2children.append(widget)
				self.e2children[0].is_child = True
			else:
				qme=QtGui.QMouseEvent( event.type(),lp,event.button(),event.buttons(),event.modifiers())
				if (self.is_child): QtCore.QCoreApplication.sendEvent(self.qwidget,qme)
				else: QtCore.QCoreApplication.sendEvent(cw,qme)
				
			self.gen_texture = True
			self.updateTexture()
		
	def mouseMoveEvent(self,event):
		#pos = QtGui.QCursor.pos()
		
		#print event.x(),pos.x(),event.y(),pos.y()
		#pos = self.parent.mapFromGlobal(pos)
		#print "now",event.x(),pos.x(),event.y(),pos.y()
	
		if event.modifiers() == Qt.ControlModifier:
			self.cam.mouseMoveEvent(event)
		else:
			if ( self.childreceiver != None ):
				# this means this class already knows that the mouse event is in the child
				# that is being displayed
				self.childreceiver.mouseMoveEvent(event)
				self.childreceiver = None
				return
			else:
				l=self.mouseinwin(event.x(),viewport_height()-event.y(),self.width(),self.height())
				cw=self.qwidget.childAt(l[0],l[1])
				self.current = cw
				if ( self.current != self.previous ):
					QtGui.QToolTip.hideText()
					if ( self.current != None ):
						qme=QtCore.QEvent(QtCore.QEvent.Enter)
						QtCore.QCoreApplication.sendEvent(self.current,qme)
						
					if ( self.previous != None ):
						qme=QtCore.QEvent(QtCore.QEvent.Leave)
						QtCore.QCoreApplication.sendEvent(self.previous,qme)
				
				self.previous = self.current
				if cw == None:
					QtGui.QToolTip.hideText()
					if ( self.previous != None ):
						qme=QtCore.QEvent(QtCore.QEvent.Leave)
						QtCore.QCoreApplication.sendEvent(self.previous,qme)
						self.gen_texture = True
						self.updateTexture()
					return
				gp=self.qwidget.mapToGlobal(QtCore.QPoint(l[0],l[1]))
				lp=cw.mapFromGlobal(gp)
				qme=QtGui.QMouseEvent(event.type(),lp,event.button(),event.buttons(),event.modifiers())
				QtCore.QCoreApplication.sendEvent(cw,qme)
			# FIXME
			# setting the gen_texture flag true here causes the texture to be regenerated
			# when the mouse moves over it, which is inefficient.
			# The fix is to only set the gen_texture flag when mouse movement
			# actually causes a change in the appearance of the widget (for instance, list boxes from comboboxes)
			self.gen_texture = True
			self.updateTexture()
	
	def keyPressEvent(self,event):
		print "we are in this place"
		if ( self.childreceiver != None ):
			# this means this class already knows that the mouse event is in the child
			# that is being displayed
			self.childreceiver.keyPressEvent(event)
			self.childreceiver = None
			return
		else:
			pos = self.parent.mapFromGlobal(QtGui.QCursor.pos())
			l=self.mouseinwin(pos.x(),viewport_height()-pos.y(),self.width(),self.height())
			cw=self.qwidget.childAt(l[0],l[1])
			self.current = cw
			print cw
			if ( self.current != self.previous ):
				QtGui.QToolTip.hideText()
				if ( self.current != None ):
					qme=QtCore.QEvent(QtCore.QEvent.Enter)
					QtCore.QCoreApplication.sendEvent(self.current,qme)
					
				if ( self.previous != None ):
					qme=QtCore.QEvent(QtCore.QEvent.Leave)
					QtCore.QCoreApplication.sendEvent(self.previous,qme)
			
			self.previous = self.current
			if cw == None:
				QtGui.QToolTip.hideText()
				if ( self.previous != None ):
					qme=QtCore.QEvent(QtCore.QEvent.Leave)
					QtCore.QCoreApplication.sendEvent(self.previous,qme)
					self.gen_texture = True
					self.updateTexture()
				return
			#gp=self.qwidget.mapToGlobal(QtCore.QPoint(l[0],l[1]))
			#lp=cw.mapFromGlobal(gp)
			#qme=QtGui.QMouseEvent(event.type(),lp,event.button(),event.buttons(),event.modifiers())
			#print cw,event
			QtCore.QCoreApplication.sendEvent(cw,event)
		# FIXME
		# setting the gen_texture flag true here causes the texture to be regenerated
		# when the mouse moves over it, which is inefficient.
		# The fix is to only set the gen_texture flag when mouse movement
		# actually causes a change in the appearance of the widget (for instance, list boxes from comboboxes)
		self.gen_texture = True
		self.updateTexture()

	def mouseReleaseEvent(self,event):
		if event.modifiers() == Qt.ControlModifier:
			self.cam.mouseReleaseEvent(event)
		else:
			if ( self.childreceiver != None ):
				# this means this class already knows that the mouse event is in the child
				# that is being displayed
				#try:
				self.childreceiver.mouseReleaseEvent(event)
				self.childreceiver = None
				self.e2children.pop()
				return
			else:
				# if we have any children (i.e. a drop down combo box) it should now disappear
				if len(self.e2children) > 0:
					self.e2children.pop()
					return
			
			l=self.mouseinwin(event.x(),viewport_height()-event.y(),self.width(),self.height())
			cw=self.qwidget.childAt(l[0],l[1])
			if cw == None: return
			gp=self.qwidget.mapToGlobal(QtCore.QPoint(l[0],l[1]))
			lp=cw.mapFromGlobal(gp)
			if (isinstance(cw,QtGui.QComboBox)):
				print "it's a combo"
				#cw.showPopup()
			else:
				qme=QtGui.QMouseEvent(event.type(),lp,event.button(),event.buttons(),event.modifiers())
				if (self.is_child):
					##print self.qwidget
					##print self.qwidget.currentIndex().row()
					##print self.widget_parent
					##print self.qwidget.rect().left(),self.qwidget.rect().right(),self.qwidget.rect().top(),self.qwidget.rect().bottom()
					##print lp.x(),lp.y()
					self.widget_parent.setCurrentIndex(self.qwidget.currentIndex().row())
					#self.widget_parent.changeEvent(QtCore.QEvent())
					#self.widget_parent.highlighted(self.qwidget.currentIndex().row())
					#self.qwidget.commitData(self.qwidget.parent())
					##print self.qwidget.currentText()
					#self.widget_parent.setVisible(True)
					#self.widget_parent.setEnabled(True)
					#self.qwidget.setVisible(True)
					#QtCore.QCoreApplication.sendEvent(self.widget_parent,qme)
					#self.qwidget.setVisible(False)
					self.widget_parent.emit(QtCore.SIGNAL("activated(QString)"),self.widget_parent.itemText(self.qwidget.currentIndex().row()))
				else:
					#self.qwidget.setVisible(True)
					QtCore.QCoreApplication.sendEvent(cw,qme)
					#self.qwidget.setVisible(False)
			
			self.gen_texture = True
			self.updateTexture()
		
	def leaveEvent(self):
		if (self.current != None) : 
			qme = QtCore.QEvent(QtCore.QEvent.Leave)
			QtCore.QCoreApplication.sendEvent(self.current,qme)
			self.current = None
			self.previouse = None
			self.gen_texture = True
			self.updateTexture()
			
	def enterEvent():
		pass
	def timerEvent(self,event=None):
		pass
		#self.cam.motion_rotate(.2,.2)


class EMFloatingWidgets(QtOpenGL.QGLWidget):
	def __init__(self,parent=None):
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True)
		# enable multisampling to combat aliasing
		fmt.setSampleBuffers(True)
		# stenciling is for object dependent shading
		fmt.setStencil(True)
		QtOpenGL.QGLWidget.__init__(self,fmt, parent)
		
		self.setMouseTracking(True)
		self.fov = 2*180*atan2(1,5)/pi
		self.zNear= 1000
		self.zFar = 3000
		self.floatwidget = EMFloatingWidgetsCore(self)
		
		self.cam = Camera()
		
		self.animatables = [] # an array of animatable objects - must have the animate(time) function
		
		self.timer = QTimer()
		QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.timeout)
		self.timer.start(10)
		
	def get_depth_for_height(self, height):
		# This function returns the width and height of the renderable 
		# area at the origin of the data volume
		depth = height/(2.0*tan(self.fov/2.0*pi/180.0))
		return depth
	
	def timeout(self):
		
		if len(self.animatables) == 0: return
		
		for i,animatable in enumerate(self.animatables):
			if not animatable.animate(time.time()):
				# this could be dangerous
				self.animatables.pop(i)
		
		self.updateGL()
	
	def register_animatable(self,animatable):
		self.animatables.append(animatable)
	
	def get_render_dims_at_depth(self, depth):
		# This function returns the width and height of the renderable 
		# area at the origin of the data volume
		height = -2*tan(self.fov/2.0*pi/180.0)*(depth)
		width = self.aspect*height
		return [width,height]

	
	def initializeGL(self):
		#print "initializeGL"
		glClearColor(0,0,0,0)
		
		glLightfv(GL_LIGHT0, GL_AMBIENT, [0.1, 0.1, 0.1, 1.0])
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_POSITION, [0.1,.1,1.,0.])
	
		glEnable(GL_LIGHTING)
		glEnable(GL_LIGHT0)
		glEnable(GL_DEPTH_TEST)
		
		glEnable(GL_NORMALIZE)
		
		glClearStencil(0)
		glEnable(GL_STENCIL_TEST)
		
		glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST)
		glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST)
		glHint(GL_TEXTURE_COMPRESSION_HINT, GL_NICEST)
	
		# enable multisampling to combat aliasing
		if ( "GL_ARB_multisample" in glGetString(GL_EXTENSIONS) ): glEnable(GL_MULTISAMPLE)
		else: glDisable(GL_MULTISAMPLE)
		
	def draw(self):
		#print "draw"
		glClear(GL_COLOR_BUFFER_BIT)
		if glIsEnabled(GL_DEPTH_TEST):
			glClear(GL_DEPTH_BUFFER_BIT)
		if glIsEnabled(GL_STENCIL_TEST):
			glClear(GL_STENCIL_BUFFER_BIT)
			
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
	
		self.floatwidget.render()
		
	def resizeGL(self, width, height):
		#print "resizeGL"
		glViewport(0,0,self.width(),self.height())
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		#glFrustum(-1.*width/height,1.*width/height, -1.,1., 5.,15.)
		
		# fov angle is the given by
		#self.fov = 2*180*atan2(1,5)/pi
		# aspect ratio is given by
		self.aspect = float(self.width())/float(self.height())
		# this is the same as the glFrustum call above
		depth = self.get_depth_for_height(height_plane)
		
		#gluPerspective(self.fov,self.aspect,depth-depth/4,depth+depth/4)
		gluPerspective(self.fov,self.aspect,self.zNear,self.zFar)
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		try: self.floatwidget.resizeEvent(width,height)
		except: print "couldn't resize floatwidget"
		
	def mousePressEvent(self, event):
		self.floatwidget.mousePressEvent(event)
	
	def mouseMoveEvent(self, event):
		self.floatwidget.mouseMoveEvent(event)
		
	def mouseReleaseEvent(self, event):
		self.floatwidget.mouseReleaseEvent(event)

	def mouseDoubleClickEvent(self, event):
		self.floatwidget.mouseDoubleClickEvent(event)

	def wheelEvent(self, event):
		self.floatwidget.wheelEvent(event)

	def toolTipEvent(self, event):
		self.floatwidget.toolTipEvent(event)
		QtGui.QToolTip.hideText()

	def dragMoveEvent(self,event):
		print "received drag move event, but I don't do anything about it :("
		
	def event(self,event):
		if event.type() == QtCore.QEvent.MouseButtonPress: 
			self.mousePressEvent(event)
			return True
		elif event.type() == QtCore.QEvent.MouseButtonRelease:
			self.mouseReleaseEvent(event)
			return True
		elif event.type() == QtCore.QEvent.MouseMove: 
			self.mouseMoveEvent(event)
			return True
		elif event.type() == QtCore.QEvent.MouseButtonDblClick: 
			self.mouseDoubleClickEvent(event)
			return True
		elif event.type() == QtCore.QEvent.Wheel: 
			self.wheelEvent(event)
			return True
		elif event.type() == QtCore.QEvent.ToolTip: 
			self.toolTipEvent(event)
			return True
		else: 
			return QtOpenGL.QGLWidget.event(self,event)

	def hoverEvent(self,event):
		self.floatwidget.hoverEvent(event)

	def get_near_plane_dims(self):
		height = 2.0*self.zNear * tan(self.fov/2.0*pi/180.0)
		width = self.aspect * height
		return [width,height]
	
	def get_start_z(self):
		return self.zNear

class EMFloatingWidgetsCore:
	"""A QT widget for rendering EMData objects. It can display single 2D or 3D images 
	or sets of 2D images.
	"""
	def __init__(self, parent=None):
		#print "init"
		self.parent = parent
	
		self.imtex=0
		self.current = None
		self.previous = None
	
		self.initFlag = True
		self.qwidgets = []
		
		self.suppressUpdateGL = False
		self.rotor = None
		#print "init done"
	
	def emit(self,signal,event,a=None,b=None):
		pass
	
	def register_animatable(self,animatable):
		self.parent.register_animatable(animatable)
		
	def get_depth_for_height(self, height):
		try: 
			return self.parent.get_depth_for_height(height)
		except:
			print "parent can't get height for depth"
			return 0

	def height(self):
		return self.parent.height()
	
	def width(self):
		return self.parent.width()

	def updateGL(self):
		if not self.suppressUpdateGL:
			try: self.parent.updateGL()
			except: pass

	def addQtWidgetDrawer(self,widget):
		w = EMGLViewQtWidget(self)
		w.setQtWidget(widget)
		self.qwidgets.append(w)
		
		#print "initializeGL done"
	def render(self):
		if ( self.initFlag == True ):
			self.initFlag = False
			self.fd = QtGui.QFileDialog(self.parent,"Open File",QtCore.QDir.currentPath(),QtCore.QString("Image files (*.img *.hed *.mrc)"))
			QtCore.QObject.connect(self.fd, QtCore.SIGNAL("finished(int)"), self.finished)
			self.fd.show()
			self.fd.hide()
			#self.qwidgets.append(EMGLViewQtWidget(self.parent))
			#self.qwidgets[0].setQtWidget(self.fd)
			##self.qwidgets[0].cam.set_cam_x(-100)
			
			
			#self.qwidgets.append(EMGLViewQtWidget(self.parent))
			#self.qwidgets[1].setQtWidget(self.fd2)
			#self.qwidgets[1].cam.set_cam_x(-200)
			
			a = EMGLViewQtWidget(self.parent)
			a.setQtWidget(self.fd)
			rotor = EMGLRotorWidget(self)
			rotor.add_widget(a)
			#rotor.add_widget(d)
			
			
			w = EMGLView2D(self,test_image(0,size=(256,256)))
			rotor.add_widget(w)
			insp = w.get_inspector()
			e = EMGLViewQtWidget(self.parent)
			e.setQtWidget(insp)
			rotor.add_widget(e)
			
			w2 = EMGLView2D(self,[test_image(0,size=(512,512)),test_image(1,size=(512,512))]*4)
			rotor.add_widget(w2)
			insp2 = w2.get_inspector()
			f = EMGLViewQtWidget(self.parent)
			f.setQtWidget(insp2)
			rotor.add_widget(f)
			
			w3 = EMGLView3D(self,test_image_3d(3))
			ww = EM3DGLWindow(self.parent,w3)
			rotor.add_widget(ww)
			insp3 = w3.get_inspector()
			g = EMGLViewQtWidget(self.parent)
			g.setQtWidget(insp3)
			rotor.add_widget(g)
		
			self.qwidgets.append(EM3DGLWindow(self,rotor))
			
			#rotor2 = EMGLRotorWidget(self,45,50,15,EMGLRotorWidget.BOTTOM_ROTARY,60)
			#rotor2.add_widget(a)
			#rotor2.add_widget(e)
			#rotor2.add_widget(g)
			#rotor2.add_widget(w2)
			
			#self.qwidgets.append(EM3DGLWindow(self,rotor2))
			
			#rotor3 = EMGLRotorWidget(self,-15,50,-15,EMGLRotorWidget.TOP_ROTARY,60)
			#rotor3.add_widget(w)
			#rotor3.add_widget(e)
			#rotor3.add_widget(f)
			##rotor3.add_widget(w2)
			
			#self.qwidgets.append(EM3DGLWindow(self,rotor3))
			
			#rotor4 = EMGLRotorWidget(self,-25,10,40,EMGLRotorWidget.LEFT_ROTARY)
			##rotor4.add_widget(w2)
			#rotor4.add_widget(e)
			#rotor4.add_widget(g)
			##rotor4.add_widget(ww)
			
			#self.qwidgets.append(EM3DGLWindow(self,rotor4))
			
		glPushMatrix()
		glTranslate(-100,0,-1250)
		for i in self.qwidgets:
			#print "getting opengl matrices"
			glPushMatrix()
			i.draw()
			glPopMatrix()
			#print "paint child done"
		
		
		glPopMatrix()
		#print "draw done"
	def finished(self,val):
		if ( val == 1 ):
			for i in self.fd.selectedFiles():
				a=EMData.read_images(str(i))
				if len(a) == 1:
					a = a[0]
					if a.get_zsize() != 1:
						w = EMGLView3D(self,a)
						self.qwidgets.append(w)
					else:
						w = EMGLView2D(self,a)
						self.qwidgets.append(w)
				else:
					w = EMGLView2D(self,a)
					self.qwidgets.append(w)
					
	def timer(self):
		pass
		#self.updateGL()
		
	def bindTexture(self,pixmap):
		return self.parent.bindTexture(pixmap)
	
	def deleteTexture(self,val):
		return self.parent.deleteTexture(val)
	
	def get_render_dims_at_depth(self, depth):
		try: return self.parent.get_render_dims_at_depth(depth)
		except:
			print "parent can't get render dims at for depth"
			return

	def resizeEvent(self, width, height):
		for i in self.qwidgets:
			i.set_update_P_inv()
	
	def mousePressEvent(self, event):
		for i in self.qwidgets:
			if ( i.isinwin(event.x(),self.height()-event.y()) ):
				i.mousePressEvent(event)
				intercepted = True
				self.updateGL()
				return
	
	def mouseMoveEvent(self, event):
		for i in self.qwidgets:
			if ( i.isinwin(event.x(),self.height()-event.y()) ):
				self.current = i
				if (self.current != self.previous ):
					if ( self.previous != None ):
						self.previous.leaveEvent()
				i.mouseMoveEvent(event)
				self.previous = i
				self.updateGL()
				return
		
	def mouseReleaseEvent(self, event):
		for i in self.qwidgets:
			if ( i.isinwin(event.x(),self.height()-event.y()) ):
				i.mouseReleaseEvent(event)
				self.updateGL()
				return
					
		
	def mouseDoubleClickEvent(self, event):
		for i in self.qwidgets:
			if ( i.isinwin(event.x(),self.height()-event.y()) ):
				i.mouseDoubleClickEvent(event)
				self.updateGL()
				return
		
		
	def wheelEvent(self, event):
		for i in self.qwidgets:
				if ( i.isinwin(event.x(),self.height()-event.y()) ):
					i.wheelEvent(event)
					self.updateGL()
					return

	def toolTipEvent(self, event):
		for i in self.qwidgets:
			if ( i.isinwin(event.x(),self.height()-event.y()) ):
				i.toolTipEvent(event)
				self.updateGL()
				return
		
		QtGui.QToolTip.hideText()
		

	def dragMoveEvent(self,event):
		print "received drag move event"
		
	def event(self,event):
		#print "event"
		#QtGui.QToolTip.hideText()
		if event.type() == QtCore.QEvent.MouseButtonPress: 
			self.mousePressEvent(event)
			return True
		elif event.type() == QtCore.QEvent.MouseButtonRelease:
			self.mouseReleaseEvent(event)
			return True
		elif event.type() == QtCore.QEvent.MouseMove: 
			self.mouseMoveEvent(event)
			return True
		elif event.type() == QtCore.QEvent.MouseButtonDblClick: 
			self.mouseDoubleClickEvent(event)
			return True
		elif event.type() == QtCore.QEvent.Wheel: 
			self.wheelEvent(event)
			return True
		elif event.type() == QtCore.QEvent.ToolTip: 
			self.toolTipEvent(event)
			return True
		else: 
			return QtOpenGL.QGLWidget.event(self,event)

	def hoverEvent(self,event):
		print "hoverEvent"
		if self.inspector :
			for i in self.qwidgets:
				if ( i.isinwin(event.x(),self.height()-event.y()) ):
					i.hoverEvent(event)
					break
		self.updateGL()
	
	def get_near_plane_dims(self):
		return self.parent.get_near_plane_dims()
	
	def get_start_z(self):
		return self.parent.get_start_z()

# This is just for testing, of course
if __name__ == '__main__':
	app = QtGui.QApplication(sys.argv)
	window = EMFloatingWidgets()
	window2 = EMParentWin(window)
	window2.show()
	
	sys.exit(app.exec_())
