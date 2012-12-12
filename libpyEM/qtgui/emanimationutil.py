#!/usr/bin/env python

#
# Author: David Woolford, 10/10/2008 (woolford@bcm.edu)
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

from time import time

from PyQt4 import QtCore
from math import sin,acos
from EMAN2 import Vec3f, Transform

class Animator:
	'''
	Register Animatables with this class
	'''
	def __init__(self):
		self.animatables = []
		self.time = -1
		self.begin_time = time()
		
		self.timer_enabled = False
		self.timer_interval =  15
		
	def time_out(self):
		self.time  = time()
		
		if len(self.animatables) != 0:
			rm = []
			for a,i in enumerate(self.animatables):
				if not i.animate(self.time): rm.append(a)
				
			rm.reverse()
			for a in rm:
				self.animatables.pop(a)			
		
			self.update()
		else:
			if not QtCore.QObject.disconnect(self.timer, QtCore.SIGNAL("timeout()"), self.time_out):
				print "failed to disconnect timer"
			
			self.timer_enabled = False
		
	def get_time(self):
		return self.time
	
	def update(self): pass
	
	def register_animatable(self,animatable):
		if not self.timer_enabled:
			self.__enable_timer()
		self.animatables.append(animatable)
	
	
	def __enable_timer(self):
		if self.timer_enabled == False:
			self.timer = QtCore.QTimer()
			QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.time_out)
			
			self.timer.start(self.timer_interval)
			self.timer_enabled = True
		else: print "timer already enabled in Animator"
		
	def animation_done_event(self,animated): raise NotImplementedError("Inheriting classes should define this function")

class Animatable:
	cache_dts = None
	def __init__(self):
		self.time = 0		# time indicates the current time used for the basis of animation.
		self.time_interval = 0.3 # 0.3 seconds for the animation to complete
		self.inverse_time_inverval = 1.0/self.time_interval
		self.time_begin = 0 # records the time at which the animation was begun
		self.animated = True
		self.n = 100
		if Animatable.cache_dts == None:
			self.init_cache_dts()
		
	def init_cache_dts(self):
		from math import tanh,sin,pi
		Animatable.cache_dts = []
		tanh_approach = False
		linear_approach = True
		for i in range(self.n):
			if tanh_approach:
				val = (1+ (tanh(-4+float(i)/(self.n-1)*8)))/2.0
				Animatable.cache_dts.append(val)
			elif linear_approach:
				#  Linear
				Animatable.cache_dts.append(float(i)/(self.n-1))
			else:
				# sine approach
				Animatable.cache_dts.append(sin(float(i)/(self.n-1)*pi/2))
		#print Animatable.cache_dts
		
	def set_animated(self,val=True):
		self.animated = val
		
	def cancel_animation(self):
		self.animation = False
		
	def is_animated(self): return self.animated
	
	def animate(self,t):
		if not self.animated: return False
		
		if self.time_begin == 0:
			self.time_begin = t
			
		self.time = t -self.time_begin

		if self.time > self.time_interval:
			self.time_begin = 0
			self.time = 0
			self.animated = False
			self.target.animation_done_event(self)
			return 0
		else:
			dt = self.get_dt()
			self.calculate_animation(dt)
			return 1
	
	def calculate_animation(self,dt): raise
	
	def get_dt(self):
		idx = int( self.n*self.time*self.inverse_time_inverval)
		#print Animatable.cache_dts
		return Animatable.cache_dts[idx]
		
		
	def transform(self): raise
	
class SingleAxisRotationAnimation(Animatable):
	def __init__(self,target,start,end,axis=[0,1,0]):
		Animatable.__init__(self)
		self.target = target
		self.start = start
		self.current = start
		self.end = end
		self.axis = axis
	
	def __del__(self):
		self.target.animation_done_event(self)
	
	def get_start(self):
		return self.start
	
	def get_end(self):
		return self.end
	
	def get_current(self):
		return self.current
	
	def calculate_animation(self,dt):
		'''
		based on the assumption that dt goes from 0 to 1
		'''
		if dt > 1: raise
		self.current = (1-dt)*self.start + dt*self.end
		self.target.set_rotation(self.current)
	#def transform(self):
		#glRotate(self.current,*self.axis)
		
class XYScaleAnimation(Animatable):
	def __init__(self,target,start,end):
		Animatable.__init__(self)
		self.start = start
		self.current = start
		self.end = end
		self.target = target
	
	def __del__(self):
		self.target.animation_done_event(self)
	
	
	def get_start(self):
		return self.start
	
	def get_end(self):
		return self.end
	
	def get_current(self):
		return self.current
	
	def calculate_animation(self,dt):
		'''
		based on the assumption that dt goes from 0 to 1
		'''
		if dt > 1: raise
		self.current = (1-dt)*self.start + dt*self.end
		self.target.set_xy_scale(self.current)
	def transform(self):
		glScale(self.current,self.current,1.0)
		
class SingleValueIncrementAnimation(Animatable):
	def __init__(self,target,start,end):
		Animatable.__init__(self)
		self.start = start
		self.current = start
		self.end = end
		self.target = target
	
	def __del__(self):
		self.target.animation_done_event(self)
	
	def get_start(self):
		return self.start
	
	def get_end(self):
		return self.end
	
	def get_current(self):
		return self.current
	
	def calculate_animation(self,dt):
		'''
		based on the assumption that dt goes from 0 to 1
		'''
		if dt > 1: raise
		self.current = (1-dt)*self.start + dt*self.end
		self.target.set_animation_increment(self.current)
		
	def transform(self):pass

class TranslationAnimation(Animatable):
	'''
	Animates from 3D oint one point to another
	'''
	def __init__(self,target,start,end):
		Animatable.__init__(self)
		self.start = start
		self.current = start
		self.end = end
		self.target = target
	
	def __del__(self):
		self.target.animation_done_event(self)
	
	def get_start(self):
		return self.start
	
	def get_end(self):
		return self.end
	
	def get_current(self):
		return self.current
	
	def calculate_animation(self,dt):
		'''
		based on the assumption that dt goes from 0 to 1
		'''
		if dt > 1: raise
		x = (1-dt)*self.start[0] + dt*self.end[0]
		y = (1-dt)*self.start[1] + dt*self.end[1]
		z = (1-dt)*self.start[2] + dt*self.end[2]
		self.target.set_translation((x,y,z))
		
	def transform(self):pass
		
	
class LineAnimation(Animatable):
	'''
	Animates from one point to another
	Should really be 2D line animation
	'''
	def __init__(self,target,start,end):
		Animatable.__init__(self)
		self.start = start
		self.current = start
		self.end = end
		self.target = target
	
	def __del__(self):
		self.target.animation_done_event(self)
	
	def get_start(self):
		return self.start
	
	def get_end(self):
		return self.end
	
	def set_end(self,end):
		self.end = end
		self.start = self.current
		self.time = 0
		self.time_begin = 0
	
	def get_current(self):
		return self.current
	
	def calculate_animation(self,dt):
		'''
		based on the assumption that dt goes from 0 to 1
		'''
		if dt > 1: raise
		x = (1-dt)*self.start[0] + dt*self.end[0]
		y = (1-dt)*self.start[1] + dt*self.end[1]
		self.current = (x,y)
		self.target.set_line_animation(x,y)
		
	def transform(self):pass
	
	
class OrientationListAnimation(Animatable):
	'''
	Animates through a list of orientations using great circles on the unit sphere
	'''
	def __init__(self,target,transforms,radius=1):
		if len(transforms) < 2: raise RuntimeError("Can run orientation list animation on a list less than length 2")
		Animatable.__init__(self)
		self.target = target
		self.transitions = len(transforms)-1
		self.transforms = transforms
		self.stage = 0
		self.radius = radius
		self.z_point = Vec3f(0,0,1)
		self.arc_points = []
		self.omega = 1
		self.rot_vec = None
		self.__establish_sphere_points(self.stage)
		self.arc_points[0].append(self.radius*self.start_point)

	def __establish_sphere_points(self,stage):
		self.start_point = self.transforms[stage].transpose()*self.z_point
		self.end_point = self.transforms[stage+1].transpose()*self.z_point
		if self.start_point != self.end_point:
			self.angle = acos(self.start_point.dot(self.end_point))
			self.omega = 1
		else:
			self.angle = 0
			normal = Vec3f(-self.start_point[2],0,-self.start_point[0])
			normal.normalize()
			
			T = Transform()
			d = {"type":"spin"}
			d["omega"] = self.omega
			d["n1"] = normal[0]
			d["n2"] = normal[1]
			d["n3"] = normal[2]
			T.set_rotation(d)
			self.rot_vec = self.start_point
			p = self.start_point
			self.start_point = T*self.rot_vec
			self.end_point = T*self.rot_vec
			self.omega += 1
			
		self.sinangle = sin(self.angle)
		self.arc_points.append([])
	def __del__(self):
		self.target.animation_done_event(self)
	
	def animate(self,t):
		if not self.animated: return False
		
		if self.time_begin == 0:
			self.time_begin = t
			
		self.time = t -self.time_begin

		if self.time > self.time_interval:
			if self.arc_points[self.stage][-1] != self.end_point:
				self.arc_points[self.stage].append(self.radius*self.end_point)
				self.target.arc_animation_update(self.arc_points)
			self.time_begin = 0
			self.time = 0
			self.stage += 1
			if self.stage < self.transitions:
				self.__establish_sphere_points(self.stage)
				return 1
			else:
				
				self.animated = False
				self.target.animation_done_event(self)
				return 0
		else:
			dt = self.get_dt()
			self.calculate_animation(dt)
			return 1
	
	def calculate_animation(self,dt):
		'''
		based on the assumption that dt goes from 0 to 1
		'''
		if dt > 1: raise
		if self.angle != 0:
			t1 = dt*self.angle
			t2 = self.angle - t1
			next = (sin(t2)*self.start_point + sin(t1)*self.end_point)/self.sinangle
		else:
			T = Transform()
			d = {"type":"spin"}
			d["omega"] = dt*360
			d["n1"] = self.rot_vec[0]
			d["n2"] = self.rot_vec[1]
			d["n3"] = self.rot_vec[2]
			T.set_rotation(d)
			
			next = T*self.start_point
			
			#next = Vec3f(self.start_point[0],self.start_point[1],self.start_point[2])
		
		self.arc_points[self.stage].append(self.radius*next)
		self.target.arc_animation_update(self.arc_points)
		
	def transform(self):pass
		