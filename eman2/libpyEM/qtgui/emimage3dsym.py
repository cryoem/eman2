#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
# and David Woolford 10/26/2007 (woolford@bcm.edu)
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

from EMAN2 import *
from OpenGL import GL, GLU, GLUT
from OpenGL.GL import *
from OpenGL.GLU import *
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from emapplication import EMApp, get_application
from emglobjects import EM3DModel, EMGLWidget, Camera, Camera2, EMViewportDepthTools2, EMGLProjectionViewMatrices, get_default_gl_colors
from emimageutil import EMTransformPanel
from math import *
from time import *
from valslider import ValSlider
from weakref import WeakKeyDictionary
import weakref

MAG_INCREMENT_FACTOR = 1.1

class Orientations:
	def __init__(self):
		self.prop = 5.0
		self.sym = "c1"
		self.sym_object=parsesym(self.sym)
		self.nomirror = True
		self.angle_label = "delta"
		self.strategy = "eman"
	
	def get_prop(self): return self.prop
	def mirror_enabled(self): return not self.nomirror
	def get_sym(self): return self.sym
	def set_prop(self,val): self.prop = float(val)
	def set_mirror(self,val): self.nomirror = not val
	def set_angle_label(self,val): self.angle_label = str(val)
	def set_strategy(self,val): self.strategy = str(val)
	def set_symmetry(self,val):
		self.sym = str(val)
		self.sym_object= Symmetries.get(self.sym)

class ColumnGraphics:
	def __init__(self):
		self.small_column_color = "dark grey"
		self.tall_column_color = "white"
		self.max_score = 0.0
		self.min_score = 0.0
		self.interval = 0.0
		if not hasattr(self,"colors"): self.colors = get_default_gl_colors()
		
		self.set_mixed_color()
		self.column_scores = None
		
	def set_column_scores(self,scores):
		
		self.column_scores = scores
		if scores != None:
			self.max_score = max(scores)
			self.min_score = min(scores)
			self.interval = self.max_score - self.min_score
		else:
			self.max_score = 0.0
			self.min_score = 0.0
			self.interval = 0.0
		
	def set_mixed_color(self):
		self.mixed_color = MixedColor(self.colors[self.small_column_color], self.colors[self.tall_column_color])
		
	def set_max_score(self,score):
		self.max_score = score
		self.interval = self.max_score - self.min_score
	def set_min_score(self,score):
		self.min_score = score
		self.interval = self.max_score - self.min_score
		
	def set_small_column_color(self,color):
		self.small_column_color = color
		self.set_mixed_color()
	
	def set_tall_column_color(self,color):
		self.tall_column_color = color
		self.set_mixed_color()
	
	def load_mixed_gl_color(self,value):
		if self.interval == 0:
			frac = 1.0
		else:
			frac = (value-self.min_score)/self.interval
		
		self.mixed_color.load_gl_color(frac)
		
	def load_basic_gl_color(self):
		self.mixed_color.load_color_2()
	
class MixedColor:
	def __init__(self,color1,color2):
		self.color1 = color1 # the color of the smallest column, for example
		self.color2 = color2 # the color of the largeset column, for example
		
		self.a1 = self.color1["ambient"]
		self.a2 = self.color2["ambient"]
		
		self.d1 = self.color1["diffuse"]
		self.d2 = self.color2["diffuse"]
		
		self.s1 = self.color1["specular"]
		self.s2 = self.color2["specular"]
		
		self.e1 = self.color1["emission"]
		self.e2 = self.color2["emission"]
		
		self.h1 = self.color1["shininess"]
		self.h2 = self.color2["shininess"]
		
	def load_gl_color(self,frac):
		'''
		frac should be in [0,1]
		'''
		ifrac = 1-frac
		r = range(4)
		r2 = range(3)
		ambient = [ ifrac*self.a1[i]+frac*self.a2[i] for i in r]
		diffuse = [ ifrac*self.d1[i]+frac*self.d2[i] for i in r]
		specular = [ ifrac*self.s1[i]+frac*self.s2[i] for i in r]
		emission = [ ifrac*self.e1[i]+frac*self.e2[i] for i in r2]
		shininess = ifrac*self.h1+frac*self.h2
		
		glColor(ambient)
		glMaterial(GL_FRONT,GL_AMBIENT,ambient)
		glMaterial(GL_FRONT,GL_DIFFUSE,diffuse)
		glMaterial(GL_FRONT,GL_SPECULAR,specular)
		glMaterial(GL_FRONT,GL_EMISSION,emission)
		glMaterial(GL_FRONT,GL_SHININESS,shininess)
	
	def load_color_2(self):
		glColor(self.a2)
		glMaterial(GL_FRONT,GL_AMBIENT,self.a2)
		glMaterial(GL_FRONT,GL_DIFFUSE,self.d2)
		glMaterial(GL_FRONT,GL_SPECULAR,self.s2)
		glMaterial(GL_FRONT,GL_EMISSION,self.e2)
		glMaterial(GL_FRONT,GL_SHININESS,self.h2)
	

class EulerData:
	'''
	A mixin for the EM3DSymModel - takes care of everything that needs to occur
	if you have supplied a list of EMData objects - the EMData's must have the
	'xform.projection' header attribute
	'''
	def __init__(self):
		self.data = None
		self.score_options = []
		self.eulers = []
		
	def __getitem__(self,idx):
		return self.data[idx]
	
	def __len__(self): return len(self.data)
	
	def set_data(self,data):
		self.data = data
		self.eulers = []
		for i in xrange(len(self.data)):
			if hasattr(self.data,"get_image_header"):
				d = self.data.get_image_header(i)
			else:
				d = self.data[i]
			try:
				self.eulers.append(d["xform.projection"])
			except:
				self.eulers.append(Transform())
#				raise RuntimeError("The data must all have the xform.projection header attribute")
		
		# get the first header and get any items that can be cast to a float
		header = self.data[0].get_attr_dict()
		self.score_options = []
		for key,value in header.items():
			try:
				# make sure we can cast the value to a float, that way in can become the height of the cylinder
				float(value)
				self.score_options.append(key)
			except: pass
		
		self.score_options.append("None")
		self.score_options.sort()
		
	def get_eulers(self): return self.eulers
	
	def get_score_options(self): return self.score_options
	
	def get_score_list(self,key,normalize=True,log_scale=False):
		l = []
		for i in range(len(self.data)):
			if hasattr(self.data,"get_image_header"):
				d = self.data.get_image_header(i)
			else:
				d = self.data[i]
			 
			if log_scale:
				try: l.append(log(float(d[key])))
				except: l.append(0)
			else:
				try: l.append(float(d[key]))
				except: l.append(0)
		if normalize: return self.normalize_float_list(l)
		else: return l
	
	def normalize_float_list(self,l):
		mn = min(l)
		mx = max(l)
		diff = float(mx-mn)
		if diff != 0: n = [ (val-mn)/(diff) for val in l ]
		else: n = [1 for val in l]
		return n
	
	
class EM3DSymModel(EM3DModel,Orientations,ColumnGraphics):
	def __init__(self, gl_widget, eulerfilename=None):
		self.arc_anim_dl = None
		self.arc_anim_points = None
		self.window_title = "EulerXplor"
		EM3DModel.__init__(self, gl_widget)
		Orientations.__init__(self)
		ColumnGraphics.__init__(self)
		
		self.eulerfilename = eulerfilename # the file name that the Euler data is read from
		self.vdtools = EMViewportDepthTools2(gl_widget)
		self.eulers = [] # will eventually store Transform objects
		#self.points = [] # will eventually store the points on the asymmetric unit
		self.point_colors = [] # will eventually store colors for the different points
		self.cam = Camera2(self) # Stores the camera position
		self.cam.allow_phi_rotations = False

		self.cube = False
		self.inspector=None
		self.mouse_mode = None # In future we might facilitate changing mouse modes
		
		self.sphere_points_dl = 0 # this will be a display list of points or cylinders on the unit sphere
		self.spheredl = 0 # This is a display list of the low resolution sphere
		self.highresspheredl = 0 # This is a display list of a high resolution sphere
		self.diskdl = 0 # This is a display list for a disk
		self.cappedcylinderdl = 0 # This is a display list for a capped cylinder (the cap is only at the outer end)
		self.great_arc_dl = 0 # This will be a display list for great arcs 
		self.asym_u_triangles_dl = 0  # This will be a display list for triangles in the asymmetric units
		self.cylinderdl = 0 # This will be a display list for a cylinder

		self.gq= None #  a glu quadric
		
		self.width_scale = 1.0 # a with scale factor for points on the unit sphere
		self.height_scale = 1.0 # a height scale factor for points on the unit sphere
		self.arc_width_scale = 0.2 # The width of the great arcs 
		self.force_update = True  # Force update everything - causes a couple of dispay lists to be regenerated
#		
		self.display_euler = True # Display sphere points flag
		self.display_tri = False # Display asymm unit triangles flag
		self.display_arc = True # Display great arcs flag
		self.display_all_syms = False
#		self.sym_object = None
		self.update_sphere_points_dl = True # A way to force the update of the display list for the points on the unit sphere
		self.flatten = 0.0	# animates transition from sphere to plane, 0-1, 0 = sphere
		
		self.radius = 50 # the radius of the sphere - important display parameter
		self.arc_segments = 16 # number of discrete segments in every great arc 

		self.eulers_specified = False # from e2au.py = if you want to specify a specific set of eulers this is required
		self.specified_eulers = None # from e2au.py = if you want to specify a specific set of eulers this is required
		self.colors_specified = False # as above
		self.specified_colors = None # as above
		
		self.file = None # For tracing - hasn't been tested in a while
		self.lr = -1 # For tracing - hasn't been tested in a while
		self.hr = -1 # For tracing - hasn't been tested in a while
		self.tracedata = [] # For tracing - hasn't been tested in a while
		self.trace_dl = 0 # For tracing - hasn't been tested in a while
		self.reduce = False  # For tracing - hasn't been tested in a while
		
		self.colors = get_default_gl_colors() # A map for OpenGL (light) colors
		self.arc_color = "emerald" # The color of the arcs, can be changed using the inspector
		
		self.euler_data = None # will be turned into an EulerData if necessary
		self.image_display_window = None # used if people are looking at image data, for example, if they opened a list of EMData's in e2.py using euler_display
		self.displayed_image_number = None # used to update the displayed image if new data is set 
		self.special_euler = None # If set, the special Euler will get the special color - should be an index (int)
		self.special_color = "gold"
		
		self.log_scale = False # Stores whether or not the cylinder height (score) is shown in log scale
		self.current_cylinder_score = None # stores which header attribute is being used to scale the heights of the cylinders
	
		self.mirror_eulers = False # If True the drawn Eulers are are also rendered on the opposite side of the sphere - see make_sym_dl_list. e2eulerxplor turns this to True. Also useful for verifying that we have accurately demarcated the non mirror redundant portion of thhe asymmetric unit

		self.initialized = True
	def __del__(self):
		#self.clear_gl_memory() # this is intentionally commented out, it makes sense to clear the memory but not here
		if self.image_display_window != None:
			w = self.image_display_window 
			self.image_display_window = None
			w.close()
			
		self.close_image_display()
	def close_image_display(self):
		if self.image_display_window:
			tmp = self.image_display_window
			self.image_display_window = None
			tmp.close()
			
	def clear_gl_memory(self):
		try: self.get_gl_widget().makeCurrent() # this is important  when you have more than one OpenGL context operating at the same time
		except: return # maybe the object is already dead in which case OpenGL probably cleaned the memory up anyway
		vals = ["sphere_points_dl","spheredl","highresspheredl","diskdl","cappedcylinderdl","great_arc_dl","asym_u_triangles_dl","cylinderdl","trace_dl" ]
		for dl in vals :
			if getattr(self,dl) != 0:
				glDeleteLists(getattr(self,dl),1)
				setattr(self,dl,0)
	
	def set_data(self,emdata_list):
		return self.set_emdata_list_as_data(emdata_list)
				
	def set_emdata_list_as_data(self,emdata_list,default="None"):
		self.euler_data = EulerData()
		self.euler_data.set_data(emdata_list)
		self.specify_eulers(self.euler_data.get_eulers())
		if default != "None":
			self.current_cylinder_score = default
			l = self.euler_data.get_score_list(default,log_scale=self.log_scale)
			self.set_column_scores(l)
		else: self.set_column_scores(None)
		self.get_inspector().set_score_options(self.euler_data.get_score_options(),default)
		self.force_update = True
		if self.image_display_window and self.displayed_image_number and len(self.euler_data) > self.displayed_image_number:
			self.image_display_window.set_data(self.euler_data[self.displayed_image_number],"Data")
			self.image_display_window.updateGL()
		else: self.displayed_image_number = None # blanket response
		
	def set_column_score_key(self,key):
		if key == "None":
			self.set_column_scores(None)
		else:
			l = self.euler_data.get_score_list(key,log_scale=self.log_scale)
			self.current_cylinder_score = key
			self.set_column_scores(l)
			
#		self.regen_dl()
	
	def set_log_scale(self,val):
		self.log_scale = val
		if str(self.current_cylinder_score) != "None":
			l = self.euler_data.get_score_list(self.current_cylinder_score,log_scale=self.log_scale)
			self.set_column_scores(l)
	
	def set_width_scale(self,val):
		self.width_scale = val
		self.regen_dl()
		
	def set_height_scale(self,val):
		self.height_scale = val
		self.regen_dl()
		
	def set_arc_width_scale(self,val):
		self.arc_width_scale = val
		self.regen_dl()
		
	def set_arc_segments(self,val):
		self.arc_segments = val
		self.regen_dl()
		
	def set_arc_color(self,s): self.arc_color = s
	
	def keyPressEvent(self,event):
		
		if event.key() == Qt.Key_F1:
			self.display_web_help("http://blake.bcm.edu/emanwiki/EMAN2/Programs/e2eulerxplor")
		elif event.key() == Qt.Key_F :
			if self.flatten>0 : self.flatten=0.0
			else: self.flatten=1.0
			self.generate_current_display_list(True)
			self.updateGL()
		else:
			EM3DModel.keyPressEvent(self,event)
			
	def object_picked(self,object_number):
		resize_necessary = False
		if self.image_display_window == None:
			from emimage2d import EMImage2DWidget
			self.image_display_window = EMImage2DWidget()
			QtCore.QObject.connect(self.image_display_window,QtCore.SIGNAL("module_closed"),self.on_image_display_window_closed)
			resize_necessary = True
				
		self.image_display_window.set_data(self.euler_data[object_number],"Data")
		self.displayed_image_number = object_number
		if resize_necessary:
			get_application().show_specific(self.image_display_window)
			self.image_display_window.optimally_resize()
		else:
			self.image_display_window.updateGL()
			
		if object_number != self.special_euler:
			self.special_euler = object_number
			self.regen_dl()
	
	def on_image_display_window_closed(self):
		self.image_display_window = None
	
	def mouseReleaseEvent(self,event):
		
		if self.euler_data != None: #and (event.modifiers()&Qt.ShiftModifier or self.mouse_mode == "pick"):
			
			v = self.vdtools.wview.tolist()
				
			glSelectBuffer(512)
			glRenderMode(GL_SELECT)
			glInitNames()
			glMatrixMode(GL_PROJECTION)
			glPushMatrix()
			glLoadIdentity()
			gluPickMatrix(event.x(),v[-1]-event.y(),2,2,v)
			self.get_gl_widget().load_perspective()
			glMatrixMode(GL_MODELVIEW)
			glInitNames()
			self.render()
			glMatrixMode(GL_PROJECTION)
			glPopMatrix()
			glMatrixMode(GL_MODELVIEW)
			glFlush()
		
			hits = list(glRenderMode(GL_RENDER))
			for hit in hits:
				a,b,c=hit
				if len(c) > 0:
					self.object_picked(int(c[0]-1))
					break
#			else:
#				if self.special_euler != None:
#					self.special_euler = None
#					self.regen_dl()
			#return
		
		
		EM3DModel.mouseReleaseEvent(self,event)
	
	def get_type(self):
		return "Symmetry Viewer"
	
	def update_data(self,data):
		pass

	def eye_coords_dif(self,x1,y1,x2,y2,mdepth=True):
		return self.vdtools.eye_coords_dif(x1,y1,x2,y2,mdepth)

	def set_radius(self,radius):
		if ( radius > 0 ):
			self.radius = radius
			self.force_update = True
		else:
			print "Error, tried to set a zero or negative radius (",radius,")"
			exit(1)
	
	def trace_great_triangles(self,inc_mirror):
		triangles = self.sym_object.get_asym_unit_triangles(inc_mirror)
		if ( self.asym_u_triangles_dl != 0 ): glDeleteLists(self.asym_u_triangles_dl, 1)
		
		self.asym_u_triangles_dl=glGenLists(1)
		
		if self.sym_object.get_name() == "d" and self.sym_object.get_nsym()% 4 != 0:
			if inc_mirror:
				pass
#				p1 = (3*triangles[0][2]+triangles[0][1])/4.0
#				p2 = (triangles[0][2]+3*triangles[0][1])/4.0
#				a = p1.normalize()
#				a = p2.normalize()
#				
#				t = []
#				t.append([triangles[0][0],triangles[0][1],p2])
#				t.append([triangles[0][0],p2, p1])
#				t.append([triangles[0][0],p1, triangles[0][2]])
#				triangles = t
			else:
				p = (triangles[0][2]+triangles[0][1])/2.0
				a = p.normalize()
				t = []
				t.append([triangles[0][0],triangles[0][1],p])
				t.append([triangles[0][0],p, triangles[0][2]])
				triangles = t
	
		glNewList(self.asym_u_triangles_dl,GL_COMPILE)
		
		glPushMatrix()
		glScalef(self.radius,self.radius,self.radius)
		if ( self.sym_object.is_h_sym() ):
			pass
		else:
			glBegin(GL_TRIANGLES)
			l = len(triangles)
			for i,t in enumerate(triangles):
				p1 = t[0]-t[1]
				p2 = t[1]-t[2]
				n = p2.cross(p1)
				n.normalize()
				glNormal(n[0],n[1],n[2])
				
				for j,p in enumerate(t):
					self.loadcolor(i,j,l)
					glVertex(p[0],p[1],p[2])
			glEnd()
		glPopMatrix()
		glEndList()
		
	def loadcolor(self,i,j,l):
		if ( l == 1):
			if (j == 0):
				self.load_gl_color("gold")
			elif (j == 1):
				self.load_gl_color("ruby")
			elif ( j == 2):
				self.load_gl_color("emerald")
		elif ( l == 2):
			if i == 0:
				if (j == 0):
					self.load_gl_color("gold")
				elif (j == 1):
					self.load_gl_color("ruby")
				elif ( j == 2):
					self.load_gl_color("emerald")
			elif i == 1:
				if (j == 0):
					self.load_gl_color("gold")
				elif (j == 1):
					self.load_gl_color("emerald")
				elif ( j == 2):
					self.load_gl_color("ruby")
		else:
			self.load_gl_color("gold")
				
	
	def load_gl_color(self,name):
		color = self.colors[name]
		glColor(color["ambient"])
		glMaterial(GL_FRONT,GL_AMBIENT,color["ambient"])
		glMaterial(GL_FRONT,GL_DIFFUSE,color["diffuse"])
		glMaterial(GL_FRONT,GL_SPECULAR,color["specular"])
		glMaterial(GL_FRONT,GL_EMISSION,color["emission"])
		glMaterial(GL_FRONT,GL_SHININESS,color["shininess"])
	
	def gl_color(self,color):
		glColor(*color)
		glMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,color)
		glMaterial(GL_FRONT,GL_SPECULAR,[0,0,0,1])
		glMaterial(GL_FRONT,GL_EMISSION,(0,0,0,1.0))
		glMaterial(GL_FRONT,GL_SHININESS,64.0)

	def trace_great_arcs(self, points):
		
		self.init_basic_shapes()
	
		if ( self.great_arc_dl != 0 ): glDeleteLists(self.great_arc_dl, 1)
		n = len(points)
		if ( n <= 1 ):
			self.great_arc_dl = 0
			return
		
		self.great_arc_dl=glGenLists(1)
		
		glNewList(self.great_arc_dl,GL_COMPILE)
		
		if ( self.sym_object.is_h_sym() ):
			glPushMatrix()
			dz = self.sym_object.get_params()["equator_range"]
			if ( dz == 0 ): dz = 5
			if (not self.nomirror):
				glTranslatef(0,0,dz)
			self.generate_arcs(points,n,1)
			glPopMatrix()
			
			glPushMatrix()
			glTranslatef(0,0,-dz)
			self.generate_arcs(points,n,1)
			glPopMatrix()
		else:
			self.generate_arcs(points,n)
		glEndList()
		
	def generate_arcs(self,points,n,halt=0):
		for i in range(0,n-halt):
			p1 = points[i]
			if ( i == n-1 ): p2 = points[0]
			else: p2 = points[i+1]
			angle = acos(p2.dot(p1))
			sinangle = sin(angle)
			prev = self.radius*Vec3f(p1[0],p1[1],p1[2])
			for t in range(1, self.arc_segments+1):
				timeangle = float(t)/float(self.arc_segments)*angle
				p1Copy = self.radius*Vec3f(p1[0],p1[1],p1[2]*(1.0-self.flatten))
				p2Copy = self.radius*Vec3f(p2[0],p2[1],p2[2]*(1.0-self.flatten))
				next = (sin(angle-timeangle)*p1Copy + sin(timeangle)*p2Copy)/sinangle
				
				self.cylinder_to_from(next,prev,self.arc_width_scale)
				prev = Vec3f(next[0],next[1],next[2])
				
	
	def cylinder_to_from(self,next,prev,scale=0.5):
		dx = next[0] - prev[0]
		dy = next[1] - prev[1]
		dz = next[2] - prev[2]
		try:
			length = sqrt(dx**2 + dy**2 + dz**2)
		except: return
		
		if length == 0: return
		
		alt = acos(dz/length)*180.0/pi
		phi = atan2(dy,dx)*180.0/pi
		
		glPushMatrix()
		glTranslatef(prev[0],prev[1],prev[2] )
		#print "positioned at", prev[0],prev[1],prev[2]
		glRotatef(90+phi,0,0,1)
		glRotatef(alt,1,0,0)
		
		glScalef(scale,scale,length)
		glCallList(self.cylinderdl)
		glPopMatrix()

	def specify_colors(self,colors):
		self.colors_specified = True
		self.specified_colors = colors

	def specify_eulers(self,eulers):
		self.eulers_specified = True
		self.specified_eulers = eulers
	
	def init_basic_shapes(self):
		self.get_gl_widget().makeCurrent()
		if self.gq == None:
			
			self.gq=gluNewQuadric() # a quadric for general use
			gluQuadricDrawStyle(self.gq,GLU_FILL)
			gluQuadricNormals(self.gq,GLU_SMOOTH)
			gluQuadricOrientation(self.gq,GLU_OUTSIDE)
			gluQuadricTexture(self.gq,GL_FALSE)
		
		if ( self.cylinderdl == 0 ):
			self.cylinderdl=glGenLists(1)
				
			glNewList(self.cylinderdl,GL_COMPILE)
			glPushMatrix()
			gluCylinder(self.gq,1.0,1.0,1.0,12,2)
			glPopMatrix()
				
			glEndList()
					
		if self.diskdl == 0:
			self.diskdl=glGenLists(1)
				
			glNewList(self.diskdl,GL_COMPILE)
			gluDisk(self.gq,0,1,12,2)
			glEndList()
		
		if self.spheredl == 0:
			self.spheredl=glGenLists(1)
				
			glNewList(self.spheredl,GL_COMPILE)
			gluSphere(self.gq,.5,4,2)
			glEndList()

		
		if self.highresspheredl == 0:
			self.highresspheredl=glGenLists(1)
				
			glNewList(self.highresspheredl,GL_COMPILE)
			gluSphere(self.gq,.5,16,16)
			glEndList()
			
		if ( self.cappedcylinderdl == 0 ):
			self.cappedcylinderdl=glGenLists(1)
			glNewList(self.cappedcylinderdl,GL_COMPILE)
			glCallList(self.cylinderdl)
			glTranslate(0,0,1)
			glCallList(self.diskdl)
				
			glEndList()
		
#	def closeEvent(self,event):
#		print "EM3DSymModel.closeEvent()!"
#		self.close_image_display()

# Added capibility to read Euler angle files, Ryan style	
	def generate_current_display_list(self,force=False): 
		self.init_basic_shapes()
		
		l = dump_orientgens_list()
		parms = []
		li = l[str(self.strategy)]
		for i in range(1,len(li),3):
			parms.append(li[i])
		
		if (self.angle_label) not in parms:
			return
		else:
			og = self.strategy + ":" + self.angle_label + "=" + str(self.prop)
		
		if self.nomirror == True : val = 0
		else: val = 1
		self.trace_great_arcs(self.sym_object.get_asym_unit_points(val))
		self.trace_great_triangles(val)
		if ('inc_mirror' in parms):
			og += ":inc_mirror=" + str(val)
					
		[og_name,og_args] = parsemodopt(str(og))
		
		if self.eulerfilename:
			filedebug = True
		else:
			filedebug = False
			
		if ( not filedebug and not self.eulers_specified):
			eulers = self.sym_object.gen_orientations(og_name, og_args)
		elif self.eulers_specified:
			eulers = self.specified_eulers
		else:
			f = file(str(self.eulerfilename))
			lines=f.readlines()
			angles=[]
			eulers = []
			for line in lines:
				angles = str.split(line)
				alt = angles[2]
				az = angles[3]
				#print alt, az
				eulers.append(Transform({"type":"eman","az":float(az),"alt":float(alt)}))
				
		self.eulers = eulers
		self.make_sym_dl_list(eulers)
	
	def make_sym_dl_list(self,eulers):
		'''
		@param eulers - a list of Transform objects storing orientations on the unit sphere
		'''
		
		if self.sphere_points_dl > 0:
			glDeleteLists(self.sphere_points_dl,1)
		self.sphere_points_dl = glGenLists(1)
		glNewList(self.sphere_points_dl,GL_COMPILE)
		for i,p in enumerate(eulers):
			if i == self.special_euler:
				self.load_gl_color(self.special_color)
			else:
				if self.column_scores != None:
					self.load_mixed_gl_color(self.column_scores[i])
				else:
					self.load_basic_gl_color()
			
			if self.flatten>0.0 :
				glPushMatrix()
				d = p.get_rotation("eman")
				# this is the tranpose, which is then converted from a right hand coordinate system (EMAN2) into a left hand coordinate system (OpenGL)
				# hence there are no negative signs
				glRotate(d["az"],0,0,1)
				glTranslate(0,-(1.0-cos(d["alt"]*pi/180.0))*self.radius,0)

				if self.column_scores != None:
					glScalef(self.width_scale*(d["alt"]+5.0)/90.0,self.width_scale*(d["alt"]+5.0)/90.0,self.height_scale*self.column_scores[i])
				else: 
					glScalef(self.width_scale*(d["alt"]+5.0)/90.0,self.width_scale*(d["alt"]+5.0)/90.0,self.height_scale)

				glPushName(i+1)
				glCallList(self.cappedcylinderdl)
				glPopName()
				glPopMatrix()
			else:
				glPushMatrix()
				d = p.get_rotation("eman")
				# this is the tranpose, which is then converted from a right hand coordinate system (EMAN2) into a left hand coordinate system (OpenGL)
				# hence there are no negative signs
				glRotate(d["az"],0,0,1)
				glRotate(d["alt"],1,0,0)
				glRotate(d["phi"],0,0,1)
				glTranslate(0,0,self.radius)
				if self.column_scores != None:
					glScalef(self.width_scale,self.width_scale,self.height_scale*self.column_scores[i])
				else: 
					glScalef(self.width_scale,self.width_scale,self.height_scale)

				glPushName(i+1)
				glCallList(self.cappedcylinderdl)
				glPopName()
				glPopMatrix()
			
			if self.mirror_eulers and not self.nomirror:
				s = p.negate()

				if self.flatten >0 :
					glPushMatrix()
					d = s.get_rotation("eman")
					glRotate(d["az"],0,0,1)
					glTranslate(0,(1.0-d["alt"]*2.0/pi)*self.radius,0)

					if self.column_scores != None:
						glScalef(self.width_scale*(d["alt"]+5.0)/90.0,self.width_scale*(d["alt"]+5.0)/90.0,self.height_scale*self.column_scores[i])
					else: 
						glScalef(self.width_scale*(d["alt"]+5.0)/90.0,self.width_scale*(d["alt"]+5.0)/90.0,self.height_scale)
		
					glPushName(i+1)
					glCallList(self.cappedcylinderdl)
					glPopName()
					glPopMatrix()
				else :
					glPushMatrix()
					d = s.get_rotation("eman")
					# this is the tranpose, which is then converted from a right hand coordinate system (EMAN2) into a left hand coordinate system (OpenGL)
					# hence there are no negative signs
					glRotate(d["az"],0,0,1)
					glRotate(d["alt"],1,0,0)
					glRotate(d["phi"],0,0,1)
					glTranslate(0,0,self.radius)
					if self.column_scores != None:
						glScalef(self.width_scale,self.width_scale,self.height_scale*self.column_scores[i])
					else: 
						glScalef(self.width_scale,self.width_scale,self.height_scale)
		
					glPushName(i+1)
					glCallList(self.cappedcylinderdl)
					glPopName()
					glPopMatrix()
				
			
		glEndList()
		

	def angular_deviation(self,t1,t2):

		v1 = Vec3f([0,0,1]*t1)
		v2 = Vec3f([0,0,1]*t2)
		t = v2.dot(v1)
		#print t
		if t > 1: 
			if t > 1.1:
				print "error, the precision is a problem, are things normalized?"
				exit(1)
			t = 1
		if t < -1:
			if t < -1.1:
				print "error, the precision is a problem, are things normalized?"
				exit(1)
			t = -1
					
		angle = acos(t)*180/pi
		
		return angle
		
	def trace_update(self,f,lr,hr):
		if not hasattr(self,"trace_file_name"):
			self.trace_file_name = None
		if f != self.trace_file_name:
			
			print "parsing file",f,self.file
			try:
				ff=file(f,'r')
			except:
				print 'couldnt read',f 
				return
			lines=ff.readlines()
			self.tracedata = []
			for line in lines:
				s = str.split(str.strip(line))
				n = len(self.tracedata)
				if s[1] == '********':
					self.tracedata.append([])
				elif s[1] == '->':
					idx = str.find(s[3],',')
					alt = float(s[3][1:idx])
					rem = s[3][idx+1:len(s[3])-1]
					idx = str.find(rem,',')
					val = rem[0:idx]
					if val[idx-1] == ')':
						val = rem[0:idx-1]
					az = float(val)
					#print [alt,az,0]
					self.tracedata[n-1].append([alt,az,0])
			
			self.file = ff
			self.trace_file_name = f
		else:
			print "that file is already loaded"
		
		self.lr = lr
		self.hr = hr
		
		
		if self.reduce:
			print "reducing"
			for k in range(lr,hr):
				particle = self.tracedata[k]
				for orient in particle:
					t = Transform({"type":"eman","az":orient[1],"alt":orient[0],"phi":orient[2]})
					t = self.sym_object.reduce(t,0)
					d = t.get_rotation("eman")
					orient[1] = d["az"]
					orient[0] = d["alt"]
					orient[2] = d["phi"]
					
			syms = self.sym_object.get_syms()
			for k in range(lr,hr):
				particle = self.tracedata[k]
				n = len(particle)
				for i in range(1,n):
					o1 = particle[i-1]
					o2 = particle[i]
					t1 = Transform({"type":"eman","az":o1[1],"alt":o1[0],"phi":o1[2]})
					t2 = Transform({"type":"eman","az":o2[1],"alt":o2[0],"phi":o2[2]})
				
					angle = self.angular_deviation(t1,t2)
					
					for t in syms:
						t2 = Transform({"type":"eman","az":o2[1],"alt":o2[0],"phi":o2[2]})*t
						
						tmp = self.angular_deviation(t1,t2)
						
						if tmp < angle:
							angle = tmp
							
							d = t2.get_rotation("eman")
							particle[i][1] = d["az"]
							particle[i][0] = d["alt"]
							particle[i][2] = d["phi"]
						
		
		self.trace_dl = glGenLists(1)
		
		glNewList( self.trace_dl,GL_COMPILE)
		for i in range(self.lr,self.hr):
			for j in range(0,len(self.tracedata[i])-1):
				alt = self.tracedata[i][j][0]
				az = self.tracedata[i][j][1]
				phi = self.tracedata[i][j][2]
				#print "a",alt,az
				T = Transform({"type":"eman","az":az,"alt":alt})
				T.invert()
				a = T*Vec3f(0,0,1)
				
				if (j == 0):
					glPushMatrix()
					d = T.get_rotation("eman")
					glRotate(-d["phi"],0,0,1)
					glRotate(-d["alt"],1,0,0)
					glRotate(-d["az"],0,0,1)
					
					glTranslate(0,0,self.radius)
					glCallList(self.spheredl)
					glPopMatrix()
				
				#print a[0],a[1],a[2]
				alt = self.tracedata[i][j+1][0]
				az = self.tracedata[i][j+1][1]
				#print "b",alt,az
				T = Transform({"type":"eman","az":az,"alt":alt})
				T.invert()
				b = T*Vec3f(0,0,1)
				#print b[0],b[1],b[2]
				##if a == b: continue
				glPushMatrix()
				self.cylinder_to_from(b*self.radius,a*self.radius)
				glPopMatrix()
				#print b,a
		glEndList()
	
		self.updateGL()
		
	def arc_animation_update(self,points):
		self.arc_anim_points = points
		
#		if self.arc_anim_dl:
#			glDeleteLists(self.arc_anim_dl,1)
#			
#		self.arc_anim_dl = glGenLists(1)
#		glNewList(self.arc_anim_dl,GL_COMPILE)
#		
#		for i in xrange(0,len(points)-2):
#			self.cylinder_to_from(points[i],points[i+1])
#		glEndList()
		self.updateGL()
		
	def set_display_all_syms(self,val):
		self.display_all_syms = val

	def set_gl_widget(self,gl_widget):
		'''
		Needed to specialized this function (from EM3DModel) because self.vdtools is absolutely required to be current and accurate
		(this object uses GLU picking and needs access to the viewport matrix) 
		'''
		self.vdtools = EMViewportDepthTools2(gl_widget)
		EM3DModel.set_gl_widget(self, gl_widget)

	def render(self):
		self.init_basic_shapes()
		if self.vdtools == None:
			self.vdtools = EMViewportDepthTools2(self.get_gl_widget())
		
		if self.update_sphere_points_dl:
			self.make_sym_dl_list(self.eulers)
			self.update_sphere_points_dl = False
		
		lighting = glIsEnabled(GL_LIGHTING)
		cull = glIsEnabled(GL_CULL_FACE)
		polygonmode = glGetIntegerv(GL_POLYGON_MODE)
		
		glEnable(GL_LIGHTING)
		glDisable(GL_CULL_FACE)
		
		glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
		
		glPushMatrix()
		self.cam.position(True)
		# the ones are dummy variables atm... they don't do anything
		self.vdtools.update(1,1)
		glPopMatrix()
		
		self.cam.position()
		self.model_matrix = glGetDoublev(GL_MODELVIEW_MATRIX) # this is for e2au.py ... there is a better solution but no time
		
		if ( self.sym_object.is_h_sym() ):
			a = {}
			a["daz"] = 60
			dz = 5
			a["tz"] = dz
			self.sym_object.insert_params(a)
		
		if ( self.sphere_points_dl == 0 or self.force_update):
			self.generate_current_display_list()
			self.force_update = False
			if ( self.sphere_points_dl == 0 ) : 
				print "error, you can't draw an empty list"
				return
			
		for i,t in enumerate(self.sym_object.get_syms()):
			if i > 0 and not self.display_all_syms: break
			d = t.get_rotation("eman")
			glPushMatrix()
			# negated because OpenGL uses the opposite handed coordinate system
			glRotate(-d["phi"],0,0,1)
			glRotate(-d["alt"],1,0,0)
			glRotate(-d["az"],0,0,1)
			
			if ( self.sym_object.is_h_sym() == False and self.asym_u_triangles_dl != 0 and self.display_tri ):
				glCallList(self.asym_u_triangles_dl)
			
			if ( self.great_arc_dl != 0 and self.display_arc ):
				self.load_gl_color(self.arc_color)
				glCallList(self.great_arc_dl)
				
			if self.display_euler and self.sym_object != None:
				glCallList(self.sphere_points_dl)
			glPopMatrix()
			
		if (self.arc_anim_points):
			self.load_gl_color("gold")
			for points in self.arc_anim_points:
				for i in range(0,len(points)-1):
					self.cylinder_to_from(points[i],points[i+1],scale=0.25)

			if len(self.arc_anim_points[-1]) > 0:
				self.cylinder_to_from(Vec3f(0,0,0),self.arc_anim_points[-1][-1],scale=0.25)
		

		if self.trace_dl != 0:
			self.load_gl_color("black")
			glStencilFunc(GL_EQUAL,self.rank,0)
			glStencilOp(GL_KEEP,GL_KEEP,GL_REPLACE)
			#print "rendering trace"
			glPushMatrix()
			glCallList(self.trace_dl)
			glPopMatrix()
		glStencilFunc(GL_EQUAL,self.rank,self.rank)
		glStencilOp(GL_KEEP,GL_KEEP,GL_KEEP)

		self.draw_bc_screen()

		glStencilFunc(GL_ALWAYS,1,1)
		if self.cube:
			glPushMatrix()
			
			self.draw_volume_bounds()
			glPopMatrix()
			
		if ( lighting ): glEnable(GL_LIGHTING)
		if ( cull ): glEnable(GL_CULL_FACE)
		
		if ( polygonmode[0] == GL_LINE ): glPolygonMode(GL_FRONT, GL_LINE)
		if ( polygonmode[1] == GL_LINE ): glPolygonMode(GL_BACK, GL_LINE)
		
	def update_inspector(self,t3d):
		if not self.inspector or self.inspector ==None:
			self.inspector=EMSymInspector(self)
		self.inspector.update_rotations(t3d)
	
	def get_inspector(self):
		if not self.inspector : self.inspector=EMSymInspector(self)
		return self.inspector
		
	def regen_dl(self, dummy=False):
		self.force_update = True
		self.updateGL()
			
	def resize(self):
		pass
	
		#self.vdtools.set_update_P_inv()
		
	def toggle_sym_display(self,bool):
		self.display_euler = bool
		
	def triangletog(self,bool):
		self.display_tri = bool
		
	def arctog(self,bool):
		self.display_arc = bool

	def reducetog(self,bool):
		self.reduce = bool

class EMSymViewerWidget(EMGLWidget, EMGLProjectionViewMatrices):
	allim=WeakKeyDictionary()
	def __init__(self, sym="c1", filename=None):
		EMSymViewerWidget.allim[self]=0
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True)
		fmt.setDepth(1)
		fmt.setSampleBuffers(True)
		
		EMGLWidget.__init__(self)
		self.setFormat(fmt)
		EMGLProjectionViewMatrices.__init__(self)

		self.filename = filename
		self.fov = 50 # field of view angle used by gluPerspective
		self.startz = 1
		self.endz = 5000
		self.cam = Camera()
		self.cam.setCamTrans("default_z",-100)
		
		sym_model = EM3DSymModel(self,eulerfilename=self.filename)
		self.model = sym_model
		sym_model.set_symmetry(sym)
		sym_model.regen_dl()
		
		self.resize(640,640)
		
#	def __del__(self):
#		print "sym viewer widget died"
	def get_near_plane_dims(self):
		height = 2.0 * self.get_start_z()*tan(self.fov/2.0*pi/180.0)
		width = self.aspect * height
		return [width,height]
	def get_render_dims_at_depth(self,depth):
		# This function returns the width and height of the renderable 
		# area at the origin of the data volume
		height = -2*tan(self.fov/2.0*pi/180.0)*(depth)
		width = self.aspect*height
		
		return [width,height]
	def get_start_z(self):
		return self.startz
	def initializeGL(self):
		glEnable(GL_NORMALIZE)
		glEnable(GL_LIGHT0)
		glEnable(GL_DEPTH_TEST)
		#print "Initializing"
#		glLightfv(GL_LIGHT0, GL_AMBIENT, [0.0, 0.0, 0.0, 1.0])
#		glLightfv(GL_LIGHT0, GL_DIFFUSE, [.8,.8,.8, 1.0])
#		glLightfv(GL_LIGHT0, GL_SPECULAR, [.2, .2, .2, 1.0])
#		glLightfv(GL_LIGHT0, GL_POSITION, [0.1,.1,1.,0.])
		glLightfv(GL_LIGHT0, GL_AMBIENT, [0.2, 0.2, 0.2, 1.0])
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [.4,.4,.4, 1.0])
		glLightfv(GL_LIGHT0, GL_SPECULAR, [.1,.1,.1, 1.0])
		glLightfv(GL_LIGHT0, GL_POSITION, [0.5,0.7,11.,0.])
		GL.glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST)
		glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,GL_TRUE)
		GL.glClearColor(0.875,0.875,0.875,1)
		#GL.glClearAccum(0,0,0,0)
	
		glShadeModel(GL_SMOOTH)
		
		#glClearStencil(0)
		#glEnable(GL_STENCIL_TEST)
	def keyPressEvent(self,event):
		if self.model:
			self.model.keyPressEvent(event)
			self.updateGL()
	def load_perspective(self):
		gluPerspective(self.fov,self.aspect,self.startz,self.endz)
	def mouseDoubleClickEvent(self,event):
		if self.model:
			self.model.mouseDoubleClickEvent(event)
			self.updateGL()
	def mouseMoveEvent(self, event):
		if self.model:
			self.model.mouseMoveEvent(event)
			self.updateGL()
	def mousePressEvent(self, event):
		if self.model:
			self.model.mousePressEvent(event)
			self.updateGL()
	def mouseReleaseEvent(self, event):
		if self.model:
			self.model.mouseReleaseEvent(event)
			self.updateGL()
	def paintGL(self):
		#glClear(GL_ACCUM_BUFFER_BIT)
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT )
		
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		self.cam.position()
		
		glPushMatrix()
		self.model.render()
		glPopMatrix()
	def resizeGL(self, width, height):
		if width<=0 or height<=0 : 
			return
		# just use the whole window for rendering
		glViewport(0,0,self.width(),self.height())
		
		# maintain the aspect ratio of the window we have
		self.aspect = float(self.width())/float(self.height())
		
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		# using gluPerspective for simplicity
		gluPerspective(self.fov,self.aspect,self.startz,self.endz)
		
		# switch back to model view mode
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		self.set_projection_view_update()
		self.model.resize()
	def show_inspector(self,force=0):
		self.model.show_inspector(self,force)
	def wheelEvent(self, event):
		if self.model:
			self.model.wheelEvent(event)
		self.updateGL()


class SparseSymChoicesWidgets:
	'''
	An encapsulation of the most basic of symmetry inspector widgits
	Used by both EMSymChoiceDialog and EMSymInspector
	'''
	def __init__(self,widget,target):
		'''
		@param widget some kind of QtWidget, used as the first argument when initializing Validators (really)
		@param target an EM3DSymModel instance
		
		The target is very important, the signals emitted by this object are connect to functions in the 
		EM3DSymViewerWidget
		'''
		self.widget = weakref.ref(widget)
		self.target = weakref.ref(target)
		self.busy = False
	def add_top_buttons(self,vbl):
		'''
		Adds 4 buttons in a grid
		| Display Eulers | Display Triangles |
		| Display Arcs   | All Syms          |
		@param vbl a QtGui.QVBoxLayout - all widgets and layouts are added to it
		'''
		self.busy = True
		self.button_hbl1 = QtGui.QHBoxLayout()
		self.symtogdisplay = QtGui.QPushButton("Display Eulers")
		self.symtogdisplay.setCheckable(1)
		self.symtogdisplay.setChecked(1)
		self.button_hbl1.addWidget(self.symtogdisplay)
		
		self.triangletog = QtGui.QPushButton("Display Triangles")
		self.triangletog.setCheckable(1)
		self.triangletog.setChecked(0)
		self.button_hbl1.addWidget(self.triangletog)
		
		vbl.addLayout(self.button_hbl1)
		
		self.button_hbl2 = QtGui.QHBoxLayout()
		
		self.arctog = QtGui.QPushButton("Display Arcs")
		self.arctog.setCheckable(1)
		self.arctog.setChecked(1)
		self.button_hbl2.addWidget(self.arctog)
		
		self.symtog = QtGui.QPushButton("All syms")
		self.symtog.setCheckable(1)
		self.button_hbl2.addWidget(self.symtog)
		vbl.addLayout(self.button_hbl2)
		
		
		QtCore.QObject.connect(self.symtog, QtCore.SIGNAL("toggled(bool)"), self.set_display_all_syms)
		QtCore.QObject.connect(self.symtogdisplay, QtCore.SIGNAL("clicked(bool)"), self.toggle_sym_display)
		QtCore.QObject.connect(self.triangletog, QtCore.SIGNAL("clicked(bool)"), self.triangle_tog)
		QtCore.QObject.connect(self.arctog, QtCore.SIGNAL("clicked(bool)"), self.arc_tog)
		self.busy = False

	def set_display_all_syms(self,val):
		'''
		This function is a slot for the signal that is emitted when the self.symtog button is checked
		'''
		if self.busy: return
		self.target().set_display_all_syms(val)
		self.target().updateGL()
		
	def toggle_sym_display(self,val):
		'''
		This function is a slot for the signal that is emitted when the self.symtogdisplay button is checked
		'''
		if self.busy: return
		self.target().toggle_sym_display(val)
		self.target().updateGL()
		
	def triangle_tog(self,val):
		'''
		This function is a slot for the signal that is emitted when the self.triangletog button is checked
		'''
		if self.busy: return
		self.target().triangletog(val)
		self.target().updateGL()
		
	def arc_tog(self,val):
		'''
		This function is a slot for the signal that is emitted when the self.arctog button is checked
		'''
		if self.busy: return
		self.target().arctog(val)
		self.target().updateGL()
	
	def add_symmetry_options(self,vbl,enable_orient_gen=True):
		'''
		Add common symmetry options to a QtGui.QVBoxLayout
		@param vbl a QtGui.QVBoxLayout - all widgets and layouts are added to it
		@param enable_orient_gen a boolean indicating whether or not the user should be permitted to change the distribution of eulers - this is False when using E2eulerxplor to examine refinement results
		
		Makes QtCore.QObject connections to functions of self.target() (see bottom of this function)
		'''
		self.busy = True
		self.maintab = QtGui.QWidget()
		maintab = self.maintab
		maintab.vbl = QtGui.QVBoxLayout(self.maintab)
		maintab.vbl.setMargin(0)
		maintab.vbl.setSpacing(6)
		maintab.vbl.setObjectName("Main")
		
		self.hbl_sym = QtGui.QHBoxLayout()
		self.hbl_sym.setMargin(0)
		self.hbl_sym.setSpacing(6)
		self.hbl_sym.setObjectName("Sym")
		maintab.vbl.addLayout(self.hbl_sym)
		
		self.sym_combo = QtGui.QComboBox(maintab)
		self.symmetries = []
		self.symmetries.append(' Icosahedral ')
		self.symmetries.append(' Octahedral ')
		self.symmetries.append(' Tetrahedral ')
		self.symmetries.append(' D ')
		self.symmetries.append(' C ')

		self.sym_map = {}
		self.sym_map[" Icosahedral "] = "icos"
		self.sym_map[" Octahedral "] = "oct"
		self.sym_map[" Tetrahedral "] = "tet"
		self.sym_map[" D "] = "d"
		self.sym_map[" C "] = "c"

		idx_default = 0
		for idx,i in enumerate(self.symmetries): 
			self.sym_combo.addItem(i)
				
		self.sym_combo.setCurrentIndex(idx_default)
		self.hbl_sym.addWidget(self.sym_combo)
		
		self.sym_label = QtGui.QLabel()
		self.sym_label.setText('C/D sym')
		self.hbl_sym.addWidget(self.sym_label)
		
		self.pos_int_validator = QtGui.QIntValidator(self.widget())
		self.pos_int_validator.setBottom(1)
		self.sym_text = QtGui.QLineEdit()
		self.sym_text.setValidator(self.pos_int_validator)
		self.sym_text.setText("7")
		self.sym_text.setFixedWidth(50)
		self.hbl_sym.addWidget(self.sym_text)
		self.sym_text.setEnabled(False)
		
		self.set_sym(self.target().get_sym())
		
		if enable_orient_gen:
			self.angle_label = QtGui.QComboBox()
			self.angle_label.addItem('Angle Based')
			self.angle_label.addItem('Number Based')
			self.hbl_sym.addWidget(self.angle_label)
			
			self.pos_double_validator = QtGui.QDoubleValidator(self.widget())
			self.pos_double_validator.setBottom(0.05)
			self.prop_text = QtGui.QLineEdit()
			self.prop_text.setValidator(self.pos_double_validator)
			self.prop_text.setText(str(self.target().get_prop()))
			self.prop_text.setFixedWidth(50)
			self.hbl_sym.addWidget(self.prop_text)
			
			self.hbl_sym2 = QtGui.QHBoxLayout()
			self.hbl_sym2.setMargin(0)
			self.hbl_sym2.setSpacing(6)
			self.hbl_sym2.setObjectName("Sym2")
			maintab.vbl.addLayout(self.hbl_sym2)
			
			self.og_label = QtGui.QLabel()
			self.og_label.setText('Strategy')
			self.hbl_sym2.addWidget(self.og_label)
			
			self.strategy_label = QtGui.QComboBox()
			l = dump_orientgens_list()
				
			n = len(l)
			for i in l:
				self.strategy_label.addItem(str(i))
				
			self.strategy_label.setCurrentIndex(n-1)
			self.hbl_sym2.addWidget(self.strategy_label)
			
			self.mirror_checkbox = QtGui.QCheckBox("Mirror")
			self.hbl_sym2.addWidget(self.mirror_checkbox)
			self.mirror_checkbox.setChecked(self.target().mirror_enabled())
		else:
			self.mirror_checkbox = QtGui.QCheckBox("Mirror")
			self.hbl_sym.addWidget(self.mirror_checkbox)
			self.mirror_checkbox.setChecked(self.target().mirror_enabled())
			
		QtCore.QObject.connect(self.sym_combo, QtCore.SIGNAL("currentIndexChanged(QString)"), self.sym_changed)
		QtCore.QObject.connect(self.sym_text, QtCore.SIGNAL("editingFinished()"), self.sym_number_changed)

		if enable_orient_gen:
			QtCore.QObject.connect(self.prop_text, QtCore.SIGNAL("editingFinished()"), self.prop_changed)
			QtCore.QObject.connect(self.angle_label, QtCore.SIGNAL("currentIndexChanged(QString)"), self.angle_label_changed)
			QtCore.QObject.connect(self.strategy_label, QtCore.SIGNAL("currentIndexChanged(QString)"), self.strategy_changed)
		QtCore.QObject.connect(self.mirror_checkbox, QtCore.SIGNAL("stateChanged(int)"), self.set_mirror)
	
		vbl.addWidget(maintab)
		self.busy = False
		
	def set_mirror(self,val):
		'''
		This function is a slot for the signal that is emitted when the self.mirror_checkbox is checked
	 	@param val the value that was delivered by Qt, indicated whether or not the checkbox is checked 
		'''
		if self.busy: return
		self.target().set_mirror(val)
		self.target().regen_dl()
	
	def prop_changed(self):
		'''
		This function is a slot for the signal that is emitted when the self.prop_text text box is altered 
		'''
		if self.busy: return
		self.target().set_prop(float(self.prop_text.text()))
		self.target().regen_dl()
		
	def angle_label_changed(self,string):
		'''
		This function is a slot for the signal that is emitted when the self.angle_label combo box has its current index altered
		@param the current text of self.angle_label
		'''
		if self.busy: return
		if string == "Angle Based": s = "delta"
		elif string == "Number Based": s = "n"
		else: raise RuntimeError
		
		self.target().set_angle_label(s)
		self.target().regen_dl()
		
	def strategy_changed(self,string):
		'''
		This function is a slot for the signal that is emitted when the self.strategy_label combo box has its current index altered
		@param the current text of self.strategy_label
		'''
		if self.busy: return
		self.target().set_strategy(str(string))
		self.target().regen_dl()
		
	def sym_number_changed(self):
		'''
		This function is a slot for the signal that is emitted when the self.sym_tex text box is altered 
		'''
		if self.busy: return
		self.target().set_symmetry(self.get_sym())
		self.target().regen_dl()
		
	def sym_changed(self, sym):
		'''
		This function is a slot for the signal that is emitted when the self.sym_label combo box has its current index altered
		@param the current text of self.sym_label
		'''
		if self.busy: return
		if sym == ' D ' or sym == ' C ' or sym == ' H ':
			self.sym_text.setEnabled(True)
		else:
			self.sym_text.setEnabled(False)
		
		self.target().set_symmetry(self.get_sym())
		self.target().regen_dl()

	def set_sym(self,sym):
		'''
		A way of setting the widgets created by this object
		@param sym a symmetry string such as "c3", "icos" etc
		'''
		s = sym.lower()
		if s == "icos":
			self.sym_combo.setCurrentIndex(0)
		elif s == "oct":
			self.sym_combo.setCurrentIndex(1)
		elif s == "tet":
			self.sym_combo.setCurrentIndex(2)
		else:
			try:
				nsym = s[1:]
			except:
				print "can't interpret",sym
				return
			self.sym_text.setEnabled(True)
			if s[0] == "d":
				self.sym_combo.setCurrentIndex(3)
			elif s[0] == "c":
				self.sym_combo.setCurrentIndex(4)
			else:
				print "can't interpret",sym
				return
			
			self.sym_text.setText(s[1:])
			
	def get_sym(self):
		'''
		A convenient way to get the symmetry currently defined by the state of the widgets
		@return a symmetry string such as "c3", "icos", "h3" etc
		'''
		sym = self.sym_map[str(self.sym_combo.currentText())]
		if sym in ['c','d']:
			sym = sym+self.sym_text.displayText()
		return str(sym)
	
	
	def get_params_dict(self):
		'''
		This function is called when the user hits ok in the EMSymChoiceDialog
		@return a dictionary that embodies the main parameters associated with this object
		'''
		d = {}
		# warning - do not change the name of these keys unless you know what you're doing. The workflow refine forms assumes these keys exist, as they are
		# see EMOrientationDistDialog in emform.py
		d["sym"] = self.get_sym()
		print self.get_sym()
		d["orientgen"] = str(self.strategy_label.currentText())
		d["approach"] = str(self.angle_label.currentText())
		d["value"] = float(self.prop_text.text())
		d["inc_mirror"] =  self.mirror_checkbox.isChecked()
		
		return d

class EMSymChoiceDialog(QtGui.QDialog):
	'''
	This is a dialog one can use to get the parameters you can use for 
	generating orientations evenly covering the asymmetric unit (etc)
	Just call exec_ and get the return value - if it's a dictionary
	then the user hit ok, if it's None then they hit cancel
	'''
	def __init__(self,sym="d7"):
		'''
		@param sym some kind of symmetry, such as "d7", "icos" etc
		'''
		QtGui.QDialog.__init__(self)		
		self.setWindowTitle("Choose Distribution Parameters")
		self.setWindowIcon(QtGui.QIcon(get_image_directory() + "eulerxplor.png"))

		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		
		self.sym_widget = EMSymViewerWidget(sym)		
		self.sym_widget.model.enable_inspector(False)
		
		self.sparse_syms_widgets = SparseSymChoicesWidgets(self,self.sym_widget.model)
		self.sparse_syms_widgets.add_top_buttons(self.vbl)
		self.sparse_syms_widgets.add_symmetry_options(self.vbl)
		self.sparse_syms_widgets.set_sym(sym)

		self.vbl.addWidget(self.sym_widget,10)
		
		self.button_hbl = QtGui.QHBoxLayout()
		self.ok = QtGui.QPushButton("Ok")
		self.ok.setDefault(True)
		self.cancel = QtGui.QPushButton("Cancel")
		self.button_hbl.addWidget(self.cancel )
		self.button_hbl.addWidget(self.ok )
		self.vbl.addLayout(self.button_hbl)
		
		self.resize(300,400)
	
		self.dialog_result = None
	
		QtCore.QObject.connect(self.ok, QtCore.SIGNAL("clicked(bool)"), self.on_ok)
		QtCore.QObject.connect(self.cancel, QtCore.SIGNAL("clicked(bool)"), self.on_cancel)
		

		
	def on_ok(self,unused=None):
		'''
		Slot which is connected to the self.ok button
		'''
		self.dialog_result = self.sparse_syms_widgets.get_params_dict()
		self.accept()
		
	def on_cancel(self,unused=None):
		'''
		Slot which is connected to the self.cancel button
		'''
		self.accept()
	
	def exec_(self):
		'''
		Customized exec_ function
		@return None if the user hit cancel or a dictionary containing important parameters if the user hit ok
		'''
		QtGui.QDialog.exec_(self)
		return self.dialog_result
		
	
class EMSymInspector(QtGui.QWidget):
	def __init__(self,target,enable_trace=True,enable_og=True) :
		self.busy = True
		QtGui.QWidget.__init__(self,None)
		self.setWindowIcon(QtGui.QIcon(get_image_directory() + "eulerxplor.png"))
		self.target=weakref.ref(target)
		
		self.score_options_hbl  = None # will eventually be a combo 
		self.rotation_sliders = EMTransformPanel(self.target(),self)
		self.enable_trace = enable_trace
		self.enable_og = enable_og
		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		
		self.sparse_syms_widgets = SparseSymChoicesWidgets(self,self.target())
		self.sparse_syms_widgets.add_top_buttons(self.vbl)
		
		self.add_symmetry_options()
		self.n3_showing = False
		
		self.tabwidget = QtGui.QTabWidget()
		self.tabwidget.addTab(self.get_display_tab(), "Display")
		self.tabwidget.addTab(self.get_transform_tab(), "Transform")
		self.vbl.addWidget(self.tabwidget)
		
		self.file = None
		self.lr = -1
		self.hr = -1
		self.busy = False

	def update_rotations(self,t3d):
		self.rotation_sliders.update_rotations(t3d)
	
	def set_scale(self,val):
		self.rotation_sliders.set_scale(val)
	
	def set_xy_trans(self, x, y):
		self.rotation_sliders.set_xy_trans(x,y)
	
	def set_xyz_trans(self,x,y,z):
		self.rotation_sliders.set_xyz_trans(x,y,z)	
		
	def get_transform_layout(self):
		return self.maintab.vbl
	
	def set_score_options(self,options,default=None):
		self.busy = True
		if options == None or len(options) == 0:
			if self.score_options_hbl != None:
				self.display_tab.vbl.removeItem(self.score_options_hbl)
				self.score_options_hbl.deleteLater()
				QtCore.QObject.disconnect(self.score_options,QtCore.SIGNAL("currentIndexChanged(int)"),self.score_option_changed)
			return

		if self.score_options_hbl == None:
			self.score_options_hbl = QtGui.QHBoxLayout()
			self.score_options = self.__get_combo(options,default)
			self.score_options_hbl.addWidget(QtGui.QLabel("Cylinder Score:",self))
			self.score_options_hbl.addWidget(self.score_options)
			self.cylinder_log = QtGui.QCheckBox("log scale")
			self.cylinder_log.setChecked(self.target().log_scale)
			self.score_options_hbl.addWidget(self.cylinder_log)
			self.display_tab.vbl.addLayout(self.score_options_hbl)
			
			QtCore.QObject.connect(self.score_options,QtCore.SIGNAL("currentIndexChanged(int)"),self.score_option_changed)
			QtCore.QObject.connect(self.cylinder_log,QtCore.SIGNAL("stateChanged(int)"),self.cylinder_log_clicked)
		else:
			self.score_options.clear()
			idx = 0
			for i,k in enumerate(options): 
				self.score_options.addItem(k)
				if k == default: idx = i
			self.score_options.setCurrentIndex(idx)
		self.busy = False
		
		
	def cylinder_log_clicked(self,val):
		self.target().set_log_scale(val)
		self.target().regen_dl()
		
	def score_option_changed(self):
		if self.busy: return
		s = str(self.score_options.currentText())
		self.target().set_column_score_key(s)
		self.target().regen_dl()
			
	def set_sym(self,sym):
		self.sparse_syms_widgets.set_sym(sym)

	def get_sym(self):
		return self.sparse_syms_widgets.get_sym()

	def trace_update(self):
		lr = int(self.lowrange.displayText())
		hr = int(self.highrange.displayText())
		file = str(self.tracefile.displayText())
		
		if ( file != self.file or lr != self.lr or hr != self.hr ):
			self.file = file
			self.hr = hr
			self.lr = lr
			
			self.target().trace_update(file,lr,hr)

	def toggle_trace(self,bool):
		self.tracefile.setEnabled(bool)
		self.lowrange.setEnabled(bool)
		self.highrange.setEnabled(bool)
		
		
	def get_transform_tab(self):
		self.transform_tab = QtGui.QWidget()
		self.transform_tab.vbl = QtGui.QVBoxLayout(self.transform_tab)
		self.rotation_sliders.addWidgets(self.transform_tab.vbl)
		
		return self.transform_tab
	
	def arc_color_changed(self):
		if self.busy: return
		s = str(self.arc_color.currentText())
		self.target().set_arc_color(s)
		self.target().regen_dl()
		
	def small_column_color_changed(self):
		if self.busy: return
		s = str(self.small_column_color.currentText())
		self.target().set_small_column_color(s)
		self.target().regen_dl()
		
	def tall_column_color_changed(self):
		if self.busy: return
		s = str(self.tall_column_color.currentText())
		self.target().set_tall_column_color(s)
		self.target().regen_dl()
	
	def get_display_tab(self):
		
		self.display_tab = QtGui.QWidget()
		self.display_tab.vbl = QtGui.QVBoxLayout(self.display_tab)
				
#		self.glcontrast = ValSlider(self.display_tab,(1.0,5.0),"GLShd:")
#		self.glcontrast.setObjectName("GLShade")
#		self.glcontrast.setValue(1.0)
#		self.display_tab.vbl.addWidget(self.glcontrast)
#		
#		self.glbrightness = ValSlider(self.display_tab,(-1.0,0.0),"GLBst:")
#		self.glbrightness.setObjectName("GLBoost")
#		self.glbrightness.setValue(0.1)
#		self.glbrightness.setValue(0.0)
#		self.display_tab.vbl.addWidget(self.glbrightness)
		
		keys = self.target().colors.keys()
		keys.sort()
		self.arc_color = self.__get_combo(keys,self.target().arc_color)
		hbl1 = QtGui.QHBoxLayout()
		hbl1.addWidget(QtGui.QLabel("Arc Color:",self))
		hbl1.addWidget(self.arc_color)
		self.display_tab.vbl.addLayout(hbl1)
		
		self.tall_column_color = self.__get_combo(keys,self.target().tall_column_color)
		hbl2 = QtGui.QHBoxLayout()
		hbl2.addWidget(QtGui.QLabel("Higher Cylinder Color:",self))
		hbl2.addWidget(self.tall_column_color)
		self.display_tab.vbl.addLayout(hbl2)
		
		
		self.small_column_color = self.__get_combo(keys,self.target().small_column_color)
		hbl3 = QtGui.QHBoxLayout()
		hbl3.addWidget(QtGui.QLabel("Lower Cylinder Color:",self))
		hbl3.addWidget(self.small_column_color)
		self.display_tab.vbl.addLayout(hbl3)
		
		self.width_scale = ValSlider(self,(0.05,10.0),"Cylinder Width:")
		self.width_scale.setValue(self.target().width_scale)
		self.display_tab.vbl.addWidget(self.width_scale)
		
		self.height_scale = ValSlider(self,(0.05,10.0),"Cylinder Height:")
		self.height_scale.setValue(self.target().height_scale)
		self.display_tab.vbl.addWidget(self.height_scale)
		
		self.arc_width_scale = ValSlider(self,(0.05,10.0),"Arc Width:")
		self.arc_width_scale.setValue(self.target().arc_width_scale)
		self.display_tab.vbl.addWidget(self.arc_width_scale)
		
		hbl_l = QtGui.QHBoxLayout()
		arc_div_label = QtGui.QLabel("Arc Segments:")
		arc_div_label.setAlignment(Qt.AlignLeft|Qt.AlignVCenter)
		hbl_l.addWidget(arc_div_label)
		self.arc_divisions = QtGui.QSpinBox(self)
		self.arc_divisions.setRange(1,1000)
		self.arc_divisions.setValue(int(self.target().arc_segments))
		hbl_l.addWidget(self.arc_divisions)
		self.display_tab.vbl.addLayout(hbl_l)
		
		QtCore.QObject.connect(self.width_scale, QtCore.SIGNAL("valueChanged"), self.target().set_width_scale)
		QtCore.QObject.connect(self.height_scale, QtCore.SIGNAL("valueChanged"), self.target().set_height_scale)
		QtCore.QObject.connect(self.arc_width_scale, QtCore.SIGNAL("valueChanged"), self.target().set_arc_width_scale)
#		QtCore.QObject.connect(self.glcontrast, QtCore.SIGNAL("valueChanged"), self.target().set_GL_contrast)
#		QtCore.QObject.connect(self.glbrightness, QtCore.SIGNAL("valueChanged"), self.target().set_GL_brightness)
		QtCore.QObject.connect(self.arc_color,QtCore.SIGNAL("currentIndexChanged(int)"),self.arc_color_changed)
		QtCore.QObject.connect(self.small_column_color,QtCore.SIGNAL("currentIndexChanged(int)"),self.small_column_color_changed)
		QtCore.QObject.connect(self.tall_column_color,QtCore.SIGNAL("currentIndexChanged(int)"),self.tall_column_color_changed)
		QtCore.QObject.connect(self.arc_divisions, QtCore.SIGNAL("valueChanged(int)"), self.target().set_arc_segments)
		
		return self.display_tab
	
	def __get_combo(self,keys,default):
		combo = QtGui.QComboBox()
		idx = 0
		for i,k in enumerate(keys): 
			combo.addItem(k)
			if k == default: idx = i
		combo.setCurrentIndex(idx)
		return combo
	
	def add_symmetry_options(self):
		
		self.sparse_syms_widgets.add_symmetry_options(self.vbl,self.enable_og)

		if self.enable_trace:		
			self.hbl_pt = QtGui.QHBoxLayout()
			self.hbl_pt.setMargin(0)
			self.hbl_pt.setSpacing(6)
			self.hbl_pt.setObjectName("Ptl Trace")
			
			
			self.tracetog = QtGui.QPushButton("Trace")
			self.tracetog.setCheckable(1)
			self.tracetog.setChecked(0)
			self.hbl_pt.addWidget(self.tracetog)
			
			self.tracefile = QtGui.QLineEdit()
			self.tracefile.setText("filename.txt")
			self.tracefile.setFixedWidth(100)
			self.hbl_pt.addWidget(self.tracefile)
			self.tracefile.setEnabled(False)
			
			self.pt_label = QtGui.QLabel()
			self.pt_label.setText('Range')
			self.hbl_pt.addWidget(self.pt_label)
			
			self.pos_int_validator2 = QtGui.QIntValidator(self)
			self.pos_int_validator2.setBottom(0)
			self.lowrange = QtGui.QLineEdit()
			self.lowrange.setValidator(self.pos_int_validator2)
			self.lowrange.setText("1")
			self.lowrange.setFixedWidth(50)
			self.hbl_pt.addWidget(self.lowrange)
			self.lowrange.setEnabled(False)
			
			self.pt_label_to = QtGui.QLabel()
			self.pt_label_to.setText('to')
			self.hbl_pt.addWidget(self.pt_label_to)
			
			self.highrange = QtGui.QLineEdit()
			self.highrange.setValidator(self.pos_int_validator2)
			self.highrange.setText("1")
			self.highrange.setFixedWidth(50)
			self.hbl_pt.addWidget(self.highrange)
			self.highrange.setEnabled(False)
			
			self.reducetog = QtGui.QPushButton("Reduce")
			self.reducetog.setCheckable(1)
			self.reducetog.setChecked(0)
			self.hbl_pt.addWidget(self.reducetog)
			
			self.vbl.addLayout(self.hbl_pt)

		if self.enable_trace:
			QtCore.QObject.connect(self.tracetog, QtCore.SIGNAL("clicked(bool)"), self.toggle_trace)
			QtCore.QObject.connect(self.reducetog, QtCore.SIGNAL("clicked(bool)"), self.target().reducetog)
			QtCore.QObject.connect(self.lowrange, QtCore.SIGNAL("editingFinished()"), self.trace_update)
			QtCore.QObject.connect(self.highrange, QtCore.SIGNAL("editingFinished()"), self.trace_update)
			QtCore.QObject.connect(self.tracefile, QtCore.SIGNAL("editingFinished()"), self.trace_update)

		
	def slider_rotate(self):
		if self.busy: return
		self.target().load_rotation(self.get_current_rotation())

	
if __name__ == '__main__':
	em_app = EMApp()
	
	#First demonstration
	dialog = EMSymChoiceDialog()
	choices_dict = dialog.exec_()
	print choices_dict
	
	#Second demonstration
	window = EMSymViewerWidget()
	
	em_app.show()
	em_app.execute()