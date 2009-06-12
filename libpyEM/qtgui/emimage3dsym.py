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



from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from OpenGL import GL,GLU,GLUT
from OpenGL.GL import *
from OpenGL.GLU import *
from valslider import ValSlider
from math import *
from EMAN2 import *
import sys
import numpy
from weakref import WeakKeyDictionary
from time import time
from PyQt4.QtCore import QTimer

from time import *

from emglobjects import EMImage3DGUIModule, Camera, Camera2, EMViewportDepthTools2,EMGLProjectionViewMatrices
from emimageutil import ImgHistogram, EMEventRerouter, EMTransformPanel
from emapplication import EMStandAloneApplication, EMQtWidgetModule, EMGUIModule
import weakref

MAG_INCREMENT_FACTOR = 1.1

class EM3DSymViewerModule(EMImage3DGUIModule):

	def get_qt_widget(self):
		if self.qt_context_parent == None:	
			self.under_qt_control = True
			from emimageutil import EMParentWin
			self.gl_context_parent = EMSymViewerWidget(self)
			self.qt_context_parent = EMParentWin(self.gl_context_parent)
			self.qt_context_parent.setWindowIcon(QtGui.QIcon(get_image_directory() + "eulerxplor.png"))
			self.qt_context_parent.setWindowTitle("Eulerxplor")
			self.gl_widget = self.gl_context_parent
		return self.qt_context_parent


	def __init__(self,application=None,inspector_go=True,ensure_gl_context=True,application_control=True):
		EMImage3DGUIModule.__init__(self,application,ensure_gl_context,application_control)
		
		
		self.eulers = [] # will eventually store Transform objects
		self.points = [] # will eventually store the points on the asymmetric unit
		self.point_colors = [] # will eventually store colors for the different points
		self.init()
		self.initialized = True
		
		if inspector_go: self.get_inspector()
		
#	def __del__(self):
#		print "sym died" 
		
	def get_inspector(self):
		if not self.inspector : self.inspector=EMSymInspector(self)
		return self.inspector
	
	def keyPressEvent(self,event):
		
		if event.key() == Qt.Key_F1:
			self.display_web_help("http://blake.bcm.edu/emanwiki/EMAN2/Programs/e2eulerxplor")
		else:
			EMImage3DGUIModule.keyPressEvent(self,event)
	
	def get_type(self):
		return "Symmetry Viewer"
	
	def update_data(self,data):
		pass

	def init(self):
		self.data=None

		self.mmode=0
		self.cam = Camera2(self)
		self.vdtools = None # EMViewportDepthTools2(self.get_gl_context_parent())
		
		self.cube = False
		self.inspector=None
		
		self.sym_dl = 0
		self.spheredl = 0
		self.highresspheredl = 0
		self.diskdl = 0

		self.glcontrast = 1.0
		self.glbrightness = 0.0
		
		self.rank = 1
		
		self.force_update = False
		self.force_force_update = False
		
		self.sym = "icos"
		self.prop = 5
		self.perturb = False
		self.nomirror = True
		self.angle_label = 'delta'
		self.strategy = 'eman'
		
		self.display_euler = True
		self.display_tri = False
		self.display_arc = True
		self.sym_object = None
		self.update_sym_dl = True
		
		self.radius = 50
		
		self.arc_t = 16
		self.arc_dl = 0
		
		self.tri_dl = 0
		
		self.cylinderdl = 0
		
		
		self.reduce = False
		
		self.eulers_specified = False # from e2au.py = if you want to specify a specific set of eulers this is required
		self.specified_eulers = None # from e2au.py = if you want to specify a specific set of eulers this is required
		self.colors_specified = False # as above
		self.specified_colors = None # as above
		
		self.file = None
		self.lr = -1
		self.hr = -1
		self.tracedata = []
		self.trace_dl = 0
		self.gq=gluNewQuadric()
		gluQuadricDrawStyle(self.gq,GLU_FILL)
		gluQuadricNormals(self.gq,GLU_SMOOTH)
		gluQuadricOrientation(self.gq,GLU_OUTSIDE)
		gluQuadricTexture(self.gq,GL_FALSE)
		
		self.set_sym(self.sym)
		self.regen_dl()
		
		
		#self.glbasicobjects = EMBasicOpenGLObjects()
		#print self.glbasicobjects.getSphereDL()
	
	def eye_coords_dif(self,x1,y1,x2,y2,mdepth=True):
		return self.vdtools.eye_coords_dif(x1,y1,x2,y2,mdepth)

	def set_radius(self,radius):
		if ( radius > 0 ):
			self.radius = radius
			self.force_force_update = True
			self.updateGL()
		else:
			print "Error, tried to set a zero or negative radius (",radius,")"
			exit(1)
	
	def trace_great_triangles(self,inc_mirror):
		triangles = self.sym_object.get_asym_unit_triangles(inc_mirror)
		if ( self.tri_dl != 0 ): glDeleteLists(self.tri_dl, 1)
		
		self.tri_dl=glGenLists(1)
		
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
	
		glNewList(self.tri_dl,GL_COMPILE)
		
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
				self.gold()
			elif (j == 1):
				self.red()
			elif ( j == 2):
				self.green()
		elif ( l == 2):
			if i == 0:
				if (j == 0):
					self.gold()
				elif (j == 1):
					self.red()
				elif ( j == 2):
					self.green()
			elif i == 1:
				if (j == 0):
					self.gold()
				elif (j == 1):
					self.green()
				elif ( j == 2):
					self.red()
		else:
			self.gold()
				
	
	def gl_color(self,color):
		glColor(*color)
		glMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,color)
		glMaterial(GL_FRONT,GL_SPECULAR,color)
		glMaterial(GL_FRONT,GL_EMISSION,(0,0,0,1.0))
		glMaterial(GL_FRONT,GL_SHININESS,64.0)
	
	def gold(self):
		glMaterial(GL_FRONT,GL_AMBIENT,(0.24725, 0.2245, 0.0645,1.0))
		glMaterial(GL_FRONT,GL_DIFFUSE,(0.34615, 0.3143, 0.0903,1.0))
		glMaterial(GL_FRONT,GL_SPECULAR,(1.000, 0.9079885, 0.26086934,1.0))
		glMaterial(GL_FRONT,GL_EMISSION,(0,0,0,1.0))
		glMaterial(GL_FRONT,GL_SHININESS,4.0)
		
	def black(self):
		glMaterial(GL_FRONT,GL_AMBIENT,(0.05375,  0.05,     0.06625 ,1.0))
		glMaterial(GL_FRONT,GL_DIFFUSE,(0.18275,  0.17,     0.22525,1.0))
		glMaterial(GL_FRONT,GL_SPECULAR,(0.66, 0.65, 0.69,1.0))
		glMaterial(GL_FRONT,GL_EMISSION,(0,0,0,1.0))
		glMaterial(GL_FRONT,GL_SHININESS, 128.0)
		
	def green(self):
		glColor(.2,.9,.2)
		# this is a nice light blue color (when lighting is on)
		# and is the default color of the frame
		glMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,(.1,.3,.1,1.0))
		glMaterial(GL_FRONT,GL_SPECULAR,(.8,.8,.8,1.0))
		glMaterial(GL_FRONT,GL_EMISSION,(0,0,0,1.0))
		glMaterial(GL_FRONT,GL_SHININESS,32.0)
		
	def red(self):
		glMaterial(GL_FRONT,GL_AMBIENT,(0.1745, 0.01175, 0.01175,1.0))
		glMaterial(GL_FRONT,GL_DIFFUSE,(0.61424, 0.04136, 0.04136,1.0))
		glMaterial(GL_FRONT,GL_SPECULAR,(0.927811, 0.826959, 0.826959,1.0))
		glMaterial(GL_FRONT,GL_EMISSION,(0,0,0,1.0))
		glMaterial(GL_FRONT,GL_SHININESS,32.0)
		
	def trace_great_arcs(self, points):
		if ( self.cylinderdl == 0 ):
			self.cylinderdl=glGenLists(1)
				
			glNewList(self.cylinderdl,GL_COMPILE)
			glPushMatrix()
			gluCylinder(self.gq,1.0,1.0,1.0,12,2)
			glPopMatrix()
				
			glEndList()
		
		if ( self.arc_dl != 0 ): glDeleteLists(self.arc_dl, 1)
		n = len(points)
		if ( n <= 1 ):
			self.arc_dl = 0
			return
		
		self.arc_dl=glGenLists(1)
		
		glNewList(self.arc_dl,GL_COMPILE)
		
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
			for t in range(1, self.arc_t+1):
				timeangle = float(t)/float(self.arc_t)*angle
				p1Copy = self.radius*Vec3f(p1[0],p1[1],p1[2])
				p2Copy = self.radius*Vec3f(p2[0],p2[1],p2[2])
				next = (sin(angle-timeangle)*p1Copy + sin(timeangle)*p2Copy)/sinangle
				
				self.cylinder_to_from(next,prev)
				prev = Vec3f(next[0],next[1],next[2])
				
	
	def cylinder_to_from(self,next,prev):
		dx = next[0] - prev[0]
		dy = next[1] - prev[1]
		dz = next[2] - prev[2]
		
		length = sqrt(dx**2 + dy**2 + dz**2)
		
		if length == 0: return
		
		alt = acos(dz/length)*180.0/pi
		phi = atan2(dy,dx)*180.0/pi
		
		glPushMatrix()
		glTranslatef(prev[0],prev[1],prev[2] )
		#print "positioned at", prev[0],prev[1],prev[2]
		glRotatef(90+phi,0,0,1)
		glRotatef(alt,1,0,0)
		
		glScalef(0.5,0.5,length)
		glCallList(self.cylinderdl)
		glPopMatrix()

	def specify_colors(self,colors):
		self.colors_specified = True
		self.specified_colors = colors

	def specify_eulers(self,eulers):
		self.eulers_specified = True
		self.specified_eulers = eulers
	
	def generate_current_display_list(self,force=False):
		if self.inspector != None:
			sym = self.inspector.get_sym()
			prop = self.inspector.get_prop()
			mirror = not self.inspector.get_mirror()
			perturb = self.inspector.get_perturb()
			angle_label = self.inspector.get_angle_label()
			strategy = str(self.inspector.get_orient_label())
		
#		if not force and self.sym == sym and self.prop == prop and self.nomirror == mirror and self.perturb == perturb and self.angle_label == angle_label and self.strategy == strategy and self.force_force_update == False: return
#		else:
			self.sym = sym
			self.prop = prop
			self.nomirror = mirror
			self.perturb = perturb
			self.angle_label = angle_label
			self.strategy = strategy
		
		if self.diskdl == 0:
			self.diskdl=glGenLists(1)
				
			glNewList(self.diskdl,GL_COMPILE)
			gluDisk(self.gq,.1,.8,8,2)
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

		if ( self.sym_dl != 0 ): glDeleteLists( self.sym_dl, 1)
		
		self.sym_dl = glGenLists(1)
	
		if (self.sym_dl == 0):
			self.sym = None
			self.prop = None
			return #OpenGL is not initialized yet
		self.sym_object = parsesym(str(self.sym))
		
		
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
			
		if ('perturb' in parms):
			if ( self.perturb == True ) : val = 1
			else: val = 0
			og += ":perturb="+ str(val)
		
		[og_name,og_args] = parsemodopt(str(og))
		
		filedebug = False
		if ( not filedebug and not self.eulers_specified):
			eulers = self.sym_object.gen_orientations(og_name, og_args)
		elif self.eulers_specified:
			eulers = self.specified_eulers
		else:
			f = file("angles.txt")
			lines=f.readlines()
			angles=[]
			eulers = []
			for line in lines:
				angles = str.split(line)
				alt = angles[0]
				az = angles[1]
				eulers.append(Transform({"type":"eman","az":az}))
		
		self.eulers = eulers
		if not self.colors_specified: self.point_colors = []
		else: self.point_colors = self.specified_colors
		self.points = []
		for i in eulers:
			p = i.transpose()*Vec3f(0,0,self.radius)
			self.points.append(p)
			if not self.colors_specified: self.point_colors.append((0.34615, 0.3143, 0.0903,1))
		self.make_sym_dl_list(self.points,self.point_colors,eulers)
	
	def set_sym(self,sym):
		'''
		Warning, this probably doesn't work in general don't use it
		'''
		self.sym = sym
		self.sym_object = parsesym(str(sym))
		self.get_inspector().set_sym(sym)
		self.force_update = True
	
	def make_sym_dl_list(self,points,point_colors,eulers):
		
		if self.sym_dl > 0:
			glDeleteLists(self.sym_dl,1)
		self.sym_dl = glGenLists(1)
		glNewList(self.sym_dl,GL_COMPILE)
		for i,p in enumerate(eulers):
			self.gl_color(point_colors[i])
			glPushMatrix()
			d = p.get_rotation("eman")
			glRotate(d["az"],0,0,1)
			glRotate(d["alt"],1,0,0)
			glRotate(d["phi"],0,0,1)
#			glTranslate(p[0],p[1],p[2])
	   	   	glTranslate(0,0,self.radius)
			glCallList(self.diskdl)
			glPopMatrix()
			
		glEndList()
		
	def set_point_colors(self,new_color_map):
		for k,v in new_color_map.items():
			try:
				self.point_colors[k] = v
			except: pass
		self.update_sym_dl = True
		#self.make_sym_dl_list(self.points,self.point_colors)	
	
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
		print "going through range",self.lr,self.hr
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
		
	def render(self):
		if self.vdtools == None:
			self.vdtools = EMViewportDepthTools2(self.get_gl_context_parent())
		
		if self.update_sym_dl:
			self.make_sym_dl_list(self.points,self.point_colors,self.eulers)
			self.update_sym_dl = False
		
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
		
		
		if ( self.sym_dl == 0 or self.force_update):
			self.generate_current_display_list()
			self.force_update = False
			if ( self.sym_dl == 0 ) : 
				print "error, you can't draw an empty list"
				return
		
		if ( self.arc_dl != 0 and self.display_arc ):
			glColor(.2,.9,.2)
			# this is a nice light blue color (when lighting is on)
			# and is the default color of the frame
			#glMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,(.1,.3,.1,1.0))
			#glMaterial(GL_FRONT,GL_SPECULAR,(.8,.8,.8,1.0))
			#glMaterial(GL_FRONT,GL_SHININESS,40.0)
			self.green()
			glCallList(self.arc_dl)
			
			if ( self.sym_object.is_h_sym() ):
				a = {}
				a["daz"] = 60
				dz = 5
				a["tz"] = dz
				self.sym_object.insert_params(a)
			if self.inspector != None and self.inspector.sym_toggled():
				for i in range(1,self.sym_object.get_nsym()):
					t = self.sym_object.get_sym(i)
					t.invert()
					d = t.get_rotation("eman")
					glPushMatrix()
					if ( self.sym_object.is_h_sym() ):
						trans = t.get_trans()
						glTranslatef(trans[0],trans[1],trans[2])
					# negated because OpenGL uses the opposite handed coordinate system
					glRotate(-d["phi"],0,0,1)
					glRotate(-d["alt"],1,0,0)
					glRotate(-d["az"],0,0,1)
					glCallList(self.arc_dl)
					glPopMatrix()
		
		if ( self.tri_dl != 0 and self.display_tri ):
			if ( self.sym_object.is_h_sym() != True ):
				glCallList(self.tri_dl)
			
				if self.inspector != None and self.inspector.sym_toggled():
					for i in range(1,self.sym_object.get_nsym()):
						t = self.sym_object.get_sym(i)
						t.invert()
						d = t.get_rotation("eman")
						glPushMatrix()
						# negated because OpenGL uses the opposite handed coordinate system
						glRotate(-d["phi"],0,0,1)
						glRotate(-d["alt"],1,0,0)
						glRotate(-d["az"],0,0,1)
						glCallList(self.tri_dl)
						glPopMatrix()
	
		if self.display_euler and self.sym_object != None:
			glColor(.9,.2,.8)
			self.gold()
			glStencilFunc(GL_EQUAL,self.rank,0)
			glStencilOp(GL_KEEP,GL_KEEP,GL_REPLACE)
			glPushMatrix()
			glCallList(self.sym_dl)
			glPopMatrix()
			if ( self.sym_object.is_h_sym() ):
				a = {}
				a["daz"] = 60
				dz = 5
				a["tz"] = dz
				self.sym_object.insert_params(a)
		
			if self.inspector != None and self.inspector.sym_toggled():
				#for i in range(1,self.sym_object.get_nsym()):
					#t = self.sym_object.get_sym(i)
				for t in self.sym_object.get_syms():
					d = t.get_rotation("eman")
					glPushMatrix()
					if ( self.sym_object.is_h_sym() ):
						trans = t.get_trans()
						glTranslatef(trans[0],trans[1],-trans[2])
						glRotate(d["phi"],0,0,1)
#						glRotate(d["alt"],1,0,0)
#						glRotate(d["az"],0,0,1)
					else:
						# negated because OpenGL uses the opposite handed coordinate system
						glRotate(-d["phi"],0,0,1)
						glRotate(-d["alt"],1,0,0)
						glRotate(-d["az"],0,0,1)
					glCallList(self.sym_dl)
					glPopMatrix()
			
			
		if self.trace_dl != 0:
			glColor(.9,.2,.8)
			# this is a nice light blue color (when lighting is on)
			# and is the default color of the frame
			self.black()
			glStencilFunc(GL_EQUAL,self.rank,0)
			glStencilOp(GL_KEEP,GL_KEEP,GL_REPLACE)
			#print "rendering trace"
			glPushMatrix()
			glCallList(self.trace_dl)
			glPopMatrix()
		glStencilFunc(GL_EQUAL,self.rank,self.rank)
		glStencilOp(GL_KEEP,GL_KEEP,GL_KEEP)
		
		#glPushMatrix()
		## FIXME the approach here is very inefficient
		#glLoadIdentity()
		#[width,height] = self.parent.get_near_plane_dims()
		#z = self.parent.get_start_z()
		#glTranslate(-width/2.0,-height/2.0,-z-0.01)
		#glScalef(width,height,1.0)
		#self.draw_bc_screen()
		#glPopMatrix()
		
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
			
	def resizeEvent(self,width=0,height=0):
		pass
	
		#self.vdtools.set_update_P_inv()
		
	def toggle_sym_display(self,bool):
		self.display_euler = bool
		self.updateGL()
		
	def triangletog(self,bool):
		self.display_tri = bool
		self.updateGL()
		
	def arctog(self,bool):
		self.display_arc = bool
		self.updateGL()

	def reducetog(self,bool):
		self.reduce = bool

class EMSymViewerWidget(QtOpenGL.QGLWidget,EMEventRerouter,EMGLProjectionViewMatrices):
	
	allim=WeakKeyDictionary()
	def __init__(self, em_slice_viwer):

		assert(isinstance(em_slice_viwer,EM3DSymViewerModule))
		
		EMSymViewerWidget.allim[self]=0
		
		fmt=QtOpenGL.QGLFormat()
		fmt.setDoubleBuffer(True)
		fmt.setDepth(1)
		fmt.setSampleBuffers(True)
		QtOpenGL.QGLWidget.__init__(self,fmt)
		EMEventRerouter.__init__(self,em_slice_viwer)
		EMGLProjectionViewMatrices.__init__(self)
		
		
		self.fov = 50 # field of view angle used by gluPerspective
		self.startz = 1
		self.endz = 5000
		self.cam = Camera()
		self.cam.setCamTrans("default_z",-100)
		
		self.resize(640,640)

	def initializeGL(self):
		glEnable(GL_NORMALIZE)
		glEnable(GL_LIGHT0)
		glEnable(GL_DEPTH_TEST)
		#print "Initializing"
		glLightfv(GL_LIGHT0, GL_AMBIENT, [0.9, 0.9, 0.9, 1.0])
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_POSITION, [0.5,0.7,11.,0.])
		GL.glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST)
		glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,GL_TRUE)
		GL.glClearColor(0,0,0,0)
		#GL.glClearAccum(0,0,0,0)
	
		glShadeModel(GL_SMOOTH)
		
		#glClearStencil(0)
		#glEnable(GL_STENCIL_TEST)
		
	def paintGL(self):
		#glClear(GL_ACCUM_BUFFER_BIT)
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT )
		
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		
		self.cam.position()
		
		glPushMatrix()
		self.target().render()
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
		self.target().resizeEvent()

	def get_start_z(self):
		return self.startz
	
	def get_near_plane_dims(self):
		height = 2.0 * self.get_start_z()*tan(self.fov/2.0*pi/180.0)
		width = self.aspect * height
		return [width,height]

	def show_inspector(self,force=0):
		self.target().show_inspector(self,force)

	def get_render_dims_at_depth(self,depth):
		# This function returns the width and height of the renderable 
		# area at the origin of the data volume
		height = -2*tan(self.fov/2.0*pi/180.0)*(depth)
		width = self.aspect*height
		
		return [width,height]
	
	
class EMSymInspector(QtGui.QWidget):
	def __init__(self,target) :
		QtGui.QWidget.__init__(self,None)
		self.setWindowIcon(QtGui.QIcon(get_image_directory() + "eulerxplor.png"))
		self.target=weakref.ref(target)
		self.rotation_sliders = EMTransformPanel(self.target(),self)
		
		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		
		self.hbl = QtGui.QHBoxLayout()
		self.hbl.setMargin(0)
		self.hbl.setSpacing(6)
		self.hbl.setObjectName("hbl")
		self.vbl.addLayout(self.hbl)
		
		
		self.vbl2 = QtGui.QVBoxLayout()
		self.vbl2.setMargin(0)
		self.vbl2.setSpacing(6)
		self.vbl2.setObjectName("vbl2")
		self.hbl.addLayout(self.vbl2)
	
		#self.cubetog = QtGui.QPushButton("Cube")
		#self.cubetog.setCheckable(1)
		#self.vbl2.addWidget(self.cubetog)
		
		self.button_hbl1 = QtGui.QHBoxLayout()
		
		self.symtogdisplay = QtGui.QPushButton("Display Eulers")
		self.symtogdisplay.setCheckable(1)
		self.symtogdisplay.setChecked(1)
		self.button_hbl1.addWidget(self.symtogdisplay)
		
		self.triangletog = QtGui.QPushButton("Display Triangles")
		self.triangletog.setCheckable(1)
		self.triangletog.setChecked(0)
		self.button_hbl1.addWidget(self.triangletog)
		
		self.vbl2.addLayout(self.button_hbl1)
		
		self.button_hbl2 = QtGui.QHBoxLayout()
		
		self.arctog = QtGui.QPushButton("Display Arcs")
		self.arctog.setCheckable(1)
		self.arctog.setChecked(1)
		self.button_hbl2.addWidget(self.arctog)
		
		self.symtog = QtGui.QPushButton("All syms")
		self.symtog.setCheckable(1)
		self.button_hbl2.addWidget(self.symtog)
		self.vbl2.addLayout(self.button_hbl2)
		
		
		self.vbl.addWidget(self.get_main_tab())
		
		self.n3_showing = False
		
		self.file = None
		self.lr = -1
		self.hr = -1
		
#		QtCore.QObject.connect(self.glcontrast, QtCore.SIGNAL("valueChanged"), target.set_GL_contrast)
#		QtCore.QObject.connect(self.glbrightness, QtCore.SIGNAL("valueChanged"), target.set_GL_brightness)
		QtCore.QObject.connect(self.sym_combo, QtCore.SIGNAL("currentIndexChanged(QString)"), self.sym_changed)
		QtCore.QObject.connect(self.sym_text, QtCore.SIGNAL("editingFinished()"), target.regen_dl)
		QtCore.QObject.connect(self.prop_text, QtCore.SIGNAL("editingFinished()"), target.regen_dl)
		QtCore.QObject.connect(self.mirror_checkbox, QtCore.SIGNAL("stateChanged(int)"), target.regen_dl)
		QtCore.QObject.connect(self.perturbtog, QtCore.SIGNAL("toggled(bool)"), self.perturbtoggled)
		#QtCore.QObject.connect(self.cubetog, QtCore.SIGNAL("toggled(bool)"), target.toggle_cube)
		QtCore.QObject.connect(self.symtog, QtCore.SIGNAL("toggled(bool)"), target.updateGL)
		QtCore.QObject.connect(self.symtogdisplay, QtCore.SIGNAL("clicked(bool)"), target.toggle_sym_display)
		QtCore.QObject.connect(self.triangletog, QtCore.SIGNAL("clicked(bool)"), target.triangletog)
		QtCore.QObject.connect(self.arctog, QtCore.SIGNAL("clicked(bool)"), target.arctog)
		QtCore.QObject.connect(self.tracetog, QtCore.SIGNAL("clicked(bool)"), self.toggle_trace)
		QtCore.QObject.connect(self.reducetog, QtCore.SIGNAL("clicked(bool)"), self.target().reducetog)
		QtCore.QObject.connect(self.lowrange, QtCore.SIGNAL("editingFinished()"), self.trace_update)
		QtCore.QObject.connect(self.highrange, QtCore.SIGNAL("editingFinished()"), self.trace_update)
		QtCore.QObject.connect(self.tracefile, QtCore.SIGNAL("editingFinished()"), self.trace_update)
		#QtCore.QObject.connect(self.cbb, QtCore.SIGNAL("currentIndexChanged(QString)"), target.setColor)
		QtCore.QObject.connect(self.angle_label, QtCore.SIGNAL("currentIndexChanged(QString)"), self.angle_label_changed)
		QtCore.QObject.connect(self.orient_label, QtCore.SIGNAL("currentIndexChanged(QString)"), self.orient_label_changed)
	
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
	
	
	def sym_toggled(self):
		return self.symtog.isChecked()
	
	def perturbtoggled(self,val):
		self.target().regen_dl()
		
	def angle_label_changed(self,string):
		self.target().regen_dl()
		
	def orient_label_changed(self,string):
		label = string
		self.angle_label.clear()
		if label == 'random':
			self.angle_label.addItem('n')
		else:
			self.angle_label.addItem('delta')
			self.angle_label.addItem('n')
			
		self.target().regen_dl()

	def sym_changed(self, sym):
		if sym == ' D ' or sym == ' C ' or sym == ' H ':
			self.sym_text.setEnabled(True)
		else:
			self.sym_text.setEnabled(False)
		
		self.target().regen_dl()

	def set_sym(self,sym):
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
			
			if s[0] == "d":
				self.sym_combo.setCurrentIndex(3)
			elif s[0] == "c":
				self.sym_combo.setCurrentIndex(4)
			elif s[0] == "h":
				self.sym_combo.setCurrentIndex(5)
			else:
				print "can't interpret",sym
				return
			
			self.sym_text.setText(s[1:])
			self.target().regen_dl()
		return
		

	def get_sym(self):
		sym = self.sym_map[str(self.sym_combo.currentText())]
		if sym in ['c','d','h']:
			sym = sym+self.sym_text.displayText()
		return sym
	
	def get_angle_label(self):
		return self.angle_label.currentText()
	
	def get_orient_label(self):
		return self.orient_label.currentText()
	
	def get_prop(self):
		return float(self.prop_text.displayText())
	
	def get_mirror(self):
		return self.mirror_checkbox.checkState() == Qt.Checked
	
	def get_perturb(self):
		return self.perturbtog.isChecked()

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

	def get_main_tab(self):
	
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
		self.symmetries.append(' H ')
		self.sym_map = {}
		self.sym_map[" Icosahedral "] = "icos"
		self.sym_map[" Octahedral "] = "oct"
		self.sym_map[" Tetrahedral "] = "tet"
		self.sym_map[" D "] = "d"
		self.sym_map[" C "] = "c"
		self.sym_map[" H "] = "h"
		idx_default = 0
		
		default_sym = self.target().sym
		
		for idx,i in enumerate(self.symmetries): 
			
			self.sym_combo.addItem(i)
			if self.sym_map[i][0] == default_sym[0]:
				idx_default = idx
				
		self.sym_combo.setCurrentIndex(idx_default)
		self.hbl_sym.addWidget(self.sym_combo)
		
		self.sym_label = QtGui.QLabel()
		self.sym_label.setText('C/D/H sym')
		self.hbl_sym.addWidget(self.sym_label)
		
		self.pos_int_validator = QtGui.QIntValidator(self)
		self.pos_int_validator.setBottom(1)
		self.sym_text = QtGui.QLineEdit(self)
		self.sym_text.setValidator(self.pos_int_validator)
		self.sym_text.setText(str(self.target().prop))
		self.sym_text.setFixedWidth(50)
		self.hbl_sym.addWidget(self.sym_text)
		self.sym_text.setEnabled(False)
		
		self.angle_label = QtGui.QComboBox(self)
		self.angle_label.addItem('delta')
		self.angle_label.addItem('n')
		self.hbl_sym.addWidget(self.angle_label)
		
		self.pos_double_validator = QtGui.QDoubleValidator(self)
		self.pos_double_validator.setBottom(0.05)
		self.prop_text = QtGui.QLineEdit(self)
		self.prop_text.setValidator(self.pos_double_validator)
		self.prop_text.setText("5.0")
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
		
		
		self.orient_label = QtGui.QComboBox(self)
		l = dump_orientgens_list()
		
		n = len(l)
		for i in l:
			self.orient_label.addItem(str(i))
		
		self.orient_label.setCurrentIndex(n-1)
		self.hbl_sym2.addWidget(self.orient_label)

		#here		
		self.hbl_pt = QtGui.QHBoxLayout()
		self.hbl_pt.setMargin(0)
		self.hbl_pt.setSpacing(6)
		self.hbl_pt.setObjectName("Ptl Trace")
		
		
		self.tracetog = QtGui.QPushButton("Trace")
		self.tracetog.setCheckable(1)
		self.tracetog.setChecked(0)
		self.hbl_pt.addWidget(self.tracetog)
		
		self.tracefile = QtGui.QLineEdit(self)
		self.tracefile.setText("filename.txt")
		self.tracefile.setFixedWidth(100)
		self.hbl_pt.addWidget(self.tracefile)
		self.tracefile.setEnabled(False)
		
		self.pt_label = QtGui.QLabel()
		self.pt_label.setText('Range')
		self.hbl_pt.addWidget(self.pt_label)
		
		self.pos_int_validator2 = QtGui.QIntValidator(self)
		self.pos_int_validator2.setBottom(0)
		self.lowrange = QtGui.QLineEdit(self)
		self.lowrange.setValidator(self.pos_int_validator2)
		self.lowrange.setText("1")
		self.lowrange.setFixedWidth(50)
		self.hbl_pt.addWidget(self.lowrange)
		self.lowrange.setEnabled(False)
		
		self.pt_label_to = QtGui.QLabel()
		self.pt_label_to.setText('to')
		self.hbl_pt.addWidget(self.pt_label_to)
		
		self.highrange = QtGui.QLineEdit(self)
		self.highrange.setValidator(self.pos_int_validator2)
		self.highrange.setText("1")
		self.highrange.setFixedWidth(50)
		self.hbl_pt.addWidget(self.highrange)
		self.highrange.setEnabled(False)
		
		self.reducetog = QtGui.QPushButton("Reduce")
		self.reducetog.setCheckable(1)
		self.reducetog.setChecked(0)
		self.hbl_pt.addWidget(self.reducetog)
		
		maintab.vbl.addLayout(self.hbl_pt)
#		#end here
#		
		self.mirror_checkbox = QtGui.QCheckBox("Mirror")
		self.hbl_sym2.addWidget(self.mirror_checkbox)
		
		self.perturbtog = QtGui.QPushButton("Perturb")
		self.perturbtog.setCheckable(1)
		self.hbl_sym2.addWidget(self.perturbtog)
		
#		self.glcontrast = ValSlider(maintab,(1.0,5.0),"GLShd:")
#		self.glcontrast.setObjectName("GLShade")
#		self.glcontrast.setValue(1.0)
#		maintab.vbl.addWidget(self.glcontrast)
#		
#		self.glbrightness = ValSlider(maintab,(-1.0,0.0),"GLBst:")
#		self.glbrightness.setObjectName("GLBoost")
#		self.glbrightness.setValue(0.1)
#		self.glbrightness.setValue(0.0)
#		maintab.vbl.addWidget(self.glbrightness)
		
		#self.cbb = QtGui.QComboBox(maintab)
		#self.vbl.addWidget(self.cbb)
		#self.cbb.deleteLater()

		self.rotation_sliders.addWidgets(maintab.vbl)
		
		return maintab

	def slider_rotate(self):
		self.target().load_rotation(self.get_current_rotation())
	
	def set_hist(self,hist,minden,maxden):
		self.hist.set_data(hist,minden,maxden)
	
if __name__ == '__main__':
	em_app = EMStandAloneApplication()
	window = EM3DSymViewerModule(application=em_app)
	em_app.show()
	em_app.execute()
