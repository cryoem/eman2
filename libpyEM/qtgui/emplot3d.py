#!/usr/bin/env python

#
# Author: David Woolford 10/26/2007 (woolford@bcm.edu)
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
from OpenGL.GLUT import *
from valslider import ValSlider
from math import *
from EMAN2 import *
import sys
import numpy
from emimageutil import ImgHistogram,EMParentWin
from weakref import WeakKeyDictionary
import weakref
from time import time
from PyQt4.QtCore import QTimer

from time import *
import copy

from emglobjects import EMImage3DGUIModule, Camera2,get_default_gl_colors,EMViewportDepthTools2,get_RGB_tab,get_gl_lights_vector,draw_volume_bounds,init_glut
from emlights import *
from emimageutil import EMTransformPanel

MAG_INCREMENT_FACTOR = 1.1

	
class EMPlot3DModule(EMLightsDrawer,EMImage3DGUIModule):
	def eye_coords_dif(self,x1,y1,x2,y2,mdepth=True):
		return self.vdtools.eye_coords_dif(x1,y1,x2,y2,mdepth)
	
	def __init__(self, application=None,parent=None):
		EMImage3DGUIModule.__init__(self,application)
		EMLightsDrawer.__init__(self)
		self.parent = parent
		
		self.mmode = 0
		self.wire = False
		self.light = True
		
		self.cam=Camera2(self)
		
		self.brightness = 0
		self.contrast = 10
		self.glcontrast = 1.0
		self.glbrightness = 0.0
		self.rank = 1
		self.inspector=None
		
		self.vdtools = EMViewportDepthTools2(self.gl_context_parent)
		
		self.gl_context_parent.cam.default_z = -25	 # this is me hacking
		self.gl_context_parent.cam.cam_z = -25 # this is me hacking
	
		self.highresspheredl = 0 # display list i
		self.draw_dl = 0
		self.clear_data() # use this so there's only one function
		self.init_plot_lights()
		glClearColor(1,1,1,1)
		
		self.draw_data_cube = True
		
		self.sphere_scale = 5
		
		self.allowable_shapes =  ["Sphere","Cube","Dodecahedron","Icosahedron","Tetrahedron","Cone","Octahedron","Teapot"]
		self.init_color_themes()
		
	def set_sphere_scale(self,value):
		if value != 0 and value != self.sphere_scale:
			self.sphere_scale = value
			self.full_refresh()
			self.updateGL()
	
	def init_plot_lights(self):
		glLightfv(GL_LIGHT0, GL_AMBIENT, [0.1, 0.1, 0.1, 1.0])
		glLightfv(GL_LIGHT0, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0])
		glLightfv(GL_LIGHT0, GL_POSITION, [0.02,-0.02,1.,0.])
		 
	def get_type(self):
		return "lights"
	
	def clear_data(self):
		self.data = {} # this is the data sets
		self.ndata = {} # this is the normalized data sets (values go from 0 to 1)
		self.visibility = {} # this is the visibility of the each of the data sets
		self.axes = {} # This stores which axes are currently being viewed, for each data set
		self.axis_colors = {} # Stores colors assigned to axes
		self.shapes = {} # Stores the shapes used draw scatter plots
		self.plot_type = {} # This stores whether the plot is being treated as 2D
		self.min = {} # Stores the minimum value in each data set
		self.max = {} # Stores the maximum value in each data set
		
	
	
	def set_data_from_file(self,filename,clear_current=False):
		"""Reads a keyed data set from a file. Automatically interpret the file contents."""
#		self.del_shapes()
		
		if clear_current: 
			self.data = {}
			self.axes = {}
			self.visibility = {}
		
		file_type = Util.get_filename_ext(filename)
		em_file_type = EMUtil.get_image_ext_type(file_type)
		
		if em_file_type != IMAGE_UNKNOWN:
			
			im=EMData.read_images(filename)
			if len(im) == 1:
				im = im[0]
				l = [[i for i in range(im.get_xsize())]]
				k = im.get_data_as_vector()
				i = 0
				while i < len(k):
					l.append(k[i:i+im.get_xsize()])
					i += im.get_xsize()
				self.set_data(filename,l)
			elif im[0].get_attr_default("isvector",0):
#				all=[i.get_data_as_vector() for i in im]

				all=[]
				for j in range(im[0].get_xsize()):
					r=[]
					for i in range(len(im)):
						r.append(im[i][j,0])
					all.append(r)
				self.set_data("Vecset",all)
			else:
				for idx,image in enumerate(im):
					l = [i for i in range(image.get_size())]
					k = image.get_data_as_vector()
					i = 0
					while i < len(k):
						l.append(k[i:i+image.get_xsize()])
						i += image.get_xsize() 
					self.set_data("Image "+str(idx),l)
				
		elif file_type == 'fp':
			fin=file(filename)
			fph=struct.unpack("120sII",fin.read(128))
			ny=fph[1]
			nx=fph[2]
			data=[]
			for i in range(nx):
				data.append(struct.unpack("%df"%ny,fin.read(4*ny)))
				
			self.set_data(filename,data)
		else:
			try:
				fin=file(filename)
				fin.seek(0)
				rdata=fin.readlines()
				rdata=[i for i in rdata if i[0]!='#']
				if ',' in rdata[0]: rdata=[[float(j) for j in i.split(',')] for i in rdata]
				else : rdata=[[float(j) for j in i.split()] for i in rdata]
				nx=len(rdata[0])
				ny=len(rdata)
				data=[[rdata[j][i] for j in range(ny)] for i in range(nx)]
				
				self.set_data(data,remove_directories_from_name(filename))
			except:
				traceback.print_exc()
				print "couldn't read",filename
				return False
				
		return True
	
	def parse_txt_file(filename):
		try:
			fin=file(filename)
			fin.seek(0)
			rdata=fin.readlines()
			rdata=[i for i in rdata if i[0]!='#']
			if ',' in rdata[0]: rdata=[[float(j) for j in i.split(',')] for i in rdata]
			else : rdata=[[float(j) for j in i.split()] for i in rdata]
			nx=len(rdata[0])
			ny=len(rdata)
			data=[[rdata[j][i] for j in range(ny)] for i in range(nx)]
			return data
		except:
			return None
		
	parse_txt_file = staticmethod(parse_txt_file)
	
	def set_data(self,data,key="data",clear_current=False,shape="Sphere"):
		'''
		
		'''
		
		if len(data) < 2:
			print "error, the length of the input data must be atleast 2"
			return
		l = -1
		for d in data:
			if l == -1: l = len(d)
			else:
				if l != len(d):
					print "error, the axis data must all be the same length"
					return
		
		if clear_current: 
			self.clear_data() # use this so there's only one function
		
		# note no checking, if the old data exists it is simply replaced
		self.data[key] = data
		self.visibility[key] = True
		if len(data) == 2:
			self.axes[key] = [0,1]
		else: # the data must have atleast 3 axes, this is guaranteed by the checking functions at the entry point into this function 
			self.axes[key] = [0,1,2]
			
			keys = self.colors.keys()
			n = len(keys)
			
			
			self.axis_colors[key] = "Color cube"
			self.shapes[key] = shape
			
		min  = []
		max = []
		for d in data:
			s = copy.copy(d)
			s.sort()
			min.append(s[0])
			max.append(s[-1])
		self.min[key] = min
		self.max[key] = max


	def init_color_themes(self):
		self.color_themes = {}
		self.color_themes["Splash"] = ["bluewhite","green","blue","cyan","purple","yellow","black","white"]
		self.color_themes["Color cube"] = ["black","blue","green","cyan","red","purple","yellow","white"]
		self.color_themes["Pastel"] = ["yellow","red","orange","pearl","green","blue","bluewhite","purple"]
		self.color_themes["Pluto"] = ["red","pearl","blue","cyan","white","black","pearl","white"]
		self.color_themes["Indigo"] = ["black","purple","blue","green","turquoise","bluewhite","cyan","white"]
		self.color_themes["Sun rise"] = ["blue","gold","orange","bluewhite","silver","red","pearl","white"]
		self.color_themes["Grey Scale"] = ["black","black","black","black","white","white","white","white"]
		self.color_themes["Red Scale"] = ["red","red","red","red","white","white","white","white"]
		self.color_themes["Blue Scale"] = ["blue","blue","blue","blue","white","white","white","white"]
		self.color_themes["Green Scale"] = ["emerald","emerald","emerald","emerald","white","white","white","white"]
		
		self.basic_color_themes = ["red","blue","green","purple","cyan","yellow","black","orange","white"]
		
		keys = self.colors.keys()
		n = len(keys)
		
		self.color_themes["Random"] = [keys[Util.get_irand(0,n-1)] for i in range(8) ]
		
		self.global_theme = "Color cube"
		self.using_global_theme = False

	def render(self):

		#if (not isinstance(self.data,EMData)): return
		lighting = glIsEnabled(GL_LIGHTING)
		cull = glIsEnabled(GL_CULL_FACE)
		depth = glIsEnabled(GL_DEPTH_TEST)
		polygonmode = glGetIntegerv(GL_POLYGON_MODE)

		glDisable(GL_CULL_FACE)
		glEnable(GL_DEPTH_TEST)
		
		glShadeModel(GL_SMOOTH)

		glStencilFunc(GL_EQUAL,self.rank,0)
		glStencilOp(GL_KEEP,GL_KEEP,GL_REPLACE)
		
		
		if ( self.wire ):
			glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
		else:
			glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
		
		if self.light:
			glEnable(GL_LIGHTING)
		else:
			glDisable(GL_LIGHTING)

		glPushMatrix()
		self.cam.position()
			
		
		glEnable(GL_NORMALIZE)
		#HERE
		
		if self.highresspheredl == 0:
			self.highresspheredl=glGenLists(1)
			glNewList(self.highresspheredl,GL_COMPILE)
			gluSphere(self.gq,.5,12,8)
			glEndList()
		
		if ( self.draw_dl == 0 ):
			self.draw_dl=glGenLists(1)
				
			glNewList(self.draw_dl,GL_COMPILE)
			self.scatter_plot()
							
			glEndList(self.draw_dl)
			
		glCallList(self.draw_dl)
		
		glPopMatrix()
		EMLightsDrawer.draw(self)
		
		self.draw_bc_screen()


		if ( lighting ): glEnable(GL_LIGHTING)
		else: glDisable(GL_LIGHTING)
		if ( cull ): glEnable(GL_CULL_FACE)
		else: glDisable(GL_CULL_FACE)
		if ( depth ): glEnable(GL_DEPTH_TEST)
		else : glDisable(GL_DEPTH_TEST)
		
		if ( polygonmode[0] == GL_LINE ): glPolygonMode(GL_FRONT, GL_LINE)
		else: glPolygonMode(GL_FRONT, GL_FILL)
		if ( polygonmode[1] == GL_LINE ): glPolygonMode(GL_BACK, GL_LINE)
		else: glPolygonMode(GL_BACK, GL_FILL)
		
	def full_refresh(self):
		glDeleteLists(self.draw_dl,1)
		self.draw_dl = 0
		
	def scatter_plot(self):
		'''
		Will do a scatter plot using the available data
		If any data is 2D it will simply ignore it
		'''
		init_glut()
	
		minx,miny,minz,maxx,maxy,maxz = None,None,None,None,None,None
		mins = [None,None,None]
		maxs = [None,None,None]
	
		keys = self.data.keys()
		vis_keys = []
		for key in keys:
			if self.visibility[key]:
				min = self.min[key]
				max = self.max[key]
				axes = self.axes[key]
				
				for i,m in enumerate(mins):
					if m == None or min[axes[i]] < m: mins[i] = min[axes[i]]
					
				for i,m in enumerate(maxs):
					if m == None or max[axes[i]] >  m: maxs[i] = max[axes[i]]
				vis_keys.append(key)
				
		if len(vis_keys) == 0: return # there are no visible plots
		
		minx,miny,minz = mins[0],mins[1],mins[2]
		maxx,maxy,maxz = maxs[0],maxs[1],maxs[2]
		
		width = maxx-minx
		height = maxy-miny
		depth = maxz-minz
		
		if width == 0: width = 1
		if height == 0: height = 1
		if depth == 0: depth = 1
		glPushMatrix()
		glTranslate(-(maxx+minx)/2,-(maxy+miny)/2,-(maxz+minz)/2) # centers the plot
		
		f = 1/(2.0*sqrt(3.0))
		g = 1/(sqrt(3.0))
		for key in vis_keys:
			d = self.data[key]
			a = self.axes[key]
			mn = self.min[key]
			mx = self.max[key]
			x = d[a[0]]
			y = d[a[1]]
			if len(a) == 3: z = d[a[2]]
			else: z = None

			if z != None:
				#
				for i in range(len(x)):
					
					if not self.using_global_theme:
						color = self.get_color_interp(self.axis_colors[key],(x[i]-mn[0])/float(mx[0]-mn[0]), (y[i]-mn[1])/float(mx[1]-mn[1]),(z[i]-mn[2])/float(mx[2]-mn[2]))
					else:
						color = self.get_color_interp(self.global_theme,(x[i]-minx)/float(width), (y[i]-miny)/float(height),(z[i]-minz)/float(depth))
			
					glMaterial(GL_FRONT, GL_AMBIENT, color["ambient"])
					glMaterial(GL_FRONT, GL_DIFFUSE, color["diffuse"])
					glMaterial(GL_FRONT, GL_SPECULAR, color["specular"])
					glMaterial(GL_FRONT, GL_EMISSION, color["emission"])
					glMaterial(GL_FRONT, GL_SHININESS, color["shininess"])
					glColor(color["ambient"])
					
					
					glPushMatrix()
					glTranslate(x[i],y[i],z[i])
					glScale(self.sphere_scale,self.sphere_scale,self.sphere_scale)
						
					shape = self.shapes[key]
						
					if shape == "Sphere":
						glCallList(self.highresspheredl)
					elif shape == "Dodecahedron":
						glScale(f,f,f)
						glutSolidDodecahedron()
					elif shape == "Icosahedron":
						glScale(0.5,0.5,0.5)
						glutSolidIcosahedron()
					elif shape ==  "Cube":
						glutSolidCube(1.0)
					elif shape == "Tetrahedron":
						glScale(g,g,g)
						glutSolidTetrahedron()
					elif shape == "Cone":
						glutSolidCone(0.5,0.5,12,6)
						glRotate(180,1,0,0)
						glutSolidCone(0.5,0.5,12,6)
					elif shape == "Octahedron":
						glScale(0.5,0.5,0.5)
						glutSolidOctahedron()
					elif shape == "Teapot":
						glutSolidTeapot(.5)
					else:
						print "Unknown shape:", shape 
					
					glPopMatrix()
					
			if self.draw_data_cube:
				
				light_on = glIsEnabled(GL_LIGHTING)
				glDisable(GL_LIGHTING)
				glPushMatrix()
				glTranslate(minx,miny,minz)
				glColor(0,0,0,1)
				draw_volume_bounds(width,height,depth,False)
				glPopMatrix()
				if light_on: glEnable(GL_LIGHTING)
							
							
		glPopMatrix()
		
	def get_color_interp(self,theme,x,y,z):
		types = ["ambient","diffuse","specular","emission","shininess"]
		p = [x,y,z]
		q = [1-c for c in p]
		
		if theme in self.colors.keys():
			return self.colors[theme]
		
		colors = [self.colors[v] for v in self.color_themes[theme]]
		#[c000,c100,c010,c110,c001,c101,c011,c111]
		factors = [q[0]*q[1]*q[2],p[0]*q[1]*q[2],q[0]*p[1]*q[2],p[0]*p[1]*q[2]]
		factors.extend([q[0]*q[1]*p[2],p[0]*q[1]*p[2],q[0]*p[1]*p[2],p[0]*p[1]*p[2]])
		
		result = {}
		for gl_color_type in types:
			r = []
			if gl_color_type != "shininess":
				n = len(colors[0][gl_color_type])
				for j in range(n):
					v = 0
					for i,color in enumerate(colors):
						v += factors[i]*color[gl_color_type][j]
					
					r.append(v)
				result[gl_color_type] = r
			else:
				v = 0
				for i,color in enumerate(colors):
					v += factors[i]*color[gl_color_type]
			
				result[gl_color_type] = v

		#print result
		return result
					
	def get_color(self,xcolor,ycolor,zcolor,xcolor_begin,ycolor_begin,zcolor_begin,x,y,z):
		
		result = {}
		t = ["ambient","diffuse","specular","emission","shininess"]
		colors = [xcolor,ycolor,zcolor]
		begin_colors = [xcolor_begin,ycolor_begin,zcolor_begin]
		coord = [x,y,z]
		inv_coord = [1-c or c in coord]
		
		
		fact = 1/3.0
		for c in t:
			r = []
			if c != "shininess":
				n = len(colors[0][c])
				for j in range(n):
					v = 0
					for i,color in enumerate(colors):
						v += coord[i]*color[c][j] +(1-coord[i])*begin_colors[i][c][j]
					
					r.append(v*fact)
				result[c] = r
			else:
				v = 0
				for i,color in enumerate(colors):
					v += (coord[i]*color[c])/float(len(colors))
			
				result[c] = v

		#print result
		return result
	
	def setInit(self):

		#self.cam.default_z = -1.25*32
		#self.cam.cam_z = -1.25*32
		
		if not self.inspector or self.inspector ==None:
			self.inspector=EMPlot3DInspector(self)
		
		self.load_colors()
	
	def toggle_wire(self,val):
		self.wire = not self.wire
		self.updateGL()
		
	def toggle_light(self,val):
		self.light = not self.light
		self.updateGL()
	
	def update_inspector(self,t3d):
		if not self.inspector or self.inspector ==None:
			self.inspector=EMPlot3DInspector(self)
		self.inspector.update_rotations(t3d)
	
	def get_inspector(self):
		if not self.inspector : 
			self.inspector=EMPlot3DInspector(self)
		return self.inspector

class EMPlot3DInspector(QtGui.QWidget,EMLightsInspectorBase):
	def __init__(self,target) :
		QtGui.QWidget.__init__(self,None)
		EMLightsInspectorBase.__init__(self)
		self.target=weakref.ref(target)
		
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
		
		self.wiretog = QtGui.QPushButton("Wire")
		self.wiretog.setCheckable(1)
		self.vbl2.addWidget(self.wiretog)
		
		self.lighttog = QtGui.QPushButton("Light")
		self.lighttog.setCheckable(1)
		self.lighttog.setChecked(True)
		self.vbl2.addWidget(self.lighttog)
		
		self.tabwidget = QtGui.QTabWidget()
		self.maintab = None
		self.tabwidget.addTab(self.get_plot_tab(), "Plot")
		self.tabwidget.addTab(self.get_main_tab(), "Transform")
		self.tabwidget.addTab(self.get_light_tab(), "Lights")
		self.tabwidget.addTab(self.get_GL_tab(),"GL")
		self.vbl.addWidget(self.tabwidget)
		self.n3_showing = False
		self.quiet = False
		
		self.data_change()

		QtCore.QObject.connect(self.wiretog, QtCore.SIGNAL("toggled(bool)"), target.toggle_wire)
		QtCore.QObject.connect(self.lighttog, QtCore.SIGNAL("toggled(bool)"), target.toggle_light)
		QtCore.QObject.connect(self.glcontrast, QtCore.SIGNAL("valueChanged"), target.set_GL_contrast)
		QtCore.QObject.connect(self.glbrightness, QtCore.SIGNAL("valueChanged"), target.set_GL_brightness)
		

	def update_rotations(self,t3d):
		self.rotation_sliders.update_rotations(t3d)
	
	def set_scale(self,val):
		self.rotation_sliders.set_scale(val)
	
	def set_xy_trans(self, x, y):
		self.rotation_sliders.set_xy_trans(x,y)
	
	def set_xyz_trans(self,x,y,z):
		self.rotation_sliders.set_xyz_trans(x,y,z)
		#return self.advanced_tab
		
	def get_plot_tab(self):
		self.plot_tab = QtGui.QWidget()
		plot_tab = self.plot_tab
		
		self.setlist=QtGui.QListWidget(self)
		self.setlist.setSizePolicy(QtGui.QSizePolicy.Preferred,QtGui.QSizePolicy.Expanding)
		
		vbl = QtGui.QVBoxLayout(self.plot_tab)
		vbl.addWidget(self.setlist)
		
		gl=QtGui.QGridLayout()
		gl.addWidget(QtGui.QLabel("X Col:",self),0,0,Qt.AlignRight)
		self.slidex=QtGui.QSpinBox(self)
		self.slidex.setRange(-1,1)
		gl.addWidget(self.slidex,0,1,Qt.AlignLeft)
		
		gl.addWidget(QtGui.QLabel("Y Col:",self),0,2,Qt.AlignRight)
		self.slidey=QtGui.QSpinBox(self)
		self.slidey.setRange(-1,1)
		gl.addWidget(self.slidey,0,3,Qt.AlignLeft)
		
		gl.addWidget(QtGui.QLabel("Z Col:",self),0,4,Qt.AlignRight)
		self.slidez=QtGui.QSpinBox(self)
		self.slidez.setRange(-1,1)
		gl.addWidget(self.slidez,0,5,Qt.AlignLeft)
		vbl.addLayout(gl)
		
		
	
		
		self.sphere_scale = ValSlider(self.plot_tab,(0.05,10.0),"Shape scale:")
		self.sphere_scale.setValue(self.target().sphere_scale)
		vbl.addWidget(self.sphere_scale)
		
		hbl1 = QtGui.QHBoxLayout()
		self.shape=QtGui.QComboBox(self)
		items = self.target().allowable_shapes
		for k in items: self.shape.addItem(k)
		
		hbl1.addWidget(QtGui.QLabel("Shape:",self))
		hbl1.addWidget(self.shape)
		vbl.addLayout(hbl1)
		
		
		hbl = QtGui.QHBoxLayout()
		self.theme=QtGui.QComboBox(self)
		keys = self.target().color_themes.keys()
		keys.sort()
		for k in keys: self.theme.addItem(k)
		self.theme.insertSeparator(len(keys))
		keys = self.target().colors.keys()
		for k in keys: self.theme.addItem(k)
		
		self.global_theme = QtGui.QCheckBox("Global theme")
		self.global_theme.setChecked(self.target().using_global_theme)
		
		
		hbl.addWidget(QtGui.QLabel("Color theme:",self))
		hbl.addWidget(self.theme)
		hbl.addWidget(self.global_theme)
		
		
		vbl.addLayout(hbl)
#
		QtCore.QObject.connect(self.slidex, QtCore.SIGNAL("valueChanged(int)"), self.new_data_cols)
		QtCore.QObject.connect(self.slidey, QtCore.SIGNAL("valueChanged(int)"), self.new_data_cols)
		QtCore.QObject.connect(self.slidez, QtCore.SIGNAL("valueChanged(int)"), self.new_data_cols)
		QtCore.QObject.connect(self.setlist,QtCore.SIGNAL("itemChanged(QListWidgetItem*)"),self.list_item_changed)
		QtCore.QObject.connect(self.setlist,QtCore.SIGNAL("currentRowChanged(int)"),self.plot_data_selected)
		QtCore.QObject.connect(self.sphere_scale, QtCore.SIGNAL("valueChanged"), self.target().set_sphere_scale)
		QtCore.QObject.connect(self.theme,QtCore.SIGNAL("currentIndexChanged(int)"),self.theme_changed)
		QtCore.QObject.connect(self.shape,QtCore.SIGNAL("currentIndexChanged(int)"),self.shape_changed)
		QtCore.QObject.connect(self.global_theme, QtCore.SIGNAL("stateChanged(int)"), self.global_theme_checked)
		
		return plot_tab
	
	def plot_data_selected(self,row):
		self.quiet=1
		try:
			i=str(self.setlist.item(row).text())
		except: 
#			print "plot error"
			return
		n = len(self.target().data[i])-1
		axes = self.target().axes[i]
		self.slidex.setRange(0,n)
		self.slidey.setRange(0,n)
		self.slidez.setRange(0,n)
		self.slidex.setValue(axes[0])
		self.slidey.setValue(axes[1])
		self.slidez.setValue(axes[2])
		
		self.__set_theme()
			
		shape = self.target().shapes[i]
		for j in range(self.shape.count()):
			if str(self.shape.itemText(j)) == shape:
				self.shape.setCurrentIndex(j)
				break

		self.quiet=0

	def __set_theme(self):
		if self.target().using_global_theme: theme = self.target().global_theme
		else: 
			items = self.setlist.selectedItems()
			if len(items) != 1:
				print "error, set list does not have a single selection. This should not happen"
			else:
				theme = self.target().axis_colors[str(items[0].text())]
		
		for j in range(self.theme.count()):
			if str(self.theme.itemText(j)) == theme:
				self.theme.setCurrentIndex(j)
				break
		
	
	def new_data_cols(self,i):
		if self.quiet: return
		if self.target(): 
			self.target().axes[str(self.setlist.currentItem().text())] = [int(self.slidex.value()),int(self.slidey.value()),int(self.slidez.value())]
			self.target().full_refresh()
			self.target().updateGL()
			
	def theme_changed(self,i):
		if self.quiet: return
		if self.target(): 
			if self.target().using_global_theme:
				self.target().global_theme = str(self.theme.itemText(i))
			else:
				self.target().axis_colors[str(self.setlist.currentItem().text())] = str(self.theme.itemText(i))
			
			self.target().full_refresh()
			self.target().updateGL()
			
	def global_theme_checked(self,i):
		self.target().using_global_theme = self.global_theme.isChecked()
		self.quiet = 1
		self.__set_theme()
		self.quiet = 0
		self.target().full_refresh()
		self.target().updateGL()
		
	
	def shape_changed(self,i):
		if self.quiet: return
		if self.target(): 
			self.target().shapes[str(self.setlist.currentItem().text())] = str(self.shape.itemText(i))
			self.target().full_refresh()
			self.target().updateGL()
			
	
	def data_change(self):
		
		self.setlist.clear()
		
		#flag1 = Qt.ItemFlags(Qt.ItemIsTristate)
		flag2 = Qt.ItemFlags(Qt.ItemIsSelectable)
		flag3 = Qt.ItemFlags(Qt.ItemIsEnabled)
		flag4 = Qt.ItemFlags(Qt.ItemIsUserCheckable)
		
		keys=self.target().data.keys()
		visible = self.target().visibility
		keys.sort()
				
		for i,j in enumerate(keys) :
			a = QtGui.QListWidgetItem(j)
			a.setFlags(flag2|flag3|flag4)
			if visible[j]: a.setCheckState(Qt.Checked)
			else: a.setCheckState(Qt.Unchecked)
			
			self.setlist.addItem(a)

		if len(keys) > 0 : self.setlist.setCurrentRow(0)
		
	def list_item_changed(self,item):
		checked = False
		if item.checkState() == Qt.Checked: checked = True
		
		name = str(item.text())
		if self.target().visibility[name] != checked:
			self.target().visibility[name] = checked
			self.target().full_refresh()
			self.target().updateGL()
	
	
	def get_GL_tab(self):
		self.gltab = QtGui.QWidget()
		gltab = self.gltab
		
		gltab.vbl = QtGui.QVBoxLayout(self.gltab )
		gltab.vbl.setMargin(0)
		gltab.vbl.setSpacing(6)
		gltab.vbl.setObjectName("Main")
		
		self.hbl_color = QtGui.QHBoxLayout()
		self.hbl_color.setMargin(0)
		self.hbl_color.setSpacing(6)
		gltab.vbl.addLayout(self.hbl_color)
		
		self.glcontrast = ValSlider(gltab,(1.0,5.0),"GLShd:")
		self.glcontrast.setObjectName("GLShade")
		self.glcontrast.setValue(1.0)
		gltab.vbl.addWidget(self.glcontrast)
		
		self.glbrightness = ValSlider(gltab,(-1.0,0.0),"GLBst:")
		self.glbrightness.setObjectName("GLBoost")
		self.glbrightness.setValue(0.1)
		self.glbrightness.setValue(0.0)
		gltab.vbl.addWidget(self.glbrightness)
	
		return gltab
	
	def get_main_tab(self):
		if ( self.maintab == None ):
			self.maintab = QtGui.QWidget()
			maintab = self.maintab
			maintab.vbl = QtGui.QVBoxLayout(self.maintab)
			maintab.vbl.setMargin(0)
			maintab.vbl.setSpacing(6)
			maintab.vbl.setObjectName("Main")
			
			self.rotation_sliders = EMTransformPanel(self.target(),self)
			self.rotation_sliders.addWidgets(maintab.vbl)
		
		return self.maintab

def get_test_data():
	'''
	Returns a test data in a format suitable for the 3D plotter
	'''
	data = []
	n = 50
	r = range(n)
	data.append([i for i in r])
	
	data.append([n/2.0*sin(i*pi/5) for i in r])
	
	#
	data.append([Util.get_frand(-n/2.0,n/2.0) for i in r])
	data.append([n/2.0*cos(i*pi/10) for i in r])
	
	data.append([n-i-1 for i in r])
	data.append([i**2 for i in r])
	data.append([sqrt(i) for i in r])
	
	return data


def get_other_test_data():
	'''
	Returns a test data in a format suitable for the 3D plotter
	'''
	data = []
	n = 50
	r = range(n)
	data.append([i+n*1.2 for i in r])
	#data.append([n-i-1 for i in r])
	#data.append([i**2 for i in r])
	data.append([n/2.0*sin(i*pi/5) for i in r])
	
	
	data.append([n/2.0*cos(i*pi/10) for i in r])
	
	data.append([sqrt(i) for i in r])
	data.append([Util.get_frand(-n/2.0,n/2.0) for i in r])
	
	return data

# This is just for testing, of course
if __name__ == '__main__':
	from emapplication import EMStandAloneApplication
	em_app = EMStandAloneApplication()
	window = EMPlot3DModule(application=em_app)
	window.set_data(get_test_data(),"test data")
	window.set_data(get_other_test_data(),"other data",shape="Cube")
	em_app.show()
	em_app.execute()
	
