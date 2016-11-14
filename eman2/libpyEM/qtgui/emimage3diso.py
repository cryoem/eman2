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
from time import time
import weakref
from time import *
from libpyGLUtils2 import GLUtil

from emglobjects import EMViewportDepthTools, Camera2, get_default_gl_colors,get_RGB_tab, EM3DModel
from emimageutil import ImgHistogram, EMTransformPanel


MAG_INCREMENT_FACTOR = 1.1

class EMIsosurfaceModel(EM3DModel):
	def eye_coords_dif(self,x1,y1,x2,y2,mdepth=True):
		return self.vdtools.eye_coords_dif(x1,y1,x2,y2,mdepth)

	def __init__(self,gl_widget, image=None,enable_file_browse=False):
		self.data = None
		EM3DModel.__init__(self, gl_widget)
		self.init()
		self.initialized = True
		
		self.cam=Camera2(self)
		self.tex_name = 0
		self.texture = False

		self.brightness = 0
		self.contrast = 10
		self.glcontrast = 1.0
		self.glbrightness = 0.0
		self.rank = 1
		self.inspector=None
		self.data_copy = None		
		self.vdtools = EMViewportDepthTools(self)
		self.enable_file_browse = enable_file_browse
		self.force_update = False
		if image :
			self.set_data(image)
		
#	def __del__(self):
#		print "iso died"
#		pass
	def set_force_update(self,val=True): self.force_update = val

	def get_type(self):
		return "Isosurface"
	
	def get_emit_signals_and_connections(self):
		return {"set_threshold":self.set_threshold}
	
	def update_data_and_texture(self):
		
		self.data_copy = self.data.copy()
		self.data_copy.add(self.brightness)
		self.data_copy.mult(self.contrast)
		
		hist = self.data_copy.calc_hist(256,self.minden,self.maxden)
		self.inspector.set_hist(hist,self.minden,self.maxden) 

		if ( self.texture ): self.gen_texture()
	
	def gen_texture(self):
		if ( self.texture == False ): return
		if ( self.tex_name != 0 ):
			glDeleteTextures(self.tex_name)
		
		if ( self.data_copy == None ):
			self.tex_name = GLUtil.gen_gl_texture(self.data)
		else:
			self.tex_name = GLUtil.gen_gl_texture(self.data_copy)
	
	def render(self):
		if (not isinstance(self.data,EMData)): return
		#a = time()
		lighting = glIsEnabled(GL_LIGHTING)
		cull = glIsEnabled(GL_CULL_FACE)
		depth = glIsEnabled(GL_DEPTH_TEST)
		polygonmode = glGetIntegerv(GL_POLYGON_MODE)
		normalize = glIsEnabled(GL_NORMALIZE)
		
		
		glEnable(GL_CULL_FACE)
		glCullFace(GL_BACK)
		#glDisable(GL_CULL_FACE)
		glEnable(GL_DEPTH_TEST)
		glEnable(GL_NORMALIZE)
		#glDisable(GL_NORMALIZE)
		if ( self.wire ):
			glPolygonMode(GL_FRONT,GL_LINE);
		else:
			glPolygonMode(GL_FRONT,GL_FILL);
		
		if self.light:
			glEnable(GL_LIGHTING)
		else:
			glDisable(GL_LIGHTING)

		
		glPushMatrix()
		self.cam.position(True)
		# the ones are dummy variables atm... they don't do anything
		self.vdtools.update(1,1)
		glPopMatrix()
		
		self.cam.position()
		glShadeModel(GL_SMOOTH)
		if ( self.isodl == 0 or self.force_update ):
			self.get_iso_dl()
			self.force_update = False
		glStencilFunc(GL_EQUAL,self.rank,0)
		glStencilOp(GL_KEEP,GL_KEEP,GL_REPLACE)
		glMaterial(GL_FRONT, GL_AMBIENT, self.colors[self.isocolor]["ambient"])
		glMaterial(GL_FRONT, GL_DIFFUSE, self.colors[self.isocolor]["diffuse"])
		glMaterial(GL_FRONT, GL_SPECULAR, self.colors[self.isocolor]["specular"])
		glMaterial(GL_FRONT, GL_SHININESS, self.colors[self.isocolor]["shininess"])
		glMaterial(GL_FRONT, GL_EMISSION, self.colors[self.isocolor]["emission"])
		glColor(self.colors[self.isocolor]["ambient"])
		glPushMatrix()
		glTranslate(-self.data.get_xsize()/2.0,-self.data.get_ysize()/2.0,-self.data.get_zsize()/2.0)
		if ( self.texture ):
			glScalef(self.data.get_xsize(),self.data.get_ysize(),self.data.get_zsize())
		glCallList(self.isodl)
		glPopMatrix()
		
		self.draw_bc_screen()
		
		glStencilFunc(GL_ALWAYS,1,1)
		if self.cube:
			glDisable(GL_LIGHTING)
			glPushMatrix()
			self.draw_volume_bounds()
			glPopMatrix()
			
		if ( lighting ): glEnable(GL_LIGHTING)
		else: glDisable(GL_LIGHTING)
		if ( not cull ): glDisable(GL_CULL_FACE)
		else: glDisable(GL_CULL_FACE)
		if ( depth ): glEnable(GL_DEPTH_TEST)
		else : glDisable(GL_DEPTH_TEST)
		
		if ( not normalize ): glDisable(GL_NORMALIZE)
		
		if ( polygonmode[0] == GL_LINE ): glPolygonMode(GL_FRONT, GL_LINE)
		else: glPolygonMode(GL_FRONT, GL_FILL)
		#if ( polygonmode[1] == GL_LINE ): glPolygonMode(GL_BACK, GL_LINE)
		#else: glPolygonMode(GL_BACK, GL_FILL)
		
		#print "total time is", time()-a
		
	def init(self):
		self.mmode = 0
		self.inspector=None
		self.isothr=0.5
		self.isorender=None
		self.isodl = 0
		self.smpval=-1
		self.griddl = 0
		self.scale = 1.0
		self.cube = False
		self.wire = False
		self.light = True
		
		
	def get_iso_dl(self):
		# create the isosurface display list
		self.isorender.set_surface_value(self.isothr)
		self.isorender.set_sampling(self.smpval)
		
		if ( self.texture ):
			if ( self.tex_name == 0 ):
				self.update_data_and_texture()
		
		face_z = False
		if self.data.get_zsize() <= 2:
			face_z = True
		
		if ( self.texture  ):
			self.isodl = GLUtil.get_isosurface_dl(self.isorender, self.tex_name,face_z)
		else:
			self.isodl = GLUtil.get_isosurface_dl(self.isorender, 0,face_z)
		#time2 = clock()
		#dt1 = time2 - time1
		#print "It took %f to render the isosurface" %dt1
	
	def update_data(self,data):
		if data==None or (isinstance(data,EMData) and data.get_zsize()<=1) :
			print "Error, tried to set data that is invalid for EMIsosurface"
			return
		self.data=data
		self.isorender=MarchingCubes(data)
		self.get_iso_dl()
		self.updateGL()
	
	def set_data(self,data):
		"""Pass in a 3D EMData object"""
		
		if data==None:
			print "Error, tried to set data that is invalid for EMIsosurface"
			return
		self.data=data
		if self.isodl != 0:
			self.force_update = True
		
		self.minden=data.get_attr("minimum")
		self.maxden=data.get_attr("maximum")
		mean=data.get_attr("mean")
		sigma=data.get_attr("sigma")
		
		
		#d = data.get_attr_dict()
		#x,y,z,act = 0,0,0,False
		#if d.has_key("origin_x"):
			
			#x =  d["origin_x"]/d["apix_x"] + data.get_xsize()/2
			#act = True
		#if d.has_key("origin_y"):
			#y =  d["origin_y"]/d["apix_y"] + data.get_ysize()/2
			#act = True
		#if d.has_key("origin_z"):
			#z =  d["origin_z"]/d["apix_z"] + data.get_zsize()/2
			#act = True

		if not self.inspector or self.inspector == None:
			self.inspector=EMIsoInspector(self)
		
		#if act:
			#self.inspector.set_xyz_trans(x,y,z)
		
		hist = data.calc_hist(256,self.minden,self.maxden)
		self.inspector.set_hist(hist,self.minden,self.maxden) 
		iso_threshold = mean+3.0*sigma
		self.inspector.set_thresholds(self.minden,self.maxden,iso_threshold)
		self.isothr = iso_threshold
		self.brightness = -self.isothr
		
		self.isorender=MarchingCubes(data)
		self.inspector.set_sampling_range(self.isorender.get_sampling_range())
#		nx,ny,nz = data.get_xsize(),data.get_ysize(),data.get_zsize()
#		if nx > 256 or ny > 256 or nz > 256:
#			self.isorender.set_sampling(2)
		
		self.load_colors()
		self.inspector.set_materials(self.colors,self.isocolor)
		
		from emglobjects import EM3DGLWidget
		if isinstance(self.get_gl_widget(),EM3DGLWidget):
			self.get_gl_widget().set_camera_defaults(self.data)
	
	def load_colors(self):
		self.colors = get_default_gl_colors()
		
		self.isocolor = "bluewhite"
	
	def get_material(self):
		return self.colors[self.isocolor]
	
	def set_threshold(self,val):
		if (self.isothr != val):
			self.isothr = val
			self.brightness = -val
			if ( self.texture ):
				self.update_data_and_texture()
			self.get_iso_dl()
		
			if self.emit_events: self.emit(QtCore.SIGNAL("set_threshold"),val)
			self.updateGL()
	
	def set_sample(self,val):
		if ( self.smpval != int(val)):
			# the minus two is here because the marching cubes thinks -1 is the high level of detail, 0 is the next best and  so forth
			# However the user wants the highest level of detail to be 1, and the next best to be 2 and then 3 etc
			self.smpval = int(val)-2
			self.get_iso_dl()
			self.updateGL()
	
	def set_material(self,val):
		#print val
		self.isocolor = str(val)
		self.updateGL()
		
	def toggle_cube(self):
		self.cube = not self.cube
		self.updateGL()
	
	def toggle_wire(self,val):
		self.wire = not self.wire
		self.updateGL()
		
	def toggle_light(self,val):
		self.light = not self.light
		self.updateGL()
	
	def toggle_texture(self):
		self.texture = not self.texture
		if ( self.texture ):
			self.update_data_and_texture()
		
		self.get_iso_dl()
		self.updateGL()
	
	def update_inspector(self,t3d):
		self.get_inspector().update_rotations(t3d)
	
	def get_inspector(self):
		if not self.inspector : self.inspector=EMIsoInspector(self,self.enable_file_browse)
		return self.inspector
		
	def set_contrast(self,val):
		self.contrast = val
		self.update_data_and_texture()
		self.updateGL()
		
	def set_brightness(self,val):
		self.brightness = val
		self.update_data_and_texture()
		self.updateGL()
		
	def resize(self):
		self.vdtools.set_update_P_inv()

	def get_mrc_file(self): #added by muthu
		return self.get_inspector().mrcfileName
		

class EMIsoInspector(QtGui.QWidget):
	def __init__(self,target,enable_browse=False) :
		QtGui.QWidget.__init__(self,None)

		self.setWindowIcon(QtGui.QIcon(get_image_directory() +"desktop.png"))
		self.target=weakref.ref(target)
		self.rotation_sliders = EMTransformPanel(target,self)
		
		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		
		self.mrcChanged = False #added by Muthu
		
		if enable_browse:
			hblbrowse = QtGui.QHBoxLayout()
			self.mrc_text = QtGui.QLineEdit()
			hblbrowse.addWidget(self.mrc_text)
			self.mrc_browse = QtGui.QPushButton("Browse")
			hblbrowse.addWidget(self.mrc_browse)
			self.vbl.addLayout(hblbrowse)

			QtCore.QObject.connect(self.mrc_text, QtCore.SIGNAL("textEdited(const QString&)"), self.on_mrc_text_change) #added by Muthu
			QtCore.QObject.connect(self.mrc_browse, QtCore.SIGNAL("clicked(bool)"), self.on_mrc_browse) # added by Muthu

		self.hbl = QtGui.QHBoxLayout()
		self.hbl.setMargin(0)
		self.hbl.setSpacing(6)
		self.hbl.setObjectName("hbl")
		self.vbl.addLayout(self.hbl)
		
		self.hist = ImgHistogram(self)
		self.hist.setObjectName("hist")
		self.hbl.addWidget(self.hist)
		
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
		self.vbl2.addWidget(self.lighttog)
		
		self.cubetog = QtGui.QPushButton("Cube")
		self.cubetog.setCheckable(1)
		self.vbl2.addWidget(self.cubetog)
		
		self.texturetog = QtGui.QPushButton("Texture")
		self.texturetog.setCheckable(1)
		self.vbl2.addWidget(self.texturetog)
		self.texture = False
		
		self.tabwidget = QtGui.QTabWidget()
		self.maintab = None
		self.tabwidget.addTab(self.get_main_tab(), "Main")
		self.texturetab = None
		self.tabwidget.addTab(self.get_GL_tab(),"GL")
		self.tabwidget.addTab(self.get_texture_tab(),"Texture")
		self.get_texture_tab().setEnabled(False)
		self.vbl.addWidget(self.tabwidget)
		self.n3_showing = False
		
		QtCore.QObject.connect(self.thr, QtCore.SIGNAL("valueChanged"), self.on_threshold_slider)
		QtCore.QObject.connect(self.contrast, QtCore.SIGNAL("valueChanged"), target.set_contrast)
		QtCore.QObject.connect(self.bright, QtCore.SIGNAL("valueChanged"), target.set_brightness)
		QtCore.QObject.connect(self.cbb, QtCore.SIGNAL("currentIndexChanged(QString)"), self.set_material)
		QtCore.QObject.connect(self.smp, QtCore.SIGNAL("valueChanged(int)"), target.set_sample)
		QtCore.QObject.connect(self.wiretog, QtCore.SIGNAL("toggled(bool)"), target.toggle_wire)
		QtCore.QObject.connect(self.lighttog, QtCore.SIGNAL("toggled(bool)"), target.toggle_light)
		QtCore.QObject.connect(self.texturetog, QtCore.SIGNAL("toggled(bool)"), self.toggle_texture)
		QtCore.QObject.connect(self.cubetog, QtCore.SIGNAL("toggled(bool)"), target.toggle_cube)
		QtCore.QObject.connect(self.glcontrast, QtCore.SIGNAL("valueChanged"), target.set_GL_contrast)
		QtCore.QObject.connect(self.glbrightness, QtCore.SIGNAL("valueChanged"), target.set_GL_brightness)
		
		QtCore.QObject.connect(self.ambient_tab.r, QtCore.SIGNAL("valueChanged"), self.update_material)
		QtCore.QObject.connect(self.ambient_tab.g, QtCore.SIGNAL("valueChanged"), self.update_material)
		QtCore.QObject.connect(self.ambient_tab.b, QtCore.SIGNAL("valueChanged"), self.update_material)
		QtCore.QObject.connect(self.diffuse_tab.r, QtCore.SIGNAL("valueChanged"), self.update_material)
		QtCore.QObject.connect(self.diffuse_tab.g, QtCore.SIGNAL("valueChanged"), self.update_material)
		QtCore.QObject.connect(self.diffuse_tab.b, QtCore.SIGNAL("valueChanged"), self.update_material)
		QtCore.QObject.connect(self.specular_tab.r, QtCore.SIGNAL("valueChanged"), self.update_material)
		QtCore.QObject.connect(self.specular_tab.g, QtCore.SIGNAL("valueChanged"), self.update_material)
		QtCore.QObject.connect(self.specular_tab.b, QtCore.SIGNAL("valueChanged"), self.update_material)
		QtCore.QObject.connect(self.emission_tab.r, QtCore.SIGNAL("valueChanged"), self.update_material)
		QtCore.QObject.connect(self.emission_tab.g, QtCore.SIGNAL("valueChanged"), self.update_material)
		QtCore.QObject.connect(self.emission_tab.b, QtCore.SIGNAL("valueChanged"), self.update_material)
		QtCore.QObject.connect(self.shininess, QtCore.SIGNAL("valueChanged"), self.update_material)



	def on_mrc_text_change(self,text): #if enable_browse, added by muthu
		print "Use the Browse button to update the mrc file"

	def on_mrc_browse(self): #if enable_browse, added by muthu
		import os
		self.mrcfileName = QtGui.QFileDialog.getOpenFileName(self, "open file", os.getcwd(), "Text files (*.mrc)")
		if (self.mrcfileName == ""): return
		mrcData = EMData(str(self.mrcfileName))
		self.target().set_data(mrcData)
		self.mrc_text.setText(self.mrcfileName) 
		self.mrcChanged = True
		self.target().updateGL()
	
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
	
	
	def update_material(self):
		self.target().isocolor = "custom"
		custom = {}
		
		custom["ambient"] = [self.ambient_tab.r.getValue(), self.ambient_tab.g.getValue(), self.ambient_tab.b.getValue(),1.0]
		custom["diffuse"] = [self.diffuse_tab.r.getValue(), self.diffuse_tab.g.getValue(), self.diffuse_tab.b.getValue(),1.0]
		custom["specular"] = [self.specular_tab.r.getValue(), self.specular_tab.g.getValue(), self.specular_tab.b.getValue(),1.0]
		custom["emission"] = [self.emission_tab.r.getValue(), self.emission_tab.g.getValue(), self.emission_tab.b.getValue(),1.0]
		custom["shininess"] = self.shininess.getValue()
		self.target().colors["custom"] = custom

		n = self.cbb.findText(QtCore.QString("custom"))
		if n < 0: return
		self.cbb.setCurrentIndex(n)
		self.target().updateGL()
	
	def set_material(self,color):
		self.target().set_material(color)
		material = self.target().get_material()
		
		self.ambient_tab.r.setValue(material["ambient"][0])
		self.ambient_tab.g.setValue(material["ambient"][1])
		self.ambient_tab.b.setValue(material["ambient"][2])
		
		self.diffuse_tab.r.setValue(material["diffuse"][0])
		self.diffuse_tab.g.setValue(material["diffuse"][1])
		self.diffuse_tab.b.setValue(material["diffuse"][2])
		
		self.specular_tab.r.setValue(material["specular"][0])
		self.specular_tab.g.setValue(material["specular"][1])
		self.specular_tab.b.setValue(material["specular"][2])
		
		self.emission_tab.r.setValue(material["emission"][0])
		self.emission_tab.g.setValue(material["emission"][1])
		self.emission_tab.b.setValue(material["emission"][2])
		
		self.shininess.setValue(material["shininess"])
	
	def get_RGB_tab(self, name=""):
		return get_RGB_tab(self,name)
		#rgbtab = QtGui.QWidget(self)
		#rgbtab.vbl = QtGui.QVBoxLayout(rgbtab)
		#rgbtab.vbl.setMargin(0)
		#rgbtab.vbl.setSpacing(6)
		#rgbtab.vbl.setObjectName(name)
		
		#rgbtab.r = ValSlider(rgbtab,(0.0,1.0),"R:")
		#rgbtab.r.setObjectName("R")
		#rgbtab.r.setValue(0.5)
		#rgbtab.vbl.addWidget(rgbtab.r)
		
		#rgbtab.g = ValSlider(rgbtab,(0.0,1.0),"G:")
		#rgbtab.g.setObjectName("G")
		#rgbtab.g.setValue(0.5)
		#rgbtab.vbl.addWidget(rgbtab.g)
		
		#rgbtab.b = ValSlider(rgbtab,(0.0,1.0),"B:")
		#rgbtab.b.setObjectName("B")
		#rgbtab.b.setValue(0.5)
		#rgbtab.vbl.addWidget(rgbtab.b)
		
		#return rgbtab
	
	def get_GL_tab(self):
		self.gltab = QtGui.QWidget()
		gltab = self.gltab
		
		gltab.vbl = QtGui.QVBoxLayout(self.gltab )
		gltab.vbl.setMargin(0)
		gltab.vbl.setSpacing(6)
		gltab.vbl.setObjectName("GL")
		
		
		self.glcontrast = ValSlider(gltab,(1.0,5.0),"GLShd:")
		self.glcontrast.setObjectName("GLShade")
		self.glcontrast.setValue(1.0)
		gltab.vbl.addWidget(self.glcontrast)
		
		self.glbrightness = ValSlider(gltab,(-1.0,0.0),"GLBst:")
		self.glbrightness.setObjectName("GLBoost")
		self.glbrightness.setValue(0.1)
		self.glbrightness.setValue(0.0)
		gltab.vbl.addWidget(self.glbrightness)
	
		self.material_tab_widget = QtGui.QTabWidget()
		self.ambient_tab = self.get_RGB_tab("ambient")
		self.material_tab_widget.addTab(self.ambient_tab, "Ambient")
		
		self.diffuse_tab = self.get_RGB_tab("diffuse")
		self.material_tab_widget.addTab(self.diffuse_tab, "Diffuse")
		
		self.specular_tab = self.get_RGB_tab("specular")
		self.material_tab_widget.addTab(self.specular_tab, "Specular")
		
		self.emission_tab = self.get_RGB_tab("emission")
		self.material_tab_widget.addTab(self.emission_tab, "Emission")
		
		gltab.vbl.addWidget(self.material_tab_widget)

		self.shininess = ValSlider(gltab,(0,128),"Shininess:")
		self.shininess.setObjectName("Shininess")
		self.shininess.setValue(64)
		gltab.vbl.addWidget(self.shininess)

		self.hbl_color = QtGui.QHBoxLayout()
		self.hbl_color.setMargin(0)
		self.hbl_color.setSpacing(6)
		self.hbl_color.setObjectName("Material")
		gltab.vbl.addLayout(self.hbl_color)
		
		self.color_label = QtGui.QLabel()
		self.color_label.setText('Material')
		self.hbl_color.addWidget(self.color_label)
		
		self.cbb = QtGui.QComboBox(gltab)
		self.hbl_color.addWidget(self.cbb)
		
		return gltab
	
	def toggle_texture(self):
		self.texture = not self.texture
		self.target().toggle_texture()
		self.get_texture_tab().setEnabled(self.texture)
	
	def get_texture_tab(self):
		if ( self.texturetab == None ):
			self.texturetab = QtGui.QWidget()
			texturetab = self.texturetab
			texturetab.vbl = QtGui.QVBoxLayout(self.texturetab)
			texturetab.vbl.setMargin(0)
			texturetab.vbl.setSpacing(6)
			texturetab.vbl.setObjectName("Main")
		
			self.contrast = ValSlider(texturetab,(0.0,20.0),"Cont:")
			self.contrast.setObjectName("contrast")
			self.contrast.setValue(10.0)
			texturetab.vbl.addWidget(self.contrast)
	
			self.bright = ValSlider(texturetab,(-5.0,5.0),"Brt:")
			self.bright.setObjectName("bright")
			self.bright.setValue(0.1)
			self.bright.setValue(0.0)
			texturetab.vbl.addWidget(self.bright)
			
			#self.glcontrast = ValSlider(texturetab,(1.0,5.0),"GLShd:")
			#self.glcontrast.setObjectName("GLShade")
			#self.glcontrast.setValue(1.0)
			#texturetab.vbl.addWidget(self.glcontrast)
			
			#self.glbrightness = ValSlider(texturetab,(-1.0,0.0),"GLBst:")
			#self.glbrightness.setObjectName("GLBoost")
			#self.glbrightness.setValue(0.1)
			#self.glbrightness.setValue(0.0)
			#texturetab.vbl.addWidget(self.glbrightness)
			
		return self.texturetab
	
	def get_main_tab(self):
		if ( self.maintab == None ):
			self.maintab = QtGui.QWidget()
			maintab = self.maintab
			maintab.vbl = QtGui.QVBoxLayout(self.maintab)
			maintab.vbl.setMargin(0)
			maintab.vbl.setSpacing(6)
			maintab.vbl.setObjectName("Main")
			
			self.thr = ValSlider(maintab,(0.0,4.0),"Thr:")
			self.thr.setObjectName("thr")
			self.thr.setValue(0.5)
			maintab.vbl.addWidget(self.thr)
			
			self.hbl_smp = QtGui.QHBoxLayout()
			self.hbl_smp.setMargin(0)
			self.hbl_smp.setSpacing(6)
			self.hbl_smp.setObjectName("Sample")
			maintab.vbl.addLayout(self.hbl_smp)
			
			self.smp_label = QtGui.QLabel()
			self.smp_label.setText('Sample Level')
			self.hbl_smp.addWidget(self.smp_label)
			
			self.smp = QtGui.QSpinBox(maintab)
			self.smp.setValue(1)
			self.hbl_smp.addWidget(self.smp)
	
			self.rotation_sliders.addWidgets(maintab.vbl)
			
		return self.maintab
	
	def set_sampling_range(self,range):
		self.smp.setMinimum(1)
		self.smp.setMaximum(1+range-1)
	
	def slider_rotate(self):
		self.target().load_rotation(self.get_current_rotation())
	
	
	def set_materials(self,colors,current_color):
		a = 0
		for i in colors:
			self.cbb.addItem(i)
			if ( i == current_color):
				self.cbb.setCurrentIndex(a)
			a += 1

	def on_threshold_slider(self,val):
		self.target().set_threshold(val)
		self.bright.setValue(-val,True)
		
	def set_thresholds(self,low,high,val):
		self.thr.setRange(low,high)
		self.thr.setValue(val, True)
		self.bright.setValue(-val,True)
	
	def set_sample(self,low,high,val):
		self.smp.setRange(int(low),int(high))
		self.smp.setValue(val, True)
		
	def set_hist(self,hist,minden,maxden):
		self.hist.set_data(hist,minden,maxden)

if __name__ == '__main__':
	from emglobjects import EM3DGLWidget
	from emapplication import EMApp
	app = EMApp()
	window = EM3DGLWidget()
	iso_model = EMIsosurfaceModel(window, test_image_3d(1,size=(64,64,64)))
	window.set_model(iso_model)
	window.updateGL()
	
	app.show()
	app.execute()