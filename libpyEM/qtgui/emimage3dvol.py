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
import weakref
from time import *

from emglobjects import EM3DModel, EMViewportDepthTools, Camera2, Camera, EMOpenGLFlagsAndTools
from emimageutil import ImgHistogram, EMTransformPanel

MAG_INCREMENT_FACTOR = 1.1


class EMVolumeModel(EM3DModel):
	def __init__(self, gl_widget, image=None):
		self.data = None
		EM3DModel.__init__(self, gl_widget)
		
		self.init()
		self.initialized = True
		
		self.initializedGL= False
		
		self.inspector=None
		
		self.tex_names_list = []		# A storage object, used to remember and later delete texture names
		
		self.axes_idx = -1
		self.axes = []
		self.axes.append(Vec3f(1,0,0))
		self.axes.append(Vec3f(0,1,0))
		self.axes.append(Vec3f(0,0,1))
		#self.axes.append( Vec3f(-1,0,0) )
		#self.axes.append( Vec3f(0,-1,0) )
		#self.add_render_axis(1,1,1)
		#self.add_render_axis(-1,1,1)
		#self.add_render_axis(-1,-1,1)
		#self.add_render_axis(1,-1,1)
		
		#self.add_render_axis(1,1,0)
		#self.add_render_axis(-1,1,0)
		#self.add_render_axis(-1,-1,0)
		#self.add_render_axis(1,-1,0)
		
		#self.add_render_axis(0,1,1)
		#self.add_render_axis(0,-1,1)
		
		#self.add_render_axis(1,0,1)
		#self.add_render_axis(-1,0,1)
			
		if image :
			self.set_data(image)
#	def __del__(self):
#		print "vol died"

	def add_render_axis(self,a,b,c):
		v = Vec3f(a,b,c);
		v.normalize()
		self.axes.append( v )
	
	def eye_coords_dif(self,x1,y1,x2,y2,mdepth=True):
		return self.vdtools.eye_coords_dif(x1,y1,x2,y2,mdepth)
	
	def get_type(self):
		return "Volume"

	def init(self):
		self.data=None

		self.mmode=0
		
		self.cam = Camera2(self)
		self.cube = False
		
		self.vdtools = EMViewportDepthTools(self)
		
		self.contrast = 1.0
		self.brightness = 0.0
		self.texsample = 1.0
		self.glcontrast = 1.0
		self.glbrightness = 0.0
		self.cube = False
		
		self.tex_name = 0
		
		self.rank = 1
		
		self.tex_dl = 0
		self.inspector=None
		
		self.force_texture_update = False

		self.glflags = EMOpenGLFlagsAndTools()		# OpenGL flags - this is a singleton convenience class for testing texture support
	
	def update_data(self,data):
		self.set_data(data)
		self.updateGL()
		
	def set_data(self,data):
		"""Pass in a 3D EMData object"""
		
		
		self.data=data
		if data==None:
			print "Error, the data is empty"
			return
		
	 	if (isinstance(data,EMData) and data.get_zsize()<=1) :
			print "Error, the data is not 3D"
			return
		
		if not self.inspector or self.inspector ==None:
			self.inspector=EMVolumeInspector(self)
		
		self.update_data_and_texture()
		
		from emglobjects import EM3DGLWidget
		if isinstance(self.get_gl_widget(),EM3DGLWidget):
			self.get_gl_widget().set_camera_defaults(self.data)
		
	def test_accum(self):
		# this code will do volume rendering using the accumulation buffer
		# I opted not to go this way because you can't retain depth in the accumulation buffer
		# Note that it only works in the z-direction
		glClear(GL_ACCUM_BUFFER_BIT)
		
		self.accum = True
		self.zsample = self.texsample*(self.data.get_zsize())
		
		if self.tex_name == 0:
			print "Error, can not render 3D texture - texture name is 0"
			return
		
		
		for z in range(0,int(self.texsample*(self.data.get_zsize()))):
			glEnable(GL_TEXTURE_3D)
			glBindTexture(GL_TEXTURE_3D, self.tex_name)
			glPixelStorei(GL_UNPACK_ALIGNMENT, 1)
			glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP)
			glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP)
			glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP)
			glTexParameter(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
			glTexParameter(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)
			glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE)

			glBegin(GL_QUADS)
			
			zz = float(z)/float(self.data.get_zsize()-1)/self.texsample
			glTexCoord3f(0,0,zz)
			glVertex3f(0,0,zz)
			
			glTexCoord3f(1,0,zz)
			glVertex3f(1,0,zz)
			
			glTexCoord3f(1,1,zz)
			glVertex3f(1,1,zz)
			
			glTexCoord3f(0,1,zz)
			glVertex3f(0,1,zz)
		
			glEnd()
			glDisable(GL_TEXTURE_3D)
		
			if ( self.accum ):
				glAccum(GL_ADD, 1.0/self.zsample*self.brightness)
				glAccum(GL_ACCUM, 1.0/self.zsample*self.contrast)
		
		
		glAccum(GL_RETURN, 1.0)
		
	def render(self):
		lighting = glIsEnabled(GL_LIGHTING)
		cull = glIsEnabled(GL_CULL_FACE)
		depth = glIsEnabled(GL_DEPTH_TEST)
		
		polygonmode = glGetIntegerv(GL_POLYGON_MODE)

		glDisable(GL_LIGHTING)
		glDisable(GL_CULL_FACE)
		glDisable(GL_DEPTH_TEST)
		
		glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
		
		glPushMatrix()
		self.cam.position(True)
		# the ones are dummy variables atm... they don't do anything
		self.vdtools.update(1,1)
		glPopMatrix()
		
		self.cam.position()
		self.vdtools.store_model()

		# here is where the correct display list (x,y or z direction) is determined
		self.texture_update_if_necessary()

		glStencilFunc(GL_EQUAL,self.rank,0)
		glStencilOp(GL_KEEP,GL_KEEP,GL_REPLACE)
		glPushMatrix()
		glTranslate(-self.data.get_xsize()/2.0,-self.data.get_ysize()/2.0,-self.data.get_zsize()/2.0)
		glScalef(self.data.get_xsize(),self.data.get_ysize(),self.data.get_zsize())
		glEnable(GL_BLEND)
		#glBlendEquation(GL_MAX)
		if self.glflags.blend_equation_supported():
			glBlendEquation(GL_FUNC_ADD)
		glDepthMask(GL_FALSE)
		glBlendFunc(GL_ONE, GL_ONE)
		glCallList(self.tex_dl)
		glDepthMask(GL_TRUE)
		glDisable(GL_BLEND)
		glPopMatrix()

		# this is the accumulation buffer version of the volume model - it was for testing purposes
		# and is left here commented out incase anyone wants to investigate it in the future
		#glPushMatrix()
		#glTranslate(-self.data.get_xsize()/2.0,-self.data.get_ysize()/2.0,-self.data.get_zsize()/2.0)
		#glScalef(self.data.get_xsize(),self.data.get_ysize(),self.data.get_zsize())
		#self.test_accum()
		#glPopMatrix()
		
		#breaks in desktop!
		#glStencilFunc(GL_EQUAL,self.rank,self.rank)
		#glStencilOp(GL_KEEP,GL_KEEP,GL_KEEP)
		#glPushMatrix()
		#glLoadIdentity()
		#[width,height] = self.parent.get_near_plane_dims()
		#z = self.parent.get_start_z()
		#glTranslate(-width/2.0,-height/2.0,-z-0.01)
		#glScalef(width,height,1.0)
		self.draw_bc_screen()
		#glPopMatrix()
		
		glStencilFunc(GL_ALWAYS,1,1)
		if self.cube:
			glPushMatrix()
			self.draw_volume_bounds()
			glPopMatrix()
			
		if ( lighting ): glEnable(GL_LIGHTING)
		if ( cull ): glEnable(GL_CULL_FACE)
		if ( depth ): glEnable(GL_DEPTH_TEST)
		
		if ( polygonmode[0] == GL_LINE ): glPolygonMode(GL_FRONT, GL_LINE)
		if ( polygonmode[1] == GL_LINE ): glPolygonMode(GL_BACK, GL_LINE)
	
	def texture_update_if_necessary(self):
		
		t3d = self.vdtools.getEmanMatrix()
		
		
		point = Vec3f(0,0,1)
		
		point = point*t3d
		
		point[0] = abs(point[0])
		point[1] = abs(point[1])
		point[2] = abs(point[2])
		#point[1] = -point[1]
		#if ( point[2] < 0 ):
			#point[2] = -point[2]
			#point[1] = -point[1]
			#point[0] = -point[0]
	
		currentaxis = self.axes_idx
		
		closest = 2*pi
		lp = point.length()
		point.normalize()
		idx = 0
		for i in self.axes:
			try:
				angle = abs(acos(point.dot(i)))
			except:
				t3d.printme()
				print 'warning, there is a bug in the volume render which may cause incorrect rendering'
				return
			if (angle < closest):
				closest = angle
				self.axes_idx = idx
			
			idx += 1

		if (currentaxis != self.axes_idx or self.force_texture_update):
			#print self.axes[self.axes_idx]
			self.gen_texture()
			
	def gen_texture(self):
		if self.glflags.threed_texturing_supported():
			self.get_3D_texture()
		else:
			self.gen_2D_texture()
			
	def get_3D_texture(self):
		if ( self.tex_dl != 0 ): glDeleteLists( self.tex_dl, 1)
		
		self.tex_dl = glGenLists(1)
		
		if self.tex_dl == 0:
			print "Error, failed to generate display list"
			return
		
		n = self.get_dimension_size()
		v = self.axes[self.axes_idx]
		[t,alt,phi] = self.get_eman_transform(v)
		v1 = t*Vec3f(-0.5,-0.5,0)
		v2 = t*Vec3f(-0.5, 0.5,0)
		v3 = t*Vec3f( 0.5, 0.5,0)
		v4 = t*Vec3f( 0.5,-0.5,0)
		vecs = [v1,v2,v3,v4]
		total = self.texsample*n
		self.data.add(self.brightness)
		if self.contrast != 0:
			self.data.mult(self.contrast*1.0/total)
		if ( self.force_texture_update ):
			if self.tex_name != 0:
				glDeleteTextures(self.tex_name)
			
			self.tex_name = self.glflags.gen_textureName(self.data)
			
			self.force_texture_update = False
			
		
		if self.tex_name == 0:
			print "Error, can not render 3D texture - texture name is 0"
			return
		
		
		glNewList(self.tex_dl,GL_COMPILE)
		glEnable(GL_TEXTURE_3D)
		glBindTexture(GL_TEXTURE_3D, self.tex_name)
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1)
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE)
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP)
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP)
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP)
		glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
		if ( not data_dims_power_of(self.data,2) and self.glflags.npt_textures_unsupported()):
			glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR)
		else:
			glTexParameterf(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR)

		glPushMatrix()
		glTranslate(0.5,0.5,0.5)
		glBegin(GL_QUADS)
		
		
		for i in range(0,int(self.texsample*n)):
			nn = float(i)/float(n)/self.texsample

			trans = (nn-0.5)*v
			
			for r in vecs:
			
				w = [r[0] + trans[0], r[1] + trans[1], r[2] + trans[2]]
				t = [w[0]+0.5,w[1]+0.5,w[2]+0.5]
				glTexCoord3fv(t)
				glVertex3fv(w)
			
		glEnd()
		glPopMatrix()
		glDisable(GL_TEXTURE_3D)
		glEndList()
		
		if self.contrast != 0:	self.data.mult(total*1.0/self.contrast)
		self.data.add(-self.brightness)
		
		
	def get_eman_transform(self,p):
		
		if ( p[2] == 0 ):
			alt = 90
		else :
			alt = acos(p[2])*180.0/pi
		
		phi = atan2(p[0],p[1])
		phi *= 180.0/pi
		
		t = Transform({"type":"eman","alt":alt,"phi":phi})
		
		return [t,alt,phi]
			
	def get_dimension_size(self):
		if ( self.axes_idx == 0 ):
			return self.data.get_xsize()
		elif ( self.axes_idx == 1 ):
			return self.data.get_ysize()
		elif ( self.axes_idx == 2 ):
			return self.data.get_zsize()
		else:
			#print "unsupported axis"
			# this is a hack and needs to be fixed eventually
			return self.data.get_xsize()
			#return 0
	def get_correct_dims_2d_emdata(self):
		if ( self.axes_idx == 0 ):
			return EMData(self.data.get_ysize(),self.data.get_zsize())
		elif ( self.axes_idx == 1 ):
			return EMData(self.data.get_xsize(),self.data.get_zsize())
		elif ( self.axes_idx == 2 ):
			return EMData(self.data.get_xsize(),self.data.get_ysize())
		else:
			#print "unsupported axis"
			# this is a hack and needs to be fixed eventually
			return EMData(self.data.get_xsize(),self.data.get_zsize())

	
	def gen_2D_texture(self):
			
		if ( self.tex_dl != 0 ): 
			glDeleteLists( self.tex_dl, 1)
		
		for i in self.tex_names_list:
			glDeleteTextures(i)
			
		self.tex_dl = glGenLists(1)
		if (self.tex_dl == 0 ):
			print "error, could not generate list"
			return

		glNewList(self.tex_dl,GL_COMPILE)
		glEnable(GL_TEXTURE_2D)

		n = self.get_dimension_size()
		v = self.axes[self.axes_idx]
		[t,alt,phi] = self.get_eman_transform(v)
		total = self.texsample*n
		for i in range(0,int(self.texsample*n)):
			nn = float(i)/float(n)/self.texsample
			
			trans = (nn-0.5)*v
			t.set_trans(2.0*int(n/2)*trans)
			
			if False and EMUtil.cuda_available(): # disable for the time being - big textures won't work on CPU
				tmp = self.data.cut_slice_cuda(t)
			else:
				tmp = self.get_correct_dims_2d_emdata() 
				tmp.cut_slice(self.data,t,True)
			tmp.add(self.brightness)
			tmp.mult(self.contrast*1.0/total)
			#tmp.write_image("tmp.img",-1)
			
			# get the texture name, store it, and bind it in OpenGL
			tex_name = self.glflags.gen_textureName(tmp)
			self.tex_names_list.append(tex_name)
			glBindTexture(GL_TEXTURE_2D, tex_name)
			
			self.loat_default_2D_texture_parms()
			
			glPushMatrix()
			
			glTranslate(trans[0]+0.5,trans[1]+0.5,trans[2]+0.5)
			glRotatef(-phi,0,0,1)
			glRotatef(-alt,1,0,0)
			glBegin(GL_QUADS)
			glTexCoord2f(0,0)
			glVertex2f(-0.5,-0.5)
			
			glTexCoord2f(1,0)
			glVertex2f( 0.5,-0.5)
			
			glTexCoord2f(1,1)
			glVertex2f( 0.5, 0.5)
			
			glTexCoord2f(0,1)
			glVertex2f(-0.5, 0.5)
			glEnd()
			glPopMatrix()
		
		glDisable(GL_TEXTURE_2D)
		glEndList()
		
		# this may have been toggled (i.e. if the image contrast or brightness changed)
		if self.force_texture_update == True:
			self.force_texture_update = False
	
	def loat_default_2D_texture_parms(self):
		glPixelStorei(GL_UNPACK_ALIGNMENT, 1)
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP)
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP)
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR)
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE)
		if ( not data_dims_power_of(self.data,2) and self.glflags.npt_textures_unsupported()):
			glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR)
		else:
			glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR) 
		
	def update_data_and_texture(self):
	
		if ( not isinstance(self.data,EMData) ): return
		
		hist = self.data.calc_hist(256,0,1.0)
		self.inspector.set_hist(hist,0,1.0)

		self.force_texture_update = True

	def set_contrast(self,val):
		self.contrast = val
		self.update_data_and_texture()
		self.updateGL()
		
	def set_brightness(self,val):
		self.brightness = val
		self.update_data_and_texture()
		self.updateGL()
		
	def set_texture_sample(self,val):
		if ( val < 0 ) :
			print "Error, cannot handle texture sample less than 0"
			return
		
		self.texsample = val
		self.force_texture_update = True
		self.updateGL()

	def update_inspector(self,t3d):
		if not self.inspector or self.inspector ==None:
			self.inspector=EMVolumeInspector(self)
		self.inspector.update_rotations(t3d)
	
	def get_inspector(self):
		if not self.inspector : self.inspector=EMVolumeInspector(self)
		return self.inspector
		
	def resize(self):
		self.vdtools.set_update_P_inv()


class EMVolumeInspector(QtGui.QWidget):
	def __init__(self,target) :
		QtGui.QWidget.__init__(self,None)
		self.target=weakref.ref(target)
		self.setWindowIcon(QtGui.QIcon(get_image_directory() +"desktop.png"))
		self.rotation_sliders = EMTransformPanel(target,self)
		
		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		
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
	
		self.cubetog = QtGui.QPushButton("Cube")
		self.cubetog.setCheckable(1)
		self.vbl2.addWidget(self.cubetog)
		
		self.defaults = QtGui.QPushButton("Defaults")
		self.vbl2.addWidget(self.defaults)
		
		self.tabwidget = QtGui.QTabWidget()
		
		self.tabwidget.addTab(self.get_main_tab(), "Main")
		self.tabwidget.addTab(self.get_GL_tab(),"GL")
		
		self.vbl.addWidget(self.tabwidget)
		
		self.n3_showing = False
		
		QtCore.QObject.connect(self.contrast, QtCore.SIGNAL("valueChanged"), target.set_contrast)
		QtCore.QObject.connect(self.glcontrast, QtCore.SIGNAL("valueChanged"), target.set_GL_contrast)
		QtCore.QObject.connect(self.glbrightness, QtCore.SIGNAL("valueChanged"), target.set_GL_brightness)
		QtCore.QObject.connect(self.bright, QtCore.SIGNAL("valueChanged"), target.set_brightness)
		QtCore.QObject.connect(self.cubetog, QtCore.SIGNAL("toggled(bool)"), target.toggle_cube)
		QtCore.QObject.connect(self.defaults, QtCore.SIGNAL("clicked(bool)"), self.set_defaults)
		QtCore.QObject.connect(self.smp, QtCore.SIGNAL("valueChanged(int)"), target.set_texture_sample)
	
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
	
	def get_GL_tab(self):
		self.gltab = QtGui.QWidget()
		gltab = self.gltab
		
		gltab.vbl = QtGui.QVBoxLayout(self.gltab )
		gltab.vbl.setMargin(0)
		gltab.vbl.setSpacing(6)
		gltab.vbl.setObjectName("Main")
		
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
	
		self.maintab = QtGui.QWidget()
		maintab = self.maintab
		maintab.vbl = QtGui.QVBoxLayout(self.maintab)
		maintab.vbl.setMargin(0)
		maintab.vbl.setSpacing(6)
		maintab.vbl.setObjectName("Main")
			
		self.contrast = ValSlider(maintab,(0.0,20.0),"Cont:")
		self.contrast.setObjectName("contrast")
		self.contrast.setValue(1.0)
		maintab.vbl.addWidget(self.contrast)

		self.bright = ValSlider(maintab,(-5.0,5.0),"Brt:")
		self.bright.setObjectName("bright")
		self.bright.setValue(0.1)
		self.bright.setValue(0.0)
		maintab.vbl.addWidget(self.bright)

		self.hbl_smp = QtGui.QHBoxLayout()
		self.hbl_smp.setMargin(0)
		self.hbl_smp.setSpacing(6)
		self.hbl_smp.setObjectName("Texture Oversampling")
		maintab.vbl.addLayout(self.hbl_smp)
		
		self.smp_label = QtGui.QLabel()
		self.smp_label.setText('Texture Oversampling')
		self.hbl_smp.addWidget(self.smp_label)
		
		self.smp = QtGui.QSpinBox(maintab)
		self.smp.setMaximum(10)
		self.smp.setMinimum(1)
		self.smp.setValue(1)
		self.hbl_smp.addWidget(self.smp)

		self.lowlim=0
		self.highlim=1.0
		self.busy=0

		self.rotation_sliders.addWidgets(maintab.vbl)

		return maintab
	
	def set_defaults(self):
		self.contrast.setValue(1.0)
		self.bright.setValue(0.0)
		self.glcontrast.setValue(1.0)
		self.glbrightness.setValue(0.0)
		self.rotation_sliders.set_defaults()
	
	def slider_rotate(self):
		self.target().load_rotation(self.get_current_rotation())
	
	def set_hist(self,hist,minden,maxden):
		self.hist.set_data(hist,minden,maxden)

	#def set_scale(self,newscale):
		#self.scale.setValue(newscale)


if __name__ == '__main__':
	from emapplication import EMApp
	from emglobjects import EM3DGLWidget
	
	em_app = EMApp()
	window = EM3DGLWidget()
	vol_model = EMVolumeModel(window)
	window.set_model(vol_model)
	
	if len(sys.argv)==1 : 
		data = []
		#for i in range(0,200):
		e = test_image_3d(1,size=(64,64,64))
		window.set_data(e)
	else :
		a=EMData(sys.argv[1])
		vol_model.set_file_name(sys.argv[1])
		window.set_data(a)

	em_app.show()
	em_app.execute()
