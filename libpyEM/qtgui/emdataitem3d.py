#!/usr/bin/env python
#
# Author: Ross Coleman (racolema@bcm.edu)
# Copyright (c) 2000-2011 Baylor College of Medicine

from OpenGL import GL
from EMAN2 import EMData, MarchingCubes
from emitem3d import EMItem3D, EMItem3DInspector


class EMDataItem3D(EMItem3D):
	def __init__(self, data, parent = None, children = set(), transform = None):
		self.data = data
		EMItem3D.__init__(self, parent, children, transform)
	def getSceneGui(self):
		if not self.widget:
			self.widget = EMDataItem3DInspector("DATA", self)
		return self.widget
	
class EMDataItem3DInspector(EMItem3DInspector):
	def __init__(self, name, item3d):
		EMItem3DInspector.__init__(self, name, item3d)


class EMIsosurface(EMItem3D):
	def __init__(self, parent, children = set(), transform = None):
		EMItem3D.__init__(self, parent, children, transform)
		
		#self.mmode = 0
		#self.inspector=None
		self.isothr=0.5
		self.isorender=None
		self.isodl = 0
		self.smpval=-1
		self.griddl = 0
		self.scale = 1.0
		self.cube = False
		self.wire = False
		self.light = True
		
		self.tex_name = 0
		self.texture = False

		self.brightness = 0
		self.contrast = 10
		self.rank = 1
		self.data_copy = None		
		#self.vdtools = EMViewportDepthTools(self)
		#self.enable_file_browse = enable_file_browse
		self.force_update = False
		
		data = self.parent.data
		assert isinstance(data, EMData)
		self.isorender = MarchingCubes(data)
		
	
	def getSceneGui(self):
		if not self.widget: #TODO: code and use EMIsosurfaceInspector
			self.widget = EMDataItem3DInspector("ISO", self)
		return self.widget

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
		
	def renderNode(self):
		if (not isinstance(self.parent.data,EMData)): return
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
		