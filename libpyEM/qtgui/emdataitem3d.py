#!/usr/bin/env python
#
# Author: Ross Coleman (racolema@bcm.edu)
# Copyright (c) 2000-2011 Baylor College of Medicine

from OpenGL import GL
from OpenGL.GL import *
from PyQt4 import QtCore, QtGui
from libpyGLUtils2 import GLUtil
from EMAN2 import EMData, MarchingCubes
from emitem3d import EMItem3D, EMItem3DInspector
from emscene3d import EMInspectorControlShape
from valslider import ValSlider

from emglobjects import get_default_gl_colors


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
	
	def addControls(self, vboxlayout):
		hblbrowse = QtGui.QHBoxLayout()
		self.file_line_edit = QtGui.QLineEdit()
		hblbrowse.addWidget(self.file_line_edit)
		self.file_browse_button = QtGui.QPushButton("Browse")
		hblbrowse.addWidget(self.file_browse_button)
		vboxlayout.addLayout(hblbrowse)
		
		self.file_browse_button.clicked.connect(self.on_file_browse)
		
	def on_file_browse(self):
		#TODO: replace this with an EMAN2 browser window once we re-write it
		file_path = QtGui.QFileDialog.getOpenFileName(self, "Open 3D Volume Map")
		self.file_line_edit.setText(file_path)
		new_emdata = EMData(file_path)
		self.item3d.data = new_emdata
		for child in self.item3d.getChildren():
			child.inspector.dataChanged()
		self.inspector.updateSceneGraph()

class EMIsosurfaceInspector(EMInspectorControlShape):
	def __init__(self, name, item3d):
		EMItem3DInspector.__init__(self, name, item3d)
		
		QtCore.QObject.connect(self.thr, QtCore.SIGNAL("valueChanged"), self.onThresholdSlider)
		QtCore.QObject.connect(self.cullbackface, QtCore.SIGNAL("stateChanged(int)"), self.onCullFaces)
		self.dataChanged()
		
	def addControls(self, vbox):
		# Perhaps we sould allow the inspector control this?
		self.cullbackface = QtGui.QCheckBox("Cull Back Face Polygons")
		self.cullbackface.setCheckState(QtCore.Qt.Checked) 
		self.thr = ValSlider(self,(0.0,4.0),"Thr:")
		self.thr.setObjectName("thr")
		self.thr.setValue(0.5)
		
		vbox.addWidget(self.cullbackface)
		vbox.addWidget(self.thr)
	
	def dataChanged(self):
		data = self.item3d.parent.data
		
		minden=data.get_attr("minimum")
		maxden=data.get_attr("maximum")
		mean=data.get_attr("mean")
		sigma=data.get_attr("sigma")
		iso_threshold = mean+3.0*sigma
		
		self.thr.setRange(minden,maxden)
		self.thr.setValue(iso_threshold, True)
		
	def onThresholdSlider(self,val):
		self.item3d.setThreshold(val)
#		self.bright.setValue(-val,True)
		self.inspector.updateSceneGraph()
		
	def onCullFaces(self):
		if self.cullbackface.checkState() == QtCore.Qt.Checked:
			self.item3d.cullbackfaces = True
		else:
			self.item3d.cullbackfaces = False
		self.inspector.updateSceneGraph()

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
		self.cullbackfaces = True
		
		self.tex_name = 0
		self.texture = False

		self.brightness = 0
		self.contrast = 10
		self.rank = 1
		self.data_copy = None		
		#self.vdtools = EMViewportDepthTools(self)
		#self.enable_file_browse = enable_file_browse
		self.force_update = False
		
		self.load_colors()
		data = self.parent.data
		assert isinstance(data, EMData)
		self.isorender = MarchingCubes(data)
		
		# color Needed for inspector to work John Flanagan
		self.diffuse = self.colors[self.isocolor]["diffuse"]
		self.specular = self.colors[self.isocolor]["specular"]
		self.ambient = self.colors[self.isocolor]["ambient"]		
		self.shininess = self.colors[self.isocolor]["shininess"]
		
	# I have added these methods so the inspector can set the color John Flanagan
	def setAmbientColor(self, red, green, blue, alpha=1.0):
		self.ambient = [red, green, blue, alpha]

	def setDiffuseColor(self, red, green, blue, alpha=1.0):
		self.diffuse = [red, green, blue, alpha]
		
	def setSpecularColor(self, red, green, blue, alpha=1.0):
		self.specular = [red, green, blue, alpha]
	
	def setShininess(self, shininess):
		self.shininess = shininess
		
	def getSceneGui(self):
		if not self.widget:
			self.widget = EMIsosurfaceInspector("ISOSURFACE", self)
		return self.widget
	
	def load_colors(self):
		self.colors = get_default_gl_colors()
		
		self.isocolor = "bluewhite"

	def get_iso_dl(self):
		# create the isosurface display list
		self.isorender.set_surface_value(self.isothr)
		self.isorender.set_sampling(self.smpval)
		
		if ( self.texture ):
			if ( self.tex_name == 0 ):
				self.update_data_and_texture()
		
		face_z = False
		if self.parent.data.get_zsize() <= 2:
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
		#lighting = glIsEnabled(GL_LIGHTING)
		cull = glIsEnabled(GL_CULL_FACE)
		#depth = glIsEnabled(GL_DEPTH_TEST)
		polygonmode = glGetIntegerv(GL_POLYGON_MODE)
		#normalize = glIsEnabled(GL_NORMALIZE)
	
		if self.cullbackfaces:
			glEnable(GL_CULL_FACE)
			glCullFace(GL_BACK)
		else:
			if not cull:
				glDisable(GL_CULL_FACE)
		# The lighting, depthtest, and normalization are controlled at the EMScene3d level. It should not be down with child widgets John Flanagan
		#glDisable(GL_CULL_FACE)
		#glEnable(GL_DEPTH_TEST)
		#glEnable(GL_NORMALIZE)
		#glDisable(GL_NORMALIZE)
		if ( self.wire ):
			glPolygonMode(GL_FRONT,GL_LINE);
		else:
			glPolygonMode(GL_FRONT,GL_FILL);
		
		#if self.light:
		#	glEnable(GL_LIGHTING)
		#else:
		#	glDisable(GL_LIGHTING)

		
#		glPushMatrix()
#		self.cam.position(True)
#		# the ones are dummy variables atm... they don't do anything
#		self.vdtools.update(1,1)
#		glPopMatrix()
		
#		self.cam.position()
		glShadeModel(GL_SMOOTH)
		if ( self.isodl == 0 or self.force_update):
			self.get_iso_dl()
			self.force_update = False
		glStencilFunc(GL_EQUAL,self.rank,0)
		glStencilOp(GL_KEEP,GL_KEEP,GL_REPLACE)
		# This is needed for the inspector to work John Flanagan	
		glMaterialfv(GL_FRONT, GL_DIFFUSE, self.diffuse)
		glMaterialfv(GL_FRONT, GL_SPECULAR, self.specular)
		glMaterialf(GL_FRONT, GL_SHININESS, self.shininess)
		glMaterialfv(GL_FRONT, GL_AMBIENT, self.ambient)
		glColor(self.ambient)
		# I have commented this out b/c it need to be controled by the inspector John Flanagan
		#glMaterial(GL_FRONT, GL_AMBIENT, self.colors[self.isocolor]["ambient"])
		#glMaterial(GL_FRONT, GL_DIFFUSE, self.colors[self.isocolor]["diffuse"])
		#glMaterial(GL_FRONT, GL_SPECULAR, self.colors[self.isocolor]["specular"])
		#glMaterial(GL_FRONT, GL_SHININESS, self.colors[self.isocolor]["shininess"])
		#glMaterial(GL_FRONT, GL_EMISSION, self.colors[self.isocolor]["emission"])
		#glColor(self.colors[self.isocolor]["ambient"])
		
		glPushMatrix()
		glTranslate(-self.parent.data.get_xsize()/2.0,-self.parent.data.get_ysize()/2.0,-self.parent.data.get_zsize()/2.0)
		if ( self.texture ):
			glScalef(self.parent.data.get_xsize(),self.parent.data.get_ysize(),self.parent.data.get_zsize())
		glCallList(self.isodl)
		glPopMatrix()
		
#		self.draw_bc_screen() #TODO: check into porting this from EM3DModel
		
		glStencilFunc(GL_ALWAYS,1,1)
		if self.cube:
			#glDisable(GL_LIGHTING)
			glPushMatrix()
			self.draw_volume_bounds()
			glPopMatrix()
		
		# Again the lighting, depthtest, and normalization are controlled at the EMScene3d level. It should not be down with child widgets John Flanagan
		#if ( lighting ): glEnable(GL_LIGHTING)
		#else: glDisable(GL_LIGHTING)
		# What is the point of this conditional testing.... It's always TRUE!!!
		if ( not cull ): glDisable(GL_CULL_FACE)
		else: glDisable(GL_CULL_FACE)
		#if ( depth ): glEnable(GL_DEPTH_TEST)
		#else : glDisable(GL_DEPTH_TEST)
		
		#if ( not normalize ): glDisable(GL_NORMALIZE)
		
		# What is the point of this conditonal testing.... It's always TRUE!!!
		if ( polygonmode[0] == GL_LINE ): glPolygonMode(GL_FRONT, GL_LINE)
		else: glPolygonMode(GL_FRONT, GL_FILL)
		#if ( polygonmode[1] == GL_LINE ): glPolygonMode(GL_BACK, GL_LINE)
		#else: glPolygonMode(GL_BACK, GL_FILL)
		
		#print "total time is", time()-a
		
	def setThreshold(self, val):
		if (self.isothr != val):
			self.isothr = val
#			self.brightness = -val
#			if ( self.texture ):
#				self.update_data_and_texture()
			self.get_iso_dl()
		
#			if self.emit_events: self.emit(QtCore.SIGNAL("set_threshold"),val)
#			self.updateGL()
		