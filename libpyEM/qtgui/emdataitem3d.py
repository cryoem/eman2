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
from emimageutil import ImgHistogram
from valslider import ValSlider
from emshapeitem3d import EMInspectorControlShape
from emglobjects import get_default_gl_colors
import os.path

class EMDataItem3D(EMItem3D):
	name = "Data"
	def __init__(self, data, parent = None, children = set(), transform = None):
		EMItem3D.__init__(self, parent, children, transform)
		self.setData(data)
		
	def getEvalString(self):
		return "EMDataItem3D(\"%s\")"%self.path
		
	def getItemInspector(self):
		if not self.item_inspector:
			self.item_inspector = EMDataItem3DInspector("DATA", self)
		return self.item_inspector
	
	def setData(self, data):
		if isinstance(data, EMData):
			self.data = data
			if data.has_attr("source_path"):
				self.path = data["source_path"]
			else:
				self.path = None 
		else:
			self.data = EMData(str(data))
			self.path = str(data)		
	
class EMDataItem3DInspector(EMItem3DInspector):
	def __init__(self, name, item3d):
		EMItem3DInspector.__init__(self, name, item3d)
	
	def updateItemControls(self):
		super(EMDataItem3DInspector, self).updateItemControls()
		# Anything that needs to be updated when the scene is rendered goes here.....
		
	def addControls(self, gridbox):
		super(EMDataItem3DInspector, self).addControls(gridbox)
		hblbrowse = QtGui.QHBoxLayout()
		self.file_line_edit = QtGui.QLineEdit()
		hblbrowse.addWidget(self.file_line_edit)
		self.file_browse_button = QtGui.QPushButton("Browse")
		hblbrowse.addWidget(self.file_browse_button)
		gridbox.addLayout(hblbrowse, 3, 0, 1, 1)
		
		self.file_browse_button.clicked.connect(self.onFileBrowse)
		
		# Set to default, but run only once and not in each base class
		if type(self) == EMDataItem3DInspector: self.updateItemControls()
		
	def onFileBrowse(self):
		#TODO: replace this with an EMAN2 browser window once we re-write it
		file_path = QtGui.QFileDialog.getOpenFileName(self, "Open 3D Volume Map")
		self.file_line_edit.setText(file_path)
		self.item3d().setData(file_path)
		for child in self.item3d().getChildren():
			child.getItemInspector().dataChanged()
		self.inspector.updateSceneGraph()

class EMIsosurfaceInspector(EMInspectorControlShape):
	def __init__(self, name, item3d):
		EMItem3DInspector.__init__(self, name, item3d, numgridcols=2)	# for the iso inspector we need two grid cols for extra space....
		
		QtCore.QObject.connect(self.thr, QtCore.SIGNAL("valueChanged"), self.onThresholdSlider)
		self.cullbackface.toggled.connect(self.onCullFaces)
		self.wireframe.toggled.connect(self.onWireframe)
		self.sampling_spinbox.valueChanged[int].connect(self.onSampling)
		self.dataChanged()
	
	def updateItemControls(self):
		super(EMIsosurfaceInspector, self).updateItemControls()
		# Anything that needs to be updated when the scene is rendered goes here.....
		
	def addControls(self, gridbox):
		super(EMIsosurfaceInspector, self).addControls(gridbox)
		self.histogram_widget = ImgHistogram(self)
		self.histogram_widget.setObjectName("hist")

		# Perhaps we should allow the inspector control this?
		isoframe = QtGui.QFrame()
		isoframe.setFrameShape(QtGui.QFrame.StyledPanel)
		isogridbox = QtGui.QGridLayout()
		
		self.cullbackface = QtGui.QCheckBox("Cull Back Face Polygons")
		self.cullbackface.setChecked(True)
		self.wireframe = QtGui.QCheckBox("Wireframe mode")
		self.wireframe.setChecked(False)
		self.thr = ValSlider(self,(0.0,4.0),"Threshold:")
		self.thr.setObjectName("thr")
		self.thr.setValue(0.5)
		self.sampling_label = QtGui.QLabel("Sample Level:")
		self.sampling_spinbox = QtGui.QSpinBox()
		self.sampling_spinbox.setValue(1)
		sampling_hbox_layout = QtGui.QHBoxLayout()
		sampling_hbox_layout.addWidget(self.sampling_label)
		sampling_hbox_layout.addWidget(self.sampling_spinbox)

		isogridbox.addWidget(self.histogram_widget, 0, 0, 1, 1)
		isogridbox.addWidget(self.cullbackface, 1, 0, 1, 1)
		isogridbox.addWidget(self.wireframe, 2,0,1,1)
		isogridbox.addLayout(sampling_hbox_layout, 3,0,1,1)
		isogridbox.setRowStretch(4,1)
		isoframe.setLayout(isogridbox)
		gridbox.addWidget(isoframe, 2, 1, 2, 1)
		gridbox.addWidget(self.thr, 4, 0, 1, 2)
		
		# Set to default, but run only once and not in each base class
		if type(self) == EMIsosurfaceInspector: self.updateItemControls()
	
	def dataChanged(self):
		data = self.item3d().parent.data
		
		minden=data.get_attr("minimum")
		maxden=data.get_attr("maximum")
		mean=data.get_attr("mean")
		sigma=data.get_attr("sigma")
		iso_threshold = mean+3.0*sigma
		
		self.thr.setRange(minden,maxden)
		self.thr.setValue(iso_threshold, True)
		
		histogram_data = data.calc_hist(256,minden,maxden)
		self.histogram_widget.set_data(histogram_data,minden,maxden) 
		
		self.item3d().force_update = True
		self.item3d().isorender = MarchingCubes(data)
		self.setSamplingRange(self.item3d().isorender.get_sampling_range())

	def onCullFaces(self):
		self.item3d().cullbackfaces = self.cullbackface.isChecked()
		self.inspector.updateSceneGraph()
		
	def onThresholdSlider(self,val):
		self.item3d().setThreshold(val)
#		self.bright.setValue(-val,True)
		self.inspector.updateSceneGraph()
	
	def onSampling(self, val):
		self.item3d().setSample(val)
		self.inspector.updateSceneGraph()
	
	def onWireframe(self):
		self.item3d().wire = self.wireframe.isChecked()
		self.inspector.updateSceneGraph()
	
	def setSamplingRange(self,range):
		self.sampling_spinbox.setMinimum(1)
		self.sampling_spinbox.setMaximum(1+range-1)

class EMIsosurface(EMItem3D):
	name = "Isosurface"
	nodetype = "DataChild" 
	def __init__(self, parent=None, children = set(), transform = None):
		EMItem3D.__init__(self, parent, children, transform)
		
		self.isothr=0.5
		self.isodl = 0
		self.smpval=-1
		self.cube = False
		self.wire = False
		self.cullbackfaces = True
		self.tex_name = 0
		self.texture = False

#		self.brightness = 0
#		self.contrast = 10
#		self.rank = 1
		self.data_copy = None
		self.force_update = False
		self.loadColors()
	
		# color Needed for inspector to work John Flanagan
		self.diffuse = self.colors[self.isocolor]["diffuse"]
		self.specular = self.colors[self.isocolor]["specular"]
		self.ambient = self.colors[self.isocolor]["ambient"]		
		self.shininess = self.colors[self.isocolor]["shininess"]
		
	def setParent(self, parent):
		"""
		For use in session restore
		"""
		self.parent = parent
		if parent:
			data = self.parent.data
			assert isinstance(data, EMData)
			self.isorender = MarchingCubes(data)
		
	def getEvalString(self):
		return "EMIsosurface()"
		
	# I have added these methods so the inspector can set the color John Flanagan
	def setAmbientColor(self, red, green, blue, alpha=1.0):
		self.ambient = [red, green, blue, alpha]

	def setDiffuseColor(self, red, green, blue, alpha=1.0):
		self.diffuse = [red, green, blue, alpha]
		
	def setSpecularColor(self, red, green, blue, alpha=1.0):
		self.specular = [red, green, blue, alpha]
	
	def setShininess(self, shininess):
		self.shininess = shininess
		
	def getItemInspector(self):
		if not self.item_inspector:
			self.item_inspector = EMIsosurfaceInspector("ISOSURFACE", self)
		return self.item_inspector
	
	def loadColors(self):
		self.colors = get_default_gl_colors()
		
		self.isocolor = "bluewhite"

	def setSample(self,val):
		if ( self.smpval != int(val)):
			# the minus two is here because the marching cubes thinks -1 is the high level of detail, 0 is the next best and  so forth
			# However the user wants the highest level of detail to be 1, and the next best to be 2 and then 3 etc
			self.smpval = int(val)-2
			self.getIsosurfaceDisplayList()
		
	def getIsosurfaceDisplayList(self):
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
		cull = glIsEnabled(GL_CULL_FACE)
		polygonmode = glGetIntegerv(GL_POLYGON_MODE)
	
		if self.cullbackfaces:
			glEnable(GL_CULL_FACE)
			glCullFace(GL_BACK)
		else:
			if not cull:
				glDisable(GL_CULL_FACE)

		if ( self.wire ):
			glPolygonMode(GL_FRONT,GL_LINE);
		else:
			glPolygonMode(GL_FRONT,GL_FILL);
		

		glShadeModel(GL_SMOOTH)
		if ( self.isodl == 0 or self.force_update):
			self.getIsosurfaceDisplayList()
			self.force_update = False
		
		# This code draws an outline around the isosurface
		if self.is_selected:
			glPushAttrib( GL_ALL_ATTRIB_BITS )
		
			# First render the cylinder, writing the outline to the stencil buffer
			glClearStencil(0)
			glClear( GL_STENCIL_BUFFER_BIT )
			glEnable( GL_STENCIL_TEST )
			glStencilFunc( GL_ALWAYS, 1, 0xFFFF )		# Write to stencil buffer
			glStencilOp( GL_KEEP, GL_KEEP, GL_REPLACE )	# Only pixels that pass the depth test are written to the stencil buffer
			if ( self.wire ):
				glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
			else:
				glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);	
			self.renderIso()
		
			# Then render the outline
			glStencilFunc( GL_NOTEQUAL, 1, 0xFFFF )		# The object itself is stenciled out
			glStencilOp( GL_KEEP, GL_KEEP, GL_REPLACE )
			glLineWidth( 4.0 )				# By increasing the line width only the outline is drawn
			glPolygonMode( GL_FRONT_AND_BACK, GL_LINE )
			glMaterialfv(GL_FRONT, GL_EMISSION, [0.0, 1.0, 0.0, 1.0])
			self.renderIso()
	
			glPopAttrib()
		else:
			self.renderIso()
		
#		self.draw_bc_screen() #TODO: check into porting this from EM3DModel
		
		if self.cube:
			#glDisable(GL_LIGHTING)
			glPushMatrix()
			self.draw_volume_bounds()
			glPopMatrix()
		
		if cull: glEnable(GL_CULL_FACE)
		else: glDisable(GL_CULL_FACE)
		
		if ( polygonmode[0] == GL_LINE ): 
			glPolygonMode(GL_FRONT, GL_LINE)
		else: 
			glPolygonMode(GL_FRONT, GL_FILL)
		
		#print "total time is", time()-a
	
	def renderIso(self):
		# This is needed for the inspector to work John Flanagan	
		glMaterialfv(GL_FRONT, GL_DIFFUSE, self.diffuse)
		glMaterialfv(GL_FRONT, GL_SPECULAR, self.specular)
		glMaterialf(GL_FRONT, GL_SHININESS, self.shininess)
		glMaterialfv(GL_FRONT, GL_AMBIENT, self.ambient)
		glColor(self.ambient)

		
		glPushMatrix()
		glTranslate(-self.parent.data.get_xsize()/2.0,-self.parent.data.get_ysize()/2.0,-self.parent.data.get_zsize()/2.0)
		if ( self.texture ):
			glScalef(self.parent.data.get_xsize(),self.parent.data.get_ysize(),self.parent.data.get_zsize())
		glCallList(self.isodl)
		glPopMatrix()
	
	def setThreshold(self, val):
		if (self.isothr != val):
			self.isothr = val
#			self.brightness = -val
#			if ( self.texture ):
#				self.update_data_and_texture()
			self.getIsosurfaceDisplayList()
		
#			if self.emit_events: self.emit(QtCore.SIGNAL("set_threshold"),val)
#			self.updateGL()
		