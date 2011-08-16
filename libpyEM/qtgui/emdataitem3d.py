#!/usr/bin/env python
#
# Author: Ross Coleman (racolema@bcm.edu)
# Copyright (c) 2000-2011 Baylor College of Medicine

from OpenGL import GL
from OpenGL.GL import *
from PyQt4 import QtCore, QtGui
from libpyGLUtils2 import GLUtil
from EMAN2 import EMData, MarchingCubes, Transform
from emitem3d import EMItem3D, EMItem3DInspector, drawBoundingBox
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
		if self.transform:
			return "EMDataItem3D(\"%s\", transform=Transform())"%self.path
		else:
			return "EMDataItem3D(\"%s\")"%self.path
		
	def getItemInspector(self):
		if not self.item_inspector:
			self.item_inspector = EMDataItem3DInspector("DATA", self)
		return self.item_inspector
	
	def getData(self):
		return self.data
	
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
		
		for child in self.getChildren():
			child.dataChanged()
	
class EMDataItem3DInspector(EMItem3DInspector):
	def __init__(self, name, item3d):
		EMItem3DInspector.__init__(self, name, item3d)
	
	def updateItemControls(self):
		super(EMDataItem3DInspector, self).updateItemControls()
		# Anything that needs to be updated when the scene is rendered goes here.....
		self.file_path_label.setText(self.item3d().path)
		
	def addControls(self, gridbox):
		super(EMDataItem3DInspector, self).addControls(gridbox)
		hblbrowse = QtGui.QHBoxLayout()
		self.file_path_label = QtGui.QLabel()
		hblbrowse.addWidget(self.file_path_label)
		self.file_browse_button = QtGui.QPushButton("Browse")
		hblbrowse.addWidget(self.file_browse_button)
		gridbox.addLayout(hblbrowse, 3, 0, 1, 1)
		
		self.file_browse_button.clicked.connect(self.onFileBrowse)
		
		# Set to default, but run only once and not in each base class
		if type(self) == EMDataItem3DInspector: self.updateItemControls()
		
	def onFileBrowse(self):
		#TODO: replace this with an EMAN2 browser window once we re-write it
		file_path = QtGui.QFileDialog.getOpenFileName(self, "Open 3D Volume Map")
		self.file_path_label.setText(file_path)
		self.item3d().setData(file_path)
		self.inspector.updateSceneGraph()

class EMIsosurfaceInspector(EMInspectorControlShape):
	def __init__(self, name, item3d):
		EMItem3DInspector.__init__(self, name, item3d, numgridcols=2)	# for the iso inspector we need two grid cols for extra space....
		self.updateItemControls()
		
		QtCore.QObject.connect(self.thr, QtCore.SIGNAL("valueChanged"), self.onThresholdSlider)
		QtCore.QObject.connect(self.histogram_widget, QtCore.SIGNAL("thresholdChanged(float)"), self.onHistogram)
		self.cullbackface.toggled.connect(self.onCullFaces)
		self.wireframe.toggled.connect(self.onWireframe)
		self.sampling_spinbox.valueChanged[int].connect(self.onSampling)
	
	def updateItemControls(self):
		super(EMIsosurfaceInspector, self).updateItemControls()
		# Anything that needs to be updated when the scene is rendered goes here.....
		
		minden = self.item3d().minden
		maxden = self.item3d().maxden
		
		self.thr.setRange(minden, maxden)
		self.thr.setValue(self.item3d().isothr, True)
		
		self.histogram_widget.set_data(self.item3d().histogram_data,minden,maxden)
		self.setSamplingRange(self.item3d().isorender.get_sampling_range())
		
	def addControls(self, gridbox):
		super(EMIsosurfaceInspector, self).addControls(gridbox)
		self.histogram_widget = ImgHistogram(self, inithreshold=0.5)
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
	
	def onCullFaces(self):
		self.item3d().cullbackfaces = self.cullbackface.isChecked()
		self.inspector.updateSceneGraph()
		
	def onThresholdSlider(self,val):
		self.item3d().setThreshold(val)
		self.histogram_widget.setProbe(val)
#		self.bright.setValue(-val,True)
		self.inspector.updateSceneGraph()
	
	def onHistogram(self, value):
		self.thr.setValue(value)
		
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
		
		self.isothr = None #Will be set in self.dataChanged()
		self.isodl = 0
		self.smpval=-1
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
		
		if parent: self.dataChanged()
		
	def setParent(self, parent):
		"""
		For use in session restore
		"""
		self.parent = parent
		if parent:
			self.dataChanged()
		
	def getEvalString(self):
		if self.transform:
			return "EMIsosurface(transform=Transform())"
		else:
			return "EMIsosurface()"
	
	def dataChanged(self):
		data = self.getParent().getData()
		
		if self.isothr: #there was data previously
			normalized_threshold = (self.isothr - self.mean)/self.sigma
		else:
			normalized_threshold = 3.0
		
		self.minden = data.get_attr("minimum")
		self.maxden = data.get_attr("maximum")
		self.mean   = data.get_attr("mean")
		self.sigma  = data.get_attr("sigma")
		
		self.isothr = self.mean+normalized_threshold*self.sigma
		self.histogram_data = data.calc_hist(256,self.minden, self.maxden)
		
		self.force_update = True
		self.isorender = MarchingCubes(data)
		
		self.getItemInspector().updateItemControls()
	
		
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
			
			#Ross: these two lines add bounding box
			data = self.getParent().data
			drawBoundingBox(data.get_xsize(), data.get_ysize(), data.get_zsize())
			
			
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
		