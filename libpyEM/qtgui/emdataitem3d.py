#!/usr/bin/env python
#
# Author: Ross Coleman (racolema@bcm.edu)
# Author: James Michael Bell, 2016 (jmbell@bcm.edu)
# Copyright (c) 2011- Baylor College of Medicine
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

from EMAN2 import *
from embrowser import EMBrowserWidget
from emglobjects import EMViewportDepthTools, Camera2, get_default_gl_colors, get_RGB_tab, EM3DModel
from emglobjects import get_default_gl_colors
from emimageutil import ImgHistogram
from emitem3d import EMItem3D, EMItem3DInspector, drawBoundingBox
from emshapeitem3d import EMInspectorControlShape
from libpyGLUtils2 import GLUtil
import math
import os.path
import sys
from valslider import ValSlider, EMLightControls, CameraControls, EMSpinWidget, EMQTColorWidget, EMANToolButton

from OpenGL import GL
from OpenGL.GL import *
from OpenGL.GLU import *
from PyQt4 import QtCore, QtGui

import numpy as np

class EMDataItem3D(EMItem3D):
	"""
	This class is the scene graph node that has a reference to an EMData object. Though it has an orientation, it can not be displayed directly.
	Instead, its children are displayable, and are sized, postioned, and oriented relative to this node.
	"""
	name = "Data"

	@staticmethod
	def getNodeDialogWidget(attribdict):
		"""
		Get Data Widget
		"""
		datawidget = QtGui.QWidget()
		grid = QtGui.QGridLayout()
		node_name_data_label = QtGui.QLabel("Data Label")
		attribdict["node_name"] = QtGui.QLineEdit()
		data_path_label = QtGui.QLabel("Data Path")
		attribdict["data_path"] = QtGui.QLineEdit()
		browse_button = QtGui.QPushButton("Browse")
		grid.addWidget(node_name_data_label, 0, 0, 1, 2)
		grid.addWidget(attribdict["node_name"], 0, 2, 1, 2)
		grid.addWidget(data_path_label, 1, 0, 1, 2)
		grid.addWidget(attribdict["data_path"], 1, 2, 1, 2)
		grid.addWidget(browse_button, 2, 0, 1, 4)
		EMItem3D.get_transformlayout(grid, 4, attribdict)
		datawidget.setLayout(grid)

		EMDataItem3D.attribdict = attribdict
		QtCore.QObject.connect(browse_button, QtCore.SIGNAL('clicked()'), EMDataItem3D._on_browse)

		return datawidget

	@staticmethod
	def _on_browse():
		filename = QtGui.QFileDialog.getOpenFileName(None, 'Get file', os.getcwd())
		if filename:
			EMDataItem3D.attribdict["data_path"].setText(filename)
			name = os.path.basename(str(filename))
			EMDataItem3D.attribdict["node_name"].setText(str(name))

	@staticmethod
	def getNodeForDialog(attribdict):
		"""
		Create a new node using a attribdict
		"""
		return EMDataItem3D(str(attribdict["data_path"].text()), transform=EMItem3D.getTransformFromDict(attribdict))

	def __init__(self, data, parent = None, children = set(), transform=None, n=0):
		if not transform: transform = Transform()	# Object initialization should not be put in the constructor. Causes issues
		EMItem3D.__init__(self, parent, children, transform=transform)
		self.setData(data,n)
		self.renderBoundingBox = False

	def setSelectedItem(self, is_selected):
		""" Set SG apix to curent selection"""
		EMItem3D.setSelectedItem(self, is_selected)
		sg = self.getRootNode()
		try:
			if is_selected:
				sg.setAPix(self.data['apix_x'])
			else:
				sg.setAPix(None)
		except:
			print "ERROR in setting APIX"

	def getRenderBoundingBox(self):
		return self.renderBoundingBox

	def setRenderBoundingBox(self, state):
		self.renderBoundingBox = state

	def getEvalString(self):
		return "EMDataItem3D(\"%s\")"%os.path.abspath(self.path)

	def getItemInspector(self):
		if not self.item_inspector:
			self.item_inspector = EMDataItem3DInspector("DATA", self)
		return self.item_inspector

	def getBoundingBoxDimensions(self):
		data = self.getData()
		return (data["nx"], data["ny"], data["nz"])

	def getData(self):
		return self.data

	def setData(self, data, n=0):
		if isinstance(data, EMData):
			self.data = data
			if data.has_attr("source_path"):
				self.path = data["source_path"]
			else:
				self.path = None
		else:
			if n>0: self.data = EMData(str(data),n)
			else: self.data = EMData(str(data))
			self.path = str(data)

		for child in self.getChildren():
			try:	# What if child ins't a data type?
				child.dataChanged()
			except:
				pass

		# Get header info, get Apix info if this is available
		try:
			self.boundingboxsize = "%dx%dx%d\tApix: %5.2f"%(self.getBoundingBoxDimensions()[0],self.getBoundingBoxDimensions()[1],self.getBoundingBoxDimensions()[2],round(data['apix_x'],3))
		except:
			self.boundingboxsize = "%dx%dx%d"%(self.getBoundingBoxDimensions()[0],self.getBoundingBoxDimensions()[1],self.getBoundingBoxDimensions()[2])

		if self.item_inspector: self.item_inspector.updateMetaData()

	def getItemDictionary(self):
		"""
		Return a dictionary of item parameters (used for restoring sessions
		"""
		return super(EMDataItem3D, self).getItemDictionary()

	def renderNode(self):
		if self.renderBoundingBox: drawBoundingBox(*self.getBoundingBoxDimensions())

class EMDataItem3DInspector(EMItem3DInspector):
	def __init__(self, name, item3d):
		EMItem3DInspector.__init__(self, name, item3d)

	def updateItemControls(self):
		""" Updates this item inspector. Function is called by the item it observes"""
		super(EMDataItem3DInspector, self).updateItemControls()
		# Anything that needs to be updated when the scene is rendered goes here.....
		if self.item3d().path: self.file_path_label.setText(self.item3d().path)
		self.data_checkbox.setChecked(self.item3d().getRenderBoundingBox())


	def addTabs(self):
		""" Add a tab for each 'column' """

		tabwidget = QtGui.QWidget()
		gridbox = QtGui.QGridLayout()
		tabwidget.setLayout(gridbox)
		self.addTab(tabwidget, "data")
		# add data tab first, then basic
		super(EMDataItem3DInspector, self).addTabs()
		EMDataItem3DInspector.addControls(self, gridbox)

	def addControls(self, gridbox):
		""" Construct all the widgets in this Item Inspector """
		dataframe = QtGui.QFrame()
		dataframe.setFrameShape(QtGui.QFrame.StyledPanel)
		lfont = QtGui.QFont()
		lfont.setBold(True)
		datagridbox = QtGui.QGridLayout()

		self.data_checkbox= QtGui.QCheckBox("Display Bounding Box")
		datagridbox.addWidget(self.data_checkbox, 0, 0)
		self.file_browse_button = QtGui.QPushButton("Set Data Source")
		datagridbox.addWidget(self.file_browse_button, 1, 0)
		dataframe.setLayout(datagridbox)
		gridbox.addWidget(dataframe, 2, 0)

		self.file_path_label = QtGui.QLabel()
		self.file_path_label.setAlignment(QtCore.Qt.AlignCenter)
		self.file_path_label.setFont(lfont)
		gridbox.addWidget(self.file_path_label, 3, 0)

		self.file_browse_button.clicked.connect(self.onFileBrowse)
		QtCore.QObject.connect(self.data_checkbox, QtCore.SIGNAL("stateChanged(int)"), self.onBBoxChange)

		# Set to default, but run only once and not in each base class
		if type(self) == EMDataItem3DInspector: self.updateItemControls()

	def onFileBrowse(self):
		#TODO: replace this with an EMAN2 browser window once we re-write it
		file_path = QtGui.QFileDialog.getOpenFileName(self, "Open 3D Volume Map")
		if file_path:
			self.file_path_label.setText(file_path)
			self.item3d().setData(file_path)
			self.inspector().updateSceneGraph()

	def onBBoxChange(self, state):
		self.item3d().setRenderBoundingBox(not self.item3d().getRenderBoundingBox())
		self.inspector().updateSceneGraph()

class EMSliceItem3D(EMItem3D):
	"""
	This displays a slice of the volume that can be oriented any direction.
	Its parent in the tree data structure that forms the scene graph must be an EMDataItem3D instance.
	"""
	name = "Slice"
	nodetype = "DataChild"

	@staticmethod
	def getNodeDialogWidget(attribdict):
		"""
		Get Slice Widget
		"""
		slicewidget = QtGui.QWidget()
		grid = QtGui.QGridLayout()
		node_name_slice_label = QtGui.QLabel("Slice Name")
		attribdict["node_name"] = QtGui.QLineEdit(str(EMSliceItem3D.name))
		grid.addWidget(node_name_slice_label, 0, 0, 1, 2)
		grid.addWidget(attribdict["node_name"], 0, 2, 1, 2)
		EMItem3D.get_transformlayout(grid, 2, attribdict)
		slicewidget.setLayout(grid)

		return slicewidget

	@staticmethod
	def getNodeForDialog(attribdict):
		"""
		Create a new node using a attribdict
		"""
		return EMSliceItem3D(attribdict["parent"], transform=EMItem3D.getTransformFromDict(attribdict))

	def __init__(self, parent=None, children = set(), transform=None):
		"""
		@param parent: should be an EMDataItem3D instance for proper functionality.
		"""
		if not transform: transform = Transform()	# Object initialization should not be put in the constructor. Causes issues
		EMItem3D.__init__(self, parent, children, transform=transform)
		self.texture2d_name = 0
		self.texture3d_name = 0
		self.use_3d_texture = False
		self.force_texture_update = True

		self.colors = get_default_gl_colors()
		self.isocolor = "bluewhite"

		# color Needed for inspector to work John Flanagan
		self.diffuse = self.colors[self.isocolor]["diffuse"]
		self.specular = self.colors[self.isocolor]["specular"]
		self.ambient = self.colors[self.isocolor]["ambient"]
		self.shininess = self.colors[self.isocolor]["shininess"]

		if parent: self.dataChanged()

	# I have added these methods so the inspector can set the color John Flanagan
	def setAmbientColor(self, red, green, blue, alpha=1.0):
		self.ambient = [red, green, blue, alpha]

	def setDiffuseColor(self, red, green, blue, alpha=1.0):
		self.diffuse = [red, green, blue, alpha]

	def setSpecularColor(self, red, green, blue, alpha=1.0):
		self.specular = [red, green, blue, alpha]

	def setShininess(self, shininess):
		self.shininess = shininess

	def useDefaultBrightnessContrast(self):
		"""
		This applies default settings for brightness and contrast.
		"""

		data = self.getParent().getData()
		min = data.get_attr("minimum")
		max = data.get_attr("maximum")
		self.brightness = -min
		if max != min:
			self.contrast = 1.0/(max-min)
		else:
			self.contrast = 1

	def gen3DTexture(self):
		"""
		If no 3D texture exists, this creates one. It always returns the number that identifies the current 3D texture.
		"""
		if self.texture3d_name == 0:
			data_copy = self.getParent().getData().copy()
			data_copy.add(self.brightness)
			data_copy.mult(self.contrast)

			if True: #MUCH faster than generating texture in Python
				self.texture3d_name = GLUtil.gen_gl_texture(data_copy, GL.GL_LUMINANCE)
			else:
				self.texture3d_name = GL.glGenTextures(1)
				GL.glBindTexture(GL.GL_TEXTURE_3D, self.texture3d_name)
				GL.glTexImage3D(GL.GL_TEXTURE_3D,0,GL.GL_LUMINANCE, data_copy["nx"], data_copy["ny"], data_copy["nz"],0, GL.GL_ALPHA, GL.GL_FLOAT, data_copy.get_data_as_vector())

			#print "Slice Node's 3D texture == ", self.texture3d_name
		return self.texture3d_name

	def dataChanged(self):
		"""
		When the EMData changes for EMDataItem3D parent node, this method is called. It is responsible for updating the state of the slice node.
		"""
		if self.texture2d_name != 0:
			GL.glDeleteTextures(self.texture2d_name)
			self.texture2d_name = 0
		if self.texture3d_name != 0:
			GL.glDeleteTextures(self.texture3d_name)
			self.texture3d_name = 0

		self.useDefaultBrightnessContrast()

	def getEvalString(self):
		return "EMSliceItem3D()"

	def getItemInspector(self):
		if not self.item_inspector:
			self.item_inspector = EMSliceInspector("SLICE", self)
		return self.item_inspector

	def getItemDictionary(self):
		"""
		Return a dictionary of item parameters (used for restoring sessions
		"""
		dictionary = super(EMSliceItem3D, self).getItemDictionary()
		dictionary.update({"COLOR":[self.ambient, self.diffuse, self.specular, self.shininess]})
		return dictionary

	def setUsingDictionary(self, dictionary):
		"""
		Set item attributes using a dictionary, used in session restoration
		"""
		super(EMSliceItem3D, self).setUsingDictionary(dictionary)
		self.setAmbientColor(dictionary["COLOR"][0][0], dictionary["COLOR"][0][1], dictionary["COLOR"][0][2], dictionary["COLOR"][0][3])
		self.setDiffuseColor(dictionary["COLOR"][1][0], dictionary["COLOR"][1][1], dictionary["COLOR"][1][2], dictionary["COLOR"][1][3])
		self.setSpecularColor(dictionary["COLOR"][2][0], dictionary["COLOR"][2][1], dictionary["COLOR"][2][2], dictionary["COLOR"][2][3])

	def renderNode(self):
		data = self.getParent().getData()

		nx = data["nx"]
		ny = data["ny"]
		nz = data["nz"]
		interior_diagonal = math.sqrt(nx**2+ny**2+nz**2) #A square with sides this big could hold any slice from the volume
		#The interior diagonal is usually too big, and OpenGL textures work best with powers of 2 so let's get the next smaller power of 2
		diag = 2**(int(math.floor( math.log(interior_diagonal)/math.log(2) ))) #next smaller power of 2
		diag2 = diag/2

		glPushAttrib( GL_ALL_ATTRIB_BITS )
		GL.glDisable(GL.GL_LIGHTING)
		GL.glColor3f(1.0,1.0,1.0)

		if not self.use_3d_texture: #Use 2D texture

			# Any time self.transform changes, a new 2D texture is REQUIRED.
			# It is easiest to create a new 2D texture every time this is called, and seems fast enough.
			# Thus, a new 2D texture is created whether it is needed or not.

			temp_data = EMData(diag, diag)
			temp_data.cut_slice(data, self.transform)
			temp_data.add(self.brightness)
			temp_data.mult(self.contrast)

			if self.texture2d_name != 0:
				GL.glDeleteTextures(self.texture2d_name)

			self.texture2d_name = GLUtil.gen_gl_texture(temp_data, GL.GL_LUMINANCE)


			#For debugging purposes, draw an outline
			GL.glBegin(GL.GL_LINE_LOOP)
			GL.glVertex3f(-diag2, -diag2, 0)
			GL.glVertex3f(-diag2, diag2, 0)
			GL.glVertex3f(diag2, diag2, 0)
			GL.glVertex3f(diag2, -diag2, 0)
			GL.glEnd()


			#Now draw the texture on another quad

			GL.glEnable(GL.GL_TEXTURE_2D)
			GL.glBindTexture(GL.GL_TEXTURE_2D, self.texture2d_name)
			GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_S, GL.GL_CLAMP)
			GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_T, GL.GL_CLAMP)
			GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_LINEAR)#GL.GL_NEAREST)
			GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_LINEAR)#GL.GL_NEAREST)

			GL.glTexEnvf(GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_REPLACE)

			GL.glBegin(GL.GL_QUADS)
			GL.glTexCoord2f(0, 0)
			GL.glVertex3f(-diag2, -diag2, 0)
			GL.glTexCoord2f(0, 1)
			GL.glVertex3f(-diag2, diag2, 0)
			GL.glTexCoord2f(1, 1)
			GL.glVertex3f(diag2, diag2, 0)
			GL.glTexCoord2f(1, 0)
			GL.glVertex3f(diag2, -diag2, 0)
			glEnd()

			GL.glDisable(GL.GL_TEXTURE_2D)

		else: #Using a 3D texture

			# Generating a new 3D texture is slower than a new 2D texture.
			# Creating a new texture is needed if brightness or contrast change or if the EMData in self.getParent().getData() has changed.

			if self.force_texture_update:
				GL.glDeleteTextures(self.texture3d_name)
				self.texture3d_name = 0
			self.gen3DTexture()

			quad_points = [(-diag2, -diag2, 0), (-diag2, diag2, 0), (diag2, diag2, 0), (diag2, -diag2, 0)]

			#For debugging purposes, draw an outline
			GL.glMatrixMode(GL.GL_MODELVIEW)
			GL.glBegin(GL.GL_LINE_LOOP)
			for i in range(4):
				GL.glVertex3f(*quad_points[i])
			GL.glEnd()

			#Now draw the texture on another quad

			GL.glEnable(GL.GL_TEXTURE_3D)
			GL.glBindTexture(GL.GL_TEXTURE_3D, self.texture3d_name)
			GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_WRAP_S, GL.GL_CLAMP_TO_BORDER)
			GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_WRAP_T, GL.GL_CLAMP_TO_BORDER)
			GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_WRAP_R, GL.GL_CLAMP_TO_BORDER)
			GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_LINEAR)#GL.GL_NEAREST)
			GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_LINEAR)#GL.GL_NEAREST)
	#		GL.glTexParameterfv(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_BORDER_COLOR, (0.0, 0.0,1.0,1.0,1.0))

			GL.glTexEnvf(GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_REPLACE)
			GL.glMatrixMode(GL.GL_TEXTURE)
			GL.glLoadIdentity()
			GL.glTranslatef(0.5, 0.5, 0.5) #Put the origin at the center of the 3D texture
			GL.glScalef(1.0/nx, 1.0/ny, 1.0/nz) #Scale to make the texture coords the same as data coords
			GLUtil.glMultMatrix(self.transform) #Make texture coords the same as EMSliceItem3D coords
			GL.glMatrixMode(GL.GL_MODELVIEW)

			GL.glBegin(GL.GL_QUADS)
			for i in range(4):
				GL.glTexCoord3f(*quad_points[i])
				GL.glVertex3f(*quad_points[i])
			glEnd()

			GL.glMatrixMode(GL.GL_TEXTURE)
			GL.glLoadIdentity()
			GL.glMatrixMode(GL.GL_MODELVIEW)

			GL.glDisable(GL.GL_TEXTURE_3D)

		GL.glEnable(GL.GL_LIGHTING)
		glPopAttrib()

class EMSliceInspector(EMInspectorControlShape):
	def __init__(self, name, item3d):
		EMInspectorControlShape.__init__(self, name, item3d)

		self.constrained_plane_combobox.currentIndexChanged.connect(self.onConstrainedOrientationChanged)
		self.use_3d_texture_checkbox.clicked.connect(self.on3DTextureCheckbox)
		QtCore.QObject.connect(self.constrained_slider, QtCore.SIGNAL("valueChanged"), self.onConstraintSlider)
		QtCore.QObject.connect(self.brightness_slider, QtCore.SIGNAL("valueChanged"), self.onBrightnessSlider)
		QtCore.QObject.connect(self.contrast_slider, QtCore.SIGNAL("valueChanged"), self.onContrastSlider)

		self.updateItemControls()

	def updateItemControls(self):
		""" Updates this item inspector. Function is called by the item it observes"""
		super(EMSliceInspector, self).updateItemControls()
		# Anything that needs to be updated when the scene is rendered goes here.....
		self.use_3d_texture_checkbox.setChecked(self.item3d().use_3d_texture)
		data = self.item3d().getParent().getData()
		min = data["minimum"]
		max = data["maximum"]
		mean = data["mean"]
		std_dev = data["sigma"]

		self.brightness_slider.setValue(self.item3d().brightness)
		self.brightness_slider.setRange(-max, -min)

		self.contrast_slider.setValue(self.item3d().contrast)
		self.contrast_slider.setRange(0.001, 1.0)

	def addTabs(self):
		""" Add a tab for each 'column' """
		tabwidget = QtGui.QWidget()
		gridbox = QtGui.QGridLayout()
		tabwidget.setLayout(gridbox)
		self.addTab(tabwidget, "slices")
		# add slices tab first then basic tab
		super(EMSliceInspector, self).addTabs()
		EMSliceInspector.addControls(self, gridbox)

	def addControls(self, gridbox):
		""" Construct all the widgets in this Item Inspector """
		sliceframe = QtGui.QFrame()
		sliceframe.setFrameShape(QtGui.QFrame.StyledPanel)
		slice_grid_layout = QtGui.QGridLayout()

		self.constrained_group_box = QtGui.QGroupBox("Constrained Slices")
		self.constrained_group_box.setCheckable(True)
		self.constrained_group_box.setChecked(False)

		self.constrained_plane_combobox = QtGui.QComboBox()
		self.constrained_plane_combobox.addItems(["XY", "YZ", "ZX"])
		self.constrained_slider = ValSlider(label="Trans:")

		constrained_layout = QtGui.QVBoxLayout()
		constrained_layout.addWidget(self.constrained_plane_combobox)
		constrained_layout.addWidget(self.constrained_slider)
		constrained_layout.addStretch()
		self.constrained_group_box.setLayout(constrained_layout)

		self.use_3d_texture_checkbox = QtGui.QCheckBox("Use 3D Texture")
		self.use_3d_texture_checkbox.setChecked(self.item3d().use_3d_texture)

		self.brightness_slider = ValSlider(label="Bright:")
		self.contrast_slider = ValSlider(label="Contr:")

		slice_grid_layout.addWidget(self.constrained_group_box, 0, 1, 2, 1)
		slice_grid_layout.addWidget(self.use_3d_texture_checkbox, 2, 1, 1, 1)
		slice_grid_layout.addWidget(self.brightness_slider, 3, 1, 1, 1)
		slice_grid_layout.addWidget(self.contrast_slider, 4, 1, 1, 1)
		slice_grid_layout.setRowStretch(5,1)
		sliceframe.setLayout(slice_grid_layout)
		gridbox.addWidget(sliceframe, 2, 0, 2, 1)

	def on3DTextureCheckbox(self):
		self.item3d().use_3d_texture = self.use_3d_texture_checkbox.isChecked()
		if self.inspector:
			self.inspector().updateSceneGraph()

	def onConstrainedOrientationChanged(self):
		self.constrained_slider.setValue(0)
		(nx, ny, nz) = self.item3d().getParent().getBoundingBoxDimensions()
		range = (0, nx)
		plane = str(self.constrained_plane_combobox.currentText())
		if plane == "XY": range = (-nz/2.0, nz/2.0)
		elif plane == "YZ": range = (-nx/2.0, nx/2.0)
		elif plane == "ZX": range = (-ny/2.0, ny/2.0)
		self.constrained_slider.setRange(*range)
		self.onConstraintSlider()

	def onConstraintSlider(self):
		value = self.constrained_slider.getValue()
		transform = self.item3d().getTransform()
		plane = str(self.constrained_plane_combobox.currentText())
		if plane == "XY":
			transform.set_rotation((0,0,1))
			transform.set_trans(0,0,value)
		elif plane == "YZ":
			transform.set_rotation((1,0,0))
			transform.set_trans(value, 0, 0)
		elif plane == "ZX":
			transform.set_rotation((0,1,0))
			transform.set_trans(0,value,0)

		if self.inspector:
			self.inspector().updateSceneGraph()

	def onBrightnessSlider(self):
		self.item3d().brightness = self.brightness_slider.getValue()
		self.item3d().force_texture_update = True
		if self.inspector:
			self.inspector().updateSceneGraph()

	def onContrastSlider(self):
		self.item3d().contrast = self.contrast_slider.getValue()
		self.item3d().force_texture_update = True
		if self.inspector:
			self.inspector().updateSceneGraph()

class EMVolumeItem3D(EMItem3D):
	"""
	This displays a partially transparent view through the EMData.
	Its parent in the tree data structure that forms the scene graph must be an EMDataItem3D instance.
	"""

	name = "Volume"
	nodetype = "DataChild"

	@staticmethod
	def getNodeDialogWidget(attribdict):
		"""
		Get Volume Widget
		"""
		volumewidget = QtGui.QWidget()
		grid = QtGui.QGridLayout()
		node_name_volume_label = QtGui.QLabel("Volume Name")
		attribdict["node_name"] = QtGui.QLineEdit(str(EMVolumeItem3D.name))
		grid.addWidget(node_name_volume_label, 0, 0, 1, 2)
		grid.addWidget(attribdict["node_name"], 0, 2, 1, 2)
		EMItem3D.get_transformlayout(grid, 2, attribdict)
		volumewidget.setLayout(grid)

		return volumewidget

	@staticmethod
	def getNodeForDialog(attribdict):
		"""
		Create a new node using a attribdict
		"""
		return EMVolumeItem3D(attribdict["parent"], transform=EMItem3D.getTransformFromDict(attribdict))

	def __init__(self, parent=None, children = set(), transform=None):
		"""
		@param parent: should be an EMDataItem3D instance for proper functionality.
		"""
		if not transform: transform = Transform()	# Object initialization should not be put in the constructor. Causes issues
		EMItem3D.__init__(self, parent, children, transform=transform)

		self.texture_name = 0
		self.colors = get_default_gl_colors()
		self.isocolor = "bluewhite"
		self.isothr = None #Will be set in self.dataChanged()

		# color Needed for inspector to work John Flanagan
		self.diffuse = self.colors[self.isocolor]["diffuse"]
		self.specular = self.colors[self.isocolor]["specular"]
		self.ambient = self.colors[self.isocolor]["ambient"]
		self.shininess = self.colors[self.isocolor]["shininess"]

		self.force_texture_update = True
		if parent: self.dataChanged()

	# I have added these methods so the inspector can set the color John Flanagan
	def setAmbientColor(self, red, green, blue, alpha=1.0):
		self.ambient = [red, green, blue, alpha]

	def setDiffuseColor(self, red, green, blue, alpha=1.0):
		self.diffuse = [red, green, blue, alpha]

	def setSpecularColor(self, red, green, blue, alpha=1.0):
		self.specular = [red, green, blue, alpha]

	def setShininess(self, shininess):
		self.shininess = shininess

	def dataChanged(self):
		"""
		When the EMData changes for EMDataItem3D parent node, this method is called. It is responsible for updating the state of the slice node.
		"""

		data = self.getParent().getData()

		self.minden = data.get_attr("minimum")
		self.maxden = data.get_attr("maximum")
		self.mean   = data.get_attr("mean")
		self.sigma  = data.get_attr("sigma")

		if self.isothr: #there was data previously
			normalized_threshold = (self.isothr - self.mean)/self.sigma
		else:
			normalized_threshold = 3.0

		self.isothr = self.mean+normalized_threshold*self.sigma

		self.histogram_data = data.calc_hist(256,self.minden, self.maxden)

		if self.texture_name != 0:
			GL.glDeleteTextures(self.texture_name)
		self.useDefaultBrightnessThreshold()

	def getEvalString(self):
		return "EMVolumeItem3D()"

	def getItemInspector(self):
		if not self.item_inspector:
			self.item_inspector = EMVolumeInspector("VOLUME", self)
		return self.item_inspector

	def getItemDictionary(self):
		"""
		Return a dictionary of item parameters (used for restoring sessions
		"""
		dictionary = super(EMVolumeItem3D, self).getItemDictionary()
		dictionary.update({"COLOR":[self.ambient, self.diffuse, self.specular, self.shininess]})
		return dictionary

	def setUsingDictionary(self, dictionary):
		"""
		Set item attributes using a dictionary, used in session restoration
		"""
		super(EMVolumeItem3D, self).setUsingDictionary(dictionary)
		self.setAmbientColor(dictionary["COLOR"][0][0], dictionary["COLOR"][0][1], dictionary["COLOR"][0][2], dictionary["COLOR"][0][3])
		self.setDiffuseColor(dictionary["COLOR"][1][0], dictionary["COLOR"][1][1], dictionary["COLOR"][1][2], dictionary["COLOR"][1][3])
		self.setSpecularColor(dictionary["COLOR"][2][0], dictionary["COLOR"][2][1], dictionary["COLOR"][2][2], dictionary["COLOR"][2][3])

	def gen3DTexture(self):
		"""
		Creates 3D texture if none exists. Returns the number that identifies the current texture.
		"""
		if self.texture_name == 0:
			data_copy = self.getParent().getData().copy()
			data_copy.add(self.brightness)
			data_copy.mult(self.contrast)

			if True: #MUCH faster than generating texture in Python
				self.texture_name = GLUtil.gen_gl_texture(data_copy, GL.GL_ALPHA)
			else:
				self.texture_name = GL.glGenTextures(1)
				GL.glBindTexture(GL.GL_TEXTURE_3D, self.texture_name)
				GL.glTexImage3D(GL.GL_TEXTURE_3D,0, GL.GL_ALPHA,data_copy["nx"], data_copy["ny"], data_copy["nz"],0, GL.GL_ALPHA, GL.GL_FLOAT, data_copy.get_data_as_vector())

			#print "Volume Node's 3D texture ==", self.texture_name
		return self.texture_name

	def useDefaultBrightnessThreshold(self):
		"""
		Applies default brightness and contrast.
		"""

		data = self.getParent().getData()
		min = data["minimum"]
		max = data["maximum"]

		self.brightness = 0 #-main
		if max != min:
			self.contrast = 1.0/(max - min)
		else:
			self.contrast = 1

	def renderNode(self):
		"""
		This creates a 3D texture if needed. Then, it uses transformed texture coordinates to apply the needed image data to a stack of rectangles that are slices of the volume.
		The rectangles are created parallel to the plane of the screen. Each one has partly transparent texture data applied to it.
		The rectangles are created from back to front (near the viewer) because a background color is needed before the newest slice can be blended with it.
		"""

		if self.force_texture_update:
			GL.glDeleteTextures(self.texture_name)
			self.texture_name = 0

		self.gen3DTexture()
		self.force_texture_update = False
		(nx, ny, nz) = self.getParent().getBoundingBoxDimensions()
		interior_diagonal = math.sqrt(nx**2+ny**2+nz**2) #A square with sides this big could hold any slice from the volume
		diag2 = interior_diagonal/2

		quad_points = [(-diag2, -diag2, 0), (-diag2, diag2, 0), (diag2, diag2, 0), (diag2, -diag2, 0)]
		transform_std_coords = self.getTransformStdCoord()


		#Find root node
		rootNode = self
		while rootNode.getParent():
			rootNode = rootNode.getParent()

		GL.glMatrixMode(GL.GL_MODELVIEW)
		GL.glPushMatrix()
		GL.glLoadIdentity() #Make the polygon be in the plane of the screen
		transf = Transform(transform_std_coords) # Get a copy, not a reference, need to use original later
		transf.set_rotation({"type":"eman"}) #Removing rotation
		GLUtil.glMultMatrix(transf)
		rootNode.camera.setCameraPosition()

		#For debugging purposes, draw an outline
		GL.glDisable(GL.GL_LIGHTING)
		GL.glColor3f(1.0,1.0,1.0)
		GL.glBegin(GL.GL_LINE_LOOP)
		for i in range(4):
			GL.glVertex3f(*quad_points[i])
		GL.glEnd()

		#Now draw the texture on another quad
		GL.glEnable(GL.GL_TEXTURE_3D)
		GL.glBindTexture(GL.GL_TEXTURE_3D, self.texture_name)
		GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_WRAP_S, GL.GL_CLAMP_TO_BORDER)
		GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_WRAP_T, GL.GL_CLAMP_TO_BORDER)
		GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_WRAP_R, GL.GL_CLAMP_TO_BORDER)
		GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_LINEAR)#GL.GL_NEAREST)
		GL.glTexParameterf(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_LINEAR)#GL.GL_NEAREST)
#		GL.glTexParameterfv(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_BORDER_COLOR, (0.0, 0.0,1.0,1.0,1.0))

		GL.glMatrixMode(GL.GL_TEXTURE)
		GL.glLoadIdentity()
		GL.glTranslatef(0.5, 0.5, 0.5) #Put the origin at the center of the 3D texture
		GL.glScalef(1.0/nx, 1.0/ny, 1.0/nz) #Scale to make the texture coords the same as data coords

		transf = Transform(transform_std_coords) # Get a copy, not a reference, in case we need to use original later
		transf.set_trans((0,0,0)) #Removing translation
		transf.set_scale(1.0) #Removing scaling
		GLUtil.glMultMatrix(transf.inverse())
		GL.glMatrixMode(GL.GL_MODELVIEW)

		# Draw the volume by drawing a stack of semi-transparent slices oriented in the plane of the screen. Order is important: back to front.
		GL.glEnable(GL.GL_BLEND)
		GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)
		GL.glTexEnvf(GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_REPLACE)

		GL.glBegin(GL.GL_QUADS)
		half_depth = int(diag2)
		for dz in range(-half_depth, half_depth+1):
			for i in range(4):
				coords = quad_points[i]
				GL.glTexCoord3f(coords[0], coords[1], coords[2] + dz)
				GL.glVertex3f(coords[0], coords[1], coords[2] + dz)
		glEnd()

		GL.glPopMatrix()

		GL.glMatrixMode(GL.GL_TEXTURE)
		GL.glLoadIdentity() #Return texture coord system to default
		GL.glMatrixMode(GL.GL_MODELVIEW)

		GL.glDisable(GL.GL_TEXTURE_3D)
		GL.glEnable(GL.GL_LIGHTING)
		GL.glDisable(GL.GL_BLEND)

class EMVolumeInspector(EMInspectorControlShape):
	def __init__(self, name, item3d):

		EMInspectorControlShape.__init__(self, name, item3d)

		QtCore.QObject.connect(self.brightness_slider, QtCore.SIGNAL("valueChanged"), self.onBrightnessSlider)
		QtCore.QObject.connect(self.contrast_slider, QtCore.SIGNAL("valueChanged"), self.onContrastSlider)

	def updateItemControls(self):
		""" Updates this item inspector. Function is called by the item it observes"""
		super(EMVolumeInspector, self).updateItemControls()
		# Anything that needs to be updated when the scene is rendered goes here.....

		data = self.item3d().getParent().getData()
		min = data["minimum"]
		max = data["maximum"]


		self.minden = self.item3d().minden
		self.maxden = self.item3d().maxden

		#self.thr.setRange(minden, maxden)
		#self.thr.setValue(self.item3d().isothr, True)

		self.histogram_widget.set_data(self.item3d().histogram_data,self.minden,self.maxden)

		self.brightness_slider.setValue(self.item3d().brightness)
		self.brightness_slider.setRange(-max, -min)
		self.contrast_slider.setValue(self.item3d().contrast)
		self.contrast_slider.setRange(0.001, 1.0)

	def addTabs(self):
		""" Add a tab for each 'column' """
		tabwidget = QtGui.QWidget()
		gridbox = QtGui.QGridLayout()
		tabwidget.setLayout(gridbox)
		self.addTab(tabwidget, "volume")
		# add volume tab then basic tab
		super(EMVolumeInspector, self).addTabs()
		EMVolumeInspector.addControls(self, gridbox)

	def addControls(self, gridbox):
		""" add controls for the volumes """
		self.histogram_widget = ImgHistogram(self)
		self.histogram_widget.setObjectName("hist")

		volframe = QtGui.QFrame()
		volframe.setFrameShape(QtGui.QFrame.StyledPanel)
		vol_grid_layout = QtGui.QGridLayout()

		probeframe = QtGui.QFrame()
		probeframe.setFrameShape(QtGui.QFrame.StyledPanel)
		probelayout = QtGui.QGridLayout()
		probelayout.setAlignment(QtCore.Qt.AlignTop)
		self.range = QtGui.QLabel("Range: %1.3f, %1.3f"%(self.item3d().minden,self.item3d().maxden))
		self.level = QtGui.QLabel("Level: %1.3f"%self.item3d().isothr)
		self.level.setFixedSize(100,40)
		self.color = QtGui.QLabel("Color:")
		self.cappingcolor = EMQTColorWidget(parent=self)

		probelayout.addWidget(self.range, 0, 0, 1, 1)
		probelayout.addWidget(self.level, 0, 1, 1, 1)
		probelayout.addWidget(self.color, 0, 2, 1, 1)
		probelayout.addWidget(self.cappingcolor, 0, 3, 1, 1)
		probeframe.setLayout(probelayout)

#		paramframe = QtGui.QFrame()
		self.brightness_slider = ValSlider(label="Bright:")
		self.contrast_slider = ValSlider(label="Contr:")

		vol_grid_layout.addWidget(self.histogram_widget, 0, 0, 1, 1)
		vol_grid_layout.addWidget(probeframe, 1,0,1,1,)
		vol_grid_layout.addWidget(self.brightness_slider, 2, 0, 1, 1)
		vol_grid_layout.addWidget(self.contrast_slider, 3, 0, 1, 1)

		self.level.setAlignment(QtCore.Qt.AlignCenter)
		vol_grid_layout.setRowStretch(4,1)
		vol_grid_layout.setColumnStretch(1,10)
		vol_grid_layout.setColumnStretch(2,10)

		volframe.setLayout(vol_grid_layout)
		gridbox.addWidget(volframe, 2, 0, 2, 1)
		if type(self) == EMVolumeInspector: self.updateItemControls()

		###
		transp=0.85
		self.probeposition=[[self.item3d().isothr,0],[(self.item3d().isothr+self.item3d().maxden)/2,transp],[self.item3d().maxden,transp]]
		self.probepresent=0;

		levelvalue =  self.probeposition[self.probepresent][0]

		if (self.maxden-self.minden) != 0:
			for i in range(3):
				self.probeposition[i][0] = int(255.0*(self.probeposition[i][0] - self.minden)/(self.maxden-self.minden))
				self.probeposition[i][1] = int(125 - 125.0*(self.probeposition[i][1]))
				if self.probeposition[i][0] > 255: self.probeposition[i][0] = 255
				if self.probeposition[i][0] < 3: self.probeposition[i][0] = 3
				if self.probeposition[i][1] > 125: self.probeposition[i][1]=125
				if self.probeposition[i][1] < 3 : self.probeposition[i][1] =3
		self.probecolor=[[255,255,255],[255,255,255],[255,255,255]]

		self.histogram_widget.setDynamicProbe(self.probeposition, self.probecolor, self.probepresent,levelvalue) # The needs to be node AFTER the data is set

		QtCore.QObject.connect(self.cappingcolor,QtCore.SIGNAL("newcolor(QColor)"),self._on_cap_color)

	def _on_cap_color(self, color):
		rgb = color.getRgb()
		self.histogram_widget.setColor(rgb)

		'''
		self.scenegraph().makeCurrent()
		self.scenegraph().camera.setCapColor(float(rgb[0])/255.0, float(rgb[1])/255.0, float(rgb[2])/255.0)
		if self.scenegraph().camera.getCappingMode(): self.updateSceneGraph()
		'''

	def onBrightnessSlider(self):
		self.item3d().brightness = self.brightness_slider.getValue()
		self.item3d().force_texture_update = True
		if self.inspector:
			self.inspector().updateSceneGraph()

	def onContrastSlider(self):
		self.item3d().contrast = self.contrast_slider.getValue()
		self.item3d().force_texture_update = True
		if self.inspector:
			self.inspector().updateSceneGraph()


class EMIsosurfaceInspector(EMInspectorControlShape):
	def __init__(self, name, item3d):

		EMInspectorControlShape.__init__(self, name, item3d)	# for the iso inspector we need two grid cols for extra space....

		QtCore.QObject.connect(self.thr, QtCore.SIGNAL("valueChanged"), self.onThresholdSlider)
		QtCore.QObject.connect(self.histogram_widget, QtCore.SIGNAL("thresholdChanged(float)"), self.onHistogram)
		self.cullbackface.toggled.connect(self.onCullFaces)
		self.wireframe.toggled.connect(self.onWireframe)
		self.colorbyradius.toggled.connect(self.onColorByRadius)
		self.colorbymap.toggled.connect(self.onColorByMap)
		self.cmapbrowse.clicked.connect(self.onFileBrowse)
		self.sampling_spinbox.valueChanged[int].connect(self.onSampling)
		QtCore.QObject.connect(self.innercolorscaling, QtCore.SIGNAL("valueChanged(int)"), self.reColorScale)
		QtCore.QObject.connect(self.outercolorscaling, QtCore.SIGNAL("valueChanged(int)"), self.reColorScale)
		QtCore.QObject.connect(self.cmapmin, QtCore.SIGNAL("valueChanged(int)"), self.reColorMapMinMax)
		QtCore.QObject.connect(self.cmapmax, QtCore.SIGNAL("valueChanged(int)"), self.reColorMapMinMax)

	def updateItemControls(self):
		""" Updates this item inspector. Function is called by the item it observes"""
		super(EMIsosurfaceInspector, self).updateItemControls()
		# Anything that needs to be updated when the scene is rendered goes here.....

		minden = self.item3d().minden
		maxden = self.item3d().maxden

		self.thr.setRange(minden, maxden)
		self.thr.setValue(self.item3d().isothr, True)

		self.histogram_widget.set_data(self.item3d().histogram_data,minden,maxden)
		self.setSamplingRange(self.item3d().isorender.get_sampling_range())

		#Set color radius
		if self.item3d().rgbmode == 1:
			self.colorbyradius.setChecked(True)
			self.colorbymap.setChecked(False)
		self.innercolorscaling.setValue(self.item3d().innerrad)
		self.outercolorscaling.setValue(self.item3d().outerrad)

		# Colormap data
		if self.item3d().rgbmode == 2:
			self.colorbyradius.setChecked(False)
			self.colorbymap.setChecked(True)
		colormapdata = self.item3d().cmapdata
		if colormapdata:
			self.cmapmin.setValue(self.item3d().cmapmin, quiet=1)
			self.cmapmax.setValue(self.item3d().cmapmax, quiet=1)
			cmrange = colormapdata.get_attr('maximum') - colormapdata.get_attr('minimum')
			self.cmapmin.setIncrement(cmrange/50.0)
			self.cmapmax.setIncrement(cmrange/50.0)
			rounding = int(math.ceil(math.fabs(math.log(cmrange/2.0)))+1)
			#print rounding
			self.cmapmin.setRounding(rounding)
			self.cmapmax.setRounding(rounding)
			self.colormap.setText(self.item3d().cmapfilename)
			if str(self.colormap.text()) != "":
				self.colorbymap.setEnabled(True)
			else:
				self.colorbymap.setEnabled(False)

	def addTabs(self):
		""" Add a tab for each 'column' """
		tabwidget = QtGui.QWidget()
		gridbox = QtGui.QGridLayout()
		tabwidget.setLayout(gridbox)
		self.addTab(tabwidget, "isosurface")
		# add isosurface tab first, then basic tab
		super(EMIsosurfaceInspector, self).addTabs()
		EMIsosurfaceInspector.addControls(self, gridbox)

	def addControls(self, gridbox):
		""" Construct all the widgets in this Item Inspector """
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
		self.thr.setValue(self.item3d().isothr)
		self.sampling_label = QtGui.QLabel("Sample Level:")
		self.sampling_spinbox = QtGui.QSpinBox()
		self.sampling_spinbox.setValue(1)

		sampling_hbox_layout = QtGui.QHBoxLayout()
		sampling_hbox_layout.addWidget(self.sampling_label)
		sampling_hbox_layout.addWidget(self.sampling_spinbox)

		# Color by radius frame
		cbrframe = QtGui.QFrame()
		cbrframe.setFrameShape(QtGui.QFrame.StyledPanel)
		cbrlayout = QtGui.QGridLayout()
		cbrlayout.setAlignment(QtCore.Qt.AlignTop)
		self.colorbyradius = QtGui.QCheckBox("Color By Radius")
		self.colorbyradius.setChecked(False)
		self.innercolorscaling = EMSpinWidget(0.0, 1.0, rounding=0)
		self.innercolorscaling.setMinimumWidth(120)
		self.outercolorscaling = EMSpinWidget(0.0, 1.0, rounding=0)
		self.outercolorscaling.setMinimumWidth(120)
		cbrlayout.addWidget(self.colorbyradius, 0, 0, 1, 2)
		cbrlabelinner = QtGui.QLabel("Inner Shell")
		cbrlabelinner.setAlignment(QtCore.Qt.AlignCenter)
		cbrlabelouter = QtGui.QLabel("Outer Shell")
		cbrlabelouter.setAlignment(QtCore.Qt.AlignCenter)
		cbrlayout.addWidget(cbrlabelinner, 1, 0, 1, 1)
		cbrlayout.addWidget(cbrlabelouter, 1, 1, 1, 1)
		cbrlayout.addWidget(self.innercolorscaling, 2, 0, 1, 1)
		cbrlayout.addWidget(self.outercolorscaling, 2, 1, 1, 1)
		cbrframe.setLayout(cbrlayout)

		# Color by map frame
		cbmframe = QtGui.QFrame()
		cbmframe.setFrameShape(QtGui.QFrame.StyledPanel)
		cbmlayout = QtGui.QGridLayout()
		cbmlayout.setAlignment(QtCore.Qt.AlignTop)
		self.colorbymap = QtGui.QCheckBox("Color By Map")
		self.colorbymap.setEnabled(False)
		self.colormap = QtGui.QLineEdit("")
		self.cmapbrowse = QtGui.QPushButton("Browse")
		self.cmapmin = EMSpinWidget(0.0, 0.1, rounding=2)
		self.cmapmax = EMSpinWidget(0.0, 0.1, rounding=2)
		cbmlayout.addWidget(self.colorbymap, 0, 0, 1, 2)
		cbmlayout.addWidget(self.colormap, 1, 0)
		cbmlayout.addWidget(self.cmapbrowse, 1, 1)
		cbmlayout.addWidget(self.cmapmin, 2, 0)
		cbmlayout.addWidget(self.cmapmax, 2, 1)
		cbmframe.setLayout(cbmlayout)


		isogridbox.addWidget(self.histogram_widget, 0, 0, 1, 1)
		isogridbox.addWidget(self.cullbackface, 1, 0, 1, 1)
		isogridbox.addWidget(self.wireframe, 2,0,1,1)
		isogridbox.addLayout(sampling_hbox_layout, 3,0,1,1)
		isogridbox.addWidget(cbrframe, 4, 0, 1, 1)
		isogridbox.addWidget(cbmframe, 5, 0, 1, 1)
		isogridbox.setRowStretch(4,1)
		isoframe.setLayout(isogridbox)
		gridbox.addWidget(isoframe, 2, 0, 2, 1)
		gridbox.addWidget(self.thr, 4, 0, 1, 2)

		# Set to default, but run only once and not in each base class
		if type(self) == EMIsosurfaceInspector: self.updateItemControls()
		self.histogram_widget.setProbe(self.item3d().isothr) # The needs to be node AFTER the data is set


	def onFileBrowse(self):
		""" Find a color map file """
		self.openbrowser = EMBrowserWidget(withmodal=True,multiselect=False)
		QtCore.QObject.connect(self.openbrowser, QtCore.SIGNAL("ok"),self._onopen_ok)
		QtCore.QObject.connect(self.openbrowser, QtCore.SIGNAL("cancel"),self._onopen_cancel)
		self.openbrowser.show()

	def _onopen_ok(self):
		""" load color map file """
		file_path = self.openbrowser.getResult()[0]
		if file_path:
			self.item3d().setCmapData(file_path)
		self.openbrowser.close()
		self.openbrowser=None

	def _onopen_cancel(self):
		""" Never mind....."""
		self.openbrowser.close()
		self.openbrowser=None

	def onColorByMap(self):
		""" Display colors map """
		if self.colorbymap.isChecked():
			self.item3d().setRGBmode(2)
		else:
			self.item3d().setRGBmode(0)

		self.inspector().updateSceneGraph()

	def reColorMapMinMax(self, val):
		self.item3d().setCmapMinMax(self.cmapmin.getValue(), self.cmapmax.getValue())
		self.inspector().updateSceneGraph()

	def onCullFaces(self):
		self.item3d().cullbackfaces = self.cullbackface.isChecked()
		self.inspector().updateSceneGraph()

	def onThresholdSlider(self,val):
		self.item3d().setThreshold(val)
		self.histogram_widget.setProbe(val)
		self.inspector().updateSceneGraph()

	def onHistogram(self, val):

		self.thr.setValue(val, quiet=1)
		self.item3d().setThreshold(val)
		self.inspector().updateSceneGraph()

	def onSampling(self, val):
		self.item3d().setSample(val)
		self.inspector().updateSceneGraph()

	def onWireframe(self):
		self.item3d().wire = self.wireframe.isChecked()
		self.inspector().updateSceneGraph()

	def onColorByRadius(self):
		if self.colorbyradius.isChecked():
			self.item3d().setRGBmode(1)
		else:
			self.item3d().setRGBmode(0)
		self.inspector().updateSceneGraph()

	def reColorScale(self):
		self.item3d().setRGBcolorScaling(self.innercolorscaling.getValue(), self.outercolorscaling.getValue())
		self.inspector().updateSceneGraph()

	def setSamplingRange(self,range):
		self.sampling_spinbox.setMinimum(1)
		self.sampling_spinbox.setMaximum(1+range-1)

class EMIsosurface(EMItem3D,EM3DModel):
	"""
	This displays an isosurface, which is a surface containing all the voxels that have the value of the given threshold.
	It must have an EMDataItem3D as its parent in the tree data structure for the scene graph.
	"""
	name = "Isosurface"
	nodetype = "DataChild"

	@staticmethod
	def getNodeDialogWidget(attribdict):
		"""
		Get Isosurface Widget
		"""
		isosurfacewidget = QtGui.QWidget()
		grid = QtGui.QGridLayout()
		node_name_data_label = QtGui.QLabel("Isosurface Name")
		attribdict["node_name"] = QtGui.QLineEdit(str(EMIsosurface.name))
		grid.addWidget(node_name_data_label, 0, 0, 1, 2)
		grid.addWidget(attribdict["node_name"], 0, 2, 1, 2)
		EMItem3D.get_transformlayout(grid, 2, attribdict)
		isosurfacewidget.setLayout(grid)

		return isosurfacewidget

	@staticmethod
	def getNodeForDialog(attribdict):

		"""
		Create a new node using a attribdict
		"""
		return EMIsosurface(attribdict["parent"], transform=EMItem3D.getTransformFromDict(attribdict))

	def __init__(self, parent=None, children = set(), transform=None):
		"""
		@param parent: should be an EMDataItem3D instance for proper functionality.
		"""
		if not transform: transform = Transform()	# Object initialization should not be put in the constructor. Causes issues
		EMItem3D.__init__(self, parent, children, transform=transform)
		self.data = None
		self.data_copy = None

		self.light= True
		self.isothr = None #Will be set in self.dataChanged()
		self.isodl = 0
		self.smpval=-1
		self.wire = False
		self.cullbackfaces = True
		self.tex_name = 0
		self.texture = False
		self.rgbmode = 0
		self.innerrad = 0.0
		self.outerrad = 0.0
		self.cmapmin = 0.0
		self.cmapmax = 0.0
		self.cmapdata = None
		self.cmapfilename = ""
		self.dxsize=0
		self.dysize=0
		self.dzsize=0

		self.cube=False
		self.cam=Camera2(self)
		self.vdtools = EMViewportDepthTools(self)

#		self.brightness = 0
#		self.contrast = 10
		self.rank = 1
		self.data_copy = None
		self.force_update = False
		self.loadColors()

		# color Needed for inspector to work John Flanagan
		self.diffuse = self.colors[self.isocolor]["diffuse"]
		self.specular = self.colors[self.isocolor]["specular"]
		self.ambient = self.colors[self.isocolor]["ambient"]
		self.shininess = self.colors[self.isocolor]["shininess"]

		if parent: self.dataChanged()

	def getEvalString(self):
		return "EMIsosurface(),getEvalString"

	def dataChanged(self):
		"""
		When the EMData changes for EMDataItem3D parent node, this method is called. It is responsible for updating the state of the slice node.
		"""
		data = self.getParent().getData()

		self.minden = data.get_attr("minimum")
		self.maxden = data.get_attr("maximum")
		self.mean   = data.get_attr("mean")
		self.sigma  = data.get_attr("sigma")

		#This computes initial threshold. Steven Murray does seem to have much success with this though.
		if self.isothr: #there was data previously
			normalized_threshold = (self.isothr - self.mean)/self.sigma
		else:
			normalized_threshold = 3.0

		self.isothr = self.mean+normalized_threshold*self.sigma
		self.histogram_data = data.calc_hist(256,self.minden, self.maxden)

		self.force_update = True
		self.isorender = MarchingCubes(data)
		self.outerrad = data.get_xsize()/2.0

		self.dxsize=data.get_xsize()/2.0
		self.dysize=data.get_ysize()/2.0
		self.dzsize=data.get_zsize()/2.0

		if ( self.texture ): self.gen_texture()
		if self.item_inspector: self.getItemInspector().updateItemControls() # The idea is to use lazy evaluation for the item inspectors. Forcing inspector creation too early causes bugs!


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

	def getItemDictionary(self):
		"""
		Return a dictionary of item parameters (used for restoring sessions
		"""
		dictionary = super(EMIsosurface, self).getItemDictionary()
		dictionary.update({"ISOPARS":[self.wire, self.cullbackfaces, self.isothr, self.rgbmode, self.innerrad, self.outerrad, self.cmapmin, self.cmapmax, self.cmapfilename],"COLOR":[self.ambient, self.diffuse, self.specular, self.shininess]})
		return dictionary

	def setUsingDictionary(self, dictionary):
		"""
		Set item attributes using a dictionary, used in session restoration
		"""
		super(EMIsosurface, self).setUsingDictionary(dictionary)
		self.setAmbientColor(dictionary["COLOR"][0][0], dictionary["COLOR"][0][1], dictionary["COLOR"][0][2], dictionary["COLOR"][0][3])
		self.setDiffuseColor(dictionary["COLOR"][1][0], dictionary["COLOR"][1][1], dictionary["COLOR"][1][2], dictionary["COLOR"][1][3])
		self.setSpecularColor(dictionary["COLOR"][2][0], dictionary["COLOR"][2][1], dictionary["COLOR"][2][2], dictionary["COLOR"][2][3])
		self.setShininess(dictionary["COLOR"][3])
		self.wire = dictionary["ISOPARS"][0]
		self.cullbackfaces = dictionary["ISOPARS"][1]
		self.setThreshold(dictionary["ISOPARS"][2])
		try:
			if dictionary["ISOPARS"][8]:
				self.setCmapData(dictionary["ISOPARS"][8])
				self.setCmapMinMax(dictionary["ISOPARS"][6], dictionary["ISOPARS"][7])
		except:
			pass
		try:
			self.setRGBmode(dictionary["ISOPARS"][3])
			self.setRGBcolorScaling(dictionary["ISOPARS"][4], dictionary["ISOPARS"][5])
		except:
			pass

	def loadColors(self):
		self.colors = get_default_gl_colors()

		self.isocolor = "bluewhite"

	def setSample(self,val):
		if ( self.smpval != int(val)):
			# the minus two is here because the marching cubes thinks -1 is the high level of detail, 0 is the next best and  so forth
			# However the user wants the highest level of detail to be 1, and the next best to be 2 and then 3 etc
			self.smpval = int(val)-2
			self.getIsosurfaceContours()

	def setRGBmode(self, mode):
		""" Set the RGB mode """
		self.rgbmode = mode
		self.isorender.set_rgb_mode(mode)

	def setCmapData(self, data):
		""" Sets the cmap data """
		if isinstance(data, EMData):
			self.cmapdata = data
			if data.has_attr('source_path'):
				self.cmapfilename = data.get_attr('source_path')
		else:
			self.cmapdata = EMData(str(data))
			self.cmapfilename = os.path.abspath(data)

		self.isorender.set_cmap_data(self.cmapdata)
		self.setCmapMinMax(self.cmapdata.get_attr('minimum'), self.cmapdata.get_attr('maximum'))

		# Update inspector
		if self.item_inspector: self.getItemInspector().updateItemControls()

	def setCmapMinMax(self, minimum, maximum):
		""" Sets the cmap min and max """
		self.cmapmin = minimum
		self.cmapmax = maximum
		self.isorender.set_cmap_minmax(minimum, maximum)

	def setRGBcolorScaling(self, inner, outer):
		self.innerrad = inner
		self.outerrad = outer
		self.isorender.set_rgb_scale(inner, outer)

	def getIsosurfaceContours(self):
		# create the isosurface display list
		self.isorender.set_surface_value(self.isothr)
		self.isorender.set_sampling(self.smpval)

		if (float(glGetString(GL_VERSION).split(".")[0])>2):
			GLUtil.contour_isosurface(self.isorender)
		else:
			if ( self.texture ):
				if ( self.tex_name == 0 ):
					self.update_data_and_texture()

			face_z = False
			if self.dzsize <= 2:
				face_z = True

			if ( self.texture  ):
				self.isodl = GLUtil.get_isosurface_dl(self.isorender, self.tex_name,face_z)
			else:
				self.isodl = GLUtil.get_isosurface_dl(self.isorender, 0,face_z)

	def renderNode(self):
		if (not isinstance(self.parent.data,EMData)): return
		#a = time()
		if (float(glGetString(GL_VERSION).split(".")[0])>2):

			scenegraph = self.getRootNode()
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
			if (self.force_update):
				self.getIsosurfaceContours()
				self.force_update = False

			if self.rgbmode != 0:
				glEnable(GL_COLOR_MATERIAL)
				glColorMaterial(GL_FRONT, GL_AMBIENT)
				glColorMaterial(GL_FRONT, GL_DIFFUSE)
			else:
				glDisable(GL_COLOR_MATERIAL)
			# This code draws an outline around the isosurface
			if (self.is_selected or self.getParent().is_selected) and glGetIntegerv(GL_RENDER_MODE) == GL_RENDER and not self.isSelectionHidded(): # No need for outlining in selection mode
				glPushAttrib( GL_ALL_ATTRIB_BITS )

				if scenegraph.camera.getCappingMode(): glDisable(GL_CULL_FACE)
				# First get a stencil of the object silluette
				glClearStencil(0)
				glClear( GL_STENCIL_BUFFER_BIT )
				glEnable( GL_STENCIL_TEST )
				glStencilFunc( GL_ALWAYS, 1, 0xFFFF )		# Write to stencil buffer
				glStencilOp( GL_KEEP, GL_KEEP, GL_REPLACE )	# Make stencil of object outline
				if ( self.wire ):
					glPolygonMode(GL_FRONT_AND_BACK,GL_LINE)
				else:
					glPolygonMode(GL_FRONT_AND_BACK,GL_FILL)
				self.renderIso()

				# Then render the outline
				glStencilFunc( GL_NOTEQUAL, 1, 0xFFFF )		# The object itself is stenciled out
				glLineWidth( 4.0 )				# By increasing the line width only the outline is drawn
				glPolygonMode( GL_FRONT_AND_BACK, GL_LINE )
				glMaterialfv(GL_FRONT, GL_EMISSION, [0.0, 1.0, 0.0, 1.0])
				self.renderIso()

				glPopAttrib()

			elif (scenegraph.camera.getCappingMode() and not scenegraph.zslicemode and glGetIntegerv(GL_RENDER_MODE) == GL_RENDER):
				# First get a stencil of the object silluette
				glPushAttrib( GL_ALL_ATTRIB_BITS )

				glDisable(GL_CULL_FACE)
				glClear( GL_STENCIL_BUFFER_BIT )
				glClearStencil(0)
				glEnable( GL_STENCIL_TEST )
				glStencilFunc( GL_ALWAYS, 1, 0xFFFF )		# Write to stencil buffer
				glStencilOp( GL_KEEP, GL_INVERT, GL_INVERT )	# Stencil buffer is 0 along clipping planes

				self.renderIso()

				glStencilFunc( GL_NOTEQUAL, 0, 0xFFFF )

				glDisable(GL_COLOR_MATERIAL)
				glMaterialfv(GL_FRONT, GL_SPECULAR, [0.0,0.0,0.0,1.0])
				glMaterialfv(GL_FRONT, GL_AMBIENT, scenegraph.camera.getCapColor())
				glMaterialfv(GL_FRONT, GL_DIFFUSE, scenegraph.camera.getCapColor())

				# Draw plane of the capping color
				glPushMatrix()
				glLoadIdentity()

				x = float(scenegraph.camera.getWidth()/2.0)
				y = float(scenegraph.camera.getHeight()/2.0)
				z = -float(scenegraph.camera.getClipNear() + 0.5)

				glBegin(GL_QUADS)
				glVertex3f(-x, -y, z)
				glVertex3f(x, -y, z)
				glVertex3f(x, y, z)
				glVertex3f(-x, y, z)
				glEnd()

				glPopMatrix()

				glPopAttrib()
			else:
				glPushAttrib( GL_ALL_ATTRIB_BITS )

				self.renderIso()

				glPopAttrib()
	#		self.draw_bc_screen() #TODO: check into porting this from EM3DModel

			if cull: glEnable(GL_CULL_FACE)
			else: glDisable(GL_CULL_FACE)

			if ( polygonmode[0] == GL_LINE ):
				glPolygonMode(GL_FRONT, GL_LINE)
			else:
				glPolygonMode(GL_FRONT, GL_FILL)

			#print "total time is", time()-a
		else:
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
				self.getIsosurfaceContours()
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
			glTranslate(-self.dxsize,-self.dysize,-self.dzsize)
			if ( self.texture ):
				glScalef(self.dxsize*2.0,self.dysize*2.0,self.dzsize*2.0)
#			print self.isodl
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

#		print "renderiso"
		GLUtil.render_using_VBOs(self.isorender, 0, 0)

		glPopMatrix()

	def setThreshold(self, val):
		if (self.isothr != val):
			self.isothr = val
#			self.brightness = -val

			self.getIsosurfaceContours()

	def gen_texture(self):
		if ( self.texture == False ): return
		if ( self.tex_name != 0 ):
			glDeleteTextures(self.tex_name)

		if ( self.data_copy == None ):
			self.tex_name = GLUtil.gen_gl_texture(self.data)
		else:
			self.tex_name = GLUtil.gen_gl_texture(self.data_copy)

	def eye_coords_dif(self,x1,y1,x2,y2,mdepth=True):
		return self.vdtools.eye_coords_dif(x1,y1,x2,y2,mdepth)
