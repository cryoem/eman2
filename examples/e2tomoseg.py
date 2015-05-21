#!/usr/bin/env python

#
# Author: James Michael Bell 5/20/2014 (jmbell@bcm.edu)
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
from emapplication import EMApp
from emimage2d import EMImage2DWidget
from emscene3d import EMScene3D, EMInspector3D, EMQTreeWidget
import os
from valslider import ValSlider, EMANToolButton, EMSpinWidget, EMQTColorWidget

from PyQt4 import QtCore
from PyQt4 import QtGui

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <image file> ...
	
	e2tomoseg.py is a simple tomogram segmentation tool with built-in automation.
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	#parser.add_argument("-t","--tomo",type=str,help="<rawptcl>,<classmx> Show particles associated class-averages")
	#parser.add_argument("-s","--seg",type=str,help="A specialized flag that disables auto contrast for the display of particles stacks and 2D images only.")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()
	
	logid=E2init(sys.argv)
	
	app = EMApp()
	
	volume = TomoSegVolumeViewer() # window to see 3D segmentation overlayed on top of current slice
	slice = TomoSegSliceViewer() # window to view current slice of tomogram
	tools = TomoSegInspector() # window to hold tools/annotation tree for segmentation
	
	volume.show() # show current slice in 2D with tomogram segmentation in 3D
	slice.show() # show slice in 2D with segmentation overlay
	tools.show() # show tools for annotation
	
	app.exec_()

	E2end(logid)
		
class TomoSegVolumeViewer(EMImage2DWidget):
	
	def __init__(self):
		super(TomoSegVolumeViewer,self).__init__()

class TomoSegSliceViewer(EMScene3D):
	
	def __init__(self):
		super(TomoSegSliceViewer,self).__init__()

class TomoSegInspector(QtGui.QWidget):
	
	def __init__(self):
		super(TomoSegInspector,self).__init__()
		
		self.mintreewidth = 250		# minimum width of the tree
		self.mincontrolwidth = 0
		
		vbox = QtGui.QVBoxLayout(self)
		self.inspectortab = QtGui.QTabWidget()
		self.inspectortab.addTab(self.getToolsWidget(), "Tools")
		self.inspectortab.addTab(self.getTreeWidget(), "Annotations")
		self.inspectortab.addTab(self.getUtilsWidget(), "Utils")
		vbox.addWidget(self.inspectortab)
		
		self.setLayout(vbox)
		self.updateGeometry()

	def getToolsWidget(self):
		tooltabs = QtGui.QTabWidget()
		tooltabs.addTab(self.getAutomaticTools(), "Automatic")
		tooltabs.addTab(self.getSemiAutomaticTools(), "Interactive")
		tooltabs.addTab(self.getManualTools(), "Manual")
		return tooltabs
	
	def getAutomaticTools(self):
		widget = QtGui.QWidget()
		
# 		font = QtGui.QFont()
# 		font.setBold(True)
# 		toollabel = QtGui.QLabel("Automatic Segmentation")
# 		toollabel.setFont(font)
		
		grid = QtGui.QGridLayout()
		self.rotatetool = EMANToolButton()
# 		self.rotatetool.setIcon(QtGui.QIcon(QtGui.QPixmap(rotateicon)))
 		self.rotatetool.setToolTip("Description")
 		self.rotatetool_label = QtGui.QLabel("Name")
		
		self.translatetool =EMANToolButton()
# 		self.translatetool.setIcon(QtGui.QIcon(QtGui.QPixmap(crosshairsicon)))
 		self.translatetool.setToolTip("Description")
 		self.translatetool_label = QtGui.QLabel("Name")
		
		self.ztranslate = EMANToolButton()
# 		self.ztranslate.setIcon(QtGui.QIcon(QtGui.QPixmap(ztransicon)))
 		self.ztranslate.setToolTip("Description")
 		self.ztranslatetool_label = QtGui.QLabel("Name")
		
		self.scaletool = EMANToolButton()
# 		self.scaletool.setIcon(QtGui.QIcon(QtGui.QPixmap(scaleicon)))
 		self.scaletool.setToolTip("Description")
 		self.scaletool_label = QtGui.QLabel("Name")
		
		self.rulertool = EMANToolButton()
# 		self.rulertool.setIcon(QtGui.QIcon(QtGui.QPixmap(rulericon)))
 		self.rulertool.setToolTip("Description")
 		self.rulertool_label = QtGui.QLabel("Name")
		
		self.selectiontool = EMANToolButton()
# 		self.selectiontool.setIcon(QtGui.QIcon(QtGui.QPixmap(selectionicon)))
 		self.selectiontool.setToolTip("Description")
 		self.selectiontool_label = QtGui.QLabel("Name")
		
		self.multiselectiontool = EMANToolButton()
# 		self.multiselectiontool.setIcon(QtGui.QIcon(QtGui.QPixmap(multiselectoricon)))
 		self.multiselectiontool.setToolTip("Description")
 		self.multiselectiontool_label = QtGui.QLabel("Name")
		
		self.linetool = EMANToolButton()
# 		self.linetool.setIcon(QtGui.QIcon(QtGui.QPixmap(lineicon)))
		self.linetool.setToolTip("Description")
		self.linetool_label = QtGui.QLabel("Name")
		
		self.cubetool = EMANToolButton()
# 		self.cubetool.setIcon(QtGui.QIcon(QtGui.QPixmap(cubeicon)))
 		self.cubetool.setToolTip("Description")
 		self.cubetool_label = QtGui.QLabel("Name")
		
		self.spheretool = EMANToolButton()
# 		self.spheretool.setIcon(QtGui.QIcon(QtGui.QPixmap(sphereicon)))
 		self.spheretool.setToolTip("Description")
 		self.spheretool_label = QtGui.QLabel("Name")
		
		self.cylindertool = EMANToolButton()
# 		self.cylindertool.setIcon(QtGui.QIcon(QtGui.QPixmap(cylindericon)))
 		self.cylindertool.setToolTip("Description")
 		self.cylindertool_label = QtGui.QLabel("Name")
		
		self.conetool = EMANToolButton()
# 		self.conetool.setIcon(QtGui.QIcon(QtGui.QPixmap(coneicon)))
 		self.conetool.setToolTip("Description")
 		self.conetool_label = QtGui.QLabel("Name")
		
		self.texttool = EMANToolButton()
# 		self.texttool.setIcon(QtGui.QIcon(QtGui.QPixmap(texticon)))
 		self.texttool.setToolTip("Description")
 		self.texttool_label = QtGui.QLabel("Name")
 		
		self.datatool = EMANToolButton()
# 		self.datatool.setIcon(QtGui.QIcon(QtGui.QPixmap(dataicon)))
 		self.datatool.setToolTip("Description")
		self.datatool_label = QtGui.QLabel("Name")
		
		self.apptool = EMANToolButton()
#		self.apptool.setIcon(QtGui.QIcon(QtGui.QPixmap(appicon)))
 		self.apptool.setToolTip("Description")
		self.apptool_label = QtGui.QLabel("Name")

		# buttons
		grid.addWidget(self.selectiontool,1,0)
		grid.addWidget(self.multiselectiontool,2,0)
		grid.addWidget(self.translatetool,3,0)
		grid.addWidget(self.ztranslate,4,0)
		grid.addWidget(self.rotatetool,5,0)
		grid.addWidget(self.scaletool,6,0)
		grid.addWidget(self.rulertool,7,0)
		grid.addWidget(self.linetool,8,0)
		grid.addWidget(self.cubetool,9,0)
		grid.addWidget(self.spheretool,10,0)
		grid.addWidget(self.cylindertool,11,0)
		grid.addWidget(self.conetool,12,0)
		grid.addWidget(self.texttool,13,0)
		grid.addWidget(self.datatool,14,0)
		grid.addWidget(self.apptool,15,0)
		grid.setAlignment(QtCore.Qt.AlignLeft)
		# labels
		#grid.addWidget(toollabel,0,0)
		grid.addWidget(self.selectiontool_label,1,1)
		grid.addWidget(self.multiselectiontool_label,2,1)
		grid.addWidget(self.translatetool_label,3,1)
		grid.addWidget(self.ztranslatetool_label,4,1)
		grid.addWidget(self.rotatetool_label,5,1)
		grid.addWidget(self.scaletool_label,6,1)
		grid.addWidget(self.rulertool_label,7,1)
		grid.addWidget(self.linetool_label,8,1)
		grid.addWidget(self.cubetool_label,9,1)
		grid.addWidget(self.spheretool_label,10,1)
		grid.addWidget(self.cylindertool_label,11,1)
		grid.addWidget(self.conetool_label,12,1)
		grid.addWidget(self.texttool_label,13,1)
		grid.addWidget(self.datatool_label,14,1)
		grid.addWidget(self.apptool_label,15,1)
		grid.setAlignment(QtCore.Qt.AlignLeft)
				
		widget.setLayout(grid)
		return widget

	def getSemiAutomaticTools(self):
		widget = QtGui.QWidget()
		grid = QtGui.QGridLayout()
		
# 		font = QtGui.QFont()
# 		font.setBold(True)
# 		toollabel = QtGui.QLabel("Semi-Automatic Segmentation")
# 		toollabel.setFont(font)
		
		self.rotatetool = EMANToolButton()
# 		self.rotatetool.setIcon(QtGui.QIcon(QtGui.QPixmap(rotateicon)))
 		self.rotatetool.setToolTip("Description")
 		self.rotatetool_label = QtGui.QLabel("Name")
		
		self.translatetool =EMANToolButton()
# 		self.translatetool.setIcon(QtGui.QIcon(QtGui.QPixmap(crosshairsicon)))
 		self.translatetool.setToolTip("Description")
 		self.translatetool_label = QtGui.QLabel("Name")
		
		self.ztranslate = EMANToolButton()
# 		self.ztranslate.setIcon(QtGui.QIcon(QtGui.QPixmap(ztransicon)))
 		self.ztranslate.setToolTip("Description")
 		self.ztranslatetool_label = QtGui.QLabel("Name")
		
		self.scaletool = EMANToolButton()
# 		self.scaletool.setIcon(QtGui.QIcon(QtGui.QPixmap(scaleicon)))
 		self.scaletool.setToolTip("Description")
 		self.scaletool_label = QtGui.QLabel("Name")
		
		self.rulertool = EMANToolButton()
# 		self.rulertool.setIcon(QtGui.QIcon(QtGui.QPixmap(rulericon)))
 		self.rulertool.setToolTip("Description")
 		self.rulertool_label = QtGui.QLabel("Name")
		
		self.selectiontool = EMANToolButton()
# 		self.selectiontool.setIcon(QtGui.QIcon(QtGui.QPixmap(selectionicon)))
 		self.selectiontool.setToolTip("Description")
 		self.selectiontool_label = QtGui.QLabel("Name")
		
		self.multiselectiontool = EMANToolButton()
# 		self.multiselectiontool.setIcon(QtGui.QIcon(QtGui.QPixmap(multiselectoricon)))
 		self.multiselectiontool.setToolTip("Description")
 		self.multiselectiontool_label = QtGui.QLabel("Name")
		
		self.linetool = EMANToolButton()
# 		self.linetool.setIcon(QtGui.QIcon(QtGui.QPixmap(lineicon)))
		self.linetool.setToolTip("Description")
		self.linetool_label = QtGui.QLabel("Name")
		
		self.cubetool = EMANToolButton()
# 		self.cubetool.setIcon(QtGui.QIcon(QtGui.QPixmap(cubeicon)))
 		self.cubetool.setToolTip("Description")
 		self.cubetool_label = QtGui.QLabel("Name")
		
		self.spheretool = EMANToolButton()
# 		self.spheretool.setIcon(QtGui.QIcon(QtGui.QPixmap(sphereicon)))
 		self.spheretool.setToolTip("Description")
 		self.spheretool_label = QtGui.QLabel("Name")
		
		self.cylindertool = EMANToolButton()
# 		self.cylindertool.setIcon(QtGui.QIcon(QtGui.QPixmap(cylindericon)))
 		self.cylindertool.setToolTip("Description")
 		self.cylindertool_label = QtGui.QLabel("Name")
		
		self.conetool = EMANToolButton()
# 		self.conetool.setIcon(QtGui.QIcon(QtGui.QPixmap(coneicon)))
 		self.conetool.setToolTip("Description")
 		self.conetool_label = QtGui.QLabel("Name")
		
		self.texttool = EMANToolButton()
# 		self.texttool.setIcon(QtGui.QIcon(QtGui.QPixmap(texticon)))
 		self.texttool.setToolTip("Description")
 		self.texttool_label = QtGui.QLabel("Name")
 		
		self.datatool = EMANToolButton()
# 		self.datatool.setIcon(QtGui.QIcon(QtGui.QPixmap(dataicon)))
 		self.datatool.setToolTip("Description")
		self.datatool_label = QtGui.QLabel("Name")
		
		self.apptool = EMANToolButton()
#		self.apptool.setIcon(QtGui.QIcon(QtGui.QPixmap(appicon)))
 		self.apptool.setToolTip("Description")
		self.apptool_label = QtGui.QLabel("Name")

		# buttons
		grid.addWidget(self.selectiontool,1,0)
		grid.addWidget(self.multiselectiontool,2,0)
		grid.addWidget(self.translatetool,3,0)
		grid.addWidget(self.ztranslate,4,0)
		grid.addWidget(self.rotatetool,5,0)
		grid.addWidget(self.scaletool,6,0)
		grid.addWidget(self.rulertool,7,0)
		grid.addWidget(self.linetool,8,0)
		grid.addWidget(self.cubetool,9,0)
		grid.addWidget(self.spheretool,10,0)
		grid.addWidget(self.cylindertool,11,0)
		grid.addWidget(self.conetool,12,0)
		grid.addWidget(self.texttool,13,0)
		grid.addWidget(self.datatool,14,0)
		grid.addWidget(self.apptool,15,0)
		grid.setAlignment(QtCore.Qt.AlignLeft)
		# labels
		#grid.addWidget(toollabel,0,0)
		grid.addWidget(self.selectiontool_label,1,1)
		grid.addWidget(self.multiselectiontool_label,2,1)
		grid.addWidget(self.translatetool_label,3,1)
		grid.addWidget(self.ztranslatetool_label,4,1)
		grid.addWidget(self.rotatetool_label,5,1)
		grid.addWidget(self.scaletool_label,6,1)
		grid.addWidget(self.rulertool_label,7,1)
		grid.addWidget(self.linetool_label,8,1)
		grid.addWidget(self.cubetool_label,9,1)
		grid.addWidget(self.spheretool_label,10,1)
		grid.addWidget(self.cylindertool_label,11,1)
		grid.addWidget(self.conetool_label,12,1)
		grid.addWidget(self.texttool_label,13,1)
		grid.addWidget(self.datatool_label,14,1)
		grid.addWidget(self.apptool_label,15,1)
		grid.setAlignment(QtCore.Qt.AlignLeft)
		
		widget.setLayout(grid)
		return widget

	def getManualTools(self):
		widget = QtGui.QWidget()
		grid = QtGui.QGridLayout()
		
# 		font = QtGui.QFont()
# 		font.setBold(True)
# 		toollabel = QtGui.QLabel("Manual Segmentation")
# 		toollabel.setFont(font)
		
		self.rotatetool = EMANToolButton()
# 		self.rotatetool.setIcon(QtGui.QIcon(QtGui.QPixmap(rotateicon)))
 		self.rotatetool.setToolTip("Description")
 		self.rotatetool_label = QtGui.QLabel("Name")
		
		self.translatetool =EMANToolButton()
# 		self.translatetool.setIcon(QtGui.QIcon(QtGui.QPixmap(crosshairsicon)))
 		self.translatetool.setToolTip("Description")
 		self.translatetool_label = QtGui.QLabel("Name")
		
		self.ztranslate = EMANToolButton()
# 		self.ztranslate.setIcon(QtGui.QIcon(QtGui.QPixmap(ztransicon)))
 		self.ztranslate.setToolTip("Description")
 		self.ztranslatetool_label = QtGui.QLabel("Name")
		
		self.scaletool = EMANToolButton()
# 		self.scaletool.setIcon(QtGui.QIcon(QtGui.QPixmap(scaleicon)))
 		self.scaletool.setToolTip("Description")
 		self.scaletool_label = QtGui.QLabel("Name")
		
		self.rulertool = EMANToolButton()
# 		self.rulertool.setIcon(QtGui.QIcon(QtGui.QPixmap(rulericon)))
 		self.rulertool.setToolTip("Description")
 		self.rulertool_label = QtGui.QLabel("Name")
		
		self.selectiontool = EMANToolButton()
# 		self.selectiontool.setIcon(QtGui.QIcon(QtGui.QPixmap(selectionicon)))
 		self.selectiontool.setToolTip("Description")
 		self.selectiontool_label = QtGui.QLabel("Name")
		
		self.multiselectiontool = EMANToolButton()
# 		self.multiselectiontool.setIcon(QtGui.QIcon(QtGui.QPixmap(multiselectoricon)))
 		self.multiselectiontool.setToolTip("Description")
 		self.multiselectiontool_label = QtGui.QLabel("Name")
		
		self.linetool = EMANToolButton()
# 		self.linetool.setIcon(QtGui.QIcon(QtGui.QPixmap(lineicon)))
		self.linetool.setToolTip("Description")
		self.linetool_label = QtGui.QLabel("Name")
		
		self.cubetool = EMANToolButton()
# 		self.cubetool.setIcon(QtGui.QIcon(QtGui.QPixmap(cubeicon)))
 		self.cubetool.setToolTip("Description")
 		self.cubetool_label = QtGui.QLabel("Name")
		
		self.spheretool = EMANToolButton()
# 		self.spheretool.setIcon(QtGui.QIcon(QtGui.QPixmap(sphereicon)))
 		self.spheretool.setToolTip("Description")
 		self.spheretool_label = QtGui.QLabel("Name")
		
		self.cylindertool = EMANToolButton()
# 		self.cylindertool.setIcon(QtGui.QIcon(QtGui.QPixmap(cylindericon)))
 		self.cylindertool.setToolTip("Description")
 		self.cylindertool_label = QtGui.QLabel("Name")
		
		self.conetool = EMANToolButton()
# 		self.conetool.setIcon(QtGui.QIcon(QtGui.QPixmap(coneicon)))
 		self.conetool.setToolTip("Description")
 		self.conetool_label = QtGui.QLabel("Name")
		
		self.texttool = EMANToolButton()
# 		self.texttool.setIcon(QtGui.QIcon(QtGui.QPixmap(texticon)))
 		self.texttool.setToolTip("Description")
 		self.texttool_label = QtGui.QLabel("Name")
 		
		self.datatool = EMANToolButton()
# 		self.datatool.setIcon(QtGui.QIcon(QtGui.QPixmap(dataicon)))
 		self.datatool.setToolTip("Description")
		self.datatool_label = QtGui.QLabel("Name")
		
		self.apptool = EMANToolButton()
#		self.apptool.setIcon(QtGui.QIcon(QtGui.QPixmap(appicon)))
 		self.apptool.setToolTip("Description")
		self.apptool_label = QtGui.QLabel("Name")

		# buttons
		grid.addWidget(self.selectiontool,1,0)
		grid.addWidget(self.multiselectiontool,2,0)
		grid.addWidget(self.translatetool,3,0)
		grid.addWidget(self.ztranslate,4,0)
		grid.addWidget(self.rotatetool,5,0)
		grid.addWidget(self.scaletool,6,0)
		grid.addWidget(self.rulertool,7,0)
		grid.addWidget(self.linetool,8,0)
		grid.addWidget(self.cubetool,9,0)
		grid.addWidget(self.spheretool,10,0)
		grid.addWidget(self.cylindertool,11,0)
		grid.addWidget(self.conetool,12,0)
		grid.addWidget(self.texttool,13,0)
		grid.addWidget(self.datatool,14,0)
		grid.addWidget(self.apptool,15,0)
		grid.setAlignment(QtCore.Qt.AlignLeft)
		# labels
		#grid.addWidget(toollabel,0,0)
		grid.addWidget(self.selectiontool_label,1,1)
		grid.addWidget(self.multiselectiontool_label,2,1)
		grid.addWidget(self.translatetool_label,3,1)
		grid.addWidget(self.ztranslatetool_label,4,1)
		grid.addWidget(self.rotatetool_label,5,1)
		grid.addWidget(self.scaletool_label,6,1)
		grid.addWidget(self.rulertool_label,7,1)
		grid.addWidget(self.linetool_label,8,1)
		grid.addWidget(self.cubetool_label,9,1)
		grid.addWidget(self.spheretool_label,10,1)
		grid.addWidget(self.cylindertool_label,11,1)
		grid.addWidget(self.conetool_label,12,1)
		grid.addWidget(self.texttool_label,13,1)
		grid.addWidget(self.datatool_label,14,1)
		grid.addWidget(self.apptool_label,15,1)
		grid.setAlignment(QtCore.Qt.AlignLeft)
				
		widget.setLayout(grid)
		return widget
	
	def getTreeWidget(self):
		"""
		This returns the treeview-control panel widget
		"""
		widget = QtGui.QWidget()
		hbox = QtGui.QHBoxLayout(widget)
		treeframe = QtGui.QFrame()
		treeframe.setFrameShape(QtGui.QFrame.StyledPanel)
		treeframe.setLayout(self._get_tree_layout(widget))
		treeframe.setMinimumWidth(self.mintreewidth)
		hbox.addWidget(treeframe)
		self.stacked_widget = QtGui.QStackedWidget()
		self.stacked_widget.setFrameShape(QtGui.QFrame.StyledPanel)
		hbox.addWidget(self.stacked_widget)
		widget.setLayout(hbox)
		
		return widget
	
	def _get_tree_layout(self, parent):
		"""
		Returns the tree layout
		"""
		tvbox = QtGui.QVBoxLayout()
		self.tree_widget = EMQTreeWidget(parent)
		self.tree_widget.setHeaderLabel("Choose a item")
		tvbox.addWidget(self.tree_widget)
		self.tree_node_button_add = QtGui.QPushButton("Add Object")
		self.tree_node_button_remove = QtGui.QPushButton("Remove Object")
		self.tree_node_slider = ValSlider(label="Seq:")
		self.tree_node_slider.setIntonly(True)
		self.tree_node_slider.setRange(0,1)
		self.tree_node_slider.setValue(0)
		tvbox.addWidget(self.tree_node_button_add)
		tvbox.addWidget(self.tree_node_button_remove)
		tvbox.addWidget(self.tree_node_slider)
		
# 		QtCore.QObject.connect(self.tree_widget, QtCore.SIGNAL("itemClicked(QTreeWidgetItem*,int)"), self._tree_widget_click)
# 		QtCore.QObject.connect(self.tree_widget, QtCore.SIGNAL("visibleItem(QTreeWidgetItem*)"), self._tree_widget_visible)
# 		QtCore.QObject.connect(self.tree_widget, QtCore.SIGNAL("editItem(QTreeWidgetItem*)"), self._tree_widget_edit)
# 		QtCore.QObject.connect(self.tree_node_button_remove, QtCore.SIGNAL("clicked()"), self._tree_widget_remove)
# 		QtCore.QObject.connect(self.tree_node_button_add, QtCore.SIGNAL("clicked()"), self._on_add_button)
# 		QtCore.QObject.connect(self.tree_node_slider, QtCore.SIGNAL("valueChanged"), self._slider_change)
		
		return tvbox

	
	def _recursiveupdatetreeselvis(self, item):
		item.setSelectionStateBox()
		item.getVisibleState()
		for childidx in xrange(item.childCount()):
			self._recursiveupdatetreeselvis(item.child(childidx))
			
	def updateTreeSelVis(self, selecteditem=None):
		"""
		Update the selection and visibility states. Makes the Sel Vis states an observer of the SG
		"""
		# Update the tree
		self._recursiveupdatetreeselvis(self.tree_widget.topLevelItem(0))
		# Set the desired item if desired
		if selecteditem:
			try:
				self.stacked_widget.setCurrentWidget(selecteditem.getItemInspector())
				self.tree_widget.setCurrentItem(selecteditem.EMQTreeWidgetItem)
				#if selecteditem: self.scenegraph().setCurrentSelection(selecteditem)
			except:
				pass
			# Unsure unqiue selection
			self.ensureUniqueTreeLevelSelection(selecteditem)
			
	def ensureUniqueTreeLevelSelection(self, item):
		"""
		Make sure that we don't select both an ancestor and child at the same time
		"""
		for ancestor in item.getSelectedAncestorNodes():
			if ancestor.EMQTreeWidgetItem:			# Not al ancestors are listed on the inspector tree (such as a data node)
				ancestor.EMQTreeWidgetItem.setSelectionState(False)
		for child in item.getAllSelectedNodes()[1:]: 	# Lop the node itself off
			child.EMQTreeWidgetItem.setSelectionState(False)
	
	def _slider_change(self):
		mdl=int(self.tree_node_slider.getValue())
# 		for i,child in enumerate(self.scenegraph().getChildren()):
# 			if i==mdl: child.setVisibleItem(True)
# 			else: child.setVisibleItem(False)
		
		#self.scenegraph().updateSG()

	def _recursiveAdd(self, parentitem, parentnode,depth=0):
		"""
		Helper function to laod the SG
		"""
		for child in parentnode.getChildren():
			if not child.getLabel(): child.setLabel(child.name)
			addeditem = self.addTreeNode(child.getLabel(), child, parentitem)
			self._recursiveAdd(addeditem, child,depth+1)
		# Expand the data items
		if parentitem.childCount() > 0: parentitem.setExpanded(True)
		self.tree_node_slider.setRange(0,len(parentnode.getChildren())-1)

		
	def loadSG(self):
		"""
		Load the SG
		"""
		pass
		#rootitem = self.addTreeNode("All Objects", self.scenegraph())
		#self._recursiveAdd(rootitem, self.scenegraph())
		
	def addTreeNode(self, name, item3d, parentitem=None, insertionindex=-1):
		"""
		Add a node (item3d) to the TreeWidget if not parent node, otherwise add a child to parent node
		We need to get a GUI for the treeitem. The treeitem and the GUI need know each other so they can talk
		The Treeitem also needs to know the node, so it can talk to the node.
		You can think of this as a three way conversation (the alterative it to use a mediator, but that is not worth it w/ only three players)
		"""
		tree_item = EMQTreeWidgetItem(QtCore.QStringList(name), item3d, parentitem)	# Make a QTreeItem widget, and let the TreeItem talk to the scenegraph node and its GUI
		item3d.setEMQTreeWidgetItem(tree_item)				# Reference to the EMQTreeWidgetItem
		item_inspector = item3d.getItemInspector()				# Get the node GUI controls 
		#return tree_item
		item_inspector.setInspector(self)					# Associate the item GUI with the inspector
		self.stacked_widget.addWidget(item_inspector)			# Add a widget to the stack
		item3d.setLabel(name)						# Set the label
		# Set icon status
		tree_item.setSelectionStateBox()
		# Set parent if one exists	
		if not parentitem:
			self.tree_widget.insertTopLevelItem(0, tree_item)
		else:
			if insertionindex >= 0:
				parentitem.insertChild(insertionindex, tree_item)
			else:
				parentitem.addChild(tree_item)
		return tree_item
	
	def removeTreeNode(self, parentitem, childindex):
		# I am using the parent item rather than the item itself b/c the stupid widget has no , remove self function...
		# Remove both the QTreeWidgetItem and the widget from the WidgetStack, otherwise we'll get memory leaks 
		if parentitem.child(childindex).item3d():
			self.stacked_widget.removeWidget(parentitem.child(childindex).item3d().getItemInspector())
		parentitem.takeChild(childindex)
	
	def clearTree(self):
		"""
		Clear the entire tree
		"""
		if self.tree_widget.topLevelItem(0):
			self.tree_widget.topLevelItem(0).removeAllChildren(self)
			self.tree_widget.takeTopLevelItem(0)
		
	def _tree_widget_click(self, item, col, quiet=False):
		"""
		When a user clicks on the selection tree check box
		"""
		self.stacked_widget.setCurrentWidget(item.item3d().getItemInspector())
		item.setSelectionState(item.checkState(0))
		# This code is to prevent both decendents and childer from being selected....
		if item.checkState(0) == QtCore.Qt.Checked: self.ensureUniqueTreeLevelSelection(item.item3d())
		if not item.item3d().isSelectedItem(): item.item3d().getItemInspector().updateItemControls() # This is too update a widget, translation and rotation may change in parent nodes change
		#self.scenegraph().setCurrentSelection(item.item3d())
		#if not quiet: self.updateSceneGraph()
		
	def _tree_widget_visible(self, item):
		"""
		When a user clicks on the visible icon
		"""
		item.toggleVisibleState()
	
	def _tree_widget_edit(self):
		"""
		When a use middle clicks
		"""
		nodedialog = NodeEditDialog(self, self.tree_widget.currentItem())
		nodedialog.exec_()
		self.activateWindow()
	
	def _on_add_button(self):
		nodedialog =  NodeDialog(self, self.tree_widget.currentItem())
		nodedialog.exec_()
		self.activateWindow()
		
	def _tree_widget_remove(self):
		"""
		When a use wants to remove a node_name
		"""
		item = self.tree_widget.currentItem()
		if item.parent:
			self.removeTreeNode(item.parent(), item.parent().indexOfChild(item)) 
			item.parent().item3d().removeChild(item.item3d())
			# In case we delete the currently selected item, we want6 to move selected item to last selection
# 			if self.scenegraph().getCurrentSelection() == item.item3d():
# 				self.scenegraph().setCurrentSelection(self.tree_widget.currentItem().item3d())
# 			self.updateSceneGraph()
		else:
			print "Error cannot remove root node!!"

	def getUtilsWidget(self):
		"""
		Return the utilites widget
		"""
		uwidget = QtGui.QWidget()
		uvbox = QtGui.QVBoxLayout()
		font = QtGui.QFont()
		font.setBold(True)
		
		self.opensession_button = QtGui.QPushButton("Open Session")
		self.savesession_button = QtGui.QPushButton("Save Session")
		self.savebutton = QtGui.QPushButton("Save Image Snapshot")
		
		self.open_tomogram_button = QtGui.QPushButton("Open Tomogram")
		self.open_segmentation_button = QtGui.QPushButton("Open Segmentation")
		self.save_segmentation_button = QtGui.QPushButton("Save Segmentation")
		
		uvbox.addWidget(self.opensession_button)
		uvbox.addWidget(self.savesession_button)
		uvbox.addWidget(self.savebutton)
		
		uvbox.addWidget(self.open_tomogram_button)
		uvbox.addWidget(self.open_segmentation_button)
		uvbox.addWidget(self.save_segmentation_button)
		uwidget.setLayout(uvbox)
		
		QtCore.QObject.connect(self.savebutton, QtCore.SIGNAL("clicked()"),self._on_save)
		QtCore.QObject.connect(self.savesession_button, QtCore.SIGNAL("clicked()"),self._on_save_session)
		QtCore.QObject.connect(self.opensession_button, QtCore.SIGNAL("clicked()"),self._on_open_session)
		
		QtCore.QObject.connect(self.open_tomogram_button, QtCore.SIGNAL("clicked()"),self._on_open_tomogram)
		QtCore.QObject.connect(self.open_segmentation_button, QtCore.SIGNAL("clicked()"),self._on_open_segmentation)
		QtCore.QObject.connect(self.save_segmentation_button, QtCore.SIGNAL("clicked()"),self._on_save_segmentation)
		
		return uwidget

	def _on_open_session(self):
		"""
		Open a session... (might want to add a warning dialog that this will close the current session)
		"""
		filename = QtGui.QFileDialog.getOpenFileName(self, 'Open Session', os.getcwd(), "*.eman")
		
	def _on_save_session(self):
		"""
		Return a list of all the child items (actually a tree of sorts)
		"""
		filename = QtGui.QFileDialog.getSaveFileName(self, 'Save Session', os.getcwd(), "*.eman")

	def _on_save(self):
		"""
		Save a snapshot of the scene
		"""
		filename = QtGui.QFileDialog.getSaveFileName(self, 'Save Image', os.getcwd(), "(*.tiff *.jpeg *.png)")

	def _on_open_tomogram(self):
		"""
		Open a session
		"""
		filename = QtGui.QFileDialog.getOpenFileName(self, 'Open Tomogram', os.getcwd(), "*.hdf,*.mrc")
		#if filename:
		#	self.scenegraph().loadSession(filename)
		
	def _on_open_segmentation(self):
		"""
		Open a session
		"""
		# Open the file
		filename = QtGui.QFileDialog.getOpenFileName(self, 'Open Segmentation', os.getcwd(), "*.hdf,*.xml")
		#if filename:
		#	self.scenegraph().loadSession(filename)
		
	def _on_save_segmentation(self):
		"""
		Save a snapshot of the scene
		"""
		filename = QtGui.QFileDialog.getSaveFileName(self, 'Save Segmentation', os.getcwd(), "(*.hdf,*.xml)")
		#if filename: # if we cancel
		#	self.scenegraph().saveSnapShot(filename)
		
	def updateInspector(self):
		"""
		Update Inspector,is called whenever the scence changes
		"""
		pass
	
	def updateTree(self, currentnode=None):
		"""
		Update the SG tree
		"""
		self.clearTree()
		self.loadSG()
		# either set the current node to the argument or set it to the one PM has selected
		if currentnode:
			self.tree_widget.setCurrentItem(currentnode.EMQTreeWidgetItem)
			idx = self.stacked_widget.indexOf(currentnode.getItemInspector())
			if idx >= 0: self.stacked_widget.setCurrentIndex(idx)
			self.scenegraph().setCurrentSelection(currentnode)
		else:
			node = self.scenegraph().getCurrentSelection()
			self.tree_widget.setCurrentItem(node.EMQTreeWidgetItem)
			idx = self.stacked_widget.indexOf(node.getItemInspector())
			if idx >= 0: self.stacked_widget.setCurrentIndex(idx)
			
	def updateSceneGraph(self):
		""" 
		Updates SG, in the near future this will be improved to allow for slow operations
		"""
		pass

if __name__ == '__main__':
	main()
