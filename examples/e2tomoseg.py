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
	
	volume = TomoSegVolViewer() # window to see 3D segmentation overlayed on top of current slice
	slice = TomoSegSliceViewer() # window to view current slice of tomogram
	tools = TomoSegInspector() # window to hold tools/annotation tree for segmentation
	
	volume.show() # show current slice in 2D with tomogram segmentation in 3D
	slice.show() # show slice in 2D with segmentation overlay
	tools.show() # show tools for annotation
	
	app.exec_()

	E2end(logid)
		
class TomoSegVolViewer(EMImage2DWidget):
	
	def __init__(self):
		super(TomoSegVolViewer,self).__init__()

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
		self.inspectortab.addTab(self.getTreeWidget(), "Labels")
		self.inspectortab.addTab(self.getUtilsWidget(), "Utils")
		#self.inspectortab.addTab(self.getFileWidget(), "File")
		vbox.addWidget(self.inspectortab)
			
		self.setLayout(vbox)
		self.updateGeometry()

	def getToolsWidget(self):
		widget = QtGui.QWidget()
		tvbox = QtGui.QVBoxLayout()
		font = QtGui.QFont()
		font.setBold(True)
		toollabel = QtGui.QLabel("Tools")
		toollabel.setFont(font)
		self.rotatetool = EMANToolButton()
		self.rotatetool.setIcon(QtGui.QIcon(QtGui.QPixmap(rotateicon)))
		self.rotatetool.setToolTip("Rotate X/Y\nMouse: Right 'n' drag\nHot Key: R")
		self.translatetool =EMANToolButton()
		self.translatetool.setIcon(QtGui.QIcon(QtGui.QPixmap(crosshairsicon)))
		self.translatetool.setToolTip("Translate X/Y\nMouse: Left 'n' drag\nHot Key: T")
		self.ztranslate = EMANToolButton()
		self.ztranslate.setIcon(QtGui.QIcon(QtGui.QPixmap(ztransicon)))
		self.ztranslate.setToolTip("Translate Z\nHot Key: Z")
		self.scaletool = EMANToolButton()
		self.scaletool.setIcon(QtGui.QIcon(QtGui.QPixmap(scaleicon)))
		self.scaletool.setToolTip("Scale\nHot Key: S")
		self.rulertool = EMANToolButton()
		self.rulertool.setIcon(QtGui.QIcon(QtGui.QPixmap(rulericon)))
		self.rulertool.setToolTip("Ruler\nHot Key: U\nDoes NOT account for scaling")
		self.selectiontool = EMANToolButton()
		self.selectiontool.setIcon(QtGui.QIcon(QtGui.QPixmap(selectionicon)))
		self.selectiontool.setToolTip("Select objects\nMouse: Left 'n' drag\nMultiple = + Shift\nHot Key: Esc")
		self.multiselectiontool = EMANToolButton()
		self.multiselectiontool.setIcon(QtGui.QIcon(QtGui.QPixmap(multiselectoricon)))
		self.multiselectiontool.setToolTip("Select multiple objects\nMouse: Left 'n' drag\nHot Key: M")
		self.linetool = EMANToolButton()
		self.linetool.setIcon(QtGui.QIcon(QtGui.QPixmap(lineicon)))
		self.linetool.setToolTip("Insert Line\nHot Key: L")
		self.cubetool = EMANToolButton()
		self.cubetool.setIcon(QtGui.QIcon(QtGui.QPixmap(cubeicon)))
		self.cubetool.setToolTip("Insert Cube\nHot Key: C")
		self.spheretool = EMANToolButton()
		self.spheretool.setIcon(QtGui.QIcon(QtGui.QPixmap(sphereicon)))
		self.spheretool.setToolTip("Insert Sphere\nHot Key: P")
		self.cylindertool = EMANToolButton()
		self.cylindertool.setIcon(QtGui.QIcon(QtGui.QPixmap(cylindericon)))
		self.cylindertool.setToolTip("Insert Cylinder\nHot Key: Y")
		self.conetool = EMANToolButton()
		self.conetool.setIcon(QtGui.QIcon(QtGui.QPixmap(coneicon)))
		self.conetool.setToolTip("Insert Cone\nHot Key: O")
		self.texttool = EMANToolButton()
		self.texttool.setIcon(QtGui.QIcon(QtGui.QPixmap(texticon)))
		self.texttool.setToolTip("Insert Text\nHot Key: X")
		self.datatool = EMANToolButton()
		self.datatool.setIcon(QtGui.QIcon(QtGui.QPixmap(dataicon)))
		self.datatool.setToolTip("Insert Data\nHot Key: E")
		self.apptool = EMANToolButton()
		self.apptool.setIcon(QtGui.QIcon(QtGui.QPixmap(appicon)))
		self.apptool.setToolTip("Application dependent\nHot Key: A")
		
		tvbox.addWidget(toollabel)
		tvbox.addWidget(self.selectiontool)
		tvbox.addWidget(self.multiselectiontool)
		tvbox.addWidget(self.translatetool)
		tvbox.addWidget(self.ztranslate)
		tvbox.addWidget(self.rotatetool)
		tvbox.addWidget(self.scaletool)
		tvbox.addWidget(self.rulertool)
		tvbox.addWidget(self.linetool)
		tvbox.addWidget(self.cubetool)
		tvbox.addWidget(self.spheretool)
		tvbox.addWidget(self.cylindertool)
		tvbox.addWidget(self.conetool)
		tvbox.addWidget(self.texttool)
		tvbox.addWidget(self.datatool)
		tvbox.addWidget(self.apptool)
		tvbox.setAlignment(QtCore.Qt.AlignLeft)
		
		widget.setLayout(tvbox)
		return widget
	
	def _rotatetool_clicked(self, state):
		self.scenegraph().setMouseMode("rotate")
		
	def _transtool_clicked(self, state):
		self.scenegraph().setMouseMode("xytranslate")
		
	def _ztranstool_clicked(self, state):
		self.scenegraph().setMouseMode("ztranslate")
		
	def _scaletool_clicked(self, state):
		self.scenegraph().setMouseMode("scale")
	
	def _rulertool_clicked(self, state):
		self.scenegraph().setMouseMode("ruler")
		
	def _seltool_clicked(self, state):
		self.scenegraph().setMouseMode("selection")
		
	def _multiseltool_clicked(self, state):
		self.scenegraph().setMouseMode("multiselection")
	
	def _linetool_clicked(self, state):
		self.scenegraph().setMouseMode("line")
		
	def _cubetool_clicked(self, state):
		self.scenegraph().setMouseMode("cube")
		
	def _spheretool_clicked(self, state):
		self.scenegraph().setMouseMode("sphere")
		
	def _cylindertool_clicked(self, state):
		self.scenegraph().setMouseMode("cylinder")
	
	def _conetool_clicked(self, state):
		self.scenegraph().setMouseMode("cone")
		
	def _texttool_clicked(self, state):
		self.scenegraph().setMouseMode("text")
		
	def _datatool_clicked(self, state):
		self.scenegraph().setMouseMode("data")
				
	def _apptool_clicked(self, state):
		self.scenegraph().setMouseMode("app")
	
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
		
		QtCore.QObject.connect(self.tree_widget, QtCore.SIGNAL("itemClicked(QTreeWidgetItem*,int)"), self._tree_widget_click)
		QtCore.QObject.connect(self.tree_widget, QtCore.SIGNAL("visibleItem(QTreeWidgetItem*)"), self._tree_widget_visible)
		QtCore.QObject.connect(self.tree_widget, QtCore.SIGNAL("editItem(QTreeWidgetItem*)"), self._tree_widget_edit)
		QtCore.QObject.connect(self.tree_node_button_remove, QtCore.SIGNAL("clicked()"), self._tree_widget_remove)
		QtCore.QObject.connect(self.tree_node_button_add, QtCore.SIGNAL("clicked()"), self._on_add_button)
		QtCore.QObject.connect(self.tree_node_slider, QtCore.SIGNAL("valueChanged"), self._slider_change)
		
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
				if selecteditem: self.scenegraph().setCurrentSelection(selecteditem)
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
		for i,child in enumerate(self.scenegraph().getChildren()):
			if i==mdl: child.setVisibleItem(True)
			else: child.setVisibleItem(False)
		
		self.scenegraph().updateSG()

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
		rootitem = self.addTreeNode("All Objects", self.scenegraph())
		self._recursiveAdd(rootitem, self.scenegraph())
		
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
		self.scenegraph().setCurrentSelection(item.item3d())
		if not quiet: self.updateSceneGraph()
		
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
			if self.scenegraph().getCurrentSelection() == item.item3d():
				self.scenegraph().setCurrentSelection(self.tree_widget.currentItem().item3d())
			self.updateSceneGraph()
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
		# Controls frame
		frame = QtGui.QFrame()
		frame.setFrameShape(QtGui.QFrame.StyledPanel)
		gridbox = QtGui.QGridLayout()
		backgroundcolor_label = QtGui.QLabel("Background Color", frame)
		backgroundcolor_label.setFont(font)
		self.backgroundcolor = EMQTColorWidget(parent=frame)
		self.hideselectionbutton = QtGui.QCheckBox("Hide Display Selections")
		self.hideselectionbutton.setMinimumHeight(100)
		self.hideselectionbutton.setFont(font)
		gridbox.addWidget(backgroundcolor_label, 0, 0)
		gridbox.addWidget(self.backgroundcolor, 0, 1)
		gridbox.addWidget(self.hideselectionbutton, 1, 0, 1, 2)
		gridbox.setAlignment(QtCore.Qt.AlignCenter)
		gridbox.setSpacing(10)
		frame.setLayout(gridbox)
		# Buttons frame
		uvbox.addWidget(frame)
		self.opensession_button = QtGui.QPushButton("Open Session")
		self.savesession_button = QtGui.QPushButton("Save Session")
		self.savebutton = QtGui.QPushButton("Save Image Snapshot")
		uvbox.addWidget(self.opensession_button)
		uvbox.addWidget(self.savesession_button)
		uvbox.addWidget(self.savebutton)
		uwidget.setLayout(uvbox)
		
		QtCore.QObject.connect(self.backgroundcolor,QtCore.SIGNAL("newcolor(QColor)"),self._on_bg_color)
		QtCore.QObject.connect(self.hideselectionbutton, QtCore.SIGNAL("clicked()"),self._on_hide)
		QtCore.QObject.connect(self.savebutton, QtCore.SIGNAL("clicked()"),self._on_save)
		QtCore.QObject.connect(self.savesession_button, QtCore.SIGNAL("clicked()"),self._on_save_session)
		QtCore.QObject.connect(self.opensession_button, QtCore.SIGNAL("clicked()"),self._on_open_session)
		
		return uwidget
	
	def _on_hide(self):
		"""
		Hide display selections
		"""
		for node in self.scenegraph().getAllNodes():
				node.setHiddenSelected(self.hideselectionbutton.isChecked())
		self.updateSceneGraph()
		
	def _on_open_session(self):
		"""
		Open a session
		"""
		# Open the file
		filename = QtGui.QFileDialog.getOpenFileName(self, 'Open Session', os.getcwd(), "*.eman")
		if filename:
			self.scenegraph().loadSession(filename)
		
	def _on_save_session(self):
		"""
		Return a list of all the child items (actually a tree of sorts)
		"""
		filename = QtGui.QFileDialog.getSaveFileName(self, 'Save Session', os.getcwd(), "*.eman")
		if filename: # if we cancel
			self.scenegraph().saveSession(filename)

	def _on_save(self):
		"""
		Save a snapshot of the scene
		"""
		filename = QtGui.QFileDialog.getSaveFileName(self, 'Save Image', os.getcwd(), "(*.tiff *.jpeg *.png)")
		if filename: # if we cancel
			self.scenegraph().saveSnapShot(filename)
	
	def _on_bg_color(self, color):
		rgb = color.getRgb()
		self.scenegraph().makeCurrent()
		self.scenegraph().setClearColor(float(rgb[0])/255.0, float(rgb[1])/255.0, float(rgb[2])/255.0)
		self.updateSceneGraph()

		
	def updateInspector(self):
		"""
		Update Inspector,is called whenever the scence changes
		"""
		#tool buttons
		if self.scenegraph().getMouseMode() == "selection": self.selectiontool.setDown(True)
		if self.scenegraph().getMouseMode() == "multiselection": self.multiselectiontool.setDown(True)
		if self.scenegraph().getMouseMode() == "rotate": self.rotatetool.setDown(True)
		if self.scenegraph().getMouseMode() == "xytranslate": self.translatetool.setDown(True)
		if self.scenegraph().getMouseMode() == "ztranslate": self.ztranslate.setDown(True)
		if self.scenegraph().getMouseMode() == "scale": self.scaletool.setDown(True)
		if self.scenegraph().getMouseMode() == "ruler": self.rulertool.setDown(True)
		if self.scenegraph().getMouseMode() == "cube": self.cubetool.setDown(True)
		if self.scenegraph().getMouseMode() == "sphere": self.spheretool.setDown(True)
		if self.scenegraph().getMouseMode() == "cylinder": self.cylindertool.setDown(True)
		if self.scenegraph().getMouseMode() == "cone": self.conetool.setDown(True)
		if self.scenegraph().getMouseMode() == "line": self.linetool.setDown(True)
		if self.scenegraph().getMouseMode() == "text": self.texttool.setDown(True)
		if self.scenegraph().getMouseMode() == "data": self.datatool.setDown(True)
		if self.scenegraph().getMouseMode() == "app": self.apptool.setDown(True)
		# Enable/Disable some tool buttons
		if self.scenegraph().getAPix():
			self.rulertool.setEnabled(True)
		else:
			self.rulertool.setEnabled(False)
			if self.scenegraph().getMouseMode() == "ruler":
				self.scenegraph().setMouseMode("rotate")	# Return to a default state
				self.rotatetool.setDown(True)
		# utils
		self.backgroundcolor.setColor(QtGui.QColor(255*self.scenegraph().clearcolor[0],255*self.scenegraph().clearcolor[1],255*self.scenegraph().clearcolor[2]))
		self.hideselectionbutton.setChecked(self.scenegraph().isSelectionHidded())
	
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
		self.scenegraph().updateSG()

# XPM format Cursors
visibleicon = [
    '16 12 3 1',
    'a c #0000ff',
    'b c #000000',
    'c c None',
    'cccccccccccccccc',
    'ccccccbbbbcccccc',
    'ccccbbbbbbbbcccc',
    'ccbbbccccccbbbcc',
    'cbbccccaaccccbbc',
    'cbccccaaaaccccbc',
    'cbccccaaaaccccbc',
    'cbbccccaaccccbbc',
    'ccbbbccccccbbbcc',
    'ccccbbbbbbbbcccc',
    'ccccccbbbbcccccc',
    'cccccccccccccccc'
]
    
invisibleicon = [
    '16 12 2 1',
    'b c #000000',
    'c c None',
    'cbbcccccccccbbcc',
    'ccbbcccccccbbccc',
    'cccbbcccccbbcccc',
    'ccccbbcccbbccccc',
    'cccccbbcbbcccccc',
    'ccccccbbbccccccc',
    'ccccccbbbccccccc',
    'cccccbbcbbcccccc',
    'ccccbbcccbbccccc',
    'cccbbcccccbbcccc',
    'ccbbcccccccbbccc',
    'cbbcccccccccbbcc'
]

zrotatecursor = [
    '15 14 2 1',
    'b c #00ff00',
    'c c None',
    'ccccccccccccccc',
    'ccccbbbbbbccbcc',
    'ccbbbbbbbbbbbbc',
    'cbbbcccccbbbbbc',
    'bbbcccccbbbbbbb',
    'bbcccccbbbbbbbc',
    'ccccccccccccccc',
    'ccccccccccccccc',
    'cbbbbbbbcccccbb',
    'bbbbbbbcccccbbb',
    'cbbbbbcccccbbbc',
    'cbbbbbbbbbbbbcc',
    'ccbccbbbbbbcccc',
    'ccccccccccccccc'
]

xyrotatecursor = [
    '14 13 2 1',
    'b c #00ff00',
    'c c None',
    'cccccccccccccc',
    'ccccbbbbbbcccc',
    'ccbbbbbbbbbbcc',
    'cbbbccccccbbbc',
    'bbbccccccccbbb',
    'bbbccccccccbbb',
    'bbbccccccccbbb',
    'bbbccccccccbbb',
    'bbbccccccccbbb',
    'cbbbccccccbbbc',
    'ccbbbbbbbbbbcc',
    'ccccbbbbbbcccc',
    'cccccccccccccc'
]

crosshairscursor = [
    '16 16 2 1',
    'b c #00ff00',
    'c c None',
    'cccccccbbcccccccc',
    'ccccccbbbbccccccc',
    'cccccbbbbbbcccccc',
    'cccccccbbcccccccc',
    'cccccccbbcccccccc',
    'ccbccccbbccccbccc',
    'cbbccccbbccccbbcc',
    'bbbbbbbbbbbbbbbbb',
    'bbbbbbbbbbbbbbbbb',
    'cbbccccbbccccbbcc',
    'ccbccccbbccccbccc',
    'cccccccbbcccccccc',
    'cccccccbbcccccccc',
    'cccccbbbbbbcccccc',
    'ccccccbbbbccccccc',
    'cccccccbbcccccccc'
]   

zhaircursor = [
    '16 16 2 1',
    'b c #00ff00',
    'c c None',
    'ccccccccccccccccc',
    'cccbccccccccccccc',
    'ccbbbcccccccccccc',
    'cbcbcbccccccccccc',
    'cccbccccccccccccc',
    'cccbccccccccccccc',
    'cccbccccccbbccccc',
    'cccbccccbbccccccc',
    'cccbccbbccccccccc',
    'cccbbbccccccccccc',
    'cccbccccccccccccc',
    'ccccbbccccccccccc',
    'ccccccbbccccccccc',
    'ccccccccbbccccccc',
    'ccccccccccbbccccc',
    'ccccccccccccccccc'
]

scalecursor = [
    '16 16 2 1',
    'b c #00ff00',
    'c c None',
    'bbbbbbbbcccccccc',
    'bccccccccccccccc',
    'bccccccccccccccc',
    'bccccccccccccccc',
    'bccccccccccccccc',
    'bccccbbbbbbccccc',
    'bccccbccccbccccc',
    'bccccbccccbccccc',
    'cccccbccccbccccb',
    'cccccbccccbccccb',
    'cccccbbbbbbccccb',
    'cccccccccccbcccb',
    'ccccccccccccbccb',
    'cccccccccccccbcb',
    'ccccccccccccccbb',
    'ccccccccbbbbbbbb'
]   

selectorcursor = [
    '16 16 2 1',
    'b c #00ff00',
    'c c None',
    'cbbbbbbbbbcccccc',
    'bcccccccccbccccc',
    'cbbbbbbbcccbcccc',
    'cccbccccccccbccc',
    'ccccbbbbccccbccc',
    'cccbccccccccbccc',
    'ccccbbbbcccbcbcc',
    'cccbccccccbcccbc',
    'ccccbbbbbbcccccb',
    'cccccccbccbcccbc',
    'ccccccccbccccbcc',
    'cccccccccbccbccc',
    'ccccccccccbbcccc',
    'cccccccccccccccc',
    'cccccccccccccccc',
    'cccccccccccccccc'
]

cubecursor = [
    '16 16 2 1',
    'b c #00ff00',
    'c c None',
    'cccccccccccccccccc',
    'cccccccccccccccccc',
    'ccccccbbbbbbbbbccc',
    'cccccbbbbbbbbbbccc',
    'ccccbbbbbbbbbbbccc',
    'cccbbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbcccc',
    'ccbbbbbbbbbbbccccc',
    'ccbbbbbbbbbbcccccc',
    'ccbbbbbbbbbccccccc',
    'cccccccccccccccccc',
    'cccccccccccccccccc'
]

spherecursor = [
    '16 16 2 1',
    'b c #00ff00',
    'c c None',
    'ccccccccccccccccc',
    'cccccccbbbccccccc',
    'cccccbbbbbbbccccc',
    'ccccbbbbbbbbbcccc',
    'cccbbbbbbbbbbbccc',
    'cccbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'cccbbbbbbbbbbbccc',
    'cccbbbbbbbbbbbccc',
    'ccccbbbbbbbbbcccc',
    'cccccbbbbbbbccccc',
    'cccccccbbbccccccc',
    'ccccccccccccccccc',
    'ccccccccccccccccc'
]

cylindercursor = [
    '16 16 2 1',
    'b c #00ff00',
    'c c None',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'cccccbbbbbbbccccc',
    'cccbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'cccbbbbbbbbbbbccc',
    'cccccbbbbbbbccccc',
    'ccccccccccccccccc'
] 

conecursor = [
    '16 16 2 1',
    'b c #00ff00',
    'c c None',
    'cccccccbbcccccccc',
    'ccccccbbbbccccccc',
    'ccccccbbbbccccccc',
    'cccccbbbbbbcccccc',
    'cccccbbbbbbcccccc',
    'ccccbbbbbbbbccccc',
    'ccccbbbbbbbbccccc',
    'cccbbbbbbbbbbcccc',
    'cccbbbbbbbbbbcccc',
    'ccbbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbccc',
    'cbbbbbbbbbbbbbbcc',
    'cbbbbbbbbbbbbbbcc',
    'cccbbbbbbbbbbcccc',
    'cccccbbbbbbcccccc',
    'ccccccccccccccccc'
]

linecursor = [
    '16 16 2 1',
    'b c #00ff00',
    'c c None',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'cccccccccccccbbcc',
    'ccccccccccccbbccc',
    'cccccccccccbbcccc',
    'ccccccccccbbccccc',
    'cccccccccbbcccccc',
    'ccccccccbbccccccc',
    'cccccccbbcccccccc',
    'ccccccbbccccccccc',
    'cccccbbcccccccccc',
    'ccccbbccccccccccc',
    'cccbbcccccccccccc',
    'ccbbccccccccccccc',
    'ccccccccccccccccc',
    'ccccccccccccccccc'
] 

datacursor = [
    '17 16 2 1',
    'b c #00ff00',
    'c c None',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'bbbbbbbcbbcccccbb',
    'bbbbbbbcbbbcccbbb',
    'bbccccccbbbbcbbbb',
    'bbccccccbbbbbbbbb',
    'bbccccccbbcbbbcbb',
    'bbbbbcccbbcccccbb',
    'bbbbbcccbbcccccbb',
    'bbccccccbbcccccbb',
    'bbccccccbbcccccbb',
    'bbccccccbbcccccbb',
    'bbbbbbbcbbcccccbb',
    'bbbbbbbcbbcccccbb',
    'ccccccccccccccccc',
    'ccccccccccccccccc'
] 

textcursor = [
    '16 16 2 1',
    'b c #00ff00',
    'c c None',
    'ccccccccccccccccc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbcccbbbcccbbcc',
    'ccbccccbbbccccbcc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'ccccccbbbbbcccccc',
    'cccccbbbbbbbccccc',
    'ccccccccccccccccc'
] 

appcursor = [
    '17 16 2 1',
    'b c #00ff00',
    'c c None',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'caaacccaaacccaaac',
    'acccacacccacaccca',
    'acccacacccacaccca',
    'acccacacccacaccca',
    'acccacacccacaccca',
    'aaaaacaaaaccaaaac',
    'acccacacccccacccc',
    'acccacacccccacccc',
    'acccacacccccacccc',
    'acccacacccccacccc',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'ccccccccccccccccc'
]

cubeicon = [
    '16 16 2 1',
    'b c #000055',
    'c c None',
    'cccccccccccccccccc',
    'cccccccccccccccccc',
    'ccccccbbbbbbbbbccc',
    'cccccbbbbbbbbbbccc',
    'ccccbbbbbbbbbbbccc',
    'cccbbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbcccc',
    'ccbbbbbbbbbbbccccc',
    'ccbbbbbbbbbbcccccc',
    'ccbbbbbbbbbccccccc',
    'cccccccccccccccccc',
    'cccccccccccccccccc'
] 

sphereicon = [
    '16 16 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccccc',
    'cccccccbbbccccccc',
    'cccccbbbbbbbccccc',
    'ccccbbbbbbbbbcccc',
    'cccbbbbbbbbbbbccc',
    'cccbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'cccbbbbbbbbbbbccc',
    'cccbbbbbbbbbbbccc',
    'ccccbbbbbbbbbcccc',
    'cccccbbbbbbbccccc',
    'cccccccbbbccccccc',
    'ccccccccccccccccc',
    'ccccccccccccccccc'
] 

cylindericon = [
    '16 16 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'cccccbbbbbbbccccc',
    'cccbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'cccbbbbbbbbbbbccc',
    'cccccbbbbbbbccccc',
    'ccccccccccccccccc'
] 

coneicon = [
    '16 16 2 1',
    'b c #000055',
    'c c None',
    'cccccccbbcccccccc',
    'ccccccbbbbccccccc',
    'ccccccbbbbccccccc',
    'cccccbbbbbbcccccc',
    'cccccbbbbbbcccccc',
    'ccccbbbbbbbbccccc',
    'ccccbbbbbbbbccccc',
    'cccbbbbbbbbbbcccc',
    'cccbbbbbbbbbbcccc',
    'ccbbbbbbbbbbbbccc',
    'ccbbbbbbbbbbbbccc',
    'cbbbbbbbbbbbbbbcc',
    'cbbbbbbbbbbbbbbcc',
    'cccbbbbbbbbbbcccc',
    'cccccbbbbbbcccccc',
    'ccccccccccccccccc'
]

lineicon = [
    '16 16 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'cccccccccccccbbcc',
    'ccccccccccccbbccc',
    'cccccccccccbbcccc',
    'ccccccccccbbccccc',
    'cccccccccbbcccccc',
    'ccccccccbbccccccc',
    'cccccccbbcccccccc',
    'ccccccbbccccccccc',
    'cccccbbcccccccccc',
    'ccccbbccccccccccc',
    'cccbbcccccccccccc',
    'ccbbccccccccccccc',
    'ccccccccccccccccc',
    'ccccccccccccccccc'
] 

texticon = [
    '16 16 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccccc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbbbbbbbbbbbbcc',
    'ccbbcccbbbcccbbcc',
    'ccbccccbbbccccbcc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'cccccccbbbccccccc',
    'ccccccbbbbbcccccc',
    'cccccbbbbbbbccccc',
    'ccccccccccccccccc'
] 

dataicon = [
    '17 16 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'bbbbbbbcbbcccccbb',
    'bbbbbbbcbbbcccbbb',
    'bbccccccbbbbcbbbb',
    'bbccccccbbbbbbbbb',
    'bbccccccbbcbbbcbb',
    'bbbbbcccbbcccccbb',
    'bbbbbcccbbcccccbb',
    'bbccccccbbcccccbb',
    'bbccccccbbcccccbb',
    'bbccccccbbcccccbb',
    'bbbbbbbcbbcccccbb',
    'bbbbbbbcbbcccccbb',
    'ccccccccccccccccc',
    'ccccccccccccccccc'
] 

appicon = [
    '17 16 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'caaacccaaacccaaac',
    'acccacacccacaccca',
    'acccacacccacaccca',
    'acccacacccacaccca',
    'acccacacccacaccca',
    'aaaaacaaaaccaaaac',
    'acccacacccccacccc',
    'acccacacccccacccc',
    'acccacacccccacccc',
    'acccacacccccacccc',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'ccccccccccccccccc'
]

rotateicon = [
    '15 14 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccc',
    'ccccbbbbbbccbcc',
    'ccbbbbbbbbbbbbc',
    'cbbbcccccbbbbbc',
    'bbbcccccbbbbbbb',
    'bbcccccbbbbbbbc',
    'ccccccccccccccc',
    'ccccccccccccccc',
    'cbbbbbbbcccccbb',
    'bbbbbbbcccccbbb',
    'cbbbbbcccccbbbc',
    'cbbbbbbbbbbbbcc',
    'ccbccbbbbbbcccc',
    'ccccccccccccccc'
]

crosshairsicon = [
    '16 16 2 1',
    'b c #000055',
    'c c None',
    'cccccccbbcccccccc',
    'ccccccbbbbccccccc',
    'cccccbbbbbbcccccc',
    'cccccccbbcccccccc',
    'cccccccbbcccccccc',
    'ccbccccbbccccbccc',
    'cbbccccbbccccbbcc',
    'bbbbbbbbbbbbbbbbb',
    'bbbbbbbbbbbbbbbbb',
    'cbbccccbbccccbbcc',
    'ccbccccbbccccbccc',
    'cccccccbbcccccccc',
    'cccccccbbcccccccc',
    'cccccbbbbbbcccccc',
    'ccccccbbbbccccccc',
    'cccccccbbcccccccc'
] 

multiselectoricon = [
    '17 16 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'ccccccccccccccccc',
    'cbcbcbcbcbcbcbccc',
    'ccccccccccccccccc',
    'cbcccccccccccbccc',
    'ccccccccccccccccc',
    'cbcccccccccccbccc',
    'ccccccccccccccccc',
    'cbcccccccccccbccc',
    'cccccccccccccbbbc',
    'cbcbcbcbcbcbcbbbc',
    'cccccccccccccbbbb',
    'ccccccccccccccccb',
    'ccccccccccccccccc'
] 

scaleicon = [
    '16 16 2 1',
    'b c #000055',
    'c c None',
    'bbbbbbbbcccccccc',
    'bccccccccccccccc',
    'bccccccccccccccc',
    'bccccccccccccccc',
    'bccccccccccccccc',
    'bccccbbbbbbccccc',
    'bccccbccccbccccc',
    'bccccbccccbccccc',
    'cccccbccccbccccb',
    'cccccbccccbccccb',
    'cccccbbbbbbccccb',
    'cccccccccccbcccb',
    'ccccccccccccbccb',
    'cccccccccccccbcb',
    'ccccccccccccccbb',
    'ccccccccbbbbbbbb'
]

ztransicon = [
    '16 16 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccccc',
    'cccbccccccccccccc',
    'ccbbbcccccccccccc',
    'cbcbcbccccccccccc',
    'cccbccccccccccccc',
    'cccbccccccccccccc',
    'cccbccccccbbccccc',
    'cccbccccbbccccccc',
    'cccbccbbccccccccc',
    'cccbbbccccccccccc',
    'cccbccccccccccccc',
    'ccccbbccccccccccc',
    'ccccccbbccccccccc',
    'ccccccccbbccccccc',
    'ccccccccccbbccccc',
    'ccccccccccccccccc'
] 

selectionicon = [
    '16 16 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccccc',
    'ccbcccccccccccccc',
    'ccbbccccccccccccc',
    'ccbcbbccccccccccc',
    'ccbcccbcccccccccc',
    'ccbccccbbcccccccc',
    'ccbccccccbccccccc',
    'ccbcccccccbbccccc',
    'ccbcccccbbbbbcccc',
    'ccbccbccbcccccccc',
    'ccbcbcbccbccccccc',
    'ccbbccbccbccccccc',
    'ccbccccbccbcccccc',
    'ccccccccbbbcccccc',
    'ccccccccccccccccc',
    'ccccccccccccccccc'
] 

rulericon = [
    '16 16 2 1',
    'b c #000055',
    'c c None',
    'cccbbbbbbbbbbbccc',
    'cccbcccccccccbccc',
    'cccbbbcccccccbccc',
    'cccbcccccccccbccc',
    'cccbbbbbcccccbccc',
    'cccbcccccccccbccc',
    'cccbbbcccccccbccc',
    'cccbcccccccccbccc',
    'cccbbbbbcccccbccc',
    'cccbcccccccccbccc',
    'cccbbbcccccccbccc',
    'cccbcccccccccbccc',
    'cccbbbbbcccccbccc',
    'cccbcccccccccbccc',
    'cccbbbbbbbbbbbccc',
    'ccccccccccccccccc'
] 

if __name__ == '__main__':
	main()
