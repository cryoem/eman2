#!/usr/bin/env python
#
# Author: Ross Coleman (racolema@bcm.edu)
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

from OpenGL import GL
from PyQt4 import QtCore, QtGui
from EMAN2 import Transform, Vec4f, Vec3f
from libpyGLUtils2 import GLUtil
from valslider import ValSlider, EMSpinWidget
import weakref
import math
import numpy

class EMVertexBuffers(object):
	""" Class to hold vertex buffers """
	vbuffers = []
	numbuffers = 2
	
	@staticmethod
	def getBuffers():
		if not len(EMVertexBuffers.vbuffers): 
			EMVertexBuffers.vbuffers =  GL.glGenBuffers(EMVertexBuffers.numbuffers)
		return EMVertexBuffers.vbuffers
	
def drawBoundingBox(width, height, depth):
	""" Draw bounding box use C side for easy of vertex arrays """
	
	GL.glPushAttrib(GL.GL_LIGHTING_BIT)
	GL.glMaterialfv(GL.GL_FRONT, GL.GL_EMISSION, [1.0,1.0,1.0,1.0])
	GLUtil.glDrawBoundingBox(width, height, depth)
	GL.glPopAttrib()
	
class EMItem3D(object): #inherit object for new-style class (new-stype classes required for super() and Python properties)
	"""
	The base class for nodes in our scene graph, which is used in our 3D scene viewer.
	In our case, the scene graph is a tree data structure.
	"""
	# Class attrib to connect openGL int identifiers to class instances
	vbuffers = 0
	selection_idx_dict = {}
	selection_recycle = []
	selection_intname = -1
	name = "General 3D Item"
	nodetype = "BaseNode"
	
	@staticmethod
	def get_transformlayout(layout, idx, attribdict):
		"""
		Adds a transform layout to a dialog
		@param layout, the layout to append to
		@param idx, the row to being appnding to
		"""
		font = QtGui.QFont()
		font.setBold(True)
		translatelabel = QtGui.QLabel("Translation")
		translatelabel.setFont(font)
		translatelabel.setAlignment(QtCore.Qt.AlignCenter)
		layout.addWidget(translatelabel, idx, 0, 1, 4)
		txlabel = QtGui.QLabel("Tx")
		tylabel = QtGui.QLabel("Ty")
		txlabel.setAlignment(QtCore.Qt.AlignRight)
		tylabel.setAlignment(QtCore.Qt.AlignRight)
		attribdict["tx"] = QtGui.QLineEdit("0.0")
		attribdict["ty"] = QtGui.QLineEdit("0.0")
		attribdict["tx"].setMinimumWidth(100.0)
		attribdict["ty"].setMinimumWidth(100.0)
		layout.addWidget(txlabel, idx+1, 0, 1, 1)
		layout.addWidget(attribdict["tx"], idx+1, 1, 1, 1)
		layout.addWidget(tylabel, idx+1, 2, 1, 1)
		layout.addWidget(attribdict["ty"], idx+1, 3, 1, 1)
		tzlabel = QtGui.QLabel("Tz")
		zoomlabel = QtGui.QLabel("Zm")
		tylabel.setAlignment(QtCore.Qt.AlignRight)
		zoomlabel.setAlignment(QtCore.Qt.AlignRight)
		attribdict["tz"] = QtGui.QLineEdit("0.0")
		attribdict["zoom"] = QtGui.QLineEdit("1.0")
		attribdict["tz"].setMinimumWidth(100.0)
		attribdict["zoom"].setMinimumWidth(100.0)
		layout.addWidget(tzlabel, idx+2, 0, 1, 1)
		layout.addWidget(attribdict["tz"], idx+2, 1, 1, 1)
		layout.addWidget(zoomlabel, idx+2, 2, 1, 1)
		layout.addWidget(attribdict["zoom"], idx+2, 3, 1, 1)
		rotatelabel = QtGui.QLabel("EMAN Rotation")
		rotatelabel.setFont(font)
		rotatelabel.setAlignment(QtCore.Qt.AlignCenter)
		layout.addWidget(rotatelabel, idx+3, 0, 1, 4)
		azlabel = QtGui.QLabel("Az")
		azlabel.setAlignment(QtCore.Qt.AlignRight)
		attribdict["az"] = QtGui.QLineEdit("0.0")
		altlabel = QtGui.QLabel("Alt")
		altlabel.setAlignment(QtCore.Qt.AlignRight)
		attribdict["alt"] = QtGui.QLineEdit("0.0")
		attribdict["az"] .setMinimumWidth(100.0)
		attribdict["alt"].setMinimumWidth(100.0)
		layout.addWidget(azlabel, idx+4, 0, 1, 1)
		layout.addWidget(attribdict["az"], idx+4, 1, 1, 1)
		layout.addWidget(altlabel, idx+4, 2, 1, 1)
		layout.addWidget(attribdict["alt"], idx+4, 3, 1, 1)
		philabel = QtGui.QLabel("Phi")
		philabel.setAlignment(QtCore.Qt.AlignRight)
		attribdict["phi"] = QtGui.QLineEdit("0.0")
		layout.addWidget(philabel, idx+5, 0, 1, 1)
		layout.addWidget(attribdict["phi"], idx+5, 1, 1, 1)
	
	@staticmethod
	def getTransformFromDict(attribdict):
		""" Return a transform using a dict created using the above function"""
		return Transform({"type":"eman","az":float(attribdict["az"].text()),"alt":float(attribdict["alt"].text()),"phi":float(attribdict["phi"].text()),"tx":float(attribdict["tx"].text()),"ty":float(attribdict["ty"].text()),"tz":float(attribdict["tz"].text()),"scale":float(attribdict["zoom"].text())})
	
	@staticmethod
	def getNodeDialogWidget(attribdict):
		"""
		Get Data Widget
		"""
		nodewidget = QtGui.QWidget()
		grid = QtGui.QGridLayout()
		node_name_label = QtGui.QLabel("Node Name")
		attribdict["node_name"] = QtGui.QLineEdit()
		grid.addWidget(node_name_label , 0, 0, 1, 2)
		grid.addWidget(attribdict["node_name"], 0, 2, 1, 2)
		EMItem3D.get_transformlayout(grid, 2, attribdict)
		nodewidget.setLayout(grid)
		
		return nodewidget
	
	@staticmethod
	def getNodeForDialog(attribdict):
		"""
		Create a new node using a attribdict
		"""
		return EMItem3D(attribdict["parent"], transform=EMItem3D.getTransformFromDict(attribdict))
		
	def __init__(self, parent = None, children = [], transform=None):
		"""
		@type parent: EMItem3D
		@param parent: the parent node to the current node or None for the root node
		@type children: list
		@param children: the child nodes
		@type transform: Transform or None
		@param transform: The transformation (rotation, scaling, translation) that should be applied before rendering this node and its children 
		"""
		self.label = None	# Customisabl label, used to label the inspector in the tree
		self.setParent(parent)
		self.setChildren(children)
		if not transform: transform = Transform()
		self.transform = transform
		self.is_visible = True 
		self.is_selected = False
		self.hide_selection = False
		self.item_inspector = None			# This is an inspector widget
		self.EMQTreeWidgetItem = None 			# This is an inspector tree item
		self.boundingboxsize = None
		self.getAndSetUniqueInteger()
	
	def getChildren(self): return self.children
	def setChildren(self, children): 
		self.children = list(children)
		for child in children:
			child.parent = self
	def getParent(self): return self.parent
	def setParent(self, parent): 
		self.parent = parent
		#if parent:	# Why is this here is causes chaos b/c things get added twice. You can only implment this if you check to ensure self is not already a child of parent!
			#parent.addChild(self)
	def isSelectedItem(self): return self.is_selected
	def setSelectedItem(self, is_selected): self.is_selected = is_selected
	def isSelectionHidded(self): return self.hide_selection
	def getHiddenSelected(self): return self.hide_selection
	def setHiddenSelected(self, hide_selection): self.hide_selection = hide_selection
	def getTransform(self): return self.transform
	def setTransform(self, transform): self.transform = transform
	def isVisibleItem(self): return self.is_visible
	def setVisibleItem(self, is_visible): self.is_visible = is_visible
	def setLabel(self, label): self.label = label
	def getLabel(self): return self.label
	
	def __del__(self):
		EMItem3D.selection_recycle.append(self.intname)
	
	def getEvalString(self):
		"""
		Retrun a string that after eval can reinstatiate the object
		"""
		return "EMItem3D()"
	
	def getAndSetUniqueInteger(self):
		"""
		Stuff for the selection mechanism, return a unique int for each instance of EMItem3D
		"""
		if len(EMItem3D.selection_recycle) > 0:
			self.intname = EMItem3D.selection_recycle.pop()
		else:
			EMItem3D.selection_intname += 1
			self.intname = EMItem3D.selection_intname
		EMItem3D.selection_idx_dict[self.intname] = weakref.ref(self)
		
	def addChild(self, node):
		"""
		Adds a child node, if not already in the list of child nodes.
		@type node: EMItem3D
		@param node: the child node to add
		"""
		self.children.append(node)
		node.parent = self
	
	def insertChild(self, node, nodeindex):
		"""
		Inserts a child node, if not already in the list of child nodes.
		@type node: EMItem3D
		@param node: the child node to add
		@param nodeindex: the index for insertion
		"""
		self.children.insert(nodeindex, node)
		node.parent = self
		
	def addChildren(self, nodes):
		"""
		Adds all the provided child nodes which are not already in the list of child nodes.
		@type nodes: an iterable collection of EMItem3D (or subclass) objects
		@param nodes: the nodes which will be added as child nodes
		"""
		for node in nodes:
			if not node in self.children:
				self.children.append(node)
				node.parent = self
			
	def hasChild(self, node):
		"""
		Tests whether the supplied node is a child of the current node. 
		@type node: EMItem3D
		@param node: test whether this node is a child node of self 
		"""
		return node in self.children
		
	def getChildIndex(self, node):
		"""
		Returns the index of child in the list of children. Raises an execption if not found
		"""
		return self.children.index(node)
	
	def displayTree(self, level = 1):
		"""
		This prints a representation of the subtree rooted at this node to standard output. Useful for debugging.
		"""
		indent = "\t"*(level-1)
		marker = "<-->" if self.parent else "-->"
		print indent, marker, self.intname
		for child in self.children:
			child.displayTree(level+1)
	
	def addParentReferences(self):
		"""
		For the subtree rooted at self, give each child node a reference to its parent node.
		"""
		for child in self.children:
			child.parent = self
			child.addParentReferences()
	
	def removeParentReferences(self):
		"""
		For the subtree rooted at self, break reference cycles by setting self.parent to None.
		"""
		self.parent = None
		for child in self.children:
			child.removeParentReferences()
	
	def removeChild(self, node):
		"""
		Remove the supplied node from the list of child nodes This also removes all its descendant nodes. 
		@type node: EMItem3D
		@param node: the node to remove
		"""
		if node in self.children:
			node.removeParentReferences()
		while node in self.children: #In case a duplicate got in, somehow
			self.children.remove(node)
	
	def getSelectedAncestorNodes(self):
		"""
		Return a list of selected node ancestors
		"""
		selected_ancestors = []
		node = self
		while node.parent:
			if node.parent.is_selected:
				selected_ancestors.append(node.parent)
			node = node.parent
			
		return selected_ancestors
	
	def getRootNode(self):
		"""
		Returns the root node
		"""
		node = self
		while node.parent:
			node = node.parent
		return node
		
	def getAllSelectedNodes(self):
		"""
		For the tree rooted at self, this recursive method returns a list of all the selected nodes.
		@return: a list of selected nodes
		"""
		selected_list = []
		if self.is_selected:
			selected_list.append(self)
		for child in self.children: #Recursion ends on leaf nodes here
			selected_list.extend(child.getAllSelectedNodes()) #Recursion
		
		return selected_list
	
	def getAllNodes(self):
		"""
		For the tree rooted at self, this recursive method returns a list of all the nodes.
		@return: a list of selected nodes
		"""
		nodelist = []
		nodelist.append(self)
		for child in self.children: #Recursion ends on leaf nodes here
			nodelist.extend(child.getAllNodes()) #Recursion
		
		return nodelist
		
		
	def getNearbySelectedNodes(self):
		"""
		For the tree rooted at self, this recursive method returns a list of the selected nodes that are near self.
		A selected node will not be in the returned list if one of its ancestor nodes is also selected. 
		@return: a list of selected nodes
		"""
		selected_list = []
		if self.is_selected:
			return [self]
		else:
			for child in self.children:
				selected_list.extend(child.getNearbySelectedNodes())
		
		return selected_list
	
	def getDistantSelectedNodes(self):
		"""
		For the tree rooted at self, this recursive method returns a list of the selected nodes that are distant from self.
		A selected node will not be in the returned list if one of its descendant nodes is also selected. 
		@return: a list of selected nodes
		"""
		selected_list = []
		for child in self.children:
			selected_list.extend(child.getDistantSelectedNodes())
		if not selected_list: #either this is a leaf node, or there are no selected nodes in the child subtrees
			if self.is_selected:
				selected_list.append(self)
		
		return selected_list

	def getItemInspector(self):
		"""
		Return a Qt widget that controls the scene item
		"""
		if not self.item_inspector: self.item_inspector = EMItem3DInspector("Node", self)
		return self.item_inspector
		
	def setEMQTreeWidgetItem(self, node):
		"""
		Relate a QtreeItem to this node
		"""
		self.EMQTreeWidgetItem = node
	
	def getTransformStdCoord(self):
		"""
		This returns the transform in standard coordinate system, not one changed by7 parent matrices
		"""
		tt = t = self.transform
		tp = self.getParentMatrixProduct()
		if tp: tt = tp*t
		return tt
		
	def getParentMatrixProduct(self):
		"""
		Get the product of all parent matrices
		This is a recursive function
		"""
		if self.parent:
			if self.parent.getParentMatrixProduct():
				return self.parent.getParentMatrixProduct()*self.parent.getTransform()
			else:
				return self.parent.getTransform()
		else:
			return None
		
	def updateMatrices(self, params, xformtype):
		"""
		The matrcies are updated in such a way that each is done in the standard cooridnate system and the std coord system is not perturbed by the others
		@type params: List
		@param params: A list defining how the transform in each active node is modified
		@type xfromtype: sting
		@param xformtype: The sort of transform we wish to do
		"""	
		if self.is_selected:
			if xformtype == "rotate":
				if self.parent:
					self.transform.rotate_origin_newbasis(self.getParentMatrixProduct(), params[0], params[1], params[2], params[3])
				else:
					self.transform.rotate_origin(Transform({"type":"spin","omega":params[0],"n1":params[1],"n2":params[2],"n3":params[3]}))
			elif xformtype == "translate":
				if self.parent: 
					self.transform.translate_newbasis(self.getParentMatrixProduct(), params[0], params[1], params[2])
				else:
					self.transform.translate(params[0], params[1], params[2])
					
			elif xformtype == "scale":
				self.transform.scale(params[0])
			else:
				raise Exception,"Invalid transformation type"
		
		# Now tell all children to update
		# TODO: we MIGHT want to get rid of the recursive part of this algorithm
		#        instead, the calling function would get a list of all selected nodes, and apply transformations to each
		for child in self.children:
			child.updateMatrices(params, xformtype) #Note: the transformation is only applied to SELECTED nodes
	
	def getItemDictionary(self):
		"""
		Return a dictionary of item parameters (used for restoring sessions
		"""
		return {"TRANSFORMATION":self.transform.get_params("eman"),"CONSTRUCTOR":self.getEvalString(),"NAME":str(self.getLabel()),"VISIBLE":self.isVisibleItem(),"SELECTED":self.isSelectedItem(),"NODETYPE":self.nodetype,"HIDDENSEL":self.getHiddenSelected()}
		
	def setUsingDictionary(self, dictionary):
		"""
		Set item attributes using a dictionary, used in session restoration
		"""
		self.setTransform(Transform(dictionary["TRANSFORMATION"]))
		self.setVisibleItem(dictionary["VISIBLE"])
		self.setSelectedItem(dictionary["SELECTED"])
		try:
			self.setHiddenSelected(dictionary["HIDDENSEL"])
		except:
			pass
		self.setLabel(dictionary["NAME"])
		
	def render(self):
		"""
		This is the method to call to render the node and its child nodes. 
		It calls self.renderNode() to render the current node. 
		Usually, this method is unchanged in subclasses. 
		"""
		if not self.is_visible:
			return #Also applies to subtree rooted at this node

		if self.transform.is_identity():
			GL.glPushName(self.intname)
			self.renderNode()
			for child in self.children:
				child.render()
			GL.glPopName()
		else:
			GL.glPushMatrix()
			GL.glPushName(self.intname)
			GLUtil.glMultMatrix(self.transform) #apply the transformation
			
			self.renderNode()
			for child in self.children:
				child.render()
			GL.glPopName()
			GL.glPopMatrix()

		if self.item_inspector != None and self.is_selected: 
				self.item_inspector.updateItemControls()

	def renderNode(self):
		"""
		This method, which is called by self.render(), renders the current node.
		It should be implemented in subclasses that represent visible objects.
		"""
		pass

	def keyPressEvent(self, event): pass
	def keyReleaseEvent(self, event): pass
	def mouseDoubleClickEvent(self, event): pass
	def mouseMoveEvent(self, event): pass
	def mousePressEvent(self, event): pass
	def mouseReleaseEvent(self, event): pass
	def wheelEvent(self, event): pass
		
class EMItem3DInspector(QtGui.QTabWidget):
	"""
	Class to make the EMItem GUI controls
	"""
	def __init__(self, name, item3d):
		QtGui.QTabWidget.__init__(self)
		self.item3d = weakref.ref(item3d)
		self.name = name
		self.inspector = None
		self.transfromboxmaxheight = 400	# This might be problematic
		
		self.addTabs()	
			
	def setInspector(self, inspector):
		""" This is a reference back to the main inspector, which holds all the item inspectors"""
		self.inspector = weakref.ref(inspector)
		
	def addTabs(self):
		""" Add a tab for each 'column' """
		tabwidget = QtGui.QWidget()
		gridbox = QtGui.QGridLayout()
		
		EMItem3DInspector.addControls(self, gridbox)
		
		tabwidget.setLayout(gridbox)
		self.addTab(tabwidget, "basic")
		
	def addControls(self, gridbox):
		""" Construct all the widgets in this Item Inspector """
		# selection box and label
		font = QtGui.QFont()
		font.setBold(True)
		label = QtGui.QLabel(self.name,self)
		label.setFont(font)
		label.setAlignment(QtCore.Qt.AlignCenter)
		gridbox.addWidget(label, 0, 0, 1, 1)
		databox = QtGui.QHBoxLayout()
		self.boundingbox = None
		if self.item3d().boundingboxsize:
			self.boundingbox = QtGui.QLabel("Size: "+self.item3d().boundingboxsize,self)
			databox.addWidget(self.boundingbox)
		gridbox.addLayout(databox, 1, 0, 1, 1)
		# angluar controls
		xformframe = QtGui.QFrame()
		xformframe.setFrameShape(QtGui.QFrame.StyledPanel)
		xformbox = QtGui.QGridLayout()
		xformlabel = QtGui.QLabel("Transformation", xformframe)
		xformlabel.setFont(font)
		xformlabel.setAlignment(QtCore.Qt.AlignCenter)
		xformbox.addWidget(xformlabel, 0, 0, 1, 2)
		# Rotations
		self.rotcombobox = QtGui.QComboBox()
		xformbox.addWidget(self.rotcombobox, 1, 0, 1, 2)
		self.rotstackedwidget = QtGui.QStackedWidget()
		self.addRotationWidgets()
		xformbox.addWidget(self.rotstackedwidget, 2, 0, 1, 2)
		#translations
		txlabel = QtGui.QLabel("TX",xformframe)
		txlabel.setAlignment(QtCore.Qt.AlignCenter)
		xformbox.addWidget(txlabel, 3, 0, 1, 1)
		tylabel = QtGui.QLabel("TY",xformframe)
		tylabel.setAlignment(QtCore.Qt.AlignCenter)
		xformbox.addWidget(tylabel, 3, 1, 1, 1)
		self.tx = EMSpinWidget(0.0, 1.0)
		self.ty = EMSpinWidget(0.0, 1.0)
		xformbox.addWidget(self.tx, 4, 0, 1, 1)
		xformbox.addWidget(self.ty, 4, 1, 1, 1)
		tzlabel = QtGui.QLabel("TZ",xformframe)
		tzlabel.setAlignment(QtCore.Qt.AlignCenter)
		xformbox.addWidget(tzlabel, 5, 0, 1, 1)
		zoomlabel = QtGui.QLabel("Zoom",xformframe)
		zoomlabel.setAlignment(QtCore.Qt.AlignCenter)
		xformbox.addWidget(zoomlabel, 5, 1, 1, 1)
		self.tz = EMSpinWidget(0.0, 1.0)
		self.zoom = EMSpinWidget(1.0, 0.1, postivemode=True, wheelstep=0.1)
		xformbox.addWidget(self.tz, 6, 0, 1, 1)
		xformbox.addWidget(self.zoom, 6, 1, 1, 1)
		self.resetbuttontx = QtGui.QPushButton("Reset Tx")
		self.resetbuttonrot = QtGui.QPushButton("Reset Rot")
		xformbox.addWidget(self.resetbuttontx, 7, 0, 1, 1)
		xformbox.addWidget(self.resetbuttonrot, 7, 1, 1, 1)
		xformframe.setLayout(xformbox)
		xformframe.setMaximumWidth(350)
				
		xformframe.setMaximumHeight(self.transfromboxmaxheight)
		xformframe.setLayout(xformbox)
		gridbox.addWidget(xformframe, 2, 0, 1, 1)
		
		# set to default, but run only as a base class
		if type(self) == EMItem3DInspector: self.updateItemControls()
		
		QtCore.QObject.connect(self.tx,QtCore.SIGNAL("valueChanged(int)"),self._on_translation)
		QtCore.QObject.connect(self.ty,QtCore.SIGNAL("valueChanged(int)"),self._on_translation)
		QtCore.QObject.connect(self.tz,QtCore.SIGNAL("valueChanged(int)"),self._on_translation)
		QtCore.QObject.connect(self.zoom,QtCore.SIGNAL("valueChanged(int)"),self._on_scale)
		QtCore.QObject.connect(self.resetbuttontx,QtCore.SIGNAL("clicked()"),self._on_resettx)
		QtCore.QObject.connect(self.resetbuttonrot,QtCore.SIGNAL("clicked()"),self._on_resetrot)
	
	def _on_translation(self, value):
		"""
		Need to contain the right coords. And do translation in the correct corrd system
		"""
		tt = t = Transform({"tx":self.tx.getValue(),"ty":self.ty.getValue(),"tz":self.tz.getValue()})
		tp = self.item3d().getParentMatrixProduct()
		if tp: tt = tp.inverse()*t
		self.item3d().getTransform().set_trans(tt.get_trans())
		self.inspector().updateSceneGraph()
		
	def _on_scale(self, value):
		self.item3d().getTransform().set_scale(self.zoom.getValue())
		self.inspector().updateSceneGraph()

	def _on_resettx(self):
		
		self.item3d().getTransform().set_trans(0.0, 0.0, 0.0)
		self.updateItemControls()
		self.inspector().updateSceneGraph()
		
	def _on_resetrot(self):
		self.item3d().getTransform().set_rotation({"type":"eman","az":0.0,"alt":0.0,"phi":0.0})
		self.updateItemControls()
		self.inspector().updateSceneGraph()
	
	def _isRotNaN(self, rot1, rot2, rot3):
		""" Better check to make sure get_rotation did not return Nan, so to prevent a crash """
		if rot1 != rot1: return True
		if rot2 != rot2: return True
		if rot3 != rot3: return True
		return False
		
	def updateItemControls(self):
		""" Updates this item inspector. Function is called by the item it observes"""
		# Translation update
		stdtransfrom = self.item3d().getTransformStdCoord()
		translation =  stdtransfrom.get_trans()
		
		self.tx.setValue(translation[0])
		self.ty.setValue(translation[1])
		self.tz.setValue(translation[2])
		# Rotation update
		rotation =  stdtransfrom.get_rotation(str(self.rotcombobox.currentText()))
		is_identity = stdtransfrom.is_rot_identity()
		comboboxidx = self.rotcombobox.currentIndex()
		if comboboxidx == 0:
			if self._isRotNaN(rotation["az"],rotation["alt"],rotation["phi"]): return
			self.emanazslider.setValue(rotation["az"], quiet=1)
			self.emanaltslider.setValue(rotation["alt"], quiet=1)
			self.emanphislider.setValue(rotation["phi"], quiet=1)
		if comboboxidx == 1:
			if self._isRotNaN(rotation["gamma"],rotation["beta"],rotation["alpha"]): return
			self.imagicgammaslider.setValue(rotation["gamma"], quiet=1)
			self.imagicbetaslider.setValue(rotation["beta"], quiet=1)
			self.imagicalphaslider.setValue(rotation["alpha"], quiet=1)
		if comboboxidx == 2:
			if self._isRotNaN(rotation["psi"],rotation["theta"],rotation["phi"]): return
			self.spiderpsislider.setValue(rotation["psi"], quiet=1)
			self.spiderthetaslider.setValue(rotation["theta"], quiet=1)
			self.spiderphislider.setValue(rotation["phi"], quiet=1)
		if comboboxidx == 3:
			if self._isRotNaN(rotation["phi"],rotation["theta"],rotation["omega"]): return
			self.mrcpsislider.setValue(rotation["phi"], quiet=1)
			self.mrcthetaslider.setValue(rotation["theta"], quiet=1)
			self.mrcomegaslider.setValue(rotation["omega"], quiet=1)
		if comboboxidx == 4:
			if self._isRotNaN(rotation["ztilt"],rotation["ytilt"],rotation["xtilt"]): return
			self.xyzzslider.setValue(rotation["ztilt"], quiet=1)
			self.xyzyslider.setValue(rotation["ytilt"], quiet=1)
			self.xyzxslider.setValue(rotation["xtilt"], quiet=1)
		if comboboxidx == 5:
			if self._isRotNaN(rotation["n1"],rotation["n2"],rotation["n3"]): return
			if is_identity and self.spinn1slider.getValue() == 0.0 and self.spinn2slider.getValue() == 0.0 and self.spinn3slider.getValue() == 0.0:
				self.spinomegaslider .setValue(0.0, quiet=1)
				self.spinn1slider.setValue(0.0, quiet=1)
				self.spinn2slider.setValue(0.0, quiet=1)
				self.spinn3slider.setValue(1.0, quiet=1)				
			else:
				self.spinomegaslider .setValue(rotation["omega"], quiet=1)
				# Don't change slider if reult is Nan
				if rotation["n1"] == rotation["n1"]: self.spinn1slider.setValue(rotation["n1"], quiet=1)
				if rotation["n2"] == rotation["n2"]: self.spinn2slider.setValue(rotation["n2"], quiet=1)
				if rotation["n3"] == rotation["n3"]: self.spinn3slider.setValue(rotation["n3"], quiet=1)
		if comboboxidx == 6:
			if self._isRotNaN(rotation["n1"],rotation["n2"],rotation["n3"]): return
			if is_identity and self.spinn1slider.getValue() == 0.0 and self.spinn2slider.getValue() == 0.0 and self.spinn3slider.getValue() == 0.0:
				self.spinomegaslider.setValue(0.0, quiet=1)
				self.sgirotn1slider.setValue(0.0, quiet=1)
				self.sgirotn2slider.setValue(0.0, quiet=1)
				self.sgirotn3slider.setValue(1.0, quiet=1)
			else:
				self.spinomegaslider.setValue(rotation["q"], quiet=1)
				# Don't change slider if reult is Nan
				if rotation["n1"] == rotation["n1"]: self.sgirotn1slider.setValue(rotation["n1"], quiet=1)
				if rotation["n2"] == rotation["n2"]: self.sgirotn2slider.setValue(rotation["n2"], quiet=1)
				if rotation["n3"] == rotation["n3"]: self.sgirotn3slider.setValue(rotation["n3"], quiet=1)
		if comboboxidx == 7:
			if self._isRotNaN(rotation["e1"],rotation["e2"],rotation["e3"]): return
			if is_identity:
				self.quaternione0slider.setValue(1.0, quiet=1)
				self.quaternione1slider.setValue(0.0, quiet=1)
				self.quaternione2slider.setValue(0.0, quiet=1)
				self.quaternione3slider.setValue(0.0, quiet=1)
			else:	
				self.quaternione0slider.setValue(rotation["e0"], quiet=1)
				self.quaternione1slider.setValue(rotation["e1"], quiet=1)
				self.quaternione2slider.setValue(rotation["e2"], quiet=1)
				self.quaternione3slider.setValue(rotation["e3"], quiet=1)
		# Scaling update
		self.zoom.setValue(self.item3d().getTransform().get_scale())
		
	def updateMetaData(self):
		"""
		I didn't want to put this in update b/c this data doesn't change very often, and I don't want to waste CPU
		Its a judgement call really, less coupling vs. more efficiency
		"""
		if self.boundingbox: self.boundingbox.setText("Size: "+self.item3d().boundingboxsize)
			
		
	def addRotationWidgets(self):
		""" Add alll the widgets for the various EMAN2 rotation conventions """
		EMANwidget = QtGui.QWidget()
		Imagicwidget = QtGui.QWidget()
		Spiderwidget = QtGui.QWidget()
		MRCwidget = QtGui.QWidget()
		XYZwidget = QtGui.QWidget()
		spinwidget = QtGui.QWidget()
		sgirotwidget = QtGui.QWidget()
		quaternionwidget = QtGui.QWidget()
		# EMAN
		emanbox = QtGui.QVBoxLayout()
		self.emanazslider = ValSlider(EMANwidget, (0.0, 360.0), "  Az", rounding = 1)
		self.emanaltslider = ValSlider(EMANwidget, (0.0, 180.0), "Alt", rounding = 1)
		self.emanphislider = ValSlider(EMANwidget, (0.0, 360.0), "Phi", rounding = 1)
		emanbox.addWidget(self.emanazslider)
		emanbox.addWidget(self.emanaltslider)
		emanbox.addWidget(self.emanphislider)
		EMANwidget.setLayout(emanbox)
		# Imagic
		imagicbox = QtGui.QVBoxLayout()
		self.imagicgammaslider = ValSlider(Imagicwidget, (0.0, 360.0), "Gamma", rounding = 1)
		self.imagicbetaslider = ValSlider(Imagicwidget, (0.0, 180.0), "     Beta", rounding = 1)
		self.imagicalphaslider = ValSlider(Imagicwidget, (0.0, 360.0), "   Alpha", rounding = 1)
		imagicbox.addWidget(self.imagicgammaslider)
		imagicbox.addWidget(self.imagicbetaslider)
		imagicbox.addWidget(self.imagicalphaslider)
		Imagicwidget.setLayout(imagicbox)
		# Spider
		spiderbox = QtGui.QVBoxLayout()
		self.spiderpsislider = ValSlider(Spiderwidget, (0.0, 360.0), "   Psi", rounding = 1)
		self.spiderthetaslider = ValSlider(Spiderwidget, (0.0, 180.0), "Theta", rounding = 1)
		self.spiderphislider = ValSlider(Spiderwidget, (0.0, 360.0), "   Phi", rounding = 1)
		spiderbox.addWidget(self.spiderpsislider)
		spiderbox.addWidget(self.spiderthetaslider)
		spiderbox.addWidget(self.spiderphislider)
		Spiderwidget.setLayout(spiderbox)
		# MRC
		mrcbox = QtGui.QVBoxLayout()
		self.mrcpsislider = ValSlider(MRCwidget, (0.0, 360.0), "      Psi", rounding = 1)
		self.mrcthetaslider = ValSlider(MRCwidget, (0.0, 180.0), "  Theta", rounding = 1)
		self.mrcomegaslider = ValSlider(MRCwidget, (0.0, 360.0), "Omega", rounding = 1)
		mrcbox.addWidget(self.mrcpsislider)
		mrcbox.addWidget(self.mrcthetaslider)
		mrcbox.addWidget(self.mrcomegaslider)
		MRCwidget.setLayout(mrcbox)
		# XYZ
		xyzbox = QtGui.QVBoxLayout()
		self.xyzzslider = ValSlider(XYZwidget, (0.0, 360.0), "Z", rounding = 1)
		self.xyzyslider = ValSlider(XYZwidget, (0.0, 180.0), "Y", rounding = 1)
		self.xyzxslider = ValSlider(XYZwidget, (0.0, 360.0), "X", rounding = 1)
		xyzbox.addWidget(self.xyzzslider)
		xyzbox.addWidget(self.xyzyslider)
		xyzbox.addWidget(self.xyzxslider)
		XYZwidget.setLayout(xyzbox)
		# spin
		spinbox = QtGui.QVBoxLayout()
		self.spinomegaslider = ValSlider(spinwidget, (0.0, 180.0), "Omega", rounding = 1)
		self.spinn1slider = ValSlider(spinwidget, (0.0, 1.0), "       N1", rounding = 4)
		self.spinn2slider = ValSlider(spinwidget, (0.0, 1.0), "       N2", rounding = 4)
		self.spinn3slider = ValSlider(spinwidget, (0.0, 1.0), "       N3", rounding = 4)
		spinbox.addWidget(self.spinomegaslider)
		spinbox.addWidget(self.spinn1slider)
		spinbox.addWidget(self.spinn2slider)
		spinbox.addWidget(self.spinn3slider)
		spinwidget.setLayout(spinbox)
		# sgirot
		sgirotbox = QtGui.QVBoxLayout()
		self.sgirotqslider = ValSlider(sgirotwidget, (0.0, 180.0), " Q", rounding = 1)
		self.sgirotn1slider = ValSlider(sgirotwidget, (0.0, 1.0), "N1", rounding = 4)
		self.sgirotn2slider = ValSlider(sgirotwidget, (0.0, 1.0), "N2", rounding = 4)
		self.sgirotn3slider = ValSlider(sgirotwidget, (0.0, 1.0), "N3", rounding = 4)
		sgirotbox.addWidget(self.sgirotqslider)
		sgirotbox.addWidget(self.sgirotn1slider)
		sgirotbox.addWidget(self.sgirotn2slider)
		sgirotbox.addWidget(self.sgirotn3slider)
		sgirotwidget.setLayout(sgirotbox)
		# quaternion
		quaternionbox = QtGui.QVBoxLayout()
		self.quaternione0slider = ValSlider(quaternionwidget, (0.0, 1.0), "E0", rounding = 4)
		self.quaternione1slider = ValSlider(quaternionwidget, (0.0, 1.0), "E1", rounding = 4)
		self.quaternione2slider = ValSlider(quaternionwidget, (0.0, 1.0), "E2", rounding = 4)
		self.quaternione3slider = ValSlider(quaternionwidget, (0.0, 1.0), "E3", rounding = 4)
		quaternionbox.addWidget(self.quaternione0slider)
		quaternionbox.addWidget(self.quaternione1slider)
		quaternionbox.addWidget(self.quaternione2slider)
		quaternionbox.addWidget(self.quaternione3slider)
		quaternionwidget.setLayout(quaternionbox)
		# Add widgets to the stack
		self.rotstackedwidget.addWidget(EMANwidget)
		self.rotstackedwidget.addWidget(Imagicwidget)
		self.rotstackedwidget.addWidget(Spiderwidget)
		self.rotstackedwidget.addWidget(MRCwidget)
		self.rotstackedwidget.addWidget(XYZwidget)
		self.rotstackedwidget.addWidget(spinwidget)
		self.rotstackedwidget.addWidget(sgirotwidget)
		self.rotstackedwidget.addWidget(quaternionwidget)
		# add choices to combobox
		self.rotcombobox.addItem("EMAN")
		self.rotcombobox.addItem("Imagic")
		self.rotcombobox.addItem("Spider")
		self.rotcombobox.addItem("MRC")
		self.rotcombobox.addItem("XYZ")
		self.rotcombobox.addItem("spin")
		self.rotcombobox.addItem("sgirot")
		self.rotcombobox.addItem("quaternion")
	
		# Signal for all sliders
		QtCore.QObject.connect(self.rotcombobox, QtCore.SIGNAL("activated(int)"), self._rotcombobox_changed)
		QtCore.QObject.connect(self.emanazslider,QtCore.SIGNAL("valueChanged"),self._on_EMAN_rotation)
		QtCore.QObject.connect(self.emanaltslider,QtCore.SIGNAL("valueChanged"),self._on_EMAN_rotation)
		QtCore.QObject.connect(self.emanphislider,QtCore.SIGNAL("valueChanged"),self._on_EMAN_rotation)
		QtCore.QObject.connect(self.imagicgammaslider,QtCore.SIGNAL("valueChanged"),self._on_Imagic_rotation)
		QtCore.QObject.connect(self.imagicbetaslider,QtCore.SIGNAL("valueChanged"),self._on_Imagic_rotation)
		QtCore.QObject.connect(self.imagicalphaslider,QtCore.SIGNAL("valueChanged"),self._on_Imagic_rotation)
		QtCore.QObject.connect(self.spiderpsislider,QtCore.SIGNAL("valueChanged"),self._on_Spider_rotation)
		QtCore.QObject.connect(self.spiderthetaslider,QtCore.SIGNAL("valueChanged"),self._on_Spider_rotation)
		QtCore.QObject.connect(self.spiderphislider,QtCore.SIGNAL("valueChanged"),self._on_Spider_rotation)
		QtCore.QObject.connect(self.mrcpsislider,QtCore.SIGNAL("valueChanged"),self._on_MRC_rotation)
		QtCore.QObject.connect(self.mrcthetaslider,QtCore.SIGNAL("valueChanged"),self._on_MRC_rotation)
		QtCore.QObject.connect(self.mrcomegaslider,QtCore.SIGNAL("valueChanged"),self._on_MRC_rotation)
		QtCore.QObject.connect(self.xyzzslider,QtCore.SIGNAL("valueChanged"),self._on_XYZ_rotation)
		QtCore.QObject.connect(self.xyzyslider,QtCore.SIGNAL("valueChanged"),self._on_XYZ_rotation)
		QtCore.QObject.connect(self.xyzxslider,QtCore.SIGNAL("valueChanged"),self._on_XYZ_rotation)
		QtCore.QObject.connect(self.spinomegaslider,QtCore.SIGNAL("valueChanged"),self._on_spin_rotation)
		QtCore.QObject.connect(self.spinn1slider,QtCore.SIGNAL("valueChanged"),self._on_spin_rotation)
		QtCore.QObject.connect(self.spinn2slider,QtCore.SIGNAL("valueChanged"),self._on_spin_rotation)
		QtCore.QObject.connect(self.spinn3slider,QtCore.SIGNAL("valueChanged"),self._on_spin_rotation)
		QtCore.QObject.connect(self.sgirotqslider,QtCore.SIGNAL("valueChanged"),self._on_sgirot_rotation)
		QtCore.QObject.connect(self.sgirotn1slider,QtCore.SIGNAL("valueChanged"),self._on_sgirot_rotation)
		QtCore.QObject.connect(self.sgirotn2slider,QtCore.SIGNAL("valueChanged"),self._on_sgirot_rotation)
		QtCore.QObject.connect(self.sgirotn3slider,QtCore.SIGNAL("valueChanged"),self._on_sgirot_rotation)
		QtCore.QObject.connect(self.quaternione0slider,QtCore.SIGNAL("valueChanged"),self._on_quaternion_rotation)
		QtCore.QObject.connect(self.quaternione1slider,QtCore.SIGNAL("valueChanged"),self._on_quaternion_rotation)
		QtCore.QObject.connect(self.quaternione2slider,QtCore.SIGNAL("valueChanged"),self._on_quaternion_rotation)
		QtCore.QObject.connect(self.quaternione3slider,QtCore.SIGNAL("valueChanged"),self._on_quaternion_rotation)	
		
	def _rotcombobox_changed(self, idx):
		self.rotstackedwidget.setCurrentIndex(idx)
		self.updateItemControls()
		
	def _on_EMAN_rotation(self, value):
		self._set_rotation_std_coords(Transform({"type":"eman","az":self.emanazslider.getValue(),"alt":self.emanaltslider.getValue(),"phi":self.emanphislider.getValue()}))
		self.inspector().updateSceneGraph()
		
	def _on_Imagic_rotation(self, value):
		self._set_rotation_std_coords(Transform({"type":"imagic","gamma":self.imagicgammaslider.getValue(),"beta":self.imagicbetaslider.getValue(),"alpha":self.imagicalphaslider.getValue()}))
		self.inspector().updateSceneGraph()
		
	def _on_Spider_rotation(self, value):
		self._set_rotation_std_coords(Transform({"type":"spider","psi":self.spiderpsislider.getValue(),"theta":self.spiderthetaslider.getValue(),"phi":self.spiderphislider.getValue()}))
		self.inspector().updateSceneGraph()
		
	def _on_MRC_rotation(self, value):
		self._set_rotation_std_coords(Transform({"type":"mrc","phi":self.mrcpsislider.getValue(),"theta":self.mrcthetaslider.getValue(),"omega":self.mrcomegaslider.getValue()}))
		self.inspector().updateSceneGraph()
		
	def _on_XYZ_rotation(self, value):
		self._set_rotation_std_coords(Transform({"type":"xyz","ztilt":self.xyzzslider.getValue(),"ytilt":self.xyzyslider.getValue(),"xtilt":self.xyzxslider.getValue()}))
		self.inspector().updateSceneGraph()
		
	def _on_spin_rotation(self, value):
		v = Vec3f(self.spinn1slider.getValue(), self.spinn2slider.getValue(), self.spinn3slider.getValue())
		v.normalize()
		self._set_rotation_std_coords(Transform({"type":"spin","omega":self.spinomegaslider.getValue(),"n1":v[0],"n2":v[1],"n3":v[2]}))
		self.inspector().updateSceneGraph()
		
	def _on_sgirot_rotation(self, value):
		v = Vec3f(self.sgirotn1slider.getValue(), self.sgirotn2slider.getValue(), self.sgirotn3slider.getValue())
		v.normalize()
		self._set_rotation_std_coords(Transform({"type":"sgirot","q":self.sgirotqslider.getValue(),"n1":v[0],"n2":v[1],"n3":v[2]}))
		self.inspector().updateSceneGraph()
		
	def _on_quaternion_rotation(self, value):
		v = Vec4f(self.quaternione0slider.getValue(), self.quaternione1slider.getValue(), self.quaternione2slider.getValue(), self.quaternione3slider.getValue())
		v.normalize()
		self._set_rotation_std_coords(Transform({"type":"quaternion","e0":v[0],"e1":v[1],"e2":v[2],"e3":v[3]}))
		self.inspector().updateSceneGraph()
		
	def _set_rotation_std_coords(self, rotation):
		""" This function sets the rotation as if there were no preceeding ones, otherwise a rot around Z could be arounf y,x, etc.
		Works by transforming local coords into global corrds"""
		tt = rotation
		tp = self.item3d().getParentMatrixProduct()
		if tp: tt = tp.inverse()*rotation
		self.item3d().getTransform().set_rotation(tt.get_rotation())
		
if __name__ == '__main__':
	#Test code
	root = EMItem3D(0)
	a = EMItem3D(root)
	b = EMItem3D(root)
	c = EMItem3D(root)
#	root.addChildren([a,b,c])
	aa = EMItem3D(a)
	ab = EMItem3D(a)
	ac = EMItem3D(a)
#	a.addChildren([aa,ab,ac])
	ba = EMItem3D(b)
	bb = EMItem3D(b)
	bc = EMItem3D(b)
#	b.addChildren([ba,bb,bc])
	ca = EMItem3D(c)
	cb = EMItem3D(c)
	cc = EMItem3D(c)
#	c.addChildren([ca,cb,cc])
	
	aaa = EMItem3D(aa)
	aab = EMItem3D(aa)
	aac = EMItem3D(aa)
#	aa.addChildren([aaa,aab,aac])
	
	a.is_selected = True
	aab.is_selected = True
	ba.is_selected = True
	bc.is_selected = True
	c.is_selected = True
	cc.is_selected = True
	
	print "getAllSelectedNodes() test: "
	print "\tpassed?  ", set(root.getAllSelectedNodes()) == set([a,aab,ba,bc,c,cc])
	print "getNearbySelectedNodes() test: "
	print "\tpassed?  ", set(root.getNearbySelectedNodes()) == set([a,ba,bc,c])
	print "getDistantSelectedNodes() test: "
	print "\tpassed?  ", set(root.getDistantSelectedNodes()) == set([aab,ba,bc,cc])
	
	root.displayTree()
	print "\n"
	a.removeParentReferences()
	root.displayTree()
	print "\n"
	root.addParentReferences()
	root.displayTree()
	print "\n"
	
	print "removing child..."
	root.removeChild(a)
	root.displayTree()
	del a, b, c, aa, ab, ac, ba, bb, bc, ca, cb, cc, aaa, aab, aac #Ensure instances are deleted before the class object is deleted, which is important in EMItem3D.__del__(self):