#!/usr/bin/env python
#
# Author: Ross Coleman (racolema@bcm.edu)
# Copyright (c) 2000-2011 Baylor College of Medicine

from OpenGL import GL
from PyQt4 import QtCore, QtGui
from EMAN2 import Transform
from libpyGLUtils2 import GLUtil
from valslider import ValSlider, EMSpinWidget
import weakref


class EMItem3D(object): #inherit object for new-style class (new-stype classes required for super() and Python properties)
	"""
	The base class for nodes in our scene graph, which is used in our 3D scene viewer.
	In our case, the scene graph is a tree data structure.
	"""
	# Class attrib to connect openGL int identifiers to class instances
	selection_idx_dict = {}
	selection_recycle = []
	selection_intname = -1
	
	def __init__(self, parent = None, children = set(), transform = None):
		"""
		@type parent: EMItem3D
		@param parent: the parent node to the current node or None for the root node
		@type children: set type (list or tuple will be converted to a set)
		@param children: the child nodes
		@type transform: Transform or None
		@param transform: The transformation (rotation, scaling, translation) that should be applied before rendering this node and its children 
		"""
		self.parent = parent		
		self.children = set(children)
		self.transform = transform
		self.is_visible = True 
		self.is_selected = False
		self.widget = None			# This is an inspector widget
		self.EMQTreeWidgetItem = None 		# This is an inspector tree item
		self.boundingboxsize = None
		self.getAndSetUniqueInteger()
	
	def getChildren(self): return self.children
	def setChildren(self, children): self.children = set(children)
	def getParent(self): return self.parent
	def setParent(self, parent): self.parent = parent
	def isSelectedItem(self): return self.is_selected
	def setSelectedItem(self, is_selected): self.is_selected = is_selected
	def getTransform(self): return self.transform
	def setTransform(self, transform): self.transform = transform
	def isVisibleItem(self): return self.is_visible
	def setVisibleItem(self, is_visible): self.is_visible = is_visible

	
	def __del__(self):
		EMItem3D.selection_recycle.append(self.intname)
	
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
		Adds a child node, if not already in the set of child nodes.
		@type node: EMItem3D
		@param node: the child node to add
		"""
		self.children.add(node)
		node.parent = self
	
	def addChildren(self, nodes):
		"""
		Adds all the provided child nodes which are not already in the set of child nodes.
		@type nodes: an iterable collection of EMItem3D (or subclass) objects
		@param nodes: the nodes which will be added as child nodes
		"""
		for node in nodes:
			self.children.add(node)
			node.parent = self
			
	def hasChild(self, node):
		"""
		Tests whether the supplied node is a child of the current node. 
		@type node: EMItem3D
		@param node: test whether this node is a child node of self 
		"""
		return node in self.children
		
	def removeChild(self, node):
		"""
		Remove the supplied node from the set of child nodes. 
		@type node: EMItem3D
		@param node: the node to remove
		"""
		self.children.remove(node)
		node.parent = None
		
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

	def getSceneGui(self):
		"""
		Return a Qt widget that controls the scene item
		"""
		raise NotImplementedError("GUI controls must be implemented in the subclasses")
	
	def setEMQTreeWidgetItem(self, node):
		"""
		Relate a QtreeItem to this node
		"""
		self.EMQTreeWidgetItem = node
		
	def update_matrices(self, params, xformtype):
		"""
		@type params: List
		@param params: A list defining how the transform in each active node is modified
		@type xfromtype: sting
		@param xformtype: The sort of transform we wish to do
		"""
		if self.is_selected:
			if xformtype == "rotate":
				self.transform.rotate_origin(Transform({"type":"spin","Omega":params[0],"n1":params[1],"n2":params[2],"n3":params[3]}))
			elif xformtype == "translate":
				self.transform.translate(params[0], params[1], params[2])
			elif xformtype == "scale":
				self.transform.scale(params[0])
			else:
				raise Exception,"Invalid transformation type"
		
		# Now tell all children to update
		for child in self.children:
			child.update_matrices(params, xformtype)
			
	def render(self):
		"""
		This is the method to call to render the node and its child nodes. 
		It calls self.renderNode() to render the current node. 
		Usually, this method is unchanged in subclasses. 
		"""
		if not self.is_visible:
			return #Also applies to subtree rooted at this node
		
		if self.transform != None:
			if self.widget != None and self.is_selected: self.widget.updateItemControls()
			GL.glPushMatrix()
			GL.glPushName(self.intname)
			GLUtil.glMultMatrix(self.transform) #apply the transformation
			
			self.renderNode()
			for child in self.children:
				child.render()
			GL.glPopName()
			GL.glPopMatrix()
			
		else:
			GL.glPushName(self.intname)
			self.renderNode()
			for child in self.children:
				child.render()
			GL.glPopName()
		

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

        
class EMInspectorControlBasic(QtGui.QWidget):
    """
    Class to make the EMItem GUI controls
    """
    def __init__(self, name, sgnode):
        QtGui.QWidget.__init__(self)
        self.sgnode = sgnode
        self.name = name
        self.inspector = None
        self.transfromboxmaxheight = 400
        
        igvbox = QtGui.QVBoxLayout()
        self.addBasicControls(igvbox)
        self.addColorControls(igvbox)
        self.addControls(igvbox)
        self.setLayout(igvbox)
    
    def setInspector(self, inspector):
        self.inspector = inspector
        
    def addBasicControls(self, igvbox):
        # selection box and label
        font = QtGui.QFont()
        font.setBold(True)
        label = QtGui.QLabel(self.name,self)
        label.setFont(font)
        label.setAlignment(QtCore.Qt.AlignCenter)
        igvbox.addWidget(label)
        databox = QtGui.QHBoxLayout()
        if self.sgnode.boundingboxsize:
            databox.addWidget(QtGui.QLabel("Size: "+str(self.sgnode.boundingboxsize)+u'\u00B3',self))
        igvbox.addLayout(databox)
        # angluar controls
        xformframe = QtGui.QFrame()
        xformframe.setFrameShape(QtGui.QFrame.StyledPanel)
        xformbox = QtGui.QVBoxLayout()
        xformlabel = QtGui.QLabel("Transformation", xformframe)
        xformlabel.setFont(font)
        xformlabel.setAlignment(QtCore.Qt.AlignCenter)
        xformbox.addWidget(xformlabel)
        # Rotations
        self.rotcombobox = QtGui.QComboBox()
        xformbox.addWidget(self.rotcombobox)
        self.rotstackedwidget = QtGui.QStackedWidget()
        self.addRotationWidgets()
        xformbox.addWidget(self.rotstackedwidget)
        #translations
        textbox = QtGui.QHBoxLayout()
        txlabel = QtGui.QLabel("TX",xformframe)
        txlabel.setAlignment(QtCore.Qt.AlignCenter)
        textbox.addWidget(txlabel)
        tylabel = QtGui.QLabel("TY",xformframe)
        tylabel.setAlignment(QtCore.Qt.AlignCenter)
        textbox.addWidget(tylabel)
        xformbox.addLayout(textbox)
        box = QtGui.QHBoxLayout()
        self.tx = EMSpinWidget(0.0, 1.0)
        self.ty = EMSpinWidget(0.0, 1.0)
        box.addWidget(self.tx)
        box.addWidget(self.ty)
        xformbox.addLayout(box)
        zoombox = QtGui.QHBoxLayout()
        tzlabel = QtGui.QLabel("TZ",xformframe)
        tzlabel.setAlignment(QtCore.Qt.AlignCenter)
        zoombox.addWidget(tzlabel)
        zoomlabel = QtGui.QLabel("Zoom",xformframe)
        zoomlabel.setAlignment(QtCore.Qt.AlignCenter)
        zoombox.addWidget(zoomlabel)
        xformbox.addLayout(zoombox)
        zoomwidgetbox = QtGui.QHBoxLayout()
        self.tz = EMSpinWidget(0.0, 1.0)
        self.zoom = EMSpinWidget(0.0, 0.1)
        zoomwidgetbox.addWidget(self.tz)
        zoomwidgetbox.addWidget(self.zoom)
        xformbox.addLayout(zoomwidgetbox)
                
        xformframe.setMaximumHeight(self.transfromboxmaxheight)
        xformframe.setLayout(xformbox)
        igvbox.addWidget(xformframe)
        
        QtCore.QObject.connect(self.tx,QtCore.SIGNAL("valueChanged(int)"),self._on_translation)
        QtCore.QObject.connect(self.ty,QtCore.SIGNAL("valueChanged(int)"),self._on_translation)
        QtCore.QObject.connect(self.tz,QtCore.SIGNAL("valueChanged(int)"),self._on_translation)
        QtCore.QObject.connect(self.zoom,QtCore.SIGNAL("valueChanged(int)"),self._on_scale)
    
    def _on_translation(self, value):
        self.sgnode.getTransform().set_trans(self.tx.getValue(), self.ty.getValue(), self.tz.getValue())
        self.inspector.updateSceneGraph()
        
    def _on_scale(self, value):
        self.sgnode.getTransform().set_scale(self.zoom.getValue())
        self.inspector.updateSceneGraph()
        
    def addColorControls(self, igvbox):
        pass
    
    def addControls(self, igvbox):
        pass
    
    def updateItemControls(self):
        # Translation update
        translation =  self.sgnode.getTransform().get_trans()
        self.tx.setValue(translation[0])
        self.ty.setValue(translation[1])
        self.tz.setValue(translation[2])
        # Rotation update
        rotation =  self.sgnode.getTransform().get_rotation(str(self.rotcombobox.currentText()))
        comboboxidx = self.rotcombobox.currentIndex()
        if comboboxidx == 0:
            self.emanazslider.setValue(rotation["az"], quiet=1)
            self.emanaltslider.setValue(rotation["alt"], quiet=1)
            self.emanphislider.setValue(rotation["phi"], quiet=1)
        if comboboxidx == 1:
            self.imagicgammaslider.setValue(rotation["gamma"], quiet=1)
            self.imagicbetaslider.setValue(rotation["beta"], quiet=1)
            self.imagicalphaslider.setValue(rotation["alpha"], quiet=1)
        if comboboxidx == 2:
            self.spiderpsislider.setValue(rotation["psi"], quiet=1)
            self.spiderthetaslider.setValue(rotation["theta"], quiet=1)
            self.spiderphislider.setValue(rotation["phi"], quiet=1)
        if comboboxidx == 3:
            self.mrcpsislider.setValue(rotation["phi"], quiet=1)
            self.mrcthetaslider.setValue(rotation["theta"], quiet=1)
            self.mrcomegaslider.setValue(rotation["omega"], quiet=1)
        if comboboxidx == 4:
            self.xyzzslider.setValue(rotation["ztilt"], quiet=1)
            self.xyzyslider.setValue(rotation["ytilt"], quiet=1)
            self.xyzxslider.setValue(rotation["xtilt"], quiet=1)
        if comboboxidx == 5:
            self.spinomegaslider .setValue(rotation["Omega"], quiet=1)
            self.spinn1slider.setValue(rotation["n1"], quiet=1)
            self.spinn2slider.setValue(rotation["n2"], quiet=1)
            self.spinn3slider.setValue(rotation["n3"], quiet=1)
        if comboboxidx == 6:
            self.spinomegaslider.setValue(rotation["q"], quiet=1)
            self.sgirotn1slider.setValue(rotation["n1"], quiet=1)
            self.sgirotn2slider.setValue(rotation["n2"], quiet=1)
            self.sgirotn3slider.setValue(rotation["n3"], quiet=1)
        if comboboxidx == 7:
            self.quaternione0slider.setValue(rotation["e0"], quiet=1)
            self.quaternione1slider.setValue(rotation["e1"], quiet=1)
            self.quaternione2slider.setValue(rotation["e2"], quiet=1)
            self.quaternione3slider.setValue(rotation["e3"], quiet=1)
        # Scaling update
        self.zoom.setValue(self.sgnode.getTransform().get_scale())
        
    def addRotationWidgets(self):
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
        self.emanazslider = ValSlider(EMANwidget, (0.0, 360.0), "  Az")
        self.emanaltslider = ValSlider(EMANwidget, (0.0, 180.0), "Alt")
        self.emanphislider = ValSlider(EMANwidget, (0.0, 360.0), "Phi")
        emanbox.addWidget(self.emanazslider)
        emanbox.addWidget(self.emanaltslider)
        emanbox.addWidget(self.emanphislider)
        EMANwidget.setLayout(emanbox)
        # Imagic
        imagicbox = QtGui.QVBoxLayout()
        self.imagicgammaslider = ValSlider(Imagicwidget, (0.0, 360.0), "Gamma")
        self.imagicbetaslider = ValSlider(Imagicwidget, (0.0, 180.0), "     Beta")
        self.imagicalphaslider = ValSlider(Imagicwidget, (0.0, 360.0), "   Alpha")
        imagicbox.addWidget(self.imagicgammaslider)
        imagicbox.addWidget(self.imagicbetaslider)
        imagicbox.addWidget(self.imagicalphaslider)
        Imagicwidget.setLayout(imagicbox)
        # Spider
        spiderbox = QtGui.QVBoxLayout()
        self.spiderpsislider = ValSlider(Spiderwidget, (0.0, 360.0), "   Psi")
        self.spiderthetaslider = ValSlider(Spiderwidget, (0.0, 180.0), "Theta")
        self.spiderphislider = ValSlider(Spiderwidget, (0.0, 360.0), "   Phi")
        spiderbox.addWidget(self.spiderpsislider)
        spiderbox.addWidget(self.spiderthetaslider)
        spiderbox.addWidget(self.spiderphislider)
        Spiderwidget.setLayout(spiderbox)
        # MRC
        mrcbox = QtGui.QVBoxLayout()
        self.mrcpsislider = ValSlider(MRCwidget, (0.0, 360.0), "      Psi")
        self.mrcthetaslider = ValSlider(MRCwidget, (0.0, 180.0), "  Theta")
        self.mrcomegaslider = ValSlider(MRCwidget, (0.0, 360.0), "Omega")
        mrcbox.addWidget(self.mrcpsislider)
        mrcbox.addWidget(self.mrcthetaslider)
        mrcbox.addWidget(self.mrcomegaslider)
        MRCwidget.setLayout(mrcbox)
        # XYZ
        xyzbox = QtGui.QVBoxLayout()
        self.xyzzslider = ValSlider(XYZwidget, (0.0, 360.0), "Z")
        self.xyzyslider = ValSlider(XYZwidget, (0.0, 180.0), "Y")
        self.xyzxslider = ValSlider(XYZwidget, (0.0, 360.0), "X")
        xyzbox.addWidget(self.xyzzslider)
        xyzbox.addWidget(self.xyzyslider)
        xyzbox.addWidget(self.xyzxslider)
        XYZwidget.setLayout(xyzbox)
        # spin
        spinbox = QtGui.QVBoxLayout()
        self.spinomegaslider = ValSlider(spinwidget, (0.0, 360.0), "Omega")
        self.spinn1slider = ValSlider(spinwidget, (0.0, 1.0), "       N1")
        self.spinn2slider = ValSlider(spinwidget, (0.0, 1.0), "       N2")
        self.spinn3slider = ValSlider(spinwidget, (0.0, 1.0), "       N3")
        spinbox.addWidget(self.spinomegaslider)
        spinbox.addWidget(self.spinn1slider)
        spinbox.addWidget(self.spinn2slider)
        spinbox.addWidget(self.spinn3slider)
        spinwidget.setLayout(spinbox)
        # sgirot
        sgirotbox = QtGui.QVBoxLayout()
        self.sgirotqslider = ValSlider(sgirotwidget, (0.0, 360.0), " Q")
        self.sgirotn1slider = ValSlider(sgirotwidget, (0.0, 1.0), "N1")
        self.sgirotn2slider = ValSlider(sgirotwidget, (0.0, 1.0), "N2")
        self.sgirotn3slider = ValSlider(sgirotwidget, (0.0, 1.0), "N3")
        sgirotbox.addWidget(self.sgirotqslider)
        sgirotbox.addWidget(self.sgirotn1slider)
        sgirotbox.addWidget(self.sgirotn2slider)
        sgirotbox.addWidget(self.sgirotn3slider)
        sgirotwidget.setLayout(sgirotbox)
        # quaternion
        quaternionbox = QtGui.QVBoxLayout()
        self.quaternione0slider = ValSlider(quaternionwidget, (0.0, 1.0), "E0")
        self.quaternione1slider = ValSlider(quaternionwidget, (0.0, 1.0), "E1")
        self.quaternione2slider = ValSlider(quaternionwidget, (0.0, 1.0), "E2")
        self.quaternione3slider = ValSlider(quaternionwidget, (0.0, 1.0), "E3")
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
        self.sgnode.getTransform().set_rotation({"type":"eman","az":self.emanazslider.getValue(),"alt":self.emanaltslider.getValue(),"phi":self.emanphislider.getValue()})
        self.inspector.updateSceneGraph()
        
    def _on_Imagic_rotation(self, value):
        self.sgnode.getTransform().set_rotation({"type":"imagic","gamma":self.imagicgammaslider.getValue(),"beta":self.imagicbetaslider.getValue(),"alpha":self.imagicalphaslider.getValue()})
        self.inspector.updateSceneGraph()
        
    def _on_Spider_rotation(self, value):
        self.sgnode.getTransform().set_rotation({"type":"spider","psi":self.spiderpsislider.getValue(),"theta":self.spiderthetaslider.getValue(),"phi":self.spiderphislider.getValue()})
        self.inspector.updateSceneGraph()
        
    def _on_MRC_rotation(self, value):
        self.sgnode.getTransform().set_rotation({"type":"mrc","phi":self.mrcpsislider.getValue(),"theta":self.mrcthetaslider.getValue(),"omega":self.mrcomegaslider.getValue()})
        self.inspector.updateSceneGraph()
        
    def _on_XYZ_rotation(self, value):
        self.sgnode.getTransform().set_rotation({"type":"xyz","ztilt":self.xyzzslider.getValue(),"ytilt":self.xyzyslider.getValue(),"xtilt":self.xyzxslider.getValue()})
        self.inspector.updateSceneGraph()
        
    def _on_spin_rotation(self, value):
        self.sgnode.getTransform().set_rotation({"type":"spin","Omega":self.spinomegaslider.getValue(),"n1":self.spinn1slider.getValue(),"n2":self.spinn2slider.getValue(),"n3":self.spinn3slider.getValue()})
        self.inspector.updateSceneGraph()
        
    def _on_sgirot_rotation(self, value):
        self.sgnode.getTransform().set_rotation({"type":"sgirot","q":self.sgirotqslider.getValue(),"n1":self.sgirotn1slider.getValue(),"n2":self.sgirotn2slider.getValue(),"n3":self.sgirotn3slider.getValue()})
        self.inspector.updateSceneGraph()
        
    def _on_quaternion_rotation(self, value):
        self.sgnode.getTransform().set_rotation({"type":"quaternion","e0":self.quaternione0slider.getValue(),"e1":self.quaternione1slider.getValue(),"e2":self.quaternione2slider.getValue(),"e3":self.quaternione3slider.getValue()})
        self.inspector.updateSceneGraph()
        




if __name__ == '__main__':
	#Test code
	root = EMItem3D(0)
	a = EMItem3D(root)
	b = EMItem3D(root)
	c = EMItem3D(root)
	root.addChildren([a,b,c])
	aa = EMItem3D()
	ab = EMItem3D()
	ac = EMItem3D()
	a.addChildren([aa,ab,ac])
	ba = EMItem3D()
	bb = EMItem3D()
	bc = EMItem3D()
	b.addChildren([ba,bb,bc])
	ca = EMItem3D()
	cb = EMItem3D()
	cc = EMItem3D()
	c.addChildren([ca,cb,cc])
	
	aaa = EMItem3D()
	aab = EMItem3D()
	aac = EMItem3D()
	aa.addChildren([aaa,aab,aac])
	
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