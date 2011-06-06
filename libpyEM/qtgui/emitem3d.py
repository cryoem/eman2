from OpenGL import GL
from EMAN2 import Transform
from libpyGLUtils2 import GLUtil
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
		#NOTE: Accessor methods are not needed; Python "properties" can be used instead if needed.
		self.parent = parent		
		self.children = set(children)
		self.transform = transform
		self.is_visible = True 
		self.is_selected = False
		self.widget = None
		self.boundingboxsize = None
		self.getnset_unique_integer()
		
	def __del__(self):
		EMItem3D.selection_recycle.append(self.intname)
	
	def getnset_unique_integer(self):
		"""
		Stuff for the selection mechanism, return a unique int for each instance of EMItem3D
		"""
		if len(EMItem3D.selection_recycle) > 0:
			self.intname = EMItem3D.selection_recycle.pop()
		else:
			EMItem3D.selection_intname += 1
			self.intname = EMItem3D.selection_intname
		EMItem3D.selection_idx_dict[self.intname] = weakref.ref(self)
		
	def add_child(self, node):
		"""
		Adds a child node, if not already in the set of child nodes.
		@type node: EMItem3D
		@param node: the child node to add
		"""
		self.children.add(node)
		node.parent = self
	
	def add_children(self, nodes):
		"""
		Adds all the provided child nodes which are not already in the set of child nodes.
		@type nodes: an iterable collection of EMItem3D (or subclass) objects
		@param nodes: the nodes which will be added as child nodes
		"""
		for node in nodes:
			self.children.add(node)
			node.parent = self
			
	def has_child(self, node):
		"""
		Tests whether the supplied node is a child of the current node. 
		@type node: EMItem3D
		@param node: test whether this node is a child node of self 
		"""
		return node in self.children
		
	def remove_child(self, node):
		"""
		Remove the supplied node from the set of child nodes. 
		@type node: EMItem3D
		@param node: the node to remove
		"""
		self.children.remove(node)
		node.parent = None
		
	def get_all_selected_nodes(self):
		"""
		For the tree rooted at self, this recursive method returns a list of all the selected nodes.
		@return: a list of selected nodes
		"""
		selected_list = []
		if self.is_selected:
			selected_list.append(self)
		for child in self.children: #Recursion ends on leaf nodes here
			selected_list.extend(child.get_all_selected_nodes()) #Recursion
		
		return selected_list
	
	def get_nearest_selected_nodes(self):
		"""
		For the tree rooted at self, this recursive method returns a list of the selected nodes that are nearest to self.
		A selected node will not be in the returned list if one of its ancestor nodes is also selected. 
		@return: a list of selected nodes
		"""
		selected_list = []
		if self.is_selected:
			return [self]
		else:
			for child in self.children:
				selected_list.extend(child.get_nearest_selected_nodes())
		
		return selected_list
	
	def get_farthest_selected_nodes(self):
		"""
		For the tree rooted at self, this recursive method returns a list of the selected nodes that are farthest from self.
		A selected node will not be in the returned list if one of its descendant nodes is also selected. 
		@return: a list of selected nodes
		"""
		selected_list = []
		for child in self.children:
			selected_list.extend(child.get_farthest_selected_nodes())
		if not selected_list: #either this is a leaf node, or there are no selected nodes in the child subtrees
			if self.is_selected:
				selected_list.append(self)
		
		return selected_list

	def get_scene_gui(self):
		"""
		Return a Qt widget that controls the scene item
		"""
		raise NotImplementedError("GUI controls must be implented in the subclasses")
		
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
		It calls self.render_node() to render the current node. 
		Usually, this method is unchanged in subclasses. 
		"""
		if not self.is_visible:
			return #Also applies to subtree rooted at this node
		
		if self.transform != None:
			GL.glPushMatrix()
			GL.glPushName(self.intname)
			GLUtil.glMultMatrix(self.transform) #apply the transformation
			
			self.render_node()
			for child in self.children:
				child.render()
			GL.glPopName()
			GL.glPopMatrix()
			
		else:
			GL.glPushName(self.intname)
			self.render_node()
			for child in self.children:
				child.render()
			GL.glPopName()
		

	def render_node(self):
		"""
		This method, which is called by self.render(), renders the current node.
		It should be implemented in subclasses that represent visible objects.
		"""
		pass

	def get_scene_gui(self):
		"""
		Return a Qt widget that controls the scene item
		"""
		return self.widget

	def keyPressEvent(self, event): pass
	def keyReleaseEvent(self, event): pass
	def mouseDoubleClickEvent(self, event): pass
	def mouseMoveEvent(self, event): pass
	def mousePressEvent(self, event): pass
	def mouseReleaseEvent(self, event): pass
	def wheelEvent(self, event): pass

if __name__ == '__main__':
	#Test code
	root = EMItem3D(0)
	a = EMItem3D(1,root)
	b = EMItem3D(2,root)
	c = EMItem3D(3,root)
	root.add_children([a,b,c])
	aa = EMItem3D(4)
	ab = EMItem3D(5)
	ac = EMItem3D(6)
	a.add_children([aa,ab,ac])
	ba = EMItem3D(7)
	bb = EMItem3D(8)
	bc = EMItem3D(9)
	b.add_children([ba,bb,bc])
	ca = EMItem3D(10)
	cb = EMItem3D(11)
	cc = EMItem3D(12)
	c.add_children([ca,cb,cc])
	
	aaa = EMItem3D(13)
	aab = EMItem3D(14)
	aac = EMItem3D(15)
	aa.add_children([aaa,aab,aac])
	
	a.is_selected = True
	aab.is_selected = True
	ba.is_selected = True
	bc.is_selected = True
	c.is_selected = True
	cc.is_selected = True
	
	print "get_all_selected_nodes() test: "
	print "\tpassed?  ", set(root.get_all_selected_nodes()) == set([a,aab,ba,bc,c,cc])
	print "get_nearest_selected_nodes() test: "
	print "\tpassed?  ", set(root.get_nearest_selected_nodes()) == set([a,ba,bc,c])
	print "get_farthest_selected_nodes() test: "
	print "\tpassed?  ", set(root.get_farthest_selected_nodes()) == set([aab,ba,bc,cc])