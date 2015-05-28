#!/usr/bin/env python

#
# Author: James Michael Bell, 2016 (jmbell@bcm.edu)
# Bond Drawing Code By: Muthu Alagappan, m.alagappan901@gmail.com, 07/22/09
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

from EMAN2 import *
from emglobjects import get_default_gl_colors
from emitem3d import EMItem3D, EMItem3DInspector, drawBoundingBox
from libpyGLUtils2 import GLUtil
import os
import sys

from OpenGL.GL import *
from OpenGL.GLU import *
from PyQt4 import QtCore, QtGui

import numpy as np

class EMPDBItem3D(EMItem3D):
	"""
	This class is the scene graph node that has a reference to a Protein Data Bank (PDB) structure. 
	Though it has an orientation, it can not be displayed directly. Instead, its children are 
	displayable, and are sized, postioned, and oriented relative to this node.
	"""
	
	name = "PDB"
	nodetype = "SceneChild"
	
	@staticmethod
	def getNodeDialogWidget(attribdict):
		"""Get PDB Widget"""
		pdbwidget = QtGui.QWidget()
		grid = QtGui.QGridLayout()
		node_name_data_label = QtGui.QLabel("PDB Model Label")
		attribdict["node_name"] = QtGui.QLineEdit()
		data_path_label = QtGui.QLabel("PDB Model Path")
		attribdict["data_path"] = QtGui.QLineEdit()
		browse_button = QtGui.QPushButton("Browse")
		grid.addWidget(node_name_data_label, 0, 0, 1, 2)
		grid.addWidget(attribdict["node_name"], 0, 2, 1, 2)
		grid.addWidget(data_path_label, 1, 0, 1, 2)
		grid.addWidget(attribdict["data_path"], 1, 2, 1, 2)
		grid.addWidget(browse_button, 2, 0, 1, 4)
		EMItem3D.get_transformlayout(grid, 4, attribdict)
		pdbwidget.setLayout(grid)
		EMPDBItem3D.attribdict = attribdict
		QtCore.QObject.connect(browse_button, QtCore.SIGNAL('clicked()'), EMPDBItem3D._on_browse)
		return pdbwidget
	
	@staticmethod
	def _on_browse():
		filename = QtGui.QFileDialog.getOpenFileName(None, 'Get file', os.getcwd())
		if filename:
			EMPDBItem3D.attribdict["data_path"].setText(filename)
			#name = os.path.basename(str(filename))
			name = str(filename).split('/')[-1].split('.')[0]
			EMPDBItem3D.attribdict["node_name"].setText(str(name))

	@staticmethod
	def getNodeForDialog(attribdict):
		"""Create a new node using a attribdict"""
		return EMPDBItem3D(str(attribdict["data_path"].text()), transform=EMItem3D.getTransformFromDict(attribdict))
	
	def __init__(self, pdb_file, parent=None, children = set(), transform=None, style='bs'):
		#if not transform: transform = Transform()	# Object initialization should not be put in the constructor. Causes issues
		EMItem3D.__init__(self, parent, children, transform=transform)
		self.setData(pdb_file)
		self.diffuse = [0.5,0.5,0.5,1.0]
		self.specular = [1.0,1.0,1.0,1.0]
		self.ambient = [1.0, 1.0, 1.0, 1.0]
		self.shininess = 25.0
		self.pdb_file = pdb_file
		self.renderBoundingBox = False
	
	def setSelectedItem(self, is_selected):
		""" Set SG apix to curent selection"""
		EMItem3D.setSelectedItem(self, is_selected)
		sg = self.getRootNode()

	def getBoundingBoxDimensions(self):
		data = self.getData()
		return np.abs(np.max(data,axis=0).astype(np.int)) + np.std(data,axis=0).astype(np.int)

	def getRenderBoundingBox(self): return self.renderBoundingBox
	def setRenderBoundingBox(self,state): self.renderBoundingBox = state
	def getEvalString(self): return "EMPDBItem3D(\"%s\")"%os.path.abspath(self.path)
	def getName(self): return self.model_name
	def getData(self): return self.data

	def setData(self, path):
		if path == None:
			self.path = str(self.attribdict['data_path'].text())
			self.model_name = str(self.attribdict['node_name'].text())
		else:
			try:
				self.path = str(path.text())
			except:
				self.path = str(path)
			self.model_name = self.path.split('/')[-1].split('.')[0]
			self.attribdict = {}
			self.attribdict['data_path'] = self.path
			self.attribdict['node_name'] = self.model_name
		
		self.parser = PDBReader()
		try:
			self.parser.read_from_pdb(self.path)
			self.natoms = self.parser.get_number_points()
			self.data = np.array(self.parser.get_points()).reshape(self.natoms,3)
			self.atom_names = self.parser.get_atomName()
		except:
			print('Could not validate this PDB file. Try another or redownload, as your copy may be corrupt.')
		
		for child in self.getChildren():
			try: child.dataChanged()
			except: pass
		
		bbsize = self.getBoundingBoxDimensions()
		self.boundingboxsize = "%dx%dx%d"%(bbsize[0],bbsize[1],bbsize[2])
		
		if self.item_inspector: self.item_inspector.updateMetaData()

	def getItemDictionary(self):
		"""Return a dictionary of item parameters (used for restoring sessions"""
		return super(EMPDBItem3D, self).getItemDictionary()
	
	def renderNode(self):
		if self.renderBoundingBox: drawBoundingBox(*self.getBoundingBoxDimensions())
	
	def getItemInspector(self):
		"""Return a Qt widget that controls the scene item"""
		if not self.item_inspector: self.item_inspector = EMPDBItem3DInspector(self.model_name, self)
		return self.item_inspector

class EMPDBItem3DInspector(EMItem3DInspector):
	
	def __init__(self, name, item3d):
		EMItem3DInspector.__init__(self, name, item3d)
	
	def updateItemControls(self):
		""" Updates this item inspector. Function is called by the item it observes"""
		super(EMPDBItem3DInspector, self).updateItemControls()
		# Anything that needs to be updated when the scene is rendered goes here.....
		if self.item3d().path:
			self.file_path_label.setText(self.item3d().path)
		self.data_checkbox.setChecked(self.item3d().getRenderBoundingBox())

	def addTabs(self):
		""" Add a tab for each 'column' """
		tabwidget = QtGui.QWidget()
		gridbox = QtGui.QGridLayout()
		tabwidget.setLayout(gridbox)
		self.addTab(tabwidget, "data")
		# add data tab first, then basic
		super(EMPDBItem3DInspector, self).addTabs()
		EMPDBItem3DInspector.addControls(self, gridbox)

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
		if type(self) == EMPDBItem3DInspector: self.updateItemControls()

	def onFileBrowse(self):
		#TODO: replace this with an EMAN2 browser window once we re-write it
		file_path = QtGui.QFileDialog.getOpenFileName(self, "Open PDB Model")
		if file_path:
			self.file_path_label.setText(file_path)
			self.item3d().setData(file_path)
			self.inspector().updateSceneGraph()
	
	def onBBoxChange(self, state):
		self.item3d().setRenderBoundingBox(not self.item3d().getRenderBoundingBox())
		self.inspector().updateSceneGraph()
	
	###########
	# John Flanagan added these methods so the inspector can set the color
	###########
	
	def setAmbientColor(self, red, green, blue, alpha=1.0):
		self.ambient = [red, green, blue, alpha]
	
	def setDiffuseColor(self, red, green, blue, alpha=1.0):
		self.diffuse = [red, green, blue, alpha]

	def setSpecularColor(self, red, green, blue, alpha=1.0):
		self.specular = [red, green, blue, alpha]

	def setShininess(self, shininess):
		self.shininess = shininess
	
	############

class EMBallStickModel(EMPDBItem3D):
	
	"""Ball and stick representation of a PDB model."""
	
	nodetype = "PDBChild"
	representation = "Ball and Stick"

	@staticmethod
	def getNodeDialogWidget(attribdict):
		"""Get Ball and Stick Model Widget"""
		ballstickwidget = QtGui.QWidget()
		grid = QtGui.QGridLayout()
		node_name_model_label = QtGui.QLabel("PDB Structure Name")
		attribdict["node_name"] = QtGui.QLineEdit(str(EMBallStickModel.representation))
		grid.addWidget(node_name_model_label, 0, 0, 1, 2)
		grid.addWidget(attribdict["node_name"], 0, 2, 1, 2)
		EMItem3D.get_transformlayout(grid, 2, attribdict)
		ballstickwidget.setLayout(grid)
		return ballstickwidget

	@staticmethod
	def getNodeForDialog(attribdict):
		"""Create a new node using a attribdict"""
		parent = attribdict["parent"]
		pdb_file = parent.attribdict['data_path']
		transform = EMItem3D.getTransformFromDict(attribdict)
		return EMBallStickModel(pdb_file=pdb_file, parent=parent, transform=transform)
	
	def __init__(self, pdb_file, parent=None, children = set(), transform=None):
		"""
		@param parent: should be an EMPDBItem3D instance for proper functionality.
		"""
		if not transform: transform = Transform()	# Object initialization should not be put in the constructor. Causes issues
		EMPDBItem3D.__init__(self, pdb_file=pdb_file, parent=parent, children=children, transform=transform)
		self.first_render_flag = True # this is used to catch the first call to the render function - so you can do an GL context sensitive initialization when you know there is a valid context
		self.gq = None # will be a glu quadric
		self.dl = None
		self.cylinderdl = 0 # will be a cylinder with no caps
		self.diskdl = 0 # this will be a flat disk
		self.spheredl = 0 # this will be a low resolution sphere
		self.highresspheredl = 0 # high resolution sphere
		self.cappedcylinderdl = 0 # a capped cylinder
		self.radius = 100
		self.colors = get_default_gl_colors()
		amino_acids_list = ["ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TYR","TRP","VAL"]
		self.side_chains = {aa:[] for aa in amino_acids_list}
	
	def load_gl_color(self,name):
		color = self.colors[name]
		glColor(color["ambient"])
		glMaterial(GL_FRONT,GL_AMBIENT,color["ambient"])
		glMaterial(GL_FRONT,GL_DIFFUSE,color["diffuse"])
		glMaterial(GL_FRONT,GL_SPECULAR,color["specular"])
		glMaterial(GL_FRONT,GL_EMISSION,color["emission"])
		glMaterial(GL_FRONT,GL_SHININESS,color["shininess"])
	
	def setLabel(self, label): self.label = label
	def current_text(self):  return self.text
	def getEvalString(self): return "EMBallStickModel()"

	def getItemInspector(self):
		if not self.item_inspector: self.item_inspector = EMBallStickModelInspector("BALL AND STICK", self)
		return self.item_inspector

	def getItemDictionary(self):
		"""Return a dictionary of item parameters (used for restoring sessions"""
		dictionary = super(EMBallStickModel, self).getItemDictionary()
		dictionary.update({"COLOR":[self.ambient, self.diffuse, self.specular, self.shininess]})
		return dictionary

	def setUsingDictionary(self, dictionary):
		"""Set item attributes using a dictionary, used in session restoration"""
		super(EMBallStickItem3D, self).setUsingDictionary(dictionary)
		self.setAmbientColor(dictionary["COLOR"][0][0], dictionary["COLOR"][0][1], dictionary["COLOR"][0][2], dictionary["COLOR"][0][3])
		self.setDiffuseColor(dictionary["COLOR"][1][0], dictionary["COLOR"][1][1], dictionary["COLOR"][1][2], dictionary["COLOR"][1][3])
		self.setSpecularColor(dictionary["COLOR"][2][0], dictionary["COLOR"][2][1], dictionary["COLOR"][2][2], dictionary["COLOR"][2][3])

	def buildResList(self): # calls PDBReader to read the given pdb file and create a list (self.allResidues) of lists (x,y,z,atom name, residue name) of lists (all the values for that residue)
		self.allResidues = []
# 		data = self.getData()
# 		data = data - np.mean(data,axis=0)
		point_x = self.parser.get_x()
		point_y = self.parser.get_y()
		point_z = self.parser.get_z()
		point_x = point_x - np.mean(point_x)
		point_y = point_y - np.mean(point_y)
		point_z = point_z - np.mean(point_z)
		point_atomName = self.parser.get_atomName()
		point_resName = self.parser.get_resName()
		point_resNum = self.parser.get_resNum()
		x =[]
		y =[]
		z =[]
		atomName =[]
		resName = []
		amino = []
		currentRes = point_resNum[0]
		for i in range(len(self.data)):
			if point_resNum[i]==currentRes:
				x.append(point_x[i])
				y.append(point_y[i])
				z.append(point_z[i])
				temp = point_atomName[i]
				temp2 = temp.strip()
				atomName.append(temp2)
				resName.append(point_resName[i])
			else:
				currentRes = point_resNum[i]
				amino.append(x[:])
				amino.append(y[:])
				amino.append(z[:])
				amino.append(atomName[:])
				amino.append(resName[:])
				self.allResidues.append(amino[:])
				del amino[:]
				del x[:]
				del y[:]
				del z[:]
				del atomName[:]
				del resName[:]
				x.append(point_x[i])
				y.append(point_y[i])
				z.append(point_z[i])
				temp = point_atomName[i]
				temp2 = temp.strip()
				atomName.append(temp2)
				resName.append(point_resName[i])
				if (i == (len(point_x)-1)): 
					amino.append(x[:])
					amino.append(y[:])
					amino.append(z[:])
					amino.append(atomName[:])
					amino.append(resName[:])
					self.allResidues.append(amino[:])
					break
	
	def draw_objects(self):
		self.init_basic_shapes() # only does something the first time you call it
		if self.dl == None: #self.dl is the display list, every time a new file is added, this is changed back to None
			self.dl=glGenLists(1)
			glNewList(self.dl,GL_COMPILE)
			self.buildResList()
			for res in self.allResidues: #goes through self.allResidues and displays a sphere for every atom in the pdb
				for i in range(len(res[0])):
					glPushMatrix()
					glTranslate(res[0][i], res[1][i], res[2][i])
					glScale(1,1,1)
					if (str(res[3][i])[0] == 'C'): self.load_gl_color("white")
					elif (str(res[3][i])[0] == 'N'): self.load_gl_color("green")
					elif (str(res[3][i])[0] == 'O'): self.load_gl_color("blue")
					elif (str(res[3][i])[0] == 'S'): self.load_gl_color("red")
					else: self.load_gl_color("silver")
					glCallList(self.highresspheredl)
					glPopMatrix()
			for k in range(len(self.allResidues)):
				res = self.allResidues[k]
				key =  res[4][0]
				if self.side_chains.has_key(key):
					self.renderResidues(res,self)
					continue
				if k !=0: #connects residues together from the nitrogen of one residue to the O of the next residue
					nt = [0,0,0]
					pt = [0,0,0]
					nt[0] = res[0][0]
					nt[1] = res[1][0]
					nt[2] = res[2][0]
					pt[0] = self.allResidues[k-1][0][2]
					pt[1] = self.allResidues[k-1][1][2]
					pt[2] = self.allResidues[k-1][2][2]
					self.cylinder_to_from(nt, pt, 0.2)
			glEndList()
		try:
			glCallList(self.dl)
		except:
			print "call list failed",self.dl
			glDeleteLists(self.dl,1)
			self.dl = None
	
	def init_basic_shapes(self):
		if self.gq == None:
			self.gq=gluNewQuadric() # a quadric for general use
			gluQuadricDrawStyle(self.gq,GLU_FILL)
			gluQuadricNormals(self.gq,GLU_SMOOTH)
			gluQuadricOrientation(self.gq,GLU_OUTSIDE)
			gluQuadricTexture(self.gq,GL_FALSE)
		if self.cylinderdl == 0:
			self.cylinderdl=glGenLists(1)
			glNewList(self.cylinderdl,GL_COMPILE)
			glPushMatrix()
			gluCylinder(self.gq,1.0,1.0,1.0,12,2)
			glPopMatrix()
			glEndList()
		if self.diskdl == 0:
			self.diskdl=glGenLists(1)
			glNewList(self.diskdl,GL_COMPILE)
			gluDisk(self.gq,0,1,12,2)
			glEndList()
		if self.spheredl == 0:
			self.spheredl=glGenLists(1)
			glNewList(self.spheredl,GL_COMPILE)
			gluSphere(self.gq,.5,4,2)
			glEndList()
		if self.highresspheredl == 0:
			self.highresspheredl=glGenLists(1)
			glNewList(self.highresspheredl,GL_COMPILE)
			gluSphere(self.gq,.5,16,16)
			glEndList()
		if self.cappedcylinderdl == 0:
			self.cappedcylinderdl=glGenLists(1)
			glNewList(self.cappedcylinderdl,GL_COMPILE)
			glCallList(self.cylinderdl)
			glPushMatrix()
			glTranslate(0,0,1)
			glCallList(self.diskdl)
			glPopMatrix()
			glPushMatrix()
			glRotate(180,0,1,0)
			glCallList(self.diskdl)
			glPopMatrix()
			glEndList()
	
	def makeStick(self, res, index1, index2): #draws a cylinder between two atoms once the index for start and stop is given
		n = [0,0,0]
		p = [0,0,0]
		p[0] = res[0][index1]
		p[1] = res[1][index1]
		p[2] = res[2][index1]
		n[0] = res[0][index2]
		n[1] = res[1][index2]
		n[2] = res[2][index2]
		self.cylinder_to_from(n, p, 0.2)
	
	def cylinder_to_from(self,next,prev,scale=0.5):
		dx = next[0] - prev[0]
		dy = next[1] - prev[1]
		dz = next[2] - prev[2]
		try: length = np.sqrt(dx**2 + dy**2 + dz**2)
		except: return
		if length == 0: return
		alt = np.arccos(dz/length)*180.0/np.pi
		phi = np.arctan2(dy,dx)*180.0/np.pi
		glPushMatrix()
		glTranslatef(prev[0], prev[1], prev[2] )
		glRotatef(90.0+phi,0,0,1)
		glRotatef(alt,1,0,0)
		glScalef(scale,scale,length)
		self.load_gl_color("silver")
		glCallList(self.cylinderdl)
		glPopMatrix()
	
	def renderNode(self):
		if self.is_selected and glGetIntegerv(GL_RENDER_MODE) == GL_RENDER and not self.isSelectionHidded(): # No need to draw outline in selection mode
			#if glGetIntegerv(GL_RENDER_MODE) == GL_RENDER: print "X"
			glPushAttrib( GL_ALL_ATTRIB_BITS )
			# First render the cylinder, writing the outline to the stencil buffer
			glClearStencil(0)
			glClear( GL_STENCIL_BUFFER_BIT )
			glEnable( GL_STENCIL_TEST )
			# Write to stencil buffer
			glStencilFunc( GL_ALWAYS, 1, 0xFFFF )
			# Only pixels that pass the depth test are written to the stencil buffer
			glStencilOp( GL_KEEP, GL_KEEP, GL_REPLACE )
			glPolygonMode( GL_FRONT_AND_BACK, GL_FILL )
			self.renderShape()
			# Then render the outline
			glStencilFunc( GL_NOTEQUAL, 1, 0xFFFF )
			glStencilOp( GL_KEEP, GL_KEEP, GL_REPLACE )
			# By increasing the line width only the outline is drawn
			glLineWidth( 4.0 )
			glPolygonMode( GL_FRONT_AND_BACK, GL_LINE )
			glMaterialfv(GL_FRONT, GL_EMISSION, [0.0, 1.0, 0.0, 1.0])
			self.renderShape()
			glPopAttrib()	
		else:
			glPushAttrib( GL_ALL_ATTRIB_BITS )
			self.renderShape()
			glPopAttrib()
	
	def renderShape(self):
		glDisable(GL_COLOR_MATERIAL)
		glMaterialfv(GL_FRONT, GL_DIFFUSE, self.diffuse)
		glMaterialfv(GL_FRONT, GL_SPECULAR, self.specular)
		glMaterialf(GL_FRONT, GL_SHININESS, self.shininess)
		glMaterialfv(GL_FRONT, GL_AMBIENT, self.ambient)
		if self.first_render_flag: self.first_render_flag = False
		glPushMatrix()
		self.draw_objects()
		glPopMatrix()
	
	@staticmethod
	def renderResidues(res,target):
		aa = res[4][0]
		if aa == "ALA":
			try: t1 = res[3].index('CB')
			except: pass
			try: target.makeStick(res, 0, 1)
			except: pass	
			try: target.makeStick(res, 1, 2)
			except: pass
			try: target.makeStick(res, 2, 3)
			except: pass
			try: target.makeStick(res, 1, t1)
			except: pass
		elif aa == "ARG":
			try: t1 = res[3].index('CB')
			except: pass
			try: t2 = res[3].index('CG')
			except: pass
			try: t3 = res[3].index('CD')
			except: pass
			try: t4 = res[3].index('NE')
			except: pass
			try: t5 = res[3].index('CZ')
			except: pass
			try: t6 = res[3].index('NH1')
			except: pass
			try: t7 = res[3].index('NH2')
			except: pass
			try: target.makeStick(res, 0, 1)
			except: pass
			try: target.makeStick(res, 1, 2)
			except: pass
			try: target.makeStick(res, 2, 3)
			except: pass
			try: target.makeStick(res, 1, t1)
			except: pass
			try: target.makeStick(res, t1, t2)
			except: pass
			try: target.makeStick(res, t2, t3)
			except: pass
			try: target.makeStick(res, t3, t4)
			except: pass
			try: target.makeStick(res, t4, t5)
			except: pass
			try: target.makeStick(res, t5, t6)
			except: pass
			try: target.makeStick(res, t5, t7)
			except: pass
		elif aa == "ASP":
			try: t1 = res[3].index('CB')
			except: pass
			try: t2 = res[3].index('CG')
			except: pass
			try: t3 = res[3].index('OD1')
			except: pass
			try: t4 = res[3].index('OD2')
			except: pass
			try: target.makeStick(res, 0, 1)
			except: pass
			try: target.makeStick(res, 1, 2)
			except: pass
			try: target.makeStick(res, 2, 3)
			except: pass
			try: target.makeStick(res, 1, t1)
			except: pass
			try: target.makeStick(res, t1, t2)
			except: pass
			try: target.makeStick(res, t2, t3)
			except: pass
			try: target.makeStick(res, t2, t4)
			except: pass
		elif aa == "ASN":
			try: t1 = res[3].index('CB')
			except: pass
			try: t2 = res[3].index('CG')
			except: pass
			try: t3 = res[3].index('OD1')
			except: pass
			try: t4 = res[3].index('ND2')
			except: pass
			try: target.makeStick(res, 0, 1)
			except: pass
			try: target.makeStick(res, 1, 2)
			except: pass
			try: target.makeStick(res, 2, 3)
			except: pass
			try: target.makeStick(res, 1, t1)
			except: pass
			try: target.makeStick(res, t1, t2)
			except: pass
			try: target.makeStick(res, t2, t3)
			except: pass
			try: target.makeStick(res, t2, t4)
			except: pass
		elif aa == "CYS":
			try: t1 = res[3].index('CB')
			except: pass
			try: t2 = res[3].index('SG')
			except: pass
			try: target.makeStick(res, 0, 1)
			except: pass
			try: target.makeStick(res, 1, 2)
			except: pass
			try: target.makeStick(res, 2, 3)
			except: pass
			try: target.makeStick(res, 1, t1)
			except: pass
			try: target.makeStick(res, t1, t2)
			except: pass
		elif aa == "GLY":
			try: target.makeStick(res, 0, 1)
			except: pass
			try: target.makeStick(res, 1, 2)
			except: pass
			try: target.makeStick(res, 2, 3)
			except: pass
		elif aa == "GLN":
			try: t1 = res[3].index('CB')
			except: pass
			try: t2 = res[3].index('CG')
			except: pass
			try: t3 = res[3].index('CD')
			except: pass
			try: t4 = res[3].index('OE1')
			except: pass
			try: t5 = res[3].index('NE2')
			except: pass
			try: target.makeStick(res, 0, 1)
			except: pass
			try: target.makeStick(res, 1, 2)
			except: pass
			try: target.makeStick(res, 2, 3)
			except: pass
			try: target.makeStick(res, 1, t1)
			except: pass
			try: target.makeStick(res, t1, t2)
			except: pass
			try: target.makeStick(res, t2, t3)
			except: pass
			try: target.makeStick(res, t3, t4)
			except: pass
			try: target.makeStick(res, t3, t5)
			except: pass
		elif aa == "GLU":
			try: t1 = res[3].index('CB')
			except: pass
			try: t2 = res[3].index('CG')
			except: pass
			try: t3 = res[3].index('CD')
			except: pass
			try: t4 = res[3].index('OE1')
			except: pass
			try: t5 = res[3].index('OE2')
			except: pass
			try: target.makeStick(res, 0, 1)
			except: pass
			try: target.makeStick(res, 1, 2)
			except: pass
			try: target.makeStick(res, 2, 3)
			except: pass
			try: target.makeStick(res, 1, t1)
			except: pass
			try: target.makeStick(res, t1, t2)
			except: pass
			try: target.makeStick(res, t2, t3)
			except: pass
			try: target.makeStick(res, t3, t4)
			except: pass
			try: target.makeStick(res, t3, t5)
			except: pass
		elif aa == "HIS":
			try: t1 = res[3].index('CB')
			except: pass
			try: t2 = res[3].index('CG')
			except: pass
			try: t3 = res[3].index('CD2')
			except: pass
			try: t4 = res[3].index('ND1')
			except: pass
			try: t5 = res[3].index('NE2')
			except: pass
			try: t6 = res[3].index('CE1')
			except: pass
			try: target.makeStick(res, 0, 1)
			except: pass
			try: target.makeStick(res, 1, 2)
			except: pass
			try: target.makeStick(res, 2, 3)
			except: pass
			try: target.makeStick(res, 1, t1)
			except: pass
			try: target.makeStick(res, t1, t2)
			except: pass
			try: target.makeStick(res, t2, t3)
			except: pass
			try: target.makeStick(res, t2, t4)
			except: pass
			try: target.makeStick(res, t3, t5)
			except: pass
			try: target.makeStick(res, t5, t6)
			except: pass
			try: target.makeStick(res, t4, t6)
			except: pass
		elif aa == "ILE":
			try: t1 = res[3].index('CB')
			except: pass
			try: t2 = res[3].index('CG1')
			except: pass
			try: t3 = res[3].index('CG2')
			except: pass
			try: t4 = res[3].index('CD1')
			except: pass
			try: target.makeStick(res, 0, 1)
			except: pass
			try: target.makeStick(res, 1, 2)
			except: pass
			try: target.makeStick(res, 2, 3)
			except: pass
			try: target.makeStick(res, 1, t1)
			except: pass
			try: target.makeStick(res, t1, t2)
			except: pass
			try: target.makeStick(res, t1, t3)
			except: pass
			try: target.makeStick(res, t2, t4)
			except: pass
		elif aa == "LEU":
			try: t1 = res[3].index('CB')
			except: pass
			try: t2 = res[3].index('CG')
			except: pass
			try: t3 = res[3].index('CD1')
			except: pass
			try: t4 = res[3].index('CD2')
			except: pass
			try: target.makeStick(res, 0, 1)
			except: pass
			try: target.makeStick(res, 1, 2)
			except: pass
			try: target.makeStick(res, 2, 3)
			except: pass
			try: target.makeStick(res, 1, t1)
			except: pass
			try: target.makeStick(res, t1, t2)
			except: pass
			try: target.makeStick(res, t2, t3)
			except: pass
			try: target.makeStick(res, t2, t4)
			except: pass
		elif aa == "LYS":
			try: t1 = res[3].index('CB')
			except: pass
			try: t2 = res[3].index('CG')
			except: pass
			try: t3 = res[3].index('CD')
			except: pass
			try: t4 = res[3].index('CE')
			except: pass
			try: t5 = res[3].index('NZ')
			except: pass
			try: target.makeStick(res, 0, 1)
			except: pass
			try: target.makeStick(res, 1, 2)
			except: pass
			try: target.makeStick(res, 2, 3)
			except: pass
			try: target.makeStick(res, 1, t1)
			except: pass
			try: target.makeStick(res, t1, t2)
			except: pass
			try: target.makeStick(res, t2, t3)
			except: pass
			try: target.makeStick(res, t3, t4)
			except: pass
			try: target.makeStick(res, t4, t5)
			except: pass
		elif aa == "MET":
			try: t1 = res[3].index('CB')
			except: pass
			try: t2 = res[3].index('CG')
			except: pass
			try: t3 = res[3].index('SD')
			except: pass
			try: t4 = res[3].index('CE')
			except: pass
			try: target.makeStick(res, 0, 1)
			except: pass
			try: target.makeStick(res, 1, 2)
			except: pass
			try: target.makeStick(res, 2, 3)
			except: pass
			try: target.makeStick(res, 1, t1)
			except: pass
			try: target.makeStick(res, t1, t2)
			except: pass
			try: target.makeStick(res, t2, t3)
			except: pass
			try: target.makeStick(res, t3, t4)
			except: pass
		elif aa == "PHE":
			try: t1 = res[3].index('CB')
			except: pass
			try: t2 = res[3].index('CG')
			except: pass
			try: t3 = res[3].index('CD1')
			except: pass
			try: t4 = res[3].index('CD2')
			except: pass
			try: t5 = res[3].index('CE1')
			except: pass
			try: t6 = res[3].index('CE2')
			except: pass
			try: t7 = res[3].index('CZ')
			except: pass
			try: target.makeStick(res, 0, 1)
			except: pass
			try: target.makeStick(res, 1, 2)
			except: pass
			try: target.makeStick(res, 2, 3)
			except: pass
			try: target.makeStick(res, 1, t1)
			except: pass
			try: target.makeStick(res, t1, t2)
			except: pass
			try: target.makeStick(res, t2, t3)
			except: pass
			try: target.makeStick(res, t2, t4)
			except: pass
			try: target.makeStick(res, t3, t5)
			except: pass
			try: target.makeStick(res, t4, t6)
			except: pass
			try: target.makeStick(res, t5, t7)
			except: pass
			try: target.makeStick(res, t6, t7)
			except: pass
		elif aa == "PRO":
			try: t1 = res[3].index('CB')
			except: pass
			try: t2 = res[3].index('CG')
			except: pass
			try: t3 = res[3].index('CD')
			except: pass
			try: t4 = res[3].index('N')
			except: pass
			try: target.makeStick(res, 0, 1)
			except: pass
			try: target.makeStick(res, 1, 2)
			except: pass
			try: target.makeStick(res, 2, 3)
			except: pass
			try: target.makeStick(res, 1, t1)
			except: pass
			try: target.makeStick(res, t1, t2)
			except: pass
			try: target.makeStick(res, t2, t3)
			except: pass
			try: target.makeStick(res, t3, t4)
			except: pass
		elif aa == "SER":
			try: t1 = res[3].index('CB')
			except: pass
			try: t2 = res[3].index('OG')
			except: pass
			try: target.makeStick(res, 0, 1)
			except: pass
			try: target.makeStick(res, 1, 2)
			except: pass
			try: target.makeStick(res, 2, 3)
			except: pass
			try: target.makeStick(res, 1, t1)
			except: pass
			try: target.makeStick(res, t1, t2)
			except: pass
		elif aa == "THR":
			try: t1 = res[3].index('CB')
			except: pass
			try: t2 = res[3].index('CG2')
			except: pass
			try: t3 = res[3].index('OG1')
			except: pass
			try: target.makeStick(res, 0, 1)
			except: pass
			try: target.makeStick(res, 1, 2)
			except: pass
			try: target.makeStick(res, 2, 3)
			except: pass
			try: target.makeStick(res, 1, t1)
			except: pass
			try: target.makeStick(res, t1, t2)
			except: pass
			try: target.makeStick(res, t1, t3)
			except: pass
		elif aa == "TRP":
			try: t1 = res[3].index('CB')
			except: pass
			try: t2 = res[3].index('CG')
			except: pass
			try: t3 = res[3].index('CD1')
			except: pass
			try: t4 = res[3].index('CD2')
			except: pass
			try: t5 = res[3].index('NE1')
			except: pass
			try: t6 = res[3].index('CE2')
			except: pass
			try: t7 = res[3].index('CE3')
			except: pass
			try: t8 = res[3].index('CZ3')
			except: pass
			try: t9 = res[3].index('CH2')
			except: pass
			try: t10 = res[3].index('CZ2')
			except: pass
			try: target.makeStick(res, 0, 1)
			except: pass
			try: target.makeStick(res, 1, 2)
			except: pass
			try: target.makeStick(res, 2, 3)
			except: pass
			try: target.makeStick(res, 1, t1)
			except: pass
			try: target.makeStick(res, t1, t2)
			except: pass
			try: target.makeStick(res, t2, t3)
			except: pass
			try: target.makeStick(res, t2, t4)
			except: pass
			try: target.makeStick(res, t3, t5)
			except: pass
			try: target.makeStick(res, t5, t6)
			except: pass
			try: target.makeStick(res, t4, t6)
			except: pass
			try: target.makeStick(res, t4, t7)
			except: pass
			try: target.makeStick(res, t7, t8)
			except: pass
			try: target.makeStick(res, t8, t9)
			except: pass
			try: target.makeStick(res, t10, t9)
			except: pass
		elif aa == "VAL":
			try: t1 = res[3].index('CB')
			except: pass
			try: t2 = res[3].index('CG2')
			except: pass
			try: t3 = res[3].index('CG1')
			except: pass
			try: target.makeStick(res, 0, 1)
			except: pass
			try: target.makeStick(res, 1, 2)
			except: pass
			try: target.makeStick(res, 2, 3)
			except: pass
			try: target.makeStick(res, 1, t1)
			except: pass
			try: target.makeStick(res, t1, t2)
			except: pass
			try: target.makeStick(res, t1, t3)
			except: pass
		elif aa == "TYR":
			try: t1 = res[3].index('CB')
			except: pass
			try: t2 = res[3].index('CG')
			except: pass
			try: t3 = res[3].index('CD1')
			except: pass
			try: t4 = res[3].index('CD2')
			except: pass
			try: t5 = res[3].index('CE1')
			except: pass
			try: t6 = res[3].index('CE2')
			except: pass
			try: t7 = res[3].index('CZ')
			except: pass
			try: t8 = res[3].index('OH')
			except: pass
			try: target.makeStick(res, 0, 1)
			except: pass
			try: target.makeStick(res, 1, 2)
			except: pass
			try: target.makeStick(res, 2, 3)
			except: pass
			try: target.makeStick(res, 1, t1)
			except: pass
			try: target.makeStick(res, t1, t2)
			except: pass
			try: target.makeStick(res, t2, t3)
			except: pass
			try: target.makeStick(res, t2, t4)
			except: pass
			try: target.makeStick(res, t3, t5)
			except: pass
			try: target.makeStick(res, t4, t6)
			except: pass
			try: target.makeStick(res, t5, t7)
			except: pass
			try: target.makeStick(res, t6, t7)
			except: pass
			try: target.makeStick(res, t7, t8)
			except: pass

class EMBallStickModelInspector(EMPDBItem3DInspector):
	
	def __init__(self, name, item3d):
		EMPDBItem3DInspector.__init__(self, name, item3d)

class EMSphereModel(EMPDBItem3D):
	
	"""Spheres representation of a PDB model."""

	nodetype = "PDBChild"
	representation = "Spheres"
		
	@staticmethod
	def getNodeDialogWidget(attribdict):
		"""Get Spheres Model Widget"""
		sphereswidget = QtGui.QWidget()
		grid = QtGui.QGridLayout()
		node_name_model_label = QtGui.QLabel("PDB Structure Name")
		attribdict["node_name"] = QtGui.QLineEdit(str(EMSphereModel.representation))
		grid.addWidget(node_name_model_label, 0, 0, 1, 2)
		grid.addWidget(attribdict["node_name"], 0, 2, 1, 2)
		EMItem3D.get_transformlayout(grid, 2, attribdict)
		sphereswidget.setLayout(grid)
		return sphereswidget

	@staticmethod
	def getNodeForDialog(attribdict):
		"""Create a new node using a attribdict"""
		parent = attribdict["parent"]
		pdb_file = parent.attribdict["data_path"]
		transform=EMItem3D.getTransformFromDict(attribdict)
		return EMSphereModel(pdb_file=pdb_file, parent=parent, transform=transform)
	
	@staticmethod
	def get_vanderwaals_radius(name):
		# returns van der waals radius in angstroms
		if 'H' in name: return 1.2
		elif 'C' in name: return 1.7
		elif 'N' in name: return 1.5
		elif 'O' in name: return 1.4
		elif 'F' in name: return 1.35
		elif 'P' in name: return 1.9
		elif 'S' in name: return 1.85
		elif 'CL' in name: return 1.8
		else: return 1.5
	
	def __init__(self, pdb_file, parent=None, children = set(), transform=None, radius=1.75):
		"""
		@param parent: should be an EMPDBItem3D instance for proper functionality.
		"""
		if not transform: transform = Transform()	# Object initialization should not be put in the constructor. Causes issues
		EMPDBItem3D.__init__(self, pdb_file=pdb_file, parent=parent, children=children, transform=transform)
		self.first_render_flag = True # this is used to catch the first call to the render function - so you can do an GL context sensitive initialization when you know there is a valid context
		self.gq = None # will be a glu quadric
		self.dl = None
# 		self.spheredl = 0 # this will be a low resolution sphere
		self.highresspheredl = 0 # high resolution sphere
		self.radius = radius
		self.colors = get_default_gl_colors()
		self.vwr = [self.get_vanderwaals_radius(an) for an in self.atom_names]
		self.coords = self.data - np.mean(self.data,axis=0)
	
	def current_text(self):
		return self.text
	
	def load_gl_color(self,name):
		color = self.colors[name]
		glColor(color["ambient"])
		glMaterial(GL_FRONT,GL_AMBIENT,color["ambient"])
		glMaterial(GL_FRONT,GL_DIFFUSE,color["diffuse"])
		glMaterial(GL_FRONT,GL_SPECULAR,color["specular"])
		glMaterial(GL_FRONT,GL_EMISSION,color["emission"])
		glMaterial(GL_FRONT,GL_SHININESS,color["shininess"])

	def getEvalString(self):
		return "EMSphereModel()"

	def getItemInspector(self):
		if not self.item_inspector: self.item_inspector = EMSphereModelInspector("SPHERES", self)
		return self.item_inspector

	def getItemDictionary(self):
		"""Return a dictionary of item parameters (used for restoring sessions"""
		dictionary = super(EMSphereModel, self).getItemDictionary()
		dictionary.update({"COLOR":[self.ambient, self.diffuse, self.specular, self.shininess]})
		return dictionary

	def setUsingDictionary(self, dictionary):
		"""Set item attributes using a dictionary, used in session restoration"""
		super(EMSphereModel, self).setUsingDictionary(dictionary)
		self.setAmbientColor(dictionary["COLOR"][0][0], dictionary["COLOR"][0][1], dictionary["COLOR"][0][2], dictionary["COLOR"][0][3])
		self.setDiffuseColor(dictionary["COLOR"][1][0], dictionary["COLOR"][1][1], dictionary["COLOR"][1][2], dictionary["COLOR"][1][3])
		self.setSpecularColor(dictionary["COLOR"][2][0], dictionary["COLOR"][2][1], dictionary["COLOR"][2][2], dictionary["COLOR"][2][3])
	
	def init_basic_shapes(self):
		if self.gq == None:
			self.gq=gluNewQuadric() # a quadric for general use
			gluQuadricDrawStyle(self.gq,GLU_FILL)
			gluQuadricNormals(self.gq,GLU_SMOOTH)
			gluQuadricOrientation(self.gq,GLU_OUTSIDE)
			gluQuadricTexture(self.gq,GL_FALSE)
		if self.highresspheredl == 0:
			self.highresspheredl=glGenLists(1)
			glNewList(self.highresspheredl,GL_COMPILE)
			gluSphere(self.gq,self.radius,20,20)
			glEndList()
	
	def renderNode(self):
		if self.is_selected and glGetIntegerv(GL_RENDER_MODE) == GL_RENDER and not self.isSelectionHidded(): # No need to draw outline in selection mode
			#if glGetIntegerv(GL_RENDER_MODE) == GL_RENDER: print "X"
			glPushAttrib( GL_ALL_ATTRIB_BITS )
			# First render the cylinder, writing the outline to the stencil buffer
			glClearStencil(0)
			glClear( GL_STENCIL_BUFFER_BIT )
			glEnable( GL_STENCIL_TEST )
			# Write to stencil buffer
			glStencilFunc( GL_ALWAYS, 1, 0xFFFF )
			# Only pixels that pass the depth test are written to the stencil buffer
			glStencilOp( GL_KEEP, GL_KEEP, GL_REPLACE )
			glPolygonMode( GL_FRONT_AND_BACK, GL_FILL )
			self.renderShape()
			# Then render the outline
			glStencilFunc( GL_NOTEQUAL, 1, 0xFFFF )
			glStencilOp( GL_KEEP, GL_KEEP, GL_REPLACE )
			# By increasing the line width only the outline is drawn
			glLineWidth( 4.0 )
			glPolygonMode( GL_FRONT_AND_BACK, GL_LINE )
			glMaterialfv(GL_FRONT, GL_EMISSION, [0.0, 1.0, 0.0, 1.0])
			self.renderShape()
			glPopAttrib()
		else:
			glPushAttrib( GL_ALL_ATTRIB_BITS )
			self.renderShape()
			glPopAttrib()
	
	def renderShape(self):
		glDisable(GL_COLOR_MATERIAL)
		glMaterialfv(GL_FRONT, GL_DIFFUSE, self.diffuse)
		glMaterialfv(GL_FRONT, GL_SPECULAR, self.specular)
		glMaterialf(GL_FRONT, GL_SHININESS, self.shininess)
		glMaterialfv(GL_FRONT, GL_AMBIENT, self.ambient)
		if self.first_render_flag: self.first_render_flag = False
		glPushMatrix()
		self.draw_objects()
		glPopMatrix()
	
	def draw_objects(self):
		self.init_basic_shapes() # only does something the first time you call it
		if self.dl == None: #self.dl is the display list, every time a new file is added, this is changed back to None
			self.dl=glGenLists(1)
			glNewList(self.dl,GL_COMPILE)
			for i in xrange(self.natoms):
				glPushMatrix()
				glTranslate(self.coords[i][0],self.coords[i][1],self.coords[i][2])
				glScale(self.vwr[i],self.vwr[i],self.vwr[i])
				if 'C' in self.atom_names[i]: self.load_gl_color("white")
				elif 'N' in self.atom_names[i]: self.load_gl_color("green")
				elif 'O' in self.atom_names[i]: self.load_gl_color("blue")
				elif 'S' in self.atom_names[i]: self.load_gl_color("red")
				elif 'H' in self.atom_names[i]: self.load_gl_color("yellow")
				else: self.load_gl_color("dark_grey")
				
				glCallList(self.highresspheredl)
				glPopMatrix()
			glEndList()
		try:
			glCallList(self.dl)
		except:
			print "call list failed",self.dl
			glDeleteLists(self.dl,1)
			self.dl = None

class EMSphereModelInspector(EMPDBItem3DInspector):
	
	def __init__(self, name, item3d):
		EMPDBItem3DInspector.__init__(self, name, item3d)

if __name__ == '__main__' :
	print("WARNING: This module is not designed to be run as a program. The browser you see is for testing purposes.")
	from emapplication import EMApp
	from embrowser import EMBrowserWidget
	app = EMApp()
	browser = EMBrowserWidget(withmodal = False, multiselect = False)
	browser.show()
	app.execute()
	