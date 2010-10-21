#!/usr/bin/env python

# Author: Muthu Alagappan, m.alagappan901@gmail.com,  07/22/09
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
#

from EMAN2 import PDBReader, EMData
from PyQt4 import QtCore, QtGui
from emimage3d import EMImage3DWidget
from emimage3diso import EMIsosurfaceModel
from empdbviewer import *

class EMPDBValTool(QtCore.QObject):
	'''
	EMPDB versus isosurface visual evaluation
	'''
	def __init__(self, application=None,ensure_gl_context=True,application_control=True):
		QtCore.QObject.__init__(self)
		self.pdb_model = None # will eventually be a EMPDBModel
		self.iso_model = None # will eventually be a EMIsosurfaceModel

		self.current_mrc = "" #holds the current mrc file name
		self.current_pdb = "" #holds the current pdb file name

		self.t = 0

		self.parent_widget = EMImage3DWidget()
		self.__init_pdb_model()
		self.__init_iso_model()
		
		self.inspector = None
		self.get_inspector()
		self.parent_widget.inspector = self.inspector #FIXME: this seems like a hack

	def __init_pdb_model(self):
		if self.pdb_model == None:
			self.pdb_model = EMPDBModel(None,False,False)
			self.parent_widget.add_model(self.pdb_model)
			self.__set_model_contexts(self.pdb_model)

			
	def __set_model_contexts(self,model):
		model.set_qt_context_parent(self.parent_widget)
		model.set_gl_context_parent(self.parent_widget)
		model.set_gl_widget(self.parent_widget)
		model.set_dont_delete_parent() # stops a RunTimeError
		model.under_qt_control = True #self.under_qt_control
		
	def __init_iso_model(self):
		if self.iso_model == None:
			self.iso_model = EMIsosurfaceModel(None,None,False,False, enable_file_browse=True)
			self.parent_widget.add_model(self.iso_model)
			self.__set_model_contexts(self.iso_model)

	
	def set_pdb_file(self,pdb_file):
		self.current_pdb = pdb_file
		print self.current_pdb
		if self.pdb_model == None:
			self.__init_pdb_model()
		self.pdb_model.set_current_text(pdb_file) #updates GL with the new pdb file

		
	def set_iso_file(self,file_name):
		self.current_mrc = file_name
		data = EMData(file_name)
		self.set_iso_data(data)
		self.iso_model.get_inspector().mrc_text.setText(file_name)

	def set_iso_data(self,data):
		if self.iso_model == None: self. __init_iso_model()
		self.iso_model.set_data(data)
		print "set force update"
		self.iso_model.set_force_update(True)
		self.iso_model.updateGL()

	def update_pdb_file(self):  #updates current_pdb with the current file
		self.current_pdb = self.pdb_model.get_pdb_file()

	def update_mrc_file(self):  #updates current_mrc with the current file
		if self.iso_model.get_inspector().mrcChanged:		
			self.current_mrc = self.iso_model.get_mrc_file()
	
	def draw_objects(self):
		if self.pdb_model == None:
			self.__init_pdb_model()
		if self.iso_model == None: 
			self. __init_iso_model()
		if self.pdb_model != None:
			glPushMatrix()
			self.pdb_model.draw_objects()
			glPopMatrix()
		if self.iso_model != None:
			glPushMatrix()
			self.iso_model.render()
			glPopMatrix()
			
	def get_inspector(self):
		if self.inspector == None:
			self.inspector = EMPDBValToolInspector(self)
		return self.inspector

	def run_validate(self): #called by the inspector to tell the mediator to validate the fit
		val1 = int(self.get_inspector().text1.text()) #val 1 = number of transformations
		val2 = float(self.get_inspector().text2.text()) #val 2 = isosurface threshold value

		self.emit(QtCore.SIGNAL("HELP!"))
		self.emit(QtCore.SIGNAL("run_validate"), str(self.current_mrc), str(self.current_pdb), val1, val2)


class EMPDBValToolInspector(QtGui.QWidget):
	def __init__(self,target,enable_advanced=False):
		QtGui.QWidget.__init__(self,None)
		self.target = weakref.ref(target) # prevent a strong cycle - this target object should be an EM3DModel, but that could change depending on who builds on this object
		self.setWindowIcon(QtGui.QIcon(get_image_directory() +"single_image_3d.png"))
		
		self.tabwidget = QtGui.QTabWidget()
		self.vbl = QtGui.QVBoxLayout(self) # this is the main vbl
		self.vbl.addWidget(self.tabwidget) 
		self.setLayout(self.vbl)

		self.advanced_tab = None
		self.text1 = "" #text line edit for "number of transformations"
		self.text2 = "" #text line edit for "isosurface threshold value"
		self.options_module = None # will eventually be the options tab
		
		self.addTab(self.target().iso_model.get_inspector(),"Isosurface")
		if enable_advanced:
			self.insert_advanced_tab()
		self.addTab(self.target().pdb_model.get_inspector(), "PDB")
		self.__init_options_module()

	def __init_options_module(self):
		if self.options_module == None:

			self.opt_tab = QtGui.QWidget()
			opt_tab = self.opt_tab

			v_main = QtGui.QVBoxLayout(self.opt_tab)	
			v_main.setMargin(0)
			v_main.setSpacing(3)
			v_main.setObjectName("Options")

			hbl1 = QtGui.QHBoxLayout()
			hbl1.setMargin(0)
			hbl1.setSpacing(6)
			text1_label = QtGui.QLabel("Number of Transformations: ")
			hbl1.addWidget(text1_label)
			self.text1 = QtGui.QLineEdit()
			self.text1.setAlignment(Qt.AlignRight)
			hbl1.addWidget(self.text1)
			v_main.addLayout(hbl1)

			hbl2 = QtGui.QHBoxLayout()
			hbl2.setMargin(0)
			hbl2.setSpacing(48)
			text2_label = QtGui.QLabel("Isosurface Threshold: ")
			hbl2.addWidget(text2_label)
			self.text2 = QtGui.QLineEdit()
			self.text2.setAlignment(Qt.AlignRight)
			hbl2.addWidget(self.text2)
			v_main.addLayout(hbl2)

			self.validate = QtGui.QPushButton("Validate")
			v_main.addWidget(self.validate)

			self.addTab(self.opt_tab,"Validate Options")

			QtCore.QObject.connect(self.validate, QtCore.SIGNAL("clicked(bool)"), self.runValidate)

	def addTab(self,widget,name):
		self.tabwidget.addTab(widget,name)
		
	def insert_advanced_tab(self):
		#FIXME: doesn't work (EMImageInspector3D vs EMPDBValToolInspector)
		if self.advanced_tab == None:
			from emimage3d import EM3DAdvancedInspector
			self.advanced_tab = EM3DAdvancedInspector(self.target(), self)
			
		self.advanced_tab.update_rotations(self.target().get_current_transform())
		self.advanced_tab.set_scale(self.target().cam.scale)
		self.tabwidget.addTab(self.advanced_tab,"Advanced")
		self.settingsrow = self.tabwidget.count()-1
		self.tabwidget.setCurrentIndex(self.settingsrow)
		
	def runValidate(self, i):
		self.target().update_pdb_file()
		self.target().update_mrc_file()
		
		if (str(self.text1.text()) == ""): self.text1.setText("20") #default value 20
		if (str(self.text2.text()) == ""): self.text2.setText("0.1") #default value 5.0

		if ((str(self.target().current_mrc) != "") and (str(self.target().current_pdb)!= "")): #only calls validate if both strings are valid
			self.target().run_validate()
		if (str(self.target().current_mrc) == ""):
			print "Please properly load an mrc file"
		if (str(self.target().current_pdb) == ""):
			print "Please properly load a pdb file"

	def set_directional_light_dir(self,d):
		if self.advanced_tab: self.advanced_tab.set_directional_light_dir(d)
		
	def set_positional_light_dir(self,d):
		if self.advanced_tab: self.advanced_tab.set_positional_light_dir(d)
			
	def set_positional_light_pos(self,d):
		if self.advanced_tab: self.advanced_tab.set_positional_light_pos(d)
			
	def set_scale(self,val):
		if self.advanced_tab: self.advanced_tab.set_scale(val)
	
	def set_xy_trans(self, x, y):
		if self.advanced_tab: self.advanced_tab.set_xy_trans(x,y)
	
	def set_xyz_trans(self,x,y,z):
		if self.advanced_tab: self.advanced_tab.set_xyz_trans(x,y,z)
		
	def update_rotations(self,t3d):
		if self.advanced_tab:
			self.advanced_tab.update_rotations(t3d)

if __name__ == '__main__':
	from emapplication import EMApp
	em_app = EMApp()
	pdb_tool = EMPDBValTool()

	pdb_tool.set_pdb_file("fh-solution-0-1UF2-T.pdb")
	pdb_tool.set_iso_file("rdv-target2.mrc")
	em_app.show()
	em_app.execute()

		
		
