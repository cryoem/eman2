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

from em3dmodule import *
from EMAN2 import PDBReader

from empdbviewer import *
from emimage3diso import *
import sys
import os


class EMPDBValTool(EM3DModule):
	'''
	EMPDB versus isosurface visual evaluation
	'''
	def __init__(self, application=None,ensure_gl_context=True,application_control=True):
		EM3DModule.__init__(self,application,ensure_gl_context=ensure_gl_context,application_control=application_control)
		self.pdb_module = None # will eventually be a EMPDBViewer
		self.iso_module = None # will eventuall be a EMImage3DModule

		self.current_mrc = "" #holds the current mrc file name
		self.current_pdb = "" #holds the current pdb file name

		self.t = 0

	def __init_pdb_module(self):
		if self.pdb_module == None:
			self.pdb_module = EMPDBViewer(None,False,False)
			self.__set_module_contexts(self.pdb_module)
			self.get_inspector().addTab(self.pdb_module.get_inspector(),"PDB")
			
	def __set_module_contexts(self,module):
		module.set_qt_context_parent(self.qt_context_parent)
		module.set_gl_context_parent(self.gl_context_parent)
		module.set_gl_widget(self.gl_context_parent)
		module.set_dont_delete_parent() # stops a RunTimeError
		module.under_qt_control = self.under_qt_control
		
	def __init_iso_module(self):
		if self.iso_module == None:

			self.iso_module = EMIsosurfaceModule(None,None,False,False, enable_file_browse=True)

			self.get_inspector().addTab(self.iso_module.get_inspector(),"Isosurface")
			self.__set_module_contexts(self.iso_module)

	
	def set_pdb_file(self,pdb_file):
		self.current_pdb = pdb_file
		print self.current_pdb
		if self.pdb_module == None:
			self.__init_pdb_module()
		self.pdb_module.set_current_text(pdb_file) #updates GL with the new pdb file

		
	def set_iso_file(self,file_name):
		self.current_mrc = file_name
		data = EMData(file_name)
		self.set_iso_data(data)
		self.iso_module.get_inspector().mrc_text.setText(file_name)

	def set_iso_data(self,data):
		if self.iso_module == None: self. __init_iso_module()
		self.iso_module.set_data(data)
		print "set force update"
		self.iso_module.set_force_update(True)
		self.iso_module.updateGL()

	def update_pdb_file(self):  #updates current_pdb with the current file
		self.current_pdb = self.pdb_module.get_pdb_file()

	def update_mrc_file(self):  #updates current_mrc with the current file
		if self.iso_module.get_inspector().mrcChanged:		
			self.current_mrc = self.iso_module.get_mrc_file()
	
	def draw_objects(self):
		
		if self.pdb_module == None:
			self.__init_pdb_module()

		if self.iso_module == None: 
			self. __init_iso_module()

		if self.pdb_module != None:
			glPushMatrix()
			self.pdb_module.draw_objects()
			glPopMatrix()
			
		if self.iso_module != None:
			glPushMatrix()
			self.iso_module.render()
			glPopMatrix()
			
	def get_inspector(self):
		if self.inspector == None:
			self.inspector = EMPDBValToolInspector(self)
		return self.inspector

	def run_validate(self): #called by the inspector to tell the mediator to validate the fit
		val1 = int(self.get_inspector().text1.text()) #val 1 = number of transformations
		val2 = float(self.get_inspector().text2.text()) #val 2 = isosurface threshold value

		self.emit(QtCore.SIGNAL("run_validate"), str(self.current_mrc), str(self.current_pdb), val1, val2)

class EMPDBValToolInspector(EM3DInspector):
	def __init__(self,target,enable_advanced=False):
		EM3DInspector.__init__(self,target,enable_advanced)

		self.text1 = "" #text line edit for "number of transformations"
		self.text2 = "" #text line edit for "isosurface threshold value"

		self.options_module = None # will eventually be the options tab
		self.__init_options_module()


		
	def addTab(self,widget,name):
		self.tabwidget.addTab(widget,name)


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

if __name__ == '__main__':
	from emapplication import EMStandAloneApplication
	em_app = EMStandAloneApplication()
	window = EMPDBValTool()
	em_app.show()

	window.set_pdb_file("fh-solution-0-1UF2-T.pdb")
	window.set_iso_file("rdv-target2.mrc")

	em_app.execute()

		
		
