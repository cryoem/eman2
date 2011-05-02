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
from emapplication import EMApp, get_application
from PyQt4 import QtCore, QtGui
from emimage3d import EMImage3DWidget
from emimage3diso import EMIsosurfaceModel
from empdbviewer import *
from emselector import EMSelectorDialog


class EMPDBValWidget(QtGui.QWidget):
	'''
	EMPDB versus isosurface visual evaluation
	'''
	def __init__(self):
		QtGui.QWidget.__init__(self)

		self.pdb_model = None # will eventually be a EMPDBModel
		self.iso_model = None # will eventually be a EMIsosurfaceModel
		self.viewer_window = EMImage3DWidget()
		self.__init_gui()

		self.validate_button.clicked.connect(self.run_validate)
		self.pdb_browse_button.clicked.connect(self.browse_pdb)
		self.pdb_line_edit.textChanged.connect(self.update_pdb_file)
		self.volume_browse_button.clicked.connect(self.browse_iso)
		self.volume_line_edit.textChanged.connect(self.update_iso_file)
		
		get_application().attach_child(self)
		
	def __init_gui(self):
		self.data_groupbox = QtGui.QGroupBox(self.tr("Data"))
		
		pdb_label = QtGui.QLabel("PDB:")
		self.pdb_line_edit = QtGui.QLineEdit()
		self.pdb_browse_button = QtGui.QPushButton(self.tr("Browse"))
		pdb_layout = QtGui.QHBoxLayout()
		pdb_layout.addWidget(pdb_label)
		pdb_layout.addWidget(self.pdb_line_edit)
		pdb_layout.addWidget(self.pdb_browse_button)		
		
		volume_label = QtGui.QLabel("Volume:")
		self.volume_line_edit = QtGui.QLineEdit()
		self.volume_browse_button = QtGui.QPushButton(self.tr("Browse"))
		volume_layout = QtGui.QHBoxLayout()
		volume_layout.addWidget(volume_label)
		volume_layout.addWidget(self.volume_line_edit)
		volume_layout.addWidget(self.volume_browse_button)
		
		data_layout = QtGui.QVBoxLayout()
		data_layout.addLayout(pdb_layout)
		data_layout.addLayout(volume_layout)
		
		self.data_groupbox.setLayout(data_layout)

		
		self.validation_groupbox = QtGui.QGroupBox(self.tr("Validation"))
		
		transformations_label = QtGui.QLabel(self.tr("&Number of Transformations"))
		self.transformations_spinbox = QtGui.QSpinBox()
		self.transformations_spinbox.setMaximum(9999)
		self.transformations_spinbox.setValue(20)
		transformations_label.setBuddy(self.transformations_spinbox)
		transformations_layout = QtGui.QHBoxLayout()
		transformations_layout.addWidget(transformations_label)
		transformations_layout.addWidget(self.transformations_spinbox)
		
		threshold_label = QtGui.QLabel(self.tr("Isosurface &Threshold"))
		self.threshold_doublespinbox = QtGui.QDoubleSpinBox()
		self.threshold_doublespinbox.setValue(0.1)
		threshold_label.setBuddy(self.threshold_doublespinbox)
		threshold_layout = QtGui.QHBoxLayout()
		threshold_layout.addWidget(threshold_label)
		threshold_layout.addWidget(self.threshold_doublespinbox)
		
		self.validate_button = QtGui.QPushButton(self.tr("&Validate"))
		
		validation_layout = QtGui.QVBoxLayout()
		validation_layout.addLayout(transformations_layout)
		validation_layout.addLayout(threshold_layout)
		validation_layout.addWidget(self.validate_button)
		self.validation_groupbox.setLayout(validation_layout)
		
		layout = QtGui.QVBoxLayout()
		layout.addWidget(self.data_groupbox)
		layout.addWidget(self.validation_groupbox)
		
		self.setLayout(layout)

	def __init_iso_model(self):
		if self.iso_model == None:
			self.viewer_window.add_isosurface()
			self.iso_model = self.viewer_window.viewables[-1]
			#self.iso_model = EMIsosurfaceModel(self.viewer_window, None, enable_file_browse=False)
			#self.viewer_window.add_model(self.iso_model)
			#self.viewer_window.get_inspector().update_selection()
			#self.__set_model_contexts(self.iso_model)

	def __init_pdb_model(self):
		if self.pdb_model == None:
			self.pdb_model = EMPDBModel(self.viewer_window)
			self.viewer_window.add_model(self.pdb_model)
			self.__set_model_contexts(self.pdb_model)

	def __set_model_contexts(self,model):
		model.set_gl_widget(self.viewer_window)
		model.set_dont_delete_parent() # stops a RunTimeError
	
	def browse_iso(self):
		em_selector = EMSelectorDialog()
		file_path = em_selector.exec_()
		get_application().detach_child(em_selector)
		self.volume_line_edit.setText(file_path)
		
	def browse_pdb(self):
		file_path = QtGui.QFileDialog.getOpenFileName(filter="Protein Data Bank (*.pdb)")
		self.pdb_line_edit.setText(file_path)
		
	def closeEvent(self, event):
		self.viewer_window.close()
		QtGui.QWidget.closeEvent(self, event)
		
	def draw_objects(self):
		if self.iso_model == None: 
			self. __init_iso_model()
		if self.pdb_model == None:
			self.__init_pdb_model()
		if self.pdb_model != None:
			glPushMatrix()
			self.pdb_model.draw_objects()
			glPopMatrix()
		if self.iso_model != None:
			glPushMatrix()
			self.iso_model.render()
			glPopMatrix()

	def run_validate(self): 
		num_transformations = self.transformations_spinbox.value()
		threshold = self.threshold_doublespinbox.value()
		current_pdb = str(self.pdb_line_edit.text())
		current_mrc = str(self.volume_line_edit.text())

		self.emit(QtCore.SIGNAL("run_validate"), current_mrc, current_pdb, num_transformations, threshold)
		
	def update_iso_file(self):
		iso_file_path = str(self.volume_line_edit.text())
		data = EMData(iso_file_path)
		self.viewer_window.set_data(data, iso_file_path)
		self.iso_model = self.viewer_window.viewables[0]
#		if self.iso_model == None: 
#			self. __init_iso_model()
#		self.iso_model.set_data(data)
		self.viewer_window.updateGL()
		
#		self.iso_model.get_inspector().mrc_text.setText(file_name) #Use self.iso_model.data["source_path"] instead

	def update_pdb_file(self):
		pdb_file = str(self.pdb_line_edit.text())
		if self.pdb_model == None:
			self.__init_pdb_model()
		self.pdb_model.set_current_text(pdb_file) #updates GL with the new pdb file


if __name__ == '__main__':
	em_app = EMApp()
	pdb_widget = EMPDBValWidget()

	pdb_widget.volume_line_edit.setText("rdv-target2.mrc")
	pdb_widget.pdb_line_edit.setText("fh-solution-0-1UF2-T.pdb")
	
	em_app.show()
	em_app.execute()
