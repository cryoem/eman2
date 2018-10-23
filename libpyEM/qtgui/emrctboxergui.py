#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
#
# Author: John Flanagan, 04/08/2011 (jfflanag@bcm.edu)
# Edited by: Stephen Murray (scmurray@bcm.edu) May 2014
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
from builtins import range
from EMAN2 import get_image_directory, dump_processors_list
from PyQt4 import QtCore, QtGui
from .emrctstrategy import Strategy2IMGMan, Strategy2IMGPair
from EMAN2jsondb import js_open_dict
from EMAN2 import *
import os

class ControlPannel(QtWidgets.QWidget):
	'''This controls the RCT boxer. Normally this will not need to be midified. If a new pair pciking strategy is to be implmented, then
	A new GUI class should be added as decribed below and the new tool needs to be added to __init__ and the functions: current_tool_combobox_changed 
	and add_picker_tools'''
	def __init__(self, mediator):
		QtWidgets.QWidget.__init__(self)
		self.mediator = mediator
		self.db = js_open_dict("info/emboxerrctgui.json")
		self.qualitydb = js_open_dict("e2boxercache/quality.json")
		self.setWindowIcon(QtGui.QIcon(get_image_directory() +"green_boxes.png"))
		self.setWindowTitle("e2RCTboxer")
		
		# Here is where additional tools can be added
		self.manual_tool = ManualPicker(self.mediator, self.db)
		self.pair_picker_tool = PairPickerTool(self.mediator, self.db)
		
		vbox = QtWidgets.QVBoxLayout(self)
		
		self.tab_widget = QtWidgets.QTabWidget()
		self.tab_widget.addTab(self.get_main_tab(),"Main")
		self.tab_widget.addTab(self.get_processor_tab(),"Processor")
		self.tab_widget.addTab(self.get_filter_tab(),"Filter")
		vbox.addWidget(self.tab_widget)
		

		self.add_controls(vbox)			# add done button, among other things
		
		self.setLayout(vbox)
		self.add_picker_tools()
		self.updateGeometry()

		# Initialize tools
		self.tools_stacked_widget.setCurrentIndex(self.db.getdefault("toolidx",dfl=0))
		self.current_tool_combobox.setCurrentIndex(self.db.getdefault("toolidx",dfl=0))
		E2loadappwin("e2rctboxer","main",self)
#		E2loadappwin("e2rctboxer","maintab",self.wimage)
#		E2loadappwin("e2rctboxer","processortab",self.wfft)
#		E2loadappwin("e2rctboxer","filtertab",self.wplot)
		
	def get_main_tab(self):
		mainwidget = QtWidgets.QWidget()
		vbox = QtWidgets.QVBoxLayout()
		
		# Make the main tools layout
		mlayout = QtWidgets.QVBoxLayout()
		self.get_main(mlayout)			# add the widgets const for all tools
		msplitter = QtWidgets.QSplitter(QtCore.Qt.Vertical)
		msplitter.setFrameShape(QtWidgets.QFrame.StyledPanel)
		msplitter.setLayout(mlayout)
		vbox.addWidget(msplitter)
		
		# Make the tools layout
		self.add_boxing_button_group(vbox)	# add tool specific widgets		

		mainwidget.setLayout(vbox)
		
		return mainwidget
	
	def get_filter_tab(self):
		
		filterwidget = QtWidgets.QWidget()
		vbox = QtWidgets.QVBoxLayout()
		
		hbl=QtWidgets.QHBoxLayout()
		flabel = QtWidgets.QLabel("Filter Type:",self)
		hbl.addWidget(flabel)
		self.filter_combobox = QtWidgets.QComboBox()
		hbl.addWidget(self.filter_combobox)
		vbox.addLayout(hbl)
		
		hbl=QtWidgets.QHBoxLayout()
		klabel = QtWidgets.QLabel("Kernal Size:",self)
		hbl.addWidget(klabel)
		self.kernel_combobox = QtWidgets.QComboBox()
		hbl.addWidget(self.kernel_combobox)
		vbox.addLayout(hbl)
		
		self.kernel_stacked_widget = QtWidgets.QStackedWidget()
		vbox.addWidget(self.kernel_stacked_widget)
		
		self.filter_but=QtWidgets.QPushButton("Filter")
		vbox.addWidget(self.filter_but)
		
		filterwidget.setLayout(vbox)
		
		self.gridkernel = []
		self.gridkernel.append(self.add_custom_kernels(3))
		self.gridkernel.append(self.add_custom_kernels(5))
		self.gridkernel.append(self.add_custom_kernels(7))
		self.add_kernel_sizes()
		
		self.filter_but.clicked[bool].connect(self.on_filter)
		self.kernel_combobox.activated[int].connect(self.kernel_combobox_changed)
		self.filter_combobox.activated[int].connect(self.filter_combobox_changed)
		
		return filterwidget
	
	def get_processor_tab(self):
		
		processorwidget = QtWidgets.QWidget()
		vboxa = QtWidgets.QVBoxLayout()
		
		vbox1 = QtWidgets.QVBoxLayout()
		hbl=QtWidgets.QHBoxLayout()
		flabel = QtWidgets.QLabel("Filter:",self)
		hbl.addWidget(flabel)
		self.processor_combobox = QtWidgets.QComboBox()
		proc_data = dump_processors_list()
		for key in list(proc_data.keys()):
			if len(key) >= 5 and key[:7] == "filter.":
				#print key
				self.processor_combobox.addItem(key)
		self.processor_combobox.setCurrentIndex(int(self.db.getdefault("processor",dfl=0)))
		hbl.addWidget(self.processor_combobox)
		vbox1.addLayout(hbl)
		
		hbl2=QtWidgets.QHBoxLayout()
		plabel = QtWidgets.QLabel("Parameters:",self)
		hbl2.addWidget(plabel)
		self.params_listbox = QtWidgets.QLineEdit(str(self.db.getdefault("processorparams",dfl="")), self)
		hbl2.addWidget(self.params_listbox)
		vbox1.addLayout(hbl2)
		vbox1.setAlignment(QtCore.Qt.AlignTop)
		vboxa.addLayout(vbox1)
		
		vbox2 = QtWidgets.QVBoxLayout()
		self.processor_but=QtWidgets.QPushButton("Filter")
		vbox2.addWidget(self.processor_but)
		vboxa.addLayout(vbox2)
		
		processorwidget.setLayout(vboxa)
		
		self.processor_but.clicked[bool].connect(self.on_processor)
		self.processor_combobox.activated[int].connect(self.processor_combobox_changed)
		self.params_listbox.editingFinished.connect(self.params_listbox_changed)
		
		return processorwidget
	
	def on_processor(self):
		filtertype = str(self.processor_combobox.currentText())
		parstring = str(self.params_listbox.text())
		pars = parstring.split(":")
		pardict = {}
		for par in pars:
			parpair = par.split("=")
			pardict[parpair[0]] = float(parpair[1])
			
		for window in self.mediator.windowlist:
			fdata = window.data.process(filtertype,pardict)
			print("filtered win")
			window.reload_image(fdata, window.filename+"_filted")
		
	def processor_combobox_changed(self, idx):
		self.db["processor"] = idx
		self.processor_combobox.setCurrentIndex(idx)
		
	def params_listbox_changed(self):
		self.db["processorparams"] = self.params_listbox.text()
		
	def on_filter(self):
		
		kernel = self.gridkernel[self.kernel_combobox.currentIndex()]
		kt = []
		for i in kernel:
			if i.text():
				kt.append(float(i.text()))
			else:
				kt.append(0.0)
		#filter each window
		for window in self.mediator.windowlist:
			fdata = window.data.process("filter.convolution.kernel",{"kernel":kt})
			print("filtered win")
			window.reload_image(fdata, window.filename+"_filted")

	
	def add_kernel_sizes(self):
		self.kernel_combobox.addItem("3x3")
		self.kernel_combobox.addItem("5x5")
		self.kernel_combobox.addItem("7x7")
		# Load data
		idx = self.db.getdefault("kernelsizeidx",dfl=0)
		self.kernel_combobox.setCurrentIndex(idx)
		self.kernel_combobox_changed(idx)
		
	def kernel_combobox_changed(self, idx):
		self.kernel_stacked_widget.setCurrentIndex(idx)
		self.db["kernelsizeidx"] = idx
		self.load_filters(idx)
	
	def load_filters(self, idx):
		# remove old items
		for i in range(self.filter_combobox.count()):
			self.filter_combobox.removeItem(i)
		# add new ones
		# 3x3
		if idx == 0:
			self.filter_combobox.addItem("Lowpass Rect")
		# 5x5
		if idx == 1:
			self.filter_combobox.addItem("Lowpass Rect")
		# 7x7
		if idx == 2:
			self.filter_combobox.addItem("Lowpass Rect")
			
		self.filter_combobox_changed(idx)
			
	def filter_combobox_changed(self, idx):
		# 3x3 Filters
		if self.kernel_combobox.currentIndex() == 0:
			if self.filter_combobox.currentText() == "Lowpass Rect":
				for i in range(9):
					self.gridkernel[0][i].setText("1")
		# 5x5 Filters
		if self.kernel_combobox.currentIndex() == 1:
			if self.filter_combobox.currentText() == "Lowpass Rect":
				for i in range(25):
					self.gridkernel[1][i].setText("1")	
		# 7x7 Filters
		if self.kernel_combobox.currentIndex() == 2:
			if self.filter_combobox.currentText() == "Lowpass Rect":
				for i in range(49):
					self.gridkernel[2][i].setText("1")	
			
	def add_custom_kernels(self, size):
		self.kernelwidget = QtWidgets.QWidget()
		grid3 = QtWidgets.QGridLayout()
		kernelwidgetidx = []
		for i in range(size):
			for j in range(size):
				kw = QtWidgets.QLineEdit("",self)
				kw.setFixedSize(40,25)	# This could be an issue......
				kernelwidgetidx.append(kw)
				grid3.addWidget(kw, i, j)
		
		grid3.setAlignment(QtCore.Qt.AlignCenter)
		self.kernelwidget.setLayout(grid3)
		self.kernel_stacked_widget.addWidget(self.kernelwidget)
		
		return kernelwidgetidx
		
	def get_main(self, layout):
		hbl=QtWidgets.QHBoxLayout()
		self.box_size_label = QtWidgets.QLabel("Box Size:",self)
		hbl.addWidget(self.box_size_label)
		
		self.pos_int_validator = QtGui.QIntValidator(2,5000, self)	#Anything bigger than 5,000 is crazy!!!!
		self.boxsize = QtWidgets.QLineEdit(str(self.mediator.boxsize),self)
		self.boxsize.setValidator(self.pos_int_validator)
		
		hbl.addWidget(self.boxsize)
		layout.addLayout(hbl)
		
		self.boxsize.editingFinished.connect(self.new_boxsize)
		
	def add_boxing_button_group(self,layout):
		self.tool_button_group_box = QtWidgets.QGroupBox("Tools")
		
		grid = QtWidgets.QGridLayout()
		self.current_tool_combobox = QtWidgets.QComboBox()
		grid.addWidget(QtWidgets.QLabel("Current Boxing Tool:"),0,0)
		grid.addWidget(self.current_tool_combobox,0,1)
		# Add stacked widget
		self.tools_stacked_widget = QtWidgets.QStackedWidget()
		grid.addWidget(self.tools_stacked_widget,1,0,1,2)
		# Add quality combobox
		self.quality = QtWidgets.QComboBox()
		for i in range(5): self.quality.addItem(str(i))
		# check full path then check basename
		if self.mediator.windowlist[0].filename not in self.qualitydb:
			self.quality.setCurrentIndex(self.qualitydb.getdefault(self.mediator.windowlist[0].filename,dfl=0))
		else:
			self.quality.setCurrentIndex(self.qualitydb.getdefault(os.path.basename(self.mediator.windowlist[0].filename),dfl=0))
		grid.addWidget(QtWidgets.QLabel("Quality Score"),2,0)
		grid.addWidget(self.quality, 2,1)
		# add to layout
		self.tool_button_group_box.setLayout(grid)
		layout.addWidget(self.tool_button_group_box,0,)
		
		self.current_tool_combobox.activated[int].connect(self.current_tool_combobox_changed)
		self.quality.activated[int].connect(self.quality_score_changed)
	
	def add_controls(self, layout):
		butbox = QtWidgets.QHBoxLayout()
		self.write_box_but=QtWidgets.QPushButton("Write Boxes")
		butbox.addWidget(self.write_box_but)
		self.write_but=QtWidgets.QPushButton("Write Ptcls")
		butbox.addWidget(self.write_but)
		layout.addLayout(butbox)
		self.done_but=QtWidgets.QPushButton("Done")
		layout.addWidget(self.done_but)
		
		self.write_box_but.clicked[bool].connect(self.on_write_box)
		self.write_but.clicked[bool].connect(self.on_write)
		self.done_but.clicked[bool].connect(self.on_done)
	
	# This function configures the tools up tool change
	def current_tool_combobox_changed(self, idx):
		self.tools_stacked_widget.setCurrentIndex(idx)
		self.db["toolidx"] = idx
		if self.current_tool_combobox.currentText() == "Manual":
			self.manual_tool.configure_widget()
			#print "Set strategy to Manual"
		if self.current_tool_combobox.currentText() == "Pair Picker":
			self.pair_picker_tool.configure_widget()
			#print "Set strategy to Pair Picker"
	
	def quality_score_changed(self, idx):
		for window in self.mediator.windowlist:
			self.qualitydb[window.filename] = idx
		
	# Here is where additional tools can be added
	def add_picker_tools(self):
		self.tools_stacked_widget.addWidget(self.manual_tool.get_widget())
		self.current_tool_combobox.addItem("Manual")
		self.tools_stacked_widget.addWidget(self.pair_picker_tool.get_widget())
		self.current_tool_combobox.addItem("Pair Picker")
	
	# This function configures the picking tools upon startup
	def configure_tools(self):
		self.manual_tool.configure_widget()
		self.pair_picker_tool.configure_widget()
		
	def new_boxsize(self):
		self.mediator.boxsize = int(self.boxsize.text())
		self.db["box_size"] = self.mediator.boxsize
		
		for window in self.mediator.windowlist:
			window.boxes.reset_images()
			window.boxes.reset_shapes()
			window.update_mainwin()
			window.update_particles()
	
	def closeEvent(self,event):
		E2saveappwin("e2rctboxer","main",self)
		self.on_done()
		
	def on_write(self):
		print("Saving Particles")
		for window in self.mediator.windowlist:
			splitpath = os.path.split(os.path.splitext(window.filename)[0])
			if splitpath[0] == '':
				window.write_particles(window.filename, ("particles/"+splitpath[1]),self.mediator.boxsize,normproc="normalize.edgemean")
			else:
				window.write_particles(window.filename, (splitpath[0]+"/particles/"+splitpath[1]),self.mediator.boxsize,normproc="normalize.edgemean")
	def on_write_box(self):
		for window in self.mediator.windowlist:
			window.write_boxes(os.path.splitext(window.filename)[0]+".box",self.mediator.boxsize)
		
	def on_done(self):
		for wid in self.mediator.widgetlist:
			if wid != None:
				wid.close()

# Current tools. Other tools can be added by simply adding a Pciker GUi and then building a 
# corresponding Strategy based by subclassing Strategy in emrctstrategy
class ManualPicker(QtWidgets.QWidget):
	def __init__(self, mediator, db):
		QtWidgets.QWidget.__init__(self)
		self.mediator = mediator
		self.db=db
		vbl = QtWidgets.QVBoxLayout()
		label = QtWidgets.QLabel("Manual Picker", self)
		boldfont = QtGui.QFont()
		boldfont.setBold(True)
		label.setFont(boldfont)
		label.setAlignment(QtCore.Qt.AlignTop)
		self.clr_but = QtWidgets.QPushButton("Clear", self)
		vbl.addWidget(label)
		vbl.addWidget(self.clr_but)
		self.setLayout(vbl)
		
		self.mpsplitter = QtWidgets.QSplitter(QtCore.Qt.Vertical)
		self.mpsplitter.setFrameShape(QtWidgets.QFrame.StyledPanel)
		self.mpsplitter.addWidget(self)
		self.clr_but.clicked[bool].connect(self.on_clear)
	
	def on_clear(self):
		for window in self.mediator.windowlist:
			window.boxes.clear_boxes()
			window.update_mainwin()
			
	def configure_widget(self):
		self.mediator.set_strategy(Strategy2IMGMan)
	
	def get_widget(self):	
		return self.mpsplitter
		
class PairPickerTool(QtWidgets.QWidget):
	def __init__(self, mediator, db):
		QtWidgets.QWidget.__init__(self)
		self.mediator = mediator
		self.db = db
		self.updateboxes = False
		self.minpp_for_xform = 3
		self.centertilts = False
		
		# GUI code below here
		ppwidget = QtWidgets.QWidget()
		
		vbl = QtWidgets.QVBoxLayout()
		label = QtWidgets.QLabel("Pair Picker", self)
		boldfont = QtGui.QFont()
		boldfont.setBold(True)
		label.setFont(boldfont)
		vbl.addWidget(label)

		self.updateboxes_cb = QtWidgets.QCheckBox("Update box positions")
		self.updateboxes_cb.setChecked(False)
		vbl.addWidget(self.updateboxes_cb)
		
		self.centertilts_cb = QtWidgets.QCheckBox("Center opposite box position")
		self.centertilts_cb.setChecked(False)
		vbl.addWidget(self.centertilts_cb)
		
		hbl = QtWidgets.QHBoxLayout()
		slabel = QtWidgets.QLabel("Min pairs for xform", self)
		hbl.addWidget(slabel)
		self.spinbox = QtWidgets.QSpinBox(self)
		self.spinbox.setMinimum(self.minpp_for_xform)
		self.spinbox.setMaximum(1000)
		hbl.addWidget(self.spinbox)
		vbl.addLayout(hbl)
		
		hta = QtWidgets.QHBoxLayout()
		tlabel = QtWidgets.QLabel("Computed tilt angle", self)
		hta.addWidget(tlabel)
		self.tiltangle = QtWidgets.QLineEdit("", self)
		self.tiltangle.setReadOnly(True)
		hta.addWidget(self.tiltangle)
		vbl.addLayout(hta)
		
		htax = QtWidgets.QHBoxLayout()
		talabel = QtWidgets.QLabel("Computed tilt axis (Y)", self)
		htax.addWidget(talabel)
		self.tiltaxis = QtWidgets.QLineEdit("", self)
		self.tiltaxis.setReadOnly(True)
		htax.addWidget(self.tiltaxis)
		vbl.addLayout(htax)
		
		hgamma = QtWidgets.QHBoxLayout()
		gammalabel = QtWidgets.QLabel("Gamma angle", self)
		hgamma.addWidget(gammalabel)
		self.gamma = QtWidgets.QLineEdit("", self)
		self.gamma.setReadOnly(True)
		hgamma.addWidget(self.gamma)
		vbl.addLayout(hgamma)
		
		hmb = QtWidgets.QHBoxLayout()
		mlabel = QtWidgets.QLabel("Mask Type", self)
		hmb.addWidget(mlabel)
		self.mask_combobox = QtWidgets.QComboBox()
		self.mask_combobox.setEnabled(False)
		hmb.addWidget(self.mask_combobox)
		vbl.addLayout(hmb)
		
		hbb = QtWidgets.QHBoxLayout()
		self.upboxes_but = QtWidgets.QPushButton("Update Boxes", self)
		self.upboxes_but.setEnabled(False)
		hbb.addWidget(self.upboxes_but)
		
		self.centerboxes_but = QtWidgets.QPushButton("Center Boxes", self)
		self.centerboxes_but.setEnabled(False)
		hbb.addWidget(self.centerboxes_but)
		vbl.addLayout(hbb)
		
		self.clr_but = QtWidgets.QPushButton("Clear", self)
		vbl.addWidget(self.clr_but)
		self.setLayout(vbl)
		
		self.ppsplitter = QtWidgets.QSplitter(QtCore.Qt.Vertical)
		self.ppsplitter.setFrameShape(QtWidgets.QFrame.StyledPanel)
		self.ppsplitter.addWidget(self)
		
		self.spinbox.valueChanged[int].connect(self.on_spinbox)
		self.updateboxes_cb.stateChanged[int].connect(self.on_updateboxes)
		self.centertilts_cb.stateChanged[int].connect(self.on_centertilts)
		self.clr_but.clicked[bool].connect(self.on_clear)
		self.centerboxes_but.clicked[bool].connect(self.on_centerboxes_but)
		self.upboxes_but.clicked[bool].connect(self.on_upboxes_but)
	
		# Initialize
		self.spinbox.setValue(self.db.getdefault("ppspinbox",dfl=self.minpp_for_xform))
		self.updateboxes_cb.setChecked(self.db.getdefault("ppcheckbox",dfl=self.updateboxes))
		self.centertilts_cb.setChecked(self.db.getdefault("tiltcheckbox",dfl=self.centertilts))
		self.addmasks()
	
	def addmasks(self):
		self.mask_combobox.addItem("None")
		self.mask_combobox.addItem("LineMask")
		self.mask_combobox.addItem("SolidMask")
		self.mask_combobox.setCurrentIndex(self.db.getdefault("masktype",dfl=0))
		
		self.mask_combobox.activated[int].connect(self.mask_combobox_changed)
		
	def mask_combobox_changed(self, idx):
		self.db["masktype"] = idx
		self.mediator.untilt_win.masktype = self.mask_combobox.currentText()
		self.mediator.tilt_win.masktype = self.mask_combobox.currentText()
		self.mediator.strategy.compute_mask()
	
	# This function is called upon widget startup
	def configure_widget(self):
		self.mediator.set_strategy(Strategy2IMGPair)
		self.mediator.strategy.initial_calculations()
		self.mask_combobox_changed(self.mask_combobox.currentIndex())
		self.mediator.configure_strategy(self)
			
	def on_spinbox(self, value):
		self.db["ppspinbox"] = value
		self.minpp_for_xform = value
		self.mediator.configure_strategy(self)
		
	def on_updateboxes(self):
		self.db["ppcheckbox"] = self.updateboxes_cb.isChecked()
		self.updateboxes = self.updateboxes_cb.isChecked()
		self.mediator.configure_strategy(self)
	
	def on_centertilts(self):
		self.db["tiltcheckbox"] = self.centertilts_cb.isChecked()
		self.centertilts = self.centertilts_cb.isChecked()
		self.mediator.configure_strategy(self)
		
	def on_clear(self):
		for window in self.mediator.windowlist:
			window.boxes.clear_boxes()
			window.update_mainwin()	
			
	def on_upboxes_but(self):
		self.mediator.handle_strategy_signal("updateboxes")
		
	def on_centerboxes_but(self):
		self.mediator.handle_strategy_signal("centerboxes")
	
	def get_widget(self):
		return self.ppsplitter
