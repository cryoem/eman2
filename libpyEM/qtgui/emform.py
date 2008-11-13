#!/usr/bin/env python

#
# Author: David Woolford 11/7/2008 (woolford@bcm.edu)
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


from emdatastorage import ParamDef
from PyQt4 import QtGui,QtCore
from PyQt4.QtCore import Qt
import os
from emselector import EMSelectorModule
from emapplication import EMQtWidgetModule
from EMAN2 import Util


class ParamTable(list):
	'''
	This is just a list of ParamDef objects
	'''
	def __init__(self,name=None,desc_short=None,desc_long=""):
		list.__init__(self)
		self.vartype = "paramtable"
		self.name = name
		self.desc_short = desc_short
		self.desc_long = desc_long
		
class EMFormModule(EMQtWidgetModule):
	'''
	params should be a list of ParamDef objects
	application should be an EMAN2 type application
	After initializing this object you should probably call application.show(this)
	or application.show_specific(this), depending on what you're doing
	'''
	def __init__(self,params,application):
		self.application = application
		self.widget = EMFormWidget(self,params)
		EMQtWidgetModule.__init__(self,self.widget,application)
		
	def get_desktop_hint(self):
		return "form"
		

class EMFormWidget(QtGui.QWidget):
	'''
	See the example in __main__ below
	If ok is clicked the "emform_ok" signal is emitted along with a dictionary containing all of the form entries
	If cancel is clicked the "emform_cancel" signal is emmitted. No extra information is sent in this case
	'''
	def __init__(self,parent,params=None):
		QtGui.QWidget.__init__(self,None)
		self.parent = parent
		self.params = params
		self.event_handlers = [] # used to keep event handlers in memory
		self.output_writers = [] # used to register output write objects, for the purpose of returning the results
		
		self.setWindowIcon(QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/eman.png"))
		
		self.__auto_incorporate = {}
		self.__auto_incorporate["float"] = self.__incorporate_float
		self.__auto_incorporate["int"] = self.__incorporate_int
		self.__auto_incorporate["url"] = self.__incorporate_url
		self.__auto_incorporate["choice"] = self.__incorporate_choice
		self.__auto_incorporate["string"] = self.__incorporate_string
		self.__auto_incorporate["text"] = self.__incorporate_text
		self.__auto_incorporate["boolean"] = self.__incorporate_boolean
		self.__auto_incorporate["stringlist"] = self.__incorporate_stringlist
		self.__auto_incorporate["intlist"] = self.__incorporate_intlist
		self.__auto_incorporate["floatlist"] = self.__incorporate_floatlist
		self.__auto_incorporate["paramtable"] = self.__incorporate_paramtable
		
		self.vbl = QtGui.QVBoxLayout(self)
		self.__incorporate_params(self.params,self.vbl)
		self.__add_ok_cancel_buttons(self.vbl)
	
	def __incorporate_params(self,params,layout):
		for param in self.params:
			try:
				if len(param) != 1 and not isinstance(param,ParamTable):
					hbl=QtGui.QHBoxLayout()
					for iparam in param:
						self.__auto_incorporate[iparam.vartype](iparam,hbl)
					layout.addLayout(hbl)
					continue
					
			except: pass
			self.__auto_incorporate[param.vartype](param,layout)
	
	def __incorporate_paramtable(self,paramtable,layout):
		
		num_choices = None
		# first check that there are no inconsistencies in the number of parameter choices
		for param in paramtable:
			if num_choices == None:
				num_choices = len(param.choices)
			else:
				if len(param.choices) != num_choices:
					print "error, the number of choices is not consistent in __incorporate_paramtable"
					return
		
		hbl=QtGui.QHBoxLayout()
		hbl.setMargin(0)
		hbl.setSpacing(2)
		
		table_widget = QtGui.QTableWidget(num_choices, len(paramtable), None)
		#self.items = []
		table_widget.setSortingEnabled(False)
		for i,param in enumerate(paramtable):
			for j,choice in enumerate(param.choices):
				item = QtGui.QTableWidgetItem(str(choice))
				flag2 = Qt.ItemFlags(Qt.ItemIsSelectable)
				flag3 = Qt.ItemFlags(Qt.ItemIsEnabled)
				flag4 = Qt.ItemFlags(Qt.ItemIsEditable)
				#flags = flags.
				if i == 0:
					item.setFlags(flag2|flag3|flag4)
				else:
					item.setFlags(flag3|flag4)
				table_widget.setItem(j, i, item)
				
				
			item = QtGui.QTableWidgetItem(param.desc_short)
			table_widget.setHorizontalHeaderItem(i,item)
			
		table_widget.setSortingEnabled(True)
		hbl.addWidget(table_widget)
		
		groupbox = QtGui.QGroupBox(paramtable.desc_short)
		groupbox.setToolTip(paramtable.desc_long)
		groupbox.setLayout(hbl)
		layout.addWidget(groupbox)
		
		self.output_writers.append(ParamTableWriter(paramtable.name,table_widget,type(paramtable[0].choices[0])))
		
	def __incorporate_floatlist(self,param,layout):
		hbl=QtGui.QHBoxLayout()
		hbl.setMargin(0)
		hbl.setSpacing(2)
		
		list_widget = QtGui.QListWidget(None)
		
		list_widget.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
		list_widget.setMouseTracking(True)	
	
		for choice in param.choices:
			a = QtGui.QListWidgetItem(str(choice),list_widget)
			if choice in param.defaultunits:
				list_widget.setItemSelected(a,True)
			
		hbl.addWidget(list_widget)
		groupbox = QtGui.QGroupBox(param.desc_short)
		groupbox.setToolTip(param.desc_long)
		groupbox.setLayout(hbl)
		
		layout.addWidget(groupbox)

		self.output_writers.append(FloatListParamWriter(param.name,list_widget))
	
	def __incorporate_stringlist(self,param,layout):
		return self.__incorporate_list(param,layout,str)
	
	def __incorporate_floatlist(self,param,layout):
		return self.__incorporate_list(param,layout,float)
	
	def __incorporate_intlist(self,param,layout):
		return self.__incorporate_list(param,layout,int)
	
	def __incorporate_list(self,param,layout,type_of):
		hbl=QtGui.QHBoxLayout()
		hbl.setMargin(0)
		hbl.setSpacing(2)
		
		list_widget = QtGui.QListWidget(None)
		
		list_widget.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
		list_widget.setMouseTracking(True)	
	
		for choice in param.choices:
			a = QtGui.QListWidgetItem(str(choice),list_widget)
			if choice in param.defaultunits:
				list_widget.setItemSelected(a,True)
			
		hbl.addWidget(list_widget)
		groupbox = QtGui.QGroupBox(param.desc_short)
		groupbox.setToolTip(param.desc_long)
		groupbox.setLayout(hbl)
		
		layout.addWidget(groupbox)

		self.output_writers.append(ListWidgetParamWriter(param.name,list_widget,type_of))
		#layout.addLayout(hbl)
	
	def __incorporate_boolean(self,param,layout):
		hbl=QtGui.QHBoxLayout()
		hbl.setMargin(0)
		hbl.setSpacing(2)
		check_box = QtGui.QCheckBox(str(param.desc_short),self)
		check_box.setChecked(bool(param.defaultunits))
		hbl.addWidget(check_box)
		layout.addLayout(hbl)
		self.output_writers.append(BoolParamWriter(param.name,check_box))
	
	def __incorporate_string(self,param,layout):
		if param.choices != None and len(param.choices) > 1:
			self.__incorporate_string_with_choices(param,layout)
		else:
			hbl=QtGui.QHBoxLayout()
			hbl.setMargin(0)
			hbl.setSpacing(2)
			label = QtGui.QLabel(param.desc_short+":",self)
			label.setToolTip(param.desc_long)
			hbl.addWidget(label)
			line_edit = QtGui.QLineEdit(str(param.defaultunits),self)
			hbl.addWidget(line_edit)
			hbl.name = param.name
			layout.addLayout(hbl)
			self.output_writers.append(StringParamWriter(param.name,line_edit))
		
	def __incorporate_float(self,param,layout):
		if param.choices != None and len(param.choices) > 1:
			self.__incorporate_float_with_choices(param,layout)
		else:
			hbl=QtGui.QHBoxLayout()
			hbl.setMargin(0)
			hbl.setSpacing(2)
			label = QtGui.QLabel(param.desc_short+":",self)
			label.setToolTip(param.desc_long)
			hbl.addWidget(label)
			double_validator = QtGui.QDoubleValidator(self)
			line_edit = QtGui.QLineEdit(str(param.defaultunits),self)
			line_edit.setValidator(double_validator)
			hbl.addWidget(line_edit)
			hbl.name = param.name
			layout.addLayout(hbl)
			self.output_writers.append(FloatParamWriter(param.name,line_edit))
		
	def __incorporate_int(self,param,layout):
		if param.choices != None and len(param.choices) > 1:
			self.__incorporate_int_with_choices(param,layout)
		else:
			hbl=QtGui.QHBoxLayout()
			hbl.setMargin(0)
			hbl.setSpacing(2)
			label = QtGui.QLabel(param.desc_short+":",self)
			label.setToolTip(param.desc_long)
			hbl.addWidget(label)
			pos_int_validator = QtGui.QIntValidator(self)
			line_edit = QtGui.QLineEdit(str(param.defaultunits),self)
			line_edit.setValidator(pos_int_validator)
			hbl.addWidget(line_edit)
			hbl.name = param.name
			layout.addLayout(hbl)
			self.output_writers.append(IntParamWriter(param.name,line_edit))
		
	def __incorporate_text(self,param,layout):
		hbl=QtGui.QHBoxLayout()
		hbl.setMargin(0)
		hbl.setSpacing(2)
				
		text_edit = QtGui.QTextEdit("",self)
		text_edit.setReadOnly(True)
		text_edit.setWordWrapMode(QtGui.QTextOption.WordWrap)
		text_edit.setText(param.defaultunits)
		hbl.addWidget(text_edit)
		
		groupbox = QtGui.QGroupBox(param.desc_short)
		groupbox.setToolTip(param.desc_long)
		groupbox.setLayout(hbl)
		
		layout.addWidget(groupbox)

		self.output_writers.append(TextParamWriter(param.name,text_edit))	
	
	def __incorporate_url(self,param,layout):
		vbl=QtGui.QVBoxLayout()
		hbl=QtGui.QHBoxLayout()
		hbl.setMargin(0)
		hbl.setSpacing(2)
		defaults = ""
		if param.defaultunits != None:
			for i,d in enumerate(param.defaultunits):
				defaults += d
				if i != (len(param.defaultunits)-1): 
					defaults += '\n'
				
		text_edit = QtGui.QTextEdit("",self)
		text_edit.setWordWrapMode(QtGui.QTextOption.NoWrap)
		text_edit.setText(defaults)
		hbl.addWidget(text_edit)
		vbl.addLayout(hbl)
		
		hbl2=QtGui.QHBoxLayout()
		hbl2.setMargin(0)
		hbl2.setSpacing(2)
		
		browse_button = QtGui.QPushButton("Browse",self)
		hbl2.addWidget(browse_button)
		clear_button = QtGui.QPushButton("Clear",self)
		hbl2.addWidget(clear_button)
		vbl.addLayout(hbl2)
		
		groupbox = QtGui.QGroupBox(param.desc_short)
		groupbox.setToolTip(param.desc_long)
		groupbox.setLayout(vbl)
		
		layout.addWidget(groupbox)
		
		self.event_handlers.append(UrlEventHandler(self,text_edit,browse_button,clear_button,self.parent.application,param.desc_short))
		self.output_writers.append(UrlParamWriter(param.name,text_edit))

	def __incorporate_choice(self,param,layout):
		
		hbl = QtGui.QHBoxLayout()
		buttons = []
		for choice in param.choices:
			button = QtGui.QRadioButton(str(choice))
			if choice == param.defaultunits: button.setChecked(True)
			hbl.addWidget( button)
			buttons.append(button)
		
		groupbox = QtGui.QGroupBox(param.desc_short)
		groupbox.setToolTip(param.desc_long)
		groupbox.setLayout(hbl)
		layout.addWidget(groupbox)
		self.output_writers.append(ChoiceParamWriter(param.name,buttons,type(param.choices[0])))
	
	def __incorporate_float_with_choices(self,param,layout):
		hbl=QtGui.QHBoxLayout()
		hbl.setMargin(0)
		hbl.setSpacing(2)
		
		hbl=QtGui.QHBoxLayout()
		label = QtGui.QLabel(param.desc_short,self)
		label.setToolTip(param.desc_long)
		hbl.addWidget(label)
				
		combo = QtGui.QComboBox(self)
		idx_default = 0
		for i,float_p in enumerate(param.choices):
			combo.addItem(str(float_p))
			
			if str(float_p) == str(param.defaultunits):
				idx_default = i
		combo.setCurrentIndex(idx_default)

		hbl.addWidget(combo)
		layout.addLayout(hbl)
		
		self.output_writers.append(FloatChoiceParamWriter(param.name,combo))
	
	def __incorporate_int_with_choices(self,param,layout):
		hbl=QtGui.QHBoxLayout()
		hbl.setMargin(0)
		hbl.setSpacing(2)
		
		hbl=QtGui.QHBoxLayout()
		label = QtGui.QLabel(param.desc_short,self)
		label.setToolTip(param.desc_long)
		hbl.addWidget(label)
				
		combo = QtGui.QComboBox(self)
		idx_default = 0
		for i,integer in enumerate(param.choices):
			combo.addItem(str(integer))
			if str(integer) == str(param.defaultunits):
				idx_default = i
		combo.setCurrentIndex(idx_default)

		hbl.addWidget(combo)
		layout.addLayout(hbl)
		
		self.output_writers.append(IntChoiceParamWriter(param.name,combo))
	
	def __incorporate_string_with_choices(self,param,layout):
		hbl=QtGui.QHBoxLayout()
		hbl.setMargin(0)
		hbl.setSpacing(2)
		
		hbl=QtGui.QHBoxLayout()
		label = QtGui.QLabel(param.desc_short,self)
		label.setToolTip(param.desc_long)
		hbl.addWidget(label)
				
		combo = QtGui.QComboBox(self)
		idx_default = 0
		for i,string in enumerate(param.choices):
			combo.addItem(string)
			if string == param.defaultunits:
				idx_default = i
		combo.setCurrentIndex(idx_default)

		hbl.addWidget(combo)
		layout.addLayout(hbl)
		
		self.output_writers.append(StringChoiceParamWriter(param.name,combo))
	
	def __add_ok_cancel_buttons(self,layout):
		hbl=QtGui.QHBoxLayout()
		label = QtGui.QLabel("Form commands:",self)
		hbl.addWidget(label)
		
		ok_button = QtGui.QPushButton("Ok",self)
		ok_button.setToolTip("When you click ok the values in the form are sent to the calling program in a dictionary")
		hbl.addWidget(ok_button)
		cancel_button = QtGui.QPushButton("Cancel",self)
		hbl.addWidget(cancel_button)
		layout.addLayout(hbl)
		QtCore.QObject.connect(ok_button,QtCore.SIGNAL("clicked(bool)"),self.ok_pressed)
		QtCore.QObject.connect(cancel_button,QtCore.SIGNAL("clicked(bool)"),self.cancel_pressed)
		
	def ok_pressed(self,bool):
		ret = {}
		for output in self.output_writers: output.write_data(ret)
		self.parent.emit(QtCore.SIGNAL("emform_ok"),ret) # getting the parent to emit ensures integration with the desktop
		
	def cancel_pressed(self,bool):
		self.parent.emit(QtCore.SIGNAL("emform_cancel")) # getting the parent to emit ensures integration with the desktop


	def update_texture(self):
		self.parent.force_texture_update()
		
	def closeEvent(self,event):
		self.emit(QtCore.SIGNAL("emform_close"))

		
class ParamTableWriter:
	def __init__(self,param_name,table_widget,type_of):
		self.param_name = param_name
		self.table_widget = table_widget
		self.type_of = type_of
		
	def write_data(self,dict):
		sel = [self.type_of(item.text()) for item in self.table_widget.selectedItems()]
		print sel,self.param_name
		print "\n\n\n"
		dict[self.param_name] = sel

class BoolParamWriter:
	def __init__(self,param_name,check_box):
		self.param_name = param_name
		self.check_box = check_box
		
	def write_data(self,dict):
		dict[self.param_name] = bool(self.check_box.isChecked())

class FloatChoiceParamWriter:
	def __init__(self,param_name,combo):
		self.param_name = param_name
		self.combo = combo
		
	def write_data(self,dict):
		dict[self.param_name] = float(self.combo.currentText())

class IntChoiceParamWriter:
	def __init__(self,param_name,combo):
		self.param_name = param_name
		self.combo = combo
		
	def write_data(self,dict):
		dict[self.param_name] = int(self.combo.currentText())

class StringChoiceParamWriter:
	def __init__(self,param_name,combo):
		self.param_name = param_name
		self.combo = combo
		
	def write_data(self,dict):
		dict[self.param_name] = str(self.combo.currentText())


class ListWidgetParamWriter:
	def __init__(self,param_name,list_widget,type_of):
		self.param_name = param_name
		self.list_widget = list_widget
		self.type_of = type_of
		
	def write_data(self,dict):
		choices = []
		for a in self.list_widget.selectedItems():
			choices.append(self.type_of(a.text()))
			
		dict[self.param_name] = choices

class StringParamWriter:
	def __init__(self,param_name,line_edit):
		self.param_name = param_name
		self.line_edit = line_edit
		
	def write_data(self,dict):
		dict[self.param_name] = str(self.line_edit.text())


class FloatParamWriter:
	def __init__(self,param_name,line_edit):
		self.param_name = param_name
		self.line_edit = line_edit
		
	def write_data(self,dict):
		dict[self.param_name] = float(self.line_edit.text())

class IntParamWriter:
	def __init__(self,param_name,line_edit):
		self.param_name = param_name
		self.line_edit = line_edit
		
	def write_data(self,dict):
		dict[self.param_name] = int(self.line_edit.text())
		
class TextParamWriter:
	def __init__(self,param_name,text_edit):
		self.param_name = param_name
		self.text_edit = text_edit
		
	def write_data(self,dict):
		dict[self.param_name] = str(self.text_edit.toPlainText())
		
class UrlParamWriter:
	def __init__(self,param_name,text_edit):
		self.param_name = param_name
		self.text_edit = text_edit
		
	def write_data(self,dict):
		strings = [i for i in str(self.text_edit.toPlainText()).split()]
		dict[self.param_name] = strings

		
class ChoiceParamWriter:
	def __init__(self,param_name,list_radio_buttons,correct_type):
		self.param_name = param_name
		self.list_radio_buttons = list_radio_buttons
		self.correct_type = correct_type
	def write_data(self,dict):
		choice = None
		for button in self.list_radio_buttons:
			if button.isChecked():
				choice = self.correct_type(str(button.text()))
		dict[self.param_name] = choice

class UrlEventHandler:
	'''
	The browse and cancel events have to be sent to the correct line edit, so this handles it
	Basically to simplify things when there is more than one url type.
	'''
	def __init__(self,target,text_edit,browse_button,clear_button,application,title=""):
		self.target = target
		self.application = application
		self.text_edit = text_edit
		self.browser = None # this will be the browser itself
		self.browser_title = title
		QtCore.QObject.connect(browse_button,QtCore.SIGNAL("clicked(bool)"),self.browse_pressed)
		QtCore.QObject.connect(clear_button,QtCore.SIGNAL("clicked(bool)"),self.clear_pressed)
		
	def browse_pressed(self,bool):
		if self.browser == None:
			self.browser = EMSelectorModule(self.application)
			self.browser.widget.desktop_hint = "form" # this is to make things work as expected in the desktop
			self.browser.setWindowTitle(self.browser_title)
			self.application.show_specific(self.browser)
			QtCore.QObject.connect(self.browser.widget,QtCore.SIGNAL("ok"),self.on_browser_ok)
			QtCore.QObject.connect(self.browser.widget,QtCore.SIGNAL("cancel"),self.on_browser_cancel)
		else:
			self.application.show_specific(self.browser)

	def on_browser_cancel(self):
		self.application.close_specific(self.browser)
		self.browser = None
	def on_browser_ok(self,stringlist):
		new_string = str(self.text_edit.toPlainText())
		present_files = new_string.split() # this is to avoid redundancies
		for i,s in enumerate(stringlist):
			if s in present_files: continue
			if len(new_string) != 0 :
				new_string += '\n'
			new_string += s

		self.text_edit.setText(new_string)
		self.target.update_texture()# in the desktop the texture would have to be updated
		self.application.close_specific(self.browser)
		self.browser = None
	def clear_pressed(self,bool):
		#self.target.update_texture()# in the desktop the texture would have to be updated
		self.text_edit.clear()

def get_example_params():
	params = []
	params.append(ParamDef(name="box size",vartype="int",desc_short="int",desc_long="An integer value",property=None,defaultunits=128,choices=[]))
	params.append(ParamDef(name="apix",vartype="float",desc_short="float",desc_long="A floating point value",property=None,defaultunits=1.0,choices=[]))
	ps = ParamDef(name="model",vartype="string",desc_short="string",desc_long="A string value",property=None,defaultunits="Three thousand series",choices=None)
	params.append(ParamDef(name="Is there a God?",vartype="boolean",desc_short="boolean",desc_long="Something that is true or false",property=None,defaultunits=Util.get_irand(0,1),choices=None))
	p1 = ParamDef(name="True or false",vartype="boolean",desc_short="boolean",desc_long="Something that is true or false",property=None,defaultunits=Util.get_irand(0,1),choices=None)
	p2 = ParamDef(name="integer length",vartype="int",desc_short="int",desc_long="An integer value",property=None,defaultunits=128,choices=[])
	params.append([p1,p2,ps])
	params.append(ParamDef(name="Normalization method",vartype="string",desc_short="string",desc_long="Choose from a list of strings",property=None,defaultunits="normalize",choices=["normalize","normalize.edgemean","normalize.other"]))
	
	pi = ParamDef(name="Int 1 to 10",vartype="int",desc_short="int",desc_long="Choose from a list of ints",property=None,defaultunits=5,choices=[1,2,3,4,5,6,7,8,9,10])
	pf = ParamDef(name="Float 1 to 10",vartype="float",desc_short="float",desc_long="Choose from a list of floats",property=None,defaultunits=2.1,choices=[1.1,2.1,3.1,4.1,5.1,6.1,7.1,8.1,9.1,10.1])
	params.append([pi,pf])
	params.append(ParamDef(name="file_names",vartype="url",desc_short="url",desc_long="This is an editable list of file names with convenient browse and clear buttons",property=None,defaultunits=["tmp.txt","other.mrc"],choices=[]))
	params.append(ParamDef(name="comparitor",vartype="choice",desc_short="choice",desc_long="This is a string choice",property=None,defaultunits="frc",choices=["frc","phase","sqeuclidean"]))
	params.append(ParamDef(name="lucky number",vartype="choice",desc_short="choice",desc_long="This is to demonstrate that the type of choice is preserved. When you click ok the value in the return dictionary corresponding to this form entry will be an integer",property=None,defaultunits=3,choices=[1,2,3,98]))
	
	params.append(ParamDef(name="song",vartype="text",desc_short="text",desc_long="A potentially very long description",property=None,defaultunits="Jack and Jill went up the hill\nTo fetch a pail of water.\nJack fell down and broke his crown,\nAnd Jill came tumbling after.",choices=None))
	
	params.append(ParamDef(name="Selected Files",vartype="stringlist",desc_short="stringlist",desc_long="Choose from a list of strings",property=None,defaultunits=["C.mrc","E.mrc"],choices=[chr(i)+".mrc" for i in range(65,91)]))
	
	pil = ParamDef(name="Int 1 to 10 from a list",vartype="intlist",desc_short="intlist",desc_long="Choose from a list of ints",property=None,defaultunits=[5],choices=[1,2,3,4,5,6,7,8,9,10])
	pfl = ParamDef(name="Float 1 to 10 from a list",vartype="floatlist",desc_short="floatlist",desc_long="Choose from a list of floats",property=None,defaultunits=[2.1],choices=[1.1,2.1,3.1,4.1,5.1,6.1,7.1,8.1,9.1,10.1])
	a = ParamTable(name="table_choice",desc_short="Please choose using this information",desc_long="The left most column is what you're choosing from, the extra columns are used only to assist in the decision making process")
	a.append(pil)
	a.append(pfl)
	params.append([pil,pfl,a])
	
	return params

def on_ok(dict):
	print "got the ok signal, the return dictionary is",dict
	
def on_cancel():
	print "got the cancel signal"
	
# This is just for testing, of course
if __name__ == '__main__':
	
	from emapplication import EMStandAloneApplication
	em_app = EMStandAloneApplication()
	window = EMFormModule(params=get_example_params(),application=em_app)
	window.setWindowTitle("A test form")
	QtCore.QObject.connect(window.widget,QtCore.SIGNAL("emform_ok"),on_ok)
	QtCore.QObject.connect(window.widget,QtCore.SIGNAL("emform_cancel"),on_cancel)
	
	em_app.show()
	em_app.execute()
	
	
