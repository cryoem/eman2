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
from emapplication import EMQtWidgetModule,get_application
from EMAN2 import Util, get_image_directory	
import weakref

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
		self.enable_multiple_selection = True
		
		
class EMFormModule(EMQtWidgetModule):
	'''
	params should be a list of ParamDef objects
	application should be an EMAN2 type application
	After initializing this object you should probably call application.show(this)
	or application.show_specific(this), depending on what you're doing
	'''
	def __init__(self,params,application):
		self.application = weakref.ref(application)
		self.widget = EMFormWidget(self,params)
		EMQtWidgetModule.__init__(self,self.widget)
		
	def get_desktop_hint(self):
		return "form"
	
	#def __del__(self):
		#EMQtWidgetModule.__del__(self)

class EMFormWidget(QtGui.QWidget):
	'''
	See the example in __main__ below
	If ok is clicked the "emform_ok" signal is emitted along with a dictionary containing all of the form entries
	If cancel is clicked the "emform_cancel" signal is emmitted. No extra information is sent in this case
	'''
	def __init__(self,parent,params=None,disable_ok_cancel=False):
		QtGui.QWidget.__init__(self,None)
		self.parent = weakref.ref(parent)
		self.params = params
		self.event_handlers = [] # used to keep event handlers in memory
		self.resize_event_handlers = [] # used to keep resize event handlers in memory
		self.output_writers = [] # used to register output write objects, for the purpose of returning the results
		self.name_widget_map = {} # used to map parameter names to qt widgets - used by booleans to automatically disable and enable widgets
		self.__init_icons()
		self.setWindowIcon(QtGui.QIcon(get_image_directory() + "/eman.png"))
		
		self.auto_incorporate = {}
		self.auto_incorporate["float"] = self.__incorporate_float
		self.auto_incorporate["int"] = self.__incorporate_int
		self.auto_incorporate["url"] = self.__incorporate_url
		self.auto_incorporate["choice"] = self.__incorporate_choice
		self.auto_incorporate["string"] = self.__incorporate_string
		self.auto_incorporate["text"] = self.__incorporate_text
		self.auto_incorporate["boolean"] = self.__incorporate_boolean
		self.auto_incorporate["stringlist"] = self.__incorporate_stringlist
		self.auto_incorporate["intlist"] = self.__incorporate_intlist
		self.auto_incorporate["floatlist"] = self.__incorporate_floatlist
		self.auto_incorporate["dict"] = self.__incorporate_dict
		self.auto_incorporate["paramtable"] = self.__incorporate_paramtable
		
		self.vbl = QtGui.QVBoxLayout()
		self.incorporate_params(self.params,self.vbl)
		if not disable_ok_cancel: self.__add_ok_cancel_buttons(self.vbl)
		self.setLayout(self.vbl)
	def __init_icons(self):
		self.emdata_icon = QtGui.QIcon(get_image_directory() + "/single_image.png")
		self.emdata_3d_icon = QtGui.QIcon(get_image_directory() + "/single_image_3d.png")
		self.emdata_matrix_icon = QtGui.QIcon(get_image_directory() + "/multiple_images.png")
		self.plot_icon = QtGui.QIcon(get_image_directory() + "/plot.png")
		
	def incorporate_params(self,params,layout):
		for param in self.params:
			try:
				if len(param) != 1 and not isinstance(param,ParamTable):
					hbl=QtGui.QHBoxLayout()
					for iparam in param:
						self.auto_incorporate[iparam.vartype](iparam,hbl)
					layout.addLayout(hbl)
					continue
					
			except: pass
			self.auto_incorporate[param.vartype](param,layout)
		
		self.enable_boolean_dependents()
	def enable_boolean_dependents(self):
		''' This just takes care of making sure any dependent widgets are correctly enabled by default
		'''
		for event_handler in self.event_handlers:
			if isinstance(event_handler,BoolDependentsEventHandler):
				event_handler.checkbox_state_changed() # the argument is essentially unused
				
	def get_ptable_icon(self,paramtable):
		if hasattr(paramtable,"icon_type"):
			icon_type = paramtable.icon_type
			if icon_type == "single_image": return self.emdata_icon
			elif icon_type == "matrix_image": return self.emdata_matrix_icon
			elif icon_type == "3d_image": return self.emdata_3d_icon
			elif icon_type == "2d_plot" : return self.plot_icon
			else: return None
			
		return None
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
		icon = self.get_ptable_icon(paramtable)
		
		if not paramtable.enable_multiple_selection:
			table_widget.setSelectionMode(QtGui.QAbstractItemView.SingleSelection)
			
		selected_items = [] # used to ensure default selection is correct
		exclusions = []
		if hasattr(paramtable,"exclusions"): exclusions = paramtable.exclusions
		table_widget.setSortingEnabled(False)
		max_len_sum = 0
		for i,param in enumerate(paramtable):
			for j,choice in enumerate(param.choices):
				if i == 0 and icon != None: item = QtGui.QTableWidgetItem(icon,str(choice))
				else: item = QtGui.QTableWidgetItem(str(choice))
				flag2 = Qt.ItemFlags(Qt.ItemIsSelectable)
				flag3 = Qt.ItemFlags(Qt.ItemIsEnabled)
				flag4 = Qt.ItemFlags(Qt.ItemIsEditable)
				if i == 0:
					if str(choice) not in exclusions:
						item.setFlags(flag2|flag3)
					else:
						# exluded items are displayed but they are not selectable
						# this was originally added for e2boxer -the write output form needs to show which images are are excluded
						item.setFlags(flag3)
						item.setTextColor(QtGui.QColor(0,128,0))
						item.setToolTip("This item is excluded")
					if param.defaultunits != None and len(param.defaultunits) > 0:
						if choice in param.defaultunits:
							selected_items.append(item)
				else:
					item.setFlags(flag3)
				item.setTextAlignment(QtCore.Qt.AlignHCenter)
				
				table_widget.setItem(j, i, item)
				
			item = QtGui.QTableWidgetItem(param.desc_short)
			item.setTextAlignment(QtCore.Qt.AlignHCenter)
			item.setToolTip(param.desc_long)
			table_widget.setHorizontalHeaderItem(i,item)
		
		table_widget.setSortingEnabled(True)
		table_widget.setToolTip(paramtable.desc_long)
		hbl.addWidget(table_widget,1)
		
		groupbox = QtGui.QGroupBox(paramtable.desc_short)
		groupbox.setToolTip(paramtable.desc_long)
		groupbox.setLayout(hbl)
		layout.addWidget(groupbox,10)
		
		optional_attr = ["convert_text","context_menu"]
		for opt in optional_attr:
			if hasattr(paramtable,opt):
				setattr(table_widget,opt,getattr(paramtable,opt))
		
		table_event_handler = ParamTableEventHandler(self,table_widget)
		self.event_handlers.append(table_event_handler)
		self.resize_event_handlers.append(table_event_handler)
		table_widget.resizeColumnsToContents()
		
		if len(paramtable[0].choices) > 0:
			type_of = type(paramtable[0].choices[0])
		else:
			type_of = type(str) # this really doesn't matter
	
		for item in selected_items: 
			item.setSelected(True)
		self.output_writers.append(ParamTableWriter(paramtable.name,table_widget,type_of))
		self.name_widget_map[paramtable.name] = groupbox
		
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
		
		layout.addWidget(groupbox,0)

		self.output_writers.append(ListWidgetParamWriter(param.name,list_widget,type_of))
		self.name_widget_map[param.name] = groupbox
		#layout.addLayout(hbl)
	
	def __incorporate_boolean(self,param,layout):
		hbl=QtGui.QHBoxLayout()
		hbl.setMargin(0)
		hbl.setSpacing(2)
		check_box = QtGui.QCheckBox(str(param.desc_short),self)
		check_box.setChecked(bool(param.defaultunits))
		check_box.setToolTip(param.desc_long)
		hbl.addWidget(check_box,0)
		layout.addLayout(hbl)
		self.output_writers.append(BoolParamWriter(param.name,check_box))
		self.name_widget_map[param.name] = check_box
		
		if hasattr(param,"dependents"):
			self.event_handlers.append(BoolDependentsEventHandler(self,check_box,param.dependents))    
	
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
			hbl.addWidget(line_edit,0)
			hbl.name = param.name
			layout.addLayout(hbl)
			self.output_writers.append(StringParamWriter(param.name,line_edit))
			self.name_widget_map[param.name] =  [line_edit,label]
		
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
			hbl.addWidget(line_edit,0)
			hbl.name = param.name
			layout.addLayout(hbl)
			self.output_writers.append(FloatParamWriter(param.name,line_edit))
			self.name_widget_map[param.name] = [line_edit,label]
		
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
			hbl.addWidget(line_edit,0)
			hbl.name = param.name
			layout.addLayout(hbl)
			self.output_writers.append(IntParamWriter(param.name,line_edit))
			self.name_widget_map[param.name] = [line_edit,label]
		
	def __incorporate_text(self,param,layout):
		hbl=QtGui.QHBoxLayout()
		hbl.setMargin(0)
		hbl.setSpacing(2)
				
		text_edit = QtGui.QTextEdit("",self)
		text_edit.setReadOnly(True)
		text_edit.setWordWrapMode(QtGui.QTextOption.WordWrap)
		text_edit.setText(param.defaultunits)
		hbl.addWidget(text_edit,0)
		
		groupbox = QtGui.QGroupBox(param.desc_short)
		groupbox.setToolTip(param.desc_long)
		groupbox.setLayout(hbl)
		
		layout.addWidget(groupbox,1)

		self.output_writers.append(TextParamWriter(param.name,text_edit))	
		self.name_widget_map[param.name] = groupbox
	
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
		hbl.addWidget(text_edit,0)
		vbl.addLayout(hbl)
		
		hbl2=QtGui.QHBoxLayout()
		hbl2.setMargin(0)
		hbl2.setSpacing(2)
		
		browse_button = QtGui.QPushButton("Browse",self)
		hbl2.addWidget(browse_button)
		clear_button = QtGui.QPushButton("Clear",self)
		hbl2.addWidget(clear_button,0)
		vbl.addLayout(hbl2)
		
		groupbox = QtGui.QGroupBox(param.desc_short)
		groupbox.setToolTip(param.desc_long)
		groupbox.setLayout(vbl)
		
		layout.addWidget(groupbox,0)
		
		self.event_handlers.append(UrlEventHandler(self,text_edit,browse_button,clear_button,self.parent().application(),param.desc_short))
		self.output_writers.append(UrlParamWriter(param.name,text_edit))
		self.name_widget_map[param.name] = groupbox

	def __incorporate_dict(self,param,layout):
		'''
		A dictionary is turned into two combo boxes - An event handler recognises when the first combo value changes (the dictionary keys), and if so
		changes the entries in the second combo to values in the dictionary corresponding to the keys
		'''
		# hbl - tht
		hbl=QtGui.QHBoxLayout()
		hbl.setMargin(0)
		hbl.setSpacing(2)
		
		keys = param.choices.keys()
		keys.sort() # yes this is somewhat restrictive but it was my only way around something
#		label = QtGui.QLabel(param.desc_short+":",self)
#		label.setToolTip(param.desc_long)
#		hbl.addWidget(label)
		
		combo = QtGui.QComboBox(self)
		idx_default = 0
		for i,k in enumerate(keys):
			combo.addItem(str(k))
			
			if str(k) == str(param.defaultunits):
				idx_default = i
				
		combo.setCurrentIndex(idx_default)
		combo_default = keys[idx_default]
		hbl.addWidget(combo)
		
		
		combo2 = QtGui.QComboBox(self)
		for v in param.choices[combo_default]:
			combo2.addItem(str(v))
			
		hbl.addWidget(combo2)
		
		groupbox = QtGui.QGroupBox(param.desc_short)
		groupbox.setToolTip(param.desc_long)
		groupbox.setLayout(hbl)
		
		layout.addWidget(groupbox)
		
		self.output_writers.append(DictParamWriter(param,combo,combo2))
		self.event_handlers.append(DictEventHandler(param.choices,combo,combo2))
		
		self.name_widget_map[param.name] = groupbox

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
		layout.addWidget(groupbox,0)
		self.output_writers.append(ChoiceParamWriter(param.name,buttons,type(param.choices[0])))
		
		self.name_widget_map[param.name] = groupbox
	
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

		hbl.addWidget(combo,0)
		layout.addLayout(hbl)
		
		self.output_writers.append(FloatChoiceParamWriter(param.name,combo))
		self.name_widget_map[param.name] = [combo,label]
	
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

		hbl.addWidget(combo,0)
		layout.addLayout(hbl)
		
		self.output_writers.append(IntChoiceParamWriter(param.name,combo))
		self.name_widget_map[param.name] = [combo,label]
	
	def __incorporate_string_with_choices(self,param,layout):
		hbl=QtGui.QHBoxLayout()
		hbl.setMargin(0)
		hbl.setSpacing(2)
		
		hbl=QtGui.QHBoxLayout()
		label = QtGui.QLabel(param.desc_short,self)
		label.setToolTip(param.desc_long)
		hbl.addWidget(label,0)
				
		combo = QtGui.QComboBox(self)
		idx_default = 0
		for i,string in enumerate(param.choices):
			combo.addItem(string)
			if string == param.defaultunits:
				idx_default = i
		combo.setCurrentIndex(idx_default)

		hbl.addWidget(combo,0)
		layout.addLayout(hbl)
		
		self.output_writers.append(StringChoiceParamWriter(param.name,combo))
		self.name_widget_map[param.name] = [combo,label]
	
	def __add_ok_cancel_buttons(self,layout):
		hbl=QtGui.QHBoxLayout()
		label = QtGui.QLabel("Form commands:",self)
		hbl.addWidget(label)
		
		ok_button = QtGui.QPushButton("Ok",self)
		ok_button.setToolTip("When you click ok the values in the form are sent to the calling program")
		hbl.addWidget(ok_button)
		cancel_button = QtGui.QPushButton("Cancel",self)
		hbl.addWidget(cancel_button,0)
		layout.addLayout(hbl)
		QtCore.QObject.connect(ok_button,QtCore.SIGNAL("clicked(bool)"),self.ok_pressed)
		QtCore.QObject.connect(cancel_button,QtCore.SIGNAL("clicked(bool)"),self.cancel_pressed)
		
	def ok_pressed(self,bool):
		ret = {}
		for output in self.output_writers: output.write_data(ret)
		self.parent().emit(QtCore.SIGNAL("emform_ok"),ret) # getting the parent to emit ensures integration with the desktop
		
	def cancel_pressed(self,bool):
		self.parent().emit(QtCore.SIGNAL("emform_cancel")) # getting the parent to emit ensures integration with the desktop


	def update_texture(self):
		self.parent().force_texture_update()
		
	def closeEvent(self,event):
		self.parent().emit(QtCore.SIGNAL("emform_close"))
		
	def resizeEvent(self,event):
		for event_handler in self.resize_event_handlers:
			event_handler.resizeEvent(event)

	def display_file(self,filename):
		self.parent().emit(QtCore.SIGNAL("display_file"),filename)

		
class ParamTableWriter:
	def __init__(self,param_name,table_widget,type_of):
		self.param_name = param_name
		self.table_widget = table_widget
		self.type_of = type_of
		
	def write_data(self,dict):
		sel = [self.type_of(item.text()) for item in self.table_widget.selectedItems()]
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
		text = self.combo.currentText()
		if len(text) != 0:
			dict[self.param_name] = float(text)
		# else the key is not written to the dictionary!

class IntChoiceParamWriter:
	def __init__(self,param_name,combo):
		self.param_name = param_name
		self.combo = combo
		
	def write_data(self,dict):
		text = self.combo.currentText()
		if len(text) != 0:
			dict[self.param_name] = int(text)
		# else the key is not written to the dictionary!

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
		text = self.line_edit.text()
		if len(text) != 0:
			dict[self.param_name] = float(text)
		# else the key is not written to the dictionary!

class IntParamWriter:
	def __init__(self,param_name,line_edit):
		self.param_name = param_name
		self.line_edit = line_edit
		
	def write_data(self,dict):
		text = self.line_edit.text()
		if len(text) != 0:
			dict[self.param_name] = int(text)
		# else the key is not written to the dictionary!
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
		strings = [i for i in str(self.text_edit.toPlainText()).split('\n') if len(i) != 0]
		rm = []
		for i,s in enumerate(strings):
			if len(s) == 0:
				rm.append(i)

				
		rm.reverse()
		for i in rm:
			strings.pop(i)
		dict[self.param_name] = strings

class DictParamWriter:
	def __init__(self,param,combo1,combo2):
		self.param = param
		self.combo1 = combo1
		self.combo2 = combo2
		
	def write_data(self,dict):
		'''
		Here I do my best to preserve the type - the combo converted the keys and values into strings. Now I try to convert them back
		Also there are two keys entered into the dictionary, these keys are extracted based on the param.name parameter, which is split - if the returning list
		has more than 1 entry than the first and the last are the keys, respectively. Else if the list has only one entry, it becomes the first key, and the string "_selection" is appended
		to this to make the second key
		
		'''
		return_keys = self.param.name.split()
		if len(return_keys) == 0: raise # there is no generic mechanism - in addition, what if there is more than one dictionary in the form?
		elif len(return_keys) == 1:
			key1 = return_keys[0]
			key2 = key1+"_selection"
		else:
			key1 = return_keys[0]
			key2 = return_keys[-1]
	
		# get value1
		idx1 = self.combo1.currentIndex()
		keys = self.param.choices.keys()
		keys.sort() # because it was sorted above 
		value1 = keys[idx1] # this preserves the type - overkill, well, foolproof, yes a bit more so
		
		# get value2
		idx2 = self.combo2.currentIndex()
		value2 = self.param.choices[value1][idx2] # type preserving again
		
		dict[key1] = value1
		dict[key2] = value2
		
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

class ParamTableEventHandler:
	'''
	handles events for param tables, atm this is only the double click event, which can
	be used to trigger image display, for example
	'''
	def __init__(self,target,table_widget):
		self.target = weakref.ref(target)
		self.table_widget = table_widget
		table_widget.contextMenuEvent = self.contextMenuEvent
				
		QtCore.QObject.connect(table_widget, QtCore.SIGNAL("itemDoubleClicked(QTableWidgetItem*)"),self.table_item_double_clicked)
		
	def table_item_double_clicked(self,item):
		if hasattr(self.table_widget,"convert_text"):
			self.target().display_file( self.table_widget.convert_text(str(item.text())))
		else: pass
	
	def contextMenuEvent(self,event):
		if hasattr(self.table_widget,"context_menu"):
			menu = QtGui.QMenu()
			for k in self.table_widget.context_menu.keys():
				menu.addAction(k)
			QtCore.QObject.connect(menu,QtCore.SIGNAL("triggered(QAction*)"),self.menu_action_triggered)
			menu.exec_(event.globalPos())
	
	def menu_action_triggered(self,action):
		items = self.table_widget.selectedItems()
		names = [self.table_widget.convert_text(str(item.text())) for item in items]
		self.table_widget.context_menu[str(action.text())](names,self.target().parent().application())
	
	def table_item_clicked(self,item):
		#
		pass
	def resizeEvent(self,event):
		return
		cols = self.table_widget.columnCount()
		cumulative_width = 0
		for i in range(cols):
			cumulative_width += self.table_widget.columnWidth(i)
		
		tab_widget_width = self.table_widget.geometry().width()
		print cumulative_width,self.table_widget.width(),tab_widget_width,self.table_widget.frameSize().width()
		print self.table_widget.getContentsMargins()
		
		if cumulative_width < self.table_widget.width():
			scale = float(self.table_widget.width())/cumulative_width
			for i in range(cols):
				w = self.table_widget.columnWidth(i)
				self.table_widget.setColumnWidth(i,w*scale)
		
		
		
class UrlEventHandler:
	'''
	The browse and cancel events have to be sent to the correct line edit, so this handles it
	'''
	def __init__(self,target,text_edit,browse_button,clear_button,application,title=""):
		self.target = weakref.ref(target)
		self.application = weakref.ref(application)
		self.text_edit = text_edit
		self.browser = None # this will be the browser itself
		self.browser_title = title
		QtCore.QObject.connect(browse_button,QtCore.SIGNAL("clicked(bool)"),self.browse_pressed)
		QtCore.QObject.connect(clear_button,QtCore.SIGNAL("clicked(bool)"),self.clear_pressed)
		
	def browse_pressed(self,bool):
		if self.browser == None:
			self.browser = EMSelectorModule(False,False)
			self.browser.widget.desktop_hint = "form" # this is to make things work as expected in the desktop
			self.browser.setWindowTitle(self.browser_title)
			get_application().show_specific(self.browser)
			QtCore.QObject.connect(self.browser,QtCore.SIGNAL("ok"),self.on_browser_ok)
			QtCore.QObject.connect(self.browser,QtCore.SIGNAL("cancel"),self.on_browser_cancel)
		else:
			get_application().show_specific(self.browser)

	def on_browser_cancel(self):
		get_application().close_specific(self.browser)
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
		self.target().update_texture()# in the desktop the texture would have to be updated
		get_application().close_specific(self.browser)
		self.browser = None
	def clear_pressed(self,bool):
		#self.target().update_texture()# in the desktop the texture would have to be updated
		self.text_edit.clear()
		
class DictEventHandler:
	'''
	Dictionaries are presented as two combo boxes - when the first combo box changes the values in the second box are updated (according to what is in the dictionary)
	'''
	def __init__(self,dict,combo1,combo2):
		self.dict = dict
		self.combo1 = combo1
		self.combo2 = combo2
		
		QtCore.QObject.connect(self.combo1, QtCore.SIGNAL("currentIndexChanged(int)"),self.combo1_index_changed)
	
	def combo1_index_changed(self,i):
		
		keys = self.dict.keys()
		keys.sort() # because the keys are sorted in the display
		key = keys[i]
		values = self.dict[key]
		
		self.combo2.clear()
		for v in values:
			self.combo2.addItem(str(v))
			
class BoolDependentsEventHandler:
	'''
	This event handler works on the assumption that a boolean type is always a checkbox.
	If the boolean types also has the "dependents" attribute, then toggling the checkbox will enable/disable the dependent widgets
	'''
	def __init__(self,target,checkbox,dependents):
		'''
		Need the target because we need access to the name_widget_map
		'''
		self.target = weakref.ref(target)
		self.checkbox = checkbox
		self.dependents = dependents # a list of depent names (ParamDef.name)
		
		QtCore.QObject.connect(self.checkbox, QtCore.SIGNAL("stateChanged(int)"),self.checkbox_state_changed)
		
	def checkbox_state_changed(self,integer=0):
		name_map = self.target().name_widget_map
		
		enabled = False
		if self.checkbox.isChecked(): enabled = True
		for name in self.dependents:
			widgets = name_map[name]
			if isinstance(widgets,list):
				for widget in widgets: widget.setEnabled(enabled)
			else: widgets.setEnabled(enabled)
		

class EMTableFormModule(EMQtWidgetModule):
	def __init__(self,params,application):
		self.application = weakref.ref(application)
		self.widget = EMTableFormWidget(self,params)
		EMQtWidgetModule.__init__(self,self.widget)
		
	def get_desktop_hint(self):
		return "form"
		
class EMTableFormWidget(EMFormWidget):
	'''
	See the example in __main__ below
	If ok is clicked the "emform_ok" signal is emitted along with a dictionary containing all of the form entries
	If cancel is clicked the "emform_cancel" signal is emmitted. No extra information is sent in this case
	'''
	def __init__(self,parent,params=None):
		EMFormWidget.__init__(self,parent,params)
		
		
	def incorporate_params(self,params,layout):
		tabwidget = QtGui.QTabWidget(self)
		
		for title,paramlist in params:
			
			widget = QtGui.QWidget(None)
			vbl =  QtGui.QVBoxLayout(widget)
			for param in paramlist:
				
				try:
					if len(param) != 1 and not isinstance(param,ParamTable):
						hbl=QtGui.QHBoxLayout()
						for iparam in param:
							self.auto_incorporate[iparam.vartype](iparam,hbl)
						vbl.addLayout(hbl)
						continue
						
				except: pass
				self.auto_incorporate[param.vartype](param,vbl)
			
			tabwidget.addTab(widget,title)
		
		layout.addWidget(tabwidget)
		
		self.enable_boolean_dependents()
				

def get_example_form_params():
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
	params.append(ParamDef(name="comparator",vartype="choice",desc_short="choice",desc_long="This is a string choice",property=None,defaultunits="frc",choices=["frc","phase","sqeuclidean"]))
	params.append(ParamDef(name="lucky number",vartype="choice",desc_short="choice",desc_long="This is to demonstrate that the type of choice is preserved. When you click ok the value in the return dictionary corresponding to this form entry will be an integer",property=None,defaultunits=3,choices=[1,2,3,98]))
	
	params.append(ParamDef(name="song",vartype="text",desc_short="text",desc_long="A potentially very long description",property=None,defaultunits="Jack and Jill went up the hill\nTo fetch a pail of water.\nJack fell down and broke his crown,\nAnd Jill came tumbling after.",choices=None))
	
	params.append(ParamDef(name="Selected Files",vartype="stringlist",desc_short="stringlist",desc_long="Choose from a list of strings",property=None,defaultunits=["C.mrc","E.mrc"],choices=[chr(i)+".mrc" for i in range(65,91)]))
	
	params.append(ParamDef(name="Combo1 and combo2", vartype="dict",desc_short="dict",desc_long="Choose from the combo on the left and the combo box on the right will be updated. The 'and' in the the name differentiates the key in the return dictionary", property=None, defaultunits="reconstructors", choices={"reconstructors":["a","b","c"],"processors":["1","2","3"]} ))
	
	pil = ParamDef(name="Int 1 to 10 from a list",vartype="intlist",desc_short="intlist",desc_long="Choose from a list of ints",property=None,defaultunits=[5],choices=[1,2,3,4,5,6,7,8,9,10])
	pfl = ParamDef(name="Float 1 to 10 from a list",vartype="floatlist",desc_short="floatlist",desc_long="Choose from a list of floats",property=None,defaultunits=[2.1],choices=[1.1,2.1,3.1,4.1,5.1,6.1,7.1,8.1,9.1,10.11111111111111111111111111111])
	a = ParamTable(name="table_choice",desc_short="Please choose using this information",desc_long="The left most column is what you're choosing from, the extra columns are used only to assist in the decision making process")
	a.append(pil)
	a.append(pfl)
	params.append([pil,pfl,a])
	
	return params
def get_example_table_form_params():
	params = get_example_form_params()
	p1 = params[0:len(params)/3]
	p2 = params[len(params)/3:2*len(params)/3]
	p3 = params[2*len(params)/3:]
	
	par =[]
	par.append(["Input", p1])
	par.append(["This too",p2])
	par.append(["And also", p3])
	return par

def on_ok(dict):
	print "got the ok signal, the return dictionary is",dict
	
def on_cancel():
	print "got the cancel signal"
	
# This is just for testing, of course
if __name__ == '__main__':
	
	from emapplication import EMStandAloneApplication
	em_app = EMStandAloneApplication()
	window = EMFormModule(params=get_example_form_params(),application=em_app)
	window.setWindowTitle("A test form")
	QtCore.QObject.connect(window,QtCore.SIGNAL("emform_ok"),on_ok)
	QtCore.QObject.connect(window,QtCore.SIGNAL("emform_cancel"),on_cancel)
	
	window2= EMTableFormModule(params=get_example_table_form_params(),application=em_app)
	window2.setWindowTitle("A test form")
	QtCore.QObject.connect(window2,QtCore.SIGNAL("emform_ok"),on_ok)
	QtCore.QObject.connect(window2,QtCore.SIGNAL("emform_cancel"),on_cancel)
	
	em_app.show()
	em_app.execute()
	
	
