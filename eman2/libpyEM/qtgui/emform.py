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
from emselector import EMSelectorDialog
from emapplication import get_application
from EMAN2 import Util, get_image_directory,file_exists,dump_aligners_list,dump_processors_list
import EMAN2
import weakref
import warnings
import time

class EMButtonDialog:
	'''
	A base class for adding a dialog to the form
	Call add_to_layout to add the button to some layout.
	'''
	def __init__(self,desc_short="Dialog Button",desc_long="",icon=None):
		'''
		@param desc_short that which will appear as the main text on the button
		@param desc_long the tool tip for the button
		@param icon if specificed this is a QtGui.QIcon, and will become the icon of the button
		'''
		self.target = None
		self.vartype = "EMButtonDialog"
		self.desc_short = desc_short
		self.desc_long = desc_long
		self.icon = icon
	def set_target(self,target): 
		'''
		@param target - an instance of an EMFormWidget
		Makes a weak reference to the target, seeing as the target has a reference to this
		'''
		self.target = weakref.ref(target)
		
	def add_to_layout(self,layout):
		'''
		Add the button to the given layout
		@param layout a QtGui.QLayout object
		'''
		if self.icon != None: self.button = QtGui.QPushButton(self.icon,self.desc_short)
		else: self.button = QtGui.QPushButton(self.desc_short)
		self.button.setToolTip(self.desc_long)
		layout.addWidget(self.button)
		QtCore.QObject.connect(self.button, QtCore.SIGNAL("clicked(bool)"), self.on_button)
	
	def on_button(self,unused=None): 
		'''
		Inheriting classes must supply this function
		The slot that is connected to the signal emitted when the button is clicked
		'''
		raise NotImplementedError("Inheriting classes are menat to supply this function")
		
class EMOrientationDistDialog(EMButtonDialog):
	'''
	This class is used by the refine form in the workflow - it runs a dialog enabling the user to select the 
	parameters that will define the distribution of projections on the unit sphere
	It is highly specialized - it assumes that the EMSymChoiceDialog returns a dictionary with certain keys in it,
	and assumes the that self.target (which is a weak ref to and EMFormWidget) has an attribute called name_widget_map
	which is a dictionary - that dictionary is assumed to have certain keys, and the values are assumed to lists which have
	widgets in an assumed order, or single objects
	In short, liable to break if someone makes changes to the code.
	'''
	def __init__(self):
		EMButtonDialog.__init__(self,desc_short="Interactive Parameters",desc_long="Interact with the unit sphere to determine your orientation distribution parameters",icon=QtGui.QIcon(get_image_directory() + "eulerxplor.png"))
		
	def on_button(self,unused=None):
		'''
		The crux of this class - dictionary keys and various attributes are assumed to exist
		'''
		name_map = self.target().name_widget_map
		symname = str(name_map["symname"][0].currentText()) # assume the first entry in the list is the combo
		symnumber = name_map["symnumber"][0] # assume the first entry in the list is a text entry
		if symnumber.isEnabled(): symname += str(symnumber.text())
		
		from emimage3dsym import EMSymChoiceDialog
		dialog = EMSymChoiceDialog(symname)
		result = dialog.exec_()
		if result != None:
			combo = name_map["orientgen"][0] # assume the first entry in the list is the combo
			s = [combo.itemText(i) for i in range(combo.count())]
			for i,text in enumerate(s):
				if text == result["orientgen"]:
					combo.setCurrentIndex(i)
					break
			else:
				raise RuntimeError("The return value from the dialog was not present in the orientgen combo. This means dependent code has evolved. Contact developers")
			
			checkbox = name_map["incmirror"]
			checkbox.setChecked(result["inc_mirror"])
			
			approach = result["approach"]
			button_group = name_map["orientopt"][1] # assume the second entry in the list is a list of buttons
			for button in button_group:
				if button.text() == approach:
					button.setChecked(True)
					break
			else:
				raise RuntimeError("The approach key was unknown or has changed. Contact developers")
			
			value = result["value"]
			float_widget = name_map["orientopt_entry"][0] # assume the first entry in the list is a text entry
			float_widget.setText(str(value))
			
			sym = result["sym"]
			if sym != symname:
				sym_n = None
				if sym not in ["icos","oct","tet"]:
					if len(sym) > 1:
						sym_n = int(sym[1:])
					sym = sym[:1]
					
				sym_combo = name_map["symname"][0] # assume the first entry in the list is the combo
				s = [sym_combo.itemText(i) for i in range(sym_combo.count())]
				for i,text in enumerate(s):
					if text == sym:
						sym_combo.setCurrentIndex(i)
						break
				else:
					raise RuntimeError("The symmetry value was uninterpretable. Contact developers.")
				
				if sym_n != None:
					symnumber.setText(str(sym_n))
	
class EMParamTable(list):
	'''
	NOTE: THE EMPARAMTABLE BECAME DEPRECATED IN FAVOR OF THE EMFILETABLE around May 2009. However, this class is still used in remote cases that need to be cleaned up
	This object is constructed as a table in the EMFormWidget
	This object list of ParamDef objects : each of these ParamDef objects are generally a list-type of some sort, such
	as intlist, floatlist, stringlist etc. They should all be the same length, but this is not tested - to do so would 
	require the list.append function to be redefined
	Presently, the user selects their options only from the first column. The other columns are useful only for displaying helpful metadata, such as image mean and image dimensions
	'''
	def __init__(self,name=None,desc_short=None,desc_long=""):
		'''
		@param name the all important name of the parameter supplied by this object, becomes a key in a dictionary.
		@param desc_short a concise (short) descriptive title. Will become a general title for the associated table in EMFormWidget
		@param desc_long will be used as a tool tip. Text that is helpful but not too long
		'''
		list.__init__(self)
		self.vartype = "EMParamTable"
		self.name = name #
		self.desc_short = desc_short # Will become a general title for the associated table
		self.desc_long = desc_long # For the tool tip
		self.enable_multiple_selection = True # use by EMFormWidget when it creates the associated table - the user can select more than one entry if this is true
		
	def custom_addition(self,layout,table_widget):
		'''
		When the associated table is being created and added in the form module, this function is called, enabling
		different things to be added to the layout in a custom fashion (such as an "Add" button).
		@param layout a Qt Layout (e.g. QVBoxLayout, QHBoxLayout - objects that support the 'addWidget' and 'addLayout' syntax
		@param table_widget - the table widget itself, which is an instance of a QtGui.QTableWidget
		'''
		pass
	
	def add_optional_table_attr(self,table_widget):
		'''
		Called by the EMFormWidget to add extra attributes to the the table widget, and this is used
		as the basis for creating context menus (context_menu attribute, which is a dictionary), and for
		converting the name in the table to the absolute file system path (convert_text attribute,
		which is a function)
		@param table_widget the QtGui.QTableWidget which will have the new attributes
		'''
		optional_attr = ["convert_text","context_menu"] # these two attributes are the only ones currently used (even for inheriting classes)
		for opt in optional_attr:
			if hasattr(self,opt):
				setattr(table_widget,opt,getattr(self,opt))
				
	def convert_text(file_name): return file_name
	
	convert_text = staticmethod(convert_text)
				
	def build_table(self,table_widget,icon):
		
		exclusions = []
		if hasattr(self,"exclusions"): exclusions = self.exclusions # exclusions are a list of highlight entries - they get colored green
		
		flag2 = Qt.ItemFlags(Qt.ItemIsSelectable)
		flag3 = Qt.ItemFlags(Qt.ItemIsEnabled)
		flag4 = Qt.ItemFlags(Qt.ItemIsEditable)
		selected_items = []
		for i,param in enumerate(self):
			for j,choice in enumerate(param.choices):
				if i == 0 and icon != None: item = QtGui.QTableWidgetItem(icon,str(choice))
				else: item = QtGui.QTableWidgetItem(str(choice))
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
			
		for item in selected_items: 
			item.setSelected(True)

class EMFileTable(QtGui.QTableWidget):
	def __init__(self,listed_names=[],name="filenames",desc_short="File Names",desc_long="A list of file names",single_selection=False,enable_save=True):
		'''
		@param listed_names The names that will be listed in the first column of the table
		@param column_data A list of EMFileTable.EMContextMenuData objects
		'''
		
		QtGui.QTableWidget.__init__(self)
		self.name = name # the name of the parameter ultimately return form the form
		self.listed_names = listed_names # listed names in first column
		self.column_data = [] # list of EMColumnData objects
		self.column_data_refs = [] # keep a reference to more complicated to column data objects - someone has to
		self.button_data = [] # extra button info which can be used to add things to the table
		self.icon = None # set be inheriting function
		self.exclusions = [] # Correspond to listed names that will be colored green and disabled
		self.default_selections = [] # Default selected names in first column
		self.desc_short = desc_short
		self.desc_long = desc_long
		self.vartype = "file_table" # This is used by the EMFormWidget to insert this object correctly into a widget
		self.name_conversions = {} # This is used to convert the displayed name to the real name of the file on the operating system
		self.context_menu_data = {} # see self.get_context_menu_dict help
		QtCore.QObject.connect(self, QtCore.SIGNAL("itemDoubleClicked(QTableWidgetItem*)"),self.table_item_double_clicked)

		if enable_save: self.context_menu_data["Save As"] = EMFileTable.save_as
		self.context_menu_refs = [] # to keep a reference to context menus related objects - somebody has to

		self.single_selection = single_selection
		if self.single_selection:
			self.setSelectionMode(QtGui.QAbstractItemView.SingleSelection)
		else:
			self.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
		
		self.timer = None
		self.timer_interval = 5000
		self.animated_columns = {}
		self.busy = 0
		self.items_selected_by_default = True
		
	def set_items_selected_by_default(self,val):
		self.items_selected_by_default = val
		
	def register_animated_column(self,column_title):
		if self.timer == None:
			self.timer = QtCore.QTimer()
			QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.time_out)
			self.timer.start(self.timer_interval)
			
		self.animated_columns[column_title] = -1 # -1 is a flag
		
	def time_out(self):
		if self.busy :
			self.timer.stop()
			from emapplication import EMErrorMessageDisplay
			EMErrorMessageDisplay.run(["Disabling updates of %s for speed" %key] )
		
		stime=time.time()	
		self.busy = 1
		for key,value in self.animated_columns.items():
			if value == -1:
				for i in xrange(0,self.columnCount()):
					if (str(self.horizontalHeaderItem(i).text())) == key:
						self.animated_columns[key] = i
						break
				else:
					from emapplication import EMErrorMessageDisplay
					EMErrorMessageDisplay.run(["Can't animate %s" %key] )
					self.animated_columns.pop(key)
					self.busy = 0
					return
				
		
		for key,value in self.animated_columns.items():
			cd = self.column_data[value-1]
			for i in xrange(0,len(self.listed_names)):
				
				item = self.item(i,value)
				item.setText(cd.function(self.convert_text(str(self.item(i,0).text()))))

		self.busy = 0
		
		if time.time()-stime>0.5 :
			self.timer.stop()
			from emapplication import EMErrorMessageDisplay
			EMErrorMessageDisplay.run(["Disabling updates of %s for speed" %key] )
		
	
	def convert_text(self,name):
		'''
		Sometimes the first column displays a shortened version of the name of a file on
		disk, but it occurs that you want to recall the full name. This function does that
		'''
		for key,value in self.name_conversions.items():
			if value == name:
				return key
		return None
	
	def get_context_menu_dict(self):
		'''
		@return a dictionary - keys are used to add context menu items, values are functions which are called
		These functions (which are the dictionary values) take two arguments, the first being a list of strings
		the second being a QtGui.QTableWidget
		'''
		return self.context_menu_data
	
	def add_context_menu_action(self,name,action):
		'''
		For these purposes, a context menu action consists of a name and an action
		The name is a string, the action is a function - this function takes two arguments,
		one being a list of strings, the other being a table_widget. Note there is not self argument
		for these "action-functions"
		'''
		self.context_menu_data[name] = action
	
	def add_context_menu_data(self,context_menu_data):
		'''
		For these purposes, a context menu action consists of a name and an action
		The name is a string, the action is a function - this function takes two arguments,
		one being a list of strings, the other being a table_widget. Note there is not self argument
		for these "action-functions"
		'''
		self.context_menu_refs.append(context_menu_data)
		for key,value in context_menu_data.items():
			self.add_context_menu_action(key,value)

	def add_column_data(self,column_data):
		'''
		@param column_data an instance of (or something that has the attributes of) EMColumnData
		Only works if you call this before the table is constructed - doesn't "update" the table  
		'''
		if isinstance(column_data,EMFileTable.EMColumnData):
			self.column_data.append(column_data)
		else:
			self.column_data_refs.append(column_data)
			for key,value in column_data.column_data.items():
				self.column_data.append(EMFileTable.EMColumnData(key,value,""))
	def remove_column_data(self,column_data_name):
		'''
		@param column_data_name the name attribute of the column data attribute that you wish to remove
		@exception RuntimeError raised if there is not column_data entry that has that specific name
		Will not work if there is more than one column data entry with the same name
		'''
		for i,column_data in enumerate(self.column_data):
			if column_data.name == column_data_name:
				self.column_data.pop(i)
				break
		else:
			raise RuntimeError("Attempt to remove a column data that didn't exist (%s)" %column_data_name)
			
		
	def insert_column_data(self,idx,column_data):
		'''
		@param column_data an instance of (or something that has the attributes of) EMColumnData
		Only works if you call this before the table is constructed - doesn't "update" the table  
		'''
		self.column_data.insert(idx,column_data)
		
	def add_button_data(self,button_data):
		'''
		Buttons are added to the layout containing the table. They are below the table
		@param button_data an instance of (or something that has the attributes of) EMButton data
		Only works if you call this before the table is constructed - doesn't "update" the table  
		'''
		self.button_data.append(button_data)
	
	def convert_name(self,name):
		'''
		This is used to get the display name of the file. For instance, to display the base name not the entire name
		@param name a file name - The file should exist on disk
		@return the desired converted name
		'''
		return EMAN2.base_name(name,bdb_keep_dir=True)
	
	def display_name(self,name):
		'''
		Get the display name for this file name
		@param name a file name - The file should exist on disk
		@return the display nane
		Internally caches the originally name in a dictionary so it can be recovered
		Redefine self.convert_name to achieve your custom-desired display name
		'''
		if not self.name_conversions.has_key(name):
			converted = self.convert_name(name) 
			self.name_conversions[name] = converted
			return converted # for efficiency
		
		return self.name_conversions[name]
	
	def build_table(self):
		'''
		Builds the table contents. Should only need to be called once, just prior to this object's
		insertion into a widget or layout
		'''
		sorting = self.isSortingEnabled()
		self.setSortingEnabled(False)
		self.setRowCount(len(self.listed_names))
		self.setColumnCount(len(self.column_data)+1)
			
		flag2 = Qt.ItemFlags(Qt.ItemIsSelectable)
		flag3 = Qt.ItemFlags(Qt.ItemIsEnabled)
		flag4 = Qt.ItemFlags(Qt.ItemIsEditable)
		
		selected_items = []

		# first step is to insert the first column, which is the file names
		for i,name in enumerate(self.listed_names):
			if self.icon != None: item = QtGui.QTableWidgetItem(self.icon,self.display_name(name))
			else: item = QtGui.QTableWidgetItem(self.display_name(name))

			if not file_exists(name):item.setTextColor(QtGui.QColor(0,128,0))
				
			if name not in self.exclusions :
				item.setFlags(flag2|flag3)
				item.setToolTip(name)
			else:
				# exluded items are displayed but they are not selectable
				# this was originally added for e2boxer -the write output form needs to show which images are are excluded
				item.setFlags(flag3)
				item.setTextColor(QtGui.QColor(0,128,0))
				item.setToolTip("This item is excluded")
			
			if name in self.default_selections:
				selected_items.append(item)
	
			item.setTextAlignment(QtCore.Qt.AlignHCenter)
				
			self.setItem(i, 0, item)
				
		item = QtGui.QTableWidgetItem(self.desc_short)
		item.setTextAlignment(QtCore.Qt.AlignHCenter)
		item.setToolTip(self.desc_long)
		self.setHorizontalHeaderItem(0,item)
		
		# second step is to add the columns
		col = 1
		for cd in self.column_data:
			item = QtGui.QTableWidgetItem(cd.name)
			item.setTextAlignment(QtCore.Qt.AlignHCenter)
			item.setToolTip(cd.tooltip)
			
			self.setHorizontalHeaderItem(col,item)
			for i in xrange(0,len(self.listed_names)):
				try : item = QtGui.QTableWidgetItem(cd.function(self.listed_names[i]))
				except : item = QtGui.QTableWidgetItem("-")
				item.setTextAlignment(QtCore.Qt.AlignHCenter)
				item.setFlags(flag3)
				if cd.lt_function: # This is how sort gets customized
					import new 
					item.__lt__ = new.instancemethod(cd.lt_function,item,QtGui.QTableWidgetItem)
				self.setItem(i, col, item)
			col += 1

		self.resizeColumnsToContents()
		self.setSortingEnabled(True)
		
	def add_entries(self,list_of_names):
		'''
		A dynamic way of adding entries to the table, for instance as a result of the user adding a file name from
		within the interface
		@param list_of_names the list of names to be added to the table
		Fills out all columns automatically
		'''
		sorting = self.isSortingEnabled()
		self.setSortingEnabled(False)
		r = self.rowCount()
		self.setRowCount(r+len(list_of_names))
		flag2 = Qt.ItemFlags(Qt.ItemIsSelectable)
		flag3 = Qt.ItemFlags(Qt.ItemIsEnabled)
		flag4 = Qt.ItemFlags(Qt.ItemIsEditable)
		new_items = []
		for i in xrange(0,len(list_of_names)):
			if self.icon != None: item = QtGui.QTableWidgetItem(self.icon,self.display_name(list_of_names[i]))
			else: item = QtGui.QTableWidgetItem(self.display_name(list_of_names[i]))
			item.setFlags(flag2|flag3)
			item.setTextAlignment(QtCore.Qt.AlignHCenter)
			self.setItem(r+i, 0, item)
			new_items.append(item)
			if self.single_selection and i == (len(list_of_names)-1):
				item.setSelected(True)
			
			for j, cd in enumerate(self.column_data):
				item = QtGui.QTableWidgetItem(cd.function(list_of_names[i]))
				item.setTextAlignment(QtCore.Qt.AlignHCenter)
				item.setFlags(flag3)
				if cd.lt_function: # This is how sort gets customized
					import new 
					item.__lt__ = new.instancemethod(cd.lt_function,item,QtGui.QTableWidgetItem)
				self.setItem(r+i,j+1, item)
			
		self.setSortingEnabled(sorting)
		if r == 0:
			self.resizeColumnsToContents()
	
		
		if self.items_selected_by_default:
			for item in new_items: 
				item.setSelected(True)

	def contextMenuEvent(self,event):
		'''
		Redefinition of QtGui.QTableWidget.contextMenuEvent
		Creates a context menu using self.context_menu_data, which is a dictionary
		@param event a QtGui.QContextMenuEvent - it is accepted
		'''
		menu = QtGui.QMenu()
		cmenu = self.context_menu_data
		for k in cmenu.keys():
			menu.addAction(k)
		QtCore.QObject.connect(menu,QtCore.SIGNAL("triggered(QAction*)"),self.menu_action_triggered)
		menu.exec_(event.globalPos())
		event.accept()
	
	def menu_action_triggered(self,action):
		'''
		Slot for context menu triggering
		@param action The context menu action that was triggered
		'''
		items = self.selectedItems()
		names = [str(item.text()) for item in items]
		self.context_menu_data[str(action.text())](names,self)

	def table_item_double_clicked(self,item):
		'''If an item is double clicked this function is called. Redefine it to 
		enable the displaying of a file, for example. You should probably check
		to make sure that the item is from the first column (column 0)
		'''
		pass
	
	def save_as(list_of_names,table_widget):
		'''
		Made static
		See the context menu dictionary in __init__, called when the user clicks "Save as"
		Iterates through the names and prompts the user for file names
		@param list_of_names a list table_widget column 0 entries, 
		@table_widget The table widget from which the entries were selected
		'''
		for name in list_of_names:
			from emsave import LightEMDataSave, save_data
			tmp = LightEMDataSave(table_widget.convert_text(name))
			val = save_data(tmp)
			if val == "":
				break
	save_as = staticmethod(save_as)
	
	def custom_addition(self,layout):
		'''
		When the associated table is being created and added in the form module, this function is called, enabling
		different things to be added to the layout in a custom fashion (such as an "Add" button).
		@param layout a Qt Layout (e.g. QVBoxLayout, QHBoxLayout - objects that support the 'addWidget' and 'addLayout' syntax
		'''
		for button_data in self.button_data:
			button = QtGui.QPushButton(button_data.name,None)
			layout.addWidget(button,0)
			QtCore.QObject.connect(button,QtCore.SIGNAL("clicked(bool)"),button_data.function)
			QtCore.QObject.connect(button,QtCore.SIGNAL("clicked(bool)"),self.sendupdate)
	def sendupdate(self):
		self.emit(QtCore.SIGNAL("updateform"))
		
	class EMColumnData:
		'''
		This class defines what's required to add column data to the EMFileTable
		'''
		def __init__(self,name,function,tooltip,lt_function=None):
			self.name = name # The name of the column of data
			self.function = function # The function which is called to populate the column with meta data - takes a file name as an argument, returns a string
			self.tooltip = tooltip # The helpful tooltip
			self.lt_function = lt_function # less than function - if specified is used as the operator< and sophisticates the sorting behavior
			
	class EMButtonData:
		'''
		This class defines what's required to add button data to the EMFileTable
		'''
		def __init__(self,name,function):
			self.name = name # The name of the button
			self.function = function # that which is called when the button is clicked (takes a single argument)
			
def float_lt(self,item2):
	'''
	This function is used to customize the sorting behavior of columns in the EMFileTable (and any subclasses)
	'''
	try:
		f = float(self.text())
	except:
		return 0
	
	try:
		f2 = float(item2.text())
	except:
		return 1
	
	return f > f2 

def int_lt(self,item2):
	'''
	This function is used to customize the sorting behavior of columns in the EMFileTable (and any subclasses)
	'''
	try:
		f = int(self.text())
	except:
		return 0
	
	try:
		f2 = int(item2.text())
	except:
		return 1
	
	return f > f2 

	

class EM2DFileTable(EMFileTable):
	def __init__(self,listed_names=[],name="filenames",desc_short="File Names",desc_long="A list of file names",single_selection=False):
		'''
		see EMFileTable for comments on parameters
		'''
		EMFileTable.__init__(self,listed_names,name,desc_short,desc_long,single_selection)
		self.icon = QtGui.QIcon(get_image_directory() + "/single_image.png")
		self.display_module = None
		self.module_events_manager = None
	
	def table_item_double_clicked(self,item):
		'''
		See EMFileTable.table_item_double_clicked for comments
		'''
		if item.column() != 0: return # only can display files from the first column
		filename = self.convert_text(str(item.text()))
		if not file_exists(filename): return # this happens sometimes when there is filtered data but no raw data
		get_application().setOverrideCursor(Qt.BusyCursor)
		if self.display_module == None:
			from emimage import EMWidgetFromFile
			self.display_module = EMWidgetFromFile(filename,get_application())
			from emapplication import ModuleEventsManager
			self.module_events_manager = ModuleEventsManager(self,self.display_module)
		else:
			from EMAN2 import EMData
			import emscene3d
			import emdataitem3d 
			
			data=emdataitem3d.EMDataItem3D(filename)
			self.display_module.insertNewNode(os.path.basename(filename), data, parentnode=self.display_module)
			isosurface = emdataitem3d.EMIsosurface(data)
			self.display_module.insertNewNode("Iso", isosurface, parentnode=data)
					
		#self.module().emit(QtCore.SIGNAL("launching_module"),"Browser",module)
		get_application().show_specific(self.display_module)
		#self.add_module([str(module),"Display",module])
		get_application().setOverrideCursor(Qt.ArrowCursor)

	def module_closed(self,module_instance):
		self.display_module = None
			

class EM3DFileTable(EM2DFileTable):
	'''
	Basically makes sure the correct icon is displayed.
	'''
	def __init__(self,listed_names=[],name="filenames",desc_short="File Names",desc_long="A list of file names",single_selection=False):
		'''
		see EMFileTable for comments on parameters
		'''
		EM2DFileTable.__init__(self,listed_names,name,desc_short,desc_long,single_selection)
		self.icon = QtGui.QIcon(get_image_directory() + "/single_image_3d.png")
	

class EMTomographicFileTable(EMFileTable):
	'''
	Basically makes sure the correct icon is displayed.
	Passes on image display - because tomographic images are huge
	'''
	def __init__(self,listed_names=[],name="filenames",desc_short="File Names",desc_long="A list of file names",single_selection=False,enable_save=False):
		'''
		see EMFileTable for comments on parameters
		'''
		EMFileTable.__init__(self,listed_names,name,desc_short,desc_long,single_selection,enable_save)
		self.icon = QtGui.QIcon(get_image_directory() + "/single_image_3d.png")
	
	def table_item_double_clicked(self,item):
		''' Just pass because it's probably a really big image'''
		pass

class EM2DStackTable(EMFileTable):
	def __init__(self,listed_names=[],name="filenames",desc_short="File Names",desc_long="A list of file names",single_selection=False,enable_save=True):
		'''
		see EMFileTable for comments on parameters
		'''
		EMFileTable.__init__(self,listed_names,name,desc_short,desc_long,single_selection,enable_save)
		self.icon = QtGui.QIcon(get_image_directory() + "/multiple_images.png")
		self.display_module = None
		self.module_events_manager = None
		if enable_save: self.context_menu_data["Save As"] = EM2DStackTable.save_as
		
#	def __del__(self):
#		print "stack table dies"
		
	def save_as(list_of_names,table_widget):
		'''
		Made static
		See the context menu dictionary in __init__, called when the user clicks "Save as"
		Iterates through the names and prompts the user for file names
		@param list_of_names a list table_widget column 0 entries, 
		@table_widget The table widget from which the entries were selected
		'''
		for name in list_of_names:
			if not file_exists(name):
				error("File %s doesn't exist" %s, "Error")
				continue
			from emimagemx import EMDataListCache
			tmp = EMDataListCache(table_widget.convert_text(name))
			from emsave import save_data
			val = save_data(tmp)
			if val == "":
				break
	save_as = staticmethod(save_as)
		
	def table_item_double_clicked(self,item):
		'''
		See EMFileTable.table_item_double_clicked for comments
		'''
		if item.column() != 0: return # only can display files from the first column
		filename = self.convert_text(str(item.text()))
		if not file_exists(filename): return # this happens sometimes when there is filtered data but no raw data
		get_application().setOverrideCursor(Qt.BusyCursor)
		if self.display_module == None:
			from emimage import EMWidgetFromFile
			self.display_module = EMWidgetFromFile(filename,get_application())
			from emapplication import ModuleEventsManager
			self.module_events_manager = ModuleEventsManager(self,self.display_module)
		else:
			from emimagemx import EMLightWeightParticleCache
			self.display_module.set_data(EMLightWeightParticleCache.from_file(filename)) #  I know this looks stupid, but c'est la vie
			self.display_module.updateGL()
					
		#self.module().emit(QtCore.SIGNAL("launching_module"),"Browser",module)
		get_application().show_specific(self.display_module)
		#self.add_module([str(module),"Display",module])
		get_application().setOverrideCursor(Qt.ArrowCursor)

	def module_closed(self,module_instance):
		self.display_module = None
		
class EMPlotTable(EMFileTable):
	def __init__(self,listed_names=[],name="filenames",desc_short="File Names",desc_long="A list of file names",single_selection=False):
		'''
		see EMFileTable for comments on parameters
		'''
		EMFileTable.__init__(self,listed_names,name,desc_short,desc_long,single_selection)
		self.icon = QtGui.QIcon(get_image_directory() + "/plot.png")
		self.display_module = None
		self.module_events_manager = None
	
#	def __del__(self):
#		print "2D table dies"
#	
	def table_item_double_clicked(self,item):
		'''
		See EMFileTable.table_item_double_clicked for comments
		'''
		if item.column() != 0: return # only can display files from the first column
		filename = self.convert_text(str(item.text()))
		if not file_exists(filename): return # this happens sometimes when there is filtered data but no raw data
		get_application().setOverrideCursor(Qt.BusyCursor)
		if self.display_module == None:
			from emimage import EMWidgetFromFile
			self.display_module = EMWidgetFromFile(filename,get_application())
			from emapplication import ModuleEventsManager
			self.module_events_manager = ModuleEventsManager(self,self.display_module)
		else:
			self.display_module.set_data_from_file(filename,True)
					
		#self.module().emit(QtCore.SIGNAL("launching_module"),"Browser",module)
		get_application().show_specific(self.display_module)
		#self.add_module([str(module),"Display",module])
		get_application().setOverrideCursor(Qt.ArrowCursor)

	def module_closed(self,module_instance):
		self.display_module = None
		
	def num_plot_entries(filename):
		try:
			f = file(filename,"r")
			lines = f.readlines()
			return str(len(lines))
		except:
			return "Error"
	num_plot_entries = staticmethod(num_plot_entries)

class EM2DStackExamineTable(EM2DStackTable):
	'''
	A specialized 2D stack table used by the Workflow 'Particle Examination and Stack Generation' stages, for
	the purpose of defining bad particles. The catch is, the chosen bad particles affect more than one particle stack,
	because we have filtered versions of image stacks. I.E. The user can define particles in any of the image stacks
	but should be able to observe the bad particles if they open the unfiltered particle stack, for example
	
	
	Also, inheriting for EM2DStackTable gives us save_as for free
	'''
	def __init__(self,listed_names=[],name="filenames",desc_short="File Names",desc_long="A list of file names",single_selection=False,name_map={}):
		'''
		see EMFileTable for comments on parameters
		'''
		EM2DStackTable.__init__(self,listed_names,name,desc_short,desc_long,single_selection)
		self.icon = QtGui.QIcon(get_image_directory() + "/multiple_images.png")
		self.display_module = None
		self.module_events_manager = None
		self.context_menu_data["Save As"] = EM2DStackTable.save_as
		self.name_map = name_map
		
	def set_name_map(self,name_map):
		self.name_map = name_map

	def table_item_double_clicked(self,item):
		'''
		See EMFileTable.table_item_double_clicked for comments
		'''
		print "X"
		if item.column() != 0: return # only can display files from the first column
		filename = self.convert_text(str(item.text()))
		if not file_exists(filename): return # this happens sometimes when there is filtered data but no raw data
		get_application().setOverrideCursor(Qt.BusyCursor)
		if self.display_module == None:
			from emimagemx import EMImageMXWidget
			self.display_module = EMImageMXWidget(None,get_application())
			from emapplication import ModuleEventsManager
			#self.module_events_manager = ModuleEventsManager(self,self.display_module)
		
		self.display_module.set_data(filename,filename) #  I know this looks stupid, but c'est la vie
		self.display_module.updateGL()
		if self.name_map.has_key(filename):
			self.display_module.set_single_active_set(self.name_map[filename])
		else:
			self.display_module.clear_sets()
		
		#self.module().emit(QtCore.SIGNAL("launching_module"),"Browser",module)
		get_application().show_specific(self.display_module)
		#self.add_module([str(module),"Display",module])
		get_application().setOverrideCursor(Qt.ArrowCursor)

	def module_closed(self,module_instance):
		self.display_module = None

class EMBrowseEventHandler:
	'''
	Base class for browse event handlers - came into existence because there are many different ways of handler the results
	of the browsing operation.
	It's a MIXIN because it actually supplies functionality
	May 2009: This class might be going out of fashion, since the SelectModule can be used as a dialog it is no longer necessary to
	intercept the "ok" and "cancel" signals.... instead just use exec_ ... 
	'''
	def __init__(self,browse_button):
		warnings.warn("EMBrowseEventHandler.__init__()", DeprecationWarning)
		self.browser = None
		self.browser_title = "Set this to be clear"
		QtCore.QObject.connect(browse_button,QtCore.SIGNAL("clicked(bool)"),self.browse_pressed)
		
	def browse_pressed(self,bool):
		if self.browser == None:
			self.browser = EMSelectorDialog(False, False)
			self.browser.setWindowTitle(self.browser_title)
			self.browser.exec_()
			QtCore.QObject.connect(self.browser,QtCore.SIGNAL("ok"),self.on_browser_ok)
			QtCore.QObject.connect(self.browser,QtCore.SIGNAL("cancel"),self.on_browser_cancel)
		else:
			self.browser.exec_()

	def on_browser_cancel(self):
		self.browser = None

	def on_browser_ok(self,string_list):
		'''
		The list is return from the browser, these are the things that user selected
		'''
		raise NotImplementedError("Browser ok not implemented. Inheriting class is supposed to do this")

		
def get_table_items_in_column(table_widget,column):
	'''
	Gets the table items from a particular column
	@param table_widget a QtGui.QTableWidget
	@param column the column from which you want to retrieve the table items
	@return a list of QtGui.QTableWidgetItemsW 
	'''
	r = table_widget.rowCount()
	entries = []
	for i in xrange(0,r):
		entries.append(table_widget.item(i,column))
		
	return entries
	
class EMEmanStrategyWidget(QtGui.QWidget):
	'''
	Something that knows how to automatically display and build parameters strings for 
	Eman strategy types, such as aligners, cmps, etc	
	'''
	def __init__(self,dumped_dict={},name="strategy",desc_short="Strategy",desc_long="Choose a strategy",defaultunits=None):
		QtGui.QWidget.__init__(self)
		self.name = name
		self.output_writers = [] # used to register output write objects, for the purpose of returning the results
		self.name_widget_map = {} # used to map parameter names to qt widgets - used by booleans to automatically disable and enable widgets
		self.dumped_dict = dumped_dict
		self.desc_short = desc_short
		self.desc_long = desc_long
		self.defaultunits = defaultunits
		self.vartype = "strategy"
		self.strategy_widget = {}
		self.strategy_output = {}
		self.current_strategy = None
		
		for key in self.dumped_dict.keys():
			self.strategy_widget[key] = None
		
		self.current_widget = None
		
		self.auto_incorporate = {}
		self.auto_incorporate["FLOAT"] = IncorpFloat()
		self.auto_incorporate["INT"] = IncorpInt()
		self.auto_incorporate["STRING"]= IncorpString()
		self.auto_incorporate["BOOL"]= IncorpBool()
		
		
	def build_widgets(self):
		
		self.vbl = QtGui.QVBoxLayout(self)
		self.dynamic_layout = QtGui.QVBoxLayout()
		groupbox = QtGui.QGroupBox(self.desc_short)
		groupbox.setToolTip(self.desc_long)
		
		self.main_combo = QtGui.QComboBox()
		start_idx = None
		dumped_dict_keys = self.dumped_dict.keys()
		dumped_dict_keys.sort()
		for i,key in enumerate(dumped_dict_keys):
			if key == self.defaultunits:
				start_idx = i
			self.main_combo.addItem(key)
		
		if start_idx == None and len(self.dumped_dict) > 0:
			start_idx = 1
			self.main_combo.setCurrentIndex(start_idx)
		
		self.dynamic_layout.addWidget(self.main_combo)
		groupbox.setLayout(self.dynamic_layout)
		
		self.vbl.addWidget(groupbox)
		
		QtCore.QObject.connect(self.main_combo, QtCore.SIGNAL("currentIndexChanged(QString)"), self.selection_changed)
		
		if start_idx != None:
			if start_idx != 0:
				self.main_combo.setCurrentIndex(start_idx) # this also serves as the purpose of populating the widgets into the interface on the first load
			else:
				self.selection_changed(self.main_combo.itemText(start_idx))
	
	def selection_changed(self,s):
		strategy = str(s)
		self.current_strategy = strategy
		if self.current_widget != None:
			self.dynamic_layout.removeWidget(self.current_widget)
			self.current_widget.hide()
			self.current_widget = None
			
		widget = self.strategy_widget[strategy]
		if widget == None:
			self.output_writers = []
			widget = self.get_strategy_widget(strategy)
			self.strategy_widget[strategy] = widget
			self.strategy_output[strategy] = self.output_writers
			
		else:widget.show()
		self.current_widget = widget
		
		self.dynamic_layout.addWidget(widget)
	
	def get_strategy_widget(self,strategy):
		data = self.dumped_dict[strategy]
		if (len(data) -1) % 3 != 0: raise RuntimeError("The format of the data is unknown") # first entry is descriptive text, then they should be in groups of threes, see dump_aligners_list, for example
		
		widget = QtGui.QWidget()
		vbl = QtGui.QVBoxLayout(widget)
		widget.setToolTip(data[0])
		params = []
		tmp_params = []
		for i in xrange(1,len(data),3):
			vartype = data[i+1]
			
			if self.auto_incorporate.has_key(vartype):
				name = data[i]
				desc_long = data[i+2]
				p = ParamDef(name=name,vartype=vartype,desc_short=name,desc_long=desc_long,defaultunits="")
				tmp_params.append(p)
			
				if len(tmp_params) == 3:
					params.append(tmp_params)
					tmp_params = []
			
			else: print "ignoring",vartype
				
		if len(tmp_params) != 0:
			params.append(tmp_params)
			tmp_params = []
		
		for param in params:
			if isinstance(param,list):
				hbl=QtGui.QHBoxLayout()
				for iparam in param:
					self.auto_incorporate[iparam.vartype](iparam,hbl,self)
				vbl.addLayout(hbl)
			else:
				self.auto_incorporate[param.vartype](param,vbl,self)
		
		return widget
	
	def write_data(self,dict):
		
		d = {}
		output_writers = self.strategy_output[self.current_strategy]
		if output_writers != None:
			for writer in output_writers:
				writer.write_data(d)
			
		
		result = self.current_strategy
		
		for key,val in d.items():
			if isinstance(val,str):
				if len(val) > 0:
					result += ":"+key+"="+val
			elif isinstance(val,bool):
				result += ":"+key+"="+str(int(val))
			else: # float, int
				result += ":"+key+"="+str(val)
				
		
		
		dict[self.name] = result


class EMFormWidget(QtGui.QWidget):
	'''
	See the example in __main__ below
	If ok is clicked the "emform_ok" signal is emitted along with a dictionary containing all of the form entries
	If cancel is clicked the "emform_cancel" signal is emmitted. No extra information is sent in this case
	'''
	def __init__(self,params=None,disable_ok_cancel=False):
		QtGui.QWidget.__init__(self,None)
		self.params = params
		self.event_handlers = [] # used to keep event handlers in memory
		self.resize_event_handlers = [] # used to keep resize event handlers in memory
		self.output_writers = [] # used to register output write objects, for the purpose of returning the results
		self.name_widget_map = {} # used to map parameter names to qt widgets - used by booleans to automatically disable and enable widgets
		self.__init_icons()
		self.setWindowIcon(QtGui.QIcon(get_image_directory() + "/desktop.png"))

		self.auto_incorporate = {}
		self.auto_incorporate["float"] = IncorpFloat()
		self.auto_incorporate["int"] = IncorpInt()
		self.auto_incorporate["url"]= IncorpUrl()
		self.auto_incorporate["choice"] = IncorpChoice()
		self.auto_incorporate["string"]= IncorpString()
		self.auto_incorporate["text"]= IncorpText()
		self.auto_incorporate["boolean"]= IncorpBool()
		self.auto_incorporate["stringlist"]= IncorpStringList()
		self.auto_incorporate["intlist"]= IncorpIntList()
		self.auto_incorporate["floatlist"]= IncorpFloatList()
		self.auto_incorporate["dict"]= IncorpDict()
		self.auto_incorporate["EMParamTable"]= IncorpParamTable()
		self.auto_incorporate["file_table"]= IncorpFileTable()
		self.auto_incorporate["EMButtonDialog"]= IncorpButtonDialog()
		self.auto_incorporate["strategy"]= IncorpStrategy()
		
		vbl = QtGui.QVBoxLayout()
		self.incorporate_params(self.params,vbl)
		if not disable_ok_cancel: self.__add_ok_cancel_buttons(vbl)
		self.setLayout(vbl)
		
		get_application().attach_child(self)

#	def __del__(self):
#		# this stops the QOBject : do not delete message
#		self.deleteLater()
		
	def __init_icons(self):
		self.emdata_icon = QtGui.QIcon(get_image_directory() + "/single_image.png")
		self.emdata_3d_icon = QtGui.QIcon(get_image_directory() + "/single_image_3d.png")
		self.emdata_matrix_icon = QtGui.QIcon(get_image_directory() + "/multiple_images.png")
		self.plot_icon = QtGui.QIcon(get_image_directory() + "/plot.png")
	
	def closeEvent(self, event):
		self.emit(QtCore.SIGNAL("emform_close"))
		QtGui.QWidget.closeEvent(self, event)
	
	def incorporate_params(self,params,layout):
		for param in params:
			act = True
			if isinstance(param,list):
				ftable = True
				for p in param:
					if not isinstance(p,EMFileTable):
						ftable = False
						break
				if ftable:
					self.incorporate_ftable_list(param,layout)
					act = False
				elif not isinstance(param,EMParamTable):
					hbl=QtGui.QHBoxLayout()
					for iparam in param:
#						if iparam.vartype == "EMButtonDialog": print "it was a button"
						self.auto_incorporate[iparam.vartype](iparam,hbl,self)
					layout.addLayout(hbl)
					act = False
					
			if act:
				self.auto_incorporate[param.vartype](param,layout,self)
		
		self.enable_boolean_dependents()
		
	def incorporate_list(self,param,layout,target,type_of):
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
	
		target.output_writers.append(ListWidgetParamWriter(param.name,list_widget,type_of))
		target.name_widget_map[param.name] = groupbox
	#layout.addLayout(hbl)
		
	def enable_boolean_dependents(self):
		''' This just takes care of making sure any dependent widgets are correctly enabled by default
		'''
		for event_handler in self.event_handlers:
			if isinstance(event_handler,BoolDependentsEventHandler):
				event_handler.checkbox_state_changed() # the argument is essentially unused
				
	def get_ptable_icon(self,EMParamTable):
		if hasattr(EMParamTable,"icon_type"):
			icon_type = EMParamTable.icon_type
			if icon_type == "single_image": return self.emdata_icon
			elif icon_type == "matrix_image": return self.emdata_matrix_icon
			elif icon_type == "3d_image": return self.emdata_3d_icon
			elif icon_type == "2d_plot" : return self.plot_icon
			else: return None
			
		return None
	
	def incorporate_ftable_list(self,file_table_list,layout):
		
		table = QtGui.QTabWidget()
		
		for paramtable in file_table_list:
			
		
			vbl=QtGui.QVBoxLayout()
			hbl=QtGui.QHBoxLayout()
			hbl.setMargin(0)
			hbl.setSpacing(2)
			
			paramtable.build_table()
			hbl.addWidget(paramtable)
			vbl.addLayout(hbl)
			
			paramtable.custom_addition(vbl)
			
#			groupbox = QtGui.QGroupBox(paramtable.desc_short)
#			groupbox.setToolTip(paramtable.desc_long)
			page = QtGui.QWidget()
			page.setLayout(vbl)
			if hasattr(paramtable,"icon"):
				table.addTab(page,paramtable.icon,paramtable.desc_short)
			else:
				table.addTab(page,paramtable.desc_short)
			
			self.output_writers.append(EMFileTableWriter(paramtable.name,paramtable,str))
			
		layout.addWidget(table,10)
	
	
	
	def incorporate_float_with_choices(self,param,layout,target):
		hbl=QtGui.QHBoxLayout()
		hbl.setMargin(0)
		hbl.setSpacing(2)
		
		hbl=QtGui.QHBoxLayout()
		label = QtGui.QLabel(param.desc_short,target)
		label.setToolTip(param.desc_long)
		hbl.addWidget(label)
				
		combo = QtGui.QComboBox(target)
		idx_default = 0
		for i,float_p in enumerate(param.choices):
			combo.addItem(str(float_p))
			
			if str(float_p) == str(param.defaultunits):
				idx_default = i
		combo.setCurrentIndex(idx_default)

		hbl.addWidget(combo,0)
		layout.addLayout(hbl)
		
		target.output_writers.append(FloatChoiceParamWriter(param.name,combo))
		target.name_widget_map[param.name] = [combo,label]
	
	def incorporate_int_with_choices(self,param,layout,target):
		hbl=QtGui.QHBoxLayout()
		hbl.setMargin(0)
		hbl.setSpacing(2)
		
		hbl=QtGui.QHBoxLayout()
		label = QtGui.QLabel(param.desc_short,target)
		label.setToolTip(param.desc_long)
		hbl.addWidget(label)
				
		combo = QtGui.QComboBox(target)
		idx_default = 0
		for i,integer in enumerate(param.choices):
			combo.addItem(str(integer))
			if str(integer) == str(param.defaultunits):
				idx_default = i
		combo.setCurrentIndex(idx_default)

		hbl.addWidget(combo,0)
		layout.addLayout(hbl)
		
		target.output_writers.append(IntChoiceParamWriter(param.name,combo))
		target.name_widget_map[param.name] = [combo,label]
	
	def incorporate_string_with_choices(self,param,layout,target):
		hbl=QtGui.QHBoxLayout()
		hbl.setMargin(0)
		hbl.setSpacing(2)
		
		hbl=QtGui.QHBoxLayout()
		label = QtGui.QLabel(param.desc_short,target)
		label.setToolTip(param.desc_long)
		hbl.addWidget(label,0)
				
		combo = QtGui.QComboBox(target)
		idx_default = 0
		for i,string in enumerate(param.choices):
			combo.addItem(string)
			if string == param.defaultunits:
				idx_default = i
		combo.setCurrentIndex(idx_default)

		hbl.addWidget(combo,0)
		layout.addLayout(hbl)
		
		target.output_writers.append(StringChoiceParamWriter(param.name,combo))
		target.name_widget_map[param.name] = [combo,label]
	
	def __add_ok_cancel_buttons(self,layout):
		hbl=QtGui.QHBoxLayout()
		label = QtGui.QLabel("Form commands:")
		hbl.addWidget(label)
		
		ok_button = QtGui.QPushButton("Ok")
		ok_button.setToolTip("When you click ok the values in the form are sent to the calling program")
		hbl.addWidget(ok_button)
		cancel_button = QtGui.QPushButton("Cancel")
		hbl.addWidget(cancel_button,0)
		layout.addLayout(hbl)
		QtCore.QObject.connect(ok_button,QtCore.SIGNAL("clicked(bool)"),self.ok_pressed)
		QtCore.QObject.connect(cancel_button,QtCore.SIGNAL("clicked(bool)"),self.cancel_pressed)
		
	def ok_pressed(self,bool):
		ret = {}
		for output in self.output_writers: output.write_data(ret)
		self.emit(QtCore.SIGNAL("emform_ok"),ret)
		
	def cancel_pressed(self,bool):
		self.emit(QtCore.SIGNAL("emform_cancel"))


	def update_texture(self):
		pass


	def display_file(self,filename):
		self.emit(QtCore.SIGNAL("display_file"),filename)


class IncorpStrategy:
	def __init__(self): pass
	def __call__(self,strategy,layout,target=None):
		num_choices = None
		strategy.build_widgets()
		layout.addWidget(strategy)
		if target != None: target.output_writers.append(strategy)
		
class IncorpButtonDialog:
	def __init__(self): pass
	def __call__(self,buttondialog,layout,target):
		buttondialog.set_target(target)
		buttondialog.add_to_layout(layout)
		target.event_handlers.append(buttondialog)    
		

class IncorpFileTable:
	def __init__(self): pass
	def __call__(self,paramtable,layout,target=None):
		num_choices = None
		# first check that there are no inconsistencies in the number of parameter choices
		
		vbl=QtGui.QVBoxLayout()
		hbl=QtGui.QHBoxLayout()
		hbl.setMargin(0)
		hbl.setSpacing(2)
		
		paramtable.build_table()
		hbl.addWidget(paramtable)
		vbl.addLayout(hbl)
		
		paramtable.custom_addition(vbl)
		
		groupbox = QtGui.QGroupBox(paramtable.desc_short)
		groupbox.setToolTip(paramtable.desc_long)
		
		#if hasattr(paramtable,"icon"):groupbox.setWindowIcon(paramtable.icon) #it would be nice if this worked
		
		
		groupbox.setLayout(vbl)
		layout.addWidget(groupbox,10)
		if target != None: target.output_writers.append(EMFileTableWriter(paramtable.name,paramtable,str))

class IncorpParamTable:
	def __init__(self): pass
	def __call__(self,paramtable,layout,target):
	
		num_choices = None
		# first check that there are no inconsistencies in the number of parameter choices
		for param in paramtable:
			if num_choices == None:
				num_choices = len(param.choices)
			else:
				if len(param.choices) != num_choices:
					print "error, the number of choices is not consistent in __incorporate_paramtable"
					return
		
		vbl=QtGui.QVBoxLayout()
		hbl=QtGui.QHBoxLayout()
		hbl.setMargin(0)
		hbl.setSpacing(2)
		
		vbl.addLayout(hbl)
		
		table_widget = QtGui.QTableWidget(num_choices, len(paramtable), None)
		icon = target.get_ptable_icon(paramtable)
		
		if not paramtable.enable_multiple_selection:
			table_widget.setSelectionMode(QtGui.QAbstractItemView.SingleSelection)
			
		selected_items = [] # used to ensure default selection is correct
		
		table_widget.setSortingEnabled(False)
		
		paramtable.add_optional_table_attr(table_widget)
		paramtable.build_table(table_widget,icon)			
		table_widget.setSortingEnabled(True)
		table_widget.setToolTip(paramtable.desc_long)
		hbl.addWidget(table_widget,1)
		
		groupbox = QtGui.QGroupBox(paramtable.desc_short)
		groupbox.setToolTip(paramtable.desc_long)
		groupbox.setLayout(vbl)
		layout.addWidget(groupbox,10)
		
		table_event_handler = EMParamTableEventHandler(target,table_widget)
		target.event_handlers.append(table_event_handler)
		target.resize_event_handlers.append(table_event_handler)
		table_widget.resizeColumnsToContents()
		
		if len(paramtable[0].choices) > 0:
			type_of = type(paramtable[0].choices[0])
		else:
			type_of = str # This case arise when the table is empty and the user had an add button.
	
#		for item in selected_items: 
#			item.setSelected(True)
		target.output_writers.append(EMParamTableWriter(paramtable.name,table_widget,type_of))
		
		paramtable.custom_addition(vbl,table_widget)
		
		target.name_widget_map[paramtable.name] = groupbox
		
class IncorpStringList:
	def __init__(self): pass
	def __call__(self,param,layout,target):
		target.incorporate_list(param,layout,target,str)

class IncorpFloatList:
	def __init__(self): pass
	def __call__(self,param,layout,target):
		target.incorporate_list(param,layout,target,float)

class IncorpIntList:
	def __init__(self): pass
	def __call__(self,param,layout,target):
		target.incorporate_list(param,layout,target,int)



class IncorpBool:
	def __init__(self): pass
	def __call__(self,param,layout,target):
		hbl=QtGui.QHBoxLayout()
		hbl.setMargin(0)
		hbl.setSpacing(2)
		check_box = QtGui.QCheckBox(str(param.desc_short),target)
		check_box.setChecked(bool(param.defaultunits))
		check_box.setToolTip(param.desc_long)
		hbl.addWidget(check_box,0)
		layout.addLayout(hbl)
		target.output_writers.append(BoolParamWriter(param.name,check_box))
		target.name_widget_map[param.name] = check_box
		
		if hasattr(param,"dependents"):
			target.event_handlers.append(BoolDependentsEventHandler(target,check_box,param.dependents,hasattr(param,"invert_logic") and param.invert_logic))    

class IncorpString:
	def __init__(self): pass
	def __call__(self,param,layout,target):
		if param.choices != None and len(param.choices) > 1:
			target.incorporate_string_with_choices(param,layout,target)
		else:
			hbl=QtGui.QHBoxLayout()
			hbl.setMargin(0)
			hbl.setSpacing(2)
			label = QtGui.QLabel(param.desc_short+":",target)
			label.setToolTip(param.desc_long)
			hbl.addWidget(label)
			line_edit = QtGui.QLineEdit(str(param.defaultunits),target)
			hbl.addWidget(line_edit,0)
			hbl.name = param.name
			layout.addLayout(hbl)
			target.output_writers.append(StringParamWriter(param.name,line_edit))
			target.name_widget_map[param.name] =  [line_edit,label]
	
class IncorpFloat:
	def __init__(self): pass
	def __call__(self,param,layout,target):
		if param.choices != None and len(param.choices) > 1:
			target.incorporate_float_with_choices(param,layout,target)
		else:
			hbl=QtGui.QHBoxLayout()
			hbl.setMargin(0)
			hbl.setSpacing(2)
			label = QtGui.QLabel(param.desc_short+":",target)
			label.setToolTip(param.desc_long)
			hbl.addWidget(label)
			double_validator = QtGui.QDoubleValidator(target)
			line_edit = QtGui.QLineEdit(str(param.defaultunits),target)
			line_edit.setValidator(double_validator)
			hbl.addWidget(line_edit,0)
			hbl.name = param.name
			layout.addLayout(hbl)
			target.output_writers.append(FloatParamWriter(param.name,line_edit))
			target.name_widget_map[param.name] = [line_edit,label]
	
class IncorpInt:
	def __init__(self): pass
	def __call__(self,param,layout,target):
		if param.choices != None and len(param.choices) > 1:
			target.incorporate_int_with_choices(param,layout,target)
		else:
			hbl=QtGui.QHBoxLayout()
			hbl.setMargin(0)
			hbl.setSpacing(2)
			label = QtGui.QLabel(param.desc_short+":",target)
			label.setToolTip(param.desc_long)
			hbl.addWidget(label)
			pos_int_validator = QtGui.QIntValidator(target)
			line_edit = QtGui.QLineEdit("",target)
			line_edit.setValidator(pos_int_validator)
			line_edit.setText(str(param.defaultunits))
			hbl.addWidget(line_edit,0)
			hbl.name = param.name
			layout.addLayout(hbl)
			target.output_writers.append(IntParamWriter(param.name,line_edit))
			target.name_widget_map[param.name] = [line_edit,label]
	
class IncorpText:
	def __init__(self): pass
	def __call__(self,param,layout,target):
#			hbl=QtGui.QHBoxLayout()
#			hbl.setMargin(0)
#			hbl.setSpacing(2)
				
		text_edit = QtGui.QTextEdit("",target)
		text_edit.setReadOnly(True)
		text_edit.setWordWrapMode(QtGui.QTextOption.WordWrap)
		text_edit.setText(param.defaultunits)
		#text_edit.setSizePolicy(QtGui.QSizePolicy.Minimum,QtGui.QSizePolicy.Minimum)
#			hbl.addWidget(text_edit,0)
		
#			groupbox = QtGui.QGroupBox(param.desc_short)
#			groupbox.setToolTip(param.desc_long)
#			groupbox.setLayout(hbl)
		
#			layout.addWidget(groupbox,0)
		layout.addWidget(text_edit,1)

		target.output_writers.append(TextParamWriter(param.name,text_edit))	
		#target.name_widget_map[param.name] = groupbox
	
class IncorpUrl:
	def __init__(self): pass
	def __call__(self,param,layout,target):
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
				
		text_edit = QtGui.QTextEdit("",target)
		text_edit.setWordWrapMode(QtGui.QTextOption.NoWrap)
		text_edit.setText(defaults)
		hbl.addWidget(text_edit,0)
		vbl.addLayout(hbl)
		
		hbl2=QtGui.QHBoxLayout()
		hbl2.setMargin(0)
		hbl2.setSpacing(2)
		
		browse_button = QtGui.QPushButton("Browse",target)
		hbl2.addWidget(browse_button)
		clear_button = QtGui.QPushButton("Clear",target)
		hbl2.addWidget(clear_button,0)
		vbl.addLayout(hbl2)
		
		groupbox = QtGui.QGroupBox(param.desc_short)
		groupbox.setToolTip(param.desc_long)
		groupbox.setLayout(vbl)
		
		layout.addWidget(groupbox,0)
		
		target.event_handlers.append(UrlEventHandler(target,text_edit,browse_button,clear_button,get_application(),param.desc_short))
		target.output_writers.append(UrlParamWriter(param.name,text_edit))
		target.name_widget_map[param.name] = groupbox

class IncorpDict:
	def __init__(self): pass
	def __call__(self,param,layout,target):
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
#		label = QtGui.QLabel(param.desc_short+":",target)
#		label.setToolTip(param.desc_long)
#		hbl.addWidget(label)
		
		combo = QtGui.QComboBox(target)
		idx_default = 0
		for i,k in enumerate(keys):
			combo.addItem(str(k))
			
			if str(k) == str(param.defaultunits):
				idx_default = i
				
		combo.setCurrentIndex(idx_default)
		combo_default = keys[idx_default]
		hbl.addWidget(combo)
		
		
		combo2 = QtGui.QComboBox(target)
		for v in param.choices[combo_default]:
			combo2.addItem(str(v))
			
		hbl.addWidget(combo2)
		
		groupbox = QtGui.QGroupBox(param.desc_short)
		groupbox.setToolTip(param.desc_long)
		groupbox.setLayout(hbl)
		
		layout.addWidget(groupbox)
		
		target.output_writers.append(DictParamWriter(param,combo,combo2))
		target.event_handlers.append(DictEventHandler(param.choices,combo,combo2))
		
		target.name_widget_map[param.name] = groupbox

class IncorpChoice:
	def __init__(self): pass
	def __call__(self,param,layout,target):
	
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
		target.output_writers.append(ChoiceParamWriter(param.name,buttons,type(param.choices[0])))
		
		target.name_widget_map[param.name] = [groupbox,buttons]

class EMParamTableWriter:
	def __init__(self,param_name,table_widget,type_of):
		self.param_name = param_name
		self.table_widget = table_widget
		self.type_of = type_of
		
	def write_data(self,dict):
		sel = [self.type_of(item.text()) for item in self.table_widget.selectedItems()]
		dict[self.param_name] = sel
		
class EMFileTableWriter:
	def __init__(self,param_name,table_widget,type_of=str):
		self.param_name = param_name
		self.table_widget = table_widget
		self.type_of = type_of
		
	def write_data(self,dict):
		sel = [self.table_widget.convert_text(self.type_of(item.text())) for item in self.table_widget.selectedItems()]
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

class EMParamTableEventHandler:
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
		self.table_widget.context_menu[str(action.text())](names,self.table_widget)
	
	def table_item_clicked(self,item):
		#
		pass
	def resize(self):
		return
		cols = self.table_widget.columnCount()
		cumulative_width = 0
		for i in range(cols):
			cumulative_width += self.table_widget.columnWidth(i)
		
		tab_widget_width = self.table_widget.geometry().width()
		#print cumulative_width,self.table_widget.width(),tab_widget_width,self.table_widget.frameSize().width()
		#print self.table_widget.getContentsMargins()
		
		if cumulative_width < self.table_widget.width():
			scale = float(self.table_widget.width())/cumulative_width
			for i in range(cols):
				w = self.table_widget.columnWidth(i)
				self.table_widget.setColumnWidth(i,w*scale)
		
		
		
class UrlEventHandler(EMBrowseEventHandler):
	'''
	The browse and cancel events have to be sent to the correct line edit, so this handles it
	'''
	def __init__(self,target,text_edit,browse_button,clear_button,application,title=""):
		self.target = weakref.ref(target)
		self.application = weakref.ref(application)
		self.text_edit = text_edit
		EMBrowseEventHandler.__init__(self,browse_button)
		self.browser_title = title
		
		QtCore.QObject.connect(clear_button,QtCore.SIGNAL("clicked(bool)"),self.clear_pressed)
		
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
	def __init__(self,target,checkbox,dependents,invert_logic=False):
		'''
		Need the target because we need access to the name_widget_map
		'''
		self.target = weakref.ref(target)
		self.checkbox = checkbox
		self.dependents = dependents # a list of depent names (ParamDef.name)
		self.invert_logic = invert_logic
		QtCore.QObject.connect(self.checkbox, QtCore.SIGNAL("stateChanged(int)"),self.checkbox_state_changed)
		
	def checkbox_state_changed(self,integer=0):
		name_map = self.target().name_widget_map
		
		enabled = False
		if self.checkbox.isChecked(): enabled = True
		
		if self.invert_logic: enabled = not enabled
		
		for name in self.dependents:
			widgets = name_map[name]
			if isinstance(widgets,list):
				for widget in widgets: widget.setEnabled(enabled)
			else: widgets.setEnabled(enabled)
		
		
class EMTableFormWidget(EMFormWidget):
	'''
	See the example in __main__ below
	If ok is clicked the "emform_ok" signal is emitted along with a dictionary containing all of the form entries
	If cancel is clicked the "emform_cancel" signal is emmitted. No extra information is sent in this case
	'''
	def __init__(self,params=None):
		EMFormWidget.__init__(self,params)
	
#	def __del__(self): print "del table form widget"
		
	def incorporate_params(self,params,layout):
		tabwidget = QtGui.QTabWidget(self)
		
		for title,paramlist in params:
			
			widget = QtGui.QWidget(None)
			vbl =  QtGui.QVBoxLayout(widget)
			#print paramlist
			EMFormWidget.incorporate_params(self,paramlist,vbl)
#			for param in paramlist:
#				EMFormWidget.incorporate_params(self, params, vbl)
#				if isinstance(param,list) and len(param) != 1:
#					hbl=QtGui.QHBoxLayout()
#					for iparam in param:
#						try:
#							self.auto_incorporate[iparam.vartype](iparam,hbl,self)
#						except:
#							if iparam.vartype == "EMButtonDialog": 
#								self.auto_incorporate[iparam.vartype](iparam,hbl,self)
#					vbl.addLayout(hbl)
#					continue
#				else:
#					self.auto_incorporate[param.vartype](param,vbl,self)
#			
			tabwidget.addTab(widget,title)
		
		layout.addWidget(tabwidget)
		
		self.enable_boolean_dependents()
				

def get_example_form_params():
	params = []
	pstrategy = EMEmanStrategyWidget(dump_aligners_list(),name="align",desc_short="Strategy",desc_long="Choose a strategy",defaultunits=None)
	params.append(pstrategy)
	pstrategy2 = EMEmanStrategyWidget(dump_processors_list(),name="proc",desc_short="Processors",desc_long="Choose a processor",defaultunits=None)
	params.append(pstrategy2)
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
	a = EMParamTable(name="table_choice",desc_short="Please choose using this information",desc_long="The left most column is what you're choosing from, the extra columns are used only to assist in the decision making process")
	a.append(pil)
	a.append(pfl)
	params.append([pil,pfl,a])
	
	return params


def get_small_example_form_params():
	params = []
	
	params.append(ParamDef(name="box size",vartype="int",desc_short="int",desc_long="An integer value",property=None,defaultunits=128,choices=[]))
	params.append(ParamDef(name="apix",vartype="float",desc_short="float",desc_long="A floating point value",property=None,defaultunits=1.0,choices=[]))
	ps = ParamDef(name="model",vartype="string",desc_short="string",desc_long="A string value",property=None,defaultunits="Three thousand series",choices=None)
	params.append(ParamDef(name="Is there a God?",vartype="boolean",desc_short="boolean",desc_long="Something that is true or false",property=None,defaultunits=Util.get_irand(0,1),choices=None))
	p1 = ParamDef(name="True or false",vartype="boolean",desc_short="boolean",desc_long="Something that is true or false",property=None,defaultunits=Util.get_irand(0,1),choices=None)
	p2 = ParamDef(name="integer length",vartype="int",desc_short="int",desc_long="An integer value",property=None,defaultunits=128,choices=[])
	#params.append([p1,p2,ps])
	#params.append(ParamDef(name="Normalization method",vartype="string",desc_short="string",desc_long="Choose from a list of strings",property=None,defaultunits="normalize",choices=["normalize","normalize.edgemean","normalize.other"]))
	
	pi = ParamDef(name="Int 1 to 10",vartype="int",desc_short="int",desc_long="Choose from a list of ints",property=None,defaultunits=5,choices=[1,2,3,4,5,6,7,8,9,10])
	pf = ParamDef(name="Float 1 to 10",vartype="float",desc_short="float",desc_long="Choose from a list of floats",property=None,defaultunits=2.1,choices=[1.1,2.1,3.1,4.1,5.1,6.1,7.1,8.1,9.1,10.1])
	#params.append([pi,pf])
	#params.append(ParamDef(name="file_names",vartype="url",desc_short="url",desc_long="This is an editable list of file names with convenient browse and clear buttons",property=None,defaultunits=["tmp.txt","other.mrc"],choices=[]))
	params.append(ParamDef(name="comparator",vartype="choice",desc_short="choice",desc_long="This is a string choice",property=None,defaultunits="frc",choices=["frc","phase","sqeuclidean"]))
	params.append(ParamDef(name="lucky number",vartype="choice",desc_short="choice",desc_long="This is to demonstrate that the type of choice is preserved. When you click ok the value in the return dictionary corresponding to this form entry will be an integer",property=None,defaultunits=3,choices=[1,2,3,98]))
	
	#params.append(ParamDef(name="song",vartype="text",desc_short="text",desc_long="A potentially very long description",property=None,defaultunits="Jack and Jill went up the hill\nTo fetch a pail of water.\nJack fell down and broke his crown,\nAnd Jill came tumbling after.",choices=None))
	
	params.append(ParamDef(name="Selected Files",vartype="stringlist",desc_short="stringlist",desc_long="Choose from a list of strings",property=None,defaultunits=["C.mrc","E.mrc"],choices=[chr(i)+".mrc" for i in range(65,91)]))
	
	params.append(ParamDef(name="Combo1 and combo2", vartype="dict",desc_short="dict",desc_long="Choose from the combo on the left and the combo box on the right will be updated. The 'and' in the the name differentiates the key in the return dictionary", property=None, defaultunits="reconstructors", choices={"reconstructors":["a","b","c"],"processors":["1","2","3"]} ))
	
	pil = ParamDef(name="Int 1 to 10 from a list",vartype="intlist",desc_short="intlist",desc_long="Choose from a list of ints",property=None,defaultunits=[5],choices=[1,2,3,4,5,6,7,8,9,10])
	pfl = ParamDef(name="Float 1 to 10 from a list",vartype="floatlist",desc_short="floatlist",desc_long="Choose from a list of floats",property=None,defaultunits=[2.1],choices=[1.1,2.1,3.1,4.1,5.1,6.1,7.1,8.1,9.1,10.11111111111111111111111111111])
	a = EMParamTable(name="table_choice",desc_short="Please choose using this information",desc_long="The left most column is what you're choosing from, the extra columns are used only to assist in the decision making process")
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
	
	from emapplication import EMApp
	em_app = EMApp()
	window = EMFormWidget(params=get_example_form_params())
	window.setWindowTitle("A test form")
	QtCore.QObject.connect(window,QtCore.SIGNAL("emform_ok"),on_ok)
	QtCore.QObject.connect(window,QtCore.SIGNAL("emform_cancel"),on_cancel)
	
	window2= EMTableFormWidget(params=get_example_table_form_params())
	window2.setWindowTitle("A test form")
	QtCore.QObject.connect(window2,QtCore.SIGNAL("emform_ok"),on_ok)
	QtCore.QObject.connect(window2,QtCore.SIGNAL("emform_cancel"),on_cancel)
	
	em_app.show()
	em_app.execute()
	
	
