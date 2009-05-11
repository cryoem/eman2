#!/usr/bin/env python

#
# Author: David Woolford (woolford@bcm.edu)
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

import PyQt4
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
import os
import re
from EMAN2 import get_image_directory,e2getcwd,get_dtag,EMData,get_files_and_directories,db_open_dict,strip_file_tag,remove_file
from EMAN2 import remove_directories_from_name,Util,EMUtil,IMAGE_UNKNOWN,get_file_tag,db_check_dict,file_exists, base_name
from emimage2d import EMImage2DModule
from emapplication import EMStandAloneApplication, EMQtWidgetModule, EMProgressDialogModule, get_application
from EMAN2db import EMAN2DB
from EMAN2db import db_convert_path
from emplot2d import EMPlot2DModule
from emapplication import EMQtWidgetModule, ModuleEventsManager
import weakref
import copy
from emsave import save_data

read_header_only = True

class EMSelectorModule(EMQtWidgetModule):
	def __init__(self,single_selection=False,save_as_mode=True):
		self.widget = EMSelectorDialog(self,single_selection,save_as_mode)
		EMQtWidgetModule.__init__(self,self.widget)
		
#	def __del__(self):
#		import sys
#		print "selector module death", sys.getrefcount(self.widget)
		
	def exec_(self):
		'''
		Wraps self.widget.exec_
		@return a list of selected filenames
		'''
		return self.widget.exec_()
		

class EMSelectorDialog(QtGui.QDialog):
	'''
	This class is somewhat of an aberration. It can be used both as a file system browser
	and as a save-file dialog. To see an example of using it as a system browser, see the
	EMBrowserModule. To see an example of using it as a save-file dialog, see emsave.py.
	Also, emsprworkflow uses this class as a browser for selecting files.  
	I wish I had time to correct the design... but that said, it's not too far from where it should
	be - the code is relatively clean, and correcting the design involves breaking some of the member
	functions into separate classes.
	'''
	def __init__(self,module,single_selection=False,save_as_mode=True):
		'''
		@param module - should be an EMSelectorModule
		@param application - should be a type of application as supplied in emapplication
		@param single_selection - should selections be limited to singles?
		@param save_as_mode - should be True if using this object as a pure dialog (using exec_).
		'''
		QtGui.QDialog.__init__(self,None)
		self.setFocusPolicy(Qt.StrongFocus)
		self.module=weakref.ref(module)
		self.desktop_hint = "dialog"
		self.single_selection = single_selection
		self.db_listing = EMDBListing(self)
		self.dir_listing = EMFSListing(self)
		
		self.hbl = QtGui.QVBoxLayout(self)
		self.hbl.setMargin(0)
		self.hbl.setSpacing(6)
		self.hbl.setObjectName("hbl")
		
		self.__init_icons()
		
		self.__init_filter_combo()
		
		self.first_list_widget = QtGui.QListWidget(None)
		self.starting_directory = e2getcwd()
		self.historical_starting_directory = e2getcwd() # just incase anyone ever needs it (True). This should never be changed
		
		self.current_force = None # use to keep track of what is currently being forced, either 2D plotting, 2D image showing, or neither
		self.selections = []
		self.current_list_widget = None
		self.lock = True
		self.list_widgets = []
		self.previews = [] # keeps track of all of the preview windows
		self.module_events = [] # used to coordinate signals from the modules, especially close events, to free memory
		self.list_widget_data= [] # entries should be tuples containing (current folder item)
		self.splitter = QtGui.QSplitter()
		self.__add_list_widget(self.first_list_widget)
		self.__add_list_widget()
		if not save_as_mode: self.__add_list_widget() # this is just to distinguish the saving mode from the browsing mode
		
		self.hbl.addWidget(self.splitter,1)
		
		self.__load_directory_data(self.starting_directory,self.first_list_widget)
		#self.first_list_widget.setCurrentRow(-1)

		self.save_as_mode = False
		if save_as_mode:
			hbl2=QtGui.QHBoxLayout()
			hbl2.setMargin(0)
			hbl2.setSpacing(2)
			label = QtGui.QLabel("Save as",self)
			hbl2.addWidget(label)
			self.save_as_line_edit = QtGui.QLineEdit("",self)
			hbl2.addWidget(self.save_as_line_edit,0)
			self.hbl.addLayout(hbl2)
			self.save_as_mode = True
			self.validator = None # an optional feature which makes the behavior of the dialog more sophisticated - see emsave.py
			self.dialog_result = ""

		self.bottom_hbl = QtGui.QHBoxLayout()
		self.bottom_hbl.addWidget(self.filter_text,0)
		self.bottom_hbl.addWidget(self.filter_combo,1)
		self.__init_buttons()
		self.bottom_hbl.addWidget(self.cancel_button,0)
		self.bottom_hbl.addWidget(self.ok_button,0)
		self.hbl.addLayout(self.bottom_hbl)
		
		if not save_as_mode:
			self.gl_image_preview = None
			
			self.bottom_hbl2 = QtGui.QHBoxLayout()
			self.__init_preview_options()
			self.bottom_hbl2.addWidget(self.preview_options,0)
			self.hbl.addLayout(self.bottom_hbl2)
			
			self.__init__force_2d_tb()
			self.__init__force_plot_tb()
			self.bottom_hbl2.addWidget(self.force_2d,0)
			self.bottom_hbl2.addWidget(self.force_plot,0)
			
			self.bottom_hbl3 = QtGui.QHBoxLayout()
			self.__init_plot_options()
			self.bottom_hbl3.addWidget(self.replace,0)
			self.bottom_hbl3.addWidget(self.include,0)
			
			#self.hbl.addLayout(self.bottom_hbl3)
			
			self.groupbox = QtGui.QGroupBox("Plot/3D options")
			self.groupbox.setLayout(self.bottom_hbl3)
			self.groupbox.setEnabled(False)
			
			self.bottom_hbl2.addWidget(self.groupbox)
		
		if not save_as_mode:
			self.resize(480,480)
		else:
			self.resize(480,400)
		
		self.lock = False
		
		self.paint_events = 0
		
		self.timer_interval = 500
		self.timer = QtCore.QTimer()
		QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.time_out)
		
		self.timer.start(self.timer_interval)
		
		self.selected_files = []
		
	def set_validator(self,validator):
		'''
		Sets the validator which is used only in save_as mode (when this is being used as a file selecting dialog)
		An example validator is an emsave.EMSaveImageValidator
		'''
		self.validator = validator
	
	def exec_(self):
		'''
		Wraps QtGui.QDialog.exec_
		@return a list of selected filenames
		'''
		QtGui.QDialog.exec_(self)
		return self.dialog_result
		
	
	def time_out(self):
		'''
		This function takes care of automatic updates - if the file system changes then so does
		the information we display
		'''
		for i,widget in enumerate(self.list_widgets):
			if hasattr(widget,"directory"):
				statinfo = os.stat(widget.directory)
				if statinfo[-2] != widget.mod_time:
					crow = widget.currentRow()
					old_row_text = widget.item(crow)
					if old_row_text != None: old_row_text = old_row_text.text()
					widget.clear()
					widget.directory_lister.load_directory_data(widget.directory_arg,widget)
					if i == 0:
							a = QtGui.QListWidgetItem(self.up_arrow_icon,"../",None)
							a.type_of_me = "go up a directory"
							widget.insertItem(0,a)
					new_row_text = widget.item(crow) # this could potentially by None
					if new_row_text != None: new_row_text = new_row_text.text()
					# if the old_row_text is is None then nothing was selected and we don't have to worry about clearing the other widgets, which may have been displaying metadata
					if old_row_text != None and new_row_text != old_row_text:
						#print "row lost"
						# if the current row was deleted then the next list widgets need to be cleared
						for j in range(i+1,len(self.list_widgets)):
							widget = self.list_widgets[j]
							widget.clear()
							self.list_widget_data[j] = None
							# get rid of these attributes, they could have been set by the listing objects
							del_attr = ["directory","directory_lister","directory_arg","mod_time"]
							for attr in del_attr:
								if hasattr(widget,attr): delattr(widget,attr)
				
					return # only do one at a time
				
	def get_desktop_hint(self):
		return self.desktop_hint
		
	def __init_icons(self):
		directory = get_image_directory()
		self.setWindowIcon(QtGui.QIcon(directory + "/display_icon.png"))
		self.folder_icon = QtGui.QIcon(directory + "/Folder.png")
		self.folder_files_icon = QtGui.QIcon(directory + "/FolderFiles.png")
		self.file_icon = QtGui.QIcon(directory + "/File.png")
		self.database_icon = QtGui.QIcon(directory + "/database.png")
		self.key_icon = QtGui.QIcon(directory + "/Key.png")
		self.basic_python_icon = QtGui.QIcon(directory + "/boxhabanosclose.png")
		self.dict_python_icon = QtGui.QIcon(directory + "/Bag.png")
		self.ab_refboxes_icon = QtGui.QIcon(directory + "/black_box.png")
		self.ab_manboxes_icon = QtGui.QIcon(directory + "/black_box.png")
		self.ab_autoboxes_icon = QtGui.QIcon(directory + "/green_boxes.png")
		self.emdata_icon = QtGui.QIcon(directory + "/single_image.png")
		self.emdata_3d_icon = QtGui.QIcon(directory + "/single_image_3d.png")
		self.emdata_matrix_icon = QtGui.QIcon(directory + "/multiple_images.png")
		self.emdata_matrix_3d_icon = QtGui.QIcon(directory + "/multiple_images_3d.png")
		self.up_arrow_icon = QtGui.QIcon(directory + "/up_arrow.png")
		self.plot_icon = QtGui.QIcon(directory + "/plot.png")

	def __init_plot_options(self):
		self.replace = QtGui.QRadioButton("Replace")
		self.include = QtGui.QRadioButton("Include")
		self.include.setChecked(True)
	
	#def __init_3d_options(self):
		#self.threed_options_label = QtGui.QLabel("3D options:",self)
		#self.replace_3d = QtGui.QCheckBox("Replace")
		#self.include_3d = QtGui.QCheckBox("Include")
		#self.include_3d.setChecked(True)

	def __init_preview_options(self):
		self.preview_options = QtGui.QComboBox(self)
		self.preview_options.addItem("No preview")
		self.preview_options.addItem("Single preview")
		self.preview_options.addItem("Multi preview")
		self.preview_options.setCurrentIndex(0)
		
		QtCore.QObject.connect(self.preview_options, QtCore.SIGNAL("currentIndexChanged(QString)"), self.preview_options_changed)
	
	def preview_options_changed(self,qstring):
		if str(qstring) == "Single preview":
			self.groupbox.setEnabled(True)
		else:
			self.groupbox.setEnabled(False)
			
		if str(qstring) == "No preview":
			self.force_2d.setEnabled(False)
			self.force_plot.setEnabled(False)
		else:
			self.force_2d.setEnabled(True)
			self.force_plot.setEnabled(True)
	
	def previews_allowed(self):
		if self.save_as_mode: return False
		
		return str(self.preview_options.currentText()) != "No preview"
	
	def single_preview_only(self):
		return str(self.preview_options.currentText()) == "Single preview"
	
	def __init_buttons(self):
		self.ok_button = QtGui.QPushButton("Ok")
		self.ok_button.adjustSize()
		
		self.cancel_button = QtGui.QPushButton("Cancel")
		self.cancel_button.adjustSize()
	
		QtCore.QObject.connect(self.ok_button, QtCore.SIGNAL("clicked(bool)"),self.ok_button_clicked)
		QtCore.QObject.connect(self.cancel_button, QtCore.SIGNAL("clicked(bool)"),self.cancel_button_clicked)
	
	def __init__force_2d_tb(self):
		self.force_2d = QtGui.QCheckBox("Force 2D")
		self.force_2d.setChecked(False)
		self.force_2d.setEnabled(False)
		
		QtCore.QObject.connect(self.force_2d, QtCore.SIGNAL("clicked(bool)"),self.force_2d_clicked)
		
	def __init__force_plot_tb(self):
		self.force_plot = QtGui.QCheckBox("Force plot")
		self.force_plot.setChecked(False)
		self.force_plot.setEnabled(False)
		
		QtCore.QObject.connect(self.force_plot, QtCore.SIGNAL("clicked(bool)"),self.force_plot_clicked)
	
	def force_2d_clicked(self):
		self.force_clicked(self.force_2d)
		
	def force_plot_clicked(self):
		self.force_clicked(self.force_plot)
	
	def force_clicked(self,f):
		if self.current_force == None:
			self.current_force = f
		elif f == self.current_force:
			self.current_force = None
			return
		else:
			self.current_force.setChecked(False)
			self.current_force = f

	def single_preview_clicked(self,bool):
		pass
		#print "not supported"
	
	def get_current_directory(self):
		directory = self.starting_directory 
		for i in range(len(self.list_widgets)-1):
			items = self.list_widgets[i].selectedItems()
			if len(items) >  0:
				directory += "/" + items[0].text()
#			elif len(items) > 1:
#				raise
			else:
				break
		
		directory = str(directory)
		if os.path.isfile(directory):
			[directory,unused] = os.path.split(directory)
		
		if not os.path.isdir(directory): # try for bdb format
			directory = directory.replace("bdb","EMAN2DB")
			if not os.path.isdir(directory): # in this case it's probably a bdb file minus the ".bdb"
				[directory,unused] = os.path.split(directory)
			if not os.path.isdir(directory):
				return None
		
		if len(directory) == 0: directory = "/"
		if directory[-1] != "/": directory += "/"
		
		return directory
	
	def cancel_button_clicked(self,bool):
		if self.save_as_mode:
			self.accept()
		else:
			self.module().emit(QtCore.SIGNAL("cancel"),self.selections)
	
	def ok_button_clicked(self,bool):
		if self.save_as_mode:
			directory = self.get_current_directory()
			if directory == None:
				msg = QtGui.QMessageBox()
				msg.setText("Can not deduce the current directory. Please update your selection")
				msg.exec_()
				return
			names = str(self.save_as_line_edit.text()).split()
			names = [name.strip(";") for name in names]
			
			if len(names)== 1:
				file = self.__convert_name_to_write_image_format(names[0],directory)
#				file = directory + names[0]
			else:
				file = [self.__convert_name_to_write_image_format(name,directory) for name in names]
			
			if self.validator == None:
				self.dialog_result = file
				self.accept()
			else:
				if not self.validator.validate_file_name(file): return
				else: 
					self.dialog_result = file
					self.accept()
			
		else:
			self.module().emit(QtCore.SIGNAL("ok"),self.selections)
			
	def __convert_name_to_write_image_format(self,name,directory):
		if len(name) > 3 and name[0:4] == "bdb:":
			if len(directory) > 0:
				last_bit = name[4:]
				dir = directory.rstrip("EMAN2DB/")
				dir = dir.rstrip("EMAN2DB")
				ret = "bdb:"+dir+"#"+last_bit
		elif directory.find("EMAN2DB/") != -1 or directory.find("EMAN2DB") != -1: # this test should be sufficient for establishing that bdb is the desired format
			ret = db_convert_path(directory+name)
		else: ret = directory + name
			
		return ret
	
	def __init_filter_combo(self):
		self.filter_text = QtGui.QLabel("Filter:",self)
		self.filter_combo = QtGui.QComboBox(None)
		self.filter_combo.addItem("EM types")
		self.filter_combo.addItem("Databases") # this doesn't really do anything
		self.filter_combo.addItem("*.spi,*.hdf,*.img, bdb:")
		self.filter_combo.addItem("*.hdf,*.img,*.hed,*.spi,bdb:,*.tif,*.mrc,*.dm3, *.pif")
		self.filter_combo.addItem("*.*")
		self.filter_combo.addItem("*")
		self.filter_combo.setEditable(True)
	
		QtCore.QObject.connect(self.filter_combo, QtCore.SIGNAL("currentIndexChanged(int)"),self.filter_index_changed)
#		QtCore.QObject.connect(self.filter_combo, QtCore.SIGNAL("currentIndexChanged(QString&)"),self.filter_index_changed)

	def filter_index_changed(self):
		self.__redo_list_widget_contents()
	
	def __redo_list_widget_contents(self):
		self.lock = True
		dtag = get_dtag()
		
		directory = self.starting_directory+dtag
		for i,data in  enumerate(self.list_widget_data):
			
			if data != None:d = str(data.text())
			old_row = self.list_widgets[i].currentRow()
			self.__load_directory_data(directory,self.list_widgets[i])
			self.list_widget_data[i] = self.list_widgets[i].item(old_row)
			if data == None: return
			else:
				directory += dtag + d
	
		self.lock = False
		
	def list_widget_context_menu_event(self,event):
		focus = None
		for l in self.list_widgets:
			if l.hasFocus(): 
				focus = l
				break
		else:
			raise
			return # No list widget has the focus even though its contextMenuEvent was triggered.
		event.accept()
		selected_items = []
		for i in range(l.count()):
			item = l.item(i)
			if item.isSelected(): selected_items.append(item)
			
		if len(selected_items) == 0: return
		
		# just for debu
		
		first_item = selected_items[0]
		
		md = first_item.get_metadata()
		nx,ny,nz=1,1,1
		if md != None:
			nx = md["nx"]
			ny = md["ny"]
			nz = md["nz"]

		multi_images_all_same_dims = False
		if md != None and len(selected_items) > 1:
			req_keys = ["nx","ny","nz"]
			gtg = True # good to go
			for k in req_keys:
				if not md.has_key(k): gtg = False
				break
			
			if gtg:
				multi_images_all_same_dims = True
				# this means I think it's an emdata. This is potentially a fatal assumption, but unlikely for the time being
				# now check to make sure all other selected items have the samve value
				
				for i in range(1,len(selected_items)):
					mdt = selected_items[i].get_metadata()
					gtgt = True # good to go tmp
					for k in req_keys:
						if not mdt.has_key(k): gtg = False
						break
					
					if not gtgt or mdt["nx"] != nx or mdt["ny"] != ny or mdt["nz"] != nz:
						multi_images_all_same_dims = False
						break

		options = first_item.context_menu_options
		if len(options) == 0: return
		options_keys = options.keys()
		
		# Make sure the options are applicable to all selected items
		for i in range(1,len(selected_items)):
			o_k = selected_items[i].context_menu_options.keys()
			rm = []
			for k in options_keys:
				if k not in o_k: rm.append(k)
			for r in rm: options_keys.remove(r)
			
			if len(options_keys) == 0: return
		
		menu = QtGui.QMenu()
		self.menu_selected_items = selected_items
		for k in options_keys:
#			action = QtGui.QAction(k,menu)
#			action.items = selected_items
			menu.addAction(k)
		
		if multi_images_all_same_dims:
			# These are customized actions that somewhat break the modularity, but I don't think it's too bad
			menu.addSeparator()
			menu.addAction("Save As Stack")
			if nz == 1: 
				menu.addAction("Preview Subset")
		
		if nz == 1 and ny > 1:	menu.addAction("Open in e2boxer")
			
		QtCore.QObject.connect(menu,QtCore.SIGNAL("triggered(QAction*)"),self.menu_action_triggered)
		self.action_list_widget = l # only set if the menu acutally triggers
		menu.exec_(event.globalPos())
	
	def menu_action_triggered(self,action):
		items = self.menu_selected_items
		
		cont = self.__check_action(action.text(),items) # this kind of breaks the OO design, but it's necessary in the current framework
		if not cont: return
		total = len(items)
		
		if action.text() == "Save As Stack":
			#data = []
			#for item in self.menu_selected_items:
				#data.append(item)
				
			save_data(self.menu_selected_items)
		elif action.text() == "Preview Subset":
			data = []
			for item in self.menu_selected_items:
				data.append(item.get_emdata())
			self.preview_data(data,"")
		elif action.text() == "Open in e2boxer":
			names = []
			for item in self.menu_selected_items:
				names.append(item.get_path())
			from e2boxer import EMBoxerModule, EmptyObject
			options = EmptyObject()
			options.filenames = names
			cwd = e2getcwd()
			options.running_mode = "gui"
			options.method = "Swarm"
			options.boxsize=128
			if db_check_dict("bdb:e2boxer.project"):
				db = db_open_dict("bdb:e2boxer.project")
				if db.has_key("working_boxsize"):
					options.boxsize = db["working_boxsize"] 
			
			self.boxer_module = EMBoxerModule(get_application(),options)
			self.boxer_module.show_guis()
				
		else:		
			items_acted_on = []
			for i,item in enumerate(items):
				if not item.context_menu_options[str(action.text())](self):
					# this is fine if the user hit cancel when they were saving a whole bunch of images
					# but it's not so fine if an error was thrown while deleting... hmmm needs some thought
					break
				else:
					items_acted_on.append(item)

		self.action_list_widget = None
		self.menu_selected_items = None
		
	def __post_action(self,action_str,items_acted_on):
		if action_str == "Delete":
			c_list_widget = self.action_list_widget
			items = [c_list_widget.item(i) for i in range(c_list_widget.count())]
			rm_indices = []
			for j,i in enumerate(items):
				if i in items_acted_on:
					rm_indices.append(j)
			
			rm_indices.reverse()
			for idx in rm_indices:
				c_list_widget.takeItem(idx)
			
			for ii,l in enumerate(self.list_widgets):
				if l == c_list_widget:
					break
			
			ii += 1
			if ii <= (len(self.list_widgets) -1):
				for i in range(ii,len(self.list_widgets)):
					self.list_widgets[i].clear()
			
	def __check_action(self,action_str,items):
		if action_str == "Delete":
			msg = QtGui.QMessageBox()
			msg.setText("Deletion will be permanent. Are you sure you want to delete the selected file(s)?")
			s = ""
			for i in items: s+=i.text()+"\n"
			msg.setInformativeText(s)
			msg.setStandardButtons(QtGui.QMessageBox.Cancel | QtGui.QMessageBox.Ok )
			msg.setDefaultButton(QtGui.QMessageBox.Cancel)
			ret = msg.exec_()
			if ret == QtGui.QMessageBox.Cancel: return False
			else: return True
		else:
			return True
		
	def __add_list_widget(self, list_widget = None):
		if list_widget == None:	list_widget = QtGui.QListWidget(None)
		
		list_widget.contextMenuEvent = self.list_widget_context_menu_event
		
		if self.single_selection:list_widget.setSelectionMode(QtGui.QAbstractItemView.SingleSelection)
		else: list_widget.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
		list_widget.setMouseTracking(True)	
		self.list_widgets.append(list_widget)
		self.splitter.addWidget(list_widget)
		self.list_widget_data.append(None)
		
		QtCore.QObject.connect(list_widget, QtCore.SIGNAL("itemDoubleClicked(QListWidgetItem*)"),self.list_widget_dclicked)
		QtCore.QObject.connect(list_widget, QtCore.SIGNAL("itemPressed(QListWidgetItem*)"),self.list_widget_clicked)
		QtCore.QObject.connect(list_widget, QtCore.SIGNAL("currentRowChanged (int)"),self.list_widget_row_changed)
		#QtCore.QObject.connect(list_widget, QtCore.SIGNAL("paintEvent (int)"),self.list_widget_row_changed)
		QtCore.QObject.connect(list_widget, QtCore.SIGNAL("itemEntered(QListWidgetItem*)"),self.list_widget_item_entered)
		#QtCore.QObject.connect(list_widget, QtCore.SIGNAL("currentItemChanged(QListWidgetItem*,QListWidgetItem*)"),self.list_widget_current_changed)
		#QtCore.QObject.connect(list_widget, QtCore.SIGNAL("itemChanged(QListWidgetItem*)"),self.list_widget_item_changed)
		#\QtCore.QObject.connect(list_widget, QtCore.SIGNAL("itemActivated(QListWidgetItem*)"),self.list_widget_item_activated)
		#QtCore.QObject.connect(list_widget, QtCore.SIGNAL("activated(QModelIndex)"),self.activated)
		QtCore.QObject.connect(list_widget, QtCore.SIGNAL("itemSelectionChanged()"),self.selection_changed)
	
	def __go_back_a_directory(self):
		self.lock = True
		dtag = get_dtag()
		
		new_dir = self.starting_directory[0:self.starting_directory.rfind(dtag)]
#		new_dir.replace("/bdb/","EMAN2DB/") # this works for the time being but needs investigating
#		if len(new_dir) > 4 and new_dir[-4:] == "/bdb": 
#			new_dir.replace("/bdb","/EMAN2DB")  # this works for the time being but needs investigating
#			print "replaced",new_dir
#		else:
#			print "no replace"
		if len(new_dir) == 0: new_dir = get_dtag()
		elif  new_dir[-1] == ":": new_dir += get_dtag() # C: becomes C:/
		elif len(new_dir) == 1 and new_dir != get_dtag(): new_dir += ":/"# C becomes C:/
		#if  self.db_listing.responsible_for(new_dir):
			#new_dir = self.db_listing.
		#print "going back a directory",new_dir
#		if not os.access(new_dir,os.R_OK):
#			print new_dir
#			print "can't go up a directory, don't have read permission"
#			self.new_dir = e2gethome() # send the user back to home
			
		self.starting_directory = new_dir
		for j in range(0,len(self.list_widgets)):
			widget = self.list_widgets[j]
			widget.clear()
			self.list_widget_data[j] = None
			# get rid of these attributes, they could have been set by the listing objects
			del_attr = ["directory","directory_lister","directory_arg","mod_time"]
			for attr in del_attr:
				if hasattr(widget,attr): delattr(widget,attr)
				
		self.__load_directory_data(self.starting_directory,self.first_list_widget)
		#self.hide_preview()
		self.lock = False
		
	def __go_forward_a_directory(self):
		self.lock = True
		self.starting_directory = self.starting_directory + get_dtag() + str(self.list_widget_data[0].text())
		dtag = get_dtag()
		directory = self.starting_directory 
		for i in range(len(self.list_widgets)-1):
			items = []
			li = self.list_widgets[i]
			lip = self.list_widgets[i+1]
			old_row = lip.currentRow()
			n = lip.count()
			for j in range(n-1,-1,-1):
				items.append(lip.takeItem(j))
				
			li.clear()	
			for k in items:
				li.insertItem(0,k)
			
			li.setCurrentRow(old_row)
			req_attr = ["mod_time","directory","directory_lister","directory_arg"]
			if  hasattr(self.list_widgets[i+1],req_attr[0]): # if it has the first it has them all
				for attr in req_attr: setattr(li,attr,getattr(lip,attr))
	   	   	else:
	   	   		# just make sure the attributes don't exist
	   	   		for attr in req_attr:
	   	   			if hasattr(li,attr): delattr(li,attr)
				
			
			self.list_widget_data[i] = li.item(old_row)
			directory += dtag + str(self.list_widget_data[i].text())
		
		
		a = QtGui.QListWidgetItem(self.up_arrow_icon,"../",None)
		a.type_of_me = "go up a directory"
		self.list_widgets[0].insertItem(0,a)
		
		self.lock = False
		self.hide_preview()
		
	def selection_changed(self):
		if not self.lock:
			self.__update_selections()
		
	def __update_selections(self):
		'''
		Makes the list of currently selected files accurate and up to date. Called when
		something has been clicked in a a list widget
		'''
		
		#print self.current_list_widget
		# get the directory 
		dtag = get_dtag()
		directory = self.starting_directory+dtag
		idx = 0
		for i,list_widget in enumerate(self.list_widgets):
			if list_widget == self.current_list_widget: 
				break
			if (self.list_widget_data[i] == None):
				# this means it could be a header attribute, for example
				return
			directory += str(self.list_widget_data[i].text()) + dtag
		else:
			print "no list widget has focus?"
			return
		
		# now make the list of selections reflect what is need to display them using
		
		self.selections = []
		items = self.current_list_widget.selectedItems()
		if len(items) > 0:
			a = items[0]
			previewable,previewer= self.__is_previewable(a)
			
			self.selections = previewer.paths(items)
			
		# Make sure the save as line is updated if we're in that mode
		if self.save_as_mode:
			text = ""
			for i,name in enumerate(self.selections):
				if i > 0:
					text+= "; "
				text += base_name(name)
				
			self.save_as_line_edit.setText(text)

	def hide_preview(self):
		if self.save_as_mode: return
		
		if self.gl_image_preview  != None:
			get_application().hide_specific(self.gl_image_preview)
	
	def list_widget_item_entered(self,item):
		self.current_list_widget = item.listWidget()

	def list_widget_row_changed(self,i):
		return
	
	def keyPressEvent(self,event):
		if event.key() == Qt.Key_F1:
			try:
				import webbrowser
				webbrowser.open("http://blake.bcm.edu/emanwiki/e2display")
			except:
				pass
	
	def list_widget_clicked(self,item):
		self.list_widget_clickerooni(item,False)
	
	def list_widget_dclicked(self,item):
		self.list_widget_clickerooni(item,True)
	
	def list_widget_clickerooni(self,item,allow_preview=True):
		if self.lock : return
		if self.current_list_widget == None: return
		if item.type_of_me == "value": return #it's just a value in the db
		self.__update_selections()
		
		if item == None: return		
		
		if item.text() == "../": 
			self.__go_back_a_directory()
			return
		
		file = self.starting_directory+"/"
		directory = self.starting_directory+"/"
		#idx = 0
		for i,list_widget in enumerate(self.list_widgets):
			if list_widget == self.current_list_widget: 
				idx = i
				file += str(item.text())
				break
			file += str(self.list_widget_data[i].text()) + "/"
			directory += str(self.list_widget_data[i].text()) + "/"
		else:
			print "no list widget has focus?"
			return
		
		if allow_preview and self.previews_allowed() and self.try_preview_item(item):
			
			for i in range(idx+1,len(self.list_widgets)):
				self.list_widgets[i].clear()
				self.list_widget_data[i] = None
				
			if not self.check_preview_item_wants_to_list(item):
				print "error 101" # facetiousness?
				return
		
	
		n = len(self.list_widgets)-1
		if self.current_list_widget  == self.list_widgets[n]:
				self.lock = True
				self.list_widget_data[n] = item
				self.__go_forward_a_directory()
				self.__load_directory_data(file+"/",self.list_widgets[n],item)
				self.lock = False
				return

		if self.__load_directory_data(file,self.list_widgets[idx+1],item):
			#if old_item != None:
				#old_item.setBackgroundColor(QtGui.QColor(255,255,255))
			##item.setBackgroundColor(QtGui.QColor(64,190,0,63))	
			self.list_widget_data[idx] = item
			self.list_widget_data[idx+1] = None
			
		for i in range(idx+2,len(self.list_widgets)):
			self.list_widgets[i].clear()
			self.list_widget_data[i] = None
	
	def try_preview_item(self,item):
		get_application().setOverrideCursor(Qt.BusyCursor)
		preview_occured = item.do_preview(self)
		get_application().setOverrideCursor(Qt.ArrowCursor)
		return preview_occured

	def check_preview_item_wants_to_list(self,item):
		if self.db_listing.preview_item_wants_to_list(item):
			return True
		elif self.dir_listing.preview_item_wants_to_list(item):
			return True
		
		return False
	
	def preview_plot(self,filename):
		if self.save_as_mode: return
		
		if self.single_preview_only():
			if not isinstance(self.gl_image_preview,EMPlot2DModule):
				if self.gl_image_preview != None: get_application().close_specific(self.gl_image_preview)
				self.gl_image_preview = EMPlot2DModule(get_application())
	
			self.gl_image_preview.set_data_from_file(filename,self.replace.isChecked())
			get_application().show_specific(self.gl_image_preview)
			self.gl_image_preview.updateGL()
			
		else:
			preview = EMPlot2DModule(get_application())
			preview.set_data_from_file(filename,self.replace.isChecked())
			get_application().show_specific(preview)
			
	def preview_plot_list(self,title,list_data):
		if self.save_as_mode: return
		
		if self.single_preview_only():
			if not isinstance(self.gl_image_preview,EMPlot2DModule):
				if self.gl_image_preview != None: get_application().close_specific(self.gl_image_preview)
				self.gl_image_preview = EMPlot2DModule(get_application())
	
			self.gl_image_preview.set_data(title,list_data,self.replace.isChecked())
			get_application().show_specific(self.gl_image_preview)
			self.gl_image_preview.updateGL()
			
		else:
			preview =EMPlot2DModule(get_application())
			preview.set_data_from_file(filename,self.replace.isChecked())
			get_application().show_specific(preview)
			
	def preview_data(self,a,filename=""):
		if self.save_as_mode: return
		
		from emimage import EMImageModule, EMModuleFromFile
		
		using_file_names_only = False
		if a == None: using_file_names_only = True # For the image matrix, you can load large image stacks if you specify only the file name
		
		f_2d = self.force_2d.isChecked()
		f_plot = self.force_plot.isChecked()
		
		if self.single_preview_only() and len(self.previews) != 0:
			old_preview = self.previews[-1] # this means we always choose the last preview if the user suddenly goes from multipreview to single preview
			if not using_file_names_only:
				preview = EMImageModule(data=a,app=get_application(),force_2d=f_2d,force_plot=f_plot,old=old_preview,filename=filename,replace=self.replace.isChecked())
			else:
				preview = EMModuleFromFile(filename,application=get_application(), force_2d=f_2d,force_plot=f_plot,old=old_preview)
			if preview != old_preview:
				self.module_closed(old_preview)
				old_preview.closeEvent(None)
				self.previews.append(preview)
				self.module_events.append(ModuleEventsManager(self,preview))
				try: preview.optimally_resize()
				except: pass
		else:
			if not using_file_names_only:
				preview = EMImageModule(data=a,app=get_application(),force_2d=f_2d,force_plot=f_plot,filename=filename)
			else:
				preview = EMModuleFromFile(filename,application=get_application(), force_2d=f_2d,force_plot=f_plot)
			self.previews.append(preview)
			self.module_events.append(ModuleEventsManager(self,preview))
			try: preview.optimally_resize()
			except: pass
					
		get_application().show_specific(preview)
				
		preview.updateGL()
			
		get_application().setOverrideCursor(Qt.ArrowCursor)	
	
	def module_closed(self,module):
		import sys
		for i,mod in enumerate(self.previews):
			if mod == module:
				p = self.previews.pop(i)
				self.module_events.pop(i)
				# right here is where the memory should be cleaned up for the module, so you could put some print statements like this if you were memory debugging:
				# print sys.getrefcount(p)
				# To date I have made sure the modules are being deleted
				return
			
		print "failed to close module?" # this shouldn't happen if I have managed everything correctly
	
	def closeEvent(self,event):
		self.module_events=[]
		for mod in self.previews:
			mod.closeEvent(None)
		
		self.module().closeEvent(event)
	
	def get_file_filter(self):
		return str(self.filter_combo.currentText())
	
	def filter_strings(self,strings):
		
		filters = str(self.filter_combo.currentText()).split(",")
		
		for j,f in enumerate(filters):
			s = f.replace("*","\w*")
			s = s.replace(".","\.")
			filters[j] = s
		
		reg_exp = []
		for f in filters:
			reg_exp.append(re.compile(f))
		
		solution = []
		for s in strings:
			for r in reg_exp:
				if len(re.findall(r,s)) != 0:
					solution.append(s)
					break
					
		
		return solution
	
	def __is_previewable(self,item):
		#if  self.db_listing.responsible_for(file_name):
		if self.db_listing.is_previewable(item): return True,self.db_listing
		else: return self.dir_listing.is_previewable(item),self.dir_listing
		
	def __is_non_empty_directory(self,s):
		'''
		Returns true if s is a non empty directory
		'''

		dirs, files  = get_files_and_directories(s)
		files = self.filter_strings(files)
		file_length = len(files)
		if file_length == 0: file_length = len(dirs)
	
		if file_length != 0: return True
		
		return False
	
	def __load_database_directory(self,database_name,list_widget):
		
		print "loading database directory",database_name
	
	
	def __refresh_directory_data(self):
		self.__load_directory_data(self.ldd_directory,self.ldd_list_widget,self.ldd_item)
	
	def __load_directory_data(self,directory,list_widget,item=None):
		self.ldd_directory = directory
		self.ldd_list_widget = list_widget
		self.ldd_item = item
		
		list_widget.clear()
		if (list_widget == self.list_widgets[0]):
			self.lock = True
			a = QtGui.QListWidgetItem(self.up_arrow_icon,"../",list_widget)
			a.type_of_me = "go up a directory"
			self.lock = False
			
		if item != None and item.type_of_me == "key_value":
			a = QtGui.QListWidgetItem(self.basic_python_icon,str(item.value),list_widget)
			a.type_of_me = "value"
			a.value = item.value
			return
		
		if self.db_listing.responsible_for(directory):
			if  self.db_listing.load_database_data(directory,list_widget): return True
			elif item != None and self.db_listing.load_metadata(item,list_widget): return True
			elif self.db_listing.load_directory_data(directory,list_widget): return True
			elif self.db_listing.load_database_interrogation(directory,list_widget): return True
			
			else: return False
		else: 
			if item != None and self.dir_listing.load_metadata(item,list_widget): return True
			else: return self.dir_listing.load_directory_data(directory,list_widget)
		
	def make_replacements(self,dirs,list_widget):
		self.db_listing.make_replacements(dirs,list_widget)



class EMFSListing:
	'''
	FS stands for File System
	'''
	def __init__(self,target):
		self.target = weakref.ref(target)
		
		self.emdata_mx = "directory_emdata_mx"
		self.emdata_mx_3d = "directory_emdata_mx_3d"
		self.emdata_3d = "directory_emdata_3d"
		self.emdata = "directory_emdata"
		self.emdata_mx_member = "directory_emdata_mx_member" # for individual images in mxs
		self.plot_data = "plot_data"
 		
		self.previewable_types = [self.emdata_mx,self.emdata_3d,self.emdata,self.emdata_mx_member,self.plot_data, self.emdata_mx_3d] # add any others as functionality grows
		self.threed_dim_limit = 128
		pass
		
	def paths(self,items):
		ret = []
		for item in items:
			if self.is_previewable(item):
			 	ret.append(item.full_path)
		return ret

	def filter_strings(self,strings):
		
		filt = self.target().get_file_filter()
		if filt == "EM types": 	
			return [i for i in strings if i[-4:]!=".hed"]
#			return strings # this is a bit of hack unfortunately
		
		filters = filt.split(",")

		for j,f in enumerate(filters):
			s = f.replace("*","\w*")
			s = s.replace(".","\.")
			filters[j] = s
		
		reg_exp = []
		for f in filters:
			reg_exp.append(re.compile(f))
		
		solution = []
		for s in strings:
			for r in reg_exp:
				if len(re.findall(r,s)) != 0:
					solution.append(s)
					break
					
		
		return solution
	
	def load_directory_data(self,directory,list_widget):
		e = EMData()
		dtag = get_dtag()
		read_header_only = True
		
		dirs, files = get_files_and_directories(directory)
		if len(files) == 0 and len(dirs) == 0:
			# something is unusual about the directory, for instance, it is a file
			return
		files = self.filter_strings(files)
		filt = self.target().get_file_filter()
		
		statinfo = os.stat(directory)
		list_widget.mod_time = statinfo[-2] # this is to record when the last modification time was, so automatic updates are possible
		list_widget.directory = directory
		list_widget.directory_lister = self
		list_widget.directory_arg = directory
		
		dirs.sort()
		files.sort()
		 
		self.target().make_replacements(dirs,list_widget)
		 
		for i in dirs:
			if i[0] == '.': continue
			
			file_length = 0
	   	   	d, f = [], []
			f = self.filter_strings(f)
			file_length = len(f)
			if file_length == 0: file_length = len(d)
			if file_length != 0:
				a = QtGui.QListWidgetItem(self.target().folder_files_icon,i,list_widget)
				b = EMGenericItem("directory_file",i)
				accrue_public_attributes(a,b)
			else:
				a = QtGui.QListWidgetItem(self.target().folder_icon,i,list_widget)
				b = EMGenericItem("directory",i)
				accrue_public_attributes(a,b)
					
		for file in files:
			if file[0] == '.': continue
			if file[-1] == '~': continue
			#print EMUtil.get_image_ext_type(Util.get_filename_ext((file)),get_file_tag(file)
			extension = Util.get_filename_ext(file)
			full_name = directory+dtag+file
			# note, if this if statement is allowed to proceed on Windows in the case of a png then the program
			# crashes. In December of 2008 I thus changed this if statement to automatically exclude unecessary files
			# such as pngs and jpges...etc.
			if EMUtil.get_image_ext_type(extension) != IMAGE_UNKNOWN and extension not in ["png","jpeg","jpg","JPG"]:
				try:
					e.read_image(full_name,0, read_header_only)
					if EMUtil.get_image_count(full_name) > 1:
						if e.get_zsize() > 1:
							a = QtGui.QListWidgetItem(self.target().emdata_matrix_3d_icon,file,list_widget)
							b = EMFS3DImageStackItem(self.emdata_mx_3d,full_name)
							accrue_public_attributes(a,b)
						else:
							a = QtGui.QListWidgetItem(self.target().emdata_matrix_icon,file,list_widget)
							b = EMFS2DImageStackItem(self.emdata_mx,full_name)
							accrue_public_attributes(a,b)
					else:
						if e.get_zsize() > 1:
							a = QtGui.QListWidgetItem(self.target().emdata_3d_icon,file,list_widget)
							b = EMFSSingleImageItem(self.emdata_3d,full_name)
						else:
							a = QtGui.QListWidgetItem(self.target().emdata_icon,file,list_widget)
							b = EMFSSingleImageItem(self.emdata_3d,full_name)
	   	   	   	   	   	accrue_public_attributes(a,b)
				except:
					a = QtGui.QListWidgetItem(self.target().file_icon,file,list_widget) # this happens when files are corrupted	
					b = EMGenericItem("regular_file",file)
					accrue_public_attributes(a,b)
			
			elif EMPlot2DModule.is_file_readable(full_name):
				a = QtGui.QListWidgetItem(self.target().plot_icon,file,list_widget)
				b = EMFSPlotItem(self.plot_data,full_name)
				accrue_public_attributes(a,b)
			else:
				if filt != "EM types":
					a = QtGui.QListWidgetItem(self.target().file_icon,file,list_widget)
					b = EMGenericItem("regular_file",file)
					accrue_public_attributes(a,b)

		return True
			
	def do_preview(self,item):
		return item.do_preview(self.target())
	
	def load_metadata(self,item,list_widget):
		if item.type_of_me in [self.emdata_3d, self.emdata,self.emdata_mx_member]:
			data = self.get_emdata_header(item)
			keys = data.keys()
			keys.sort() # alphabetical order
			for k in keys:
				v = data[k]
				a = QtGui.QListWidgetItem(self.target().basic_python_icon,str(k)+":"+str(v),list_widget)
				b = EMKeyValueItem("key_value",k,v)
				accrue_public_attributes(a,b)

			return True
		elif item.type_of_me == self.emdata_mx:
			for i in range(EMUtil.get_image_count(item.full_path)):
				a = QtGui.QListWidgetItem(self.target().emdata_icon,str(i),list_widget)
				b = EMFSImageStackMemberItem(self.emdata_mx_member,item.full_path,i)
				accrue_public_attributes(a,b)
				
			return True
		
		elif item.type_of_me == self.emdata_mx_3d:
			for i in range(EMUtil.get_image_count(item.full_path)):
				a = QtGui.QListWidgetItem(self.target().emdata_3d_icon,str(i),list_widget)
				b = EMFSImageStackMemberItem(self.emdata_mx_member,item.full_path,i)
				accrue_public_attributes(a,b)
				
			return True
				
				
		return False
	
	def get_emdata_header(self,item):
		return item.get_metadata()
			
	def is_previewable(self,item):
		return item.type_of_me in self.previewable_types
		
	def preview_item_wants_to_list(self,item):
		return item.is_previewable()

	#def is_previewable(self,file_or_folder):
		## this may not be strictly correct, seeing as if it's a file it will return true
		#return os.path.isfile(file_or_folder)

class EMListingItem:
	'''
	Base class definition providing the pubic interface of list widget items as 
	required by the EMSelectorDialog, public interface required by the EMListing objects (DB and FS)
	'''
	def __init__(self):
		self.context_menu_options = {} # this is used for running context menu actions
		pass
	def do_preview(self,target):
		'''
		A function that is supposed to call the targets preview_data or preview_plot function,
		then return True
		'''
		return False
	
	def is_previewable(self):
		'''
		Essentially asking if the subclassing object provides its own definition
		of the do_preview functoin
		''' 
		return False

	def get_metadata(self):
		'''
		Returns a dictionary if defined
		'''
		return None
	
	def get_attr_dict(self):
		'''
		A wrapper for get_metadata - originally added to allow
		a generic interface in EMStackSaveDialog
		'''
		return self.get_metadata()
	
	def get_path(self):
		'''
		Returns the file name path
		DB items have to use a special syntax
		'''
		return None
	
	def get_emdata(self):
		'''
		supposed to return an EMData object, if possible
		'''
		return None

class EMFSListingItem(EMListingItem):
	def __init__(self,type,full_name):
		EMListingItem.__init__(self)
		self.type_of_me = type
		self.full_path = full_name
		self.context_menu_options["Delete"] = self.delete
		self.context_menu_options["Save As"] = self.save_as
		
	def delete(self,target):
		try:
			remove_file(self.full_path)
			return True
		except: return False
		
	def save_as(self,target):
		save_data(self.get_emdata())
		#return save_as(self.full_path,target.application())
	
	def get_path(self):
		return self.full_path
	
	def get_emdata(self):
		e = EMData()
		e.read_image(self.full_path)
		return e
	
class EMFS2DImageStackItem(EMFSListingItem):
	def __init__(self,type,full_name):
		EMFSListingItem.__init__(self,type,full_name)
		
	def do_preview(self,target):
		target.preview_data(None,self.full_path)
		return True
		
	def is_previewable(self): return True
	
	def get_emdata(self):
		'''This one returns a list'''
		from emimagemx import EMDataListCache
		return EMDataListCache(self.full_path)
	
class EMFS3DImageStackItem(EMFSListingItem):
	def __init__(self,type,full_name):
		EMFSListingItem.__init__(self,type,full_name)

	# no preview for this item as of Feb 2009
	def get_emdata(self):
		'''This one returns a list'''
		e = EMData().read_images(self.full_path)
		return e
	
class EMFSImageStackMemberItem(EMFSListingItem):
	def __init__(self,type,full_name,idx):
		EMFSListingItem.__init__(self,type,full_name)
		self.idx = idx
		
	def do_preview(self,target):
		a=EMData(self.full_path,self.idx)
		target.preview_data(a,self.full_path)
		return True
		
	def is_previewable(self): return True

	def get_metadata(self):
		e = EMData()
		e.read_image(self.full_path,self.idx, read_header_only)
		return e.get_attr_dict()
	
	def get_emdata(self):
		e = EMData()
		e.read_image(self.full_path,self.idx)
		return e
	
class EMFSSingleImageItem(EMFSListingItem):
	def __init__(self,type,full_name):
		EMFSListingItem.__init__(self,type,full_name)

	def do_preview(self,target): 
		target.preview_data(EMData(self.full_path),self.full_path)
		return True
		
	def is_previewable(self): return True

	def get_metadata(self):
		e = EMData()
		e.read_image(self.full_path,0, read_header_only)
		return e.get_attr_dict()
	
	def get_emdata(self):
		e = EMData()
		e.read_image(self.full_path,0)
		return e

class EMFSPlotItem(EMFSListingItem):
	def __init__(self,type,full_name):
		EMFSListingItem.__init__(self,type,full_name)
		
	def do_preview(self,target):
		target.preview_plot(self.full_path)
		return True
		
	def is_previewable(self): return True
	
class EMDBListing:
	def __init__(self,target):
		self.target = weakref.ref(target)
		self.directory_replacements = {"EMAN2DB":"bdb"}
	
		self.db_mx = "database_emdata_mx"

		self.emdata_3d_entry = "database_emdata_3d_entry"
		self.emdata_entry = "database_emdata_entry"
		self.db_dict_emdata_entry = "database_dictionary_emdata"
		self.db_list_plot = "list_plot"
		
		self.previewable_types = [self.db_mx,self.emdata_3d_entry,self.emdata_entry,self.db_dict_emdata_entry,self.db_list_plot] # add any others as functionality grows
	
	def paths(self,items):
		if len(items) > 1:
			# the point here is to see if the selected images are in the same bdb matrix directory
			# if so we can exploit the "?" terminology
			# Now the selector doesn't facilitate choosing images from seperate directories
			# so if the first two items are from the same matrix database, then they all are
			# if this wasn't the case then we could iterate through the items and group them, but
			# currently isn't so hence this simplified approach
			# FEB 2009 - I'm a bit confused about what I was doing, I added a try block because I was getting an error sometimes. (d.woolford)			
			item0 = items[0]
			item1 = items[1]
			try:
			# check to see if the first two items are in the same directory
				if (item0.database_directory+item0.database) == (item1.database_directory+item1.database):
					# we know assume they are all in the same directory
					return_val = "bdb:"+item0.database_directory+'#'+item0.database + "?"
					for i,item in enumerate(items):
						if i != 0:
							return_val += ","
						return_val += str(item.database_key)
						
					return [return_val]
			except: pass
					
		# if the return didn't happen then this generic loop will work -
		ret = []
		for item in items:
			if item.type_of_me == self.db_mx:
				ret.append("bdb:"+item.database_directory+'#'+item.database)
			elif item.type_of_me in [self.emdata_3d_entry,self.emdata_entry]:
				# there is one problem with this approach - if the user is in a matrix
				# directory then choosing the 0th image won't really work as was intended
				if item.database_key != 0:
					ret.append("bdb:"+item.database_directory+'#'+item.database+"?"+str(item.database_key))
				else:
					ret.append("bdb:"+item.database_directory+'#'+item.database)
			else:
				pass #it's not viewable
		return ret

	def responsible_for(self,file_or_folder):
		real_name = self.convert_to_absolute_path(file_or_folder)
		#print file_or_folder,real_name
		dtag = get_dtag()
		split = real_name.split(dtag)
		split.reverse() # this probably makes things faster
		for s in split:
			if s in self.directory_replacements.keys() or (len(s) > 4 and s[-4:] == ".bdb"): return True 
	
		return False
	
	def make_replacements(self,dirs,list_widget):
		rm = []
		for i,directory in enumerate(dirs):
			d = remove_directories_from_name(directory)
			
			if d in self.directory_replacements.keys():
				a = QtGui.QListWidgetItem(self.target().database_icon,self.directory_replacements[d],list_widget)
				rm.append(i)
				b = EMGenericItem("regular",None)
				accrue_public_attributes(a,b)
		
		rm.reverse()
		for r in rm:
			dirs.pop(r)
	
	def convert_to_absolute_path(self,file_or_folder):
		dtag = get_dtag()
		ret = file_or_folder
		found = False
		for dir_rep in self.directory_replacements.items():
			if ret.find(dtag+dir_rep[1]) != -1:
				ret = ret.replace(dtag+dir_rep[1],dtag+dir_rep[0])
				found = True
		if not found: return ret
		if (not os.path.isdir(ret)) and (not os.path.isfile(ret)):
			if ret[-1] == dtag: ret = ret[:-1]
			ret += ".bdb"
			
		return ret
			
	def is_database_directory(self,directory):
		if remove_directories_from_name(directory) in self.directory_replacements.values(): return True
		else: return False

	def load_directory_data(self,directory,list_widget):
		
		'''
		Displays the file/folder information in the directory /home/someonone/data/EMAN2DB
		this will typically consist of .bdb (database) files, but may contain folders and other
		EMAN2DB directories.
		
		At the moment I have only written the code so that it supports the interrogation of the .bdb
		files, and am displaying the other folders only as a I reminder that they need to be dealt with
		'''
		if not remove_directories_from_name(directory) in self.directory_replacements.values():
			 return False

		dtag = get_dtag()
		real_directory = self.convert_to_absolute_path(directory)
		dirs,files = get_files_and_directories(real_directory)
		files.sort()
		dirs.sort()
		
		statinfo = os.stat(real_directory)
		list_widget.mod_time = statinfo[-2] # this is to record when the last modification time was, so automatic updates are possible
		list_widget.directory = real_directory
		list_widget.directory_lister = self
		list_widget.directory_arg = directory
		
		
		for i in dirs:
			if i[0] == '.': continue
			
			if i == "EMAN2DB":
				a = QtGui.QListWidgetItem(self.target().database_icon,"bdb",list_widget) # this is something we do not wish to support
				b = EMGenericItem("unwanted","bdb")  # really haven't accommodated for this...
				accrue_public_attributes(a,b)
				continue
		
			file_length = 0
			d,f = [], []
			file_length = len(f)
			if file_length == 0: file_length = len(d)
				
			if file_length != 0:
				a = QtGui.QListWidgetItem(self.target().folder_files_icon,i,list_widget)
				b = EMGenericItem("directory_file",i)
				accrue_public_attributes(a,b)
			else:
				a = QtGui.QListWidgetItem(self.target().folder_icon,i,list_widget)
				b = EMDBDirectoryItem("directory",real_directory)
				accrue_public_attributes(a,b)
			
		for file in files:
			if file[len(file)-3:] == "bdb":
				f = file.rpartition(".bdb")
				db_directory = self.get_emdatabase_directory(real_directory)
				
				db_name = "bdb:"+db_directory+"#"+f[0]
				db = db_open_dict(db_name,ro=True)
				
				if db and db.has_key("maxrec"):
					#n = DB[f[0]]["maxrec"]
					n = db["maxrec"]
					if n >= 1:
						a = QtGui.QListWidgetItem(self.target().emdata_matrix_icon,f[0],list_widget)
						b = EMDB2DImageStackItem(self.db_mx,db_directory,f[0])
						accrue_public_attributes(a,b)
					elif n == 0:
						#d = DB[f[0]].get_header(0)
						d = db.get_header(0)
						if d["nz"] <= 1:
							a = QtGui.QListWidgetItem(self.target().emdata_icon,f[0],list_widget)
							type = self.emdata_entry
						else:
							a = QtGui.QListWidgetItem(self.target().emdata_3d_icon,f[0],list_widget)
							type = self.emdata_3d_entry
						
						b = EMDBSingleImageItem(type,db_directory,f[0],0)
						accrue_public_attributes(a,b)				
					else:
						a = QtGui.QListWidgetItem(self.target().database_icon,f[0],list_widget)
						b = EMGenericItem("regular",f[0])
						accrue_public_attributes(a,b)
				else:
					a = QtGui.QListWidgetItem(self.target().database_icon,f[0],list_widget)
					b = EMGenericItem("regular",f[0])
					accrue_public_attributes(a,b)
					
		return True

	def is_database_file(self,file_name):
		file = self.convert_to_absolute_path(file_name)
		if len(file) > 4 and file[-4:] == ".bdb":
			if self.get_last_directory(file) == "EMAN2DB":
				if file_exists(file):
					return True
			
		return False

	def load_database_data(self,directory,list_widget):
		
		if not self.is_database_file(directory): 
			return False
		
		file = self.convert_to_absolute_path(directory)
 
		db_directory = self.get_emdatabase_directory(file)

		key = remove_directories_from_name(file)
		key = strip_file_tag(key)
		
		db_name = "bdb:"+db_directory+"#"+key
		db = db_open_dict(db_name,ro=True)
		
		list_widget.clear()
		#items = DB[key] # NOTE items should be called "db" or something else
		items = db
		keys = items.keys()
		keys.sort() # puts them alphabetical order
		for k in keys:
			if k == '': continue
			_type =db.item_type(k)
			if _type == dict:
				a = QtGui.QListWidgetItem(self.target().database_icon,str(k),list_widget)
				b = EMGenericItem("database_dict",str(k))  # really haven't accommodated for this...
				accrue_public_attributes(a,b)
				#a.type_of_me = "database_dict"
			elif _type == EMData:
				data = db.get_header(k)
				if data["nz"] > 1:
					a = QtGui.QListWidgetItem(self.target().emdata_3d_icon,str(k),list_widget)
					b = EMDBSingleImageItem(self.emdata_3d_entry,db_directory,key,k)
					accrue_public_attributes(a,b)
				else:
					a = QtGui.QListWidgetItem(self.target().emdata_icon,str(k),list_widget)
					b = EMDBSingleImageItem(self.emdata_entry,db_directory,key,k)
					accrue_public_attributes(a,b)
			elif _type == list and len(db[k]) == 2:
				if 1:
					if (isinstance(db[k][0][0],float) or isinstance(db[k][0][0],int)) and (isinstance(db[k][1][0],float) or isinstance(db[k][0][0],int)):
						v = db[k]
				
						a = QtGui.QListWidgetItem(self.target().plot_icon,str(k)+":"+str(v),list_widget)
						b = EMDBPlotItem(self.db_list_plot,db_directory,key,k)
						accrue_public_attributes(a,b)
					else: raise
				else:
					# yes redundant but no time
					v = db[k]
					a = QtGui.QListWidgetItem(self.target().basic_python_icon,str(k)+":"+str(v),list_widget)
					b = EMKeyValueItem("key_value",k,v)
					accrue_public_attributes(a,b)
					
			else:
				#if type(i) in [str,float,int,tuple,list,bool]:
				v = db[k]
				a = QtGui.QListWidgetItem(self.target().basic_python_icon,str(k)+":"+str(v),list_widget)
				b = EMKeyValueItem("key_value",k,v)
				accrue_public_attributes(a,b)
			
		return True
				
	
	def get_last_directory(self,file):
		dtag = get_dtag()
		idx1 = file.rfind(dtag)
		if idx1 > 0:
			ret = file[0:idx1]
		else: return ret
		
		idx2 = ret.rfind(dtag)
		if idx2 > 0:
			ret = ret[idx2+1:]
		
		return ret
		
	def get_emdatabase_directory(self,file):
		'''
		Get the database where EMAN2DB should be opening in order to open the given file
		e.g. if db path is /home/someone/work/EMAN2DB/data.bdb will return /home/someone/work
		'''
		
		idx1 = file.find("EMAN2DB")
		if idx1 > 0:
			return file[0:idx1-1]
		else: return None
		
	
	def load_metadata(self,item,list_widget):
		
		data = item.get_metadata()
		if data != None:
			keys = data.keys()
			keys.sort() # alphabetical order
			for k in keys:
				v = data[k]
				a = QtGui.QListWidgetItem(self.target().basic_python_icon,str(k)+":"+str(v),list_widget)
				b = EMKeyValueItem("key_value",k,v)
				accrue_public_attributes(a,b)
			return True
			
		return False
		
	
	def load_database_interrogation(self,file_name,list_widget):
		dtag = get_dtag()
		split = file_name.split(dtag)
		
		rm = []
		for i,s in enumerate(split):
			if len(s) == 0: rm.append(i)
		
		rm.reverse()
		for k in rm: split.pop(k)
		
		if len(split) > 2 : # must atleast have EMAN2DB/something.bdb/dictionary
			split.reverse() # this probably makes things faster
			for j in range(2,len(split)):
				if split[j] in self.directory_replacements.values():
					break
			else:
				return False
			
			real_name = self.convert_to_absolute_path(file_name)
			db_directory = self.get_emdatabase_directory(real_name)
			
			key = split[j-1]
			item_key = split[j-2]
			
			db_name = "bdb:"+db_directory+"#"+key
			db = db_open_dict(db_name,ro=True)
			try:
				for ii in range(j-2,-1,-1):
					db = db[split[ii]]
			except: 
				print 0	
				return False
			
			if type(db) == dict:
				keys = db.keys()
				keys.sort() # puts them alphabetical order
				for k in keys:
					i = db[k]
					if k == "auto_boxes":
						a = QtGui.QListWidgetItem(self.target().ab_autoboxes_icon,str(k),list_widget)
						b = EMGenericItem("auto_boxes",str(k))
						accrue_public_attributes(a,b)
					elif k == "reference_boxes":
						a = QtGui.QListWidgetItem(self.target().ab_autoboxes_icon,str(k),list_widget)
						b = EMGenericItem("reference_boxes",str(k))
						accrue_public_attributes(a,b)
					elif k == "manual_boxes":
						a = QtGui.QListWidgetItem(self.target().ab_autoboxes_icon,str(k),list_widget)
						b = EMGenericItem("manual_boxes",str(k))
						accrue_public_attributes(a,b)
					elif type(i) in [str,float,int,tuple,list,bool]:
						a = QtGui.QListWidgetItem(self.target().basic_python_icon,str(k),list_widget)
						b = EMGenericItem("value",str(k))
						accrue_public_attributes(a,b)
					elif type(i) == dict:
						a = QtGui.QListWidgetItem(self.target().dict_python_icon,str(k),list_widget)
						b = EMGenericItem("python_dict",str(k))
						accrue_public_attributes(a,b)
					elif type(i) == EMData:
						a = QtGui.QListWidgetItem(self.target().emdata_icon,str(k),list_widget)
						database_dictionary_keys = [split[jj] for jj in range(j-2,-1,-1)]
						database_dictionary_keys.append(k)
						b = EMDBDictSingleImageItem(self.db_dict_emdata_entry,db_directory,key,database_dictionary_keys)
						accrue_public_attributes(a,b)
					elif type(i) == list and len(i) >= 2:
						# warning - I haven't tested this
						# double warning - it is really really untested only leaving for a short while
						a = QtGui.QListWidgetItem(self.target().plot_icon,str(k),list_widget)
						b = EMDBPlotItem(self.db_list_plot,db_directory,key,k)
						accrue_public_attributes(a,b)
					else:
						a = QtGui.QListWidgetItem(self.target().basic_python_icon,str(k),list_widget)
						b = EMGenericItem("value",str(k))
						accrue_public_attributes(a,b)
			elif isinstance(db,EMData):
				print "this shouldn't happen"
				#self.target().preview_data(db)
				return False
			else:
				a = QtGui.QListWidgetItem(self.target().basic_python_icon,str(db),list_widget)
				b = EMGenericItem("value",str(k))
				accrue_public_attributes(a,b)
			return True
				
		else: return False
	
	def is_previewable(self,item):
		return item.type_of_me in self.previewable_types

	def preview_item_wants_to_list(self,item):
		'''
		Sometimes a previeable item will still be able to list information
		'''
		return item.is_previewable()
	
	def do_preview(self,item):
		#if not os.path.isfile(item.full_path): return False
		return item.do_preview(self.target())
	
	def get_emdata_header(self,item):
		return item.get_metadata()

	def load_database_variables(self,directory,list_widget):
		pass

def accrue_public_attributes(to_this,from_this):
	'''
	Accrues the public attributes of the from_this object to to_this.
	By public it is meant that all attributes beginning with "__" are ignored
	'''
	for attr in dir(from_this):
		if len(attr) > 2 and attr[:2] == "__": continue
		setattr(to_this,attr,getattr(from_this,attr))
	
#	
#class MultiDelegate(object):
#	'''
#	This class was the basis of the idea behind the accrue_public_attributes function
#	'''
#	def __init__(self, *args):
#	    object.__setattr__(self,"_objects",args)
#	def __getattr__(self, attr):
#	    for obj in self._objects:
#	         if attr in dir(obj):
#	             return getattr(obj, attr)
#	     types = ",".join(obj.__class__.__name__ for obj in self._objects)
#	 raise AttributeError, "%s object has no attribute '%s'" % (types, attr)
#	 def __setattr__(self, attr, value):
#	     # but you could do something more useful here too
#	 raise TypeError, "Can't set attributes of MultiDelegate"


class EMDB2DImageStackItem(EMListingItem):
	def __init__(self,type,db_directory,db):
		EMListingItem.__init__(self)
		self.type_of_me = type
		self.database_directory = db_directory
		self.database = db
		
		self.context_menu_options["Delete"] = self.delete
		self.context_menu_options["Save As"] = self.save_as
		
	def do_preview(self,target):
		db_name = self.get_path()
		target.preview_data(None,db_name)	
		return True
	
	def is_previewable(self): return True

	def delete(self,target):
		db_name = self.get_path()
		try:
			remove_file(db_name)
			return True
		except: return False
		
	def save_as(self,target):
		save_data(self.get_emdata())
#		db_name = self.get_path()
#		return save_as(db_name,target.application())
#	
	def get_emdata(self):
		from emimagemx import EMDataListCache
		return EMDataListCache(self.get_path())
	
	def get_path(self):
		return "bdb:"+self.database_directory+"#"+self.database

class EMDBDictSingleImageItem(EMListingItem):
	def __init__(self,type,db_directory,key,dictionary_keys):
		EMListingItem.__init__(self)
		self.type_of_me = type
		self.database_directory = db_directory
		self.database = key
		self.database_dictionary_keys = database_dictionary_keys
	
	def do_preview(self,target):
		db_name = self.get_path()
		db = db_open_dict(db_name)

		for key in self.database_dictionary_keys:
			db = db[key]
		# db should now be an EMData
		target.preview_data(db)
		return True
	
	def is_previewable(self): return True
	
	def get_metadata(self):
		db_name = self.get_path()
		db = db_open_dict(db_name,ro=True)
		for key in self.database_dictionary_keys:
			db = db[key]
		
		return db.get_attr_dict()
	
	def get_emdata(self):
		e = EMData()
		e.read_image(self.get_path(),0)
		return e

	def get_path(self):
		return "bdb:"+self.database_directory+"#"+self.database

class EMDBSingleImageItem(EMListingItem):
	def __init__(self,type,db_directory,key,k):
		EMListingItem.__init__(self)
		self.type_of_me = type
		self.database_directory = db_directory
		self.database = key
		self.database_key = k
		self.context_menu_options["Save As"] = self.save_as

	def do_preview(self,target):
		db_name = self.get_path()
		db = db_open_dict(db_name,ro=True)
		
		data = db[self.database_key]
		target.preview_data(data)
		return True
	
	def is_previewable(self): return True

	def get_metadata(self):
		db_name = self.get_path()
		db = db_open_dict(db_name,ro=True)
		data = db.get_header(self.database_key)
		return data
	
	def get_emdata(self):
		db_name = self.get_path()
		return EMData(db_name,self.database_key)
		
	def save_as(self,target):
		db_name = self.get_path()
		return save_data(self.get_emdata())
	
	def get_path(self):
		return "bdb:"+self.database_directory+"#"+self.database

class EMDBPlotItem(EMListingItem):
	def __init__(self,type,db_directory,key,k):
		EMListingItem.__init__(self)
		self.type_of_me = type
		self.db_dir = db_directory
		self.db = key
		self.db_key = k
		
	def do_preview(self,target):
		db_name = "bdb:"+self.db_dir+"#"+self.db
		db = db_open_dict(db_name)
		plot = db[self.db_key]
		target.preview_plot_list(self.db_key,plot)
		return True
	
	def is_previewable(self): return True
		
class EMKeyValueItem(EMListingItem):
	def __init__(self,type,k,v):
		EMListingItem.__init__(self)
		self.type_of_me = type
		self.key = k
		self.value = v

class EMDBDirectoryItem(EMListingItem):
	def __init__(self,type,real_directory):
		EMListingItem.__init__(self)
		self.type_of_me = type
		self.full_path = real_directory

class EMGenericItem(EMListingItem):
	'''
	This one is just where things end up if they have no special treatment, and/or only the type is relevant
	'''
	def __init__(self,type,key):
		EMListingItem.__init__(self)
		self.type_of_me = type
		self.key = key
		
	def do_preview(self,target): pass
	

		
class EMBrowserDialog(EMSelectorDialog):
	def __init__(self,target):
		EMSelectorDialog.__init__(self,target,False,False)
		self.preview_options.setCurrentIndex(1)
		self.preview_options_changed(self.preview_options.currentText())
		self.ok_button.setEnabled(False)
		self.cancel_button.setEnabled(False)

class EMBrowserModule(EMQtWidgetModule):
	def __init__(self):
		self.widget = EMBrowserDialog(self)
		EMQtWidgetModule.__init__(self,self.widget)

	def __del__(self):
#		import sys
#		print "browser module death", sys.getrefcount(self.widget)
		self.widget.destroy()

app = None
def on_done(string_list):
	if len(string_list) != 0:
		for s in string_list:
			print s,
		print
	app.quit()

def on_cancel(string_list):
	app.quit()


if __name__ == '__main__':
	em_app = EMStandAloneApplication()
	#dialog = EMSelectorDialog(None,em_app)
	em_qt_widget = EMSelectorModule()
	QtCore.QObject.connect(em_qt_widget,QtCore.SIGNAL("ok"),on_done)
	QtCore.QObject.connect(em_qt_widget,QtCore.SIGNAL("cancel"),on_cancel)
	em_app.show()
	em_app.execute()


