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
EMAN2DB = "EMAN2DB"

class EMSelectorModule(EMQtWidgetModule):
	def __init__(self,single_selection=False,save_as_mode=True):
		self.widget = EMSelectorTemplate(QtGui.QDialog)(self,single_selection)
		EMQtWidgetModule.__init__(self,self.widget)
		self.widget.setWindowTitle("EMAN2 Selector")
		
#	def __del__(self):
#		import sys
#		print "selector module death", sys.getrefcount(self.widget)
		
	def exec_(self):
		'''
		Wraps self.widget.exec_
		@return a list of selected filenames
		'''
		return self.widget.exec_()
		

def EMSelectorTemplate(Type):
	''' 
	This is just an example of using nested scope inheritance to achieve templated inheritance, ala what you'd do in C++
	It solved a problem for me, but the best solution may be to break the functionality into two (written, not templated) classes
	Valid Types are QtGui.Dialog and QtGui.QWidget
	'''
	class EMSelector(Type):
		'''
		This class is hybrid. It can be used both as a file system browser
		and as a save-file dialog. To see an example of using it as a system browser, see the
		EMBrowserModule. To see an example of using it as a save-file dialog, see emsave.py.
		Also, emsprworkflow uses this class as a browser for selecting files.  
		I wish I had time to correct the design... but that said, it's not too far from where it should
		be - the code is relatively clean, and correcting the design involves breaking some of the member
		functions into separate classes.
		'''
		def __init__(self,module,single_selection=False):
			'''
			@param module - should be an EMSelectorModule
			@param application - should be a type of application as supplied in emapplication
			@param single_selection - should selections be limited to singles?
			@param save_as_mode - should be True if using this object as a pure dialog (using exec_).
			'''
			Type.__init__(self,None)
			self.setFocusPolicy(Qt.StrongFocus)
			self.module=weakref.ref(module) # Avoid strong cycle
			self.desktop_hint = "dialog" # So the desktop knows how to display this
			self.single_selection = single_selection # Flag indicating single selection in interface
			self.browse_delegates = [EMBDBReader(self), EMFileSystemReader(self)] # Object capable of returning listed items based on url
			
			self.hbl = QtGui.QVBoxLayout(self)
			self.hbl.setMargin(0)
			self.hbl.setSpacing(6)
			self.hbl.setObjectName("hbl")
			
			self.__init_icons()
			
			self.__init_filter_combo()
			
			self.current_force = None # use to keep track of what is currently being forced, either 2D plotting, 2D image showing, or neither
			self.selections = []
			self.current_list_widget = None
			self.lock = True
			self.list_widgets = []
			self.previews = [] # keeps track of all of the preview windows
			self.module_events = [] # used to coordinate signals from the modules, especially close events, to free memory
			self.list_widget_data= [] # entries should be tuples containing (current folder item)
			self.splitter = QtGui.QSplitter(self)
			self.splitter.setChildrenCollapsible(False)
			self.starting_directory = e2getcwd()
			self.historical_starting_directory = e2getcwd() # just incase anyone ever needs it (True). This should never be changed

			self.__add_list_widget()
			self.__add_list_widget()
			if Type != QtGui.QDialog: self.__add_list_widget() # this is just to distinguish the saving mode from the browsing mode
			
			self.hbl.addWidget(self.splitter,1)
			
			self.__load_url(self.starting_directory,self.list_widgets[0])
	
			self.dialog_mode = False
			if Type == QtGui.QDialog:
				hbl2=QtGui.QHBoxLayout()
				hbl2.setMargin(0)
				hbl2.setSpacing(2)
				self.selection_label = QtGui.QLabel("Save As",self)
				hbl2.addWidget(self.selection_label)
				self.save_as_line_edit = QtGui.QLineEdit("",self)
				hbl2.addWidget(self.save_as_line_edit,0)
				self.hbl.addLayout(hbl2)
				self.dialog_mode = True
				self.validator = None # an optional feature which makes the behavior of the dialog more sophisticated - see emsave.py
				self.dialog_result = ""
	
			self.bottom_hbl = QtGui.QHBoxLayout()
			self.bottom_hbl.addWidget(self.filter_text,0)
			self.bottom_hbl.addWidget(self.filter_combo,1)
			self.__init_buttons()
			self.bottom_hbl.addWidget(self.cancel_button,0)
			self.bottom_hbl.addWidget(self.ok_button,0)
			self.hbl.addLayout(self.bottom_hbl)
			
			if Type != QtGui.QDialog: 
				self.gl_image_preview = None
				
				self.bottom_hbl2 = QtGui.QHBoxLayout()
				self.__init_preview_options()
				self.bottom_hbl2.addWidget(self.preview_options,0)
				self.hbl.addLayout(self.bottom_hbl2)
				
				self.bottom_hbl3 = QtGui.QHBoxLayout()
				self.__init_plot_options()
				self.bottom_hbl3.addWidget(self.replace,0)
				self.bottom_hbl3.addWidget(self.include,0)
				
				self.groupbox = QtGui.QGroupBox("Plot/3D options")
				self.groupbox.setLayout(self.bottom_hbl3)
				self.groupbox.setEnabled(False)
				
				self.bottom_hbl2.addWidget(self.groupbox)
			else:
				self.ok_button.setDefault(True)
			
			if Type != QtGui.QDialog:
				self.resize(480,480)
			else:
				self.resize(480,400)
			
			self.lock = False
			
			self.paint_events = 0
			
			self.timer_interval = 500 # half a second
			self.timer = QtCore.QTimer()
			QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.time_out) # for auto refresh
			
			self.timer.start(self.timer_interval)
			
			self.selected_files = []
			
		def __del__(self):
			pass
		
		def set_selection_text(self,text):
			'''
			Selection label is a QLabel, by default its text is "Save As", but you can change it
			'''
			self.selection_label.setText(text)
		
		def set_validator(self,validator):
			'''
			Sets the validator which is used only in save_as mode (when this is being used as a file selecting dialog)
			An example validator is an emsave.EMSaveImageValidator
			'''
			self.validator = validator
		
		if Type == QtGui.QDialog:
			def exec_(self):
				'''
				Wraps QtGui.QDialog.exec_
				@return a list of selected filenames
				'''
				QtGui.QDialog.exec_(self)
				return self.dialog_result
		else:
			def exec_(self):
				print "exec performs no function - inherit from a QtGui.QDialog instead"
				return ""
		
		def time_out(self):
			'''
			This function takes care of automatic updates - if the file system changes then so does
			the information we display
			'''
			if self.lock: return
			self.lock = True
			for idx,widget in enumerate(self.list_widgets):
				if widget.get_mod_time() != None and widget.get_delegate() != None:
					mod_time = widget.get_delegate().url_mod_time(widget.get_url())
					if mod_time != widget.get_mod_time():
						old_items = [widget.item(i) for i in range(widget.count())]
						new_items = widget.get_delegate().get_items(widget.get_url())
						old_set = set([str(item.text()) for item in old_items])
						new_set = set([str(item.text()) for item in new_items])
						
						added = new_set - old_set
						removed = old_set - new_set
						
						if len(added) > 0:
							for k in added:
								new_item = (item for item in new_items if item.text() == k).next()
								widget.addItem(new_item)

						if len(removed) > 0:
							for k in removed:
								old_item = (item for item in old_items if item.text() == k).next()
								if old_item.isSelected():
									for j in range(idx+1,len(self.list_widgets)):
										self.list_widgets[j].clear()
								widget.takeItem(widget.row(old_item))

						widget.set_mod_time(mod_time)
						self.lock = False
						return
							
			self.lock = False
					
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
			self.euler_icon = QtGui.QIcon(directory + "eulerxplor.png")
			
		def __init_plot_options(self):
			self.replace = QtGui.QRadioButton("Replace")
			self.include = QtGui.QRadioButton("Include")
			self.include.setChecked(True)
	
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
		
		def previews_allowed(self):
			if self.dialog_mode: return False
			
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
			if self.dialog_mode:
				self.accept()
			else:
				self.module().emit(QtCore.SIGNAL("cancel"),self.selections)
		
		def ok_button_clicked(self,bool):
			if self.dialog_mode:
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
					v = directory
					# strip out the EMAN2DB
					# assumption that only one "EMAN2DB" exists in the string
					for cand in ["EMAN2DB/","EMAN2DB"]:
						l = len(cand)
						n = v.find(cand)
						if  n != -1:
							v = v[:n]+v[n+l:]
							break # if we find one the first then there is no need t
					
					if len(v) > 0 and v[-1] == "/":
						v = v[:-1]
					
					ret = "bdb:"+v+"#"+last_bit
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
			self.filter_combo.addItem("*.hdf,*.img,*.hed,*.spi,bdb:,*.tif,*.mrc,*.dm3, *.pif, *.rec")
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
				self.__load_url(directory,self.list_widgets[i])
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
				if nz == 1 and not self.dialog_mode: 
					menu.addAction("Preview Subset")
					if len(selected_items) < 256:
						menu.addAction("View Subset In 3D")
			elif not self.dialog_mode: # dialog mode stops all previews
				if md != None and len(selected_items) == 1:
					if nz > 1:
						menu.addSeparator()
						menu2 = QtGui.QMenu("Open With")
						menu2.addAction(self.emdata_3d_icon,"Volume Viewer")
						menu2.addAction(self.emdata_3d_icon,"Slice Viewer")
						menu2.addAction(self.emdata_icon,"Single 2D")
						menu2.addAction(self.emdata_matrix_icon,"Multi 2D")
						menu2.addAction(self.plot_icon,"Plot 2D")
						menu.addMenu(menu2)
					elif nz == 1:
						menu.addSeparator()
						menu2 = QtGui.QMenu("Open With")
						menu2.addAction(self.plot_icon,"2D Plot")
						if ny < 1025 and nx < 1025:
							menu2.addAction(self.emdata_3d_icon,"3D Viewer")
							menu2.addAction(self.plot_icon,"3D Plot")
							menu.addMenu(menu2)
				elif len(selected_items) == 1 and selected_items[0].is_2d_stack():
					menu.addSeparator()
					menu2 = QtGui.QMenu("Open With")
					menu2.addAction(self.emdata_icon,"Single 2D")
					if EMUtil.get_image_count(selected_items[0].get_url()) < 1001:
						menu2.addAction(self.plot_icon,"2D Plot")
						menu2.addAction(self.emdata_3d_icon,"3D Viewer")
					md = EMData()
					md.read_image(selected_items[0].get_url(),0,True)
					md = md.get_attr_dict()
					if md.has_key("xform.projection"):
						# the assumption is that they all have the xform.projection header attribute, which could fail
						menu2.addAction(self.euler_icon,"Euler Viewer")
					menu.addMenu(menu2)
				
				elif selected_items[0].is_2d_plot():
					menu.addSeparator()
					menu2 = QtGui.QMenu("Open With")
					menu2.addAction(self.plot_icon,"3D Plot")
					menu.addMenu(menu2)
					
				
	#		if nz == 1 and ny > 1:	menu.addAction("Open in e2boxer")
				
			QtCore.QObject.connect(menu,QtCore.SIGNAL("triggered(QAction*)"),self.menu_action_triggered)
			self.action_list_widget = l # only set if the menu acutally triggers
			menu.exec_(event.globalPos())
		
		def menu_action_triggered(self,action):
			items = self.menu_selected_items
			
			cont = self.__check_action(action.text(),items) # this kind of breaks the OO design, but it's necessary in the current framework
			if not cont: return
			total = len(items)
			
			if action.text() == "Save As Stack":
				save_data(self.menu_selected_items)
			elif action.text() == "Preview Subset":
				data = []
				for item in self.menu_selected_items:
					data.append(item.get_emdata())
				self.preview_data(data,"")
			elif action.text() == "Euler Viewer":
				get_application().setOverrideCursor(Qt.BusyCursor)
				self.preview_euler_view( self.menu_selected_items[0].get_url())
				get_application().setOverrideCursor(Qt.ArrowCursor)
			elif action.text() == "2D Plot" or action.text() == "Plot 2D":
				get_application().setOverrideCursor(Qt.BusyCursor)
				self.preview_plot( self.menu_selected_items[0].get_url())
				get_application().setOverrideCursor(Qt.ArrowCursor)
			elif action.text() == "3D Plot":
				from emplot3d import EMPlot3DModule
				get_application().setOverrideCursor(Qt.BusyCursor)
				if self.menu_selected_items[0].is_2d_plot():
					data =	EMPlot3DModule.parse_txt_file(self.menu_selected_items[0].get_url())
					d = [ i*.01 for i in xrange(0,len(data[0])) ]
					self.preview_item_with_module([d,data[1],data[1]],EMPlot3DModule)
				else:
					data = self.menu_selected_items[0].get_emdata()
					v = data.get_data_as_vector()
					x = [i for i in range(data.get_xsize()) for j in range(data.get_ysize())]
					y = [j for i in range(data.get_xsize()) for j in range(data.get_ysize())]
					self.preview_item_with_module([x,y,v],EMPlot3DModule)
				get_application().setOverrideCursor(Qt.ArrowCursor)
			elif action.text() == "View Subset In 3D":
				data1 = self.menu_selected_items[0].get_emdata()
				new_data = EMData(data1.get_xsize(),data1.get_ysize(),len(self.menu_selected_items))
				new_data.insert_clip(data1,[0,0,0])
				for i in xrange(1,len(self.menu_selected_items)):
					new_data.insert_clip(self.menu_selected_items[i].get_emdata(),[0,0,i])
					
				from emimage3d import EMImage3DModule
				self.preview_item_with_module(new_data,EMImage3DModule)
				get_application().setOverrideCursor(Qt.ArrowCursor)
			elif action.text() == "Multi 2D":
				from emimagemx import EMImageMXModule
				item = self.menu_selected_items[0]
				get_application().setOverrideCursor(Qt.BusyCursor)
				data = item.get_url()
				self.preview_item_with_module(data,EMImageMXModule)
				get_application().setOverrideCursor(Qt.ArrowCursor)
			elif action.text() == "Single 2D":
				from emimage2d import EMImage2DModule
				item = self.menu_selected_items[0]
				get_application().setOverrideCursor(Qt.BusyCursor)
				data = item.get_emdata()
				self.preview_item_with_module(data,EMImage2DModule)
				get_application().setOverrideCursor(Qt.ArrowCursor)
			elif action.text() == "Volume Viewer":
				from emimage3dvol import EMVolumeModule
				item = self.menu_selected_items[0]
				get_application().setOverrideCursor(Qt.BusyCursor)
				data = item.get_emdata()
				data.process_inplace("normalize")
				self.preview_item_with_module(data,EMVolumeModule)
				get_application().setOverrideCursor(Qt.ArrowCursor)
			elif action.text() == "Slice Viewer":
				from emimage3dslice import EM3DSliceViewerModule
				item = self.menu_selected_items[0]
				get_application().setOverrideCursor(Qt.BusyCursor)
				data = item.get_emdata()
				data.process_inplace("normalize")
				self.preview_item_with_module(data,EM3DSliceViewerModule)
				get_application().setOverrideCursor(Qt.ArrowCursor)
			elif action.text() == "3D Viewer":
				from emimage3d import EMImage3DModule
				from emimagemx import EMDataListCache
				item = self.menu_selected_items[0]
				get_application().setOverrideCursor(Qt.BusyCursor)
				data = item.get_emdata()
				if isinstance(data,EMDataListCache):
					data0 = data[0]
					new_data = EMData(data0.get_xsize(),data0.get_ysize(), len(data))
					new_data.insert_clip(data0,[0,0,0])
					for i in xrange(1,len(data)):
						new_data.insert_clip(data[i],[0,0,i])
					data = new_data
				elif data.get_zsize() == 1:
					new_data = EMData(data.get_xsize(),data.get_ysize(),2)
					new_data.insert_clip(data,[0,0,0])
					new_data.insert_clip(data,[0,0,1])
					data = new_data
				 
				self.preview_item_with_module(data,EMImage3DModule)
				get_application().setOverrideCursor(Qt.ArrowCursor)
					
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
#			
#		def __post_action(self,action_str,items_acted_on):
#			if action_str == "Delete":
#				c_list_widget = self.action_list_widget
#				items = [c_list_widget.item(i) for i in range(c_list_widget.count())]
#				rm_indices = []
#				for j,i in enumerate(items):
#					if i in items_acted_on:
#						rm_indices.append(j)
#				
#				rm_indices.reverse()
#				for idx in rm_indices:
#					c_list_widget.takeItem(idx)
#				
#				for ii,l in enumerate(self.list_widgets):
#					if l == c_list_widget:
#						break
#				
#				ii += 1
#				if ii <= (len(self.list_widgets) -1):
#					for i in range(ii,len(self.list_widgets)):
#						self.list_widgets[i].clear()
				
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
			if list_widget == None:	list_widget = EMListWidget()
			
			list_widget.contextMenuEvent = self.list_widget_context_menu_event
			
			if self.single_selection:list_widget.setSelectionMode(QtGui.QAbstractItemView.SingleSelection)
			else: list_widget.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
			list_widget.setMouseTracking(True)	
			self.list_widgets.append(list_widget)
			self.splitter.addWidget(list_widget)
			
			self.list_widget_data.append(None)
			
			QtCore.QObject.connect(list_widget, QtCore.SIGNAL("itemDoubleClicked(QListWidgetItem*)"),self.list_widget_dclicked)
			#QtCore.QObject.connect(list_widget, QtCore.SIGNAL("itemPressed(QListWidgetItem*)"),self.list_widget_clicked)
			#QtCore.QObject.connect(list_widget, QtCore.SIGNAL("currentRowChanged (int)"),self.list_widget_row_changed)
			#QtCore.QObject.connect(list_widget, QtCore.SIGNAL("paintEvent (int)"),self.list_widget_row_changed)
			QtCore.QObject.connect(list_widget, QtCore.SIGNAL("itemEntered(QListWidgetItem*)"),self.list_widget_item_entered)
			#QtCore.QObject.connect(list_widget, QtCore.SIGNAL("currentItemChanged(QListWidgetItem*,QListWidgetItem*)"),self.list_widget_current_changed)
			#QtCore.QObject.connect(list_widget, QtCore.SIGNAL("itemChanged(QListWidgetItem*)"),self.list_widget_item_changed)
			#\QtCore.QObject.connect(list_widget, QtCore.SIGNAL("itemActivated(QListWidgetItem*)"),self.list_widget_item_activated)
			#QtCore.QObject.connect(list_widget, QtCore.SIGNAL("activated(QModelIndex)"),self.activated)
			#QtCore.QObject.connect(list_widget, QtCore.SIGNAL("itemSelectionChanged()"),self.selection_changed)
			
#		def list_widget_current_changed(self,new,prev):
#			if new != None:
#				if self.current_list_widget == new.listWidget():
#					self.list_item_selected(new,False)
		
		def __go_back_a_directory(self):
			self.lock = True
			dtag = get_dtag()
			
			new_dir = self.starting_directory[0:self.starting_directory.rfind(dtag)]

			if len(new_dir) == 0: new_dir = get_dtag()
			elif  new_dir[-1] == ":": new_dir += get_dtag() # C: becomes C:/
			elif len(new_dir) == 1 and new_dir != get_dtag(): new_dir += ":/"# C becomes C:/

				
			self.starting_directory = new_dir
			for j in range(0,len(self.list_widgets)):
				widget = self.list_widgets[j]
				widget.clear()
				self.list_widget_data[j] = None
					
			self.__load_url(self.starting_directory,self.list_widgets[0])
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
				li.set_mod_time(lip.get_mod_time())
				li.set_url(lip.get_url())
				li.set_delegate(lip.get_delegate())

				self.list_widget_data[i] = li.item(old_row)
				directory += dtag + str(self.list_widget_data[i].text())
			
			# The last widget must have these attributes removed if they existed,
			# or else the time_out function will produce incorrect results
			last_widget = self.list_widgets[-1]

			a = EMUpArrowItem(None,"../")
			self.list_widgets[0].insertItem(0,a)
			
			self.lock = False
			self.hide_preview()
			
		def __update_selections(self):
			'''
			Makes the list of currently selected files accurate and up to date. Called when
			something has been clicked in a a list widget
			'''
			
			#print self.current_list_widget
			# get the directory 
#			dtag = get_dtag()
#			directory = self.starting_directory+dtag
#			idx = 0
#			for i,list_widget in enumerate(self.list_widgets):
#				if list_widget == self.current_list_widget: 
#					break
#				if (self.list_widget_data[i] == None):
#					# this means it could be a header attribute, for example
#					return
#				directory += str(self.list_widget_data[i].text()) + dtag
#			else:
#				print "no list widget has focus?"
#				return
			
			# now make the list of selections reflect what is need to display them using
			
			self.selections = []
			items = self.current_list_widget.selectedItems()
			if len(items) > 0:
				a = items[0]
#				previewable,previewer= a.is_previewable(),a.get_delegate()
#				if previewable:
				self.selections = self.previewable_full_paths(items)
			
			# Make sure the save as line is updated if we're in that mode
			if self.dialog_mode:
				text = ""
				for i,name in enumerate(self.selections):
					if i > 0:
						text+= "; "
					text += base_name(name)
					
				self.save_as_line_edit.setText(text)
		
		def previewable_full_paths(self,items):
			return [item.get_url() for item in items if item.is_previewable()]
		
		def hide_preview(self):
			if self.dialog_mode: return
			
			if self.gl_image_preview  != None:
				get_application().close_specific(self.gl_image_preview)
		
		def list_widget_item_entered(self,item):
			list_widget = item.listWidget()
			if list_widget != self.current_list_widget:
				if self.current_list_widget != None:
					QtCore.QObject.disconnect(self.current_list_widget,QtCore.SIGNAL("itemSelectionChanged()"), self.current_item_changed)
				self.current_list_widget = item.listWidget()
				QtCore.QObject.connect(self.current_list_widget,QtCore.SIGNAL("itemSelectionChanged()"), self.current_item_changed)
#				
		def current_item_changed(self):
			'''
			This function handles any change in current item
			'''
			if self.lock: return
			if self.current_list_widget == None: return 
			item = self.current_list_widget.currentItem()
			if item != None:
				self.list_item_selected(item)
				
		def keyPressEvent(self,event):
			if event.key() == Qt.Key_F1:
				try:
					import webbrowser
					webbrowser.open("http://blake.bcm.edu/emanwiki/EMAN2/Programs/emselector")
				except:
					pass
		
#		def list_widget_clicked(self,item):
#			self.list_item_selected(item)
#		
		def list_widget_dclicked(self,item):
			self.preview_item(item)

		def list_item_selected(self,item):
			if self.lock : return
			#if self.current_list_widget == None: return
			if item.get_name() == EMGenericItem.NAME: return #it's just a value in the db
			
			self.__update_selections()
			if item == None: return		
			
			if item.get_name() == EMUpArrowItem.NAME: 
				self.__go_back_a_directory()
				return
			
			self.current_list_widget = item.listWidget()
			if self.current_list_widget  == self.list_widgets[-1]:
					self.lock = True
					self.list_widget_data[-1] = item
					self.__go_forward_a_directory()
					self.__load_url_from_item(self.list_widgets[-1],item)
					self.lock = False
					return
			
			idx = (i for i in range(len(self.list_widgets)) if self.list_widgets[i] == self.current_list_widget).next()
			if self.__load_url_from_item(self.list_widgets[idx+1],item):
				self.list_widget_data[idx] = item
				self.list_widget_data[idx+1] = None
				
			for i in range(idx+2,len(self.list_widgets)):
				self.list_widgets[i].clear()
				self.list_widget_data[i] = None
		
		def preview_item(self,item):
			'''
			previews the item (loads an appropriate display module) if possible
			Returns True if the item was loaded into a display module, otherwise returns False
			'''
			if item.is_previewable():
				get_application().setOverrideCursor(Qt.BusyCursor)
				preview_occured = item.do_preview(self)
				get_application().setOverrideCursor(Qt.ArrowCursor)
				return preview_occured
			return False
				
		def preview_euler_view(self, full_path):
			if self.dialog_mode: return
			from emimagemx import EMLightWeightParticleCache
			from emimage3dsym import EM3DSymViewerModule
			if self.single_preview_only():
				if not isinstance(self.gl_image_preview,EM3DSymViewerModule):
					if self.gl_image_preview != None: get_application().close_specific(self.gl_image_preview)
					self.gl_image_preview = EM3DSymViewerModule()
					QtCore.QObject.connect(self.gl_image_preview.emitter(), QtCore.SIGNAL("module_closed"),self.preview_module_closed)
				
				self.gl_image_preview.set_emdata_list_as_data(EMLightWeightParticleCache.from_file(full_path))
				get_application().show_specific(self.gl_image_preview)
				self.gl_image_preview.updateGL()
				
			else:
				preview = EM3DSymViewerModule(get_application())
				preview.set_emdata_list_as_data(EMLightWeightParticleCache.from_file(full_path))
				get_application().show_specific(preview)
				
		def preview_plot(self,filename):
			if self.dialog_mode: return
			
			if self.single_preview_only():
				if not isinstance(self.gl_image_preview,EMPlot2DModule):
					if self.gl_image_preview != None: get_application().close_specific(self.gl_image_preview)
					self.gl_image_preview = EMPlot2DModule(get_application())
					QtCore.QObject.connect(self.gl_image_preview.emitter(), QtCore.SIGNAL("module_closed"),self.preview_module_closed)
		
				self.gl_image_preview.set_data_from_file(filename,self.replace.isChecked())
				get_application().show_specific(self.gl_image_preview)
				self.gl_image_preview.updateGL()
				
			else:
				preview = EMPlot2DModule(get_application())
				preview.set_data_from_file(filename,self.replace.isChecked())
				get_application().show_specific(preview)
				
		def preview_plot_list(self,title,list_data):
			if self.dialog_mode: return
			
			if self.single_preview_only():
				if not isinstance(self.gl_image_preview,EMPlot2DModule):
					if self.gl_image_preview != None: get_application().close_specific(self.gl_image_preview)
					self.gl_image_preview = EMPlot2DModule(get_application())
					QtCore.QObject.connect(self.gl_image_preview.emitter(), QtCore.SIGNAL("module_closed"),self.preview_module_closed)
		
				self.gl_image_preview.set_data(title,list_data,self.replace.isChecked())
				get_application().show_specific(self.gl_image_preview)
				self.gl_image_preview.updateGL()
				
			else:
				preview =EMPlot2DModule(get_application())
				preview.set_data_from_file(filename,self.replace.isChecked())
				get_application().show_specific(preview)
				
		def preview_data(self,a,filename=""):
			if self.dialog_mode: return
			
			from emimage import EMImageModule, EMModuleFromFile
			
			using_file_names_only = False
			if a == None: using_file_names_only = True # For the image matrix, you can load large image stacks if you specify only the file name
			
			if self.single_preview_only() and len(self.previews) != 0:
				old_preview = self.previews[-1] # this means we always choose the last preview if the user suddenly goes from multipreview to single preview
				if not using_file_names_only:
					preview = EMImageModule(data=a,app=get_application(),old=old_preview,filename=filename,replace=self.replace.isChecked())
				else:
					preview = EMModuleFromFile(filename,application=get_application(),old=old_preview)
				if preview != old_preview:
					self.module_closed(old_preview)
					old_preview.closeEvent(None)
					self.previews.append(preview)
					self.module_events.append(ModuleEventsManager(self,preview))
					try: preview.optimally_resize()
					except: pass
			else:
				if not using_file_names_only:
					preview = EMImageModule(data=a,app=get_application(),filename=filename)
				else:
					preview = EMModuleFromFile(filename,application=get_application())
				self.previews.append(preview)
				self.module_events.append(ModuleEventsManager(self,preview))
				try: preview.optimally_resize()
				except: pass
						
			get_application().show_specific(preview)
					
			preview.updateGL()
				
			get_application().setOverrideCursor(Qt.ArrowCursor)
			
		def preview_item_with_module(self,data,module_type):
			get_application().setOverrideCursor(Qt.BusyCursor)
			if self.single_preview_only() and len(self.previews) != 0:
				old_preview = self.previews[-1]
				if isinstance(old_preview,module_type):
					old_preview.set_data(data)
					old_preview.updateGL()
					return
				else:
					self.module_closed(old_preview)
					old_preview.closeEvent(None)
		
			preview = module_type()
			preview.set_data(data)
			preview.updateGL()
			self.previews.append(preview)
			self.module_events.append(ModuleEventsManager(self,preview))
			get_application().show_specific(preview)
			get_application().setOverrideCursor(Qt.ArrowCursor)
		
		def preview_module_closed(self):
			'''
			Slot for signal that is emitted when a module is closed
			'''
			self.gl_image_preview = None
		
		def module_closed(self,module):
			import sys
			for i,mod in enumerate(self.previews):
				if mod == module:
					p = self.previews.pop(i)
					mod = self.module_events.pop(i)
					mod.disconnect_object()
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
						
#		def __is_non_empty_directory(self,s):
#			'''
#			Returns true if s is a non empty directory
#			'''
#	
#			dirs, files  = get_files_and_directories(s)
#			files = self.filter_strings(files)
#			file_length = len(files)
#			if file_length == 0: file_length = len(dirs)
#		
#			if file_length != 0: return True
#			
#			return False
		
#		def __load_database_directory(self,database_name,list_widget):
#			
#			print "loading database directory",database_name
#		
		
#		def __refresh_directory_data(self):
#			#self.__load_url(self.ldd_directory,self.ldd_list_widget,self.ldd_item)
#			self.__load_url_from_item(self.ldd_list_widget,self.ldd_item)
#		
		def __load_url(self,url,list_widget):
			'''
			A way to load items in the list_widget using a url
			especially at start up
			'''
			get_application().setOverrideCursor(Qt.BusyCursor)
			
			list_widget.clear()
			if (list_widget == self.list_widgets[0]):
				self.lock = True
				a = EMUpArrowItem(None,"../")
				list_widget.addItem(a)
				self.lock = False

			for delegate in self.browse_delegates:
				if delegate.handles_url(url): 
					items = delegate.get_items(url)
					for item in items: list_widget.addItem(item)
					ret = True
					list_widget.set_url(url)
					list_widget.set_mod_time(delegate.url_mod_time(url))
					list_widget.set_delegate(delegate)
					break
			else:
				print "unknown url",url
				ret = False
				
			get_application().setOverrideCursor(Qt.ArrowCursor)
			return ret
		
		def __load_url_from_item(self,list_widget,item):
			'''
			A way to load items in the list_widget using data from the item
			i.e. user clicked on the item
			'''
			
			url = item.get_url()
			list_widget.clear()
			
			ret = False
			get_application().setOverrideCursor(Qt.BusyCursor)
			if item.loads_metadata():
				items = item.get_metadata_items()
				for item in items: list_widget.addItem(item)
				
				ret = True
			else:
				for delegate in self.browse_delegates:
					if delegate.handles_url(url): 
						items = delegate.get_items(url)
						for item in items: list_widget.addItem(item)
						ret = True
						list_widget.set_url(url)
						list_widget.set_mod_time(delegate.url_mod_time(url))
						list_widget.set_delegate(delegate)
						break
				else:
					print "unknown url",url
					ret = False
			
			get_application().setOverrideCursor(Qt.ArrowCursor)
			return ret	
					
	
	return EMSelector

class EMListWidget(QtGui.QListWidget):
	'''
	Customized ListWidget as displayed in the browser
	
	'''
	def __init__(self,*args):
		QtGui.QListWidget.__init__(self,*args)
		self.reset_vars()
	
	def clear(self):
		self.reset_vars()
		QtGui.QListWidget.clear(self)
		
	def reset_vars(self):
		'''
		Can't use reset as a name because the QListWidget has it already, and it's vital
		'''
		self.mod_time = None # keeps track of mod time, used to provoke auto list widget repopulation
		self.url = None # keep track of the url that was used to populate the list widget
		self.delegate = None # the delegate the populated the list widget
		
	def set_mod_time(self,mod_time): self.mod_time = mod_time
	def get_mod_time(self): return self.mod_time
	
	def set_url(self,url): self.url = url
	def get_url(self): return self.url
	
	def get_delegate(self): return self.delegate
	def set_delegate(self,delegate): self.delegate = delegate

class EMBrowseDelegate:
	'''
	Base class for objects that can read urls and return lists of ListWidgetItems
	to the EMSelector
	'''
	def __init__(self): pass
	
	def handles_url(self,url):
		'''
		Definitive way of testing whether or not the object handles the url
		@param url a string, e.g. "/home/d.woolford/", "sludtke@10.10.10.10:~/test/"
		'''
		raise NotImplementedException("Inheriting classes must supply this functionality")
	
	def get_items(self,url):
		'''
		Get a list of EMListItems for display in the browser
		@return a list of EMListItems
		'''
		raise NotImplementedException("Inheriting classes must supply this functionality")
	
	def get_emdata(self,full_path,idx=0):
		'''
		EMListItems call this function
		Get a fully loaded EMData.
		This is so the read routine is in a single location and makes the EMListItems generic
		'''
		raise NotImplementedException("Inheriting classes must supply this functionality")
		
	def get_metadata(self,full_path,idx=0):
		'''
		EMListItems call this function
		Get the metadata, typically EMData header information. Must return a dict, or similar
		This is so the read routine is in a single location and makes the EMListItems generic
		'''
		raise NotImplementedException("Inheriting classes must supply this functionality")
	
	def get_stack_data(self,full_path):
		'''
		EMListItems call this function
		Might return a cache, of if that's not possible, a list of EMData objects, for example.
		The return value should 'look' like a list
		'''
		#from emimagemx import EMLightWeightParticleCache
		#return EMLightWeightParticleCache.from_file(full_path)
		raise NotImplementedException("Inheriting classes must supply this functionality")
	
	def url_mod_time(self,url):
		'''
		Get the last time the url was modified
		eg. return os.stat(url)[-2]
		May return None to indicate the call is not valid/supported - this will mean the corresponding list widget will not be automatically updated
		'''
		raise NotImplementedException("Inheriting classes must supply this functionality")
	
class EMFileSystemReader(EMBrowseDelegate):
	'''
	The Delegate for reading file system contents 
	Main function is to return the the correct list items to the EMSelector 
	'''
	def __init__(self,target):
		self.target = weakref.ref(target)
		self.threed_dim_limit = 128
		
	def get_emdata(self,full_path,idx=0):
		'''
		All items delegate by this should call this function to get a fully loaded EMData
		That way the read routine is in the one location
		'''
		e = EMData()
		e.read_image(full_path,idx)
		return e
		
	def get_metadata(self,full_path,idx=0):
		'''
		All items that load metadata using EMData io call this function
		'''
		e = EMData()
		e.read_image(full_path,idx, True)
		return e.get_attr_dict()
	
	def get_stack_data(self,full_path):
		'''
		This function is called by EM2DStackItem
		Return is a cache that can be treated like a list of EMData objects
		Can give speed ups
		'''
		from emimagemx import EMLightWeightParticleCache
		return EMLightWeightParticleCache.from_file(full_path)
	
	def url_mod_time(self,url):
		'''
		Get the last time the url was modified
		May return None to indicate the call is not valid/supported - this will mean the corresponding list widget will not be automatically updated
		'''
		if self.handles_url(url): return os.stat(url)[-2]
		else: return None
	
	def handles_url(self,url):
		'''
		Expected interface Delegate::handles_url
		'''
		if os.path.exists(url):
			if url.endswith("EMAN2DB") or url.endswith("EMAN2DB/"): return False # File System reader does not handle directories that begin with EMAN2DB
			else: return True
		
		return False

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
	
	def get_items(self,url):
		'''
		'''
		dirs, files = get_files_and_directories(url)
		if len(files) == 0 and len(dirs) == 0:
			# something is unusual about the directory, for instance, it is a file
			return []
		
		dirs.sort()
		files.sort()
		 
		return_items = []
		
		for i in dirs:
			if i[0] == '.': continue

			d = remove_directories_from_name(i)
			
			if d == "EMAN2DB": return_items.append(EMBDBFolderItem(self,"bdb",url+"/EMAN2DB")) # the EMBDBReader will know how handle this url
			else: return_items.append(EMFSFolderItem(self,i,url+"/"+i))
			
		for file in files:
			a = self.__get_file_item(file,url)
			if a != None: return_items.append(a)
	
		return return_items
	
	def __get_file_item(self,file,directory):
		'''
		Called internally
		@returns the correct item for the given file in the given directory, returns None in exceptional circumstances
		'''
		dtag = get_dtag()
		filt = self.target().get_file_filter()
		if file[0] == '.': return None
		if file[-1] == '~': return None
		
		extension = Util.get_filename_ext(file)
		full_name = directory+dtag+file
		
		e = EMData()
		item = None
		# note, if this if statement is allowed to proceed on Windows in the case of a png then the program
		# crashes. In December of 2008 I thus changed this if statement to automatically exclude unecessary files
		# such as pngs and jpges...etc.
		if EMUtil.get_image_ext_type(extension) != IMAGE_UNKNOWN and extension not in ["png","jpeg","jpg","JPG"]:
			try:
				e.read_image(full_name,0, read_header_only)
				if EMUtil.get_image_count(full_name) > 1:
					if e.get_zsize() > 1: item = EM3DStackItem(self,file,full_name)
					else: item = EM2DStackItem(self,file,full_name)
				else:
					if e.get_zsize() > 1: item = EM3DImageItem(self,file,full_name)
					else: item = EM2DImageItem(self,file,full_name)
			except:
				item = EMGenericItem(self,file,file)
		
		elif EMPlot2DModule.is_file_readable(full_name):
			item = EMFSPlotItem(self,file,full_name)
		else:
			if filt != "EM types":
				item = EMGenericItem(self,file,file)

		return item

class EMListItem(QtGui.QListWidgetItem):
	'''
	Base class definition providing the pubic interface of list widget items as 
	required by the EMSelector
	'''
	ICON = None
	def __init__(self,delegate=None,text=""):
		'''
		@param delegate an instance of an EMBrowseDelegate - a strong reference is made to this
		@param text the string that will be displayed in the QtGui.QListWidgetItem
		'''
		QtGui.QListWidgetItem.__init__(self,self.get_icon(),text)
		self.delegate = delegate
		self.context_menu_options = {} # this is used for running context menu actions
		self.icon = None
		self.metadata = None # subclass objects can use this to cache metadata
		self.data = None # subclassing objects can use this to cache data - actually this isn't used at the momemt for fear of memory hogging
		
	def get_delegate(self): return self.delegate
	
	def get_name(self):
		'''	Must return a unique name'''
		raise NotImplementedError("Inheriting classes must supply this function (%s)" %self)
	
	def get_icon(self):
		'''Supply your own Icon'''
		if EMListItem.ICON == None:
			EMListItem.ICON = QtGui.QIcon(get_image_directory() + "/File.png")
		return EMListItem.ICON
	
	def do_preview(self,emselector):
		'''
		The item must call the correct the function in the EMSelector to display itself.
		I contemplated doing it the other way (EMSelector instead knows how to display every item) but
		decided against thinking it would be too finicky to manage
		return True to indicate that displaying data was successful, False otherwise
		'''
		return False
	
	def is_previewable(self):
		'''
		Essentially asking if the subclassing object provides a definition
		of the do_preview function (that also does something useful, such as launch a display module)
		''' 
		return False

	def get_attr_dict(self):
		'''
		A wrapper for get_metadata - originally added to allow
		a generic interface in EMStackSaveDialog
		Unfortunately needs to stay for the time being
		'''
		return self.get_metadata()
	
	def get_url(self):
		'''	Returns the file name path, for BDB stuff this should in full bdb syntax'''
		return None
	
	def get_emdata(self):
		'''
		supposed to return an EMData object, if possible.
		Is customized to return a list of EMDatas, or a Cache, in some instances
		'''
		return None
	
	def is_2d_stack(self):
		'''Ask if a particular item is a 2D stack'''
		return False
	
	def is_2d_plot(self):
		'''Ask if a particular item is a 2D plot'''
		return False
	
	def get_metadata(self):
		'''
		Returns a dictionary if defined
		'''
		return None
	
	def loads_metadata(self):
		'''
		Ask if the item is capable of loading metadata
		This means the item can supply its own list of ListItems detailing metadata info
		See get_metadata_items
		'''
		return self.get_metadata() != None
	
	def get_metadata_items(self):
		'''
		Get the items that will be displayed if the item is capable if supplying metadata
		'''
		
		data = self.get_metadata()
		if data != None:
			return_items = []
			keys = data.keys()
			keys.sort() # alphabetical order
			for k in keys:
				v = data[k]
				return_items.append(EMKeyValueItem(self,str(k),k,v))
			
			return return_items
		return None
	
class EMUpArrowItem(EMListItem):
	ICON = None
	NAME = "go up a directory"
	def get_icon(self):
		'''
		Supply your own Icon
		'''
		if EMUpArrowItem.ICON == None:
			EMUpArrowItem.ICON = QtGui.QIcon(get_image_directory() + "/up_arrow.png")
		return EMUpArrowItem.ICON
	
	def get_name(self): return EMUpArrowItem.NAME
	
class EMDataListItem(EMListItem):
	'''
	Objects of this type are items that list EMData types -
	they happen to have the fact that they can be deleted and saved in common, that is all
	ABSTRACT
	'''
	def __init__(self,delegate=None,text="",full_name=""):
		EMListItem.__init__(self,delegate,text)
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
	
	def get_url(self):
		return self.full_path
	
class EM2DStackItem(EMDataListItem):
	ICON = None
	NAME = "2D stack"
	def get_name(self): return EM2DStackItem.NAME
	
	def do_preview(self,target):
		target.preview_data(None,self.full_path)
		return True
	
	def get_icon(self):
		'''Supply your own Icon	'''
		if EM2DStackItem.ICON == None:
			EM2DStackItem.ICON = QtGui.QIcon(get_image_directory() + "/multiple_images.png")
		return EM2DStackItem.ICON
	
	def is_previewable(self): return True
	
	def get_emdata(self):
		'''This one returns a list'''
		return self.delegate.get_stack_data(self.full_path)

	def is_2d_stack(self): return True
	
	def loads_metadata(self): return True
	
	def get_metadata_items(self):
		return [EM2DImageItem(self.delegate,str(i),self.full_path,i) for i in xrange(0,EMUtil.get_image_count(self.full_path))]
	
class EM3DStackItem(EMDataListItem):
	ICON = None
	NAME = "3D stack"
	def get_name(self): return EM3DStackItem.NAME
	
	def get_icon(self):
		'''
		Supply your own Icon
		'''
		if EM3DStackItem.ICON == None:
			EM3DStackItem.ICON = QtGui.QIcon(get_image_directory() + "/multiple_images_3d.png")
		return EM3DStackItem.ICON
	# no preview for this item as of Feb 2009
	def get_emdata(self):
		return self.delegate.get_stack_data(self.full_path)
	
	def loads_metadata(self): return True
	
	def is_previewable(self): return True # this is a hack - object's that are previewable should supply a useful do_preview function. But this doesn't actually break anything. But it's a hack :(
		
	def get_metadata_items(self):
		return [EM3DImageItem(self.delegate,str(i),self.full_path,i) for i in xrange(0,EMUtil.get_image_count(self.full_path))]

class EM2DImageItem(EMDataListItem):
	ICON = None
	NAME = "2D image"
	def __init__(self,delegate=None,text="",full_name="",idx=0):
		'''
		Have to supply init because we need the idx
		'''
		EMDataListItem.__init__(self,delegate,text,full_name)
		self.idx = idx
	
	def get_name(self): return EM2DImageItem.NAME
	
	def get_icon(self):
		'''
		Supply your own Icon
		'''
		if EM2DImageItem.ICON == None:
			EM2DImageItem.ICON = QtGui.QIcon(get_image_directory() + "/single_image.png")
		return EM2DImageItem.ICON

	def do_preview(self,target): 
		target.preview_data(self.get_emdata(),self.full_path)
		return True
		
	def is_previewable(self): return True

	def get_metadata(self):
		if self.metadata == None:
			self.metadata = self.delegate.get_metadata(self.full_path,self.idx)
		return self.metadata
	
	def get_emdata(self):
		return self.delegate.get_emdata(self.full_path,self.idx)

class EM3DImageItem(EM2DImageItem):
	ICON = None
	NAME = "3D image"
	def get_name(self): return EM3DImageItem.NAME
	def get_icon(self):
		'''
		Supply your own Icon
		'''
		if EM3DImageItem.ICON == None:
			EM3DImageItem.ICON = QtGui.QIcon(get_image_directory() + "/single_image_3d.png")
		return EM3DImageItem.ICON

class EMFSPlotItem(EMListItem):
	ICON = None
	NAME = "fs plot"
	
	def __init__(self,delegate=None,text="",full_name=""):
		EMListItem.__init__(self,delegate,text)
		self.full_path = full_name
	
	def get_name(self): return EMFSPlotItem.NAME
	
	def get_icon(self):
		'''
		Supply your own Icon
		'''
		if EMFSPlotItem.ICON == None:
			EMFSPlotItem.ICON = QtGui.QIcon(get_image_directory() + "/plot.png")
		return EMFSPlotItem.ICON
	
	def do_preview(self,target):
		target.preview_plot(self.full_path)
		return True
	
	def is_2d_plot(self): return True
		
	def is_previewable(self): return True
	
	def get_url(self): return self.full_path
	
class EMFSFolderItem(EMListItem):
	ICON = None
	NAME = "fs folder"
	
	def __init__(self,delegate=None,text="",full_name=""):
		EMListItem.__init__(self,delegate,text)
		self.full_path = full_name
	
	def get_name(self): return EMFSFolderItem.NAME
	
	def get_icon(self):
		'''
		Supply your own Icon
		'''
		if EMFSFolderItem.ICON == None:
			EMFSFolderItem.ICON = QtGui.QIcon(get_image_directory() + "/Folder.png")
		return EMFSFolderItem.ICON
	
	def get_url(self): return self.full_path

class EMBDBReader(EMBrowseDelegate):
	def __init__(self,target):
		self.target = weakref.ref(target)
		self.directory_replacements = {"EMAN2DB":"bdb"}
	
	def get_emdata(self,full_path,idx=0):
		'''
		All items delegate by this should call this function to get a fully loaded EMData
		That way the read routine is in the one location
		'''
		e = EMData()
		e.read_image(full_path,idx)
		return e
		
	def get_metadata(self,full_path,idx=0):
		'''
		All items that load metadata using EMData io call this function
		'''
		db_name =  full_path
		db = db_open_dict(db_name,ro=True)
		data = db.get_header(idx)
		return data
	
	def get_stack_data(self,full_path):
		'''
		This function is called by EM2DStackItem
		Return is a cache that can be treated like a list of EMData objects
		'''
		from emimagemx import EMLightWeightParticleCache
		return EMLightWeightParticleCache.from_file(full_path)
	
	def url_mod_time(self,url):
		'''
		Get the last time the url was modified
		May return None to indicate the call is not valid/supported - this will mean the corresponding list widget will not be automatically updated
		'''
		if self.handles_url(url): return os.stat(url)[-2]
		else: return None
	
	def handles_url(self,url):
		
		if url.endswith("EMAN2DB") or url.endswith("EMAN2DB/"): return True		
		return False
		
	
	def get_items(self,url):
#		if self.__is_database_file(url): 
#			list_widget.clear()
#			return self.__get_database_data_items(url)
		if self.handles_url(url): #os.path.exists(url) and (url.endswith("EMAN2DB") or url.endswith("EMAN2DB/")):
			return self.__get_bdb_directory_items(url)
		else: raise RuntimeError("Unknown url %s" %url)
	
	
#	def __get_database_data_items(self,url):
#		print "getting database items"
#		
#		if not self.__is_database_file(url): 
#			raise RuntimeError("Unknown url type %s" %url)
#		
#		file = self.__convert_to_absolute_path(url)
# 
#		db_directory = self.__get_emdatabase_directory(file)
#
#		key = remove_directories_from_name(file)
#		key = strip_file_tag(key)
#		
#		db_name = "bdb:"+db_directory+"#"+key
#		db = db_open_dict(db_name,ro=True)
#		
#		#items = DB[key] # NOTE items should be called "db" or something else
#		items = db
#		keys = items.keys()
#		keys.sort() # puts them alphabetical order
#		
#		return_items = []
#		for k in keys:
#			if k == '': continue
#			_type =db.item_type(k)
#			a = None
#			if _type == dict: a = EMBDBDictItem(self,str(k),db_directory,str(k))
#			elif _type == EMData:
#				data = db.get_header(k)
#				if data["nz"] > 1: a = EM3DImageItem(self,str(k), "bdb:"+db_directory+"#"+key,k)
#				else: a = EM2DImageItem(self,str(k), "bdb:"+db_directory+"#"+key,k)
#			elif _type == list and len(db[k]) == 2:
#				try:
#					if (isinstance(db[k][0][0],float) or isinstance(db[k][0][0],int)) and (isinstance(db[k][1][0],float) or isinstance(db[k][0][0],int)):
#						v = db[k]
#						a = EMBDBPlotItem(self,str(k)+":"+str(v),db_directory,key,k)
#					else: 
#						pass
#				except:
#					v = db[k]
#					a = EMKeyValueItem(self,str(k),k,v)
#			else:
#				v = db[k]
#				a = EMKeyValueItem(self,str(k),k,v)
#
#			if a != None: return_items.append(a)
#		
#		return return_items
	
	def __get_bdb_directory_items(self,url):
		
		'''
		Displays the file/folder information in the directory /home/someonone/data/EMAN2DB
		this will typically consist of .bdb (database) files, but may contain folders and other
		EMAN2DB directories.
		
		At the moment I have only written the code so that it supports the interrogation of the .bdb
		files, and am displaying the other folders only as a I reminder that they need to be dealt with
		'''
		#if not (os.path.exists(url) and (url.endswith("EMAN2DB") or url.endswith("EMAN2DB/"))):\
		if not self.handles_url(url): raise RuntimeError("Unknown url %s" %url)

		dtag = get_dtag()
		#real_directory = self.__convert_to_absolute_path(url)
		real_directory = url
		dirs,files = get_files_and_directories(real_directory)
		files.sort()
		dirs.sort()
		
		return_items = []
		for i in dirs:
			if i[0] == '.': continue
			
			if i == "EMAN2DB":
				b = EMGenericItem(self,"bdb","unwanted")  # really haven't accommodated for this...
				continue

			a = EMBDBDirectoryItem(self,i,url+"/"+i)
			return_items.append(a)
			
		for file in files:
			a = self.__get_bdb_file_item(file,real_directory)
			if a != None: return_items.append(a)
		return return_items

	def __get_bdb_file_item(self,file,real_directory):
		if not file[-3:] == "bdb": return None
		f = file.rpartition(".bdb")
		db_directory = self.__get_emdatabase_directory(real_directory)
		
		db_name = "bdb:"+db_directory+"#"+f[0]
		db = db_open_dict(db_name,ro=True)
		
		try:
			db.has_key("maxrec")
		except:
			# sometimes when the browser is updating in real time a database file is 
			# created, however only one of the two files exists (one is ptcl.bdb,
			# the other something like ptcl_200x200x1.bdb (etc), even though the other
			# is just about to be written... so I wait for 2 seconds and try a second time
			import time
			time.sleep(1)
			db = db_open_dict(db_name,ro=True)
			try:
				db.has_key("maxrec")
			except:
#					from emsprworkflow import EMErrorMessageDisplay
#					EMErrorMessageDisplay.run(["Warning: the %s database might be corrupted." %db_name], "Data loss" )
					
				return False
		
		if db and len(db) > 0:
			#n = DB[f[0]]["maxrec"]
			n = len(db)
			if n >= 1:
				d = db.get_header(0)
				if d["nz"] == 1: a = EM2DStackItem(self,f[0],"bdb:"+db_directory+"#"+f[0])
				elif d["nz"] > 1: a = EM3DStackItem(self,f[0],"bdb:"+db_directory+"#"+f[0])
			elif n == 0:
				d = db.get_header(0)
				if d["nz"] <= 1: a = EM2DImageItem(self,f[0], "bdb:"+db_directory+"#"+f[0],0)
				else: a = EM3DImageItem(self,f[0], "bdb:"+db_directory+"#"+f[0],0)
			else: a = EMBDBDictItem(self,f[0],db_directory,f[0])
		else: a = EMBDBDictItem(self,f[0],db_directory,f[0])
			
		a.file_name = file
		return a

#	def __convert_to_absolute_path(self,file_or_folder):
#		print "converting",file_or_folder,
#		dtag = get_dtag()
#		ret = file_or_folder
#		found = False
#		for dir_rep in self.directory_replacements.items():
#			if ret.find(dtag+dir_rep[1]) != -1:
#				ret = ret.replace(dtag+dir_rep[1],dtag+dir_rep[0])
#				found = True
#		if not found: return ret
#		if (not os.path.isdir(ret)) and (not os.path.isfile(ret)):
#			if ret[-1] == dtag: ret = ret[:-1]
#			ret += ".bdb"
#		
#		print ret
#		return ret
#
#	def __is_database_file(self,file_name):
#		file = self.__convert_to_absolute_path(file_name)
#		if len(file) > 4 and file[-4:] == ".bdb":
#			if self.__get_last_directory(file) == "EMAN2DB":
#				if file_exists(file):
#					return True
#			
#		return False

#	def __get_last_directory(self,file):
#		dtag = get_dtag()
#		idx1 = file.rfind(dtag)
#		if idx1 > 0:
#			ret = file[0:idx1]
#		else: return ret
#		
#		idx2 = ret.rfind(dtag)
#		if idx2 > 0:
#			ret = ret[idx2+1:]
#		
#		return ret
		
	def __get_emdatabase_directory(self,file):
		'''
		Get the database where EMAN2DB should be opening in order to open the given file
		e.g. if db path is /home/someone/work/EMAN2DB/data.bdb will return /home/someone/work
		'''
		idx1 = file.find("EMAN2DB")
		if idx1 > 0:
			return file[0:idx1-1]
		else: return None

class EMBDBEntryItem(EMListItem):
	'''
	Base items for BDB entries
	Stores the database directory and the name of the database
	'''
	def __init__(self,delegate=None,text="",db_directory="",db=""):
		EMListItem.__init__(self,delegate,text)
		self.database_directory = db_directory
		self.database = db
		
	def get_url(self):
		return "bdb:"+self.database_directory+"#"+self.database
	
class EMBDBPlotItem(EMBDBEntryItem):
	ICON = None
	NAME = "db plot"
	def __init__(self,delegate,text,db_directory,db,k):
		EMBDBEntryItem.__init__(self,delegate,text,db_directory,db)
#		self.database_directory = db_directory
#		self.database = db
		self.db_key = k
		
	def get_name(self): return EMBDBPlotItem.NAME
	
	def get_icon(self):
		'''
		Supply your own Icon
		'''
		if EMBDBPlotItem.ICON == None:
			EMBDBPlotItem.ICON = QtGui.QIcon(get_image_directory() + "/plot.png")
		return EMBDBPlotItem.ICON
	
		
	def do_preview(self,target):
		db_name = "bdb:"+self.database_directory+"#"+self.database
		db = db_open_dict(db_name)
		plot = db[self.db_key]
		target.preview_plot_list(self.db_key,plot)
		return True
	
	def is_previewable(self): return True
	
class EMBDBFolderItem(EMListItem):
	ICON = None
	NAME = "db folder"
	def __init__(self,delegate=None,text="",real_directory=""):
		EMListItem.__init__(self,delegate,text)
		self.full_path = real_directory
	
	def get_name(self): return EMBDBFolderItem.NAME
	
	def get_icon(self):
		'''
		Supply your own Icon
		'''
		if EMBDBFolderItem.ICON == None:
			EMBDBFolderItem.ICON = QtGui.QIcon(get_image_directory() + "/database.png")

		return EMBDBFolderItem.ICON
	
	def get_url(self):
		return self.full_path
		
class EMKeyValueItem(EMListItem):
	def __init__(self,delegate=None,text="",k="",v=None):
		EMListItem.__init__(self,delegate,text)
		self.key = k
		self.value = v
		
	def get_name(self): return "key_value"

	def loads_metadata(self): return True
	
	def get_metadata_items(self):
		return [ EMGenericItem(self.delegate,str(self.value)) ]

class EMBDBDirectoryItem(EMListItem):
	NAME = "bdb directory"
	ICON = None
	def __init__(self,delegate,text,real_directory):
		EMListItem.__init__(self,delegate,text)
		self.full_path = real_directory
	
	def get_name(self): return EMBDBDirectoryItem.NAME
	
	def get_icon(self):
		'''
		Supply your own Icon
		'''
		if EMBDBDirectoryItem.ICON == None:
			EMBDBDirectoryItem.ICON = QtGui.QIcon(get_image_directory() + "/Folder.png")

		return EMBDBDirectoryItem.ICON
	
	def get_url(self): return self.full_path

class EMBDBDictItem(EMBDBEntryItem):
	ICON = None
	NAME = "bdb dict"
	def __init__(self,delegate,text,db_directory,db):
		EMBDBEntryItem.__init__(self,delegate,text,db_directory,db)

	def get_name(self): return EMBDBDictItem.NAME
	
	def get_icon(self):
		'''
		Supply your own Icon
		'''
		if EMBDBDictItem.ICON == None:
			EMBDBDictItem.ICON = QtGui.QIcon(get_image_directory() + "/database.png")

		return EMBDBDictItem.ICON
	
	def loads_metadata(self): return True
	
	def get_metadata(self):
		if self.metadata == None:
			self.metadata = db_open_dict("bdb:"+self.database_directory+"#"+self.database,ro=True)
		return self.metadata
	
class EMGenericItem(EMListItem):
	'''
	A dead end item, displays a value. Has no metadata and there is no consequence for clicking on it. 
	'''
	NAME = "generic"
	def __init__(self,delegate=None,text="",key=None):
		EMListItem.__init__(self,delegate,text)
		self.key = key
		
	def do_preview(self,target): pass
	
	def get_name(self): return EMGenericItem.NAME

class EMBrowserDialog(object):
	def __new__(self,target):
		selector = EMSelectorTemplate(QtGui.QWidget)(target,False)
		selector.setWindowTitle("EMAN2 Browser")
		selector.preview_options.setCurrentIndex(1)
		selector.preview_options_changed(selector.preview_options.currentText())
		selector.ok_button.setEnabled(False)
		selector.cancel_button.setEnabled(False)
		return selector

class EMBrowserModule(EMQtWidgetModule):
	def __init__(self):
		self.widget = EMBrowserDialog(self)
		EMQtWidgetModule.__init__(self,self.widget)

	

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
	#dialog = EMSelector(None,em_app)
	em_qt_widget = EMSelectorModule()
	QtCore.QObject.connect(em_qt_widget.emitter(),QtCore.SIGNAL("ok"),on_done)
	QtCore.QObject.connect(em_qt_widget.emitter(),QtCore.SIGNAL("cancel"),on_cancel)
	em_app.show()
	em_app.execute()


