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
from EMAN2 import EMData,Region, gimme_image_dimensions3D, EMFunctor, get_file_tag, name_has_no_tag, remove_directories_from_name, file_exists, strip_file_tag, EMData
from emimage2d import EMImage2DModule
from emapplication import EMStandAloneApplication, EMQtWidgetModule
from EMAN2db import EMAN2DB
#from boxertools import TrimSwarmAutoBoxer


class EMSelectorDialog(QtGui.QDialog):
	def __init__(self,target,application):
		QtGui.QDialog.__init__(self,None)
		self.setFocusPolicy(Qt.StrongFocus)
		self.application=application
		self.target=target
		
		self.db_listing = EMBDBListing(self)
		self.dir_listing = EMDirectoryListing(self)
		
		self.hbl = QtGui.QVBoxLayout(self)
		self.hbl.setMargin(0)
		self.hbl.setSpacing(6)
		self.hbl.setObjectName("hbl")
	
		self.list_hdl = QtGui.QHBoxLayout()
		
		self.__init_icons()
		
		self.__init_filter_combo()
		
		self.first_list_widget = QtGui.QListWidget(None)
		self.starting_directory = os.getcwd()
		
		self.selections = []
		self.current_list_widget = None
		self.lock = True
		self.list_widgets = []
		self.list_widget_data= [] # entries should be tuples containing (current folder item)
		self.__add_list_widget(self.first_list_widget)
		self.__add_list_widget()
		self.__add_list_widget()
		
		self.__load_directory_data(self.starting_directory,self.first_list_widget)
		#self.first_list_widget.setCurrentRow(-1)
		
		self.hbl.addLayout(self.list_hdl)

		self.bottom_hbl = self.list_hdl = QtGui.QHBoxLayout()
		self.bottom_hbl.addWidget(self.filter_combo,1)
		self.__init_done_button()
		self.bottom_hbl.addWidget(self.done_button,0)
		self.__init__single_preview_tb()
		self.bottom_hbl.addWidget(self.single_preview,0)
		self.hbl.addLayout(self.bottom_hbl)
		self.gl_image_preview = None
		
		self.resize(480,480)
		
		self.lock = False
		
		self.paint_events = 0
		
		
		
	def get_desktop_hint(self):
		return "dialog"
		
	def set_application(self,app):
		self.application = app
		
	def __init_icons(self):
		self.setWindowIcon(QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/eman.png"))
		self.folder_icon = QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/Folder.png")
		self.folder_files_icon = QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/FolderFiles.png")
		self.file_icon = QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/File.png")
		self.database_icon = QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/database.png")
		self.key_icon = QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/Key.png")
		self.basic_python_icon = QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/boxhabanosclose.png")
		self.dict_python_icon = QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/Bag.png")
		self.ab_refboxes_icon = QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/black_box.png")
		self.ab_manboxes_icon = QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/black_box.png")
		self.ab_autoboxes_icon = QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/green_boxes.png")
		self.emdata_icon = QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/boxer_erase.png");
		
	
	def __init__single_preview_tb(self):
		self.single_preview = QtGui.QCheckBox("Single preview")
		self.single_preview.setChecked(True)
		
		QtCore.QObject.connect(self.single_preview, QtCore.SIGNAL("clicked(bool)"),self.single_preview_clicked)
	
	def __init_done_button(self):
		self.done_button = QtGui.QPushButton("Done")
		self.done_button.adjustSize()
	
		QtCore.QObject.connect(self.done_button, QtCore.SIGNAL("clicked(bool)"),self.done_button_clicked)
	
	def single_preview_clicked(self,bool):
		pass
		#print "not supported"
	
	def done_button_clicked(self,bool):
		self.emit(QtCore.SIGNAL("done"),self.selections)
		
	def __init_filter_combo(self):
		self.filter_combo = QtGui.QComboBox(None)
		self.filter_combo.addItem("*.mrc,*.hdf,*.img")
		self.filter_combo.addItem("Databases")
		self.filter_combo.addItem("*.*")
		self.filter_combo.setEditable(True)
	
		QtCore.QObject.connect(self.filter_combo, QtCore.SIGNAL("currentIndexChanged(int)"),self.filter_index_changed)
		QtCore.QObject.connect(self.filter_combo, QtCore.SIGNAL("currentIndexChanged(QString&)"),self.filter_index_changed)

	def filter_index_changed(self):
		self.__redo_list_widget_contents()
	
	def __redo_list_widget_contents(self):
		self.lock = True
		directory = self.starting_directory+"/"
		for i,data in  enumerate(self.list_widget_data):
			
			if data != None:d = str(data.text())
			old_row = self.list_widgets[i].currentRow()
			self.__load_directory_data(directory,self.list_widgets[i])
			self.list_widget_data[i] = self.list_widgets[i].item(old_row)
			if data == None: return
			else:
				directory += '/' + d
	
		self.lock = False
	def __add_list_widget(self, list_widget = None):
		if list_widget == None:	list_widget = QtGui.QListWidget(None)
		
		list_widget.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
		list_widget.setMouseTracking(True)	
		self.list_widgets.append(list_widget)
		self.list_hdl.addWidget(list_widget)
		self.list_widget_data.append(None)
		
		QtCore.QObject.connect(list_widget, QtCore.SIGNAL("itemDoubleClicked(QListWidgetItem*)"),self.list_widget_dclicked)
		QtCore.QObject.connect(list_widget, QtCore.SIGNAL("itemPressed(QListWidgetItem*)"),self.list_widget_clicked)
		QtCore.QObject.connect(list_widget, QtCore.SIGNAL("currentRowChanged (int)"),self.list_widget_row_changed)
		#QtCore.QObject.connect(list_widget, QtCore.SIGNAL("paintEvent (int)"),self.list_widget_row_changed)
		QtCore.QObject.connect(list_widget, QtCore.SIGNAL("itemEntered(QListWidgetItem*)"),self.list_widget_item_entered)
		QtCore.QObject.connect(list_widget, QtCore.SIGNAL("currentItemChanged(QListWidgetItem*,QListWidgetItem*)"),self.list_widget_current_changed)
		QtCore.QObject.connect(list_widget, QtCore.SIGNAL("itemChanged(QListWidgetItem*)"),self.list_widget_item_changed)
		QtCore.QObject.connect(list_widget, QtCore.SIGNAL("itemActivated(QListWidgetItem*)"),self.list_widget_item_activated)
		#QtCore.QObject.connect(list_widget, QtCore.SIGNAL("activated(QModelIndex)"),self.activated)
		QtCore.QObject.connect(list_widget, QtCore.SIGNAL("selectionChanged(QItemSelection,QItemSelection)"),self.selection_changed)
		
	def __go_back_a_directory(self):
		self.lock = True
		new_dir = self.starting_directory[0:self.starting_directory.rfind('/')]
		#if  self.db_listing.responsible_for(new_dir):
			#new_dir = self.db_listing.
		#print "going back a directory",new_dir
		#if not os.access(new_dir,os.R_OK):
			#print "can't go up a directory, don't have read permission"
			#return
		self.starting_directory = new_dir
		for j in range(0,len(self.list_widgets)):
			self.list_widgets[j].clear()
			self.list_widget_data[j] = None
		self.__load_directory_data(self.starting_directory,self.first_list_widget)
		#self.hide_preview()
		self.lock = False
		
	def __go_forward_a_directory(self):
		self.lock = True
		self.starting_directory = self.starting_directory + '/' + str(self.list_widget_data[0].text())
		
		directory = self.starting_directory 
		for i in range(len(self.list_widgets)-1):
			items = []
			old_row = self.list_widgets[i+1].currentRow()
			n = self.list_widgets[i+1].count()
			for j in range(n-1,-1,-1):
				items.append(self.list_widgets[i+1].takeItem(j))
				
			self.list_widgets[i].clear()	
			for k in items:
				self.list_widgets[i].insertItem(0,k)
			
			self.list_widgets[i].setCurrentRow(old_row)
			
			self.list_widget_data[i] = self.list_widgets[i].item(old_row)
			directory += '/' + str(self.list_widget_data[i].text())
		
		self.list_widgets[0].insertItem(0,"../")
		self.lock = False
		#self.hide_preview()
		
	def selection_changed(self,item1,item2):
		pass
		
	def __update_selections(self):
		'''
		Makes the list of currently selected files accurate and up to date. Called when
		something has been clicked in a a list widget
		'''
		
		# get the directory 
		directory = self.starting_directory+"/"
		idx = 0
		for i,list_widget in enumerate(self.list_widgets):
			if list_widget == self.current_list_widget: 
				break
			directory += str(self.list_widget_data[i].text()) + "/"
		else:
			print "no list widget has focus?"
			return
		
		# now make the list of selections accurate
		self.selections = []
		for a in self.current_list_widget.selectedIndexes():
			file = directory+str(a.data().toString())
			if self.__is_previewable(file): self.selections.append(file)
		
		# if there are no selections then close the preview
		#if len(self.selections) == 0:
			#self.hide_preview()
	
	def hide_preview(self):
		if self.gl_image_preview  != None:
			self.application.hide_specific(self.gl_image_preview)
	
	def list_widget_item_activated(self,item):
		pass
		
	
	def list_widget_item_changed(self,item):
		pass
	
	def list_widget_current_changed(self,item1,item2):
		pass
		
	def list_widget_item_entered(self,item):
		self.current_list_widget = item.listWidget()

	def list_widget_row_changed(self,i):
		return
	
	def keyPressEvent(self,event):
		pass
	
	def list_widget_clicked(self,item):
		if self.lock : return
		if self.current_list_widget == None: return

		self.__update_selections()

		if item == None: return
		
		if item.text() == "../": 
			self.__go_back_a_directory()
			return
		
		file = self.starting_directory+"/"
		directory = self.starting_directory+"/"
		idx = 0
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


		if self.set_preview(file):
			for i in range(idx+1,len(self.list_widgets)):
				self.list_widgets[i].clear()
				self.list_widget_data[i] = None
			return
		
	
		n = len(self.list_widgets)-1
		if self.current_list_widget  == self.list_widgets[n]:
				self.list_widget_data[n] = item
				self.__go_forward_a_directory()
				self.__load_directory_data(file+'/',self.list_widgets[n])
				return

		if self.__load_directory_data(file,self.list_widgets[idx+1]):
			#if old_item != None:
				#old_item.setBackgroundColor(QtGui.QColor(255,255,255))
			##item.setBackgroundColor(QtGui.QColor(64,190,0,63))	
			self.list_widget_data[idx] = item
			self.list_widget_data[idx+1] = None
			
		for i in range(idx+2,len(self.list_widgets)):
			self.list_widgets[i].clear()
			self.list_widget_data[i] = None
		
	def list_widget_dclicked(self,item):
		
		file = self.starting_directory+"/"
		for i,list_widget in enumerate(self.list_widgets):
			if item.listWidget() == list_widget:
				file += str(item.text())
				break
			file += str(self.list_widget_data[i].text()) + "/"
	
		
		if self.__is_previewable(file):
			print "Opening ", file
	
	def set_preview(self,filename):
		
		if self.db_listing.do_preview(filename):
			return True
		elif self.dir_listing.do_preview(filename):
			return True
		
		return False
		
	def preview_data(self,a,filename=""):
		self.application.setOverrideCursor(Qt.BusyCursor)
		
		if type(a) == list and len(a) == 1:
			a = a[0]
			data = []
			if a.get_zsize() != 1:
				for z in range(a.get_zsize()):
					image = a.get_clip(Region(0,0,z,a.get_xsize(),a.get_ysize(),1))
					data.append(image)
				a = data
		
		if self.single_preview.isChecked():
			if self.gl_image_preview == None:
				self.gl_image_preview = EMImage2DModule(application=self.application)
				self.gl_image_preview.set_data(a,filename)
			else:
				self.gl_image_preview.set_data(a,filename,retain_current_settings=True)
				
			#self.gl_image_preview.set_file_name(f)
			self.application.show_specific(self.gl_image_preview)
			self.gl_image_preview.updateGL()
		else:
			preview = EMImage2DModule(application=self.application)
			preview.set_data(a,filename)
			self.application.show_specific(preview)
			preview.updateGL()
			
		self.application.setOverrideCursor(Qt.ArrowCursor)
		
	
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
	
	def __is_previewable(self,s):
		if  self.db_listing.responsible_for(s):
			if self.db_listing.is_previewable(s): return True 
			else: return False
		else: return self.dir_listing.is_previewable(s)
		
	def __is_non_empty_directory(self,s):
		'''
		Returns true if s is a non empty directory
		'''

		for root, dirs, files in os.walk(s):
			files = self.filter_strings(files)
			file_length = len(files)
			if file_length == 0: file_length = len(dirs)
		
			if file_length != 0: return True
			
			return False
	
	def __load_database_directory(self,database_name,list_widget):
		
		print "loading database directory",database_name
	
	
	def __load_directory_data(self,directory,list_widget):
		
		list_widget.clear()
		if (list_widget == self.list_widgets[0]):
			self.lock = True
			QtGui.QListWidgetItem("../",list_widget)
			self.lock = False
		
		if self.db_listing.responsible_for(directory):
			if  self.db_listing.load_database_data(directory,list_widget): return True
			elif self.db_listing.load_directory_data(directory,list_widget): return True
			elif self.db_listing.load_database_interrogation(directory,list_widget): return True
			else: return False
		else: return self.dir_listing.load_directory_data(directory,list_widget)
		
	def make_replacements(self,dirs,list_widget):
		self.db_listing.make_replacements(dirs,list_widget)


class EMDirectoryListing:
	def __init__(self,target):
		self.target = target
		pass
	
	def load_directory_data(self,directory,list_widget):
		for root, dirs, files in os.walk(directory):
			files = self.target.filter_strings(files)
			
			dirs.sort()
			files.sort()
			 
			self.target.make_replacements(dirs,list_widget)
			 
			for i in dirs:
				if i[0] == '.': continue
				
				file_length = 0
				for r, d, f in os.walk(directory+"/"+i):
					f = self.target.filter_strings(f)
					file_length = len(f)
					if file_length == 0: file_length = len(d)
					break
				if file_length != 0:
					a = QtGui.QListWidgetItem(self.target.folder_files_icon,i,list_widget)
				else:
					a = QtGui.QListWidgetItem(self.target.folder_icon,i,list_widget)
				
			#for i,file in enumerate(files):
				#if file[0] == '.': continue
				#if get_file_tag(file) == "bdb":
					#a = QtGui.QListWidgetItem(self.database_icon,file,list_widget)
					#files.pop(i)
					
			for file in files:
				a = QtGui.QListWidgetItem(self.target.file_icon,file,list_widget)

			return True
			
		return False
	
	def do_preview(self,file_or_folder):
		if not os.path.isfile(file_or_folder): return False
		
		filename = file_or_folder
		try: 
			dims = gimme_image_dimensions3D(filename)
			if dims[2] != 1:
				for d in dims: 
					if d > 128:
						print "3D image too big, no preview available"
						return False
			elif dims[1] != 1:
				for d in [dims[0],dims[1]]:
					if d > 1024:
						print "2D image too big, no preview available"
						return False
			a=EMData.read_images(filename)
		except: 
			return False
		
		self.target.preview_data(a,filename)
		return True

	def do_preview_in_browse(self,file_or_folder):
		if not os.path.isfile(file_or_folder): return False
		
		filename = file_or_folder
		try: 
			dims = gimme_image_dimensions3D(filename)
			if dims[2] != 1:
				for d in dims: 
					if d > 128:
						print "3D image too big, no preview available"
						return False
			elif dims[1] != 1:
				for d in [dims[0],dims[1]]:
					if d > 1024:
						print "2D image too big, no preview available"
						return False
			a=EMData.read_images(filename)
		except: 
			return False
		
		self.target.preview_data(a,filename)
		return True

	def is_previewable(self,file_or_folder):
		# this may not be strictly correct, seeing as if it's a file it will return true
		return os.path.isfile(file_or_folder)
		

class EMBDBListing:
	def __init__(self,target):
		self.target = target
		self.directory_replacements = {"EMAN2DB":"bdb"}
	
	def responsible_for(self,file_or_folder):
		real_name = self.convert_to_absolute_path(file_or_folder)
		split = real_name.split('/')
		split.reverse() # this probably makes things faster
		for s in split:
			if s in self.directory_replacements.keys() or (len(s) > 4 and s[-4:] == ".bdb"): return True 
	
		return False
	
	def make_replacements(self,dirs,list_widget):
		rm = []
		for i,directory in enumerate(dirs):
			d = remove_directories_from_name(directory)
			
			if d in self.directory_replacements.keys():
				a = QtGui.QListWidgetItem(self.target.database_icon,self.directory_replacements[d],list_widget)
				rm.append(i)
		
		rm.reverse()
		for r in rm:
			dirs.pop(r)
	
	def convert_to_absolute_path(self,file_or_folder):
		ret = file_or_folder
		for dir_rep in self.directory_replacements.items():
			ret = ret.replace('/'+dir_rep[1],'/'+dir_rep[0])
		if (not os.path.isdir(ret)) and (not os.path.isfile(ret)):
			if ret[-1] == "/": ret = ret[:-1]
			ret += ".bdb"
			
		return ret
			
	def is_database_directory(self,directory):
		if remove_directories_from_name(directory) in self.directory_replacements.values(): return True
		else: return False

	def load_directory_data(self,directory,list_widget):
		if not remove_directories_from_name(directory) in self.directory_replacements.values():
			 return False

		real_directory = self.convert_to_absolute_path(directory)
		for root, dirs, files in os.walk(real_directory):
			files.sort()
			dirs.sort()
			
			for i in dirs:
				if i[0] == '.': continue
				
				if i == "EMAN2DB":
					a = QtGui.QListWidgetItem(self.target.database_icon,"bdb",list_widget)
					continue
			
				file_length = 0
				for r, d, f in os.walk(real_directory+"/"+i):
					file_length = len(f)
					if file_length == 0: file_length = len(d)
					break
				if file_length != 0:
					a = QtGui.QListWidgetItem(self.target.folder_files_icon,i,list_widget)
				else:
					a = QtGui.QListWidgetItem(self.target.folder_icon,i,list_widget)
				
			for file in files:
				if file[len(file)-3:] == "bdb":
					f = file.rpartition(".bdb")
					a = QtGui.QListWidgetItem(self.target.database_icon,f[0],list_widget)
				#else:
					#a = QtGui.QListWidgetItem(self.target.key_icon,file,list_widget)
				
			return True
				
		return False

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
		DB = EMAN2DB.open_db(db_directory)
		
		
		#print "the database directory is",db_directory
		key = remove_directories_from_name(file)
		key = strip_file_tag(key)
		DB.open_dict(key)
		
		list_widget.clear()
		item = DB[key]
		keys = item.keys()
		keys.sort() # puts them alphabetical order
		for k in keys:
			i = item[k]
			if k == '': continue
			if type(i) == dict:
				a = QtGui.QListWidgetItem(self.target.database_icon,str(k),list_widget)
			elif type(i) == EMData:
				a = QtGui.QListWidgetItem(self.target.emdata_icon,str(k),list_widget)
			else:
				#if type(i) in [str,float,int,tuple,list,bool]:
				a = QtGui.QListWidgetItem(self.target.basic_python_icon,str(k),list_widget)	
		return True
				
	
	def get_last_directory(self,file):
		idx1 = file.rfind('/')
		if idx1 > 0:
			ret = file[0:idx1]
		else: return ret
		
		idx2 = ret.rfind('/')
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
			return file[0:idx1]
		else: return None
		
	
	def load_database_interrogation(self,file_name,list_widget):
		split = file_name.split('/')
		
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

			DB = EMAN2DB.open_db(db_directory)
			
			key = split[j-1]
			item_key = split[j-2]
			
			DB.open_dict(key)
			item = DB[key]
			
			#item = item[item_key]
			for ii in range(j-2,-1,-1):
				item = item[split[ii]]
			
			if type(item) == dict:
				keys = item.keys()
				keys.sort() # puts them alphabetical order
				for k in keys:
					i = item[k]
					if k == "auto_boxes":
						a = QtGui.QListWidgetItem(self.target.ab_autoboxes_icon,str(k),list_widget)
					elif k == "reference_boxes":
						a = QtGui.QListWidgetItem(self.target.ab_autoboxes_icon,str(k),list_widget)
					elif k == "manual_boxes":
						a = QtGui.QListWidgetItem(self.target.ab_autoboxes_icon,str(k),list_widget)
					elif type(i) in [str,float,int,tuple,list,bool]:
						a = QtGui.QListWidgetItem(self.target.basic_python_icon,str(k),list_widget)
					elif type(i) == dict:
						a = QtGui.QListWidgetItem(self.target.dict_python_icon,str(k),list_widget)
					elif type(i) == EMData:
						a = QtGui.QListWidgetItem(self.target.emdata_icon,str(k),list_widget)
					else:
						a = QtGui.QListWidgetItem(self.target.basic_python_icon,str(k),list_widget)
			elif isinstance(item,EMData):
				print "this shouldn't happen"
				self.target.preview_data(item)
				return False
			else:
				a = QtGui.QListWidgetItem(self.target.basic_python_icon,str(item),list_widget)
			
			return True
				
		else: return False 
	
	def is_previewable(self,file_name):
		return self.do_preview(file_name,fake_it=True)
	
	def do_preview(self,file_name,fake_it=False):
		split = file_name.split('/')
		
		rm = []
		for i,s in enumerate(split):
			if len(s) == 0: rm.append(i)
		
		rm.reverse()
		for k in rm: split.pop(k)
		
		if len(split) > 1 : # must atleast have EMAN2DB/something.bdb/dictionary
			split.reverse() # this probably makes things faster
			for j in range(1,len(split)):
				if split[j] in self.directory_replacements.values():
					break
			else: return False
			
			real_name = self.convert_to_absolute_path(file_name)
			db_directory = self.get_emdatabase_directory(real_name)
			DB = EMAN2DB.open_db(db_directory)
			
			key = split[j-1]
			#item_key = split[j-2]
			
			
			#print key,item_key
			DB.open_dict(key)
			item = DB[key]
			for ii in range(j-2,-1,-1):
				for t in [type(split[ii]),float,int]:
					try:
						key = t(split[ii])
						if item[key] != None: 
							item = item[key]
							break
					except:
						pass
			
			if isinstance(item,EMData):
				if not fake_it: self.target.preview_data(item)
				return True
			
		return False
			
	def load_database_variables(self,directory,list_widget):
		pass
		
app = None
def on_done(string_list):
	print "on done"
	if len(string_list) != 0:
		for s in string_list:
			print s,
		print
	app.quit()


if __name__ == '__main__':
	em_app = EMStandAloneApplication()
	app = em_app
	dialog = EMSelectorDialog(None,em_app)
	em_qt_widget = EMQtWidgetModule(dialog,em_app)
	QtCore.QObject.connect(dialog,QtCore.SIGNAL("done"),on_done)
	em_app.show()
	em_app.execute()




