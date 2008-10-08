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
from EMAN2 import EMData,Region, gimme_image_dimensions3D, EMFunctor
from emimage2d import EMImage2DModule
from emapplication import EMStandAloneApplication, EMQtWidgetModule

class EMSelectorDialog(QtGui.QDialog):
	def __init__(self,target,application):
		QtGui.QDialog.__init__(self,None)
		#self.setMouseTracking(True)
		self.application=application
		self.target=target
		
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
		self.hbl.addLayout(self.bottom_hbl)
		self.gl_image_preview = None
		
		self.resize(480,480)
		
		self.lock = False
		
		self.paint_events = 0
	def set_application(self,app):
		self.application = app
		
	def __init_icons(self):
		self.setWindowIcon(QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/eman.png"))
		self.folder_icon = QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/Folder.png")
		self.folder_files_icon = QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/FolderFiles.png")
		self.file_icon = QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/File.png")
		
	
	def __init_done_button(self):
		self.done_button = QtGui.QPushButton("Done")
		self.done_button.adjustSize()
	
		QtCore.QObject.connect(self.done_button, QtCore.SIGNAL("clicked(bool)"),self.done_button_clicked)
		
	def done_button_clicked(self,bool):
		self.emit(QtCore.SIGNAL("done"),self.selections)
		
	def __init_filter_combo(self):
		self.filter_combo = QtGui.QComboBox(None)
		self.filter_combo.addItem("*.mrc,*.hdf,*.img")
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
		QtCore.QObject.connect(list_widget, QtCore.SIGNAL("activated(QModelIndex)"),self.activated)
		QtCore.QObject.connect(list_widget, QtCore.SIGNAL("selectionChanged(QItemSelection,QItemSelection)"),self.selection_changed)
		
	def __go_back_a_directory(self):
		self.lock = True
		new_dir = self.starting_directory[0:self.starting_directory.rfind('/')]
		if not os.access(new_dir,os.R_OK): 
			print "can't go up a directory, don't have read permission"
			return
		self.starting_directory = new_dir
		for j in range(0,len(self.list_widgets)):
			self.list_widgets[j].clear()
			self.list_widget_data[j] = None
		self.__load_directory_data(self.starting_directory,self.first_list_widget)
		self.hide_preview()
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
		self.hide_preview()
	def selection_changed(self,item1,item2):
		pass

	def activated(self,index):
		directory = self.starting_directory+"/"
		idx = 0
		for i,list_widget in enumerate(self.list_widgets):
			if list_widget == self.current_list_widget: 
				break
			directory += str(self.list_widget_data[i].text()) + "/"
		else:
			print "no list widget has focus?"
			return
		
		self.selections = []
		for a in self.current_list_widget.selectedIndexes():
			file = directory+str(a.data().toString())
			if self.__is_file(file): self.selections.append( file)
			
		if len(self.selections) == 0:
			self.hide_preview()
	
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
		#print self.current_list_widget
	
	def list_widget_row_changed(self,i):
		return
	
	def list_widget_clicked(self,item):
		if self.lock : return
		if self.current_list_widget == None: return
		#if self.paint_events < 2 :
			 #self.current_list_widget.setCurrentRow(-1)
			 #return
		
		#item = self.current_list_widget.currentItem()
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

		
		n = len(self.list_widgets)-1
		if self.current_list_widget  == self.list_widgets[n] and not self.__is_file(file):
			self.list_widget_data[n] = item
			self.__go_forward_a_directory()
			self.__load_directory_data(file+'/',self.list_widgets[n])
			return
		
		

		# check to see if we have to go back a directory
	
		
		if self.__is_file(file):
			self.__set_preview(file)
			for i in range(idx+1,len(self.list_widgets)):
				self.list_widgets[i].clear()
				self.list_widget_data[i] = None
			return
	
		
		old_item = self.list_widget_data[idx]
		
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
	
		
		if self.__is_file(file):
			print "Opening ", file
	
	def __set_preview(self,filename):
		self.application.setOverrideCursor(Qt.BusyCursor)
		try: 
			dims = gimme_image_dimensions3D(filename)
			if dims[2] != 1:
				for d in dims: 
					if d > 128:
						print "3D image too big, no preview available"
						self.application.setOverrideCursor(Qt.ArrowCursor)
						return
			elif dims[1] != 1:
				for d in [dims[0],dims[1]]:
					if d > 1024:
						print "2D image too big, no preview available"
						self.application.setOverrideCursor(Qt.ArrowCursor)
						return
			a=EMData.read_images(filename)
		except: 
			self.application.setOverrideCursor(Qt.ArrowCursor)
			return
		
		if len(a) == 1:
			a = a[0]
			data = []
			if a.get_zsize() != 1:
				for z in range(a.get_zsize()):
					image = a.get_clip(Region(0,0,z,a.get_xsize(),a.get_ysize(),1))
					data.append(image)
				a = data
		
		if self.gl_image_preview == None:
			self.gl_image_preview = EMImage2DModule(application=self.application)
			
		self.gl_image_preview.set_data(a,filename)
		#self.gl_image_preview.set_file_name(f)
		self.application.show_specific(self.gl_image_preview)
		self.gl_image_preview.updateGL()
		
		self.application.setOverrideCursor(Qt.ArrowCursor)
		
		self.activateWindow()

	
	def __filter_strings(self,strings):
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
	
	def __is_file(self,s):
		for root, dirs, files in os.walk(s):
			return False
		
		return True
		

	def __is_non_empty_directory(self,s):
		'''
		Returns true if s is a non empty directory
		'''

		for root, dirs, files in os.walk(s):
			files = self.__filter_strings(files)
			file_length = len(files)
			if file_length == 0: file_length = len(dirs)
		
			if file_length != 0: return True
			
			return False
	
	
	def __load_directory_data(self,directory,list_widget):
		
		for root, dirs, files in os.walk(directory):
			files = self.__filter_strings(files)
			
			list_widget.clear()
			dirs.sort()
			files.sort()
			
			if (list_widget == self.list_widgets[0]):
				self.lock = True
				QtGui.QListWidgetItem("../",list_widget)
				self.lock = False
				 
			for i in dirs:
				if i[0] == '.': continue
				file_length = 0
				for r, d, f in os.walk(directory+"/"+i):
					f = self.__filter_strings(f)
					file_length = len(f)
					if file_length == 0: file_length = len(d)
					break
				if file_length != 0:
					a = QtGui.QListWidgetItem(self.folder_files_icon,i,list_widget)
				else:
					a = QtGui.QListWidgetItem(self.folder_icon,i,list_widget)
				
			for i in files:
				if i[0] == '.': continue
				a = QtGui.QListWidgetItem(self.file_icon,i,list_widget)

			return True
			
		return False

	#def enterEvent(self,event):
		#self.lock = False
		
		
	#def leaveEvent(self,event):
		#self.lock = True
		
	#def paintEvent(self,event):
		#self.paint_events += 1 # i had to do this to avoid unwanted behavior

app = None
def on_done(string_list):
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




