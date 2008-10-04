#!/usr/bin/env python

#
# Author: David Woolford (sludtke@bcm.edu)
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

from emapplication import EMStandAloneApplication, EMQtWidgetModule

class EMSelectorDialog(QtGui.QDialog):
	def __init__(self,target) :
		QtGui.QDialog.__init__(self,None)
		self.target=target
		
		self.hbl = QtGui.QVBoxLayout(self)
		self.hbl.setMargin(0)
		self.hbl.setSpacing(6)
		self.hbl.setObjectName("hbl")
		
		self.list_vbl = QtGui.QHBoxLayout()
		
		self.setWindowIcon(QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/eman.png"))
		
		self.folder_icon = QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/Folder.png")
		self.folder_files_icon = QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/FolderFiles.png")
		self.file_icon = QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/File.png")
		self.first_list_widget = QtGui.QListWidget(None)
		self.starting_directory = os.getcwd()
		
		
		#for i in self.starting_cwd_contents:
			#print i,os.stat(os.getcwd()+"/"+i)
			#folder_icon = QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/Folder.png")
			#file_icon = QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/File.png")
			#a = QtGui.QListWidgetItem(file_icon,i,self.list_view)
			##self.list_view.insertItem(a)

		self.list_widgets = []
		self.list_widget_data= [] # entries should be tuples containing (current folder item)
		self.__add_list_widget(self.first_list_widget)
		self.__add_list_widget()
		self.__add_list_widget()
		
		self.__load_directory_data(self.starting_directory,self.first_list_widget)
		
		self.hbl.addLayout(self.list_vbl)
	
	def mousePressEvent(self,event):
		print "pressed"
		
	def mouseReleaseEvent(self,event):
		print "released"
		
	def __add_list_widget(self, list_widget = None):
		if list_widget == None:	list_widget = QtGui.QListWidget(None)
			
		self.list_widgets.append(list_widget)
		self.list_vbl.addWidget(list_widget)
		self.list_widget_data.append(None)
		
		QtCore.QObject.connect(list_widget, QtCore.SIGNAL("itemDoubleClicked(QListWidgetItem*)"),self.list_widget_dclicked)
		QtCore.QObject.connect(list_widget, QtCore.SIGNAL("itemPressed(QListWidgetItem*)"),self.list_widget_clicked)
	
	def __go_back_a_directory(self):
		self.starting_directory = self.starting_directory[0:self.starting_directory.rfind('/')]
		
		for j in range(0,len(self.list_widgets)):
			self.list_widgets[j].clear()
			self.list_widget_data[j] = None
		self.__load_directory_data(self.starting_directory,self.first_list_widget)
		print "a",self.starting_directory
		
	def list_widget_clicked(self,item):
		
		if item.text() == "../": 
			self.__go_back_a_directory()
			return
		
		val = self.starting_directory+"/"
		for i,list_widget in enumerate(self.list_widgets):
			if item.listWidget() == list_widget:
				act = False
				if i == (len(self.list_widgets)-1): 
					self.__add_list_widget()
					act = True
				else:
					for j in range(i+1,len(self.list_widgets)): 
						self.list_widgets[j].clear()
						self.list_widget_data[j] = None
						
				old_item = self.list_widget_data[i]
				val += str(item.text())
				if self.__load_directory_data(val,self.list_widgets[i+1]):
					if old_item != None:
						old_item.setBackgroundColor(QtGui.QColor(255,255,255))
					item.setBackgroundColor(QtGui.QColor(64,190,0,63))	
					self.list_widget_data[i] = item
				
				if act:
					widget = self.list_widgets.pop(0)
					self.list_vbl.removeWidget(widget)
					#self.list_vbl.update()
					#widget.destroy()
					#widget.deleteLater()
					old = self.list_widget_data.pop(0)
					self.list_widgets[0].insertItem(0,"../")
					self.starting_directory += "/"+str(old.text())
					print "b",self.starting_directory
				
				return
			val += str(self.list_widget_data[i].text()) + "/"
	def list_widget_dclicked(self,item):
		pass

			
	
	def __load_directory_data(self,directory,list_widget):
		
		for root, dirs, files in os.walk(directory):
	
			list_widget.clear()
			dirs.sort()
			files.sort()
			
			if (list_widget == self.list_widgets[0]):
				"it's the first widget"
				QtGui.QListWidgetItem("../",list_widget)
				 
			for i in dirs:
				if i[0] == '.': continue
				file_length = 0
				for r, d, f in os.walk(directory+"/"+i):
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
		


if __name__ == '__main__':
	em_app = EMStandAloneApplication()
	dialog = EMSelectorDialog(None)
	em_qt_widget = EMQtWidgetModule(dialog,em_app)
	em_app.show()
	em_app.execute()

