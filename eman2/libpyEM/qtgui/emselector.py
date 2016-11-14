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

from EMAN2 import get_image_directory, get_dtag, EMData, \
	get_files_and_directories, db_open_dict, strip_file_tag, remove_file, \
	remove_directories_from_name, Util, EMUtil, IMAGE_UNKNOWN, base_name, \
	file_exists, base_name
from EMAN2db import EMAN2DB, db_convert_path, db_open_dict, db_check_dict, e2getcwd
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from emapplication import ModuleEventsManager, EMApp, get_application
from emimage2d import EMImage2DWidget
from emimagemx import EMImageMXWidget
from emimage3diso import EMIsosurfaceModel
from emimage3dslice import EM3DSliceModel
from emimage3dsym import EM3DSymModel
from emimage3dvol import EMVolumeModel
from emimageutil import EMTransformPanel
from emplot2d import EMPlot2DWidget
#from e2simmxxplor import EMSimmxExplorer
from emsave import save_data
import PyQt4
import math
import os
import re
import weakref

read_header_only = True
EMAN2DB = "EMAN2DB"

MDS = "%" # metadata separator

class EMActionDelegate:
	'''
	interface for action delegates - they are notified when the widget owning them is closed
	'''
	def closeEvent(self,event): pass
	
class EMItemAction:
	'''
	interface for single item actions
	'''
	
	def item_action(self,item,target): raise NotImplementedException
	
class EMMultiItemAction:
	'''
	interface for multiple item actions
	'''
	
	def multi_item_action(self,items,target): raise NotImplementedException
	
class EMSaveItemAction(EMItemAction,EMMultiItemAction,EMActionDelegate):
	
	def item_action(self,item,target):
		data = item.get_data()
		if data != None: save_data(data)
		
	def multi_item_action(self,items,target):
		for item in items:
			data = item.get_data()
			if data != None: 
				result = save_data(data)
				if not result: return
		
class EMSaveStackSaveAction(EMMultiItemAction,EMActionDelegate):
	
	def multi_item_action(self,items,target):
		save_data(items) # save_data knows how to deal with items

class EMDeleteItemAction(EMItemAction,EMMultiItemAction,EMActionDelegate):
	
	def item_action(self,item,target):
		self.__delete_items( [item] )

	def __delete_items(self,items):
		msg = QtGui.QMessageBox()
		msg.setText("Deletion will be permanent. Are you sure you want to delete the selected file(s)?")
		s = ""
		for i in items: s+=i.text()+"\n"
		msg.setInformativeText(s)
		msg.setStandardButtons(QtGui.QMessageBox.Cancel | QtGui.QMessageBox.Ok )
		msg.setDefaultButton(QtGui.QMessageBox.Cancel)
		ret = msg.exec_()
		if ret == QtGui.QMessageBox.Cancel: return False
		
		for item in items:
			delegate= item.get_delegate()
			delegate.delete_url(item.get_url())
			
	def multi_item_action(self,items,target):
		self.__delete_items(items)
			

			
def DataDisplayModuleTemplate(Type,get_data_attr="get_data",data_functors=[],usescenegraph=False):
	'''
	This would be similar to policy based design using templated inheritance
	Type is a type of EMAN2 display module. Search for DataDisplayModuleTemplate
	'''
	class DataDisplayModule(EMItemAction,EMActionDelegate):
		def __init__(self):
			self.module_type = Type
			self.display_modules = [] # a list of display modules
			self.module_events = [] # a list of ModuleEventsManagers
			self.get_data_attr = get_data_attr
			self.data_functors = data_functors # functors that can be called once the data is acquired
			
		def item_action(self,item,target):
			from EMAN2 import Transform
			from emdataitem3d import EMDataItem3D, EMIsosurface
			from emshapeitem3d import EMCube
			
			name = os.path.basename(str(item.get_url()))
			single_mode = target.single_preview_only()
			if single_mode and len(self.display_modules) != 0:
				old_module = self.display_modules[-1]
				data = getattr(item,self.get_data_attr)()
				for funct in self.data_functors: 
					funct(data)
				if self.module_type == EMPlot2DWidget: # aka template specialization
					old_module.set_data(data,item.get_url())
				elif self.module_type == EM3DSymModel:
					old_module.model.set_data(data)
				else: 
					if usescenegraph:
						emdata3d = EMDataItem3D(data, transform=Transform())
						#emdata3d.setSelectedItem(True)
						isosurface = EMIsosurface(emdata3d, transform=Transform())
						old_module.insertNewNode(name, emdata3d, parentnode=old_module)
						old_module.insertNewNode("Iso", isosurface, parentnode=emdata3d)
						old_module.initialViewportDims(emdata3d.getData().get_xsize())
						old_module.setCurrentSelection(isosurface)	# Set isosurface to display upon inspector loading
					else:
						old_module.set_data(data)
				old_module.setWindowTitle(item.get_url())
				old_module.show()
				old_module.updateGL()
				return
			
			from e2simmxxplor import EMSimmxExplorer

			if self.module_type == EM3DSymModel: #TODO: get correct symmetry or switch to e2eulerxplor.py
				from emimage3dsym import EMSymViewerWidget
				widget = EMSymViewerWidget()
			elif self.module_type in (EMIsosurfaceModel, EMVolumeModel, EM3DSliceModel, EMSimmxExplorer):
				from emglobjects import EM3DGLWidget
				widget = EM3DGLWidget()
				model = self.module_type(widget)
				widget.set_model(model)
			else: 
				widget = self.module_type()
				
			data = getattr(item,self.get_data_attr)()
			for funct in self.data_functors: 
				funct(data)
			
			if self.module_type == EM3DSymModel:
				widget.model.set_data(data)
			elif self.module_type == EMPlot2DWidget:
				widget.set_data(data,item.get_url())
			else:
				if usescenegraph:
					emdata3d = EMDataItem3D(data, transform=Transform())
					#emdata3d.setSelectedItem(True)
					isosurface = EMIsosurface(emdata3d, transform=Transform())
					widget.insertNewNode(name, emdata3d, parentnode=widget)
					widget.insertNewNode("Iso", isosurface, parentnode=emdata3d)
					widget.initialViewportDims(emdata3d.getData().get_xsize())
					widget.setCurrentSelection(isosurface)	# Set isosurface to display upon inspector loading
				else:
					widget.set_data(data)
			self.display_modules.append(widget)
			self.module_events.append(ModuleEventsManager(self,widget))
			widget.show()
			widget.updateGL()
			widget.setWindowTitle(item.get_url())
			
		def module_closed(self,module):
			for i,mod in enumerate(self.display_modules):
				if mod == module:
					p = self.display_modules.pop(i)
					mod = self.module_events.pop(i)
					mod.disconnect_object()
					return
				
			print "failed to close module?" # this shouldn't happen if I have managed everything correctly
		
		def closeEvent(self,event):
			self.module_events = [] # this should take care of complications that arise if self.module_closed is called
			for module in self.display_modules:
				module.closeEvent(event)
				
			self.display_modules = []
			
	return DataDisplayModule

class EM2DStackPreviewAction(DataDisplayModuleTemplate(EMImageMXWidget,"get_2d_stack"),EMMultiItemAction):
	'''
	This is like a template specialization of the DataDisplayModuleTemplate in the case of
	using an EMImageMXWidget. The reason is because we support a special "Preview Subset" action.
	'''
	def multi_item_action(self,items,target):
		single_mode = target.single_preview_only()
		data = []
		from emimagemx import ApplyAttribute
		for item in items:
			data.append([item.image_path(),item.get_idx(),[ApplyAttribute("Img #",item.get_idx())]])
		
		from emimagemx import EMLightWeightParticleCache
		data = EMLightWeightParticleCache(data)
		
		if single_mode and len(self.display_modules) != 0:
			old_module = self.display_modules[-1]
			old_module.set_data(data)
			old_module.setWindowTitle(item.get_url())
			old_module.updateGL()
			return
		
		module = self.module_type()
		module.set_data(data)

		self.display_modules.append(module)
		self.module_events.append(ModuleEventsManager(self,module))
		module.show()
		module.updateGL()
		module.setWindowTitle(item.get_url())

# Global location for common and unique strings
DELETE = "Delete"
SAVE_AS ="Save As"
VIEWER_3D = "3D Viewer"
VOLUME_VIEWER =  "Volume Viewer"
SLICE_VIEWER = "Slice Viewer"
SINGLE_2D_VIEWER = "Single 2D"
MULTI_2D_VIEWER = "Multi 2D"
PLOT_2D_VIEWER = "Plot 2D"
PLOT_3D_VIEWER = "Plot 3D"
EULER_VIEWER = "Euler View"
SIMMX_EULER_VIEWER = "Simmx Euler View"
PREVIEW_SUBSET = "Preview Subset"
SAVE_SUBSET = "Save Subset"


def EMSelectorBaseTemplate(Type):
	'''
	This is templated inheritance. I need the selector to be a dialog, and I need it to be a normal widget.
	See the EMSelectorDialogType and EMBrowserType
	Types currently in use are the QtGui.QWidget and the QtGui.QDialog
	'''
	class EMSelectorBase(Type):
		def __init__(self, single_selection=False):
			'''
			@param single_selection - should selections be limited to singles?
			'''
			Type.__init__(self,None)
			self.setFocusPolicy(Qt.StrongFocus)
#			self.module=weakref.ref(module) # Avoid strong cycle
			self.single_selection = single_selection # Flag indicating single selection in interface
			self.browse_delegates = [EMBDBDelegate(self), EMFileSystemDelegate(self)] # Object capable of returning listed items based on url- Add your own
			
			self.hbl = QtGui.QVBoxLayout(self)
			self.hbl.setMargin(0)
			self.hbl.setSpacing(6)
			self.hbl.setObjectName("hbl")
			
			self.setWindowIcon(QtGui.QIcon(get_image_directory() + "/display_icon.png"))
			
			self.__init_filter_combo()
			
			self.current_force = None # use to keep track of what is currently being forced, either 2D plotting, 2D image showing, or neither
			self.selections = []
			self.current_list_widget = None
			self.lock = True
			self.list_widgets = []
			self.previews = [] # keeps track of all of the preview windows
#			self.module_events = [] # used to coordinate signals from the modules, especially close events, to free memory
			self.list_widget_data= [] # entries should be tuples containing (current folder item)
			self.splitter = QtGui.QSplitter(self)
			self.splitter.setChildrenCollapsible(False)

			self.add_list_widget()
			self.add_list_widget()
			
			self.hbl.addWidget(self.splitter,1)
			
			self.__load_url(e2getcwd(),self.list_widgets[0])
	
			self.bottom_hbl = QtGui.QHBoxLayout()
			self.bottom_hbl.addWidget(self.filter_text,0)
			self.bottom_hbl.addWidget(self.filter_combo,1)
			self.__init_buttons()
			self.bottom_hbl.addWidget(self.cancel_button,0)
			self.bottom_hbl.addWidget(self.ok_button,0)
			self.hbl.addLayout(self.bottom_hbl)
			
			self.lock = False
			
			self.paint_events = 0
			
			self.timer_interval = 500 # half a second
			self.timer = QtCore.QTimer()
			QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.time_out) # for auto refresh
			
			self.timer.start(self.timer_interval)
			
			self.selected_files = []
			
			get_application().attach_child(self)
			
		def __init_buttons(self):
			self.ok_button = QtGui.QPushButton("Ok")
			self.ok_button.adjustSize()
			
			self.cancel_button = QtGui.QPushButton("Cancel")
			self.cancel_button.adjustSize()
		
			QtCore.QObject.connect(self.ok_button, QtCore.SIGNAL("clicked(bool)"),self.ok_button_clicked)
			QtCore.QObject.connect(self.cancel_button, QtCore.SIGNAL("clicked(bool)"),self.cancel_button_clicked)
		
		def ok_button_clicked(self,bool):
			''' Slot for OK button '''
			#print "EMSelectorBase.ok_button_clicked"
			self.emit(QtCore.SIGNAL("ok"),self.selections)
			self.emit(QtCore.SIGNAL("oky"))
		
		def cancel_button_clicked(self,bool):
			''' Slot for Cancel button '''
			#print "EMSelectorBase.cancel_button_clicked"
			self.emit(QtCore.SIGNAL("cancel"),self.selections)
		
		
		def __del__(self):
			pass
		
		def set_selection_text(self,text):
			'''
			Selection label is a QLabel, by default its text is "Save As", but you can change it
			Called in emsprworkflow
			'''
			self.selection_label.setText(text)
		
		def time_out(self):
			'''
			This function takes care of automatic updates - if a url is detected to have changed
			then so does the information we display
			This is somewhat inefficient but it works
			'''
			if self.lock: return
			self.lock = True
			for idx,widget in enumerate(self.list_widgets):
				if widget.get_mod_time() != None and widget.get_delegate() != None:
					mod_time = widget.get_delegate().url_mod_time(widget.get_url())
					if mod_time != widget.get_mod_time():
						old_items = [widget.item(i) for i in range(widget.count())]
						new_items = widget.get_delegate().get_items(widget.get_url())
						old_set = set([str(item.text()) for item in old_items if item.get_name() != EMUpArrowItem.NAME ]) # forget about the up arrow, it just confuses things
						new_set = set([str(item.text()) for item in new_items])
						
						added = new_set - old_set
						removed = old_set - new_set
						
						if len(added) > 0:
							for k in added:
								new_item = (item for item in new_items if item.text() == k).next()
								widget.addItem(new_item)

						if len(removed) > 0:
							rm = []
							clear_flag=False
							for k in removed:
								old_item = (item for item in old_items if item.text() == k).next()
								if old_item.isSelected():clear_flag = True
								rm.append(widget.row(old_item))
							
							rm.sort()
							rm.reverse()
							for rm_idx in rm: widget.takeItem(rm_idx)
							
							if clear_flag:
								for j in range(idx+1,len(self.list_widgets)):
									self.list_widgets[j].clear()
						widget.set_mod_time(mod_time)
						self.lock = False
						return
							
			self.lock = False
					
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
			'''
			Clears the widgets and refreshes widget 0
			This is not the best way of doing it
			'''
			self.lock = True
			directory = folderize(self.list_widgets[0].get_url())
			for list_widget in self.list_widgets: list_widget.clear()
			self.__load_url(directory,self.list_widgets[0])
			self.lock = False
	
		def add_list_widget(self):
			list_widget = EMListWidget(self)
			
			
			#list_widget.contextMenuEvent = self.list_widget_context_menu_event
			
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

		def list_widget_dclicked(self,item):
			'''
			Inheriting class should supply specialize this
			'''
			pass
		
		def list_widget_context_menu_event(self,event):
			'''
			Inheriting class should supply specialize this
			'''
			pass
		
		def __go_up_a_directory(self,item=None):
			
			if item != None:
				#url = item.get_url()
				url = self.list_widgets[0].get_url()
				for delegate in self.browse_delegates:
					if delegate.handles_url(url): 
						new_dir = delegate.parent_url(url)
						break
				else:
					raise RuntimeError("Unknown url %s" %url)
			else:
				raise
			
#			if len(new_dir) == 0: new_dir = get_dtag()
#			elif  new_dir[-1] == ":": new_dir += get_dtag() # C: becomes C:/
#			elif len(new_dir) == 1 and new_dir != get_dtag(): new_dir += ":/"# C becomes C:/

			self.lock = True	
			for j in range(0,len(self.list_widgets)):
				widget = self.list_widgets[j]
				widget.clear()
				self.list_widget_data[j] = None
					
			self.__load_url(new_dir,self.list_widgets[0])
			self.lock = False
			
		def __go_forward_a_directory(self):
			self.lock = True
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
				#directory = folderize(directory) + str(self.list_widget_data[i].text())
				#directory += dtag + str(self.list_widget_data[i].text())
			
			# The last widget must have these attributes removed if they existed,
			# or else the time_out function will produce incorrect results
			last_widget = self.list_widgets[-1]

			a = EMUpArrowItem(None,"../",self.list_widgets[0].get_url())
			self.list_widgets[0].insertItem(0,a)
			
			self.lock = False
			
		def update_selections(self):
			'''
			Makes the list of currently selected files accurate and up to date. Called when
			something has been clicked in a a list widget
			'''
			self.selections = []
			items = self.current_list_widget.selectedItems()
			if len(items) > 0:
				a = items[0]
				self.selections = self.valid_emdata_urls(items)
					
		def valid_emdata_urls(self,items):
			'''
			Get the valid EMData urls from the list of items. Useful in EMAN2 when selecting images, etc
			'''
			# This ALWAYS required save files to be valid EMData objects, which was not always desirable...
			#return [item.emdata_save_as_url() for item in items if item.emdata_save_as_url() != None]
		
			ret=[]
			for item in items:
				if item.emdata_save_as_url() != None : ret.append(item.emdata_save_as_url())
				elif item.get_url() != None: ret.append(item.get_url())

			return ret
		
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

		def list_item_selected(self,item):
			if self.lock : return
			#if self.current_list_widget == None: return
			if item.get_name() == EMGenericItem.NAME: return #it's just a value in the db
			
			self.update_selections()
			if item == None: return		
			
			if item.get_name() == EMUpArrowItem.NAME: 
				self.__go_up_a_directory(item)
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
						
		def __load_url(self,url,list_widget):
			'''
			A way to load items in the list_widget using a url
			especially at start up
			'''
			get_application().setOverrideCursor(Qt.BusyCursor)
			
			list_widget.clear()
			if (list_widget == self.list_widgets[0]):
				self.lock = True
				a = EMUpArrowItem(None,"../",url)
				list_widget.addItem(a)
				self.lock = False

#			print self.browse_delegates
#			self.browse_delegates.sort(key=fspsort)

			for delegate in self.browse_delegates:
				if delegate.handles_url(url): 
					items = delegate.get_items(url)
					for item in items:
						if item.get_url()!=None and ".hed" in item.get_url(): continue
						list_widget.addItem(item)
					ret = True
					list_widget.set_url(url)
					list_widget.set_mod_time(delegate.url_mod_time(url))
					list_widget.set_delegate(delegate)
					break
			else:
				raise RuntimeError("unknown url %s" %url)
				ret = False
				
			get_application().setOverrideCursor(Qt.ArrowCursor)
			return ret
		
		def __load_url_from_item(self,list_widget,item):
			'''
			A way to load items in the list_widget using data from the item
			i.e. user clicked on the item
			'''
			return self.__load_url(item.get_url(),list_widget)
			
	
	return EMSelectorBase

def fspsort(x):
#	print x
	if x[:4].lower()=="bdb:" : return x
	if "." not in x : return x
	y=x.rsplit(".",1)
	return y[1]+"."+y[0]

EMBrowserType = EMSelectorBaseTemplate(QtGui.QWidget)
class EMBrowser(EMBrowserType):
	def __init__(self, single_selection=False, usescenegraph=False):
		EMBrowserType.__init__(self,single_selection)
		
		self.add_list_widget() # 3 panels in browsing mode
		self.usescenegraph = usescenegraph
		
		self.__init_action_delegates()
		
		bottom_hbl2 = QtGui.QHBoxLayout()
		self.__init_preview_options()
		bottom_hbl2.addWidget(self.preview_options,0)
		self.hbl.addLayout(bottom_hbl2)
				
		bottom_hbl3 = QtGui.QHBoxLayout()
		self.__init_plot_options()
		bottom_hbl3.addWidget(self.replace,0)
		bottom_hbl3.addWidget(self.include,0)
		
		self.groupbox = QtGui.QGroupBox("Plot/3D options")
		self.groupbox.setLayout(bottom_hbl3)
		self.groupbox.setEnabled(False)
		
		bottom_hbl2.addWidget(self.groupbox)
		self.resize(480,480)
	
		#Below: from old EMBrowserModule
		self.setWindowTitle("EMAN2 Browser")
		self.preview_options.setCurrentIndex(0)
		self.preview_options_changed(self.preview_options.currentText())
		self.ok_button.setEnabled(False)
		self.cancel_button.setEnabled(False)
		#End of old EMBrowserModule code
	
	def __del__(self):
		pass
	
	def closeEvent(self,event):
		for delegate in self.action_delegates.values():
			delegate.closeEvent(event)
		EMBrowserType.closeEvent(self, event)
	
	def __init_action_delegates(self):
		'''
		All of the available actions for the context menu (right click)
		'''
		self.action_delegates = {}
		if self.usescenegraph:
			from emscene3d import EMScene3D
			self.action_delegates[VIEWER_3D] = DataDisplayModuleTemplate(EMScene3D, usescenegraph=self.usescenegraph)()
		else:
			from emimage3d import EMImage3DWidget
			self.action_delegates[VIEWER_3D] = DataDisplayModuleTemplate(EMImage3DWidget)()
		from emimagemx import ApplyProcessor
		self.action_delegates[VOLUME_VIEWER] = DataDisplayModuleTemplate(EMVolumeModel,data_functors=[ApplyProcessor("normalize",{})])()
		self.action_delegates[SLICE_VIEWER] = DataDisplayModuleTemplate(EM3DSliceModel,data_functors=[ApplyProcessor("normalize",{})])()
		self.action_delegates[SINGLE_2D_VIEWER] = DataDisplayModuleTemplate(EMImage2DWidget)()
		stack_action = EM2DStackPreviewAction()
		self.action_delegates[MULTI_2D_VIEWER] = stack_action #DataDisplayModuleTemplate(EMImageMXWidget,"get_2d_stack")
		self.action_delegates[PLOT_2D_VIEWER] = DataDisplayModuleTemplate(EMPlot2DWidget)()
		from emplot3d import EMPlot3DWidget
		self.action_delegates[PLOT_3D_VIEWER] = DataDisplayModuleTemplate(EMPlot3DWidget)()
		self.action_delegates[EULER_VIEWER] = DataDisplayModuleTemplate(EM3DSymModel)()
		from e2simmxxplor import EMSimmxExplorer
		self.action_delegates[SIMMX_EULER_VIEWER] = DataDisplayModuleTemplate(EMSimmxExplorer,"get_url")()
		
		self.action_delegates[SAVE_AS] = EMSaveItemAction()
		self.action_delegates[DELETE] = EMDeleteItemAction()
		self.action_delegates[PREVIEW_SUBSET] = stack_action # can use the same object
		self.action_delegates[SAVE_SUBSET] = EMSaveStackSaveAction()
		
	def __init_plot_options(self):
		self.replace = QtGui.QRadioButton("Replace")
		self.include = QtGui.QRadioButton("Include")
		self.include.setChecked(True)

	def __init_preview_options(self):
		self.preview_options = QtGui.QComboBox(self)
		#self.preview_options.addItem("No preview")
		self.preview_options.addItem("Single preview")
		self.preview_options.addItem("Multi preview")
		#self.preview_options.setCurrentIndex(0)
		
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
	
	def list_widget_dclicked(self,item):
		self.preview_item(item)

	def preview_item(self,item):
		'''
		previews the item (loads an appropriate display module) if possible
		Returns True if the item was loaded into a display module, otherwise returns False
		'''
		preview_occured = False
		get_application().setOverrideCursor(Qt.BusyCursor)
		view_action = item.default_view_action()
		if view_action != None:
			self.action_delegates[view_action].item_action(item,self)
			preview_occured = True
		get_application().setOverrideCursor(Qt.ArrowCursor)
			
		return preview_occured
	
	def list_widget_context_menu_event(self,event):
		'''
		Dynamic menu creation depending on item type
		'''
		event.accept()
		focus = self.current_list_widget
		l = focus
		if focus == None: return
		selected_items = l.selectedItems()
		if len(selected_items) == 0: return
		
		menu = QtGui.QMenu()
		self.menu_selected_items = selected_items
		if len(selected_items) == 1:
			first_item = selected_items[0]
			actions = first_item.actions()
			if actions == None: return
			for action in actions: menu.addAction(action)
		else:
			multi_actions = []
			common = [SAVE_AS,DELETE]
			for item in selected_items:
				actions = item.actions()
				if actions == None: 
					common = [] # the items has no actions - abort
					break 
				rm = []
				for i,c in enumerate(common):
					if c not in actions:
						rm.append(i)
				rm.reverse()
				for idx in rm: common.pop(idx)
						
				if len(common) == 0: break
				
			for action in common: menu.addAction(action)
			
			all_2d = True
			
			# some specialized behavior for 2D images with the same dimensions
			for item in selected_items:
				if not item.get_name() == EM2DImageItem.NAME:
					all_2d = False
					break

			if all_2d:
				all_same_dims = True
				nx = -1
				ny = -1
				for item in selected_items:
					md = item.get_metadata()
					if nx == -1 and ny == -1:
						nx = md["nx"]
						ny = md["ny"]
					else:
						if nx != md["nx"] or ny != md["ny"]:
							all_same_dims = False
							
				if all_same_dims:
					menu.addAction(PREVIEW_SUBSET)
					menu.addAction(SAVE_SUBSET)
				

		QtCore.QObject.connect(menu,QtCore.SIGNAL("triggered(QAction*)"),self.menu_action_triggered)
		self.action_list_widget = l # only set if the menu acutally triggers
		menu.exec_(event.globalPos())
		
	def menu_action_triggered(self,action):
		'''
		Slot for right click menu
		'''
		items = self.menu_selected_items
		
		total = len(items)
		if total == 0: return
		if total == 1 and self.action_delegates.has_key(str(action.text())):
			get_application().setOverrideCursor(Qt.BusyCursor)
			self.action_delegates[str(action.text())].item_action(items[0],self)
			get_application().setOverrideCursor(Qt.ArrowCursor)
			return
		else:
			get_application().setOverrideCursor(Qt.BusyCursor)
			self.action_delegates[str(action.text())].multi_item_action(items,self)
			get_application().setOverrideCursor(Qt.ArrowCursor)
			return

EMSelectorDialogType = EMSelectorBaseTemplate(QtGui.QDialog)
class EMSelectorDialog(EMSelectorDialogType):
	def __init__(self,single_selection=False,save_as_mode=True): #TODO: figure out whether save_as_mode is needed (unused)
		EMSelectorDialogType.__init__(self,single_selection)	

		hbl2=QtGui.QHBoxLayout()
		hbl2.setMargin(0)
		hbl2.setSpacing(2)
		self.selection_label = QtGui.QLabel(SAVE_AS,self)
		hbl2.addWidget(self.selection_label)
		self.save_as_line_edit = QtGui.QLineEdit("",self)
		hbl2.addWidget(self.save_as_line_edit,0)
		self.hbl.insertLayout(1,hbl2)
		self.dialog_mode = True
		self.validator = None # an optional feature which makes the behavior of the dialog more sophisticated - see emsave.py
		self.dialog_result = ""
		
		self.ok_button.setDefault(True)
			
		self.resize(480,400)
		self.setWindowTitle("EMAN2 Selector")
		
	def exec_(self):
		'''
		Wraps QtGui.QDialog.exec_
		@return a list of selected filenames
		'''
		QtGui.QDialog.exec_(self)
		return self.dialog_result
	
	def set_validator(self,validator):
		'''
		Sets the validator which is used only in save_as mode (when this is being used as a file selecting dialog)
		An example validator is an emsave.EMSaveImageValidator
		'''
		self.validator = validator
		
	def get_current_directory(self):
		directory = self.list_widgets[0].get_url()
		for i in range(len(self.list_widgets)-1):
			items = self.list_widgets[i].selectedItems()
			if len(items) >  0:
				directory = folderize(directory) + items[0].text()
				#directory += "/" + items[0].text()
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
		directory = folderize(directory)
		#if directory[-1] != "/": directory += "/"
		
		return directory
	
	def update_selections(self):
		'''
		Makes the list of currently selected files accurate and up to date. Called when
		something has been clicked in a a list widget
		'''
		# now make the list of selections reflect what is need to display them using
		EMSelectorDialogType.update_selections(self)
		
		text = ""
		for i,name in enumerate(self.selections):
			if i > 0:
				text+= "; "
			text += base_name(name)
			
		self.save_as_line_edit.setText(text)
	
	def cancel_button_clicked(self,bool):
		EMSelectorDialogType.cancel_button_clicked(self, bool)
		
		self.accept()
	
	def ok_button_clicked(self,bool):
		EMSelectorDialogType.ok_button_clicked(self, bool)
		
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
	
	
class EMListWidget(QtGui.QListWidget):
	'''
	Customized ListWidget as displayed in the browser
	'''
	def __init__(self,target,*args):
		self.target = weakref.ref(target)
		QtGui.QListWidget.__init__(self,*args)
		self.reset_vars()
	
	def clear(self):
		self.reset_vars()
		QtGui.QListWidget.clear(self)
		
	def contextMenuEvent(self,event):
		self.target().list_widget_context_menu_event(event)
		
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
	
	def get_data(self,full_path,idx=0):
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
	
	def parent_url(self,url):
		'''
		Get the parent url - for example, if the up arrow is hit
		'''
		raise NotImplementedException("Inheriting classes must supply this functionality")
	
	def delete_url(self,url):
		'''
		delete the url
		'''
		raise NotImplementedException("Inheriting classes must supply this functionality")
	
	def emdata_save_as_url(self,url):
		'''
		What would be the name of the EMData if this url where being used to call EMData.write_image
		May return None to indicate that the function call makes no sense
		'''
		return None

def folderize(url):
	'''
	makes sure the last char in the url is "/" 
	'''
	if url[-1] == "/": return url
	else: return url + "/"

class EMFileSystemDelegate(EMBrowseDelegate):
	'''
	The Delegate for reading file system contents 
	Main function is to return the the correct list items to the EMSelector 
	'''
	
	def __init__(self,target):
		self.target = weakref.ref(target)
		self.threed_dim_limit = 128
		
	def get_data(self,full_path,idx=0):
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
		This function is called by EM2DStackItem and EM3DImageItem
		Return is a cache that can be treated like a list of EMData objects
		Can give speed ups
		'''
		from emimagemx import EMLightWeightParticleCache,EM3DDataListCache
		md = self.get_metadata(full_path,0)
		if md["nz"] > 1: return EM3DDataListCache(full_path)
		else: return EMLightWeightParticleCache.from_file(full_path)
	
	def url_mod_time(self,url):
		'''
		Get the last time the url was modified
		May return None to indicate the call is not valid/supported - this will mean the corresponding list widget will not be automatically updated
		'''
		if os.path.exists(url): return os.stat(url)[-2]
		else: return None
	
	def emdata_save_as_url(self,url):
		return url
	
	def handles_url(self,url):
		'''
		Expected interface Delegate::handles_url
		'''
		if "EMAN2DB" in url: return False
		
		vals = url.split(MDS)
		if os.path.exists(vals[0]): return True
		
		return False
		
	def delete_url(self,url):
		remove_file(url)
		
	def parent_url(self,url):
		vals = url.split(MDS)
		if len(vals) == 1:
			if os.path.isfile(url):
				return os.path.dirname(url)
			elif os.path.isdir(url):
				if len(url) == 1: return url # probably at the root
				vals = url.split("/")
				l = 0
				idx = -1
				if len(vals[-1]) == 0: # '/' at the end
					idx = -2
					l = 1
				
				l += len(vals[idx])
				return url[:-l]
			else:
				raise RuntimeError("Unknown url %s" %url)
			
		else:
			# metadata
			l = 0
			idx = -1
			if len(vals[-1]) == 0: 
				l = 1
				idx = -2
				
			if len(vals) > math.fabs(idx):
				l += len(vals[idx]) + 1 # plus one for the split
			
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
		if os.path.isdir(url): return self.__get_dir_items(url)
		elif os.path.isfile(url): return self.__get_file_metadata_items(url)
		else: return self.__get_metadata_items(url)
		
	def __get_dir_items(self,url):
		'''
		This function called when the url is a directory
		'''
		if not os.path.isdir(url): raise RuntimeError("%s is not a directory" %url)
		
		dirs, files = get_files_and_directories(url)
		if len(files) == 0 and len(dirs) == 0:
			# something is unusual about the directory, for instance, it is a file
			return []
		
		dirs.sort()
		files.sort(key=fspsort)
		 
		return_items = []
		
		for i in dirs:
			if i[0] == '.': continue

			d = remove_directories_from_name(i)
			
			if d == "EMAN2DB": return_items.append(EMBDBFolderItem(self,"bdb",folderize(url)+"EMAN2DB")) # the EMBDBDelegate will know how handle this url
			else: return_items.append(EMFSFolderItem(self,i,folderize(url)+i))
			
		for file in files:
			a = self.__get_file_item(file,url)
			if a != None: return_items.append(a)
	
		return return_items
	
	
	
	def __get_file_metadata_items(self,url):
		'''
		This function called when the url is a file
		'''
		if not os.path.isfile(url): raise RuntimeError("%s is not a file" %url)
		
		return_items = []
		item = None
		extension = Util.get_filename_ext(url)
		if EMUtil.get_image_ext_type(extension) != IMAGE_UNKNOWN and extension not in ["png","jpeg","jpg","JPG"]:
			e = EMData()
			e.read_image(url,0, True)
			if EMUtil.get_image_count(url) > 1:
				
			
				if e.get_zsize() > 1:
					return_items = [EM3DMetaImageItem(self,str(i),url,i) for i in xrange(0,EMUtil.get_image_count(url))]
				else:
					return_items = [EM2DMetaImageItem(self,str(i),url,i) for i in xrange(0,EMUtil.get_image_count(url))]
			else:
				d = e.get_attr_dict()
				keys = d.keys()
				keys.sort() #alphabetical order
				return_items = [EMDataHeaderItem(self,str(k)+" : "+str(d[k]),url,k,d[k]) for k in keys]
			

		return return_items
	
	def __get_metadata_items(self,url):
		vals = url.split(MDS)
		if not os.path.isfile(vals[0]): raise RuntimeError("% is not a valid file system url" %url)
		
		return_items = []
		item = None
		extension = Util.get_filename_ext(vals[0])
		if EMUtil.get_image_ext_type(extension) != IMAGE_UNKNOWN and extension not in ["png","jpeg","jpg","JPG"]:
			e = EMData()
			e.read_image(vals[0],0, True)
			val_idx = 1
			if EMUtil.get_image_count(vals[0]) > 1:
				e.read_image(vals[0],int(vals[1]), True)
				val_idx += 1
			d = e.get_attr_dict()
				
			if len(vals) == val_idx:
				keys = d.keys()
				keys.sort() #alphabetical order
				return_items = [EMDataHeaderItem(self,str(k)+" : "+str(d[k]),url,k,d[k]) for k in keys]
			elif len(vals) == val_idx+1:
				val = d[vals[val_idx]]
				return_items = [EMGenericItem(self,str(val))]
			else: pass

		return return_items
		
		
	def __get_file_item(self,file,directory):
		'''
		Called internally
		@returns the correct item for the given file in the given directory, returns None in exceptional circumstances
		'''
		filt = self.target().get_file_filter()
		if file[0] == '.': return None
		if file[-1] == '~': return None
		
		extension = Util.get_filename_ext(file)
		full_name = folderize(directory)+file
		
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
				item = EMGenericFileItem(self,file,full_name)
		
		elif EMPlot2DWidget.is_file_readable(full_name):
			item = EMFSPlotItem(self,file,full_name)
		else:
			if filt != "EM types":
				item = EMGenericFileItem(self,file,full_name)

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

	def get_attr_dict(self):
		'''
		A wrapper for get_metadata - originally added to allow
		a generic interface in EMStackSaveDialog
		Unfortunately needs to stay for the time being
		'''
		return self.get_metadata()
	
	def emdata_save_as_url(self):
		'''
		What would be the name of the EMData if this url where being used to call EMData.write_image
		May return None to indicate that the function call makes no sense
		'''
		return None
	
	def get_url(self):
		'''	Returns the file name path, for BDB stuff this should in full bdb syntax'''
		return None
	
	def get_data(self):
		'''
		Generally returns an EMData objec. May return a list of floats, for example.
		Is customized to return a list of EMDatas, or a Cache, in some instances
		'''
		return None
	
	def get_metadata(self):
		'''
		Returns a dictionary if defined
		'''
		return None
	def actions(self):
		'''
		Return a list of strings defining feasible view actions
		'''
		return None
	
	def default_view_action(self):
		'''
		Return the default view action
		'''
		return None
	
class EMUpArrowItem(EMListItem):
	ICON = None
	NAME = "go up a directory"
	
	def __init__(self,delegate,text,url):
		EMListItem.__init__(self,delegate,text)
		self.url = url # this should be the url that can be handed to EMBrowseDelegate.parent_url

	def get_icon(self):
		'''
		Supply your own Icon
		'''
		if EMUpArrowItem.ICON == None:
			EMUpArrowItem.ICON = QtGui.QIcon(get_image_directory() + "/up_arrow.png")
		return EMUpArrowItem.ICON
	
	def get_name(self): return EMUpArrowItem.NAME
	
	
	def get_url(self): return self.url
	
class EMDataListItem(EMListItem):
	'''
	Objects of this type are items that list EMData types -
	they happen to have the fact that they can be deleted and saved in common, that is all
	ABSTRACT
	'''
	def __init__(self,delegate=None,text="",full_name=""):
		EMListItem.__init__(self,delegate,text)
		self.full_path = full_name
		self.context_menu_options[DELETE] = self.delete
		self.context_menu_options[SAVE_AS] = self.save_as
		
	def delete(self,target):
		try:
			remove_file(self.full_path)
			return True
		except: return False
		
	def save_as(self,target):
		save_data(self.get_data())
	
	def get_url(self):
		return self.full_path

	
	def actions(self):
		'''
		Return a list of strings defining feasible view actions
		'''
		return [DELETE,SAVE_AS]

class EMDataHeaderItem(EMListItem):
	NAME = "EMData Header"
	def __init__(self,delegate=None,text="",url="",key="",value=None):
		EMListItem.__init__(self,delegate,text)
		self.url = url
		self.key = key
		self.value = value
		
	def get_name(self): return EMDataHeaderItem.NAME

	def get_url(self):
		return self.url +MDS+str(self.key)

class EMStack2DCapableMixin:
	'''
	a 2D stack capable item is something that knows how to supply
	data to the EMImageMX set_data function.
	It's a mixin
	'''
	def get_2d_stack(self):
		raise NotImplementedException("Inheriting classes must supply this function")

class EM2DStackItem(EMDataListItem,EMStack2DCapableMixin):
	ICON = None
	NAME = "2D stack"
	def get_name(self): return EM2DStackItem.NAME
	
	def get_icon(self):
		'''Supply your own Icon	'''
		if EM2DStackItem.ICON == None:
			EM2DStackItem.ICON = QtGui.QIcon(get_image_directory() + "/multiple_images.png")
		return EM2DStackItem.ICON
	
	def emdata_save_as_url(self): return self.delegate.emdata_save_as_url(self.get_url())
	
	def get_data(self):
		'''This one returns a list'''
		return self.delegate.get_stack_data(self.full_path)
	
	def get_2d_stack(self):
		return self.delegate.get_stack_data(self.full_path)
	
	def get_metadata(self):
		if self.metadata == None:
			self.metadata = self.delegate.get_metadata(self.full_path,0)
		return self.metadata
	
	def default_view_action(self): return MULTI_2D_VIEWER
	
	def actions(self):
		from e2simmx import PROJ_FILE_ATTR,PART_FILE_ATTR
		ret = [DELETE,SAVE_AS]
		ret.append(SINGLE_2D_VIEWER)
		ret.append(MULTI_2D_VIEWER)
		
		md = self.get_metadata()
		if md.has_key("xform.projection"):
			# the assumption is that they all have the xform.projection header attribute, which could be fatal
			ret.append(EULER_VIEWER)
		if md.has_key(PROJ_FILE_ATTR) and md.has_key(PART_FILE_ATTR): 
			ret.append(SIMMX_EULER_VIEWER)
		return  ret
	
class EM3DStackItem(EMDataListItem):
	ICON = None
	NAME = "3D stack"
	def get_name(self): return EM3DStackItem.NAME
	
	def emdata_save_as_url(self): return self.delegate.emdata_save_as_url(self.get_url())
	
	def get_icon(self):
		'''
		Supply your own Icon
		'''
		if EM3DStackItem.ICON == None:
			EM3DStackItem.ICON = QtGui.QIcon(get_image_directory() + "/multiple_images_3d.png")
		return EM3DStackItem.ICON
	# no preview for this item as of Feb 2009
	def get_data(self):
		return self.delegate.get_stack_data(self.full_path)
	
class EM2DImageItem(EMDataListItem):
	ICON = None
	NAME = "2D image"
	def __init__(self,delegate=None,text="",full_name="",idx=0):
		'''
		Have to supply init because we need the idx
		'''
		EMDataListItem.__init__(self,delegate,text,full_name)
		self.idx = idx
		
	def emdata_save_as_url(self): return self.delegate.emdata_save_as_url(self.get_url())
	
	def image_path(self): return self.full_path
	def get_idx(self): return self.idx
	
	def get_name(self): return EM2DImageItem.NAME
	
	def get_icon(self):
		'''
		Supply your own Icon
		'''
		if EM2DImageItem.ICON == None:
			EM2DImageItem.ICON = QtGui.QIcon(get_image_directory() + "/single_image.png")
		return EM2DImageItem.ICON

	def get_metadata(self):
		if self.metadata == None:
			self.metadata = self.delegate.get_metadata(self.full_path,self.idx)
		return self.metadata
	
	def get_data(self):
		return self.delegate.get_data(self.full_path,self.idx)
	
	def default_view_action(self): return SINGLE_2D_VIEWER
	
	def actions(self):
		ret = [DELETE,SAVE_AS]
		ret.append(SINGLE_2D_VIEWER)
		
		md = self.get_metadata()
		ret.append(PLOT_2D_VIEWER)
		if md["ny"] < 1025 and md["nx"] < 1025:
			ret.append(VIEWER_3D)
			ret.append(PLOT_3D_VIEWER)
			
		return ret
	
class EM3DImageItem(EM2DImageItem,EMStack2DCapableMixin):
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
	
	def default_view_action(self):
		md = self.get_metadata()
		if md["nx"] > 256 or md["ny"] > 256 or md["nz"] > 256:
			return MULTI_2D_VIEWER
		else: return VIEWER_3D
	
	def actions(self):
		ret = [DELETE,SAVE_AS]
		ret.append(VIEWER_3D)
		ret.append(VOLUME_VIEWER)
		ret.append(SLICE_VIEWER)
		ret.append(SINGLE_2D_VIEWER)
		ret.append(MULTI_2D_VIEWER)
		return ret
	
	def get_2d_stack(self):
		return self.delegate.get_stack_data(self.full_path)

class EM2DMetaImageItem(EM2DImageItem):
	'''
	This is a 2D Image in a stack
	'''
#	NAME = "2D meta image"
#	def get_name(self): return EM2DMetaImageItem.NAME
	
	def get_url(self):
		return self.full_path+MDS+str(self.idx)
	
	def actions(self):
		'''
		No DELETE
		'''
		ret = [SAVE_AS]
		ret.append(SINGLE_2D_VIEWER)
		
		md = self.get_metadata()
		ret.append(PLOT_2D_VIEWER)
		if md["ny"] < 1025 and md["nx"] < 1025:
			ret.append(VIEWER_3D)
			ret.append(PLOT_3D_VIEWER)
			
		return ret

class EM3DMetaImageItem(EM3DImageItem):
	'''
	This is a 3D Image in a stack
	'''
#	NAME = "3D meta image"
#	def get_name(self): return EM3DMetaImageItem.NAME
	
	def get_url(self):
		return self.full_path+MDS+str(self.idx)
	
	def actions(self):
		'''
		No Delete
		'''
		ret = [SAVE_AS]
		ret.append(VIEWER_3D)
		ret.append(VOLUME_VIEWER)
		ret.append(SLICE_VIEWER)
		ret.append(SINGLE_2D_VIEWER)
		ret.append(MULTI_2D_VIEWER)
		return ret
	
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
	
	def get_data(self):
		return EMPlot2DWidget.get_data_from_file(self.get_url())
		
	def get_url(self): return self.full_path
	
	def default_view_action(self): return PLOT_2D_VIEWER
	
	def actions(self):
		ret = [DELETE]
		ret.append(PLOT_2D_VIEWER)
		return  ret
	
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
	
	def actions(self):
		ret = [DELETE]
		return  ret

class EMBDBDelegate(EMBrowseDelegate):
	def __init__(self,target):
		self.target = weakref.ref(target)
		self.directory_replacements = {"EMAN2DB":"bdb"}
	
	def get_data(self,full_path,idx=0):
		'''
		All items delegate by this should call this function to get a fully loaded EMData
		That way the read routine is in the one location
		'''
		e = EMData()
		e.read_image(db_convert_path(full_path),idx)
		return e
		
	def get_metadata(self,full_path,idx=0):
		'''
		All items that load metadata using EMData io call this function
		'''
		db_name =  db_convert_path(full_path)
		db = db_open_dict(db_name,ro=True)
		data = db.get_header(idx)
		return data
	
	def get_stack_data(self,full_path):
		'''
		This function is called by EM2DStackItem and EM3DImageItem
		Return is a cache that can be treated like a list of EMData objects
		'''
		from emimagemx import EMLightWeightParticleCache,EM3DDataListCache
		md = self.get_metadata(full_path,0)
		if md["nz"] > 1: return EM3DDataListCache(db_convert_path(full_path))
		else: return EMLightWeightParticleCache.from_file(db_convert_path(full_path))
	
	def url_mod_time(self,url):
		'''
		Get the last time the url was modified
		May return None to indicate the call is not valid/supported - this will mean the corresponding list widget will not be automatically updated
		'''
		if os.path.exists(url): return os.stat(url)[-2]
		else: return None
	
	def emdata_save_as_url(self,url):
		return db_convert_path(url)
	
	def delete_url(self,url):
		remove_file(db_convert_path(url))
	
	def handles_url(self,url):
		
		return "EMAN2DB" in url		
	
	def parent_url(self,url):
		
		vals = url.split(MDS)
		if len(vals) == 1:
			#if not os.path.isdir(url): raise RuntimeError("Unknown url %s" %url)

			vals = url.split("/")
			if len(url) == 1: return url # probably at the root
			idx = -1
			l = 0
			if len(vals[-1]) == 0: # '/' at the end
				idx = -2
				l = 1
			
			l += len(vals[idx]) + 1 # plus one for the '/'
			return url[:-l]
		else:
			# metadata
			l = 0
			idx = -1
			if len(vals[-1]) == 0: 
				l = 1
				idx = -2
				
			if len(vals) > math.fabs(idx):
				l += len(vals[idx]) + 1 # plus one for the split
				
			return url[:-l]
	
	def get_items(self,url):
#		if self.__is_database_file(url): 
#			list_widget.clear()
#			return self.__get_database_data_items(url)
		if self.handles_url(url): #os.path.exists(url) and (url.endswith("EMAN2DB") or url.endswith("EMAN2DB/")):
			if os.path.isdir(url): return self.__get_bdb_directory_items(url)
			elif self.__is_image_url(url):return self.__get_bdb_image_items(url)
			else:return self.__get_bdb_type_items(url)
		else: raise RuntimeError("Unknown url %s" %url)
	
	def __is_image_url(self,url):
		vals = url.split(MDS)
		is_image = False
		try:
			db_url = db_convert_path(vals[0])
			e = EMData()
			e.read_image(db_url,0,False)
			is_image = True
		except: pass
		return is_image
	
	def __get_bdb_image_items(self,url):
		if not self.__is_image_url(url): raise RuntimeError("Unknown image url %s" %url)
		vals = url.split(MDS)
		db_url = db_convert_path(vals[0])
		
		return_items = []
		
		db = db_open_dict(db_url,ro=True)
		n = len(db) 
		if len(vals) == 1:
			for i in range(n):
				d = db.get_header(i)
				if d!=None and d.has_key("nz") : break
			if  n > 1:
				if d["nz"] > 1: return_items = [EM3DMetaImageItem(self,str(i),url,i) for i in xrange(0,n)]
				else: return_items = [EM2DMetaImageItem(self,str(i),url,i) for i in xrange(0,n)]
			else: 
				keys = d.keys()
				keys.sort() #alphabetical order
				return_items = [EMDataHeaderItem(self,str(k)+" : "+str(d[k]),url,k,d[k]) for k in keys]
		else:
			val_idx = 1
			if len(vals) > 1 and n > 1:
				val_idx += 1
				d = db.get_header(int(vals[1]))
			else:
				d = db.get_header(0)
				
			if len(vals) == val_idx:
				keys = d.keys()
				keys.sort() #alphabetical order
				return_items = [EMDataHeaderItem(self,str(k)+" : "+str(d[k]),url,k,d[k]) for k in keys]
			elif len(vals) == val_idx+1:
				val = d[vals[val_idx]]
				return_items = [EMGenericItem(self,str(val))]
			else:
				raise RuntimeError("Unknown url %s" %url) # an EMData header only goes so deep
			
		return return_items
	
	def __get_bdb_type_items(self,url):
		
		vals = url.split(MDS)
		
		return_items = []
		if len(vals) == 1:
			return []
		else:
			#if vals[0][-1] != "/": vals[0] += "/" # this is so the db_convert_path function works
			vals[0] = folderize(vals[0])
			db_url = db_convert_path(vals[0])
			db = db_open_dict(db_url+vals[1])
			
			for db_key in vals[2:]:db = db[db_key]
			
			try:
				for k,val in db.items():
					if isinstance(val,dict):
						return_items.append(EMBDBDictItem(self,str(k),url,str(k)))
					else:
						return_items.append(EMBDBKeyValueItem(self,str(k),url,str(k)))
			except:
				return_items = [EMGenericItem(self,str(db))]
			#print db.keys()
		
		return return_items
		
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

			a = EMBDBDirectoryItem(self,i,folderize(url)+i)
			return_items.append(a)
			
		for file in files:
			a = self.__get_bdb_file_item(file,real_directory)
			if a != None: return_items.append(a)
		return return_items

	def __get_bdb_file_item(self,file,real_directory):
		if not file[-3:] == "bdb": return None
		f = file.rpartition(".bdb")
		db_directory = self.__get_database_directory(real_directory)
		
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
#					from emapplication import EMErrorMessageDisplay
#					EMErrorMessageDisplay.run(["Warning: the %s database might be corrupted." %db_name], "Data loss" )
					
				return None
		
		if db and len(db) > 0:
			#n = DB[f[0]]["maxrec"]
			n = len(db)
			if n > 1:
				for i in range(n):
					d = db.get_header(i)
					try:
						if d!=None and d.has_key("nz") : break
					except: 
						return EMBDBItem(self,f[0],real_directory+MDS+f[0])
				if d["nz"] == 1:
					#a = EM2DStackItem(self,f[0],"bdb:"+db_directory+"#"+f[0])
					a = EM2DStackItem(self,f[0],folderize(real_directory)+f[0])
				elif d["nz"] > 1:
					#a = EM3DStackItem(self,f[0],"bdb:"+db_directory+"#"+f[0])
					a = EM3DStackItem(self,f[0],folderize(real_directory)+f[0])
			elif n == 1:
				d = db.get_header(0)
				if d["nz"] <= 1:
					#a = EM2DImageItem(self,f[0], "bdb:"+db_directory+"#"+f[0],0)
					a = EM2DImageItem(self,f[0],folderize(real_directory)+f[0],0)
				else:
					#a = EM3DImageItem(self,f[0], "bdb:"+db_directory+"#"+f[0],0)
					a = EM3DImageItem(self,f[0],folderize(real_directory)+f[0],0)
					
			else:
				a = EMBDBItem(self,f[0],real_directory+MDS+f[0])
		else:
			a = EMBDBItem(self,f[0],real_directory+MDS+f[0])
			
		a.file_name = file
		return a
		
	def __get_database_directory(self,file):
		'''
		Get the database where EMAN2DB should be opening in order to open the given file
		e.g. if db path is /home/someone/work/EMAN2DB/data.bdb will return /home/someone/work
		'''
		idx1 = file.find("EMAN2DB")
		if idx1 > 0:
			return file[0:idx1-1]
		else: return None
	
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

class EMBDBItem(EMListItem):
	ICON = None
	NAME = "bdb dict"
	def __init__(self,delegate,text,db_url):
		EMListItem.__init__(self,delegate,text)
		self.full_path = db_url
		
	def get_name(self): return EMBDBItem.NAME
	
	def get_icon(self):
		'''
		Supply your own Icon
		'''
		if EMBDBItem.ICON == None:
			EMBDBItem.ICON = QtGui.QIcon(get_image_directory() + "/database.png")

		return EMBDBItem.ICON

	def get_url(self): return self.full_path

class EMBDBKeyValueItem(EMListItem):
	'''
	Most basic BDB item, basic string display
	'''
	NAME = "BDB key value"
	def __init__(self,delegate=None,text="",url="",key=""):
		EMListItem.__init__(self,delegate,text)
		self.url = url
		self.key = key
		
	def get_name(self): return EMBDBKeyValueItem.NAME

	def get_url(self):
		return self.url +MDS+str(self.key)

class EMBDBDictItem(EMListItem):
	'''
	Most basic BDB item, basic string display
	'''
	NAME = "BDB dict"
	ICON = None
	def __init__(self,delegate=None,text="",url="",key=""):
		EMListItem.__init__(self,delegate,text)
		self.url = url
		self.key = key
		
	def get_name(self): return EMBDBDictItem.NAME

	def get_url(self):
		return self.url +MDS+str(self.key)
	
	def get_icon(self):
		'''
		Supply your own Icon
		'''
		if EMBDBDictItem.ICON == None:
			EMBDBDictItem.ICON = QtGui.QIcon(get_image_directory() + "/Bag.png")

		return EMBDBDictItem.ICON

class EMGenericItem(EMListItem):
	'''
	A dead end item, displays a value. Has no metadata and there is no consequence for clicking on it. 
	'''
	NAME = "generic"
	def __init__(self,delegate=None,text="",key=None):
		EMListItem.__init__(self,delegate,text)
		self.key = key
		
	def get_name(self): return EMGenericItem.NAME

class EMGenericFileItem(EMGenericItem):
	'''
	A generic item that can be deleted
	'''
	def actions(self):
		ret = [DELETE]
		return  ret
	
	def get_url(self):
		return str(self.key)


app = None


if __name__ == '__main__':
	em_app = EMApp()
#	em_qt_widget = EMSelectorModule(save_as_mode=False)
#    em_app.show()
	
	em_selector = EMSelectorDialog()
	files = em_selector.exec_()
	#print files
	#print "Press Ctrl-C to exit" #FIXME: figure out why Ctrl-C is required to terminate the program
	em_app.exit(0)
	em_app.execute() 
