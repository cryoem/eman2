#!/usr/bin/env python
#
# Author: David Woolford 11/10/08 (woolford@bcm.edu)
# Copyright (c) 2000-2008 Baylor College of Medicine
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


from PyQt4 import QtGui,QtCore
import os
from emapplication import EMQtWidgetModule,ModuleEventsManager
from emsprworkflow import *
from emselector import EMBrowserModule
from e2boxer import EMBoxerModule



class EMTaskMonitorModule(object):
	def __new__(cls,application):
		widget = TaskMonitorWidget(application)
		widget.resize(200,200)
		module = EMQtWidgetModule(widget,application)
		application.show_specific(module)
		return module

class TaskMonitorWidget(QtGui.QWidget):
	def get_desktop_hint(self):
		return "task_related"
	
	def __init__(self,application):
		QtGui.QWidget.__init__(self,None)
		
		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		
		self.vbl2 = QtGui.QVBoxLayout()
		self.vbl2.setMargin(0)
		self.vbl2.setSpacing(6)
		self.vbl2.setObjectName("vbl2")
		
		self.num_rows = 2
		self.num_cols = 2
		
		self.list_widget = QtGui.QListWidget(None)
		QtCore.QObject.connect(self.list_widget,QtCore.SIGNAL("itemClicked(QListWidgetItem *)"),self.list_widget_item_clicked)
			
		self.vbl2.addWidget(self.list_widget) 
		
		self.groupbox = QtGui.QGroupBox("Running tasks")
		self.groupbox.setToolTip("This is a list of currently running tasks.")
		self.groupbox.setLayout(self.vbl2)
		self.vbl.addWidget(self.groupbox)
	
	def set_entries(self,entries_dict):
		self.entries_dict=entries_dict
		self.list_widget.clear()
		for val in entries_dict:
			a = QtGui.QListWidgetItem(val,self.list_widget)
			
	def list_widget_item_clicked(self,item):
		print "item clicked",item
	

class EMWorkFlowSelector(object):
	def __new__(cls,parent,application):
		widget = EMWorkFlowSelectorWidget(parent,application)
		widget.resize(200,200)
		#gl_view = EMQtGLView(EMDesktop.main_widget,widget)
		module = EMQtWidgetModule(widget,application)
		application.show_specific(module)
		#desktop_task_widget = EM2DGLWindow(gl_view)
		return module
	
class EMWorkFlowSelectorWidget(QtGui.QWidget):
	def get_desktop_hint(self):
		return "inspector"
	
	def __init__(self,target,application):
		QtGui.QWidget.__init__(self,None)
		self.target=target
		self.application = application
		
		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		
		self.hbl_buttons2 = QtGui.QHBoxLayout()
		
		self.tree_widget = QtGui.QTreeWidget(self)
		self.tree_widget_entries = []
		self.launchers = {} # convenience launchers, to implement a pseudo factory
		self.tasks = {} # a current list of tasks that are running. Helps to avoid running the same task more than once
		self.event_managers = {} # event managers, coordinating signals and slots, same as keys in self.tasks
		self.gui_modules = [] # keeps track of all gui modules that are running
		self.module_events = [] # keeks track of gui module events managers
		
		spr = QtGui.QTreeWidgetItem(QtCore.QStringList("SPR"))
		self.launchers["SPR"] = self.launch_spr
		
		self.tree_widget_entries.append(spr)
		self.tree_widget_entries.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Browse")))
		self.launchers["Browse"] = self.launch_browser
		self.tree_widget_entries.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Boxer")))
		self.launchers["Boxer"] = self.launch_boxer_general
		self.tree_widget.insertTopLevelItems(0,self.tree_widget_entries)

		
		spr_list = []
		
		rd = QtGui.QTreeWidgetItem(QtCore.QStringList("Micrograph/CCD"))
		self.launchers["Micrograph/CCD"] = self.launch_mic_ccd
		spr_list.append(rd)
		
		ap = QtGui.QTreeWidgetItem(QtCore.QStringList("Particles"))
		self.launchers["Particles"] = self.launch_particle_report
		spr_list.append(ap)
		ctf = QtGui.QTreeWidgetItem(QtCore.QStringList("CTF"))
		self.launchers["CTF"] = self.launch_ctf_report
		spr_list.append(ctf)
		#spr_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Initial model")))
		#spr_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Refinement")))
		spr.addChildren(spr_list)
		
		
		ctf_list = []
		ctf_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Automatic fitting - e2ctf")))
		ctf_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Interactive tuning - e2ctf")))
		ctf_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Generate output - e2ctf")))
		ctf.addChildren(ctf_list)
		#self.launchers["e2ctf"] = self.launch_e2ctf_management
		
		
		ap_list = []
		ap_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Interactive boxing - e2boxer")))
		self.launchers["Interactive boxing - e2boxer"] = self.launch_e2boxer_management
		ap_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Auto boxing - e2boxer")))
		ap_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Particle import")))
		self.launchers["Particle import"] = self.launch_particle_import
		ap.addChildren(ap_list)
		
		rd_list = []
		rd_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Import micrograph/ccd")))
		self.launchers["Import micrograph/ccd"] = self.launch_import_mic_ccd
		rd.addChildren(rd_list)
		
		self.tree_widget.setHeaderLabel("Choose a task")
		
		self.hbl_buttons2.addWidget(self.tree_widget)
		
		self.close = QtGui.QPushButton("Close")
		
		self.vbl.addLayout(self.hbl_buttons2)
		self.vbl.addWidget(self.close)
		
		#QtCore.QObject.connect(self.tree_widget, QtCore.SIGNAL("itemDoubleClicked(QTreeWidgetItem*,int)"), self.tree_widget_double_click)
		QtCore.QObject.connect(self.tree_widget, QtCore.SIGNAL("itemClicked(QTreeWidgetItem*,int)"), self.tree_widget_click)
		QtCore.QObject.connect(self.close, QtCore.SIGNAL("clicked()"), self.target.close)
	
	def launch_browser(self):
		module = EMBrowserModule(self.application)
		self.application.show_specific(module)
		self.add_module(module)
		
	def launch_boxer_general(self):
		module = EMBoxerModule(self.application,None)
#		self.application.show_specific(module) # it doesn't work this way for boxer ... hmmmmm ...
		self.add_module(module)
		
	def add_module(self,module):
		self.gui_modules.append(module)
		self.module_events.append(ModuleEventsManager(self,module))
		self.update_task_list()
	
	def module_idle(self,module):
		# yes this is just the same as module_closed... thinking in progress
		for i,mod in enumerate(self.gui_modules):
			if mod == module:
				self.gui_modules.pop(i)
				self.module_events.pop(i)
				self.update_task_list()
				return
			
		print "failed to handle idle?" # this shouldn't happen if I have managed everything correctly
	
	def module_closed(self,module):
		for i,mod in enumerate(self.gui_modules):
			if mod == module:
				self.gui_modules.pop(i)
				self.module_events.pop(i)
				self.update_task_list()
				return
			
		print "failed to close module?" # this shouldn't happen if I have managed everything correctly
	
	def launch_e2ctf_management(self):
		self.launch_task(E2CTFTask,"e2ctf")
	
	def launch_ctf_report(self):
		self.launch_task(CTFReportTask,"ctf")
	
	def launch_particle_report(self):
		self.launch_task(ParticleReportTask,"particles")
	
	def launch_particle_import(self):
		self.launch_task(ParticleImportTask,"particle_import")
		
	def launch_e2boxer_management(self):
		self.launch_task(E2BoxerTask,"e2boxer")
	
	def launch_spr(self):
		self.launch_task(SPRInitTask,"spr")
			
	def launch_mic_ccd(self):
		self.launch_task(MicrographCCDTask,"mic_ccd")
		
	def launch_import_mic_ccd(self):
		self.launch_task(MicrographCCDImportTask,"import_mic_ccd")
	
	def launch_task(self,task_type,task_unique_identifier):
		'''
		You can only have one of the task forms showing at any time
		'''
		if len(self.tasks) > 0: 
			self.__clear_tasks()
		
		if not self.tasks.has_key(task_unique_identifier):
			task = task_type(self.application)
			
			task.run_form()
			self.tasks[task_unique_identifier] = task
			self.event_managers[task_unique_identifier] = TaskEventsManager(task,self,task_unique_identifier)
		else:
			self.application.show_specific(self.tasks[task_unique_identifier])
			
		self.update_task_list()
	
	def update_task_list(self):
		tasks = []
		tasks_dict = {}
		
		if len(self.tasks) > 1:
			print "error, I haven't been written to handle more than one task"
			# solution is to turn the tasks into a module and add more modules
			return
		
		for val in self.tasks.keys():
			tasks.append(str(val))
			#tasks_dict["workflow"] = 
		
		for val in self.gui_modules:
			tasks.append(str(val))
		
		self.emit(QtCore.SIGNAL("tasks_updated"),tasks)
	
	def __clear_tasks(self):
		for v in self.tasks.values():
			v.closeEvent(None)
			#self.application.close_specific(v)
		self.tasks = {}
		self.event_managers = {}
	
	def pop_task_event_pair(self,task):
		self.tasks.pop(task)
		self.event_managers.pop(task)
		self.update_task_list()
		
	def gui_running(self,task_key):
		'''
		Tells this to take the task out of the dictionary of current tasks and put in another dictionary.
		Essentially saying, "make this task disappear, but keep a reference to it until the associated gui ends"
		'''
		module = self.tasks.pop(task_key)
		self.gui_modules.append(module)
		self.module_events.append(None) # this makes the handling of the gui_modules and module_events more general, see code
		self.update_task_list()
	def gui_exit(self,module):
		for i,mod in enumerate(self.gui_modules):
			if mod == module:
				self.gui_modules.pop(i)
				self.module_events.pop(i)
				self.update_task_list()
				return
			
		print "gui exit failed" #this shouldn't happen if I have managed everything correctly
	
	def tree_widget_click(self,tree_item,i):
		task = str(tree_item.text(0))
		if task in self.launchers.keys():
			self.launchers[task]()
		else:
			if task == "Browse":
				print "hi"
	
#	def tree_widget_double_click(self,tree_item,i):
#		task = tree_item.text(0)
#		if task == "Browse":
#			self.target.add_browser_frame()
#		if task == "Thumb":
#			self.target.add_selector_frame()
#		elif task == "Boxer":
#			self.target.add_boxer_frame()


class TaskEventsManager:
	'''
	This object coordinates the signals of the WorkFlowTask objects with the EMSelectorWorkFlowWidget.
	When the WorkFlowTask is idle this object tells the SelectorWFW to release its reference, freeign memory etc
	Occasionally the WorkFlowTask will launch a gui, in which case the task is idle but it must remain in memory in order for
	signals of the spawned gui to connect properly. In which case the SelectorWFW can be told that it needs to keep a reference
	to the task.
	'''
	def __init__(self,task,selector,key):
		self.task = task
		self.selector = selector
		self.key = key
		QtCore.QObject.connect(self.task, QtCore.SIGNAL("task_idle"), self.on_task_idle) # this typically gets emitted when the user hits ok or cancel on the 
		QtCore.QObject.connect(self.task, QtCore.SIGNAL("gui_running"),self.on_gui_running) # this one 
		QtCore.QObject.connect(self.task, QtCore.SIGNAL("gui_exit"),self.on_gui_exit)
#		
	def on_gui_running(self):
		''' 
		
		'''
		self.selector.gui_running(self.key)
		
	def on_gui_exit(self):
		self.selector.gui_exit(self.task)
		
	def on_task_idle(self):
		self.selector.pop_task_event_pair(self.key)
		

class EMWorkFlowManager:
	def __init__(self,application):
		self.application = application
		self.selector = EMWorkFlowSelector(self,application)
	
		self.task_monitor = EMTaskMonitorModule(em_app)
		QtCore.QObject.connect(self.selector.qt_widget, QtCore.SIGNAL("tasks_updated"),self.task_monitor.qt_widget.set_entries)
	#window = sprinit.run_form() 
	
	def close(self):
		self.application.close_specific(self.selector)
		self.application.close_specific(self.task_monitor)

if __name__ == '__main__':
	
	from emapplication import EMStandAloneApplication
	em_app = EMStandAloneApplication()
	sprinit = EMWorkFlowManager(em_app)
	
	
	
	#em_app.show()
	em_app.execute()	