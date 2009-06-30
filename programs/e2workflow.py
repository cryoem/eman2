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
from PyQt4.QtCore import Qt
import os
from emapplication import EMQtWidgetModule,ModuleEventsManager
from emsprworkflow import *
from emtprworkflow import *
from emselector import EMBrowserModule
from e2boxer import EMBoxerModule
from EMAN2 import process_running,kill_process
from emanimationutil import Animator
from emimage import EMModuleFromFile

import time
import weakref

class EMTaskMonitorModule(EMQtWidgetModule):
	def __init__(self,application=None):
		self.application = application
		self.widget = EMTaskMonitorWidget(self,application)
		self.widget.resize(256,200)
		EMQtWidgetModule.__init__(self,self.widget)

class EMTaskMonitorWidget(QtGui.QWidget,Animator):
	def get_desktop_hint(self):
		return "workflow"

	def __init__(self,module,application):
		QtGui.QWidget.__init__(self,None)
		self.setWindowIcon(QtGui.QIcon(get_image_directory() + "workflow.png"))
		self.setWindowTitle("Running Tasks")
		Animator.__init__(self)
		self.module = weakref.ref(module)
		# animation specific
		self.timer_interval = 500 # half a second time interval
		self.register_animatable(self)
		# remember the number of entries in the history db when this was initialized
		# so that keeping track of processes using pids is more efficient
		global HOMEDB
		HOMEDB=EMAN2db.EMAN2DB.open_db()
		HOMEDB.open_dict("history")
		try: self.init_history_db_entries = HOMEDB.history["count"]
		except: self.init_history_db_entries = 0
		if self.init_history_db_entries == None: self.init_history_db_entries = 0
		#print self.init_history_db_entries
		if db_check_dict("bdb:project"):
			project_db = db_open_dict("bdb:project")
			self.current_process_info = project_db.get("workflow.process_ids",{})

			# If a project is copied from one machine to another while a job is running,
			# the job table will be corrupt. If we run into this, we clear the table
			for k in self.current_process_info.keys(): 
				try: x=HOMEDB.history[k]["pid"]
				except: 
					project_db["workflow.process_ids"]={}
					self.current_process_info = {}
					break
		else:
			self.current_process_info = {}
		self.entries_dict = {}
		self.history_check = []
		
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
		self.kill = QtGui.QPushButton("Kill")
		
		self.vbl.addWidget(self.kill)
		
		self.set_entries(self.entries_dict)
	
		QtCore.QObject.connect(self.kill, QtCore.SIGNAL("clicked()"), self.on_kill)
		
	def animate(self,time):
		global HOMEDB
		HOMEDB=EMAN2db.EMAN2DB.open_db()
		HOMEDB.open_dict("history")
		
		if not HOMEDB.history: return
		
		for i,pid in enumerate(self.history_check):
			try: total = HOMEDB.history["count"]+1
			except: return True
			for j in range(self.init_history_db_entries+1,HOMEDB.history["count"]+1):
				try:
					if HOMEDB.history[j]["pid"] == pid:
						project_db = db_open_dict("bdb:project")
						#print "adding",j
						self.current_process_info[j] = HOMEDB.history[j]
						self.history_check.pop(i) # ok what if two are ended at the same time? oh well there is time lag...
						project_db ["workflow.process_ids"] = self.current_process_info
						db_close_dict("bdb:project")
						self.set_entries(self.entries_dict)
						return True
				except:
					continue
		
		
		for key in self.current_process_info.keys():
			#if not process_running(HOMEDB.history[key]["pid"]):
			if not process_running(HOMEDB.history[key]["pid"]):
				project_db = db_open_dict("bdb:project")
				self.current_process_info.pop(key)
				self.set_entries(self.entries_dict)
				project_db ["workflow.process_ids"] = self.current_process_info
				db_close_dict("bdb:project")
				self.set_entries(self.entries_dict)
				return True
			else:
				self.current_process_info[key] = HOMEDB.history[key]
			#else:
				#print "process running",HOMEDB.history[key]["pid"]
		
		# kill defunct child processes
		# i.e. If I don't do this defunct processes will still show up as running
		if len(self.current_process_info) > 0: 
			try:
				pid,stat = os.waitpid(0,os.WNOHANG)
			except: pass 
			
		self.set_entries(self.entries_dict)
		return True
				
			
	def set_entries(self,entries_dict):
		selected_items = self.list_widget.selectedItems() # need to preserve the selection
		
		s_text = None
		if len(selected_items) == 1 :
			s = str(selected_items[0].text())
			if len(s) >= 4:
				s_text= s[:4]
		elif len(selected_items) > 1:
			print "can't handle more than one selection in the task manager" # yes this should work, no time
		
		
		self.entries_dict=entries_dict
		self.list_widget.clear()
		for k,val in entries_dict.items():
			a = QtGui.QListWidgetItem(val,self.list_widget)
			if s_text != None and len(val) >=4:
				if val[:4] == s_text:
					a.setSelected(True)
			a.module = k
		self.list_processes(s_text)
		
	def list_processes(self,s_text):
		
		for key in self.current_process_info.keys():
			d = self.current_process_info[key]
			pid = str(d["pid"])
			prog = get_file_tag(d["args"][0])
			progress = "0%"
			try:
				progress = str(float(d["progress"])*100)
				if len(progress) > 4:
					progress = progress[:4]
				progress += "%"
			except: pass
			t = str(time.ctime(d["start"]))
			s = pid + "\t" +progress+"\t"+ prog + "\t" + t
			
			a =  QtGui.QListWidgetItem(s,self.list_widget)
			
			if s_text != None and len(s) >=4:
				if s[:4] == s_text:
					a.setSelected(True)
					
			a.module = "process"
			a.pid = pid
	
	def accrue_process_info(self):
		
		for i in range(len(self.current_process_info),len(self.current_processes)):
			for j in range(self.init_history_db_entries,HOMEDB.history["count"]):
				print db[j]["pid"],self.current_processes[i]
				
	
	def add_process(self,pid):
		self.history_check.append(pid) # the program hasn't added the entry to the db (There is lag)
	
	def list_widget_item_clicked(self,item):
		#val = None
		#for k,v in  self.entries_dict.items():
			#if v == str(item.text()):
				#val = k
				#break
			
		#if val == None:
			#print "error, this shouldn't happen"
		#print item.module
		self.module().emit(QtCore.SIGNAL("task_selected"),str(item.text()),item.module)
	
	def on_kill(self):
		selected_items = self.list_widget.selectedItems()
		if len(selected_items) == 0:
			return
		
		# not sure if there will ever be more than one selected but oh well let's support closing them all
		
		for item in selected_items:
			if item.module == "process":
				if not kill_process(int(item.pid)):
					print "kill process failed for some reason"
			else:
				#print self.entries_dict
				self.entries_dict.pop(item.module)
				self.module().emit(QtCore.SIGNAL("task_killed"),str(item.text()),item.module)
				self.set_entries(self.entries_dict)
			
	def keyPressEvent(self,event):
		if event.key() == Qt.Key_F1:
			try:
				import webbrowser
				webbrowser.open("http://blake.bcm.edu/emanwiki/e2workflow")
				return
			except: pass

			try: from PyQt4 import QtWebKit
			except: return
			try:
				try:
					test = self.browser
				except: 
					self.browser = QtWebKit.QWebView()
					self.browser.load(QtCore.QUrl("http://blake.bcm.edu/emanwiki/e2workflow"))
					self.browser.resize(800,800)
				
				if not self.browser.isVisible(): self.browser.show()
			except: pass		
		
		#self.list_widget.
			#print "removing item",item
			#self.list_widget.removeItemWidget(item)

class EMWorkFlowSelector(EMQtWidgetModule):
	def __init__(self,application,task_monitor):
		self.application = application
		self.widget = EMWorkFlowSelectorWidget(self,application,task_monitor)
		#self.widget.resize(256,200)
		EMQtWidgetModule.__init__(self,self.widget)


#class EMWorkFlowSelector(object):
	#def __new__(cls,parent,application,task_monitor):
		#widget = EMWorkFlowSelectorWidget(parent,application,task_monitor)
		#widget.resize(200,200)
		##gl_view = EMQtGLView(EMDesktop.main_widget,widget)
		#module = EMQtWidgetModule(widget,application)
		#application.show_specific(module)
		##desktop_task_widget = EM2DGLWindow(gl_view)
		#return module
	
class EMWorkFlowSelectorWidget(QtGui.QWidget):
	def get_desktop_hint(self):
		return "workflow"
	
	def __init__(self,module,application,task_monitor):
		'''
		Target is not really currently used
		application is required
		task_monitor should probably be supplied, it should be an instance of a EMTaskMonitorWidget
		'''
		QtGui.QWidget.__init__(self,None)
		self.__init_icons()
		self.setWindowIcon(self.icons["desktop"])
		self.setWindowTitle("EMAN2 Tasks")
		self.application = weakref.ref(application)
		self.module = weakref.ref(module)
		self.task_monitor = task_monitor
		
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
		
		spr = self.__get_spr_tree(self.tree_widget)
		self.tree_widget_entries.append(spr)

		self.tree_widget_entries.append(self.__get_tomo_part_recon_tree(self.tree_widget))
		self.tree_widget_entries.append(self.__get_utilities_tree(self.tree_widget))
		self.tree_widget_entries.append(self.__get_generic_interfaces_tree(self.tree_widget))
		
		browser_entry = QtGui.QTreeWidgetItem(QtCore.QStringList("Browse"))
		browser_entry.setIcon(0,self.icons["display_icon"])
		self.tree_widget_entries.append(browser_entry)
		self.launchers["Browse"] = self.launch_browser

		
		self.tree_widget.insertTopLevelItems(0,self.tree_widget_entries)
		self.tree_widget.setHeaderLabel("Choose a task")

		self.hbl_buttons2.addWidget(self.tree_widget)
		
		self.vbl.addLayout(self.hbl_buttons2)
		#self.vbl.addWidget(self.close)
		
		#QtCore.QObject.connect(self.tree_widget, QtCore.SIGNAL("itemDoubleClicked(QTreeWidgetItem*,int)"), self.tree_widget_double_click)
		QtCore.QObject.connect(self.tree_widget, QtCore.SIGNAL("itemClicked(QTreeWidgetItem*,int)"), self.tree_widget_click)
		#QtCore.QObject.connect(self.close, QtCore.SIGNAL("clicked()"), self.target.close)

	def __init_icons(self):
		self.icons = {}
		self.icons["single_image"] = QtGui.QIcon(get_image_directory() + "single_image.png")
		self.icons["single_image_3d"] = QtGui.QIcon(get_image_directory() + "single_image_3d.png")
		self.icons["eulerxplor"] = QtGui.QIcon(get_image_directory() + "eulerxplor.png")
		self.icons["ctf"] = QtGui.QIcon(get_image_directory() + "ctf.png")
		self.icons["green_boxes"] = QtGui.QIcon(get_image_directory() + "green_boxes.png")
		self.icons["desktop"] = QtGui.QIcon(get_image_directory() + "desktop.png")
		self.icons["eman"] = QtGui.QIcon(get_image_directory() + "eman.png")
		self.icons["multiple_images"] = QtGui.QIcon(get_image_directory() + "multiple_images.png")
		self.icons["multiple_images_3d"] = QtGui.QIcon(get_image_directory() + "multiple_images_3d.png")
		self.icons["black_box"] = QtGui.QIcon(get_image_directory() + "black_box.png")
		self.icons["refine"] = QtGui.QIcon(get_image_directory() + "refine.png")
		self.icons["plot"] = QtGui.QIcon(get_image_directory() + "plot.png")
		self.icons["display_icon"] = QtGui.QIcon(get_image_directory() + "display_icon.png")		
		self.icons["feather"] = QtGui.QIcon(get_image_directory() + "feather.png")
		self.icons["classes"] = QtGui.QIcon(get_image_directory() + "classes.png")
		self.icons["tomo_hunter"] = QtGui.QIcon(get_image_directory() + "tomo_hunter.png")
		
	def __get_spr_tree(self,tree_widget):
		
		spr = QtGui.QTreeWidgetItem(QtCore.QStringList("Single Particle Reconstruction"))
		self.launchers["Single Particle Reconstruction"] = self.launch_spr
		
		
		spr_list = []
		
		rd = QtGui.QTreeWidgetItem(QtCore.QStringList("Raw Data"))
		rd.setIcon(0,self.icons["single_image"])
		self.launchers["Raw Data"] = self.launch_mic_ccd_report#EMRawDataReportTask
		
		rd_list = []
		rd_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Filter Raw Data")))
		rd_list[-1].setIcon(0,self.icons["single_image"])
		self.launchers["Filter Raw Data"] = self.launch_filter_mic
		rd.addChildren(rd_list)
 
		spr_list.append(rd)
		
		ap = QtGui.QTreeWidgetItem(QtCore.QStringList("Particles"))
		ap.setIcon(0,self.icons["green_boxes"])
		self.launchers["Particles"] = self.launch_particle_report
		spr_list.append(ap)
		ctf = QtGui.QTreeWidgetItem(QtCore.QStringList("CTF"))
		ctf.setIcon(0,self.icons["ctf"])
		self.launchers["CTF"] = self.launch_ctf_report
		spr_list.append(ctf)
		mis = QtGui.QTreeWidgetItem(QtCore.QStringList("Particle Sets"))
		mis.setIcon(0,self.icons["multiple_images"])
		spr_list.append(mis)
		self.launchers["Particle Sets"] = self.launch_mis_report
		refine2d = QtGui.QTreeWidgetItem(QtCore.QStringList("Reference Free Class Averages"))
		self.launchers["Reference Free Class Averages"] = self.launch_refine2d_report
		refine2d.setIcon(0,self.icons["classes"])
		spr_list.append(refine2d)
		init_model = QtGui.QTreeWidgetItem(QtCore.QStringList("Initial Model"))
		init_model.setIcon(0,self.icons["single_image_3d"])
		self.launchers["Initial Model"] = self.launch_initmodel_report
		spr_list.append(init_model)
		refinement = QtGui.QTreeWidgetItem(QtCore.QStringList("3D Refinement"))
		refinement.setIcon(0,self.icons["refine"])
		self.launchers["3D Refinement"] = self.launch_refinement_report
		spr_list.append(refinement)
		spr_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Eulers")))
		spr_list[-1].setIcon(0,self.icons["eulerxplor"])
		self.launchers["Eulers"] = self.launch_eulers
		resolution = QtGui.QTreeWidgetItem(QtCore.QStringList("Resolution"))
		resolution.setIcon(0,self.icons["plot"])
		self.launchers["Resolution"] = self.launch_resolution_report
		spr_list.append(resolution)
		
		#spr_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Initial model")))
		#spr_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Refinement")))
		spr.addChildren(spr_list)
		
		
		ctf_list = []
		ctf_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Automated Fitting - e2ctf")))
		self.launchers["Automated Fitting - e2ctf"] = self.launch_e2ctf_auto_ft
		ctf_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Interactive Tuning - e2ctf")))
		self.launchers["Interactive Tuning - e2ctf"] = self.launch_e2ctf_tune
		ctf_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Generate Output - e2ctf")))
		self.launchers["Generate Output - e2ctf"] =  self.launch_e2ctf_write_ouptut
		ctf_list[0].setIcon(0,self.icons["ctf"])
		ctf_list[1].setIcon(0,self.icons["ctf"])
		ctf_list[2].setIcon(0,self.icons["ctf"])
		ctf.addChildren(ctf_list)
		#self.launchers["e2ctf"] = self.launch_e2ctf_management
		
		
		ap_list = []
		ap_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Interactive Boxing - e2boxer")))
		self.launchers["Interactive Boxing - e2boxer"] = self.launch_e2boxer_gui
		ap_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Auto Boxing - e2boxer")))
		self.launchers["Auto Boxing - e2boxer"] = self.launch_e2boxer_auto
		ap_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Generate Output - e2boxer")))
		self.launchers["Generate Output - e2boxer"] = self.launch_e2boxer_output
		ap_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Particle Import")))
		self.launchers["Particle Import"] = self.launch_particle_import
		ap_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Coordinate Import")))
		self.launchers["Coordinate Import"] = self.launch_ptcl_coord_import
		ap_list[0].setIcon(0,self.icons["green_boxes"])
		ap_list[1].setIcon(0,self.icons["green_boxes"])
		ap_list[2].setIcon(0,self.icons["green_boxes"])
		ap_list[3].setIcon(0,self.icons["green_boxes"])
		ap_list[4].setIcon(0,self.icons["green_boxes"])
		ap.addChildren(ap_list)
		
		mis_list = []
#		mis_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Examine Particles")))
#		self.launchers["Examine Particles"] = self.launch_examine_particle_stacks
		mis_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Build Particle Sets")))
		self.launchers["Build Particle Sets"] = self.launch_make_ptcl_set
		mis_list[0].setIcon(0,self.icons["multiple_images"])
#		mis_list[1].setIcon(0,self.icons["multiple_images"])
		mis.addChildren(mis_list)
		
		refine2d_list = []
		refine2d_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Generate Classes (Sets) - e2refine2d")))
		refine2d_list[-1].setIcon(0,self.icons["classes"])
		self.launchers["Generate Classes (Sets) - e2refine2d"] = self.launch_refine2d_choose_stacks
#		refine2d_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Generate Classes (Particles) - e2refine2d")))
#		self.launchers["Generate Classes (Particles) - e2refine2d"] = self.launch_refine2d_choose_ptcls
#		refine2d_list[-1].setIcon(0,self.icons["classes"])
		refine2d.addChildren(refine2d_list)
		
		init_model_list = []
		init_model_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Make Model- e2initialmodel")))
		init_model_list[-1].setIcon(0,self.icons["single_image_3d"])
		self.launchers["Make Model- e2initialmodel"] = self.launch_e2makeinitial
#		init_model_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Import model")))
#		self.launchers["Import model"] = self.launch_import_initial_model
		init_model.addChildren(init_model_list)
		
		refine_list = []
		refine_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Run (Sets) - e2refine")))
		refine_list[-1].setIcon(0,self.icons["refine"])
		self.launchers["Run (Sets) - e2refine"] = self.launch_e2refine_sets
#		refine_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Run (Particles) e2refine")))
#		refine_list[-1].setIcon(0,self.icons["refine"])
#		self.launchers["Run (Particles) e2refine"] = self.launch_e2refine
		
		refinement.addChildren(refine_list)
		
		resolution_list = []
		resolution_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Run e2eotest")))
		resolution_list[-1].setIcon(0,self.icons["plot"])
		self.launchers["Run e2eotest"] = self.launch_e2eotest
		resolution_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Run e2resolution")))
		resolution_list[-1].setIcon(0,self.icons["plot"])
		self.launchers["Run e2resolution"] = self.launch_e2resolution
		resolution.addChildren(resolution_list)
		
		tree_widget.expandItem(spr)
		tree_widget.expandItem(rd)
		tree_widget.expandItem(ap)
		tree_widget.expandItem(mis)
		tree_widget.expandItem(ctf)	
		tree_widget.expandItem(refine2d)
		tree_widget.expandItem(init_model)
		tree_widget.expandItem(refinement)
		tree_widget.expandItem(resolution)
		
		return spr
	
	def __get_utilities_tree(self,tree_widget):
		
		util = QtGui.QTreeWidgetItem(QtCore.QStringList("Utilities"))

		util_list = []
		util_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Eulers Tool")))
		util_list[-1].setIcon(0,self.icons["eulerxplor"])
		self.launchers["Eulers Tool"] = self.launch_asym_unit
		history = QtGui.QTreeWidgetItem(QtCore.QStringList("History"))
		history.setIcon(0,self.icons["feather"])
		util_list.append(history)
		self.launchers["History"] = self.launch_view_history
		self.directory = QtGui.QTreeWidgetItem(QtCore.QStringList("Working directory"))
		self.directory.setIcon(0,self.icons["desktop"])
		self.directory.setToolTip(0,os.getcwd())
		util_list.append(self.directory)
		preferences = QtGui.QTreeWidgetItem(QtCore.QStringList("Preferences"))
		preferences.setIcon(0,self.icons["desktop"])
		util_list.append(preferences)
		self.launchers["Preferences"] = self.launch_view_preferences
		self.launchers["Working directory"] = self.launch_change_directory
		lights = QtGui.QTreeWidgetItem(QtCore.QStringList("3D Lights"))
		lights.setIcon(0,self.icons["desktop"])
		util_list.append(lights)
		self.launchers["3D Lights"] = self.launch_lights_tool
		fonts = QtGui.QTreeWidgetItem(QtCore.QStringList("3D Fonts"))
		fonts.setIcon(0,self.icons["desktop"])
		util_list.append(fonts)
		self.launchers["3D Fonts"] = self.launch_fonts_tool
		plot3d = QtGui.QTreeWidgetItem(QtCore.QStringList("3D Plot"))
		plot3d.setIcon(0,self.icons["plot"])
		util_list.append(plot3d)
		self.launchers["3D Plot"] = self.launch_3d_plot_tool
		self.launchers["3D Volume Viewer"] = self.launch_3d_volume_viewer
		desktop = QtGui.QTreeWidgetItem(QtCore.QStringList("3D Desktop"))
		desktop.setIcon(0,self.icons["desktop"])
		util_list.append(desktop)
		self.launchers["3D Desktop"] = self.launch_desktop
		iso3d = QtGui.QTreeWidgetItem(QtCore.QStringList("3D Isosurface Viewer"))
		iso3d.setIcon(0,self.icons["single_image_3d"])
		util_list.append(iso3d)
		vol3d = QtGui.QTreeWidgetItem(QtCore.QStringList("3D Volume Viewer"))
		vol3d.setIcon(0,self.icons["single_image_3d"])
		util_list.append(vol3d)
		self.launchers["3D Isosurface Viewer"] = self.launch_3d_iso_viewer
		slice3d = QtGui.QTreeWidgetItem(QtCore.QStringList("3D Slice Viewer"))
		slice3d.setIcon(0,self.icons["single_image_3d"])
		util_list.append(slice3d)
		self.launchers["3D Slice Viewer"] = self.launch_3d_slice_viewer
		image3d = QtGui.QTreeWidgetItem(QtCore.QStringList("3D Image Viewer"))
		image3d.setIcon(0,self.icons["single_image_3d"])
		util_list.append(image3d)
		self.launchers["3D Image Viewer"] = self.launch_3d_image_viewer
		image2d = QtGui.QTreeWidgetItem(QtCore.QStringList("2D Image Viewer"))
		image2d.setIcon(0,self.icons["single_image"])
		util_list.append(image2d)
		self.launchers["2D Image Viewer"] = self.launch_2d_image_viewer
		stack2d = QtGui.QTreeWidgetItem(QtCore.QStringList("2D Set Viewer"))
		stack2d.setIcon(0,self.icons["multiple_images"])
		util_list.append(stack2d)
		self.launchers["2D Set Viewer"] = self.launch_2d_stack_viewer
		
		flatform = QtGui.QTreeWidgetItem(QtCore.QStringList("Flat Form"))
		flatform.setIcon(0,self.icons["eman"])
		util_list.append(flatform)
		self.launchers["Flat Form"] = self.launch_flat_form
		
		tabform = QtGui.QTreeWidgetItem(QtCore.QStringList("Tab Form"))
		tabform.setIcon(0,self.icons["eman"])
		util_list.append(tabform)
		self.launchers["Tab Form"] = self.launch_table_form
		
		
		util.addChildren(util_list)
		return util
	
	def __get_generic_interfaces_tree(self,tree_widget):
		
		gi = QtGui.QTreeWidgetItem(QtCore.QStringList("Generic Interfaces"))

		gi_list = []
		gi_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Boxer")))
		gi_list[-1].setIcon(0,self.icons["green_boxes"])
		self.launchers["Boxer"] = self.launch_boxer_general
		gi_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("CTF ")))
		self.launchers["CTF "] = self.launch_ctf_general
		gi_list[-1].setIcon(0,self.icons["ctf"])
		gi_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Refine2D ")))
		self.launchers["Refine2D "] = self.launch_refine2d_general

		gi.addChildren(gi_list)
		return gi
	
	def __get_tomo_part_recon_tree(self,tree_widget):
		
		tomo = QtGui.QTreeWidgetItem(QtCore.QStringList("Tomographic Particle Reconstruction"))

		self.launchers["Tomographic Particle Reconstruction"] = self.launch_tomography
		
		tomo_list = []
		
		rd = QtGui.QTreeWidgetItem(QtCore.QStringList("Raw Tomogram Files"))
		tomo_list.append(rd)
		tomo_list[-1].setIcon(0,self.icons["single_image_3d"])
		
		self.launchers["Raw Tomogram Files"] = self.launch_tomo_raw_files
		
		ap = QtGui.QTreeWidgetItem(QtCore.QStringList("Tomogram Particles"))
		ap.setIcon(0,self.icons["green_boxes"])
		self.launchers["Tomogram Particles"] = self.launch_tomo_particle_report
		tomo_list.append(ap)
		
		filt = QtGui.QTreeWidgetItem(QtCore.QStringList("Filtered Particles"))
		filt.setIcon(0,self.icons["green_boxes"])
		self.launchers["Filtered Particles"] = self.launch_tomo_filt_ptcl_report
		tomo_list.append(filt)
		
		
		probes = QtGui.QTreeWidgetItem(QtCore.QStringList("Tomogram Probes"))
		tomo_list.append(probes)
		tomo_list[-1].setIcon(0,self.icons["black_box"])
		self.launchers["Tomogram Probes"] = self.launch_tomo_probe_report
		
		ali = QtGui.QTreeWidgetItem(QtCore.QStringList("Particle Alignment Report"))
		self.launchers["Particle Alignment Report"] = self.launch_tomo_ptcl_ali_report
		ali.setIcon(0,self.icons["tomo_hunter"])
		tomo_list.append(ali)
		
		aves = QtGui.QTreeWidgetItem(QtCore.QStringList("Particle Averages"))
		tomo_list.append(aves)
		tomo_list[-1].setIcon(0,self.icons["refine"])
		self.launchers["Particle Averages"] = self.launch_tomo_ave_report
		
		
		ap_list = []
		ap_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Interactive Boxing - e2tomoboxer")))
		self.launchers["Interactive Boxing - e2tomoboxer"] = self.launch_e2tomoboxer_gui
		ap_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Generate Output - e2tomoboxer")))
		self.launchers["Generate Output - e2tomoboxer"] = self.launch_e2tomoboxer_output
		ap_list[0].setIcon(0,self.icons["green_boxes"])
		ap_list[1].setIcon(0,self.icons["green_boxes"])
		ap.addChildren(ap_list)
		
		
		filt_list = []
		filt_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Filter/Rotate Tomo Particles")))
		filt_list[0].setIcon(0,self.icons["green_boxes"])
		self.launchers["Filter/Rotate Tomo Particles"] = self.launch_tomo_filter_ptcls
		filt_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Scale/Clip Tomo Particles")))
		filt_list[1].setIcon(0,self.icons["green_boxes"])
		self.launchers["Scale/Clip Tomo Particles"] = self.launch_tomo_scale_clip_ptcls
		filt_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Mask Particles")))
		filt_list[2].setIcon(0,self.icons["green_boxes"])
		self.launchers["Mask Particles"] = self.launch_tomo_mask_ptcls
		filt_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Norm/Invert Particles")))
		filt_list[3].setIcon(0,self.icons["green_boxes"])
		self.launchers["Norm/Invert Particles"] = self.launch_tomo_norm_inv_ptcls
		filt.addChildren(filt_list)
		
		rd_list = []
		rd_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Reconstruct ALI File")))
		self.launchers["Reconstruct ALI File"] = self.launch_reconstruct_ali
		rd_list[0].setIcon(0,self.icons["single_image_3d"])
		rd.addChildren(rd_list)
		
		probe_list = []
		probe_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Generate Bootstrapped Probe")))
		probe_list[-1].setIcon(0,self.icons["black_box"])
		self.launchers["Generate Bootstrapped Probe"] = self.launch_tomo_boostrap_probe
		probes.addChildren(probe_list)
		
		ali_list = []
		ali_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Tomohunter")))
		ali_list[-1].setIcon(0,self.icons["tomo_hunter"])
		self.launchers["Tomohunter"] = self.launch_tomohunter
		ali.addChildren(ali_list)
		
	
		
		aves_list = []
		aves_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Generate Average")))
		self.launchers["Generate Average"] = self.launch_gen_ave_tomo
		aves_list[0].setIcon(0,self.icons["refine"])
		aves.addChildren(aves_list)
		
		
		
		tomo.addChildren(tomo_list)
		
		tree_widget.expandItem(rd)
		tree_widget.expandItem(ap)
		tree_widget.expandItem(ali)
		tree_widget.expandItem(probes)
		tree_widget.expandItem(aves)
		tree_widget.expandItem(filt)
		return tomo
	
	def launch_reconstruct_ali(self):
		 self.launch_task(EMReconstructAliFile(),"Reconstruct ALI File")
		
	def task_killed(self,module_string,module):
		print module
		module.closeEvent(None)
		self.module().emit(QtCore.SIGNAL("module_closed"),"module_string",module)
	
	def display_file(self,filename):
		get_application().setOverrideCursor(Qt.BusyCursor)
		
		if len(filename) > 18 and filename[-19:] == "convergence.results":
			
			db = db_open_dict(filename,ro=True)
			keys = db.keys()
			res = get_e2resolution_results_list(keys)
			eo = get_e2eotest_results_list(keys)
			conv = get_convergence_results_list(keys)
			from emplot2d import EMPlot2DModule,colortypes
			module = EMPlot2DModule(get_application())
			i = 0
			max = len(colortypes)
			
			for k in conv:
				module.set_data(k,db[k],color=(i%max),linewidth=1) # there are only a ceratin number of  colors
				i += 1
			
			for plot in [eo,res]:
				for k in plot:
					module.set_data(k,db[k],color=(i%max),linewidth=3) # there are only a ceratin number of  colors
					i += 1
					
			get_application().show_specific(module)
			self.add_module([str(module),"Display",module])
		else:
			module = EMModuleFromFile(filename,get_application())
			if module != None:
				self.module().emit(QtCore.SIGNAL("launching_module"),"Browser",module)
				get_application().show_specific(module)
				self.add_module([str(module),"Display",module])
		get_application().setOverrideCursor(Qt.ArrowCursor)
		
	def launch_asym_unit(self):
		from e2eulerxplor import EMAsymmetricUnitViewer
		get_application().setOverrideCursor(Qt.BusyCursor)
		module = EMAsymmetricUnitViewer(get_application(),auto=False)
		self.module().emit(QtCore.SIGNAL("launching_module"),"Eulers Tool",module)
		get_application().show_specific(module)
		self.add_module([str(module),"Eulerxplor",module])
		get_application().setOverrideCursor(Qt.ArrowCursor)
		
	def launch_eulers(self):
		from e2eulerxplor import EMAsymmetricUnitViewer
		get_application().setOverrideCursor(Qt.BusyCursor)
		module = EMAsymmetricUnitViewer(get_application(),auto=True)
		self.module().emit(QtCore.SIGNAL("launching_module"),"Eulers",module)
		get_application().show_specific(module)
		self.add_module([str(module),"Eulerxplor",module])
		get_application().setOverrideCursor(Qt.ArrowCursor)
		
	def launch_lights_tool(self):
		from emlights import EMLights
		get_application().setOverrideCursor(Qt.BusyCursor)
		module = EMLights(application=em_app)
		self.module().emit(QtCore.SIGNAL("launching_module"),"EMLights",module)
		get_application().show_specific(module)
		self.add_module([str(module),"EMLights",module])
		get_application().setOverrideCursor(Qt.ArrowCursor)
	
	def launch_fonts_tool(self):
		from em3Dfonts import EM3DFontWidget
		get_application().setOverrideCursor(Qt.BusyCursor)
		module = EM3DFontWidget()
		self.module().emit(QtCore.SIGNAL("launching_module"),"EMLights",module)
		get_application().show_specific(module)
		self.add_module([str(module),"EMLights",module])
		get_application().setOverrideCursor(Qt.ArrowCursor)
	
	def launch_desktop(self):
		os.system("e2desktop.py")
	
	def launch_3d_volume_viewer(self):
		from emimage3dvol import EMVolumeModule
		self.launch_3d_viewer(EMVolumeModule(application=em_app),"3D Volumer Viewer")
	
	def launch_3d_iso_viewer(self):
		from emimage3diso import EMIsosurfaceModule
		self.launch_3d_viewer(EMIsosurfaceModule(application=em_app),"3D Isosurface Viewer")
		
	def launch_3d_slice_viewer(self):
		from emimage3dslice import EM3DSliceViewerModule
		self.launch_3d_viewer(EM3DSliceViewerModule(application=em_app),"3D Slice Viewer")
	
	def launch_3d_image_viewer(self):
		from emimage3d import EMImage3DModule
		self.launch_3d_viewer(EMImage3DModule(application=em_app),"3D Image Viewer")

	def launch_3d_viewer(self,module,name):
		from emimage3dvol import EMVolumeModule
		
		get_application().setOverrideCursor(Qt.BusyCursor)
		#module = EMVolumeModule(application=em_app)
		module.set_data(test_image_3d(1,size=(64,64,64)))
		self.module().emit(QtCore.SIGNAL("launching_module"),name,module)
		get_application().show_specific(module)
		self.add_module([str(module),name,module])
		get_application().setOverrideCursor(Qt.ArrowCursor)

	def launch_2d_image_viewer(self):
		get_application().setOverrideCursor(Qt.BusyCursor)
		from emimage2d import EMImage2DModule
		module = EMImage2DModule(application=em_app)
		module.set_data(test_image(Util.get_irand(0,10),size=(256,256)))
		self.module().emit(QtCore.SIGNAL("launching_module"),"2D Image Viewer",module)
		get_application().show_specific(module)
		self.add_module([str(module),"2D Image Viewer",module])
		get_application().setOverrideCursor(Qt.ArrowCursor)
		module.optimally_resize()
		
	def launch_2d_stack_viewer(self):
		get_application().setOverrideCursor(Qt.BusyCursor)
		from emimagemx import EMImageMXModule
		module = EMImageMXModule(application=em_app)
		data = [test_image(Util.get_irand(0,9)) for i in range(64)]
		module.set_data(data)
		self.module().emit(QtCore.SIGNAL("launching_module"),"2D Set Viewer",module)
		get_application().show_specific(module)
		self.add_module([str(module),"2D Set Viewer",module])
		get_application().setOverrideCursor(Qt.ArrowCursor)
		module.optimally_resize()
		
	def launch_flat_form(self):
		get_application().setOverrideCursor(Qt.BusyCursor)
		from emform import get_small_example_form_params,EMFormModule
		module = EMFormModule(params=get_small_example_form_params(),application=em_app)
		self.module().emit(QtCore.SIGNAL("launching_module"),"Flat Form",module)
		get_application().show_specific(module)
		self.add_module([str(module),"Flat Form",module])
		get_application().setOverrideCursor(Qt.ArrowCursor)
		
	def launch_table_form(self):
		get_application().setOverrideCursor(Qt.BusyCursor)
		from emform import get_example_table_form_params,EMTableFormModule
		module = EMTableFormModule(params=get_example_table_form_params(),application=em_app)
		self.module().emit(QtCore.SIGNAL("launching_module"),"Tabular Form",module)
		get_application().show_specific(module)
		self.add_module([str(module),"Tabular Form",module])
		get_application().setOverrideCursor(Qt.ArrowCursor)
	  	
	def launch_3d_plot_tool(self):
		from emplot3d import EMPlot3DModule,get_test_data,get_other_test_data
		get_application().setOverrideCursor(Qt.BusyCursor)
		module = EMPlot3DModule(application=em_app)
		module.set_data("test data",get_test_data())
		module.set_data("other data",get_other_test_data(),shape="Cube")
		self.module().emit(QtCore.SIGNAL("launching_module"),"Plot3D",module)
		get_application().show_specific(module)
		self.add_module([str(module),"Plot3D",module])
		get_application().setOverrideCursor(Qt.ArrowCursor)
	
	def launch_browser(self):
		import subprocess
		subprocess.Popen("e2display.py")
#		get_application().setOverrideCursor(Qt.BusyCursor)
#		module = EMBrowserModule()
#		self.module().emit(QtCore.SIGNAL("launching_module"),"Browser",module)
#		get_application().show_specific(module)
#		self.add_module([str(module),"Browse",module])
#		get_application().setOverrideCursor(Qt.ArrowCursor)
	
	def launch_refine2d_general(self):
		self.launch_task(E2Refine2DWithGenericTask(False),"e2refine2d Generic")
	
	def launch_ctf_general(self):
		self.launch_task(E2CTFGenericTask(),"e2ctf Generic")
	
	def launch_e2boxer_output(self):
		self.launch_task(E2BoxerOutputTask(),"e2boxer Output")
	
	def launch_boxer_general(self):
		self.launch_task(E2BoxerGenericTask(),"e2boxer General")
#		module = EMBoxerModule(get_application(),None)
#		self.module.emit(QtCore.SIGNAL("launching_module"),"Boxer",module)
#		module.show_guis()
#		self.add_module([str(module),"Boxer",module])
		
	def add_module(self,module):
		self.gui_modules.append(module)
		self.module_events.append(ModuleEventsManager(self,module[2]))
		self.update_task_list()
	
	def module_idle(self,module):
		# yes this is just the same as module_closed... thinking in progress
		for i,mod in enumerate(self.gui_modules):
			if mod[2] == module:
				self.gui_modules.pop(i)
				self.module_events.pop(i)
				self.update_task_list()
				return
			
		print "failed to handle idle?" # this shouldn't happen if I have managed everything correctly
	
	def module_closed(self,module):
		for i,mod in enumerate(self.gui_modules):
			if mod[2] == module:
				self.gui_modules.pop(i)
				self.module_events.pop(i)
				self.update_task_list()
				return
			
		print "failed to close module?" # this shouldn't happen if I have managed everything correctly

	def launch_e2resolution(self): self.launch_task(E2ResolutionTask(),"Run e2resolution")
	def launch_e2eotest(self): self.launch_task(E2EotestTask(),"Run e2eotest")
	def launch_resolution_report(self): self.launch_task(ResolutionReportTask(),"Resolution Report")
	def launch_e2refine_sets(self): self.launch_task(E2RefineChooseSetsTask(),"Choose Set For e2refine")
	def launch_e2refine(self): self.launch_task(E2RefineChooseParticlesTask(),"Choose Particles For e2refine")
	def launch_refinement_report(self): self.launch_task(RefinementReportTask(),"Refinement Report")
#	def launch_import_initial_model(self): self.launch_task(ImportInitialModels(),"import initial models")
	def launch_e2makeinitial(self): self.launch_task(E2InitialModel(),"e2makeinitialmodel")
	def launch_initmodel_report(self): self.launch_task(EMInitialModelReportTask(),"Initial Model Report")
	def launch_refine2d_report(self): self.launch_task(E2Refine2DReportTask(),"e2refine2d Report")
	def launch_refine2d_exec(self): self.launch_task(E2Refine2DRunTask(),"e2refine2d Parameters")
	def	launch_refine2d_choose_stacks(self): self.launch_task(E2Refine2DChooseSetsTask(),"Choose Sets For e2refine2d")
	def launch_refine2d_choose_ptcls(self): self.launch_task(E2Refine2DChooseParticlesTask(),"Choose Particles For e2refine2d")	
	def launch_e2ctf_write_ouptut(self): self.launch_task(E2CTFOutputTask(),"e2ctf Write Output")
	def launch_e2ctf_tune(self): self.launch_task(E2CTFGuiTask(),"e2ctf Intreface")
	def launch_e2ctf_auto_ft(self): self.launch_task(E2CTFAutoFitTask(),"e2ctf Auto Fitting")
	def launch_ctf_report(self):self.launch_task(CTFReportTask(),"CTF Report")
	def launch_mis_report(self): self.launch_task(EMSetReportTask(),"Project Sets")
	def launch_make_ptcl_set(self):	self.launch_task(E2MakeSetChooseDataTask(),"Build Particle Set")
	def launch_examine_particle_stacks(self): self.launch_task(E2ParticleExamineChooseDataTask(),"Examine Particles")
	def launch_particle_report(self): self.launch_task(EMParticleReportTask(),"Particle Report")
	def launch_tomo_particle_report(self): self.launch_task(EMTomoParticleReportTask(),"Tomogram Particles")
	def launch_tomo_probe_report(self): self.launch_task(EMTomoProbeReportTask(),"Tomogram Probes")
	def launch_tomo_ptcl_ali_report(self): self.launch_task(EMTomoChooseAlignedSetTask(),"Tomo Particle Alignment Report Report")
	def launch_tomo_ave_report(self): self.launch_task(EMTomoAveragesReportTask(),"Tomo Averages Report")
	def launch_tomo_boostrap_probe(self): error("This form is yet to be implemented")
	def launch_gen_ave_tomo(self): self.launch_task(EMTomoGenerateAverageChooseDataTask(),"Choose Alignment Set")
	def launch_tomo_filt_ptcl_report(self): self.launch_task(EMTomoChooseFilteredPtclsTask(),"Filtered Tomo Particles")
	def launch_tomo_filter_ptcls(self): self.launch_task(EMTomoChooseFilteredPtclsForFiltTask(),"Filter/Rotate Tomo Particles")
	def launch_tomo_scale_clip_ptcls(self):
		task = 	EMTomoChooseFilteredPtclsForFiltTask(E2TomoScaleClipTask)
		self.launch_task(task,"Scale/Clip Tomo Particles")
	def launch_tomo_mask_ptcls(self):
		task = EMTomoChooseFilteredPtclsForFiltTask( E2TomoMaskTask)
		self.launch_task(task,"Mask Tomo Particles")
	def launch_tomo_norm_inv_ptcls(self):
		task = EMTomoChooseFilteredPtclsForFiltTask( E2TomoNormalizeInvertTask)
		self.launch_task(task,"Norm/Invert Tomo Particles")
		
		
	def launch_particle_import(self):self.launch_task(EMParticleImportTask(),"Import Particles")
	
	def launch_ptcl_coord_import(self):self.launch_task(EMParticleCoordImportTask(),"Import Coordinate Data")
		
	def launch_mic_ccd_report(self): self.launch_task(EMRawDataReportTask(),"Raw Data")
		
	def launch_e2boxer_auto(self):
		self.launch_task(E2BoxerAutoTask(),"e2boxer Automated Boxing")
		
	def launch_e2boxer_gui(self):
		self.launch_task(E2BoxerGuiTask(),"e2boxer Interface")
	
	def launch_e2tomoboxer_gui(self):
		self.launch_task(E2TomoBoxerGuiTask(),"e2tomoboxer Interface")
	
	def launch_e2tomoboxer_output(self):
		self.launch_task(E2TomoBoxerOutputTask(),"Generate e2tomoboxer output")
	
	def launch_spr(self):
		self.launch_task(SPRInitTask(),"Single Particle Reconstruction")

	def launch_tomography(self):
		self.launch_task(EMTomoRawDataReportTask(),"Tomo Raw Data")
	def launch_tomohunter(self):
		self.launch_task(EMTomoChooseTomoHunterPtclsTask(),"Tomohunter")
	def launch_tomo_raw_files(self):
		self.launch_task(EMTomoRawDataReportTask(),"Tomo Raw Data")
			
#	def launch_mic_ccd(self):
#		self.launch_task(MicrographCCDTask(),"Micrograph/CCD report")
		
	def launch_view_history(self):
		self.launch_task(HistoryTask(),"History")
		
	def launch_view_preferences(self):
		from e2preferences import EMPreferencesTask
		self.launch_task(EMPreferencesTask(),"Preferences")
	
	def launch_filter_mic(self):
		self.launch_task(EMFilterRawDataTask(),"Filter Micrographs")
	
	def launch_change_directory(self):
		if len(self.tasks) > 0: 
			self.__clear_tasks()
			
		self.module().emit(QtCore.SIGNAL("task_selected"),"Workflow","Workflow")
		if not self.tasks.has_key("Change directory"):
			task = ChangeDirectoryTask()
			
			wd = task.run_form()
			if wd != None: self.directory.setToolTip(0,wd)
			
	def launch_task(self,task,task_unique_identifier):
		'''
		You can only have one of the task forms showing at any time
		'''
		get_application().setOverrideCursor(Qt.BusyCursor)
		if len(self.tasks) > 0: 
			self.__clear_tasks()
		
		self.module().emit(QtCore.SIGNAL("task_selected"),"Workflow","Workflow")
		if not self.tasks.has_key(task_unique_identifier):
			#task = task_type()
			
			task.run_form()
			self.tasks[task_unique_identifier] = task
			self.event_managers[task_unique_identifier] = TaskEventsManager(task,self,task_unique_identifier)
		else:
			get_application().show_specific(self.tasks[task_unique_identifier])
			
		self.update_task_list()
		get_application().setOverrideCursor(Qt.ArrowCursor)
	def update_task_list(self):
		tasks = []
		tasks_dict = {}
		
		if len(self.tasks) > 1:
			print "error, I haven't been written to handle more than one task"
			# solution is to turn the tasks into a module and add more modules
			return
		
		for val in self.tasks.keys():
			tasks.append(str(val))
			tasks_dict["Workflow"] = val
			#tasks_dict["workflow"] = 
		
		for val in self.gui_modules:
			tasks_dict[val[2]] = val[1]
			tasks.append(str(val[1]))
		
		self.module().emit(QtCore.SIGNAL("tasks_updated"),tasks_dict)
	
	def __clear_tasks(self):
		for v in self.tasks.values():
			v.closeEvent(None)
			#get_application().close_specific(v)
		self.tasks = {}
		self.event_managers = {}
	
	
	def on_replace_task(self,old_task,module_task,task_name):
		self.pop_task_event_pair(old_task)
		self.launch_task(module_task,task_name)
		
	
	def pop_task_event_pair(self,task):
		old_task = self.tasks.pop(task)
		self.event_managers.pop(task)
		self.update_task_list()
		
	def gui_running(self,task_key,module_string_name,module):
		'''
		Tells this to take the task out of the dictionary of current tasks and put in another dictionary.
		Essentially saying, "make this task disappear, but keep a reference to it until the associated gui ends"
		'''
		
		self.module().emit(QtCore.SIGNAL("launching_module"),module_string_name,module) # for the desktop to prepare a cube
		module = [self.tasks.pop(task_key),module_string_name,module]
		self.gui_modules.append(module)
		self.module_events.append(None) # this makes the handling of the gui_modules and module_events more general, see code
		self.update_task_list()
	
	def gui_exit(self,module):
		for i,mod in enumerate(self.gui_modules):
			if mod[0] == module:
				#print "gui exit succeeded",module
				self.gui_modules.pop(i)
				self.module_events.pop(i)
				self.update_task_list()
				return
			
		print "gui exit failed",module #this shouldn't happen if I have managed everything correctly
	
	def tree_widget_click(self,tree_item,i):
		task = str(tree_item.text(0))
		if task in self.launchers.keys():
			self.launchers[task]()
		else:
			if task == "Browse":
				print "hi"
				
	def keyPressEvent(self,event):
		if event.key() == Qt.Key_F1:
			try:
				import webbrowser
				webbrowser.open("http://blake.bcm.edu/emanwiki/e2workflow")
				return
			except: pass
			
			try: from PyQt4 import QtWebKit
			except: return
			try:
				try:
					test = self.browser
				except: 
					self.browser = QtWebKit.QWebView()
					self.browser.load(QtCore.QUrl("http://blake.bcm.edu/emanwiki/e2workflow"))
					self.browser.resize(800,800)
				
				if not self.browser.isVisible(): self.browser.show()
			except: pass
	
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
		self.selector = weakref.ref(selector)
		self.key = key
		QtCore.QObject.connect(self.task.emitter(), QtCore.SIGNAL("task_idle"), self.on_task_idle) # this typically gets emitted when the user hits ok or cancel on the 
		QtCore.QObject.connect(self.task.emitter(), QtCore.SIGNAL("gui_running"),self.on_gui_running) # this one 
		QtCore.QObject.connect(self.task.emitter(), QtCore.SIGNAL("replace_task"),self.on_replace_task)
		QtCore.QObject.connect(self.task.emitter(), QtCore.SIGNAL("gui_exit"),self.on_gui_exit)
		QtCore.QObject.connect(self.task.emitter(), QtCore.SIGNAL("process_started"), self.on_process_started)
		QtCore.QObject.connect(self.task.emitter(), QtCore.SIGNAL("display_file"), self.selector().display_file)
#		
	def on_gui_running(self,module_string_name,module):
		''' 
		
		'''
		self.selector().gui_running(self.key,module_string_name,module)
		
	def on_gui_exit(self):
		self.selector().gui_exit(self.task)
		
	def on_task_idle(self):
		self.selector().pop_task_event_pair(self.key)
		
	def on_process_started(self,pid):
		self.selector().task_monitor.add_process(pid)
	
	def on_replace_task(self,task,task_name):
		self.selector().on_replace_task(self.key,task,task_name)

class EMWorkFlowManager:
	def __init__(self,application):
		self.application = weakref.ref(application)
		
	
		self.task_monitor = EMTaskMonitorModule(get_application())
		self.selector = EMWorkFlowSelector(application,self.task_monitor.qt_widget)
		self.selector.qt_widget.resize(300,540)
		QtCore.QObject.connect(self.selector.emitter(), QtCore.SIGNAL("tasks_updated"),self.task_monitor.qt_widget.set_entries)
		QtCore.QObject.connect(self.task_monitor.emitter(), QtCore.SIGNAL("task_killed"),self.selector.qt_widget.task_killed)
		
	
	
	def close(self):
		self.selector().closeEvent(None)
		self.task_monitor.closeEvent(None)
		#get_application().close_specific(self.selector)
		#get_application().close_specific(self.task_monitor)

if __name__ == '__main__':
	logid=E2init(sys.argv)
	from emapplication import EMStandAloneApplication
	em_app = EMStandAloneApplication()
	sprinit = EMWorkFlowManager(em_app)
	
	em_app.show()
	em_app.execute()
	E2end(logid)