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
from emapplication import ModuleEventsManager
from emsprworkflow import *
from emtprworkflow import *
from emselector import EMBrowser
from e2boxer import EMBoxerModule
from EMAN2 import process_running,kill_process
from EMAN2db import db_open_dict, EMTask
from emanimationutil import Animator
from emglobjects import EM3DGLWidget
from emimage3d import EMImage3DWidget
from emimage import EMWidgetFromFile

import time
import weakref
import traceback

class EMTaskMonitorWidget(QtGui.QWidget):
	def __init__(self):
		QtGui.QWidget.__init__(self,None)
		self.setWindowIcon(QtGui.QIcon(get_image_directory() + "workflow.png"))
		self.setWindowTitle("Running Tasks")
		
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
		
#		self.vbl.addWidget(self.kill)
		
#		self.set_entries(self.entries_dict)
	
#		QtCore.QObject.connect(self.kill, QtCore.SIGNAL("clicked()"), self.on_kill)
		
		self.resize(256,200)
		get_application().attach_child(self)
		
		self.tasks=None
		
		# A timer for updates
		self.timer = QtCore.QTimer(self);
 		QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.update_list)
 		self.timer.start(3000)
	
	def check_task(self,fin,ptsk):
		"""Note that this modifies ptsk in-place"""
		fin.seek(ptsk[0])
		try : t=fin.readline()
		except: return False
		
		tsk=t.split("\t")
		ptsk[3]=tsk[1]
		if self.is_running(tsk) : return True
		
		return False
	
	def is_running(self,tsk):
		"""Check to see if the 'tsk' string is actively running"""
		
		try:
			# This means the task must be done
			if tsk[1][4]=="/" : return False
			
			# or not running
			if '/' in tsk[2] : pid=int(tsk[2].split("/")[0])
			else : pid=int(tsk[2])
			
			if os.name=="posix":
				# don't do this on windows (may actually kill the process)
				os.kill(pid,0)			# raises an exception if doesn't exist
			else :
				# not a good approach, but the best we have right now on windows
				tm=time.localtime()
				mon=int(tsk[0][5:7])
				yr=int(tsk[0][:4])
				day=int(tsk[0][8:10])
				if mon!=tm[1] or yr!=tm[0] or tm[2]-day>1 : raise Exception 
				
			if os.name=="nt": tsk[4]=tsk[4].convert("\\","/")
			command = tsk[4].split()[0].split("/")[-1]
			if command=="e2workflow.py" : raise Exception
		except:
#			traceback.print_exc()
			return False

		return True
	
	def parse_new_commands(self,fin):
		self.last_update=time.time()
		while 1:
			loc=fin.tell()
			try: 
				t=fin.readline()
				if len(t)==0 : raise Exception
			except: 
				self.lastpos=loc
				break
				
			tsk=t.split("\t")
			if not self.is_running(tsk) : continue
			
			# if we get here, then we have a (probably) active task
			
			# This contains (file location,process id, status, command name)
			if '/' in tsk[2] : pid=int(tsk[2].split("/")[0])
			else : pid=int(tsk[2])
			if os.name=="nt": tsk[4]=tsk[4].convert("\\","/")
			command = tsk[4].split()[0].split("/")[-1]
			self.tasks.append([loc,pid,tsk[0],tsk[1],command])
		
	
	def update_list(self):
		
		if self.tasks==None:
			try: fin=file(".eman2log.txt","r")
			except: return
			
			self.tasks=[]
			self.parse_new_commands(fin)
			if len(self.tasks)==0 :
				
				#print "no commands"
				self.tasks=None
				return
		
		else:
			# only check the file when it changes
			try: 
				if os.stat(".eman2log.txt").st_mtime < self.last_update : return
			except:
				print "Error, couldn't stat(.eman2log.txt)."
				return
				
			try: fin=file(".eman2log.txt","r")
			except: return
			
			# Go through existing entries and see if they've finished
			self.tasks=[i for i in self.tasks if self.check_task(fin,i)]
			
			# now check for new entries
			fin.seek(self.lastpos)
			self.parse_new_commands(fin)
				
		self.list_widget.clear()
		items=["%s  %s(%s)  %s"%(t[2][5:16],t[3][:4],t[1],t[4]) for t in self.tasks]
		self.list_widget.addItems(QtCore.QStringList(items))
	
	def animate(self,time):
		print "animate"
		return True
				
			
	def set_entries(self,entries_dict):
		print "set_entries",entries_dict
			
	def list_processes(self,s_text):
		print "list_processes ",s_text

	def accrue_process_info(self):
		print "accrue_process_info"
	
	def add_process(self,pid):
		print "add_process ",pid
	
	def list_widget_item_clicked(self,item):
		print "list_widget_item_clicked ",item.text()
#		self.emit(QtCore.SIGNAL("task_selected"),str(item.text()),item.module)
	
	def on_kill(self):
		print "on_kill"
		selected_items = self.list_widget.selectedItems()
		if len(selected_items) == 0:
			return
		
		# not sure if there will ever be more than one selected but oh well let's support closing them all
		
		#for item in selected_items:
			#if item.module == "process":
				#if not kill_process(int(item.pid)):
					#print "kill process failed for some reason"
			#else:
				#self.emit(QtCore.SIGNAL("task_killed"),str(item.text()),item.module)
	
class EMWorkFlowSelectorWidget(QtGui.QWidget):
	def __init__(self,task_monitor):
		'''
		task_monitor should probably be supplied, it should be an instance of a EMTaskMonitorWidget
		'''
		QtGui.QWidget.__init__(self,None)
		self.__init_icons()
		self.setWindowIcon(self.icons["desktop"])
		self.setWindowTitle("EMAN2 Tasks")
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

		get_application().attach_child(self)
	def __init_icons(self):
		self.icons = {}
		print get_image_directory()
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
		ctf_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Generate Structure Factor - e2ctf")))
		self.launchers["Generate Structure Factor - e2ctf"] =  self.launch_e2ctf_write_sf
		ctf_list[0].setIcon(0,self.icons["ctf"])
		ctf_list[1].setIcon(0,self.icons["ctf"])
		ctf_list[2].setIcon(0,self.icons["ctf"])
		ctf_list[3].setIcon(0,self.icons["ctf"])
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
		mis_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Evaluate Particle Sets")))
		self.launchers["Evaluate Particle Sets"] = self.launch_evaluate_ptcl_set
		mis_list[0].setIcon(0,self.icons["multiple_images"])
		mis_list[1].setIcon(0,self.icons["multiple_images"])
		mis.addChildren(mis_list)
		
		refine2d_list = []
		refine2d_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Generate Classes - e2refine2d")))
		refine2d_list[-1].setIcon(0,self.icons["classes"])
		self.launchers["Generate Classes - e2refine2d"] = self.launch_refine2d_choose_stacks
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
		refine_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Run - e2refine")))
		refine_list[-1].setIcon(0,self.icons["refine"])
		self.launchers["Run - e2refine"] = self.launch_e2refine_sets
		refine_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Run - e2refinemulti")))
		refine_list[-1].setIcon(0,self.icons["refine"])
		self.launchers["Run - e2refinemulti"] = self.launch_e2refinemulti_sets
		
		freealign_list = []
		freealign = QtGui.QTreeWidgetItem(QtCore.QStringList("Frealign"))
		freealign.setIcon(0,self.icons["refine"])
		freealign_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Run - e2refinetofrealign")))
		freealign_list[-1].setIcon(0,self.icons["refine"])
		self.launchers["Run - e2refinetofrealign"] = self.launch_e2refinetofrealign
		freealign_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Run - e2runfrealign")))
		freealign_list[-1].setIcon(0,self.icons["refine"])
		self.launchers["Run - e2runfrealign"] = self.launch_e2runfrealign
		freealign.addChildren(freealign_list)
		freealign_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Run - e2refinefromfrealign")))
		freealign_list[-1].setIcon(0,self.icons["refine"])
		self.launchers["Run - e2refinefromfrealign"] = self.launch_e2refinefromfrelaign
		freealign.addChildren(freealign_list)
		
		refine_list.append(freealign)
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
		rd.setIcon(0,self.icons["single_image_3d"])
		self.launchers["Raw Tomogram Files"] = self.launch_tomography
		tomo_list.append(rd)
		
		ac = QtGui.QTreeWidgetItem(QtCore.QStringList("Box Tomogram Particles"))
		ac.setIcon(0,self.icons["green_boxes"])
		self.launchers["Box Tomogram Particles"] = self.launch_tomo_boxer
		tomo_list.append(ac)
		
		ali = QtGui.QTreeWidgetItem(QtCore.QStringList("Tomogram Alignment and Averaging"))
		ali.setIcon(0,self.icons["tomo_hunter"])
		self.launchers["Tomogram Alignment and Averaging"] = self.launch_tomo_hunter
		tomo_list.append(ali)
		
		tomo.addChildren(tomo_list)
		
		tree_widget.expandItem(ac)
		tree_widget.expandItem(ali)
		
		return tomo
	
	def launch_reconstruct_ali(self):
		 self.launch_task(EMReconstructAliFile(),"Reconstruct ALI File")
		
	def task_killed(self,module_string,module):
		print module
		module.close()
		self.emit(QtCore.SIGNAL("module_closed"),"module_string",module)
	
	def display_file(self,filename):
		get_application().setOverrideCursor(Qt.BusyCursor)
		
		if len(filename) > 18 and filename[-19:] == "convergence.results":
			
			db = db_open_dict(filename,ro=True)
			keys = db.keys()
			res = get_e2resolution_results_list(keys)
			eo = get_e2eotest_results_list(keys)
			conv = get_convergence_results_list(keys)
			from emplot2d import EMPlot2DWidget,colortypes
			module = EMPlot2DWidget(get_application())
			i = 0
			max = len(colortypes)
			
			for k in conv:
				module.set_data(db[k],k,color=(i%max),linewidth=1) # there are only a ceratin number of  colors
				i += 1
			
			for plot in [eo,res]:
				for k in plot:
					module.set_data(db[k],k,color=(i%max),linewidth=3) # there are only a ceratin number of  colors
					i += 1
					
			module.show()
			self.add_module([str(module),"Display",module])
		else:
			module = EMWidgetFromFile(filename,get_application())
			if module != None:
				self.emit(QtCore.SIGNAL("launching_module"),"Browser",module)
				get_application().show_specific(module)
				self.add_module([str(module),"Display",module])
		get_application().setOverrideCursor(Qt.ArrowCursor)
		
	def launch_asym_unit(self):
		from emimage3dsym import EMSymViewerWidget
		get_application().setOverrideCursor(Qt.BusyCursor)
		widget = EMSymViewerWidget()
		self.emit(QtCore.SIGNAL("launching_module"),"Eulers Tool",widget)
		get_application().show_specific(widget)
		self.add_module([str(widget),"Eulerxplor",widget])
		get_application().setOverrideCursor(Qt.ArrowCursor)
		
	def launch_eulers(self):
		from e2eulerxplor import EMEulerWidget
		get_application().setOverrideCursor(Qt.BusyCursor)
		widget = EMEulerWidget()
		self.emit(QtCore.SIGNAL("launching_module"),"Eulers",widget)
		get_application().show_specific(widget)
		self.add_module([str(widget),"Eulerxplor",widget])
		get_application().setOverrideCursor(Qt.ArrowCursor)
		
	def launch_lights_tool(self):
		from emlights import EMLights
		from emglobjects import EM3DGLWidget
		get_application().setOverrideCursor(Qt.BusyCursor)
		window = EM3DGLWidget()
		em_lights = EMLights(window)
		window.set_model(em_lights)
		self.emit(QtCore.SIGNAL("launching_module"),"EMLights",window)
		get_application().show_specific(window)
		self.add_module([str(window),"EMLights",window])
		get_application().setOverrideCursor(Qt.ArrowCursor)
	
	def launch_fonts_tool(self):
		from em3Dfonts import EM3DFontWidget
		get_application().setOverrideCursor(Qt.BusyCursor)
		module = EM3DFontWidget()
		self.emit(QtCore.SIGNAL("launching_module"),"EMLights",module)
		get_application().show_specific(module)
		self.add_module([str(module),"EMLights",module])
		get_application().setOverrideCursor(Qt.ArrowCursor)
	
	def launch_3d_volume_viewer(self):
		from emimage3dvol import EMVolumeModel
		widget = EM3DGLWidget()
		self.launch_3d_viewer(EMVolumeModel(widget),"3D Volumer Viewer")
	
	def launch_3d_iso_viewer(self):
		from emimage3diso import EMIsosurfaceModel
		widget = EM3DGLWidget()
		self.launch_3d_viewer(EMIsosurfaceModel(widget),"3D Isosurface Viewer")
		
	def launch_3d_slice_viewer(self):
		from emimage3dslice import EM3DSliceModel
		widget = EM3DGLWidget()
		self.launch_3d_viewer(EM3DSliceModel(widget),"3D Slice Viewer")
	
	def launch_3d_image_viewer(self):
		widget = EMImage3DWidget(application=em_app)
		
		get_application().setOverrideCursor(Qt.BusyCursor)
		data = test_image_3d(1,size=(64,64,64))
		widget.set_data(data)
		name = "3D Image Viewer"
		self.emit(QtCore.SIGNAL("launching_module"), name, widget)
		get_application().show_specific(widget)
		self.add_module([str(widget),name,widget])
		get_application().setOverrideCursor(Qt.ArrowCursor)

	def launch_3d_viewer(self,model,name):
		get_application().setOverrideCursor(Qt.BusyCursor)
		widget = model.get_gl_widget()
		widget.set_model(model)
		data = test_image_3d(1,size=(64,64,64))
		model.set_data(data)
		widget.set_data(data)
		self.emit(QtCore.SIGNAL("launching_module"),name,widget)
		get_application().show_specific(widget)
		self.add_module([str(widget),name,widget])
		get_application().setOverrideCursor(Qt.ArrowCursor)

	def launch_2d_image_viewer(self):
		get_application().setOverrideCursor(Qt.BusyCursor)
		from emimage2d import EMImage2DWidget
		module = EMImage2DWidget(application=em_app)
		module.set_data(test_image(Util.get_irand(0,10),size=(256,256)))
		self.emit(QtCore.SIGNAL("launching_module"),"2D Image Viewer",module)
		get_application().show_specific(module)
		self.add_module([str(module),"2D Image Viewer",module])
		get_application().setOverrideCursor(Qt.ArrowCursor)
		module.optimally_resize()
		
	def launch_2d_stack_viewer(self):
		get_application().setOverrideCursor(Qt.BusyCursor)
		from emimagemx import EMImageMXWidget
		module = EMImageMXWidget(application=em_app)
		data = [test_image(Util.get_irand(0,9)) for i in range(64)]
		module.set_data(data)
		self.emit(QtCore.SIGNAL("launching_module"),"2D Set Viewer",module)
		get_application().show_specific(module)
		self.add_module([str(module),"2D Set Viewer",module])
		get_application().setOverrideCursor(Qt.ArrowCursor)
		module.optimally_resize()
		
	def launch_flat_form(self):
		get_application().setOverrideCursor(Qt.BusyCursor)
		from emform import get_small_example_form_params,EMFormWidget
		module = EMFormWidget(params=get_small_example_form_params())
		self.emit(QtCore.SIGNAL("launching_module"),"Flat Form",module)
		get_application().show_specific(module)
		self.add_module([str(module),"Flat Form",module])
		get_application().setOverrideCursor(Qt.ArrowCursor)
		
	def launch_table_form(self):
		get_application().setOverrideCursor(Qt.BusyCursor)
		from emform import get_example_table_form_params,EMTableFormWidget
		module = EMTableFormWidget(params=get_example_table_form_params())
		self.emit(QtCore.SIGNAL("launching_module"),"Tabular Form",module)
		get_application().show_specific(module)
		self.add_module([str(module),"Tabular Form",module])
		get_application().setOverrideCursor(Qt.ArrowCursor)
	  	
	def launch_3d_plot_tool(self):
		from emplot3d import EMPlot3DWidget,get_test_data,get_other_test_data
		get_application().setOverrideCursor(Qt.BusyCursor)
		module = EMPlot3DWidget()
		module.set_data(get_test_data(),"test data")
		module.set_data(get_other_test_data(),"other data",shape="Cube")
		#module.set_data("test data",get_test_data())
		#module.set_data("other data",get_other_test_data(),shape="Cube")
		self.emit(QtCore.SIGNAL("launching_module"),"Plot3D",module)
		get_application().show_specific(module)
		self.add_module([str(module),"Plot3D",module])
		get_application().setOverrideCursor(Qt.ArrowCursor)
	
	def launch_browser(self):
#		import subprocess
#		subprocess.Popen("e2display.py")
		get_application().setOverrideCursor(Qt.BusyCursor)
		module = EMBrowser()
		self.emit(QtCore.SIGNAL("launching_module"),"Browser",module)
		get_application().show_specific(module)
		self.add_module([str(module),"Browse",module])
		get_application().setOverrideCursor(Qt.ArrowCursor)
	
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
	def launch_e2refinemulti_sets(self): self.launch_task(E2RefineMultiChooseSetsTask(),"Choose Set For e2refinemulti")
	def launch_e2refinetofrealign(self): self.launch_task(E2RefineToFreeAlign(),"Run RefineToFreAlign")
	def launch_e2runfrealign(self): self.launch_task(E2RunFreAlign(), "Run RunFreAlign")
	def launch_e2refinefromfrelaign(self): self.launch_task(E2RefineFromFreAlign(), "Run RefineFromFreAlign")
	def launch_e2refine(self): self.launch_task(E2RefineChooseParticlesTask(),"Choose Particles For e2refine")
	def launch_refinement_report(self): self.launch_task(RefinementReportTask(),"Refinement Report")
	def launch_import_initial_model(self): self.launch_task(ImportInitialModels(),"import initial models")
	def launch_e2makeinitial(self): self.launch_task(E2InitialModel(),"e2makeinitialmodel")
	def launch_initmodel_report(self): self.launch_task(EMInitialModelReportTask(),"Initial Model Report")
	def launch_refine2d_report(self): self.launch_task(E2Refine2DReportTask(),"e2refine2d Report")
	def launch_refine2d_exec(self): self.launch_task(E2Refine2DRunTask(),"e2refine2d Parameters")
	def	launch_refine2d_choose_stacks(self): self.launch_task(E2Refine2DChooseSetsTask(),"Choose Sets For e2refine2d")
	def launch_refine2d_choose_ptcls(self): self.launch_task(E2Refine2DChooseParticlesTask(),"Choose Particles For e2refine2d")	
	def launch_e2ctf_write_ouptut(self): self.launch_task(E2CTFOutputTask(),"e2ctf Write Output")
	def launch_e2ctf_write_sf(self): self.launch_task(E2CTFSFOutputTask(),"e2ctf Structure Factor")
	def launch_e2ctf_tune(self): self.launch_task(E2CTFGuiTask(),"e2ctf Intreface")
	def launch_e2ctf_auto_ft(self): self.launch_task(E2CTFAutoFitTask(),"e2ctf Auto Fitting")
	def launch_ctf_report(self):self.launch_task(CTFReportTask(),"CTF Report")
	def launch_mis_report(self): self.launch_task(EMSetReportTask(),"Project Sets")
	def launch_make_ptcl_set(self):	self.launch_task(E2MakeSetChooseDataTask(),"Build Particle Set")
	def launch_evaluate_ptcl_set(self):	self.launch_task(E2EvaluateSetTask(),"Build Particle Set")
	def launch_examine_particle_stacks(self): self.launch_task(E2ParticleExamineChooseDataTask(),"Examine Particles")
	def launch_particle_report(self): self.launch_task(EMParticleReportTask(),"Particle Report")	
	def launch_particle_import(self):self.launch_task(EMParticleImportTask(),"Import Particles")
	
	def launch_tomography(self):
		self.launch_task(EMTomoRawDataReportTask(),"Tomo Raw Data")
	def launch_tomo_boxer(self):
		self.launch_task(E2TomoBoxerGuiTask(),"Box Tomograms")
	def launch_tomo_hunter(self): 
		self.launch_task(EMTomoBootstrapTask(),"Tomo Hunter")
		
	def launch_ptcl_coord_import(self):self.launch_task(EMParticleCoordImportTask(),"Import Coordinate Data")
		
	def launch_mic_ccd_report(self): self.launch_task(EMRawDataReportTask(),"Raw Data")
		
	def launch_e2boxer_auto(self):
		self.launch_task(E2BoxerAutoTask(),"e2boxer Automated Boxing")
		
	def launch_e2boxer_gui(self):
		self.launch_task(E2BoxerGuiTask(),"e2boxer Interface")
	
	def launch_spr(self):
		self.launch_task(SPRInitTask(),"Single Particle Reconstruction")

	def launch_mic_ccd(self):
		self.launch_task(MicrographCCDTask(),"Micrograph/CCD report")
		
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
			
		self.emit(QtCore.SIGNAL("task_selected"),"Workflow","Workflow")
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
		
		self.emit(QtCore.SIGNAL("task_selected"),"Workflow","Workflow")
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
		
		for (task_id, task) in self.tasks.items():
			tasks.append(str(task_id))
			#tasks_dict["Workflow"] = task_id
			tasks_dict[task] = task_id
		
		for val in self.gui_modules:
			tasks_dict[val[2]] = val[1]
			tasks.append(str(val[1]))
		
		self.emit(QtCore.SIGNAL("tasks_updated"),tasks_dict)
	
	def __clear_tasks(self):
		for v in self.tasks.values():
			v.close()
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
		
		self.emit(QtCore.SIGNAL("launching_module"),module_string_name,module) # for the desktop to prepare a cube
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
		
	
		self.task_monitor = EMTaskMonitorWidget()
		self.selector = EMWorkFlowSelectorWidget(self.task_monitor)
		self.selector.resize(300,540)
		#QtCore.QObject.connect(self.selector, QtCore.SIGNAL("tasks_updated"),self.task_monitor.set_entries) # I don't like all this output!
		QtCore.QObject.connect(self.task_monitor, QtCore.SIGNAL("task_killed"),self.selector.task_killed)
		
	
	
	def close(self):
		self.selector().close()
		self.task_monitor.close()
		#get_application().close_specific(self.selector)
		#get_application().close_specific(self.task_monitor)

if __name__ == '__main__':
	logid=E2init(sys.argv)
	from emapplication import EMApp
	em_app = EMApp()
	sprinit = EMWorkFlowManager(em_app)
	
	em_app.show()
	em_app.execute()
	E2end(logid)
