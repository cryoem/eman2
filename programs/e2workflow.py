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
from emapplication import EMQtWidgetModule
from emsprworkflow import *

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
		
		spr = QtGui.QTreeWidgetItem(QtCore.QStringList("SPR"))
		self.launchers["SPR"] = self.launch_spr
		
		self.tree_widget_entries.append(spr)
		self.tree_widget_entries.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Browse")))
		self.tree_widget_entries.append(QtGui.QTreeWidgetItem(QtCore.QStringList("Boxer")))
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
		ctf_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("e2ctf")))
		ctf.addChildren(ctf_list)
		self.launchers["e2ctf"] = self.launch_e2ctf_management
		
		
		ap_list = []
		ap_list.append(QtGui.QTreeWidgetItem(QtCore.QStringList("e2boxer")))
		self.launchers["e2boxer"] = self.launch_e2boxer_management
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
		
		QtCore.QObject.connect(self.tree_widget, QtCore.SIGNAL("itemDoubleClicked(QTreeWidgetItem*,int)"), self.tree_widget_double_click)
		QtCore.QObject.connect(self.tree_widget, QtCore.SIGNAL("itemClicked(QTreeWidgetItem*,int)"), self.tree_widget_click)
		QtCore.QObject.connect(self.close, QtCore.SIGNAL("clicked()"), self.target.close)
	
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
		if len(self.tasks) > 0: 
			self.__clear_tasks()
		
		if not self.tasks.has_key(task_unique_identifier):
			task = task_type(self.application)
			
			task.run_form()
			self.tasks[task_unique_identifier] = task
			self.event_managers[task_unique_identifier] = TaskEventsManager(task,self,task_unique_identifier)
		else:
			self.application.show_specific(self.tasks[task_unique_identifier])
	
	def __clear_tasks(self):
		for v in self.tasks.values():
			v.closeEvent(None)
			#self.application.close_specific(v)
		self.tasks = {}
		self.event_managers = {}
	
	def pop_task_event_pair(self,task):
		self.tasks.pop(task)
		self.event_managers.pop(task)
		
	def tree_widget_click(self,tree_item,i):
		task = str(tree_item.text(0))
		if task in self.launchers.keys():
			self.launchers[task]()
		else:
			if task == "Browse":
				print "hi"
	
	def tree_widget_double_click(self,tree_item,i):
		task = tree_item.text(0)
		if task == "Browse":
			self.target.add_browser_frame()
		if task == "Thumb":
			self.target.add_selector_frame()
		elif task == "Boxer":
			self.target.add_boxer_frame()


class TaskEventsManager:
	def __init__(self,task,selector,key):
		self.task = task
		self.selector = selector
		self.key = key
		QtCore.QObject.connect(self.task, QtCore.SIGNAL("task_idle"), self.on_task_idle)
	
	def on_task_idle(self):
		print "recieved task idle"
		self.selector.pop_task_event_pair(self.key)
		

class EMWorkFlowManager:
	def __init__(self,application):
		self.application = application
		self.selector = EMWorkFlowSelector(self,application)
	
	def close(self):
		self.application.close_specific(self.selector)

if __name__ == '__main__':
	
	from emapplication import EMStandAloneApplication
	em_app = EMStandAloneApplication()
	sprinit = EMWorkFlowManager(em_app)
	#window = sprinit.run_form() 
	
	#em_app.show()
	em_app.execute()	