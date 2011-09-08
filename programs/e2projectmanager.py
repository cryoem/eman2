#!/usr/bin/env python
#
# Author: John Flanagan (jfflanag@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine


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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import Qt
import os, json
from EMAN2db import db_open_dict

class EMProjectManager(QtGui.QMainWindow):
	def __init__(self):
		QtGui.QMainWindow.__init__(self)
		# default PM attributes
		self.pm_icon = "../images/EMAN2Icon.tif"
		self.pm_cwd = os.getcwd()
		self.usingEMEN = False
		# Load the project DataBase
		self.loadPMdb()
		# Load icons
		self._load_icons()
		
		# Make the project manager layout
		self.makeMenues()
		font = QtGui.QFont()
		font.setBold(True)
		centralwidget = QtGui.QWidget()
		grid = QtGui.QGridLayout()
		grid.addWidget(self.makeTilteBarWidget(), 0, 0, 1, 2)
		workflowcontrollabel = QtGui.QLabel("Workflow Control", centralwidget)
		workflowcontrollabel.setFont(font)
		grid.addWidget(workflowcontrollabel, 1,0)
		guilabel = QtGui.QLabel("EMAN2 program GUI", centralwidget)
		guilabel.setFont(font)
		grid.addWidget(guilabel, 1,1)
		grid.addWidget(self.makeStackedWidget(),2,0)
		grid.addWidget(self.makeStackedGUIwidget(),2,1,2,1)
		# Make Browse button
		self.browsebutton = QtGui.QPushButton("Browse")
		grid.addWidget(self.browsebutton, 3, 0)
		# Make status bar
		self.statusbar = EMAN2StatusBar("Welcome to the EMAN2 Project Manager")
		grid.addWidget(self.statusbar, 4,0,1,2)
		centralwidget.setLayout(grid)
		# Make status bar
		self.setCentralWidget(centralwidget)
		# Make logbook
		self.logbook = LogBook()
		self.logbook.setGeometry(400,400,400,400)
		self.logbook.show()
		
		QtCore.QObject.connect(self.browsebutton,QtCore.SIGNAL("clicked()"),self._on_browse)
		
	def _on_browse(self):
		print "Browse"
		#self.pm_icon = "/home/john/Desktop/MmCpnValidate.ti"
		
	def closeEvent(self, event):
		self.logbook.close()
	
	def loadPMdb(self):
		# This probably wont work on WINDOWS.....
		self.pm_projects_db = db_open_dict("bdb:"+os.environ['HOME']+"#:pm_projects")
		# Load the project if we are already in the right dir
		self.pn_project_name = "Unknown"
		for project in self.pm_projects_db.keys():
			if self.pm_projects_db[project]["CWD"] == self.pm_cwd:
				self.pn_project_name = project
		
	def makeMenues(self):
		"""
		Make the menus
		"""
		menubar = self.menuBar()
		
		# File menu
		filemenu = menubar.addMenu('&File')
		# exit
		exit = QtGui.QAction('Exit', self)
		exit.setShortcut('Ctrl+Q')
		exit.setStatusTip('Exit application')
		self.connect(exit, QtCore.SIGNAL('triggered()'), QtCore.SLOT('close()'))
		filemenu.addAction(exit)
		
		# Project
		projectmenu = menubar.addMenu('&Project')
		newproject = QtGui.QAction('New Project', self)
		newproject.setShortcut('Ctrl+N')
		newproject.setStatusTip('Create New Project')
		self.connect(newproject, QtCore.SIGNAL('triggered()'), self._on_newproject)
		projectmenu.addAction(newproject)
		newproject.setStatusTip('New Project')
		openproject = QtGui.QAction('Open Project', self)
		openproject.setShortcut('Ctrl+O')
		openproject.setStatusTip('Open Project')
		self.connect(openproject, QtCore.SIGNAL('triggered()'), self._on_openproject)
		projectmenu.addAction(openproject)
		editproject = QtGui.QAction('Edit Project', self)
		editproject.setShortcut('Ctrl+E')
		editproject.setStatusTip('Edit Project')
		projectmenu.addAction(editproject)
		importdataproject = QtGui.QAction('Import data', self)
		importdataproject.setShortcut('Ctrl+I')
		importdataproject.setStatusTip('Import data from EMEN or Disk')
		projectmenu.addAction(importdataproject)
		
		# Options
		#optionsmenu = menubar.addMenu('&Options')
		
		# Mode
		modemenu = menubar.addMenu('&Mode')
		# modemenu
		sprmode = QtGui.QAction('SPR mode', self)
		sprmode.setShortcut('Ctrl+S')
		sprmode.setStatusTip('SPR mode')
		self.connect(sprmode, QtCore.SIGNAL('triggered()'), self._on_sprmode)
		modemenu.addAction(sprmode)
		tomomode = QtGui.QAction('TOMO mode', self)
		tomomode.setShortcut('Ctrl+T')
		tomomode.setStatusTip('TOMO mode')
		self.connect(tomomode, QtCore.SIGNAL('triggered()'), self._on_tomomode)
		modemenu.addAction(tomomode)
		
		# Utils
		utilsmenu = menubar.addMenu('&Utilities')
		filebrowser = QtGui.QAction('File Browser', self)
		filebrowser.setShortcut('Ctrl+F')
		filebrowser.setStatusTip('File Browser')
		configureEMEN = QtGui.QAction('Configure EMEN', self)
		configureEMEN.setShortcut('Ctrl+C')
		configureEMEN.setStatusTip('Configure EMEN')
		utilsmenu.addAction(filebrowser)
		utilsmenu.addAction(configureEMEN)
	
		# Help
		helpmenu = menubar.addMenu('&Help')
		about = QtGui.QAction('About', self)
		about.setStatusTip('About')
		helpmenu.addAction(about)
		helpdoc = QtGui.QAction('Help', self)
		helpdoc.setStatusTip('Help')
		helpmenu.addAction(helpdoc)
	
	def _on_sprmode(self):
		self.tree_stacked_widget.setCurrentIndex(0)
		
	def _on_tomomode(self):
		self.tree_stacked_widget.setCurrentIndex(1)
		
	def _on_newproject(self):
		np = DialogNewProject()
		np.exec_()
		self.activateWindow()
		
	def _on_openproject(self):
		np = DialogOpenProject(self)
		np.exec_()
		self.activateWindow()
		
	def makeTilteBarWidget(self):
		"""
		Make the title bar widget (ICON + label)
		"""
		tbwidget = QtGui.QFrame()
		tbwidget.setFrameShape(QtGui.QFrame.StyledPanel)
		grid = QtGui.QGridLayout()
		self.PMIcon = PMIcon(self.pm_icon, tbwidget)
		grid.addWidget(self.PMIcon,0 , 0, 2, 1)
		self.PMTitle = QtGui.QLabel("EMAN2 Project Manager ")
		self.PMProjectNameBanner = QtGui.QLabel("Project Name: "+self.pn_project_name)
		self.PMProjectNameBanner.setAlignment(QtCore.Qt.AlignCenter)
		titlefont = QtGui.QFont()
		titlefont.setPointSize(30)
		titlefont.setBold(True)
		self.PMTitle.setFont(titlefont)
		self.PMTitle.setAlignment(QtCore.Qt.AlignBottom)
		pmnamefont = QtGui.QFont()
		pmnamefont.setBold(True)
		self.PMProjectNameBanner.setFont(pmnamefont)
		grid.addWidget(self.PMTitle, 0, 1)
		grid.addWidget(self.PMProjectNameBanner, 1, 1)
		tbwidget.setLayout(grid)
		
		return tbwidget
	
	def _load_icons(self):
		self.icons = {}
		self.icons["single_image"] = QtGui.QIcon("../images/macimages/single_image.png")
		self.icons["multiple_images"] = QtGui.QIcon("../images/macimages/multiple_images.png")
		self.icons["green_boxes"] = QtGui.QIcon("../images/macimages/green_boxes.png")
		self.icons["ctf"] = QtGui.QIcon("../images/macimages/ctf.png")
		self.icons["web"] = QtGui.QIcon("../images/macimages/classes.png")
		self.icons["single_image_3d"] = QtGui.QIcon("../images/macimages/single_image_3d.png")
		self.icons["refine"] = QtGui.QIcon("../images/macimages/refine.png")
		self.icons["eulers"] = QtGui.QIcon("../images/macimages/eulerxplor.png")
		self.icons["resolution"] = QtGui.QIcon("../images/macimages/plot.png")
		self.icons["tomo_hunter"] = QtGui.QIcon("../images/macimages/tomo_hunter.png")
	
	def makeStackedWidget(self):
		"""
		This is the stacked widget to manage the tree types
		"""
		self.tree_stacked_widget = QtGui.QStackedWidget()
		self.tree_stacked_widget.setMaximumWidth(300)
		self.tree_stacked_widget.addWidget(self.makeSPRTreeWidget())
		self.tree_stacked_widget.addWidget(self.makeTomoTreeWidget())
		
		return self.tree_stacked_widget
	
	def makeStackedGUIwidget(self):
		"""
		Make a stacked widget
		When a python script is called for the first time a GUI widget is made and added to the stack
		"""
		self.gui_stacked_widget = QtGui.QStackedWidget()
		self.gui_stacked_widget.setFrameShape(QtGui.QFrame.StyledPanel)
		
		return self.gui_stacked_widget
	
	def _add_children(self, toplevel, widgetitem, treewidgetdict):
		""" recursive hlper function for loadTree"""
		for child in toplevel["CHILDREN"]:
			treewidgetdict[child["NAME"]] = QtGui.QTreeWidgetItem(QtCore.QStringList(child["NAME"]))
			treewidgetdict[child["NAME"]].setIcon(0, self.icons[child["ICON"]])
			self._add_children(child, treewidgetdict[child["NAME"]], treewidgetdict)
			widgetitem.addChild(treewidgetdict[child["NAME"]])
			
		
	def loadTree(self, filename, treename, treewidgetdict):
		"""
		Load a workflow tree from a JSON file
		@param filename The JSON filename
		@param treename The name of the QTreWidget
		@param treewidgetdict the dictionary containing the QtreeWidgets
		"""
		jsonfile = open(filename, 'r')
		tree = json.load(jsonfile)
		jsonfile.close()
		
		QTree = QtGui.QTreeWidget()
		QTree.setMinimumHeight(400)
		QTree.setHeaderLabel(treename)
		
		for toplevel in tree:
			treewidgetdict[toplevel["NAME"]] = QtGui.QTreeWidgetItem(QtCore.QStringList(toplevel["NAME"]))
			treewidgetdict[toplevel["NAME"]].setIcon(0, self.icons[toplevel["ICON"]])
			self._add_children(toplevel, treewidgetdict[toplevel["NAME"]], treewidgetdict)
			QTree.addTopLevelItem(treewidgetdict[toplevel["NAME"]])
			
		return QTree
		
	def makeSPRTreeWidget(self):
		"""
		Make the SPR tree
		"""
		self.SPRtreewidgetdict = {}
		self.SPRtree = self.loadTree('spr.json', 'SPR', self.SPRtreewidgetdict)
		return self.SPRtree
	
	def makeTomoTreeWidget(self):
		"""
		Make the tomo tree widget
		"""
		self.TOMOtreewidgetdict = {}
		self.TOMOtree = self.loadTree('tomo.json', 'Tomography', self.TOMOtreewidgetdict)
		return self.TOMOtree
			
class EMAN2StatusBar(QtGui.QLabel):
	"""
	The Stats bar for PM
	"""
	def __init__(self, text):
		QtGui.QLabel.__init__(self, text)
		self.setFrameShape(QtGui.QFrame.Panel | QtGui.QFrame.Sunken)
		self.setLineWidth(2)
		self.setMargin(4)
		
	def setMessage(self, text):
		self.setText(text)

class PMIcon(QtGui.QLabel):
	"""
	The Icon manager for PM
	"""
	def __init__(self, image, parent=None):
		QtGui.QLabel.__init__(self, ("<img src=\"%s\" />")%image, parent)
	
	def setIcon(self, image):
		self.setText(("<img src=\"%s\" />")%image)
		
class LogBook(QtGui.QWidget):
	"""
	The Logbook for PM
	"""
	def __init__(self):
		QtGui.QWidget.__init__(self)
		self.setWindowTitle('LogBook')
		grid = QtGui.QGridLayout()
		font = QtGui.QFont()
		font.setBold(True)
		textlabel = QtGui.QLabel("EMAN2 LogBook")
		textlabel.setFont(font)
		self.texteditbox = QtGui.QTextEdit()
		grid.addWidget(textlabel,0,0)
		grid.addWidget(self.texteditbox,1,0,1,2)
		self.savepb = QtGui.QPushButton("Save")
		self.closepb = QtGui.QPushButton("Close")
		grid.addWidget(self.savepb, 2,0)
		grid.addWidget(self.closepb, 2,1)
		self.setLayout(grid)

class DialogNewProject(QtGui.QDialog):
	"""
	Generate the new projects dialog
	"""
	def __init__(self):
		QtGui.QDialog.__init__(self, pm)
		self.setWindowTitle('New Project')
		self.pm = pm
		
		# Dialog line edit fields
		minwidth = 8*max(len(self.pm.pm_cwd),len(self.pm.pm_icon))
		frame = QtGui.QFrame()
		frame.setFrameStyle(QtGui.QFrame.StyledPanel)
		grid = QtGui.QGridLayout()
		project_name_label = QtGui.QLabel("Project Name")
		self.project_name = QtGui.QLineEdit("New Project")
		self.project_name.setMinimumWidth(minwidth)
		grid.addWidget(project_name_label, 0, 0)
		grid.addWidget(self.project_name, 0, 1)
		project_directory_label = QtGui.QLabel("Project Directory")
		self.project_directory = QtGui.QLineEdit(self.pm.pm_cwd)
		self.project_directory.setMinimumWidth(minwidth)
		grid.addWidget(project_directory_label, 1, 0)
		grid.addWidget(self.project_directory, 1, 1)
		icon_path_label = QtGui.QLabel("Project Icon")
		self.icon_path = QtGui.QLineEdit(self.pm.pm_icon)
		self.icon_path.setMinimumWidth(minwidth)
		grid.addWidget(icon_path_label, 2, 0)
		grid.addWidget(self.icon_path, 2, 1)
		# Microscope parameters
		if self.pm.usingEMEN == False:
			micrscope_cs_label = QtGui.QLabel("Microscope CS")
			self.micrscope_cs = QtGui.QLineEdit()
			microscope_voltage_label = QtGui.QLabel("Microscope Voltage")
			self.microscope_voltage = QtGui.QLineEdit()
			microscope_apix_label = QtGui.QLabel("Microscope apix")
			self.microscope_apix = QtGui.QLineEdit()
			grid.addWidget(micrscope_cs_label, 3, 0)
			grid.addWidget(self.micrscope_cs, 3, 1)
			grid.addWidget(microscope_voltage_label, 4, 0)
			grid.addWidget(self.microscope_voltage, 4, 1)
			grid.addWidget(microscope_apix_label, 5, 0)
			grid.addWidget(self.microscope_apix, 5, 1)
		else:
			pass
			# Get these data from EMEN
			
		frame.setLayout(grid)
		
		# Ok, cancel buttons
		makeproject_pb = QtGui.QPushButton("Make Project")
		cancel_pb = QtGui.QPushButton("Cancel")
		sgrid = QtGui.QGridLayout()
		sgrid.addWidget(frame,0,0,1,2)
		sgrid.addWidget(makeproject_pb,1,0)
		sgrid.addWidget(cancel_pb,1,1)
		self.setLayout(sgrid)
		
		self.connect(makeproject_pb, QtCore.SIGNAL('clicked()'), self._on_makeproject)
		self.connect(cancel_pb, QtCore.SIGNAL('clicked()'), self._on_cancel)
		
	def _on_makeproject(self):
		for project in self.pm.pm_projects_db.keys():
			if project == self.project_name.text():
				self.pm.statusbar.setMessage("Project Name is alreay in use!!!")
				return
			if self.pm.pm_projects_db[project]["CWD"] == self.project_directory.text():
				self.pm.statusbar.setMessage("This project directory is already in use by another project!!!")
				return
		self.pm.pm_projects_db[self.project_name.text()] = {"CWD":self.project_directory.text(),"ICON":self.icon_path.text(),"CS":self.micrscope_cs.text(),"VOLTAGE":self.microscope_voltage.text(),"APIX":self.microscope_apix.text()}
		self.pm.statusbar.setMessage("Project: "+self.project_name.text()+" created :)")
		self.done(0)
		
	def _on_cancel(self):
		self.done(1)
		
class DialogOpenProject(QtGui.QDialog):
	"""
	Generate the new projects dialog
	"""
	def __init__(self, pm):
		QtGui.QDialog.__init__(self)
		self.pm = pm
		self.setWindowTitle('Open Project')
		grid = QtGui.QGridLayout()
		self.list_widget = QtGui.QListWidget()
		# Populate list
		for project in self.pm.pm_projects_db.keys():
			item = QtGui.QListWidgetItem(project)
			self.list_widget.addItem(item) 
			if self.pm.pn_project_name == project:
				self.list_widget.setCurrentItem(item)
		self.list_widget.sortItems()
		
		grid.addWidget(self.list_widget, 0, 0, 1, 2)
		open_pb = QtGui.QPushButton("Open Project")
		cancel_pb = QtGui.QPushButton("Cancel")
		grid.addWidget(open_pb)
		grid.addWidget(cancel_pb)
		self.setLayout(grid)
		
		self.connect(open_pb, QtCore.SIGNAL('clicked()'), self._on_openproject)
		self.connect(cancel_pb, QtCore.SIGNAL('clicked()'), self._on_cancel)
		
	def _on_openproject(self):
		self.pm.pn_project_name = self.list_widget.currentItem().text()
		self.pm.statusbar.setMessage("Opened Project: "+self.pm.pn_project_name)
		self.pm.PMProjectNameBanner.setText("Project Name: "+self.pm.pn_project_name)
		self.PMIcon.setIcon(self.pm.pm_projects_db[self.pm.pn_project_name]['ICON'])
		self.done(0)
		# do some stuff......
		
	def _on_cancel(self):
		self.done(1)
		
if __name__ == "__main__":
	import sys
	app = QtGui.QApplication(sys.argv)
	pm = EMProjectManager()
	pm.show()
	app.exec_()