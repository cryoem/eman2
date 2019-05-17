#!/usr/bin/env python
from __future__ import print_function
from __future__ import division
#
# Author: John Flanagan Oct 20th 2011 (jfflanag@bcm.edu)
# Copyright (c) 2000-2011 Baylor College of Medicine

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

from past.utils import old_div
from builtins import range
from EMAN2 import *
from PyQt5 import QtCore, QtGui, QtWidgets
from eman2_gui.pmicons import *
import os, json, re, glob, signal
import subprocess
from eman2_gui.empmwidgets import *
from eman2_gui.valslider import EMQTColorWidget
from eman2_gui.embrowser import EMBrowserWidget

class EMProjectManager(QtWidgets.QMainWindow):
	""" The EM Project Manager is a QT application to provide a GUI for EMAN2 job management.
	See the wiki for more details. For documentation see the EMAN2 WIKI """
	def __init__(self):
		QtWidgets.QMainWindow.__init__(self)
		# Default PM attributes
		self.pm_cwd = os.getcwd()
		self.pn_project_name_default='Unknown'
		self.pm_icon = self.pm_icon_default = get_image_directory() + "EMAN2Icon.png"
		self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(pmicon)))

		# Set defaults
		self.notebook = None
		self.taskmanager = None
		self.wikipage = None
		self.thehelp = None
		self.table = None

		# Load the project DataBase
		self.loadPMdb()
		# Load icons
		self._load_icons()

		# Make the project manager layout
		self.makeMenues()
		font = QtGui.QFont()
		font.setBold(True)
		centralwidget = QtWidgets.QWidget()
		vsplitter = QtWidgets.QSplitter(QtCore.Qt.Vertical)

		# Make the titlebars
		grid = QtWidgets.QGridLayout()
		grid.addWidget(self.makeTilteBarWidget(), 0, 0, 1, 3)
		#grid.addWidget(workflowcontrollabel, 1,0)
		grid.addWidget(self.makeModeWidget(font))
		guilabel = QtWidgets.QLabel("EMAN2 Program Interface", centralwidget)
		guilabel.setFont(font)
		guilabel.setMaximumHeight(20)
		grid.addWidget(guilabel, 1,1)

		# Make the GUI, Tree and ToolBar
		grid.addWidget(self.makeStackedWidget()   ,2,0,2,1)
		grid.addWidget(self.makeStackedGUIwidget(),2,1,1,1)
		grid.addWidget(self.makeGUIToolButtons()  ,2,2,1,1)
		grid.addWidget(self.makeCMDButtonsWidget(),3,1,1,1)

		grid.setColumnStretch(0,0)
		grid.setColumnStretch(1,1)
		grid.setColumnStretch(2,0)

		# Make status bar
		self.statusbar = EMAN2StatusBar("Welcome to the EMAN2 Project Manager","font-weight:bold;")
		centralwidget.setLayout(grid)
		# Add splitter to adjust status bar
		vsplitter.addWidget(centralwidget)
		vsplitter.addWidget(self.statusbar)
		vsplitter.setSizes([1000,100])

		self.setCentralWidget(vsplitter)
		E2loadprojtype("e2projectmanager","main",self)
		E2loadappwin("e2projectmanager","main",self)
		# Update the project are construction
		self.updateProject()

	def mousePressEvent(self, event):
		self.setFocus()
	
	def closeEvent(self, event):
		""" Upon PM close, close the taskmanager and the logbook """
		E2saveappwin("e2projectmanager","main",self)
		E2saveprojtype("e2projectmanager","main",self)
		if self.notebook: self.notebook.close()
		if self.taskmanager: self.taskmanager.close()
		if self.thehelp: self.thehelp.close()

	def loadPMdb(self):
		"""
		Load the PM database. This is a global database. Each user on each machine has one
		"""
		self.pm_projects_db = js_open_dict("info/project.json")
		self.pn_project_name = self.pm_projects_db.setdefault("project_name",dfl=self.pn_project_name_default)

	def _load_icons(self):
		"""
		Load icons used for the tree. Additional icons can be added using icons.json
		"""
		self.icons = {}
		EMAN2DIR = e2getinstalldir()

		jsonfile = open(EMAN2DIR+'/lib/pmconfig/icons.json', 'r')
		data = jsonfile.read()
		data = self.json_strip_comments(data)
		tree = json.loads(data)
		jsonfile.close()
		icons=tree[0]

		for icon in icons:
			self.icons[icon] = QtGui.QIcon(EMAN2DIR+icons[icon])

	def makeMenues(self):
		"""
		Make the menus
		"""
		menubar = self.menuBar()

		# File menu
		filemenu = menubar.addMenu('&File')
		# exit
		exit = QtWidgets.QAction('Exit', self)
		exit.setShortcut('Ctrl+Q')
		exit.setStatusTip('Exit application')
		exit.triggered.connect(self.close)
		filemenu.addAction(exit)

		# Project
		projectmenu = menubar.addMenu('&Project')
		openproject = QtWidgets.QAction('Open Project', self)
		openproject.setShortcut('Ctrl+O')
		openproject.setStatusTip('Open Project')
		openproject.triggered.connect(self._on_openproject)
		projectmenu.addAction(openproject)
		editproject = QtWidgets.QAction('Edit Project', self)
		editproject.setShortcut('Ctrl+E')
		editproject.setStatusTip('Edit Project')
		editproject.triggered.connect(self._on_editproject)
		projectmenu.addAction(editproject)
		# Options
		# optionsmenu = menubar.addMenu('&Options')

		# Utils
		utilsmenu = menubar.addMenu('&Utilities')
		filebrowser = QtWidgets.QAction('File Browser', self)
		filebrowser.setShortcut('Ctrl+F')
		filebrowser.setStatusTip('File Browser')
		utilsmenu.addAction(filebrowser)
		utilsmenu.addSeparator()
		filebrowser.triggered.connect(self._on_browse)
		self.dumpterminal = QtWidgets.QAction('Dump Terminal', self)
		self.dumpterminal.setCheckable(True)
		self.dumpterminal.setChecked(False)
		utilsmenu.addAction(self.dumpterminal)



		# Help
		helpmenu = menubar.addMenu('&Help')
		about = QtWidgets.QAction('About', self)
		about.setStatusTip('About')
		helpmenu.addAction(about)
		helpdoc = QtWidgets.QAction('Help', self)
		helpdoc.setStatusTip('Help')
		helpmenu.addAction(helpdoc)

	def makeModeWidget(self, font):
		""" Return a mode control widget """
		widget = QtWidgets.QWidget()
		box = QtWidgets.QHBoxLayout()
		box.setContentsMargins(0,0,0,0)
		workflowcontrollabel = QtWidgets.QLabel("Workflow Mode", widget)
		workflowcontrollabel.setFont(font)
		workflowcontrollabel.setMaximumHeight(20)
		self.modeCB = QtWidgets.QComboBox()
		# To add a new mode add an item to the list, and then add the json file in fuction: makeStackedWidget
		self.modeCB.addItem("SPR")
		self.modeCB.addItem("Tomo")
		self.modeCB.addItem("SPT (Legacy)")

		box.addWidget(workflowcontrollabel)
		box.addWidget(self.modeCB)
		widget.setLayout(box)

		self.modeCB.activated[int].connect(self._onModeChange)

		return widget

	def _onModeChange(self, idx):
		self.tree_stacked_widget.setCurrentIndex(idx)
		if idx == 1 or idx == 2: self.pm_icon = get_image_directory() + "tomoseg.png"
		else: self.pm_icon = get_image_directory() + "EMAN2Icon.png"
		
		self.pm_projects_db["project_icon"] = self.pm_icon

		self.updateProject()
		self._on_cmd_cancel()

	def _on_openproject(self):
		self.openbrowser = EMBrowserWidget(withmodal=True,multiselect=False)
		self.openbrowser.ok.connect(self._onopen_ok)
		self.openbrowser.cancel.connect(self._onopen_cancel)
		self.openbrowser.show()
		self.activateWindow()

	def _onopen_cancel(self):
		self.openbrowser.close()

	def _onopen_ok(self):
		# Here we need to get a new dir and change to it.

		self.changeDirectory(self.openbrowser.getCWD())
		self.loadPMdb()
		self.updateProject()
		self.statusbar.setMessage("Project Name is Now: '%s'"%self.pn_project_name, "font-weight:bold;")

	def _on_editproject(self):
		""" Open edit dialog """
		np = ProjectDialog(self)
		np.exec_()
		self.activateWindow()

	def makeTilteBarWidget(self):
		"""
		Make the title bar widget (Project ICON + label)
		"""
		tbwidget = QtWidgets.QFrame()
		tbwidget.setFrameShape(QtWidgets.QFrame.StyledPanel)
		grid = QtWidgets.QGridLayout()
		self.PMIcon = PMIcon(self.pm_icon, tbwidget)
		self.PMIcon.setAlignment(QtCore.Qt.AlignLeft)
		grid.addWidget(self.PMIcon,0 , 0, 2, 1)
		self.PMTitle = QtWidgets.QLabel("EMAN2 Project Manager ")
		self.PMTitle.setAlignment(QtCore.Qt.AlignCenter)
		self.PMProjectNameBanner = QtWidgets.QLabel("Project Name: "+self.pn_project_name)
		self.PMProjectNameBanner.setAlignment(QtCore.Qt.AlignCenter)
		titlefont = QtGui.QFont()
		titlefont.setPointSize(30)
		titlefont.setBold(True)
		titlefont.setItalic(True)
		self.PMTitle.setFont(titlefont)
		self.PMTitle.setAlignment(QtCore.Qt.AlignBottom)
		pmnamefont = QtGui.QFont()
		pmnamefont.setBold(True)
		self.PMProjectNameBanner.setFont(pmnamefont)
		grid.addWidget(self.PMTitle, 0, 1, 1, 2)
		grid.addWidget(self.PMProjectNameBanner, 1, 1, 1, 2)
		tbwidget.setLayout(grid)

		return tbwidget

	def makeStackedWidget(self):
		"""
		This is the stacked widget to manage the tree types. To add modes, do so here. Be sure to add the mode to the combo box in function: makeModeWidget
		"""
		self.tree_stacked_widget = QtWidgets.QStackedWidget()
		self.tree_stacked_widget.setMinimumWidth(300)
		self.tree_stacked_widget.addWidget(self.makeTreeWidget(e2getinstalldir()+'/lib/pmconfig/spr.json', 'Single Particle Refinement'))
		self.tree_stacked_widget.addWidget(self.makeTreeWidget(e2getinstalldir()+'/lib/pmconfig/tomo.json', 'Tomography'))
		self.tree_stacked_widget.addWidget(self.makeTreeWidget(e2getinstalldir()+'/lib/pmconfig/tomo_legacy.json', 'SPT (Legacy)'))

		return self.tree_stacked_widget

	def makeStackedGUIwidget(self):
		"""
		Make a stacked widget
		When a python script is called for the first time a GUI widget is made and added to the stack
		The first widget on the stack is the blank widget, the rest are e2program widgets
		"""
		self.gui_stacked_widget = QtWidgets.QStackedWidget()
		# Set the initial height of the browser
		#self.gui_stacked_widget.setMinimumHeight(250)
		self.gui_stacked_widget.setFrameShape(QtWidgets.QFrame.StyledPanel)
		# Blank screen widget
		self.gui_stacked_widget.addWidget(QtWidgets.QWidget())
		self.stackedWidgetHash = {}

		return self.gui_stacked_widget

	def makeGUIToolButtons(self):
		"""
		Get ToolBar widget
		"""
		toolwidget = QtWidgets.QFrame()
		# toolwidget.setFrameShape(QtWidgets.QFrame.StyledPanel)
		tbox = QtWidgets.QVBoxLayout()
		self.browsebutton = QtWidgets.QToolButton()
		self.browsebutton.setIcon(QtGui.QIcon(QtGui.QPixmap(browseicon)))
		self.browsebutton.setToolTip("Browse button")
		self.browsebutton.setMinimumWidth(30)
		self.browsebutton.setMinimumHeight(30)
		self.helpbutton = PMToolButton()
		self.helpbutton.setIcon(QtGui.QIcon(QtGui.QPixmap(helpicon)))
		self.helpbutton.setToolTip("Help button")
		self.logbutton = PMToolButton()
		self.logbutton.setToolTip("Display the note book")
		self.logbutton.setIcon(QtGui.QIcon(QtGui.QPixmap(noteicon)))
		self.taskmanagerbutton = PMToolButton()
		self.taskmanagerbutton.setToolTip("Display the task manager")
		self.taskmanagerbutton.setIcon(QtGui.QIcon(QtGui.QPixmap(taskicon)))
		tbox.addWidget(self.browsebutton)
		tbox.addWidget(self.helpbutton)
		tbox.addWidget(self.logbutton)
		tbox.addWidget(self.taskmanagerbutton)
		self.makeProgramToolButtons(tbox)
		tbox.setContentsMargins(0,0,0,0)
		tbox.setAlignment(QtCore.Qt.AlignTop)
		toolwidget.setLayout(tbox)

		self.browsebutton.clicked.connect(self._on_browse)
		self.helpbutton.stateChanged[bool].connect(self._on_helpbutton)
		self.logbutton.stateChanged[bool].connect(self._on_logbutton)
		self.taskmanagerbutton.stateChanged[bool].connect(self._on_taskmgrbutton)

		return toolwidget

	def _on_browse(self):
		"""
		Launch the browser
		"""
		win=EMBrowserWidget(withmodal=False,multiselect=False)
		win.show()
		win.setAttribute(QtCore.Qt.WA_DeleteOnClose)
	#self.window = EMBrowserWidget(withmodal=False,multiselect=False)
	#self.window.show()

	def _on_helpbutton(self, state):
		"""
		Load the help info
		"""
		if state:
			self.loadTheHelp()
		else:
			self.thehelp.hide()

	def _on_logbutton(self, state):
		"""
		Load the log book
		"""
		if state:
			self.loadNoteBook()
		else:
			self.notebook.hide()

	def _on_wikibutton(self):
		"""
		Load the wiki help
		"""
		self.loadWiki()

	def _on_wizardbutton(self):
		"""
		Load the wizardbutton
		"""
		self.loadWizard()

	def _on_taskmgrbutton(self, state):
		"""
		Load the log book
		"""
		if state:
			self.loadTaskManager()
		else:
			self.taskmanager.hide()

	def loadTheHelp(self):
		"""
		Make the help
		"""
		if not self.thehelp:
			self.thehelp = TheHelp(self)
		self.thehelp.show()

	def loadNoteBook(self):
		"""
		Make logbook
		"""
		if not self.notebook:
			self.notebook = NoteBook(self)
		self.notebook.show()

	def loadTaskManager(self):
		"""
		Make task manager
		"""
		if not self.taskmanager:
			self.taskmanager = TaskManager(self)
		self.taskmanager.show()

	def makeProgramToolButtons(self, tbox):
		"""
		Make the browse button
		"""
		#programtoolwidget = QtWidgets.QFrame()
		#programtoolwidget.setFrameShape(QtWidgets.QFrame.StyledPanel)
		#tbox = QtWidgets.QVBoxLayout()
		self.wikibutton = QtWidgets.QToolButton()
		self.wikibutton.setIcon(QtGui.QIcon(QtGui.QPixmap(wikiicon)))
		self.wikibutton.setToolTip("Show Wiki button")
		self.wikibutton.setMinimumWidth(30)
		self.wikibutton.setMinimumHeight(30)
		self.wikibutton.setEnabled(False)
		self.expertbutton = PMToolButton()
		self.expertbutton.setDown(False, quiet=True)
		self.expertbutton.setIcon(QtGui.QIcon(QtGui.QPixmap(experticon)))
		self.expertbutton.setToolTip("ExpertMode")
		self.expertbutton.setEnabled(False)
		self.wizardbutton = QtWidgets.QToolButton()
		self.wizardbutton.setIcon(QtGui.QIcon(QtGui.QPixmap(wizardicon)))
		self.wizardbutton.setToolTip("Form Wizard")
		self.wizardbutton.setMinimumWidth(30)
		self.wizardbutton.setMinimumHeight(30)
		self.wizardbutton.setEnabled(False)
		tbox.addWidget(self.wikibutton)
		tbox.addWidget(self.wizardbutton)
		tbox.addWidget(self.expertbutton)
		#programtoolwidget.setLayout(tbox)
		#tbox.setContentsMargins(0,0,0,0)
		#tbox.setAlignment(QtCore.Qt.AlignTop)

		self.wikibutton.clicked.connect(self._on_wikibutton)
		self.wizardbutton.clicked.connect(self._on_wizardbutton)
		self.expertbutton.stateChanged[bool].connect(self._on_expertmodechanged)

	#return programtoolwidget

	def _on_expertmodechanged(self, state):
		"""
		Change the GUI upon expert mode toggling
		"""
		if self.gui_stacked_widget.currentIndex() >= 1:	# First widget in the stack is the blank widget
			self.setProgramExpertMode(int(state)+1)	# If we are able to set the state, then button is available (state = 0), hence we toggle between '1', available and OFF or '2', available and ON
			self._set_GUI(self.getProgram(), self.getProgramMode())

	def loadWiki(self):
		"""
		Make wiki
		"""
		import webbrowser
		webbrowser.open(self.getProgramWikiPage())

	def loadWizard(self):
		"""
		Load the wizard from a filebox
		"""

		jsonfile = open(e2getinstalldir()+self.getProgramWizardFile(), 'r')
		data = jsonfile.read()
		data = self.json_strip_comments(data)
		wizarddata = json.loads(data)
		jsonfile.close()

		# Now load the wizard (wizard is modal)
		self.wizard = EMWizard(wizarddata, self.gui_stacked_widget.currentWidget().guiwidget)
		self.wizard.show()

	def makeCMDButtonsWidget(self):
		"""
		Get the e2program interface command buttons widget
		"""
		cmdwidget = QtWidgets.QFrame()
		cmdwidget.setFrameShape(QtWidgets.QFrame.StyledPanel)
		cmdwidget.setMaximumHeight(40)
		hbox = QtWidgets.QHBoxLayout()
		self.cancelbutton = QtWidgets.QPushButton("Cancel")
		self.launchbutton = QtWidgets.QPushButton("Launch")
		self.launchbutton.setFocusPolicy(Qt.StrongFocus)
		hbox.addWidget(self.cancelbutton)
		hbox.addWidget(self.launchbutton)
		hbox.setContentsMargins(4,4,4,4)
		cmdwidget.setLayout(hbox)

		self.cancelbutton.clicked.connect(self._on_cmd_cancel)
		self.launchbutton.clicked.connect(self._on_cmd_launch)

		return cmdwidget

	def clearE2Interface(self):
		"""
		Clear the e2program interface
		"""
		self.getCurrentTree().clearSelection()
		self.tree_stacked_widget.currentWidget().setCurrentItem(None)
		self.gui_stacked_widget.setCurrentIndex(0)	# blank screen
		self.updateProject()

	def _on_cmd_cancel(self):
		"""
		Cancel the present form. This just pulls up the blank screen and updates the project
		"""
		self.clearE2Interface()

	def _on_cmd_launch(self):
		"""
		Launch the command from the GUI
		"""
		if self.gui_stacked_widget.currentIndex() == 0: return	# No cmd to run

		# Else run the command
		cmd = self.gui_stacked_widget.currentWidget().getCommand()
		if self.launchScript(cmd):
			self.clearE2Interface()

	def launchScript(self, cmd):
		"""
		Start the script running
		"""
		if not cmd: return False	# Don't excecute a broken script
		# --ipd=-2 ; tell .eman2log.txt that this job is already in the pm note book
		if self.dumpterminal.isChecked():
			# This wont work on Windows(or for multiple e2projectmanger instances)....... Fix later
			# Also this code is experimental, it needs to be finished....
			try : os.unlink("%s/pmfifo"%os.getcwd())
			except : pass
			os.mkfifo("%s/pmfifo"%os.getcwd())
			stdoutpipe = open("%s/pmfifo"%os.getcwd(),"w+",0)
		else:
			stdoutpipe = None

		if self.notebook and self.getProgramNoteLevel() > 0:
			# Only take notes if note level is greater than 0
			self.notebook.insertNewJob(cmd,local_datetime())
			self.notebook.writeNotes()
			child = EMPopen((str(cmd)+" --ppid=-2"), shell=True, cwd=self.pm_cwd, stdout=stdoutpipe, bufsize=1)
		else:
			if self.getProgramNoteLevel() > 0:
				print("NOT Writing notes, ppid=-1")
				child = EMPopen(str(cmd), shell=True, cwd=self.pm_cwd, stdout=stdoutpipe, bufsize=1)
			else:
				print("NOT Writing notes, ppid=-2")
				child = EMPopen((str(cmd)+" --ppid=-2"), shell=True, cwd=self.pm_cwd, stdout=stdoutpipe, bufsize=1)
		# Dump terminal stdout if desired
		if self.dumpterminal.isChecked(): child.realTimeCommunicate(self.statusbar)

		self.statusbar.setMessage("Program %s Launched!!!!"%str(cmd).split()[0],"font-weight:bold;")

		return True

	def loadUsage(self, program):
		"""
		Read in and return the usage from an e2program
		"""
		try:
			f = open(e2getinstalldir()+"/bin/"+program,"r")
		except:
			self.statusbar.setMessage("Can't open usage file '%s'"%program,"color:red;")
			return
		begin = False
		helpstr = ""
		bregex = re.compile('usage\s*=\s*"""')
		eregex = re.compile('"""')
		for line in f:
			if re.search(eregex, line) and begin:
				line = re.sub(eregex, "", line)
				helpstr = helpstr + line.strip()
				break
			if re.search(bregex, line):
				begin = True
				line = re.sub(bregex, "", line)
			if begin:
				helpstr = helpstr + line.strip() + " "
		f.close()

		return helpstr

	def _add_children(self, toplevel, widgetitem):
		"""
		Recursive helper function for loadTree
		"""
		for child in toplevel["CHILDREN"]:
			qtreewidget = PMQTreeWidgetItem([child["NAME"]])
			if "ICON" in  child: qtreewidget.setIcon(0, self.icons[child["ICON"]])
			# Optional program
			if "PROGRAM" in child: qtreewidget.setProgram(child["PROGRAM"])
			# Optional table to display rather than a program
			if "TABLE" in child: qtreewidget.setTable(child["TABLE"])
			# Optional mode for the program to run in. The default is to have no mode
			if "MODE" in child: qtreewidget.setMode(child["MODE"])
			# Option note level and note elevel > 0 means add job to notebook. Good to prevent a lot of crap from piling up!
			if "NOTELEVEL" in child: qtreewidget.setNoteLevel(child["NOTELEVEL"])
			# Optional wikipage, to load info about the program from the wiki
			if "WIKIPAGE" in child: qtreewidget.setWikiPage(child["WIKIPAGE"])
			# Optional wizard file to add a wizard for this e2program
			if "WIZARD" in child: qtreewidget.setWizardFile(child["WIZARD"])
			# Optional expertmode, to enable expert mode GUI widgets
			if "EXPERT" in child:  qtreewidget.setExpertMode(child["EXPERT"])
			self._add_children(child, qtreewidget)
			widgetitem.addChild(qtreewidget)


	def loadTree(self, filename, treename):
		"""
		Load a workflow tree from a JSON file
		@param filename The JSON filename
		@param treename The name of the QTreeWidget
		@param treewidgetdict the dictionary containing the QTreeWidgets
		"""
		try:
			jsonfile = open(filename, 'r')
		except:
			self.statusbar.setMessage("Can't open configuration file '%s'"%filename,"color:red;")
			return
		data = jsonfile.read()
		data = self.json_strip_comments(data)
		tree = json.loads(data)
		jsonfile.close()

		QTree = QtWidgets.QTreeWidget()
		QTree.setHeaderLabel(treename)

		for toplevel in tree:
			qtreewidget = PMQTreeWidgetItem([toplevel["NAME"]])
			if "ICON" in toplevel: qtreewidget.setIcon(0, self.icons[toplevel["ICON"]])
			# Optional program
			if "PROGRAM" in toplevel: qtreewidget.setProgram(toplevel["PROGRAM"])
			# Optional table to display rather than a program
			if "TABLE" in toplevel: qtreewidget.setTable(toplevel["TABLE"])
			# Optional mode for the program to run in. The default is to have no mode
			if "MODE" in toplevel: qtreewidget.setMode(toplevel["MODE"])
			# Option note level and note elevel > 0 means add job to notebook. Good to prevent a lot of crap from piling up!
			if "NOTELEVEL" in toplevel: qtreewidget.setNoteLevel(toplevel["NOTELEVEL"])
			# Optional wikipage, to load info about the program from the wiki
			if "WIKIPAGE" in toplevel: qtreewidget.setWikiPage(toplevel["WIKIPAGE"])
			# Optional wizardfile, add a wizard to use with an e2program
			if "WIZARD" in toplevel: qtreewidget.setWizardFile(toplevel["WIZARD"])
			# Optional expertmode, to enable expert mode GUI widgets
			if "EXPERT" in toplevel:  qtreewidget.setExpertMode(toplevel["EXPERT"])
			self._add_children(toplevel, qtreewidget)
			QTree.addTopLevelItem(qtreewidget)

		QTree.itemClicked[QtWidgets.QTreeWidgetItem, int].connect(self._tree_widget_click)

		return QTree

	def json_strip_comments(self, data):
		"""This method takes a JSON-serialized string and removes
		JavaScript-style comments. These include // and /* */"""
		r = re.compile('/\\*.*\\*/', flags=re.M|re.S)
		data = r.sub("", data)
		data = re.sub("\s//.*\n", "", data)
		return data

	def _tree_widget_click(self, item, col):
		# Display the program GUI
		if item.getProgram():
			self._set_GUI(item.getProgram(), item.getMode())
			self.updateProject()
		elif item.getTable():
			# Only one table is allowed to be displayed at a time
			if self.table: self.table.close()
			self.table = eval(item.getTable())
			self.table.show()
		else:
			self.gui_stacked_widget.setCurrentIndex(0)
			self.clearE2Interface()


	def _set_GUI(self, program, mode):
		"""
		Set the current GUI widget
		"""

		programfile = program
		# Need a widget for each program for each mode
		if self.getProgramExpertMode() == 2:
			program = "Expert"+program+mode # There is a separate widget for the advanced mode
		else:
			program = program+mode
		# Only generate the GUI widget once.....
		if program in self.stackedWidgetHash:
			# load the widget using value from the correct DB
			self.gui_stacked_widget.setCurrentIndex(self.stackedWidgetHash[program])
			self.gui_stacked_widget.widget(self.stackedWidgetHash[program]).updateWidget()
		else:
			# Or make the widget
			self.stackedWidgetHash[program] = self.gui_stacked_widget.count()
			guioptions = self._read_e2program(programfile, mode)
			# Now actually make the widget
			widget = PMProgramWidget(guioptions, programfile, self, mode)
			self.gui_stacked_widget.addWidget(widget)
			self.gui_stacked_widget.setCurrentIndex(self.stackedWidgetHash[program])

	def _read_e2program(self, e2program, mode):
		"""
		This a pseudo import function to load parser info
		"""
		parser = EMArgumentParser()
		try:
			f = open(e2getinstalldir()+"/bin/"+e2program,"r")
		except:
			self.statusbar.setMessage("Can't open file '%s'"%e2program,"color:red;")
			return
		# Regex objects
		lineregex = re.compile("^\s*parser\.add_",flags=re.I) # eval parser.add_ lines, which are  not commented out.
		moderegex = re.compile("mode\s*=\s*[\"'].*%s[{0,1}.*]{0,1}.*[\"']"%mode,flags=re.I)	# If the program has a mode only eval lines with the right mode.
		modedefre = re.compile("%s\[([\w\.\(\)\']*)\]"%mode,flags=re.I)
		defaultre = re.compile("default\s*=\s*[^,]*")

		# Read line and do preprocessing (set mode defaults if desired)
		for line in f:
			if mode:
				if not re.search(moderegex, line): continue	# If we are running the program in a mode, then only eval mode lines
				string = re.findall(modedefre, re.findall(moderegex, line)[0])
				if string:
					default = re.findall(defaultre, line)
					if default:
						repl = "default=%s"%string[0]
						line = re.sub(defaultre, repl, line)
			if re.search(lineregex, line):
				eval(line)
				continue
			if 'parser.parse_args()' in line:
				break
		f.close()
		return parser.getGUIOptions()

	def makeTreeWidget(self, jsonfile, treename):
		"""
		Make a tree for a mode
		"""
		return self.loadTree(jsonfile, treename)

	def getCurrentTree(self):
		"""
		Return current tree
		"""
		return self.tree_stacked_widget.currentWidget()

	def getProgram(self):
		"""
		Return the current program
		"""
		try:
			return self.tree_stacked_widget.currentWidget().currentItem().getProgram()
		except:
			return None

	def getProgramMode(self):
		"""
		Return the current mode of the current program
		"""
		try:
			return self.tree_stacked_widget.currentWidget().currentItem().getMode()
		except:
			return None

	def getProgramNoteLevel(self):
		"""
		Return the program note level, used for noting
		"""
		try:
			return self.tree_stacked_widget.currentWidget().currentItem().getNoteLevel()
		except:
			return None

	def getProgramWikiPage(self):
		"""
		Return the program wiki page
		"""
		try:
			return self.tree_stacked_widget.currentWidget().currentItem().getWikiPage()
		except:
			return None

	def getProgramWizardFile(self):
		"""
		Return the program wiki page
		"""
		try:
			return self.tree_stacked_widget.currentWidget().currentItem().getWizardFile()
		except:
			return None

	def getProgramExpertMode(self):
		"""
		Return expert mode 0 = not available, 1 = available but not used, 2 = available and used
		"""
		try:
			return self.tree_stacked_widget.currentWidget().currentItem().getExpertMode()
		except:
			return None

	def setProgramExpertMode(self, state):
		"""
		Set expert mode
		"""
		self.tree_stacked_widget.currentWidget().currentItem().setExpertMode(state)

	def getCS(self):
		""" Return the project CS """
		try:
			return self.pm_projects_db.setdefault("global.microscope_cs",dfl=2.0)
		except:
			return ""

	def getInvar(self):
		""" Return the project invariant type """
		try:
			return self.pm_projects_db.setdefault("global.invariant_type",dfl="bispec")
		except:
			return ""

	def getVoltage(self):
		""" Return the project Voltage """
		try:
			return self.pm_projects_db.setdefault("global.microscope_voltage",dfl=200.0)
		except:
			return ""

	def getAPIX(self):
		""" Return the project Apix """
		try:
			return self.pm_projects_db.setdefault("global.apix",dfl=1.0)
		except:
			return ""

	def getMass(self):
		""" Return the particle mass """
		try:
			return self.pm_projects_db.setdefault("global.particle_mass", dfl=500.0)
		except:
			return ""

	def getPMCWD(self):
		""" Return the CWD that the pm is working in """
		return self.pm_cwd

	def changeDirectory(self, directory):
		"""
		Change directory we are working in
		"""
		os.chdir(directory)
		self.pm_cwd = directory
		# Update widget to reflect new CWD
		if self.gui_stacked_widget.currentIndex() > 0:
			self.gui_stacked_widget.currentWidget().updateWidget()

	def updateProject(self):
		"""
		Update the project. Make sure buttons, DB, title and icon observe the PM state
		"""
		# Set name, use db value, otherwise whatever the default is set to
		self.pn_project_name = self.pm_projects_db.setdefault("project_name", dfl=self.pn_project_name_default)
		self.PMProjectNameBanner.setText("Project Name: "+self.pn_project_name)

		# Set icon, use db value, otherwise whatever the default is set to
		self.pm_icon = self.pm_projects_db.setdefault("project_icon", dfl=self.pm_icon_default)
		self.PMIcon.setIcon(self.pm_icon)

		# Set buttons
		self.logbutton.setDown(bool(self.notebook) and self.notebook.isVisible(), quiet=True)
		self.taskmanagerbutton.setDown(bool(self.taskmanager) and self.taskmanager.isVisible(), quiet=True)
		self.helpbutton.setDown(bool(self.thehelp) and self.thehelp.isVisible(), quiet=True)

		# Wikipage
		if self.getProgramWikiPage():
			self.wikibutton.setEnabled(True)
		else:
			self.wikibutton.setEnabled(False)

		# Wizardpage
		if self.getProgramWizardFile():
			self.wizardbutton.setEnabled(True)
		else:
			self.wizardbutton.setEnabled(False)

		# Set Expert mode button
		if self.getProgramExpertMode() > 0:
			self.expertbutton.setEnabled(True)
			if self.getProgramExpertMode() == 2:
				self.expertbutton.setDown(True)
			else:
				self.expertbutton.setDown(False)
		else:
			self.expertbutton.setEnabled(False)
			
		

class EMPopen(subprocess.Popen):
	"""
	This was originally done to implement realtime STDOUT feed to the PM statusbar. Unfortunately this caused Broken Pipe errors causing PM to crash. I am experimenting with named pipes
	but this is still a work in progress
	"""
	def __init__(self, args, bufsize=0, executable=None,stdin=None, stdout=None, stderr=None,preexec_fn=None, close_fds=False, shell=False,cwd=None, env=None, universal_newlines=False,startupinfo=None, creationflags=0):
		subprocess.Popen.__init__(self, args, bufsize=bufsize, executable=executable, stdin=stdin, stdout=stdout, stderr=stderr, preexec_fn=preexec_fn, close_fds=close_fds, shell=shell, cwd=cwd, env=env, universal_newlines=universal_newlines, startupinfo=startupinfo, creationflags=creationflags)

	def realTimeCommunicate(self, msgbox):
		self.pmfifo = open("%s/pmfifo"%os.getcwd(),"r+")
		self.msgbox = msgbox
		rtcom = threading.Thread(target=self.realTimeChatter)
		rtcom.start()

	def realTimeChatter(self):
		while 1:
			if self.pmfifo:
				#pass
				stdout = self.pmfifo.readline()
				if stdout:
					self.msgbox.setMessage(stdout)
			else:
				break

			if self.poll() !=None:
				print("\n\nDONE\n\n")
				break

			# A HACK to prevent broken pipes (this may be a bit buggy), but it fixes an apparent bug in subprocess module
			# Some sort of syncronization issue
			time.sleep(0.1)

class EMAN2StatusBar(QtWidgets.QTextEdit):
	"""
	The status bar for PM
	"""
	def __init__(self, text, style):
		QtWidgets.QTextEdit.__init__(self)
		self.setFrameShape(QtWidgets.QFrame.Panel | QtWidgets.QFrame.Sunken)
		self.setLineWidth(2)
		#self.setContentsMargins(4, 4, 4, 4)
		self.setTextInteractionFlags(QtCore.Qt.NoTextInteraction)
		self.viewport().setCursor(QtCore.Qt.ArrowCursor)
		self.setMessage(text, style)
		self.setToolTip("This is the Status Bar")

	def setMessage(self, text, style=""):
		textcursor = self.textCursor()
		textcursor.movePosition(QtGui.QTextCursor.End)
		if textcursor.columnNumber() !=0:
			self.insertHtml("<br><span style=%s>%s</span>"%(style, text))
		else:
			self.insertHtml("<span style=%s>%s</span>"%(style,text))
		self.verticalScrollBar().setValue(self.verticalScrollBar().maximum())



class PMIcon(QtWidgets.QLabel):
	"""
	The icon manager for PM
	"""
	def __init__(self, image, parent=None):
		QtWidgets.QLabel.__init__(self, ("<img src=\"%s\" />")%image, parent)

	def setIcon(self, image):
		self.setText(("<img src=\"%s\" />")%image)

class EMWizard(QtWidgets.QWizard):
	"""
	This creates a wizard for filling out EM program forms.
	wizarddata is a list of dicts
	e2gui is a reference to a PMProgramWidget which is what this wizard is to fill out
	"""
	def __init__(self, wizarddata, e2gui):
		QtWidgets.QWizard.__init__(self)
		#self.setModal(True)
		self.setWindowTitle("e2program wizard")
		self.e2gui = e2gui

		for page in wizarddata:
			self.addPage(EMWizardPage(page, self.e2gui))


class EMWizardPage(QtWidgets.QWizardPage):
	"""
	Make a page for the wizard, uses input validation. If input is not valid send a message to the PM status bar and you will not be allowed to continue to the next wizard page
	"""
	def __init__(self, page, e2gui):
		QtWidgets.QWizardPage.__init__(self)
		self.e2gui = e2gui
		self.widgetlist = []

		#Set Hep info
		self.setTitle(page["TITLE"])
		label = QtWidgets.QLabel(page["INST"])
		label.setWordWrap(True)
		# add to grid
		grid = QtWidgets.QGridLayout()
		grid.addWidget(label,0,0)

		# Add widgets stuff
		frame = QtWidgets.QFrame()
		frame.setFrameStyle(QtWidgets.QFrame.StyledPanel)
		framegrid = QtWidgets.QGridLayout()
		framegrid.setContentsMargins(6,0,6,0)

		# Insanely if I were to just addWidget using self.e2gui.widgethash[widget], it steals it from the PMGUI layout, so I need to make a copy !!!
		for i, widget in enumerate(page["WIDGETS"]):
			widgetcopy = type(self.e2gui.widgethash[widget]).copyWidget(self.e2gui.widgethash[widget])
			self.widgetlist.append((self.e2gui.widgethash[widget],widgetcopy))
			framegrid.addWidget(widgetcopy,i,0)

		frame.setLayout(framegrid)
		grid.addWidget(frame,1,0)
		self.setLayout(grid)

	def validatePage(self):
		""" Validate input and update GUI"""
		for widget in self.widgetlist:
			if widget[1].getErrorMessage():
				#print widget[1].getErrorMessage()
				self.e2gui.pm().statusbar.setMessage(widget[1].getErrorMessage(),"color:red;")
				return False
			widget[0].setValue(widget[1].getValue())
		return True

class TheHelp(QtWidgets.QWidget):
	"""
	A little widget to aid in the daily chores. Good help is hard to find these days.....
	This is a GUI to e2help.py. Documentation is autogenerated from the EMAN2 C++ core
	"""
	def __init__(self, pm=None):
		QtWidgets.QWidget.__init__(self)
		if pm:
			self.pm = weakref.ref(pm)
		else:
			self.pm = None
		self.helptopics = []
		self.widgetgeometry = None

		self.setWindowTitle('The Help')
		grid = QtWidgets.QGridLayout()
		grid.addWidget(self.getToolBar(), 0, 0)
		self.textbox = QtWidgets.QTextEdit()
		self.textbox.setReadOnly(True)
		grid.addWidget(self.textbox, 1, 0)
		self.setLayout(grid)

		self.setMinimumWidth(600)
		self._helpchange(0)

	def getToolBar(self):
		""" Return the toolbar widget """
		tbwidget = QtWidgets.QWidget()
		grid = QtWidgets.QGridLayout()
		grid.setContentsMargins(6,0,6,0)

		font = QtGui.QFont()
		font.setBold(True)
		helplabel = QtWidgets.QLabel("EMAN2 topic:")
		helplabel.setFont(font)
		helplabel.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)

		hbox = QtWidgets.QHBoxLayout()
		hbox.addWidget(PMIcon(get_image_directory() + "SirEMAN2.png"))
		hbox.addWidget(helplabel)
		hbox.setContentsMargins(0,0,0,0)
		grid.addLayout(hbox,0, 0)
		self.helpcb = QtWidgets.QComboBox()
		grid.addWidget(self.helpcb, 0, 1, 1, 2)

		self.helpcb.addItem("aligners")
		self.helpcb.addItem("averagers")
		self.helpcb.addItem("analyzers")
		self.helpcb.addItem("cmps")
		self.helpcb.addItem("orientgens")
		self.helpcb.addItem("processors")
		self.helpcb.addItem("projectors")
		self.helpcb.addItem("reconstuctors")
		self.helpcb.addItem("symmetries")
		self.helptopics.append(["aligners", dump_aligners_list()])
		self.helptopics.append(["averagers", dump_averagers_list()])
		self.helptopics.append(["analyzers", dump_analyzers_list()])
		self.helptopics.append(["cmps", dump_cmps_list()])
		self.helptopics.append(["orientgens", dump_orientgens_list()])
		self.helptopics.append(["processors", dump_processors_list()])
		self.helptopics.append(["projectors", dump_projectors_list()])
		self.helptopics.append(["reconstuctors", dump_reconstructors_list()])
		self.helptopics.append(["symmetries", dump_symmetries_list()])


		self.helpcb.activated[int].connect(self._helpchange)

		tbwidget.setLayout(grid)
		return tbwidget

	def _helpchange(self, idx):
		""" Read in documenation info from the dump_.*_list() functions. Should reflect what is in the C dicts """
		self.helpcb.setCurrentIndex(idx)
		helpdict =  self.helptopics[idx][1]
		helpdoc = "<B><H3>Listed below is a list of EMAN2 <I>%s</I></H3></B><BR>"%self.helptopics[idx][0]

		keys = list(helpdict.keys())
		keys.sort()
		for key in keys:
			helpdoc += "<B>%s</B>"%(key)
			eman2item = helpdict[key]
			helpdoc += "<UL><LI><I>Description:</I> %s</LI>"%eman2item[0]
			for param in range(old_div((len(eman2item)-1),3)):
				helpdoc += "<LI><I>Parameter:</I> &nbsp;<B>%s(</B><SPAN style='color:red;'>%s</SPAN><B>)</B>, %s</LI>"%(eman2item[param*3 +1],eman2item[param*3 +2],eman2item[param*3 +3])
			helpdoc += "</UL>"

		self.textbox.setHtml(helpdoc)

	def hideEvent(self, event):
		""" This remembers the geometry when we hide the widget """
		QtWidgets.QWidget.hideEvent(self, event)
		self.widgetgeometry = self.geometry()

	def showEvent(self, event):
		""" This recalls the geometry when we show the widget """
		QtWidgets.QWidget.showEvent(self, event)
		if self.widgetgeometry: self.setGeometry(self.widgetgeometry)

	def closeEvent(self, event):
		if self.pm: self.pm().thehelp = None
		if self.pm: self.pm().updateProject()

class NoteBook(QtWidgets.QWidget):
	"""
	The Notebook for PM. The note book will reflect top levels jobs run, even if they were run on the command line. This class needs some work, b/c the notebook exhibits some funny font and color
	behaviour. This is probably because the cursor is not a strict observer of the NoteBook state.
	When jobs are launched the command is recorded in this widget. In addition this widget saves its text as an HTML file named. pmnotes.html in the local project directory
	"""
	def __init__(self, pm):
		QtWidgets.QWidget.__init__(self)
		self.pm = weakref.ref(pm)
		self.donotsave = False
		self.widgetgeometry = None

		self.setWindowTitle('NoteBook')
		grid = QtWidgets.QGridLayout()
		font = QtGui.QFont()
		font.setBold(True)
		textlabel = QtWidgets.QLabel("EMAN2 NoteBook")
		textlabel.setFont(font)
		self.texteditbox = PMTextEdit(self)
		grid.addWidget(textlabel,0,0)
		grid.addWidget(self.getToolBar(),1,0,1,2)
		grid.addWidget(self.texteditbox,2,0,1,2)
		self.savepb = QtWidgets.QPushButton("Save")
		self.closepb = QtWidgets.QPushButton("Close")
		grid.addWidget(self.savepb, 3,0)
		grid.addWidget(self.closepb, 3,1)
		self.setLayout(grid)
		self.setMinimumWidth(600)

		# Load the logbook and update
		self.loadNoteBook()
		self.checkEMAN2LogFile()
		self.writeNotes() # save changes from checkEMAN2LogFile()

		self.savepb.clicked.connect(self._on_save)
		self.closepb.clicked.connect(self._on_close)

	def getToolBar(self):
		""" Return the toolbar widget """
		tbwidget = QtWidgets.QWidget()
		hbox = QtWidgets.QHBoxLayout()
		self.dbdict = js_open_dict(self.pm().getPMCWD()+"/info/notebook.json")

		# font type
		self.fontdb = QtGui.QFontDatabase()
		self.fontfamily = QtWidgets.QComboBox()
		self.fontfamily.addItems(self.fontdb.families())

		# font size
		self.fontsizecb = QtWidgets.QComboBox()

		# Bold italic, underline
		self.boldbutton = PMToolButton()
		self.italicbutton = PMToolButton()
		self.underlinebutton = PMToolButton()
		self.boldbutton.setIcon(QtGui.QIcon(QtGui.QPixmap(boldicon)))
		self.italicbutton.setIcon(QtGui.QIcon(QtGui.QPixmap(italicicon)))
		self.underlinebutton.setIcon(QtGui.QIcon(QtGui.QPixmap(underlineicon)))
		self.fontcolor = EMQTColorWidget(red=0,green=0,blue=0)

		# Add widgets
		hbox.addWidget(self.fontfamily)
		hbox.addWidget(self.fontsizecb)
		hbox.addWidget(self.boldbutton)
		hbox.addWidget(self.italicbutton)
		hbox.addWidget(self.underlinebutton)
		hbox.addWidget(self.fontcolor)
		tbwidget.setLayout(hbox)

		# Set defaults
		if 'FONTFAMILY' in self.dbdict:
			idx = self.fontfamily.findText(self.dbdict['FONTFAMILY'])
			self.fontfamily.setCurrentIndex(idx)
			self._load_fontsizes()
			self.fontsizecb.setCurrentIndex(self.fontsizecb.findText(self.dbdict['FONTSIZE']))
			self.boldbutton.setDown(bool(self.dbdict['BOLDFONT']), quiet=1)
			self.italicbutton.setDown(bool(self.dbdict['ITALICFONT']), quiet=1)
			self.underlinebutton.setDown(bool(self.dbdict['UNDERLINE']), quiet=1)
			self.fontcolor.setColor(QtGui.QColor(self.dbdict['FONTCOLOR'][0],self.dbdict['FONTCOLOR'][1],self.dbdict['FONTCOLOR'][2]))
			self.texteditbox.updateFont()
		else:
			self._load_fontsizes()

		# Connect signals
		self.fontfamily.activated[int].connect(self._fontfamilychange)
		self.fontsizecb.activated[int].connect(self._fontchange)
		self.boldbutton.stateChanged[bool].connect(self._fontchange)
		self.italicbutton.stateChanged[bool].connect(self._fontchange)
		self.underlinebutton.stateChanged[bool].connect(self._fontchange)
		self.fontcolor.newcolor[QtGui.QColor].connect(self._fontchange)

		return tbwidget

	def _load_fontsizes(self):
		for i in self.fontdb.pointSizes(self.fontfamily.currentText()):
			self.fontsizecb.addItem(str(i))

	def _fontfamilychange(self):
		self._load_fontsizes()
		self._fontchange()

	def _fontchange(self):
		self.dbdict['FONTFAMILY'] = self.fontfamily.currentText()
		self.dbdict['FONTSIZE'] = self.fontsizecb.currentText()
		self.dbdict['BOLDFONT'] = self.boldbutton.isDown()
		self.dbdict['ITALICFONT'] = self.italicbutton.isDown()
		self.dbdict['UNDERLINE'] = self.underlinebutton.isDown()
		self.dbdict['FONTCOLOR'] = self.fontcolor.getColor().getRgb()
		self.texteditbox.updateFont()

	def loadNoteBook(self):
		""" Load the current log book """
		try:
			fin=open(self.pm().getPMCWD()+"/pmnotes.html","r")
			self.texteditbox.setHtml(fin.read())
			fin.close()
		except:
			self.pm().statusbar.setMessage("No pmnotes file found, starting a new one...","font-weight:bold;")
			self.texteditbox.insertHtml("<b>Project Name: %s<br>EMAN Version: %s<br>Path: %s<br></b>"%(self.pm().pn_project_name,EMANVERSION,self.pm().getPMCWD()))

	def checkEMAN2LogFile(self):
		""" Check the log file for any changes """
		try:
			fin=open(self.pm().getPMCWD()+"/.eman2log.txt","r+")
		except:
			return
		lastlineloc = 0
		while True:
			try:
				line = fin.readline()
				if len(line) == 0: raise
			except:
				break
			fields = line.split('\t')
			# For top level jobs ensure that they are in the pmnotes file
			if "-1" in fields[2]:
				self._checklogforjob(fields[4].strip(), fields[0].strip())
				# We replace -1 with -2 to indicate that this event has been recorded in the pm log book. So that we don't check over and over again...
				inipos = fin.tell()
				fin.seek(lastlineloc+51)
				fin.write("-2")
				fin.seek(inipos)

			lastlineloc = fin.tell()
		fin.close()

	def insertNewJob(self, job, time):
		""" Insert a new job into the Html """
		textcursor = self.texteditbox.textCursor()
		textcursor.movePosition(QtGui.QTextCursor.End)
		textcursor.insertText('\n')
		textcursor.insertHtml("<b>Job Launched on <i>%s</i>:<font color=blue> %s</font></b>"%(time, job))

	def writeNotes(self):
		""" Write the log file to disk """
		fout=open(self.pm().getPMCWD()+"/pmnotes.html","w")
		fout.write(self.texteditbox.toHtml())
		fout.close()

	def _checklogforjob(self, job, time):
		""" Check to see if this job is already in the html """
		if time in self.texteditbox.toHtml():
			return True
		else:
			self.insertNewJob(job, time)
		return False

	def _on_save(self):
		self.writeNotes()

	def _on_close(self):
		self.donotsave = True
		self.close()

	def hideEvent(self, event):
		""" This remembers the geometry when we hide the widget """
		QtWidgets.QWidget.hideEvent(self, event)
		self.widgetgeometry = self.geometry()

	def showEvent(self, event):
		""" This recalls the geometry when we show the widget """
		QtWidgets.QWidget.showEvent(self, event)
		if self.widgetgeometry: self.setGeometry(self.widgetgeometry)

	def closeEvent(self, event):
		if not self.donotsave: self.writeNotes()
		self.pm().notebook = None
		self.pm().updateProject()

class PMTextEdit(QtWidgets.QTextEdit):
	""" Sub class of the QTextEdit widget to observe the PMNoteBook """
	def __init__(self, parent):
		QtWidgets.QTextEdit.__init__(self)
		self.parent = weakref.ref(parent)

		self.setPMFontWeight(QtGui.QFont.Normal)
		self.setPMTextColor(QtGui.QColor(0,0,0))

	def mousePressEvent(self,event):
		""" Stupidly, TextEdit resets the font based on its context, which I find undesireable """
		QtWidgets.QTextEdit.mousePressEvent(self, event)
		self.setFontWeight(self.pmfontweight)
		self.setTextColor(self.textcolor)

	def setPMFontWeight(self, weight):
		self.pmfontweight = weight
		self.setFontWeight(weight)

	def setPMTextColor(self, color):
		self.textcolor = color
		self.setTextColor(color)

	def setPMFontSize(self, size):
		self.fontsize = size
		self.setFontPointSize(size)

	def setPMFontFamily(self, family):
		self.fontfamily = family
		self.setFontFamily(family)

	def setPMFontItalic(self, italic):
		self.italicfont = italic
		self.setFontItalic(italic)

	def setPMFontUnderline(self, underline):
		self.underlinefont = underline
		self.setFontUnderline(underline)

	def updateFont(self):
		""" When the NoteBook state is changed this function is called """
		if self.parent().boldbutton.isDown():
			self.setPMFontWeight(QtGui.QFont.Bold)
		else:
			self.setPMFontWeight(QtGui.QFont.Normal)
		try : self.setPMFontSize(int(self.parent().fontsizecb.currentText()))
		except: self.setPMFontSize(10)
		self.setPMTextColor(self.parent().fontcolor.getColor())
		self.setPMFontFamily(self.parent().fontfamily.currentText())
		self.setPMFontItalic(self.parent().italicbutton.isDown())
		self.setPMFontUnderline(self.parent().underlinebutton.isDown())

class TaskManager(QtWidgets.QWidget):
	"""
	The Task manager for PM, unlike e2workflow.py, this actualy works. Well sort of... It uses a SIGTERM soft kill, which when it works is nice because EMAN2 jobs are
	processed by EMAN2 eveent handlers, thus producing clean kills. Often though the job completely ignores SIGTERM and keeps on humming. Perhaps SIGTERM should
	be replaced with SIGKILL to hard kill it. Downside is that the kill could be dirty and corrupt. On the up side your jobs will actually be killed!
	"""
	def __init__(self, pm):
		QtWidgets.QWidget.__init__(self)
		self.pm = weakref.ref(pm)
		self.widgetgeometry = None
		self.setWindowTitle('Task Manager')

		grid = QtWidgets.QGridLayout()
		font = QtGui.QFont()
		font.setBold(True)
		textlabel = QtWidgets.QLabel("Running Tasks")
		textlabel.setFont(font)
		grid.addWidget(textlabel,0,0)
		self.list_widget = QtWidgets.QListWidget()
		grid.addWidget(self.list_widget,1,0,1,2)
		self.killpb = QtWidgets.QPushButton("Kill")
		self.closepb = QtWidgets.QPushButton("Close")
		grid.addWidget(self.killpb, 2,0)
		grid.addWidget(self.closepb, 2,1)
		self.setLayout(grid)
		self.resize(400, 200)


		self.closepb.clicked.connect(self._on_close)
		self.killpb.clicked.connect(self._on_kill)
		self.tasks=None
		self.update_tasks()

		# A timer for updates
		self.timer = QtCore.QTimer(self);
		self.timer.timeout.connect(self.update_tasks)
		self.timer.start(2000)

	def check_task(self,fin,ptsk):
		"""Note that this modifies ptsk in-place"""
		fin.seek(ptsk[0])
		try : t=fin.readline()
		except: return False

		tsk=t.split("\t")
		ptsk[3]=tsk[1]
		if self.is_running(tsk, ptsk[0]) : return True

		return False

	def is_running(self,tsk,loc=None):
		"""Check to see if the 'tsk' string is actively running"""
		try:
			# This means the task must be done
			if tsk[1][4]=="/" or "crashed" in tsk[1]: return False

			# Or not running
			if '/' in tsk[2] : pid=int(tsk[2].split("/")[0])
			else : pid=int(tsk[2])

			if os.name=="posix":
				# Don't do this on windows (may actually kill the process)
				os.kill(pid,0)			# raises an exception if doesn't exist
			else :
				# Not a good approach, but the best we have right now on windows
				tm=time.localtime()
				mon=int(tsk[0][5:7])
				yr=int(tsk[0][:4])
				day=int(tsk[0][8:10])
				if mon!=tm[1] or yr!=tm[0] or tm[2]-day>1 : raise Exception

			if os.name=="nt": tsk[4]=tsk[4].replace("\\","/")
			command = tsk[4].split()[0].split("/")[-1]
			if command=="e2projectmanager.py" : raise Exception
		except:
			# If the process is not running nor complete then list it as crashed
			if  loc:
				fout=open(self.pm().getPMCWD()+"/.eman2log.txt","r+")
				fout.seek(loc+20)
				fout.write("crashed            ")
				fout.close()
			return False

		return True

	def parse_new_commands(self,fin):
		self.last_update=time.time()
		while True:
			loc=fin.tell()
			try:
				t=fin.readline()
				if len(t)==0 : raise Exception
			except:
				self.lastpos=loc
				break

			tsk=t.split("\t")

			if not self.is_running(tsk,loc) : continue

			# if we get here, then we have a (probably) active task

			# This contains (file location, process id, status, command name)
			ppid = None
			if '/' in tsk[2]:
				ids = tsk[2].split("/")
				pid=int(ids[0])
				if len(ids) == 2:
					ppid = int(ids[1].strip())
			else : pid=int(tsk[2])
			if os.name=="nt": tsk[4]=tsk[4].replace("\\","/")
			command = tsk[4].split()[0].split("/")[-1]
			self.tasks.append([loc,pid,tsk[0],tsk[1],command,ppid])

	def update_tasks(self):
		if self.tasks==None:
			try: fin=open(self.pm().getPMCWD()+"/.eman2log.txt","r")
			except: return
			self.prevfin = fin
			self.tasks=[]
			self.parse_new_commands(fin)
			if len(self.tasks)==0 :
				self.tasks=None
				return
		else:
			try:
				if os.stat(self.pm().getPMCWD()+"/.eman2log.txt").st_mtime < self.last_update:
					# See if any jobs have crashed, if not then do nothing
					self.tasks=[i for i in self.tasks if self.check_task(self.prevfin,i)]
					self.update_list()
					return

			except:
				print("Error, couldn't stat(.eman2log.txt).")
				return

			# Load in the file
			try: fin=open(self.pm().getPMCWD()+"/.eman2log.txt","r")
			except: return
			self.prevfin = fin

			# Go through existing entries and see if they've finished
			self.tasks=[i for i in self.tasks if self.check_task(fin,i)]

			# Now check for new entries
			fin.seek(self.lastpos)
			self.parse_new_commands(fin)
		# Update the list
		self.update_list()

	def update_list(self):
		#self.list_widget.clear()
		listitems = self.getListItems()
		for t in self.tasks:
			if t[1] in list(listitems.keys()):
				listitems[t[1]].setText("%s  %s (%s)  %s"%(t[2][5:16],t[3][:4],t[1],t[4]))
				del(listitems[t[1]])
				continue
			listwigetitem = PMQListWidgetItem("%s  %s (%s)  %s"%(t[2][5:16],t[3][:4],t[1],t[4]))
			listwigetitem.setPID(t[1])
			listwigetitem.setPPID(t[5])
			listwigetitem.setProgramName(t[4])
			self.list_widget.addItem(listwigetitem)
		# Remove items that have stopped running
		for item in list(listitems.values()):
			self.list_widget.takeItem(self.list_widget.row(item))

	def getListItems(self):
		itemdict = {}
		for i in range(self.list_widget.count()):
			itemdict[self.list_widget.item(i).getPID()] = self.list_widget.item(i)
		return itemdict

	def reset(self):
		self.tasks = None
		self.update_tasks()

	def _on_close(self):
		self.close()

	def _on_kill(self):
		killsig=signal.SIGTERM
		modifiers = QtWidgets.QApplication.keyboardModifiers()
		if modifiers == QtCore.Qt.ShiftModifier:
			print("Shift held. Will force kill processes")
			killsig=signal.SIGKILL

		selitems = self.list_widget.selectedItems()
		self.update_tasks()
		for item in selitems:
			curpid=item.getPID()
			curppid=item.getPPID()
			print("Killing process PID: {}, PPID: {}".format(curpid, curppid))
			if os.name == 'posix': # Kill by this method will only work for unix/linux type OS
				
				if curppid<0:
					#### this is a parent process. kill itself then all its children
					os.kill(curpid,killsig)
					self._recursivekill(curpid)
				else:
					#### this is a child process. kill all its siblings then parent
					for task in self.tasks:
						if curppid == task[5]: 
							#### processes share the same parent
							print("Killing self process", task[1])
							os.kill(task[1], killsig)
							self._recursivekill(task[1], killsig)
					## now parent
					print("Killing parent")
					os.kill(curppid,killsig)
					
				## The PID Mafia occurs below: # Kill all children and siblings
				#for task in self.tasks:
					#if item.getPPID() == task[5]:
						#print("killing self process", task[1])
						#os.kill(task[1],signal.SIGTERM)
						#self._recursivekill(task[1])
				## kill parent (top level item)
				#if item.getPPID() > 0:
					#print("KIlling parent")
					#os.kill(item.getPPID(),signal.SIGTERM)
			else:
				# Windows kill
				pass

	def _recursivekill(self, pid, sig=signal.SIGTERM):
		for task in self.tasks:
			if pid == task[5]:
				print("Killing child process", task[1])
				os.kill(task[1],sig)
				self._recursivekill(task[1])

	def hideEvent(self, event):
		""" This remembers the geometry when we hide the widget """
		QtWidgets.QWidget.hideEvent(self, event)
		self.widgetgeometry = self.geometry()

	def showEvent(self, event):
		""" This recalls the geometry when we show the widget """
		QtWidgets.QWidget.showEvent(self, event)
		if self.widgetgeometry: self.setGeometry(self.widgetgeometry)

	def closeEvent(self, event):
		self.pm().taskmanager = None
		self.pm().updateProject()

class PMQListWidgetItem(QtWidgets.QListWidgetItem):
	""" A subclass of the QListWidgetItem for use in the task manger listwidget"""
	def __init__(self, text):
		QtWidgets.QListWidgetItem.__init__(self, text)
		self.pid = None
		self.programname = None

	def setPID(self, pid):
		self.pid = pid

	def getPID(self):
		return self.pid

	def setPPID(self, ppid):
		self.ppid = ppid

	def getPPID(self):
		return self.ppid

	def setProgramName(self, name):
		self.programname = name

	def getProgramName(self):
		return self.programname

class PMProgramWidget(QtWidgets.QTabWidget):
	"""
	Creates a program interface for each e2 program, for each mode, etc. This is a tab widget and there are three tabs, a GUI tab, a command line tab and a help tab
	"""
	def __init__(self, options, program, pm, mode):
		QtWidgets.QTabWidget.__init__(self)
		self.pm = weakref.ref(pm)
		self.setMinimumHeight(210) # Size of the tool bar

		# Add GUI tab
		self.guiwidget = PMGUIWidget(options, program, pm, mode)
		self.addTab(self.guiwidget, "GUI")

		# Add command tab
		self.guitexteditbox = QtWidgets.QTextEdit("")
		self.guitexteditbox.setWordWrapMode(QtGui.QTextOption.WrapAnywhere)
		self.addTab(self.guitexteditbox, "Command")

		# Add help tab
		self.helptexteditbox = QtWidgets.QTextEdit("")
		self.helptexteditbox.setWordWrapMode(QtGui.QTextOption.WordWrap)
		self.helptexteditbox.setReadOnly(True)
		self.helptexteditbox.viewport().setCursor(QtCore.Qt.ArrowCursor)
		self.helptexteditbox.setText(self.pm().loadUsage(self.pm().getProgram()))
		self.addTab(self.helptexteditbox, "Help")

		self.previoustab = 0

		self.currentChanged[int].connect(self._on_tabchange)

	def updateWidget(self):
		""" Delegate to guiwidget """
		self.guiwidget.updateWidget()

	def getCommand(self):
		""" Delegate to gui widget """
		if self.currentIndex() == 0:
			return self.guiwidget.getCommand()
		elif self.currentIndex() == 1:
			if not self.guiwidget.errorstate:
				return self.guitexteditbox.toPlainText()
			else:
				self.pm().statusbar.setMessage("Error(s) in GUI parameters. Please fix....","color:red;")
				return None
		else:
			self.pm().statusbar.setMessage("Error: Can't launch from help menu","color:red;")
			return None

	def _on_tabchange(self, idx):
		""" Implements cross talk between GUI and Comnmand widgets """
		if idx == 1:
			errormsg = self.guiwidget.getErrorMessage()
			if errormsg:
				self.guitexteditbox.setHtml("Error in parameters<br>"+errormsg)
			else:
				self.guitexteditbox.setHtml(self.guiwidget.getCommand())
		if idx == 0 and self.previoustab == 1:
			self.guiwidget.updateGUIFromCmd(str(self.guitexteditbox.toPlainText()))

		self.previoustab = idx

class PMGUIWidget(QtWidgets.QScrollArea):
	"""
	Creates a GUI widget using a dict derived from the e2program options. Instances of this widget are added to the QStackedWidget on the right hand side of the PM.
	When the user clicks on a leaf node in the workflow tree an instance of this class is created, if it doesn't already exist, and added to the stackedwidget. If it
	already exists, then it is rendered visible.
	"""
	def __init__(self, options, program, pm, mode):
		QtWidgets.QScrollArea.__init__(self)
		self.errorstate = False
		# I need both an ordered list and an associative means of accessing the widgets
		self.widgetlist = []
		self.widgethash = {}
		self.cwd  = pm.pm_cwd	# The working directory that were in
		self.program = program
		self.db = js_open_dict("{}/info/pm/{}.json".format(str(self.cwd),self.program))
		self.pm = weakref.ref(pm)
		self.mode = mode

		# Loop through options (a list of dicts) and generate the GUI widget
		gridbox = QtWidgets.QGridLayout()
		for option in options:
			"""create the correct widget type"""
			if ('expert' in  option) and self.pm().getProgramExpertMode() < 2: continue	# Do no add expertmode options if not in expert mode (defined as 2, available AND turned on) Should have used ENUM here

			if option['guitype'] == 'header':
				widget = PMHeaderWidget(option['name'], option['title'])
			if option['guitype'] == 'filebox':
				widget = PMFileNameWidget(option['name'], self.getDefault(option), self.getSharingMode(option), self.getBrowser(option), postional=self.getPositional(option), initdefault=self.getDefault(option, nodb=True),checkfileexist=self.getFileCheck(option), infolabels=self.getInfoLabels(option))
			if option['guitype'] == 'dirbox':
				widget = PMDirectoryWidget(option['name'], option['dirbasename'], self.getDefault(option), self.getSharingMode(option), postional=self.getPositional(option), initdefault=self.getDefault(option, nodb=True))
			if option['guitype'] == 'symbox':
				widget = PMSymWidget(option['name'], self.getDefault(option), self.getSharingMode(option), initdefault=self.getDefault(option, nodb=True))
			if option['guitype'] == 'intbox':
				widget = PMIntEntryWidget(option['name'], self.getDefault(option), self.getSharingMode(option), self.getLRange(option), self.getURange(option), postional=self.getPositional(option), initdefault=self.getDefault(option, nodb=True))
			if option['guitype'] == 'shrinkbox':
				widget = PMShrinkEntryWidget(option['name'], self.getDefault(option), self.getSharingMode(option), 2, postional=self.getPositional(option), initdefault=self.getDefault(option, nodb=True))
			if option['guitype'] == 'floatbox':
				widget = PMFloatEntryWidget(option['name'], self.getDefault(option), self.getSharingMode(option), self.getLRange(option), self.getURange(option), postional=self.getPositional(option), initdefault=self.getDefault(option, nodb=True))
			if option['guitype'] == 'boolbox':
				widget = PMBoolWidget(option['name'], self.getDefault(option), self.getSharingMode(option), initdefault=self.getDefault(option, nodb=True))
			if option['guitype'] == 'strbox':
				widget = PMStringEntryWidget(option['name'], self.getDefault(option), self.getSharingMode(option), postional=self.getPositional(option), initdefault=self.getDefault(option, nodb=True),returnNone=self.getreturnNone(option))
			if option['guitype'] == 'comboparambox':
				widget = PMComboParamsWidget(option['name'], self.getChoices(option), self.getDefault(option), self.getSharingMode(option), postional=self.getPositional(option), initdefault=self.getDefault(option, nodb=True),returnNone=self.getreturnNone(option))
			if option['guitype'] == 'combobox':
				widget = PMComboWidget(option['name'], self.getChoices(option), self.getDefault(option), self.getSharingMode(option), postional=self.getPositional(option), initdefault=self.getDefault(option, nodb=True),returnNone=self.getreturnNone(option))
			if option['guitype'] == 'automask3d':
				widget = PMAutoMask3DWidget(option['name'], self.getDefault(option), self.getSharingMode(option), initdefault=self.getDefault(option, nodb=True))
			if option['guitype'] == 'fsctable':
				widget = PMFSCTableWidget(option['name'], self.getDefault(option), self.getSharingMode(option), postional=self.getPositional(option), initdefault=self.getDefault(option, nodb=True))

			# Setup each widget
			widget.pmmessage[str].connect(self._on_message)
			widget.setToolTip(option['help'])
			self.widgethash[option['name']] = widget
			self.widgetlist.append(widget)
			gridbox.addWidget(widget, option['row'], option['col'], self.getRowSpan(option), self.getColSpan(option))

		# Now make a widget and add it to the scroll area
		gridbox.setContentsMargins(0,0,0,0)
		self.scwidget = QtWidgets.QWidget()
		self.scwidget.setLayout(gridbox)
		self.scwidget.setMinimumWidth(self.width()-1.5*self.verticalScrollBar().width())
		self.scwidget.setMaximumWidth(self.width()-1.5*self.verticalScrollBar().width())
		self.setWidget(self.scwidget)
		self.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
		self.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
		self.setSizePolicy(QtWidgets.QSizePolicy.Fixed,QtWidgets.QSizePolicy.Fixed)

	def getRowSpan(self, option):
		"""Return the rowspan"""
		rowspan = 1
		if 'rowspan' in option: rowspan = option['rowspan']
		return rowspan

	def getColSpan(self, option):
		"""Return the rowspan"""
		colspan = 1
		if 'colspan' in option: colspan = option['colspan']
		return colspan

	def getDefault(self, option, nodb = False):
		""" Return the default value according to the following rules"""
		# If there is a DB and its usage is desired the default will be the DB value
		k=option['name']+self.getSharingMode(option)
		if not nodb and k in self.db: return self.db[k]	# Return the default if it exists in the DB
		default = ""
		if 'default' in option: default = option['default']
		if type(default) == str and "self.pm()" in default: default = eval(default)	# eval CS, apix, voltage, etc, a bit of a HACK, but it works
		return default

	def getLRange(self, option):
		""" Return the lower range if it exists"""
		lrange = None
		if 'lrange' in option: lrange = option['lrange']
		return lrange

	def getURange(self, option):
		""" Return the upper range if it exists"""
		urange = None
		if 'urange' in option: urange = option['urange']
		return urange

	def getChoices(self, option):
		""" Get a list of choices for combo boxes """
		choices = []
		if 'choicelist' in option:
			choices = eval(option['choicelist'])
			# If it is a dict, get the keys
			if type(choices) == type({}):
				choices = list(choices.keys())
		return choices

	def getPositional(self, option):
		# See if the arugment is positional or not
		positional = False	# Default if not provided
		if 'positional' in option: positional = option['positional']
		return positional

	def getFileCheck(self, option):
		filecheck = True
		if 'filecheck' in option: filecheck = option['filecheck']
		return filecheck
	
	
	def getInfoLabels(self, option):
		infolabel = False
		if 'infolabel' in option: infolabel = option['infolabel']
		return infolabel
	

	def getSharingMode(self, option):
		""" Returns whether or not the widget desires DB sharing with other modes """
		if 'nosharedb' in option:
			return self.mode
		else:
			return ""

	def getBrowser(self, option):
		""" Returns the sort of browser to use """
		browser = "EMBrowserWidget(withmodal=True,multiselect=True)"
		if 'browser' in option:
			browser = option['browser']
		return browser

	def getreturnNone(self, option):
		""" Sets whether or not we will actually return None as an argument or leave it blank. Only some args can accept None, others crash if this is input """
		if 'returnNone' in option:
			return True
		else:
			return False

	def updateWidget(self):
		# Reload the DB if necessary (when projects are changed)
		thiscwd = self.pm().pm_cwd
		if thiscwd != self.cwd:
			self.cwd = thiscwd
			self.db = js_open_dict("{}/info/pm/{}.json".format(str(self.cwd),self.program))
		# Might enable serval dbs to be loaded, but we will implement this later
		for widget in self.widgetlist:
			# If this is not a value holding widget continue (bool widget, herader, etc)
			if widget.getArgument() == None: continue
			try: widget.setValue(self.db[widget.getName()+widget.getMode()], quiet=True)
			except: widget.setValue(widget.initdefault, quiet=True)

	def updateGUIFromCmd(self, cmd):
		""" Update the GUI from a command """
		if self.errorstate or not cmd: return
		options = re.findall('\s--\S*',cmd)
		args = re.findall('\s[^-]{2}\S*',cmd)
		widgethash = dict(self.widgethash)	# Make a copy of the widget hash

		# Process
		for option in options:
			ov = option.split('=', 1)
			if ov[0][3:] not in self.widgethash:
				self.pm().statusbar.setMessage("Rubbish!!! Option '%s' not found."%ov[0][3:],"color:red;")
				continue
			if len(ov) == 2:
				self._setValueJournaling(self.widgethash[ov[0][3:]], ov[1])
			#self.widgethash[ov[0][2:]].setValue(ov[1])
			else:
				self._setValueJournaling(self.widgethash[ov[0][3:]], True)
			#self.widgethash[ov[0][2:]].setValue(True)
			# Pop the widget off a copy of the hash
			del(widgethash[ov[0][3:]])

		posargs = []
		# Now do the widgets which are not listed in the above list
		for name,widget in list(widgethash.items()):
			if isinstance(widget, PMHeaderWidget):
				continue
			if isinstance(widget, PMBoolWidget):
				self._setValueJournaling(widget, False)
				continue
			# Process arguments, if positional widget
			if widget.getPositional():
				posargs.append(widget)
				continue

			# Otherwise set to none	because the user deleted it
			self._setValueJournaling(widget, '')

		# Loop through pos args. The standard is that either there is one arg per widget or one widget and many args
		if len(posargs) == len(args):
			for posarg in posargs:
				self._setValueJournaling(posarg, (args.pop())[1:])
		else:
			argsstring = ''.join(str(n ) for n in args)
			self._setValueJournaling(posargs[0], argsstring)

	def _setValueJournaling(self, widget, value):
		""" Set a value, but if an error occurs then revert to the old value """
		origval = widget.getValue()
		widget.setValue(value)
		if widget.getErrorMessage():
			widget.setValue(origval)


	def getCommand(self):
		""" Loop and check for errors and set the DB. If errors are encountered, then return None """
		self.errorstate = False
		args = self.program
		for widget in self.widgetlist:
			# If this is not a value holding widget continue
			if widget.getArgument() == None: continue
			# Check for errors before we launch script
			if widget.getErrorMessage():
				self.pm().statusbar.setMessage(self.getErrorMessage(),"color:red;")
				self.errorstate = True
				return None
			# Save the value
			self.db[widget.getName()+widget.getMode()] = widget.getValue()
			args += " "+widget.getArgument()

		return args

	def getErrorMessage(self):
		""" Check for any error messages """
		errormsg = ""
		self.errorstate = False
		for widget in list(self.widgethash.values()):
			if widget.getErrorMessage():
				self.errorstate = True
				errormsg += (widget.getErrorMessage()+"<br>")
		return errormsg

	def _on_message(self, s):
		self.pm().statusbar.setMessage(str(s),"color:red;")

class PMQTreeWidgetItem(QtWidgets.QTreeWidgetItem):
	"""
	Custom QTreeWidget for PM, holds a bunch of properties relating to the e2program (or table) it represents.
	"""
	def __init__(self, qstring):
		QtWidgets.QTreeWidgetItem.__init__(self, qstring)
		self.program = None
		self.table = None
		self.mode = ""
		self.notelevel = 0
		self.wikipage = None
		self.wizardfile = None
		self.exmodestate = False

	def setProgram(self, program):
		""" The name of the program the tree widget is supposed to run """
		self.program = program

	def getProgram(self):
		return self.program

	def setTable(self, table):
		""" Set the table to display, as an alternative to a program """
		self.table = table

	def getTable(self):
		return self.table

	def setMode(self, mode):
		""" The mode the program is supposed to run in """
		self.mode = mode

	def getMode(self):
		return self.mode

	def setNoteLevel(self, notelevel):
		""" The note level of a program """
		self.notelevel = int(notelevel)

	def getNoteLevel(self):
		return self.notelevel

	def setWikiPage(self, page):
		""" The wiki for a program """
		self.wikipage = page

	def getWikiPage(self):
		return self.wikipage

	def setWizardFile(self, wizardfile):
		""" Set the wizard file """
		self.wizardfile = wizardfile

	def getWizardFile(self):
		return self.wizardfile

	def setExpertMode(self, state):
		"""
		Sets the expert mode state
		0 = not available
		1 = available, but turned off
		2 = available and turned on
		"""
		self.exmodestate = state

	def getExpertMode(self):
		return self.exmodestate

class PMToolButton(QtWidgets.QToolButton):
	""" Create a toggle button """
	stateChanged = QtCore.pyqtSignal(bool)

	def __init__(self):
		QtWidgets.QToolButton.__init__(self)
		self.setMinimumWidth(30)
		self.setMinimumHeight(30)

	def setDown(self, state, quiet=False):
		QtWidgets.QToolButton.setDown(self, state)
		if not quiet: self.stateChanged.emit(state)

	def mousePressEvent(self, event):
		self.setDown(not self.isDown())

	def mouseReleaseEvent(self, event):
		pass

class ProjectDialog(QtWidgets.QDialog):
	"""
	Class for the Project New and Edit dialogs
	"""
	def __init__(self, pm):
		QtWidgets.QDialog.__init__(self)
		self.pm = pm

		# Dialog line edit fields
		minwidth = 8*max(len(self.pm.pm_cwd),len(self.pm.pm_icon))
		frame = QtWidgets.QFrame()
		frame.setFrameStyle(QtWidgets.QFrame.StyledPanel)
		grid = QtWidgets.QGridLayout()
		# Add intro
		textbox = QtWidgets.QTextEdit("")
		textbox.setHtml("Welcome to the EMAN2 project manager. Please add project specific parameters below. For questions email: <a href='mailto:sludtke@bcm.edu'>sludtke@bcm.edu<\a>")
		textbox.setMaximumHeight(66)
		textbox.setReadOnly(True)
		textbox.viewport().setCursor(QtCore.Qt.ArrowCursor)

		grid.addWidget(textbox, 0, 0, 1, 2)
		# Add pm name and icon
		project_name_label = QtWidgets.QLabel("Project Name")
		self.project_name = QtWidgets.QLineEdit()
		self.project_name.setMinimumWidth(minwidth)
		grid.addWidget(project_name_label, 1, 0)
		grid.addWidget(self.project_name, 1, 1)
		icon_path_label = QtWidgets.QLabel("Project Icon")
		self.icon_path = QtWidgets.QLineEdit()
		self.icon_path.setMinimumWidth(minwidth)
		grid.addWidget(icon_path_label, 2, 0)
		grid.addWidget(self.icon_path, 2, 1)
		# Mass
		particle_mass_label = QtWidgets.QLabel("Particle Mass (kDa)")
		self.particle_mass = QtWidgets.QLineEdit()
		grid.addWidget(particle_mass_label, 3, 0)
		grid.addWidget(self.particle_mass, 3, 1)
		# Scope pars
		micrscope_cs_label = QtWidgets.QLabel("Microscope CS (mm)")
		self.micrscope_cs = QtWidgets.QLineEdit()
		microscope_voltage_label = QtWidgets.QLabel("Microscope Voltage")
		self.microscope_voltage = QtWidgets.QLineEdit()
		microscope_apix_label = QtWidgets.QLabel("Microscope apix")
		self.microscope_apix = QtWidgets.QLineEdit()
		self.invar_type = QtWidgets.QComboBox()
		self.invar_type.addItem("bispec")
		self.invar_type.addItem("harmonic")
		
		grid.addWidget(micrscope_cs_label, 4, 0)
		grid.addWidget(self.micrscope_cs, 4, 1)
		grid.addWidget(microscope_voltage_label, 5, 0)
		grid.addWidget(self.microscope_voltage, 5, 1)
		grid.addWidget(microscope_apix_label, 6, 0)
		grid.addWidget(self.microscope_apix, 6, 1)
		grid.addWidget(self.invar_type, 7, 1)

		frame.setLayout(grid)

		# Ok, cancel buttons
		done_pb = QtWidgets.QPushButton("Ok")
		cancel_pb = QtWidgets.QPushButton("Cancel")
		sgrid = QtWidgets.QGridLayout()
		sgrid.addWidget(frame,0,0,1,2)
		sgrid.addWidget(done_pb,1,0)
		sgrid.addWidget(cancel_pb,1,1)
		self.setLayout(sgrid)

		done_pb.clicked.connect(self._on_done)
		cancel_pb.clicked.connect(self._on_cancel)

		# Set values
		self.fillFields()

	def fillFields(self):
		self.project_name.setText(self.pm.pn_project_name)
		self.icon_path.setText(self.pm.pm_icon)
		self.micrscope_cs.setText(str(self.pm.getCS()))
		self.microscope_voltage.setText(str(self.pm.getVoltage()))
		self.microscope_apix.setText(str(self.pm.getAPIX()))
		self.particle_mass.setText(str(self.pm.getMass()))
		self.invar_type.setCurrentText(str(self.pm.getInvar()))

	def _on_done(self):
		self.pm.pm_projects_db['global.invariant_type'] = str(self.invar_type.currentText())
		self.pm.pm_projects_db['global.particle_mass'] = str(self.particle_mass.text())
		self.pm.pm_projects_db['global.microscope_cs'] = str(self.micrscope_cs.text())
		self.pm.pm_projects_db['global.microscope_voltage'] = str(self.microscope_voltage.text())
		self.pm.pm_projects_db['global.apix'] = str(self.microscope_apix.text())
		self.pm.pm_projects_db['project_icon'] = str(self.icon_path.text())
		self.pm.pm_projects_db['project_name'] = str(self.project_name.text())
		self.pm.statusbar.setMessage("Project: "+self.project_name.text()+" edited!","font-weight:bold;")
		self.pm.updateProject()
		self.done(0)

	def _on_cancel(self):
		self.done(1)

if __name__ == "__main__":
	import sys

	if os.path.isdir("EMAN2DB") and os.path.isfile("EMAN2DB/project.bdb") :
		print("""ERROR: This appears to be an EMAN2.0x project directory. EMAN2.1 uses a number
of different conventions. To use e2projectmanager with an older project, you will need to
first upgrade the project with e2projectupdate21.py. You can still use the e2display.py
GUI directly to browse the contents of old-style projects.""")
		sys.exit(1)

	try: sets=os.listdir("sets")
	except: sets=[]
	
	# Fix naming convention for 2.3
	bis=[i for i in sets if "_bispec" in i]
	inv=[i for i in sets if "_invar" in i]
	if len(bis)>0 and len(bis)>len(inv):
		print("""Warning: detected a pre EMAN2.3 project (_bispec files are now named _invar). The project is being updated
for the new convention in a semi-backwards compatible way. If there are any issues, you may want to run e2ctf_auto
with output_only, and regenerate any sets/""")
		for ii in bis:
			try:
				rpl=ii.replace("_bispec","_invar")
				i="sets/"+ii
				if not os.path.exists("sets/"+rpl) :
#					print("rename",i,"sets/"+rpl)
#					print("link",rpl,i)			# may fail in some cases
					os.rename(i,"sets/"+rpl)
					os.symlink(rpl,i)			# may fail in some cases
			except: pass
		
		# All _bispec particles
		bis=[("particles/"+i,i.replace("_bispec","_invar")) for i in os.listdir("particles") if "_bispec" in i]
		for i in bis:
			try:
				if not os.path.exists("particles/"+i[1]):
#					print("rename",i[0],"particles/"+i[1])
#					print("link",i[1],i[0])			# may fail in some cases
					os.rename(i[0],"particles/"+i[1])
					os.symlink(i[1],i[0])			# may fail in some cases
			except: pass

	from eman2_gui.emapplication import EMApp
	app = EMApp()
	#app = QtWidgets.QApplication(sys.argv)
	pm = EMProjectManager()
	pm.show()
	try: pm.raise_()
	except: pass
	app.exec_()
