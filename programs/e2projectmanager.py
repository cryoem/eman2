#!/usr/bin/env python
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

wikiicon = [
    '15 14 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccc',
    'ccccccccccccccc',
    'bcccbcbcbccbccb',
    'bcccbcccbccbccc',
    'bcccbcbcbcbcccb',
    'bcccbcbcbcbcccb',
    'bcbcbcbcbbccccb',
    'bcbcbcbcbbccccb',
    'bcbcbcbcbcbcccb',
    'bcbcbcbcbcbcccb',
    'bcbcbcbcbccbccb',
    'cbcbccbcbccbccb',
    'ccccccccccccccc',
    'ccccccccccccccc'
]

helpicon = [
    '15 14 2 1',
    'b c #000055',
    'c c None',
    'cccccbbbbbccccc',
    'ccccbbbbbbbcccc',
    'cccbbcccccbbccc',
    'ccbbcccccccbbcc',
    'ccbbccccccbbccc',
    'cccccccccbbcccc',
    'ccccccccbbccccc',
    'cccccccbbcccccc',
    'ccccccbbccccccc',
    'ccccccbbccccccc',
    'ccccccbbccccccc',
    'ccccccccccccccc',
    'ccccccbbccccccc',
    'ccccccbbccccccc'
]

logicon = [
    '15 14 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccc',
    'ccccccccccccccc',
    'bccccbbcccbbbcc',
    'bcccbccbcbcccbc',
    'bcccbccbcbccccc',
    'bcccbccbcbccccc',
    'bcccbccbcbccccc',
    'bcccbccbcbccbbc',
    'bcccbccbcbcccbc',
    'bcccbccbcbcccbc',
    'bbbbcbbcccbbbcc',
    'ccccccccccccccc',
    'ccccccccccccccc',
    'ccccccccccccccc'
]

noteicon = [
    '15 14 2 1',
    'b c #000055',
    'c c None',
    'ccccccccbbccccc',
    'cccccccbbbbbccc',
    'cccccccbbbbbbcc',
    'ccccccbbbbbbccc',
    'ccccccbbbbbbccc',
    'cccccbbbbbbcccc',
    'ccccbbbbbbccccc',
    'ccccbbbbbbccccc',
    'cccbbbbbbcccccc',
    'ccbbbbbbccccccc',
    'ccbccbbbccccccc',
    'ccbcccbcccccccc',
    'cccbbbccccccccc',
    'ccccccccccccccc'
]

taskicon = [
    '15 14 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccc',
    'ccccccccccccccc',
    'ccccccccccccccc',
    'ccccccccccccccc',
    'ccccccccccccccc',
    'bbbbbbbbbbbbbbb',
    'bbbbbbbbbcccccb',
    'bbbbbbbbbbbbbbb',
    'ccccccccccccccc',
    'ccccccccccccccc',
    'ccccccccccccccc',
    'ccccccccccccccc',
    'ccccccccccccccc',
    'ccccccccccccccc'
]

wizardicon = [
    '15 14 2 1',
    'b c #000055',
    'c c None',
    'cccccccccbbbccc',
    'ccccccccbbbbccc',
    'cccccccbbbccccc',
    'ccccccbbbbccccc',
    'cccccbbbbbccccc',
    'cccccbbbbbccccc',
    'ccccbbbbbbccccc',
    'ccccbbbbbbbcccc',
    'cccbbbbbbbbcccc',
    'cccbbcccbbbcccc',
    'ccbbbcbcbbbbccc',
    'bbbbbcccbbbbbbc',
    'cbbbbbbbbbbbbcc',
    'cccbbbbbbbbcccc'
]

experticon = [
    '15 14 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccc',
    'bbbbbcbbccccbbc',
    'bbbbbcbbccccbbc',
    'bbcccccbbccbbcc',
    'bbcccccbbccbbcc',
    'bbbbccccbbbbccc',
    'bbbbcccccbbcccc',
    'bbcccccccbbcccc',
    'bbccccccbbbbccc',
    'bbcccccbbccbbcc',
    'bbcccccbbccbbcc',
    'bbbbbcbbccccbbcc',
    'bbbbbcbbccccbbc',
    'ccccccccccccccc'
]

boldicon = [
    '15 14 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccc',
    'ccbbbbbbbbbcccc',
    'ccbbbbbbbbbbccc',
    'ccbbccccccbbccc',
    'ccbbccccccbbccc',
    'ccbbccccccbbccc',
    'ccbbbbbbbbbcccc',
    'ccbbbbbbbbbcccc',
    'ccbbccccccbbccc',
    'ccbbccccccbbccc',
    'ccbbccccccbbccc',
    'ccbbbbbbbbbbccc',
    'ccbbbbbbbbbcccc',
    'ccccccccccccccc'
]

italicicon = [
    '15 14 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccc',
    'cccbbbbbbbbbbbc',
    'cccbbbbbbbbbbbc',
    'cccccccbbbccccc',
    'cccccccbbbccccc',
    'ccccccbbbcccccc',
    'ccccccbbbcccccc',
    'ccccccbbbcccccc',
    'cccccbbbccccccc',
    'cccccbbbccccccc',
    'cccccbbbccccccc',
    'cbbbbbbbbbbbccc',
    'cbbbbbbbbbbbccc',
    'ccccccccccccccc'
]

underlineicon = [
    '15 14 2 1',
    'b c #000055',
    'c c None',
    'ccccccccccccccc',
    'ccbbcccccccbbcc',
    'ccbbcccccccbbcc',
    'ccbbcccccccbbcc',
    'ccbbcccccccbbcc',
    'ccbbcccccccbbcc',
    'ccbbcccccccbbcc',
    'ccbbcccccccbbcc',
    'ccbbcccccccbbcc',
    'cccbbbbbbbbbccc',
    'ccccbbbbbbbcccc',
    'ccccccccccccccc',
    'bbbbbbbbbbbbbbb',
    'bbbbbbbbbbbbbbb'
]

pmicon = [
    '15 14 2 1',
    'b c #000055',
    'c c None',
    'cccccccbccccccc',
    'cccccbbbbbccccc',
    'cccccbbbbbccccc',
    'cccccbbbbbccccc',
    'cccccbbbbbccccc',
    'cccccbbbbbbbccc',
    'cccccbbbbbbbccc',
    'cccccbbbbbccccc',
    'cccccbbbbbccccc',
    'cccccbbbbbccccc',
    'ccccbbbcbbbcccc',
    'cccbbbcccbbbccc',
    'cccbbbcccbbbccc',
    'cccbbbbbbbbbccc'
]

from EMAN2 import *
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import Qt
import os, json, re, glob, signal
import subprocess
from EMAN2db import db_open_dict
from empmwidgets import *
from valslider import EMQTColorWidget
from embrowser import EMBrowserWidget

class EMProjectManager(QtGui.QMainWindow):
	""" The EM Project Manager is a QT application to provide a GUI for EMAN2 job managment. 
	See the wiki for more details """
	def __init__(self):
		QtGui.QMainWindow.__init__(self)
		# default PM attributes
		self.pm_cwd = os.getcwd()
		self.pm_icon = os.getenv("EMAN2DIR")+"/images/EMAN2Icon.png"
		self.setWindowIcon(QtGui.QIcon(QtGui.QPixmap(pmicon)))
		
		# Set Defaults
		self.usingEMEN = False
		self.expertmode = False
		self.notebook = None
		self.taskmanager = None
		self.wikipage = None
		
		# Load the project DataBase
		self.loadPMdb()
		# Load icons
		self._load_icons()
		
		# Make the project manager layout
		self.makeMenues()
		font = QtGui.QFont()
		font.setBold(True)
		centralwidget = QtGui.QWidget()
		vsplitter = QtGui.QSplitter(QtCore.Qt.Vertical)
		
		# Make the tiltebars
		grid = QtGui.QGridLayout()
		grid.addWidget(self.makeTilteBarWidget(), 0, 0, 1, 3)
		workflowcontrollabel = QtGui.QLabel("Workflow Control", centralwidget)
		workflowcontrollabel.setFont(font)
		grid.addWidget(workflowcontrollabel, 1,0)
		guilabel = QtGui.QLabel("EMAN2 program GUI", centralwidget)
		guilabel.setFont(font)
		grid.addWidget(guilabel, 1,1)
				
		# Make the GUI, Tree and ToolBar
		grid.addWidget(self.makeStackedWidget(),2,0,1,1)
		grid.addWidget(self.makeBrowseButtonWidget(),3,0,1,1)
		grid.addWidget(self.makeStackedGUIwidget(),2,1,1,1)
		grid.addWidget(self.makeGUIToolButtons(),2,2,1,1)
		grid.addWidget(self.makeCMDButtonsWidget(),3,1,1,1)
		
		# Make status bar
		self.statusbar = EMAN2StatusBar("Welcome to the EMAN2 Project Manager")
		centralwidget.setLayout(grid)
		# Add splitter to adjust statsus bar
		vsplitter.addWidget(centralwidget)
		vsplitter.addWidget(self.statusbar)
		vsplitter.setSizes([10,1])
		
		self.setCentralWidget(vsplitter)
		
		#Update the project are construction
		self.updateProject()
		
	def closeEvent(self, event):
		""" Upon PM close, close the taskmanager and the logbook """
		if self.notebook: self.notebook.close()
		if self.taskmanager: self.taskmanager.close()
	
	def loadPMdb(self):
		"""
		Load the PM database. This is a global database. Each user on each machine has one
		"""
		# This probably wont work on WINDOWS.....
		self.pm_projects_db = db_open_dict("bdb:"+os.environ['HOME']+"#:pm_projects")
		# Load the project if we are already in the right dir
		self.pn_project_name = "Unknown"
		for project in self.pm_projects_db.keys():
			if project == "Unknown": continue	# Don't load the default project
			if self.pm_projects_db[project]["CWD"] == self.pm_cwd:
				self.pn_project_name = project
				return
		# Default, if nothing is found
		self.pm_projects_db["Unknown"] = {"CWD":self.pm_cwd,"ICON":self.pm_icon}
	
	def _load_icons(self):
		"""
		Load icons used for the tree. Additonal icons can be added using icsons.json
		"""
		self.icons = {}
		EMAN2DIR = os.getenv("EMAN2DIR")

		jsonfile = open(os.getenv("EMAN2DIR")+'/lib/pmconfig/icons.json', 'r')
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
		self.connect(editproject, QtCore.SIGNAL('triggered()'), self._on_editproject)
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
		self.connect(filebrowser, QtCore.SIGNAL('triggered()'), self._on_browse)
		utilsmenu.addAction(filebrowser)
	
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
		np = DialogNewProject(self)
		np.exec_()
		self.activateWindow()
		
	def _on_openproject(self):
		np = DialogOpenProject(self)
		np.exec_()
		self.activateWindow()
		
	def _on_editproject(self):
		np = DialogEditProject(self)
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
		self.PMIcon.setAlignment(QtCore.Qt.AlignLeft)
		grid.addWidget(self.PMIcon,0 , 0, 2, 1)
		self.PMTitle = QtGui.QLabel("EMAN2 Project Manager ")
		self.PMTitle.setAlignment(QtCore.Qt.AlignCenter)
		self.PMProjectNameBanner = QtGui.QLabel("Project Name: "+self.pn_project_name)
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
		This is the stacked widget to manage the tree types
		"""
		self.tree_stacked_widget = QtGui.QStackedWidget()
		self.tree_stacked_widget.setMinimumWidth(300)
		self.tree_stacked_widget.addWidget(self.makeSPRTreeWidget())
		self.tree_stacked_widget.addWidget(self.makeTomoTreeWidget())
		
		return self.tree_stacked_widget
	
	def makeStackedGUIwidget(self):
		"""
		Make a stacked widget
		When a python script is called for the first time a GUI widget is made and added to the stack
		The First Widget on the stack is the blank widget
		The Second Widget on the stack is the text edit widget for displaying comand line info and help
		"""
		self.gui_stacked_widget = QtGui.QStackedWidget()
		# Set the initial height of the browser
		self.gui_stacked_widget.setMinimumHeight(250)
		self.gui_stacked_widget.setFrameShape(QtGui.QFrame.StyledPanel)
		# Blank screen widget
		self.gui_stacked_widget.addWidget(QtGui.QWidget())
		# Textbox widget
		self.gui_stacked_widget.addWidget(self.getCMDTextEdit())
		self.stackedWidgetHash = {}
		
		return self.gui_stacked_widget
	
	def getCMDTextEdit(self):
		"""
		Get the Command Text edit Widget. This is incorperated into the GUI Widget Stack 
		"""
		self.guitexteditbox = QtGui.QTextEdit("")
		self.guitexteditbox.setWordWrapMode(QtGui.QTextOption.WrapAnywhere)
		return self.guitexteditbox
		
	def makeGUIToolButtons(self):
		"""
		Get get ToolBar widget
		"""
		toolwidget = QtGui.QFrame()
		#toolwidget.setFrameShape(QtGui.QFrame.StyledPanel)
		tbox = QtGui.QVBoxLayout()
		self.wikibutton = PMToolButton()
		self.wikibutton.setIcon(QtGui.QIcon(QtGui.QPixmap(wikiicon)))
		self.wikibutton.setToolTip("Show Wiki button")
		self.helpbutton = PMToolButton()
		self.helpbutton.setIcon(QtGui.QIcon(QtGui.QPixmap(helpicon)))
		self.helpbutton.setToolTip("Help button")
		self.expertbutton = PMToolButton()
		self.expertbutton.setDown(self.expertmode, quiet=True)
		self.expertbutton.setIcon(QtGui.QIcon(QtGui.QPixmap(experticon)))
		self.expertbutton.setToolTip("ExpertMode")
		self.wizardbutton = QtGui.QToolButton()
		self.wizardbutton.setIcon(QtGui.QIcon(QtGui.QPixmap(wizardicon)))
		self.wizardbutton.setToolTip("Form Wizard")
		self.wizardbutton.setMinimumWidth(30)
		self.wizardbutton.setMinimumHeight(30)
		self.logbutton = PMToolButton()
		self.logbutton.setToolTip("Display the note book")
		self.logbutton.setIcon(QtGui.QIcon(QtGui.QPixmap(noteicon)))
		self.taskmanagerbutton = PMToolButton()
		self.taskmanagerbutton.setToolTip("Diaplay the task manager")
		self.taskmanagerbutton.setIcon(QtGui.QIcon(QtGui.QPixmap(taskicon)))
		tbox.addWidget(self.wikibutton)
		tbox.addWidget(self.helpbutton)
		tbox.addWidget(self.logbutton)
		tbox.addWidget(self.taskmanagerbutton)
		tbox.addWidget(self.wizardbutton)
		tbox.addWidget(self.expertbutton)
		tbox.setContentsMargins(0,0,0,0)
		tbox.setAlignment(QtCore.Qt.AlignTop)
		toolwidget.setLayout(tbox)
		
		QtCore.QObject.connect(self.expertbutton,QtCore.SIGNAL("stateChanged(bool)"),self._on_expertmodechanged)
		QtCore.QObject.connect(self.helpbutton,QtCore.SIGNAL("stateChanged(bool)"),self._on_helpbutton)
		QtCore.QObject.connect(self.wikibutton,QtCore.SIGNAL("stateChanged(bool)"),self._on_wikibutton)
		QtCore.QObject.connect(self.logbutton,QtCore.SIGNAL("stateChanged(bool)"),self._on_logbutton)
		QtCore.QObject.connect(self.taskmanagerbutton,QtCore.SIGNAL("stateChanged(bool)"),self._on_taskmgrbutton)
		
		return toolwidget

	def _on_expertmodechanged(self, state):
		""" 
		Change the GUI upon expert mode
		"""
		self.expertmode = state
		if self.gui_stacked_widget.currentIndex() >= 2:	# First two widgets in the stack are the blank and textbox widgets
			self._set_GUI(self.getProgram(), self.getProgramMode())
		
	def _on_helpbutton(self, state):
		"""
		Load help info for the current GUI widgetitem
		"""
		# If there is wizard help, the use it otherwise load the usage info
		if self.getProgram() == "programName": return
		if state:
			# Set the stacked widget to the help info
			self.guitexteditbox.setWordWrapMode(QtGui.QTextOption.WordWrap)
			self.guitexteditbox.setText(self.loadUsage(self.getProgram()))
			self.gui_stacked_widget.setCurrentIndex(1)
		else:
			# Set the stacked widget to the GUI
			self._set_GUI(self.getProgram(), self.getProgramMode())
	
	def _on_logbutton(self, state):
		"""Load the log book
		"""
		if state:
			self.loadNoteBook()
		else:
			self.notebook.hide()
	
	def _on_wikibutton(self, state):
		""" Load the wiki help """
		if state:
			self.loadWiki()
		else:
			pass
			
	def _on_taskmgrbutton(self, state):
		"""Load the log book
		"""
		if state:
			self.loadTaskManager()
		else:
			self.taskmanager.hide()
	
	def loadUsage(self, program):
		"""
		Read in and return the usage from an e2 program
		"""
		try:
			f = open(os.getenv("EMAN2DIR")+"/bin/"+program,"r")
		except:
			self.statusbar.setMessage("Can't open usage file '%s'"%program)
			return
		begin = False
		helpstr = ""
		bregex = re.compile('usage\s*=\s*"""')
		eregex = re.compile('"""')
		for line in f.xreadlines():
			if re.search(eregex, line) and begin:
				helpstr = helpstr + line.strip()
				break
			if re.search(bregex, line):
				begin = True
			if begin:
				helpstr = helpstr + line.strip() + " "
		f.close()
		
		return helpstr
		
	def loadNoteBook(self):
		"""
		Make logbook
		"""
		if not self.notebook:
			self.notebook = NoteBook(self)
		self.notebook.show()
			
	def loadTaskManager(self):
		"""
		Make notebook
		"""
		if not self.taskmanager:
			self.taskmanager = TaskManager(self)
		self.taskmanager.show()
	
	def loadWiki(self):
		"""
		Make wiki
		"""
		if not self.wikipage:
			self.wikipage = WikiPage(self)
		self.wikipage.update()
		
	def makeBrowseButtonWidget(self):
		# Make the browse buttonv
		browsewidget = QtGui.QFrame()
		browsewidget.setFrameShape(QtGui.QFrame.StyledPanel)
		hbox = QtGui.QHBoxLayout()
		self.browsebutton = QtGui.QPushButton('Browse')
		hbox.addWidget(self.browsebutton)
		browsewidget.setLayout(hbox)
		hbox.setContentsMargins(4,4,4,4)
		
		QtCore.QObject.connect(self.browsebutton,QtCore.SIGNAL("clicked()"),self._on_browse)
		
		return browsewidget
		
	def _on_browse(self):
		self.window = EMBrowserWidget(withmodal=False,multiselect=False)
		self.window.show()
	
	def makeCMDButtonsWidget(self):
		"""
		Get the GUI command buttons widget
		"""
		cmdwidget = QtGui.QFrame()
		cmdwidget.setFrameShape(QtGui.QFrame.StyledPanel)
		hbox = QtGui.QHBoxLayout()
		self.cancelbutton = QtGui.QPushButton("Cancel")
		self.cmdlinebutton = QtGui.QPushButton("Get CMD")
		self.launchbutton = QtGui.QPushButton("Launch")
		hbox.addWidget(self.cancelbutton)
		hbox.addWidget(self.cmdlinebutton)
		hbox.addWidget(self.launchbutton)
		hbox.setContentsMargins(4,4,4,4)
		cmdwidget.setLayout(hbox)
		
		QtCore.QObject.connect(self.cancelbutton,QtCore.SIGNAL("clicked()"),self._on_cmd_cancel)
		QtCore.QObject.connect(self.cmdlinebutton,QtCore.SIGNAL("clicked()"),self._on_cmd_getcmd)
		QtCore.QObject.connect(self.launchbutton,QtCore.SIGNAL("clicked()"),self._on_cmd_launch)
		
		return cmdwidget
	
	def _on_cmd_cancel(self):
		"""
		'cancel' the present form'. This just pulls up the blank screen and updates the project
		"""
		self.gui_stacked_widget.setCurrentIndex(0)
		self.updateProject()
		
	def _on_cmd_getcmd(self):
		""" 
		Dispaly the command generated from the GUI 
		"""
		if self.gui_stacked_widget.currentIndex() <= 1: return	# Obviously we dnot want to run get command on the blank of cmd widgets
		cmd = self.gui_stacked_widget.currentWidget().getCommand()
		if not cmd: return					# The command has some errors
		self.guitexteditbox.setWordWrapMode(QtGui.QTextOption.WrapAnywhere)
		self.guitexteditbox.setText(cmd)
		self.gui_stacked_widget.setCurrentIndex(1)

	
	def _on_cmd_launch(self):
		""" 
		Launch the command from the GUI 
		"""
		if self.gui_stacked_widget.currentIndex() == 0: return	# No cmd to run
		# Else Run the command
		cmd = None
		if self.gui_stacked_widget.currentIndex() == 1:	# Use the command in the textbox
			if self.helpbutton.isDown(): return	# Don't launch using the help info
			cmd = self.gui_stacked_widget.widget(1).toPlainText()
			
		else:						# Get the command from the GUI
			cmd = self.gui_stacked_widget.currentWidget().getCommand()
		if self.launchScript(cmd): 
			self.gui_stacked_widget.setCurrentIndex(0)
			self.updateProject()	
	
	def launchScript(self, cmd):
		"""
		Start the script running 
		"""
		if not cmd: return False	# Don't excecute a broken script
		# --ipd=-2 tells the pm log book that this job is already in the pm
		if self.notebook and self.getProgramNoteLevel() > 0:
			# Only take notes if note level is greater than 0
			self.notebook.insertNewJob(cmd,local_datetime())
			self.notebook.writeNotes()
			child = subprocess.Popen((str(cmd)+" --ppid=-2"), shell=True, cwd=self.pm_projects_db[self.pn_project_name]["CWD"])
		else:
			if self.getProgramNoteLevel() > 0:
				child = subprocess.Popen(str(cmd), shell=True, cwd=self.pm_projects_db[self.pn_project_name]["CWD"])
			else:
				child = subprocess.Popen((str(cmd)+" --ppid=-2"), shell=True, cwd=self.pm_projects_db[self.pn_project_name]["CWD"])
		self.statusbar.setMessage("Program %s Launched!!!!"%str(cmd).split()[0])
		
		return True
		
	def _add_children(self, toplevel, widgetitem):
		""" 
		recursive hlper function for loadTree
		"""
		for child in toplevel["CHILDREN"]:
			qtreewidget = PMQTreeWidgetItem(QtCore.QStringList(child["NAME"]))
			qtreewidget.setIcon(0, self.icons[child["ICON"]])
			qtreewidget.setProgram(child["PROGRAM"])
			# Optional mode for the program to run in. The default is to have no mode
			if "MODE" in child: qtreewidget.setMode(child["MODE"])
			# Option note level and note elevel > 0 means add job to notebook. Good to prevent a lot of crap from piling up!
			if "NOTELEVEL" in child: qtreewidget.setNoteLevel(child["NOTELEVEL"])
			self._add_children(child, qtreewidget)
			widgetitem.addChild(qtreewidget)
			
		
	def loadTree(self, filename, treename):
		"""
		Load a workflow tree from a JSON file
		@param filename The JSON filename
		@param treename The name of the QTreWidget
		@param treewidgetdict the dictionary containing the QtreeWidgets
		"""
		try:
			jsonfile = open(filename, 'r')
		except:
			self.statusbar.setMessage("Can't open configureation file '%s'"%filename)
			return
		data = jsonfile.read()
		data = self.json_strip_comments(data)
		tree = json.loads(data)
		jsonfile.close()
		
		QTree = QtGui.QTreeWidget()
		QTree.setHeaderLabel(treename)
		
		for toplevel in tree:
			qtreewidget = PMQTreeWidgetItem(QtCore.QStringList(toplevel["NAME"]))
			qtreewidget.setIcon(0, self.icons[toplevel["ICON"]])
			qtreewidget.setProgram(toplevel["PROGRAM"])
			# Optional mode for the program to run in. The default is to have no mode
			if "MODE" in toplevel: qtreewidget.setMode(toplevel["MODE"])
			# Option note level and note elevel > 0 means add job to notebook. Good to prevent a lot of crap from piling up!
			if "NOTELEVEL" in toplevel: qtreewidget.setNoteLevel(toplevel["NOTELEVEL"])
			self._add_children(toplevel, qtreewidget)
			QTree.addTopLevelItem(qtreewidget)
		
		QtCore.QObject.connect(QTree, QtCore.SIGNAL("itemClicked(QTreeWidgetItem*,int)"), self._tree_widget_click)
		
		return QTree
	
	def json_strip_comments(self, data):
		"""This method takes a JSON-serialized string and removes
		JavaScript-style comments. These include // and /* */"""
		r = re.compile('/\\*.*\\*/', flags=re.M|re.S)
		data = r.sub("", data)
		data = re.sub("\s//.*\n", "", data)
		return data
        
        def _tree_widget_click(self, item, col):
		# Display the progrma GUI
		if item.getProgram() != "programName":
			self._set_GUI(item.getProgram(), item.getMode())
			self.updateProject()
		else:
			self.gui_stacked_widget.setCurrentIndex(0)

	
	def _set_GUI(self, program, mode):
		""" 
		Set the current GUI widget
		"""
		
		programfile = program
		# Need a widget for each program for each mode
		if self.expertmode: 
			program = "Expert"+program+mode# There is a separate widget for the advanced mode
		else:
			program = program+mode
		# Only generate the GUI widget once.....
		if program in self.stackedWidgetHash:
			# load the widget using value from the correct DB
			self.gui_stacked_widget.widget(self.stackedWidgetHash[program]).updateWidget()
			self.gui_stacked_widget.setCurrentIndex(self.stackedWidgetHash[program])
		else:
			# OR make the widget
			self.stackedWidgetHash[program] = self.gui_stacked_widget.count()
			guioptions = self._read_e2program(programfile, mode)
			# Now actually make the widget
			widget = PMGUIWidget(guioptions, programfile, self, mode)
			self.gui_stacked_widget.addWidget(widget)
			self.gui_stacked_widget.setCurrentIndex(self.stackedWidgetHash[program])
				
	def _read_e2program(self, e2program, mode):
		"""
		This a a pseudo import function to load paser info
		"""
		parser = EMArgumentParser()
		try: 
			f = open(os.getenv("EMAN2DIR")+"/bin/"+e2program,"r")
		except:
			self.statusbar.setMessage("Can't open file '%s'"%e2program)
			return
		# Regex objects
		lineregex = re.compile("^\s*parser\.add_",flags=re.I) # eval parser.add_ lines, which are  not commented out.
		moderegex = re.compile("mode\s*=\s*[\"'].*%s[{0,1}.*]{0,1}.*[\"']"%mode,flags=re.I)	# If the program has a mode only eval lines with the right mode.
		modedefre = re.compile("%s\[([\w\.\(\)\']*)\]"%mode,flags=re.I)
		defaultre = re.compile("default\s*=\s*[^,]*")
		
		# Read line and do preprocessing(set mode defaults if desired)
		for line in f.xreadlines():
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
		
	def makeSPRTreeWidget(self):
		"""
		Make the SPR tree
		"""
		self.SPRtree = self.loadTree(os.getenv("EMAN2DIR")+'/lib/pmconfig/spr.json', 'SPR')
		return self.SPRtree
	
	def makeTomoTreeWidget(self):
		"""
		Make the tomo tree widget
		"""
		self.TOMOtree = self.loadTree(os.getenv("EMAN2DIR")+'/lib/pmconfig/tomo.json', 'Tomography')
		return self.TOMOtree
		
	def getProgram(self):
		"""
		return the current program
		"""
		return self.tree_stacked_widget.currentWidget().currentItem().getProgram()
		
	def getProgramMode(self):
		"""
		return the current mode of the current program
		"""
		return self.tree_stacked_widget.currentWidget().currentItem().getMode()
		
	def getProgramNoteLevel(self):
		"""
		return the program note level, used for noting
		"""
		return self.tree_stacked_widget.currentWidget().currentItem().getNoteLevel()
		
	def getCS(self):
		""" Return the project CS """
		try:
			return self.pm_projects_db[self.pn_project_name]['CS']
		except:
			return ""
		
	def getVoltage(self):
		""" Return the project Voltage """
		try:
			return self.pm_projects_db[self.pn_project_name]['VOLTAGE']
		except:
			return ""
		
	def getAPIX(self):
		""" Return the project Apix """
		try:
			return self.pm_projects_db[self.pn_project_name]['APIX']
		except:
			return ""
	
	def getPMCWD(self):
		""" return the CWD that the pm is working in """
		return self.pm_projects_db[self.pn_project_name]['CWD']
		
	def updateProject(self):
		"""
		Update the Project. Make sure buttons, DB, title and icon observe the state of the PM
		"""
		self.PMProjectNameBanner.setText("Project Name: "+self.pn_project_name)
		# Icon loading may be a bit slow
		if self.pm_icon != self.pm_projects_db[self.pn_project_name]['ICON']:
			self.PMIcon.setIcon(self.pm_projects_db[self.pn_project_name]['ICON'])
		self.pm_icon = self.pm_projects_db[self.pn_project_name]['ICON']
		if self.pm_cwd != self.pm_projects_db[self.pn_project_name]["CWD"]:
			self.pm_cwd = self.pm_projects_db[self.pn_project_name]["CWD"]
			os.chdir(self.pm_cwd)
			# update widget to reflect new CWD
			if self.gui_stacked_widget.currentIndex() > 1:
				self.gui_stacked_widget.currentWidget().updateWidget()
		self.logbutton.setDown(bool(self.notebook), quiet=True)
		self.taskmanagerbutton.setDown(bool(self.taskmanager), quiet=True)
		# Help button should only be down if the textbox widget is displayed
		if self.gui_stacked_widget.currentIndex() != 1:
			self.helpbutton.setDown(False, quiet=True)
	
class EMAN2StatusBar(QtGui.QTextEdit):
	"""
	The Stats bar for PM
	"""
	def __init__(self, text):
		QtGui.QTextEdit.__init__(self)
		self.setFrameShape(QtGui.QFrame.Panel | QtGui.QFrame.Sunken)
		self.setLineWidth(2)
		#self.setMargin(4)
		self.setTextInteractionFlags(QtCore.Qt.NoTextInteraction)
		self.viewport().setCursor(QtCore.Qt.ArrowCursor)
		self.setMessage(text)
		self.setToolTip("This is the Status Bar")
		
	def setMessage(self, text):
		textcursor = self.textCursor()
		textcursor.movePosition(QtGui.QTextCursor.End)
		if textcursor.columnNumber() !=0:
			self.insertPlainText("\n"+text)
		else:
			self.insertPlainText(text)
		self.verticalScrollBar().setValue(self.verticalScrollBar().maximum())
		

		
class PMIcon(QtGui.QLabel):
	"""
	The Icon manager for PM
	"""
	def __init__(self, image, parent=None):
		QtGui.QLabel.__init__(self, ("<img src=\"%s\" />")%image, parent)
	
	def setIcon(self, image):
		self.setText(("<img src=\"%s\" />")%image)

class WikiPage():
	"""
	Load the wikipage
	"""
	def __init__(self, pm):
		self.pm = weakref.ref(pm)
	
		
		self.downloadHTML('xx')
	
	def downloadHTML(self, page):
		import webbrowser
		webbrowser.open('http://blake.bcm.tmc.edu/emanwiki/EMAN2/Programs/e2proc2d')
		#self.wikibox.setHtml(wiki.read())
		
	def update(self):
		pass
	
	def _on_close(self):
		self.close()
		
class NoteBook(QtGui.QWidget):
	"""
	The Logbook for PM. The note book will reflect top levels jobs run, even if they were run on the command line
	"""
	def __init__(self, pm):
		QtGui.QWidget.__init__(self)
		self.pm = weakref.ref(pm)
		self.donotsave = False
		
		self.setWindowTitle('NoteBook')
		grid = QtGui.QGridLayout()
		font = QtGui.QFont()
		font.setBold(True)
		textlabel = QtGui.QLabel("EMAN2 NoteBook")
		textlabel.setFont(font)
		self.texteditbox = PMTextEdit(self)
		grid.addWidget(textlabel,0,0)
		grid.addWidget(self.getToolBar(),1,0,1,2)
		grid.addWidget(self.texteditbox,2,0,1,2)
		self.savepb = QtGui.QPushButton("Save")
		self.closepb = QtGui.QPushButton("Close")
		grid.addWidget(self.savepb, 3,0)
		grid.addWidget(self.closepb, 3,1)
		self.setLayout(grid)
		self.setMinimumWidth(600)
		
		# Load the logbook and update
		self.loadNoteBook()
		self.checkEMAN2LogFile()
		
		self.connect(self.savepb, QtCore.SIGNAL('clicked()'), self._on_save)
		self.connect(self.closepb, QtCore.SIGNAL('clicked()'), self._on_close)
	
	def getToolBar(self):
		""" Return the toolbar widget """
		tbwidget = QtGui.QWidget()
		hbox = QtGui.QHBoxLayout()
		self.dbdict = db_open_dict("bdb:"+self.pm().getPMCWD()+"#notebook")
		
		# font type
		self.fontdb = QtGui.QFontDatabase()
		self.fontfamily = QtGui.QComboBox()
		self.fontfamily.addItems(self.fontdb.families())	
		
		# font size
		self.fontsizecb = QtGui.QComboBox()
		
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
		self.connect(self.fontfamily, QtCore.SIGNAL("activated(int)"), self._fontfamilychange)
		self.connect(self.fontsizecb, QtCore.SIGNAL("activated(int)"), self._fontchange)
		self.connect(self.boldbutton, QtCore.SIGNAL("stateChanged(bool)"), self._fontchange)
		self.connect(self.italicbutton, QtCore.SIGNAL("stateChanged(bool)"), self._fontchange)
		self.connect(self.underlinebutton, QtCore.SIGNAL("stateChanged(bool)"), self._fontchange)
		self.connect(self.fontcolor, QtCore.SIGNAL("newcolor(QColor)"), self._fontchange)
		
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
			self.pm().statusbar.setMessage("No pmnotes file found, starting a new one...")
			self.texteditbox.insertHtml("<b>Project Name: %s<br>EMAN Version: %s<br>Path: %s<br></b>"%(self.pm().pn_project_name,EMANVERSION,self.pm().pm_projects_db[self.pm().pn_project_name]["CWD"]))
			
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
		
	def closeEvent(self, event):
		if not self.donotsave: self.writeNotes()
		self.pm().notebook = None
		self.pm().updateProject()

class PMTextEdit(QtGui.QTextEdit):
	""" Sub class of the QTextEdit widget to observe the PMNoteBook """
	def __init__(self, parent):
		QtGui.QTextEdit.__init__(self)
		self.parent = weakref.ref(parent)
		
		self.setPMFontWeight(QtGui.QFont.Normal)
		self.setPMTextColor(QtGui.QColor(0,0,0))
		
	def mousePressEvent(self,event):
		""" Stupidly, TextEdit resets the font based on its context, which I find undesireable """
		QtGui.QTextEdit.mousePressEvent(self, event)
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
		self.setPMFontSize(int(self.parent().fontsizecb.currentText()))
		self.setPMTextColor(self.parent().fontcolor.getColor())
		self.setPMFontFamily(self.parent().fontfamily.currentText())
		self.setPMFontItalic(self.parent().italicbutton.isDown())
		self.setPMFontUnderline(self.parent().underlinebutton.isDown())
		
class TaskManager(QtGui.QWidget):
	"""
	The Logbook for PM
	"""
	def __init__(self, pm):
		QtGui.QWidget.__init__(self)
		self.pm = weakref.ref(pm)
		self.setWindowTitle('Task Manager')
		grid = QtGui.QGridLayout()
		font = QtGui.QFont()
		font.setBold(True)
		textlabel = QtGui.QLabel("Running Tasks")
		textlabel.setFont(font)
		grid.addWidget(textlabel,0,0)
		self.list_widget = QtGui.QListWidget()
		grid.addWidget(self.list_widget,1,0,1,2)
		self.killpb = QtGui.QPushButton("Kill")
		self.closepb = QtGui.QPushButton("Close")
		grid.addWidget(self.killpb, 2,0)
		grid.addWidget(self.closepb, 2,1)
		self.setLayout(grid)
		self.resize(400, 200)
		
		
		self.connect(self.closepb, QtCore.SIGNAL('clicked()'), self._on_close)
		self.connect(self.killpb, QtCore.SIGNAL('clicked()'), self._on_kill)
		self.tasks=None
		self.update_tasks()
		
		# A timer for updates
		self.timer = QtCore.QTimer(self);
 		QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.update_tasks)
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
			
			# This contains (file location,process id, status, command name)
			ppid = None
			if '/' in tsk[2]: 
				ids = tsk[2].split("/")
				pid=int(ids[0])
				if len(ids) == 2:
					ppid = int(ids[1].strip())
			else : pid=int(tsk[2])
			if os.name=="nt": tsk[4]=tsk[4].convert("\\","/")
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
					# see if any jobs have crashed, if not then do nothing
					self.tasks=[i for i in self.tasks if self.check_task(self.prevfin,i)]
					self.update_list()
					return
				
			except:
				print "Error, couldn't stat(.eman2log.txt)."
				return
			
			# Load in the file
			try: fin=open(self.pm().getPMCWD()+"/.eman2log.txt","r")
			except: return
			self.prevfin = fin
			
			# Go through existing entries and see if they've finished
			self.tasks=[i for i in self.tasks if self.check_task(fin,i)]
			
			# now check for new entries
			fin.seek(self.lastpos)
			self.parse_new_commands(fin)
		# Update the list
		self.update_list()

	def update_list(self):
		#self.list_widget.clear()
		listitems = self.getListItems()
		for t in self.tasks:
			if t[1] in listitems.keys():
				del(listitems[t[1]])
				continue
			listwigetitem = PMQListWidgetItem("%s  %s(%s)  %s"%(t[2][5:16],t[3][:4],t[1],t[4]))
			listwigetitem.setPID(t[1])
			listwigetitem.setPPID(t[5])
			listwigetitem.setProgramName(t[4])
			self.list_widget.addItem(listwigetitem)
		# Remove items that have stopped running
		for item in listitems.values():
			self.list_widget.takeItem(self.list_widget.row(item))
	
	def getListItems(self):
		itemdict = {}
		for i in xrange(self.list_widget.count()):
			itemdict[self.list_widget.item(i).getPID()] = self.list_widget.item(i)
		return itemdict
			
	def reset(self):
		self.tasks = None
		self.update_tasks()
		
	def _on_close(self):
		self.close()	
	
	def _on_kill(self):
		selitems = self.list_widget.selectedItems()
		self.update_tasks()
		for item in selitems:
			if os.name == 'posix': # Kill by this method will only work for unix/linux type OS
				# The PID Mafia occurs below: # Kill all children and siblings
				for task in self.tasks:
					if item.getPPID() == task[5]:
						print "killing self process", task[1]
						os.kill(task[1],signal.SIGTERM)
						self._recursivekill(task[1])
				# kill parent (top level item)
				if item.getPPID() > 0:
					print "KIlling parent"
					os.kill(item.getPPID(),signal.SIGTERM)
			else:
				# Windows kill
				pass
	
	def _recursivekill(self, pid):
		for task in self.tasks:
			if pid == task[5]:
				print "Killing child process", task[1]
				os.kill(task[1],signal.SIGTERM)
				self._recursivekill(task[1])
		
	def closeEvent(self, event):
		self.pm().taskmanager = None
		self.pm().updateProject()

class PMQListWidgetItem(QtGui.QListWidgetItem):
	""" A subclass of the QListWidgetItem for use in the task manger """
	def __init__(self, text):
		QtGui.QListWidgetItem.__init__(self, text)
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
		
class PMGUIWidget(QtGui.QScrollArea):
	"""
	Creates a GUI widget using a dict derived from the e2program options
	"""
	def __init__(self, options, program, pm, mode):
		QtGui.QScrollArea.__init__(self)
		self.widgetlist = []
		self.cwd  = pm.pm_projects_db[pm.pn_project_name]['CWD']	# The working directory that were in
		self.program = program
		self.db = db_open_dict("bdb:"+str(self.cwd)+"#"+program)
		self.pm = weakref.ref(pm)
		self.mode = mode
		
		gridbox = QtGui.QGridLayout()
		for option in options:
			"""create the correct widget type"""
			if ('expert' in  option) and not self.pm().expertmode: continue	# Do no add expertmode optionsif not in expert mode

			if option['guitype'] == 'header': 
				widget = PMHeaderWidget(option['name'], option['title'])
			if option['guitype'] == 'filebox':
				widget = PMFileNameWidget(option['name'], self.getDefault(option), self.getSharingMode(option), self.getBrowser(option), postional=self.getPositional(option), initdefault=self.getDefault(option, nodb=True),checkfileexist=self.getFileCheck(option))
				fileboxwidget = widget
			if option['guitype'] == 'dirbox':
				widget = PMDirectoryWidget(option['name'], option['dirbasename'], self.getDefault(option), self.getSharingMode(option), postional=self.getPositional(option), initdefault=self.getDefault(option, nodb=True))
			if option['guitype'] == 'symbox':
				widget = PMSymWidget(option['name'], self.getDefault(option), self.getSharingMode(option), initdefault=self.getDefault(option, nodb=True))
			if option['guitype'] == 'multisymbox':
				widget = PMMultiSymWidget(option['name'], self.getSharingMode(option), initdefault=self.getDefault(option, nodb=True))
				self.connect(fileboxwidget,QtCore.SIGNAL("pmfilename(QString)"),widget.update)
				widget.update(fileboxwidget.getValue())
				widget.setValue(self.getDefault(option))
			if option['guitype'] == 'intbox':
				widget = PMIntEntryWidget(option['name'], self.getDefault(option), self.getSharingMode(option), self.getLRange(option), self.getURange(option), postional=self.getPositional(option), initdefault=self.getDefault(option, nodb=True))
			if option['guitype'] == 'floatbox':
				widget = PMFloatEntryWidget(option['name'], self.getDefault(option), self.getSharingMode(option), self.getLRange(option), self.getURange(option), postional=self.getPositional(option), initdefault=self.getDefault(option, nodb=True))
			if option['guitype'] == 'boolbox':
				widget = PMBoolWidget(option['name'], self.getDefault(option), self.getSharingMode(option), initdefault=self.getDefault(option, nodb=True))
			if option['guitype'] == 'strbox':
				widget = PMStringEntryWidget(option['name'], self.getDefault(option), self.getSharingMode(option), postional=self.getPositional(option), initdefault=self.getDefault(option, nodb=True))
			if option['guitype'] == 'comboparambox':
				widget = PMComboParamsWidget(option['name'], self.getChoices(option), self.getDefault(option), self.getSharingMode(option), postional=self.getPositional(option), initdefault=self.getDefault(option, nodb=True))
			if option['guitype'] == 'combobox':
				widget = PMComboWidget(option['name'], self.getChoices(option), self.getDefault(option), self.getSharingMode(option), postional=self.getPositional(option), initdefault=self.getDefault(option, nodb=True))
			if option['guitype'] == 'automask3d':	
				widget = PMAutoMask3DWidget(option['name'], self.getDefault(option), self.getSharingMode(option), initdefault=self.getDefault(option, nodb=True))
			
			# Setup each widget
			self.connect(widget,QtCore.SIGNAL("pmmessage(QString)"),self._on_message)
			widget.setToolTip(option['help'])
			self.widgetlist.append(widget)
			gridbox.addWidget(widget, option['row'], option['col'], self.getRowSpan(option), self.getColSpan(option))
		
		# Now make a widget and add it to the scroll area
		gridbox.setContentsMargins(0,0,0,0)
		self.scwidget = QtGui.QWidget()
		self.scwidget.setLayout(gridbox)
		self.scwidget.setMinimumWidth(self.width()-1.5*self.verticalScrollBar().width())
		self.scwidget.setMaximumWidth(self.width()-1.5*self.verticalScrollBar().width())
		self.setWidget(self.scwidget)
		self.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
		self.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOn)
		self.setSizePolicy(QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed)
		
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
		""" return the default value according to the folowing rules"""
		# If there is a DB and its usage is desired the default will be the DB value
		if not nodb and self.db[option['name']+self.getSharingMode(option)]: return self.db[option['name']+self.getSharingMode(option)]	# Return the default if it exists in the DB
		default = ""
		if 'default' in option: default = option['default']
		if type(default) == str and "self.pm()" in default: default = eval(default)	# eval CS, apix, voltage, etc, a bit of a HACK, but it works
		return default
		
	def getLRange(self, option):
		""" Return the Lower range if it exists"""
		lrange = None
		if 'lrange' in option: lrange = option['lrange']
		return lrange
		
	def getURange(self, option):
		""" Return the Upper range if it exists"""
		urange = None
		if 'urange' in option: urange = option['urange']
		return urange
	
	def getChoices(self, option):
		choices = []
		if 'choicelist' in option:
			choices = eval(option['choicelist'])
			# If it is a dict, get the keys
			if type(choices) == type({}):
				choices = choices.keys()
		return choices
	
	def getPositional(self, option):
		# See if the arugment is positional or not
		positional = False	# Defult if not provided
		if 'positional' in option: positional = option['positional']
		return positional
	
	def getFileCheck(self, option):
		filecheck = True
		if 'filecheck' in option: filecheck = option['filecheck']
		return filecheck
	
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
		
	def updateWidget(self):
		# reload the DB if necessary (when projects are changed)
		thiscwd = self.pm().pm_projects_db[self.pm().pn_project_name]['CWD']
		if thiscwd != self.cwd:
			self.cwd = thiscwd
			self.db = db_open_dict("bdb:"+str(self.cwd)+"#"+self.program)
		# Might enable serval dbs to be loaded, but we will implment this later
		for widget in self.widgetlist:
			# If this is not a value holding widget continue
			if widget.getArgument() == None: continue
			if self.db[widget.getName()+widget.getMode()] == None:
				widget.setValue(widget.initdefault)
			else:
				widget.setValue(self.db[widget.getName()+widget.getMode()])
		
	def getCommand(self):
		#Loop and check for errors and set the DB
		args = self.program
		for widget in self.widgetlist:
			# If this is not a value holding widget continue
			if widget.getArgument() == None: continue
			# Check for errors before we launch script
			if widget.getErrorMessage():
				self.pm().statusbar.setMessage(widget.getErrorMessage())
				return None
			# Save the value
			self.db[widget.getName()+widget.getMode()] = widget.getValue()
			args += " "+widget.getArgument()
			
		self.pm().statusbar.setMessage("")	# Blank Status bar
		return args
			
	def _on_message(self, QString):
		self.pm().statusbar.setMessage(str(QString))
		
class PMQTreeWidgetItem(QtGui.QTreeWidgetItem):
	"""
	Custon QTreeWidget for PM
	"""
	def __init__(self, qstring):
		QtGui.QTreeWidgetItem.__init__(self, qstring)
		self.program = None
		self.mode = ""
		self.notelevel = 0
		
	def setProgram(self, program):
		""" The name of the program the tree widget is supposed to run """
		self.program = program
		
	def getProgram(self):
		return self.program
		
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

class PMToolButton(QtGui.QToolButton):
	""" Create a toogle button """
	def __init__(self):
		QtGui.QToolButton.__init__(self)
		self.setMinimumWidth(30)
		self.setMinimumHeight(30)
		
	def setDown(self, state, quiet=False):
		QtGui.QToolButton.setDown(self, state)
		if not quiet: self.emit(QtCore.SIGNAL("stateChanged(bool)"), state)
		
	def mousePressEvent(self, event):
		self.setDown(not self.isDown())
	
	def mouseReleaseEvent(self, event):
		pass

class ProjectDialogBase(QtGui.QDialog):
	"""
	Base class for the Project New and Edit dialogs
	"""
	def __init__(self, pm):
		QtGui.QDialog.__init__(self)
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
		
		# Scope pars
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
		
		frame.setLayout(grid)
		
		# Ok, cancel buttons
		done_pb = QtGui.QPushButton("Ok")
		cancel_pb = QtGui.QPushButton("Cancel")
		sgrid = QtGui.QGridLayout()
		sgrid.addWidget(frame,0,0,1,2)
		sgrid.addWidget(done_pb,1,0)
		sgrid.addWidget(cancel_pb,1,1)
		self.setLayout(sgrid)
		
		self.connect(done_pb, QtCore.SIGNAL('clicked()'), self._on_done)
		self.connect(cancel_pb, QtCore.SIGNAL('clicked()'), self._on_cancel)
		
	def _on_cancel(self):
		self.done(1)
		
class DialogNewProject(ProjectDialogBase):
	"""
	Generate the new projects dialog
	"""
	def __init__(self, pm):
		ProjectDialogBase.__init__(self, pm)
		self.setWindowTitle('New Project')
		
		
	def _on_done(self):
		for project in self.pm.pm_projects_db.keys():
			if project == "Unknown": continue
			if project == self.project_name.text():
				self.pm.statusbar.setMessage("Project Name is alreay in use!!!")
				return
			if self.pm.pm_projects_db[project]["CWD"] == self.project_directory.text():
				self.pm.statusbar.setMessage("This project directory is already in use by another project!!!")
				return
		self.pm.pm_projects_db[str(self.project_name.text())] = {"CWD":str(self.project_directory.text()),"ICON":str(self.icon_path.text()),"CS":str(self.micrscope_cs.text()),"VOLTAGE":str(self.microscope_voltage.text()),"APIX":str(self.microscope_apix.text())}
		self.pm.pn_project_name = str(self.project_name.text())
		
		self.pm.statusbar.setMessage("Project: "+self.project_name.text()+" created :)")
		self.pm.updateProject()
		self.done(0)
	
class DialogEditProject(ProjectDialogBase):
	"""
	Edit an existing project
	"""
	def __init__(self, pm):
		ProjectDialogBase.__init__(self, pm)
		self.setWindowTitle('Edit Project')
		self.project_name.setText(self.pm.pn_project_name)
		self.project_name.setReadOnly(True)
		self.project_directory.setText(self.pm.pm_cwd)
		
		# CS, VOLTAGE, APIX
		db = self.pm.pm_projects_db[self.pm.pn_project_name]
		self.micrscope_cs.setText(db['CS'])
		self.microscope_voltage.setText(db['VOLTAGE'])
		self.microscope_apix.setText(db['APIX'])
		
	def _on_done(self):
		self.pm.pm_projects_db[self.pm.pn_project_name] = {"CWD":str(self.project_directory.text()),"ICON":str(self.icon_path.text()),"CS":str(self.micrscope_cs.text()),"VOLTAGE":str(self.microscope_voltage.text()),"APIX":str(self.microscope_apix.text())}
		self.pm.statusbar.setMessage("Project: "+self.project_name.text()+" edited!")
		self.pm.updateProject()
		self.done(0)
	
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
			if project == "Unknown": continue
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
		self.pm.pn_project_name = str(self.list_widget.currentItem().text())
		self.pm.statusbar.setMessage("Opened Project: "+self.pm.pn_project_name)
		self.pm.updateProject()
		if self.pm.gui_stacked_widget.currentIndex() >= 2:
			self.pm.gui_stacked_widget.currentWidget().updateWidget()
		else:
			self.pm.gui_stacked_widget.setCurrentIndex(0)
		self.done(0)
		
	def _on_cancel(self):
		self.done(1)
		
if __name__ == "__main__":
	import sys
	from emapplication import EMApp
	app = EMApp()
	#app = QtGui.QApplication(sys.argv)
	pm = EMProjectManager()
	pm.show()
	app.exec_()