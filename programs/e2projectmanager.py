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

from EMAN2 import *
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import Qt
import os, json, re, glob
import subprocess
from EMAN2db import db_open_dict
from pmwidgets import *

class EMProjectManager(QtGui.QMainWindow):
	def __init__(self):
		QtGui.QMainWindow.__init__(self)
		# default PM attributes
		self.pm_cwd = os.getcwd()
		self.pm_icon = os.getenv("EMAN2DIR")+"/images/EMAN2Icon.png"
		
		# Set Defaults
		self.usingEMEN = False
		self.expertmode = False
		self.logbook = None
		self.taskmanager = None
		
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
		grid.addWidget(self.makeTilteBarWidget(), 0, 0, 1, 3)
		workflowcontrollabel = QtGui.QLabel("Workflow Control", centralwidget)
		workflowcontrollabel.setFont(font)
		grid.addWidget(workflowcontrollabel, 1,0)
		guilabel = QtGui.QLabel("EMAN2 program GUI", centralwidget)
		guilabel.setFont(font)
		grid.addWidget(guilabel, 1,1)
		grid.addWidget(self.makeStackedWidget(),2,0,2,1)
		grid.addWidget(self.makeStackedGUIwidget(),2,1,1,1)
		grid.addWidget(self.makeGUIToolButtons(),2,2,1,1)
		grid.addWidget(self.makeCMDButtonsWidget(),3,1,1,1)

		# Make status bar
		self.statusbar = EMAN2StatusBar("Welcome to the EMAN2 Project Manager")
		grid.addWidget(self.statusbar, 4,0,1,3)
		centralwidget.setLayout(grid)
		# Make status bar
		self.setCentralWidget(centralwidget)
		
		#Update the project are construction
		self.updateProject()
		
	def closeEvent(self, event):
		if self.logbook: self.logbook.close()
	
	def loadPMdb(self):
		"""
		Load the PM database
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
		Load icons used for the tree 
		"""
		self.icons = {}
		EMAN2DIR = os.getenv("EMAN2DIR")
		self.icons["single_image"] = QtGui.QIcon(EMAN2DIR+"/images/macimages/single_image.png")
		self.icons["multiple_images"] = QtGui.QIcon(EMAN2DIR+"/images/macimages/multiple_images.png")
		self.icons["green_boxes"] = QtGui.QIcon(EMAN2DIR+"/images/macimages/green_boxes.png")
		self.icons["ctf"] = QtGui.QIcon(EMAN2DIR+"/images/macimages/ctf.png")
		self.icons["web"] = QtGui.QIcon(EMAN2DIR+"/images/macimages/classes.png")
		self.icons["single_image_3d"] = QtGui.QIcon(EMAN2DIR+"/images/macimages/single_image_3d.png")
		self.icons["refine"] = QtGui.QIcon(EMAN2DIR+"/images/macimages/refine.png")
		self.icons["eulers"] = QtGui.QIcon(EMAN2DIR+"/images/macimages/eulerxplor.png")
		self.icons["resolution"] = QtGui.QIcon(EMAN2DIR+"/images/macimages/plot.png")
		self.icons["tomo_hunter"] = QtGui.QIcon(EMAN2DIR+"/images/macimages/tomo_hunter.png")
		
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
	
	def makeStackedWidget(self):
		"""
		This is the stacked widget to manage the tree types
		"""
		self.tree_stacked_widget = QtGui.QStackedWidget()
		self.tree_stacked_widget.setMinimumWidth(300)
		self.tree_stacked_widget.setMaximumWidth(300)
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
		self.gui_stacked_widget.setFrameShape(QtGui.QFrame.StyledPanel)
		self.gui_stacked_widget.addWidget(QtGui.QWidget())
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
		self.logbutton = PMToolButton()
		self.logbutton.setToolTip("Display the log book")
		self.logbutton.setIcon(QtGui.QIcon(QtGui.QPixmap(logicon)))
		self.taskmanagerbutton = PMToolButton()
		self.taskmanagerbutton.setToolTip("Diaplay the task manager")
		self.taskmanagerbutton.setIcon(QtGui.QIcon(QtGui.QPixmap(taskicon)))
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
			self.loadLogBook()
		else:
			self.logbook.close()
	
	def _on_taskmgrbutton(self, state):
		"""Load the log book
		"""
		if state:
			self.loadTaskManager()
		else:
			self.taskmanager.close()
	
	def loadUsage(self, program):
		"""
		Read in and return the usage from an e2 program
		"""
		f = open(os.getenv("EMAN2DIR")+"/bin/"+program,"r")
		begin = False
		helpstr = ""
		for line in f.xreadlines():
			if '"""' in line and begin:
				helpstr = helpstr + line.strip() + " "
				break
			if 'usage = """' in line:
				begin = True
			if begin:
				helpstr = helpstr + line.strip() + " "
		f.close()
		
		return helpstr
		
	def loadLogBook(self):
		"""
		Make logbook
		"""
		if not self.logbook:
			self.logbook = LogBook(self)
			self.logbook.show()
			
	def loadTaskManager(self):
		"""
		Make logbook
		"""
		if not self.taskmanager:
			self.taskmanager = TaskManager(self)
			self.taskmanager.show()		
		
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
		child = subprocess.Popen(str(cmd), shell=True, cwd=self.pm_projects_db[self.pn_project_name]["CWD"])
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
			self._add_children(child, qtreewidget)
			widgetitem.addChild(qtreewidget)
			
		
	def loadTree(self, filename, treename):
		"""
		Load a workflow tree from a JSON file
		@param filename The JSON filename
		@param treename The name of the QTreWidget
		@param treewidgetdict the dictionary containing the QtreeWidgets
		"""
		jsonfile = open(filename, 'r')
		data = jsonfile.read()
		data = self.json_strip_comments(data)
		tree = json.loads(data)
		jsonfile.close()
		
		QTree = QtGui.QTreeWidget()
		QTree.setMinimumHeight(400)
		QTree.setHeaderLabel(treename)
		
		for toplevel in tree:
			qtreewidget = PMQTreeWidgetItem(QtCore.QStringList(toplevel["NAME"]))
			qtreewidget.setIcon(0, self.icons[toplevel["ICON"]])
			qtreewidget.setProgram(toplevel["PROGRAM"])
			# Optional mode for the program to run in. The default is to have no mode
			if "MODE" in toplevel: qtreewidget.setMode(toplevel["MODE"])
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
			self.statusbar.setMessage("GUI '%s' displayed"%str(item.getProgram()))	# Blank Status bar
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
			program = "Expert"+program+mode	# There is a separate widget for the advanced mode
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
			widget = PMGUIWidget(guioptions, programfile, self)
			self.gui_stacked_widget.addWidget(widget)
			self.gui_stacked_widget.setCurrentIndex(self.stackedWidgetHash[program])
				
	def _read_e2program(self, e2program, mode):
		"""
		This a a pseudo import function to load paser info
		"""
		parser = EMArgumentParser()
		f = open(os.getenv("EMAN2DIR")+"/bin/"+e2program,"r")
		lineregex = re.compile("^\s*parser\.add_",flags=re.I) # eval parser.add_ lines, which are  not commented out.
		moderegex = re.compile("mode\s*=\s*[\"'].*%s.*[\"']"%mode,flags=re.I)	# If the program has a mode only eval lines with the right mode.
		for line in f.xreadlines():
			if mode and not re.search(moderegex, line): continue	# If we are running the program in a mode, then only eval mode lines
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
		self.logbutton.setDown(bool(self.logbook), quiet=True)
		self.taskmanagerbutton.setDown(bool(self.taskmanager), quiet=True)
		# Help button should only be down if the textbox widget is displayed
		if self.gui_stacked_widget.currentIndex() != 1:
			self.helpbutton.setDown(False, quiet=True)
		
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
	def __init__(self, pm):
		QtGui.QWidget.__init__(self)
		self.pm = weakref.ref(pm)
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
		
		self.connect(self.closepb, QtCore.SIGNAL('clicked()'), self._on_close)
		
	def _on_close(self):
		self.close()
		
	def closeEvent(self, event):
		self.pm().logbook = None
		self.pm().updateProject()

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
		textlabel = QtGui.QLabel("Task Manager")
		textlabel.setFont(font)
		grid.addWidget(textlabel,0,0)
		self.killpb = QtGui.QPushButton("Kill")
		self.closepb = QtGui.QPushButton("Close")
		grid.addWidget(self.killpb, 1,0)
		grid.addWidget(self.closepb, 1,1)
		self.setLayout(grid)
			
		self.connect(self.closepb, QtCore.SIGNAL('clicked()'), self._on_close)
		
	def _on_close(self):
		self.close()	
		
	def closeEvent(self, event):
		self.pm().taskmanager = None
		self.pm().updateProject()
		
class PMGUIWidget(QtGui.QScrollArea):
	"""
	Creates a GUI widget using a dict derived from the e2program options
	"""
	def __init__(self, options, program, pm):
		QtGui.QScrollArea.__init__(self)
		self.widgetlist = []
		self.cwd  = pm.pm_projects_db[pm.pn_project_name]['CWD']	# The working directory that were in
		self.program = program
		self.db = db_open_dict("bdb:"+str(self.cwd)+"#"+program)
		self.pm = weakref.ref(pm)
		
		gridbox = QtGui.QGridLayout()
		for option in options:
			"""create the correct widget type"""
			if ('expert' in  option) and not self.pm().expertmode: continue	# Do no add expertmode optionsif not in expert mode

			if option['guitype'] == 'header': 
				widget = PMHeaderWidget(option['name'], option['title'])
			if option['guitype'] == 'filebox':
				widget = PMFileNameWidget(option['name'], self.getDefault(option), postional=self.getPositional(option), initdefault=self.getDefault(option, nodb=True),checkfileexist=self.getFileCheck(option))
				fileboxwidget = widget
			if option['guitype'] == 'symbox':
				widget = PMSymWidget(option['name'], self.getDefault(option), initdefault=self.getDefault(option, nodb=True))
			if option['guitype'] == 'multisymbox':
				widget = PMMultiSymWidget(option['name'], initdefault=self.getDefault(option, nodb=True))
				self.connect(fileboxwidget,QtCore.SIGNAL("pmfilename(QString)"),widget.update)
				widget.update(fileboxwidget.getValue())
				widget.setValue(self.getDefault(option))
			if option['guitype'] == 'intbox':
				widget = PMIntEntryWidget(option['name'], self.getDefault(option), self.getLRange(option), self.getURange(option), postional=self.getPositional(option), initdefault=self.getDefault(option, nodb=True))
			if option['guitype'] == 'floatbox':
				widget = PMFloatEntryWidget(option['name'], self.getDefault(option), self.getLRange(option), self.getURange(option), postional=self.getPositional(option), initdefault=self.getDefault(option, nodb=True))
			if option['guitype'] == 'boolbox':
				widget = PMBoolWidget(option['name'], self.getDefault(option), initdefault=self.getDefault(option, nodb=True))
			if option['guitype'] == 'strbox':
				widget = PMStringEntryWidget(option['name'], self.getDefault(option), postional=self.getPositional(option), initdefault=self.getDefault(option, nodb=True))
			if option['guitype'] == 'comboparambox':
				widget = PMComboParamsWidget(option['name'], self.getChoices(option), self.getDefault(option), postional=self.getPositional(option), initdefault=self.getDefault(option, nodb=True))
			if option['guitype'] == 'combobox':
				widget = PMComboWidget(option['name'], self.getChoices(option), self.getDefault(option), postional=self.getPositional(option), initdefault=self.getDefault(option, nodb=True))
			if option['guitype'] == 'automask3d':	
				widget = PMAutoMask3DWidget(option['name'], self.getDefault(option), initdefault=self.getDefault(option, nodb=True))
			
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
		default = ""
		# Set default to GUI default if available otherwise set to default if available
		if 'guidefault' in option:
			default = option['guidefault']
		else:
			if 'default' in option: default = option['default']
		# If there is no DataBase or it isn't desired return
		if nodb: return default
		# If there is a DB and its usage is desired the default will be the DB value
		if self.db[option['name']]: return self.db[option['name']]	# Return the default if it exists in the DB
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
			if self.db[widget.getName()] == None:
				widget.setValue(widget.initdefault)
			else:
				widget.setValue(self.db[widget.getName()])
		
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
			self.db[widget.getName()] = widget.getValue()
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

class PMToolButton(QtGui.QToolButton):
	""" Create a toogle button """
	def __init__(self):
		QtGui.QToolButton.__init__(self)
		
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
	logid=E2init(sys.argv)
	app = EMApp()	# God only knows what this does.... the current GUI scheme is overly complex!
	#app = QtGui.QApplication(sys.argv)
	pm = EMProjectManager()
	pm.show()
	app.exec_()
	E2end(logid)