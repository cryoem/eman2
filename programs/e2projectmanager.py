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

class EMProjectManager(QtGui.QMainWindow):
	def __init__(self):
		QtGui.QMainWindow.__init__(self)
		self.makemenues()
		
		font = QtGui.QFont()
		font.setBold(True)
		
		centralwidget = QtGui.QWidget()
		grid = QtGui.QGridLayout()
		grid.addWidget(self.maketiltebarwidget(), 0, 0, 1, 2)
		workflowcontrollabel = QtGui.QLabel("Workflow Control", centralwidget)
		workflowcontrollabel.setFont(font)
		grid.addWidget(workflowcontrollabel, 1,0)
		logbooklabel = QtGui.QLabel("LogBook", centralwidget)
		logbooklabel.setFont(font)
		grid.addWidget(logbooklabel, 1,1)
		grid.addWidget(self.makestackedwidget(),2,0)
		grid.addWidget(self.maketexteditwidget(),2,1)
		centralwidget.setLayout(grid)

		self.setCentralWidget(centralwidget)
		self.makestatusbar()
		
	def makemenues(self):
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
		projectmenu.addAction(newproject)
		newproject.setStatusTip('New Project')
		openproject = QtGui.QAction('Open Project', self)
		openproject.setShortcut('Ctrl+O')
		openproject.setStatusTip('Open Project')
		projectmenu.addAction(openproject)
		editproject = QtGui.QAction('Edit Project', self)
		editproject.setShortcut('Ctrl+E')
		editproject.setStatusTip('Edit Project')
		projectmenu.addAction(editproject)
		
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
		utilsmenu.addAction(filebrowser)
	
		# Help
		helpmenu = menubar.addMenu('&Help')
		about = QtGui.QAction('About', self)
		about.setStatusTip('About')
		helpmenu.addAction(about)
		helpdoc = QtGui.QAction('Help', self)
		helpdoc.setStatusTip('Help')
		helpmenu.addAction(helpdoc)
		
	def maketiltebarwidget(self):
		tbwidget = QtGui.QSplitter(QtCore.Qt.Horizontal)
		tbwidget.setFrameShape(QtGui.QFrame.StyledPanel)
		hbox = QtGui.QHBoxLayout()
		self.PMIcon = QtGui.QLabel("<img src=\"../images/EMAN2Icon.tif\" />", tbwidget)
		hbox.addWidget(self.PMIcon)
		self.PMTitle = QtGui.QLabel("EMAN2 Project Manager ")
		titlefont = QtGui.QFont()
		titlefont.setPointSize(30)
		titlefont.setBold(True)
		self.PMTitle.setFont(titlefont)
		hbox.addWidget(self.PMTitle)
		tbwidget.setLayout(hbox)
		
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
	
	def makestackedwidget(self):
		self.tree_stacked_widget = QtGui.QStackedWidget()
		self.tree_stacked_widget.setMaximumWidth(300)
		self.tree_stacked_widget.addWidget(self.maketreewidget())
		self.tree_stacked_widget.addWidget(self.maketomotreewidget())
		
		return self.tree_stacked_widget
	
	def _on_sprmode(self):
		self.tree_stacked_widget.setCurrentIndex(0)
		
	def _on_tomomode(self):
		self.tree_stacked_widget.setCurrentIndex(1)
		
	def maketomotreewidget(self):
		self.TomoTree = QtGui.QTreeWidget()
		self.TomoTree.setHeaderLabel("Tomography")
		
		self.tomogrpahy = QtGui.QTreeWidgetItem(QtCore.QStringList("Tomographic Particle Reconstruction"))
		self.tomogrpahy.setIcon(0, self.icons['single_image_3d'])
		self.raw_tomographic_files =  QtGui.QTreeWidgetItem(QtCore.QStringList("Raw Tomogram Files"))
		self.raw_tomographic_files.setIcon(0, self.icons['single_image_3d'])
		self.tomogrpahy.addChild(self.raw_tomographic_files)
		self.box_tomogram_particles = QtGui.QTreeWidgetItem(QtCore.QStringList("Box Tomogram Files"))
		self.box_tomogram_particles.setIcon(0, self.icons['green_boxes'])
		self.tomogrpahy.addChild(self.box_tomogram_particles)
		self.tomogram_alignment_averaging = QtGui.QTreeWidgetItem(QtCore.QStringList("Tomogram Alignment and Averaging"))
		self.tomogram_alignment_averaging.setIcon(0, self.icons['tomo_hunter'])
		self.tomogrpahy.addChild(self.tomogram_alignment_averaging)
		
		self.TomoTree.addTopLevelItem(self.tomogrpahy)
		
		return self.TomoTree
		
	def maketreewidget(self):
		self.SPRtree = QtGui.QTreeWidget()
		self.SPRtree.setMinimumHeight(400)
		self.SPRtree.setHeaderLabel("SPR")
		self._load_icons()
		
		# Raw Data
		self.rawdata = QtGui.QTreeWidgetItem(QtCore.QStringList("Raw Data"))
		self.rawdata.setIcon(0, self.icons['single_image'])
		self.filter_raw_data = QtGui.QTreeWidgetItem(QtCore.QStringList("Filter Raw Data"))
		self.filter_raw_data.setIcon(0, self.icons['single_image'])
		self.rawdata.addChild(self.filter_raw_data)
		# Particles
		self.particles = QtGui.QTreeWidgetItem(QtCore.QStringList("Particles"))
		self.particles.setIcon(0, self.icons['green_boxes'])
		self.interactive_boxing = QtGui.QTreeWidgetItem(QtCore.QStringList("Interactive Boxing-  e2boxer"))
		self.interactive_boxing.setIcon(0, self.icons['green_boxes'])
		self.particles.addChild(self.interactive_boxing)
		self.auto_boxing = QtGui.QTreeWidgetItem(QtCore.QStringList("Auto Boxing -e2boxer"))
		self.auto_boxing.setIcon(0, self.icons['green_boxes'])
		self.particles.addChild(self.auto_boxing)
		self.generate_output = QtGui.QTreeWidgetItem(QtCore.QStringList("Generate Output -e2boxer"))
		self.generate_output.setIcon(0, self.icons['green_boxes'])
		self.particles.addChild(self.generate_output)
		self.particle_import = QtGui.QTreeWidgetItem(QtCore.QStringList("Particle Import"))
		self.particle_import.setIcon(0, self.icons['green_boxes'])
		self.particles.addChild(self.particle_import)
		self.coordinate_import = QtGui.QTreeWidgetItem(QtCore.QStringList("Coordinate Import"))
		self.coordinate_import.setIcon(0, self.icons['green_boxes'])
		self.particles.addChild(self.coordinate_import)
		# CTF
		self.ctf = QtGui.QTreeWidgetItem(QtCore.QStringList("CTF"))
		self.ctf.setIcon(0, self.icons['ctf'])
		self.automated_fitted = QtGui.QTreeWidgetItem(QtCore.QStringList("Automated Fitting -e2ctf"))
		self.automated_fitted.setIcon(0, self.icons['ctf'])
		self.ctf.addChild(self.automated_fitted)
		self.interactive_tuning = QtGui.QTreeWidgetItem(QtCore.QStringList("Interactive Tuning -e2ctf"))
		self.interactive_tuning.setIcon(0, self.icons['ctf'])
		self.ctf.addChild(self.interactive_tuning)
		self.generate_output = QtGui.QTreeWidgetItem(QtCore.QStringList("Generate Output - e2ctf"))
		self.generate_output.setIcon(0, self.icons['ctf'])
		self.ctf.addChild(self.generate_output)
		self.generate_structure_factor = QtGui.QTreeWidgetItem(QtCore.QStringList("Generate Structure Factor - e2ctf"))
		self.generate_structure_factor.setIcon(0, self.icons['ctf'])
		self.ctf.addChild(self.generate_structure_factor)
		# Particle Sets
		self.particle_sets = QtGui.QTreeWidgetItem(QtCore.QStringList("Particle Sets"))
		self.particle_sets.setIcon(0, self.icons['multiple_images'])
		self.build_particle_sets = QtGui.QTreeWidgetItem(QtCore.QStringList("Build Particle Sets"))
		self.build_particle_sets.setIcon(0, self.icons['multiple_images'])
		self.particle_sets.addChild(self.build_particle_sets)
		self.evalute_particle_sets = QtGui.QTreeWidgetItem(QtCore.QStringList("Evaluate Particle Sets"))
		self.evalute_particle_sets.setIcon(0, self.icons['multiple_images'])
		self.particle_sets.addChild(self.evalute_particle_sets)
		# Reference Free Class Averages
		self.reference_free_cas = QtGui.QTreeWidgetItem(QtCore.QStringList("Reference Free Class Averages"))
		self.reference_free_cas.setIcon(0, self.icons['web'])
		self.generate_classes = QtGui.QTreeWidgetItem(QtCore.QStringList("Generate Classes - e2refine2d"))
		self.generate_classes.setIcon(0, self.icons['web'])
		self.reference_free_cas.addChild(self.generate_classes)
		# Initial Model
		self.initial_model = QtGui.QTreeWidgetItem(QtCore.QStringList("Initial Model"))
		self.initial_model.setIcon(0, self.icons['single_image_3d'])
		self.make_initialmodel_random = QtGui.QTreeWidgetItem(QtCore.QStringList("Make Model - e2initialmodel"))
		self.make_initialmodel_random.setIcon(0, self.icons['single_image_3d'])
		self.initial_model.addChild(self.make_initialmodel_random)
		# 3D Refinement
		self.threed_refinement = QtGui.QTreeWidgetItem(QtCore.QStringList("3D Refinement"))
		self.threed_refinement.setIcon(0, self.icons['refine'])
		self.e2refine = QtGui.QTreeWidgetItem(QtCore.QStringList("Run e2refine"))
		self.e2refine.setIcon(0, self.icons['refine'])
		self.threed_refinement.addChild(self.e2refine)
		self.e2refinemulti = QtGui.QTreeWidgetItem(QtCore.QStringList("Run e2refinemulti"))
		self.e2refinemulti.setIcon(0, self.icons['refine'])
		self.threed_refinement.addChild(self.e2refinemulti)
		self.frealign = QtGui.QTreeWidgetItem(QtCore.QStringList("Frealign"))
		self.frealign.setIcon(0, self.icons['refine'])
		self.threed_refinement.addChild(self.frealign)
		self.e2refinetofrealign = QtGui.QTreeWidgetItem(QtCore.QStringList("Run e2refinetofrealign"))
		self.e2refinetofrealign.setIcon(0, self.icons['refine'])
		self.frealign.addChild(self.e2refinetofrealign)
		self.e2runfrealign = QtGui.QTreeWidgetItem(QtCore.QStringList("Run e2runfrealign"))
		self.e2runfrealign.setIcon(0, self.icons['refine'])
		self.frealign.addChild(self.e2runfrealign)
		self.e2refinefromfrealign = QtGui.QTreeWidgetItem(QtCore.QStringList("Run e2refinefromfrealign"))
		self.e2refinefromfrealign.setIcon(0, self.icons['refine'])
		self.frealign.addChild(self.e2refinefromfrealign)
		self.eulers = QtGui.QTreeWidgetItem(QtCore.QStringList("Eulers"))
		self.eulers.setIcon(0, self.icons['eulers'])
		# Resolution
		self.resolution =  QtGui.QTreeWidgetItem(QtCore.QStringList("Resolution"))
		self.resolution.setIcon(0, self.icons['resolution'])
		self.e2eotest = QtGui.QTreeWidgetItem(QtCore.QStringList("Run e2eotest"))
		self.e2eotest.setIcon(0, self.icons['resolution'])
		self.resolution.addChild(self.e2eotest)
		self.e2resolution = QtGui.QTreeWidgetItem(QtCore.QStringList("Run e2resolution"))
		self.e2resolution.setIcon(0, self.icons['resolution'])
		self.resolution.addChild(self.e2resolution)
		
		self.SPRtree.addTopLevelItem(self.rawdata)
		self.SPRtree.addTopLevelItem(self.particles)
		self.SPRtree.addTopLevelItem(self.ctf)
		self.SPRtree.addTopLevelItem(self.particle_sets)
		self.SPRtree.addTopLevelItem(self.reference_free_cas)
		self.SPRtree.addTopLevelItem(self.initial_model)
		self.SPRtree.addTopLevelItem(self.threed_refinement)
		self.SPRtree.addTopLevelItem(self.eulers)
		self.SPRtree.addTopLevelItem(self.resolution)
		
		return self.SPRtree
	
	def maketexteditwidget(self):
		self.texteditbox = QtGui.QTextEdit()
		return self.texteditbox
		
	def makestatusbar(self):
		self.statusbar = QtGui.QStatusBar(self)
		self.statusbar.showMessage("Welcome to the EMAN2 Project Manager", 0)
		#self.statusbar.addWidget(QtGui.QLabel("XXXX"))
		self.setStatusBar(self.statusbar)


		
if __name__ == "__main__":
	import sys
	app = QtGui.QApplication(sys.argv)
	pm = EMProjectManager()
	pm.show()
	app.exec_()