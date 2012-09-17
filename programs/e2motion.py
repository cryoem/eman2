#!/usr/bin/env python

#
# Author: Steven Ludtke  9/14/2012 
# Copyright (c) 2012- Baylor College of Medicine
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

import sys
import os
import weakref
from sys import argv
from EMAN2 import *
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import Qt
from emapplication import get_application, EMApp
from emimage2d import EMImage2DWidget
from emimagemx import EMImageMXWidget
from valslider import *
import embrowser

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options]

	WARNING: This program still under development.
	
	This program is used to interactively generate sequences of class-averages from sets of particles. It can be
	used for many purposes, but is primarily intended at studies of macromolecular dynamics and variability. A stack
	of particles ostensibly in the same 3-D (but not 2-D) orientation are read in, then alignment, classification
	and averaging is performed to produce pseudo time-series animations detailing some aspect of the structure's
	variability.
	
	This program is NOT designed for use with stacks of tens of thousands of images. All images are kept in system
	memory. All images are theoretically supposed to share a common overall orientation. That is, for a normal single
	particle project, you would first subclassify the overall image stack used for 3-D using e2refine2d.py or
	e2refine.py, then run this program only on a subset of particles.
	
	If an existing project is specified with --path, previous results will be re-opened"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--path",type=str,default=None,help="Path for the refinement, default=auto")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	global options
	(options, args) = parser.parse_args()
	
	if options.path==None:
		options.path=numbered_path("motion",True)
		os.makedirs(options.path)

	pid=E2init(argv)
	
	app = EMApp()
	motion=EMMotion(app,options.path)
	motion.show()
	app.execute()
	
	E2end(pid)

class EMMotion(QtGui.QMainWindow):
	"""This is the main window for the EMMotion application"""
	
	def __init__(self,application,path=None):
		"""application is an QApplication instance. ptclstack is the path to the file containing the particles to analyze. path is the path for ouput files""" 
		QtGui.QWidget.__init__(self)

		self.app=weakref.ref(application)
		self.path=path

		self.setWindowTitle("Main Window (e2spt_boxer.py)")

#		self.setWindowTitle("e2spt_boxer.py")
		
		# Menu Bar
		self.mfile=self.menuBar().addMenu("File")
		self.mfileopen=self.mfile.addAction("Select Particles")
		self.mfileopencls=self.mfile.addAction("Particles from Classes")
		self.mfileopencls.setEnabled(False)
		self.mfilequit=self.mfile.addAction("Quit")

		#self.mwin=self.menuBar().addMenu("Window")
		#self.mwin_boxes=self.mwin.addAction("Particles")
		#self.mwin_single=self.mwin.addAction("Single Particle")
		#self.mwin_average=self.mwin.addAction("Averaging")


		self.setCentralWidget(QtGui.QWidget())
		self.gbl = QtGui.QGridLayout(self.centralWidget())
		cen=self.centralWidget()
		
		###### Alignment Mask
		# widget for editing the alignment mask
		self.wlalimaskdraw=QtGui.QLabel("<big><pre>Edit</pre></big>")
		self.wlalimaskdraw.setAlignment(Qt.AlignHCenter)
		self.gbl.addWidget(self.wlalimaskdraw,0,1)
		
		self.wlalimaskdraw2=QtGui.QLabel("<big><pre>A\nl\ni\ng\nn</pre></big>")
		self.gbl.addWidget(self.wlalimaskdraw2,1,0)
		
		self.w2dalimaskdraw=EMImage2DWidget()
		self.gbl.addWidget(self.w2dalimaskdraw,1,1)
		
		# Buttons for controlling mask
		self.hbl1=QtGui.QHBoxLayout()
		self.gbl.addLayout(self.hbl1,2,1)
		self.hbl1.addStretch(5)
		
		self.wbautoali=QtGui.QPushButton("Auto")
		self.hbl1.addWidget(self.wbautoali)
		
		self.wbresetali=QtGui.QPushButton("Reset")
		self.hbl1.addWidget(self.wbresetali)

		self.hbl1.addStretch(5)

		# Widget for setting alignment mask blur and base level
		self.vbl1=QtGui.QVBoxLayout()
		self.gbl.addLayout(self.vbl1,1,2)
		self.vbl1.addStretch(5)
		
		self.wlalimaskblur=QtGui.QLabel("Blur")
		self.vbl1.addWidget(self.wlalimaskblur)
		
		self.wsbalimaskblur=QtGui.QSpinBox()
		self.wsbalimaskblur.setRange(0,25)
		self.vbl1.addWidget(self.wsbalimaskblur)
		
		self.vbl1.addSpacing(16)
		
		self.wlalimaskbase=QtGui.QLabel("Base")
		self.vbl1.addWidget(self.wlalimaskbase)
		
		self.wsbalimaskbase=QtGui.QSpinBox()
		self.wsbalimaskbase.setRange(0,100)
		self.vbl1.addWidget(self.wsbalimaskbase)
		
		self.vbl1.addSpacing(16)
		
		self.wbaligo=QtGui.QPushButton(QtCore.QChar(0x2192))
		self.vbl1.addWidget(self.wbaligo)
		
		self.vbl1.addStretch(5)
		
		# widget for displaying the masked alignment reference
		self.wlalimask=QtGui.QLabel("<big><pre>Reference</pre></big>")
		self.wlalimask.setAlignment(Qt.AlignHCenter)
		self.gbl.addWidget(self.wlalimask,0,3)
		
		self.w2dalimask=EMImage2DWidget()
		self.gbl.addWidget(self.w2dalimask,1,3)
		
		self.hbl1a=QtGui.QHBoxLayout()
		self.gbl.addLayout(self.hbl1a,2,3)
		self.hbl1a.addStretch(5)
		
		self.wbrecalcref=QtGui.QPushButton("Recompute")
		self.hbl1a.addWidget(self.wbrecalcref)
		
		self.hbl1a.addStretch(5)


		###### ROI Mask
		# widget for editing the ROI mask
		self.wlroimaskdraw=QtGui.QLabel("<big><pre>R\nO\nI</pre></big>")
		self.gbl.addWidget(self.wlroimaskdraw,4,0)
		
		self.w2clsmaskdraw=EMImage2DWidget()
		self.gbl.addWidget(self.w2clsmaskdraw,4,1)

		# Buttons for controlling mask
		self.hbl2=QtGui.QHBoxLayout()
		self.gbl.addLayout(self.hbl2,5,1)
		self.hbl2.addStretch(5)
		
		self.wbautoroi=QtGui.QPushButton("Auto")
		self.hbl2.addWidget(self.wbautoroi)
		
		self.wbresetroi=QtGui.QPushButton("Reset")
		self.hbl2.addWidget(self.wbresetroi)
		
		self.hbl2.addStretch(5)

		# Widget for setting alignment mask blur and base level
		self.vbl2=QtGui.QVBoxLayout()
		self.gbl.addLayout(self.vbl2,4,2)
		self.vbl2.addStretch(5)
		
		self.wlclsmaskblur=QtGui.QLabel("Blur")
		self.vbl2.addWidget(self.wlclsmaskblur)
		
		self.wsbclsmaskblur=QtGui.QSpinBox()
		self.wsbclsmaskblur.setRange(0,25)
		self.vbl2.addWidget(self.wsbclsmaskblur)

		self.vbl2.addSpacing(16)
		
		self.wbroigo=QtGui.QPushButton(QtCore.QChar(0x2192))
		self.vbl2.addWidget(self.wbroigo)


		# widget for displaying the masked ROI
		self.w2clsmask=EMImage2DWidget()
		self.gbl.addWidget(self.w2clsmask,4,3)
		
		self.vbl2.addStretch(5)

		self.wlarrow1=QtGui.QLabel(QtCore.QChar(0x2192))
		self.gbl.addWidget(self.wlarrow1,2,4)

		###### Results
		# Widget showing lists of different result sets
		self.vbl3=QtGui.QVBoxLayout()
		self.gbl.addLayout(self.vbl3,1,6,5,1)
		
		self.wllistresult=QtGui.QLabel("Results")
#		self.wllistresult.setAlignment(Qt.AlignHCenter)
		self.vbl3.addWidget(self.wllistresult)
		
		self.wlistresult=QtGui.QListWidget()
		self.vbl3.addWidget(self.wlistresult)

		###### Parameters for processing
		self.vgb1=QtGui.QGroupBox("Launch Job")
		self.vbl3.addWidget(self.vgb1)
		
		self.vbl3a=QtGui.QVBoxLayout()
		self.vgb1.setLayout(self.vbl3a)
		
		self.wvbclasses=ValBox(None,(0,256),"# Classes",32)
		self.wvbclasses.setIntonly(True)
		self.vbl3a.addWidget(self.wvbclasses)
		
		self.wvbcores=ValBox(None,(0,256),"# Threads",2)
		self.wvbcores.setIntonly(True)
		self.vbl3a.addWidget(self.wvbcores)
		
		self.wcbprocmode=QtGui.QComboBox()
		self.wcbprocmode.addItem("PCA / k-means")
		self.wcbprocmode.addItem("Average Density")
		self.vbl3a.addWidget(self.wcbprocmode)
		
		self.wpbprogress=QtGui.QProgressBar()
		self.wpbprogress.setEnabled(False)
		self.wpbprogress.setMinimum(0)
		self.wpbprogress.setMaximum(100)
		self.wpbprogress.setValue(0)
		self.vbl3a.addWidget(self.wpbprogress)
		
		# doubles as a cancel button
		self.wbcompute=QtGui.QPushButton("Compute")
		self.vbl3a.addWidget(self.wbcompute)

		self.wlarrow2=QtGui.QLabel(QtCore.QChar(0x2192))
		self.gbl.addWidget(self.wlarrow2,2,7)


		###### Output widgets
		# Class-averages
		self.wlclasses=QtGui.QLabel("<big><pre>Classes</pre></big>")
		self.wlclasses.setAlignment(Qt.AlignHCenter)
		self.gbl.addWidget(self.wlclasses,0,9)
		
		self.w2dclasses=EMImage2DWidget()
		self.gbl.addWidget(self.w2dclasses,1,9)
		
		## Buttons for controlling mask
		#self.hbl1=QtGui.QHBoxLayout()
		#self.gbl.addLayout(self.hbl1,2,1)
		#self.hbl1.addStretch(5)
		
		#self.wbautoali=QtGui.QPushButton("Auto")
		#self.hbl1.addWidget(self.wbautoali)
		
		#self.wbresetali=QtGui.QPushButton("Reset")
		#self.hbl1.addWidget(self.wbresetali)"bdb:%s"

		#self.hbl1.addStretch(5)

		QtCore.QObject.connect(self.wbautoali,QtCore.SIGNAL("clicked(bool)"),self.aliAutoPress)
		QtCore.QObject.connect(self.wbresetali,QtCore.SIGNAL("clicked(bool)"),self.aliResetPress)
		QtCore.QObject.connect(self.wbaligo,QtCore.SIGNAL("clicked(bool)"),self.aliGoPress)
		QtCore.QObject.connect(self.wbrecalcref,QtCore.SIGNAL("clicked(bool)"),self.aliRecalcRefPress)
		QtCore.QObject.connect(self.wbautoroi,QtCore.SIGNAL("clicked(bool)"),self.roiAutoPress)
		QtCore.QObject.connect(self.wbresetroi,QtCore.SIGNAL("clicked(bool)"),self.roiResetPress)
		QtCore.QObject.connect(self.wbroigo,QtCore.SIGNAL("clicked(bool)"),self.roiGoPress)
		QtCore.QObject.connect(self.wbcompute,QtCore.SIGNAL("clicked(bool)"),self.doCompute)

		QtCore.QObject.connect(self.mfileopen,QtCore.SIGNAL("triggered(bool)")  ,self.menuFileOpen  )


		# set up draw mode
		insp=self.w2dalimaskdraw.get_inspector()
		insp.hide()
		insp.mmtab.setCurrentIndex(3)
		insp.dtpenv.setText("0.0")
		
		insp=self.w2dalimaskdraw.get_inspector()
		insp.hide()
		insp.mmtab.setCurrentIndex(3)
		insp.dtpenv.setText("0.0")
		

		if path!=None : self.initPath(path)

	def initPath(self,path):
		if path[:4].lower()!="bdb:" : path="bdb:"+path
		self.path=path
		dicts=db_list_dicts(self.path)
		
		# If particles exists then we can fully initialize
		if db_check_dict("%s#particles"%self.path) :
			self.particles=EMData.read_images("%s#particles"%self.path,0,-1)	# read in the entire particle stack
			
			cl=[i for i in dicts if "classes" in i]
			cl.sort()
			self.wllistresult.addItems(QtCore.QStringList(cl))
			
			self.bootstrap()
		else :
			self.particles=None
			
		return

	def menuFileOpen(self,x):
		if self.particles!=None:
			QtGui.QMessageBox.warning(None,"Error","%s already contains a stack of particles. A new folder is required to start with a new stack of particles. Rerun without --path option."%self.path)
			return

		self.dialog = embrowser.EMBrowserWidget(withmodal=False,multiselect=False)
		QtCore.QObject.connect(dialog,QtCore.SIGNAL("ok"),self.gotPath)
		QtCore.QObject.connect(dialog,QtCore.SIGNAL("cancel"),self.gotPath)
	
	def gotPath(self):
		ptcl=self.dialog.getResult()
		self.dialog=None
		if ptcl == None : return		# user pressed cancel, but we still need to clean up
		
		n=EMUtil.get_image_count(ptcl)
		sz=EMData(ptcl,0,True)["nx"]
		
		# We don't want the user to overwhelm the system
		if sz*sz*4*n > 5.0e8 :
			nmax=5.0e8/(sz*sz*4)
			r=QtGui.QMessageBox.question(None,"Are you sure ?","WARNING: This full particle set will require %d+ gb memory to process. Select Yes to use only the first %d particles, No to use the entire stack, or Cancel to abort. You may also consider using e2proc2d.py --meanshrink= to produce a downsampled stack for processing."%(int(sz*sz*9*n/1.0e9+.5),nmax),QtGui.QMessageBox.Yes|QtGui.QMessageBox.No|QtGui.QMessageBox.Cancel)
			if r==QtGui.QMessageBox.Cancel : return
			if r==QtGui.QMessageBox.Yes : n=nmax
			
		launch_childprocess("e2bdb.py %s  --makevstack=%s#particles --step=0,1,%d"%(path,self.path,n))
		
		self.initPath(self.path)
	
	def bootstrap(self):
		"""This will create an initial alignment for the particles when no reference is available"""
		
		
		
	def aliAutoPress(self,x):
		pass
	
	def aliResetPress(self,x):
		pass
	
	def aliGoPress(self,x):
		pass
	
	def aliRecalcRefPress(self,x):
		pass
	
	def roiAutoPress(self,x):
		pass
	
	def roiResetPress(self,x):
		pass
	
	def roiGoPress(self,x):
		pass
	
	def doCompute(self,x):
		pass
	
	

if __name__ == "__main__":
	main()
