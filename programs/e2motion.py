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
		self.mfile_open=self.mfile.addAction("Open")
		self.mfile_quit=self.mfile.addAction("Quit")

		self.mwin=self.menuBar().addMenu("Window")
		self.mwin_boxes=self.mwin.addAction("Particles")
		self.mwin_single=self.mwin.addAction("Single Particle")
		self.mwin_average=self.mwin.addAction("Averaging")


		self.setCentralWidget(QtGui.QWidget())
		self.gbl = QtGui.QGridLayout(self.centralWidget())
		cen=self.centralWidget()
		
		# widget for editing the alignment mask
		self.wlalimaskdraw=QtGui.QLabel("<big><b>Edit</b></big>")
		self.wlalimaskdraw.setAlignment(Qt.AlignHCenter)
		self.gbl.addWidget(self.wlalimaskdraw,0,1)
		
		self.wlalimaskdraw2=QtGui.QLabel("A\nl\ni\ng\nn")
		self.gbl.addWidget(self.wlalimaskdraw2,1,0)
		
		self.w2dalimaskdraw=EMImage2DWidget()
		self.gbl.addWidget(self.w2dalimaskdraw,1,1)
		
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
		
		self.vbl1.addStretch(5)
		
		# widget for displaying the masked alignment reference
		self.wlalimask=QtGui.QLabel("Ref")
		self.wlalimask.setAlignment(Qt.AlignHCenter)
		self.gbl.addWidget(self.wlalimask,0,3)
		
		self.w2dalimask=EMImage2DWidget()
		self.gbl.addWidget(self.w2dalimask,1,3)

		# widget for editing the ROI mask
		self.wlroimaskdraw=QtGui.QLabel("R\nO\nI")
		self.gbl.addWidget(self.wlroimaskdraw,4,0)
		
		self.w2clsmaskdraw=EMImage2DWidget()
		self.gbl.addWidget(self.w2clsmaskdraw,4,1)

		# Widget for setting alignment mask blur and base level
		self.vbl2=QtGui.QVBoxLayout()
		self.gbl.addLayout(self.vbl2,4,2)
		self.vbl2.addStretch(5)
		
		self.wlclsmaskblur=QtGui.QLabel("Blur")
		self.vbl2.addWidget(self.wlclsmaskblur)
		
		self.wsbclsmaskblur=QtGui.QSpinBox()
		self.wsbclsmaskblur.setRange(0,25)
		self.vbl2.addWidget(self.wsbclsmaskblur)

		# widget for displaying the masked ROI
		self.w2clsmask=EMImage2DWidget()
		self.gbl.addWidget(self.w2clsmask,4,3)
		
		self.vbl2.addStretch(5)

		# Widget showing lists of different result sets
		self.wlistresult=QtGui.QListWidget()
		self.gbl.addWidget(self.wlistresult,1,6,4,1)

		if path!=None : self.initPath(path)

	def initPath(self,path):
		return

if __name__ == "__main__":
	main()
