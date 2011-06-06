#!/usr/bin/env python

#
# Author: Steven Ludtke, 06/06/2011 
# Copyright (c) 2011- Baylor College of Medicine
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#

#
from EMAN2 import *
from emimagemx import EMImageMXWidget

import sys
from optparse import OptionParser
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
#from OpenGL import GL,GLU,GLUT
from emapplication import EMApp
import os
from EMAN2db import *
from valslider import *

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog 
	
	This program provides tools for evaluating particle data in various ways. For example it will allow you to select class-averages
	containing bad (or good) particles and manipulate the project to in/exclude them.
	
"""

	parser = OptionParser(usage=usage,version=EMANVERSION)
	
	parser.add_option("--gui",action="store_true",help="Start the GUI for interactive use (default=True)",default=True)
	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n", type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()
	
	#logid=E2init(sys.argv)
	
	app = EMApp()
	control=EMEvalPtclTool(verbose=options.verbose)
	control.show()
	app.execute()

#	E2end(logid)

class EMClassPtclTool(QtGui.QWidget):
	"""This class is a tab widget for inspecting particles within class-averages"""
	
	def __init__(self):
		QtGui.QWidget.__init__(self)
		self.vbl = QtGui.QVBoxLayout(self)
	
		# A listwidget for selecting which class-average file we're looking at
		self.wfilesel=QtGui.QListWidget()
		self.vbl.addWidget(self.wfilesel)
		
		# Mark particles in selected classes as bad
		self.wmarkbad=QtGui.QGroupBox("Mark Bad",self)
		self.wmarkbad.setFlat(False)
		self.vbl.addWidget(self.wmarkbad)
		
		self.vbl2=QtGui.QVBoxLayout(self.wmarkbad)
		
		self.wmarkused=CheckBox(self.wmarkbad,"Included Ptcls",1,100)
		self.vbl2.addWidget(self.wmarkused)
		
		self.wmarkunused=CheckBox(self.wmarkbad,"Excluded Ptcls",1,100)
		self.vbl2.addWidget(self.wmarkunused)
		
		self.wmarkbut=QtGui.QPushButton("Mark as Bad",self.wmarkbad)
		self.vbl2.addWidget(self.wmarkbut)
		
		# Make a new set from selected classes
		self.wmakeset=QtGui.QGroupBox("Make Set",self)
		self.wmakeset.setFlat(False)
		self.vbl.addWidget(self.wmakeset)
		
		self.vbl3=QtGui.QVBoxLayout(self.wmakeset)
		
		self.wmarkused=CheckBox(self.wmakeset,"Included Ptcls",1,100)
		self.vbl3.addWidget(self.wmarkused)
		
		self.wmarkunused=CheckBox(self.wmakeset,"Excluded Ptcls",1,100)
		self.vbl3.addWidget(self.wmarkunused)
		
		self.wsetname=StringBox(self.wmakeset,"Set Name","set-",100)
		self.vbl3.addWidget(self.wsetname)
		
		self.wmakebut=QtGui.QPushButton("Make New Set",self.wmakeset)
		self.vbl3.addWidget(self.wmakebut)
	
		
		# View windows, one for class-averages, one for good particles and one for bad particles
		self.vclasses=None
		self.vgoodptcl=None
		self.vbadptcl=None
		
		self.updatefiles()
		
	def updatefiles(self):
		subdir=sorted([i for i in os.listdir(e2getcwd()) if "r2d_" in i or "refine_" in i])
		for d in subdir:
			dbs=db_list_dicts("bdb:"+d)
			dbs.sort()
			for db in dbs:
				if "classes_" in db or "allrefs_" in db :
					self.wfilesel.addItem("%s/%s"%(d,db))

class EMEvalPtclTool(QtGui.QMainWindow):
	"""This class represents the EMTomoBoxer application instance.  """
	
	def __init__(self,verbose=0):
		QtGui.QMainWindow.__init__(self)
		
		app=QtGui.qApp
		self.setWindowTitle("e2evalparticles.py")
		
		# Menu Bar
		self.mfile=self.menuBar().addMenu("File")
#		self.mfile_save_processed=self.mfile.addAction("Save processed data")
		self.mfile_quit=self.mfile.addAction("Quit")

		self.wtabs=QtGui.QTabWidget()
		self.setCentralWidget(self.wtabs)
		
		self.wclasstab=EMClassPtclTool()
		self.wtabs.addTab(self.wclasstab,"Classes")

		
		# file menu
		QtCore.QObject.connect(self.mfile_quit,QtCore.SIGNAL("triggered(bool)")  ,self.menu_file_quit)

	def menu_file_quit(self):
		self.close()

if __name__ == "__main__":
	main()
