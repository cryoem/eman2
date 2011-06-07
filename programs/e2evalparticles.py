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
		
		# A widget containing the current particle filename, editable by the user
		# If edited it will also impact set generation !
		self.wptclfile=QtGui.QComboBox(self)
		self.vbl.addWidget(self.wptclfile)

		# Selection tools
		self.wselectg=QtGui.QGroupBox("Select",self)
		self.wselectg.setFlat(False)
		self.vbl.addWidget(self.wselectg)
		
		self.hbl0=QtGui.QHBoxLayout(self.wselectg)

		self.wselallb=QtGui.QPushButton("All")
		self.hbl0.addWidget(self.wselallb)
		
		self.wselnoneb=QtGui.QPushButton("None")
		self.hbl0.addWidget(self.wselnoneb)
		
		self.wsel3db=QtGui.QPushButton("From 3D")
		self.hbl0.addWidget(self.wsel3db)

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
	
		
		QtCore.QObject.connect(self.wfilesel,QtCore.SIGNAL("itemSelectionChanged()"),self.fileUpdate)
		QtCore.QObject.connect(self.wptclfile,QtCore.SIGNAL("currentIndexChanged(int)"),self.ptclChange)
		QtCore.QObject.connect(self.wselallb,QtCore.SIGNAL("clicked(bool)"),self.selAllClasses)
		QtCore.QObject.connect(self.wselnoneb,QtCore.SIGNAL("clicked(bool)"),self.selNoClasses)
		QtCore.QObject.connect(self.wsel3db,QtCore.SIGNAL("clicked(bool)"),self.sel3DClasses)
		
		# View windows, one for class-averages, one for good particles and one for bad particles
		self.vclasses=None
		self.vgoodptcl=None
		self.vbadptcl=None
		
		self.updateFiles()
	
	def ptclChange(self,value):
		self.vgoodptcl.set_data(None)
		self.vbadptcl.set_data(None)
		self.vgoodptcl.setWindowTitle("Included Particles (%s)"%ptclfile)
		self.vbadptcl.setWindowTitle("Excluded Particles (%s)"%ptclfile)
	
	def updateFiles(self):
		subdir=sorted([i for i in os.listdir(e2getcwd()) if "r2d_" in i or "refine_" in i])
		for d in subdir:
			dbs=db_list_dicts("bdb:"+d)
			dbs.sort()
			for db in dbs:
				if "classes_" in db or "allrefs_" in db :
					self.wfilesel.addItem("%s/%s"%(d,db))

		dbs=db_list_dicts("bdb:sets")
		dbs.sort()
		for db in dbs:
			self.wptclfile.addItem("sets#"+db)

	def curFile(self):
		"return the currently selected file as a readable path"
		return "bdb:"+str(self.wfilesel.item(self.wfilesel.currentRow()).text())		# text of the currently selected item

	def curPtclFile(self):
		"return the particle file associated with the currently selected classes file"
		return "bdb:"+str(self.wptclfile.currentText())		# text of the currently selected item
		

	def fileUpdate(self):
		"Called when the user selects a file from the list or need to completely refresh display"
		
		QtGui.qApp.setOverrideCursor(Qt.BusyCursor)
		
		if self.vclasses==None :
			self.vclasses=EMImageMXWidget()
			self.vclasses.setWindowTitle("Classes (%s)"%self.curFile())
			self.vclasses.set_mouse_mode("App")
			QtCore.QObject.connect(self.vclasses,QtCore.SIGNAL("mx_image_selected"),self.classSelect)
			QtCore.QObject.connect(self.vclasses,QtCore.SIGNAL("mx_image_double"),self.classDouble)


		self.classes=EMData.read_images(self.curFile())
		self.vclasses.set_data(self.classes)
		self.vclasses.set_single_active_set(self.curFile()[4:])
		self.vclasses.set_mouse_mode("App")

		# This makes sure the particle file is in the list of choices and is selected
		try:
			ptclfile=self.classes[0]["class_ptcl_src"]
			if ptclfile.lower()[:4]=="bdb:" : ptclfile=ptclfile[4:]
			i=self.wptclfile.findText(ptclfile)
			if i==-1 : 
				self.wptclfile.insertItem(0,ptclfile)
				self.wptclfile.setCurrentIndex(0)
			else:
				self.wptclfile.setCurrentIndex(i)
		except:
			QtGui.QMessageBox.warning(self,"Error !","This image does not appear to be a class average. (No class_ptcl_src, etc.)")
			ptclfile="None"
			
		
		# Make sure our display widgets exist
		if self.vgoodptcl==None :
			self.vgoodptcl=EMImageMXWidget()
		self.vgoodptcl.setWindowTitle("Included Particles (%s)"%self.curPtclFile())

		if self.vbadptcl==None :
			self.vbadptcl=EMImageMXWidget()
		self.vbadptcl.setWindowTitle("Excluded Particles (%s)"%self.curPtclFile())

		self.vclasses.show()
		self.vgoodptcl.show()
		self.vbadptcl.show()
		
		QtGui.qApp.setOverrideCursor(Qt.ArrowCursor)

	def classSelect(self,event,lc):
		"Single clicked class particle. lc=(img#,x,y,image_dict)"
		
		ptclfile=self.curPtclFile()
		try:
			ptclgood=lc[3]["class_ptcl_idxs"]
			self.vgoodptcl.set_data(EMData.read_images(ptclfile,ptclgood))
		except:
			QtGui.QMessageBox.warning(self,"Error !","This image does not appear to be a class average. (No class_ptcl_src, etc.)")
			return
		try:
			ptclbad=lc[3]["exc_class_ptcl_idxs"]
			self.vbadptcl.set_data(EMData.read_images(ptclfile,ptclbad))
		except:
			ptclbad=[]
			self.vbadptcl.set_data(None)
		
		self.vgoodptcl.show()
		self.vbadptcl.show()
		
	def classDouble(self,event,lc):
		self.vclasses.image_set_associate(lc[0],update_gl=True)

	def closeEvent(self,event):
		try : self.vclasses.close()
		except: pass
		try : self.vgoodptcl.close()
		except: pass
		try : self.vbadptcl.close()
		except: pass
		
		QtGui.QWidget.closeEvent(self, event)		

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
		
	def closeEvent(self,event):
		self.wclasstab.close()
		QtGui.QWidget.closeEvent(self, event)		

if __name__ == "__main__":
	main()
