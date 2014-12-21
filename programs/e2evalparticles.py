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
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
#from OpenGL import GL,GLU,GLUT
from emapplication import EMApp
import os
from EMAN2db import *
from valslider import *
import traceback

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [classfile]

	This program provides tools for evaluating particle data in various ways. For example it will allow you to select class-averages
	containing bad (or good) particles and manipulate the project to in/exclude them. It will locate files with class-averages automatically,
	but you can specify additional files at the command-line.

"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_header(name="runeval", help='Click Launch to launch the particle evaluation interface', title="### Click Launch to run e2evalparticles ###", row=0, col=0, rowspan=1, colspan=1)
	parser.add_argument("--gui",action="store_true",help="Start the GUI for interactive use (default=True)",default=True)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()

	#logid=E2init(sys.argv, options.ppid)

	app = EMApp()
	control=EMEvalPtclTool(args,verbose=options.verbose)
	control.show()
	app.execute()

#	E2end(logid)

class EMClassPtclTool(QtGui.QWidget):
	"""This class is a tab widget for inspecting particles within class-averages"""

	def __init__(self,extrafiles=None):
		QtGui.QWidget.__init__(self)
		self.vbl = QtGui.QVBoxLayout(self)

		self.extrafiles=extrafiles

		# A listwidget for selecting which class-average file we're looking at
		self.wclassfilel=QtGui.QLabel("Class-average File:")
		self.vbl.addWidget(self.wclassfilel)

		self.wfilesel=QtGui.QListWidget()
		self.vbl.addWidget(self.wfilesel)
		self.vbl.addSpacing(5)

		# A widget containing the current particle filename, editable by the user
		# If edited it will also impact set generation !
		self.wptclfilel=QtGui.QLabel("Particle Data File:")
		self.vbl.addWidget(self.wptclfilel)

		self.wptclfile=QtGui.QComboBox(self)
		self.vbl.addWidget(self.wptclfile)
		self.vbl.addSpacing(5)

		# Selection tools
		self.wselectg=QtGui.QGroupBox("Select",self)
		self.wselectg.setFlat(False)
		self.vbl.addWidget(self.wselectg)
		self.vbl.addSpacing(5)

		self.gbl0=QtGui.QGridLayout(self.wselectg)

		self.wselallb=QtGui.QPushButton("All")
		self.gbl0.addWidget(self.wselallb,0,0)

		self.wselnoneb=QtGui.QPushButton("None")
		self.gbl0.addWidget(self.wselnoneb,0,1)

		self.wselrangeb=QtGui.QPushButton("Range")
		self.gbl0.addWidget(self.wselrangeb,1,0)

		self.wselinvertb=QtGui.QPushButton("Invert")
		self.gbl0.addWidget(self.wselinvertb,0,2)

		self.wsel3db=QtGui.QPushButton("From 3D")
		self.gbl0.addWidget(self.wsel3db,1,2)

		self.wprocessg=QtGui.QGroupBox("Process results",self)
		self.wprocessg.setFlat(False)
		self.vbl.addWidget(self.wprocessg)

		self.vbl2=QtGui.QVBoxLayout(self.wprocessg)

		self.wselused=CheckBox(None,"Included Ptcls",1,100)
		self.vbl2.addWidget(self.wselused)

		self.wselunused=CheckBox(None,"Excluded Ptcls",1,100)
		self.vbl2.addWidget(self.wselunused)

		# Mark particles in selected classes as bad
		self.wmarkbut=QtGui.QPushButton("Mark as Bad")
		self.vbl2.addWidget(self.wmarkbut)

		# Mark particles in selected classes as good
		self.wmarkgoodbut=QtGui.QPushButton("Mark as Good")
		self.vbl2.addWidget(self.wmarkgoodbut)

		# Make a new set from selected classes
		self.wmakebut=QtGui.QPushButton("Make New Set")
		self.vbl2.addWidget(self.wmakebut)
#		self.wmakebut.setEnabled(False)

		# Save list
		self.wsavebut=QtGui.QPushButton("Save Particle List")
		self.vbl2.addWidget(self.wsavebut)

		# Save micrograph dereferenced lists
		self.wsaveorigbut=QtGui.QPushButton("Save CCD-based List")
		self.vbl2.addWidget(self.wsaveorigbut)


		QtCore.QObject.connect(self.wfilesel,QtCore.SIGNAL("itemSelectionChanged()"),self.fileUpdate)
		QtCore.QObject.connect(self.wptclfile,QtCore.SIGNAL("currentIndexChanged(int)"),self.ptclChange)
		QtCore.QObject.connect(self.wselallb,QtCore.SIGNAL("clicked(bool)"),self.selAllClasses)
		QtCore.QObject.connect(self.wselnoneb,QtCore.SIGNAL("clicked(bool)"),self.selNoClasses)
		QtCore.QObject.connect(self.wselrangeb,QtCore.SIGNAL("clicked(bool)"),self.selRangeClasses)
		QtCore.QObject.connect(self.wselinvertb,QtCore.SIGNAL("clicked(bool)"),self.selInvertClasses)
		QtCore.QObject.connect(self.wsel3db,QtCore.SIGNAL("clicked(bool)"),self.sel3DClasses)
		QtCore.QObject.connect(self.wmakebut,QtCore.SIGNAL("clicked(bool)"),self.makeNewSet)
		QtCore.QObject.connect(self.wmarkbut,QtCore.SIGNAL("clicked(bool)"),self.markBadPtcl)
		QtCore.QObject.connect(self.wmarkgoodbut,QtCore.SIGNAL("clicked(bool)"),self.markGoodPtcl)
		QtCore.QObject.connect(self.wsavebut,QtCore.SIGNAL("clicked(bool)"),self.savePtclNum)
		QtCore.QObject.connect(self.wsaveorigbut,QtCore.SIGNAL("clicked(bool)"),self.saveOrigPtclNum)

		# View windows, one for class-averages, one for good particles and one for bad particles
		self.vclasses=None
		self.vgoodptcl=None
		self.vbadptcl=None

		self.updateFiles()

	def makeNewSet(self,x):
		"Makes a new particle set based on the selected class-averages"
		setname=QtGui.QInputDialog.getText(None,"Set Name","Please specify the name for the set. If you specify an existing set, new particles will be added to the end")
		if setname[1]==False : return
		else: setname=setname[0]
		if setname[-4:]!=".lst" : setname=setname+".lst"
		if not "/" in setname : setname="sets/"+setname

		lst=LSXFile(self.curPtclFile())		# lst file for dereferenceing
		lstout=LSXFile(setname)
		include=[]
		# iterate over each particle from each marked class-average
		for n in self.curPtclIter(self.wselused.getValue(),self.wselunused.getValue()):
			try :
				orign,origfile,comment=lst.read(n)			# the original file/number dereferenced from the LST file
			except:
				QtGui.QMessageBox.warning(self,"Error !","The data_source '%s' does not follow EMAN2.1 project conventions. Cannot find raw particles for set."%srcfile)
				return

			include.append((origfile,orign,comment))		# build a list so we can sort by frame
		
		# write the new set
		for i in sorted(include) : lstout.write(-1,i[1],i[0],i[2])

	def markBadPtcl(self,x):
		"Mark particles from the selected class-averages as bad in the set interface"

		r=QtGui.QMessageBox.question(None,"Are you sure ?","WARNING: There is no undo for this operation. It will  mark all particles associated with the selected class-averages as bad. Are you sure you want to proceed ?",QtGui.QMessageBox.Yes|QtGui.QMessageBox.Cancel)
		if r==QtGui.QMessageBox.Cancel : return

		lst=LSXFile(self.curPtclFile())		# lst file for dereferenceing
		ptcls={}						# dictionary keyed by original frame filename with list of selected particle #s
		# iterate over each particle from each marked class-average
		for n in self.curPtclIter(self.wselused.getValue(),self.wselunused.getValue()):
			try :
				orign,origfile,comment=lst.read(n)
			except:
				QtGui.QMessageBox.warning(self,"Error !","The data_source '%s' does not follow EMAN2.1 project conventions. Cannot find raw particles for set."%srcfile)
				return

			try: ptcls[origfile].append(orign)		# try to add to a list for an existing filename
			except: ptcls[origfile]=[orign]			# creates the list for this filename if it's new

		#now mark the particles as bad
		newbad=0
		totbad=0
		for origfile in ptcls:
			js=js_open_dict(info_name(origfile))	# get the info dict for this file

			try: sets=js["sets"]
			except: sets={"bad_particles":[]}
			try: badset=set(sets["bad_particles"])
			except: badset=set()

			try:
				newset=list(set(ptcls[origfile])|badset)
				sets["bad_particles"]=newset	# update the set of bad particles for this image file
				js["sets"]=sets
				totbad+=len(badset)
				newbad+=len(newset)-len(badset)
			except:
				print "Error setting bad particles in ",origfile
				
		print newbad, " new particles marked as bad. Total of ",totbad," in affected micrographs"

	def markGoodPtcl(self,x):
		"Mark particles from the selected class-averages as good in the set interface"

		r=QtGui.QMessageBox.question(None,"Are you sure ?","WARNING: There is no undo for this operation. It will un-mark all particles associated with the selected class-averages as bad. Are you sure you want to proceed ?",QtGui.QMessageBox.Yes|QtGui.QMessageBox.Cancel)
		if r==QtGui.QMessageBox.Cancel : return

		lst=LSXFile(self.curPtclFile())		# lst file for dereferenceing
		ptcls={}						# dictionary keyed by original frame filename with list of selected particle #s
		# iterate over each particle from each marked class-average
		for n in self.curPtclIter(self.wselused.getValue(),self.wselunused.getValue()):
			try :
				orign,origfile,comment=lst.read(n)
			except:
				QtGui.QMessageBox.warning(self,"Error !","The data_source '%s' does not follow EMAN2.1 project conventions. Cannot find raw particles for set."%srcfile)
				return

			try: ptcls[origfile].append(orign)		# try to add to a list for an existing filename
			except: ptcls[origfile]=[orign]			# creates the list for this filename if it's new

		#now mark the particles as good
		badafter=0
		badbefore=0
		for origfile in ptcls:
			js=js_open_dict(info_name(origfile))	# get the info dict for this file
			try:
				badset=set(js["sets"]["bad_particles"])
				js["sets"]["bad_particles"]=list(badset-set(ptcls[origfile]))	# update the set of bad particles for this image file
			except:
				pass		# since marking as good is the same as removing from the bad list, if there is no bad list, there is nothing to do

			try: sets=js["sets"]
			except: continue	# if no current bad particles, nothing to mark good
			try: badset=sets["bad_particles"]
			except: continue

			try:
				newset=list(badset-set(ptcls[origfile]))
				sets["bad_particles"]=newset	# update the set of bad particles for this image file
				js["sets"]=sets
				badbefore+=len(badset)
				badafter+=len(newset)
			except:
				continue
		
		print badbefore," bad particles before processing, now ",badafter

	def savePtclNum(self,x):
		"Saves a list of particles from marked classes into a text file"

		filename=QtGui.QInputDialog.getText(None,"Filename","Please enter a filename for the particle list. The file will contain the particle number (within the particle file) for each particle associated with a selected class-average.")
		if filename[1]==False or filename[0]=="" : return

		out=file(filename[0],"w")
		for i in self.curPtclIter(self.wselused.getValue(),self.wselunused.getValue()): out.write("%d\n"%i)
		out.close()

	def saveOrigPtclNum(self,x):
		"Saves a file containing micrograph-dereferenced particles"
		filename=QtGui.QInputDialog.getText(None,"Filename","Please enter a filename for the particle list. The file will contain particle number and image file, one per line. Image files will be referenced back to the original per-CCD frame stacks.")
		if filename[1]==False or filename[0]=="" : return

		lst=LSXFile(self.curPtclFile())		# lst file for dereferenceing
		include=[]
		# iterate over each particle from each marked class-average
		for n in self.curPtclIter(self.wselused.getValue(),self.wselunused.getValue()):
			try :
				orign,origfile,comment=lst.read(n)			# the original file/number dereferenced from the LST file
			except:
				QtGui.QMessageBox.warning(self,"Error !","The data_source '%s' does not follow EMAN2.1 project conventions. Cannot find raw particles for set."%srcfile)
				return

			include.append((origfile,orign,comment))		# build a list so we can sort by frame
		
		# write the output file
		out=file(filename,"w")
		for i in sorted(include) : out.write("{}\t{}\n".format(i[1],i[0]))
		out=None

	def selAllClasses(self,x):
		"Mark all classes as selected"
		self.vclasses.all_set()

	def selNoClasses(self,x):
		"Clear selection"
		self.vclasses.clear_set()

	def selRangeClasses(self,x):
		"Select a range of images (ask the user for the range)"
		rng=QtGui.QInputDialog.getText(None,"Select Range","Enter the range of particle values as first-last (inclusive). Merges with existing selection.")
		if rng[1]==False : return

		try:
			x0,x1=rng[0].split("-")
			x0=int(x0)
			x1=int(x1)+1
		except:
			QtGui.QMessageBox.warning(self,"Error !","Invalid range specified. Use: min-max")
			return

		self.vclasses.subset_set(range(x0,x1))

	def selInvertClasses(self,x):
		"Inverts the current selection set"
		self.vclasses.invert_set()

	def sel3DClasses(self,x):
		"Select a range of images based on those used in a 3-D reconstruction associated with this classes file. Removes current selection first."

		f=self.curFile()
		if not '#classes_' in f :
			QtGui.QMessageBox.warning(self,"Error !","A classes_xx file from a refine_xx directory is not currently selected")
			return

		# construct the path to the threed_xx file
		num=f.split("_")[-1]
		pre=f.split("#")[0]
		d3path="%s#threed_%s"%(pre,num)
		try:
			a=EMData(d3path,0,True)
			goodptcl=a["threed_ptcl_idxs"]
		except:
			QtGui.QMessageBox.warning(self,"Error !","Cannot read classes from "+d3path)
			return

		self.vclasses.clear_set()
		self.vclasses.subset_set(goodptcl)

	def ptclChange(self,value):
		"Called when a change of particle data file occurs to zero out the display"
		try:
			self.vgoodptcl.set_data(None)
			self.vbadptcl.set_data(None)
		except:
			pass

	def updateFiles(self):
		"Updates the list of classes files"
		subdir=sorted([i for i in os.listdir(e2getcwd()) if "r2d_" in i or "relion2d_" in i or "refine_" in i or "multi_" in i])
		for d in subdir:
			try: dbs=os.listdir(d)
			except: continue
			dbs.sort()
			for db in dbs:
				if "classes_" in db or "allrefs_" in db :
					self.wfilesel.addItem("%s/%s"%(d,db))

		for f in self.extrafiles:
			self.wfilesel.addItem(f)

		dbs=os.listdir("sets")
		dbs.sort()
		for db in dbs:
			self.wptclfile.addItem("sets/"+db)

	def curPtclIter(self,included=True,excluded=True):
		"This is a generator function which yields n (in self.curPtclFile()) for all particles associated with selected classes"
		for ci in self.curSet():
			try :
				c=EMData(self.curFile(),ci,True)		# read header for current class average
				if included :
					incl=c["class_ptcl_idxs"]
					if isinstance(incl,int) : incl=[incl]		# This should not happen, but seems to sometimes for some reason
					for i in incl:
						yield(i)
				if excluded and c.has_attr("exc_class_ptcl_idxs"):
					excl=c["exc_class_ptcl_idxs"]
					if isinstance(excl,int) : excl=[excl]		# This should not happen, but seems to sometimes for some reason
					for i in excl:
						yield(i)
			except:
				print "Problem with class %d (%s). Skipping"%(ci,self.curFile())
				traceback.print_exc()
				continue

	def curFile(self):
		"return the currently selected file as a readable path"
		return str(self.wfilesel.item(self.wfilesel.currentRow()).text())		# text of the currently selected item

	def curSet(self):
		"return a set (integers) of the currently selected class-averages"
		
		return self.vclasses.get_set("evalptcl")

	def curPtclFile(self):
		"return the particle file associated with the currently selected classes file"
		return str(self.wptclfile.currentText())		# text of the currently selected item


	def fileUpdate(self):
		"Called when the user selects a file from the list or need to completely refresh display"

		QtGui.qApp.setOverrideCursor(Qt.BusyCursor)

		if self.vclasses==None :
			self.vclasses=EMImageMXWidget()
			self.vclasses.set_mouse_mode("App")
			QtCore.QObject.connect(self.vclasses,QtCore.SIGNAL("mx_image_selected"),self.classSelect)
			QtCore.QObject.connect(self.vclasses,QtCore.SIGNAL("mx_image_double"),self.classDouble)

		self.vclasses.set_title("Classes")


#		self.classes=EMData.read_images(self.curFile())
		self.vclasses.set_data(self.curFile(),self.curFile())
#		self.vclasses.set_single_active_set("selected")		# This makes the 'set' representing the selected class-averages current
		self.vclasses.set_mouse_mode("App")
		self.vclasses.enable_set("evalptcl",[])

		# This makes sure the particle file is in the list of choices and is selected
		try:
			ptclfile=EMData(self.curFile(),0,True)["class_ptcl_src"]
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
		self.vgoodptcl.set_title("Included Particles")

		if self.vbadptcl==None :
			self.vbadptcl=EMImageMXWidget()
		self.vbadptcl.set_title("Excluded Particles")

		self.vclasses.show()
		self.vgoodptcl.show()
		self.vbadptcl.show()

		QtGui.qApp.setOverrideCursor(Qt.ArrowCursor)

	def classSelect(self,event,lc):
		"Single clicked class particle. lc=(img#,x,y,image_dict)"

		QtGui.qApp.setOverrideCursor(Qt.BusyCursor)
		ptclfile=self.curPtclFile()
		try:
			ptclgood=lc[3]["class_ptcl_idxs"]
			self.vgoodptcl.set_data(EMData.read_images(ptclfile,ptclgood))
		except:
			QtGui.QMessageBox.warning(self,"Error !","This image does not appear to be a class average. (No class_ptcl_src, etc.)")
			QtGui.qApp.setOverrideCursor(Qt.ArrowCursor)
			return
		try:
			ptclbad=lc[3]["exc_class_ptcl_idxs"]
			self.vbadptcl.set_data(EMData.read_images(ptclfile,ptclbad))
		except:
			ptclbad=[]
			self.vbadptcl.set_data(None)

		self.vgoodptcl.show()
		self.vbadptcl.show()
		QtGui.qApp.setOverrideCursor(Qt.ArrowCursor)

	def classDouble(self,event,lc):
		self.vclasses.image_set_associate(lc[0],update_gl=True)

	def closeEvent(self,event):
		try :
			self.vclasses.commit_sets()
			self.vclasses.close()
		except: pass
		try : self.vgoodptcl.close()
		except: pass
		try : self.vbadptcl.close()
		except: pass

		QtGui.QWidget.closeEvent(self, event)

class EMEvalPtclTool(QtGui.QMainWindow):
	"""This class represents the EMTomoBoxer application instance.  """

	def __init__(self,extrafiles=None,verbose=0):
		QtGui.QMainWindow.__init__(self)

		app=QtGui.qApp
		self.setWindowTitle("e2evalparticles")

		# Menu Bar
		self.mfile=self.menuBar().addMenu("File")
#		self.mfile_save_processed=self.mfile.addAction("Save processed data")
		self.mfile_quit=self.mfile.addAction("Quit")

		self.wtabs=QtGui.QTabWidget()
		self.setCentralWidget(self.wtabs)

		self.wclasstab=EMClassPtclTool(extrafiles)
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
