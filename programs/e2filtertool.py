#!/usr/bin/env python

#
# Author: Steven Ludtke  3/4/2011
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA
#

from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import Qt
from PyQt4.QtCore import QTimer

import sys
import os
import weakref
import threading

from EMAN2 import *
from emapplication import get_application, EMApp
from emimage2d import EMImage2DWidget
from emimagemx import EMImageMXWidget
from emscene3d import EMScene3D
from emdataitem3d import EMDataItem3D, EMIsosurface
from emshape import EMShape
from valslider import *
import traceback

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <image file>

	Provides a GUI interface for applying a sequence of processors to an image, stack of images or a volume."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

#	parser.add_argument("--boxsize","-B",type=int,help="Box size in pixels",default=64)
#	parser.add_argument("--shrink",type=int,help="Shrink factor for full-frame view, default=0 (auto)",default=0)
	parser.add_argument("--apix",type=float,help="Override the A/pix value stored in the file header",default=0.0)
#	parser.add_argument("--force2d",action="store_true",help="Display 3-D data as 2-D slices",default=False)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()

	if len(args) != 1: parser.error("You must specify a single data file on the command-line.")
	if not file_exists(args[0]): parser.error("%s does not exist" %args[0])
#	if options.boxsize < 2: parser.error("The boxsize you specified is too small")
#	# The program will not run very rapidly at such large box sizes anyhow
#	if options.boxsize > 2048: parser.error("The boxsize you specified is too large.\nCurrently there is a hard coded max which is 2048.\nPlease contact developers if this is a problem.")

#	logid=E2init(sys.argv,options.ppid)

	app = EMApp()
	pix_init()
	control=EMFilterTool(datafile=args[0],apix=options.apix,force2d=False,verbose=options.verbose)
#	control=EMFilterTool(datafile=args[0],apix=options.apix,force2d=options.force2d,verbose=options.verbose)
	control.show()
	try: control.raise_()
	except: pass

	app.execute()
#	E2end(logid)

def filtchange(name,value):
	return {}

class EMProcessorWidget(QtGui.QWidget):
	"""A single processor with parameters"""

	plist=dump_processors_list()

	# Sorted list of the stuff before the first '.'
	cats=set([i.split(".")[0] for i in plist.keys()])
	cats=list(cats)
	cats.sort()

	# For some parameters we will try to assign sensible default values and ranges. This information should probably be a part
	# of the processor objects themselves, but isn't at present.
	# (start enabled, range, default value, change callback)
	parmdefault= {
		"cutoff_abs":(1,(0,.5),.1,filtchange),
		"cutoff_freq":(0,(0,1.0),None,filtchange),
		"cutoff_pixels":(0,(0,128),None,filtchange),
		"cutoff_resolv":(0,(0,1.0),None,filtchange),
		"sigma":(0,(0,3.0),.5,None),
		"apix":(0,(0.2,10.0),2.0,None)
	}

	def __init__(self,parent=None,tag=None):
		app=QtGui.qApp

		QtGui.QWidget.__init__(self,parent)
		self.gbl = QtGui.QGridLayout(self)
		self.gbl.setColumnStretch(0,0)
		self.gbl.setColumnStretch(1,0)
		self.gbl.setColumnStretch(2,1)
		self.gbl.setColumnStretch(3,3)

		# Enable checkbox
		self.wenable=QtGui.QCheckBox(self)
		self.wenable.setChecked(False)			# disable new processors by default to permit their values to be set
		self.gbl.addWidget(self.wenable,0,1)

		# List of processor categories
		self.wcat=QtGui.QComboBox(self)
		self.wcat.addItem("")
		for i in self.cats: self.wcat.addItem(i)
#		self.wcat.setCurrentindex(self.wcat.findText("processor"))
		self.gbl.addWidget(self.wcat,0,2)

		# List of processor subcategories
		self.wsubcat=QtGui.QComboBox(self)
		self.gbl.addWidget(self.wsubcat,0,3)
#		self.update_subcat()

		#button grid
		self.gbl2=QtGui.QGridLayout()
		self.gbl.addLayout(self.gbl2,0,0,1,1)
		if get_platform().lower()=="darwin": self.gbl2.setSpacing(10)
		else: self.gbl2.setSpacing(1)

#		self.gbl2.setColumnStretch(0,1)
#		self.gbl2.setColumnStretch(1,1)

#		self.wup = QtGui.QPushButton(app.style().standardIcon(QtGui.QStyle.SP_ArrowUp),"")
		self.wup = QtGui.QPushButton(pix_up,"",self)
		self.wup.setMaximumSize(QtCore.QSize(17, 17))
#		self.wup.setSizePolicy(QtGui.QSizePolicy.Fixed,QtGui.QSizePolicy.Fixed)
		self.gbl2.addWidget(self.wup,0,0)

		self.wdown = QtGui.QPushButton(pix_down,"",self)
		self.wdown.setMaximumSize(QtCore.QSize(17, 17))
		self.gbl2.addWidget(self.wdown,1,0)

		self.wplus = QtGui.QPushButton(pix_plus,"",self)
		self.wplus.setMaximumSize(QtCore.QSize(17, 17))
		self.gbl2.addWidget(self.wplus,1,1)

		self.wminus= QtGui.QPushButton(pix_minus,"",self)
		self.wminus.setMaximumSize(QtCore.QSize(17, 17))
		self.gbl2.addWidget(self.wminus,0,1)


		QtCore.QObject.connect(self.wcat,QtCore.SIGNAL("currentIndexChanged(int)"),self.eventCatSel)
		QtCore.QObject.connect(self.wsubcat,QtCore.SIGNAL("currentIndexChanged(int)"),self.eventSubcatSel)
		QtCore.QObject.connect(self.wup,QtCore.SIGNAL("clicked(bool)"),self.butUp)
		QtCore.QObject.connect(self.wdown,QtCore.SIGNAL("clicked(bool)"),self.butDown)
		QtCore.QObject.connect(self.wplus,QtCore.SIGNAL("clicked(bool)"),self.butPlus)
		QtCore.QObject.connect(self.wminus,QtCore.SIGNAL("clicked(bool)"),self.butminus)
		QtCore.QObject.connect(self.wenable,QtCore.SIGNAL("clicked(bool)"),self.updateFilt)

		self.parmw=[]

		self.setTag(tag)
#		print "Alloc processor ",self.tag

#	def __del__(self):
#		print "Free processor ",self.tag
#		QtGui.QWidget.__del__(self)

	def __getstate__(self):
		"used when pickling"
		proc=self.processorName()
		if not self.plist.has_key(proc) : return None		# invalid processor selection
		if self.wenable.isChecked() : proc=(proc,True)		# disabled, so we return None
		else: proc=(proc,False)

		parms={}
		for w in self.parmw:
			parms[w.getLabel()]=(w.getValue(),w.getEnabled())

		return (proc,parms)

# TODO - incomplete
	def __setstate__(self,state):
		"used when unpickling"
		self.init()
		if state==None : return

		proc,parms=state
		proc,enabled=proc

		try: proc_cat,proc_scat=proc.split(".",1)
		except:
			proc_cat=proc
			proc_scat="---"

		self.wcat.setCurrentIndex(self.wcat.findText(proc_cat))
		self.wsubcat.setCurrentIndex(self.wsubcat.findText(proc_scat))

	def getAsProc(self):
		"Returns the currently defined processor as a --process string for use in commands"

		if not self.wenable.isChecked() : return ""

		proc=self.processorName()
		if not self.plist.has_key(proc) : return ""		# invalid processor selection

		enabled=[]
		for w in self.parmw:
			if w.getEnabled() : enabled.append("%s=%s"%(w.getLabel(),str(w.getValue())))

		return "--process=%s:%s "%(self.processorName(),":".join(enabled))

	def getAsText(self):
		"Returns the currently defined processor as a 3 line string for persistence"

		proc=self.processorName()
		if not self.plist.has_key(proc) : return None		# invalid processor selection

		if self.wenable.isChecked() : ret="#$ enabled\n"
		else: ret="#$ disabled\n"

		enabled=[]
		disabled=[]
		for w in self.parmw:
			if w.getEnabled() : enabled.append("%s=%s"%(w.getLabel(),str(w.getValue())))
			else : disabled.append("%s=%s"%(w.getLabel(),str(w.getValue())))

		ret+="# %s\n--process=%s:%s\n\n"%(":".join(disabled),self.processorName(),":".join(enabled))

		return ret


	def setFromText(self,text):
		"""Sets the GUI from the text written to disk for persistent storage
		text should contain a 3-tuple of lines
		#$ enabled or disabled
		# name=value name=value ...
		--process=processor:name=value:name=value:..."""

#		print "set"

		if len(text)!=3 :
			raise Exception,"setFromText requires 3-tuple of strings"
		if text[0][0]!="#" or text[1][0]!="#" or text[2][0]!="-" :
			raise Exception,"Problem unpacking '%s' from file"%text

		disabled=parsemodopt("X:"+text[1][1:].strip())[1]			# dictionary of disabled values
		proc,enabled=parsemodopt(text[2].split("=",1)[1])	# dictionary of enabled values

		try: proc_cat,proc_scat=proc.split(".",1)
		except:
			proc_cat=proc
			proc_scat="---"

		self.wcat.setCurrentIndex(self.wcat.findText(proc_cat))
#		self.eventCatSel(0)
		self.wsubcat.setCurrentIndex(self.wsubcat.findText(proc_scat))
#		self.eventSubcatSel(0)


		for w in self.parmw:
			lbl=w.getLabel()
			if lbl in enabled:
				w.setValue(enabled[lbl])
				w.setEnabled(1)
			elif lbl in disabled:
				w.setValue(disabled[lbl])
				w.setEnabled(0)

		if "enabled" in text[0] : self.wenable.setChecked(True)
		else : self.wenable.setChecked(False)

#		print proc_cat,proc_scat,enabled,disabled
		return

	def setTag(self,tag):
		self.tag=tag

	def butUp(self):
		self.emit(QtCore.SIGNAL("upPress"),self.tag)

	def butDown(self):
		self.emit(QtCore.SIGNAL("downPress"),self.tag)

	def butPlus(self):
		self.emit(QtCore.SIGNAL("plusPress"),self.tag)

	def butminus(self):
		self.emit(QtCore.SIGNAL("minusPress"),self.tag)


	def eventCatSel(self,idx):
		"Called when the user selects a processor category"
		cat=str(self.wcat.currentText())
		#print "catsel ",cat
		#traceback.print_stack(limit=2)
		scats=['.'.join(i.split('.',1)[1:]) for i in self.plist if i.split(".")[0]==cat]
		scats.sort()
		for i in xrange(len(scats)):
			if scats[i]=="" : scats[i]="---"
		self.wsubcat.clear()
		self.wsubcat.addItem("")
		if len(scats)>0 : self.wsubcat.addItems(scats)

	def eventSubcatSel(self,idx):
		"Called when the user selects a processor subcategory. Updates the parameter widget set."
		proc=self.processorName()
		#print "subcatsel ",proc
		#traceback.print_stack(limit=2)

		try: parms=self.plist[proc]
		except: parms=["Please select a processor"]

		self.wsubcat.setToolTip(parms[0])

		# Remove old parameter widgets
		for w in self.parmw:
			self.gbl.removeWidget(w)
			w.close()				# get it off the screen
			w.setParent(None)		# do this before we can actually delete the widget
			del(w)					# finally, get it deleted (may not take effect until parmw=[] below

		# Iterate over the parameters
		# parms[i] - name of the parameter
		# parms[i+1] - type of parameter
		# parms[i+2] - helpstring for parameter
		self.parmw=[]
		self.ninput=0
		for i in range(1,len(parms),3):
			self.ninput+=1
			try: dflt=self.parmdefault[parms[i]]		# contains (start enabled, range, default value, change callback)
			except: dflt=(1,(0,5.0),None,None)			# default parameter settings if we don't have anything predefined

			if parms[i+1] in ("FLOAT","INT"):
				self.parmw.append(ValSlider(self,dflt[1],parms[i],dflt[2],100,dflt[0]))
				if parms[i+1]=="INT" :
					self.parmw[-1].setIntonly(1)
	#			self.parmw[-1].hboxlayout.setContentsMargins ( 11.0,5.0,5.0,5.0 )
			elif parms[i+1]=="BOOL" :
				self.parmw.append(CheckBox(self,parms[i],dflt[2],100,dflt[0]))

			elif parms[i+1]=="STRING" :
				self.parmw.append(StringBox(self,parms[i],dflt[2],100,dflt[0]))

			elif parms[i+1]=="TRANSFORM" :
				self.parmw.append(StringBox(self,parms[i],dflt[2],100,dflt[0]))

			elif parms[i+1]=="EMDATA" :
				self.parmw.append(StringBox(self,parms[i],dflt[2],100,dflt[0]))

			elif parms[i+1]=="FLOATARRAY" :
				self.parmw.append(StringBox(self,parms[i],dflt[2],100,dflt[0]))

			elif parms[i+1]=="XYDATA" :
				self.parmw.append(StringBox(self,parms[i],dflt[2],100,dflt[0]))

			else: print "Unknown parameter type",parms[i+1],parms

			self.parmw[-1].setToolTip(parms[i+2])
			self.gbl.addWidget(self.parmw[-1],self.ninput,1,1,4)
			QtCore.QObject.connect(self.parmw[-1], QtCore.SIGNAL("valueChanged"), self.updateFilt)
			QtCore.QObject.connect(self.parmw[-1], QtCore.SIGNAL("enableChanged"), self.updateFilt)

		self.updateFilt()

	def updateFilt(self,val=None):
		"Called whenever the processor changes"
		#if self.wenable.isChecked() :
		self.emit(QtCore.SIGNAL("processorChanged"),self.tag)

	def processorName(self):
		"Returns the name of the currently selected processor"
		cat=str(self.wcat.currentText())
		scat=str(self.wsubcat.currentText())
		if scat=="" :  return cat
		else : return cat+"."+scat

	def processorParms(self):
		"Returns the currently defined processor as (name,dict) or None if not set to a valid state"
		if not self.wenable.isChecked() : return None		# disabled, so we return None

		proc=self.processorName()
		if not self.plist.has_key(proc) : return None		# invalid processor selection


		parms={}
		for w in self.parmw:
			if w.getEnabled() : parms[w.getLabel()]=w.getValue()

		return (proc,parms)

class EMFilterTool(QtGui.QMainWindow):
	"""This class represents the EMFilterTool application instance.  """

	def __init__(self,datafile=None,apix=0.0,force2d=False,verbose=0):
		QtGui.QMainWindow.__init__(self)

		app=QtGui.qApp
		self.apix=apix
		self.force2d=force2d
		self.setWindowTitle("e2filtertool.py")

		# Menu Bar
		self.mfile=self.menuBar().addMenu("File")
#		self.mfile_save_processor=self.mfile.addAction("Save Processor Param")
		self.mfile_save_stack=self.mfile.addAction("Save Processed Stack")
		self.mfile_save_map=self.mfile.addAction("Save Processed Map")
		self.mfile_quit=self.mfile.addAction("Quit")

		self.setCentralWidget(QtGui.QWidget())
		self.vblm = QtGui.QVBoxLayout(self.centralWidget())		# The contents of the main window

		# List of processor sets
		self.wsetname=QtGui.QComboBox()
		self.wsetname.setEditable(True)
		psetnames=[i.split("_",1)[1][:-4].replace("_"," ") for i in os.listdir(".") if i[:11]=="filtertool_"]
		try: psetnames.remove("default")  # remove default if it exists
		except: pass
		psetnames.insert(0,"default")     # add it back at the top of the list
		for i in psetnames : self.wsetname.addItem(i)
		self.vblm.addWidget(self.wsetname)

		# scrollarea for processor widget
		self.processorsa=QtGui.QScrollArea()
		self.vblm.addWidget(self.processorsa)

		# Actual widget contianing processors being scrolled
		self.processorpanel=QtGui.QWidget()
		self.processorsa.setWidget(self.processorpanel)
		self.processorsa.setWidgetResizable(True)
		self.vbl = QtGui.QVBoxLayout(self.processorpanel)

		self.processorlist=[]
		self.addProcessor()

		# file menu
#		QtCore.QObject.connect(self.mfile_save_processor,QtCore.SIGNAL("triggered(bool)")  ,self.menu_file_save_processor  )
		QtCore.QObject.connect(self.mfile_save_stack,QtCore.SIGNAL("triggered(bool)")  ,self.menu_file_save_stack  )
		QtCore.QObject.connect(self.mfile_save_map,QtCore.SIGNAL("triggered(bool)")  ,self.menu_file_save_map  )
		QtCore.QObject.connect(self.mfile_quit,QtCore.SIGNAL("triggered(bool)")  ,self.menu_file_quit)

		QtCore.QObject.connect(self.wsetname,QtCore.SIGNAL("currentIndexChanged(int)"),self.setChange)


		self.viewer=None			# viewer window for data
		self.processors=[]			# list of processor specifications (name,dict)
		self.origdata=[]
		self.filtdata=[]
		self.nx=0
		self.ny=0
		self.nz=0
		self.busy=False		# used to prevent multiple updates during restore
		self.oldset=None	# currently selected processor set name

		if datafile!=None : self.setData(datafile)

		self.needupdate=1
		self.needredisp=0
		self.lastredisp=1
		self.procthread=None
		self.errors=None		# used to communicate errors back from the reprocessing thread

		self.restore_processorset("default")

		self.timer=QTimer()
		QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.timeOut)
		self.timer.start(100)
		E2loadappwin("e2filtertool","main",self)

#		QtCore.QObject.connect(self.boxesviewer,QtCore.SIGNAL("mx_image_selected"),self.img_selected)


	def setChange(self,line):
		"""When the user selects a new set or hits enter after a new name"""
		#print "setchange ",line
		newset=str(self.wsetname.currentText())

		if newset==self.oldset : return

		self.save_current_processorset(self.oldset)
		self.restore_processorset(newset)

	def addProcessor(self,after=-1):
		#print "addProc ",after,len(self.processorlist)
		after+=1
		if after==0 : after=len(self.processorlist)
		epw=EMProcessorWidget(self.processorpanel,tag=after)
		self.processorlist.insert(after,epw)
		self.vbl.insertWidget(after,epw)
		QtCore.QObject.connect(epw,QtCore.SIGNAL("upPress"),self.upPress)
		QtCore.QObject.connect(epw,QtCore.SIGNAL("downPress"),self.downPress)
		QtCore.QObject.connect(epw,QtCore.SIGNAL("plusPress"),self.plusPress)
		QtCore.QObject.connect(epw,QtCore.SIGNAL("minusPress"),self.minusPress)
		QtCore.QObject.connect(epw,QtCore.SIGNAL("processorChanged"),self.procChange)

		# Make sure all the tags are correct
		for i in range(len(self.processorlist)): self.processorlist[i].setTag(i)

	def delProcessor(self,idx):
		#print "delProc ",idx,len(self.processorlist)
		self.vbl.removeWidget(self.processorlist[idx])
		self.processorlist[idx].close()
		del self.processorlist[idx]

		# Make sure all the tags are correct
		for i in range(len(self.processorlist)): self.processorlist[i].setTag(i)

	def swapProcessors(self,tag):
		"This swaps 2 adjacent processors, tag and tag+1"
		w=self.processorlist[tag]
		self.processorlist[tag-1],self.processorlist[tag]=self.processorlist[tag],self.processorlist[tag-1]
		self.vbl.removeWidget(w)
		self.vbl.insertWidget(tag-1,w)

		# Make sure all the tags are correct
		for i in range(len(self.processorlist)): self.processorlist[i].setTag(i)

	def upPress(self,tag):
		if tag==0 : return
		self.swapProcessors(tag)

	def downPress(self,tag):
		if tag==len(self.processorlist)-1 : return

		self.swapProcessors(tag+1)

	def plusPress(self,tag):
		self.addProcessor(tag)

	def minusPress(self,tag):
		if len(self.processorlist)==1 : return		# Can't delete the last processor
		self.delProcessor(tag)

	def timeOut(self):
		if self.busy : return

		# Spawn a thread to reprocess the data
		if self.needupdate and self.procthread==None:
			self.procthread=threading.Thread(target=self.reprocess)
			self.procthread.start()

		if self.errors:
			QtGui.QMessageBox.warning(None,"Error","The following processors encountered errors during processing of 1 or more images:"+"\n".join(self.errors))
			self.errors=None

		# When reprocessing is done, we want to redisplay from the main thread
		if self.needredisp :
			self.needredisp=0
			self.viewer.show()
			if self.nz==1 or self.force2d: self.viewer.set_data(self.procdata)
			else :
				self.sgdata.setData(self.procdata[0])
				self.viewer.updateSG()

	def procChange(self,tag):
#		print "change"
		self.needupdate=1

	def reprocess(self):
		"Called when something about a processor changes (or when the data changes)"

		if self.busy: return
		self.needupdate=0		# we set this immediately so we reprocess again if an update happens while we're processing
		self.procdata=[im.copy() for im in self.origdata]

		needredisp=0
		errors=[]
		for p in self.processorlist:
			pp=p.processorParms()				# processor parameters
			if pp==None: continue				# disabled processor
			proc=Processors.get(pp[0],pp[1])	# the actual processor object
#			print pp

			errflag=False
			for im in self.procdata:
				try: proc.process_inplace(im)
				except: errflag=True

			if errflag: errors.append(str(pp))
			needredisp=1					# if all processors are disabled, we don't want a redisplay
			self.lastredisp=1

		if len(errors)>0 :
			self.errors=errors

		self.procthread=None					# This is a reference to ourselves (the thread doing the processing). we reset it ourselves before returning
		self.needredisp=max(needredisp,self.lastredisp)
		self.lastredisp=needredisp

	def setData(self,data):
		if data==None :
			self.data=None
			return

		elif isinstance(data,str) :
			self.datafile=data
			self.nimg=EMUtil.get_image_count(data)

			self.origdata=EMData(data,0)

			if self.origdata["nz"]==1:
				if self.nimg>20 :
					self.origdata=EMData.read_images(data,range(0,self.nimg,self.nimg/20))		# read regularly separated images from the file totalling ~20
				elif self.nimg>1 :
					self.origdata=EMData.read_images(data,range(self.nimg))
				else: self.origdata=[self.origdata]
			else :
				self.origdata=[self.origdata]

		else :
			self.datafile=None
			if isinstance(data,EMData) : self.origdata=[data]
			else : self.origdata=data

		self.nx=self.origdata[0]["nx"]
		self.ny=self.origdata[0]["ny"]
		self.nz=self.origdata[0]["nz"]
		if self.apix<=0.0 : self.apix=self.origdata[0]["apix_x"]
		EMProcessorWidget.parmdefault["apix"]=(0,(0.2,10.0),self.apix,None)

		if self.viewer!=None : self.viewer.close()
		if self.nz==1 or self.force2d:
			if len(self.origdata)>1 :
				self.viewer=EMImageMXWidget()
				self.mfile_save_stack.setEnabled(True)
				self.mfile_save_map.setEnabled(False)
			else :
				self.viewer=EMImage2DWidget()
				self.mfile_save_stack.setEnabled(False)
				self.mfile_save_map.setEnabled(True)
		else :
			self.mfile_save_stack.setEnabled(False)
			self.mfile_save_map.setEnabled(True)
			self.viewer = EMScene3D()
			self.sgdata = EMDataItem3D(test_image_3d(3), transform=Transform())
			isosurface = EMIsosurface(self.sgdata, transform=Transform())
			self.viewer.insertNewNode('Data', self.sgdata, parentnode=self.viewer)
			self.viewer.insertNewNode("Iso", isosurface, parentnode=self.sgdata)

		E2loadappwin("e2filtertool","image",self.viewer.qt_parent)

		self.procChange(-1)

	def save_current_processorset(self,name):
		"""Saves the current processor and parameters to a text file"""
#		print "saveset ",name

		try: out=file("filtertool_%s.txt"%(name.replace(" ","_")),"w")	# overwrite the current contents
		except:
			traceback.print_exc()
			print "No permission to store processorset info"
			return

		out.write("# This file contains the parameters for the processor set named\n# %s\n# Each of the --process lines below is in the correct syntax for use with e2proc2d.py or e2proc3d.py.\n# Use the full set sequentially in a single command to replicate the processor set\n\n"%name)

		for filt in self.processorlist:
			ftext=filt.getAsText()
			if ftext!=None : out.write(ftext)

		out.close()

	def restore_processorset(self,name):
		"""This restores a processorset that has been written to a text file"""
		# Erase all existing processors
#		print "loadset ",name

		self.oldset=name

		self.busy=True
		while len(self.processorlist)>0 : self.delProcessor(0)

		# Open the file
		try: infile=file("filtertool_%s.txt"%(name.replace(" ","_")),"r")
		except:
			self.addProcessor()
			self.busy=False
			return		# file can't be read, must be a new set (or we'll assume that)

		while 1:
			l=infile.readline()
			if l=="" : break
			if l[:2]=="#$" :								# This indicates the start of a new processor
				l=[l.strip(),infile.readline().strip(),infile.readline().strip()]
				self.addProcessor()
				self.processorlist[-1].setFromText(l)

		if len(self.processorlist)==0 : self.addProcessor()

		self.busy=False
		self.needupdate=1

	def menu_file_save_processor(self):
		"Saves the processor in a usable form to a text file"
		#out=file("processor.txt","a")
		#out.write("Below you will find the processor options in sequence. They can be passed together in order\nas options to a single e2proc2d.py or e2proc3d.py command to achieve the\nsame results as in e2filtertool.py\n")
		#for p in self.processorlist:
			#pp=p.processorParms()				# processor parameters
			#if pp==None : continue
			#out.write("--process=%s"%pp[0])
			#for k in pp[1]:
				#out.write(":%s=%s"%(k,str(pp[1][k])))
			#out.write("\n")
		#out.write("\n")
		#out.close()

		#QtGui.QMessageBox.warning(None,"Saved","The processor parameters have been added to the end of 'processor.txt'")

		self.save_current_processorset(str(self.wsetname.currentText()))

	def menu_file_save_stack(self):
		"Processes the entire current stack, and saves as a new name"

		name=QtGui.QInputDialog.getText(None,"Enter Filename","Enter an output filename for the entire processed particle stack (not just the displayed images).")
		if not name[1] : return		# canceled

		allfilt=" ".join([i.getAsProc() for i in self.processorlist])

		n=EMUtil.get_image_count(self.datafile)
		from PyQt4.QtGui import QProgressDialog
		progressdialog=QProgressDialog("Processing Images","Abort",0,n,self)
		progressdialog.setMinimumDuration(1000)

		e=E2init(["e2proc2d.py",self.datafile,str(name[0]),allfilt])	# we don't actually run this program, since we couldn't have a progress dialog easo;y then

		pp=[i.processorParms() for i in self.processorlist]

		for i in xrange(n):
			im=EMData(self.datafile,i)
			QtGui.qApp.processEvents()
			for p in pp: im.process_inplace(p[0],p[1])
			im.write_image(str(name[0]),i)
			progressdialog.setValue(i+1)
			if progressdialog.wasCanceled() :
				print "Processing Cancelled"
				break

		progressdialog.setValue(n)

		E2end(e)

	def menu_file_save_map(self):
		"saves the current processored map"

		if len(self.procdata)==1 and self.procdata[0]["nz"]>1 :
			try: os.unlink("processed_map.hdf")
			except : pass
			self.procdata[0].write_image("processed_map.hdf",0)
			QtGui.QMessageBox.warning(None,"Saved","The processed map has been saved as processed_map.hdf")
		else :
			try: os.unlink("processed_images.hdf")
			except: pass
			for i in self.procdata: i.write_image("processed_images.hdf",-1)
			QtGui.QMessageBox.warning(None,"Saved","The processed image(s) has been saved as processed_images.hdf. WARNING: this will include only be a subset of the images in a large image stack. To process the full stack, use e2proc2d.py with the options in filtertool_<filtername>.txt")

	def menu_file_quit(self):
		self.close()

	def closeEvent(self,event):
		E2saveappwin("e2filtertool","main",self)
		self.save_current_processorset(str(self.wsetname.currentText()))

#		print "Exiting"
		if self.viewer!=None :
			E2saveappwin("e2filtertool","image",self.viewer.qt_parent)
			self.viewer.close()
		event.accept()
		#self.app().close_specific(self)
		self.emit(QtCore.SIGNAL("module_closed")) # this signal is important when e2ctf is being used by a program running its own event loop

	#def closeEvent(self,event):
		#self.target().done()


pix_plus=None
pix_up=None
pix_down=None
pix_minus=None

def pix_init():
	global pix_plus,pix_minus,pix_up,pix_down

	pix_plus=QtGui.QIcon(QtGui.QPixmap(["15 15 3 1",
	" 	c None",
	".	c black",
	"X	c grey",
	"               ",
	"               ",
	"               ",
	"               ",
	"       .       ",
	"       .X      ",
	"       .X      ",
	"    .......    ",
	"     XX.XXXX   ",
	"       .X      ",
	"       .X      ",
	"        X      ",
	"               ",
	"               ",
	"               "
	]))

	pix_up=QtGui.QIcon(QtGui.QPixmap(["15 15 3 1",
	" 	c None",
	".	c black",
	"X	c grey",
	"               ",
	"               ",
	"               ",
	"       .       ",
	"      ...      ",
	"     .....     ",
	"    .......    ",
	"    XXX.XXXX   ",
	"       .X      ",
	"       .X      ",
	"       .X      ",
	"       .X      ",
	"        X      ",
	"               ",
	"               "
	]))

	pix_down=QtGui.QIcon(QtGui.QPixmap(["15 15 3 1",
	" 	c None",
	".	c black",
	"X	c grey",
	"               ",
	"               ",
	"               ",
	"       .       ",
	"       .X      ",
	"       .X      ",
	"       .X      ",
	"       .X      ",
	"    .......    ",
	"     .....X    ",
	"      ...X     ",
	"       .X      ",
	"               ",
	"               ",
	"               "
	]))

	pix_minus=QtGui.QIcon(QtGui.QPixmap(["15 15 3 1",
	" 	c None",
	".	c black",
	"X	c grey",
	"               ",
	"               ",
	"               ",
	"               ",
	"               ",
	"               ",
	"               ",
	"    .......    ",
	"     XXXXXXX   ",
	"               ",
	"               ",
	"               ",
	"               ",
	"               ",
	"               "
	]))

if __name__ == "__main__":
	main()



