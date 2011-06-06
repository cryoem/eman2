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

from optparse import OptionParser
import sys
import os
import weakref
import threading

from EMAN2 import *
from emapplication import get_application, EMApp
from emimage2d import EMImage2DWidget
from emimagemx import EMImageMXWidget
from emimage3d import EMImage3DWidget
from emshape import EMShape
from valslider import *

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <image file>

	Provides a GUI interface for applying a sequence of processors to an image, stack of images or a volume."""

	parser = OptionParser(usage=usage,version=EMANVERSION)

#	parser.add_option("--boxsize","-B",type="int",help="Box size in pixels",default=64)
#	parser.add_option("--shrink",type="int",help="Shrink factor for full-frame view, default=0 (auto)",default=0)
	parser.add_option("--apix",type="float",help="Override the A/pix value stored in the file header",default=0.0)
	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n", type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()
		
	if len(args) != 1: parser.error("You must specify a single data file on the command-line.")
	if not file_exists(args[0]): parser.error("%s does not exist" %args[0])
#	if options.boxsize < 2: parser.error("The boxsize you specified is too small")
#	# The program will not run very rapidly at such large box sizes anyhow
#	if options.boxsize > 2048: parser.error("The boxsize you specified is too large.\nCurrently there is a hard coded max which is 2048.\nPlease contact developers if this is a problem.")
	
#	logid=E2init(sys.argv)
	
	app = EMApp()
	pix_init()
	control=EMProcessorTool(datafile=args[0],apix=options.apix,verbose=options.verbose)
	control.show()
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
		self.wenable.setChecked(True)
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
		self.gbl2.setSpacing(1)
		
		self.gbl2.setColumnStretch(0,1)
		self.gbl2.setColumnStretch(1,1)
		
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
#		print "Alloc filter ",self.tag
	
#	def __del__(self):
#		print "Free filter ",self.tag
#		QtGui.QWidget.__del__(self)

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
		scats=[i.split('.',1)[1] for i in self.plist if i.split(".")[0]==cat]
		scats.sort()
		self.wsubcat.clear()
		self.wsubcat.addItem("")
		if len(scats)>0 : self.wsubcat.addItems(scats)

	def eventSubcatSel(self,idx):
		"Called when the user selects a processor subcategory. Updates the parameter widget set."
		proc=self.processorName()
		
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

			else: print "Unknown parameter type",parms[i+1],parms
			
			self.parmw[-1].setToolTip(parms[i+2])
			self.gbl.addWidget(self.parmw[-1],self.ninput,1,1,4)
			QtCore.QObject.connect(self.parmw[-1], QtCore.SIGNAL("valueChanged"), self.updateFilt)
			QtCore.QObject.connect(self.parmw[-1], QtCore.SIGNAL("enableChanged"), self.updateFilt)
		
		self.updateFilt()
	
	def updateFilt(self,val=None):
		"Called whenever the processor changes"
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
	
class EMProcessorTool(QtGui.QMainWindow):
	"""This class represents the EMTomoBoxer application instance.  """
	
	def __init__(self,datafile=None,apix=0.0,verbose=0):
		QtGui.QMainWindow.__init__(self)
		
		app=QtGui.qApp
		self.apix=apix
		self.setWindowTitle("e2processortool.py")
		
		# Menu Bar
		self.mfile=self.menuBar().addMenu("File")
		self.mfile_save_processed=self.mfile.addAction("Save processed data")
		self.mfile_quit=self.mfile.addAction("Quit")

		self.setCentralWidget(QtGui.QScrollArea())
		self.filterpanel=QtGui.QWidget()
		self.centralWidget().setWidget(self.filterpanel)
		self.centralWidget().setWidgetResizable(True)
		self.vbl = QtGui.QVBoxLayout(self.filterpanel)

		self.filterlist=[]
		self.addFilter()
		
		# file menu
		QtCore.QObject.connect(self.mfile_save_processed,QtCore.SIGNAL("triggered(bool)")  ,self.menu_file_save_processed  )
		QtCore.QObject.connect(self.mfile_quit,QtCore.SIGNAL("triggered(bool)")  ,self.menu_file_quit)

		self.viewer=None			# viewer window for data
		self.processors=[]				# list of processor specifications (name,dict)
		self.origdata=[]
		self.filtdata=[]
		self.nx=0
		self.ny=0
		self.nz=0

		if datafile!=None : self.setData(datafile)
		
		self.needupdate=1
		self.needredisp=0
		self.procthread=None
		
		self.timer=QTimer()
		QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.timeOut)
		self.timer.start(0.1)
#		QtCore.QObject.connect(self.boxesviewer,QtCore.SIGNAL("mx_image_selected"),self.img_selected)

	def addFilter(self,after=-1):
		after+=1
		if after==0 : after=len(self.filterlist)
		epw=EMProcessorWidget(self.filterpanel,tag=after)
		self.filterlist.insert(after,epw)
		self.vbl.insertWidget(after,epw)
		QtCore.QObject.connect(epw,QtCore.SIGNAL("upPress"),self.upPress)
		QtCore.QObject.connect(epw,QtCore.SIGNAL("downPress"),self.downPress)
		QtCore.QObject.connect(epw,QtCore.SIGNAL("plusPress"),self.plusPress)
		QtCore.QObject.connect(epw,QtCore.SIGNAL("minusPress"),self.minusPress)
		QtCore.QObject.connect(epw,QtCore.SIGNAL("processorChanged"),self.procChange)

		# Make sure all the tags are correct
		for i in range(len(self.filterlist)): self.filterlist[i].setTag(i)

	def delFilter(self,idx):
		self.vbl.removeWidget(self.filterlist[idx])
		del self.filterlist[idx]
		
		# Make sure all the tags are correct
		for i in range(len(self.filterlist)): self.filterlist[i].setTag(i)

	def swapFilters(self,tag):
		"This swaps 2 adjacent filters, tag and tag+1"
		w=self.filterlist[tag]
		self.filterlist[tag-1:tag+1]=self.filterlist[tag:tag-2:-1]		#swap the entries in the list
		w=self.vbl.removeWidget(w)
		self.vbl.insertWidget(tag-1,w)
		
		# Make sure all the tags are correct
		for i in range(len(self.filterlist)): self.filterlist[i].setTag(i)

	def upPress(self,tag):
		if tag==0 : return
		self.swapfilters(tag)
	
	def downPress(self,tag):
		if tag==len(self.filterlist)-1 : return

		self.swapFilters(tag+1)
	
	def plusPress(self,tag):
		self.addFilter(tag)
		
	def minusPress(self,tag):
		if len(self.filterlist)==1 : return		# Can't delete the last filter
		self.delFilter(tag)
		
	def timeOut(self):
		# Spawn a thread to reprocess the data
		if self.needupdate and self.procthread==None: 
			self.procthread=threading.Thread(target=self.reprocess)
			self.procthread.start()
		
		# When reprocessing is done, we want to redisplay from the main thread
		if self.needredisp :
			self.needredisp=0
			self.viewer.show()
			if self.nz==1: self.viewer.set_data(self.procdata)
			else : self.viewer.set_data(self.procdata[0])
		
		
	def procChange(self,tag):
		self.needupdate=1
		
	def reprocess(self):
		"Called when something about a processor changes (or when the data changes)"
		
		self.needupdate=0		# we set this immediately so we reprocess again if an update happens while we're processing
		self.procdata=[im.copy() for im in self.origdata]
		
		for p in self.filterlist:
			pp=p.processorParms()				# processor parameters
			if pp==None: continue				# disabled processor
			proc=Processors.get(pp[0],pp[1])	# the actual processor object
			
			for im in self.procdata:
				proc.process_inplace(im)
				
		self.needredisp=1
		self.procthread=None					# This is a reference to ourselves (the thread doing the processing). we reset it ourselves before returning
		
	
	def setData(self,data):
		if data==None :
			self.data=None
			return
		
		elif isinstance(data,str) :
			self.datafile=data
			self.nimg=EMUtil.get_image_count(data)
			
			self.origdata=EMData(data,0)
			
			if self.nz==1:
				if self.n>20 : 
					self.origdata=EMData.read_images(data,range(20))
				elif self.n>1 :
					self.origdata=EMData.read_images(data,range(self.n))
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

		if self.viewer!=None : self.viewer.close()
		if self.nz==1 : self.viewer=EMImageMXWidget()
		else : self.viewer=EMImage3DWidget()
		
		self.procChange(-1)

	def menu_file_save_processed(self):
		"incomplete"
		
	def menu_file_quit(self):
		self.close()
		
	def closeEvent(self,event):
		print "Exiting"
		if self.viewer!=None : self.viewer.close()
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
		
		
		
