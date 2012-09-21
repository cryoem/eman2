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
import threading
import time
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
#		os.makedirs(options.path)

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

		self.wbdrawali=QtGui.QPushButton("Draw")
		self.hbl1.addWidget(self.wbdrawali)
		self.wbdrawali.hide()						# this functionality won't work with the current widget

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
		self.wsbalimaskbase.setValue(10)
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
		
		self.wbrecalcref=QtGui.QPushButton("Realign")
		self.hbl1a.addWidget(self.wbrecalcref)
		
		self.wbrrecalcref=QtGui.QPushButton("Rerefine")
		self.hbl1a.addWidget(self.wbrrecalcref)
		
		self.hbl1a.addStretch(5)


		###### ROI Mask
		# widget for editing the ROI mask
		self.wlroimaskdraw=QtGui.QLabel("<big><pre>R\nO\nI</pre></big>")
		self.gbl.addWidget(self.wlroimaskdraw,4,0)
		
		self.w2droimaskdraw=EMImage2DWidget()
		self.gbl.addWidget(self.w2droimaskdraw,4,1)

		# Buttons for controlling mask
		self.hbl2=QtGui.QHBoxLayout()
		self.gbl.addLayout(self.hbl2,5,1)
		self.hbl2.addStretch(5)
		
		self.wbdrawroi=QtGui.QPushButton("Draw")
		self.hbl1.addWidget(self.wbdrawroi)
		self.wbdrawroi.hide()							# this button won't work right for now

		self.wbautoroi=QtGui.QPushButton("Auto")
		self.hbl2.addWidget(self.wbautoroi)
		
		self.wbresetroi=QtGui.QPushButton("Reset")
		self.hbl2.addWidget(self.wbresetroi)
		
		self.hbl2.addStretch(5)

		# Widget for setting alignment mask blur and base level
		self.vbl2=QtGui.QVBoxLayout()
		self.gbl.addLayout(self.vbl2,4,2)
		self.vbl2.addStretch(5)
		
		self.wlroimaskblur=QtGui.QLabel("Blur")
		self.vbl2.addWidget(self.wlroimaskblur)
		
		self.wsbroimaskblur=QtGui.QSpinBox()
		self.wsbroimaskblur.setRange(0,25)
		self.vbl2.addWidget(self.wsbroimaskblur)

		self.vbl2.addSpacing(16)
		
		self.wbroigo=QtGui.QPushButton(QtCore.QChar(0x2192))
		self.vbl2.addWidget(self.wbroigo)


		# widget for displaying the masked ROI
		self.w2droimask=EMImage2DWidget()
		self.gbl.addWidget(self.w2droimask,4,3)
		
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
		
		try :
			cores=num_cpus()
		except:
			cores=2
		if cores==1 : cores=2
		self.wvbcores=ValBox(None,(0,256),"# Threads",cores)
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
		self.wpbprogress.reset()
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

		self.w2dptcl=EMImage2DWidget()
		self.gbl.addWidget(self.w2dptcl,4,9)
		
		## Buttons for controlling mask
		#self.hbl1=QtGui.QHBoxLayout()
		#self.gbl.addLayout(self.hbl1,2,1)
		#self.hbl1.addStretch(5)
		
		#self.wbautoali=QtGui.QPushButton("Auto")
		#self.hbl1.addWidget(self.wbautoali)
		
		#self.wbresetali=QtGui.QPushButton("Reset")
		#self.hbl1.addWidget(self.wbresetali)"bdb:%s"

		#self.hbl1.addStretch(5)

		QtCore.QObject.connect(self.wbdrawali,QtCore.SIGNAL("clicked(bool)"),self.aliDrawMode)
		QtCore.QObject.connect(self.wbautoali,QtCore.SIGNAL("clicked(bool)"),self.aliAutoPress)
		QtCore.QObject.connect(self.wbresetali,QtCore.SIGNAL("clicked(bool)"),self.aliResetPress)
		QtCore.QObject.connect(self.wbaligo,QtCore.SIGNAL("clicked(bool)"),self.aliGoPress)
		QtCore.QObject.connect(self.wbrecalcref,QtCore.SIGNAL("clicked(bool)"),self.aliRecalcRefPress)
		QtCore.QObject.connect(self.wbrrecalcref,QtCore.SIGNAL("clicked(bool)"),self.aliRRecalcRefPress)
		QtCore.QObject.connect(self.wbdrawroi,QtCore.SIGNAL("clicked(bool)"),self.roiDrawMode)
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
		
		insp=self.w2droimaskdraw.get_inspector()
		insp.hide()
		insp.mmtab.setCurrentIndex(3)
		insp.dtpenv.setText("0.0")
		
		self.aliimg=None		# This is the unmasked alignment reference image
		self.alisig=None		# This is the standard deviation of the alignment reference
		self.alimask=None		# This is the drawn-upon version of aliimg used to generate the mask
		self.alimasked=None		# This is the masked alignment reference used for the actual alignments
		self.roimask=None		# Region of interest mask drawing image
		self.roimasked=None		# Region of interest display widget

		self.particles=None
		self.particles_ali=None

		self.path=path
		QtCore.QTimer.singleShot(200,self.afterStart)

	def afterStart(self):
		"""This gets called once the event loop is running"""
		print "App running, initializing"
		self.initPath(self.path)

	def initPath(self,path):
		if path[:4].lower()!="bdb:" : path="bdb:"+path
		self.path=path
		dicts=db_list_dicts(self.path)
		
		# If particles exists then we can fully initialize
		if db_check_dict("%s#particles"%self.path) :
			self.particles=EMData.read_images("%s#particles"%self.path)	# read in the entire particle stack
			for p in self.particles: p.process_inplace("normalize.edgemean")
			
			cl=[i for i in dicts if "classes" in i]
			cl.sort()
			self.wlistresult.addItems(QtCore.QStringList(cl))
			
			self.bootstrap()
		else :
			self.particles=None
			
		return

	def menuFileOpen(self,x):
		if self.particles!=None:
			QtGui.QMessageBox.warning(None,"Error","%s already contains a stack of particles. A new folder is required to start with a new stack of particles. Rerun without --path option."%self.path)
			return

		self.dialog = embrowser.EMBrowserWidget(withmodal=True,multiselect=False)
		QtCore.QObject.connect(self.dialog,QtCore.SIGNAL("ok"),self.gotPath)
		QtCore.QObject.connect(self.dialog,QtCore.SIGNAL("cancel"),self.gotPath)
		self.dialog.show()
	
	def gotPath(self):
		ptcl=self.dialog.getResult()
		self.dialog=None
		if ptcl == None : return		# user pressed cancel, but we still need to clean up
		ptcl=ptcl[0]
		
		n=EMUtil.get_image_count(ptcl)
		sz=EMData(ptcl,0,True)["nx"]
		
		# We don't want the user to overwhelm the system
		if sz*sz*4*n > 5.0e8 :
			nmax=5.0e8/(sz*sz*4)
			r=QtGui.QMessageBox.question(None,"Are you sure ?","WARNING: This full particle set will require %d+ gb memory to process. Select Yes to use only the first %d particles, No to use the entire stack, or Cancel to abort. You may also consider using e2proc2d.py --meanshrink= to produce a downsampled stack for processing."%(int(sz*sz*9*n/1.0e9+.5),nmax),QtGui.QMessageBox.Yes|QtGui.QMessageBox.No|QtGui.QMessageBox.Cancel)
			if r==QtGui.QMessageBox.Cancel : return
			if r==QtGui.QMessageBox.Yes : n=nmax
		
		if ptcl[:4].lower()=="bdb:" : task="e2bdb.py %s  --makevstack=%s#particles --step=0,1,%d"%(ptcl,self.path,n)
		else : task="e2proc2d.py %s %s#particles --last=%d"%(ptcl,self.path,n)
		print task
		launch_childprocess(task)
		
		self.initPath(self.path)
	
	def threadAlign(self,stack,ref,outstack):
		"""This method is designed to run in a thread and perform alignments of a stack of inputs to a single reference"""
	
		for p in stack:
			p2=p.process("filter.lowpass.gauss",{"cutoff_abs":0.1})
			a=p2.align("rotate_translate_flip",ref)
#			a=p2.copy()
			a["match_qual"]=a.cmp("frc",ref,{"sweight":0})		# compute similarity to unmasked reference
			outstack.append((a["match_qual"],a))

	def threadTreeAlign(self,tree,tree2,lock):
		"""Averages pairs in tree to produce 1/2 as many averages in tree2. Competes with other threads for elements to work on."""
		
		while len(tree)>0:
			# we need to pull a pair from 'tree'. We use a lock to make sure we don't end up with a bunch of singletons at the end
			lock.acquire()
			try:
				a=tree.pop()
			except:
				lock.release()
				break
			
			try: b=tree.pop()
			except:
				tree2.append(a)
				lock.release()
				break
			lock.release()
			
			# now we do the alignment and dump the result into the shared output list
			c=a.align("rotate_translate_flip",b)
			c.add(b)
			c.process_inplace("xform.centerofmass",{"threshold":0.5})
			tree2.append(c)

	def bootstrap(self):
		"""This will create an initial alignment reference for the particles when no reference is available"""

		self.wpbprogress.setEnabled(True)
		self.wpbprogress.reset()
		QtGui.qApp.processEvents()
		nthr=int(self.wvbcores.getValue())		# number of threads to use for faster alignments
		print nthr

		tree=[i.process("math.meanshrink",{"n":2}) for i in self.particles]
		for j,i in enumerate(tree): 
			self.wpbprogress.setValue(int(j*100/len(tree)))
			QtGui.qApp.processEvents()
			i.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.1})
		
		lock=threading.Lock()

		print "b1"
		while len(tree)>1:
			maxt=len(tree)
			tree2=[]
			
			# launch nthr threads to do the alignments
			thrs=[]
			for i in range(nthr):
				print i
				thrs.append(threading.Thread(target=self.threadTreeAlign, args=(tree,tree2,lock)))
				thrs[-1].start()
				
			# Wait for threads to finish
			while 1:
				time.sleep(0.2)
				print len(tree)
				self.wpbprogress.setValue(int((maxt-len(tree))*100/maxt))
				QtGui.qApp.processEvents()

				# If any threads are alive, it breaks out of the inner loop, if none are alive, the else block breaks out of the outer loop
				for t in thrs:
					if t.isAlive() : break
				else:
					break

			tree=tree2
		
#		tree[0]=tree[0].get_clip(Region(-tree[0]["nx"]/2,-tree[0]["ny"]/2,tree[0]["nx"]*2,tree[0]["ny"]*2))
		tree[0].process_inplace("xform.scale",{"clip":tree[0]["nx"]*2,"scale":2.0})
		tree[0].process_inplace("xform.centerofmass",{"threshold":0.5})
		tree[0].process_inplace("normalize.edgemean")

		# We look for the most 'featurefull' direction to put on Y
		f=tree[0].do_fft()
		f.ri2ap()
		pl=f.calc_az_dist(90,-90.0,2.0,3,16)
#		for i in range(90): print "%d\t%f"%(i,pl[i])
		ml=2.0*pl.index(max(pl))
		tree[0].rotate(ml,0,0)				# put the initial reference in the preferred orientation
		tree[0].write_image("zzz.hdf",0)
		
		# One final alignment pass
		self.particles_ali=[]
		self.alisig=EMData(self.particles[0]["nx"],self.particles[0]["ny"],1)
		avgr=Averagers.get("mean",{"sigma":self.alisig})
		thrs=[]
		
		print "b2"
		# launch nthr threads to do the alignments
		for i in range(nthr):
			print i
			thrs.append(threading.Thread(target=self.threadAlign, args=(self.particles[i::nthr],tree[0],self.particles_ali)))
			thrs[-1].start()

		# wait for the threads to finish running
		while (1):
			time.sleep(0.2)
			print len(self.particles_ali)
			self.wpbprogress.setValue(int(len(self.particles_ali)*100/len(self.particles)))
			QtGui.qApp.processEvents()
			
			# If any threads are alive, it breaks out of the inner loop, if none are alive, the else block breaks out of the outer loop
			for t in thrs:
				if t.isAlive() : break
			else:
				break
		
		print "b3"
		self.particles_ali.sort()
		for i,a in self.particles_ali[:len(self.particles)/2]: avgr.add_image(a)
		
		self.aliimg=avgr.finish()

		self.w2dptcl.set_data([i for j,i in self.particles_ali])
		

		self.aliimg.write_image("zzz.hdf",1)
		self.alisig.write_image("zzz.hdf",2)
		
		self.setAliRef(self.aliimg)
		
		self.wpbprogress.reset()
		self.wpbprogress.setEnabled(False)
	
	def setAliRef(self,img):
		"""Sets a new alignment reference, creating a new empty mask if none is present"""
		
		self.aliimg=img
		
		# If there is an existing mask drawn, we copy it out before overwriting
		mask=None
		if self.alimask!=None:
			mask=self.alimask
			mask.process_inplace("threshold.binary",{"value":0.001})
		
		self.alimask=img.copy()							# we need a copy of the pristine alignment reference to mask later
		self.alimask.add(-self.alimask["minimum"]+.01)		# The minimum value in the image is >0, so when we draw with 0 'color' we can extract the mask info
		if mask!=None : self.alimask.mult(mask)
		
		self.w2dalimaskdraw.set_data(self.alimask)
		self.aliGoPress()								# regenerate the actual alignment reference image (alimasked)
		
		rmask=None
		if self.roimask!=None:
			rmask=self.roimask
			rmask.process_inplace("threshold.binary",{"value":0.001})
			
		self.roimask=img.copy()
		self.roimask.add(-self.roimask["minimum"]+.01)		# The minimum value in the image is >0, so when we draw with 0 'color' we can extract the mask info
		if rmask!=None: self.roimask.mult(rmask)
		
		self.w2droimaskdraw.set_data(self.roimask)
		self.roiGoPress()								# regenerate the region of interest display

	def aliDrawMode(self,x=False):
		insp=self.w2dalimaskdraw.inspector
		if self.wbdrawali.text()=="Draw":
			insp.dtpenv.setText("0.0")

	def aliAutoPress(self,x=False):
		pass
	
	def aliResetPress(self,x=False):
		self.alimask=self.aliimg.copy()
		self.alimask.add(-self.alimask["minimum"]+.01)		# The minimum value in the image is >0, so when we draw with 0 'color' we can extract the mask info
		self.w2dalimaskdraw.set_data(self.alimask)
		self.aliGoPress()

	def aliGoPress(self,x=False):
		self.alimasked=self.aliimg.copy()
		
		base=self.wsbalimaskbase.value()
		blur=self.wsbalimaskblur.value()
		
		mask=self.alimask.process("threshold.binary",{"value":0.001})			# binarize the drawn mask
		mask.process_inplace("math.linear",{"scale":-1.0,"shift":1.0+base/100.0})		# invert the mask (user selects the region to include, not exclude)
		mask.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.5/(blur+.01)})
		
		self.alimasked.mult(mask)
		self.w2dalimask.set_data(self.alimasked)

	def rAlignToRef(self):
		"""realigns particles to reference, but does not compute a new reference"""
		self.wpbprogress.setEnabled(True)
		self.wpbprogress.reset()
		
		newali=[]
		for i,p in enumerate(self.particles):
			self.wpbprogress.setValue(int(i*100/len(self.particles)))
			QtGui.qApp.processEvents()
#			a=p.align("rotate_translate_flip",self.alimasked)
			p2=p.process("filter.lowpass.gauss",{"cutoff_abs":0.1})
			a=p2.align("refine",self.alimasked,{"verbose":0,"xform.align2d":self.particles_ali[i][1]["xform.align2d"]},"ccc",{})		# doesn't really seem to make an improvement
#			a["match_qual"]=a.cmp("ccc",self.aliimg)		# compute similarity to unmasked reference
			a["match_qual"]=a.cmp("frc",self.aliimg,{"sweight":0})		# compute similarity to unmasked reference
			newali.append((a["match_qual"],a))

		newali.sort()
		self.particles_ali=newali
		self.w2dptcl.set_data([i for j,i in self.particles_ali])

		self.wpbprogress.reset()
		self.wpbprogress.setEnabled(False)

	def alignToRef(self):
		"""realigns particles to reference, but does not compute a new reference"""
		self.wpbprogress.setEnabled(True)
		self.wpbprogress.reset()
		
		self.particles_ali=[]
		for i,p in enumerate(self.particles):
			self.wpbprogress.setValue(int(i*100/len(self.particles)))
			QtGui.qApp.processEvents()
			p2=p.process("filter.lowpass.gauss",{"cutoff_abs":0.1})
			a=p2.align("rotate_translate_flip",self.alimasked)
#			a=p.align("refine",tree[0],{"verbose":0,"xform.align2d":a["xform.align2d"]},"ccc",{})		# doesn't really seem to make an improvement
#			a["match_qual"]=a.cmp("ccc",self.aliimg)		# compute similarity to unmasked reference
			a["match_qual"]=a.cmp("frc",self.aliimg,{"sweight":0})		# compute similarity to unmasked reference
			self.particles_ali.append((a["match_qual"],a))
		self.particles_ali.sort()
		self.w2dptcl.set_data([i for j,i in self.particles_ali])

		self.wpbprogress.reset()
		self.wpbprogress.setEnabled(False)

	def aliRecalcRefPress(self,x=False):
		if len(self.particles)==0 : return
		
		self.alignToRef()			# realign particles to current masked reference
		
		# Compute the new average
		self.alisig=EMData(self.particles[0]["nx"],self.particles[0]["ny"],1)
		avgr=Averagers.get("mean",{"sigma":self.alisig})
		for q,i in self.particles_ali[:len(self.particles)/2] : avgr.add_image(i)
		self.aliimg=avgr.finish()
		
		self.setAliRef(self.aliimg)

	def aliRRecalcRefPress(self,x=False):
		if len(self.particles)==0 : return
		
		self.rAlignToRef()			# realign particles to current masked reference
		
		# Compute the new average
		self.alisig=EMData(self.particles[0]["nx"],self.particles[0]["ny"],1)
		avgr=Averagers.get("mean",{"sigma":self.alisig})
		for q,i in self.particles_ali[:len(self.particles)/2] : avgr.add_image(i)
		self.aliimg=avgr.finish()
		
		self.setAliRef(self.aliimg)
	
	def roiDrawMode(self,x=False):
		pass
	
	def roiAutoPress(self,x=False):
		pass
	
	def roiResetPress(self,x=False):
		self.roimask=self.aliimg.copy()
		self.roimask.add(-self.roimask["minimum"]+.01)		# The minimum value in the image is >0, so when we draw with 0 'color' we can extract the mask info
		self.w2droimaskdraw.set_data(self.roimask)
		self.roiGoPress()
	
	def roiGoPress(self,x=False):
		self.roimasked=self.aliimg.copy()
		
		blur=self.wsbroimaskblur.value()
		
		mask=self.roimask.process("threshold.binary",{"value":0.001})			# binarize the drawn mask
		mask.process_inplace("math.linear",{"scale":-1.0,"shift":1.0})		# invert the mask (user selects the region to include, not exclude)
		mask.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.5/(blur+.01)})
	
		self.roimasked.mult(mask)
		self.w2droimask.set_data(self.roimasked)
	
	def doCompute(self,x=False):
		pass
	
	

if __name__ == "__main__":
	main()
