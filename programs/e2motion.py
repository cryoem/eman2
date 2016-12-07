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
import Queue
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

	parser.add_argument("--threads", default=0,type=int,help="Number of alignment threads to run in parallel on a single computer. This is the only parallelism supported by e2spt_align at present.")
	parser.add_argument("--path",type=str,default=None,help="Path for the refinement, default=auto")
	parser.add_argument("--iter",type=int,help="Iteration number within path. Default = start a new iteration",default=0)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	global options
	(options, args) = parser.parse_args()
	
	if options.path==None:
		options.path=numbered_path("m2d",True)
#		os.makedirs(options.path)

	if options.threads<1 : options.threads=num_cpus()

	if options.path==None:
		fls=[int(i[-2:]) for i in os.listdir(".") if i[:4]=="m2d_" and len(i)==6 and str.isdigit(i[-2:])]
		if len(fls)==0 : fls=[1]
		options.path = "m2d_{:02d}".format(max(fls))
		if options.verbose : print "Using --path ",options.path
		
	if not os.path.exists(options.path) :
		os.mkdir(options.path)
		if itr==0 : itr=1
		
	parms=js_open_dict("{}/0_a2d_parms.json".format(options.path))
	
	if not parms.has_key(options.iter) :
		try: options.iter=max([int(i) for i in parms.keys()])
		except: options.iter=0
		print "Iteration: ",options.iter

	pid=E2init(argv)
	
	app = EMApp()
	motion=EMMotion(app,options.path,options.iter,options.threads)
	motion.show()
	app.execute()
	
	E2end(pid)

class EMMotion(QtGui.QMainWindow):
	"""This is the main window for the EMMotion application"""
	
	def __init__(self,application,path=None,piter=None,threads=4):
		"""application is an QApplication instance. ptclstack is the path to the file containing the particles to analyze. path is the path for ouput files""" 
		QtGui.QWidget.__init__(self)

		self.aliimg=None		# This is the unmasked alignment reference image
		self.alisig=None		# This is the standard deviation of the alignment reference
		self.alimask=None		# This is the drawn-upon version of aliimg used to generate the mask
		self.alimasked=None		# This is the masked alignment reference used for the actual alignments
		self.roidrawmask=None	# Region of interest mask drawing image
		self.roimask=None		# Region of interest actual mask
		self.roimasked=None		# Region of interest display widget

		self.particles=None

		self.app=weakref.ref(application)
		self.path=path
		self.iter=piter

		self.setWindowTitle("Main Window (e2motion.py)")

#		self.setWindowTitle("e2motion.py")
		
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
		
		######
		# Folder parameters
		self.vgb0=QtGui.QGroupBox("Particle Data")
		self.gbl.addWidget(self.vgb0,0,0,1,4)
		
		self.gbl2=QtGui.QGridLayout(self.vgb0)
		self.wlpath=QtGui.QLabel("Path: {}".format(self.path))
		self.gbl2.addWidget(self.wlpath,0,0)
		self.gbl2.setColumnStretch(0,1)
		
		self.wvbiter=ValBox(label="Iter:",value=self.iter)
		self.wvbiter.setIntonly(True)
		self.gbl2.addWidget(self.wvbiter,0,2)
		self.gbl2.setColumnStretch(2,0)
		
		self.wvsnum=ValSlider(rng=(-.2,0),label="Nptcl:",value=250)
		self.wvsnum.setIntonly(True)
		self.gbl2.addWidget(self.wvsnum,0,3)
		self.gbl2.setColumnStretch(3,4)
		
		self.wlnptcl=QtGui.QLabel(" ")
		self.gbl2.addWidget(self.wlnptcl,0,5)
		self.gbl2.setColumnStretch(5,2)
		
		self.wbdoavg=QtGui.QPushButton("Make Avg")
		self.gbl2.addWidget(self.wbdoavg,0,8)
		
		###### Alignment Mask
		# widget for editing the alignment mask
		self.wlalimaskdraw=QtGui.QLabel("<big><pre>Edit</pre></big>")
		self.wlalimaskdraw.setAlignment(Qt.AlignHCenter)
		self.gbl.addWidget(self.wlalimaskdraw,2,1)
		
		self.wlalimaskdraw2=QtGui.QLabel("<big><pre>A\nl\ni\ng\nn</pre></big>")
		self.gbl.addWidget(self.wlalimaskdraw2,3,0)
		
		self.w2dalimaskdraw=EMImage2DWidget()
		self.gbl.addWidget(self.w2dalimaskdraw,3,1)
		
		# Buttons for controlling mask
		self.hbl1=QtGui.QHBoxLayout()
		self.gbl.addLayout(self.hbl1,4,1)
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
		self.gbl.addLayout(self.vbl1,3,2)
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

		self.wlalimaskrot=QtGui.QLabel("Rot")
		self.vbl1.addWidget(self.wlalimaskrot)
		
		self.wsbalimaskrot=QtGui.QSpinBox()
		self.wsbalimaskrot.setRange(0,360)
		self.wsbalimaskrot.setValue(0)
		self.vbl1.addWidget(self.wsbalimaskrot)
		
		self.vbl1.addSpacing(16)
		
		self.wbaligo=QtGui.QPushButton(QtCore.QChar(0x2192))
		self.vbl1.addWidget(self.wbaligo)

		self.vbl1.addStretch(5)
		
		# widget for displaying the masked alignment reference
		self.wlalimask=QtGui.QLabel("<big><pre>Reference</pre></big>")
		self.wlalimask.setAlignment(Qt.AlignHCenter)
		self.gbl.addWidget(self.wlalimask,2,3)
		
		self.w2dalimask=EMImage2DWidget()
		self.gbl.addWidget(self.w2dalimask,3,3)
		
		self.hbl1a=QtGui.QHBoxLayout()
		self.gbl.addLayout(self.hbl1a,4,3)
		self.hbl1a.addStretch(5)
		
		self.wbrecalcref=QtGui.QPushButton("Realign")
		self.hbl1a.addWidget(self.wbrecalcref)
		
		self.wbrrecalcref=QtGui.QPushButton("Rerefine")
		self.hbl1a.addWidget(self.wbrrecalcref)
		
		self.hbl1a.addStretch(5)


		###### ROI Mask
		# widget for editing the ROI mask
		self.wlroimaskdraw=QtGui.QLabel("<big><pre>R\nO\nI</pre></big>")
		self.gbl.addWidget(self.wlroimaskdraw,6,0)
		
		self.w2droimaskdraw=EMImage2DWidget()
		self.gbl.addWidget(self.w2droimaskdraw,6,1)

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
		self.gbl.addLayout(self.vbl2,6,2)
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
		self.gbl.addWidget(self.w2droimask,6,3)
		
		self.vbl2.addStretch(5)

		self.wlarrow1=QtGui.QLabel(QtCore.QChar(0x2192))
		self.gbl.addWidget(self.wlarrow1,4,4)

		###### Results
		# Widget showing lists of different result sets
		self.vbl3=QtGui.QVBoxLayout()
		self.gbl.addLayout(self.vbl3,3,6,5,1)
		
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

		self.wvbnbasis=ValBox(None,(0,64),"# PCA Vec",8)
		self.wvbnbasis.setIntonly(True)
		self.vbl3a.addWidget(self.wvbnbasis)

		#self.wvbptclpct=ValBox(None,(0,100),"% Ptcl Incl",60)
		#self.wvbptclpct.setIntonly(True)
		#self.vbl3a.addWidget(self.wvbptclpct)

		
		## fill in a default value for number of threads
		#try :
			#cores=num_cpus()
		#except:
			#cores=2
		#if cores==1 : cores=2
		cores=threads+1		# one thread for GUI
		
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
		self.gbl.addWidget(self.wlarrow2,4,7)


		###### Output widgets
		# Class-averages
		self.wlclasses=QtGui.QLabel("<big><pre>Classes</pre></big>")
		self.wlclasses.setAlignment(Qt.AlignHCenter)
		self.gbl.addWidget(self.wlclasses,2,9)
		
		self.w2dclasses=EMImage2DWidget()
		self.gbl.addWidget(self.w2dclasses,3,9)

		self.wbshowptcl=QtGui.QPushButton(QtCore.QChar(0x2193))
		self.gbl.addWidget(self.wbshowptcl,4,9)

		self.w2dptcl=EMImage2DWidget()
		self.gbl.addWidget(self.w2dptcl,6,9)
		
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
		QtCore.QObject.connect(self.wbshowptcl,QtCore.SIGNAL("clicked(bool)"),self.showParticles)
		QtCore.QObject.connect(self.wvbiter,QtCore.SIGNAL("valueChanged"),self.newIter)
		QtCore.QObject.connect(self.wvsnum,QtCore.SIGNAL("valueChanged"),self.newThresh)
		QtCore.QObject.connect(self.wbdoavg,QtCore.SIGNAL("clicked(bool)"),self.avgPress)

		QtCore.QObject.connect(self.mfileopen,QtCore.SIGNAL("triggered(bool)")  ,self.menuFileOpen  )


		# set up draw mode
		insp=self.w2dalimaskdraw.get_inspector()
		insp.hide()
		insp.mmtab.setCurrentIndex(5)
		insp.dtpenv.setText("0.0")
		
		insp=self.w2droimaskdraw.get_inspector()
		insp.hide()
		insp.mmtab.setCurrentIndex(5)
		insp.dtpenv.setText("0.0")
		

		self.path=path
		QtCore.QTimer.singleShot(500,self.afterStart)

	def afterStart(self):
		"""This gets called once the event loop is running"""
#		print "App running, initializing"
		self.initPath(self.path,self.iter)

	def initPath(self,path,itr):
		
		parms=js_open_dict("{}/0_a2d_parms.json".format(options.path))

		self.wvbiter.setValue(itr)
		self.newIter()
		
		return

	def newIter(self,x=0):
		itr=int(self.wvbiter.getValue())
		
		try: 
			dct=js_open_dict("{}/particle_parms_{:02d}.json".format(self.path,itr))
			self.particles=[(j["score"],j["xform.align2d"],eval(i)[0],int(eval(i)[1])) for i,j in dct.items()]
			self.particles.sort()
			if len(self.particles)==0 : raise Exception
		except:
			self.particles=None
			self.wlnptcl.setText("No Data")
			print "Warning: no particle alignment data found for iter=",itr
			return
		
		m=0
		s=0
		for i in self.particles:
			m+=i[0]
			s+=i[0]**2
		m/=len(self.particles)
		s=sqrt(s/len(self.particles)-m**2)
		self.ptclmean=m
		self.ptclsigma=s
		
		self.wvsnum.setRange(1,len(self.particles)/10)
				
		self.newThresh()

	def newThresh(self,x=0):
		if self.particles==None or len(self.particles)<3:
			self.wlnptcl.setText("No Data")
			self.setAliRef(None)
			return
		
		n2use=self.wvsnum.getValue()
		if n2use>len(self.particles): 
			n2use=len(self.particles)
			self.wvsnum.setValue(n2use)
			self.wvsnum.setRange(1,n2use)
		self.wlnptcl.setText("/{:1d}    {:1.3f} sigma, {:1.1f} %".format(len(self.particles),(self.particles[n2use-1][0]-self.ptclmean)/self.ptclsigma,100.0*float(n2use)/len(self.particles)))

	def avgPress(self,x=0):
		if self.particles==None or len(self.particles)<3:
			self.wlnptcl.setText("No Data")
			self.setAliRef(None)
			return
		
		n2use=self.wvsnum.getValue()
		avg=Averagers.get("mean")
		for i in self.particles[:n2use]:
			img=EMData(i[2],i[3]).process("xform",{"transform":i[1]})
			avg.add_image(img)
		
		self.aliimg=avg.finish()
		
		self.setAliRef(self.aliimg)
		
#		self.wlnptcl.setText("{:1.3f} *sigma, {:1.1f} %".format((self.particles[n2use][0]-self.ptclmean)/self.ptclsigma,100.0*float(n2use)/len(self.particles)))
		self.wlnptcl.setText("/{:1d}    {:1.3f} sigma, {:1.1f} %".format(len(self.particles),(self.particles[n2use][0]-self.ptclmean)/self.ptclsigma,100.0*float(n2use)/len(self.particles)))
		

	def menuFileOpen(self,x):
		if self.particles!=None:
			QtGui.QMessageBox.warning(None,"Error","%s already contains a stack of particles. A new folder is required to start with a new stack of particles. Rerun without --path option."%self.path)
			return

		self.dialog = embrowser.EMBrowserWidget(withmodal=True,multiselect=False)
		QtCore.QObject.connect(self.dialog,QtCore.SIGNAL("ok"),self.gotPath)
		QtCore.QObject.connect(self.dialog,QtCore.SIGNAL("cancel"),self.gotPath)
		self.dialog.show()
	
	def gotPath(self):
		#ptcl=self.dialog.getResult()
		#self.dialog=None
		#if ptcl == None : return		# user pressed cancel, but we still need to clean up
		#ptcl=ptcl[0]
		
		#n=EMUtil.get_image_count(ptcl)
		#sz=EMData(ptcl,0,True)["nx"]
		
		## We don't want the user to overwhelm the system
		#if sz*sz*4*n > 5.0e8 :
			#nmax=5.0e8/(sz*sz*4)
			#r=QtGui.QMessageBox.question(None,"Are you sure ?","WARNING: This full particle set will require %d+ gb memory to process. Select Yes to use only the first %d particles, No to use the entire stack, or Cancel to abort. You may also consider using e2proc2d.py --meanshrink= to produce a downsampled stack for processing."%(int(sz*sz*12*n/1.0e9+.5),nmax),QtGui.QMessageBox.Yes|QtGui.QMessageBox.No|QtGui.QMessageBox.Cancel)
			#if r==QtGui.QMessageBox.Cancel : return
			#if r==QtGui.QMessageBox.Yes : n=nmax
		
		#task="e2proclst.py %s --create %s/particles.lst --range=0,%d,1"%(ptcl,self.path,n)
		#print task
		#launch_childprocess(task)
		
		#self.initPath(self.path)
		print "does nothing"
	
	def threadAlign(self,curlist,ref,refnomask,outstack):
		"""This method is designed to run in a thread and perform alignments of a stack of inputs to a single reference"""
	
		for p in curlist:
#			p2=p.process("filter.lowpass.gauss",{"cutoff_abs":0.1})
			p2=EMData(p[2],p[3])
			a=p2.align("rotate_translate_tree",ref)
#			a=p2.copy()
#			a["match_qual"]=a.cmp("frc",ref,{"sweight":0})		# compute similarity to unmasked reference
			a["match_qual"]=a.cmp("ccc",refnomask)		# compute similarity to unmasked reference
			outstack.put((a["match_qual"],a["xform.align2d"],p[2],p[3]))

	def threadRAlign(self,curlist,ref,refnomask,outstack):
		"""This method is designed to run in a thread and perform alignments of a stack of inputs to a single reference."""
	
		for i,p in enumerate(curlist):
#			p2=p.process("filter.lowpass.gauss",{"cutoff_abs":0.1})
			p2=EMData(p[2],p[3])
			a=p2.align("refine",ref,{"verbose":0,"xform.align2d":p[1]},"ccc",{})		# doesn't really seem to make an improvement
#			a["match_qual"]=a.cmp("frc",ref,{"sweight":0})		# compute similarity to unmasked reference
			a["match_qual"]=a.cmp("ccc",refnomask)		# compute similarity to unmasked reference
			outstack.put((a["match_qual"],a["xform.align2d"],p[2],p[3]))

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
			c=a.align("rotate_translate_tree",b)
			c.add(b)
			c.process_inplace("xform.centerofmass",{"threshold":c["mean"]+c["sigma"]})
#			c.process_inplace("xform.center")
			tree2.append(c)

	def threadFilt(self,ptclstack,outstack):
		"""Allows particles to be shrunk/filtered in parallel"""
		
		for i in ptclstack:
			j=i.process("math.meanshrink",{"n":2})
			j.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.1})
			outstack.append(j)

	def waitThreads(self,thrs,lst,maxl):
		
		# Wait for threads to finish
		while 1:
			time.sleep(0.2)
			self.wpbprogress.setValue(int(lst.qsize()*100/maxl))
			QtGui.qApp.processEvents()

			# If any threads are alive, it breaks out of the inner loop, if none are alive, the else block breaks out of the outer loop
			for t in thrs:
				if t.isAlive() : break
			else:
				break
				
		self.wpbprogress.reset()


	def bootstrap(self):
		"""This will create an initial alignment reference for the particles when no reference is available"""

		self.wpbprogress.setEnabled(True)
		self.wpbprogress.reset()
		QtGui.qApp.processEvents()
		nthr=int(self.wvbcores.getValue())		# number of threads to use for faster alignments

		print "bs1"
		tree=[]
		# launch nthr threads to do the alignments
		thrs=[]
		for i in range(nthr):
			thrs.append(threading.Thread(target=self.threadFilt, args=(self.particles,tree)))
			thrs[-1].start()
		
		
		self.waitThreads(thrs,tree,len(self.particles))
		
		print "bs2"
		lock=threading.Lock()
		
		while len(tree)>1:
			maxt=len(tree)
			tree2=[]
			
			# launch nthr threads to do the alignments
			thrs=[]
			for i in range(nthr):
				thrs.append(threading.Thread(target=self.threadTreeAlign, args=(tree,tree2,lock)))
				thrs[-1].start()
				
			self.waitThreads(thrs,tree2,nthr)

			tree=tree2
		
#		tree[0]=tree[0].get_clip(Region(-tree[0]["nx"]/2,-tree[0]["ny"]/2,tree[0]["nx"]*2,tree[0]["ny"]*2))
		tree[0].process_inplace("xform.scale",{"clip":tree[0]["nx"]*2,"scale":2.0})
		tree[0].process_inplace("normalize.edgemean")
#		tree[0].process_inplace("xform.centerofmass",{"threshold":0.5})
#		tree[0].process_inplace("normalize.edgemean")

		# We look for the most 'featurefull' direction to put on Y
		f=tree[0].do_fft()
		f.ri2ap()
		pl=f.calc_az_dist(90,-90.0,2.0,3,16)
#		for i in range(90): print "%d\t%f"%(i,pl[i])
		ml=2.0*pl.index(max(pl))
		tree[0].rotate(ml,0,0)				# put the initial reference in the preferred orientation
		#tree[0].write_image("zzz.hdf",0)
		
		print "bs3"
		# One final alignment pass
		self.particles_ali=[]
		self.alisig=EMData(self.particles[0]["nx"],self.particles[0]["ny"],1)
		avgr=Averagers.get("mean",{"sigma":self.alisig})
		thrs=[]
		
		# launch nthr threads to do the alignments
		for i in range(nthr):
			thrs.append(threading.Thread(target=self.threadAlign, args=(self.particles[i::nthr],tree[0],tree[0],self.particles_ali)))
			thrs[-1].start()

		self.waitThreads(thrs,self.particles_ali,len(self.particles))
		
		print "bs4"
		self.particles_ali.sort()
		
		#out=file("x.txt","w")
		#for i in self.particles_ali: out.write("%f\n"%i[1]["match_qual"])
		#out.close()
		
#		ptclfrac=self.wvbptclpct.getValue()/100.0
		for i,a in self.particles_ali[:int(len(self.particles)*ptclfrac)]: avgr.add_image(a)
		
		self.aliimg=avgr.finish()

		self.w2dptcl.set_data([i for j,i in self.particles_ali])
		

		#self.aliimg.write_image("zzz.hdf",1)
		#self.alisig.write_image("zzz.hdf",2)
		
		print "bs5"
		self.setAliRef(self.aliimg)
		
		self.wpbprogress.reset()
		self.wpbprogress.setEnabled(False)
		
	
	def setAliRef(self,img):
		"""Sets a new alignment reference, creating a new empty mask if none is present"""
		
		if img==None : return
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
		if self.roidrawmask!=None:
			rmask=self.roidrawmask
			rmask.process_inplace("threshold.binary",{"value":0.001})
			
		self.roidrawmask=img.copy()
		self.roidrawmask.add(-self.roidrawmask["minimum"]+.01)		# The minimum value in the image is >0, so when we draw with 0 'color' we can extract the mask info
		if rmask!=None: self.roidrawmask.mult(rmask)
		
		self.w2droimaskdraw.set_data(self.roidrawmask)
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
		rot=self.wsbalimaskrot.value()
		
		mask=self.alimask.process("threshold.binary",{"value":0.001})			# binarize the drawn mask
		mask.process_inplace("math.linear",{"scale":-1.0,"shift":1.0+base/100.0})		# invert the mask (user selects the region to include, not exclude)
		mask.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.5/(blur+.01)})
		
		self.alimasked.mult(mask)
		if rot!=0 :
			self.alimasked.process_inplace("xform",{"transform":Transform({"type":"2d","alpha":rot})})

		self.w2dalimask.set_data(self.alimasked)

	def rAlignToRef(self):
		"""realigns particles to reference, but does not compute a new reference"""
		self.wpbprogress.setEnabled(True)
		self.wpbprogress.reset()
		
		nthr=int(self.wvbcores.getValue())		# number of threads to use for faster alignments
		
		jsd=Queue.Queue(0)
		self.particles_ali=[]
		thrs=[]
		# launch nthr threads to do the alignments
		for i in range(nthr):
			thrs.append(threading.Thread(target=self.threadRAlign, args=(self.particles[i::nthr],self.alimasked,self.aliimg,jsd)))
			thrs[-1].start()
			
		self.waitThreads(thrs,jsd,len(self.particles))

		# for refinement we use the current iteration
		itr=self.wvbiter.getValue()
		
		dct=js_open_dict("{}/particle_parms_{:02d}.json".format(self.path,itr))
		while not jsd.empty():
			i=jsd.get()
			dct[(i[2],i[3])]={"xform.align2d":i[1],"score":i[0]}
		
		self.wpbprogress.reset()
		self.wpbprogress.setEnabled(False)
		
			
	def alignToRef(self):
		"""realigns particles to reference, but does not compute a new reference"""
		self.wpbprogress.setEnabled(True)
		self.wpbprogress.reset()
		
		nthr=int(self.wvbcores.getValue())		# number of threads to use for faster alignments

		jsd=Queue.Queue(0)
		n2use=self.wvsnum.getValue()
		thrs=[]
		# launch nthr threads to do the alignments
		for i in range(nthr):
			thrs.append(threading.Thread(target=self.threadAlign, args=(self.particles[i:n2use:nthr],self.alimasked,self.aliimg,jsd)))
			thrs[-1].start()
			
		self.waitThreads(thrs,jsd,len(self.particles))

		# we find a new iteration to use for the new alignment
		itr=1
		while os.path.exists("{}/particle_parms_{:02d}.json".format(self.path,itr)): itr+=1
		
		dct=js_open_dict("{}/particle_parms_{:02d}.json".format(self.path,itr))
		while not jsd.empty():
			i=jsd.get()
			dct[(i[2],i[3])]={"xform.align2d":i[1],"score":i[0]}
			
		self.wvbiter.setValue(itr)
		
		
		self.wpbprogress.reset()
		self.wpbprogress.setEnabled(False)

	def aliRecalcRefPress(self,x=False):
		if len(self.particles)==0 : return
		
		self.alignToRef()			# realign particles to current masked reference
		
		## Compute the new average
		#self.alisig=EMData(self.particles[0]["nx"],self.particles[0]["ny"],1)
		#avgr=Averagers.get("mean",{"sigma":self.alisig})
		#ptclfrac=self.wvbptclpct.getValue()/100.0
		#for q,i in self.particles_ali[:int(len(self.particles)*ptclfrac)] : avgr.add_image(i)
		#self.aliimg=avgr.finish()
		
		#self.setAliRef(self.aliimg)

	def aliRRecalcRefPress(self,x=False):
		if len(self.particles)==0 : return
		
		self.rAlignToRef()			# realign particles to current masked reference
		
		## Compute the new average
		#self.alisig=EMData(self.particles[0]["nx"],self.particles[0]["ny"],1)
		#avgr=Averagers.get("mean",{"sigma":self.alisig})
		#ptclfrac=self.wvbptclpct.getValue()/100.0
		#for q,i in self.particles_ali[:int(len(self.particles)*ptclfrac)] : avgr.add_image(i)
		#self.aliimg=avgr.finish()
		
		#self.setAliRef(self.aliimg)
	
	def roiDrawMode(self,x=False):
		pass
	
	def roiAutoPress(self,x=False):
		pass
	
	def roiResetPress(self,x=False):
		self.roidrawmask=self.aliimg.copy()
		self.roidrawmask.add(-self.roidrawmask["minimum"]+.01)		# The minimum value in the image is >0, so when we draw with 0 'color' we can extract the mask info
		self.w2droimaskdraw.set_data(self.roidrawmask)
		self.roiGoPress()
	
	def roiGoPress(self,x=False):
		self.roimasked=self.aliimg.copy()
		
		blur=self.wsbroimaskblur.value()
		
		mask=self.roidrawmask.process("threshold.binary",{"value":0.001})			# binarize the drawn mask
		mask.process_inplace("math.linear",{"scale":-1.0,"shift":1.0})		# invert the mask (user selects the region to include, not exclude)
		mask.process_inplace("filter.lowpass.gauss",{"cutoff_abs":0.5/(blur+.01)})
#		mask.process_inplace("threshold.belowtozero",{"minvalue":0.05})				# this limits the range of the mask to improve PCA performance
	
		self.roimask=mask
		self.roimasked.mult(mask)
		self.w2droimask.set_data(self.roimasked)
	
	def showParticles(self):
		pass
	
	def doCompute(self,x=False):
		mode=self.wcbprocmode.currentIndex()
		
		if mode==0 : self.doComputePCA()
		elif mode==1 : self.doComputeAvD()
		else : QtGui.QMessageBox.warning("Unknown mode %d"%mode)
		
		return
		
	def doComputePCA(self):
		"""Particle classification by MSA"""
		
		n2use=self.wvsnum.getValue()
		toclass=self.particles[:n2use]

		nclasses=self.wvbclasses.getValue()
		nbasis=self.wvbnbasis.getValue()
		
		
		# Find a class number in the current iteration
		clnums=[i.split("_")[-1][:2] for i in os.listdir(self.path) if "classes_{:02d}".format(self.iter) in i]
		for i in xrange(len(clnums)):
			try: clnums[i]=int(clnums[i])
			except: clnums[i]=0
		clnums.append(0)
		clnum=max(clnums)+1		# number of the next free class averagee
		
		self.wpbprogress.setValue(10.0)
		
		# compute PCA
		ptcl4an=[EMData(name,i).process("xform",{"transform":xform}).process("normalize.edgemean") for score,xform,name,i in toclass]
		pca=Analyzers.get("pca_large",{"mask":self.roimask,"nvec":nbasis})
		for p in ptcl4an:
			p2=p.copy()
			p2.mult(self.roimask)
			p2.process_inplace("normalize.unitlen")
			pca.insert_image(p2)
			
		basis=pca.analyze()
		for b in basis: b.write_image("{}/basis_{:02d}_{:02d}.hdf".format(self.path,self.iter,clnum),-1)
		self.roimask.write_image("{}/maskroi_{:02d}_{:02d}.hdf".format(self.path,self.iter,clnum),0)

		self.wpbprogress.setValue(33.0)
		
		# project each particle into the basis subspace
#		proj=[[p.cmp("dot",b,{"normalize":0,"negative":0} for b in basis] for p in ptcl4an]
		projs=[]					# 1 EMData per particle, each is a basis subspace projection (1-D image)
		for p in ptcl4an:
			proj=EMData(nbasis-1,1,1)
			projs.append(proj)
			for i,b in enumerate(basis[1:]):
				proj.set_value_at(i,0,p.cmp("dot",b,{"normalize":0,"negative":0}))
				
		# k-means classification
		an=Analyzers.get("kmeans",{"ncls":nclasses,"minchange":n2use//100+1,"slowseed":0,"mininclass":min((n2use//(nclasses*4)+1),10)})
		
		an.insert_images_list(projs)
		centers=an.analyze()
		self.wpbprogress.setValue(66.0)

		classes=[None for i in range(nclasses)]
		classlst=[[] for i in range(nclasses)]
		for n,i in enumerate(ptcl4an):
			try: classes[projs[n]["class_id"]].add(i)
			except: classes[projs[n]["class_id"]]=i.copy()
			classlst[projs[n]["class_id"]].append(n)
		
		
		# Make class-averages
		fsp="{}/classes_{:02d}_{:02d}.hdf".format(self.path,self.iter,clnum)
		print fsp
		for i,c in enumerate(classes): 
			self.wpbprogress.setValue(66+33*i/len(classes))
			c.process_inplace("normalize.edgemean")
			c["class_ptcl_idxs"]=classlst[i]
#			c["exc_class_ptcl_idxs"]=[]
			c["class_ptcl_src"]=toclass[0][2]
			c.write_image(fsp,-1)
		
		
			
		self.classes=classes
		self.w2dclasses.set_data(self.classes)
		
		self.wpbprogress.reset()

		
	def doComputeAvD(self):
		"""Particle classification based on average density within the mask"""
		
		n2use=self.wvsnum.getValue()
		toclass=self.particles[:n2use]
		
		tosort=[]
		for score,xform,name,i in toclass:
			img=EMData(name,i).process("xform",{"transform":xform}).process("normalize")
			imgm=img.copy()
			imgm.mult(self.roimask)
			tosort.append((imgm["mean"],img))
		
		nclasses=self.wvbclasses.getValue()

#		c.process_inplace("normalize.toimage",{"to":self.alimasked,"ignore_zero":1})
		
		# sort by density
		tosort.sort()
		
		# Make the class-averages
		self.classes=[]
		clssz=len(tosort)/nclasses		# particles per class
		
		for cl in xrange(nclasses):
			avgr=Averagers.get("mean")
			for i in tosort[cl*clssz:(cl+1)*clssz]: avgr.add_image(i[1])
			self.classes.append(avgr.finish())

		self.w2dclasses.set_data(self.classes)
	

if __name__ == "__main__":
	main()
