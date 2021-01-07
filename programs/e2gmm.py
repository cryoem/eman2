#!/usr/bin/env python
#
# Author: Steven Ludtke  9/28/2020 
# Copyright (c) 2020- Baylor College of Medicine
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
# This is a basic GUI for Muyuan's GMM processing

from past.utils import old_div
from future import standard_library
standard_library.install_aliases()
from builtins import range
import sys
import os
import weakref
import threading
import time
from sys import argv
from EMAN2 import *
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt
from eman2_gui.emapplication import get_application, EMApp
from eman2_gui.emplot2d import EMPlot2DWidget
from eman2_gui.emimage2d import EMImage2DWidget
from eman2_gui.emscene3d import EMScene3D
from eman2_gui.emshapeitem3d import EMScatterPlot3D 
from eman2_gui.valslider import *
import queue
from eman2_gui import embrowser
import sklearn.decomposition as skdc


os.environ["CUDA_VISIBLE_DEVICES"]='0' 
os.environ["TF_FORCE_GPU_ALLOW_GROWTH"]='true' 
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' #### reduce log output

import traceback
import tensorflow as tf
import tensorflow.keras.models
import numpy as np

def butval(widg):
	if widg.isChecked(): return "1"
	return "0"

# shortcut, if we had a central GUI file like EMAN2.py, it could go there...
def showerror(msg,parent=None):
	QtWidgets.QMessageBox.warning(parent,"ERROR",msg)

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options]

	WARNING: This program still under development.
	
	This program is designed to work interactively with Gaussian mixture models, automating as much 
	as possible of the tasks associated with studying dynamics using GMMs. It requires either a 
	refine_XX folder or a set of particles with known orientations as a starting point."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

#	parser.add_argument("--path",type=str,default=None,help="Path for the gmm_XX folder to work in, by default a new folder will be created")
	parser.add_argument("--threads", default=-1,type=int,help="Number of alignment threads to run in parallel on a single computer. This is the only parallelism supported by e2spt_align at present.")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	global options
	(options, args) = parser.parse_args()
	
	#if options.path==None:
		#options.path=numbered_path("gmm",True)
#		os.makedirs(options.path)

	if options.threads<1 : options.threads=num_cpus()
		
	#if not os.path.exists(options.path) :
		#os.mkdir(options.path)
		
#	parms=js_open_dict("{}/0_a2d_parms.json".format(options.path))
	

	pid=E2init(argv)
	
	app = EMApp()
	emgmm=EMGMM(app,options)
	emgmm.show()
	emgmm.raise_()
	app.execute()
	
	E2end(pid)

class EMGMM(QtWidgets.QMainWindow):
	"""This is the main window for the e2gmm application"""
	
	def __init__(self,application,opt):
		self.options=opt
		"""application is an QApplication instance. path is the path for ouput files""" 
		QtWidgets.QWidget.__init__(self)

		self.particles=None

		self.app=weakref.ref(application)

		self.setWindowTitle("Main Window (e2gmm.py)")

		# Menu Bar
		#self.mfile=self.menuBar().addMenu("File")
		#self.mfileopen=self.mfile.addAction("Select Particles")
		#self.mfileopencls=self.mfile.addAction("Particles from Classes")
		#self.mfileopencls.setEnabled(False)
		#self.mfilequit=self.mfile.addAction("Quit")

		self.setCentralWidget(QtWidgets.QWidget())
		self.gbl = QtWidgets.QGridLayout(self.centralWidget())
		cen=self.centralWidget()
		
		self.gbl.setColumnStretch(0,0)
		self.gbl.setColumnStretch(2,4)
		self.gbl.setColumnStretch(4,5)
		
		# Plot
		self.wplot2d = EMPlot2DWidget()
		self.gbl.addWidget(self.wplot2d,0,2,3,1)
		self.wplot2d.set_mouse_emit(True)
		
		# 3-D View
		self.wview3d = EMScene3D()
		self.gbl.addWidget(self.wview3d,0,4,4,1)
		
		# Left pane for GMM folder and parameters
		self.gbll = QtWidgets.QGridLayout()
		self.gbl.addLayout(self.gbll,0,0,4,1)		# 4 rows tall, 1 column wide
		
		# gmm_XX folder selection
		self.gblfld = QtWidgets.QGridLayout()
		self.gbll.addLayout(self.gblfld,0,0)
		
		self.wlistgmm = QtWidgets.QListWidget()
		self.wlistgmm.setSizePolicy(QtWidgets.QSizePolicy.Preferred,QtWidgets.QSizePolicy.Expanding)
		self.update_gmms()
		self.gblfld.addWidget(self.wlistgmm,0,0,1,2)
		
		self.wbutnewgmm=QtWidgets.QPushButton("New GMM")
		self.gblfld.addWidget(self.wbutnewgmm,1,0)
		
		self.wbutrefine = QtWidgets.QPushButton("refine_XX")
		self.gblfld.addWidget(self.wbutrefine,1,1)
		
		# run selection
		self.gblrun = QtWidgets.QGridLayout()
		self.gbll.addLayout(self.gblrun,1,0)
		
		self.wlistrun = QtWidgets.QListWidget()
		self.wlistrun.setSizePolicy(QtWidgets.QSizePolicy.Preferred,QtWidgets.QSizePolicy.Expanding)
#		self.update_runs()
		self.gblrun.addWidget(self.wlistrun,0,0,1,2)
		
		self.wbutrerun=QtWidgets.QPushButton("Rerun")
		self.gblrun.addWidget(self.wbutrerun,1,0)
		
		self.wbutnewrun=QtWidgets.QPushButton("New Run")
		self.gblrun.addWidget(self.wbutnewrun,1,1)
		
		# The form with details about the selected gmm_XX folder
		self.gflparm = QtWidgets.QFormLayout()
		self.gbll.addLayout(self.gflparm,2,0)
		
		
		self.wedbox = QtWidgets.QLineEdit("256")
		self.wedbox.setToolTip("Box size of input particles in pixels")
		self.gflparm.addRow("Box Size:",self.wedbox)
		
		self.wedapix = QtWidgets.QLineEdit("1.0")
		self.wedapix.setToolTip("A/pix of input particles")
		self.gflparm.addRow("A/pix:",self.wedapix)

		self.wedmask = QtWidgets.QLineEdit("")
		self.wedmask.setToolTip("3-D volume mask")
		self.gflparm.addRow("Mask:",self.wedmask)
		
		self.wedngauss = QtWidgets.QLineEdit("256")
		self.wedngauss.setToolTip("Number of Gaussians")
		self.gflparm.addRow("N Gauss:",self.wedngauss)
		
		self.weddim = QtWidgets.QLineEdit("4")
		self.weddim.setToolTip("Number of dimensions in the latent space (middle network layer)")
		self.gflparm.addRow("Latent Dim:",self.weddim)
		
		self.wedtrainiter = QtWidgets.QLineEdit("10")
		self.wedtrainiter.setToolTip("Training iterations per stage")
		self.gflparm.addRow("Train iter:",self.wedtrainiter)
		
		self.wbutpos = QtWidgets.QPushButton("Position")
		self.wbutpos.setCheckable(True)
		self.wbutpos.setChecked(True)
		self.wbutpos.setToolTip("Include changes of position in the GMM (motion)")
		self.gflparm.addRow("Parameters",self.wbutpos)
		
		self.wbutamp = QtWidgets.QPushButton("Amplitude")
		self.wbutamp.setCheckable(True)
		self.wbutamp.setChecked(False)
		self.wbutamp.setToolTip("Include changes of amplitude in the GMM (ligand binding)")
		self.gflparm.addRow(" ",self.wbutamp)
		
		self.wbutsig = QtWidgets.QPushButton("Sigma")
		self.wbutsig.setCheckable(True)
		self.wbutsig.setChecked(False)
		self.wbutsig.setToolTip("Include changes of Gaussian Width in the GMM (rarely useful)")
		self.gflparm.addRow(" ",self.wbutsig)
		
		self.wlabruntime = QtWidgets.QLabel("-")
		self.gflparm.addRow("Run:",self.wlabruntime)
		
		# Widgets below plot
		self.gblpltctl = QtWidgets.QGridLayout()
		self.gbl.addLayout(self.gblpltctl,3,2)
		
		self.gblpltctl.addWidget(QtWidgets.QLabel("X Col:",self),0,0,Qt.AlignRight)
		self.wsbxcol=QtWidgets.QSpinBox(self)
		self.wsbxcol.setRange(0,10)
		self.gblpltctl.addWidget(self.wsbxcol,0,1,Qt.AlignLeft)
		self.wsbxcol.setValue(0)

		self.gblpltctl.addWidget(QtWidgets.QLabel("Y Col:",self),1,0,Qt.AlignRight)
		self.wsbycol=QtWidgets.QSpinBox(self)
		self.wsbycol.setRange(0,10)
		self.gblpltctl.addWidget(self.wsbycol,1,1,Qt.AlignLeft)
		self.wsbycol.setValue(1)
		
		self.wbutdrgrp=QtWidgets.QButtonGroup()
		
		self.wbutdrmid=QtWidgets.QPushButton("Net Mid")
		self.wbutdrmid.setCheckable(True)
		self.wbutdrmid.setChecked(True)
		self.gblpltctl.addWidget(self.wbutdrmid,0,2)
		self.wbutdrgrp.addButton(self.wbutdrmid,0)
		
		self.wbutdrpca=QtWidgets.QPushButton("PCA")
		self.wbutdrpca.setCheckable(True)
		self.gblpltctl.addWidget(self.wbutdrpca,0,3)
		self.wbutdrgrp.addButton(self.wbutdrpca,1)
		
		self.wbutdrpca=QtWidgets.QPushButton("ICA")
		self.wbutdrpca.setCheckable(True)
		self.gblpltctl.addWidget(self.wbutdrpca,0,4)
		self.wbutdrgrp.addButton(self.wbutdrpca,2)
		
		self.gblpltctl.addWidget(QtWidgets.QLabel("New Dim:",self),1,3,Qt.AlignRight)
		self.wsbnewdim=QtWidgets.QSpinBox(self)
		self.wsbnewdim.setRange(2,10)
		self.gblpltctl.addWidget(self.wsbnewdim,1,4,Qt.AlignLeft)
		
		self.wcbpntpln=QtWidgets.QComboBox()
		self.wcbpntpln.addItem("Plane")
		self.wcbpntpln.addItem("Point")
		self.wcbpntpln.addItem("Region")
		self.gblpltctl.addWidget(self.wcbpntpln,1,2)
		
		# Connections
		self.wlistgmm.currentRowChanged[int].connect(self.sel_gmm)
		self.wbutnewgmm.clicked[bool].connect(self.add_gmm)
		self.wbutrefine.clicked[bool].connect(self.setgmm_refine)
		self.wlistrun.currentRowChanged[int].connect(self.sel_run)
		self.wbutnewrun.clicked[bool].connect(self.new_run)
		self.wbutrerun.clicked[bool].connect(self.do_run)
		self.wbutdrgrp.buttonClicked[QtWidgets.QAbstractButton].connect(self.plot_mode_sel)
		self.wsbxcol.valueChanged[int].connect(self.wplot2d.setXAxisAll)
		self.wsbycol.valueChanged[int].connect(self.wplot2d.setYAxisAll)
		self.wplot2d.mouseDown[QtGui.QMouseEvent,tuple].connect(self.plot_mouse)
		self.wplot2d.mouseUp[QtGui.QMouseEvent,tuple].connect(self.plot_mouse)
		self.wplot2d.mouseDrag[QtGui.QMouseEvent,tuple].connect(self.plot_mouse)
		E2loadappwin("e2gmm","main",self)

		self.gaussplot=EMScatterPlot3D()
		self.wview3d.insertNewNode("Gauss Model",self.gaussplot)

		#QtCore.QTimer.singleShot(500,self.afterStart)

	def plot_mouse(self,event,loc):
		mmode=str(self.wcbpntpln.currentText())
		dim=self.currun.get("dim",4)
		latent=np.zeros(dim)
		if mmode=="Plane":
			if self.plotmode==0:
				latent[self.wsbxcol.value()]=loc[0]
				latent[self.wsbycol.value()]=loc[1]
			if self.plotmode==1:
				newdim=self.wsbnewdim.value()
				sel=np.zeros(newdim)
				sel[self.wsbxcol.value()]=loc[0]
				sel[self.wsbycol.value()]=loc[1]
				latent=self.decomp.inverse_transform(sel)
		elif mmode=="Point":
			try: sel=self.wplot2d.selected[0]
			except: return						# no new selected point
			latent=self.midresult[:,sel]
		elif mmode=="Region":
			return
		else: print("mode error")
		
		# run the current latent vector through the decoder then pull the result out in a useful 'shape'
		gauss=np.array(self.decoder(latent[None,...]))[0].transpose()
		box=int(self.wedbox.text())
		gauss[:3]*=box
		self.gaussplot.setData(gauss,gauss[4][0])
		self.wview3d.update()

	def plot_mode_sel(self,but):
		"""Plot mode selected"""
		
		self.plotmode=self.wbutdrgrp.id(but)
		newdim=self.wsbnewdim.value()
		
		if self.plotmode==0:
#			print("norm")
			self.data=self.midresult
		elif self.plotmode==1:
#			print("pca")
			self.decomp=skdc.PCA(n_components=newdim)
			self.data=self.decomp.fit_transform(self.midresult.transpose()).transpose()
		elif self.plotmode==2:
			pass
			
		self.wsbxcol.setRange(0,len(self.data)-1)
		self.wsbxcol.setValue(0)
		self.wsbycol.setRange(0,len(self.data)-1)
		self.wsbycol.setValue(1)
		
		self.wplot2d.set_data(self.data,"map")

	def new_run(self,clk=False):
		"""Create a new run and run() it"""
		name=str(QtWidgets.QInputDialog.getText(self,"Run Name","Enter a name for this run. Current parameters will be used.")[0])
		if not self.jsparm.has_key("run_"+name) : self.wlistrun.addItem(name)
		self.currunkey=name
		self.do_run()

	def do_run(self,clk=False):
		"""Run the current job with current parameters"""
		self.currun={}
		self.currun["ngauss"]=int(self.wedngauss.text())
		self.currun["dim"]=int(self.weddim.text())
		self.currun["mask"]=str(self.wedmask.text())
		self.currun["trainiter"]=int(self.wedtrainiter.text())
		self.currun["pas"]=butval(self.wbutpos)+butval(self.wbutamp)+butval(self.wbutsig)
		self.currun["time"]=local_datetime()
		self.jsparm["run_"+self.currunkey]=self.currun
		
		modelout=f"{self.gmm}/{self.currunkey}_model_gmm.txt"
		
		prog=QtWidgets.QProgressDialog("Running networks. Progress updates here are limited. See the Console for detailed output.","Abort",0,9)
		prog.show()
		# First step with very coarse model, gradually increasing size improves convergence
		run(f"e2gmm_refine.py --projs {self.gmm}/proj_in.hdf --npt {self.currun['ngauss']} --maxboxsz 24 --modelout {modelout} --niter 20 --mask {self.currun['mask']} --nmid {self.currun['dim']}")
		prog.setValue(1)
		self.app().processEvents()
		
		box=24
		n=2
		while box<self.jsparm["boxsize"]:
			box=good_size(box*2)
			# in the last iteration we do some additional things
			if box>=self.jsparm["boxsize"]:
				box=self.jsparm["boxsize"]
				s=f"--evalmodel {self.gmm}/{self.currunkey}_model_projs.hdf --evalsize {self.jsparm['boxsize']}"		# in the final iteration
			else: s=""
				
			# iterate until box size is the full size of the particles
			run(f"e2gmm_refine.py --projs {self.gmm}/proj_in.hdf --npt {self.currun['ngauss']} --maxboxsz {box} --model {self.gmm}/{self.currunkey}_model_gmm.txt --modelout {modelout} --niter {self.currun['trainiter']} --mask {self.currun['mask']} --nmid {self.currun['dim']} {s}")
			prog.setValue(n)
			self.app().processEvents()
			n+=1

		# make3d on gaussian output for comparison
		run(f"e2make3dpar.py --input {self.gmm}/{self.currunkey}_model_projs.hdf --output {self.gmm}/{self.currunkey}_model_recon.hdf --pad {good_size(self.jsparm['boxsize']*1.25)} --mode trilinear --keep 1 --threads {self.options.threads}")
		prog.setValue(8)
		self.app().processEvents()

		# heterogeneity analysis
		run(f"e2gmm_refine.py --model {modelout} --ptclsin {self.gmm}/particles.lst --heter --maxboxsz {self.jsparm['boxsize']} --gradout {self.gmm}/{self.currunkey}_grads.hdf --mask {self.currun['mask']} --midout {self.gmm}/{self.currunkey}_mid.txt --decoderout {self.gmm}/{self.currunkey}_decoder.h5 --pas {self.currun['pas']}")
		prog.setValue(9)
		self.app().processEvents()

	def update_gmms(self):
		"""Updates the display of gmm_XX folders"""
		#self.gmm=str(self.wlistgmm.currentItem().text())
		self.gmms=[i for i in sorted(os.listdir(".")) if i[:4]=="gmm_"]
		self.wlistgmm.clear()
		for i in self.gmms:
			self.wlistgmm.addItem(i)
	
	def sel_gmm(self,line):
		"""Called when the user selects a new GMM from the list. Should not be called directly by the user without updating wlistgmm"""

#		print("sel gmm",line)
		if line<0 : 
			self.wedbox.setText("")
			self.wedapix.setText("")
			self.wlistrun.clear()
			return
		self.gmm=str(self.wlistgmm.item(line).text())
		self.jsparm=js_open_dict(f"{self.gmm}/0_gmm_parms.json")
		# These are associated with the whole GMM
		self.wedbox.setText(f'{self.jsparm.getdefault("boxsize",128)}')
		self.wedapix.setText(f'{self.jsparm.getdefault("apix","")}')
		
		# these may vary from run to run
		self.wlistrun.clear()
		for k in self.jsparm.keys():
			if k[:4]!="run_": continue
			self.wlistrun.addItem(k[4:])
		self.sel_run(self.wlistrun.count()-1)
			
	def sel_run(self,line):
		"""Called when the user selects a new run from the list"""
		
#		print("sel_run",line)
		if line<0: return
		self.currunkey=str(self.wlistrun.item(line).text())
		self.currun=self.jsparm.getdefault("run_"+self.currunkey,{"dim":4,"mask":f"{self.gmm}/mask.hdf","trainiter":10,"pas":"100","time":"-"})
		self.wedngauss.setText(f'{self.currun.get("ngauss",256)}')
		self.weddim.setText(f'{self.currun.get("dim",4)}')
		self.wedmask.setText(f'{self.currun.get("mask",self.jsparm.getdefault("mask",f"{self.gmm}/mask.hdf"))}')
		self.wedtrainiter.setText(f'{self.currun.get("trainiter",10)}')
		pas=self.currun.get("pas","100")
		self.wbutpos.setChecked(int(pas[0]))
		self.wbutamp.setChecked(int(pas[1]))
		self.wbutsig.setChecked(int(pas[2]))
		self.wlabruntime.setText(self.currun.get("time","-"))
		
		# Decoder model for generating Gaussians
		try: self.decoder = tf.keras.models.load_model(f"{self.gmm}/{self.currunkey}_decoder.h5")
		except: 
			showerror(f"Run {self.currunkey} results incomplete. No stored decoder found.",self)
			return

		# Middle layer for every particle
		
		self.midresult=np.loadtxt(f"{self.gmm}/{self.currunkey}_mid.txt")[:,1:].transpose()
		self.wbutdrmid.click()
		

	def add_gmm(self,clk=False):
		"""Creates a new numbered gmm_XX folder"""
		self.update_gmms()
		if len(self.gmms)==0 : crt=1
		else: crt=int(self.gmms[-1][-2:])+1
		newgmm=f"gmm_{crt:02d}"
		try: os.mkdir(newgmm)
		except:
			showerror(f"Cannot create {newgmm}",self)
			return
		self.gmm=newgmm
		self.setgmm_refine()
		self.wlistgmm.addItem(newgmm)
		
	def setgmm_refine(self,clk=False):
		"""Allows the user to base the refine_XX folder for the currently selected gmm_XX folder""" 
		self.jsparm=js_open_dict(f"{self.gmm}/0_gmm_parms.json")
		
		# double check if the user will be destroying results
		if self.jsparm.has_key("refinepath"):
			ans=QtWidgets.QMessageBox.question(self,"Are you sure?",f"{self.gmm} has already been configured to work on {self.jsparm['refinepath']}. Continuing may invalidate current results. Proceed?")
			if ans==QtWidgets.QMessageBox.No: return
			
		# Get the name of an existing refinement
		rpath=str(QtWidgets.QFileDialog.getExistingDirectory(self,"Please select an existing refine_xx folder to seed the analysis"))
		if not os.path.isdir(rpath) : 
			showerror("Invalid path")
			return
		self.jsparm["refinepath"]=rpath
		
		### setup the folder
		try:
			itr=max([int(i.split("_")[1]) for i in os.listdir(rpath) if i[:12]=="projections_" and i.split("_")[1].isdigit()])
		except:
			showerror("No projections in refine folder")
			return
		
		# Copy projections from refine folder
		eprj=f"{rpath}/projections_{itr:02d}_even.hdf"
		oprj=f"{rpath}/projections_{itr:02d}_odd.hdf"
		n=EMUtil.get_image_count(eprj)
		for i in range(n):
			a1=EMData(eprj,i)
			a2=EMData(oprj,i)
			a=a1+a2
			a.mult(0.5)
			a.write_compressed(f"{self.gmm}/proj_in.hdf",i,10)
		self.jsparm["boxsize"]=a["nx"]
		self.jsparm["apix"]=a["apix_x"]
#		self.jsparm["res"]=file_resolution(f"{rpath}/fsc_maskedtight_{itr:02d}.txt")

		# Copy mask from refine folder
		a=EMData(f"{rpath}/mask_tight.hdf")
		a.write_compressed(f"{self.gmm}/mask.hdf",0,8)
		self.jsparm["mask"]=f"{self.gmm}/mask.hdf"

		# Extract particles from refine folder
		run(f"e2evalrefine.py {rpath} --extractorientptcl {self.gmm}/particles.lst")
		
	
	def closeEvent(self,event):
		E2saveappwin("e2gmm","main",self)


if __name__ == "__main__":
	main()
