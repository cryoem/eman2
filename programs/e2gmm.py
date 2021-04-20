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
from eman2_gui.emdataitem3d import EMDataItem3D,EMIsosurface
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
	parser.add_argument("--threads", default=-1,type=int,help="Number of alignment threads to run in parallel on a single computer.")
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
	QtWidgets.QMessageBox.warning(None,"Warning","""This program is still experimental. While functional
many capabilities of the underlying e2gmm_refine program
are not yet available through this interface. Matching number
of gaussians to resolution/volume is key in obtaining good
distributions.""")
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
		self.gbl.addWidget(self.wview3d,0,4,3,1)
		
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
		self.gblrun.addWidget(self.wbutrerun,1,1)
		
		self.wbutnewrun=QtWidgets.QPushButton("New Run")
		self.gblrun.addWidget(self.wbutnewrun,1,0)
		
		# The form with details about the selected gmm_XX folder
		self.gflparm = QtWidgets.QFormLayout()
		self.gbll.addLayout(self.gflparm,2,0)
		
		self.wlpath = QtWidgets.QLabel("-")
		self.wlpath.setToolTip("Path this GMM is based on")
		self.gflparm.addRow("Path:",self.wlpath)
		
		self.wedbox = QtWidgets.QLineEdit("256")
		self.wedbox.setToolTip("Box size of input particles in pixels")
		self.gflparm.addRow("Box Size:",self.wedbox)

		self.wedres = QtWidgets.QLineEdit("25")
		self.wedres.setToolTip("Maximum resolution to use for gaussian fitting")
		self.gflparm.addRow("Target Res:",self.wedres)

		self.wedsym = QtWidgets.QLineEdit("c1")
		self.wedsym.setToolTip("Symmetry used during refinement")
		self.gflparm.addRow("Symmetry:",self.wedsym)
		
		self.wedapix = QtWidgets.QLineEdit("1.0")
		self.wedapix.setToolTip("A/pix of input particles")
		self.gflparm.addRow("A/pix:",self.wedapix)

		self.wedmask = QtWidgets.QLineEdit("")
		self.wedmask.setToolTip("3-D volume mask")
		self.gflparm.addRow("Mask:",self.wedmask)
		
		self.wedngauss = QtWidgets.QLineEdit("64")
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

		# Widgets below 3D
		self.gbl3dctl = QtWidgets.QGridLayout()
		self.gbl.addLayout(self.gbl3dctl,3,4)
		
		#self.wbutmap=QtWidgets.QPushButton("Map")
		#self.wbutmap.setCheckable(True)
		#self.wbutmap.setChecked(True)
		#self.gbl3dctl.addWidget(self.wbutmap,0,0)

		#self.wbutspheres=QtWidgets.QPushButton("Sphere Mdl")
		#self.wbutspheres.setCheckable(True)
		#self.wbutspheres.setChecked(True)
		#self.gbl3dctl.addWidget(self.wbutspheres,0,1)

		self.wvssphsz=ValSlider(self,(1,50),"Size:",3.0,90)
		self.gbl3dctl.addWidget(self.wvssphsz,0,2)
		
		# Connections
		self.wlistgmm.currentRowChanged[int].connect(self.sel_gmm)
		#self.wbutspheres.clicked[bool].connect(self.new_3d_opt)
		#self.wbutmap.clicked[bool].connect(self.new_3d_opt)
		self.wvssphsz.valueChanged.connect(self.new_sph_size)
		self.wbutnewgmm.clicked[bool].connect(self.add_gmm)
		self.wbutrefine.clicked[bool].connect(self.setgmm_refine)
		self.wlistrun.currentRowChanged[int].connect(self.sel_run)
		self.wbutnewrun.clicked[bool].connect(self.new_run)
		self.wbutrerun.clicked[bool].connect(self.do_run)
		self.wbutdrgrp.buttonClicked[QtWidgets.QAbstractButton].connect(self.plot_mode_sel)
		self.wsbxcol.valueChanged[int].connect(self.wplot2d.setXAxisAll)
		self.wsbycol.valueChanged[int].connect(self.wplot2d.setYAxisAll)
		self.wplot2d.mousedown[QtGui.QMouseEvent,tuple].connect(self.plot_mouse)
		self.wplot2d.mouseup[QtGui.QMouseEvent,tuple].connect(self.plot_mouse)
		self.wplot2d.mousedrag[QtGui.QMouseEvent,tuple].connect(self.plot_mouse)
		E2loadappwin("e2gmm","main",self)

		self.gaussplot=EMScatterPlot3D()
		self.mapdataitem=EMDataItem3D(None)
		self.mapiso=EMIsosurface(self.mapdataitem)
		self.wview3d.insertNewNode("Neutral Map",self.mapdataitem)
		self.wview3d.insertNewNode("Isosurface",self.mapiso,parentnode=self.mapdataitem)
		self.wview3d.insertNewNode("Gauss Model",self.gaussplot)

		#QtCore.QTimer.singleShot(500,self.afterStart)

	def do_events(self,delay=0.1):
		"""process the event loop with a small delay to allow user abort, etc."""
		t=time.time()
		while (time.time()-t<delay): 
			self.app().processEvents()
	
	def new_sph_size(self,newval=10):
		self.gaussplot.setPointSize(newval)
		self.wview3d.update()
	
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
		self.gaussplot.setData(gauss,self.wvssphsz.value)
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

	def new_3d_opt(self,clk=False):
		"""When the user changes selections for the 3-D display"""
		self.gaussplot.setVisibleItem(butval(self.wbutspheres))
		self.mapdataitem.setVisibleItem(butval(self.wbutmap))
		print(self.gaussplot.isVisibleItem(),self.mapdataitem.isVisibleItem())

	def new_run(self,clk=False):
		"""Create a new run and run() it"""
		name=str(QtWidgets.QInputDialog.getText(self,"Run Name","Enter a name for this run. Current parameters will be used.")[0])
		if not self.jsparm.has_key("run_"+name) : self.wlistrun.addItem(name)
		self.currunkey=name
		self.do_run()

	def do_run(self,clk=False):
		"""Run the current job with current parameters"""
		self.currun={}
		self.currun["boxsize"]=int(self.wedbox.text())
		self.currun["targres"]=float(self.wedres.text())
		self.currun["apix"]=float(self.wedapix.text())
		self.currun["sym"]=str(self.wedsym.text())
		self.currun["ngauss"]=int(self.wedngauss.text())
		self.currun["dim"]=int(self.weddim.text())
		self.currun["mask"]=str(self.wedmask.text())
		self.currun["trainiter"]=int(self.wedtrainiter.text())
		self.currun["pas"]=butval(self.wbutpos)+butval(self.wbutamp)+butval(self.wbutsig)
		self.currun["time"]=local_datetime()
		self.jsparm["run_"+self.currunkey]=self.currun
		
		maxbox=(int(self.currun["boxsize"]*(2*self.currun["apix"])/self.currun["targres"])//2)*2
		print(f"Target res {self.currun['targres']} -> max box size {maxbox}")
		modelout=f"{self.gmm}/{self.currunkey}_model_gmm.txt"
		modelseg=f"{self.gmm}/{self.currunkey}_model_seg.txt"
		
		sym=self.currun["sym"]
		prog=QtWidgets.QProgressDialog("Running networks. Progress updates here are limited. See the Console for detailed output.","Abort",0,3)
		prog.show()
		self.do_events(1)
		curngauss=self.currun["ngauss"]//2**(int(log(maxbox/16)/log(2.0))+1)
		
		#### Original method, pure network approach
		## First step with very coarse model, gradually increasing size improves convergence
		#run(f"e2gmm_refine.py --projs {self.gmm}/proj_in.hdf --npt {curngauss} --sym {sym} --maxboxsz 16 --modelout {modelout} --niter {self.currun['trainiter']*2} --mask {self.currun['mask']} --nmid {self.currun['dim']}")
		#prog.setValue(1)
		#self.do_events()
		#if prog.wasCanceled() : return
		
		#box=16
		#n=2
		#while box<maxbox:
			#box=good_size(box*2)
			#curngauss*=2
			#cungauss=min(curngauss,self.currun["ngauss"])

			## in the last iteration we do some additional things
			#if box>=maxbox:
				#box=maxbox
				#s=f"--evalmodel {self.gmm}/{self.currunkey}_model_projs.hdf --evalsize {self.jsparm['boxsize']}"		# in the final iteration
			#else: s=""
				
			## iterate until box size is the full size of the particles
			#er=run(f"e2gmm_refine.py --projs {self.gmm}/proj_in.hdf --npt {curngauss} --sym {sym} --maxboxsz {box} --model {modelout} --modelout {modelout} --niter {self.currun['trainiter']} --mask {self.currun['mask']} --nmid {self.currun['dim']} {s}")
			#if er :
				#showerror("Error running e2gmm_refine, see console for details. GPU memory exhaustion is a common issue. Consider reducing the target resolution.")
				#return
			#prog.setValue(n)
			#self.do_events()
			#if prog.wasCanceled() : return
			#n+=1
		#### New method using segmentation for initial seed
		print("segmentation to initialize")
		inmap=EMData(f"{self.gmm}/input_map.hdf")
		inmask=EMData(self.currun["mask"])
		inmap.mult(inmask)
		std=inmap["sigma_nonzero"]
		#seg=inmap.process("segment.kmeans",{"nseg":self.currun["ngauss"]/4,"thr":std,"maxiter":40,"verbose":1})
		apix=self.currun["apix"]
		seg=inmap.process("segment.distance",{"minsegsep":self.currun['targres']/apix*0.7,"maxsegsep":self.currun['targres']/apix*1.3,"thr":std})
		nx=seg["nx"]
		coord=seg["segment_centers"]
		ng=len(coord)//3
		with open(modelseg,"w") as mdl:
			#for i in range(self.currun["ngauss"]//4):
			for i in range(ng):
				mdl.write(f"{coord[i*3]/nx-0.5:0.4f}\t{coord[i*3+1]/nx-0.5:0.4f}\t{coord[i*3+2]/nx-0.5:0.4f}\t1.0\t1.0\n")
			
			coord=None
			seg=None
		
		print(ng," Gaussian seeds")
		prog.setValue(1)
		if prog.wasCanceled() : return

		er=run(f"e2gmm_refine.py --projs {self.gmm}/proj_in.hdf --npt {max(self.currun['ngauss'],ng)} --sym {sym} --maxboxsz {maxbox} --model {modelseg} --modelout {modelout} --niter {self.currun['trainiter']} --mask {self.currun['mask']} --nmid {self.currun['dim']} --evalmodel {self.gmm}/{self.currunkey}_model_projs.hdf --evalsize {self.jsparm['boxsize']}")
		if er :
			showerror("Error running e2gmm_refine, see console for details. GPU memory exhaustion is a common issue. Consider reducing the target resolution.")
			return


		# make3d on gaussian output for comparison
		er=run(f"e2make3dpar.py --input {self.gmm}/{self.currunkey}_model_projs.hdf --output {self.gmm}/{self.currunkey}_model_recon.hdf --pad {good_size(self.jsparm['boxsize']*1.25)} --mode trilinear --keep 1 --threads {self.options.threads}")
		prog.setValue(2)
		self.do_events()
		if prog.wasCanceled() : return

		# heterogeneity analysis
		er=run(f"e2gmm_refine.py --model {modelout} --ptclsin {self.gmm}/particles.lst --heter --sym {sym} --maxboxsz {maxbox} --gradout {self.gmm}/{self.currunkey}_grads.hdf --mask {self.currun['mask']} --nmid {self.currun['dim']} --midout {self.gmm}/{self.currunkey}_mid.txt --decoderout {self.gmm}/{self.currunkey}_decoder.h5 --pas {self.currun['pas']}")
		if er :
			showerror("Error running e2gmm_refine, see console for details. Memory is a common issue. Consider reducing the target resolution.")
			return
		prog.setValue(3)
		self.do_events()
		
		self.sel_run(0)

	def update_gmms(self):
		"""Updates the display of gmm_XX folders"""
		#self.gmm=str(self.wlistgmm.currentItem().text())
		self.gmms=[i for i in sorted(os.listdir(".")) if i[:4]=="gmm_" and os.path.isdir(i)]
		self.wlistgmm.clear()
		for i in self.gmms:
			self.wlistgmm.addItem(i)
	
	def sel_gmm(self,line):
		"""Called when the user selects a new GMM from the list. Should not be called directly by the user without updating wlistgmm"""

#		print("sel gmm",line)
		if line<0 : 
			self.wedbox.setText("")
			self.wedapix.setText("")
			self.wedsym.setText("c1")
			self.wlistrun.clear()
			return
		self.gmm=str(self.wlistgmm.item(line).text())
		self.jsparm=js_open_dict(f"{self.gmm}/0_gmm_parms.json")
		# These are associated with the whole GMM
		self.wlpath.setText(f'{self.jsparm.getdefault("refinepath","-")}')
		self.wedbox.setText(f'{self.jsparm.getdefault("boxsize",128)}')
		self.wedapix.setText(f'{self.jsparm.getdefault("apix","")}')
		self.wedsym.setText(f'{self.jsparm.getdefault("sym","c1")}')
		self.wedmask.setText(f'{self.jsparm.getdefault("mask",f"{self.gmm}/mask.hdf")}')
		
		# these may vary from run to run
		self.wlistrun.clear()
		for k in self.jsparm.keys():
			if k[:4]!="run_": continue
			self.wlistrun.addItem(k[4:])
		self.sel_run(self.wlistrun.count()-1)
		
		try: self.mapdataitem.setData(EMData(f'{self.gmm}/input_map.hdf'))
		except: print("Error: input_map.hdf missing")
			
	def sel_run(self,line):
		"""Called when the user selects a new run from the list"""
		
#		print("sel_run",line)
		if line<0: return
		self.currunkey=str(self.wlistrun.item(line).text())
		self.currun=self.jsparm.getdefault("run_"+self.currunkey,{"dim":4,"mask":f"{self.gmm}/mask.hdf","trainiter":10,"pas":"100","time":"-"})
		self.wedres.setText(f'{self.currun.get("targres",20)}')
		self.wedapix.setText(f'{self.currun.get("apix",self.jsparm.getdefault("apix",""))}')
		self.wedngauss.setText(f'{self.currun.get("ngauss",64)}')
		self.weddim.setText(f'{self.currun.get("dim",4)}')
		self.wedsym.setText(f'{self.currun.get("sym","c1")}')
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
		
		self.plot_mouse(None,(0,0))
		

	def add_gmm(self,clk=False):
		"""Creates a new numbered gmm_XX folder"""
		try: newgmm=num_path_new("gmm")
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
		try:
			rpath=os.path.relpath(str(QtWidgets.QFileDialog.getExistingDirectory(self,"Please select an existing refine_xx folder to seed the analysis")))
		except: return
		if not os.path.isdir(rpath) : 
			showerror("Invalid path")
			return
		self.jsparm["refinepath"]=rpath
		
		### setup the folder
		try:
			itr=max([int(i.split("_")[1]) for i in os.listdir(rpath) if i[:7]=="threed_" and i.split("_")[1].isdigit()])
		except:
			showerror("No projections in refine folder")
			return

		self.app().setOverrideCursor(Qt.BusyCursor)
		rparm=js_open_dict(f"{rpath}/0_refine_parms.json")
		if rparm["breaksym"] : self.jsparm["sym"]="c1"
		else: self.jsparm["sym"]=rparm["sym"]
		
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

		# Copy map from refine folder
		a=EMData(f"{rpath}/threed_{itr:02d}.hdf")
		a.write_compressed(f"{self.gmm}/input_map.hdf",0,12)
		self.jsparm["source_map"]=f"{rpath}/threed_{itr:02d}.hdf"
		
		# Copy mask from refine folder
		a=EMData(f"{rpath}/mask_tight.hdf")
		a.write_compressed(f"{self.gmm}/mask.hdf",0,8)
		self.jsparm["mask"]=f"{self.gmm}/mask.hdf"

		# Extract particles from refine folder
		run(f"e2evalrefine.py {rpath} --extractorientptcl {self.gmm}/particles.lst")
		self.app().setOverrideCursor(Qt.ArrowCursor)		
		
	
	def closeEvent(self,event):
		E2saveappwin("e2gmm","main",self)


if __name__ == "__main__":
	main()
