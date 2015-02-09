#!/usr/bin/env python

#
# Author: Steven Ludtke, 07/28/2011 (sludtke@bcm.edu)
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
from EMAN2db import db_open_dict, db_close_dict, db_check_dict, db_list_dicts
from OpenGL import GL,GLUT
from math import *
import os
import sys
import weakref
import e2ctf
import threading
from numpy import array,arange
import traceback

try:
	from PyQt4 import QtCore, QtGui, QtOpenGL
	from PyQt4.QtCore import Qt
	from PyQt4.QtCore import QTimer
	from emshape import *
	from valslider import *
except:
	print "Warning: PyQt4 must be installed"
	sys.exit(1)

from Simplex import Simplex



debug=False
logid=None

sfcurve=None		# This will store a global structure factor curve if specified
sfcurve2=None
envelopes=[]		# simplex minimizer needs to use a global at the moment

def main():
	global debug,logid
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <single image file> ...

This program will allow you to evaluate an individual scanned micrograph or CCD frame, by looking at its
power spectrum in various ways."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_pos_argument(name="image",help="Image to process with e2evalimage.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True,multiselect=True)",  row=0, col=0,rowspan=1, colspan=2, mode="eval")
	parser.add_header(name="evalimageheader", help='Options below this label are specific to e2evalimage', title="### e2evalimage options ###", default=None, row=1, col=0, rowspan=1, colspan=2, mode="eval")

	parser.add_argument("--gui",action="store_true",help="This is a GUI-only program. This option is provided for self-consistency",default=True)
	parser.add_argument("--apix",type=float,help="Angstroms per pixel for all images",default=None, guitype='floatbox', row=3, col=0, rowspan=1, colspan=1, mode="eval['self.pm().getAPIX()']")
	parser.add_argument("--constbfactor",type=float,help="Set B-factor to fixed specified value, negative value autofits",default=-1.0, guitype='floatbox', row=8, col=0, rowspan=1, colspan=1, mode='eval[-1.0]')
	parser.add_argument("--voltage",type=float,help="Microscope voltage in KV",default=None, guitype='floatbox', row=3, col=1, rowspan=1, colspan=1, mode="eval['self.pm().getVoltage()']")
	parser.add_argument("--cs",type=float,help="Microscope Cs (spherical aberation)",default=None, guitype='floatbox', row=4, col=0, rowspan=1, colspan=1, mode="eval['self.pm().getCS()']")
	parser.add_argument("--ac",type=float,help="Amplitude contrast (percentage, default=10)",default=10, guitype='floatbox', row=4, col=1, rowspan=1, colspan=1, mode="eval")
	parser.add_argument("--astigmatism",action="store_true",help="Includes astigmatism in automatic fitting",default=False, guitype='boolbox', row=8, col=1, rowspan=1, colspan=1, mode='eval')
	parser.add_argument("--box",type=int,help="Forced box size in grid mode. Overrides any previous setting. ",default=-1, guitype='intbox', row=5, col=0, rowspan=1, colspan=1, mode="eval")
	parser.add_argument("--usefoldername",action="store_true",help="If you have the same image filename in multiple folders, and need to import into the same project, this will prepend the folder name on each image name",default=False,guitype='boolbox',row=5, col=1, rowspan=1, colspan=1, mode="eval")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()

	logid=E2init(sys.argv,options.ppid)

	from emapplication import EMApp
	app=EMApp()
	gui=GUIEvalImage(args,options.voltage,options.apix,options.cs,options.ac,options.box,options.usefoldername,options.constbfactor,options.astigmatism)
	gui.show()

	try:
		gui.wimage.raise_()
		gui.wfft.raise_()
		gui.wplot.raise_()
		gui.raise_()
	except: pass

# 	try: gui.raise_()
# 	except: pass
	app.execute()

	E2end(logid)



class GUIEvalImage(QtGui.QWidget):
	def __init__(self,images,voltage=None,apix=None,cs=None,ac=10.0,box=512,usefoldername=False,constbfactor=-1,fitastig=False):
		"""Implements the CTF fitting dialog using various EMImage and EMPlot2D widgets
		'data' is a list of (filename,ctf,im_1d,bg_1d,quality)
		'parms' is [box size,ctf,box coord,set of excluded boxnums,quality,oversampling]
		"""
		try:
			from emimage2d import EMImage2DWidget
		except:
			print "Cannot import EMAN image GUI objects (EMImage2DWidget)"
			sys.exit(1)
		try:
			from emplot2d import EMPlot2DWidget
		except:
			print "Cannot import EMAN plot GUI objects (is matplotlib installed?)"
			sys.exit(1)

		QtGui.QWidget.__init__(self,None)
		self.setWindowIcon(QtGui.QIcon(get_image_directory() + "ctf.png"))

		self.fitastig=fitastig
		self.nodir=not usefoldername
		self.data=None
		self.curset=0
		self.calcmode=1
		self.f2dmode=0
		self.f2danmode=0
		self.plotmode=0
		self.curfilename = None

		self.ringrad=1
		self.xpos1=(10,0)
		self.xpos2=(0,10)
#		self.db = db_open_dict('bdb:mgquality')

		self.defaultvoltage=voltage
		self.defaultapix=apix
		self.defaultcs=cs
		self.defaultac=ac
		self.constbfactor=constbfactor

		# Per image parameters to keep track of
		# for each image [box size,ctf,box coord,set of excluded boxnums,quality,oversampling]
		self.parms=[]
#		db_fparms=db_open_dict("bdb:e2ctf.frameparms",True)

		# This little block of code deals with Z stacks of images and multiple image files,
		# rewriting the name(s) with an enumeration including ;number or ,number
		newimages=[]
		for i in images:
			n=EMUtil.get_image_count(i)
			if n!=1:
				for j in xrange(n): newimages.append(i+",%d"%j)
			else:
				h=EMData(i,0,True)		# read header
				n=h["nz"]
				if n!=1:
					for j in xrange(n): newimages.append(i+";%d"%j)
				else : newimages.append(i)
		images=newimages

		# Now we try to restore old image information
		for i in images:
			try:
				pd=js_open_dict(info_name(i,nodir=self.nodir))
				parms=pd["ctf_frame"]
				parms[4]=pd["quality"]
				if box>=64 and parms[0]!=box:
					parms[0]=box
					parms[3]=set()
#				if parms==None : raise Exception
			except:
				ctf = EMAN2Ctf()
				ctf.from_dict({'defocus':0.0,'dfdiff':0.0,'dfang':0.0,'bfactor':200.0,'ampcont':self.defaultac,'voltage':200.0,'cs':4.1,'apix':1.0,'dsbg':-1})
				if self.defaultvoltage!=None : ctf.voltage=self.defaultvoltage
				if self.defaultcs!=None : ctf.cs=self.defaultcs
				if self.defaultapix!=None : ctf.apix=self.defaultapix
				parms=[int(box),ctf,(256,256),set(),5,1]
				print "Initialize new parms for: ",base_name(i)

			if self.constbfactor>0 : parms[1].bfactor=self.constbfactor

			if parms[0]<64 :
				parms[0]=512
				parms[3]=set()

			self.parms.append(parms)

		self.wimage=EMImage2DWidget()
		self.wimage.setWindowTitle("e2evalimage - Micrograph")

		self.wfft=EMImage2DWidget()
		self.wfft.setWindowTitle("e2evalimage - 2D FFT")

		self.wplot=EMPlot2DWidget()
		self.wplot.setWindowTitle("e2evalimage - Plot")


		self.wimage.connect(self.wimage,QtCore.SIGNAL("mousedown"),self.imgmousedown)
		self.wimage.connect(self.wimage,QtCore.SIGNAL("mousedrag"),self.imgmousedrag)
		self.wimage.connect(self.wimage,QtCore.SIGNAL("mouseup")  ,self.imgmouseup)
		self.wfft.connect(self.wfft,QtCore.SIGNAL("mousedown"),self.fftmousedown)
		self.wfft.connect(self.wfft,QtCore.SIGNAL("mousedrag"),self.fftmousedrag)
		self.wfft.connect(self.wfft,QtCore.SIGNAL("mouseup")  ,self.fftmouseup)
		self.wplot.connect(self.wplot,QtCore.SIGNAL("mousedown"),self.plotmousedown)

		self.wimage.mmode="app"
		self.wfft.mmode="app"

		# This object is itself a widget we need to set up
		self.gbl = QtGui.QGridLayout(self)
		self.gbl.setMargin(8)
		self.gbl.setSpacing(6)

		# plot list and plot mode combobox
		self.setlist=e2ctf.MyListWidget(self)
		self.setlist.setSizePolicy(QtGui.QSizePolicy.Preferred,QtGui.QSizePolicy.Expanding)
		for i in images:
			self.setlist.addItem(i)
		self.gbl.addWidget(self.setlist,0,0,10,2)

		self.lcalcmode=QtGui.QLabel("Region:",self)
		self.gbl.addWidget(self.lcalcmode,10,0)

		self.scalcmode=QtGui.QComboBox(self)
		self.scalcmode.addItem("Single Region")
		self.scalcmode.addItem("Tiled Boxes")
		self.scalcmode.setCurrentIndex(1)
		self.gbl.addWidget(self.scalcmode,10,1)


		self.lcalcmode=QtGui.QLabel("2D FFT:",self)
		self.gbl.addWidget(self.lcalcmode,11,0)

		self.s2dmode=QtGui.QComboBox(self)
		self.s2dmode.addItem("Power Spectrum")
		self.s2dmode.addItem("Bg Subtracted")
		self.s2dmode.addItem("Background")
		self.gbl.addWidget(self.s2dmode,11,1)

		self.lcalcmode=QtGui.QLabel("Annotate:",self)
		self.gbl.addWidget(self.lcalcmode,12,0)

		self.s2danmode=QtGui.QComboBox(self)
		self.s2danmode.addItem("Ctf Zeroes")
		self.s2danmode.addItem("Resolution Ring")
		self.s2danmode.addItem("2-D Xtal")
		self.s2danmode.addItem("None")
		self.gbl.addWidget(self.s2danmode,12,1)


		self.lcalcmode=QtGui.QLabel("Plot:",self)
		self.gbl.addWidget(self.lcalcmode,13,0)

		self.splotmode=QtGui.QComboBox(self)
		self.splotmode.addItem("Bgsub and Fit")
		self.splotmode.addItem("Fg and Bg")
		self.splotmode.addItem("Bgsub, 45 deg slices")
		self.splotmode.addItem("Fg, 45 deg slices")
		self.splotmode.addItem("SSNR (NOT BY PTCL !)")
		self.gbl.addWidget(self.splotmode,13,1)

		# ValSliders for CTF parameters

		#self.samp = ValSlider(self,(0,5.0),"Amp:",0)
		#self.vbl.addWidget(self.samp)

		self.sdefocus=ValSlider(self,(0,5),"Defocus:",0.0,90)
		self.gbl.addWidget(self.sdefocus,0,2,1,3)

		self.sbfactor=ValSlider(self,(0,1600),"B factor:",200.0,90)
		self.gbl.addWidget(self.sbfactor,1,2,1,3)

		self.sdfdiff=ValSlider(self,(0,1),"DF Diff:",0,90)
		self.gbl.addWidget(self.sdfdiff,2,2,1,3)

		self.sdfang=ValSlider(self,(0,180),"Df Angle:",0,90)
		self.gbl.addWidget(self.sdfang,3,2,1,3)

		self.sampcont=ValSlider(self,(0,100),"% AC",10.0,90)
		if self.defaultac!=None : self.sampcont.setValue(self.defaultac)
		self.gbl.addWidget(self.sampcont,4,2,1,3)

		self.sang45=ValSlider(self,(-22.5,22.5),"45 mode ang",0.0,90)
		self.gbl.addWidget(self.sang45,5,2,1,3)

		self.squality=ValSlider(self,(0,9),"Quality (0-9):",0,90)
		self.squality.setIntonly(True)
		self.gbl.addWidget(self.squality,6,2,1,3)

		self.brefit=QtGui.QPushButton("Refit")
		self.gbl.addWidget(self.brefit,7,2)

		self.cbgadj=CheckBox(None,"CTF BG Adj",1)
		self.gbl.addWidget(self.cbgadj,7,3)


#		self.sapix=ValSlider(self,(.2,10),"A/Pix:",2,90)
#		self.vbl.addWidget(self.sapix)

		self.sapix=ValBox(self,(0,500),"A/pix:",1.0,90)
		if self.defaultapix!=None : self.sapix.setValue(self.defaultapix)
		self.gbl.addWidget(self.sapix,10,2)

		self.svoltage=ValBox(self,(0,500),"Voltage (kV):",200,90)
		if self.defaultvoltage!=None : self.svoltage.setValue(self.defaultvoltage)
		self.gbl.addWidget(self.svoltage,11,2)

		self.scs=ValBox(self,(0,5),"Cs (mm):",4.1,90)
		if self.defaultcs!=None : self.scs.setValue(self.defaultcs)
		self.gbl.addWidget(self.scs,12,2)

		self.sboxsize=ValBox(self,(128,1024),"Box Size:",512,90)
		self.sboxsize.intonly=True
		self.gbl.addWidget(self.sboxsize,13,2)

#		self.soversamp=ValBox(self,(0,500),"Oversample:",1,90)
#		self.soversamp.intonly=True
#		self.gbl.addWidget(self.soversamp,14,2)					# disable oversampling for now, not properly integrated with save/recall mechanism, and not clear if it's really useful


		# this is just a spacer
		self.gbl.setColumnStretch(3,2)

		self.bbox=QtGui.QGroupBox("Project")
		self.gbl.addWidget(self.bbox,10,4,4,1)

		self.bvbl=QtGui.QVBoxLayout()
		self.bbox.setLayout(self.bvbl)

		self.bimport=QtGui.QPushButton("Import")
		self.bvbl.addWidget(self.bimport)

		self.cinvert=CheckBox(None,"Invert")
		self.bvbl.addWidget(self.cinvert)

		self.cxray=CheckBox(None,"X-ray Pixels")
		self.bvbl.addWidget(self.cxray)

		QtCore.QObject.connect(self.bimport, QtCore.SIGNAL("clicked(bool)"),self.doImport)
		QtCore.QObject.connect(self.brefit, QtCore.SIGNAL("clicked(bool)"),self.doRefit)
		QtCore.QObject.connect(self.cbgadj, QtCore.SIGNAL("valueChanged"),self.bgAdj)
		QtCore.QObject.connect(self.sdefocus, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.sbfactor, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.sdfdiff, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.sdfang, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.sapix, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.sampcont, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.svoltage, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.scs, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.sboxsize, QtCore.SIGNAL("valueChanged"), self.newBox)
#		QtCore.QObject.connect(self.soversamp, QtCore.SIGNAL("valueChanged"), self.newBox)
		QtCore.QObject.connect(self.sang45, QtCore.SIGNAL("valueChanged"), self.recalc_real)
		QtCore.QObject.connect(self.squality,QtCore.SIGNAL("valueChanged"),self.newQualityFactor)
		QtCore.QObject.connect(self.setlist,QtCore.SIGNAL("currentRowChanged(int)"),self.newSet)
		QtCore.QObject.connect(self.setlist,QtCore.SIGNAL("keypress"),self.listkey)
		QtCore.QObject.connect(self.scalcmode,QtCore.SIGNAL("currentIndexChanged(int)"),self.newCalcMode)
		QtCore.QObject.connect(self.s2dmode,QtCore.SIGNAL("currentIndexChanged(int)"),self.new2DMode)
		QtCore.QObject.connect(self.s2danmode,QtCore.SIGNAL("currentIndexChanged(int)"),self.new2DAnMode)
		QtCore.QObject.connect(self.splotmode,QtCore.SIGNAL("currentIndexChanged(int)"),self.newPlotMode)

	   	#QtCore.QObject.connect(self.saveparms,QtCore.SIGNAL("clicked(bool)"),self.on_save_params)
		#QtCore.QObject.connect(self.recallparms,QtCore.SIGNAL("clicked(bool)"),self.on_recall_params)
		#QtCore.QObject.connect(self.refit,QtCore.SIGNAL("clicked(bool)"),self.on_refit)
		#QtCore.QObject.connect(self.output,QtCore.SIGNAL("clicked(bool)"),self.on_output)


		self.resize(720,380) # figured these values out by printing the width and height in resize event

		### This section is responsible for background updates
		self.busy=False
		self.needupdate=True
		self.needredisp=False
		self.procthread=None
		self.errors=None		# used to communicate errors back from the reprocessing thread

		self.timer=QTimer()
		QtCore.QObject.connect(self.timer, QtCore.SIGNAL("timeout()"), self.timeOut)
		self.timer.start(100)

		self.setWindowTitle("e2evalimage - Control Panel")

		self.wimage.show()
		self.wfft.show()
		self.wplot.show()
		E2loadappwin("e2evalimage","main",self)
		E2loadappwin("e2evalimage","image",self.wimage.qt_parent)
		E2loadappwin("e2evalimage","fft",self.wfft.qt_parent)
		E2loadappwin("e2evalimage","plot",self.wplot.qt_parent)
#		self.recalc()

	def listkey(self,event):

		if event.key()>=Qt.Key_0 and event.key()<=Qt.Key_9 :
			q=int(event.key())-Qt.Key_0
			self.squality.setValue(q)
		elif event.key() == Qt.Key_Left:
			self.sdefocus.setValue(self.sdefocus.getValue()-0.03)
		elif event.key() == Qt.Key_Right:
			self.sdefocus.setValue(self.sdefocus.getValue()+0.03)
		elif event.key()==Qt.Key_I :
			self.doImport()
		elif event.key()==Qt.Key_U :
			self.unImport()
		elif event.key()==Qt.Key_F :
			self.doRefit()


	def closeEvent(self,event):
#		QtGui.QWidget.closeEvent(self,event)
		E2saveappwin("e2evalimage","main",self)
		E2saveappwin("e2evalimage","image",self.wimage.qt_parent)
		E2saveappwin("e2evalimage","fft",self.wfft.qt_parent)
		E2saveappwin("e2evalimage","plot",self.wplot.qt_parent)

		self.writeCurParm()
		event.accept()
		QtGui.qApp.exit(0)
		#app=QtGui.qApp
		#if self.wimage != None:
			#app.close_specific(self.wimage)
			#self.wimage = None
		#if self.wfft != None:
			#app.close_specific(self.wfft)
		#if self.wplot != None:
			#app.close_specific(self.wplot)
		#app.close_specific(self)
#		self.emit(QtCore.SIGNAL("module_closed")) # this signal is important when e2ctf is being used by a program running its own event loop

	def update_plot(self):
#		if self.wplot == None: return # it's closed/not visible

		if self.incalc : return		# no plot updates during a recomputation

		parms=self.parms[self.curset]
		apix=self.sapix.getValue()
		ds=1.0/(apix*parms[0]*parms[5])
		ctf=parms[1]
		bg1d=array(ctf.background)
		r=len(ctf.background)
		s=arange(0,ds*r,ds)

		# This updates the FFT CTF-zero circles
		if self.f2danmode==0 :
			fit=ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP)
			shp={}
			nz=0
			for i in range(1,len(fit)):
				if fit[i-1]*fit[i]<=0.0:
					nz+=1
					shp["z%d"%i]=EMShape(("circle",0.0,0.0,1.0/nz,r,r,i,1.0))

			self.wfft.del_shapes()
			self.wfft.add_shapes(shp)
			self.wfft.updateGL()
		# Single measurement circle mode
		elif self.f2danmode==1 :
			self.wfft.del_shapes()
			if self.ringrad==0: self.ringrad=1.0
			self.wfft.add_shape("ring",EMShape(("circle",0.2,1.0,0.2,r,r,self.ringrad,1.0)))
			self.wfft.add_shape("ringlbl",EMShape(("scrlabel",0.2,1.0,0.2,10,10,"r=%d pix -> 1/%1.2f 1/A (%1.4f)"%(self.ringrad,1.0/(self.ringrad*ds),self.ringrad*ds),24.0,2.0)))
			self.wfft.updateGL()
		# 2-D Crystal mode
		elif self.f2danmode==2 :
			shp={}
			for a in range(-5,6):
				for b in range(-5,6):
					shp["m%d%d"%(a,b)]=EMShape(("circle",1.0,0.0,0.0,a*self.xpos1[0]+b*self.xpos2[0]+self.fft["nx"]/2-1,a*self.xpos1[1]+b*self.xpos2[1]+self.fft["ny"]/2,3,1.0))

			self.wfft.del_shapes()
			self.wfft.add_shapes(shp)
			self.wfft.add_shape("xtllbl",EMShape(("scrlabel",1.0,0.3,0.3,10,10,"Unit Cell: %1.2f,%1.2f"%(1.0/(hypot(*self.xpos1)*ds),1.0/(hypot(*self.xpos2)*ds)),60.0,2.0)))
#			except: pass
			self.wfft.updateGL()
		else:
			self.wfft.del_shapes()
			self.wfft.updateGL()


		# Now update the plots for the correct plot mode
		if self.plotmode==0:
			try: bgsub=self.fft1d-bg1d
			except:
				print "Error computing bgsub on this image"
				return
			self.wplot.set_data((s,bgsub),"fg-bg",quiet=True,color=0)

			fit=array(ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP))		# The fit curve
			fit=fit*fit			# squared

			# auto-amplitude for b-factor adjustment
			rto,nrto=0,0
			for i in range(int(.04/ds)+1,min(int(0.15/ds),len(s)-1)):
				if bgsub[i]>0 :
					rto+=fit[i]
					nrto+=fabs(bgsub[i])
			if nrto==0 : rto=1.0
			else : rto/=nrto
			fit=[fit[i]/rto for i in range(len(s))]

#			print ctf_cmp((self.sdefocus.value,self.sbfactor.value,rto),(ctf,bgsub,int(.04/ds)+1,min(int(0.15/ds),len(s)-1),ds,self.sdefocus.value))

			self.wplot.set_data((s,fit),"fit",color=1)
			self.wplot.setAxisParms("s (1/"+ "$\AA$" +")","Intensity (a.u)")
		elif self.plotmode==1:
			self.wplot.set_data((s[1:],self.fft1d[1:]),"fg",quiet=True,color=1)
			self.wplot.set_data((s[1:],bg1d[1:]),"bg",color=0)
			self.wplot.setAxisParms("s (1/"+ "$\AA$" +")","Intensity (a.u)")
		elif self.plotmode==2:
			if self.fft1dang==None: self.recalc_real()
			bgsub=self.fft1d-bg1d
			bgsuba=[array(self.fft1dang[i])-bg1d for i in xrange(4)]
					# Write the current image parameters to the database

#			for i in xrange(4): bgsuba[i][0]=0
			self.wplot.set_data((s,bgsub),"fg",quiet=True,color=0)
			self.wplot.set_data((s[3:],bgsuba[0][3:]),"fg 0-45",quiet=True,color=2)
			self.wplot.set_data((s[3:],bgsuba[1][3:]),"fg 45-90",quiet=True,color=3)
			self.wplot.set_data((s[3:],bgsuba[2][3:]),"fg 90-135",quiet=True,color=4)
			self.wplot.set_data((s[3:],bgsuba[3][3:]),"fg 135-180",quiet=True,color=6)

			fit=array(ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP))		# The fit curve
			fit=fit*fit			# squared

			# auto-amplitude for b-factor adjustment
			rto,nrto=0,0
			for i in range(int(.04/ds)+1,min(int(0.15/ds),len(s)-1)):
				if bgsub[i]>0 :
					rto+=fit[i]
					nrto+=fabs(bgsub[i])
			if nrto==0 : rto=1.0
			else : rto/=nrto
			fit/=rto

			self.wplot.set_data((s,fit),"fit",color=1)
			self.wplot.setAxisParms("s (1/"+ "$\AA$" + ")","Intensity (a.u)")

		elif self.plotmode==3:
			if self.fft1dang==None: self.recalc_real()
			#bgsub=self.fft1d-bg1d
			#bgsuba=[array(self.fft1dang[i])-bg1d for i in xrange(4)]
			fg=self.fft1d
			fga=[array(self.fft1dang[i]) for i in xrange(4)]

			for i in xrange(4): fga[i][0]=0
			self.wplot.set_data((s,fg),"fg",quiet=True,color=0)
			self.wplot.set_data((s,fga[0]),"fg 0-45",quiet=True,color=2)
			self.wplot.set_data((s,fga[1]),"fg 45-90",quiet=True,color=3)
			self.wplot.set_data((s,fga[2]),"fg 90-135",quiet=True,color=4)
			self.wplot.set_data((s,fga[3]),"fg 135-180",quiet=True,color=6)

			#fit=array(ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP))		# The fit curve
			#fit=fit*fit			# squared

			## auto-amplitude for b-factor adjustment
			#rto,nrto=0,0
			#for i in range(int(.04/ds)+1,min(int(0.15/ds),len(s)-1)):
				#if bgsub[i]>0 :
					#rto+=fit[i]
					#nrto+=fabs(bgsub[i])
			#if nrto==0 : rto=1.0
			#else : rto/=nrto
			#fit/=rto

			#self.wplot.set_data((s,fit),"fit",color=1)
			#self.wplot.setAxisParms("s (1/"+ "$\AA$" + ")","Intensity (a.u)")
		if self.plotmode==4:
			if min(bg1d)<=0.0 : bg1d+=min(bg1d)+max(bg1d)/10000.0
			ssnr=(self.fft1d-bg1d)/bg1d
			self.wplot.set_data((s,ssnr),"SSNR",quiet=True,color=0)

			#fit=array(ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP))		# The fit curve
			#fit=fit*fit			# squared

			## auto-amplitude for b-factor adjustment
			#rto,nrto=0,0
			#for i in range(int(.04/ds)+1,min(int(0.15/ds),len(s)-1)):
				#if bgsub[i]>0 :
					#rto+=fit[i]
					#nrto+=fabs(bgsub[i])
			#if nrto==0 : rto=1.0
			#else : rto/=nrto
			#fit=[fit[i]/rto for i in range(len(s))]

##			print ctf_cmp((self.sdefocus.value,self.sbfactor.value,rto),(ctf,bgsub,int(.04/ds)+1,min(int(0.15/ds),len(s)-1),ds,self.sdefocus.value))

			#self.wplot.set_data((s,fit),"fit",color=1)
			self.wplot.setAxisParms("s (1/"+ "$\AA$" + ")","Est. SSNR")


	def timeOut(self):
		if self.busy : return

		# Redisplay before spawning thread for more interactive display
		if self.needredisp :
			try: self.redisplay()
			except: pass

		# Spawn a thread to reprocess the data
		if self.needupdate and self.procthread==None:
			self.procthread=threading.Thread(target=self.recalc_real)
			self.procthread.start()

		if self.errors:
			QtGui.QMessageBox.warning(None,"Error","The following processors encountered errors during processing of 1 or more images:"+"\n".join(self.errors))
			self.errors=None

	def doRefit(self):
		parms=self.parms[self.curset]
		apix=self.sapix.getValue()
		ds=1.0/(apix*parms[0]*parms[5])

		try:
			parms[1]=e2ctf.ctf_fit(self.fft1d,parms[1].background,parms[1].background,self.fft,self.fftbg,parms[1].voltage,parms[1].cs,parms[1].ampcont,apix,bgadj=False,autohp=True,verbose=1)
			if self.cbgadj.getValue() :
				self.bgAdj()
				parms[1]=e2ctf.ctf_fit(self.fft1d,parms[1].background,parms[1].background,self.fft,self.fftbg,parms[1].voltage,parms[1].cs,parms[1].ampcont,apix,bgadj=False,autohp=True,verbose=1)
				if self.fitastig : e2ctf.ctf_fit_stig(self.fft,self.fftbg,parms[1],verbose=1)
				self.bgAdj()
			else:
				if self.fitastig : e2ctf.ctf_fit_stig(self.fft,self.fftbg,parms[1],verbose=1)

		except:
			print "CTF Autofit Failed"
			traceback.print_exc()
			parms[1].defocus=1.0

		if self.constbfactor>0 : parms[1].bfactor=self.constbfactor
		self.sdefocus.setValue(parms[1].defocus,True)
		self.sbfactor.setValue(parms[1].bfactor,True)
		self.sampcont.setValue(parms[1].ampcont,True)

		self.update_plot()


	def unImport(self,val=None):
		print "unimport ",base_name(self.setlist.item(self.curset).text(),nodir=self.nodir)
		item=base_name(self.setlist.item(self.curset).text(),nodir=self.nodir)
		try: os.unlink("micrographs/%s.hdf"%item)
		except: print "Couldn't delete micrographs/%s.hdf"%item

	def doImport(self,val=None):
		"""Imports the currently selected image into a project"""
		print "import ",base_name(self.setlist.item(self.curset).text(),nodir=self.nodir)

		# This is just the (presumably) unique portion of the filename
		item=base_name(self.setlist.item(self.curset).text(),nodir=self.nodir)

		# create directory if necessary
		if not os.access("micrographs",os.R_OK) :
			try : os.mkdir("micrographs")
			except:
				QtGui.QMessageBox.warning(self,"Error !","Cannot create micrographs directory")
				return

		#db=db_open_dict("bdb:micrographs#%s"%item)
		self.data["ctf"]=self.parms[self.curset][1]

		if self.cinvert.getValue()!=0 : self.data.mult(-1)
		if self.cxray.getValue() : self.data.process_inplace("threshold.clampminmax.nsigma",{"nsigma":4,"tomean":1})
		self.data.write_image("micrographs/%s.hdf"%item)
		self.writeCurParm()

#		js_open_dict(info_name(item))["ctf"]=[self.parms[val][1],None,None,None,None]
# 		db_parms=db_open_dict("bdb:e2ctf.parms")
# 		db_parms[item]=[self.parms[val][1].to_string(),self.fft1d,self.parms[val][1].background,self.parms[val][4]]

	def writeCurParm(self):
		"Called to store the current parameters for this image to the frameparms database"
		parms=self.parms[self.curset]
		js=js_open_dict(info_name(self.curfilename,nodir=self.nodir))
		js.setval("ctf_frame",parms,True)
		js.setval("quality",parms[4])
# 		db_fparms=db_open_dict("bdb:e2ctf.frameparms")
# 		curtag=item_name(str(self.setlist.item(self.curset).text()))
# 		db_fparms[curtag]=self.parms[self.curset]

	def newSet(self,val):
		"called when a new data set is selected from the list"

		# Write the current image parameters to the database
		if val!=self.curset : self.writeCurParm()

		# now set the new item
		self.curset=val

		# This deals with Z stacks and multi image files
		fsp=str(self.setlist.item(val).text())


		if "," in fsp :
			fsp,n=fsp.split(",")
			self.data=EMData(fsp,int(n))	# read the image from disk
		elif ";" in fsp :
			fsp,n=fsp.split(";")
			hdr=EMData(fsp,0,True)
			self.data=EMData(fsp,0,False,Region(0,0,int(n),hdr["nx"],hdr["ny"],1))	# read the image from disk
		else :
			self.data=EMData(fsp,0)	# read the image from disk

		self.wimage.setWindowTitle("e2evalimage - " + fsp.split("/")[-1])
		self.wfft.setWindowTitle("e2evalimage - 2D FFT - "+fsp.split("/")[-1])
		self.wplot.setWindowTitle("e2evalimage - Plot - "+fsp.split("/")[-1])


		if self.defaultapix!=None : self.data["apix_x"]=self.defaultapix
		self.wimage.set_data(self.data)
		self.curfilename = str(self.setlist.item(val).text())

		ctf=self.parms[val][1]
		# if voltage is 0, we need to initialize
		if ctf.voltage==0:
			try: ctf.voltage=self.data["microscope_voltage"]
			except: pass
			if ctf.voltage==0 : ctf.voltage=self.defaultvoltage
			try: ctf.cs=self.data["microscope_cs"]
			except: pass
			if ctf.cs==0 : ctf.cs=self.defaultcs
			ctf.apix=self.data["apix_x"]
			ctf.defocus=0.0		#triggers fitting
			ctf.bfactor=200.0
			ctf.ampcont=self.defaultac

		if self.constbfactor>0 : ctf.bfactor=self.constbfactor

		self.sdefocus.setValue(ctf.defocus,True)
		self.sbfactor.setValue(ctf.bfactor,True)
		self.sapix.setValue(ctf.apix,True)
		self.sampcont.setValue(ctf.ampcont,True)
		self.svoltage.setValue(ctf.voltage,True)
		self.scs.setValue(ctf.cs,True)
		self.sboxsize.setValue(self.parms[val][0],True)
		self.squality.setValue(self.parms[val][4],True)

		#if self.guiim != None:
##			print self.data
			#self.guiim.set_data(self.data[val][4])
			#if self.guiiminit:
				#self.guiim.optimally_resize()
				#self.guiiminit = False
			#self.guiim.updateGL()
		self.recalc()

	def recalc(self):
		self.needupdate=True

	def recalc_real(self):
		"Called to recompute the power spectra, also updates plot"

		self.needupdate=False

		if self.data==None :
			self.procthread=None
			return

		self.incalc=True	# to avoid incorrect plot updates

		# To simplify expressions
		parms=self.parms[self.curset]
		apix=self.sapix.getValue()
		if len(parms)==5 : parms.append(1)		# for old projects where there was no oversampling specification
		else: parms[5]=max(1,int(parms[5]))
		ds=1.0/(apix*parms[0]*parms[5])

		# Mode where user drags the box around the parent image
		if self.calcmode==0:

			# extract the data and do an fft
			clip=self.data.get_clip(Region(parms[2][0],parms[2][1],parms[0],parms[0]))
			clip.process_inplace("normalize.edgemean")

			if parms[5]>1 :
				clip=clip.get_clip(Region(0,0,parms[0]*parms[5],parms[0]*parms[5]))		# since we aren't using phases, doesn't matter if we center it or not
			self.fft=clip.do_fft()
#			self.fft.mult(1.0/parms[0]**2)
			self.fft.mult(1.0/parms[0])

		# mode where user selects/deselcts tiled image set
		elif self.calcmode==1:
			# update the box display on the image
			nx=self.data["nx"]/parms[0]-1
			self.fft=None
			nbx=0
			for x in range(nx):
				for y in range(self.data["ny"]/parms[0]-1):
					# User deselected this one
					if int(x+y*nx) in parms[3] : continue

					# read the data and make the FFT
					clip=self.data.get_clip(Region(x*parms[0]+parms[0]/2,y*parms[0]+parms[0]/2,parms[0],parms[0]))
					clip.process_inplace("normalize.edgemean")
					if parms[5]>1 :
						clip=clip.get_clip(Region(0,0,parms[0]*parms[5],parms[0]*parms[5]))		# since we aren't using phases, doesn't matter if we center it or not
					fft=clip.do_fft()
#					fft.mult(parms[0])
					fft.ri2inten()
					if self.fft==None: self.fft=fft
					else: self.fft+=fft
					nbx+=1

			self.fft.mult(1.0/(nbx*parms[0]**2))
			self.fft.process_inplace("math.sqrt")
			self.fft["is_intensity"]=0				# These 2 steps are done so the 2-D display of the FFT looks better. Things would still work properly in 1-D without it
#			self.fft.mult(1.0/(nbx*parms[0]**2))

		self.fftbg=self.fft.process("math.nonconvex")
		self.fft1d=self.fft.calc_radial_dist(self.fft.get_ysize()/2,0.0,1.0,1)	# note that this handles the ri2inten averages properly
		if self.plotmode==2 or self.plotmode==3:
			self.fft1dang=array(self.fft.calc_radial_dist(self.fft.get_ysize()/2,0.0,1.0,4,self.sang45.getValue()*.017453292,1))	# This form generates 4 sequential power spectra representing angular ranges
			self.fft1dang=self.fft1dang.reshape((4,self.fft.get_ysize()/2))
		else:
			self.fft1dang=None

		# Compute 1-D curve and background
		bg_1d=e2ctf.low_bg_curve(self.fft1d,ds)
		parms[1].background=bg_1d
		parms[1].dsbg=ds

		self.fft1d=array(self.fft1d)

		self.needredisp=True
		self.incalc=False
		time.sleep(.2)			# help make sure update has a chance
		self.procthread=None

		self.bgAdj()
#		dbquality = self.db[os.path.basename(self.curfilename)]
#		print dbquality
# 		item=item_name(self.curfilename)
# 		db=db_open_dict("bdb:e2ctf.parms")
# 		if db[item]==None or db[item[3]]==None : db[item]=[parms[1].to_string(),self.fft1d,parms[1].background,5]
# #		print item,db[item][3]
# 		try: self.squality.setValue(int(db[item[3]]))
# 		except: self.squality.setValue(5)

	def redisplay(self):

		if self.incalc: return

		self.needredisp=False
		self.busy=True
		parms=self.parms[self.curset]
		apix=self.sapix.getValue()
		ds=1.0/(apix*parms[0]*parms[5])


		# Fitting not done yet. Need to make 2D background somehow
		if parms[1].defocus==0:
			self.doRefit()

		self.wimage.show()
		self.wfft.show()
		self.wplot.show()
		try:
			self.wimage.raise_()
			self.wfft.raise_()
			self.wplot.raise_()
		except: pass

		self.update_plot()

		# To simplify expressions

		if self.calcmode==0:
			# update the box display on the image
			self.wimage.del_shapes()
			self.wimage.add_shape("box",EMShape(("rect",.3,.9,.3,parms[2][0],parms[2][1],parms[2][0]+parms[0],parms[2][1]+parms[0],1)))
			self.wimage.updateGL()
		elif self.calcmode==1:
			# update the box display on the image
			nx=self.data["nx"]/parms[0]-1
			shp={}
			for x in range(nx):
				for y in range(self.data["ny"]/parms[0]-1):
					# User deselected this one
					if int(x+y*nx) in parms[3] : continue

					# Make a shape for this box
					shp["box%02d%02d"%(x,y)]=EMShape(("rect",.3,.9,.3,(x+.5)*parms[0],(y+.5)*parms[0],(x+1.5)*parms[0],(y+1.5)*parms[0],1))

			self.wimage.del_shapes()
			self.wimage.add_shapes(shp)
			self.wimage.updateGL()

		if self.f2dmode>0 :
			if self.f2dmode==1 : self.wfft.set_data(self.fft-self.fftbg)
			else : self.wfft.set_data(self.fftbg)

		else :
			self.wfft.set_data(self.fft)
		self.busy=False

	def newCalcMode(self,mode):
		self.calcmode=mode
		self.recalc()

	def new2DMode(self,mode):
		self.f2dmode=mode
		self.recalc()

	def new2DAnMode(self,mode):
		self.f2danmode=mode
		self.needredisp=True

	def newPlotMode(self,mode):
		self.plotmode=mode
		self.wplot.set_data(None,replace=True,quiet=True)	# clear the data so plots are properly redisplayed, but don't update the display
		self.needredisp=True

	def newBox(self):
		parms=self.parms[self.curset]
		parms[0]=self.sboxsize.value
#		parms[5]=self.soversamp.value
		parms[5]=1
		parms[3]=set()
		self.recalc()

	def newQualityFactor(self):
		parms=self.parms[self.curset]
		parms[4]=self.squality.value

	def bgAdj(self):

		if self.cbgadj.getValue() :
			parms=self.parms[self.curset]
			apix=self.sapix.getValue()
#			ds=1.0/(apix*parms[0]*parms[5])
			ds=parms[1].dsbg
			ctf=parms[1]
#			bg_1d=e2ctf.low_bg_curve(self.fft1d,ds)
			bg_1d=list(parms[1].background)

			xyd=XYData()
#			lz=int(ctf.zero(0)/ds)

			# Find the minimum value near the origin, which we'll use as a zero (though it likely should not be)
			mv=(self.fft1d[1],1)
			fz=int(ctf.zero(0)/(ds*2))
			for lz in xrange(1,fz):
				mv=min(mv,(self.fft1d[lz],lz))

#			print mv,int(ctf.zero(0)/(ds*2)),min(self.fft1d[1:int(ctf.zero(0)/(ds*2))])
			xyd.insort(mv[1],mv[0])

			# now we add all of the zero locations to our XYData object
			for i in xrange(100):
				z=int(ctf.zero(i)/ds)
				if z>=len(bg_1d)-1: break
				if self.fft1d[z-1]<self.fft1d[z] and self.fft1d[z-1]<self.fft1d[z+1]: mv=(z-1,self.fft1d[z-1])
				elif self.fft1d[z]<self.fft1d[z+1] : mv=(z,self.fft1d[z])
				else : mv=(z+1,self.fft1d[z+1])
				xyd.insort(mv[0],mv[1])

			# new background is interpolated XYData
			parms[1].background=[xyd.get_yatx_smooth(i,1) for i in xrange(len(bg_1d))]

			# if our first point (between the origin and the first 0) is too high, we readjust it once
			bs=[self.fft1d[i]-parms[1].background[i] for i in xrange(fz)]
			if min(bs)<0 :
				mv=(bs[0],self.fft1d[0],0)
				for i in xrange(1,fz): mv=min(mv,(bs[i],self.fft1d[i],i))
				xyd.set_x(0,mv[2])
				xyd.set_y(0,mv[1])
				
				parms[1].background=[xyd.get_yatx_smooth(i,1) for i in xrange(len(bg_1d))]

		self.needredisp=True


	def newCTF(self) :
		parms=self.parms[self.curset]
		parms[1].defocus=self.sdefocus.value
		parms[1].bfactor=self.sbfactor.value
		parms[1].dfdiff=self.sdfdiff.value
		parms[1].dfang=self.sdfang.value
		parms[1].apix=self.sapix.value
		parms[1].ampcont=self.sampcont.value
		parms[1].voltage=self.svoltage.value
		parms[1].cs=self.scs.value
		self.bgAdj()

		self.needredisp=True


	def imgmousedown(self,event) :
		if self.calcmode==0:
			m=self.wimage.scr_to_img((event.x(),event.y()))
			parms=self.parms[self.curset]
			parms[2]=(m[0]-parms[0]/2,m[1]-parms[0]/2)
			self.recalc()
			self.needredisp=True
		#self.guiim.add_shape("cen",["rect",.9,.9,.4,x0,y0,x0+2,y0+2,1.0])

	def imgmousedrag(self,event) :
		if self.calcmode==0:
			m=self.wimage.scr_to_img((event.x(),event.y()))
			parms=self.parms[self.curset]
			parms[2]=(m[0]-parms[0]/2,m[1]-parms[0]/2)
			self.needredisp=True
			self.recalc()

		# box deletion when shift held down
		#if event.modifiers()&Qt.ShiftModifier:
			#for i,j in enumerate(self.boxes):

	def imgmouseup(self,event) :
		m=self.wimage.scr_to_img((event.x(),event.y()))
		if self.calcmode==1:
			parms=self.parms[self.curset]
			nx=self.data["nx"]/parms[0]-1
			grid=int((m[0]-parms[0]/2)/parms[0])+int((m[1]-parms[0]/2)/parms[0])*nx
			if grid in parms[3] : parms[3].remove(grid)
			else: parms[3].add(grid)
			self.needredisp=True
			self.recalc()


	def fftmousedown(self,event,m) :
		#m=self.wfft.scr_to_img((event.x(),event.y()))

		if self.f2danmode==1:
			self.ringrad=hypot(m[0]-self.fft["nx"]/2,m[1]-self.fft["ny"]/2)
			self.needredisp=True
		elif self.f2danmode==2:
			if (event.modifiers()&Qt.ControlModifier): self.xpos2=((m[0]-self.fft["nx"]/2)/3.0,(m[1]-self.fft["ny"]/2)/3.0)
			else: self.xpos1=((m[0]-self.fft["nx"]/2)/3.0,(m[1]-self.fft["ny"]/2)/3.0)
			self.needredisp=True


		#self.guiim.add_shape("cen",["rect",.9,.9,.4,x0,y0,x0+2,y0+2,1.0])

	def fftmousedrag(self,event,m) :
		#m=self.wfft.scr_to_img((event.x(),event.y()))

		if self.f2danmode==1:
			self.ringrad=hypot(m[0]-self.fft["nx"]/2,m[1]-self.fft["ny"]/2)
			self.needredisp=True
		elif self.f2danmode==2:
			if (event.modifiers()&Qt.ControlModifier): self.xpos2=((m[0]-self.fft["nx"]/2)/3.0,(m[1]-self.fft["ny"]/2)/3.0)
			else: self.xpos1=((m[0]-self.fft["nx"]/2)/3.0,(m[1]-self.fft["ny"]/2)/3.0)
			self.needredisp=True
		# box deletion when shift held down
		#if event.modifiers()&Qt.ShiftModifier:
			#for i,j in enumerate(self.boxes):

	def fftmouseup(self,event,m) :
		"up"
		#m=self.wfft.scr_to_img((event.x(),event.y()))


	def plotmousedown(self,event) :
		"mousedown in plot"
#		m=self.guiim.scr_to_img((event.x(),event.y()))

if __name__ == "__main__":
	main()
