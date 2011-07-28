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
from optparse import OptionParser
from OpenGL import GL,GLUT
from math import *
import os
import sys
import weakref

from Simplex import Simplex

debug=False
logid=None

sfcurve=None		# This will store a global structure factor curve if specified
sfcurve2=None
envelopes=[]		# simplex minimizer needs to use a global at the moment

def main():
	global debug,logid
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <single image file> ...

This program will allow you to evaluate an individual scanned micrograph or CCD frame, by looking at its
power spectrum in various ways."""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--gui",action="store_true",help="This is a GUI-only program. This option is provided for self-consistency",default=True)
	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n", type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()
	
	logid=E2init(sys.argv)

	from emapplication import EMApp
	app=EMApp()
	gui=GUIEvalImage(args)
	gui.show()
	app.execute()

	E2end(logid)


		
class GUIEvalImage(QtGui.QWidget):
	def __init__(self,images):
		"""Implements the CTF fitting dialog using various EMImage and EMPlot2D widgets
		'data' is a list of (filename,ctf,im_1d,bg_1d,im_2d,bg_2d)
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
		
		self.data=data
		self.curset=0
		self.plotmode=0
		
		self.wimage=EMImage2DWidget()
		
		self.wfft=EMImage2DWidget()
		self.guiiminit = True # a flag that's used to auto resize the first time the gui's set_data function is called
		self.wplot=EMPlot2DWidget()
		
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
		self.hbl = QtGui.QHBoxLayout(self)
		self.hbl.setMargin(0)
		self.hbl.setSpacing(6)
		self.hbl.setObjectName("hbl")
		
		# plot list and plot mode combobox
		self.vbl2 = QtGui.QVBoxLayout()
		self.setlist=QtGui.QListWidget(self)
		self.setlist.setSizePolicy(QtGui.QSizePolicy.Preferred,QtGui.QSizePolicy.Expanding)
		for i in images:
			self.setlist.addItem(i)
		self.vbl2.addWidget(self.setlist)
		
		self.splotmode=QtGui.QComboBox(self)
		self.splotmode.addItem("Live Box")
		self.splotmode.addItem("Tiled Boxes")
		self.vbl2.addWidget(self.splotmode)
		self.hbl.addLayout(self.vbl2)
		
		# ValSliders for CTF parameters
		self.vbl = QtGui.QVBoxLayout()
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		self.hbl.addLayout(self.vbl)
		
		#self.samp = ValSlider(self,(0,5.0),"Amp:",0)
		#self.vbl.addWidget(self.samp)
				
		self.sdefocus=ValSlider(self,(0,5),"Defocus:",1.0,90)
		self.vbl.addWidget(self.sdefocus)
		
		self.sbfactor=ValSlider(self,(0,1600),"B factor:",200.0,90)
		self.vbl.addWidget(self.sbfactor)
		
		self.sampcont=ValSlider(self,(0,100),"% AC",10.0,90)
		self.vbl.addWidget(self.sampcont)
		
#		self.sapix=ValSlider(self,(.2,10),"A/Pix:",2,90)
#		self.vbl.addWidget(self.sapix)

		self.hbl2 = QtGui.QHBoxLayout()

		self.sboxsize=ValBox(self,(0,500),"Box Size:",256,90)
		self.sboxsize.intonly=True
		self.hbl2.addWidget(self.sboxsize)

		self.sapix=ValBox(self,(0,500),"A/pix:",1.0,90)
		self.hbl2.addWidget(self.sapix)

		self.svoltage=ValBox(self,(0,500),"Voltage (kV):",200,90)
		self.hbl2.addWidget(self.svoltage)
		
		self.scs=ValBox(self,(0,5),"Cs (mm):",4.1,90)
		self.hbl2.addWidget(self.scs)		
		
		self.vbl.addLayout(self.hbl2)
		
		#self.hbl_buttons = QtGui.QHBoxLayout()
		#self.saveparms = QtGui.QPushButton("Save parms")
		#self.recallparms = QtGui.QPushButton("Recall")
		#self.refit = QtGui.QPushButton("Refit")
		#self.output = QtGui.QPushButton("Output")
		#self.hbl_buttons.addWidget(self.refit)
		#self.hbl_buttons.addWidget(self.saveparms)
		#self.hbl_buttons.addWidget(self.recallparms)
		#self.hbl_buttons2 = QtGui.QHBoxLayout()
		#self.hbl_buttons2.addWidget(self.output)
		#self.vbl.addLayout(self.hbl_buttons)
		#self.vbl.addLayout(self.hbl_buttons2)
		
		QtCore.QObject.connect(self.sdefocus, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.sbfactor, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.sapix, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.sampcont, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.svoltage, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.scs, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.boxsize, QtCore.SIGNAL("valueChanged"), self.newBox)
		QtCore.QObject.connect(self.setlist,QtCore.SIGNAL("currentRowChanged(int)"),self.newSet)
		QtCore.QObject.connect(self.splotmode,QtCore.SIGNAL("currentIndexChanged(int)"),self.newPlotMode)

	   	QtCore.QObject.connect(self.saveparms,QtCore.SIGNAL("clicked(bool)"),self.on_save_params)
		QtCore.QObject.connect(self.recallparms,QtCore.SIGNAL("clicked(bool)"),self.on_recall_params)
		QtCore.QObject.connect(self.refit,QtCore.SIGNAL("clicked(bool)"),self.on_refit)
		QtCore.QObject.connect(self.output,QtCore.SIGNAL("clicked(bool)"),self.on_output)
		
		self.update_data()
		
		self.resize(720,380) # figured these values out by printing the width and height in resize event
		
		
		self.setWindowTitle("CTF")

		
	def closeEvent(self,event):
#		QtGui.QWidget.closeEvent(self,event)
#		self.app.app.closeAllWindows()
		if self.guiim != None:
			self.app().close_specific(self.guiim)
			self.guiim = None 
		if self.guiplot != None:
			self.app().close_specific(self.guiplot)
		event.accept()
		self.app().close_specific(self)
		self.emit(QtCore.SIGNAL("module_closed")) # this signal is important when e2ctf is being used by a program running its own event loop

	def update_plot(self):
		if self.guiplot == None: return # it's closed/not visible
		val=self.curset
		ctf=self.data[val][1]
		ds=self.data[val][1].dsbg
		r=len(ctf.background)
		s=[ds*i for i in range(r)]
		
		# This updates the image circles
		fit=ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP)
		shp={}
		nz=0
		for i in range(1,len(fit)):
			if fit[i-1]*fit[i]<=0.0: 
				nz+=1
				shp["z%d"%i]=EMShape(("circle",0.0,0.0,1.0/nz,r,r,i,1.0))
#				if nz==1: print ("circle",0.0,0.0,1.0/nz,r,r,i,1.0)
		
		self.guiim.del_shapes()
		self.guiim.add_shapes(shp)
		self.guiim.updateGL()
		
		if self.plotmode==1:
			self.guiplot.set_data((s,self.data[val][2]),"fg",True,True,color=1)
			self.guiplot.set_data((s,self.data[val][7]),"bg(concave)",quiet=True,color=0,linetype=2)
			self.guiplot.set_data((s,self.data[val][3]),"bg",color=0)
			self.guiplot.setAxisParms("s (1/"+ u"\u212B" +")","Intensity (a.u)")
		elif self.plotmode==0: 
			bgsub=[self.data[val][2][i]-self.data[val][3][i] for i in range(len(self.data[val][2]))]
			self.guiplot.set_data((s,bgsub),"fg-bg",True,True,color=0)
			
			lowbg=self.data[val][7]
			lowbgsub=[self.data[val][2][i]-lowbg[i] for i in range(len(self.data[val][2]))]
			self.guiplot.set_data((s,lowbgsub),"fg-lowbg",False,True,color=0,linetype=2)

			
			fit=ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP)		# The fit curve
			fit=[sfact2(s[i])*fit[i]**2 for i in range(len(s))]		# squared * structure factor

			# auto-amplitude for b-factor adjustment
			rto,nrto=0,0
			for i in range(int(.04/ds)+1,min(int(0.15/ds),len(s)-1)): 
				if bgsub[i]>0 : 
					#rto+=fit[i]**2/fabs(bgsub[i])
					#nrto+=fit[i]
					#rto+=fit[i]**2
					#nrto+=bgsub[i]**2
					rto+=fit[i]
					nrto+=fabs(bgsub[i])
			if nrto==0 : rto=1.0
			else : rto/=nrto
			fit=[fit[i]/rto for i in range(len(s))]
			
#			print ctf_cmp((self.sdefocus.value,self.sbfactor.value,rto),(ctf,bgsub,int(.04/ds)+1,min(int(0.15/ds),len(s)-1),ds,self.sdefocus.value))
			
			self.guiplot.set_data((s,fit),"fit",color=1)
			self.guiplot.setAxisParms("s (1/"+ u"\u212B" + ")","Intensity (a.u)")


	def newSet(self,val):
		"called when a new data set is selected from the list"
		self.curset=val

		self.sdefocus.setValue(self.data[val][1].defocus,True)
		self.sbfactor.setValue(self.data[val][1].bfactor,True)
#		self.sapix.setValue(self.data[val][1].apix)
		self.sampcont.setValue(self.data[val][1].ampcont,True)
		self.svoltage.setValue(self.data[val][1].voltage,True)
		self.scs.setValue(self.data[val][1].cs,True)
		self.sdfdiff.setValue(self.data[val][1].dfdiff,True)
		self.sdfang.setValue(self.data[val][1].dfang,True)
		self.squality.setValue(self.data[val][6],True)

		try : ptcl=str(self.data[val][4]["ptcl_repr"])
		except: ptcl="?"
		try: 
			ctf=self.data[val][1]
			ssnr="%1.4f"%(sum(ctf.snr)/float(len(ctf.snr)))
		except:
			ssnr="?"
		self.imginfo.setText("%s particles     SNR = %s"%(ptcl,ssnr))
		
		if self.guiim != None: 
#			print self.data
			self.guiim.set_data(self.data[val][4])
			if self.guiiminit:
				self.guiim.optimally_resize()
				self.guiiminit = False
			self.guiim.updateGL()
		self.update_plot()

	def newPlotMode(self,mode):
		self.plotmode=mode
		self.update_plot()

	def newCTF(self) :
		self.data[self.curset][1].defocus=self.sdefocus.value
		self.data[self.curset][1].bfactor=self.sbfactor.value
		self.data[self.curset][1].dfdiff=self.sdfdiff.value
		self.data[self.curset][1].dfang=self.sdfang.value
#		self.data[self.curset][1].apix=self.sapix.value
		self.data[self.curset][1].ampcont=self.sampcont.value
		self.data[self.curset][1].voltage=self.svoltage.value
		self.data[self.curset][1].cs=self.scs.value
		self.update_plot()
		

	def imgmousedown(self,event) :
		m=self.guiim.scr_to_img((event.x(),event.y()))
		#self.guiim.add_shape("cen",["rect",.9,.9,.4,x0,y0,x0+2,y0+2,1.0])
		
	def imgmousedrag(self,event) :
		m=self.guiim.scr_to_img((event.x(),event.y()))
		
		# box deletion when shift held down
		#if event.modifiers()&Qt.ShiftModifier:
			#for i,j in enumerate(self.boxes):
		
	def imgmouseup(self,event) :
		m=self.guiim.scr_to_img((event.x(),event.y()))

	def fftmousedown(self,event) :
		m=self.guiim.scr_to_img((event.x(),event.y()))
		#self.guiim.add_shape("cen",["rect",.9,.9,.4,x0,y0,x0+2,y0+2,1.0])
		
	def fftmousedrag(self,event) :
		m=self.guiim.scr_to_img((event.x(),event.y()))
		
		# box deletion when shift held down
		#if event.modifiers()&Qt.ShiftModifier:
			#for i,j in enumerate(self.boxes):
		
	def fftmouseup(self,event) :
		m=self.guiim.scr_to_img((event.x(),event.y()))


	def plotmousedown(self,event) :
		m=self.guiim.scr_to_img((event.x(),event.y()))
	
	def run(self):
		"""If you make your own application outside of this object, you are free to use
		your own local app.exec_(). This is a convenience for ctf-only programs."""
		self.app.exec_()
		
#		E2saveappwin("boxer","imagegeom",self.guiim)
#		try:
#			E2setappval("boxer","imcontrol",self.guiim.inspector.isVisible())
#			if self.guiim.inspector.isVisible() : E2saveappwin("boxer","imcontrolgeom",self.guiim.inspector)
#		except : E2setappval("boxer","imcontrol",False)
		
		return


if __name__ == "__main__":
	main()
