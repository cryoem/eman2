#!/usr/bin/env python

#
# Author: Steven Ludtke, 09/04/2014 (sludtke@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
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

# e2ctf.py  09/04/2014 Steven Ludtke
# This is a program for determining CTF parameters and (optionally) phase flipping images

from EMAN2 import *
from EMAN2db import db_open_dict, db_close_dict, db_check_dict, db_list_dicts
from optparse import OptionParser
from OpenGL import GL,GLUT
from math import *
import os
import sys
import weakref
import traceback
from numpy import array,arange


def main():
	progname = os.path.basename(sys.argv[0])

	usage = """prog [options]
A simple CTF simulation program. Doesn't read or process data. Just does mathematical simulations.
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--apix",type=float,help="Angstroms per pixel for all images",default=1.0, guitype='floatbox', row=4, col=0, rowspan=1, colspan=1, mode="autofit['self.pm().getAPIX()']")
	parser.add_argument("--voltage",type=float,help="Microscope voltage in KV",default=300.0, guitype='floatbox', row=4, col=1, rowspan=1, colspan=1, mode="autofit['self.pm().getVoltage()']")
	parser.add_argument("--cs",type=float,help="Microscope Cs (spherical aberation)",default=4.1, guitype='floatbox', row=5, col=0, rowspan=1, colspan=1, mode="autofit['self.pm().getCS()']")
	parser.add_argument("--ac",type=float,help="Amplitude contrast (percentage, default=10)",default=10, guitype='floatbox', row=5, col=1, rowspan=1, colspan=1, mode='autofit')
	parser.add_argument("--samples",type=int,help="Number of samples in the plotted curve",default=256)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()

	from emapplication import EMApp
	app=EMApp()
	gui=GUIctfsim(app,options.apix,options.voltage,options.cs,options.ac,options.samples)
	gui.show_guis()
	app.exec_()

#		print "done execution"


try:
	from PyQt4 import QtCore, QtGui, QtOpenGL
	from PyQt4.QtCore import Qt
	from emshape import *
	from valslider import ValSlider
except:
	print "Warning: PyQt4 must be installed to use the --gui option"
	class dummy:
		pass
	class QWidget:
		"A dummy class for use when Qt not installed"
		def __init__(self,parent):
			print "Qt4 has not been loaded"
	QtGui=dummy()
	QtGui.QWidget=QWidget


class MyListWidget(QtGui.QListWidget):
	"""Exactly like a normal list widget but intercepts a few keyboard events"""

	def keyPressEvent(self,event):

		if event.key() in (Qt.Key_Up,Qt.Key_Down) :
			QtGui.QListWidget.keyPressEvent(self,event)
			return

		self.emit(QtCore.SIGNAL("keypress"),event)
#		event.key()==Qt.Key_I


class GUIctfsim(QtGui.QWidget):
	def __init__(self,application,apix=1.0,voltage=300.0,cs=4.1,ac=10.0,samples=256):
		"""CTF simulation dialog
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

		self.app = weakref.ref(application)

		self.df_voltage=voltage
		self.df_apix=apix
		self.df_cs=cs
		self.df_ac=ac
		self.df_samples=samples
		self.img=None

		QtGui.QWidget.__init__(self,None)
		self.setWindowIcon(QtGui.QIcon(get_image_directory() + "ctf.png"))

		self.data=[]
		self.curset=0
		self.plotmode=0

		self.guiim=EMImage2DWidget(application=self.app())
		self.guiiminit = True # a flag that's used to auto resize the first time the gui's set_data function is called
		self.guiplot=EMPlot2DWidget(application=self.app())
#		self.guirealim=EMImage2DWidget(application=self.app())	# This will show the original particle images

#		self.guirealim.connect(self.guirealim,QtCore.SIGNAL("keypress"),self.realimgkey)
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mousedown"),self.imgmousedown)
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mousedrag"),self.imgmousedrag)
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mouseup")  ,self.imgmouseup)
		self.guiplot.connect(self.guiplot,QtCore.SIGNAL("mousedown"),self.plotmousedown)

		self.guiim.mmode="app"

		# This object is itself a widget we need to set up
		self.hbl = QtGui.QHBoxLayout(self)
		self.hbl.setMargin(0)
		self.hbl.setSpacing(6)
		self.hbl.setObjectName("hbl")

		# plot list and plot mode combobox
		self.vbl2 = QtGui.QVBoxLayout()
		self.setlist=MyListWidget(self)
		self.setlist.setSizePolicy(QtGui.QSizePolicy.Preferred,QtGui.QSizePolicy.Expanding)
		self.vbl2.addWidget(self.setlist)

		self.splotmode=QtGui.QComboBox(self)
		self.splotmode.addItem("Amplitude")
		self.splotmode.addItem("Intensity")
		self.splotmode.addItem("Int w sum")
		self.splotmode.addItem("Amp w sum")
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

		self.imginfo=QtGui.QLabel("Info",self)
		self.vbl.addWidget(self.imginfo)

		self.sdefocus=ValSlider(self,(0,5),"Defocus:",0,90)
		self.vbl.addWidget(self.sdefocus)

		self.sbfactor=ValSlider(self,(0,1600),"B factor:",0,90)
		self.vbl.addWidget(self.sbfactor)

		self.sdfdiff=ValSlider(self,(0,1),"DF Diff:",0,90)
		self.vbl.addWidget(self.sdfdiff)

		self.sdfang=ValSlider(self,(0,180),"Df Angle:",0,90)
		self.vbl.addWidget(self.sdfang)

		self.sampcont=ValSlider(self,(0,100),"% AC",0,90)
		self.vbl.addWidget(self.sampcont)

		self.sapix=ValSlider(self,(.2,10),"A/Pix:",2,90)
		self.vbl.addWidget(self.sapix)

		self.svoltage=ValSlider(self,(0,1000),"Voltage (kV):",0,90)
		self.vbl.addWidget(self.svoltage)

		self.scs=ValSlider(self,(0,5),"Cs (mm):",0,90)
		self.vbl.addWidget(self.scs)

		self.ssamples=ValSlider(self,(32,1024),"# Samples:",0,90)
		self.ssamples.setIntonly(True)
		self.vbl.addWidget(self.ssamples)


		self.hbl_buttons = QtGui.QHBoxLayout()
		self.newbut = QtGui.QPushButton("New")
		self.hbl_buttons.addWidget(self.newbut)
		self.vbl.addLayout(self.hbl_buttons)

		self.on_new_but()

		QtCore.QObject.connect(self.sdefocus, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.sbfactor, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.sdfdiff, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.sdfang, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.sapix, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.sampcont, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.svoltage, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.scs, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.ssamples, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.setlist,QtCore.SIGNAL("currentRowChanged(int)"),self.newSet)
		QtCore.QObject.connect(self.setlist,QtCore.SIGNAL("keypress"),self.listkey)
		QtCore.QObject.connect(self.splotmode,QtCore.SIGNAL("currentIndexChanged(int)"),self.newPlotMode)

	   	QtCore.QObject.connect(self.newbut,QtCore.SIGNAL("clicked(bool)"),self.on_new_but)


		self.resize(720,380) # figured these values out by printing the width and height in resize event


		E2loadappwin("e2ctfsim","main",self)
		E2loadappwin("e2ctfsim","image",self.guiim.qt_parent)
#		E2loadappwin("e2ctf","realimage",self.guirealim.qt_parent)
		E2loadappwin("e2ctfsim","plot",self.guiplot.qt_parent)

		self.setWindowTitle("CTF")

	def listkey(self,event):

		if event.key()>=Qt.Key_0 and event.key()<=Qt.Key_9 :
			q=int(event.key())-Qt.Key_0
			self.squality.setValue(q)
		elif event.key() == Qt.Key_Left:
			self.sdefocus.setValue(self.sdefocus.getValue()-0.01)
		elif event.key() == Qt.Key_Right:
			self.sdefocus.setValue(self.sdefocus.getValue()+0.01)
		elif event.key()==Qt.Key_R :
			self.on_recall_params()


	def on_new_but(self):
		ctf=EMAN2Ctf()
		ctf.defocus=1.0
		ctf.voltage=self.df_voltage
		ctf.apix=self.df_apix
		ctf.cs=self.df_cs
		ctf.ac=self.df_ac
		ctf.samples=self.df_samples
		self.data.append((str(len(self.setlist)+1),ctf))
		self.curset=len(self.data)
		self.update_data()
		
	def show_guis(self):
		if self.guiim != None:
			self.app().show_specific(self.guiim)
		if self.guiplot != None:
			self.app().show_specific(self.guiplot)
		#if self.guirealim != None:
			#self.app().show_specific(self.guirealim)

		self.show()

	def closeEvent(self,event):
#		QtGui.QWidget.closeEvent(self,event)
#		self.app.app.closeAllWindows()
		E2saveappwin("e2ctf","main",self)

		if self.guiim != None:
			E2saveappwin("e2ctf","image",self.guiim.qt_parent)
			self.app().close_specific(self.guiim)
			self.guiim = None
		if self.guiplot != None:
			E2saveappwin("e2ctf","plot",self.guiplot.qt_parent)
			self.app().close_specific(self.guiplot)
		#if self.guirealim != None:
			#E2saveappwin("e2ctf","realimage",self.guirealim.qt_parent)
			#self.app().close_specific(self.guirealim)

		event.accept()
		self.app().close_specific(self)
		self.emit(QtCore.SIGNAL("module_closed")) # this signal is important when e2ctf is being used by a program running its own event loop

	def update_data(self):
		"""This will make sure the various widgets properly show the current data sets"""
		self.setlist.clear()
		for i,j in enumerate(self.data):
			self.setlist.addItem(j[0])
		self.setlist.setCurrentRow(self.curset)

	def update_plot(self):
		if self.guiplot == None: return # it's closed/not visible

		for d in xrange(len(self.data)):
			ctf=self.data[d][1]
			ds=1.0/(ctf.apix*2.0*ctf.samples)
			s=arange(0,ds*ctf.samples,ds)
			
			curve=array(ctf.compute_1d(len(s)*2,ds,Ctf.CtfType.CTF_AMP))
			if self.plotmode==1 or self.plotmode==2:
				curve=curve**2
			
			if self.plotmode==2 or self.plotmode==3:
				if d==0 : avg=curve[:]
				else:
					if len(curve)!=len(avg) :
						print "Number of samples must be fixed to compute an average ({})".format(d+1)
					else:
						avg+=curve
			
			self.guiplot.set_data((s,curve),self.data[d][0],d==0,True,color=d+1)

		if self.plotmode in (2,3) :
			self.guiplot.set_data((s,avg),"Sum",False,True,color=0)
			
		self.guiplot.setAxisParms("s (1/$\AA$)","CTF")

		ctf.compute_2d_complex(self.img,Ctf.CtfType.CTF_AMP,None)
		self.guiim.set_data(self.img)

	def newSet(self,val=0):
		"called when a new data set is selected from the list"
		self.curset=val

		self.sdefocus.setValue(self.data[val][1].defocus,True)
		self.sbfactor.setValue(self.data[val][1].bfactor,True)
		self.sapix.setValue(self.data[val][1].apix,True)
		self.sampcont.setValue(self.data[val][1].ampcont,True)
		self.svoltage.setValue(self.data[val][1].voltage,True)
		self.scs.setValue(self.data[val][1].cs,True)
		self.sdfdiff.setValue(self.data[val][1].dfdiff,True)
		self.sdfang.setValue(self.data[val][1].dfang,True)
		self.ssamples.setValue(self.data[val][1].samples,True)

		# make new image if necessary
		if self.img==None or self.img["ny"]!=self.data[val][1].samples :
			self.img=EMData(self.data[val][1].samples+2,self.data[val][1].samples)
			self.img.to_zero()
			self.img.set_complex(1)
			self.guiim.set_data(self.img)
#		self.imginfo.setText("%s particles     SNR = %s"%(ptcl,ssnr))

		#if self.guiim != None:
##			print self.data
			#self.guiim.set_data(self.data[val][4])
			#if self.guiiminit:
				#self.guiim.optimally_resize()
				#self.guiiminit = False
			#self.guiim.updateGL()
		#self.update_plot()

#		print "self.data[val]=",self.data[val][0].split('#')[-1]


		self.guiim.qt_parent.setWindowTitle("e2ctfsim - 2D FFT - "+self.data[val][0])
#		self.guirealim.qt_parent.setWindowTitle("e2ctf - "+self.data[val][0].split('#')[-1])
		self.guiplot.qt_parent.setWindowTitle("e2ctfsim - Plot ")

		#n=EMUtil.get_image_count(self.data[val][0])
		#if n>1:
			#self.ptcldata=EMData.read_images(self.data[val][0],range(0,min(20,n)))
			#im=sum(self.ptcldata)
			#im.mult(1.0/len(self.ptcldata))
			#self.ptcldata.insert(0,im)
			#self.guirealim.set_data(self.ptcldata)
		#else : self.guirealim.set_data([EMData()])
		self.update_plot()

	def newPlotMode(self,mode):
		self.plotmode=mode
		self.update_plot()

	def newCTF(self) :
#		print traceback.print_stack()
		self.data[self.curset][1].defocus=self.sdefocus.value
		self.data[self.curset][1].bfactor=self.sbfactor.value
		self.data[self.curset][1].dfdiff=self.sdfdiff.value
		self.data[self.curset][1].dfang=self.sdfang.value
		self.data[self.curset][1].apix=self.sapix.value
		self.data[self.curset][1].ampcont=self.sampcont.value
		self.data[self.curset][1].voltage=self.svoltage.value
		self.data[self.curset][1].cs=self.scs.value
		self.data[self.curset][1].samples=self.ssamples.value
		
		if self.img==None or self.img["ny"]!=self.ssamples.value :
			self.img=EMData(self.ssamples.value+2,self.ssamples.value)
			self.img.to_zero()
			self.img.set_complex(1)
			self.guiim.set_data(self.img)
		self.update_plot()

	def realimgkey(self,event):
		"""Keypress in the image display window"""

		if event.key()==Qt.Key_I:			# if user presses I in this window we invert the stack on disk
			fsp=self.data[self.curset][0]
			n=EMUtil.get_image_count(fsp)
			print "Inverting images in %s"%fsp
			for i in xrange(n):
				img=EMData(fsp,i)
				img.mult(-1.0)
				img.write_image(fsp,i)

			#self.ptcldata=EMData.read_images(fsp,range(0,20))
			#self.guirealim.set_data(self.ptcldata)


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
