#!/usr/bin/env python

#
# Author: Steven Ludtke, 10/17/2007 (sludtke@bcm.edu)
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

# e2boxer.py  10/17/2007  Steven Ludtke
# This is a program for performing various CTF related operations on
# images and sets of images

from EMAN2 import *
from optparse import OptionParser
from math import *
import time
import os
import sys

pl=()

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <input stack/image> ...
	
Various CTF-related operations on images."""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--gui",action="store_true",help="Start the GUI for interactive fitting",default=False)
	parser.add_option("--powspec",action="store_true",help="Compute the power spectrum of the input image(s)",default=False)
	
	#parser.add_option("--boxsize","-B",type="int",help="Box size in pixels",default=-1)
	#parser.add_option("--dbin","-D",type="string",help="Filename to read an existing box database from",default=None)
	#parser.add_option("--auto","-A",type="string",action="append",help="Autobox using specified method: circle, ref, grid, pspec",default=[])
	#parser.add_option("--threshold","-T",type="float",help="(auto:ref) Threshold for keeping particles. 0-4, 0 excludes all, 4 keeps all.",default=2.0)
	#parser.add_option("--refptcl","-R",type="string",help="(auto:ref) A stack of reference images. Must have the same scale as the image being boxed.",default=None)
	#parser.add_option("--nretest",type="int",help="(auto:ref) Number of reference images (starting with the first) to use in the final test for particle quality.",default=-1)
	#parser.add_option("--retestlist",type="string",help="(auto:ref) Comma separated list of image numbers for retest cycle",default="")
	#parser.add_option("--ptclsize","-P",type="int",help="(auto:circle) Approximate size (diameter) of the particle in pixels. Not required if reference particles are provided.",default=0)
	#parser.add_option("--overlap",type="int",help="(auto:grid) number of pixels of overlap between boxes. May be negative.")
	#parser.add_option("--farfocus",type="string",help="filename or 'next', name of an aligned far from focus image for preliminary boxing",default=None)
	#parser.add_option("--dbout",type="string",help="filename to write EMAN1 style box database file to",default=None)
	#parser.add_option("--ptclout",type="string",help="filename to write boxed out particles to",default=None)
	#parser.add_option("--savealiref",action="store_true",help="Stores intermediate aligned particle images in boxali.hdf. Mainly for debugging.",default=False)
	
	(options, args) = parser.parse_args()
	if len(args)<1 : parser.error("Input image required")

	logid=E2init(sys.argv)
	
	ps2d=[]
	for i in args:
		ps2d.append(powspec(i))
	
	ps1d=[i.calc_rad_dist(i.get_ysize()/2,0.0,1.0,1) for i in ps2d]
	
	
def powspec(stackfile):
	"""This routine will read the images from the specified file, and compute the average
	2-D power spectrum for the stack. Results returned as a 2-D FFT intensity/0 image"""
	
	n=EMUtil.get_image_count(stackfile)
	
	for i in range(n):
		im=EMImage(stackfile,i)
		imf=im.do_fft()
		imf.ri2inten()
		if i==0: av=imf
		else: av+=imf
	
	av/=(float(n)*av.get_xsize()*av.get_ysize())
	
	return av


try:
	from PyQt4 import QtCore, QtGui, QtOpenGL
	from PyQt4.QtCore import Qt
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


class GUIctf:
	def __init__(self,imagefsp,boxes,thr,boxsize=-1):
		"""Implements the CTF fitting dialog using various EMImage and EMPlot2D widgets"""
		try:
			from emimage import EMImage,get_app
		except:
			print "Cannot import EMAN image GUI objects (emimage,etc.)"
			sys.exit(1)
		try:
			from emplot2d import EMPlot2D
		except:
			print "Cannot import EMAN plot GUI objects (is matplotlib installed?)"
			sys.exit(1)
		
		self.app=get_app()
		
		self.guiim=EMImage()
		self.guiplot=EMPlot2D()
		
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mousedown"),self.imgmousedown)
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mousedrag"),self.imgmousedrag)
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mouseup")  ,self.imgmouseup)
		self.guimx.connect(self.guiplot,QtCore.SIGNAL("mousedown"),self.plotmousedown)
		
		self.guiim.mmode="app"
		self.guictl=GUIctfPanel(self)
		
		#try:
			#E2loadappwin("boxer","imagegeom",self.guiim)
			
			#if E2getappval("boxer","imcontrol") :
				#self.guiim.showInspector(True)
				#E2loadappwin("boxer","imcontrolgeom",self.guiim.inspector)
		#except:
			#pass
		
		self.guiim.show()
		self.guimx.show()
		self.guictl.show()
		
		self.boxupdate()

	def imgmousedown(self,event) :
		m=self.guiim.scrtoimg((event.x(),event.y()))
		#self.guiim.addShape("cen",["rect",.9,.9,.4,x0,y0,x0+2,y0+2,1.0])
		
	def imgmousedrag(self,event) :
		m=self.guiim.scrtoimg((event.x(),event.y()))
		
		# box deletion when shift held down
		#if event.modifiers()&Qt.ShiftModifier:
			#for i,j in enumerate(self.boxes):
		
	def imgmouseup(self,event) :
		m=self.guiim.scrtoimg((event.x(),event.y()))
	
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

class GUIctfPanel(QtGui.QWidget):
	def __init__(self,target) :
		
		QtGui.QWidget.__init__(self,None)
		self.target=target
		
		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		
		#self.info = QtGui.QLabel("%d Boxes"%len(target.boxes),self)
		#self.vbl.addWidget(self.info)

		#self.thr = ValSlider(self,(0.0,3.0),"Threshold:")
		#self.thr.setValue(target.threshold)
		#self.vbl.addWidget(self.thr)

		#self.hbl1=QtGui.QHBoxLayout()
		#self.hbl1.setMargin(0)
		#self.hbl1.setSpacing(2)
		#self.vbl.addLayout(self.hbl1)
		
		#self.lblbs=QtGui.QLabel("Box Size:",self)
		#self.hbl1.addWidget(self.lblbs)
		
		#self.bs = QtGui.QLineEdit(str(target.boxsize),self)
		#self.hbl1.addWidget(self.bs)

		#self.hbl2=QtGui.QHBoxLayout()
		#self.hbl2.setMargin(0)
		#self.hbl2.setSpacing(2)
		#self.vbl.addLayout(self.hbl2)

		#self.done=QtGui.QPushButton("Done")
		#self.vbl.addWidget(self.done)
		
		#self.connect(self.bs,QtCore.SIGNAL("editingFinished()"),self.newBoxSize)
		#self.connect(self.thr,QtCore.SIGNAL("valueChanged"),self.newThresh)
		#self.connect(self.done,QtCore.SIGNAL("clicked(bool)"),self.target.app.quit)
##		self.target.connect(self.target,QtCore.SIGNAL("nboxes"),self.nboxesChanged)
		


if __name__ == "__main__":
	main()
