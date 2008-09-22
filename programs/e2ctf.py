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
	parser.add_option("--dbin",type="string",help="Box locations used when input is a whole micrograph")
	parser.add_option("--powspec",action="store_true",help="Compute the power spectrum of the input image(s)",default=False)
	parser.add_option("--bgmask",type="int",help="Compute the background power spectrum from the edge of the image, specify a mask radius in pixels which will largely mask out the particles",default=0)
	parser.add_option("--runav",type="int",help="Apply a running average to the background with +-runav pixels",default=0)
	parser.add_option("--apix",type="float",help="Angstroms per pixel for all images",default=None)
	parser.add_option("--nonorm",action="store_true",help="Suppress per image real-space normalization",default=False)
	
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
	
	ps1d=[]
	ps2d=[]
	# This is for reading an entire micrograph
	# and processing it two different ways from the box db
	if options.dbin:
		im=EMData(args[0],0)
		im.process_inplace("normalize.edgemean")
		
		# x,y,xsize,ysize,quality,changed
		boxes=[[int(j) for j in i.split()] for i in file(options.dbin,"r").readlines() if i[0]!="#"]	# this reads the whole box db file
		
		ps2d.append(powspecdb(im,boxes))
		i3,i2=powspecbg(im,boxes[0][2])
		ps2d.append(i3)
		
		names=args[:]
		names.append(names[0]+" BG")
		
		ps1d=[ps2d[0].calc_radial_dist(ps2d[0].get_ysize()/2,0.0,1.0,1),i2]
		norm=sum(ps1d[0][-4:-1])/sum(ps1d[1][-4:-1])
		for i in range(len(ps1d[1])): ps1d[1][i]*=norm
		
	# This reads already boxed images
	else :
		names=[]
		for i in args:
			im=EMData(i,0)
			ys=im.get_ysize()
			if options.bgmask:
				mask1=EMData(ys,ys,1)
				mask1.to_one()
				mask1.process_inplace("mask.sharp",{"outer_radius":options.bgmask})
				mask2=mask1.copy()*-1+1
				ps2d.append(powspec(i,not options.nonorm,mask1))
			else : ps2d.append(powspec(i,not options.nonorm))
			ps1d.append(ps2d[-1].calc_radial_dist(ps2d[-1].get_ysize()/2,0.0,1.0,1))
			ratio=pi*(options.bgmask)**2
			ratio=ratio/(ys*ys-ratio)
#			ps1d.append(bgedge1d(i,-options.bgedge,not options.nonorm))
			out=file("fg1d.txt","w")
			for a,b in enumerate(ps1d[-1]): out.write("%d\t%1.5f\n"%(a,b))
			out.close()
			names.append(i)
			if options.bgmask>0 :
				bg2d=powspec(i,not options.nonorm,mask2)
#				ps2d.append(bg2d)
				bg=bg2d.calc_radial_dist(bg2d.get_ysize()/2,0.0,1.0,1)
				bg=[y*ratio for y in bg]
				if options.runav:
					ra=options.runav
					bg=bg[:options.runav]+[sum(bg[v-ra:v+ra+1])/(ra*2+1) for v in range(ra,len(bg)-ra)]+bg[len(bg)-ra:]
				for v in range(len(bg)): ps1d[-1][v]-=bg[v]
#				ps1d.append(bg)
#				names.append(i+"(bg)")

				out=file("bg1d.txt","w")
#				for a,b in enumerate(ps1d[-1]): out.write("%d\t%1.5f\n"%(a,b))
				for a,b in enumerate(bg): out.write("%d\t%1.5f\n"%(a,b))
				out.close()

			#if options.bgedge>0 :
				#bg=bgedge2d(i,options.bgedge)
				#ps2d.append(bg)
				#ps1d.append(ps2d[-1].calc_radial_dist(ps2d[-1].get_ysize()/2,0.0,1.0,1))
				#names.append(i+"(bg)")
			
#			if options.bgedge>0 :
#				bg=bgedge1d(i,options.bgedge,not options.nonorm)
##				for ii in range(1,len(bg)): bg[ii]/=sqrt(ii)
##				bgs=[ps1d[-1][ii]/bg[ii] for ii in range(2,len(bg)/2)]
##				mins=min(bgs)
##				bg=[ii*mins for ii in bg]
#				ps1d.append(bg)
#				ps2d.append(ps2d[-1])
#				names.append(i+"(bg)")
#				out=file("bg1d.txt","w")
#				for a,b in enumerate(ps1d[-1]): out.write("%d\t%1.5f\n"%(a,b))
#				out.close()

			

	if options.gui : 
		gui=GUIctf(names,ps1d,ps2d)
		gui.sapix.setValue(options.apix)
		gui.run()

def powspecbg(image,size):
	"""This routine will 'gridbox' the entire image, and compute a power spectrum consisting
	of the minimum value in each pixel. Hopefully this will approximate the background."""
	
	avgr=Averagers.get("minmax",{"max":0})
#	avgr=Averagers.get("image")
	allav=[]	# used for development, can be removed in production
	
	norm=size*size
	n=0
	for y in range(size/2,image.get_ysize()-size*3/2,size/2):
		for x in range(size/2,image.get_xsize()-size*3/2,size/2):
			b=image.get_clip(Region(x,y,size,size))
			imf=b.do_fft()
			imf.ri2inten()
			i2=imf.calc_radial_dist(imf.get_ysize()/2,0.0,1.0,1)
			if n==0 : i2a=i2[:]
			else :
				for i,j in enumerate(i2):
					i2a[i]=min(j,i2a[i])
			avgr.add_image(imf)
			n+=1
			allav.append(i2)	# store the individual power spectra
	
	av=avgr.finish()
	av/=norm
	av.set_value_at(0,0,0.0)
	
	i2a=[i/norm for i in i2a]
	
	av.set_complex(1)
	av.set_attr("is_intensity", 1)
	
	# SVD on individual spectra
	allav2=[]
	for i in allav:
		im=EMData(len(i),1,1)
		im.to_zero()
		for j in range(1,len(i)): im.set_value_at(j,0,0,i[j])
		allav2.append(im) 
	#s=Util.svdcmp(allav2,10)
	#for i in s: i.write_image("basis.hed",-1)
	
	# PCA on individual spectra
	#msk=allav2[0].copy()
	#msk.to_one()
	#an=Analyzers.get("pca",{"mask":msk ,"nvec":20})
	#an.insert_images_list(allav2)
	#s=an.analyze()
	for i in allav2: i.write_image("basis.hed",-1)
	
	
	# this writes a 2D image containing all of the individual power spectra
	aa=EMData(len(allav[0]),len(allav))
	aa.to_zero()
	for y,a in enumerate(allav):
		for x in range(1,len(a)):
			aa.set_value_at(x,y,a[x])
	aa.write_image("psp.mrc")
	
	return av,i2a

def powspecdb(image,boxes):
	"""This routine will read the images from the specified file, and compute the average
	2-D power spectrum for the stack. Results returned as a 2-D FFT intensity/0 image"""
	
	
	for j,i in enumerate(boxes):
		b=image.get_clip(Region(i[0],i[1],i[2],i[3]))
		imf=b.do_fft()
		imf.ri2inten()
		if j==0: av=imf
		else: av+=imf
	
	av/=(float(len(boxes))*av.get_xsize()*av.get_ysize())
	av.set_value_at(0,0,0.0)

	av.set_complex(1)
	av.set_attr("is_intensity", 1)
	return av

def powspec(stackfile,edgenorm=True,mask=None):
	"""This routine will read the images from the specified file, optionally edgenormalize,
	optionally apply a mask then compute the average
	2-D power spectrum for the stack. Results returned as a 2-D FFT intensity/0 image"""
	
	n=EMUtil.get_image_count(stackfile)
	
	for i in range(n):
		im=EMData(stackfile,i)
		if edgenorm : im.process_inplace("normalize.edgemean")
		if mask : im*=mask
		imf=im.do_fft()
		imf.ri2inten()
		if i==0: av=imf
		else: av+=imf
	
	av/=(float(n)*av.get_ysize()*av.get_ysize())
	av.set_value_at(0,0,0.0)
#	av.process_inplace("xform.fourierorigin.tocenter")
	
	av.set_complex(1)
	av.set_attr("is_intensity", 1)
	return av

def bgedge1d(stackfile,width,edgenorm=True):
	"""This routine will read the images from the specified file, and compute the average
	1-D power spectrum computed using lines taken from the edge of the image. Returns the
	1-D power spectrum as a list of floats. Negative widths work from the center out."""
	
	n=EMUtil.get_image_count(stackfile)
	
	for i in range(n):
		im=EMData(stackfile,i)
		if edgenorm : im.process_inplace("normalize.edgemean")
		if i==0 : 
			xs=im.get_xsize()
		
			# values to cover a "width" border around the image
			if width>0:
				pos=range(width)+range(xs-width,xs)
			else: pos=range(-width,xs+width)
		
			# avg will store the running sum of the 1D spectra
			avg=EMData(xs+2,1,1)
			avg.to_zero()
			avg.set_complex(1)
			avg.set_ri(1)
		
		for j in pos:
			c=im.get_clip(Region(0,j,xs,1))
			cf=c.do_fft()
			cf.ri2inten()
			avg+=cf
			
			c=im.get_clip(Region(j,0,1,xs))
			c.set_size(xs,1,1)
			cf=c.do_fft()
			cf.ri2inten()
			avg+=cf
			
#	avg.ri2inten()
	avg/=len(pos)*2*n
	avg/=xs
#	avg/=(xs+2)/4			# didn't work the math out for this, but seems to be correct
	ret=[avg.get(v,0,0) for v in range(0,xs+2,2)]
	
	return ret

def bgedge2d(stackfile,width):
	"""This routine will read the images from the specified file, and compute the average
	2-D power spectrum computed using boxes taken from the edge of the image. Returns the
	1-D power spectrum as a list of floats"""
	
	n=EMUtil.get_image_count(stackfile)
	av=None
	
	for i in range(n):
		im=EMData(stackfile,i)
		
		xs=im.get_xsize()		# x size of image
		xst=int(floor(xs/ceil(xs/width)))	# step to use so we cover xs with width sized blocks
		
		# Build a list of all boxes around the edge
		boxl=[]
		for x in range(0,xs-xst/2,xst): 
			boxl.append((x,0))
			boxl.append((x,xs-xst))
		for y in range(xst,xs-3*xst/2,xst):
			boxl.append((0,y))
			boxl.append((xs-xst,y))
			
		for b in boxl:
			r=im.get_clip(Region(b[0],b[1],width,width))
			imf=r.do_fft()
			imf.ri2inten()
			if av : av+=imf
			else: av=imf
	
	av/=(n*len(boxl)*width*width)
	av.set_value_at(0,0,0.0)

	av.set_complex(1)
	av.set_attr("is_intensity", 1)
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


class GUIctf(QtGui.QWidget):
	def __init__(self,names,pow1d,pow2d):
		"""Implements the CTF fitting dialog using various EMImage and EMPlot2D widgets
		names is a list of strings with the names for each power spectrum
		pow1d is a list of the 1-D power spectra of the images
		pow2d is a list of EMData objects with the 2-D power spectra
		"""
		try:
			from emimage import EMImage,get_app
		except:
			print "Cannot import EMAN image GUI objects (emimage,etc.)"
			sys.exit(1)
#		try:
		from emplot2d import EMPlot2D,NewPlot2DWin
		#except:
			#print "Cannot import EMAN plot GUI objects (is matplotlib installed?)"
			#sys.exit(1)
		
		self.app=get_app()
		
		QtGui.QWidget.__init__(self,None)
		
		self.names=names
		self.pow1d=pow1d
		self.pow2d=pow2d
		if names and (len(names)!=len(pow1d) or len(names)!=len(pow2d)) :
			raise Exception,"Uneven number of data sets in GUIctf (%d,%d,%d)"%(len(names),len(pow1d),len(pow2d))
		
		if not names :
			self.names=[]
			self.pow1d=[]
			self.pow2d=[]
		
		try: self.guiimw=EMImage(pow2d[0])
		except: self.guiimw=EMImage()
		self.guiim=self.guiimw.child
		self.guiplotw=NewPlot2DWin()
		self.guiplot=self.guiplotw.child
		
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mousedown"),self.imgmousedown)
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mousedrag"),self.imgmousedrag)
		self.guiim.connect(self.guiim,QtCore.SIGNAL("mouseup")  ,self.imgmouseup)
		self.guiplot.connect(self.guiplot,QtCore.SIGNAL("mousedown"),self.plotmousedown)
		
		self.guiim.mmode="app"
		
		#try:
			#E2loadappwin("boxer","imagegeom",self.guiim)
			
			#if E2getappval("boxer","imcontrol") :
				#self.guiim.showInspector(True)
				#E2loadappwin("boxer","imcontrolgeom",self.guiim.inspector)
		#except:
			#pass
		

		# This object is itself a widget we need to set up
		self.hbl = QtGui.QHBoxLayout(self)
		self.hbl.setMargin(0)
		self.hbl.setSpacing(6)
		self.hbl.setObjectName("hbl")
		
		# plot list
		self.setlist=QtGui.QListWidget(self)
		self.setlist.setSizePolicy(QtGui.QSizePolicy.Preferred,QtGui.QSizePolicy.Expanding)
		self.hbl.addWidget(self.setlist)
		
		self.vbl = QtGui.QVBoxLayout()
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		self.hbl.addLayout(self.vbl)
		
		self.samp = ValSlider(self,(0,5.0),"Amp:",0)
		self.vbl.addWidget(self.samp)
		
		self.sdefocus=ValSlider(self,(0,5.0),"Defocus:",0)
		self.vbl.addWidget(self.sdefocus)
		
		self.sbfactor=ValSlider(self,(0,500),"B factor:",0)
		self.vbl.addWidget(self.sbfactor)
		
		self.sapix=ValSlider(self,(.2,10),"A/Pix:",2)
		self.vbl.addWidget(self.sapix)

		
		QtCore.QObject.connect(self.samp, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.sdefocus, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.sbfactor, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.sapix, QtCore.SIGNAL("valueChanged"), self.newCTF)
		QtCore.QObject.connect(self.setlist,QtCore.SIGNAL("currentRowChanged(int)"),self.newSet)

		self.updateData()
		
		self.guiimw.show()
		self.guiplotw.show()
		self.show()
		
	def newData(self,name,p1d,p2d):
		self.names.append(name)
		self.pow1d.append(p1d)
		self.pow2d.append(p2d)
		self.updateData()
		
	def updateData(self):
		"""This will make sure the various widgets properly show the current data sets"""
		self.setlist.clear()
		for i,j in enumerate(self.names):
			self.setlist.addItem(j)
			l=len(self.pow1d[i])
			apix=self.sapix.value
			self.guiplot.setData(j,[[ii/(l*2.0*apix) for ii in range(l)],self.pow1d[i]])


	def newSet(self,val):
		self.guiim.setData(self.pow2d[val])

	def newCTF(self) :
		df=self.sdefocus.value
		self.updateData()

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
	
	def plotmousedown(self,event) :
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

#class GUIctfPanel(QtGui.QWidget):
	#def __init__(self,target) :
		
		#QtGui.QWidget.__init__(self,None)
		#self.target=target
		
		#self.vbl = QtGui.QVBoxLayout(self)
		#self.vbl.setMargin(0)
		#self.vbl.setSpacing(6)
		#self.vbl.setObjectName("vbl")
		
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
