#!/usr/bin/env python

#
# Author: Steven Ludtke 11/09/2006 (sludtke@bcm.edu)
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA
#
#

from EMAN2 import *
from emimage import EMImage
try: from emplot2d import EMPlot2D,NewPlot2DWin
except: pass

from emimageutil import EMParentWin
import sys
from optparse import OptionParser
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from OpenGL import GL,GLU,GLUT
#from valslider import ValSlider
#from math import *
#import numpy
#from emimageutil import ImgHistogram
#from weakref import WeakKeyDictionary

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <image file> ...
	
"""
	global app,win,options

	parser = OptionParser(usage=usage,version=EMANVERSION)

#	parser.add_option("--gui",action="store_true",help="Start the GUI for interactive boxing",default=False)
#	parser.add_option("--auto","-A",type="string",action="append",help="Autobox using specified method: circle, ref, grid",default=[])
#	parser.add_option("--threshold","-T",type="float",help="Threshold for keeping particles. 0-4, 0 excludes all, 4 keeps all.",default=2.0)
#	parser.add_option("--maxbad","-M",type="int",help="Maximumum number of unassigned helices",default=2)
#	parser.add_option("--minhelix","-H",type="int",help="Minimum residues in a helix",default=6)
#	parser.add_option("--apix","-P",type="float",help="A/Pixel",default=1.0)
	parser.add_option("--classmx",type="string",help="<classmx>,<#> Show particles in one class from a classification matrix. Pass raw particle file as first argument to command.")
	parser.add_option("--classes",type="string",help="<rawptcl>,<classmx> Show particles associated class-averages")
	parser.add_option("--plot",action="store_true",default=False,help="Data file(s) should be plotted rather than displayed in 2-D")
	
	(options, args) = parser.parse_args()
	if len(args)<1 : parser.error("Input image required")

	logid=E2init(sys.argv)
#        GLUT.glutInit(sys.argv)

	app = QtGui.QApplication(sys.argv)
	win=[]

	if options.plot:
		plot(args)
	elif options.classes:
		options.classes=options.classes.split(",")
		imgs=EMData.read_images(args[0])
		display(imgs,args[0])
		
		QtCore.QObject.connect(win[0].child,QtCore.SIGNAL("mousedown"),lambda a,b:selectclass(options.classes[0],options.classes[1],a,b))
		try:
			out=file("selected.lst","w")
			out.write("#LST\n")
			out.close()
		except: pass
	elif options.classmx:
		options.classmx=options.classmx.split(",")
		clsnum=int(options.classmx[1])
		imgs=getmxim(args[0],options.classmx[0],clsnum)
		display(imgs,args[0])
	else:
		for i in args:
			a=EMData.read_images(i)
			display(a,i)
			
	sys.exit(app.exec_())

	E2end(logid)

def selectclass(rawfsp,mxfsp,event,lc):
	"""Used in 'classes' mode when the user selects a particular class-average"""
	global win
	
	clsnum=lc[0]
	
	if event.modifiers()==Qt.ShiftModifier :
		ptcls=getmxinfo(rawfsp,mxfsp,clsnum)
		out=file("selected.lst","a")
		for i in ptcls: i.write("%d\t%s\n"%(i[1],i[0]))
		out.close()
	else:	
		ptcls=getmxim(rawfsp,mxfsp,clsnum)
		if len(win)==1:
			display(ptcls)
		else:
			win[1].child.setData(ptcls)

def getmxinfo(fsp,fsp2,clsnum):
	"""reads particle references associated with a particular class.
	fsp2 is the matrix file and fsp is the raw particle file"""
	mx=EMData(fsp2,0)
	imgs=[(fsp,i) for i in range(mx.get_ysize()) if mx.get(0,i)==clsnum]
	return imgs

def getmxim(fsp,fsp2,clsnum):
	"""reads the raw particles associated with a particular class.
	fsp2 is the matrix file and fsp is the raw particle file"""
	mx=EMData(fsp2,0)
	dx=EMData(fsp2,2)
	dy=EMData(fsp2,3)
	da=EMData(fsp2,4)
	imgs=[(EMData(fsp,i),dx.get(0,i),dy.get(0,i),da.get(0,i)) for i in range(mx.get_ysize()) if mx.get(0,i)==clsnum]
	for i in imgs :
		print i 
		i[0].rotate_translate(i[3],0,0,i[1],i[2],0)
	imgs=[i[0] for i in imgs]
	return imgs

def display(img,title="EMImage"):
#	print img
	
	if len(img)==1 : img=img[0]
	w=EMImage(img)
	w.setWindowTitle(title)
	try:
		if file_exists(title):
			w.child.setImageFileName(title)
	except: pass
	w.show()
	win.append(w)

def plot(files):
	plotw=EMPlot2D()
	for f in files:
		plotw.setDataFromFile(f,f)
		
	w=EMParentWin(plotw)
	w.setWindowTitle("EMPlot2D")
	w.show()
	win.append(w)
	

# If executed as a program
if __name__ == '__main__':
	main()	
