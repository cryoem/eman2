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

from EMAN2 import EMANVERSION,E2init,E2end,EMData,get_file_tag,file_exists
from emimage import EMImageModule,EMModuleFromFile

from emimageutil import EMParentWin
import sys
from optparse import OptionParser
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
from OpenGL import GL,GLU,GLUT
from emapplication import EMStandAloneApplication, EMProgressDialogModule
from emselector import EMBrowserModule
import os


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <image file> ...
	
	Can be run without any arguments, in which case a browsing interface pops up.
	
"""
	global app,win,options

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--classmx",type="string",help="<classmx>,<#> Show particles in one class from a classification matrix. Pass raw particle file as first argument to command.")
	parser.add_option("--classes",type="string",help="<rawptcl>,<classmx> Show particles associated class-averages")
	parser.add_option("--plot",action="store_true",default=False,help="Data file(s) should be plotted rather than displayed in 2-D")
	parser.add_option("--plot3",action="store_true",default=False,help="Data file(s) should be plotted rather than displayed in 3-D")
	
	(options, args) = parser.parse_args()
	

	logid=E2init(sys.argv)
#        GLUT.glutInit(sys.argv)

	app = EMStandAloneApplication()
	gapp = app
	#QtGui.QApplication(sys.argv)
	win=[]

	if len(args)<1 :
		dialog = EMBrowserModule()
		QtCore.QObject.connect(dialog,QtCore.SIGNAL("ok"),on_browser_done)
		QtCore.QObject.connect(dialog,QtCore.SIGNAL("cancel"),on_browser_cancel)
		app.show()
	elif options.plot:
		plot(args,app)
	elif options.plot3:
		plot_3d(args,app)
	elif options.classes:
		options.classes=options.classes.split(",")
		imgs=EMData.read_images(args[0])
		display(imgs,app,args[0])
		
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
		display(imgs,app,args[0])
	else:
		for i in args:
			if not file_exists(i):
				print "file doesn't exist:",i
				sys.exit(1)
			display_file(i,app)
	
	app.exec_()

	E2end(logid)
	
def on_browser_done(string_list):
	if len(string_list) != 0:
		for s in string_list:
			print s,
		print
def on_browser_cancel():
	pass

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
			win[1].child.set_data(ptcls)

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

def display_file(filename,app):
	w = EMModuleFromFile(filename,application=app)
	w.get_qt_widget().setWindowTitle(get_file_tag(filename))
	
	app.show_specific(w)
	try: w.optimally_resize()
	except: pass
	return w


def display(img,app,title="EMAN2 image"):
	if len(img)==1 : img=img[0]
	w=EMImageModule(data=img,old=None,app=app)
	w.get_qt_widget().setWindowTitle(title)
	try:
		if file_exists(title):
			w.set_file_name(title)
	except: pass
	app.show_specific(w)
	try: w.optimally_resize()
	except: pass
	return w

def plot(files,app):
	from emplot2d import EMPlot2DModule
	plotw=EMPlot2DModule(application=app)
	for f in files:
		plotw.set_data_from_file(f)
	plotw.get_qt_widget().setWindowTitle("2D Plot")
	app.show_specific(plotw)
	return plotw

def plot_3d(files,app):
	from emplot3d import EMPlot3DModule
	plotw=EMPlot3DModule(application=app)
	for f in files:
		plotw.set_data_from_file(f)
	plotw.get_qt_widget().setWindowTitle("3D Plot")
	app.show_specific(plotw)
	return plotw


# If executed as a program
if __name__ == '__main__':
	main()	
