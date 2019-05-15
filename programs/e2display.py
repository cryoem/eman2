#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

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

from builtins import range
from EMAN2 import EMANVERSION, E2init, E2end, EMData, base_name, file_exists, EMArgumentParser
import EMAN2db
from eman2_gui.emapplication import EMApp
from eman2_gui import embrowser
from eman2_gui.emimage import EMImageWidget, EMWidgetFromFile
from eman2_gui.emscene3d import EMScene3D
import os
import sys

import OpenGL
OpenGL.ERROR_CHECKING = False
from OpenGL import GL, GLU, GLUT
from PyQt5.QtCore import Qt

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <image file> ...

	This program can be used to visualize most files used in EMAN2. Running it without arguments
	will open a browser window with more flexible functionality than the command-line.
	
	"""
	global app,win,options

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--classmx",type=str,help="<classmx>,<#> Show particles in one class from a classification matrix. Pass raw particle file as first argument to command.")
	parser.add_argument("--classes",type=str,help="<rawptcl>,<classmx> Show particles associated class-averages")
	parser.add_argument("--pdb",type=str,help="<pdb file> Show PDB structure.")
	parser.add_argument("--singleimage",action="store_true",default=False,help="Display a stack in a single image view")
	parser.add_argument("--plot",action="store_true",default=False,help="Data file(s) should be plotted rather than displayed in 2-D")
	parser.add_argument("--hist",action="store_true",default=False,help="Data file(s) should be plotted as a histogram rather than displayed in 2-D.")
	parser.add_argument("--plot3d",action="store_true",default=False,help="Data file(s) should be plotted rather than displayed in 3-D")
	parser.add_argument("--fullrange",action="store_true",default=False,help="A specialized flag that disables auto contrast for the display of particles stacks and 2D images only.")
	parser.add_argument("--newwidget",action="store_true",default=False,help="Use the new 3D widgetD. Highly recommended!!!!")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higher number means higher level of verboseness")

	(options, args) = parser.parse_args()

#	logid=E2init(sys.argv)

	app = EMApp()
	#gapp = app
	#QtWidgets.QApplication(sys.argv)
	win=[]
	if options.fullrange:
		print("""The --fullrange option has been removed, and replaced with an option in user preferences.
To set your preferences for full-range 2-D display, please run:
e2procjson.py --setoption display2d.autocontrast:true
""")
		sys.exit(0)
	
	if len(args) < 1:
		global dialog
		file_list = []
		dialog = embrowser.EMBrowserWidget(withmodal=False,multiselect=False)
		dialog.show()
		try: dialog.raise_()
# 			QtCore.QObject.connect(dialog,QtCore.SIGNAL("ok"),on_browser_done)
# 			QtCore.QObject.connect(dialog,QtCore.SIGNAL("cancel"),on_browser_cancel)
		except: pass
	
	elif options.pdb:
		load_pdb(args,app)
	
	elif options.plot:
		plot(args,app)
		
	elif options.hist:
		hist(args,app)
	
	elif options.plot3d:
		plot_3d(args,app)
	
	elif options.classes:
		options.classes=options.classes.split(",")
		imgs=EMData.read_images(args[0])
		display(imgs,app,args[0])

		win[0].child.mousedown.connect(lambda a,b:selectclass(options.classes[0],options.classes[1],a,b))
		try:
			out=open("selected.lst","w")
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
				print("%s doesn't exist" %i)
				sys.exit(1)
			display_file(i,app,options.singleimage,usescenegraph=options.newwidget)

	app.exec_()

#	E2end(logid)

def on_browser_done():
	pass

def on_browser_cancel():
	pass

def selectclass(rawfsp,mxfsp,event,lc):
	"""Used in 'classes' mode when the user selects a particular class-average"""
	global win

	clsnum=lc[0]

	if event.modifiers()==Qt.ShiftModifier :
		ptcls=getmxinfo(rawfsp,mxfsp,clsnum)
		out=open("selected.lst","a")
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
		print(i)
		i[0].rotate_translate(i[3],0,0,i[1],i[2],0)
	imgs=[i[0] for i in imgs]
	return imgs

def display_file(filename,app,force_2d=False,usescenegraph=False):
	w = EMWidgetFromFile(filename,application=app,force_2d=force_2d)
	w.setWindowTitle(base_name(filename))
	app.show_specific(w)
	try: w.optimally_resize()
	except: pass
	try: w.raise_()
	except: pass
	return w

def display(img,app,title="EMAN2 image"):
	if len(img)==1 : img=img[0]
	w=EMImageWidget(data=img,old=None,app=app)
	w.setWindowTitle(title)
	try:
		if file_exists(title):
			w.set_file_name(title)
	except: pass
	app.show_specific(w)
	try: w.optimally_resize()
	except: pass
	return w

def plot(files,app):
	from eman2_gui.emplot2d import EMPlot2DWidget
	plotw=EMPlot2DWidget(application=app)
	for f in files:
		plotw.set_data_from_file(f,quiet=True)
	plotw.setWindowTitle("2D Plot")
	app.show_specific(plotw)
	return plotw

def hist(files,app):
	from eman2_gui.emhist import EMHistogramWidget
	histw=EMHistogramWidget(application=app)
	for f in files:
		histw.set_data_from_file(f,quiet=True)
	histw.setWindowTitle(f)
	app.show_specific(histw)
	return histw

def plot_3d(files,app):
	from eman2_gui.emplot3d import EMPlot3DWidgetNew
	plotw=EMPlot3DWidgetNew(application=app)
	for f in files:
		plotw.set_data_from_file(f)
	plotw.setWindowTitle("3D Plot")
	app.show_specific(plotw)
	return plotw

def load_pdb(files,app):
	from eman2_gui.empdbitem3D import EMPDBItem3D, EMBallStickModel
	scene=EMScene3D()
	title = []
	for f in files:
		pdb_model = EMPDBItem3D(f)
		scene.insertNewNode(f.split("/")[-1],pdb_model)
		modeltype = EMBallStickModel(f) 
		scene.insertNewNode(modeltype.representation, modeltype, parentnode = pdb_model)
		scene.addChild(EMPDBItem3D(f))
		title.append(pdb_model.getName())
	scene.setWindowTitle(" ".join(title))
	scene.show()
	return scene

# If executed as a program
if __name__ == '__main__':
	main()
