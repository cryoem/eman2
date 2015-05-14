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

from EMAN2 import EMANVERSION, E2init, E2end, EMData, base_name, file_exists, EMArgumentParser
import EMAN2db
from emapplication import EMApp
import embrowser
from emdataitem3d import EMStructureItem3D
from emimage import EMImageWidget, EMWidgetFromFile
import os
import sys
from OpenGL import GL, GLU, GLUT
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt


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
	parser.add_argument("--singleimage",action="store_true",default=False,help="Display a stack in a single image view")
	parser.add_argument("--pdb",type=str,nargs='*',help="Specify the location of one or more PDB files you wish to inspect.")
	parser.add_argument("--plot",action="store_true",default=False,help="Data file(s) should be plotted rather than displayed in 2-D")
	parser.add_argument("--plot3",action="store_true",default=False,help="Data file(s) should be plotted rather than displayed in 3-D")
	parser.add_argument("--fullrange",action="store_true",default=False,help="A specialized flag that disables auto contrast for the display of particles stacks and 2D images only.")
	parser.add_argument("--newwidget",action="store_true",default=False,help="Use the new 3D widgetD. Highly recommended!!!!")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()

#	logid=E2init(sys.argv)

	app = EMApp()
	gapp = app
	#QtGui.QApplication(sys.argv)
	win=[]
	if options.fullrange:
		fullrangeparms = set_full_range()

	if len(args) < 1 and not options.pdb:
		dialog = embrowser.EMBrowserWidget(withmodal=False,multiselect=False)
		dialog.show()
		try: dialog.raise_()
		except: pass
		#QtCore.QObject.connect(dialog,QtCore.SIGNAL("ok"),on_browser_done)
		#QtCore.QObject.connect(dialog,QtCore.SIGNAL("cancel"),on_browser_cancel)
		dialog.show()
	
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
		
	elif options.pdb:
		load_pdb(options.pdb)
		
	else:
		for i in args:
			if not file_exists(i):
				print "%s doesn't exist" %i
				sys.exit(1)
			display_file(i,app,options.singleimage,usescenegraph=options.newwidget)

	if options.fullrange:
		revert_full_range(fullrangeparms)

	app.exec_()

#	E2end(logid)

def set_full_range():
	'''
	Turns all auto contrasting flags to False etc.
	This is just a convenience function for users who do no want to use e2preferences.
	This makes sense if the user would like auto contrasting to be on or off at the same time (i.e. at one moment they
	want it off, and at the next they want it on, on a regular basis).
	Note if auto contrast is on that it is not difficult to manipulate the contrast settings manually anyway.
	Regular users are advised to just use e2preferences.
	@return the current settings - so the calling function can call revert_full_range
	'''
	current_settings = {}
	#global HOMEDB
	#HOMEDB=EMAN2db.EMAN2DB.open_db()
	#HOMEDB.open_dict("display_preferences")
	#db = HOMEDB.display_preferences
	#auto_contrast = db.get("display_2d_auto_contrast",dfl=True)
	#db["display_2d_auto_contrast"] = False
	#current_settings["display_2d_auto_contrast"] = auto_contrast
	current_settings["display_2d_auto_contrast"] = True

	#stack_auto_contrast = db.get("display_stack_auto_contrast",dfl=True)
	#stack_np_for_auto = db.get("display_stack_np_for_auto",dfl=5)

	#db["display_stack_auto_contrast"] = False
	#current_settings["display_stack_auto_contrast"] = stack_auto_contrast
	current_settings["display_stack_auto_contrast"] = True

	#db["display_stack_np_for_auto"] = -1
	#current_settings["display_stack_np_for_auto"] = stack_np_for_auto
	current_settings["display_stack_np_for_auto"] = -1

	return current_settings

def revert_full_range(d):
	'''
	Reverts the call to set_full_range.
	@param d - that which was returned by set_full_range
	'''
	#global HOMEDB
	#HOMEDB=EMAN2db.EMAN2DB.open_db()
	#HOMEDB.open_dict("display_preferences")
	#db = HOMEDB.display_preferences

	#for key,value in d.items():
		#db[key] = d[key]

	pass

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
	from emplot2d import EMPlot2DWidget
	plotw=EMPlot2DWidget(application=app)
	for f in files:
		plotw.set_data_from_file(f)
	plotw.setWindowTitle("2D Plot")
	app.show_specific(plotw)
	return plotw

def plot_3d(files,app):
	from emplot3d import EMPlot3DWidget
	plotw=EMPlot3DWidget()
	for f in files:
		plotw.set_data_from_file(f)
	plotw.setWindowTitle("3D Plot")
	app.show_specific(plotw)
	return plotw

def load_pdb(pdbfs):
		from emscene3d import EMScene3D
		viewer = EMScene3D()
		models = [EMStructureItem3D(pdbf) for pdbf in pdbfs]
		viewer.addChildren(models)
		viewer.show()

# If executed as a program
if __name__ == '__main__':
	main()
