#!/usr/bin/env python

#
# Author: Steven Ludtke, 06/06/2011
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
from emimagemx import EMImageMXWidget

import sys
from PyQt4 import QtCore, QtGui, QtOpenGL
from PyQt4.QtCore import Qt
#from OpenGL import GL,GLU,GLUT
from emapplication import EMApp
import os
from EMAN2db import *
from valslider import *
import traceback

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [classfile]

	THIS PROGRAM IS NOT YET COMPLETE

	This program allows you to manually mark bad particles via a graphical interface.

"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_pos_argument(name="particles",help="List the file to process with e2ctf here.", default="", guitype='filebox', browser="EMCTFParticlesTable(withmodal=True,multiselect=True)",  filecheck=False, row=0, col=0,rowspan=1, colspan=2, mode='autofit,tuning,genoutp,gensf')
	parser.add_argument("--allparticles",action="store_true",help="Will process all particle stacks stored in the particles subdirectory (no list of files required)",default=False, guitype='boolbox',row=1, col=0, mode='autofit,tuning,genoutp,gensf')
	parser.add_argument("--minptcl",type=int,help="Files with fewer than the specified number of particles will be skipped",default=0,guitype='intbox', row=2, col=0, mode='autofit,tuning,genoutp,gensf')
	parser.add_argument("--minqual",type=int,help="Files with a quality value lower than specified will be skipped",default=0,guitype='intbox', row=2, col=1, mode='autofit,tuning,genoutp,gensf')
	parser.add_argument("--gui",action="store_true",help="Start the GUI for interactive use (default=True)",default=True)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()

	if options.allparticles:
		args=["particles/"+i for i in os.listdir("particles") if "__" not in i and i[0]!="." and ".hed" not in i ]
		args.sort()
		if options.verbose : print "%d particle stacks identified"%len(args)

	# remove any files that don't have enough particles from the list
	if options.minptcl>0 :
		args=[i for i in args if imcount(i)>=options.minptcl]
		if options.verbose: print "{} stacks after minptcl filter".format(len(args))


	# remove files with quality too low
	if options.minqual>0 :
		outargs=[]
		for i in args:
			try:
				if js_open_dict(info_name(i))["quality"]>=options.minqual : outargs.append(i)
			except:
#				traceback.print_exc()
				print "Unknown quality for {}, including it".format(info_name(i))
				outargs.append(i)

		args=outargs

		if options.verbose: print "{} stacks after quality filter".format(len(args))
	#logid=E2init(sys.argv, options.ppid)

	app = EMApp()
	control=EMMarkPtclTool(args,verbose=options.verbose)
	control.show()
	app.execute()

#	E2end(logid)

class EMMarkPtclTool(QtGui.QMainWindow):
	"""This is a tool for marking bad particles"""

	def __init__(self,extrafiles=None,verbose=0):
		QtGui.QMainWindow.__init__(self)

		app=QtGui.qApp
		self.setWindowTitle("e2markbadparticles")

		# Menu Bar
		self.mfile=self.menuBar().addMenu("File")
#		self.mfile_save_processed=self.mfile.addAction("Save processed data")
		self.mfile_quit=self.mfile.addAction("Quit")

		self.wtabs=QtGui.QTabWidget()
		self.setCentralWidget(self.wtabs)

		self.wclasstab=EMClassPtclTool(extrafiles)
		self.wtabs.addTab(self.wclasstab,"Classes")

		self.vbl2 = QtGui.QVBoxLayout()
		self.setlist=MyListWidget(self)
		self.setlist.setSizePolicy(QtGui.QSizePolicy.Preferred,QtGui.QSizePolicy.Expanding)
		self.vbl2.addWidget(self.setlist)
		

		# file menu
		QtCore.QObject.connect(self.mfile_quit,QtCore.SIGNAL("triggered(bool)")  ,self.menu_file_quit)
		QtCore.QObject.connect(self.setlist,QtCore.SIGNAL("currentRowChanged(int)"),self.newSet)
		QtCore.QObject.connect(self.setlist,QtCore.SIGNAL("keypress"),self.listkey)

	def menu_file_quit(self):
		self.close()

	def closeEvent(self,event):
		self.wclasstab.close()
		QtGui.QWidget.closeEvent(self, event)

if __name__ == "__main__":
	main()
