#!/usr/bin/env python

#
# Author: Steven Ludtke (sludtke@bcm.edu)
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA
#
#

import PyQt4
from PyQt4 import QtCore, QtGui,Qt
from emapplication import EMApp
from EMAN2 import *


class EMBrowserWidget(QtGui.QWidget):
	"""This widget is a file browser for EMAN2. In addition to being a regular file browser, it supports:
	- getting information about recognized data types
	- embedding BDB: databases into the observed filesystem
	- remote database access (EMEN2)
	"""
	
	def __init__(self,parent=None,withmodal=False):
		QtGui.QWidget.__init__(self,parent)
		
		self.vbl = QtGui.QVBoxLayout(self)
		
		# Top Toolbar area
		self.wtools=QtGui.QWidget()
		self.vbl.addWidget(self.wtools)
		
		# Central verticalregion has bookmarks and tree
		self.hbl1 = QtGui.QHBoxLayout()
		
		# Bookmarks implemented with a toolbar
		self.wbookmarks = QtGui.QToolBar()
		self.wbookmarks.setOrientation(2)
		self.wbookmarks.addAction("EMEN2",self.BookmarkPress)
		self.wbookmarks.addSeparator()
		self.wbookmarks.addAction("Root",self.BookmarkPress)
		self.wbookmarks.addAction("Current",self.BookmarkPress)
		self.wbookmarks.addAction("Home",self.BookmarkPress)
		self.hbl1.addWidget(self.wbookmarks)
		
		self.wtree = QtGui.QTreeView()
		self.hbl1.addWidget(self.wtree)
		
		self.vbl.addLayout(self.hbl1)
		
		
		self.hbl2 = QtGui.QHBoxLayout()

		self.wbutshow = QtGui.QPushButton("Show")
		self.hbl2.addWidget(self.wbutshow)
		
		self.wbutadd = QtGui.QPushButton("Add")
		self.hbl2.addWidget(self.wbutadd)
		
		self.wbutnew = QtGui.QPushButton("New")
		self.hbl2.addWidget(self.wbutnew)

		# buttons for selector use
		self.wspace1=QtGui.QSpacerItem(100,10,QtGui.QSizePolicy.MinimumExpanding)
		self.hbl2.addSpacerItem(self.wspace1)
		
		self.wbutcancel=QtGui.QPushButton("Cancel")
		self.hbl2.addWidget(self.wbutcancel)
		
		self.wbutok=QtGui.QPushButton("OK")
		self.hbl2.addWidget(self.wbutok)
		
		self.vbl.addLayout(self.hbl2)
		
	def BookmarkPress(self,action):
		""
		
# This is just for testing, of course
if __name__ == '__main__':
	em_app = EMApp()
	window = EMBrowserWidget(withmodal=True)
		
	window.show()
	sys.exit(em_app.exec_())
		