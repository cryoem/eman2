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
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import Qt
from emapplication import EMApp
from EMAN2 import *
import os.path
import traceback

class EMDirEntry:
	"""Represents a directory entry in the filesystem"""
	
	# list of lambda functions to extract column values for sorting
	col=(lambda x:x.name,lambda x:x.filetype,lambda x:x.size,lambda x:x.dim,lambda x:x.nimg,lambda x:x.date)
#	classcount=0
	
	def __init__(self,root,name,parent=None):
		self.__parent=parent	# single parent
		self.__children=None	# ordered list of children, None indicates no check has been made, empty list means no children, otherwise list of names or list of EMDirEntrys
		self.root=str(root)		# Path prefixing name
		self.name=str(name)		# name of this path element (string)
		#self.seq=EMDirEntry.classcount
		#EMDirEntry.classcount+=1
		
		if name[:4].lower=="bdb:":
			self.isbdb=True
			self.name=name[4:]
		else:
			self.isbdb=False
			
		if self.isbdb :
			self.filepath=os.path.join(self.root,"EMAN2DB",self.name+".bdb")
		else : self.filepath=os.path.join(self.root,self.name)
		stat=os.stat(self.filepath)
		
		self.size=stat[6]		# file size (integer, bytes)
		self.date=local_datetime(stat[8])	# modification date (string: yyyy/mm/dd hh:mm:ss)

		# These can be expensive so we only get them on request, or if they are fast 
		self.dim=None			# dimensions (string)
		self.filetype=None		# file type (string, "Folder" for directory)
		self.nimg=None			# number of images in file (int)

		# Directories don't really have extra info
		if os.path.isdir(self.filepath) :
			self.filetype="Folder"
			self.dim=""
			self.nimg=""
			self.size=0			# more convenient to show as zero
		# BDB details are cached and can be retrieved quickly
		elif self.isbdb:
			self.filetype="BDB"
			try:
				d=db_open_dict("bdb:%s#00image_counts",True)
				p=d[name]
				self.nimg=p[1]
				if p[2][1]==1 : self.dim=str(p[2][0])
				elif p[2][2]==1 : self.dim="%dx%d"%(p[2][0],p[2][1])
				else : self.dim="%dx%dx%d"%(p[2][0],p[2][1],p[2][2])
			except:
				self.nimg=-1
				self.dim="-"
		
#		print "Init DirEntry ",self,self.__dict__
		
	#def __repr__(self):
		#return "<EMDirEntry %d>"%self.seq

	#def __str__(self):
		#return "<EMDirEntry %d>"%self.seq
	
	def sort(self,column,order):
		"Recursive sorting"
		if self.__children==None or len(self.__children)==0 or isinstance(self.__children[0],str): return
		
		self.__children.sort(key=EMDirEntry[column])
		
	def parent(self):
		"""Return the parent"""
		return self.__parent

	def child(self,n):
		"""Returns nth child or None"""
		self.fillChildEntries()
		return self.__children[n]
	
	def nChildren(self):
		"""Count of children"""
		self.fillChildNames()
#		print "EMDirEntry.nChildren(%s) = %d"%(self.filepath,len(self.__children))
		return len(self.__children)
	
	def fillChildNames(self):
		"""Makes sure that __children contains at LEAST a list of names"""
		if self.__children==None:
			
			if not os.path.isdir(self.filepath) :
				self.__children=[]
				return
			
			self.__children=os.listdir(self.filepath)
			if "EMAN2DIR" in self.__children :
				self.__children.remove("EMAN2DIR")
				
				t=["bdb:"+i for i in db_list_dicts(self.filepath)]
				self.__children.extend(t)
				
			self.__children.sort()
				
	def fillChildEntries(self):
		"""Makes sure that __children have been filled with EMDirEntries when appropriate"""
		if self.__children == None : self.fillChildNames()
		if len(self.__children)==0 : return		# nothing to do. No children
		if not isinstance(self.__children[0],str) : return 		# this implies entries have already been filled
		
		for i,n in enumerate(self.__children):
			self.__children[i]=EMDirEntry(self.filepath,n,self)
	
	def fillDetails(self):
		"""Fills in the expensive metadata about this entry"""
		if self.filetype!=None : return
		
		# FIXME - finish this
		
class EMFileItemModel(QtCore.QAbstractItemModel):
	"""This ItemModel represents the local filesystem. We don't use the normal filesystem item model because we want
	to provide more info on images, and we need to merge BDB: files into the file view."""
	
	headers=("Name","Type","Size","Dim","N","Date")
	
	def __init__(self,startpath=None):
		QtCore.QAbstractItemModel.__init__(self)
		self.root=EMDirEntry(startpath,"")					# EMDirEntry as a parent for the root path
		self.rootpath=startpath			# root path for current browser
		self.last=(0,0)
#		print "Init FileItemModel ",self,self.__dict__

	def canFetchMore(self,idx):
		return False
		
	def columnCount(self,parent):
		#print "EMFileItemModel.columnCount()=6"
		return 6
		
	def rowCount(self,parent):
		if parent.isValid() : 
#			print "rowCount(%s) = %d"%(str(parent),parent.internalPointer().nChildren())
			return parent.internalPointer().nChildren()
			
#		print "rowCount(root) = %d"%self.root.nChildren()
		return self.root.nChildren()
		
	def data(self,index,role):
		"QModelIndex, Qt::DisplayRole"
		
		if not index.isValid() : return None
		if role!=Qt.DisplayRole : return None
		
		data=index.internalPointer()
		#if index.column()==0 : print "EMFileItemModel.data(%d %d %s)=%s"%(index.row(),index.column(),index.parent(),str(data.__dict__))
		col=index.column()
		if col==0 : return str(data.name)
		elif col==1 : return str(data.filetype)
		elif col==2 : return str(data.size)
		elif col==3 : return str(data.dim)
		elif col==4 : return str(data.nimg)
		elif col==5 : return str(data.date)
		
	def headerData(self,sec,orient,role):
		if orient==Qt.Horizontal:
			if role==Qt.DisplayRole :
				return EMFileItemModel.headers[sec]
			elif role==Qt.ToolTipRole:
				return None								# May fill this in later
			
		return None
			
	def hasChildren(self,parent):
		#print "EMFileItemModel.hasChildren(%d,%d,%s)"%(parent.row(),parent.column(),str(parent.internalPointer()))
		try: 
			if parent.isValid():
				if parent.internalPointer().nChildren()>0 : return True
				else : return False
			return True
		except: return False
		
	def hasIndex(self,row,col,parent):
		try:
			if parent.isValid(): 
				data=parent.internalPointer().child(row)
			else: data=self.root.children(row)
		except: return False
		return True
		
	def index(self,row,column,parent):
		try:
			if parent.isValid(): 
				data=parent.internalPointer().child(row)
			else: data=self.root.child(row)
		except:
			traceback.print_exc()
			return QtCore.QModelIndex()			# No data, return invalid
		return self.createIndex(row,column,data)
		
	def parent(self,index):
		"qmodelindex"
		
		if index.isValid(): 
			data=index.internalPointer().parent()
		else: return QtCore.QModelIndex()
		if data==None : return QtCore.QModelIndex()			# No data, return invalid
		
		# Ok, this following statement appears to be crazy. When you select a whole row, the TreeView demands
		# that all selected indices have the same parent. The trick is that it doesn't just require that the
		# parent indexes point to the same parent, but that the index objects themselves be identical.
		# this little hack is intended to return sequential requests for the same linked parent object
		# as a single index, rather than making a new one each time. If this isn't done, you get a lot of
		# errors about different parents
		if index.internalPointer()==self.last[0]: return self.last[1]
		
		self.last=(index.internalPointer(),self.createIndex(index.row(),index.column(),data))
		return self.last[1]

	def sort(self,column,order):
		"Trigger recursive sorting"
		self.root.sort(column,order)
	
	
class myQItemSelection(QtGui.QItemSelectionModel):
	
	def select(self,tl,br):
		print tl.indexes()[0].row(),tl.indexes()[0].column(),int(br)
		QtGui.QItemSelectionModel.select(self,tl,QtGui.QItemSelectionModel.SelectionFlags(QtGui.QItemSelectionModel.ClearAndSelect+QtGui.QItemSelectionModel.Rows))
		

class EMBrowserWidget(QtGui.QWidget):
	"""This widget is a file browser for EMAN2. In addition to being a regular file browser, it supports:
	- getting information about recognized data types
	- embedding BDB: databases into the observed filesystem
	- remote database access (EMEN2)
	"""
	
	def __init__(self,parent=None,withmodal=False):
		QtGui.QWidget.__init__(self,parent)
		
		self.gbl = QtGui.QGridLayout(self)
		
		2# Top Toolbar area
		self.wtools=QtGui.QWidget()
		self.wtoolhbl=QtGui.QHBoxLayout(self.wtools)
		self.wtoolhbl.setContentsMargins(0,0,0,0)
		
		self.wbutback=QtGui.QPushButton("<-")
		self.wbutback.setMaximumWidth(42)
		self.wtoolhbl.addWidget(self.wbutback,0)

		self.wbutup=QtGui.QPushButton("^")
		self.wbutup.setMaximumWidth(42)
		self.wtoolhbl.addWidget(self.wbutup,0)

		self.wbutfwd=QtGui.QPushButton("->")
		self.wbutfwd.setMaximumWidth(42)
		self.wtoolhbl.addWidget(self.wbutfwd,0)

		# Text line for showing (or editing) full path
		self.wpath = QtGui.QLineEdit()
		self.wtoolhbl.addWidget(self.wpath,5)
				
		#self.wspacet1=QtGui.QSpacerItem(100,10,QtGui.QSizePolicy.MinimumExpanding)
		#self.wtoolhbl.addSpacerItem(self.wspacet1)


		self.wbutinfo=QtGui.QPushButton("Info")
		self.wbutinfo.setCheckable(True)
		self.wtoolhbl.addWidget(self.wbutinfo,1)

		self.gbl.addWidget(self.wtools,0,0,1,2)
		
		
		### Central verticalregion has bookmarks and tree
		# Bookmarks implemented with a toolbar in a frame
		self.wbookmarkfr = QtGui.QFrame()
		self.wbookmarkfr.setFrameStyle(QtGui.QFrame.StyledPanel|QtGui.QFrame.Raised)
		self.wbmfrbl=QtGui.QVBoxLayout(self.wbookmarkfr)
		
		self.wbookmarks = QtGui.QToolBar()	
		#self.wbookmarks.setAutoFillBackground(True)
		#self.wbookmarks.setBackgroundRole(QtGui.QPalette.Dark)
		self.wbookmarks.setOrientation(2)
		self.addBookmark("EMEN2","emen2")
		self.wbookmarks.addSeparator()
		self.addBookmark("Root","/")
		self.addBookmark("Current",".")
		self.addBookmark("Home",e2gethome())
		self.wbmfrbl.addWidget(self.wbookmarks)
		
		self.gbl.addWidget(self.wbookmarkfr,1,0)
		
		self.wtree = QtGui.QTreeView()
		self.wtree.setSelectionMode(1)			# single selection
		self.wtree.setSelectionBehavior(1)		# select rows
		self.wtree.setAllColumnsShowFocus(True)
		self.gbl.addWidget(self.wtree,1,1)
		
		# Lower region has buttons for actions
		self.hbl2 = QtGui.QHBoxLayout()

		self.wbutshow = QtGui.QPushButton("Show")
		self.hbl2.addWidget(self.wbutshow)
		
		self.wbutadd = QtGui.QPushButton("Add")
		self.hbl2.addWidget(self.wbutadd)
		
		self.wbutnew = QtGui.QPushButton("New")
		self.hbl2.addWidget(self.wbutnew)

		self.wbutsave = QtGui.QPushButton("Save")
		self.hbl2.addWidget(self.wbutsave)


		# buttons for selector use
		if withmodal :
			self.wspace1=QtGui.QSpacerItem(100,10,QtGui.QSizePolicy.MinimumExpanding)
			self.hbl2.addSpacerItem(self.wspace1)
			
			self.wbutcancel=QtGui.QPushButton("Cancel")
			self.hbl2.addWidget(self.wbutcancel)
			
			self.wbutok=QtGui.QPushButton("OK")
			self.hbl2.addWidget(self.wbutok)

			QtCore.QObject.connect(self.wbutcancel, QtCore.SIGNAL('clicked(bool)'), self.buttonCancel)
			QtCore.QObject.connect(self.wbutok, QtCore.SIGNAL('clicked(bool)'), self.buttonOk)

		self.gbl.addLayout(self.hbl2,3,1)

		QtCore.QObject.connect(self.wbutshow, QtCore.SIGNAL('clicked(bool)'), self.buttonShow)
		QtCore.QObject.connect(self.wbutadd, QtCore.SIGNAL('clicked(bool)'), self.buttonAdd)
		QtCore.QObject.connect(self.wbutnew, QtCore.SIGNAL('clicked(bool)'), self.buttonNew)
		QtCore.QObject.connect(self.wbutsave, QtCore.SIGNAL('clicked(bool)'), self.buttonSave)
		QtCore.QObject.connect(self.wbutback, QtCore.SIGNAL('clicked(bool)'), self.buttonBack)
		QtCore.QObject.connect(self.wbutfwd, QtCore.SIGNAL('clicked(bool)'), self.buttonFwd)
		QtCore.QObject.connect(self.wbutup, QtCore.SIGNAL('clicked(bool)'), self.buttonUp)
		QtCore.QObject.connect(self.wbutinfo, QtCore.SIGNAL('clicked(bool)'), self.buttonInfo)
		QtCore.QObject.connect(self.wtree, QtCore.SIGNAL('clicked(const QModelIndex)'), self.itemSel)
		QtCore.QObject.connect(self.wpath, QtCore.SIGNAL('returnPressed()'), self.editPath)
		QtCore.QObject.connect(self.wbookmarks, QtCore.SIGNAL('actionTriggered(QAction*)'), self.bookmarkPress)

		self.curmodel=None
		self.curpath=None

	def editPath(self):
		print "Return pressed in path editor"

	def itemSel(self,qmi):
		print "Item selected ",qmi.row(),qmi.column()
		qism=self.wtree.selectionModel()
		print qism.selectedRows()

	def buttonOk(self,tog):
		"Button press"
		
	def buttonCancel(self,tog):
		"Button press"
		
	def buttonShow(self,tog):
		"Button press"
		
	def buttonNew(self,tog):
		"Button press"

	def buttonSave(self,tog):
		"Button press"
		
	def buttonAdd(self,tog):
		"Button press"
		
	def buttonBack(self,tog):
		"Button press"
		
	def buttonFwd(self,tog):
		"Button press"
		
	def buttonUp(self,tog):
		"Button press"
		
	def buttonInfo(self,tog):
		"Button press"

	def addBookmark(self,label,path):
		"""Add a new bookmark"""
		act=self.wbookmarks.addAction(label)
		act.setData(path)
	
	def bookmarkPress(self,action):
		""
		print "Got action ",action.text(),action.data().toString()
		
		newpath=action.data().toString()
		self.curpath=newpath
		self.curmodel=EMFileItemModel(newpath)
		self.wpath.setText(newpath)
		self.wtree.setModel(self.curmodel)
#		self.wtree.setSelectionModel(myQItemSelection(self.curmodel))
		
# This is just for testing, of course
if __name__ == '__main__':
	em_app = EMApp()
	window = EMBrowserWidget(withmodal=True)
		
	window.show()
	sys.exit(em_app.exec_())
		