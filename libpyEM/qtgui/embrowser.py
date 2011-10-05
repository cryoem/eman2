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

class EMFileType:
	"""This is an abstract base class for handling interaction with files of different type"""

	# A class dictionary keyed by EMDirEntry filetype string with value beign a single subclass of EMFileType. filetype strings are unique
	# When you add a new EMFiletype subclass, it must be added to this dictionary to be functional
	typesbyft = {}
	
	# a class dictionary keyed by file extension with values being either a single subclass or a list of possible subclasses for that extension
	# When you add a new EMFiletype subclass, it must be added to this dictionary to be functional
	typesbyext = {}

	# A list of all types that need to be checked when the file extension can't be interpreted
	alltocheck=()

	def __init__(self,path,header):
		self.path=None			# the current path this FileType is representing

	def setFile(self,path):
		"""Represent a new file. Will update inspector if open. Assumes isValid already checked !"""
		self.path=path

	@staticmethod
	def name():
		"The unique name of this FileType. Stored in EMDirEntry.filetype for each file."
		return None

	@staticmethod
	def isValid(path,header):
		"Returns (size,n,dim) if the referenced path is a file of this type, false if not valid. The first 4k block of data from the file is provided as well to avoid unnecesary file access."
		return False

	def menuItems(self):
		"Returns a list of (name,callback) tuples detailing the operations the user can call on the current file"
		return []
		
def isprint(s):
	"returns True if the string contains only printable ascii characters"
	
	# Seems like no isprint() in python, this does basically the same thing
	mpd=s.translate("AAAAAAAAAABAABAAAAAAAAAAAAAAAAAA"+"B"*95+"A"*129)
	if "A" in mpd : return False
	return True

class EMTextFileType(EMFileType):
	"""FileType for files containing normal ASCII text"""
	
	@staticmethod
	def name():
		"The unique name of this FileType. Stored in EMDirEntry.filetype for each file."
		return "Text"

	@staticmethod
	def isValid(path,header):
		"Returns (size,n,dim) if the referenced path is a file of this type, None if not valid. The first 4k block of data from the file is provided as well to avoid unnecesary file access."

		if not isprint(header) : return False			# demand printable Ascii. FIXME: what about unicode ?

		try: size=os.stat(path)[6]
		except: return False
		
		if size>5000000 : dim="big"
		else :
			f=file(path,"r").read()
			lns=max(f.count("\n"),f.count("\r"))
			dim="%d ln"%lns

		return (size,"-",lns)

	def menuItems(self):
		"Returns a list of (name,callback) tuples detailing the operations the user can call on the current file"
		return []
		

class EMPlotFileType(EMFileType):
	"""FileType for files containing normal ASCII text"""
	
	@staticmethod
	def name():
		"The unique name of this FileType. Stored in EMDirEntry.filetype for each file."
		return "Plot"

	@staticmethod
	def isValid(path,header):
		"Returns (size,n,dim) if the referenced path is a file of this type, None if not valid. The first 4k block of data from the file is provided as well to avoid unnecesary file access."

		if not isprint(header): return False
		
		hdr=header.splitlines()
		for l in hdr:
			if l[0]=="#" : continue		# comment lines ok
			
		# TODO : Finish this
		
		return True

	def menuItems(self):
		"Returns a list of (name,callback) tuples detailing the operations the user can call on the current file"
		return []
		

class EMFolderFileType(EMFileType):
	"""FileType for Folders"""

	@staticmethod
	def name():
		"The unique name of this FileType. Stored in EMDirEntry.filetype for each file."
		return "Folder"
	@staticmethod
	def isValid(path,header):
		"Returns (size,n,dim) if the referenced path is a file of this type, None if not valid. The first 4k block of data from the file is provided as well to avoid unnecesary file access."
		return False

	def menuItems(self):
		"Returns a list of (name,callback) tuples detailing the operations the user can call on the current file"
		return []

class EMBdbFileType(EMFileType):
	"""FileType for Folders"""

	@staticmethod
	def name():
		"The unique name of this FileType. Stored in EMDirEntry.filetype for each file."
		return "BDB"
		
	@staticmethod
	def isValid(path,header):
		"Returns (size,n,dim) if the referenced path is a file of this type, None if not valid. The first 4k block of data from the file is provided as well to avoid unnecesary file access."
		return False

	def menuItems(self):
		"Returns a list of (name,callback) tuples detailing the operations the user can call on the current file"
		return []

class EMImageFileType(EMFileType):
	"""FileType for files containing a single 2-D image"""

	@staticmethod
	def name():
		"The unique name of this FileType. Stored in EMDirEntry.filetype for each file."
		return "Image"
		
	@staticmethod
	def isValid(path,header):
		"Returns (size,n,dim) if the referenced path is a file of this type, None if not valid. The first 4k block of data from the file is provided as well to avoid unnecesary file access."
		return False

	def menuItems(self):
		"Returns a list of (name,callback) tuples detailing the operations the user can call on the current file"
		return []
		

class EMStackFileType(EMFileType):
	"""FileType for files containing a set of 2-D images"""

	@staticmethod
	def name():
		"The unique name of this FileType. Stored in EMDirEntry.filetype for each file."
		return "Text"
	@staticmethod
	def isValid(path,header):
		"Returns (size,n,dim) if the referenced path is a file of this type, None if not valid. The first 4k block of data from the file is provided as well to avoid unnecesary file access."
		return False

	def menuItems(self):
		"Returns a list of (name,callback) tuples detailing the operations the user can call on the current file"
		return []
		

# These are set all together at the end rather than after each class for efficiency
EMFileType.typesbyft = {
	"Folder":EMFolderFileType,
	"BDB":EMBdbFileType,
	"Text":EMTextFileType,
	"Plot":EMPlotFileType,
	"Image":EMImageFileType,
	"Image Stack":EMStackFileType
}

# Note that image types are not included here, and are handled via a separate mechanism
# note that order is important in the tuples. The most specific filetype should go first, and
# the most general last (they will be checked in order)
EMFileType.extbyft = {
	".txt":(EMPlotFileType,EMTextFileType),
	
}

# Default Choices when extension doesn't work
# We don't need to test for things like Images because they are fully tested outside this mechanism
EMFileType.alltocheck = (EMPlotFileType, EMTextFileType)


class EMDirEntry:
	"""Represents a directory entry in the filesystem"""
	
	# list of lambda functions to extract column values for sorting
	col=(lambda x:x.name,lambda x:x.filetype,lambda x:x.size,lambda x:x.dim,lambda x:x.nimg,lambda x:x.date)
#	classcount=0
	
	def __init__(self,root,name,parent=None,hidedot=True):
		"""The path for this item is root/name. 
		Parent (EMDirEntry) must be specified if it exists.
		hidedot will cause hidden files (starting with .) to be excluded"""
		self.__parent=parent	# single parent
		self.__children=None	# ordered list of children, None indicates no check has been made, empty list means no children, otherwise list of names or list of EMDirEntrys
		self.root=str(root)		# Path prefixing name
		self.name=str(name)		# name of this path element (string)
		self.hidedot=hidedot	# If set children beginning with . will be hidden
		if self.root[-1]=="/" or self.root[-1]=="\\" : self.root=self.root[:-1]
		#self.seq=EMDirEntry.classcount
		#EMDirEntry.classcount+=1
		
		if name[:4].lower()=="bdb:":
			self.isbdb=True
			self.name=name[4:]
		else:
			self.isbdb=False
			
		if self.isbdb :
			self.filepath=os.path.join(self.root,"EMAN2DB",self.name+".bdb")
		else : self.filepath=os.path.join(self.root,self.name)
		
		try: stat=os.stat(self.filepath)
		except : stat=(0,0,0,0,0,0,0,0,0)
		
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
			
		# BDB details are already cached and can be retrieved quickly
		elif self.isbdb:
			self.filetype="BDB"
			try:
				info=db_get_image_info(self.path())
				self.nimg=info[0]
				if self.nimg>0:
					if info[1][1]==1 : self.dim=str(info[1][0])
					elif info[1][2]==1 : self.dim="%dx%d"%(info[1][0],info[1][1])
					else : self.dim="%dx%dx%d"%(info[1][0],info[1][1],info[1][2])
					self.size=info[1][0]*info[1][1]*info[1][2]*4*self.nimg
				else:
					self.dim="-"

				#d=db_open_dict("bdb:%s#00image_counts"%root,True)
				#p=d[self.name]
				#print self.name,p,root
				#self.nimg=p[1]
				#if p[2][1]==1 : self.dim=str(p[2][0])
				#elif p[2][2]==1 : self.dim="%dx%d"%(p[2][0],p[2][1])
				#else : self.dim="%dx%dx%d"%(p[2][0],p[2][1],p[2][2])
			except:
				traceback.print_exc()
				self.nimg=-1
				self.dim="-"
		
#		print "Init DirEntry ",self,self.__dict__
		
	def __repr__(self):
		return "<EMDirEntry %s>"%self.path()

	#def __str__(self):
		#return "<EMDirEntry %d>"%self.seq
	
	def path(self):
		"""The full path of the current item"""
		if self.isbdb: return "bdb:%s#%s"%(self.root,self.name)
		return os.path.join(self.root,self.name)
	
	def fileTypeClass(self):
		"Returns the FileType class corresponding to the named filetype if it exists. None otherwise"
		try: return EMFileType.typesbyft[self.filetype]
		except: return None
	
	def sort(self,column,order):
		"Recursive sorting"
		if self.__children==None or len(self.__children)==0 or isinstance(self.__children[0],str): return
		
		self.__children.sort(key=EMDirEntry.col[column],reverse=order)
		
	def parent(self):
		"""Return the parent"""
		return self.__parent

	def child(self,n):
		"""Returns nth child or None"""
		self.fillChildEntries()
		try: return self.__children[n]
		except:
			print "Request for child %d of children %s (%d)"%(n,self.__children,len(self.__children))
			traceback.print_stack()
			raise Exception
	
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
			
			# read the child filenames
			if self.hidedot : self.__children=[i for i in os.listdir(self.filepath) if i[0]!='.']
			else : self.__children=os.listdir(self.filepath)
			
			if "EMAN2DB" in self.__children :
				self.__children.remove("EMAN2DB")
				
				t=["bdb:"+i for i in db_list_dicts("bdb:"+self.filepath)]
				self.__children.extend(t)
				
			self.__children.sort()
			
#			print self.path(),self.__children
			
	def fillChildEntries(self):
		"""Makes sure that __children have been filled with EMDirEntries when appropriate"""
		if self.__children == None : self.fillChildNames()
		if len(self.__children)==0 : return		# nothing to do. No children
		if not isinstance(self.__children[0],str) : return 		# this implies entries have already been filled
		
		for i,n in enumerate(self.__children):
			self.__children[i]=EMDirEntry(self.filepath,n,self)
	
	def fillDetails(self):
		"""Fills in the expensive metadata about this entry. Returns False if no update was necessary."""
		if self.filetype!=None : return False		# must all ready be filled in
		
		try: self.nimg=EMUtil.get_image_count(self.path())
		except: self.nimg=0
			
		# we have an image file
		if self.nimg>0 :
			try: tmp=EMData(self.path(),0,True)		# try to read an image header for the file
			except : print "Error : no first image in %s."%self.path()
			
			if tmp[ny]==1 : self.dim=str(tmp[ny])
			elif info[1][2]==1 : self.dim="%dx%d"%(info[1][0],info[1][1])
			else : self.dim="%dx%dx%d"%(info[1][0],info[1][1],info[1][2])
			
		# Ok, we need to try to figure out what kind of file this is
		else:
			head = file(self.path(),"rb").read(4096)		# Most FileTypes should be able to identify themselves using the first 4K block of a file
			
			try: guesses=EMFileType.extbyft[os.path.splitext(self.path())[1]]		# This will get us a list of possible FileTypes for this extension
			except: guesses=EMFileType.alltocheck
			
			for guess in guesses:
				try : size,n,dim=guess.isValid(self.path(),head)		# This will raise an exception if isValid returns False
				except: continue
				
				# If we got here, we found a match
				self.filetype=guess.name()
				self.dim=dim
				self.nimg=n
				self.size=size
				break	
			else:		# this only happens if no match was found
				self.filetype="-"
				self.dim="-"
				self.nimg="-"
				
		return True

def nonone(val):
	"Returns '-' for None, otherwise the string representation of the passed value"
	try : 
		if val!=None: return str(val)
		return "-"
	except: return "X"

def humansize(val):
	"Representation of an integer in readable form"
	try : val=int(val)
	except: return val
	
	if val>1000000000 : return "%d g"%(val/1000000000)
	elif val>1000000 : return "%d m"%(val/1000000)
	elif val>1000 : return "%d k"%(val/1000)
	return str(val)

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
		"""The data is loaded lazily for children already. We don't generally 
	need to worry about there being SO many file that this is necessary, so it always
	returns False."""
		return False
		
	def columnCount(self,parent):
		"Always 6 columns"
		#print "EMFileItemModel.columnCount()=6"
		return 6
		
	def rowCount(self,parent):
		"Returns the number of children for a given parent"
#		if parent.column() !=0 : return 0

		if parent.isValid() : 
#			print "rowCount(%s) = %d"%(str(parent),parent.internalPointer().nChildren())
			return parent.internalPointer().nChildren()
			
#		print "rowCount(root) = %d"%self.root.nChildren()
		return self.root.nChildren()
	
	def data(self,index,role):
		"Returns the data for a specific location as a string"
		
		if not index.isValid() : return None
		if role!=Qt.DisplayRole : return None
		
		data=index.internalPointer()
		if data==None : 
			print "Error with index ",index.row(),index.column()
			return "XXX"
		#if index.column()==0 : print "EMFileItemModel.data(%d %d %s)=%s"%(index.row(),index.column(),index.parent(),str(data.__dict__))

		col=index.column()
		if col==0 : 
			if data.isbdb : return "bdb:"+data.name
			return nonone(data.name)
		elif col==1 : return nonone(data.filetype)
		elif col==2 : return humansize(data.size)
		elif col==3 :
			if data.dim==0 : return "-"
			return nonone(data.dim)
		elif col==4 : return nonone(data.nimg)
		elif col==5 : return nonone(data.date)
		
	def headerData(self,sec,orient,role):
		if orient==Qt.Horizontal:
			if role==Qt.DisplayRole :
				return EMFileItemModel.headers[sec]
			elif role==Qt.ToolTipRole:
				return None								# May fill this in later
			
		return None
			
	def hasChildren(self,parent):
		"Returns whether the index 'parent' has any child items"
		#print "EMFileItemModel.hasChildren(%d,%d,%s)"%(parent.row(),parent.column(),str(parent.internalPointer()))
#		if parent.column()!=0 : return False
		try: 
			if parent.isValid():
				if parent.internalPointer().nChildren()>0 : return True
				else : return False
			return True
		except: return False
		
	def hasIndex(self,row,col,parent):
		"Test if the specified index would exist"
#		print "EMFileItemModel.hasIndex(%d,%d,%s)"%(row,column,parent.internalPointer())
		try:
			if parent.isValid(): 
				data=parent.internalPointer().child(row)
			else: data=self.root.child(row)
		except:
			traceback.print_exc()
			return False
		return True
		
	def index(self,row,column,parent):
		"produces a new QModelIndex for the specified item"
#		if column==0 : print "Index :",row,column,parent.internalPointer(),
		try:
			if parent.isValid(): 
				data=parent.internalPointer().child(row)
			else: 
				data=self.root.child(row)
		except:
			traceback.print_exc()
#			print "None"
			return QtCore.QModelIndex()			# No data, return invalid
#		if column==0 :print data
		return self.createIndex(row,column,data)
		
	def parent(self,index):
		"Returns the parent of the specified index"
		
		if index.isValid(): 
			try: data=index.internalPointer().parent()
			except:
				print "Parent index error: ",str(index.__dict__)
			
		else: return QtCore.QModelIndex()
		if data==None : return QtCore.QModelIndex()			# No data, return invalid
				
		return self.createIndex(index.row(),0,data)		# parent is always column 0

	def sort(self,column,order):
		"Trigger recursive sorting"
		if column<0 : return
		self.root.sort(column,order)
#		self.emit(QtCore.SIGNAL("layoutChanged()"))
		self.layoutChanged.emit()
		
	def details(self,index):
		"""This will trigger loading the (expensive) details about the specified index, and update the display"""
		if not index.isValid(): return
		
		if index.internalPointer().fillDetails() : 
			self.dataChanged.emit(index, self.createIndex(index.row(),5,index.internalPointer()))
		
			
	
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
		self.wtree.sortByColumn(-1,0)			# start unsorted
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
		QtCore.QObject.connect(self.wbutup, QtCore.SIGNAL('hasclicked(bool)'), self.buttonUp)
		QtCore.QObject.connect(self.wbutinfo, QtCore.SIGNAL('clicked(bool)'), self.buttonInfo)
		QtCore.QObject.connect(self.wtree, QtCore.SIGNAL('clicked(const QModelIndex)'), self.itemSel)
		QtCore.QObject.connect(self.wtree, QtCore.SIGNAL('activated(const QModelIndex)'), self.itemActivate)
		QtCore.QObject.connect(self.wpath, QtCore.SIGNAL('returnPressed()'), self.editPath)
		QtCore.QObject.connect(self.wbookmarks, QtCore.SIGNAL('actionTriggered(QAction*)'), self.bookmarkPress)

		self.curmodel=None	# The current data model displayed in the tree
		self.curpath=None	# The path represented by the current data model
		self.curft=None		# a fileType instance for the currently hilighted object
		self.models={}

	def editPath(self):
		print "Return pressed in path editor"

	def itemSel(self,qmi):
#		print "Item selected",qmi.row(),qmi.column(),qmi.internalPointer().path()
		qism=self.wtree.selectionModel().selectedRows()
		if len(qism)>1 : self.wpath.setText("<multiple select>")
		elif len(qism)==1 : 
			obj=qism[0].internalPointer()
			self.wpath.setText(obj.path())
			self.curmodel.details(qism[0])
			
			# This makes an instance of a FileType for the selected object
			ftc=obj.fileTypeClass() 
			if ftc!=None: 
				self.curft=ftc(obj.path())
				# TODO continue here
			
		
	def itemActivate(self,qmi):
#		print "Item activated",qmi.row(),qmi.column()
		itm=qmi.internalPointer()
		if itm.nChildren()>0:
			self.setPath(itm.path())


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

	def setPath(self,path):
		"""Sets the current root path for the browser window"""
		self.curpath=path
		self.wpath.setText(path)

		if path in self.models :
			self.curmodel=self.models[path]
		else : 
			self.curmodel=EMFileItemModel(path)
			self.models[self.curpath]=self.curmodel

		self.wtree.setSortingEnabled(False)
		self.wtree.setModel(self.curmodel)
		self.wtree.setSortingEnabled(True)
		
	def bookmarkPress(self,action):
		""
		print "Got action ",action.text(),action.data().toString()
		
		self.setPath(action.data().toString())
#		self.wtree.setSelectionModel(myQItemSelection(self.curmodel))
		
# This is just for testing, of course
if __name__ == '__main__':
	em_app = EMApp()
	window = EMBrowserWidget(withmodal=True)
		
	window.show()
	sys.exit(em_app.exec_())
		