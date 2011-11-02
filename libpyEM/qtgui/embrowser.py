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
from emimage2d import *
from emimagemx import *
from emplot2d import *
from emplot3d import *
from emscene3d import *
from emdataitem3d import *
import re
import threading

# This is a floating point number-finding regular expression
renumfind=re.compile(r"-?[0-9]+\.*[0-9]*[eE]?[-+]?[0-9]*")


def isprint(s):
	"returns True if the string contains only printable ascii characters"
	
	# Seems like no isprint() in python, this does basically the same thing
	mpd=s.translate("AAAAAAAAABBAABAAAAAAAAAAAAAAAAAA"+"B"*95+"A"*129)
	if "A" in mpd : 
		ind=mpd.index("A")
		print "bad chr %d at %d"%(ord(s[ind]),ind)
		return False
	return True

class EMFileType:
	"""This is a base class for handling interaction with files of different type. It includes a number of excution methods common to
	several different subclasses"""

	# A class dictionary keyed by EMDirEntry filetype string with value beign a single subclass of EMFileType. filetype strings are unique
	# When you add a new EMFiletype subclass, it must be added to this dictionary to be functional
	typesbyft = {}
	
	# a class dictionary keyed by file extension with values being either a single subclass or a list of possible subclasses for that extension
	# When you add a new EMFiletype subclass, it must be added to this dictionary to be functional
	typesbyext = {}

	# A list of all types that need to be checked when the file extension can't be interpreted
	alltocheck=()

	def __init__(self,path):
		self.path=path			# the current path this FileType is representing

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

	def actions(self):
		"""Returns a list of (name,help,callback) tuples detailing the operations the user can call on the current file.
		callbacks will also be passed a reference to the browser object."""
		return []

	def plot2dApp(self,brws):
		"Append self to current plot"
		brws.busy()

		data=EMData(self.path)

		try: 
			target=brws.viewplot2d[-1]
			target.set_data(data,self.path.split("/")[-1].split("#")[-1])
		except: 
			target=EMPlot2DWidget()
			brws.viewplot2d.append(target)
			target.set_data(data,self.path.split("/")[-1].split("#")[-1])

		brws.notbusy()
		target.show()
		target.raise_()
		
	def plot2dNew(self,brws):
		"Make a new plot"
		brws.busy()

		data=EMData(self.path)
		
		target=EMPlot2DWidget()
		brws.viewplot2d.append(target)
		target.set_data(data,self.path.split("/")[-1].split("#")[-1])

		brws.notbusy()
		target.show()
		target.raise_()
		

	def show3dApp(self,brws):
		"Add to current 3-D window"
		brws.busy()

		data=EMDataItem3D(self.path)

		try: 
			target=brws.view3d[-1]
		except: 
			target=EMScene3D()
			brws.view3d.append(target)

		target.insertNewNode(self.path.split("/")[-1].split("#")[-1],data)
		iso = EMIsosurface(data)
		target.insertNewNode('Isosurface', iso, parentnode=data)
		brws.notbusy()
		target.show()
		target.raise_()

	def show3DNew(self,brws):
		"New 3-D window"
		brws.busy()

		data=EMDataItem3D(self.path)

		target=EMScene3D()
		brws.view3ds.append(target)
		
		target.insertNewNode(self.path.split("/")[-1].split("#")[-1],data)
		iso = EMIsosurface(data)
		target.insertNewNode('Isosurface', iso, parentnode=data)
		brws.notbusy()
		
		target.show()
		target.raise_()
		
	def show2dStack(self,brws):
		"A set of 2-D images together in an existing window"
		brws.busy()
		if self.dim[2]>1:
			data=[]
			for z in range(self.dim[2]):
				data.append(EMData(self.path,0,False,Region(0,0,z,self.dim[0],self.dim[1],1)))
		else : data=EMData.read_images(self.path)
		
		try: 
			target=brws.view2ds[-1]
			target.set_data(data)
		except: 
			target=EMImageMXWidget()
			target.set_data(data)
			brws.view2ds.append(target)
			
		brws.notbusy()
		target.show()
		target.raise_()
		
	def show2dStackNew(self,brws):
		"A set of 2-D images together in a new window"
		brws.busy()
		if self.dim[2]>1:
			data=[]
			for z in range(self.dim[2]):
				data.append(EMData(self.path,0,False,Region(0,0,z,self.dim[0],self.dim[1],1)))
		else : data=EMData.read_images(self.path)
		
		target=EMImageMXWidget()
		target.set_data(data)
		brws.view2ds.append(target)

		brws.notbusy()
		target.show()
		target.raise_()
		
	def show2dSingle(self,brws):
		"Show a single 2-D image"
		brws.busy()
		if self.nimg>1 : data=EMData.read_images(self.path)
		else : data=EMData(self.path)
		
		try: 
			target=brws.view2d[-1]
			target.set_data(data)
		except: 
			target=EMImage2DWidget(data)
			brws.view2d.append(target)

		brws.notbusy()
		target.show()
		target.raise_()

	def show2dSingleNew(self,brws):
		"Show a single 2-D image"
		brws.busy()
		if self.nimg>1 : data=EMData.read_images(self.path)
		else : data=EMData(self.path)
		
		target=EMImage2DWidget(data)
		brws.view2d.append(target)

		brws.notbusy()
		target.show()
		target.raise_()
		
		
	def showFilterTool(self,brws):
		"Open in e2filtertool.py"
		
		os.system("e2filtertool.py %s &"%self.path)


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

		return (size,"-",dim)

	def actions(self):
		"No actions other than the inspector for text files"
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
		
		# We need to try to count the columns in the file
		hdr=header.splitlines()
		numc=0
		for l in hdr:
			if l[0]=="#" : continue		# comment lines ok
			
			try: numc=len([float(i) for i in renumfind.findall(l)])		# number of numeric columns
			except: 
				return False		# shouldn't really happen...
				
			if numc>0 : break			# just finding the number of columns
		
		# If we couldn't find at least one valid line with some numbers, give up
		if numc==0 : return False
		
		try: size=os.stat(path)[6]
		except: return False

		# Make sure all of the lines have the same number of columns
		fin=file(path,"r")
		numr=0
		for l in fin:
			if l[0]=="#" : continue
			
			lnumc=len([float(i) for i in renumfind.findall(l)])
			if lnumc!=0 and lnumc!=numc : return False				# 0 means the line contains no numbers, we'll live with that, but if there are numbers, it needs to match
			if lnumc!=0 : numr+=1
			
		return (size,"-","%d x %d"%(numr,numc))

	def __init__(self,path):
		self.path=path			# the current path this FileType is representing

		# Make sure all of the lines have the same number of columns
		fin=file(path,"r")
		numr=0
		numc=0
		for l in fin:
			if l[0]=="#" : continue
			
			lnumc=len([float(i) for i in renumfind.findall(l)])
			if lnumc!=0 and numc==0 : numc=lnumc
			elif lnumc!=0 and lnumc!=numc : 
				print "Error: invalid Plot file :",path
				self.numr=0
				self.numc=0
				return
			elif lnumc!=0 : numr+=1
			
			self.numr=numr
			self.numc=numc
			
		
		
	def actions(self):
		"""Returns a list of (name,help,callback) tuples detailing the operations the user can call on the current file.
		callbacks will also be passed a reference to the browser object."""
		
		if self.numc>2 : return [("Plot 2D+","Make new plot",self.plot2dNew),("Plot 2D","Add to current plot",self.plot2dApp),
			("Plot 3D+","Make new 3-D plot",self.plot3dNew),("Plot 3D","Add to current 3-D plot",self.plot3dApp)]
		return [("Plot 2D+","Make new plot",self.plot2dNew),("Plot 2D","Add to current plot",self.plot2dApp)]
		
	def plot2dApp(self,brws):
		"Append self to current plot"
		brws.busy()

		data1=[]
		fin=file(self.path,"r")
		numr=0
		for l in fin:
			if l[0]=="#" : continue
			data1.append([float(i) for i in renumfind.findall(l)])
		
		data=[]
		for c in xrange(self.numc):
			data.append([i[c] for i in data1])
		
		try: 
			target=brws.viewplot2d[-1]
			target.set_data(data,self.path.split("/")[-1].split("#")[-1])
		except: 
			target=EMPlot2DWidget()
			brws.viewplot2d.append(target)
			target.set_data(data,self.path.split("/")[-1].split("#")[-1])

		brws.notbusy()
		target.show()
		target.raise_()
		
	def plot2dNew(self,brws):
		"Make a new plot"
		brws.busy()

		data1=[]
		fin=file(self.path,"r")
		numr=0
		for l in fin:
			if l[0]=="#" : continue
			data1.append([float(i) for i in renumfind.findall(l)])
		
		data=[]
		for c in xrange(self.numc):
			data.append([i[c] for i in data1])
		
		target=EMPlot2DWidget()
		brws.viewplot2d.append(target)
		target.set_data(data,self.path.split("/")[-1].split("#")[-1])

		brws.notbusy()
		target.show()
		target.raise_()
		
	def plot3dApp(self,brws):
		"Append self to current 3-D plot"
		
	def plot3dNew(self,brws):
		"Make a new 3-D plot"
		

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

	def actions(self):
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

	def __init__(self,path):
		self.path=path			# the current path this FileType is representing
		self.bdb=db_open_dict(path,ro=True)
		
		# here we assume the bdb either contains numbered images OR key/value pairs. Not strictly valid,
		# but typically true
		self.nimg=len(self.bdb)
		if self.nimg==0 : self.keys=self.bdb.keys()
		else: self.keys=None
		
		if self.nimg>0:
			im0=EMData(path,0,True)
			self.dim=(im0["nx"],im0["ny"],im0["nz"])
		else : self.dim=(0,0,0)
			
	def actions(self):
		"Returns a list of (name,help,callback) tuples detailing the operations the user can call on the current file"
		
		# single 3-D
		if self.nimg==1 and self.dim[2]>1 :
			return [("Show 3D","Add to 3D window",self.show3dApp),("Show 3D+","New 3D Window",self.show3DNew),("Show Stack","Show as set of 2-D Z slices",self.show2dStack),
				("Show Stack+","Show all images together in a new window",self.show2dStackNew),("Show 2D","Show in a scrollable 2D image window",self.show2dSingle),
				("Show 2D+","Show all images, one at a time in a new window",self.show2dSingleNew),("Chimera","Open in chimera (if installed)",self.showChimera),
				("FilterTool","Open in e2filtertool.py",self.showFilterTool)]
		# single 2-D
		elif self.nimg==1 and self.dim[1]>1 :
			return [("Show 2D","Show in a 2D single image display",self.show2dSingle),("Show 2D+","Show in new 2D single image display",self.show2dSingleNew),("FilterTool","Open in e2filtertool.py",self.showFilterTool)]
		# single 1-D
		elif self.nimg==1:
			return [("Plot 2D","Add to current plot",self.plot2dApp),("Plot 2D+","Make new plot",self.plot2dNew),
				("Show 2D","Replace in 2D single image display",self.show2dSingle),("Show 2D+","New 2D single image display",self.show2dSingleNew)]
		# 3-D stack
		elif self.nimg>1 and self.dim[2]>1 :
			return [("Show 3D","Show all in a single 3D window",self.show3DNew),("Chimera","Open in chimera (if installed)",self.showChimera)]
		# 2-D stack
		elif self.nimg>1 and self.dim[1]>1 :
			return [("Show Stack","Show all images together in one window",self.show2dStack),("Show Stack+","Show all images together in a new window",self.show2dStackNew),
				("Show 2D","Show all images, one at a time in current window",self.show2dSingle),("Show 2D+","Show all images, one at a time in a new window",self.show2dSingleNew),("FilterTool","Open in e2filtertool.py",self.showFilterTool)]
		# 1-D stack
		elif self.nimg>0:
			return [("Plot 2D","Plot all on a single 2-D plot",self.plot2dNew)]
		
		return []

	def showChimera(self,brws):
		"Open in Chimera"
		
		if get_platform()=="Linux":
			os.system("e2proc3d.py %s /tmp/vol.hdf"%self.path)		# Probably not a good hack to use, but it will do for now...
			os.system("chimera /tmp/vol.hdf&")
		else : print "Sorry, I don't know how to run Chimera on this platform"
		
		
class EMImageFileType(EMFileType):
	"""FileType for files containing a single 2-D image"""

	def __init__(self,path):
		self.path=path			# the current path this FileType is representing
		self.nimg=EMUtil.get_image_count(path)
		im0=EMData(path,0,True)
		self.dim=(im0["nx"],im0["ny"],im0["nz"])
		
	@staticmethod
	def name():
		"The unique name of this FileType. Stored in EMDirEntry.filetype for each file."
		return "Image"
		
	@staticmethod
	def isValid(path,header):
		"Returns (size,n,dim) if the referenced path is a file of this type, None if not valid. The first 4k block of data from the file is provided as well to avoid unnecesary file access."
		return False

	def actions(self):
		"Returns a list of (name,callback) tuples detailing the operations the user can call on the current file"
		# single 3-D
		if  self.dim[2]>1 :
			return [("Show 3D","Add to 3D window",self.show3dApp),("Show 3D+","New 3D Window",self.show3DNew),("Show Stack","Show as set of 2-D Z slices",self.show2dStack),
				("Show Stack+","Show all images together in a new window",self.show2dStackNew),("Show 2D","Show in a scrollable 2D image window",self.show2dSingle),
				("Show 2D+","Show all images, one at a time in a new window",self.show2dSingleNew),("Chimera","Open in chimera (if installed)",self.showChimera),
				("FilterTool","Open in e2filtertool.py",self.showFilterTool)]
		# single 2-D
		elif  self.dim[1]>1 :
			return [("Show 2D","Show in a 2D single image display",self.show2dSingle),("Show 2D+","Show in new 2D single image display",self.show2dSingleNew),("FilterTool","Open in e2filtertool.py",self.showFilterTool)]
		# single 1-D
		else:
			return [("Plot 2D","Add to current plot",self.plot2dApp),("Plot 2D+","Make new plot",self.plot2dNew),
				("Show 2D","Replace in 2D single image display",self.show2dSingle),("Show 2D+","New 2D single image display",self.show2dSingleNew)]
		
	def showChimera(self,brws):
		"Open in Chimera"
		
		if get_platform()=="Linux":
			# these types are supported natively in Chimera
			if EMUtil.get_image_type("tst.hdf") in (IMAGE_HDF,IMAGE_MRC,IMAGE_SPIDER,IMAGE_SINGLE_SPIDER):
				os.system("chimera %s &"%self.path)
			else :
				os.system("e2proc3d.py %s /tmp/vol.hdf"%self.path)		# Probably not a good hack to use, but it will do for now...
				os.system("chimera /tmp/vol.hdf&")
		else : print "Sorry, I don't know how to run Chimera on this platform"
		

class EMStackFileType(EMFileType):
	"""FileType for files containing a set of 1-3D images"""

	@staticmethod
	def name():
		"The unique name of this FileType. Stored in EMDirEntry.filetype for each file."
		return "Text"
	@staticmethod
	def isValid(path,header):
		"Returns (size,n,dim) if the referenced path is a file of this type, None if not valid. The first 4k block of data from the file is provided as well to avoid unnecesary file access."
		return False

	def __init__(self,path):
		self.path=path			# the current path this FileType is representing
		self.nimg=EMUtil.get_image_count(path)
		im0=EMData(path,0,True)
		self.dim=(im0["nx"],im0["ny"],im0["nz"])
		
	def actions(self):
		"Returns a list of (name,callback) tuples detailing the operations the user can call on the current file"
		# 3-D stack
		if self.nimg>1 and self.dim[2]>1 :
			return [("Show 3D","Show all in a single 3D window",self.show3DNew),("Chimera","Open in chimera (if installed)",self.showChimera)]
		# 2-D stack
		elif self.nimg>1 and self.dim[1]>1 :
			return [("Show Stack","Show all images together in one window",self.show2dStack),("Show Stack+","Show all images together in a new window",self.show2dStackNew),
				("Show 2D","Show all images, one at a time in current window",self.show2dSingle),("Show 2D+","Show all images, one at a time in a new window",self.show2dSingleNew),("FilterTool","Open in e2filtertool.py",self.showFilterTool)]
		# 1-D stack
		elif self.nimg>1:
			return [("Plot 2D","Plot all on a single 2-D plot",self.plot2dNew)]
		else : print "Error: stackfile with <2 images ? (%s)"%self.path
		
		return []
		
	def showChimera(self,brws):
		"Open in Chimera"
		
		if get_platform()=="Linux":
			# these types are supported natively in Chimera
			if EMUtil.get_image_type("tst.hdf") in (IMAGE_HDF,IMAGE_MRC,IMAGE_SPIDER,IMAGE_SINGLE_SPIDER):
				os.system("chimera %s &"%self.path)
			else :
				os.system("e2proc3d.py %s /tmp/vol.hdf"%self.path)		# Probably not a good hack to use, but it will do for now...
				os.system("chimera /tmp/vol.hdf&")
		else : print "Sorry, I don't know how to run Chimera on this platform"

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
		
		if not self.isbdb : self.size=stat[6]		# file size (integer, bytes)
		else: self.size="-"
		self.date=local_datetime(stat[8])			# modification date (string: yyyy/mm/dd hh:mm:ss)

		# These can be expensive so we only get them on request, or if they are fast 
		self.dim=None			# dimensions (string)
		self.filetype=None		# file type (string, "Folder" for directory)
		self.nimg=None			# number of images in file (int)

		# Directories don't really have extra info
		if os.path.isdir(self.filepath) :
			self.filetype="Folder"
			self.dim=""
			self.nimg=""
			self.size=""
			
		
#		print "Init DirEntry ",self,self.__dict__
		
	def __repr__(self):
		return "<EMDirEntry %s>"%self.path()

	#def __str__(self):
		#return "<EMDirEntry %d>"%self.1seq
	
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
		
		# BDB details are already cached and can often be retrieved quickly
		if self.isbdb:
			self.filetype="BDB"
			try:
				info=db_get_image_info(self.path())
				self.nimg=info[0]
				if self.nimg>0:
					if info[1][1]==1 : self.dim=str(info[1][0])
					elif info[1][2]==1 : self.dim="%d x %d"%(info[1][0],info[1][1])
					else : self.dim="%d x %d x %d"%(info[1][0],info[1][1],info[1][2])
					self.size=info[1][0]*info[1][1]*info[1][2]*4*self.nimg
				else:
					self.dim="-"

			except:
				traceback.print_exc()
				self.nimg=-1
				self.dim="-"

			return True
		
		# we do this this way because there are so many possible image file exensions, and sometimes
		# people use a non-standard one (esp for MRC files)
		try: self.nimg=EMUtil.get_image_count(self.path())
		except: self.nimg=0
			
		# we have an image file
		if self.nimg>0 :
			try: tmp=EMData(self.path(),0,True)		# try to read an image header for the file
			except : print "Error : no first image in %s."%self.path()
			
			if tmp["ny"]==1 : self.dim=str(tmp["nx"])
			elif tmp["nz"]==1 : self.dim="%d x %d"%(tmp["nx"],tmp["ny"])
			else : self.dim="%d x %d x %d"%(tmp["nx"],tmp["ny"],tmp["nz"])
			if self.nimg==1 : self.filetype="Image"
			else : self.filetype="Image Stack"
			
		# Ok, we need to try to figure out what kind of file this is
		else:
			head = file(self.path(),"rb").read(4096)		# Most FileTypes should be able to identify themselves using the first 4K block of a file
			
			try: guesses=EMFileType.extbyft[os.path.splitext(self.path())[1]]		# This will get us a list of possible FileTypes for this extension
			except: guesses=EMFileType.alltocheck
			
#			print "-------\n",guesses
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

		if parent!=None and parent.isValid() : 
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
			if parent!=None and parent.isValid():
				if parent.internalPointer().nChildren()>0 : return True
				else : return False
			return True
		except: return False
		
	def hasIndex(self,row,col,parent):
		"Test if the specified index would exist"
#		print "EMFileItemModel.hasIndex(%d,%d,%s)"%(row,column,parent.internalPointer())
		try:
			if parent!=None and parent.isValid(): 
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
			if parent!=None and parent.isValid(): 
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
	"""For debugging"""
	
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
		
		self.resize(780,580)
		self.gbl = QtGui.QGridLayout(self)
		
		# Top Toolbar area
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
		self.hbl2 = QtGui.QGridLayout()

		self.wbutmisc=[]
	
		# 10 buttons for context-dependent actions
		self.hbl2.setRowStretch(0,1)
		self.hbl2.setRowStretch(1,1)
		
		for i in range(5):
			self.hbl2.setColumnStretch(i,2)
			for j in range(2):
				self.wbutmisc.append(QtGui.QPushButton(""))
				self.hbl2.addWidget(self.wbutmisc[-1],j,i)
				self.wbutmisc[-1].hide()
#				self.wbutmisc[-1].setEnabled(False)
				QtCore.QObject.connect(self.wbutmisc[-1], QtCore.SIGNAL('clicked(bool)'), lambda x,v=i*2+j:self.buttonMisc(v))

		self.wbutxx=QtGui.QLabel("")
		self.wbutxx.setMaximumHeight(12)
		self.hbl2.addWidget(self.wbutxx,0,6)
		self.wbutyy=QtGui.QLabel("")
		self.wbutyy.setMaximumHeight(12)
		self.hbl2.addWidget(self.wbutyy,1,6)

		# buttons for selector use
		if withmodal :
#			self.wspace1=QtGui.QSpacerItem(100,10,QtGui.QSizePolicy.MinimumExpanding)
#			self.hbl2.addSpacerItem(self.wspace1)
			
			self.wbutcancel=QtGui.QPushButton("Cancel")
			self.hbl2.addWidget(self.wbutcancel,1,7)
			
			self.wbutok=QtGui.QPushButton("OK")
			self.hbl2.addWidget(self.wbutok,1,8)

			self.hbl2.setColumnStretch(6,1)
			self.hbl2.setColumnStretch(7,1)
			self.hbl2.setColumnStretch(8,1)
			
			QtCore.QObject.connect(self.wbutcancel, QtCore.SIGNAL('clicked(bool)'), self.buttonCancel)
			QtCore.QObject.connect(self.wbutok, QtCore.SIGNAL('clicked(bool)'), self.buttonOk)

		self.gbl.addLayout(self.hbl2,3,1)

		QtCore.QObject.connect(self.wbutback, QtCore.SIGNAL('clicked(bool)'), self.buttonBack)
		QtCore.QObject.connect(self.wbutfwd, QtCore.SIGNAL('clicked(bool)'), self.buttonFwd)
		QtCore.QObject.connect(self.wbutup, QtCore.SIGNAL('hasclicked(bool)'), self.buttonUp)
		QtCore.QObject.connect(self.wbutinfo, QtCore.SIGNAL('clicked(bool)'), self.buttonInfo)
		QtCore.QObject.connect(self.wtree, QtCore.SIGNAL('clicked(const QModelIndex)'), self.itemSel)
		QtCore.QObject.connect(self.wtree, QtCore.SIGNAL('activated(const QModelIndex)'), self.itemActivate)
		QtCore.QObject.connect(self.wtree, QtCore.SIGNAL('expanded(const QModelIndex)'), self.itemExpand)
		QtCore.QObject.connect(self.wpath, QtCore.SIGNAL('returnPressed()'), self.editPath)
		QtCore.QObject.connect(self.wbookmarks, QtCore.SIGNAL('actionTriggered(QAction*)'), self.bookmarkPress)

		self.curmodel=None	# The current data model displayed in the tree
		self.curpath=None	# The path represented by the current data model
		self.curft=None		# a fileType instance for the currently hilighted object
		self.curactions=[]	# actions returned by the filtetype. Cached for speed
		self.models={}		# Cached models to avoid a lot of rereading (not sure if this is really worthwhile)

		# Each of these contains a list of open windows displaying different types of content
		# The last item in the list is considered the 'current' item for any append operations
		self.view2d=[]
		self.view2ds=[]
		self.view3d=[]
		self.viewplot2d=[]
		self.viewplot3d=[]

		# These items are used to do gradually filling in of file details for better interactivity
		self.updtimer=QTimer()		# This causes the actual display updates, which can't be done from a python thread
		QtCore.QObject.connect(self.updtimer, QtCore.SIGNAL('timeout()'), self.updateDetailsDisplay)
		self.updthreadexit=False		# when set, this triggers the update thread to exit
		self.updthread=threading.Thread(target=self.updateDetails)	# The actual thread
		self.updlist=[]				# List of QModelIndex items in need of updating
		self.redrawlist=[]			# List of QModelIndex items in need of redisplay
		self.expanded=set()			# We get multiple expand events for each path element, so we need to keep track of which ones we've updated

		self.setPath(".")	# start in the local directory
		self.updthread.start()
		self.updtimer.start(300)

	def busy(self):
		"display a busy cursor"
		QtGui.qApp.setOverrideCursor(Qt.BusyCursor)

	def notbusy(self):
		"normal arrow cursor"
		QtGui.qApp.setOverrideCursor(Qt.ArrowCursor)

	def updateDetails(self):
		"""This is spawned as a thread to gradually fill in file details in the background"""
		
		while 1:
			if self.updthreadexit : break
			if len(self.updlist)==0 :
				time.sleep(1.0)				# If there is nothing to update at the moment, we don't need to spin our wheels as much
			else:
				de=self.updlist.pop()
				if de.internalPointer().fillDetails() : 
					self.redrawlist.append(de)		# if the update changed anything, we trigger a redisplay of this entry
					time.sleep(0.07)					# prevents updates from happening too fast and slowing the machine down
#				print "### ",de.internalPointer().path()

	def updateDetailsDisplay(self):
		"""Since we can't do GUI updates from a thread, this is a timer event to update the display after the beckground thread
		gets the details for each item"""
		
		if len(self.redrawlist)==0 : return
		
		# we emit a datachanged event for each item
		for i in self.redrawlist : self.curmodel.dataChanged.emit(i,self.curmodel.createIndex(i.row(),5,i.internalPointer()))

	def editPath(self):
		"Set a new path"
		self.setPath(str(self.wpath.text()))

	def itemSel(self,qmi):
#		print "Item selected",qmi.row(),qmi.column(),qmi.internalPointer().path()
		qism=self.wtree.selectionModel().selectedRows()
		if len(qism)>1 : self.wpath.setText("<multiple select>")
		elif len(qism)==1 : 
			obj=qism[0].internalPointer()
			self.wpath.setText(obj.path())
			self.curmodel.details(qism[0])
			self.wtree.resizeColumnToContents(2)
			self.wtree.resizeColumnToContents(3)
			
			# This makes an instance of a FileType for the selected object
			ftc=obj.fileTypeClass() 
			if ftc!=None: 
				self.curft=ftc(obj.path())
				
				self.curactions=self.curft.actions()
#				print actions
				for i,b in enumerate(self.wbutmisc):
					try:
						b.setText(self.curactions[i][0])
						b.setToolTip(self.curactions[i][1])
						b.show()
#						b.setEnabled(True)
					except:
						b.hide()
#						b.setEnabled(False)
			
		
	def itemActivate(self,qmi):
#		print "Item activated",qmi.row(),qmi.column()
		itm=qmi.internalPointer()
		if itm.nChildren()>0:
			self.setPath(itm.path())
	
	def itemExpand(self,qmi):
		"Called when an item is expanded"
	
		if qmi.internalPointer().path() in self.expanded: return
		self.expanded.add(qmi.internalPointer().path())
#		print "expand ",qmi.internalPointer().path()
	
		# Otherwise we get expand events on a single-click
		if qmi.internalPointer().filetype!="Folder" : return
		
		# we add the child items to the list needing updates
		for i in xrange(self.curmodel.rowCount(qmi)-1,-1,-1):
			self.updlist.append(self.curmodel.index(i,0,qmi))

	def buttonMisc(self,num):
		"Misc Button press"
		
		print "press ",self.curactions[num][0]
		
		self.curactions[num][2](self)				# This calls the action method

	def buttonOk(self,tog):
		"Button press"
		
	def buttonCancel(self,tog):
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
		
		self.updlist=[]
		
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
		self.wtree.resizeColumnToContents(0)
		self.wtree.resizeColumnToContents(2)
		self.wtree.resizeColumnToContents(3)

		self.expanded=set()
		# we add the child items to the list needing updates
		for i in xrange(self.curmodel.rowCount(None)-1,-1,-1):
			self.updlist.append(self.curmodel.index(i,0,None))


	def bookmarkPress(self,action):
		""
		print "Got action ",action.text(),action.data().toString()
		
		self.setPath(action.data().toString())
#		self.wtree.setSelectionModel(myQItemSelection(self.curmodel))
	
	def closeEvent(self,event):
		print "Exiting"
		try: window.updthreadexit=True
		except:pass
		for w in self.view2d+self.view2ds+self.view3d+self.viewplot2d+self.viewplot3d:
			w.close()

		event.accept()
		#self.app().close_specific(self)
		self.emit(QtCore.SIGNAL("module_closed")) # this signal is important when e2ctf is being used by a program running its own eve

# This is just for testing, of course
if __name__ == '__main__':
	em_app = EMApp()
	window = EMBrowserWidget(withmodal=True)
		
	window.show()
	ret=em_app.exec_()
	try: window.updthreadexit=True
	except:pass
	sys.exit(ret)
		