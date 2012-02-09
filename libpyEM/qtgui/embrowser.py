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
import os
import traceback
from emimage2d import *
from emimagemx import *
from emplot2d import *
from emplot3d import *
from valslider import StringBox
import re
import threading
import time

# This is a floating point number-finding regular expression
renumfind=re.compile(r"-?[0-9]+\.*[0-9]*[eE]?[-+]?[0-9]*")


def isprint(s):
	"returns True if the string contains only printable ascii characters"
	
	# Seems like no isprint() in python, this does basically the same thing
	mpd=s.translate("AAAAAAAAABBAABAAAAAAAAAAAAAAAAAA"+"B"*95+"A"*129)
	if "A" in mpd : 
		ind=mpd.index("A")
#		print "bad chr %d at %d"%(ord(s[ind]),ind)
		return False
	return True

class EMFileType(object):
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

	setsmode = False
	
	def __init__(self,path):
		self.path=path			# the current path this FileType is representing
		self.setsdb=None		# The bad particles DB 

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

	@staticmethod
	def infoClass():
		"Returns a reference to the QWidget subclass for displaying information about this file"
		return EMInfoPane
		
	def setSetsDB(self, db_name):
		"Sets the emmxwidget to sets dbname"
		if db_name[:4].lower()!="bdb:": 
			self.setsdb = os.path.splitext(os.path.basename(db_name))[0]
		else:
			self.setsdb = db_parse_path(db_name)[1]
		
	def getSetsDB(self):
		"Returns the sets mode"
		return self.setsdb

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

		data=emdataitem3d.EMDataItem3D(self.path)

		try: 
			target=brws.view3d[-1]
		except: 
			target=emscene3d.EMScene3D()
			brws.view3d.append(target)

		target.insertNewNode(self.path.split("/")[-1].split("#")[-1],data,parentnode=target)
		iso = emdataitem3d.EMIsosurface(data)
		target.insertNewNode('Isosurface', iso, parentnode=data)
		brws.notbusy()
		target.show()
		target.raise_()

	def show3DNew(self,brws):
		"New 3-D window"
		brws.busy()

		data=emdataitem3d.EMDataItem3D(self.path)

		target=emscene3d.EMScene3D()
		brws.view3d.append(target)
		
		target.insertNewNode(self.path.split("/")[-1].split("#")[-1],data)
		iso = emdataitem3d.EMIsosurface(data)
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
			if self.getSetsDB(): target.set_single_active_set(self.getSetsDB())
		except: 
			target=EMImageMXWidget()
			target.set_data(data)
			if self.getSetsDB(): target.set_single_active_set(self.getSetsDB())
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
		if self.getSetsDB(): target.set_single_active_set(self.getSetsDB())
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

	@staticmethod
	def infoClass():
		"Returns a reference to the QWidget subclass for displaying information about this file"
		return EMTextInfoPane
		
	def actions(self):
		"No actions other than the inspector for text files"
		return []
		
class EMHTMLFileType(EMFileType):
	"""FileType for files containing HTML text"""
	
	@staticmethod
	def name():
		"The unique name of this FileType. Stored in EMDirEntry.filetype for each file."
		return "HTML"

	@staticmethod
	def isValid(path,header):
		"Returns (size,n,dim) if the referenced path is a file of this type, None if not valid. The first 4k block of data from the file is provided as well to avoid unnecesary file access."

		if not isprint(header) : return False			# demand printable Ascii. FIXME: what about unicode ?
		if not "<html>" in header.lower() : return False # For the moment, we demand an <html> tag somewhere in the first 4k

		try: 
			size=os.stat(path)[6]
		except: return False
		
		if size>5000000 : dim="big"
		else :
			f=file(path,"r").read()
			lns=max(f.count("\n"),f.count("\r"))
			dim="%d ln"%lns

		return (size,"-",dim)

	@staticmethod
	def infoClass():
		"Returns a reference to the QWidget subclass for displaying information about this file"
		return EMHTMLInfoPane
		
	def actions(self):
		"No actions other than the inspector for HTML files"
		return [("Firefox","Open in Firefox",self.showFirefox)]

	def showFirefox(self,brws):
		"""Try to open file in firefox"""
		
		if get_platform()=="Linux":
				os.system("firefox -new-tab file://%s"%self.path)

		else : print "Sorry, I don't know how to run Firefox on this platform"


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
				
			if l[0]=="#" or len(l)<2 or "nan" in l : continue
			
			lnumc=len([float(i) for i in renumfind.findall(l)])
			if lnumc!=0 and lnumc!=numc : return False				# 0 means the line contains no numbers, we'll live with that, but if there are numbers, it needs to match
			if lnumc!=0 : numr+=1
		
		return (size,"-","%d x %d"%(numr,numc))
		
	@staticmethod
	def infoClass():
		"Returns a reference to the QWidget subclass for displaying information about this file"
		return EMPlotInfoPane

	def __init__(self,path):
		EMFileType.__init__(self,path)	# the current path this FileType is representing

		# Make sure all of the lines have the same number of columns
		fin=file(path,"r")
		numr=0
		numc=0
		for l in fin:
			if "nan" in l : 
				print "Warning, NaN present in file"
				continue
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
			if l[0]=="#" or "nan" in l: continue
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
			if l[0]=="#" or "nan" in l: continue
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

	@staticmethod
	def infoClass():
		"Returns a reference to the QWidget subclass for displaying information about this file"
		return EMFolderInfoPane


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

	@staticmethod
	def infoClass():
		"Returns a reference to the QWidget subclass for displaying information about this file"
		return EMBDBInfoPane

	def __init__(self,path):
		EMFileType.__init__(self,path)	# the current path this FileType is representing
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
		EMFileType.__init__(self,path)	# the current path this FileType is representing
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

	@staticmethod
	def infoClass():
		"Returns a reference to the QWidget subclass for displaying information about this file"
		return EMImageInfoPane

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
			if EMUtil.get_image_type(self.path) in (IMAGE_HDF,IMAGE_MRC,IMAGE_SPIDER,IMAGE_SINGLE_SPIDER):
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

	@staticmethod
	def infoClass():
		"Returns a reference to the QWidget subclass for displaying information about this file"
		return EMStackInfoPane

	def __init__(self,path):
		EMFileType.__init__(self,path)	# the current path this FileType is representing
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
	"Image Stack":EMStackFileType,
	"HTML":EMHTMLFileType
}

# Note that image types are not included here, and are handled via a separate mechanism
# note that order is important in the tuples. The most specific filetype should go first, and
# the most general last (they will be checked in order)
EMFileType.extbyft = {
	".txt":(EMPlotFileType,EMTextFileType),
	".htm":(EMHTMLFileType,EMTextFileType),
	".html":(EMHTMLFileType,EMTextFileType),
}

# Default Choices when extension doesn't work
# We don't need to test for things like Images because they are fully tested outside this mechanism
EMFileType.alltocheck = (EMPlotFileType, EMTextFileType)


class EMDirEntry(object):
	"""Represents a directory entry in the filesystem"""
	
	# list of lambda functions to extract column values for sorting
	col=(lambda x:int(x.index),lambda x:x.name,lambda x:x.filetype,lambda x:x.size,lambda x:x.dim,lambda x:x.nimg,lambda x:x.date)
#	classcount=0
	
	def __init__(self,root,name,index,parent=None,hidedot=True):
		"""The path for this item is root/name. 
		Parent (EMDirEntry) must be specified if it exists.
		hidedot will cause hidden files (starting with .) to be excluded"""
		self.__parent=parent	# single parent
		self.__children=None	# ordered list of children, None indicates no check has been made, empty list means no children, otherwise list of names or list of EMDirEntrys
		self.root=str(root)		# Path prefixing name
		self.name=str(name)		# name of this path element (string)
		self.index=str(index)
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
		return "<%s %s>"%(self.__class__.__name__ ,self.path())
	
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
		self.__children.sort(key=self.__class__.col[column],reverse=order)
		
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
			self.__children[i]=self.__class__(self.filepath,n,i,self)
	
	def checkCache(self, db, name=""):
		""" 
		Returns a dict from the cache if it exists AND self.path() has an access time equal to cache.
		If the cache is used by more than one browser, then a name for the nth browser must be given otherwise 
		cache info may be out of date. 'name' is needed because the same file might be involved in several MV widgets
		"""
		# modify time is used rather than access time b/c access time only has 24 hour resolution.
		# This means that ONLY file metadata should be cached and NOT associated data in DBs which will
		# not necessarily modify the original file!!!!!
		if db.has_key(self.path()) and db[self.path()+name+'lastaccesstime'] == self.statFile(self.path())[8]: 
			return db[self.path()]
		else:
			return {}
			
	def setCache(self, db, dbdict, name=""):
		""" 
		Sets a cache 
		'name' is needed because the same file might be involved in several MV widgets
		"""
		db[self.path()] = dbdict
		db[self.path()+name+'lastaccesstime'] = self.statFile(self.path())[8]
	
	def cacheMiss(self, cache, *args):
		""" Check the cahce generate a miss is some data is missing """
		if not cache: return True
		miss = False	# Basicaly a dirty bit
		for arg in args:
			if cache.has_key(arg): 
				self.__setattr__(arg, cache[arg])
			else:
				miss = True
		return miss
	
	def updateCache(self, db, cache, name="", *args):
		""" update Cache """
		dir(self)
		for arg in args:
			cache[arg] = self.__getattribute__(arg)
		self.setCache(db, cache, name)
			
	def statFile(self, filename):
		""" Stat either a file or BDB """
		if filename[:4].lower()=="bdb:":
			path,dictname,keys=db_parse_path(filename)
			path=path+"/EMAN2DB/"+dictname+".bdb"
			return os.stat(path)
		else:
			return os.stat(filename)
			
	def fillDetails(self, db):
		"""Fills in the expensive metadata about this entry. Returns False if no update was necessary."""
		if self.filetype!=None : return False		# must all ready be filled in
		
		# Check the cache for metadata
		name='browser'
		if os.access("EMAN2DB",os.R_OK):
			cache = self.checkCache(db, name=name)
			if not self.cacheMiss(cache,'filetype','nimg','dim','size'): return 
		
		# BDB details are already cached and can often be retrieved quickly
		if self.isbdb:
			self.filetype="BDB"
			try:
				info=db_get_image_info(self.path())
#				print self.path,info
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
			
			# Cache if EMAN2DB already exists
			if os.access("EMAN2DB",os.R_OK):
				self.updateCache(db, cache, name, 'filetype', 'nimg', 'dim', 'size')
			
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
		
		# Cache if EMAN2DB already exists
		if os.access("EMAN2DB",os.R_OK):
			self.updateCache(db, cache, name, 'filetype', 'nimg', 'dim', 'size')
			
		return True

	def getBaseName(self, name, extension=False):
		""" return a sensible basename """
		if name[:4].lower()!="bdb:":
			if extension:
				return os.path.basename(name)
			else:
				return os.path.splitext(os.path.basename(name))[0]
		else:
			return db_parse_path(name)[1]
			
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
	
	headers=("Row","Name","Type","Size","Dim","N","Date")
	
	def __init__(self,startpath=None,direntryclass=EMDirEntry):
		QtCore.QAbstractItemModel.__init__(self)
		self.root=direntryclass(startpath,"", 0)				# EMDirEntry as a parent for the root path
		self.rootpath=startpath							# root path for current browser
		self.last=(0,0)
		self.db=None
#		print "Init FileItemModel ",self,self.__dict__

	def canFetchMore(self,idx):
		"""The data is loaded lazily for children already. We don't generally 
	need to worry about there being SO many file that this is necessary, so it always
	returns False."""
		return False
		
	def columnCount(self,parent):
		"Always 7 columns"
		#print "EMFileItemModel.columnCount()=6"
		return 7
		
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
		if col==0:
			return nonone(data.index)
		elif col==1 : 
			if data.isbdb : return "bdb:"+data.name
			return nonone(data.name)
		elif col==2 : return nonone(data.filetype)
		elif col==3 : return humansize(data.size)
		elif col==4 :
			if data.dim==0 : return "-"
			return nonone(data.dim)
		elif col==5 : return nonone(data.nimg)
		elif col==6 : return nonone(data.date)
		
	def headerData(self,sec,orient,role):
		# In case we use the QTableViews
		if orient==Qt.Vertical:
			if role==Qt.DisplayRole :
				return str(sec)
			elif role==Qt.ToolTipRole:
				return None	
		# This works for all Q*views
		if orient==Qt.Horizontal:
			if role==Qt.DisplayRole :
				return self.__class__.headers[sec]
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
	
	def getCacheDB(self):
		if not self.db:
			self.db=db_open_dict("bdb:browsercache")
		return self.db
					
	def details(self,index):
		"""This will trigger loading the (expensive) details about the specified index, and update the display"""
		if not index.isValid(): return
		
		if index.internalPointer().fillDetails(self.getCacheDB()) : 
			self.dataChanged.emit(index, self.createIndex(index.row(),5,index.internalPointer()))
		
			
class myQItemSelection(QtGui.QItemSelectionModel):
	"""For debugging"""
	
	def select(self,tl,br):
		print tl.indexes()[0].row(),tl.indexes()[0].column(),int(br)
		QtGui.QItemSelectionModel.select(self,tl,QtGui.QItemSelectionModel.SelectionFlags(QtGui.QItemSelectionModel.ClearAndSelect+QtGui.QItemSelectionModel.Rows))

class EMInfoPane(QtGui.QWidget):
	"""Subclasses of this class will be used to display information about specific files. Each EMFileType class will return the
	pointer to the appropriate infoPane subclass for displaying information about the file it represents. The subclass instances
	are allocated by the infoWin class"""
	
	def __init__(self,parent=None):
		"Set our GUI up"
		QtGui.QWidget.__init__(self,parent)
		
		# Root class represents no target
		self.hbl=QtGui.QHBoxLayout(self)
		self.lbl=QtGui.QLabel("No Information Available")
		self.hbl.addWidget(self.lbl)

	def display(self,target):
		"display information for the target EMDirEntry with EMFileType ftype"
		
		self.target=target
		
		
		return
		
class EMTextInfoPane(EMInfoPane):
	
	def __init__(self,parent=None):
		QtGui.QWidget.__init__(self,parent)
		
		self.vbl=QtGui.QVBoxLayout(self)
		
		# text editing widget
		self.text=QtGui.QTextEdit()
		self.text.setAcceptRichText(False)
		self.text.setReadOnly(True)
		self.vbl.addWidget(self.text)
		
		# Find box
		self.wfind=StringBox(label="Find:")
		self.vbl.addWidget(self.wfind)
		
		# Buttons
		self.hbl=QtGui.QHBoxLayout()
		
		self.wbutedit=QtGui.QPushButton("Edit")
		self.hbl.addWidget(self.wbutedit)
		
		self.wbutcancel=QtGui.QPushButton("Revert")
		self.wbutcancel.setEnabled(False)
		self.hbl.addWidget(self.wbutcancel)
		
		self.wbutok=QtGui.QPushButton("Save")
		self.wbutok.setEnabled(False)
		self.hbl.addWidget(self.wbutok)
		
		self.vbl.addLayout(self.hbl)
		
		QtCore.QObject.connect(self.wfind, QtCore.SIGNAL("valueChanged"),self.find)
		QtCore.QObject.connect(self.wbutedit, QtCore.SIGNAL('clicked(bool)'), self.buttonEdit)
		QtCore.QObject.connect(self.wbutcancel, QtCore.SIGNAL('clicked(bool)'), self.buttonCancel)
		QtCore.QObject.connect(self.wbutok, QtCore.SIGNAL('clicked(bool)'), self.buttonOk)
		
	def display(self,data):
		"display information for the target EMDirEntry"
		self.target=data
		
		self.text.setPlainText(file(self.target.path(),"r").read())
		
		self.text.setReadOnly(True)
		self.wbutedit.setEnabled(True)
		self.wbutcancel.setEnabled(False)
		self.wbutok.setEnabled(False)

	def find(self,value):
		"Find a string"
		
		if not self.text.find(value):
			# this implements wrapping
			self.text.moveCursor(1,0)
			self.text.find(value)

	def buttonEdit(self,tog):
		self.text.setReadOnly(False)
		self.wbutedit.setEnabled(False)
		self.wbutcancel.setEnabled(True)
		self.wbutok.setEnabled(True)
		
	def buttonCancel(self,tog):
		self.display(self.target)

	def buttonOk(self,tog):
		try: file(self.target.path(),"w").write(str(self.text.toPlainText()))
		except: QtGui.QMessageBox.warning(self,"Error !","File write failed")

class EMHTMLInfoPane(EMInfoPane):
	
	def __init__(self,parent=None):
		QtGui.QWidget.__init__(self,parent)
		
		self.vbl=QtGui.QVBoxLayout(self)
		
		# text editing widget
		self.text=QtGui.QTextEdit()
		self.text.setAcceptRichText(True)
		self.text.setReadOnly(True)
		self.vbl.addWidget(self.text)
		
		# Find box
		self.wfind=StringBox(label="Find:")
		self.vbl.addWidget(self.wfind)
		
		# Buttons
		self.hbl=QtGui.QHBoxLayout()
		
		self.wbutedit=QtGui.QPushButton("Edit")
		self.hbl.addWidget(self.wbutedit)
		
		self.wbutcancel=QtGui.QPushButton("Revert")
		self.wbutcancel.setEnabled(False)
		self.hbl.addWidget(self.wbutcancel)
		
		self.wbutok=QtGui.QPushButton("Save")
		self.wbutok.setEnabled(False)
		self.hbl.addWidget(self.wbutok)
		
		self.vbl.addLayout(self.hbl)
		
		QtCore.QObject.connect(self.wfind, QtCore.SIGNAL("valueChanged"),self.find)
		QtCore.QObject.connect(self.wbutedit, QtCore.SIGNAL('clicked(bool)'), self.buttonEdit)
		QtCore.QObject.connect(self.wbutcancel, QtCore.SIGNAL('clicked(bool)'), self.buttonCancel)
		QtCore.QObject.connect(self.wbutok, QtCore.SIGNAL('clicked(bool)'), self.buttonOk)
		
	def display(self,data):
		"display information for the target EMDirEntry"
		self.target=data
		
		self.text.setHtml(file(self.target.path(),"r").read())
		
		self.text.setReadOnly(True)
		self.wbutedit.setEnabled(True)
		self.wbutcancel.setEnabled(False)
		self.wbutok.setEnabled(False)

	def find(self,value):
		"Find a string"
		
		if not self.text.find(value):
			# this implements wrapping
			self.text.moveCursor(1,0)
			self.text.find(value)

	def buttonEdit(self,tog):
		self.text.setReadOnly(False)
		self.wbutedit.setEnabled(False)
		self.wbutcancel.setEnabled(True)
		self.wbutok.setEnabled(True)
		
	def buttonCancel(self,tog):
		self.display(self.target)

	def buttonOk(self,tog):
		try: file(self.target.path(),"w").write(str(self.text.toHtml()))
		except: QtGui.QMessageBox.warning(self,"Error !","File write failed")

class EMPlotInfoPane(EMInfoPane):
	
	def __init__(self,parent=None):
		QtGui.QWidget.__init__(self,parent)
		
		self.gbl=QtGui.QGridLayout(self)
		
		# List as alternate mechanism for selecting image number(s)
		self.plotdata=QtGui.QTableWidget()
		self.gbl.addWidget(self.plotdata,0,0)

	def display(self,target):
		"display information for the target EMDirEntry"
		self.target=target
		self.plotdata.clear()


		# read the data into a list of lists
		numc=0
		data=[]
		for l in file(target.path(),"r"):
			if l[0]=="#" : continue
			
			vals=[float(i) for i in renumfind.findall(l)]
			if len(vals)==0 : continue 
			
			if numc==0 : numc=len(vals)
			elif numc!=len(vals) : break
			data.append(vals)
			
			if len(data)==2500 : break			# if the table is too big, we just do a ...

		if len(data)==2500: self.plotdata.setRowCount(2501)
		else : self.plotdata.setRowCount(len(data))
		self.plotdata.setColumnCount(numc)
		self.plotdata.setVerticalHeaderLabels([str(i) for i in xrange(len(data))])
		self.plotdata.setHorizontalHeaderLabels([str(i) for i in xrange(numc)])
		
		for r in xrange(len(data)):
			for c in xrange(numc):
				self.plotdata.setItem(r,c,QtGui.QTableWidgetItem("%1.4g"%data[r][c]))

		if len(data)==2500:
			self.plotdata.setVerticalHeaderItem(2500,QtGui.QTableWidgetItem("..."))

class EMFolderInfoPane(EMInfoPane):
	
	def __init__(self,parent=None):
		QtGui.QWidget.__init__(self,parent)
		
		self.vbl=QtGui.QVBoxLayout(self)

	def display(self,target):
		"display information for the target EMDirEntry"
		self.target=target

class EMBDBInfoPane(EMInfoPane):
	
	maxim=500
	def __init__(self,parent=None):
		QtGui.QWidget.__init__(self,parent)
		
		self.gbl=QtGui.QGridLayout(self)
		
		# Spinbox for selecting image number
		self.wimnum=QtGui.QSpinBox()
		self.wimnum.setRange(0,0)
		self.gbl.addWidget(self.wimnum,0,0)
		
		# List as alternate mechanism for selecting image number(s)
		self.wimlist=QtGui.QListWidget()
		self.gbl.addWidget(self.wimlist,1,0)
		
		# Actual header contents
		self.wheadtree=QtGui.QTreeWidget()
		self.wheadtree.setColumnCount(2)
		self.wheadtree.setHeaderLabels(["Item","Value"])
		self.gbl.addWidget(self.wheadtree,0,1,2,1)
		
		self.gbl.setColumnStretch(0,1)
		self.gbl.setColumnStretch(1,4)
		
		QtCore.QObject.connect(self.wimnum, QtCore.SIGNAL("valueChanged(int)"),self.imNumChange)
		QtCore.QObject.connect(self.wimlist, QtCore.SIGNAL("itemSelectionChanged()"),self.imSelChange)
#		QtCore.QObject.connect(self.wbutedit, QtCore.SIGNAL('clicked(bool)'), self.buttonEdit)
		

	def display(self,target):
		"display information for the target EMDirEntry"
		self.target=target
		self.bdb=db_open_dict(self.target.path())
		
		# Set up image selectors for stacks
		if target.nimg==0 :
			self.wimnum.hide()
			k=self.bdb.keys()
			k.sort()
			self.wimlist.addItems(k)
			self.wimlist.show()
			self.curim=0
		else:
			self.wimnum.setRange(0,target.nimg)
			self.wimlist.clear()
			self.wimlist.addItems([str(i) for i in range(0,min(target.nimg,self.maxim))])
			if target.nimg>self.maxim: self.wimlist.addItem("...")
			self.wimnum.show()
			self.wimlist.show()

		self.wheadtree.clear()

	def imNumChange(self,num):
		"New image number"
		if num<500 : self.wimlist.setCurrentRow(num)
		else : self.showItem(num)
		
	def imSelChange(self):
		"New image selection"
		
		val=self.wimlist.currentItem().text()
		try:
			val=int(val)
			self.wimnum.setValue(val)
		except:
			val=str(val)
		
		self.showItem(val)
		
	def showItem(self,key):
		"""Shows header information for the selected item"""
		self.wheadtree.clear()
		trg=self.bdb.get_header(key)
		
		if trg==None :
			print "Warning: tried to read unavailable key: %s"%key
			#print self.bdb.keys()
			return
		
		self.addTreeItem(trg)
		
	def addTreeItem(self,trg,parent=None):
		"""(recursively) add an item to the tree"""
		itms=[]
		# Dictionaries may require recursion
		if isinstance(trg,dict):
			for k in trg.keys():
				itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((str(k),str(trg[k])))))
				if isinstance(trg[k],list) or isinstance(trg[k],tuple) or isinstance(trg[k],set) or isinstance(trg[k],dict):
					self.addTreeItem(trg[k],itms[-1])

		elif isinstance(trg,list) or isinstance(trg,tuple) or isinstance(trg,set):
			for k in trg:
				if isinstance(k,list) or isinstance(k,tuple) or isinstance(k,set) or isinstance(k,dict):
					try: itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((k.__class__.__name__,""))))
					except: itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList(("??",""))))
					self.addTreeItem(k,itms[-1])
				else:
					itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((str(k),""))))
			
		else:
			itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((str(trg),""))))

		if parent==None: 
			self.wheadtree.addTopLevelItems(itms)
			self.wheadtree.resizeColumnToContents(0)
		else : parent.addChildren(itms)


class EMImageInfoPane(EMInfoPane):
	
	maxim=500
	def __init__(self,parent=None):
		QtGui.QWidget.__init__(self,parent)
		
		self.gbl=QtGui.QGridLayout(self)
		
		# Actual header contents
		self.wheadtree=QtGui.QTreeWidget()
		self.wheadtree.setColumnCount(2)
		self.wheadtree.setHeaderLabels(["Item","Value"])
		self.gbl.addWidget(self.wheadtree,0,0)
		
		

	def display(self,target):
		"display information for the target EMDirEntry"
		self.target=target
		
		self.wheadtree.clear()
		try: trg=EMData(self.target.path(),0,True).get_attr_dict()		# read the header only, discard the emdata object
		except:
			print "Error reading: ",self.target.path(),key
			
		if trg==None :
			print "Warning: tried to read unavailable key: %s"%key
			#print self.bdb.keys()
			return
		
		self.addTreeItem(trg)
		
	def addTreeItem(self,trg,parent=None):
		"""(recursively) add an item to the tree"""
		itms=[]
		# Dictionaries may require recursion
		if isinstance(trg,dict):
			for k in trg.keys():
				itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((str(k),str(trg[k])))))
				if isinstance(trg[k],list) or isinstance(trg[k],tuple) or isinstance(trg[k],set) or isinstance(trg[k],dict):
					self.addTreeItem(trg[k],itms[-1])

		elif isinstance(trg,list) or isinstance(trg,tuple) or isinstance(trg,set):
			for k in trg:
				if isinstance(k,list) or isinstance(k,tuple) or isinstance(k,set) or isinstance(k,dict):
					try: itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((k.__class__.__name__,""))))
					except: itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList(("??",""))))
					self.addTreeItem(k,itms[-1])
				else:
					itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((str(k),""))))
			
		else:
			itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((str(trg),""))))

		if parent==None: 
			self.wheadtree.addTopLevelItems(itms)
			self.wheadtree.resizeColumnToContents(0)
		else : parent.addChildren(itms)

class EMStackInfoPane(EMInfoPane):
	
	maxim=500
	def __init__(self,parent=None):
		QtGui.QWidget.__init__(self,parent)
		
		self.gbl=QtGui.QGridLayout(self)
		
		# Spinbox for selecting image number
		self.wimnum=QtGui.QSpinBox()
		self.wimnum.setRange(0,0)
		self.gbl.addWidget(self.wimnum,0,0)
		
		# List as alternate mechanism for selecting image number(s)
		self.wimlist=QtGui.QListWidget()
		self.gbl.addWidget(self.wimlist,1,0)
		
		# Actual header contents
		self.wheadtree=QtGui.QTreeWidget()
		self.wheadtree.setColumnCount(2)
		self.wheadtree.setHeaderLabels(["Item","Value"])
		self.gbl.addWidget(self.wheadtree,0,1,2,1)
		
		self.gbl.setColumnStretch(0,1)
		self.gbl.setColumnStretch(1,4)
		
		QtCore.QObject.connect(self.wimnum, QtCore.SIGNAL("valueChanged(int)"),self.imNumChange)
		QtCore.QObject.connect(self.wimlist, QtCore.SIGNAL("itemSelectionChanged()"),self.imSelChange)
#		QtCore.QObject.connect(self.wbutedit, QtCore.SIGNAL('clicked(bool)'), self.buttonEdit)
		

	def display(self,target):
		"display information for the target EMDirEntry"
		self.target=target
		
		# Set up image selectors for stacks
		self.wimnum.setRange(0,target.nimg)
		self.wimlist.clear()
		self.wimlist.addItems([str(i) for i in range(0,min(target.nimg,self.maxim))])
		if target.nimg>self.maxim: self.wimlist.addItem("...")
		self.wimnum.show()
		self.wimlist.show()

		self.wheadtree.clear()

	def imNumChange(self,num):
		"New image number"
		if num<500 : self.wimlist.setCurrentRow(num)
		else : self.showItem(num)
		
	def imSelChange(self):
		"New image selection"
		
		val=self.wimlist.currentItem().text()
		try:
			val=int(val)
			self.wimnum.setValue(val)
		except:
			print "Error with key :",val
			return
		
		self.showItem(val)
		
	def showItem(self,key):
		"""Shows header information for the selected item"""
		self.wheadtree.clear()
		try: trg=EMData(self.target.path(),key,True).get_attr_dict()		# read the header only, discard the emdata object
		except:
			print "Error reading: ",self.target.path(),key
			
		if trg==None :
			print "Warning: tried to read unavailable key: %s"%key
			#print self.bdb.keys()
			return
		
		self.addTreeItem(trg)
		
	def addTreeItem(self,trg,parent=None):
		"""(recursively) add an item to the tree"""
		itms=[]
		# Dictionaries may require recursion
		if isinstance(trg,dict):
			for k in trg.keys():
				itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((str(k),str(trg[k])))))
				if isinstance(trg[k],list) or isinstance(trg[k],tuple) or isinstance(trg[k],set) or isinstance(trg[k],dict):
					self.addTreeItem(trg[k],itms[-1])

		elif isinstance(trg,list) or isinstance(trg,tuple) or isinstance(trg,set):
			for k in trg:
				if isinstance(k,list) or isinstance(k,tuple) or isinstance(k,set) or isinstance(k,dict):
					try: itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((k.__class__.__name__,""))))
					except: itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList(("??",""))))
					self.addTreeItem(k,itms[-1])
				else:
					itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((str(k),""))))
			
		else:
			itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((str(trg),""))))

		if parent==None: 
			self.wheadtree.addTopLevelItems(itms)
			self.wheadtree.resizeColumnToContents(0)
		else : parent.addChildren(itms)

class EMInfoWin(QtGui.QWidget):
	"""The info window"""
	
	def __init__(self,parent=None):
		QtGui.QWidget.__init__(self,parent)
		
		self.target=None
		
		self.stack=QtGui.QStackedLayout(self)
		
		# We add one instance of 'infoPane' parent class to represent nothing
		
		self.stack.addWidget(EMInfoPane())
		
	def set_target(self,target,ftype):
		"""Display the info pane for target EMDirEntry with EMFileType instance ftype"""
		
		self.target=target
		self.ftype=ftype
		
		if target==None :
			self.stack.setCurrentIndex(0)
			return
		
		cls=ftype.infoClass()
		
		for i in range(self.stack.count()):
			if isinstance(self.stack.itemAt(i),cls) :
				self.stack.setCurrentIndex(i)
				pane=self.stack.itemAt(i)
				pane.display(target)
				break
		else:
			# If we got here, then we need to make a new instance of the appropriate pane
			if cls==None : print "No class ! (%s)"%str(ftype)
			pane=cls()
			i=self.stack.addWidget(pane)		# add the new pane and get its index
			pane.display(target)
			self.stack.setCurrentIndex(i)		# put the new pane on top
			

class EMBrowserWidget(QtGui.QWidget):
	"""This widget is a file browser for EMAN2. In addition to being a regular file browser, it supports:
	- getting information about recognized data types
	- embedding BDB: databases into the observed filesystem
	- remote database access (EMEN2)
	"""
	
	def __init__(self,parent=None,withmodal=False,multiselect=False,startpath=".",setsmode=None):
		# although this looks dumb it is necessary to break Python's issue with circular imports(a major weakness of Python IMO)
		global emscene3d, emdataitem3d
		import emscene3d
		import emdataitem3d 
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
		self.addBookmark("EMEN2","emen2:")
		self.wbookmarks.addSeparator()
		self.addBookmark("SSH","ssh:")
		self.wbookmarks.addSeparator()
		self.addBookmark("Root","/")
		self.addBookmark("Current",os.getcwd())
		self.addBookmark("Home",e2gethome())
		self.wbmfrbl.addWidget(self.wbookmarks)
		
		self.gbl.addWidget(self.wbookmarkfr,1,0)
		
		# This is the main window listing files and metadata
		self.wtree = QtGui.QTreeView()
		if multiselect: self.wtree.setSelectionMode(3)	# extended selection
		else : self.wtree.setSelectionMode(1)			# single selection
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
		QtCore.QObject.connect(self.wbutup, QtCore.SIGNAL('clicked(bool)'), self.buttonUp)
		QtCore.QObject.connect(self.wbutinfo, QtCore.SIGNAL('clicked(bool)'), self.buttonInfo)
		QtCore.QObject.connect(self.wtree, QtCore.SIGNAL('clicked(const QModelIndex)'), self.itemSel)
		QtCore.QObject.connect(self.wtree, QtCore.SIGNAL('activated(const QModelIndex)'), self.itemActivate)
		QtCore.QObject.connect(self.wtree, QtCore.SIGNAL('expanded(const QModelIndex)'), self.itemExpand)
		QtCore.QObject.connect(self.wpath, QtCore.SIGNAL('returnPressed()'), self.editPath)
		QtCore.QObject.connect(self.wbookmarks, QtCore.SIGNAL('actionTriggered(QAction*)'), self.bookmarkPress)

		self.setsmode=setsmode	# The sets mode is used when selecting bad particles 
		self.curmodel=None	# The current data model displayed in the tree
		self.curpath=None	# The path represented by the current data model
		self.curft=None		# a fileType instance for the currently hilighted object
		self.curactions=[]	# actions returned by the filtetype. Cached for speed
		self.models={}		# Cached models to avoid a lot of rereading (not sure if this is really worthwhile)
		self.pathstack=[]	# A stack of previous viewed paths
		self.infowin=None	# The 'info' window instance

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

		self.setPath(startpath)	# start in the local directory
		self.updthread.start()
		self.updtimer.start(300)

		self.result=None			# used in modal mode. Holds final selection

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
				if de.internalPointer().fillDetails(self.curmodel.getCacheDB()) : 
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

# 	def keyPressEvent(self,event):
# 		"""Make sure we update selection when keyboard is pressed"""
# 		print "key",event.__dict__
# 		QtGui.QTreeView.keyPressEvent(self.wtree,event)
# 		self.itemSel(None)

	def itemSel(self,qmi):
#		print "Item selected",qmi.row(),qmi.column(),qmi.internalPointer().path()
		qism=self.wtree.selectionModel().selectedRows()
		if len(qism)>1 : 
			self.wpath.setText("<multiple select>")
			if self.infowin!=None and not self.infowin.isHidden() :
				self.infowin.set_target(None)
			
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
				if self.setsmode: self.curft.setSetsDB(obj.path())	# If we want to enable bad particel picking
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
			# Bug fix by JFF (What if filetype is None????)
			else:
				self.curft=None
				self.curactions=[]
				for i,b in enumerate(self.wbutmisc):
					try:
						b.hide()
					except:
						pass
			
			if self.infowin!=None and not self.infowin.isHidden() :
				self.infowin.set_target(obj,ftc)
		
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
		"One of the programmable action buttons was pressed"
		
#		print "press ",self.curactions[num][0]
		
		self.curactions[num][2](self)				# This calls the action method

	def buttonOk(self,tog):
		"Button press"
		qism=self.wtree.selectionModel().selectedRows()
		self.result=[i.internalPointer().path().replace(os.getcwd(),".") for i in qism]
		self.emit(QtCore.SIGNAL("ok")) # this signal is important when e2ctf is being used by a program running its own eve
		
	def buttonCancel(self,tog):
		"Button press"
		self.result=[]
		self.emit(QtCore.SIGNAL("cancel")) # this signal is important when e2ctf is being used by a program running its own eve
		
	def buttonBack(self,tog):
		"Button press"
		try : 
			l=self.pathstack.index(self.curpath)
			if l==0 : raise Exception
		except : return
		
		self.setPath(self.pathstack[l-1],True)
		
	def buttonFwd(self,tog):
		"Button press"
		try : 
			l=self.pathstack.index(self.curpath)
			if l==len(self.pathstack)-1 : raise Exception
		except : return
		
		self.setPath(self.pathstack[l+1],True)
		
	def buttonUp(self,tog):
		"Button press"

		if "/" in self.curpath : newpath=self.curpath.rsplit("/",1)[0]
		else: newpath=os.path.realpath(self.curpath).rsplit("/",1)[0]

		print "Newpath: ",newpath
		if len(newpath)>1 : self.setPath(newpath)
		
	def buttonInfo(self,tog):
		if tog :
			if self.infowin==None :
				self.infowin=EMInfoWin()
				self.infowin.resize(500,600)
			self.infowin.show()
			self.infowin.raise_()
			qism=self.wtree.selectionModel().selectedRows()
			if len(qism)==1 : 
				self.infowin.set_target(qism[0].internalPointer(),self.curft)
			else : self.infowin.set_target(None,None)
		else :
			if self.infowin!=None:
				self.infowin.hide()

	def getResult(self):
		"""In modal mode, this will return the result of a selection. Returns None before
		ok/cancel have been pressed. If cancel is pressed or ok is pressed with no selection,
		returns an empty list []. With a valid selection, returns a list of path strings.
		When called after ok/cancel, also closes the dialog."""

		if self.result==None: return None
		
		self.close()
		return self.result

	def getCWD(self):
		""" In modal mode, this will return the directory the browser is in. This is useful for 
		using the browser to select a directory of interest. """
		
		# If a directory is selected, return this
		if self.result and os.path.isdir(self.result[0]):
			self.close()
			return self.result[0]
		# If there is no current path and no result dir
		if self.curpath==None: return None
		# Return current path
		self.close()
		return self.curpath
		
	def addBookmark(self,label,path):
		"""Add a new bookmark"""
		act=self.wbookmarks.addAction(label)
		act.setData(path)

	def setPath(self,path,silent=False,inimodel=EMFileItemModel):
		"""Sets the current root path for the browser window. If silent is true,
		the path history will not be updated."""
		
		self.updlist=[]
		
		self.curpath=str(path)
		self.wpath.setText(path)

		if path in self.models :
			self.curmodel=self.models[path]
		else : 
			self.curmodel=inimodel(path)
			self.models[self.curpath]=self.curmodel

		self.wtree.setSortingEnabled(False)
		self.wtree.setModel(self.curmodel)
		self.wtree.setSortingEnabled(True)
		self.wtree.resizeColumnToContents(0)
		self.wtree.resizeColumnToContents(1)
		self.wtree.resizeColumnToContents(3)
		self.wtree.resizeColumnToContents(4)

		self.expanded=set()
		# we add the child items to the list needing updates
		for i in xrange(self.curmodel.rowCount(None)-1,-1,-1):
			self.updlist.append(self.curmodel.index(i,0,None))

		if not silent:
			try: self.pathstack.remove(self.curpath)
			except: pass
			self.pathstack.append(self.curpath)

	def bookmarkPress(self,action):
		""
#		print "Got action ",action.text(),action.data().toString()
		
		self.setPath(action.data().toString())
#		self.wtree.setSelectionModel(myQItemSelection(self.curmodel))
	
	def closeEvent(self,event):
		self.updthreadexit=True
		for w in self.view2d+self.view2ds+self.view3d+self.viewplot2d+self.viewplot3d:
			w.close()

		if self.infowin!=None:
			self.infowin.close()

		event.accept()
		#self.app().close_specific(self)
		self.emit(QtCore.SIGNAL("module_closed")) # this signal is important when e2ctf is being used by a program running its own eve

# This is just for testing, of course
def test_result():
	global window
	print "Returned"
	print window.getResult()

if __name__ == '__main__':
	em_app = EMApp()
	window = EMBrowserWidget(withmodal=True,multiselect=True)
	QtCore.QObject.connect(window, QtCore.SIGNAL("ok"),test_result)
	QtCore.QObject.connect(window, QtCore.SIGNAL("cancel"),test_result)

	window.show()
	ret=em_app.exec_()
	try: window.updthreadexit=True
	except:pass
	sys.exit(ret)

	