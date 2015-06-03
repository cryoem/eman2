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

from EMAN2 import *
from EMAN2jsondb import js_open_dict
from emapplication import EMApp
from emimage2d import *
from emimagemx import *
from empdbitem3d import *
from emplot2d import *
from emplot3d import *
from expand_string import expand_string
from libpyUtils2 import EMUtil
from matching import matches_pats
import os
import re
from string import lower
import threading
import time
import traceback
from valslider import StringBox
import weakref

from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import Qt, QString, QChar


#---------------------------------------------------------------------------
def display_error(msg) :
	"""Displays an error message, in gui and on terminal."""

	print msg
	sys.stdout.flush()
	QtGui.QMessageBox.warning(None, "Error", msg)

# This is a floating point number-finding regular expression

renumfind = re.compile(r"-?[0-9]+\.*[0-9]*[eE]?[-+]?[0-9]*")

# We need to sort ints and floats as themselves, not string, John Flanagan

def safe_int(v) :
	"""Performs a safe conversion from a string to an int. If a non int is presented we return the lowest possible value"""

	try :
		return int(v)
	except (ValueError, TypeError) :
		return -sys.maxint-1

def safe_float(v) :
	"""Performs a safe conversion from a string to a float. If a non float is presented we return the lowest possible value"""

	try :
		return float(v)
	except :
		try : return float(v.split()[0])
		except :
			return sys.float_info.min

def isprint(s) :
	"""returns True if the string contains only printable ascii characters"""

	# Seems like no isprint() in python, this does basically the same thing

	mpd = s.translate("AAAAAAAAABBAABAAAAAAAAAAAAAAAAAA"+"B"*95+"A"*129)

	if "A" in mpd :
		ind = mpd.index("A")
#		print "bad chr %d at %d"%(ord(s[ind]), ind)
		return False

	return True

def askFileExists() :
	"""Opens a dialog and asks the user what to do if a file to be written to already exists"""

	box = QtGui.QMessageBox(4, "File Exists", "File already exists. What would you like to do ?")	# 4 means we're asking a question
	b1 = box.addButton("Append", QtGui.QMessageBox.AcceptRole)
	b2 = box.addButton("Overwrite", QtGui.QMessageBox.AcceptRole)
	b3 = box.addButton("Cancel", QtGui.QMessageBox.AcceptRole)

	box.exec_()

	if box.clickedButton() == b1 : return "append"
	elif box.clickedButton() == b2 : return "overwrite"
	else : return "cancel"

#---------------------------------------------------------------------------

class EMFileType(object) :
	"""This is a base class for handling interaction with files of different types. It includes a number of excution methods common to
	several different subclasses"""

	# A class dictionary keyed by EMDirEntry filetype string with value beign a single subclass of EMFileType. filetype strings are unique
	# When you add a new EMFiletype subclass, it must be added to this dictionary to be functional

	typesbyft = {}

	# a class dictionary keyed by file extension with values being either a single subclass or a list of possible subclasses for that extension
	# When you add a new EMFiletype subclass, it must be added to this dictionary to be functional

	typesbyext = {}

	# A list of all types that need to be checked when the file extension can't be interpreted

	alltocheck = ()

	setsmode = False

	def __init__(self, path) :
		if path[:2] == "./" : path = path[2:]
		self.path = path			# the current path this FileType is representing
		# self.setsdb = None		# sets associated with this file
		self.n = -1				# used to look at one image from a stack. -1 means not set.

	def setFile(self, path) :
		"""Represent a new file. Will update inspector if open. Assumes isValid already checked !"""

		if path[:2] == "./" : path = path[2:]
		self.path = path

	def setN(self, n) :
		"""Change the specific image methods should target"""

		self.n = n

	@staticmethod
	def name() :
		"""The unique name of this FileType. Stored in EMDirEntry.filetype for each file."""

		return None

	@staticmethod
	def isValid(path, header) :
		"""Returns (size, n, dim) if the referenced path is a file of this type, false if not valid. The first 4k block of data from the file is provided as well to avoid unnecesary file access."""

		return False

	@staticmethod
	def infoClass() :
		"""Returns a reference to the QWidget subclass for displaying information about this file"""

		return EMInfoPane

	#def setSetsDB(self, db_name) :
		#"Sets the emmxwidget to sets dbname"
		# self.setsdb = base_name(db_name)

	#def getSetsDB(self) :
		#"Returns the sets mode"
		#return self.setsdb

	def actions(self) :
		"""Returns a list of (name, help, callback) tuples detailing the operations the user can call on the current file.
		callbacks will also be passed a reference to the browser object."""

		return []

	def saveAs(self, brws) :
		"""Save an image file/stack to a new file"""

		outpath = QtGui.QInputDialog.getText(None, "Save Filename", "Filename to save to (type determined by extension)", 0, self.path)

		if outpath[1] != True : return

		outpath = str(outpath[0])

		if outpath == "" :		
			display_error ("No file name to save to entered.")
			return

		not_writeable_extensions = \
			[".dm2", ".dm3", ".dm4", ".v4l", ".sal", ".fts"] + \
			[".DM2", ".DM3", ".DM4", ".V4L", ".SAL", ".FTS"]

		for ext in not_writeable_extensions :
			if outpath.endswith(ext) :
				display_error ("Image file type '" + ext + "' is not writeable.")
				return

		if not EMUtil.is_valid_filename(outpath) :
			display_error ("Output file name '" + outpath + \
				"' does not have a valid image file extension.")
			return

		action = "overwrite"

		if file_exists(outpath) :
			action = askFileExists()
			if action == "cancel" : return

			if action == "overwrite" :
				remove_image(outpath)

		brws.busy()

		if self.n >= 0 : ns = [self.n]
		else : ns = xrange(EMUtil.get_image_count(self.path))

		for i in ns :
			im = EMData(self.path, i)
			im.write_image(outpath, -1)

		brws.notbusy()

	def plot2dApp(self, brws) :
		"""Append self to current plot"""

		brws.busy()

		if self.n >= 0 : data = EMData(self.path)
		else : data = EMData(self.path, self.n)

		try :
			target = brws.viewplot2d[-1]
			target.set_data(data, self.path.split("/")[-1].split("#")[-1])
		except :
			target = EMPlot2DWidget()
			brws.viewplot2d.append(target)
			target.set_data(data, self.path.split("/")[-1].split("#")[-1])

		target.qt_parent.setWindowTitle(self.path.split('/')[-1])

		brws.notbusy()
		target.show()
		target.raise_()

	def plot2dNew(self, brws) :
		"""Make a new plot"""

		brws.busy()

		if self.n >= 0 : data = EMData(self.path)
		else : data = EMData(self.path, self.n)

		target = EMPlot2DWidget()
		brws.viewplot2d.append(target)
		target.set_data(data, self.path.split("/")[-1].split("#")[-1])

		target.qt_parent.setWindowTitle(self.path.split('/')[-1])

		brws.notbusy()
		target.show()
		target.raise_()

	def show3dApp(self, brws) :
		"""Add to current 3-D window"""

		brws.busy()
		
		if self.n >= 0 : data = emdataitem3d.EMDataItem3D(self.path, n = self.n)
		else : data = emdataitem3d.EMDataItem3D(self.path)

		try :
			target = brws.view3d[-1]
		except :
			target = emscene3d.EMScene3D()
			brws.view3d.append(target)

		target.insertNewNode(self.path.split("/")[-1].split("#")[-1], data, parentnode = target)
		iso = emdataitem3d.EMIsosurface(data)
		target.insertNewNode('Isosurface', iso, parentnode = data)
		target.initialViewportDims(data.getData().get_xsize())	# Scale viewport to object size
		target.setCurrentSelection(iso)				# Set isosurface to display upon inspector loading
		target.updateSG()	# this is needed because this might just be an addition to the SG rather than initialization

		target.setWindowTitle(self.path.split('/')[-1])

		brws.notbusy()
		target.show()
		target.raise_()

	def show3DNew(self, brws) :
		"""New 3-D window"""

		brws.busy()

		if self.n >= 0 : data = emdataitem3d.EMDataItem3D(self.path, n = self.n)
		else : data = emdataitem3d.EMDataItem3D(self.path)

		target = emscene3d.EMScene3D()
		brws.view3d.append(target)

		target.insertNewNode(self.path.split("/")[-1].split("#")[-1], data)
		iso = emdataitem3d.EMIsosurface(data)
		target.insertNewNode('Isosurface', iso, parentnode = data)
		target.initialViewportDims(data.getData().get_xsize())	# Scale viewport to object size
		target.setCurrentSelection(iso)				# Set isosurface to display upon inspector loading
		brws.notbusy()
		target.setWindowTitle(self.path.split('/')[-1])

		target.show()
		target.raise_()

	def show3DAll(self, brws) :
		"""All in new 3-D window (3-D stacks)"""

		brws.busy()

		target = emscene3d.EMScene3D()
		brws.view3d.append(target)

		for n in xrange(self.nimg) :
			data = emdataitem3d.EMDataItem3D(self.path, n = n)
			target.insertNewNode("{} #{}".format(self.path.split("/")[-1], n), data)
			iso = emdataitem3d.EMIsosurface(data)
			target.insertNewNode('Isosurface', iso, parentnode = data)

		target.initialViewportDims(data.getData().get_xsize())	# Scale viewport to object size
		target.setCurrentSelection(iso)				# Set isosurface to display upon inspector loading
		brws.notbusy()
		target.setWindowTitle(self.path.split('/')[-1])

		target.show()
		target.raise_()

	def show2dStack(self, brws) :
		"""A set of 2-D images together in an existing window"""
	
		brws.busy()

		# if self.dim[2] > 1 :
			# data = []
			# for z in range(self.dim[2]) :
				# data.append(EMData(self.path, 0, False, Region(0, 0, z, self.dim[0], self.dim[1], 1)))
		# else : data = EMData.read_images(self.path)

		try :
			target = brws.view2ds[-1]
			target.set_data(self.path, self.path)
			#if self.getSetsDB() : target.set_single_active_set(self.getSetsDB())
		except :
			target = EMImageMXWidget()
			target.set_data(self.path, self.path)
			QtCore.QObject.connect(target, QtCore.SIGNAL("mx_image_double"), target.mouse_double_click)		# this makes class average viewing work in app mode
			# if self.getSetsDB() : target.set_single_active_set(self.getSetsDB())
			brws.view2ds.append(target)

		target.qt_parent.setWindowTitle("Stack - "+self.path.split('/')[-1])

		brws.notbusy()
		target.show()
		target.raise_()

	def show2dStackNew(self, brws) :
		"""A set of 2-D images together in a new window"""

		brws.busy()

		# if self.dim[2] > 1 :
			# data = []
			# for z in range(self.dim[2]) :
				# data.append(EMData(self.path, 0, False, Region(0, 0, z, self.dim[0], self.dim[1], 1)))
		# else : data = EMData.read_images(self.path)

		target = EMImageMXWidget()
		target.set_data(self.path, self.path)
		QtCore.QObject.connect(target, QtCore.SIGNAL("mx_image_double"), target.mouse_double_click)
		# if self.getSetsDB() : target.set_single_active_set(self.getSetsDB())
		brws.view2ds.append(target)

		target.qt_parent.setWindowTitle("Stack - "+self.path.split('/')[-1])

		brws.notbusy()
		target.show()
		target.raise_()

	def show2dSingle(self, brws) :
		"""Show a single 2-D image"""

		brws.busy()

		if self.nimg > 1 :
			if self.n >= 0 : data = EMData(self.path, self.n)
			else : data = EMData.read_images(self.path)
		else : data = EMData(self.path)

		try :
			target = brws.view2d[-1]
			target.set_data(data)
		except :
			target = EMImage2DWidget(data)
			brws.view2d.append(target)

		target.setWindowTitle(self.path.split('/')[-1])

		brws.notbusy()
		target.show()
		target.raise_()

	def show2dSingleNew(self, brws) :
		"""Show a single 2-D image"""

		brws.busy()

		if self.nimg > 1 :
			if self.n >= 0 : data = EMData(self.path, self.n)
			else : data = EMData.read_images(self.path)
		else : data = EMData(self.path)

		target = EMImage2DWidget(data)
		brws.view2d.append(target)

		target.qt_parent.setWindowTitle(self.path.split('/')[-1])

		brws.notbusy()
		target.show()
		target.raise_()

	def showFilterTool(self, brws) :
		"""Open in e2filtertool.py"""

		os.system("e2filtertool.py %s &"%self.path)

#---------------------------------------------------------------------------

class EMTextFileType(EMFileType) :
	"""FileType for files containing normal ASCII text"""

	@staticmethod
	def name() :
		"""The unique name of this FileType. Stored in EMDirEntry.filetype for each file."""

		return "Text"

	@staticmethod
	def isValid(path, header) :
		"""Returns (size, n, dim) if the referenced path is a file of this type, None if not valid. The first 4k block of data from the file is provided as well to avoid unnecesary file access."""

		if not isprint(header) : return False			# demand printable Ascii. FIXME: what about unicode ?

		try : size = os.stat(path)[6]
		except : return False

		if size > 5000000 : dim = "big"
		else :
			f = file(path, "r").read()
			lns = max(f.count("\n"), f.count("\r"))
			dim = "%d ln"%lns

		return (size, "-", dim)

	@staticmethod
	def infoClass() :
		"""Returns a reference to the QWidget subclass for displaying information about this file"""

		return EMTextInfoPane

	def actions(self) :
		"""No actions other than the inspector for text files"""

		return []

class EMHTMLFileType(EMFileType) :
	"""FileType for files containing HTML text"""

	@staticmethod
	def name() :
		"""The unique name of this FileType. Stored in EMDirEntry.filetype for each file."""

		return "HTML"

	@staticmethod
	def isValid(path, header) :
		"""Returns (size, n, dim) if the referenced path is a file of this type, None if not valid. The first 4k block of data from the file is provided as well to avoid unnecesary file access."""

		if not isprint(header) : return False			# demand printable Ascii. FIXME: what about unicode ?
		if not "<html>" in header.lower() : return False # For the moment, we demand an <html> tag somewhere in the first 4k

		try :
			size = os.stat(path)[6]
		except : return False

		if size > 5000000 : dim = "big"
		else :
			f = file(path, "r").read()
			lns = max(f.count("\n"), f.count("\r"))
			dim = "%d ln"%lns

		return (size, "-", dim)

	@staticmethod
	def infoClass() :
		"""Returns a reference to the QWidget subclass for displaying information about this file"""

		return EMHTMLInfoPane

	def actions(self) :
		"""No actions other than the inspector for HTML files"""

		return [("Firefox", "Open in Firefox", self.showFirefox)]

	def showFirefox(self, brws) :
		"""Try to open file in firefox"""

		if get_platform() == "Linux" :
			os.system("firefox -new-tab file://%s"%os.path.abspath(self.path))
		elif get_platform() == "Darwin" :
			os.system("open {}".format(os.path.abspath(self.path)))		# uses the default browser
		else : 
			print "Sorry, I don't know how to run Firefox on this platform"

#---------------------------------------------------------------------------

class EMPlotFileType(EMFileType) :
	"""FileType for files containing normal ASCII text"""

	@staticmethod
	def name() :
		"""The unique name of this FileType. Stored in EMDirEntry.filetype for each file."""

		return "Plot"

	@staticmethod
	def isValid(path, header) :
		"""Returns (size, n, dim) if the referenced path is a file of this type, None if not valid. The first 4k block of data from the file is provided as well to avoid unnecesary file access."""

		if not isprint(header) : return False

		# We need to try to count the columns in the file

		hdr = header.splitlines()
		numc = 0

		for l in hdr :
			if l[0] == "#" : continue		# comment lines ok

			try : numc = len([float(i) for i in renumfind.findall(l)])		# number of numeric columns
			except :
				return False		# shouldn't really happen...

			if numc > 0 : break			# just finding the number of columns

		# If we couldn't find at least one valid line with some numbers, give up

		if numc == 0 : return False

		try : size = os.stat(path)[6]
		except : return False

		# Make sure all of the lines have the same number of columns

		fin = file(path, "r")
		numr = 0

		for l in fin :
			if l[0] == "#" or len(l) < 2 or "nan" in l : continue

			lnumc = len([float(i) for i in renumfind.findall(l)])
			if lnumc != 0 and lnumc != numc : return False				# 0 means the line contains no numbers, we'll live with that, but if there are numbers, it needs to match
			if lnumc != 0 : numr += 1

		return (size, "-", "%d x %d"%(numr, numc))

	@staticmethod
	def infoClass() :
		"""Returns a reference to the QWidget subclass for displaying information about this file"""

		return EMPlotInfoPane

	def __init__(self, path) :
		if path[:2] == "./" : path = path[2:]
		EMFileType.__init__(self, path)	# the current path this FileType is representing

		# Make sure all of the lines have the same number of columns

		fin = file(path, "r")
		numr = 0
		numc = 0

		for l in fin :
			if "nan" in l :
				print "Warning, NaN present in file"
				continue

			if l[0] == "#" : continue

			lnumc = len([float(i) for i in renumfind.findall(l)])

			if lnumc != 0 and numc == 0 : numc = lnumc
			elif lnumc != 0 and lnumc != numc :
				print "Error: invalid Plot file:", path
				self.numr = 0
				self.numc = 0
				return
			elif lnumc != 0 : numr += 1

			self.numr = numr
			self.numc = numc

	def actions(self) :
		"""Returns a list of (name, help, callback) tuples detailing the operations the user can call on the current file.
		callbacks will also be passed a reference to the browser object."""

		if self.numc > 2 : return [("Plot 2D", "Add to current plot", self.plot2dApp), ("Plot 2D+", "Make new plot", self.plot2dNew), 
			("Plot 3D+", "Make new 3-D plot", self.plot3dNew), ("Plot 3D", "Add to current 3-D plot", self.plot3dApp)]

		return [("Plot 2D", "Add to current plot", self.plot2dApp), ("Plot 2D+", "Make new plot", self.plot2dNew)]

	def plot2dApp(self, brws) :
		"""Append self to current plot"""

		brws.busy()

		data1 = []
		fin = file(self.path, "r")
		numr = 0

		for l in fin :
			if l[0] == "#" or "nan" in l : continue
			data1.append([float(i) for i in renumfind.findall(l)])

		data = []

		for c in xrange(self.numc) :
			data.append([i[c] for i in data1])

		try :
			target = brws.viewplot2d[-1]
			target.set_data(data, remove_directories_from_name(self.path, 1))
		except :
			target = EMPlot2DWidget()
			brws.viewplot2d.append(target)
			target.set_data(data, remove_directories_from_name(self.path, 1))

		target.qt_parent.setWindowTitle(self.path.split('/')[-1])

		brws.notbusy()
		target.show()
		target.raise_()

	def plot2dNew(self, brws) :
		"""Make a new plot"""

		brws.busy()

		data1 = []
		fin = file(self.path, "r")
		numr = 0

		for l in fin :
			if l[0] == "#" or "nan" in l : continue
			data1.append([float(i) for i in renumfind.findall(l)])

		data = []

		for c in xrange(self.numc) :
			data.append([i[c] for i in data1])

		target = EMPlot2DWidget()
		brws.viewplot2d.append(target)
		target.set_data(data, remove_directories_from_name(self.path, 1))

		target.qt_parent.setWindowTitle(self.path.split('/')[-1])

		brws.notbusy()
		target.show()
		target.raise_()

	def plot3dApp(self, brws) :
		"""Append self to current 3-D plot"""

	def plot3dNew(self, brws) :
		"""Make a new 3-D plot"""

#---------------------------------------------------------------------------

class EMFolderFileType(EMFileType) :
	"""FileType for Folders"""

	@staticmethod
	def name() :
		"""The unique name of this FileType. Stored in EMDirEntry.filetype for each file."""

		return "Folder"
	@staticmethod
	def isValid(path, header) :
		"""Returns (size, n, dim) if the referenced path is a file of this type, None if not valid. The first 4k block of data from the file is provided as well to avoid unnecesary file access."""

		return False

	@staticmethod
	def infoClass() :
		"""Returns a reference to the QWidget subclass for displaying information about this file"""

		return EMFolderInfoPane

	def actions(self) :
		"""Returns a list of (name, callback) tuples detailing the operations the user can call on the current file"""

		return []

#---------------------------------------------------------------------------

class EMJSONFileType(EMFileType) :
	"""FileType for JSON files"""

	@staticmethod
	def name() :
		"""The unique name of this FileType. Stored in EMDirEntry.filetype for each file."""

		return "JSON"

	@staticmethod
	def isValid(path, header) :
		"""Returns (size, n, dim) if the referenced path is a file of this type, None if not valid. The first 4k block of data from the file is provided as well to avoid unnecesary file access."""

		if path[-5:] == ".json" and header.strip()[0] == "{" : return (humansize(os.stat(path).st_size), "-", "-")
		else : return None
			# sz = len(js_open_dict(path).keys())
			# return ("-", sz, "-")

	@staticmethod
	def infoClass() :
		"""Returns a reference to the QWidget subclass for displaying information about this file"""

		return EMJSONInfoPane

	def __init__(self, path) :
		if path[:2] == "./" : path = path[2:]

		EMFileType.__init__(self, path)	# the current path this FileType is representing

		self.js = js_open_dict(path)
		self.keys = self.js.keys()
		self.dim = (0, 0, 0)

	def __del__(self) :
		try : self.js.close()
		except : pass

	def actions(self) :
		"""Returns a list of (name, help, callback) tuples detailing the operations the user can call on the current file"""

		# No actions for now...

		return []

#---------------------------------------------------------------------------

class EMBdbFileType(EMFileType) :
	"""FileType for Folders"""

	@staticmethod
	def name() :
		"""The unique name of this FileType. Stored in EMDirEntry.filetype for each file."""

		return "BDB"

	@staticmethod
	def isValid(path, header) :
		"""Returns (size, n, dim) if the referenced path is a file of this type, None if not valid. The first 4k block of data from the file is provided as well to avoid unnecesary file access."""

		return False

	@staticmethod
	def infoClass() :
		"""Returns a reference to the QWidget subclass for displaying information about this file"""

		return EMBDBInfoPane

	def __init__(self, path) :
		if path[:2] == "./" : path = path[2:]

		EMFileType.__init__(self, path)	# the current path this FileType is representing
		self.bdb = db_open_dict(path, ro = True)

		# here we assume the bdb either contains numbered images OR key/value pairs. Not strictly valid,
		# but typically true

		self.nimg = len(self.bdb)

		if self.nimg == 0 : self.keys = self.bdb.keys()
		else : self.keys = None

		if self.nimg > 0 :
			im0 = EMData(path, 0, True)
			self.dim = (im0["nx"], im0["ny"], im0["nz"])
		else : self.dim = (0, 0, 0)

	def actions(self) :
		"""Returns a list of (name, help, callback) tuples detailing the operations the user can call on the current file"""

		# single 3-D
		if self.nimg == 1 and self.dim[2] > 1 :
			return [("Show 3D", "Add to 3D window", self.show3dApp), ("Show 3D+", "New 3D Window", self.show3DNew), ("Show Stack", "Show as set of 2-D Z slices", self.show2dStack), 
				("Show Stack+", "Show all images together in a new window", self.show2dStackNew), ("Show 2D", "Show in a scrollable 2D image window", self.show2dSingle), 
				("Show 2D+", "Show all images, one at a time in a new window", self.show2dSingleNew), ("Chimera", "Open in chimera (if installed)", self.showChimera), 
				("FilterTool", "Open in e2filtertool.py", self.showFilterTool), ("ProjXYZ", "Make projections along Z,Y,X", self.showProjXYZ ),("Save As", "Saves images in new file format", self.saveAs)]
		# single 2-D
		elif self.nimg == 1 and self.dim[1] > 1 :
			return [("Show 2D", "Show in a 2D single image display", self.show2dSingle), ("Show 2D+", "Show in new 2D single image display", self.show2dSingleNew), ("FilterTool", "Open in e2filtertool.py", self.showFilterTool), ("Save As", "Saves images in new file format", self.saveAs)]
		# single 1-D
		elif self.nimg == 1 :
			return [("Plot 2D", "Add to current plot", self.plot2dApp), ("Plot 2D+", "Make new plot", self.plot2dNew), 
				("Show 2D", "Replace in 2D single image display", self.show2dSingle), ("Show 2D+", "New 2D single image display", self.show2dSingleNew), ("Save As", "Saves images in new file format", self.saveAs)]
		# 3-D stack
		elif self.nimg > 1 and self.dim[2] > 1 :
			return [("Show 3D", "Show all in a single 3D window", self.show3DNew), ("Chimera", "Open in chimera (if installed)", self.showChimera), ("Save As", "Saves images in new file format", self.saveAs)]
		# 2-D stack
		elif self.nimg > 1 and self.dim[1] > 1 :
			return [("Show Stack", "Show all images together in one window", self.show2dStack), ("Show Stack+", "Show all images together in a new window", self.show2dStackNew), 
				("Show 2D", "Show all images, one at a time in current window", self.show2dSingle), ("Show 2D+", "Show all images, one at a time in a new window", self.show2dSingleNew), ("FilterTool", "Open in e2filtertool.py", self.showFilterTool), ("Save As", "Saves images in new file format", self.saveAs)]
		# 1-D stack
		elif self.nimg > 0 :
			return [("Plot 2D", "Plot all on a single 2-D plot", self.plot2dNew), ("Save As", "Saves images in new file format", self.saveAs)]

		return []

	def showProjXYZ(self,brws) :
		"""Show XYZ projections of 3-D volume"""

		brws.busy()

		tmp=EMData(self.path, self.n)
		data=[tmp.process("misc.directional_sum",{"axis":axis}) for axis in "zyx"]
	
		target = EMImage2DWidget(data)
		brws.view2d.append(target)

		target.qt_parent.setWindowTitle(self.path.split('/')[-1])

		brws.notbusy()
		target.show()
		target.raise_()

	def showChimera(self, brws) :
		"""Open in Chimera"""

		if get_platform() == "Linux" :
			os.system("e2proc3d.py %s /tmp/vol.hdf"%self.path)		# Probably not a good hack to use, but it will do for now...
			os.system("chimera /tmp/vol.hdf&")
		elif get_platform() == "Darwin" :
			os.system("e2proc3d.py %s /tmp/vol.hdf"%self.path)		# Probably not a good hack to use, but it will do for now...
			os.system("/Applications/Chimera.app/Contents/MacOS/chimera /tmp/vol.hdf&")
		else : print "Sorry, I don't know how to run Chimera on this platform"

#---------------------------------------------------------------------------

class EMImageFileType(EMFileType) :
	"""FileType for files containing a single 2-D image"""

	def __init__(self, path) :
		if path[:2] == "./" : path = path[2:]

		EMFileType.__init__(self, path)	# the current path this FileType is representing
		self.nimg = EMUtil.get_image_count(path)
		im0 = EMData(path, 0, True)
		self.dim = (im0["nx"], im0["ny"], im0["nz"])

	def closeEvent(self, event) :
#		E2saveappwin("e2display", "main", self)

		self.updthreadexit = True

		for w in self.view2d+self.view2ds+self.view3d+self.viewplot2d+self.viewplot3d :
			w.close()

		if self.infowin != None :
			self.infowin.close()

		event.accept()
		# self.app().close_specific(self)
		self.emit(QtCore.SIGNAL("module_closed"))

	@staticmethod
	def name() :
		"""The unique name of this FileType. Stored in EMDirEntry.filetype for each file."""

		return "Image"

	@staticmethod
	def isValid(path, header) :
		"""Returns (size, n, dim) if the referenced path is a file of this type, None if not valid. The first 4k block of data from the file is provided as well to avoid unnecesary file access."""

		return False

	@staticmethod
	def infoClass() :
		"""Returns a reference to the QWidget subclass for displaying information about this file"""

		return EMImageInfoPane

	def actions(self) :
		"""Returns a list of (name, callback) tuples detailing the operations the user can call on the current file"""

		# single 3-D
		if  self.dim[2] > 1 :
			return [("Show 3D", "Add to 3D window", self.show3dApp), ("Show 3D+", "New 3D Window", self.show3DNew), ("Show Stack", "Show as set of 2-D Z slices", self.show2dStack), 
				("Show Stack+", "Show all images together in a new window", self.show2dStackNew), ("Show 2D", "Show in a scrollable 2D image window", self.show2dSingle), 
				("Show 2D+", "Show all images, one at a time in a new window", self.show2dSingleNew), ("Chimera", "Open in chimera (if installed)", self.showChimera), 
				("FilterTool", "Open in e2filtertool.py", self.showFilterTool), ("ProjXYZ", "Make projections along Z,Y,X", self.showProjXYZ ), ("Save As", "Saves images in new file format", self.saveAs)]
		## 2-D stack, STEVE: THIS SHOULD NOT BE HERE
		# elif self.nimg > 1 :
			# return [("Show Stack", "Show as set of 2-D Z slices", self.show2dStack), ("Show Stack+", "Show all images together in a new window", self.show2dStackNew), ("Show 2D", "Show in a scrollable 2D image window", self.show2dSingle), 
				# ("Show 2D+", "Show all images, one at a time in a new window", self.show2dSingleNew), ("FilterTool", "Open in e2filtertool.py", self.showFilterTool), ("Save As", "Saves images in new file format", self.saveAs)]
		elif  self.dim[1] > 1 :
			return [("Show 2D", "Show in a 2D single image display", self.show2dSingle), ("Show 2D+", "Show in new 2D single image display", self.show2dSingleNew), ("FilterTool", "Open in e2filtertool.py", self.showFilterTool), ("Save As", "Saves images in new file format", self.saveAs)]
		# single 1-D
		else :
			return [("Plot 2D", "Add to current plot", self.plot2dApp), ("Plot 2D+", "Make new plot", self.plot2dNew), 
				("Show 2D", "Replace in 2D single image display", self.show2dSingle), ("Show 2D+", "New 2D single image display", self.show2dSingleNew), ("Save As", "Saves images in new file format", self.saveAs)]

	def showProjXYZ(self,brws) :
		"""Show XYZ projections of 3-D volume"""

		brws.busy()

		if self.n>=0 : tmp=EMData(self.path, self.n)
		else: tmp=EMData(self.path, 0)
		data=[tmp.process("misc.directional_sum",{"axis":axis}) for axis in "zyx"]
	
		target = EMImage2DWidget(data)
		brws.view2d.append(target)

		target.qt_parent.setWindowTitle(self.path.split('/')[-1])

		brws.notbusy()
		target.show()
		target.raise_()

	def showChimera(self, brws) :
		"""Open in Chimera"""

		if get_platform() == "Linux" :
			# these types are supported natively in Chimera
			if EMUtil.get_image_type(self.path) in (IMAGE_HDF, IMAGE_MRC, IMAGE_SPIDER, IMAGE_SINGLE_SPIDER) :
				os.system("chimera %s &"%self.path)
			else :
				os.system("e2proc3d.py %s /tmp/vol.hdf"%self.path)		# Probably not a good hack to use, but it will do for now...
				os.system("chimera /tmp/vol.hdf&")
		elif get_platform() == "Darwin" :
			os.system("e2proc3d.py %s /tmp/vol.hdf"%self.path)		# Probably not a good hack to use, but it will do for now...
			os.system("/Applications/Chimera.app/Contents/MacOS/chimera /tmp/vol.hdf&")
		else : print "Sorry, I don't know how to run Chimera on this platform"

#---------------------------------------------------------------------------

class EMStackFileType(EMFileType) :
	"""FileType for files containing a set of 1-3D images"""

	@staticmethod
	def name() :
		"""The unique name of this FileType. Stored in EMDirEntry.filetype for each file."""""

		return "Image Stack"

	@staticmethod
	def isValid(path, header) :
		"""Returns (size, n, dim) if the referenced path is a file of this type, None if not valid. The first 4k block of data from the file is provided as well to avoid unnecesary file access."""

		return False

	@staticmethod
	def infoClass() :
		"""Returns a reference to the QWidget subclass for displaying information about this file"""

		return EMStackInfoPane

	def __init__(self, path) :
		if path[:2] == "./" : path = path[2:]

		EMFileType.__init__(self, path)	# the current path this FileType is representing
		self.nimg = EMUtil.get_image_count(path)

		try : im0 = EMData(path, 0, True)
		except :
			for i in xrange(1, 10) :
				try : im0 = EMData(path, i, True)
				except : continue
				break

		try : self.dim = (im0["nx"], im0["ny"], im0["nz"])
		except :
			print "First 10 images all missing in ", path
			self.dim = "?"

	def actions(self) :
		"""Returns a list of (name, callback) tuples detailing the operations the user can call on the current file"""
		# 3-D stack
		if self.nimg > 1 and self.dim[2] > 1:
			return [("Show all 3D", "Show all in a single 3D window", self.show3DAll), ("Show 1st 3D", "Show only the first volume", self.show3DNew), ("Chimera", "Open in chimera (if installed)", self.showChimera), ("Save As", "Saves images in new file format", self.saveAs)]
		# 2-D stack
		elif self.nimg > 1 and self.dim[1] > 1:
			return [("Show Stack", "Show all images together in one window", self.show2dStack), ("Show Stack+", "Show all images together in a new window", self.show2dStackNew), 
				("Show 2D", "Show all images, one at a time in current window", self.show2dSingle), ("Show 2D+", "Show all images, one at a time in a new window", self.show2dSingleNew), ("FilterTool", "Open in e2filtertool.py", self.showFilterTool), ("Save As", "Saves images in new file format", self.saveAs)]
		# 1-D stack
		elif self.nimg > 1:
			return [("Plot 2D", "Plot all on a single 2-D plot", self.plot2dNew), ("Save As", "Saves images in new file format", self.saveAs)]
		else : print "Error: stackfile with < 2 images ? (%s)"%self.path

		return []

	def showChimera(self, brws):
		"""Open in Chimera"""
		if get_platform() == "Linux":
			# these types are supported natively in Chimera
			if EMUtil.get_image_type("tst.hdf") in (IMAGE_HDF, IMAGE_MRC, IMAGE_SPIDER, IMAGE_SINGLE_SPIDER) :
				os.system("chimera %s &"%self.path)
			else:
				os.system("e2proc3d.py %s /tmp/vol.hdf"%self.path)		# Probably not a good hack to use, but it will do for now...
				os.system("chimera /tmp/vol.hdf&")
		elif get_platform() == "Darwin" :
			os.system("e2proc3d.py %s /tmp/vol.hdf"%self.path)		# Probably not a good hack to use, but it will do for now...
			os.system("/Applications/Chimera.app/Contents/MacOS/chimera /tmp/vol.hdf&")
		else : print "Sorry, I don't know how to run Chimera on this platform"

#---------------------------------------------------------------------------

class EMPDBFileType(EMFileType):
	
	"""FileType for files with original pdb format"""
	
	def __init__(self, path) :
		if path[:2] == "./" : path = path[2:]
		EMFileType.__init__(self, path)	# the current path this FileType is representing
	
	@staticmethod
	def name():
		"""
		The unique name of this FileType. Stored in EMDirEntry.filetype for each file.
		"""
		return "PDB"

	@staticmethod
	def infoClass():
		"""
		Returns a reference to the QWidget subclass for displaying information about this file
		"""
		return EMPDBInfoPane

	@staticmethod
	def isValid(path, header) :
		"""
		Returns (size, n, dim) if the referenced path is a file of this type, None if not valid. 
		The first 4k block of data from the file is provided as well to avoid unnecesary file access.
		"""
		proper_exts = ['pdb','ent']
		ext = os.path.basename(path).split('.')[-1]
		if ext not in proper_exts: return False
		
		if not isprint(header) : return False			# demand printable Ascii. FIXME: what about unicode ?
		
		try : size = os.stat(path)[6]
		except : return False
		
		if size > 5000000: dim = "big"
		else :
			f = file(path, "r").read()
			lns = max(f.count("\n"), f.count("\r"))
			dim = "%d ln"%lns
		
		return (size, "-", dim)

	def actions(self):
		"""
		Returns a list of (name, callback) tuples detailing the operations the user can call on the current file
		"""
		return [("Show Ball and Stick", "Show ball and stick representation of this PDB model in a new 3D window", self.showBallStick3DNew), ("Show Ball and Stick +", "Show ball and stick representation of this PDB model in the current 3D window", self.showBallStick3dApp), ("Show Spheres", "Show spheres representation of this PDB model in a new 3D window", self.showSpheres3DNew), ("Show Spheres +", "Show spheres representation of this PDB model in the current 3D window", self.showSpheres3dApp), ("Chimera", "Open this PDB file in chimera (if installed)", self.showChimera), ("Save As", "Saves a copy of the selected PDB file", self.saveAs)]

	def showSpheres3DApp(self, brws):
		"""New 3-D window"""
		brws.busy()
		pdb_model = EMPDBItem3D(self.path)
		target = emscene3d.EMScene3D()
		brws.view3d.append(target)
		target.insertNewNode(self.path.split("/")[-1].split("#")[-1],pdb_model)
		modeltype = EMSphereModel(self.path)
		target.insertNewNode(modeltype.representation, modeltype, parentnode = pdb_model)
		target.initialViewportDims(pdb_model.getBoundingBoxDimensions()[0])	# Scale viewport to object size
		target.setCurrentSelection(modeltype)	# Set style to display upon inspector loading
		brws.notbusy()
		target.setWindowTitle(pdb_model.getName())
		target.show()
		target.raise_()

	def showSpheres3dNew(self, brws):
		"""Add to current 3-D window"""
		brws.busy()
		pdb_model = EMPDBItem3D(self.path)
		try: target = brws.view3d[-1]
		except:
			target = emscene3d.EMScene3D()
			brws.view3d.append(target)
		target.insertNewNode(self.path.split("/")[-1].split("#")[-1],pdb_model, parentnode = target)
		modeltype = EMSphereModel(self.path)
		target.insertNewNode(modeltype.representation, modeltype, parentnode = pdb_model)
		target.initialViewportDims(pdb_model.getBoundingBoxDimensions()[0])	# Scale viewport to object size
		target.setCurrentSelection(modeltype)	# Set style to display upon inspector loading
		#target.updateSG()	# this is needed because this might just be an addition to the SG rather than initialization
		target.setWindowTitle(pdb_model.getName())
		brws.notbusy()
		target.show()
		target.raise_()

	def showBallStick3DApp(self, brws):
		"""New 3-D window"""
		brws.busy()
		pdb_model = EMPDBItem3D(self.path)
		target = emscene3d.EMScene3D()
		brws.view3d.append(target)
		target.insertNewNode(self.path.split("/")[-1].split("#")[-1],pdb_model)
		modeltype = EMBallStickModel(self.path) #parent=pdb_model)
		target.insertNewNode(modeltype.representation, modeltype, parentnode = pdb_model)
		target.initialViewportDims(pdb_model.getBoundingBoxDimensions()[0])	# Scale viewport to object size
		target.setCurrentSelection(modeltype)	# Set style to display upon inspector loading
		brws.notbusy()
		target.setWindowTitle(pdb_model.getName())
		target.show()
		target.raise_()

	def showBallStick3dNew(self, brws):
		"""Add to current 3-D window"""
		brws.busy()
		pdb_model = EMPDBItem3D(self.path)
		try: target = brws.view3d[-1]
		except:
			target = emscene3d.EMScene3D()
			brws.view3d.append(target)
		target.insertNewNode(self.path.split("/")[-1].split("#")[-1],pdb_model, parentnode = target)
		modeltype = EMBallStickModel(self.path)#parent=pdb_model)
		target.insertNewNode(modeltype.representation, modeltype, parentnode = pdb_model)
		target.initialViewportDims(pdb_model.getBoundingBoxDimensions()[0])	# Scale viewport to object size
		target.setCurrentSelection(modeltype)	# Set style to display upon inspector loading
		#target.updateSG()	# this is needed because this might just be an addition to the SG rather than initialization
		target.setWindowTitle(pdb_model.getName())
		brws.notbusy()
		target.show()
		target.raise_()

	def showChimera(self):
		"""Open in Chimera"""
		if get_platform() == "Linux":
			os.system("chimera %s &" % self.path)
		elif get_platform() == "Darwin":
			os.system("/Applications/Chimera.app/Contents/MacOS/chimera %s &" % self.path)
		else:
			print "Sorry, I don't know how to run Chimera on this platform"

#---------------------------------------------------------------------------

# These are set all together at the end rather than after each class for efficiency

EMFileType.typesbyft = {
	"Folder"      : EMFolderFileType,
	"JSON"        : EMJSONFileType,
	"BDB"         : EMBdbFileType,
	"PDB"         : EMPDBFileType,
	"Text"        : EMTextFileType,
	"Plot"        : EMPlotFileType,
	"Image"       : EMImageFileType,
	"Image Stack" : EMStackFileType,
	"HTML"        : EMHTMLFileType
}

# Note that image types are not included here, and are handled via a separate mechanism
# note that order is important in the tuples. The most specific filetype should go first, and
# the most general last (they will be checked in order)

EMFileType.extbyft = {
	".json" : (EMJSONFileType, EMTextFileType),
	".pdb"  : (EMPDBFileType,  EMTextFileType),
	".ent"  : (EMPDBFileType,  EMTextFileType),
	".txt"  : (EMPlotFileType, EMTextFileType),
	".htm"  : (EMHTMLFileType, EMTextFileType),
	".html" : (EMHTMLFileType, EMTextFileType)
}

# Default Choices when extension doesn't work
# We don't need to test for things like Images because they are fully tested outside this mechanism

EMFileType.alltocheck = (EMPlotFileType, EMPDBFileType, EMTextFileType)

#---------------------------------------------------------------------------

class EMDirEntry(object) :
	"""Represents a directory entry in the filesystem"""

	# list of lambda functions to extract column values for sorting
	col = (lambda x:int(x.index), lambda x:x.name, lambda x:x.filetype, lambda x:x.size, lambda x:x.dim, lambda x:safe_int(x.nimg), lambda x:x.date)
#	classcount = 0

	def __init__(self, root, name, index, parent = None, hidedot = True, dirregex = None) :
		"""The path for this item is root/name.
		Parent (EMDirEntry) must be specified if it exists.
		hidedot will cause hidden files (starting with .) to be excluded"""

		self.__parent = parent	# single parent
		self.__children = None	# ordered list of children, None indicates no check has been made, empty list means no children, otherwise list of names or list of EMDirEntrys
		self.dirregex = dirregex	# only list files using regex
		self.root = str(root)	# Path prefixing name
		self.name = str(name)	# name of this path element (string)
		self.index = str(index)
		self.hidedot = hidedot	# If set children beginning with . will be hidden
		
		if self.root[-1] == "/" or self.root[-1] == "\\" : self.root = self.root[:-1]

		# self.seq = EMDirEntry.classcount
		#EMDirEntry.classcount += 1

		if name[:4].lower() == "bdb:" :
			self.isbdb = True
			self.name = name[4:]
		else :
			self.isbdb = False

		if self.isbdb :
			self.filepath = os.path.join(self.root, "EMAN2DB", self.name+".bdb")
		else : self.filepath = os.path.join(self.root, self.name)

		try : stat = os.stat(self.filepath)
		except : stat = (0, 0, 0, 0, 0, 0, 0, 0, 0)

		if not self.isbdb : self.size = stat[6]		# file size (integer, bytes)
		else : self.size = "-"
		self.date = local_datetime(stat[8])			# modification date (string: yyyy/mm/dd hh:mm:ss)

		# These can be expensive so we only get them on request, or if they are fast

		self.dim = None			# dimensions (string)
		self.filetype = None		# file type (string, "Folder" for directory)
		self.nimg = None			# number of images in file (int)

		# Directories don't really have extra info

		if os.path.isdir(self.filepath) :
			self.filetype = "Folder"
			self.dim = ""
			self.nimg = ""
			self.size = ""
		
#		print "Init DirEntry ", self, self.__dict__

	def __repr__(self) :
		return "<%s %s>"%(self.__class__.__name__ , self.path())

	def path(self) :
		"""The full path of the current item"""

		if self.isbdb : return "bdb:%s#%s"%(self.root, self.name)

		return os.path.join(self.root, self.name).replace("\\", "/")

	def truepath(self) :
		"""The actual path to the file for the current item"""

		if self.isbdb : return "%s/EMAN2DB/%s.bdb"%(self.root, self.name)

		return os.path.join(self.root, self.name).replace("\\", "/")

	def fileTypeClass(self) :
		"""Returns the FileType class corresponding to the named filetype if it exists. None otherwise"""
		try:
			filetype = EMFileType.typesbyft[self.filetype]
			return filetype
		except:
			return None

	def sort(self, column, order) :
		"""Recursive sorting"""

		if self.__children == None or len(self.__children) == 0 or isinstance(self.__children[0], str) : return

		self.__children.sort(key = self.__class__.col[column], reverse = order)

		# This keeps the row numbers consistent

		for i, child in enumerate(self.__children) :
			child.index = i

	def findSelected(self, ret) :
		"""Used to retain selection during sorting. Returns a list of (parent, row) pairs."""

		if self.__children == None or len(self.__children) == 0 or isinstance(self.__children[0], str) : return

#		if len(ret) == 0 : print "findsel"

		for i, child in enumerate(self.__children) :
			try :
				if child.sel :
#					print "sel ->", child.name
					child.sel = False
					ret.append((self, i))
			except :
#				print "not ->", child.name
				pass

			child.findSelected(ret)

	def parent(self) :
		"""Return the parent"""

		return self.__parent

	def child(self, n) :
		"""Returns nth child or None"""

		self.fillChildEntries()

		try : return self.__children[n]
		except :
			print "Request for child %d of children %s (%d)"%(n, self.__children, len(self.__children))
			traceback.print_stack()
			raise Exception

	def nChildren(self) :
		"""Count of children"""

		self.fillChildNames()

#		print "EMDirEntry.nChildren(%s) = %d"%(self.filepath, len(self.__children))

		return len(self.__children)

	def fillChildNames(self) :
		"""Makes sure that __children contains at LEAST a list of names. This function needs to reimplmented to make derived browsers, 
		NOTE!!!! You must have nimg implmented in your reimplmentation (I know this is a bad design....)"""

		if self.__children == None :
			if not os.path.isdir(self.filepath) :
				self.__children = []
				return

			# read the child filenames

			if self.hidedot :
				filelist = [i for i in os.listdir(self.filepath) if i[0] != '.']
			else :
				filelist = os.listdir(self.filepath)

			# Weed out undesirable files

			self.__children = []

			if self.dirregex != None :
				for child in filelist :
#					print child, self.dirregex.search (child)

					ctt=self.filepath + "/" + child
					have_dir = os.path.isdir(ctt)
					if not (have_dir or os.path.isfile(ctt) or os.path.islink(ctt)) : continue

					if isinstance(self.dirregex, str) :
						if have_dir :
							chl = child + ".dir"
						else :
							chl = child
						matching = matches_pats(chl, self.dirregex)
						have_dir = False
					else :
						matching = self.dirregex.match(child) != None

					if have_dir or matching :
						self.__children.append(child)
			# elif self.regex :
				# for child in filelist :
					# if not self.regex.search(child) or child == "EMAN2DB" :
						#s elf.__children.append(child)
			else :
				self.__children = filelist

			if "EMAN2DB" in self.__children :
				self.__children.remove("EMAN2DB")

				if self.dirregex != None :
					if isinstance(self.dirregex, str) :
						t = ["bdb:"+i for i in db_list_dicts("bdb:"+self.filepath) if matches_pats(i, self.dirregex)]
					else :
						t = ["bdb:"+i for i in db_list_dicts("bdb:"+self.filepath) if self.dirregex.match(i) != None]

#					for i in db_list_dicts("bdb:"+self.filepath) : print i, self.dirregex.search (i)
				else :
					t = ["bdb:"+i for i in db_list_dicts("bdb:"+self.filepath)]

				self.__children.extend(t)

			self.__children.sort( )

#			print self.path( ), self.__children

	def fillChildEntries(self) :
		"""Makes sure that __children have been filled with EMDirEntries when appropriate"""

		if self.__children == None : self.fillChildNames()
		if len(self.__children) == 0 : return		# nothing to do. No children
		if not isinstance(self.__children[0], str) : return 		# this implies entries have already been filled

		for i, n in enumerate(self.__children) :
			self.__children[i] = self.__class__(self.filepath, n, i, self)

	def statFile(self, filename) :
		"""Stat either a file or BDB"""

		if filename[:4].lower() == "bdb:" :
			path, dictname, keys = db_parse_path(filename)
			path = path + "/EMAN2DB/" + dictname + ".bdb"

			return os.stat(path)
		else :
			try :
				return os.stat(filename)
			except :
				return os.stat(".")				# this is just a failsafe for things like broken links

	def cache_old(self, check_info = True) :
		"""Makes sure that the cache has been updated since the last file change or info file change.
		Returns 0 if the cache should be up to date.
		Returns 1 if the image file has changed
		Returns 2 if the info file has changed
		Returns 3 if both have changed

		This method is largely called in subclasses."""

		tp = self.truepath()

		if self.updtime >= os.stat(tp).st_mtime : change = 0
		else : change = 1

		if check_info :
			try :
				if self.updtime < os.stat(info_name(tp)).st_mtime : change += 2
			except : pass		# If there is no info file, we don't have to worry about it

		return change

	def fillDetails(self) :
		"""Fills in the expensive metadata about this entry. Returns False if no update was necessary.
		Returns 0 if nothing was done
       		 1 if the metadata was filled in by probing the file
       		 2 if the metadata was filled in from the cache"""

		if self.filetype != None : return 0		# must all ready be filled in

#		print "X %s\t%s\t%s"%(self.root, self.name, self.path())

		cache = None
		cachename = self.name+"!main"

		try :
			cache = js_open_dict(self.root+"/.browsercache.json")
			self.updtime, self.dim, self.filetype, self.nimg, self.size = cache[cachename]		# try to read the cache for the current file

			if self.cache_old(False) == 0 : return 2 		# current cache, no further update necessary
		except :
			pass

		# Check the cache for metadata

		if not (os.path.isfile(self.path()) or os.path.islink(self.path())) : 
			self.filetype="SPECIAL"
			return 0

		name = 'browser'

		# BDB details are already cached and can often be retrieved quickly

		if self.isbdb :
			self.filetype = "BDB"

			try :
				info = db_get_image_info(self.path())
#				print self.path, info
				self.nimg = info[0]

				if self.nimg > 0 :
					if info[1][1] == 1 : self.dim = str(info[1][0])
					elif info[1][2] == 1 : self.dim = "%d x %d"%(info[1][0], info[1][1])
					else : self.dim = "%d x %d x %d"%(info[1][0], info[1][1], info[1][2])

					self.size = info[1][0]*info[1][1]*info[1][2]*4*self.nimg
				else :
					self.dim = "-"
			except :
				traceback.print_exc()
				self.nimg = -1
				self.dim = "-"

			if cache != None : cache.setval(cachename, (time.time(), self.dim, self.filetype, self.nimg, self.size), True)

			return 1

		# we do this this way because there are so many possible image file exensions, and sometimes
		# people use a non-standard one (esp for MRC files)

		try : self.nimg = EMUtil.get_image_count(self.path())
		except : self.nimg = 0

		# we have an image file

		if self.nimg > 0 :
			try : tmp = EMData(self.path(), 0, True)		# try to read an image header for the file
			except :
				for i in xrange(1, 10) :
					try : tmp = EMData(self.path(), i, True)
					except : continue
					break
				if i == 9 :
					print "Error: all of the first 10 images are missing ! : ",self.path()
					return

			if tmp["ny"] == 1 : self.dim = str(tmp["nx"])
			elif tmp["nz"] == 1 : self.dim = "%d x %d"%(tmp["nx"], tmp["ny"])
			else : self.dim = "%d x %d x %d"%(tmp["nx"], tmp["ny"], tmp["nz"])

			if self.nimg == 1 : self.filetype = "Image"
			else : self.filetype = "Image Stack"

		# Ok, we need to try to figure out what kind of file this is

		else :
			try :
				head = file(self.path(), "rb").read(4096)		# Most FileTypes should be able to identify themselves using the first 4K block of a file

				try : guesses = EMFileType.extbyft[os.path.splitext(self.path())[1]]		# This will get us a list of possible FileTypes for this extension
				except : guesses = EMFileType.alltocheck

	#			print "-------\n", guesses

				for guess in guesses :
					try : size, n, dim = guess.isValid(self.path(), head)		# This will raise an exception if isValid returns False
					except : continue

					# If we got here, we found a match
					self.filetype = guess.name()
					self.dim = dim
					self.nimg = n
					self.size = size

					break
				else :		# this only happens if no match was found
					self.filetype = "-"
					self.dim = "-"
					self.nimg = "-"
			except :
				self.filetype = "-"
				self.dim = "-"
				self.nimg = "-"

		if cache != None :
			try : cache.setval(cachename, (time.time(), self.dim, self.filetype, self.nimg, self.size), True)
			except : pass
		return 1

#---------------------------------------------------------------------------

def nonone(val) :
	"""Returns '-' for None, otherwise the string representation of the passed value"""

	try :
		if val != None : return str(val)
		return "-"
	except : return "X"

def humansize(val) :
	"""Representation of an integer in readable form"""

	try : val = int(val)
	except : return val

	if val > 1000000000 : return "%d g"%(val/1000000000)
	elif val > 1000000 : return "%d m"%(val/1000000)
	elif val > 1000 : return "%d k"%(val/1000)
	return str(val)

#---------------------------------------------------------------------------

class EMFileItemModel(QtCore.QAbstractItemModel) :
	"""This ItemModel represents the local filesystem. We don't use the normal filesystem item model because we want
	to provide more info on images, and we need to merge BDB: files into the file view."""

	headers = ("Row", "Name", "Type", "Size", "Dim", "N", "Date")

	def __init__(self, startpath = None, direntryclass = EMDirEntry, dirregex = None) :
		QtCore.QAbstractItemModel.__init__(self)
		if startpath[:2] == "./" : startpath = startpath[2:]
		self.root = direntryclass(startpath, "", 0, dirregex = dirregex)				# EMDirEntry as a parent for the root path
		self.rootpath = startpath							# root path for current browser
		self.last = (0, 0)
		self.db = None
#		print "Init FileItemModel ", self, self.__dict__

	def canFetchMore(self, idx) :
		"""The data is loaded lazily for children already. We don't generally
		need to worry about there being SO many file that this is necessary, so it always 
		returns False."""

		return False

	def columnCount(self, parent) :
		"""Always 7 columns"""

		#print "EMFileItemModel.columnCount() = 6"

		return 7

	def rowCount(self, parent) :
		"""Returns the number of children for a given parent"""

#		if parent.column() != 0 : return 0

		if parent != None and parent.isValid() :
#			print "rowCount(%s) = %d"%(str(parent), parent.internalPointer().nChildren())

			return parent.internalPointer().nChildren()

#		print "rowCount(root) = %d"%self.root.nChildren()

		return self.root.nChildren()

	def data(self, index, role) :
		"""Returns the data for a specific location as a string"""

		if not index.isValid() : return None
		if role != Qt.DisplayRole : return None

		data = index.internalPointer()

		if data == None :
			print "Error with index ", index.row(), index.column()
			return "XXX"

		#if index.column() == 0 : print "EMFileItemModel.data(%d %d %s) = %s"%(index.row(), index.column(), index.parent(), str(data.__dict__))

		col = index.column()

		if col == 0 :
			return nonone(data.index)
		elif col == 1 :
			if data.isbdb : return "bdb:"+data.name
			return nonone(data.name)
		elif col == 2 : return nonone(data.filetype)
		elif col == 3 : return humansize(data.size)
		elif col == 4 :
			if data.dim == 0 : return "-"
			return nonone(data.dim)
		elif col == 5 : return nonone(data.nimg)
		elif col == 6 : return nonone(data.date)

	def headerData(self, sec, orient, role) :
		# In case we use the QTableViews

		if orient == Qt.Vertical :
			if role == Qt.DisplayRole :
				return str(sec)
			elif role == Qt.ToolTipRole :
				return None

		# This works for all Q*views

		if orient == Qt.Horizontal :
			if role == Qt.DisplayRole :
				return self.__class__.headers[sec]
			elif role == Qt.ToolTipRole :
				return None								# May fill this in later

		return None

	def hasChildren(self, parent) :
		"""Returns whether the index 'parent' has any child items"""

		#print "EMFileItemModel.hasChildren(%d, %d, %s)"%(parent.row(), parent.column(), str(parent.internalPointer()))
#		if parent.column() != 0 : return False

		try :
			if parent != None and parent.isValid() :
				if parent.internalPointer().nChildren() > 0 : return True
				else : return False
			return True
		except : return False

	def hasIndex(self, row, col, parent) :
		"""Test if the specified index would exist"""

#		print "EMFileItemModel.hasIndex(%d, %d, %s)"%(row, column, parent.internalPointer())

		try :
			if parent != None and parent.isValid() :
				data = parent.internalPointer().child(row)
			else : data = self.root.child(row)
		except :
			traceback.print_exc()
			return False
		return True

	def index(self, row, column, parent) :
		"""produces a new QModelIndex for the specified item"""

#		if column == 0 : print "Index :", row, column, parent.internalPointer(),

		try :
			if parent != None and parent.isValid() :
				data = parent.internalPointer().child(row)
			else :
				data = self.root.child(row)
		except :
			traceback.print_exc()
#			print "None"
			return QtCore.QModelIndex()			# No data, return invalid

#		if column == 0 : print data

		return self.createIndex(row, column, data)

	def parent(self, index) :
		"""Returns the parent of the specified index"""

#		print "parent ", index.row()

		if index.isValid() :
			try : data = index.internalPointer().parent()
			except :
				print "Parent index error:", str(index.__dict__)
		else : return QtCore.QModelIndex()

		if data == None : return QtCore.QModelIndex()			# No data, return invalid

		return self.createIndex(index.row(), 0, data)		# parent is always column 0

	def sort(self, column, order) :
		"""Trigger recursive sorting"""

		if column < 0 : return

		self.root.sort(column, order)
#		self.emit(QtCore.SIGNAL("layoutChanged()"))
		self.layoutChanged.emit()

	def findSelected(self, toplevel = True) :
		"""Makes a list of QModelIndex items for all items marked as selected"""

#		print "findsel"

		sel = []
		self.root.findSelected(sel)

		if toplevel : return [self.createIndex(i[1], 0, i[0]) for i in sel if i[0] == self.root]

		return [self.createIndex(i[1], 0, i[0]) for i in sel]

	def details(self, index) :
		"""This will trigger loading the (expensive) details about the specified index, and update the display"""

		if not index.isValid() : return

		if index.internalPointer().fillDetails() :
			self.dataChanged.emit(index, self.createIndex(index.row(), 5, index.internalPointer()))

#---------------------------------------------------------------------------

class myQItemSelection(QtGui.QItemSelectionModel) :
	"""For debugging"""

	def select(self, tl, br) :
		print tl.indexes()[0].row(), tl.indexes()[0].column(), int(br)
		QtGui.QItemSelectionModel.select(self, tl, QtGui.QItemSelectionModel.SelectionFlags(QtGui.QItemSelectionModel.ClearAndSelect+QtGui.QItemSelectionModel.Rows))

#---------------------------------------------------------------------------

class EMInfoPane(QtGui.QWidget) :
	"""Subclasses of this class will be used to display information about specific files. Each EMFileType class will return the
	pointer to the appropriate infoPane subclass for displaying information about the file it represents. The subclass instances
	are allocated by the infoWin class"""

	def __init__(self, parent = None) :
		"""Set our GUI up"""

		QtGui.QWidget.__init__(self, parent)

		# self.setTitle("e2dispaly.py Information Pane")

		self.setWindowTitle("e2display.py Information Pane") # Jesus

		# Root class represents no target

		self.hbl = QtGui.QHBoxLayout(self)
		self.lbl = QtGui.QLabel("No Information Available")
		self.hbl.addWidget(self.lbl)

	def display(self, target) :
		"""display information for the target EMDirEntry with EMFileType ftype"""

		# self.setTitle("e2dispaly.py Information Pane")

		self.target = target
		self.setWindowTitle("e2display.py Information Pane") # Jesus

		return

	def busy(self) :
		pass

	def notbusy(self) :
		pass

class EMTextInfoPane(EMInfoPane) :
	def __init__(self, parent = None) :
		QtGui.QWidget.__init__(self, parent)

		self.vbl = QtGui.QVBoxLayout(self)

		# text editing widget

		self.text = QtGui.QTextEdit()
		self.text.setAcceptRichText(False)
		self.text.setReadOnly(True)
		self.vbl.addWidget(self.text)

		# Find box

		self.wfind = StringBox(label = "Find:")
		self.vbl.addWidget(self.wfind)

		# Buttons

		self.hbl = QtGui.QHBoxLayout()

		self.wbutedit = QtGui.QPushButton("Edit")
		self.hbl.addWidget(self.wbutedit)

		self.wbutcancel = QtGui.QPushButton("Revert")
		self.wbutcancel.setEnabled(False)
		self.hbl.addWidget(self.wbutcancel)

		self.wbutok = QtGui.QPushButton("Save")
		self.wbutok.setEnabled(False)
		self.hbl.addWidget(self.wbutok)

		self.vbl.addLayout(self.hbl)

		QtCore.QObject.connect(self.wfind, QtCore.SIGNAL("valueChanged"), self.find)
		QtCore.QObject.connect(self.wbutedit, QtCore.SIGNAL('clicked(bool)'), self.buttonEdit)
		QtCore.QObject.connect(self.wbutcancel, QtCore.SIGNAL('clicked(bool)'), self.buttonCancel)
		QtCore.QObject.connect(self.wbutok, QtCore.SIGNAL('clicked(bool)'), self.buttonOk)

	def display(self, data) :
		"""display information for the target EMDirEntry"""

		self.target = data

		self.text.setPlainText(file(self.target.path(), "r").read())

		self.text.setReadOnly(True)
		self.wbutedit.setEnabled(True)
		self.wbutcancel.setEnabled(False)
		self.wbutok.setEnabled(False)

	def find(self, value) :
		"""Find a string"""

		if not self.text.find(value) :
			# this implements wrapping
			self.text.moveCursor(1, 0)
			self.text.find(value)

	def buttonEdit(self, tog) :
		self.text.setReadOnly(False)
		self.wbutedit.setEnabled(False)
		self.wbutcancel.setEnabled(True)
		self.wbutok.setEnabled(True)

	def buttonCancel(self, tog) :
		self.display(self.target)

	def buttonOk(self, tog) :
		try : file(self.target.path(), "w").write(str(self.text.toPlainText()))
		except : QtGui.QMessageBox.warning(self, "Error !", "File write failed")

#---------------------------------------------------------------------------

class EMHTMLInfoPane(EMInfoPane) :
	def __init__(self, parent = None) :
		QtGui.QWidget.__init__(self, parent)

		self.vbl = QtGui.QVBoxLayout(self)

		# text editing widget

		self.text = QtGui.QTextEdit()
		self.text.setAcceptRichText(True)
		self.text.setReadOnly(True)
		self.vbl.addWidget(self.text)

		# Find box

		self.wfind = StringBox(label = "Find:")
		self.vbl.addWidget(self.wfind)

		# Buttons

		self.hbl = QtGui.QHBoxLayout()

		self.wbutedit = QtGui.QPushButton("Edit")
		self.hbl.addWidget(self.wbutedit)

		self.wbutcancel = QtGui.QPushButton("Revert")
		self.wbutcancel.setEnabled(False)
		self.hbl.addWidget(self.wbutcancel)

		self.wbutok = QtGui.QPushButton("Save")
		self.wbutok.setEnabled(False)
		self.hbl.addWidget(self.wbutok)

		self.vbl.addLayout(self.hbl)

		QtCore.QObject.connect(self.wfind, QtCore.SIGNAL("valueChanged"), self.find)
		QtCore.QObject.connect(self.wbutedit, QtCore.SIGNAL('clicked(bool)'), self.buttonEdit)
		QtCore.QObject.connect(self.wbutcancel, QtCore.SIGNAL('clicked(bool)'), self.buttonCancel)
		QtCore.QObject.connect(self.wbutok, QtCore.SIGNAL('clicked(bool)'), self.buttonOk)

	def display(self, data) :
		"""display information for the target EMDirEntry"""

		self.target = data
		self.text.setHtml(file(self.target.path(), "r").read())
		self.text.setReadOnly(True)
		self.wbutedit.setEnabled(True)
		self.wbutcancel.setEnabled(False)
		self.wbutok.setEnabled(False)

	def find(self, value) :
		"""Find a string"""

		if not self.text.find(value) :
			# this implements wrapping
			self.text.moveCursor(1, 0)
			self.text.find(value)

	def buttonEdit(self, tog) :
		self.text.setReadOnly(False)
		self.wbutedit.setEnabled(False)
		self.wbutcancel.setEnabled(True)
		self.wbutok.setEnabled(True)

	def buttonCancel(self, tog) :
		self.display(self.target)

	def buttonOk(self, tog) :
		try : file(self.target.path(), "w").write(str(self.text.toHtml()))
		except : QtGui.QMessageBox.warning(self, "Error !", "File write failed")

#---------------------------------------------------------------------------

class EMPDBInfoPane(EMInfoPane) :
	def __init__(self, parent = None) :
		QtGui.QWidget.__init__(self, parent)
		self.vbl = QtGui.QVBoxLayout(self)
		# text editing widget
		self.text = QtGui.QTextEdit()
		self.text.setAcceptRichText(False)
		self.text.setReadOnly(True)
		self.vbl.addWidget(self.text)
		# Find box
		self.wfind = StringBox(label = "Find:")
		self.vbl.addWidget(self.wfind)
		# Buttons
		self.hbl = QtGui.QHBoxLayout()
		self.wbutedit = QtGui.QPushButton("Edit")
		self.hbl.addWidget(self.wbutedit)
		self.wbutcancel = QtGui.QPushButton("Revert")
		self.wbutcancel.setEnabled(False)
		self.hbl.addWidget(self.wbutcancel)
		self.wbutok = QtGui.QPushButton("Save")
		self.wbutok.setEnabled(False)
		self.hbl.addWidget(self.wbutok)
		self.vbl.addLayout(self.hbl)
		QtCore.QObject.connect(self.wfind, QtCore.SIGNAL("valueChanged"), self.find)
		QtCore.QObject.connect(self.wbutedit, QtCore.SIGNAL('clicked(bool)'), self.buttonEdit)
		QtCore.QObject.connect(self.wbutcancel, QtCore.SIGNAL('clicked(bool)'), self.buttonCancel)
		QtCore.QObject.connect(self.wbutok, QtCore.SIGNAL('clicked(bool)'), self.buttonOk)

	def display(self, data) :
		"""display information for the target EMDirEntry"""
		self.target = data
		self.text.setPlainText(file(self.target.path(), "r").read())
		self.text.setReadOnly(True)
		self.wbutedit.setEnabled(True)
		self.wbutcancel.setEnabled(False)
		self.wbutok.setEnabled(False)

	def find(self, value) :
		"""Find a string"""
		if not self.text.find(value) :
			# this implements wrapping
			self.text.moveCursor(1, 0)
			self.text.find(value)

	def buttonEdit(self, tog) :
		self.text.setReadOnly(False)
		self.wbutedit.setEnabled(False)
		self.wbutcancel.setEnabled(True)
		self.wbutok.setEnabled(True)

	def buttonCancel(self, tog) :
		self.display(self.target)

	def buttonOk(self, tog) :
		try : file(self.target.path(), "w").write(str(self.text.toPlainText()))
		except : QtGui.QMessageBox.warning(self, "Error !", "File write failed")

#---------------------------------------------------------------------------

class EMPlotInfoPane(EMInfoPane) :
	def __init__(self, parent = None) :
		QtGui.QWidget.__init__(self, parent)

		self.gbl = QtGui.QGridLayout(self)

		# List as alternate mechanism for selecting image number(s)

		self.plotdata = QtGui.QTableWidget()
		self.gbl.addWidget(self.plotdata, 0, 0)

	def display(self, target) :
		"""display information for the target EMDirEntry"""

		self.target = target
		self.plotdata.clear()

		# read the data into a list of lists

		numc = 0
		data = []

		for l in file(target.path(), "r") :
			if l[0] == "#" : continue

			vals = [float(i) for i in renumfind.findall(l)]

			if len(vals) == 0 : continue

			if numc == 0 : numc = len(vals)
			elif numc != len(vals) : break
			data.append(vals)

			if len(data) == 2500 : break			# if the table is too big, we just do a ...

		if len(data) == 2500 : self.plotdata.setRowCount(2501)
		else : self.plotdata.setRowCount(len(data))

		self.plotdata.setColumnCount(numc)
		self.plotdata.setVerticalHeaderLabels([str(i) for i in xrange(len(data))])
		self.plotdata.setHorizontalHeaderLabels([str(i) for i in xrange(numc)])

		for r in xrange(len(data)) :
			for c in xrange(numc) :
				self.plotdata.setItem(r, c, QtGui.QTableWidgetItem("%1.4g"%data[r][c]))

		if len(data) == 2500 :
			self.plotdata.setVerticalHeaderItem(2500, QtGui.QTableWidgetItem("..."))

#---------------------------------------------------------------------------

class EMFolderInfoPane(EMInfoPane) :
	def __init__(self, parent = None) :
		QtGui.QWidget.__init__(self, parent)

		self.vbl = QtGui.QVBoxLayout(self)

	def display(self, target) :
		"""display information for the target EMDirEntry"""

		self.target = target

#---------------------------------------------------------------------------

class EMBDBInfoPane(EMInfoPane) :
	maxim = 500

	def __init__(self, parent = None) :
		QtGui.QWidget.__init__(self, parent)

		self.gbl = QtGui.QGridLayout(self)

		# Spinbox for selecting image number

		self.wimnum = QtGui.QSpinBox()
		self.wimnum.setRange(0, 0)
		self.gbl.addWidget(self.wimnum, 0, 0)

		# List as alternate mechanism for selecting image number(s)

		self.wimlist = QtGui.QListWidget()
		self.gbl.addWidget(self.wimlist, 1, 0)

		# Actual header contents

		self.wheadtree = QtGui.QTreeWidget()
		self.wheadtree.setColumnCount(2)
		self.wheadtree.setHeaderLabels(["Item", "Value"])
		self.gbl.addWidget(self.wheadtree, 0, 1, 2, 1)

		self.gbl.setColumnStretch(0, 1)
		self.gbl.setColumnStretch(1, 4)

		# Lower region has buttons for actions

		self.hbl2 = QtGui.QGridLayout()

		self.wbutmisc = []

		# 10 buttons for context-dependent actions

		self.hbl2.setRowStretch(0, 1)
		self.hbl2.setRowStretch(1, 1)

		for i in range(5) :
			self.hbl2.setColumnStretch(i, 2)
	
			for j in range(2) :
				self.wbutmisc.append(QtGui.QPushButton(""))
				self.hbl2.addWidget(self.wbutmisc[-1], j, i)
				self.wbutmisc[-1].hide()
				QtCore.QObject.connect(self.wbutmisc[-1], QtCore.SIGNAL('clicked(bool)'), lambda x, v = i*2+j :self.buttonMisc(v))

		# These just clean up the layout a bit

		self.wbutxx = QtGui.QLabel("")
		self.wbutxx.setMaximumHeight(12)
		self.hbl2.addWidget(self.wbutxx, 0, 6)
		self.wbutyy = QtGui.QLabel("")
		self.wbutyy.setMaximumHeight(12)
		self.hbl2.addWidget(self.wbutyy, 1, 6)

		self.gbl.addLayout(self.hbl2, 2, 0, 1, 2)

		QtCore.QObject.connect(self.wimnum, QtCore.SIGNAL("valueChanged(int)"), self.imNumChange)
		QtCore.QObject.connect(self.wimlist, QtCore.SIGNAL("itemSelectionChanged()"), self.imSelChange)
##		QtCore.QObject.connect(self.wbutedit, QtCore.SIGNAL('clicked(bool)'), self.buttonEdit)

		self.view2d = []
		self.view3d = []
		self.view2ds = []
		self.viewplot2d = []

	def hideEvent(self, event) :
		"""If this pane is no longer visible close any child views"""

		for v in self.view2d : v.close()
		for v in self.view3d : v.close()
		for v in self.view2ds : v.close()
		event.accept()

	def buttonMisc(self, but) :
		"""Context sensitive button press"""

		val = self.wimlist.currentItem().text()
		try :
			val = int(val)
		except :
			QtGui.QMessageBox.warning(self, "Error", "Sorry, cannot display string-keyed images")
			return

		self.curft.setN(val)
		self.curactions[but][2](self)				# This calls the action method

	def display(self, target) :
		"""display information for the target EMDirEntry"""

		self.target = target
		self.bdb = db_open_dict(self.target.path())

		# Set up image selectors for stacks

		if target.nimg == 0 :
			self.wimnum.hide()
			k = self.bdb.keys()
			k.sort()
			self.wimlist.addItems(k)
			self.wimlist.show()
			self.curim = 0
		else :
			self.wimnum.setRange(0, target.nimg)
			self.wimlist.clear()
			self.wimlist.addItems([str(i) for i in range(0, min(target.nimg, self.maxim))])

			if target.nimg > self.maxim : self.wimlist.addItem("...")

			self.wimnum.show()
			self.wimlist.show()

		self.wheadtree.clear()

		# This sets up action buttons which can be used on individual images in a stack

		try :
			self.curft = EMImageFileType(target.path())
			self.curactions = self.curft.actions()

			for i, b in enumerate(self.wbutmisc) :
				try :
					b.setText(self.curactions[i][0])
					b.setToolTip(self.curactions[i][1])
					b.show()
				except :
					b.hide()
		except :
			# Not a readable image or volume

			pass

	def imNumChange(self, num) :
		"""New image number"""

		if num < 500 : self.wimlist.setCurrentRow(num)
		else : self.showItem(num)

	def imSelChange(self) :
		"""New image selection"""

		val = self.wimlist.currentItem().text()

		try :
			val = int(val)
			self.wimnum.setValue(val)
		except :
			val = str(val)

		self.showItem(val)

	def showItem(self, key) :
		"""Shows header information for the selected item"""

		self.wheadtree.clear()
		trg = self.bdb.get_header(key)

		# print "Warning: tried to read unavailable key: %s"%key

		if trg == None and key != "..." :
			#print self.bdb.keys()
			self.addTreeItem("*None*")
			return

		self.addTreeItem(trg)

	def addTreeItem(self, trg, parent = None) :
		"""(recursively) add an item to the tree"""

		itms = []

		# Dictionaries may require recursion

		if isinstance(trg, dict) :
			for k in sorted(trg.keys()) :
				itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((str(k), str(trg[k])))))
				if isinstance(trg[k], list) or isinstance(trg[k], tuple) or isinstance(trg[k], set) or isinstance(trg[k], dict) :
					self.addTreeItem(trg[k], itms[-1])
		elif isinstance(trg, list) or isinstance(trg, tuple) or isinstance(trg, set) :
			for k in trg :
				if isinstance(k, list) or isinstance(k, tuple) or isinstance(k, set) or isinstance(k, dict) :
					try : itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((k.__class__.__name__, ""))))
					except : itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList(("??", ""))))
					self.addTreeItem(k, itms[-1])
				else :
					itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((str(k), ""))))
		else :
			itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((str(trg), ""))))

		if parent == None :
			self.wheadtree.addTopLevelItems(itms)
			self.wheadtree.resizeColumnToContents(0)
		else : parent.addChildren(itms)

#---------------------------------------------------------------------------

class EMJSONInfoPane(EMInfoPane) :
	def __init__(self, parent = None) :
		QtGui.QWidget.__init__(self, parent)

		self.gbl = QtGui.QGridLayout(self)

		# List of keys

		self.wkeylist = QtGui.QListWidget()
		self.gbl.addWidget(self.wkeylist, 1, 0)

		# contents of a single key

		self.wheadtree = QtGui.QTreeWidget()
		self.wheadtree.setColumnCount(2)
		self.wheadtree.setHeaderLabels(["Key/#", "Value"])
		self.gbl.addWidget(self.wheadtree, 0, 1, 2, 1)

		self.gbl.setColumnStretch(0, 1)
		self.gbl.setColumnStretch(1, 4)

		# Lower region has buttons for actions

		self.hbl2 = QtGui.QGridLayout()

		self.wbutmisc = []

		# 10 buttons for context-dependent actions

		self.hbl2.setRowStretch(0, 1)
		self.hbl2.setRowStretch(1, 1)

		for i in range(5) :
			self.hbl2.setColumnStretch(i, 2)

			for j in range(2) :
				self.wbutmisc.append(QtGui.QPushButton(""))
				self.hbl2.addWidget(self.wbutmisc[-1], j, i)
				self.wbutmisc[-1].hide()
				QtCore.QObject.connect(self.wbutmisc[-1], QtCore.SIGNAL('clicked(bool)'), lambda x, v = i*2+j :self.buttonMisc(v))

		# These just clean up the layout a bit

		self.wbutxx = QtGui.QLabel("")
		self.wbutxx.setMaximumHeight(12)
		self.hbl2.addWidget(self.wbutxx, 0, 6)
		self.wbutyy = QtGui.QLabel("")
		self.wbutyy.setMaximumHeight(12)
		self.hbl2.addWidget(self.wbutyy, 1, 6)

		self.gbl.addLayout(self.hbl2, 2, 0, 1, 2)

		QtCore.QObject.connect(self.wkeylist, QtCore.SIGNAL("itemSelectionChanged()"), self.imSelChange)
		QtCore.QObject.connect(self.wheadtree, QtCore.SIGNAL("itemExpanded(QTreeWidgetItem*)"), self.treeExp)
		QtCore.QObject.connect(self.wheadtree, QtCore.SIGNAL("itemCollapsed(QTreeWidgetItem*)"), self.treeExp)
		QtCore.QObject.connect(self.wheadtree, QtCore.SIGNAL("itemSelectionChanged()"), self.treeSel)
		QtCore.QObject.connect(self.wheadtree, QtCore.SIGNAL("itemActivated(QTreeWidgetItem*, int)"), self.treeAct)
##		QtCore.QObject.connect(self.wbutedit, QtCore.SIGNAL('clicked(bool)'), self.buttonEdit)
		self.view2d = []
		self.view3d = []
		self.view2ds = []

	def hideEvent(self, event) :
		"""If this pane is no longer visible close any child views"""

		for v in self.view2d : v.close()
		for v in self.view3d : v.close()
		for v in self.view2ds : v.close()
		event.accept()

	def buttonMisc(self, but) :
		"""Context sensitive button press"""

		val = self.wkeylist.currentItem().text()

		# self.curft.setN(val)
		# self.curactions[but][2](self)				# This calls the action method

	def display(self, target) :
		"""display information for the target EMDirEntry"""

		self.target = target
		self.js = js_open_dict(self.target.path())

		# Set up image selectors for stacks

		k = self.js.keys()
		k.sort()
		self.wkeylist.addItems(k)
		self.wkeylist.show()
		self.curim = 0

		self.wheadtree.clear()

		# This sets up action buttons which can be used on individual images in a stack

	def imSelChange(self) :
		"""New image selection"""

		val = str(self.wkeylist.currentItem().text())

		self.showItem(val)

	def showItem(self, key) :
		"""Shows header information for the selected item"""

		self.wheadtree.clear()
		trg = self.js[key]

		if trg == None :
#			print "Warning: tried to read unavailable key: %s"%key

			self.addTreeItem("*None*")
			return

		self.addTreeItem(trg)
		self.wheadtree.resizeColumnToContents(0)
		self.wheadtree.resizeColumnToContents(1)

	def treeExp(self, item) :
		"""Make sure the tree columns get resized when the user expands/contracts the content"""

		self.wheadtree.resizeColumnToContents(1)
		self.wheadtree.resizeColumnToContents(0)

	def treeSel(self) :
		"""When the selection changes in the tree"""

		#FIXME - Not implemented yet

	def treeAct(self, item, col) :
		"""When a tree item is 'activated' """

		#FIXME - Not implemented yet

	def addTreeItem(self, trg, parent = None) :
		"""(recursively) add an item to the tree"""

		itms = []

		# Dictionaries may require recursion

		if isinstance(trg, dict) :
			for k in sorted(trg.keys()) :
				if isinstance(trg[k], (list, tuple, set, dict, EMAN2Ctf)) :
					itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((str(k), ""))))
					self.addTreeItem(trg[k], itms[-1])
				else : itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((str(k), str(trg[k])))))
		elif isinstance(trg, (list, tuple, set)) :
			if isinstance(trg, set) : trg = sorted(trg)		# make a list temporarily
			if len(trg) > 120 : vals = range(0, 50)+[-1]+range(len(trg)-50, len(trg))
			else : vals = xrange(len(trg))
			for k in vals :
				if k == -1 : itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList(("...", "..."))))
				else :
					v = trg[k]
					if isinstance(v, (list, tuple, set, dict, EMAN2Ctf)) :
						itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((str(k), ""))))
						self.addTreeItem(v, itms[-1])
					else : itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((str(k), str(v)))))
		elif isinstance(trg, EMAN2Ctf) :
			itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList(("EMAN2Ctf", ""))))
			subitms = []
			for k, v in trg.to_dict().items() :
				if isinstance(v, (list, tuple)) :
					v = ["%1.3g"%i for i in v]
					subitms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((str(k), ", ".join(v)))))
				else : subitms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((str(k), str(v)))))
			itms[-1].addChildren(subitms)
		elif isinstance(trg, EMData) :
			itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList(("EMData", ""))))
		else :
			itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((str(trg), ""))))

		if parent == None :
			self.wheadtree.addTopLevelItems(itms)
			self.wheadtree.resizeColumnToContents(0)
		else : parent.addChildren(itms)

#---------------------------------------------------------------------------

class EMImageInfoPane(EMInfoPane) :
	maxim = 500

	def __init__(self, parent = None) :
		QtGui.QWidget.__init__(self, parent)

		self.gbl = QtGui.QGridLayout(self)

		# Actual header contents

		self.wheadtree = QtGui.QTreeWidget()
		self.wheadtree.setColumnCount(2)
		self.wheadtree.setHeaderLabels(["Item", "Value"])
		self.gbl.addWidget(self.wheadtree, 0, 0)

	def display(self, target) :
		"""display information for the target EMDirEntry"""

		self.target = target

		self.wheadtree.clear()
		try : trg = EMData(self.target.path(), 0, True).get_attr_dict()		# read the header only, discard the emdata object
		except :
			print "Error reading:", self.target.path(), key

		if trg == None :
			#print "Warning: tried to read unavailable key: %s"%key
	
			self.addTreeItem("*None*")
	
			#print self.bdb.keys()

			return

		self.addTreeItem(trg)

	def addTreeItem(self, trg, parent = None) :
		"""(recursively) add an item to the tree"""

		itms = []

		# Dictionaries may require recursion

		if isinstance(trg, dict) :
			for k in sorted(trg.keys()) :
				itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((str(k), str(trg[k])))))
				if isinstance(trg[k], list) or isinstance(trg[k], tuple) or isinstance(trg[k], set) or isinstance(trg[k], dict) :
					self.addTreeItem(trg[k], itms[-1])
		elif isinstance(trg, list) or isinstance(trg, tuple) or isinstance(trg, set) :
			for k in trg :
				if isinstance(k, list) or isinstance(k, tuple) or isinstance(k, set) or isinstance(k, dict) :
					try : itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((k.__class__.__name__, ""))))
					except : itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList(("??", ""))))
					self.addTreeItem(k, itms[-1])
				else :
					itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((str(k), ""))))
		else :
			itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((str(trg), ""))))

		if parent == None :
			self.wheadtree.addTopLevelItems(itms)
			self.wheadtree.resizeColumnToContents(0)
		else : parent.addChildren(itms)

#---------------------------------------------------------------------------

class EMStackInfoPane(EMInfoPane) :
	maxim = 500

	def __init__(self, parent = None) :
		QtGui.QWidget.__init__(self, parent)

		# self.setWindowTitle("e2display.py Information Pane") # Jesus
		# self.setTitle("e2dispaly.py Information Pane")

		self.gbl = QtGui.QGridLayout(self)

		self.gbl.label1 = QtGui.QLabel("Images") # Jesus
		self.gbl.addWidget(self.gbl.label1, 0, 0) # Jesus

		self.gbl.label2 = QtGui.QLabel("Header Info") # Jesus
		self.gbl.addWidget(self.gbl.label2, 0, 1) # Jesus

		'''Spinbox for selecting image number'''

		self.wimnum = QtGui.QSpinBox()

		# self.wimnum.setRange(0, 0) # JOHN
		# self.gbl.addWidget(self.wimnum, 0, 0) # JOHN

		self.wimnum.setRange(1, 0) # Jesus
		self.gbl.addWidget(self.wimnum, 1, 0) # Jesus

		'''List as alternate mechanism for selecting image number(s)'''

		self.wimlist = QtGui.QListWidget()

		# self.gbl.addWidget(self.wimlist, 1, 0) # JOHN

		self.gbl.addWidget(self.wimlist, 2, 0) # Jesus

		'''Actual header contents'''
		self.wheadtree = QtGui.QTreeWidget()

		# self.wheadtree.setColumnCount(2) #
		self.wheadtree.setHeaderLabels(["Item", "Value"])
		# self.gbl.addWidget(self.wheadtree, 0, 1, 2, 1) # JOHN
		self.gbl.addWidget(self.wheadtree, 1, 1, 2, 1) # Jesus

		self.gbl.setColumnStretch(0, 1)
		self.gbl.setColumnStretch(1, 4)

		'''Lower region has buttons for actions'''

		self.hbl2 = QtGui.QGridLayout()

		self.wbutmisc = []

		'''10 buttons for context-dependent actions'''

		# self.hbl2.setRowStretch(0, 1) # JOHN
		# self.hbl2.setRowStretch(1, 1) # JOHN

		self.hbl2.setRowStretch(1, 1) # Jesus
		self.hbl2.setRowStretch(2, 1) # Jesus

		for i in xrange(5) :
			self.hbl2.setColumnStretch(i, 2)
	
			for j in range(2) :
				self.wbutmisc.append(QtGui.QPushButton(""))
				self.hbl2.addWidget(self.wbutmisc[-1], j, i)
				self.wbutmisc[-1].hide()
				QtCore.QObject.connect(self.wbutmisc[-1], QtCore.SIGNAL('clicked(bool)'), lambda x, v = i*2+j :self.buttonMisc(v))

		# These just clean up the layout a bit

		self.wbutxx = QtGui.QLabel("")
		self.wbutxx.setMaximumHeight(12)
		# self.hbl2.addWidget(self.wbutxx, 0, 6) # JOHN
		self.hbl2.addWidget(self.wbutxx, 1, 6) # Jesus

		self.wbutyy = QtGui.QLabel("")
		self.wbutyy.setMaximumHeight(12)

		# self.hbl2.addWidget(self.wbutyy, 1, 6) # JOHN
		self.hbl2.addWidget(self.wbutyy, 2, 6) # Jesus

		# self.gbl.addLayout(self.hbl2, 2, 0, 1, 2) # JOHN
		self.gbl.addLayout(self.hbl2, 3, 0, 1, 2) # Jesus

		QtCore.QObject.connect(self.wimnum, QtCore.SIGNAL("valueChanged(int)"), self.imNumChange)
		QtCore.QObject.connect(self.wimlist, QtCore.SIGNAL("itemSelectionChanged()"), self.imSelChange)
#		QtCore.QObject.connect(self.wbutedit, QtCore.SIGNAL('clicked(bool)'), self.buttonEdit)
		self.view2d = []
		self.view3d = []
		self.view2ds = []
		self.viewplot2d = []

	def closeEvent(self, event) :
#		E2saveappwin("e2display", "main", self)

		self.updthreadexit = True

		for w in self.view2d+self.view2ds+self.view3d+self.viewplot2d :
			w.close()

		if self.infowin != None :
			self.infowin.close()

		event.accept()

		# self.app().close_specific(self)
		self.emit(QtCore.SIGNAL("module_closed"))

	def hideEvent(self, event) :
		"""If this pane is no longer visible close any child views"""

		for v in self.view2d : v.close()
		event.accept()

	def buttonMisc(self, but) :
		"""Context sensitive button press"""

		self.curft.setN(int(self.wimnum.value()))
		self.curactions[but][2](self)				# This calls the action method

	def display(self, target) :
		"""display information for the target EMDirEntry"""

		self.target = target

		# Set up image selectors for stacks

		self.wimnum.setRange(0, target.nimg-1)
		self.wimlist.clear()
		self.wimlist.addItems([str(i) for i in range(0, min(target.nimg, self.maxim))])

		if target.nimg > self.maxim : self.wimlist.addItem("...")

		self.wimnum.show()
		self.wimlist.show()

		self.wheadtree.clear()

		# This sets up action buttons which can be used on individual images in a stack

		self.curft = EMImageFileType(target.path())
		self.curactions = self.curft.actions()

		for i, b in enumerate(self.wbutmisc) :
			try :
				b.setText(self.curactions[i][0])
				b.setToolTip(self.curactions[i][1])
				b.show()
			except :
				b.hide()

	def imNumChange(self, num) :
		"""New image number"""

		if num < 500 : self.wimlist.setCurrentRow(num)
		else : self.showItem(num)

	def imSelChange(self) :
		"""New image selection"""

		try :
			val = self.wimlist.currentItem().text()
			val = int(val)
			self.wimnum.setValue(val)
		except :
			try :
				val = int(self.wimnum.value())
			except :
				print "Error with key:", val
				return

		self.showItem(val)

	def showItem(self, key) :
		"""Shows header information for the selected item"""

		self.wheadtree.clear()

		try : trg = EMData(self.target.path(), key, True).get_attr_dict()		# read the header only, discard the emdata object
		except :
			print "Error reading:", self.target.path(), key

		if trg == None :
#			print "Warning: tried to read unavailable key: %s"%key

			self.addTreeItem("*None*")

			# print self.bdb.keys()

			return

		self.addTreeItem(trg)

	def addTreeItem(self, trg, parent = None) :
		"""(recursively) add an item to the tree"""

		itms = []

		# Dictionaries may require recursion

		if isinstance(trg, dict) :
			for k in sorted(trg.keys()) :
				itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((str(k), str(trg[k])))))
				if isinstance(trg[k], list) or isinstance(trg[k], tuple) or isinstance(trg[k], set) or isinstance(trg[k], dict) :
					self.addTreeItem(trg[k], itms[-1])
		elif isinstance(trg, list) or isinstance(trg, tuple) or isinstance(trg, set) :
			for k in trg :
				if isinstance(k, list) or isinstance(k, tuple) or isinstance(k, set) or isinstance(k, dict) :
					try : itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((k.__class__.__name__, ""))))
					except : itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList(("??", ""))))
					self.addTreeItem(k, itms[-1])
				else :
					itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((str(k), ""))))
		else :
			itms.append(QtGui.QTreeWidgetItem(QtCore.QStringList((str(trg), ""))))

		if parent == None :
			self.wheadtree.addTopLevelItems(itms)
			self.wheadtree.resizeColumnToContents(0)
		else : parent.addChildren(itms)

#---------------------------------------------------------------------------

class EMInfoWin(QtGui.QWidget) :
	"""The info window"""

	def __init__(self, parent = None) :
		QtGui.QWidget.__init__(self, parent)

		self.target = None
		self.stack = QtGui.QStackedLayout(self)

		# We add one instance of 'infoPane' parent class to represent nothing

		self.stack.addWidget(EMInfoPane())

	def set_target(self, target, ftype) :
		"""Display the info pane for target EMDirEntry with EMFileType instance ftype"""

		self.target = target
		self.ftype = ftype

		if target == None :
			self.stack.setCurrentIndex(0)
			return

		cls = ftype.infoClass()

		for i in range(self.stack.count()) :
			if isinstance(self.stack.itemAt(i), cls) :
				self.stack.setCurrentIndex(i)
				pane = self.stack.itemAt(i)
				pane.display(target)
				break
		else :
			# If we got here, then we need to make a new instance of the appropriate pane

			if cls == None : print "No class ! (%s)"%str(ftype)
			pane = cls()
			i = self.stack.addWidget(pane)		# add the new pane and get its index
			pane.display(target)
			self.stack.setCurrentIndex(i)		# put the new pane on top

	def closeEvent(self, event) :
		QtGui.QWidget.closeEvent(self, event)
		self.emit(QtCore.SIGNAL("winclosed()"))

class SortSelTree(QtGui.QTreeView) :
	"""This is a subclass of QtGui.QTreeView. It is almost identical but implements selection processing with sorting.
	The correct way of doing this in QT4.2 is to use a QSortFilterProxy object, but that won't work properly in this case."""

	def __init__(self, parent = None) :
		QtGui.QTreeView.__init__(self, parent)
		self.header().setClickable(True)
		self.connect(self.header(), QtCore.SIGNAL("sectionClicked(int)"), self.colclick)
		self.scol = -1
		self.sdir = 1

	def setSortingEnabled(self, enable) :
#		print "enable ", enable
		return		# always enabled

	def colclick(self, col) :
		if col == self.scol : self.sdir ^= 1
		else :
			self.scol = col
			self.sdir = 1

		self.header().setSortIndicator(self.scol, self.sdir)
		self.header().setSortIndicatorShown(True)

		self.sortByColumn(self.scol, self.sdir)

	def sortByColumn(self, col, ascend) :
		if col == -1 : return

		# we mark all selected records

		try :
			for s in self.selectedIndexes() :
				if s.column() == 0 :
					s.internalPointer().sel = True
#					print s.row()
		except :
			pass
	
#			print "no model"

		# then do the actual sort

		QtGui.QTreeView.sortByColumn(self, col, ascend)

		# then set a new selection list

		sel = self.model().findSelected()
		if len(sel) == 0 :return

		qis = QtGui.QItemSelection()
		for i in sel : qis.select(i, i)
		self.selectionModel().select(qis, QtGui.QItemSelectionModel.ClearAndSelect|QtGui.QItemSelectionModel.Rows)

#		for i in sel : self.selectionModel().select(i, QtGui.QItemSelectionModel.ClearAndSelect)
#		self.update()

#---------------------------------------------------------------------------

class EMBrowserWidget(QtGui.QWidget) :
	"""This widget is a file browser for EMAN2. In addition to being a regular file browser, it supports:
	- getting information about recognized data types
	- embedding BDB: databases into the observed filesystem
	- remote database access (EMEN2)*
	"""

	def __init__(self, parent = None, withmodal = False, multiselect = False, startpath = ".", setsmode = None) :
		"""withmodal - if specified will have ok/cancel buttons, and provide a mechanism for a return value (not truly modal)
		multiselect - if True, multiple files can be simultaneously selected
		startpath - default "."
		setsmode - Used during bad particle marking
		dirregex - default "", a regular expression for filtering filenames (directory names not filtered)
		"""

		# although this looks dumb it is necessary to break Python's issue with circular imports(a major weakness of Python IMO)

		global emscene3d, emdataitem3d
		import emscene3d
		import emdataitem3d

		QtGui.QWidget.__init__(self, parent)

		self.setWindowTitle("e2display.py Browser") # Jesus

		# label = QtGui.QLabel(self);
      # label.setText("Window Title");
      # self.setWindowTitle("Window Title");

		self.withmodal = withmodal
		self.multiselect = multiselect

		self.resize(780, 580)
		self.gbl = QtGui.QGridLayout(self)

		# Top Toolbar area

		self.wtoolhbl = QtGui.QHBoxLayout()
		self.wtoolhbl.setContentsMargins(0, 0, 0, 0)

		self.wbutback = QtGui.QPushButton(QString(QChar(0x2190)))
		self.wbutback.setMaximumWidth(36)
		self.wbutback.setEnabled(False)
		self.wtoolhbl.addWidget(self.wbutback, 0)

		self.wbutfwd = QtGui.QPushButton(QString(QChar(0x2192)))
		self.wbutfwd.setMaximumWidth(36)
		self.wbutfwd.setEnabled(False)
		self.wtoolhbl.addWidget(self.wbutfwd, 0)

		# Text line for showing (or editing) full path

		self.lpath = QtGui.QLabel("  Path:")
		self.wtoolhbl.addWidget(self.lpath)

		self.wpath = QtGui.QLineEdit()
		self.wtoolhbl.addWidget(self.wpath, 5)

		# self.wspacet1 = QtGui.QSpacerItem(100, 10, QtGui.QSizePolicy.MinimumExpanding)
		# self.wtoolhbl.addSpacerItem(self.wspacet1)

		self.wbutinfo = QtGui.QPushButton("Info")
		self.wbutinfo.setCheckable(True)
		self.wtoolhbl.addWidget(self.wbutinfo, 1)

		self.gbl.addLayout(self.wtoolhbl, 0, 0, 1, 2)

		# 2nd Top Toolbar area

		self.wtoolhbl2 = QtGui.QHBoxLayout()
		self.wtoolhbl2.setContentsMargins(0, 0, 0, 0)

		self.wbutup = QtGui.QPushButton(QString(QChar(0x2191)))
		self.wbutup.setMaximumWidth(36)
		self.wtoolhbl2.addWidget(self.wbutup, 0)

		self.wbutrefresh = QtGui.QPushButton(QString(QChar(0x21ba)))
		self.wbutrefresh.setMaximumWidth(36)
		self.wtoolhbl2.addWidget(self.wbutrefresh, 0)

		# Text line for showing (or editing) full path

		self.lfilter = QtGui.QLabel("Filter:")
		self.wtoolhbl2.addWidget(self.lfilter)

		self.wfilter = QtGui.QComboBox()
		self.wfilter.setEditable(True)
		self.wfilter.setInsertPolicy(QtGui.QComboBox.InsertAtBottom)
		self.wfilter.addItem("")
		self.wfilter.addItem("(.(?!_ctf))*$")
		self.wfilter.addItem(".*\.img")
		self.wfilter.addItem(".*\.box")
		self.wfilter.addItem(".*\.hdf")
		self.wfilter.addItem(".*_ptcls$")
		self.wfilter.addItem(".*\.mrc")
		self.wfilter.addItem(".*\.tif")
		self.wfilter.addItem(".*\.pdb")
		self.wfilter.addItem("help")
		self.wtoolhbl2.addWidget(self.wfilter, 5)

		# self.wspacet1 = QtGui.QSpacerItem(100, 10, QtGui.QSizePolicy.MinimumExpanding)
		# self.wtoolhbl.addSpacerItem(self.wspacet1)

		self.selectall = QtGui.QPushButton("Sel All")
		self.wtoolhbl2.addWidget(self.selectall, 1)
		self.selectall.setEnabled(withmodal)

		self.gbl.addLayout(self.wtoolhbl2, 1, 0, 1, 2)

		### Central verticalregion has bookmarks and tree
		# Bookmarks implemented with a toolbar in a frame

		self.wbookmarkfr = QtGui.QFrame()
		self.wbookmarkfr.setFrameStyle(QtGui.QFrame.StyledPanel|QtGui.QFrame.Raised)
		self.wbmfrbl = QtGui.QVBoxLayout(self.wbookmarkfr)

		self.wbookmarks = QtGui.QToolBar()
		# self.wbookmarks.setAutoFillBackground(True)
		# self.wbookmarks.setBackgroundRole(QtGui.QPalette.Dark)
		self.wbookmarks.setOrientation(2)
		self.addBookmark("EMEN2", "emen2:")
		self.wbookmarks.addSeparator()
		self.addBookmark("SSH", "ssh:")
		self.wbookmarks.addSeparator()
		self.addBookmark("Root", "/")
		self.addBookmark("Current", os.getcwd())
		self.addBookmark("Home", e2gethome())
		self.wbmfrbl.addWidget(self.wbookmarks)

		self.gbl.addWidget(self.wbookmarkfr, 2, 0)

		# This is the main window listing files and metadata

		self.wtree = SortSelTree()

		if multiselect : self.wtree.setSelectionMode(3)	# extended selection
		else : self.wtree.setSelectionMode(1)			# single selection

		self.wtree.setSelectionBehavior(1)		# select rows
		self.wtree.setAllColumnsShowFocus(True)
		self.wtree.sortByColumn(-1, 0)			# start unsorted
		self.gbl.addWidget(self.wtree, 2, 1)

		# Lower region has buttons for actions

		self.hbl2 = QtGui.QGridLayout()

		self.wbutmisc = []

		# 10 buttons for context-dependent actions

		self.hbl2.setRowStretch(0, 1)
		self.hbl2.setRowStretch(1, 1)

		for i in range(5) :
			self.hbl2.setColumnStretch(i, 2)
	
			for j in range(2) :
				self.wbutmisc.append(QtGui.QPushButton(""))
				self.hbl2.addWidget(self.wbutmisc[-1], j, i)
				self.wbutmisc[-1].hide()
#				self.wbutmisc[-1].setEnabled(False)
				QtCore.QObject.connect(self.wbutmisc[-1], QtCore.SIGNAL('clicked(bool)'), lambda x, v = i*2+j :self.buttonMisc(v))

		self.wbutxx = QtGui.QLabel("")
		self.wbutxx.setMaximumHeight(12)
		self.hbl2.addWidget(self.wbutxx, 0, 6)
		self.wbutyy = QtGui.QLabel("")
		self.wbutyy.setMaximumHeight(12)
		self.hbl2.addWidget(self.wbutyy, 1, 6)

		# buttons for selector use

		if withmodal :
#			self.wspace1 = QtGui.QSpacerItem(100, 10, QtGui.QSizePolicy.MinimumExpanding)
#			self.hbl2.addSpacerItem(self.wspace1)

			self.wbutcancel = QtGui.QPushButton("Cancel")
			self.hbl2.addWidget(self.wbutcancel, 1, 7)

			self.wbutok = QtGui.QPushButton("OK")
			self.hbl2.addWidget(self.wbutok, 1, 8)

			self.hbl2.setColumnStretch(6, 1)
			self.hbl2.setColumnStretch(7, 1)
			self.hbl2.setColumnStretch(8, 1)

			QtCore.QObject.connect(self.wbutcancel, QtCore.SIGNAL('clicked(bool)'), self.buttonCancel)
			QtCore.QObject.connect(self.wbutok, QtCore.SIGNAL('clicked(bool)'), self.buttonOk)

		self.gbl.addLayout(self.hbl2, 4, 1)

		QtCore.QObject.connect(self.wbutback, QtCore.SIGNAL('clicked(bool)'), self.buttonBack)
		QtCore.QObject.connect(self.wbutfwd, QtCore.SIGNAL('clicked(bool)'), self.buttonFwd)
		QtCore.QObject.connect(self.wbutup, QtCore.SIGNAL('clicked(bool)'), self.buttonUp)
		QtCore.QObject.connect(self.wbutrefresh, QtCore.SIGNAL('clicked(bool)'), self.buttonRefresh)
		QtCore.QObject.connect(self.wbutinfo, QtCore.SIGNAL('clicked(bool)'), self.buttonInfo)
		QtCore.QObject.connect(self.selectall, QtCore.SIGNAL('clicked(bool)'), self.selectAll)
		QtCore.QObject.connect(self.wtree, QtCore.SIGNAL('clicked(const QModelIndex)'), self.itemSel)
		QtCore.QObject.connect(self.wtree, QtCore.SIGNAL('activated(const QModelIndex)'), self.itemActivate)
		QtCore.QObject.connect(self.wtree, QtCore.SIGNAL('doubleClicked(const QModelIndex)'), self.itemDoubleClick)
		QtCore.QObject.connect(self.wtree, QtCore.SIGNAL('expanded(const QModelIndex)'), self.itemExpand)
		QtCore.QObject.connect(self.wpath, QtCore.SIGNAL('returnPressed()'), self.editPath)
		QtCore.QObject.connect(self.wbookmarks, QtCore.SIGNAL('actionTriggered(QAction*)'), self.bookmarkPress)
		QtCore.QObject.connect(self.wfilter, QtCore.SIGNAL('currentIndexChanged(int)'), self.editFilter)

		self.setsmode = setsmode	# The sets mode is used when selecting bad particles
		self.curmodel = None	# The current data model displayed in the tree
		self.curpath = None	# The path represented by the current data model
		self.curft = None		# a fileType instance for the currently hilighted object
		self.curactions = []	# actions returned by the filtetype. Cached for speed
#		self.models = {}		# Cached models to avoid a lot of rereading (not sure if this is really worthwhile)
		self.pathstack = []	# A stack of previous viewed paths
		self.infowin = None	# The 'info' window instance

		# Each of these contains a list of open windows displaying different types of content
		# The last item in the list is considered the 'current' item for any append operations

		self.view2d = []
		self.view2ds = []
		self.view3d = []
		self.viewplot2d = []
		self.viewplot3d = []

		# These items are used to do gradually filling in of file details for better interactivity

		self.updtimer = QTimer()		# This causes the actual display updates, which can't be done from a python thread
		QtCore.QObject.connect(self.updtimer, QtCore.SIGNAL('timeout()'), self.updateDetailsDisplay)
		self.updthreadexit = False		# when set, this triggers the update thread to exit
		self.updthread = threading.Thread(target = self.updateDetails)	# The actual thread
		self.updlist = []				# List of QModelIndex items in need of updating
		self.redrawlist = []			# List of QModelIndex items in need of redisplay
		self.needresize = 0			# Used to resize column widths occaisonally
		self.expanded = set()			# We get multiple expand events for each path element, so we need to keep track of which ones we've updated

		self.setPath(startpath)	# start in the local directory
		self.updthread.start()
		self.updtimer.start(200)

		self.result = None			# used in modal mode. Holds final selection

		E2loadappwin("e2display", "main", self)

	def busy(self) :
		"""display a busy cursor"""

		QtGui.qApp.setOverrideCursor(Qt.BusyCursor)

	def notbusy(self) :
		"""normal arrow cursor"""

		QtGui.qApp.setOverrideCursor(Qt.ArrowCursor)

	def updateDetails(self) :
		"""This is spawned as a thread to gradually fill in file details in the background"""

		while 1 :
			if self.updthreadexit : break

			if len(self.updlist) == 0 :
				time.sleep(1.0)				# If there is nothing to update at the moment, we don't need to spin our wheels as much
			else :
				de = self.updlist.pop()

# 				print de.internalPointer().truepath()

				r = de.internalPointer().fillDetails()
				if r == 1 :
					self.redrawlist.append(de)		# if the update changed anything, we trigger a redisplay of this entry
					time.sleep(0.01)			# prevents updates from happening too fast and slowing the machine down
				if r == 2 :						# This means we're reading from a cache, and we should probably update as fast as possible
					self.redrawlist.append(de)
# 				print "### ", de.internalPointer().path()

	def updateDetailsDisplay(self) :
		"""Since we can't do GUI updates from a thread, this is a timer event to update the display after the beckground thread
		gets the details for each item"""

		if self.needresize > 0 :
			self.needresize -= 1
			self.wtree.resizeColumnToContents(3)
			self.wtree.resizeColumnToContents(4)

		if len(self.redrawlist) == 0 :
			return

		# we emit a datachanged event
		rdr = []
		while len(self.redrawlist) > 0 :
			i = self.redrawlist.pop()
			rdr.append((i.row(), i.internalPointer()))		# due to threads, we do it this way to make sure we don't miss any

		# We emit only a single event here for efficiency
		self.curmodel.dataChanged.emit(self.curmodel.createIndex(min(rdr)[0], 0, min(rdr)[1]), self.curmodel.createIndex(max(rdr)[0], 5, max(rdr)[1]))

		self.needresize = 2

	def editFilter(self, newfilt) :
		"""Sets a new filter. Requires reloading the current directory."""

		self.setPath(str(self.wpath.text()))

	def editPath(self) :
		"""Set a new path"""

		self.setPath(str(self.wpath.text()))

# 	def keyPressEvent(self, event) :
# 		"""Make sure we update selection when keyboard is pressed"""
#
# 		print "key", event.__dict__
# 		QtGui.QTreeView.keyPressEvent(self.wtree, event)
# 		self.itemSel(None)

	def itemSel(self, qmi) :
#		print "Item selected", qmi.row(), qmi.column(), qmi.parent(), qmi.internalPointer().path()

		qism = self.wtree.selectionModel().selectedRows()

		if len(qism) > 1 :
			self.wpath.setText("<multiple select>")

			if self.infowin != None and not self.infowin.isHidden() :
				self.infowin.set_target(None)
		elif len(qism) == 1 :
			obj = qism[0].internalPointer()
			self.wpath.setText(obj.path())
			self.curmodel.details(qism[0])
			self.wtree.resizeColumnToContents(2)
			self.wtree.resizeColumnToContents(3)

			# This makes an instance of a FileType for the selected object

			ftc = obj.fileTypeClass()
	
			if ftc != None :
				self.curft = ftc(obj.path())
				
#				if self.setsmode : self.curft.setSetsDB(re.sub(r'_ctf_flip$|_ctf_wiener$', '', obj.path()))	# If we want to enable bad particle picking (treat ctf and raw data bads as same)

				self.curactions = self.curft.actions()

#				print actions

				for i, b in enumerate(self.wbutmisc) :
					try :
						b.setText(self.curactions[i][0])
						b.setToolTip(self.curactions[i][1])
						b.show()
#						b.setEnabled(True)
					except :
						b.hide()
#						b.setEnabled(False)
			# Bug fix by JFF (What if filetype is None????)
			else :
				self.curft = None
				self.curactions = []

				for i, b in enumerate(self.wbutmisc) :
					try :
						b.hide()
					except :
						pass

			if self.infowin != None and not self.infowin.isHidden() :
				self.infowin.set_target(obj, ftc)

	def itemActivate(self, qmi) :
#		print "Item activated", qmi.row(), qmi.column()

		itm = qmi.internalPointer()

		#if itm.nChildren() > 0 :
			# self.setPath(itm.path())
		#else :
			#if self.withmodal and not self.multiselect :
				# self.buttonOk(True)
				#return

			#try :
				# self.curactions[0][2](self)
			#except :
				#pass

	def itemDoubleClick(self, qmi) :
#		print "Item activated", qmi.row(), qmi.column()

		itm = qmi.internalPointer()

		if itm.nChildren() > 0 :
			self.setPath(itm.path())
		else :
			if self.withmodal and not self.multiselect :
				self.buttonOk(True)
				return

			try :
				self.curactions[0][2](self)
			except :
				pass

	def itemExpand(self, qmi) :
		"""Called when an item is expanded"""

		if qmi.internalPointer().path() in self.expanded : return
		self.expanded.add(qmi.internalPointer().path())

#		print "expand", qmi.internalPointer().path()

		# Otherwise we get expand events on a single-click

		if qmi.internalPointer().filetype != "Folder" : return

		# we add the child items to the list needing updates

		for i in xrange(self.curmodel.rowCount(qmi)-1, -1, -1) :
			self.updlist.append(self.curmodel.index(i, 0, qmi))

	def buttonMisc(self, num) :
		"""One of the programmable action buttons was pressed"""

#		print "press ", self.curactions[num][0]

		self.curactions[num][2](self)				# This calls the action method

	def buttonOk(self, tog) :
		"""When the OK button is pressed, this will emit a signal. The receiver should call the getResult method (once) to get the list of paths"""

		qism = self.wtree.selectionModel().selectedRows()
		self.result = [i.internalPointer().path().replace(os.getcwd(), ".") for i in qism]
		self.updtimer.stop()
 		self.emit(QtCore.SIGNAL("ok")) # this signal is important when e2ctf is being used by a program running its own eve

	def buttonCancel(self, tog) :
		"""When the Cancel button is pressed, a signal is emitted, but getResult should not be called."""

		self.result = []
		self.updtimer.stop()
 		self.emit(QtCore.SIGNAL("cancel")) # this signal is important when e2ctf is being used by a program running its own eve
		self.close()

	def selectAll(self) :
		self.wtree.selectAll()

	def buttonBack(self, tog) :
		"""Button press"""

		# I don't like the stack idea, it's annoying, so I am using a circular array instead John F

		l = self.pathstack.index(self.curpath)
		self.setPath(self.pathstack[(l-1)], True)

		if l == 1 :
			self.wbutback.setEnabled(False)
			self.wbutfwd.setEnabled(True)

	def buttonFwd(self, tog) :
		"""Button press"""

		# I don't like the stack idea, it's annoying, so I am using a circular array instead John F

		l = self.pathstack.index(self.curpath)
		self.setPath(self.pathstack[(l+1)], True)

		if l == len(self.pathstack) - 2 :
			self.wbutback.setEnabled(True)
			self.wbutfwd.setEnabled(False)

	def buttonUp(self, tog) :
		"""Button press"""

		if "/" in self.curpath : newpath = self.curpath.rsplit("/", 1)[0]
		else : newpath = os.path.realpath(self.curpath).rsplit("/", 1)[0]

		print "Newpath:", newpath

		#if len(newpath) > 1 : self.setPath(newpath)	# What if we want to return to CWD, '.' # John F

		if len(newpath) > 0 : self.setPath(newpath)

	def buttonRefresh(self, tog) :
		"""Button press"""

		self.setPath(self.curpath)

	def infowinClosed(self) :
		self.wbutinfo.setChecked(False)

	def buttonInfo(self, tog) :
		if tog :
			if self.infowin == None :
				self.infowin = EMInfoWin()
				self.infowin.resize(500, 600)
			self.infowin.show()
			self.infowin.raise_()
			qism = self.wtree.selectionModel().selectedRows()
			if len(qism) == 1 :
				self.infowin.set_target(qism[0].internalPointer(), self.curft)
			else : self.infowin.set_target(None, None)
			QtCore.QObject.connect(self.infowin, QtCore.SIGNAL('winclosed()'), self.infowinClosed)
		else :
			if self.infowin != None :
				self.infowin.hide()

	def getResult(self) :
		"""In modal mode, this will return the result of a selection. Returns None before
		ok/cancel have been pressed. If cancel is pressed or ok is pressed with no selection, 
		returns an empty list []. With a valid selection, returns a list of path strings.
		When called after ok/cancel, also closes the dialog."""

		if self.result == None : return None

		self.close()

		for i in xrange(len(self.result)) :
			if self.result[i][:2] == "./" : self.result[i] = self.result[i][2:]
		return self.result

	def getCWD(self) :
		""" In modal mode, this will return the directory the browser is in. This is useful for
		using the browser to select a directory of interest. """

		# If a directory is selected, return this

		if self.result and os.path.isdir(self.result[0]) :
			self.close()
			return self.result[0]

		# If there is no current path and no result dir

		if self.curpath == None : return None

		# Return current path

		self.close()

		return self.curpath

	def addBookmark(self, label, path) :
		"""Add a new bookmark"""

		act = self.wbookmarks.addAction(label)
		act.setData(path)

	def setPath(self, path, silent = False, inimodel = EMFileItemModel) :
		"""Sets the current root path for the browser window. If silent is true, 
		the path history will not be updated."""

		if path != ""  and  path[0] == ":" :
			os.system(path[1:])
			return

		path = expand_string(path)

		if path == "" :
			path = "."

		path = path.replace("\\", "/")
		if path[:2] == "./" : path = path[2:]

		self.updlist = []
		self.redrawlist = []

		filt = str(self.wfilter.currentText()).strip()

		if filt != ""  and  not os.path.isdir(path) :
			path = os.path.dirname(path)
			if path == "" :
				path = "."

		self.curpath = str(path)
		self.wpath.setText(path)

		if filt == "" :
			filt = None
		elif filt == "?"  or  lower(filt) == "help" :
			filt = None
			hlp  = \
			"Enter a regular expression to filter files to see, or\n" + \
			"enter wildcard file name patterns as in Linux 'ls' or DOS 'dir',\n" + \
			"where you may use wildcards *, ?, %, #, or !.\n" + \
			"    * matches any number of characters, including none\n" + \
			"    ? matches 1 character, or none\n" + \
			"    % matches 1 character\n" + \
			"    # matches 1 decimal digit, 0 to 9\n" + \
			"    ! matches 1 upper or lower case letter\n" + \
			"The syntax is\n" + \
			"    { pattern } [ 'not' { pattern } ]\n" + \
			"Any patterns after the optional 'not' exclude matching files.\n" + \
			"Folders have an implicit '.dir' extension, to match or exclude.\n" + \
			"If a wildcard pattern fails to work, it is probably interpreted\n" + \
			"as a regular expression, so append ' [' to it to prevent it.\n" + \
			"Examples:\n" + \
			"1. *.dir         - find all directories (folders)\n" + \
			"2. * not *.dir   - find all non-directories\n" + \
			"3. *.txt *.tiff  - find all text files or tiff files\n" + \
			"4. *             - find all files"

			QtGui.QMessageBox.warning(None, "Info", hlp)
		else :
			try :
				flt = re.compile(filt)
				filt = flt
			except :
				filt = filt
#				QtGui.QMessageBox.warning(self, "Error", "Bad filter expression")

		# if path in self.models :
			# self.curmodel = self.models[path]
		# else :
			# self.curmodel = inimodel(path)
			# self.models[self.curpath] = self.curmodel

		if filt != None and filt != "" :
			try :
				self.curmodel = inimodel(path, dirregex = filt)
			except :
				self.curmodel = inimodel(path)
 				filt = None
#				QtGui.QMessageBox.warning(None, "Error", "Filtering not allowed.")
				print "Filtering is not implemented in this instance of the file browser."
		else :
			self.curmodel = inimodel(path)

		self.wtree.setSortingEnabled(False)
		self.wtree.setModel(self.curmodel)
		self.wtree.setSortingEnabled(True)
		self.wtree.resizeColumnToContents(0)
		self.wtree.resizeColumnToContents(1)
		self.wtree.resizeColumnToContents(3)
		self.wtree.resizeColumnToContents(4)

		self.expanded = set()

		# we add the child items to the list needing updates

		for i in xrange(self.curmodel.rowCount(None)-1, -1, -1) :
			self.updlist.append(self.curmodel.index(i, 0, None))

		if not silent :
			try : self.pathstack.remove(self.curpath)
			except : pass
			self.pathstack.append(self.curpath)
			if len(self.pathstack) > 1 :
				self.wbutback.setEnabled(True)
			else :
				self.wbutback.setEnabled(False)
				self.wbutfwd.setEnabled(False)

	def bookmarkPress(self, action) :
		""""""

#		print "Got action ", action.text(), action.data().toString()

		self.setPath(action.data().toString())
#		self.wtree.setSelectionModel(myQItemSelection(self.curmodel))

	def closeEvent(self, event) :
		E2saveappwin("e2display", "main", self)
		self.updthreadexit = True

		for w in self.view2d+self.view2ds+self.view3d+self.viewplot2d+self.viewplot3d :
			w.close()

		if self.infowin != None :
			self.infowin.close()

		event.accept()
		self.updtimer.stop()
		# self.app().close_specific(self)
		self.emit(QtCore.SIGNAL("module_closed"))

# This is just for testing, of course

def test_result() :
	global window
	print "Returned"
	print window.getResult()

if __name__ == '__main__' :
	em_app = EMApp()
	window = EMBrowserWidget(withmodal = True, multiselect = True)
	QtCore.QObject.connect(window, QtCore.SIGNAL("ok"), test_result)
	QtCore.QObject.connect(window, QtCore.SIGNAL("cancel"), test_result)

	window.show()
	ret = em_app.exec_()
	try : window.updthreadexit = True
	except : pass
	sys.exit(ret)
