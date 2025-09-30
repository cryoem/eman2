###!/usr/bin/env python
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

from past.utils import old_div
from builtins import range
from builtins import object
from EMAN2 import *
from EMAN2jsondb import js_open_dict
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt, QTimer
from .emapplication import EMApp
from .emimage2d import *
from .emimagemx import *
from .empdbitem3d import *
from .emplot2d import *
from .emhist import *
from .emplot3d import *
from libpyUtils2 import EMUtil
from .matching import matches_pats
from .valslider import StringBox
import os
import re
import threading
import time
import traceback
import weakref
import random


def display_error(msg) :
	"""Displays an error message, in gui and on terminal."""
	print(msg)
	sys.stdout.flush()
	QtWidgets.QMessageBox.warning(None, "Error", msg)

# We need to sort ints and floats as themselves, not string, John Flanagan
def safe_int(v) :
	"""Performs a safe conversion from a string to an int. If a non int is presented we return the lowest possible value"""
	try :
		return int(v)
	except (ValueError, TypeError) :
		return -sys.maxsize-1

def safe_float(v) :
	"""Performs a safe conversion from a string to a float. If a non float is presented we return the lowest possible value"""
	try :
		return float(v)
	except :
		try : return float(v.split()[0])
		except :
			return sys.float_info.min

def size_sortable(s):
	try: return int(s)
	except: return 0
#	try:
	#if s.endswith("g"): return int(s[:-1])*1000000000
	#elif s.endswith("m"): return int(s[:-1])*1000000
	#elif s.endswith("k"): return int(s[:-1])*1000
	#else: return int(s)
#	except: return 0

def isprint(s) :
	"""returns True if the string contains only printable ascii characters"""
	# Seems like no isprint() in python, this does basically the same thing
	try: s=s.decode("utf-8")
	except: return False
	mpd = s.translate("AAAAAAAAABBAABAAAAAAAAAAAAAAAAAA"+"B"*95+"A"*129)

	if "A" in mpd :
		ind = mpd.index("A")
#		print "bad chr %d at %d"%(ord(s[ind]), ind)
		return False

	return True

def askFileExists() :
	"""Opens a dialog and asks the user what to do if a file to be written to already exists"""
	box = QtWidgets.QMessageBox(4, "File Exists", "File already exists. What would you like to do ?")	# 4 means we're asking a question
	b1 = box.addButton("Append", QtWidgets.QMessageBox.AcceptRole)
	b2 = box.addButton("Overwrite", QtWidgets.QMessageBox.AcceptRole)
	b3 = box.addButton("Cancel", QtWidgets.QMessageBox.AcceptRole)

	box.exec_()

	if box.clickedButton() == b1 : return "append"
	elif box.clickedButton() == b2 : return "overwrite"
	else : return "cancel"

def makeOrthoProj(ptcl,layers,center,highpass,lowpass,stkout):
	"""makes restricted orthogonal projections from a 3-D volume and returns as a single 2d image
	ptcl - 3D input ovlume
	layers - +- layer range about center for integration, -1 is full projection
	center - center location, 0 being the middle of the box
	highpass - optional high-pass filter in A (<0 disables)
	lowpass - optional low-pass filter in A (<0 disables)
	"""
	
	# these are the range limited orthogonal projections, with optional filtration
	if layers>=0:
		first=ptcl["nx"]/2+center-layers
		last=ptcl["nx"]/2+center+layers+1
	else:
		first=0
		last=-1
	x=ptcl.process("misc.directional_sum",{"axis":"x","first":first,"last":last})
	if lowpass>0 : x.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/lowpass})
	if highpass>0 : x.process_inplace("filter.highpass.gauss",{"cutoff_freq":1.0/highpass})
	get_application().processEvents()	# keeps the GUI happy

	y=ptcl.process("misc.directional_sum",{"axis":"y","first":first,"last":last})
	if lowpass>0 : y.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/lowpass})
	if highpass>0 : y.process_inplace("filter.highpass.gauss",{"cutoff_freq":1.0/highpass})
	get_application().processEvents()

	z=ptcl.process("misc.directional_sum",{"axis":"z","first":first,"last":last})
	if lowpass>0 : z.process_inplace("filter.lowpass.gauss",{"cutoff_freq":1.0/lowpass})
	if highpass>0 : z.process_inplace("filter.highpass.gauss",{"cutoff_freq":1.0/highpass})

	# different directions sometimes have vastly different standard deviations, independent normalization may help balance
	x.process_inplace("normalize")
	y.process_inplace("normalize")
	z.process_inplace("normalize")
	
	get_application().processEvents()
	
	if stkout:
		hall=[x,y,z]
		for h in hall:
			h.set_attr_dict(ptcl.get_attr_dict())
	else:
		# we pack the 3 projections into a single 2D image
		hall=EMData(x["nx"]*3,x["ny"],1)
		hall.insert_clip(x,(0,0))
		hall.insert_clip(y,(x["nx"],0))
		hall.insert_clip(z,(x["nx"]*2,0))
		hall.set_attr_dict(ptcl.get_attr_dict())

	return hall

class EMFileType(object) :
	"""This is a base class for handling interaction with files of different types. It includes a number of execution methods common to
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
		"""Returns (size, n, dim) if the referenced path is a file of this type, false if not valid. The first 4k block of data from the file is provided as well to avoid unnecessary file access."""
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
		outpath = QtWidgets.QInputDialog.getText(None, "Save Filename", "Filename to save to (type determined by extension)", 0, self.path)

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
		else : ns = list(range(EMUtil.get_image_count(self.path)))

		for i in ns :
			im = EMData(self.path, i)
			im.write_image(outpath, -1)

		brws.notbusy()

	def plot2dApp(self, brws) :
		"""Append self to current plot"""
		brws.busy()

		if self.n <= 0 : data = EMData(self.path)
		else : data = EMData(self.path, self.n)

		try :
			target = brws.viewplot2d[-1]
			if target.closed : 
				brws.viewplot2d.pop()
				raise Exception
			target.set_data(data, display_path(self.path))
		except :
			target = EMPlot2DWidget()
			brws.viewplot2d.append(target)
			target.set_data(data, display_path(self.path))

		target.qt_parent.setWindowTitle(display_path(self.path))

		brws.notbusy()
		target.show()
		target.raise_()

	def plot2dNew(self, brws) :
		"""Make a new plot"""
		brws.busy()

		if self.n <= 0 : data = EMData(self.path)
		else : data = EMData(self.path, self.n)

		target = EMPlot2DWidget()
		brws.viewplot2d.append(target)
		target.set_data(data, display_path(self.path))

		target.qt_parent.setWindowTitle(display_path(self.path))

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

		target.insertNewNode(display_path(self.path), data, parentnode = target)
		iso = emdataitem3d.EMIsosurface(data)
		self.loadIsoColor(iso)
		target.insertNewNode('Isosurface', iso, parentnode = data)
		target.initialViewportDims(data.getData().get_xsize())	# Scale viewport to object size
		target.setCurrentSelection(iso)				# Set isosurface to display upon inspector loading
		target.updateSG()	# this is needed because this might just be an addition to the SG rather than initialization

		target.setWindowTitle(display_path(self.path))

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

		target.insertNewNode(display_path(self.path), data)
		iso = emdataitem3d.EMIsosurface(data)
		self.loadIsoColor(iso)
		target.insertNewNode('Isosurface', iso, parentnode = data)
		target.initialViewportDims(data.getData().get_xsize())	# Scale viewport to object size
		target.setCurrentSelection(iso)				# Set isosurface to display upon inspector loading
		
		brws.notbusy()
		target.setWindowTitle(display_path(self.path))

		target.show()
		target.raise_()

	def loadIsoColor(self,iso):
		"""Common routine to look for a fscvol file we can associate with the loaded volume to use as a colormap"""
		if not E2getappval("e2display","isosurface_autocolor_res",True) : return
	
		base,fsp=os.path.split(self.path)
		if base=="": base="."
		n=fsp.rsplit("_",1)[-1][:-4]
		if not n.isdigit(): return
		resfsp=f"{base}/fscvol_{n}.hdf"
		if resfsp==self.path or fsp[:6]=="fscvol": return
		if not os.path.isfile(resfsp): 
			#print("doesn't exist",resfsp)
			return
		print("loading ",resfsp)
		iso.setCmapData(resfsp)
		iso.setRGBmode(2)

	def show3DAll(self, brws) :
		"""All in new 3-D window (3-D stacks)"""
		brws.busy()

		target = emscene3d.EMScene3D()
		brws.view3d.append(target)

		for n in range(self.nimg) :
			data = emdataitem3d.EMDataItem3D(self.path, n = n)
			target.insertNewNode("{} #{}".format(display_path(self.path), n), data)
			iso = emdataitem3d.EMIsosurface(data)
#			self.loadIsoColor(iso)
			target.insertNewNode('Isosurface', iso, parentnode = data)

		target.initialViewportDims(data.getData().get_xsize())	# Scale viewport to object size
		target.setCurrentSelection(iso)				# Set isosurface to display upon inspector loading
		brws.notbusy()
		target.setWindowTitle(display_path(self.path))

		target.show()
		target.raise_()

	def show2dAvg(self, brws, new=False) :
		"""Show the unaligned average of all images in a stack"""
		brws.busy()

		# this averages images in chunks of 1000 (when more than 1000 images in the file)
		avg=sum([sum(EMData.read_images(self.path,list(range(i,min(i+1000,self.nimg))))) for i in range(0,self.nimg,1000)])

		if not new:
			try :
				target = brws.view2d[-1]
				if target.closed : 
					brws.view2d.pop()
					raise Exception

				target.set_data(avg)
			except :
				new=True
				
		if new:
			target = EMImage2DWidget()
			target.set_data(avg)
			brws.view2d.append(target)

		target.setWindowTitle("sum of"+display_path(self.path))

		brws.notbusy()
		target.show()
		target.raise_()

	def show2dAvgRnd(self, brws) :
		"""Show averages of random subsets of 1/4 or 1000 images from the stack"""

		# minimum number of images required
		if self.nimg<8 : return
		brws.busy()

		avgs=[]
		for i in range(10):
			imns=[random.randrange(self.nimg) for j in range(min(self.nimg//4,1000))]
			imns.sort()
			ims=EMData.read_images(self.path,imns)
			avgs.append(sum(ims))
			avgs[-1].mult(1.0/len(imns))

		try :
			target = brws.view2d[-1]
			if target.closed : 
				brws.view2d.pop()
				raise Exception

			target.set_data(avgs)
			#if self.getSetsDB() : target.set_single_active_set(self.getSetsDB())
		except :
			target = EMImage2DWidget()
			target.set_data(avgs)
#			target.mx_image_double.connect(target.mouse_double_click)		# this makes class average viewing work in app mode
			brws.view2ds.append(target)

		target.qt_parent.setWindowTitle("Random Avg Stack - "+display_path(self.path))

		brws.notbusy()
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
			if target.closed : 
				brws.view2ds.pop()
				raise Exception
			target.set_data(self.path, self.path)
			#if self.getSetsDB() : target.set_single_active_set(self.getSetsDB())
		except :
			target = EMImageMXWidget()
			target.set_data(self.path, self.path)
			target.mx_image_double.connect(target.mouse_double_click)		# this makes class average viewing work in app mode
			# if self.getSetsDB() : target.set_single_active_set(self.getSetsDB())
			brws.view2ds.append(target)

		target.qt_parent.setWindowTitle("Stack - "+display_path(self.path))

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
		target.mx_image_double.connect(target.mouse_double_click)
		# if self.getSetsDB() : target.set_single_active_set(self.getSetsDB())
		brws.view2ds.append(target)

		target.qt_parent.setWindowTitle("Stack - "+display_path(self.path))

		brws.notbusy()
		target.show()
		target.raise_()

	def show2dStack3z(self, brws) :
		"""A set of 2-D images derived from a stack of 3-D Volumes"""
		brws.busy()

		# if self.dim[2] > 1 :
			# data = []
			# for z in range(self.dim[2]) :
				# data.append(EMData(self.path, 0, False, Region(0, 0, z, self.dim[0], self.dim[1], 1)))
		# else : data = EMData.read_images(self.path)
		modifiers = QtWidgets.QApplication.keyboardModifiers()
		if modifiers == QtCore.Qt.ShiftModifier:
			#print("rotate x")
			xyz="y"
		elif modifiers == QtCore.Qt.ControlModifier:
			#print("rotate y")
			xyz="x"
		else:
			xyz="z"

		data=[EMData(self.path,i).process("misc.directional_sum",{"axis":xyz}) for i in range(self.nimg)]

		try :
			target = brws.view2ds[-1]
			if target.closed : 
				brws.view2ds.pop()
				raise Exception
			target.set_data(data,self.path)
			#if self.getSetsDB() : target.set_single_active_set(self.getSetsDB())
		except :
			target = EMImageMXWidget()
			target.set_data(data,self.path)
			#target.mx_image_double.connect(target.mouse_double_click)		# this makes class average viewing work in app mode
			# if self.getSetsDB() : target.set_single_active_set(self.getSetsDB())
			brws.view2ds.append(target)

		target.qt_parent.setWindowTitle("Stack - "+display_path(self.path))

		brws.notbusy()
		target.show()
		target.raise_()

	def show2dStack3sec(self, brws) :
		"""A set of 2-D images derived from a stack of 3-D Volumes"""
		
		modifiers = QtWidgets.QApplication.keyboardModifiers()
		
		if modifiers == QtCore.Qt.ShiftModifier:
			self.showProjXYZ(brws)
			return
		
		try:
			ret=self.secparm.exec_()
		except:
			self.secparm=EMSliceParamDialog(brws,self.nimg)
			ret=self.secparm.exec_()
		
		if not ret: return	# cancel
	
		# these won't be available if nimg==1
		try:
			img0=self.secparm.wspinmin.value()
			img0=max(img0,0)
			img1=self.secparm.wspinmax.value()
			if img1<0 or img1>self.nimg: img1=self.nimg
			imgstep=self.secparm.wspinstep.value()
		except:
			img0=0
			img1=1
			imgstep=1
		layers=self.secparm.wspinlayers.value()
		center=self.secparm.wspincenter.value()
		applyxf=self.secparm.wcheckxf.checkState()
		applysym=str(self.secparm.wlesym.text())
		stkout=self.secparm.wcheckstk.checkState()
		oldwin=self.secparm.wcheckoldwin.checkState()
		highpass=float(self.secparm.wlehp.text())
		lowpass=float(self.secparm.wlelp.text())
		
		maskfsp=str(self.secparm.wlemask.text())
		mask=None
		if len(maskfsp)>4 :
			try: mask=EMData(maskfsp,0)
			except: print("ERROR: unable to read mask file: ",maskfsp)

		reffsp=str(self.secparm.wleref.text())
		ref=None
		if len(reffsp)>4 :
			try: ref=EMData(reffsp,0)
			except: print("ERROR: unable to read ref file: ",maskfsp)
		
		brws.busy()

		# if self.dim[2] > 1 :
			# data = []
			# for z in range(self.dim[2]) :
				# data.append(EMData(self.path, 0, False, Region(0, 0, z, self.dim[0], self.dim[1], 1)))
		# else : data = EMData.read_images(self.path)
				
		progress = QtWidgets.QProgressDialog("Generating projections", "Cancel", 0, img1-img0-1,None)
#		progress.show()
		
		# make an empty array so we can easily interleave
		data=[None for i in range(img0,img1,imgstep)]
		nd=len(data)
		if stkout: data=data*3
		
		# reference goes last!
		if ref!=None:
			if mask!=None: ref.mult(mask)
			hall=makeOrthoProj(ref,layers,center,highpass,lowpass,stkout)
			if isinstance(hall,EMData): data.append(hall)
			else: data.extend(hall)
		
		c=0
		for i in range(img0,img1,imgstep):
			try: ptcl=EMData(self.path,i)
			except:
				print(f"Error reading {self.path} {i}")
				continue
			
			try: xf=ptcl["xform.align3d"]
			except: xf=Transform()
			
			if applyxf: ptcl.process_inplace("xform",{"transform":xf})
			if applysym is not None and applysym!="" and applysym.lower()!="c1": ptcl.process_inplace("xform.applysym",{"sym":applysym})
			if mask!=None : ptcl.mult(mask)
			
			time.sleep(0.001)
			get_application().processEvents()
			
			hall=makeOrthoProj(ptcl,layers,center,highpass,lowpass,stkout)
			
			if isinstance(hall,EMData):
				for k in ["ptcl_repr","class_ptcl_idxs","class_ptlc_src","orig_file","orig_n","source_path","source_n"]:
					try: hall[k]=ptcl[k]
					except: pass
				data[c]=hall
			else:
				data[c]=hall[0]
				data[c+nd]=hall[1]
				data[c+2*nd]=hall[2]
			
			c+=1
			progress.setValue(i-img0)
			
			if progress.wasCanceled():
#				progress.close()
				return

		if self.nimg==1 or stkout:
			if oldwin : 
				try: 
					target=brws.view2d[-1]
					if target.closed : 
						brws.view2d.pop()
						raise Exception
				except:
					target = EMImage2DWidget()
					brws.view2d.append(target)
				old=target.get_data(True)
				if isinstance(old,list) : old.extend(data)
				else: 
					if old!=None: print("embrowser.py error: old is type ",type(old))
					old=data
				target.set_data(old,self.path)
			else:
				target = EMImage2DWidget()
				target.set_data(data,self.path)
				brws.view2d.append(target)
		else:
			target = EMImageMXWidget()
			target.set_data(data,self.path)
			target.mx_image_double.connect(target.mouse_double_click)		# this makes class average viewing work in app mode
			brws.view2ds.append(target)

		target.qt_parent.setWindowTitle("Stack - "+display_path(self.path))

		brws.notbusy()
		target.show()
		target.raise_()

	def show2dSingle(self, brws, new=False) :
		"""Show a single 2-D image"""
		brws.busy()

		if self.nimg > 1 :
			if self.n >= 0 : data = EMData(self.path, self.n)
			else : data = EMData.read_images(self.path)
		else : data = EMData(self.path)
		
		#### allow view from x/y/z axis
		modifiers = QtWidgets.QApplication.keyboardModifiers()
		
		if modifiers == QtCore.Qt.ShiftModifier:
			#print("rotate x")
			xyz=0
		elif modifiers == QtCore.Qt.ControlModifier:
			#print("rotate y")
			xyz=1
		else:
			xyz=-1
			
		if new==False:

			try :
				target = brws.view2d[-1]
				if target.closed : 
					brws.view2d.pop()
					raise Exception
				target.set_data(data, xyz=xyz)
			except :
				new=True
		
		if new:
			target = EMImage2DWidget()
			target.set_data(data, xyz=xyz)
			brws.view2d.append(target)

		target.setWindowTitle(display_path(self.path))

		brws.notbusy()
		target.show()
		target.raise_()

	def show2dSingle30(self, brws) :
		"""Show a single 3-D volume as a 2-D image"""
		brws.busy()

		data = EMData(self.path,0)

		try :
			target = brws.view2d[-1]
			if target.closed : 
				brws.view2d.pop()
				raise Exception
			target.set_data(data)
		except :
			target = EMImage2DWidget(data)
			brws.view2d.append(target)

		target.setWindowTitle(display_path(self.path))

		brws.notbusy()
		target.show()
		target.raise_()
	
	def show2dSingle31(self, brws) :
		"""Show a single 3-D volume as a 2-D image"""
		brws.busy()
		#print(self.path)

		data = EMData(self.path,1)

		try :
			target = brws.view2d[-1]
			if target.closed : 
				brws.view2d.pop()
				raise Exception
			target.set_data(data)
		except :
			target = EMImage2DWidget(data)
			brws.view2d.append(target)

		target.setWindowTitle(display_path(self.path))

		brws.notbusy()
		target.show()
		target.raise_()

	def show2dSingleNew(self, brws) :
		"""Show a single 2-D image"""
		self.show2dSingle(brws, new=True)

		#brws.busy()

		#if self.nimg > 1 :
			#if self.n >= 0 : data = EMData(self.path, self.n)
			#else : data = EMData.read_images(self.path)
		#else : data = EMData(self.path)

		#modifiers = QtWidgets.QApplication.keyboardModifiers()
		#if modifiers == QtCore.Qt.ShiftModifier:
			#print("rotate x")
			#target = EMImage2DWidget()
			#target.set_data(data, xyz=0)
		#if modifiers == QtCore.Qt.ControlModifier:
			#print("rotate y")
			#target = EMImage2DWidget()
			#target.set_data(data, xyz=1)
		#else:
			#target = EMImage2DWidget(data)
		#brws.view2d.append(target)

		#target.qt_parent.setWindowTitle(display_path(self.path))

		#brws.notbusy()
		#target.show()
		#target.raise_()

	def showFilterTool(self, brws) :
		"""Open in e2filtertool.py"""
		modifiers = QtWidgets.QApplication.keyboardModifiers()
		cmd="e2filtertool.py {}".format(self.path)
		if modifiers == QtCore.Qt.ShiftModifier:
			#print("Running filter tool in safe mode...")
			cmd+=" --safemode"
		
		if self.n>=0:
			cmd+=" --idx {:d} ".format(self.n)
			
		os.system(cmd+"&");


class EMTextFileType(EMFileType) :
	"""FileType for files containing normal ASCII text"""

	@staticmethod
	def name() :
		"""The unique name of this FileType. Stored in EMDirEntry.filetype for each file."""
		return "Text"

	@staticmethod
	def isValid(path, header) :
		"""Returns (size, n, dim) if the referenced path is a file of this type, None if not valid. The first 4k block of data from the file is provided as well to avoid unnecessary file access."""
#		if not isprint(header) : return False			# demand printable Ascii. FIXME: what about unicode ?
		try: s=header.decode("utf-8")
		except: return False

		try : size = os.stat(path)[6]
		except : return False

		if size > 5000000 : dim = "big"
		else :
			f = open(path, "r").read()
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
		"""Returns (size, n, dim) if the referenced path is a file of this type, None if not valid. The first 4k block of data from the file is provided as well to avoid unnecessary file access."""

		try: s=header.decode("utf-8")
		except: return False
#		if not isprint(header) : return False			# demand printable Ascii. FIXME: what about unicode ?
		if isinstance(header, bytes):
			header=header.decode("utf-8")
		if not "<html>" in header.lower() : return False # For the moment, we demand an <html> tag somewhere in the first 4k

		try :
			size = os.stat(path)[6]
		except : return False

		if size > 5000000 : dim = "big"
		else :
			f = open(path, "r").read()
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
			print("Sorry, I don't know how to run Firefox on this platform")


class EMPDFFileType(EMFileType) :
	"""FileType for PDF files"""

	@staticmethod
	def name() :
		"""The unique name of this FileType. Stored in EMDirEntry.filetype for each file."""
		return "PDF"

	@staticmethod
	def isValid(path, header) :
		"""Returns (size, n, dim) if the referenced path is a file of this type, None if not valid. The first 4k block of data from the file is provided as well to avoid unnecessary file access."""
		if header[:4]!=b"%PDF": return False

		try :
			size = os.stat(path)[6]
		except : return False

		return (f"{size}", "-", "-")

	@staticmethod
	def infoClass() :
		"""Returns a reference to the QWidget subclass for displaying information about this file"""
		return EMPDFInfoPane

	def actions(self) :
		"""Few actions"""
		if get_platform() == "Linux" :
			return [("gv", "Open using gv", self.showPlatform),("Firefox", "Open in Firefox", self.showFirefox)]
		elif get_platform() == "Darwin" :
			return [("Open", "Default Open", self.showPlatform)]
		else : 
			return []

	def showPlatform(self, brws) :
		"""Try to open file in firefox"""
		if get_platform() == "Linux" :
			os.system("gv %s"%os.path.abspath(self.path))
		elif get_platform() == "Darwin" :
			os.system("open {}".format(os.path.abspath(self.path)))		# uses the default browser
		else : 
			print("Sorry, I don't know how to display PDF files on this platform")

	def showFirefox(self, brws) :
		"""Try to open file in firefox"""
		if get_platform() == "Linux" :
			os.system("firefox -new-tab file://%s"%os.path.abspath(self.path))
		elif get_platform() == "Darwin" :
			os.system("open {}".format(os.path.abspath(self.path)))		# uses the default browser
		else : 
			print("Sorry, I don't know how to run Firefox on this platform")


class EMPlotFileType(EMFileType) :
	"""FileType for files containing normal ASCII text"""

	@staticmethod
	def name() :
		"""The unique name of this FileType. Stored in EMDirEntry.filetype for each file."""
		return "Plot"

	@staticmethod
	def isValid(path, header) :
		"""Returns (size, n, dim) if the referenced path is a file of this type, None if not valid. The first 4k block of data from the file is provided as well to avoid unnecessary file access."""
		try: s=header.decode("utf-8")
		except: return False

		# We need to try to count the columns in the file
		header=header.decode("utf-8")
		hdr = header.splitlines()
		numc = 0
		#print(path)
		for l in hdr :
			if len(l)==0: continue
			if l[0] == "#" : continue		# comment lines ok
			if len(l)>4000:
				return (os.stat(path)[6], '-', 'big')
			numc=renumfind.findall(l)

			try:
				numc = len([float(i) for i in numc])		# number of numeric columns
			except:
				numc=0
				break

			if numc > 0 : break			# just finding the number of columns

		# If we couldn't find at least one valid line with some numbers, give up
		if numc == 0 : return False

		size = os.stat(path)[6]

		# Make sure all of the lines have the same number of columns
		fin = open(path, "r")
		numr = 0

		for l in fin :
			if l[0] == "#" or len(l) < 2 or "nan" in l : continue
			try: lnumc = len([float(i) for i in renumfind.findall(l)])
			except: continue
			if lnumc != 0 and lnumc != numc : return False				# 0 means the line contains no numbers, we'll live with that, but if there are numbers, it needs to match
			if lnumc != 0 : numr += 1
		
		if numr<3 and numc<3: return False

		return (size, "-", "%d x %d"%(numr, numc))

	@staticmethod
	def infoClass() :
		"""Returns a reference to the QWidget subclass for displaying information about this file"""
		return EMPlotInfoPane

	def __init__(self, path) :
		if path[:2] == "./" : path = path[2:]
		EMFileType.__init__(self, path)	# the current path this FileType is representing

		# Make sure all of the lines have the same number of columns
		fin = open(path, "r")
		numr = 0
		numc = 0

		for l in fin :
			if "nan" in l :
				print("Warning, NaN present in file")
				continue

			if l[0] == "#" : continue
			lnumc = len([float(i) for i in renumfind.findall(l)])

			if lnumc != 0 and numc == 0 : numc = lnumc
			elif lnumc != 0 and lnumc != numc :
				print("Error: invalid Plot file:", path)
				self.numr = 0
				self.numc = 0
				return
			elif lnumc != 0 : numr += 1

			self.numr = numr
			self.numc = numc

	def actions(self) :
		"""Returns a list of (name, help, callback) tuples detailing the operations the user can call on the current file.
		callbacks will also be passed a reference to the browser object."""

		if self.numc > 2 : rtr=[("Plot 2D", "Add to current plot", self.plot2dApp), ("Plot 2D+", "Make new plot", self.plot2dNew),
			("Histogram", "Add to current histogram", self.histApp),("Histogram +", "Make new histogram", self.histNew),
			("Plot 3D", "Add to current 3-D plot", self.plot3dApp),("Plot 3D+", "Make new 3-D plot", self.plot3dNew)]
		else: rtr=[("Plot 2D", "Add to current plot", self.plot2dApp), ("Plot 2D+", "Make new plot", self.plot2dNew),("Hist 2D", "Add to current histogram", self.histApp),("Hist 2D+", "Make new histogram", self.histNew)]
		if self.numc>=3 and self.numc<=5: rtr.append(("Spheres","Each X line is X-Y-Z[-A[-S]]. Show as spheres in 3-D",self.showSpheres))

		return rtr

	def showSpheres(self,brws):
		"""New 3-D window"""
		brws.busy()

		target = emscene3d.EMScene3D()
		brws.view3d.append(target)

		data=np.loadtxt(self.path).transpose()
		gaussplots=[emshapeitem3d.EMScatterPlot3D()]

		target.insertNewNode(self.path.split("/")[-1],gaussplots[0])
		data[:3]*=256.0
		gaussplots[0].setData(data)

		# target.initialViewportDims(data.getData().get_xsize())	# Scale viewport to object size
		# target.setCurrentSelection(iso)				# Set isosurface to display upon inspector loading

		brws.notbusy()
		target.setWindowTitle(display_path(self.path))

		target.show()
		target.raise_()

	def plot2dApp(self, brws) :
		"""Append self to current plot"""
		brws.busy()

		#data1 = []
		#fin = file(self.path, "r")
		#numr = 0

		#for l in fin :
			#if l[0] == "#" or "nan" in l : continue
			#data1.append([float(i) for i in renumfind.findall(l)])

		#data = []

		#for c in xrange(self.numc) :
			#data.append([i[c] for i in data1])

		try :
			target = brws.viewplot2d[-1]
			if target.closed : 
				brws.viewplot2d.pop()
				raise Exception
			target.set_data_from_file(self.path)
			#target.set_data(data, remove_directories_from_name(self.path, 1))
		except :
			target = EMPlot2DWidget()
			brws.viewplot2d.append(target)
			target.set_data_from_file(self.path)

		target.qt_parent.setWindowTitle(display_path(self.path))

		brws.notbusy()
		target.show()
		target.raise_()

	def plot2dNew(self, brws) :
		"""Make a new plot"""
		brws.busy()

		#data1 = []
		#fin = file(self.path, "r")
		#numr = 0

		#for l in fin :
			#if l[0] == "#" or "nan" in l : continue
			#data1.append([float(i) for i in renumfind.findall(l)])

		#data = []

		#for c in xrange(self.numc) :
			#data.append([i[c] for i in data1])

		target = EMPlot2DWidget()
		brws.viewplot2d.append(target)
		target.set_data_from_file(self.path)
		#target.set_data(data, remove_directories_from_name(self.path, 1))

		target.qt_parent.setWindowTitle(display_path(self.path))

		brws.notbusy()
		target.show()
		target.raise_()
	
	def histNew(self, brws) :
		"""Make a new plot"""
		brws.busy()

		data1 = []
		fin = open(self.path, "r")
		numr = 0

		for l in fin :
			if l[0] == "#" or "nan" in l : continue
			data1.append([float(i) for i in renumfind.findall(l)])

		data = []

		for c in range(self.numc) :
			data.append([i[c] for i in data1])

		target = EMHistogramWidget()
		brws.viewhist.append(target)
		target.set_data(data, remove_directories_from_name(self.path, 1))

		target.qt_parent.setWindowTitle(display_path(self.path))

		brws.notbusy()
		target.show()
		target.raise_()

	def histApp(self, brws) :
		"""Append self to current plot"""
		brws.busy()

		data1 = []
		fin = open(self.path, "r")
		numr = 0

		for l in fin :
			if l[0] == "#" or "nan" in l : continue
			data1.append([float(i) for i in renumfind.findall(l)])

		data = []

		for c in range(self.numc) :
			data.append([i[c] for i in data1])

		try :
			target = brws.viewhist[-1]
			if target.closed : 
				brws.viewhist.pop()
				raise Exception
			target.set_data(data, remove_directories_from_name(self.path, 1))
		except :
			target = EMHistogramWidget()
			brws.viewhist.append(target)
			target.set_data(data, remove_directories_from_name(self.path, 1))
		
		target.qt_parent.setWindowTitle(display_path(self.path))

		brws.notbusy()
		target.show()
		target.raise_()

	def plot3dApp(self, brws) :
		"""Append self to current plot"""
		brws.busy()

		try :
			target = brws.viewplot3d[-1]
			if target.closed : 
				brws.viewplot3d.pop()
				raise Exception
			target.set_data_from_file(self.path)
			#target.set_data(data, remove_directories_from_name(self.path, 1))
		except :
			target = EMPlot3DWidget()
			brws.viewplot3d.append(target)
			target.set_data_from_file(self.path)

		target.qt_parent.setWindowTitle(display_path(self.path))

		brws.notbusy()
		target.show()
		target.raise_()

	def plot3dNew(self, brws) :
		"""Make a new plot"""
		brws.busy()

		target = EMPlot3DWidget()
		brws.viewplot3d.append(target)
		target.set_data_from_file(self.path)
		#target.set_data(data, remove_directories_from_name(self.path, 1))

		target.qt_parent.setWindowTitle(display_path(self.path))

		brws.notbusy()
		target.show()
		target.raise_()


class EMFolderFileType(EMFileType) :
	"""FileType for Folders"""

	@staticmethod
	def name() :
		"""The unique name of this FileType. Stored in EMDirEntry.filetype for each file."""
		return "Folder"

	@staticmethod
	def isValid(path, header) :
		"""Returns (size, n, dim) if the referenced path is a file of this type, None if not valid. The first 4k block of data from the file is provided as well to avoid unnecessary file access."""
		return False

	@staticmethod
	def infoClass() :
		"""Returns a reference to the QWidget subclass for displaying information about this file"""
		return EMFolderInfoPane

	def actions(self) :
		"""Returns a list of (name, callback) tuples detailing the operations the user can call on the current file"""
		return []


class EMJSONFileType(EMFileType) :
	"""FileType for JSON files"""

	@staticmethod
	def name() :
		"""The unique name of this FileType. Stored in EMDirEntry.filetype for each file."""
		return "JSON"

	@staticmethod
	def isValid(path, header) :
		"""Returns (size, n, dim) if the referenced path is a file of this type, None if not valid. The first 4k block of data from the file is provided as well to avoid unnecessary file access."""

		header=header.decode("utf-8")
		if path[-5:] == ".json" and header.strip()[0] == "{" : 
			sz = len(js_open_dict(path).keys())
			return (humansize(os.stat(path).st_size), sz, "-")
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
		self.keys = list(self.js.keys())
		self.dim = (0, 0, 0)

	def __del__(self) :
		try : self.js.close()
		except : pass

	def actions(self) :
		"""Returns a list of (name, help, callback) tuples detailing the operations the user can call on the current file"""
		# No actions for now...
		if len(self.js.values())==0:
			return []
		v=self.js.values()[0]
		if type(v)==dict:
			if "xform.align3d" in v:
				ret=[
					("Plot 2D", "plot xform", self.plot2dApp),
					("Plot 2D+", "plot xform", self.plot2dNew),
					("Histogram", "histogram xform", self.histApp),
					("Histogram+", "histogram xform", self.histNew)
					]
				try:
					fsp,n=eval(self.js.keys()[0])
					tmp=EMData(fsp,n)
					ret.append(("All XYZ", "Show restricted XYZ projections of all", self.show2dStack3sec))
				except:
					pass
			
				return ret

		return []
	
	def plot2dNew(self, brws):
		self.plot2dApp(brws, True)
		
	def plot2dApp(self, brws, new=False) :
		"""Append self to current plot"""
		brws.busy()
		
		modifiers = QtWidgets.QApplication.keyboardModifiers()
		if modifiers == QtCore.Qt.ShiftModifier:
			inv=False
		else:
			inv=True

		rows = []
		for k in self.keys:
			dct = self.js[k]
			#### inverse since we want the transform from reference to particle
			tf = Transform(dct[u'xform.align3d'])
			if inv: 
				tf=tf.inverse()
			t = tf.get_trans()
			r = tf.get_rotation()
			row = [r["az"],r["alt"],r["phi"],t[0],t[1],t[2]]
			if "score" in dct:
				row.append(dct["score"])

			rows.append([float(x) for x in row])
		
		data=np.array(rows).T.tolist()
		
		if new:
			target = EMPlot2DWidget()
			brws.viewplot2d.append(target)
		else:
			try :
				target = brws.viewplot2d[-1]
				if target.closed : 
					brws.viewplot2d.pop()
					raise Exception
				#target.set_data(data, remove_directories_from_name(self.path, 1))
			except :
				target = EMPlot2DWidget()
				brws.viewplot2d.append(target)
		
		target.set_data(data, display_path(self.path))

		target.qt_parent.setWindowTitle(display_path(self.path))

		brws.notbusy()
		target.show()
		target.raise_()
	
	def histNew(self, brws):
		self.histApp(brws, True)
	
	def histApp(self, brws, new=False) :
		"""Make a new plot"""
		modifiers = QtWidgets.QApplication.keyboardModifiers()
		if modifiers == QtCore.Qt.ShiftModifier:
			inv=False
		else:
			inv=True

		rows = []
		for k in self.keys:
			dct = self.js[k]
			#### inverse since we want the transform from reference to particle
			tf = Transform(dct[u'xform.align3d'])
			if inv: 
				tf=tf.inverse()
			t = tf.get_trans()
			r = tf.get_rotation()
			row = [r["az"],r["alt"],r["phi"],t[0],t[1],t[2]]
			if "score" in dct:
				row.append(dct["score"])

			rows.append([float(x) for x in row])
		
		data=np.array(rows).T.tolist()
		if new==False and len(brws.viewhist)>0:
			target=brws.viewhist[-1]
		
		else:
			target = EMHistogramWidget()
			brws.viewhist.append(target)
			
		target.set_data(data, display_path(self.path))

		target.qt_parent.setWindowTitle(display_path(self.path))

		brws.notbusy()
		target.show()
		target.raise_()

	def show2dStack3sec(self, brws) :
		"""A set of 2-D images derived from a stack of 3-D Volumes referenced from a JSON file"""
		try:
			ret=self.secparm.exec_()
		except:
			self.secparm=EMSliceParamDialog(brws,len(self.js))
			ret=self.secparm.exec_()
		
		if not ret: return	# cancel
	
		# with JSON files, img0 and 1 are particle number in the referenced key
		img0=self.secparm.wspinmin.value()
		img0=max(img0,0)
		img1=self.secparm.wspinmax.value()
		if img1<0 or img1>len(self.js): img1=len(self.js)
		if img1<=img0 : img1+=1
		imgstep=self.secparm.wspinstep.value()
		layers=self.secparm.wspinlayers.value()
		center=self.secparm.wspincenter.value()
		lowpass=float(self.secparm.wlelp.text())
		highpass=float(self.secparm.wlehp.text())
		applyxf=self.secparm.wcheckxf.checkState()
		applysym=str(self.secparm.wlesym.text())
		stkout=self.secparm.wcheckstk.checkState()
		
		maskfsp=str(self.secparm.wlemask.text())
		mask=None
		if len(maskfsp)>4 :
			try: mask=EMData(maskfsp,0)
			except: print("ERROR: unable to read mask file: ",maskfsp)

		reffsp=str(self.secparm.wleref.text())
		ref=None
		if len(reffsp)>4 :
			try: ref=EMData(reffsp,0)
			except: print("ERROR: unable to read ref file: ",maskfsp)
		
#		print(img0,img1,ungstep,layers)
		
		brws.busy()
				
		progress = QtWidgets.QProgressDialog("Generating projections", "Cancel", 0, (img1-img0-1)/imgstep,None)
		
		data=[]
		if ref!=None:
			if mask!=None: ref.mult(mask)
			hall=makeOrthoProj(ref,layers,center,highpass,lowpass,stkout)
			if isinstance(hall,EMData): data.append(hall)
			else: data.extend(hall)
			
		i=0
		for k in self.keys:
			try:
				fsp,n=eval(k)
			except:
				print("Key error: ",k)
				continue
			
			if n<img0 or n>=img1 or (n-img0)%imgstep!=0 : 
				continue
		
			try: ptcl=EMData(fsp,n)
			except:
				print(f"Error reading {self.path} {i}")
				continue
			
			# We use the xform.align3d if it's in the JSON file, otherwise try the header
			try: xf=self.js[k]["xform.align3d"]
			except: 
				try: xf=ptcl["xform.align3d"]
				except: xf=Transform()
			if applyxf: ptcl.process_inplace("xform",{"transform":xf})
			if applysym is not None and applysym!="" and applysym.lower()!="c1": ptcl.process_inplace("xform.applysym",{"sym":applysym})
			if mask!=None : ptcl.mult(mask)
						
			time.sleep(0.001)
			get_application().processEvents()

			hall=makeOrthoProj(ptcl,layers,center,highpass,lowpass,stkout)
			
			if isinstance(hall,EMData):
				# copy some parameters from the original image header
				for kh in ["ptcl_repr","class_ptcl_idxs","class_ptlc_src","orig_file","orig_n","source_path","source_n"]:
					try: hall[kh]=ptcl[kh]
					except: pass
				data.append(hall)
			else:
				data.extend(hall)
			
			# copy some parameters from the JSON file to the final image for possible display
			hall["xform.align3d"]=xf
			try: hall["score"]=self.js[k]["score"]
			except: pass
			try: hall["coverage"]=self.js[k]["coverage"]
			except: pass
			
			progress.setValue(i-img0)
			i+=1
			
			if progress.wasCanceled():
#				progress.close()
				return

		if stkout:
			target = EMImage2DWidget()
			target.set_data(data,self.path)
			brws.view2d.append(target)
		else:
			target = EMImageMXWidget()
			target.set_data(data,self.path)
			target.mx_image_double.connect(target.mouse_double_click)		# this makes class average viewing work in app mode
			brws.view2ds.append(target)

		target.qt_parent.setWindowTitle("Stack - "+display_path(self.path))

		brws.notbusy()
		target.show()
		target.raise_()


class EMBdbFileType(EMFileType) :
	"""FileType for Folders"""

	@staticmethod
	def name() :
		"""The unique name of this FileType. Stored in EMDirEntry.filetype for each file."""
		return "BDB"

	@staticmethod
	def isValid(path, header) :
		"""Returns (size, n, dim) if the referenced path is a file of this type, None if not valid. The first 4k block of data from the file is provided as well to avoid unnecessary file access."""
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

		if self.nimg == 0 : self.keys = list(self.bdb.keys())
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
				("Show 2D", "Show all images, one at a time in current window", self.show2dSingle), ("Show 2D+", "Show all images, one at a time in a new window", self.show2dSingleNew), 
				("Avg All", "Unaligned average of entire stack",self.show2dAvg),("Avg Sample","Averages random min(1/4 of images,1000) multiple times",self.show2dAvgRnd),
				("FilterTool", "Open in e2filtertool.py", self.showFilterTool), ("Save As", "Saves images in new file format", self.saveAs)]
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

		target.qt_parent.setWindowTitle(display_path(self.path))

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
		else : print("Sorry, I don't know how to run Chimera on this platform")


class EMImageFileType(EMFileType) :
	"""FileType for files containing a single 2-D image"""
	module_closed = QtCore.pyqtSignal()

	def __init__(self, path) :
		if path[:2] == "./" : path = path[2:]

		EMFileType.__init__(self, path)	# the current path this FileType is representing
		self.nimg = EMUtil.get_image_count(path)
		im0 = EMData(path, 0, True)
		self.dim = (im0["nx"], im0["ny"], im0["nz"])

	def closeEvent(self, event) :
#		E2saveappwin("e2display", "main", self)

		self.updthreadexit = True

		for w in self.view2d+self.view2ds+self.view3d+self.viewplot2d+self.viewplot3d+self.viewhist :
			w.close()

		if self.infowin != None :
			self.infowin.close()

		event.accept()
		# self.app().close_specific(self)
		self.module_closed.emit()

	@staticmethod
	def name() :
		"""The unique name of this FileType. Stored in EMDirEntry.filetype for each file."""
		return "Image"

	@staticmethod
	def isValid(path, header) :
		"""Returns (size, n, dim) if the referenced path is a file of this type, None if not valid. The first 4k block of data from the file is provided as well to avoid unnecessary file access."""
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
				("FilterTool", "Open in e2filtertool.py", self.showFilterTool), ("Rng XYZ", "Show restricted XYZ projection", self.show2dStack3sec), ("Save As", "Saves images in new file format", self.saveAs)]
#				("FilterTool", "Open in e2filtertool.py", self.showFilterTool), ("ProjXYZ", "Make projections along Z,Y,X", self.showProjXYZ ), ("Save As", "Saves images in new file format", self.saveAs)]
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

		target.qt_parent.setWindowTitle(display_path(self.path))

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
		else : print("Sorry, I don't know how to run Chimera on this platform")


class EMStackFileType(EMFileType) :
	"""FileType for files containing a set of 1-3D images"""

	@staticmethod
	def name() :
		"""The unique name of this FileType. Stored in EMDirEntry.filetype for each file."""""
		return "Image Stack"

	@staticmethod
	def isValid(path, header) :
		"""Returns (size, n, dim) if the referenced path is a file of this type, None if not valid. The first 4k block of data from the file is provided as well to avoid unnecessary file access."""
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
			for i in range(1, 10) :
				try : im0 = EMData(path, i, True)
				except : continue
				break

		try : self.dim = (im0["nx"], im0["ny"], im0["nz"])
		except :
			print("First 10 images all missing in ", path)
			self.dim = "?"
		
		self.xfparms=False
		if path.endswith(".lst"):
			info=load_lst_params(path, [0])[0]
			if ("xform.align3d" in info) or ("xform.projection" in info):
				self.xfparms=True
		
	def actions(self) :
		"""Returns a list of (name, callback) tuples detailing the operations the user can call on the current file"""
		# 3-D stack
		if self.nimg > 1 and self.dim[2] > 1:
			rtr= [("Show all 3D", "Show all in a single 3D window", self.show3DAll), ("Show 1st 3D", "Show only the first volume", self.show3DNew),("Show 1st 2D", "Show first volume as 2D stack", self.show2dSingle30),("Show 2nd 2D", "Show second volume as 2D stack", self.show2dSingle31),("Show All Zproj", "Show Z projection of all volumes", self.show2dStack3z),("All XYZ", "Show restricted XYZ projections of all", self.show2dStack3sec), ("Chimera", "Open in chimera (if installed)", self.showChimera), ("Save As", "Saves images in new file format", self.saveAs)]
		# 2-D stack
		elif self.nimg > 1 and self.dim[1] > 1:
			rtr=[("Show Stack", "Show all images together in one window", self.show2dStack), ("Show Stack+", "Show all images together in a new window", self.show2dStackNew), 
				("Show 2D", "Show all images, one at a time in current window", self.show2dSingle), ("Show 2D+", "Show all images, one at a time in a new window", self.show2dSingleNew), 
				("Avg All", "Unaligned average of entire stack",self.show2dAvg),("Avg Rnd Subset","Averages random min(1/4 of images,1000) multiple times",self.show2dAvgRnd),
				("FilterTool", "Open in e2filtertool.py", self.showFilterTool), ("Save As", "Saves images in new file format", self.saveAs)]
			if self.dim[0]>=3 and self.dim[0]<=5: rtr.append(("Spheres","Each X line is X-Y-Z[-A[-S]]. Show as spheres in 3-D",self.showSpheres))
			else: print("Nope ",self.dim)
			
		# 1-D stack
		elif self.nimg > 1:
			rtr= [("Plot 2D", "Plot all on a single 2-D plot", self.plot2dNew), ("Save As", "Saves images in new file format", self.saveAs)]
		else : 
			rtr=[]
			print("Error: stackfile with < 2 images ? (%s)"%self.path)
			
		if self.xfparms:
			rtr.extend([("Plot 2D", "Plot xform", self.plot2dLstApp),("Plot 2D+", "plot xform in new window", self.plot2dLstNew)])

		return rtr

	def showSpheres(self,brws):
		"""New 3-D window"""
		brws.busy()

		target = emscene3d.EMScene3D()
		brws.view3d.append(target)

		data=EMData.read_images(self.path)
		gaussplots=[emshapeitem3d.EMScatterPlot3D() for i in data]
		for i,d in enumerate(data):
			target.insertNewNode(f"Iter_{i}",gaussplots[i])
			if i>0: gaussplots[i].setVisibleItem(False)
			cp=data[i].numpy().copy().transpose()
			cp[:3]*=256.0
			gaussplots[i].setData(cp)

		# target.initialViewportDims(data.getData().get_xsize())	# Scale viewport to object size
		# target.setCurrentSelection(iso)				# Set isosurface to display upon inspector loading

		brws.notbusy()
		target.setWindowTitle(display_path(self.path))

		target.show()
		target.raise_()


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
		else : print("Sorry, I don't know how to run Chimera on this platform")

	def plot2dNew(self, brws) :
		self.plot2dApp(brws, True)

	def plot2dApp(self, brws, new=False) :
		"""Append self to current plot"""
		brws.busy()

		data = EMData.read_images(self.path)

		if new:
			target = EMPlot2DWidget()
			brws.viewplot2d.append(target)
			
		else:
			try :
				target = brws.viewplot2d[-1]
			except :
				target = EMPlot2DWidget()
				brws.viewplot2d.append(target)
				
		for i,d in enumerate(data): target.set_data(d, f"{i},{display_path(self.path)}")
		
		target.qt_parent.setWindowTitle(display_path(self.path))

		brws.notbusy()
		target.show()
		target.raise_()

	def plot2dLstNew(self, brws):
		self.plot2dLstApp(brws, True)
		
	def plot2dLstApp(self, brws, new=False) :
		"""Append self to current plot"""
		brws.busy()
		rows = []
		print("Reading from {}...".format(self.path))
		params=load_lst_params(self.path)
		keys=list(params[0].keys())
		if "xform.align3d" in keys:
			xfs=[p["xform.align3d"].inverse().get_params("eman") for p in params]
		elif "xform.projection" in keys:
			xfs=[p["xform.projection"].get_params("eman") for p in params]
			
		xfkeys=["az", "alt", "phi", "tx", "ty", "tz"]
		
		data=[[x[k] for k in xfkeys] for x in xfs]
		data=np.array(data)
		
		for etc in ["score", "class"]:
			if etc in keys:
				xfkeys+=[etc]
				scr=np.array([p[etc] for p in params])
				data=np.hstack([data, scr[:,None]])
		
		data=data.T.tolist()
		print("The {} columns are {}".format(len(xfkeys), xfkeys))
		
		if new:
			target = EMPlot2DWidget()
			brws.viewplot2d.append(target)
		else:
			try :
				target = brws.viewplot2d[-1]
				if target.closed : 
					brws.viewplot2d.pop()
					raise Exception
				#target.set_data(data, remove_directories_from_name(self.path, 1))
			except :
				target = EMPlot2DWidget()
				brws.viewplot2d.append(target)
		
		target.set_data(data, display_path(self.path))

		target.qt_parent.setWindowTitle(display_path(self.path))

		brws.notbusy()
		target.show()
		target.raise_()


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
		The first 4k block of data from the file is provided as well to avoid unnecessary file access.
		"""
		proper_exts = ['pdb','ent']
		ext = os.path.basename(path).split('.')[-1]
		if ext not in proper_exts: return False
		
		try: s=header.decode("utf-8")
		except: return False
#		if not isprint(header) : return False			# demand printable Ascii. FIXME: what about unicode ?
		
		try : size = os.stat(path)[6]
		except : return False
		
		if size > 5000000: dim = "big"
		else :
			f = open(path, "r").read()
			lns = max(f.count("\n"), f.count("\r"))
			dim = "%d ln"%lns
		
		return (size, "-", dim)

	def actions(self):
		"""
		Returns a list of (name, callback) tuples detailing the operations the user can call on the current file
		"""
		return [("Show Ball and Stick", "Show ball and stick representation of this PDB model in a new 3D window", self.showBallStick3dNew), ("Show Ball and Stick +", "Show ball and stick representation of this PDB model in the current 3D window", self.showBallStick3dApp), ("Show Spheres", "Show spheres representation of this PDB model in a new 3D window", self.showSpheres3dNew), ("Show Spheres +", "Show spheres representation of this PDB model in the current 3D window", self.showSpheres3dApp), ("Chimera", "Open this PDB file in chimera (if installed)", self.showChimera), ("Save As", "Saves a copy of the selected PDB file", self.saveAs)]

	def showSpheres3dApp(self, brws):
		"""New 3-D window"""
		brws.busy()
		pdb_model = EMPDBItem3D(self.path)
		target = emscene3d.EMScene3D()
		brws.view3d.append(target)
		target.insertNewNode(display_path(self.path),pdb_model)
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
		target.insertNewNode(display_path(self.path),pdb_model, parentnode = target)
		modeltype = EMSphereModel(self.path)
		target.insertNewNode(modeltype.representation, modeltype, parentnode = pdb_model)
		target.initialViewportDims(pdb_model.getBoundingBoxDimensions()[0])	# Scale viewport to object size
		target.setCurrentSelection(modeltype)	# Set style to display upon inspector loading
		#target.updateSG()	# this is needed because this might just be an addition to the SG rather than initialization
		target.setWindowTitle(pdb_model.getName())
		brws.notbusy()
		target.show()
		target.raise_()

	def showBallStick3dApp(self, brws):
		"""New 3-D window"""
		brws.busy()
		pdb_model = EMPDBItem3D(self.path)
		target = emscene3d.EMScene3D()
		brws.view3d.append(target)
		target.insertNewNode(display_path(self.path),pdb_model)
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
		target.insertNewNode(display_path(self.path),pdb_model, parentnode = target)
		modeltype = EMBallStickModel(self.path)#parent=pdb_model)
		target.insertNewNode(modeltype.representation, modeltype, parentnode = pdb_model)
		target.initialViewportDims(pdb_model.getBoundingBoxDimensions()[0])	# Scale viewport to object size
		target.setCurrentSelection(modeltype)	# Set style to display upon inspector loading
		#target.updateSG()	# this is needed because this might just be an addition to the SG rather than initialization
		target.setWindowTitle(pdb_model.getName())
		brws.notbusy()
		target.show()
		target.raise_()

	def showChimera(self,brws):
		"""Open in Chimera"""
		if get_platform() == "Linux":
			os.system("chimera %s &" % self.path)
		elif get_platform() == "Darwin":
			os.system("/Applications/Chimera.app/Contents/MacOS/chimera %s &" % self.path)
		else:
			print("Sorry, I don't know how to run Chimera on this platform")


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
	"HTML"        : EMHTMLFileType,
	"PDF"         : EMPDFFileType
}

# Note that image types are not included here, and are handled via a separate mechanism
# note that order is important in the tuples. The most specific filetype should go first, and
# the most general last (they will be checked in order)
EMFileType.extbyft = {
	".json" : (EMJSONFileType, EMTextFileType),
	".pdb"  : (EMPDBFileType,  EMTextFileType),
	".ent"  : (EMPDBFileType,  EMTextFileType),
	".txt"  : (EMPlotFileType, EMTextFileType),
	".mdoc"  : (EMTextFileType,),
	".rawtlt": (EMTextFileType,),
	".pdf"  : (EMPDFFileType,),
	".htm"  : (EMHTMLFileType, EMTextFileType),
	".html" : (EMHTMLFileType, EMTextFileType)
}

# Default Choices when extension doesn't work
# We don't need to test for things like Images because they are fully tested outside this mechanism
EMFileType.alltocheck = (EMPlotFileType, EMPDBFileType, EMTextFileType)

BDBWARN=False
class EMDirEntry(object) :
	"""Represents a directory entry in the filesystem"""

	# list of lambda functions to extract column values for sorting
	col = (lambda x:int(x.index), lambda x:x.name, lambda x:x.filetype if x.filetype!=None else "", lambda x:size_sortable(x.size), lambda x:str(x.dim), lambda x:safe_int(x.nimg), lambda x:x.date)
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
		if self.filetype in EMFileType.typesbyft:
			filetype = EMFileType.typesbyft[self.filetype]
		else:
			filetype=None
		return filetype

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
			print("Request for child %d of children %s (%d)"%(n, self.__children, len(self.__children)))
			traceback.print_stack()
			raise Exception

	def nChildren(self) :
		"""Count of children"""
		self.fillChildNames()

#		print "EMDirEntry.nChildren(%s) = %d"%(self.filepath, len(self.__children))

		return len(self.__children)

	def fillChildNames(self) :
		"""Makes sure that __children contains at LEAST a list of names. This function needs to reimplemented to make derived browsers, 
		NOTE!!!! You must have nimg implemented in your reimplementation (I know this is a bad design....)"""
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
				global BDBWARN
				self.__children.remove("EMAN2DB")
				if not BDBWARN : print("WARNING: BDB (EMAN2DB/) no longer supported. EMAN2.91 or earlier required to view.")
				BDBWARN=True

				#if self.dirregex != None :
					#if isinstance(self.dirregex, str) :
						#t = ["bdb:"+i for i in db_list_dicts("bdb:"+self.filepath) if matches_pats(i, self.dirregex)]
					#else :
						#t = ["bdb:"+i for i in db_list_dicts("bdb:"+self.filepath) if self.dirregex.match(i) != None]

##					for i in db_list_dicts("bdb:"+self.filepath) : print i, self.dirregex.search (i)
				#else :
					#t = ["bdb:"+i for i in db_list_dicts("bdb:"+self.filepath)]

				#self.__children.extend(t)

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
		if self.path()[:4].lower()!="bdb:" and not (os.path.isfile(self.path()) or os.path.islink(self.path())) : 
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

		# we do this this way because there are so many possible image file extensions, and sometimes
		# people use a non-standard one (esp for MRC files)
		try : self.nimg = EMUtil.get_image_count(self.path())
		except : self.nimg = 0

		# we have an image file
		if self.nimg > 0 :
			try : 
				tmp = EMData(self.path(), 0, True)		# try to read an image header for the file
			except :
				#print("Error: first image missing ! : ",self.path())
				return 0
				#for i in range(1, 10) :
					#try : tmp = EMData(self.path(), i, True)
					#except : continue
					#break
				#if i == 9 :
					#print("Error: all of the first 10 images are missing ! : ",self.path())
					#return 0

			if tmp["ny"] == 1 : self.dim = str(tmp["nx"])
			elif tmp["nz"] == 1 : self.dim = "%d x %d"%(tmp["nx"], tmp["ny"])
			else : self.dim = "%d x %d x %d"%(tmp["nx"], tmp["ny"], tmp["nz"])

			if self.nimg == 1 : self.filetype = "Image"
			else : self.filetype = "Image Stack"

		# Ok, we need to try to figure out what kind of file this is
		else :
			head = open(self.path(), "rb").read(16384)		# Most FileTypes should be able to identify themselves using the first 4K block of a file
			ext=os.path.splitext(self.path())[1]
			if ext not in EMFileType.extbyft:
				guesses=EMFileType.alltocheck
			else: guesses = EMFileType.extbyft[ext]		# This will get us a list of possible FileTypes for this extension

	#			print "-------\n", guesses

			for guess in guesses :
				ret=guess.isValid(self.path(), head)
				if ret==False: continue
				size, n, dim = ret

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

		if cache != None :
			try : cache.setval(cachename, (time.time(), self.dim, self.filetype, self.nimg, self.size), True)
			except : pass
		return 1

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

	if val > 2**30 : return "%d g"%(val/(2**30))
	elif val > 2**20 : return "%d m"%(val/(2**20))
	elif val > 2**10 : return "%d k"%(val/(2**10))
	return str(val)

	#if val > 1000000000 : return "%d g"%(old_div(val,1000000000))
	#elif val > 1000000 : return "%d m"%(old_div(val,1000000))
	#elif val > 1000 : return "%d k"%(old_div(val,1000))
	#return str(val)

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
		return 7

	def rowCount(self, parent) :
		"""Returns the number of children for a given parent"""
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
			print("Error with index ", index.row(), index.column())
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
				print("Parent index error:", str(index.__dict__))
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


class myQItemSelection(QtCore.QItemSelectionModel) :
	"""For debugging"""

	def select(self, tl, br) :
		print(tl.indexes()[0].row(), tl.indexes()[0].column(), int(br))
		QtCore.QItemSelectionModel.select(self, tl, QtCore.QItemSelectionModel.SelectionFlags(QtCore.QItemSelectionModel.ClearAndSelect+QtCore.QItemSelectionModel.Rows))


class EMInfoPane(QtWidgets.QWidget) :
	"""Subclasses of this class will be used to display information about specific files. Each EMFileType class will return the
	pointer to the appropriate infoPane subclass for displaying information about the file it represents. The subclass instances
	are allocated by the infoWin class"""

	def __init__(self, parent = None) :
		"""Set our GUI up"""
		QtWidgets.QWidget.__init__(self, parent)

		# self.setTitle("e2display.py Information Pane")

		self.setWindowTitle("e2display.py Information Pane") # Jesus

		# Root class represents no target
		self.hbl = QtWidgets.QHBoxLayout(self)
		self.lbl = QtWidgets.QLabel("No Information Available")
		self.hbl.addWidget(self.lbl)

	def display(self, target) :
		"""display information for the target EMDirEntry with EMFileType ftype"""
		# self.setTitle("e2display.py Information Pane")

		self.target = target
		self.setWindowTitle("e2display.py Information Pane") # Jesus

		return

	def busy(self) :
		pass

	def notbusy(self) :
		pass


class EMTextInfoPane(EMInfoPane) :
	def __init__(self, parent = None) :
		QtWidgets.QWidget.__init__(self, parent)

		self.vbl = QtWidgets.QVBoxLayout(self)

		# text editing widget
		self.text = QtWidgets.QTextEdit()
		self.text.setAcceptRichText(False)
		self.text.setReadOnly(True)
		self.vbl.addWidget(self.text)

		# Find box
		self.wfind = StringBox(label = "Find:")
		self.vbl.addWidget(self.wfind)

		# Buttons
		self.hbl = QtWidgets.QHBoxLayout()

		self.wbutedit = QtWidgets.QPushButton("Edit")
		self.hbl.addWidget(self.wbutedit)

		self.wbutcancel = QtWidgets.QPushButton("Revert")
		self.wbutcancel.setEnabled(False)
		self.hbl.addWidget(self.wbutcancel)

		self.wbutok = QtWidgets.QPushButton("Save")
		self.wbutok.setEnabled(False)
		self.hbl.addWidget(self.wbutok)

		self.vbl.addLayout(self.hbl)

		self.wfind.valueChanged.connect(self.find)
		self.wbutedit.clicked.connect(self.buttonEdit)
		self.wbutcancel.clicked.connect(self.buttonCancel)
		self.wbutok.clicked.connect(self.buttonOk)

	def display(self, data) :
		"""display information for the target EMDirEntry"""
		self.target = data

		self.text.setPlainText(open(self.target.path(), "r").read())

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

	def buttonEdit(self, tog=None) :
		self.text.setReadOnly(False)
		self.wbutedit.setEnabled(False)
		self.wbutcancel.setEnabled(True)
		self.wbutok.setEnabled(True)

	def buttonCancel(self, tog=None) :
		self.display(self.target)

	def buttonOk(self, tog=None) :
		try : open(self.target.path(), "w").write(str(self.text.toPlainText()))
		except : QtWidgets.QMessageBox.warning(self, "Error !", "File write failed")


class EMHTMLInfoPane(EMInfoPane) :
	def __init__(self, parent = None) :
		QtWidgets.QWidget.__init__(self, parent)

		self.vbl = QtWidgets.QVBoxLayout(self)

		# text editing widget
		self.text = QtWidgets.QTextEdit()
		self.text.setAcceptRichText(True)
		self.text.setReadOnly(True)
		self.vbl.addWidget(self.text)

		# Find box
		self.wfind = StringBox(label = "Find:")
		self.vbl.addWidget(self.wfind)

		# Buttons
		self.hbl = QtWidgets.QHBoxLayout()

		self.wbutedit = QtWidgets.QPushButton("Edit")
		self.hbl.addWidget(self.wbutedit)

		self.wbutcancel = QtWidgets.QPushButton("Revert")
		self.wbutcancel.setEnabled(False)
		self.hbl.addWidget(self.wbutcancel)

		self.wbutok = QtWidgets.QPushButton("Save")
		self.wbutok.setEnabled(False)
		self.hbl.addWidget(self.wbutok)

		self.vbl.addLayout(self.hbl)

		self.wfind.valueChanged.connect(self.find)
		self.wbutedit.clicked.connect(self.buttonEdit)
		self.wbutcancel.clicked.connect(self.buttonCancel)
		self.wbutok.clicked.connect(self.buttonOk)

	def display(self, data) :
		"""display information for the target EMDirEntry"""
		self.target = data
		self.text.setHtml(open(self.target.path(), "r").read())
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

	def buttonEdit(self,tog=None) :
		self.text.setReadOnly(False)
		self.wbutedit.setEnabled(False)
		self.wbutcancel.setEnabled(True)
		self.wbutok.setEnabled(True)

	def buttonCancel(self,tog=None) :
		self.display(self.target)

	def buttonOk(self,tog=None) :
		try : open(self.target.path(), "w").write(str(self.text.toHtml()))
		except : QtWidgets.QMessageBox.warning(self, "Error !", "File write failed")


class EMPDBInfoPane(EMInfoPane) :
	def __init__(self, parent = None) :
		QtWidgets.QWidget.__init__(self, parent)
		self.vbl = QtWidgets.QVBoxLayout(self)
		# text editing widget
		self.text = QtWidgets.QTextEdit()
		self.text.setAcceptRichText(False)
		self.text.setReadOnly(True)
		self.vbl.addWidget(self.text)
		# Find box
		self.wfind = StringBox(label = "Find:")
		self.vbl.addWidget(self.wfind)
		# Buttons
		self.hbl = QtWidgets.QHBoxLayout()
		self.wbutedit = QtWidgets.QPushButton("Edit")
		self.hbl.addWidget(self.wbutedit)
		self.wbutcancel = QtWidgets.QPushButton("Revert")
		self.wbutcancel.setEnabled(False)
		self.hbl.addWidget(self.wbutcancel)
		self.wbutok = QtWidgets.QPushButton("Save")
		self.wbutok.setEnabled(False)
		self.hbl.addWidget(self.wbutok)
		self.vbl.addLayout(self.hbl)
		self.wfind.valueChanged.connect(self.find)
		self.wbutedit.clicked.connect(self.buttonEdit)
		self.wbutcancel.clicked.connect(self.buttonCancel)
		self.wbutok.clicked.connect(self.buttonOk)

	def display(self, data) :
		"""display information for the target EMDirEntry"""
		self.target = data
		self.text.setPlainText(open(self.target.path(), "r").read())
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

	def buttonEdit(self,tog=None) :
		self.text.setReadOnly(False)
		self.wbutedit.setEnabled(False)
		self.wbutcancel.setEnabled(True)
		self.wbutok.setEnabled(True)

	def buttonCancel(self,tog=None) :
		self.display(self.target)

	def buttonOk(self,tog=None) :
		try : open(self.target.path(), "w").write(str(self.text.toPlainText()))
		except : QtWidgets.QMessageBox.warning(self, "Error !", "File write failed")


class EMPlotInfoPane(EMInfoPane) :
	def __init__(self, parent = None) :
		QtWidgets.QWidget.__init__(self, parent)

		self.gbl = QtWidgets.QGridLayout(self)

		# List as alternate mechanism for selecting image number(s)
		self.plotdata = QtWidgets.QTableWidget()
		self.gbl.addWidget(self.plotdata, 0, 0)

	def display(self, target) :
		"""display information for the target EMDirEntry"""
		self.target = target
		self.plotdata.clear()

		# read the data into a list of lists
		numc = 0
		data = []

		for l in open(target.path(), "r") :
			if l[0] == "#" : continue

			vals = [float(i) for i in renumfind.findall(l)]

			if len(vals) == 0 : continue

			if numc == 0 : numc = len(vals)
			elif numc != len(vals) : break
			data.append(vals)

			if len(data) == 2500 : break			# if the table is too big, we just do a ...

		if len(data) == 2500 : self.plotdata.setRowCount(2501)
		else : self.plotdata.setRowCount(len(data))

		v=np.array(data)
		mean=np.mean(data, axis=0)
		std=np.std(data, axis=0)
		self.plotdata.setColumnCount(numc)
		self.plotdata.setVerticalHeaderLabels(["Mean", "Std"]+[str(i) for i in range(len(data))])
		self.plotdata.setHorizontalHeaderLabels([str(i) for i in range(numc)])
		for c in range(numc) :
			self.plotdata.setItem(0, c, QtWidgets.QTableWidgetItem("%1.4g"%mean[c]))
			self.plotdata.setItem(1, c, QtWidgets.QTableWidgetItem("%1.4g"%std[c]))
		if len(data)==2500: data=data[:-2]

		for r in range(len(data)) :
			for c in range(numc) :
				self.plotdata.setItem(r+2, c, QtWidgets.QTableWidgetItem("%1.4g"%data[r][c]))

		if len(data) == 2500 :
			self.plotdata.setVerticalHeaderItem(2500, QtWidgets.QTableWidgetItem("..."))


class EMFolderInfoPane(EMInfoPane) :
	def __init__(self, parent = None) :
		QtWidgets.QWidget.__init__(self, parent)

		self.vbl = QtWidgets.QVBoxLayout(self)

	def display(self, target) :
		"""display information for the target EMDirEntry"""
		self.target = target


class EMBDBInfoPane(EMInfoPane) :
	maxim = 1000

	def __init__(self, parent = None) :
		QtWidgets.QWidget.__init__(self, parent)

		self.gbl = QtWidgets.QGridLayout(self)

		# Spinbox for selecting image number
		self.wimnum = QtWidgets.QSpinBox()
		self.wimnum.setRange(0, 0)
		self.gbl.addWidget(self.wimnum, 0, 0)

		# List as alternate mechanism for selecting image number(s)
		self.wimlist = QtWidgets.QListWidget()
		self.gbl.addWidget(self.wimlist, 1, 0)

		# Actual header contents
		self.wheadtree = QtWidgets.QTreeWidget()
		self.wheadtree.setColumnCount(2)
		self.wheadtree.setHeaderLabels(["Item", "Value"])
		self.gbl.addWidget(self.wheadtree, 0, 1, 2, 1)

		self.gbl.setColumnStretch(0, 1)
		self.gbl.setColumnStretch(1, 4)

		# Lower region has buttons for actions
		self.hbl2 = QtWidgets.QGridLayout()

		self.wbutmisc = []

		# 10 buttons for context-dependent actions
		self.hbl2.setRowStretch(0, 1)
		self.hbl2.setRowStretch(1, 1)

		for i in range(5) :
			self.hbl2.setColumnStretch(i, 2)
	
			for j in range(2) :
				self.wbutmisc.append(QtWidgets.QPushButton("-"))
				self.hbl2.addWidget(self.wbutmisc[-1], j, i)
				self.wbutmisc[-1].setEnabled(False)
				self.wbutmisc[-1].clicked.connect(lambda x, v = i*2+j :self.buttonMisc(v))

		# These just clean up the layout a bit
		self.wbutxx = QtWidgets.QLabel("")
		self.wbutxx.setMaximumHeight(12)
		self.hbl2.addWidget(self.wbutxx, 0, 6)
		self.wbutyy = QtWidgets.QLabel("")
		self.wbutyy.setMaximumHeight(12)
		self.hbl2.addWidget(self.wbutyy, 1, 6)

		self.gbl.addLayout(self.hbl2, 2, 0, 1, 2)

		self.wimnum.valueChanged[int].connect(self.imNumChange)
		self.wimlist.itemSelectionChanged.connect(self.imSelChange)
##		QtCore.QObject.connect(self.wbutedit, QtCore.SIGNAL('clicked(bool)'), self.buttonEdit)

		self.view2d = []
		self.view3d = []
		self.view2ds = []
		self.viewplot2d = []
		self.viewhist = []

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
			QtWidgets.QMessageBox.warning(self, "Error", "Sorry, cannot display string-keyed images")
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
			k = list(self.bdb.keys())
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
					b.setEnabled(True)
				except :
					b.setText("-")
					b.setToolTip("")
					b.setEnabled(False)
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
				itms.append(QtWidgets.QTreeWidgetItem(list((str(k), str(trg[k])))))
				if isinstance(trg[k], list) or isinstance(trg[k], tuple) or isinstance(trg[k], set) or isinstance(trg[k], dict) :
					self.addTreeItem(trg[k], itms[-1])
		elif isinstance(trg, list) or isinstance(trg, tuple) or isinstance(trg, set) :
			for k in trg :
				if isinstance(k, list) or isinstance(k, tuple) or isinstance(k, set) or isinstance(k, dict) :
					try : itms.append(QtWidgets.QTreeWidgetItem(list((k.__class__.__name__, ""))))
					except : itms.append(QtWidgets.QTreeWidgetItem(list(("??", ""))))
					self.addTreeItem(k, itms[-1])
				else :
					itms.append(QtWidgets.QTreeWidgetItem(list((str(k), ""))))
		else :
			itms.append(QtWidgets.QTreeWidgetItem(list((str(trg), ""))))

		if parent == None :
			self.wheadtree.addTopLevelItems(itms)
			self.wheadtree.resizeColumnToContents(0)
		else : parent.addChildren(itms)


class EMJSONInfoPane(EMInfoPane) :
	def __init__(self, parent = None) :
		QtWidgets.QWidget.__init__(self, parent)

		self.gbl = QtWidgets.QGridLayout(self)

		# List of keys
		self.wkeylist = QtWidgets.QListWidget()
		self.gbl.addWidget(self.wkeylist, 1, 0)

		# contents of a single key
		self.wheadtree = QtWidgets.QTreeWidget()
		self.wheadtree.setColumnCount(2)
		self.wheadtree.setHeaderLabels(["Key/#", "Value"])
		self.gbl.addWidget(self.wheadtree, 0, 1, 2, 1)

		self.gbl.setColumnStretch(0, 1)
		self.gbl.setColumnStretch(1, 4)

		# Lower region has buttons for actions
		self.hbl2 = QtWidgets.QGridLayout()

		self.wbutmisc = []

		# 10 buttons for context-dependent actions
		self.hbl2.setRowStretch(0, 1)
		self.hbl2.setRowStretch(1, 1)

		for i in range(5) :
			self.hbl2.setColumnStretch(i, 2)

			for j in range(2) :
				self.wbutmisc.append(QtWidgets.QPushButton("-"))
				self.hbl2.addWidget(self.wbutmisc[-1], j, i)
				self.wbutmisc[-1].setEnabled(False)
				self.wbutmisc[-1].clicked.connect(lambda x, v = i*2+j :self.buttonMisc(v))

		# These just clean up the layout a bit
		self.wbutxx = QtWidgets.QLabel("")
		self.wbutxx.setMaximumHeight(12)
		self.hbl2.addWidget(self.wbutxx, 0, 6)
		self.wbutyy = QtWidgets.QLabel("")
		self.wbutyy.setMaximumHeight(12)
		self.hbl2.addWidget(self.wbutyy, 1, 6)

		self.gbl.addLayout(self.hbl2, 2, 0, 1, 2)

		self.wkeylist.itemSelectionChanged.connect(self.imSelChange)
		self.wheadtree.itemExpanded[QtWidgets.QTreeWidgetItem].connect(self.treeExp)
		self.wheadtree.itemCollapsed[QtWidgets.QTreeWidgetItem].connect(self.treeExp)
		self.wheadtree.itemSelectionChanged.connect(self.treeSel)
		self.wheadtree.itemActivated[QtWidgets.QTreeWidgetItem, int].connect(self.treeAct)
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
		k = list(self.js.keys())
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
					itms.append(QtWidgets.QTreeWidgetItem(list((str(k), ""))))
					self.addTreeItem(trg[k], itms[-1])
				else : itms.append(QtWidgets.QTreeWidgetItem(list((str(k), str(trg[k])))))
		elif isinstance(trg, (list, tuple, set)) :
			if isinstance(trg, set) : trg = sorted(trg)		# make a list temporarily
			if len(trg) > 120 : vals = list(range(0, 50))+[-1]+list(range(len(trg)-50, len(trg)))
			else : vals = list(range(len(trg)))
			for k in vals :
				if k == -1 : itms.append(QtWidgets.QTreeWidgetItem(list(("...", "..."))))
				else :
					v = trg[k]
					if isinstance(v, (list, tuple, set, dict, EMAN2Ctf)) :
						itms.append(QtWidgets.QTreeWidgetItem(list((str(k), ""))))
						self.addTreeItem(v, itms[-1])
					else : itms.append(QtWidgets.QTreeWidgetItem(list((str(k), str(v)))))
		elif isinstance(trg, EMAN2Ctf) :
			itms.append(QtWidgets.QTreeWidgetItem(list(("EMAN2Ctf", ""))))
			subitms = []
			for k, v in list(trg.to_dict().items()) :
				if isinstance(v, (list, tuple)) :
					v = ["%1.3g"%i for i in v]
					subitms.append(QtWidgets.QTreeWidgetItem(list((str(k), ", ".join(v)))))
				else : subitms.append(QtWidgets.QTreeWidgetItem(list((str(k), str(v)))))
			itms[-1].addChildren(subitms)
		elif isinstance(trg, EMData) :
			itms.append(QtWidgets.QTreeWidgetItem(list(("EMData", ""))))
		else :
			itms.append(QtWidgets.QTreeWidgetItem(list((str(trg), ""))))

		if parent == None :
			self.wheadtree.addTopLevelItems(itms)
			self.wheadtree.resizeColumnToContents(0)
		else : parent.addChildren(itms)


class EMImageInfoPane(EMInfoPane) :
	maxim = 500

	def __init__(self, parent = None) :
		QtWidgets.QWidget.__init__(self, parent)

		self.gbl = QtWidgets.QGridLayout(self)

		# Actual header contents
		self.wheadtree = QtWidgets.QTreeWidget()
		self.wheadtree.setColumnCount(2)
		self.wheadtree.setHeaderLabels(["Item", "Value"])
		self.gbl.addWidget(self.wheadtree, 0, 0)

	def display(self, target) :
		"""display information for the target EMDirEntry"""
		self.target = target

		self.wheadtree.clear()
		try : trg = EMData(self.target.path(), 0, True).get_attr_dict()		# read the header only, discard the emdata object
		except :
			print("Error reading:", self.target.path(), key)

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
				itms.append(QtWidgets.QTreeWidgetItem(list((str(k), str(trg[k])))))
				if isinstance(trg[k], list) or isinstance(trg[k], tuple) or isinstance(trg[k], set) or isinstance(trg[k], dict) :
					self.addTreeItem(trg[k], itms[-1])
		elif isinstance(trg, list) or isinstance(trg, tuple) or isinstance(trg, set) :
			for k in trg :
				if isinstance(k, list) or isinstance(k, tuple) or isinstance(k, set) or isinstance(k, dict) :
					try : itms.append(QtWidgets.QTreeWidgetItem(list((k.__class__.__name__, ""))))
					except : itms.append(QtWidgets.QTreeWidgetItem(list(("??", ""))))
					self.addTreeItem(k, itms[-1])
				else :
					itms.append(QtWidgets.QTreeWidgetItem(list((str(k), ""))))
		else :
			itms.append(QtWidgets.QTreeWidgetItem(list((str(trg), ""))))

		if parent == None :
			self.wheadtree.addTopLevelItems(itms)
			self.wheadtree.resizeColumnToContents(0)
		else : parent.addChildren(itms)


class EMStackInfoPane(EMInfoPane) :
	module_closed = QtCore.pyqtSignal()
	maxim = 500

	def __init__(self, parent = None) :
		QtWidgets.QWidget.__init__(self, parent)

		# self.setWindowTitle("e2display.py Information Pane") # Jesus
		# self.setTitle("e2display.py Information Pane")

		self.gbl = QtWidgets.QGridLayout(self)

		self.gbl.label1 = QtWidgets.QLabel("Images") # Jesus
		self.gbl.addWidget(self.gbl.label1, 0, 0) # Jesus

		self.gbl.label2 = QtWidgets.QLabel("Header Info") # Jesus
		self.gbl.addWidget(self.gbl.label2, 0, 1) # Jesus

		'''Spinbox for selecting image number'''
		self.wimnum = QtWidgets.QSpinBox()

		# self.wimnum.setRange(0, 0) # JOHN
		# self.gbl.addWidget(self.wimnum, 0, 0) # JOHN

		self.wimnum.setRange(1, 0) # Jesus
		self.gbl.addWidget(self.wimnum, 1, 0) # Jesus

		'''List as alternate mechanism for selecting image number(s)'''
		self.wimlist = QtWidgets.QListWidget()

		# self.gbl.addWidget(self.wimlist, 1, 0) # JOHN

		self.gbl.addWidget(self.wimlist, 2, 0) # Jesus

		'''Actual header contents'''
		self.wheadtree = QtWidgets.QTreeWidget()

		# self.wheadtree.setColumnCount(2) #
		self.wheadtree.setHeaderLabels(["Item", "Value"])
		# self.gbl.addWidget(self.wheadtree, 0, 1, 2, 1) # JOHN
		self.gbl.addWidget(self.wheadtree, 1, 1, 2, 1) # Jesus

		self.gbl.setColumnStretch(0, 1)
		self.gbl.setColumnStretch(1, 4)

		'''Lower region has buttons for actions'''
		self.hbl2 = QtWidgets.QGridLayout()

		self.wbutmisc = []

		'''10 buttons for context-dependent actions'''
		# self.hbl2.setRowStretch(0, 1) # JOHN
		# self.hbl2.setRowStretch(1, 1) # JOHN

		self.hbl2.setRowStretch(1, 1) # Jesus
		self.hbl2.setRowStretch(2, 1) # Jesus

		for i in range(5) :
			self.hbl2.setColumnStretch(i, 2)
	
			for j in range(2) :
				self.wbutmisc.append(QtWidgets.QPushButton("-"))
				self.hbl2.addWidget(self.wbutmisc[-1], j, i)
				self.wbutmisc[-1].setEnabled(False)
				self.wbutmisc[-1].clicked.connect(lambda x, v = i*2+j :self.buttonMisc(v))

		# These just clean up the layout a bit
		self.wbutxx = QtWidgets.QLabel("")
		self.wbutxx.setMaximumHeight(12)
		# self.hbl2.addWidget(self.wbutxx, 0, 6) # JOHN
		self.hbl2.addWidget(self.wbutxx, 1, 6) # Jesus

		self.wbutyy = QtWidgets.QLabel("")
		self.wbutyy.setMaximumHeight(12)

		# self.hbl2.addWidget(self.wbutyy, 1, 6) # JOHN
		self.hbl2.addWidget(self.wbutyy, 2, 6) # Jesus

		# self.gbl.addLayout(self.hbl2, 2, 0, 1, 2) # JOHN
		self.gbl.addLayout(self.hbl2, 3, 0, 1, 2) # Jesus

		self.wimnum.valueChanged[int].connect(self.imNumChange)
		self.wimlist.itemSelectionChanged.connect(self.imSelChange)
#		QtCore.QObject.connect(self.wbutedit, QtCore.SIGNAL('clicked(bool)'), self.buttonEdit)
		self.view2d = []
		self.view3d = []
		self.view2ds = []
		self.viewplot2d = []
		self.viewhist = []

	def closeEvent(self, event) :
#		E2saveappwin("e2display", "main", self)

		self.updthreadexit = True

		for w in self.view2d+self.view2ds+self.view3d+self.viewplot2d+self.viewhist:
			w.close()

		if self.infowin != None :
			self.infowin.close()

		event.accept()

		# self.app().close_specific(self)
		self.module_closed.emit()

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
				b.setEnabled(True)
			except :
				b.setText("-")
				b.setToolTip("")
				b.setEnabled(False)

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
				print("Error with key:", val)
				return

		self.showItem(val)

	def showItem(self, key) :
		"""Shows header information for the selected item"""
		self.wheadtree.clear()

		try : trg = EMData(self.target.path(), key, True).get_attr_dict()		# read the header only, discard the emdata object
		except :
			print("Error reading:", self.target.path(), key)

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
				itms.append(QtWidgets.QTreeWidgetItem(list((str(k), str(trg[k])))))
				if isinstance(trg[k], list) or isinstance(trg[k], tuple) or isinstance(trg[k], set) or isinstance(trg[k], dict) :
					self.addTreeItem(trg[k], itms[-1])
		elif isinstance(trg, list) or isinstance(trg, tuple) or isinstance(trg, set) :
			for k in trg :
				if isinstance(k, list) or isinstance(k, tuple) or isinstance(k, set) or isinstance(k, dict) :
					try : itms.append(QtWidgets.QTreeWidgetItem(list((k.__class__.__name__, ""))))
					except : itms.append(QtWidgets.QTreeWidgetItem(list(("??", ""))))
					self.addTreeItem(k, itms[-1])
				else :
					itms.append(QtWidgets.QTreeWidgetItem(list((str(k), ""))))
		else :
			itms.append(QtWidgets.QTreeWidgetItem(list((str(trg), ""))))

		if parent == None :
			self.wheadtree.addTopLevelItems(itms)
			self.wheadtree.resizeColumnToContents(0)
		else : parent.addChildren(itms)


class EMInfoWin(QtWidgets.QWidget) :
	"""The info window"""
	winclosed = QtCore.pyqtSignal()

	def __init__(self, parent = None) :
		QtWidgets.QWidget.__init__(self, parent)

		self.target = None
		self.stack = QtWidgets.QStackedLayout(self)

		# We add one instance of 'infoPane' parent class to represent nothing
		self.stack.addWidget(EMInfoPane())

	def set_target(self, target, ftype) :
		"""Display the info pane for target EMDirEntry with EMFileType instance ftype"""
		self.target = target
		self.ftype = ftype

		if target == None :
			self.stack.setCurrentIndex(0)
			return

		if hasattr(ftype, "infoClass"):
			cls = ftype.infoClass()
		else:
			return

		for i in range(self.stack.count()) :
			if isinstance(self.stack.itemAt(i), cls) :
				self.stack.setCurrentIndex(i)
				pane = self.stack.itemAt(i)
				pane.display(target)
				break
		else :
			# If we got here, then we need to make a new instance of the appropriate pane
			if cls == None : print("No class ! (%s)"%str(ftype))
			#self.winclosed = QtCore.pyqtSignal()
			pane = cls()
			i = self.stack.addWidget(pane)		# add the new pane and get its index
			pane.display(target)
			self.stack.setCurrentIndex(i)		# put the new pane on top

	def closeEvent(self, event) :
		QtWidgets.QWidget.closeEvent(self, event)
		self.winclosed.emit()


class SortSelTree(QtWidgets.QTreeView) :
	"""This is a subclass of QtWidgets.QTreeView. It is almost identical but implements selection processing with sorting.
	The correct way of doing this in QT4.2 is to use a QSortFilterProxy object, but that won't work properly in this case."""

	def __init__(self, parent = None) :
		QtWidgets.QTreeView.__init__(self, parent)
		self.header().setSectionsClickable(True)
		self.header().sectionClicked[int].connect(self.colclick)
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
		QtWidgets.QTreeView.sortByColumn(self, col, ascend)

		# then set a new selection list
		sel = self.model().findSelected()
		if len(sel) == 0 :return

		qis = QtCore.QItemSelection()
		for i in sel : qis.select(i, i)
		self.selectionModel().select(qis, QtCore.QItemSelectionModel.ClearAndSelect|QtCore.QItemSelectionModel.Rows)

#		for i in sel : self.selectionModel().select(i, QtCore.QItemSelectionModel.ClearAndSelect)
#		self.update()

class EMSliceParamDialog(QtWidgets.QDialog):
	"""This modal dialog asks the user for parameters required for the XYZ 3-D stack viewer"""
	dlay=-1
	dcnt=0
	dlp="-1"
	dhp="-1"
	dmask=""
	sym="c1"
	oldcheck=0
	
	
	def __init__(self, parent = None,nimg = 1) :
		QtWidgets.QDialog.__init__(self,parent)
		nimg=int(nimg)
		
		self.setWindowTitle("Slice Parameters")

#		self.resize(780, 580)
		self.vbl = QtWidgets.QVBoxLayout(self)
		self.fol = QtWidgets.QFormLayout()
		self.vbl.addLayout(self.fol)

		if nimg>1:
			self.wspinmin=QtWidgets.QSpinBox()
			self.wspinmin.setRange(-1,nimg-1)
			self.wspinmin.setValue(-1)
			self.wspinmin.setToolTip("Start of range of volumes to display, -1 for all")
			self.fol.addRow("First Image #:",self.wspinmin)
			
			self.wspinmax=QtWidgets.QSpinBox()
			self.wspinmax.setRange(-1,nimg-1)
			self.wspinmax.setValue(-1)
			self.wspinmax.setToolTip("End (exclusive) of range of volumes to display, -1 for all")
			self.fol.addRow("Last Image #:",self.wspinmax)

			self.wspinstep=QtWidgets.QSpinBox()
			self.wspinstep.setRange(1,nimg//2)
			self.wspinstep.setValue(1)
			self.wspinstep.setToolTip("Step for range of volumes to display (partial display of file)")
			self.fol.addRow("Step:",self.wspinstep)

		self.wspinlayers=QtWidgets.QSpinBox()
		self.wspinlayers.setRange(-1,256)
		self.wspinlayers.setValue(self.dlay)
		self.wspinlayers.setToolTip("Projection about center +- selected number of layers, eg- 0 -> central section only, -1 full projection")
		self.fol.addRow("Sum Layers (0->1, 1->3, 2->5, ...):",self.wspinlayers)

		self.wspincenter=QtWidgets.QSpinBox()
		self.wspincenter.setRange(-256,256)
		self.wspincenter.setValue(self.dcnt)
		self.wspincenter.setToolTip("Center about which sum is generated, 0=center of volume")
		self.fol.addRow("Center for sum on Z:",self.wspincenter)

		self.wlelp=QtWidgets.QLineEdit(self.dlp)
		self.wlelp.setToolTip("If >0 applies a low-pass filter in 2D. Specify in A.")
		self.fol.addRow("Lowpass (A):",self.wlelp)

		self.wlehp=QtWidgets.QLineEdit(self.dhp)
		self.wlehp.setToolTip("If >0 applies a high-pass filter in 2D. Specify in A.")
		self.fol.addRow("Highpass (A):",self.wlehp)

		self.wlemask=QtWidgets.QLineEdit(self.dmask)
		self.wlemask.setToolTip("Optional filename of a mask volume (same dimensions)")
		self.fol.addRow("Mask volume:",self.wlemask)

		self.wleref=QtWidgets.QLineEdit("")
		self.wleref.setToolTip("Optional filename of a reference volume (same dimensions)")
		self.fol.addRow("Reference volume:",self.wleref)

		self.wcheckxf=QtWidgets.QCheckBox("enable")
		self.wcheckxf.setChecked(0)
		self.wcheckxf.setToolTip("If set, applies the xform from the JSON/lst file or image header before making projections")
		self.fol.addRow("Apply xform.align3d from image:",self.wcheckxf)

		self.wlesym=QtWidgets.QLineEdit(self.sym)
		self.wlesym.setToolTip("Applies the specified symmetry to each particle after xform but before projection")
		self.fol.addRow("Symmetry after xform:",self.wlesym)

		self.wcheckstk=QtWidgets.QCheckBox("enable")
		self.wcheckstk.setChecked(0)
		self.wcheckstk.setToolTip("If set, makes 3 square images instead of a single rectangular image. Good for FFTs.")
		self.fol.addRow("Stack output:",self.wcheckstk)

		self.wcheckoldwin=QtWidgets.QCheckBox("enable")
		self.wcheckoldwin.setChecked(self.oldcheck)
		self.wcheckoldwin.setToolTip("If set, adds the new projection set to the last existing window")
		self.fol.addRow("Same window:",self.wcheckoldwin)


		self.bhb = QtWidgets.QHBoxLayout()
		self.vbl.addLayout(self.bhb)
		self.wbutok = QtWidgets.QPushButton("OK")
		self.bhb.addWidget(self.wbutok)

		self.wbutcancel = QtWidgets.QPushButton("Cancel")
		self.bhb.addWidget(self.wbutcancel)
	
		self.wbutok.clicked.connect(self.okpress)
		self.wbutcancel.clicked.connect(self.reject)
		self.wbutok.setDefault(1)
		
	def okpress(self,state):
		EMSliceParamDialog.dlay=self.wspinlayers.value()
		EMSliceParamDialog.dcnt=self.wspincenter.value()
		EMSliceParamDialog.dlp=self.wlelp.text()
		EMSliceParamDialog.dhp=self.wlehp.text()
		EMSliceParamDialog.dmask=self.wlemask.text()
		EMSliceParamDialog.sym=self.wlesym.text()
		EMSliceParamDialog.oldcheck=self.wcheckoldwin.checkState()
		self.accept()

class EMBrowserWidget(QtWidgets.QWidget) :
	"""This widget is a file browser for EMAN2. In addition to being a regular file browser, it supports:
	- getting information about recognized data types
	- embedding BDB: databases into the observed filesystem
	- remote database access (EMEN2)*
	"""
	ok = QtCore.pyqtSignal()
	cancel = QtCore.pyqtSignal()
	module_closed = QtCore.pyqtSignal()

	def __init__(self, parent = None, withmodal = False, multiselect = False, startpath = ".", setsmode = None, dirregex="") :
		"""withmodal - if specified will have ok/cancel buttons, and provide a mechanism for a return value (not truly modal)
		multiselect - if True, multiple files can be simultaneously selected
		startpath - default "."
		setsmode - Used during bad particle marking
		dirregex - default "", a regular expression for filtering filenames (directory names not filtered)
		"""
		# although this looks dumb it is necessary to break Python's issue with circular imports(a major weakness of Python IMO)
		global emscene3d, emdataitem3d, emshapeitem3d
		from . import emscene3d
		from . import emdataitem3d
		from . import emshapeitem3d

		QtWidgets.QWidget.__init__(self, parent)
		
		cwd=os.path.basename(os.getcwd())
		self.setWindowTitle(f"{cwd} - e2display") # Jesus

		# label = QtWidgets.QLabel(self);
      # label.setText("Window Title");
      # self.setWindowTitle("Window Title");

		self.withmodal = withmodal
		self.multiselect = multiselect

		self.resize(780, 580)
		self.gbl = QtWidgets.QGridLayout(self)

		# Top Toolbar area
		self.wtoolhbl = QtWidgets.QHBoxLayout()
		self.wtoolhbl.setContentsMargins(0, 0, 0, 0)

		self.wbutback = QtWidgets.QPushButton(chr(0x2190))
		self.wbutback.setMaximumWidth(36)
		self.wbutback.setEnabled(False)
		self.wtoolhbl.addWidget(self.wbutback, 0)

		self.wbutfwd = QtWidgets.QPushButton(chr(0x2192))
		self.wbutfwd.setMaximumWidth(36)
		self.wbutfwd.setEnabled(False)
		self.wtoolhbl.addWidget(self.wbutfwd, 0)

		# Text line for showing (or editing) full path
		self.lpath = QtWidgets.QLabel("  Path:")
		self.wtoolhbl.addWidget(self.lpath)

		self.wpath = QtWidgets.QLineEdit()
		self.wtoolhbl.addWidget(self.wpath, 5)

		# self.wspacet1 = QtWidgets.QSpacerItem(100, 10, QtWidgets.QSizePolicy.MinimumExpanding)
		# self.wtoolhbl.addSpacerItem(self.wspacet1)

		self.wbutinfo = QtWidgets.QPushButton("Info")
		self.wbutinfo.setCheckable(True)
		self.wtoolhbl.addWidget(self.wbutinfo, 1)

		self.gbl.addLayout(self.wtoolhbl, 0, 0, 1, 2)

		# 2nd Top Toolbar area
		self.wtoolhbl2 = QtWidgets.QHBoxLayout()
		self.wtoolhbl2.setContentsMargins(0, 0, 0, 0)

		self.wbutup = QtWidgets.QPushButton(chr(0x2191))
		self.wbutup.setMaximumWidth(36)
		self.wtoolhbl2.addWidget(self.wbutup, 0)

		self.wbutrefresh = QtWidgets.QPushButton(chr(0x21ba))
		self.wbutrefresh.setMaximumWidth(36)
		self.wtoolhbl2.addWidget(self.wbutrefresh, 0)

		# Text line for showing (or editing) full path
		self.lfilter = QtWidgets.QLabel("Filter:")
		self.wtoolhbl2.addWidget(self.lfilter)

		self.wfilter = QtWidgets.QComboBox()
		self.wfilter.setEditable(True)
		self.wfilter.setInsertPolicy(QtWidgets.QComboBox.InsertAtBottom)
		self.wfilter.addItem("")
		self.wfilter.addItem(r"(.(?!_ctf))*$")
		self.wfilter.addItem(r".*\.img")
		self.wfilter.addItem(r".*\.box")
		self.wfilter.addItem(r".*\.hdf")
		self.wfilter.addItem(r".*_ptcls$")
		self.wfilter.addItem(r".*\.mrc")
		self.wfilter.addItem(r".*\.tif")
		self.wfilter.addItem(r".*\.pdb")
		self.wfilter.addItem("help")
		self.wtoolhbl2.addWidget(self.wfilter, 5)
		if dirregex!="":
			self.wfilter.setEditText(dirregex)
		# self.wspacet1 = QtWidgets.QSpacerItem(100, 10, QtWidgets.QSizePolicy.MinimumExpanding)
		# self.wtoolhbl.addSpacerItem(self.wspacet1)

		self.selectall = QtWidgets.QPushButton("Sel All")
		self.wtoolhbl2.addWidget(self.selectall, 1)
		self.selectall.setEnabled(withmodal)

		self.gbl.addLayout(self.wtoolhbl2, 1, 0, 1, 2)

		### Central vertical region has bookmarks and tree
		# Bookmarks implemented with a toolbar in a frame

		self.wbookmarkfr = QtWidgets.QFrame()
		self.wbookmarkfr.setFrameStyle(QtWidgets.QFrame.StyledPanel|QtWidgets.QFrame.Raised)
		self.wbmfrbl = QtWidgets.QVBoxLayout(self.wbookmarkfr)

		self.wbookmarks = QtWidgets.QToolBar()
		# self.wbookmarks.setAutoFillBackground(True)
		# self.wbookmarks.setBackgroundRole(QtGui.QPalette.Dark)
		self.wbookmarks.setOrientation(2)
		#self.addBookmark("EMEN2", "emen2:")
		#self.wbookmarks.addSeparator()
		#self.addBookmark("SSH", "ssh:")
		#self.wbookmarks.addSeparator()
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
		self.hbl2 = QtWidgets.QGridLayout()

		self.wbutmisc = []

		# 10 buttons for context-dependent actions
		self.hbl2.setRowStretch(0, 1)
		self.hbl2.setRowStretch(1, 1)

		for i in range(5) :
			self.hbl2.setColumnStretch(i, 2)
	
			for j in range(2) :
				self.wbutmisc.append(QtWidgets.QPushButton("-"))
				self.hbl2.addWidget(self.wbutmisc[-1], j, i)
				self.wbutmisc[-1].setEnabled(False)
#				self.wbutmisc[-1].setEnabled(False)
				self.wbutmisc[-1].clicked.connect(lambda x, v = i*2+j :self.buttonMisc(v))

		self.wbutxx = QtWidgets.QLabel("")
		self.wbutxx.setMaximumHeight(12)
		self.hbl2.addWidget(self.wbutxx, 0, 6)
		self.wbutyy = QtWidgets.QLabel("")
		self.wbutyy.setMaximumHeight(12)
		self.hbl2.addWidget(self.wbutyy, 1, 6)

		# buttons for selector use
		if withmodal :
#			self.wspace1 = QtWidgets.QSpacerItem(100, 10, QtWidgets.QSizePolicy.MinimumExpanding)
#			self.hbl2.addSpacerItem(self.wspace1)

			self.wbutcancel = QtWidgets.QPushButton("Cancel")
			self.hbl2.addWidget(self.wbutcancel, 1, 7)

			self.wbutok = QtWidgets.QPushButton("OK")
			self.hbl2.addWidget(self.wbutok, 1, 8)

			self.hbl2.setColumnStretch(6, 1)
			self.hbl2.setColumnStretch(7, 1)
			self.hbl2.setColumnStretch(8, 1)

			self.wbutcancel.clicked.connect(self.buttonCancel)
			self.wbutok.clicked.connect(self.buttonOk)

		self.gbl.addLayout(self.hbl2, 4, 1)

		self.wbutback.clicked.connect(self.buttonBack)
		self.wbutfwd.clicked.connect(self.buttonFwd)
		self.wbutup.clicked.connect(self.buttonUp)
		self.wbutrefresh.clicked.connect(self.buttonRefresh)
		self.wbutinfo.clicked.connect(self.buttonInfo)
		self.selectall.clicked.connect(self.selectAll)
		self.wtree.clicked[QtCore.QModelIndex].connect(self.itemSel)
		self.wtree.activated[QtCore.QModelIndex].connect(self.itemActivate)
		self.wtree.doubleClicked[QtCore.QModelIndex].connect(self.itemDoubleClick)
		self.wtree.expanded[QtCore.QModelIndex].connect(self.itemExpand)
		self.wpath.returnPressed.connect(self.editPath)
		self.wbookmarks.actionTriggered[QtWidgets.QAction].connect(self.bookmarkPress)
		self.wfilter.currentIndexChanged[int].connect(self.editFilter)

		self.setsmode = setsmode	# The sets mode is used when selecting bad particles
		self.curmodel = None	# The current data model displayed in the tree
		self.curpath = None	# The path represented by the current data model
		self.curft = None		# a fileType instance for the currently highlighted object
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
		self.viewhist = []

		# These items are used to do gradually filling in of file details for better interactivity
		self.updtimer = QTimer()		# This causes the actual display updates, which can't be done from a python thread
		self.updtimer.timeout.connect(self.updateDetailsDisplay)
		self.updthreadexit = False		# when set, this triggers the update thread to exit
		self.updthread = threading.Thread(target = self.updateDetails)	# The actual thread
		self.updlist = []				# List of QModelIndex items in need of updating
		self.redrawlist = []			# List of QModelIndex items in need of redisplay
		self.needresize = 0			# Used to resize column widths occasionally
		self.expanded = set()			# We get multiple expand events for each path element, so we need to keep track of which ones we've updated

		self.setPath(startpath)	# start in the local directory
		self.updthread.start()
		self.updtimer.start(200)

		self.result = None			# used in modal mode. Holds final selection

		E2loadappwin("e2display", "main", self)

	def busy(self) :
		"""display a busy cursor"""
		QtWidgets.qApp.setOverrideCursor(Qt.BusyCursor)

	def notbusy(self) :
		"""normal arrow cursor"""
		QtWidgets.qApp.setOverrideCursor(Qt.ArrowCursor)

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
		"""Since we can't do GUI updates from a thread, this is a timer event to update the display after the background thread
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
		rid=[(r[0],i) for i,r in enumerate(rdr)]
		rmin=rdr[min(rid)[1]]; rmax=rdr[max(rid)[1]]
		self.curmodel.dataChanged.emit(self.curmodel.createIndex(rmin[0], 0, rmin[1]),
				 self.curmodel.createIndex(rmax[0], 5, rmax[1]))

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
# 		QtWidgets.QTreeView.keyPressEvent(self.wtree, event)
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
				if os.path.exists(obj.truepath()):
					self.curft = ftc(obj.path())
				else:
					print("Error: file {} does not exist...".format(obj.path()))
					return
				
#				if self.setsmode : self.curft.setSetsDB(re.sub(r'_ctf_flip$|_ctf_wiener$', '', obj.path()))	# If we want to enable bad particle picking (treat ctf and raw data bads as same)

				self.curactions = self.curft.actions()

#				print actions

				for i, b in enumerate(self.wbutmisc) :
					try :
						b.setText(self.curactions[i][0])
						b.setToolTip(self.curactions[i][1])
#						b.show()
						b.setEnabled(True)
					except :
						b.setText("-")
						b.setToolTip("")
						b.setEnabled(False)
			else :
				self.curft = None
				self.curactions = []

				for i, b in enumerate(self.wbutmisc) :
					b.setText("-")
					b.setToolTip("")
					b.setEnabled(False)

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
		for i in range(self.curmodel.rowCount(qmi)-1, -1, -1) :
			self.updlist.append(self.curmodel.index(i, 0, qmi))

	def buttonMisc(self, num) :
		"""One of the programmable action buttons was pressed"""
#		print "press ", self.curactions[num][0]

		try: self.curactions[num][2](self)				# This calls the action method
		except: 
			traceback.print_exc()
			# sometimes we are missing an action

	def buttonOk(self,tog=None) :
		"""When the OK button is pressed, this will emit a signal. The receiver should call the getResult method (once) to get the list of paths"""
		qism = self.wtree.selectionModel().selectedRows()
		self.result = [i.internalPointer().path().replace(os.getcwd(), ".") for i in qism]
		self.updtimer.stop()
		self.ok.emit() # this signal is important when e2ctf is being used by a program running its own eve

	def buttonCancel(self,tog=None) :
		"""When the Cancel button is pressed, a signal is emitted, but getResult should not be called."""
		self.result = []
		self.updtimer.stop()
		self.cancel.emit() # this signal is important when e2ctf is being used by a program running its own eve
		self.close()

	def selectAll(self) :
		self.wtree.selectAll()

	def buttonBack(self,tog=None) :
		"""Button press"""
		# I don't like the stack idea, it's annoying, so I am using a circular array instead John F
		l = self.pathstack.index(self.curpath)
		self.setPath(self.pathstack[(l-1)], True)

		if l == 1 :
			self.wbutback.setEnabled(False)
			self.wbutfwd.setEnabled(True)

	def buttonFwd(self,tog=None) :
		"""Button press"""
		# I don't like the stack idea, it's annoying, so I am using a circular array instead John F
		l = self.pathstack.index(self.curpath)
		self.setPath(self.pathstack[(l+1)], True)

		if l == len(self.pathstack) - 2 :
			self.wbutback.setEnabled(True)
			self.wbutfwd.setEnabled(False)

	def buttonUp(self,tog=None) :
		"""Button press"""
		if "/" in self.curpath : newpath = self.curpath.rsplit("/", 1)[0]
		else : newpath = os.path.realpath(self.curpath).rsplit("/", 1)[0]

		print("Newpath:", newpath)

		#if len(newpath) > 1 : self.setPath(newpath)	# What if we want to return to CWD, '.' # John F

		if len(newpath) > 0 : self.setPath(newpath)

	def buttonRefresh(self,tog=None) :
		"""Button press"""
		self.setPath(self.curpath)

	def infowinClosed(self) :
		self.wbutinfo.setChecked(False)

	def buttonInfo(self,tog=None) :
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
			self.infowin.winclosed.connect(self.infowinClosed)
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

		for i in range(len(self.result)) :
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

		path = os.path.expandvars(path)

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
		elif filt == "?"  or  filt.lower() == "help" :
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

			QtWidgets.QMessageBox.warning(None, "Info", hlp)
		else :
			try :
				flt = re.compile(filt)
				filt = flt
			except :
				filt = filt
#				QtWidgets.QMessageBox.warning(self, "Error", "Bad filter expression")

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
#				QtWidgets.QMessageBox.warning(None, "Error", "Filtering not allowed.")
				print("Filtering is not implemented in this instance of the file browser.")
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
		for i in range(self.curmodel.rowCount(None)-1, -1, -1) :
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
		self.setPath(action.data())
#		self.wtree.setSelectionModel(myQItemSelection(self.curmodel))

	def closeEvent(self, event) :
		E2saveappwin("e2display", "main", self)
		self.updthreadexit = True

		for w in self.view2d+self.view2ds+self.view3d+self.viewplot2d+self.viewplot3d+self.viewhist :
			w.close()

		if self.infowin != None :
			self.infowin.close()

		event.accept()
		self.updtimer.stop()
		# self.app().close_specific(self)
		self.module_closed.emit()

# This is just for testing, of course
def main():
	#em_app = EMApp()
	window = EMBrowserWidget(withmodal = True, multiselect = True)

	window.show()
	em_app.execute()


if __name__ == '__main__' :
	main()
