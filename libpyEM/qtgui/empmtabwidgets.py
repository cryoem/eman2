#!/usr/bin/env python
#
# Author: John Flanagan Dec 1st 2011 (jfflanag@bcm.edu)
# Copyright (c) 2000-2011 Baylor College of Medicine
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

# These classes are subclasses of EMBrowserWidget to provide additonal models to represent various types of files in the GUI.

from EMAN2 import *
import os
import re
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import Qt
from embrowser import EMBrowserWidget, EMFileItemModel, EMDirEntry, nonone, safe_int,safe_float


class EMRefine2dTable(EMBrowserWidget):
	""" Widgewt to display refinement directories """
	def __init__(self, withmodal=False, multiselect=False):
		EMBrowserWidget.__init__(self, withmodal=withmodal, multiselect=multiselect, startpath=".")

	def setPath(self,path,silent=False):
		super(EMRefine2dTable, self).setPath(path,silent=silent,inimodel=EMRefine2dModel)

class EMRefine2dModel(EMFileItemModel):
	""" Item model for the refinement data """
	def __init__(self, startpath=None, dirregex=None):
		regex = re.compile('^r2d')
		EMFileItemModel.__init__(self, startpath=startpath, direntryclass=EMDirEntry, dirregex=regex)

########################################################################################################################

class EMValidateTable(EMBrowserWidget):
	""" Widgewt to display refinement directories """
	def __init__(self, withmodal=False, multiselect=False):
		EMBrowserWidget.__init__(self, withmodal=withmodal, multiselect=multiselect, startpath=".")

	def setPath(self,path,silent=False):
		super(EMValidateTable, self).setPath(path,silent=silent,inimodel=EMValidateModel)

class EMValidateModel(EMFileItemModel):
	""" Item model for the refinement data """
	def __init__(self, startpath=None, dirregex=None):
		regex = re.compile('^TiltValidate')
		EMFileItemModel.__init__(self, startpath=startpath, direntryclass=EMDirEntry, dirregex=regex)

########################################################################################################################

class EMRefineTable(EMBrowserWidget):
	""" Widgewt to display refinement directories """
	def __init__(self, withmodal=False, multiselect=False):
		EMBrowserWidget.__init__(self, withmodal=withmodal, multiselect=multiselect, startpath=".")

	def setPath(self,path,silent=False):
		super(EMRefineTable, self).setPath(path,silent=silent,inimodel=EMRefineModel)

class EMRefineModel(EMFileItemModel):
	""" Item model for the refinement data """
	def __init__(self, startpath=None, dirregex=None):
		regex = re.compile('^refine|^frealign')
		EMFileItemModel.__init__(self, startpath=startpath, direntryclass=EMDirEntry, dirregex=regex)

########################################################################################################################

class EMModelsTable(EMBrowserWidget):
	""" Widget to display junk from e2boxercache """
	def __init__(self, withmodal=False, multiselect=False):
		EMBrowserWidget.__init__(self, withmodal=withmodal, multiselect=multiselect, startpath="./initial_models")

	def setPath(self,path,silent=False):
		super(EMModelsTable, self).setPath(path,silent=silent,inimodel=EMModelsModel)


class EMModelsModel(EMFileItemModel):
	""" Item model for the raw data """

	headers=("Row","3D Models","Quality", "Dims")

	def __init__(self, startpath=None, dirregex=None):
		EMFileItemModel.__init__(self, startpath=startpath, direntryclass=EMModelsEntry, dirregex=dirregex)

	def columnCount(self,parent):
		"Always 4 columns"
		#print "EMFileItemModel.columnCount()=6"
		return 4

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
		elif col==2 :
			if data.quality==0 : return "-"
			return nonone(data.quality)
		elif col==3 :
			if data.dims==0 : return "-"
			return nonone(data.dims)

class EMModelsEntry(EMDirEntry):
	""" Subclassing of EMDirEntry to provide functionality"""

	col=(lambda x:int(x.index), lambda x:x.name,lambda x:safe_int(x.quality), lambda x:x.dims)

	def __init__(self,root,name,i,parent=None,hidedot=True,dirregex=None):
		EMDirEntry.__init__(self,root,name,i,parent=parent,hidedot=hidedot,dirregex=dirregex)
		self.quality=None
		self.dims=None

	def fillDetails(self):
		"""Fills in the expensive metadata about this entry. Returns False if no update was necessary."""
		if self.filetype!=None : return 0		# must all ready be filled in

		cache=None
		cachename=self.name+"!models"
		try:
			cache=js_open_dict(self.root+"/.browsercache.json")
			self.updtime,self.dims,self.filetype,self.nimg,self.quality=cache[cachename]		# try to read the cache for the current file
			if self.updtime>=os.stat(self.truepath()).st_mtime : return 2 		# current cache, no further update necessary
		except:
			pass

		# Should only be this:
		self.filetype="Image"

		# get image counts Needed for get info. Requiring this in in a reimplemeted function is a bad design... I din't do this :)
		try:
			self.nimg = EMUtil.get_image_count(self.path())
		except:
			pass

		# Get particle stack headers
		a = None
		try:
			a = EMData(self.path(),0,True)
		except:
			pass
		if a:
			self.dims = "%dx%dx%d"%(a.get_xsize(),a.get_ysize(),a.get_zsize())
			try:
				self.quality = a['quality']
			except:
				pass

			if cache!=None : cache.setval(cachename,(time.time(),self.dims,self.filetype,self.nimg,self.quality),True)

		return 1

################################################################################################################

class EMSetsTable(EMBrowserWidget):
	""" Widget to display junk from e2boxercache """
	def __init__(self, withmodal=False, multiselect=False):
		EMBrowserWidget.__init__(self, withmodal=withmodal, multiselect=multiselect, startpath="./sets")

	def setPath(self,path,silent=False):
		super(EMSetsTable, self).setPath(path,silent=silent,inimodel=EMSetsModel)


class EMSetsModel(EMFileItemModel):
	""" Item model for the raw data """

	headers=("Row","Name","Num Particles", "Dims")

	def __init__(self, startpath=None, dirregex=None):
		EMFileItemModel.__init__(self, startpath=startpath, direntryclass=EMSetsEntry, dirregex=dirregex)

	def columnCount(self,parent):
		"Always 4 columns"
		#print "EMFileItemModel.columnCount()=6"
		return 4

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
		elif col==2 :
			if data.nimg==0 : return "-"
			return nonone(data.nimg)
		elif col==3 :
			if data.dims==0 : return "-"
			return nonone(data.dims)

class EMSetsEntry(EMDirEntry):
	""" Subclassing of EMDirEntry to provide functionality"""

	col=(lambda x:int(x.index),lambda x:x.name,lambda x:safe_int(x.partcount), lambda x:x.dims)

	def __init__(self,root,name,i,parent=None,hidedot=True,dirregex=None):
		EMDirEntry.__init__(self,root,name,i,parent=parent,hidedot=hidedot,dirregex=dirregex)
		self.dims=None

	def fillDetails(self):
		"""Fills in the expensive metadata about this entry. Returns False if no update was necessary."""
		if self.filetype!=None : return False		# must all ready be filled in

		cache=None
		cachename=self.name+"!sets"
		try:
			cache=js_open_dict(self.root+"/.browsercache.json")
			self.updtime,self.dims,self.filetype,self.nimg=cache[cachename]		# try to read the cache for the current file
			if self.updtime>=os.stat(self.truepath()).st_mtime : return 2 		# current cache, no further update necessary
		except:
			pass

		# get image counts
		try:
			self.nimg = EMUtil.get_image_count(self.path())
			if self.nimg == 1 : self.filetype="Image"
			if int(self.nimg) > 1 : self.filetype="Image Stack"
		except:
			self.filetype="-"

		# Get particle stack headers
		a = None
		try:
			a = EMData(self.path(),0,True)
		except:
			pass
		if a:
			self.dims = "%dx%dx%d"%(a.get_xsize(),a.get_ysize(),a.get_zsize())

			if cache!=None : cache.setval(cachename,(time.time(),self.dims,self.filetype,self.nimg),True)


		return True

###########################################################################################################################

class EMParticlesTable(EMBrowserWidget):
	""" Widget to display junk from e2boxercache """
	def __init__(self, withmodal=False, multiselect=False):
		EMBrowserWidget.__init__(self, withmodal=withmodal, multiselect=multiselect, startpath="./particles")

	def setPath(self,path,silent=False):
		super(EMParticlesTable, self).setPath(path,silent=silent,inimodel=EMParticlesModel)

class EMParticlesModel(EMFileItemModel):
	""" Item model for the raw data """

	headers=("Row","Raw Data Files","Type", "Num Particles", "Particle Dims", "Quality")

	def __init__(self, startpath=None, dirregex=None):
		EMFileItemModel.__init__(self, startpath=startpath, direntryclass=EMParticlesEntry, dirregex=dirregex)

	def columnCount(self,parent):
		"Always 6 columns"
		#print "EMFileItemModel.columnCount()=6"
		return 6

	def data(self,index,role):
		"Returns the data for a specific location as a string"

		if not index.isValid() : return None
		if role!=Qt.DisplayRole : return None

		data=index.internalPointer()
		if data==None :
			print "Error with index ",index.row(),index.column()
			return "Rubbish"
		#if index.column()==0 : print "EMFileItemModel.data(%d %d %s)=%s"%(index.row(),index.column(),index.parent(),str(data.__dict__))

		col=index.column()
		if col==0:
			return nonone(data.index)
		elif col==1 :
			if data.isbdb : return "bdb:"+data.name
			return nonone(data.name)
		elif col==2 :
			if data.typ==0 : return "-"
			return nonone(data.typ)
		elif col==3 :
			if data.nimg==0 : return "-"
			return nonone(data.nimg)
		elif col==4 :
			if data.particledim==0 : return "-"
			return nonone(data.particledim)
		elif col==5 :
			if data.quality==0 : return "-"
			return nonone(data.quality)

class EMParticlesEntry(EMDirEntry):
	""" Subclassing of EMDirEntry to provide functionality"""

	col=(lambda x:int(x.index),lambda x:x.name,lambda x:x.typ,lambda x:safe_int(x.nimg), lambda x:x.particledim, lambda x:safe_int(x.quality))

	def __init__(self,root,name,i,parent=None,hidedot=True,dirregex=None):
		EMDirEntry.__init__(self,root,name,i,parent=parent,hidedot=hidedot,dirregex=dirregex)
		self.typ = None
		self.particledim=None
		self.defocus=None
		self.bfactor=None
		self.snr=None
		self.quality=None
		self.sampling=None

	def fillDetails(self):
		"""Fills in the expensive metadata about this entry. Returns False if no update was necessary."""
		if self.filetype!=None : return False		# must all ready be filled in

		# get quality


		# Check the cache for metadata
		cache=None
		cachename=self.name+"!particles"
		try:
			cache=js_open_dict(self.root+"/.browsercache.json")
			self.updtime,self.particledim,self.filetype,self.nimg,self.typ,self.quality=cache[cachename]		# try to read the cache for the current file
			old=self.cache_old()
			if old==0 : return 2 		# current cache, no further update necessary
		except:
			old=3

		if old & 2 :
			try:
				qdb=js_open_dict(info_name(self.path()))
				self.quality=qdb["quality"]
			except: self.quality="-"

		if old & 1:
		# get image counts
			try:
				self.nimg = EMUtil.get_image_count(self.path())
				if self.nimg==1 : self.filetype="Image"
				if self.nimg > 1 : self.filetype="Image Stack"

			except:
				pass

			# Get particle stack headers
			a = None
			try:
				a = EMData(self.path(),0,True)
			except:
				pass
			if a:
				self.particledim = a.get_xsize()

			# get partivle set type
			try:
				self.typ = str(self.path().split('_ctf_')[1])
			except:
				pass

		if cache!=None : cache.setval(cachename,(time.time(),self.particledim,self.filetype,self.nimg,self.typ,self.quality),True)

		return 1

###########################################################################################################################

class EMCTFParticlesTable(EMBrowserWidget):
	""" Widget to display junk from e2boxercache """
	def __init__(self, withmodal=False, multiselect=False):
		EMBrowserWidget.__init__(self, withmodal=withmodal, multiselect=multiselect, startpath="./particles")

	def setPath(self,path,silent=False):
		super(EMCTFParticlesTable, self).setPath(path,silent=silent,inimodel=EMCTFParticlesModel)

class EMCTFParticlesModel(EMFileItemModel):
	""" Item model for the raw data """

	headers=("Row","Raw Data Files","Type", "Num Particles", "Particle Dims", "Defocus", "B Factor", "SNR-lo","SNR-hi", "Quality", "Sampling")

	def __init__(self,startpath=None, direntryclass=None, dirregex=None):
		if not direntryclass:
			EMFileItemModel.__init__(self, startpath=startpath, direntryclass=EMCTFParticlesEntry, dirregex=dirregex)
		else:
			EMFileItemModel.__init__(self, startpath=startpath, direntryclass=direntryclass, dirregex=dirregex)

	def columnCount(self,parent):
		"Always 11 columns"
		#print "EMFileItemModel.columnCount()=6"
		return 11

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
		elif col==2 :
			if data.typ==0 : return "-"
			return nonone(data.typ)
		elif col==3 :
			if data.nimg==0 : return "-"
			return nonone(data.nimg)
		elif col==4 :
			if data.particledim==0 : return "-"
			return nonone(data.particledim)
		elif col==5 :
			if data.defocus==0 : return "-"
			return nonone(data.defocus)
		elif col==6 :
			if data.bfactor==0 : return "-"
			return nonone(data.bfactor)
		elif col==7 :
			if data.snr==0 : return "-"
			return nonone(data.snr)
		elif col==8 :
			if data.snrhi==0 : return "-"
			return nonone(data.snrhi)
		elif col==9 :
			if data.quality==0 : return "-"
			return nonone(data.quality)
		elif col==10 :
			if data.sampling==0 : return "-"
			return nonone(data.sampling)

class EMCTFParticlesEntry(EMDirEntry):
	""" Subclassing of EMDirEntry to provide functionality"""

	col=(lambda x:int(x.index),lambda x:x.name,lambda x:x.typ,lambda x:safe_int(x.nimg), lambda x:x.particledim, lambda x:safe_float(x.defocus), lambda x:safe_float(x.bfactor), lambda x:safe_float(x.snr), lambda x:safe_float(x.snrhi),lambda x:safe_int(x.quality), lambda x:x.sampling)

	def __init__(self,root,name,i,parent=None,hidedot=True,dirregex=None):
		EMDirEntry.__init__(self,root,name,i,parent=parent,hidedot=hidedot,dirregex=dirregex)
		self.typ = None
		self.particledim=None
		self.defocus=None
		self.bfactor=None
		self.snr=None
		self.snrhi=None
		self.quality=None
		self.sampling=None

	def fillDetails(self):
		"""Fills in the expensive metadata about this entry. Returns False if no update was necessary."""
		if self.filetype!=None : return False		# must all ready be filled in

		# Check the cache for metadata
		cache=None
		self.updtime=0
		cachename=self.name+"!ctf"
		try:
			cache=js_open_dict(self.root+"/.browsercache.json")
			self.updtime,self.particledim,self.filetype,self.nimg,self.typ,self.defocus,self.bfactor,self.sampling,self.snr,self.snrhi,self.quality,self.badparticlecount=cache[cachename]		# try to read the cache for the current file
			old=self.cache_old()
			if old==0 : return 2 		# current cache, no further update necessary
		except:
			old=3

		if old&2:
			try:
				db=js_open_dict(info_name(self.truepath()))
				self.quality=db["quality"]
			except:
				self.quality="-"
			try:
				self.badparticlecount=len(db["select_bad"])		# Note that this is used in a subclass, but not actually used here
				print self.badparticlecount
			except:
				self.badparticlecount=0
			try:
				ctf=db["ctf"][0]
				if ctf.dfdiff==0 : self.defocus = "%.3f" % ctf.defocus
				else : self.defocus = "%.3f (%.3f, %.1f)" % (ctf.defocus,ctf.dfdiff,ctf.dfang)
				self.bfactor = "%.3f" % ctf.bfactor
				self.sampling = str(len(ctf.snr))
				s0=max(2,int(1.0/(200.0*ctf.dsbg)))
				s1=max(2,int(1.0/(20.0*ctf.dsbg)))
				s2=int(1.0/(10.0*ctf.dsbg))
				s3=int(1.0/(4.0*ctf.dsbg))
				self.snr = "%.3f" %  (sum(ctf.snr[s0:s1])/(s1-s0))
				self.snrhi = "%.3f" %  (sum(ctf.snr[s2:s3])/(s3-s2))
			except:
				self.defocus="-"
				self.bfactor="-"
				self.sampling="-"
				self.snr="-"
				self.snrhi="-"

		if old&1 :
			# get image counts
			try:
				self.nimg = EMUtil.get_image_count(self.path())
				if self.nimg==1 : self.filetype="Image"
				if self.nimg > 1 : self.filetype="Image Stack"

			except:
				pass

			# Get particle stack headers
			a = None
			try:
				a = EMData(self.path(),0,True)
			except:
				pass
			if a:
				self.particledim = a.get_xsize()

			# get partivle set type
			try:
				self.typ = str(self.path().split('_ctf_')[1])
			except:
				pass

		# Update the cache
		if cache!=None and old!=0: cache.setval(cachename,(time.time(),self.particledim,self.filetype,self.nimg,self.typ,self.defocus,self.bfactor,self.sampling,self.snr,self.snrhi,self.quality,self.badparticlecount),True)

		return True

#####################################################################################################################

class EMCTFcorrectedParticlesTable(EMCTFParticlesTable):
	""" Widget to display junk from e2boxercache. Same as EMCTFParticlesTable, but only displays CTF corrected  """
	def __init__(self, withmodal=False, multiselect=False):
		EMCTFParticlesTable.__init__(self, withmodal=withmodal, multiselect=multiselect)

	def setPath(self,path,silent=False):
		super(EMCTFParticlesTable, self).setPath(path,silent=silent,inimodel=EMCTFcorrectedParticlesModel)

class EMCTFcorrectedParticlesModel(EMCTFParticlesModel):
	""" Item model for the raw data """
	def __init__(self, startpath=None, dirregex=None):
		EMCTFParticlesModel.__init__(self, startpath=startpath, direntryclass=EMCTFcorrectedParticlesEntry, dirregex=dirregex)

class EMCTFcorrectedParticlesEntry(EMCTFParticlesEntry):
	""" Subclassing of EMDirEntry to provide functionality"""
	def __init__(self,root,name,i,parent=None,hidedot=True,dirregex=None):
		EMCTFParticlesEntry.__init__(self,root,name,i,parent=None,hidedot=True,dirregex=dirregex)

######################################################################################################################

class EMParticlesEditTable(EMBrowserWidget):
	""" Widget to display junk from e2boxercache """
	def __init__(self, withmodal=False, multiselect=False):
		EMBrowserWidget.__init__(self, withmodal=withmodal, multiselect=multiselect, startpath="./particles",setsmode=True)
		self.wfilter.setCurrentIndex(1)

	def setPath(self,path,silent=False):
		super(EMParticlesEditTable, self).setPath(path,silent=silent,inimodel=EMParticlesEditModel)

class EMParticlesEditModel(EMFileItemModel):
	""" Item model for the raw data """

	headers=("Row","Filename", "# Ptcl", "Bad Ptcl", "Defocus", "B Factor", "SNR-lo", "SNR-hi", "Quality", "Sampling", "Size")

	def __init__(self,startpath=None,dirregex=None):
		EMFileItemModel.__init__(self, startpath=startpath, direntryclass=EMParticlesEditEntry,dirregex=dirregex)

	def columnCount(self,parent):
		"Always 11 columns"
		#print "EMFileItemModel.columnCount()=6"
		return 11

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
		elif col==2 :
			if data.nimg==0 : return "-"
			return nonone(data.nimg)
		elif col==3 :
			if data.badparticlecount==0 : return "-"
			return nonone(data.badparticlecount)
		elif col==4 :
			if data.defocus==0 : return "-"
			return nonone(data.defocus)
		elif col==5 :
			if data.bfactor==0 : return "-"
			return nonone(data.bfactor)
		elif col==6 :
			if data.snr==0 : return "-"
			return nonone(data.snr)
		elif col==7 :
			if data.snrhi==0 : return "-"
			return nonone(data.snrhi)
		elif col==8 :
			if data.quality==0 : return "-"
			return nonone(data.quality)
		elif col==9 :
			if data.sampling==0 : return "-"
			return nonone(data.sampling)
		elif col==10 :
			if data.particledim==0 : return "-"
			return nonone(data.particledim)

class EMParticlesEditEntry(EMCTFParticlesEntry):
	""" Subclassing of EMDirEntry to provide functionality"""

	col=(lambda x:int(x.index),lambda x:x.name,lambda x:safe_int(x.nimg), lambda x:safe_int(x.badparticlecount), lambda x:safe_float(x.defocus), lambda x:safe_float(x.bfactor), lambda x:safe_float(x.snr), lambda x:safe_float(x.snrhi), lambda x:safe_int(x.quality), lambda x:x.sampling, lambda x:x.particledim)

	def __init__(self,root,name,i,parent=None,hidedot=True,dirregex=None):
		EMCTFParticlesEntry.__init__(self,root=root,name=name,i=i,parent=parent,hidedot=hidedot,dirregex=dirregex)
		self.badparticlecount = None

	def fillDetails(self):
		super(EMParticlesEditEntry, self).fillDetails()

#######################################################################################################################

class EMBoxesTable(EMBrowserWidget):
	""" Widget to display junk from e2boxercache """
	def __init__(self, withmodal=False, multiselect=False):
		EMBrowserWidget.__init__(self, withmodal=withmodal, multiselect=multiselect, startpath="./micrographs")

	def setPath(self,path,silent=False):
		super(EMBoxesTable, self).setPath(path,silent=silent,inimodel=EMBoxesModel)


class EMBoxesModel(EMFileItemModel):
	""" Item model for the raw data """

	headers=("Row","Raw Data Files","Stored Boxes", "Quality")

	def __init__(self,startpath=None,dirregex=None):
		EMFileItemModel.__init__(self, startpath=startpath, direntryclass=EMBoxesEntry,dirregex=dirregex)

	def columnCount(self,parent):
		"Always 4 columns"
		#print "EMFileItemModel.columnCount()=6"
		return 4

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
		elif col==2 :
			if data.boxcount==0 : return "-"
			return nonone(data.boxcount)
		elif col==3 :
			if data.quality==0 : return "-"
			return nonone(data.quality)

class EMBoxesEntry(EMDirEntry):
	""" Subclassing of EMDirEntry to provide functionality"""

	col=(lambda x:int(x.index),lambda x:x.name,lambda x:safe_int(x.boxcount), lambda x:safe_int(x.quality), lambda x:safe_int(x.mgquality))

	def __init__(self,root,name,i,parent=None,hidedot=True,dirregex=None):
		EMDirEntry.__init__(self,root,name,i,parent=parent,hidedot=hidedot,dirregex=dirregex)
		self.boxcount=None
		self.quality=None
		self.mgquality=None

	def fillDetails(self):
		"""Fills in the expensive metadata about this entry. Returns False if no update was necessary."""
		if self.filetype!=None : return 0		# must all ready be filled in


		# Check the cache for metadata
		cache=None
		cachename=self.name+"!boxes"
		try:
			cache=js_open_dict(self.root+"/.browsercache.json")
			self.updtime,self.boxcount,self.filetype,self.quality=cache[cachename]		# try to read the cache for the current file
			old=self.cache_old()
			if old==0 : return 2 		# current cache, no further update necessary
		except:
			old=3

		if old&2 :
			try:
				db=js_open_dict(info_name(self.path()))
				self.quality=db["quality"]
			except: self.quality="-"

			try:
				boxdb=db["boxes"]
				self.boxcount=len(boxdb)
			except:
				self.boxcount="-"

		if old&1 :
			# get image counts
			try:
				self.nimg = EMUtil.get_image_count(self.path())
				if self.nimg==1 : self.filetype="Image"
			except:
				self.filetype="-"

		# Set cache
		if cache!=None : cache.setval(cachename,(time.time(),self.boxcount,self.filetype,self.quality),True)

		return 1

#############################################################################################################################

class EMRCTBoxesTable(EMBrowserWidget):
	""" Widget to display junk from e2boxercache """
	def __init__(self, withmodal=False, multiselect=False):
		EMBrowserWidget.__init__(self, withmodal=withmodal, multiselect=multiselect, startpath="./micrographs")

	def setPath(self,path,silent=False):
		super(EMRCTBoxesTable, self).setPath(path,silent=silent,inimodel=EMRCTBoxesModel)


class EMRCTBoxesModel(EMFileItemModel):
	""" Item model for the raw data """

	headers=("Row","Raw Data Files","Stored Boxes", "Quality")

	def __init__(self,startpath=None,dirregex=None):
		EMFileItemModel.__init__(self, startpath=startpath, direntryclass=EMRCTBoxesEntry,dirregex=dirregex)

	def columnCount(self,parent):
		"Always 4 columns"
		#print "EMFileItemModel.columnCount()=6"
		return 4

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
		elif col==2 :
			if data.boxcount==0 : return "-"
			return nonone(data.boxcount)
		elif col==3 :
			if data.quality==0 : return "-"
			return nonone(data.quality)

class EMRCTBoxesEntry(EMDirEntry):
	""" Subclassing of EMDirEntry to provide functionality"""

	col=(lambda x:int(x.index),lambda x:x.name,lambda x:safe_int(x.boxcount), lambda x:safe_int(x.quality), lambda x:safe_int(x.mgquality))

	def __init__(self,root,name,i,parent=None,hidedot=True,dirregex=None):
		EMDirEntry.__init__(self,root,name,i,parent=parent,hidedot=hidedot,dirregex=dirregex)
		self.boxcount=None
		self.quality=None
		self.mgquality=None

	def fillDetails(self):
		"""Fills in the expensive metadata about this entry. Returns False if no update was necessary."""
		if self.filetype!=None : return 0		# must all ready be filled in


		# Check the cache for metadata
		cache=None
		cachename=self.name+"!rct"
		try:
			cache=js_open_dict(self.root+"/.browsercache.json")
			self.updtime,self.boxcount,self.filetype,self.quality=cache[cachename]		# try to read the cache for the current file
			old=self.cache_old()
			if old==0 : return 2 		# current cache, no further update necessary
		except:
			old=3

		if old&2 :
			try:
				db=js_open_dict(info_name(self.path()))
				self.quality=db["quality"]
			except: self.quality="-"

			try:
				boxdb=db["boxes_tilted"]
				self.boxcount=len(boxdb)
				self.filetype="Tilted"
			except:
				try:
					boxdb=db["boxes_untilted"]
					self.boxcount=len(boxdb)
					self.filetype="Untilted"
				except:
					self.boxcount="-"


		if old&1 :
			# get image counts
			if self.filetype in (None,"-") :
				try:
					self.nimg = EMUtil.get_image_count(self.path())
					if self.nimg==1 : self.filetype="Image"
				except:
					self.filetype="-"

		# Set cache
		if cache!=None : cache.setval(cachename,(time.time(),self.boxcount,self.filetype,self.quality),True)

		return 1

#################################################################################################################################

class EMSubTomosTable(EMBrowserWidget):
	""" Widget to display Raw Data """
	def __init__(self, withmodal=False, multiselect=False):
		EMBrowserWidget.__init__(self, withmodal=withmodal, multiselect=multiselect, startpath="./subtomograms")

	def setPath(self,path,silent=False):
		super(EMSubTomosTable, self).setPath(path,silent=silent,inimodel=EMSubTomosModel)

class EMSubTomosModel(EMFileItemModel):
	""" Item model for the raw data """

	headers=("Row","Subtomograms","Num Subtomos", "Dims")

	def __init__(self,startpath=None,dirregex=None):
		EMFileItemModel.__init__(self, startpath=startpath, direntryclass=EMSubTomosEntry,dirregex=dirregex)

	def columnCount(self,parent):
		"Always 4 columns"
		#print "EMFileItemModel.columnCount()=6"
		return 4

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
		elif col==2 :
			if data.nimg==0 : return "-"
			return nonone(data.nimg)
		elif col==3 :
			if data.dim==0 : return "-"
			return nonone(data.dim)

class EMSubTomosEntry(EMDirEntry):
	""" Subclassing of EMDirEntry to provide functionality"""

	col=(lambda x:int(x.index),lambda x:x.name,lambda x:safe_int(x.nimg), lambda x:x.dim)

	def __init__(self,root,name,i,parent=None,hidedot=True,dirregex=None):
		EMDirEntry.__init__(self,root,name,i,parent=parent,hidedot=hidedot,dirregex=dirregex)

	def fillDetails(self):
		# Maybe add code to cache results.....
		super(EMSubTomosEntry, self).fillDetails()

#################################################################################################################################

class EMRawDataTable(EMBrowserWidget):
	""" Widget to display Raw Data """
	def __init__(self, withmodal=False, multiselect=False, startpath="./micrographs"):
		EMBrowserWidget.__init__(self, withmodal=withmodal, multiselect=multiselect, startpath=startpath)

	def setPath(self,path,silent=False):
		super(EMRawDataTable, self).setPath(path,silent=silent,inimodel=EMRawDataModel)


class EMRawDataModel(EMFileItemModel):
	""" Item model for the raw data """

	headers=("Row","Raw Data Files","Dimensions", "Quality")

	def __init__(self,startpath=None,dirregex=None):
		EMFileItemModel.__init__(self, startpath=startpath, direntryclass=EMRawDataEntry,dirregex=dirregex)

	def columnCount(self,parent):
		"Always 4 columns"
		#print "EMFileItemModel.columnCount()=6"
		return 4

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
		elif col==2 :
			if data.dim==0 : return "-"
			return nonone(data.dim)
		elif col==3 :
			if data.quality==0 : return "-"
			return nonone(data.quality)

class EMRawDataEntry(EMDirEntry):
	""" Subclassing of EMDirEntry to provide functionality"""

	col=(lambda x:int(x.index),lambda x:x.name,lambda x:x.dim, lambda x:safe_int(x.quality))

	def __init__(self,root,name,i,parent=None,hidedot=True,dirregex=None):
		EMDirEntry.__init__(self,root,name,i,parent=parent,hidedot=hidedot,dirregex=dirregex)
		self.quality = None # init
		self.mgquality = None

	def fillDetails(self):
#		super(EMRawDataEntry, self).fillDetails(db)
		if self.filetype!=None : return False		# must all ready be filled in

		# Check the cache for metadata
		cache=None
		cachename=self.name+"!raw"
		try:
			cache=js_open_dict(self.root+"/.browsercache.json")
			self.updtime,self.filetype,self.dim,self.quality=cache[cachename]		# try to read the cache for the current file
			old=self.cache_old()
			if old==0 : return 2 		# current cache, no further update necessary
		except:
			old=3

		if old & 2 :
			try:
				qdb=js_open_dict(info_name(self.path()))
				self.quality=qdb["quality"]
			except: self.quality="-"

		if old & 1:
			try:
				tmp=EMData(self.path(),0,True)		# try to read an image header for the file
				if tmp["ny"]==1 : self.dim=str(tmp["nx"])
				elif tmp["nz"]==1 : self.dim="%d x %d"%(tmp["nx"],tmp["ny"])
				else : self.dim="%d x %d x %d"%(tmp["nx"],tmp["ny"],tmp["nz"])
				self.filetype="Image"
			except :
				self.filetype="-"

		if cache!=None : cache.setval(cachename,(time.time(),self.filetype,self.dim,self.quality),True)

		return 1


#################################################################################################################################

class EMTomoDataTable(EMBrowserWidget):
	""" Widget to display Raw Data """
	def __init__(self, withmodal=False, multiselect=False):
		EMBrowserWidget.__init__(self, withmodal=withmodal, multiselect=multiselect, startpath="./rawtomograms")

	def setPath(self,path,silent=False):
		super(EMTomoDataTable, self).setPath(path,silent=silent,inimodel=EMTomoDataModel)

class EMTomoDataModel(EMFileItemModel):
	""" Item model for the raw data """

	headers=("Row","Raw Data Files","Dimensions")

	def __init__(self,startpath=None,dirregex=None):
		EMFileItemModel.__init__(self, startpath=startpath, direntryclass=EMTomoDataEntry,dirregex=dirregex)

	def columnCount(self,parent):
		"Always 3 columns"
		#print "EMFileItemModel.columnCount()=6"
		return 3

	def data(self,index,role):
		"Returns the data for a specific location as a string"

		if not index.isValid() : return None
		if role!=Qt.DisplayRole : return None

		data=index.internalPointer()
		if data==None :
			print "Error with index ",index.row(),index.column()
			return "XXX"

		col=index.column()
		if col==0:
			return nonone(data.index)
		elif col==1 :
			if data.isbdb : return "bdb:"+data.name
			return nonone(data.name)
		elif col==2 :
			if data.dim==0 : return "-"
			return nonone(data.dim)

class EMTomoDataEntry(EMDirEntry):
	""" Subclassing of EMDirEntry to provide functionality"""

	col=(lambda x:int(x.index),lambda x:x.name,lambda x:x.dim)

	def __init__(self,root,name,i,parent=None,hidedot=True,dirregex=None):
		EMDirEntry.__init__(self,root,name,i,parent=parent,hidedot=hidedot,dirregex=dirregex)
		self.mgquality = None

	def fillDetails(self):
		# Maybe add code to cache results.....
		super(EMTomoDataEntry, self).fillDetails()



