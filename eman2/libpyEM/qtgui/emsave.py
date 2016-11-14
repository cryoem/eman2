#!/usr/bin/env python

#
# Author: David Woolford (woolford@bcm.edu) April 2009
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

from PyQt4.QtCore import Qt
from PyQt4 import QtGui,QtCore
from EMAN2 import EMData, file_exists, gimme_image_dimensions3D,get_image_directory,EMUtil,base_name,gm_time_string
from EMAN2db import db_check_dict, db_remove_dict
import os
# For example usage see http://blake.bcm.edu/emanwiki/EMAN2ImageFormats#SavingEMDatafromPython

class EMFileTypeValidator:
	'''
	A basic validator class - checks to make sure the file name is valid using the file type
	'''
	def __init__(self,type="img"):
		'''
		@param type a file type string such as "tlt","hdf" etc
		'''
		self.type = type
	
	def validate_file_name(self,file_name):
		'''
		Validate the file name 
		Checks to make sure the file name is valid
		@return 0 if something went wrong, for example, the file is not valid
		'''
		
		msg = QtGui.QMessageBox()
		vals = file_name.split(".")
		if len(vals) < 2 or vals[-1] != self.type:
			# error is therefore a string
			msg.setText("%s is not of type %s" %(file_name,self.type))
			msg.exec_()
			return 0
			
		return 1
	

class EMCoordFileValidator:
	'''
	checks to make sure the file name supplied is readable as a box database, in the traditional EMAN1 sense 
	'''
	def __init__(self): pass
	
	def validate_file_name(self,file_name):
		'''
		Validate the file name 
		Checks to make sure the first two non commented lines of the file defines at least 4 values that can be converted to an int
		@return 0 if something went wrong, for example, the file is not valid
		'''
		
		if not isinstance(file_name,list):
			lst = [file_name]
		else:
			lst = file_name
		for f in lst:
			try:
				fin=file(f)
				fin.seek(0)
				rdata = []
				while (len(rdata) < 2):
					line = fin.readline()
					if line == '': # this is equivalent to EOF
						break
					elif line[0] == '#': continue
					else: rdata.append(line)
					
				#print rdata
				if len(rdata) < 2: raise RuntimeError()
	
				rdata=[[float(j) for j in i.split()] for i in rdata]
				
				if len(rdata[0]) == 0 or len(rdata[1]) == 0: raise RuntimeError()
				if  len(rdata[0]) != len(rdata[1]) : raise RuntimeError()
	
				return 1
				
			except:
				msg = QtGui.QMessageBox()
				msg.setText("%s is not a valid coordinate file" %file_name)
				msg.exec_()
				return 0


def save_data(item_object):
	'''
	Intended to be the only function you need ever call if you want to save EM image data, be it one or many.
	@param item_object can be an EMData, a list of EMDatas, an EMDataListCache, or a list of EMListingItems in the selector
	@return the name of the file that was written to disk. May be and empty string ("") if no writing occurred
	'''
	saver = None
	if EMSingleImageSaveDialog.validate_save_argument(item_object):
		saver = EMSingleImageSaveDialog()
	elif EMStackSaveDialog.validate_save_argument(item_object):
		saver = EMStackSaveDialog()
	else:
		raise NotImplementedError(repr(item_object)+" is not supported by this function")
	
	saver.save(item_object)
	
	return saver.get_file_name()

class LightEMDataSave:
	'''
	Used for file io - never reads the image until you actually call write_image, so if the user hits cancel
	they will not experience any time lags due to file reading (which would have then been in vain)
	'''
	def __init__(self,file_name,idx=0):
		'''
		@exception RuntimeError raised if the file_name does not exist
		'''
		if not file_exists(file_name): raise RuntimeError("%s does not exist" %file_name)
		self.file_name = file_name
		self.idx = idx
		
	def write_image(self,*args):
		a = EMData()
		a.read_image(self.file_name,self.idx)
		a.write_image(*args)

	def get_attr_dict(self):
		a = EMData()
		a.read_image(self.file_name,self.idx,True)
		return a.get_attr_dict()
	
class EMFileSaver:
	'''
	Base class for file savers. This function is tightly linked to the save_data function in this
	file.
	'''
	def __init__(self): pass
	def validate_file_name(self,file_name):
		'''
		Validate the file name and do any internal preparation for the impending call to
		save.
		Preparation includes checking whether any appending or overwriting needs to occur
		@return 1 if ready for save to be called, 0 if not
		'''
		raise NotImplementedError("Inheriting classes must supply the get_file_name function")
	
	def validate_save_argument(self,item_object):
		'''
		A function for verifying that the argument to the save function is acceptable. Used both
		internally in the save function and externally by the file scope save_data function
		MUST be made static by adding this line 
		validate_save_argument = Callable(validate_save_argument)
		@param item_object the argument that will ultimately be given to the save function
		'''
		raise NotImplementedError("Inheriting classes must supply the validate_save_argument function, and also make it static")
	
	def save(self,item_object):
		'''
		Save EMAN2 style image data to disk
		This function is the main saving interface, The item_object can be anything, but
		it should first have been validated by a call to validate_save_argument
		@param item_object an object that has all of the attributes required to save data (this is intentionally a broad definition)
		'''
		raise NotImplementedError("Inheriting classes must supply the save function")
	
	def get_file_name(self):
		'''
		This function should return the name of the file that was saved. If no
		file was saved, it should return ""
		@return the name of the file that was saved
		'''
		raise NotImplementedError("Inheriting classes must supply the get_file_name function")
	
class EMSingleImageSaveDialog(EMFileSaver):
	'''
	A class used for saving single images. Will append if the file exists and the type allows it
	'''
	def __init__(self):
		EMFileSaver.__init__(self)

		self.__item = None # this will be a reference to the item
		self.file_name_used = "" # ultimately stores the name of the file that was stored to disk, if it actually occurred
		self.validator = None
	def validate_save_argument(item):
		'''
		Verify that the function argument will work if supplied to the save function
		@param item the item you want to test
		@return true or False
		'''
		if not isinstance(item,EMData) and not isinstance(item,LightEMDataSave):return False
		return True
	
	# mimic a static function
	validate_save_argument = staticmethod(validate_save_argument)
	
	def save(self,item):
		'''
		Runs a file saving dialog and prompts the user to make sure nothing goes awry
		@param item an EMData object
		@exception RuntimeError raised if item is not of type EMData
		'''
		fine = EMSingleImageSaveDialog.validate_save_argument(item)
		if not fine: raise RuntimeError("item is not an EMData instance")
		
		self.validator = EMSaveImageValidator([item])
		self.__item = item
		from emselector import EMSelectorDialog
		selector = EMSelectorDialog(True,True)
		selector.set_validator(self.validator)
		file = selector.exec_()
		if file != "":
			self.__save_file(file)

	def __save_file(self,file):
		'''
		Called internally to save the file. Checks to make sure it's a valid file name
		and may run a dialog asking if the user wants to overwrite or append data.
		@param file the file name
		@return 0 if there was an error, 1 if there was not
		@exception NotImplementedError raised when the EMFileExistsDialog returns an unknown code
		'''
		
		tmp_file_object = EMDummyTmpFileHandle(file)
		if self.validator.is_overwriting():
			tmp_file_object = EMTmpFileHandle(file)

		out_file = tmp_file_object.get_tmp_file_name()
		try:
			self.__item.write_image(out_file, -1)
		except:
			msg = QtGui.QMessageBox()
			msg.setText("An exception occured while writing %s, please try again" %out_file)
			msg.exec_()
			tmp_file_object.remove_tmp_file()
			return 1
		
		tmp_file_object.finalize_renaming()		
		self.file_name_used = tmp_file_object.get_final_file_name()	
		return 1
	
	def get_file_name(self):
		'''
		Get the file name that was used to write to disk
		Can be empty if file save was cancelled
		'''
		return self.file_name_used

class EMSaveImageValidator:
	'''
	A validator class - checks to make sure the file name is valid
	and stores an overwrite flag. Will trigger a Qt message saying the file name is invalid.
	Will also trigger and overwrite/append dialog if appropriate.
	'''
	def __init__(self,item_list):
		'''
		@param item_list should be a list of EMDatas, and EMDataListCache, or a list of EMListingWidgets
		If you have only one EMData just put it in a list
		'''
		self.overwrite = False # overwriting flag
		self.append = False # appending flag - mostly just for book keeping
		
		# FIXME - Supported types should retrieved from EMUtil (in progress)
		if len(item_list) == 1:
			self.__file_filt = [".hdf", ".img", ".hed",".spi","bdb:",".tif",".mrc",".dm3",".pif"] # allowable file types, be good if this was automatic
		else:
			if item_list[0].get_attr_dict()["nz"] == 1:
				self.__file_filt=[".hdf", ".img",".hed",".spi","bdb:"]
			else:
				self.__file_filt =[".hdf"] # 3D only works for hdf (and bdb)
				
		self.__item_list = item_list
	def is_overwriting(self):
		'''
		Ask whether or not overwriting has been chosen as the write method
		@return True or False
		'''
		return self.overwrite
	
	def __validate_file_name(self,file_name):
		'''
		Checks to see if the file_name is valid in terms of what EMAN2 knows it can write
		@param file_name the file_name that is to be checked
		@return an error string if an error was encountered. Otherwise return None
		'''
		if len(file_name) > 4 and file_name[:4] == "bdb:":
			return None # The user has specified a file in bdb format
			
		splt = file_name.split(".")
		if len(splt) < 2:
			return "The file name you specified is invalid - can't determine the file type"
			
		# the length of splt is at least 2
		if ("."+splt[-1]) not in self.__file_filt:
			 return "The file type you specified: %s, is not valid in this context" %splt[-1]
			
		return None
		
	def validate_file_name(self,file_name):
		'''
		Validate the file name and do any internal preparation for the impending call to
		save.
		Checks to make sure the file name is valid
		Confirms overwrite and asks to append
		@return 0 if something went wrong (e.g. the user cancelled), 1 if it's okay to call save now
		'''
		
		msg = QtGui.QMessageBox()
		error = self.__validate_file_name(file_name)
		if error != None:
			# error is therefore a string
			msg.setText(error)
			msg.exec_()
			return 0
		
		# If we make it here the file name is fine, but it may exist
		self.overwrite = False
		self.append = False
		if file_exists(file_name):
			file_exists_dialog = EMFileExistsDialog(file_name,self.__item_list)
			code = file_exists_dialog.exec_()
			if code == 0:
				return 0
			elif code == 1:
				self.overwrite = True
			elif code == 2:
				self.append = True
			else: raise NotImplementedError("A developer has modified the code and this has caused an error. Please notify us")
			
		return 1

class EMStackSaveDialog(EMFileSaver):
	'''
	A class used for saving stacks, should work for stacks of all dimensions
	'''
	def __init__(self):
		EMFileSaver.__init__(self)

		self.__item_list = None # this will be a reference to the item list'
		self.file_name_used = "" #ultimately used to store the filename that was written to disk - if it occurred
		self.overwrite = False
		self.append = False
		
	def validate_save_argument(item_object):
		'''
		Verify that the function argument will work if supplied to the save function.
		This function could be updated if anyone ever makes it possible to save more objects
		@param item_object the item you want to test
		@return true or False
		'''
		
		#from emimagemx import EMMXDataCache
		fine = False
		if hasattr(item_object,"get_item_from_emsave"): # this is True EMMXDataCache
			fine = True
		elif isinstance(item_object,list):
		    if hasattr(item_object[0],"get_attr_dict"):
		    	fine = True

		return fine
	# mimic a static function
	validate_save_argument = staticmethod(validate_save_argument)
	
	def save(self,item_list):
		'''
		The main function
		@param item_list a list of items - will change to become more generic
		@raise RuntimeError if the the function argument is not acceptable
		'''
		fine = EMStackSaveDialog.validate_save_argument(item_list)
		if not fine:
			raise RuntimeError("item_list must be a list of EMData instances, a list of ListWidgetItems, or an EMDataListCache")

		self.__item_list = item_list
		from emselector import EMSelectorDialog
		selector = EMSelectorDialog(True,True)
		self.validator = EMSaveImageValidator(item_list)
		selector.set_validator(self.validator)
		file = selector.exec_()
		if file != "":
			self.__save_file(str(file))

	def __save_file(self,file):
		'''
		Called internally to save the file. Checks to make sure it's a valid file name
		and may run a dialog asking if the user wants to overwrite or append data.
		@param file the file name
		@return 0 if there was an error, 1 if there was not
		@exception NotImplementedError raised when the EMFileExistsDialog returns an unknown code
		'''
		msg = QtGui.QMessageBox()
			
		tmp_file_object = EMDummyTmpFileHandle(file)
		if self.validator.is_overwriting():
			tmp_file_object = EMTmpFileHandle(file)
		
		total_images = len(self.__item_list)

		out_file = tmp_file_object.get_tmp_file_name()
		progress = QtGui.QProgressDialog("Writing files", "abort", 0, 2*total_images,None)
		progress.show()
		tally = 0
		exc_list = None
		if hasattr(self.__item_list, "excluded_list"):
			# this is a stop gap - there is probably a better way to do it...
			exc_list = self.__item_list.excluded_list
		for i in range(total_images):
			if exc_list:
				try:
					(j for j in exc_list if j == i ).next() 
					# it's in the exc_list
					continue
				except: pass
			
			if isinstance(self.__item_list[i],EMData):
				d = self.__item_list[i]
			else:
				try:
					d = self.__item_list[i].get_data() # this will be case from the selector
				except:
					print "unknown situation" # contact David Woolford 
					# this might be redundant now
					#d = self.__item_list.get_item_from_emsave(i) # this will be the case from emimagemx
					#if d == None: continue # this will be the case if the image is shown as deleted in the emimagemx interface
			
			progress.setValue(tally)
			tally += 1
			
			
			# this is necessary to cause the progress bar to update on Mac, not sure about windows
			QtCore.QCoreApplication.instance().processEvents()
			
			try:
				d.write_image(out_file,-1)
			except:
				msg.setText("An exception occured while writing %s, please try again" %out_file)
				msg.exec_()
				tmp_file_object.remove_tmp_file()
				progress.close()
				return 1
			#else if d == None this is the equivalent of the particle being deleted, which makes sense for the EMDataListCache
				
				
			progress.setValue(tally)
			tally += 1
			# this is necessary to cause the progress bar to update on Mac, not sure about windows
			QtCore.QCoreApplication.instance().processEvents()
			if progress.wasCanceled():
				tmp_file_object.remove_tmp_file()
				progress.close()
				return 1
			
		
		tmp_file_object.finalize_renaming()			
		progress.close()
		self.file_name_used = tmp_file_object.get_final_file_name()		
		return 1
	
	def get_file_name(self):
		'''
		Get the file name that was used to write to disk
		Can be empty if file save was cancelled
		'''
		return self.file_name_used

class EMFileExistsDialog(QtGui.QDialog):
	'''
	Runs a dialog asking if the user wants to overwrite,append to, or cancel the operation.
	Appending may not be possible in which case nly the Overwrite and Cancel buttons are available.
	Return code is important - see exec_ documentation
	@code
	file_exists_dialog = EMFileExistsDialog(file,self.__item_list)
	return_code = file_exists_dialog.exec_()
	if return_code == 0:
		return 1 # operation was cancelled
	elif code == 1:
		overwrite = True # user chose overwrite
	elif code == 2: 
		append = True # user chose append
	else: raise
	@endcode
	
	'''
	def __init__(self,filename,item_list):
		'''
		@param filename the name of the file that is being overwritten
		@param item_list a list - the first object must have supply the get_attr_dict function, the keys "nx", "ny" and "nz therein
		If you want to use this function for a single EMData object then just put it in a list
		'''
		QtGui.QDialog.__init__(self,None)
		self.resize(480,320)
		self.setWindowIcon(QtGui.QIcon(get_image_directory() + "/eman.png"))
		self.setWindowTitle("File already exists")
		self.appendable_types = ["hed","img","spi","hdf"] #image types that can be stacks - TODO bdb ?
		
		append_enable = False 
		
		splt = filename.split(".")
		if splt[-1] in self.appendable_types or (len(filename) > 3 and filename[:4] == "bdb:"):
			(nx,ny,nz) = gimme_image_dimensions3D(filename)
			d = item_list[0].get_attr_dict()
			nz1 = d["nz"]
			ny1 = d["ny"]
			nx1 = d["nx"]
			# can only append if the dimensions are the same
			if nx == nx1 and ny == ny1 and nz == nz1:
				append_enable=True
		
		# some widgets
		vbl = QtGui.QVBoxLayout(self)
		hbl = QtGui.QHBoxLayout()
		overwrite = QtGui.QPushButton("Overwrite")
		cancel = QtGui.QPushButton("Cancel")
		
		# add Widgets to Layouts
		hbl.addWidget(cancel)
		if append_enable:
			append = QtGui.QPushButton("Append")
			hbl.addWidget(append)
		hbl.addWidget(overwrite)
		if (len(filename) > 3 and filename[:4] == "bdb:"):
			overwrite.setEnabled(False)
			overwrite.setToolTip("Overwriting bdb files is currently disabled.")
		
		# Text to alert the user
		hbl2 = QtGui.QHBoxLayout()
		text_edit = QtGui.QTextEdit("",self)
		text_edit.setReadOnly(True)
		text_edit.setWordWrapMode(QtGui.QTextOption.WordWrap)
		if (filename == ""):
			if append_enable:
				help = "The file already exists. You can choose to append to it, to overwrite it, or to cancel this operation."
			else:
				help = "The file already exists. You can choose to overwrite it, or to cancel this operation."			
		else:
			if append_enable:
				help = "The file %s already exists. You can choose to append to it, to overwrite it, or to cancel this operation." %filename
			else:
				help = "The file %s already exists. You can choose to to overwrite it, or to cancel this operation." %filename
		
		help += "\n\nBe careful about overwriting data as it may cause errors when using EMAN2 virtual stacks."
		
		text_edit.setText(help)
		hbl2.addWidget(text_edit,0)
		
		groupbox = QtGui.QGroupBox("Warning")
		groupbox.setLayout(hbl2)
		vbl.addWidget(groupbox)
		
		vbl.addLayout(hbl)
		
		if append_enable:
			QtCore.QObject.connect(append, QtCore.SIGNAL("clicked(bool)"), self.append_clicked)
		QtCore.QObject.connect(overwrite, QtCore.SIGNAL("clicked(bool)"), self.overwrite_clicked)
		QtCore.QObject.connect(cancel, QtCore.SIGNAL("clicked(bool)"), self.cancel_clicked)
		
		self.__result = 0
		
	def overwrite_clicked(self,bool):
		'''
		QtCore.SLOT
		calls accept to end execution
		'''
		self.__result = 1
		self.accept()
	
	def append_clicked(self,bool):
		'''
		QtCore.SLOT
		calls accept to end execution
		'''
		self.__result = 2
		self.accept()
		
	def cancel_clicked(self,bool):
		'''
		QtCore.SLOT
		calls accept to end execution
		'''
		self.__result = 0
		self.accept()
	
	def exec_(self):
		'''
		Wraps QtGui.QDialog.exec_ but returns a custom return value
		@return important integer code
		The return integer codes are as follows:
		0 - The user hit cancel
		1 - The user hit overwrite
		2 - The user hit append 
		'''
		QtGui.QDialog.exec_(self)
		return self.__result
	
class EMTmpFileHandle(object):
	'''
	A factory for getting a EMTmpFileHandle
	In the original design the type (e.g. hdf, img) of the file_name is critically important,
	seeing as it will be eventually renamed
	@param file_name the file_name used to deduce a temporary file name.
	@exception raised if the file_name does not exist on the file system
	@exception raised if the file_name has a unrecognized (or no) type
	'''
	def __new__(cls, file_name):
		if not file_exists(file_name): raise
		splt = file_name.split(".")
		if len(file_name) > 4 and file_name[:4] == "bdb:": pass
		elif EMUtil.get_image_ext_type(splt[-1]) == EMUtil.ImageType.IMAGE_UNKNOWN:
			raise
		
		if splt[-1] in ["img","hed"]: return EMImagicTmpFileHandle(file_name)
		elif file_name[:4] == "bdb:": return EMDBTmpFileHandle(file_name)
		else: return EMGeneralTmpFileHandle(file_name)
		# okay lots of tests, now return the right one

class EMTmpFileHandleBase:
	'''
	This class originally added to deal with issues that arise when users overwrite data on disk
	using file saving dialogs. You want to write to a temporary file, and when writing is finished,
	rename the temporary file to the original file (reasons not discussed).
	Hence this class will determine a temporary file name, supply it, and then rename
	it to the original file name once you have completed writing. Recovery functionality is also
	incorporated in that if there is a progress dialog and the user hits cancel, then you
	may call remove_tmp_file to clean up.
	The need for several classes arises from the existence of the Imagic and BDB formats
	'''
	
	def __init__(self):
		pass
	def get_tmp_file_name(self): raise NotImplementedException("Inheriting functions must supply this function")
	def get_final_file_name(self): raise NotImplementedException("Inheriting functions must supply this function")
	def finalize_renaming(self): raise NotImplementedException("Inheriting functions must supply this function")
	def remove_tmp_file(self): raise NotImplementedException("Inheriting functions must supply this function")
	
class EMDummyTmpFileHandle(EMTmpFileHandleBase):
	'''
	This class is for convenience in emselector.py - in the case where
	no tmp file is needed but the code wants to act in a uniform fashion, then
	I just use this object to mimic the function calls that otherwise happen
	'''
	def __init__(self,file_name):
		EMTmpFileHandleBase.__init__(self)
		self.__file_name = file_name
	
	def get_tmp_file_name(self): return self.__file_name
	def get_final_file_name(self): return self.__file_name
	def finalize_renaming(self):pass
	def remove_tmp_file(self): pass
	
	
class EMGeneralTmpFileHandle(EMTmpFileHandleBase):
	
	'''
	This tmp file handle works for most EM image types, such as mrc, hdf, spi, etc - e.g.
	all file types that have a single entry on disk.
	See EMTmpFileHandleBase for more information
	'''
	def __init__(self,file_name):
		'''
		@param file_name the name of the file that which is being overwritten
		@exception raised if the file_name does not exist on the file system
		@exception raised if the file_name has a unrecognized (or no) type
		'''
		if not os.path.exists(file_name): raise
		splt = file_name.split(".")
		if len(splt) < 2: raise
		if EMUtil.get_image_ext_type(splt[-1]) == EMUtil.ImageType.IMAGE_UNKNOWN:
			raise
		self.__file_ext = "."+splt[-1]
		
		EMTmpFileHandleBase.__init__(self)
		self.__file_name = base_name(file_name)
		self.__full_file_name = file_name
		self.__tmp_file_name = None
		self.__establish_tmp_file_name()

	def get_final_file_name(self): 
		'''
		Get the final file name of the file that is to be written to disk
		Returns the name regardless of whether or not it was actually written
		@return the the final file name
		'''
		return self.__full_file_name
	def get_tmp_file_name(self):
		'''
		Get the name of the temporary file 
		'''
		return self.__tmp_file_name

	def __gen_tmp_file_name(self):
		'''
		Called internally - attempts to make a unique file name using gm_time_string
		'''
		new_ending = gm_time_string()
		new_ending = new_ending.replace(" ","_") #just because 
		new_ending = new_ending.replace(":","_") #just because
		self.__tmp_file_name = self.__file_name + new_ending + self.__file_ext
	def __establish_tmp_file_name(self):
		'''
		establish self.__tmp_file_name
		'''
		self.__gen_tmp_file_name()
		while os.path.exists(self.__tmp_file_name):
			self.__gen_tmp_file_name()
			
	def remove_tmp_file(self):
		'''
		removes the temporary file
		'''
		if os.path.exists(self.__tmp_file_name):
			os.remove(self.__tmp_file_name)
		else:
			pass
				
	def finalize_renaming(self):
		'''
		Overwrite the original file, called when writing is done 
		'''
		os.rename(self.__tmp_file_name,self.__full_file_name)
	
class EMImagicTmpFileHandle(EMTmpFileHandleBase):
	'''
	This tmp file handle works for Imagic file types.
	See EMTmpFileHandleBase for more information
	'''
	def __init__(self,file_name):
		'''
		@param file_name the name of the file that which is being overwritten
		@exception if either the .hed or .img file corrsponding to the old_file_name does not exist on the file system
		@exception if the old_file_name is not an imagic file
		'''
		
		splt = file_name.split(".")
		if len(splt) < 2 or splt[-1] not in ["img","hed"]: raise
		
		self.__file_name_root = base_name(file_name)
		
		if not os.path.exists(self.__file_name_root+".img"): raise
		if not os.path.exists(self.__file_name_root+".hed"): raise
		
		EMTmpFileHandleBase.__init__(self)
		self.__file_name_root_hed = self.__file_name_root+".hed"
		self.__file_name_root_img = self.__file_name_root+".img"
		self.__establish_tmp_file_name()
		
	def get_tmp_file_name(self):
		'''
		Get the name of the temporary file - 
		@return the filename that ends in an "img" 
		'''
		return self.__tmp_file_name_img
	
	def __gen_tmp_file_names(self):
		'''
		Called internally - attempts to make a unique file name using gm_time_string
		'''
		new_ending = gm_time_string()
		new_ending = new_ending.replace(" ","_") #just because 
		new_ending = new_ending.replace(":","_") #just because
		self.__tmp_file_name_hed = self.__file_name_root + new_ending + ".hed"
		self.__tmp_file_name_img = self.__file_name_root + new_ending + ".img"
		
	def __establish_tmp_file_name(self):
		'''
		Called internally to establish temporary file name member variables
		'''
		self.__gen_tmp_file_names()
		while (os.path.exists(self.__tmp_file_name_hed) or os.path.exists(self.__tmp_file_name_img)):
			self.__gen_tmp_file_names()

	def remove_tmp_file(self):
		'''
		removes the temporary file
		'''
		if (os.path.exists(self.__tmp_file_name_hed) and os.path.exists(self.__tmp_file_name_img)):
			os.remove(self.__tmp_file_name_hed)
			os.remove(self.__tmp_file_name_img)
		else:
			pass
	
	def finalize_renaming(self):
		'''
		Overwrite the original file, called when writing is done. Overwrites the img and hed files
		'''
		os.rename(self.__tmp_file_name_hed,self.__file_name_root_hed)
		os.rename(self.__tmp_file_name_img,self.__file_name_root_img)

	def get_final_file_name(self): return self.__file_name_root_img


class EMDBTmpFileHandle(EMTmpFileHandleBase):
	'''
	This tmp file handle works for Imagic file types.
	See EMTmpFileHandleBase for more information
	'''
	def __init__(self,file_name):
		'''
		@param file_name the name of the file that which is being overwritten
		@exception RuntimeError raised of the file_name is not a valid bdb style database handle
		'''
		
		if not db_check_dict(file_name): raise RuntimeError("%s is not a valid database name" %file_name)
		self.__orig_db_name = file_name
		EMTmpFileHandleBase.__init__(self)
		self.__tmp_db_name = ""
		self.__establish_tmp_file_name()
	def get_tmp_file_name(self):
		'''
		Get the name of the temporary file - 
		@return the filename that ends in an "img" 
		'''
		return self.__tmp_db_name
	
	def __gen_tmp_file_name(self):
		'''
		Called internally - attempts to make a unique file name using gm_time_string
		'''
		new_ending = gm_time_string()
		new_ending = new_ending.replace(" ","_") #just because 
		new_ending = new_ending.replace(":","_") #just because
		self.__tmp_db_name = self.__orig_db_name + new_ending
		
	def __establish_tmp_file_name(self):
		'''
		Called internally to establish temporary file name member variables
		'''
		self.__gen_tmp_file_name()
		while (db_check_dict(self.__tmp_db_name)):
			self.__gen_tmp_file_name()

	def remove_tmp_file(self):
		'''
		removes the temporary file
		'''
		if (db_check_dict(self.__tmp_db_name)):
			db_remove_dict(self.__tmp_db_name)
		else:
			pass
	
	def finalize_renaming(self):
		'''
		Overwrite the original file, called when writing is done. Overwrites the img and hed files
		'''
		raise NotImplementedError("Woops waiting on an email")

	def get_final_file_name(self): return self.__orig_db_name
	
