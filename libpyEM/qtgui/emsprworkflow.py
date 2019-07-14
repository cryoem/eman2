#!/usr/bin/env python
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
#
# Author: David Woolford 11/10/08 (woolford@bcm.edu)
# Copyright (c) 2000-2008 Baylor College of Medicine
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

from past.utils import old_div
from builtins import range
from builtins import object
from .emform import EMFormWidget,EMParamTable,EMTableFormWidget
from .emdatastorage import ParamDef
from PyQt5 import QtGui, QtWidgets,QtCore
from PyQt5.QtCore import Qt
from EMAN2db import db_check_dict, db_open_dict,db_remove_dict,db_list_dicts,db_close_dict, e2getcwd
from EMAN2 import *
import os
import copy
from .emapplication import EMProgressDialog, get_application, EMErrorMessageDisplay, error
from e2ctf import pspec_and_ctf_fit,GUIctf,write_e2ctf_output,get_gui_arg_img_sets,init_sfcurve
import subprocess
import weakref
import time
from .emsave import save_data
from .emimagemx import EMDataListCache
import traceback

USING_RELATIVE_DIRS = True # used by database infrastructure for recording file names

def workflow_path(path="",dir=None):
	'''
	Makes paths relative if possible.
	'''
	
	if USING_RELATIVE_DIRS:
		if dir==None: dir = e2getcwd()
		from .emselector import folderize
		dir = folderize(dir)
		name = path
		if dir in name:
			if name[:4].lower()=="bdb:" :
				name="bdb:"+name[len(dir)+4:]
			else: name = name[len(dir):]
		
		return name
	
	return path

spr_raw_data_dict = "global.spr_raw_data_dict"
spr_filt_ptcls_map = "global.spr_filt_ptcls_map"# this is for recovery and back compatibility
spr_ptcls_dict = "global.spr_ptcls_dict"
spr_rfca_dict = "global.spr_rfca_dict"
spr_stacks_map = "global.spr_stacks_map"
spr_sets_dict = "global.spr_sets_dict"
spr_init_models_dict = "global.spr_init_models_dictx" #my DB has been corrupted!!!!!
spr_freealign_models_dict = "global.spr_freealign_models_dict"
spr_freealign_dirs_dict = "global.spr_freealign_dirs_dict"

#These are actually used
tpr_subtomo_stacks = "global.tpr_subtomo_stacks"
tpr_raw_data_dict = "global.tpr_raw_data_dict"
tpr_subtomo_ref = "global.tpr_subtomo_ref"

#Currently these are garbage!!
tpr_filt_ptcls_map = "global.tpr_filt_ptcls_map"
tpr_ptcls_dict = "global.tpr_ptcls_dict"
tpr_probes_dict = "global.tpr_probes_dict"
tpr_ptcl_ave_dict = "global.tpr_ptcl_ave_dict"
tpr_ptcls_ali_dict = "global.tpr_ptcls_ali_dict"

class EmptyObject(object):
	'''
	This just because I need an object I can assign attributes to, and object() doesn't seem to work
	'''
	def __init__(self):
		pass

class WorkFlowTask(object):
	display_file = QtCore.pyqtSignal()
	task_idle = QtCore.pyqtSignal()
	process_started = QtCore.pyqtSignal()

	def __init__(self):
		self.window_title = "Set me please" # inheriting classes should set this
		self.preferred_size = (480,640) # inheriting classes can change this if they choose
		self.form_db_name = None # specify this to make use of automated parameter storage (see write_db_entries(self,...) ) - don't forget the "bdb:"
		self.project_db_entries = ["global.num_cpus","global.apix","global.microscope_voltage","global.microscope_cs","global.memory_available","global.particle_mass"] # used to write entries to a specific db
		self.core_object = QtCore.QObject()
	
	def emitter(self):
		return self.core_object
	
	
	def __del__(self):
		self.core_object.deleteLater()
		
	def run_form(self):
		self.form = EMFormWidget(self.get_params())
		self.form.resize(*self.preferred_size)
		self.form.setWindowTitle(self.window_title)
		get_application().show_specific(self.form)
		self.form.emform_ok.connect(self.on_form_ok)
		self.form.emform_cancel.connect(self.on_form_cancel)
		self.form.emform_close.connect(self.on_form_close)
		self.form.display_file.connect(self.on_display_file)
	
	def get_params(self): raise NotImplementedError
	
	def on_display_file(self,filename):
		self.display_file.emit(filename)
		
	def on_form_ok(self,params):
		for k,v in list(params.items()):
			self.write_db_entry(k,v)
		
		self.disconnect_form()
		self.form.close()
		self.form = None
	
		self.task_idle.emit()
		
	def on_form_cancel(self):
		self.disconnect_form()
		self.form.close()
		self.form = None
		self.task_idle.emit()
	
	def disconnect_form(self):
		self.form.emform_ok.disconnect(self.on_form_ok)
		self.form.emform_cancel.disconnect(self.on_form_cancel)
		self.form.emform_close.disconnect(self.on_form_close)
		self.form.display_file.disconnect(self.on_display_file)
	
	
	def emit(self,*args,**kargs):
		self.core_object.emit(*args,**kargs)
		
	def on_form_close(self):
		self.disconnect_form()
		#self.task_idle.emit()

	def close(self):
		if self.form != None: 
			self.form.close()
		self.task_idle.emit()
        
	def closeEvent(self,event):
		self.close()
		
	def write_db_entries(self,dictionary):
		'''
		Write the dictionary key/entries into the database using self.form_db_name
		Writes all keys except for "blurb" - note the the "blurb" key is mostly used in the context
		of these forms to display helpful information to the user - it doesn't need to be stored in the
		database 
		'''
		if self.form_db_name != None: db = db_open_dict(self.form_db_name)
		else: db = None
		
		project_db = db_open_dict("bdb:project")
		for k,v in list(dictionary.items()):
			if k == "blurb": continue
			
			if k in self.project_db_entries: project_db[k] = v
			else:
				if db != None: db[k] = v
		
		#if self.form_db_name != None: db_close_dict(self.form_db_name)
		#db_close_dict("bdb:project")
		
	def get_default_filenames_from_form_db(self,key="filenames"):
		'''
		Opens the self.form_db_name database and retrieves the filenames entry
		Returns None if self.form_db_name is None
		Returns an empty list if the "filenames" entry doesn't exist in the existing database
		'''
		default_selections = None
		if self.form_db_name != None:
			db = db_open_dict(self.form_db_name)
			default_selections = db.get(key,dfl=[])
			#db_close_dict(self.form_db_name)
			
		return default_selections
	
	def write_db_entry(self,key,value):
		'''
		This function is becoming deprecated, used write_db_entries instead
		'''
		db = db_open_dict("bdb:project")
		if len(key) > 5 and key[:6] == "global":
			db[key] = value
		else:
			pass
		
		#db_close_dict("bdb:project")
	def get_wd(self):
		'''
		Get the working directory, originally introduced to provide a centralized mechanism for accessing the working directory,
		specifically for the purpose of spawning processes. Could be used more generally, however.
		'''
		return e2getcwd()

	def spawn_task(self,program,options,string_args,bool_args,additional_args=[],temp_file_name="e2workflow_tmp.txt"):
		'''
		splits the task over the available processors
		example-
		program="e2ctf.py"
		options is an object with the filenames, all string args and all bool_args as attributes 
		string_args=["
		bool_args=["
		additional_args=["--auto_db,--autofit"]
		temp_file_name = "etctf_auto_tmp.txt"
		'''
		project_db = db_open_dict("bdb:project")
		ncpu = project_db.get("global.num_cpus",dfl=num_cpus())
		cf = old_div(float(len(options.filenames)),float(ncpu)) # common factor
		
		files = []
		for n in range(ncpu):
			files.append([])
		
		if len(temp_file_name) > 3 and temp_file_name[-4] == ".":
			temp_fname_root = temp_file_name[:-4]
			temp_fname_end = temp_file_name[-4:]
		else:
			temp_fname_root = temp_file_name
			temp_fname_end = ""
		# distribute the names into bins
		for i,f in enumerate(options.filenames):
			idx = i % ncpu
			files[idx].append(f) 
		
		for n in range(ncpu):
			#print "n"
			filenames = files[n]
			if len(filenames) == 0: continue # maybe there are more CPUS then filenames
								
			args = [e2getinstalldir()+"/bin/"+program]
	
			for name in filenames:
				args.append(name)
			
			for string in string_args:
				args.append("--"+string+"="+str(getattr(options,string)))

			# okay the user can't currently change these, but in future the option might be there
			for string in bool_args:
				# these are all booleans so the following works:
				if getattr(options,string):
					args.append("--"+string)
					
			for arg in additional_args:
				args.append(arg)
#			print "command is ",program
			for i in args: print(i, end=' ')
			print()
			
			#print args
			#fname = temp_fname_root +"_"+str(n)+temp_fname_end
			#file = open(fname,"w+")
			
			cmdstr = ' '.join(args)
#			process = subprocess.Popen(args_adjusted,stdout=file,stderr=subprocess.STDOUT)
			process = subprocess.Popen(cmdstr, shell=True)
			print("started process",process.pid)
			self.process_started.emit(process.pid)
			
		#db_close_dict("bdb:project")
	
	def spawn_single_task(self,program,options,string_args,bool_args,additional_args=[],temp_file_name="e2workflow_tmp.txt"):
		'''
		runs a single job
		example-
		program="e2ctf.py"
		options is an object with the filenames, all string args and all bool_args as attributes 
		string_args=["
		bool_args=["
		additional_args=["--auto_db,--autofit"]
		temp_file_name = "etctf_auto_tmp.txt"
		'''
		project_db = db_open_dict("bdb:project")	
								
		#args = [program]
		args = [e2getinstalldir()+"/bin/"+program]
		if hasattr(options,"filenames"):
			for name in options.filenames:
				args.append(name)

		for string in string_args:
			args.append("--"+string+"="+str(getattr(options,string)))

		# okay the user can't currently change these, but in future the option might be there
		for string in bool_args:
			# these are all booleans so the following works:
			if getattr(options,string):
				args.append("--"+string)
				
		for arg in additional_args:
			args.append(arg)
		
		# this prints the full command
		for i in args: print(i, end=' ')
#		print ""
#		
		#print args
		#file = open(temp_file_name,"w+")

		cmdstr = ' '.join(args)
#		process = subprocess.Popen(args_adjusted,stdout=file,stderr=subprocess.STDOUT)
		process = subprocess.Popen(cmdstr, shell=True)
		print("started process",process.pid)
		self.process_started.emit(process.pid)
		
		#db_close_dict("bdb:project")
		
	def run_select_files_msg(self):
		'''
		Runs a QMessageBox asking for the user to select files for processing
		'''
		msg = QtWidgets.QMessageBox()
		msg.setWindowTitle("Almost")
		msg.setText("Please select files for processing")
		msg.exec_()

	def show_error_message(self,error_message):
		'''
		error_message is a list of error messages
		'''
		EMErrorMessageDisplay.run(error_message)
		
	def get_latest_r2d_classes(self):
		dirs = get_numbered_directories("r2d_")
		# all right everything left in dirs is "r2d_??" where the ?? is castable to an int, so we should be safe now
		class_files = []
		class_dims = []
		class_ptcls = []
		for dir in dirs:
			classes_db = None
			# check for 00 to 09 but 00 is replaced with "init"
			db_first_part = "bdb:"+dir+"#classes_"
			cont = True
			for i in range(0,10):
				for j in range(0,10):
					if i == 0 and j == 0:
						db_name = db_first_part+"init"
					else:
						db_name = db_first_part+str(i)+str(j)
						
					if db_check_dict(db_name):
						classes_db = db_name
					else:
						if i != 0 or j != 0:
							cont = False
							break
						#else just check for 01 in case the user has specified the --initial argument
				if not cont:
					break
				
			if classes_db != None:
				class_files.append(classes_db)
				cl_db = db_open_dict(classes_db,ro=True)
				if "maxrec" in cl_db:
					class_ptcls.append(cl_db["maxrec"]+1)
					hdr = cl_db.get_header(0)
					class_dims.append(str(hdr["nx"])+'x'+str(hdr["ny"])+'x'+str(hdr["nz"]))
				else:
					class_ptcls.append("")
					class_dims.append("")
					
				#db_close_dict(classes_db)
					
		return class_files,class_ptcls,class_dims
	
	def check_sym(self,params,options):
		error_message = []
		for key in params:
			if key[:7] == "symname":	
				if params[key] in ["c","d","h"]:
					n = params["symnumber"]
					fail = False
					if len(n) == 0: fail = True
					try: int(n)
					except: fail = True
			
					if not fail:
						if int(n) < 1:
							fail = True
					
					if fail:
						error_message.append("The symmetry number must be specified for c,d, and h.")
					else:
						options.sym=params["symname"]+n
				elif len(params["symnumber"]) != 0:
					error_message.append("There is something entered in the symmetry number box but you have not specified c, d or h symmetry.")
				else:
					options.sym = params["symname"]	
			
		return error_message
	
	def get_cmps_list(self):
		return list(dump_cmps_list().keys())
	
	def get_aligners_list(self):
		return list(dump_aligners_list().keys())
	
	def get_projectors_list(self):
		return list(dump_projectors_list().keys())
	
	def get_orientgens_list(self):
		return list(dump_orientgens_list().keys())
		
	def get_averagers_list(self):
		return list(dump_averagers_list().keys())
		
#		cmps.append("None") I think this is necessary
		
class EMProjectDataDict(object):
	''' This class encapsulate the common routines used to get, add to and remove data dictionaries from the database.
	These data dictionaries are used for storing almost all types of data in the workflow, and are persistently located
	on disk using Berkeley DBs. The keys of these dictionaries are unique identifiers, such as the full name of a
	file on disk. The values are dictionaries that enumerate any filtered data that are associated with this particular key.
	It became necessary for us to make it so that one "these filtered data entries" is actually a named version of the
	original data - it's key will be "Original Data". Here's what this dict looks like:
	
	{"bdb:particles#1002_ptcls":{"Original Data":"bdb:particles#1002_ptcls", "Phase flipped":"bdb:particles#1002_ptcls_phase"}}
	'''
	original_data = "Original Data"
	recovery_list_to_dict = None
	recovery_dict_to_dict = None
	def __init__(self,data_dict_name="global.spr_raw_ptcls_map",db_name="bdb:project"):
		'''
		@param data_dict_name the key to a database entry which holds a dictionary, the keys of which are strings the values of which are dictionaries. See http://blake.bcm.edu/emanwiki/Eman2AppMetadata
		@param db_name the name of the database opened by the called to db_open_dict. Should probably be bdb:project
		@param recovery_list_name in June 2009 we converted list style storage to the more versatile dictionary style storage. This is the old list name. For back compatibility
		@param recovery_dict_name possibly a recovery mechanism if we ever change the name of the data_dict_name. Not currently used
		
		Terminology : I will call big dictionary the 'data dictionary' and the value of this which is a dictionary the 'filter dictionary' 
		'''
		self.db_name = db_name # the name of the database
		self.data_dict_name = data_dict_name # the key of the database dictionary
		self.data_dict = None # created the first time it's need
		
		if EMProjectDataDict.recovery_list_to_dict == None:
			self.__init_recovery_list_to_dict()
		try:
			self.recovery_list_name =  EMProjectDataDict.recovery_list_to_dict[self.data_dict_name]
		except:
			self.recovery_list_name = None
			
		if EMProjectDataDict.recovery_dict_to_dict == None:
			self.__init_recovery_dict_to_dict()
		try:
			self.recovery_dict_name =  EMProjectDataDict.recovery_dict_to_dict[self.data_dict_name]
		except:
			self.recovery_dict_name = None
		
		
	def __init_recovery_list_to_dict(self):
		'''
		Initializes EMProjectDataDict.recovery_list_to_dict
		'''
		d = {}
		d[spr_ptcls_dict] = "global.spr_ptcls"
		d[spr_raw_data_dict] = "global.spr_raw_file_names"
		d[spr_rfca_dict] = "global.spr_ref_free_class_aves"
		d[spr_init_models_dict] = "global.spr_init_models"
		d[spr_sets_dict] = "global.spr_stacks"
		d[tpr_raw_data_dict] = "global.tpr_tomograms"
		d[tpr_ptcls_dict] = "global.tpr_ptcls"
		d[tpr_probes_dict] = "global.tpr_probes"
		d[tpr_ptcl_ave_dict] = "global.tpr_ptcl_averages"
		EMProjectDataDict.recovery_list_to_dict = d
		
	def __init_recovery_dict_to_dict(self):
		'''
		Initializes EMProjectDataDict.recovery_dict_to_dict
		'''
		d = {}
		d[spr_ptcls_dict] = "global.spr_filt_ptcls_map"
		d[spr_sets_dict] = "global.spr_stacks_map"
		d[tpr_ptcls_dict] = "global.tpr_filt_ptcls_map"
		EMProjectDataDict.recovery_dict_to_dict = d
		
	def keys(self):
		'''
		@return the keys of the data dictionary
		'''
		dict = self.get_data_dict()
		return list(dict.keys())
		
	def get_names(self,filt=original_data):
		'''
		@param filt the key to retrieve from the filter dictionary
		@return a list of names that have entries corresponding to this filter entry
		'''
		dict = self.get_data_dict()
		ret = []
		for tag,dict2 in list(dict.items()):
			for key,name in list(dict2.items()):
				if key == filt:
					ret.append(name)
		return ret
	
	def add_names(self,list_of_names,filt=original_data,use_file_tag=False):
		'''
		@param list_of_names a list of names - should be full path file names
		@param filt the key to use when inserted the names in the filter dictionary
		@param use_file_tag - if true the keys in the main dictionary are file tags, not full names - this is useful for particles, for example
		The names in the list of names are used as keys in the data dictionary, and an 
		entry is added to the filter dictionary using the filt parameters (will be the key)
		and the name (will be the value)
		'''
		dict = self.get_data_dict()

		for name in list_of_names:
			
			n = workflow_path(name)

			tmp = {}
			tmp[filt] = n
			if use_file_tag: dict[base_name(n)] = tmp
			else: dict[n] = tmp
			
			
		self.__write_db(dict)
		
	def update(self,more_data_dict):
		'''
		@param more_data_dict a dictionary that will be added to the data dictionary using dict.update
		'''
		dict = self.get_data_dict()
		dict.update(more_data_dict)
		self.__write_db(dict)
		
		
	def add_name(self,name):
		'''Convenience function see add_names'''
		self.add_names([name])
			
	def remove_names(self,list_of_names):
		'''
		@param list_of_names a list of names to be removed from this data dictionary
		The names should be keys in the data dictionary. If not you get an error message
		'''
		dict = self.get_data_dict()
		err = []
		for name in list_of_names:
			try:
				dict.pop(name)
			except:
				err.append("%s is not a key in the %s dictionary" %(name,self.data_dict_name))
		
		if len(err) > 0:
			error(err)
		
		self.__write_db(dict)
	
	def remove_name(self,name): 
		'''Convenience function see remove_names'''
		self.remove_names([name])
		
	def get_data_dict(self):
		'''
		Called internally to get the self.data_dict attribute - creates it the first time the function
		@return the database dictionary itself
		'''
		if self.data_dict == None:
			
			if self.recovery_dict_name != None:
				# also uncomment this line in about a year - i.e. around June 2010
				self.__recover_dict_to_dict_from_old_name()
			elif self.recovery_list_name != None:
#				# uncomment this line in about a year - i.e. around June 2010
				self.__recover_list_to_dict_from_old_name()
				
			# also uncomment this line in about a year - i.e. around June 2010
			self.__recover_missing_key_dict()
			
			# remove the dictionary if any data was lost - do not delete
			self.__remove_dict_if_files_lost()
			
			# this recovery function removed on August 14th
#			if self.data_dict_name == spr_ptcls_dict:
#				self.__convert_ptcls_dict_to_tags()
			
#			print self.db_name
			project_db = db_open_dict(self.db_name)
#			print project_db
			self.data_dict = project_db.get(self.data_dict_name,dfl={})
		
		return self.data_dict
	
	def __write_db(self,dict):
		'''
		@param dict what can be considered the data dictionary corresponding to the state of this class object
		Write the data to disk
		'''
		project_db = db_open_dict(self.db_name)
		project_db[self.data_dict_name] = dict
		
#	def __convert_ptcls_dict_to_tags(self):
#		'''
#		Converts the ptcls dict so that its keys are file tags, not full names
#		This was deemed necessary when we ran into problems using the full names for particle files,
#		specifically at the point of building sets and removing bad particles
#		This function is only called if the self.data_dict_name is spr_ptcls_dict (in self.get_data_dict)
#		'''
#		
#		project_db = db_open_dict(self.db_name)
#		data_dict = project_db.get(spr_ptcls_dict,dfl={})
#		acted = False
#		for key, value in data_dict.items():
#			tag = base_name(key)
#			print tag,key
#			if tag != key:
#				acted = True
#				data_dict.pop(key)
#				data_dict[tag] = value
#		
#		if acted:
#			error("Converted global.spr_ptcls_dict to using file tags","Warning")
#			project_db[spr_ptcls_dict] = data_dict
	
	def __remove_dict_if_files_lost(self):
		'''
		Looks at all of the data in the filter dictionary, if any of that data does not exist the
		entry is automatically removed from the data dictionary
		This function can probably be removed a year from now, somewhere around June 2010
		'''
		project_db = db_open_dict(self.db_name)
		dict = project_db.get(self.data_dict_name,dfl={})
		rem = []
		for name, map in list(dict.items()):			
			for key,image_name in list(map.items()):
				try:
					if not file_exists(image_name) and not os.path.isdir(image_name[:image_name.rindex("/")]):
						map.pop(key)
						rem.append("%s data no longer exists for %s. Automatically removing from project" %(key,image_name))
				except: pass
			if len(map) == 0:
				dict.pop(name)
				rem.append("%s no longer exists. Automatically removing from project" %(name))
		
		if len(rem) != 0:
			EMErrorMessageDisplay.run(rem,"Warning")
			project_db[self.data_dict_name] = dict
				

	def __recover_dict_to_dict_from_old_name(self):
		'''
		Recovers dicts name in old style format - I had to do this to correct a naming inconsistency
		This function can probably be removed a year from now, somewhere around June 2010
		'''
		project_db = db_open_dict(self.db_name)
		acted = False
		if self.data_dict_name not in project_db or project_db[self.data_dict_name]==None or len(project_db[self.data_dict_name]) == 0:
			if self.recovery_dict_name in project_db:
				dict = project_db[self.recovery_dict_name]
				project_db[self.data_dict_name] = dict
	#			project_db[self.recovery_dict_name] = None # the only problem is if someone started using this name for something else # not necessary
				acted = True
			
		if acted:
			EMErrorMessageDisplay.run("Automatically renamed old database %s to %s " %(self.recovery_dict_name,self.data_dict_name),"Warning")
			
	def __recover_missing_key_dict(self):
		'''
		Detects if the EMProjectDataDict.original_data is present, if not inserts it, but only inserts it
		if the file exists. This will preserve any sets, for example, that were created from flat files.
		A highly specialized function that works a tight rope.
		It's for back compatibility. It should be fine.
		'''
		project_db = db_open_dict(self.db_name)
		dict = project_db.get(self.data_dict_name,dfl={})
		update = False
		for tag,dict2 in list(dict.items()):
			if EMProjectDataDict.original_data not in dict2:
				if file_exists(tag):
					dict2[EMProjectDataDict.original_data] = tag
					update = True
		if update:
			EMErrorMessageDisplay.run("Automatically converted an old dict for the %s database" %self.data_dict_name,"Warning")
			project_db[self.data_dict_name] = dict
		
			
	def __recover_list_to_dict_from_old_name(self):
		'''
		Recovers an old-style list if it existed.
		This function can probably be removed a year from now, somewhere around June 2010
		'''
		project_db = db_open_dict("bdb:project")
		if self.recovery_list_name == None or self.recovery_list_name not in project_db: return
		
		project_data  = project_db.get(self.data_dict_name,dfl={})
		old_project_data  = project_db.get(self.recovery_list_name,dfl=[])
		if len(old_project_data) > 0:
			EMErrorMessageDisplay.run("Recovering %s database (for back compatibility)" %self.recovery_list_name, "Warning")
			for name in old_project_data:
				if name not in project_data:
					project_data[name] = {}
			project_db[self.data_dict_name] = project_data
			
		project_db[ self.recovery_list_name] = None

	
class EMRawDataReportTask(WorkFlowTask):	
	'''This form displays raw data that are associated with the project. You browse to add raw data, or right click and choose Add.''' 
	documentation_string = "This page shows raw micrographs/ccd frames currently associated with the project. It is possible to add additional images directly on this panel, which will \
leave them in-place and not copy them into the project database. This will limit some later operations and leave the project with less metadata at the end, but will save disk space. \
Note that the data cannot be filtered unless it is imported."
	warning_string = "\n\n\nNOTE: There are no images currently associated with the project. Please associate or import images"
	def __init__(self):
		WorkFlowTask.__init__(self)
		self.window_title = "Micrographs In Project"
		self.project_data_at_init = None # stores the known project files when the form is created and shown - that way if cancel is hit we can restore the original parameters
		self.project_list = spr_raw_data_dict
		
	def get_image_dimensions(file_name):
		'''
		A static function for getting the dimensions of a file as a string
		'''
		try:
			
			nx,ny,nz = gimme_image_dimensions3D(file_name)
			return "%ix%ix%i" %(nx,ny,nz)
		except:
			return "Error"
	
	get_image_dimensions = staticmethod(get_image_dimensions)

	def get_raw_data_table(self):
		'''
		Gets an EM2DFileTable - this is type of class that the emform knows how to handle 
		'''
		
		data_dict = EMProjectDataDict(self.project_list)
		self.project_data_at_init = data_dict.get_data_dict() # so if the user hits cancel this can be reset
		print(self.project_data_at_init)
		project_names = list(data_dict.keys())
		
		from .emform import EM2DFileTable,EMFileTable
		table = EM2DFileTable(project_names,desc_short="Raw Data Files",desc_long="")
		context_menu_data = EMRawDataReportTask.ProjectListContextMenu(self.project_list)
		table.add_context_menu_data(context_menu_data)
		table.add_button_data(EMRawDataReportTask.ProjectAddRawDataButton(table,context_menu_data))
	
		#p.append(pdims) # don't think this is really necessary
		return table,len(project_names)
	
	def get_raw_data_table_custom(self):
		'''
		Calls get_raw_data_table and then adds the Dimensions column
		'''
		table,n = self.get_raw_data_table()
		from .emform import EMFileTable
		table.add_column_data(EMFileTable.EMColumnData("Dimensions",EMRawDataReportTask.get_image_dimensions,"The dimensions of the file on disk"))
		return table,n
	
	class ProjectAddRawDataButton(object):
		def __init__(self,table_widget,context_menu_data):
			self.table_widget = weakref.ref(table_widget)
			self.context_menu_data = context_menu_data
			self.name = "Browse To Add"
			
		def function(self,bool):
			self.context_menu_data.context_menu["Add"]([],self.table_widget())
			

	class ProjectListContextMenu(object):
		def __init__(self,project_list=spr_raw_data_dict,remove_only=False,using_file_tags=False):
			self.project_list = project_list
			self.validator = AddFilesToProjectValidator(self.project_list)
			self.context_menu = {}
			self.context_menu["Remove"] = EMRawDataReportTask.ProjectListContextMenu.RemoveFilesFromProject(self.project_list,using_file_tags)
			if not remove_only: self.context_menu["Add"] = EMRawDataReportTask.ProjectListContextMenu.AddFilesToProjectViaContext(self.project_list)
		
		def items(self):
			return list(self.context_menu.items())
		
		
		class RemoveFilesFromProject(object):
			def __init__(self,project_list,using_file_tags=False):
				self.project_list = project_list
				self.using_file_tags = using_file_tags
			def __call__(self,names,table_widget):
				if len(names) == 0: return # nothing happened
			
				from .emform import get_table_items_in_column
				entries = get_table_items_in_column(table_widget,0)
				text_entries = [table_widget.convert_text(str(i.text())) for i in entries]
				
				data_dict = EMProjectDataDict(self.project_list)
				project_names = list(data_dict.keys())
				
				full_names = [table_widget.convert_text(name) for name in names]
				db_full_names = [table_widget.convert_text(name) for name in names]
				if self.using_file_tags:
					db_full_names = [base_name(table_widget.convert_text(name)) for name in names]
					
				#print names
				#print project_names
				#print full_names
				#print db_full_names
				
				for name in db_full_names:
					if name not in project_names:
						# this shouldn't happen, but luckily it was here because it alerted me to an error once
						EMErrorMessageDisplay.run(["%s is not in the project list" %name])
						return
				
				indices = [ text_entries.index(name) for name in full_names]
				indices.sort()
				indices.reverse()
				for idx in indices:
					table_widget.removeRow(idx)
					
				data_dict.remove_names(db_full_names)
				
		class AddFilesToProject(object):
			def __init__(self,project_list):
				self.project_list = project_list
				
			def __call__(self,list_of_names,table_widget):
		
				data_dict = EMProjectDataDict(self.project_list)
				project_names = list(data_dict.keys())
				
				for name in list_of_names:
					if not file_exists(name) and not os.path.isdir(name[:name.rfind("/")]):
						EMErrorMessageDisplay.run(["%s does not exists" %name])
						return
				
				for name in list_of_names:
					if name in project_names:
						EMErrorMessageDisplay.run(["%s is already in the project" %name])
						return
						
				# if we make it here we're good
				# first add entries to the table
				table_widget.add_entries(list_of_names)
				data_dict.add_names(list_of_names)
		
		class AddFilesToProjectViaContext(object):
			task_idle = QtCore.pyqtSignal()

			def __init__(self,project_list):
				self.project_list = project_list
				self.validator = AddFilesToProjectValidator(self.project_list)
				
			def __call__(self,list_of_names,table_widget):
			
		#def add_files_from_context_menu(self,list_of_names,table_widget):
				from .emselector import EMSelectorDialog
				selector = EMSelectorDialog(save_as_mode=False)
				
				selector.set_selection_text("Selection(s)")
				selector.set_validator(self.validator)
				files = selector.exec_()
				selector.close()

				if files != "":
					if isinstance(files,str): files = [files]
					
					from .emform import get_table_items_in_column
					entries = get_table_items_in_column(table_widget,0)
					strentires = [str(i.text())[4:] for i in entries]
					strfiles = [i[5+len(os.getcwd()):] for i in files]
					error_messages = []
					for idx,tag in enumerate(strfiles):
						if tag in strentires:
							error_messages.append("%s is already listed" %files[idx])
					
				
					if len(error_messages) > 0:
						EMErrorMessageDisplay.run(error_messages)
						return
				a = EMRawDataReportTask.ProjectListContextMenu.AddFilesToProject(self.project_list)
				a(files,table_widget)
				#self.add_files_to_project(files,table_widget)
			
	def get_params(self):
		
		project_db = db_open_dict("bdb:project")
		
		params = []
		
		p,n = self.get_raw_data_table_custom()

		params.append(ParamDef(name="blurb",vartype="text",desc_short="Files",desc_long="",property=None,defaultunits=self.__doc__,choices=None))
		params.append(p)
			
		return params
	
	def on_form_cancel(self):
		self.recover_original_raw_data_list()
		
		self.form.close()
		self.form = None
		self.task_idle.emit()
	
	def recover_original_raw_data_list(self):
		'''
		Called if the user hits cancel - if they removed some files or added files the changes
		are not saved unless the user hits ok
		'''
		project_db = db_open_dict("bdb:project")
		project_db[self.project_list] = self.project_data_at_init
	
	
	def on_form_ok(self,params):
		self.form.close()
		self.form = None
		

		
class AddFilesToProjectValidator(object):
	def __init__(self,project_list=spr_raw_data_dict):
		self.project_list = project_list
	def validate_file_name(self,list_of_names):
		'''
		a validator for the select module
		@exception RuntimeError thrown if list_of_names is not a list
		@return 0 if something went wrong , 1 if it's okay to call save now
		'''
		if not isinstance(list_of_names,list):
			if isinstance(list_of_names,str): list_of_names = [list_of_names]
			else: raise RuntimeError("Files needs to be a list")
		
		data_dict = EMProjectDataDict(self.project_list)
		project_names = list(data_dict.keys())
		
		for name in list_of_names:
			if not file_exists(name) and not os.path.isdir(name[:name.rfind("/")]): #See if its a dir, which is ok:
				EMErrorMessageDisplay.run(["%s does not exists" %name])
				return 0
		
		for name in list_of_names:
			if name in project_names:
				EMErrorMessageDisplay.run(["%s is already in the project" %name])
				return 0
			
		return 1

class ParticleWorkFlowTask(WorkFlowTask):
	'''
	Encapsulates some functionality  common to the particle based work flow tasks
	Such tasks should inherit from this class, not the the WorkFlowTask
	'''
	def __init__(self):
		WorkFlowTask.__init__(self)

	def get_particle_selection_table(self,ptcl_list,table=None,single_selection=False,enable_ctf=True):
		'''
		
		'''
		from .emform import EM2DStackTable,EMFileTable,float_lt,int_lt
		if table==None: # you can hand in your own table (it's done in E2ParticleExamineTask)
			table = EM2DStackTable(ptcl_list,desc_short="Particles",desc_long="",single_selection=single_selection)
		if len(ptcl_list) != 0:
			a = EMData()
			a.read_image(ptcl_list[0],0,True)
			d = a.get_attr_dict()
			if enable_ctf and "ctf" in d:
				self.column_data = CTFColumns()
				table.add_column_data(EMFileTable.EMColumnData("Defocus",self.column_data.get_defocus,"The estimated defocus",float_lt))
				table.add_column_data(EMFileTable.EMColumnData("B Factor",self.column_data.get_bfactor,"The estimated B factor, note this is ~4x greater than in EMAN1",float_lt))
				table.add_column_data(EMFileTable.EMColumnData("SNR",self.column_data.get_snr,"The averaged SNR",float_lt))
				table.add_column_data(EMFileTable.EMColumnData("Quality",self.column_data.get_quality,"The quality of the fit as judged by e2ctf",int_lt))
				table.add_column_data(EMFileTable.EMColumnData("Sampling",self.column_data.get_sampling,"The amount of sampling used for generating CTF parameters",int_lt))
		else:
			context_menu_data = ParticleWorkFlowTask.DataContextMenu()
			table.add_context_menu_data(context_menu_data)
			table.add_button_data(ParticleWorkFlowTask.AddDataButton(table,context_menu_data))
		#table.insert_column_data(1,EMFileTable.EMColumnData("Particles On Disk",EMParticleReportTask.get_num_ptcls,"Particles currently stored on disk that are associated with this image"))
		#table.insert_column_data(2,EMFileTable.EMColumnData("Particle Dims",EMParticleReportTask.get_particle_dims,"The dimensions of the particles that are stored on disk"))
	
				
		table.add_column_data(EMFileTable.EMColumnData("Particles On Disk",EMParticleReportTask.get_num_ptcls,"Particles currently stored on disk that are associated with this image",int_lt))
		table.add_column_data(EMFileTable.EMColumnData("Particle Dims",EMParticleReportTask.get_particle_dims,"The dimensions of the particles that are stored on disk"))

		return table, len(ptcl_list)

	
	class AddDataButton(object):
		def __init__(self,table_widget,context_menu_data):
			self.table_widget = weakref.ref(table_widget)
			self.context_menu_data = context_menu_data
			self.name = "Browse To Add"
			
		def function(self,bool):
			self.context_menu_data.context_menu["Add"]([],self.table_widget())
			

	class DataContextMenu(object):
		def __init__(self,validator=None):
			
			self.validator = validator
			self.context_menu = {}
			self.context_menu["Remove"] = ParticleWorkFlowTask.RemoveDataFromTable()
			self.context_menu["Add"] = ParticleWorkFlowTask.AddDataToTable(validator)
		
		def items(self):
			return list(self.context_menu.items())
		
		
	class RemoveDataFromTable(object):
		def __call__(self,names,table_widget):
			if len(names) == 0: return # nothing happened
		
			from .emform import get_table_items_in_column
			entries = get_table_items_in_column(table_widget,0)
			text_entries = [table_widget.convert_text(str(i.text())) for i in entries]

			full_names = [table_widget.convert_text(name) for name in names]
	
			indices = [ text_entries.index(name) for name in full_names]
			indices.sort()
			indices.reverse()
			for idx in indices:
				table_widget.removeRow(idx)
				
	class AddDataToTable(object):
		def __init__(self,validator=None):
			self.validator = validator
			
		def __call__(self,list_of_names,table_widget):
		
			from .emselector import EMSelectorDialog
			selector = EMSelectorDialog(save_as_mode=False)
			
			if self.validator != None: 
				selector.set_validator(self.validator)
			files = selector.exec_()
			selector.close()
			
			if files != "":
				if isinstance(files,str): files = [files]
				
				from .emform import get_table_items_in_column
				entries = get_table_items_in_column(table_widget,0)
				entrie_tags = [base_name(str(i.text())) for i in entries]
				file_tags = [base_name(i) for i in files]
				error_messages = []
				for idx,tag in enumerate(file_tags):
					if tag in entrie_tags:
						error_messages.append("%s is already listed" %files[idx])
				
			
				if len(error_messages) > 0:
					EMErrorMessageDisplay.run(error_messages)
					return
				
				table_widget.add_entries(files)
		
	def get_quality_score(image_name):
		'''
		Used by the initial models table to get a quality score
		'''
		a = EMData()
		a.read_image(image_name,0,True)
		d = a.get_attr_dict()
		if "quality" in d: return "%.3f" %(d["quality"])
		else: return "-"
	
	get_quality_score = staticmethod(get_quality_score)
	
class CTFColumns(object):
	'''
	Basically some functions with a cache - the cache is to avoid
	re-reading stuff from disk multiple times
	'''
	def __init__(self):
		self.ctf_cache = {}
					
#		def __del__(self):
#			print "CTF columns dies"
		
	
	def get_ctf(self,name):
		ctf = None
		if name in self.ctf_cache:
			ctf = self.ctf_cache[name]
		else:
			try:
				a = EMData(name,0,True)
				d = a.get_attr_dict()
				ctf = d["ctf"]
				self.ctf_cache[name] = ctf
			except:
				self.ctf_cache[name] = None
		
		return ctf
	
	def get_defocus(self,name):
		ctf = self.get_ctf(name)
		try: return "%.3f" %ctf.defocus
		except: return ""
	
	def get_bfactor(self,name):
		ctf = self.get_ctf(name)
		try: return "%.3f" %ctf.bfactor
		except: return ""
			
	def get_sampling(self,name):
		ctf = self.get_ctf(name)
		try: return str(len(ctf.background))
		except: return ""
			
	def get_snr(self,name):
		ctf = self.get_ctf(name)
		if ctf != None:
			snr = 0
			try: snr = old_div(sum(ctf.snr),len(ctf.snr))
			except: pass
			return "%.3f" %snr
		else: return ""
		
	def get_quality(self,name):
		'''
		Quality is only applicable in the instance of the database
		'''
#		print name
		if db_check_dict("bdb:e2ctf.parms"):
			ctf_db = db_open_dict("bdb:e2ctf.parms",ro=True)
			try:
#				print name,base_name(name),ctf_db[base_name(name)]
				quality = ctf_db[base_name(name).split("_ctf")[0]][3]
				return "%d" %quality
			except:
				pass
		return "-"


class CTFDBColumns(CTFColumns):
	def __init__(self):
		CTFColumns.__init__(self)

	def get_ctf(self,name):
				ctf = None
				if name in self.ctf_cache:
					ctf = self.ctf_cache[name]
				else:
					if db_check_dict("bdb:e2ctf.parms"):
						ctf_db = db_open_dict("bdb:e2ctf.parms",ro=False)
						try:
							vals = ctf_db[base_name(name)][0]
							ctf = EMAN2Ctf()
							ctf.from_string(vals)
							self.ctf_cache[name] = ctf
						except:
							pass
				return ctf
			
	def get_quality(self,name):
		'''
		Quality is only applicable in the instance of the database
		'''
		if db_check_dict("bdb:e2ctf.parms"):
			ctf_db = db_open_dict("bdb:e2ctf.parms",ro=False)
			try:
				quality = ctf_db[base_name(name)][3]
				return "%d" %quality
			except:
				pass
		return "-"

class EMParticleReportTask(ParticleWorkFlowTask):
	'''This tool is for displaying the particles that are currently associated with this project.'''
	task_idle = QtCore.pyqtSignal()

	def __init__(self):
		ParticleWorkFlowTask.__init__(self)
		self.window_title = "Project Particles"
		self.project_list = spr_ptcls_dict
		self.project_data_at_init = None

	def get_project_particle_table(self):
		data_dict = EMProjectDataDict(spr_ptcls_dict)
		particle_data = data_dict.get_data_dict()
		particle_names = data_dict.get_names()
		
		self.project_data_at_init = particle_data # so if the user hits cancel this can be reset

		from .emform import EM2DStackTable,EMFileTable,int_lt
		table = EM2DStackTable(particle_names,desc_short="Project Particle Sets",desc_long="")
		context_menu_data = EMRawDataReportTask.ProjectListContextMenu(self.project_list,using_file_tags=False)		# STEVE, changed from true
		table.add_context_menu_data(context_menu_data)
		table.add_button_data(EMRawDataReportTask.ProjectAddRawDataButton(table,context_menu_data))
		table.insert_column_data(1,EMFileTable.EMColumnData("Particles On Disk",EMParticleReportTask.get_num_ptcls,"Particles currently stored on disk that are associated with this image",int_lt))
		table.insert_column_data(2,EMFileTable.EMColumnData("Particle Dims",EMParticleReportTask.get_particle_dims,"The dimensions of the particles that are stored on disk"))
		
		return table
	
	def on_form_cancel(self):
		self.recover_original_raw_data_list()
		
		self.form.close()
		self.form = None
		self.task_idle.emit()
	
	def recover_original_raw_data_list(self):
		'''
		Called if the user hits cancel - if they removed some files or added files the changes
		are not saved unless the user hits ok
		'''
		project_db = db_open_dict("bdb:project")
		project_db[self.project_list] = self.project_data_at_init
	
	def get_particle_dims(file_name):
		try:
			nx,ny,nz = gimme_image_dimensions3D(file_name)
			return "%ix%ix%i" %(nx,ny,nz)
		except:
			return "Error"
	
	def get_num_ptcls(file_name):
		return str(EMUtil.get_image_count(file_name))
	
	get_particle_dims = staticmethod(get_particle_dims)
	get_num_ptcls = staticmethod(get_num_ptcls)
	
	def get_params(self):
		params = []
		
	
		table = self.get_project_particle_table()
		
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=self.__doc__,choices=None))
		params.append(table)  
		
		return params
	
	
class E2CTFWorkFlowTask(EMParticleReportTask):
	'''
	Common functionality for E2CTF Work flow taskss
	'''
	def __init__(self):
		EMParticleReportTask.__init__(self)
		self.form_db_name = "bdb:emform.e2ctf"
	
	def get_ctf_param_table(self):
		'''
		
		'''		
		table = self.get_project_particle_table()
		
		from .emform import EMFileTable,float_lt,int_lt
		self.column_data = CTFDBColumns()
		table.add_column_data(EMFileTable.EMColumnData("Defocus",self.column_data.get_defocus,"The estimated defocus",float_lt))
		table.add_column_data(EMFileTable.EMColumnData("B Factor",self.column_data.get_bfactor,"The estimated B factor, note this is ~4x greater than in EMAN1",float_lt))
		table.add_column_data(EMFileTable.EMColumnData("SNR",self.column_data.get_snr,"The averaged SNR",float_lt))
		table.add_column_data(EMFileTable.EMColumnData("Quality",self.column_data.get_quality,"The quality of the fit as judged by e2ctf",int_lt))
		table.add_column_data(EMFileTable.EMColumnData("Sampling",self.column_data.get_sampling,"The amount of sampling used for generating CTF parameters",int_lt))
		return table,0
	
	def get_full_ctf_table(self,project_names=None,no_particles=False):
		'''
		Gets the ctf param table but also adds information about the wiener and phase flipped
		particles on disk (number, dimensions)
		'''
		table,n = self.get_ctf_param_table()
	
		self.other_column_data = E2CTFWorkFlowTask.MoreCTFColumns()
		from .emform import EMFileTable,int_lt
		table.add_column_data(EMFileTable.EMColumnData("Phase flip",self.other_column_data.get_num_phase_flipped,"The number of phase flipped particles on disk",int_lt))
		table.add_column_data(EMFileTable.EMColumnData("Phase flip dims",self.other_column_data.phase_flipped_dim,"The dimensions of the phase flippped particles"))
		table.add_column_data(EMFileTable.EMColumnData("Phase flip hp",self.other_column_data.get_num_phase_flipped_hp,"The number of phase flipped high pass filtered particles on disk",int_lt))
#		table.add_column_data(EMFileTable.EMColumnData("Phase flip hp dims",self.other_column_data.phase_flipped_dim,"The dimensions of the phase flippped particles"))
		table.add_column_data(EMFileTable.EMColumnData("Wiener filt",self.other_column_data.get_num_wein_filt,"The number of Wiener filtered particles on disk",int_lt))
#		table.add_column_data(EMFileTable.EMColumnData("Wiener filt dims",self.other_column_data.wien_filt_dim,"The dimensions of the Wiener filtered particles"))
		return table, n
	
	class MoreCTFColumns(object):
		'''
		Basically some functions with a cache - the cache is to avoid
		re-reading stuff from disk multiple times
		'''
		def __init__(self):
			
			data_dict = EMProjectDataDict(spr_ptcls_dict)
			self.db_map = data_dict.get_data_dict()
#			db = db_open_dict("bdb:project",ro=True)
#			self.db_map = db.get(,dfl={})
#			update = False
#			for name,map in self.db_map.items():
#				for key,image_name in map.items():
#					if not file_exists(image_name):
#						map.pop(key)
#						update = True
#				if len(map) == 0:
#					self.db_map.pop(name)
#					
#			if update:
#				EMErrorMessageDisplay.run("Warning, filtered particle data was lost.")
#				db["global.spr_filt_ptcls_map"] = self.db_map					
						
#		def __del__(self):
#			print "CTF columns dies"
		
		def __get_num_filtered(self,name,filt):
			
			tag = base_name(name)
			if tag in self.db_map:
				val = self.db_map[tag]
				if filt in val:
					file_name = val[filt]
					return str(EMUtil.get_image_count(file_name))
				
			return ""
		
		def __get_dim_filtered(self,name,filt):
			tag = base_name(name)
			if tag in self.db_map:
				val = self.db_map[tag]
				if filt in val:
					file_name = val[filt]
					nx,ny,nz = gimme_image_dimensions3D(file_name)
					return "%ix%ix%i" %(nx,ny,nz)

			return ""
		
		def get_num_phase_flipped(self,name):
			return self.__get_num_filtered(name,"Phase flipped")
				
		def phase_flipped_dim(self,name):
			return self.__get_dim_filtered(name,"Phase flipped")
		
		def get_num_phase_flipped_hp(self,name):
			return self.__get_num_filtered(name,"Phase flipped hp")
				
		def get_num_wein_filt(self,name):
			return self.__get_num_filtered(name,"Wiener filtered")
				
		def wien_filt_dim(self,name):
			return self.__get_dim_filtered(name,"Wiener filtered")
		
		
				
	def get_names_with_ctf_params(self):
		'''
		opens the e2ctf.parms directory and returns all a list of lists like this:
		[[db_name_key, real_image_name],....]
		eg
		[[neg_001,neg_001.hdf],[ptcls_01,bdb:particles#ptcls_01_ptcls],...] etc
		e2ctf is responsible for making sure the last data entry for each image is the original image name (this was first enforced by d.woolford)
		'''
		if not db_check_dict("bdb:e2ctf.parms"): return None
		parms_db = db_open_dict("bdb:e2ctf.parms",ro=True)
		
		ret = []
		for key,data in list(parms_db.items()):
			if data == None:
				print("error?",key)
				continue
			ret.append([key,data[-1]]) # parms[-1] should be the original filename
		#db_close_dict("bdb:e2ctf.parms")
		return ret

class E2CTFOutputTask(E2CTFWorkFlowTask):	
	"""Select the particle data for which you wish to generate phase flipped and/or Wiener filtered output and hit OK.\nThis will cause the workflow to spawn processes based on the available CPUs that write the output into a predefined location in the EMAN2 database.\nNote that the Wiener filtered output images are also phase flipped."""
	task_idle = QtCore.pyqtSignal()

	def __init__(self):
		E2CTFWorkFlowTask.__init__(self)
		self.window_title = "Generate e2ctf Output"

	def get_params(self):
		params = []		

		p,n = self.get_full_ctf_table()
		db = db_open_dict(self.form_db_name)
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=self.__doc__,choices=None))
		params.append(p)
		pos = ParamDef(name="oversamp",vartype="int",desc_short="Oversampling",desc_long="If greater than 1, data will be oversampled during phase flipping. Using this will make the flipping process irreversible and visibly depress the power spectra near zeroes.",property=None,defaultunits=db.get("oversampout",dfl=1),choices=None)
		pwiener = ParamDef(name="wiener",vartype="boolean",desc_short="Wiener",desc_long="Wiener filter your particle images using parameters in the database. Phase flipping will also occur",property=None,defaultunits=db.get("wiener",dfl=True),choices=None)
		pphase = ParamDef(name="phaseflip",vartype="boolean",desc_short="Phase flip",desc_long="Phase flip your particle images using parameters in the database",property=None,defaultunits=db.get("phaseflip",dfl=True),choices=None)
		pphasehp = ParamDef(name="phasefliphp",vartype="boolean",desc_short="Phase flip hp",desc_long="Phase flip your particle images and apply an automatic high-pass filter",property=None,defaultunits=db.get("phasefliphp",dfl=False),choices=None)
		params.append(pos)
		params.append([pphase,pphasehp,pwiener])
			#db_close_dict(self.form_db_name)
		
		return params

	def get_default_ctf_options(self,params):
		'''
		These are the options required to run pspec_and_ctf_fit in e2ctf.py, works in e2workflow
		'''
		
		error_message = []	
		
		options = EmptyObject()
	
#		filenames = params["filenames"]
##
#		db_file_names = []
#		for i,name in enumerate(filenames):
#			db_name=name
#			db_file_names.append(db_name)
#			if not file_exists(db_name):
#				print "error, particle entry doesn't exist for",name,"aborting."
#				return None
			
		options.filenames = params["filenames"]
		options.wiener = params["wiener"]
		options.phaseflip = params["phaseflip"]
		options.phasefliphp = params["phasefliphp"]
		options.oversamp = params["oversamp"]
		if not options.wiener and not options.phaseflip and not options.phasefliphp:
			error_message.append("Please choose at atleast one of the phaseflip or Wiener options.")
			
		if options.oversamp < 1:
			error_message.append("The oversampling factor must be atleast 1")
			
		if len(error_message) > 0:
			self.show_error_message(error_message)
			return None
#		
		return options
	
	def on_form_ok(self,params):
		
		if  "filenames" not in params or ("filenames" in params and len(params["filenames"]) == 0):
			error("Please select files to process")
			return

		options = self.get_default_ctf_options(params)
		if options != None and len(options.filenames) > 0 and (options.wiener or options.phaseflip or options.phasefliphp):
			self.write_db_entries(params)
			string_args = []
			bool_args = ["wiener","phaseflip","phasefliphp"]
#			additional_args = ["--dbds=%s"  %spr_ptcls_dict,"--computesf"] # don't want automatic computesf any more
			additional_args = ["--dbds=%s"  %spr_ptcls_dict]
			temp_file_name = "e2ctf_output_stdout.txt"
			self.spawn_single_task("e2ctf.py",options,string_args,bool_args,additional_args,temp_file_name)
			# Steve directed that the output writing task should use a single thead on July 3rd 2009
			#self.spawn_single_task("e2ctf.py",options,string_args,bool_args,additional_args,temp_file_name)
			

			self.form.close()
			self.task_idle.emit()
		else:
			return
		
class E2CTFOutputTaskGeneral(E2CTFOutputTask):
	''' Use this form for generating CTF-related output. 
	'''
	task_idle = QtCore.pyqtSignal()
	warning_string = "\n\n\nNOTE: There are no CTF parameters currently stored for any images in the local database. You can change this by running automated fitting with e2ctf."
	
	def __init__(self):
		E2CTFOutputTask.__init__(self)
		self.window_title = "Generate e2ctf Output"
		self.names = None
		
	def set_names(self,names):
		self.names = names

	def get_custom_table(self,names):
		
		from .emform import EM2DStackTable,EMFileTable,float_lt,int_lt
		table = EM2DStackTable(names,desc_short="Particle Images",desc_long="")
		context_menu_data = EMRawDataReportTask.ProjectListContextMenu(names,using_file_tags=True)
		table.add_context_menu_data(context_menu_data)
		table.add_button_data(EMRawDataReportTask.ProjectAddRawDataButton(table,context_menu_data))
	#	table.insert_column_data(1,EMFileTable.EMColumnData("Particles On Disk",EMParticleReportTask.get_num_ptcls,"Particles currently stored on disk that are associated with this image"))
		table.insert_column_data(2,EMFileTable.EMColumnData("Particle Dims",EMParticleReportTask.get_particle_dims,"The dimensions of the particles that are stored on disk"))
		self.column_data = CTFDBColumns()
		table.add_column_data(EMFileTable.EMColumnData("Defocus",self.column_data.get_defocus,"The estimated defocus",float_lt))
		table.add_column_data(EMFileTable.EMColumnData("B Factor",self.column_data.get_bfactor,"The estimated B factor, note this is ~4x greater than in EMAN1",float_lt))
		table.add_column_data(EMFileTable.EMColumnData("SNR",self.column_data.get_snr,"The averaged SNR",float_lt))
		table.add_column_data(EMFileTable.EMColumnData("Quality",self.column_data.get_quality,"The quality of the fit as judged by e2ctf",int_lt))
		table.add_column_data(EMFileTable.EMColumnData("Sampling",self.column_data.get_sampling,"The amount of sampling used for generating CTF parameters",int_lt))
		return table,len(names)
	

	def get_params(self):
		params = []
		
		if self.names != None: 
			names = self.names
			p,num = self.get_custom_table(names)
		else:
			p,num = self.get_ctf_param_table()
		
		
		if num == 0:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=self.__doc__,choices=None))
		else:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=self.__doc__,choices=None))
			params.append(p)
			pwiener = ParamDef(name="wiener",vartype="boolean",desc_short="Wiener",desc_long="Wiener filter your particle images using parameters in the database. Phase flipping will also occur",property=None,defaultunits=False,choices=None)
			pphase = ParamDef(name="phaseflip",vartype="boolean",desc_short="Phase flip",desc_long="Phase flip your particle images using parameters in the database",property=None,defaultunits=False,choices=None)
			pphasehp = ParamDef(name="phasefliphp",vartype="boolean",desc_short="Phase flip hp",desc_long="Phase flip your particle images and apply an automatic high-pass filter",property=None,defaultunits=False,choices=None)
		
			params.append([pphase,pphasehp,pwiener])
			
		return params
	
	def get_ctf_options(self,params):
		'''
		This is a way to get the ctf optiosn if one is using the "alternate" path, which means just in the context of general use of e2ctf (not the workflow)
		'''
		options = EmptyObject()
		
		if  "filenames" in params and len(params["filenames"]) == 0:
			return None # this is fine, for example, there were no files to work

		selected_filenames = params["filenames"]
		
		filenames = []
		names = self.get_names_with_ctf_params()
		for name in names:
			if name[0] in selected_filenames:
				filenames.append(name[1])
		options.filenames = filenames

		options.wiener = params["wiener"]
		options.phaseflip = params["phaseflip"]
#		
		return options

	def on_form_ok(self,params):
		for k,v in list(params.items()):
			self.write_db_entry(k,v)

		options = self.get_ctf_options(params)
		if options != None and len(options.filenames) > 0 and (options.wiener or options.phaseflip):
			
			string_args = []
			bool_args = ["wiener","phaseflip"]
			additional_args = []
			temp_file_name = "e2ctf_output_stdout.txt"
			self.spawn_single_task("e2ctf.py",options,string_args,bool_args,additional_args,temp_file_name)
			

			self.form.close()
			self.task_idle.emit()
		else:
			self.form.close()
			self.task_idle.emit()
			
	def on_form_cancel(self):
		
		self.form.close()
		self.form = None
		self.task_idle.emit()
	
class EMPartSetOptions(object):
	def __init__(self,data_dict_name,bdb_only=False):
		self.data_dict_name  = data_dict_name
		self.bdb_only = bdb_only # restricts returned options the image sets that exist only in the database
		self.image_count = False # if False merely counts the number of stacks - this is useful for the set choosing form
		
	def get_particle_options(self):
		'''
		@return stacks map - key is the name of the filtered set, value is a list of names in the that filtered set
		@return stacks_name_map - basic map - maps the selection string (e.g. "Original Data (2)" ) to the real name of the filtered set (e.g. "Original Data")
		@return choices a list of choices such as [ "Original Data (2)","Wiener Filtered (2)" ]
		@return name_map a way to map whichever name to the main key in the data base dictionary.e.g. key="bdb:particles#1100_ptcls_wienered", value="bdb:particles#1100_ptcls"
		'''
		data_dict = EMProjectDataDict(self.data_dict_name)

		stacks_map = {} # used to cache data in get_params
		stacks_name_map = {} # used to recall which selection was taken	
		db = db_open_dict("bdb:project")

		name_map = {} # a way to map a filtered image to its originating image
		filter_opts = {} # key is the filter type, value is the number of images with this filter type
		main_dict = data_dict.get_data_dict()

		for name,d in list(main_dict.items()):
			for filt,ptcl_name in list(d.items()):
				if self.bdb_only != False and ( len(ptcl_name) > 3 and ptcl_name[:4] != "bdb:"): 
					continue
				name_map[ptcl_name] = name
				if filt in filter_opts:
					if self.image_count: filter_opts[filt] += EMUtil.get_image_count(ptcl_name)
					else: filter_opts[filt] += 1
					stacks_map[filt].append(ptcl_name)
				else:
					if self.image_count: filter_opts[filt] = EMUtil.get_image_count(ptcl_name)
					else: filter_opts[filt] = 1
					stacks_map[filt] = [ptcl_name]
		
		choices = []
		for filt,num in list(filter_opts.items()):
			name = filt+" ("+str(num)+")"
			choices.append( name )
			stacks_name_map[name] = filt
		
		return stacks_map, stacks_name_map,choices,name_map
