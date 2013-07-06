#!/usr/bin/env python
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

from emform import EMFormWidget,EMParamTable,EMTableFormWidget
from emdatastorage import ParamDef
from PyQt4 import QtGui,QtCore
from PyQt4.QtCore import Qt
from EMAN2db import db_check_dict, db_open_dict,db_remove_dict,db_list_dicts,db_close_dict, e2getcwd
from EMAN2 import *
import os
import copy
from emapplication import EMProgressDialog, get_application, EMErrorMessageDisplay, error
from e2ctf import pspec_and_ctf_fit,GUIctf,write_e2ctf_output,get_gui_arg_img_sets,init_sfcurve
import subprocess
import weakref
from e2history import HistoryForm
import time
from emsave import save_data
from emimagemx import EMDataListCache
import traceback

USING_RELATIVE_DIRS = True # used by database infrastructure for recording file names

def workflow_path(path="",dir=None):
	'''
	Makes paths relative if possible.
	'''
	
	if USING_RELATIVE_DIRS:
		if dir==None: dir = e2getcwd()
		from emselector import folderize
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

class EmptyObject:
	'''
	This just because I need an object I can assign attributes to, and object() doesn't seem to work
	'''
	def __init__(self):
		pass

class WorkFlowTask:
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
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_ok"),self.on_form_ok)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_cancel"),self.on_form_cancel)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_close"),self.on_form_close)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("display_file"),self.on_display_file)
	
	def get_params(self): raise NotImplementedError
	
	def on_display_file(self,filename):
		self.emit(QtCore.SIGNAL("display_file"),filename)	
		
	def on_form_ok(self,params):
		for k,v in params.items():
			self.write_db_entry(k,v)
		
		self.disconnect_form()
		self.form.close()
		self.form = None
	
		self.emit(QtCore.SIGNAL("task_idle"))
		
	def on_form_cancel(self):
		self.disconnect_form()
		self.form.close()
		self.form = None
		self.emit(QtCore.SIGNAL("task_idle"))
	
	def disconnect_form(self):
		QtCore.QObject.disconnect(self.form,QtCore.SIGNAL("emform_ok"),self.on_form_ok)
		QtCore.QObject.disconnect(self.form,QtCore.SIGNAL("emform_cancel"),self.on_form_cancel)
		QtCore.QObject.disconnect(self.form,QtCore.SIGNAL("emform_close"),self.on_form_close)
		QtCore.QObject.disconnect(self.form,QtCore.SIGNAL("display_file"),self.on_display_file)
	
	
	def emit(self,*args,**kargs):
		self.core_object.emit(*args,**kargs)
		
	def on_form_close(self):
		self.disconnect_form()
		self.emit(QtCore.SIGNAL("task_idle"))

	def close(self):
		if self.form != None: 
			self.form.close()
		self.emit(QtCore.SIGNAL("task_idle"))
        
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
		for k,v in dictionary.items():
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
		specificially for the purpose of spawning processes. Could be used more generally, however.
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
		cf = float(len(options.filenames))/float(ncpu) # common factor
		
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
			for i in args: print i,
			print
			
			#print args
			#fname = temp_fname_root +"_"+str(n)+temp_fname_end
			#file = open(fname,"w+")
			
			cmdstr = ' '.join(args)
#			process = subprocess.Popen(args_adjusted,stdout=file,stderr=subprocess.STDOUT)
			process = subprocess.Popen(cmdstr, shell=True)
			print "started process",process.pid
			self.emit(QtCore.SIGNAL("process_started"),process.pid)
			
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
		for i in args: print i,
#		print ""
#		
		#print args
		#file = open(temp_file_name,"w+")

		cmdstr = ' '.join(args)
#		process = subprocess.Popen(args_adjusted,stdout=file,stderr=subprocess.STDOUT)
		process = subprocess.Popen(cmdstr, shell=True)
		print "started process",process.pid
		self.emit(QtCore.SIGNAL("process_started"),process.pid)
		
		#db_close_dict("bdb:project")
		
	def run_select_files_msg(self):
		'''
		Runs a QMessageBox asking for the user to select files for processing
		'''
		msg = QtGui.QMessageBox()
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
		# allright everything left in dirs is "r2d_??" where the ?? is castable to an int, so we should be safe now
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
						#else just check for 01 incase the user has specified the --initial arugment
				if not cont:
					break
				
			if classes_db != None:
				class_files.append(classes_db)
				cl_db = db_open_dict(classes_db,ro=True)
				if cl_db.has_key("maxrec"):
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
		return dump_cmps_list().keys()
	
	def get_aligners_list(self):
		return dump_aligners_list().keys()
	
	def get_projectors_list(self):
		return dump_projectors_list().keys()
	
	def get_orientgens_list(self):
		return dump_orientgens_list().keys()
		
	def get_averagers_list(self):
		return dump_averagers_list().keys()
		
#		cmps.append("None") I think this is necessary
		
class HistoryTask(WorkFlowTask,HistoryForm):
	def __init__(self):
		WorkFlowTask.__init__(self)
		# don't need HistoryFrom init
		self.wd = os.getcwd()
		self.window_title = "History"
	
	def run_form(self):	
		self.form = EMFormWidget(self.get_history_table())
		self.form.resize(*self.preferred_size)
		self.form.setWindowTitle(self.window_title)
		self.form.setWindowIcon(QtGui.QIcon(os.getenv("EMAN2DIR")+"/images/feather.png"))
		get_application().show_specific(self.form)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_ok"),self.on_form_ok)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_cancel"),self.on_form_cancel)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_close"),self.on_form_close)
		
class ChangeDirectoryTask(WorkFlowTask):
	def __init__(self):
		WorkFlowTask.__init__(self)
		self.window_title = "Change Project Directory"
	
	def run_form(self):	
		
		fsp=QtGui.QFileDialog.getExistingDirectory(None, "Choose a directory")
		fsp = str(fsp)
		if os.path.exists(fsp):
			os.chdir(fsp)
			return fsp
		else: return None
		
	def closeEvent(self,event):
		pass
		#self.form.close()

class EMProjectDataDict:
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
		return dict.keys()
		
	def get_names(self,filt=original_data):
		'''
		@param filt the key to retrieve from the filter dictionary
		@return a list of names that have entries corresponding to this filter entry
		'''
		dict = self.get_data_dict()
		ret = []
		for tag,dict2 in dict.items():
			for key,name in dict2.items():
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
		for name, map in dict.items():			
			for key,image_name in map.items():
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
		if not project_db.has_key(self.data_dict_name) or project_db[self.data_dict_name]==None or len(project_db[self.data_dict_name]) == 0:
			if project_db.has_key(self.recovery_dict_name):
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
		for tag,dict2 in dict.items():
			if not dict2.has_key(EMProjectDataDict.original_data):
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
		if self.recovery_list_name == None or not project_db.has_key(self.recovery_list_name): return
		
		project_data  = project_db.get(self.data_dict_name,dfl={})
		old_project_data  = project_db.get(self.recovery_list_name,dfl=[])
		if len(old_project_data) > 0:
			EMErrorMessageDisplay.run("Recovering %s database (for back compatibility)" %self.recovery_list_name, "Warning")
			for name in old_project_data:
				if not project_data.has_key(name):
					project_data[name] = {}
			project_db[self.data_dict_name] = project_data
			
		project_db[ self.recovery_list_name] = None

	
class SPRInitTask(WorkFlowTask):
	'''Welcome to the EMAN2 workflow. Use this tool to step through and manage the process of generating single particle reconstructions. Get started by entering what you can of the parameters in this form and then proceed to the next step task in the workflow'''
	def __init__(self):
		WorkFlowTask.__init__(self)
		self.window_title = "Project Settings"
	def get_params(self):
		params = []
		project_db = db_open_dict("bdb:project")
		params.append(ParamDef(name="blurb",vartype="text",desc_short="SPR",desc_long="Information regarding this task",property=None,defaultunits=self.__doc__,choices=None))
		
		papix = ParamDef(name="global.apix",vartype="float",desc_short="A/pix for project",desc_long="The physical distance represented by the pixel spacing",property=None,defaultunits=project_db.get("global.apix",dfl=1.1),choices=None)
		pmass = ParamDef(name="global.particle_mass",vartype="float",desc_short="Particle mass (kda)",desc_long="The mass of the particle in kilodaltons",property=None,defaultunits=project_db.get("global.particle_mass",dfl=800),choices=None)
		
		pvolt = ParamDef(name="global.microscope_voltage",vartype="float",desc_short="Microscope voltage",desc_long="The operating voltage of the microscope in kilo volts",property=None,defaultunits=project_db.get("global.microscope_voltage",dfl=300),choices=None)
		pcs = ParamDef(name="global.microscope_cs",vartype="float",desc_short="Microscope Cs",desc_long="Microscope spherical aberration constant",property=None,defaultunits=project_db.get("global.microscope_cs",dfl=2.0),choices=None)
		pncp = ParamDef(name="global.num_cpus",vartype="int",desc_short="Number of CPUs",desc_long="Number of CPUS available for the project to use",property=None,defaultunits=project_db.get("global.num_cpus",dfl=num_cpus()),choices=None)
		mem = memory_stats()
		pmem = ParamDef(name="global.memory_available",vartype="float",desc_short="Memory usage (%.2f Gb total)" %mem[0],desc_long="The total amount of system memory you want to make available to the project in gigabytes",property=None,defaultunits=project_db.get("global.memory_available",dfl=mem[1]),choices=None)
		params.append(pmass)
		params.append(papix)
		params.append(pvolt)
		params.append(pcs)
		params.append(pncp)
		params.append(pmem)
		#db_close_dict("bdb:project")
		return params

	def write_db_entry(self,key,value):
		WorkFlowTask.write_db_entry(self,key,value)		


class EMRawDataReportTask(WorkFlowTask):	
	'''This form displays raw data that are associated with the project. You browse to add raw data, or right click and choose Add.''' 
	documentation_string = "This page shows raw micrographs/ccd frames currently associated with the project. It is possible to add additional images directly on this panel, which will \
leave them in-place and not copy them into the project database. This will limit some later operations and leave the project with less metadata at the end, but will save disk space. \
Note that the data cannot be filtered unless it is imported."
	warning_string = "\n\n\nNOTE: There are no images currenty associated with the project. Please associate or import images"
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
		print self.project_data_at_init
		project_names = data_dict.keys()
		
		from emform import EM2DFileTable,EMFileTable
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
		from emform import EMFileTable
		table.add_column_data(EMFileTable.EMColumnData("Dimensions",EMRawDataReportTask.get_image_dimensions,"The dimensions of the file on disk"))
		return table,n
	
	class ProjectAddRawDataButton:
		def __init__(self,table_widget,context_menu_data):
			self.table_widget = weakref.ref(table_widget)
			self.context_menu_data = context_menu_data
			self.name = "Browse To Add"
			
		def function(self,bool):
			self.context_menu_data.context_menu["Add"]([],self.table_widget())
			

	class ProjectListContextMenu:
		def __init__(self,project_list=spr_raw_data_dict,remove_only=False,using_file_tags=False):
			self.project_list = project_list
			self.validator = AddFilesToProjectValidator(self.project_list)
			self.context_menu = {}
			self.context_menu["Remove"] = EMRawDataReportTask.ProjectListContextMenu.RemoveFilesFromProject(self.project_list,using_file_tags)
			if not remove_only: self.context_menu["Add"] = EMRawDataReportTask.ProjectListContextMenu.AddFilesToProjectViaContext(self.project_list)
		
		def items(self):
			return self.context_menu.items()
		
		
		class RemoveFilesFromProject:
			def __init__(self,project_list,using_file_tags=False):
				self.project_list = project_list
				self.using_file_tags = using_file_tags
			def __call__(self,names,table_widget):
				if len(names) == 0: return # nothing happened
			
				from emform import get_table_items_in_column
				entries = get_table_items_in_column(table_widget,0)
				text_entries = [table_widget.convert_text(str(i.text())) for i in entries]
				
				data_dict = EMProjectDataDict(self.project_list)
				project_names = data_dict.keys()
				
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
				
		class AddFilesToProject:
			def __init__(self,project_list):
				self.project_list = project_list
				
			def __call__(self,list_of_names,table_widget):
		
				data_dict = EMProjectDataDict(self.project_list)
				project_names = data_dict.keys()
				
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
		
		class AddFilesToProjectViaContext:
			def __init__(self,project_list):
				self.project_list = project_list
				self.validator = AddFilesToProjectValidator(self.project_list)
				
			def __call__(self,list_of_names,table_widget):
			
		#def add_files_from_context_menu(self,list_of_names,table_widget):
				from emselector import EMSelectorDialog
				selector = EMSelectorDialog(save_as_mode=False)
				
				selector.set_selection_text("Selection(s)")
				selector.set_validator(self.validator)
				files = selector.exec_()
				selector.close()

				if files != "":
					if isinstance(files,str): files = [files]
					
					from emform import get_table_items_in_column
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
		self.emit(QtCore.SIGNAL("task_idle"))
	
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
		

		
class AddFilesToProjectValidator:
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
		project_names = data_dict.keys()
		
		for name in list_of_names:
			if not file_exists(name) and not os.path.isdir(name[:name.rfind("/")]): #See if its a dir, which is ok:
				EMErrorMessageDisplay.run(["%s does not exists" %name])
				return 0
		
		for name in list_of_names:
			if name in project_names:
				EMErrorMessageDisplay.run(["%s is already in the project" %name])
				return 0
			
		return 1

class EMFilterRawDataTask(WorkFlowTask):	
	documentation_string = """This tool allows you to import micrographs/ccd frames into the project. This copies the files into the internal \
project database, and gives an opportunity to apply a number of common filters to the data before importing:
- Invert - EMAN2 expects particle data to be positive, ie - particles should appear white. If particles are dark, select this.
- Filter Xray Pixels - Use this for CCD frames. It will remove over-bright individual pixels
- Associate with project - always select this
- In-place processing - technically allows processing without importing to the database. Not suggested.s"""
	def __init__(self):
		WorkFlowTask.__init__(self)
		self.window_title = "Filter Raw Data"
		self.thumb_shrink = -1
		self.output_names = [] # eventually populated with output names. Could be the names of the files themselves if processing is in place
		self.form_db_name = "bdb:emform.filter_raw_data"
		
	def get_table(self,ptcl_list=[]):
		from emform import EM2DFileTable,EMFileTable,float_lt
		table = EM2DFileTable(ptcl_list,desc_short="Raw Data",desc_long="")
		context_menu_data = ParticleWorkFlowTask.DataContextMenu()
		table.add_context_menu_data(context_menu_data)
		table.add_button_data(ParticleWorkFlowTask.AddDataButton(table,context_menu_data))
		table.add_column_data(EMFileTable.EMColumnData("Dimensions",EMRawDataReportTask.get_image_dimensions,"The dimensions of the file on disk"))
		table.add_column_data(EMFileTable.EMColumnData("Mean",EMFilterRawDataTask.get_mean,"The mean pixel value",float_lt))
		table.add_column_data(EMFileTable.EMColumnData("Sigma",EMFilterRawDataTask.get_sigma,"The standard deviation of the pixel values",float_lt))
	
		return table, len(ptcl_list)
	
	def get_mean(file_name,idx=0):
		'''
		A static function for getting the dimensions of a file as a string
		'''
		a = EMData()
		a.read_image(file_name,idx,True)
		try: val = a.get_attr("mean")
		except: return ""
		if val != None:	return str(val)
		else: return ""
	
	def get_sigma(file_name,idx=0):
		'''
		A static function for getting the dimensions of a file as a string
		'''
		a = EMData()
		a.read_image(file_name,idx,True)
		try: val = a.get_attr("sigma")
		except: return ""
		if val != None:	return str(val)
		else: return ""
	
	get_mean = staticmethod(get_mean)
	get_sigma = staticmethod(get_sigma)

	def get_params(self):
#		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2CTFGenericTask.documentation_string,choices=None))
		
		db = db_open_dict(self.form_db_name)
		
		p,n = self.get_table([])
		params = []
		project_db = db_open_dict("bdb:project")
		params.append(ParamDef(name="blurb",vartype="text",desc_short="Filtering Raw Data",desc_long="",property=None,defaultunits=EMFilterRawDataTask.documentation_string,choices=None))
		pinvert = ParamDef(name="invert",vartype="boolean",desc_short="Invert",desc_long="Invert pixel intensities",property=None,defaultunits=db.get("invert",dfl=False),choices=None)
		pxray = ParamDef(name="xraypixel",vartype="boolean",desc_short="Filter X-ray pixels",desc_long="Filter X-ray pixels (4*sigma)",property=None,defaultunits=db.get("xraypixel",dfl=True),choices=None)
		pnorm = ParamDef(name="norm.edgemean",vartype="boolean",desc_short="Edge norm",desc_long="Normalize using the normalize.edgemean processor",property=None,defaultunits=db.get("norm.edgemean",dfl=True),choices=None)
		pthumbnail = ParamDef(name="thumbs",vartype="boolean",desc_short="Generate thumbnail",desc_long="Generate thumbnails for e2boxer",property=None,defaultunits=db.get("thumbs",dfl=True),choices=None)
		passociate = ParamDef(name="project_associate",vartype="boolean",desc_short="Associate with project",desc_long="Associate filtered images with the project",property=None,defaultunits=db.get("project_associate",dfl=True),choices=None)
		pinplace = ParamDef(name="inplace",vartype="boolean",desc_short="Inplace processing",desc_long="Process images inplace to save disk space",property=None,defaultunits=db.get("inplace",dfl=False),choices=None)
		psuffix = ParamDef(name="suffix",vartype="string",desc_short="Output Suffix", desc_long="This text will be appended to the names of the output files",property=None,defaultunits=db.get("suffix",dfl=""),choices=None )
	
		pinplace.dependents = ["suffix"] # these are things that become disabled when the pwb checkbox is unchecked etc
		pinplace.invert_logic = True
		params.append(p)
		params.append([pinvert,pxray,pnorm])
		params.append([pthumbnail,passociate])
		params.append([pinplace,psuffix])
		
		#db_close_dict("bdb:project")
		return params
	
	def on_form_ok(self,params):
		
		error_message = self.check_params(params)
		if len(error_message):
			self.show_error_message(error_message)
			return
 		
		self.do_filtering(params,self.output_names)

		self.form.close()
		self.form = None
	
		self.emit(QtCore.SIGNAL("task_idle"))
		
		self.write_db_entries(params)
		
	def check_params(self,params):
		error_message = []
		filenames = params["filenames"]
		
		if len(filenames) == 0:
			error_message.append("Please specify files to filter.")
			return error_message
		
		import copy
		self.output_names = copy.deepcopy(filenames)
		
		num_proc_ops = 0
		num_thm_ops = 0
		if params["invert"]: num_proc_ops += 1
		if params["xraypixel"]: num_proc_ops += 1
		if params["thumbs"]:num_thm_ops += 1
		if params["norm.edgemean"]:num_proc_ops += 1
		
		if num_thm_ops+ num_proc_ops== 0: error_message.append("Please choose an operation to perform")
		
		if num_proc_ops == 0 and params["inplace"]:
			error_message.append(["Inplace filtering is only valid if you supply a processing operation"])
		
		if not params["inplace"]:
			reps = []
			for file in filenames:
				if len(file) > 3 and file[:4] == "bdb:":
					reps.append(file+params["suffix"])
				else:
					idx = file.rfind(".")
					if idx == -1: error_message.append("Couldn't interpret %s" %file)
					else:
						type_ = file[idx+1:]
						if type_ not in get_supported_2d_write_formats(): type_ = "mrc"
						reps.append(file[:idx]+params["suffix"]+"."+type_)
			self.output_names = reps
			
			for name in self.output_names:
				if file_exists(name):
					error_message.append("The output file %s already exists. Please remove it or change your suffix." %name)
		else:
			for file in self.output_names:
				if len(file) > 3 and file[:4] == "bdb:": continue
				else:
					idx = file.rfind(".")
					if idx == -1: error_message.append("Couldn't interpret %s" %file)
					else:
						type_ = file[idx+1:]
						if type_ not in get_supported_2d_write_formats():
							error_message.append("Can not inplace filter %s files (%s)" %(type_,file))		
						
		for name in filenames:
			if len(name) == 0: continue
			if not file_exists(name):
				error_message.append("File %s doesn't exist." %name)
			else:
				try:
					e = EMData()
					e.read_image(name,0,1)
					n = EMUtil.get_image_count(name)
					if n > 1:
						error_message.append("File %s contains more than 1 image." %name)
				except:
					error_message.append("File %s is not a valid EM image." %name)
				
		data_dict = EMProjectDataDict(spr_raw_data_dict)
		project_names = data_dict.keys()
		if params["project_associate"]:
			for name in self.output_names:
				# this will change one day - the imported name will be altered
				if name in project_names:
					error_message.append("%s is already associated with the project" %name)
					
		return error_message

	def do_filtering(self,params,output_names):
		
		filenames = params["filenames"]
		
		# get the number of process operation - the progress dialog reflects image copying and image processing operations
		num_processing_operations = 2 # there is atleast a copy and a disk write
		if params["invert"]: num_processing_operations += 1
		if params["xraypixel"]: num_processing_operations += 1
		if params["thumbs"]:num_processing_operations += 1
		if params["norm.edgemean"]:num_processing_operations += 1
		
		# now add the files to db (if they don't already exist
		progress = EMProgressDialog("Filtering Raw Files", "Abort", 0, len(filenames)*num_processing_operations,None)
		progress.show()
		i = 0
		cancelled = False # if the user cancels the import then we must act
		cancelled_writes = []
		get_application().setOverrideCursor(Qt.BusyCursor)
		for j in xrange(0,len(filenames)):
			name = filenames[j]
			outname = output_names[j]
			
			e = EMData()
			e.read_image(name,0)
			i += 1
			progress.setValue(i)	
			get_application().processEvents()
			e.set_attr("disk_file_name",name)
			
			write_large = False
			if params["norm.edgemean"]:
				e.process_inplace("normalize.edgemean")
				i += 1
				progress.setValue(i)
				get_application().processEvents()
				write_large = True
			
			if params["invert"]:
				e.mult(-1)
				i += 1
				progress.setValue(i)
				get_application().processEvents()
				write_large = True
			
			if params["xraypixel"]:
				e.process_inplace("threshold.clampminmax.nsigma",{"nsigma":4,"tomean":True})
				i += 1
				progress.setValue(i)
				get_application().processEvents()
				write_large = True
			
			if write_large:	e.write_image(outname,0)
			#db_close_dict(db_name)
			cancelled_writes.append(outname)
			i += 1
			progress.setValue(i)
			get_application().processEvents()
				
			if params["thumbs"]:
				shrink = self.get_thumb_shrink(e.get_xsize(),e.get_ysize())
				thumb = e.process("math.meanshrink",{"n":shrink})
				thumb.process_inplace("normalize.edgemean")
				from emboxerbase import set_idd_image_entry
				if not write_large: outname = name
				set_idd_image_entry(outname,"image_thumb",thumb) # boxer uses the full name
				i += 1
				progress.setValue(i)
				get_application().processEvents()
				
			if progress.wasCanceled():
				cancelled = True
				break
			
		progress.setValue(len(filenames))
		progress.close()
		get_application().setOverrideCursor(Qt.ArrowCursor)
		
		if not cancelled and params["project_associate"]:
			data_dict = EMProjectDataDict(spr_raw_data_dict)
			data_dict.add_names(output_names)
		
	def get_thumb_shrink(self,nx,ny):
		if self.thumb_shrink == -1:
			shrink = 1
			inx =  nx/2
			iny =  ny/2
			while ( inx >= 128 and iny >= 128):
				inx /= 2
				iny /= 2
				shrink *= 2
		
			self.thumb_shrink=shrink
		
		return self.thumb_shrink
			
	def on_import_cancel(self):
		print "canceled"
		
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
		from emform import EM2DStackTable,EMFileTable,float_lt,int_lt
		if table==None: # you can hand in your own table (it's done in E2ParticleExamineTask)
			table = EM2DStackTable(ptcl_list,desc_short="Particles",desc_long="",single_selection=single_selection)
		if len(ptcl_list) != 0:
			a = EMData()
			a.read_image(ptcl_list[0],0,True)
			d = a.get_attr_dict()
			if enable_ctf and d.has_key("ctf"):
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

	
	class AddDataButton:
		def __init__(self,table_widget,context_menu_data):
			self.table_widget = weakref.ref(table_widget)
			self.context_menu_data = context_menu_data
			self.name = "Browse To Add"
			
		def function(self,bool):
			self.context_menu_data.context_menu["Add"]([],self.table_widget())
			

	class DataContextMenu:
		def __init__(self,validator=None):
			
			self.validator = validator
			self.context_menu = {}
			self.context_menu["Remove"] = ParticleWorkFlowTask.RemoveDataFromTable()
			self.context_menu["Add"] = ParticleWorkFlowTask.AddDataToTable(validator)
		
		def items(self):
			return self.context_menu.items()
		
		
	class RemoveDataFromTable:
		def __call__(self,names,table_widget):
			if len(names) == 0: return # nothing happened
		
			from emform import get_table_items_in_column
			entries = get_table_items_in_column(table_widget,0)
			text_entries = [table_widget.convert_text(str(i.text())) for i in entries]

			full_names = [table_widget.convert_text(name) for name in names]
	
			indices = [ text_entries.index(name) for name in full_names]
			indices.sort()
			indices.reverse()
			for idx in indices:
				table_widget.removeRow(idx)
				
	class AddDataToTable:
		def __init__(self,validator=None):
			self.validator = validator
			
		def __call__(self,list_of_names,table_widget):
		
			from emselector import EMSelectorDialog
			selector = EMSelectorDialog(save_as_mode=False)
			
			if self.validator != None: 
				selector.set_validator(self.validator)
			files = selector.exec_()
			selector.close()
			
			if files != "":
				if isinstance(files,str): files = [files]
				
				from emform import get_table_items_in_column
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
		if d.has_key("quality"): return "%.3f" %(d["quality"])
		else: return "-"
	
	get_quality_score = staticmethod(get_quality_score)
	
class CTFColumns:
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
		if self.ctf_cache.has_key(name):
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
			try: snr = sum(ctf.snr)/len(ctf.snr)
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
				if self.ctf_cache.has_key(name):
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

def ptable_convert_2(text):
	'''
	This is needed in one location - very close to removal
	'''
	return text

class EMParticleReportTask(ParticleWorkFlowTask):
	'''This tool is for displaying the particles that are currently associated with this project.'''
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

		from emform import EM2DStackTable,EMFileTable,int_lt
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
		self.emit(QtCore.SIGNAL("task_idle"))
	
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
	
	
class EMParticleImportTask(ParticleWorkFlowTask):
	'''This task is for importing particle into the project. The data will be copied into the particles directory. This is essential if you wish to use the data to generating sets'''
	def __init__(self):
		ParticleWorkFlowTask.__init__(self)
		self.window_title = "Import Particles"

	def get_project_particle_table(self):
		data_dict = EMProjectDataDict(spr_ptcls_dict)
		particle_data = data_dict.get_data_dict() # this is for back compatibility only - it cleans up old databases
		
		from emform import EM2DStackTable,EMFileTable,int_lt
		table = EM2DStackTable([],desc_short="Particles",desc_long="")
		context_menu_data = EMParticleImportTask.ContextMenu(spr_ptcls_dict)
		table.add_context_menu_data(context_menu_data)
		table.add_button_data(EMParticleImportTask.ProjectAddRawDataButton(table,context_menu_data))
		table.insert_column_data(1,EMFileTable.EMColumnData("Particles On Disk",EMParticleReportTask.get_num_ptcls,"Particles currently stored on disk that are associated with this image",int_lt))
		table.insert_column_data(2,EMFileTable.EMColumnData("Particle Dims",EMParticleReportTask.get_particle_dims,"The dimensions of the particles that are stored on disk"))
		
		return table

	def get_params(self):
		params = []
		
		table = self.get_project_particle_table()
		
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=self.__doc__,choices=None))
		params.append(table)  
		
		return params
	
	def update_params(self,params):
		warning = []
		name_map = {}
		for name in params["filenames"]:
			# Changed by Steve on 9/14/10. Used to just make a reference. Changed to actually copy data.
			#if len(name) > 3 and name[:4] == "bdb:":
				#warning.append(name)
				#name_map[name] = name
			#else:
			
			new_name = "bdb:particles#"+base_name(name)
			if file_exists(new_name):
				i = 0
				while 1:
					c_name = new_name + "_" + str(i)
					if not file_exists(c_name):
						new_name = c_name
						break
						
			name_map[name] = new_name
		params["name_map"] = name_map
		
		if len(warning) > 0:
			s = "The files: "
			for w in warning:
				s += w
				if w != warning[-1]: s+= ", "
			s += " are already in database format: they will merely be associated with the project"
			error(s,"Warning")
		
		print params["name_map"]
	def on_form_ok(self,params):
		
		if  params.has_key("filenames") and len(params["filenames"]) == 0:
			self.run_select_files_msg()
			return
		
		self.update_params(params)
		success, cmd = self.import_data(params)
		if not success:
			s = "The command: " + cmd, "failed to complete"
			error(s)
			return
		else:
			data_dict = EMProjectDataDict(spr_ptcls_dict)
			print params["name_map"].values()
			data_dict.add_names(params["name_map"].values(),use_file_tag=True)
		
		self.emit(QtCore.SIGNAL("task_idle"))
		self.form.close()
		self.form = None
		
	def import_data(self,params):

		get_application().setOverrideCursor(Qt.BusyCursor)
		progress = QtGui.QProgressDialog("Importing files into database...", "Abort import", 0, len(params["name_map"]),None)
		progress.show()
		i = 0
		progress.setValue(i)
		get_application().processEvents()

		
		for infile,output in params["name_map"].items():
			#if len(input) > 3 and infile[:4] == "bdb:": 
				#i += 1
				#progress.setValue(i)
 				#get_application().processEvents()
				#continue
			
			cmd = "e2proc2d.py"
 			cmd += " "+infile
 			cmd += " "+output
 			success = (os.system(cmd) in (0,12))
 			if not success or progress.wasCanceled():
 				progress.close()
 				return False,cmd
 			else:
 				i += 1
 				progress.setValue(i)
 				get_application().processEvents()
		
		progress.close()
		
		get_application().setOverrideCursor(Qt.ArrowCursor)
		return True,"success"
	
	
	class ProjectAddRawDataButton:
		def __init__(self,table_widget,context_menu_data):
			self.table_widget = weakref.ref(table_widget)
			self.context_menu_data = context_menu_data
			self.name = "Browse To Add"
			
		def function(self,bool):
			self.context_menu_data.context_menu["Add"]([],self.table_widget())
			

	class ContextMenu:
		def __init__(self,project_list):
			self.project_list = project_list
			self.validator = AddFilesToProjectValidator(self.project_list)
			self.context_menu = {}
			self.context_menu["Remove"] = EMParticleImportTask.ContextMenu.RemoveFiles()
			self.context_menu["Add"] = EMParticleImportTask.ContextMenu.AddFilesViaContext(self.project_list)
		
		def items(self):
			return self.context_menu.items()
		
		
		class RemoveFiles:
			def __init__(self): pass
			def __call__(self,names,table_widget):
				if len(names) == 0: return # nothing happened
			
				from emform import get_table_items_in_column
				entries = get_table_items_in_column(table_widget,0)
				text_entries = [table_widget.convert_text(str(i.text())) for i in entries]
				
				indices = [ text_entries.index(name) for name in full_names]
				indices.sort()
				indices.reverse()
				for idx in indices:
					table_widget.removeRow(idx)
					
		class AddFiles:
			def __init__(self,project_list):
				self.project_list = project_list
				
			def __call__(self,list_of_names,table_widget):
#				print list_of_names
				project_db = db_open_dict("bdb:project")
				project_names = project_db.get(self.project_list,dfl=[])
				project_name_tags = [base_name(name) for name in project_names]
#				print project_names,project_name_tags

				for name in list_of_names:
					if not file_exists(name):
						EMErrorMessageDisplay.run(["%s does not exists" %name])
						return
				
#				for name in list_of_names:
#					if base_name(name) in project_name_tags:
#						EMErrorMessageDisplay.run(["%s is already in the project" %name])
#						return

				err_list=[i for i in list_of_names if (base_name(i) in project_name_tags)]
				if len(err_list)>0 : EMErrorMessageDisplay.run(["%s are already in the project"%(",".join(err_list))])

				# only include files not already in the project
				list_of_names=[i for i in list_of_names if not (base_name(i) in project_name_tags)]
						
				# if we make it here we're good
				# first add entries to the table
				table_widget.add_entries(list_of_names)
				
		
		class AddFilesViaContext:
			def __init__(self,project_list):
				self.project_list = project_list
				self.validator = AddFilesToProjectValidator(self.project_list)
				
			def __call__(self,list_of_names,table_widget):
			
		#def add_files_from_context_menu(self,list_of_names,table_widget):
				from emselector import EMSelectorDialog
				selector = EMSelectorDialog(save_as_mode=False)
				
				selector.set_selection_text("Selection(s)")
				selector.set_validator(self.validator)
				files = selector.exec_()
				selector.close()
#				print files
				if files != "":
					if isinstance(files,str): files = [files]
					
					from emform import get_table_items_in_column
					entries = get_table_items_in_column(table_widget,0)
					entrie_tags = [base_name(str(i.text())) for i in entries]
					file_tags = [base_name(i) for i in files]
#					print file_tags
					error_messages = []
					for idx,tag in enumerate(file_tags):
						if tag in entrie_tags:
							error_messages.append("%s is already listed" %files[idx])
					
				
					if len(error_messages) > 0:
						EMErrorMessageDisplay.run(error_messages)
						return
				a = EMParticleImportTask.ContextMenu.AddFiles(self.project_list)
				a(files,table_widget)


class EMParticleCoordImportTask(WorkFlowTask):
	"""Use this task for importing coordinate files directly into the project database. Generally you do this after you have associated raw data with the project. Coordinate files are matched with input (raw data) files according to name. For example, if your micrograph is called "321a.mrc" then your coordinate file should be called "321a.box" or "321a.data" (both will work). After you import the coordinates you can run the "Generate Output - e2boxer" task to write particle image data."""
	def __init__(self):
		WorkFlowTask.__init__(self)
		self.window_title = "Import Particle Coordinates"
		
	def get_params(self):
		params = []
		from emform import EMPlotTable
		from emsave import EMCoordFileValidator
		coord_table = EMPlotTable(name="coordfiles",desc_short="Coord File",desc_long="Use this tool to browse for your coordinate files")
		context_menu_data = ParticleWorkFlowTask.DataContextMenu(EMCoordFileValidator())
		coord_table.add_context_menu_data(context_menu_data)
		coord_table.add_button_data(ParticleWorkFlowTask.AddDataButton(coord_table,context_menu_data))
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=self.__doc__,choices=None))
		params.append(coord_table)  
		
		return params

	
	def on_form_ok(self,params):
		if not params.has_key("coordfiles") or len(params["coordfiles"]) == 0:
			error("Please choose files to import")
			return
		
#		from e2boxer import merge_boxes_as_manual_to_db
#		merge_boxes_as_manual_to_db(params["coordfiles"])
		pdb=db_open_dict("bdb:.#project")
		img_list=pdb.get("global.spr_raw_data_dict",dfl={}).keys()
		if len(img_list)==0 :
			print "No image files in project !!!"
			return
		bdb=db_open_dict("bdb:e2boxercache#boxes")
		for ff in params["coordfiles"]:
			f=ff.rsplit("/",1)[-1]
			fsp=best_match(f,img_list)
			lns=file(ff,"r").readlines()
			lns=[l.split() for l in lns]
			bxs=[[int(l[0])+int(l[2])/2,int(l[1])+int(l[3])/2,"manual"] for l in lns]
			bdb[fsp]=bxs
			print "imported ",f," into ",fsp
								
		self.disconnect_form()
		self.form.close()
		self.form = None
		self.emit(QtCore.SIGNAL("task_idle"))
		
from os.path import commonprefix
def best_match(name,lst):
	ml=[(len(commonprefix((name,i.split('/')[-1]))),n) for n,i in enumerate(lst)]
#	print name,lst,ml
	if max(ml)[0]==0 : return None
	return lst[max(ml)[1]]
		
class E2BoxerTask(ParticleWorkFlowTask):
	'''
	Provides some common functions for the e2boxer tasks
	'''
	def __init__(self):
		ParticleWorkFlowTask.__init__(self)
		self.form_db_name = "bdb:emform.e2boxer"
		self.report_task = None  #will eventually store a EMRawDataReportTask
		
		data_dict = EMProjectDataDict(spr_ptcls_dict)
		dict = data_dict.get_data_dict() # this is to protect against back compatibility problems
	
	
	def get_boxes_in_database(file_name):
		
#		print "GDB",file_name
		db_name = "bdb:e2boxercache#boxes"
		#db_name = "bdb:e2boxer.cache"		
		box_maps = {}
		nbox = 0
		if db_check_dict(db_name):
			e2boxer_db = db_open_dict(db_name,ro=True)
			if e2boxer_db.has_key(file_name):
				return str(len(e2boxer_db[file_name]))
		return "-"

	get_boxes_in_database = staticmethod(get_boxes_in_database)
#	get_num_particles_project = staticmethod(get_num_particles_project)
#	get_particle_dims_project = staticmethod(get_particle_dims_project)
	
	class ParticleColumns:
		def __init__(self,project_dict=spr_ptcls_dict):
			self.header_cache = {}
			self.translation_cache = {}
			self.project_dict = project_dict
			data_dict = EMProjectDataDict(self.project_dict)
			particle_names = data_dict.get_names()
			
			# cache everything at the start, it's the most efficient approach, in this strategy
			for name in particle_names:
				a = EMData()
				try: a.read_image(name,0,True) # header only
				except: a=EMData(1,1)
				d = a.get_attr_dict()
				if d.has_key("ptcl_source_image"):
					file_name = d["ptcl_source_image"]
					self.header_cache[name] = d
					self.translation_cache[file_name] = name
			
		def get_num_particles_project(self,file_name):
			'''
			Get the particles in the project that are associated with a specific file name
			This is useful for the e2boxer forms, which are used to take raw data files and
			produce boxed output - i.e. if the user wants to know if they've already
			written boxed output for a given raw file
			@param file_name a file name, should be a file that's in global.spr_raw_data_dict - this is not checked though
			Note that the only thing that defines the relationship is whether or not the particle's 
			'''
			if self.translation_cache.has_key(file_name):
				name = self.translation_cache[file_name]
				d = self.header_cache[name]
				if d.has_key("ptcl_source_image"):
					if d["ptcl_source_image"] == file_name:
						return str(EMUtil.get_image_count(name))
			
			
			data_dict = EMProjectDataDict(self.project_dict)
			particle_names = data_dict.get_names()
			if len(particle_names) == 0: return "-"
			

			for name in particle_names:
				a = EMData(1,1)
				try : a.read_image(name,0,True) # header only
				except: pass
				d = a.get_attr_dict()
				if d.has_key("ptcl_source_image"):
					if d["ptcl_source_image"] == file_name:
						self.header_cache[name] = d
						self.translation_cache[file_name] = name
						return str(EMUtil.get_image_count(name))
					
			return "0"
		
		def get_particle_dims_project(self,file_name):
			'''
			'''
			if self.translation_cache.has_key(file_name):
				name = self.translation_cache[file_name]
				d = self.header_cache[name]
				if d.has_key("ptcl_source_image"):
					if d["ptcl_source_image"] == file_name:
						nx,ny,nz = gimme_image_dimensions3D(name)
						return "%ix%ix%i" %(nx,ny,nz)
			
			data_dict = EMProjectDataDict(self.project_dict)
			particle_names = data_dict.get_names()
			if len(particle_names) == 0: return "-"
			
			for name in particle_names:
				a = EMData(1,1)
				try: a.read_image(name,0,True) # header only
				except: pass
				d = a.get_attr_dict()
				#print d
				if d.has_key("ptcl_source_image"):
					if d["ptcl_source_image"] == file_name:
						nx,ny,nz = gimme_image_dimensions3D(name)
						self.header_cache[name] = d
						self.translation_cache[file_name] = name
						return "%ix%ix%i" %(nx,ny,nz)
			return ""
			
	
	def get_boxer_basic_table(self):
		'''
		
		Returns a table like this:
		
		|| Project image name || Boxes in e2boxer db ||
		
		Returns the table, and the the number of entries (p,n)
		if n is zero there are no entries in the table and the calling function can act appropriately
		'''
		
		self.report_task = EMRawDataReportTask()
		table,n = self.report_task.get_raw_data_table()
		from emform import EMFileTable,int_lt
		table.insert_column_data(0,EMFileTable.EMColumnData("Stored Boxes",E2BoxerTask.get_boxes_in_database,"Boxes currently stored in the EMAN2 database",int_lt))
		table.add_column_data(EMFileTable.EMColumnData("Quality",E2BoxerTask.get_quality,"Quality metadata score stored in local database",int_lt))
	
		return table, n
	
	def get_quality(file_name):
		'''
		A static function for getting the number of boxes associated with each file
		'''
		from emboxerbase import get_database_entry
		val = get_database_entry(file_name,"quality")
		
		if val == None: return "-"
		else: return str(val)
		
	get_quality = staticmethod(get_quality)
		
	
	def get_project_files_that_have_db_boxes_in_table(self):
		
		self.report_task = EMRawDataReportTask()
		table,n = self.report_task.get_raw_data_table()
		
		data_dict = EMProjectDataDict(spr_ptcls_dict)
		dict = data_dict.get_data_dict() # this is to protect against back compatibility problems. This is necessary for the columns_object to operate without throwing (in unusual circumstances the user deletes the particles, and this accomodates for it)
	
		from emform import EMFileTable,int_lt
		table.insert_column_data(0,EMFileTable.EMColumnData("Stored Boxes",E2BoxerTask.get_boxes_in_database,"Boxes currently stored in the EMAN2 database",int_lt))
		table.insert_column_data(1,EMFileTable.EMColumnData("Quality",E2BoxerTask.get_quality,"Quality metadata score stored in local database",int_lt))

		self.columns_object = E2BoxerTask.ParticleColumns()
		table.insert_column_data(2,EMFileTable.EMColumnData("Particles On Disk",self.columns_object.get_num_particles_project,"Particles currently stored on disk that are associated with this image",int_lt))
		table.insert_column_data(3,EMFileTable.EMColumnData("Particle Dims",self.columns_object.get_particle_dims_project,"The dimensions of the particles that are stored on disk"))
		#self.tmp = E2BoxerTask.Tmp()
		return table, n
		

#		
	def __get_e2boxer_data(self,project_names):
		
		db_name = "bdb:e2boxer.cache"		
		box_maps = {}
		if db_check_dict(db_name):
			e2boxer_db = db_open_dict(db_name,ro=True)
			for name in e2boxer_db.keys():
				d = e2boxer_db[name]
				if not isinstance(d,dict): continue
				if not d.has_key("e2boxer_image_name"): # this is the test, if something else has this key then we're screwed.
					continue
				name = d["e2boxer_image_name"]
				if not name in project_names: continue
				dim = ""
				nbox = 0
				for key in ["auto_boxes","manual_boxes","reference_boxes"]:
					if d.has_key(key):
						boxes = d[key]
						if boxes != None:
							nbox += len(boxes)
							if dim == "" and len(boxes) > 0:
								box = boxes[0]
								dim = str(box.xsize) + "x"+str(box.ysize)
				
				if nbox == 0: nbox = "" # just so it appears as nothin in the interface			
				box_maps[name] = [dim,nbox]
		
		nboxes = []
		dimensions = []
		for name in project_names:
			if box_maps.has_key(name):
				dimensions.append(box_maps[name][0])
				nboxes.append(box_maps[name][1])
			else:
				dimensions.append("")
				nboxes.append("")
				
		return nboxes,dimensions
	
	def on_form_cancel(self):
		if self.report_task != None:
			self.report_task.recover_original_raw_data_list()
		self.form.close()
		self.form = None
		self.emit(QtCore.SIGNAL("task_idle"))


class E2BoxerGenericTask(ParticleWorkFlowTask):
	documentation_string = "Fill me in"
	def __init__(self):
		ParticleWorkFlowTask.__init__(self)
		self.window_title = "e2boxer"
		self.preferred_size = (480,200)
		self.form_db_name = "bdb:emform.e2boxer"
		
	def get_params(self):
		params = []		
#		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2CTFGenericTask.documentation_string,choices=None))
		db = db_open_dict(self.form_db_name)
		params.append(ParamDef(name="running_mode",vartype="choice",desc_short="Choose your running mode",desc_long="There are three Boxer related task which are generally run in order",property=None,defaultunits=db.get("running_mode",dfl="Interactive boxing"),choices=["Interactive boxing", "Autoboxing", "Write output"]))
		#db_close_dict(self.form_db_name)
		return params

	def on_form_ok(self,params):
		if params["running_mode"] == "Interactive boxing":
			self.emit(QtCore.SIGNAL("replace_task"),E2BoxerGuiTaskGeneral(),"e2boxer interface launcher")
			self.form.close()
			self.form = None
		elif params["running_mode"] == "Autoboxing":
			self.emit(QtCore.SIGNAL("replace_task"),E2BoxerAutoTaskGeneral(),"e2boxer automated boxing")
			self.form.close()
			self.form = None
		elif params["running_mode"] == "Write output":
			self.emit(QtCore.SIGNAL("replace_task"),E2BoxerOutputTaskGeneral(),"e2boxer write output")
			self.form.close()
			self.form = None	
		else:
			self.form.close()
			self.form = None
			self.emit(QtCore.SIGNAL("task_idle"))
			return
		
		self.write_db_entries(params)
	def write_db_entry(self,key,value):
		pass


class E2BoxerAutoTask(E2BoxerTask):
	'''
	DO NOT USE THIS !
	'''
	
	def __init__(self):
		E2BoxerTask.__init__(self)
		self.window_title = "e2boxer Autoboxing"
		self.boxer_module = None # this will actually point to an EMBoxerModule, potentially
		recover_old_boxer_database()
	def get_params(self):
		params = []
		
		params.append(ParamDef(name="blurb",vartype="text",desc_short="Using e2boxer",desc_long="",property=None,defaultunits=self.__doc__,choices=None))
		
		db_name = "bdb:e2boxercache#swarm_boxers"
		if db_check_dict(db_name):
			db = db_open_dict(db_name)
			project_list = spr_raw_data_dict
		
			data_dict = EMProjectDataDict(project_list)
			self.project_data_at_init = data_dict.get_data_dict() # so if the user hits cancel this can be reset
			project_names = data_dict.keys()
			
			culled_names = []
			for name in project_names:
				if db.has_key(name) and len(db[name]) > 0:
					culled_names.append(name)
			
			from emform import EM2DFileTable,EMFileTable
			table = EM2DFileTable(culled_names,desc_short="Choose Autoboxer",desc_long="",single_selection=True,name="autoboxer")
			table.add_column_data(EMFileTable.EMColumnData("Autoboxer ID",E2BoxerAutoTask.auto_boxer_id,"Autoboxer time stamp indicating last alteration"))
			table.add_column_data(EMFileTable.EMColumnData("Mode",E2BoxerAutoTask.swarm_mode,"Swarm picking mode"))
			table.add_column_data(EMFileTable.EMColumnData("Proximity Threshold",E2BoxerAutoTask.proximity_threshold,"Swarm autoboxer proximity threshold"))
			table.add_column_data(EMFileTable.EMColumnData("Threshold",E2BoxerAutoTask.swarm_threshold,"Swarm correlation threshold"))
			table.add_column_data(EMFileTable.EMColumnData("References",E2BoxerAutoTask.swarm_num_refs,"The number of references used by the Swarm autoboxer"))
			#table.add_column_data(EMFileTable.EMColumnData("Profile",E2BoxerAutoTask.swarm_profile,"Swarm picking profile"))
			
			data_table,n = self.get_boxer_basic_table() # note n is unused, it's a refactoring residual		
			
			params.append([table,data_table])
			
		return params
	
	def auto_boxer_id(filename):
		db_name = "bdb:e2boxercache#swarm_boxers"
		if db_check_dict(db_name):
			db = db_open_dict(db_name)
			data = db[filename]
			return str(data[-1][9])
		return "-"
	auto_boxer_id = staticmethod(auto_boxer_id)
		
	def proximity_threshold(filename):
		db_name = "bdb:e2boxercache#swarm_boxers"
		if db_check_dict(db_name):
			db = db_open_dict(db_name)
			data = db[filename]
			return str(data[-1][10])
		return "-"	
	proximity_threshold = staticmethod(proximity_threshold)
	
	def swarm_mode(filename):
		db_name = "bdb:e2boxercache#swarm_boxers"
		if db_check_dict(db_name):
			db = db_open_dict(db_name)
			data = db[filename]
			return str(data[-1][7])
		return "-"
	
	swarm_mode = staticmethod(swarm_mode)
	
	def swarm_threshold(filename):
		db_name = "bdb:e2boxercache#swarm_boxers"
		if db_check_dict(db_name):
			db = db_open_dict(db_name)
			data = db[filename]
			return str(data[-1][3])
		return "-"
		
	swarm_threshold = staticmethod(swarm_threshold)
	
	
	def swarm_profile(filename):
		db_name = "bdb:e2boxercache#swarm_boxers"
		if db_check_dict(db_name):
			db = db_open_dict(db_name)
			data = db[filename]
			return str(data[-1][2])
		return "-"
		
	swarm_profile = staticmethod(swarm_profile)
	
	def swarm_num_refs(filename):
		db_name = "bdb:e2boxercache#swarm_boxers"
		if db_check_dict(db_name):
			db = db_open_dict(db_name)
			data = db[filename]
			return str(len(data[-1][0]))
		return "-"
	swarm_num_refs = staticmethod(swarm_num_refs)
			
	def on_form_ok(self,params): 
		if  params.has_key("filenames") and len(params["filenames"]) == 0:
			self.run_select_files_msg()
			return
		
		if  params.has_key("autoboxer") and len(params["autoboxer"]) == 0:
			error("Please choose autoboxer data")
			return
	
		else:
			self.write_db_entries(params) # will only write filenames
			options = EmptyObject()
			for k,v in params.items():
				setattr(options,k,v)
			
			options.autoboxer = params["autoboxer"][0]
			string_args = ["autoboxer"]
			additional_args = []
			bool_args = []
			temp_file_name = "e2boxer_autobox_stdout.txt"
			self.spawn_task("e2boxer.py",options,string_args,bool_args,additional_args,temp_file_name)
			self.emit(QtCore.SIGNAL("task_idle"))
			self.form.close()
			self.form = None

class E2BoxerAutoTaskGeneral(E2BoxerAutoTask):
	def __init__(self):
		E2BoxerAutoTask.__init__(self)
		self.window_title = "e2boxer Autoboxing"
		self.boxer_module = None # this will actually point to an EMBoxerModule, potentially

	def get_params(self):
		params = []
		
		
		db_name = "bdb:e2boxer.cache"
		fail = False
		if db_check_dict(db_name):
			db = db_open_dict(db_name,ro=True)
			if not db.has_key("current_autoboxer"): fail = True
		else:
			fail = True
			
		if fail:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="Using e2boxer",desc_long="",property=None,defaultunits=E2BoxerAutoTask.documentation_string+E2BoxerAutoTask.warning_string,choices=None))
		else:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="Using e2boxer",desc_long="",property=None,defaultunits=E2BoxerAutoTask.documentation_string,choices=None))
			params.append(ParamDef(name="filenames",vartype="url",desc_short="File Names",desc_long="The names of the particle files you want to interactively box using e2boxer",property=None,defaultunits=[],choices=[]))
		
		return params


class OldBoxerRecoveryDialog(QtGui.QDialog):
	def __init__(self):
		'''
		@param sym some kind of symmetry, such as "d7", "icos" etc
		'''
		QtGui.QDialog.__init__(self)
		self.setWindowTitle("Old Boxer Recovery")
		self.setWindowIcon(QtGui.QIcon(get_image_directory() + "green_boxes.png"))

		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		
		text_edit = QtGui.QTextEdit("")
		text_edit.setReadOnly(True)
		text_edit.setWordWrapMode(QtGui.QTextOption.WordWrap)
		text_edit.setText("The workflow has detected you have data stored in the local database that was generated with an old version of e2boxer. You can recover it (recommended), in which case the old data is converted so it can be interpreted within the current framework. Alternatively you can just delete it, which means the box coordinates will be lost forever.")
		self.vbl.addWidget(text_edit)
		self.button_hbl = QtGui.QHBoxLayout()
		self.recover = QtGui.QPushButton("Recover")
		self.recover.setToolTip("The old database will be converted to a format recognized by the new boxer. The old database will then be deleted.")
		self.recover.setDefault(True)
		self.remove = QtGui.QPushButton("Remove")
		self.remove.setToolTip("The old database will be removed and all previous boxing results will be deleted from disk.")
		self.cancel = QtGui.QPushButton("Cancel")
		self.cancel.setToolTip("The operation will be cancelled.")
		self.button_hbl.addWidget(self.cancel )
		self.button_hbl.addWidget(self.remove )
		self.button_hbl.addWidget(self.recover )
		self.vbl.addLayout(self.button_hbl)

	
		QtCore.QObject.connect(self.recover, QtCore.SIGNAL("clicked(bool)"), self.on_recover)
		QtCore.QObject.connect(self.remove, QtCore.SIGNAL("clicked(bool)"), self.on_remove)
		QtCore.QObject.connect(self.cancel, QtCore.SIGNAL("clicked(bool)"), self.on_cancel)
		self.ret_code = 0
	def on_cancel(self,int):
		self.ret_code = 0
		self.accept()
		
	def on_remove(self,int):
		self.ret_code = 1
		self.accept()
		
	def on_recover(self,int):
		self.ret_code = 2
		self.accept()
	
	def exec_(self):
		'''
		Customized exec_ function
		@return None if the user hit cancel or a dictionary containing important parameters if the user hit ok
		'''
		QtGui.QDialog.exec_(self)
		return self.ret_code

def recover_old_boxer_database():
	'''
	This function will become redundant, probably by the beginning of 2010
	'''
	old_boxer_database = "bdb:e2boxer.cache"
	if db_check_dict(old_boxer_database):
		recovery_items = []
		db = db_open_dict(old_boxer_database)
		for key,value in db.items():
			if isinstance(key,str) and len(key) > 2 and key[-3:] == "_DD":
				if isinstance(value,dict):
					if value.has_key("reference_boxes") or value.has_key("auto_boxes") or value.has_key("manual_boxes"):
						if value.has_key("e2boxer_image_name"):
							recovery_items.append([key,value])
		if len(recovery_items) > 0:
			dialog = OldBoxerRecoveryDialog()
			code = dialog.exec_()
			if code == 0:
				return 0
			else:
				if code == 1:
					db_remove_dict(old_boxer_database)
				else: # code == 2
					new_db_name = "bdb:e2boxercache#boxes"
					db = db_open_dict(new_db_name)
					for item in recovery_items:
						from pyemtbx.boxertools import TrimBox
						value = item[1]
						name = value["e2boxer_image_name"]
						auto = None
						ref = None
						man = None
						new_boxes = []
						for box_name in ["reference_boxes","auto_boxes","manual_boxes"]:
							if value.has_key(box_name):
								old_boxes =  value[box_name]
								for box in old_boxes:
									x = box.xcorner+box.xsize/2
									y = box.ycorner+box.ysize/2
	#									score = box.correlation_score
									new_boxes.append([x,y,"manual"])
						
						print name,len(new_boxes)
						if db.has_key(name):
							new_boxes = new_boxes.extend(db[name])
						
						db[name] = new_boxes
					
					db_remove_dict(old_boxer_database)
			return 1
		
	return 1
			

class E2BoxerGuiTask(E2BoxerTask):	
	documentation_string = """Select the frames you want to select boxes from, enter your boxsize, and hit OK.
NOTE - SELECTING A GOOD BOX SIZE IS CRITICAL. See the wiki for a list of good sizes. Make sure the box is ~2x the size of your particle. \
Changing box size later can be very painful, and CTF correction relies on a sufficiently large box (larger than used with EMAN1.

This will lauch e2boxer and automatically load the selected images for boxing. \
Generally you don't want to work with more than ~10 at a time. To autobox, make sure you select SWARM mode before selecting boxes. """
	
	warning_string = "\n\n\nNOTE: There are no images currenty associated with the project. Please import or specify which images you want as part of this project in step 1 of the workflow and try again."
	
	def __init__(self):
		E2BoxerTask.__init__(self)
		self.window_title = "Launch e2boxer Interface"
		self.boxer_module = None # this will actually point to an EMBoxerModule, potentially
		recover_old_boxer_database()
	def get_params(self):
		params = []
		
		p,n = self.get_boxer_basic_table() # note n is unused, it's a refactoring residual		
		params.append(ParamDef(name="blurb",vartype="text",desc_short="Interactive use of e2boxer",desc_long="",property=None,defaultunits=E2BoxerGuiTask.documentation_string,choices=None))
		params.append(p)
		db = db_open_dict(self.form_db_name)
		params.append(ParamDef(name="interface_boxsize",vartype="int",desc_short="Box size",desc_long="An integer value",property=None,defaultunits=db.get("interface_boxsize",dfl=128),choices=[]))
			#db_close_dict(self.form_db_name)
		return params
			
	def on_form_ok(self,params):
		
		if not params.has_key("filenames"): return
		
		if  params.has_key("filenames") and len(params["filenames"]) == 0:
			self.run_select_files_msg()
			return

		if  params.has_key("interface_boxsize") and params["interface_boxsize"] < 1:
			self.show_error_message(["Must specify a positive, non zero boxsize."])
			return
		else:
			self.write_db_entries(params)
			options = EmptyObject()
			for key in params.keys():
				setattr(options,key,params[key])
			options.boxsize = params["interface_boxsize"]
			options.running_mode = "gui"
			options.method = "Swarm"
			
			from emboxerbase import EMBoxerModule
			from e2boxer import  SwarmTool
			self.boxer_module = EMBoxerModule(params["filenames"],params["interface_boxsize"])
															
			
			# this is an example of how to add your own custom tools:
			self.boxer_module.add_tool(SwarmTool,particle_diameter=options.boxsize)
			
			
#			from e2boxer import EMBoxerModule
#			self.boxer_module = EMBoxerModule(get_application(),options)
			self.emit(QtCore.SIGNAL("gui_running"),"Boxer",self.boxer_module) # The controlled program should intercept this signal and keep the E2BoxerTask instance in memory, else signals emitted internally in boxer won't work
			
			QtCore.QObject.connect(self.boxer_module, QtCore.SIGNAL("module_closed"), self.on_boxer_closed)
			self.form.close()
#			self.boxer_module.show_guis()
			self.boxer_module.show_interfaces()
			self.form = None
			
	def on_form_close(self):
		# this is to avoid a task_idle signal, which would be incorrect if e2boxer is running
		if self.boxer_module == None:
			self.emit(QtCore.SIGNAL("task_idle"))
		else: pass
	
	def on_boxer_closed(self):
		if self.boxer_module != None:
			self.boxer_module = None
			self.emit(QtCore.SIGNAL("gui_exit"))
	

class E2BoxerGuiTaskGeneral(E2BoxerGuiTask):	
	def __init__(self):
		E2BoxerTask.__init__(self)
		self.window_title = "Launch e2boxer Interface"
		self.boxer_module = None # this will actually point to an EMBoxerModule, potentially

	def get_params(self):
		params = []
		params.append(ParamDef(name="blurb",vartype="text",desc_short="Using e2boxer",desc_long="",property=None,defaultunits=E2BoxerGuiTask.documentation_string,choices=None))
		params.append(ParamDef(name="filenames",vartype="url",desc_short="File Names",desc_long="The names of the particle files you want to interactively box using e2boxer",property=None,defaultunits=[],choices=[]))
		db = db_open_dict(self.form_db_name)
		params.append(ParamDef(name="interface_boxsize",vartype="int",desc_short="Box size",desc_long="An integer value",property=None,defaultunits=db.get("interface_boxsize",dfl=128),choices=[]))
		#db_close_dict(self.form_db_name)
		return params
	
	
class E2BoxerOutputTask(E2BoxerTask):	
	"""This task will write the selected particles to output files. The default format is 'BDB', which should be used if you plan to continue processing using the workflow interface. HDF is the only other format which will preserve full metadata. img (IMAGIC) and spi (SPIDER) will lose metadata if used.
	
	Select the images you wish to generate output for, enter the box size and normalization etc, and then hit OK."""
	def __init__(self):
		E2BoxerTask.__init__(self)
		self.window_title = "Generate e2boxer Output"
		self.output_formats = ["bdb","hdf","img","spi"] # disable img from the workflow because in EMAN2 we want to store more metadata in the header
		recover_old_boxer_database()
#	def __del__(self):
#		print "output task dies"
	
	def get_params(self):
		params = []
		
		p,n = self.get_project_files_that_have_db_boxes_in_table()

		params.append(ParamDef(name="blurb",vartype="text",desc_short="Using e2boxer",desc_long="",property=None,defaultunits=self.__doc__,choices=None))
		params.append(p)
		self.add_general_params(params)

		return params
	
	def add_general_params(self,params):
		'''
		Functionality used in several places
		'''
		db = db_open_dict(self.form_db_name)
		pbox = ParamDef(name="output_boxsize",vartype="int",desc_short="Box size",desc_long="An integer value",property=None,defaultunits=db.get("output_boxsize",dfl=128),choices=[])	
		pfo = ParamDef(name="force",vartype="boolean",desc_short="Force overwrite",desc_long="Whether or not to force overwrite files that already exist",property=None,defaultunits=db.get("force",dfl=False),choices=None)
		pwc = ParamDef(name="write_dbbox",vartype="boolean",desc_short="Write box coord files",desc_long="Whether or not box db files should be written",property=None,defaultunits=db.get("write_dbbox",dfl=False),choices=None)
		pwb = ParamDef(name="write_ptcls",vartype="boolean",desc_short="Write box image files",desc_long="Whether or not box images should be written",property=None,defaultunits=db.get("write_ptcls",dfl=True),choices=None)
		pinv = ParamDef(name="invert",vartype="boolean",desc_short="Invert",desc_long="Do you want the pixel intensities in the output inverted?",property=None,defaultunits=db.get("invert",dfl=False),choices=None)
		pn =  ParamDef(name="norm",vartype="string",desc_short="Normalize images",desc_long="How the output box images should be normalized",property=None,defaultunits=db.get("norm",dfl="normalize.edgemean"),choices=["normalize","normalize.edgemean","none"])
		pop = ParamDef(name="format",vartype="string",desc_short="Output image format",desc_long="The format of the output box images",property=None,defaultunits=db.get("format",dfl="bdb"),choices=self.output_formats)
		
		#db_close_dict(self.form_db_name)
		pwb.dependents = ["invert","norm","format"] # these are things that become disabled when the pwb checkbox is unchecked etc
		
		params.append([pbox,pfo])
		params.append([pwc,pwb,pinv])
		params.append(pn)
		params.append(pop)
		
	
	def check_params(self,params):
		
		error_message = []
		if params["output_boxsize"] < 1: error_message.append("Boxsize must be greater than 0.")
		if not params["write_ptcls"] and not params["write_dbbox"]: error_message.append("You must choose at least one of the write_coords/write_box_images options")
	
		return error_message
	
	def on_form_ok(self,params):	
		if  params.has_key("filenames") and len(params["filenames"]) == 0:
			self.run_select_files_msg()
			return
		
		error_message = self.check_params(params)
		if len(error_message) >0: 
			self.show_error_message(error_message)
			return
		
		else:
			self.write_db_entries(params)
			options = EmptyObject()
			for k,v in params.items():
				setattr(options,k,v)	
			options.boxsize = params["output_boxsize"]
			
			options.just_output=True # this is implicit, it has to happen
			
			string_args = ["norm","format","boxsize"]
			bool_args = ["force","write_dbbox","write_ptcls","invert"]
			additional_args = ["--dbls=%s" %spr_ptcls_dict]
			temp_file_name = "e2boxer_autobox_stdout.txt"
			self.spawn_single_task("e2boxer.py",options,string_args,bool_args,additional_args,temp_file_name)
			self.emit(QtCore.SIGNAL("task_idle"))
			self.form.close()
			self.form = None

class E2BoxerOutputTaskGeneral(E2BoxerOutputTask):
	documentation_string = "Write me"
	def __init__(self):
		E2BoxerOutputTask.__init__(self)
		self.window_title = "Generate e2boxer Output"
		
	def get_params(self):
		params = []
		params.append(ParamDef(name="blurb",vartype="text",desc_short="Using e2boxer",desc_long="",property=None,defaultunits=E2BoxerOutputTaskGeneral.documentation_string,choices=None))
		
		p = self.get_e2boxer_boxes_table(project_check=False)
		params.append(p)
		
		self.add_general_params(params)
	
#		boxer_project_db = db_open_dict("bdb:e2boxer.project")
#		params.append(ParamDef(name="boxsize",vartype="int",desc_short="Box size",desc_long="An integer value",property=None,defaultunits=boxer_project_db.get("interface_boxsize",dfl=128),choices=[]))
		return params
	
	def get_e2boxer_boxes_table(self,project_check=True):
		db_name = "bdb:e2boxer.cache"
		p = EMParamTable(name="filenames",desc_short="Current boxes generated by e2boxer",desc_long="")
		names = []
		nboxes = []
		dimensions = []
		
		if project_check:
			project_db = db_open_dict("bdb:project")
			project_data = project_db.get(spr_raw_data_dict,dfl={})
			project_names = project_data.keys()
		
		if db_check_dict(db_name):
			e2boxer_db = db_open_dict(db_name,ro=True)
			for name in e2boxer_db.keys():
				d = e2boxer_db[name]
				if not isinstance(d,dict): continue
				if not d.has_key("e2boxer_image_name"): # this is the test, if something else has this key then we're screwed.
					continue

				name = d["e2boxer_image_name"]
				if project_check:
					if not name in project_names: continue
				names.append(name)
				
				dim = ""
				nbox = 0
				for key in ["auto_boxes","manual_boxes","reference_boxes"]:
					if d.has_key(key):
						boxes = d[key]
						nbox += len(boxes)
						if dim == "" and len(boxes) > 0:
							box = boxes[0]
							dim = str(box.xsize) + "x"+str(box.ysize)
							
							
				nboxes.append(nbox)
				dimensions.append(dim)
			
		pnames = ParamDef(name="Filenames",vartype="stringlist",desc_short="File Names",desc_long="The filenames",property=None,defaultunits=None,choices=names)
		pboxes = ParamDef(name="Num boxes",vartype="intlist",desc_short="Boxes in DB",desc_long="The number of boxes stored for this image in the database",property=None,defaultunits=None,choices=nboxes)
		pdims = ParamDef(name="Dimensions",vartype="stringlist",desc_short="Dimensions",desc_long="The dimensions boxes",property=None,defaultunits=None,choices=dimensions)
		
		p = EMParamTable(name="filenames",desc_short="Choose a subset of these images",desc_long="")
		p.append(pnames)
		p.append(pboxes)
		p.append(pdims)
		return p
	
class E2BoxerProgramOutputTask(E2BoxerOutputTask):
	'''
	This task is called from e2boxer itself. Not from the workflow
	'''
	documentation_string = "Use this form to write output file from within the e2boxer interface.\nYou can choose to write image files in a number of formats. The bdb file format is most useful if you are using EMAN2. If you plan to use your data with other programs, including EMAN1, you must choose either the hdf or img output formats.\nYou can also choose to write EMAN1 style .box files"
	def __init__(self,application,filenames,target,exclusions=[]):
		E2BoxerOutputTask.__init__(self)
		self.window_title = "Generate e2boxer Output"
		self.filenames = filenames
		self.target = weakref.ref(target)
		self.exclusions = exclusions
		self.output_formats = ["bdb","img","hdf"]
		
	def get_params(self):
		params = []
		params.append(ParamDef(name="blurb",vartype="text",desc_short="E2Boxer output form",desc_long="",property=None,defaultunits=E2BoxerProgramOutputTask.documentation_string,choices=None))
		
		p = EMParamTable(name="filenames",desc_short="Choose a subset of these images",desc_long="")
		pnames = ParamDef(name="Filenames",vartype="stringlist",desc_short="File Names",desc_long="The filenames",property=None,defaultunits=None,choices=self.filenames)
		p.append(pnames)
		setattr(p,"convert_text", ptable_convert_2)
		setattr(p,"icon_type","single_image")
		setattr(p,"exclusions",self.exclusions)
		
		params.append(p)
		
		self.add_general_params(params)
		return params
	
	def on_form_ok(self,params):

		if  params.has_key("filenames") and len(params["filenames"]) == 0:
			self.run_select_files_msg()
			return
		
		error_message = self.check_params(params)
		if len(error_message) > 0: 
			self.show_error_message(error_message)
			return
		else:
			if params["write_dbbox"]:
				self.target().write_coord_files(params["filenames"],params["output_boxsize"],params["force"])
			if params["write_ptcls"]:
				normproc = False
				if params["norm"] != "none":
					normproc=True
				self.target().write_box_image_files(params["filenames"],params["output_boxsize"],params["force"],params["format"],normproc,params["norm"],params["invert"])
				
			self.emit(QtCore.SIGNAL("task_idle"))
			self.form.close()
			self.form = None

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
		
		from emform import EMFileTable,float_lt,int_lt
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
		from emform import EMFileTable,int_lt
		table.add_column_data(EMFileTable.EMColumnData("Phase flip",self.other_column_data.get_num_phase_flipped,"The number of phase flipped particles on disk",int_lt))
		table.add_column_data(EMFileTable.EMColumnData("Phase flip dims",self.other_column_data.phase_flipped_dim,"The dimensions of the phase flippped particles"))
		table.add_column_data(EMFileTable.EMColumnData("Phase flip hp",self.other_column_data.get_num_phase_flipped_hp,"The number of phase flipped high pass filtered particles on disk",int_lt))
#		table.add_column_data(EMFileTable.EMColumnData("Phase flip hp dims",self.other_column_data.phase_flipped_dim,"The dimensions of the phase flippped particles"))
		table.add_column_data(EMFileTable.EMColumnData("Wiener filt",self.other_column_data.get_num_wein_filt,"The number of Wiener filtered particles on disk",int_lt))
#		table.add_column_data(EMFileTable.EMColumnData("Wiener filt dims",self.other_column_data.wien_filt_dim,"The dimensions of the Wiener filtered particles"))
		return table, n
	
	class MoreCTFColumns:
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
			if self.db_map.has_key(tag):
				val = self.db_map[tag]
				if val.has_key(filt):
					file_name = val[filt]
					return str(EMUtil.get_image_count(file_name))
				
			return ""
		
		def __get_dim_filtered(self,name,filt):
			tag = base_name(name)
			if self.db_map.has_key(tag):
				val = self.db_map[tag]
				if val.has_key(filt):
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
		for key,data in parms_db.items():
			if data == None:
				print "error?",key
				continue
			ret.append([key,data[-1]]) # parms[-1] should be the original filename
		#db_close_dict("bdb:e2ctf.parms")
		return ret

class CTFReportTask(E2CTFWorkFlowTask):
	
	documentation_string = "This tool is for displaying the currently determined CTF parameters for the particles associated with the project. It also displays \
the number of phase flipped and/or wiener filtered images corresponding to each particle set. Normally you don't do anything on this page other than get a \
status report on current processing."
	warning_string = "\n\n\nNOTE: There are no particles currently associated with the project. Please go to the \"Particles\" task and import/box particles first."
	def __init__(self):
		E2CTFWorkFlowTask.__init__(self)
		self.window_title = "CTF Parameters And Images"

	def get_params(self):
		params = []
		p,n = self.get_full_ctf_table()

		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=CTFReportTask.documentation_string,choices=None))
		params.append(p)
		return params

	def write_db_entry(self,key,value):
		pass		

class E2CTFGenericTask(ParticleWorkFlowTask):
	documentation_string = "Fill me in"
	def __init__(self):
		ParticleWorkFlowTask.__init__(self)
		self.window_title = "Generic e2ctf Task"
		self.preferred_size = (480,200)
		self.form_db_name = "bdb:emform.e2ctf"
		
		
	def get_params(self):
		params = []		
#		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2CTFGenericTask.documentation_string,choices=None))
		db = db_open_dict(self.form_db_name)
		params.append(ParamDef(name="running_mode",vartype="choice",desc_short="Choose your running mode",desc_long="There are three CTF related task which are generally run in order",property=None,defaultunits=db.get("running_mode",dfl="auto params"),choices=["auto params", "interactively fine tune", "write output"]))
		#db_close_dict(self.form_db_name)
		return params

	def on_form_ok(self,params):
		if params["running_mode"] == "auto params":
			self.emit(QtCore.SIGNAL("replace_task"),E2CTFAutoFitTaskGeneral(),"ctf auto fit")
			self.form.close()
			self.form = None
		elif params["running_mode"] == "interactively fine tune":
			self.emit(QtCore.SIGNAL("replace_task"),E2CTFGuiTaskGeneral(),"fine tune ctf")
			self.form.close()
			self.form = None
		elif params["running_mode"] == "write output":
			self.emit(QtCore.SIGNAL("replace_task"),E2CTFOutputTaskGeneral(),"ctf output")
			self.form.close()
			self.form = None	
		else:
			self.form.close()
			self.form = None
			self.emit(QtCore.SIGNAL("task_idle"))
			return
		
		self.write_db_entries(params)
			
class E2CTFAutoFitTask(E2CTFWorkFlowTask):	
	documentation_string = """Select the particles you wish to generate CTF parameters for, enter the appropriate parameters such as microscope voltage etc, and hit OK.\nThis will cause the \
workflow to spawn processes based on the available CPUs. Once finished the automatically determined CTF parameters will be stored in the EMAN2 database. It is often worthwhile to follow \
the following pattern:
1 - autofit
2 - manually fine-tune parameters for a few sets at different defocuses
3 - generate a structure factor using these sets
4 - rerun autofit
5 - manually check the fitting results. Correct any major errors.
"""
	warning_string = "\n\n\nNOTE: There are no particles currently associated with the project. Please go to the \"Particles\" task and import/box particles first."

	def __init__(self):
		E2CTFWorkFlowTask.__init__(self)
		self.window_title = "e2ctf Auto Fitting"

	def get_params(self):
		params = []
		p,n= self.get_ctf_param_table()
#		if n == 0 and False:
#			params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2CTFAutoFitTask.documentation_string+E2CTFAutoFitTask.warning_string,choices=None))
#		else:
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2CTFAutoFitTask.documentation_string,choices=None))
		params.append(p)
		self.add_general_params(params)

		return params
	
	def add_general_params(self,params):
		project_db = db_open_dict("bdb:project")
		db = db_open_dict(self.form_db_name)
		papix = ParamDef(name="global.apix",vartype="float",desc_short="A/pix for project",desc_long="The physical distance represented by the pixel spacing",property=None,defaultunits=project_db.get("global.apix",dfl=1.1),choices=None)
		pvolt = ParamDef(name="global.microscope_voltage",vartype="float",desc_short="Microscope voltage",desc_long="The operating voltage of the microscope in kilo volts",property=None,defaultunits=project_db.get("global.microscope_voltage",dfl=300),choices=None)
		pcs = ParamDef(name="global.microscope_cs",vartype="float",desc_short="Microscope Cs",desc_long="Microscope spherical aberration constant",property=None,defaultunits=project_db.get("global.microscope_cs",dfl=2.0),choices=None)
		pac = ParamDef(name="ac",vartype="float",desc_short="Amplitude contrast",desc_long="The amplitude contrast constant. It is recommended that this value is identical in all of your images.",property=None,defaultunits=db.get("ac",dfl=10),choices=None)
		pos = ParamDef(name="oversamp",vartype="int",desc_short="Oversampling",desc_long="If greater than 1, oversampling by this amount will be used when fitting CTF. It is not automatically applied when phase flipping particles since it would make the process irreversible.",property=None,defaultunits=db.get("oversamp",dfl=1),choices=None)
		pncp = ParamDef(name="global.num_cpus",vartype="int",desc_short="Number of CPUs",desc_long="Number of CPUS available for the project to use",property=None,defaultunits=project_db.get("global.num_cpus",dfl=num_cpus()),choices=None)
		pahp = ParamDef(name="autohp",vartype="boolean",desc_short="Auto high pass",desc_long="Automatic high pass filter of the SNR only to remove initial sharp peak, phase-flipped data is not directly affected (default false)",property=None,defaultunits=db.get("autohp",dfl=False),choices=None)
		pns = ParamDef(name="nosmooth",vartype="boolean",desc_short="No smoothing",desc_long="Disable smoothing of the background (running-average of the log with adjustment at the zeroes of the CTF)",property=None,defaultunits=db.get("nosmooth",dfl=False),choices=None)
		#db_close_dict(self.form_db_name)
		
		params.append([papix,pvolt])
		params.append([pcs,pac])
		params.append([pos,pncp])
		params.append([pahp,pns])

		#db_close_dict("bdb:project")
	
	def get_default_ctf_options(self,params):
		'''
		These are the options required to run pspec_and_ctf_fit in e2ctf.py
		'''
		
		error_message = []
		
		if not params.has_key("filenames"): return None # this is fine
		
		if  params.has_key("filenames") and len(params["filenames"]) == 0:
			self.run_select_files_msg()
			return None
		
		filenames = params["filenames"]
		boxsize = None
		db_file_names = []
		for i,name in enumerate(filenames):
			db_name= name
			db_file_names.append(db_name)
			if not file_exists(db_name):
				error_message.append("error, %s does not exist." %name)
				continue
			
			if boxsize == None:
				a = EMData()
				a.read_image(db_name,0,True)
				hdr = a.get_attr_dict()
				boxsize = hdr["nx"] # no consideration is given for non square images
				#db_close_dict(db_name)
			else:
				a = EMData()
				a.read_image(db_name,0,True)
				hdr = a.get_attr_dict()
				#db_close_dict(db_name)
				if boxsize != hdr["nx"]: # no consideration is given for non square images
					error_message.append("Can't run e2ctf on images with different box sizes.")
		
		if boxsize == None or boxsize < 2:
			error_message.append("error, boxsize is less than 2.")
		
		
		if boxsize == None:
			error_message.append("Something went wrong in the form")
			self.show_error_message(error_message)
			return None
		
		options = EmptyObject()
#		options.bgmask = boxsize/2.6
		options.filenames = db_file_names
		error_message.extend(self.append_general_options(options,params))
		
		if len(error_message) > 0:
			self.show_error_message(error_message)
			return None
	
		return options
	
	def append_general_options(self,options,params):
		'''
		This is done in more than one place hence this function
		'''
		
		options.nosmooth = params["nosmooth"]
		options.nonorm = False
		options.autohp = params["autohp"]
		options.invert = False
		options.oversamp = params["oversamp"]
		options.ac = params["ac"]
		options.apix = params["global.apix"]
		options.cs = params["global.microscope_cs"]
		options.voltage = params["global.microscope_voltage"]
		
		
		error_message = []
		if options.oversamp < 1:
			error_message.append("You must specify a value for oversamp that is atleast 1.")
		if options.apix <= 0:
			error_message.append("You must specify a positive non zero value for the apix.")
		if options.voltage <= 0:
			error_message.append("You must specify a positive non zero value for the voltage.")
		
		return error_message
		
	def on_form_ok(self,params):
		
		options = self.get_default_ctf_options(params)
		if options != None:
			self.write_db_entries(params)
			
			string_args = ["oversamp","ac","apix","cs","voltage"]
			bool_args = ["nosmooth","nonorm","autohp","invert"]
			additional_args = ["--autofit","--curdefocushint"]
			temp_file_name = "e2ctf_autofit_stdout.txt"
			self.spawn_task("e2ctf.py",options,string_args,bool_args,additional_args,temp_file_name)
			
			self.form.close()
			self.emit(QtCore.SIGNAL("task_idle"))
			
		else:
			return
		
	def on_form_close(self):
		# this is to avoid a task_idle signal, which would be incorrect if e2boxer is running
		self.emit(QtCore.SIGNAL("task_idle"))
			
class E2CTFAutoFitTaskGeneral(E2CTFAutoFitTask):
	'''
	This one has a generic url browser to get the input file names as opposed to the more rigid project based one
	'''
	documentation_string = "Use this tool to use e2ctf to generate ctf parameters for the particles located anywhere on disk. Use the browse buttons to locate the files on disk, enter the fitting parameters such as microscope voltage and angstrom per pixel, and hit OK. \nThis will cause the workflow to spawn processes based on the available CPUs. Once finished the automatically determined CTF parameters will be stored in the EMAN2 database."
	
	def __init__(self):
		E2CTFAutoFitTask.__init__(self)
		self.window_title = "e2ctf Auto Fitting"
		self.file_check = True
	def get_params(self):
		params = []		
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2CTFAutoFitTask.documentation_string,choices=None))
		
		params.append(ParamDef(name="filenames",vartype="url",desc_short="File Names",desc_long="The names of the particle files you want to generate automated ctf parameters for",property=None,defaultunits=[],choices=[]))
		
		self.add_general_params(params)

		return params
	
	def get_default_ctf_options(self,params):
		'''
		These are the options required to run pspec_and_ctf_fit in e2ctf.py
		'''
		
		if not params.has_key("filenames"): return None # this is fine
		
		if  params.has_key("filenames") and len(params["filenames"]) == 0:
			return None # this is fine
		
		
		filenames = params["filenames"]
		fine,message = check_files_are_em_images(filenames)
		
		if not fine:
			print message
			return None
		
		boxsize = None
		db_file_names = []
		for i,name in enumerate(filenames):
			a = EMData()
			a.read_image(name,0,True)
			if boxsize == None:
				boxsize = a.get_attr("nx") # no consideration is given for non square images
			elif boxsize != a.get_attr("nx"): # no consideration is given for non square images
					print "error, can't run e2ctf on images with different box sizes." # Specifically, I can not deduce the bgmask option for the group"
					return None
		
		if boxsize == None or boxsize < 2:
			print "error, boxsize is less than 2"
			return None
		
		options = EmptyObject()
#		options.bgmask = boxsize/2.6
		options.filenames = filenames
		self.append_general_options(options,params)
		return options

class E2CTFOutputTask(E2CTFWorkFlowTask):	
	"""Select the particle data for which you wish to generate phase flipped and/or Wiener filtered output and hit OK.\nThis will cause the workflow to spawn processes based on the available CPUs that write the output into a predefined location in the EMAN2 database.\nNote that the Wiener filtered output images are also phase flipped."""
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
		
		if  not params.has_key("filenames") or (params.has_key("filenames") and len(params["filenames"]) == 0):
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
			self.emit(QtCore.SIGNAL("task_idle"))
		else:
			return
		
class E2CTFSFOutputTask(E2CTFWorkFlowTask):	
	"""This task is for generating a structure factor file using e2ctf.py. Select the files you wish to use as the basis of the the structure factor approximation and hit OK. This will spawn e2ctf as a process, the output of which will be a file called 'strucfac.txt' """
	def __init__(self):
		E2CTFWorkFlowTask.__init__(self)
		self.window_title = "Generate e2ctf Structure Factor"

	def get_params(self):
		params = []		

		p,n = self.get_full_ctf_table()
		db = db_open_dict(self.form_db_name)
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=self.__doc__,choices=None))
		params.append(p)
		pos = ParamDef(name="oversamp",vartype="int",desc_short="Oversampling",desc_long="If greater than 1, oversampling by this amount will be used when fitting CTF. It is not automatically applied when phase flipping particles since it would make the process irreversible.",property=None,defaultunits=db.get("oversamp",dfl=1),choices=None)
		params.append(pos)
		
		return params
	
	def on_form_ok(self,params):

				
		if  not params.has_key("filenames") or  len(params["filenames"]) == 0:
			error("Please select files to process")
			return

		options = EmptyObject()
		options.filenames = params["filenames"]
		options.oversamp = params["oversamp"]
		
		string_args = []
		bool_args = []
		additional_args = ["--computesf"]
		temp_file_name = "e2ctf_sfoutput_stdout.txt"
		self.spawn_single_task("e2ctf.py",options,string_args,bool_args,additional_args,temp_file_name)
		
		self.form.close()
		self.emit(QtCore.SIGNAL("task_idle"))
		self.write_db_entries(params)


class E2CTFOutputTaskGeneral(E2CTFOutputTask):
	''' Use this form for generating CTF-related output. 
	'''
	warning_string = "\n\n\nNOTE: There are no CTF parameters currently stored for any images in the local database. You can change this by running automated fitting with e2ctf."
	
	def __init__(self):
		E2CTFOutputTask.__init__(self)
		self.window_title = "Generate e2ctf Output"
		self.names = None
		
	def set_names(self,names):
		self.names = names

	def get_custom_table(self,names):
		
		from emform import EM2DStackTable,EMFileTable,float_lt,int_lt
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
		
		if  params.has_key("filenames") and len(params["filenames"]) == 0:
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
		for k,v in params.items():
			self.write_db_entry(k,v)

		options = self.get_ctf_options(params)
		if options != None and len(options.filenames) > 0 and (options.wiener or options.phaseflip):
			
			string_args = []
			bool_args = ["wiener","phaseflip"]
			additional_args = []
			temp_file_name = "e2ctf_output_stdout.txt"
			self.spawn_single_task("e2ctf.py",options,string_args,bool_args,additional_args,temp_file_name)
			

			self.form.close()
			self.emit(QtCore.SIGNAL("task_idle"))
		else:
			self.form.close()
			self.emit(QtCore.SIGNAL("task_idle"))
			
	def on_form_cancel(self):
		
		self.form.close()
		self.form = None
		self.emit(QtCore.SIGNAL("task_idle"))
	
class E2CTFGuiTask(E2CTFWorkFlowTask):	
	documentation_string = "Autofitting tends to either work very well or get the defocus completely wrong. It is wise to at least quickly go through the data and insure that \
defocus values are reasonable. If not, roughly adjust the defocus and press the refit button. If you manually vary parameters, press save for each set, or your changes will \
be lost. B-factors are not as important as in EMAN1, and use the X-ray convention (4x the EMAN1 values). Try to get them in a reasonable range, at least. This is particularly \
important when manually fitting before determining a structure factor."
	warning_string = "\n\n\nNOTE: There are no particles associated with the project and/or there are no previously generated CTF parameters for these particles. To establish project particles go to the \"Particles\" task. To generate CTF parameters go to the \"Automated fitting - e2ctf\" task" 
	def __init__(self):
		E2CTFWorkFlowTask.__init__(self)
		self.window_title = "Launch e2ctf Interface"
		self.gui = None # will eventually be a e2ctf gui

	def get_params(self):

#		ptcl_names = self.get_particle_db_names(strip_ptcls=False) # particles in the project directory
#		if ptcl_names != None and len(ptcl_names) != 0: 
		p,n = self.get_ctf_param_table()
#		else:
#			n = 0
		params = []		
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2CTFGuiTask.documentation_string,choices=None))
		params.append(p)
		return params
	
	def get_default_ctf_options(self,params):
		'''
		These are the options required to run pspec_and_ctf_fit in e2ctf.py
		'''
		
		if  params.has_key("filenames") and len(params["filenames"]) == 0:
			self.run_select_files_msg()
			return None

		options = EmptyObject()
		filenames = params["filenames"]
#
		db_file_names = []
		for i,name in enumerate(filenames):
			db_name=name
			db_file_names.append(db_name)
			if not file_exists(db_name):
				print "No project particles entry exists for",name,"aborting."
				return None
		options.filenames = db_file_names
#		
		return options

	
	def on_form_ok(self,params):
		for k,v in params.items():
			self.write_db_entry(k,v)

		options = self.get_default_ctf_options(params)
		if options != None and len(options.filenames) > 0:
			
			img_sets = get_gui_arg_img_sets(options.filenames)
		
			init_sfcurve("auto")
			self.gui=GUIctf(get_application(),img_sets)
			self.emit(QtCore.SIGNAL("gui_running"), "CTF", self.gui) # so the desktop can prepare some space!
			self.form.close()
			QtCore.QObject.connect(self.gui,QtCore.SIGNAL("module_closed"), self.on_ctf_closed)
			self.gui.show_guis()
		else:
			return
	
	def on_ctf_closed(self):
		if self.gui != None:
			self.gui = None
			self.emit(QtCore.SIGNAL("gui_exit")) #
		
	def on_form_close(self):
		# this is to avoid a task_idle signal, which would be incorrect if e2boxer is running
		if self.gui == None:
			self.emit(QtCore.SIGNAL("task_idle"))
		else: pass
		
class E2CTFGuiTaskGeneral(E2CTFGuiTask):
	''' This one uses the names in the e2ctf.parms to generate it's table of options, not the particles in the particles directory
	'''
	warning_string = "\n\n\nNOTE: There are there are no previously generated CTF parameters. Please run automated fitting using e2ctf first." 
	
	documentation_string = "Write me"
	def __init__(self):
		E2CTFGuiTask.__init__(self)

	def get_params(self):
		params = []		
		
		names = self.get_names_with_ctf_params()
		n = [l[0] for l in names]
		p,num = self.get_ctf_param_table()
		
		if num == 0:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2CTFGuiTask.documentation_string+E2CTFGuiTaskGeneral.warning_string,choices=None))
		else:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2CTFGuiTaskGeneral.documentation_string,choices=None))
			params.append(p)
		return params
	
	def get_default_ctf_options(self,params):
		'''
		These are the options required to run pspec_and_ctf_fit in e2ctf.py
		'''
	
		if  params.has_key("filenames") and len(params["filenames"]) == 0:
			self.run_select_files_msg()
			return None # this is fine
		
		
		options = EmptyObject()

		selected_filenames = params["filenames"]
		
		filenames = []
		names = self.get_names_with_ctf_params()
		for name in names:
			if name[0] in selected_filenames:
				filenames.append(name[1])
		options.filenames = filenames
#		
		return options


class EMPartSetOptions:
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

		for name,d in main_dict.items():
			for filt,ptcl_name in d.items():
				if self.bdb_only != False and ( len(ptcl_name) > 3 and ptcl_name[:4] != "bdb:"): 
					continue
				name_map[ptcl_name] = name
				if filter_opts.has_key(filt):
					if self.image_count: filter_opts[filt] += EMUtil.get_image_count(ptcl_name)
					else: filter_opts[filt] += 1
					stacks_map[filt].append(ptcl_name)
				else:
					if self.image_count: filter_opts[filt] = EMUtil.get_image_count(ptcl_name)
					else: filter_opts[filt] = 1
					stacks_map[filt] = [ptcl_name]
		
		choices = []
		for filt,num in filter_opts.items():
			name = filt+" ("+str(num)+")"
			choices.append( name )
			stacks_name_map[name] = filt
		
		return stacks_map, stacks_name_map,choices,name_map
	
class EMParticleOptions(EMPartSetOptions):
	''' 
	e2refine2d and e2refine, from with the worklfow setting, both want to know the same 
	thing, basically, what are the available particles (including filtered options)
	'''
	
	def __init__(self,bdb_only=False):
		EMPartSetOptions.__init__(self,spr_ptcls_dict,bdb_only)
	
 		 	
class E2ParticleExamineChooseDataTask(ParticleWorkFlowTask):
	"""Choose the particle data you wish to examine. This will pop a second form listing the particles stacks along with other relevant information"""
	documentation_string = "On the next screen you will have the opportunity to eliminate 'bad' particles. This page selects which type of particles you \
want to look at for doing this. Normally Wiener filtered particles are used for this purpose. The selection here has no impact on the output generated from \
this stage. It is used for display purposes only !" 
	def __init__(self):
		ParticleWorkFlowTask.__init__(self)
		self.window_title = "Choose Data For Particle Examination"
		self.preferred_size = (480,300)
		self.form_db_name ="bdb:emform.workflow_particle_examine"

	def get_params(self):
		ptcl_opts = EMParticleOptions(bdb_only=True)
		self.particles_map, self.particles_name_map,choices,self.name_map = ptcl_opts.get_particle_options()
			
		params = []		
		
		#if as_string:
		#params.append(ParamDef(name="particle_set_choice",vartype="string",desc_long="Choose the particle data set you wish to use to generate a starting data for e2refine2d",desc_short=title,property=None,defaultunits=db.get("particle_set_choice",dfl=""),choices=choices))
		db = db_open_dict(self.form_db_name)
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=str(self.__doc__),choices=None))
		
		if len(self.particles_map) > 0 and len(choices) > 0:
			params.append(ParamDef(name="particle_set_choice",vartype="choice",desc_long="Choose the particle data set you wish to use to generate a starting data for e2refine2d",desc_short="Choose data",property=None,defaultunits=db.get("particle_set_choice",dfl=""),choices=choices))
		return params

	def on_form_ok(self,params):
		if params.has_key("particle_set_choice"):

			choice = params["particle_set_choice"]
			
#			name_map = {}
#			particles = self.particles_map[self.particles_name_map[choice]]
#			for name in particles:
#				if self.name_map.has_key(name):
#					name_map[name] = base_name(self.name_map[name])
#				else:
#					name_map[name] = base_name(name)
					
			self.emit(QtCore.SIGNAL("replace_task"),E2ParticleExamineTask(particles,self.name_map),"Particle Set Examination")
		
		self.form.close()
		self.form = None
		
		self.write_db_entries(params)
		

		
class E2ParticleExamineTask(E2CTFWorkFlowTask):
	'''This task is for defining bad particles. Double click on the particle stack you want to examine and then (once the stack viewer has loaded) click on particles to define them as bad. You can come back and review your bad particle selections at any time. You can write stacks that exclude particles you've defined as bad using 'Make Particle Set' tool.'''
	def __init__(self,particle_stacks=[],name_map={}):
		E2CTFWorkFlowTask.__init__(self)
		self.window_title = "Particle Set Examination"
		self.particle_stacks = particle_stacks
		self.name_map = name_map

	def get_params(self):
		params = []
		from emform import EM2DStackExamineTable,EMFileTable,int_lt
		table = EM2DStackExamineTable(self.particle_stacks,desc_short="Particles",desc_long="",name_map=self.name_map)
		bpc = E2ParticleExamineTask.BadParticlesColumn(self.name_map)
		# think bpc will be devoured by the garbage collector? Think again! the table has a reference to one of bpc's functions, and this is enough to prevent death
		table.insert_column_data(1,EMFileTable.EMColumnData("Bad Particles",bpc.num_bad_particles,"The number of particles you have identified as bad",int_lt))
		table.register_animated_column("Bad Particles")
		
		p,n = self.get_particle_selection_table(self.particle_stacks,table) # I wish this came from a "has a" not an "is a" relationship. But too late to chaneg
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=self.__doc__,choices=None))
		params.append(p)

		return params
	
	class BadParticlesColumn:
		def __init__(self,name_map):
			self.name_map = name_map
	
		def num_bad_particles(self,filename):
			if not db_check_dict("bdb:select"): return "0"
			if not self.name_map.has_key(filename):
				EMErrorMessageDisplay.run(["The name map doesn't have the %s key" %filename])
				return "0"
		
			db_name = self.name_map[filename]
			db = db_open_dict("bdb:select",ro=True)
			if not db.has_key(db_name): return "0"
			set_list = db[db_name]
			if set_list == None: return "0"
			else: return str(len(set_list))
	
	
	def on_form_ok(self,params):
		self.form.close()
		self.form = None
		self.emit(QtCore.SIGNAL("task_idle"))
		
	def on_form_close(self):
		# this is to avoid a task_idle signal, which would be incorrect if e2boxer is running
		self.emit(QtCore.SIGNAL("task_idle"))
		
	def on_form_cancel(self):
		
		self.form.close()
		self.form = None
		self.emit(QtCore.SIGNAL("task_idle"))

 		 	
class E2MakeSetChooseDataTask(E2ParticleExamineChooseDataTask):
	"""Choose the particle data you wish to examine. This will pop a second form listing the particles stacks along with other relevant information"""
	documentation_string = "Choose the data for particle examination" 
	def __init__(self):
		E2ParticleExamineChooseDataTask.__init__(self)
		self.window_title = "Choose Data For Making Particle Sets"
		self.preferred_size = (480,300)
		self.form_db_name ="bdb:emform.make_ptcl_set"

	def on_form_ok(self,params):
		if params.has_key("particle_set_choice"):

			choice = params["particle_set_choice"]
			
			name_map = {}
			for filt_name, particles in self.particles_map.items():
				for name in particles:
					if self.name_map.has_key(name):
						name_map[name] = base_name(self.name_map[name])
					else:
						name_map[name] = base_name(name)
			particles = self.particles_map[self.particles_name_map[choice]]
			
			self.emit(QtCore.SIGNAL("replace_task"),E2MakeSetTask(particles,name_map,self.name_map),"Make Particle Sets")
		
		self.form.close()
		self.form = None
		
		self.write_db_entries(params)
		
class E2MakeSetTask(E2ParticleExamineTask):
	'''This task is for making agglomerated particle sets which function as input for the 2D and 3D refinement tasks. Sets for any linked particle sets (such as phase flipped or Wiener filtered data) are automatically generated. You can choose to exclude any bad particles that you have defined.'''
	def __init__(self,particle_stacks=[],name_map={},data_name_map={}):
		E2ParticleExamineTask.__init__(self,particle_stacks,data_name_map)
		self.window_title = "Make Particle Sets"
		self.form_db_name = "bdb:emform.make_ptcl_set"
		self.data_name_map = data_name_map
		
	def get_params(self):
		params = E2ParticleExamineTask.get_params(self)
		
		db = db_open_dict(self.form_db_name) # see eman wiki for a list of what args are kept in this db
	   	pstack_name = ParamDef(name="stack_name",vartype="string",desc_short="Set Name",desc_long="What you want to call this stack. Leave out file types (such as hdf,bdb etc). Names for any linked particle sets will be generated automatically.",property=None,defaultunits=db.get("stack_name",dfl="set"),choices=None)
		pexclude_bad = ParamDef(name="exclude_bad",vartype="boolean",desc_short="Exclude Bad Particles",desc_long="Exclude bad particles from generated stack(s)",property=None,defaultunits=db.get("exclude_bad",dfl=True),choices=None)
 	
 		params.append(pstack_name)
 		params.append(pexclude_bad)
		return params
	
	def check_params(self,params):
		error_message = []
		if len(params["filenames"]) == 0:
			error_message.append("Please select files for processing")
		
		if len(params["stack_name"]) == 0:
			error_message.append("Please supply a stack name")
		else:
			pass
		
		if len(error_message) == 0:
			return True
		else:
			EMErrorMessageDisplay.run(error_message)
			return False
		
	def execute(self,params):
		db = db_open_dict("bdb:project")
		if not db.has_key(spr_ptcls_dict):
			EMErrorMessageDisplay.run("Something is wrong, the global.spr_ptcls_dict is supposed to exist")
			return False
		
		base_stack_root = params["stack_name"]
		test_stack_name = "bdb:sets#"+base_stack_root
		if file_exists(test_stack_name):
			EMErrorMessageDisplay.run("The %s stack already exists. Remove it, or try a different name" %test_stack_name)
			return False

		data_dict = EMProjectDataDict(spr_ptcls_dict)
		project_data = data_dict.get_data_dict()
		
		output_stacks = {}
		output_stacks[""] = []
		stack_type_map = {}
		
		print params["filenames"]
		for name in params["filenames"]:
			root_name = self.data_name_map[name]
			dict = project_data[root_name]
			for filt,name in dict.items():
#				f = "_"+filt.split()[0].lower()			# This produces redundant names sometimes
				f ="_"+filt.lower().replace(" ","_")	# This doesn't

				if not stack_type_map.has_key(f):
					stack_type_map[f] = filt
				
				if output_stacks.has_key(f): output_stacks[f].append(name)
				else: output_stacks[f] = [name]
		
		for f in output_stacks.keys():
			test_stack_name = "bdb:sets#"+base_stack_root+f
			if file_exists(test_stack_name):
				EMErrorMessageDisplay.run("The %s stack already exists. Remove it, or try a different name" %test_stack_name)
				return False
			
		stacks_map = {}
	 	progress = QtGui.QProgressDialog("Making virtual stacks...", "Abort import", 0, 1000,None)
	 	progress.setWindowIcon(QtGui.QIcon(get_image_directory() + "/eman.png"))
		progress.show()
		n=-1
		print output_stacks
		for key,filenames in output_stacks.items():
			n+=1
			if len(filenames) == 0:
				print "Warning, there were no files in the list"
				continue
			if len(filenames[0]) > 3 and filenames[0][:4] == "bdb:":
				success,cmd = self.make_v_stack(filenames,base_stack_root+key,"sets",params["exclude_bad"],progress,n,len(output_stacks.items()))
				
			else:
				EMErrorMessageDisplay.run("The generation of stacks for flat (non database) files is currently disabled. A particle set (%s) is being ignored" %stack_type_map[key], "Warning")
				continue
				# This will work eventually - the only problem is that the exclusioons are not factored in, that's all
				#success,cmd = self.make_stack(filenames,base_stack_root+key,"stacks/")
				
			if not success:
				EMErrorMessageDisplay.run("The command %s failed. Please try again" %cmd, "Failure")
				return False
			else:
				if key != "": stacks_map[stack_type_map[key]] ="bdb:sets#"+ base_stack_root+key
			
		progress.close()
		stack_root = params["stack_name"]	
	
		stacks_data_dict = EMProjectDataDict(spr_sets_dict)
		update_map = {}
		update_map["bdb:sets#"+base_stack_root] = stacks_map
		stacks_data_dict.update(update_map)

		
		return True
		
	def make_v_stack(self,filenames,out_name,path,exclude_bad,progress=None,N=0,M=1):
	 	
		if db_check_dict("bdb:select"):
			select_db = db_open_dict("bdb:select")
		else:
			select_db = None
		cmd = ""
		
		get_application().setOverrideCursor(Qt.BusyCursor)
		for i,name in enumerate(filenames):
			cmd = "e2bdb.py"
			no_exc = True
			if exclude_bad and select_db != None:
				exc_db = self.name_map[name]
				if select_db.has_key(exc_db):
					cmd += " "+name+"?exclude."+exc_db
					no_exc=False
			if no_exc:
				cmd += " "+name
				
			cmd += " --appendvstack=bdb:"+path+"#"+out_name
			print cmd
#			print N,i,len(filenames),M
			success = (os.system(cmd) in (0,11,12))
			if not success:
				get_application().setOverrideCursor(Qt.ArrowCursor)
				return False,cmd
			else:
				get_application().setOverrideCursor(Qt.ArrowCursor)
				progress.setValue(1000*(N+(i+1.0)/len(filenames))/M)
				get_application().processEvents()
				
			try:
				if progress.wasCanceled():
					get_application().setOverrideCursor(Qt.ArrowCursor)
					return False,"did not finish. The operation was cancelled"
			except: 
				pass
		
		get_application().setOverrideCursor(Qt.ArrowCursor)
		return True,cmd

	def make_stack(self,filenames,out_name,path):
		'''
		This one's a bit more specialized to handle flat files and avoid massive copying
		'''

			
		progress = QtGui.QProgressDialog("Importing files into database...", "Abort import", 0, len(filenames),None)
		progress.setWindowIcon(QtGui.QIcon(get_image_directory() + "/eman.png"))
		progress.show()

		for i,name in enumerate(filenames):
			cmd = "e2proc2d.py"
 			cmd += " "+name
 			cmd += " "+path+out_name
 			success = (os.system(cmd) in (0,12))
 			if not success:
 				progress.close()
 				return False,cmd
 			else:
 				progress.setValue(i+1)
 				get_application().processEvents()
		
		progress.close()
		
		return True,cmd

			
	def on_form_ok(self,params):
		proceed = self.check_params(params)
		if not proceed: return
		
		
		proceed = self.execute(params)
		if not proceed: return
		
		self.form.close()
		self.form = None
		self.emit(QtCore.SIGNAL("task_idle"))
	
		self.write_db_entries(params)

class EMSetReportTask(ParticleWorkFlowTask):
	'''This form displays a list of the stacks you currently have associated with the project
	'''
	def __init__(self):
		ParticleWorkFlowTask.__init__(self)
		self.window_title = "Project Sets"
		self.stack_data_at_init = []

	def get_project_particle_table(self):
		data_dict = EMProjectDataDict(spr_sets_dict)
		stack_data = data_dict.get_data_dict()
		stack_names = stack_data.keys()
		
		self.stack_data_at_init = stack_data # so if the user hits cancel this can be reset
		
		init_stacks_map = {}
		
		from emform import EM2DStackTable,EMFileTable,int_lt
		table = EM2DStackTable(stack_names,desc_short="Project Particle Sets",desc_long="",enable_save=False)
#		context_menu_data = EMRawDataReportTask.ProjectListContextMenu(spr_sets_dict,remove_only=True)
#		table.add_context_menu_data(context_menu_data)
		#table.add_button_data(EMRawDataReportTask.ProjectAddRawDataButton(table,context_menu_data))
#		table.insert_column_data(1,EMFileTable.EMColumnData("Particles On Disk",EMParticleReportTask.get_num_ptcls,"Particles in this image"))
#		table.insert_column_data(2,EMFileTable.EMColumnData("Particle Dims",EMParticleReportTask.get_particle_dims,"The dimensions of the particles that are stored on disk"))
#		
#		project_db = db_open_dict("bdb:project",ro=True)
#		if project_db.has_key(spr_sets_dict):
			
		if len(stack_names) != 0:
			
			filt_options = []
			for key,value in stack_data.items():
				for filt,name in value.items():
					if filt not in filt_options:
						filt_options.append(filt)
			
			for filt in filt_options:
				fsi = EMSetReportTask.FilteredSetInfo(filt)
				table.add_column_data(EMFileTable.EMColumnData(filt,fsi.num_particles_in_filt_stack,"The number of particles in the associated stack (if it exists)",int_lt))

		return table
	
	class FilteredSetInfo:
		def __init__(self,filt):
			self.filt = filt
			data_dict = EMProjectDataDict(spr_sets_dict)
			self.stacks_map = data_dict.get_data_dict()
			
		def num_particles_in_filt_stack(self,filename):
			d = self.stacks_map[filename]
			if d.has_key(self.filt):
				return str(EMUtil.get_image_count(d[self.filt]))
			else: return "-"
			
	
	def on_form_cancel(self):
		self.recover_original_raw_data_list()
		
		self.form.close()
		self.form = None
		self.emit(QtCore.SIGNAL("task_idle"))
	
	def recover_original_raw_data_list(self):
		'''
		Called if the user hits cancel - if they removed some files or added files the changes
		are not saved unless the user hits ok
		'''
		project_db = db_open_dict("bdb:project")
		project_db[spr_sets_dict] = self.stack_data_at_init

	def get_params(self):
		params = []
		
	
		table = self.get_project_particle_table()
		
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=self.__doc__,choices=None))
		params.append(table)  
		
		return params

class E2Refine2DReportTask(ParticleWorkFlowTask):
	"""This form displays reference free class averages that are associated with the project."""
	
	def __init__(self):
		ParticleWorkFlowTask.__init__(self)
		self.window_title = "2D Class Averages "
		
		
	def get_params(self):
		params = []		
		
		self.rfcat = E2RefFreeClassAveTool()
		p,n = self.rfcat.get_reference_free_class_averages_table()

		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=self.__doc__,choices=None))
		params.append(p)
		return params
	
	def on_form_cancel(self):
		self.rfcat.recover_original_raw_data_list()
		
		self.form.close()
		self.form = None
		self.emit(QtCore.SIGNAL("task_idle"))


class EMClassificationTools(ParticleWorkFlowTask):
	'''
	Encapsulation of common functionality.
	Specifically - e2classaverage and e2simmx pages - both of which are used by e2refine and e2refine2d
	'''
	def __init__(self):
		ParticleWorkFlowTask.__init__(self)
		
	def get_simmx_page(self):
		
		params = []
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2RefineParticlesTask.simmx_documentation,choices=None))

		db = db_open_dict(self.form_db_name)
		pshrink = ParamDef(name="shrink",vartype="int",desc_short="Shrink",desc_long="Shrink the data at various stages in refinement, for speed purposes",property=None,defaultunits=db.get("shrink",dfl=2),choices=[])
		ptwostage = ParamDef(name="twostage",vartype="int",desc_short="2 Stage Simmx",desc_long="This will determine the orientation in 2 stages, usually 5-10x faster, number specifies shrink factor for first stage, 0 disables",property=None,defaultunits=db.get("twostage",dfl=0),choices=[])
		pprefilt = ParamDef(name="prefilt",vartype="boolean",desc_short="PS Match Ref",desc_long="Filter references to match particles before alignment. Works best with usefilt -> Wiener filtered particles",property=None,defaultunits=db.get("prefilt",dfl=0),choices=[])

		
		params.append([pshrink,ptwostage,pprefilt])
		params.extend(self.get_cls_simmx_params(parameter_prefix="sim"))
		
		#db_close_dict(self.form_db_name)
		return ["Simmx",params]
	

	def get_classaverage_page(self,include_sep=True):
		params = []
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2RefineParticlesTask.class_documentation,choices=None))

		db = db_open_dict(self.form_db_name)
		
		psep = ParamDef(name="sep",vartype="int",desc_short="Class separation",desc_long="The number of classes a particle can contribute towards",property=None,defaultunits=db.get("sep",dfl=1),choices=[])
		piter = ParamDef(name="classiter",vartype="int",desc_short="Averaging iterations",desc_long="The number of class averaging iterations",property=None,defaultunits=db.get("classiter",dfl=2),choices=[])
		
		averagers = self.get_averagers_list()
		averagers.sort()
		paverager =  ParamDef("classaverager",vartype="string",desc_short="Averager",desc_long="The method used for generating class averages",property=None,defaultunits=db.get("classaverager",dfl="mean"),choices=averagers)
		
		pkeep = ParamDef(name="classkeep",vartype="float",desc_short="keep",desc_long="The fraction of particles to keep in each class average. If sigma based is checked this value is interpreted in standard deviations from the mean instead",property=None,defaultunits=db.get("classkeep",dfl=0.8),choices=[])
		pkeepsig = ParamDef(name="classkeepsig",vartype="boolean",desc_short="Sigma based",desc_long="If checked the keep value is interpreted in standard deviations from the mean instead of basic ratio",property=None,defaultunits=db.get("classkeepsig",dfl=True),choices=[])
		prefsf = ParamDef(name="classrefsf",vartype="boolean",desc_short="Set proj SF on avg",desc_long="If checked each class average will be filtered to match the 1-D structure factor of the reference projection",property=None,defaultunits=db.get("classrefsf",dfl=False),choices=[])
		
		pnormproc =  ParamDef("classnormproc",vartype="string",desc_short="Normalization processor",desc_long="The normalization method applied to the particles before averaging",property=None,defaultunits=db.get("classnormproc",dfl="normalize.edgemean"),choices=["normalize","normalize.edgemean","None"])
		
		#db_close_dict(self.form_db_name)
		
		if include_sep: 
			params.append([piter,psep])
			params.append([pkeep,pkeepsig])
		else: # this happens in the eotest
			# I made the next row longer because it seemed like there was room
			params.append([piter,pkeep,pkeepsig])
		params.append([pnormproc])
		params.append([paverager,prefsf])
		params.extend(self.get_cls_simmx_params(parameter_prefix="class"))

		return ["Class averaging",params]
	
	def get_cls_simmx_params(self,parameter_prefix=""):
		params = []
		aligners = self.get_aligners_list()
		cmps = self.get_cmps_list()
		
		aligners.sort()
		cmps.sort()
		
		db = db_open_dict(self.form_db_name)
		
		palign =  ParamDef(name=parameter_prefix+"align",vartype="string",desc_short="Aligner",desc_long="The aligner being used",property=None,defaultunits=db.get(parameter_prefix+"align",dfl="rotate_translate_flip"),choices=aligners)
		palignargs =  ParamDef(name=parameter_prefix+"alignargs",vartype="string",desc_short="params",desc_long="Parameters for the aligner, see \"e2help.py aligners\"",property=None,defaultunits=db.get(parameter_prefix+"alignargs",dfl=""),choices=[])
		
		paligncmp =  ParamDef(name=parameter_prefix+"aligncmp",vartype="string",desc_short="Align comparator",desc_long="The comparator being used",property=None,defaultunits=db.get(parameter_prefix+"aligncmp",dfl="ccc"),choices=cmps)
		paligncmpargs =  ParamDef(name=parameter_prefix+"aligncmpargs",vartype="string",desc_short="params",desc_long="Parameters for this comparator, see \"e2help.py cmps\"",property=None,defaultunits=db.get(parameter_prefix+"aligncmpargs",dfl=""),choices=[])	
		
		
		pralign =  ParamDef(name=parameter_prefix+"ralign",vartype="string",desc_short="Refine aligner",desc_long="The refine aligner being used",property=None,defaultunits=db.get(parameter_prefix+"ralign",dfl="None"),choices=["None","refine"])
		pralignargs =  ParamDef(name=parameter_prefix+"ralignargs",vartype="string",desc_short="params",desc_long="Parameters for this aligner, see \"e2help.py aligners\"",property=None,defaultunits=db.get(parameter_prefix+"ralignargs",dfl=""),choices=[])
		
		praligncmp =  ParamDef(name=parameter_prefix+"raligncmp",vartype="string",desc_short="Refine align comparator",desc_long="The comparator being used for refine alignment",property=None,defaultunits=db.get(parameter_prefix+"raligncmp",dfl="ccc"),choices=cmps)
		praligncmpargs =  ParamDef(name=parameter_prefix+"raligncmpargs",vartype="string",desc_short="params",desc_long="Parameters for thos comparator, see \"e2help.py cmps\"",property=None,defaultunits=db.get(parameter_prefix+"raligncmpargs",dfl=""),choices=[])	
		
		pcmp  =  ParamDef(name=parameter_prefix+"cmp",vartype="string",desc_short="Main comparator",desc_long="The comparator to determine the final quality metric",defaultunits=db.get(parameter_prefix+"cmp",dfl="ccc"),choices=cmps)
		pcmpargs =  ParamDef(name=parameter_prefix+"cmpargs",vartype="string",desc_short="params",desc_long="Parameters for the this comparator, see \"e2help.py cmps\"",property=None,defaultunits=db.get(parameter_prefix+"cmpargs",dfl=""),choices=[])	
	
		#db_close_dict(self.form_db_name)

		params.append([pcmp,pcmpargs])
		params.append([palign,palignargs])
		params.append([paligncmp,paligncmpargs])
		params.append([pralign,pralignargs])
		params.append([praligncmp,praligncmpargs])
		
		return params
	
	def check_aligners_and_cmps(self,params,options,parameter_prefix="class",page="Class averaging"):
		
		error_message = []
		vals = []
		
		vals.append(["align","alignargs"])
		vals.append(["ralign","ralignargs"])
		vals.append(["aligncmp","aligncmpargs"])
		vals.append(["raligncmp","raligncmpargs"])
		vals.append(["cmp","cmpargs"])
		
		
		for v in vals:
			v[0] = parameter_prefix + v[0]
			v[1] = parameter_prefix + v[1]
		
		for v in vals:
			setattr(options,v[0],params[v[0]])
			setattr(options,v[1],params[v[1]])
		
		for v in vals:
			if getattr(options,v[0]) == "None": setattr(options,v[0],None)
			elif len(getattr(options,v[1])) != 0: setattr(options,v[0],getattr(options,v[0])+":"+getattr(options,v[1]))
		
		
		for i,v in enumerate(vals):
			arg = getattr(options,v[0])
			if arg != None:
				if i > 1: # its a cmp, yes a hack but I have no time
					if not check_eman2_type(arg,Cmps,"Cmp",False): error_message.append("There is problem with the " +v[0]+ " comparator argument in the "+page+" page.")
				else:
					if not check_eman2_type(arg,Aligners,"Aligner",False): error_message.append("There is problem with the " +v[0]+ " aligner argument in the "+page+" page.")
  	
  		return error_message
  	
  	def add_classaverage_args(self,options,string_args,bool_args,additional_args,include_sep=True):
		
		optionals = ["classcmp","classalign","classaligncmp","classralign","classraligncmp"]
		for opt in optionals:
			if getattr(options,opt) != None: string_args.append(opt)
			
		string_args.extend(["classiter","classkeep","classnormproc","classaverager"])
		if include_sep: string_args.append("sep")
		bool_args.append("classkeepsig")
		bool_args.append("classrefsf")
	
	def check_classaverage_page(self,params,options):
		error_message = []
		
		if params.has_key("sep") and params["sep"] <= 0: # sometimes this key is absent (from the e2eotest form)
			error_message.append("The separation argument in the Class average page must be atleast 1")
		
		if params["classiter"] < 0:
			error_message.append("The number of class averaging iterations must be atleast 0")
		
		if params["classkeepsig"] == False:
			if params["classkeep"] < 0 or params["classkeep"] > 1:
				error_message.append("The keep parameter in the Class average page must be between 0 and 1. This does not hold if the \'Sigma based\' option is selected.")
				
		error_message.extend(self.check_aligners_and_cmps(params,options,"class", "Class average"))
		
		if len(error_message) > 0: return error_message # calling program should act and discontinue
		
		if params.has_key("sep"): options.sep = params["sep"] # sometimes this key is absent (from the e2eotest form)
		options.classkeep = params["classkeep"]
		options.classkeepsig = params["classkeepsig"]
		options.classnormproc = params["classnormproc"]
		options.classiter = params["classiter"]
		options.classrefsf=params["classrefsf"]
		
		options.classaverager = params["classaverager"] # at the moment there are no extra averager parameter, but if that changes then the parameteres would have to be checked
		
		return error_message

	
	def add_simmx_args(self,options,string_args,bool_args,additional_args,include_shrink=True):
		
		optionals = ["simcmp","simalign","simaligncmp","simralign","simraligncmp"]
		for opt in optionals:
			if getattr(options,opt) != None: string_args.append(opt)
		
		if include_shrink and options.shrink > 1: string_args.append("shrink") # e2simmx doesn't like it if shrink is 1
		if options.twostage>0 : string_args.append("twostage")
		try: 
			if options.prefilt : string_args.append("prefilt")
		except :
			pass
#			traceback.print_exc()
#			print "Error setting prefilt option"
	
	def check_simmx_page(self,params,options):
		error_message = []
		
		if not params.has_key("shrink"): params["shrink"] = 1
		if not params.has_key("twostage"): params["twostage"] = 0
		
		if params["shrink"] <= 0:
			error_message.append("The shrink argument in the simmx page must be atleast 1")

		if params["twostage"] < 0:
			error_message.append("The shrink argument in the simmx page must be atleast 0")
			
		options.shrink=params["shrink"]
		options.twostage=params["twostage"]
		
		error_message.extend(self.check_aligners_and_cmps(params,options,"sim","Simmx"))
		
		return error_message

class E2Refine2DTask(EMClassificationTools):
	'''
	Common e2refine2D functionality
	'''
	documentation_string = "This form is a way for the user to supply arguments to and execute e2refine2d.py"
	def __init__(self):
		EMClassificationTools.__init__(self)
		self.window_title = "Run e2refine2d"
		self.form_db_name = "bdb:emform.e2refine2d"
		self.workflow_setting = True # you can turn this off if you don't want results stored workflow list of reference free class averages
		
	def get_general_params(self):
		db = db_open_dict(self.form_db_name)
		
		params = []		
		piter = ParamDef(name="iter",vartype="int",desc_short="Refinement iterations",desc_long="The number of times the e2refine2d svd-based class averaging procedure is iterated",property=None,defaultunits=db.get("iter",dfl=8),choices=[])
		pnaliref = ParamDef(name="naliref",vartype="int",desc_short="# alignment references",desc_long="The number of alignment references to use when determining particle orientations",property=None,defaultunits=db.get("naliref",dfl=6),choices=[])
		pnbasisfp = ParamDef(name="nbasisfp",vartype="int",desc_short="# basis fp",desc_long="The number of MSA basis vectors to use when classifying",property=None,defaultunits=db.get("nbasisfp",dfl=8),choices=[])
		pncls = ParamDef(name="ncls",vartype="int",desc_short="# classes",desc_long="The number of classes to produce",property=None,defaultunits=db.get("ncls",dfl=32),choices=[])

		pnp = ParamDef(name="normproj",vartype="boolean",desc_short="Normalize projection vectors",desc_long="Normalizes each projected vector into the MSA subspace",property=None,defaultunits=db.get("normproj",dfl=True),choices=None)
		
		project_db = db_open_dict("bdb:project")
		#pncp = ParamDef(name="global.num_cpus",vartype="int",desc_short="Number of CPUs",desc_long="Number of CPUS available for e2refine2d to use",property=None,defaultunits=project_db.get("global.num_cpus",dfl=num_cpus()),choices=None)
		pinitclasses =  ParamDef(name="initial",vartype="string",desc_short="Initial class averages",desc_long="A file (full path) containing starting class averages. If note specificed will generate starting class averages automatically.",property=None,defaultunits=db.get("initial",dfl=""),choices=[])	
		
		#db_close_dict(self.form_db_name)
		#db_close_dict("bdb:project")
		
		params.append([pncls,pnp])
		params.append(piter)
		params.append([pnaliref,pnbasisfp])
		params.append([pinitclasses])
		
		pparallel = ParamDef(name="parallel",vartype="string",desc_short="Parallel",desc_long="Parallel arguments (advanced). Leave blank if unsure",property=None,defaultunits=db.get("parallel",dfl=""),choices=None)
		params.append(pparallel)
		return params

	def check_main_page(self,params,options):
		error_message = []
		interested_in = ["initial","iter","naliref","nbasisfp","ncls","filenames","normproj"] #including filenames could be a hack
		for val in interested_in:
			setattr(options,val,params[val])
		
		#options.parallel = params["global.num_cpus"]
		if options.initial == "" or str(options.initial) == "None" or str(options.initial) == "none":
			options.initial = None
		
		if params.has_key("parallel") and len(params["parallel"]) > 0:
			options.parallel = params["parallel"]
		
		if options.initial != None and len(options.initial) > 0:
	 		if not file_exists(options.initial):
	 			error_message.append("The initial class averages file you specified (%s) does not exist." %(options.initial))
	 		
	 	if options.iter < 0: error_message.append("The number of e2refine2d iterations must be atleast 0.")
		if options.naliref < 1:	error_message.append("The number alignment references must be atleast 1.")
	  	if options.nbasisfp < 1: error_message.append("The number of MSA basis vectors must be atleast 1.")
	  	if options.ncls < 2: error_message.append("The number of classes must be atleast 2.")
	  	#if options.parallel < 1: error_message.append("The number CPUs availables must be atleast 1.")
  		
  		if len(error_message) != 0:
 			return error_message
	 	
	 	# if we make it here we are almost definitely good to go, the only thing that can fail is the e2bdb or e2proc2d commands
	 	options.path = numbered_path("r2d",True)
	 	if len(options.filenames) > 0 and options.filenames[0][:4] == "bdb:": 
	 		if len(options.filenames) > 1:
	 			bdb_success, bdb_cmd = self.make_v_stack(options)
	 			input = "bdb:"+options.path+"#all"
	 		else:
	 			input = options.filenames[0]
	 			bdb_success = True
	 			
		 	if bdb_success:
		 		if options.shrink > 1:
		 			cmd = "e2proc2d.py"
		 			cmd += " " + input
		 			options.input =  "bdb:"+options.path+"#all"+str(options.shrink)
		 			cmd += " "+options.input
		 			cmd += " --process=math.meanshrink:n="+str(options.shrink)
		 			
		 			get_application().setOverrideCursor(Qt.BusyCursor)
		 			success = (os.system(cmd) in (0,12))
		 			get_application().setOverrideCursor(Qt.ArrowCursor)
		 			
		 			if not success:
		 				return ["e2proc2d.py shrinking command failed. This command was\n" + cmd +"\nTry again please. If the failure occurs a second time please contact developers."]

		 		else: options.input = input
		 		
		 		options.filenames = [] # this is so spawn_task doesn't supply args to e2refine2d.py
		 		return [] # THIS IS IT, THE POINT OF SUCCESS - returning this means e2refine2d.py is good to go using the given parameters
		 	else:
		 		return ["e2bdb.py command failed. The command was\n" + bdb_cmd +"\nTry again please. If the failure occurs a second time please contact developers"]
		else:
			return self.process_specified_files(options,params)
  		
  		return error_message
  	
  	def add_main_args(self,options,string_args,bool_args,additional_args,):
  		string_args.extend( ["iter","naliref","nbasisfp","path","input","ncls"] ) # May 2009, took "parallel" out
		bool_args.append("normproj")
		optionals = ["initial","parallel"]
		for opt in optionals:
			if hasattr(options,opt):
				if getattr(options,opt) != None: string_args.append(opt)
	 		  
	def make_v_stack(self,options):
	 	
	 	cmd = "e2bdb.py"

	 	for name in options.filenames:
	 		cmd += " "+name
	 	
	 	cmd += " --makevstack=bdb:"+options.path+"#all"
	 	
	 	get_application().setOverrideCursor(Qt.BusyCursor)
	 	success = os.system(cmd)
		print "mvs ",success
		success = (success in (0,12))
	 	get_application().setOverrideCursor(Qt.ArrowCursor)
	 	return success,cmd

	def process_specified_files(self,options,params):
		error_message = []
		if not params.has_key("filenames") or len(params["filenames"]) == 0:
			return ["Please specify files to process"]
#			error_message.append("Please specify files to process")
#			self.show_error_message(error_message)
#			return None
		
		rx = -1
		a = EMData()
		for file in params["filenames"]:
			if not file_exists(file): error_message.append("%s does not exist" %s)
			else:
				try:
					nx,ny,nz = gimme_image_dimensions3D(file)
				except:
					error_message.append("%s is not a valid EM image" %file)
					continue
				
				
				if nz > 1:
					error_message.append("%s is 3D, refine2d works on 2D images only" &file)
					break
				else:
					if nx != ny:
						error_message.append("%s contains images that are not square. Refine2d requires square images" %file)
						break
					if rx == -1:
						rx = nx
					elif nx != rx:
						error_message.append("The dimensions of the image files you specified are not uniform")
						break
						
		if rx == -1:
			error_message.append("Couldn't determine image dimensions")
		
		if len(error_message) != 0:
			return error_message
			#self.show_error_message(error_message)
			#return None
		
		# if we make it here we are good to go
		options.input = None
		options.filenames = [] # important - makes the spawn_process task work
		if len(params["filenames"]) > 1 or options.shrink > 1:
			progress = EMProgressDialog("Processing files...", "Abort import", 0, len(params["filenames"]),None)
			progress.show()
		
			for i,file in enumerate(params["filenames"]):
				cmd = "e2proc2d.py"
				cmd += " %s" %file # the input name
				
				if options.shrink > 1: out_name = "bdb:"+options.path+"#all"+str(options.shrink)
				else: out_name = "bdb:"+options.path+"#all"
				
				if options.input == None: options.input = out_name
	 			
	 			cmd += " %s" %out_name # the output name
	 			
	 			if options.shrink > 1: 
	 				cmd += " --process=math.meanshrink:n="+str(options.shrink)
	 			
	 			get_application().setOverrideCursor(Qt.BusyCursor)
	 			success = (os.system(cmd) in (0,12))
	 			get_application().setOverrideCursor(Qt.ArrowCursor)
	 			if not success:
	 				return ["Command %s failed" %cmd]
	 			
	 			
	 			progress.setValue(i+1)
	 			get_application().processEvents()
	 			if progress.wasCanceled():
	 				db_remove_dict(options.input)
	 				progress.close()
	 				return ["Processing was cancelled"]
	 			
	 		try: db_close_dict(out_name)
	 		except:
	 			print "db close dict failed",out_name
	 			
	 		progress.close()
	 		
	 	else:
	 		options.input = params["filenames"][0]
	 	
	 	return []

class EMSetsOptions(EMPartSetOptions):
	''' 
	e2refine2d and e2refine, from with the worklfow setting, both want to know the same 
	thing, basically, what are the available stacks (including filtered options)
	'''
	
	def __init__(self,bdb_only=False):
		EMPartSetOptions.__init__(self,spr_sets_dict,bdb_only)
		self.image_count = False

class E2Refine2DChooseParticlesTask(ParticleWorkFlowTask):
	documentation_string = "Select the type of particles you wish to run e2refine2d on. Typically the phase_flipped_hp particles will give the \
 best results, followed by the phase_flipped particles. Wiener filtered particles have strong contrast, but tend to produce blurry class-averages, \
 since it is the particles that are filtered, not the averages. After you hit OK, you will be prompted for specific refinement options." 
	def __init__(self):
		ParticleWorkFlowTask.__init__(self)
		self.window_title = "Choose Particles For e2refine2d"
		self.preferred_size = (480,300)
		self.form_db_name ="bdb:emform.e2refine2d"
	
	def get_params(self):
		ptcl_opts = EMParticleOptions()
		self.particles_map, self.particles_name_map,choices,unused = ptcl_opts.get_particle_options()
		choices.append("Specify")
			
		params = []		
		#if as_string:
		#params.append(ParamDef(name="particle_set_choice",vartype="string",desc_long="Choose the particle data set you wish to use to generate a starting data for e2refine2d",desc_short=title,property=None,defaultunits=db.get("particle_set_choice",dfl=""),choices=choices))
		db = db_open_dict(self.form_db_name)
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2Refine2DChooseParticlesTask.documentation_string,choices=None))
		params.append(ParamDef(name="particle_set_choice",vartype="choice",desc_long="Choose the particle data set you wish to use to generate a starting data for e2refine2d",desc_short="Choose data",property=None,defaultunits=db.get("particle_set_choice",dfl=""),choices=choices))
		return params

	def on_form_ok(self,params):
		if not params.has_key("particle_set_choice") or params["particle_set_choice"] == None:
			error("Please choose some data")
			return
		
		choice = params["particle_set_choice"]
		
		if choice == "Specify":
			self.emit(QtCore.SIGNAL("replace_task"),E2Refine2DWithGenericTask(),"e2refine2d arguments")
		else:
			self.emit(QtCore.SIGNAL("replace_task"),E2Refine2DRunTask(self.particles_map[self.particles_name_map[choice]]),"e2refine2d arguments")
		self.form.close()
		self.form = None
		
		self.write_db_entries(params)


		 	
class E2Refine2DChooseSetsTask(ParticleWorkFlowTask):
	documentation_string = "Choose the data you wish to use for use for running e2refine2d from the list of options below and hit OK. This will pop up a second form asking you to fill in more details.\n\nNote that usually you should have 4 options to choose from below. If you are not seeing all 4 options it means you should go back in the work flow, import particles, and generate phase flipped and Wiener filtered output." 
	def __init__(self):
		ParticleWorkFlowTask.__init__(self)
		self.window_title = "Choose Sets For e2refine2d"
		self.preferred_size = (480,300)
		self.form_db_name ="bdb:emform.e2refine"
	
	def get_params(self):
		ptcl_opts = EMSetsOptions(bdb_only=True)
		self.particles_map, self.particles_name_map,choices,unused = ptcl_opts.get_particle_options()
		choices.append("Specify")
			
		params = []
		
		#if as_string:
		#params.append(ParamDef(name="particle_set_choice",vartype="string",desc_long="Choose the particle data set you wish to use to generate a starting data for e2refine2d",desc_short=title,property=None,defaultunits=db.get("particle_set_choice",dfl=""),choices=choices))
		db = db_open_dict(self.form_db_name)
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2Refine2DChooseParticlesTask.documentation_string,choices=None))
		params.append(ParamDef(name="particle_set_choice_stack",vartype="choice",desc_long="Choose the particle data set you wish to use to generate a starting data for e2refine2d",desc_short="Choose data",property=None,defaultunits=db.get("particle_set_choice_stack",dfl=""),choices=choices))
		return params

	def on_form_ok(self,params):
		if not params.has_key("particle_set_choice_stack") or params["particle_set_choice_stack"] == None:
			error("Please choose some data")
			return
		choice = params["particle_set_choice_stack"]
		
		if choice == "Specify":
			self.emit(QtCore.SIGNAL("replace_task"),E2Refine2DWithGenericTask(),"e2refine2d arguments")
		else:
			task = E2Refine2DRunTask(self.particles_map[self.particles_name_map[choice]])
			task.single_selection = True
			self.emit(QtCore.SIGNAL("replace_task"),task,"e2refine2d arguments")
		self.form.close()
		self.form = None
		
		self.write_db_entries(params)
		
class E2RefFreeClassAveTool:
	def __init__(self):
		self.project_data_at_init = None
	
	def get_reference_free_class_averages_table(self):
		data_dict = EMProjectDataDict(spr_rfca_dict)
		data = data_dict.get_data_dict()
		self.project_data_at_init = data # so if the user hits cancel this can be reset
		names = data.keys()
		from emform import EM2DStackTable,EMFileTable,int_lt
		table = EM2DStackTable(names,desc_short="Class Averages",desc_long="")
		context_menu_data = EMRawDataReportTask.ProjectListContextMenu(spr_rfca_dict)
		table.add_context_menu_data(context_menu_data)
		table.add_button_data(EMRawDataReportTask.ProjectAddRawDataButton(table,context_menu_data))
		table.insert_column_data(1,EMFileTable.EMColumnData("Particles On Disk",EMParticleReportTask.get_num_ptcls,"Particles currently stored on disk that are associated with this image",int_lt))
		table.insert_column_data(2,EMFileTable.EMColumnData("Particle Dims",EMParticleReportTask.get_particle_dims,"The dimensions of the particles that are stored on disk"))
		
		return table, len(names)
	
	def recover_original_raw_data_list(self):
		'''
		Called if the user hits cancel - if they removed some files or added files the changes
		are not saved unless the user hits ok
		'''
		project_db = db_open_dict("bdb:project")
		project_db[spr_rfca_dict] = self.project_data_at_init
		

class E2Refine2DRunTask(E2Refine2DTask):
	documentation_string = "Select a particle set to use for 2-D refinement.  Select how many classes you want to generate. The other options on this page will generally produce \
reasonably good results with default settings. Normalize projection vectors can make a dramatic improvement on some data sets, but is bad for others. See the e2refine2d.py \
documentation for other details. Don't forget the other tabs !\
	"
	warning_string = "\n\nThere are no files" 
	def __init__(self,particles,single_selection = False):
		E2Refine2DTask.__init__(self)
		self.window_title = "Run e2refine2d"
		self.end_tag = "_ptcls"
		self.particles = particles
		self.single_selection = single_selection
		
	
	def run_form(self):
		self.form = EMTableFormWidget(self.get_params())
		self.form.resize(*self.preferred_size)
		self.form.setWindowTitle(self.window_title)
		get_application().show_specific(self.form)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_ok"),self.on_form_ok)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_cancel"),self.on_form_cancel)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_close"),self.on_form_close)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("display_file"),self.on_display_file)
		
	def get_params(self):
	 	params = []
		
		params.append(self.get_main_params())
		params.append(self.get_simmx_page())
		params.append(self.get_classaverage_page(include_sep=False))
		
		return params
	
	def get_main_params(self):
		params = []
		
#		if self.end_tag != "generic":
#			p,n = self.get_particle_selection_table(tag=self.end_tag)
#		else:
#			p = ParamDef(name="filenames",vartype="url",desc_short="Input file name(s)",desc_long="The names of the particle files you want to use as in the input data for e2refine2d.py",property=None,defaultunits=[],choices=[])
#			n = 1 # just to fool the next bit, that's all
#
	   	p,n = self.get_particle_selection_table(self.particles,None,self.single_selection,enable_ctf=False)
#		if n == 0:
#			params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2Refine2DRunTask.documentation_string+E2Refine2DRunTask.warning_string,choices=None))
#		else:
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2Refine2DRunTask.documentation_string,choices=None))
		params.append(p)
	
		other_params = self.get_general_params()
		params.extend(other_params)

		return ["General",params]

	def on_form_ok(self,params):
		self.write_db_entries(params) # I wrote the entries here so that the filenames entry is not altered, allowing the form parameters to be correctly memorized. This breaks a general pattern I am following. It probably doesn't matter at all.
		
#		if self.end_tag != "generic":
#			names = ["bdb:particles#"+name for name in params["filenames"]]
#			params["filenames"] = names 
			
		options = EmptyObject()
		for checker in [self.check_classaverage_page,self.check_simmx_page,self.check_main_page]: # check main page needs the shrink parameter to be checked first
			error_message = checker(params,options)
			if len(error_message) > 0 :
				self.show_error_message(error_message)
#				self.display_errors(error_message)
				return
		# w'oh if we make it here a lot of checking has occured. Now get the args in order to spawn_single_task
		string_args = []
		bool_args = []
		additional_args = []
		
		for get_args in [self.add_main_args,self.add_simmx_args,self.add_classaverage_args]:
			if get_args == self.add_classaverage_args:
				error = get_args(options,string_args,bool_args,additional_args,include_sep=False)
			elif get_args == self.add_simmx_args:
				error = get_args(options,string_args,bool_args,additional_args,include_shrink=False)
			else:
				error = get_args(options,string_args,bool_args,additional_args)
		  	#error = get_args(options,string_args,bool_args,additional_args)
		
			if error != None: # not too fast, something still could have gone wrong
				#self.display_errors([error])
				self.show_error_message([error])
				return
		if self.workflow_setting: 
			additional_args = ["--dbls=%s" %spr_rfca_dict]
		else: 
			additional_args = []
		temp_file_name = "e2refine2d_stdout.txt"
		self.spawn_single_task("e2refine2d.py",options,string_args,bool_args,additional_args,temp_file_name)
		self.emit(QtCore.SIGNAL("task_idle"))
		self.form.close()
		self.form = None

class E2Refine2DWithGenericTask(E2Refine2DRunTask):
	def __init__(self,workflow_setting=True):
		E2Refine2DRunTask.__init__(self,[])
		self.end_tag = "generic"
		self.workflow_setting = workflow_setting
			
							
class E2InitialModelsTool:
	def __init__(self):
		self.project_data_at_init = None
		
	def get_initial_models_table(self, makebutton=1):
		data_dict = EMProjectDataDict(spr_init_models_dict)
		init_model_data = data_dict.get_data_dict()
		self.project_data_at_init = init_model_data # so if the user hits cancel this can be reset
		init_model_names = init_model_data.keys()

		from emform import EM3DFileTable,EMFileTable,float_lt
		table = EM3DFileTable(init_model_names,name="model",desc_short="Initial Models",desc_long="")
		context_menu_data = EMRawDataReportTask.ProjectListContextMenu(spr_init_models_dict)
		table.add_context_menu_data(context_menu_data)
		if makebutton:
			table.add_button_data(EMRawDataReportTask.ProjectAddRawDataButton(table,context_menu_data))
		table.add_column_data(EMFileTable.EMColumnData("Quality",E2InitialModelsTool.get_quality_score,"This the quality score as determined by e2initialmodel.py",float_lt))

		table.add_column_data(EMFileTable.EMColumnData("Dimensions",EMRawDataReportTask.get_image_dimensions,"The dimensions of the file on disk"))
		#p.append(pdims) # don't think this is really necessary
		return table,len(init_model_names)
	
	def get_quality_score(image_name):
		'''
		Used by the initial models table to get a quality score
		'''
		a = EMData()
		a.read_image(image_name,0,True)
		d = a.get_attr_dict()
		if d.has_key("quality"): return "%.3f" %(d["quality"])
		else: return "-"
	
	get_quality_score = staticmethod(get_quality_score)
	
	def recover_original_raw_data_list(self):
		'''
		Called if the user hits cancel - if they removed some files or added files the changes
		are not saved unless the user hits ok
		'''
		project_db = db_open_dict("bdb:project")
		project_db[spr_init_models_dict] = self.project_data_at_init
	
class EMInitialModelReportTask(ParticleWorkFlowTask):
	"""This form displays the initial models currently associated with the project. You can associate initial models with the project using e2makeinitialmodel or by importing them directly, see the options below."""
	warning_string = "\n\n\nNOTE: There are no initial models currently associated with the project."
	def __init__(self):
		ParticleWorkFlowTask.__init__(self)
		self.window_title = "Project Initial Models"
		self.imt = None # will eventually become an E2IntialModelsTool
	
	def get_params(self):
		params = []
		
		self.imt = E2InitialModelsTool()
		p,n = self.imt.get_initial_models_table()
		
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=self.__doc__,choices=None))
		params.append(p)
		return params

	def on_form_cancel(self):
		self.imt.recover_original_raw_data_list()
		
		self.form.close()
		self.form = None
		self.emit(QtCore.SIGNAL("task_idle"))
		

class E2InitialModel(ParticleWorkFlowTask):
	"""Use this form for creating initial model with e2initialmodel."""
	def __init__(self):
		ParticleWorkFlowTask.__init__(self)
		self.window_title = "Run e2initialmodel"
		self.form_db_name = "bdb:emform.e2initialmodel"
		
	def get_params(self):
		params = []
		
		self.rfcat = E2RefFreeClassAveTool()
		p,n = self.rfcat.get_reference_free_class_averages_table()

		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=self.__doc__,choices=None))
		
		db = db_open_dict(self.form_db_name)
		piter = ParamDef(name="iter",vartype="int",desc_short="Iterations",desc_long="The number of times each 3D is iteratively refined",property=None,defaultunits=db.get("iter",dfl=8),choices=[])
		ptries = ParamDef(name="tries",vartype="int",desc_short="Tries",desc_long="The number of 3D models to generate",property=None,defaultunits=db.get("tries",dfl=10),choices=[])
		syms = ["icos","oct","tet","c","d","h"]
		psym =  ParamDef(name="symname",vartype="string",desc_short="Symmetry",desc_long="Symmetry to be imposed during refinement",property=None,defaultunits=db.get("symname",dfl="c"),choices=syms)
		psymnum = ParamDef(name="symnumber",vartype="string",desc_short="Symmetry number",desc_long="In C,D and H symmetry, this is the symmetry number",property=None,defaultunits=db.get("symnumber",dfl=""),choices=None)
		#db_close_dict(self.form_db_name)
		orientgens = self.get_orientgens_list()
		porientgens =  ParamDef(name="orientgen",vartype="string",desc_short="Orientation generator",desc_long="The method of orientation generation",property=None,defaultunits=db.get("orientgen",dfl="eman"),choices=orientgens)

		p.enable_multiple_selection = False
		params.append(p)
		params.append([piter,ptries,porientgens])
		params.append([psym,psymnum])
	
		return params
	
	def on_form_ok(self,params):
		
		error_message = []
		
		if not params.has_key("filenames"):
			# THERE ARE no classes to choose from and the user has hit ok
			self.emit(QtCore.SIGNAL("task_idle"))
			self.form.close()
			self.form = None
			return
		
		error_message = []
		if len(params["filenames"]) != 1: error_message.append("Please choose a single file to proceed")
		
		if params["iter"] < 0: error_message.append("Please specify an iter value of atleast 0")
		if params["tries"] < 1: error_message.append("Please specify a number of tries that is atleast 1")
		
		
		# copied from e2refine
		
		options = EmptyObject()
		
		error_message.extend(self.check_sym(params,options))
		options.orientgen = params["orientgen"]
		
		if len(error_message) > 0:
			self.show_error_message(error_message)
			return	
		else:
			# obviously if we make it here the parameters are sensible so we can store them (and later use them again for the user's convenience)
			self.write_db_entries(params)
			
			options.input = params["filenames"][0]
			options.iter = params["iter"]
			options.tries = params["tries"]
			options.filenames = []
			#options.sym - taken care of by check_sum
			string_args = ["iter","input","tries","sym","orientgen"]
			bool_args = []
			additional_args = ["--dbls=%s" %spr_init_models_dict]
			temp_file_name = "e2initialmodel_stdout.txt"
			self.spawn_single_task("e2initialmodel.py",options,string_args,bool_args,additional_args,temp_file_name)
			self.emit(QtCore.SIGNAL("task_idle"))
			self.form.close()
			self.form = None
	def on_form_cancel(self):
		self.rfcat.recover_original_raw_data_list()
		
		self.form.close()
		self.form = None
		self.emit(QtCore.SIGNAL("task_idle"))

class RefinementReportTask(ParticleWorkFlowTask):
	documentation_string = "This form displays the models produced at the end of each refinement."
	warning_string = "\n\n\nNOTE: There are no results available."
	def __init__(self):
		ParticleWorkFlowTask.__init__(self)
		self.window_title = "Refinement Results"
	
	def get_params(self):
		params = []
		
		p,n = self.get_last_refinement_models_table()
		
		if n == 0:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=RefinementReportTask.documentation_string+RefinementReportTask.warning_string,choices=None))
		else:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=RefinementReportTask.documentation_string,choices=None))
			params.append(p)
		return params
	
	def get_rmodels_list(self):
		'''
		Get the most recent refinement models in each refinement directory
		Has a preference for listing the threed_filt_?? models as opposed
		to the threed_?? models, seeing as they are generally better.
		@return a list of strings
		'''
#		dirs = get_numbered_directories("refine_")
		dirs = get_prefixed_directories("refine_")
		threed_files = []
		for dir in dirs:
			threed_db = None
			# check for 00 to 09 but 00 is replaced with "init"
			db_first_part = "bdb:"+dir+"#threed_"
			cont = True
			threed_file = None
			for i in range(0,10):
				if not cont: break
				for j in range(0,10):
					if not cont: break
					
					db_name_filt = db_first_part+"filt_"+str(i)+str(j)
					if db_check_dict(db_name_filt):
						threed_file = db_name_filt
					else:
						db_name = db_first_part+str(i)+str(j)
						if db_check_dict(db_name):
							threed_file = db_name
						else:
							cont = False
							
			if threed_file != None:
				threed_files.append(threed_file)
				
		return threed_files
	
	def get_last_refinement_models_table(self):

		from emform import EM3DFileTable,EMFileTable
		names = self.get_rmodels_list()
		table = EM3DFileTable(names,name="model",desc_short="Most Recent Models",desc_long="")

#		table.add_column_data(EMFileTable.EMColumnData("Mean",EMFilterRawDataTask.get_mean,"The mean pixel value"))
#		table.add_column_data(EMFileTable.EMColumnData("Sigma",EMFilterRawDataTask.get_sigma,"The standard deviation of the pixel values"))
		table.add_column_data(EMFileTable.EMColumnData("Dimensions",EMRawDataReportTask.get_image_dimensions,"The dimensions of the file on disk"))

		return table,len(names)
	
	
			
class E2Make3DTools:
	'''
	e2eotest and e2refine tasks both need the functionality embodied here
	'''
	
	preprocessor_cache = None
	
	def __init__(self):pass
	
	def add_make3d_args(self,options,string_args,bool_args,additional_args):
		
		string_args.extend(["m3diter","m3dkeep","recon"])
		bool_args.append("m3dkeepsig")
		bool_args.append("m3dsetsf")
		if hasattr(options,"m3dpreprocess"): string_args.append("m3dpreprocess")
		if hasattr(options,"m3dpostprocess"): string_args.append("m3dpostprocess")
		if hasattr(options,"m3dpostprocess2"): string_args.append("m3dpostprocess2")
		if hasattr(options,"pad"): string_args.append("pad")

	def check_make3d_page(self,params,options):
		error_message = []
		
		if params["m3diter"] < 0:
			error_message.append("The number of make3d iterations must be atleast 0")
		
		if params["m3dkeepsig"] == False:
			if params["m3dkeep"] < 0 or params["m3dkeep"] > 1:
				error_message.append("The keep parameter in the Make3D page must be between 0 and 1. This does not hold if the \'Sigma based\' option is selected.")
		
		if len(error_message) > 0 : return error_message # calling program should discontinue
		
		
		
		if len(params["pad"]) > 0:
			try: int(params["pad"])
			except: error_message.append("The value you entered for padding is nonsensical")
			pad = int(params["pad"])
			options.pad = int(params["pad"])
			if params.has_key("filenames") and params["filenames"] > 0:
				nx,ny = gimme_image_dimensions2D(params["filenames"][0])
				if nx >= pad or ny >= pad:
					error_message.append("You must specify a value for padding that is larger than the image dimensions - the image dimensions are %i x %i and your pad value is %i" %(nx,ny,pad))				
#				else:
#					options.pad = int(params["pad"])
#				except:
#					error_message.append("Can't get the dimensions of the first image???")
			else:
				pass # the user not entering filenames is an error, so after they've correct that we'll with the issues here
		
		options.m3diter = params["m3diter"]
		options.m3dkeep = params["m3dkeep"]
		options.m3dkeepsig = params["m3dkeepsig"]
		options.m3dsetsf = params["m3dsetsf"]
		
		options.recon = params["recon"]
		
		if params["m3dpreprocess"] != "None":
			options.m3dpreprocess = params["m3dpreprocess"]
		
		if options.m3dsetsf :
			try:
				db_misc=db_open_dict("bdb:e2ctf.misc",True)
				m=db_misc["strucfac"]
			except:
				error_message.append("You must determine a structure factor before you can use the setsf option")
		
		if params["m3dpostprocess"] != "None":
			if len(params["m3dpostprocessargs"]) == 0:
				error_message.append("If you supply a post processor for make3d, you have to supply one of the cutoff_abs, cutoff_freq, or cutoff_pixels parameters")
			else:
				s = params["m3dpostprocess"] + ":"
				s += params["m3dpostprocessargs"]
				p = parsemodopt(s)
				if p[0] == None:
					error_message.append("Error can't interpret the make3d post processor string (%s)" %(s))
				else:
					try:
						Processors.get(p[0], p[1])
						options.m3dpostprocess = s
					except:
						error_message.append("Error, can't interpret parameters for the make3d post processor (%s)" %(s))
						vals = dump_processors_list()
						values = vals[p[0]]
						s = "The parameters for the %s processor are:"  %p[0]
						
						for i in xrange(1,len(values),3):
							s += " " + values[i] +","
						s = s[:-1] # get rid of the last column
						error_message.append(s)
		return error_message

		if params["m3dpostprocess2"] != "None":
			if len(params["m3dpostprocessargs2"]) == 0:
				error_message.append("If you supply a second post processor for make3d, you have to supply one of the cutoff_abs, cutoff_freq, or cutoff_pixels parameters")
			else:
				s = params["m3dpostprocess2"] + ":"
				s += params["m3dpostprocessargs2"]
				p = parsemodopt(s)
				if p[0] == None:
					error_message.append("Error can't interpret the make3d post processor string (%s)" %(s))
				else:
					try:
						Processors.get(p[0], p[1])
						options.m3dpostprocess2 = s
					except:
						error_message.append("Error, can't interpret parameters for the make3d post processor (%s)" %(s))
						vals = dump_processors_list()
						values = vals[p[0]]
						s = "The parameters for the %s processor are:"  %p[0]
						
						for i in xrange(1,len(values),3):
							s += " " + values[i] +","
						s = s[:-1] # get rid of the last column
						error_message.append(s)
		return error_message

	def get_make3d_page(self):
		
		db = db_open_dict(self.form_db_name)
		pkeep = ParamDef(name="m3dkeep",vartype="float",desc_short="keep",desc_long="The fraction of particles to keep in each class average. If sigma based is checked this value is interpreted in standard deviations from the mean instead",property=None,defaultunits=db.get("m3dkeep",dfl=0.8),choices=[])
		pkeepsig = ParamDef(name="m3dkeepsig",vartype="boolean",desc_short="Sigma based",desc_long="If checked the keep value is interpreted in standard deviations from the mean instead of basic ratio",property=None,defaultunits=db.get("m3dkeepsig",dfl=True),choices=[])
		
		piter = ParamDef(name="m3diter",vartype="int",desc_short="Reconstruction iterations",desc_long="The number of times the reconstruction algorithm is iterated",property=None,defaultunits=db.get("m3diter",dfl=3),choices=[])

		pkeepsig = ParamDef(name="m3dkeepsig",vartype="boolean",desc_short="Sigma based",desc_long="If checked the keep value is interpreted in standard deviations from the mean instead of basic ratio",property=None,defaultunits=db.get("m3dkeepsig",dfl=True),choices=[])

		pnormproc =  ParamDef("m3dpreprocess",vartype="string",desc_short="Normalization processor",desc_long="The normalization method applied to the class averages",property=None,defaultunits=db.get("m3dpreprocess",dfl="normalize.edgemean"),choices=["normalize","normalize.edgemean","None"])
		
		psetsf = ParamDef(name="m3dsetsf",vartype="boolean",desc_short="Set SF",desc_long="If checked, this will impose the precomputed structure factor on the model before postprocessing",property=None,defaultunits=db.get("m3dsetsf",dfl=False),choices=[])
		
		ppostproc =  ParamDef("m3dpostprocess",vartype="string",desc_short="Post processor",desc_long="A post processor applied to the reconstructed model",property=None,defaultunits=db.get("m3dpostprocess",dfl="None"),choices=self.get_postprocess_filt_options())
		ppostprocargs =  ParamDef(name="m3dpostprocessargs",vartype="string",desc_short="params",desc_long="Parameters for the post processor see \"e2help.py processors\"",property=None,defaultunits=db.get("m3dpostprocessargs",dfl=""),choices=[])	

		ppostproc2 =  ParamDef("m3dpostprocess2",vartype="string",desc_short="Post processor 2",desc_long="A second post processor applied to the reconstructed model",property=None,defaultunits=db.get("m3dpostprocess2",dfl="None"),choices=self.get_postprocess_filt_options())
		ppostprocargs2 =  ParamDef(name="m3dpostprocessargs2",vartype="string",desc_short="params 2",desc_long="Parameters for the post processor see \"e2help.py processors\"",property=None,defaultunits=db.get("m3dpostprocessargs2",dfl=""),choices=[])	

		precon = ParamDef("recon",vartype="string",desc_short="Reconstructor",desc_long="The method used to perform 3D reconstruction",property=None,defaultunits=db.get("recon",dfl="fourier"),choices=["fourier","wiener_fourier","back_projection"])
		ppad = ParamDef("pad",vartype="string",desc_short="Pad to",desc_long="The amount to which you want to pad the 3D volume when Fourier inversion is being used. At least 25% is recommended", defaultunits=db.get("pad",dfl=""),choices=[])
		params = []
		
		#db_close_dict(self.form_db_name)
		
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2RefineParticlesTask.make3d_documentation,choices=None))

		
		params.append([precon,piter])
		params.append([pnormproc,ppad])
		params.append([pkeep,pkeepsig,psetsf])
		params.append([ppostproc,ppostprocargs])
#		params.append([ppostproc2,ppostprocargs2])
		return ["Make3d", params]
		
	
	def get_postprocess_filt_options(self):
		if E2Make3DTools.preprocessor_cache == None:
			a = dump_processors_list()
			l = ["None"]
			for key in a.keys():
				if len(key) > 5 and key[:6] == "filter":
					vals = key.split(".")
					if len(vals) > 1:
						if vals[1] in ["lowpass","highpass"]:
							l.append(key)
							
			E2Make3DTools.preprocessor_cache = l
			 
		return E2Make3DTools.preprocessor_cache
		
class E2RefineParticlesTaskBase(EMClassificationTools, E2Make3DTools):
	'''
	This class is a base class for eman refinements, to avoid massive code duplicantions
	'''

	def __init__(self,ptcls_list,usefilt_ptcls_list):
		'''
		@param ptcls_list the list of particle files that will form the primary input to e2refine
		@param usefilt_ptcls_list the list of usefilt particles corresponding in length to the ptcls_list, or None
		'''
		self.ptcls = ptcls_list
		self.usefilt_ptcls = usefilt_ptcls_list
		self.imt = None # will eventually become an E2IntialModelsTool
	 	EMClassificationTools.__init__(self)
	 	E2Make3DTools.__init__(self)
	 	
	 	self.window_title = "Run e2refine"
	 	self.form_db_name = "bdb:emform.e2refine"
	 	self.single_selection = False

	class UsefiltColumn:
		def __init__(self,ptcls,usefilt_ptcls):
			if len(ptcls) != len(usefilt_ptcls):
				raise RuntimeError("The usefilt and raw particle lists must be the same length")
			
			self.filt_map = {}
			for i in xrange(0,len(ptcls)):
				self.filt_map[ptcls[i]] = usefilt_ptcls[i]
				
		def get_usefilt_name(self,name):
			return self.filt_map[name]
		
	def run_form(self):
		self.form = EMTableFormWidget(self.get_params())
		self.form.resize(*self.preferred_size)
		self.form.setWindowTitle(self.window_title)
		get_application().show_specific(self.form)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_ok"),self.on_form_ok)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_cancel"),self.on_form_cancel)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_close"),self.on_form_close)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("display_file"),self.on_display_file)
		
	def get_params(self):
	 	params = []
		
	#	params.append(self.get_intro_params())
		params.append(self.get_main_params())
		params.append(self.get_main_params_2())
		params.append(self.get_project3d_page())
		params.append(self.get_simmx_page())
		params.append(self.get_classaverage_page())
		params.append(self.get_make3d_page())
		
		return params
	
	def on_form_cancel(self):
		if self.imt != None: self.imt.recover_original_raw_data_list()
		
		self.form.close()
		self.form = None
		self.emit(QtCore.SIGNAL("task_idle"))
	
	def on_form_ok(self,params):
		
		options = EmptyObject()
#		error_message = self.check_main_page(params,options)

		for checker in [self.check_main_page,self.check_project3d_page,self.check_simmx_page,self.check_classaverage_page,self.check_make3d_page]:
			error_message = checker(params,options)
			if len(error_message) > 0 :
				self.display_errors(error_message)
				return	
		self.write_db_entries(params)
		# w'oh if we make it here a lot of checking has occured. Now get the args in order to spawn_single_task
		string_args = []
		bool_args = []
		
		additional_args = []
		
		for get_args in [self.add_general_args,self.add_project3d_args,self.add_simmx_args,self.add_classaverage_args,self.add_make3d_args]:
		  	error = get_args(options,string_args,bool_args,additional_args)
		
			if error != None: # not too fast, something still could have gone wrong
				self.display_errors([error])
				return
		
		temp_file_name = "e2refine_stdout.txt"
		
		# Steve is rethinking how we remember programs arguments
		#self.write_db_parms(options,string_args,bool_args)

		self.spawn_single_task("e2refine.py",options,string_args,bool_args,additional_args,temp_file_name)
		
		self.emit(QtCore.SIGNAL("task_idle"))
		self.form.close()
		self.form = None

# This functionality is being redesigned and pends a discussion with Steve ludtke with respect to the history mechanism
#	def write_db_parms(self,options,string_args,bool_args):
#		db = db_open_dict("bdb:e2refine.args")
#		
#		for string in string_args:
#			db[string] = getattr(options,string)
#			
#		for string in bool_args:
#			db[string] = getattr(options,string)
#			
#		db_close_dict("bdb:e2refine.args")
#		
	
	def display_errors(self,error_message):
		'''
		error_message is a list of strings
		'''
		
		if len(error_message) > 0:
			self.show_error_message(error_message)
	 		
	
	def get_main_params(self):
		'''
		General/broad refine params
		'''
		params = []
#		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits="Choose the particles you wish to refine",choices=None))

		
		if self.ptcls == None:
			p,n = self.get_particle_selection_table([],None,self.single_selection,enable_ctf=False)
		else:
			p,n = self.get_particle_selection_table(self.ptcls,None,self.single_selection,enable_ctf=False)
			if self.usefilt_ptcls != None and len(self.usefilt_ptcls) > 0:
				from emform import EMFileTable
				self.column_data = E2RefineParticlesTask.UsefiltColumn(self.ptcls,self.usefilt_ptcls)
				p.add_column_data(EMFileTable.EMColumnData("Usefilt data",self.column_data.get_usefilt_name,"The usefilt data"))
				
#			p = ParamDef(name="filenames",vartype="url",desc_short="Input file name(s)",desc_long="The names of the particle files you want to use as in the input data for e2refine2d.py",property=None,defaultunits=[],choices=[])
#			n = 1
		
		# I could check to see if the database exists but it seems unnecessary
		# In the event that the database doesn't exist it is created and 
		# a new entry is created on disk. The only inconvenient aspect of this comes
		# if the user hits cancel - then there is a file on disk even though
		# the user never agreed to anything
		db = db_open_dict(self.form_db_name) # see eman wiki for a list of what args are kept in this db
		project_db = db_open_dict("bdb:project")
		
		params.append(p)
		
	   	pmass = ParamDef(name="global.particle_mass",vartype="float",desc_short="Particle mass (kda)",desc_long="The mass of the particle in kilodaltons. Leave blank if unknown",property=None,defaultunits=project_db.get("global.particle_mass",dfl=800),choices=None)
		papix = ParamDef(name="global.apix",vartype="float",desc_short="Angtsrom per pixel",desc_long="The physical distance represented by the pixel spacing",property=None,defaultunits=project_db.get("global.apix",dfl=1.1),choices=None)
		
		params.append([papix,pmass])
		
		piter = ParamDef(name="iter",vartype="int",desc_short="Refinement iterations",desc_long="The number of times 3D refinement should be iterated",property=None,defaultunits=db.get("iter",dfl=3),choices=[])
		plowmem = ParamDef(name="lowmem",vartype="boolean",desc_short="Low mem",desc_long="Causes various programs to restrict memory usage but results in increased CPU time.",property=None,defaultunits=db.get("lowmem",dfl=False),choices=None)

	   	params.append([piter,plowmem])
	   	
	   	pparallel = ParamDef(name="parallel",vartype="string",desc_short="Parallel",desc_long="Parallel arguments (advanced). Leave blank if unsure",property=None,defaultunits=db.get("parallel",dfl=""),choices=None)
		
		params.append(pparallel)
	
		return ["Particles",params]
		 
	
	def get_main_params_2(self):
		
		params = []
#		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits="",choices=None))

		# I could check to see if the database exists but it seems unnecessary
		# In the event that the database doesn't exist it is created and 
		# a new entry is created on disk. The only inconvenient aspect of this comes
		# if the user hits cancel - then there is a file on disk even though
		# the user never agreed to anything
		db = db_open_dict(self.form_db_name) # see eman wiki for a list of what args are kept in this db
		project_db = db_open_dict("bdb:project")
		
		self.imt = E2InitialModelsTool()
		p1,n1 = self.imt.get_initial_models_table()
		p1.enable_multiple_selection = False
		params.append(p1)
			
		pautomask = ParamDef(name="automask3d",vartype="boolean",desc_short="Auto mask 3D",desc_long="Causes automasking of the 3D volume to occur at the end of each iteration",property=None,defaultunits=db.get("automask3d",dfl=True),choices=None)
		
		params.append(pautomask)
		
		def_rad = 30
		def_mask_dltn = 5
		if len(self.ptcls) > 0:
			nx,ny,nz = gimme_image_dimensions3D(self.ptcls[0])
			def_rad = int(nx/8.0)
			def_mask_dltn = int(nx/20.0)
		
		pamthreshold =  ParamDef(name="amthreshold",vartype="float",desc_short="Threshold",desc_long="An isosurface threshold that well defines your structure.",property=None,defaultunits=db.get("amthreshold",dfl=0.7),choices=None)
		pamradius =  ParamDef(name="amradius",vartype="int",desc_short="Radius",desc_long="The radius of a sphere at the the origin which contains seeding points for the flood file operation using the given threshold",property=None,defaultunits=db.get("amradius",dfl=def_rad),choices=None)
		pamnmax =  ParamDef(name="amnmax",vartype="int",desc_short="NMax",desc_long="Number of highest value voxels to use as a seed",property=None,defaultunits=db.get("amradius",dfl=def_rad),choices=None)
		pamnshells =  ParamDef(name="amnshells",vartype="int",desc_short="Mask dilations",desc_long="The number of dilations to apply to the mask after the flood fill operation has finished. Suggest 5% of the boxsize",property=None,defaultunits=db.get("amnshells",dfl=def_mask_dltn),choices=None)
		pamngaussshells =  ParamDef(name="amnshellsgauss",vartype="int",desc_short="Post Gaussian dilations",desc_long="The number of dilations to apply to the dilated mask, using a gaussian fall off. Suggest 5% of the boxsize",property=None,defaultunits=db.get("amnshellsgauss",dfl=def_mask_dltn),choices=None)
		
		pautomask.dependents = ["amthreshold","amradius","amnshells","amnshellsgauss"] # these are things that become disabled when the pautomask checkbox is checked etc
		
		params.append([pamthreshold,pamnmax,pamradius])
		params.append([pamnshells,pamngaussshells])

		#db_close_dict(self.form_db_name)
		#db_close_dict("bdb:project")
		
		return ["Model",params]

	def add_general_args(self,options,string_args,bool_args,additional_args):
		
		options.path = numbered_path("refine",True)
		
		if len(options.filenames) > 1 :
			if options.filenames[0][:4] == "bdb:":
				success,cmd = self.make_v_stack(options.filenames,"all",options,"input")
				if not success:
					return cmd + " failed"
			else:
				success, cmd = self.make_stack(options.filenames, "all",options,"input")
				if not success:
					return cmd + " failed"
		else:
			options.input = options.filenames[0]
			
		string_args.append("input")
		
		if hasattr(options,"usefilt_names"):
			if len(options.usefilt_names) > 1 :
				if options.usefilt_names[0][:4] == "bdb:":
						success,cmd = self.make_v_stack(options.usefilt_names,"usefilt",options,"usefilt")# sets the attribute for us
						if not success:
							return cmd + " failed"
				else:
					success, cmd = self.make_stack(options.usefilt_names,"usefilt",options,"usefilt")
					if not success:
						return cmd + " failed"
			else:
				options.usefilt = options.usefilt_names[0]
			
			string_args.append("usefilt")
		
		if hasattr(options,"parallel"):string_args.append("parallel")
		
		error = self.check_model(options)
		
		if error != None:
			return error
		
		opt_attr = ["mass","apix","automask3d"] # these does not necessarily have to be specified
		for attr in opt_attr:
			if hasattr(options,attr): string_args.append(attr) 
		
		options.filenames = [] # important for this to happen so the argument doesn't have all the filenames as args
		string_args.extend(["iter","sym","model","path"])
		bool_args.append("lowmem")
		
		return None # returning None is good
	
	def check_model(self,options):
		
		model = options.model[0] # options.model is a list
		if not file_exists(model): # why did I do this? Oh well doesn't hurt # Retrospective note - it was useful as I did Steve's second round of alterations
			return "the initial model %s does not exist" %model
		
		nx,ny = gimme_image_dimensions2D(options.input)
		if nx != ny:
			return "input images aren't square"
		
		
		x,y,z = gimme_image_dimensions3D(model)
		
		
		if x != y or z != y:
			return "initial model isn't square"
		
		if nx != x:
			scale = float(nx)/x
			new_model = "bdb:"+options.path + "#initial_model"
			
			image = EMData()
			image.read_image(model,0)
			start = (x-nx)/2
			if scale > 1:
				image.clip_inplace(Region(start,start,start,nx,nx,nx))
				t = Transform()
				t.set_scale(scale)
				image.transform(t)	
			else:
				t = Transform()
				t.set_scale(scale)
				image.transform(t)
				image.clip_inplace(Region(start,start,start,nx,nx,nx))
				
			image.write_image(new_model,0) # db got opened here

		 	options.model = new_model
		 	db_close_dict(new_model) # force synchronization so e2refine.py will definetely run -
		else:
			options.model = model # all good
			
		return None
		
		
	def make_stack(self,filenames,out_name,options,attr):
		'''
		This one's a bit more specialized to handle flat files and avoid massive copying
		'''
		if len(filenames) == 1:
			setattr(options,attr,filenames[0])
			return True,""
		else:
			fail = False
			# check if they're all bdb files, in which case we can make a v stack
			for name in filenames:
				if name[0:4] != "bdb:":
					fail = True
					break
				
			if fail: # we can't make a vstack
				# potentially lots of e2proc2d
				progress = QtGui.QProgressDialog("Importing files into database...", "Abort import", 0, len(filenames),None)
				progress.show()
	
		  	   	i = 0
		  	   	setattr(options,attr, "bdb:"+options.path+"#"+out_name)
				for i,name in enumerate(filenames):
					cmd = "e2proc2d.py"
		 			cmd += " "+name
		 			cmd += " "+getattr(options,attr)
		 			success = (os.system(cmd) in (0,12))
		 			if not success:
		 				progress.close()
		 				return False,cmd
		 			else:
		 				progress.setValue(i+1)
		 				get_application().processEvents()
				
				progress.close()
				
				return True,cmd
			else:
				# This actually never happens the code needs changing
				raise NotImplementedException

	
	def make_v_stack(self,filenames,out_name,options,attr):
	 	
	 	cmd = "e2bdb.py"
	 	for name in filenames:
	 		cmd += " "+name
	 	
	 	cmd += " --makevstack=bdb:"+options.path+"#"+out_name
	 	
	 	print "executing cmd", cmd
	 	
	 	get_application().setOverrideCursor(Qt.BusyCursor)
	 	success = os.system(cmd)
	 	success = (success in (0,11,12))
	 	get_application().setOverrideCursor(Qt.ArrowCursor)
	 	
	 	setattr(options,attr,"bdb:"+options.path+"#"+out_name) # Note important 
	 	
	 	return success,cmd
	 
	
	def check_main_page(self,params,options):
		'''
		Called internally to check that the user has entered correct parameters in the main page
		returns a potentially empty list of error messages, if it is empty it means there are no errors
		Also sets argument attributes of the options object, killing two birds with one stone
		'''
		error_message = []
		#filenames
		if len(params["filenames"]) == 0:
			error_message.append("Please choose files to form the input data set.")
			
		if len(params["model"]) == 0:
			error_message.append("Please choose a starting model.")
			
		options.filenames = params["filenames"]
			#print options.filenames
		#usefilt
		
		if self.usefilt_ptcls != None and len(self.usefilt_ptcls) > 0:
			usefilt_names = [self.column_data.get_usefilt_name(name) for name in params["filenames"]]
			options.usefilt_names = usefilt_names

		
		if params.has_key("global.particle_mass"): 
			if params["global.particle_mass"] <= 0:
				error_message.append("The particle mass must be greater than 0")
			else:
				options.mass = params["global.particle_mass"]
			
		if params.has_key("global.apix"):
			if params["global.apix"] <= 0:
				error_message.append("The angstrom per pixel must be greater than  0")
			else:
				options.apix = params["global.apix"]
				
		if params.has_key("parallel") and len(params["parallel"]) > 0: 
			options.parallel = params["parallel"]
				
		if params["automask3d"]:
			# the user wants to do automasking
			names = ["amthreshold","amradius","amnshells","amnshellsgauss","amnmax"]
			arg = ""
			for i,name in enumerate(names):
				if not params.has_key(name):
					error_message.append("Missing automask parameter %s" %name[2:])
					continue
				#elif i == 1:
					#if params[name] <=0:
						#error_message.append("The automask radius parameter must be greater than 0")
						#continue
				elif i in [2,3]:
					if params[name] < 0:
						error_message.append("The automask dilation parameters must be atleast 0")
						continue
				# if we make it here than no error conditions were encountered, so we're safe to just append the argument
				if i != 0:
					arg +=","
				arg+= str(params[name])
		
			options.automask3d=arg
				
		#symmetry
		error_message.extend(self.check_sym(params,options))
		
		# iterations
		if params["iter"] < 1:
			error_message.append("The number of refinement iterations must be atleast 1.")
		else:
			options.iter = params["iter"]
			
		options.lowmem = params["lowmem"]
		options.model = params["model"] # can't get this one wrong
		
		return error_message
	
	def get_usefilt_options(self):
		if self.end_tag != "generic":
			
			n = self.get_total_particles(self.end_tag)
			
			available_filt_files = []
			number = []
			for i,tag in enumerate(self.usefilt_tags):
				if tag != self.end_tag:
					n_ = self.get_total_particles(tag=tag)
					if n_ > 0 and n_ == n:
						available_filt_files.append(self.usefilt_display_names[i])
					
			return available_filt_files
		else:
			return []
		
	def add_project3d_args(self,options,string_args,bool_args,additional_args):
		
		string_args.extend(["orientgen","projector"])
		# sym is already taken care of in the main args
		return None # no error to report
		
	def get_project3d_page(self):
		params = []
		
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2RefineParticlesTask.project3d_documentation,choices=None))

		db = db_open_dict(self.form_db_name)
		
		projectors = self.get_projectors_list()
		orientgens = self.get_orientgens_list()
			
		pprojector =  ParamDef(name="projector",vartype="string",desc_short="Projector",desc_long="The method used to generate projections",property=None,defaultunits=db.get("projector",dfl="standard"),choices=projectors)
		
		porientgens =  ParamDef(name="orientgen",vartype="string",desc_short="Orientation generator",desc_long="The method of orientation generation",property=None,defaultunits=db.get("orientgen",dfl="eman"),choices=orientgens)
		
		pmirror = ParamDef(name="incmirror",vartype="boolean",desc_short="Include mirror",desc_long="Include the mirror portion of the asymmetric uni",property=None,defaultunits=db.get("incmirror",False),choices=[])
		
		
		orient_options = ["Angle Based", "Number Based"]
		porientoptions = ParamDef(name="orientopt",vartype="choice",desc_short="Method of generating orientation distribution",desc_long="Choose whether you want the orientations generating based on an angle or based on a total number of orientations desired",property=None,defaultunits=db.get("orientopt",dfl=orient_options[0]),choices=orient_options)
		porientoptionsentry  =  ParamDef(name="orientopt_entry",vartype="float",desc_short="value",desc_long="Specify the value corresponding to your choice",property=None,defaultunits=db.get("orientopt_entry",dfl=5),choices=[])
		from emform import EMOrientationDistDialog
		buttondialog = EMOrientationDistDialog()
		params.append([pprojector,porientgens])
		
		syms = ["icos","oct","tet","d","c","h"]
		psym =  ParamDef(name="symname",vartype="string",desc_short="Symmetry",desc_long="Symmetry to be imposed during refinement",property=None,defaultunits=db.get("symname",dfl="c"),choices=syms)
		psymnum = ParamDef(name="symnumber",vartype="string",desc_short="Symmetry number",desc_long="In C,D and H symmetry, this is the symmetry number",property=None,defaultunits=db.get("symnumber",dfl="1"),choices=None)
		params.append([psym,psymnum])
			
		params.append([porientoptions,porientoptionsentry])
		params.append([pmirror,buttondialog])
		
		#db_close_dict(self.form_db_name)
		
		return ["Project 3D",params]
	
	def check_project3d_page(self,params,options):
		
		error_message = []
		if params["orientopt_entry"] < 0:
			error_message.append("Please enter a positive non zero value for the angle/number of projections in the Project3D settings")

		if params["orientgen"] == "rand" and params["orientopt"] == "Angle Based":
			error_message.append("The random orientation generator doesn't work with the \'Angle Based\' argument, please choose \'Number Based\' instead") 
		
		if int(params["orientopt_entry"]) !=  params["orientopt_entry"] and  params["orientopt"] == "Number Based":
			error_message.append("In project3d - for the Number Based orientation method the number must be an integer")
		
		options.orientgen = params["orientgen"]
		if params["orientopt"] == "Angle Based":
			options.orientgen += ":delta="
		else:
			options.orientgen += ":n="
		
		options.orientgen += str(params["orientopt_entry"])	
		
		if params["incmirror"]:
			options.orientgen += ":inc_mirror=1"
		else:
			options.orientgen += ":inc_mirror=0"
			
		
		options.projector = params["projector"]
		
		return error_message

class E2RefineParticlesTask(E2RefineParticlesTaskBase, EMClassificationTools, E2Make3DTools):
	'''
	This task will harness the parameters for, and launch, e2refine.py
	'''
	 
	general_documentation = "These are the general parameters for 3D refinement in EMAN2. Please select which particles you wish to use as part of this process, specify your starting model, and fill in other parameters such as symmetry and whether or not the usefilt option should be used."
	project3d_documentation = "These  parameters are used by e2project3d. Several orientation generation techniques provide alternative methods for distributing orientations in the asymmetric unit. Orientations can be generated based on your desired angular spacing, or alternatively on the desired total number of projections. In the latter case EMAN2 will generate a number as close as possible to the specified number, but note that there is no guarantee of a perfect match. You can also vary the method by which projections are generated. If you check the \'include mirror\' option you should be sure to use aligners to that do not perform mirror alignment."
	simmx_documentation = """These  parameters are used by e2simmx, a program that compares each particle to each projection and records quality scores. \
To do this the particles must first be aligned to the projections using the aligners you specify. Once aligned the \'Main comparator\' is used to record \
the quality score. These quality values are recorded to an image matrix on handed on to the next stage in the refinement process.
- The shrink parameter causes all projections and particles to be shrunken by the given amount prior to comparison. This can provide a significant time advantage, though at the expense
of resolution. Note however that the class averaging stage, which can involve iterative alignment, does not use shrunken data.
- 2 stage simmx is still experimental. If set to 2 instead of zero, classification will be performed in two stages resulting in a 5-25x speedup, but with a potential decrease in accuracy.
- PS match ref will force the power spectra of the particle and reference to be the same before comparison. Necessary for some comparators.
- Main comparator is used to decide which reference a particle most looks-like (e2help.py cmps -v2)
- Aligner - use default
- Align comparator and refine align comparator allow you to select which comparators are used for particle->reference alignment. In most cases ccc is adequate, but sometimes you may wish to match the main comparator.
- Refine align - if set to 'refine', alignments will be more accurate, and thus classification will be more accurate. Severe speed penalty.

For comparators here are some possible choices:

ccc (no options) - Simple dot product. Fast, can work well, but in some situations will cause a deterministic orientation bias (like GroEL side views which end up tilted). Works poorly for very noisy data unless usefilt particles are used for alignment.

frc zeromask=1:snrweight=1 - Fourier Ring Correlation with signal to noise ratio weighting and reference based masks. Works poorly without SNR weighting. Masking is optional, but a good idea.

phase zeromask=1:snrweight=1 - Mean phase error. same options as for frc. Do NOT use phase without snrweight=1

sqeuclidean normto=1:zeromask=1 - similar to ccc, but with additional options to better match densities. Only works well in conjunction with PS match ref, and usefilt with Wiener filtered particles.

"""
	class_documentation = """These parameters address how class-averages are made. For the comparators see the previous tab:
Averaging iterations - Use 6 or 7 when doing intial refinements to eliminate model bias. Use 2 (maybe 1) when pushing resolution
Class separation - puts each particle in the best 'n' classes, a larger number here combined with a small angular step can somewhat mimic maximum liklihood methods (smoother map at the cost of resolution)
keep - determines how many 'bad' particles are thrown away either in terms of sigma, or an absolute value (1 keeps 100%, .8 keeps 80%)
averager - either ctf.auto for ctf amplitude correction or mean for no CTF amplitude correction
set sf proj - this will filter the class-averages to match the radial power spectrum of the projections. Not good for initial rounds of refinement, but may be useful later.
"""
	make3d_documentation = """Parameters for 3D reconstruction:
- Use the default 'fourier' reconstructor
pad to - should be some number a bit larger than your box size. This should be a 'good' box size as well (see wiki)
keep - similar to keep in class-averaging, but determines how many class averages are excluded from the reconstruction
set SF - This will force the reconstruction to be filtered to match the structure factor determined during CTF correction. If used it should be combined with a gaussian lowpass filter at the targeted resolution
post-process - This is an optional filter to apply to the model as a final step, filter.lowpass.gauss with 'cutoff_freq=<1/resolution>' is good with set SF. If set SF is not used, note that the model will already \
 be somewhat filtered even without this."""
 
	def __init__(self,ptcls_list,usefilt_ptcls_list):
		E2RefineParticlesTaskBase.__init__(self,ptcls_list,usefilt_ptcls_list) 

class E2RefineMultiParticlesTask(E2RefineParticlesTaskBase, EMClassificationTools, E2Make3DTools):
	'''
	This task will harness the parameters for, and launch, e2refine.py
	'''
	 
	general_documentation = "These are the general parameters for 3D refinement in EMAN2. Please select which particles you wish to use as part of this process, specify your starting model, and fill in other parameters such as symmetry and whether or not the usefilt option should be used."
	project3d_documentation = "These  parameters are used by e2project3d. Several orientation generation techniques provide alternative methods for distributing orientations in the asymmetric unit. Orientations can be generated based on your desired angular spacing, or alternatively on the desired total number of projections. In the latter case EMAN2 will generate a number as close as possible to the specified number, but note that there is no guarantee of a perfect match. You can also vary the method by which projections are generated. If you check the \'include mirror\' option you should be sure to use aligners to that do not perform mirror alignment."
	simmx_documentation = """These  parameters are used by e2simmx, a program that compares each particle to each projection and records quality scores. \
To do this the particles must first be aligned to the projections using the aligners you specify. Once aligned the \'Main comparator\' is used to record \
the quality score. These quality values are recorded to an image matrix on handed on to the next stage in the refinement process.
- The shrink parameter causes all projections and particles to be shrunken by the given amount prior to comparison. This can provide a significant time advantage, though at the expense
of resolution. Note however that the class averaging stage, which can involve iterative alignment, does not use shrunken data.
- 2 stage simmx is still experimental. If set to 2 instead of zero, classification will be performed in two stages resulting in a 5-25x speedup, but with a potential decrease in accuracy.
- PS match ref will force the power spectra of the particle and reference to be the same before comparison. Necessary for some comparators.
- Main comparator is used to decide which reference a particle most looks-like (e2help.py cmpRefineMultiChooseSetsTasks -v2)
- Aligner - use default
- Align comparator and refine align comparator allow you to select which comparators are used for particle->reference alignment. In most cases ccc is adequate, but sometimes you may wish to match the main comparator.
- Refine align - if set to 'refine', alignments will be more accurate, and thus classification will be more accurate. Severe speed penalty.

For comparators here are some possible choices:

ccc (no options) - Simple dot product. Fast, can work well, but in some situations will cause a deterministic orientation bias (like GroEL side views which end up tilted). Works poorly for very noisy data unless usefilt particles are used for alignment.

frc zeromask=1:snrweight=1 - Fourier Ring Correlation with signal to noise ratio weighting and reference based masks. Works poorly without SNR weighting. Masking is optional, but a good idea.

phase zeromask=1:snrweight=1 - Mean phase error. same options as for frc. Do NOT use phase without snrweight=1

sqeuclidean normto=1:zeromask=1 - similar to ccc, but with additional options to better match densities. Only works well in conjunction with PS match ref, and usefilt with Wiener filtered particles.

"""
	class_documentation = """These parameters address how class-averages are made. For the comparators see the previous tab:
Averaging iterations - Use 6 or 7 when doing intial refinements to eliminate model bias. Use 2 (maybe 1) when pushing resolution
Class separation - puts each particle in the best 'n' classes, a larger number here combined with a small angular step can somewhat mimic maximum liklihood methods (smoother map at the cost of resolution)
keep - determines how many 'bad' particles are thrown away either in terms of sigma, or an absolute value (1 keeps 100%, .8 keeps 80%)
averager - either ctf.auto for ctf amplitude correction or mean for no CTF amplitude correction
set sf proj - this will filter the class-averages to match the radial power spectrum of the projections. Not good for initial rounds of refinement, but may be useful later.
"""
	make3d_documentation = """Parameters for 3D reconstruction:
- Use the default 'fourier' reconstructor
pad to - should be some number a bit larger than your box size. This should be a 'good' box size as well (see wiki)
keep - similar to keep in class-averaging, but determines how many class averages are excluded from the reconstruction
set SF - This will force the reconstruction to be filtered to match the structure factor determined during CTF correction. If used it should be combined with a gaussian lowpass filter at the targeted resolution
post-process - This is an optional filter to apply to the model as a final step, filter.lowpass.gauss with 'cutoff_freq=<1/resolution>' is good with set SF. If set SF is not used, note that the model will already \
 be somewhat filtered even without this."""
 
	def __init__(self,ptcls_list,usefilt_ptcls_list):
		E2RefineParticlesTaskBase.__init__(self,ptcls_list,usefilt_ptcls_list) 
		
	def on_form_ok(self,params):
		
		options = EmptyObject()
#		error_message = self.check_main_page(params,options)
		params["symname"]="c"
		params["symnumber"]="1"

		for checker in [self.check_main_page,self.check_project3d_page,self.check_simmx_page,self.check_classaverage_page,self.check_make3d_page]:
			error_message = checker(params,options)
			if len(error_message) > 0 :
				self.display_errors(error_message)
				return	
		self.write_db_entries(params)
		# w'oh if we make it here a lot of checking has occured. Now get the args in order to spawn_single_task	
		string_args = []
		bool_args = []
		
		additional_args = []
		
		for get_args in [self.add_general_args,self.add_project3d_args,self.add_simmx_args,self.add_classaverage_args,self.add_make3d_args]:
		  	error = get_args(options,string_args,bool_args,additional_args)
		
			if error != None: # not too fast, something still could have gone wrong
				self.display_errors([error])
				return
		
		temp_file_name = "e2refine_stdout.txt"
		
		# Steve is rethinking how we remember programs arguments
		#self.write_db_parms(options,string_args,bool_args)
		
		# Hack to enable multimodels
		if len(params["model"]) > 0:
			mlist = ""
			tsym = ""
			for i,m in enumerate(params["model"]):
				mlist += m+","
				tsym += params["symname"+str(i)]+params["symnumber"+str(i)]+","
			mlist = mlist[:-1]
			tsym = tsym[:-1]
			#string_args[string_args.index("model")] = mlist
			setattr(options, "models", mlist)
			setattr(options, "sym", tsym)

		midx = string_args.index("model")
		string_args.remove("model")
		string_args.insert(midx, "models")
		
		# Set the path, this a a hack b/c I can't decipher Davids code.
		dirnum = 1
		while os.access("refinemulti_%02d"%dirnum,os.F_OK):
			dirnum += 1
		setattr(options, 'path', ("refinemulti_%02d"%dirnum))
		os.mkdir("refinemulti_%02d"%dirnum)
		# end of hack
		
		self.spawn_single_task("e2refinemulti.py",options,string_args,bool_args,additional_args,temp_file_name)
		
		self.emit(QtCore.SIGNAL("task_idle"))
		self.form.close()
		self.form = None
		
	def get_main_params_2(self):
		
		params = []
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits="Select initial models for multireference alignment. To select click one one, and then click on another while holding down control. To add additional models to the list,use the Initial Model tab in the previous section",choices=None))

		# I could check to see if the database exists but it seems unnecessary
		# In the event that the database doesn't exist it is created and 
		# a new entry is created on disk. The only inconvenient aspect of this comes
		# if the user hits cancel - then there is a file on disk even though
		# the user never agreed to anything
		db = db_open_dict(self.form_db_name) # see eman wiki for a list of what args are kept in this db
		project_db = db_open_dict("bdb:project")
		
		self.imt = E2InitialModelsTool()
		p1,n1 = self.imt.get_initial_models_table(makebutton=0) # Do not add a Browse to Add button, it just created a headache
		QtCore.QObject.connect(p1,QtCore.SIGNAL("updateform"),self.updateform)
		
		p1.enable_multiple_selection = False
		params.append(p1)
			
		pautomask = ParamDef(name="automask3d",vartype="boolean",desc_short="Auto mask 3D",desc_long="Causes automasking of the 3D volume to occur at the end of each iteration",property=None,defaultunits=db.get("automask3d",dfl=True),choices=None)
		
		params.append(pautomask)
		
		def_rad = 30
		def_mask_dltn = 5
		if len(self.ptcls) > 0:
			nx,ny,nz = gimme_image_dimensions3D(self.ptcls[0])
			def_rad = int(nx/8.0)
			def_mask_dltn = int(nx/20.0)
		
		pamthreshold =  ParamDef(name="amthreshold",vartype="float",desc_short="Threshold",desc_long="An isosurface threshold that well defines your structure.",property=None,defaultunits=db.get("amthreshold",dfl=0.7),choices=None)
		pamradius =  ParamDef(name="amradius",vartype="int",desc_short="Radius",desc_long="The radius of a sphere at the the origin which contains seeding points for the flood file operation using the given threshold",property=None,defaultunits=db.get("amradius",dfl=def_rad),choices=None)
		pamnmax =  ParamDef(name="amnmax",vartype="int",desc_short="NMax",desc_long="Number of highest value voxels to use as a seed",property=None,defaultunits=db.get("amradius",dfl=def_rad),choices=None)
		pamnshells =  ParamDef(name="amnshells",vartype="int",desc_short="Mask dilations",desc_long="The number of dilations to apply to the mask after the flood fill operation has finished. Suggest 5% of the boxsize",property=None,defaultunits=db.get("amnshells",dfl=def_mask_dltn),choices=None)
		pamngaussshells =  ParamDef(name="amnshellsgauss",vartype="int",desc_short="Post Gaussian dilations",desc_long="The number of dilations to apply to the dilated mask, using a gaussian fall off. Suggest 5% of the boxsize",property=None,defaultunits=db.get("amnshellsgauss",dfl=def_mask_dltn),choices=None)
		
		pautomask.dependents = ["amthreshold","amradius","amnshells","amnshellsgauss"] # these are things that become disabled when the pautomask checkbox is checked etc
		params.append([pamthreshold,pamnmax,pamradius])
		params.append([pamnshells,pamngaussshells])
	
		syms = ["icos","oct","tet","d","c","h"]
		for s in xrange(n1):
			psym =  ParamDef(name="symname"+str(s),vartype="string",desc_short="Model"+str(s+1)+" Symmetry",desc_long="Symmetry to be imposed during refinement",property=None,defaultunits=db.get("symname"+str(s),dfl="c"),choices=syms)
			psymnum = ParamDef(name="symnumber"+str(s),vartype="string",desc_short="Symmetry number",desc_long="In C,D and H symmetry, this is the symmetry number",property=None,defaultunits=db.get("symnumber"+str(s),dfl="1"),choices=None)
			params.append([psym,psymnum])

		#db_close_dict(self.form_db_name)
		#db_close_dict("bdb:project")
		
		return ["Model",params]
	
	def updateform(self):
		#print "CLICKED"
		#self.on_form_cancel()
		#self.run_form()
		pass
		
	def get_project3d_page(self):
		params = []
		
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2RefineParticlesTask.project3d_documentation,choices=None))

		db = db_open_dict(self.form_db_name)
		
		projectors = self.get_projectors_list()
		orientgens = self.get_orientgens_list()
			
		pprojector =  ParamDef(name="projector",vartype="string",desc_short="Projector",desc_long="The method used to generate projections",property=None,defaultunits=db.get("projector",dfl="standard"),choices=projectors)
		
		porientgens =  ParamDef(name="orientgen",vartype="string",desc_short="Orientation generator",desc_long="The method of orientation generation",property=None,defaultunits=db.get("orientgen",dfl="eman"),choices=orientgens)
		
		pmirror = ParamDef(name="incmirror",vartype="boolean",desc_short="Include mirror",desc_long="Include the mirror portion of the asymmetric uni",property=None,defaultunits=db.get("incmirror",False),choices=[])
		
		
		orient_options = ["Angle Based", "Number Based"]
		porientoptions = ParamDef(name="orientopt",vartype="choice",desc_short="Method of generating orientation distribution",desc_long="Choose whether you want the orientations generating based on an angle or based on a total number of orientations desired",property=None,defaultunits=db.get("orientopt",dfl=orient_options[0]),choices=orient_options)
		porientoptionsentry  =  ParamDef(name="orientopt_entry",vartype="float",desc_short="value",desc_long="Specify the value corresponding to your choice",property=None,defaultunits=db.get("orientopt_entry",dfl=5),choices=[])
		from emform import EMOrientationDistDialog
		buttondialog = EMOrientationDistDialog()
		params.append([pprojector,porientgens])
		
		params.append([porientoptions,porientoptionsentry])
		params.append([pmirror,buttondialog])
		
		#db_close_dict(self.form_db_name)
		
		return ["Project 3D",params]
		
class E2EvaluateSetTask(WorkFlowTask):
	"""Launch Steven's Evaluate particles task. Just click 'ok' """
	def __init__(self):
		WorkFlowTask.__init__(self)
		
	def get_params(self):
	 	params = []
	 	params.append(ParamDef(name="blurb",vartype="text",desc_short="Select directories for FreAlign",desc_long="",property=None,defaultunits=self.__doc__,choices=None))

		return params
		
	def on_form_ok(self,params):
		
		processlist = "e2evalparticles.py"
		
		child = subprocess.Popen(processlist, shell=True)
		
		self.form.close()
		self.form = None
		
class E2RefineFromFreAlign(WorkFlowTask):
	"""Select the directory you want to use"""
	def __init__(self):
		WorkFlowTask.__init__(self)
		self.window_title = "Launch e2refinefromfrealign"
		self.report_task = None

	def get_ref_table(self):
		data_dict = EMProjectDataDict(spr_freealign_dirs_dict)
		init_model_data = data_dict.get_data_dict()
		self.project_data_at_init = init_model_data # so if the user hits cancel this can be reset
		init_model_names = init_model_data.keys()

		from emform import EM3DFileTable,EMFileTable,float_lt
		table = EM3DFileTable(init_model_names,name="dirs",desc_short="FreAlign Dirs",desc_long="")
		context_menu_data = EMRawDataReportTask.ProjectListContextMenu(spr_freealign_dirs_dict)
		table.add_context_menu_data(context_menu_data)
		table.add_button_data(EMRawDataReportTask.ProjectAddRawDataButton(table,context_menu_data))
		
		return table,len(init_model_names)
		
	def get_params(self):
	 	params = []
	 	
	 	r,rn = self.get_ref_table()
	 	params.append(ParamDef(name="blurb",vartype="text",desc_short="Select directories for FreAlign",desc_long="",property=None,defaultunits=self.__doc__,choices=None))
		params.append(r)
		return params
		
	def on_form_ok(self,params):
		
		if not params.has_key("dirs"):
			EMErrorMessageDisplay.run(["Please select dirs for processing"])
			return
		if  params.has_key("dirs") and len(params["dirs"]) == 0:
			EMErrorMessageDisplay.run(["Please select dirs for processing"])
			return
			
		self.write_db_entries(params)
		
		
		edir = params["dirs"][0]
		edir = edir[:edir.rindex("/")]
		
		prflist = "e2refinefromfrealign.py"
		prflist += " "+edir

		child = subprocess.Popen(prflist, shell=True)
		
		self.form.close()
		self.form = None
		
class E2RunFreAlign(WorkFlowTask):
	"""Select the directory you want to use for frealign"""
	def __init__(self):
		WorkFlowTask.__init__(self)
		self.window_title = "Launch e2runfrealign"
		self.report_task = None

	def get_ref_table(self):
		data_dict = EMProjectDataDict(spr_freealign_dirs_dict)
		init_model_data = data_dict.get_data_dict()
		self.project_data_at_init = init_model_data # so if the user hits cancel this can be reset
		init_model_names = init_model_data.keys()

		from emform import EM3DFileTable,EMFileTable,float_lt
		table = EM3DFileTable(init_model_names,name="dirs",desc_short="FreAlign Dirs",desc_long="")
		context_menu_data = EMRawDataReportTask.ProjectListContextMenu(spr_freealign_dirs_dict)
		table.add_context_menu_data(context_menu_data)
		table.add_button_data(EMRawDataReportTask.ProjectAddRawDataButton(table,context_menu_data))
		
		return table,len(init_model_names)
		
	def get_params(self):
	 	params = []
	 	
	 	r,rn = self.get_ref_table()
	 	params.append(ParamDef(name="blurb",vartype="text",desc_short="Select directories for FreAlign",desc_long="",property=None,defaultunits=self.__doc__,choices=None))
		params.append(r)
		#piter = ParamDef(name="iter",vartype="int",desc_short="iterations", desc_long="number of FreAlign iterations",property=None,defaultunits=5,choices=None)
		#params.append(piter)
		return params
		
	def on_form_ok(self,params):
		
		if not params.has_key("dirs"):
			EMErrorMessageDisplay.run(["Please select dirs for processing"])
			return
		if  params.has_key("dirs") and len(params["dirs"]) == 0:
			EMErrorMessageDisplay.run(["Please select dirs for processing"])
			return
			
		self.write_db_entries(params)
		
		
		edir = params["dirs"][0]
		edir = edir[:edir.rindex("/")]
		
		prflist = "e2runfrealign.py"
		
		child = subprocess.Popen(prflist, cwd=edir, shell=True)
		
		self.form.close()
		self.form = None
		

		
class E2RefineToFreeAlign(WorkFlowTask):
	"""Select the file you want to use for frealign

fbeaut: Apply extra real space symmetry averaging and masking to beautify the final map just prior to output
fcref:  Apply an FOM filter to the final reconstruction using the function SQRT(2.0*FSC/(1.0+FSC))
fstat:  Calculate additional statistics in the resolution table at the end. Using this option requires about 50% more memory
reslow: Resolution of the data in the refinement. This is the low resolution number
reshigh:Resolution of the data in the refinement. This is the high resolution number
rrec:   Resolution to which the reconstruction is calculated
thresh: Phase residual cutoff. Any particles with a higher phase residual will not be included in the reconstruction

"""
	def __init__(self):
		WorkFlowTask.__init__(self)
		self.window_title = "Launch e2refinetofrealign"
		self.form_db_name = "bdb:emform.refinetofrealign"
		self.report_task = None

	def get_ref_table(self):
		data_dict = EMProjectDataDict(spr_freealign_models_dict)
		init_model_data = data_dict.get_data_dict()
		self.project_data_at_init = init_model_data # so if the user hits cancel this can be reset
		init_model_names = init_model_data.keys()

		from emform import EM3DFileTable,EMFileTable,float_lt
		table = EM3DFileTable(init_model_names,name="model",desc_short="Starting Models",desc_long="")
		context_menu_data = EMRawDataReportTask.ProjectListContextMenu(spr_freealign_models_dict)
		table.add_context_menu_data(context_menu_data)
		table.add_button_data(EMRawDataReportTask.ProjectAddRawDataButton(table,context_menu_data))
		table.add_column_data(EMFileTable.EMColumnData("Quality",E2InitialModelsTool.get_quality_score,"This the quality score as determined by e2initialmodel.py",float_lt))

		table.add_column_data(EMFileTable.EMColumnData("Dimensions",EMRawDataReportTask.get_image_dimensions,"The dimensions of the file on disk"))
		#p.append(pdims) # don't think this is really necessary
		return table,len(init_model_names)
		
	def get_params(self):
	 	params = []
	 	db = db_open_dict(self.form_db_name)
	 	
	 	r,rn = self.get_ref_table()
	 	params.append(ParamDef(name="blurb",vartype="text",desc_short="Interactive use of tomohunter",desc_long="",property=None,defaultunits=self.__doc__,choices=None))
		params.append(r)
		pfbeaut = ParamDef(name="fbeaut",vartype="boolean",desc_short="fbeaut",desc_long="FreeAlign fbeaut",property=None,defaultunits=db.get("fbeaut",dfl=False),choices=None)
		pfcref = ParamDef(name="fcref",vartype="boolean",desc_short="fcref",desc_long="FreeAlign cref",property=None,defaultunits=db.get("fcref",dfl=False),choices=None)
		pfstat = ParamDef(name="fstat",vartype="boolean",desc_short="fstat",desc_long="FreeAlign fstat",property=None,defaultunits=db.get("fstat",dfl=False),choices=None)
		params.append([pfbeaut, pfcref, pfstat])
		preslow = ParamDef(name="reslow",vartype="float",desc_short="reslow", desc_long="low resolution",property=None,defaultunits=db.get("reslow",dfl=200.0),choices=None)
		preshigh = ParamDef(name="reshigh",vartype="float",desc_short="reshigh", desc_long="high resolution",property=None,defaultunits=db.get("reshigh",dfl=25.0),choices=None)
		params.append([preslow, preshigh])
		prrec = ParamDef(name="rrec",vartype="float",desc_short="rrec", desc_long="rrec for FreeAlign",property=None,defaultunits=db.get("rrec",dfl=10.0),choices=None)
		pthresh = ParamDef(name="thresh",vartype="float",desc_short="thresh", desc_long="FreeAlign threshold",property=None,defaultunits=db.get("thresh",dfl=90.0),choices=None)
		params.append([prrec, pthresh])
		
		return params
		
	def on_form_ok(self,params):
		
		if not params.has_key("model"):
			EMErrorMessageDisplay.run(["Please select files for processing"])
			return
		if  params.has_key("model") and len(params["model"]) == 0:
			EMErrorMessageDisplay.run(["Please select files for processing"])
			return
			
		self.write_db_entries(params)
		
		redata = params["model"][0]
		try:
			directory = re.compile(r':.+#').findall(redata)[0][1:-1]
		except:
			raise RuntimeError("The directory %s does not conform to the format: bdb:refinedir#threed_??", redata)
		try:
			iteration = re.compile(r'#threed.*_[0-9]+').findall(redata)[0][-2:]
		except:
			raise RuntimeError("The model %s does not conform to the format: bdb:refinedir#threed_??+", redata)
		iteration = str(int(iteration)) # This is a hack, to remove preceding zeros
		
		e2falist = "e2refinetofrealign.py"
		e2falist += " "+directory
		e2falist += " "+iteration
		if(params["fbeaut"]):
			e2falist += " --fbeaut"
		if(params["fcref"]):
			e2falist += " --fcref"
		if(params["fstat"]):
			e2falist += " --fstat"
		e2falist += (" --rrec="+str(params["rrec"]))
		e2falist += " --reslow="+str(params["reslow"])
		e2falist += " --reshigh="+str(params["reshigh"])
		e2falist += " --thresh="+str(params["thresh"])
		
		print e2falist
		child = subprocess.Popen(e2falist, shell=True)
		
		self.form.close()
		self.form = None
		
		
class EMChooseTask(ParticleWorkFlowTask):
	'''Fill me in'''
	documentation_string = "This selects the type of data to use for the refinement. Normally phase_flipped or phase_flipped_hp will be selected in the top \
section. Usefilt permits you to use one set of data for alignment and classification, and the set above when actually building a model. This can provide significant \
benefits in some cases. Normally the wiener_filtered set or None will be selected here. If you are using the 'prefilt' option in the refinement parameters, you \
should use Wiener filtered particles here. If not, you may experiment both ways. After completing this form, you will be presented with a set of forms with \
detailed refinement options." 
	def __init__(self,filt_name,ptcl_opts,choice_subscript=""):
		ParticleWorkFlowTask.__init__(self)
		self.window_title = "Choose Data For e2refine"
		self.preferred_size = (480,300)
		self.form_db_name = "bdb:emform.e2refine"
		self.filt_name = filt_name
		self.ptcl_opts = ptcl_opts
		self.choice_subscript = choice_subscript # for data base uniqueness
		self.ptcl_choice = "ptcl_choice"+choice_subscript
		self.usefilt_ptcl_choice = "usefilt_ptcl_choice"+choice_subscript
		self.single_selection = False
		
		self.name_map = None
		self.particles_map = None
		self.particles_name_map = None
		
	def get_params(self):
		self.particles_map, self.particles_name_map,choices,self.name_map = self.ptcl_opts.get_particle_options()
		choices.append("Specify")
		
		db = db_open_dict(self.form_db_name)
		ppart = ParamDef(name=self.ptcl_choice,vartype="choice",desc_long="Choose the particle data set you wish to use to generate a starting data for e2refine",desc_short="Particle data",property=None,defaultunits=db.get(self.ptcl_choice,dfl=""),choices=choices)
		import copy
		usefilt_choices = copy.deepcopy(choices)
		usefilt_choices.append("None")
		usefilt_choices.remove("Specify")
		pusefilt = ParamDef(name=self.usefilt_ptcl_choice,vartype="choice",desc_long="Choose the particle data set you wish to use to generate a starting data for e2refine",desc_short="Usefilt data",property=None,defaultunits=db.get(self.usefilt_ptcl_choice,dfl="None"),choices=usefilt_choices)
		params = []		
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=self.__doc__,choices=None))
		params.append(ppart)
		params.append(pusefilt)
		return params
	
	def on_form_ok(self,params):
		self.write_db_entries(params)
		if not params.has_key(self.ptcl_choice) or params[self.ptcl_choice] == None:
			error("Please choose particle data")
			return
		if not params.has_key(self.usefilt_ptcl_choice) or params[self.usefilt_ptcl_choice] == None:
			error("Please choose usefilt data")
			return
		choice = params[self.ptcl_choice]
		usefilt_choice = params[self.usefilt_ptcl_choice]
		
		if choice == usefilt_choice:
			EMErrorMessageDisplay.run(["The Particle and Usefilt data must be different.\nIf you don't want to include Usefilt data choose None"])
			return
		
		if choice == "Specify" or usefilt_choice =="Specify":
			EMErrorMessageDisplay.run(["Specify not currently supported"])
			return
		else:
			err = []
			bad = []
			filt_name = self.particles_name_map[choice]
			ptcls = self.particles_map[filt_name]
			db = db_open_dict("bdb:project")
			
			data_dict = EMProjectDataDict(self.filt_name)
			db_map = data_dict.get_data_dict()
			ptcl_keys = db_map.keys()

			if usefilt_choice != "None":
				usefilt_name = self.particles_name_map[usefilt_choice]
				
				intersection_ptcls = []
				intersection_usefilt_ptcls = []
				
				for name in ptcls:
					d = self.name_map[name]
					filt_map = db_map[d]
					if filt_map.has_key(usefilt_name):
						name2 = filt_map[usefilt_name]
						if EMUtil.get_image_count(name) != EMUtil.get_image_count(name2):
							err.append("The number of images in the usefilt file %s is not consistent with the particles file %s" %(name,name2))
							bad.append(name)
							continue
						intersection_ptcls.append(name)
						intersection_usefilt_ptcls.append(name2)
				
			else:
				intersection_ptcls = ptcls
				intersection_usefilt_ptcls = None
				
			if len(err) > 0:
				e = "\n\n All of these sets will now be automatically removed from the list of sets in the project."
				err.append(e)
				error(err)
				
				for b in bad: db_map.pop(b)
				db[self.filt_name] = db_map 
			
			task = self.task_to_use(intersection_ptcls,intersection_usefilt_ptcls)
			task.single_selection = self.single_selection
			self.emit(QtCore.SIGNAL("replace_task"),task,"e2refine2d arguments")
			self.form.close()
			self.form = None
		
		self.write_db_entries(params)

	def replace_task(self): raise NotImplementedException("Inheriting class needs to implement this")
	
	def task_to_use(self, intersection_ptcls, intersection_usefilt_ptcls): raise NotImplementedException("Inheriting class needs to implement this")


class E2RefineChooseParticlesTask(EMChooseTask):
	'''This form is for choosing e2refine input and usefilt data. After you hit ok a second form will appear asking for more parameters.'''
	documentation_string = "This form is for choosing e2refine input and usefilt data. After you hit ok a second form will appear asking for more parameters." 
	def __init__(self):
		EMChooseTask.__init__(self,spr_ptcls_dict,EMParticleOptions(),"")
		self.window_title = "Choose Data For e2refine"
		self.preferred_size = (480,300)
		self.form_db_name = "bdb:emform.e2refine"
		
	# Lets you swap out one refine form for another (Template method pattern, is used here)
	def task_to_use(self, intersection_ptcls, intersection_usefilt_ptcls):
		task = E2RefineParticlesTask(intersection_ptcls,intersection_usefilt_ptcls)
		return task

class E2RefineChooseSetsTask(EMChooseTask):
	'''This form is for choosing e2refine input and usefilt data. After you hit ok a second form will appear asking for more parameters.'''
	documentation_string = "This form is for choosing e2refine input and usefilt data. After you hit ok a second form will appear asking for more parameters." 
	def __init__(self):
		EMChooseTask.__init__(self,spr_sets_dict,EMSetsOptions(),"_stack")
		self.window_title = "Choose Data For e2refine"
		self.preferred_size = (480,300)
		self.form_db_name = "bdb:emform.e2refine"
		self.single_selection = True
		
	# Lets you swap out one refine form for another (Template method pattern, is used here)
	def task_to_use(self, intersection_ptcls, intersection_usefilt_ptcls):
		task = E2RefineParticlesTask(intersection_ptcls,intersection_usefilt_ptcls)
		return task
		
class E2RefineMultiChooseSetsTask(EMChooseTask):
	'''This form is for choosing e2refine input and usefilt data. After you hit ok a second form will appear asking for more parameters.'''
	documentation_string = "This form is for choosing e2refinemulti input and usefilt data. After you hit ok a second form will appear asking for more parameters." 
	def __init__(self):
		EMChooseTask.__init__(self,spr_sets_dict,EMSetsOptions(),"_stack")
		self.window_title = "Choose Data For e2refinemulti"
		self.preferred_size = (480,300)
		self.form_db_name = "bdb:emform.e2refinemulti"
		self.single_selection = True
		
	# Lets you swap out one refine form for another (Template method pattern, is used here)
	def task_to_use(self, intersection_ptcls, intersection_usefilt_ptcls):
		task = E2RefineMultiParticlesTask(intersection_ptcls,intersection_usefilt_ptcls)
		return task	
		
class ResolutionReportTask(ParticleWorkFlowTask):
	documentation_string = "This form displays information related to the estimated resolution of refinement results.\nIf you double click on any of the entries \
 in the table you will see the convergence plot (thin lines, NOT useful in measuring resolution) and any associated resolution curve (thick line for resolution \
 determination, must run eotest first)."
	warning_string = "\n\n\nNOTE: There are no results available."
	def __init__(self):
		ParticleWorkFlowTask.__init__(self)
		self.window_title = "Project Resolution Report"
	
	def get_params(self):
		params = []
		
		p,n = self.get_resolution_table()
		
		if n == 0:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=ResolutionReportTask.documentation_string+ResolutionReportTask.warning_string,choices=None))
		else:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=ResolutionReportTask.documentation_string,choices=None))
			params.append(p)
		return params
	
	def find_first_point_5_crossing(self,xaxis,yaxis):
		'''
		Find the first 0.5 crossing in the FSC - interpolate and try to return an accurate estimate
		Disregards the first five entries
		Assumes the Nyquist frequency is correct
		'''
		idx = 0
		if len(yaxis) > 5:
			idx = 6
		
		soln = -1
		while (idx < len(yaxis)-1 ):
			if yaxis[idx] >= 0.5 and yaxis[idx+1] <= 0.5:
				v1 = yaxis[idx]
				v2 = yaxis[idx+1]
				if v1 != v2:
					d = v1-v2
					offset = v1-0.5
					interp = offset/d
					soln = idx+interp
				else: soln = idx
				break
			else:
				idx += 1
		
		try:
			if soln == -1:
				return "???"
			elif int(soln) == soln:
				return "%.1f" %(1.0/xaxis(soln))
			else:
				# interpolated frequency
				return "%.1f" %(1.0/(soln/len(yaxis)*xaxis[-1]))
		except:
			return "invalid"
		
				
	def get_resolution_table(self):
		'''
		Looks for bdb:r2d_??#classes_?? and the bdb:r2d_??#classes_init file, finds the most recent one, then fills in the number of particles in
		in the class average file and also its dimensions.
		'''
		dirs = get_prefixed_directories("refine_") + get_prefixed_directories("frealign_") + get_prefixed_directories("refinemulti_")
		
		available_dirs = []
		total_iterations = []
		eotest_res = []
		e2resolution_res = []
		for dir in dirs:
			db_name = "bdb:"+dir+"#convergence.results"
			if db_check_dict(db_name):
				db = db_open_dict(db_name,ro=True)
				keys = db.keys()
				if len(keys) > 0:
					available_dirs.append(dir)
					
					res = get_e2resolution_results_list(keys)
					eo = get_e2eotest_results_list(keys)
					conv = get_convergence_results_list(keys)
					total_iterations.append(len(conv))
					
					if len(res) > 0:
						# get the latest one, this will be the last as guaranteed by sorted results
						last_res = res[-1]
						[xaxis,yaxis] = db[last_res]
						resolution = self.find_first_point_5_crossing(xaxis,yaxis)
						e2resolution_res.append(resolution)
					else:
						e2resolution_res.append("")
						
					if len(eo) > 0:
						last_res = eo[-1]
						[xaxis,yaxis] = db[last_res]
						resolution = self.find_first_point_5_crossing(xaxis,yaxis)
						eotest_res.append(resolution)
					else:
						eotest_res.append("")
				
				#db_close_dict(db_name)

		if len(available_dirs) > 0:
			p = EMParamTable(name="filenames",desc_short="Resolution Evaluation",desc_long="")
			pnames = ParamDef(name="Dirs",vartype="stringlist",desc_short="Refinement Directory", desc_long="EMAN2 refinement directory", property=None,defaultunits=None,choices=available_dirs)
			piter = ParamDef(name="Iterations",vartype="intlist",desc_short="Total Iterations",desc_long="The number of 3D refinement iterations that have occured in this directory",property=None,defaultunits=None,choices=total_iterations)
			peo = ParamDef(name="eo",vartype="stringlist",desc_short="e2eotest",desc_long="0.5 e2eotest resolution estimate",property=None,defaultunits=None,choices=eotest_res)
			pres = ParamDef(name="res",vartype="stringlist",desc_short="e2resolution",desc_long="0.5 e2resolution resolution estimate",property=None,defaultunits=None,choices=e2resolution_res)
			
			p.append(pnames)
			p.append(piter)
			p.append(peo)
			p.append(pres)
			
			setattr(p,"convert_text", resolution_display_convert)
			setattr(p,"icon_type","2d_plot")
			
			return p,len(available_dirs)
		else:
			return None,0
		
def resolution_display_convert(dir):
	'''
	This "display_convert" function breaks the mold a little in that it returns the name of a database that has plots in it. It's handled in e2workflow.py
	'''
	return "bdb:"+dir+"#convergence.results"


def get_e2resolution_results_list(keys):
		'''
		Extract the names from the keys that match the e2resolution.py output naming convention
		(keys is a list of keys in the convergence.results dictionary, in a refinement directory)
		'''
		solns = []
		for k in keys:
			if len(k) > 6 and k[-7:] == "res_fsc":
				solns.append(k)
		solns.sort()
		return solns
	
def get_e2eotest_results_list(keys):
	'''
	Extract the names from the keys that match the e2eotest.py output naming convention
	(keys is a list of keys in the convergence.results dictionary, in a refinement directory)
	'''
	solns = []
	for k in keys:
		if len(k) > 7 and k[0:8] == "even_odd":
			solns.append(k)
	solns.sort()
	return solns

def get_convergence_results_list(keys):
	'''
	Extract the names from the keys that match the e2refine.py convergence plot output naming convention
	(keys is a list of keys in the convergence.results dictionary, in a refinement directory)
	'''
	solns = []
	if "init_00_fsc" in keys:
		solns.append("init_00_fsc")
		
	i = 0
	while True:
		s1 = str(i)
		s2 = str(i+1)
		if len(s1) == 1: s1 = "0"+s1
		if len(s2) == 1: s2 = "0"+s2
		k = s1+"_"+s2+"_fsc"
		if k in keys:
			solns.append(k)
		else:
			break

		i += 1
	return solns

class E2EotestTask(EMClassificationTools,E2Make3DTools):
	'''
	Run e2eotest from the workflow setting
	Inherits from E2RefineParticlesTask because it uses some forms that are/almost identical
	'''
	 
	general_documentation = "These are parameters required to run an even-odd test in EMAN2"
	documentation_string = "This form allows you to specify parameters for the even-odd resolution test. Normally these parameters should match \
those used during refinement."
	warning_string = "\n\n\nThere are no refinement results available to use as the basis of running e2eotest"
	def __init__(self):
	 	EMClassificationTools.__init__(self)
	 	E2Make3DTools.__init__(self)
	 	self.window_title = "Run e2eotest"
	 	self.form_db_name = "bdb:emform.e2eotest"
	 	self.dir_and_iter = {} # will eventually be useful information about directories that will work for e2eotest
	 		
	def run_form(self):
		self.form = EMTableFormWidget(self.get_params())
		self.form.resize(*self.preferred_size)
		self.form.setWindowTitle(self.window_title)
		get_application().show_specific(self.form)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_ok"),self.on_form_ok)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_cancel"),self.on_form_cancel)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_close"),self.on_form_close)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("display_file"),self.on_display_file)
		
	
	def get_params(self):
	 	params = []
		
		# do this so that we have 
		self.__set_available_iteration_data()
		
		
		if len(self.dir_and_iter) == 0:
			params.append(self.get_cant_proceed_page())
			return params
		
		params.append(self.get_main_params())
		params.append(self.get_classaverage_page(include_sep=False))
		params.append(self.get_make3d_page())
		
		return params
	
	def get_cant_proceed_page(self):
		params = []		

		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2EotestTask.documentation_string+E2EotestTask.warning_string,choices=None))

		return ["Can't proceed",params]
	def __set_available_iteration_data(self):
		'''
		This function is called in get_params to accrue the directory data
		'''
		dirs = get_prefixed_directories("refine_") + get_prefixed_directories("frealign_") + get_prefixed_directories("refinemulti_")
		dirs.sort()
		
		nec_files = [ "classes_", "classify_","projections_"]
		
		self.dir_and_iter = {}
		for dir in dirs:
			fail = False
			available_iters = []
			
			register_db_name = "bdb:"+dir+"#register"
				
			# needs to exist
			if not db_check_dict(register_db_name):
				continue
			# cmd dictionary needs to be stored
			db = db_open_dict(register_db_name,ro=True)
			if not db.has_key("cmd_dict") or db["cmd_dict"]==None:
				continue
	
			cmd = db["cmd_dict"]
			# need to be able to get the input data
			if not cmd.has_key("input"):
				continue
			
			for i in range(0,10):
				for j in range(0,10):
					end = str(i) + str(j)
					for file in nec_files:
						db_first_part = "bdb:"+dir+"#" + file
						db_name = db_first_part + end
						if not db_check_dict(db_name):
							fail = True
							break
					if not fail:
						available_iters.append(end)
					else:
						break
				if fail: break
			
			# might be some empyt ones so just forget those
			if len(available_iters) > 0: 
				available_iters.reverse()
				self.dir_and_iter[dir] = available_iters
	
	def get_main_params(self):
		'''
		General/broad refine params
		'''
		# have to get the directories 
		params = []
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2EotestTask.general_documentation,choices=None))
		
		# I could check to see if the database exists but it seems unnecessary
		# In the event that the database doesn't exist it is created and 
		# a new entry is created on disk. The only inconvenient aspect of this comes
		# if the user hits cancel - then there is a file on disk even though
		# the user never agreed to anything
		db = db_open_dict(self.form_db_name) # see eman wiki for a list of what args are kept in this db
		
		params.append(ParamDef(name="path and iteration", vartype="dict",desc_short="Directory and iteration",desc_long="Select the directory containing the refinement and the iteration you wish to use as the input", property=None, defaultunits="",choices=self.dir_and_iter  ))
	
		
		plowmem = ParamDef(name="lowmem",vartype="boolean",desc_short="Low mem",desc_long="Causes various programs to restrict memory usage but results in increased CPU time.",property=None,defaultunits=db.get("lowmem",dfl=False),choices=None)
		pusefilt = ParamDef(name="usefilt",vartype="boolean",desc_short="Usefilt",desc_long="Will use the 'usefilt' data for class alignment if it exists in the refinement directory",property=None,defaultunits=db.get("usefilt",dfl=False),choices=None)
		
		syms = ["icos","oct","tet","d","c","h"]
		
		psym =  ParamDef(name="symname",vartype="string",desc_short="Symmetry",desc_long="Symmetry to be imposed during refinement",property=None,defaultunits=db.get("symname",dfl="c"),choices=syms)
		psymnum = ParamDef(name="symnumber",vartype="string",desc_short="Symmetry number",desc_long="In C,D and H symmetry, this is the symmetry number",property=None,defaultunits=db.get("symnumber",dfl="1"),choices=None)
		

		params.append([plowmem,pusefilt])
		params.append([psym,psymnum])
		
		
		#db_close_dict(self.form_db_name)
		
		return ["General",params]
	
	
	def check_main_page(self,params,options):
		error_message = []
		
		options.path = params["path"]
		options.iteration = params["iteration"]
		options.lowmem = params["lowmem"]
		
		if params["usefilt"] == True:
			
			register_db_name = "bdb:"+options.path+"#register"
			db = db_open_dict(register_db_name)
			cmd = db["cmd_dict"]
			# need to be able to get the input data
			if not cmd.has_key("usefilt") or str(cmd["usefilt"]) == "None":
				error_message.append("You have checked usefilt but no usefilt data was used for the refinement in the given directory")
			else:
				options.usefilt = file = cmd["usefilt"]

		
		error_message.extend(self.check_sym(params,options))
		
		return error_message
		
	def add_general_args(self,options,string_args,bool_args,additional_args):
		
		
		if hasattr(options,"usefilt"): string_args.append("usefilt")
		string_args.extend(["path","iteration","sym"])
		bool_args.append("lowmem")
		
		return None # returning None is good
		
	
	def on_form_ok(self,params):

		if len(self.dir_and_iter) == 0: return # The user has the can't proceed page

		options = EmptyObject()
		
		for checker in [self.check_main_page,self.check_classaverage_page,self.check_make3d_page]:
			error_message = checker(params,options)
			if len(error_message) > 0 :
				EMErrorMessageDisplay.run(error_message)
				return
			

		self.write_db_entries(params)
#		# w'oh if we make it here a lot of checking has occured. Now get the args in order to spawn_single_task
		string_args = []
		bool_args = []
		
		additional_args = ["--force"]
#		
		for get_args in [self.add_general_args,self.add_classaverage_args,self.add_make3d_args]:
			if get_args == self.add_classaverage_args:
				error = get_args(options,string_args,bool_args,additional_args,include_sep=False)
			else:
				error = get_args(options,string_args,bool_args,additional_args)
		
			if error != None: # not too fast, something still could have gone wrong
				self.display_errors([error])
				return
			
		temp_file_name = "e2eotest_stdout.txt"
#		
#		self.write_db_parms(options,string_args,bool_args)
#		
	   	options.filenames = [] # spawn single task expects a filenames attribute
		self.spawn_single_task("e2eotest.py",options,string_args,bool_args,additional_args,temp_file_name)
		self.emit(QtCore.SIGNAL("task_idle"))
		self.form.close()
		self.form = None
#	
class E2ResolutionTask(WorkFlowTask):
	'''

	'''
	 
	general_documentation = "These are parameters required to run an e2resolution.py."
	warning_string = "\n\n\nThere are no refinement results available to use as the basis of running e2resolution"
	def __init__(self):
	 	WorkFlowTask.__init__(self)
	 	self.window_title = "Run e2resolution"
	 	self.dir_and_iter = {} # will eventually be useful information about directories that will work for e2eotest
	 	
		
	def get_params(self):
		'''
		General/broad refine params
		'''
		# have to get the directories 
		dirs = get_prefixed_directories("refine_") + get_prefixed_directories("frealign_") + get_prefixed_directories("refinemulti_")
		dirs.sort()
		
		nec_files = ["threed_filt_","threed_mask_"]
		
		self.dir_and_iter = {}
		for dir in dirs:
			fail = False
			available_iters = []
			for i in range(0,10):
				for j in range(0,10):
					end = str(i) + str(j)
					for file in nec_files:
						db_first_part = "bdb:"+dir+"#" + file
						db_name = db_first_part + end
						if not db_check_dict(db_name):
							fail = True
							break
					if not fail:
						available_iters.append(end)
					else:
						break
				if fail: break
			
			# might be some empyt ones so just forget those
			if len(available_iters) > 0: 
				available_iters.reverse()
				self.dir_and_iter[dir] = available_iters
	
		params = []
		if len(self.dir_and_iter) == 0:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2ResolutionTask.general_documentation+ E2ResolutionTask.warning_string,choices=None))
			return params
		
		params = []
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2ResolutionTask.general_documentation,choices=None))
		
	
		params.append(ParamDef(name="path and iteration", vartype="dict",desc_short="Directory and iteration",desc_long="Select the directory containing the refinement and the iteration you wish to use as the input", property=None, defaultunits="",choices=self.dir_and_iter  ))
		
		project_db = db_open_dict("bdb:project")
		papix = ParamDef(name="global.apix",vartype="float",desc_short="Angtsrom per pixel",desc_long="The physical distance represented by the pixel spacing",property=None,defaultunits=project_db.get("global.apix",dfl=1.1),choices=None)
		params.append(papix)
		
		#db_close_dict("bdb:project")
		
		return params
	
	def on_form_ok(self,params):

		if len(self.dir_and_iter) == 0: return # the user has been issued a warning about data lacking, ok does nothing

		if  params.has_key("global.apix") and params["global.apix"] <= 0:
			self.show_error_message(["Apix must be greater than  zero"])
			return
		else:
			self.write_db_entries(params)

		options = EmptyObject()			
		options.path = params["path"]
		options.apix = params["global.apix"]
		
		string_args = ["path","apix"]
		bool_args = []
		
		image = "bdb:"+params["path"] + "#threed_filt_"+params["iteration"]
		mask = "bdb:"+params["path"] + "#threed_mask_"+params["iteration"]
		output = "tmp"
		additional_args = [image,mask,output]

		temp_file_name = "e2resolution_stdout.txt"
	   	options.filenames = [] # spawn single task expects a filenames attribute
		self.spawn_single_task("e2resolution.py",options,string_args,bool_args,additional_args,temp_file_name)
		self.emit(QtCore.SIGNAL("task_idle"))
		self.form.close()
		self.form = None
#	

	
if __name__ == '__main__':
	
	from emapplication import EMApp
	em_app = EMApp()
	sprinit = SPRInitTask()
	window = sprinit.run_form() 
	#em_app.show()
	em_app.execute()	
