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

from emform import EMFormModule,ParamTable
from emdatastorage import ParamDef
from PyQt4 import QtGui,QtCore
from EMAN2db import db_check_dict, db_open_dict,db_remove_dict,db_list_dicts,db_close_dict
from EMAN2 import EMData,get_file_tag,EMAN2Ctf,num_cpus,memory_stats,check_files_are_2d_images,check_files_are_em_images,numbered_path,dump_aligners_list,dump_cmps_list
import os
import copy
from emapplication import EMProgressDialogModule
from e2boxer import EMBoxerModule
from e2ctf import pspec_and_ctf_fit,GUIctfModule,write_e2ctf_output,get_gui_arg_img_sets
import subprocess
from pyemtbx.boxertools import set_idd_image_entry, TrimBox
import weakref

class EmptyObject:
	'''
	This just because I need an object I can assign attributes to, and object() doesn't seem to work
	'''
	def __init__(self):
		pass
	
	

class WorkFlowTask(QtCore.QObject):
	def __init__(self,application):
		QtCore.QObject.__init__(self)
		self.application = weakref.ref(application)
		self.window_title = "Set me please"
		self.preferred_size = (480,640)
	
	def run_form(self):
		self.form = EMFormModule(self.get_params(),self.application())
		self.form.qt_widget.resize(*self.preferred_size)
		self.form.setWindowTitle(self.window_title)
		self.application().show_specific(self.form)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_ok"),self.on_form_ok)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_cancel"),self.on_form_cancel)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_close"),self.on_form_close)
		
	def on_form_ok(self,params):
		for k,v in params.items():
			self.write_db_entry(k,v)
			
		self.application().close_specific(self.form)
		self.form = None
	
		self.emit(QtCore.SIGNAL("task_idle"))
		
	def on_form_cancel(self):
		self.application().close_specific(self.form)
		self.form = None

		self.emit(QtCore.SIGNAL("task_idle"))
	
	def on_form_close(self):
		self.emit(QtCore.SIGNAL("task_idle"))

	def closeEvent(self,event):
		self.application().close_specific(self.form)
		#self.emit(QtCore.SIGNAL("task_idle")
		
	def write_db_entry(self,key,value):
		'''
		Call this function if you need to
		'''
		db = db_open_dict("bdb:project")
		if key == "global.apix":
			db["global.apix"] = value
		elif key == "global.microscope_voltage":
			db["global.microscope_voltage"] = value
		elif key == "global.microscope_cs":
			db["global.microscope_cs"] = value
		elif key == "global.memory_available":
			#mem = memory_stats()
			#if value > mem[0]:
				#print "error, memory usage is beyond the total amount available"
			#else:
			# we're avoiding validation because users have peculiar reasons some times
			db["global.memory_available"] = value
		elif key == "global.num_cpus":
			#n = num_cpus()
			#if value > n:
				#print "error, num_cpus more than the available cpus"
			#else:
			#  we're avoiding validation because users have peculiar reasons some times
			db["global.num_cpus"] = value
		else:
			pass
		
		db_close_dict("bdb:project")
	def get_wd(self):
		'''
		Get the working directory, originally introduced to provide a centralized mechanism for accessing the working directory,
		specificially for the purpose of spawning processes. Could be used more generally, however.
		'''
		return os.getcwd()

	def run_task(self,program,options,string_args,bool_args,additional_args=[],temp_file_name="e2workflow_tmp.txt"):
		'''
		splits the task over the available processors
		example-
		program="e2ctf.py"
		options is an object with the filenames, all string args and all bool_args as attributes 
		string_args=["
		bool_args=["
		additional_args=["--auto_db,--auto_fit"]
		temp_file_name = "etctf_auto_tmp.txt"
		'''
		project_db = db_open_dict("bdb:project")
		ncpu = project_db.get("global.num_cpus",dfl=num_cpus())
		cf = float(len(options.filenames))/float(ncpu) # common factor
		for n in range(ncpu):
			b = int(n*cf)
			t = int(n+1*cf)
			if n == (ncpu-1):
				t = len(options.filenames) # just make sure of it, round off error could 
			
			if b == t:
				continue # it's okay this happens when there are more cpus than there are filenames	
			filenames = options.filenames[b:t]
								
			args = [program]
	
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
			#print "command is ",program
			#for i in args: print i
			
			#print args
			file = open(temp_file_name,"w+")
			process = subprocess.Popen(args,stdout=file,stderr=subprocess.STDOUT)
			print "started process",process.pid
			self.emit(QtCore.SIGNAL("process_started"),process.pid)
			
		db_close_dict("bdb:project")


class SPRInitTask(WorkFlowTask):
	'''
	A class that manages the initialization component of a Single Particle
	Reconstruction workflow
	'''
	
	# stolen from wikipedia
	documentation_string = "In physics, in the area of microscopy, single particle reconstruction is a technique in which large numbers of images (10,000 - 1,000,000) of ostensibly identical individual molecules or macromolecular assemblies are combined to produce a 3 dimensional reconstruction. This is a complementary technique to crystallography of biological molecules. As molecules/assembies become larger, it becomes more difficult to prepare high resolution crystals. For single particle reconstruction, the opposite is true. Larger objects actually improve the resolution of the final structure. In single particle reconstruction, the molecules/assemblies in solution are prepared in a thin layer of vitreous (glassy) ice, then imaged on an electron cryomicroscope (see Transmission electron microscopy). Images of individual molecules/assemblies are then selected from the micrograph and then a complex series of algorithms is applied to produce a full volumetric reconstruction of the molecule/assembly. In the 1990s this technique was limited to roughly 2 nm resolution, providing only gross features of the objects being studied. However, recent improvements in both microscope technology as well as available computational capabilities now make 0.5 nm resolution possible."
	
	def __init__(self,application):
		WorkFlowTask.__init__(self,application)
		self.window_title = "Project information"
	def get_params(self):
		params = []
		project_db = db_open_dict("bdb:project")
		#params.append(ParamDef(name="global.micrograph_ccd_filenames",vartype="url",desc_short="File names",desc_long="The raw data from which particles will be extracted and ultimately refined to produce a reconstruction",property=None,defaultunits=db_entry("global.micrograph_ccd_filenames","bdb:project",[]),choices=[]))
		params.append(ParamDef(name="blurb",vartype="text",desc_short="SPR",desc_long="Information regarding this task",property=None,defaultunits=SPRInitTask.documentation_string,choices=None))
		
		papix = ParamDef(name="global.apix",vartype="float",desc_short="A/pix for project",desc_long="The physical distance represented by the pixel spacing",property=None,defaultunits=project_db.get("global.apix",dfl=1.1),choices=None)
		pvolt = ParamDef(name="global.microscope_voltage",vartype="float",desc_short="Microscope voltage",desc_long="The operating voltage of the microscope in kilo volts",property=None,defaultunits=project_db.get("global.microscope_voltage",dfl=300),choices=None)
		pcs = ParamDef(name="global.microscope_cs",vartype="float",desc_short="Microscope Cs",desc_long="Microscope spherical aberration constant",property=None,defaultunits=project_db.get("global.microscope_cs",dfl=2.0),choices=None)
		pncp = ParamDef(name="global.num_cpus",vartype="int",desc_short="Number of CPUs",desc_long="Number of CPUS available for the project to use",property=None,defaultunits=project_db.get("global.num_cpus",dfl=num_cpus()),choices=None)
		mem = memory_stats()
		pmem = ParamDef(name="global.memory_available",vartype="float",desc_short="Memory usage ("+str(mem[0])+ " Gb total)",desc_long="The total amount of system memory you want to make available to the project in gigabytes",property=None,defaultunits=project_db.get("global.memory_available",dfl=mem[1]),choices=None)
		params.append(papix)
		params.append(pvolt)
		params.append(pcs)
		params.append(pncp)
		params.append(pmem)
		db_close_dict("bdb:project")
		return params

	def write_db_entry(self,key,value):
		WorkFlowTask.write_db_entry(self,key,value)		


class MicrographCCDImportTask(WorkFlowTask):	
	documentation_string = "Use this tool for importing flat files into the raw_data directory in the project database. Files that you import in this way will be automatically added the list of files in the project."
	
	def __init__(self,application):
		WorkFlowTask.__init__(self,application)
		self.window_title = "Import micrographs"
		self.thumb_shrink = -1
	def get_params(self):
		params = []
		project_db = db_open_dict("bdb:project")
		params.append(ParamDef(name="blurb",vartype="text",desc_short="Importing image data",desc_long="",property=None,defaultunits=MicrographCCDImportTask.documentation_string,choices=None))
		params.append(ParamDef(name="import_micrograph_ccd_files",vartype="url",desc_short="File names",desc_long="The raw data from which particles will be extracted and ultimately refined to produce a reconstruction",property=None,defaultunits=[],choices=[]))
		pinvert = ParamDef(name="invert",vartype="boolean",desc_short="Invert",desc_long="Tick this if you want eman2 to invert your images while importing",property=None,defaultunits=False,choices=None)
		pxray = ParamDef(name="xraypixel",vartype="boolean",desc_short="X-ray pixel",desc_long="Tick this if you want eman2 to automatically filter out X-ray pixels while importing",property=None,defaultunits=False,choices=None)
		pthumbnail = ParamDef(name="thumbs",vartype="boolean",desc_short="Thumbnails",desc_long="Tick this if you want eman2 to automatically generate thumbnails for your images. This will save time at later stages in the project",property=None,defaultunits=True,choices=None)
		
		params.append([pinvert,pxray,pthumbnail])
		
		db_close_dict("bdb:project")
		return params
	
	def on_form_ok(self,params):
		for k,v in params.items():
			if k == "import_micrograph_ccd_files":
				self.do_import(params)
			else:
				self.write_db_entry(k,v)
		
		
		
		
		self.application().close_specific(self.form)
		self.form = None
	
		self.emit(QtCore.SIGNAL("task_idle"))

	def do_import(self,params):
		filenames = params["import_micrograph_ccd_files"]
		
			
		no_dir_names = [get_file_tag(name) for name in filenames]
		
		for name in filenames:
			if name.find("bdb:rawdata#") != -1:
				print "you can't import files that are already in the project raw data directory,",name,"is invalid"
				return
		
		for name in no_dir_names:
			if no_dir_names.count(name) > 1:
				print "you can't use images with the same name (",name,")"
				return
		
		project_db = db_open_dict("bdb:project")
		
		current_project_files = project_db.get("global.micrograph_ccd_filenames",dfl=[])
		cpft = [get_file_tag(file) for file in current_project_files]
		
		# get the number of process operation - the progress dialog reflects image copying and image processing operations
		num_processing_operations = 2 # there is atleast a copy and a disk write
		if params["invert"]: num_processing_operations += 1
		if params["xraypixel"]: num_processing_operations += 1
		if params["thumbs"]:num_processing_operations += 1
		
		
		# now add the files to db (if they don't already exist
		progress = EMProgressDialogModule(self.application(),"Importing files into database...", "Abort import", 0, len(filenames)*num_processing_operations,None)
		progress.qt_widget.show()
		i = 0
		cancelled = False # if the user cancels the import then we must act
		cancelled_dbs = []
		for name in filenames:
			
			tag = get_file_tag(name)
			if tag in cpft:
				print "can't import images have identical tags to those already in the database"
				continue
			
			
			db_name = "bdb:raw_data#"+tag
			if db_check_dict(db_name):
				print "there is already a raw_data database entry for",tag
				continue
			else:
				e = EMData()
				e.read_image(name,0)
				i += 1
				progress.qt_widget.setValue(i)	
				self.application().processEvents()
				e.set_attr("disk_file_name",name)
				
				if params["invert"]:
					e.mult(-1)
					i += 1
					progress.qt_widget.setValue(i)
					self.application().processEvents()
				
				if params["xraypixel"]:
					e.process_inplace("threshold.clampminmax.nsigma",{"nsigma":4,"tomean":True})
					i += 1
					progress.qt_widget.setValue(i)
					self.application().processEvents()
				
				e.write_image(db_name,0)
				cancelled_dbs.append(db_name)
				i += 1
				progress.qt_widget.setValue(i)
				self.application().processEvents()
				current_project_files.append(db_name)
				
				if params["thumbs"]:
					shrink = self.get_thumb_shrink(e.get_xsize(),e.get_ysize())
					thumb = e.process("math.meanshrink",{"n":shrink})
					set_idd_image_entry(db_name,"image_thumb",thumb) # boxer uses the full name
					i += 1
					progress.qt_widget.setValue(i)
					self.application().processEvents()
					
				if progress.qt_widget.wasCanceled():
					cancelled = True
					for data_db in cancelled_dbs: # policy here is to remove only the raw data dbs - the e2boxer thumbnails are tiny and I don't have time...
						db_remove_dict(data_db)
					break
			
				
		progress.qt_widget.setValue(len(filenames))
		progress.qt_widget.close()
		
		if not cancelled:
			project_db["global.micrograph_ccd_filenames"] = current_project_files
			db_close_dict("bdb:project")
		
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
		
	
class MicrographCCDTask(WorkFlowTask):
	
	documentation_string = "This is a list of micrographs or CCD frames that you choose to associate with this project. You can add and remove file names by editing the text entries directly and or by using the browse and clear buttons. In addition to being able to specify images that are stored on your hard drive in the usual way, you can also choose images from EMAN2 style databases."
	
	
	def __init__(self,application):
		WorkFlowTask.__init__(self,application)
		self.window_title = "Project micrographs"

	def get_params(self):
		params = []
		project_db = db_open_dict("bdb:project")
		params.append(ParamDef(name="blurb",vartype="text",desc_short="Raw image data",desc_long="",property=None,defaultunits=MicrographCCDTask.documentation_string,choices=None))
		params.append(ParamDef(name="global.micrograph_ccd_filenames",vartype="url",desc_short="File names",desc_long="The raw data from which particles will be extracted and ultimately refined to produce a reconstruction",property=None,defaultunits=project_db.get("global.micrograph_ccd_filenames",dfl=[]),choices=[]))
		
		db_close_dict("bdb:project")
		return params

	def write_db_entry(self,key,value):
		if key == "global.micrograph_ccd_filenames":
			if value != None:
				stipped_input = [get_file_tag(name) for name in value]
				for name in stipped_input:
					if stipped_input.count(name) > 1:
						print "you can't use images with the same file tag (",name,")"
						return

				new_names = []
		
				# now just update the static list of project file names
				e = EMData()
				read_header_only = True
				for i,name in enumerate(value):
					
					if os.path.exists(name) or db_check_dict(name):
						try:
							a = e.read_image(name,0,read_header_only)
						except: continue
						
						new_names.append(name)
						
				project_db = db_open_dict("bdb:project")
				project_db["global.micrograph_ccd_filenames"] = new_names
				db_close_dict("bdb:project")
		else:  pass


class ParticleWorkFlowTask(WorkFlowTask):
	'''
	Encapsulates some functionality  common to the particle based work flow tasks
	Such tasks should inherit from this class, not the the WorkFlowTask
	'''
	def __init__(self,application):
		WorkFlowTask.__init__(self,application)
	
	
	def get_particle_db_names_versatile(self,end_string="_ptcls",strip_end=False):
		'''
		
		'''
		if os.path.exists("particles/EMAN2DB"):
			dbs = db_list_dicts("bdb:particles#")
			ptcl_dbs = []
			for db in dbs:
				db_strip = get_file_tag(db)
				if len(db_strip) > len(end_string)-1:
					if db_strip[-len(end_string):] == end_string:
						if strip_end:
							ptcl_dbs.append(db_strip[:-len(end_string)])
						else:
							ptcl_dbs.append(db_strip)
			return ptcl_dbs
		else: return []
	
	def get_particle_db_names(self,strip_ptcls=True):
		'''
		returns a list of particle databases in the project
		for instance if these dbs files exist in the particle bdb directory:
		image_001_ptcls.bdb
		image_1000_a_ptcls.bdb
		image_001_ctf.bdbs
		then this function will return [ image_001, image_1000_a]
		Note the if strip_ptcls is False then the return list would instead be [ image_001_ptcls, image_1000_a_ptcls]
		If there are no such dbs then an empty list is returned
		'''
		return self.get_particle_db_names_versatile(end_string="_ptcls",strip_end=strip_ptcls)

	def get_particle_and_project_names(self):
		'''
		
		'''
		ptcl_dbs = self.get_particle_db_names()
		
		project_db = db_open_dict("bdb:project")
		project_files = project_db.get("global.micrograph_ccd_filenames",dfl=[])
		for name in project_files:
			stripped = get_file_tag(name)
			if stripped not in ptcl_dbs:
				ptcl_dbs.append(stripped)
		db_close_dict("bdb:project")
		return ptcl_dbs
	
	def get_num_particles(self,project_files_names):
		ptcl_dbs = self.get_particle_db_names()
		vals = []
		for name in project_files_names:
			name_strip = get_file_tag(name)
			if name_strip in ptcl_dbs:
				db_name = "bdb:particles#"+name_strip+"_ptcls"
				if db_check_dict(db_name):
					db = db_open_dict(db_name)
					vals.append(db["maxrec"]+1)
					db_close_dict(db_name)
				else:
					vals.append("")
			else:
				vals.append("")
		return vals
	
	def get_particle_dims(self,project_names):
		if not os.path.exists("particles/EMAN2DB"):
			vals = [ "" for i in range(len(project_names))]
			return vals
		else:
			vals = []
			for name in project_names:
				name_strip = get_file_tag(name)
				db_name = "bdb:particles#"+name_strip+"_ptcls"
				if db_check_dict(db_name):
					db = db_open_dict(db_name)
					hdr = db.get_header(0)
					vals.append(str(hdr["nx"])+'x'+str(hdr["ny"])+'x'+str(hdr["nz"]))
					db_close_dict(db_name)
				else: vals.append("")
			return vals
		
	def get_particle_dims_direct(self,file_names):
		'''
		Iterates through the given names opening the bdb:particles#name
		and returning the dimensions as string e.g. 128x128x1, in a list
		You can use this function if the file_names are exactly those that exist in the particles BDB
		otherwise use get_particle_dims
		'''
		vals = []
		for name in file_names:
			db_name = "bdb:particles#"+name
			if db_check_dict(db_name):
				db = db_open_dict(db_name)
				hdr = db.get_header(0)
				vals.append(str(hdr["nx"])+'x'+str(hdr["ny"])+'x'+str(hdr["nz"]))
				db_close_dict(db_name)
			else: vals.append("No data")
		return vals

	def get_num_particles_direct(self,file_names):
		'''
		Iterates through the given names opening the bdb:particles#name
		and returning the number of particles in each db, in a list. If for some reason the dictionary doesn't exist will
		return a -1 in the list instead
		You can use this function if the file_names are exactly those that exist in the particles BDB
		otherwise use get_num_particles
		'''
		vals = []
		for name in file_names:
			db_name = "bdb:particles#"+name
			if db_check_dict(db_name):
				db = db_open_dict(db_name)
				vals.append(db["maxrec"]+1)
				db_close_dict(db_name)
			else:
				vals.append(-1)
		return vals

	def get_project_particle_param_table(self):
		'''
		Use the names in the global.micrograph_ccd_filenames to build a table showing the corresponding and  current number of boxed particles in the particles directory, and also lists their dimensions
		Puts the information in a ParamTable.
		'''
		project_db = db_open_dict("bdb:project")	
		project_names = project_db.get("global.micrograph_ccd_filenames",dfl=[])
		
		return self.__make_particle_param_table(project_names)
	
	def get_particle_param_table(self):
		'''
		Inspects the particle databases in the particles directory, gathering their names, number of particles and dimenions. Puts the information in a ParamTable.
		'''
		project_names = self.get_particle_db_names(strip_ptcls=True)
		
		return self.__make_particle_param_table(project_names)
	
	def __make_particle_param_table(self,project_names):
		'''
		Functionality used in two places (directly above)
		'''
		
		num_boxes = self.get_num_particles(project_names)
		dimensions = self.get_particle_dims(project_names)
		
		pnames = ParamDef(name="global.micrograph_ccd_filenames",vartype="stringlist",desc_short="File names",desc_long="The raw data from which particles will be extracted and ultimately refined to produce a reconstruction",property=None,defaultunits=None,choices=project_names)
		pboxes = ParamDef(name="Num boxes",vartype="intlist",desc_short="Particles on disk",desc_long="The number of box images stored for this image in the database",property=None,defaultunits=None,choices=num_boxes)
		pdims = ParamDef(name="Dimensions",vartype="stringlist",desc_short="Particle dims",desc_long="The dimensions of the particle images",property=None,defaultunits=None,choices=dimensions)
		
		p = ParamTable(name="filenames",desc_short="Choose a subset of these images",desc_long="")
		p.append(pnames)
		p.append(pboxes)
		p.append(pdims)
		
		return p, len(num_boxes)
	
	def get_total_particles_project_exclusive(self,tag="_ptcls"):
		'''
		A way to get the total number of particles that have a certain 
		'''
		project_db = db_open_dict("bdb:project")	
		project_names = project_db.get("global.micrograph_ccd_filenames",dfl=[])
		stripped_project_names = [get_file_tag(name) for name in project_names ]
		particle_names = self.get_particle_db_names_versatile(tag,strip_end=True)
		
		good_names = []
		for name in particle_names:
			if name in stripped_project_names:
				good_names.append(name+tag)
		
		n = 0
		for name in good_names:
			db_name = "bdb:particles#"+name
			if db_check_dict(db_name):
				pt_db = db_open_dict(db_name)
				if pt_db.has_key("maxrec"):
					val = pt_db["maxrec"]
					if val != None:
						n += val+1 # maxrec is always one less than the actual number stored
						
		return n
	
	def get_total_particles(self,tag="_ptcls"):
		'''
		A way to get the total number of particles that have a certain ending
		tag = "_ptcls"
		tag = "_ptcls_ctf_wiener"
		tag = "_ptcls_ctf_phase" ALL WORK
		'''
		particle_names = self.get_particle_db_names_versatile(tag,strip_end=False)
		
		n = 0
		for name in particle_names:
			db_name = "bdb:particles#"+name
			if db_check_dict(db_name):
				pt_db = db_open_dict(db_name)
				if pt_db.has_key("maxrec"):
					val = pt_db["maxrec"]
					if val != None:
						n += val+1 # maxrec is always one less than the actual number stored
						
		return n
	
	def get_particle_selection_table(self,tag="_ptcls"):
		'''
		tag = "_ptcls"
		tag = "_ptcls_ctf_wiener"
		tag = "_ptcls_ctf_phase" ALL WORK
		'''
		particle_names = self.get_particle_db_names_versatile(tag,strip_end=False)
		n = []
		dims = []
		for name in particle_names:
			db_name = "bdb:particles#"+name
			act = True
			if db_check_dict(db_name):
				pt_db = db_open_dict(db_name)
				if pt_db.has_key("maxrec"):
					val = pt_db["maxrec"]
					if val != None:
						n.append(val)
						hdr = pt_db.get_header(0)
						dims.append(str(hdr["nx"])+'x'+str(hdr["ny"])+'x'+str(hdr["nz"]))
						act = False
						
			if act:
				n.append("")
				dims.append("")

		p = ParamTable(name="filenames",desc_short="Choose a subset of these images",desc_long="")
			
		pnames = ParamDef(name="names",vartype="stringlist",desc_short="File names",desc_long="The particles that will be used",property=None,defaultunits=None,choices=particle_names)
		pboxes = ParamDef(name="Num boxes",vartype="intlist",desc_short="Particles on disk",desc_long="The number of box images stored for this image in the database",property=None,defaultunits=None,choices=n)
		pdims = ParamDef(name="Dimensions",vartype="stringlist",desc_short="Particle dims",desc_long="The dimensions of the particle images",property=None,defaultunits=None,choices=dims)
		
	
		
		p.append(pnames)
		p.append(pboxes)
		p.append(pdims)
		
		return p,len(pnames)
		
		
class ParticleReportTask(ParticleWorkFlowTask):
	'''
	
	Reports the current status of the particles in the projec
	'''
	
	documentation_string = "This tool is for displaying the particles that are currently associated with this project. This list is generating by inspecting the contents of the project particles directory."
	warning_string = "\n\n\nNOTE: There are no particles currently associated with the project. You can add particles to the project using e2boxer or by importing them directly - see from the list of options associated with this task." 
	def __init__(self,application):
		ParticleWorkFlowTask.__init__(self,application)
		self.window_title = "Project particles"

	def get_params(self):
		params = []
		
		p,n = self.get_particle_param_table()
		if n == 0:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=ParticleReportTask.documentation_string+ParticleReportTask.warning_string,choices=None))
		else:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=ParticleReportTask.documentation_string,choices=None))
			params.append(p)  
		
		return params


class ParticleImportTask(ParticleWorkFlowTask):	
	'''
	
	A task for importing particles into the project particles directory
	'''
	documentation_string = "Use this tool for importing particle images into the particle directory in the project database. Files that you import in this way will be automatically added the list of particle files in the project."
	
	def __init__(self,application):
		ParticleWorkFlowTask.__init__(self,application)
		self.window_title = "Import particles"

	def get_params(self):
		params = []
		project_db = db_open_dict("bdb:project")
		params.append(ParamDef(name="blurb",vartype="text",desc_short="Description",desc_long="",property=None,defaultunits=ParticleImportTask.documentation_string,choices=None))
		params.append(ParamDef(name="import_particle_files",vartype="url",desc_short="Imported particles",desc_long="A list of particle files that have been or will be imported into this project",property=None,defaultunits=[],choices=[]))
		db_close_dict("bdb:project")
		return params


	def write_db_entry(self,k,v):
		if k == "import_particle_files":
			# first make sure that the file tags that will correspond to the imported file names do not conflict with names in the project files
			project_db = db_open_dict("bdb:project")
			project_file_names = project_db.get("global.micrograph_ccd_filenames",dfl=[])
			pfnt = [ get_file_tag(name) for name in project_file_names]
			
			for name in v:
				if get_file_tag(name) in pfnt:
					print "error, you can't import particles that have the same file name as one of the project files - the problem is with:",get_file_tag(name)
					db_close_dict("bdb:project")
					return
			# now check to see if there are duplicated filetags in the incoming list
			# Actually the url form widget may have already dealt with this?
			nt = [ get_file_tag(name) for name in v]
			for name in v:
				if nt.count(get_file_tag(name)) > 1:
					print "error, you can't import particles that have the same file names! The problem is with:",get_file_tag(name)
					db_close_dict("bdb:project")
					return
			
			# now check to see if there isn't already an entry in the particle directory that corresponds to this name
			particle_dbs = self.get_particle_db_names()
			for name in v:
				if get_file_tag(name) in particle_dbs:
					print "error, you can't import particles that already have entries in the particle database! The problem is with:",get_file_tag(name)
					db_close_dict("bdb:project")
					return
			# okay if we make it here we're fine, import that particles
			progress = EMProgressDialogModule(self.application(),"Importing files into database...", "Abort import", 0, len(v),None)
			progress.qt_widget.show()
#			self.application().show_specific(progress)
			read_header_only = True
			a= EMData()
			cancelled_dbs = []
			for i,name in enumerate(v):
				progress.qt_widget.setValue(i)
				self.application().processEvents()
				if (os.path.exists(name) or db_check_dict(name)) and is_2d_image_mx(name):
					
					try:
						a.read_image(name,0,read_header_only)
					except:
						print "error,",name,"is not a valid image. Ignoring"
						continue
					
					tag = get_file_tag(name)
					db_name = "bdb:particles#"+tag+"_ptcls"
					imgs = EMData().read_images(name)
					cancelled_dbs.append(db_name)
					for img in imgs: img.write_image(db_name,-1)
				else:
					print "error,",name,"doesn't exist. Ignoring"
					
				if progress.qt_widget.wasCanceled():
					cancelled = True
					for db_name in cancelled_dbs:
						db_remove_dict(db_name)
					break
			
					
			progress.qt_widget.setValue(len(v))
			#self.application().close_specific(progress)
			progress.qt_widget.close()
				
			
			db_close_dict("bdb:project")

		

class E2BoxerTask(ParticleWorkFlowTask):
	'''
	Provides some common functions for the e2boxer tasks
	'''
	def __init__(self,application):
		ParticleWorkFlowTask.__init__(self,application)
		
	def get_e2boxer_boxes_and_project_particles_table(self):
		'''
		
		Returns a table like this:
		
		|| Project image name || Particles currently in desktop || Dimensions of Particles || Boxes in e2boxer db || Dims of boxes in e2boxer db|
		
		Returns the table, and the the number of entries (p,n)
		if n is zero there are no entries in the table and the calling function can act appropriately
		'''
		
		p,n = self.get_project_particle_param_table() # now p is a ParamTable with rows for as many files as there in the project
		# also, p contains columns with filename | particle number | particle dimensions
		 
		project_db = db_open_dict("bdb:project")	
		project_names = project_db.get("global.micrograph_ccd_filenames",dfl=[])
		
		
		nboxes,dimensions = self.__get_e2boxer_data(project_names)
		pboxes = ParamDef(name="Num boxes",vartype="intlist",desc_short="Boxes in DB",desc_long="The number of boxes stored for this image in the database",property=None,defaultunits=None,choices=nboxes)
		pdims = ParamDef(name="DB Box Dims",vartype="stringlist",desc_short="Dims in DB",desc_long="The dimensions boxes",property=None,defaultunits=None,choices=dimensions)
		
		
		p_reordered = ParamTable(name="filenames",desc_short="Choose a subset of these images",desc_long="") # because I want the boxes in db to come first
		p_reordered.append(p[0])
		p_reordered.append(pboxes)
		p_reordered.extend(p[1:])
		
		
		#p.append(pdims) # don't think this is really necessary
		return p_reordered,len(nboxes)
	
	def get_project_files_that_have_db_boxes_in_table(self):
		project_db = db_open_dict("bdb:project")	
		project_names = project_db.get("global.micrograph_ccd_filenames",dfl=[])
		
		flat_boxes = self.get_num_particles(project_names)
		flat_dims = self.get_particle_dims(project_names)
		
		
		db_boxes,db_dims = self.__get_e2boxer_data(project_names)
		
		
		
		for i in range(len(db_dims)-1,-1,-1):
			if db_boxes[i] == "":
				for data in [project_names,flat_boxes,flat_dims,db_boxes,db_dims]:
					data.pop(i)
		
		
		pnames = ParamDef(name="global.micrograph_ccd_filenames",vartype="stringlist",desc_short="File names",desc_long="The raw data from which particles will be extracted and ultimately refined to produce a reconstruction",property=None,defaultunits=None,choices=project_names)
		pboxes = ParamDef(name="Particles on disk",vartype="intlist",desc_short="Particles on disk",desc_long="The number of box images stored for this image in the database",property=None,defaultunits=None,choices=flat_boxes)
		pdims = ParamDef(name="Particle dimensions",vartype="stringlist",desc_short="Particle dims",desc_long="The dimensions of the particle images",property=None,defaultunits=None,choices=flat_dims)
		pdbboxes = ParamDef(name="Boxes in DB",vartype="intlist",desc_short="Boxes in DB",desc_long="The number of boxes stored for this image in the database",property=None,defaultunits=None,choices=db_boxes)
		pdvdims = ParamDef(name="DB Box Dims",vartype="stringlist",desc_short="Dims in DB",desc_long="The dimensions boxes",property=None,defaultunits=None,choices=db_dims)
		
		p = ParamTable(name="filenames",desc_short="Choose a subset of these images",desc_long="")
		p.append(pnames)
		p.append(pdbboxes)
		p.append(pboxes)
		p.append(pdims)
		
		#p.append(pdvdims) # decided this wasn't necessary
		
		return p,len(db_dims)
		
	def __get_e2boxer_data(self,project_names):
		
		db_name = "bdb:e2boxer.cache"		
		box_maps = {}
		if db_check_dict(db_name):
			e2boxer_db = db_open_dict(db_name)
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


class E2BoxerGenericTask(ParticleWorkFlowTask):
	documentation_string = "Fill me in"
	def __init__(self,application):
		ParticleWorkFlowTask.__init__(self,application)
		self.window_title = "e2ctf"
		self.preferred_size = (480,200)
		
		
	def get_params(self):
		params = []		
#		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2CTFGenericTask.documentation_string,choices=None))
		
		params.append(ParamDef(name="running_mode",vartype="choice",desc_short="Choose your running mode",desc_long="There are three Boxer related task which are generally run in order",property=None,defaultunits="auto",choices=["Interactive boxing", "Autoboxing", "Write output"]))
		
		return params

	def on_form_ok(self,params):
		if params["running_mode"] == "Interactive boxing":
			self.emit(QtCore.SIGNAL("replace_task"),E2BoxerGuiTaskGeneral,"e2boxer interface launcher")
			self.application().close_specific(self.form)
			self.form = None
		elif params["running_mode"] == "Autoboxing":
			self.emit(QtCore.SIGNAL("replace_task"),E2BoxerAutoTaskGeneral,"e2boxer automated boxing")
			self.application().close_specific(self.form)
			self.form = None
		elif params["running_mode"] == "Write output":
			self.emit(QtCore.SIGNAL("replace_task"),E2BoxerOutputTaskGeneral,"e2boxer write output")
			self.application().close_specific(self.form)
			self.form = None	
		else:
			self.application().close_specific(self.form)
			self.form = None
			self.emit(QtCore.SIGNAL("task_idle"))
			
	def write_db_entry(self,key,value):
		pass


class E2BoxerAutoTask(E2BoxerTask):
	'''
	A task for running automated boxing in the project context
	'''
	documentation_string = "Select the images you wish to run autoboxing on and hit OK.\nThis will cause the workflow to spawn processes based on the available CPUs.\nData will be autoboxed using the current autoboxer in the database, which is placed there by e2boxer."
	warning_string = "\n\n\nNOTE: This feature is currently disabled as there is no autoboxing information in the EMAN2 database. You can fix this situation by using e2boxer to interactively box a few images first"
	
	def __init__(self,application):
		E2BoxerTask.__init__(self,application)
		self.window_title = "e2boxer auto"
		self.boxer_module = None # this will actually point to an EMBoxerModule, potentially

	def get_params(self):
		params = []
		
		
		db_name = "bdb:e2boxer.cache"
		fail = False
		if db_check_dict(db_name):
			db = db_open_dict(db_name)
			if not db.has_key("current_autoboxer"): fail = True
		else:
			fail = True
			
		if fail:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="Using e2boxer",desc_long="",property=None,defaultunits=E2BoxerAutoTask.documentation_string+E2BoxerAutoTask.warning_string,choices=None))
		else:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="Using e2boxer",desc_long="",property=None,defaultunits=E2BoxerAutoTask.documentation_string,choices=None))
			p,n = self.get_e2boxer_boxes_and_project_particles_table()
			params.append(p)
	
		return params
			
	def on_form_ok(self,params): 
		for k,v in params.items():
			self.write_db_entry(k,v)
			
		if not params.has_key("filenames") or len(params["filenames"]) == 0:
			self.emit(QtCore.SIGNAL("task_idle"))
			self.application().close_specific(self.form)
			self.form = None
			return
	
		else:
			options = EmptyObject()
			for k,v in params.items():
				setattr(options,k,v)
			
			string_args = []
			bool_args = []
			additional_args = ["--method=Swarm", "--auto=db"]
			temp_file_name = "e2boxer_autobox_stdout.txt"
			self.run_task("e2boxer.py",options,string_args,bool_args,additional_args,temp_file_name)
			self.emit(QtCore.SIGNAL("task_idle"))
			self.application().close_specific(self.form)
			self.form = None

	def write_db_entry(self,key,value):
		if key == "boxsize":
			boxer_project_db = db_open_dict("bdb:e2boxer.project")
			boxer_project_db["working_boxsize"] = value
			db_close_dict("bdb:e2boxer.project")
		else:
			pass
	
class E2BoxerAutoTaskGeneral(E2BoxerAutoTask):
	def __init__(self,application):
		E2BoxerAutoTask.__init__(self,application)
		self.window_title = "e2boxer auto"
		self.boxer_module = None # this will actually point to an EMBoxerModule, potentially

	def get_params(self):
		params = []
		
		
		db_name = "bdb:e2boxer.cache"
		fail = False
		if db_check_dict(db_name):
			db = db_open_dict(db_name)
			if not db.has_key("current_autoboxer"): fail = True
		else:
			fail = True
			
		if fail:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="Using e2boxer",desc_long="",property=None,defaultunits=E2BoxerAutoTask.documentation_string+E2BoxerAutoTask.warning_string,choices=None))
		else:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="Using e2boxer",desc_long="",property=None,defaultunits=E2BoxerAutoTask.documentation_string,choices=None))
			params.append(ParamDef(name="filenames",vartype="url",desc_short="File names",desc_long="The names of the particle files you want to interactively box using e2boxer",property=None,defaultunits=[],choices=[]))
		
		return params

class E2BoxerGuiTask(E2BoxerTask):	
	documentation_string = "Select the images you want to box, enter your boxsize, and hit OK. This will lauch e2boxer and automatically load the selected images for boxing."
	
	warning_string = "\n\n\nNOTE: There are no images currenty associated with the project. Please import or specify which images you want as part of this project in step 1 of the workflow and try again."
	
	def __init__(self,application):
		E2BoxerTask.__init__(self,application)
		self.window_title = "e2boxer interface"
		self.boxer_module = None # this will actually point to an EMBoxerModule, potentially

	def get_params(self):
		params = []
		
		p,n = self.get_e2boxer_boxes_and_project_particles_table()
		
		if n == 0:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="Interactive use of e2boxer",desc_long="",property=None,defaultunits=E2BoxerGuiTask.documentation_string+E2BoxerGuiTask.warning_string,choices=None))
		else:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="Interactive use of e2boxer",desc_long="",property=None,defaultunits=E2BoxerGuiTask.documentation_string,choices=None))
			params.append(p)
			boxer_project_db = db_open_dict("bdb:e2boxer.project")
			params.append(ParamDef(name="boxsize",vartype="int",desc_short="Box size",desc_long="An integer value",property=None,defaultunits=boxer_project_db.get("working_boxsize",dfl=128),choices=[]))
		return params
			
	def on_form_ok(self,params):
		for k,v in params.items():
			self.write_db_entry(k,v)
			
		if not params.has_key("filenames") or len(params["filenames"]) == 0:
			self.emit(QtCore.SIGNAL("task_idle"))
			self.application().close_specific(self.form)
			self.form = None
			return

		else:
			options = EmptyObject()
			for key in params.keys():
				setattr(options,key,params[key])
			options.running_mode = "gui"
			options.method = "Swarm"
			
			
			self.boxer_module = EMBoxerModule(self.application(),options)
			self.emit(QtCore.SIGNAL("gui_running"),"Boxer",self.boxer_module) # The controlled program should intercept this signal and keep the E2BoxerTask instance in memory, else signals emitted internally in boxer won't work
			
			QtCore.QObject.connect(self.boxer_module, QtCore.SIGNAL("module_idle"), self.on_boxer_idle)
			QtCore.QObject.connect(self.boxer_module, QtCore.SIGNAL("module_closed"), self.on_boxer_closed)
			self.application().close_specific(self.form)
			self.boxer_module.show_guis()
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
	
	def on_boxer_idle(self):
		'''
		Presently this means boxer did stuff but never opened any guis, so it's safe just to emit the signal
		'''
		self.boxer_module = None
		self.emit(QtCore.SIGNAL("gui_exit"))

	def write_db_entry(self,key,value):
		if key == "boxsize":
			boxer_project_db = db_open_dict("bdb:e2boxer.project")
			boxer_project_db["working_boxsize"] = value
			db_close_dict("bdb:e2boxer.project")
		else:
			pass

class E2BoxerGuiTaskGeneral(E2BoxerGuiTask):	
	def __init__(self,application):
		E2BoxerTask.__init__(self,application)
		self.window_title = "e2boxer interface"
		self.boxer_module = None # this will actually point to an EMBoxerModule, potentially

	def get_params(self):
		params = []
		params.append(ParamDef(name="blurb",vartype="text",desc_short="Using e2boxer",desc_long="",property=None,defaultunits=E2BoxerGuiTask.documentation_string,choices=None))
		params.append(ParamDef(name="filenames",vartype="url",desc_short="File names",desc_long="The names of the particle files you want to interactively box using e2boxer",property=None,defaultunits=[],choices=[]))
		boxer_project_db = db_open_dict("bdb:e2boxer.project")
		params.append(ParamDef(name="boxsize",vartype="int",desc_short="Box size",desc_long="An integer value",property=None,defaultunits=boxer_project_db.get("working_boxsize",dfl=128),choices=[]))
		
		return params
	
	
class E2BoxerOutputTask(E2BoxerTask):	
	documentation_string = "Select the images you wish to generate output for, enter the box size and normalization etc, and then hit OK.\nThis will cause the workflow to spawn output writing processes using the available CPUs. Note that the bdb option is the preferred output format, in this mode output particles are written directly to the EMAN project database."
	warning_string = "\n\n\nNOTE: There are no boxes currently stored in the database. To rectify this situation use e2boxer to interactively box your images, or alternatively used autoboxing information stored in the database to autobox your images."
	def __init__(self,application):
		E2BoxerTask.__init__(self,application)
		self.window_title = "e2boxer output"
	
	def get_params(self):
		params = []
		
		p,n = self.get_project_files_that_have_db_boxes_in_table()
		if n == 0:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="Using e2boxer",desc_long="",property=None,defaultunits=E2BoxerOutputTask.documentation_string+E2BoxerOutputTask.warning_string,choices=None))
		else:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="Using e2boxer",desc_long="",property=None,defaultunits=E2BoxerOutputTask.documentation_string,choices=None))
			params.append(p)
			self.add_general_params(params)

		return params
	
	def add_general_params(self,params):
		'''
		Functionality used in several places
		'''
		boxer_project_db = db_open_dict("bdb:e2boxer.project")
		pbox = ParamDef(name="boxsize",vartype="int",desc_short="Box size",desc_long="An integer value",property=None,defaultunits=boxer_project_db.get("working_boxsize",dfl=128),choices=[])	
		pfo = ParamDef(name="force",vartype="boolean",desc_short="Force overwrite",desc_long="Whether or not to force overwrite files that already exist",property=None,defaultunits=False,choices=None)
		pwc = ParamDef(name="write_coord_files",vartype="boolean",desc_short="Write box db files",desc_long="Whether or not box db files should be written",property=None,defaultunits=False,choices=None)
		pwb = ParamDef(name="write_box_images",vartype="boolean",desc_short="Write box image files",desc_long="Whether or not box images should be written",property=None,defaultunits=True,choices=None)
		pn =  ParamDef(name="normproc",vartype="string",desc_short="Normalize images",desc_long="How the output box images should be normalized",property=None,defaultunits="normalize.edgemean",choices=["normalize","normalize.edgemean","none"])
		pop = ParamDef(name="outformat",vartype="string",desc_short="Output image format",desc_long="The format of the output box images",property=None,defaultunits="bdb",choices=["bdb","img","hdf"])
		params.append([pbox,pfo])
		params.append([pwc,pwb])
		params.append(pn)
		params.append(pop)
		
			
	def on_form_ok(self,params):
		for k,v in params.items():
			self.write_db_entry(k,v)
			
		if not params.has_key("filenames") or len(params["filenames"]) == 0:
			self.emit(QtCore.SIGNAL("task_idle"))
			self.application().close_specific(self.form)
			self.form = None
			return
	
		else:
			options = EmptyObject()
			for k,v in params.items():
				setattr(options,k,v)
				
			options.just_output=True # this is implicit, it has to happen
			
			string_args = ["normproc","outformat","boxsize"]
			bool_args = ["force","write_coord_files","write_box_images","just_output"]
			additional_args = ["--method=Swarm", "--auto=db"]
			temp_file_name = "e2boxer_autobox_stdout.txt"
			self.run_task("e2boxer.py",options,string_args,bool_args,additional_args,temp_file_name)
			self.emit(QtCore.SIGNAL("task_idle"))
			self.application().close_specific(self.form)
			self.form = None
#	
	def write_db_entry(self,key,value):
		if key == "boxsize":
			boxer_project_db = db_open_dict("bdb:e2boxer.project")
			boxer_project_db["working_boxsize"] = value
			db_close_dict("bdb:e2boxer.project")
		else:
			pass

class E2BoxerOutputTaskGeneral(E2BoxerOutputTask):
	documentation_string = "Write me"
	def __init__(self,application):
		E2BoxerOutputTask.__init__(self,application)
		
	def get_params(self):
		params = []
		params.append(ParamDef(name="blurb",vartype="text",desc_short="Using e2boxer",desc_long="",property=None,defaultunits=E2BoxerOutputTaskGeneral.documentation_string,choices=None))
		
		p = self.get_e2boxer_boxes_table(project_check=False)
		params.append(p)
		
		self.add_general_params(params)
	
#		boxer_project_db = db_open_dict("bdb:e2boxer.project")
#		params.append(ParamDef(name="boxsize",vartype="int",desc_short="Box size",desc_long="An integer value",property=None,defaultunits=boxer_project_db.get("working_boxsize",dfl=128),choices=[]))
		return params
	
	def get_e2boxer_boxes_table(self,project_check=True):
		db_name = "bdb:e2boxer.cache"
		p = ParamTable(name="filenames",desc_short="Current boxes generated by e2boxer",desc_long="")
		names = []
		nboxes = []
		dimensions = []
		
		if project_check:
			project_db = db_open_dict("bdb:project")	
			project_names = project_db.get("global.micrograph_ccd_filenames",dfl=[])
		
		if db_check_dict(db_name):
			e2boxer_db = db_open_dict(db_name)
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
			
		pnames = ParamDef(name="Filenames",vartype="stringlist",desc_short="File names",desc_long="The filenames",property=None,defaultunits=None,choices=names)
		pboxes = ParamDef(name="Num boxes",vartype="intlist",desc_short="Boxes in DB",desc_long="The number of boxes stored for this image in the database",property=None,defaultunits=None,choices=nboxes)
		pdims = ParamDef(name="Dimensions",vartype="stringlist",desc_short="Dimensions",desc_long="The dimensions boxes",property=None,defaultunits=None,choices=dimensions)
		
		p = ParamTable(name="filenames",desc_short="Choose a subset of these images",desc_long="")
		p.append(pnames)
		p.append(pboxes)
		p.append(pdims)
		return p

class E2CTFWorkFlowTask(ParticleWorkFlowTask):
	'''
	Common functionality for E2CTF Work flow taskss
	'''
	def __init__(self,application):
		ParticleWorkFlowTask.__init__(self,application)

	def get_ctf_param_table(self,project_names=None,no_particles=False):
		'''
		particle_files_names should be the return variable of self.get_particle_db_names(strip_ptcls=False), or get_ctf_project_names
		'''
		if project_names == None: # just left this here in case anyone is wandering what to do
			project_names = self.get_particle_db_names(strip_ptcls=False)
		
		defocus,dfdiff,dfang,bfactor,noise = self.get_ctf_info(project_names)
		
		pnames = ParamDef(name="micrograph_ccd_filenames",vartype="stringlist",desc_short="File names",desc_long="The raw data from which particles will be extracted and ultimately refined to produce a reconstruction",property=None,defaultunits=None,choices=project_names)
		pdefocus = ParamDef(name="Defocus",vartype="floatlist",desc_short="Defocus",desc_long="Estimated defocus of the microscope",property=None,defaultunits=None,choices=defocus)
		pbfactor = ParamDef(name="Bfactor",vartype="floatlist",desc_short="B factor",desc_long="Estimated B factor of the microscope",property=None,defaultunits=None,choices=bfactor)
		pnoise = ParamDef(name="Noise",vartype="intlist",desc_short="Noise profile length",desc_long="The number of entries in the noise profile",property=None,defaultunits=None,choices=noise)
		
		num_phase,num_wiener,num_particles,phase_dims,wiener_dims,particle_dims = self.get_ctf_particle_info(project_names)
		pboxes = ParamDef(name="Num boxes",vartype="intlist",desc_short="Particles on disk",desc_long="The number of particles stored for this image in the database",property=None,defaultunits=None,choices=num_particles)
		pphase = ParamDef(name="Num phase",vartype="intlist",desc_short="Phase flipped",desc_long="The number of Wiener filter particles stored for this image in the database",property=None,defaultunits=None,choices=num_phase)
		pwiener = ParamDef(name="Num wiener",vartype="intlist",desc_short="Wienered",desc_long="The number of phase flipped particles stored for this image in the database",property=None,defaultunits=None,choices=num_wiener)
			
		pphasedim = ParamDef(name="phase dim",vartype="stringlist",desc_short="Phase ptcl dims",desc_long="The dimensions of the phase flipped images",property=None,defaultunits=None,choices=phase_dims)
		pwienerdim = ParamDef(name="wiener dim",vartype="stringlist",desc_short="Wiener ptcl dims",desc_long="The dimensions of the Wiener filtered images",property=None,defaultunits=None,choices=wiener_dims)
		pparticledim = ParamDef(name="particle dim",vartype="stringlist",desc_short="Particle dims",desc_long="The dimensions of the particle images",property=None,defaultunits=None,choices=particle_dims)
		
		#print len(num_phase),len(num_wiener),len(num_particles),len(phase_dims),len(wiener_dims),len(particle_dims)
		
		p = ParamTable(name="filenames",desc_short="Current CTF parameters",desc_long="")
		
		p.append(pnames)
		p.append(pdefocus)
		p.append(pbfactor)
		p.append(pnoise)
		if not no_particles: 
			p.append(pboxes)
			p.append(pparticledim)
			p.append(pphase)
			p.append(pphasedim)
			p.append(pwiener)
			p.append(pwienerdim)
		
		
		return p,len(defocus)
	
	def get_ctf_info(self,project_names):
		if db_check_dict("bdb:e2ctf.parms"):
			defocus = []
			dfdiff = []
			dfang = []
			bfactor = []
			noise_profile=[]
			ctf_db = db_open_dict("bdb:e2ctf.parms")
			for name in project_names:
				try:
					vals = ctf_db[name][0]
					ctf = EMAN2Ctf()
					ctf.from_string(vals)
					vals = [ctf.defocus,ctf.dfdiff,ctf.dfang,ctf.bfactor,len(ctf.background)]
				except:
					vals = ["","","","",'']  # only need 5 at the moment
					
				defocus.append(vals[0])
				dfdiff.append(vals[1])
				dfang.append(vals[2])
				bfactor.append(vals[3])
				noise_profile.append(vals[4])
				
			db_close_dict("bdb:e2ctf.parms")

			return defocus,dfdiff,dfang,bfactor,noise_profile
	
		else:
			dummy = ["" for i in range(len(project_names))]
			return dummy,dummy,dummy,dummy,dummy
	def get_ctf_particle_info(self,project_names):
		'''
		Returns three string lists containing entries correpsonding to the number of phase corrected, wiener correct, and regular particles in the database
		
		You could call this function with the names of the particle images in the bdb particles directory, i.e. all pariticle file names would
		look like
		img_001_ptcls
		img_003_ptcls
		And these would correspond to entries in the particles database, ie the following bdb entries would exist and contain particles
		bdb:particles#:img_001_ptcls
		bdb:particles#:img_003_ptcls
		
		But you might also call it with the filenames from the e2ctf.parms directory
		
		'''
		phase = []
		wiener = []
		particles = []
		
		phase_dims = []
		wiener_dims = []
		particle_dims = []
		
		for name in project_names:
			particle_db_name = "bdb:particles#"+name
			wiener_ptcl_db_name = particle_db_name + "_ctf_wiener"
			flip_ptcl_db_name = particle_db_name + "_ctf_flip"
			
			d = {}
			d[particle_db_name] = particles
			d[wiener_ptcl_db_name] = wiener
			d[flip_ptcl_db_name] = phase
			
			for db_name,data_list in d.items():
				if db_check_dict(db_name):
					particle_db = db_open_dict(db_name)
					if particle_db.has_key("maxrec"): 
						data_list.append(particle_db["maxrec"]+1)
					else: data_list.append("")
					db_close_dict(db_name)
				else: data_list.append("")
			
			d = {}
			d[particle_db_name] = particle_dims
			d[wiener_ptcl_db_name] = wiener_dims
			d[flip_ptcl_db_name] = phase_dims
			
			for db_name,data_list in d.items():
				if db_check_dict(db_name):
					particle_db = db_open_dict(db_name)
					try:
						hdr = particle_db.get_header(0)
						data_list.append(str(hdr["nx"])+'x'+str(hdr["ny"])+'x'+str(hdr["nz"]))
					except:
						data_list.append("")
				else: data_list.append("")
					
			
		
		return phase,wiener,particles,phase_dims,wiener_dims,particle_dims
		
	def get_ctf_project_names(self):
		'''
		Inspects the e2ctf.parms database for any images that have ctf parameters
		'''
		if db_check_dict("bdb:e2ctf.parms"):
			ctf_parms_db = db_open_dict("bdb:e2ctf.parms")
			vals = ctf_parms_db.keys()
			#print vals
			if vals == None: vals = []
			db_close_dict("bdb:e2ctf.parms")
			return vals
		else:
			return []
		
	def get_project_particle_names_with_ctf(self):
		'''
		Get the names of particles in the project db that also have an entry in the e2ctf.parms db
		This is useful for the E2CTFGui and Output task
		'''
		
		ptcl_names = self.get_particle_db_names(strip_ptcls=False) # particles in the project directory
		ctf_names = self.get_names_with_ctf_params()
		ctf_ptcl_names = [l[0] for l in ctf_names]
		
		interactable_names = []
		for name in ptcl_names:
			if name in ctf_ptcl_names:
				interactable_names.append(name)
				
		return interactable_names
		
	def get_names_with_ctf_params(self):
		'''
		opens the e2ctf.parms directory and returns all a list of lists like this:
		[[db_name_key, real_image_name],....]
		eg
		[[neg_001,neg_001.hdf],[ptcls_01,bdb:particles#ptcls_01_ptcls],...] etc
		e2ctf is responsible for making sure the last data entry for each image is the original image name (this was first enforced by d.woolford)
		'''
		parms_db = db_open_dict("bdb:e2ctf.parms")
		
		ret = []
		for key,data in parms_db.items():
			ret.append([key,data[-1]]) # parms[-1] should be the original filename
				
		return ret

class CTFReportTask(E2CTFWorkFlowTask):
	
	documentation_string = "This tool is for displaying the currently determined CTF parameters for the particles associated with the project. It also displays the number of phase flipped and/or wiener filtered images corresponding to each particle set."
	warning_string = "\n\n\nNOTE: There are no particles currently associated with the project. Please go to the \"Particles\" task and import/box particles first."
	def __init__(self,application):
		E2CTFWorkFlowTask.__init__(self,application)
		self.window_title = "Project particles"

	def get_params(self):
		params = []
		p,n = self.get_ctf_param_table(None)
		
		if n == 0:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=CTFReportTask.documentation_string+CTFReportTask.warning_string,choices=None))
		else:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=CTFReportTask.documentation_string,choices=None))
			params.append(p)
		return params

	def write_db_entry(self,key,value):
		pass
		
class E2CTFGenericTask(ParticleWorkFlowTask):
	documentation_string = "Fill me in"
	def __init__(self,application):
		ParticleWorkFlowTask.__init__(self,application)
		self.window_title = "e2ctf"
		self.preferred_size = (480,200)
		
		
	def get_params(self):
		params = []		
#		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2CTFGenericTask.documentation_string,choices=None))
		
		params.append(ParamDef(name="running_mode",vartype="choice",desc_short="Choose your running mode",desc_long="There are three CTF related task which are generally run in order",property=None,defaultunits="auto",choices=["auto params", "interactively fine tune", "write output"]))
		
		return params

	def on_form_ok(self,params):
		if params["running_mode"] == "auto params":
			self.emit(QtCore.SIGNAL("replace_task"),E2CTFAutoFitTaskGeneral,"ctf auto fit")
			self.application().close_specific(self.form)
			self.form = None
		elif params["running_mode"] == "interactively fine tune":
			self.emit(QtCore.SIGNAL("replace_task"),E2CTFGuiTaskGeneral,"fine tune ctf")
			self.application().close_specific(self.form)
			self.form = None
		elif params["running_mode"] == "write output":
			self.emit(QtCore.SIGNAL("replace_task"),E2CTFOutputTaskGeneral,"ctf output")
			self.application().close_specific(self.form)
			self.form = None	
		else:
			self.application().close_specific(self.form)
			self.form = None
			self.emit(QtCore.SIGNAL("task_idle"))
			
	def write_db_entry(self,key,value):
		pass

class E2CTFAutoFitTask(E2CTFWorkFlowTask):	
	documentation_string = "Select the particles you wish to generate CTF parameters for, enter the appropriate parameters such as microscope voltage etc, and hit OK.\nThis will cause the workflow to spawn processes based on the available CPUs. Once finished the automatically determined CTF parameters will be stored in the EMAN2 database."
	warning_string = "\n\n\nNOTE: There are no particles currently associated with the project. Please go to the \"Particles\" task and import/box particles first."

	def __init__(self,application):
		E2CTFWorkFlowTask.__init__(self,application)
		self.window_title = "e2ctf auto fit"

	def get_params(self):
		params = []
		p,n= self.get_ctf_param_table(self.get_particle_db_names(strip_ptcls=False),no_particles=True)
		
		if n == 0:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2CTFAutoFitTask.documentation_string+E2CTFAutoFitTask.warning_string,choices=None))
		else:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2CTFAutoFitTask.documentation_string,choices=None))
			params.append(p)
			self.add_general_params(params)

		return params
	
	def add_general_params(self,params):
		project_db = db_open_dict("bdb:project")
		ctf_misc_db = db_open_dict("bdb:e2ctf.misc")
		papix = ParamDef(name="global.apix",vartype="float",desc_short="A/pix for project",desc_long="The physical distance represented by the pixel spacing",property=None,defaultunits=project_db.get("global.apix",dfl=1.1),choices=None)
		pvolt = ParamDef(name="global.microscope_voltage",vartype="float",desc_short="Microscope voltage",desc_long="The operating voltage of the microscope in kilo volts",property=None,defaultunits=project_db.get("global.microscope_voltage",dfl=300),choices=None)
		pcs = ParamDef(name="global.microscope_cs",vartype="float",desc_short="Microscope Cs",desc_long="Microscope spherical aberration constant",property=None,defaultunits=project_db.get("global.microscope_cs",dfl=2.0),choices=None)
		pac = ParamDef(name="working_ac",vartype="float",desc_short="Amplitude contrast",desc_long="The amplitude contrast constant. It is recommended that this value is identical in all of your images.",property=None,defaultunits=ctf_misc_db.get("working_ac",dfl=10),choices=None)
		pos = ParamDef(name="working_oversamp",vartype="int",desc_short="Oversampling",desc_long="If greater than 1, oversampling by this amount will be used when images are being phase flipped and Wiener filtered.",property=None,defaultunits=ctf_misc_db.get("working_oversamp",dfl=1),choices=None)
		pncp = ParamDef(name="global.num_cpus",vartype="int",desc_short="Number of CPUs",desc_long="Number of CPUS available for the project to use",property=None,defaultunits=project_db.get("global.num_cpus",dfl=num_cpus()),choices=None)
		mem = memory_stats()
		
		params.append([papix,pvolt])
		params.append([pcs,pac])
		params.append([pos,pncp])

		db_close_dict("bdb:project")
	
	def get_default_ctf_options(self,params):
		'''
		These are the options required to run pspec_and_ctf_fit in e2ctf.py
		'''
		
		if not params.has_key("filenames") or len(params["filenames"]) == 0:
			return None # this is fine there are no particles assocaited with the project
		
		filenames = params["filenames"]
		
		boxsize = None
		db_file_names = []
		for i,name in enumerate(filenames):
			db_name="bdb:particles#"+name
			db_file_names.append(db_name)
			if not db_check_dict(db_name):
				print "error, can't particle entry doesn't exist for",name,"aborting."
				return None
			
			if boxsize == None:
				db = db_open_dict(db_name)
				hdr = db.get_header(0)
				boxsize = hdr["nx"] # no consideration is given for non square images
				db_close_dict(db_name)
			else:
				db = db_open_dict(db_name)
				hdr = db.get_header(0)
				db_close_dict(db_name)
				if boxsize != hdr["nx"]: # no consideration is given for non square images
					print "error, can't run e2ctf on images with different box sizes. Specifically, I can not deduce the bgmask option for the group"
					return None
		
		if boxsize == None or boxsize < 2:
			print "error, boxsize is less than 2"
			return None
		
		
		options = EmptyObject()
		options.bgmask = boxsize/2
		options.filenames = db_file_names
		self.append_general_options(options,params)
	
		return options
	
	def append_general_options(self,options,params):
		'''
		This is done in more than one place hence this function
		'''
		
		options.nosmooth = False
		options.nonorm = False
		options.autohp = False
		options.invert = False
		options.oversamp = params["working_oversamp"]
		options.ac = params["working_ac"]
		options.apix = params["global.apix"]
		options.cs = params["global.microscope_cs"]
		options.voltage = params["global.microscope_voltage"]
		
	def on_form_ok(self,params):
		for k,v in params.items():
			self.write_db_entry(k,v)
		
		options = self.get_default_ctf_options(params)
		if options != None:
			
			string_args = ["bgmask","oversamp","ac","apix","cs","voltage"]
			bool_args = ["nosmooth","nonorm","autohp","invert"]
			additional_args = ["--auto_fit"]
			temp_file_name = "e2ctf_autofit_stdout.txt"
			self.run_task("e2ctf.py",options,string_args,bool_args,additional_args,temp_file_name)
			
			self.application().close_specific(self.form)
			self.emit(QtCore.SIGNAL("task_idle"))
			
		else:
			self.application().close_specific(self.form)
			self.emit(QtCore.SIGNAL("task_idle"))
		
	def on_form_close(self):
		# this is to avoid a task_idle signal, which would be incorrect if e2boxer is running
		self.emit(QtCore.SIGNAL("task_idle"))

	def write_db_entry(self,key,value):
		if key == "working_ac":
			ctf_misc_db = db_open_dict("bdb:e2ctf.misc")
			ctf_misc_db["working_ac"] = value
			db_close_dict("bdb:e2ctf.misc")
		elif key == "working_oversamp":
			ctf_misc_db = db_open_dict("bdb:e2ctf.misc")
			ctf_misc_db["working_oversamp"] = value
			db_close_dict("bdb:e2ctf.misc")
		else:
			# there are some general parameters that need writing:
			WorkFlowTask.write_db_entry(self,key,value)
			
class E2CTFAutoFitTaskGeneral(E2CTFAutoFitTask):
	'''
	This one has a generic url browser to get the input file names as opposed to the more rigid project based one
	'''
	documentation_string = "Use this tool to use e2ctf to generate ctf parameters for the particles located anywhere on disk. Use the browse buttons to locate the files on disk, enter the fitting parameters such as microscope voltage and angstrom per pixel, and hit OK. \nThis will cause the workflow to spawn processes based on the available CPUs. Once finished the automatically determined CTF parameters will be stored in the EMAN2 database."
	
	def __init__(self,application):
		E2CTFAutoFitTask.__init__(self,application)
		self.window_title = "e2ctf auto fit"
		self.file_check = True
	def get_params(self):
		params = []		
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2CTFAutoFitTask.documentation_string,choices=None))
		
		params.append(ParamDef(name="filenames",vartype="url",desc_short="File names",desc_long="The names of the particle files you want to generate automated ctf parameters for",property=None,defaultunits=[],choices=[]))
		
		self.add_general_params(params)

		return params
	
	def get_default_ctf_options(self,params):
		'''
		These are the options required to run pspec_and_ctf_fit in e2ctf.py
		'''
		
		if not params.has_key("filenames") or len(params["filenames"]) == 0:
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
		options.bgmask = boxsize/2
		options.filenames = filenames
		self.append_general_options(options,params)
		return options

class E2CTFOutputTask(E2CTFWorkFlowTask):	
	documentation_string = "Select the particle data for which you wish to generate phase flipped and/or Wiener filtered output and hit OK.\nThis will cause the workflow to spawn processes based on the available CPUs that write the output into a predefined location in the EMAN2 database.\nNote that the Wiener filtered output images are also phase flipped."
	warning_string = "\n\n\nNOTE: There are no particles associated with the project and/or there are no previously generated CTF parameters for these particles. To establish project particles go to the \"Particles\" task. To generate CTF parameters go to the \"Automated fitting - e2ctf\" task" 
	def __init__(self,application):
		E2CTFWorkFlowTask.__init__(self,application)
		self.window_title = "e2ctf management"

	def get_params(self):
		params = []		
		
		p,n = self.get_ctf_param_table(self.get_project_particle_names_with_ctf())
		if n == 0:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2CTFOutputTask.documentation_string+E2CTFOutputTask.warning_string,choices=None))
		else:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2CTFOutputTask.documentation_string,choices=None))
			params.append(p)
			pwiener = ParamDef(name="wiener",vartype="boolean",desc_short="Wiener",desc_long="Wiener filter your particle images using parameters in the database. Phase flipping will also occur",property=None,defaultunits=False,choices=None)
			pphase = ParamDef(name="phaseflip",vartype="boolean",desc_short="Phase flip",desc_long="Phase flip your particle images using parameters in the database",property=None,defaultunits=False,choices=None)
			
			params.append([pphase,pwiener])
		return params

	def get_default_ctf_options(self,params):
		'''
		These are the options required to run pspec_and_ctf_fit in e2ctf.py, works in e2workflow
		'''
		
		if not params.has_key("filenames") or len(params["filenames"]) == 0:
			return None # this is fine, for example, there were no files to work
		
		
		options = EmptyObject()
	
		filenames = params["filenames"]
#
		db_file_names = []
		for i,name in enumerate(filenames):
			db_name="bdb:particles#"+name
			db_file_names.append(db_name)
			if not db_check_dict(db_name):
				print "error, can't particle entry doesn't exist for",name,"aborting."
				return None
		options.filenames = db_file_names
		options.wiener = params["wiener"]
		options.phaseflip = params["phaseflip"]
#		
		return options
	
	def on_form_ok(self,params):
		for k,v in params.items():
			self.write_db_entry(k,v)

		options = self.get_default_ctf_options(params)
		if options != None and len(options.filenames) > 0 and (options.wiener or options.phaseflip):

			string_args = []
			bool_args = ["wiener","phaseflip"]
			additional_args = []
			temp_file_name = "e2ctf_output_stdout.txt"
			self.run_task("e2ctf.py",options,string_args,bool_args,additional_args,temp_file_name)
			

			self.application().close_specific(self.form)
			self.emit(QtCore.SIGNAL("task_idle"))
		else:
			self.application().close_specific(self.form)
			self.emit(QtCore.SIGNAL("task_idle"))
	
	def on_ctf_closed(self):
		self.emit(QtCore.SIGNAL("gui_exit")) #
		
	def on_form_close(self):
		# this is to avoid a task_idle signal, which would be incorrect if e2boxer is running
		self.emit(QtCore.SIGNAL("task_idle"))

class E2CTFOutputTaskGeneral(E2CTFOutputTask):
	''' This one uses the names in the e2ctf.parms to generate it's table of options, not the particles in the particles directory
	'''
	warning_string = "\n\n\nNOTE: There are no CTF parameters currently stored for any images in the local database. You can change this by running automated fitting with e2ctf."
	
	def __init__(self,application):
		E2CTFOutputTask.__init__(self,application)
		self.window_title = "e2ctf management"

	def get_params(self):
		params = []
		names = self.get_names_with_ctf_params()
		n = [l[0] for l in names]
		p,num = self.get_ctf_param_table(n)
		
		if num == 0:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2CTFOutputTask.documentation_string+E2CTFOutputTaskGeneral.warning_string,choices=None))
		else:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2CTFOutputTask.documentation_string,choices=None))
			params.append(p)
			pwiener = ParamDef(name="wiener",vartype="boolean",desc_short="Wiener",desc_long="Wiener filter your particle images using parameters in the database. Phase flipping will also occur",property=None,defaultunits=False,choices=None)
			pphase = ParamDef(name="phaseflip",vartype="boolean",desc_short="Phase flip",desc_long="Phase flip your particle images using parameters in the database",property=None,defaultunits=False,choices=None)
		
			params.append([pphase,pwiener])
			
		return params
	
	def get_ctf_options(self,params):
		'''
		This is a way to get the ctf optiosn if one is using the "alternate" path, which means just in the context of general use of e2ctf (not the workflow)
		'''
		options = EmptyObject()
		
		if not params.has_key("filenames") or len(params["filenames"]) == 0:
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
			self.run_task("e2ctf.py",options,string_args,bool_args,additional_args,temp_file_name)
			

			self.application().close_specific(self.form)
			self.emit(QtCore.SIGNAL("task_idle"))
		else:
			self.application().close_specific(self.form)
			self.emit(QtCore.SIGNAL("task_idle"))
	
class E2CTFGuiTask(E2CTFWorkFlowTask):	
	documentation_string = "Select the particle data you wish to evaluate/tweak in the e2ctf interactive interface and hit OK. This will launch e2ctf and the selected images will automatically be loaded for viewing. Once inside the e2ctf interface you can save your tweaked parameters to the database using the Save button."
	warning_string = "\n\n\nNOTE: There are no particles associated with the project and/or there are no previously generated CTF parameters for these particles. To establish project particles go to the \"Particles\" task. To generate CTF parameters go to the \"Automated fitting - e2ctf\" task" 
	def __init__(self,application):
		E2CTFWorkFlowTask.__init__(self,application)
		self.window_title = "e2ctf management"
		self.gui = None # will eventually be a e2ctf gui

	def get_params(self):
		
		
		ptcl_names = self.get_particle_db_names(strip_ptcls=False) # particles in the project directory
		ctf_names = self.get_names_with_ctf_params()
		ctf_ptcl_names = [l[0] for l in ctf_names]
		
		p,n = self.get_ctf_param_table(self.get_project_particle_names_with_ctf(),no_particles=True)
		
		params = []		
		if n == 0:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2CTFGuiTask.documentation_string+E2CTFGuiTask.warning_string,choices=None))
		else:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2CTFGuiTask.documentation_string,choices=None))
		  	params.append(p)
		return params
	
	def get_default_ctf_options(self,params):
		'''
		These are the options required to run pspec_and_ctf_fit in e2ctf.py
		'''
		
		if not params.has_key("filenames") or len(params["filenames"]) == 0:
			#print "There are no filenames" # this shouldn't happen
			return None
		
		
		options = EmptyObject()
		filenames = params["filenames"]
#
		db_file_names = []
		for i,name in enumerate(filenames):
			db_name="bdb:particles#"+name
			db_file_names.append(db_name)
			if not db_check_dict(db_name):
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
		
			
			self.gui=GUIctfModule(self.application(),img_sets)
			self.emit(QtCore.SIGNAL("gui_running"), "CTF", self.gui) # so the desktop can prepare some space!
			self.application().close_specific(self.form)
			QtCore.QObject.connect(self.gui,QtCore.SIGNAL("module_closed"), self.on_ctf_closed)
			self.gui.show_guis()
		else:
			self.application().close_specific(self.form)
			self.emit(QtCore.SIGNAL("task_idle"))
	
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
	def __init__(self,application):
		E2CTFGuiTask.__init__(self,application)

	def get_params(self):
		params = []		
		
		names = self.get_names_with_ctf_params()
		n = [l[0] for l in names]
		p,num = self.get_ctf_param_table(n)
		
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
		
		if not params.has_key("filenames") or len(params["filenames"]) == 0:
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


class E2Refine2DReportTask(ParticleWorkFlowTask):
	documentation_string = "This form displays the current sets of reference free class averages currently associated with the project"
	def __init__(self,application):
		ParticleWorkFlowTask.__init__(self,application)
		self.window_title = "Refine 2D class averages "
		
		
	def get_params(self):
		params = []		
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2Refine2DReportTask.documentation_string,choices=None))
		
		return params
			
	def write_db_entry(self,key,value):
		pass
	
class E2Refine2DTask(ParticleWorkFlowTask):
	'''
	Common e2refine2D functionality
	'''
	documentation_string = "This form is a way for the user to supply arguments to and execute e2refine2d.py"
	def __init__(self,application):
		ParticleWorkFlowTask.__init__(self,application)
		self.window_title = "Run e2refine2d"		
		
	def get_general_params(self):
		params = []		
		piter = ParamDef(name="iter",vartype="int",desc_short="Refinement iterations",desc_long="The number of times the e2refine2d svd-based class averaging procedure is iterated",property=None,defaultunits=0,choices=[])
		piterclassav = ParamDef(name="iterclassav",vartype="int",desc_short="Class averaging iterations",desc_long="The number of iterative class average steps that occur in e2classaverage",property=None,defaultunits=2,choices=[])
		pnaliref = ParamDef(name="naliref",vartype="int",desc_short="# alignment references",desc_long="The number of alignment references to use when determining particle orientations",property=None,defaultunits=8,choices=[])
		pnbasisfp = ParamDef(name="nbasisfp",vartype="int",desc_short="# basis fp",desc_long="The number of MSA basis vectors to use when classifying",property=None,defaultunits=5,choices=[])
		
		aligners = dump_aligners_list().keys()
		aligners.append("None")
		cmps = dump_cmps_list().keys()
		cmps.append("None")
		
		psimalign =  ParamDef(name="simalign",vartype="string",desc_short="Aligner",desc_long="The aligner being used",property=None,defaultunits="rotate_translate_flip",choices=aligners)
		psimalignargs =  ParamDef(name="simalignargs",vartype="string",desc_short="",desc_long="Parameters for the aligner, see \"e2help.py aligners\"",property=None,defaultunits="",choices=[])
		
		psimaligncmp =  ParamDef(name="simaligncmp",vartype="string",desc_short="Comparitor",desc_long="The comparitor being used",property=None,defaultunits="dot",choices=cmps)
		psimaligncmpsargs =  ParamDef(name="simaligncmpargs",vartype="string",desc_short="",desc_long="Parameters for the comparitor, see \"e2help.py cmps\"",property=None,defaultunits="",choices=[])	
		
		pnp = ParamDef(name="normproj",vartype="boolean",desc_short="Normalize projection vectors",desc_long="Normalizes each projected vector into the MSA subspace",property=None,defaultunits=False,choices=None)
		
		
		project_db = db_open_dict("bdb:project")
		pncp = ParamDef(name="parallel",vartype="int",desc_short="Number of CPUs",desc_long="Number of CPUS available for e2refine2d to use",property=None,defaultunits=project_db.get("global.num_cpus",dfl=num_cpus()),choices=None)

		
		params.append([piter,piterclassav])
		params.append([pnaliref,pnbasisfp])
		params.append([psimalign,psimalignargs])
		params.append([psimaligncmp,psimaligncmpsargs])
		params.append([pncp,pnp])
		
		return params
			
	def write_db_entry(self,key,value):
		pass
		
	
class E2Refine2DCreateDataSetTask(ParticleWorkFlowTask):
	documentation_string = "Choose the data you wish to use for use for running e2refine2d from the list of options below and hit OK. This will pop up a second form asking you to fill in more details.\n\nNote that usually you should have 4 options to choose from below. If you are not seeing all 4 options it means you should go back in the work flow, import particles, and generate phase flipped and Wiener filtered output." 
	def __init__(self,application):
		ParticleWorkFlowTask.__init__(self,application)
		self.window_title = "e2refine2d - getting starting"
		self.preferred_size = (480,300)
		
	def get_params(self):
		params = []		
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2Refine2DCreateDataSetTask.documentation_string,choices=None))
		
		
		n1 = self.get_total_particles(tag="_ptcls")
		n2 = self.get_total_particles(tag="_ptcls_ctf_flip")
		n3 = self.get_total_particles(tag="_ptcls_ctf_wiener")
		
		choices = []
		if n1 > 0:
			choices.append("Particles ("+str(n1)+")")
		if n2 > 0:
			choices.append("Phase flipped ptcls ("+str(n2)+")")
		if n3 > 0:
			choices.append("Wiener ptcls ("+str(n3)+")")
			
		choices.append("Specify files")
		
		
		params.append(ParamDef(name="particle_set_choice",vartype="choice",desc_long="Choose the particle data set you wish to use to generate a starting data for e2refine2d",desc_short="Choose data",property=None,defaultunits=None,choices=choices))
			
		return params
	
	def on_form_ok(self,params):
		for k,v in params.items():
			self.write_db_entry(k,v)

		choice = params["particle_set_choice"]
		
		if choice[:9] == "Particles":
			self.emit(QtCore.SIGNAL("replace_task"),E2Refine2DCreateParticleSetTask,"e2refine2d options")
			self.application().close_specific(self.form)
			self.form = None
		elif choice[:5] == "Phase":
			self.emit(QtCore.SIGNAL("replace_task"),E2Refine2DWithPhasePtclsTask,"e2refine2d options")
			self.application().close_specific(self.form)
			self.form = None
		elif choice[:6] == "Wiener":
			self.emit(QtCore.SIGNAL("replace_task"),E2Refine2DWithWienPtclsTask,"e2refine2d options")
			self.application().close_specific(self.form)
			self.form = None
		elif choice == "Specify files":
			self.emit(QtCore.SIGNAL("replace_task"),E2Refine2DWithGenericTask,"e2refine2d options")
			self.application().close_specific(self.form)
			self.form = None
		else:
			print "the code has changed since the original author wrote it, something is wrong!!"

	
			
	def write_db_entry(self,key,value):
		pass

class E2Refine2DCreateParticleSetTask(E2Refine2DTask):
	documentation_string = "This form enables the user to create starting data sets for e2refine2d.\nChoose from the list of options below in terms of which data you wish to create the initial data set from."
	def __init__(self,application):
		E2Refine2DTask.__init__(self,application)
		self.window_title = "Create refine 2D starting data set"
		self.end_tag = "_ptcls"
	def get_params(self):
		params = []
		
		if self.end_tag != "generic":
			p,n = self.get_particle_selection_table(tag=self.end_tag)
		else:
			p = ParamDef(name="filenames",vartype="url",desc_short="File names",desc_long="The names of the particle files you want to use as in the input data for e2refine2d.py",property=None,defaultunits=[],choices=[])
			n = 1 # just to fool the next bit, that's all
			
		if n == 0:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2Refine2DCreateParticleSetTask.documentation_string+E2Refine2DCreateParticleSetTask.warning_string,choices=None))
		else:
			params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2Refine2DCreateParticleSetTask.documentation_string,choices=None))
			params.append(p)  
			params.append(ParamDef(name="shrink",vartype="int",desc_short="Shrink factor",desc_long="The the downsampling rate used to shrink the data",property=None,defaultunits=1,choices=[]))
		
		other_params = self.get_general_params()
		
		params.extend(other_params)
		return params

	def on_form_ok(self,params):
		for k,v in params.items():
			self.write_db_entry(k,v)
			
		
		
		
		if not params.has_key("filenames") or len(params["filenames"]) == 0 or params["shrink"] < 1:
			print "there was an error, either you didn't choose any filenames or shrink was less than 1. Try again please"
			self.application().close_specific(self.form)
			self.emit(QtCore.SIGNAL("task_idle"))
			return
		
		if params["shrink"] == 1:
			# if shrink is one then we can just make a vstack with e2bdb
			options = EmptyObject()
			names = ["bdb:particles#"+name for name in params["filenames"]]
			options.filenames = names
			string_args = []
			bool_args = []
			additional_args = ["--makevstack=bdb:r2d#"+numbered_path("start_data",True)]
			temp_file_name = "e2bdb_stdout.txt"
			#self.run_task("e2bdb.py",options,string_args,bool_args,additional_args,temp_file_name)
			self.emit(QtCore.SIGNAL("task_idle"))
			self.application().close_specific(self.form)
			self.form = None
			self.emit(QtCore.SIGNAL("task_idle"))
			self.application().close_specific(self.form)
			self.form = None
		else:
			pass
		
		
		
	def write_db_entry(self,key,value):
		pass

class E2Refine2DWithPhasePtclsTask(E2Refine2DCreateParticleSetTask):
	def __init__(self,application):
		E2Refine2DCreateParticleSetTask.__init__(self,application)
		self.end_tag = "_ptcls_ctf_flip"
		
class E2Refine2DWithWienPtclsTask(E2Refine2DCreateParticleSetTask):
	def __init__(self,application):
		E2Refine2DCreateParticleSetTask.__init__(self,application)
		self.end_tag = "_ptcls_ctf_wiener"
		
class E2Refine2DWithGenericTask(E2Refine2DCreateParticleSetTask):
	def __init__(self,application):
		E2Refine2DCreateParticleSetTask.__init__(self,application)
		self.end_tag = "generic"


if __name__ == '__main__':
	
	from emapplication import EMStandAloneApplication
	em_app = EMStandAloneApplication()
	sprinit = SPRInitTask(em_app)
	window = sprinit.run_form() 
	
	#em_app.show()
	em_app.execute()	