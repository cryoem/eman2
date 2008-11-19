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
from EMAN2db import db_check_dict, db_open_dict,db_remove_dict,db_list_dicts
from EMAN2 import EMData,get_file_tag,EMAN2Ctf,num_cpus,memory_stats
import os
import copy
from emapplication import EMProgressDialogModule
from e2boxer import EMBoxerModule
from e2ctf import pspec_and_ctf_fit,GUIctfModule,write_e2ctf_output,get_gui_arg_img_sets
import subprocess
class EmptyObject:
	'''
	This just because I need an object I can assign attributes to, and object() doesn't seem to work
	'''
	def __init__(self):
		pass

class WorkFlowTask(QtCore.QObject):
	def __init__(self,application):
		QtCore.QObject.__init__(self)
		self.application = application
		self.window_title = "Set me please"
	
	def run_form(self):
		self.form = EMFormModule(self.get_params(),self.application)
		self.form.qt_widget.resize(640,640)
		self.form.setWindowTitle(self.window_title)
		self.application.show_specific(self.form)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_ok"),self.on_form_ok)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_cancel"),self.on_form_cancel)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_close"),self.on_form_close)
		
	def on_form_ok(self,params):
		for k,v in params.items():
			if k != "blurb": self.write_db_entry(k,v)
			
		self.application.close_specific(self.form)
		self.form = None
	
		self.emit(QtCore.SIGNAL("task_idle"))
		
	def on_form_cancel(self):
		self.application.close_specific(self.form)
		self.form = None

		self.emit(QtCore.SIGNAL("task_idle"))
	
	def on_form_close(self):
		self.emit(QtCore.SIGNAL("task_idle"))

	def closeEvent(self,event):
		self.application.close_specific(self.form)
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
	def get_wd(self):
		'''
		Get the working directory, originally introduced to provide a centralized mechanism for accessing the working directory,
		specificially for the purpose of spawning processes. Could be used more generally, however.
		'''
		return os.getcwd()

class ParticleWorkFlowTask(WorkFlowTask):
	'''
	Enncapsulates some functionality  common to the particle based work flow tasks
	Such task should inherit from this class, not the the WorkFlowTask
	'''
	def __init__(self,application):
		WorkFlowTask.__init__(self,application)
		
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
		if os.path.exists("particles/EMAN2DB"):
			dbs = db_list_dicts("bdb:particles#")
			ptcl_dbs = []
			for db in dbs:
				db_strip = get_file_tag(db)
				if len(db_strip) > 5:
					if db_strip[-6:] == "_ptcls":
						if strip_ptcls:
							ptcl_dbs.append(db_strip[:-6])
						else:
							ptcl_dbs.append(db_strip)
			return ptcl_dbs
		else: return []

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
				else:
					vals.append(0)
			else:
				vals.append(0)
		return vals
	
	def get_particle_dims(self,particle_file_names):
		if not os.path.exists("particles/EMAN2DB"):
			vals = [ "" for i in range(len(particle_file_names))]
			return vals
		else:
			vals = []
			for name in particle_file_names:
				name_strip = get_file_tag(name)
				db_name = "bdb:particles#"+name_strip+"_ptcls"
				if db_check_dict(db_name):
					db = db_open_dict(db_name)
					hdr = db.get_header(0)
					vals.append(str(hdr["nx"])+'x'+str(hdr["ny"])+'x'+str(hdr["nz"]))
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
			else:
				vals.append(-1)
		return vals


class ParticleReportTask(ParticleWorkFlowTask):
	
	documentation_string = "This tool is for displaying the particles that are currently associated with this project. You can add particles to the project by importing them or by using e2boxer."
	
	def __init__(self,application):
		ParticleWorkFlowTask.__init__(self,application)
		self.window_title = "Project particles"

	def get_params(self):
		params = []
		
		
		#particle_file_names = self.get_particle_and_project_names()
		particle_file_names = self.get_particle_db_names(strip_ptcls=False)
		
		num_boxes = self.get_num_particles_direct(particle_file_names)
		dimensions = self.get_particle_dims_direct(particle_file_names)
		
		
		pnames = ParamDef(name="global.micrograph_ccd_filenames",vartype="stringlist",desc_short="File names",desc_long="The raw data from which particles will be extracted and ultimately refined to produce a reconstruction",property=None,defaultunits=None,choices=particle_file_names)
		pboxes = ParamDef(name="Num boxes",vartype="intlist",desc_short="Number of boxes",desc_long="The number of box images stored for this image in the database",property=None,defaultunits=None,choices=num_boxes)
		pdims = ParamDef(name="Dimensions",vartype="stringlist",desc_short="Dimensions",desc_long="The dimensions of the particle images",property=None,defaultunits=None,choices=dimensions)
		
		
		p = ParamTable(name="filenames",desc_short="Choose images to box",desc_long="")
		p.append(pnames)
		p.append(pboxes)
		p.append(pdims)
		
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=ParticleReportTask.documentation_string,choices=None))
		
		params.append(p)
		
		#boxer_project_db = db_open_dict("bdb:e2boxer.project")
		#params.append(ParamDef(name="boxsize",vartype="int",desc_short="Box size",desc_long="An integer value",property=None,defaultunits=boxer_project_db.get("working_boxsize",dfl=128),choices=[]))
		#params.append(ParamDef(name="global.imported_micrograph_ccd_files",vartype="url",desc_short="File names",desc_long="The raw data from which particles will be extracted and ultimately refined to produce a reconstruction",property=None,defaultunits=project_db.get("global.imported_micrograph_ccd_files",dfl=[]),choices=[]))
		return params

	
	
#	def get_num_particles(self,particle_file_names):
#		if not os.path.exists("particles/EMAN2DB"):
#			vals = [ 0 for i in range(len(particle_file_names))]
#			return vals
#		else:
#			vals = []
#			for name in particle_file_names:
#				try:
#					db = db_open_dict("bdb:particles#"+name+"_ptcls")
#					vals.append(db["maxrec"]+1)
#				except:	vals.append(0)
#			return vals

	def write_db_entry(self,key,value):
		pass


class ParticleImportTask(ParticleWorkFlowTask):	
	documentation_string = "Use this tool for importing particle images into the particle directory in the project database. Files that you import in this way will be automatically added the list of particle files in the project."
	
	def __init__(self,application):
		ParticleWorkFlowTask.__init__(self,application)
		self.window_title = "Import particles"

	def get_params(self):
		params = []
		project_db = db_open_dict("bdb:project")
		params.append(ParamDef(name="blurb",vartype="text",desc_short="Description",desc_long="",property=None,defaultunits=ParticleImportTask.documentation_string,choices=None))
		params.append(ParamDef(name="import_particle_files",vartype="url",desc_short="Imported particles",desc_long="A list of particle files that have been or will be imported into this project",property=None,defaultunits=[],choices=[]))
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
					return
			# now check to see if there are duplicated filetags in the incoming list
			# Actually the url form widget may have already dealt with this?
			nt = [ get_file_tag(name) for name in v]
			for name in v:
				if nt.count(get_file_tag(name)) > 1:
					print "error, you can't import particles that have the same file names! The problem is with:",get_file_tag(name)
					return
			
			# now check to see if there isn't already an entry in the particle directory that corresponds to this name
			particle_dbs = self.get_particle_db_names()
			for name in v:
				if get_file_tag(name) in particle_dbs:
					print "error, you can't import particles that already have entries in the particle database! The problem is with:",get_file_tag(name)
					return
			# okay if we make it here we're fine, import that particles
			progress = EMProgressDialogModule(self.application,"Importing files into database...", "Abort import", 0, len(v),None)
			self.application.show_specific(progress)
			read_header_only = True
			a= EMData()
			for i,name in enumerate(v):
				progress.qt_widget.setValue(i)
				if os.path.exists(name) or db_check_dict(name):
					
					try:
						a.read_image(name,0,read_header_only)
					except:
						print "error,",name,"is not a valid image. Ignoring"
						continue
					
					tag = get_file_tag(name)
					db_name = "bdb:particles#"+tag+"_ptcls"
					imgs = EMData().read_images(name)
					for img in imgs: img.write_image(db_name,-1)
				else:
					print "error,",name,"doesn't exist. Ignoring"
					
			progress.qt_widget.setValue(len(v))
			self.application.close_specific(progress)
			
				
				# why doesn't this work :(
		

class E2BoxerGuiTask(ParticleWorkFlowTask):	
	documentation_string = "Use this tool to box the selected images using e2boxer. Choose from the list of images and hit ok. This will lauch e2boxer and automatically load the selected images for boxing. Alternatively you may choose the auto_db option to execute automated boxing using information stored in the e2boxer database."
	
	def __init__(self,application):
		ParticleWorkFlowTask.__init__(self,application)
		self.window_title = "e2boxer interface"
		self.boxer_module = None # this will actually point to an EMBoxerModule, potentially

	def get_params(self):
		params = []
		project_db = db_open_dict("bdb:project")
		params.append(ParamDef(name="blurb",vartype="text",desc_short="Using e2boxer",desc_long="",property=None,defaultunits=E2BoxerGuiTask.documentation_string,choices=None))
		
		project_file_names = project_db.get("global.micrograph_ccd_filenames",dfl=[])
		
		num_boxes = self.get_num_particles(project_file_names)
#		if len(num_boxes) > 0:
		pnames = ParamDef(name="micrograph_ccd_filenames",vartype="stringlist",desc_short="File names",desc_long="The raw data from which particles will be extracted and ultimately refined to produce a reconstruction",property=None,defaultunits=None,choices=project_file_names)
		pboxes = ParamDef(name="Num boxes",vartype="intlist",desc_short="Number of boxes",desc_long="The number of box images stored for this image in the database",property=None,defaultunits=None,choices=num_boxes)
		
		
		p = ParamTable(name="filenames",desc_short="Choose images to box",desc_long="")
		p.append(pnames)
		p.append(pboxes)
		params.append(p)
	
		boxer_project_db = db_open_dict("bdb:e2boxer.project")
		params.append(ParamDef(name="boxsize",vartype="int",desc_short="Box size",desc_long="An integer value",property=None,defaultunits=boxer_project_db.get("working_boxsize",dfl=128),choices=[]))
		
		return params
			
	def on_form_ok(self,params):
		for k,v in params.items():
			if k != "blurb": self.write_db_entry(k,v)
			
		if not params.has_key("filenames") or len(params["filenames"]) == 0:
			self.emit(QtCore.SIGNAL("task_idle"))
			self.application.close_specific(self.form)
			self.form = None
			return

		else:
			options = EmptyObject()
			for key in params.keys():
				setattr(options,key,params[key])
			options.running_mode = "gui"
			options.method = "Swarm"
				
			self.boxer_module = EMBoxerModule(self.application,options)
			QtCore.QObject.connect(self.boxer_module, QtCore.SIGNAL("module_idle"), self.on_boxer_idle)
			QtCore.QObject.connect(self.boxer_module, QtCore.SIGNAL("module_closed"), self.on_boxer_closed)
			self.application.close_specific(self.form)
			self.form = None
			self.emit(QtCore.SIGNAL("gui_running")) # The controlled program should intercept this signal and keep the E2BoxerTask instance in memory, else signals emitted internally in boxer won't work
			
	def on_form_close(self):
		# this is to avoid a task_idle signal, which would be incorrect if e2boxer is running
		if self.boxer_module == None:
			self.emit(QtCore.SIGNAL("task_idle"))
		else: pass
	
	def on_boxer_closed(self):
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

		else:
			pass

class CTFWorkFlowTask(ParticleWorkFlowTask):
	def __init__(self,application):
		ParticleWorkFlowTask.__init__(self,application)

	def get_ctf_param_table(self,particle_file_names=None):
		'''
		particle_files_names should be the return variable of self.get_particle_db_names(strip_ptcls=False), or get_ctf_project_names
		'''
		if particle_file_names == None: # just left this here in case anyone is wandering what to do
			particle_file_names = self.get_particle_db_names(strip_ptcls=False)
		
		defocus,dfdiff,dfang,bfactor = self.get_ctf_info(particle_file_names)
		
		pnames = ParamDef(name="micrograph_ccd_filenames",vartype="stringlist",desc_short="File names",desc_long="The raw data from which particles will be extracted and ultimately refined to produce a reconstruction",property=None,defaultunits=None,choices=particle_file_names)
		pdefocus = ParamDef(name="Defocus",vartype="floatlist",desc_short="Defocus",desc_long="Estimated defocus of the microscope",property=None,defaultunits=None,choices=defocus)
		pbfactor = ParamDef(name="Bfactor",vartype="floatlist",desc_short="B factor",desc_long="Estimated B factor of the microscope",property=None,defaultunits=None,choices=bfactor)
		num_phase,num_wiener,num_particles = self.get_ctf_particle_info(particle_file_names)
		pboxes = ParamDef(name="Num boxes",vartype="intlist",desc_short="Particles",desc_long="The number of particles stored for this image in the database",property=None,defaultunits=None,choices=num_particles)
		pwiener = ParamDef(name="Num wiener",vartype="intlist",desc_short="Phase flipped",desc_long="The number of Wiener filter particles stored for this image in the database",property=None,defaultunits=None,choices=num_phase)
		pphase = ParamDef(name="Num phase",vartype="intlist",desc_short="Wienered",desc_long="The number of phase flipped particles stored for this image in the database",property=None,defaultunits=None,choices=num_wiener)
		
		p = ParamTable(name="filenames",desc_short="Current CTF parameters",desc_long="")
		
		p.append(pnames)
		p.append(pboxes)
		p.append(pwiener)
		p.append(pphase)
		p.append(pdefocus)
		p.append(pbfactor)
		
		return p
	
	def get_ctf_info(self,particle_file_names):
		if db_check_dict("bdb:e2ctf.parms"):
			defocus = []
			dfdiff = []
			dfang = []
			bfactor = []
			ctf_db = db_open_dict("bdb:e2ctf.parms")
			for name in particle_file_names:
				try:
					vals = ctf_db[name][0]
					ctf = EMAN2Ctf()
					ctf.from_string(vals)
					vals = [ctf.defocus,ctf.dfdiff,ctf.dfang,ctf.bfactor]
				except:
					vals = [0,0,0,0]  # only need 4 at the moment
					
				defocus.append(vals[0])
				dfdiff.append(vals[1])
				dfang.append(vals[2])
				bfactor.append(vals[3])

			return defocus,dfdiff,dfang,bfactor
	
		else:
			dummy = [0 for i in range(len(particle_file_names))]
			return dummy,dummy,dummy,dummy
	def get_ctf_particle_info(self,particle_file_names):
		'''
		Returns three string lists containing entries that are either "yes" or "no"
		The first question asked is do phase corrected images exist for this file name in the particles db?
		The second question asked is do wiener filtered iamges exist for this file name in the particles db?
		The third question is, do particle images exist for the given file name in the particles db?
		A fourth question could be, do all of the existing particle files have the same number of images in them? That might be something that is 
		useful, at some stage
		
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
		
		for name in particle_file_names:
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
				else: data_list.append("")
		
		return phase,wiener,particles
#		else:
#			print "no db for particles"
#			dummy = [-1 for i in range(len(particle_file_names))]
#			return dummy,dummy,dummy
		
		
		
	def get_ctf_project_names(self):
		'''
		Inspects the e2ctf.parms database for any images that have ctf parameters
		'''
		if db_check_dict("bdb:e2ctf.parms"):
			ctf_parms_db = db_open_dict("bdb:e2ctf.parms")
			vals = ctf_parms_db.keys()
			print vals
			if vals == None: vals = []
			return vals
		else:
			print "empty, empty"
			return []

class E2CTFAutoFitTask(CTFWorkFlowTask):	
	documentation_string = "Use this tool to use e2ctf to generate ctf parameters for the particles located in the project particle directory"
	
	def __init__(self,application):
		CTFWorkFlowTask.__init__(self,application)
		self.window_title = "e2ctf auto fit"
		self.options = None # will enventually store e2ctf options
		self.gui = None # will eventually be a e2ctf gui

	def get_params(self):
		params = []		
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2CTFAutoFitTask.documentation_string,choices=None))
		
		params.append(self.get_ctf_param_table(self.get_particle_db_names(strip_ptcls=False)))
		
		project_db = db_open_dict("bdb:project")
		ctf_misc_db = db_open_dict("bdb:e2ctf.misc")
		papix = ParamDef(name="global.apix",vartype="float",desc_short="A/pix for project",desc_long="The physical distance represented by the pixel spacing",property=None,defaultunits=project_db.get("global.apix",dfl=1.1),choices=None)
		pvolt = ParamDef(name="global.microscope_voltage",vartype="float",desc_short="Microscope voltage",desc_long="The operating voltage of the microscope",property=None,defaultunits=project_db.get("global.microscope_voltage",dfl=300),choices=None)
		pcs = ParamDef(name="global.microscope_cs",vartype="float",desc_short="Microscope Cs",desc_long="Microscope spherical aberration constant",property=None,defaultunits=project_db.get("global.microscope_cs",dfl=2.0),choices=None)
		pac = ParamDef(name="working_ac",vartype="float",desc_short="Amplitude contrast",desc_long="The amplitude contrast constant. It is recommended that this value is identical in all of your images.",property=None,defaultunits=ctf_misc_db.get("working_ac",dfl=10),choices=None)
		pos = ParamDef(name="working_oversamp",vartype="int",desc_short="Oversampling",desc_long="If greater than 1, oversampling by this amount will be used when images are being phase flipped and Wiener filtered.",property=None,defaultunits=ctf_misc_db.get("working_oversamp",dfl=1),choices=None)
		pncp = ParamDef(name="global.num_cpus",vartype="int",desc_short="Number of CPUs",desc_long="Number of CPUS available for the project to use",property=None,defaultunits=project_db.get("global.num_cpus",dfl=num_cpus()),choices=None)
		mem = memory_stats()
		
		params.append([papix,pvolt,pcs])
		params.append([pac,pos,pncp])

		return params
	
	def get_default_ctf_options(self,params):
		'''
		These are the options required to run pspec_and_ctf_fit in e2ctf.py
		'''
		
		if not params.has_key("filenames") and len(params["filenames"]) == 0:
			print "there is an internal error. You are asking for the default ctf options but there are no filenames to deduce the bgmask option from" # this shouldn't happen
			return None
		
		
		options = EmptyObject()
		options.nosmooth = False
		options.nonorm = False
		options.autohp = False
		options.invert = False
		
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
			else:
				db = db_open_dict(db_name)
				hdr = db.get_header(0)
				if boxsize != hdr["nx"]: # no consideration is given for non square images
					print "error, can't run e2ctf on images with different box sizes. Specifically, I can not deduce the bgmask option for the group"
					return None
		
		if boxsize == None or boxsize < 2:
			print "error, boxsize is less than 2"
			return None
		
		options.bgmask = boxsize/2
		options.oversamp = params["working_oversamp"]
		options.ac = params["working_ac"]
		options.apix = params["global.apix"]
		options.cs = params["global.microscope_cs"]
		options.voltage = params["global.microscope_voltage"]
		options.filenames = db_file_names
		
		return options

	
	def on_form_ok(self,params):
		for k,v in params.items():
			if k != "blurb": self.write_db_entry(k,v)
		
		self.options = self.get_default_ctf_options(params)
		if self.options != None:
			options = self.options
			project_db = db_open_dict("bdb:project")
			ncpu = project_db.get("global.num_cpus",dfl=num_cpus())
			cf = float(len(options.filenames))/float(ncpu) # common factor
			for n in range(ncpu):
				b = int(n*cf)
				t = int(n+1*cf)
				if n == (ncpu-1):
					t = len(options.filenames) # just make sure of it, round off error could 
				
				if b == t:
					print "hmmmm b equals t"
					continue # it's okay this happens when there are more cpus than there are filenames	
				filenames = options.filenames[b:t]
									
				args = ["e2ctf.py"]
		
				for name in filenames:
					args.append(name)
				
				for string in ["bgmask","oversamp","ac","apix","cs","voltage"]:
					args.append("--"+string+"="+str(getattr(options,string)))
	
				# okay the user can't currently change these, but in future the option might be there
				for string in ["nosmooth","nonorm","autohp","invert"]:
					# these are all booleans so the following works:
					if getattr(options,string):
						args.append("--"+string)
				args.append("--auto_fit")
				#print "command is ",program
				#for i in args: print i
				
				#print args
				file = open("e2ctf_autofit_stdout.txt","w+")
				process = subprocess.Popen(args,stdout=file,stderr=subprocess.STDOUT)
				print "started",process.pid
				self.emit(QtCore.SIGNAL("process_started"),process.pid)
			
			self.application.close_specific(self.form)
			self.emit(QtCore.SIGNAL("task_idle"))
			
		else:
			self.application.close_specific(self.form)
			self.emit(QtCore.SIGNAL("task_idle"))
			
	def on_ctf_closed(self):
		self.gui = None
		write_e2ctf_output(self.options)
		self.options = None
		self.emit(QtCore.SIGNAL("gui_exit")) #
		
	def on_form_close(self):
		# this is to avoid a task_idle signal, which would be incorrect if e2boxer is running
		if self.gui == None:
			self.emit(QtCore.SIGNAL("task_idle"))
		else: pass

	def write_db_entry(self,key,value):
		if key == "working_ac":
			ctf_misc_db = db_open_dict("bdb:e2ctf.misc")
			ctf_misc_db["working_ac"] = value
		elif key == "working_oversamp":
			ctf_misc_db = db_open_dict("bdb:e2ctf.misc")
			ctf_misc_db["working_oversamp"] = value
		else:
			# there are some general parameters that need writing:
			WorkFlowTask.write_db_entry(self,key,value)

class E2CTFOutputTask(CTFWorkFlowTask):	
	documentation_string = "Use this tool to use e2ctf to generate ctf parameters for the particles located in the project particle directory"
	
	def __init__(self,application):
		CTFWorkFlowTask.__init__(self,application)
		self.window_title = "e2ctf management"
		self.options = None # will enventually store e2ctf options
		self.gui = None # will eventually be a e2ctf gui

	def get_params(self):
		params = []		
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2CTFOutputTask.documentation_string,choices=None))
		
		params.append(self.get_ctf_param_table(self.get_particle_db_names(strip_ptcls=False)))
		
		pwiener = ParamDef(name="wiener",vartype="boolean",desc_short="Wiener + phase flip",desc_long="Wiener filter your particle images using parameters in the database. Phase flipping will also occur",property=None,defaultunits=False,choices=None)
		pphase = ParamDef(name="phaseflip",vartype="boolean",desc_short="Phase flip",desc_long="Phase flip your particle images using parameters in the database",property=None,defaultunits=False,choices=None)
		
		params.append(pwiener)
		params.append(pphase)
		return params
	
	def get_default_ctf_options(self,params):
		'''
		These are the options required to run pspec_and_ctf_fit in e2ctf.py
		'''
		
		if not params.has_key("filenames") and len(params["filenames"]) == 0:
			print "there is an internal error. You are asking for the default ctf options but there are no filenames to deduce the bgmask option from" # this shouldn't happen
			return None
		
		
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
			if k != "blurb": self.write_db_entry(k,v)

		options = self.get_default_ctf_options(params)
		if options != None and len(options.filenames) > 0 and (options.wiener or options.phaseflip):
			
			project_db = db_open_dict("bdb:project")
			ncpu = project_db.get("global.num_cpus",dfl=num_cpus())
			cf = float(len(options.filenames))/float(ncpu) # common factor
			for n in range(ncpu):
				b = int(n*cf)
				t = int(n+1*cf)
				if n == (ncpu-1):
					t = len(options.filenames) # just make sure of it, round off error could be problematic
				
				if b == t:
					continue # it's okay this happens when there are more cpus than there are filenames	
				filenames = options.filenames[b:t]
									
				args = ["e2ctf.py"]
		
				for name in filenames:
					args.append(name)
					
				if options.wiener:
					args.append("--wiener")
					
				if options.phaseflip:
					args.append("--phaseflip")

				file = open("e2ctf_output_stdout.txt","w+")
				process = subprocess.Popen(args,stdout=file,stderr=subprocess.STDOUT)
				print "started",process.pid
				self.emit(QtCore.SIGNAL("process_started"),process.pid)
			
			self.application.close_specific(self.form)
			self.emit(QtCore.SIGNAL("task_idle"))
		else:
			self.application.close_specific(self.form)
			self.emit(QtCore.SIGNAL("task_idle"))
	
	def on_ctf_closed(self):
		self.gui = None
		self.options = None
		self.emit(QtCore.SIGNAL("gui_exit")) #
		
	def on_form_close(self):
		# this is to avoid a task_idle signal, which would be incorrect if e2boxer is running
		if self.gui == None:
			self.emit(QtCore.SIGNAL("task_idle"))
		else: pass

	def write_db_entry(self,key,value):
		if key == "working_ac":
			ctf_misc_db = db_open_dict("bdb:e2ctf.misc")
			ctf_misc_db["working_ac"] = value
		elif key == "working_oversamp":
			ctf_misc_db = db_open_dict("bdb:e2ctf.misc")
			ctf_misc_db["working_oversamp"] = value
		else:
			# there are some general parameters that need writing:
			WorkFlowTask.write_db_entry(self,key,value)
	

class E2CTFGuiTask(CTFWorkFlowTask):	
	documentation_string = "Use this tool to use e2ctf to generate ctf parameters for the particles located in the project particle directory"
	
	def __init__(self,application):
		CTFWorkFlowTask.__init__(self,application)
		self.window_title = "e2ctf management"
		self.options = None # will enventually store e2ctf options
		self.gui = None # will eventually be a e2ctf gui

	def get_params(self):
		params = []		
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2CTFGuiTask.documentation_string,choices=None))
		
		params.append(self.get_ctf_param_table(self.get_particle_db_names(strip_ptcls=False)))
		return params
	
	def get_default_ctf_options(self,params):
		'''
		These are the options required to run pspec_and_ctf_fit in e2ctf.py
		'''
		
		if not params.has_key("filenames") and len(params["filenames"]) == 0:
			print "there is an internal error. You are asking for the default ctf options but there are no filenames to deduce the bgmask option from" # this shouldn't happen
			return None
		
		
		options = EmptyObject()
#		options.nosmooth = False
#		options.nonorm = False
#		options.autohp = False
#		options.invert = False
#		
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
#		
		return options

	
	def on_form_ok(self,params):
		for k,v in params.items():
			if k != "blurb": self.write_db_entry(k,v)

		self.options = self.get_default_ctf_options(params)
		if self.options != None and len(self.options.filenames) > 0:
			
			img_sets = get_gui_arg_img_sets(self.options.filenames)
		
			self.gui=GUIctfModule(self.application,img_sets)
			self.application.show_specific(self.gui)
			self.application.close_specific(self.form)
			QtCore.QObject.connect(self.gui.qt_widget,QtCore.SIGNAL("module_closed"), self.on_ctf_closed())
			self.emit(QtCore.SIGNAL("gui_running")) # 
		else:
			self.application.close_specific(self.form)
			self.emit(QtCore.SIGNAL("task_idle"))
	
	def on_ctf_closed(self):
		self.gui = None
		self.options = None
		self.emit(QtCore.SIGNAL("gui_exit")) #
		
	def on_form_close(self):
		# this is to avoid a task_idle signal, which would be incorrect if e2boxer is running
		if self.gui == None:
			self.emit(QtCore.SIGNAL("task_idle"))
		else: pass

	def write_db_entry(self,key,value):
		if key == "working_ac":
			ctf_misc_db = db_open_dict("bdb:e2ctf.misc")
			ctf_misc_db["working_ac"] = value
		elif key == "working_oversamp":
			ctf_misc_db = db_open_dict("bdb:e2ctf.misc")
			ctf_misc_db["working_oversamp"] = value
		else:
			# there are some general parameters that need writing:
			WorkFlowTask.write_db_entry(self,key,value)
	
class CTFReportTask(CTFWorkFlowTask):
	
	documentation_string = "This tool is for displaying the currently determined CTF parameters for the particles located in the project particle directory."
	
	def __init__(self,application):
		CTFWorkFlowTask.__init__(self,application)
		self.window_title = "Project particles"

	def get_params(self):
		params = []
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=CTFReportTask.documentation_string,choices=None))
		
		params.append(self.get_ctf_param_table(self.get_ctf_project_names()))
		return params

	def write_db_entry(self,key,value):
		pass


class MicrographCCDImportTask(WorkFlowTask):	
	documentation_string = "Use this tool for importing flat files into the raw_data directory in the project database. Files that you import in this way will be automatically added the list of files in the project."
	
	def __init__(self,application):
		WorkFlowTask.__init__(self,application)
		self.window_title = "Import micrographs"

	def get_params(self):
		params = []
		project_db = db_open_dict("bdb:project")
		params.append(ParamDef(name="blurb",vartype="text",desc_short="Importing image data",desc_long="",property=None,defaultunits=MicrographCCDImportTask.documentation_string,choices=None))
		params.append(ParamDef(name="import_micrograph_ccd_files",vartype="url",desc_short="File names",desc_long="The raw data from which particles will be extracted and ultimately refined to produce a reconstruction",property=None,defaultunits=[],choices=[]))
		pinvert = ParamDef(name="invert",vartype="boolean",desc_short="Invert",desc_long="Tick this if you want eman2 to invert your images while importing",property=None,defaultunits=False,choices=None)
		pxray = ParamDef(name="xraypixel",vartype="boolean",desc_short="X-ray pixel",desc_long="Tick this if you want eman2 to automatically filter out X-ray pixels while importing",property=None,defaultunits=False,choices=None)
		pthumbnail = ParamDef(name="thumbs",vartype="boolean",desc_short="Thumbnails",desc_long="Tick this if you want eman2 to automatically generate thumbnails for your images. This will save time at later stages in the project",property=None,defaultunits=True,choices=None)
		
		params.append([pinvert,pxray,pthumbnail])
		return params

	def write_db_entry(self,key,value):
		if key == "import_micrograph_ccd_files":
			
			no_dir_names = [get_file_tag(name) for name in value]
			
			for name in value:
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

			# now add the files to db (if they don't already exist
			progress = EMProgressDialogModule(self.application,"Importing files into database...", "Abort import", 0, len(value),None)
			self.application.show_specific(progress)
			for i,name in enumerate(value):
				progress.qt_widget.setValue(i)
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
					e.set_attr("disk_file_name",name)
					e.write_image(db_name,0)
					raw_data_db = db_open_dict(db_name)
					current_project_files.append(db_name)
				
				# why doesn't this work :(
				#print progress.qt_widget.wasCanceled()
				#if progress.qt_widget.wasCanceled():
					#print "it was cancelled"
			progress.qt_widget.setValue(len(value))
			self.application.close_specific(progress)
			
			project_db["global.micrograph_ccd_filenames"] = current_project_files
			
		else:
			print "unknown key:",key,"this object is",self
			
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
		else:
			print "unknown key:",key,"this object is",self
					
	

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
		pvolt = ParamDef(name="global.microscope_voltage",vartype="float",desc_short="Microscope voltage",desc_long="The operating voltage of the microscope",property=None,defaultunits=project_db.get("global.microscope_voltage",dfl=300),choices=None)
		pcs = ParamDef(name="global.microscope_cs",vartype="float",desc_short="Microscope Cs",desc_long="Microscope spherical aberration constant",property=None,defaultunits=project_db.get("global.microscope_cs",dfl=2.0),choices=None)
		pncp = ParamDef(name="global.num_cpus",vartype="int",desc_short="Number of CPUs",desc_long="Number of CPUS available for the project to use",property=None,defaultunits=project_db.get("global.num_cpus",dfl=num_cpus()),choices=None)
		mem = memory_stats()
		pmem = ParamDef(name="global.memory_available",vartype="float",desc_short="Memory usage ("+str(mem[0])+ " Gb total)",desc_long="The total amount of system memory you want to make available to the project in gigabytes",property=None,defaultunits=project_db.get("global.memory_available",dfl=mem[1]),choices=None)
		params.append(papix)
		params.append(pvolt)
		params.append(pcs)
		params.append(pncp)
		params.append(pmem)
		
		return params

	def write_db_entry(self,key,value):
		WorkFlowTask.write_db_entry(self,key,value)
		
		
if __name__ == '__main__':
	
	from emapplication import EMStandAloneApplication
	em_app = EMStandAloneApplication()
	sprinit = SPRInitTask(em_app)
	window = sprinit.run_form() 
	
	#em_app.show()
	em_app.execute()	