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
from EMAN2 import EMData,get_file_tag
import os
import copy
from emapplication import EMProgressDialogModule
from e2boxer import EMBoxerModule

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
		self.form.setWindowTitle(self.window_title)
		self.application.show_specific(self.form)
		QtCore.QObject.connect(self.form.widget,QtCore.SIGNAL("emform_ok"),self.on_form_ok)
		QtCore.QObject.connect(self.form.widget,QtCore.SIGNAL("emform_cancel"),self.on_form_cancel)
		QtCore.QObject.connect(self.form.widget,QtCore.SIGNAL("emform_close"),self.on_form_close)
		
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
		#self.emit(QtCore.SIGNAL("task_idle"))

class ParticleImportTask(WorkFlowTask):	
	documentation_string = "Use this tool for importing particle image into the project database under the particle directory. Files that you import in this way will be automatically added the list of particle files in the project."
	
	def __init__(self,application):
		WorkFlowTask.__init__(self,application)
		self.window_title = "Import particles"

	def get_params(self):
		params = []
		project_db = db_open_dict("bdb:project")
		params.append(ParamDef(name="blurb",vartype="text",desc_short="Description",desc_long="",property=None,defaultunits=ParticleImportTask.documentation_string,choices=None))
		params.append(ParamDef(name="global.imported_particle_files",vartype="url",desc_short="Imported particles",desc_long="A list of particle files that have been or will be imported into this project",property=None,defaultunits=project_db.get("global.imported_particle_files",dfl=[]),choices=[]))
		return params


	def write_db_entry(self,k,v):
		pass
		

class E2BoxerTask(WorkFlowTask):	
	documentation_string = "Use this tool to review which images in the project have been boxed, and to identify which images need boxing. You can choose from the list of images and hit ok which will lauch e2boxer and automatically load the selected images for boxing"
	
	def __init__(self,application):
		WorkFlowTask.__init__(self,application)
		self.window_title = "e2boxer management"
		self.boxer_module = None # this will actually point to an EMBoxerModule, potentially

	def get_params(self):
		params = []
		project_db = db_open_dict("bdb:project")
		params.append(ParamDef(name="blurb",vartype="text",desc_short="Using e2boxer",desc_long="",property=None,defaultunits=E2BoxerTask.documentation_string,choices=None))
		
		project_file_names = project_db.get("global.micrograph_ccd_filenames",dfl=[])
		num_boxes = self.get_num_particles(project_file_names)
		pnames = ParamDef(name="global.micrograph_ccd_filenames",vartype="stringlist",desc_short="File names",desc_long="The raw data from which particles will be extracted and ultimately refined to produce a reconstruction",property=None,defaultunits=None,choices=project_file_names)
		pboxes = ParamDef(name="Num boxes",vartype="intlist",desc_short="Number of boxes",desc_long="The number of box images stored for this image in the database",property=None,defaultunits=None,choices=num_boxes)
		
		
		p = ParamTable(name="filenames",desc_short="Choose images to box",desc_long="")
		p.append(pnames)
		p.append(pboxes)
		params.append(p)
		
		boxer_project_db = db_open_dict("bdb:e2boxer.project")
		params.append(ParamDef(name="boxsize",vartype="int",desc_short="Box size",desc_long="An integer value",property=None,defaultunits=boxer_project_db.get("working_boxsize",dfl=128),choices=[]))
		params.append(ParamDef(name="method",vartype="choice",desc_short="Boxing mode",desc_long="Currently only one mode is supported, but this could change",property=None,defaultunits="Swarm",choices=["Swarm"]))
		params.append(ParamDef(name="running_mode",vartype="choice",desc_short="Boxing mode",desc_long="Whether to load the GUI or run automatic boxing based on information stored in the database",property=None,defaultunits="gui",choices=["gui","auto_db"]))
		#params.append(ParamDef(name="global.imported_micrograph_ccd_files",vartype="url",desc_short="File names",desc_long="The raw data from which particles will be extracted and ultimately refined to produce a reconstruction",property=None,defaultunits=project_db.get("global.imported_micrograph_ccd_files",dfl=[]),choices=[]))
		return params
	
	def get_num_particles(self,project_files_names):
		if not os.path.exists("particles/EMAN2DB"):
			vals = [ 0 for i in range(len(project_files_names))]
			return vals
		else:
			#particle_db = project_db = db_open_dict("bdb:particles#")
			dbs = db_list_dicts("bdb:particles")
			ptcl_dbs = []
			for db in dbs:
				db_strip = get_file_tag(db)
				if len(db_strip) > 5:
					if db_strip[-6:] == "_ptcls":
						ptcl_dbs.append(db_strip[:-6])
			
			
			vals = []
			for name in project_files_names:
				name_strip = get_file_tag(name)
				if name_strip in ptcl_dbs:
					try:
						db = db_open_dict("bdb:particles#"+name_strip+"_ptcls")
						vals.append(db["maxrec"]+1)
					except:
						vals.append(0)
				else:
					vals.append(0)
			return vals
			
	def on_form_ok(self,params):
		for k,v in params.items():
			if k != "blurb": self.write_db_entry(k,v)
			
		if len(params["filenames"]) == 0:
			self.emit(QtCore.SIGNAL("task_idle"))
			self.application.close_specific(self.form)
			self.form = None
			return

		else:
		
			options = EmptyObject()
			for key in params.keys():
				setattr(options,key,params[key])
				
			self.boxer_module = EMBoxerModule(self.application,options)
			QtCore.QObject.connect(self.boxer_module, QtCore.SIGNAL("e2boxer_idle"), self.on_boxer_idle)
			self.application.close_specific(self.form)
			self.form = None
			
	def on_form_close(self):
		# this is to avoid a task_idle signal, which would be incorrect if e2boxer is running
		if self.boxer_module == None:
			self.emit(QtCore.SIGNAL("task_idle"))
		else: pass
	
	def on_boxer_idle(self):
		self.boxer_module = None
		print "done"
		self.emit(QtCore.SIGNAL("task_idle"))

	def write_db_entry(self,key,value):
		if key == "boxsize":
			boxer_project_db = db_open_dict("bdb:e2boxer.project")
			boxer_project_db["working_boxsize"] = value

		else:
			pass	

class ParticlesTask(WorkFlowTask):
	
	documentation_string = "This is a list of particles images that are associated with this project. You can add and remove file names by editing the text entries directly and or by using the browse and clear buttons."
	
	
	def __init__(self,application):
		WorkFlowTask.__init__(self,application)
		self.window_title = "Project particles"

	def get_params(self):
		params = []
		
		
		particle_file_names = self.get_project_particle_file_names()
		num_boxes = self.get_num_particles(particle_file_names)
		dimensions = self.get_particle_dims(particle_file_names)
		
		
		pnames = ParamDef(name="global.micrograph_ccd_filenames",vartype="stringlist",desc_short="File names",desc_long="The raw data from which particles will be extracted and ultimately refined to produce a reconstruction",property=None,defaultunits=None,choices=particle_file_names)
		pboxes = ParamDef(name="Num boxes",vartype="intlist",desc_short="Number of boxes",desc_long="The number of box images stored for this image in the database",property=None,defaultunits=None,choices=num_boxes)
		pdims = ParamDef(name="Dimensions",vartype="stringlist",desc_short="Dimensions",desc_long="The dimensions of the particle images",property=None,defaultunits=None,choices=dimensions)
		
		
		p = ParamTable(name="filenames",desc_short="Choose images to box",desc_long="")
		p.append(pnames)
		p.append(pboxes)
		p.append(pdims)
		params.append(p)
		
		#boxer_project_db = db_open_dict("bdb:e2boxer.project")
		#params.append(ParamDef(name="boxsize",vartype="int",desc_short="Box size",desc_long="An integer value",property=None,defaultunits=boxer_project_db.get("working_boxsize",dfl=128),choices=[]))
		#params.append(ParamDef(name="global.imported_micrograph_ccd_files",vartype="url",desc_short="File names",desc_long="The raw data from which particles will be extracted and ultimately refined to produce a reconstruction",property=None,defaultunits=project_db.get("global.imported_micrograph_ccd_files",dfl=[]),choices=[]))
		return params
	
	def get_project_particle_file_names(self):
		project_db = db_open_dict("bdb:project")
		project_file_names = project_db.get("global.micrograph_ccd_filenames",dfl=[])
		result = [get_file_tag(name)+"_ptcls" for name in project_file_names]
		
		dbs = db_list_dicts("bdb:particles")
		
		for db in dbs:
			db_strip = get_file_tag(db)
			if db_strip not in result:
				result.append(db)
				
		return result
					
	def get_particle_dims(self,particle_file_names):
		if not os.path.exists("particles/EMAN2DB"):
			vals = [ "" for i in range(len(particle_file_names))]
			return vals
		else:
			vals = []
			for name in particle_file_names:
				try:
					db = db_open_dict("bdb:particles#"+name)
					hdr = db.get_header(0)
					#print hdr
					vals.append(str(hdr["nx"])+'x'+str(hdr["ny"])+'x'+str(hdr["nz"]))
				except:	vals.append("")
			return vals
	
	def get_num_particles(self,particle_file_names):
		if not os.path.exists("particles/EMAN2DB"):
			vals = [ 0 for i in range(len(particle_file_names))]
			return vals
		else:
			vals = []
			for name in particle_file_names:
				try:
					db = db_open_dict("bdb:particles#"+name)
					vals.append(db["maxrec"]+1)
				except:	vals.append(0)
			return vals

	def write_db_entry(self,key,value):
		pass


class MicrographCCDImportTask(WorkFlowTask):	
	documentation_string = "Use this tool for importing flat files into the project database under the raw_data directory. Files that you import in this way will be automatically added the list of files in the project."
	
	def __init__(self,application):
		WorkFlowTask.__init__(self,application)
		self.window_title = "Import micrographs"

	def get_params(self):
		params = []
		project_db = db_open_dict("bdb:project")
		params.append(ParamDef(name="blurb",vartype="text",desc_short="Importing image data",desc_long="",property=None,defaultunits=MicrographCCDImportTask.documentation_string,choices=None))
		params.append(ParamDef(name="import_micrograph_ccd_files",vartype="url",desc_short="File names",desc_long="The raw data from which particles will be extracted and ultimately refined to produce a reconstruction",property=None,defaultunits=[],choices=[]))
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
		print "cancelled"
		
	
class MicrographCCDTask(WorkFlowTask):
	
	documentation_string = "This is a list of micrographs or CCD frames that are associated with this project. You can add and remove file names by editing the text entries directly and or by using the browse and clear buttons. In addition to being able to specify images that are stored on your hard drive in the usual way, you can also choose images from EMAN2 style databases."
	
	
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
		params.append(papix)
		params.append(pvolt)
		params.append(pcs)
		
		return params

	def write_db_entry(self,key,value):
	
		if key == "global.apix":
			db = db_open_dict("bdb:project")
			db["global.apix"] = value
		elif key == "global.microscope_voltage":
			db = db_open_dict("bdb:project")
			db["global.microscope_voltage"] = value
		elif key == "global.microscope_cs":
			db = db_open_dict("bdb:project")
			db["global.microscope_cs"] = value
		else:
			print "unknown key:",key,"this object is",self
		
		
if __name__ == '__main__':
	
	from emapplication import EMStandAloneApplication
	em_app = EMStandAloneApplication()
	sprinit = SPRInitTask(em_app)
	window = sprinit.run_form() 
	
	#em_app.show()
	em_app.execute()	