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

from emform import EMFormModule
from emdatastorage import ParamDef
from PyQt4 import QtGui,QtCore
from EMAN2db import db_check_dict, db_open_dict,db_remove_dict
from EMAN2 import EMData,remove_directories_from_name
import os
import copy


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

	
	
class MicrographCCDImportTask(WorkFlowTask):	
	documentation_string = "Use this tool to manage the data that you have copied (or will copy) into into the project database. If you have done this than you do not necessarily have to keep a copy of the data in the local directory (or elsewhere). The files in this list are automatically included in the project. You can remove entries from list and they will be removed from the project database."
	
	def __init__(self,application):
		WorkFlowTask.__init__(self,application)
		self.window_title = "Import micrographs"

	def get_params(self):
		params = []
		project_db = db_open_dict("bdb:project")
		params.append(ParamDef(name="blurb",vartype="text",desc_short="Importing image data",desc_long="",property=None,defaultunits=MicrographCCDImportTask.documentation_string,choices=None))
		params.append(ParamDef(name="global.imported_micrograph_ccd_files",vartype="url",desc_short="File names",desc_long="The raw data from which particles will be extracted and ultimately refined to produce a reconstruction",property=None,defaultunits=project_db.get("global.imported_micrograph_ccd_files",dfl=[]),choices=[]))
		return params

	def write_db_entry(self,key,value):
		if key == "global.imported_micrograph_ccd_files":
			project_db = db_open_dict("bdb:project")
			remove_files = project_db.get("global.imported_micrograph_ccd_files",dfl=[])
			orig_names = copy.copy(remove_files)
			orig_bdb_names = ["bdb:raw_data#"+remove_directories_from_name(name) for name in orig_names]
			
			no_dir_names = [remove_directories_from_name(name) for name in value]
			
			for name in no_dir_names:
				if no_dir_names.count(name) > 1:
					print "you can't use images with the same name (",name,")"
					return
				
			
			new_names = []
			
			for name in value:
				if os.path.exists(name): # if the name exists 
					new_names.append(name) # store it as a name corresponding to an image that should be stored
					
					for i in range(len(remove_files)-1,-1,-1):
						if remove_files[i] == name:
							remove_files.pop(i) # we don't need to remove it
							cont = False
			
			# now remove the files that the user no longer wants to exist
			for rm in remove_files:
				db_remove_dict("bdb:raw_data#"+remove_directories_from_name(rm))
			
			# now add the files to db (if they don't already exist
			progress = QtGui.QProgressDialog("Importing files into database...", "Abort Copy", 0, len(new_names),None)
			progress.show()
			for i,name in enumerate(new_names):
				progress.setValue(i)
				db_name = "bdb:raw_data#"+remove_directories_from_name(name)
				if name in orig_names: # this means the entry
					if not db_check_dict(db_name):
						print "there is an internal error. The database entry for",name,"doesn't exist"
					else:
						#it's the same image, continue
						continue
				else:
					e = EMData(name)
					e.set_attr("disk_file_name",name)
					e.write_image(db_name,0)
					raw_data_db = db_open_dict(db_name)
					
			new_bdb_names = ["bdb:raw_data#"+remove_directories_from_name(name) for name in new_names]
			
			project_files = project_db.get("global.micrograph_ccd_filenames",dfl=[])
			new_project_files = [ file for file in project_files if file not in orig_bdb_names]
			new_project_files.extend(new_bdb_names)
			project_db["global.micrograph_ccd_filenames"] = new_project_files
			project_db["global.imported_micrograph_ccd_files"] = new_names
			
		else:
			print "unknown key:",key,"this object is",self
			
						
	
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
			project_db = db_open_dict("bdb:project")
			current_project_files = project_db.get("global.micrograph_ccd_filenames",dfl=[])
			
			imported_names = project_db.get("global.imported_micrograph_ccd_files",dfl=[])
			imported_db_names = ["bdb:raw_data#"+remove_directories_from_name(name) for name in imported_names]
			
			
			# Make sure the list of imported names is update to date, and any imported data
			# that has been removed is removed from the database
			keep_imported_db_names = []
			rm_imported_db_names = []
			for name in value:
				if name in imported_db_names:
					file_name = db_open_dict(name).get_header(0)["disk_file_name"]
					keep_imported_db_names.append(file_name)
				else:
					# pass
					db_remove_dict(name) # or we could just pass and leave the data there
		
			project_db["global.imported_micrograph_ccd_files"] = keep_imported_db_names
			
			
			new_names = []
	
			# now just update the static list of project file names
			if value != None:
				e = EMData()
				for i,name in enumerate(value):
					#if len(name) >= 13 and name[:13] == "bdb:raw_data#":
						
					if os.path.exists(name) or db_check_dict(name):
						cont = True
						for i in range(len(current_project_files)-1,-1,-1):
							if current_project_files[i] == name:
								new_names.append(current_project_files.pop(i))
								cont = False
								
						if not cont: continue # the image is already stored in the database
						
						cool_to_go = True
						read_header_only = True
						try:
							a = e.read_image(name,0,read_header_only)
						except: cool_to_go = False
						
						new_names.append(name)
				
				db = db_open_dict("bdb:project")
				db["global.micrograph_ccd_filenames"] = new_names
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