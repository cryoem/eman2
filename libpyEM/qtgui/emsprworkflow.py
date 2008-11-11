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

def db_entry(entry,db_path,alternate):
	# am only doing it this way because there is no current strategy for doing it more generically
	if db_check_dict(db_path):
		db = db_open_dict(db_path)
		ret = db[entry]
		if ret != None: return ret
		else: return alternate
	else: return alternate
	
def write_db_entry(key,value):
	if key == "global.project_files":
		pre_existing_files = db_entry("global.project_files","bdb:project",[])
		new_names = []
		if value != None:
			e = EMData()
			for name in value:
				if os.path.exists(name):
					cont = True
					for i in range(len(pre_existing_files)-1,-1,-1):
						if pre_existing_files[i] == name:
							new_names.append(pre_existing_files.pop(i))
							cont = False
							
					if not cont: continue # the image is already stored in the database
					
					cool_to_go = True
					read_header_only = True
					try:
						a = e.read_image(name,0,read_header_only)
					except: cool_to_go = False
					
					if cool_to_go:
						db_path = "bdb:raw_data#"+remove_directories_from_name(name)
						b = EMData(name)
						b.write_image(db_path,0)
					else: continue
					
					
					new_names.append(name)
			
			print "new names ", new_names
			db = db_open_dict("bdb:project")
			db["global.project_files"] = new_names
			
			print "pre",pre_existing_files
			for pre in pre_existing_files:
				db_remove_dict("bdb:raw_data#"+remove_directories_from_name(pre))
				
	elif key == "global.apix":
		db = db_open_dict("bdb:project")
		db["global.apix"] = value
	elif key == "global.microscope_voltage":
		db = db_open_dict("bdb:project")
		db["global.microscope_voltage"] = value
	elif key == "global.microscope_cs":
		db = db_open_dict("bdb:project")
		db["global.microscope_cs"] = value

	
		
	else:
		pass
	
	
class SPRInitModule:
	'''
	A class that manages the initialization component of a Single Particle
	Reconstruction workflow
	'''
	
	def __init__(self,application):
		self.form = None # will potentially reference an EMFormModule
		self.application = application
		pass
	
	def run_form(self):
		self.form = EMFormModule(self.get_params(),self.application)
		self.application.show_specific(self.form)
		QtCore.QObject.connect(self.form.widget,QtCore.SIGNAL("emform_ok"),self.on_form_ok)
		QtCore.QObject.connect(self.form.widget,QtCore.SIGNAL("emform_cancel"),self.on_form_cancel)
	
	def on_form_ok(self,params):
		for k,v in params.items():
			write_db_entry(k,v)
			
		self.application.close_specific(self.form)
		
	def on_form_cancel(self):
		self.application.close_specific(self.form)
	
	def get_params(self):
		params = []
		params.append(ParamDef(name="global.project_files",vartype="url",desc_short="File names",desc_long="The raw data from which particles will be extracted and ultimately refined to produce a reconstruction",property=None,defaultunits=db_entry("global.project_files","bdb:project",[]),choices=[]))
		papix = ParamDef(name="global.apix",vartype="float",desc_short="A/pix for project",desc_long="The physical distance represented by the pixel spacing",property=None,defaultunits=db_entry("global.apix","bdb:project",1.0),choices=None)
		pvolt = ParamDef(name="global.microscope_voltage",vartype="float",desc_short="Microscope voltage",desc_long="The operating voltage of the microscope",property=None,defaultunits=db_entry("global.microscope_voltage","bdb:project",300.0),choices=None)
		pcs = ParamDef(name="global.microscope_cs",vartype="float",desc_short="Microscope Cs",desc_long="Microscope spherical aberration constant",property=None,defaultunits=db_entry("global.microscope_cs","bdb:project",2.0),choices=None)
		params.append([papix,pvolt,pcs])
		
		return params
		
if __name__ == '__main__':
	
	from emapplication import EMStandAloneApplication
	em_app = EMStandAloneApplication()
	sprinit = SPRInitModule(em_app)
	window = sprinit.run_form() 
	
	#em_app.show()
	em_app.execute()	