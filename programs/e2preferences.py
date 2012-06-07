#!/usr/bin/env python

#
# Author: David Woolford (woolford@bcm.edu)
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
import EMAN2db
from emdatastorage import ParamDef 

from emsprworkflow import WorkFlowTask
from emapplication import get_application

global HOMEDB
HOMEDB=EMAN2db.EMAN2DB.open_db()

class EMPreferencesTask(WorkFlowTask):
	def __init__(self,application=None):
		WorkFlowTask.__init__(self)
		self.preferred_size = [240,240]
		self.window_title = "EMAN2 preferences"
		
	def get_params(self):
		params = []
		params.append(["e2display",self.__get_display_preference_params()])
		params.append(["e2spt_boxer",self.__get_tomoboxer_preference_params()])
		return params
		
	def __get_display_preference_params(self):
		HOMEDB.open_dict("display_preferences")
		db = HOMEDB.display_preferences
		p2d_auto_contrast = ParamDef(name="display_2d_auto_contrast",vartype="boolean",desc_short="2D image auto contrast",desc_long="Should the 2D image display module adjust contrast settings automatically?",property=None,defaultunits=db.get("display_2d_auto_contrast",dfl=True),choices=None)
		p2d_stack_auto_contrast = ParamDef(name="display_stack_auto_contrast",vartype="boolean",desc_short="Stack (2D) - auto contrast",desc_long="Should the stack display module adjust contrast settings automatically?",property=None,defaultunits=db.get("display_stack_auto_contrast",dfl=True),choices=None)
		p2d_stack_n = ParamDef(name="display_stack_np_for_auto",vartype="int",desc_short="Stack (2D) - # particles used for contrast settings",desc_long="When the stack viewer starts up it investigates the parameters of the first n particles to determine contrast settings. Specify -1 to force the stack viewer to investigate all particles.",property=None,defaultunits=db.get("display_stack_np_for_auto",dfl=5),choices=None)

		params = []
		params.append(p2d_auto_contrast)
		params.append(p2d_stack_auto_contrast)
		params.append(p2d_stack_n)
		
		self.__display_entries = ["display_2d_auto_contrast","display_stack_auto_contrast","display_stack_np_for_auto"]
		
		return params
	
	def __get_tomoboxer_preference_params(self):
		HOMEDB.open_dict("e2tomoboxer_preferences")
		db = HOMEDB.e2tomoboxer_preferences
		plargest_dim = ParamDef(name="largest_allowable_dimension",vartype="int",desc_short="Shrink to this",desc_long="The largest permissible image dimension of a tomogram after shrinking",property=None,defaultunits=db.get("largest_allowable_dimension",dfl=1024),choices=None)

		params = []
		params.append(plargest_dim)
		
		self.__tomoboxer_entries = ["largest_allowable_dimension"]
		return params

	def on_form_ok(self, params):
		self.write_db_entries(params)
		self.form.close()
		self.form = None
		from PyQt4 import QtCore
		self.emit(QtCore.SIGNAL("task_idle"))
	
	def run_form(self):
		from emform import EMTableFormWidget
		self.form = EMTableFormWidget(self.get_params())
		self.form.resize(*self.preferred_size)
		self.form.setWindowTitle(self.window_title)
		get_application().show_specific(self.form)
		
		from PyQt4 import QtCore
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_ok"),self.on_form_ok)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_cancel"),self.on_form_cancel)
		QtCore.QObject.connect(self.form,QtCore.SIGNAL("emform_close"),self.on_form_close)
		
	def write_db_entries(self,params):
		HOMEDB.open_dict("display_preferences")
		db_display = HOMEDB.display_preferences
		HOMEDB.open_dict("e2tomoboxer_preferences")
		db_tomo = HOMEDB.e2tomoboxer_preferences
		tmp_db = None
		for key,value in params.items():
			if key in self.__tomoboxer_entries:
				tmp_db = db_tomo
			elif key in self.__display_entries:
				tmp_db = db_display
			else: raise NotImplementedError("A developer has probably forgotten something")
			tmp_db[key] = value
			
		
def main():
	from emapplication import EMApp

	em_app = EMApp()
	
	pref_task = EMPreferencesTask()
	pref_task.run_form()
	em_app.execute()


if __name__ == "__main__":
	main()