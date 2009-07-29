#!/usr/bin/env python

# Author: Muthu Alagappan, 07/22/09
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

from EMAN2 import *
from PyQt4 import QtCore
from e2foldhunterstat import *
from emapplication import EMStandAloneApplication,get_application
from empdbvaltool import EMPDBValTool
from emplot3d import *

class E2ValidateMed():
	
	def __init__(self):
		self.em_val = None
		self.plot3d =  None
		self.fh_stat = E2FoldHunterStat()

	def start(self):
		if self.em_val == None:
			self.__init_em_val()
				
		get_application().show_specific(self.em_val)

	def __init_em_val(self):
		if self.em_val == None: 
			self.em_val = EMPDBValTool()
			QtCore.QObject.connect(self.em_val.emitter(), QtCore.SIGNAL("run_validate"),self.on_em_validate_requested)
			QtCore.QObject.connect(self.em_val.emitter(), QtCore.SIGNAL("module_closed"),self.on_em_val_closed)
	
	def on_em_val_closed(self):
		self.em_val = None

	def __init_plot3d(self):
		if self.plot3d == None: 
			self.plot3d = EMPlot3DModule()
			QtCore.QObject.connect(self.plot3d.emitter(), QtCore.SIGNAL("view_transform"),self.on_transform_requested)
			QtCore.QObject.connect(self.plot3d.emitter(), QtCore.SIGNAL("module_closed"),self.on_plot3d_closed)

	def on_transform_requested(self, new_pdb_file):
		if self.em_val == None:
			self.__init_em_val()
			get_application().show_specific(self.em_val)
		self.em_val.set_pdb_file(str(new_pdb_file), shouldDelete=True)

	def on_plot3d_closed(self):
		self.plot3d = None

	def on_em_validate_requested(self, mrc_file, pdb_file, trans, iso_thresh):

		vals = {}
		rotList = []
		data = []
		initPoint = []

		vals, rotList, b, data, initPoint = self.fh_stat.gen_data(mrc_file, pdb_file, trans, iso_thresh)
		get_application().close_specific(self.plot3d)
		self.plot3d = None

		if self.plot3d == None:
			self.__init_plot3d()

		self.plot3d.set_Vals(vals)
		self.plot3d.set_Rotations(rotList)
		self.plot3d.set_Probe(b)
		self.plot3d.set_data(data,"Rotation Angles with Z Score")
		self.plot3d.set_data(initPoint, "Original Probe", shape = "Cube")
		get_application().show_specific(self.plot3d)

if __name__ == '__main__':
	from emapplication import EMStandAloneApplication
	em_app = EMStandAloneApplication()
	window = E2ValidateMed()
	window.start()
	#em_app.show()
	em_app.execute()

		
		
