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

from em3dmodule import *
from EMAN2 import PDBReader

from empdbviewer import *
from emimage3diso import *


class EMPDBValTool(EM3DModule):
	'''
	EMPDB versus isosurface visual evaluation
	'''
	def __init__(self, application=None,ensure_gl_context=True,application_control=True):
		EM3DModule.__init__(self,application,ensure_gl_context=ensure_gl_context,application_control=application_control)
		#self.fName = raw_input ("Enter the file name of a pdb file: ")
		self.pdb_module = None # will eventually be a EMPDBViewer
		self.iso_module = None # will eventuall be a EMImage3DModule
	
		self.t = 0
		
	def __init_pdb_module(self):
		if self.pdb_module == None:
			self.pdb_module = EMPDBViewer(None,False,False)
			self.__set_module_contexts(self.pdb_module)
			self.get_inspector().addTab(self.pdb_module.get_inspector(),"PDB")
			
	def __set_module_contexts(self,module):
		module.set_qt_context_parent(self.qt_context_parent)
		module.set_gl_context_parent(self.gl_context_parent)
		module.set_gl_widget(self.gl_context_parent)
		module.set_dont_delete_parent() # stops a RunTimeError
		module.under_qt_control = self.under_qt_control
		
	def __init_iso_module(self):
		if self.iso_module == None:
			self.iso_module = EMIsosurfaceModule(None,None,False,False,enable_file_browse=True)
			self.get_inspector().addTab(self.iso_module.get_inspector(),"Isosurface")
			self.__set_module_contexts(self.iso_module)
	
	def set_pdb_file(self,pdb_file):
		if self.pdb_module == None:
			self.__init_pdb_module()
		self.pdb_module.set_current_text(pdb_file)
		
	def set_iso_file(self,file_name):
		data = EMData(file_name)
		self.set_iso_data(data)
	
	def set_iso_data(self,data):
		if self.iso_module == None: self. __init_iso_module()
		self.iso_module.set_data(data)
	
	def draw_objects(self):
		
		if self.pdb_module == None:
			self.__init_pdb_module()
		if self.iso_module == None:
			self.__init_iso_module()


		if self.pdb_module != None:
			glPushMatrix()
			self.pdb_module.draw_objects()
			glPopMatrix()
			
		if self.iso_module != None:
			glPushMatrix()
			self.iso_module.render()
			glPopMatrix()
			
	def get_inspector(self):
		if self.inspector == None:
			self.inspector = EMPDBValToolInspector(self)
		return self.inspector
	

class EMPDBValToolInspector(EM3DInspector):
	def __init__(self,target,enable_advanced=False):
		EM3DInspector.__init__(self,target,enable_advanced)
		
	def addTab(self,widget,name):
		self.tabwidget.addTab(widget,name)

	
if __name__ == '__main__':
	from emapplication import EMStandAloneApplication
	em_app = EMStandAloneApplication()
	window = EMPDBValTool()
	em_app.show()
#	window.set_pdb_file("fh-solution-0-1UF2-T.pdb")
#	window.set_iso_file("rdv-target2.mrc")
	em_app.execute()

		
		