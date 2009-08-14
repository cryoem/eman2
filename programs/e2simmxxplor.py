#!/usr/bin/env python

#
# Author: David Woolford 08/13/08 (woolford@bcm.edu
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

from emapplication import EMStandAloneApplication, get_application
from emimage3dsym import EM3DSymViewerModule,EMSymInspector
from e2eulerxplor import InputEventsManager
import os,sys
from EMAN2 import *
from optparse import OptionParser
	
def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog  <simmx file> <projection file>  <particles file>
	
Simmx xplor for comparator evaluation
"""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	(options, args) = parser.parse_args()
	
	
	logid=E2init(sys.argv)
	
	em_app = EMStandAloneApplication()
	window = EMSimmxExplorer(application=em_app)
	window.set_simmx_file(args[0])
	window.set_projection_file(args[1])
	window.set_particle_file(args[2])

	em_app.show()
	em_app.execute()
	
	E2end(logid)


class EMSimmxExplorer(EM3DSymViewerModule):
	def get_desktop_hint(self): return "image"
	
	def __init__(self,application=None,ensure_gl_context=True,application_control=True,projection_file="",simmx_file="",particle_file=""):
		self.init_lock = True # a lock indicated that we are still in the __init__ function
		self.au_data = None # This will be a dictionary, keys will be refinement directories, values will be something like available iterations for visual study	
		EM3DSymViewerModule.__init__(self,application,ensure_gl_context=ensure_gl_context,application_control=application_control)
		#InputEventsManager.__init__(self)
	
		self.projection_file = projection_file # a projection file produced by e2project3d
		self.simmx_file = simmx_file # a similarity matrix produced by e2simmx or by e2classaverage - should have the same number of columns (nx) as the number of projections
		self.particle_file = particle_file # a similarity matrix produced by e2simmx or by e2classaverage - should have the same number of columns (nx) as the number of projections
		self.num_particles = None # keep track of the total number of particles ( the y dimension of the simmx file)
		self.current_particle = None # keep track of the current particle
		self.current_projection = None # keep track of the current projection
		self.projections = None # a list of the projections - the initial design keeps them all in memory - this could be overcome
		self.mx_display = None # mx display module for displaying projection and aligned particle
		
		
	def set_projection_file(self,file): self.projection_file = file
	def set_simmx_file(self,file): self.simmx_file = file
	def set_particle_file(self,file): self.particle_file = file
	
	def init_vitals(self):
		'''
		Vital information
		'''
		if not file_exists(self.projection_file): raise RuntimeError("%s does not exist, this is vital" %self.projection_file)
		if not file_exists(self.simmx_file): raise RuntimeError("%s does not exist, this is vital" %self.simmx_file)
		if not file_exists(self.particle_file): raise RuntimeError("%s does not exist, this is vital" %self.particle_file)
		
		n = EMUtil.get_image_count(self.projection_file)
		nx,ny =  gimme_image_dimensions2D(self.simmx_file)

		if nx != n: raise RuntimeError("Error: the number of projections %d does not match the x dimension %d of the simmx file" %(n,nxs))
	
		self.num_particles = ny
		
		from e2eulerxplor import get_eulers_from
		eulers = get_eulers_from(self.projection_file)
		
		self.specify_eulers(eulers)
		
		self.projections = EMData().read_images(self.projection_file)
		
		self.current_particle = 0
		e = EMData()
		for i in range(len(self.projections)):
			r = Region(self.current_particle,i,1,1)
			e.read_image(self.simmx_file,0,False,r)
			self.projections[i].set_attr("cmp",e.get(0))
		
		self.set_emdata_list_as_data(self.projections,"cmp")
	
	def render(self):
		if self.projections == None:
			self.init_vitals()
			
		EM3DSymViewerModule.render(self)
	
	def object_picked(self,object_number):
		if object_number == self.current_projection: return
		self.current_projection = object_number
		resize_necessary = False
		if self.mx_display == None:
			from emimagemx import EMImageMXModule
			self.mx_display = EMImageMXModule()
			from PyQt4 import QtCore
			QtCore.QObject.connect(self.mx_display.emitter(),QtCore.SIGNAL("module_closed"),self.on_mx_display_closed)
			resize_necessary = True
		
		data = []
		projection = EMData()
		projection.read_image(self.projection_file,self.current_projection)
		
		r = Region(self.current_particle,self.current_projection,1,1)
		e = EMData()
		d = {}
		e.read_image(self.simmx_file,1,False,r)
		d["tx"] = e.get(0)
		e.read_image(self.simmx_file,2,False,r)
		d["ty"] =  e.get(0)
		e.read_image(self.simmx_file,3,False,r)
		d["alpha"] = e.get(0)
		e.read_image(self.simmx_file,4,False,r)
		d["mirror"] = bool(e.get(0))
		d["type"] = "2d"
		t = Transform(d)
		
		particle = EMData()
		n = EMUtil.get_image_count(self.particle_file)
		particle.read_image(self.particle_file,self.current_particle)
		particle.transform(t)
		
		self.mx_display.set_data([projection,particle])
		if resize_necessary:
			get_application().show_specific(self.mx_display)
			self.mx_display.optimally_resize()
		else:
			self.mx_display.updateGL()
			
		if object_number != self.special_euler:
			self.special_euler = object_number
			self.regen_dl()
			
	def on_mx_display_closed(self):
		self.mx_display = None
#		
#class EMSimmxXplorInspector(EMSymInspector):
#	def get_desktop_hint(self):
#		return "inspector"
#	
#	def __init__(self,target,enable_trace=False,enable_og=False) :
#		EMSymInspector.__init__(self,target,enable_trace=enable_trace,enable_og=enable_og)
#		
#	def add_au_table(self):
#		
#		self.au_tab= QtGui.QWidget()
#		self.au_tab.vbl = QtGui.QVBoxLayout(self.au_tab)
#		
#		self.au_data = self.target().au_data
#		combo_entries = self.au_data.keys()
#		combo_entries.sort()
#		combo_entries.reverse()
#		self.combo = QtGui.QComboBox(self)
#		for e in combo_entries:
#			self.combo.addItem(e)
#			
#		self.connect(self.combo,QtCore.SIGNAL("currentIndexChanged(QString&)"),self.on_combo_change)
#		self.connect(self.combo,QtCore.SIGNAL("currentIndexChanged(const QString&)"),self.on_combo_change)
#			
#		self.au_tab.vbl.addWidget(self.combo)
#		self.refine_dir = combo_entries[0]
#		
#		self.list_widget = QtGui.QListWidget(None)
#		
#		self.list_widget.setSelectionMode(QtGui.QAbstractItemView.SingleSelection)
#		self.list_widget.setMouseTracking(True)
#		QtCore.QObject.connect(self.list_widget,QtCore.SIGNAL("itemClicked(QListWidgetItem *)"),self.list_widget_item_clicked)
#	
#		self.update_classes_list(first_time=True)
#		self.au_tab.vbl.addWidget(self.list_widget)
#		self.tabwidget.insertTab(0,self.au_tab,"Refinement")
#		self.tabwidget.setCurrentIndex(0)
#	
#	def on_combo_change(self,s):
#		self.refine_dir = str(s)
#		self.update_classes_list()
#	
#	def update_classes_list(self,first_time=False):
#		selected_items = self.list_widget.selectedItems() # need to preserve the selection
#		
#		s_text = None
#		if len(selected_items) == 1 :
#			s_text = str(selected_items[0].text())
#			if len(s_text) > 4: s_text = s_text[-4:] 
#			
#		self.list_widget.clear()
#		for i,vals in enumerate(self.au_data[self.refine_dir]):
#			choice = vals[0]
#			
#			a = QtGui.QListWidgetItem(str(choice),self.list_widget)
#			if first_time and i == 0:
#				self.list_widget.setItemSelected(a,True)
#			elif len(choice) > 4 and (choice[-4:] == s_text):
#				self.list_widget.setItemSelected(a,True)
#				
#		selected_items = self.list_widget.selectedItems() # need to preserve the selection
#		if len(selected_items) == 1:
#			self.emit(QtCore.SIGNAL("au_selected"),self.refine_dir,str(selected_items[0].text()))
#	
#	def list_widget_item_clicked(self,item):
#		self.emit(QtCore.SIGNAL("au_selected"),self.refine_dir,str(item.text()))
#	
#	def on_mouse_mode_clicked(self,bool):
#		for button in self.mouse_mode_buttons:
#			if button.isChecked():
#				self.target().set_events_mode(str(button.text()))
	
if __name__ == '__main__':
	main()