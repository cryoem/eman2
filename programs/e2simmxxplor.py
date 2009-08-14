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
	if len(args) > 2:window.set_particle_file(args[2])

	em_app.show()
	em_app.execute()
	
	E2end(logid)


class EMSimmxExplorer(EM3DSymViewerModule):
	def get_desktop_hint(self): return "image"
	
	def __init__(self,application=None,ensure_gl_context=True,application_control=True,projection_file=None,simmx_file=None,particle_file=None):
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
		self.mirror_eulers = True # If True the drawn Eulers are are also rendered on the opposite side of the sphere - see EM3DSymViewerModule.make_sym_dl_lis
		
		
	def set_projection_file(self,file): self.projection_file = file
	def set_simmx_file(self,file): self.simmx_file = file
	def set_particle_file(self,file): self.particle_file = file
	
	def init_vitals(self):
		'''
		Vital information
		'''
		if not file_exists(self.projection_file): raise RuntimeError("%s does not exist, this is vital" %self.projection_file)
		if not file_exists(self.simmx_file): raise RuntimeError("%s does not exist, this is vital" %self.simmx_file)
		if self.particle_file != None and not file_exists(self.particle_file): raise RuntimeError("%s does not exist, this is vital" %self.particle_file)
		
		n = EMUtil.get_image_count(self.projection_file)
		nx,ny =  gimme_image_dimensions2D(self.simmx_file)

		if nx != n: raise RuntimeError("Error: the number of projections %d does not match the x dimension %d of the simmx file" %(n,nxs))
	
		self.num_particles = ny
		
		from e2eulerxplor import get_eulers_from
		eulers = get_eulers_from(self.projection_file)
		
		self.specify_eulers(eulers)
		
		self.projections = EMData().read_images(self.projection_file)
		
		self.current_particle = 0
		self.__update_cmp_attribute(False)
#		e = EMData()
#		for i in range(len(self.projections)):
#			r = Region(self.current_particle,i,1,1)
#			e.read_image(self.simmx_file,0,False,r)
#			self.projections[i].set_attr("cmp",e.get(0))
		
		self.set_emdata_list_as_data(self.projections,"cmp")
		
	def __update_cmp_attribute(self,update=True):
		'''
		Uses self.current_particle to update the cmp attribute in the projections
		'''
		e = EMData()
		for i in range(len(self.projections)):
			r = Region(i,self.current_particle,1,1)
			e.read_image(self.simmx_file,0,False,r)
			self.projections[i].set_attr("cmp",e.get(0))
		
		self.regen_dl()
		
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
		
		self.__update_display(False)

		if resize_necessary:
			get_application().show_specific(self.mx_display)
			self.mx_display.optimally_resize()
		else:
			self.mx_display.updateGL()
			
		if object_number != self.special_euler:
			self.special_euler = object_number
			self.regen_dl()
			
	def __update_display(self,update=True):
		'''
		Uses self.current_particle and self.current_projection to udpate the self.mx_display
		'''
		if self.mx_display == None: return
		
		data = []
		projection = EMData()
		projection.read_image(self.projection_file,self.current_projection)
		
		if self.particle_file != None:
			r = Region(self.current_projection,self.current_particle,1,1)
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
			particle.process_inplace("normalize.toimage.lsq",{"to":projection})
			
			self.mx_display.set_data([projection,particle])
		else:
			self.mx_display.set_data([projection])
		if update: self.mx_display.updateGL()
			
	def on_mx_display_closed(self):
		self.mx_display = None
		
	def get_inspector(self):
		if not self.inspector : 
			self.inspector=EMSimmxXplorInspector(self)
		return self.inspector
	
	def set_ptcl_idx(self,idx):
		if self.current_particle != idx:
			self.current_particle = idx
			self.__update_cmp_attribute()
			self.__update_display(True)
		
class EMSimmxXplorInspector(EMSymInspector):
	def get_desktop_hint(self):
		return "inspector"
	
	def __init__(self,target,enable_trace=False,enable_og=False) :
		EMSymInspector.__init__(self,target,enable_trace=enable_trace,enable_og=enable_og)
		self.add_simmx_options()
		
	def add_simmx_options(self):
		from PyQt4 import QtGui,QtCore
		self.simmx_tab= QtGui.QWidget()
		vbl = QtGui.QVBoxLayout(self.simmx_tab)
		
		from valslider import ValSlider
		self.ptcl_slider = ValSlider(None,(0,self.target().num_particles-1),"Particle:")
		self.ptcl_slider.setValue(self.target().current_particle)
		self.ptcl_slider.setIntonly(True)
		vbl.addWidget(self.ptcl_slider)
		
		self.connect(self.ptcl_slider, QtCore.SIGNAL("valueChanged"), self.set_ptcl_idx)
			
		self.tabwidget.insertTab(0,self.simmx_tab,"Simmx")
		self.tabwidget.setCurrentIndex(0)
		
	def set_ptcl_idx(self,val):
		self.target().set_ptcl_idx(val)
		self.score_option_changed()
	
	
	
if __name__ == '__main__':
	main()