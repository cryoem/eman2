#!/usr/bin/env python

#
# Author: Steve Ludtke 02/05/10
# Copyright (c) 2010 Baylor College of Medicine
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
from PyQt4 import QtGui,QtCore
from e2eulerxplor import get_eulers_from
from emimagemx import EMImageMXModule
from emplot2d import EMPlot2DModule
from valslider import ValSlider

	
def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog  <projection file>  <particles file>
	
Compare different similarity metrics for a given particle stack. Note that all particles and projections are
read into memory. Do not use it on large sets of particles !!!
"""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	(options, args) = parser.parse_args()
	
	if len(args)<2 :
		print "Error, please specify projection file and particles file"
		sys.exit(1)
	
	logid=E2init(sys.argv)
	
	em_app = EMStandAloneApplication()
	window = EMCmpExplorer(application=em_app)
	window.set_data(args[0],args[1])

	em_app.show()
	em_app.execute()
	
	E2end(logid)


class EMCmpExplorer(EM3DSymViewerModule):
	def get_desktop_hint(self): return "image"
	
	def __init__(self,application=None,ensure_gl_context=True,application_control=True,projection_file=None,simmx_file=None,particle_file=None):
		self.init_lock = True # a lock indicated that we are still in the __init__ function
		self.au_data = None # This will be a dictionary, keys will be refinement directories, values will be something like available iterations for visual study	
		EM3DSymViewerModule.__init__(self,application,ensure_gl_context=ensure_gl_context,application_control=application_control)
		self.window_title = "SimmxXplor"
		#InputEventsManager.__init__(self)
	
		self.projection_file = projection_file 	# a projection file produced by e2project3d
		self.particle_file = particle_file 		# A file containing particles to be examined

		self.current_particle = 0 				# keep track of the current particle
		self.current_projection = None 			# keep track of the current projection

		self.ptcl_display = None			# display all particles
		self.mx_display = None 					# mx display module for displaying projection and aligned particle
		self.lay = None 						# 2d plot for displaying comparison between particle and projection
		
	def set_data(self,projections,particles):
		'''
		Initialize data
		'''
		if not file_exists(projections): raise RuntimeError("%s does not exist" %self.projection_file)
		if not file_exists(particles): raise RuntimeError("%s does not exist" %self.particle_file)
		
		
		# Deal with particles
		self.ptcl_data=EMData.read_images(particles,range(200))
		if self.ptcl_display==None : 
			self.ptcl_display = EMImageMXModule()
			self.ptcl_display.set_mouse_mode("App")
		self.ptcl_display.set_data(self.ptcl_data)
		QtCore.QObject.connect(self.ptcl_display.emitter(),QtCore.SIGNAL("mx_image_selected"),self.ptcl_selected)		
		QtCore.QObject.connect(self.ptcl_display.emitter(),QtCore.SIGNAL("module_closed"),self.on_mx_display_closed)

		# deal with projections
		self.proj_data=EMData.read_images(projections)

		eulers = [i["xform.projection"] for i in self.proj_data]
		self.specify_eulers(eulers)
		
		for i in self.proj_data : i["cmp"]=0.0
		self.set_emdata_list_as_data(self.proj_data,"cmp")

	def get_num_particles(self):
		if self.ptcl_data==None : return 0
		return len(self.ptcl_data)
		
	def render(self):
		if self.inspector == None: self.get_inspector()
			
		EM3DSymViewerModule.render(self)
	
	def object_picked(self,object_number):
		if object_number == self.current_projection: return
		self.current_projection = object_number
		resize_necessary = False
		if self.mx_display == None:
			self.mx_display = EMImageMXModule()
			QtCore.QObject.connect(self.mx_display.emitter(),QtCore.SIGNAL("module_closed"),self.on_mx_display_closed)
			resize_necessary = True

		#if self.frc_display == None:
			#self.frc_display = EMPlot2DModule()
#			QtCore.QObject.connect(self.frc_display.emitter(),QtCore.SIGNAL("module_closed"),self.on_frc_display_closed)

		self.update_display(False)

		if resize_necessary:
			get_application().show_specific(self.mx_display)
			self.mx_display.optimally_resize()
			get_application().show_specific(self.frc_display)
#			self.frc_display.optimally_resize()
		else:
			self.mx_display.updateGL()
			self.frc_display.updateGL()
			
		if object_number != self.special_euler:
			self.special_euler = object_number
			self.regen_dl()
			
	def update_display(self,update=True):
		'''
		Uses self.current_particle and self.current_projection to udpate the self.mx_display
		'''
		if self.mx_display == None : return
		
		data = []
		projection = EMData()
		projection.read_image(self.projection_file,self.current_projection)
		
		if self.particle_file != None:
			r = Region(self.current_projection,self.current_particle,1,1)
			e = EMData()
			d = {}
			e.read_image(self.simmx_file,0,False,r)
			simval=e.get(0)
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
			#n = EMUtil.get_image_count(self.particle_file)
			particle.read_image(self.particle_file,self.current_particle)
			particle.transform(t)
			particle.process_inplace("normalize.toimage",{"to":projection})
			if particle["norm_mult"]<0 : particle*=-1.0
			
			particle_masked=particle.copy()
			tmp=projection.process("threshold.notzero")
			particle_masked.mult(tmp)
			
			projection["cmp"]=0
			particle["cmp"]=simval
			particle_masked["cmp"]=0
			self.mx_display.set_data([projection,particle,particle_masked])
			
			if self.frc_display != None :
				try : apix=particle["ctf"].apix
				except: apix=particle["apix_x"]
				
				frc=projection.calc_fourier_shell_correlation(particle)
				frcm=projection.calc_fourier_shell_correlation(particle_masked)
				nf=len(frc)/3
				self.frc_display.set_data(([i/apix for i in frc[:nf]],frc[nf:nf*2]),"frc")
				self.frc_display.set_data(([i/apix for i in frcm[:nf]],frcm[nf:nf*2]),"frcm")
				
				try : 
					ctf=particle["ctf"]
					ds=1.0/(ctf.apix*particle["ny"])
					snr=ctf.compute_1d(particle["ny"],ds,Ctf.CtfType.CTF_SNR)
					ses=[i*ds for i in range(particle["ny"]/2)]
					self.frc_display.set_data((ses,snr),"SNR",color=1)

					amp=ctf.compute_1d(particle["ny"],ds,Ctf.CtfType.CTF_AMP)
					self.frc_display.set_data((ses,amp),"CTF",color=2)
				except : ctf=None

		else:
			self.mx_display.set_data([projection])
		if update: self.mx_display.updateGL()
			
	def on_mx_display_closed(self):
		self.mx_display = None
		
	def get_inspector(self):
		if not self.inspector : 
			self.inspector=EMSimmxXplorInspector(self)
		return self.inspector
	
	def ptcl_selected(self,event,lc):
		"""slot for image selection events from image mx"""
		self.set_ptcl_idx(lc[0])
	
	def set_ptcl_idx(self,idx):
		"""Select the index of the current particle to use for comparisons"""
		if self.current_particle != idx:
			self.current_particle = idx

			progress = QtGui.QProgressDialog("Computing alignments", "Abort", 0, len(self.proj_data),None)
			progress.show()
			# redetermines particle alignments
			# then we can quickly compute a series of different similarity values
			ptcl=self.ptcl_data[idx]
			for i,p in enumerate(self.proj_data):
				ali=p.align("rotate_translate_flip",ptcl,{},"dot",{})
				ali=p.align("refine",ptcl,{"xform.align2d":ali["xform.align2d"]},"dot",{})
				p["ptcl.align2d"]=ali["xform.align2d"]
				progress.setValue(i)
				QtCore.QCoreApplication.instance().processEvents()
			
			progress.close()
			self.set_cmp("dot:normalize=1")
			self.update_display(True)
	
	def set_cmp(self,cmpstring):
		"""Select the comparator. Passed as a standard name:attr=value:attr=value string"""
		cmpopt=parsemodopt(str(cmpstring))
		
		progress = QtGui.QProgressDialog("Computing similarities", "Abort", 0, len(self.proj_data),None)
		progress.show()
		for i,p in enumerate(self.proj_data):
			ali=p.copy()
			ali.transform(p["ptcl.align2d"])
			p["cmp"]=-self.ptcl_data[self.current_particle].cmp(cmpopt[0],ali,cmpopt[1])
			progress.setValue(i)
			QtCore.QCoreApplication.instance().processEvents()
			
		progress.close()
		self.regen_dl(True)
		EM3DSymViewerModule.render(self)
	
class EMSimmxXplorInspector(EMSymInspector):
	def get_desktop_hint(self):
		return "inspector"
	
	def __init__(self,target,enable_trace=False,enable_og=False) :		
		EMSymInspector.__init__(self,target,enable_trace=enable_trace,enable_og=enable_og)
		self.add_cmp_options()
		
#	def __del__(self):
#		print "simmx xplor died"
		
	def add_cmp_options(self):
		self.cmp_tab= QtGui.QWidget()
		vbl = QtGui.QVBoxLayout(self.cmp_tab)
		
		self.cmp_combo=QtGui.QComboBox()
		self.cmp_combo.setEditable(True)
		self.cmp_combo.setInsertPolicy(self.cmp_combo.InsertAlphabetically)
		self.cmp_combo.addItem("dot:normalize=1")
		self.cmp_combo.addItem("phase:snrweight=1")
		self.cmp_combo.addItem("phase:snrweight=1:minres=200:maxres=5")
		self.cmp_combo.addItem("phase:snrweight=1:minres=200:maxres=10")
		self.cmp_combo.addItem("phase:snrweight=1:minres=200:maxres=15")
		self.cmp_combo.addItem("phase:snrweight=1:minres=100:maxres=15")
		self.cmp_combo.addItem("phase:snrweight=1:minres=50:maxres=15")
		self.cmp_combo.addItem("frc:snrweight=1")
		self.cmp_combo.addItem("frc:snrweight=1:minres=200:maxres=5")
		self.cmp_combo.addItem("frc:snrweight=1:minres=200:maxres=10")
		self.cmp_combo.addItem("frc:snrweight=1:minres=200:maxres=15")
		vbl.addWidget(self.cmp_combo)
		
		self.tabwidget.insertTab(0,self.cmp_tab,"Cmp")
		self.tabwidget.setCurrentIndex(0)
		
		self.connect(self.cmp_combo, QtCore.SIGNAL("currentIndexChanged(QString)"), self.cmp_changed)
		
	def __init_ptcl_slider(self,layout):
		if self.combo == None: n = self.target().get_num_particles()
		else:
			dir = str(self.combo.currentText())
			data = (d for d in self.data if d[0] == dir).next()
			n = EMUtil.get_image_count(data[1])
		
		if n == None: n = 0
		else: n -= 1
		#self.ptcl_slider = ValSlider(None,(0,n),"Particle:")
		#self.ptcl_slider.setValue(self.target().current_particle)
		#self.ptcl_slider.setIntonly(True)
		#layout.addWidget(self.ptcl_slider)
		#self.connect(self.ptcl_slider, QtCore.SIGNAL("valueChanged"), self.set_ptcl_idx)
		self.ptcl_slider=QtGui.QSpinBox()
		self.ptcl_slider.setRange(0,1000)
		self.ptcl_slider.setSingleStep(1)
		self.ptcl_slider.setValue(0)
		layout.addWidget(self.ptcl_slider)
		self.connect(self.ptcl_slider, QtCore.SIGNAL("valueChanged(int)"), self.on_)
		
		
	def set_ptcl_idx(self,val):
		self.target().set_ptcl_idx(val)
		self.score_option_changed() # this might be a hack - it forces the cylinder heights to update.
	
		
	def simmx_selection(self):
		selected_items = self.list_widget.selectedItems() # need to preserve the selection
		if len(selected_items) == 0: return None
		else: return str(selected_items[0].text())

	def cmp_changed(self,s):
		"When a new comapartor string is selected"
		self.target().set_cmp(s)
#		self.target.regen_dl()
			

	
if __name__ == '__main__':
	main()