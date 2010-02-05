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
	
	
	
	logid=E2init(sys.argv)
	
	em_app = EMStandAloneApplication()
	window = EMSimmxExplorer(application=em_app)
	if len(args) > 0: window.set_projection_file(args[0])
	if len(args) > 1: window.set_particle_file(args[1])

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

		self.ptcl_display = EMImageMx				# display all particles
		self.mx_display = None 					# mx display module for displaying projection and aligned particle
		self.lay = None 						# 2d plot for displaying comparison between particle and projection
		
	def set_data(self,projections,particles):
		'''
		Initialize data
		'''
		
		
		# Deal with particles
		self.ptcl_data=EMData.read_images(particles,range(200))
		self.ptcl_display = EMImageMXModule()
		self.ptcl_display.set_data(self.ptcl_data)
		QtCore.QObject.connect(self.mx_display.emitter(),QtCore.SIGNAL("mx_image_selected"),self.ptcl_selected)		
		QtCore.QObject.connect(self.mx_display.emitter(),QtCore.SIGNAL("module_closed"),self.on_mx_display_closed)
			
	def get_num_particles(self):
		if self.ptcl_data==None : return 0
		return len(self.ptcl_data)
		
	def init_vitals(self):
		'''
		Vital information
		'''
		if not file_exists(self.projection_file): raise RuntimeError("%s does not exist, this is vital" %self.projection_file)
		if not file_exists(self.particle_file): raise RuntimeError("%s does not exist, this is vital" %self.particle_file)

		# deal with projections
		self.proj_data=EMData.read_images(projections)


		self.num_particles = ny
		
		eulers = get_eulers_from(self.projection_file)
		
		self.specify_eulers(eulers)
		
		self.projections = EMData().read_images(self.projection_file)
		
		self.__update_cmp_attribute(False)
		
		self.set_emdata_list_as_data(self.projections,"cmp")
		
	def __update_cmp_attribute(self,update=True):
		'''
		Uses self.current_particle to update the cmp attribute in the projections
		'''
		e = EMData()
		for i in range(len(self.projections)):
			r = Region(i,self.current_particle,1,1)
			e.read_image(self.simmx_file,0,False,r)
			self.projections[i].set_attr("cmp",-e.get(0))		# We plot -1* values so taller peaks are better...
		
		self.regen_dl()
		
	def render(self):
		if self.inspector == None: self.get_inspector()
		if self.projections == None:
			try:self.init_vitals()
			except: return
			#except:
			#	self.get_inspector()
			#	print "failed to initialize vital information"
			#	return
			
		EM3DSymViewerModule.render(self)
	
	def object_picked(self,object_number):
		if object_number == self.current_projection: return
		self.current_projection = object_number
		resize_necessary = False
		if self.mx_display == None:
			self.mx_display = EMImageMXModule()
			QtCore.QObject.connect(self.mx_display.emitter(),QtCore.SIGNAL("module_closed"),self.on_mx_display_closed)
			resize_necessary = True

		if self.frc_display == None:
			self.frc_display = EMPlot2DModule()
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
	
	def set_ptcl_idx(self,idx):
		if self.current_particle != idx:
			self.current_particle = idx
			self.__update_cmp_attribute()
			self.update_display(True)
		
class EMSimmxXplorInspector(EMSymInspector):
	def get_desktop_hint(self):
		return "inspector"
	
	def __init__(self,target,enable_trace=False,enable_og=False) :
		self.ptcl_slider = None
		self.combo = None
		
		EMSymInspector.__init__(self,target,enable_trace=enable_trace,enable_og=enable_og)
		if self.target().projection_file == None:
			self.add_simmx_dir_data_widgets()
		else: self.add_simmx_options()
		
#	def __del__(self):
#		print "simmx xplor died"
		
	def add_simmx_options(self):
		self.simmx_tab= QtGui.QWidget()
		vbl = QtGui.QVBoxLayout(self.simmx_tab)
		
		self.__init_ptcl_slider(vbl)	
		self.tabwidget.insertTab(0,self.simmx_tab,"Simmx")
		self.tabwidget.setCurrentIndex(0)
		
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
		self.connect(self.ptcl_slider, QtCore.SIGNAL("valueChanged(int)"), self.set_ptcl_idx)
		
		
	def set_ptcl_idx(self,val):
		self.target().set_ptcl_idx(val)
		self.score_option_changed() # this might be a hack - it forces the cylinder heights to update.
	
	def add_simmx_dir_data_widgets(self):
		self.data = simmx_xplore_dir_data()
		if len(self.data) == 0: raise RuntimeError("There is no simmx refinement data in the current directory")
		
		self.simmx_dir_tab= QtGui.QWidget()
		vbl = QtGui.QVBoxLayout(self.simmx_dir_tab)
		
		# This is the combo-box with the list of refine_* directories
		combo_entries = [d[0] for d in self.data]
		combo_entries.sort()
		combo_entries.reverse()
		self.combo = QtGui.QComboBox(self)
		for e in combo_entries: self.combo.addItem(e)
			
		self.connect(self.combo,QtCore.SIGNAL("currentIndexChanged(QString&)"),self.on_combo_change)
		self.connect(self.combo,QtCore.SIGNAL("currentIndexChanged(const QString&)"),self.on_combo_change)
			
		vbl.addWidget(self.combo)
		
		self.list_widget = QtGui.QListWidget(None)
		
		self.list_widget.setSelectionMode(QtGui.QAbstractItemView.SingleSelection)
		self.list_widget.setMouseTracking(True)
		QtCore.QObject.connect(self.list_widget,QtCore.SIGNAL("itemClicked(QListWidgetItem *)"),self.list_widget_item_clicked)
	
		self.update_simmx_list(True)
		vbl.addWidget(self.list_widget)
		self.tabwidget.insertTab(0,self.simmx_dir_tab,"Simmx")
		self.tabwidget.setCurrentIndex(0)
		
		self.update_target()
		
		if self.ptcl_slider == None:
			self.__init_ptcl_slider(vbl)
		else:
			dir = str(self.combo.currentText())
			data = (d for d in self.data if d[0] == dir).next()
			n = EMUtil.get_image_count(data[1])
			if n == None: n = 0
			else: n -= 1
			self.ptcl_slider.setRange(0,n)
			
		
	def update_target(self,update=False):
		simmx_file = self.simmx_selection()
		if simmx_file == None:
			return
		
		dir = str(self.combo.currentText())
		data = (d for d in self.data if d[0] == dir).next()
		
		self.target().set_particle_file(data[1])
		
		try:
			idx = ( i for i in range(len(data[3])) if data[3][i] == simmx_file ).next()
		except:
			return
		projection_file = data[2][idx]
		
		self.target().set_simmx_file(simmx_file)
		self.target().set_projection_file(projection_file)
		if update: 
			self.target().init_vitals()
			self.target().update_display()
			self.target().updateGL()
		
	def simmx_selection(self):
		selected_items = self.list_widget.selectedItems() # need to preserve the selection
		if len(selected_items) == 0: return None
		else: return str(selected_items[0].text())
		
	def update_simmx_list(self,first_time=False):
		selected_items = self.list_widget.selectedItems() # need to preserve the selection
		
		dir = str(self.combo.currentText())
		
		data = (d for d in self.data if d[0] == dir).next()
#	
#		
		s_text = None
		if len(selected_items) == 1 :
			s_text = str(selected_items[0].text())
			
		self.list_widget.clear()
		for i,vals in enumerate(data[3]):
			choice = vals
			
			a = QtGui.QListWidgetItem(str(choice),self.list_widget)
			if first_time and i == 0:
				self.list_widget.setItemSelected(a,True)
			elif choice == s_text:
				self.list_widget.setItemSelected(a,True)
				
				
	def on_combo_change(self,s):
		self.update_simmx_list()
		self.update_target(True)
		
	def list_widget_item_clicked(self,item):
		self.update_target(True)
		

def simmx_xplore_dir_data():
	'''
	Returns a list containing the names of files that can be used by e2simmxplor
	Format - [refine_dir, particle_file,[projection_00,projection_01..],[simmx_00, simmx_01..]]
	'''
	dirs,files = get_files_and_directories()
	
	dirs.sort()
	for i in range(len(dirs)-1,-1,-1):
		if len(dirs[i]) != 9:
			dirs.pop(i)
		elif dirs[i][:7] != "refine_":
			dirs.pop(i)
		else:
			try: int(dirs[i][7:])
			except: dirs.pop(i)
	
	
	ret = []
	
	for dir in dirs:
		# Get the name of the input 3-D Model
		register_db_name = "bdb:"+dir+"#register"
		if not db_check_dict(register_db_name): continue
		db = db_open_dict(register_db_name,ro=True)
		if not db.has_key("cmd_dict"): continue
		# need to be able to get the input data
		cmd = db["cmd_dict"]
		if cmd==None or not cmd.has_key("input"): continue
		inp = cmd["input"]
		
		dcts=db_list_dicts("bdb:"+dir)		# list all databases
		projs = []
		simxs  = []
		for name in dcts:
			if "simmx" in name :
				num=name.split("_")[1]		# we assume names of the form simmx_12_ext
				name="bdb:%s#%s"%(dir,name)
				namep="bdb:%s#projections_%s"%(dir,num)
				
				if db_check_dict(namep) :
					projs.append(namep)
					simxs.append(name)
				
		if len(projs) > 0:
			projs.reverse()
			simxs.reverse()
			ret.append([dir,inp,projs,simxs])
	
	return ret
	
if __name__ == '__main__':
	main()