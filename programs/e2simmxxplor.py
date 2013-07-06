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


import os,sys
from PyQt4 import QtGui,QtCore
from valslider import ValSlider

from e2eulerxplor import get_eulers_from
from emimagemx import EMImageMXWidget
from emimage2d import EMImage2DWidget
from emplot2d import EMPlot2DWidget
from EMAN2db import db_convert_path, db_open_dict, db_check_dict, db_list_dicts
from emapplication import EMApp, get_application
from emglobjects import EM3DGLWidget
from emimage3dsym import EM3DSymModel,EMSymInspector, EMSymViewerWidget
from EMAN2 import *

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog  <simmx file> <projection file>  <particles file>

	This program allows you to look at per-particle classification based on a pre-computed similarity
	matrix. If a particle seems to be mis-classified, you can use this program to help figure out
	why. It will display a single asymmetric triangle on the unit sphere, with cylinders representing
	the similarity value for the selected particle vs. each possible reference projection. Use the
	normal middle-click on the asymmetric unit viewer for more options.
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()



	logid=E2init(sys.argv,options.ppid)

	em_app = EMApp()

	window = EMSymViewerWidget()
	explorer = EMSimmxExplorer(window)
	window.model = explorer

	if len(args) > 0: explorer.set_simmx_file(args[0])
	if len(args) > 1: explorer.set_projection_file(args[1])
	if len(args) > 2: explorer.set_particle_file(args[2])

	em_app.show()
	em_app.execute()

	E2end(logid)


class EMSimmxExplorer(EM3DSymModel):
	def __init__(self, gl_widget, projection_file=None,simmx_file=None,particle_file=None):
		self.init_lock = True # a lock indicated that we are still in the __init__ function
		self.au_data = None # This will be a dictionary, keys will be refinement directories, values will be something like available iterations for visual study
		EM3DSymModel.__init__(self, gl_widget)
		self.window_title = "SimmxXplor"

		self.projection_file = projection_file # a projection file produced by e2project3d
		self.simmx_file = simmx_file # a similarity matrix produced by e2simmx or by e2classaverage - should have the same number of columns (nx) as the number of projections
		self.particle_file = particle_file # a similarity matrix produced by e2simmx or by e2classaverage - should have the same number of columns (nx) as the number of projections
		self.num_particles = None # keep track of the total number of particles ( the y dimension of the simmx file)
		self.current_particle = 0 # keep track of the current particle
		self.current_projection = None # keep track of the current projection
		self.projections = None # a list of the projections - the initial design keeps them all in memory - this could be overcome
		self.mx_display = None # mx display widget for displaying projection and aligned particle
		self.frc_display = None # 2d plot for displaying comparison between particle and projection
		self.mirror_eulers = True # If True the drawn Eulers are are also rendered on the opposite side of the sphere - see EM3DSymModel.make_sym_dl_lis

	def set_data(self,simmx_file):
		'''
		Can only set data using the name of a simmx file
		'''
		if "EMAN2DB" in simmx_file: simmx_file = db_convert_path(simmx_file)
		self.simmx_file = simmx_file
		e = EMData()
		e.read_image(simmx_file,0,True)
		self.projection_file = e["projection_file"]
		self.particle_file = e["particle_file"]
		self.init_vitals()
		self.regen_dl()

	def set_projection_file(self,file): self.projection_file = file
	def set_simmx_file(self,file): self.simmx_file = file
	def set_particle_file(self,file): self.particle_file = file

	def get_num_particles(self):
		if self.simmx_file == None: return None
		nx,ny =  gimme_image_dimensions2D(self.simmx_file)
		return ny
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
		allcmp=[]
		for i in range(len(self.projections)):
			r = Region(i,self.current_particle,1,1)
			e.read_image(self.simmx_file,0,False,r)
			allcmp.append(e.get(0))

		if min(allcmp)>0 :
			for i in range(len(self.projections)):
				self.projections[i].set_attr("cmp",1.0/allcmp[i])		# We plot reciprocalk values so taller peaks are better...
		else :
			for i in range(len(self.projections)):
				self.projections[i].set_attr("cmp",-allcmp[i])		# We plot -1* values so taller peaks are better...

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

		EM3DSymModel.render(self)

	def object_picked(self,object_number):
		if object_number == self.current_projection: return
		self.current_projection = object_number
		resize_necessary = False
		if self.mx_display == None:
			self.mx_display = EMImage2DWidget()
#			self.mx_display = EMImageMXWidget()
			QtCore.QObject.connect(self.mx_display,QtCore.SIGNAL("module_closed"),self.on_mx_display_closed)
			resize_necessary = True

		if self.frc_display == None:
			self.frc_display = EMPlot2DWidget()
#			QtCore.QObject.connect(self.frc_display,QtCore.SIGNAL("module_closed"),self.on_frc_display_closed)

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
			particle["xform.align2d"]=t
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

					awfrcm=[frcm[nf+i]*snr[i] for i in range(len(snr))]
					self.frc_display.set_data((ses,awfrcm),"frcm*SNR",color=3)

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
	dirs=[i for i in dirs if "refine_" in i]

	ret = []

	for dir in dirs:
			dl=os.listdir(dir)
			ptcl=js_open_dict("{path}/0_refine_parms.json".format(path=dir))["input"]
			pfx=dir+"/"

			prjs=[pfx+i for i in dl if "projections_" in i and i[-4:]==".hdf" and "_even" in i]
			simmx=[pfx+i for i in dl if "simmx_" == i[:6] and "stg1_" not in i and "_even" in i]

			prjs+=[pfx+i for i in dl if "proj_stg1" in i and i[-4:]==".hdf" and "_even" in i]
			simmx+=[pfx+i for i in dl if "simmx_stg1" == i[:10] and "_even" in i]

			if len(prjs)==len(simmx) : ret.append([dir,str(ptcl[0]),prjs,simmx])
			else : print "Mismatch in :",zip(prjs,simmx)

			prjs=[pfx+i for i in dl if "projections_" in i and i[-4:]==".hdf" and "_odd" in i]
			simmx=[pfx+i for i in dl if "simmx_" == i[:6] and "stg1_" not in i and "_odd" in i]

			prjs+=[pfx+i for i in dl if "proj_stg1" in i and i[-4:]==".hdf" and "_odd" in i]
			simmx+=[pfx+i for i in dl if "simmx_stg1" == i[:10] and "_odd" in i]

			if len(prjs)==len(simmx) : ret.append([dir,str(ptcl[1]),prjs,simmx])
			else : print "Mismatch in :",zip(prjs,simmx)

	print ret

	return ret

if __name__ == '__main__':
	main()
