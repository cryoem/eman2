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

from emapplication import EMApp, get_application
from emimage3dsym import EM3DSymModel,EMSymInspector
from e2eulerxplor import InputEventsManager
import os,sys
from EMAN2 import *
from PyQt4 import QtGui,QtCore
from e2eulerxplor import get_eulers_from
from emimagemx import EMImageMXModule
from emplot2d import EMPlot2DModule
from valslider import ValSlider

	
def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog  <projection file>  <particles file>
	
Compare different similarity metrics for a given particle stack. Note that all particles and projections are
read into memory. Do not use it on large sets of particles !!!
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	
	if len(args)<2 :
		print "Error, please specify projection file and particles file"
		sys.exit(1)
	
	logid=E2init(sys.argv,options.ppid)
	
	em_app = EMApp()
	window = EM3DGLWidget() #TODO: see if this should be a subclass of EMSymViewerWidget instead
	explorer = EMCmpExplorer(window)
	window.set_model(explorer)
	explorer.set_data(args[0],args[1])

	em_app.show()
	em_app.execute()
	
	E2end(logid)


class EMCmpExplorer(EM3DSymModel):
	def __init__(self, gl_widget, projection_file=None,simmx_file=None,particle_file=None):
		self.init_lock = True # a lock indicated that we are still in the __init__ function
		self.au_data = None # This will be a dictionary, keys will be refinement directories, values will be something like available iterations for visual study	
		EM3DSymModel.__init__(self,gl_widget)
		self.window_title = "SimmxXplor"
		#InputEventsManager.__init__(self)
	
		self.projection_file = projection_file 	# a projection file produced by e2project3d
		self.particle_file = particle_file 		# A file containing particles to be examined

		self.current_particle = -1 				# keep track of the current particle
		self.current_projection = None 			# keep track of the current projection

		self.ptcl_display = None			# display all particles
		self.mx_display = None 					# mx display module for displaying projection and aligned particle
		self.lay = None 						# 2d plot for displaying comparison between particle and projection
		
		self.simcmp = "dot:normalize=1"
		self.align = "rotate_translate_flip"
		self.aligncmp = "dot"
		self.refine = "refine"
		self.refinecmp = "dot:normalize=1"
		self.shrink=1
		
		
	def set_data(self,projections,particles):
		'''
		Initialize data
		'''
		if not file_exists(projections): raise RuntimeError("%s does not exist" %self.projection_file)
		if not file_exists(particles): raise RuntimeError("%s does not exist" %self.particle_file)
		
		self.projection_file = projections 	# a projection file produced by e2project3d
		self.particle_file = particles 		# A file containing particles to be examined
		self.set_shrink(self.shrink)
		
	def set_shrink(self,shrink):
		"""This actually loads the data ..."""
		
		self.shrink=shrink
		# Deal with particles
		n=min(EMUtil.get_image_count(self.particle_file),800)
		self.ptcl_data=[i for i in EMData.read_images(self.particle_file,range(n)) if i!=None]
		if self.shrink>1 :
			for i in self.ptcl_data : i.process_inplace("math.meanshrink",{"n":self.shrink})
		for i in self.ptcl_data : i.process_inplace("normalize.edgemean",{})

		if self.ptcl_display==None : 
			self.ptcl_display = EMImageMXModule()
			self.ptcl_display.set_mouse_mode("App")
			QtCore.QObject.connect(self.ptcl_display,QtCore.SIGNAL("mx_image_selected"),self.ptcl_selected)		
			QtCore.QObject.connect(self.ptcl_display,QtCore.SIGNAL("module_closed"),self.on_mx_display_closed)
		self.ptcl_display.set_data(self.ptcl_data)

		# deal with projections
		self.proj_data=EMData.read_images(self.projection_file)
		if self.shrink>1 :
			for i in self.proj_data : i.process_inplace("math.meanshrink",{"n":self.shrink})
		for i in self.proj_data : i.process_inplace("normalize.edgemean",{})

		eulers = [i["xform.projection"] for i in self.proj_data]
		self.specify_eulers(eulers)
		
		for i in self.proj_data : i["cmp"]=0.0
		self.set_emdata_list_as_data(self.proj_data,"cmp")

	def get_num_particles(self):
		if self.ptcl_data==None : return 0
		return len(self.ptcl_data)
		
	def render(self):
		if self.inspector == None: self.get_inspector()
			
		EM3DSymModel.render(self)
	
	def object_picked(self,object_number):
		if object_number == self.current_projection: return
		self.current_projection = object_number
		resize_necessary = False
		if self.mx_display == None:
			self.mx_display = EMImageMXModule()
			QtCore.QObject.connect(self.mx_display,QtCore.SIGNAL("module_closed"),self.on_mx_display_closed)
			resize_necessary = True

		#if self.frc_display == None:
			#self.frc_display = EMPlot2DModule()
#			QtCore.QObject.connect(self.frc_display,QtCore.SIGNAL("module_closed"),self.on_frc_display_closed)

		self.update_display(False)

		if resize_necessary:
			get_application().show_specific(self.mx_display)
			self.mx_display.optimally_resize()
#			get_application().show_specific(self.frc_display)
#			self.frc_display.optimally_resize()
		else:
			self.mx_display.updateGL()
#			self.frc_display.updateGL()
			
		if object_number != self.special_euler:
			self.special_euler = object_number
			self.regen_dl()
			
	def update_display(self,update=True):
		'''
		Uses self.current_particle and self.current_projection to udpate the self.mx_display
		'''
		if self.mx_display == None : return
		
		if self.current_particle<0 or self.current_projection==None : return
		
		dlist=[]
		dlist.append(self.proj_data[self.current_projection].copy())	# aligned projection
		dlist[0].transform(dlist[0]["ptcl.align2d"])					
		tmp=dlist[0].process("threshold.notzero")
		dlist.append(self.ptcl_data[self.current_particle].copy())		# original particle
		dlist[1].process_inplace("normalize.toimage",{"to":dlist[0]})
		dlist.append(self.proj_data[self.current_projection].copy())	# filtered projection
		dlist[2].process_inplace("filter.matchto",{"to":dlist[1]})
		dlist[2].mult(tmp)
		dlist[2].process_inplace("normalize.toimage",{"to":dlist[1]})
		dlist.append(dlist[2].copy())									# particle with projection subtracted
		dlist[3].sub(dlist[1])
		
		#dlist.append(self.ptcl_data[self.current_particle].copy())		# same as 1 and 2 above, but with a mask
		#tmp=dlist[0].process("threshold.notzero")
		#dlist[4].mult(tmp)
		#dlist[4].process_inplace("filter.matchto",{"to":dlist[0]})
		#dlist[4].mult(tmp)
		#dlist[4].process_inplace("normalize.toimage",{"to":dlist[0]})
		#dlist.append(dlist[3].copy())
		#dlist[5].sub(dlist[0])
		
		self.mx_display.set_data(dlist)
		

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
	
	def set_alignment(self,align,aligncmp,refine,refinecmp):
		"""sets alignment algorithms and recomputes"""
		self.align = str(align)
		self.aligncmp = str(aligncmp)
		self.refine = str(refine)
		self.refinecmp = str(refinecmp)

		self.update_align()



	def set_ptcl_idx(self,idx):
		"""Select the index of the current particle to use for comparisons"""
		if self.current_particle != idx:
			self.current_particle = idx
			self.update_align()

	def update_align(self):
		if self.current_particle<0 : return
		ptcl=self.ptcl_data[self.current_particle]
		
		progress = QtGui.QProgressDialog("Computing alignments", "Abort", 0, len(self.proj_data),None)
		progress.show()
		# redetermines particle alignments
		# then we can quickly compute a series of different similarity values
		aopt=parsemodopt(self.align)
		acmp=parsemodopt(self.aligncmp)
		ropt=parsemodopt(self.refine)
		rcmp=parsemodopt(self.refinecmp)
		for i,p in enumerate(self.proj_data):
			try:
				ali=p.align(aopt[0],ptcl,aopt[1],acmp[0],acmp[1])
				if self.refine!="" :
					ropt[1]["xform.align2d"]=ali["xform.align2d"]
					ali=p.align(ropt[0],ptcl,ropt[1],rcmp[0],rcmp[1])
			except:
				print traceback.print_exc()
				QtGui.QMessageBox.warning(None,"Error","Problem with alignment parameters")
				progress.close()
				return
			p["ptcl.align2d"]=ali["xform.align2d"]
			progress.setValue(i)
			QtCore.QCoreApplication.instance().processEvents()
		
		progress.close()
		self.update_cmp()
#		self.update_display(True)
	
	def set_cmp(self,cmpstring):
		"""Select the comparator. Passed as a standard name:attr=value:attr=value string"""
		self.simcmp=str(cmpstring)
		self.update_cmp()
		
	def update_cmp(self):
		cmpopt=parsemodopt(self.simcmp)
		
		progress = QtGui.QProgressDialog("Computing similarities", "Abort", 0, len(self.proj_data),None)
		progress.show()
		ptcl=self.ptcl_data[self.current_particle]
		for i,p in enumerate(self.proj_data):
			ali=p.copy()
			ali.transform(p["ptcl.align2d"])
			try : p["cmp"]=-ptcl.cmp(cmpopt[0],ali,cmpopt[1])
			except:
				print traceback.print_exc()
				QtGui.QMessageBox.warning(None,"Error","Invalid similarity metric string, or other comparison error")
				progress.close()
				return
			progress.setValue(i)
			QtGui.qApp.processEvents()
			
		progress.close()
		self.set_emdata_list_as_data(self.proj_data,"cmp")
#		self.regen_dl(True)
		EM3DSymModel.render(self)
	
class EMSimmxXplorInspector(EMSymInspector):
	def __init__(self,target,enable_trace=False,enable_og=False) :		
		EMSymInspector.__init__(self,target,enable_trace=enable_trace,enable_og=enable_og)
		self.add_cmp_options()
		
#	def __del__(self):
#		print "simmx xplor died"
		
	def add_cmp_options(self):
		self.cmp_tab= QtGui.QWidget()
		gridl = QtGui.QGridLayout(self.cmp_tab)
		
		self.cmp_shrinkl=QtGui.QLabel("Shrink:")
		gridl.addWidget(self.cmp_shrinkl,0,1)
		
		self.cmp_shrink=QtGui.QSpinBox()
		self.cmp_shrink.setRange(1,5)
		self.cmp_shrink.setValue(1)
		gridl.addWidget(self.cmp_shrink,0,2)
		
		self.cmp_ali=QtGui.QLineEdit("rotate_translate_flip")
		gridl.addWidget(self.cmp_ali,1,0,1,2)
		
		self.cmp_alicmp=QtGui.QLineEdit("dot")
		gridl.addWidget(self.cmp_alicmp,1,2,1,2)
		
		self.cmp_refine=QtGui.QLineEdit("refine")
		gridl.addWidget(self.cmp_refine,2,0,1,2)
		
		self.cmp_refinecmp=QtGui.QLineEdit("dot:normalize=1")
		gridl.addWidget(self.cmp_refinecmp,2,2,1,2)
		
		self.cmp_realignb=QtGui.QPushButton("Set Alignment")
		gridl.addWidget(self.cmp_realignb,3,2)
		
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
		gridl.addWidget(self.cmp_combo,4,0,1,4)
		
		self.tabwidget.insertTab(0,self.cmp_tab,"Cmp")
		self.tabwidget.setCurrentIndex(0)
		
		self.connect(self.cmp_combo, QtCore.SIGNAL("currentIndexChanged(QString)"), self.cmp_changed)
		self.connect(self.cmp_realignb, QtCore.SIGNAL("clicked(bool)"), self.cmp_realign)
#		self.connect(self.cmp_shrink, QtCore.SIGNAL("valueChanged(int)"), self.ali_changed)

		
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
	
	def cmp_realign(self,b):
		"Redo alignment of current particle with new parameters"
		shrink=self.cmp_shrink.value()
		if self.target().shrink!=shrink : self.target().set_shrink(shrink)
		self.target().set_alignment(str(self.cmp_ali.text()),str(self.cmp_alicmp.text()),str(self.cmp_refine.text()),str(self.cmp_refinecmp.text()))

	
if __name__ == '__main__':
	main()