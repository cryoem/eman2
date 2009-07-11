#!/usr/bin/env python
#
# Author: David Woolford (woolford@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine


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

from emboxerbase import *
from valslider import ValSlider
from EMAN2 import BoxingTools,gm_time_string,Transform
from pyemtbx.boxertools import CoarsenedFlattenedImageCache,FLCFImageCache

SWARM_TEMPLATE_MIN = TEMPLATE_MIN # this comes from emboxerbase

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <image> <image2>....

a refactoring of e2boxer
"""

	parser = OptionParser(usage=usage,version=EMANVERSION)
	parser.add_option("--boxsize","-B",type="int",help="Box size in pixels",default=128)
	
	(options, args) = parser.parse_args()
	
	error_message = check(options,args)
	if len(error_message) > 0:
		error = "\n"
		for e in error_message:
			error += "Error: "+e +"\n"
		parser.error(error)
		
	args = [abs_path(arg) for arg in args] # always try to use full file names 

	application = EMStandAloneApplication()
#	QtCore.QObject.connect(gui, QtCore.SIGNAL("module_idle"), on_idle)
	
	module = EMBoxerModule(args,options.boxsize)
	module.show_interfaces()
	# this is an example of how to add your own custom tools:
	module.add_2d_window_mouse_tool(SwarmEventHandling,SwarmPanel,particle_diameter=options.boxsize)
	application.execute()


class SwarmPanel:
	def __init__(self,target,particle_diameter=128):
		self.busy = True
		self.particle_diameter = particle_diameter
		self.target = weakref.ref(target)
		self.widget = None
		self.busy = False
		
	def icon(self):
		from PyQt4 import QtGui
		return QtGui.QIcon(get_image_directory() + "swarm_icon.png")
		
	def get_widget(self):
		if self.widget == None:
			from PyQt4 import QtCore, QtGui, Qt
			self.widget = QtGui.QWidget()
			vbl = QtGui.QVBoxLayout(self.widget)
			vbl.setMargin(0)
			vbl.setSpacing(6)
			vbl.setObjectName("vbl")
			
			hbl = QtGui.QHBoxLayout()
			hbl.addWidget(QtGui.QLabel("Particle Diameter:"))
			
			self.ptcl_diam_edit = QtGui.QLineEdit(str(self.particle_diameter))
			hbl.addWidget(self.ptcl_diam_edit)
			vbl.addLayout(hbl)
			
			self.thrbut = QtGui.QRadioButton("Threshold")
			self.selbut = QtGui.QRadioButton("Selective")
			self.selbut.setChecked(True)
			self.morselbut = QtGui.QRadioButton("More Selective")
			
			self.methodhbox = QtGui.QHBoxLayout()
			self.methodhbox.addWidget(self.thrbut)
			self.methodhbox.addWidget(self.selbut)
			self.methodhbox.addWidget(self.morselbut)
			
			self.groupbox = QtGui.QGroupBox("Auto Box Method")
			self.groupbox.setLayout(self.methodhbox)
			vbl.addWidget(self.groupbox)
			
			self.advanced_hbl2 = QtGui.QHBoxLayout()
			self.enable_interactive_params  = QtGui.QCheckBox("Interactive Threshold")
			self.enable_interactive_params.setChecked(False) #FIXME should this be related to the state of the autoboxer?
			self.thr = ValSlider(None,(0.0,3.0),"")
			self.thr.setValue(1.0)
			self.thr.setEnabled(False)
			self.advanced_hbl2.addWidget(self.enable_interactive_params)
			self.advanced_hbl2.addWidget(self.thr)
			vbl.addLayout(self.advanced_hbl2)
			
			hbl_ww = QtGui.QHBoxLayout()
			self.update_template = QtGui.QCheckBox("Update Template")
			self.update_template.setChecked(True)
			hbl_ww.addWidget(self.update_template)
			
			
			self.view_template=QtGui.QPushButton("View Template")
			hbl_ww.addWidget(self.view_template)
			self.clear=QtGui.QPushButton("Clear")
			hbl_ww.addWidget(self.clear)
		
			vbl.addLayout(hbl_ww)
		
			QtCore.QObject.connect(self.ptcl_diam_edit,QtCore.SIGNAL("editingFinished()"),self.new_ptcl_diam)
			QtCore.QObject.connect(self.update_template,QtCore.SIGNAL("clicked(bool)"),self.update_template_checked)
			QtCore.QObject.connect(self.clear, QtCore.SIGNAL("clicked(bool)"), self.clear_clicked)
			
		return self.widget
	
	def new_ptcl_diam(self):
		if self.busy: return
		print "yo"
		
	def update_template_checked(self,val):
		if self.busy: return
		print "yo"
		
	def hide(self):
		if self.widget != None:
			self.widget.hide()
	
	def clear_clicked(self,val):
		self.target().clear_all()

class SwarmShrunkenImageMediator:
	def __init__(self,template_radius,subsample_rate):
		self.template_radius = template_radius
		self.subsample_rate = subsample_rate
	
	def get_template_radius(self): return self.template_radius
	def get_subsample_rate(self): return self.subsample_rate

def gen_rot_ave_template(image_name,ref_boxes,shrink,box_size,iter=4):
	
	ptcls = []
	mediator = SwarmShrunkenImageMediator(box_size/(2*shrink),shrink)
	
	real_box_size = box_size/shrink
	averages = []
	for box in ref_boxes:
		if box.in_template == False: continue # it's excluded from the template
		
		xc = box.x/shrink-real_box_size/2
		yc = box.y/shrink-real_box_size/2
		r = Region(xc,yc,real_box_size,real_box_size)
		image = CoarsenedFlattenedImageCache.get_image(box.image_name,mediator)
		particle = image.get_clip(r)
		
		particle.process_inplace("normalize.edgemean")
		ptcls.append(particle)
		
	if len(ptcls) == 0:
		raise RuntimeError("No boxes are flagged as contributing to the template")
	
	# bootstrap the original average 
	ave = ptcls[0].copy()
	for i in range(1,len(ptcls)):
		ta = ptcls[i].align("translational",ave)
		ave.add(ta)
	
	ave.process_inplace("xform.centeracf")
	ave.process_inplace("math.radialaverage")
	ave.process_inplace("normalize.edgemean")
	ave.process_inplace("mask.sharp",{'outer_radius':ave.get_xsize()/2})
	averages.append(ave)
	averages[-1].set_attr("creation_time_stamp", gm_time_string())
	# 4 is a magic number
	for n in range(0,iter):
		t = []
		for idx,particle in enumerate(ptcls):
			ta = particle.align("translational",ave)
			t.append(ta)
			
	
		ave = t[0].copy()
		for i in range(1,len(ptcls)):
			ave.add(t[i])
			
		ave.process_inplace("xform.centeracf")
		ave.process_inplace("math.radialaverage")
		ave.process_inplace("normalize.edgemean")
		
		# edge normalize here SL before
		ave.process_inplace("mask.sharp",{'outer_radius':ave.get_xsize()/2})
		averages.append(ave)
		averages[-1].set_attr("creation_time_stamp", gm_time_string())
		
	return averages

class SwarmFLCFImageMediator:
	def __init__(self,template_radius,subsample_rate,template_image):
		self.template_radius = template_radius
		self.subsample_rate = subsample_rate
		self.template_object = SwarmFLCFImageMediator.TemplateObject(template_image)
	
	def get_template_radius(self): return self.template_radius
	def get_subsample_rate(self): return self.subsample_rate
	def get_template_object(self): return self.template_object
	
	class TemplateObject:
		def __init__(self,template_image):
			self.template_image = template_image
			
		def get_template_ts(self): return self.template_image["creation_time_stamp"]
		
		def get_template(self): return self.template_image

		def get_radius(self): return 
class SwarmBox:
	def __init__(self,x,y,image_name,in_template=True,profile=None):
		self.x = x
		self.y = y
		self.in_template = in_template
		self.image_name = image_name
		self.profile=profile
		self.peak_x = None
		self.peak_y = None
		self.peak_score = None
		
	def update_picking_data(self,mediator):
		
		shrink = mediator.get_subsample_rate()
		x = int(self.x/shrink)
		y = int(self.y/shrink)
		search_radius = mediator.get_template_radius()
		correlation = FLCFImageCache.get_image(self.image_name,mediator)
		
		peak_location = BoxingTools.find_radial_max(correlation,x,y,search_radius )
		peak_location2 = BoxingTools.find_radial_max(correlation,peak_location[0],peak_location[1],search_radius )
		if (peak_location != peak_location2):
			# there is no local max, this particle will be ignorned
			self.profile = None
			self.peak_score = None
			self.peak_x = None
			self.peak_y = None
			return
	
		# store the peak location
		self.peak_x = peak_location[0]
		self.peak_y = peak_location[1]
	
		# store the correlation value at the correlation max
		self.peak_score = correlation.get(self.peak_x,self.peak_y)
		
		if self.peak_score != None:
			# store the profile
			self.profile = BoxingTools.get_min_delta_profile(correlation,self.peak_x,self.peak_y, mediator.get_template_radius() )
		

def compare_box_correlation(box1,box2):
	c1 = box1[3]
	c2 = box2[3]
	if c1 > c2: return -1
	elif c1 == c2: return 0
	else: return 1


class SwarmBoxer:
	THRESHOLD = "Threshold"
	SELECTIVE = "Selective"
	MORESELECTIVE = "More Selective"
	def __init__(self,particle_diameter=128):
		self.particle_diameter = particle_diameter
		self.update_template = True
		self.pick_mode = SwarmBoxer.SELECTIVE	# the autobox method - see EMData::BoxingTools for more details
		self.profile_mode = BoxingTools.CmpMode.SWARM_AVERAGE_RATIO # the way a profile is generated
		BoxingTools.set_mode(self.profile_mode)
		self.interactive_threshold_enabled = False
		self.interactive_threshold = None
		self.profile_max = 0.8 # this is a percentage - it stops the profile trough point from going to the end
		
		self.reset()
		
	
	def reset(self):
		self.ref_boxes = []
		self.template_change = False
		self.templates = None
		self.profile = None
		self.peak_score = None
		self.profile_trough_point = None
		self.profile_max = 0.8 # this is a percentage - it stops the profile trough point from going to the end
		
	def add_ref(self,x,y,image_name):
		self.ref_boxes.append(SwarmBox(x,y,image_name,self.update_template))
		if not self.template_change: self.template_change = self.update_template # this is safe
	
	def get_subsample_rate(self): 
		return int(math.ceil(float(self.particle_diameter)/float(SWARM_TEMPLATE_MIN)))
	
	def auto_box(self,image_name,exclusion_image):
		
		self.target().get_2d_window().updateGL()
		shrink = self.get_subsample_rate()
		if self.template_change:
			self.templates = gen_rot_ave_template(image_name,self.ref_boxes,self.get_subsample_rate(),self.particle_diameter)
			self.template_change = False
#			for template in self.templates:
#				template.write_image("template.hdf",-1)
	   	  
		
		mediator = SwarmFLCFImageMediator(self.particle_diameter/(shrink*2),shrink,self.templates[-1])
		correlation_image = FLCFImageCache.get_image(image_name,mediator)
		
		exclusion_shrink = exclusion_image.get_attr("shrink")
		
		if shrink != exclusion_shrink:
			# to do: test this
			print "weird",shrink,exclusion_shrink,self.particle_diameter, SWARM_TEMPLATE_MIN,TEMPLATE_MIN
			rescale = float(exclusion_shrink)/shrink
			oldx = exclusion_image.get_xsize()
			oldy = exclusion_image.get_ysize()
			newx = correlation_image.get_xsize()
			newy = correlation_image.get_ysize()
			r = Region((oldx-newx)/2,(oldy-newy)/2,newx,newy)
			t = Transform()
			t.set_scale(rescale)
			if rescale > 1.0:
				exclusion_image.clip_inplace(r)
				exclusion_image.transform(t)
			else:
				exclusion_image.transform(t)
				exclusion_image.clip_inplace(r)
				
		for box in self.ref_boxes:
			box.update_picking_data(mediator)
				
		self.update_opt_picking_data()
		
		if self.pick_mode == SwarmBoxer.THRESHOLD: mode = 0
		elif self.pick_mode == SwarmBoxer.SELECTIVE: mode = 1
		elif self.pick_mode == SwarmBoxer.MORESELECTIVE: mode = 2
				
		searchradius = self.templates[-1].get_xsize()/2
		correlation = FLCFImageCache.get_image(image_name,mediator)
		soln = BoxingTools.auto_correlation_pick(correlation,self.peak_score,searchradius,self.profile,exclusion_image,self.profile_lower_point,mode)
		
		scaled_template = self.templates[-1].process("math.transform.scale",{"scale":shrink,"clip":self.particle_diameter})
		scaled_template.process_inplace("xform.centeracf")

		boxes = []
		for b in soln:
			x = b[0]
			y = b[1]
			xx = int(x*shrink)
			yy = int(y*shrink)
			type = "auto"
			peak_score = correlation.get(x,y)
			box = [xx,yy,type,peak_score]
			self.center_propagate(box,image_name,scaled_template,self.particle_diameter)
			boxes.append(box)

	   	boxes.sort(compare_box_correlation) # sorting like this will often put large ice contaminations in a group, thanks Pawel Penczek
		self.target().add_boxes(boxes)
	
	def update_opt_picking_data(self):
		self.profile = None
		self.peak_score = None
		self.profile_lower_point = None
		
		for box in self.ref_boxes:
			if box.peak_score == None: continue
			
			if self.peak_score == None or box.peak_score < self.peak_score:
				self.peak_score = box.peak_score
		
			from copy import copy
			if self.profile == None: self.profile = copy(box.profile)
			else:
				if len(self.profile) != len(box.profile): raise RuntimeError("This should never happen")
				
				profile = box.profile
				for j in xrange(0,len(self.profile)):
					if profile[j] < self.profile[j]: self.profile[j] = profile[j]
				
		max_radius =  int(len(self.profile)*self.profile_max)
		tmp = self.profile[0]
		self.profile_lower_point = 0
		for i in range(1,max_radius):
			if self.profile[i] > tmp and tmp > 0:
				tmp = self.profile[i]
				self.profile_lower_point = i
		
	def center_propagate(self,box,image_name,template,box_size):
	  	global BigImageCache
	  	image = BigImageCache.get_image_directly(image_name)
			 
		xc = box[0]-box_size/2
		yc = box[1]-box_size/2
		r = Region(xc,yc,box_size,box_size)
		particle = image.get_clip(r)
		ccf  = particle.calc_ccf(template)
		trans = ccf.calc_max_location_wrap(particle.get_xsize()/2,particle.get_ysize()/2,0)
		dx = trans[0]
		dy = trans[1]
		box[0] += trans[0]
		box[1] += trans[1]
					
class SwarmEventHandling(SwarmBoxer):
	'''
	A class that knows how to handle mouse erase events for a GUIBox
	'''	
	def __init__(self,target,panel_object=None,particle_diameter=128):
		SwarmBoxer.__init__(self,particle_diameter)
		self.target = weakref.ref(target)
		
		self.ref_box_name = "ref"
		self.auto_box_name = "auto"
		
	def unique_name(self): return "Swarm"
	
	def get_2d_window(self): return self.target().get_2d_window()
	def set_panel_object(self,panel): self.panel_object = panel
		
	def mouse_move(self,event):
		pass
		
	def mouse_wheel(self,event):
		pass
	def mouse_down(self,event):
		m = self.get_2d_window().scr_to_img((event.x(),event.y()))
		box_num = self.target().detect_box_collision(m)
		from PyQt4.QtCore import Qt
		if box_num == -1:
			if event.modifiers()&Qt.ShiftModifier : return # the user tried to delete nothing
			self.target().add_box(m[0],m[1],type="ref")
			self.target().clear_boxes(["auto"])
			self.get_2d_window().updateGL()
			self.add_ref(m[0],m[1],self.target().current_file())
			self.auto_box(self.target().current_file(),self.target().get_exclusion_image(mark_boxes=True))
		else:
			pass
#			box_type = self.target().get_box_type(box_num)
#			if box_type in [self.auto_box_name,self.ref_box_name]:
#				from PyQt4.QtCore import Qt
#				if event.modifiers()&Qt.ShiftModifier :
#					# remove the box
#					self.target().remove_box(box_num)
#				else:
#					# if we make it here than the we're moving a box
#					self.moving=[m,box_num]
#					self.target().moving_box_established(box_num)
#			else:
#				print 'should change mouse handler now'
			
			
	def mouse_drag(self,event) :
		pass
	def mouse_up(self,event) :
		pass
	
	def clear_all(self):
		self.target().clear_boxes(["ref","auto"])
		self.reset()

if __name__ == "__main__":
	main()
