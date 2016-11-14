#!/usr/bin/env python

#
# Authors: James Michael Bell, 06/03/2015
# Copyright (c) 2015 Baylor College of Medicine
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA	2111-1307 USA
#

from EMAN2 import *
import sys
from emboxerbase import *
import subprocess as sp

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <image> <image2>....
	
	Auto Boxing strategy making use of mathematical morphology. It is still in need of some work.
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_pos_argument(name="micrographs",help="List the file to process with e2boxer here.", default="", guitype='filebox', browser="EMBoxesTable(withmodal=True,multiselect=True)",row=0,col=0,rowspan=1,colspan=3,mode="boxing,extraction")
	parser.add_header(name="boxerheader", help='Options below this label are specific to e2boxer', title="### e2boxer options ###", row=1, col=0, rowspan=1, colspan=3, mode="boxing,extraction")
	parser.add_argument("--boxsize","-b",type=int,help="Box size in pixels",default=-1, guitype='intbox', row=2, col=0, rowspan=1, colspan=3, mode="boxing,extraction")
	parser.add_argument("--xmin",type=int,default=0)
	parser.add_argument("--xmax",type=int,default=-1)
	parser.add_argument("--xstep",type=int,default=16)
	parser.add_argument("--ymin",type=int,default=0)
	parser.add_argument("--ymax",type=int,default=-1)
	parser.add_argument("--ystep",type=int,default=16)
	(options, args) = parser.parse_args()
	
	app = EMApp()
	boxer = AutoBoxer(args,options)
	boxer.show_interfaces()
	app.execute()


class AutoBoxer(EMBoxerModule):
	
	"""
	Class to automatically pick particles from a micrograph
	"""
	
	def __init__(self,micrographs,options,type="morphological",fx=0,fy=0,lx=-1,ly=-1):
		super(AutoBoxer,self).__init__(micrographs,options.boxsize)
		self.box_list = MorphBoxList(self)
		self.add_tool(MorphBoxingTool)
		for i in xrange(len(self.file_names)):
			self.set_current_file_by_idx(i)
			f = self.current_file()
			if self.get_num_boxes(f) == 0:
				hdr = EMData(f,0,True).get_attr_dict()
				if options.xmin == 0: fx = int(options.boxsize)
				if options.ymin == 0: fy = int(options.boxsize)
				if options.xmax == -1: lx = int(hdr['nx']-options.boxsize)
				if options.ymax == -1: ly = int(hdr['ny']-options.boxsize)
				boxes = []
				for y in xrange(fy,ly,options.ystep):
					for x in xrange(fx,lx,options.xstep):
						boxes.append([x,y,type])
				self.add_boxes(boxes)
	
	def morph_tiles(self,options):
		# apply filters in filter test to each box individually
		pass
	
	def classify_tiles(self,options):
		# perform a classification to distinguish particles from nonsense
		pass
	
	def remove_nonsense(self,options):
		# removes bad particles (maybe store in a separate file for easy reference?)
		pass
	
	def adjust_filters(self):
		# Just like the filtertool. Allow user to adjust all parameters in the morph panel.
		pass
	
	@staticmethod
	def get_num_boxes(file_name):
		"""
		A static function for getting the number of boxes associated with each file
		"""
		box_list = EMBoxList()
		box_list.load_boxes_from_database(file_name)
		return len(box_list)


class MorphBoxList(EMBoxList):

	"""
	Class to show filtered particles. Needs an on/off switch!
	"""

	def __init__(self,target):
		super(MorphBoxList,self).__init__(target)
	
	def get_particle_images(self,image_name,box_size):
		return [self.process_box(box.get_image(image_name,box_size,"normalize.edgemean")) for box in self.boxes]
	
	def process_box(self,image):
		img = image.copy()
		img.process_inplace('math.meanshrink',{'n':2})
		img.process_inplace('filter.highpass.gauss',{'sigma':0.015})
		#img.process_inplace('normalize.edgemean')
		#img.process_inplace('math.sigma',{'value1':15.0,'value2':0.0})
		#img.process_inplace('normalize.edgemean')
		img.process_inplace('filter.lowpass.gauss',{'sigma':0.12})
		#img.process_inplace('normalize.edgemean')
		#img.process_inplace('threshold.belowtozero',{'minval':0.098})
		#img.process_inplace('math.gradient.magnitude')
		img.process_inplace('math.gradient.magnitude')
		img.process_inplace('morph.object.density',{'thresh':5.0})
		#img.process_inplace('morph.majority',{'nmaj':1,'thresh':0.0})
		#img.process_inplace('morph.object.density',{'thresh':0.0})
		#img.process_inplace('mask.addshells.multilevel',{'nshells':3})
		#img.process_inplace('threshold.belowtozero',{'minval':250})
		#img.process_inplace('histogram.bin',{'nbins':3})
		#img.process_inplace('normalize.maxmin')
		img.process_inplace('normalize.edgemean')
		return self.center_particle_in_box(img)
	
	def center_particle_in_box(self,img):
		# will iteratively move box to center on peak of density
		# possibly a good opportunity to use the image gradient? (contour plot)
		
		#img.process_inplace('xform.center')
		
		return img


class MorphBoxingTool(EMBoxingTool):

	BOX_TYPE = "morphological"
	EMBox.set_box_color(BOX_TYPE,[1,1,1])

	def __init__(self,target):
		super(MorphBoxingTool,self).__init__(target)
		self.target = weakref.ref(target)
		self.moving = None
		self.panel_object = None
		self.moving_data = None

	def get_widget(self):
		if self.panel_object == None:
			self.panel_object = MorphBoxingPanel(self)
		return self.panel_object.get_widget()

	def icon(self):
		from PyQt4 import QtGui
		return QtGui.QIcon(get_image_directory() + "white_box.png")

	def set_panel_object(self,panel): self.panel_object = panel
	def unique_name(self): return MorphBoxingTool.BOX_TYPE

	def set_current_file(self,file_name,active_tool=False):
		"""
		If the behavior of this Handler does not if the file changes, but the function needs to be supplied
		"""
		pass

	def get_2d_window(self): return self.target().get_2d_window()

	def mouse_down(self,event) :
		m = self.get_2d_window().scr_to_img((event.x(),event.y()))
		box_num = self.target().detect_box_collision(m)
		from PyQt4.QtCore import Qt
		if box_num == -1:
			if event.modifiers()&Qt.ShiftModifier : return # the user tried to delete nothing
			box_num = self.target().add_box(m[0],m[1],MorphBoxingTool.BOX_TYPE)
			if self.panel_object.auto_center_checkbox.isChecked():
				self.try_to_center_ref(box_num)
			self.moving=[m,box_num]
		else:
			box = self.target().get_box(box_num)
			if box.type == MorphBoxingTool.BOX_TYPE:
		 		if event.modifiers()&Qt.ShiftModifier :
					self.target().remove_box(box_num)
				else:
					self.moving=[m,box_num]
			else:
				raise EMUnknownBoxType,box.type

	def mouse_drag(self,event) :
		m=self.get_2d_window().scr_to_img((event.x(),event.y()))
		from PyQt4.QtCore import Qt
		if event.modifiers()&Qt.ShiftModifier:
			box_num = self.target().detect_box_collision(m)
			if ( box_num != -1):
				box = self.target().get_box(box_num)
				if box.type ==  MorphBoxingTool.BOX_TYPE:
					self.target().remove_box(box_num)
				else:
					raise EMUnknownBoxType,box.type

		elif self.moving != None:
			oldm = self.moving[0]
			self.target().move_box(self.moving[1],m[0]-oldm[0],m[1]-oldm[1])
			self.moving[0] = m

	def mouse_up(self,event) :
		if self.moving != None:
			self.target().box_released(self.moving[1])
		self.moving=None

	def mouse_move(self,event): pass

	def clear_all(self):
		self.target().clear_boxes([MorphBoxingTool.BOX_TYPE],cache=True)

	def moving_ptcl_established(self,box_num,x,y):
		box = self.target().get_box(box_num)
		if box.type != MorphBoxingTool.BOX_TYPE:
			raise EMUnknownBoxType,box.type
		self.moving_data = [x,y,box_num]

	def move_ptcl(self,box_num,x,y,scale):
		if self.moving_data == None: return
		dx = self.moving_data[0] - x
		dy = y - self.moving_data[1]
		self.target().move_box(self.moving_data[2],dx,dy)
		self.moving_data = [x,y,self.moving_data[2]]

	def release_moving_ptcl(self,box_num,x,y):
		if self.moving_data == None: return
		self.target().box_placement_update_exclusion_image_n(box_num)
		self.moving_data = None

	def delete_ptcl(self,box_num):
		box = self.target().get_box(box_num)
		if box.type != MorphBoxingTool.BOX_TYPE:
			raise EMUnknownBoxType,box.type
		self.target().remove_box(box_num)

	def get_unique_box_types(self):
		return [MorphBoxingTool.BOX_TYPE]

	def boxes_erased(self,list_of_boxes):
		"""
		No need to act here for the morph boxing tool - everything is fine?
		"""
		pass

	def xform_center_propagate(self,box,image_name,template,box_size):
		"""
		Centers a box that was generated in a shrunken image by getting the 'real particle' out of the large
		image on disk and doing a ccf with the template - then I just find the peak and use that to center
		@param box a list like [x,y,type,float] - only x and y are used
		@param image_name the name of the image we're operating on
		@param template the template correctly scaled to the have the same angstrom per pixel as the image (named image_name) stored on disk
		@param box_size the size of the box used to center
		Returns the dx and dy parameters, i.e. does not actually alter the box
		"""
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
		return dx,dy

	def try_to_center_ref(self,box_num): #Modified from that in SwarmBoxer
		box = self.target().get_box(box_num)
		box_size = self.target().get_box_size()
		img_filename = self.target().current_file()
		ptcl = box.get_image(img_filename, box_size)
		centered_ptcl = ptcl.process("xform.centeracf")
		dx,dy = self.xform_center_propagate([box.x,box.y],img_filename,centered_ptcl,box_size)
		self.target().move_box(box_num, dx, dy)


class MorphBoxingPanel:
	
	def __init__(self,target):
		self.target = weakref.ref(target)
		self.widget = None
	
	def get_widget(self):
		if self.widget == None:
			from PyQt4 import QtCore, QtGui, Qt
			self.widget = QtGui.QWidget()
			vbl = QtGui.QVBoxLayout(self.widget)
			vbl.setMargin(0)
			vbl.setSpacing(6)
			vbl.setObjectName("vbl")
			self.auto_center_checkbox = QtGui.QCheckBox("Auto-center")
			self.clear=QtGui.QPushButton("Clear")
			vbl.addWidget(self.auto_center_checkbox)
			vbl.addWidget(self.clear)
			QtCore.QObject.connect(self.clear, QtCore.SIGNAL("clicked(bool)"), self.clear_clicked)
		return self.widget
	
	def clear_clicked(self,val):
		self.target().clear_all()


class ErasingPanel: # copied for ideas for the morph panel

	def __init__(self,target,erase_radius=128):
		self.busy = True
		self.erase_radius = erase_radius
		self.target = weakref.ref(target)
		self.erase_rad_edit = None
		self.widget = None
		self.busy = False

	def set_erase_radius(self, erase_rad_edit):
		self.busy=True
		self.erase_radius = erase_rad_edit
		if self.erase_rad_edit != None: self.erase_rad_edit.setValue(erase_rad_edit)
		self.busy=False

	def get_widget(self):
		if self.widget == None:
			from PyQt4 import QtCore, QtGui, Qt
			self.widget = QtGui.QWidget()
			vbl = QtGui.QVBoxLayout(self.widget)
			vbl.setMargin(0)
			vbl.setSpacing(6)
			vbl.setObjectName("vbl")

			hbl = QtGui.QHBoxLayout()
			hbl.addWidget(QtGui.QLabel("Erase Radius:"))
			from valslider import ValSlider
			self.erase_rad_edit = ValSlider(None,(0.0,1000.0),"")
			self.erase_rad_edit.setValue(int(self.erase_radius))
			self.erase_rad_edit.setEnabled(True)
			hbl.addWidget(self.erase_rad_edit)

			self.unerase = QtGui.QCheckBox("Unerase")
			self.unerase.setChecked(False)

			vbl.addLayout(hbl)
			vbl.addWidget(self.unerase)
			QtCore.QObject.connect(self.erase_rad_edit,QtCore.SIGNAL("sliderReleased"),self.new_erase_radius)
			QtCore.QObject.connect(self.unerase,QtCore.SIGNAL("clicked(bool)"),self.unerase_checked)

		return self.widget

	def new_erase_radius(self, erase_rad_edit):
		if self.busy: return
		self.target().set_erase_radius(erase_rad_edit)

	def unerase_checked(self,val):
		if self.busy: return
		self.target().toggle_unerase(val)


if __name__ == "__main__":
	main()
