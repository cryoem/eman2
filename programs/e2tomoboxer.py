#!/usr/bin/env python

#
# Author: David Woolford 04/16/2009 (woolford@bcm.edu)
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

from optparse import OptionParser
import sys
import os
from EMAN2 import file_exists, gimme_image_dimensions3D,EMANVERSION,EMData,Region,Transform
from pyemtbx.boxertools import Cache,get_idd_image_entry,set_idd_image_entry, Box
from emapplication import get_application, EMStandAloneApplication,EMQtWidgetModule
from emimage2d import EMImage2DModule
import weakref
from emshape import EMShape
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import Qt

tomo_db_name = "bdb:e2tomoboxercache#"

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <image>
	
Manual particle selection from tomograms. This version is specifically aimed at cubic boxes
for tomographic analysis."""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--gui",action="store_true",help="Start the GUI for interactive boxing",default=False)
	parser.add_option("--boxsize","-B",type="int",help="Box size in pixels",default=-1)
	parser.add_option("--auto","-A",type="string",action="append",help="Autobox using specified method: ref, grid, db, cmd",default=[])

	(options, args) = parser.parse_args()
	
	
	application = EMStandAloneApplication()
	tomo_mediator= EMTomoBoxerModule(args[0])
	
	tomo_mediator.show_guis()
#	QtCore.QObject.connect(gui, QtCore.SIGNAL("module_idle"), on_idle)
	application.execute()
	
#	print args[0]
#	a = TomogramZProjectionImage(args[0])
#	a.get_image()


class TomogramZProjectionImage:
	'''
	Sums the tomograph (3D image) along z, stores and returns the result
	'''
	def __init__(self,file_name):
		'''
		@param file_name the name of the file on disk - should be a 3D tomogram
		@exception RuntimeError raised if the file_name is not a file on disk
		@exception RuntimeError raised if the file on disk is not 3D (nz > 1)
		'''
		if not file_exists(file_name): raise RuntimeError("The file %s does not exist" %file_name)
		nx,ny,nz = gimme_image_dimensions3D(file_name)
		if nz <= 1: raise RuntimeError("The file must be 3D (nz>1)")
		self.file_name = file_name
		# this will be either None or an EMData object
		self.image = get_idd_image_entry(self.file_name,"tomogram_z_projection",tomo_db_name)
	def get_image_name(self):
		return self.file_name
	
	def get_construction_argument(self):
		return self.file_name
	
	def get_image(self):
		if self.image == None:
			tmp = EMData(self.file_name)
			self.image = tmp.project("standard",Transform())
			self.image.write_image("result.mrc")
#			nx,ny,nz = gimme_image_dimensions3D(self.file_name)
#			self.image = EMData(nx,ny)
#			self.image.to_zero()
#			tmp = EMData()
#			for i in range(nz):
#				print i
#				tmp.read_image(self.file_name,0,False,Region(0,0,i,nx,ny,i+1))
#				tmp.set_size(nx,ny)
#				self.image = self.image + tmp
#			
			# cache to database
			set_idd_image_entry(self.file_name,"tomogram_z_projection",self.image,tomo_db_name)	
		return self.image
	
ZProjectionCache = Cache(TomogramZProjectionImage)


class EMTBBoxManipulations():
	'''
	A class that responds to added, moved and removed box signals emitted
	by the EMImage2DModule
	'''
	def __init__(self,target,main_2d_window):
		self.target = weakref.ref(target) # prevent a strong cycle
		self.main_2d_window = main_2d_window
		self.moving = None
		
		self.__connect_signals_to_slots()
	def __connect_signals_to_slots(self):
		
		'''
		connects the signals of the main 2D window to the slots of this object
		'''
		qt_target = get_application().get_qt_emitter(self.main_2d_window)
		
		QtCore.QObject.connect(qt_target,QtCore.SIGNAL("mousedown"),self.mouse_down)
		QtCore.QObject.connect(qt_target,QtCore.SIGNAL("mousedrag"),self.mouse_drag)
		QtCore.QObject.connect(qt_target,QtCore.SIGNAL("mouseup")  ,self.mouse_up  )
		QtCore.QObject.connect(qt_target,QtCore.SIGNAL("keypress"),self.key_press)
		QtCore.QObject.connect(qt_target,QtCore.SIGNAL("mousewheel"),self.mouse_wheel)

	def __box_collision_detect(self,coord_data):
		'''
		A wrapper function that could easily disappear
		@param coord_data a coordinate interms of the main image's pixel coordinate system
		'''
		return self.target().detect_box_collision(coord_data)
	
	def mouse_down(self,event) :
		'''
		@param a QtGui.QMouseEvent sent from the EMImage2DModule
		'''
		coord_data = self.main_2d_window.scr_to_img((event.x(),event.y()))
		box_num =  self.__box_collision_detect(coord_data)
		if box_num == -1:
			if not self.main_2d_window.coords_within_image_bounds(coord_data):	return
			
			if event.modifiers() & Qt.ShiftModifier : return # the user tried to delete nothing
			
			# If we get here, we need to add a new box
			box_size = self.target().get_box_size()
			lbx = coord_data[0]-box_size/2
			lby = coord_data[1]-box_size/2
			box = [lbx,lby,lbx+box_size,lby+box_size]
			
			self.target().add_box(box)
			self.moving=[coord_data,box_num]
			
			self.main_2d_window.set_active(box_num,.9,.9,.2)
			self.main_2d_window.updateGL()
			
		
		elif event.modifiers()&Qt.ShiftModifier : # box removal
			self.target().remove_box(box_num)
		else:
			# if we make it here than the we're moving a box
			box = self.target().get_box(box_num)
			self.moving=[coord_data,box_num]
			self.target().set_active_box(box_num)
							
	def mouse_drag(self,event) :
		'''
		@param a QtGui.QMouseEvent sent from the EMImage2DModule
		'''
		
		coord_data = self.main_2d_window.scr_to_img((event.x(),event.y()))
		
		if event.modifiers()&Qt.ShiftModifier:
			box_num = self.__box_collision_detect(coord_data)
			if ( box_num != -1):
				self.target().remove_box(box_num)
				#self.target().mouse_click_update_ppc()
			
		elif self.moving != None:
			old_coord = self.moving[0]
			
			self.target().move_box(self.moving[1],coord_data[0]-old_coord[0],coord_data[1]-old_coord[1])
			self.moving[0] = coord_data
	
	def mouse_up(self,event) :
		'''
		@param a QtGui.QMouseEvent sent from the EMImage2DModule
		'''
		if self.moving != None:
			self.moving=None
		
	def key_press(self,event):
		'''
		@param a QtGui.QKeyEvent sent from the EMImage2DModule
		'''
		pass 
	def mouse_wheel(self,event):
		'''
		@param a QtGui.QMouseEvent sent from the EMImage2DModule
		'''
		pass
	

class EMBoxList:
	'''
	Meant to store the primary data associated with a box, namely its bottom left and top right
	coordinates. Stores coordinates as a list of length 4, e.g.
	[bottom left x, bottom left y, top right x, top right y]
	Provides additional collision detection functionality
	Provides additional box size changing and box moving functionality
	Supports list-style iterator support 
	'''
	def __init__(self):
		self.boxes = [] # main data container
		self.current_iter = 0 # iterator support
		self.max_idx = 0 # iterator support
		
	def append(self,box_coords):
		'''
		@param box_coords should be a list like [bottom left x, bottom left y, top right x, top right y]
		@exception RuntimeError raised if box_coords is not a list
		@exception RuntimeError raised if box_coords is not of length 4
		'''
		if not isinstance(box_coords,list): raise RuntimeError("box_coords should be a list")
		if len(box_coords) != 4: raise RuntimeError("box_coords should have precisely 4 entries")
		
		self.boxes.append(box_coords)
	
	def detect_collision(self,coords):
		'''
		@param coords a coordinate list, the first two values are the coordinates
		@exception RuntimeError raised if coords is not a list
		@exception RuntimeError raised if coords is not at least length 2
		The last exception is an aberration because we already have a framework where coordinate
		data is embedded in a list that is longer than 2 in length 
		@return the list-coordinate of the first box that collided with the coordinate, else -1 
		'''
		if not isinstance(coords,list) and not isinstance(coords,tuple): raise RuntimeError("coords should be a list")
		if len(coords) < 2: raise RuntimeError("coords should have 2 entries")
		
		for box_num,box in enumerate(self.boxes):
			if coords[0]<box[0] or coords[0]>box[2] or coords[1]<box[1] or coords[1]>box[3]:
				# no collision
				continue
			# if we make it here there has been a collision, the box already exists
			return box_num
		
		return -1
	
	def change_box_size(self,new_box_size):
		'''
		@param the box size - this is much more like diameter than it is radius
		@exception RuntimeError raised if the box_size is less than 1
		#exception RuntimeError raised if the first internally stored box is not a square - to alert developers if things change
		'''
		if new_box_size < 1: raise RuntimeError("The boxsize must be at least 1")
		if len(self.boxes) == 0: return # this is fine, the calling program might not know and there is not box_size attribute
		
		box_0 = self.boxes[0]
		old_box_size = box_0[2]-box_0[0] # note no consideration 
		check_box_size = box_0[3]-box_0[1]
		if old_box_size != check_box_size:
			raise NotImplementedError("This code does not currently handle uneven box dimensions - developer alert") 
	
		box_shift = (new_box_size - old_box_size)/2 # this should be fine for odd sizes
		new_boxes = []
		for box in self.boxes:
			lbx = box[0]-box_shift
			lby = box[1]-box_shift
			new_boxes.append([lbx,lby,lbx+new_box_size,lby+new_box_size])
		
		self.boxes = new_boxes
			
		
	
	def move_box(self,idx,dx,dy):
		'''
		move box at idx in the x and y direction
		@param idx the index of the object to be moved
		@param dx the amount by which to move in the x direction
		@param dy the amount to move in the y direction
		'''
		box = self.boxes[idx]
		box[0] += dx
		box[1] += dy
		box[2] += dx
		box[3] += dy
	
	def __getitem__(self,idx):
		'''
		Support for [] syntax
		@idx the index of the object which is to be retrieved
		'''
		return self.boxes[idx]
	
	def __len__(self):
		'''
		support for standard len
		'''
		return len(self.boxes)
	
	def __iter__(self):
		''' 
		Iteration support - initialization 
		'''
		self.current_iter = 0
		self.max_idx = len(self.boxes)
		return self
	
	def next(self):
		''' 
		Iteration support
		@return the next element in the iteration
		@raise StopIteration exception when the end of the container is reached
		'''
		if self.current_iter >= self.max_idx:
			raise StopIteration
		else:
			self.current_iter += 1
			return self.boxes[self.current_iter-1]
		
	def pop(self,idx):
		'''
		@param idx the index of object that is to be popped
		@return the popped object
		'''
		return self.boxes.pop(idx)
	
class EMBoxDisplayShapes:
	def __init__(self): pass
		
	def gen_shapes(box_list,shape_string="rectpoint"):
		'''
		@param box_list should be an instance of an EMBoxList
		@exception RuntimeError raise if box_list if not an EMBoxList
		
		'''
		shapes = {}
		
		for i,box in enumerate(box_list):
			shape = EMShape([shape_string,0,0.9,0.4,box[0],box[1],box[2],box[3],2.0])
			shapes[i] = shape
			
		return shapes
	gen_shapes = staticmethod(gen_shapes)
	
class EMTomoBoxerModule:
	'''
	This module is essentially a Mediator (see GoF) - it coordinates the activities of several EMAN2 modules
	that would otherwise not necessary interact. For tomographic boxing we need a main display (for
	these purposes that is a projection along z) and a display for the isolated ("boxed") objects
	'''
	def __init__(self,file_name):
		'''
		@file_name the name of a file on disk
		@exception RuntimeError raised if the file does not exist
		'''
		if not file_exists(file_name): raise RuntimeError("The file (%s) does not exist" %file_name)
		self.file_name = file_name
		self.box_size = 256
		self.box_list = EMBoxList()
		self.active_box = -1 # an index storing which box the user had clicked on
		
		# initialize mouse  handlers first
		self.__init_signal_handlers()
		# initialize windows
		self.__init_main_2d_window()
		# initialize the inspector
		self.__init_inspector()
		
	def __init_signal_handlers(self):
		'''
		Initialize mouse event handlers - "events" are actually pseudo-events because
		this module is not a QtWidget, the "events" are actually triggered from signals
		that are emitted by member objects and intercepted by slots which are member functions.
		'''
		self.signal_slot_handlers = {}
		

	def __init_inspector(self):
		self.inspector = EMTomoBoxerInspectorModule(self)
		
	def __init_main_2d_window(self):
		'''
		Called internally to initialize the main 2D image viewing window.
		Also connects signals to their corresponding slots
		'''
		from PyQt4 import QtCore
		
		self.main_2d_window = EMImage2DModule()
		self.main_2d_window.set_data(ZProjectionCache.get_image_directly(self.file_name),self.file_name)
		get_application().show_specific(self.main_2d_window)
		
		self.main_2d_window.setWindowTitle(self.file_name)
		
		self.signal_slot_handlers["box_manipulation"] = EMTBBoxManipulations(self,self.main_2d_window)
		
				
	def show_guis(self):
		'''
		Called in the context of an application - generally one creates and application, then
		creates an EMTomoBoxerModule, then tells the EMTomoBoxerModule to "show_guis", meaning
		all of the interfaces should appear
		'''
		get_application().show_specific(self.main_2d_window)
		self.main_2d_window.optimally_resize()
		
		get_application().show_specific(self.inspector)
		
	def set_box_size(self,box_size):
		'''
		@param the box size - this is much more like diameter than it is radius
		@exception RuntimeError raised if the box_size is less than 1
		@exception RuntimeError raised if the the box_size is bigger than the current image's dimensions
		'''
		if box_size < 1: raise RuntimeError("The boxsize must be at least 1")
		image = ZProjectionCache.get_image_directly(self.file_name)
		if box_size > image.get_xsize() or box_size > image.get_ysize():
			raise RuntimeError("The box size can not be larger than any of the image's dimensions")
	
		self.box_list.change_box_size(box_size)
		self.box_size = box_size
		self.__refresh_box_display()
	
	def get_box_size(self):
		''' @return the current box size'''
		return self.box_size
	def get_current_image_name(self):
		''' @return the name of the current file '''
		return self.file_name
	def get_shape_string(self):
		''' @return the string name used to generate EMShapes '''
		return "rectpoint"
		
	def detect_box_collision(self,coords):
		'''
		See EMBoxList.detect_collision help
		'''
		return self.box_list.detect_collision(coords)
	
	def add_box(self,box):
		'''
		See EMBoxList.append help
		'''
		self.box_list.append(box)
		self.active_box = len(self.box_list)-1
		self.__refresh_box_display()
	
	def remove_box(self,box_idx):
		'''
		See EMBoxList.pop help
		'''
		self.box_list.pop(box_idx)
		self.__refresh_box_display()
		
	def move_box(self,box_idx,dx,dy):
		'''
		See EMBoxList.move_box for help
		'''
		self.box_list.move_box(box_idx,dx,dy)
		self.__refresh_box_display()
		
	def get_box(self,box_idx):
		'''
		See EMBoxList.__getitem__ for help
		'''
		return self.box_list[box_idx]
	
	def set_active_box(self,box_idx):
		self.active_box = box_idx
		self.main_2d_window.set_active(self.active_box,.9,.9,.1)
		self.main_2d_window.shapechange = 1
		self.main_2d_window.updateGL()
	
	def __refresh_box_display(self):
		'''
		Refreshes the display of boxes. Typically called when a box has been added or moved.
		Tells concerned interfaces to update their displays
		'''
		shapes = EMBoxDisplayShapes.gen_shapes(self.box_list,self.get_shape_string())
		self.main_2d_window.set_shapes(shapes,1.0)
		self.main_2d_window.set_active(self.active_box,.9,.9,.1)
		print self.active_box
		self.main_2d_window.updateGL()
		

class EMTomoBoxerInspectorModule(EMQtWidgetModule):
	'''
	'''
	def __init__(self,target):
		self.application = weakref.ref(get_application())
		self.widget = EMTomoBoxerInspector(target)
		EMQtWidgetModule.__init__(self,self.widget)
		
	def get_desktop_hint(self):
		return "inspector"
		
class EMTomoBoxerInspector(QtGui.QWidget):
	def get_desktop_hint(self):
		return "inspector"
	
	def __init__(self,target) :
		
		QtGui.QWidget.__init__(self,None)
		self.target=weakref.ref(target)
		
		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		
		box_size_hbl=QtGui.QHBoxLayout()
		box_size_hbl.setMargin(0)
		box_size_hbl.setSpacing(2)
		
		self.box_size_label = QtGui.QLabel("Box Size:",self)
		box_size_hbl.addWidget(self.box_size_label)
		self.pos_int_validator = QtGui.QIntValidator(self)
		self.pos_int_validator.setBottom(1)
		self.box_size = QtGui.QLineEdit(str(self.target().get_box_size()),self)
		self.box_size.setValidator(self.pos_int_validator)
		box_size_hbl.addWidget(self.box_size)
		
		self.vbl.addLayout(box_size_hbl)
		
		self.connect(self.box_size,QtCore.SIGNAL("editingFinished()"),self.new_box_size)
		
	def new_box_size(self):
		box_size=int(self.box_size.text())
		self.target().set_box_size(box_size)


if __name__ == "__main__":
	main()
		
		
		