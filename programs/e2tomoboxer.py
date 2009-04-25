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
from EMAN2 import file_exists, gimme_image_dimensions3D,EMANVERSION,EMData,Region,Transform,get_image_directory,db_check_dict,db_open_dict,db_close_dict
from pyemtbx.boxertools import Cache,get_idd_image_entry,set_idd_image_entry, Box
from emapplication import get_application, EMStandAloneApplication,EMQtWidgetModule
from emimage2d import EMImage2DModule
import weakref
from emshape import EMShape
from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import Qt
import math
import EMAN2db

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
#	a = TomogramProjection(args[0])
#	a.get_image()

class ShrunkenTomogram:
	'''
	Takes care of the automatic shrinking of tomogram, if necessary
	It may be true that the shrink factor was 1, in which case shrink didn't really occur.
	Could be made generic for other applications if the need arises
	'''
	db_name = "shrunken_image"
	def __init__(self,file_name):
		'''
		@param file_name the name of the file, which is a tomographic reconstruction sitting on disk
		@exception RuntimeError raised if the file_name is not a file on disk
		@exception RuntimeError raised if the file on disk is not 3D (nx,ny,nz > 1)
		'''
		if not file_exists(file_name): raise RuntimeError("The file %s does not exist" %file_name)
		nx,ny,nz = gimme_image_dimensions3D(file_name)
		if nx <=1 or ny <= 1 or nz <= 1: raise RuntimeError("The file must be 3D (nx,ny,nz>1)")
		
		self.file_name = file_name
		self.image = get_idd_image_entry(self.file_name,ShrunkenTomogram.db_name,tomo_db_name)
		
	def get_shrunken_resident_name(file_name):
		'''
		asks whether or not this file name has a corresponding shrunken image in the
		e2tomoboxer cache, and if so, for its name (embeded in list - [db_name, image_index] ([string,int])
		return None if there is no shrunken cached image, or [string,int], which can be used as the arguments to EMData.read_image(string,int)
		'''

		db_name =  tomo_db_name+ShrunkenTomogram.db_name
		if not db_check_dict(db_name): 
		#		print "db failed",db_name
			return [None,None] #
		
		db = db_open_dict(db_name)
		
		if db.has_key("maxrec"):
			for idx in range(0,db["maxrec"]+1):
				header = db.get_header(idx)
				try:
					if header["shrunken_from"] == file_name:
						return [db_name,idx]
				except:
					pass
				
		db_close_dict(db_name)
		return [None,None]
	
	get_shrunken_resident_name = staticmethod(get_shrunken_resident_name)
		
	def get_construction_argument(self):
		'''
		Used by the cache for retrieving the correct object
		'''
		return self.file_name
	
	def get_image(self):
		if self.image == None:
			global HOMEDB
			HOMEDB=EMAN2db.EMAN2DB.open_db()
			HOMEDB.open_dict("e2tomoboxer_preferences")
			db = HOMEDB.e2tomoboxer_preferences	
			target_max_dim = db.get("largest_allowable_dimension",dfl=1024)
			nx,ny,nz = gimme_image_dimensions3D(self.file_name)
			shrink = 1
			for val in [nx,ny,nz]:
				s = int(math.ceil(float(val)/target_max_dim))
				if s > shrink: shrink = s
			
			if shrink > 1:
				shrink_thread = ShrinkThread(self.file_name,shrink)
				progress = QtGui.QProgressDialog(shrink_thread.get_action_string(), "abort", 0, shrink_thread.get_total_steps(),None)
				progress.setWindowTitle("Shrinking tomogram by %i" %shrink)
				progress.show()
				tally = -1
				shrink_thread.start()
				while 1:
					get_application().processEvents()
					if progress.wasCanceled():
						shrink_thread.quit()
						progress.close()
						# probably the calling program should do something about this
						# the image returned is going to be None
						break
					
					if shrink_thread.get_progress() > tally:
						tally = shrink_thread.get_progress()
						progress.setValue(tally)
						progress.setLabelText(shrink_thread.get_action_string())
											
					if shrink_thread.get_progress() == shrink_thread.get_total_steps()-1:
						progress.close()
						self.image = shrink_thread.get_image()
						self.image.set_attr("shrunken_by",shrink)
						self.image.set_attr("shrunken_from",self.file_name)
						set_idd_image_entry(self.file_name,ShrunkenTomogram.db_name,self.image,tomo_db_name) #cache to disk
						break
			else:
				self.image = EMData()
				self.image.read_image(self.file_name,0)
		
		return self.image

class ShrinkThread(QtCore.QThread):
	'''
	A thread for running the shrink operation - so the user can cancel if it's taking too long or
	they think twice.
	Designed to be used with a progress dialog - you can get the progress (get_progress) to update the
	progress bar, and get the action string (get_action_string) to update the label on the progress dialog.
	See ShrunkenTomogram for an example of how this object is used
	Once the job is done call get_image to get the result 
	'''
	def __init__(self,file_name,shrink):
		'''
		@param file_name the name of the file on disk, or in bdb format
		@param shrink the shrink factor, should be an integer greater than 1
		@exception RuntimeError raised if the file_name is not a file on disk
		@exception RuntimeError raised if shrink is not an int
		@exception RuntimeError raised if shrink is not greater than 1
		'''
		if not file_exists(file_name): raise RuntimError("The file %s does not exist" %file_name)
		if not isinstance(shrink,int) or shrink < 2: raise RuntimeError("The shrink factor must be an integer greater than 1")
		
		QtCore.QThread.__init__(self)
		self.file_name = file_name # store the file name
		self.shrink = shrink # store the shrink
		self.current_step = -1 # so the calling function can ask about the progress and update a progress dialog
		self.total_steps = 2 # the total number of steps I chose this object to have (read,shrink,write)
		self.image = None # the final image - so the calling object can retrieve it
		self.action_string = "Reading image" # so the progress dialog can display a useful string
		
	def get_total_steps(self):
		''' @return the total number of operations this object will perform, i.e. for initializing a progress dialgo
		'''
		return self.total_steps
	
	def get_progress(self):
		''' @return the current progress as an int. Before run is called this is -1.'''
		return self.current_step
	
	def get_action_string(self):
		''' @return a useful string saying what this is currently doing, i.e. for a progress dialog'''
		return self.action_string
	
	def get_image(self):
		''' @return the shrunken image, or None if for example the thread was terminatea'''
		return self.image
	
	def run(self):
		''' run - overwrites QtThread.run - all work is done here
		'''
		try:
			self.image = EMData()
			self.image.read_image(self.file_name,0)
			print "read image"
			self.current_step = 0
			self.action_string = "Shrinking by %i" %self.shrink
			self.image.process_inplace("math.meanshrink",{"n":self.shrink})
			print "shrunk the image"
			self.current_step = 1
		except: pass
		self.exit() # is this necessary
		
		
class TomogramProjection:
	'''
	Sums a tomograph (3D image) along the narrowest dimension, this is not necessarily the z dimension
	Caches this compressed image on disk and in memory. Disk cache is to prevent the projection operation from
	occurring more than once.
	Designed to be used the the boxertools.Cache object - in particular using the get_image_directly interface 
	Stores the dimension that was compressed as a header attribute ("sum_direction")
	Stores the number of pixels along the compressed dimension as a header attribute ("sum_pixels")
	'''
	def __init__(self,file_name):
		'''
		@param file_name the name of the file on disk - should be a 3D tomogram
		@exception RuntimeError raised if the file_name is not a file on disk
		@exception RuntimeError raised if the file on disk is not 3D (nx,ny,nz > 1)
		'''
		if not file_exists(file_name): raise RuntimeError("The file %s does not exist" %file_name)
		nx,ny,nz = gimme_image_dimensions3D(file_name)
		if nx <=1 or ny <= 1 or nz <= 1: raise RuntimeError("The file must be 3D (nx,ny,nz>1)")
		self.file_name = file_name
		# this will be either None or an EMData object
		self.image = get_idd_image_entry(self.file_name,"tomogram_z_projection",tomo_db_name)
	def get_image_name(self):
		'''
		@return the file name of the image 
		'''
		return self.file_name

	def get_construction_argument(self):
		'''
		@return the file name of the image 
		'''
		return self.file_name
	
	def get_image(self):
		'''
		@return a version of the tomogram that has been compressed/projected along the narrowest dimension 
		If this projection did not previously exist then it is created, and the "sum_direction" and "sum_pixels" 
		header attributes are created so that the original dimensions can be recalled if necessary.
		'''
		if self.image == None:
			
			st = ShrunkenTomogram(self.file_name)
		
			tmp = st.get_image()
			if tmp == None: # this means the shrink operation was cancelled
				return None # the calling function should be sophisticated to handle this!
		
			nx,ny,nz = tmp["nx"],tmp["ny"],tmp["nz"]
			if ny < nz:
				self.image = tmp.process("misc.directional_sum",{"direction":"y"})
				self.image.set_size(nx,nz,1)
				self.image.set_attr("sum_pixels",tmp.get_ysize())
				self.image.set_attr("sum_direction","y")
				
			elif nx < nz:
				self.image = tmp.process("misc.directional_sum",{"direction":"x"})
				self.image.set_size(ny,nz,1)
				self.image.set_attr("sum_pixels",tmp.get_xsize())
				self.image.set_attr("sum_direction","x")
			else:
				self.image = tmp.process("misc.directional_sum",{"direction":"z"})
				self.image.set_size(nx,ny,1)
				self.image.set_attr("sum_pixels",tmp.get_zsize())
				self.image.set_attr("sum_direction","z")
				
			set_idd_image_entry(self.file_name,"tomogram_z_projection",self.image,tomo_db_name)	
		return self.image
	
TomoProjectionCache = Cache(TomogramProjection)


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
		coord_data = list(self.main_2d_window.scr_to_img((event.x(),event.y())))
		coord_data[0] = int(coord_data[0])
		coord_data[1] = int(coord_data[1])
		box_num =  self.__box_collision_detect(coord_data)
		if box_num == -1:
			if not self.main_2d_window.coords_within_image_bounds(coord_data):	return
			
			if event.modifiers() & Qt.ShiftModifier : return # the user tried to delete nothing
			
			self.target().add_box(coord_data)
			self.moving=[coord_data,box_num]
		
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
		
		coord_data = list(self.main_2d_window.scr_to_img((event.x(),event.y())))
		coord_data[0] = int(coord_data[0])
		coord_data[1] = int(coord_data[1])
		
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
	
class EMTBPositionBoxManipulations():
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
		return self.target().detect_side_view_box_collision(coord_data)
	
	def mouse_down(self,event) :
		'''
		@param a QtGui.QMouseEvent sent from the EMImage2DModule
		'''
		coord_data = list(self.main_2d_window.scr_to_img((event.x(),event.y())))
		coord_data[0] = int(coord_data[0])
		coord_data[1] = int(coord_data[1])
		box_num =  self.__box_collision_detect(coord_data)
		if box_num != -1:

#			# if we make it here than the we're moving a box
			box = self.target().get_box(box_num)
			self.moving=[coord_data,box_num]
#			self.target().set_active_box(box_num)
							
	def mouse_drag(self,event) :
		'''
		@param a QtGui.QMouseEvent sent from the EMImage2DModule
		'''
		
		coord_data = list(self.main_2d_window.scr_to_img((event.x(),event.y())))
		coord_data[0] = int(coord_data[0])
		coord_data[1] = int(coord_data[1])
#		print "mouse drag"#
		if self.moving != None:
			old_coord = self.moving[0]
			
			self.target().move_box(self.moving[1],coord_data[0]-old_coord[0],0,coord_data[1]-old_coord[1])
			self.moving[0] = coord_data
#	
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
	
class EMCoordList:
	'''
	Meant to store the primary data associated with a box, namely its bottom left and top right
	coordinates. Stores central coordinates as a list of length 3, e.g. [center x, center y, center z]
	Provides additional collision detection functionality
	Provides additional box size changing and box moving functionality
	Supports list-style iterator support 
	'''
	def __init__(self):
		self.coords = [] # main data container, entries are lists of length 4
		self.current_iter = 0 # iterator support
		self.max_idx = 0 # iterator support
		
	def append(self,box_coords):
		'''
		@param box_coords should be a list like [center_x,center_y] or [center_x,center_y,center_z] - if unsure of z don't specify it, will set it -1 internally
		@exception RuntimeError raised if box_coords is not a list
		@exception RuntimeError raised if box_coords is not of length 2 or 3
		'''
		if not isinstance(box_coords,list): raise RuntimeError("box_coords should be a list")
		if len(box_coords) not in [2,3]: raise RuntimeError("box_coords should have 2 or 3 entries")
		
		if len(box_coords) == 2:
			box_coords.append(-1)
			
		self.coords.append(box_coords)
	
	def detect_collision(self,coords,box_size):
		'''
		@param coords a coordinate list, the first two values are the coordinates
		@param boxsize the size of the boxed centered on the coordinate used as the basis of collision detection
		@exception RuntimeError raised if coords is not a list or a tuple
		@exception RuntimeError raised if coords is not at least length 2
		The last exception is an aberration because we already have a framework where coordinate
		data is embedded in a list that is longer than 2 in length 
		@return the list-coordinate of the first box that collided with the coordinate, else -1 
		'''
		if not isinstance(coords,list) and not isinstance(coords,tuple): raise RuntimeError("coords should be a list")
		if len(coords) < 2: raise RuntimeError("coords should have 2 entries")
		
		for box_num,box in enumerate(self.coords):
			lbx = box[0]-box_size/2
			lby = box[1]-box_size/2
			rtx = lbx+box_size
			rty = lby+box_size
			if coords[0]<lbx or coords[0]>rtx or coords[1]<lby or coords[1]>rty:
				# no collision
				continue
			# if we make it here there has been a collision, the box already exists
			return box_num
		
		return -1			
		
	def move_box(self,idx,dx,dy,dz=0):
		'''
		move box at idx in the x and y and z direction
		@param idx the index of the object to be moved
		@param dx the amount by which to move in the x direction
		@param dy the amount to move in the y direction
		@param dz the amount to move in the z direction
		'''
		box = self.coords[idx]
		box[0] += dx
		box[1] += dy
		box[2] += dz
#		box[3] += dy
	
	def __getitem__(self,idx):
		'''
		Support for [] syntax
		@idx the index of the object which is to be retrieved
		'''
		return self.coords[idx]
	
	def __len__(self):
		'''
		support for standard len
		'''
		return len(self.coords)
	
	def __iter__(self):
		''' 
		Iteration support - initialization 
		'''
		self.current_iter = 0
		self.max_idx = len(self.coords)
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
			return self.coords[self.current_iter-1]
		
	def pop(self,idx):
		'''
		@param idx the index of object that is to be popped
		@return the popped object
		'''
		return self.coords.pop(idx)
	
class EMBoxDisplayShapes:
	def __init__(self): pass
		
	def gen_shapes(box_list,box_size,shape_string="rectpoint"):
		'''
		@param box_list should be an instance of an EMCoordList
		@exception RuntimeError raise if box_list if not an EMCoordList
		
		'''
		shapes = {}
		
		for i,box in enumerate(box_list):
			lbx = int(box[0])-box_size/2+(box_size+1)%2
			lby = int(box[1])-box_size/2+(box_size+1)%2
			rtx = lbx+box_size
			rty = lby+box_size
			shape = EMShape([shape_string,0,0.9,0.4,lbx,lby,rtx,rty,2.0])
			shapes[i] = shape
			
		return shapes
	gen_shapes = staticmethod(gen_shapes)
	
class EMTomoBoxerModule:
	'''
	This module is essentially a Mediator (see Design Patterns) - it coordinates the activities of several EMAN2 modules
	that would otherwise not necessary interact. For tomographic boxing we need a main display (for
	these purposes that is a projection along z) and a display for the isolated ("boxed") objects. The user can manipulate the
	location of the boxed sub tomograms and the display in both windows is updated
	'''
	def __init__(self,file_name):
		'''
		@file_name the name of a file on disk
		@exception RuntimeError raised if the file does not exist
		'''
		if not file_exists(file_name): raise RuntimeError("The file (%s) does not exist" %file_name)
		self.file_name = file_name
		self.box_size = 64
		self.coord_list = EMCoordList()
		self.active_box = -1 # an index storing which box the user had clicked on
		self.positioning_2d_window = None # will eventually be used to position boxed images in the z direction
		
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
		Also creates a signal handler for managing box addition, movement and removal
		'''
		
		self.main_2d_window = EMImage2DModule()
		image = TomoProjectionCache.get_image_directly(self.file_name)
		if image == None: # the user cancelled the shrink operation
			sys.exit(1)
		self.main_2d_window.set_data(image,self.file_name)
		self.main_2d_window.setWindowTitle(self.file_name)
		
		self.signal_slot_handlers["box_manipulation"] = EMTBBoxManipulations(self,self.main_2d_window)
	
	def __init_positioning_2d_window(self):
		'''
		Called internally to intiialize the 2D window used to center the box in the z direction
		Also connects signals to their corresponding 
		'''
		self.positioning_2d_window = EMImage2DModule()
		self.positioning_2d_window.setWindowTitle("Z positioning")
		
		self.signal_slot_handlers["position_manipulation"] = EMTBPositionBoxManipulations(self,self.positioning_2d_window)
		
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
		image = TomoProjectionCache.get_image_directly(self.file_name)
		if box_size > image.get_xsize() or box_size > image.get_ysize():
			raise RuntimeError("The box size can not be larger than any of the image's dimensions")
	
		self.box_size = box_size
		self.__refresh_box_display()
	
	def get_box_size(self):
		''' @return the current box size'''
		return self.box_size
	
	def get_anticipated_output_box_size(self):
		''' Look at the images, see if there's a large unshrunken one, use the scaling factors to guess what a good output is'''
		[file_name,idx] = ShrunkenTomogram.get_shrunken_resident_name(self.file_name)
		if file_name ==  None: 
			scale = 1
		else:
			db = db_open_dict(file_name)
			hdr = db.get_header(idx)
			scale = hdr["shrunken_by"]
		
		return scale*self.box_size
	
	def get_current_image_name(self):
		''' @return the name of the current file '''
		return self.file_name
	def get_shape_string(self):
		''' @return the string name used to generate EMShapes '''
		return "rectpoint"
		
	def detect_box_collision(self,coords):
		'''
		See EMCoordList.detect_collision help
		'''
		return self.coord_list.detect_collision(coords,self.box_size)
	
	def detect_side_view_box_collision(self,coords):
		if self.active_box == None: return
		
		active_box =  self.coord_list[self.active_box]
		
		by = active_box[2]-self.box_size/2
		ty = by+self.box_size
		if coords[1]<by or coords[1]>ty: return -1
		else: return self.active_box
	
	def add_box(self,coord):
		'''
		See EMCoordList.append help
		'''
		self.coord_list.append(coord)
		self.active_box = len(self.coord_list)-1
		self.__refresh_box_display()
	
	def remove_box(self,box_idx):
		'''
		See EMCoordList.pop help
		'''
		self.coord_list.pop(box_idx)
		self.__refresh_box_display()
		
	def move_box(self,box_idx,dx,dy,dz=0):
		'''
		See EMCoordList.move_box for help
		'''
		self.coord_list.move_box(box_idx,dx,dy,dz)
		self.__refresh_box_display()
		
	def get_box(self,box_idx):
		'''
		See EMCoordList.__getitem__ for help
		'''
		return self.coord_list[box_idx]
	
	def set_active_box(self,box_idx):
		self.active_box = box_idx
		self.main_2d_window.set_active(self.active_box,.9,.9,.1)
		self.main_2d_window.shapechange = 1
		self.main_2d_window.updateGL()
		self.__refresh_positioner()
	
	def __refresh_box_display(self):
		'''
		Refreshes the display of boxes. Typically called when a box has been added or moved.
		Tells concerned interfaces to update their displays
		'''
		shapes = EMBoxDisplayShapes.gen_shapes(self.coord_list,self.box_size,self.get_shape_string())
		self.main_2d_window.set_shapes(shapes,1.0)
		self.main_2d_window.set_active(self.active_box,.9,.9,.1)
		self.main_2d_window.updateGL()
		
		self.__refresh_positioner()
		
	def __refresh_positioner(self):
		if self.active_box != -1:
			coord = self.coord_list[self.active_box]
			image = TomoProjectionCache.get_image_directly(self.file_name)
			a = EMData()
			direction = image.get_attr("sum_direction")
			half_box = self.box_size/2
			
			[file_name,idx] = ShrunkenTomogram.get_shrunken_resident_name(self.file_name)
			if file_name ==  None: 
				file_name = self.file_name
				idx = 0
			if direction == "z":
				nz = image.get_attr("sum_pixels")
				nx = self.box_size
				ny = self.box_size
				region = Region(coord[0]-half_box,coord[1]-half_box,0,nx,ny,nz)
				
				a.read_image(file_name,idx,False,region)
				b = a.process("misc.directional_sum",{"direction":"y"})
				b.set_size(int(nx),int(nz),1)
			elif direction == "y":
				ny = image.get_attr("sum_pixels")
				nx = self.box_size
				nz = self.box_size
				region = Region(coord[0]-half_box,0,coord[1]-half_box,nx,ny,nz)
				a.read_image(file_name,idx,False,region)
				b = a.process("misc.directional_sum",{"direction":"z"})
				b.set_size(int(nx),int(ny),1)
			elif direction == "x":
				nx = image.get_attr("sum_pixels")
				ny = self.box_size
				nz = self.box_size
				region = Region(0,coord[0]-half_box,coord[1]-half_box,nx,ny,nz)
				a.read_image(file_name,idx,False,region)
				b = a.process("misc.directional_sum",{"direction":"z"})
				b.set_size(int(nx),int(ny),1)
			
			if self.positioning_2d_window == None:
				self.__init_positioning_2d_window()
				self.positioning_2d_window.set_data(b)
				get_application().show_specific(self.positioning_2d_window)
				self.positioning_2d_window.optimally_resize()
			else:
				self.positioning_2d_window.set_data(b)
			
			
			if coord[2] == -1:
				com = b.calc_center_of_mass()
				coord[2] = b.get_ysize() - com[1]
				
			l = int(coord[2])  - self.box_size/2+(self.box_size+1)%2
			r = l + self.box_size
			shape = EMShape(["rectpoint",0,0.9,0.4,0,l,self.box_size,r,2.0])
			shapes = {0:shape}
			self.positioning_2d_window.set_shapes(shapes,1)
			self.positioning_2d_window.updateGL()
	
	def write_stack(self,out_file_name, input_file_name,xout,yout,zout):
		[file_name,idx] = ShrunkenTomogram.get_shrunken_resident_name(self.file_name)
		if file_name ==  None: 
			file_name = self.file_name
			idx = 0
		
		out_idx = -1
		image = TomoProjectionCache.get_image_directly(self.file_name)
		direction = image.get_attr("sum_direction")
		if input_file_name == file_name: # then we do not need to scale
			half_box = self.box_size/2
			box = self.box_size
			fname = file_name
			scale = 1.0
		else: # we need to scale
			db = db_open_dict(file_name)
			hdr = db.get_header(idx)
			orig_file = hdr["shrunken_from"]
			scale = hdr["shrunken_by"]
			half_box = self.box_size/2*scale
			box = self.box_size*scale
			fname = orig_file
		
		progress = QtGui.QProgressDialog("Writing_files", "abort", 0, 2*len(self.coord_list),None)
		progress.show()
		tally = -1
		for i,coord in enumerate(self.coord_list):	
			a = EMData()
			if direction == "z":
				region = Region(scale*coord[0]-xout/2,scale*coord[1]--yout/2,scale*coord[2]-zout/2,xout,yout,zout)
				a.read_image(fname,idx,False,region)
			elif direction == "y":
				region = Region(scale*coord[0]-xout/2,scale*coord[2]-zout/2,scale*coord[1]-yout/2,xout,zout,yout)
				a.read_image(fname,idx,False,region)
			elif direction == "x":
				region = Region(scale*coord[2]-zout/2,scale*coord[0]-xout/2,scale*coord[1]-yout/2,zout,xout,yout)
				a.read_image(fname,idx,False,region)
			else:
				progress.close()
				raise NotImplementedException
			
			if progress.wasCanceled():
				# Should probably delete the file(s)
				break
			
			tally += 1
			progress.setValue(tally)

			if isinstance(out_file_name,list):
				print out_file_name[i]
				a.write_image(out_file_name[i],0)
			else:
				a.write_image(out_file_name,-1)
			
			tally += 1
			progress.setValue(tally)
			
		progress.close()
			
		
	def write_images(self,out_file_names, input_file_name,xout,yout,zout):
		self.write_stack(out_file_names, input_file_name,xout,yout,zout)
	
	def run_output_dialog(self):
		if self.get_num_boxes() == 0: return
		
		[file_name,idx] = ShrunkenTomogram.get_shrunken_resident_name(self.file_name)
		file_names = [self.file_name]
		file_indices = [0]
		if file_name !=  None: 
			file_names.append(file_name)
			file_indices.append(idx)
		
		self.form = TomoBoxerOutputForm(self,file_names,file_indices)
		self.form.run_form()
		
	def get_num_boxes(self):
		''' @return the number of boxes currently in the image '''
		return len(self.coord_list)

from emsprworkflow import EMFormTask
from emdatastorage import ParamDef
class TomoBoxerOutputForm(EMFormTask):
	def __init__(self,target,file_names,file_indices):
		EMFormTask.__init__(self)
		self.window_title = "e2tomoboxer write output"
		self.preferred_size = [240,240]
		self.target = weakref.ref(target)
		self.file_names = file_names # a list of file_names
		self.file_indices = file_indices # list the file indices, import if accessing images from stacks
		self.output_file_names = None
	def get_params(self):
		params = []
		box_size = self.target().get_anticipated_output_box_size()
		xout = ParamDef(name="xout",vartype="int",desc_short="Output nx",desc_long="The number of pixels of the output sub regions will have in the x direction",property=None,defaultunits=box_size,choices=None)
		yout = ParamDef(name="yout",vartype="int",desc_short="Output ny",desc_long="The number of pixels of the output sub regions will have in the y direction",property=None,defaultunits=box_size,choices=None)
		zout = ParamDef(name="zout",vartype="int",desc_short="Output nz",desc_long="The number of pixels of the output sub regions will have in the z direction",property=None,defaultunits=box_size,choices=None)
		
		pdatafile = ParamDef(name="data_file",vartype="string",desc_short="Extract from",desc_long="The name of the file that will be used for the subvolume extraction.",property=None,defaultunits=[],choices=self.file_names)
		pwrite_stack = ParamDef(name="write_stack",vartype="boolean",desc_short="Write to 3D stack",desc_long="Tick if you want the output files written to a stack.",property=None,defaultunits=False,choices=None)
		poutname = ParamDef(name="out_file",vartype="string",desc_short="Output file name",desc_long="If writing to a stack this is literally the name of the output file.\nIf writing many output image names the name you specify will be altered, e.g. by adding '00' for the first image and so on.\nYou will be told if any files already exist.\nEMAN2 style bdb syntax is supported.",property=None,defaultunits="something.hdf",choices=None)
	
		params.append([xout,yout,zout])
		params.append(pdatafile)
		params.append(pwrite_stack)
		params.append(poutname)
		return params
	
	def on_form_ok(self,params):
		from EMAN2 import get_supported_3d_stack_formats
		# check file_names
		error_message = []
		error_message.extend(self.__out_file_name_ok(params["out_file"],params["write_stack"]))
		error_message.extend(self.__out_sizes_ok(params["xout"],params["yout"],params["zout"]))
						
		if len(error_message) != 0:
			from emsprworkflow import EMErrorMessageDisplay
			EMErrorMessageDisplay.run(error_message)
			return
	
		# if we make it here we can go ahead and write
		if params["write_stack"]:
			self.target().write_stack(self.output_file_names,params["data_file"],params["xout"],params["yout"],params["zout"])
		else:
			self.target().write_images(self.output_file_names,params["data_file"],params["xout"],params["yout"],params["zout"])
			
		self.form.closeEvent(None)
		self.form = None
		self.emit(QtCore.SIGNAL("task_idle"))
	
	def __out_sizes_ok(self,nx,ny,nz):
		'''
		Check to make sure nx,ny and nz are greater than 1
		'''
		l1_warning = False
		g2048_warning = False
		error_messages = []
		for i in nx,ny,nz:
			if i < 1 and not l1_warning:
				error_messages.append("dimensions must be greater than 1")
				l1_warning = True
			if i > 2048 and not g2048_warning:
				error_messages.append("dimensions must be less than 2048. If this is a problem contact developers")
				g2048_warning = True
		
		return error_messages

	def __out_file_name_ok(self,out_file,write_stack):
		'''
		checks to make sure the out_file_name is ok
		returns a list of error messages. If the list is empty everything is fine
		'''
		error_messages = []
		fine = True
		
		if len(out_file) > 4 and out_file[:4] == "bdb:":
			if db_check_dict(out_file):
				error_messages.append("The bdb file name you specified exists")
				fine = False
		elif file_exists(out_file):
			# should prompt the user to see if they want to append or overwrite
			error_messages.append("The file name you specified exists")
			fine = False
		else:
			d = out_file.split('.')
			if len(d) < 2: 
				error_messages.append("The file name has no type")
				fine = False
			else:
				from EMAN2 import get_supported_3d_stack_formats,get_supported_3d_formats
				if write_stack:
					fmts = get_supported_3d_stack_formats()
				else:
					fmts = get_supported_3d_formats()
			
				if d[-1] not in fmts:
					e = "The format you specifed is not supported. We support ["
					for f in fmts:
						if f != fmts[0]:
							e += ", "
						e+= f
					e+=", bdb:]"
					error_messages.append(e)
					fine = False
					
		if not write_stack:
			self.output_file_names = []
			bdb_type = False
			if len(out_file) > 4 and out_file[:4] == "bdb:": 
				bdb_type = True
				name = out_file_name
				ext = ''
			else:
				# this is guaranteed to be safe
				d = out_file.split('.')
				ext = "."+d[-1]
				from EMAN2 import get_file_tag
				name = get_file_tag(out_file)
			
			n = self.target().get_num_boxes()
			digits = len(str(n))
			for i in range(self.target().get_num_boxes()):
				si = str(i)
				this_digits = len(si)
				s = '_'
				for j in range(digits-this_digits): s+= '0'
				s += si
				this_name =  name+s+ext
				if file_exists(this_name):
					error_messages.append("%s exists" %this_name)
				else:
					self.output_file_names.append(this_name)
		else:
			self.output_file_names = out_file
					
		return error_messages
				
			
				
			

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
		self.setWindowIcon(QtGui.QIcon(get_image_directory() +"eman.png"))
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
		
		self.gen_output_but=QtGui.QPushButton("Write output")
		self.vbl.addWidget(self.gen_output_but)
		
		
		self.connect(self.box_size,QtCore.SIGNAL("editingFinished()"),self.new_box_size)
		self.connect(self.gen_output_but,QtCore.SIGNAL("clicked(bool)"),self.target().run_output_dialog)

	def new_box_size(self):
		box_size=int(self.box_size.text())
		self.target().set_box_size(box_size)


if __name__ == "__main__":
	main()
		
		
		