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


'''
This is a a really basic version of boxer that can be copied and used as the basis of developing something more advanced
Design is meant to be granular and logical, i.e. easy to add to.

To see how it works just run it from the commandline:
./e2boxerbase.py test.mrc --boxsize=128 # single image handling
-or-
./e2boxerbase.py test.mrc test2.mrc test3.mrc --boxsize=128 # multiple image handling enabled

The EMBoxerModule is basically the epicenter of everything: functions like "add_box" and "move_box" are probably good starting
points in terms of figuring out how to adapt this code to application specific needs
'''

from emapplication import EMStandAloneApplication,get_application
from pyemtbx.boxertools import BigImageCache,BinaryCircleImageCache,Cache
from EMAN2 import file_exists,EMANVERSION,gimme_image_dimensions2D,EMData,get_file_tag,get_image_directory,Region,file_exists,gimme_image_dimensions3D,abs_path
from EMAN2db import db_open_dict,db_check_dict,db_close_dict

import os,sys,weakref,math
from optparse import OptionParser


TEMPLATE_MIN = 30

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <image> <image2>....

This is a a really basic version of boxer that can be copied 
and used as the basis of developing something more advanced.
Design is meant to be granular and logical, i.e. easy to add to.

To see how it works just run it from the commandline:
./e2boxerbase.py test.mrc --boxsize=128 # single image handling
-or-
./e2boxerbase.py test.mrc test2.mrc --boxsize=128 # multiple images

The EMBoxerModule class is the epicenter of everything: 
functions like "add_box" and "move_box" are probably good starting
points in terms of figuring out how to adapt this code to your
specific needs/algorithms
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
	#this is an example of how to add your own custom tools:
	#module.add_2d_window_mouse_tool(EraseEventHandling,ErasingPanel,erase_radius=2*options.boxsize)
	application.execute()
	
def check(options,args):
	error_message = []
	if options.boxsize < 0:
		error_message.append("Box size muse be greater than 0")
		
	if options.boxsize > 4096:
		error_message.append("Box size must be less than 4096")
		
	if len(args) == 0:
		error_message.append("Please specify files")
	else:
		for file in args:
			if not file_exists(file):
				error_message.append("%s does not exist" %file)
			else:
				nx,ny,nz = gimme_image_dimensions3D(file)
				if nz > 1 or ny == 1 or nx == 1:
					error_message.append("%s is not 2D" %file)
					
	return error_message

def get_database_entry(image_name,key,database="bdb:e2boxercache",dfl=None):
	if not db_check_dict(database+"#"+key) and dfl==None:  return None
	
	db = db_open_dict(database+"#"+key)
	
	if db.has_key(image_name):
		return db[image_name]
	elif dfl != None:
		db[image_name] = dfl
		return dfl
	else: return None
	
def set_database_entry(image_name,key,value,database="bdb:e2boxercache"):
	'''
	write a key/object pair to the Image Database Dictionary associat
	'''
	db = db_open_dict(database+"#"+key)
	db[image_name] = value
	
def set_idd_image_entry(image_name,key,image,db_title="bdb:e2boxercache#"):
	'''
	Using EMAN2 style image dbs has efficiency payoffs in various ways... 
	'''
	image.set_attr("src_image",image_name)
	db_name =  db_title+key
	# first have to make sure it's not already there
	db = db_open_dict(db_name)
	
	if db.has_key("maxrec"):
		for idx in range(0,db["maxrec"]+1):
			header = db.get_header(idx)
			try:
				if header["src_image"] == image_name:
					db[idx] = image
					#image.write_image("bdb:e2boxercache#"+key,idx)
					db_close_dict(db_name)
					return
			except:
				pass

	# if we make it here then go ahead and append
	image.write_image(db_title+key,-1)
	db_close_dict(db_name)

def get_idd_image_entry(image_name,key,db_title="bdb:e2boxercache#"):
	'''
	Using EMAN2 style image dbs has efficiency payoffs in various ways... 
	'''
	db_name =  db_title+key
	if not db_check_dict(db_name): 
#		print "db failed",db_name
		return None # t
	
	db = db_open_dict(db_name)
	
	if db.has_key("maxrec"):
		for idx in range(0,db["maxrec"]+1):
			header = db.get_header(idx)
			try:
				if header["src_image"] == image_name:
					e = db[idx]
					db_close_dict(db_name)
					return e
			except:
				pass

	db_close_dict(db_name)
	return None
	

class BoxerThumbsWindowEventHandling:
	'''
	
	'''
	def __init__(self,target,thumbs_window):
		self.target = weakref.ref(target)
		self.thumbs_window = weakref.ref(thumbs_window)
		
		self.__connect_signals_to_slots()
		
	def __connect_signals_to_slots(self):
		'''
		connects the signals of the main 2D window to the slots of this object
		'''
		from PyQt4 import QtCore
		QtCore.QObject.connect(self.thumbs_window().emitter(),QtCore.SIGNAL("mx_mouseup"),self.thumb_image_selected)
		QtCore.QObject.connect(self.thumbs_window().emitter(),QtCore.SIGNAL("module_closed"),self.module_closed)
		
	def thumb_image_selected(self,event,lc):
		im=lc[0]
		self.target().set_current_file_by_idx(im)
		self.target().get_2d_window().updateGL()
	
	def module_closed(self):
		self.target().thumbs_window_closed()

class ScaledExclusionImage:
	database_name = "boxer_exclusion_image" # named it this to avoid conflicting with ExclusionImage
	def __init__(self,image_name):
		self.image_name = image_name
		self.image = None
		
		try:
			self.image = get_idd_image_entry(self.image_name,ScaledExclusionImage.database_name)
		except:
			pass
		
	def get_image_name(self):
		return self.image_name
		
	def get_shrink(self):
		if self.image != None: return self.image.get_attr("shrink")
		else: return 0
		
	def __update_image(self,shrink):
		'''
		Updates the image using the function arguments
		If they match current parameters than nothing happens - the correct image is already cached
		'''
		nx,ny,nz = gimme_image_dimensions3D(self.image_name)
		xsize = int(nx/shrink)
		ysize = int(ny/shrink)
		if self.image == None:
			self.image = EMData(xsize,ysize)
			self.image.to_zero()
			self.image.set_attr("shrink",int(shrink))
		else:
			# if the image already exists then we must retain the information in it by scaling and resizing it
			old_shrink = self.get_shrink()
			oldxsize = int(nx/old_shrink)
			oldysize = int(ny/old_shrink)
			r = Region( (oldxsize-xsize)/2, (oldysize-ysize)/2,xsize,ysize )
			#print "clipping to",(oldxsize-xsize)/2, (oldysize-ysize)/2,xsize,ysize
			scale = float(xsize)/float(oldxsize)
			
			# the order in which you clip and scale is dependent on whether or not scale is > 1
			if scale > 1:
				# if it's greater than one than clip (enlargen the image) first
				self.image.clip_inplace(r)
				# then scale the pixels
				self.image.scale(float(xsize)/float(oldxsize))
			else:
				# if it's less than one scale first so that we retain the maximum amount of the pixel information
				self.image.scale(float(xsize)/float(oldxsize))
				self.image.clip_inplace(r)
			
			self.image.set_attr("shrink",int(shrink))
			#set_idd_key_entry(self.image_name,"exclusion_image",self.image) # not sure if this is necessary
			
		#else:
			#print "doing nothing to currently stored small image in CoarsenedFlattenedImage"
			
	def get_image(self):
		'''
		Should only be called if you know the stored image is up to date
		'''
		return self.image
	
	
	def get_image_carefully(self,shrink):
		
		if self.image == None or not (shrink == self.get_shrink()):
			self.__update_image(shrink)
		
		return self.get_image()

ScaledExclusionImageCache = Cache(ScaledExclusionImage)

class ErasingPanel:
	def __init__(self,target,erase_radius=128):
		self.busy = True
		self.erase_radius = erase_radius
		self.target = weakref.ref(target)
		self.erase_rad_edit = None
		self.widget = None
		self.busy = False
		
	def icon(self):
		from PyQt4 import QtGui
		return QtGui.QIcon(get_image_directory() + "boxer_erase.png")
		
	def set_erase_radius(self,val):
		self.busy=True
		self.erase_radius = val
		if self.erase_rad_edit != None: self.erase_rad_edit.setText(str(val))
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
			
			self.erase_rad_edit = QtGui.QLineEdit(str(self.erase_radius))
			hbl.addWidget(self.erase_rad_edit)
			
			self.unerase = QtGui.QCheckBox("Unerase")
			self.unerase.setChecked(False)
			
			vbl.addLayout(hbl)
			vbl.addWidget(self.unerase)
			QtCore.QObject.connect(self.erase_rad_edit,QtCore.SIGNAL("editingFinished()"),self.new_erase_radius)
			QtCore.QObject.connect(self.unerase,QtCore.SIGNAL("clicked(bool)"),self.unerase_checked)
			
		return self.widget
	
	def new_erase_radius(self):
		if self.busy: return
		self.target().set_erase_radius(float(self.erase_rad_edit.text()))
		
	def unerase_checked(self,val):
		if self.busy: return
		self.target().toggle_unerase(val)
	def hide(self):
		if self.widget != None:
			self.widget.hide()
			
class ManualBoxingPanel:
	def __init__(self,target):
		self.target = weakref.ref(target)
		self.widget = None
		
	def icon(self):
		from PyQt4 import QtGui
		return QtGui.QIcon(get_image_directory() + "white_box.png")

	def get_widget(self):
		if self.widget == None:
			from PyQt4 import QtCore, QtGui, Qt
			self.widget = QtGui.QWidget()
			vbl = QtGui.QVBoxLayout(self.widget)
			vbl.setMargin(0)
			vbl.setSpacing(6)
			vbl.setObjectName("vbl")
			self.clear=QtGui.QPushButton("Clear")
			vbl.addWidget(self.clear)
			
			QtCore.QObject.connect(self.clear, QtCore.SIGNAL("clicked(bool)"), self.clear_clicked)
		return self.widget
	def hide(self):
		if self.widget != None:
			self.widget.hide()

	def clear_clicked(self,val):
		self.target().clear_all()

class EraseEventHandling:
	'''
	A class that knows how to handle mouse erase events for a GUIBox
	'''
	
	def __init__(self,target,panel_object=None,erase_radius=128):
		self.target = weakref.ref(target)
		self.panel_object = panel_object
		self.erase_value = 0.1			# erase mode can be either Boxable.ERASE or Boxable.UNERASE
		self.erase_radius = erase_radius
	
	def unique_name(self): return "Erase"
	
	def set_panel_object(self,panel): self.panel_object = panel
		
	def set_erase_radius(self,val): self.erase_radius = val
	
	def set_current_file(self,file_name):
		'''
		If the behavior of this Handler does not if the file changes, but the function needs to be supplied 
		'''
		pass
	
	def toggle_unerase(self,val):
		if val: self.erase_value = 0.0
		else: self.erase_value = 0.1
	
	def get_2d_window(self): return self.target().get_2d_window()
		
	def mouse_move(self,event):
		from emshape import EMShape
		m = self.get_2d_window().scr_to_img((event.x(),event.y()))
		self.get_2d_window().add_eraser_shape("eraser",EMShape(["circle",.1,.1,.1,m[0],m[1],self.erase_radius,3]))
		self.get_2d_window().updateGL()
		
	def adjust_erase_rad(self,delta):
		v = self.erase_radius
		if delta > 0:
			v = 1.1*v
		if delta < 0:
			v = 0.9*v
		self.erase_radius = v
		self.panel_object.set_erase_radius(v)
	
	def mouse_wheel(self,event):
		from PyQt4.QtCore import Qt
		if event.modifiers()&Qt.ShiftModifier:
			from emshape import EMShape
			self.adjust_erase_rad(event.delta())
			m= self.get_2d_window().scr_to_img((event.x(),event.y()))
			self.get_2d_window().add_eraser_shape("eraser",EMShape(["circle",.1,.1,.1,m[0],m[1],self.erase_radius,3]))
			self.get_2d_window().updateGL()
	
	def mouse_down(self,event) :
		from emshape import EMShape
		m=self.get_2d_window().scr_to_img((event.x(),event.y()))
		#self.boxable.add_exclusion_area("circle",m[0],m[1],self.erase_radius)
		self.get_2d_window().add_eraser_shape("eraser",EMShape(["circle",.9,.9,.9,m[0],m[1],self.erase_radius,3]))
		self.target().exclusion_area_added("circle",m[0],m[1],self.erase_radius,self.erase_value)	

	def mouse_drag(self,event) :
		from emshape import EMShape
		m=self.get_2d_window().scr_to_img((event.x(),event.y()))
		self.get_2d_window().add_eraser_shape("eraser",EMShape(["circle",.9,.9,.9,m[0],m[1],self.erase_radius,3]))
		self.target().exclusion_area_added("circle",m[0],m[1],self.erase_radius,self.erase_value)
		# exclusion_area_added does the OpenGL update calls, so there is no need to do so here
		
	def mouse_up(self,event) :
		# we have finished erasing
		
		# make the eraser shape non visible
		self.get_2d_window().add_eraser_shape("None",None)
		self.target().erasing_done(self.erase_value)
		
class ManualBoxingEventHandling:
	'''
	A class that knows how to add, move and remove reference and non reference boxes 
	'''
	def __init__(self,target,panel_object=None):
		self.target = weakref.ref(target)
		self.panel_object = panel_object
		self.moving = None
		
	def set_panel_object(self,panel): self.panel_object = panel
	def unique_name(self): return "Manual"
	
	def set_current_file(self,file_name):
		'''
		If the behavior of this Handler does not if the file changes, but the function needs to be supplied 
		'''
		pass
		
	def get_2d_window(self): return self.target().get_2d_window()
		
	def mouse_down(self,event) :
		m = self.get_2d_window().scr_to_img((event.x(),event.y()))
		box_num = self.target().detect_box_collision(m)
		from PyQt4.QtCore import Qt
		if box_num == -1:
			if event.modifiers()&Qt.ShiftModifier : return # the user tried to delete nothing
			self.target().add_box(m[0],m[1])
		
		elif event.modifiers()&Qt.ShiftModifier :
			# remove the box
			self.target().remove_box(box_num)
		else:
			# if we make it here than the we're moving a box
			self.moving=[m,box_num]
			self.target().moving_box_established(box_num)

	def mouse_drag(self,event) :
		
		m=self.get_2d_window().scr_to_img((event.x(),event.y()))
		from PyQt4.QtCore import Qt
		if event.modifiers()&Qt.ShiftModifier:
			box_num = self.target().detect_box_collision(m)
			if ( box_num != -1):
				self.target().remove_box(box_num)
			
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
		self.target().clear_boxes(["manual"])

class Main2DWindowEventHandling:
	'''
	A class that responds to added, moved and removed box signals emitted
	by the by the main 2d image in the EMBoxerModule - this is the image
	display that shows the image that is currently being boxed along
	with any boxed regions
	'''
	def __init__(self,target,main_2d_window):
		self.target = weakref.ref(target) # prevent a strong cycle
		self.main_2d_window = main_2d_window
		
		self.mouse_handlers = {}  #stores mouse handlers
		self.mouse_handler = None  # the current mouse handler
		
		self.__connect_signals_to_slots()
		
	def set_mouse_mode(self,name):
		try:
			self.mouse_handler = self.mouse_handlers[name]
		except:
			print "unknown handler",name
	
	def set_current_file(self,file_name):
		for value in self.mouse_handlers.values(): value.set_current_file(file_name)
			
	def add_mouse_handler(self,handler,name):
		self.mouse_handlers[name] = handler
		if self.mouse_handler == None: self.mouse_handler = handler
		
	def __connect_signals_to_slots(self):
		'''
		connects the signals of the main 2D window to the slots of this object
		'''
		from PyQt4 import QtCore
		QtCore.QObject.connect(self.main_2d_window.emitter(),QtCore.SIGNAL("mousedown"),self.mouse_down)
		QtCore.QObject.connect(self.main_2d_window.emitter(),QtCore.SIGNAL("mousedrag"),self.mouse_drag)
		QtCore.QObject.connect(self.main_2d_window.emitter(),QtCore.SIGNAL("mouseup")  ,self.mouse_up  )
		QtCore.QObject.connect(self.main_2d_window.emitter(),QtCore.SIGNAL("keypress"),self.key_press)
		QtCore.QObject.connect(self.main_2d_window.emitter(),QtCore.SIGNAL("mousewheel"),self.mouse_wheel)
		QtCore.QObject.connect(self.main_2d_window.emitter(),QtCore.SIGNAL("mousemove"),self.mouse_move)
		QtCore.QObject.connect(self.main_2d_window.emitter(),QtCore.SIGNAL("module_closed"),self.module_closed)
	
	def mouse_down(self,event) :
		'''
		@param a QtGui.QMouseEvent sent from the EMImage2DModule
		'''
		if self.mouse_handler == None: return
		self.mouse_handler.mouse_down(event)
	
	def mouse_drag(self,event) :
		'''
		@param a QtGui.QMouseEvent sent from the EMImage2DModule
		'''
		if self.mouse_handler == None: return
		self.mouse_handler.mouse_drag(event)
	
	def mouse_up(self,event) :
		'''
		@param a QtGui.QMouseEvent sent from the EMImage2DModule
		'''
		if self.mouse_handler == None: return
		self.mouse_handler.mouse_up(event)
		
	def key_press(self,event):
		'''
		@param a QtGui.QKeyEvent sent from the EMImage2DModule
		'''
		if self.mouse_handler == None: return
		print "2d window key press"
		
	def mouse_wheel(self,event):
		'''
		@param a QtGui.QMouseEvent sent from the EMImage2DModule
		'''
		if self.mouse_handler == None: return
		self.mouse_handler.mouse_wheel(event)
		
	def mouse_move(self,event):
		'''
		@param a QtGui.QMouseEvent sent from the EMImage2DModule
		'''
		if self.mouse_handler == None: return
		self.mouse_handler.mouse_move(event)
		
	def module_closed(self):
		'''
		'''
		self.target().main_2d_window_closed()
	
class ParticlesWindowEventHandling:
	def __init__(self,target,particle_window):
		self.target = weakref.ref(target) # prevent a strong cycle
		self.particle_window = particle_window
		
		self.__connect_signals_to_slots()
		self.moving_box_data = None
		self.first_clicked = None
	def __connect_signals_to_slots(self):
		'''
		connects the signals of the main 2D window to the slots of this object
		'''
		from PyQt4 import QtCore
		QtCore.QObject.connect(self.particle_window.emitter(),QtCore.SIGNAL("mx_image_selected"),self.box_selected)
		QtCore.QObject.connect(self.particle_window.emitter(),QtCore.SIGNAL("mx_mousedrag"),self.box_moved)
		QtCore.QObject.connect(self.particle_window.emitter(),QtCore.SIGNAL("mx_mouseup"),self.box_released)
		QtCore.QObject.connect(self.particle_window.emitter(),QtCore.SIGNAL("mx_boxdeleted"),self.box_image_deleted)
		QtCore.QObject.connect(self.particle_window.emitter(),QtCore.SIGNAL("module_closed"),self.module_closed)
	
	def box_selected(self,event,lc):
		im=lc[0]
		if im == None: return
		self.moving_box_data = [event.x(),event.y(),im]
		self.first_clicked = im
		self.target().moving_box_established(im)
		#self.target().get_2d_window().set_active(im,.9,.9,.4)
		self.target().get_2d_window().updateGL()
		
	def box_moved(self,event,scale):
		if self.moving_box_data:
			dx = (self.moving_box_data[0] - event.x())/scale
			dy = (event.y() - self.moving_box_data[1])/scale
			self.target().move_box(self.moving_box_data[2],dx,dy)
			
			self.moving_box_data = [event.x(),event.y(),self.moving_box_data[2]]

	def box_released(self,event,lc):
		if lc == None: return
		
		im = lc[0]
		if im == None: return
		if im == self.first_clicked: self.target().particle_selected(im)
		else: self.target().box_released(self.first_clicked)
		self.moving_box_data = None
		self.first_clicked = None
		
	def box_image_deleted(self,event,lc,force_image_mx_remove=True):
		if lc[0] == None: return
		box = self.target().remove_box(lc[0])
	
	def module_closed(self):
		'''
		'''
		self.target().particles_window_closed()
	


class EMThumbsTools:
	
	def gen_thumbs(image_names=[],shrink=None):
		'''
		@param image_names a list of image names
		@param shrink the shrink factor, should be an int 
		@return a list of very small images
		'''
		from emapplication import EMProgressDialogModule
		
		if shrink == None: shrink = EMThumbsTools.get_image_thumb_shrink(image_names[0])
		
		application = get_application()
		nim = len(image_names)
		thumbs = [None for i in range(nim)]
		progress = EMProgressDialogModule(application,"Generating Thumbnails", "Abort", 0, nim,None)
		progress.qt_widget.show()
		prog = 0
		for i in range(nim):
	#				thumb = self.get_image_thumb(i)
			thumb = get_idd_image_entry(image_names[i],"image_thumb")
			if thumb == None:
				thumb=EMData(image_names[i],0)
				thumb.process_inplace("math.meanshrink",{"n":shrink})
				thumb.process_inplace("normalize.edgemean") # if there are lots than they =should all have the same contrast
				thumb.set_attr("display_name",get_file_tag(image_names[i]))
				set_idd_image_entry(image_names[i],"image_thumb",thumb)

			prog += 1
			progress.qt_widget.setValue(prog)
			application.processEvents()
		
			if progress.qt_widget.wasCanceled():
				progress.qt_widget.setValue(nim)
				progress.qt_widget.close()
				return None
			
		progress.qt_widget.setValue(nim)
		progress.qt_widget.close()
		
		for i in range(nim):
			# warning : DO NOT move this into the loop above - for some reason if you store the thumbs in the list while they're being generated there is a HUGE memory leak
			# this solves it
			thumbs[i] = get_idd_image_entry(image_names[i],"image_thumb")
				
		
		return thumbs
	
	gen_thumbs = staticmethod(gen_thumbs)
	
	def get_image_thumb_shrink(image_name):
		shrink = 1
		inx,iny =  gimme_image_dimensions2D(image_name)
		inx /= 2
		iny /= 2
		while ( inx >= 128 and iny >= 128):
			inx /= 2
			iny /= 2
			shrink *= 2
	
		return shrink
	
	get_image_thumb_shrink = staticmethod(get_image_thumb_shrink)


class EMBox:
	def __init__(self,x,y,type="manual",score=0.0):
		self.x = x # central x coordinate
		self.y = y # central y coordinate
		self.type = type # type can be customized
		self.score = score # can be some kind of score, such as correlation
		self.image = None # an image
		
	def move(self,dx,dy):
		self.x += dx
		self.y += dy
		self.image = None
		
	def get_image(self,image_name,box_size):
		
		if self.image == None or self.image.get_xsize() != box_size or self.image.get_ysize() != box_size:
			global BigImageCache
			data=BigImageCache.get_object(image_name).get_image(use_alternate=True)
			r = Region(self.x-box_size/2,self.y-box_size/2,box_size,box_size)
			self.image = data.get_clip(r)
			
		return self.image
	
	def reset_image(self): self.image = None
		
	def get_shape(self,shape_string,box_size):
		if self.type == "manual":
			r,g,b = 1,1,1
		elif self.type == "ref":
			r,g,b = 0,0,0
		elif self.type == "auto":
			r,g,b = 0.4,0.9,0.4 # historical green, the original color
		else: raise RuntimeError("Unknown type")
		from emshape import EMShape
		shape = EMShape([shape_string,r,g,b,self.x-box_size/2,self.y-box_size/2,self.x+box_size/2,self.y+box_size/2,2.0])
		return shape
	
	def collision(self,x,y,box_size):
		
		if x-box_size/2 < self.x and x+box_size/2 > self.x and y-box_size/2 < self.y and y+box_size/2 > self.y: return True
		else: return False
	

class EMBoxList:
	'''
	A list of boxes
	'''
	CACHE_SIZE = 10
	def __init__(self,target=None):
		if target: self.target=weakref.ref(target)
		else: self.target = None # someone is using this incognito
		self.current_iter = 0 # iterator support
		self.max_idx = 0 # iterator support
		self.boxes = []
		self.shapes = []
		self.shape_string = "rectpoint" # for making shapes
		
		self.undo_cache = []
		self.redo_cache = []	
	
	def get_boxes(self): return self.boxes
	
	def redo(self):
		if len(self.redo_cache) == 0: 
			print "can't redo, nothing stored"
			return False
		
		self.undo_cache.append([[box.x,box.y,box.type] for box in self.boxes])
		cache = self.redo_cache[-1]
		self.boxes = [EMBox(x,y,type) for x,y,type in cache]
		self.redo_cache.pop(-1)
		self.reset_shapes()
		return True
	def is_redoable(self): return len(self.redo_cache) > 0
	
	def undo(self):
		if len(self.undo_cache) == 0: 
			print "can't undo, nothing stored"
			return False
		
		self.redo_cache.append([[box.x,box.y,box.type] for box in self.boxes])
		cache = self.undo_cache[-1]
		self.boxes = [EMBox(x,y,type) for x,y,type in cache]
		self.undo_cache.pop(-1)
		self.reset_shapes()
		return True
	
	def is_undoable(self): return len(self.undo_cache) > 0
	
	def cache_old_boxes(self):
		if len(self.undo_cache) >= EMBoxList.CACHE_SIZE: self.undo_cache.pop(0)
		self.undo_cache.append([[box.x,box.y,box.type] for box in self.boxes])
		self.redo_cache = []
	
	def clear_boxes(self,types):
		self.cache_old_boxes()
		
		for i in xrange(len(self.boxes)-1,-1,-1):
			if self.boxes[i].type in types:
				self.boxes.pop(i)
				self.shapes.pop(i)
						
	def add_box(self,x,y,type="manual",score=0.0):
		self.cache_old_boxes()
		
		self.boxes.append(EMBox(x,y,type,score))
		self.shapes.append(None)
		
	def add_boxes(self,boxes):
		'''
		boxes should be a list like [[x,y,type],[x,y,type],....[int,int,string]]
		'''
		
		self.cache_old_boxes()
		
		for box in boxes:
			self.boxes.append(EMBox(*box))
			self.shapes.append(None)
	
	def reset_images(self):
		for box in self.boxes: box.reset_image()
		
	def reset_shapes(self):
		self.shapes = [None for i in range(len(self.boxes))]
	
	def move_box(self,i,dx,dy):
		self.boxes[i].move(dx,dy)
		self.shapes[i] = None
		
	def get_shape(self,i,box_size):
		if self.shapes[i] == None:
			shape = self.boxes[i].get_shape(self.shape_string,box_size)
			self.shapes[i] = shape
			return shape
		return self.shapes[i]
	
	def get_box_type(self,box_number):
		'''
		@param box_number the number of the box for which you want to get the type i.e. that which was returned from the detect_box_collision
		'''
		return self.boxes[box_number].type
	
	def detect_collision(self,x,y,box_size):
		for i,box in enumerate(self.boxes):
			if box.collision(x,y,box_size):
				return i
		
		return -1
		
	def get_box_shapes(self,box_size):
		
		d = {}
		for i in range(len(self.boxes)):
			if self.shapes[i] == None:
				shape = self.boxes[i].get_shape(self.shape_string,box_size)
				self.shapes[i] = shape
			
			d[i] = self.shapes[i]
			
		return d

	def get_particle_images(self,image_name,box_size):
		return [box.get_image(image_name,box_size) for box in self.boxes]
	
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
		self.shapes.pop(idx)
		return self.boxes.pop(idx)
	
	def remove_box(self,idx):
		self.cache_old_boxes()
		self.shapes.pop(idx)
		return self.boxes.pop(idx)
		
	
	def get_boxes_for_database(self):
		return [[box.x,box.y,box.type] for box in self.boxes]
			
	def save_boxes_to_database(self,image_name):

		set_database_entry(image_name,"boxes",self.get_boxes_for_database())
			
		set_database_entry(image_name,"undo_cache",self.undo_cache)
		set_database_entry(image_name,"redo_cache",self.redo_cache)

	def load_boxes_from_database(self,image_name,reset=True,load_caches=True):
		if reset:
			self.boxes = []
			self.shapes = []
			self.undo_cache = []
		
		data = get_database_entry(image_name,"boxes")
		if data != None:
			for x,y,type in data:
				self.add_box(x,y,type)
		
		if load_caches:
			self.undo_cache = get_database_entry(image_name,"undo_cache",dfl=[])
			self.redo_cache = get_database_entry(image_name,"redo_cache",dfl=[])

	def exclude_from_scaled_image(self,exclusion_image,subsample_rate):
		action = False
		for i in xrange(len(self.boxes)-1,-1,-1):
			box = self.boxes[i]
			x = int(box.x/subsample_rate)
			y = int(box.y/subsample_rate)
			if exclusion_image.get(x,y):
				if action == False: self.cache_old_boxes()
				self.pop(i)
				action = True
			
		return action
	
	def write_particles(self,input_file_name,out_file_name,box_size,invert=False,normproc=None):
		print invert,normproc
		for i,box in enumerate(self.boxes):
			image = box.get_image(input_file_name,box_size)
			if invert: image.mult(-1)
			if str(normproc) != "None": image.process_inplace(normproc)
			image.write_image(out_file_name,i)
			
	
	def write_coordinates(self,input_file_name,out_file_name,box_size):
		f=file(out_file_name,'w')
		for box in self.boxes:
			xc = box.x-box_size/2
			yc = box.y-box_size/2
			f.write(str(int(xc))+'\t'+str(int(yc))+'\t'+str(box_size)+'\t'+str(box_size)+'\n')
		f.close()

class EMBoxerModule:
	'''
	This module is essentially a Mediator (see Design Patterns) - it coordinates the activities of several EMAN2 modules
	that would otherwise not necessary interact.
	'''
	def __init__(self,file_names=[],box_size=128):
		'''
		@file_name the name of a file on disk
		@exception RuntimeError raised if the file does not exist
		'''
		self.file_names = file_names # a list of file names
		self.current_idx = None # an index into self.file_names
		self.box_size = box_size # the current box size
		
		self.signal_slot_handlers = {} # this is a dictionary, keys are (somewhat random) names, values are event handlers such as Main2DWindowEventHandling. This dict has the only reference to the event handlers
		self.inspector = None # this will be a Qt style inspector
		self.inspector_module = None # the wrapping object of self.inspector
		self.main_2d_window = None # this will be the main 2D image display, showing boxed particles etc 
		self.particles_window = None # this will be the window displaying the picked particles
		self.thumbs_window = None # this will be the window showing the thumbnails, enabling the user to change between 2D raw data
		self.image_thumbs = None # image_thumbs is a list of thumbnail images	
		self.box_list = EMBoxList(self)
		self.moving_box = None
		self.output_task = None # will be an EMAN2 style form for writing output
		# initialized the 2D window
		self.__init_main_2d_window()
		if len(self.file_names) > 1: self.__init_thumbs_window()
		
		# initialize the inspector
		self.__init_inspector()
		
		# this is an example of how to add your own custom tools:
		self.add_2d_window_mouse_tool(ManualBoxingEventHandling,ManualBoxingPanel)
		self.add_2d_window_mouse_tool(EraseEventHandling,ErasingPanel,erase_radius=2*box_size)

	
	def redo_boxes(self):
		if not self.box_list.redo(): return
		else: self.full_box_update()
		
		if self.inspector: 
			self.inspector.enable_undo_redo(self.box_list.is_undoable(),self.box_list.is_redoable())
	
	def undo_boxes(self):
		if not self.box_list.undo(): return
		else: self.full_box_update()
		if self.inspector: 
			self.inspector.enable_undo_redo(self.box_list.is_undoable(),self.box_list.is_redoable())
	
	def moving_box_established(self,box_number):
		if self.moving_box != None: raise RuntimeError("establishing a moving box that was already established?")
		self.moving_box = box_number
	
	def clear_boxes(self,type):
		self.box_list.clear_boxes(type)
		if self.particles_window:
			self.particles_window.set_data(self.box_list.get_particle_images(self.current_file(), self.box_size))
			self.particles_window.updateGL()
		if self.main_2d_window:
			self.main_2d_window.set_shapes(self.box_list.get_box_shapes(self.box_size))
			self.main_2d_window.updateGL()
			
		if self.inspector: 
			self.inspector.enable_undo_redo(self.box_list.is_undoable(),self.box_list.is_redoable())
	
	def has_thumbs(self):
		return self.image_thumbs != None
	
	def get_subsample_rate(self): 
		'''
		
		'''
		return int(math.ceil(float(self.box_size)/float(TEMPLATE_MIN)))

	
	def get_box_type(self,box_number):
		'''
		@param box_number the number of the box for which you want to get the type i.e. that which was returned from the detect_box_collision
		'''
		return self.box_list.get_box_type(box_number)
	
	def detect_box_collision(self,data):
		return self.box_list.detect_collision(data[0], data[1], self.box_size)
	
	def particle_selected(self,box_number):
		if self.moving_box != None: self.moving_box = None
		box= self.box_list[box_number]
		self.box_placement_update_exclusion_image(box.x,box.y)
		if self.main_2d_window: self.main_2d_window.register_scroll_motion(box.x,box.y)
	
	def box_released(self,box_number):
		if self.moving_box != None: self.moving_box = None
		box= self.box_list[box_number]
		self.box_placement_update_exclusion_image(box.x,box.y)
	
	def add_boxes(self,boxes):
		'''
		boxes should be a list like [[x,y,type],[x,y,type],....[int,int,string]]
		'''
		if get_database_entry(self.current_file(),"frozen",dfl=False): return
		if self.particles_window == None:
			self.__init_particles_window()
			get_application().show_specific(self.particles_window)
			
		self.box_list.add_boxes(boxes)
		self.box_list.save_boxes_to_database(self.current_file())
		if self.particles_window:
			self.particles_window.set_data(self.box_list.get_particle_images(self.current_file(), self.box_size))
			self.particles_window.updateGL()
		if self.main_2d_window:
			self.main_2d_window.update_shapes(self.box_list.get_box_shapes(self.box_size))
			self.main_2d_window.updateGL()
		
		if self.inspector: 
			self.inspector.enable_undo_redo(self.box_list.is_undoable(),self.box_list.is_redoable())
		
	
	def add_box(self,x,y,type="manual"):
		if get_database_entry(self.current_file(),"frozen",dfl=False): return
		if self.particles_window == None:
			self.__init_particles_window()
			get_application().show_specific(self.particles_window)

		self.box_placement_update_exclusion_image(x,y)
		self.box_list.add_box(x,y,type=type)
		self.box_list.save_boxes_to_database(self.current_file())
		if self.particles_window:
			self.particles_window.set_data(self.box_list.get_particle_images(self.current_file(), self.box_size))
			self.particles_window.updateGL()
		if self.main_2d_window:
			self.main_2d_window.update_shapes(self.box_list.get_box_shapes(self.box_size))
			self.main_2d_window.updateGL()
		
		if self.inspector: 
			self.inspector.enable_undo_redo(self.box_list.is_undoable(),self.box_list.is_redoable())
		
	def box_placement_update_exclusion_image(self,x,y):
		exclusion_image = self.get_exclusion_image()
		if exclusion_image != None:
			sr = self.get_subsample_rate()
			xx,yy = int(x/sr),int(y/sr)
			if exclusion_image.get(xx,yy):
				from EMAN2 import BoxingTools
				global BinaryCircleImageCache
				mask = BinaryCircleImageCache.get_image_directly(int(self.box_size/(2*sr)))
				BoxingTools.set_region(self.get_exclusion_image(),mask,xx,yy,0.0)
				set_idd_image_entry(self.current_file(),ScaledExclusionImage.database_name,self.get_exclusion_image())
				if self.main_2d_window:
					self.main_2d_window.set_other_data(self.get_exclusion_image(),self.get_subsample_rate(),True)
					self.main_2d_window.updateGL()
		
	def remove_box(self,box_number):
		if get_database_entry(self.current_file(),"frozen",dfl=False): return
		self.box_list.remove_box(box_number)
		self.box_list.save_boxes_to_database(self.current_file())
		self.full_box_update()
		if self.inspector: 
			self.inspector.enable_undo_redo(self.box_list.is_undoable(),self.box_list.is_redoable())
		
	def full_box_update(self):
		if self.particles_window != None:
			self.particles_window.set_data(self.box_list.get_particle_images(self.current_file(), self.box_size))
			self.particles_window.updateGL()
		if self.main_2d_window != None:
			self.main_2d_window.set_shapes(self.box_list.get_box_shapes(self.box_size))
			self.main_2d_window.updateGL()
	
	def move_box(self,box_number,dx,dy):
		if self.moving_box != None:
			self.box_list.cache_old_boxes()
			self.moving_box = None
		self.box_list.move_box(box_number,dx,dy)
		self.box_list.save_boxes_to_database(self.current_file())
		
		if self.particles_window:
			self.particles_window.set_data(self.box_list.get_particle_images(self.current_file(), self.box_size))
			self.particles_window.updateGL()
		if self.main_2d_window:
			self.main_2d_window.add_shape(box_number,self.box_list.get_shape(box_number,self.box_size))
			self.main_2d_window.updateGL()
			
		if self.inspector: 
			self.inspector.enable_undo_redo(self.box_list.is_undoable(),self.box_list.is_redoable())
		
	def get_exclusion_image(self,mark_boxes=False):
		'''
		@mark_boxes if true the exclusion image is copied and the locations of the current boxes are painted in as excluded regions
		This is useful for autoboxers - they  obviously dont want to box any region that already has a box in it (such as a manual box,
		or a previously autoboxed box)
		'''
		exc_image = ScaledExclusionImageCache.get_image(self.current_file(),self.get_subsample_rate())
		if not mark_boxes: return exc_image
		else:
			image = exc_image.copy()
			boxes = self.box_list.get_boxes()
			if len(boxes) > 0:
				sr = self.get_subsample_rate()
				global BinaryCircleImageCache
				mask = BinaryCircleImageCache.get_image_directly(int(self.box_size/(2*sr)))
				for box in self.box_list.get_boxes():
					x,y = int(box.x/sr),int(box.y/sr)
					from EMAN2 import BoxingTools
					BoxingTools.set_region(image,mask,x,y,0.1) # 0.1 is also the value set by the eraser - all that matters is that it's zon_zero
			
			return image
	
	def exclusion_area_added(self,typeofexclusion,x,y,radius,val):
		xx = int(x/self.get_subsample_rate())
		yy = int(y/self.get_subsample_rate())
		
		rr = int(radius/self.get_subsample_rate())
		global BinaryCircleImageCache
		mask = BinaryCircleImageCache.get_image_directly(rr)

		from EMAN2 import BoxingTools
		BoxingTools.set_region(self.get_exclusion_image(),mask,xx,yy,val)
		
		if self.main_2d_window:
			self.main_2d_window.set_other_data(self.get_exclusion_image(),self.get_subsample_rate(),True)
			self.main_2d_window.updateGL()
	
	def erasing_done(self,erase_mode):
		set_idd_image_entry(self.current_file(),ScaledExclusionImage.database_name,self.get_exclusion_image())
		act = self.box_list.exclude_from_scaled_image(self.get_exclusion_image(),self.get_subsample_rate())
		if act:
			self.box_list.save_boxes_to_database(self.current_file())
			self.full_box_update()
			
	def __init_main_2d_window(self):
		from emimage2d import EMImage2DModule
		
		if self.main_2d_window == None:
			
			self.main_2d_window= EMImage2DModule(application=get_application())
	
			self.main_2d_window.set_mouse_mode(0)
					
			self.signal_slot_handlers["2d_window"] = Main2DWindowEventHandling(self,self.main_2d_window)
			
			get_application().show_specific(self.main_2d_window)
			
	def get_2d_window(self): return self.main_2d_window
	
	def main_2d_window_closed(self):
		self.main_2d_window = None
		if self.inspector:
			self.inspector.set_2d_window_visible(False)
			
	def add_2d_window_mouse_tool(self,events_manager,panel,**kargs):
		em = events_manager(self,**kargs)
		panel = panel(em,**kargs)
		em.set_panel_object(panel)
		if self.current_idx != None: em.set_current_file(self.current_file())
		self.signal_slot_handlers["2d_window"].add_mouse_handler(em,em.unique_name())
		self.inspector.add_mouse_tool(panel,em.unique_name())
		
	
	def set_current_file_by_idx(self,idx):
		if len(self.file_names) <= idx: raise RuntimeError("The index is beyond the length of the file names list")
		
		if idx != self.current_idx:
			self.current_idx = idx
			self.set_current_file(self.file_names[idx])
	
	def set_current_file(self,file_name):
		from PyQt4 import QtCore
		get_application().setOverrideCursor(QtCore.Qt.BusyCursor)
		
		if not file_exists(file_name): raise RuntimeError("The file %s does not exist" %file_name)
		
		self.signal_slot_handlers["2d_window"].set_current_file(file_name)
		
		if self.main_2d_window != None:

	   	   	self.load_2d_window_display(file_name)
			if self.inspector != None: 
				self.inspector.set_frozen(get_database_entry(file_name,"frozen",dfl=False))
				self.inspector.set_image_quality(get_database_entry(file_name,"quality",dfl=2))
						
			# the boxes should be loaded from the database, if possible
			self.box_list.load_boxes_from_database(file_name)
			if self.inspector: 
				self.inspector.enable_undo_redo(self.box_list.is_undoable(),self.box_list.is_redoable())
			self.main_2d_window.set_shapes(self.box_list.get_box_shapes(self.box_size))
			self.main_2d_window.updateGL()
			
			particles = self.box_list.get_particle_images(self.current_file(),self.box_size)
			if len(particles) > 0:
				if self.particles_window == None: self.__init_particles_window()
				self.particles_window.set_data(particles)
				self.particles_window.updateGL()
			
			
		get_application().setOverrideCursor(QtCore.Qt.ArrowCursor)
	
	def load_2d_window_display(self,file_name):
		global BigImageCache
		data=BigImageCache.get_object(file_name).get_image(use_alternate=True)
		
		if get_idd_image_entry(file_name,ScaledExclusionImage.database_name) != None:
			self.main_2d_window.set_other_data(self.get_exclusion_image(),self.get_subsample_rate(),True)
		else:
			self.main_2d_window.set_other_data(None,self.get_subsample_rate(),True)

		self.main_2d_window.set_data(data,file_name)
		
		frozen = get_database_entry(file_name,"frozen",dfl=False)
		if frozen == None:
			set_database_entry(file_name,"frozen",False)
			frozen = False 
		self.main_2d_window.set_frozen(frozen)
		
	def set_image_quality(self,val):
		set_database_entry(self.current_file(),"quality",val)
		
	def set_frozen(self,val):
		set_database_entry(self.current_file(),"frozen",val)
		self.main_2d_window.set_frozen(val)
		self.main_2d_window.updateGL()
		
	def current_file(self):
		return self.file_names[self.current_idx]
			
	def __init_thumbs_window(self,redo_thumbs=False):
		if len(self.file_names) == 0: raise RuntimeError("Will not make a thumbs window if the number of images is zero")
		
		if self.thumbs_window == None:
			from PyQt4 import QtCore
			get_application().setOverrideCursor(QtCore.Qt.BusyCursor)
			
			
			if self.image_thumbs == None or redo_thumbs:
				self.image_thumbs = EMThumbsTools.gen_thumbs(self.file_names)
			if self.image_thumbs == None:
				sys.exit(1)
			
			from emimagemx import EMImageMXModule
			self.thumbs_window=EMImageMXModule(application=get_application())
			self.thumbs_window.desktop_hint = "rotor" # this is to make it work in the desktop
				
			self.thumbs_window.set_data(self.image_thumbs,soft_delete=True)
			self.thumbs_window.set_mouse_mode("app")
			self.thumbs_window.update_window_title("Thumbnails")
			self.signal_slot_handlers["thumbs_window"] = BoxerThumbsWindowEventHandling(self,self.thumbs_window)
			
			get_application().setOverrideCursor(QtCore.Qt.ArrowCursor)
	
	def thumbs_window_closed(self):
		self.thumbs_window = None
		if self.inspector:
			self.inspector.set_thumbs_visible(False)
	
	def __init_inspector(self):
		if self.inspector == None:
			self.inspector_module = EMBoxerInspectorModule(self)
			self.inspector = self.inspector_module.widget
			self.inspector.set_box_size(self.box_size)
			
			try:
				#self.current_file() might fail
				frozen = get_database_entry(self.current_file(),"frozen",dfl=False)
				self.inspector.set_frozen(frozen)
			except: pass # inspector sets the frozen button to false by default
			
	def get_inspector(self): return self.inspector
	
	def __init_particles_window(self):
		if self.particles_window == None:
			from emimagemx import EMImageMXModule
			self.particles_window=EMImageMXModule(application=get_application())
			self.particles_window.desktop_hint = "rotor" # this is to make it work in the desktop
				
			self.particles_window.set_mouse_mode("app")
			self.particles_window.update_window_title("Particles")
			self.signal_slot_handlers["particles_window"] = ParticlesWindowEventHandling(self,self.particles_window)
			
	
	def particles_window_closed(self):
		self.particles_window = None
		if self.inspector:
			self.inspector.set_particles_visible(False)
		
	def show_thumbs_window(self,bool):
		print self.thumbs_window
		if self.thumbs_window == None: 
			self.__init_thumbs_window()
			
		print self.thumbs_window,"now"
		if bool:
			get_application().show_specific(self.thumbs_window)
		else:
			get_application().hide_specific(self.thumbs_window)
		
	def show_2d_window(self,bool):
		resize = False
		if self.main_2d_window == None:
			resize = True 
			self.__init_main_2d_window()
			self.load_2d_window_display(self.current_file())
			self.main_2d_window.set_shapes(self.box_list.get_box_shapes(self.box_size))
		if bool:
			get_application().show_specific(self.main_2d_window)
		else:
			get_application().hide_specific(self.main_2d_window)
		
		if resize:
			self.main_2d_window.optimally_resize()
			
	def show_particles_window(self,bool):
		resize = False
		if self.particles_window == None:
			resize = True 
			self.__init_particles_window()
			self.particles_window.set_data(self.box_list.get_particle_images(self.current_file(), self.box_size))
		if bool:
			get_application().show_specific(self.particles_window)
		else:
			get_application().hide_specific(self.particles_window)
			
		if resize: self.particles_window.optimally_resize()
	def show_interfaces(self):
		if len(self.file_names) > 0:	self.set_current_file_by_idx(0)
		
		if self.main_2d_window != None:
			get_application().show_specific(self.main_2d_window)
			self.main_2d_window.optimally_resize()
		if self.thumbs_window != None: 
			get_application().show_specific(self.thumbs_window)
			self.thumbs_window.optimally_resize()
		if self.inspector_module != None: 
			get_application().show_specific(self.inspector_module)
		if self.particles_window != None: 
			get_application().show_specific(self.particles_window)
			self.particles_window.optimally_resize()
		
	def get_box_size(self): return self.box_size
	def set_box_size(self,box_size):
		self.box_size = box_size
		self.box_list.reset_images()
		self.box_list.reset_shapes()
		self.full_box_update()
		
		
	def set_main_2d_mouse_mode(self,mode):
		if self.main_2d_window != None:
			self.main_2d_window.add_eraser_shape("None",None)
			self.main_2d_window.updateGL()
		self.signal_slot_handlers["2d_window"].set_mouse_mode(mode)
		
	def done(self):
		for module in [self.main_2d_window, self.thumbs_window,self.particles_window ]:
			if module != None: module.closeEvent(None)
				
	def run_output_dialog(self):
		from emsprworkflow import E2BoxerProgramOutputTask
		if self.output_task != None: return
		from PyQt4 import QtCore
		self.output_task = EMBoxerWriteOutputTask(self.file_names)
		QtCore.QObject.connect(self.output_task.emitter(),QtCore.SIGNAL("task_idle"),self.on_output_task_idle)
		self.output_task.run_form()
		
	def on_output_task_idle(self):
		self.output_task = None

from emsprworkflow import WorkFlowTask,error
class EMBoxerWriteOutputTask(WorkFlowTask):	
	"""Use this form for writing boxed particles and/or coordinate files to disk."""
	def __init__(self,file_names=[],output_formats=["hdf","spi","img","bdb"]):
		WorkFlowTask.__init__(self)
		self.window_title = "Write Particle Output"
		self.form_db_name = "bdb:e2boxerbase"
		self.file_names = file_names
		self.output_formats = output_formats
		
	def get_table(self):
		from emform import EM2DFileTable,EMFileTable,int_lt
		table = EM2DFileTable(self.file_names,desc_short="Raw Data",desc_long="")
		table.add_column_data(EMFileTable.EMColumnData("Stored Boxes",EMBoxerWriteOutputTask.get_num_boxes,"The number of stored boxes",int_lt))
		table.add_column_data(EMFileTable.EMColumnData("Quality",EMBoxerWriteOutputTask.get_quality,"Quality metadata score stored in local database",int_lt))
	
		return table
	
	def get_quality(file_name):
		'''
		A static function for getting the number of boxes associated with each file
		'''
		val = get_database_entry(file_name,"quality")
		
		if val == None: return "-"
		else: return str(val)
		
	
	
	def get_num_boxes(file_name):
		'''
		A static function for getting the number of boxes associated with each file
		'''
		box_list = EMBoxList()
		box_list.load_boxes_from_database(file_name,load_caches=False)
		return str(len(box_list))
	
	get_quality = staticmethod(get_quality)
	get_num_boxes = staticmethod(get_num_boxes)
	
	def get_params(self):
#		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2CTFGenericTask.documentation_string,choices=None))
		from emdatastorage import ParamDef
		db = db_open_dict(self.form_db_name)
		
		params = []
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=self.__doc__,choices=None))
		params.append(self.get_table())
		
		pbox = ParamDef(name="output_boxsize",vartype="int",desc_short="Box Size",desc_long="An integer value",property=None,defaultunits=db.get("output_boxsize",dfl=128),choices=[])	
		pfo = ParamDef(name="force",vartype="boolean",desc_short="Force Overwrite",desc_long="Whether or not to force overwrite files that already exist",property=None,defaultunits=db.get("force",dfl=False),choices=None)
		pwc = ParamDef(name="write_coords",vartype="boolean",desc_short="Write Coordinates",desc_long="Whether or not to write .box files",property=None,defaultunits=db.get("write_coords",dfl=False),choices=None)
		psuffix = ParamDef(name="suffix",vartype="string",desc_short="Output Suffix", desc_long="This text will be appended to the names of the output files",property=None,defaultunits=db.get("suffix",dfl="_ptcls"),choices=None )
		pwb = ParamDef(name="write_particles",vartype="boolean",desc_short="Write Particles",desc_long="Whether or not box images should be written",property=None,defaultunits=db.get("write_particles",dfl=True),choices=None)
		pinv = ParamDef(name="invert_output",vartype="boolean",desc_short="Invert Pixels",desc_long="Do you want the pixel intensities in the output inverted?",property=None,defaultunits=db.get("invert_output",dfl=False),choices=None)
		pn =  ParamDef(name="normproc",vartype="string",desc_short="Normalize Images",desc_long="How the output box images should be normalized",property=None,defaultunits=db.get("normproc",dfl="normalize.edgemean"),choices=["normalize","normalize.edgemean","normalize.ramp.normvar","None"])
		pop = ParamDef(name="outformat",vartype="string",desc_short="Output Image Format",desc_long="The format of the output box images",property=None,defaultunits=db.get("outformat",dfl="bdb"),choices=self.output_formats)


	   	pwb.dependents = ["invert_output","normproc","outformat","suffix"] # these are things that become disabled when the pwb checkbox is unchecked etc
		
		params.append([pbox,pfo])
		params.append([pwc,pwb])
		params.append([psuffix,pinv])
		params.append(pn)
		params.append(pop)

		return params
	
	def check_params(self,params):
		error_message = []
		if params["output_boxsize"] < 1: error_message.append("Boxsize must be greater than 0.")
		if not params["write_coords"] and not params["write_particles"]: error_message.append("You must choose at least one of the write_coords/write_box_images options")
	
		return error_message
	
	def on_form_ok(self,params):	
		if  params.has_key("filenames") and len(params["filenames"]) == 0:
			self.run_select_files_msg()
			return
		
		error_message = self.check_params(params)
		if len(error_message) >0: 
			self.show_error_message(error_message)
			return
		
		particle_output_names = []
		if params["write_particles"]: particle_output_names = self.get_particle_outnames(params)

		coord_output_names = []
		if params["write_coords"]: coord_output_names = self.get_coord_outnames(params)
		
		if not params["force"]:
			warning = []
			for name in particle_output_names: 
				if file_exists(name): warning.append("%s" %name)
			for name in coord_output_names: 
				if file_exists(name):  warning.append("%s" %name)
			if len(warning) > 0:
				warning.insert(0,"The following files exist and will not be written over")
				warning.append("You can remove them or change the suffix")
				error(warning,"Error")
				return
		
		if len(particle_output_names) > 0:
			self.write_output(params["filenames"],particle_output_names,EMBoxList.write_particles,params["output_boxsize"],"Writing Particles", [params["invert_output"], params["normproc"]])
		if len(coord_output_names) > 0:
			self.write_output(params["filenames"],coord_output_names,EMBoxList.write_coordinates,params["output_boxsize"],"Writing Coordinates")

		from PyQt4 import QtCore
		self.emit(QtCore.SIGNAL("task_idle"))
		self.form.closeEvent(None)
		self.form = None
		self.write_db_entries(params)

	def write_output(self,input_names,output_names,box_list_function,box_size,msg="Writing Output",extra_args=[]):
		n = len(input_names)
		from emapplication import EMProgressDialogModule
		progress = EMProgressDialogModule(get_application(),msg, "Cancel", 0,n,None)
		progress.qt_widget.show()
		prog = 0	

		files_written = []
		for i,output in enumerate(output_names):
			input = input_names[i]
			box_list = EMBoxList()
			box_list.load_boxes_from_database(input)
			if len(extra_args) > 0:box_list_function(box_list,input,output,box_size,*extra_args)
			else:box_list_function(box_list,input,output,box_size)
			files_written.append(output)
			prog += 1
			progress.qt_widget.setValue(prog)
			get_application().processEvents()
			
			if progress.qt_widget.wasCanceled():
				from EMAN2 import remove_file
				for file in files_written: remove_file(file)
				progress.qt_widget.setValue(nim)
				progress.qt_widget.close()
				return
	
		progress.qt_widget.setValue(n)
		progress.qt_widget.close()
			
	def get_particle_outnames(self,params):
		input = params["filenames"]
		outformat = params["outformat"]
		output = []
		for name in input:
			if outformat == "bdb":
				out = "bdb:particles#" + get_file_tag(name) + params["suffix"]
			else:
				out = get_file_tag(name)+ params["suffix"]+"."+outformat
			output.append(out)
		return output
	
	def get_coord_outnames(self,params):
		input = params["filenames"]
		output = []
		for name in input:
			output.append(get_file_tag(name)+".box")
		return output

from emapplication import EMQtWidgetModule
class EMBoxerInspectorModule(EMQtWidgetModule):
	'''
	'''
	def __init__(self,target):
		self.widget = EMBoxerInspector(target)
		EMQtWidgetModule.__init__(self,self.widget)
		
	def get_desktop_hint(self):
		return "inspector"
	
from PyQt4 import QtGui
class EMBoxerInspector(QtGui.QWidget):
	def __init__(self,target) :
		from PyQt4 import QtCore, QtGui
		self.busy = True
		self.tool_dynamic_vbl = None # this will be used to dynamic add widgets as the buttons are changed
		self.dynamic_box_button_widget = None # this will be used to dynamic add widgets as the buttons are changed
		 
		QtGui.QWidget.__init__(self,None)
		self.setWindowIcon(QtGui.QIcon(get_image_directory() +"eman.png"))
		self.target=weakref.ref(target)
		
		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		
		self.tab_widget = QtGui.QTabWidget()
		self.tab_widget.addTab(self.get_main_tab(),"Main")
		self.tab_widget.addTab(self.get_display_tab(),"Display")
		self.vbl.addWidget(self.tab_widget)
		
		self.gen_output_but=QtGui.QPushButton("Write output")
		self.vbl.addWidget(self.gen_output_but)
		
		self.done_but=QtGui.QPushButton("Done")
		self.vbl.addWidget(self.done_but)
					
		
		self.connect(self.done_but,QtCore.SIGNAL("clicked(bool)"),self.on_done)
		self.connect(self.gen_output_but,QtCore.SIGNAL("clicked(bool)"),self.write_output_clicked)
		self.busy = False
	
	def on_done(self):
		self.close()
	
	def closeEvent(self,event):
		self.target().done()
		
	def write_output_clicked(self,val):
		self.target().run_output_dialog()
	
	def get_display_tab(self):
		from PyQt4 import QtCore, QtGui, Qt
		widget = QtGui.QWidget()
		vbl =  QtGui.QVBoxLayout(widget)
		
		#  Insert the plot widget
		viewhbl = QtGui.QVBoxLayout()
		
		self.viewboxes = QtGui.QCheckBox("Particle Window")
		self.viewboxes.setChecked(True)
		self.viewimage = QtGui.QCheckBox("2D Image Window")
		self.viewimage.setChecked(True)
		
		viewhbl.addWidget(self.viewboxes)
		viewhbl.addWidget(self.viewimage)
	
		if self.target().has_thumbs():
			self.viewthumbs = QtGui.QCheckBox("Thumbnails Window")
			self.viewthumbs.setChecked(True)
			viewhbl.addWidget(self.viewthumbs)
	
		viewmanagement = QtGui.QGroupBox("Displayed Windows")
		viewmanagement.setLayout(viewhbl)
		vbl.addWidget(viewmanagement)
	
		
		viewhbl2 = QtGui.QHBoxLayout()
		self.boxformats = QtGui.QComboBox(self)
		self.boxformats.addItem("square with central dot")
		self.boxformats.addItem("square")
		self.boxformats.addItem("circle with central dot")
		self.boxformats.addItem("circle")
		self.boxformats.setEnabled(False)
		viewhbl2.addWidget(self.boxformats)
		
		displayboxes = QtGui.QGroupBox("Displayed Boxes")
		displayboxes.setLayout(viewhbl2)
		vbl.addWidget(displayboxes)
	
	
		self.connect(self.viewboxes,QtCore.SIGNAL("clicked(bool)"),self.view_particles_clicked)
		self.connect(self.viewimage,QtCore.SIGNAL("clicked(bool)"),self.view_2d_window_clicked)
		if self.target().has_thumbs():
			self.connect(self.viewthumbs,QtCore.SIGNAL("clicked(bool)"),self.view_thumbs_clicked)
		
		QtCore.QObject.connect(self.boxformats, QtCore.SIGNAL("currentIndexChanged(QString)"), self.box_format_changed)

		return widget
	def view_particles_clicked(self,val):
		if self.busy: return
		self.target().show_particles_window(val)
	
	def view_2d_window_clicked(self,val):
		if self.busy: return
		self.target().show_2d_window(val)
		
	def view_thumbs_clicked(self,val):
		if self.busy: return
		self.target().show_thumbs_window(val)
		
	def set_thumbs_visible(self,val=True):
		self.busy = True
		self.viewthumbs.setChecked(val)
		self.busy = False

	def set_particles_visible(self,val=True):
		self.busy = True
		self.viewboxes.setChecked(val)
		self.busy = False

	def set_2d_window_visible(self,val=True):
		self.busy = True
		self.viewimage.setChecked(val)
		self.busy = False
		
	def box_format_changed(self,val):
		print "pass"
	
	def get_main_tab(self):
		from PyQt4 import QtCore, QtGui, Qt
		widget = QtGui.QWidget()
		vbl = QtGui.QVBoxLayout(widget)
		vbl.setMargin(0)
		vbl.setSpacing(6)
		
		self.add_top_buttons(vbl)
		
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
		
		vbl.addLayout(box_size_hbl)
	
		self.add_boxing_button_group(vbl)
		
		self.add_bottom_buttons(vbl)
		
		self.connect(self.box_size,QtCore.SIGNAL("editingFinished()"),self.new_box_size)
		return widget
	
	def add_top_buttons(self,layout):
		from PyQt4 import QtCore, QtGui, Qt
		hbl_t=QtGui.QHBoxLayout()
		self.togfreeze=QtGui.QPushButton(QtGui.QIcon(get_image_directory() + "freeze_swirl.png"),"Freeze")
		self.togfreeze.setCheckable(1)
		self.togfreeze.setChecked(False)
		hbl_t.addWidget(self.togfreeze)
			
		layout.addLayout(hbl_t)
		
		QtCore.QObject.connect(self.togfreeze, QtCore.SIGNAL("clicked(bool)"), self.toggle_freeze_clicked)
	
	def set_frozen(self,val):
		self.busy = True
		self.togfreeze.setChecked(int(val))
		self.busy = False
	
	def toggle_freeze_clicked(self,val):
		if self.busy: return
		self.target().set_frozen(val)
		
	def clear_clicked(self,val):
		if self.busy: return
		self.target().clear_boxes()
	
	def add_bottom_buttons(self,layout):
		from PyQt4 import QtCore, QtGui, Qt
		hbl_t=QtGui.QHBoxLayout()
		self.undo_boxes=QtGui.QPushButton("Undo Boxes")
		self.undo_boxes.setToolTip("Recall the last state of the boxes")
		self.undo_boxes.setEnabled(False)
		hbl_t.addWidget(self.undo_boxes)
		self.redo_boxes=QtGui.QPushButton("Redo Boxes")
		self.redo_boxes.setToolTip("Redo the last undo")
		self.redo_boxes.setEnabled(False)
		hbl_t.addWidget(self.redo_boxes)
		layout.addLayout(hbl_t)
		
		hbl_q=QtGui.QHBoxLayout()
		self.quality=QtGui.QLabel("Image Quality:")
		qual_tt = "Assign a quality number to the image. This acts as metadata for your convenience and is displayed in eman2 forms when possible."
		self.quality.setToolTip(qual_tt)
		hbl_q.addWidget(self.quality)
		
		self.image_qualities = QtGui.QComboBox()
		for i in range(0,5):
			self.image_qualities.addItem(str(i))
		self.image_qualities.setCurrentIndex(2)
		self.image_qualities.setToolTip(qual_tt)
		hbl_q.addWidget(self.image_qualities)
		layout.addLayout(hbl_q)
		
		QtCore.QObject.connect(self.redo_boxes, QtCore.SIGNAL("clicked(bool)"), self.redo_boxes_clicked)
		QtCore.QObject.connect(self.undo_boxes, QtCore.SIGNAL("clicked(bool)"), self.undo_boxes_clicked)
		QtCore.QObject.connect(self.image_qualities, QtCore.SIGNAL("currentIndexChanged(QString)"), self.image_quality_changed)
	
	def image_quality_changed(self,val):
		if self.busy: return
		self.target().set_image_quality(int(val))
		
	def set_image_quality(self,val):
		self.busy = True
		for i in range(self.image_qualities.count()):
			if int(self.image_qualities.itemText(i)) == val:
				self.image_qualities.setCurrentIndex(i)
				break
		else:
			raise RuntimeError("Unknow quality: " + str(val))
		self.busy = False
		
	def redo_boxes_clicked(self):
		self.target().redo_boxes()
	
	def undo_boxes_clicked(self,val):
		self.target().undo_boxes()
		
	def enable_undo_redo(self, enable_undo, enable_redo):
		self.undo_boxes.setEnabled(enable_undo)
		self.redo_boxes.setEnabled(enable_redo)
	
	def add_boxing_button_group(self,layout):
		from PyQt4 import QtCore, QtGui, Qt
		
		self.tool_button_group_box = QtGui.QGroupBox("Tools")
		self.tool_button_group_box_vbl = QtGui.QVBoxLayout(self.tool_button_group_box)
#		
		self.tool_dynamic_vbl = QtGui.QVBoxLayout()
		self.tool_tabs = QtGui.QTabWidget()
		self.tool_dynamic_vbl.addWidget(self.tool_tabs)
		self.tool_button_group_box_vbl.addLayout(self.tool_dynamic_vbl,1)
		layout.addWidget(self.tool_button_group_box,0,)
			
		QtCore.QObject.connect(self.tool_tabs,QtCore.SIGNAL("currentChanged(int)"),self.tool_tab_changed)
	
	def add_mouse_tool(self,mouse_tool,name):
		widget = mouse_tool.get_widget()
		self.tool_tabs.addTab(widget,mouse_tool.icon(),name)		
		self.update()
		
	def tool_tab_changed(self,idx):
		self.target().set_main_2d_mouse_mode(str(self.tool_tabs.tabText(idx)))

	def get_desktop_hint(self):
		return "inspector"
	
	def set_box_size(self,value):
		self.busy = True
		self.box_size.setText(str(value))
		self.busy = False

	def new_box_size(self):
		if self.busy: return
		box_size=int(self.box_size.text())
		self.target().set_box_size(box_size)
	
	def keyPressEvent(self,event):
		from PyQt4 import QtCore
		if event.key() == QtCore.Qt.Key_F1:
			try:
				import webbrowser
				webbrowser.open("http://blake.bcm.edu/emanwiki/e2boxer")
				return
			except: pass
			
			try: from PyQt4 import QtWebKit
			except: return
			try:
				try:
					test = self.browser
				except: 
					self.browser = QtWebKit.QWebView()
					self.browser.load(QtCore.QUrl("http://blake.bcm.edu/emanwiki/e2boxer"))
					self.browser.resize(800,800)
				
				if not self.browser.isVisible(): self.browser.show()
			except: pass
			
if __name__ == "__main__":
	main()

	
	