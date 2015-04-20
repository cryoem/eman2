#!/home/ahmad/EMAN2/python/Python-2.5.4-ucs4/bin/python
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
./emboxerbase.py test.mrc --boxsize=128 # single image handling
-or-
./emboxerbase.py test.mrc test2.mrc test3.mrc --boxsize=128 # multiple image handling enabled

The EMBoxerModule is basically the epicenter of everything: functions like "add_box" and "move_box" are probably good starting
points in terms of figuring out how to adapt this code to application specific needs
'''
from optparse import OptionParser
from emapplication import EMApp,get_application
from pyemtbx.boxertools import BigImageCache,BinaryCircleImageCache,Cache
from EMAN2 import file_exists,EMANVERSION,gimme_image_dimensions2D,EMData,get_image_directory,Region,file_exists,gimme_image_dimensions3D,abs_path,get_platform,base_name
from EMAN2db import db_open_dict,db_check_dict,db_close_dict
from EMAN2jsondb import *
from emsprworkflow import workflow_path
from EMAN2 import *

import os,sys,weakref,math, json


TEMPLATE_MIN = 30

EMBOXERBASE_DB = "e2boxercache/base.json"

def my_main():
	progname = os.path.basename(sys.argv[0])
	usage = '''%prog [options] <image> <image2>....

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
	'''

	parser = OptionParser(usage=usage,version=EMANVERSION)
	parser.add_option("--boxsize",type="int",help="Box size in pixels",default=-1)
	(options, args) = parser.parse_args()

	db = js_open_dict(EMBOXERBASE_DB)
	cache_box_size = True
	if options.boxsize == -1:
		cache_box_size = False
		options.boxsize = db.setdefault("box_size",128)

	error_message = check(options,args)
	if len(error_message) > 0:
		error = "\n"
		for e in error_message:
			error += "Error: "+e +"\n"
		parser.error(error)

	if cache_box_size: db["box_size"] = options.boxsize

	application = EMApp()
#	QtCore.QObject.connect(gui, QtCore.SIGNAL("module_idle"), on_idle)

	module = EMBoxerModule(args,options.boxsize)
	module.show_interfaces()
	#this is an example of how to add your own custom tools:
	#module.add_tool(EraseTool,ErasingPanel,erase_radius=2*options.boxsize)
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

def get_database_entry(image_name,key,dfl=None):
	db = js_open_dict(info_name(image_name))
	# First check the full path, then check for basename
	if db.has_key(key):
		return db[key]
	elif dfl != None:
		db[key] = dfl
		return dfl
	else: return None

def set_database_entry(image_name,key,value):
	'''
	write a key/object pair to the Image Database Dictionary associat
	'''
	db = js_open_dict(info_name(image_name))
	db[key] = value

def set_idd_image_entry(image_name,key,image,db_title="e2boxercache/"):
	'''
	caches images associated with autoboxing. Technically 'image' can be any object.
	'''
	db=js_open_dict("{}{}.json".format(db_title,base_name(image_name)))
	db[key]=image

def get_idd_image_entry(image_name,key,db_title="e2boxercache/",dfl=None):
	'''
	Using EMAN2 style image dbs has efficiency payoffs in various ways...
	'''
	db=js_open_dict("{}{}.json".format(db_title,base_name(image_name)))
	try: return db[key]
	except : return dfl

class ThumbsEventHandler:
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
		QtCore.QObject.connect(self.thumbs_window(),QtCore.SIGNAL("mx_mouseup"),self.thumb_image_selected)
		QtCore.QObject.connect(self.thumbs_window(),QtCore.SIGNAL("module_closed"),self.module_closed)

	def thumb_image_selected(self,event,lc):
		if lc == None: return
		im=lc[0]
		self.target().set_current_file_by_idx(im)
		try:
			self.target().get_2d_window().updateGL()
		except: pass #window is closed

	def module_closed(self):

		self.target().thumbs_window_closed()

	def add_mouse_handler(self,handler):
		pass

	def set_mouse_mode(self,name):
		pass

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


class EMBox:
	'''
	A basic encapsulation of a box - it has a central coordinate, a type attribute which can be
	customized for specific boxes, and a score attribute, which could be useful to a particular
	tool.
	Also has convenient functions for moving, getting the associated box image, etc
	'''
	BOX_COLORS = {}
	def __init__(self,x,y,type,score=0.0):
		self.x = x # central x coordinate
		self.y = y # central y coordinate
		self.type = type # type can be customized
		self.score = score # can be some kind of score, such as correlation
		self.image = None # an image


	def __getitem__(self,idx):
		'''
		A highly specialized implementation of __getitem__ - makes it so you can treat an EMBox
		as a list - throughout the emboxer code boxes are stored either as lists in the format
		[x,y,type,score], or are encapsulated as an EMBox. So if you wrote a function
		for dealing with a box stored as a list, this function helps your code to work for
		boxes stored as EMBoxes too
		'''
		if idx == 0: return self.x
		elif idx == 1: return self.y
		elif idx == 2: return self.type
		else: raise RuntimeError("Out of bounds request in EMBox")

	def to_list(self):
		'''
		A way of getting the most important information of an EMBox in a
		@return a list: [self.x,self.y,self.type]
		'''
		return [self.x,self.y,self.type]

	def set_box_color(box_type,box_color,force=False):
		'''
		static - use this function to register a box color with a particular EMBox.type attribute
		This is critical - if you don't register your unique box type using this function you'll
		get a runttime error
		@param box_type a string such as "manual" or "swarm_auto", or "swarm_ref", etc
		@param box_color an RGB list [R,G,B] (floats)
		@param force something you'd set to True if you want to force the overwrite of the old color (previously stored)
		'''
		if not force and EMBox.BOX_COLORS.has_key(box_type):
			# this is just to make sure there are no conflicts - if someone is resetting a color they
			# should know what they're doing
			raise RuntimeError("Error, attempt to set a color key (%s) that already existed" %box_type)
		EMBox.BOX_COLORS[box_type] = box_color

	set_box_color = staticmethod(set_box_color)

	def move(self,dx,dy):
		self.x += dx
		self.y += dy
		self.image = None

	def get_image(self,image_name,box_size,norm=None):

		if self.image == None or self.image.get_xsize() != box_size or self.image.get_ysize() != box_size:
			global BigImageCache
			data=BigImageCache.get_object(image_name).get_image(use_alternate=True,norm=norm) # use alternate is a red herring
			r = Region(self.x-box_size/2,self.y-box_size/2,box_size,box_size)
			self.image = data.get_clip(r)
			if norm != None and str(norm) != "None":
				self.image.process_inplace(norm)

			self.image.set_attr("ptcl_source_coord",[self.x,self.y])
			self.image.set_attr("ptcl_source_image",image_name)

		return self.image

	def reset_image(self): self.image = None

	def get_shape(self,shape_string,box_size):
		if EMBox.BOX_COLORS.has_key(self.type):
			r,g,b = EMBox.BOX_COLORS[self.type]
		else:
			r,g,b = 1.0,0.42,0.71 # hot pint, apparently ;)
		from emshape import EMShape
		shape = EMShape([shape_string,r,g,b,self.x-box_size/2,self.y-box_size/2,self.x+box_size/2,self.y+box_size/2,2.0])
		return shape

	def collision(self,x,y,box_size):

		if x-box_size/2 < self.x and x+box_size/2 > self.x and y-box_size/2 < self.y and y+box_size/2 > self.y: return True
		else: return False


class EMBoxingTool:
	'''
	This class defines the interface necessary for integration of a custom tool, such as an automatic boxer,
	into EMBoxerModule.

	EMBoxingTools need to supply a widget (which is automatically added to a tab widget in the inspector of the main module),
	when the uses activates the tab corresponding to the concrete EMBoxingTool instance then it becomes active, and
	mouse events from the various windows are rerouted to it.Thus you have to make your EMBoxingTool
	capable of reacting to the events that occur in the 2D window and particle stack viewer windows.
	(For events that happen in the 2D window see the Main2DWindowEventHandler, for events that happen
	in the particle stack viewer see the ParticlesWindowEventHandler). The expected interface for acquiring the widget
	and mouse-related functions is defined in this object.

	I (David Woolford) prefer to encapsulate all things concerned with the accompanying widget into a separate class, see
	for example the ErasingPanel and the ManualBoxingPanel, which supply the widgets for the EraseTool and ManualBoxingTool
	objects. That way the EMBoxingTool is the only object that needs to be explicitly aware of the Widget supplying
	class, and it nicely layers the design of everything.
	----
	'''

	def __init__(self,target,**kargs):
		'''
		Any such object will always be given an instance of an EMBoxerModule as the target argument - you
		should keep a weak reference to it. Additionally you can supply key word arguments
		@param target an EMBoxerModule instance
		@param kargs keyword argument initialization support
		called in EMBoxerModule.add_tool
		'''
		pass

	def icon(self):
		'''
		An icon makes your tool easy to identify and the interface look slick. Returning None
		is safe too, which just means the return of the self.unique_name function is all that appears to the user
		@return hopefully a QtGui.QIcon, but you may also return None
		called in EMBoxInspector
		'''
		raise NotImplementedException("Inheriting classes must supply this function")

	def get_widget(self):
		'''
		This function should return a QtGui.QWidget - you can put in it whatever you like. Generally
		you make your widget in this function and also make all of the signal-slot connections -
		@return a QtGui.QWidget with your widgets in it
		called in EMBoxInspector
		'''
		raise NotImplementedException("Inheriting classes must supply this function")

	def unique_name(self):
		'''
		Get a unique name - the name should be information yet succinct. Names currently uses are ManualBoxingTool.BOX_TYPE, "Erase", and "Swarm"
		The name is used as the title in a a Tab Widget, and also is used internally in the code for various purposes (such as
		retrieval)
		@return a string
		called at various locations
		'''
		raise  NotImplementedException("Inheriting classes must supply this function")

	def mouse_move(self,event):
		''''
		required by the Main2DWindowEventHandler
		How shall you respond to the mouse move event?
		@param event a QtGui.QEvent
		'''
		raise  NotImplementedException("Inheriting classes must supply this function")

	def mouse_wheel(self,event):
		'''
		required by the Main2DWindowEventHandler
		How shall you respond to the mouse wheel event?
		@param event a QtGui.QEvent
		'''
		raise  NotImplementedException("Inheriting classes must supply this function")

	def mouse_down(self,event) :
		'''
		required by the Main2DWindowEventHandler
		How shall you respond to the mouse down event?
		@param event a QtGui.QEvent
		'''
		raise  NotImplementedException("Inheriting classes must supply this function")


	def mouse_drag(self,event) :
		'''
		required by the Main2DWindowEventHandlers
		How shall you respond to the mouse drag event?
		@param event a QtGui.QEvent
		'''
		raise  NotImplementedException("Inheriting classes must supply this function")


	def mouse_up(self,event) :
		'''
		required by the Main2DWindowEventHandler
		How shall you respond to the mouse up event?
		@param event a QtGui.QEvent
		'''
		raise  NotImplementedException("Inheriting classes must supply this function")


	def key_press(self,event):
		'''
		required by the Main2DWindowEventHandler
		How shall you respond to the key press event?
		@param event a QtGui.QEvent
		'''
		raise NotImplementedException("Inheriting classes must supply this function")

	def moving_ptcl_established(self,box_num,x,y):
		'''
		This is required by the ParticlesWindowEventHandler
		This function is called when a particle is first clicked on
		@param box_num the box number (int)
		@param x the x coordinate of the mouse
		@param y the y coordinate of the mouse
		'''
		raise NotImplementedException("Inheriting classes must supply this function")

	def move_ptcl(self,box_num,x,y,scale):
		'''
		This is required by the ParticlesWindowEventHandler
		This function is called when a particle is moved
		@param box_num the box number (int)
		@param x the x coordinate of the mouse
		@param y the y coordinate of the mouse
		@param scale the current scale of the particle viewer window, for imposing the correct amount of movement
		'''
		raise NotImplementedException("Inheriting classes must supply this function")

	def release_moving_ptcl(self,box_num,x,y):
		'''
		This is required by the ParticlesWindowEventHandler
		This function is called when a particle is released
		@param box_num the box number (int)
		@param x the x coordinate of the mouse
		@param y the y coordinate of the mouse
		'''
		raise NotImplementedException("Inheriting classes must supply this function")

	def delete_ptcl(self,box_num):
		'''
		This is required by the ParticlesWindowEventHandler
		This function is called when a particle is deleted. It is the responsible of the EMBoxingTool
		to handle particle deletion.
		@param box_num the box number (int)
		'''
		raise NotImplementedException("Inheriting classes must supply this function")

	def get_unique_box_types(self):
		'''
		Should return a list containing all of the types of boxes this object 'owns'
		so for the ManualBoxingTool object this is ["Manual"] and for the SwarmTool it's
		["swarm_auto","swarm_ref","swarm_weak_ref"] etc, This is important for the
		EMBoxerModule - if you click on a box that originated from a different Mouse Mode,
		then we can automatically switch, etc
		@return a list of strings
		'''
		return [] # a default is supplied because some tools may not be concerned with this feature, such as the EraseTool

	def boxes_erased(self,list_of_boxes):
		'''
		When the eraser completes its task a list of boxes is generated.
		The Main2DWindowEventHandler gets this list and then splits
		of the boxes according to which Tool they come from, and hands them on
		to the actual tool using this function. This function exists because some
		tools need to know when associated boxes have been deleted. You can ignore
		this function if you like
		'''
		pass

	def set_current_file(self,file_name,active_tool=False):
		'''
		Called whenever the current file being studied is changed, the active_tool parameter tells the
		tool (this object) whether or not it is the current active tool
		@param file_name a string which is a full file name
		@param avtive_tool a bool indicating whether or not this tool is the active one
		'''
		pass

class EMUnknownBoxType:
	'''
	Error which is thrown when a mouse tool figures out that its operating on box that it doesn't 'own'
	'''
	def __init__(self,type):
		self.type = type

class ErasingPanel:
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
			QtCore.QObject.connect(self.erase_rad_edit,QtCore.SIGNAL("sliderReleased"),self.new_erase_radius) #"editingFinished()"
			QtCore.QObject.connect(self.unerase,QtCore.SIGNAL("clicked(bool)"),self.unerase_checked)

		return self.widget

	def new_erase_radius(self, erase_rad_edit):
		if self.busy: return
		self.target().set_erase_radius(erase_rad_edit)

	def unerase_checked(self,val):
		if self.busy: return
		self.target().toggle_unerase(val)

class ManualBoxingPanel:
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

class EraseTool(EMBoxingTool):
	'''
	A class that knows how to handle mouse erase events for a GUIBox
	'''

	def __init__(self,target,erase_radius=128):
		self.target = weakref.ref(target)
		self.erase_value = 0.1			# erase mode can be either Boxable.ERASE or Boxable.UNERASE
		self.erase_radius = erase_radius
		self.panel_object = None

	def unique_name(self): return "Erase"

	def icon(self):
		from PyQt4 import QtGui
		return QtGui.QIcon(get_image_directory() + "boxer_erase.png")

	def get_widget(self):
		if self.panel_object == None:
			self.panel_object = ErasingPanel(self,self.erase_radius)
		return self.panel_object.get_widget()

	def set_erase_radius(self,val): self.erase_radius = val

	def set_current_file(self,file_name,active_tool=False):
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
		self.get_2d_window().add_eraser_shape("eraser",["circle",.1,.1,.1,m[0],m[1],self.erase_radius,3])
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
			self.get_2d_window().add_eraser_shape("eraser",["circle",.1,.1,.1,m[0],m[1],self.erase_radius,3])
			self.get_2d_window().updateGL()

	def mouse_down(self,event) :
		from emshape import EMShape
		m=self.get_2d_window().scr_to_img((event.x(),event.y()))
		#self.boxable.add_exclusion_area("circle",m[0],m[1],self.erase_radius)
		self.get_2d_window().add_eraser_shape("eraser",["circle",.9,.9,.9,m[0],m[1],self.erase_radius,3])
		self.target().exclusion_area_added("circle",m[0],m[1],self.erase_radius,self.erase_value)

	def mouse_drag(self,event) :
		from emshape import EMShape
		m=self.get_2d_window().scr_to_img((event.x(),event.y()))
		self.get_2d_window().add_eraser_shape("eraser",["circle",.9,.9,.9,m[0],m[1],self.erase_radius,3])
		self.target().exclusion_area_added("circle",m[0],m[1],self.erase_radius,self.erase_value)
		# exclusion_area_added does the OpenGL update calls, so there is no need to do so here

	def mouse_up(self,event) :
		# we have finished erasing

		# make the eraser shape non visible
		self.get_2d_window().add_eraser_shape("None",None)
		self.target().erasing_done(self.erase_value)

	def moving_ptcl_established(self,box_num,x,y):
		box = self.target().get_box(box_num)
		raise EMUnknownBoxType,box.type # this causes the mouse mode to be changed

	def move_ptcl(self,box_num,x,y,scale):
		box = self.target().get_box(box_num)
		raise EMUnknownBoxType,box.type # this causes the mouse mode to be changed

	def release_moving_ptcl(self,box_num,x,y):
		box = self.target().get_box(box_num)
		raise EMUnknownBoxType,box.type # this causes the mouse mode to be changed

	def delete_ptcl(self,box_num):
		box = self.target().get_box(box_num)
		raise EMUnknownBoxType,box.type # this causes the mouse mode to be changed

class ManualBoxingTool:
	'''
	A class that knows how to add, move and remove reference and non reference boxes
	'''
#	SET_BOX_COLOR = True
	BOX_TYPE = "manual"
	EMBox.set_box_color(BOX_TYPE,[1,1,1])
	def __init__(self,target):
		self.target = weakref.ref(target)
		self.moving = None
		self.panel_object = None
		self.moving_data = None
#		if ManualBoxingTool.SET_BOX_COLOR:

#			ManualBoxingTool.SET_BOX_COLOR = False

	def get_widget(self):
		if self.panel_object == None:
			self.panel_object = ManualBoxingPanel(self)
		return self.panel_object.get_widget()


	def icon(self):
		from PyQt4 import QtGui
		return QtGui.QIcon(get_image_directory() + "white_box.png")


	def set_panel_object(self,panel): self.panel_object = panel
	def unique_name(self): return ManualBoxingTool.BOX_TYPE

	def set_current_file(self,file_name,active_tool=False):
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
			box_num = self.target().add_box(m[0],m[1],ManualBoxingTool.BOX_TYPE)
			if self.panel_object.auto_center_checkbox.isChecked():
				self.try_to_center_ref(box_num)

			self.moving=[m,box_num]
		else:
			box = self.target().get_box(box_num)
			if box.type == ManualBoxingTool.BOX_TYPE:
		 		if event.modifiers()&Qt.ShiftModifier :
					self.target().remove_box(box_num)
				else:
					# if we make it here than the we're moving a box
					self.moving=[m,box_num]
					#self.target().moving_box_established(box_num)
			else:
				raise EMUnknownBoxType,box.type

	def mouse_drag(self,event) :
		m=self.get_2d_window().scr_to_img((event.x(),event.y()))
		from PyQt4.QtCore import Qt
		if event.modifiers()&Qt.ShiftModifier:
			box_num = self.target().detect_box_collision(m)
			if ( box_num != -1):
				box = self.target().get_box(box_num)
				if box.type ==  ManualBoxingTool.BOX_TYPE:
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
		self.target().clear_boxes([ManualBoxingTool.BOX_TYPE],cache=True)

	def moving_ptcl_established(self,box_num,x,y):
		box = self.target().get_box(box_num)
		if box.type != ManualBoxingTool.BOX_TYPE:
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
		if box.type != ManualBoxingTool.BOX_TYPE:
			raise EMUnknownBoxType,box.type
		self.target().remove_box(box_num)

	def get_unique_box_types(self):
		return [ManualBoxingTool.BOX_TYPE]

	def boxes_erased(self,list_of_boxes):
		'''
		No need to act here for the manual boxing tool - everything is fine
		'''
		pass



	#TODO: better code reuse, not copy and paste, here
	#COPIED FROM e2boxer's SwarmBoxer !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	def xform_center_propagate(self,box,image_name,template,box_size):
		'''
		Centers a box that was generated in a shrunken image by getting the 'real particle' out of the large
		image on disk and doing a ccf with the template - then I just find the peak and use that to center
		@param box a list like [x,y,type,float] - only x and y are used
		@param image_name the name of the image we're operating on
		@param template the template correctly scaled to the have the same angstrom per pixel as the image (named image_name) stored on disk
		@param box_size the size of the box used to center
		Returns the dx and dy parameters, i.e. does not actually alter the box
		'''
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
	# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

class BoxEventsHandler:
	def __init__(self, target):
		self.target = weakref.ref(target)
		self.box_to_tool_dict = {} # this help automatic changing from one mouse tool to another when a user selects a box of a certain type
		self.busy = False
		self.mouse_handlers = {}  #stores mouse handlers
		self.mouse_handler = None  # the current mouse handler

	def set_mouse_mode(self,name):
		if self.busy: return
		self.mouse_handler = self.mouse_handlers[name]

	def change_event_handler(self,box_type):
		'''
		Called internally in when the EMUnknownBoxType error is caught
		Sets the mouse mode using the new name, and makes sure the inspector is update as well
		'''
		name = self.box_to_tool_dict[box_type]
		self.set_mouse_mode(name)
		self.busy=True
		self.target().set_inspector_tool_mode(name)
		self.busy=False

	def add_mouse_handler(self,handler):
		name = handler.unique_name()
		self.mouse_handlers[name] = handler
		if self.mouse_handler == None: self.mouse_handler = handler

		name = handler.unique_name()
		for box_type in handler.get_unique_box_types():
			if self.box_to_tool_dict.has_key(box_type):
				raise RuntimeError("Error - some EMBoxingTools are using the same box_type name, or the same tool has been added twice (%s)" %box_type)
			self.box_to_tool_dict[box_type] = name

class Main2DWindowEventHandler(BoxEventsHandler):
	'''
	A class that responds to added, moved and removed box signals emitted
	by the by the main 2d image in the EMBoxerModule - this is the image
	display that shows the image that is currently being boxed along
	with any boxed regions
	'''
	def __init__(self,target,main_2d_window):
		BoxEventsHandler.__init__(self,target)
		self.target = weakref.ref(target) # prevent a strong cycle
		self.main_2d_window = main_2d_window

		self.__connect_signals_to_slots()

#	def set_current_file(self,file_name):
#		for value in self.mouse_handlers.values(): value.set_current_file(file_name)

	def __connect_signals_to_slots(self):
		'''
		connects the signals of the main 2D window to the slots of this object
		'''
		from PyQt4 import QtCore
		QtCore.QObject.connect(self.main_2d_window,QtCore.SIGNAL("mousedown"),self.mouse_down)
		QtCore.QObject.connect(self.main_2d_window,QtCore.SIGNAL("mousedrag"),self.mouse_drag)
		QtCore.QObject.connect(self.main_2d_window,QtCore.SIGNAL("mouseup")  ,self.mouse_up  )
		QtCore.QObject.connect(self.main_2d_window,QtCore.SIGNAL("keypress"),self.key_press)
		QtCore.QObject.connect(self.main_2d_window,QtCore.SIGNAL("mousewheel"),self.mouse_wheel)
		QtCore.QObject.connect(self.main_2d_window,QtCore.SIGNAL("mousemove"),self.mouse_move)
		QtCore.QObject.connect(self.main_2d_window,QtCore.SIGNAL("module_closed"),self.module_closed)

	def boxes_erased(self,rm_boxes):
		'''
		When the eraser completes its task a list of boxes is generated.
		They are given to this function, which splits up the boxes according to which Tool
		created them, and then hands on theses on the actual tool. This function exists because some
		tools need to know when associated boxes have been deleted.
		@param rm_boxes a list of EMBox objects
		'''
		name_box_map = {}
		for name,handler in self.mouse_handlers.items():
			name_box_map[name] = handler.get_unique_box_types()

		rm_box_map = {}
		for box in rm_boxes:
			for name,box_types in name_box_map.items():
				if box.type in box_types:
					if rm_box_map.has_key(name):
						rm_box_map[name].append(box)
					else:
						rm_box_map[name] = [box]


		for name,boxes in rm_box_map.items():
			self.mouse_handlers[name].boxes_erased(boxes)

	def mouse_down(self,event) :
		'''
		@param a QtGui.QMouseEvent sent from the EMImage2DWidget
		'''
		if self.mouse_handler == None: return
		try: self.mouse_handler.mouse_down(event)
		except EMUnknownBoxType,data:
			self.change_event_handler(self.box_to_tool_dict[data.type])
			self.mouse_handler.mouse_down(event)

	def mouse_drag(self,event) :
		'''
		@param a QtGui.QMouseEvent sent from the EMImage2DWidget
		'''
		if self.mouse_handler == None: return
		try: self.mouse_handler.mouse_drag(event)
		except EMUnknownBoxType,data:
			self.change_event_handler(self.box_to_tool_dict[data.type])
			self.mouse_handler.mouse_drag(event)

	def change_event_handler(self,name):
		'''
		Called internally in when the EMUnknownBoxType error is caught
		Sets the mouse mode using the new name, and makes sure the inspector is update as well
		'''
		self.set_mouse_mode(name)
		self.target().set_inspector_tool_mode(name)

	def mouse_up(self,event) :
		'''
		@param a QtGui.QMouseEvent sent from the EMImage2DWidget
		'''
		if self.mouse_handler == None: return
		try: self.mouse_handler.mouse_up(event)
		except EMUnknownBoxType,data:
			self.change_event_handler(self.box_to_tool_dict[data.type])
			self.mouse_handler.mouse_up(event)

	def key_press(self,event):
		'''
		@param a QtGui.QKeyEvent sent from the EMImage2DWidget
		'''
		if self.mouse_handler == None: return
#		print "2d window key press"

	def mouse_wheel(self,event):
		'''
		@param a QtGui.QMouseEvent sent from the EMImage2DWidget
		'''
		if self.mouse_handler == None: return
		try: self.mouse_handler.mouse_wheel(event)
		except EMUnknownBoxType,data:
			self.change_event_handler(self.box_to_tool_dict[data.type])
			self.mouse_handler.mouse_wheel(event)

	def mouse_move(self,event):
		'''
		@param a QtGui.QMouseEvent sent from the EMImage2DWidget
		'''
		if self.mouse_handler == None: return
		try: self.mouse_handler.mouse_move(event)
		except EMUnknownBoxType,data:
			self.change_event_handler(self.box_to_tool_dict[data.type])
			self.mouse_handler.mouse_move(event)

	def module_closed(self):
		'''
		Slot that is called when the main 2d window is closed
		'''
		self.target().main_2d_window_closed()

class ParticlesWindowEventHandler(BoxEventsHandler):
	def __init__(self,target,particle_window):
		BoxEventsHandler.__init__(self,target)
		self.target = weakref.ref(target) # prevent a strong cycle
		self.particle_window = particle_window
		self.particle_window.set_reroute_delete(True)

		self.__connect_signals_to_slots()
		self.moving_box_data = None
		self.first_clicked = None

	def __connect_signals_to_slots(self):
		'''
		connects the signals of the main 2D window to the slots of this object
		'''
		from PyQt4 import QtCore
		QtCore.QObject.connect(self.particle_window,QtCore.SIGNAL("mx_image_selected"),self.box_selected)
		QtCore.QObject.connect(self.particle_window,QtCore.SIGNAL("mx_mousedrag"),self.box_moved)
		QtCore.QObject.connect(self.particle_window,QtCore.SIGNAL("mx_mouseup"),self.box_released)
		QtCore.QObject.connect(self.particle_window,QtCore.SIGNAL("mx_boxdeleted"),self.box_image_deleted)
		QtCore.QObject.connect(self.particle_window,QtCore.SIGNAL("module_closed"),self.module_closed)

	def box_selected(self,event,lc):
		if self.mouse_handler == None: return

		if lc == None or lc[0] == None: return
		im=lc[0]
		self.moving_box_data = [event.x(),event.y(),im]
		self.first_clicked = im

		try: self.mouse_handler.moving_ptcl_established(im,event.x(),event.y())
		except EMUnknownBoxType,data:
			self.change_event_handler(data.type)
			self.mouse_handler.moving_ptcl_established(im,event.x(),event.y())
		#self.target().moving_ptcl_established(im,event.x(),event.y())
		#self.target().get_2d_window().set_active(im,.9,.9,.4)
		try:
			self.target().get_2d_window().updateGL()
		except: pass #window is closed

	def box_moved(self,event,scale):
		if self.mouse_handler == None: return

		if self.moving_box_data:
			try: self.mouse_handler.move_ptcl(self.moving_box_data[2],event.x(),event.y(),scale)
			except EMUnknownBoxType,data:
				self.change_event_handler(self.box_to_tool_dict[data.type])
				self.mouse_handler.move_ptcl(self.moving_box_data[2],event.x(),event.y(),scale)
			#self.target().move_ptcl(self.moving_box_data[2],event.x(),event.y(),scale)

#			self.moving_box_data = [event.x(),event.y(),self.moving_box_data[2]]

	def box_released(self,event,lc):
		if lc == None or lc[0] == None: return

		if event.modifiers()&PyQt4.QtCore.Qt.ShiftModifier:
			self.particle_window.remove_particle_image(lc[0],event,True)
			self.particle_window.force_display_update()
			return

		if self.mouse_handler == None: return

		try: self.mouse_handler.release_moving_ptcl(self.first_clicked,event.x(),event.y())
		except EMUnknownBoxType,data:
			self.change_event_handler(self.box_to_tool_dict[data.type])
			self.mouse_handler.move_ptcl(self.moving_box_data[2],event.x(),event.y(),scale)
		#self.target().release_moving_ptcl(self.first_clicked,event.x(),event.y())
		if self.moving_box_data:
			if event.x() == self.moving_box_data[0] and event.y() == self.moving_box_data[1]:
				self.target().scroll_2d_window_to_box(self.first_clicked)

		self.first_clicked = None
		self.moving_box_data = None

	def box_image_deleted(self,event,lc,force_image_mx_remove=True):
		if lc == None or lc[0] == None: return

		box_num = lc[0]
		try: self.mouse_handler.delete_ptcl(box_num)
		except EMUnknownBoxType,data:
			self.change_event_handler(self.box_to_tool_dict[data.type])
			self.mouse_handler.delete_ptcl(box_num)

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
		from emapplication import EMProgressDialog

		if shrink == None or shrink<1.5 : shrink = EMThumbsTools.get_image_thumb_shrink(image_names[0])

		application = get_application()
		nim = len(image_names)
		thumbs = [None for i in range(nim)]
		progress = EMProgressDialog("Generating Thumbnails", "Abort", 0, nim,None)
		progress.show()
		prog = 0
		for i in range(nim):
	#				thumb = self.get_image_thumb(i)
			thumb = get_idd_image_entry(image_names[i],"image_thumb")
			if thumb == None:
				thumb=EMData(image_names[i],0)
				thumb.process_inplace("math.meanshrink",{"n":shrink})
				thumb.process_inplace("normalize.edgemean") # if there are lots than they =should all have the same contrast
				thumb.set_attr("display_name",base_name(image_names[i]))
				set_idd_image_entry(image_names[i],"image_thumb",thumb)

			prog += 1
			progress.setValue(prog)
			application.processEvents()

			if progress.wasCanceled():
				progress.setValue(nim)
				progress.close()
				return None

		progress.setValue(nim)
		progress.close()

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
		if shrink==1 : print "WARNING: input images are too small. Something is likely wrong with your file selection"

		return shrink

	get_image_thumb_shrink = staticmethod(get_image_thumb_shrink)



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

		db = js_open_dict(EMBOXERBASE_DB)
		self.shape_string = EMBoxerInspector.PTCL_SHAPE_MAP[db.setdefault("ptcl_display_shape","square")]

		self.box_color_dict = {}

	def set_shape(self,shape):
		if shape != self.shape_string:
			self.shape_string = shape
			self.reset_shapes()
			return True

		return False

	def get_boxes(self,as_dict=False):
		'''
		Get the list of boxes, optionally get them as a dict - idx is the key, box itself is the value
		@param as_dict force the return object to be a dictionary, not a list
		'''
		if not as_dict: return self.boxes
		else:
			ret = {}
			for i,box in enumerate(self.boxes): ret[i] = box
			return ret

	def get_boxes_filt(self,filt,as_dict=False):
		'''
		A way of getting all of the boxes of a certain type, for example
		self.get_boxes_filt("swarm_auto")
		@param a filter, will be compared against the EMBox.type attribute, e.g. "manual","swarm_auto","swarm_ref"
		@param as_dict - results are returned as a dict, key is the box number, value is the box itself
		@return a list of boxes that have the type filt - or - if as_dict is supplied a dict is returned and the keys are box numbers, the values are the boxes
		'''
		if as_dict: ret = {}
		else: ret = []

		for i,box in enumerate(self.boxes):
			if filt == box.type:
				if as_dict: ret[i] = box
				else: ret.append(box)

		return ret

	def clear_boxes(self,types,cache=False):
		for i in xrange(len(self.boxes)-1,-1,-1):
			if self.boxes[i].type in types:
				self.boxes.pop(i)
				self.shapes.pop(i)

		if cache:
			self.save_boxes_to_database(self.target().current_file())

	def add_box(self,x,y,type=ManualBoxingTool.BOX_TYPE,score=0.0):
		'''
		Appends a box to the end of the list of the boxes
		@return the index of the box that was just added
		'''
		self.boxes.append(EMBox(x,y,type,score))
		self.shapes.append(None)
		return len(self.boxes)-1

	def add_boxes(self,boxes):
		'''
		boxes should be a list like [[x,y,type],[x,y,type],....[int,int,string]]
		'''

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

	def get_box(self,box_number):
		'''
		@param box_number the number of the box for which you want to get
		'''
		return self.boxes[box_number]

	def set_box(self,box,box_number):
		'''
		Overwrites the box at the given index
		'''
		self.boxes[box_number] = box
		self.shapes[box_number] = None

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
		return [box.get_image(image_name,box_size,"normalize.edgemean") for box in self.boxes]

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
#		self.cache_old_boxes()
		self.shapes.pop(idx)
		return self.boxes.pop(idx)


	def get_boxes_for_database(self):
		return [[box.x,box.y,box.type] for box in self.boxes]

	def save_boxes_to_database(self,image_name):

		set_database_entry(image_name,"boxes",self.get_boxes_for_database())

	def load_boxes_from_database(self,image_name,reset=True):
		if reset:
			self.boxes = []
			self.shapes = []
			self.undo_cache = []

		data = get_database_entry(image_name,"boxes")
		if data != None:
			for x,y,type in data:
				self.add_box(x,y,type)

	def exclude_from_scaled_image(self,exclusion_image,subsample_rate):
		action = False
		for i in xrange(len(self.boxes)-1,-1,-1):
			box = self.boxes[i]
			x = int(box.x/subsample_rate)
			y = int(box.y/subsample_rate)
			if exclusion_image.get(x,y):
				self.pop(i)
				action = True

		return action

	def remove_boxes(self,l):
		'''
		@param l a list of indices IN REVERSE ORDER
		'''
		for i in l: self.pop(i)

	def write_particles(self,input_file_name,out_file_name,box_size,invert=False,normproc=None,noedges=False):
		hdr=EMData(input_file_name,0,True)
		xsize=hdr["nx"]
		ysize=hdr["ny"]
		db = js_open_dict('info/project.json')
		try: apix = db.get('global.apix')
		except:
			try:
				apix=hdr["apix_x"]
			except:
				apix=1.0
				
		# If someone did a generate output with one box size, then later increased it and repeated the process_events
		# with noedges set, then there could be some residual images with the wrong size left at the end of the file
		try: os.unlink(out_file_name)
		except: pass
		
		j=0
		for i,box in enumerate(self.boxes):
			if noedges :
				if box[0]-box_size/2<0 or box[1]-box_size/2<0 or box[0]+box_size/2>=xsize or box[1]+box_size/2>=ysize : continue
				
			image = box.get_image(input_file_name,box_size,norm=normproc)
			if invert: image.mult(-1)
			if str(normproc) != "None": image.process_inplace(normproc)
			try:
				if apix != None:
					image.set_attr('apix_x', float(apix))
					image.set_attr('apix_y', float(apix))
					image.set_attr('apix_z', float(apix))
			except: pass
			image.write_image(out_file_name,j)
			j+=1


	def write_coordinates(self,input_file_name,out_file_name,box_size):
		f = open(out_file_name,'w')
		if out_file_name.endswith('json'):
			# Write .json file
			coord_str = 'coordinates'
			data = {input_file_name: {coord_str : [], 'box_size': box_size}}
			for box in self.boxes:
				data[input_file_name][coord_str].append((box.x, box.y))
			f.write(json.dumps(data, sort_keys=True, indent=2, separators=(',', ': ')))
		else:
			# Write .box file
			for box in self.boxes:
				xc = box.x-box_size/2
				yc = box.y-box_size/2
				f.write(str(int(xc))+'\t'+str(int(yc))+'\t'+str(box_size)+'\t'+str(box_size)+'\n')
		f.close()


class EMBoxerModuleVitals(object):
	'''
	If you want to run autoboxing without the gui, then the SwarmBoxer still needs
	mediator functionality
	'''

	def __init__(self, file_names=[], box_size=128):
		'''
		@file_name the name of a file on disk
		@exception RuntimeError raised if the file does not exist
		'''
		self.file_names = file_names # a list of file names
		self.current_idx = None # an index into self.file_names
		self.box_size = box_size # the current box size
		self.box_list = EMBoxList(self)

	def set_status_message(self, mesg, timeout=5000, process_events=False):
		print mesg

	def load_default_status_msg(self):
		pass


	def clear_boxes(self, type, cache=False):
		self.box_list.clear_boxes(type,cache=cache)


	def get_subsample_rate(self):
		'''

		'''
		return int(math.ceil(float(self.box_size)/float(TEMPLATE_MIN)))

	def get_box_type(self, box_number):
		'''
		@param box_number the number of the box for which you want to get the type i.e. that which was returned from the detect_box_collision
		'''
		return self.box_list.get_box_type(box_number)

	def get_box(self, box_number):
		'''
		@param box_number the number of the box for which you want to get
		'''
		return self.box_list.get_box(box_number)

	def get_boxes_filt(self, filt, as_dict=False):
		'''
		A way of getting all of the boxes of a certain type, for example
		self.get_boxes_filt("swarm_auto")
		@param a filter, will be compared against the EMBox.type attribute, e.g. "manual","swarm_auto","swarm_ref"
		@param as_dict - results are returned as a dict, key is the box number, value is the box itself
		@return a list of boxes that have the type filt - or - if as_dict is supplied a dict is returned and the keys are box numbers, the values are the boxes
		'''
		return self.box_list.get_boxes_filt(filt,as_dict)

	def get_boxes(self, as_dict=False):
		'''
		A way of getting all of the boxes as a list or a dict
		@param as_dict - results are returned as a dict, key is the box number, value is the box itself
		@return a list of boxes  - or - if as_dict is supplied a dict is returned and the keys are box numbers, the values are the boxes
		'''
		return self.box_list.get_boxes(as_dict)

	def set_box(self, box, box_number, update_display=False):
		'''
		@param box_number the number of the box for which you want to get
		'''
		self.box_list.set_box(box,box_number)
		if update_display:
			self.full_box_update()

	def detect_box_collision(self, data):
		return self.box_list.detect_collision(data[0], data[1], self.box_size)


	def particle_selected(self, box_number):
		box = self.box_list[box_number]
		self.box_placement_update_exclusion_image(box.x,box.y)
		return box

	def box_released(self, box_number):
		box = self.box_list[box_number]
		self.box_placement_update_exclusion_image(box.x,box.y)
		return box

	def add_boxes(self, boxes, update_gl=True):
		'''
		boxes should be a list like [[x,y,type],[x,y,type],....[int,int,string]]
		'''
		self.box_list.add_boxes(boxes)
		self.box_list.save_boxes_to_database(self.current_file())

	def add_box(self, x, y, type=ManualBoxingTool.BOX_TYPE):
		self.box_placement_update_exclusion_image(x,y)
		box_num = self.box_list.add_box(x,y,type=type)
		self.box_list.save_boxes_to_database(self.current_file())
		return box_num

	def box_placement_update_exclusion_image_n(self, box_num, val=0.0, force=False):
		box = self.box_list.get_box(box_num)
		self.box_placement_update_exclusion_image(box.x,box.y,val,force)

	def box_placement_update_exclusion_image(self, x, y, val=0.0, force=False):
		exclusion_image = self.get_exclusion_image()
		if exclusion_image != None:
			sr = self.get_subsample_rate()
			xx,yy = int(x/sr),int(y/sr)
			if force or exclusion_image.get(xx,yy):
				from EMAN2 import BoxingTools
				global BinaryCircleImageCache
				mask = BinaryCircleImageCache.get_image_directly(int(self.box_size/(2*sr)))
				BoxingTools.set_region(self.get_exclusion_image(),mask,xx,yy,val)
				set_idd_image_entry(self.current_file(),ScaledExclusionImage.database_name,self.get_exclusion_image())
				if self.main_2d_window:
					self.main_2d_window.set_other_data(self.get_exclusion_image(),self.get_subsample_rate(),True)
					self.main_2d_window.updateGL()

	def remove_boxes(self, box_numbers, update_gl=True):
		'''
		Removes a list of box numbers from the display and also those that are stored in the local database
		@param box_numbers a list of integer box numbers
		'''
		self.box_list.remove_boxes(box_numbers)
		self.box_list.save_boxes_to_database(self.current_file())
		self.full_box_update(update_gl)
		self.load_default_status_msg()

	def remove_box(self, box_number, exclude_region=False):
		'''
		Removes the box from those that are stored and those that are displayed. Optionally adds the area
		defined by the removed box into the exclusion image
		@param box_number int, a box number, in terms of the boxes stored by self.box_list
		@param exclude_region bool, if True causes a small region in defined by the given box to painted into the exclusion image
		'''
		box = self.box_list.remove_box(box_number)
		self.box_list.save_boxes_to_database(self.current_file())

		if exclude_region:
			self.box_placement_update_exclusion_image(box.x,box.y,0.1,force=True)

		self.full_box_update()
		self.load_default_status_msg()

	def full_box_update(self, update_gl=True):
		pass

	def move_box(self, box_number, dx, dy):
		self.box_list.move_box(box_number,dx,dy)
		self.box_list.save_boxes_to_database(self.current_file())

	def get_exclusion_image(self, mark_boxes=False):
		'''
		@mark_boxes if true the exclusion image is copied and the locations of the current boxes are painted in as excluded regions
		This is useful for autoboxers - they  obviously dont want to box any region that already has a box in it (such as a manual box,
		or a previously autoboxed box)
		'''
		exc_image = ScaledExclusionImageCache.get_image(self.current_file(),self.get_subsample_rate())
		if not mark_boxes:
			return exc_image
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

	def exclusion_area_added(self, typeofexclusion, x, y, radius, val):
		xx = int(x/self.get_subsample_rate())
		yy = int(y/self.get_subsample_rate())

		rr = int(radius/self.get_subsample_rate())
		global BinaryCircleImageCache
		mask = BinaryCircleImageCache.get_image_directly(rr)

		from EMAN2 import BoxingTools
		BoxingTools.set_region(self.get_exclusion_image(),mask,xx,yy,val)

	def set_current_file_by_idx(self, idx):
		if len(self.file_names) <= idx:
			raise RuntimeError("The index is beyond the length of the file names list")

		if idx != self.current_idx:
			self.current_idx = idx
			self.set_current_file(self.file_names[idx])

	def set_image_quality(self, val):
		set_database_entry(self.current_file(),"quality",val)

	def current_file(self):
		return self.file_names[self.current_idx]

	def get_box_size(self):
		return self.box_size

	def set_box_size(self,box_size):
		self.box_size = box_size
		self.box_list.reset_images()
		self.box_list.reset_shapes()
		self.full_box_update()

import PyQt4
class EMBoxerModule(EMBoxerModuleVitals, PyQt4.QtCore.QObject):
	'''
	The EMBoxerModule is like a coordinator. It has 4 widgets: 1 inspector, 1 2D window viewer, and 2 particle
	stack viewers (one for viewing boxed particles, one for viewing thumbnails).
	This module is essentially a Mediator (see Design Patterns) - it coordinates the activities of several EMAN2 modules
	that would otherwise not necessary interact. Overall the interactions can be complicated and this class is an
	attempt to correctly granulate the overall design and the complexity of the classes involved.
	'''
	def __init__(self,file_names=[],box_size=128):
		'''
		@file_name the name of a file on disk
		@exception RuntimeError raised if the file does not exist
		'''
		EMBoxerModuleVitals.__init__(self, file_names=file_names, box_size=box_size)
		PyQt4.QtCore.QObject.__init__(self)

		self.signal_slot_handlers = {} # this is a dictionary, keys are (somewhat random) names, values are event handlers such as Main2DWindowEventHandler. This dict has the only reference to the event handlers
		self.tools = {} # this is just to keep track of all the tools that have been added
		self.current_tool = None # stores the name of the current tool
		self.inspector = None # this will be a Qt style inspector
#		self.inspector_module = None # the wrapping object of self.inspector
		self.main_2d_window = None # this will be the main 2D image display, showing boxed particles etc
		self.particles_window = None # this will be the window displaying the picked particles
		self.thumbs_window = None # this will be the window showing the thumbnails, enabling the user to change between 2D raw data
		self.image_thumbs = None # image_thumbs is a list of thumbnail images
		# self.moving_box = None
		self.output_task = None # will be an EMAN2 style form for writing output
		# initialized the 2D window
		self.__init_main_2d_window()
		if len(self.file_names) > 1:
			self.__init_thumbs_window()
		# initialize the inspector
		self.__init_inspector()
		# this is an example of how to add your own custom tools:
		self.add_tool(ManualBoxingTool)
		self.add_tool(EraseTool,erase_radius=2*box_size)

	# Method overrides
	def set_status_message(self,mesg,timeout=5000,process_events=False):
		if self.inspector != None:
			self.inspector.set_status_message(mesg,timeout)
			if process_events: get_application().processEvents()

	def load_default_status_msg(self):
		self.set_status_message("%d Boxes" %(len(self.box_list)), 0, False)

	def exclusion_area_added(self,typeofexclusion,x,y,radius,val):
		EMBoxerModuleVitals.exclusion_area_added(self,typeofexclusion,x,y,radius,val)

		if self.main_2d_window:
			self.main_2d_window.set_other_data(self.get_exclusion_image(),self.get_subsample_rate(),True)
			self.main_2d_window.updateGL()

	def move_box(self,box_number,dx,dy):
		EMBoxerModuleVitals.move_box(self,box_number,dx,dy)

		if self.particles_window:
			self.particles_window.set_data(self.box_list.get_particle_images(self.current_file(), self.box_size))
			self.particles_window.updateGL()
		if self.main_2d_window:
			self.main_2d_window.add_shape(box_number,self.box_list.get_shape(box_number,self.box_size))
			self.main_2d_window.updateGL()

	def full_box_update(self,update_gl=True):
		EMBoxerModuleVitals.full_box_update(self)

		if self.particles_window != None:
			self.particles_window.set_data(self.box_list.get_particle_images(self.current_file(), self.box_size))
			if update_gl:self.particles_window.updateGL()
		if self.main_2d_window != None:
			self.main_2d_window.set_shapes(self.box_list.get_box_shapes(self.box_size))
			if update_gl:self.main_2d_window.updateGL()

	def add_boxes(self,boxes,update_gl=True):
		'''
		boxes should be a list like [[x,y,type],[x,y,type],....[int,int,string]]
		'''
		if self.particles_window == None:
			self.__init_particles_window()
			get_application().show_specific(self.particles_window)

		EMBoxerModuleVitals.add_boxes(self, boxes)

		if self.particles_window:
			self.particles_window.set_data(self.box_list.get_particle_images(self.current_file(), self.box_size))
			if update_gl: self.particles_window.updateGL()
		if self.main_2d_window:
			self.main_2d_window.update_shapes(self.box_list.get_box_shapes(self.box_size))
			if update_gl: self.main_2d_window.updateGL()
		self.load_default_status_msg()

	def add_box(self,x,y,type=ManualBoxingTool.BOX_TYPE):
		if self.particles_window == None:
			self.__init_particles_window()
			get_application().show_specific(self.particles_window)
		box_num = EMBoxerModuleVitals.add_box(self, x, y,type=type)
		if self.particles_window:
			self.particles_window.set_data(self.box_list.get_particle_images(self.current_file(), self.box_size))
			self.particles_window.updateGL()
		if self.main_2d_window:
			self.main_2d_window.update_shapes(self.box_list.get_box_shapes(self.box_size))
			self.main_2d_window.updateGL()
		self.load_default_status_msg()
		return box_num

	def clear_boxes(self, type, cache=False):
		EMBoxerModuleVitals.clear_boxes(self, type, cache=cache)

		if self.particles_window:
			self.particles_window.set_data(self.box_list.get_particle_images(self.current_file(), self.box_size))
			self.particles_window.updateGL()
		if self.main_2d_window:
			self.main_2d_window.set_shapes(self.box_list.get_box_shapes(self.box_size))
			self.main_2d_window.updateGL()
		self.load_default_status_msg()

	def particle_selected(self,box_number):
		box = EMBoxerModuleVitals.particle_selected(self,box_number)
		if self.main_2d_window: self.main_2d_window.register_scroll_motion(box.x,box.y)

	# subclass methods
	def has_thumbs(self):
		return self.image_thumbs != None

	def erasing_done(self,erase_mode):
		set_idd_image_entry(self.current_file(),ScaledExclusionImage.database_name,self.get_exclusion_image())

		exclusion_image = self.get_exclusion_image()
		subsample_rate = self.get_subsample_rate()

		action = False
		rm_idxs = []
		rm_boxes = []
		for i in xrange(len(self.box_list)-1,-1,-1):
			box = self.box_list.get_box(i)
			x = int(box.x/subsample_rate)
			y = int(box.y/subsample_rate)
			if exclusion_image.get(x,y):
				rm_idxs.append(i)
				rm_boxes.append(box)
				action = True

		if action:
			self.box_list.remove_boxes(rm_idxs)
			self.signal_slot_handlers["2d_window"].boxes_erased(rm_boxes)
			self.full_box_update()

		self.load_default_status_msg()


	def set_ptcl_display_shape(self,shape):
		if self.box_list.set_shape(shape):
			if self.main_2d_window:
				self.main_2d_window.set_shapes(self.box_list.get_box_shapes(self.box_size))
				self.main_2d_window.updateGL()


	def scroll_2d_window_to_box(self,box_number):
		if self.main_2d_window:
			box = self.box_list.get_box(box_number)
			self.main_2d_window.register_scroll_motion(box.x,box.y)

	def set_main_2d_mouse_mode(self,mode):
		self.current_tool = mode
		if self.main_2d_window != None:
			self.main_2d_window.add_eraser_shape("None",None)
			self.main_2d_window.updateGL()

		for mouse_handler in self.signal_slot_handlers.values():
			mouse_handler.set_mouse_mode(mode)


	def set_inspector_tool_mode(self,mode):
		self.current_tool = mode
		self.inspector.set_tool_mode(mode)

		for mouse_handler in self.signal_slot_handlers.values():
			mouse_handler.set_mouse_mode(mode)


	def done(self):
		if self.main_2d_window != None:
			E2saveappwin("e2boxer","image",self.main_2d_window.qt_parent)
			self.main_2d_window.close()

		if self.thumbs_window != None:
			E2saveappwin("e2boxer","thumbs",self.thumbs_window.qt_parent)
			self.thumbs_window.close()

		if self.particles_window != None:
			E2saveappwin("e2boxer","particles",self.particles_window.qt_parent)
			self.particles_window.close()
		self.emit(PyQt4.QtCore.SIGNAL("module_closed"))

	def run_output_dialog(self):
		if self.current_tool=='Gauss':
			print "\n\nThis operation has been deactivated for Gauss mode.\n\nPlease use sxwindow.py for windowing!\n\n"
			error("This operation has been deactivated for Gauss mode.\n\nPlease use sxwindow.py for windowing!","Error")
			return
		from emsprworkflow import E2BoxerProgramOutputTask
		if self.output_task != None: return
		from PyQt4 import QtCore
		self.output_task = EMBoxerWriteOutputTask(self.file_names, dfl_boxsize=self.box_size, current_tool=self.current_tool)
		QtCore.QObject.connect(self.output_task.emitter(),QtCore.SIGNAL("task_idle"),self.on_output_task_idle)
		self.output_task.run_form()

	def on_output_task_idle(self):
		self.output_task = None

	def __init_thumbs_window(self,redo_thumbs=False):
		if len(self.file_names) == 0: raise RuntimeError("Will not make a thumbs window if the number of images is zero")

		if self.thumbs_window == None:
			from PyQt4 import QtCore
			get_application().setOverrideCursor(QtCore.Qt.BusyCursor)


			if self.image_thumbs == None or redo_thumbs:
				self.image_thumbs = EMThumbsTools.gen_thumbs(self.file_names)
			if self.image_thumbs == None:
				sys.exit(1)

			from emimagemx import EMImageMXWidget
			self.thumbs_window=EMImageMXWidget(application=get_application())

			self.thumbs_window.set_data(self.image_thumbs,soft_delete=True)
			self.thumbs_window.set_mouse_mode("App")
			self.thumbs_window.setWindowTitle("Thumbnails")
			self.signal_slot_handlers["thumbs_window"] = ThumbsEventHandler(self,self.thumbs_window)
			for tool in self.tools.values():
				self.signal_slot_handlers["thumbs_window"].add_mouse_handler(tool)

			get_application().setOverrideCursor(QtCore.Qt.ArrowCursor)

	def thumbs_window_closed(self):
		self.thumbs_window = None
		if self.inspector:
			self.inspector.set_thumbs_visible(False)

	def __init_inspector(self):
		if self.inspector == None:
			self.inspector = EMBoxerInspector(self)
			self.inspector.set_box_size(self.box_size)

			# try:
			# 	#self.current_file() might fail
			# 	frozen = get_database_entry(self.current_file(),"frozen",dfl=False)
			# 	self.inspector.set_frozen(frozen)
			# except: pass # inspector sets the frozen button to false by default

	def get_inspector(self):
		return self.inspector

	def __init_particles_window(self):
		if self.particles_window == None:
			from emimagemx import EMImageMXWidget
			self.particles_window=EMImageMXWidget(application=get_application())

			self.particles_window.set_mouse_mode("App")
			self.particles_window.setWindowTitle("Particles")
			self.signal_slot_handlers["particles_window"] = ParticlesWindowEventHandler(self,self.particles_window)
			for tool in self.tools.values():
				self.signal_slot_handlers["particles_window"].add_mouse_handler(tool)

	def particles_window_closed(self):
		self.particles_window = None
		if self.inspector:
			self.inspector.set_particles_visible(False)

	def show_thumbs_window(self,bool):
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
			self.__update_2d_window(self.current_file())
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
		if len(self.file_names) > 0:
			self.set_current_file_by_idx(0)

		if self.main_2d_window != None:
			get_application().show_specific(self.main_2d_window)
			self.main_2d_window.optimally_resize()
			E2loadappwin("e2boxer","image",self.main_2d_window.qt_parent)

		if self.thumbs_window != None:
			get_application().show_specific(self.thumbs_window)
			self.thumbs_window.optimally_resize()
			E2loadappwin("e2boxer","thumbs",self.thumbs_window.qt_parent)

		if self.inspector != None:
			get_application().show_specific(self.inspector)
			E2loadappwin("e2boxer","main",self.inspector)

		if self.particles_window != None:
			get_application().show_specific(self.particles_window)
			self.particles_window.optimally_resize()
			E2loadappwin("e2boxer","particles",self.particles_window.qt_parent)

	def __update_2d_window(self,file_name):
		self.set_status_message("Reading %s..." %file_name,0,True)
		global BigImageCache
		data=BigImageCache.get_object(file_name).get_image(use_alternate=True)

		if get_idd_image_entry(file_name,ScaledExclusionImage.database_name) != None:
			self.main_2d_window.set_other_data(self.get_exclusion_image(),self.get_subsample_rate(),True)
		else:
			self.main_2d_window.set_other_data(None,self.get_subsample_rate(),True)

		self.main_2d_window.set_data(data,file_name)
		self.main_2d_window.force_display_update()
		self.set_status_message("Read Image Done",1000,True)
		# frozen = get_database_entry(file_name,"frozen",dfl=False)
		# if frozen == None:
		# 	set_database_entry(file_name,"frozen",False)
		# 	frozen = False
		# self.main_2d_window.set_frozen(frozen)
	# def set_frozen(self,val):
	# 	set_database_entry(self.current_file(),"frozen",val)
	# 	self.main_2d_window.set_frozen(val)
	# 	self.main_2d_window.updateGL()

	def set_current_file(self,file_name):
		from PyQt4 import QtCore
		get_application().setOverrideCursor(QtCore.Qt.BusyCursor)

		if not file_exists(file_name): raise RuntimeError("The file %s does not exist" %file_name)

		if self.main_2d_window != None:
	   	   	self.__update_2d_window(file_name)
			if self.inspector != None:
				#self.inspector.set_frozen(get_database_entry(file_name,"frozen",dfl=False))
				self.inspector.set_image_quality(get_database_entry(file_name,"quality",dfl=2))

			# the boxes should be loaded from the database, if possible
			self.box_list.load_boxes_from_database(file_name)
			# if self.inspector:
			# 	self.inspector.enable_undo_redo(self.box_list.is_undoable(),self.box_list.is_redoable())
			self.main_2d_window.set_shapes(self.box_list.get_box_shapes(self.box_size))
			self.main_2d_window.updateGL()

			if self.particles_window == None: self.__init_particles_window()
			particles = self.box_list.get_particle_images(self.current_file(),self.box_size)
			self.particles_window.set_data(particles)
			self.particles_window.updateGL()

			self.load_default_status_msg()

		for name, mouse_handler in self.tools.items():
			mouse_handler.set_current_file(file_name,name==self.current_tool)

		get_application().setOverrideCursor(QtCore.Qt.ArrowCursor)

	def __init_main_2d_window(self):
		from emimage2d import EMImage2DWidget
		if self.main_2d_window == None:

			self.main_2d_window= EMImage2DWidget(application=get_application())

			self.main_2d_window.set_mouse_mode(0)

			self.signal_slot_handlers["2d_window"] = Main2DWindowEventHandler(self,self.main_2d_window)
			for tool in self.tools.values():
				self.signal_slot_handlers["2d_window"].add_mouse_handler(tool)

			get_application().show_specific(self.main_2d_window)

	def get_2d_window(self):
		return self.main_2d_window

	def main_2d_window_closed(self):
		self.main_2d_window = None
		if self.inspector:
			self.inspector.set_2d_window_visible(False)

	def add_tool(self,event_tool_class,**kargs):
		event_tool = event_tool_class(self,**kargs)
		if self.current_idx != None:
			event_tool.set_current_file(self.current_file(),False)

		name = event_tool.unique_name()
		if self.current_tool == None: self.current_tool = name
		self.tools[name] = event_tool

		for mouse_handler in self.signal_slot_handlers.values():
			mouse_handler.add_mouse_handler(event_tool)

		self.inspector.add_mouse_tool(event_tool)

from emsprworkflow import WorkFlowTask
from emapplication import error
class EMBoxerWriteOutputTask(WorkFlowTask):
	"""Use this form for writing boxed particles and/or coordinate files to disk."""
	def __init__(self,file_names=[],output_formats=["hdf","spi","img","bdb"],dfl_boxsize=128, current_tool=None):
		WorkFlowTask.__init__(self)
		self.window_title = "Write Particle Output"
		self.form_db_name = EMBOXERBASE_DB
		self.file_names = file_names
		self.output_formats = output_formats
		self.dfl_boxsize = dfl_boxsize
		self.current_tool = current_tool

	def get_table(self):
		from emform import EM2DFileTable,EMFileTable,int_lt
		table = EM2DFileTable(self.file_names,desc_short="Raw Data",desc_long="")
		table.add_column_data(EMFileTable.EMColumnData("Stored Boxes",EMBoxerWriteOutputTask.get_num_boxes,"The number of stored boxes",int_lt))
		table.add_column_data(EMFileTable.EMColumnData("Quality",EMBoxerWriteOutputTask.get_quality,"Quality metadata score stored in local database",int_lt))

		return table

	@staticmethod
	def get_quality(file_name):
		'''
		A static function for getting the number of boxes associated with each file
		'''
		val = get_database_entry(file_name,"quality")

		if val == None: return "-"
		else: return str(val)

	@staticmethod
	def get_num_boxes(file_name):
		'''
		A static function for getting the number of boxes associated with each file
		'''
		box_list = EMBoxList()
		box_list.load_boxes_from_database(file_name)
#		print "lbd ",file_name
		return str(len(box_list))

	def get_params(self):
#		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=E2CTFGenericTask.documentation_string,choices=None))
		from emdatastorage import ParamDef
		db = js_open_dict(self.form_db_name)
		is_gauss = self.current_tool == 'Gauss'

		# params = []
		# params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=self.__doc__,choices=None))
		# params.append(self.get_table())
		#
		# pbox = ParamDef(name="output_boxsize",vartype="int",desc_short="Box Size",desc_long="An integer value",property=None,defaultunits=db.setdefault("output_boxsize",self.dfl_boxsize),choices=[])
		# pfo = ParamDef(name="force",vartype="boolean",desc_short="Force Overwrite",desc_long="Whether or not to force overwrite files that already exist",property=None,defaultunits=db.setdefault("force",False),choices=None)
		# params.append([pbox,pfo])
		#
		# pwc = ParamDef(name="write_coords",vartype="boolean",desc_short="Write Coordinates",desc_long="Whether or not to write .box files",property=None,defaultunits=db.setdefault("write_coords", is_gauss),choices=None)
		# if not is_gauss:
		# 	pwb = ParamDef(name="write_particles",vartype="boolean",desc_short="Write Particles",desc_long="Whether or not box images should be written",property=None,defaultunits=db.setdefault("write_particles", True),choices=None)
		# 	pwb.dependents = ["invert","normproc","format","suffix"] # these are things that become disabled when the pwb checkbox is unchecked etc
		# 	params.append([pwc, pwb])
		#
		# 	psuffix = ParamDef(name="suffix",vartype="string",desc_short="Output Suffix", desc_long="This text will be appended to the names of the output files",property=None,defaultunits=db.setdefault("suffix","_ptcls"),choices=None )
		# 	pinv = ParamDef(name="invert",vartype="boolean",desc_short="Invert Pixels",desc_long="Do you want the pixel intensities in the output inverted?",property=None,defaultunits=db.setdefault("invert",False),choices=None)
		# 	pn =  ParamDef(name="normproc",vartype="string",desc_short="Normalize Images",desc_long="How the output box images should be normalized",property=None,defaultunits=db.setdefault("normproc","normalize.edgemean"),choices=["normalize","normalize.edgemean","normalize.ramp.normvar","None"])
		# 	pop = ParamDef(name="format",vartype="string",desc_short="Output Image Format",desc_long="The format of the output box images",property=None,defaultunits=db.setdefault("format","bdb"),choices=self.output_formats)
		# 	params.append([psuffix,pinv])
		# 	params.append(pn)
		# 	params.append(pop)
		# else:
		# 	params.append(pwc)
		# return params
		params = []
		params.append(ParamDef(name="blurb",vartype="text",desc_short="",desc_long="",property=None,defaultunits=self.__doc__,choices=None))
		params.append(self.get_table())

		pbox = ParamDef(name="output_boxsize",vartype="int",desc_short="Box Size",desc_long="An integer value",property=None,defaultunits=db.setdefault("output_boxsize",self.dfl_boxsize),choices=[])
		pfo = ParamDef(name="force",vartype="boolean",desc_short="Force Overwrite",desc_long="Whether or not to force overwrite files that already exist",property=None,defaultunits=db.setdefault("force",False),choices=None)
		params.append([pbox,pfo])

		pwc = ParamDef(name="write_coords",vartype="boolean",desc_short="Write Coordinates",desc_long="Whether or not to write .box files",property=None,defaultunits=db.setdefault("write_coords", is_gauss),choices=None)
		pwb = ParamDef(name="write_particles",vartype="boolean",desc_short="Write Particles",desc_long="Whether or not box images should be written",property=None,defaultunits=db.setdefault("write_particles", True),choices=None)
		pwb.dependents = ["invert","normproc","format","suffix"] # these are things that become disabled when the pwb checkbox is unchecked etc
		params.append([pwc, pwb])

		psuffix = ParamDef(name="suffix",vartype="string",desc_short="Output Suffix", desc_long="This text will be appended to the names of the output files",property=None,defaultunits=db.setdefault("suffix","_ptcls"),choices=None )
		pinv = ParamDef(name="invert",vartype="boolean",desc_short="Invert Pixels",desc_long="Do you want the pixel intensities in the output inverted?",property=None,defaultunits=db.setdefault("invert",False),choices=None)
		pn =  ParamDef(name="normproc",vartype="string",desc_short="Normalize Images",desc_long="How the output box images should be normalized",property=None,defaultunits=db.setdefault("normproc","normalize.edgemean"),choices=["normalize","normalize.edgemean","normalize.ramp.normvar","None"])
		pop = ParamDef(name="format",vartype="string",desc_short="Output Image Format",desc_long="The format of the output box images",property=None,defaultunits=db.setdefault("format","hdf"),choices=self.output_formats)
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
		if params.get("write_particles"): particle_output_names = get_particle_outnames(params)

		coord_output_names = []
		if params["write_coords"]: coord_output_names =  get_coord_outnames(params)

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
			self.write_output(params["filenames"],particle_output_names,EMBoxList.write_particles,params["output_boxsize"],"Writing Particles", [params["invert"], params["normproc"]])
		if len(coord_output_names) > 0:
			self.write_output(params["filenames"],coord_output_names,EMBoxList.write_coordinates,params["output_boxsize"],"Writing Coordinates")

		from PyQt4 import QtCore
		self.emit(QtCore.SIGNAL("task_idle"))
		self.form.close()
		self.form = None
		#self.write_db_entries(params)

	def write_output(self,input_names,output_names,box_list_function,box_size,msg="Writing Output",extra_args=[]):
		n = len(input_names)
		from emapplication import EMProgressDialog
		progress = EMProgressDialog(msg, "Cancel", 0,n,None)
		progress.show()
		prog = 0

		files_written = []
		for i,output in enumerate(output_names):
			infile = input_names[i]
			box_list = EMBoxList()
			box_list.load_boxes_from_database(infile)
			if len(extra_args) > 0:box_list_function(box_list,infile,output,box_size,*extra_args)
			else:box_list_function(box_list,infile,output,box_size)
			files_written.append(output)
			prog += 1
			progress.setValue(prog)
			get_application().processEvents()

			if progress.wasCanceled():
				from EMAN2 import remove_file
				for file in files_written: remove_file(file)
				progress.setValue(nim)
				progress.close()
				return

		progress.setValue(n)
		progress.close()

def get_particle_outnames(params):
	input = params["filenames"]

	if get_platform() == 'Windows':
		input2 = []
		for name in input:
			name = name.replace('\\', '/')
			input2.append(name)
		input = input2

	format = params["format"]
	output = []
	for name in input:
		if format == "bdb":
			out = "bdb:particles#" + base_name(name) + params["suffix"]
		else:
			out = "particles/"+base_name(name)+ params["suffix"]+"."+format
		output.append(out)
	
	try: os.mkdir("particles")
	except: pass
	return output

def get_coord_outnames(params):
	input = params["filenames"]
	output = []
	for name in input:
		if params.get('format') and params['format'].lower() == 'json':
			output.append(base_name(name) + '.json')
		else:
			output.append(base_name(name)+ '.box')
	return output

from PyQt4 import QtGui
class EMBoxerInspector(QtGui.QWidget):

	PTCL_SHAPE_MAP = {}
	PTCL_SHAPE_MAP["none"] = "hidden"
	PTCL_SHAPE_MAP["square"] = "rect"
	PTCL_SHAPE_MAP["square with marker"] = "rectpoint"
	PTCL_SHAPE_MAP["circle"] = "rcircle"
	PTCL_SHAPE_MAP["circle with marker"] = "rcirclepoint"

	def __init__(self,target) :
		from PyQt4 import QtCore, QtGui
		self.busy = True
		self.tool_dynamic_vbl = None # this will be used to dynamic add widgets as the buttons are changed
		self.dynamic_box_button_widget = None # this will be used to dynamic add widgets as the buttons are changed
		self.ptcl_display_dict = None # this will be a dict mapping the names in the
		QtGui.QWidget.__init__(self,None)
		self.setWindowIcon(QtGui.QIcon(get_image_directory() +"green_boxes.png"))
		self.setWindowTitle("e2boxer")
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

		self.status_bar = QtGui.QStatusBar()
		self.vbl.addWidget(self.status_bar)
		self.status_bar.showMessage("Ready",10000)

		self.connect(self.status_bar,QtCore.SIGNAL("messageChanged(const QString&)"),self.on_status_msg_change)
		self.connect(self.done_but,QtCore.SIGNAL("clicked(bool)"),self.on_done)
		self.connect(self.gen_output_but,QtCore.SIGNAL("clicked(bool)"),self.write_output_clicked)
		self.busy = False

	def on_status_msg_change(self,s):
		if self.busy: return
		if len(s) == 0:
			if self.target() != None:
				self.target().load_default_status_msg()

	def set_status_message(self,mesg,timeout=1000):
		self.busy = True
		self.status_bar.showMessage(mesg,timeout)
		self.busy = False
	def on_done(self):
		self.close()

	def closeEvent(self,event):
		E2saveappwin("e2boxer","main",self)
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
		for val in EMBoxerInspector.PTCL_SHAPE_MAP.keys():
			self.boxformats.addItem(val)

		viewhbl2.addWidget(self.boxformats)

		# set the default
		db = js_open_dict(EMBOXERBASE_DB)
		val = db.setdefault("ptcl_display_shape","square")
		for n in range(self.boxformats.count()):
			if str(self.boxformats.itemText(n)) == val:
				self.boxformats.setCurrentIndex(n)
				break
		else:
			raise RuntimeError("Unknown ptcl display shape %s" %val)


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
		db = js_open_dict(EMBOXERBASE_DB)
		val = db["ptcl_display_shape"] = str(val)
		self.target().set_ptcl_display_shape(EMBoxerInspector.PTCL_SHAPE_MAP[val])


	def get_main_tab(self):
		from PyQt4 import QtCore, QtGui, Qt
		widget = QtGui.QWidget()
		vbl = QtGui.QVBoxLayout(widget)
		vbl.setMargin(0)
		vbl.setSpacing(6)

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

	def add_bottom_buttons(self,layout):
		from PyQt4 import QtCore, QtGui, Qt
		hbl_t=QtGui.QHBoxLayout()

		hbl_q=QtGui.QHBoxLayout()
		self.quality=QtGui.QLabel("Image Quality:")
		qual_tt = "Assign a quality number to the image. This acts as metadata for your convenience and is displayed in eman2 forms when possible."
		self.quality.setToolTip(qual_tt)
		hbl_q.addWidget(self.quality)

		self.image_qualities = QtGui.QComboBox()
		for i in range(10):
			self.image_qualities.addItem(str(i))
		self.image_qualities.setCurrentIndex(2)
		self.image_qualities.setToolTip(qual_tt)
		hbl_q.addWidget(self.image_qualities)
		layout.addLayout(hbl_q)

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
			raise RuntimeError("Unknown quality: " + str(val))
		self.busy = False

	def add_boxing_button_group(self,layout):
		from PyQt4 import QtCore, QtGui, Qt

		self.tool_button_group_box = QtGui.QGroupBox("Tools")
		self.tool_button_group_box_vbl = QtGui.QVBoxLayout(self.tool_button_group_box)
		self.tool_dynamic_vbl = QtGui.QVBoxLayout()

		hbl = QtGui.QHBoxLayout()
		current_tool_label = QtGui.QLabel("Current Boxing Tool:")
		self.current_tool_combobox = QtGui.QComboBox()
		hbl.addWidget(current_tool_label)
		hbl.addWidget(self.current_tool_combobox)

		self.tools_stacked_widget = QtGui.QStackedWidget()
		self.tool_dynamic_vbl.addLayout(hbl)
		self.tool_dynamic_vbl.addWidget(self.tools_stacked_widget)
		self.tool_button_group_box_vbl.addLayout(self.tool_dynamic_vbl,1)
		layout.addWidget(self.tool_button_group_box,0,)

		QtCore.QObject.connect(self.current_tool_combobox, QtCore.SIGNAL("activated(int)"), self.current_tool_combobox_changed)

	def add_mouse_tool(self,mouse_tool,):
#		icon = mouse_tool.icon()
#		if icon != None:
#			self.tool_tabs.addTab(mouse_tool.get_widget(),icon,mouse_tool.unique_name())
#		else:
#			self.tool_tabs.addTab(mouse_tool.get_widget(),mouse_tool.unique_name())
		self.tools_stacked_widget.addWidget(mouse_tool.get_widget())
		self.current_tool_combobox.addItem(mouse_tool.unique_name()) #Ross
		self.update()

	def set_tool_mode(self,name):
		self.busy = True
		for i in range(self.current_tool_combobox.count()):
			if str(self.current_tool_combobox.itemText(i)) == name:
				self.current_tool_combobox.setCurrentIndex(i)
				self.tools_stacked_widget.setCurrentIndex(i)
				break
		else:
			raise RuntimeError("Don't know the tab %s" %name)
		self.busy = False

	def current_tool_combobox_changed(self, idx):
		if self.busy: return
		self.tools_stacked_widget.setCurrentIndex(idx)
		self.target().set_main_2d_mouse_mode(str( self.current_tool_combobox.currentText() ))

	def set_box_size(self,value):
		self.busy = True
		self.box_size.setText(str(value))
		self.busy = False

	def new_box_size(self):
		if self.busy: return
		box_size=int(self.box_size.text())
		db = js_open_dict(EMBOXERBASE_DB)
		db["box_size"] = box_size
		if self.target(): self.target().set_box_size(box_size)

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
	my_main()




