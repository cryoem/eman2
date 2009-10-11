#!/usr/bin/env python
#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu) and David Woolford (woolford@bcm.edu)
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

# e2boxer.py  07/27/2004  Steven Ludtke
# This program is used to box out particles from micrographs/CCD frames

import PyQt4
from PyQt4 import QtCore, QtGui, Qt
from PyQt4.QtCore import Qt
from EMAN2 import *
from pyemtbx.boxertools import *
from optparse import OptionParser
from emshape import EMShape
from emimagemx import EMImageMXModule
from emimage2d import EMImage2DModule
from emimagerotor import EMImageRotorModule
from emimagemxrotor import EMImageMXRotorModule
from emrotor import EMRotorModule
from valslider import *
from math import *
from time import *
import os
import sys
import signal
from copy import *
#from OpenGL import contextdata
from emglplot import *
from emapplication import EMProgressDialogModule,get_application
# import SPARX definitions
import global_def
from global_def import *
# end import SPARX definitions


#from emglplot import *

from time import time,sleep

from sys import getrefcount

from emglobjects import EMOpenGLFlagsAndTools
from emapplication import EMStandAloneApplication,EMQtWidgetModule
import weakref

if os.name == 'nt':
	def kill(pid):
		"""kill function for Win32"""
		import win32api
		handle = win32api.OpenProcess(1, 0, pid)
		return (0 != win32api.TerminateProcess(handle, 0))
else:
	from os import kill

def get_out_file( f ):
	import os
	dir = os.path.dirname(f)
	base = os.path.basename( f )
	(name, ext) = os.path.splitext( base )
	outf = os.path.join( dir, ("particles_" + name + ".hdf") )
	return outf
pl=()


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <image>
	
Automatic and manual particle selection. This version is specifically aimed at square boxes
for single particle analysis."""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--gui",       action="store_true",help="Start the GUI for interactive boxing",default=False)
	parser.add_option("--boxsize","-B",type="int",help="Box size in pixels",default=-1)
	parser.add_option("--auto","-A",   type="string",action="append",help="Autobox using specified method: swarm, gauss",default=[])
	parser.add_option("--write_coord_files",action="store_true",help="Write data box files",default=False)
	parser.add_option("--write_box_images",action="store_true",help="Write data box files",default=False)
	parser.add_option("--force","-f",action="store_true",help="Force overwrites old files",default=False)
	parser.add_option("--overlap",type="int",help="(auto:grid) number of pixels of overlap between boxes. May be negative.")
	#parser.add_option("--nretest",type="int",help="(auto:ref) Number of reference images (starting with the first) to use in the final test for particle quality.",default=-1)
	#parser.add_option("--retestlist",type="string",help="(auto:ref) Comma separated list of image numbers for retest cycle",default="")
	#parser.add_option("--farfocus",type="string",help="filename or 'next', name of an aligned far from focus image for preliminary boxing",default=None)
	parser.add_option("--merge_boxes_to_db",action="store_true",help="A special argument, if true all input arguments are considered to be box files and they are merged into the project database as manually selected particles",default=False)
	parser.add_option("--subsample_method",help="The method used to subsample images prior to generation of the correlation image. Available methods are standard,careful",default="standard")	
	parser.add_option("--method",          help="boxer method, Swarm or Gauss", default="Gauss")
	parser.add_option("--outformat",       help="Format of the output particles images, should be bdb,img, or hdf", default="bdb")
	parser.add_option("--just_output",action="store_true", help="Applicable if doing auto boxing using the database. Bypasses autoboxing and just writes boxes that are currently stored in the database. Useful for changing the boxsize, for example", default=False)
	parser.add_option("--normproc",        help="Normalization processor to apply to particle images. Should be normalize, normalize.edgemean or none", default="normalize.edgemean")
	parser.add_option("--invert_output",action="store_true",help="If writing output only, this will invert the pixel intensities of the boxed files",default=False)
	parser.add_option("--dbls",            type="string",help="data base list storage, used by the workflow. You can ignore this argument.",default=None)
	
	# options added for cmdline calling of screening with gauss convolution. parameters for Gauss are passed as
	#    commandline arguments and will have to be parsed only if auto="cmd" option is specified. note that these parameters
	#    are applicable only to Gauss...
	# jl 10-28-08.
	parser.add_option("--ccf_lo",      type="float", default=-1.0,help="lower CCF threshold")
	parser.add_option("--ccf_hi",      type="float", default=1.0,help="upper CCF threshold")
	
	parser.add_option("--pix_in",      type="float", default=1.0,help="input pixel size in Angstrom")
	parser.add_option("--pix_out",     type="float", default=1.0,help="output pixel size in Angstrom")

	parser.add_option("--width",       default=False,help="width of the Gaussian")
	parser.add_option("--var",         action="store_true",default=False,help="Use variance flag (default is False)")
	parser.add_option("--inv",         action="store_true",default=False,help="Invert image flag (default is False)")

	parser.add_option("--do_ctf",      action="store_true",default=False,help="determine CTF (default is False)")
	parser.add_option("--ctf_fstart",  type="float", default=0.02,help="Start frequency for CTF determination")
	parser.add_option("--ctf_fstop",   type="float", default=0.5,help="Stop frequency for CTF determination")
	parser.add_option("--ctf_window",  type="int",   default=512,help="Window size for CTF determination")
	parser.add_option("--ctf_edge",    type="int",   default=0,help="Edge for CTF determination")
	parser.add_option("--ctf_overlap", type="int",   default=50,help="Overlap value for CTF determination")
	parser.add_option("--ctf_ampcont", type="float", default=10.0,help="Amplitude contrast for CTF")
	parser.add_option("--ctf_volt",    type="float", default=200.0,help="Voltage for CTF")
	parser.add_option("--ctf_cs",      type="float", default=2.0,help="Cs for CTF")
	parser.add_option("--out_file",    default=False,help="File to write particles to")
	parser.add_option("--out_dir",     default=False,help="Directory to write particle files to")

	global options
	(options, args) = parser.parse_args()
	if options.method == "Swarm":
		print "Note: Please consider switching to e2boxer2.py"

	filenames = []
	error_message = []
	for arg in args:
		if file_exists(arg):
			if options.merge_boxes_to_db == True:
				filenames.append(arg)
			else:

				nx,ny,nz = gimme_image_dimensions3D(arg)
				if nz != 1:
					error_message.append("%s is 3D!" %arg)
				elif EMUtil.get_image_count(arg) > 1:
					error_message.append("%s contains more than one image" %arg)
				else:
					filenames.append(arg)
		else:
			error_message.append("file %s does not exist" %arg)

	options.filenames = filenames
	
	if options.merge_boxes_to_db == True:
		#The user wants to add some boxes to the database
		merge_boxes_as_manual_to_db(filenames)
		sys.exit(1)
		
	if options.boxsize > 0 and  len(args) > 0:
		nx,ny = gimme_image_dimensions2D(arg)
		if options.boxsize > nx/2 or options.boxsize > ny/2:
			error_message.append("boxsize can not be greater than half of either of the image dimensions")


	if len(error_message) > 0:
		for error in error_message:
			print error

		sys.exit(1)

	if options.gui: options.running_mode = "gui"
	elif options.auto: options.running_mode = "auto_db"
	else:
		options.running_mode = None
		#print "unknown running mode"
		# in the new framework there is not need to tell the options parser there is an error, because the boxer module is smart enough to handle it

	logid=E2init(sys.argv)
	boxes=[]
	if len(options.auto)>0:
		if "cmd" in options.auto:
			print "Autobox mode ",options.auto[0]
			print "commandline version"
			
			do_gauss_cmd_line_boxing(options)
			print "cmdline autoboxer exiting"

			sys.exit(1)

		elif "grid" in options.auto:
			image_size=gimme_image_dimensions2D(filenames[0])
			try:
				dx=-options.overlap
				if dx+options.boxsize<=0 : dx=0.0
				dy=dx
			except:
				dy=(image_size[1]%options.boxsize)*options.boxsize/image_size[1]-1
				dx=(image_size[0]%options.boxsize)*options.boxsize/image_size[0]-1
				if dy<=0 : dy=((image_size[1]-1)%options.boxsize)*options.boxsize/image_size[1]-1
				if dx<=0 : dx=((image_size[0]-1)%options.boxsize)*options.boxsize/image_size[0]-1
			
	#		print image_size,dx,dy,options.boxsize
			for y in range(options.boxsize/2,image_size[1]-options.boxsize,dy+options.boxsize):
				for x in range(options.boxsize/2,image_size[0]-options.boxsize,dx+options.boxsize):
					boxes.append([x,y,options.boxsize,options.boxsize,0.0,1])
		else:
			# this is fine the boxer module should be smart enough to handle it
			pass


	if options.running_mode == "auto_db":
		dab = RawDatabaseAutoBoxer(logid)
		dab.go(options)
	else:
		application = EMStandAloneApplication()
		options.boxes = boxes
		options.logid = logid
		gui=EMBoxerModule(application,options)
		gui.show_guis()
	#	QtCore.QObject.connect(gui, QtCore.SIGNAL("module_idle"), on_idle)
		application.execute()

	E2end(logid)
	
	print "Exiting e2boxer"

def on_idle():
	# I may need this yet
	pass

def do_gauss_cmd_line_boxing(options):
	
	# commands to execute follow autbox_multi, except for getting parameters
	#    from the db
	project_db = EMProjectDB()

	# set up a parameter dict for passing arguments to the autoboxer object.
	#    dict keys will have to follow variable names in autoboxer.
	parm_dict = {}

	# parse cmd arguments
	if (options.ccf_lo):
		try:
			parm_dict["thr_low"] = float(options.ccf_lo)
		except ValueError:
			print "could not convert ccf_lo value. bad value",options.ccf_lo,". exiting!"
			sys.exit(1)

	if (options.ccf_hi):
		try:
			parm_dict["thr_hgh"] = float(options.ccf_hi)
		except ValueError:
			print "could not convert ccf_hi value. bad value",options.ccf_hi,". exiting!"
			sys.exit(1)

	if (options.pix_in):
		try:
			parm_dict["pixel_input"] = float(options.pix_in)
		except ValueError:
			print "could not convert pix_in value. bad value",options.pix_in,". exiting!"
			sys.exit(1)

	if (options.pix_out):
		try:
			parm_dict["pixel_output"] = float(options.pix_out)
		except ValueError:
			print "could not convert pix_out value. bad value",options.pix_out,". exiting!"
			sys.exit(1)

	if (options.width):
		try:
			parm_dict["gauss_width"] = float(options.width)
		except ValueError:
			print "could not convert gaussian width. bad value",options.width,". exiting!"
			sys.exit(1)

	if not( -1 == options.boxsize):
		try:
			parm_dict["box_size"] = int(options.boxsize)
		except ValueError:
			print "could not convert boxsize. bad value",options.width,". exiting!"
			sys.exit(1)
	# this is necessary, since default boxsize is -1 (which segfaults peak_ccf) and is not overwritten
	#    anywhere else....
	else:
		print "boxsize not set! exiting!"
		sys.exit(1)
		
	if (options.var):
		parm_dict["use_variance"] = True
	else:
		parm_dict["use_variance"] = False
		
	if (options.inv):
		parm_dict["invert"] = True
	else:
		parm_dict["invert"] = False


	# check ctf cmdline flag. if set, try to read f_start and f_stop parameters
	if (options.do_ctf):
		if (options.ctf_fstart):
			try:
				parm_dict["ctf_fstart"] = int(options.ctf_fstart)
			except ValueError:
				print "could not convert fstart value. bad value",options.ctf_fstart,". exiting!"
				sys.exit(1)
		if (options.ctf_fstop):
			try:
				parm_dict["ctf_fstop"] = int(options.ctf_fstop)
			except ValueError:
				print "could not convert fstop value. bad value",options.ctf_fstop,". exiting!"
				sys.exit(1)

		if (options.ctf_window):
			try:
				parm_dict["ctf_window"] = int(options.ctf_window)
			except ValueError:
				print "could not convert window value. bad value",options.ctf_window,". exiting!"
				sys.exit(1)
		if (options.ctf_edge):
			try:
				parm_dict["ctf_edge"] = int(options.ctf_edge)
			except ValueError:
				print "could not convert edge value. bad value",options.ctf_edge,". exiting!"
				sys.exit(1)
		if (options.ctf_overlap):
			try:
				parm_dict["ctf_overlap"] = int(options.ctf_overlap)
			except ValueError:
				print "could not convert overlap value. bad value",options.ctf_overlap,". exiting!"
				sys.exit(1)

		if (options.ctf_ampcont):
			try:
				parm_dict["ctf_ampcont"] = float(options.ctf_ampcont)
			except ValueError:
				print "could not convert amplitude contrast. bad value",options.ctf_ampcont,". exiting!"
				sys.exit(1)

		if (options.ctf_cs):
			try:
				parm_dict["ctf_Cs"] = float(options.ctf_cs)
			except ValueError:
				print "could not convert Cs factor. bad value",options.ctf_cs,". exiting!"
				sys.exit(1)

		if (options.ctf_volt):
			try:
				parm_dict["ctf_volt"] = int(options.ctf_volt)
			except ValueError:
				print "could not convert voltage value. bad value",options.ctf_volt,". exiting!"
				sys.exit(1)

	if (options.out_file):
		try:
			if ("bdb" == options.outformat):
				if not(options.out_file.startswith("bdb:")):
					print "outfile format does not match outformat! assuming correct outfile"
					parm_dict["out_file"] = options.out_file
				else:
					parm_dict["out_file"] = options.out_file
			elif ("hdf" == options.outformat):
				if not(options.out_file.endswith(".hdf")):
					print "outfile format does not match outformat! assuming correct outfile"
					parm_dict["out_file"] = options.out_file
				else:
					parm_dict["out_file"] = options.out_file
			else:
				parm_dict["out_file"] = options.out_file
				print "warning: may have unrecognized file format for output file!"
				
		except:
			print "could not set output file. exiting!"
			sys.exit(1)
		
	if (options.out_dir):
		try:
			import os.path, os
			
			if (os.path.lexists(os.path.abspath(options.out_dir))):
				parm_dict["out_dir"] = os.path.abspath(options.out_dir)+os.sep
			else:
				# make the dir. if this fails, we go to default cwd in the except clause
				os.mkdir(os.path.abspath(options.out_dir))
				parm_dict["out_dir"] = os.path.abspath(options.out_dir)+os.sep
			    
		except:
			print "could not set output directory. using cwd!"
			parm_dict["out_dir"] = os.getcwd() + os.sep

	# PawelAutoBoxer is changed to allow passing in of a parameter dictionary
	#    as additional argument....

	autoboxer = PawelAutoBoxer(None,parm_dict)

	# check whether to merge particles from several images into one file. this has to be
	#    done manually at the very end of the routine; to do this, set up a list of filenames
	#    and copy particles from those at the end. create an empty list here...
	if (options.out_file):
		image_list = []

	normalize=True
	norm_method="normalize.edgemean"

	for image_name in options.filenames:
		print "cmd autoboxing",image_name
		boxable = Boxable(image_name,None,autoboxer)
		
		if boxable.is_excluded():
			print "Image",image_name,"is excluded and being ignored"
			continue

		autoboxer.set_mode_explicit(SwarmAutoBoxer.COMMANDLINE)

		if (options.do_ctf):
			# new method to determine ctf...
			print "starting ctf determination"
			this_ctf = autoboxer.auto_ctf(boxable)
		else:
			# create empty so that del later doesn't raise exceptions
			this_ctf = None

		# Tell the boxer to delete non refs - FIXME - the uniform appraoch needs to occur - see SwarmAutoBoxer.auto_box
		print "starting autoboxer"
		autoboxer.auto_box(boxable,False)

		if options.write_coord_files:
			print "write coords"
			boxable.write_coord_file(box_size=-1,force=options.force,imageformat=options.outformat)
		if options.write_box_images:
			# we don't want to use boxable.write_box_images, since these don't store all information.
			#    do all this manually, then, but the code follows Boxable.write_box_images
			from utilities import set_params2D
			#img_name = boxable.get_image_file_name(options.outformat)
			if (options.out_dir):
				img_name = parm_dict["out_dir"] + boxable.get_image_file_name(options.outformat)
			else:
				from os import getcwd,mkdir,sep
				from os.path import abspath, lexists
				img_name = abspath(getcwd()) + sep + "particles" + sep + boxable.get_image_file_name(options.outformat)
				if not(lexists(getcwd()+sep+"particles")):
				       mkdir(getcwd()+sep+"particles")
				
			#print "img_name:",img_name
			if ("bdb" == options.outformat):
				img_name = boxable.get_image_file_name(options.outformat)
				if db_check_dict(img_name):
					if not(options.force):
						print "db",img_name,"already exists, doing nothing. Use force to override this behavior"
						return
					else:
						db_remove_dict(img_name)
						db_open_dict(img_name)
			elif ("hdf" == options.outformat):	
				if file_exists(img_name):
					if not(options.force):
						print "warning, file already exists - ", img_name, " doing nothing. Use force to override this behavior"
						return
					else:
						remove_file(img_name)
				
			print "writing",boxable.num_boxes(),"boxed images to", img_name
			for single_box in boxable.boxes:
				img = single_box.get_box_image(normalize,norm_method)
				# set all necessary attributes....
				#img.set_attr( "nx" , img.get_xsize())
				#img.set_attr( "ny" , img.get_ysize())
				#img.set_attr( "nz" , img.get_zsize())
				img.set_attr( "ctf" , this_ctf)
				img.set_attr( "Pixel_size", autoboxer.pixel_output )
				img.set_attr( "Micrograph", image_name )
				img.set_attr( "Score", single_box.correlation_score )
				img.set_attr( "ctf_applied" ,0 )
				img.set_attr( "active", 1 )
				set_params2D(img, [0.0, 0.0, 0.0, 0, 1.0])
				img.write_image(img_name,-1)

			if ("bdb" == options.outformat):
				db_close_dict(img_name)
			
			#print "write box images"
			#boxable.write_box_images(box_size=-1,force=options.force,imageformat=options.outformat)

		# check whether we will need to merge particles from a single image into
		#    one file. if so, append image name to a list
		if (options.out_file):
			# need to use get_image_file_name instead of raw name....
			image_list.append(boxable.get_image_file_name(options.outformat))

		del boxable,this_ctf
	
	# now check again for output file
	if (options.out_file):
		print "merging particles...."
		# check if target file exists already. if so, and if force is set, remove it
		if parm_dict["out_file"].startswith("bdb:"):
			# check for existing dict
			if (db_check_dict(parm_dict["out_file"]) and (options.force)):
				# remove and create new
				db_remove_dict(parm_dict["out_file"])
				db_open_dict(parm_dict["out_file"])
			# force isn't set, so we can't write....
			elif (db_check_dict(parm_dict["out_file"])):
				# set image list to empty, so that nothing is written but program
				#    exits normally...
				print "outfile exists. use --force to overwrite!"
				image_list = [] 

		elif parm_dict["out_file"].endswith(".hdf"):
			# check for existing file
			if (file_exists(parm_dict["out_file"]) and options.force):	
				# remove 
				remove_file(parm_dict["out_file"])
			# force isn't set, so we can't write....
			elif (file_exists(parm_dict["out_file"])):
				# set image list to empty, so that nothing is written but program
				#    exits normally...
				print "outfile exists. use --force to overwrite!"
				image_list = [] 
							
		image_object=EMData()
		for single_filename in image_list:
			# XXX: we loop over index to save space. reading the whole list may be faster, though...
			for image_index in xrange(EMUtil.get_image_count(single_filename)):
				image_object.read_image(single_filename,image_index)
				image_object.write_image(options.out_file,-1)

	project_db.close()
	#done

def merge_boxes_as_manual_to_db(filenames, use_progress=False):
	'''
	Merges a set of .box files into the local database - stores them as manual boxes
	'''
	for filename in filenames:
		f=file(filename,'r')
		lines=f.readlines()
		boxes = []
		for line in lines:
			data = str.split(line)
			b = Box(int(data[0]),int(data[1]),int(data[2]),int(data[3]),False)
			b.ismanual = True
			boxes.append(TrimBox(b))
	
		try:
			manualboxes = get_idd_key_entry(filename,"manual_boxes")
		except:
			manualboxes = []
	
		if manualboxes == None: manualboxes = []
		manualboxes.extend(boxes)
		set_idd_key_entry(filename,"manual_boxes",manualboxes)
		set_idd_key_entry(filename,"e2boxer_image_name",get_file_tag(filename))	
		
	#def process_start(self):
		#print "process started"
		

class EMBoxerModuleMouseEventsObject:
	'''
	A base class for objects that handle mouse events in the EMBoxerModule
	
	Inheriting objects are concerned with supplying their own definitions of 
	the functions mouse_up, mouse_down, mouse_drag, mouse_move and mouse_wheel. 
	They do not have to supply their own definition for all of these functions,
	but it would not make sense to inherit from this class unless the child class
	did not atleast supply one of them.
	
	Also stores a reference to the mediator object which coordinates the requests and signals
	of subclassing objects to the EMBoxerModule class. Modify the mediator class if you ever need
	to extend the scope of the messaging/signaling/requesting system between EMBoxerModuleMouseEventsObjects
	and the EMBoxerModule.
	
	The mediator object could be replaced with a direct reference to the EMBoxerModule class, but
	this is not recommended - for documentation and clarity of code purposes.
	'''
	def __init__(self,mediator):
		'''
		Stores only a reference to the mediator
		'''
		if not isinstance(mediator,EMBoxerModuleEventsMediator):
			print "error, the mediator should be a EMBoxerModuleEventsMediator"
			return
		
		self.mediator = mediator
		
	def get_2d_gui_image(self):
		'''
		Ask the mediator to for the EMImage2D object
		'''
		return self.mediator.get_2d_gui_image()
	
	def get_gui_ctl(self):
		'''
		Ask the mediator to for the main GUI controller
		'''
		return self.mediator.get_gui_ctl()
	
	def get_mx_gui_image(self):
		'''
		Ask the mediator to for the EMImageMX object
		'''
		return self.mediator.get_mx_gui_image()
	
	def mouse_up(self,event):
		'''
		Inheriting classes to potentially define this function
		'''
		pass

	def mouse_down(self,event):
		'''
		Inheriting classes to potentially define this function
		'''
		pass
		
	def mouse_drag(self,event):
		'''
		Inheriting classes to potentially define this function
		'''
		pass
	
	def mouse_move(self,event):
		'''
		Inheriting classes to potentially define this function
		'''
		pass

	def mouse_wheel(self,event):
		'''
		Inheriting classes to potentially define this function
		'''
		pass
	

class EMBoxerModuleMouseEraseEvents(EMBoxerModuleMouseEventsObject):
	'''
	A class that knows how to handle mouse erase events for a GUIBox
	'''
	def __init__(self,mediator,eraseradius=-1):
		EMBoxerModuleMouseEventsObject.__init__(self,mediator)
		
		self.eraseradius=eraseradius	# This is a circular radius 
		self.erasemode = None			# erase mode can be either Boxable.ERASE or Boxable.UNERASE
		
	def set_mode(self,mode):
		self.erasemode = mode
		
	def set_erase_radius(self,radius):
		self.eraseradius = radius
	
	def mouse_move(self,event):
		m = self.get_2d_gui_image().scr_to_img((event.x(),event.y()))
		self.get_2d_gui_image().add_eraser_shape("eraser",EMShape(["circle",.1,.1,.1,m[0],m[1],self.eraseradius,3]))
		self.mediator.update_image_display()
		
	def mouse_wheel(self,event):
		if event.modifiers()&Qt.ShiftModifier:
			self.get_gui_ctl().adjust_erase_rad(event.delta())
			m= self.get_2d_gui_image().scr_to_img((event.x(),event.y()))
			self.get_2d_gui_image().add_eraser_shape("eraser",EMShape(["circle",.1,.1,.1,m[0],m[1],self.eraseradius,3]))
			self.mediator.update_image_display()
	
	def mouse_down(self,event) :
		m=self.get_2d_gui_image().scr_to_img((event.x(),event.y()))
		#self.boxable.add_exclusion_area("circle",m[0],m[1],self.eraseradius)
		self.get_2d_gui_image().add_eraser_shape("eraser",EMShape(["circle",.9,.9,.9,m[0],m[1],self.eraseradius,3]))
		self.mediator.exclusion_area_added("circle",m[0],m[1],self.eraseradius,self.erasemode)	

	def mouse_drag(self,event) :
		m=self.get_2d_gui_image().scr_to_img((event.x(),event.y()))
		self.get_2d_gui_image().add_eraser_shape("eraser",EMShape(["circle",.9,.9,.9,m[0],m[1],self.eraseradius,3]))
		self.mediator.exclusion_area_added("circle",m[0],m[1],self.eraseradius,self.erasemode)
		# exclusion_area_added does the OpenGL update calls, so there is no need to do so here
		
	def mouse_up(self,event) :
		# we have finished erasing
		
		# make the eraser shape non visible
		self.get_2d_gui_image().add_eraser_shape("None",None)
		self.mediator.erasing_done(self.erasemode)
	
class EMBoxerModuleParticleManipEvents(EMBoxerModuleMouseEventsObject):
	'''
	A class that knows how to add, move and remove reference and non reference boxes 
	'''
	def __init__(self,mediator):
		EMBoxerModuleMouseEventsObject.__init__(self,mediator)
		self.mode =  EMBoxerModule.REFERENCE_ADDING
		self.moving = None
		self.dynapix = False

	def set_mode(self,mode):
		if mode not in [EMBoxerModule.REFERENCE_ADDING,EMBoxerModule.MANUALLY_ADDING]:
			print 'error, that is  an illegal mode'
			return
		
		self.mode = mode
		
	def mouse_down(self,event) :
		m = self.get_2d_gui_image().scr_to_img((event.x(),event.y()))
		box_num = self.mediator.detect_box_collision(m)
		if box_num == -1:
			#if not self.mediator.within_main_image_bounds(m):	return
			#if we make it here, that means the user has clicked on an area that is not in any box
			
			if event.modifiers()&Qt.ShiftModifier : return # the user tried to delete nothing
			
			# If we get here, we need to add a new reference
			box_size = self.mediator.get_box_size()
			
			box = Box(m[0]-box_size/2,m[1]-box_size/2,box_size,box_size,True)
			box.set_image_name(self.mediator.get_current_image_name())

			box.changed = True # this is so image2D nows to repaint the shape
			
			if self.mode == EMBoxerModule.REFERENCE_ADDING:
				box.isref = True
				box.ismanual = False
				
			elif self.mode == EMBoxerModule.MANUALLY_ADDING:
				box.isref = False
				box.ismanual = True
			else:
				print 'error, unknown error in mouse_down, boxing mode'
			
			x0=box.xcorner+box.xsize/2-1
			y0=box.ycorner+box.ysize/2-1
			self.get_2d_gui_image().add_shape("cen",EMShape([self.mediator.get_shape_string(),.9,.9,.4,x0,y0,x0+2,y0+2,1.0]))
			
			self.mediator.add_box(box)
			#self.moving=[box,m,box_num]
			#self.mediator.mouse_click_update_ppc()
		
		elif event.modifiers()&Qt.ShiftModifier :
			# remove the box
			self.mediator.remove_box(box_num)
			#self.mediator.mouse_click_update_ppc()
			
		else:
			# if we make it here than the we're moving a box
			box = self.mediator.get_box(box_num)
			self.moving=[box,m,box_num]
			self.get_2d_gui_image().set_active(box_num,.9,.9,.4)
				
			x0=box.xcorner+box.xsize/2-1
			y0=box.ycorner+box.ysize/2-1
			self.get_2d_gui_image().add_shape("cen",EMShape([self.mediator.get_shape_string(),.9,.9,.4,x0,y0,x0+2,y0+2,1.0]))
			object = self.get_mx_gui_image()
			if object.is_visible(box_num) or True : self.get_mx_gui_image().set_selected([box_num],True)
			self.mediator.update_all_image_displays()
			

	def mouse_drag(self,event) :
		
		m=self.get_2d_gui_image().scr_to_img((event.x(),event.y()))
		
		if event.modifiers()&Qt.ShiftModifier:
			box_num = self.mediator.detect_box_collision(m)
			if ( box_num != -1):
				self.mediator.remove_box(box_num)
				#self.mediator.mouse_click_update_ppc()
			
		elif self.moving != None:
			# self.moving[0] is the box, self.moving[1] are the mouse coordinates
			box = self.moving[0]
			# the old m in in self.moving[2]
			oldm = self.moving[1]
			
			self.mediator.move_box(self.moving[2],m[0]-oldm[0],m[1]-oldm[1])
			self.moving[1] = m
	
	def mouse_up(self,event) :
		if self.moving != None:
			box = self.moving[0]
			if box.isref:
				self.mediator.reference_moved(box)
			self.moving=None
		
		m = self.get_2d_gui_image().scr_to_img((event.x(),event.y()))
		box_num = self.mediator.detect_box_collision(m)
		if box_num != -1 and not event.modifiers()&Qt.ShiftModifier:
			object = self.get_mx_gui_image()
			if not object.is_visible(box_num) : object.scroll_to(box_num,True)
			self.get_mx_gui_image().set_selected([box_num],True)
			self.mediator.update_all_image_displays()

class EMBoxerModuleEventsMediator:
	'''
	This class could just as easily not exist - however it remains in use for documentation
	purposes. If it was removed, the EMBoxerModuleMouseEventsObjects could have a reference to the EMBoxerModule
	class instead of this one and the code would still work. But without this class it would be difficult to
	disentangle the relationship between the EMBoxerModuleMouseEventsObjects and the EMBoxerModule.
	
	This class
	supplies a bunch of public functions that accept the requests and signals from the EMBoxerModuleMouseEventsObject
	and send them on to the 'slots' and getter functions of the the EMBoxerModule - using basic function interfacing. This
	is also motivated by the Mediator concept in the Gang of Four.
	
	All things considered, the class remains for documentation purposes. It should only be removed if 
	it poses a significant performance hit
	
	
	NOVEMBER d.woolford says just use inheritance instead, then there isn't function call expenses but the documentation aspect remains
	'''
	def __init__(self,parent):
		'''
		Stores only a reference to the parent
		'''
		if not isinstance(parent,EMBoxerModule):
			print "error, the parent of a EMBoxerModuleMouseEraseEvents must be a EMBoxerModule"
			return
		
		self.parent = weakref.ref(parent)	# need a referene to the parent to send it events, what's more it needs to be week or we get cyclic regferences

	def get_2d_gui_image(self):
		'''
		Return the parent's EMImage2D object
		'''
		return self.parent().get_2d_gui_image()
	
	def get_gui_ctl(self):
		'''
		Return the parent's Controller widgit object
		'''
		return self.parent().get_gui_ctl()
	
	def get_mx_gui_image(self):
		'''
		Return the parent's EMImageMX object
		'''
		return self.parent().get_mx_gui_image()
	
	def update_image_display(self):
		'''
		Send an event to the parent that the EMImage2D should update its display
		'''
		self.parent().update_image_display()
		
	def update_all_image_displays(self):
		'''
		Send an event to the parent that the EMImage2D and EMImageMX objects should
		update their displays
		'''
		self.parent().update_all_image_displays()
		
	def exclusion_area_added(self,typeofexclusion,x,y,radius,mode):
		'''
		Send an event to the parent that an exclusion area was added.
		The parameters define the type of exclusion area. In future the exclusion area
		should probably just be its own class.
		'''
		self.parent().exclusion_area_added(typeofexclusion,x,y,radius,mode)

	def erasing_done(self,erase_mode):
		'''
		Send an event to the parent letting it know that the user has stopped adding
		erased area
		'''
		self.parent().erasing_done(erase_mode)
		
	def detect_box_collision(self,coords):
		'''
		Ask the parent to detect a collision between a resident box and the given coordinates.
		This is in terms of the EMImage2D 
		'''
		return self.parent().detect_box_collision(coords)
	
	def get_current_image_name(self):
		'''
		Ask the parent for the current name of the image in the large image view...
		'''
		return self.parent().get_current_image_name()

	def get_box_size(self):
		'''
		Ask the parent for the current project box_size
		'''
		return self.parent().get_box_size()
	
	def add_box(self,box):
		'''
		Tell the parent to store a box
		'''
		self.parent().add_box(box)
	
	def box_display_update(self):
		'''
		Tell the parent the a general box display update needs to occur
		'''
		self.parent().box_display_update()
		
	#def mouse_click_update_ppc(self):
		#'''
		#Tell the parent that a mouse click occured and that the PPC metric should be updated
		#'''
		#self.parent.mouse_click_update_ppc()
	
	def remove_box(self,box_num):
		'''
		Tell the parent to remove a box, as given by the box_num
		'''
		self.parent().remove_box(box_num)
		
	def get_box(self,box_num):
		'''
		Ask the parent for the box, as given by the box_num
		'''
		return self.parent().get_box(box_num)
	
	def move_box(self,box_num,dx,dy):
		'''
		Tell the parent to handle the movement of a box
		'''
		self.parent().move_box(box_num,dx,dy)
	
	def reference_moved(self,box):
		'''
		Tell the parent that a reference was moved - this could trigger automatic boxing
		'''
		self.parent().reference_moved(box)
		
	def within_main_image_bounds(self,coords):
		'''
		Ask the parent to determine if the coords are within the currently display image boundaries
		'''
		return self.parent().within_main_image_bounds(coords)
	
	def get_shape_string(self):
		'''
		Gets the shape string currently used for creating shapes for the 2D image
		'''
		return self.parent().get_shape_string()


class RawDatabaseAutoBoxer:
	def __init__(self,logid=None):
		self.required_options = ["boxsize","write_coord_files","write_box_images","force","normproc","outformat","just_output","invert_output"]
		self.logid = logid
		
	def go(self,options):
		options_ready = True
		for req_opt in self.required_options:
			if not hasattr(options,req_opt):
				if not self.application:
					print "there are insufficient parameters to run autoboxing"
					print "required option:", req_opt, "is missing"
					return
			
		if options_ready:
			self.autobox_images(options)
				
	def autobox_images(self,options):
		image_names = options.filenames
		project_db = EMProjectDB()
		
		
		for i,image_name in enumerate(image_names):
			
			try:
				data = project_db[get_idd_key(image_name)]
				
				trim_autoboxer = project_db[data["autoboxer_unique_id"]]["autoboxer"]
				autoboxer = SwarmAutoBoxer(None)
				autoboxer.become(trim_autoboxer)
				print 'using cached autoboxer db'
			except:
				try:
					print "using most recent autoboxer"
					if project_db["current_autoboxer_type"]=="Gauss":
						trim_autoboxer = project_db["current_autoboxer"]
						autoboxer = PawelAutoBoxer(None)
						autoboxer.become(trim_autoboxer)
					else:
						trim_autoboxer = project_db["current_autoboxer"]
						autoboxer = SwarmAutoBoxer(None)
						autoboxer.become(trim_autoboxer)
				except:
					print "Error - there seems to be no autoboxing information in the database - autobox interactively first - bailing"
					if self.logid:  E2progress(self.logid,1.0)
					return
			autoboxer.set_mode_explicit(SwarmAutoBoxer.COMMANDLINE)
		
			boxable = Boxable(image_name,None,autoboxer)
			boxable.set_autoboxer(autoboxer)
			if boxable.is_excluded():
				print "Image",image_name,"is excluded and being ignored"
				continue
			
			
			# Tell the boxer to delete non refs - FIXME - the uniform appraoch needs to occur - see SwarmAutoBoxer.auto_box
			if not options.just_output: # This is useful if you want to just change the boxsize, or the normalization method
				autoboxer.auto_box(boxable,False)
			if options.write_coord_files:
				boxable.write_coord_file(options.boxsize,options.force)
			if options.write_box_images:
				if options.normproc == "none":normalize=False
				else: normalize=True
				boxable.write_box_images(options.boxsize,options.force,imageformat=options.outformat,normalize=normalize,norm_method=options.normproc,invert=options.invert_output)
				if options.dbls != None and len(boxable.boxes)>0:
					
					from EMAN2 import get_file_tag
					db = db_open_dict("bdb:project")	
					particle_data =  db.get(options.dbls,dfl={})
					if isinstance(particle_data,list): # this is for back compatibility - in June 2009 we transitioned from lists to dicts
						d = {}
						for name in particle_data:
							s={}
							s["Original Data"] = name
							d[get_file_tag(name)] = s
						particle_data = d

					out_name = boxable.get_image_file_name(imageformat=options.outformat)
					s={}
					s["Original Data"] = out_name
					particle_data[get_file_tag(out_name)] = s
					db[options.dbls] = particle_data
		
			if self.logid:  E2progress(self.logid,float(i+1)/len(image_names))
		project_db.close()



def gen_thumbs(image_names,n):

	application = get_application()
	nim = len(image_names)
	thumbs = [None for i in range(nim)]
	progress = EMProgressDialogModule(application,"Generating thumbnails", "Abort", 0, nim,None)
	progress.qt_widget.show()
	prog = 0
	for i in range(nim):
#				thumb = self.get_image_thumb(i)
		thumb = get_idd_image_entry(image_names[i],"image_thumb")
		if thumb == None:
			
			#global BigImageCache
			thumb=EMData(image_names[i],0)
			#while n > 1:
				#image = image.process("math.meanshrink",{"n":2})
				#n /= 2
			
			thumb.process_inplace("math.meanshrink",{"n":n})
   	   	   	#thumb.process_inplace("testimage.noise.uniform.rand")
			thumb.process_inplace("normalize.edgemean") # if there are lots than they =should all have the same contrast
			set_idd_image_entry(image_names[i],"image_thumb",thumb)
#				image = None
#					print sys.getrefcount(image)
#				return thumb
   	   	#gc.collect()
   	   	#print "collected"
		prog += 1
		progress.qt_widget.setValue(prog)
		application.processEvents()
		#print "got thumb",i
		#thumbs[i] = thumb
			
		if progress.qt_widget.wasCanceled():
			progress.qt_widget.setValue(nim)
			progress.qt_widget.close()
			return -1 
		
	for i in range(nim):
		thumbs[i] = get_idd_image_entry(image_names[i],"image_thumb")
	
	progress.qt_widget.setValue(nim)
	progress.qt_widget.close()
	
	return thumbs
			

class DatabaseAutoBoxer(QtCore.QObject,RawDatabaseAutoBoxer):
	'''
	Initialize this with the application
	Then call the member function go (with the options)
	When this object is done it will emit "db_auto_boxing_done"
	'''
	def __init__(self,application,logid=None):
		QtCore.QObject.__init__(self)
		RawDatabaseAutoBoxer.__init__(self,logid)
		self.application = weakref.ref(application)
		

	def go(self,options):
		options_ready = True
		for req_opt in self.required_options:
			if not hasattr(options,req_opt):
				if not self.application:
					print "there are insufficient parameters to run autoboxing"
					print "required option:", req_opt, "is missing"
					return
				
				options_ready = False
				self.__run_form_initialization(options)
				return
			
		if options_ready:
			self.autobox_images(options)

	def __run_form_initialization(self,options):
		from emform import EMFormModule
		self.form = EMFormModule(self.get_params(options),get_application())
		self.form.setWindowTitle("Auto boxing parameters")
		get_application().show_specific(self.form)
		#print emitter
		QtCore.QObject.connect(self.form.emitter(),QtCore.SIGNAL("emform_ok"),self.on_form_ok)
		QtCore.QObject.connect(self.form.emitter(),QtCore.SIGNAL("emform_cancel"),self.on_form_cancel)

	def on_form_ok(self,params):
		options = EmptyObject()
		for key in params.keys():
			setattr(options,key,params[key])
			
		#print params
		#print options
		
		options_ready = True
		for req_opt in self.required_options:
			if not hasattr(options,req_opt):
				print "a parameters is missing"
				return
		
		get_application().close_specific(self.form)
		self.autobox_images(options)
		
	def on_form_cancel(self):
		get_application().close_specific(self.form)
	
	def get_params(self,options):
		from emdatastorage import ParamDef
		
		try: filenames = options.filenames
		except: filenames = []
		try: boxsize = options.boxsize
		except: boxsize = 128
		try: fo = options.force
		except: fo = False
		try: wc = options.write_coord_files
		except: wc = False
		try: wb = options.write_box_images
		except: wb = True
		try: norm = options.normproc
		except: norm = "normalize.edgemean"
		try: output = options.outformat
		except: output = "hdf"
		try: jo = options.just_output
		except: jo = False
	
		params = []
		params.append(ParamDef(name="filenames",vartype="url",desc_short="File names",desc_long="The files you wish to box",property=None,defaultunits=filenames,choices=[]))
		pbox = ParamDef(name="boxsize",vartype="int",desc_short="Box size",desc_long="The output box size",property=None,defaultunits=boxsize,choices=None)
		pfo = ParamDef(name="force",vartype="boolean",desc_short="Force overwrite",desc_long="Whether or not to force overwrite files that already exist",property=None,defaultunits=fo,choices=None)
		pjo = ParamDef(name="just_output",vartype="boolean",desc_short="Just output",desc_long="Bypass autoboxing and just use the boxes that are currently stored in the database",property=None,defaultunits=jo,choices=None)
		pwc = ParamDef(name="write_coord_files",vartype="boolean",desc_short="Write box db files",desc_long="Whether or not box db files should be written",property=None,defaultunits=wc,choices=None)
		pwb = ParamDef(name="write_box_images",vartype="boolean",desc_short="Write box image files",desc_long="Whether or not box images should be written",property=None,defaultunits=wb,choices=None)
		pn =  ParamDef(name="normproc",vartype="string",desc_short="Normalize images",desc_long="How the output box images should be normalized",property=None,defaultunits="normalize.edgmean",choices=["normalize","normalize.edgemean","none"])
		pop = ParamDef(name="outformat",vartype="string",desc_short="Output image format",desc_long="The format of the output box images",property=None,defaultunits="bdb",choices=["bdb","img","hdf"])
		params.append([pbox,pfo,pjo])
		params.append([pwc,pwb])
		params.append(pn)
		params.append(pop)
		
		return params

	def close(self):
		if self.form != None:
			get_application().close_specific(self.form)
			self.form = None

class EmptyObject:
	'''
	This just because I need an object I can assign attributes to, and object() doesn't seem to work
	'''
	def __init__(self):
		pass
	
# FIXME check for pyqt sometime
#try:
#	from PyQt4 import QtCore, QtGui, QtOpenGL
#	from PyQt4.QtCore import Qt
#	from valslider import ValSlider
#except:
#	print "Warning: PyQt4 must be installed to use the --gui option"
#	class dummy:
#		pass
#	class QWidget:
#		"A dummy class for use when Qt not installed"
#		def __init__(self,parent):
#			print "Qt4 has not been loaded"
#	QtGui=dummy()
#	QtGui.QWidget=QWidget
	
class EMBoxerModule(QtCore.QObject):
	'''
	'''
	REFERENCE_ADDING = 0
	ERASING = 1
	MANUALLY_ADDING = 2
	FANCY_MODE = 'fancy'
	PLAIN_MODE = 'plain'
	def __init__(self,application,options=None):
		QtCore.QObject.__init__(self)
		"""Implements the 'boxer' GUI."""
		
		self.required_options = ["filenames","running_mode"]
		
		self.boxable = None
		self.guictl_module = None
		self.guiim = None
		self.guimx = None
		self.guimxit = None
		self.dab = None
		self.form = None
		self.output_task = None # will be an EMAN2 style form
		
		form_init = False
		for s in self.required_options: 
			if not hasattr(options,s):
				form_init = True
				break

		if form_init or len(options.filenames) == 0 or options.running_mode not in ["gui","auto_db"]:
			self.__run_form_initialization(options)
		elif options.running_mode == "gui":
			self.__gui_init(options)
		elif options.running_mode == "auto_db":
			self.__auto_box_from_db(options)
		else: # if the first if statement is false then one of the second must be true
			print "this shouldn't happen"
	
	def show_guis(self):
		'''
		Generally called after initialization. 
		Will emit the module_idle signal if no guis/forms currently exist, this can be useful for closing the module.
		'''
		something_shown = False
		# guiim should be last so that it's at the front
		# guimx should come after guimit in the desktop, just because it's a better setup
		for gui in [self.form,self.guimxit,self.guimx,self.guictl_module,self.guiim]:
			if gui != None:
				something_shown = True
				get_application().show_specific(gui)
				
		if not something_shown:
			self.emit(QtCore.SIGNAL("module_idle"))
	
	def __run_form_initialization(self,options):
		from emform import EMFormModule
		self.form = EMFormModule(self.get_params(options),get_application())
		self.form.setWindowTitle("Boxer input variables")
		QtCore.QObject.connect(self.form.emitter(),QtCore.SIGNAL("emform_ok"),self.on_form_ok)
		QtCore.QObject.connect(self.form.emitter(),QtCore.SIGNAL("emform_cancel"),self.on_form_cancel)
		QtCore.QObject.connect(self.form.emitter(),QtCore.SIGNAL("emform_close"),self.on_form_cancel)
		
	def on_form_ok(self,params):
		options = EmptyObject()
		for key in params.keys():
			setattr(options,key,params[key])
		
		if len(options.filenames) == 0 or options.running_mode not in ["gui","auto_db"]:
			from emsprworkflow import error
			error("Please select the files you want to process","")
			return
		elif options.running_mode == "gui":
			self.__disconnect_form_signals()
			get_application().close_specific(self.form)
			self.form = None
			self.__gui_init(options)
			self.show_guis()
		elif options.running_mode == "auto_db":
			self.__disconnect_form_signals()
			get_application().close_specific(self.form)
			self.form = None
			self.__auto_box_from_db(options)
			
	def __disconnect_form_signals(self):
		QtCore.QObject.disconnect(self.form.emitter(),QtCore.SIGNAL("emform_ok"),self.on_form_ok)
		QtCore.QObject.disconnect(self.form.emitter(),QtCore.SIGNAL("emform_cancel"),self.on_form_cancel)
		QtCore.QObject.disconnect(self.form.emitter(),QtCore.SIGNAL("emform_close"),self.on_form_cancel)	
		
	def on_form_cancel(self):
		# this means e2boxer isn't doing anything. The application should probably be told to close the EMBoxerModule
#		self.__disconnect_form_signals()
		get_application().close_specific(self.form)
		self.form = None
		self.emit(QtCore.SIGNAL("module_idle"))
	
	def __auto_box_from_db(self,options):
		print "auto data base boxing"
		
		if not hasattr(options,"logid"): options.logid =None
		self.dab = DatabaseAutoBoxer(get_application(),options.logid)
		QtCore.QObject.connect(self.dab,QtCore.SIGNAL("db_auto_boxing_done"),self.on_db_autoboxing_done)
		
		self.dab.go(options)

	def on_db_autoboxing_done(self):
		self.emit(QtCore.SIGNAL("module_idle"))
	
	def __gui_init(self,options):
		
		if not hasattr(options,"boxes"): options.boxes = [] # this is a temporary workaround
		self.__initialize_from_parms(options.filenames,options.boxes,options.boxsize,options.method)
		
		self.__init_guiim() # initialise the 2D image display
		#self.__init_guimx() # intialize the matrix display
		self.__init_guimx_thumbs() # initialize the thumbnail diplsy
		self.ab_sel_mediator = AutoBoxerSelectionsMediator(self)
		self.__init_guictl()
		if self.fancy_mode == EMBoxerModule.FANCY_MODE:
			self.__init_ctl_rotor()

		if isinstance(self.autoboxer,SwarmAutoBoxer):
			self.autoboxer.auto_box(self.boxable,False) # Do the automatic autoboxing - this makes the user see results immediately
		self.box_display_update(stop_animation=True) # update displays to show boxes etc

	def get_params(self,options=None):
		'''
		Gets the params in order to run the Form module. This happens when no filenames have been specified, for example. Only works for the Swarm mode, because the other mode has broken with conventions etc.
		'''

		if options == None:
			filenames = None
			box_size = 128
			running_mode = "gui"
			mode = options.method #"Swarm"
		else:
			filenames = options.filenames
			box_size = options.boxsize
			if box_size == -1: box_size = 128
			running_mode = options.running_mode
			if running_mode == None:
				running_mode = "gui"
			mode = options.method #"Swarm"

		from emdatastorage import ParamDef
		from emsprworkflow import ParticleWorkFlowTask,EMRawDataReportTask
		from emform import EM2DFileTable,EMFileTable
		table = EM2DFileTable([],desc_short="Input Files",desc_long="")
		context_menu_data = ParticleWorkFlowTask.DataContextMenu()
		table.add_context_menu_data(context_menu_data)
		table.add_button_data(ParticleWorkFlowTask.AddDataButton(table,context_menu_data))
		table.add_column_data(EMFileTable.EMColumnData("Dimensions",EMRawDataReportTask.get_image_dimensions,"The dimensions of the file on disk"))

		params = []
		#params.append(ParamDef(name="filenames",vartype="url",desc_short="File names",desc_long="The files you wish to box",property=None,defaultunits=filenames,choices=[]))
		params.append(table)
		params.append(ParamDef(name="boxsize",vartype="int",desc_short="Box size",desc_long="An integer value",property=None,defaultunits=box_size,choices=[]))
		params.append(ParamDef(name="method",vartype="choice",desc_short="Boxing mode",desc_long="Currently only one mode is supported, but this could change",property=None,defaultunits=mode,choices=[options.method]))
		params.append(ParamDef(name="running_mode",vartype="choice",desc_short="Boxing mode",desc_long="Whether to load the GUI or run automatic boxing based on information stored in the database",property=None,defaultunits=running_mode,choices=["gui","auto_db"]))
		return params

	def __initialize_from_parms(self,image_names,boxes=[],box_size=-1, default_method="Swarm"):
		
		# initialize important autoboxer related variables
		self.dynapix = False
		self.image_names = image_names
		self.current_image_idx = 0
		self.box_size = box_size
		self.fancy_mode = EMBoxerModule.PLAIN_MODE # for now just make sure the fancy stuff isn't used
		# set self.autoboxer
		if len(self.image_names) != 0:
			self.set_autoboxer(self.image_names[self.current_image_idx], default_method)

		self.eraseradius = 2*self.box_size # this happens after the autoboxer has been loaded, because the boxsize can change
		self.erasemode = None #stores the erase mode
		self.shape_string = "rectpoint" # the shape of the picked particles
		self.in_display_limbo = False	# a flag I am using to solve a problem

		# A boxable is just a class that manages boxes in terms of images
		if len(self.image_names) != 0:
			#print self.autoboxer
			self.boxable = Boxable(self.image_names[0],self,self.autoboxer)
			self.boxable.add_non_refs(boxes)
			self.boxable.set_box_size(self.box_size)
		else: self.boxable = None
		
		self.initialize_mouse_event_handlers() # initialize the mouse event handlers

		self.ptcl=[] # list of actual boxed out EMImages. This may be redundant I am working on a better solution.
		self.moving_box_data = None # a vector storing [mouse x, mouse y, box idx]
		self.moving=None # Used during a user box drag. May be redudant could potentially just use self.moving_box_data. FIXME
	
		self.ab_sel_mediator = AutoBoxerSelectionsMediator(self)
		self.guimx = None # widget for displaying matrix of smaller imagespaugay

	def __init_ctl_rotor(self):
		self.ctl_rotor = EMRotor()
		self.ctl_rotor.get_core_object().add_qt_widget(self.guictl)
		self.guictlrotor = EMParentWin(self.ctl_rotor)
		self.guictlrotor.setWindowTitle("e2boxer Controllers")

	def __init_guictl(self):
		self.guictl_module = EMBoxerModulePanelModule(get_application(),self,self.ab_sel_mediator)
		self.guictl = self.guictl_module.qt_widget
		self.guictl.set_image_quality(self.boxable.get_quality())
		self.guictl.setWindowTitle("e2boxer Controller")
		self.guictl.set_dynapix(self.dynapix)
		#if self.fancy_mode == EMBoxerModule.FANCY_MODE: self.guictl.hide()
		#get_application().show_specific(self.guictl_module)
		if isinstance(self.autoboxer,PawelAutoBoxer):
			#print "Setting GUI for Gauss boxing method"
			gauss_method_id = 1
			self.guictl.method.setCurrentIndex(gauss_method_id)
			self.guictl.method_changed( gauss_method_id ) 
			self.autoboxer.set_params_of_gui(self.boxable)

	def __init_guimx_thumbs(self):
		self.itshrink = -1 # image thumb shrink. Default value of -1 means it has to be calculated when it's first needed
		self.imagethumbs = None # image thumbs - will be a list of tiny images
		self.guimxitp = None
		self.guimxit = None
		self.__gen_image_thumbnails_widget()

		if self.guimxit != None:
			if isinstance(self.guimxit,EMImageRotorModule):
				self.guimxit.optimally_resize()
				QtCore.QObject.connect(self.guimxit.emitter(),QtCore.SIGNAL("image_selected"),self.image_selected)
			else:
				QtCore.QObject.connect(self.guimxit.emitter(),QtCore.SIGNAL("mx_mouseup"),self.image_selected)
				QtCore.QObject.connect(self.guimxit.emitter(),QtCore.SIGNAL("module_closed"),self.guimxit_closed)
		
			if isinstance(self.guimxit,EMImageRotorModule):
				self.guimxit.set_frozen(self.boxable.is_frozen(),self.current_image_idx)

	def guimxit_closed(self):
		self.guimxit = None
		self.guictl_module.widget.set_guimxit_visible(False)

	def __init_guimx(self):
		glflags = EMOpenGLFlagsAndTools()
		emftgl_supported = True
		try: a = EMFTGL()
		except: emftgl_supported = False
		#if not glflags.npt_textures_unsupported() and emftgl_supported:
			#self.guimx=EMImageMXRotorModule(application=self.application)# widget for displaying image thumbs
			#self.guimx.disable_mx_zoom()
			##self.guimx.allow_camera_rotations(False)
			##self.guimx.disable_mx_translate()
			##self.fancy_mode = EMBoxerModule.FANCY_MODE
			
		#else:
		self.guimx=EMImageMXModule(application=get_application())
		self.guimx.desktop_hint = "rotor"
		self.fancy_mode = EMBoxerModule.PLAIN_MODE
		
		self.guimx.set_mouse_mode("App")
		
		qt_target = get_application().get_qt_emitter(self.guimx)

		#self.guimx.connect(self.guimx,QtCore.SIGNAL("removeshape"),self.removeshape)
		QtCore.QObject.connect(self.guimx.emitter(),QtCore.SIGNAL("mx_image_selected"),self.box_selected)
		QtCore.QObject.connect(self.guimx.emitter(),QtCore.SIGNAL("mx_mousedrag"),self.box_moved)
		QtCore.QObject.connect(self.guimx.emitter(),QtCore.SIGNAL("mx_mouseup"),self.box_released)
		QtCore.QObject.connect(self.guimx.emitter(),QtCore.SIGNAL("mx_boxdeleted"),self.box_image_deleted)
		QtCore.QObject.connect(self.guimx.emitter(),QtCore.SIGNAL("module_closed"),self.guimx_closed)
		if self.fancy_mode == EMBoxerModule.FANCY_MODE:
			QtCore.QObject.connect(self.guimx.emitter(),QtCore.SIGNAL("inspector_shown"),self.guimx_inspector_requested)
	
	def guimx_closed(self):
		self.guimx = None
		self.guictl_module.widget.set_guimx_visible(False)
	
	def __init_guiim(self, image=None, imagename=None):
		if image == None:
			imagename = self.image_names[self.current_image_idx]
			global BigImageCache
			image=BigImageCache.get_object(imagename).get_image(use_alternate=True)
		
		self.guiim= EMImage2DModule(application=get_application())
		self.guiim.set_data(image,imagename)
		self.guiim.force_display_update()
		#get_application().show_specific(self.guiim)

		self.__update_guiim_states()
		self.guiim.set_mouse_mode(0)
		
		self.guiim.setWindowTitle(imagename)
		
		QtCore.QObject.connect(self.guiim.emitter(),QtCore.SIGNAL("mousedown"),self.mouse_down)
		QtCore.QObject.connect(self.guiim.emitter(),QtCore.SIGNAL("mousedrag"),self.mouse_drag)
		QtCore.QObject.connect(self.guiim.emitter(),QtCore.SIGNAL("mouseup")  ,self.mouse_up  )
		QtCore.QObject.connect(self.guiim.emitter(),QtCore.SIGNAL("keypress"),self.keypress)
		QtCore.QObject.connect(self.guiim.emitter(),QtCore.SIGNAL("mousewheel"),self.mouse_wheel)
		QtCore.QObject.connect(self.guiim.emitter(),QtCore.SIGNAL("mousemove"),self.mouse_move)
		QtCore.QObject.connect(self.guiim.emitter(),QtCore.SIGNAL("module_closed"),self.guiim_closed)
	
	def guiim_closed(self):
		self.guiim = None
		self.guictl_module.widget.set_guiim_visible(False)

	def __update_guiim_states(self):
		if self.guiim == None:
			current_name = self.image_names[self.current_image_idx]
			global BigImageCache
			self.__init_guiim(BigImageCache.get_image_directly(current_name),current_name)
			
		self.guiim.set_other_data(self.boxable.get_exclusion_image(False),self.autoboxer.get_subsample_rate(),True)
		
		self.guiim.set_frozen(self.boxable.is_frozen())
		self.guiim.set_excluded(self.boxable.is_excluded())
		self.guiim.set_file_name(self.image_names[self.current_image_idx])
		
	def image_selected(self,event,lc):
		get_application().setOverrideCursor(Qt.BusyCursor)
		#try:
		im=lc[0]
		#try:
		debug = False
		if im != self.current_image_idx:
			
			global BigImageCache
			image=BigImageCache.get_image_directly(self.image_names[im])
			
			if self.guiim == None:
				current_name = self.image_names[self.current_image_idx]
				self.__init_guiim(image,current_name)

			
			try: 
				self.guiim.setWindowTitle(self.image_names[im])
			except:
				print "set window title failed"
				
			self.guiim.set_data(image)
			self.boxable.cache_exc_to_db()
			self.boxable = Boxable(self.image_names[im],self,self.autoboxer)
			self.ptcl = []
			self.guiim.del_shapes()
			self.guiim.force_display_update()
			self.in_display_limbo = True
			project_db = EMProjectDB()
			data = project_db[get_idd_key(self.image_names[im])]
			ab_failure = self.autoboxer
			
			if data != None:
				try:
					autoboxer_id = data["autoboxer_unique_id"]
					trim_autoboxer = project_db[autoboxer_id]["autoboxer"]
					self.autoboxer_name = autoboxer_id
					self.autoboxer = SwarmAutoBoxer(self)
					self.autoboxer.become(trim_autoboxer)
					self.dynapix = self.autoboxer.dynapix_on()
					self.guictl.set_dynapix(self.dynapix)
				except:
					try:
						trim_autoboxer = project_db["current_autoboxer"]
						self.autoboxer = SwarmAutoBoxer(self)
						self.autoboxer.become(trim_autoboxer)
						self.autoboxer_name = self.autoboxer.get_unique_stamp()
						self.dynapix = self.autoboxer.dynapix_on()
						self.guictl.set_dynapix(self.dynapix)
					except:
						self.autoboxer = ab_failure
			
			if self.dynapix:
				self.autoboxer.auto_box(self.boxable,False)
				
			self.boxable.set_autoboxer(self.autoboxer)
			
			self.autoboxer_db_changed()
			
			if self.box_size != self.autoboxer.get_box_size():
				self.update_box_size(self.autoboxer.get_box_size())

			self.in_display_limbo = False
			
			for box in self.boxable.boxes: box.changed = True
			
			self.current_image_idx = im
			
			self.__update_guiim_states()
			if isinstance(self.guimxit,EMImageRotorModule):
				self.guimxit.set_frozen(self.boxable.is_frozen(),self.current_image_idx)
			self.guictl.set_image_quality(self.boxable.get_quality())
			
			self.box_display_update()
		#except: pass

		get_application().setOverrideCursor(Qt.ArrowCursor)
	
	def set_autoboxer(self,imagename, default_method="Swarm"):
		'''
		This function sets the self.autoboxer variable. There are 3 cases:
		
		1. The image already has an autoboxer in the database which becomes self.autoboxer
		2. There is a current autoboxer in the database in which case this becomes self.autoboxer
		3. Neither 1 nor 2 - then an empty autoboxer is created
		
		'''
		try:
			project_db = EMProjectDB()
			data = project_db[get_idd_key(imagename)]
			autoboxer_id = data["autoboxer_unique_id"]
			type_autoboxer = project_db[autoboxer_id]["autoboxer_type"]
			trim_autoboxer = project_db[autoboxer_id]["autoboxer"]
			self.autoboxer_name = autoboxer_id

			if type_autoboxer=="Swarm":
				self.autoboxer = SwarmAutoBoxer(self)
				self.autoboxer.become(trim_autoboxer)
				self.dynapix = self.autoboxer.dynapix_on()
			else:
				self.autoboxer = PawelAutoBoxer(self)
				self.autoboxer.become(trim_autoboxer)
				self.autoboxer.source = "loaded"
			if self.box_size==-1: self.box_size = self.autoboxer.get_box_size()
			else: self.autoboxer.set_box_size(self.box_size,self.image_names)
		except Exception, inst:
#			print "exception happen when load image's autoboxer: ", inst
			try:
				type_autoboxer = project_db["current_autoboxer_type"]
				trim_autoboxer = project_db["current_autoboxer"]

				#print "Current Autoboxer Type: ", type_autoboxer

				if type_autoboxer=="Swarm":
					self.autoboxer = SwarmAutoBoxer(self)
					self.autoboxer.become(trim_autoboxer)
				elif type_autoboxer=="Gauss":
					self.autoboxer = PawelAutoBoxer(self)
					self.autoboxer.become(trim_autoboxer)
					self.autoboxer.source = "loaded"
				else:
					raise Exception("Current autoboxer not set")

				self.autoboxer_name = self.autoboxer.get_unique_stamp()
				self.dynapix = self.autoboxer.dynapix_on()
				if self.box_size==-1: self.box_size = self.autoboxer.get_box_size()
				else: self.autoboxer.set_box_size(self.box_size,self.image_names)
				
			except Exception, inst:
				#print "exception happend when load current autoboxer: ", inst

				
				if self.box_size == -1: self.box_size = 128
				
				if default_method=="Swarm":
					self.autoboxer = SwarmAutoBoxer(self)
					self.autoboxer.set_box_size_explicit(self.box_size)
					self.autoboxer.set_interactive_mode(self.dynapix)
				else:
					self.autoboxer = PawelAutoBoxer(self)
					self.autoboxer.source = "new"
					
				self.autoboxer.set_box_size(self.box_size,self.image_names)

	def guiim_inspector_requested(self,event):
		self.ctl_rotor.get_core_object().add_qt_widget(self.guiim.get_core_object().get_inspector())
	
	def guimx_inspector_requested(self,event):
		if self.guimx == None: self.__init_guiimx()
		self.ctl_rotor.get_core_object().add_qt_widget(self.guimx.get_inspector())
	
	def initialize_mouse_event_handlers(self):
		self.eventsmediator = EMBoxerModuleEventsMediator(self)
		self.mouse_handlers = {}
		self.mouse_handlers["boxing"] = EMBoxerModuleParticleManipEvents(self.eventsmediator)
		self.mouse_handlers["erasing"] = EMBoxerModuleMouseEraseEvents(self.eventsmediator,self.eraseradius)
		self.mousehandler = self.mouse_handlers["boxing"]
	
	def get_shape_string(self):
		return self.shape_string
	
	def has_thumbnails(self):
		return self.guimxit != None
	
	def view_boxes_clicked(self,bool):
		if self.guimx == None: self.__init_guiimx()
		if bool:
			get_application().show_specific(self.guimx)
		else:
			get_application().hide_specific(self.guimx)
			
	def view_image_clicked(self,bool):
		if self.guiim == None:
			current_name = self.image_names[self.current_image_idx]
			global BigImageCache
			self.__init_guiim(BigImageCache.get_image_directly(current_name),current_name)

		if bool:
			get_application().show_specific(self.guiim)
		else:
			get_application().hide_specific(self.guiim)
			
	def view_thumbs_clicked(self,bool):
		if  self.guimxit== None: return
		if bool:
			get_application().show_specific(self.guimxit)
		else:
			get_application().hide_specific(self.guimxit)
			
	def within_main_image_bounds(self,coords):
		x = coords[0]
		y = coords[1]
		if x < 0 or y < 0:
			return False
		
		main_image = self.get_current_image()
		
		if x >= main_image.get_xsize() or y >= main_image.get_ysize():
			return False
		
		return True
		
	def change_current_autoboxer(self, autoboxer_id,autobox=True):
		#print "change current autoboxer"
		self.autoboxer.write_to_db()
		project_db = EMProjectDB()
		if project_db[autoboxer_id] == None:
			project_db[autoboxer_id] = {}
		trim_autoboxer = project_db[autoboxer_id]["autoboxer"]
		#print "changing autoboxer to autoboxer_",timestamp,"and its stamp in the db is",trim_autoboxer.get_creation_ts()
		self.autoboxer = SwarmAutoBoxer(self)
		self.autoboxer.become(trim_autoboxer)
		self.autoboxer.write_image_specific_references_to_db(self.boxable.get_image_name())
		self.dynapix = self.autoboxer.dynapix_on()
		self.guictl.set_dynapix(self.dynapix)
		
		self.boxable.set_autoboxer(self.autoboxer)
		
		if self.box_size != self.autoboxer.get_box_size():
			self.update_box_size(self.autoboxer.get_box_size())
			#print "box display update"
		
		if autobox:
			self.ptcl = []
			self.guiim.del_shapes()
			self.boxable.clear_and_reload_images()
			self.in_display_limbo = True
			#self.autoboxer.regressiveflag = True
			self.autoboxer.auto_box(self.boxable)
			self.in_display_limbo = False
		
		self.box_display_update()	
		
	def update_box_size(self,box_size,mode=0):
		if box_size != self.box_size:
			if mode == 0:
				self.box_size = box_size
				self.boxable.change_box_size(self.box_size)
				self.guictl.set_box_size(self.box_size)
				self.boxable.get_exclusion_image(True)
				self.boxable.reload_boxes() # this may be inefficient
				#print "clearing displays"
				self.guiim.set_other_data(self.boxable.get_exclusion_image(False),self.autoboxer.get_subsample_rate(),True)
				self.clear_displays()
			elif mode == 1:
				self.box_size = box_size
				self.autoboxer.set_box_size(self.box_size,self.image_names)
				self.boxable.reload_boxes() # this may be inefficient
				#print "clearing displays"
				self.guiim.set_other_data(self.boxable.get_exclusion_image(False),self.autoboxer.get_subsample_rate(),True)
				self.clear_displays()	
			else:
				print "error, unknown mode in update_box_size"
		
	def get_2d_gui_image(self):
		if self.guiim == None:
			current_name = self.image_names[self.current_image_idx]
			global BigImageCache
			self.__init_guiim(BigImageCache.get_image_directly(current_name),current_name)
		return self.guiim
	
	def get_gui_ctl(self):
		if self.guictl == None: self.__init_guictl()
		return self.guictl
	
	def get_mx_gui_image(self):
		if self.guimx == None: self.__init_guimx()
		return self.guimx
	
	def get_boxable(self):
		return self.boxable
		
	def exclusion_area_added(self,typeofexclusion,x,y,radius,mode):
		self.boxable.add_exclusion_area(typeofexclusion,x,y,radius,mode)
		self.guiim.set_other_data(self.boxable.get_exclusion_image(False),self.autoboxer.get_subsample_rate(),True)
		self.update_image_display()
		
	def erasing_done(self,erase_mode):
		'''
		Call this function after erasing has occured to remove all boxes in the
		erased regions
		'''
		
		if erase_mode == Boxable.ERASE:
			# tell the boxable to remove boxes (and return their number)
			[lostboxes,refboxes] = self.boxable.update_excluded_boxes()
			# after this, make sure the display is correct.
			
			# If any of the boxes are references the autoboxer needs to know about it...
			# this could potentially trigger auto boxer 
			if len(lostboxes) != 0:
				self.delete_display_boxes(lostboxes)
				#self.mouseclicks += 1
				#self.update_ppc()
				self.update_all_image_displays()
				
			else: self.update_image_display()
			
			if len(refboxes) != 0: 
				val = self.autoboxer.remove_reference(refboxes)
				if val == 2:
					self.boxable.clear_and_cache(True)
					self.clear_displays()
					return # avoid unecessary display update below
			
			self.box_display_update()
			
		elif erase_mode == Boxable.UNERASE:
			[added_boxes,added_refboxes] = self.boxable.update_included_boxes()
			
			# len(added_boxes) is always <= len(added_refboxes)
			if len(added_boxes) != 0:
				self.box_display_update()
				
				if len(added_refboxes) != 0:
					val = self.autoboxer.add_reference(add_refboxes)

				
	def detect_box_collision(self,coords):
		'''
		Detects a collision of the coordinates with any of the boxes in the current image
		stored in guiim. coords need to be in the image coordinates
		'''
		return self.collision_detect(coords,self.get_boxes())
	
	def get_current_autoboxer_ts(self):
		return self.autoboxer.get_creation_ts()
	
	def get_current_image_name(self):
		if len(self.image_names) != 0:
			return self.image_names[self.current_image_idx]
		else: return None
	
	def get_current_image(self):
		global BigImageCache
		return BigImageCache.get_image_directly(self.get_current_image_name())
	
	def get_image_names(self):
		return self.image_names
	
	def get_box_size(self):
		return self.box_size
	
	def add_box(self,box):
		if not self.boxable.is_interactive():
			return
		
		self.boxable.add_box(box)
		
		# Should this be here?
		box_num = len(self.get_boxes())
		#self.guimx.get_core_object().set_selected([box_num])
		
		# autoboxer will autobox depending on the state of its mode
		if box.isref : self.autoboxer.add_reference(box)
		
		self.box_display_update()
		
	def get_box(self,box_num):
		boxes = self.get_boxes()
		return boxes[box_num]
		
	def remove_boxes(self,box_nums):
		if not self.boxable.is_interactive(): return # what does this do again?
		
		for box_num in box_nums:
			box = self.delete_box(box_num)
		if not (box.isref or box.ismanual):
			self.boxable.add_exclusion_particle(box)
	
		self.guiim.set_other_data(self.boxable.get_exclusion_image(False),self.autoboxer.get_subsample_rate(),True)
		self.box_display_update()
	
	def remove_box(self,box_num):
		if not self.boxable.is_interactive(): return
		
		box = self.delete_box(box_num)
		if not (box.isref or box.ismanual):
			self.boxable.add_exclusion_particle(box)
		self.guiim.set_other_data(self.boxable.get_exclusion_image(False),self.autoboxer.get_subsample_rate(),True)
		self.box_display_update()
		
		return box
	def reference_moved(self,box):
		if (self.autoboxer.reference_moved(box)):
			self.box_display_update()
		
	#def mouse_click_update_ppc(self):
		#self.mouseclicks += 1
		#self.update_ppc()
	
	def set_ptcl_mx_data(self,data=None):
		'''
		Call this to set the Ptcl Mx data 
		'''
		if data != None:
			if self.guimx == None: self.__init_guimx()
			if len(data) != 0:
				#get_application().show_specific(self.guimx) #NEED TO CHANGE THIS IN THE MIDDLE OF REFURBISHMENT
				self.guimx.set_data(data)
			
				#qt_target = self.guimx.get_parent()
				#if not qt_target.isVisible():
					#get_application().show_specific(self.guimx)
	
	def clear_displays(self):
		self.ptcl = []
		if self.guiim != None: self.guiim.del_shapes()
		if self.guimx != None: self.guimx.set_data([])
		self.box_display_update() # - the user may still have some manual boxes...
	
	def big_image_change(self):
		global BigImageCache
		image=BigImageCache.get_object(self.get_current_image_name()).get_image(use_alternate=True)
		if self.guiim != None:
			self.guiim.set_data(image)
			self.guiim.del_shapes()
			self.guiim.force_display_update()
		self.box_display_update()
			

	def gen_thumbs(self):
		self.imagethumbs = gen_thumbs(self.image_names,self.get_image_thumb_shrink())
		if self.imagethumbs == -1:
			print "thumbnail process was aborted"
			sys.exit(1)

	def __gen_image_thumbnails_widget(self):
		'''
		Generates image thumbnails for a single image name
		if there is only one image in the image file on disk
		this function returns 0 and does nothing.
		Else thumbnails are generated, and image matrix is initialized (but not shown),
		and 1 is returned
		'''
		# warnilg here
		try: nim = len(self.image_names)
		except: 
			# warning - bad hacking going on
			print "the image name ", image_name, "probably doesn't exist"
			raise Exception
	
		if (nim == 1): return 0
		

		get_application().setOverrideCursor(QtCore.Qt.BusyCursor)
		
		self.gen_thumbs()
		

		self.guimxit=EMImageMXModule(application=get_application())
		self.guimxit.desktop_hint = "rotor"
			
		try:
			self.guimxit.setWindowTitle("Image Thumbs")
		except:
			pass
		#get_application().show_specific(self.guimxit)
		self.guimxit.set_data(self.imagethumbs,soft_delete=True)
		
		try:
			for i in range(0,nim):
				frozen = get_idd_key_entry(self.image_names[i],"frozen_state")
				if frozen != None:
					self.guimxit.set_frozen(frozen,i)
		except: 
			pass # this will happen if fancy widgets are being used
		
		self.guimxit.set_mouse_mode("App")
			#app.setOverrideCursor(QtCore.Qt.BusyCursor)
#		except: 
#			app.setOverrideCursor(QtCore.Qt.ArrowCursor)
#			return 0
#			pass
		
		get_application().setOverrideCursor(QtCore.Qt.ArrowCursor)
		return 1
		
	def get_image_thumb(self,i):
		
		thumb = get_idd_image_entry(self.image_names[i],"image_thumb")
		if thumb == None:
			n = self.get_image_thumb_shrink()
			#global BigImageCache
			image=EMData(self.image_names[i],0)
			#while n > 1:
				#image = image.process("math.meanshrink",{"n":2})
				#n /= 2
			thumb = image.process("math.meanshrink",{"n":n})
			thumb.process_inplace("normalize.edgemean") # if there are lots than they =should all have the same contrast
			set_idd_image_entry(self.image_names[i],"image_thumb",thumb)
			image = None
			print sys.getrefcount(image)
		return thumb
		
	def get_image_thumb_shrink(self):
		if self.itshrink == -1:
			shrink = 1
			inx,iny =  gimme_image_dimensions2D(self.image_names[self.current_image_idx])
			inx /= 2
			iny /= 2
			while ( inx >= 128 and iny >= 128):
				inx /= 2
				iny /= 2
				shrink *= 2
		
			self.itshrink=shrink
		
		return self.itshrink
		
	def box_moved(self,event,scale):
		try:
			dx = (self.moving_box_data[0] - event.x())/scale
			dy = (event.y() - self.moving_box_data[1])/scale
			self.move_box(self.moving_box_data[2],dx,dy)
			self.moving_box_data[0] = event.x()
			self.moving_box_data[1] = event.y()
			self.update_all_image_displays()
		except: pass
		
	def box_released(self,event,lc):
		if lc == None: return # this can happen in strange circumstances
		boxes = self.get_boxes()
		#print boxes
		try: box =  boxes[self.moving_box_data[2]]
		except: return # The user clicked on black
		
		if box.isref :
			self.reference_moved(box)

		im=lc[0]
		self.moving_box_data = [event.x(),event.y(),im]
		if self.guiim != None:
			self.guiim.set_active(im,.9,.9,.4)
			#self.guimx.get_core_object.set_selected(im)
			boxes = self.get_boxes()
			self.guiim.register_scroll_motion(boxes[im].xcorner+boxes[im].xsize/2,boxes[im].ycorner+boxes[im].ysize/2)
		
		self.moving_box_data = None
		
	def box_selected(self,event,lc):
		im=lc[0]
		self.moving_box_data = [event.x(),event.y(),im]
		self.guiim.set_active(im,.9,.9,.4)
		#object = self.guimx.get_core_object()
		#if not object.is_visible(lc[0]) : object.scroll_to(lc[0],True)
		#self.get_mx_gui_image().get_core_object().set_selected([lc[0]],True)
		#		boxes = self.get_boxes()
		#self.guiim.register_scroll_motion(boxes[im].xcorner+boxes[im].xsize/2,boxes[im].ycorner+boxes[im].ysize/2)
		#try:
			##self.guiim.scroll_to(boxes[im].xcorner+boxes[im].xsize/2,boxes[im].ycorner+boxes[im].ysize/2)
			#pass
			
		#except: print "box_selected() scrolling error"

	def get_boxes(self):
		return self.boxable.boxes
	
	def mouse_move(self,event):
		self.mousehandler.mouse_move(event)

	def mouse_wheel(self,event):
		self.mousehandler.mouse_wheel(event)
		
	def mouse_down(self,event) :
		self.mousehandler.mouse_down(event)
		
	def mouse_up(self,event) :
		self.mousehandler.mouse_up(event)
	
	def mouse_drag(self,event) :
		self.mousehandler.mouse_drag(event)
	
	def collision_detect(self,m,boxes):
			
		for box_num,box in enumerate(boxes):
			if m[0]<box.xcorner or m[0]>(box.xcorner +box.xsize) or m[1]<box.ycorner or m[1]>(box.ycorner +box.ysize) :
				# no collision
				continue
			# if we make it here there has been a collision, the box already exists
			return box_num
		
		return -1

	def move_box(self,box_num,dx,dy):
		box = self.get_boxes()[box_num]
		if not self.boxable.is_interactive():
			return
		
		self.boxable.move_box(box,dx,dy,box_num)
			# we have to update the reference also
		self.ptcl[box_num] = box.get_box_image()
			
		x0=box.xcorner+box.xsize/2-1
		y0=box.ycorner+box.ysize/2-1
		self.guiim.add_shape("cen",EMShape([self.shape_string,.9,.9,.4,x0,y0,x0+2,y0+2,1.0]))
		box.shape = EMShape([self.shape_string,box.r,box.g,box.b,box.xcorner,box.ycorner,box.xcorner+box.xsize,box.ycorner+box.ysize,2.0])
		self.guiim.add_shape(box_num,box.shape)
		self.box_display_update_specific(box_num)

	def change_shapes(self,shape_string):
		if shape_string in ["rectpoint","rect","rcircle","rcirclepoint"]:
			self.shape_string = shape_string
			self.box_display_update(True)
		else:
			print "unknown shape string", shapestring
	
	def update_box_colors(self,classify):
		sh=self.guiim.get_shapes()
		for i in classify.items():
			color = BoxingTools.get_color(i[1])
			
			sh[int(i[0])].shape[1] = color[0]
			sh[int(i[0])].shape[2] = color[1]
			sh[int(i[0])].shape[3] = color[2]
			sh[int(i[0])].changed=True
			
		self.box_display_update()

	def keypress(self,event):
		if event.key() == QtCore.Qt.Key_Tab:
			pass
			#self.rendcor = not self.rendcor
			#if ( self.correlation != None ):
				#if self.rendcor == True:
					#self.image.insert_clip(self.correlation,self.correlationcoords)
				#else:
					#self.image.insert_clip(self.correlationsection,self.correlationcoords)
		else: pass
	
	def update_all_image_displays(self):
		self.update_image_display()
		self.update_mx_display()
		if self.guimxit != None: self.guimxit.updateGL()
	
		#context = contextdata.getContext(None)
		#print context
		
	def update_image_display(self):
		if self.guiim != None: self.guiim.updateGL()
		
		#if self.guimxit != None: self.guimxit.get_parent()
		
	def update_mx_display(self):
		#get_application().get_qt_gl_updategl_target(self.guimx).updateGL()
		if self.guimx != None: self.guimx.updateGL()

	def delete_display_boxes(self,numbers):
		'''
		Warning - this won't work unless the numbers go from greatest to smallest- i.e. they are in reverse order
		Deletes shapes displayed by the 2D image viewer
		Pops boxed particles from the list used by the matrix image viewer (for boxes)
		'''
		#print "called delete display shapesS"
		if self.in_display_limbo: return
		
		sh=self.guiim.get_shapes()
		
		for num in numbers:
			#sh=self.guiim.get_shapes()
			k=sh.keys()
			k.sort()
			#del sh[int(num)]
			sh.pop(int(num))
			for j in k:
				if isinstance(j,int):
					if j>num :
						sh[j-1]=sh[j]
						sh.pop(j)
						#del sh[j]
			self.ptcl.pop(num)
						
			
		self.guiim.del_shapes()
		self.guiim.add_shapes(sh)

		#self.guiim.del_shapes()
		#self.guiim.add_shapes(sh)
		#print "now there are",len(sh),"shapes"
		
	def box_image_deleted(self,event,lc,force_image_mx_remove=True):
		box = self.delete_box(lc[0],force_image_mx_remove)
		self.boxable.add_exclusion_particle(box)
		self.guiim.set_other_data(self.boxable.get_exclusion_image(False),self.autoboxer.get_subsample_rate(),True)
		self.guictl.num_boxes_changed(len(self.ptcl))
		self.update_all_image_displays()
		#self.update_ppc()
		
	def delete_box(self,box_num,force_image_mx_remove=True):
		"""
		Deletes the numbered box completely
		Should only be called in the instance where a single box is being deleting - NOT when
		you are deleting a list of boxes sequentially (for that you should use delete_display_boxes
		and something to pop the box from the Boxable. See examples in this code)
		"""
		sh=self.guiim.get_shapes()
		k=sh.keys()
		k.sort()
		del sh[int(box_num)]
		for j in k:
			if isinstance(j,int):
				if j>box_num :
					sh[j-1]=sh[j]
					del sh[j]
					
		if self.guiim != None:
			self.guiim.del_shapes()
			self.guiim.add_shapes(sh)
			self.guiim.set_active(None,.9,.9,.4)
		
		if force_image_mx_remove: 
			#self.ptcl.pop(box_num)
			#self.guimx.set_data(self.ptcl)
			if self.guimx != None:self.guimx.pop_box_image(box_num)

		box = self.boxable.boxes[box_num]
		
		self.boxable.delete_box(box_num)
		# if the boxable was a reference then the autoboxer needs to be told. It will remove
		# it from its own list and potentially do autoboxing
		if box.isref:
			val = self.autoboxer.remove_reference(box)
			if val == 2:
				self.boxable.clear_and_cache(True)
				self.clear_displays()
			self.box_display_update()
				
		return box
	
	
	def box_display_update_specific(self, box_num,force=False):
		ns = {}
		box = self.get_boxes()[box_num]
		
		im=box.get_box_image()
		box.shape = EMShape([self.shape_string,box.r,box.g,box.b,box.xcorner,box.ycorner,box.xcorner+box.xsize,box.ycorner+box.ysize,2.0])
		
		#if not force and not box.isref and not box.ismanual:
			#box.shape.isanimated = True
			#box.shape.blend = 0
		
		ns[box_num]=box.shape
		if box_num>=len(self.ptcl) : self.ptcl.append(im)
		else : self.ptcl[box_num]=im
			
		self.guiim.add_shapes(ns)
		self.set_ptcl_mx_data(self.ptcl)
		
		self.guictl.num_boxes_changed(len(self.ptcl))
		
		if isinstance(self.guimxit,EMImageRotorModule):
			self.guimxit.set_extra_hud_data([len(self.ptcl)])
	
		if isinstance(self.guimxit,EMImageRotorModule):
			othershapes = {}
			for shape in ns:
				s = ns[shape].getShape()
				othershapes[shape] = EMShape(["point",s[1],s[2],s[3],(s[4]+s[6])/2,(s[5]+s[7])/2,2])
	
			self.guimxit.set_shapes(othershapes,self.get_image_thumb_shrink(),self.current_image_idx)

		self.update_all_image_displays()
		
		
	def box_display_update(self,force=False,stop_animation=False):
		
		ns = {}
		idx = 0
		#self.ptcl = []
		# get the boxes
		boxes =self.get_boxes()
		register_animation = False
		for j,box in enumerate(boxes):
			if not box.changed and not force:
				idx += 1
				continue
			
			box.changed=False
		
			im=box.get_box_image()
			box.shape = EMShape([self.shape_string,box.r,box.g,box.b,box.xcorner,box.ycorner,box.xcorner+box.xsize,box.ycorner+box.ysize,2.0])
			
			if not force and not box.isref and not box.ismanual:
				if not stop_animation:
					box.shape.isanimated = True
					box.shape.set_blend(0)
					register_animation = True
			ns[idx]=box.shape
			if idx>=len(self.ptcl) : self.ptcl.append(im)
			else : self.ptcl[idx]=im
			idx += 1
		
		if stop_animation: register_animation=False
		self.guiim.add_shapes(ns,register_animation)
		self.set_ptcl_mx_data(self.ptcl)
		
		self.guictl.num_boxes_changed(len(self.ptcl))
		if isinstance(self.guimxit,EMImageRotorModule):
			self.guimxit.set_extra_hud_data([len(self.ptcl)])
			
		if isinstance(self.guimxit,EMImageRotorModule):
			othershapes = {}
			for shape in ns:
				s = ns[shape].getShape()
				othershapes[shape] = EMShape(["point",s[1],s[2],s[3],(s[4]+s[6])/2,(s[5]+s[7])/2,2])
	
			self.guimxit.set_shapes(othershapes,self.get_image_thumb_shrink(),self.current_image_idx)


		#if self.guimxitp !=None:
			#othershapes = {}
			#for shape in ns:
				#s = ns[shape].getShape()
				#othershapes[shape] = EMShape(["point",s[1],s[2],s[3],(s[4]+s[6])/2,(s[5]+s[7])/2,2])
	
			#if isinstance(self.guimxit,EMImageRotorModule): self.guimxit.set_shapes(othershapes,self.get_image_thumb_shrink(),self.current_image_idx)
		self.update_all_image_displays()

	def force_autobox(self,bool):
		'''
		like hitting refresh - everything is updated
		'''
		self.boxable.clear_and_cache(True)
		self.autoboxer.write_image_specific_references_to_db(self.boxable.get_image_name())
		self.boxable.get_references_from_db()
		self.autoboxer.auto_box(self.boxable, True,True)
		self.box_display_update()

	def toggle_dynapix(self,bool):
		self.dynapix = bool
		self.autoboxer.set_interactive_mode(self.dynapix)
		
	def done(self):
		if self.boxable != None: self.boxable.cache_exc_to_db()
		#print "here we
		if self.guictl_module != None: self.guictl_module.closeEvent(None)
		if self.guiim != None: self.guiim.closeEvent(None)
		if self.guimx != None: self.guimx.closeEvent(None)
		if self.guimxit != None: self.guimxit.closeEvent(None)
		if self.dab != None : self.dab.close()
		if self.form != None : self.form.closeEvent(None)
		
		BigImageCache.clear_cache()
		ExclusionImageCache.clear_cache()
		SincBlackmanSubsampleCache.clear_cache()
		CoarsenedFlattenedImageCache.clear_cache()
		SubsamplerCache.clear_cache()
		InverseSigmaImageCache.clear_cache()
		BinaryCircleImageCache.clear_cache()
		
		
		self.emit(QtCore.SIGNAL("module_closed"))

	def closeEvent(self,event):
		self.done()

	def try_data(self,data,thr):
		print data, thr
		print 'try that was pressed, this feature is currently disabled'
		
	def opt_params_updated(self,thresh,profile,radius):
		self.guictl.update_data(thresh,profile,radius)
		
	def set_selection_mode(self,selmode):
		if self.autoboxer.set_selection_mode(selmode):
			self.box_display_update()
		
	def set_profile_comparitor(self,cmp_mode):
		if self.autoboxer.set_cmp_mode(cmp_mode):
			self.box_display_update()
			
	def classify(self,bool):
		self.boxable.classify()
		
	def erase_toggled(self,bool):
		# for the time being there are only two mouse modes
		
		if bool == True:
			self.mouse_handlers["erasing"].set_mode(Boxable.ERASE)
			self.mousehandler = self.mouse_handlers["erasing"]
		else:
			self.guiim.add_eraser_shape("None",None)
			self.update_image_display()
			self.mousehandler = self.mouse_handlers["boxing"]
			
	
	def unerase_toggled(self,bool):
		if bool == True:
			self.mouse_handlers["erasing"].set_mode(Boxable.UNERASE)
			self.mousehandler = self.mouse_handlers["erasing"]
		else:
			self.guiim.add_eraser_shape("None",None)
			self.update_image_display()
			self.mousehandler = self.mouse_handlers["boxing"]
	
	def update_erase_rad(self,rad):
		self.mouse_handlers["erasing"].set_erase_radius(rad)
	
	def set_dummy_box(self,box):
		self.autoboxer.set_dummy_box(box)
		self.box_display_update()
		
	def remove_dummy(self):
		self.autoboxer.set_dummy_box(None)
		self.box_display_update()
	
	def run_output_dialog(self):
		from emsprworkflow import E2BoxerProgramOutputTask
		exclusions = []
		for name in self.image_names:
			if get_idd_key_entry(name,"quality") == 0:
				exclusions.append(name)
		if self.output_task != None: return
		self.output_task = E2BoxerProgramOutputTask(get_application(),self.image_names,self,exclusions)
		QtCore.QObject.connect(self.output_task.emitter(),QtCore.SIGNAL("task_idle"),self.on_output_task_idle)
		self.output_task.run_form()

	def on_output_task_idle(self):
		self.output_task = None
	
	def write_all_box_image_files(self,box_size,forceoverwrite=False,imageformat="hdf",normalize=True,norm_method="normalize.edgemean"):
		self.write_box_image_files(self.image_names,box_size,forceoverwrite,imageformat,normalize,norm_method)
		
	def write_box_image_files(self,image_names,box_size,forceoverwrite=False,imageformat="hdf",normalize=True,norm_method="normalize.edgemean",invert=False):
		self.boxable.cache_exc_to_db()
		progress = EMProgressDialogModule(get_application(),"Writing boxed images", "Abort", 0, len(image_names),None)
		progress.qt_widget.show()
		get_application().setOverrideCursor(Qt.BusyCursor)
		for i,image_name in enumerate(image_names):
			try:
				project_db = EMProjectDB()
				data = project_db[get_idd_key(image_name)]
				trim_autoboxer = project_db[data["autoboxer_unique_id"]]["autoboxer"]
				autoboxer = SwarmAutoBoxer(self)
				autoboxer.become(trim_autoboxer)
				#print "writing box images for",image_name,"using",data["autoboxer_unique_id"]
			except:
				autoboxer = self.autoboxer
				#print "writing box images or",image_name,"using currently stored autoboxer"
			
			
			if isinstance(autoboxer,SwarmAutoBoxer):
				boxable = Boxable(image_name,self,autoboxer)
				#boxable.clear_and_cache() ATTENTION
				if boxable.is_excluded():
					print "Image",image_name,"is excluded and being ignored"
					continue
				
				mode = self.autoboxer.get_mode()
				self.autoboxer.set_mode_explicit(SwarmAutoBoxer.COMMANDLINE)
				self.autoboxer.auto_box(boxable,False)
				self.autoboxer.set_mode_explicit(mode)
				
				boxable.write_box_images(box_size,forceoverwrite,imageformat,normalize,norm_method,invert)
		
			else: 
				self.autoboxer.write_box_images(self.boxable, normalize, norm_method)
				
			if progress.qt_widget.wasCanceled():
				# yes we could probably clean up all of the images that were written to disk but not time...
				break
				
			progress.qt_widget.setValue(i+1)
			get_application().processEvents()
		get_application().setOverrideCursor(Qt.ArrowCursor)
		progress.qt_widget.close()
 
	def write_all_coord_files(self,box_size,forceoverwrite=False):
		self.write_coord_files(self.image_names,box_size,forceoverwrite)

	def write_coord_files(self,image_names,box_size,forceoverwrite=False):
		self.boxable.cache_exc_to_db()
		progress = EMProgressDialogModule(get_application(),"Writing Coordinate  Files", "Abort", 0, len(image_names),None)
		progress.qt_widget.show()
		get_application().setOverrideCursor(Qt.BusyCursor)
		for i,image_name in enumerate(image_names):
			
			try:
				project_db = EMProjectDB()
				data = project_db[get_idd_key(image_name)]
				trim_autoboxer = project_db[data["autoboxer_unique_id"]]["autoboxer"]
				autoboxer = SwarmAutoBoxer(self)
				autoboxer.become(trim_autoboxer)
				#print "writing box coordinates for",image_name,"using",data["autoboxer_unique_id"]
			except:
				autoboxer = self.autoboxer
				#print "writing box coordinates for",image_name,"using currently stored autoboxer"
			if isinstance(autoboxer,SwarmAutoBoxer):
				boxable = Boxable(image_name,self,autoboxer)
			
				if boxable.is_excluded():
					print "Image  ",image_name,"   is excluded and being ignored"
					continue
			
				mode = autoboxer.get_mode()
				autoboxer.set_mode_explicit(SwarmAutoBoxer.COMMANDLINE)
				autoboxer.auto_box(boxable,False)
				autoboxer.set_mode_explicit(mode)
				boxable.write_coord_file(box_size,forceoverwrite)
			else:
				self.autoboxer.write_box_coords( boxable )
				
			if progress.qt_widget.wasCanceled():
				# yes we could probably clean up all of the images that were written to disk but not time...
				break;
			
				
			progress.qt_widget.setValue(i+1)
			get_application().processEvents()
		get_application().setOverrideCursor(Qt.ArrowCursor)
		progress.qt_widget.close()
 
	
	def center(self,technique):
		
		if self.boxable.center(technique):
			self.box_display_update()
		else:
			print 'technique',technique,'is unsupported - check back tomorrow'
			
	def toggle_frozen(self):
		if self.boxable.is_excluded() : return
		self.boxable.toggle_frozen()
		self.boxable.write_to_db()
	
		self.guiim.set_frozen(self.boxable.is_frozen())
		
		if isinstance(self.guimxit,EMImageRotorModule): self.guimxit.set_frozen(self.boxable.is_frozen(),self.current_image_idx)
		if not self.boxable.is_frozen():
			self.change_current_autoboxer(self.boxable.get_autoboxer_id(),False)
		else:	
			self.update_image_display()
		
	def change_image_quality(self,val):
		self.boxable.set_quality(val)
		if val == Boxable.EXCLUDE:
			self.boxable.set_frozen(False) # If it was already frozen then setting to excluded overrides this
			self.guiim.set_excluded(True)
			self.guiim.set_frozen(False)
			if isinstance(self.guimxit,EMImageRotorModule): self.guimxit.set_frozen(False,self.current_image_idx)
			self.boxable.clear_and_cache(True) # tell boxable to clear its autoboxes and references -
			self.clear_displays() # tell the display to clear itself
		elif self.guiim.set_excluded(False):
			self.update_image_display()
			
		self.boxable.write_to_db() # make sure the infromation changes that just occured are written to the DB
		

	def set_boxing_method(self,ref,manual):
		'''
		Okay could do it with one argument but leaving it this way makes it obvious
		'''
		
		if ref and manual:
			print 'error, you cant set both ref and manual'
			
		if ref:
			self.mouse_handlers["boxing"].set_mode(EMBoxerModule.REFERENCE_ADDING)
		elif manual:
			self.mouse_handlers["boxing"].set_mode(EMBoxerModule.MANUALLY_ADDING)
	
	def autoboxer_db_changed(self):
		self.guictl.update_ab_table()

	def add_new_autoboxer_db(self, n):
		if not self.boxable.is_interactive():
			return None
		autoboxer = SwarmAutoBoxer(self)
		autoboxer.set_box_size_explicit(self.box_size)
		autoboxer.set_interactive_mode(self.dynapix)
		autoboxer_db_string = "autoboxer_"+autoboxer.get_creation_ts()
		trim_autoboxer = TrimSwarmAutoBoxer(autoboxer)
		convenience_name = "New " + str(n)
		trim_autoboxer.set_convenience_name(convenience_name)
		trim_autoboxer.write_to_db()
		return convenience_name
	
	def add_copy_autoboxer_db(self,autoboxer_id,n):
		#print "adding a copy of the ab with id is",autoboxer_id
		project_db = EMProjectDB()
		trim_autoboxer = copy(project_db[autoboxer_id]["autoboxer"])
		trim_autoboxer.set_creation_ts(gm_time_string())
		convenience_name = "Copy " + str(n)
		trim_autoboxer.set_convenience_name(convenience_name)
		trim_autoboxer.write_to_db()
		return convenience_name
		
		#print "done"

class AutoBoxerSelectionsMediator:
	'''
	A class for coordinating the EMBoxerModulePanel and the the EMBoxerModule in relation
	to adding and removing AutoBoxers, and changing which Boxables use which
	AutoBoxer etc
	'''
	def get_desktop_hint(self): return "inspector"
	def __init__(self,parent):
		if not isinstance(parent,EMBoxerModule):
			print "error, the AutoBoxerSelectionsMediator must be initialized with a EMBoxerModule type as its first constructor argument"
			return
		self.parent=weakref.ref(parent)
		self.current_image_idx_name = None
		self.image_names = []
		self.name_map = {}
		self.dict_data = {}
		self.set_current_image_name(parent.get_current_image_name())
		self.set_image_names(parent.get_image_names())
		
	def set_current_image_name(self,image_name):
		self.current_image_idx_name = image_name
		
	def set_image_names(self,image_names):
		self.image_names = image_names
		
	def get_autoboxer_data(self):
		project_db = EMProjectDB()
		self.dict_data = {}
		self.name_map = {}
		for i in project_db.keys():
			try:
				if i[0:10] == "autoboxer_":
					tag = project_db[i]["convenience_name"]
					self.dict_data[tag] = []
					self.name_map[i] = tag
			except:
				print "error couldn't handle",i
		
		for image_name in self.image_names:
			found = False
			try:
				data = project_db[get_idd_key(image_name)]
				#trim_autoboxer = project_db[data["autoboxer_unique_id"]]
				if data != None: found = True
			except: pass
			
			if found:
				try:
					self.dict_data[self.name_map[data["autoboxer_unique_id"]]].append(strip_after_dot(image_name))
				except: pass
					#print data
					#try:
						#print "error, an autoboxer has been lost, its stamp is",data["autoboxer_unique_id"]
					#except: pass
		return self.dict_data
	
	def get_current_autoboxer_ts(self):
		try:
			return self.name_map["autoboxer_"+self.parent().get_current_autoboxer_ts()]
		except:
			return None

	def add_new_autoboxer(self):
		return self.parent().add_new_autoboxer_db(self.get_total_autoboxers())

	def add_copy_autoboxer(self,tag):
		autoboxer_id = self.__get_autoboxer_id_from_tag(tag)
		if autoboxer_id != None:
			return self.parent().add_copy_autoboxer_db(autoboxer_id,self.get_total_autoboxers())
		else:
			print "error, couldn't find autoboxer from tag",tag
			return None
	
	def change_current_autoboxer(self,tag):
		# FIXME - is there any way to use a bidirectional map?pyt
		autoboxer_id = self.__get_autoboxer_id_from_tag(tag)
		if autoboxer_id != None:
			return self.parent().change_current_autoboxer(autoboxer_id)
		else:
			print "error, couldn't get autoboxer_id"
			return None
				
	def update_db_convenience_name(self,newname,oldname):
		'''
		Updates the unique name of the autoboxer in the DB
		if the name is already used then False is returned
		and the calling function should act on this to stop the
		name change, for example within a widget
		'''
		project_db = EMProjectDB()
		autoboxer_id = self.__get_autoboxer_id_from_tag(oldname)
		if autoboxer_id != None:
			self.name_map.pop(autoboxer_id)
			if self.__name_not_already_present(newname):
				self.name_map[autoboxer_id] = newname
				autoboxer = project_db[autoboxer_id]["autoboxer"]
				autoboxer.set_convenience_name(newname)
				autoboxer.write_to_db()
				return True
			else:
				self.name_map[autoboxer_id]["convenience_name"] = oldname
				return False
					
	def get_total_autoboxers(self):
		return len(self.name_map)
	
	def associated_images(self,tag):
		return self.dict_data[tag]
	
	def remove(self,tag):
		autoboxer_id = self.__get_autoboxer_id_from_tag(tag)
		db = db_open_dict("bdb:e2boxer.cache")
		del(db[autoboxer_id])
#		project_db = EMProjectDB()
#		project_db.pop(autoboxer_id)
		self.dict_data.pop(tag)
	
	def __name_not_already_present(self,name):
		for names in self.name_map.items():
			if names[1] == name:
			 	return False
		
		return True
	
	def __get_autoboxer_id_from_tag(self,tag):
		for names in self.name_map.items():
			if names[1] == tag:
				return names[0]
			
		print "error, couldn't find",tag,"in the namemap"
		return None
	
	def toggle_frozen(self,tag,bool):
		#frozen = boxable.is_frozen()
		#if frozen:
			#new_name = self.add_copy_autoboxer(tag)
			#self.get_autoboxer_data() # FIXME this is inefficient, could just add to self.dict_data etc
			#autoboxer_id = self.__get_autoboxer_id_from_tag(new_name)
			
			#boxable.set_autoboxer_id(autoboxer_id)
			#boxable.write_to_db()
		#boxable.write_to_db()
		self.parent().toggle_frozen()
		
	def clear_current(self):
		new_name = self.add_new_autoboxer()
		
		# if the new_name is none then the current Boxable is frozen!
		if new_name == None: return new_name
		
		self.get_autoboxer_data() # FIXME this is inefficient, could just add to self.dict_data etc
		autoboxer_id = self.__get_autoboxer_id_from_tag(new_name)
		boxable = self.parent().get_boxable()
		boxable.set_autoboxer_id(autoboxer_id)
		boxable.write_to_db()
		boxable.clear_and_cache(True)
		self.parent().clear_displays()
		return new_name

	def closeEvent(self,event=None):
		self.boxable.cache_exc_to_db()
		project_db = EMProjectDB()
		project_db.close()
		

		E2saveappwin("boxer","imagegeom",self.guiim)
		E2saveappwin("boxer","matrixgeom",self.guimx)
		E2saveappwin("boxer","controlgeom",self.guictl)
		#E2setappval("boxer","matrixnperrow",self.guimx.nperrow)
		try:
			E2setappval("boxer","imcontrol",self.guiim.inspector.isVisible())
			if self.guiim.inspector.isVisible() : E2saveappwin("boxer","imcontrolgeom",self.guiim.inspector)
		except : E2setappval("boxer","imcontrol",False)
		try:
			E2setappval("boxer","mxcontrol",self.guimx.inspector.isVisible())
			if self.guimx.inspector.isVisible() : E2saveappwin("boxer","mxcontrolgeom",self.guimx.inspector)
		except : E2setappval("boxer","mxcontrol",False)
		

		
def histogram1d( data, nbin, presize=0 ) :
	fmax = max( data )
	fmin = min( data )
	binsize = (fmax - fmin)/(nbin-2*presize)
	start = fmin - binsize*presize
	region = [None]*nbin
	hist = [None]*nbin
	for i in xrange(nbin):
		region[i] = start + (i+0.5)*binsize
		hist[i] = 0

	for d in data:
		id = int( (d-start)/binsize )
		hist[id]+=1

	return region,hist

class CcfHistogram(QtGui.QWidget):

	def __init__(self, parent):	
		QtGui.QWidget.__init__(self,parent)
		self.parent= weakref.ref(parent)  # this needs to be a weakref ask David Woolford for details, but otherwise just call self.parent() in place of self.parent
		self.data=None
                self.setMinimumSize(QtCore.QSize(256,256))
		self.PRESIZE = 28

		self.cached_boxes = []
		
	def clear( self ):
		self.ccfs = None
		self.data = None
		self.parent().ccf_range.setText( "Range: (N/A, N/A)" )
		self.parent().threshold_low.setText( "N/A" )
		self.parent().threshold_hgh.setText( "N/A" )

	def set_data( self, data ):
		self.ccfs = data
		self.nbin = self.width()
                self.data = histogram1d( data, self.nbin, self.PRESIZE )

		hmin = self.data[0][0]
		hmax = self.data[0][-1]

		info = "Range: (%8.4f, %8.4f)" % (hmin, hmax)

		self.parent().ccf_range.setText( info )
		self.parent().threshold_low.setText( str(hmin) )
		self.parent().threshold_hgh.setText( str(hmax) )
		
		self.tickers =[1, self.width()-2]

	def paintEvent( self, event ):
		p=QtGui.QPainter()
		p.begin(self)
		p.setBackground(QtGui.QColor(16,16,16))
		p.eraseRect(0,0,self.width(),self.height())
		p.setPen(QtCore.Qt.darkGray)

		if self.data is None:
			return

		hmax = max( self.data[1] )
                for i in xrange( len(self.data[1]) ):
			h = self.data[1][i]
                        p.drawLine(i, self.height(), i, int(self.height()*(1-0.9*h/hmax)) )
		
		self.drawTicker( self.tickers[0], p )
		self.drawTicker( self.tickers[1], p )
		p.end()

	def mousePressEvent(self, event):
		if event.button()==Qt.LeftButton:
			x = event.x()
			dis1 = abs( x - self.tickers[0] )
			dis2 = abs( x - self.tickers[1] )

			if dis1 < dis2:
				self.cur_ticker = 0
			else:
				self.cur_ticker = 1

			if not hasattr( self, "shapes" ):
				self.shapes = self.parent().target.guiim.get_shapes().copy()

	def mouseMoveEvent(self, event):
		if event.buttons()&Qt.LeftButton and event.x() > 0 :
			self.tickers[self.cur_ticker] = event.x()
			self.repaint()

			x = event.x()
			if x < 0 :
				thr = self.data[0][0]
			elif x > len(self.data[0]) : 
				thr = self.data[0][-1]
			else :
				thr = self.data[0][x-1]

			if self.cur_ticker==0:
				self.parent().threshold_low.setText( str(thr) )
			else:
				self.parent().threshold_hgh.setText( str(thr) )


			thr_low = float( self.parent().threshold_low.text() )
			thr_hgh = float( self.parent().threshold_hgh.text() )
	
			#guiim = self.parent.target.guiim
			target = self.parent().get_target() # using getter functions is much preferable. Add support for them yourself if you need them
			boxable = target.get_boxable()
			
			[added_boxes,added_ref_boxes] = boxable.update_included_boxes_hist(thr_low,thr_hgh)
			[lost_boxes,refs] = boxable.update_excluded_boxes_hist(thr_low,thr_hgh)
			
			if len(lost_boxes) != 0:
				target.delete_display_boxes(lost_boxes)
			
				
			if len(added_boxes) != 0 or len(lost_boxes) != 0 :
				target.box_display_update()
				target.update_all_image_displays()

	def drawTicker( self, newpos, p ) :
		p.setPen(Qt.yellow)
		for i in xrange( newpos, newpos+2):
			p.drawLine( i, self.height(), i, int(0.2*self.height()) )

class EMBoxerModulePanelModule(EMQtWidgetModule):
	'''
	params should be a list of ParamDef objects
	application should be an EMAN2 type application
	After initializing this object you should probably call application.show(this)
	or application.show_specific(this), depending on what you're doing
	'''
	def __init__(self,application,target,ab_sel_mediator):
		self.application = weakref.ref(application)
		self.widget = EMBoxerModulePanel(self,target,ab_sel_mediator)
		EMQtWidgetModule.__init__(self,self.widget)
		
	def get_desktop_hint(self):
		return "inspector"

class EMBoxerModulePanel(QtGui.QWidget):
	def get_desktop_hint(self):
		return "inspector"

	def __init__(self,module,target,ab_sel_mediator) :
		
		QtGui.QWidget.__init__(self,None)
		self.target=weakref.ref(target)
		self.ab_sel_mediator = ab_sel_mediator
		self.module = module # please set this to be a EMBoxerModulePanelModule
		
		self.vbl = QtGui.QVBoxLayout(self)
		self.vbl.setMargin(0)
		self.vbl.setSpacing(6)
		self.vbl.setObjectName("vbl")
		
		self.tabwidget = QtGui.QTabWidget(self)
		self.insert_main_tab()
		self.insert_advanced_tab()
		self.insert_view_tab()
		self.vbl.addWidget(self.tabwidget)
		
		self.dummybox = Box()
		self.dummybox.isdummy = True
		self.currentlyselected = -1 # used in the ab_table
		
		self.lock = False
		self.connect(self.bs,QtCore.SIGNAL("editingFinished()"),self.new_box_size)
	
		self.connect(self.thr,QtCore.SIGNAL("textChanged"),self.new_threshold)
		self.threshold_press = -1
		self.connect(self.thr,QtCore.SIGNAL("sliderReleased"),self.new_threshold_release)
		self.connect(self.thr,QtCore.SIGNAL("sliderPressed"),self.new_threshold_press)
		
#		self.connect(self.trythat,QtCore.SIGNAL("clicked(bool)"),self.try_dummy_parameters)
#		self.connect(self.reset,QtCore.SIGNAL("clicked(bool)"),self.target.remove_dummy)
		self.connect(self.thrbut, QtCore.SIGNAL("clicked(bool)"), self.selection_mode_changed)
	
	def new_threshold_press(self,value):
		self.threshold_press = value
    	
	def new_threshold_release(self,value):
		if ( value != self.threshold_press ):
			self.threshold_press = value
			self.new_threshold(value)

	def keyPressEvent(self,event):
		if event.key() == Qt.Key_F1:
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
	#def centerpushed(self,unused):
		#self.target.center(str(self.centerooptions.currentText()))
	def get_target(self): return self.target()

	def insert_main_tab(self):
		# this is the box layout that will store everything
		self.main_inspector = QtGui.QWidget()
		self.main_vbl =  QtGui.QVBoxLayout(self.main_inspector)

		self.boxingvbl = QtGui.QVBoxLayout()
		
		self.boxinghbl1=QtGui.QHBoxLayout()
		self.boxinghbl1.setMargin(0)
		self.boxinghbl1.setSpacing(2)
		
		self.refbutton=QtGui.QPushButton( QtGui.QIcon(get_image_directory() + "black_box.png"), "Reference")
		self.refbutton.setCheckable(1)
		self.refbutton.setChecked(True)
		self.boxinghbl1.addWidget(self.refbutton)
		
		self.manualbutton=QtGui.QPushButton(QtGui.QIcon(get_image_directory() + "white_box.png"), "Manual")
		self.manualbutton.setCheckable(1)
		self.manualbutton.setChecked(False)
		self.boxinghbl1.addWidget(self.manualbutton)
		
		self.lblbs=QtGui.QLabel("Box Size:",self)
		self.boxinghbl1.addWidget(self.lblbs)
		self.pos_int_validator = QtGui.QIntValidator(self)
		self.pos_int_validator.setBottom(1)
		self.bs = QtGui.QLineEdit(str(self.target().box_size),self)
		self.bs.setValidator(self.pos_int_validator)
		self.boxinghbl1.addWidget(self.bs)
		
		self.boxingvbl.addLayout(self.boxinghbl1)
	
		self.boxinghbl3=QtGui.QHBoxLayout()
		self.dynapix = QtGui.QCheckBox("Dynapix")
		self.dynapix.setChecked(self.target().dynapix)
		self.boxinghbl3.addWidget(self.dynapix)

		self.method=QtGui.QComboBox()
		self.swarm_icon = QtGui.QIcon(get_image_directory() + "swarm_icon.png")
		self.method.addItem( self.swarm_icon, "Swarm" )
		self.setWindowIcon( QtGui.QIcon(get_image_directory() + "green_boxes.png") )
		self.pp_icon = QtGui.QIcon(get_image_directory() + "pp_boxer_icon.png");
		self.method.addItem( self.pp_icon,"Gauss Conv" )
		self.boxinghbl3.addWidget( self.method )

		self.autobox=QtGui.QPushButton(QtGui.QIcon(get_image_directory() + "green_boxes.png"), "Autobox")
		self.boxinghbl3.addWidget(self.autobox)
		self.boxingvbl.addLayout(self.boxinghbl3)

		self.boxinghbl2=QtGui.QHBoxLayout()
		self.boxinghbl2.setMargin(2)
		self.boxinghbl2.setSpacing(6)
		#self.vbl.addLayout(self.hbl1)
		
		self.erasepic = QtGui.QIcon(get_image_directory() + "boxer_erase.png");
		self.erase=QtGui.QPushButton(self.erasepic,Boxable.ERASE)
		self.erase.setCheckable(1)
		self.boxinghbl2.addWidget(self.erase)
		
		self.unerasepic = QtGui.QIcon(get_image_directory() + "boxer_unerase.png");
		self.unerase=QtGui.QPushButton(self.unerasepic,Boxable.UNERASE)
		self.unerase.setCheckable(1)
		self.boxinghbl2.addWidget(self.unerase)
		
		self.eraseradtext=QtGui.QLabel("Erase Radius:",self)
		self.boxinghbl2.addWidget(self.eraseradtext)
		
		self.eraserad = QtGui.QLineEdit(str(self.target().eraseradius),self)
		self.boxinghbl2.addWidget(self.eraserad)
		self.eraserad.setEnabled(False)

		self.boxingvbl.addLayout(self.boxinghbl2)

		self.invert_contrast_mic = QtGui.QCheckBox("Invert Contrast")
		self.invert_contrast_mic.setChecked(False)
		self.boxingvbl.addWidget(self.invert_contrast_mic,0, Qt.AlignLeft)
		self.invert_contrast_mic.setEnabled(False) # this presents a problem with the Swarm approach
		self.connect(self.invert_contrast_mic,QtCore.SIGNAL("clicked(bool)"),self.invert_contrast_mic_toggled)

		self.boxinghbl4=QtGui.QHBoxLayout()
		self.togfreeze=QtGui.QPushButton(QtGui.QIcon(get_image_directory() + "freeze_swirl.png"),"Toggle Freeze")
		self.boxinghbl4.addWidget(self.togfreeze)
		self.clear=QtGui.QPushButton("Clear")
		self.boxinghbl4.addWidget(self.clear)
		
		self.imagequality=QtGui.QLabel("Image Quality:",self)
		self.boxinghbl4.addWidget(self.imagequality)
		
		self.imagequalities = QtGui.QComboBox(self)
		for metadata in Boxable.QUALITY_META_DATA:
			self.imagequalities.addItem(metadata)
		#self.imagequalities.setSelection(Boxable.AVERAGE)
		self.imagequalities.setCurrentIndex(Boxable.QUALITY_META_DATA_MAP[Boxable.AVERAGE])
		self.boxinghbl4.addWidget(self.imagequalities)

		self.boxingvbl.addLayout(self.boxinghbl4)

		self.infohbl = QtGui.QHBoxLayout()
		self.info = QtGui.QLabel("%d Boxes"%len(self.target().get_boxes()),self)
		#self.ppc = QtGui.QLabel("%f particles per click"%0,self)
		self.infohbl.addWidget(self.info)
		self.boxingvbl.addLayout(self.infohbl)
		
		self.interactiveboxing = QtGui.QGroupBox("Interactive Boxing")
		self.interactiveboxing.setLayout(self.boxingvbl)
		self.main_vbl.addWidget(self.interactiveboxing)

		self.gen_output_but=QtGui.QPushButton("Write output")
		self.main_vbl.addWidget(self.gen_output_but)

		self.classifybut=QtGui.QPushButton("Classify")
		self.main_vbl.addWidget(self.classifybut)
		
		
		self.done=QtGui.QPushButton("Done")
		self.main_vbl.addWidget(self.done)
		
		self.tabwidget.addTab(self.main_inspector,"Main")
		
		self.connect(self.eraserad,QtCore.SIGNAL("editingFinished()"),self.update_erase_rad)
		self.connect(self.erase, QtCore.SIGNAL("clicked(bool)"), self.erase_toggled)
		self.connect(self.unerase, QtCore.SIGNAL("clicked(bool)"), self.unerase_toggled)
		self.connect(self.autobox,QtCore.SIGNAL("clicked(bool)"),self.target().force_autobox)
		self.connect(self.togfreeze,QtCore.SIGNAL("clicked(bool)"),self.toggle_frozen)
		self.connect(self.clear,QtCore.SIGNAL("clicked(bool)"),self.clear_current)
		self.connect(self.dynapix,QtCore.SIGNAL("clicked(bool)"),self.dynapix_toggled)
		
		self.connect(self.refbutton, QtCore.SIGNAL("clicked(bool)"), self.ref_button_toggled)
		self.connect(self.manualbutton, QtCore.SIGNAL("clicked(bool)"), self.manual_button_toggled)
		
		self.connect(self.method, QtCore.SIGNAL("activated(int)"), self.method_changed)

		QtCore.QObject.connect(self.imagequalities, QtCore.SIGNAL("currentIndexChanged(QString)"), self.image_quality_changed)
		
		self.connect(self.done,QtCore.SIGNAL("clicked(bool)"),self.target().done)
		self.connect(self.classifybut,QtCore.SIGNAL("clicked(bool)"),self.target().classify)
		self.connect(self.gen_output_but,QtCore.SIGNAL("clicked(bool)"),self.target().run_output_dialog)

	def invert_contrast_mic_toggled(self):
		from EMAN2 import Util
		image_name = self.target().boxable.get_image_name()
		img = BigImageCache.get_image_directly( image_name )

		[avg,sigma,fmin,fmax] = Util.infomask( img, None, True )

		print "before invert, info: ", Util.infomask( img, None, True )
		img -= avg
		img *= -1
		img += avg
		BigImageCache.get_object(image_name).register_alternate(img)
		self.target().big_image_change()
		print " after invert, info: ", Util.infomask( img, None, True )

	def set_method( self, name ):
		if name[0:5] == "Swarm":
			tabid = self.tabwidget.indexOf( self.pawel_option )
			if tabid != -1:
				self.tabwidget.removeTab( tabid )
				self.tabwidget.insertTab( tabid, self.david_option, "Swarm Advanced" )
				self.autobox.setText("Autobox")
				self.setWindowIcon( QtGui.QIcon(get_image_directory() + "green_boxes.png") )
				self.invert_contrast_mic.setEnabled(False)
		else:
			assert name[0:5]=="Gauss"
			tabid = self.tabwidget.indexOf( self.david_option )
			if tabid != -1:
				self.tabwidget.removeTab( tabid )
				self.tabwidget.insertTab( tabid, self.pawel_option, "Gauss Advanced" )
				self.autobox.setText("Run")
				self.setWindowIcon( self.pp_icon )
				self.invert_contrast_mic.setEnabled(True)

	def method_changed(self, methodid):

		name = self.method.itemText( methodid )

		self.set_method( name )

		self.target().set_autoboxer(self.target().image_names[0], name[0:5])

	def set_image_quality(self,integer):
		self.lock = True
		self.imagequalities.setCurrentIndex(int(integer))
		self.lock = False
		
	def image_quality_changed(self,val):
		if self.lock == False:
			self.target().change_image_quality(str(val))

	def clear_current(self,unused):
		self.lock = True
		new_name = self.ab_sel_mediator.clear_current()
		self.lock = False
		if new_name != None: # if the boxable wasn't frozen...
			self.update_ab_table()
			self.setChecked(new_name)
		
	def toggle_frozen(self,bool):
		self.lock = True
		self.ab_sel_mediator.toggle_frozen(self.col1[self.currentlyselected].text(),bool)
		self.lock = False

	def setChecked(self,tag):
		
		for i,col in enumerate(self.col1):
			if str(col.text()) == tag:
				col.setCheckState(Qt.Checked)
				self.currentlychecked = i
				break
	
	def ab_table_cell_changed(self,i,j):
		if i >= len(self.col1): return
		if self.lock:
			return
		#data = self.ab_sel_mediator.get_autoboxer_data()
		#data = data.items()
		if str(self.col1[i].text()) != self.colnames[i]:
			if not self.ab_sel_mediator.update_db_convenience_name(str(self.col1[i].text()),self.colnames[i]):
				self.col1[i].setText(self.colnames[i])
			self.lock = False
			return
		
		self.lock = False
		
		if j == 0: # we're in the first row, something happened, maybe a check change
			self.lock = True
			#try:
			if i == self.currentlychecked:
				# just make sure the check stays on
				self.col1[self.currentlychecked].setCheckState(Qt.Checked)
			elif self.col1[i].checkState() == Qt.Checked:
				# uncheck the previously selected one
				self.col1[self.currentlychecked].setCheckState(Qt.Unchecked)
				self.col1[i].setCheckState(Qt.Checked)
				self.currentlychecked = i
				self.ab_sel_mediator.change_current_autoboxer(str(self.col1[i].text()))
				self.update_ab_table()
			else:
				print "error, unforeseen checkstate circumstance. Nothing done"
			#except: pass
		
		self.lock = False
				
	def ab_table_item_changed(self,item):
		print "item changed"
	
	def update_ab_table(self):
		
		data = self.ab_sel_mediator.get_autoboxer_data()
		self.col1 = []
		self.col2 = []
		self.colnames = []
		if len(data) == 0:
			return
			
		self.lock = True
		idx = 0
		for d in data.items():
			if idx >= self.ab_table.rowCount():
				self.ab_table.insertRow(idx)
			col1 = QtGui.QTableWidgetItem(d[0])
			qstr =''
			for i,s in enumerate(d[1]):
				qstr += s
				if i != len(d[1])-1:
					qstr += ', '
			col2 = QtGui.QTableWidgetItem(qstr)
			
			flag1 = Qt.ItemFlags(Qt.ItemIsUserCheckable)
			flag2 = Qt.ItemFlags(Qt.ItemIsSelectable)
			flag3 = Qt.ItemFlags(Qt.ItemIsEnabled)
			flag4 = Qt.ItemFlags(Qt.ItemIsEditable)
			#flags = flags.
			col1.setFlags(flag1|flag2|flag3|flag4) #Qt.ItemIsEnabled+Qt.ItemIsUserCheckable)) #&Qt.ItemIsSelectable))
			col1.setCheckState( Qt.Unchecked)
			col2.setFlags(flag3)
			self.ab_table.setItem(idx,0,col1)
			self.ab_table.setItem(idx,1,col2)
			self.col1.append(col1)
			self.col2.append(col2)
			self.colnames.append(col1.text())
			idx += 1
		
		# remove any overhanging columns if they already existed
		while (len(data) < self.ab_table.rowCount()):
			self.ab_table.removeRow(self.ab_table.rowCount()-1)
	
		currentsel = self.ab_sel_mediator.get_current_autoboxer_ts()
		self.currentlychecked = -1
		
		self.setChecked(currentsel)
		self.lock = False
	
	def delete_autoboxer(self,unused):
		items = self.ab_table.selectedItems()
		data = self.ab_sel_mediator.get_autoboxer_data()
		update = False
		for item in items:
			if item.column() == 0:
				if len(self.ab_sel_mediator.associated_images(str(item.text()))) != 0:
					print "can't delete an autoboxer unless it has no images associated with it"
				else:
					self.ab_sel_mediator.remove(str(item.text()))
					update = True
					
		if update:
			self.update_ab_table()
	
	def add_new_autoboxer(self,bool):
		self.ab_sel_mediator.add_new_autoboxer()
		self.update_ab_table()
		
	def add_copy_autoboxer(self,bool):
		numsel = 0
		selected = ""
		for col in self.col1:
			if self.ab_table.isItemSelected(col):
				numsel += 1
				selected = str(col.text())
				if numsel > 1:
					print "error, more than one autoboxer is selected. Please choose only one"
					return
		if numsel == 0:
			print "no autoboxers were selected, doing nothing"
			return
		else:
			self.ab_sel_mediator.add_copy_autoboxer(selected)
			self.lock=False
			self.update_ab_table()
			
	def insert_view_tab(self):
		# this is the box layout that will store everything
		self.view_inspector = QtGui.QWidget()
		self.view_vbl =  QtGui.QVBoxLayout(self.view_inspector)
		
		#  Insert the plot widget
		self.viewhbl = QtGui.QHBoxLayout()
		
		self.viewboxes = QtGui.QCheckBox("Boxed Particle Window")
		self.viewboxes.setChecked(True)
		self.viewimage = QtGui.QCheckBox("Main Image Window")
		self.viewimage.setChecked(True)
		
		self.viewhbl.addWidget(self.viewboxes)
		self.viewhbl.addWidget(self.viewimage)
		
		if self.target().has_thumbnails():
			self.viewthumbs = QtGui.QCheckBox("Image Thumbnails Window")
			self.viewthumbs.setChecked(True)
			self.viewhbl.addWidget(self.viewthumbs)
		
		self.viewmanagement = QtGui.QGroupBox("Displayed Windows")
		self.viewmanagement.setLayout(self.viewhbl)
		self.view_vbl.addWidget(self.viewmanagement)
		
			
		self.viewhbl2 = QtGui.QHBoxLayout()
		self.boxdisplay =QtGui.QLabel("Box Display Object:",self)
		self.boxdisplay.setEnabled(False)
		self.viewhbl2.addWidget(self.boxdisplay)
		
		self.boxformats = QtGui.QComboBox(self)
		self.boxformats.addItem("square with central dot")
		self.boxformats.addItem("square")
		self.boxformats.addItem("circle with central dot")
		self.boxformats.addItem("circle")
		self.boxformats.setEnabled(False)
		self.viewhbl2.addWidget(self.boxformats)
		
		self.displayboxes = QtGui.QGroupBox("Displayed Boxes")
		self.displayboxes.setLayout(self.viewhbl2)
		self.view_vbl.addWidget(self.displayboxes)
		
		self.tabwidget.addTab(self.view_inspector,"Display Options")
		
		self.connect(self.viewboxes,QtCore.SIGNAL("clicked(bool)"),self.target().view_boxes_clicked)
		self.connect(self.viewimage,QtCore.SIGNAL("clicked(bool)"),self.target().view_image_clicked)
		if self.target().has_thumbnails():
			self.connect(self.viewthumbs,QtCore.SIGNAL("clicked(bool)"),self.target().view_thumbs_clicked)
			
		QtCore.QObject.connect(self.boxformats, QtCore.SIGNAL("currentIndexChanged(QString)"), self.box_format_changed)
	
	def set_guimxit_visible(self,val=True):
		self.lock = True
		self.viewthumbs.setChecked(val)
		self.lock = False
	
	def set_guimx_visible(self,val=True):
		self.lock = True
		self.viewboxes.setChecked(val)
		self.lock = False
	
	def set_guiim_visible(self,val=True):
		self.lock = True
		self.viewimage.setChecked(val)
		self.lock = False
	
	def box_format_changed(self,new_format):
		format = str(new_format)
		if format == "square with central dot": format = "rectpoint"
		elif format == "square": format = "rect"
		elif format == "circle with central dot": format = "rcirclepoint"
		elif format == "circle": format = "rcircle"
		else: 
			print "errror, unknown format" 
			return
		self.target().change_shapes(format)

	
	def insert_advanced_tab(self):
		# this is the box layout that will store everything
		self.adv_inspector = QtGui.QWidget()
		self.advanced_vbl =  QtGui.QVBoxLayout(self.adv_inspector)
		
		#Insert the plot widget
#		self.plothbl = QtGui.QHBoxLayout()
		
		# This was commented out because the extra GL context crashes the desktop
#		self.window = EMGLPlotWidget(self)
#		self.window.setInit()
#		self.window.resize(100,100)
#		self.window2=EMParentWin(self.window)
#		self.window2.resize(100,100)
#		
#		self.plothbl.addWidget(self.window)
		
#		self.plotbuttonvbl = QtGui.QVBoxLayout()
#		
#		self.trythat=QtGui.QPushButton("Try That")
#		self.plotbuttonvbl.addWidget(self.trythat)
#		
#		self.reset=QtGui.QPushButton("Reset")
#		self.plotbuttonvbl.addWidget(self.reset)
#		
#		self.plothbl.addLayout(self.plotbuttonvbl)
#		
		self.advanced_hbl2 = QtGui.QHBoxLayout()
		
#		self.advanced_vbl2.addLayout(self.plothbl)

		self.enable_interactive_params  = QtGui.QCheckBox("Enable")
		self.enable_interactive_params.setChecked(False) #FIXME should this be related to the state of the autoboxer?
		self.thr = ValSlider(self,(0.0,3.0),"Threshold:")
		self.thr.setValue(1.0)
		self.thr.setEnabled(False)
		self.advanced_hbl2.addWidget(self.enable_interactive_params)
		self.advanced_hbl2.addWidget(self.thr)


		self.interbox = QtGui.QGroupBox("Interactive Parameters")
		self.interbox.setLayout(self.advanced_hbl2)
		self.advanced_vbl.addWidget(self.interbox)
		
		self.thrbut = QtGui.QRadioButton(SwarmAutoBoxer.THRESHOLD)
		self.selbut = QtGui.QRadioButton(SwarmAutoBoxer.SELECTIVE)
		self.selbut.setChecked(True)
		self.morselbut = QtGui.QRadioButton(SwarmAutoBoxer.MORESELECTIVE)

		self.methodhbox = QtGui.QHBoxLayout()
		self.methodhbox.addWidget(self.thrbut)
		self.methodhbox.addWidget(self.selbut)
		self.methodhbox.addWidget(self.morselbut)

		self.groupbox = QtGui.QGroupBox("Auto Box Method")
		self.groupbox.setLayout(self.methodhbox)

		self.advanced_vbl.addWidget(self.groupbox)

		self.ratiobut = QtGui.QRadioButton("Ratio")

		self.difbut = QtGui.QRadioButton("Difference")
		self.ratio_average_but = QtGui.QRadioButton("Average Ratio")
		self.ratio_average_but.setChecked(True)

		self.cmpmethodhbox = QtGui.QHBoxLayout()
		self.cmpmethodhbox.addWidget(self.ratiobut)
		self.cmpmethodhbox.addWidget(self.difbut)
		self.cmpmethodhbox.addWidget(self.ratio_average_but)

		self.cmpgroupbox = QtGui.QGroupBox("Peak Profile Comparator")
		self.cmpgroupbox.setLayout(self.cmpmethodhbox)


		self.advanced_vbl.addWidget(self.cmpgroupbox)

		self.lock = True
		self.autoboxerhdbl = QtGui.QHBoxLayout()
		# ab means autoboxer
		self.ab_table = QtGui.QTableWidget(1,2,self)
		self.ab_table.setColumnWidth(1,150)
		self.abcol0title = QtGui.QTableWidgetItem("Autoboxer ID")
		self.abcol1title = QtGui.QTableWidgetItem("Associated Images")
		self.update_ab_table()
		self.lock = True
		self.ab_table.setHorizontalHeaderItem(0,self.abcol0title)
		self.ab_table.setHorizontalHeaderItem(1,self.abcol1title)
		self.autoboxerhdbl.addWidget(self.ab_table)
		self.lock = False

		self.autoboxervbl1 = QtGui.QVBoxLayout()
		self.abcopy = QtGui.QPushButton("Copy")
		self.autoboxervbl1.addWidget(self.abcopy)
		self.abnew = QtGui.QPushButton("New")
		self.autoboxervbl1.addWidget(self.abnew)
		self.abdelete = QtGui.QPushButton("Delete")
		self.autoboxervbl1.addWidget(self.abdelete)
		self.autoboxerhdbl.addLayout(self.autoboxervbl1)

		self.abmanagement = QtGui.QGroupBox("Auto Boxer Management")
		self.abmanagement.setLayout(self.autoboxerhdbl)
		self.advanced_vbl.addWidget(self.abmanagement)


		self.david_option = self.adv_inspector

		self.pawel_option = QtGui.QWidget()
		self.pawel_option_vbox = QtGui.QVBoxLayout(self.pawel_option)
		self.pawel_option_vbox.addWidget(QtGui.QLabel("Gauss Conv's Parameters") )
		pawel_grid1 = QtGui.QGridLayout( )
		self.pawel_option_vbox.addLayout(pawel_grid1)

		pawel_grid1.addWidget( QtGui.QLabel("Input Pixel Size:") , 0, 0 )
		pawel_grid1.addWidget( QtGui.QLabel("Output Pixel Size:"), 1, 0 )
		pawel_grid1.addWidget( QtGui.QLabel("Gauss Width Adjust:"), 2, 0 )
		pawel_grid1.addWidget( QtGui.QLabel("Angstrom"), 0, 2 )
		pawel_grid1.addWidget( QtGui.QLabel("Angstrom"), 1, 2 )

		self.use_variance = QtGui.QCheckBox("Use Variance Image")
		self.use_variance.setChecked(True)
		pawel_grid1.addWidget( self.use_variance, 3, 0)
		print " aaaa ",options.pix_in,options.pix_out,"%6.3f"%options.pix_out
		self.input_pixel_size  = QtGui.QLineEdit("%6.3f"%options.pix_in, self)
		self.output_pixel_size = QtGui.QLineEdit("%6.3f"%options.pix_out, self)
		pawel_grid1.addWidget( self.input_pixel_size, 0, 1 )
		pawel_grid1.addWidget( self.output_pixel_size, 1, 1 )
		print  " BB ",self.input_pixel_size.text(),self.output_pixel_size.text()

		self.gauss_width = QtGui.QLineEdit("1.0", self)
		self.gauss_width_slider = QtGui.QSlider(QtCore.Qt.Horizontal, self)
		pawel_grid1.addWidget( self.gauss_width_slider, 2, 1 )
		pawel_grid1.addWidget( self.gauss_width, 2, 2 )

		self.gauss_width_slider.setRange( -100, 100 )
		self.gauss_width_slider.setValue( 0 )


		# --------------------
		# XXX CTF button
		self.ctf_data = None
		self.i_start = None
		self.i_stop = None

		pawel_grid1.addWidget( QtGui.QLabel("Window size:"), 4,0 )
		self.ctf_window_size = QtGui.QLineEdit("%d"%options.ctf_window, self)
		pawel_grid1.addWidget( self.ctf_window_size , 4,1 )		

		pawel_grid1.addWidget( QtGui.QLabel("Edge size:"), 5,0 )
		self.ctf_edge_size = QtGui.QLineEdit("%d"%options.ctf_edge, self)	
		pawel_grid1.addWidget( self.ctf_edge_size , 5,1 )		

		pawel_grid1.addWidget( QtGui.QLabel("Overlap:"), 6,0 )
		self.ctf_overlap_size = QtGui.QLineEdit("%d"%options.ctf_overlap, self)
		pawel_grid1.addWidget( self.ctf_overlap_size , 6,1 )

		pawel_grid1.addWidget( QtGui.QLabel("F_stop:"), 8,0 )
		self.ctf_f_stop = QtGui.QLineEdit("%5.3f"%options.ctf_fstop, self)
		pawel_grid1.addWidget( self.ctf_f_stop , 8,1 )

		pawel_grid1.addWidget( QtGui.QLabel("F_start:"), 7,0 )
		self.ctf_f_start = QtGui.QLineEdit("%5.3f"%options.ctf_fstart, self)
		pawel_grid1.addWidget( self.ctf_f_start , 7,1 )

		pawel_grid1.addWidget( QtGui.QLabel("Cs:"), 4,2 )
		self.ctf_cs = QtGui.QLineEdit("%4.1f"%options.ctf_cs, self)
		pawel_grid1.addWidget( self.ctf_cs , 4,3 )	

		pawel_grid1.addWidget( QtGui.QLabel("Voltage:"), 5,2 )
		self.ctf_volt = QtGui.QLineEdit("%5.1f"%options.ctf_volt, self)
		pawel_grid1.addWidget( self.ctf_volt , 5,3 )		

		pawel_grid1.addWidget( QtGui.QLabel("Amplitude Contrast:"), 6,2 )
		self.ctf_ampcont = QtGui.QLineEdit("%5.1f"%options.ctf_ampcont, self)
		pawel_grid1.addWidget( self.ctf_ampcont , 6,3 )		

		self.ctf_button =  QtGui.QPushButton("Estimate CTF")
		pawel_grid1.addWidget( self.ctf_button, 9, 0)
		self.connect(self.ctf_button,QtCore.SIGNAL("clicked(bool)"), self.calc_ctf)

		self.inspect_button =  QtGui.QPushButton("Inspect CTF")
		pawel_grid1.addWidget( self.inspect_button, 9, 1)
		self.connect(self.inspect_button,QtCore.SIGNAL("clicked(bool)"), self.inspect_ctf)
		self.ctf_inspector = None
		self.ctf_inspector_gone = True
		# XXX
		# --------------------


		self.pawel_option_vbox.addWidget( QtGui.QLabel("CCF Histogram") )
		self.pawel_histogram = CcfHistogram( self )
		self.pawel_option_vbox.addWidget( self.pawel_histogram )

		self.ccf_range = QtGui.QLabel("Range: (N/A, N/A)")
		self.pawel_option_vbox.addWidget( self.ccf_range )

		pawel_grid2 = QtGui.QGridLayout()
		self.pawel_option_vbox.addLayout( pawel_grid2 )

		pawel_grid2.addWidget( QtGui.QLabel("Threshold Low:"),  0, 0 )
		pawel_grid2.addWidget( QtGui.QLabel("Threshold High:"), 1, 0 )
		self.threshold_low = QtGui.QLineEdit()
		self.threshold_hgh = QtGui.QLineEdit()
		pawel_grid2.addWidget( self.threshold_low, 0, 1 )
		pawel_grid2.addWidget( self.threshold_hgh, 1, 1 )

		#self.pawel_table.setVerticalHeaderLabels( ["Input Pixel Size  (Angstrom): ", "Output Pixel size (Angstrom): "] )
		#self.pawel_table.horizontalHeader().hide()
		#self.pawel_table.setItem(0, 0, QtGui.QTableWidgetItem("1.0") )
		#self.pawel_table.setItem(1, 0, QtGui.QTableWidgetItem("1.0") )
		#self.connect(self.pawel_table, QtCore.SIGNAL("cellChanged(int,int)"), self.pawel_parm_changed)


		self.tabwidget.addTab(self.david_option,"Swarm Advanced")
		self.connect(self.abnew, QtCore.SIGNAL("clicked(bool)"), self.add_new_autoboxer)
		self.connect(self.abcopy, QtCore.SIGNAL("clicked(bool)"), self.add_copy_autoboxer)
		self.connect(self.abdelete, QtCore.SIGNAL("clicked(bool)"), self.delete_autoboxer)
		self.connect(self.ab_table, QtCore.SIGNAL("itemChanged(QtGui.QTableWidgetItem)"), self.ab_table_item_changed)
		self.connect(self.ab_table, QtCore.SIGNAL("cellChanged(int,int)"), self.ab_table_cell_changed)
		self.connect(self.selbut, QtCore.SIGNAL("clicked(bool)"), self.selection_mode_changed)
		self.connect(self.morselbut, QtCore.SIGNAL("clicked(bool)"), self.selection_mode_changed)
		self.connect(self.ratiobut, QtCore.SIGNAL("clicked(bool)"), self.cmp_box_changed)
		self.connect(self.ratio_average_but, QtCore.SIGNAL("clicked(bool)"), self.cmp_box_changed)
		#self.connect(self.centerbutton,QtCore.SIGNAL("clicked(bool)"),self.centerpushed)
		self.connect(self.difbut, QtCore.SIGNAL("clicked(bool)"), self.cmp_box_changed)
		self.connect(self.gauss_width_slider, QtCore.SIGNAL("valueChanged(int)"), self.gauss_width_changed)
		self.connect(self.gauss_width, QtCore.SIGNAL("editingFinished()"), self.gauss_width_edited)
		self.connect(self.enable_interactive_params,QtCore.SIGNAL("clicked(bool)"),self.enable_interactive_params_toggled)



	# --------------------
	# XXX: calculate ctf and defocus from current image		
	def inspect_ctf(self):
		#display(self.ctf_data)
		
		if not(self.ctf_inspector):
			self.ctf_inspector = CTFInspector(self,self.ctf_data)
			self.ctf_inspector.show()
			self.ctf_inspector_gone=False
		else:
			if (self.ctf_inspector_gone):
				self.ctf_inspector.show()
				self.ctf_inspector_gone=False
			else:
				pass
	
	def calc_ctf(self):
		# calculate power spectrum of image with welch method (welch_pw2)
		# calculate rotational average of power spectrum (rot_avg_table)
		# calculate ctf values with ctf_get
		print "starting CTF estimation"

		# get the current image
		image_name = self.target().boxable.get_image_name()
		img = BigImageCache.get_image_directly( image_name )

		# conversion from text necessary
		try:
			ctf_window_size  = int(self.ctf_window_size.text())
			ctf_edge_size    = int(self.ctf_edge_size.text())
			ctf_overlap_size = int(self.ctf_overlap_size.text())
			ctf_f_start      = float(self.ctf_f_start.text())
			ctf_f_stop       = float(self.ctf_f_stop.text())
			ctf_volt         = float(self.ctf_volt.text())
			ctf_cs           = float(self.ctf_cs.text())
			ctf_ampcont      = float(self.ctf_ampcont.text())
		except ValueError,extras:
			# conversion of a value failed!
			print "integer conversion failed."
			if not(extras.args is None):
				print extras.args[0]
			return
		except:
			print "error"
			print self.ctf_window_size.text()
			print self.ctf_overlap_size.text()
			print self.ctf_edge_size.text()
			return

		# print "determine power spectrum"
		from fundamentals import welch_pw2
		# XXX: check image dimensions, especially box size for welch_pw2!
		power_sp = welch_pw2(img, win_size=ctf_window_size, overlp_x=ctf_overlap_size, overlp_y=ctf_overlap_size,
				     edge_x=ctf_edge_size, edge_y=ctf_edge_size)

		# print "averaging power spectrum"
		from fundamentals import rot_avg_table
		avg_sp = rot_avg_table(power_sp)
		del power_sp

		# print "determine ctf"
		from morphology import defocus_gett

		px_size = float(self.output_pixel_size.text())
		print "  pixel   ",px_size

		defocus = defocus_gett(avg_sp, voltage=ctf_volt, Pixel_size=px_size, Cs=ctf_cs, wgh=ctf_cs,
				       f_start=ctf_f_start, f_stop=ctf_f_stop, parent=self)
		del avg_sp

		print "CTF estimation done"
		print "Defocus: ",defocus

		# update ctf inspector values
		if (self.ctf_inspector is not None):
			self.ctf_inspector.setData(self.ctf_data)
			self.ctf_inspector.i_start = self.i_start
			self.ctf_inspector.i_stop = self.i_stop
			if not(self.ctf_inspector_gone):
				self.ctf_inspector.update()

		# XXX: wgh?? amp_cont static to 0?
		# set image properties, in order to save ctf values
		from utilities import set_ctf
		set_ctf(img,[defocus,ctf_cs,ctf_volt,px_size,0,ctf_ampcont])
		# and rewrite image 
		img.write_image(image_name)
		del img,image_name
		
	def closeEvent(self,event):
		
		if ((self.ctf_inspector is not None) and not(self.ctf_inspector_gone)):
			del self.ctf_inspector
			self.ctf_inspector = None
		# someone overwrote my closeEvent! oh well gotta do this:
		self.target().done()
		event.accept()
		
	# XXX
	# --------------------


	def pawel_parm_changed(self, row, col ):
		from string import atof

		assert col==0

		t = self.pawel_table.item(row, col).text()

		print 'row, col, text: ', row, col, t

		if row==0:
			self.target().autoboxer.pixel_input = atof( t )
		else:
			assert row==1
			self.target().autoboxer.pixel_output = atof( t )		

	def gauss_width_edited(self):
		from string import atof
		from math import log10
		text = self.gauss_width.text()
		v = int( log10(atof(text)) * 100)
		self.gauss_width_slider.setValue( v )

	def gauss_width_changed(self, v):
		from math import pow
		s = "%.3f" % pow(10.0, v*0.01)
		self.gauss_width.setText( s )
		

	def set_dynapix(self,bool):
		self.dynapix.setChecked(bool)

	def ref_button_toggled(self,bool):
		
		if self.refbutton.isChecked():
			self.manualbutton.setChecked(False)
		
		if not self.refbutton.isChecked():
			self.manualbutton.setChecked(True)

		self.target().set_boxing_method(self.refbutton.isChecked(),self.manualbutton.isChecked())
		
	def manual_button_toggled(self,bool):
		if self.manualbutton.isChecked():
			self.refbutton.setChecked(False)
		
		if not self.manualbutton.isChecked():
			self.refbutton.setChecked(True)
			
		self.target().set_boxing_method(self.refbutton.isChecked(),self.manualbutton.isChecked())

	def erase_toggled(self,bool):
		self.unerase.setChecked(False)
		self.eraserad.setEnabled(bool)
		#self.target.get_2d_gui_image().get_qt_context_parent().setMouseTracking(bool)
		self.target().erase_toggled(bool)
		
	def unerase_toggled(self,bool):
		self.erase.setChecked(False)
		self.eraserad.setEnabled(bool)
		#self.target.get_2d_gui_image().get_qt_context_parent().setMouseTracking(bool)
		self.target().unerase_toggled(bool)
		
	def dynapix_toggled(self,bool):
		self.target().toggle_dynapix(bool)
	
	def cmp_box_changed(self,unusedbool):
		if self.ratiobut.isChecked():
			s = BoxingTools.CmpMode.SWARM_RATIO
		elif self.difbut.isChecked():
			s = BoxingTools.CmpMode.SWARM_DIFFERENCE
		elif self.ratio_average_but.isChecked():
			s = BoxingTools.CmpMode.SWARM_AVERAGE_RATIO
		else:
			print "Bug intercepted in e2boxer.py. Please email the development team."
			
		self.target().set_profile_comparitor(s)
	
	def selection_mode_changed(self,unusedbool):
		if self.thrbut.isChecked():
			s = self.thrbut.text()
		elif self.selbut.isChecked():
			s = self.selbut.text()
		elif self.morselbut.isChecked():
			s = self.morselbut.text()
		else:
			print "Bug intercepted in e2boxer.py. Please email the development team."
			
		self.target().set_selection_mode(str(s))
	
	def try_dummy_parameters(self):
		#self.dummybox.set_opt_profile(self.window.getData())
		self.dummybox.set_correlation_score(float(self.thr.getValue()))
		self.target().set_dummy_box(self.dummybox)
		
	def enable_interactive_params_toggled(self):
		enabled = self.enable_interactive_params.isChecked()
		if enabled:
			self.thr.setEnabled(True)
			self.try_dummy_parameters()
		else:
			self.thr.setEnabled(False)
			self.target().set_dummy_box(None)
			
			
	def update_data(self,thresh,data,datar):
		#print data
#		self.window.set_data(data,datar)
		self.thr.setValue(thresh,True)
#		self.resize(self.width(),self.height())
		#self.window.resizeGL(self.window.width(),self.window.height())
		#self.window.updateGL()
	def num_boxes_changed(self,n):
		self.info.setText("%d Boxes"%n)
		
	#def ppc_changed(self,f):
		#self.ppc.setText("%f ppc"%f)
	
	def adjust_erase_rad(self,delta):
		v = float(self.eraserad.text())
		if delta > 0:
			v = 1.1*v
		if delta < 0:
			v = 0.9*v
			
		self.eraserad.setText(str(int(v)))
		# this makes sure the target updates itself 
		# there may be a better approach, seeing as
		# the target called this function
		self.update_erase_rad()
		
	def update_erase_rad(self):
		v = int(self.eraserad.text())
		if ( v < 1 ): raise Exception
		self.target().update_erase_rad(v)
	
	def new_box_size(self):
		try:
			v=int(self.bs.text())
			if v<12 : raise Exception
		except:
			self.bs.setText(str(self.target().box_size))
			return
		
		get_application().setOverrideCursor(Qt.BusyCursor)
		self.target().update_box_size(v,1)
		get_application().setOverrideCursor(Qt.ArrowCursor)
	
	def set_box_size(self,box_size):
		self.bs.setText(str(box_size))
	
	def new_threshold(self,val):
		self.try_dummy_parameters()


# --------------------
# XXX: class for popup display of CTF. this is a rough prototype meant to be extended for full
#    display of values and CTF.

class CTFInspector(QtGui.QWidget):

	def __init__(self,parent,data=None) :
		QtGui.QWidget.__init__(self) 
		# we need to keep track of our parent to signal when we are gone again....
		self.parent = weakref.ref(parent) # this needs to be a weakref ask David Woolford for details, but otherwise just call self.parent() in place of self.parent
		self.setGeometry(300, 300, 250, 150)
		self.setWindowTitle("CTF Inspector")

		self.i_start = None
		self.i_stop = None

		if (data is None):
			# test data, to ensure something is displayed even if no data is set yet. this is
			#    for development only and can be removed later.....
			self.data = [[80,20,10,9,8,7,6,5,4,3,2,1,0,0,0,0,0],]
		else:
			# assume we got a triple of lists. assign it for now.
			self.data=data


	def setData(self,data):
		# check data type is a list and break, if not
		if not(type(data) is list):
			return False
		
		# free previous and reset our data to the passed list
		del self.data
		self.data = data
		# return success
		return True

	def update(self):
		QtGui.QWidget.update(self) #self.paintEvent(None)
		# print "update..."
		
	def paintEvent(self,event):
		h=self.height()
		w=self.width()

		hborder = ( min((h / 15.0),20.0))
		wborder = ( min((w / 15.0),20.0))

		# accessible height and width....
		ah = int(h-2*hborder)
		aw = int(w-2*wborder)

		p=QtGui.QPainter()
		p.begin(self)
		p.setBackground(QtGui.QColor(16,16,16))
		p.eraseRect(0,0,self.width(),self.height())
		p.setPen(Qt.white)

		# labels
		# spectrum
		# background
		# ctf

		# draw axes
		p.drawLine(int(wborder),int(hborder),int(wborder),int(h-hborder))
		p.drawLine(int(wborder),int(h-hborder),int(w-wborder),int(h-hborder))

		color = [Qt.yellow,Qt.red,Qt.blue]
		labels= ["Roo","Back","CTF"]

		if (not(self.data == []) and not(self.data is None)):

			# scaling factors in x and y, respectively. margins are left around the plot,
			#    stepw along x and 10% of height in y... explicit conversion is necessary,
			#    since we are working with ints....			
			if ((self.i_start is not None) and (self.i_stop is not None)):
				sizew = self.i_stop - self.i_start + 1
			else:
				sizew = max([len(i) for i in self.data])
				self.i_start = 0
				self.i_stop = sizew-1
				sizew=float(sizew)

			# print "range: ",self.i_start," - ",self.i_stop
				
			stepw = float(w-2*wborder) / float(sizew)

			if ((self.i_start is not None) and (self.i_stop is not None)):
				sizeh = max([max(self.data[i][self.i_start:self.i_stop]) for i in xrange(len(self.data))])
			else:
				sizeh = max([max(self.data[i]) for i in xrange(len(self.data))])
				
			sizeh = float(sizeh)
			steph = float(h-2*hborder) / float(sizeh) 

			for list_index in xrange(len(self.data)):
				
				p.setPen(color[list_index])
				metrics = p.fontMetrics()
				fw = metrics.width(str(labels[list_index]))
				fh = metrics.height()+4
				p.drawText(w-wborder-fw/2, hborder+(list_index)*fh, str(labels[list_index]))
				
				for index in xrange(self.i_start,self.i_stop):
					p.setPen(color[list_index])
					# skip first point, since there is no previous point to connect to
					if (0 == index):
						continue
					else:
						# x is normal, y is flipped (i.e. top left is (0,0))
						#oldx = int(wborder+ (stepw*(index-1)))
						oldx=int(wborder + ((w-2*wborder) / sizew * (index-1-self.i_start)))
						#newx = int(wborder+ (stepw*(index)))
						newx=int(wborder + ((w-2*wborder) / sizew * (index-self.i_start)))
						
						#oldy = int(h-hborder-steph*self.data[list_index][index-1])
						oldy=int(h-hborder-(h-2*hborder)*self.data[list_index][index-1]/sizeh)
						#newy = int(h-hborder-steph*self.data[list_index][index])
						newy=int(h-hborder-(h-2*hborder)*self.data[list_index][index]/sizeh)
						p.drawLine(oldx,oldy,newx,newy)
					if (len(self.data)-1 == list_index):
						p.setPen(Qt.white)
						p.drawLine(newx, h-hborder, newx, h-hborder-5)
						metrics = p.fontMetrics()
						fw = metrics.width(str(index))
						p.drawText(newx-fw/2, h-hborder-7, str(index))
						
		p.end()

	# closing the window is tricky: we need to notify the parent window we are gone, but
	#    cannot set parent.ctf_inspector directly, since that would destroy ourselves in
	#    the middle of the event handler, prompting an error. instead, we set a flag in
	#    the parent object and let it handle destroying, resetting or updating when
	#    it becomes necessary....
	def closeEvent(self,event):
		# set the flag of our parent object
		self.parent().ctf_inspector_gone=True
		# and close ourselves by accepting the event....
		event.accept()

# XXX
# --------------------

if __name__ == "__main__":
	main()
