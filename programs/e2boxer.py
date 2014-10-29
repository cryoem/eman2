#!/usr/bin/env python
#
# Author: David Woolford (woolford@bcm.edu)
# Copyright (c) 2000- Baylor College of Medicine


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

from EMAN2 import BoxingTools,gm_time_string,Transform, E2init, E2end, E2progress,db_open_dict,EMArgumentParser
from EMAN2db import db_check_dict
from EMAN2jsondb import *
from pyemtbx.boxertools import CoarsenedFlattenedImageCache,FLCFImageCache
from copy import deepcopy
from EMAN2 import *
from emboxerbase import *
import os

SWARM_TEMPLATE_MIN = TEMPLATE_MIN # this comes from emboxerbase

def e2boxer_check(options,args):
	error_message = check(options,args)

	if options.autoboxer:
		if not js_check_dict("e2boxercache/swarm.json"): error_message.append("There is no autoboxing information present in the current directory")
		else:
			db = js_open_dict("e2boxercache/swarm.json")
			if not db.has_key(options.autoboxer):
				s = "There is no autoboxing information present for %s." %options.autoboxer
				if len(db.keys()) > 0:
					s+= ("Choose from")
					for k in db.keys():
						try: s+=" "+str(k)+";"
						except: pass
				error_message.append(s)
			else: pass #print "yaya"
#				print db.keys(),options.autoboxer
#				for k,v in db.items(): print k,v
#				print db[options.autoboxer]
#				print db[str(options.autoboxer)]

	return error_message

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <image> <image2>....

The 'new' version of e2boxer. While it is quite usable, it is still in need of some work.

For example:

e2boxer.py ????.mrc --boxsize=256
"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_pos_argument(name="micrographs",help="List the file to process with e2boxer here.", default="", guitype='filebox', browser="EMBoxesTable(withmodal=True,multiselect=True)",  row=0, col=0,rowspan=1, colspan=3, mode="boxing,extraction")
	parser.add_header(name="boxerheader", help='Options below this label are specific to e2boxer', title="### e2boxer options ###", row=1, col=0, rowspan=1, colspan=3, mode="boxing,extraction")
	parser.add_argument("--boxsize","-B",type=int,help="Box size in pixels",default=-1, guitype='intbox', row=2, col=0, rowspan=1, colspan=3, mode="boxing,extraction")
	parser.add_argument("--write_dbbox",action="store_true",help="Write coordinate file (eman1 dbbox) files",default=False, guitype='boolbox', row=3, col=0, rowspan=1, colspan=1, mode="extraction")
	parser.add_argument("--write_ptcls",action="store_true",help="Write particles to disk",default=False, guitype='boolbox', row=3, col=1, rowspan=1, colspan=1, mode="extraction[True]")
	parser.add_argument("--exclude_edges",action="store_true",help="Don't generate output for any particles extending outside the micrograph",default=False, guitype='boolbox', row=4, col=1, rowspan=1, colspan=1, mode="extraction[True]")
	parser.add_argument("--force","-f",action="store_true",help="Force overwrite",default=False, guitype='boolbox', row=4, col=0, rowspan=1, colspan=1, mode="extraction")
	parser.add_argument("--format", help="Format of the output particle images. For EMAN2 refinement must be HDF.", default="hdf", guitype='combobox', choicelist="['hdf','img','spi']", row=6, col=0, rowspan=1, colspan=2, mode="extraction")
	parser.add_argument("--norm", type=str,help="Normalization processor to apply to written particle images. Should be normalize, normalize.edgemean,etc.Specifc \"None\" to turn this off", default="normalize.edgemean", guitype='combobox', choicelist='re_filter_list(dump_processors_list(),\'normalize\')', row=5, col=0, rowspan=1, colspan=2, mode="extraction")
	parser.add_argument("--invert",action="store_true",help="If writing outputt inverts pixel intensities",default=False, guitype='boolbox', row=3, col=2, rowspan=1, colspan=1, mode="extraction")
	parser.add_argument("--suffix",type=str,help="suffix which is appended to the names of output particle and coordinate files",default="_ptcls", guitype='strbox', expert=True, row=4, col=1, rowspan=1, colspan=2, mode="extraction")
	parser.add_argument("--dbls",type=str,help="data base list storage, used by the workflow. You can ignore this argument.",default=None)
	parser.add_argument("--autoboxer",type=str,help="A key of the swarm_boxers dict in the local directory, used by the workflow.",default=None)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--gui", action="store_true", default=True, help="Dummy option; used in older version of e2boxer")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
 	parser.add_argument("--gauss_autoboxer",type=str,help="Name of autoboxed file whose autoboxing parameters (obtained via some previous run of Gauss autoboxer via the GUI) should be used for automatic boxing.",default=None)
	parser.add_argument("--do_ctf",type=str,help="Name of file whose ctf estimation parameters (obtained via some previous run of Gauss autoboxer via the GUI) should be used for automatic ctf estimation.",default=None)

	# ctf estimation using cter
	parser.add_argument("--cter",               action="store_true",     default=False,                  help="CTF estimation using cter")
	parser.add_argument("--indir",              type=str,				 default= ".",     				 help="Directory containing micrographs to be processed.")
	parser.add_argument("--nameroot",         	type=str,		 		 default= "",     				 help="Prefix of micrographs to be processed.")
	parser.add_argument("--micsuffix",          type=str,                default="",                     help="A string denoting micrograph type. For example 'mrc', 'hdf', 'ser' ...")
	parser.add_argument("--wn",				  	type=int,		 		 default=256, 					 help="Size of window to use (should be slightly larger than particle box size)")

	parser.add_argument("--Cs",               	type=float,	 			 default= 2.0,               	 help="Microscope Cs (spherical aberation)")
	parser.add_argument("--voltage",		  	type=float,	 		     default=300.0, 			     help="Microscope voltage in KV")
	parser.add_argument("--ac",				  	type=float,	 		     default=10.0, 					 help="Amplitude contrast (percentage, default=10)")
	parser.add_argument("--kboot",			  	type=int,			     default=16, 					 help="kboot")
	parser.add_argument("--MPI",                action="store_true",   	 default=False,              	 help="use MPI version")
	parser.add_argument("--debug",              action="store_true",   	 default=False,              	 help="debug mode")
	parser.add_argument("--apix",               type=float,				 default= -1.0,                  help="pixel size in Angstroms")

	(options, args) = parser.parse_args()

	if options.cter:
		from global_def import SPARXVERSION
		import global_def
		if len(args) <2 or len(args) > 3:
			print "see usage"
			sys.exit()

		stack = None

		if len(args) == 3:
			if options.MPI:
				print "Please use single processor version if specifying a stack."
				sys.exit()
			stack = args[0]
			out1 = args[1]
			out2 = args[2]

		if len(args) == 2:
			out1 = args[0]
			out2 = args[1]

		if options.apix < 0:
			print "Enter pixel size"
			sys.exit()

		if options.MPI:
			from mpi import mpi_init, mpi_finalize
			sys.argv = mpi_init(len(sys.argv), sys.argv)

		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()

		from morphology import cter
		global_def.BATCH = True
		cter(stack, out1, out2, options.indir, options.nameroot, options.micsuffix, options.wn, \
			voltage=300.0, Pixel_size=options.apix, Cs = options.Cs, wgh=options.ac, kboot=options.kboot, MPI=options.MPI, DEBug = options.debug)
		global_def.BATCH = False

		if options.MPI:
			from mpi import mpi_finalize
			mpi_finalize()

		return

	logid=E2init(sys.argv,options.ppid)
	db = js_open_dict(EMBOXERBASE_DB)
	cache_box_size = True
	if options.boxsize == -1:
		cache_box_size = False
		options.boxsize = db.setdefault("box_size",128)

	error_message = e2boxer_check(options,args)
	if len(error_message) > 0:
		error = "\n"
		for e in error_message:
			error += "Error: "+e +"\n"
		parser.error(error)

	if cache_box_size: db["box_size"] = options.boxsize

	if options.autoboxer:
		autobox(args,options,logid)
		return

	if options.gauss_autoboxer or options.do_ctf:
		err = gauss_cmd_line_autobox(args, options,logid)
		if len(err) > 0:
			print err
			return
		if options.write_ptcls or options.write_dbbox:
			write_output(args,options,logid)
		return

	if options.write_dbbox or options.write_ptcls:
		write_output(args,options,logid)
	else:
		application = EMApp()
		module = EMBoxerModule(args,options.boxsize)


		# this is an example of how to add your own custom tools:
		module.add_tool(SwarmTool,particle_diameter=options.boxsize)

		#Make the Swarm Tool the default
		swarm_tool_instance = SwarmTool(object)
		swarm_tool_name = swarm_tool_instance.unique_name()
		module.set_inspector_tool_mode(swarm_tool_name) #Set current tool to Swarm mode


		module.add_tool(GaussTool,particle_diameter=options.boxsize)
		gauss_tool_instance = GaussTool(object)
		gauss_tool_name = gauss_tool_instance.unique_name()
		(module.tools[gauss_tool_name]).get_alternate(args[0])

		module.show_interfaces()

		application.execute(logid)

def gauss_cmd_line_autobox(args,options,logid):
	err = ''
	boxkey = options.gauss_autoboxer
	ctfkey = options.do_ctf
	if boxkey == None:
		boxkey = ctfkey
	if ctfkey == None:
		ctfkey = boxkey
	gdb_name = GaussPanel.GDB_NAME
	if not js_check_dict(gdb_name):
		err = "There is no gausse method autoboxing information present in the current directory"
		return err
	gbdb = js_open_dict(gdb_name)
	if (not gbdb.has_key(boxkey)):
		err = "There is no autoboxing information present for %s." %boxkey
		return err

	gboxer = GaussBoxer()
	if options.boxsize == -1:
		print "box size must be set! Exiting ..."
		sys.exit(1)
	boxsize=int(options.boxsize)
	do_autobox=False

	if options.gauss_autoboxer:
		do_autobox=True
		try:
			autoboxdict = gbdb[boxkey]
			gboxer.set_pixel_input(float(autoboxdict['pixel_input']))
			gboxer.set_pixel_output(float(autoboxdict['pixel_output']))
			gboxer.set_invert_contrast(autoboxdict['invert_contrast'])
			gboxer.set_use_variance(autoboxdict['use_variance'])
			gboxer.set_gauss_width(float(autoboxdict['gauss_width']))
			gboxer.set_thr_low(float(autoboxdict['thr_low']))
			gboxer.set_thr_hi(float(autoboxdict['thr_hi']))
		except:
			print 'no autoboxing parameters were found in database...will not autobox'
			do_autobox=False # maybe user just want to estimate ctf in batch mode...

	do_ctfest = False
	if options.do_ctf:
		do_ctfest = True
		ctfdict = gbdb[ctfkey]
		try:
			ctf_params={'pixel_input': float(ctfdict['pixel_input']), 'pixel_output': float(ctfdict['pixel_output']),'ctf_fstart':float(ctfdict['ctf_fstart']),'ctf_fstop':float(ctfdict['ctf_fstop']), 'ctf_window':int(ctfdict['ctf_window']),'ctf_edge':int(ctfdict['ctf_edge']),'ctf_overlap':int(ctfdict['ctf_overlap']),'ctf_ampcont':float(ctfdict['ctf_ampcont']),'ctf_volt':float(ctfdict['ctf_volt']),'ctf_cs':float(ctfdict['ctf_cs'])}
		except:
			print 'no ctf parameters stored in current directory! Estimate ctf in the GUI in gauss mode first!'
			do_ctfest = False

	for i,arg in enumerate(args):
		boxer_vitals = EMBoxerModuleVitals(file_names=[arg], box_size=boxsize)
		boxer_vitals.current_idx=0
		gboxer.target = boxer_vitals
		# Estimate ctf first and then autobox so ctf parametes need not be separately set in alternate image (i.e., the filtered and possibly downsampled image)
		if do_ctfest:
			gboxer.auto_ctf(arg,ctf_params)
		if do_autobox:
			gboxer.auto_box_cmdline(arg,boxsize=boxsize)
		E2progress(logid,float(i+1)/len(args))
	return err

def autobox(args,options,logid):

	boxer = SwarmBoxer()
	boxer.handle_file_change(options.autoboxer)
	for i,arg in enumerate(args):
		boxer_vitals = EMBoxerModuleVitals(file_names=[arg])
		boxer_vitals.current_idx=0
		boxer.target = weakref.ref(boxer_vitals)
		boxer.auto_box(arg, False, True, True)
		E2progress(logid,float(i+1)/len(args))

def write_output(args,options,logid, database="e2boxercache"):
	params = {}
	params["filenames"] = args
	params["suffix"] = options.suffix
	params["format"] = options.format
	db = js_open_dict(database+"/quality.json")
	db['suffix'] = options.suffix
	db['extension'] = os.path.splitext(args[0])[-1]

	total_progress = 0
	if options.write_ptcls:total_progress += len(args)
	if options.write_dbbox:total_progress += len(args)
	progress = 0.0
	E2progress(logid,0.0)

	if options.write_ptcls:
		if options.gauss_autoboxer != None:
			print "\n\n --write_ptcl option has been deactivated for Gauss mode.\n\nPlease use sxwindow.py for windowing!\n\n"
		else:
			names = get_particle_outnames(params)
			for i,output in enumerate(names):
				input = args[i]
				box_list = EMBoxList()
				box_list.load_boxes_from_database(input)
	
				# if box type is GaussBoxer.AUTO_NAME, the pre-process and possibly decimate image using params in db
				# only need to do this if write_ptcls is called on its own
				if (len(box_list) > 0) and (options.gauss_autoboxer == None):
					bx = box_list[0]
					if bx.type == GaussBoxer.AUTO_NAME and options.gauss_autoboxer == None:
						print 'For automatic boxing using Gauss Boxer, particles should be written at the time of autoboxing either by 1) activating "Write Output" button in GUI mode, or 2) adding --write_ptcls when autoboxing in command line mode'
						print 'I will continue and write particles, but results may not be as expected.'
	
				# if box type is GaussBoxer.AUTO_NAME, the pre-process and possibly decimate image using params in db
				# only need to do this if write_ptcls is called on its own
				if (len(box_list) > 0) and (options.gauss_autoboxer != None):
					bx = box_list[0]
					boxkey = options.gauss_autoboxer
					if bx.type == GaussBoxer.AUTO_NAME:
						# these were boxed in gauss mode, so we have to process and possibly decimate image before writing
	
						gdb_name = GaussPanel.GDB_NAME
						if not js_check_dict(gdb_name):
							print "no gauss mode autoboxing parameters were found in database...this should not happen"
							print 'exiting...'
							return
						gbdb = js_open_dict(gdb_name)
						if not gbdb.has_key(boxkey):
							print "no gauss mode autoboxing parameters were found in database for %s...this should not happen"%boxkey
							print 'exiting...'
							return
						autoboxdict = gbdb[boxkey]
	
						gboxer = GaussBoxer()
						boxsize=int(options.boxsize)
						try:
							gboxer.set_pixel_input(float(autoboxdict['pixel_input']))
							gboxer.set_pixel_output(float(autoboxdict['pixel_output']))
							gboxer.set_invert_contrast(autoboxdict['invert_contrast'])
							gboxer.set_use_variance(autoboxdict['use_variance'])
							gboxer.set_gauss_width(float(autoboxdict['gauss_width']))
							gboxer.set_thr_low(float(autoboxdict['thr_low']))
							gboxer.set_thr_hi(float(autoboxdict['thr_hi']))
						except:
							print 'no gauss mode autoboxing parameters were found in database...this should not happen'
							print 'exiting...'
							return
	
						gboxer.get_small_image(input,modecmd=True,boxsize=boxsize,ret_dummy=True)
	
				box_list.write_particles(input,output,options.boxsize,options.invert,options.norm,options.exclude_edges)
	
	
				progress += 1.0
				E2progress(logid,progress/total_progress)

	if options.write_dbbox:
		names = get_coord_outnames(params)

		for i,output in enumerate(names):
			input = args[i]
			box_list = EMBoxList()
			box_list.load_boxes_from_database(input)
			box_list.write_coordinates(input,output,options.boxsize) # input is redundant but it makes output interfaces generic

			progress += 1.0
			E2progress(logid,progress/total_progress)

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
	ave = ptcls[0].process("xform.centeracf")
	for i in range(1,len(ptcls)):
		ave.add(ptcls[i].process("xform.centeracf"))
#		ta = ptcls[i].align("translational",ave)
#		ave.add(ta)

	ave.process_inplace("xform.centeracf")
	ave.process_inplace("math.rotationalaverage")
	ave.process_inplace("normalize.edgemean")
	ave.process_inplace("mask.sharp",{'outer_radius':ave.get_xsize()/2})
	averages.append(ave)
	averages[-1].set_attr("creation_time_stamp", gm_time_string())
	for n in range(0,iter):
		t = []
		for idx,particle in enumerate(ptcls):
			ta = particle.align("translational",ave)
			t.append(ta)

		ave = t[0].copy()
		for i in range(1,len(ptcls)):
			ave.add(t[i])

		ave.process_inplace("xform.centeracf")
		ave.process_inplace("math.rotationalaverage")
		ave.process_inplace("normalize.edgemean")

		# edge normalize here SL before
		ave.process_inplace("mask.sharp",{'outer_radius':ave.get_xsize()/2})
		averages.append(ave)
		averages[-1].set_attr("creation_time_stamp", gm_time_string())

	return averages


class SwarmShrunkenImageMediator:
	def __init__(self,template_radius,subsample_rate):
		self.template_radius = template_radius
		self.subsample_rate = subsample_rate

	def get_template_radius(self): return self.template_radius
	def get_subsample_rate(self): return self.subsample_rate

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
		self.ever_moved = False # used to record if the user ever moved the ref after they first added it

	def type_name(self):
		if self.in_template: return SwarmBoxer.REF_NAME
		else: return SwarmBoxer.WEAK_REF_NAME

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


class SwarmPanel:
	DB_NAME = "e2boxercache/swarm_panel.json"

	def __init__(self,target,particle_diameter=128):
		self.busy = True
		self.particle_diameter = particle_diameter
		self.target = weakref.ref(target)
		self.widget = None
		self.busy = False

	def get_widget(self):
		if self.widget == None:
			from PyQt4 import QtCore, QtGui, Qt
			self.widget = QtGui.QWidget()
			vbl = QtGui.QVBoxLayout(self.widget)
			vbl.setMargin(0)
			vbl.setSpacing(6)
			vbl.setObjectName("vbl")

			db = js_open_dict(SwarmPanel.DB_NAME)

			hbl = QtGui.QHBoxLayout()
			ptcl_diam_label = QtGui.QLabel("Particle Diameter:")
			ptcl_diam_label.setToolTip("An estimate of the particle diameter - this used by Swarm for automatically shrinking and for determining automatic picking parameters.\nA value that is slightly larger than the particle is generally good. Err on the side of being too large, not too small.")
			hbl.addWidget(ptcl_diam_label)

			self.ptcl_diam_edit = QtGui.QLineEdit(str(self.particle_diameter))
			hbl.addWidget(self.ptcl_diam_edit)
			self.clear=QtGui.QPushButton("Clear Boxes")
			self.clear.setToolTip("Clear boxes generated by the Swarm mode.")
			hbl.addWidget(self.clear)

			vbl.addLayout(hbl)

			self.thrbut = QtGui.QRadioButton(SwarmBoxer.THRESHOLD)
			self.selbut = QtGui.QRadioButton(SwarmBoxer.SELECTIVE)
			self.selbut.setChecked(True)
			self.morselbut = QtGui.QRadioButton(SwarmBoxer.MORESELECTIVE)

			self.method_group = QtGui.QButtonGroup()
			self.method_group.addButton(self.thrbut)
			self.method_group.addButton(self.selbut)
			self.method_group.addButton(self.morselbut)

			self.methodhbox = QtGui.QHBoxLayout()
			self.methodhbox.addWidget(self.thrbut)
			self.methodhbox.addWidget(self.selbut)
			self.methodhbox.addWidget(self.morselbut)

			self.groupbox = QtGui.QGroupBox("Auto Box Method")
			self.groupbox.setToolTip("Tell Swarm what criteria to use for selecting particles.")
			self.groupbox.setLayout(self.methodhbox)
			vbl.addWidget(self.groupbox)

			hbl_ww = QtGui.QHBoxLayout()
			self.view_template=QtGui.QPushButton( QtGui.QIcon(get_image_directory() + "pp_boxer_icon.png"),"View Template")
			self.view_template.setEnabled(False)
			self.view_template.setToolTip("View the template that is being used (if more than one template is shown, it is the last one).")
			hbl_ww.addWidget(self.view_template)

			self.autobox=QtGui.QPushButton(QtGui.QIcon(get_image_directory() + "green_boxes.png"),"Autobox")
			self.autobox.setEnabled(False)
			self.autobox.setToolTip("Force Swarm to autobox the current image")
			hbl_ww.addWidget(self.autobox)

			vbl.addLayout(hbl_ww)

			hbl_aa = QtGui.QHBoxLayout()
			self.update_template = QtGui.QCheckBox("Refresh Template")
			self.update_template.setToolTip("Whether or not the act of adding a reference should force an update of the template being used by Swarm.\nOnce you have an adequate template you can turn this off and interactive picking will be faster.")
			self.update_template.setChecked(True)
			hbl_aa.addWidget(self.update_template)


			self.auto_update = QtGui.QCheckBox("Auto Update")
			self.auto_update.setToolTip("Whether or not autoboxing should occur every time you change a parameter or select a different image. This is the old dynapix button.")
			self.auto_update.setChecked(db.setdefault("auto_update",True))
			hbl_aa.addWidget(self.auto_update)
			vbl.addLayout(hbl_aa)


			self.advanced_hbl2 = QtGui.QHBoxLayout()
			self.enable_interactive_threshold  = QtGui.QCheckBox("Interactive Threshold")
			self.enable_interactive_threshold.setToolTip("Tweak the correlation threshold that is used to select particles.")
			self.enable_interactive_threshold.setChecked(False)
			from valslider import ValSlider
			self.thr = ValSlider(None,(0.0,3.0),"")
			self.thr.setValue(1.0)
			self.thr.setEnabled(False)
			self.advanced_hbl2.addWidget(self.enable_interactive_threshold)
			self.advanced_hbl2.addWidget(self.thr)
			vbl.addLayout(self.advanced_hbl2)

			self.overlap_hbl = QtGui.QHBoxLayout()
			self.enable_overlap_removal  = QtGui.QCheckBox("Proximity Threshold")
			self.enable_overlap_removal.setToolTip("Remove closely positioned particles.")
			self.enable_overlap_removal.setChecked(False)
			self.proximity_thr = ValSlider(None,(0,self.particle_diameter*2),"")
			self.proximity_thr.setValue(0.001)
			self.proximity_thr.setEnabled(False)
			self.overlap_hbl.addWidget(self.enable_overlap_removal)
			self.overlap_hbl.addWidget(self.proximity_thr)
			vbl.addLayout(self.overlap_hbl)

			hbl_bb = QtGui.QHBoxLayout()
			self.step_back=QtGui.QPushButton("Step Back")
			self.step_back.setToolTip("Recall the previous Swarm states")
			self.step_back.setEnabled(False)
			hbl_bb.addWidget(self.step_back)
			self.step_forward=QtGui.QPushButton("Step Forward")
			self.step_forward.setToolTip("Undo a step back")
			self.step_forward.setEnabled(False)
			hbl_bb.addWidget(self.step_forward)
			vbl.addLayout(hbl_bb)


			QtCore.QObject.connect(self.ptcl_diam_edit,QtCore.SIGNAL("editingFinished()"),self.new_ptcl_diam)
			QtCore.QObject.connect(self.update_template,QtCore.SIGNAL("clicked(bool)"),self.update_template_checked)
			QtCore.QObject.connect(self.auto_update,QtCore.SIGNAL("clicked(bool)"),self.auto_update_checked)
			QtCore.QObject.connect(self.clear, QtCore.SIGNAL("clicked(bool)"), self.clear_clicked)
			QtCore.QObject.connect(self.view_template, QtCore.SIGNAL("clicked(bool)"), self.view_template_clicked)
			QtCore.QObject.connect(self.autobox, QtCore.SIGNAL("clicked(bool)"), self.auto_box_clicked)
			QtCore.QObject.connect(self.method_group,QtCore.SIGNAL("buttonClicked (QAbstractButton *)"),self.method_group_clicked)
			QtCore.QObject.connect(self.enable_interactive_threshold, QtCore.SIGNAL("clicked(bool)"), self.interact_thresh_clicked)
			QtCore.QObject.connect(self.thr,QtCore.SIGNAL("sliderReleased"),self.new_threshold_release)
			QtCore.QObject.connect(self.thr,QtCore.SIGNAL("textChanged"),self.new_threshold_text_changed)
			QtCore.QObject.connect(self.step_back, QtCore.SIGNAL("clicked(bool)"), self.step_back_clicked)
			QtCore.QObject.connect(self.step_forward, QtCore.SIGNAL("clicked(bool)"), self.step_forward_clicked)
			QtCore.QObject.connect(self.proximity_thr,QtCore.SIGNAL("sliderReleased"),self.proximity_threshold_release)
			QtCore.QObject.connect(self.proximity_thr,QtCore.SIGNAL("textChanged"),self.proximity_threshold_text_changed)
			QtCore.QObject.connect(self.enable_overlap_removal, QtCore.SIGNAL("clicked(bool)"), self.enable_overlap_removal_clicked)
		return self.widget

	def update_states(self,swarm_boxer):
		self.busy = True
		self.ptcl_diam_edit.setText(str(swarm_boxer.particle_diameter))
		mode = swarm_boxer.pick_mode
		if mode == SwarmBoxer.THRESHOLD: self.thrbut.setChecked(True)
		elif mode == SwarmBoxer.SELECTIVE: self.selbut.setChecked(True)
		elif mode == SwarmBoxer.MORESELECTIVE: self.morselbut.setChecked(True)
		else: raise RuntimeError("This shouldn't happen")

		if swarm_boxer.proximity_threshold != None:
			self.enable_overlap_removal.setChecked(True)
			self.proximity_thr.setEnabled(True)
			self.proximity_thr.setValue(swarm_boxer.proximity_threshold)
		else:
			self.enable_overlap_removal.setChecked(False)
			self.proximity_thr.setEnabled(False)

		#self.auto_update.setChecked(swarm_boxer.auto_update)
		t = swarm_boxer.peak_score
		if swarm_boxer.interactive_threshold != None:
			t = swarm_boxer.interactive_threshold
			self.enable_interactive_threshold.setChecked(True)
			self.thr.setEnabled(True)
		else:
			self.enable_interactive_threshold.setChecked(False)
			self.thr.setEnabled(False)
		self.set_picking_data(t,swarm_boxer.profile,swarm_boxer.profile_trough_point)

		if swarm_boxer.templates != None and len(swarm_boxer.templates) > 0:
			self.enable_update_template(True)
			self.update_template.setChecked(swarm_boxer.update_template)
		else:
			self.set_update_template(True,False)

		self.enable_view_template(swarm_boxer.templates != None)

		#self.busy = True # BEWARE self.set_picking_data set it False!

		self.busy = False


	def enable_overlap_removal_clicked(self,val):
		if self.busy: return
		self.proximity_thr.setEnabled(val)
		self.target().set_proximity_removal_enabled(val)

	def proximity_threshold_text_changed(self):
		if self.busy: return
		val = self.proximity_thr.getValue()
		self.target().set_proximity_threshold(val)

	def proximal_threshold(self):
		'''
		Gets the value stored by the proximity threshold slider
		'''
		return self.proximity_thr.getValue()

	def proximity_threshold_release(self,val):
		if self.busy: return
		val = float(val)
		self.target().set_proximity_threshold(val)

	def new_threshold_text_changed(self):
		if self.busy: return
		val = self.thr.getValue()
		self.target().set_interactive_threshold(val)

	def new_threshold_release(self,val):
		if self.busy: return
		val = float(val)
		self.target().set_interactive_threshold(val)

	def interact_thresh_clicked(self,val):
		self.thr.setEnabled(val)
		if val == False:
			self.target().disable_interactive_threshold()

	def set_picking_data(self,threshold,profile=[],trough_point=None):
		self.busy = True
		if threshold == None: self.thr.setValue(0)
		else: self.thr.setValue(threshold)
		help = ""
		if profile != None and len(profile) > 0:
			help += "The Swarm profile is ["
			for i,p in enumerate(profile):
				if i != 0: help += " "
				help += "%.3f" %p

			help += "]\n"

		if trough_point:
			help += "The trough point is " + str(trough_point) + "."

		self.thr.setToolTip(help)
		self.busy = False

	def new_ptcl_diam(self):
		if self.busy: return
		self.target().set_particle_diameter(int(self.ptcl_diam_edit.text()))

	def update_template_checked(self,val):
		if self.busy: return
		self.target().set_update_template(val)

	def clear_clicked(self,val):
		self.target().clear_all()

	def method_group_clicked(self,button):
		if self.busy: return
		self.target().set_pick_mode(str(button.text()))

	def view_template_clicked(self,val):
		self.target().view_template_clicked()

	def auto_box_clicked(self,val):
		self.target().auto_box_clicked()

	def enable_view_template(self,val):
		self.view_template.setEnabled(val)

	def auto_update_checked(self,val):
		if self.busy: return
		db = js_open_dict(SwarmPanel.DB_NAME)
		db["auto_update"] = val
		self.target().set_auto_update(val)

	def set_update_template(self,val,enable=True):
		self.busy = True
		self.update_template.setChecked(val)
		self.update_template.setEnabled(enable)
		self.busy = False

	def enable_update_template(self,enable=True):
		self.busy = True
		self.update_template.setEnabled(enable)
		self.busy = False

	def enable_auto_box(self,val):
		self.autobox.setEnabled(val)

	def	step_back_clicked(self,val):
		self.target().step_back()

	def step_forward_clicked(self,val):
		self.target().step_forward()

	def enable_step_back(self,val,total=None):
		self.step_back.setEnabled(val)
		if total != None: self.step_back.setToolTip("%d backward steps available" %total)
		else: self.step_back.setToolTip("")

	def enable_step_forward(self,val,total=None):
		if self.widget == None: self.get_widget()
		self.step_forward.setEnabled(val)
		if total != None: self.step_forward.setToolTip("%d forward steps available" %total)
		else: self.step_forward.setToolTip("")


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
	CACHE_MAX = 10 # Each image has its last CACHE_MAX SwarmBoxer instance stored (or retrievable) automatically
	PROFILE_MAX = 0.8 # this is a percentage - it stops the profile trough point from going to the end
	REF_NAME = "swarm_ref"
	AUTO_NAME = "swarm_auto"
	WEAK_REF_NAME = "swarm_weak_ref"
	INIT = True
	MVT_THRESHOLD = 200 # a squared distance - used by the movement cache to automatically update using previously supplied user movement data
	SWARM_BOXERS = "swarm_boxers"
	SWARM_FWD_BOXERS = "swarm_fwd_boxers"
	SWARM_USER_MVT = "swarm_user_mvt"
	def __init__(self,particle_diameter=128):
		if SwarmBoxer.INIT: # this solved a strange problem with databases
			SwarmBoxer.INIT = False
			EMBox.set_box_color(SwarmBoxer.REF_NAME,[0,0,0])
			EMBox.set_box_color(SwarmBoxer.WEAK_REF_NAME,[0.2,0.2,0.4])
			EMBox.set_box_color(SwarmBoxer.AUTO_NAME,[0.4,.9,.4]) # Classic green, always was this one ;)

		self.panel_object = None # maybe it doesn't exist
		self.particle_diameter = particle_diameter
		self.update_template = True # Tied to the SwarmPanel - Whether or not the act of adding a reference should force an update of the template
		self.pick_mode = SwarmBoxer.SELECTIVE	# the autobox method - see EMData::BoxingTools for more details
		self.interactive_threshold = None  # Tied to the SwarmPanel -  this is a float when the user is playing with the threshold, but is None elsewise
		self.reset() # the first time this is called this establishes attribute variables - seeing as the function is required more than once it makes sense to do this
		self.auto_update = True # Tied to the SwarmPanel - this is the historical 'dynapix' button. It means that if any picking parameter is altered then autoboxing will be automatically triggered
		self.signal_template_update = False

		BoxingTools.set_mode(BoxingTools.CmpMode.SWARM_AVERAGE_RATIO) # this probably never needs to change - this mode has the best statistics

		self.template_viewer = None

		self.step_back_cache = [] # this will stored lists that store the parameters necessary to reconstruct a swarm boxer
		self.step_fwd_cache = [] # this will store the redo lists
		self.mvt_cache = [] # we have to remember if any of the auto selected boxes were moved, so if the user reboxes then the movements they previously supplied will be applied

		self.gui_mode = False # set this to False to stop any calls to Qt - such as the act of making the cursor busy...

	def __del__(self):
		'''
		Closes the template viewer if it exists
		'''
		if self.template_viewer != None:
			self.template_viewer.close()
			self.template_viewer = None

	def reset(self):
		'''
		Sets the self.ref_boxes, self.templates, self.profile, self.peak_score and self.profile_trough_point to safe defaults.
		Called in a couple of locations internally (notable, in the __init__ function)
		'''
		self.ref_boxes = []
		self.templates = None
		self.profile = None
		self.peak_score = None
		self.profile_trough_point = None
		self.proximity_threshold = None
		self.proximal_boxes = [] # this is like a temporary cache that stores boxes removed  by the application of the proximity threshold. If the user shortens the proximity threshold then previously removed boxes can be recalled.

	def to_list(self):
		'''
		Stores the vital attributes of this object in a list - returns a deep copy. The list contains these entries:
		-----------------------------
		Index   Variable			Description
		0		self.ref_boxes		List of SwarmBox object - this information could technically be used to regenerate the template and the picking parameters
		1		self.templates		List of EMData objects which are iteratively honed templates - the last is the one that was used for the correlation image
		2		self.profile		List of floats - The Swarm profile
		3		self.peak_score		Float - The correlation threshold being used to select particles
		4		self.profile_trough_point	Int - The trough point, used in the selective mode of Swarm picking
		5		self.particle_diameter		Int- User supplied estimate of the particle diameter
		6		self.update_template		Bool - Flag indicating whether adding a references should cause an update of the template
		7		self.pick_mode				Int - either SwarmBoxer.THRESHOLD, SwarmBoxer.SELECTIVE, or SwarmBoxer.MORESELECTIVE
		8		self.interactive_threshol	Float or None - If not None, this is a theshold value that overrides self.peak_score
		9		unique id					String - that which is return by gm_time_string - this can be used to decide what is the best autoboxer (i.e. is not really related to interactive boxing)
		10		self.proximity_threshold	Float or None - If not None, this is a threshold used to remove auto-selected particles that are too close to each other
		-----------------
		Note that it would be possible to append entries, but that if you disrupt the order which is specified above you will destroy back compatibility
		'''
		l = [self.ref_boxes,self.templates,self.profile,self.peak_score,self.profile_trough_point,self.particle_diameter,self.update_template]
		l.extend([self.pick_mode, self.interactive_threshold, gm_time_string(),self.proximity_threshold])
		return deepcopy(l)

	def load_from_list(self,l):
		'''
		Sets almost every variable attributes to a new value by assuming a specific ordering of the incoming list.
		The incoming list should almost certainly have been generated, at some point in time, but the call to self.to_list -
		See the comments in self.to_list for more details.
		@param l a list that was almost certainly generated, at some point in time, but a call to self.to_list
		'''
		self.ref_boxes = l[0]
		self.templates = l[1]
		self.profile = l[2]
		self.peak_score = l[3]
		self.profile_trough_point = l[4]
		self.particle_diameter = l[5]
		self.update_template = l[6]
		self.pick_mode = l[7]
		self.interactive_threshold = l[8]
		# entry 9 is not required, it is the creation stamp used to ascertain if common boxer instances are being used by different images, it has no relevance to interactive boxing
		try: # this try/except is for back compatibility only

			self.proximity_threshold = l[10]
		except:
			self.proximity_threshold = None

	def step_back(self):
		'''
		That which is called by the SwarmPanel when "Step Back" is clicked
		Manages all the various scenarios
		'''

		if len(self.step_back_cache) == 0: raise RuntimeError("Step backward cache has no entries")
		l = self.step_back_cache.pop(-1)
		self.step_fwd_cache.append(l)
		self.panel_object.enable_step_forward(True,len(self.step_fwd_cache))

		if len(self.step_back_cache) > 0:
			self.panel_object.enable_step_back(True,len(self.step_back_cache))
			self.load_from_last_state()
		else:
			self.target().clear_boxes([SwarmBoxer.REF_NAME,SwarmBoxer.AUTO_NAME,SwarmBoxer.WEAK_REF_NAME])
			self.reset()
			self.panel_object.enable_step_back(False)
			self.target().get_2d_window().updateGL()

		self.panel_object.update_states(self)
		if self.template_viewer != None:
			self.template_viewer.set_data(self.templates,soft_delete=True)
			self.template_viewer.updateGL()

		if self.templates != None:
			self.panel_object.enable_view_template(True)

#		self.debug_back("back")
#
#	def debug_back(self,s):
#		for i,c in enumerate(self.step_back_cache):
#			for j,b in enumerate(c[0]):
#				print s,i,j,b.in_template
#		print "---"

	def load_from_last_state(self):
		'''
		A function called in more than one location. Gets the last entry in the self.step_back_cache and harnesses
		the relevant parameters from it, then clears the EMBoxerModule of relevant boxes, and then an autobox is
		forced. Making the autoboxing sensitive to the self.auto_update parameter would make it messy and tricky to handle.
		'''
		if len(self.step_back_cache) == 0: raise RuntimeError("can't proceed there are no entries in the step_back cache")
		self.target().clear_boxes([SwarmBoxer.REF_NAME,SwarmBoxer.AUTO_NAME,SwarmBoxer.WEAK_REF_NAME])
#		self.debug_back("hb a")
		self.load_from_list(deepcopy(self.step_back_cache[-1]))
		boxes = [ [box.x, box.y, box.type_name(), box.peak_score] for box in self.ref_boxes if box.image_name == self.current_file ]
#		self.debug_back("hb b")
		self.target().add_boxes(boxes)
		self.target().get_2d_window().updateGL()
		if self.templates != None:
			self.auto_box(self.target().current_file(),parameter_update=False,cache=False)
		else:
			# it's an empty boxer... that's fine
			pass
			#print "handle the case of no templates ?"

	def step_forward(self):
		'''
		Manages the case when the user has clicked "Step Forward" in the SwarmPanel
		In this case the last entry in self.step_fwd_cache is used to set the states of
		the relevant variable attributes of this object. After this we have to clear the
		EMBoxerModule and redo autoboxing. Note that I made this function force an autobox,
		simply because the it would get messy otherwise.
		Handles any necessary updates of the SwarmPanel
		'''
		if len(self.step_fwd_cache) == 0: raise RuntimeError("Step forward cache has no entries")
		l = self.step_fwd_cache.pop(-1)

		if len(self.step_fwd_cache) == 0: self.panel_object.enable_step_forward(False)
		else: self.panel_object.enable_step_forward(True,len(self.step_fwd_cache)) # solely for the tooltip

		self.panel_object.enable_step_back(True,len(self.step_back_cache))
		self.step_back_cache.append(l)

		self.target().clear_boxes([SwarmBoxer.REF_NAME,SwarmBoxer.AUTO_NAME,SwarmBoxer.WEAK_REF_NAME])
		self.load_from_list(deepcopy(self.step_back_cache[-1]))
		boxes = [ [box.x, box.y, box.type_name(), box.peak_score] for box in self.ref_boxes if box.image_name == self.current_file]
		self.target().add_boxes(boxes)
		self.target().get_2d_window().updateGL()
		if len(self.ref_boxes) > 0:
			self.auto_box(self.target().current_file(),cache=False)

		self.panel_object.update_states(self)
		if self.template_viewer != None:
			self.template_viewer.set_data(self.templates,soft_delete=True)
			self.template_viewer.updateGL()

		if self.templates != None:
			self.panel_object.enable_view_template(True)

	def cache_current_state(self):
		'''
		Adds the current state to self.step_back_cache
		As a result resets self.step_fwd_cache (think in terms of undo and redo)
		Also updates the SwarmPanel
		'''
		self.step_back_cache.append(self.to_list())
		self.step_fwd_cache = []
		if self.gui_mode: self.panel_object.enable_step_forward(False)
		if self.gui_mode: self.panel_object.enable_step_back(True,len(self.step_back_cache))
		if len(self.step_back_cache) >= SwarmBoxer.CACHE_MAX: self.step_back_cache.pop(0)

	def cache_to_database(self):
		'''
		Stores the SwarmBoxer.SWARM_BOXERS and SwarmBoxer.SWARM_FWD_BOXER entries on the local database
		'''
		set_idd_image_entry(self.target().current_file(),SwarmBoxer.SWARM_BOXERS,self.step_back_cache)
		set_idd_image_entry(self.target().current_file(),SwarmBoxer.SWARM_FWD_BOXERS,self.step_fwd_cache)

	def handle_file_change(self,file_name,active_tool=False):
		'''
		Handles the situation when the user changes the current image being studied in the main interface
		in the EMBoxerModule.
		An autobox will occur if the image being changed to is in an empty state this the current image
		has Swarm Data that is useful
		Specifically, the last boxer from the previous image (i.e. self.step_back_cache[-1]) is pushed to the top
		of the boxer stack for the new image and auto boxing occurs. This is easily undone if the results
		are unwanted, i.e. by pressing "Step Back" in the SwarmPanel

		A variety of basic things happen too. And cached boxer data are restored from the local database,
		allowing the user to step back and step forward (as in the buttons in the Swarm Panel), etc.

		'''
		l = None
		#if active_tool:
			#if self.auto_update and len(self.step_back_cache) > 0:
		if len(self.step_back_cache) > 0:
			l = deepcopy(self.step_back_cache[-1])
			if l[1] == None:
				# there is no template, it's an empty autoboxer
				l = None

		self.reset()
		self.mvt_cache = get_database_entry(file_name,SwarmBoxer.SWARM_USER_MVT,dfl=[])

		self.step_fwd_cache = get_idd_image_entry(file_name,SwarmBoxer.SWARM_FWD_BOXERS,dfl=[])
		if len(self.step_fwd_cache) > 0: self.panel_object.enable_step_forward(True,len(self.step_fwd_cache))
		else:
			if self.panel_object: self.panel_object.enable_step_forward(False)

		self.step_back_cache = get_idd_image_entry(file_name,SwarmBoxer.SWARM_BOXERS,dfl=[])
		if self.step_back_cache == None: self.step_back_cache = []
		if len(self.step_back_cache) > 0 and l != None:
#			print self.step_back_cache[-1][9],l[9]
#			if self.step_back_cache[-1][9] == l[9]: # this is the time stamp - if they match we definitely shouldn't push on to the stacks
#				print "saved a redundant step"
#				l = None
			if self.step_back_cache[-1][1] != None: # this means the handing on of parameters only ever happens if the current state is empty
				l = None # there is a template, it's an non-trivial autoboxer
			else:
				pass
				# this means we're in the clear in terms of automatically boxing this image

		if l != None:
			self.step_back_cache.append(l)
			self.cache_to_database()

		if len(self.step_back_cache) > 0:
			if self.panel_object: self.panel_object.enable_step_back(True,len(self.step_back_cache))
		else:
			if self.panel_object: self.panel_object.enable_step_back(False)

		if l != None:
			self.load_from_last_state()
		else:
			if len(self.step_back_cache) > 0: self.load_from_list(deepcopy(self.step_back_cache[-1]))
			# else we're probably establishing an empty slate

		if self.panel_object: self.panel_object.update_states(self)
		if self.template_viewer != None:
			self.template_viewer.set_data(self.templates,soft_delete=True)
			self.template_viewer.updateGL()

	def get_proximal_boxes(self,boxes):
		'''
		very inefficient, but does the job
		'''
		if self.proximity_threshold == None: return [] # you probably should not have called this function if the proximity threshold is None anyway
		if len(boxes) < 2: return [] # can't remove overlapping in this case

		return_idxs = []

		nearness_sq = self.proximity_threshold**2 # avoid use of sqrt

		if isinstance(boxes,dict): keys = boxes.keys()
		elif isinstance(boxes,list): keys = [i for i in range(len(boxes))]

		for i,key in enumerate(keys):
			box1 = boxes[key]
			collision = False
			for j in range(i+1,len(keys)):
				key2 = keys[j]
				box2 = boxes[key2]
				if ((box1[0]-box2[0])**2 + (box1[1]-box2[1])**2) < nearness_sq:
					collision = True
					if box2[2] == SwarmBoxer.AUTO_NAME and return_idxs.count(key2) == 0: return_idxs.append(key2)

			if collision:
				# for efficiency
				if box1[2] == SwarmBoxer.AUTO_NAME and return_idxs.count(key) == 0: return_idxs.append(key)

		return return_idxs

	def set_proximity_removal_enabled(self,val):
		'''
		'''
		if val == False:
			self.proximity_threshold = None
			if len(self.proximal_boxes) != 0:
				self.target().add_boxes(self.proximal_boxes)
				self.proximal_boxes = []
		else:
			if self.proximity_threshold == None:
				self.proximity_threshold = self.panel_object.proximal_threshold()
				self.__remove_proximal_particles()
	def set_proximity_threshold(self,val):
		'''
		The SwarmPanel call this function when the user changes the proximity threshold
		'''
		if self.proximity_threshold == val: return



		if self.proximity_threshold == None or self.proximity_threshold < val or len(self.proximal_boxes) == 0:
			self.proximity_threshold = val
			self.__remove_proximal_particles()
		else:
			from PyQt4 import QtCore
			get_application().setOverrideCursor(QtCore.Qt.BusyCursor)
			self.proximity_threshold = val
			if len(self.proximal_boxes) == 0: return

			add_idxs = self.check_proximity_add_boxes(self.proximal_boxes)

			if len(add_idxs) > 0:
				add_idxs.sort()
				add_idxs.reverse()
				# already in reverse order
				boxes = []
				for idx in add_idxs:
					boxes.append(self.proximal_boxes.pop(idx))

				self.target().add_boxes(boxes)
			get_application().setOverrideCursor(QtCore.Qt.ArrowCursor)

		cache = get_idd_image_entry(self.target().current_file(),SwarmBoxer.SWARM_BOXERS,dfl=[])
		if len(cache) > 0:
			cache[-1][10] = self.proximity_threshold # just store the new proximity threshold
			get_idd_image_entry(self.target().current_file(),SwarmBoxer.SWARM_BOXERS,cache)



#	def __remove_proximal_particles(self):
#		from PyQt4 import QtCore
#		get_application().setOverrideCursor(QtCore.Qt.BusyCursor)
#		boxes = self.target().get_boxes_filt(SwarmBoxer.AUTO_NAME,as_dict=True)
#		proximal_boxes_idxs = self.get_proximal_boxes(boxes)
#		proximal_boxes_idxs.sort()
#		proximal_boxes_idxs.reverse()
#		proximal_boxes = [boxes[idx] for idx in proximal_boxes_idxs]
#
#		conv = [box.to_list() for box in proximal_boxes]
#		self.proximal_boxes.extend(conv)
#		self.target().remove_boxes(proximal_boxes_idxs)
#		get_application().setOverrideCursor(QtCore.Qt.ArrowCursor)
#
	def __remove_proximal_particles(self):
		if self.gui_mode:
			from PyQt4 import QtCore
			get_application().setOverrideCursor(QtCore.Qt.BusyCursor)
		boxes = self.target().get_boxes()
		proximal_boxes_idxs = self.get_proximal_boxes(boxes)
		proximal_boxes_idxs.sort()
		proximal_boxes_idxs.reverse()
		proximal_boxes = [boxes[idx] for idx in proximal_boxes_idxs]

		conv = [box.to_list() for box in proximal_boxes]
		self.proximal_boxes.extend(conv)
		self.target().remove_boxes(proximal_boxes_idxs)
		if self.gui_mode:
			from PyQt4 import QtCore
			get_application().setOverrideCursor(QtCore.Qt.ArrowCursor)

	def check_proximity_add_boxes(self,boxes):
		'''
		boxes is always a list
		'''

		keys = [i for i in range(len(boxes))]

		add_idxs = []
		nearness_sq = self.proximity_threshold**2 # avoid use of sqrt
		for i,key in enumerate(keys):
			box1 = boxes[key]
			nearness = None
			near_idx = None
			for j in range(0,len(keys)):
				if j == i: continue

				key2 = keys[j]
				box2 = boxes[key2]
				dist = ((box1[0]-box2[0])**2 + (box1[1]-box2[1])**2)
				if nearness == None or dist < nearness:
					nearness = dist
					near_idx = j

			if nearness > nearness_sq:
				add_idxs.append(i)

		return add_idxs

	def set_interactive_threshold(self,val):
		'''
		The SwarmPanel calls this function when the user changes the interactive threshold
		@param val the new interactive threshold to set
		'''
		if self.interactive_threshold == val: return

		self.interactive_threshold = val
		if self.auto_update and self.templates != None:
			self.auto_box(self.target().current_file(), parameter_update=True, force_remove_auto_boxes=True)
		elif not self.auto_update:
			self.panel_object.enable_auto_box(True)

	def disable_interactive_threshold(self):
		'''
		The SwarmPanel calls this function when the user clicks to disable the interactive threshold
		'''

		if self.interactive_threshold == None:
			# The user clicked "Interactive Threshold" and then clicked it again - i.e. nothing was altered
			return

		self.interactive_threshold = None
		if self.auto_update and self.templates != None:
			self.auto_box(self.target().current_file(), parameter_update=True, force_remove_auto_boxes=True)
		elif not self.auto_update:
			self.panel_object.enable_auto_box(True)

	def view_template_clicked(self):
		'''
		The SwarmPanel calls this function.
		Loads an EMImageMXWidget for viewing (if it doesn't already exist)
		Brings the template viewer to the foreground
		'''
		if self.template_viewer == None:
			from emimagemx import EMImageMXWidget
			self.template_viewer = EMImageMXWidget()

			self.template_viewer.set_data(self.templates,soft_delete=True) # should work if self.templates is None
			self.template_viewer.setWindowTitle("Templates")
			from PyQt4 import QtCore
			QtCore.QObject.connect(self.template_viewer,QtCore.SIGNAL("module_closed"),self.template_viewer_closed)

		get_application().show_specific(self.template_viewer)

	def template_viewer_closed(self):
		'''
		A slot for the module_closed signal from self.template_viewer
		'''
		self.template_viewer = None

	def set_update_template(self,val):
		'''
		Whether or not the act of adding a reference should for an update of the template
		@param val a boolena
		'''
		self.update_template = val

	def set_auto_update(self,val):
		'''
		Whether or not picked boxes should be updated when parameters are changed
		@param val a boolean
		'''
		self.auto_update = val

	def set_particle_diameter(self,diameter):
		'''
		set the particle diameter, which should be an int
		@param diameter the approximate diameter of the particle
		'''
		if diameter != self.particle_diameter:
			self.particle_diameter = diameter
			self.signal_template_update = True
			if self.auto_update and self.templates != None:
				self.auto_box(self.target().current_file())
			elif not self.auto_update: self.panel_object.enable_auto_box(True)

	def set_pick_mode(self,mode):
		'''
		Sets self.pick_mode
		@param mode should be in SwarmBoxer.THRESHOLD, SwarmBoxer.SELECTIVE, SwarmBoxer.MORESELECTIVE
		'''
		if mode not in [SwarmBoxer.THRESHOLD, SwarmBoxer.SELECTIVE, SwarmBoxer.MORESELECTIVE ]: raise RuntimeError("%s is an unknown SwarmBoxer mode" %mode)
		if mode != self.pick_mode:
			self.pick_mode = mode
			if self.auto_update and self.templates != None:
				self.auto_box(self.target().current_file(), parameter_update=False, force_remove_auto_boxes=True)
			elif not self.auto_update:
				self.panel_object.enable_auto_box(True)

	def move_auto_box(self,box_number,image_name,dx,dy):
		'''
		This is the solution to a conundrum - what if the user moves an auto box (a green one) but then reboxes, subsequently
		losing this "centering" metadata? The solution implemented here is to store the movement metadata and perform collision
		detection whenever autoboxing happens again.

		Here we update the "movement cache" to reflect the new metadata
		@param box_number the index corresponding to the stored box of the EMBoxList in the EMBoxerModule
		@param image_name the name of the image we're currently operating on - this is important because its used to store data in the local database
		@param dx the movement in the x direction (float)
		@param dy the movement in the y direction
		'''
		box = self.target().get_box(box_number)
		for i,[old,new] in enumerate(self.mvt_cache):
			dist = (new[0]-box.x)**2 + (new[1]-box.y)**2
			if dist < SwarmBoxer.MVT_THRESHOLD:
				self.mvt_cache[i][1] = [box.x+dx,box.y+dy]
				break
		else:
			self.mvt_cache.append([[box.x,box.y],[box.x+dx,box.y+dy]])

		set_database_entry(image_name,SwarmBoxer.SWARM_USER_MVT,self.mvt_cache)
		self.target().move_box(box_number,dx,dy)


	def update_auto_box_user_mvt(self,boxes):
		'''
		This is the solution to a conundrum - what if the user moves an auto box (a green one) but then reboxes, subsequently
		losing this "centering" metadata? The solution implemented here is to store the movement metadata and perform collision
		detection whenever autoboxing happens again.
		@param boxes a list of lists box data e.g. [[x,y,type,score],[x,y,type,score],....[int,int,str,float])
		'''
		for box in boxes:
			for i,[old,new] in enumerate(self.mvt_cache):
				dist = (old[0]-box[0])**2 + (old[1]-box[1])**2
				if dist < SwarmBoxer.MVT_THRESHOLD:
					box[0] = new[0]
					box[1] = new[1]

	def move_ref(self,box_number,image_name,dx,dy,set_moved=True,allow_template_update=True):
		'''
		Moves the reference stored internally and updates the EMBoxerModule so that it displays and stores the exact same thing.
		Separating what this object stores from what the EMBoxerModule stores decouple them and disentangles the overall design,
		but also makes some tasks 'fiddly'.

		This function takes care of telling the EMBoxerModule to correct what it's displaying - so the calling function should not.

		@param box_number the index corresponding to the stored box of the EMBoxList in the EMBoxerModule
		@param image_name the name of the image we're currently operating on
		@param dx the movement in the x direction (float)
		@param dy the movement in the y direction
		@param set_moved a boolean indicating whether the SwarmBox.ever_moved attribute should be set to True. Search for ever_moved to see where it's used
		'''
		box = self.target().get_box(box_number)
		for i,rbox in enumerate(self.ref_boxes):
			if rbox.x == box.x and rbox.y == box.y and os.path.basename(rbox.image_name) == os.path.basename(image_name):
				rbox.x += dx
				rbox.y += dy
				if set_moved:rbox.ever_moved = True
				new_box = EMBox(rbox.x,rbox.y,rbox.type_name(),rbox.peak_score)
				self.target().set_box(new_box,box_number,update_display=True)
				if box.type == SwarmBoxer.REF_NAME and allow_template_update:
					self.signal_template_update = True
				break
		else:
			raise RuntimeError("Attempt to move a reference failed")

	def remove_ref(self,box,image_name):
		'''
		Remove the reference, i.e. in response to a mouse delete event

		Does not tell the EMBoxModule to removed the displayed box - the calling function should do that.

		@param box an EMBox, as returned by the call to EMBoxerModule.get_box
		@param image_name the name of the image we're operating on
		'''
		for i,rbox in enumerate(self.ref_boxes):
			if rbox.x == box.x and rbox.y == box.y and os.path.basename(rbox.image_name) == os.path.basename(image_name):
				b = self.ref_boxes.pop(i)
				if self.auto_update:
					if len(self.ref_boxes) > 0:
						if box.type == SwarmBoxer.REF_NAME:
							self.signal_template_update=True
						self.auto_box(image_name,force_remove_auto_boxes=True)
					else:
						self.clear_all()
				else:
					self.panel_object.enable_auto_box(True)
				break
		else:
			raise RuntimeError("Attempt to remove a reference failed")


	def add_ref(self,x,y,image_name):
		'''
		Add a reference at the given coordinate and from the given image

		Does not tell the EMBoxerModule to add the box - the calling function must do this

		@param x the x coordinate of the box
		@param y the y coordinate of the box
		@parm image_name the name of the image we're operating on
		'''
		new_box = SwarmBox(x,y,image_name,self.update_template)
		self.ref_boxes.append(new_box)
		if self.update_template:
			self.signal_template_update = True
			box_num = self.target().add_box(x,y,type=SwarmBoxer.REF_NAME)
		else:
			box_num = self.target().add_box(x,y,type=SwarmBoxer.WEAK_REF_NAME)

		self.get_2d_window().updateGL()

		return box_num


	def ref_released(self,image_name,box_num):
		'''
		This function called when a reference box is released, i.e. by the mouse
		This can cause an autoboxing update if the self.auto_update parameter is true
		@param image_name the name of the image that is being operated on
		'''
		if self.auto_update:
			self.auto_box(image_name)
		else:
			self.try_to_center_ref(box_num)
			self.panel_object.enable_auto_box(True)

#
	def try_to_center_ref(self,box_num):
		if self.templates:
			shrink = self.get_subsample_rate()
			scaled_template = self.templates[-1].process("xform.scale",{"scale":shrink,"clip":self.particle_diameter})
			scaled_template.process_inplace("xform.centeracf")
			box = self.target().get_box(box_num)
			dx,dy = self.xform_center_propagate([box.x,box.y],self.target().current_file(),scaled_template,self.particle_diameter)
			self.move_ref(box_num,self.target().current_file(),dx,dy,set_moved=False)  # set moved = False so that it's also updated the next time through self.auto_box - this is probably unecessary but is not harmful
#

	def get_subsample_rate(self):
		'''
		Get the subsample rate advised by the SwarmBoxer, as based on self.particle_diameter and SWARM_TEMPLATE_MIN
		'''
		return int(math.ceil(float(self.particle_diameter)/float(SWARM_TEMPLATE_MIN)))

	def auto_box_clicked(self):
		'''
		When the autobox button is clicked then we force an autobox.
		'''
		self.auto_box(self.target().current_file(),force_remove_auto_boxes=True)

	def auto_box(self, image_name, parameter_update=True, force_remove_auto_boxes=False, cache=True):
		'''
		The main autoboxing function, this has a lot in it but it's not that complicated
		@param image_name the image that we're boxing
		@param parameter_update, should generally be True,but may be False if the current parameters are known to be current (see self.load_from_last_state)
		@param force_remove_auto_boxes if True all of the autoboxed boxes in the EMBoxerModule are removed and the 'entire' image is autoboxed again. This might be False if you know the template has not changed
		@param cache whether or not the newly establish state, i.e. at the end of this function, should be cached to the database and internally. Generally True but sometimes False (see self.load_from_last_state) .
		'''
		self.proximal_boxes = [] # this is always res
		if self.gui_mode:
			from PyQt4 import QtCore
			get_application().setOverrideCursor(QtCore.Qt.BusyCursor)

		if self.signal_template_update or force_remove_auto_boxes:
			self.target().clear_boxes([SwarmBoxer.AUTO_NAME])
			if self.gui_mode:
				self.get_2d_window().updateGL()


		if self.signal_template_update or not self.templates:
			if self.gui_mode:
				self.target().set_status_message("Updating Swarm Template",0)

			self.templates = gen_rot_ave_template(image_name,self.ref_boxes,self.get_subsample_rate(),self.particle_diameter)
			if self.gui_mode:
				self.panel_object.enable_update_template(True)
				self.target().set_status_message("Swarm Template Update Done",1000)
				self.panel_object.enable_view_template(True)
			self.signal_template_update = False
			if self.template_viewer != None:
				self.template_viewer.set_data(self.templates,soft_delete=True)
				self.template_viewer.updateGL()

		shrink = self.get_subsample_rate()

		exclusion_image = self.target().get_exclusion_image(mark_boxes=True)

		mediator = SwarmFLCFImageMediator(self.particle_diameter/(shrink*2), shrink, self.templates[-1])
		if self.gui_mode: self.target().set_status_message("Getting Correlation Image",0)
		correlation_image = FLCFImageCache.get_image(image_name,mediator)
		if self.gui_mode: self.target().set_status_message("Correlation Image Done",1000)

		exclusion_shrink = exclusion_image.get_attr("shrink")

		if shrink != exclusion_shrink:
			# the amount by which the exclusion is shrunken does not match the amount by which the SwarmBoxer shrinks - so we have to scale
			# to do: test this
			#print "shrink changed does this work?",shrink,exclusion_shrink,self.particle_diameter, SWARM_TEMPLATE_MIN,TEMPLATE_MIN
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


		if parameter_update:
			if not self.update_opt_picking_data():
				if self.gui_mode:
					from PyQt4 import QtCore
					get_application().setOverrideCursor(QtCore.Qt.ArrowCursor)
				print "funny error"
				return

		if self.pick_mode == SwarmBoxer.THRESHOLD: mode = 0
		elif self.pick_mode == SwarmBoxer.SELECTIVE: mode = 1
		elif self.pick_mode == SwarmBoxer.MORESELECTIVE: mode = 2

		searchradius = self.templates[-1].get_xsize()/2
		correlation = FLCFImageCache.get_image(image_name,mediator)
		# print "Correlation img is %s"%correlation

		if self.gui_mode: self.target().set_status_message("Autoboxing ....",0)

		# print "Using %s references"%(len(self.ref_boxes))

		# print "===========\n\nBoxingTools.auto_correlation_pick"
		# print "correlation: ", correlation
		# print "peak_score: ", self.peak_score
		# print "searchradius: ", searchradius
		# print "profile: ", self.profile
		# print "exclusion_image: ", exclusion_image
		# print "profile_trough_point: ", self.profile_trough_point
		# print "mode: ", mode
		# print "------"

		soln = BoxingTools.auto_correlation_pick(correlation, self.peak_score, searchradius, self.profile, exclusion_image, self.profile_trough_point, mode)

		# print "soln:", soln
		# print "==========="

		if self.gui_mode: self.target().set_status_message("Auboxing Done",1000)

		scaled_template = self.templates[-1].process("xform.scale",{"scale":shrink,"clip":self.particle_diameter})
		scaled_template.process_inplace("xform.centeracf")
		boxes = []
		for b in soln:
			x = b[0]
			y = b[1]
			xx = int(x*shrink)
			yy = int(y*shrink)
			type = SwarmBoxer.AUTO_NAME
			peak_score = correlation.get(x,y)
			box = [xx,yy,type,peak_score]
			self.center_propagate(box,image_name,scaled_template,self.particle_diameter)

			exc_x = box[0]/exclusion_shrink
			exc_y = box[1]/exclusion_shrink
			if exc_x >= exclusion_image.get_xsize() or exc_y > exclusion_image.get_ysize():
				print "Box position (%i,%i) was outside exclusion image boundaries (%i,%i)... ignoring (email this to woolford@bcm.edu)" %(exc_x,exc_y,exclusion_image.get_xsize(),exclusion_image.get_ysize())
				continue
			if exclusion_image.get(exc_x,exc_y) != 0: boxes.append(box)
			#else the particle was re-centered on an excluded region!

		if self.gui_mode: self.target().set_status_message("Updating Positions",0)
		self.update_auto_box_user_mvt(boxes)
		if self.gui_mode: self.target().set_status_message("Done",1000)

		for box in self.ref_boxes:
			if not box.ever_moved and box.image_name == self.target().current_file():
				box_num = self.target().detect_box_collision([box.x,box.y])
				if box_num == -1:
					raise RuntimeError("could not find reference")
				else:
					box = self.target().get_box(box_num)
					if box.type not in [SwarmBoxer.REF_NAME,SwarmBoxer.WEAK_REF_NAME]:
						raise RuntimeError("Did not get a reference when doing collision detection")

					dx,dy = self.xform_center_propagate([box.x,box.y],image_name,scaled_template,self.particle_diameter)
					self.move_ref(box_num,image_name,dx,dy,allow_template_update=False)


	   	boxes.sort(compare_box_correlation) # sorting like this will often put large ice contaminations in a group, thanks Pawel Penczek
		self.target().add_boxes(boxes, self.proximity_threshold == None)

		if self.proximity_threshold != None:
			self.__remove_proximal_particles()

		if cache:
			self.cache_current_state()
			self.cache_to_database()

		if self.gui_mode:
			from PyQt4 import QtCore
			self.panel_object.enable_auto_box(False)
			get_application().setOverrideCursor(QtCore.Qt.ArrowCursor)
			self.target().set_status_message("Autoboxed %d Particles" %len(boxes), 10000)
		else:
			print "Autoboxed %d Particles" %len(boxes)


		return boxes


	def update_opt_picking_data(self):
		'''
		This is the function that decides on the picking parameters of the SwarmBoxer, based on the reference
		boxes stored in self.ref_boxes
		'''
		self.profile = None
		self.peak_score = None
		self.profile_trough_point = None

		for box in self.ref_boxes:
			if self.interactive_threshold != None:
				self.peak_score = self.interactive_threshold
			else:
				if box.peak_score == None: continue

				if self.peak_score == None or box.peak_score < self.peak_score:
					self.peak_score = box.peak_score

			if self.profile == None: self.profile = deepcopy(box.profile)  # or else I'd alter it and muck everything up
			else:
				if len(self.profile) != len(box.profile): raise RuntimeError("This should never happen")

				profile = box.profile
				for j in xrange(0,len(self.profile)):
					if profile[j] < self.profile[j]: self.profile[j] = profile[j]

		if self.profile == None:
			return False

		max_radius =  int(len(self.profile)*SwarmBoxer.PROFILE_MAX)
		tmp = self.profile[0]
		self.profile_trough_point = 0
		for i in range(1,max_radius):
			if self.profile[i] > tmp and tmp > 0:
				tmp = self.profile[i]
				self.profile_trough_point = i

		if self.panel_object:
			self.panel_object.set_picking_data(self.peak_score, self.profile, self.profile_trough_point)

		return True

	def clear_all(self):
		'''
		Clears all associated boxes from the EMBoxerModule and internally, establishing a clean, blank state
		'''
		empty = (self.templates == None or len(self.templates) == 0)
		self.target().clear_boxes([SwarmBoxer.REF_NAME,SwarmBoxer.AUTO_NAME,SwarmBoxer.WEAK_REF_NAME],cache=True)
		self.reset()

		if self.panel_object:
			self.panel_object.enable_view_template(False)
			self.panel_object.enable_auto_box(False)

		if self.template_viewer != None:
			self.template_viewer.set_data(None)

		if not empty:
			self.cache_current_state()
			self.cache_to_database()
		# else:
		# 	print "avoided a redundant clear"

		self.update_template = True

		if self.panel_object:
			self.panel_object.set_update_template(True,False)


	def center_propagate(self,box,image_name,template,box_size):
		'''
		Centers the box argument in place
		@param box a list like [x,y,type,float] - only x and y are used
		@param image_name the name of the image we're operating on
		@param template the template correctly scaled to the have the same angstrom per pixel as the image (named image_name) stored on disk
		@param box_size the size of the box used to center - see xform_center_propagate
		'''
		dx,dy = self.xform_center_propagate(box,image_name,template,box_size)
		box[0] += dx
		box[1] += dy

#
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

	def boxes_erased(self,list_of_boxes,image_name):
		auto_box = False
		remove_happened = False
		for box in list_of_boxes:
			if box.type in [SwarmBoxer.REF_NAME,SwarmBoxer.WEAK_REF_NAME]:
				for i,rbox in enumerate(self.ref_boxes):
					if rbox.x == box.x and rbox.y == box.y and os.path.basename(rbox.image_name) == os.path.basename(image_name):
						remove_happened = True
						if self.auto_update: auto_box = True
						if box.type == SwarmBoxer.REF_NAME: self.signal_template_update=True
						b = self.ref_boxes.pop(i)
						break
		if auto_box:
			if len(self.ref_boxes) > 0:
				self.auto_box(image_name,force_remove_auto_boxes=True)
			else: self.clear_all()
		elif remove_happened:
			if len(self.ref_boxes) > 0:
				self.panel_object.enable_auto_box(True)
			else: self.clear_all()

class SwarmTool(SwarmBoxer,EMBoxingTool):
	'''
	A class that knows how to handle mouse erase events for a GUIBox
	'''

	def __init__(self,target,particle_diameter=128):
		SwarmBoxer.__init__(self,particle_diameter)
		self.target = weakref.ref(target)
		self.panel_object = SwarmPanel(self,self.particle_diameter)
		self.current_file = None # the name of the file that is being studied in the main viewer
		self.moving = None # keeps track of the moving box, will be a list in the format [[int,int],int,str] = [[x,y],box_number,box_type]
		self.ptcl_moving_data = None # keeps track of movement that's occuring in the particles (mx) viewer
		self.gui_mode = True

	def unique_name(self): return "Swarm"

	def get_widget(self):
		if self.panel_object == None:
			self.panel_object = SwarmPanel(self,self.particle_diameter)
		return self.panel_object.get_widget()

	def icon(self):
		from PyQt4 import QtGui
		return QtGui.QIcon(get_image_directory() + "swarm_icon.png")


	def set_current_file(self,file_name,active_tool=False):
		'''
		If the behavior of this Handler does not if the file changes, but the function needs to be supplied
		'''
		if self.current_file != file_name:
			self.current_file = file_name
			self.handle_file_change(file_name,active_tool)


	def get_2d_window(self): return self.target().get_2d_window()

	def mouse_move(self,event):
		pass

	def mouse_wheel(self,event):
		pass

	def mouse_down(self,event):
		m = self.get_2d_window().scr_to_img((event.x(),event.y()))
		box_num = self.target().detect_box_collision(m)
		from PyQt4.QtCore import Qt
		if box_num == -1:
			if event.modifiers()&Qt.ShiftModifier :
				return # the user tried to delete nothing
			box_num = self.add_ref(m[0],m[1],self.target().current_file())
			b = self.target().get_box(box_num)
			self.moving=[m,box_num,b.type]
		else:
			box = self.target().get_box(box_num)
			if box.type in [SwarmBoxer.REF_NAME,SwarmBoxer.AUTO_NAME,SwarmBoxer.WEAK_REF_NAME]:
				from PyQt4.QtCore import Qt
				if event.modifiers()&Qt.ShiftModifier :
					self.handle_box_delete(box,box_num)
				else:
					# if we make it here than the we're moving a box
					self.moving=[m,box_num,box.type]
#					self.target().moving_box_established(box_num)
			else:
				raise EMUnknownBoxType,box.type



	def handle_box_delete(self,box,box_num):
		if box.type == SwarmBoxer.AUTO_NAME:
			self.target().remove_box(box_num,exclude_region=True)
		elif box.type == SwarmBoxer.REF_NAME or box.type == SwarmBoxer.WEAK_REF_NAME:
			self.target().remove_box(box_num)
			self.remove_ref(box,self.target().current_file())
		else:
			raise EMUnknownBoxType,box.type

	def mouse_drag(self,event) :
		m=self.get_2d_window().scr_to_img((event.x(),event.y()))
		from PyQt4.QtCore import Qt
		if event.modifiers()&Qt.ShiftModifier:
			box_num = self.target().detect_box_collision(m)
			box = self.target().get_box(box_num)
			if ( box_num != -1):
				if box.type in [SwarmBoxer.REF_NAME,SwarmBoxer.AUTO_NAME,SwarmBoxer.WEAK_REF_NAME]:
					self.handle_box_delete(box,box_num)
		elif self.moving != None:
			oldm = self.moving[0]
			dx = m[0]-oldm[0]
			dy = m[1]-oldm[1]

			if self.moving[2] in [SwarmBoxer.REF_NAME,SwarmBoxer.WEAK_REF_NAME]:
				self.move_ref(self.moving[1],self.target().current_file(),dx,dy)
			else:
				self.move_auto_box(self.moving[1],self.target().current_file(),dx,dy)

			self.moving[0] = m

	def mouse_up(self,event) :
		if self.moving != None:
			self.target().box_released(self.moving[1])
			if self.moving[2] in [SwarmBoxer.REF_NAME,SwarmBoxer.WEAK_REF_NAME]:
				self.ref_released(self.target().current_file(),self.moving[1])

		self.moving= None

	def moving_ptcl_established(self,box_num,x,y):
		box = self.target().get_box(box_num)
		if box.type not in [SwarmBoxer.REF_NAME,SwarmBoxer.AUTO_NAME,SwarmBoxer.WEAK_REF_NAME]:
			raise EMUnknownBoxType,box.type

		self.ptcl_moving_data = [x,y,box_num]

	def move_ptcl(self,box_num,x,y,ptcls):
		if self.ptcl_moving_data == None: return

		dx = self.ptcl_moving_data[0] - x
		dy = y - self.ptcl_moving_data[1]
		box = self.target().get_box(box_num)

		if box.type in [SwarmBoxer.REF_NAME,SwarmBoxer.WEAK_REF_NAME]:
			self.move_ref(box_num,self.target().current_file(),dx,dy)
		else:
			self.move_auto_box(box_num,self.target().current_file(),dx,dy)

		self.ptcl_moving_data = [x,y,self.ptcl_moving_data[2]]

	def release_moving_ptcl(self,box_num,x,y):
		if self.ptcl_moving_data == None: return
		self.target().box_placement_update_exclusion_image_n(box_num)
		box = self.target().get_box(box_num)
		if box.type in [SwarmBoxer.REF_NAME,SwarmBoxer.WEAK_REF_NAME]:
			self.ref_released(self.target().current_file(),box_num)

		self.ptcl_moving_data = None

	def delete_ptcl(self,box_num):
		box = self.target().get_box(box_num)
		if box.type not in [SwarmBoxer.REF_NAME,SwarmBoxer.AUTO_NAME,SwarmBoxer.WEAK_REF_NAME]:
			raise EMUnknownBoxType,box.type

		self.handle_box_delete(self.target().get_box(box_num),box_num)

	def get_unique_box_types(self):
		return [SwarmBoxer.REF_NAME,SwarmBoxer.AUTO_NAME,SwarmBoxer.WEAK_REF_NAME]

	def boxes_erased(self,list_of_boxes):
		SwarmBoxer.boxes_erased(self,list_of_boxes,self.target().current_file())



############################################################################################################################################
# Add path for boxing using gauss convolution method (in pawelautoboxer used in sxboxer.py)

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

class GaussPanel:
	DB_NAME = "e2boxercache/gauss_panel.json"
	GDB_NAME = 'e2boxercache/gauss_box_DB.json' # cache for putting params related to gauss method autoboxer

	def __init__(self,target,particle_diameter=128):
		self.busy = True
		self.target = weakref.ref(target)
		self.widget = None
		self.busy = False
		self.ccfs = None
		self.data = None
		self.PRESIZE = 28
		self.THRNA='N/A'
		self.INVCONT=False
		self.USEVAR=True
		self.SLVAL=0
		self.GW = "1.0"
		self.ctf_inspector = None
		self.ctf_inspector_gone = True
		self.setgwbox =False

	def set_data( self, data ):
		self.ccfs = data
		#self.nbin = self.width()
		# hardcode nbin to 256 for now, which is the hardcoded width of the ccf histogram widget in sxboxer...
		self.nbin = 256
                self.data = histogram1d( data, self.nbin, self.PRESIZE )

		hmin = self.data[0][0]
		hmax = self.data[0][-1]

		#self.thr_low_edit.setText(str(hmin))
		#self.thr_hi_edit.setText(str(hmax))
		#self.new_thr_low()
		#self.new_thr_hi()

	def get_widget(self):
		if self.widget == None:

			gbdb = js_open_dict(GaussPanel.GDB_NAME)

			from PyQt4 import QtCore, QtGui, Qt
			self.widget = QtGui.QWidget()
			vbl = QtGui.QVBoxLayout(self.widget)
			vbl.setMargin(0)
			vbl.setSpacing(6)
			vbl.setObjectName("vbl")

			hgc = QtGui.QHBoxLayout()
			gconvheader = QtGui.QLabel("<b>Parameters of Gauss convolution</b>")
			hgc.addWidget(gconvheader)
			vbl.addLayout(hgc)

			hbl = QtGui.QHBoxLayout()
			pixel_input_label = QtGui.QLabel("Input Pixel Size:")
			pixel_input_label.setToolTip("Input pixel size")
			hbl.addWidget(pixel_input_label)

			pixin = gbdb.setdefault('pixel_input',None)
			if pixin == None:
				self.pixel_input_edit = QtGui.QLineEdit('1.0')
			else:
				self.pixel_input_edit = QtGui.QLineEdit(str(gbdb['pixel_input']))
			hbl.addWidget(self.pixel_input_edit)

			#vbl.addLayout(hbl)

			pixel_output_label = QtGui.QLabel("Output Pixel Size:")
			pixel_output_label.setToolTip("Output pixel size")
			hbl.addWidget(pixel_output_label)

			pixout_cache = gbdb.setdefault('pixel_output',None)
			if pixout_cache == None:
				self.pixel_output_edit = QtGui.QLineEdit('1.0')
			else:
				self.pixel_output_edit = QtGui.QLineEdit(str(gbdb['pixel_output']))
			hbl.addWidget(self.pixel_output_edit)
			self.new_pixel_output()
			self.new_pixel_input()
			vbl.addLayout(hbl)

			hbl_invcont = QtGui.QHBoxLayout()
			self.invert_contrast_chk = QtGui.QCheckBox("Invert Contrast")
			self.invert_contrast_chk.setToolTip("Invert contrast")
			invert_cache = gbdb.setdefault('invert_contrast',None)
			if invert_cache == None:
				self.invert_contrast_chk.setChecked(self.INVCONT)
				self.invert_contrast_checked(self.INVCONT)
			else:
				self.invert_contrast_chk.setChecked(invert_cache)
				self.invert_contrast_checked(invert_cache)
			hbl_invcont.addWidget(self.invert_contrast_chk)

			self.use_variance_chk = QtGui.QCheckBox("Use Variance")
			self.use_variance_chk.setToolTip("Use the variance image")
			use_variance_cache = gbdb.setdefault('use_variance',None)
			if use_variance_cache == None:
				self.use_variance_chk.setChecked(self.USEVAR)
				self.use_variance_checked(self.USEVAR)
			else:
				self.use_variance_chk.setChecked(use_variance_cache)
				self.use_variance_checked(use_variance_cache)
			hbl_invcont.addWidget(self.use_variance_chk)

			vbl.addLayout(hbl_invcont)

			hbl_gwidth = QtGui.QHBoxLayout()
			self.gauss_width_slider = QtGui.QSlider(QtCore.Qt.Horizontal)
			self.gauss_width_slider.setRange( -100, 100 )
			self.gauss_width_slider.setValue( self.SLVAL )
			hbl_gwidth.addWidget( self.gauss_width_slider)
			hbl_gwidth.addWidget(QtGui.QLabel("Gauss Width Adjust:"))
			self.gauss_width = QtGui.QLineEdit(self.GW)
			gauss_width_cache = gbdb.setdefault('gauss_width',None)
			if not(gauss_width_cache == None):
				self.gauss_width = QtGui.QLineEdit(str(gauss_width_cache))
			hbl_gwidth.addWidget( self.gauss_width)


			self.gauss_width_changed(self.SLVAL)
			self.gauss_width_edited()
			vbl.addLayout(hbl_gwidth)

			hbl_thr = QtGui.QHBoxLayout()
			thr_low_label = QtGui.QLabel("Threshold Low:")
			hbl_thr.addWidget(thr_low_label)
			self.thr_low_edit = QtGui.QLineEdit(self.THRNA)
			thrlow_cache = gbdb.setdefault('thr_low',None)
			if not(thrlow_cache == None):
				self.thr_low_edit = QtGui.QLineEdit(str(thrlow_cache))
			self.new_thr_low()
			hbl_thr.addWidget(self.thr_low_edit)
			thr_hi_label = QtGui.QLabel("Threshold High:")
			hbl_thr.addWidget(thr_hi_label)

			thrhi_cache = gbdb.setdefault('thr_hi',None)
			if thrhi_cache == None:
				self.thr_hi_edit = QtGui.QLineEdit(self.THRNA)
			else:
				self.thr_hi_edit = QtGui.QLineEdit(str(thrhi_cache))
			self.new_thr_hi()
			hbl_thr.addWidget(self.thr_hi_edit)
			vbl.addLayout(hbl_thr)

			hbl_ww = QtGui.QHBoxLayout()
			self.clear=QtGui.QPushButton("Clear Boxes")
			self.clear.setToolTip("Clear boxes generated by the Gauss mode.")
			hbl_ww.addWidget(self.clear)

			self.autobox=QtGui.QPushButton(QtGui.QIcon(get_image_directory() + "green_boxes.png"),"Run")
			self.autobox.setEnabled(True)
			self.autobox.setToolTip("Autobox using Gauss method")
			hbl_ww.addWidget(self.autobox)

			vbl.addLayout(hbl_ww)

			# add input fields for CTF estimation
			hgctf = QtGui.QHBoxLayout()
			ctftitle = QtGui.QLabel("<b>Parameters of CTF estimation</b>")
			hgctf.addWidget(ctftitle)
			vbl.addLayout(hgctf)

			hbl_wscs = QtGui.QHBoxLayout()

			window_size_label = QtGui.QLabel("Window size:")
			hbl_wscs.addWidget(window_size_label)
			self.ctf_window_size = QtGui.QLineEdit(str(gbdb.setdefault('ctf_window',"512")))

			hbl_wscs.addWidget(self.ctf_window_size)

			cs_label = QtGui.QLabel("Cs:")
			hbl_wscs.addWidget(cs_label)
			self.ctf_cs = QtGui.QLineEdit(str(gbdb.setdefault('ctf_cs',"2.0")))

			hbl_wscs.addWidget(self.ctf_cs)

			vbl.addLayout(hbl_wscs)

			hbl_esv = QtGui.QHBoxLayout()

			edge_size_label = QtGui.QLabel("Edge size:")
			hbl_esv.addWidget(edge_size_label)
			self.ctf_edge_size = QtGui.QLineEdit(str(gbdb.setdefault('ctf_edge',"0")))
			hbl_esv.addWidget(self.ctf_edge_size)

			voltage_label = QtGui.QLabel("Voltage:")
			hbl_esv.addWidget(voltage_label)
			self.ctf_volt = QtGui.QLineEdit(str(gbdb.setdefault('ctf_volt',"200.0")))
			hbl_esv.addWidget(self.ctf_volt)

			vbl.addLayout(hbl_esv)

			hbl_oac = QtGui.QHBoxLayout()

			overlap_label = QtGui.QLabel("Overlap:")
			hbl_oac.addWidget(overlap_label)
			self.ctf_overlap_size = QtGui.QLineEdit(str(gbdb.setdefault('ctf_overlap',"50")))
			hbl_oac.addWidget(self.ctf_overlap_size)

			amplitude_contrast_label = QtGui.QLabel("Amplitude Contrast:")
			hbl_oac.addWidget(amplitude_contrast_label)
			self.ctf_ampcont = QtGui.QLineEdit(str(gbdb.setdefault('ctf_ampcont',"10.0")))
			hbl_oac.addWidget(self.ctf_ampcont)

			vbl.addLayout(hbl_oac)

			# cter kboot
			hbl_kboot = QtGui.QHBoxLayout()
			kboot_label = QtGui.QLabel("kboot:")
			hbl_kboot.addWidget(kboot_label)

			self.ctf_kboot = QtGui.QLineEdit(str(gbdb.setdefault('ctf_kboot',"16")))
			hbl_kboot.addWidget(self.ctf_kboot)
			vbl.addLayout(hbl_kboot)

			hbl_estdef = QtGui.QHBoxLayout()
			hbl_fed = QtGui.QHBoxLayout()

			fstart_label = QtGui.QLabel("F_start:")
			hbl_fed.addWidget(fstart_label)
			self.ctf_f_start = QtGui.QLineEdit(str(gbdb.setdefault('ctf_fstart',"0.020")))
			hbl_fed.addWidget(self.ctf_f_start)

			estimated_defocus_label = QtGui.QLabel("Estimated defocus:")
			hbl_estdef.addWidget(estimated_defocus_label)
			self.estdef = QtGui.QLineEdit('')
			hbl_estdef.addWidget(self.estdef)
			vbl.addLayout(hbl_estdef)

			hbl_estdeferr = QtGui.QHBoxLayout()
			deferr_label = QtGui.QLabel("Estimated defocus error:")
			hbl_estdeferr.addWidget(deferr_label)
			self.deferr = QtGui.QLineEdit('')
			hbl_estdeferr.addWidget(self.deferr)
			vbl.addLayout(hbl_estdeferr)

			hbl_astamp = QtGui.QHBoxLayout()
			astig_amp_label = QtGui.QLabel("Estimated astigmatism \namplitude:")
			hbl_astamp.addWidget(astig_amp_label)
			self.astamp = QtGui.QLineEdit('')
			hbl_astamp.addWidget(self.astamp)
			vbl.addLayout(hbl_astamp)

			hbl_astamperr = QtGui.QHBoxLayout()
			astamperr_label = QtGui.QLabel("Estimated astigmatism \namplitude error:")
			hbl_astamperr.addWidget(astamperr_label)
			self.astamperr = QtGui.QLineEdit('')
			hbl_astamperr.addWidget(self.astamperr)
			vbl.addLayout(hbl_astamperr)

			hbl_astagl = QtGui.QHBoxLayout()
			astig_angle_label = QtGui.QLabel("Estimated astigmatism \nangle")
			hbl_astagl.addWidget(astig_angle_label)
			self.astagl = QtGui.QLineEdit('')
			hbl_astagl.addWidget(self.astagl)
			vbl.addLayout(hbl_astagl)

			hbl_astaglerr = QtGui.QHBoxLayout()
			astaglerr_label = QtGui.QLabel("Estimated astigmatism \nangle error:")
			hbl_astaglerr.addWidget(astaglerr_label)
			self.astaglerr = QtGui.QLineEdit('')
			hbl_astaglerr.addWidget(self.astaglerr)
			vbl.addLayout(hbl_astaglerr)

			hbl_ctf_cter = QtGui.QHBoxLayout()
			self.estimate_ctf_cter =QtGui.QPushButton("Estimate CTF using CTER")
			hbl_ctf_cter.addWidget(self.estimate_ctf_cter)
			vbl.addLayout(hbl_ctf_cter)



			#hbl_fed.addWidget(self.estdef)

			#vbl.addLayout(hbl_fed)

			#hbl_fs = QtGui.QHBoxLayout()

			#fstop_label = QtGui.QLabel("F_stop:")
			#hbl_fs.addWidget(fstop_label)
			#self.ctf_f_stop = QtGui.QLineEdit(str(gbdb.setdefault('ctf_fstop',"0.500")))
			#hbl_fs.addWidget(self.ctf_f_stop)

			#vbl.addLayout(hbl_fs)
			#hbl_ctf = QtGui.QHBoxLayout()
			#self.estimate_ctf=QtGui.QPushButton("Estimate CTF")
			#hbl_ctf.addWidget(self.estimate_ctf)

			#self.inspect_button=QtGui.QPushButton("Inspect CTF")
			#hbl_ctf.addWidget(self.inspect_button)

			#vbl.addLayout(hbl_ctf)


			QtCore.QObject.connect(self.pixel_input_edit,QtCore.SIGNAL("editingFinished()"),self.new_pixel_input)
			QtCore.QObject.connect(self.pixel_output_edit,QtCore.SIGNAL("editingFinished()"),self.new_pixel_output)
			QtCore.QObject.connect(self.autobox, QtCore.SIGNAL("clicked(bool)"), self.auto_box_clicked)
			QtCore.QObject.connect(self.clear, QtCore.SIGNAL("clicked(bool)"), self.clear_clicked)
			QtCore.QObject.connect(self.invert_contrast_chk,QtCore.SIGNAL("clicked(bool)"),self.invert_contrast_checked)
			QtCore.QObject.connect(self.use_variance_chk,QtCore.SIGNAL("clicked(bool)"),self.use_variance_checked)
			QtCore.QObject.connect(self.gauss_width_slider, QtCore.SIGNAL("valueChanged(int)"), self.gauss_width_changed)
			QtCore.QObject.connect(self.gauss_width, QtCore.SIGNAL("editingFinished()"), self.gauss_width_edited)
			QtCore.QObject.connect(self.thr_low_edit,QtCore.SIGNAL("editingFinished()"),self.new_thr_low)
			QtCore.QObject.connect(self.thr_hi_edit,QtCore.SIGNAL("editingFinished()"),self.new_thr_hi)
#			QtCore.QObject.connect(self.estimate_ctf,QtCore.SIGNAL("clicked(bool)"), self.calc_ctf)
#			QtCore.QObject.connect(self.inspect_button,QtCore.SIGNAL("clicked(bool)"), self.inspect_ctf)
			QtCore.QObject.connect(self.ctf_window_size,QtCore.SIGNAL("editingFinished()"),self.new_ctf_window)
			QtCore.QObject.connect(self.ctf_cs,QtCore.SIGNAL("editingFinished()"),self.new_ctf_cs)
			QtCore.QObject.connect(self.ctf_edge_size,QtCore.SIGNAL("editingFinished()"),self.new_ctf_edge)
			QtCore.QObject.connect(self.ctf_volt,QtCore.SIGNAL("editingFinished()"),self.new_ctf_volt)
			QtCore.QObject.connect(self.ctf_overlap_size,QtCore.SIGNAL("editingFinished()"),self.new_ctf_overlap_size)
			QtCore.QObject.connect(self.ctf_ampcont,QtCore.SIGNAL("editingFinished()"),self.new_ctf_ampcont)
			QtCore.QObject.connect(self.ctf_kboot,QtCore.SIGNAL("editingFinished()"),self.new_ctf_kboot)

			QtCore.QObject.connect(self.estimate_ctf_cter,QtCore.SIGNAL("clicked(bool)"), self.calc_ctf_cter)

		return self.widget

	def gauss_width_changed(self, v):
		from math import pow
		s = "%.3f" % pow(10.0, v*0.01)
		if self.setgwbox:
		        self.gauss_width.setText( s )
		        self.target().set_gauss_width(float(self.gauss_width.text()))

	def gauss_width_edited(self):
		from string import atof
		from math import log10
		text = self.gauss_width.text()
		v = int( log10(atof(text)) * 100)
		self.setgwbox = False
		self.gauss_width_slider.setValue( v )
		self.setgwbox = True
		self.target().set_gauss_width(float(text))

		gbdb = js_open_dict(GaussPanel.GDB_NAME)
		gbdb['gauss_width']=float(text)

	def update_states(self,gauss_boxer):
		self.busy = True
		self.pixel_input_edit.setText(str(gauss_boxer.pixel_input))

		self.busy = False

	def new_pixel_input(self):
		if self.busy: return
		pixin = float(self.pixel_input_edit.text())
		pixout=float(self.pixel_output_edit.text())
		if pixout < pixin:
			print "output pixel size cannot be smaller than input pixel size"
			self.pixel_input_edit.setText(str(self.target().pixel_input))
			return
		self.target().set_pixel_input(pixin)
		gbdb = js_open_dict(GaussPanel.GDB_NAME)
		gbdb['pixel_input']=pixin

	def new_pixel_output(self):
		if self.busy: return
		pixout=float(self.pixel_output_edit.text())
		pixin = float(self.pixel_input_edit.text())
		if pixout < pixin:
			self.pixel_output_edit.setText(str(self.target().pixel_output))
			print "output pixel size cannot be smaller than input pixel size"
			return
		self.target().set_pixel_output(pixout)
		gbdb = js_open_dict(GaussPanel.GDB_NAME)
		gbdb['pixel_output']=pixout

	def new_thr_low(self):
		if self.busy: return
		thrlow=self.thr_low_edit.text()
		if thrlow != self.THRNA:
			thrlow= float(thrlow)
		else:
			thrlow = None
		self.target().set_thr_low(thrlow)
		gbdb = js_open_dict(GaussPanel.GDB_NAME)
		gbdb['thr_low']=thrlow


	def new_thr_hi(self):
		if self.busy: return
		thrhi=self.thr_hi_edit.text()
		if thrhi != self.THRNA:
			thrhi=float(thrhi)
		else:
			thrhi=None
		self.target().set_thr_hi(thrhi)
		gbdb = js_open_dict(GaussPanel.GDB_NAME)
		gbdb['thr_hi']=thrhi

	def auto_box_clicked(self,val):
		self.target().auto_box_clicked()

	def clear_clicked(self,val):
		self.target().clear_all()

	def enable_auto_box(self,val):
		self.autobox.setEnabled(val)

	def invert_contrast_checked(self,val):
		if self.busy: return
		gbdb = js_open_dict(GaussPanel.GDB_NAME)
		gbdb["invert_contrast"] = val
		self.target().set_invert_contrast(val)

	def use_variance_checked(self,val):
		if self.busy: return
		gbdb = js_open_dict(GaussPanel.GDB_NAME)
		gbdb["use_variance"] = val
		self.target().set_use_variance(val)


	def new_ctf_window(self):
		if self.busy: return
		winsize=self.ctf_window_size.text()
		gbdb = js_open_dict(GaussPanel.GDB_NAME)
		gbdb['ctf_window']=int(winsize)

	def new_ctf_cs(self):
		if self.busy: return
		cs=self.ctf_cs.text()
		gbdb = js_open_dict(GaussPanel.GDB_NAME)
		gbdb['ctf_cs']=float(cs)

	def new_ctf_edge(self):
		if self.busy: return
		edge=self.ctf_edge_size.text()
		gbdb = js_open_dict(GaussPanel.GDB_NAME)
		gbdb['ctf_edge']=int(edge)

	def new_ctf_volt(self):
		if self.busy: return
		volt=self.ctf_volt.text()
		gbdb = js_open_dict(GaussPanel.GDB_NAME)
		gbdb['ctf_volt']=float(volt)

	def new_ctf_kboot(self):
		if self.busy: return
		kboot=self.ctf_kboot.text()
		gbdb = db_open_dict(GaussPanel.GDB_NAME)
		gbdb['ctf_kboot']=float(kboot)

	def new_ctf_overlap_size(self):
		if self.busy: return
		ov=self.ctf_overlap_size.text()
		gbdb = js_open_dict(GaussPanel.GDB_NAME)
		gbdb['ctf_overlap']=int(ov)

	def new_ctf_ampcont(self):
		if self.busy: return
		ac=self.ctf_ampcont.text()
		gbdb = js_open_dict(GaussPanel.GDB_NAME)
		gbdb['ctf_ampcont']=float(ac)

	def new_ctf_f_start(self):
		if self.busy: return
		fstart=self.ctf_f_start.text()
		gbdb = js_open_dict(GaussPanel.GDB_NAME)
		gbdb['ctf_fstart']=float(fstart)

	def new_ctf_f_stop(self):
		if self.busy: return
		fstop=self.ctf_f_stop.text()
		gbdb = js_open_dict(GaussPanel.GDB_NAME)
		gbdb['ctf_fstop']=float(fstop)

	# sxboxer's calc_ctf (from class EMBoxerModulePanel)
	# ctf is always calculated from original input micrograph
	def calc_ctf(self):
		# calculate power spectrum of image with welch method (welch_pw2)
		# calculate rotational average of power spectrum (rot_avg_table)
		# calculate ctf values with ctf_get
		#print "starting CTF estimation"
		# get the current image
		from utilities import get_im
		#image_name = self.target().boxable.get_image_name()
		#img = BigImageCache.get_image_directly( image_name )
		image_name = self.target().target().file_names[0]
		img = get_im(image_name)

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
		from fundamentals import rot_avg_table
		avg_sp = rot_avg_table(power_sp)
		del power_sp

		# print "determine ctf"
		from morphology import defocus_gett


		input_pixel_size = float(self.pixel_input_edit.text())
		output_pixel_size = float(self.pixel_output_edit.text())
		#print "Input pixel size: ", input_pixel_size
		#print "Output pixel size: ", output_pixel_size

		# this is wrong from sxboxer. wgh should be amplitude contrast
		#defocus = defocus_gett(avg_sp, voltage=ctf_volt, Pixel_size=input_pixel_size, Cs=ctf_cs, wgh=ctf_cs,f_start=ctf_f_start, f_stop=ctf_f_stop, parent=self)
		defocus = defocus_gett(avg_sp, voltage=ctf_volt, Pixel_size=input_pixel_size, Cs=ctf_cs, wgh=ctf_ampcont,f_start=ctf_f_start, f_stop=ctf_f_stop, parent=self)
		self.estdef.setText(str(defocus/10000.0))
		self.estdef.setEnabled(False)


		# update ctf inspector values

		if (self.ctf_inspector is not None):
			self.ctf_inspector.setData(self.ctf_data)
			self.ctf_inspector.i_start = self.i_start
			self.ctf_inspector.i_stop = self.i_stop
			if not(self.ctf_inspector_gone):
				self.ctf_inspector.update()
		else:
			global i_start_initial
			global i_stop_initial
			i_start_initial = self.i_start
			i_stop_initial = self.i_stop

		# XXX: wgh?? amp_cont static to 0?
		# set image properties, in order to save ctf values
		from utilities import set_ctf
		set_ctf(img, [defocus, ctf_cs, ctf_volt, input_pixel_size, 0, ctf_ampcont])
		# and rewrite image
		img.write_image(image_name)
		print [defocus, ctf_cs, ctf_volt, input_pixel_size, 0, ctf_ampcont]
		# get alternate, and set its ctf
		altimg=BigImageCache.get_object(image_name).get_image(use_alternate=True)
		set_ctf(altimg, [defocus, ctf_cs, ctf_volt, output_pixel_size, 0, ctf_ampcont])
		BigImageCache.get_object(image_name).register_alternate(altimg)
		print [defocus, ctf_cs, ctf_volt, output_pixel_size, 0, ctf_ampcont]
 		print "CTF estimation done."
 		#print "Estimated defocus value: ", defocus

		##############################################################################
		#### save ctf estimation parameters to db for command line batch processing
		gbdb = js_open_dict(GaussPanel.GDB_NAME)
		ctfdict = {'pixel_input':input_pixel_size,'pixel_output':output_pixel_size,'ctf_fstart':ctf_f_start,'ctf_fstop':ctf_f_stop, 'ctf_window':ctf_window_size,'ctf_edge':ctf_edge_size,'ctf_overlap':ctf_overlap_size,'ctf_ampcont':ctf_ampcont,'ctf_volt':ctf_volt,'ctf_cs':ctf_cs}
		#print "calc_ctf image_name: ", image_name
		if gbdb.has_key(image_name):
			olddict=gbdb[image_name]
			gbdb[image_name] = dict((olddict).items() + ctfdict.items() ) # merge the two dictionaries with conflict resolution resolved in favorr of the latest ctf parameters
		else:
			gbdb[image_name]=ctfdict

		del img
		del altimg

	def inspect_ctf(self):
		#display(self.ctf_data)

		if not(self.ctf_inspector):
			self.ctf_inspector = CTFInspectorWidget(self,self.ctf_data)
			self.ctf_inspector.show()
			self.ctf_inspector_gone=False
		else:
			if (self.ctf_inspector_gone):
				self.ctf_inspector.show()
				self.ctf_inspector_gone=False
			else:
				pass

	def calc_ctf_cter(self):
		# calculate ctf of ORIGINAL micrograph using cter in gui mode
		# this must mean cter is being calculated on a single micrograph!

		print "Starting CTER"
		# get the current image
		from utilities import get_im
		#image_name = self.target().boxable.get_image_name()
		#img = BigImageCache.get_image_directly( image_name )
		image_name = self.target().target().file_names[0]
		img = get_im(image_name)

		# conversion from text necessary
		try:
			ctf_window_size  = int(self.ctf_window_size.text())
			input_pixel_size = float(self.pixel_input_edit.text())
			output_pixel_size = float(self.pixel_output_edit.text())
			ctf_edge_size    = int(self.ctf_edge_size.text())
			ctf_overlap_size = int(self.ctf_overlap_size.text())
			ctf_volt         = float(self.ctf_volt.text())
			ctf_cs           = float(self.ctf_cs.text())
			ctf_ampcont      = float(self.ctf_ampcont.text())
			ctf_kboot        = int(self.ctf_kboot.text())

		except ValueError,extras:
			# conversion of a value failed!
			print "integer conversion failed."
			if not(extras.args is None):
				print extras.args[0]
			return
		except:
			print "error"
			return

		fname, fext = os.path.splitext(image_name)
		outpwrot = 'pwrot_%s'%fname
		outpartres = 'partres_%s'%fname

		if os.path.exists(outpwrot) or os.path.exists(outpartres):
			print "Please remove or rename %s and or %s"%(outpwrot,outpartres)
			return

		from morphology import cter
		defocus, ast_amp, ast_agl, error_defocus, error_astamp, error_astagl = cter(None, outpwrot, outpartres, None, None, ctf_window_size, voltage=ctf_volt, Pixel_size=input_pixel_size, Cs = ctf_cs, wgh=ctf_ampcont, kboot=ctf_kboot, MPI=False, DEBug= False, overlap_x = ctf_overlap_size, overlap_y = ctf_overlap_size, edge_x = ctf_edge_size, edge_y = ctf_edge_size, guimic=image_name)

		self.estdef.setText(str(defocus))
		self.estdef.setEnabled(True)

		self.astamp.setText(str(ast_amp))
		self.astamp.setEnabled(True)

		self.astagl.setText(str(ast_agl))
		self.astagl.setEnabled(True)

		self.deferr.setText(str(error_defocus))
		self.deferr.setEnabled(True)

		self.astamperr.setText(str(error_astamp))
		self.astamperr.setEnabled(True)

		self.astaglerr.setText(str(error_astagl))
		self.astaglerr.setEnabled(True)

		# XXX: wgh?? amp_cont static to 0?
		# set image properties, in order to save ctf values
		from utilities import set_ctf
		set_ctf(img, [defocus, ctf_cs, ctf_volt, input_pixel_size, 0, ctf_ampcont, ast_amp, ast_agl])
		# and rewrite image
		img.write_image(image_name)
		print [defocus, ctf_cs, ctf_volt, input_pixel_size, 0, ctf_ampcont, ast_amp, ast_agl]
		# get alternate, and set its ctf
		altimg=BigImageCache.get_object(image_name).get_image(use_alternate=True)
		set_ctf(altimg, [defocus, ctf_cs, ctf_volt, input_pixel_size, 0, ctf_ampcont, ast_amp, ast_agl])
		BigImageCache.get_object(image_name).register_alternate(altimg)
		print [defocus, ctf_cs, ctf_volt, input_pixel_size, 0, ctf_ampcont, ast_amp, ast_agl]
 		print "CTF estimation using CTER done."
 		#print "Estimated defocus value: ", defocus

		##############################################################################
		#### save ctf estimation parameters to db for command line batch processing
		gbdb = db_open_dict(GaussPanel.GDB_NAME)
		ctfdict = {'pixel_input':input_pixel_size,'pixel_output':output_pixel_size,'ctf_window':ctf_window_size,'ctf_edge':ctf_edge_size,'ctf_overlap':ctf_overlap_size,'ctf_ampcont':ctf_ampcont,'ctf_volt':ctf_volt,'ctf_cs':ctf_cs, 'ctf_kboot':ctf_kboot}
		#print "calc_ctf image_name: ", image_name
		if gbdb.has_key(image_name):
			olddict=gbdb[image_name]
			gbdb[image_name] = dict((olddict).items() + ctfdict.items() ) # merge the two dictionaries with conflict resolution resolved in favorr of the latest ctf parameters
		else:
			gbdb[image_name]=ctfdict

		del img
		del altimg



class GaussBoxer:

	CACHE_MAX = 10 # Each image has its last CACHE_MAX SwarmBoxer instance stored (or retrievable) automatically
	PROFILE_MAX = 0.8 # this is a percentage - it stops the profile trough point from going to the end
	REF_NAME = "gauss_ref"
	AUTO_NAME = "gauss_auto"
	WEAK_REF_NAME = "gauss_weak_ref"
	INIT = True
	MVT_THRESHOLD = 200 # a squared distance - used by the movement cache to automatically update using previously supplied user movement data
	GAUSS_BOXERS = "gauss_boxers"
	GAUSS_FWD_BOXERS = "gauss_fwd_boxers"
	GAUSS_USER_MVT = "gauss_user_mvt"
	def __init__(self,particle_diameter=128):
		if GaussBoxer.INIT: # this solved a strange problem with databases
			GaussBoxer.INIT = False
			EMBox.set_box_color(GaussBoxer.REF_NAME,[0,0,0])
			EMBox.set_box_color(GaussBoxer.WEAK_REF_NAME,[0.2,0.2,0.4])
			EMBox.set_box_color(GaussBoxer.AUTO_NAME,[0.4,.9,.4]) # Classic green, always was this one ;)

		self.panel_object = None # maybe it doesn't exist

                self.pixel_input = 1.0
                self.pixel_output = 1.0
                self.frequency_cutoff = 0
                self.window_size_min = 15
                self.gauss_width = 1.0
                self.use_variance = True
                self.invert = False
                self.thr_low = None
                self.thr_hgh = None
		self.gui_mode = False # set this to False to stop any calls to Qt - such as the act of making the cursor busy...

		self.mvt_cache = [] # we have to remember if any of the auto selected boxes were moved, so if the user reboxes then the movements they previously supplied will be applied


	def __del__(self):
		'''
		Closes the template viewer if it exists
		'''
		#if self.template_viewer != None:
		#	self.template_viewer.close()
		#	self.template_viewer = None

	def set_invert_contrast(self,val):
		'''
		Whether or not to invert contrast
		@param val a boolean
		'''
		self.invert = val

	def set_use_variance(self,val):
		'''
		Whether or not to use variance image
		@param val a boolean
		'''
		self.use_variance = val

	def set_gauss_width(self,gwidth):
		'''
		set gauss width, which should be a float
		@param gwidth the gauss width
		'''
		if gwidth != self.gauss_width:
			self.gauss_width = gwidth


	def set_pixel_input(self,inpix):
		'''
		set the input pixel size, which should be a float
		@param inpix the input pixel size
		'''
		if inpix != self.pixel_input:
			self.pixel_input = inpix

	def set_pixel_output(self,outpix):
		'''
		set the output pixel size, which should be a float
		@param inpix the output pixel size
		'''
		if outpix != self.pixel_output:
			self.pixel_output = outpix

	def set_thr_low(self,thrlow):
		'''
		set the low threshold, which should be a float
		@param thrlow the low threshold
		'''
		if thrlow != self.thr_low:
			self.thr_low = thrlow

	def set_thr_hi(self,thrhi):
		'''
		set the high threshold, which should be a float
		@param thrhi the hi threshold
		'''
		if thrhi != self.thr_hgh:
			self.thr_hgh = thrhi

	def auto_box_clicked(self):
		'''
		When the autobox button is clicked then we force an autobox.
		'''
		print 'file to be processed: ', self.target().current_file()
		self.auto_box(self.target().current_file(),force_remove_auto_boxes=True)

	def auto_box(self, imgname, parameter_update=True, force_remove_auto_boxes=False, cache=True):
		'''
		Autoboxing method using gauss convolution method (also used in sxboxer). Core of the implementation is taken from PawelAutoBoxer and modified to fit in the new framework.
		@param image_name the image that we're boxing
		@param parameter_update, should generally be True,but may be False if the current parameters are known to be current (see self.load_from_last_state)
		@param force_remove_auto_boxes if True all of the autoboxed boxes in the EMBoxerModule are removed and the 'entire' image is autoboxed again. This might be False if you know the template has not changed
		@param cache whether or not the newly establish state, i.e. at the end of this function, should be cached to the database and internally. Generally True but sometimes False (see self.load_from_last_state) .
		'''
		from sparx import filt_gaussl
		print "Gauss method............start auto boxing"
		if self.gui_mode:
			from PyQt4 import QtCore
			get_application().setOverrideCursor(QtCore.Qt.BusyCursor)

		# user pawelautoboxer (gauss method) to compute soln
		# downsample input image.
		small_img = self.get_small_image(imgname)
		#set_idd_image_entry(imgname,'subsampled_image',self.small_img)
		BigImageCache.get_object(imgname).register_alternate(small_img)
		#this causes the alternate image to be picked up and update micrograph window using the small image
		self.target().set_current_file(imgname)
		[avg,sigma,fmin,fmax] = Util.infomask(small_img, None, True )
		small_img -= avg
		small_img /= sigma

		if(self.use_variance):
			from morphology import power
			small_img = power(small_img, 2.0)
			print "using variance"

		boxsize = self.target().get_box_size()
		ccf = filt_gaussl( small_img, self.gauss_width/boxsize )
		del small_img
		peaks = ccf.peak_ccf( boxsize/2-1)
		del ccf
		npeak = len(peaks)/3
		print "npeak: ", npeak
		boxes = []
		ccfs = [] # ccfs are used to set threshold_low adn threshold_high after the particles have been picked. see set_data in CcfHistogram in sxboxer and set_params_of_gui in pawelautoboxer in boxertools.py
		print "thr low: ", self.thr_low
		print "thr hi: ", self.thr_hgh
		print "pixel_output: ", self.pixel_output
		print "pixel input: ", self.pixel_input
		print "invert: ", self.invert
		print "gauss width: ", self.gauss_width
		print "variance: ", self.use_variance
		for i in xrange(npeak):
			cx = peaks[3*i+1]
			cy = peaks[3*i+2]

			corr_score= peaks[3*i]
			skip = False
			if not(self.thr_low is None) and corr_score < self.thr_low:
				skip = True

			if not(self.thr_hgh is None) and corr_score > self.thr_hgh:
				skip = True

			if not skip:
				ccfs.append( peaks[3*i] )
				type = GaussBoxer.AUTO_NAME
				box = [cx,cy,type,corr_score]
				boxes.append(box)

		del peaks

		# Need to do: handle npeak when npeak = 1
		if npeak > 1:
			self.panel_object.set_data(ccfs)

		del ccfs

		if self.gui_mode: self.target().set_status_message("Auboxing Done",1000)

		if self.gui_mode: self.target().set_status_message("Updating Positions",0)
		if self.gui_mode: self.target().set_status_message("Done",1000)

		#self.target().set_box_size(boxsize/shrinkby)
		self.target().add_boxes(boxes, True)

		if self.gui_mode:
			from PyQt4 import QtCore
			get_application().setOverrideCursor(QtCore.Qt.ArrowCursor)
			self.target().set_status_message("Autoboxed %d Particles" %len(boxes), 10000)
		else:
			print "Autoboxed %d Particles" %len(boxes)

		self.panel_object.enable_auto_box(False)
		gbdb = js_open_dict(GaussPanel.GDB_NAME)
		gbdb['clear']=False
		gbdb['boxsize']=boxsize
		gbdb['filename']=imgname
		#######################################################
		# save parameters in case user later want to do cmdline autoboxing

		autoboxdict = {'boxsize':boxsize, 'pixel_input':self.pixel_input, 'pixel_output':self.pixel_output, 'invert_contrast':self.invert, 'use_variance':self.use_variance, 'gauss_width':self.gauss_width,'thr_low':self.thr_low,'thr_hi':self.thr_hgh}

		if gbdb.has_key(imgname):
			oldautoboxdict = gbdb[imgname]
			gbdb[imgname] = dict(oldautoboxdict.items() + autoboxdict.items()) # resolve conflicts in favor of new autoboxdict
		else:
			gbdb[imgname] = autoboxdict
		#######################################################
		return boxes


	def clear_all(self):
		'''
		Clears all boxes calculated in Gauss mode and start from clean slate by setting the current file (i.e., the micrograph being picked) to the commandline one. Otherwise, the already normalized micrograph will get normalized again.
		'''

		self.target().clear_boxes([GaussBoxer.REF_NAME,GaussBoxer.AUTO_NAME,GaussBoxer.WEAK_REF_NAME],cache=True)

		self.panel_object.pixel_output_edit.setText(str(self.pixel_output))
		self.panel_object.new_pixel_output()
		self.panel_object.pixel_input_edit.setText(str(self.pixel_input))
		self.panel_object.new_pixel_input()
		self.panel_object.gauss_width.setText(str(self.gauss_width))
		self.panel_object.gauss_width_edited()
		self.panel_object.setgwbox=False
		self.panel_object.gauss_width_slider.setValue(int( log10(self.gauss_width) * 100) )
		self.panel_object.setgwbox=True
		tlowtext = self.panel_object.THRNA
		if self.thr_low != None:
		        tlowtext = str(self.thr_low )
		self.panel_object.thr_low_edit.setText(tlowtext)
		self.panel_object.new_thr_low()
		thitext = self.panel_object.THRNA
		if self.thr_hgh != None:
		        thitext = str(self.thr_hgh )
                self.panel_object.thr_hi_edit.setText(thitext)
                self.panel_object.new_thr_hi()

		self.panel_object.use_variance_chk.setChecked(self.use_variance)
		self.panel_object.use_variance_checked( self.use_variance)
		self.panel_object.invert_contrast_chk.setChecked(self.invert)
		self.panel_object.invert_contrast_checked(self.invert)

		BigImageCache.get_object(self.target().file_names[0]).register_alternate(None)
		self.target().set_current_file_by_idx(0)
		self.target().set_current_file(self.target().file_names[0])
		self.panel_object.enable_auto_box(True)

		gbdb = js_open_dict(GaussPanel.GDB_NAME)
		# Don't set gbdb to None but just set it's 'clear' flag to true.
		# If GUI is invoked next time, alternate image will NOT be generated adn GUI will start from clean slate. However, if gauss mode autoboxing or ctf is invoked via command line, the paramters used for boxing will still be there for autoboxing and ctf est to work.
		gbdb['clear']=True

	def get_small_image(self,imgname,modecmd=False,boxsize=128,ret_dummy=False):

		from sparx import get_im, filt_gaussl, filt_gaussh
		subsample_rate = self.get_subsample_rate()
		frequency_cutoff = self.get_frequency_cutoff()
		template_min = self.get_window_size_min()
		gaussh_param = self.get_gaussh_param(modecmd=modecmd,boxsize=boxsize)
		invert = self.get_invert()

		img = get_im( imgname )

		# first invert image if invert is true. code taken from invert_contrast_mic_toggled in sxboxer.py
		if invert:
			[avg,sigma,fmin,fmax] = Util.infomask( img, None, True )
			img -= avg
			img *= -1
			img += avg

		img_filt = filt_gaussh( img, gaussh_param )

		if subsample_rate != 1.0:
			print "Generating downsampled image\n"
			sb = Util.sincBlackman(template_min, frequency_cutoff,1999) # 1999 taken directly from util_sparx.h
			small_img = img_filt.downsample(sb,subsample_rate)
			del sb
		else:
			small_img = img_filt.copy()

		del img_filt

		small_img.set_attr("invert", invert)
		small_img.set_attr("gaussh_param", gaussh_param)
		small_img.set_attr("subsample_rate",subsample_rate)
		small_img.set_attr("frequency_cutoff",frequency_cutoff)
		small_img.set_attr("template_min",template_min)
		from utilities import generate_ctf
		try:
			ctf_dict = img.get_attr("ctf")
			ctf_dict.apix = self.pixel_output
			small_img.set_attr("ctf",ctf_dict)
		except:
			pass

		del img

		if ret_dummy:
			BigImageCache.get_object(imgname).register_alternate(small_img)
			del small_img
			return None

		return small_img
	#############################################################################
	# parameter access functions from pawelautoboxer class
	def get_subsample_rate(self):
		return self.pixel_input/self.pixel_output

	def get_window_size_min(self):
		return 15

	def get_frequency_cutoff(self):
		return 0.5*self.get_subsample_rate()

	def get_gaussh_param(self,modecmd=False,boxsize=128):
		ratio = self.pixel_input/self.pixel_output
		if modecmd:
			return ratio/boxsize
		return ratio/self.target().get_box_size()

	def get_invert(self):
		return self.invert
	#############################################################################

	def move_auto_box(self,box_number,image_name,dx,dy):
		'''
		This is the solution to a conundrum - what if the user moves an auto box (a green one) but then reboxes, subsequently
		losing this "centering" metadata? The solution implemented here is to store the movement metadata and perform collision
		detection whenever autoboxing happens again.

		Here we update the "movement cache" to reflect the new metadata
		@param box_number the index corresponding to the stored box of the EMBoxList in the EMBoxerModule
		@param image_name the name of the image we're currently operating on - this is important because its used to store data in the local database
		@param dx the movement in the x direction (float)
		@param dy the movement in the y direction
		'''
		box = self.target().get_box(box_number)
		for i,[old,new] in enumerate(self.mvt_cache):
			dist = (new[0]-box.x)**2 + (new[1]-box.y)**2
			if dist < GaussBoxer.MVT_THRESHOLD:
				self.mvt_cache[i][1] = [box.x+dx,box.y+dy]
				break
		else:
			self.mvt_cache.append([[box.x,box.y],[box.x+dx,box.y+dy]])

		set_database_entry(image_name,GaussBoxer.GAUSS_USER_MVT,self.mvt_cache)
		self.target().move_box(box_number,dx,dy)

	def get_alternate(self,filename):
		# if there is a subsampled image in cache then user was probably using that the last time
		# if we use subsampled image from cache then also have to reload all the previous parameters, mosti mportantly box size and in/output pixel size

		#self.small_img = get_idd_image_entry(filename,'subsampled_image')

		# look at parameters stored in cache to determine if users ended previous session not from a blank slate
		# if that is the case then pre-process/downsample input micrograph and set it to alternate

		gbdb = js_open_dict(GaussPanel.GDB_NAME)

		try:
			if gbdb['clear'] == False and gbdb['filename']==filename:
				print "Restoring micrograph and particles window for Gauss mode..."
				# set parameters from cache so get_small image work from the correct parameters
				self.pixel_input = gbdb['pixel_input']
				self.pixel_output = gbdb['pixel_output']
				self.invert = gbdb['invert_contrast']
				self.use_variance = gbdb['use_variance']
				small_img = self.get_small_image(filename)
				BigImageCache.get_object(filename).register_alternate(small_img)
				del small_img
				#this causes the alternate image to be picked up by micrograph window
				self.target().current_idx=0
				self.target().set_current_file(filename)
				self.panel_object.enable_auto_box(False)
		except: pass

	# this is for automated particle picking from command line
	def auto_box_cmdline(self, imgname,boxsize=128,norm="normalize.edgemean"):

		from sparx import filt_gaussl
		print "Gauss method............start command line auto boxing"
		# user pawelautoboxer (gauss method) to compute soln
		# downsample input image.
		small_img = self.get_small_image(imgname,modecmd=True,boxsize=boxsize)
		BigImageCache.get_object(imgname).register_alternate(small_img)
		[avg,sigma,fmin,fmax] = Util.infomask( small_img, None, True )
		small_img -= avg
		small_img /= sigma

		if(self.use_variance):
			from morphology import power
			small_img = power(small_img, 2.0)
			print "using variance"

		ccf = filt_gaussl( small_img, self.gauss_width/boxsize )
		del small_img
		peaks = ccf.peak_ccf( boxsize/2-1)
		del ccf
		npeak = len(peaks)/3
		print "npeak: ", npeak
		boxes = []
		for i in xrange(npeak):
			cx = peaks[3*i+1]
			cy = peaks[3*i+2]

			corr_score= peaks[3*i]
			skip = False
			if not(self.thr_low is None) and corr_score < self.thr_low:
				skip = True

			if not(self.thr_hgh is None) and corr_score > self.thr_hgh:
				skip = True

			if not skip:
				type = GaussBoxer.AUTO_NAME
				box = [cx,cy,type,corr_score]
				boxes.append(box)
		del peaks
		# adds boxes and write to database
		self.target.add_boxes(boxes, True)
		print "Autoboxed %d Particles" %len(boxes)

	# auto_ctf is meant to be called for batch only...
	# take care of case where estimated ctf is saved into downsampled micrograph from which particles are picked.
	def auto_ctf(self,image_name,ctf_params):

		from utilities import get_im
		img = get_im(image_name)
		ctf_volt = ctf_params['ctf_volt']
		ctf_window = ctf_params['ctf_window']
		ctf_overlap = ctf_params['ctf_overlap']
		ctf_edge = ctf_params['ctf_edge']
		ctf_Cs = ctf_params['ctf_cs']
		ctf_ampcont = ctf_params['ctf_ampcont']
		ctf_fstart = ctf_params['ctf_fstart']
		ctf_fstop = ctf_params['ctf_fstop']
		self.pixel_input = ctf_params['pixel_input']
		self.pixel_output = ctf_params['pixel_output']

		from fundamentals import welch_pw2
		# XXX: check image dimensions, especially box size for welch_pw2!
		power_sp = welch_pw2(img,win_size=ctf_window,overlp_x=ctf_overlap,overlp_y=ctf_overlap,edge_x=ctf_edge,edge_y=ctf_edge)

		from fundamentals import rot_avg_table
		avg_sp = rot_avg_table(power_sp)
		del power_sp

		from morphology import defocus_gett
		defocus = defocus_gett(avg_sp,voltage=ctf_volt,Pixel_size=self.pixel_input,Cs=ctf_Cs,wgh=ctf_ampcont, f_start=ctf_fstart,f_stop=ctf_fstop)

		# set image properties, in order to save ctf values
		from utilities import set_ctf, generate_ctf
		ctf_tuple = [defocus,ctf_Cs,ctf_volt,self.pixel_output,0,ctf_ampcont]
		set_ctf(img, ctf_tuple)
		img.write_image(image_name, 0)
		print "        CTF parameters for original micrograph %s:"%image_name, ctf_tuple
		del img

class GaussTool(GaussBoxer,EMBoxingTool):
	'''
	A class that knows how to handle mouse erase events for a GUIBox
	'''

	def __init__(self,target,particle_diameter=128):
		GaussBoxer.__init__(self,particle_diameter)
		self.target = weakref.ref(target) # module can be accessed through self.target. So self.target.get_main_2D_window can be used to get the micrograph window
		self.panel_object = GaussPanel(self,self.pixel_input)
		self.current_file = None # the name of the file that is being studied in the main viewer
		self.moving = None # keeps track of the moving box, will be a list in the format [[int,int],int,str] = [[x,y],box_number,box_type]
		self.ptcl_moving_data = None # keeps track of movement that's occuring in the particles (mx) viewer
		self.gui_mode = True

	def unique_name(self): return "Gauss"

	def get_widget(self):
		if self.panel_object == None:
			self.panel_object = GaussPanel(self,self.particle_diameter)
		return self.panel_object.get_widget()

	def icon(self):
		from PyQt4 import QtGui
		return QtGui.QIcon(get_image_directory() + "swarm_icon.png")

	def get_2d_window(self): return self.target().get_2d_window()

	def mouse_move(self,event):
		pass

	def mouse_wheel(self,event):
		pass

	def mouse_down(self,event):
		# User clicking on micrograph does nothing (currently) in Gauss mode
		return

	def handle_box_delete(self,box,box_num):
		if box.type == GaussBoxer.AUTO_NAME:
			self.target().remove_box(box_num,exclude_region=True)
		elif box.type == GaussBoxer.REF_NAME or box.type == GaussBoxer.WEAK_REF_NAME:
			self.target().remove_box(box_num)
			self.remove_ref(box,self.target().current_file())
		else:
			raise EMUnknownBoxType,box.type

	def mouse_drag(self,event) :
		m=self.get_2d_window().scr_to_img((event.x(),event.y()))
		from PyQt4.QtCore import Qt
		if event.modifiers()&Qt.ShiftModifier:
			box_num = self.target().detect_box_collision(m)
			box = self.target().get_box(box_num)
			if ( box_num != -1):
				if box.type in [GaussBoxer.REF_NAME,GaussBoxer.AUTO_NAME,GaussBoxer.WEAK_REF_NAME]:
					self.handle_box_delete(box,box_num)
		elif self.moving != None:
			oldm = self.moving[0]
			dx = m[0]-oldm[0]
			dy = m[1]-oldm[1]

			if self.moving[2] in [GaussBoxer.REF_NAME,GaussBoxer.WEAK_REF_NAME]:
				self.move_ref(self.moving[1],self.target().current_file(),dx,dy)
			else:
				self.move_auto_box(self.moving[1],self.target().current_file(),dx,dy)

			self.moving[0] = m

	def mouse_up(self,event) :
		if self.moving != None:
			self.target().box_released(self.moving[1])
			if self.moving[2] in [GaussBoxer.REF_NAME,GaussBoxer.WEAK_REF_NAME]:
				self.ref_released(self.target().current_file(),self.moving[1])

		self.moving= None

	def moving_ptcl_established(self,box_num,x,y):
		box = self.target().get_box(box_num)
		if box.type not in [GaussBoxer.REF_NAME,GaussBoxer.AUTO_NAME,GaussBoxer.WEAK_REF_NAME]:
			raise EMUnknownBoxType,box.type

		self.ptcl_moving_data = [x,y,box_num]

	def move_ptcl(self,box_num,x,y,ptcls):
		if self.ptcl_moving_data == None: return

		dx = self.ptcl_moving_data[0] - x
		dy = y - self.ptcl_moving_data[1]
		box = self.target().get_box(box_num)

		if box.type in [GaussBoxer.REF_NAME,GaussBoxer.WEAK_REF_NAME]:
			self.move_ref(box_num,self.target().current_file(),dx,dy)
		else:
			self.move_auto_box(box_num,self.target().current_file(),dx,dy)

		self.ptcl_moving_data = [x,y,self.ptcl_moving_data[2]]

	def release_moving_ptcl(self,box_num,x,y):
		if self.ptcl_moving_data == None: return
		self.target().box_placement_update_exclusion_image_n(box_num)
		box = self.target().get_box(box_num)
		if box.type in [GaussBoxer.REF_NAME,GaussBoxer.WEAK_REF_NAME]:
			self.ref_released(self.target().current_file(),box_num)

		self.ptcl_moving_data = None

	def delete_ptcl(self,box_num):
		box = self.target().get_box(box_num)
		if box.type not in [GaussBoxer.REF_NAME,GaussBoxer.AUTO_NAME,GaussBoxer.WEAK_REF_NAME]:
			raise EMUnknownBoxType,box.type

		self.handle_box_delete(self.target().get_box(box_num),box_num)

	def get_unique_box_types(self):
		return [GaussBoxer.REF_NAME,GaussBoxer.AUTO_NAME,GaussBoxer.WEAK_REF_NAME]

	def boxes_erased(self,list_of_boxes):
		GaussBoxer.boxes_erased(self,list_of_boxes,self.target().current_file())

# this is class CTFInspector from sxboxer.py with very slight modifications
class CTFInspectorWidget(QtGui.QWidget):

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
		from PyQt4 import QtCore
		from PyQt4.QtCore import Qt
		if (self.i_start is None and (i_start_initial > 0)):
			self.i_start = i_start_initial
		if (self.i_stop is None and (i_stop_initial > 0)):
			self.i_stop = i_stop_initial


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
				sizeh = max([max(self.data[i][self.i_start:self.i_stop]) for i in xrange(len(self.data)-1)])
			else:
				sizeh = max([max(self.data[i]) for i in xrange(len(self.data)-1)])


			sizeh = float(sizeh)
			steph = float(h-2*hborder) / float(sizeh)

			import math
			from utilities import read_text_file
			ctfdata2 = read_text_file("procpw.txt",3)

			if ((self.i_start is not None) and (self.i_stop is not None)):
				sizehctf = max(ctfdata2[self.i_start:self.i_stop])

			else:
				sizehctf = max(ctfdata2)

			sizehctf = float(sizehctf)

			tickspacing = min(int(sizew/30)+1, 5)

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
						if (list_index == 2):
							oldy=int(h-hborder-(h-2*hborder)*ctfdata2[index-1]/sizehctf)
						else:
							oldy=int(h-hborder-(h-2*hborder)*self.data[list_index][index-1]/sizeh)
						#newy = int(h-hborder-steph*self.data[list_index][index])
						if (list_index == 2):
							newy=int(h-hborder-(h-2*hborder)*ctfdata2[index]/sizehctf)
						else:
							newy=int(h-hborder-(h-2*hborder)*self.data[list_index][index]/sizeh)

						p.drawLine(oldx,oldy,newx,newy)


					if (len(self.data)-1 == list_index):
						if index % tickspacing == 0:
							p.setPen(Qt.white)
							p.setFont(QtGui.QFont('Times',10))
							p.drawLine(newx, h-hborder, newx, h-hborder+5)
							metrics = p.fontMetrics()
							fw = metrics.width(str(index))
							p.drawText(newx-fw/2, h-hborder+14, str(index))

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

if __name__ == "__main__":
	main()
