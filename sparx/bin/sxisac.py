#!/usr/bin/env python

#
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#

import os
import sys
import random

import global_def
from   global_def import *
from optparse import OptionParser, SUPPRESS_HELP
import ConfigParser
from inspect import currentframe, getframeinfo

# from __future__ import print_function
from EMAN2 import *
from sparx import *
from logger import Logger, BaseLogger_Files
import global_def

from mpi   import  *
from math  import  *

from utilities import send_string_to_all

from applications import  ali2d_base

from utilities import program_state_stack

NAME_OF_JSON_STATE_FILE = "my_state.json"
NAME_OF_ORIGINAL_IMAGE_INDEX = "originalid"
NAME_OF_RUN_DIR = "run"
NAME_OF_MAIN_DIR = "generation_"
DIR_DELIM = os.sep

def main(args):
	progname = os.path.basename(sys.argv[0])
	usage = ( progname + " stack_file  output_directory --radius=particle_radius --img_per_grp=img_per_grp --CTF --restart_section<The remaining parameters are optional --ir=ir --rs=rs --xr=xr --yr=yr --ts=ts --maxit=maxit --dst=dst --FL=FL --FH=FH --FF=FF --init_iter=init_iter --main_maxit=main_iter" +
			" --iter_reali=iter_reali --match_first=match_first --max_round=max_round --match_second=match_second --stab_ali=stab_ali --thld_err=thld_err --indep_run=indep_run --thld_grp=thld_grp" +
			"  --generation=generation  --rand_seed=rand_seed>" )
	
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--radius",                type="int",           help="particle radius: there is no default, a sensible number has to be provided, units - pixels (default required int)")
	parser.add_option("--img_per_grp",           type="int",           default=100,        help="number of images per class: in the ideal case (essentially maximum size of class) (default 100)")
	parser.add_option("--CTF",                   action="store_true",  default=False,      help="apply phase-flip for CTF correction: if set the data will be phase-flipped using CTF information included in image headers (default False)")
	parser.add_option("--ir",                    type="int",           default=1,          help="inner ring: of the resampling to polar coordinates. units - pixels (default 1)")
	parser.add_option("--rs",                    type="int",           default=1,          help="ring step: of the resampling to polar coordinates. units - pixels (default 1)")
	parser.add_option("--xr",                    type="int",           default=-1,         help="x range: of translational search. By default, set by the program. (default -1)")
	parser.add_option("--yr",                    type="int",           default=-1,         help="y range: of translational search. By default, same as xr. (default -1)")
	parser.add_option("--ts",                    type="float",         default=1.0,        help="search step: of translational search: units - pixels (default 1.0)")
	parser.add_option("--maxit",                 type="int",           default=30,         help="number of iterations for reference-free alignment: (default 30)")
	#parser.add_option("--snr",            type="float",        default=1.0,     help="signal-to-noise ratio (only meaningful when CTF is enabled, currently not supported)")
	parser.add_option("--center_method",         type="int",           default=7,          help="method for centering: of global 2D average during initial prealignment of data (0 : no centering; -1 : average shift method; please see center_2D in utilities.py for methods 1-7) (default 7)")
	parser.add_option("--dst",                   type="float",         default=90.0,       help="discrete angle used in within group alignment: (default 90.0)")
	parser.add_option("--FL",                    type="float",         default=0.2,        help="lowest stopband: frequency used in the tangent filter (default 0.2)")
	parser.add_option("--FH",                    type="float",         default=0.3,        help="highest stopband: frequency used in the tangent filter (default 0.3)")
	parser.add_option("--FF",                    type="float",         default=0.2,        help="fall-off of the tangent filter: (default 0.2)")
	parser.add_option("--init_iter",             type="int",           default=3,          help="SAC initialization iterations: number of runs of ab-initio within-cluster alignment for stability evaluation in SAC initialization (default 3)")
	parser.add_option("--main_iter",             type="int",           default=3,          help="SAC main iterations: number of runs of ab-initio within-cluster alignment for stability evaluation in SAC (default 3)")
	parser.add_option("--iter_reali",            type="int",           default=1,          help="SAC stability check interval: every iter_reali iterations of SAC stability checking is performed (default 1)")
	parser.add_option("--match_first",           type="int",           default=1,          help="number of iterations to run 2-way matching in the first phase: (default 1)")
	parser.add_option("--max_round",             type="int",           default=20,         help="maximum rounds: of generating candidate class averages in the first phase (default 20)")
	parser.add_option("--match_second",          type="int",           default=5,          help="number of iterations to run 2-way (or 3-way) matching in the second phase: (default 5)")
	parser.add_option("--stab_ali",              type="int",           default=5,          help="number of alignments when checking stability: (default 5)")
	parser.add_option("--thld_err",              type="float",         default=0.7,        help="threshold of pixel error when checking stability: equals root mean square of distances between corresponding pixels from set of found transformations and theirs average transformation, depends linearly on square of radius (parameter ou). units - pixels. (default 0.7)")
	parser.add_option("--indep_run",             type="int",           default=4,          help="level of m-way matching for reproducibility tests: By default, perform full ISAC to 4-way matching. Value indep_run=2 will restrict ISAC to 2-way matching and 3 to 3-way matching.  Note the number of used MPI processes requested in mpirun must be a multiplicity of indep_run. (default 4)")
	parser.add_option("--thld_grp",              type="int",           default=10,         help="threshold of the size of reproducible class: essentially minimum size of class (default 10)")
	parser.add_option("--n_generations",         type="int",           default=100,        help="maximum number of generations: program stops when reaching this total number of generations: (default 100)")
	#parser.add_option("--candidatesexist",action="store_true", default=False,   help="Candidate class averages exist use them (default False)")
	parser.add_option("--rand_seed",             type="int",           help="random seed set before calculations: useful for testing purposes (default total randomness - type int)")
	parser.add_option("--new",                   action="store_true",  default=False,      help="use new code: (default False)")
	parser.add_option("--debug",                 action="store_true",  default=False,      help="debug info printout: (default False)")

	# must be switched off in production
	parser.add_option("--use_latest_master_directory",action="store_true",  default=False,      help="use latest master directory: when active, the program looks for the latest directory that starts with the word 'master', so the user does not need to provide a directory name. (default False)")
	
	parser.add_option("--restart_section",       type="string",        default=' ',        help="restart section: each generation (iteration) contains three sections: 'restart', 'candidate_class_averages', and 'reproducible_class_averages'. To restart from a particular step, for example, generation 4 and section 'candidate_class_averages' the following option is needed: '--restart_section=candidate_class_averages,4'. The option requires no white space before or after the comma. The default behavior is to restart execution from where it stopped intentionally or unintentionally. For default restart, it is assumed that the name of the directory is provided as argument. Alternatively, the '--use_latest_master_directory' option can be used. (default ' ')")
	parser.add_option("--stop_after_candidates", action="store_true",  default=False,      help="stop after candidates: stops after the 'candidate_class_averages' section. (default False)")

	##### XXXXXXXXXXXXXXXXXXXXXX option does not exist in docs XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	parser.add_option("--return_options", action="store_true", dest="return_options", default=False, help = SUPPRESS_HELP)
	parser.add_option("--skip_alignment",        action="store_true",  default=False,      help="skip alignment step: to be used if images are already aligned. (default False)")

	required_option_list = ['radius']
	(options, args) = parser.parse_args(args)

	if options.return_options:
		return parser
	
	if len(args) > 2:
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
		sys.exit()
	
	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()
	
	from isac import iter_isac
	global_def.BATCH = True

	global_def.BATCH = True
	
	command_line_provided_stack_filename = args[0]
	global_def.BATCH = True

	main_node = 0
	mpi_init(0, [])
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	nproc = mpi_comm_size(MPI_COMM_WORLD)

	# Making sure all required options appeared.
	for required_option in required_option_list:
		if not options.__dict__[required_option]:
			print "\n ==%s== mandatory option is missing.\n"%required_option
			print "Please run '" + progname + " -h' for detailed options"
			return 1

	radi  = options.radius
	center_method  = options.center_method
	if(radi < 1):  ERROR("Particle radius has to be provided!","sxisac",1,myid)

	
	use_latest_master_directory = options.use_latest_master_directory
	stop_after_candidates = options.stop_after_candidates
	# program_state_stack.restart_location_title_from_command_line = options.restart_section
	
	from utilities import qw
	program_state_stack.PROGRAM_STATE_VARIABLES = set(qw("""
		isac_generation
	"""))

	# create or reuse master directory
	masterdir = ""
	stack_processed_by_ali2d_base__filename = ""
	stack_processed_by_ali2d_base__filename__without_master_dir = ""
	error_status = 0
	if len(args) == 2:
		masterdir = args[1]
	elif len(args) == 1:
		if use_latest_master_directory:
			all_dirs = [d for d in os.listdir(".") if os.path.isdir(d)]
			import re; r = re.compile("^master.*$")
			all_dirs = filter(r.match, all_dirs)
			if len(all_dirs)>0:
				# all_dirs = max(all_dirs, key=os.path.getctime)
				masterdir = max(all_dirs, key=os.path.getmtime)
				
	#Create folder for all results or check if there is one created already
	if(myid == main_node):
		if( masterdir == ""):
			timestring = strftime("%Y_%m_%d__%H_%M_%S" + DIR_DELIM, localtime())
			masterdir = "master"+timestring
			cmd = "{} {}".format("mkdir", masterdir)
			cmdexecute(cmd)
		elif not os.path.exists(masterdir):
			# os.path.exists(masterdir) does not exist
			masterdir = args[1]
			cmd = "{} {}".format("mkdir", masterdir)
			cmdexecute(cmd)

		if(args[0][:4] == "bdb:"): filename = args[0][4:]
		else:                      filename = args[0][:-4]
		filename = os.path.basename(filename)
		stack_processed_by_ali2d_base__filename  = "bdb:" + os.path.join(masterdir, filename )
		stack_processed_by_ali2d_base__filename__without_master_dir  = "bdb:" + filename
	if_error_all_processes_quit_program(error_status)

	# send masterdir to all processes
	masterdir = send_string_to_all(masterdir)

	if myid == 0:
		if options.restart_section != " ":
			if os.path.exists(os.path.join(masterdir,NAME_OF_JSON_STATE_FILE)):
				stored_stack, stored_state = restore_program_stack_and_state(os.path.join(masterdir,NAME_OF_JSON_STATE_FILE))
				import re
				if "," in options.restart_section:
					parsed_restart_section_option = options.restart_section.split(",")
					stored_state[-1]["location_in_program"] = re.sub(r"___.*$", "___%s"%parsed_restart_section_option[0], stored_state[-1]["location_in_program"])
					generation_str_format = parsed_restart_section_option[1]
					if generation_str_format != "":
						isac_generation_from_command_line = int(generation_str_format)
						stored_state[-1]["isac_generation"] = isac_generation_from_command_line 
					else:
						isac_generation_from_command_line = 1
						if "isac_generation" in stored_state[-1]:
							del stored_state[-1]["isac_generation"]
				else:
					isac_generation_from_command_line = -1
					stored_state[-1]["location_in_program"] = re.sub(r"___.*$", "___%s"%options.restart_section, stored_state[-1]["location_in_program"])
					if "isac_generation" in stored_state[-1]:
						del stored_state[-1]["isac_generation"]
				store_program_state(os.path.join(masterdir,NAME_OF_JSON_STATE_FILE), stored_state, stored_stack)
			else:
				print "Please remove the restart_section option from the command line. The program must be started from the beginning."			
				mpi_finalize()
				sys.exit()
		else:
			isac_generation_from_command_line = -1
	
	program_state_stack(locals(), getframeinfo(currentframe()), os.path.join(masterdir,NAME_OF_JSON_STATE_FILE))	

	stack_processed_by_ali2d_base__filename = send_string_to_all(stack_processed_by_ali2d_base__filename)
	stack_processed_by_ali2d_base__filename__without_master_dir = \
		send_string_to_all(stack_processed_by_ali2d_base__filename__without_master_dir)

	#  PARAMETERS OF THE PROCEDURE
	if( options.xr == -1 ):
		#  Default values
		target_nx = 76
		target_radius = 29
		target_xr = 1
	else:  #  nx//2
		#  Check below!
		target_xr = options.xr
		target_nx = 76 + target_xr - 1 # subtract one, which is default
		target_radius = 29

	mpi_barrier(MPI_COMM_WORLD)

	# Initialization of stacks
	if(myid == main_node):
		number_of_images_in_stack = EMUtil.get_image_count(command_line_provided_stack_filename)
	else:
		number_of_images_in_stack = 0

	number_of_images_in_stack = bcast_number_to_all(number_of_images_in_stack, source_node = main_node)
	
	nxrsteps = 4
	
	init2dir = os.path.join(masterdir,"2dalignment")

	if(myid == 0):
		import subprocess
		from logger import Logger, BaseLogger_Files
		#  Create output directory
		log2d = Logger(BaseLogger_Files())
		log2d.prefix = os.path.join(init2dir)
		cmd = "mkdir -p "+log2d.prefix
		outcome = subprocess.call(cmd, shell=True)
		log2d.prefix += "/"
		# outcome = subprocess.call("sxheader.py  "+command_line_provided_stack_filename+"   --params=xform.align2d  --zero", shell=True)
	else:
		outcome = 0
		log2d = None

	if(myid == main_node):
		a = get_im(command_line_provided_stack_filename)
		nnxo = a.get_xsize()
	else:
		nnxo = 0
	nnxo = bcast_number_to_all(nnxo, source_node = main_node)

	txrm = (nnxo - 2*(radi+1))//2
	if(txrm < 0):  			ERROR( "ERROR!!   Radius of the structure larger than the window data size permits   %d"%(radi), "sxisac",1, myid)
	if(txrm/nxrsteps>0):
		tss = ""
		txr = ""
		while(txrm/nxrsteps>0):
			tts=txrm/nxrsteps
			tss += "  %d"%tts
			txr += "  %d"%(tts*nxrsteps)
			txrm =txrm//2
	else:
		tss = "1"
		txr = "%d"%txrm

	# section ali2d_base

	params2d, aligned_images = ali2d_base(command_line_provided_stack_filename, init2dir, None, 1, radi, 1, txr, txr, tss, \
		False, 90.0, center_method, 14, options.CTF, 1.0, False, \
		"ref_ali2d", "", log2d, nproc, myid, main_node, MPI_COMM_WORLD, write_headers = False, skip_alignment = options.skip_alignment)

	mpi_barrier(MPI_COMM_WORLD)
	if( myid == main_node ):
		if options.skip_alignment:
			print "========================================="
			print "Even though there is no alignment step, '%s' params are set to zero for later use."%os.path.join(init2dir, "initial2Dparams.txt")
			print "========================================="
		write_text_row(params2d,os.path.join(init2dir, "initial2Dparams.txt"))
	del params2d
	mpi_barrier(MPI_COMM_WORLD)

	#  We assume the target image size will be target_nx, radius will be 29, and xr = 1.  
	#  Note images can be also padded, in which case shrink_ratio > 1.
	shrink_ratio = float(target_radius)/float(radi)
	nx = aligned_images[0].get_xsize()
	nima = len(aligned_images)
	newx = int(nx*shrink_ratio + 0.5)

	from fundamentals import rot_shift2D, resample
	from utilities import pad, combine_params2
	if(shrink_ratio < 1.0):
		if    newx > target_nx  :
			msk = model_circle(target_radius, target_nx, target_nx)
			for im in xrange(nima):
				#  Here we should use only shifts
				alpha, sx, sy, mirror, scale = get_params2D(aligned_images[im])
				alpha, sx, sy, mirror = combine_params2(0, sx,sy, 0, -alpha, 0, 0, 0)
				aligned_images[im] = rot_shift2D(aligned_images[im], 0, sx, sy, 0)
				aligned_images[im]  = resample(aligned_images[im], shrink_ratio)
				aligned_images[im] = Util.window(aligned_images[im], target_nx, target_nx, 1)
				p = Util.infomask(aligned_images[im], msk, False)
				aligned_images[im] -= p[0]
				p = Util.infomask(aligned_images[im], msk, True)
				aligned_images[im] /= p[1]
		elif  newx == target_nx :
			msk = model_circle(target_radius, target_nx, target_nx)
			for im in xrange(nima):
				#  Here we should use only shifts
				alpha, sx, sy, mirror, scale = get_params2D(aligned_images[im])
				alpha, sx, sy, mirror = combine_params2(0, sx,sy, 0, -alpha, 0, 0, 0)
				aligned_images[im] = rot_shift2D(aligned_images[im], 0, sx, sy, 0)
				aligned_images[im]  = resample(aligned_images[im], shrink_ratio)
				p = Util.infomask(aligned_images[im], msk, False)
				aligned_images[im] -= p[0]
				p = Util.infomask(aligned_images[im], msk, True)
				aligned_images[im] /= p[1]
		elif  newx < target_nx  :	
			msk = model_circle(newx//2-2, newx,  newx)
			for im in xrange(nima):
				#  Here we should use only shifts
				alpha, sx, sy, mirror, scale = get_params2D(aligned_images[im])
				alpha, sx, sy, mirror = combine_params2(0, sx,sy, 0, -alpha, 0, 0, 0)
				aligned_images[im] = rot_shift2D(aligned_images[im], 0, sx, sy, 0)
				aligned_images[im]  = resample(aligned_images[im], shrink_ratio)
				p = Util.infomask(aligned_images[im], msk, False)
				aligned_images[im] -= p[0]
				p = Util.infomask(aligned_images[im], msk, True)
				aligned_images[im] /= p[1]
				aligned_images[im] = pad(aligned_images[im], target_nx, target_nx, 1, 0.0)
	elif(shrink_ratio == 1.0):
		if    newx > target_nx  :
			msk = model_circle(target_radius, target_nx, target_nx)
			for im in xrange(nima):
				#  Here we should use only shifts
				alpha, sx, sy, mirror, scale = get_params2D(aligned_images[im])
				alpha, sx, sy, mirror = combine_params2(0, sx,sy, 0, -alpha, 0, 0, 0)
				aligned_images[im] = rot_shift2D(aligned_images[im], 0, sx, sy, 0)
				aligned_images[im] = Util.window(aligned_images[im], target_nx, target_nx, 1)
				p = Util.infomask(aligned_images[im], msk, False)
				aligned_images[im] -= p[0]
				p = Util.infomask(aligned_images[im], msk, True)
				aligned_images[im] /= p[1]
		elif  newx == target_nx :
			msk = model_circle(target_radius, target_nx, target_nx)
			for im in xrange(nima):
				#  Here we should use only shifts
				alpha, sx, sy, mirror, scale = get_params2D(aligned_images[im])
				alpha, sx, sy, mirror = combine_params2(0, sx,sy, 0, -alpha, 0, 0, 0)
				aligned_images[im] = rot_shift2D(aligned_images[im], 0, sx, sy, 0)
				p = Util.infomask(aligned_images[im], msk, False)
				aligned_images[im] -= p[0]
				p = Util.infomask(aligned_images[im], msk, True)
				aligned_images[im] /= p[1]
		elif  newx < target_nx  :			
			msk = model_circle(newx//2-2, newx,  newx)
			for im in xrange(nima):
				#  Here we should use only shifts
				alpha, sx, sy, mirror, scale = get_params2D(aligned_images[im])
				alpha, sx, sy, mirror = combine_params2(0, sx,sy, 0, -alpha, 0, 0, 0)
				aligned_images[im] = rot_shift2D(aligned_images[im], 0, sx, sy, 0)
				#aligned_images[im]  = resample(aligned_images[im], shrink_ratio)
				p = Util.infomask(aligned_images[im], msk, False)
				aligned_images[im] -= p[0]
				p = Util.infomask(aligned_images[im], msk, True)
				aligned_images[im] /= p[1]
				aligned_images[im] = pad(aligned_images[im], target_nx, target_nx, 1, 0.0)
	elif(shrink_ratio > 1.0):
		target_radius = radi
		msk = model_circle(target_radius, nx, nx)
		for im in xrange(nima):
			#  Here we should use only shifts
			alpha, sx, sy, mirror, scale = get_params2D(aligned_images[im])
			alpha, sx, sy, mirror = combine_params2(0, sx,sy, 0, -alpha, 0, 0, 0)
			aligned_images[im] = rot_shift2D(aligned_images[im], 0, sx, sy, 0)
			p = Util.infomask(aligned_images[im], msk, False)
			aligned_images[im] -= p[0]
			p = Util.infomask(aligned_images[im], msk, True)
			aligned_images[im] /= p[1]
			aligned_images[im] = pad(aligned_images[im], target_nx, target_nx, 1, 0.0)
	del msk

	gather_compacted_EMData_to_root(number_of_images_in_stack, aligned_images, myid)
	number_of_images_in_stack = bcast_number_to_all(number_of_images_in_stack, source_node = main_node)

	if( myid == main_node ):
		for i in range(number_of_images_in_stack):  aligned_images[i].write_image(stack_processed_by_ali2d_base__filename,i)
		#  It has to be explicitly closed
		from EMAN2db import db_open_dict
		DB = db_open_dict(stack_processed_by_ali2d_base__filename)
		DB.close()
		
		fp = open("README_shrink_ratio.txt", "w")
		output_text = """
		Since, for processing purposes, isac changes the image dimensions,
		adjustment of pixel size needs to be made in subsequent steps, (e.g.
		running sxviper.py). The shrink ratio for this particular isac run is
		--------
		%.5f
		--------
		To get the pixel size for the isac output the user needs to divide
		the original pixel size by the above value. This info is saved in
		the following file: README_shrink_ratio.txt
		"""%shrink_ratio
		fp.write(output_text); fp.flush() ;fp.close()
		print output_text

	mpi_barrier(MPI_COMM_WORLD)

	global_def.BATCH = True

	os.chdir(masterdir)

	if program_state_stack(locals(), getframeinfo(currentframe())):
	# if 1:
		pass
		if (myid == main_node):
			cmdexecute("sxheader.py  --consecutive  --params=originalid   %s"%stack_processed_by_ali2d_base__filename__without_master_dir)
			cmdexecute("e2bdb.py %s --makevstack=%s_000"%(stack_processed_by_ali2d_base__filename__without_master_dir, stack_processed_by_ali2d_base__filename__without_master_dir))

	if (myid == main_node):
		main_dir_no = get_latest_directory_increment_value("./", NAME_OF_MAIN_DIR, myformat="%04d")
		print "isac_generation_from_command_line", isac_generation_from_command_line, main_dir_no
		if isac_generation_from_command_line < 0:
			if os.path.exists(NAME_OF_JSON_STATE_FILE):
				stored_stack, stored_state = restore_program_stack_and_state(NAME_OF_JSON_STATE_FILE)
				if "isac_generation" in stored_state[-1]:
					isac_generation_from_command_line = stored_state[-1]["isac_generation"]
				else:
					isac_generation_from_command_line = -1
		if isac_generation_from_command_line >= 0 and isac_generation_from_command_line <= main_dir_no: 
			for i in xrange(isac_generation_from_command_line+1, main_dir_no + 1):
				if i == isac_generation_from_command_line+1:
					backup_dir_no = get_nonexistent_directory_increment_value("./", "000_backup", myformat="%05d", start_value=1)
					cmdexecute("mkdir -p " + "000_backup" + "%05d"%backup_dir_no)
				cmdexecute("mv  " + NAME_OF_MAIN_DIR + "%04d"%i +  " 000_backup" + "%05d"%backup_dir_no)
				cmdexecute("rm  " + "EMAN2DB/"+stack_processed_by_ali2d_base__filename__without_master_dir[4:]+"_%03d.bdb"%i)
				
			# it includes both command line and json file
			my_restart_section = stored_state[-1]["location_in_program"].split("___")[-1]
			if "restart" in my_restart_section:
				if "backup_dir_no" not in locals():
					backup_dir_no = get_nonexistent_directory_increment_value("./", "000_backup", myformat="%05d", start_value=1)
					cmdexecute("mkdir -p " + "000_backup" + "%05d"%backup_dir_no)
				cmdexecute("mv  " + NAME_OF_MAIN_DIR + "%04d"%isac_generation_from_command_line +  " 000_backup" + "%05d"%backup_dir_no)
				cmdexecute("rm  " + "EMAN2DB/"+stack_processed_by_ali2d_base__filename__without_master_dir[4:]+"_%03d.bdb"%isac_generation_from_command_line )
			elif "candidate_class_averages" in my_restart_section:
				if "backup_dir_no" not in locals():
					backup_dir_no = get_nonexistent_directory_increment_value("./", "000_backup", myformat="%05d", start_value=1)
					cmdexecute("mkdir -p " + "000_backup" + "%05d"%backup_dir_no)
				cmdexecute("mv  " + NAME_OF_MAIN_DIR + "%04d"%isac_generation_from_command_line +  " 000_backup" + "%05d"%backup_dir_no)
				cmdexecute("mkdir -p " + NAME_OF_MAIN_DIR + "%04d"%isac_generation_from_command_line)
				# cmdexecute("rm -f " + NAME_OF_MAIN_DIR + "%04d/class_averages_candidate*"%isac_generation_from_command_line)
			elif "reproducible_class_averages" in my_restart_section:
				cmdexecute("rm -rf " + NAME_OF_MAIN_DIR + "%04d/ali_params_generation_*"%isac_generation_from_command_line)
				cmdexecute("rm -f " + NAME_OF_MAIN_DIR + "%04d/class_averages_generation*"%isac_generation_from_command_line)
		else:
			if os.path.exists(NAME_OF_JSON_STATE_FILE):
				stored_stack, stored_state = restore_program_stack_and_state(NAME_OF_JSON_STATE_FILE)
				if "isac_generation" in stored_state[-1]:
					isac_generation_from_command_line = stored_state[-1]["isac_generation"]
				else:
					isac_generation_from_command_line = 1
			else:
				isac_generation_from_command_line = 1
	else:
		isac_generation_from_command_line = 0
		
		
		
	isac_generation_from_command_line = mpi_bcast(isac_generation_from_command_line, 1, MPI_INT, 0, MPI_COMM_WORLD)[0]
	isac_generation = isac_generation_from_command_line - 1
	
	if (myid == main_node):
		if isac_generation == 0:
			cmdexecute("mkdir -p " + NAME_OF_MAIN_DIR + "%04d"%isac_generation)
			write_text_file([1], os.path.join(NAME_OF_MAIN_DIR + "%04d"%isac_generation, "generation_%d_accounted.txt"%isac_generation))
			write_text_file(range(number_of_images_in_stack), os.path.join(NAME_OF_MAIN_DIR + "%04d"%isac_generation, "generation_%d_unaccounted.txt"%isac_generation))

	#  Stopping criterion should be inside the program.
	while True:
		isac_generation += 1
		if isac_generation > options.n_generations:
			break

		data64_stack_current = "bdb:../"+stack_processed_by_ali2d_base__filename__without_master_dir[4:]+"_%03d"%isac_generation

		if(myid == main_node):
			accounted_images = read_text_file(os.path.join(NAME_OF_MAIN_DIR + "%04d"%(isac_generation - 1),"generation_%d_accounted.txt"%(isac_generation - 1)))
			number_of_accounted_images = len(accounted_images)
			# unaccounted_images = read_text_file(os.path.join(NAME_OF_MAIN_DIR + "%04d"%(isac_generation - 1),"generation_%d_unaccounted.txt"%(isac_generation - 1)))
			# number_of_unaccounted_images = len(unaccounted_images)
		else:
			number_of_accounted_images = 0

		number_of_accounted_images = int(mpi_bcast(number_of_accounted_images, 1, MPI_INT, 0, MPI_COMM_WORLD)[0])
		if number_of_accounted_images == 0:
			os.chdir("..")
			break

		program_state_stack.restart_location_title = "restart"
		if program_state_stack(locals(), getframeinfo(currentframe())):
			if (myid == main_node):
				cmdexecute("mkdir -p " + NAME_OF_MAIN_DIR + "%04d"%isac_generation)
				# reference the original stack
				list_file = os.path.join(NAME_OF_MAIN_DIR + "%04d"%(isac_generation - 1), "generation_%d_unaccounted.txt"%(isac_generation - 1))
				cmdexecute("e2bdb.py %s --makevstack=%s --list=%s"%(stack_processed_by_ali2d_base__filename__without_master_dir,\
						stack_processed_by_ali2d_base__filename__without_master_dir + "_%03d"%isac_generation, list_file))
			mpi_barrier(MPI_COMM_WORLD)

		os.chdir(NAME_OF_MAIN_DIR + "%04d"%isac_generation)

		program_state_stack.restart_location_title = "candidate_class_averages"
		if program_state_stack(locals(), getframeinfo(currentframe())):

			iter_isac(data64_stack_current, options.ir, target_radius, options.rs, target_xr, target_xr, options.ts, options.maxit, False, 1.0,\
				options.dst, options.FL, options.FH, options.FF, options.init_iter, options.main_iter, options.iter_reali, options.match_first, \
				options.max_round, options.match_second, options.stab_ali, options.thld_err, options.indep_run, options.thld_grp, \
				options.img_per_grp, isac_generation, False, random_seed=options.rand_seed, new=False)#options.new)

		# program_state_stack.restart_location_title = "stopped_program1"
		# program_state_stack(locals(), getframeinfo(currentframe()))
		
		program_state_stack.restart_location_title = "stop_after_candidates"
		program_state_stack(locals(), getframeinfo(currentframe()))
		if stop_after_candidates:
			mpi_finalize()
			sys.exit()

		exit_program = 0
		if(myid == main_node):
			if not os.path.exists("class_averages_candidate_generation_%d.hdf"%isac_generation):
				print "This generation (%d) no class averages were generated!"%isac_generation
				exit_program = 1
		exit_program = int(mpi_bcast(exit_program, 1, MPI_INT, 0, MPI_COMM_WORLD)[0])
		if exit_program:
			os.chdir("..")
			break

		program_state_stack.restart_location_title = "reproducible_class_averages"
		if program_state_stack(locals(), getframeinfo(currentframe())):


			iter_isac(data64_stack_current, options.ir, target_radius, options.rs, target_xr, target_xr, options.ts, options.maxit, False, 1.0,\
				options.dst, options.FL, options.FH, options.FF, options.init_iter, options.main_iter, options.iter_reali, options.match_first, \
				options.max_round, options.match_second, options.stab_ali, options.thld_err, options.indep_run, options.thld_grp, \
				options.img_per_grp, isac_generation, True, random_seed=options.rand_seed, new=False)#options.new)
			pass

		os.chdir("..")

		if (myid == main_node):
			cmdexecute("rm -f class_averages.hdf")
			cpy(["generation_%04d/class_averages_generation_%d.hdf"%(i,i) for i in xrange(1, isac_generation)], "class_averages.hdf")

		# program_state_stack.restart_location_title = "stopped_program2"
		# program_state_stack(locals(), getframeinfo(currentframe()))

	program_state_stack(locals(), getframeinfo(currentframe()), last_call="__LastCall")


	mpi_finalize()

if __name__=="__main__":
	main(sys.argv[1:])


