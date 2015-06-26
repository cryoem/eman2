#!/usr/bin/env python

import os
import sys

import global_def
from   global_def import *
from   optparse import OptionParser
import sys, ConfigParser
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

import string

NAME_OF_ORIGINAL_IMAGE_INDEX = "originalid"
NAME_OF_RUN_DIR = "run"
NAME_OF_MAIN_DIR = "generation_"
DIR_DELIM = os.sep

class ali3d_options:
	ir     = 1
	rs     = 1
	ou     = -1
	xr     = "-1"
	yr     = "-1"
	ts     = "1"
	an     = "-1"
	sym    = "d2"
	delta  = "2"
	npad   = 2
	center = 0
	CTF    = True
	ref_a  = "S"
	snr    = 1.0
	mask3D = "startm.hdf"
	fl     = 0.4
	aa     = 0.1
	initfl = 0.4
	pwreference = "rotpw3i3.txt"


def preparing_test_data():
	
	location_of_pdb_files = "/home/hvoicu/Analysis/rrviper2/single_part_demo/demo_with_cft/"
	location_of_pdb_files = "../../generate_test_data/"
	box_size = 90
	# box_size = 64
	
			
	os.system("mkdir -p generate_test_data")
	os.chdir("generate_test_data")
	for i in range(5):
		os.system("sxpdb2em.py %sgrp%d.pdb tmp%d.hdf --apix=5.2 --box=%d"%(location_of_pdb_files, i, i, box_size))
	os.system("sxcpy.py tmp*.hdf  ttt.hdf")
	os.system("rm tmp*")
	os.system("e2proc3d.py ttt.hdf fmodel_structure.hdf --process=filter.lowpass.tanh:cutoff_abs=0.15:fall_off=0.2")
	os.system("e2proc3d.py ttt.hdf model_structure.hdf --process=filter.lowpass.tanh:cutoff_abs=0.45:fall_off=0.1")
	os.system("rm ttt.hdf")
	os.system("sxprocess.py model_structure.hdf particles mic --generate_projections apix=5.2:CTF=True:boxsize=%d"%box_size)
	os.system("rm -rf ../EMAN2DB")
	os.system("cp -r EMAN2DB ../EMAN2DB")
	os.chdir("..")
	
	

	# os.system("mpirun -np 5 sxcter.py pwrot partres --indir=. --nameroot=mic --micsuffix=hdf --wn=512 --apix=5.2 --Cs=2.0 --voltage=120 --ac=10.0 --f_start=0.02 --f_stop=0.1 --MPI")
		
		


# def program_state_stack(current_state, frameinfo, file_name_of_saved_state=None, last_call="", force_starting_execution = False):
# 	# first call needs to contain file_name_of_saved_state
# 	# then call is with program_state_stack(locals(), getframeinfo(currentframe()))
# 	# needs from inspect import currentframe, getframeinfo
# 	
# 	from mpi import mpi_comm_rank, mpi_bcast, MPI_COMM_WORLD, mpi_finalize, MPI_INT
# 	from utilities import store_value_of_simple_vars_in_json_file, if_error_all_processes_quit_program
# 	
# 	location_in_program = frameinfo.filename + "_" + str(frameinfo.lineno) + "_" + last_call
# 
# 	error_status = 0
# 
# 
# 	if mpi_comm_rank(MPI_COMM_WORLD) == 0:
# 		if "file_name_of_saved_state" not in program_state_stack.__dict__:
# 			if type(file_name_of_saved_state) != type(""):
# 				print "Must provide the file name of saved state as a string in the first call of the function!"
# 				error_status = 1
# 				
# 			program_state_stack.file_name_of_saved_state = file_name_of_saved_state
# 			program_state_stack.counter = 0
# 			
# 			
# 			if (os.path.exists(file_name_of_saved_state)):
# 				import json; f = open(file_name_of_saved_state, 'r')
# 				program_state_stack.saved_state = json.load(f); f.close()
# 				program_state_stack.start_executing = 0
# 			else:
# 				# check to see if file can be created
# 				f = open(file_name_of_saved_state, "w")
# 				f.close()
# 				program_state_stack.start_executing = 1
# 		else:
# 			program_state_stack.counter += 1
# 			current_state["location_in_program"] = location_in_program
# 			if program_state_stack.start_executing == 1 or last_call != "" or force_starting_execution:
# 				# import json; f = open(program_state_stack.file_name_of_saved_state, 'w')
# 				# json.dump(current_state,f); f.close()
# 				store_value_of_simple_vars_in_json_file(program_state_stack.file_name_of_saved_state, current_state,
# 				exclude_list_of_vars = ["MPI_COMM_WORLD", "myid"])
# 				program_state_stack.start_executing = 1
# 			else:
# 				for key in set(program_state_stack.saved_state) & set(current_state):
# 					# print key, current_state[key], program_state_stack.saved_state[key]
# 					if current_state[key] != program_state_stack.saved_state[key]:
# 						break
# 				else:
# 					program_state_stack.start_executing = 1
# 					print "////////////////////////////" 
# 					print "Start executing: ", location_in_program
# 					print "////////////////////////////"
# 		
# 		store_value_of_simple_vars_in_json_file("/home/hvoicu/Analysis/sxisac_scripts/test043/ccc%04d"%program_state_stack.counter, current_state,
# 				exclude_list_of_vars = ["MPI_COMM_WORLD", "myid"])
# 		
# 	else:
# 		program_state_stack.start_executing = 0
# 		
# 	if_error_all_processes_quit_program(error_status)	
# 		
# 	program_state_stack.start_executing = mpi_bcast(program_state_stack.start_executing, 1, MPI_INT, 0, MPI_COMM_WORLD)
# 	program_state_stack.start_executing = int(program_state_stack.start_executing[0])
# 
# 	# print "program_state_stack.start_executing ", program_state_stack.start_executing
# 
# 	return program_state_stack.start_executing




def program_state_stack(full_current_state, frameinfo, file_name_of_saved_state=None, last_call="", force_starting_execution = False):
	# first call needs to contain file_name_of_saved_state
	# then call is with program_state_stack(locals(), getframeinfo(currentframe()))
	# needs: from inspect import currentframe, getframeinfo
	
	from traceback import extract_stack
	from mpi import mpi_comm_rank, mpi_bcast, MPI_COMM_WORLD, mpi_finalize, MPI_INT
	from utilities import if_error_all_processes_quit_program

	def store_program_state(filename, state, stack):
		import json
		with open(filename, "w") as fp:
			json.dump(zip(stack, state), fp, indent = 2)
		fp.close()

	def restore_program_stack_and_state(file_name_of_saved_state):
		import json; f = open(file_name_of_saved_state, 'r')
		saved_state_and_stack = json.load(f); f.close()
		return list(zip(*saved_state_and_stack)[0]), list(zip(*saved_state_and_stack)[1])
	
	def get_current_stack_info():
		return [[x[0], x[2]] for x in extract_stack()[:-2]]

	PROGRAM_STATE_VARIABLES = {"isac_generation", "i", "j"}
	START_EXECUTING_FALSE = 0
	START_EXECUTING_TRUE = 1
	START_EXECUTING_ONLY_ONE_TIME_THEN_REVERT = 2
	
	# error_status = 1
	# if_error_all_processes_quit_program(error_status)
	
	location_in_program = frameinfo.filename + "_" + str(frameinfo.lineno) + "_" + last_call
	
	current_state = {"location_in_program" : location_in_program}
	for var in PROGRAM_STATE_VARIABLES & set(full_current_state) :
		current_state[var] =  full_current_state[var]
	
	current_stack = get_current_stack_info()

	error_status = 0

	# not a real while, an if with the possibility of jumping with break
	while mpi_comm_rank(MPI_COMM_WORLD) == 0:
		if "file_name_of_saved_state" not in program_state_stack.__dict__:
			if type(file_name_of_saved_state) != type(""):
				print "Must provide the file name of saved state as a string in the first call of the function!"
				error_status = 1
				break

			program_state_stack.file_name_of_saved_state = os.getcwd() + os.sep + file_name_of_saved_state
			program_state_stack.counter = 0
			program_state_stack.track_stack = get_current_stack_info()
			program_state_stack.track_state = [dict() for i in xrange(len(program_state_stack.track_stack))]
			program_state_stack.track_state[-1] = current_state

			if (os.path.exists(file_name_of_saved_state)):
				program_state_stack.saved_stack, \
				program_state_stack.saved_state = restore_program_stack_and_state(file_name_of_saved_state)
				program_state_stack.start_executing = START_EXECUTING_FALSE
			else:
				# check to see if file can be created
				f = open(file_name_of_saved_state, "w"); f.close()
				program_state_stack.start_executing = START_EXECUTING_TRUE
		else:
			program_state_stack.counter += 1
			# print "counter: ", program_state_stack.counter
			# if program_state_stack.counter == program_state_stack.CCC:
			# 	error_status = 1
			# 	break

						
			if program_state_stack.start_executing == START_EXECUTING_ONLY_ONE_TIME_THEN_REVERT:
				program_state_stack.start_executing = START_EXECUTING_FALSE
			
			# correct track_state to reflect track_stack 
			for i in xrange(len(current_stack)):
				if i < len(program_state_stack.track_state):
					if program_state_stack.track_stack[i] != current_stack[i]:
						program_state_stack.track_state[i] = dict()
				else:
					# print "i:", i, len(program_state_stack.track_state), len(current_stack), current_stack
					program_state_stack.track_state.append(dict())
			program_state_stack.track_state[i] = current_state
			
			# correct track_stack to reflect current_stack
			program_state_stack.track_stack = current_stack
			
			# if program_state_stack.counter == 68:
			# 	print range(len(current_stack), len(program_state_stack.track_state))
				
			# delete additional elements in track_state so that size of track_state is the same as current_stack  				
			program_state_stack.track_state[len(current_stack):len(program_state_stack.track_state)] = []
			
			if program_state_stack.start_executing == START_EXECUTING_TRUE or last_call != "" or force_starting_execution:
				store_program_state(program_state_stack.file_name_of_saved_state, program_state_stack.track_state, current_stack)
				program_state_stack.start_executing = START_EXECUTING_TRUE
			else:
				if len(program_state_stack.saved_state) >= len(current_stack):
					for i in range(len(program_state_stack.saved_state)):
						if i < len(current_stack):
							if program_state_stack.track_stack[i] == current_stack[i]:
								if program_state_stack.track_state[i] == program_state_stack.saved_state[i]:
									continue
							break
						else:
							program_state_stack.start_executing = START_EXECUTING_ONLY_ONE_TIME_THEN_REVERT
							# print "////////////////////////////" 
							# print "Entering function: ", location_in_program
							# print "////////////////////////////"
							break
					else:
						program_state_stack.start_executing = START_EXECUTING_TRUE
						# print "////////////////////////////" 
						# print "Start executing: ", location_in_program
						# print "////////////////////////////"
		break
	else:
		program_state_stack.start_executing = START_EXECUTING_FALSE
		
	if_error_all_processes_quit_program(error_status)	
		
	program_state_stack.start_executing = mpi_bcast(program_state_stack.start_executing, 1, MPI_INT, 0, MPI_COMM_WORLD)
	program_state_stack.start_executing = int(program_state_stack.start_executing[0])

	# print "program_state_stack.start_executing ", program_state_stack.start_executing

	return program_state_stack.start_executing





def main():
	progname = os.path.basename(sys.argv[0])
	usage = ( progname + " stack_file --ir=ir --ou=ou --rs=rs --xr=xr --yr=yr --ts=ts --maxit=maxit --dst=dst --FL=FL --FH=FH --FF=FF --init_iter=init_iter --main_maxit=main_iter" +
			" --iter_reali=iter_reali --match_first=match_first --max_round=max_round --match_second=match_second --stab_ali=stab_ali --thld_err=thld_err --indep_run=indep_run --thld_grp=thld_grp" +
			" --img_per_grp=img_per_grp --generation=generation --candidatesexist --rand_seed=rand_seed" )
	
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--ir",             type="int",          default=1,       help="inner ring of the resampling to polar coordinates (1)")
	parser.add_option("--ou",             type="int",          default=-1,      help="outer ring of the resampling to polar coordinates (max)")
	parser.add_option("--rs",             type="int",          default=1,       help="ring step of the resampling to polar coordinates (1)")
	parser.add_option("--xr",             type="float",        default=1.0,     help="x range of translational search (1.0)")
	parser.add_option("--yr",             type="float",        default=1.0,     help="y range of translational search (1.0)")
	parser.add_option("--ts",             type="float",        default=1.0,     help="search step of translational search (1.0)")
	parser.add_option("--maxit",          type="int",          default=30,      help="number of iterations for reference-free alignment (30)")
	parser.add_option("--CTF",            action="store_true", default=False,   help="whether to use CTF information (default=False, currently True is not supported)")
	parser.add_option("--snr",            type="float",        default=1.0,     help="signal-to-noise ratio (only meaningful when CTF is enabled, currently not supported)")
	parser.add_option("--dst",            type="float",        default=90.0,    help="discrete angle used in within group alignment ")
	parser.add_option("--FL",             type="float",        default=0.1,     help="lowest stopband frequency used in the tangent filter (0.1)")
	parser.add_option("--FH",             type="float",        default=0.3,     help="highest stopband frequency used in the tangent filter (0.3)")
	parser.add_option("--FF",             type="float",        default=0.2,     help="fall-off of the tangent filter (0.2)")
	parser.add_option("--init_iter",      type="int",          default=3,       help="number of iterations of ISAC program in initialization (3)")
	parser.add_option("--main_iter",      type="int",          default=3,       help="number of iterations of ISAC program in main part (3)")
	parser.add_option("--iter_reali",     type="int",          default=1,       help="number of iterations in ISAC before checking stability (1)")
	parser.add_option("--match_first",    type="int",          default=1,       help="number of iterations to run 2-way matching in the first phase (1)")
	parser.add_option("--max_round",      type="int",          default=20,      help="maximum rounds of generating candidate averages in the first phase (20)")
	parser.add_option("--match_second",   type="int",          default=5,       help="number of iterations to run 2-way (or 3-way) matching in the second phase (5)")
	parser.add_option("--stab_ali",       type="int",          default=5,       help="number of alignments when checking stability (5)")
	parser.add_option("--thld_err",       type="float",        default=0.7,     help="the threshold of pixel error when checking stability (0.7)")
	parser.add_option("--indep_run",      type="int",          default=4,       help="number of indepentdent runs for reproducibility (default=4, only values 2, 3 and 4 are supported (4)")
	parser.add_option("--thld_grp",       type="int",          default=10,      help="the threshold of size of reproducible class (essentially minimum size of class) (10)")
	parser.add_option("--img_per_grp",    type="int",          default=100,     help="number of images per group in the ideal case (essentially maximum size of class) (100)")
	parser.add_option("--generation",     type="int",          default=1,       help="current generation number (1)")
	parser.add_option("--candidatesexist",action="store_true", default=False,   help="Candidate class averages exist use them (default False)")
	parser.add_option("--rand_seed",      type="int",          default=None,    help="random seed set before calculations, useful for testing purposes (default None - total randomness)")
	#parser.add_option("--new",            action="store_true", default=False,   help="use new code (default = False)")

	#options introduced from sxmeridien
	# parser.add_option("--ir",      		type= "int",   default= 1,			help="inner radius for rotational correlation > 0 (set to 1)")
	# parser.add_option("--ou",      		type= "int",   default= -1,			help="outer radius for rotational correlation < int(nx/2)-1 (set to the radius of the particle)")
	# parser.add_option("--rs",      		type= "int",   default= 1,			help="step between rings in rotational correlation >0  (set to 1)" ) 
	# parser.add_option("--xr",      		type="string", default= "-1",		help="range for translation search in x direction, search is +/xr (default 0)")
	# parser.add_option("--yr",      		type="string", default= "-1",		help="range for translation search in y direction, search is +/yr (default = same as xr)")
	# parser.add_option("--ts",      		type="string", default= "1",		help="step size of the translation search in both directions, search is -xr, -xr+ts, 0, xr-ts, xr, can be fractional")
	parser.add_option("--delta",   		type="string", default= "-1",		help="angular step of reference projections during initialization step (default automatically selected based on radius of the structure.)")
	parser.add_option("--an",      		type="string", default= "-1",		help="angular neighborhood for local searches (phi and theta) (Default exhaustive searches)")
	parser.add_option("--center",  		type="float",  default= -1,			help="-1: average shift method; 0: no centering; 1: center of gravity (default=-1)")
	# parser.add_option("--maxit",   		type="int",  	default= 400,		help="maximum number of iterations performed for the GA part (set to 400) ")
	parser.add_option("--outlier_percentile",type="float",    default= 95,	help="percentile above which outliers are removed every iteration")
	parser.add_option("--iteration_start",type="int",    default= 0,		help="starting iteration for rviper, 0 means go to the most recent one (default).")
	# parser.add_option("--CTF",     		action="store_true", default=False,	help="Use CTF (Default no CTF correction)")
	# parser.add_option("--snr",     		type="float",  default= 1.0,		help="Signal-to-Noise Ratio of the data (default 1.0)")
	parser.add_option("--ref_a",   		type="string", default= "S",		help="method for generating the quasi-uniformly distributed projection directions (default S)")
	parser.add_option("--sym",     		type="string", default= "c1",		help="symmetry of the refined structure")
	parser.add_option("--npad",    		type="int",    default= 2,			help="padding size for 3D reconstruction (default=2)")
	parser.add_option("--nsoft",    	type="int",    default= 1,			help="Use SHC in first phase of refinement iteration (default=1, to turn it off set to 0)")
	parser.add_option("--startangles",  action="store_true", default=False,	help="Use orientation parameters in the input file header to jumpstart the procedure")

	#options introduced for the do_volume function
	parser.add_option("--fl",			type="float",	default=0.12,		help="cut-off frequency of hyperbolic tangent low-pass Fourier filte (default 0.12)")
	parser.add_option("--aa",			type="float",	default=0.1,		help="fall-off of hyperbolic tangent low-pass Fourier filter (default 0.1)")
	parser.add_option("--pwreference",	type="string",	default="",			help="text file with a reference power spectrum (default no power spectrum adjustment)")
	parser.add_option("--mask3D",		type="string",	default=None,		help="3D mask file (default a sphere  WHAT RADIUS??)")

	parser.add_option("--use_latest_master_directory", action="store_true", dest="use_latest_master_directory", default=True)


	# options found only in sx_meridien
	# set(['outlier_percentile', 'aa', 'ref_a', 'mask3D', 'pwreference', 'iteration_start', 'npad', 'sym', 'startangles', 'nsoft', 'delta', 'an', 'fl', 'center'])



	(options, args) = parser.parse_args()
	
	if len(args) != 1:
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
		sys.exit()
	
	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()
	
	from isac import iter_isac
	global_def.BATCH = True

	# orgstack = args[0]
	command_line_provided_stack_filename = args[0]

	radi  = options.ou
	global_def.BATCH = True
	ali3d_options.ir     = options.ir
	ali3d_options.rs     = options.rs
	ali3d_options.ou     = options.ou
	ali3d_options.xr     = options.xr
	ali3d_options.yr     = options.yr
	ali3d_options.ts     = options.ts
	ali3d_options.an     = "-1"
	ali3d_options.sym    = options.sym
	ali3d_options.delta  = options.delta
	ali3d_options.npad   = options.npad
	ali3d_options.center = options.center
	ali3d_options.CTF    = options.CTF
	ali3d_options.ref_a  = options.ref_a
	ali3d_options.snr    = options.snr
	ali3d_options.mask3D = options.mask3D
	ali3d_options.pwreference = ""  #   It will have to be turned on after exhaustive done by setting to options.pwreference
	ali3d_options.fl     = 0.4
	ali3d_options.initfl = 0.4
	ali3d_options.aa     = 0.1

	if( ali3d_options.xr == "-1" ):  ali3d_options.xr = "2"


	from mpi import mpi_init, mpi_comm_rank, MPI_COMM_WORLD
	main_node = 0
	mpi_init(0, [])
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	nproc = mpi_comm_size(MPI_COMM_WORLD)
	
	
	use_latest_master_directory = options.use_latest_master_directory	

	# plan for sxisac_script
	"""

	(*)
	sxprocess.py bdb:data bdb:flip_data --phase_flip

	sxali2d get it from meridien
	
	
	applyparameters rotshift2d on all of them
	
	reduce to 64 
	
	make sure it is zoomed in
	
	
	particles on different proc, brinbg them to zero and write to disk.
	
	


	(*)
	Reduce images to 64X64
	sxprocess.py bdb:orgstack bdb:stack --changesize --ratio=????
	
	(*)
	sxheader.py bdb:flip_data --params=xform.align2d --zero
	
	(*)
	Pre-align the particles in the stack and apply the resulting parameters to create a pre-aligned stack
	
	"""
	
	# mpirun -np 1  sxisac_script.py
	
	program_state_stack(locals(), getframeinfo(currentframe()), "my_state.json")

	# if program_state_stack(locals(), getframeinfo(currentframe())):
	# 	if(myid == main_node):
	# 		preparing_test_data()
	# 
	# mpi_barrier(MPI_COMM_WORLD)

	# if program_state_stack(locals(), getframeinfo(currentframe()), force_starting_execution = True):
	# # if 1:		
	# 	pass

	# bdb_stack_location = send_string_to_all(bdb_stack_location)
	# bdb_stack_location_for_ali2d = send_string_to_all(bdb_stack_location_for_ali2d)
		
		
	# global_def.LOGFILE =  os.path.join(masterdir, global_def.LOGFILE)
	print_program_start_information()
	
	# mpi_barrier(mpi_comm)
	# from mpi import mpi_finalize
	# mpi_finalize()
	# print "mpi finalize"
	# from sys import exit
	# exit()
	
	# if program_state_stack(locals(), getframeinfo(currentframe())):
	# 	os.system("sxheader.py  bdb:particles --params=xform.align2d  --zero")
	# 	os.system("sxheader.py  bdb:particles --consecutive  --params=originalid")



	# create or reuse master directory
	masterdir = ""
	stack_processed_by_ali2d_base__filename = ""
	stack_processed_by_ali2d_base__filename__without_master_dir = ""
	error_status = 0
	if len(args) == 2:
		masterdir = args[1]
		if masterdir[-1] != DIR_DELIM:
			masterdir += DIR_DELIM
	elif len(args) == 1:
		if use_latest_master_directory:
			all_dirs = [d for d in os.listdir(".") if os.path.isdir(d)]
			import re; r = re.compile("^master.*$")
			all_dirs = filter(r.match, all_dirs)
			if len(all_dirs)>0:
				# all_dirs = max(all_dirs, key=os.path.getctime)
				masterdir = max(all_dirs, key=os.path.getmtime)
				masterdir += DIR_DELIM
				
	#Create folder for all results or check if there is one created already
	if(myid == main_node):
		if( masterdir == ""):
			timestring = strftime("%Y_%m_%d__%H_%M_%S" + DIR_DELIM, localtime())
			masterdir = "master"+timestring
			cmd = "{} {}".format("mkdir", masterdir)
			cmdexecute(cmd)
		if os.path.exists(masterdir):
			# bdb_stack_location = args[0]
			# org_stack_location = args[0]

			if ':' in args[0]:
				stack_processed_by_ali2d_base__filename = args[0].split(":")[0] + ":" + masterdir + args[0].split(":")[1]
				stack_processed_by_ali2d_base__filename__without_master_dir = args[0].split(":")[0] + ":" + args[0].split(":")[1]
			else:
				filename = os.path.basename(args[0])
				stack_processed_by_ali2d_base__filename  = "bdb:" + masterdir + os.path.splitext(filename)[0]
				stack_processed_by_ali2d_base__filename__without_master_dir  = "bdb:" + os.path.splitext(filename)[0]
		else:
			# os.path.exists(masterdir) does not exist
			ERROR('Output directory does not exist, please change the name and restart the program', "sxrviper", 1)
			error_status = 1

	if_error_all_processes_quit_program(error_status)
				
	# send masterdir to all processes
	masterdir = send_string_to_all(masterdir)
	if masterdir[-1] != DIR_DELIM:
		masterdir += DIR_DELIM

	stack_processed_by_ali2d_base__filename = send_string_to_all(stack_processed_by_ali2d_base__filename)
	stack_processed_by_ali2d_base__filename__without_master_dir = \
		send_string_to_all(stack_processed_by_ali2d_base__filename__without_master_dir)

	

	#  INPUT PARAMETERS
	radi  = options.ou
	global_def.BATCH = True
	ali3d_options.ir     = options.ir
	ali3d_options.rs     = options.rs
	ali3d_options.ou     = options.ou
	ali3d_options.xr     = options.xr
	ali3d_options.yr     = options.yr
	ali3d_options.ts     = options.ts
	ali3d_options.an     = "-1"
	ali3d_options.sym    = options.sym
	ali3d_options.delta  = options.delta
	ali3d_options.npad   = options.npad
	ali3d_options.center = options.center
	ali3d_options.CTF    = options.CTF
	ali3d_options.ref_a  = options.ref_a
	ali3d_options.snr    = options.snr
	ali3d_options.mask3D = options.mask3D
	ali3d_options.pwreference = ""  #   It will have to be turned on after exhaustive done by setting to options.pwreference
	ali3d_options.fl     = 0.4
	ali3d_options.initfl = 0.4
	ali3d_options.aa     = 0.1

	if( ali3d_options.xr == "-1" ):  ali3d_options.xr = "2"
	

	# mpi_barrier(MPI_COMM_WORLD)
	# from mpi import mpi_finalize
	# mpi_finalize()
	# print "mpi finalize"
	# from sys import exit
	# exit()


	# ali2d_base
	if program_state_stack(locals(), getframeinfo(currentframe())):
	# if 1:		

		nproc     = mpi_comm_size(MPI_COMM_WORLD)
		myid      = mpi_comm_rank(MPI_COMM_WORLD)
		main_node = 0
		
		mpi_barrier(MPI_COMM_WORLD)

		# preparing variables	
		nxinit = 64  #int(280*0.3*2)
		cushion = 8  #  the window size has to be at least 8 pixels larger than what would follow from resolution
		nxstep = 4
		projdata = [[model_blank(1,1)],[model_blank(1,1)]]
		
		mempernode = 4.0e9
		
		#  PARAMETERS OF THE PROCEDURE 
		#  threshold error
		thresherr = 0
		fq = 50 # low-freq resolution to which fuse ref volumes. [A]
		
		# Get the pixel size, if none set to 1.0, and the original image size
		if(myid == main_node):
			a = get_im(command_line_provided_stack_filename)
			nnxo = a.get_xsize()
			if( nnxo%2 == 1 ):
				ERROR("Only even-dimensioned data allowed","sxmeridien",1)
				nnxo = -1
			elif( nxinit > nnxo ):
				ERROR("Image size less than minimum permitted $d"%nxinit,"sxmeridien",1)
				nnxo = -1
			else:
				if ali3d_options.CTF:
					i = a.get_attr('ctf')
					pixel_size = i.apix
					fq = pixel_size/fq
				else:
					pixel_size = 1.0
					#  No pixel size, fusing computed as 5 Fourier pixels
					fq = 5.0/nnxo
				del a
		else:
			nnxo = 0
			pixel_size = 1.0
		nnxo = bcast_number_to_all(nnxo, source_node = main_node)
		if( nnxo < 0 ):
			mpi_finalize()
			exit()
		pixel_size = bcast_number_to_all(pixel_size, source_node = main_node)
		fq   = bcast_number_to_all(fq, source_node = main_node)
	
	
		if(radi < 1):  radi = nnxo//2-2
		elif((2*radi+2)>nnxo):  ERROR("Particle radius set too large!","sxmeridien",1,myid)
		ali3d_options.ou = radi


		#  create a vstack from input stack to the local stack in masterdir
		#  Stack name set to default
		# stack = "bdb:"+masterdir+"/rdata"
		# stack = bdb_stack_location + "_rdata"
		# stack_for_ali2d = bdb_stack_location_for_ali2d + "_rdata"
		# stack_processed_by_ali2d_base__filename = bdb_stack_location + "_adata"
		# Initialization of stacks
		if(myid == main_node):
			# if not (os.path.exists(os.path.join(masterdir,"EMAN2DB/rdata.bdb"))):
			# 	if(orgstack[:4] == "bdb:"):	cmd = "{} {} {}".format("e2bdb.py", orgstack,"--makevstack="+stack)
			# 	else:  cmd = "{} {} {}".format("sxcpy.py", orgstack, stack)
			# 	cmdexecute(cmd)
			# 	cmd = "{} {}".format("sxheader.py  --consecutive  --params=originalid", stack)
			# 	cmdexecute(cmd)
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
			cmd = "mkdir "+log2d.prefix
			outcome = subprocess.call(cmd, shell=True)
			log2d.prefix += "/"
			# outcome = subprocess.call("sxheader.py  "+command_line_provided_stack_filename+"   --params=xform.align2d  --zero", shell=True)
		else:
			outcome = 0
			log2d = None
		txrm = (nnxo - 2*(radi+1))//2
		if(txrm < 0):  			ERROR( "ERROR!!   Radius of the structure larger than the window data size permits   %d"%(radi), "sxmeridien",1, myid)
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

		
		params2d, aligned_images = ali2d_base(command_line_provided_stack_filename, init2dir, None, 1, radi, 1, txr, txr, tss, \
					False, 90.0, -1, 14, options.CTF, 1.0, False, \
					"ref_ali2d", "", log2d, nproc, myid, main_node, MPI_COMM_WORLD, write_headers = False)
	
		gather_compacted_EMData_to_root(number_of_images_in_stack, aligned_images, myid)

		if( myid == main_node ):
		
			for i in range(number_of_images_in_stack):
				# aligned_images[i].set_attr("originalid", i)
				aligned_images[i].write_image(stack_processed_by_ali2d_base__filename, i)

			write_text_row(params2d,os.path.join(init2dir, "initial2Dparams.txt"))		
			# outcome = subprocess.call("sxheader.py  "+stack+"   --params=xform.align2d  --import="+os.path.join(init2dir, "initial2Dparams.txt"), shell=True)


		mpi_barrier(MPI_COMM_WORLD)


	# # ali2d_base
	# if program_state_stack(locals(), getframeinfo(currentframe())):
	# # if 1:
	# 	pass
	# 
	# 
	# 
	# 
	# program_state_stack(locals(), getframeinfo(currentframe()), last_call="LastCall")
	# 
	# from mpi import mpi_finalize
	# mpi_finalize()
	# # import  sys
	# sys.exit()


	# do not need to use this
	# # create or reuse master directory
	# masterdir = ""
	# bdb_stack_location = ""
	# bdb_stack_location_for_ali2d = ""
	# error_status = 0
	# if len(args) == 2:
	# 	masterdir = args[1]
	# 	if masterdir[-1] != DIR_DELIM:
	# 		masterdir += DIR_DELIM
	# elif len(args) == 1:
	# 	if use_latest_master_directory:
	# 		all_dirs = [d for d in os.listdir(".") if os.path.isdir(d)]
	# 		import re; r = re.compile("^master.*$")
	# 		all_dirs = filter(r.match, all_dirs)
	# 		if len(all_dirs)>0:
	# 			# all_dirs = max(all_dirs, key=os.path.getctime)
	# 			masterdir = max(all_dirs, key=os.path.getmtime)
	# 			masterdir += DIR_DELIM
	# 			
	# #Create folder for all results or check if there is one created already
	# if(myid == main_node):
	# 	if( masterdir == ""):
	# 		timestring = strftime("%Y_%m_%d__%H_%M_%S" + DIR_DELIM, localtime())
	# 		masterdir = "master"+timestring
	# 		cmd = "{} {}".format("mkdir", masterdir)
	# 		cmdexecute(cmd)
	# 	if os.path.exists(masterdir):
	# 		if ':' in args[0]:
	# 			bdb_stack_location = args[0].split(":")[0] + ":" + masterdir + args[0].split(":")[1]
	# 			bdb_stack_location_for_ali2d = args[0].split(":")[0] + ":" + args[0].split(":")[1]
	# 			org_stack_location = args[0]
	# 
	# 			if(not os.path.exists(os.path.join(masterdir,"EMAN2DB" + DIR_DELIM))):
	# 				# cmd = "{} {}".format("cp -rp EMAN2DB", masterdir, "EMAN2DB" DIR_DELIM)
	# 				# cmdexecute(cmd)
	# 				cmd = "{} {} {}".format("e2bdb.py", org_stack_location,"--makevstack=" + bdb_stack_location + "_rdata")
	# 				cmdexecute(cmd)
	# 			
	# 				from applications import header
	# 				try:
	# 					header(bdb_stack_location + "_rdata", params=NAME_OF_ORIGINAL_IMAGE_INDEX, fprint=True)
	# 					print "Images were already indexed!"
	# 				except KeyError:
	# 					print "Indexing images"
	# 					header(bdb_stack_location + "_rdata", params=NAME_OF_ORIGINAL_IMAGE_INDEX, consecutive=True)
	# 		else:
	# 			filename = os.path.basename(args[0])
	# 			bdb_stack_location = "bdb:" + masterdir + os.path.splitext(filename)[0]
	# 			bdb_stack_location_for_ali2d = "bdb:" + os.path.splitext(filename)[0]
	# 			if(not os.path.exists(os.path.join(masterdir,"EMAN2DB" + DIR_DELIM))):
	# 				cmd = "{} {} {}".format("sxcpy.py  ", args[0], bdb_stack_location + "_rdata")
	# 				cmdexecute(cmd)
	# 			
	# 				from applications import header
	# 				try:
	# 					header(bdb_stack_location + "_rdata", params=NAME_OF_ORIGINAL_IMAGE_INDEX, fprint=True)
	# 					print "Images were already indexed!"
	# 				except KeyError:
	# 					print "Indexing images"
	# 					header(bdb_stack_location + "_rdata", params=NAME_OF_ORIGINAL_IMAGE_INDEX, consecutive=True)
	# 			
	# 			else:
	# 				ERROR('Conflicting information: EMAN2DB exists, but provided *.hdf file', "sxrviper", 1)
	# 				error_status = 1
	# 
	# 	else:
	# 		# os.path.exists(masterdir) does not exist
	# 		ERROR('Output directory does not exist, please change the name and restart the program', "sxrviper", 1)
	# 		error_status = 1
	# 
	# if_error_all_processes_quit_program(error_status)
	# 			
	# 
	# # send masterdir to all processes
	# masterdir = send_string_to_all(masterdir)
	# if masterdir[-1] != DIR_DELIM:
	# 	masterdir += DIR_DELIM
	# 


	global_def.BATCH = True
	
	if program_state_stack(locals(), getframeinfo(currentframe())):
	# if 1:
		pass
		if (myid == main_node):
			cmdexecute("sxheader.py  --consecutive  --params=originalid %s"%stack_processed_by_ali2d_base__filename)
			cmdexecute("e2bdb.py %s --makevstack=%s_001"%(stack_processed_by_ali2d_base__filename, stack_processed_by_ali2d_base__filename))
	
	if program_state_stack(locals(), getframeinfo(currentframe())):
	# if 1:
		pass

	os.chdir(masterdir + DIR_DELIM)

	# for isac_generation in range(1,10):
	for isac_generation in range(1,10):
		
		data64_stack_current = stack_processed_by_ali2d_base__filename__without_master_dir.split(":")[0]
		data64_stack_current += ":../" + stack_processed_by_ali2d_base__filename__without_master_dir.split(":")[1]
		data64_stack_current += "_%03d"%isac_generation    

		data64_stack_next = stack_processed_by_ali2d_base__filename__without_master_dir.split(":")[0]
		data64_stack_next += ":../" + stack_processed_by_ali2d_base__filename__without_master_dir.split(":")[1]
		data64_stack_next += "_%03d"%(isac_generation + 1)
		
		
		if (myid == main_node):
			cmdexecute("mkdir -p " + NAME_OF_MAIN_DIR + "%04d"%isac_generation + DIR_DELIM)
			
		mpi_barrier(MPI_COMM_WORLD)
		os.chdir(NAME_OF_MAIN_DIR + "%04d"%isac_generation + DIR_DELIM)

		if (myid == main_node):
			print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
			print "isac_generation: ", isac_generation
			print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

		if program_state_stack(locals(), getframeinfo(currentframe())):
		# if 1:
			pass

			iter_isac(data64_stack_current, options.ir, options.ou, options.rs, options.xr, options.yr, options.ts, options.maxit, False, 1.0,\
				#options.CTF, options.snr, \
				options.dst, options.FL, options.FH, options.FF, options.init_iter, options.main_iter, options.iter_reali, options.match_first, \
				options.max_round, options.match_second, options.stab_ali, options.thld_err, options.indep_run, options.thld_grp, \
				options.img_per_grp, isac_generation, options.candidatesexist, random_seed=options.rand_seed, new=False)#options.new)
	
		error_status = 0
		if program_state_stack(locals(), getframeinfo(currentframe())):
		# if 1:
			pass

			while (myid == main_node):
				
				# number_of_accounted_images = sum(1 for line in open("generation_%04d/generation_%d_accounted.txt"%(isac_generation, isac_generation))
				# number_of_unaccounted_images = sum(1 for line in open("generation_%04d/generation_%d_unaccounted.txt"%(isac_generation, isac_generation))
				number_of_accounted_images = sum(1 for line in open("generation_%d_accounted.txt"%(isac_generation)))
				number_of_unaccounted_images = sum(1 for line in open("generation_%d_unaccounted.txt"%(isac_generation)))
				
				if number_of_accounted_images == 0:
					error_status = 1
					break
					
				if number_of_unaccounted_images < 2*options.img_per_grp:
					error_status = 1
					break
			
				# cmdexecute("e2bdb.py %s --makevstack=%s --list=%s%04d/generation_%d_unaccounted.txt"%
				# 		   (data64_stack_current, data64_stack_next, NAME_OF_MAIN_DIR, isac_generation, isac_generation))
				cmdexecute("e2bdb.py %s --makevstack=%s --list=this_generation_%d_unaccounted.txt"%
						   (data64_stack_current, data64_stack_next, isac_generation))
			
				break

			if_error_all_processes_quit_program(error_status)

			os.chdir("..")
		
	global_def.BATCH = False

	if program_state_stack(locals(), getframeinfo(currentframe())):
	# if 1:
		pass

	program_state_stack(locals(), getframeinfo(currentframe()), last_call="LastCall")
	


	# import time
	# 
	# program_state_stack(locals(), getframeinfo(currentframe()), "my_state.json")
	# 
	# for i in range(4):
	# 	for j in range(4):
	# 		if program_state_stack(locals(), getframeinfo(currentframe())):
	# 			time.sleep(1)
	# 			f = open("1_%d%d_%d.txt"%(i,j, myid), "w")
	# 			f.close()
	# 		if program_state_stack(locals(), getframeinfo(currentframe())):
	# 			time.sleep(1)
	# 			f = open("2_%d%d_%d.txt"%(i,j, myid), "w")
	# 			f.close()
	# 		if program_state_stack(locals(), getframeinfo(currentframe())):
	# 			time.sleep(1)
	# 			f = open("3_%d%d_%d.txt"%(i,j, myid), "w")
	# 			f.close()
	# 		if program_state_stack(locals(), getframeinfo(currentframe())):
	# 			time.sleep(1)
	# 			f = open("4_%d%d_%d.txt"%(i,j, myid), "w")
	# 			f.close()
	# 		program_state_stack(locals(), getframeinfo(currentframe()))

	# iter_isac(args[0], options.ir, options.ou, options.rs, options.xr, options.yr, options.ts, options.maxit, False, 1.0,\
	# 	#options.CTF, options.snr, \
	# 	options.dst, options.FL, options.FH, options.FF, options.init_iter, options.main_iter, options.iter_reali, options.match_first, \
	# 	options.max_round, options.match_second, options.stab_ali, options.thld_err, options.indep_run, options.thld_grp, \
	# 	options.img_per_grp, options.generation, options.candidatesexist, random_seed=options.rand_seed, new=False)#options.new)
	# global_def.BATCH = False

	from mpi import mpi_finalize
	mpi_finalize()

if __name__ == "__main__":
	main()



# class A:
# 	a = 1
# 	b = 9
# 	
# # import copy	
# # x = copy.deepcopy(A)
# 
# x = A()
# y = A()
# A.a = 100
# 
# x.a = 888
# x.x = 888
# 
# y.a = 5555
# 
# A.a = 11111
# 
# print x.a
# print y.a





	