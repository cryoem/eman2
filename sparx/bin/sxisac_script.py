#!/usr/bin/env python

import os
import global_def
from   global_def import *
from   optparse import OptionParser
import sys, ConfigParser
from inspect import currentframe, getframeinfo


def program_state(current_state, frameinfo, file_name_of_saved_state=None):
	# first call needs to contain file_name_of_saved_state
	# then call is with program_state(locals(), getframeinfo(currentframe()))
	# needs from inspect import currentframe, getframeinfo
	
	from mpi import mpi_comm_rank, mpi_bcast, MPI_COMM_WORLD, mpi_finalize, MPI_INT 
	
	location_in_program = frameinfo.filename + "_" + str(frameinfo.lineno)
	
	if mpi_comm_rank(MPI_COMM_WORLD) == 0:
		if "file_name_of_saved_state" not in program_state.__dict__:
			if file_name_of_saved_state == None:
				print "Must provide the file name of saved state in the first call of the function"
				mpi_finalize()
				sys.exit()
			program_state.file_name_of_saved_state = file_name_of_saved_state
			
			if (os.path.exists(file_name_of_saved_state)):
				import json; f = open(file_name_of_saved_state, 'r')
				program_state.saved_state = json.load(f); f.close()
				program_state.start_executing = 0
			else:
				program_state.start_executing = 1
		else:
			current_state["location_in_program"] = location_in_program
			if program_state.start_executing == 1:
				# import json; f = open(program_state.file_name_of_saved_state, 'w')
				# json.dump(current_state,f); f.close()
				from utilities import store_value_of_simple_vars_in_json_file
				store_value_of_simple_vars_in_json_file(program_state.file_name_of_saved_state, current_state,
				exclude_list_of_vars = ["MPI_COMM_WORLD", "myid"])
			else:
				for key in set(program_state.saved_state) & set(current_state):
					# print key, current_state[key], program_state.saved_state[key]
					if current_state[key] != program_state.saved_state[key]:
						break
				else:
					program_state.start_executing = 1
	else:
		program_state.start_executing = 0
		
	program_state.start_executing = mpi_bcast(program_state.start_executing, 1, MPI_INT, 0, MPI_COMM_WORLD)
	program_state.start_executing = int(program_state.start_executing[0])

	return program_state.start_executing


def main():
	# progname = os.path.basename(sys.argv[0])
	# usage = ( progname + " stack_file --ir=ir --ou=ou --rs=rs --xr=xr --yr=yr --ts=ts --maxit=maxit --dst=dst --FL=FL --FH=FH --FF=FF --init_iter=init_iter --main_maxit=main_iter" +
	# 		" --iter_reali=iter_reali --match_first=match_first --max_round=max_round --match_second=match_second --stab_ali=stab_ali --thld_err=thld_err --indep_run=indep_run --thld_grp=thld_grp" +
	# 		" --img_per_grp=img_per_grp --generation=generation --candidatesexist --rand_seed=rand_seed" )
	# 
	# parser = OptionParser(usage,version=SPARXVERSION)
	# parser.add_option("--ir",             type="int",          default=1,       help="inner ring of the resampling to polar coordinates (1)")
	# parser.add_option("--ou",             type="int",          default=-1,      help="outer ring of the resampling to polar coordinates (max)")
	# parser.add_option("--rs",             type="int",          default=1,       help="ring step of the resampling to polar coordinates (1)")
	# parser.add_option("--xr",             type="float",        default=1.0,     help="x range of translational search (1.0)")
	# parser.add_option("--yr",             type="float",        default=1.0,     help="y range of translational search (1.0)")
	# parser.add_option("--ts",             type="float",        default=1.0,     help="search step of translational search (1.0)")
	# parser.add_option("--maxit",          type="int",          default=30,      help="number of iterations for reference-free alignment (30)")
	# #parser.add_option("--CTF",            action="store_true", default=False,   help="whether to use CTF information (default=False, currently True is not supported)")
	# #parser.add_option("--snr",            type="float",        default=1.0,     help="signal-to-noise ratio (only meaningful when CTF is enabled, currently not supported)")
	# parser.add_option("--dst",            type="float",        default=90.0,    help="discrete angle used in within group alignment ")
	# parser.add_option("--FL",             type="float",        default=0.1,     help="lowest stopband frequency used in the tangent filter (0.1)")
	# parser.add_option("--FH",             type="float",        default=0.3,     help="highest stopband frequency used in the tangent filter (0.3)")
	# parser.add_option("--FF",             type="float",        default=0.2,     help="fall-off of the tangent filter (0.2)")
	# parser.add_option("--init_iter",      type="int",          default=3,       help="number of iterations of ISAC program in initialization (3)")
	# parser.add_option("--main_iter",      type="int",          default=3,       help="number of iterations of ISAC program in main part (3)")
	# parser.add_option("--iter_reali",     type="int",          default=1,       help="number of iterations in ISAC before checking stability (1)")
	# parser.add_option("--match_first",    type="int",          default=1,       help="number of iterations to run 2-way matching in the first phase (1)")
	# parser.add_option("--max_round",      type="int",          default=20,      help="maximum rounds of generating candidate averages in the first phase (20)")
	# parser.add_option("--match_second",   type="int",          default=5,       help="number of iterations to run 2-way (or 3-way) matching in the second phase (5)")
	# parser.add_option("--stab_ali",       type="int",          default=5,       help="number of alignments when checking stability (5)")
	# parser.add_option("--thld_err",       type="float",        default=0.7,     help="the threshold of pixel error when checking stability (0.7)")
	# parser.add_option("--indep_run",      type="int",          default=4,       help="number of indepentdent runs for reproducibility (default=4, only values 2, 3 and 4 are supported (4)")
	# parser.add_option("--thld_grp",       type="int",          default=10,      help="the threshold of size of reproducible class (essentially minimum size of class) (10)")
	# parser.add_option("--img_per_grp",    type="int",          default=100,     help="number of images per group in the ideal case (essentially maximum size of class) (100)")
	# parser.add_option("--generation",     type="int",          default=1,       help="current generation number (1)")
	# parser.add_option("--candidatesexist",action="store_true", default=False,   help="Candidate class averages exist use them (default False)")
	# parser.add_option("--rand_seed",      type="int",          default=None,    help="random seed set before calculations, useful for testing purposes (default None - total randomness)")
	# #parser.add_option("--new",            action="store_true", default=False,   help="use new code (default = False)")
	# (options, args) = parser.parse_args()
	# 
	# if len(args) != 1:
	# 	print "usage: " + usage
	# 	print "Please run '" + progname + " -h' for detailed options"
	# 	sys.exit()
	# 
	# if global_def.CACHE_DISABLE:
	# 	from utilities import disable_bdb_cache
	# 	disable_bdb_cache()
	# 
	# #  The code is exclusively MPI
	# from mpi import mpi_init
	# sys.argv = mpi_init(len(sys.argv),sys.argv)
	# 
	# from isac import iter_isac
	# global_def.BATCH = True

	from mpi import mpi_init, mpi_comm_rank, MPI_COMM_WORLD
	mpi_init(0, [])
	myid = mpi_comm_rank(MPI_COMM_WORLD)


	"""
	(*)
	Reduce images to 64X64
	sxprocess.py bdb:orgstack bdb:stack --changesize --ratio=????
	
	(*)
	sxprocess.py bdb:data bdb:flip_data --phase_flip
	
	(*)
	sxheader.py bdb:flip_data --params=xform.align2d --zero
	
	(*)
	Pre-align the particles in the stack and apply the resulting parameters to create a pre-aligned stack
	
	
	"""
	
	import time

	program_state(locals(), getframeinfo(currentframe()), "my_state.json")

	for i in range(4):
		for j in range(4):
			if program_state(locals(), getframeinfo(currentframe())):
				time.sleep(1)
				f = open("1_%d%d_%d.txt"%(i,j, myid), "w")
				f.close()
			if program_state(locals(), getframeinfo(currentframe())):
				time.sleep(1)
				f = open("2_%d%d_%d.txt"%(i,j, myid), "w")
				f.close()
			if program_state(locals(), getframeinfo(currentframe())):
				time.sleep(1)
				f = open("3_%d%d_%d.txt"%(i,j, myid), "w")
				f.close()
			if program_state(locals(), getframeinfo(currentframe())):
				time.sleep(1)
				f = open("4_%d%d_%d.txt"%(i,j, myid), "w")
				f.close()
			program_state(locals(), getframeinfo(currentframe()))

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
