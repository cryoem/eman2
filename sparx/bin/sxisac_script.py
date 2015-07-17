#!/usr/bin/env python
import os
import sys

import global_def
from   global_def import *
from   optparse import OptionParser
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

import string

NAME_OF_ORIGINAL_IMAGE_INDEX = "originalid"
NAME_OF_RUN_DIR = "run"
NAME_OF_MAIN_DIR = "generation_"
DIR_DELIM = os.sep

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



def main():
	import sys
	progname = os.path.basename(sys.argv[0])
	usage = ( progname + " stack_file  <output_directory> --radius=particle_radius --img_per_grp=img_per_grp --CTF <The remaining parameters are optional --ir=ir --rs=rs --xr=xr --yr=yr --ts=ts --maxit=maxit --dst=dst --FL=FL --FH=FH --FF=FF --init_iter=init_iter --main_maxit=main_iter" +
			" --iter_reali=iter_reali --match_first=match_first --max_round=max_round --match_second=match_second --stab_ali=stab_ali --thld_err=thld_err --indep_run=indep_run --thld_grp=thld_grp" +
			"  --generation=generation --candidatesexist --rand_seed=rand_seed>" )
	
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--radius",         type="int",          default=-1,      help="Particle radius, it has to be provided.")
	parser.add_option("--img_per_grp",    type="int",          default=100,     help="number of images per group in the ideal case (essentially maximum size of class) (100)")
	parser.add_option("--CTF",            action="store_true", default=False,   help="CTF flag, if set the data will be phase-flipped")
	parser.add_option("--ir",             type="int",          default=1,       help="inner ring of the resampling to polar coordinates (1)")
	parser.add_option("--rs",             type="int",          default=1,       help="ring step of the resampling to polar coordinates (1)")
	parser.add_option("--xr",             type="int",          default=-1,      help="x range of translational search (By default set by the program)")
	parser.add_option("--yr",             type="int",          default=-1,      help="y range of translational search (same as xr")
	parser.add_option("--ts",             type="float",        default=1.0,     help="search step of translational search (1.0)")
	parser.add_option("--maxit",          type="int",          default=30,      help="number of iterations for reference-free alignment (30)")
	#parser.add_option("--snr",            type="float",        default=1.0,     help="signal-to-noise ratio (only meaningful when CTF is enabled, currently not supported)")
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
	parser.add_option("--generation",     type="int",          default=1,       help="current generation number (1)")
	parser.add_option("--candidatesexist",action="store_true", default=False,   help="Candidate class averages exist use them (default False)")
	parser.add_option("--rand_seed",      type="int",          default=None,    help="random seed set before calculations, useful for testing purposes (default None - total randomness)")
	parser.add_option("--new",            action="store_true", default=False,   help="use new code (default = False)")
	parser.add_option("--debug",          action="store_true", default=False,   help="debug info printout (default = False)")


	parser.add_option("--use_latest_master_directory", action="store_true", dest="use_latest_master_directory", default=True)


	# options found only in sx_meridien
	# set(['outlier_percentile', 'aa', 'ref_a', 'mask3D', 'pwreference', 'iteration_start', 'npad', 'sym', 'startangles', 'nsoft', 'delta', 'an', 'fl', 'center'])



	(options, args) = parser.parse_args()
	
	if len(args) > 2:
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
		sys.exit()
	
	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()
	
	from isac import iter_isac
	global_def.BATCH = True

	import cProfile
	import profile
	
	
	global_def.BATCH = True
	
	
	# orgstack = args[0]
	command_line_provided_stack_filename = args[0]

	radi  = options.radius
	
	if(radi < 1):  ERROR("Particle radius has to be provided!","sxisac",1,myid)

	global_def.BATCH = True

	from mpi import mpi_init, mpi_comm_rank, MPI_COMM_WORLD
	main_node = 0
	mpi_init(0, [])
	myid = mpi_comm_rank(MPI_COMM_WORLD)
	nproc = mpi_comm_size(MPI_COMM_WORLD)
	
	use_latest_master_directory = options.use_latest_master_directory


	from utilities import qw
	program_state_stack.PROGRAM_STATE_VARIABLES = set(qw("""
		isac_generation
		i j
		generation
		
		mloop
	"""))
	program_state_stack(locals(), getframeinfo(currentframe()), "my_state.json")


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
				masterdir += DIR_DELIM
				
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
	if_error_all_processes_quit_program(error_status, report_program_state=True)

	# send masterdir to all processes
	masterdir = send_string_to_all(masterdir)

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

	# mpi_barrier(MPI_COMM_WORLD)
	# from mpi import mpi_finalize
	# mpi_finalize()
	# print "mpi finalize"
	# from sys import exit
	# exit()

	# ali2d_base
	program_state_stack.restart_location_title = "ali2d_base"
	if program_state_stack(locals(), getframeinfo(currentframe())):
	# if 1:		


		nproc     = mpi_comm_size(MPI_COMM_WORLD)
		myid      = mpi_comm_rank(MPI_COMM_WORLD)
		main_node = 0
	
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
			cmd = "mkdir "+log2d.prefix
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

		#  centering method is set to #7
		params2d, aligned_images = ali2d_base(command_line_provided_stack_filename, init2dir, None, 1, radi, 1, txr, txr, tss, \
					False, 90.0, 7, 14, options.CTF, 1.0, False, \
					"ref_ali2d", "", log2d, nproc, myid, main_node, MPI_COMM_WORLD, write_headers = False)

		if( myid == main_node ):
			write_text_row(params2d,os.path.join(init2dir, "initial2Dparams.txt"))

		#  We assume the target image size will be target_nx, radius will be 29, and xr = 1.  Note images can be also padded, in which case shrink_ratio > 1.
		shrink_ratio = float(target_radius)/float(radi)
		nx = aligned_images[0].get_xsize()
		nima = len(aligned_images)
		newx = int(nx*shrink_ratio + 0.5)
		if    newx > target_nx  : opt = 1
		elif  newx == target_nx : opt = 0
		else                    : opt = -1

		if opt == -1 :   msk = model_circle(nx//2-2, nx,  nx)
		else:            msk = model_circle(target_radius, target_nx, target_nx)

		from fundamentals import rot_shift2D, resample
		from utilities import pad, combine_params2
		for im in xrange(nima):
			#  Here we should use only shifts
			alpha, sx, sy, mirror = combine_params2(0, params2d[im][1], params2d[im][2], 0, -params2d[im][0], 0, 0, 0)
			aligned_images[im] = rot_shift2D(aligned_images[im], 0, sx, sy, 0)
			if shrink_ratio < 1.0:
				aligned_images[im]  = resample(aligned_images[im], shrink_ratio)
			if   opt == 1 or opt == 0:
				if opt == 1:  aligned_images[im] = Util.window(aligned_images[im], target_nx, target_nx, 1)
				p = Util.infomask(aligned_images[im], msk, False)
				aligned_images[im] -= p[0]
				p = Util.infomask(aligned_images[im], msk, True)
				aligned_images[im] /= p[1]
			elif opt == -1:
				#  Different mask!
				p = Util.infomask(aligned_images[im], msk, False)
				aligned_images[im] -= p[0]
				p = Util.infomask(aligned_images[im], msk, True)
				aligned_images[im] /= p[1]					
				aligned_images[im] = pad(aligned_images[im], target_nx, target_nx, 1, 0.0)
				
				
		del msk, params2d

		gather_compacted_EMData_to_root(number_of_images_in_stack, aligned_images, myid)

		if( myid == main_node ):
			for i in range(number_of_images_in_stack):  aligned_images[i].write_image(stack_processed_by_ali2d_base__filename,i)

	
		mpi_barrier(MPI_COMM_WORLD)

	"""
	from mpi import mpi_finalize
	mpi_finalize()
	import  sys
	sys.exit()
	"""

	global_def.BATCH = True
	
	if program_state_stack(locals(), getframeinfo(currentframe())):
	# if 1:
		pass
		if (myid == main_node):
			cmdexecute("sxheader.py  --consecutive  --params=originalid   %s"%stack_processed_by_ali2d_base__filename)
			cmdexecute("e2bdb.py %s --makevstack=%s_001"%(stack_processed_by_ali2d_base__filename, stack_processed_by_ali2d_base__filename))

	if program_state_stack(locals(), getframeinfo(currentframe()), force_starting_execution = True):
	# if 1:
		if(myid == 0): print  "  STARTED OK"
		pass

	if(myid == 0): print "Location: A   " + os.getcwd()

	os.chdir(masterdir)
	if(myid == 0): print "Location: B   " + os.getcwd()

	# import random
	# image_list  = random.shuffle(range(number_of_images_in_stack))
	
	# for isac_generation in range(1,10):
	isac_generation = 0
	#  Stopping criterion should be inside the program.
	for q12345 in xrange(10):
		isac_generation += 1

		data64_stack_current = "bdb:../"+stack_processed_by_ali2d_base__filename__without_master_dir[4:]+"_%03d"%isac_generation

		data64_stack_next    = "bdb:../"+stack_processed_by_ali2d_base__filename__without_master_dir[4:]+"_%03d"%(isac_generation + 1)
		
		
		if (myid == main_node):
			if os.path.exists(NAME_OF_MAIN_DIR + "%04d"%isac_generation):
				backup_dir_no = get_latest_directory_increment_value(".", "000_backup", myformat="%05d") + 1
				main_dir_no = get_latest_directory_increment_value(".", NAME_OF_MAIN_DIR, myformat="%04d")
				cmdexecute("mkdir -p " + "000_backup" + "%05d"%backup_dir_no)
				for i in xrange(main_dir_no):
					cmdexecute("mv  " + NAME_OF_MAIN_DIR + "%04d"%i +  " 000_backup" + "%05d"%backup_dir_no)
			else:
				cmdexecute("mkdir -p " + NAME_OF_MAIN_DIR + "%04d"%isac_generation)
			
		mpi_barrier(MPI_COMM_WORLD)
		os.chdir(NAME_OF_MAIN_DIR + "%04d"%isac_generation)

		if (myid == main_node):
			print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
			print "isac_generation: ", isac_generation
			print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

		if program_state_stack(locals(), getframeinfo(currentframe())):
			iter_isac(data64_stack_current, options.ir, target_radius, options.rs, target_xr, target_xr, options.ts, options.maxit, False, 1.0,\
				options.dst, options.FL, options.FH, options.FF, options.init_iter, options.main_iter, options.iter_reali, options.match_first, \
				options.max_round, options.match_second, options.stab_ali, options.thld_err, options.indep_run, options.thld_grp, \
				options.img_per_grp, isac_generation, options.candidatesexist, random_seed=options.rand_seed, new=False)#options.new)
			pass
			
			# this_generation_members_acc = image_list[:len(image_list)/5]
			# this_generation_members_unacc = image_list[len(image_list)/5:]

			# write_text_file(this_generation_members_acc, "this_generation_%d_accounted.txt"%isac_generation)
			# write_text_file(this_generation_members_unacc, "this_generation_%d_unaccounted.txt"%isac_generation)
			
			# image_list  = random.shuffle(range(this_generation_members_unacc))
			
	
		error_status = 0
		if program_state_stack(locals(), getframeinfo(currentframe())):
			while(myid == main_node):

				# number_of_accounted_images = sum(1 for line in open("generation_%04d/generation_%d_accounted.txt"%(isac_generation, isac_generation))
				# number_of_unaccounted_images = sum(1 for line in open("generation_%04d/generation_%d_unaccounted.txt"%(isac_generation, isac_generation))
				number_of_accounted_images = sum(1 for line in open("this_generation_%d_accounted.txt"%(isac_generation)))
				number_of_unaccounted_images = sum(1 for line in open("this_generation_%d_unaccounted.txt"%(isac_generation)))
				
				if number_of_accounted_images == 0:
					error_status = 1
					pr.disable()
					s = StringIO.StringIO()
					sortby = 'cumulative'
					ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
					ps.print_stats()
					print s.getvalue()
					break
					
				if number_of_unaccounted_images < 2*options.img_per_grp:
					error_status = 1
					pr.disable()
					s = StringIO.StringIO()
					sortby = 'cumulative'
					ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
					ps.print_stats()
					print s.getvalue()
					break
			
				# cmdexecute("e2bdb.py %s --makevstack=%s --list=%s%04d/generation_%d_unaccounted.txt"%
				# 		   (data64_stack_current, data64_stack_next, NAME_OF_MAIN_DIR, isac_generation, isac_generation))
				cmdexecute("e2bdb.py %s --makevstack=%s --list=this_generation_%d_unaccounted.txt"%
						   (data64_stack_current, data64_stack_next, isac_generation))
				break

		if_error_all_processes_quit_program(error_status, report_program_state=True)

		os.chdir("..")

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

