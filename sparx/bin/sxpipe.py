#!/usr/bin/env python
from __future__ import print_function
#
# Author: Toshio Moriya 02/15/2017 (toshio.moriya@mpi-dortmund.mpg.de)
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete SPHIRE and EMAN2 software packages have some GPL dependencies,
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
# ========================================================================================
# Imports
# ========================================================================================
# Python Standard Libraries
import sys
import os
import argparse

# SPHIRE/EMAN2 Libraries
from EMAN2 import *
from sparx import *
import global_def
from global_def import *

# ========================================================================================
# Helper Functions
# ========================================================================================
# ----------------------------------------------------------------------------------------
# Generate command line
# ----------------------------------------------------------------------------------------
def get_cmd_line():
	cmd_line = ""
	for arg in sys.argv:
		cmd_line += arg + "  "
	cmd_line = "Shell line command: " + cmd_line
	return cmd_line

# ----------------------------------------------------------------------------------------
# Print progress message with time stamp
# ----------------------------------------------------------------------------------------
def print_progress(message):
	from time import strftime, localtime
	time_stamp = strftime("%Y-%m-%d_%H:%M:%S", localtime()) + " =>"
	print(time_stamp, message)

# ----------------------------------------------------------------------------------------
# Get suffix of current time stamp
# ----------------------------------------------------------------------------------------
def get_time_stamp_suffix():
	from time import strftime, localtime
	time_stamp_suffix = strftime("%Y%m%d_%H%M%S", localtime())
	return time_stamp_suffix

# ----------------------------------------------------------------------------------------
# MPI run class
# ----------------------------------------------------------------------------------------
class SXmpi_run(object):
	# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	# static class variables
	# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	RUNNING_UNDER_MPI = False
	main_mpi_proc = 0
	my_mpi_proc_id = 0
	n_mpi_procs = 1
	# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	
	# Set up MPI related variables
	@staticmethod
	def setup():
		# Detect if program is running under MPI
		SXmpi_run.RUNNING_UNDER_MPI = "OMPI_COMM_WORLD_SIZE" in os.environ
	
		SXmpi_run.main_mpi_proc = 0
		if SXmpi_run.RUNNING_UNDER_MPI:
			from mpi import mpi_init
			from mpi import MPI_COMM_WORLD, mpi_comm_rank, mpi_comm_size
		
			mpi_init(0, [])
			SXmpi_run.my_mpi_proc_id = mpi_comm_rank(MPI_COMM_WORLD)
			SXmpi_run.n_mpi_procs = mpi_comm_size(MPI_COMM_WORLD)
		else:
			assert (SXmpi_run.my_mpi_proc_id == 0)
			assert (SXmpi_run.n_mpi_procs == 1)
	
	@staticmethod
	def cleanup():
		if SXmpi_run.RUNNING_UNDER_MPI:
			from mpi import MPI_COMM_WORLD, mpi_barrier, mpi_finalize
			mpi_barrier(MPI_COMM_WORLD)
			mpi_finalize()
	
	@staticmethod
	def is_main_proc():
		return (SXmpi_run.my_mpi_proc_id == SXmpi_run.main_mpi_proc)

# ========================================================================================
# Subcommand functions
# ========================================================================================
# ----------------------------------------------------------------------------------------
# TEST COMMAND
# cd /home/moriya/mrk_qa/mrktest_pipeline
# 
# rm -r debug_mrkout_pipe03o04_sxpipe_isac_substack_isac; sxpipe.py isac_substack 'bdb:mrkout_pipe02_sxwindow#data_isac' 'mrkout_pipe03_sxisac' 'debug_mrkout_pipe03o04_sxpipe_isac_substack_isac' --isac_class_avgs_path='mrkout_pipe03_sxisac/class_averages.hdf'
# e2iminfo.py 'bdb:debug_mrkout_pipe03o04_sxpipe_isac_substack_isac#isac_substack' --header
# sxtransform2d.py 'bdb:debug_mrkout_pipe03o04_sxpipe_isac_substack_isac#isac_substack' 'bdb:debug_mrkout_pipe03o04_sxpipe_isac_substack_isac#isac_substack_shift_applied' --shift --ignore_mirror
# cd debug_mrkout_pipe03o04_sxpipe_isac_substack_isac
# e2display.py &
# cd ..
# 
# rm -r debug_mrkout_pipe03o05_sxpipe_isac_substack_beautifier; sxpipe.py isac_substack 'bdb:mrkout_pipe02_sxwindow#data_isac' 'debug_mrkout_pipe03o01_sxbeautifier' 'debug_mrkout_pipe03o05_sxpipe_isac_substack_beautifier' --isac_class_avgs_path='mrkout_pipe03_sxisac/class_averages.hdf'
# e2iminfo.py 'bdb:debug_mrkout_pipe03o05_sxpipe_isac_substack_beautifier#isac_substack' --header
# 
# rm -r debug_mrkout_pipe03o06_sxpipe_isac_substack_isac; sxpipe.py isac_substack 'bdb:mrkout_pipe02_sxwindow#data_isac' 'mrkout_pipe03_sxisac' 'debug_mrkout_pipe03o06_sxpipe_isac_substack_isac'
# e2iminfo.py 'bdb:debug_mrkout_pipe03o06_sxpipe_isac_substack_isac#isac_substack' --header
# 
# rm -r debug_mrkout_pipe03o07_sxpipe_isac_substack_beautifier; sxpipe.py isac_substack 'bdb:mrkout_pipe02_sxwindow#data_isac' 'debug_mrkout_pipe03o01_sxbeautifier' 'debug_mrkout_pipe03o07_sxpipe_isac_substack_beautifier'
# e2iminfo.py 'bdb:debug_mrkout_pipe03o07_sxpipe_isac_substack_beautifier#isac_substack' --header
# 
# rm -r debug_mrkout_pipe03o08_sxpipe_isac_substack_isac; sxpipe.py isac_substack 'bdb:mrkout_pipe02_sxwindow#data_isac' 'mrkout_pipe03_sxisac' 'debug_mrkout_pipe03o08_sxpipe_isac_substack_isac' --isac_class_avgs_path='mrkout_pipe03_sxisac/ordered_class_averages_two.hdf'
# e2iminfo.py 'bdb:debug_mrkout_pipe03o08_sxpipe_isac_substack_isac#isac_substack' --header
# sxtransform2d.py 'bdb:debug_mrkout_pipe03o08_sxpipe_isac_substack_isac#isac_substack' 'bdb:debug_mrkout_pipe03o08_sxpipe_isac_substack_isac#isac_substack_shift_applied' --shift --ignore_mirror
# e2iminfo.py 'bdb:debug_mrkout_pipe03o08_sxpipe_isac_substack_isac#isac_substack_shift_applied'
# cd debug_mrkout_pipe03o08_sxpipe_isac_substack_isac
# e2display.py &
# cd ..
# 
# ----------------------------------------------------------------------------------------
def isac_substack(args):
	from utilities import get_im, read_text_row, write_text_row, write_text_file, combine_params2, cmdexecute
	from EMAN2db import db_check_dict
	# from EMAN2db import db_open_dict, db_check_dict
	# from e2bdb import makerelpath
	
	# Check MPI execution
	if SXmpi_run.n_mpi_procs > 1:
		assert (SXmpi_run.RUNNING_UNDER_MPI)
		error_status = ("The {} subcommand supports only a single process.".format(subcommand_name), getframeinfo(currentframe()))
		if_error_then_all_processes_exit_program(error_status)
	
	# To make the execution exit upon fatal error by ERROR in global_def.py
	global_def.BATCH = True 
	
	# Check error conditions of arguments
	# subcommand_name = "isac_substack"
	command_script_basename = os.path.basename(sys.argv[0])
	subcommand_name = "{} {}".format(command_script_basename, args.subcommand)
	
	args.input_bdb_stack_path = args.input_bdb_stack_path.strip()
	if not db_check_dict(args.input_bdb_stack_path, readonly=True):
		ERROR("Input BDB image stack file does not exist. Please check the file path and restart the program.", subcommand_name) # action=1 - fatal error, exit
	args.input_run_dir = args.input_run_dir.strip()
	if not os.path.exists(args.input_run_dir):
		ERROR("ISAC or Beautifier run output directory does not exist. Please check the directory path and restart the program.", subcommand_name) # action=1 - fatal error, exit
	args.output_directory = args.output_directory.strip()
	if os.path.exists(args.output_directory):
		ERROR("Output directory exists. Please change the name and restart the program.", subcommand_name) # action=1 - fatal error, exit
	
	assert (db_check_dict(args.input_bdb_stack_path, readonly=True))
	assert (os.path.exists(args.input_run_dir))
	assert (not os.path.exists(args.output_directory))
	
	# Check error conditions of options
	defalut_isac_class_avgs_path = os.path.join(args.input_run_dir, "ordered_class_averages.hdf")
	if args.isac_class_avgs_path != "": # User provided name
		args.isac_class_avgs_path = args.isac_class_avgs_path.strip()
		if not os.path.exists(args.isac_class_avgs_path):
			ERROR("The specifed ISAC class average stack file does not exist. Please check the file path and restart the program.", subcommand_name) # action=1 - fatal error, exit
	else: # Default name of ISAC or Beautifier
		args.isac_class_avgs_path = defalut_isac_class_avgs_path
		if not os.path.exists(args.isac_class_avgs_path):
			ERROR("ISAC or Beautifier run output directory does not contain the default ISAC class average stack file ({}). Please check the directory path or specify ISAC class average stack file, then restart the program.".format(args.isac_class_avgs_path), subcommand_name) # action=1 - fatal error, exit
	args.substack_basename = args.substack_basename.strip()
	if args.substack_basename == "":
		ERROR("Substack basename cannot be empty string or only white spaces.", subcommand_name) # action=1 - fatal error, exit
	
	assert (os.path.exists(args.isac_class_avgs_path))
	assert (args.substack_basename != "")
	
	# Create output directory
	print(" ")
	print_progress("Creating output directory {}.".format(args.output_directory))
	os.mkdir(args.output_directory)
	
	# Extract the number of images in the input BDB stack
	n_fullstack_img = EMUtil.get_image_count(args.input_bdb_stack_path)
	if n_fullstack_img == 0:
		ERROR("Input BDB image stack file does contain no images. Please check the file path and restart the program.", subcommand_name) # action=1 - fatal error, exit
	assert (n_fullstack_img > 0)
	
	# ------------------------------------------------------------------------------------
	# Find out if the input ISAC directory is one of ISAC run or Beautifier run
	# ------------------------------------------------------------------------------------
	# Define enumerators for indices of items in ISAC and Beautifier path
	i_enum = -1
	i_enum += 1; idx_path_item_type = i_enum
	i_enum += 1; idx_path_item_path = i_enum
	i_enum += 1; n_idx_path_item    = i_enum

	# Define enumerators for indices of ISAC run subdirecty/file paths necessary for this command.
	i_enum = -1
	i_enum += 1; idx_isac_path_subdir_align2d                     = i_enum
	i_enum += 1; idx_isac_path_file_shrink                        = i_enum
	i_enum += 1; idx_isac_path_file_fullstack_prealign2d          = i_enum
	i_enum += 1; idx_isac_path_file_fullstack_shrunk_core_align2d = i_enum
	i_enum += 1; n_idx_isac_path                                  = i_enum

	# Define a list of ISAC paths necessary for this command.
	# 
	# The 2dalignment/initial2Dparams.txt contains the scaled-back 2D ISAC pre-alignment parameters (only shift parameters are valid and without scale) of all particles in fullstack. 
	# In the 1st generation, the ISAC core 2D alignment read this file and uses the parameters as an initial starting point (applies the alignment to internal shrunk image).
	# Note each alignment parameters does not include the ISAC core 2D alignment results.
	# 
	# The all_parameters.txt contains the ISAC core 2D alignment parameters (without scale) of all particles in fullstack.
	# From 2nd generation, the ISAC core 2D alignment read this file and uses the parameters as an initial starting point of each generation.
	# Note each alignment parameters does not include the 2D ISAC pre-alignment results, 
	# because the ISAC core alignment is done with prealignment-applied shrunk images.
	# 
	isac_path_list = [None] * n_idx_isac_path
	isac_path_list[idx_isac_path_subdir_align2d]                     = ["Subdirectory", os.path.join(args.input_run_dir, "2dalignment")]
	isac_path_list[idx_isac_path_file_shrink]                        = ["File", os.path.join(args.input_run_dir, "README_shrink_ratio.txt")]
	isac_path_list[idx_isac_path_file_fullstack_prealign2d]          = ["File", os.path.join(args.input_run_dir, "2dalignment", "initial2Dparams.txt")]
	isac_path_list[idx_isac_path_file_fullstack_shrunk_core_align2d] = ["File", os.path.join(args.input_run_dir, "all_parameters.txt")]
	assert (len(isac_path_list) == n_idx_isac_path)
	
	# Check if the contents of run output directory match with ISAC
	isac_missing_path_list = []
	for isac_path in isac_path_list:
		if not os.path.exists(isac_path[idx_path_item_path]):
			isac_missing_path_list.append(isac_path[idx_path_item_path])
	
	# Define enumerators for indices of Beautifier run subdirecty/file paths necessary for this command.
	i_enum = -1
	i_enum += 1; idx_beautifier_path_file_accounted_isac_total_align2d  = i_enum
	i_enum += 1; idx_beautifier_path_file_accounted_local_total_align2d = i_enum
	i_enum += 1; n_idx_beautifier_path                        = i_enum

	# Define a list of Beautifier paths necessary for this command.
	# 
	# The init_isac_params.txt contains the scaled-back totally-combined ISAC 2D alignment parameters of all accounted particles with fullset particle ID. 
	# The 2D alignment parameters (without fullset particle ID) of individual classes are separately stored in params_avg/params_avg_*.txt files. 
	# Each entry is only scaled-back totally-combined ISAC alignment which is used as initial starting point of Beautifier local alignment (see ali2d_single_iter() in alignment.py).
	# 
	# The ali2d_local_params.txt contains the Beautifier 2D local alignment parameters of all accounted particles with fullset particle ID. 
	# The 2D alignment parameters (without fullset particle ID) of individual classes are separately stored in ali2d_local_params_avg/ali2d_local_params_avg_*.txt files. 
	# Each entry is totally-combined alignement of scaled-back totally-combined ISAC alignment and Beautifier local alignment (see ali2d_single_iter() in alignment.py).
	# 
	beautifier_path_list = [None] * n_idx_beautifier_path
	beautifier_path_list[idx_beautifier_path_file_accounted_isac_total_align2d]  = ["File", os.path.join(args.input_run_dir, "init_isac_params.txt")]   
	beautifier_path_list[idx_beautifier_path_file_accounted_local_total_align2d] = ["File", os.path.join(args.input_run_dir, "ali2d_local_params.txt")] 
	assert (len(beautifier_path_list) == n_idx_beautifier_path)
	
	# Check if the contents of run output directory match with Beautifier
	beautifier_missing_path_list = []
	for beautifier_path in beautifier_path_list:
		if not os.path.exists(beautifier_path[idx_path_item_path]):
			beautifier_missing_path_list.append(isac_path[idx_path_item_path])

	# Detect the type of run. (ISAC run directory is prioritized over Beautifier run.)
	if len(isac_missing_path_list) > 0 and len(beautifier_missing_path_list) > 0:
		error_message = "The provided ISAC or Beautifier run output directory is not valid...\n"
		error_message = "    For ISAC, the directory must contains:\n"
		for isac_missing_path in isac_missing_path_list:
			error_message += "    {} {}\n".format(isac_missing_path[idx_path_item_path], isac_missing_path[idx_path_item_type])
		error_message = "    For Beautifier, the directory must contains:\n"
		for beautifier_missing_path in beautifier_missing_path_list:
			error_message += "    {} {}\n".format(beautifier_missing_path[idx_path_item_path], beautifier_missing_path[idx_path_item_type])
		ERROR(error_message, subcommand_name) # action=1 - fatal error, exit
	assert (len(isac_missing_path_list) == 0 or len(beautifier_missing_path_list) == 0)
	
	# Define enumerators for indices of 2D alignment parameters header entry (xform.align2d)
	i_enum = -1
	i_enum += 1; idx_header_align2d_alpha  = i_enum # 2D rotation (in-plane rotation)
	i_enum += 1; idx_header_align2d_tx     = i_enum # x-translation or x-shift
	i_enum += 1; idx_header_align2d_ty     = i_enum # y-translation or y-shift
	i_enum += 1; idx_header_align2d_mirror = i_enum # mirror
	i_enum += 1; idx_header_align2d_scale  = i_enum # scale
	i_enum += 1; n_idx_header_align2d      = i_enum

	# Define enumerators for indices of ISAC 2D alignment parameter
	# Note: ISAC does not stores "scale" to output files
	i_enum = -1
	i_enum += 1; idx_isac_align2d_alpha  = i_enum # 2D rotation (in-plane rotation)
	i_enum += 1; idx_isac_align2d_tx     = i_enum # x-translation or x-shift
	i_enum += 1; idx_isac_align2d_ty     = i_enum # y-translation or y-shift
	i_enum += 1; idx_isac_align2d_mirror = i_enum # mirror
	i_enum += 1; n_idx_isac_align2d      = i_enum

	# Define enumerators for indices of Beautifier 2D alignment parameter
	# Note: Beautifier additionally stores particle image ID to output files
	i_enum = -1
	i_enum += 1; idx_beautifier_align2d_fullstack_img_id  = i_enum # fullstack particle images ID
	i_enum += 1; idx_beautifier_align2d_alpha             = i_enum # 2D rotation (in-plane rotation)
	i_enum += 1; idx_beautifier_align2d_tx                = i_enum # x-translation or x-shift
	i_enum += 1; idx_beautifier_align2d_ty                = i_enum # y-translation or y-shift
	i_enum += 1; idx_beautifier_align2d_mirror            = i_enum # mirror
	i_enum += 1; idx_beautifier_align2d_scale             = i_enum # scale
	i_enum += 1; n_idx_beautifier_align2d                 = i_enum

	# ------------------------------------------------------------------------------------
	# Find out if the input ISAC directory is one of ISAC run or Beautifier run
	# ------------------------------------------------------------------------------------
	# Initialize fullset list with invalid 2D alignment parameters
	invalid_alpha  = 0.0
	invalid_sx     = 0.0
	invalid_sy     = 0.0
	invalid_mirror = -1
	fullstack_total_align2d_list = [[invalid_alpha, invalid_sx, invalid_sy, invalid_mirror]] * n_fullstack_img
	n_accounted_img = 0
	accounted_total_align2d_list = []
	subdir_path = None
	align2d_avg_basename = None
	if len(isac_missing_path_list) == 0:
		assert(len(beautifier_missing_path_list) > 0)
		# The specified run directory is ISAC. (ISAC run directory is prioritized over Beautifier run.)
		print(" ")
		print_progress("ISAC run output directory is specified. The program assumes the ISAC class averages are shrunk and not beautified with the original image size.")
		print_progress("Extracting the shrink ratio and 2D alingment parameters of this ISAC run...")
		
		# shrink ratio
		isac_shrink_path = isac_path_list[idx_isac_path_file_shrink][idx_path_item_path]
		assert (os.path.exists(isac_shrink_path))
		isac_shrink_file = open(isac_path_list[idx_isac_path_file_shrink][idx_path_item_path], "r")
		isac_shrink_lines = isac_shrink_file.readlines()
		isac_shrink_ratio = float(isac_shrink_lines[5])  # 6th line: shrink ratio (= [target particle radius]/[particle radius]) used in the ISAC run
		isac_radius = float(isac_shrink_lines[6])        # 7th line: particle radius at original pixel size used in the ISAC run
		isac_shrink_file.close()
		print_progress("Extracted parameter values...")
		print_progress("  ISAC shrink ratio    : {}".format(isac_shrink_ratio))
		print_progress("  ISAC particle radius : {}".format(isac_radius))
		
		# Pre-alignment (initial 2D alignment) parameters
		fullstack_prealign2d_path = isac_path_list[idx_isac_path_file_fullstack_prealign2d][idx_path_item_path]
		assert (os.path.exists(fullstack_prealign2d_path))
		fullstack_prealign2d_list = read_text_row(fullstack_prealign2d_path)
		print(" ")
		print_progress("Found {} entries in {}.".format(len(fullstack_prealign2d_list), fullstack_prealign2d_path))
		if len(fullstack_prealign2d_list) != n_fullstack_img:
			ERROR("The number of entries in {} is not consistent with {}. Please check the consistency of input datasets.".format(fullstack_prealign2d_path, args.input_bdb_stack_path), subcommand_name) # action=1 - fatal error, exit
		
		# ISAC 2D alignment parameters
		fullstack_shrunk_core_align2d_path = isac_path_list[idx_isac_path_file_fullstack_shrunk_core_align2d][idx_path_item_path]
		assert (os.path.exists(fullstack_shrunk_core_align2d_path))
		fullstack_shrunk_core_align2d_list = read_text_row(fullstack_shrunk_core_align2d_path)
		print(" ")
		print_progress("Found {} entries in {}.".format(len(fullstack_shrunk_core_align2d_list), fullstack_shrunk_core_align2d_path))
		if len(fullstack_shrunk_core_align2d_list) != n_fullstack_img:
			ERROR("The number of entries in {} is not consistent with {}. Please check the consistency of input datasets.".format(fullstack_shrunk_core_align2d_path, args.input_bdb_stack_path), subcommand_name) # action=1 - fatal error, exit
		
		assert (len(fullstack_prealign2d_list) == len(fullstack_shrunk_core_align2d_list))
		
		# For each entry (2D alignment parameters of particle image), register sxcaled back and combined 2D alignment parameters of this ISAC run to the lists
		print(" ")
		print_progress("Registering scaled back and combined 2D alignment parameters of this ISAC run...")
		for fullstack_img_id in xrange(n_fullstack_img):
			prealign2d = fullstack_prealign2d_list[fullstack_img_id]
			if len(prealign2d) != n_idx_isac_align2d:
				ERROR("Invalid number of columns {} at entry #{} in {}. It should be {}. The parameter file might be corrupted. Please consider to rerun ISAC.".format(len(prealign2d), fullstack_img_id, fullstack_prealign2d_path, n_idx_isac_align2d), subcommand_name) # action=1 - fatal error, exit
			shrunk_core_align2d = fullstack_shrunk_core_align2d_list[fullstack_img_id]
			if len(shrunk_core_align2d) != n_idx_isac_align2d:
				ERROR("Invalid number of columns {} at entry #{} in {}. It should be {}. The parameter file might be corrupted. Please consider to rerun ISAC.".format(len(shrunk_core_align2d), fullstack_img_id, fullstack_shrunk_core_align2d_path, n_idx_isac_align2d), subcommand_name) # action=1 - fatal error, exit
			if shrunk_core_align2d[idx_isac_align2d_mirror] != -1: # An accounted particle
				alpha1  = float(prealign2d[idx_isac_align2d_alpha])
				sx1     = float(prealign2d[idx_isac_align2d_tx])
				sy1     = float(prealign2d[idx_isac_align2d_ty])
				mirror1 = int(prealign2d[idx_isac_align2d_mirror])
				alpha2  = float(shrunk_core_align2d[idx_isac_align2d_alpha])
				sx2     = float(shrunk_core_align2d[idx_isac_align2d_tx])/isac_shrink_ratio # Need to apply the shrink ratio to ISAC x-shift
				sy2     = float(shrunk_core_align2d[idx_isac_align2d_ty])/isac_shrink_ratio # Need to apply the shrink ratio to ISAC y-shift
				mirror2 = int(shrunk_core_align2d[idx_isac_align2d_mirror])
				isac_total_align2d = list(combine_params2(alpha1, sx1, sy1, mirror1, alpha2, sx2, sy2, mirror2)) # return value is tuple type but we want to list! 
				
				fullstack_total_align2d_list[fullstack_img_id] = isac_total_align2d
				assert(len(fullstack_total_align2d_list[fullstack_img_id]) == n_idx_isac_align2d)
				scale = 1.0 # because this 2D alignment parameters are scaled back!
				accounted_total_align2d_list.append([fullstack_img_id, isac_total_align2d[idx_isac_align2d_alpha], isac_total_align2d[idx_isac_align2d_tx], isac_total_align2d[idx_isac_align2d_ty], isac_total_align2d[idx_isac_align2d_mirror], scale])
				assert(len(accounted_total_align2d_list[-1]) == n_idx_beautifier_align2d)
			else: # An unaccounted particle
				assert (shrunk_core_align2d[idx_isac_align2d_mirror] == -1)
				assert (len(fullstack_total_align2d_list[fullstack_img_id]) == n_idx_isac_align2d)
				assert (fullstack_total_align2d_list[fullstack_img_id][idx_isac_align2d_mirror] == -1) # Leave this entry to the initialisation value of invalid 2D alignment parameters in fullset list
		
		# Set the number of accounted images
		n_accounted_img = len(accounted_total_align2d_list)
		
		# Set subdirectory name for ISAC run case.
		# Use the corresponding subdirectory name corresponding to Beautifier output directory structure which stores the same information.
		subdir_path = os.path.join(args.output_directory, "params_avg")
		align2d_avg_basename = "params_avg"
		
	else:
		# The specified run directory is Beautifier.
		assert (len(beautifier_missing_path_list) == 0 and len(isac_missing_path_list) > 0)
		print(" ")
		print_progress("Beautifier run output directory is specified. The program assumes the ISAC class averages are beautified with the original image size.")
		print_progress("Extracting the 2D alingment parameters of this Beautifier run...")
	
		# local alignment parameters
		accounted_local_total_align2d_path = beautifier_path_list[idx_beautifier_path_file_accounted_local_total_align2d][idx_path_item_path]
		assert (os.path.exists(accounted_local_total_align2d_path))
		accounted_local_total_align2d_list = read_text_row(accounted_local_total_align2d_path)
		n_accounted_img = len(accounted_local_total_align2d_list)
		print(" ")
		print_progress("Found {} entries in {}.".format(n_accounted_img, accounted_local_total_align2d_path))
		if n_accounted_img > n_fullstack_img:
			ERROR("The number of entries in {} is not consistent with {} (the number of accounted particles is larger than ones of particles in the original fullstack). Please check the consistency of input datasets.".format(accounted_local_total_align2d_path, args.input_bdb_stack_path), subcommand_name) # action=1 - fatal error, exit
		assert(n_accounted_img <= n_fullstack_img)
		
		# For each entry (2D alignment parameters of accounted particle image), register 2D alignment parameters of this Beautifier run to the lists
		print(" ")
		print_progress("Registering 2D alignment parameters of this Beautifier run...")
		for accounted_img_id in xrange(n_accounted_img):
			local_total_param2d = accounted_local_total_align2d_list[accounted_img_id]
			if len(local_total_param2d) != n_idx_beautifier_align2d:
				ERROR("Invalid number of columns {} at entry #{} in {}. It should be {}. The parameter file might be corrupted. Please consider to rerun ISAC.".format(len(local_total_param2d), accounted_img_id, accounted_local_total_align2d_path, n_idx_beautifier_align2d), subcommand_name) # action=1 - fatal error, exit
			if local_total_param2d[idx_beautifier_align2d_mirror] == -1: # An Unaccounted Particle
				ERROR("Invalid alignment parameters of an unaccounted particle is detected at entry #{} in {}. The parameter files might be corrupted. Please consider to rerun Beautifier.".format(accounted_img_id, accounted_local_total_align2d_path), subcommand_name) # action=1 - fatal error, exit
			assert (local_total_param2d[idx_beautifier_align2d_mirror] != -1)
			
			fullstack_img_id  = int(local_total_param2d[idx_beautifier_align2d_fullstack_img_id])
			alpha             = float(local_total_param2d[idx_beautifier_align2d_alpha])
			sx                = float(local_total_param2d[idx_beautifier_align2d_tx])
			sy                = float(local_total_param2d[idx_beautifier_align2d_ty])
			mirror            = int(local_total_param2d[idx_beautifier_align2d_mirror])
			scale             = float(local_total_param2d[idx_beautifier_align2d_scale])
			
			fullstack_total_align2d_list[fullstack_img_id] = [alpha, sx, sy, mirror]
			assert(len(fullstack_total_align2d_list[fullstack_img_id]) == n_idx_isac_align2d)
			accounted_total_align2d_list.append([fullstack_img_id, alpha, sx, sy, mirror, scale])
			assert(len(accounted_total_align2d_list[-1]) == n_idx_beautifier_align2d)
			
		# Set subdirectory name for Beautifier run case.
		# Use the corresponding subdirectory name corresponding to Beautifier output directory structure which stores the same information.
		subdir_path = os.path.join(args.output_directory, "ali2d_local_params_avg")
		align2d_avg_basename = "ali2d_local_params_avg"
		
	if len(fullstack_total_align2d_list) == 0:
		ERROR("No alignment parameters are detected. Please check the contents of ISAC or Beautifier run output directory.", subcommand_name) # action=1 - fatal error, exit
	assert (len(fullstack_total_align2d_list) != 0 and len(fullstack_total_align2d_list) == n_fullstack_img)

	if len(accounted_total_align2d_list) == 0:
		ERROR("No alignment parameters of accounted particles are detected. Please check the contents of ISAC or Beautifier run output directory.", subcommand_name) # action=1 - fatal error, exit
	assert (len(accounted_total_align2d_list) != 0 and len(accounted_total_align2d_list) == n_accounted_img and n_accounted_img <= n_fullstack_img)

	# Save the 2D alignment parameters of all particles to file, using the same format as ISAC 2D alignment file (all_parameters.txt)
	fullset_total_align2d_path = os.path.join(args.output_directory, "scaled_all_parameters.txt")
	write_text_row(fullstack_total_align2d_list, fullset_total_align2d_path)
	print(" ")
	print_progress("Saved the total 2D alignment parameters of all particles in original fullstack to {}, using the same format as ISAC 2D alignment file.".format(fullset_total_align2d_path))
	
	# Save the 2D alignment parameters of all accounted particles to file, using the same format as Beautifier 2D alignment file
	accounted_total_align2d_path = os.path.join(args.output_directory, "init_isac_params.txt")
	write_text_row(accounted_total_align2d_list, accounted_total_align2d_path)
	print(" ")
	print_progress("Saved the total 2D alignment parameters of all accounted particles to {}, using the same format as Beautifier 2D alignment file.".format(accounted_total_align2d_path))
	
	# Create subdirectory
	assert (subdir_path is not None)
	if not os.path.exists(subdir_path): 
		print(" ")
		print_progress("Creating output subdirectory {}.".format(subdir_path))
		os.mkdir(subdir_path)
	
	# Check the number of default ISAC class averages in ISAC or Beautifier run
	print(" ")
	print_progress("Checking the number of default ISAC class averages {} in ISAC or Beautifier run output directory...".format(defalut_isac_class_avgs_path))
	n_default_class_avg = 0
	if os.path.exists(defalut_isac_class_avgs_path):
		n_default_class_avg = EMUtil.get_image_count(defalut_isac_class_avgs_path)
	else: 
		print_progress("WARNING! The default ISAC class averages file does not exist.")
		assert (n_default_class_avg == 0)
	print(" ")
	print_progress("Detected {} default ISAC class averages in {}".format(n_default_class_avg, defalut_isac_class_avgs_path))
	
	# Retrieve original fullstack particle IDs of members listed in ISAC class averages
	print(" ")
	print_progress("Extracting original fullstack particle IDs of members listed in ISAC class averages...")
	n_class_avg = EMUtil.get_image_count(args.isac_class_avgs_path)
	print(" ")
	print_progress("Detected {} ISAC class averages in {}".format(n_class_avg, args.isac_class_avgs_path))
	fullstack_img_id_list_of_isac_substack = []
	for class_avg_id in xrange(n_class_avg):
		fullstack_img_id_list_of_isac_class = []
		fullstack_img_id_list_of_isac_class = get_im(args.isac_class_avgs_path, class_avg_id).get_attr("members")
		fullstack_img_id_list_of_isac_class.sort()
		total_align2d_list_of_isac_class = []
		for fullstack_img_id in fullstack_img_id_list_of_isac_class:
			total_align2d = fullstack_total_align2d_list[fullstack_img_id]
			if total_align2d[idx_isac_align2d_mirror] == -1:
				ERROR("The member with original fullstack particle ID {} listed in ISAC class averages {} has the invalid 2D alignment parameters for ISAC unaccounted particle. Please check the consistency of input datasets. Worse yet, the input datasets might be corrupted. In this case, please consider to rerun ISAC.".format(fullstack_img_id, args.isac_class_avgs_path), subcommand_name) # action=1 - fatal error, exit
			assert (total_align2d[idx_isac_align2d_mirror] != -1) # all class member particles should be accounted!!!
			scale = 1.0 # because this 2D alignment parameters are scaled back!
			total_align2d_list_of_isac_class.append([total_align2d[idx_isac_align2d_alpha], total_align2d[idx_isac_align2d_tx], total_align2d[idx_isac_align2d_ty], total_align2d[idx_isac_align2d_mirror], scale])
		assert(len(total_align2d_list_of_isac_class) == len(fullstack_img_id_list_of_isac_class))
		align2d_avg_path = os.path.join(subdir_path, "%s_%03d.txt"%(align2d_avg_basename, class_avg_id))
		write_text_row(total_align2d_list_of_isac_class, align2d_avg_path)
		
		# Append class particle ID list to substack particle ID list
		fullstack_img_id_list_of_isac_substack += fullstack_img_id_list_of_isac_class
		
	# Sort the substack particle id list
	fullstack_img_id_list_of_isac_substack.sort()
	
	n_isac_substack_img = len(fullstack_img_id_list_of_isac_substack)
	print(" ")
	print_progress("Extracted {} ISAC class members from {}".format(n_isac_substack_img, args.isac_class_avgs_path))
	assert (n_accounted_img <= n_fullstack_img)
	if not n_isac_substack_img <= n_accounted_img:
		ERROR("Invalid number of ISAC class members {}. It must be smaller than or equal to the total number of ISAC accounted particles {}. The stack header might be corrupted. Please consider to rerun ISAC.".format(n_isac_substack_img, n_accounted_img), subcommand_name) # action=1 - fatal error, exit
	
	# Save the substack particle id list
	fullstack_img_id_path_of_isac_substack = os.path.join(args.output_directory, "{}_particle_id_list.txt".format(args.substack_basename))
	write_text_file(fullstack_img_id_list_of_isac_substack, fullstack_img_id_path_of_isac_substack)
	print(" ")
	print_progress("Saved original fullstack particle IDs of all members listed in ISAC class averages to {}.".format(fullstack_img_id_path_of_isac_substack))
	
	print(" ")
	print_progress("Converting 2D alignment parameters of all members listed in ISAC class averages to 2D alignment parameters header entry format...")
	isac_substack_total_header_align2d_list = []
	for isac_substack_img_id in xrange(n_isac_substack_img):
		fullstack_img_id = fullstack_img_id_list_of_isac_substack[isac_substack_img_id]
		# Get 2D alignment parameters associated with this particle and conver to 3D alignment parameters
		total_align2d = fullstack_total_align2d_list[fullstack_img_id]
		assert (total_align2d[idx_isac_align2d_mirror] != -1) # all class member particles should be accounted!!!
		# Register total_align2d to the list in xform.align2d format
		scale = 1.0 # because this 2D alignment parameters are scaled back!
		isac_substack_total_header_align2d_list.append([total_align2d[idx_isac_align2d_alpha], total_align2d[idx_isac_align2d_tx], total_align2d[idx_isac_align2d_ty], total_align2d[idx_isac_align2d_mirror], scale])
		assert (len(isac_substack_total_header_align2d_list[-1]) == n_idx_header_align2d)
	assert(len(isac_substack_total_header_align2d_list) == n_isac_substack_img)
	
	# Save the 2D alignment parameters of all members listed in ISAC class averages to file, using the xform.align2d header entry format.
	isac_substack_total_header_align2d_path = os.path.join(args.output_directory, "{}_header_align2d.txt".format(args.substack_basename))
	write_text_row(isac_substack_total_header_align2d_list, isac_substack_total_header_align2d_path)
	print(" ")
	print_progress("Saved the converted 2D alignment parameters to {}.".format(isac_substack_total_header_align2d_path))
	
	# Create virtual stack for ISAC substack
	assert (args.output_directory != "")
	print(" ")
	print_progress("Creating ISAC substack as a virtual stack...")
	virtual_bdb_substack_path = "bdb:{}#{}".format(args.output_directory, args.substack_basename)
	cmd_line = "e2bdb.py {} --makevstack={} --list={}".format(args.input_bdb_stack_path, virtual_bdb_substack_path, fullstack_img_id_path_of_isac_substack)
	status = cmdexecute(cmd_line)
	if status == 0: ERROR("\"{}\" execution failed. Exiting...".format(cmd_line), subcommand_name) # action=1 - fatal error, exit
	assert (EMUtil.get_image_count(virtual_bdb_substack_path) == n_isac_substack_img)
	
	# Import the total 2D alignment parameters to xform.align2d
	print(" ")
	print_progress("Importing the total 2D alignment parameters in the original scale to the header entry...")
	cmd_line = "sxheader.py {} --import={} --params=xform.align2d".format(virtual_bdb_substack_path, isac_substack_total_header_align2d_path)
	status = cmdexecute(cmd_line)
	if status == 0: ERROR("\"{}\" execution failed. Exiting...".format(cmd_line), subcommand_name) # action=1 - fatal error, exit
	
	# Transform xform.align2d to xform.projection
	print(" ")
	print_progress("Creating projection parameters header entry from imported 2D alignment parameters using 2D-to-3D transformation...")
	cmd_line = "sxparams_2D_to_3D.py {}".format(virtual_bdb_substack_path)
	status = cmdexecute(cmd_line)
	if status == 0: ERROR("\"{}\" execution failed. Exiting...".format(cmd_line), subcommand_name) # action=1 - fatal error, exit
	
	# Export projection parameters from xform.projection
	print(" ")
	print_progress("Exporting projection parameters from the header entry...")
	isac_substack_total_header_projection_path = os.path.join(args.output_directory, "{}_header_projection.txt".format(args.substack_basename))
	cmd_line = "sxheader.py {} --export={} --params=xform.projection".format(virtual_bdb_substack_path, isac_substack_total_header_projection_path)
	status = cmdexecute(cmd_line)
	if status == 0: ERROR("\"{}\" execution failed. Exiting...".format(cmd_line), subcommand_name) # action=1 - fatal error, exit
	
	# Print summary of processing
	print(" ")
	print_progress("Summary of processing...")
	print_progress("  Particles in fullstack  : %6d"%(n_fullstack_img)) 
	print_progress("  Accounted particles     : %6d"%(n_accounted_img))
	print_progress("  Defalut class averages  : %6d"%(n_default_class_avg)) 
	print_progress("  Provided class averages : %6d"%(n_class_avg)) 
	print_progress("  Extracted class members : %6d"%(n_isac_substack_img))
	print_progress("  ISAC substack size      : %6d"%(EMUtil.get_image_count(virtual_bdb_substack_path)))
	print(" ")

# ----------------------------------------------------------------------------------------
# TEST COMMAND
#
# cd /home/moriya/mrk_qa/mrktest_pipeline
#
# sxpipe.py resample_micrographs --help
# 
# rm -rf debug_mrkout_sxpipe_resample_micographs; time mpirun -np 10 sxpipe.py resample_micrographs 'mrkout_pipe00_inputs/micrographs/corrsum_dose_filtered/TcdA1-sialic-*_frames_sum.mrc' 'debug_mrkout_sxpipe_resample_micographs' --resample_ratio=0.5 --selection_list='debug_micrographs.txt' --check_consistency
#
# ----------------------------------------------------------------------------------------
def resample_micrographs(args):
	import glob
	import shutil
	from applications import MPI_start_end
	from inspect import currentframe, getframeinfo

	# ====================================================================================
	# Prepare processing
	# ====================================================================================
	# ------------------------------------------------------------------------------------
	# Check MPI execution
	# ------------------------------------------------------------------------------------
	if SXmpi_run.RUNNING_UNDER_MPI:
		from mpi import MPI_COMM_WORLD, mpi_barrier, mpi_reduce, MPI_INT, MPI_SUM

	# ------------------------------------------------------------------------------------
	# Set up SPHIRE global definitions
	# ------------------------------------------------------------------------------------
	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()
	
	command_script_basename = os.path.basename(sys.argv[0])
	program_name = "{} {}".format(command_script_basename, args.subcommand)
	
	# Change the name log file for error message
	original_logfilename = global_def.LOGFILE
	# global_def.LOGFILE = os.path.splitext(program_name)[0] + '_' + original_logfilename + '.txt'
	global_def.LOGFILE = os.path.splitext(command_script_basename)[0] + args.subcommand + '_' + original_logfilename + '.txt'
	
	# # To make the execution exit upon fatal error by ERROR in global_def.py
	# global_def.BATCH = True 
	
	# ------------------------------------------------------------------------------------
	# Check error conditions of arguments and options, then prepare variables for arguments
	# ------------------------------------------------------------------------------------
	mic_pattern = None
	root_out_dir = None
	# Not a real while, each "if" statement has the opportunity to use break when errors need to be reported
	error_status = None
	while True:
		mic_pattern = args.input_micrograph_pattern
		root_out_dir = args.output_directory
		
		# --------------------------------------------------------------------------------
		# Check error conditions of arguments
		# --------------------------------------------------------------------------------
		if error_status is None and mic_pattern is None:
			error_status = ("Missing required argument input_micrograph_pattern. Please run %s -h for help." % (program_name), getframeinfo(currentframe()))
			break

		if error_status is None and mic_pattern[:len("bdb:")].lower() == "bdb":
			error_status = ("BDB file can not be selected as input micrographs. Please convert the format, and restart the program. Run %s -h for help." % (program_name), getframeinfo(currentframe()))
			break
		
		if error_status is None and mic_pattern.find("*") == -1:
			error_status = ("Input micrograph file name pattern must contain wild card (*). Please check input_micrograph_pattern argument. Run %s -h for help." % (program_name), getframeinfo(currentframe()))
			break
		
		if error_status is None and root_out_dir is None:
			error_status = ("Missing required argument output_directory. Please run %s -h for help." % (program_name), getframeinfo(currentframe()))
			break

		if error_status is None and os.path.exists(root_out_dir):
			error_status = ("Output directory exists. Please change the name and restart the program.", getframeinfo(currentframe()))
			break
		
		# --------------------------------------------------------------------------------
		# Check error conditions of options
		# --------------------------------------------------------------------------------
		if error_status is None and args.resample_ratio is None:
			error_status = ("Missing required option --resample_ratio. Please run %s -h for help." % (program_name), getframeinfo(currentframe()))
			break
		
		if error_status is None and (args.resample_ratio <= 0.0 or args.resample_ratio >= 1.0):
			error_status = ("Invalid option value: --resample_ratio=%s. Please run %s -h for help." % (args.resample_ratio, program_name), getframeinfo(currentframe()))
			break
		
		if args.selection_list != None:
			if error_status is None and not os.path.exists(args.selection_list): 
				error_status = ("File specified by --selection_list option does not exists. Please check --selection_list option. Run %s -h for help." % (program_name), getframeinfo(currentframe()))
				break
		
		break
	if_error_then_all_processes_exit_program(error_status)
	assert (mic_pattern != None)
	assert (root_out_dir != None)
	
	# ------------------------------------------------------------------------------------
	# Prepare the variables for all sections
	# ------------------------------------------------------------------------------------
	# Micrograph basename pattern (directory path is removed from micrograph path pattern)
	mic_basename_pattern = os.path.basename(mic_pattern)
	
	# Global entry dictionary (all possible entries from all lists) for all mic id substring
	global_entry_dict = {} # mic id substring is the key
	subkey_input_mic_path = "Input Micrograph Path"
	subkey_selected_mic_basename = "Selected Micrograph Basename"
	
	# List keeps only id substrings of micrographs whose all necessary information are available
	valid_mic_id_substr_list = [] 
	
	# ====================================================================================
	# Obtain the list of micrograph id sustrings using a single CPU (i.e. main mpi process)
	# ====================================================================================
	# NOTE: Toshio Moriya 2016/10/24
	# The below is not a real while.  
	# It gives if-statements an opportunity to use break when errors need to be reported
	# However, more elegant way is to use 'raise' statement of exception mechanism...
	# 
	error_status = None
	while SXmpi_run.is_main_proc():
		# --------------------------------------------------------------------------------
		# Prepare variables for this section
		# --------------------------------------------------------------------------------
		# Prefix and suffix of micrograph basename pattern 
		# to find the head/tail indices of micrograph id substring
		mic_basename_tokens = mic_basename_pattern.split('*')
		assert (len(mic_basename_tokens) == 2)
		# Find head index of micrograph id substring
		mic_id_substr_head_idx = len(mic_basename_tokens[0])
		
		# --------------------------------------------------------------------------------
		# Register micrograph id substrings found in the input directory (specified by micrograph path pattern)
		# to the global entry dictionary
		# --------------------------------------------------------------------------------
		# Generate the list of micrograph paths in the input directory
		print(" ")
		print("Checking the input directory...")
		input_mic_path_list = glob.glob(mic_pattern)
		# Check error condition of input micrograph file path list
		print("Found %d microgarphs in %s." % (len(input_mic_path_list), os.path.dirname(mic_pattern)))
		if error_status is None and len(input_mic_path_list) == 0:
			error_status = ("No micrograph files are found in the directory specified by micrograph path pattern (%s). Please check input_micrograph_pattern argument. Run %s -h for help." % (os.path.dirname(mic_pattern), program_name), getframeinfo(currentframe()))
			break
		assert (len(input_mic_path_list) > 0)
		
		# Register micrograph id substrings to the global entry dictionary
		for input_mic_path in input_mic_path_list:
			# Find tail index of micrograph id substring and extract the substring from the micrograph name
			input_mic_basename = os.path.basename(input_mic_path)
			mic_id_substr_tail_idx = input_mic_basename.index(mic_basename_tokens[1])
			mic_id_substr = input_mic_basename[mic_id_substr_head_idx:mic_id_substr_tail_idx]
			assert (input_mic_path == mic_pattern.replace("*", mic_id_substr))
			if not mic_id_substr in global_entry_dict:
				# print("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from input_mic_path_list " % (mic_id_substr))
				global_entry_dict[mic_id_substr] = {}
			assert (mic_id_substr in global_entry_dict)
			global_entry_dict[mic_id_substr][subkey_input_mic_path] = input_mic_path
		assert (len(global_entry_dict) > 0)
		
		# --------------------------------------------------------------------------------
		# Register micrograph id substrings found in the selection list
		# to the global entry dictionary
		# --------------------------------------------------------------------------------
		# Generate the list of selected micrograph paths in the selection file
		selected_mic_path_list = []
		# Generate micrograph lists according to the execution mode
		if args.selection_list == None:
			print(" ")
			print("----- Running with All Micrographs Mode -----")
			# Treat all micrographs in the input directory as selected ones
			selected_mic_path_list = input_mic_path_list
		else:
			assert (args.selection_list != None)
			if os.path.splitext(args.selection_list)[1] == ".txt":
				print(" ")
				print("----- Running with Selected Micrographs Mode -----")
				print(" ")
				print("Checking the selection list...")
				assert (os.path.exists(args.selection_list))
				selected_mic_path_list = read_text_file(args.selection_list)
				
				# Check error condition of micrograph entry lists
				print("Found %d microgarph entries in %s." % (len(selected_mic_path_list), args.selection_list))
				if error_status is None and len(selected_mic_path_list) == 0:
					error_status = ("No micrograph entries are found in the selection list file. Please check selection_list option. Run %s -h for help." % (program_name), getframeinfo(currentframe()))
					break
				assert (len(selected_mic_path_list) > 1)
				if error_status is None and not isinstance(selected_mic_path_list[0], basestring):
					error_status = ("Invalid format of the selection list file. The first column must contain micrograph paths in string type. Please check selection_list option. Run %s -h for help." % (program_name), getframeinfo(currentframe()))
					break
			else:
				print(" ")
				print("----- Running with Single Micrograph Mode -----")
				print(" ")
				print("Processing a single micorgprah: %s..." % (args.selection_list))
				selected_mic_path_list = [args.selection_list]
			assert (len(selected_mic_path_list) > 0)
			
			selected_mic_directory = os.path.dirname(selected_mic_path_list[0])
			if selected_mic_directory != "":
				print("    NOTE: Program disregards the directory paths in the selection list (%s)." % (selected_mic_directory))
			
		assert (len(selected_mic_path_list) > 0)
		
		# Register micrograph id substrings to the global entry dictionary
		for selected_mic_path in selected_mic_path_list:
			# Find tail index of micrograph id substring and extract the substring from the micrograph name
			selected_mic_basename = os.path.basename(selected_mic_path)
			mic_id_substr_tail_idx = selected_mic_basename.index(mic_basename_tokens[1])
			mic_id_substr = selected_mic_basename[mic_id_substr_head_idx:mic_id_substr_tail_idx]
			if error_status is None and selected_mic_basename != mic_basename_pattern.replace("*", mic_id_substr):
				error_status = ("A micrograph name (%s) in the input directory (%s) does not match with input micrograph basename pattern (%s) (The wild card replacement with \'%s\' resulted in \'%s\'). Please correct input micrograph path pattern. Run %s -h for help." % (selected_mic_basename, os.path.dirname(mic_pattern), mic_basename_pattern, mic_id_substr, mic_basename_pattern.replace("*", mic_id_substr), program_name), getframeinfo(currentframe()))
				break
			if not mic_id_substr in global_entry_dict:
				# print("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from selected_mic_path_list " % (mic_id_substr))
				global_entry_dict[mic_id_substr] = {}
			assert (mic_id_substr in global_entry_dict)
			global_entry_dict[mic_id_substr][subkey_selected_mic_basename] = selected_mic_basename
		assert (len(global_entry_dict) > 0)
		
		del selected_mic_path_list # Do not need this anymore
		del input_mic_path_list # Do not need this anymore
		
		# --------------------------------------------------------------------------------
		# Clean up variables related to registration to the global entry dictionary
		# --------------------------------------------------------------------------------
		del mic_basename_tokens
		del mic_id_substr_head_idx
		
		# --------------------------------------------------------------------------------
		# Create the list containing only valid micrograph id substrings
		# --------------------------------------------------------------------------------
		# Prepare lists to keep track of invalid (rejected) micrographs 
		no_input_mic_id_substr_list = []
		
		print(" ")
		print("Checking consistency of the provided dataset ...")
		
		# Loop over substring id list
		for mic_id_substr in global_entry_dict:
			mic_id_entry = global_entry_dict[mic_id_substr]
			
			warinnig_messages = []
			# selected micrograph basename must have been registed always .
			if subkey_selected_mic_basename in mic_id_entry: 
				# Check if associated input micrograph exists
				if not subkey_input_mic_path in mic_id_entry:
					input_mic_path = mic_pattern.replace("*", mic_id_substr)
					warinnig_messages.append("    associated input micrograph %s." % (input_mic_path))
					no_input_mic_id_substr_list.append(mic_id_substr)
				
				if len(warinnig_messages) > 0:
					print("WARNING!!! Micrograph ID %s does not have:" % (mic_id_substr))
					for warinnig_message in warinnig_messages:
						print(warinnig_message)
					print("    Ignores this as an invalid entry.")
				else:
					# print("MRK_DEBUG: adding mic_id_substr := ", mic_id_substr)
					valid_mic_id_substr_list.append(mic_id_substr)
			# else:
			# 	assert (not subkey_selected_mic_basename in mic_id_entry)
			# 	# This entry is not in the selection list. Do nothing
			
		# Check the input dataset consistency and save the result to a text file, if necessary.
		if args.check_consistency:
			# Create output directory
			assert (not os.path.exists(root_out_dir))
			os.mkdir(root_out_dir)
			
			# Open the consistency check file
			mic_consistency_check_info_path = os.path.join(root_out_dir, "mic_consistency_check_info_%s.txt"%(get_time_stamp_suffix()))
			print(" ")
			print("Generating consistency report of the provided dataset in %s..."%(mic_consistency_check_info_path))
			mic_consistency_check_info_file = open(mic_consistency_check_info_path, "w")
			mic_consistency_check_info_file.write("# The consistency information about micrograph IDs that might have problmes with consistency among the provided dataset.\n")
			mic_consistency_check_info_file.write("# \n")
			# Loop over substring id list
			for mic_id_substr in global_entry_dict:
				mic_id_entry = global_entry_dict[mic_id_substr]
				
				consistency_messages = []
				# Check if associated input micrograph path exists
				if not subkey_input_mic_path in mic_id_entry:
					input_mic_path = mic_pattern.replace("*", mic_id_substr)
					consistency_messages.append("    associated input micrograph %s." % (input_mic_path))
				
				# Check if associated selected micrograph basename exists
				if not subkey_selected_mic_basename in mic_id_entry:
					input_mic_path = mic_pattern.replace("*", mic_id_substr)
					consistency_messages.append("    associated selected micrograph %s." % (input_mic_path))
				
				if len(consistency_messages) > 0:
					mic_consistency_check_info_file.write("Micrograph ID %s might have problems with consistency among the provided dataset:\n"%(mic_id_substr))
					for consistency_message in consistency_messages:
						mic_consistency_check_info_file.write(consistency_message)
						mic_consistency_check_info_file.write("\n")
			
			# Close the consistency check file, if necessary
			mic_consistency_check_info_file.flush()
			mic_consistency_check_info_file.close()
			
		# Since mic_id_substr is once stored as the key of global_entry_dict and extracted with the key order
		# we need sort the valid_mic_id_substr_list here
		# print("MRK_DEBUG: before sort, valid_mic_id_substr_list := ", valid_mic_id_substr_list)
		valid_mic_id_substr_list.sort()
		# print("MRK_DEBUG: after sort, valid_mic_id_substr_list := ", valid_mic_id_substr_list)
		
		# --------------------------------------------------------------------------------
		# Print out the summary of input consistency
		# --------------------------------------------------------------------------------
		print(" ")
		print("Summary of dataset consistency check...")
		print("Detected                           : %6d" % (len(global_entry_dict)))
		print("Valid                              : %6d" % (len(valid_mic_id_substr_list)))
			
		# --------------------------------------------------------------------------------
		# Clean up variables related to tracking of invalid (rejected) micrographs 
		# --------------------------------------------------------------------------------
		del no_input_mic_id_substr_list
		
		# --------------------------------------------------------------------------------
		# Check MPI error condition
		# --------------------------------------------------------------------------------
		if error_status is None and len(valid_mic_id_substr_list) < SXmpi_run.n_mpi_procs:
			error_status = ("Number of MPI processes (%d) supplied by --np in mpirun cannot be greater than %d (number of valid micrographs that satisfy all criteria to be processed)." % (SXmpi_run.n_mpi_procs, len(valid_mic_id_substr_list)), getframeinfo(currentframe()))
			break
		
		break
	# 
	# NOTE: Toshio Moriya 2016/10/24
	# The following function takes care of the case when an if-statement uses break for occurence of an error.
	# However, more elegant way is to use 'exception' statement of exception mechanism...
	# 
	if_error_then_all_processes_exit_program(error_status)
	
	# ====================================================================================
	# Obtain the list of micrograph id sustrings
	# ====================================================================================
	# --------------------------------------------------------------------------------
	# Prepare variables for this section
	# --------------------------------------------------------------------------------
	# Prepare variables related to options
	resample_ratio = args.resample_ratio
	
	# Micrograph baseroot pattern (extension are removed from micrograph basename pattern)
	# for substack file names
	mic_baseroot_pattern = os.path.splitext(mic_basename_pattern)[0]  
	
	# Prepare the counters for the global summary of micrographs
	n_mic_process = 0
	
	# keep a copy of the root output directory where the final bdb will be created
	unsliced_valid_serial_id_list = valid_mic_id_substr_list
	if SXmpi_run.RUNNING_UNDER_MPI:
		mpi_barrier(MPI_COMM_WORLD)
		# All mpi processes should know global entry directory and valid micrograph id substring list
		global_entry_dict = wrap_mpi_bcast(global_entry_dict, SXmpi_run.main_mpi_proc)
		valid_mic_id_substr_list = wrap_mpi_bcast(valid_mic_id_substr_list, SXmpi_run.main_mpi_proc)
		
		# Slice the list of valid micrograph id substrings for this mpi process
		mic_start, mic_end = MPI_start_end(len(valid_mic_id_substr_list), SXmpi_run.n_mpi_procs, SXmpi_run.my_mpi_proc_id)
		valid_mic_id_substr_list = valid_mic_id_substr_list[mic_start:mic_end]
		
	if SXmpi_run.is_main_proc():
		print(" ")
		print("Micrographs processed by main process (including percent of progress):")
		progress_percent_step = len(valid_mic_id_substr_list)/100.0 # the number of micrograms for main node divided by 100
		
		# Create output directory
		# 
		# NOTE: Toshio Moriya 2017/12/12
		# This might not be necessary since particle_img.write_image() will automatically create all directory tree necessary to save the file.
		# However, it is side-effect of the function, so we will explicitly make root output directory here.
		# 
		if not os.path.exists(root_out_dir):
			os.mkdir(root_out_dir)

	# All node should wait for main node to create root output directory
	if SXmpi_run.RUNNING_UNDER_MPI:
		mpi_barrier(MPI_COMM_WORLD) # all MPI processes should wait until the directory is created by main process
		# 
		# NOTE: 2017/12/12 Toshio Moriya
		# To walk-around synchronisation problem between all MPI nodes and a file server,
		# we use exception to assert the existence of directory.
		# 
		try: 
			os.mkdir(root_out_dir)
			assert False, "Unreachable code..."
		except OSError as err:
			pass
	else: 
		assert os.path.exists(root_out_dir)
	
	# ------------------------------------------------------------------------------------
	# Starting main parallel execution
	# ------------------------------------------------------------------------------------
	for mic_id_substr_idx, mic_id_substr in enumerate(valid_mic_id_substr_list):
		
		# --------------------------------------------------------------------------------
		# Print out progress if necessary
		# --------------------------------------------------------------------------------
		mic_basename = global_entry_dict[mic_id_substr][subkey_selected_mic_basename]
		assert (mic_basename == mic_basename_pattern.replace("*", mic_id_substr))
		if SXmpi_run.is_main_proc():
			print("%s ---> % 2.2f%%" % (mic_basename, mic_id_substr_idx / progress_percent_step))
		
		# --------------------------------------------------------------------------------
		# Read micrograph
		# --------------------------------------------------------------------------------
		mic_path = global_entry_dict[mic_id_substr][subkey_input_mic_path]
		assert (mic_path == mic_pattern.replace("*", mic_id_substr))
		try:
			mic_img = get_im(mic_path)
		except:
			print("Failed to read the associate micrograph %s for %s. The file might be corrupted. Skipping..." % (mic_path, mic_basename))
			continue
		
		# --------------------------------------------------------------------------------
		# Resample micrograph, map coordinates, and window segments from resampled micrograph using new coordinates
		# after resampling by resample_ratio, resampled pixel size = src_pixel_size/resample_ratio
		# --------------------------------------------------------------------------------
		# NOTE: 2015/04/13 Toshio Moriya
		# resample() efficiently takes care of the case resample_ratio = 1.0 but
		# it does not set apix_*. Even though it sets apix_* when resample_ratio < 1.0...
		mic_img = resample(mic_img, resample_ratio)
		
		# --------------------------------------------------------------------------------
		# Generate the output file path of particle stack for this mpi process
		# --------------------------------------------------------------------------------
		mic_basename = mic_basename_pattern.replace("*", mic_id_substr)
		output_mic_path = os.path.join(root_out_dir, mic_basename)
		mic_img.write_image(output_mic_path) 
		
		# --------------------------------------------------------------------------------
		# Prepare coordinates loop variables
		# --------------------------------------------------------------------------------
		# Update the counters for the global summary of micrographs
		n_mic_process += 1
	
	# ------------------------------------------------------------------------------------
	# Print out summary of processing
	# ------------------------------------------------------------------------------------
	if SXmpi_run.RUNNING_UNDER_MPI:
		n_mic_process = mpi_reduce(n_mic_process, 1, MPI_INT, MPI_SUM, SXmpi_run.main_mpi_proc, MPI_COMM_WORLD)
	
	# Print out the summary of all micrographs
	if SXmpi_run.is_main_proc():
		print(" ")
		print("Summary of process...")
		print("Valid                              : %6d" % (len(unsliced_valid_serial_id_list)))
		print("Processed                          : %6d" % (n_mic_process))
	
	if SXmpi_run.RUNNING_UNDER_MPI:
		mpi_barrier(MPI_COMM_WORLD)
	
	if SXmpi_run.is_main_proc():
		print(" ")
		print("DONE!!!")
		print(" ")
	
	# ====================================================================================
	# Clean up
	# ====================================================================================
	# ------------------------------------------------------------------------------------
	# Reset SPHIRE global definitions
	# ------------------------------------------------------------------------------------
	global_def.LOGFILE = original_logfilename
	
	sys.stdout.flush()


# ----------------------------------------------------------------------------------------
# TEST COMMAND
# cd /home/moriya/SPHIRE-demo/Precalculated-Results
# 
# ls -l CorrectedSums/MRK_DISCARDED
# rm -r CorrectedSums/MRK_DISCARDED
# 
# ls -l CorrectedSums/MRK_DISCARDED_DUPLICATED
# rm -r CorrectedSums/MRK_DISCARDED_DUPLICATED
# 
# sxpipe.py organize_micrographs 'CorrectedSums/corrsum/TcdA1-*_frames_sum.mrc' 'CTFest/Tutorial_micrographs_discard.txt' 'CorrectedSums/MRK_DISCARDED' --check_consistency 2>&1 | tee sxpipe_organize_micrographs01.log
# 
# sxpipe.py organize_micrographs 'CorrectedSums/corrsum/TcdA1-*_frames_sum.mrc' 'CTFest/Tutorial_micrographs_discard.txt' 'CorrectedSums/MRK_DISCARDED' --reverse --check_consistency 2>&1 | tee sxpipe_organize_micrographs02.log
# 
# sxpipe.py organize_micrographs 'CorrectedSums/corrsum/TcdA1-*_frames_sum.mrc' 'CTFest/Tutorial_micrographs_discard.txt' 'CorrectedSums/MRK_DISCARDED' --check_consistency 2>&1 | tee sxpipe_organize_micrographs03.log
# 
# cp -r CorrectedSums/MRK_DISCARDED CorrectedSums/MRK_DISCARDED_DUPLICATED
# ls -l CorrectedSums/MRK_DISCARDED_DUPLICATED
# 
# sxpipe.py organize_micrographs 'CorrectedSums/corrsum/TcdA1-*_frames_sum.mrc' 'CTFest/Tutorial_micrographs_discard.txt' 'CorrectedSums/MRK_DISCARDED' --reverse --check_consistency 2>&1 | tee sxpipe_organize_micrographs04.log
# 
# sxpipe.py organize_micrographs 'CorrectedSums/corrsum/TcdA1-*_frames_sum.mrc' 'CTFest/Tutorial_micrographs_discard.txt' 'CorrectedSums/MRK_DISCARDED_DUPLICATED' --check_consistency 2>&1 | tee sxpipe_organize_micrographs05.log
# 
# sxpipe.py organize_micrographs 'CorrectedSums/corrsum/TcdA1-*_frames_sum.mrc' 'CTFest/Tutorial_micrographs_discard.txt' 'CorrectedSums/MRK_DISCARDED_DUPLICATED' --reverse  --check_consistency 2>&1 | tee sxpipe_organize_micrographs06.log
# 
# ----------------------------------------------------------------------------------------
def organize_micrographs(args):
	import glob
	import shutil
	from utilities import read_text_file
	
	# Check MPI execution
	if SXmpi_run.n_mpi_procs > 1:
		assert (SXmpi_run.RUNNING_UNDER_MPI)
		error_status = ("The {} subcommand supports only a single process.".format(subcommand_name), getframeinfo(currentframe()))
		if_error_then_all_processes_exit_program(error_status)
	
	# To make the execution exit upon fatal error by ERROR in global_def.py
	global_def.BATCH = True 
	
	# ------------------------------------------------------------------------------------
	# Prepare the variables for all sections
	# ------------------------------------------------------------------------------------
	# Use short names for arguments and options
	src_mic_pattern = args.source_micrograph_pattern
	select_list_path = args.selection_list
	dst_dir = args.destination_directory

	# ------------------------------------------------------------------------------------
	# Check error conditions
	# ------------------------------------------------------------------------------------
	# subcommand_name = "organize_micrographs"
	command_script_basename = os.path.basename(sys.argv[0])
	subcommand_name = "{} {}".format(command_script_basename, args.subcommand)
	
	if src_mic_pattern.find("*") == -1:
		ERROR("The source micrograph path pattern must contain wild card (*). Please correct source_micrograph_pattern argument and restart the program.", subcommand_name) # action=1 - fatal error, exit
	
	if os.path.splitext(select_list_path)[1] != ".txt":
		ERROR("The extension of source micrograph selecting list file must \'.txt\'. Please choose a correct file path or change the file extension, then restart the program.", subcommand_name) # action=1 - fatal error, exit
	
	if not os.path.exists(select_list_path):
		ERROR("The micrograph selecting list file does not exist. Please choose a correct file path and restart the program.", subcommand_name) # action=1 - fatal error, exit
	assert (os.path.exists(select_list_path))
	
	# ------------------------------------------------------------------------------------
	# Define operation mode information
	# ------------------------------------------------------------------------------------
	# Micrograph basename pattern (directory path is removed from micrograph path pattern)
	mic_basename_pattern = os.path.basename(src_mic_pattern)
	src_dir = os.path.dirname(src_mic_pattern)
	record_dir = dst_dir # always use the original output directory for recording generated information
	
	# Swap input directory and output directory if necessary
	if not args.reverse:
		print(" ")
		print_progress("Running with Normal Operation Mode... ")
	else:
		assert (args.reverse)
		print(" ")
		print_progress("Running with Reverse Operation Mode... ")
		dst_dir = src_dir
		src_dir = record_dir
		src_mic_pattern = os.path.join(src_dir, mic_basename_pattern)
	
	print_progress("Source micrograph basename pattern : %s"%(src_mic_pattern))
	print_progress("Source directory                   : %s"%(src_dir))
	print_progress("Destination directory              : %s"%(dst_dir))
	print_progress("Recording directory                : %s"%(record_dir))
	print(" ")
	

	# --------------------------------------------------------------------------------
	# Prepare variables
	# --------------------------------------------------------------------------------
	# Define indices of selection list parameters
	i_enum = -1
	i_enum += 1; idx_mic_list_mic_path   = i_enum # The name or path of micrographs
	i_enum += 1; n_idx_mic_list          = i_enum

	# Global entry dictionary (all possible entries from all lists) for all mic id substring
	global_entry_dict = {} # mic id substring is the key
	subkey_src_mic_path = "Source Micrograph Path"
	subkey_dst_mic_path = "Destination Micrograph Path"
	subkey_select_mic_basename = "Selected Micrograph Basename"
	
	# List keeps only id substrings of micrographs whose all necessary information are available
	valid_mic_id_substr_list = [] 
	
	# Prefix and suffix of micrograph basename pattern 
	# to find the head/tail indices of micrograph id substring
	mic_basename_tokens = mic_basename_pattern.split("*")
	assert (len(mic_basename_tokens) == 2)
	# Find head index of micrograph id substring
	mic_id_substr_head_idx = len(mic_basename_tokens[0])

	# Set up output directory 
	dst_mic_pattern = None
	if os.path.exists(dst_dir):
		print(" ")
		print_progress("The destination directory (%s) already exists. "%(dst_dir))
		dst_mic_pattern = os.path.join(dst_dir, mic_basename_pattern)
	
	# --------------------------------------------------------------------------------
	# Register micrograph id substrings found in source directory (specified by source micrograph path pattern)
	# and associated source micrograph path to the global entry dictionary
	# --------------------------------------------------------------------------------
	# Generate the list of micrograph paths in the source directory
	print(" ")
	print_progress("Checking the source directory...")
	src_mic_path_list = glob.glob(src_mic_pattern)
	# Check error condition of source micrograph file path list
	print_progress("Found %d microgarphs in %s."%(len(src_mic_path_list), src_dir))
	if len(src_mic_path_list) == 0:
		ERROR("No micrograph files are found in the directory specified by the micrograph path pattern (%s). Please check source_micrograph_pattern argument and restart the program."%(src_dir), subcommand_name) # action=1 - fatal error, exit
	assert (len(src_mic_path_list) > 0)
	
	# Register micrograph id substrings to the global entry dictionary
	for src_mic_path in src_mic_path_list:
		# Find tail index of micrograph id substring and extract the substring from the micrograph name
		src_mic_basename = os.path.basename(src_mic_path)
		mic_id_substr_tail_idx = src_mic_basename.index(mic_basename_tokens[1])
		mic_id_substr = src_mic_basename[mic_id_substr_head_idx:mic_id_substr_tail_idx]
		assert (src_mic_path == src_mic_pattern.replace("*", mic_id_substr))
		if not mic_id_substr in global_entry_dict:
			# print("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from src_mic_path_list "%(mic_id_substr))
			global_entry_dict[mic_id_substr] = {}
		assert (mic_id_substr in global_entry_dict)
		global_entry_dict[mic_id_substr][subkey_src_mic_path] = src_mic_path
	assert (len(global_entry_dict) > 0)
	
	# Clean up variables which won't be used anymore
	del src_mic_path_list

	# --------------------------------------------------------------------------------
	# Register micrograph id substrings found in destination directory if any
	# and associated source micrograph path to the global entry dictionary
	# --------------------------------------------------------------------------------
	if dst_mic_pattern is not None:
		assert (os.path.exists(dst_dir))
		dst_mic_pattern = os.path.join(dst_dir, mic_basename_pattern)
		# Generate the list of micrograph paths in the output directory
		print(" ")
		print_progress("Checking the destination directory...")
		dst_mic_path_list = glob.glob(dst_mic_pattern)
		# Check error condition of destination micrograph file path list
		print_progress("Found %d microgarphs in %s."%(len(dst_mic_path_list), dst_dir))
		
		# Register micrograph id substrings to the global entry dictionary
		for dst_mic_path in dst_mic_path_list:
			# Find tail index of micrograph id substring and extract the substring from the micrograph name
			dst_mic_basename = os.path.basename(dst_mic_path)
			mic_id_substr_tail_idx = dst_mic_basename.index(mic_basename_tokens[1])
			mic_id_substr = dst_mic_basename[mic_id_substr_head_idx:mic_id_substr_tail_idx]
			assert (dst_mic_path == dst_mic_pattern.replace("*", mic_id_substr))
			if not mic_id_substr in global_entry_dict:
				# print("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from dst_mic_path_list "%(mic_id_substr))
				global_entry_dict[mic_id_substr] = {}
			assert (mic_id_substr in global_entry_dict)
			global_entry_dict[mic_id_substr][subkey_dst_mic_path] = dst_mic_path
		assert (len(global_entry_dict) > 0)
	
		# Clean up variables which won't be used anymore
		del dst_mic_path_list

	# --------------------------------------------------------------------------------
	# Register micrograph id substrings found in the selection list
	# and associated micrograph basename to the global entry dictionary
	# --------------------------------------------------------------------------------
	# Generate the list of select micrograph paths in the selection file
	select_mic_path_list = []
	# Generate micrograph lists according to the execution mode
	print(" ")
	print_progress("Checking the selection list...")
	select_mic_path_list = read_text_file(select_list_path)
	
	# Check error condition of micrograph entry lists
	print_progress("Found %d microgarph entries in %s."%(len(select_mic_path_list), select_list_path))
	if len(select_mic_path_list) == 0:
		ERROR("No micrograph entries are found in the selection list file (%s). Please correct selection_list option and restart the program."%(select_list_path), subcommand_name) # action=1 - fatal error, exit
	assert (len(select_mic_path_list) > 0)
	
	select_mic_dir = os.path.dirname(select_mic_path_list[0])
	if select_mic_dir != "":
		print_progress("    NOTE: Program disregards the directory paths in the source selection list (%s)."%(select_mic_dir))

	# Register micrograph id substrings to the global entry dictionary
	for select_mic_path in select_mic_path_list:
		# Find tail index of micrograph id substring and extract the substring from the micrograph name
		select_mic_basename = os.path.basename(select_mic_path)
		mic_id_substr_tail_idx = select_mic_basename.index(mic_basename_tokens[1])
		mic_id_substr = select_mic_basename[mic_id_substr_head_idx:mic_id_substr_tail_idx]
		assert (select_mic_basename == mic_basename_pattern.replace("*", mic_id_substr))
		if not mic_id_substr in global_entry_dict:
			# print("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from select_mic_path_list "%(mic_id_substr))
			global_entry_dict[mic_id_substr] = {}
		assert (mic_id_substr in global_entry_dict)
		global_entry_dict[mic_id_substr][subkey_select_mic_basename] = select_mic_basename
	assert (len(global_entry_dict) > 0)
	
	# Clean up variables which won't be used anymore
	del select_mic_path_list
	
	# --------------------------------------------------------------------------------
	# Clean up variables related to registration to the global entry dictionary
	# --------------------------------------------------------------------------------
	del mic_basename_tokens
	del mic_id_substr_head_idx

	# --------------------------------------------------------------------------------
	# Create the list containing only valid micrograph id substrings
	# --------------------------------------------------------------------------------
	print(" ")
	print_progress("Checking consistency of the provided dataset ...")

	if dst_mic_pattern is None:
		assert (not os.path.exists(dst_dir))
		# Prepare lists to keep track of invalid (rejected) micrographs
		no_src_mic_id_substr_list = []
		
		# Loop over substring id list
		for mic_id_substr in global_entry_dict:
			mic_id_entry = global_entry_dict[mic_id_substr]
		
			warinnig_messages = []
			# selected micrograph basename must have been registed always .
			if subkey_select_mic_basename in mic_id_entry: 
				# Check if associated input micrograph exists
				if not subkey_src_mic_path in mic_id_entry:
					mic_basename = mic_basename_pattern.replace("*", mic_id_substr)
					warinnig_messages.append("    associated micrograph (%s) does not exist in the source directory (%s)."%(mic_basename, src_dir))
					no_src_mic_id_substr_list.append(mic_id_substr)
				
				if len(warinnig_messages) > 0:
					print_progress("WARNING!!! Micrograph ID %s has problems with consistency among the provided dataset:"%(mic_id_substr))
					for warinnig_message in warinnig_messages:
						print_progress(warinnig_message)
					print_progress("    Ignores this as an invalid entry.")
				else:
					# print("MRK_DEBUG: adding mic_id_substr := ", mic_id_substr)
					valid_mic_id_substr_list.append(mic_id_substr)
			# else:
			# 	assert (not subkey_select_mic_basename in mic_id_entry)
			# 	# This entry is not in the selection list. Do nothing
			
		# Check the input dataset consistency and save the result to a text file, if necessary.
		if args.check_consistency:
			# Create destination directory
			assert (not os.path.exists(dst_dir))
			os.mkdir(dst_dir)
			assert (os.path.exists(dst_dir))
			
			# Open the consistency check file
			mic_consistency_check_info_path = os.path.join(record_dir, "mic_consistency_check_info_%s.txt"%(get_time_stamp_suffix()))
			print(" ")
			print_progress("Generating consistency report of the provided dataset in %s..."%(mic_consistency_check_info_path))
			mic_consistency_check_info_file = open(mic_consistency_check_info_path, "w")
			mic_consistency_check_info_file.write("# The consistency information about micrograph IDs that might have problmes with consistency among the provided dataset.\n")
			mic_consistency_check_info_file.write("# \n")
			# Loop over substring id list
			for mic_id_substr in global_entry_dict:
				mic_id_entry = global_entry_dict[mic_id_substr]
			
				consistency_messages = []
				# Check if associated micrograph path exists in source directory
				if not subkey_src_mic_path in mic_id_entry:
					mic_basename = mic_basename_pattern.replace("*", mic_id_substr)
					consistency_messages.append("    associated micrograph (%s) does not exist in the source directory (%s)."%(mic_basename, src_dir))
			
				# Check if associated micrograph basename exists in selection list
				if not subkey_select_mic_basename in mic_id_entry:
					mic_basename = mic_basename_pattern.replace("*", mic_id_substr)
					consistency_messages.append("    associated micrograph (%s) is not in the selection list (%s)."%(mic_basename, select_list_path))
			
				if len(consistency_messages) > 0:
					mic_consistency_check_info_file.write("Micrograph ID %s might have problems with consistency among the provided dataset:\n"%(mic_id_substr))
					for consistency_message in consistency_messages:
						mic_consistency_check_info_file.write(consistency_message)
						mic_consistency_check_info_file.write("\n")
		
			# Close the consistency check file, if necessary
			mic_consistency_check_info_file.flush()
			mic_consistency_check_info_file.close()
		
		# Since mic_id_substr is once stored as the key of global_entry_dict and extracted with the key order
		# we need sort the valid_mic_id_substr_list here
		# print("MRK_DEBUG: before sort, valid_mic_id_substr_list := ", valid_mic_id_substr_list)
		valid_mic_id_substr_list.sort()
		# print("MRK_DEBUG: after sort, valid_mic_id_substr_list := ", valid_mic_id_substr_list)
	
		# --------------------------------------------------------------------------------
		# Print out the summary of input consistency
		# --------------------------------------------------------------------------------
		print(" ")
		print_progress("Summary of consistency check for provided dataset...")
		print_progress("Detected                           : %6d"%(len(global_entry_dict)))
		print_progress("Valid                              : %6d"%(len(valid_mic_id_substr_list)))
		print_progress("Rejected by no source micrographs  : %6d"%(len(no_src_mic_id_substr_list)))
		print(" ")
		
		# --------------------------------------------------------------------------------
		# Clean up variables related to tracking of invalid (rejected) micrographs 
		# --------------------------------------------------------------------------------
		del no_src_mic_id_substr_list

	else:
		assert (dst_mic_pattern is not None)
		assert (os.path.exists(dst_dir))
		# Prepare lists to keep track of invalid (rejected) micrographs
		no_mic_in_both_dirs_id_substr_list = []
		already_in_dst_dir_mic_id_substr_list = []
		duplicated_in_dst_dir_mic_id_substr_list = []

		# Loop over substring id list
		for mic_id_substr in global_entry_dict:
			mic_id_entry = global_entry_dict[mic_id_substr]
		
			warinnig_messages = []
			# selected micrograph basename must have been registed always .
			if subkey_select_mic_basename in mic_id_entry: 
				# Check if associated input micrograph exists
				if not subkey_src_mic_path in mic_id_entry:
					mic_basename = mic_basename_pattern.replace("*", mic_id_substr)
					if not subkey_dst_mic_path in mic_id_entry:
						warinnig_messages.append("    associated micrograph (%s) does not exist neither in the source directory (%s) nor in the destination directory (%s)."%(mic_basename, src_dir, dst_dir))
						no_mic_in_both_dirs_id_substr_list.append(mic_id_substr)
					else:
						assert (subkey_dst_mic_path in mic_id_entry)
						warinnig_messages.append("    associated micrograph (%s) exists only in the destination directory (%s), but not in the source directory (%s)."%(mic_basename, dst_dir, src_dir))
						already_in_dst_dir_mic_id_substr_list.append(mic_id_substr)
				else: 
					assert (subkey_src_mic_path in mic_id_entry)
					if subkey_dst_mic_path in mic_id_entry:
						mic_basename = mic_basename_pattern.replace("*", mic_id_substr)
						warinnig_messages.append("    associated micrograph (%s) exist both in the source directory (%s) and in the destination directory (%s)."%(mic_basename, src_dir, dst_dir))
						duplicated_in_dst_dir_mic_id_substr_list.append(mic_id_substr)
					# else:
					#	# This should most typical case!
					#	assert (not subkey_dst_mic_path in mic_id_entry)
				if len(warinnig_messages) > 0:
					print_progress("WARNING!!! Micrograph ID %s has problems with consistency among the provided dataset:"%(mic_id_substr))
					for warinnig_message in warinnig_messages:
						print_progress(warinnig_message)
					print_progress("    Ignores this as an invalid entry.")
				else:
					# print("MRK_DEBUG: adding mic_id_substr := ", mic_id_substr)
					valid_mic_id_substr_list.append(mic_id_substr)
			# else:
			# 	assert (not subkey_select_mic_basename in mic_id_entry)
			# 	# This entry is not in the selection list. Do nothing
			
		# Check the input dataset consistency and save the result to a text file, if necessary.
		if args.check_consistency:
			assert (os.path.exists(dst_dir))
			
			# Open the consistency check file
			mic_consistency_check_info_path = os.path.join(record_dir, "mic_consistency_check_info_%s.txt"%(get_time_stamp_suffix()))
			print(" ")
			print_progress("Generating consistency report of the provided dataset in %s..."%(mic_consistency_check_info_path))
			mic_consistency_check_info_file = open(mic_consistency_check_info_path, "w")
			mic_consistency_check_info_file.write("# The consistency information about micrograph IDs that might have problmes with consistency among the provided dataset.\n")
			mic_consistency_check_info_file.write("# \n")
			# Loop over substring id list
			for mic_id_substr in global_entry_dict:
				mic_id_entry = global_entry_dict[mic_id_substr]
				
				consistency_messages = []
				# Check if associated micrograph path exists in source directory
				if not subkey_src_mic_path in mic_id_entry:
					mic_basename = mic_basename_pattern.replace("*", mic_id_substr)
					consistency_messages.append("    associated micrograph (%s) does not exist in the source directory (%s)."%(mic_basename, src_dir))
			
				# Check if associated micrograph basename exists in selection list
				if not subkey_select_mic_basename in mic_id_entry:
					mic_basename = mic_basename_pattern.replace("*", mic_id_substr)
					consistency_messages.append("    associated micrograph (%s) is not in the selection list (%s)."%(mic_basename, select_list_path))
			
				# Check if associated micrograph path does not exist in destination directory
				if subkey_dst_mic_path in mic_id_entry:
					mic_basename = mic_basename_pattern.replace("*", mic_id_substr)
					consistency_messages.append("    associated micrograph (%s) already exist in the destination directory (%s)."%(mic_basename, dst_dir))
			
				if len(consistency_messages) > 0:
					mic_consistency_check_info_file.write("Micrograph ID %s have inconsistency among provided dataset:\n"%(mic_id_substr))
					for consistency_message in consistency_messages:
						mic_consistency_check_info_file.write(consistency_message)
						mic_consistency_check_info_file.write("\n")
		
			# Close the consistency check file, if necessary
			mic_consistency_check_info_file.flush()
			mic_consistency_check_info_file.close()
		
		# Since mic_id_substr is once stored as the key of global_entry_dict and extracted with the key order
		# we need sort the valid_mic_id_substr_list here
		# print("MRK_DEBUG: before sort, valid_mic_id_substr_list := ", valid_mic_id_substr_list)
		valid_mic_id_substr_list.sort()
		# print("MRK_DEBUG: after sort, valid_mic_id_substr_list := ", valid_mic_id_substr_list)
		
		# --------------------------------------------------------------------------------
		# Print out the summary of input consistency
		# --------------------------------------------------------------------------------
		print(" ")
		print_progress("Summary of dataset consistency check...")
		print_progress("Detected                           : %6d"%(len(global_entry_dict)))
		print_progress("Valid                              : %6d"%(len(valid_mic_id_substr_list)))
		print_progress("Rejected by not found in both dirs : %6d"%(len(no_mic_in_both_dirs_id_substr_list)))
		print_progress("Rejected by already in dst dir     : %6d"%(len(already_in_dst_dir_mic_id_substr_list)))
		print_progress("Rejected by duplicated in dst dir  : %6d"%(len(duplicated_in_dst_dir_mic_id_substr_list)))
		print(" ")
		
		# --------------------------------------------------------------------------------
		# Save the list of duplicated_micrographs in duplicated_micrographs_DATE_TIME.txt 
		# under destination directory if necessary
		# --------------------------------------------------------------------------------
		if len(duplicated_in_dst_dir_mic_id_substr_list) > 0:
			duplicated_mic_list_path = os.path.join(record_dir, "duplicated_micrographs_%s.txt"%(get_time_stamp_suffix()))
			print_progress("Storing the list of duplicated micrographs in %s."%(duplicated_mic_list_path))
			print(" ")
			
			# Open the duplicated micrograph list file
			duplicated_mic_list_file = open(duplicated_mic_list_path, "w")
			for mic_id_substr in duplicated_in_dst_dir_mic_id_substr_list:
				duplicated_mic_basename = mic_basename_pattern.replace("*", mic_id_substr)
				duplicated_mic_list_file.write(duplicated_mic_basename)
				duplicated_mic_list_file.write("\n")
			# Close the duplicated micrograph list file
			duplicated_mic_list_file.flush()
			duplicated_mic_list_file.close()
		
		# --------------------------------------------------------------------------------
		# Clean up variables related to tracking of invalid (rejected) micrographs 
		# --------------------------------------------------------------------------------
		del no_mic_in_both_dirs_id_substr_list
		del already_in_dst_dir_mic_id_substr_list
		del duplicated_in_dst_dir_mic_id_substr_list

	# --------------------------------------------------------------------------------
	# Create destination directory
	# --------------------------------------------------------------------------------
	if not os.path.exists(dst_dir):
		print(" ")
		print_progress("Creating the destination directory (%s)..."%(dst_dir))
		os.mkdir(dst_dir)
	assert (os.path.exists(dst_dir))

	# --------------------------------------------------------------------------------
	# Move micrographs in selecting list form source directory to destination directory
	# --------------------------------------------------------------------------------
	# Prepare the counters for the global summary of micrographs
	n_moved_mics = 0
	
	if len(valid_mic_id_substr_list) > 0:
		print(" ")
		print_progress("Moving micrographs in the selecting list (%s) from the source directory (%s) to the destination directory (%s)..."%(select_list_path, src_dir, dst_dir))
		### print("Micrographs processed (including percent of progress):")
		### progress_percent_step = len(valid_mic_id_substr_list)*0.1 # Report every 10% of the number of micrograms

	# Loop over substring id list
	for mic_id_substr_idx, mic_id_substr in enumerate(valid_mic_id_substr_list):
		mic_id_entry = global_entry_dict[mic_id_substr]
		mic_basename = mic_id_entry[subkey_select_mic_basename]
		assert (mic_basename == mic_basename_pattern.replace("*", mic_id_substr))
		
		### # Print out progress if necessary
		### print("%s ---> % 2.2f%%" % (mic_basename, mic_id_substr_idx / progress_percent_step))
		
		# At this point, this micrograph
		# - must exist in source directory
		# - must NOT exist in destination directory
		# because of the consistency check above
		assert (subkey_src_mic_path in mic_id_entry)
		assert (os.path.exists(mic_id_entry[subkey_src_mic_path]))
		assert (not os.path.exists(os.path.join(dst_dir, mic_basename)))
		
		# Move this micrograph from input directory to output directory
		src_mic_path = mic_id_entry[subkey_src_mic_path]
		shutil.move(src_mic_path, dst_dir)
		n_moved_mics += 1
	
	# Print summary of processing
	print(" ")
	print_progress("Summary of processing...")
	print_progress("Moved      : %6d"%(n_moved_mics))
	print(" ")
	
### # ----------------------------------------------------------------------------------------
###	# NOTE: Toshio Moriya 2018/03/05
### # "reboxing" subcommand became obsolete because of "restacking" subcommand
### #
### # ----------------------------------------------------------------------------------------
### # Author #1: Christos Gatsogiannis 12/23/2015 (Christos.Gatsogiannis@mpi-dortmund.mpg.de)
### # Author #2: Toshio Moriya 02/20/2018 (toshio.moriya@mpi-dortmund.mpg.de)
### # 
### # Generate coordinates files and micrograph selection text file
### # based on information in the image headers of SPHIRE particle stack file
### # 
### # This command executes the following processes:
### # (1) Extract the following information stored in the headers of each particle image
### #     - source micrograph path 
### #     - box center coordinates within the micrograph
### #     - projection parameters 
### # (2) Convert the center coordinates to EMAN1 box coordinates format
### #     and save the results to output files.
### # (3) Transform the coordinates based on the projection parameters and user-provided 3D shift
### #     and save the results to output files.
### # (4) Save the list of extracted micrograph names to an output file.
### # 
### # ----------------------------------------------------------------------------------------
### # TEST COMMAND
### # cd /home/moriya/mrk_qa/mrktest_pipeline
### # rm -rf debug_mrkout_sxpipe_reboxing_fullset; sxpipe.py reboxing 'bdb:mrkout_pipe02_sxwindow#data_meridien' 'debug_mrkout_sxpipe_reboxing_fullset' --box_size=352
### # rm -rf debug_mrkout_sxpipe_reboxing_isac_substack; sxpipe.py reboxing 'bdb:mrkout_pipe03_sxisac/mrkout_sxpipe_isac_substack#isac_substack' 'debug_mrkout_sxpipe_reboxing_isac_substack' --box_size=352
### # rm -rf debug_mrkout_sxpipe_reboxing_sort3d_substack; sxpipe.py reboxing 'bdb:mrkout_pipe09_sxsort3d_depth_isac_subset_c5/mrkout_pipe09o02_e2bdb_makevstack#sort3d_depth_substack' 'debug_mrkout_sxpipe_reboxing_sort3d_substack' --box_size=352
### # ----------------------------------------------------------------------------------------
### def reboxing(args):
### 	# from sys import  *
### 	# import csv
### 	# import glob
### 	# import traceback
### 	# import math
### 	from EMAN2db   import db_check_dict
### 	from utilities import get_im, get_params_proj
### 	
### 	# ========================================================================================
### 	class coordinates_entry(object):
### 		def __init__(self, mic_basename):
### 			# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
### 			# class variables
### 			self.mic_basename = mic_basename
### 			self.eman1_original = []
### 			self.eman1_centered = []
### 			# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
### 
### 	# To make the execution exit upon fatal error by ERROR in global_def.py
### 	global_def.BATCH = True 
### 	
### 	# Define the name of this subcommand
### 	subcommand_name = "reboxing"
### 	
### 	# Check error conditions of arguments
### 	if not db_check_dict(args.input_stack_path, readonly=True) and not os.path.exists(args.input_stack_path):
### 		ERROR("Input stack file does not exist. Please check the input stack path and restart the program.", subcommand_name) # action=1 - fatal error, exit
### 	if os.path.exists(args.output_directory):
### 		ERROR("Output directory exists. Please change the name and restart the program.", subcommand_name) # action=1 - fatal error, exit
### 	
### 	# Check error conditions of options
### 	if args.box_size <= 0:
### 		ERROR("Invalid box size {}. Box size must be larger than zero. Please provide valid box size and restart the program.".format(args.box_size), subcommand_name) # action=1 - fatal error, exit
### 	if abs(args.shift3d_x) >= args.box_size//2:
### 		ERROR("Invalid 3D x-shift {}. 3D x-shift must be smaller than half of box size {}. Please provide valid box size and restart the program.".format(args.shift3d_x, args.box_size//2), subcommand_name) # action=1 - fatal error, exit
### 	if abs(args.shift3d_y) >= args.box_size//2:
### 		ERROR("Invalid 3D y-shift {}. 3D y-shift must be smaller than half of box size {}. Please provide valid box size and restart the program.".format(args.shift3d_y, args.box_size//2), subcommand_name) # action=1 - fatal error, exit
### 	if abs(args.shift3d_z) >= args.box_size//2:
### 		ERROR("Invalid 3D z-shift {}. 3D z-shift must be smaller than half of box size {}. Please provide valid box size and restart the program.".format(args.shift3d_z, args.box_size//2), subcommand_name) # action=1 - fatal error, exit
### 	
### 	# args.scale=1.0
### 	
### 	n_img = EMUtil.get_image_count(args.input_stack_path)
### 	# print(" ")
### 	print_progress("Found {} particle images in the input stack {}".format(n_img, args.input_stack_path))
### 
### 	# Define variables and constants used in the loop
### 	global_coordinates_dict = {} # mic basename is the key. contains original and centered box coordinates in EMAN1 format
### 	eman1_dummy = -1             # For 5th column of EMAN1 boxer format
### 	img = EMData()
### 	for img_id in xrange(n_img):
### 		# Load images 
### 		# img = get_im(args.input_stack_path, img_id)
### 		# Load only image header 
### 		img.read_image(args.input_stack_path, img_id, True)
### 		# Extract associated source micrograph path name from the image header
### 		mic_path = str(img.get_attr("ptcl_source_image"))
### 		mic_basename = os.path.basename(mic_path)
### 		if mic_basename not in global_coordinates_dict:
### 			global_coordinates_dict[mic_basename] = coordinates_entry(mic_basename)
### 			assert (global_coordinates_dict[mic_basename].mic_basename == mic_basename)
### 			# print_progress("Found new micrograph {}. Detected {} micrographs so far...".format(mic_basename, len(global_coordinates_dict)))
### 		
### 		# Extract the associated coordinates from the image header
### 		center_coordinate_x, center_coordinate_y = img.get_attr("ptcl_source_coord")
### 		# Compute the left bottom coordinates of box (EMAN1 box file format)
### 		assert(args.box_size > 0)
### 		# eman1_original_coordinate_x = center_coordinate_x - (args.box_size//2+1)
### 		# eman1_original_coordinate_y = center_coordinate_y - (args.box_size//2+1)
### 		# NOTE: 2018/02/21 Toshio Moriya
### 		# Currently, the following the way e2boxer.py calculates EMAN1 box format from particle center coordinates.
### 		eman1_original_coordinate_x = center_coordinate_x - (args.box_size//2)
### 		eman1_original_coordinate_y = center_coordinate_y - (args.box_size//2)
### 		global_coordinates_dict[mic_basename].eman1_original.append("{:6d} {:6d} {:6d} {:6d} {:6d}\n".format(eman1_original_coordinate_x, eman1_original_coordinate_y, args.box_size, args.box_size, eman1_dummy))
### 		
### 		# Extract the projection parameters from the image header
### 		proj_phi, proj_theta, proj_psi, proj_tx, proj_ty = get_params_proj(img)
### 		# Transform the coordinates according to projection parameters and user-provided 3D shift (corresponding to shifting the 3D volume)
### 		trans3x3 = Transform({"phi":float(proj_phi), "theta":float(proj_theta), "psi":float(proj_psi), "tx":float(proj_tx), "ty":float(proj_ty), "tz":0.0, "type":"spider"})
### 		origin_vec3d = Vec3f(float(args.shift3d_x), float(args.shift3d_y), float(args.shift3d_z))
### 		transformed_vec3d = trans3x3 * origin_vec3d
### 		shift2d_x = -1 * transformed_vec3d[0]
### 		shift2d_y = -1 * transformed_vec3d[1]
### 		# Transform and center the coordinates according to projection parameters and user-provided 3D shift (corresponding to shifting the 3D volume)
### 		eman1_centered_coordinate_x = int(round(eman1_original_coordinate_x + shift2d_x))
### 		eman1_centered_coordinate_y = int(round(eman1_original_coordinate_y + shift2d_y))
### 		global_coordinates_dict[mic_basename].eman1_centered.append("{:6d} {:6d} {:6d} {:6d} {:6d}\n".format(eman1_centered_coordinate_x, eman1_centered_coordinate_y, args.box_size, args.box_size, eman1_dummy))
### 	
### 	print(" ")
### 	print_progress("Found total of {} assocaited micrographs in the input stack {}".format(len(global_coordinates_dict), args.input_stack_path))
### 	
### 	assert (not os.path.exists(args.output_directory))
### 	os.mkdir(args.output_directory)
### 	
### 	mic_list_file_name = "micrographs.txt"
### 	mic_list_file_path = os.path.join(args.output_directory, mic_list_file_name)
### 	mic_list_file = open(mic_list_file_path, "w")
### 
### 	eman1_original_coordinates_subdir = "original"
### 	eman1_original_coordinates_suffix = '_original.box'
### 	os.mkdir(os.path.join(args.output_directory, eman1_original_coordinates_subdir))
### 	
### 	eman1_centered_coordinates_subdir = "centered"
### 	eman1_centered_coordinates_suffix = '_centered.box'
### 	os.mkdir(os.path.join(args.output_directory, eman1_centered_coordinates_subdir))
### 
### 	global_coordinates_list = sorted(global_coordinates_dict) # sort entries according to keys (i.e. micrograph basenames)
### 	assert (len(global_coordinates_list) == len(global_coordinates_dict))
### 	
### 	print(" ")
### 	for mic_basename in global_coordinates_list:
### 		# Write the mic base name to output file; micrograph selection text file
### 		mic_list_file.write("{}\n".format(mic_basename))
### 		
### 		mic_rootname, mic_extension = os.path.splitext(mic_basename)
### 		mic_entry = global_coordinates_dict[mic_basename]
### 		assert (len(mic_entry.eman1_original) == len(mic_entry.eman1_centered))
### 		
### 		# Save the original coordinates to output file; original EMAN1 box coordinate file or this micrograph
### 		eman1_original_coordinates_path = os.path.join(args.output_directory, eman1_original_coordinates_subdir, "{}{}".format(mic_rootname, eman1_original_coordinates_suffix))
### 		eman1_original_coordinates_file = open(eman1_original_coordinates_path, "w")
### 		for original_particle_coordinates in mic_entry.eman1_original:
### 			eman1_original_coordinates_file.write(original_particle_coordinates)
### 		eman1_original_coordinates_file.close()
### 		
### 		# Save the centered coordinates to output file; centered EMAN1 box coordinate file or this micrograph
### 		eman1_centered_coordinates_path = os.path.join(args.output_directory, eman1_centered_coordinates_subdir, "{}{}".format(mic_rootname, eman1_centered_coordinates_suffix))
### 		eman1_centered_coordinates_file = open(eman1_centered_coordinates_path, "w")
### 		for centered_particle_coordinates in mic_entry.eman1_centered:
### 			eman1_centered_coordinates_file.write(centered_particle_coordinates)
### 		eman1_centered_coordinates_file.close()
### 		
### 		# print(" ")
### 		# print_progress("Micrograph summary...")
### 		# print_progress("  Micrograph Name                : {}".format(mic_basename))
### 		# print_progress("  Extracted original coordinates : {:6d}".format(len(mic_entry.eman1_original)))
### 		# print_progress("  Extracted centered coordinates : {:6d}".format(len(mic_entry.eman1_centered)))
### 		# print_progress("  Saved original coordinates to  : {}".format(eman1_original_coordinates_path))
### 		# print_progress("  Saved centered coordinates to  : {}".format(eman1_centered_coordinates_path))
### 		print_progress(" {:6d} particle coordinates for {}...".format(len(mic_entry.eman1_original), mic_basename))
### 	
### 	mic_list_file.close()
### 	
### 	print(" ")
### 	print_progress("Global summary of processing...")
### 	print_progress("Processed particles                : {:6d}".format(n_img))
### 	print_progress("Extracted micrographs              : {:6d}".format(len(global_coordinates_list)))
### 	print_progress("Saved extracted micrograph list to : {}".format(mic_list_file_path))

# ----------------------------------------------------------------------------------------
# Author 1: Christos Gatsogiannis 12/23/2015 (Christos.Gatsogiannis@mpi-dortmund.mpg.de)
# Author 2: Toshio Moriya 03/02/2018 (toshio.moriya@mpi-dortmund.mpg.de)
# 
# --- Restacking ---
# Generate all necessary information to restack the input stack 
# (i.e. particle image ID list, CTF parameters list, projection parameters list) 
# while appling micrograph selection list. 
# Optionally, the command can directly output the virtual stack.  
# In addition, this command can be used to generate all parameters files for reboxing 
# (i.e. original/centered particle coordinates list files, CTF parameters list, 
# original/centered projection parameters list as well as micrograph selection file). 
# Optionally, user can provided a 3D shift to recenter the projection parameters and so the particle coordinates.
# 
# This command executes the following processes:
#  1. Extract the following information stored in the header of each particle image.
#     - source micrograph path (ptcl_source_image).
#     - CTF parameters if exist (ctf).
#     - projection parameters if exist (xform.projection).
#     - box center coordinates within the micrograph (ptcl_source_coord).
#  2. Save the list of extracted micrograph names to an output file.
#  3. If provided, apply the selection list to extracted micrographs.
#  4. Save the list of selected micrograph names to an output file.
#  5. Extract only particle image associating to selected micrographs.
#  6. Save the list of selected particle image IDs to an output file.
#  7. Save the CTF parameters of selected particle images to output file.
#  8. Save the original projection parameters of selected particle images to output file.
#  9. Transform the projection parameters of selected particle images  based on user-provided 3D shift, and then save the results to output files.
# 10. Convert the center coordinates to EMAN1 box coordinates format, and then save the results to output files.
# 11. Transform the coordinates based on the projection parameters and user-provided 3D shift, and then save the results to output files.
# 12. Create the output virtual stack if necessary
# 
# ----------------------------------------------------------------------------------------
# TEST COMMAND
# cd /home/moriya/mrk_qa/mrktest_pipeline
#
# sxpipe.py restacking --help
# 
# rm -rf debug_mrkout_sxpipe_restacking_fullset; sxpipe.py restacking 'bdb:mrkout_pipe02_sxwindow#data_meridien' 'debug_mrkout_sxpipe_restacking_fullset' --reboxing --rb_box_size=352
# rm -rf debug_mrkout_sxpipe_restacking_isac_substack; sxpipe.py restacking 'bdb:mrkout_pipe03_sxisac/mrkout_sxpipe_isac_substack#isac_substack' 'debug_mrkout_sxpipe_restacking_isac_substack' --reboxing --rb_box_size=352
# rm -rf debug_mrkout_sxpipe_restacking_sort3d_substack; sxpipe.py restacking 'bdb:mrkout_pipe09_sxsort3d_depth_isac_subset_c5/mrkout_pipe09o02_e2bdb_makevstack#sort3d_depth_substack' 'debug_mrkout_sxpipe_restacking_sort3d_substack' --reboxing --rb_box_size=352
# 
# rm -rf debug_mrkout_sxpipe_restacking_fullset_vstack; sxpipe.py restacking 'bdb:mrkout_pipe02_sxwindow#data_meridien' 'debug_mrkout_sxpipe_restacking_fullset_vstack' --reboxing --rb_box_size=352 --save_vstack 
# e2iminfo.py bdb:debug_mrkout_sxpipe_restacking_fullset_vstack#vstack
# rm -rf debug_mrkout_sxpipe_restacking_fullset_vstack; sxpipe.py restacking 'bdb:mrkout_pipe02_sxwindow#data_meridien' 'debug_mrkout_sxpipe_restacking_fullset_vstack' --reboxing --rb_box_size=352 --save_vstack  --sv_vstack_basename=mrkdebug_vstack
# e2iminfo.py bdb:debug_mrkout_sxpipe_restacking_fullset_vstack#mrkdebug_vstack
# 
# rm -rf debug_mrkout_sxpipe_restacking_fullset_selection; sxpipe.py restacking 'bdb:mrkout_pipe02_sxwindow#data_meridien' 'debug_mrkout_sxpipe_restacking_fullset_selection' --selection_list='debug_micrographs.txt' --reboxing --rb_box_size=352 --save_vstack  --sv_vstack_basename=mrkdebug_vstack
# e2iminfo.py bdb:debug_mrkout_sxpipe_restacking_fullset_selection#mrkdebug_vstack
# 
# rm -rf debug_mrkout_sxpipe_restacking_fullset_shift3d; sxpipe.py restacking 'bdb:mrkout_pipe02_sxwindow#data_meridien' 'debug_mrkout_sxpipe_restacking_fullset_shift3d' --selection_list='debug_micrographs.txt' --reboxing --rb_box_size=352 --save_vstack  --sv_vstack_basename=mrkdebug_vstack --shift3d_z=20
# e2iminfo.py bdb:debug_mrkout_sxpipe_restacking_fullset_shift3d#mrkdebug_vstack
# 
# rm -rf debug_mrkout_sxpipe_restacking_fullset_invalid_entry; sxpipe.py restacking 'bdb:mrkout_pipe02_sxwindow#data_meridien' 'debug_mrkout_sxpipe_restacking_fullset_invalid_entry' --selection_list='debug_micrographs_invalid_entry.txt' --reboxing --rb_box_size=352 --save_vstack  --sv_vstack_basename=mrkdebug_vstack
# e2iminfo.py bdb:debug_mrkout_sxpipe_restacking_fullset_vstack#mrkdebug_vstack
#
# ----------------------------------------------------------------------------------------
def restacking(args):
	# from sys import  *
	# import csv
	# import glob
	# import traceback
	# import math
	from EMAN2db   import db_check_dict
	from utilities import get_im, get_params_proj
	
	# ========================================================================================
	class SX_mic_entry(object):
		def __init__(self, mic_basename):
			# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
			# class variables
			self.mic_basename = mic_basename
			self.is_in_stack = False
			self.is_in_list = False
			self.img_id_list = []
			self.ctf_params_list = []
			self.original_proj_params_list = []
			self.original_coords_list = []
			self.centered_proj_params_list = []
			self.centered_coords_list = []
			# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><

	# Check MPI execution
	if SXmpi_run.n_mpi_procs > 1:
		assert (SXmpi_run.RUNNING_UNDER_MPI)
		error_status = ("The {} subcommand supports only a single process.".format(subcommand_name), getframeinfo(currentframe()))
		if_error_then_all_processes_exit_program(error_status)
	
	# To make the execution exit upon fatal error by ERROR in global_def.py
	global_def.BATCH = True 
	
	# Define the name of this subcommand
	# subcommand_name = "restacking"
	command_script_basename = os.path.basename(sys.argv[0])
	subcommand_name = "{} {}".format(command_script_basename, args.subcommand)
	
	# Check error conditions of arguments
	if not db_check_dict(args.input_bdb_stack_path, readonly=True):
		ERROR("Input BDB image stack file does not exist. Please check the input stack path and restart the program.", subcommand_name) # action=1 - fatal error, exit
	if os.path.exists(args.output_directory):
		ERROR("Output directory exists. Please change the name and restart the program.", subcommand_name) # action=1 - fatal error, exit
	
	# Check error conditions of options
	if args.selection_list is not None:
		if not os.path.exists(args.selection_list):
			ERROR("Micrograph selecting list does not exist. Please check the micrograph selecting list path and restart the program.", subcommand_name) # action=1 - fatal error, exit
		if os.path.splitext(args.selection_list)[1] != ".txt":
			ERROR("The extention of micrograph selecting list file must be .txt. Please check the micrograph selecting list path and restart the program.", subcommand_name) # action=1 - fatal error, exit
		args.sv_vstack_basename = args.sv_vstack_basename.strip()
		if args.sv_vstack_basename == "":
			ERROR("Output virtual stack basename cannot be empty string or only white spaces.", subcommand_name) # action=1 - fatal error, exit
	if not args.reboxing:
		img = EMData()
		img.read_image(args.input_bdb_stack_path, 0, True)
		img_size = img.get_xsize()
		assert (img.get_xsize() == img.get_ysize())
		if abs(args.shift3d_x) >= img_size//2:
			ERROR("Invalid 3D x-shift {}. 3D x-shift must be smaller than half of image size {}. Please provide valid values for x-shift and/or box size. Then, restart the program.".format(args.shift3d_x, args.img_size//2), subcommand_name) # action=1 - fatal error, exit
		if abs(args.shift3d_y) >= img_size//2:
			ERROR("Invalid 3D y-shift {}. 3D y-shift must be smaller than half of image size {}. Please provide valid values for y-shift and/or box size. Then, restart the program.".format(args.shift3d_y, args.img_size//2), subcommand_name) # action=1 - fatal error, exit
		if abs(args.shift3d_z) >= img_size//2:
			ERROR("Invalid 3D z-shift {}. 3D z-shift must be smaller than half of image size {}. Please provide valid values for z-shift and/or box size. Then, restart the program.".format(args.shift3d_z, args.img_size//2), subcommand_name) # action=1 - fatal error, exit
	if args.reboxing:
		if args.rb_box_size <= 0:
			ERROR("Invalid box size {}. Box size must be larger than zero. Please provide valid box size and restart the program.".format(args.rb_box_size), subcommand_name) # action=1 - fatal error, exit
		if abs(args.shift3d_x) >= args.rb_box_size//2:
			ERROR("Invalid 3D x-shift {}. 3D x-shift must be smaller than half of box size {}. Please provide valid values for x-shift and/or box size. Then, restart the program.".format(args.shift3d_x, args.rb_box_size//2), subcommand_name) # action=1 - fatal error, exit
		if abs(args.shift3d_y) >= args.rb_box_size//2:
			ERROR("Invalid 3D y-shift {}. 3D y-shift must be smaller than half of box size {}. Please provide valid values for y-shift and/or box size. Then, restart the program.".format(args.shift3d_y, args.rb_box_size//2), subcommand_name) # action=1 - fatal error, exit
		if abs(args.shift3d_z) >= args.rb_box_size//2:
			ERROR("Invalid 3D z-shift {}. 3D z-shift must be smaller than half of box size {}. Please provide valid values for z-shift and/or box size. Then, restart the program.".format(args.shift3d_z, args.rb_box_size//2), subcommand_name) # action=1 - fatal error, exit
	
	# args.scale=1.0
	
	# --------------------------------------------------------------------------------
	# Register micrograph basename found in the selection list
	# to the global entry dictionary
	# --------------------------------------------------------------------------------
	global_mic_dict = {}   # mic basename is the key. Organize the particle image information according to associated micrographs
	selected_mic_path_list = []
	# Generate micrograph lists according to the execution mode
	if args.selection_list is not None:
		assert (os.path.exists(args.selection_list))
		assert (os.path.splitext(args.selection_list)[1] == ".txt")
		
		# Generate the list of selected micrograph paths in the selection file
		print(" ")
		print_progress("----- Running with Selected Micrographs Mode -----")
		print(" ")
		print_progress("Checking the selection list {}...".format(args.selection_list))
		selected_mic_path_list = read_text_file(args.selection_list)
		
		# Check error condition of micrograph entry lists
		print("Found %d microgarph entries in %s." % (len(selected_mic_path_list), args.selection_list))
		if len(selected_mic_path_list) == 0:
			ERROR("No micrograph entries are found in the selection list file. Please check the micrograph selecting list and restart the program.", subcommand_name) # action=1 - fatal error, exit
		assert (len(selected_mic_path_list) > 0)
		if not isinstance(selected_mic_path_list[0], basestring):
			ERROR("Invalid format of the selection list file. The first column must contain micrograph paths in string type. Please check the micrograph selecting list and restart the program.", subcommand_name) # action=1 - fatal error, exit
		
		assert (len(selected_mic_path_list) > 0)
		selected_mic_directory = os.path.dirname(selected_mic_path_list[0])
		if selected_mic_directory != "":
			print_progress("    NOTE: Program disregards the directory paths in the selection list ({}).".format(selected_mic_directory))
		
		# Register micrograph basename found in the selection list
		# to the global entry dictionary
		print(" ")
		print_progress("Registering all micrographs in the selection list {}...".format(args.selection_list))
		for mic_path in selected_mic_path_list:
			mic_basename = os.path.basename(mic_path)
			if mic_basename in global_mic_dict:
				print_progress("WARNING!!! Micrograph {} is duplicated in the selection list {}. Ignoring this duplicated entry...".format(mic_path, args.selection_list))
				continue
			
			assert (mic_basename not in global_mic_dict)
			global_mic_dict[mic_basename] = SX_mic_entry(mic_basename)
			global_mic_dict[mic_basename].is_in_list = True
			assert (global_mic_dict[mic_basename].mic_basename == mic_basename)
			assert (global_mic_dict[mic_basename].is_in_stack == False)
			assert (global_mic_dict[mic_basename].is_in_list == True)
			# print_progress("Found new micrograph {}. Detected {} micrographs so far...".format(mic_basename, len(global_mic_dict)))
	else:
		assert (args.selection_list is None)
		print(" ")
		print_progress("----- Running with All Micrographs Mode -----")

	# --------------------------------------------------------------------------------
	# Register micrograph basename found in the selection list
	# to the global entry dictionary
	# --------------------------------------------------------------------------------
	print(" ")
	print_progress("Checking the input stack {}...".format(args.input_bdb_stack_path))
	n_img = EMUtil.get_image_count(args.input_bdb_stack_path)
	# print(" ")
	print_progress("Found {} particle images in the input stack {}".format(n_img, args.input_bdb_stack_path))
	
	# --------------------------------------------------------------------------------
	# Register micrograph basenames found in the stack to the global micrographs dictionary
	# --------------------------------------------------------------------------------
	# Define variables and constants used in the loop
	eman1_dummy = -1       # For 5th column of EMAN1 boxer format
	missing_ctf_params_counter = 0
	missing_proj_params_counter = 0
	img = EMData()
	for img_id in xrange(n_img):
		# Load images 
		# img = get_im(args.input_bdb_stack_path, img_id)
		img.read_image(args.input_bdb_stack_path, img_id, True)
		# Extract associated source micrograph path name from the image header
		mic_path = str(img.get_attr("ptcl_source_image"))
		mic_basename = os.path.basename(mic_path)
		# Case 1: This micrograph has been not registered yet.
		if mic_basename not in global_mic_dict:
			mic_entry = SX_mic_entry(mic_basename)
			assert (mic_entry.mic_basename == mic_basename)
			assert (not mic_entry.is_in_stack)
			mic_entry.is_in_stack = True
			assert (mic_entry.is_in_stack)
			if args.selection_list is None:
				# If selection list is not provided, process as if all micrographs in the input stack are selected
				assert (not mic_entry.is_in_list)
				mic_entry.is_in_list = True
				assert (mic_entry.is_in_list)
			global_mic_dict[mic_basename] = mic_entry
			assert (global_mic_dict[mic_basename].mic_basename == mic_basename)
			# print_progress("Found new micrograph {}. Detected {} micrographs so far...".format(mic_basename, len(global_mic_dict)))
		# Case 2: This micrograph has been registered already. (1) This is in the selection list or (2) Not first incidence in the input stack
		else:
			assert (mic_basename in global_mic_dict)
			mic_entry = global_mic_dict[mic_basename]
			if not mic_entry.is_in_stack:
				mic_entry.is_in_stack = True
			assert (mic_entry.is_in_stack)
			
		global_mic_dict[mic_basename].img_id_list.append(img_id)
		
		# Initialise CTF parameters same way as the dummy CTF of sxwindow.py
		defocus = 0.0
		cs      = 0.0
		voltage = 300.0
		apix    = 1.0
		bfactor = 0.0
		ampcont = 100.0
		dfdiff  = 0.0
		dfang   = 0.0
		# Extract the CTF parameters from the image header
		if img.has_attr("ctf"):
			defocus, cs, voltage, apix, bfactor, ampcont, dfdiff, dfang = get_ctf(img)
		else:
			missing_ctf_params_counter += 1
		
		# Register CTF parameters to global micrograph dictionary
		global_mic_dict[mic_basename].ctf_params_list.append("{:15.5f} {:15.5f} {:15.5f} {:15.5f} {:15.5f} {:15.5f} {:15.5f} {:15.5f}".format(defocus, cs, voltage, apix, bfactor, ampcont, dfdiff, dfang))
		
		# Extract the projection parameters from the image header
		proj_phi   = 0.0
		proj_theta = 0.0
		proj_psi   = 0.0
		proj_tx    = 0.0
		proj_ty    = 0.0
		if img.has_attr("xform.projection"):
			proj_phi, proj_theta, proj_psi, proj_tx, proj_ty = get_params_proj(img)
		else:
			missing_proj_params_counter += 1
			
		# Register original projection parameters to global micrograph dictionary
		global_mic_dict[mic_basename].original_proj_params_list.append("{:15.5f} {:15.5f} {:15.5f} {:15.5f} {:15.5f} ".format(proj_phi, proj_theta, proj_psi, proj_tx, proj_ty))
		
		# Transform the coordinates according to projection parameters and user-provided 3D shift (corresponding to shifting the 3D volume)
		trans3x3 = Transform({"phi":float(proj_phi), "theta":float(proj_theta), "psi":float(proj_psi), "tx":float(proj_tx), "ty":float(proj_ty), "tz":0.0, "type":"spider"})
		origin_vec3d = Vec3f(float(args.shift3d_x), float(args.shift3d_y), float(args.shift3d_z))
		transformed_vec3d = trans3x3 * origin_vec3d
		shift2d_x = -1 * transformed_vec3d[0]
		shift2d_y = -1 * transformed_vec3d[1]
		
		# Register centered projection parameters to global micrograph dictionary
		global_mic_dict[mic_basename].centered_proj_params_list.append("{:15.5f} {:15.5f} {:15.5f} {:15.5f} {:15.5f} ".format(proj_phi, proj_theta, proj_psi, 0.0, 0.0))
		
		if args.reboxing:
			# Extract the associated coordinates from the image header
			center_coordinate_x, center_coordinate_y = img.get_attr("ptcl_source_coord")
			# Compute the left bottom coordinates of box (EMAN1 box file format)
			assert(args.rb_box_size > 0)
			# original_coordinate_x = center_coordinate_x - (args.rb_box_size//2+1)
			# original_coordinate_y = center_coordinate_y - (args.rb_box_size//2+1)
			# NOTE: 2018/02/21 Toshio Moriya
			# Currently, the following the way e2boxer.py calculates EMAN1 box format from particle center coordinates.
			original_coordinate_x = center_coordinate_x - (args.rb_box_size//2)
			original_coordinate_y = center_coordinate_y - (args.rb_box_size//2)
			global_mic_dict[mic_basename].original_coords_list.append("{:6d} {:6d} {:6d} {:6d} {:6d}\n".format(original_coordinate_x, original_coordinate_y, args.rb_box_size, args.rb_box_size, eman1_dummy))
			
			# Transform and center the coordinates according to projection parameters and user-provided 3D shift (corresponding to shifting the 3D volume)
			centered_coordinate_x = int(round(original_coordinate_x + shift2d_x))
			centered_coordinate_y = int(round(original_coordinate_y + shift2d_y))
			global_mic_dict[mic_basename].centered_coords_list.append("{:6d} {:6d} {:6d} {:6d} {:6d}\n".format(centered_coordinate_x, centered_coordinate_y, args.rb_box_size, args.rb_box_size, eman1_dummy))
			
#	print(" ")
#	print_progress("Found total of {} assocaited micrographs in the input stack {}.".format(len(global_mic_dict), args.input_bdb_stack_path))
	
	if missing_ctf_params_counter > 0:
		print(" ")
		print_progress("WARNING!!! The CTF parameters (ctf header entry) are missing from {} out of {} particle images in the input stack {}.".format(missing_proj_params_counter, n_img, args.input_bdb_stack_path))
		print_progress("           The program automatically sets dummy CTF parameters for these particle images.")
		
	if missing_proj_params_counter > 0:
		print(" ")
		print_progress("WARNING!!! The projection parameters (xform.projection header entry) are missing from {} out of {} particle images in the input stack {}.".format(missing_proj_params_counter, n_img, args.input_bdb_stack_path))
		print_progress("           The program automtically set the prjection parameters to all zeros (null alignment).")

	mic_basename_list_of_input_stack = []
	mic_basename_list_of_output_stack = []
	
	print(" ")
	print_progress("Checking consistency of the provided dataset ...")
		
	# Loop over all registed micrograph basename
	for mic_basename in global_mic_dict:
		mic_entry = global_mic_dict[mic_basename]
		assert (mic_entry.mic_basename == mic_basename)
		
		if mic_entry.is_in_stack:
			mic_basename_list_of_input_stack.append(mic_basename)
			if mic_entry.is_in_list:
				# This is only condition (expected typical case) where we have to output the info of this micrograph
				mic_basename_list_of_output_stack.append(mic_basename)
			else:
				assert (not mic_entry.is_in_list)
				print_progress("    Micrograph {} is in the stack but not in the selection list.".format(mic_basename))
		else:
			assert (not mic_entry.is_in_stack)
			if mic_entry.is_in_list:
				print_progress("    Micrograph {} is in the selection list but not in the stack.".format(mic_basename))
			else:
				assert (not mic_entry.is_in_list)
				assert ("Unreachable!!! This condition should never happen!")
	
	assert (not os.path.exists(args.output_directory))
	os.mkdir(args.output_directory)

	print(" ")
	print_progress("Found total of {} micrographs in the input stack {}.".format(len(mic_basename_list_of_input_stack), args.input_bdb_stack_path))
	
	mic_basename_list_of_input_stack.sort()
	mic_basename_list_of_output_stack.sort()
	
	assert (len(mic_basename_list_of_input_stack) >= len(mic_basename_list_of_output_stack))
	assert (len(global_mic_dict) >= len(mic_basename_list_of_input_stack))  # Extra micrographs might be only selection file
	assert (len(global_mic_dict) >= len(mic_basename_list_of_output_stack))
	
	if args.selection_list is not None:
		print(" ")
		print_progress("Found total of {} valid micrographs for the output.".format(len(mic_basename_list_of_output_stack)))
	else: 
		assert (args.selection_list is None)
		assert (len(mic_basename_list_of_input_stack) == len(mic_basename_list_of_output_stack))
	
	print(" ")
	input_mic_list_file_name = "micrographs_in_input_stack.txt"
	input_mic_list_file_path = os.path.join(args.output_directory, input_mic_list_file_name)
	print_progress("Saving the list of micrographs found in the input stack {} to {}...".format(args.input_bdb_stack_path, input_mic_list_file_path))
	input_mic_list_file = open(input_mic_list_file_path, "w")
	for mic_basename in mic_basename_list_of_input_stack:
		# Write the mic base name to output file; micrograph selection text file
		input_mic_list_file.write("{}\n".format(mic_basename))
	input_mic_list_file.close()

	print(" ")
	output_mic_list_file_name = "micrographs_in_output_dataset.txt"
	output_mic_list_file_path = os.path.join(args.output_directory, output_mic_list_file_name)
	print_progress("Saving the list of valid micrograph names for the output to {}...".format(output_mic_list_file_path))
	output_mic_list_file = open(output_mic_list_file_path, "w")
	for mic_basename in mic_basename_list_of_output_stack:
		# Write the mic base name to output file; micrograph selection text file
		output_mic_list_file.write("{}\n".format(mic_basename))
	output_mic_list_file.close()

	ctf_params_list_file_name = "ctf_params_for_output_dataset.txt"
	ctf_params_list_file_path = os.path.join(args.output_directory, ctf_params_list_file_name)
	ctf_params_list_file = open(ctf_params_list_file_path, "w")

	original_proj_params_list_file_name = "original_proj_params_for_output_dataset.txt"
	original_proj_params_list_file_path = os.path.join(args.output_directory, original_proj_params_list_file_name)
	original_proj_params_list_file = open(original_proj_params_list_file_path, "w")

	centered_proj_params_list_file_name = "centered_proj_params_for_output_dataset.txt"
	centered_proj_params_list_file_path = os.path.join(args.output_directory, centered_proj_params_list_file_name)
	centered_proj_params_list_file = open(centered_proj_params_list_file_path, "w")

	if args.reboxing:
		original_coords_list_subdir = "original"
		original_coords_list_suffix = '_original.box'
		os.mkdir(os.path.join(args.output_directory, original_coords_list_subdir))
	
		centered_coords_list_subdir = "centered"
		centered_coords_list_suffix = '_centered.box'
		os.mkdir(os.path.join(args.output_directory, centered_coords_list_subdir))

	global_output_image_id_list = []
	global_ctf_params_counters = 0
	global_original_proj_params_counters = 0
	global_centered_proj_params_counters = 0
	global_original_coords_counters = 0
	global_centered_coords_counters = 0
	print(" ")
	for mic_basename in mic_basename_list_of_output_stack:
		mic_entry = global_mic_dict[mic_basename]
		assert (mic_basename == mic_entry.mic_basename)
		
		# Append particle ID list to global output stack particle ID list
		global_output_image_id_list += mic_entry.img_id_list
		
		# Count up total number of CTF parameters
		assert (len(mic_entry.ctf_params_list) == len(mic_entry.img_id_list))
		global_ctf_params_counters += len(mic_entry.ctf_params_list)
		assert (global_ctf_params_counters == len(global_output_image_id_list))
		
		# Count up total number of project parameters
		assert (len(mic_entry.original_proj_params_list) == len(mic_entry.img_id_list))
		assert (len(mic_entry.centered_proj_params_list) == len(mic_entry.img_id_list))
		global_original_proj_params_counters += len(mic_entry.original_proj_params_list)
		global_centered_proj_params_counters += len(mic_entry.centered_proj_params_list)
		assert (global_original_proj_params_counters == len(global_output_image_id_list))
		assert (global_centered_proj_params_counters == len(global_output_image_id_list))
		
		for ctf_params in mic_entry.ctf_params_list:
			ctf_params_list_file.write("{}\n".format(ctf_params))
		
		for original_proj_params in mic_entry.original_proj_params_list:
			original_proj_params_list_file.write("{}\n".format(original_proj_params))
		
		for centered_proj_params in mic_entry.centered_proj_params_list:
			centered_proj_params_list_file.write("{}\n".format(centered_proj_params))

		if args.reboxing:
			mic_rootname, mic_extension = os.path.splitext(mic_basename)
			assert (len(mic_entry.original_coords_list) == len(mic_entry.img_id_list))
			assert (len(mic_entry.centered_coords_list) == len(mic_entry.img_id_list))
			# Count up total number of coordinates
			global_original_coords_counters += len(mic_entry.original_coords_list)
			global_centered_coords_counters += len(mic_entry.centered_coords_list)
			assert (global_original_coords_counters == len(global_output_image_id_list))
			assert (global_centered_coords_counters == len(global_output_image_id_list))
		
			# Save the original coordinates to output file; original EMAN1 box coordinate file or this micrograph
			original_coords_list_path = os.path.join(args.output_directory, original_coords_list_subdir, "{}{}".format(mic_rootname, original_coords_list_suffix))
			original_coords_list_file = open(original_coords_list_path, "w")
			for original_coords in mic_entry.original_coords_list:
				original_coords_list_file.write(original_coords)
			original_coords_list_file.close()
		
			# Save the centered coordinates to output file; centered EMAN1 box coordinate file or this micrograph
			centered_coords_list_path = os.path.join(args.output_directory, centered_coords_list_subdir, "{}{}".format(mic_rootname, centered_coords_list_suffix))
			centered_coords_list_file = open(centered_coords_list_path, "w")
			for centered_particle_coordinates in mic_entry.centered_coords_list:
				centered_coords_list_file.write(centered_particle_coordinates)
			centered_coords_list_file.close()
		
		# print(" ")
		# print_progress("Micrograph summary...")
		# print_progress("  Micrograph Name                : {}".format(mic_basename))
		# print_progress("  Extracted particle image ID    : {:6d}".format(len(mic_entry.img_id_list)))
		# print_progress("  Original projection parameters : {:6d}".format(len(mic_entry.original_proj_params_list)))
		# print_progress("  Centered projection parameters : {:6d}".format(len(mic_entry.centered_proj_params_list)))
		print_progress(" {:6d} particles in {}...".format(len(mic_entry.img_id_list), mic_basename))
		if args.reboxing:
		#	print_progress("  Extracted original coordinates : {:6d}".format(len(mic_entry.original_coords_list)))
		#	print_progress("  Extracted centered coordinates : {:6d}".format(len(mic_entry.centered_coords_list)))
		#	print_progress("  Saved original coordinates to  : {}".format(original_coords_list_path))
		#	print_progress("  Saved centered coordinates to  : {}".format(centered_coords_list_path))
		#	print_progress(" {:6d} particle coordinates for {}...".format(len(mic_entry.original_coords_list), mic_basename))
			assert (len(mic_entry.original_coords_list) == len(mic_entry.img_id_list))
	
	ctf_params_list_file.close()
	original_proj_params_list_file.close()
	centered_proj_params_list_file.close()
		
	print(" ")
	output_particle_id_list_file_name = "input_stack_particle_id_for_output_dataset.txt"
	output_particle_id_list_file_path = os.path.join(args.output_directory, output_particle_id_list_file_name)
	print_progress("Saving the list of input stack particle IDs for the output dataset to {}...".format(output_particle_id_list_file_path))
	output_particle_id_list_file = open(output_particle_id_list_file_path, "w")
	for output_image_id_list in global_output_image_id_list:
		# Write the mic base name to output file; micrograph selection text file
		output_particle_id_list_file.write("{}\n".format(output_image_id_list))
	output_particle_id_list_file.flush()
	output_particle_id_list_file.close()

	if args.save_vstack:
		# Create virtual stack for output stack
		assert (args.output_directory != "")
		print(" ")
		print_progress("Creating output stack as a virtual stack...")
		virtual_bdb_stack_path = "bdb:{}#{}".format(args.output_directory, args.sv_vstack_basename)
		cmd_line = "e2bdb.py {} --makevstack={} --list={}".format(args.input_bdb_stack_path, virtual_bdb_stack_path, output_particle_id_list_file_path)
		status = cmdexecute(cmd_line)
		if status == 0: ERROR("\"{}\" execution failed. Exiting...".format(cmd_line), subcommand_name) # action=1 - fatal error, exit
		assert (EMUtil.get_image_count(virtual_bdb_stack_path) == len(global_output_image_id_list))

	print(" ")
	print_progress("Global summary of processing...")
	print_progress("Num. of extracted micrographs in selection list       : {:6d}".format(len(selected_mic_path_list)))
	print_progress("Num. of extracted micrographs in input stack          : {:6d}".format(len(mic_basename_list_of_input_stack)))
	print_progress("Saved input stack micrograph list to                  : {}".format(input_mic_list_file_path))
	print_progress("Num. of valid micrographs for output dataset          : {:6d}".format(len(mic_basename_list_of_output_stack)))
	print_progress("Saved output dataset micrograph list to               : {}".format(output_mic_list_file_path))
	print_progress(" ")
	print_progress("Num. of detected input stack particles                : {:6d}".format(n_img))
	print_progress("Num. of particle IDs in output dataset                : {:6d}".format(len(global_output_image_id_list)))
	print_progress("Saved particle ID list of output dataset to           : {}".format(output_particle_id_list_file_path))
	print_progress("Num. of CTF params in output dataset                  : {:6d}".format(global_ctf_params_counters))
	print_progress("Saved CTF params list of output dataset to            : {}".format(ctf_params_list_file_path))
	print_progress("Num. of original proj. params in output dataset       : {:6d}".format(global_original_proj_params_counters))
	print_progress("Saved original proj. params list of output dataset to : {}".format(original_proj_params_list_file_path))
	print_progress("Num. of ceneterd proj. params in output dataset       : {:6d}".format(global_centered_proj_params_counters))
	print_progress("Saved ceneterd proj. params list of output dataset to : {}".format(centered_proj_params_list_file_path))
	if args.reboxing:
		print_progress("Num. of original coordinates in output dataset        : {:6d}".format(global_original_coords_counters))
		print_progress("Saved original coordinates files in                   : {}".format(os.path.join(args.output_directory, original_coords_list_subdir)))
		print_progress("Num. of centered coordinates in output dataset        : {:6d}".format(global_centered_coords_counters))
		print_progress("Saved centered coordinates files in                   : {}".format(os.path.join(args.output_directory, centered_coords_list_subdir)))
	if args.save_vstack:
		print_progress("Save output stack as                                  : {}".format(virtual_bdb_stack_path))

# ========================================================================================
# Main function
# ========================================================================================
def main():
	# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	# Set up argument parser (supports subcommand)
	# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	parser = argparse.ArgumentParser(description="The collection of SPHIRE small pipleline tools.")
	parser.add_argument("--version", action="version", version=SPARXVERSION)
	# subparsers = parser.add_subparsers(title="subcommands", description="valid subcommands", help="additional help")
	subparsers = parser.add_subparsers(help="sub-command help", dest="subcommand")

	# create the subparser for the "isac_substack" subcommand
	parser_isac_subset = subparsers.add_parser("isac_substack", help="ISAC2 Stack Subset: Create virtual subset stack consisting from ISAC2 accounted particles by retrieving particle numbers associated with the ISAC2 or Beautifier class averages. The command also saves a list text file containing the retrieved original image numbers and 2D alingment parameters. In addition, it stores the 2D alingment parameters to stack header.")
	parser_isac_subset.add_argument("input_bdb_stack_path",          type=str,                            help="Input BDB image stack: Specify the same BDB image stack used for the associated ISAC2 run. (default required string)")
	parser_isac_subset.add_argument("input_run_dir",                 type=str,                            help="ISAC2 or Beautifier run output directory: Specify output directory of an ISAC2 or Beautifier run as an input to this command. From this directory, the program extracts the shrink ratio and 2D alingment parameters of the ISAC2 run or local 2D alingment parameters of the Beautifier run. (default required string)")
	parser_isac_subset.add_argument("output_directory",              type=str,                            help="Output directory: The results will be written here. This directory will be created automatically and it must not exist previously. (default required string)")
	parser_isac_subset.add_argument("--isac_class_avgs_path",        type=str,  default="",               help="ISAC2 or Beautifier class averages path: Specify path to a file containg ISAC2 or Beautifier class averages. The calss averages can be fullset or selected subset, as long as they are associated with the input BDB image stack and contain class member information stored in the headers. By default, the program uses the same deafult name of ordered class averages in ISAC2 or Beautifier (i.e. ordered_class_averages.hdf). (default none)")
	parser_isac_subset.add_argument("--substack_basename",           type=str,  default="isac_substack",  help="Stack subset basename: Specify the basename of ISAC2 stack subset file. It cannot be empty string or only white spaces. (default isac_substack)")
	### 
	### NOTE: Toshio Moriya 2018/01/13
	### The following options are not implemented yet.
	### parser_isac_subset.add_argument("--isac_class_id",               type=int,             default=-1,     help="ISAC class average ID: Retrieve only particle members of the specifed ISAC class. By default, retrieve from all classes. (default -1)")
	### parser_isac_subset.add_argument("--no_virtual_stack",            action="store_true",  default=False,  help="Do not create virtual stack: Use this option to create only the particle ID list text file associated with the ISAC class averages. (default False)")
	### parser_isac_subset.add_argument("--no_import_align2d",           action="store_true",  default=False,  help="Do not import alignment:  (default False)")
	parser_isac_subset.set_defaults(func=isac_substack)
	
	# create the subparser for the "resample_micrographs" subcommand
	parser_resample_micrographs = subparsers.add_parser("resample_micrographs", help="Resample Micrographs: Resample micrographs in input directory specified by input micrograph pattern with user-specified ratio. This operation changes the image dimensitions and the pixel size.")
	parser_resample_micrographs.add_argument("input_micrograph_pattern",  type=str,                             help="Input micrograph path pattern: Specify path pattern of input micrographs with a wild card (*). Use the wild card to indicate the place of variable part of the file names (e.g. serial number, time stamp, and etc). The path pattern must be enclosed by single quotes (\') or double quotes (\"). (Note: sxgui.py automatically adds single quotes (\')). The substring at the variable part must be same between the associated pair of input micrograph and coordinates file. bdb files can not be selected as input micrographs. (default required string)")
	parser_resample_micrographs.add_argument("output_directory",          type=str,                             help="Output directory: The results will be written here. This directory will be created automatically and it must not exist previously. (default required string)")
	parser_resample_micrographs.add_argument("--resample_ratio",          type=float,           default=None,   help="Resampling ratio: Specify ratio between new and original pixel size. Use a value between 0.0 and 1.0 (exclusive both ends). (default required float)")
	parser_resample_micrographs.add_argument("--selection_list",          type=str,             default=None,   help="Micrograph selecting list: Specify a name of micrograph selection list text file for Selected Micrographs Mode. The file extension must be \'.txt\'. Alternatively, the file name of a single micrograph can be specified for Single Micrograph Mode. (default none)")
	parser_resample_micrographs.add_argument("--check_consistency",       action="store_true",  default=False,  help="Check consistency of dataset: Create a text file containing the list of Micrograph ID entries might have inconsitency among the provided dataset. (i.e. mic_consistency_check_info_TIMESTAMP.txt). (default False)")
	parser_resample_micrographs.set_defaults(func=resample_micrographs)

	# create the subparser for the "organize_micrographs" subcommand
	parser_organize_micrographs = subparsers.add_parser("organize_micrographs", help="Organize Micrographs/Movies: Organize micrographs/movies by moving micrographs/movies in a selecting file from a source directory (specified by source micrographs/movies pattern) to a destination directory.")
	parser_organize_micrographs.add_argument("source_micrograph_pattern",    type=str,                                help="Source micrograph/movies path pattern: Specify path pattern of source micrographs/movies with a wild card (*). Use the wild card to indicate the place of variable part of the file names (e.g. serial number, time stamp, and etc). The path pattern must be enclosed by single quotes (\') or double quotes (\"). (Note: sxgui.py automatically adds single quotes (\')). The substring at the variable part must be same between each associated pair of micrograph/movie names. bdb files can not be selected as source micrographs/movies. (default required string)")
	parser_organize_micrographs.add_argument("selection_list",               type=str,                                help="Micrograph/Movie selection file: Specify a path of text file containing a list of selected micrograph/movie names or paths. The file extension must be \'.txt\'. The directory path of each entry will be ignored if there are any. (default required string)")
	parser_organize_micrographs.add_argument("destination_directory",        type=str,                                help="Destination directory: The micrographs/movies in selecting list will be moved to this directory. This directory will be created automatically if it does not exist. (default required string)")
	parser_organize_micrographs.add_argument("--reverse",                    action="store_true",  default=False,     help="Reverse operation: Move back micrographs/movies from the destination directory to the source directory. Please use this option to restore the previously-moved micrographs/movies. (default False)")
	parser_organize_micrographs.add_argument("--check_consistency",          action="store_true",  default=False,     help="Check consistency of dataset: Create a text file containing the list of micrograph/movie ID entries might have inconsitency among the provided dataset. (i.e. mic_consistency_check_info.txt). (default False)")
	parser_organize_micrographs.set_defaults(func=organize_micrographs)
	
###	# NOTE: Toshio Moriya 2018/03/05
### # "reboxing" subcommand became obsolete because of "restacking" subcommand
### #
###	# create the subparser for the "reboxing" subcommand
###	parser_reboxing = subparsers.add_parser("reboxing", help="Reboxing: Extract coordinates from the input stack, then center them according to projection parameters in the header and user-provided 3D shift")
###	parser_reboxing.add_argument("input_stack_path",    type=str,                  help="Input image stack: Specify path to input particle stack. (default required string)")
###	parser_reboxing.add_argument("output_directory",    type=str,                  help="Output directory: The results will be written here. This directory will be created automatically and it must not exist previously. (default required string)")
###	parser_reboxing.add_argument("--box_size",          type=int,    default=0,    help="Particle box size [Pixels]: The x and y dimensions of square area to be windowed. (default 0)")
###	parser_reboxing.add_argument("--shift3d_x",         type=int,    default=0,    help="3D x-shift [Pixels]: User-provided 3D x-shift corresponding to shifting the 3D volume along x-axis. (default 0)")
###	parser_reboxing.add_argument("--shift3d_y",         type=int,    default=0,    help="3D y-shift [Pixels]: User-provided 3D y-shift corresponding to shifting the 3D volume along y-axis. (default 0)")
###	parser_reboxing.add_argument("--shift3d_z",         type=int,    default=0,    help="3D z-shift [Pixels]: User-provided 3D z-shift corresponding to shifting the 3D volume along z-axis. (default 0)")
###	parser_reboxing.set_defaults(func=reboxing)

	# create the subparser for the "restacking" subcommand
	parser_restacking = subparsers.add_parser("restacking", help="Restacking: Generate all necessary information to restack the input stack (i.e. particle image ID list and projection parameters list) while appling micrograph selection list. Optinally, the command can directly output the virtual stack.  Also, this command can be used to generate all parameters files for reboxing (i.e. original/centered particle coordinates list files, original/centered projection parameters list as well as micrograph selection file). Optionally, user can provided a 3D shift to recenter the projection parameters and so the particle coordinates.")
	parser_restacking.add_argument("input_bdb_stack_path",    type=str,                                help="Input BDB image stack: Specify the input BDB image stack. (default required string)")
	parser_restacking.add_argument("output_directory",        type=str,                                help="Output directory: The results will be written here. This directory will be created automatically and it must not exist previously. (default required string)")
	parser_restacking.add_argument("--selection_list",        type=str,             default=None,      help="Micrograph/Movie selection file: Specify path to text file containing a list of selected micrograph/movie names or paths. The particles associated with the micrographs/movies in this list will be processed. The file extension must be \'.txt\'. The directory path of each entry will be ignored if there are any. (default none)")
	parser_restacking.add_argument("--shift3d_x",             type=int,             default=0,         help="3D x-shift [Pixels]: Provide 3D x-shift corresponding to shifting the 3D volume along x-axis. (default 0)")
	parser_restacking.add_argument("--shift3d_y",             type=int,             default=0,         help="3D y-shift [Pixels]: Provide 3D y-shift corresponding to shifting the 3D volume along y-axis. (default 0)")
	parser_restacking.add_argument("--shift3d_z",             type=int,             default=0,         help="3D z-shift [Pixels]: Provide 3D z-shift corresponding to shifting the 3D volume along z-axis. (default 0)")
	parser_restacking.add_argument("--save_vstack",           action="store_true",  default=False,     help="Save virtual stack: Use this option to save the virtual stack. By default, the virtual stack will not be generated. (default False)")
	parser_restacking.add_argument("--sv_vstack_basename",    type=str,             default="vstack",  help="Virtual stack basename: For --save_vstack, specify the basename of output virtual stack file. It cannot be empty string or only white spaces. (default vstack)")
	parser_restacking.add_argument("--reboxing",              action="store_true",  default=False,     help="Generate reboxing information: Prepare reboxing by extracting coordinates from the input stack headers, then center them according to projection parameters in the header and user-provided 3D shift. If the headers do not contain projection parameters, the program assumes the prjection parameters are all zeros (null alignment). (default False)")
	parser_restacking.add_argument("--rb_box_size",           type=int,             default=0,         help="Particle box size [Pixels]: For --reboxing option, specify the x and y dimensions of square area to be windowed. (default 0)")

	parser_restacking.set_defaults(func=restacking)

	# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	# Run specified subcommand
	# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	args = parser.parse_args() # Get namespace object from parser
	# args_dict = vars(parser.parse_args()) # convert it to dictionary object
	# print (args_dict)
	
	# ------------------------------------------------------------------------------------
	# Set up MPI related variables
	# ------------------------------------------------------------------------------------
	# Detect if program is running under MPI
	SXmpi_run.setup()
	
	# ------------------------------------------------------------------------------------
	# Execute command
	# ------------------------------------------------------------------------------------
	# Print command line
	if SXmpi_run.is_main_proc():
		print(" ")
		print_progress(get_cmd_line())
		print(" ")
	
	# Call the associated function of the specified subcommand
	args.func(args)

	if SXmpi_run.is_main_proc():
		print(" ")
		print_progress("DONE!!!")
		print(" ")

	# ------------------------------------------------------------------------------------
	# Clean up MPI related variables
	# ------------------------------------------------------------------------------------
	SXmpi_run.cleanup()
	
# ----------------------------------------------------------------------------------------
if __name__ == "__main__":
	main()

# ========================================================================================
# END OF SCRIPT
# ========================================================================================
