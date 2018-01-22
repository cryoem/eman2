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
	
	# To make the execution exit upon fatal error by ERROR in global_def.py
	global_def.BATCH = True 
	
	# Check error conditions of arguments
	subcommand_name = "isac_substack"
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
				alpha1  = prealign2d[idx_isac_align2d_alpha]
				sx1     = prealign2d[idx_isac_align2d_tx]
				sy1     = prealign2d[idx_isac_align2d_ty]
				mirror1 = prealign2d[idx_isac_align2d_mirror]
				alpha2  = shrunk_core_align2d[idx_isac_align2d_alpha]
				sx2     = shrunk_core_align2d[idx_isac_align2d_tx]/isac_shrink_ratio # Need to apply the shrink ratio to ISAC x-shift
				sy2     = shrunk_core_align2d[idx_isac_align2d_ty]/isac_shrink_ratio # Need to apply the shrink ratio to ISAC y-shift
				mirror2 = shrunk_core_align2d[idx_isac_align2d_mirror]
				isac_total_align2d = combine_params2(alpha1, sx1, sy1, mirror1, alpha2, sx2, sy2, mirror2)
				
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
			
			fullstack_img_id  = local_total_param2d[idx_beautifier_align2d_fullstack_img_id]
			alpha             = local_total_param2d[idx_beautifier_align2d_alpha]
			sx                = local_total_param2d[idx_beautifier_align2d_tx]
			sy                = local_total_param2d[idx_beautifier_align2d_ty]
			mirror            = local_total_param2d[idx_beautifier_align2d_mirror]
			scale             = local_total_param2d[idx_beautifier_align2d_scale]
			
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
	subcommand_name = "organize_micrographs"
	
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
		print_progress("Creating the destination directory (%)..."%(dst_dir))
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
	subparsers = parser.add_subparsers(help="sub-command help")	

	# create the parser for the "isac_substack" command
	parser_isac_subset = subparsers.add_parser("isac_substack", help="Create Stack Subset: Create virtual subset stack consisting from ISAC accounted particles by retrieving particle numbers associated with the ISAC or Beautifier class averages. The command also saves a list text file containing the retrieved original image numbers and 2D alingment parameters. In addition, it stores the 2D alingment parameters to stack header.")
	parser_isac_subset.add_argument("input_bdb_stack_path",          type=str,                            help="Input BDB image stack: Specify the same BDB image stack used for the associated ISAC run. (default required string)")
	parser_isac_subset.add_argument("input_run_dir",                 type=str,                            help="ISAC or Beautifier run output directory: Specify output directory of an ISAC or Beautifier run as an input to this command. From this directory, the program extracts the shrink ratio and 2D alingment parameters of the ISAC run or local 2D alingment parameters of the Beautifier run. (default required string)")
	parser_isac_subset.add_argument("output_directory",              type=str,                            help="Output directory: The results will be written here. This directory will be created automatically and it must not exist previously. (default required string)")
	parser_isac_subset.add_argument("--isac_class_avgs_path",        type=str,  default="",               help="ISAC or Beautifier class averages path: Specify path to a file containg ISAC or Beautifier class averages. The calss averages can be fullset or selected subset, as long as they are associated with the input BDB image stack and contain class member information stored in the headers. By default, the program uses the same deafult name of ordered class averages in ISAC or Beautifier (i.e. ordered_class_averages.hdf). (default none)")
	parser_isac_subset.add_argument("--substack_basename",           type=str,  default="isac_substack",  help="Substack basename: Specify the basename of ISAC substack file. It cannot be empty string or only white spaces. (default isac_substack)")
	### 
	### NOTE: Toshio Moriya 2018/01/13
	### The following options are not implemented yet.
	### parser_isac_subset.add_argument("--isac_class_id",               type=int,             default=-1,     help="ISAC class average ID: Retrieve only particle members of the specifed ISAC class. By default, retrieve from all classes. (default -1)")
	### parser_isac_subset.add_argument("--no_virtual_stack",            action="store_true",  default=False,  help="Do not create virtual stack: Use this option to create only the particle ID list text file associated with the ISAC class averages. (default False)")
	### parser_isac_subset.add_argument("--no_import_align2d",           action="store_true",  default=False,  help="Do not import alignment:  (default False)")
	parser_isac_subset.set_defaults(func=isac_substack)
	
	# create the parser for the "organize_micrographs" command
	parser_organize_micrographs = subparsers.add_parser("organize_micrographs", help="Organize micrographs: Organize micrographs by moving micrographs in a selecting file from a source directory (specified by source micrographs pattern) to a destination directory.")
	parser_organize_micrographs.add_argument("source_micrograph_pattern",    type=str,                                help="Source micrograph path pattern: Specify path pattern of source micrographs with a wild card (*). Use the wild card to indicate the place of variable part of the file names (e.g. serial number, time stamp, and etc). The path pattern must be enclosed by single quotes (\') or double quotes (\"). (Note: sxgui.py automatically adds single quotes (\')). The substring at the variable part must be same between each associated pair of micrograph names. bdb files can not be selected as source micrographs. (default required string)")
	parser_organize_micrographs.add_argument("selection_list",               type=str,                                help="Micrograph selecting list: Specify a name of text file containing a list of selected micrograph names or paths. The file extension must be \'.txt\'. The directory path of each entry will be ignored if there are any. (default required string)")
	parser_organize_micrographs.add_argument("destination_directory",        type=str,                                help="Destination directory: The micrographs in selecting list will be moved to this directory. This directory will be created automatically if it does not exist. (default required string)")
	parser_organize_micrographs.add_argument("--reverse",                    action="store_true",  default=False,     help="Reverse operation: Move back micrographs from the destination directory to the source directory. Please use this option to restore the previously-moved micrographs. (default False)")
	parser_organize_micrographs.add_argument("--check_consistency",          action="store_true",  default=False,     help="Check consistency of dataset: Create a text file containing the list of Micrograph ID entries might have inconsitency among the provided dataset. (i.e. mic_consistency_check_info_TIMESTAMP.txt). (default False)")
	parser_organize_micrographs.set_defaults(func=organize_micrographs)
	
	# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	# Run specified subcommand
	# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	args = parser.parse_args() # Get namespace object from parser
	# args_dict = vars(parser.parse_args()) # convert it to dictionary object
	# print (args_dict)
	
	print(" ")
	print_progress(get_cmd_line())
	print(" ")
	
	# Call the associated function of the specified subcommand
	args.func(args)

	print(" ")
	print_progress("DONE!!!")
	print(" ")

# ----------------------------------------------------------------------------------------
if __name__ == "__main__":
	main()

# ========================================================================================
# END OF SCRIPT
# ========================================================================================
