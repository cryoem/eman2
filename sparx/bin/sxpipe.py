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
# cd /home/moriya/mrk_develop/sxdemo/sxdemo07_20160908/mpi_bdb_ctf
# rm -r mrkout_sxpipe_isac_substack; sxpipe.py isac_substack bdb:beta20161216_pa03b_sxwindow01#data beta20161216_pa04a_sxisac01/class_averages.hdf mrkout_sxpipe_isac_substack
# rm -r mrkout_sxpipe_isac_substack_02; sxpipe.py isac_substack bdb:beta20161216_pa03b_sxwindow01#data beta20161216_pa04a_sxisac01/class_averages.hdf mrkout_sxpipe_isac_substack_02 --substack_basename=mrk_substack
# 
# ----------------------------------------------------------------------------------------
def isac_substack(args):
	from utilities import get_im, write_text_file
	from EMAN2db import db_open_dict
	from e2bdb import makerelpath
	
	# To make the execution exit upon fatal error by ERROR in global_def.py
	global_def.BATCH = True 
	
	# Check error conditions
	subcommand_name = "isac_substack"
	if not os.path.exists(args.input_isac_class_avgs_path):
		ERROR("Input ISAC class average stack file does not exist. Please check the file path and restart the program.", subcommand_name) # action=1 - fatal error, exit
	if os.path.exists(args.output_directory):
		ERROR("Output directory exists. Please change the name and restart the program.", subcommand_name) # action=1 - fatal error, exit
	if args.substack_basename.strip() == "":
		ERROR("Substack basename cannot be empty string or only white spaces.", subcommand_name) # action=1 - fatal error, exit
	
	assert (os.path.exists(args.input_isac_class_avgs_path))
	assert (not os.path.exists(args.output_directory))
	assert (args.substack_basename.strip() != "")
	
	# Create output directory
	os.mkdir(args.output_directory)
	
	# Retrieve original particle IDs of member particles listed in ISAC class average stack
	n_img_processed = EMUtil.get_image_count(args.input_isac_class_avgs_path)
	isac_substack_particle_id_list = []
	for i_img in xrange(n_img_processed):
		isac_substack_particle_id_list += get_im(args.input_isac_class_avgs_path, i_img).get_attr("members")
	isac_substack_particle_id_list.sort()
	
	# Save the substack particle id list
	isac_substack_particle_id_list_file_path = os.path.join(args.output_directory, "{0}_particle_id_list.txt".format(args.substack_basename))
	write_text_file(isac_substack_particle_id_list, isac_substack_particle_id_list_file_path)
	
	# Open the output BDB dictionary
	assert (args.output_directory != "")
	output_virtual_bdb_stack_real_path = "bdb:{0}#{1}".format(args.output_directory,args.substack_basename )
	output_virtual_bdb_stack = db_open_dict(output_virtual_bdb_stack_real_path)
	
	# Convert an absolute path to the actual output data to a relative path by eliminating any symbolic links 
	output_virtual_bdb_stack_real_path=os.path.realpath(output_virtual_bdb_stack.path)+"/"
		
	# Open the input BDB dictionary
	input_bdb_stack = db_open_dict(args.input_bdb_stack_path, ro=True) # Read only
	
	# Copy the header from input to output BDB dictionary
	n_img_detected = len(isac_substack_particle_id_list)
	print(" ")
	print_progress("Detected %d ISAC validated particles in %s"%(n_img_detected, args.input_isac_class_avgs_path))
	
	# Loop through all ISAC validated particles
	print(" ")
	n_img_processed = 0
	n_img_of_10_percent = n_img_detected // 10
	for i_img_detected, isac_substack_particle_id in enumerate(isac_substack_particle_id_list):
		# Print progress
		if i_img_detected % n_img_of_10_percent == 0:
			try:
				print_progress("Progress %5.2f%%: Processing %6dth entry (Particle ID %6d)."%(float(i_img_detected)/n_img_detected*100.0, i_img_detected, isac_substack_particle_id))
				sys.stdout.flush()
			except:
				pass
		
		# Read a particle image header from input bdb stack
		try: 
			img_header = input_bdb_stack.get(isac_substack_particle_id, nodata=1).get_attr_dict() # Need only header information
		except:
			ERROR("Failed to read image header of particle #%d from %s. Skipping this image..."%(isac_substack_particle_id, args.input_bdb_stack_path), subcommand_name, action = 0) # action = 0 - non-fatal, print a warning;
			continue
		
		# Convert an absolute path to the actual input data to a relative path by eliminating any symbolic links 
		try:
			input_bdb_stack_real_path = os.path.realpath(input_bdb_stack.get_data_path(isac_substack_particle_id))
			# Conver the path to OS specific format
			if os.name == "nt":
				output_virtual_bdb_stack_real_path = output_virtual_bdb_stack_real_path.replace("\\", "/")
				input_bdb_stack_real_path = input_bdb_stack_real_path.replace("\\", "/")
			# Takes a pair of paths /a/b/c/d and /a/b/e/f/g and returns a relative path to b from a, ../../e/f/g
			common_relative_path = makerelpath(output_virtual_bdb_stack_real_path, input_bdb_stack_real_path)
		except:
			ERROR("Failure to find common relative data path for particle image #%d. Skipping this image..."%(isac_substack_particle_id), subcommand_name, action = 0) # action = 0 - non-fatal, print a warning;
			continue
		
		# Update the image header for output
		img_header["data_path"]    = common_relative_path
		img_header["data_n"]       = isac_substack_particle_id
		img_header["data_source"]  = args.input_bdb_stack_path
		
		# Register the image header to output virtual bdb stack
		output_virtual_bdb_stack[n_img_processed] = img_header
		
		# Increment process image counts
		n_img_processed += 1
	
	# Close input and output bdb stacks
	output_virtual_bdb_stack.close()
	input_bdb_stack.close()

	# Print summary of processing
	print(" ")
	print_progress("Summary of processing...")
	print_progress("Detected  : %6d"%(n_img_detected))
	print_progress("Processed : %6d"%(n_img_processed))
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
	parser_isac_subset = subparsers.add_parser("isac_substack", help="Create Stack Subset: Create virtual subset stack consisting from ISAC accounted particles by retrieving particle numbers associated with the class averages. The command also saves a list text file containing the retrieved original image numbers.")
	parser_isac_subset.add_argument("input_bdb_stack_path",          type=str,                              help="Input BDB image stack: Specify the same BDB image stack used for the associated ISAC run. (default required string)")
	parser_isac_subset.add_argument("input_isac_class_avgs_path",    type=str,                              help="ISAC class average file path: Input ISAC class average file path. (default required string)")
	parser_isac_subset.add_argument("output_directory",              type=str,                              help="Output directory: The results will be written here. This directory will be created automatically and it must not exist previously. (default required string)")
	parser_isac_subset.add_argument("--substack_basename",           type=str,  default="isac_substack",    help="Substack basename: Specify the basename of ISAC substack file.  It cannot be empty string or only white spaces. (default isac_substack)")
	### 
	### NOTE: Toshio Moriya 2017/11/16
	### The following options are not implemented yet.
	### parser_isac_subset.add_argument("--isac_class_id",               type=int,             default=-1,     help="ISAC class average ID: Retrieve only particle members of the specifed ISAC class. By default, retrieve from all classes. (default -1)")
	### parser_isac_subset.add_argument("--no_virtual_stack",            action="store_true",  default=False,  help="Do not create virtual stack: Use this option to create only the particle ID list text file associated with the ISAC class averages. (default False)")
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
