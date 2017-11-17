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
from __future__ import print_function
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

# ========================================================================================
# Subcommand functions
# ========================================================================================
# ----------------------------------------------------------------------------------------
# TEST COMMAND
# cd /home/moriya/mrk_develop/sxdemo/sxdemo07_20160908/mpi_bdb_ctf
# rm -r mrkout_sxpipe_isac_substack; sxpipe.py isac_substack bdb:beta20161216_pa03b_sxwindow01#data beta20161216_pa04a_sxisac01/class_averages.hdf mrkout_sxpipe_isac_substack
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
	
	assert (os.path.exists(args.input_isac_class_avgs_path))
	assert (not os.path.exists(args.output_directory))
	
	# Create output directory
	os.mkdir(args.output_directory)
	
	# Retrieve original particle IDs of member particles listed in ISAC class average stack
	n_img_processed = EMUtil.get_image_count(args.input_isac_class_avgs_path)
	isac_substack_particle_id_list = []
	for i_img in xrange(n_img_processed):
		isac_substack_particle_id_list += get_im(args.input_isac_class_avgs_path, i_img).get_attr("members")
	isac_substack_particle_id_list.sort()
	
	# Save the substack particle id list
	isac_substack_particle_id_list_file_path = os.path.join(args.output_directory, "isac_substack_particle_id_list.txt")
	write_text_file(isac_substack_particle_id_list, isac_substack_particle_id_list_file_path)
	
	# Open the output BDB dictionary
	assert (args.output_directory != "")
	output_virtual_bdb_stack_real_path = "bdb:%s#isac_substack"%(args.output_directory)
	output_virtual_bdb_stack = db_open_dict(output_virtual_bdb_stack_real_path)
	
	# Convert an absolute path to the actual output data to a relative path by eliminating any symbolic links 
	output_virtual_bdb_stack_real_path=os.path.realpath(output_virtual_bdb_stack.path)+"/"
		
	# Open the input BDB dictionary
	input_bdb_stack = db_open_dict(args.input_bdb_stack_path, ro=True) # Read only
	
	# Copy the header from input to output BDB dictionary
	n_img_detected = len(isac_substack_particle_id_list)
	print_progress("Detected %d ISAC validated particles in %s"%(n_img_detected, args.input_isac_class_avgs_path))
	print(" ")
	
	# Loop through all ISAC validated particles
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
# sxpipe.py organize_micrographs 'CorrectedSums/corrsum/TcdA1-*_frames_sum.mrc' 'CorrectedSums/discard_corrsum' 'CTFest/Tutorial_micrographs_discard.txt' --check_consistency
# 
# ----------------------------------------------------------------------------------------
def organize_micrographs(args):
	import glob
	from utilities import read_text_file
	
	# To make the execution exit upon fatal error by ERROR in global_def.py
	global_def.BATCH = True 
	
	# ------------------------------------------------------------------------------------
	# Prepare the variables for all sections
	# ------------------------------------------------------------------------------------
	# Use short names for arguments and options
	input_mic_pattern = args.input_micrograph_pattern
	output_dir = args.output_directory
	
	# Check error conditions
	subcommand_name = "move_micrographs"
	if input_mic_pattern.find("*") == -1:
		ERROR("Input micrograph path pattern must contain wild card (*). Please check input_micrograph_pattern argument. Please correct input_micrograph_pattern argument and restart the program.", subcommand_name) # action=1 - fatal error, exit
	
	if not os.path.exists(args.selection_list):
		ERROR("Input micrograph selecting list file does not exist. Please correct the file path and restart the program.", subcommand_name) # action=1 - fatal error, exit
	assert (os.path.exists(args.selection_list))

	# --------------------------------------------------------------------------------
	# Prepare variables
	# --------------------------------------------------------------------------------
	# Define indices of selection list parameters
	i_enum = -1
	i_enum += 1; idx_mic_list_mic_path   = i_enum # The name or path of micrographs
	i_enum += 1; n_idx_mic_list          = i_enum

	# Global entry dictionary (all possible entries from all lists) for all mic id substring
	global_entry_dict = {} # mic id substring is the key
	subkey_input_mic_path = "Input Micrograph Path"
	subkey_output_mic_path = "Output Micrograph Path"
	subkey_selected_mic_basename = "Selected Micrograph Basename"
	
	# List keeps only id substrings of micrographs whose all necessary information are available
	valid_mic_id_substr_list = [] 

	# Micrograph basename pattern (directory path is removed from micrograph path pattern)
	mic_basename_pattern = os.path.basename(input_mic_pattern)
	
	# Prefix and suffix of micrograph basename pattern 
	# to find the head/tail indices of micrograph id substring
	mic_basename_tokens = mic_basename_pattern.split("*")
	assert (len(mic_basename_tokens) == 2)
	# Find head index of micrograph id substring
	mic_id_substr_head_idx = len(mic_basename_tokens[0])

	# Set up output directory 
	output_mic_pattern = None
	if os.path.exists(output_dir):
		print(" ")
		print_progress("Output directory (%s) exists. The program will move duplicated micrographs in output directoy to \'duplicates_DATE_TIME\' directory if there are any."%(output_dir))
		output_mic_pattern = os.path.join(output_dir, mic_basename_pattern)
	
	# --------------------------------------------------------------------------------
	# Register micrograph id substrings found in input directory (specified by micrograph path pattern)
	# and associated input micrograph path to the global entry dictionary
	# --------------------------------------------------------------------------------
	# Generate the list of micrograph paths in the input directory
	print(" ")
	print_progress("Checking the input directory...")
	input_mic_path_list = glob.glob(input_mic_pattern)
	# Check error condition of input micrograph file path list
	print("Found %d microgarphs in %s."%(len(input_mic_path_list), os.path.dirname(input_mic_pattern)))
	if len(input_mic_path_list) == 0:
		ERROR("No micrograph files are found in the directory specified by micrograph path pattern (%s). Please check input_micrograph_pattern argument and restart the program."%(os.path.dirname(input_mic_pattern))) # action=1 - fatal error, exit
	assert (len(input_mic_path_list) > 0)
	
	# Register micrograph id substrings to the global entry dictionary
	for input_mic_path in input_mic_path_list:
		# Find tail index of micrograph id substring and extract the substring from the micrograph name
		input_mic_basename = os.path.basename(input_mic_path)
		mic_id_substr_tail_idx = input_mic_basename.index(mic_basename_tokens[1])
		mic_id_substr = input_mic_basename[mic_id_substr_head_idx:mic_id_substr_tail_idx]
		assert (input_mic_path == input_mic_pattern.replace("*", mic_id_substr))
		if not mic_id_substr in global_entry_dict:
			# print("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from input_mic_path_list "%(mic_id_substr))
			global_entry_dict[mic_id_substr] = {}
		assert (mic_id_substr in global_entry_dict)
		global_entry_dict[mic_id_substr][subkey_input_mic_path] = input_mic_path
	assert (len(global_entry_dict) > 0)
	
	# Clean up variables which won't be used anymore
	del input_mic_path_list

	# --------------------------------------------------------------------------------
	# Register micrograph id substrings found in output directory if any
	# and associated input micrograph path to the global entry dictionary
	# --------------------------------------------------------------------------------
	if output_mic_pattern is not None:
		assert (os.path.exists(output_dir))
		output_mic_pattern = os.path.join(output_dir, mic_basename_pattern)
		# Generate the list of micrograph paths in the output directory
		print(" ")
		print_progress("Checking the output directory...")
		output_mic_path_list = glob.glob(output_mic_pattern)
		# Check error condition of input micrograph file path list
		print("Found %d microgarphs in %s."%(len(output_mic_path_list), output_dir))
		
		# Register micrograph id substrings to the global entry dictionary
		for output_mic_path in output_mic_path_list:
			# Find tail index of micrograph id substring and extract the substring from the micrograph name
			output_mic_basename = os.path.basename(output_mic_path)
			mic_id_substr_tail_idx = output_mic_basename.index(mic_basename_tokens[1])
			mic_id_substr = output_mic_basename[mic_id_substr_head_idx:mic_id_substr_tail_idx]
			assert (output_mic_path == output_mic_pattern.replace("*", mic_id_substr))
			if not mic_id_substr in global_entry_dict:
				# print("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from output_mic_path_list "%(mic_id_substr))
				global_entry_dict[mic_id_substr] = {}
			assert (mic_id_substr in global_entry_dict)
			global_entry_dict[mic_id_substr][subkey_output_mic_path] = output_mic_path
		assert (len(global_entry_dict) > 0)
	
		# Clean up variables which won't be used anymore
		del output_mic_path_list

	# --------------------------------------------------------------------------------
	# Register micrograph id substrings found in the selection list
	# and associated micrograph basename to the global entry dictionary
	# --------------------------------------------------------------------------------
	# Generate the list of selected micrograph paths in the selection file
	selected_mic_path_list = []
	# Generate micrograph lists according to the execution mode
	print("Checking the selection list...")
	selected_mic_path_list = read_text_file(args.selection_list)
	
	# Check error condition of micrograph entry lists
	print("Found %d microgarph entries in %s."%(len(selected_mic_path_list), args.selection_list))
	if len(selected_mic_path_list) == 0:
		ERROR("No micrograph entries are found in the selection list file. Please check selection_list option and restart the program."%(os.path.dirname(input_mic_pattern))) # action=1 - fatal error, exit
	assert (len(selected_mic_path_list) > 0)
	
	selected_mic_directory = os.path.dirname(selected_mic_path_list[0])
	if selected_mic_directory != "":
		print("    NOTE: Program disregards the directory paths in the selection list (%s)."%(selected_mic_directory))

	# Register micrograph id substrings to the global entry dictionary
	for selected_mic_path in selected_mic_path_list:
		# Find tail index of micrograph id substring and extract the substring from the micrograph name
		selected_mic_basename = os.path.basename(selected_mic_path)
		mic_id_substr_tail_idx = selected_mic_basename.index(mic_basename_tokens[1])
		mic_id_substr = selected_mic_basename[mic_id_substr_head_idx:mic_id_substr_tail_idx]
		assert (selected_mic_basename == mic_basename_pattern.replace("*", mic_id_substr))
		if not mic_id_substr in global_entry_dict:
			# print("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from selected_mic_path_list "%(mic_id_substr))
			global_entry_dict[mic_id_substr] = {}
		assert (mic_id_substr in global_entry_dict)
		global_entry_dict[mic_id_substr][subkey_selected_mic_basename] = selected_mic_basename
	assert (len(global_entry_dict) > 0)
	
	# Clean up variables which won't be used anymore
	del selected_mic_path_list
	
	# --------------------------------------------------------------------------------
	# Clean up variables related to registration to the global entry dictionary
	# --------------------------------------------------------------------------------
	del mic_basename_tokens
	del mic_id_substr_head_idx

	# --------------------------------------------------------------------------------
	# Create the list containing only valid micrograph id substrings
	# --------------------------------------------------------------------------------
	print(" ")
	print("Checking the input datasets consistency...")

	if output_mic_pattern is None:
		assert (not os.path.exists(output_dir))
		# Prepare lists to keep track of invalid (rejected) micrographs
		no_input_mic_id_substr_list = []
		
		# Loop over substring id list
		for mic_id_substr in global_entry_dict:
			mic_id_entry = global_entry_dict[mic_id_substr]
		
			warinnig_messages = []
			# selected micrograph basename must have been registed always .
			if subkey_selected_mic_basename in mic_id_entry: 
				# Check if associated input micrograph exists
				if not subkey_input_mic_path in mic_id_entry:
					input_mic_path = input_mic_pattern.replace("*", mic_id_substr)
					warinnig_messages.append("    associated input micrograph %s does not exist."%(input_mic_path))
					no_input_mic_id_substr_list.append(mic_id_substr)
			
				if len(warinnig_messages) > 0:
					print("WARNING!!! Micrograph ID %s have inconsistency amoung provided information:"%(mic_id_substr))
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
			assert (not os.path.exists(output_dir))
			os.mkdir(output_dir)
			assert (os.path.exists(output_dir))
			
			# Open the consistency check file
			inconsist_mic_list_path = os.path.join(output_dir,"inconsist_mic_id_file.txt")
			print(" ")
			print("Generating the provided information consistency report in %s..."%(inconsist_mic_list_path))
			inconsist_mic_list_file = open(inconsist_mic_list_path, "w")
			inconsist_mic_list_file.write("# The information about inconsistent micrograph IDs\n")
			# Loop over substring id list
			for mic_id_substr in global_entry_dict:
				mic_id_entry = global_entry_dict[mic_id_substr]
			
				consistency_messages = []
				# Check if associated input micrograph path exists
				if not subkey_input_mic_path in mic_id_entry:
					input_mic_path = input_mic_pattern.replace("*", mic_id_substr)
					consistency_messages.append("    associated input micrograph %s is missing."%(input_mic_path))
			
				# Check if associated selected micrograph basename exists
				if not subkey_selected_mic_basename in mic_id_entry:
					input_mic_path = input_mic_pattern.replace("*", mic_id_substr)
					consistency_messages.append("    associated selected micrograph %s is missing."%(input_mic_path))
			
				if len(consistency_messages) > 0:
					inconsist_mic_list_file.write("Micrograph ID %s have inconsistency amoung provided information:\n"%(mic_id_substr))
					for consistency_message in consistency_messages:
						inconsist_mic_list_file.write(consistency_message)
						inconsist_mic_list_file.write("\n")
		
			# Close the consistency check file, if necessary
			inconsist_mic_list_file.flush()
			inconsist_mic_list_file.close()
		
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
		print("Detected                           : %6d"%(len(global_entry_dict)))
		print("Valid                              : %6d"%(len(valid_mic_id_substr_list)))
		print("Rejected by no input micrograph    : %6d"%(len(no_input_mic_id_substr_list)))
		
		# --------------------------------------------------------------------------------
		# Clean up variables related to tracking of invalid (rejected) micrographs 
		# --------------------------------------------------------------------------------
		del no_input_mic_id_substr_list

	else:
		assert (output_mic_pattern is not None)
		assert (os.path.exists(output_dir))
		# Prepare lists to keep track of invalid (rejected) micrographs
		no_mic_in_both_dir_id_substr_list = []
		already_in_output_dir_mic_id_substr_list = []
		duplicated_in_output_dir_mic_id_substr_list = []

		# Loop over substring id list
		for mic_id_substr in global_entry_dict:
			mic_id_entry = global_entry_dict[mic_id_substr]
		
			warinnig_messages = []
			# selected micrograph basename must have been registed always .
			if subkey_selected_mic_basename in mic_id_entry: 
				# Check if associated input micrograph exists
				if not subkey_input_mic_path in mic_id_entry:
					input_mic_path = input_mic_pattern.replace("*", mic_id_substr)
					output_mic_path = output_mic_pattern.replace("*", mic_id_substr)
					if not subkey_output_mic_path in mic_id_entry:
						warinnig_messages.append("    associated micrographs do not exist neither in input and output directories (%s and %s)."%(input_mic_path, output_mic_path))
						no_mic_in_both_dir_id_substr_list.append(mic_id_substr)
					else:
						assert (subkey_output_mic_path in mic_id_entry)
						warinnig_messages.append("    associated micrographs already has been moved to output directory (%s)."%(output_mic_path))
						already_in_output_dir_mic_id_substr_list.append(mic_id_substr)
				else: 
					assert (subkey_input_mic_path in mic_id_entry)
					if subkey_output_mic_path in mic_id_entry:
						input_mic_path = input_mic_pattern.replace("*", mic_id_substr)
						output_mic_path = output_mic_pattern.replace("*", mic_id_substr)
						warinnig_messages.append("    associated micrographs are dupliicated in both input and output directories (%s and %s)."%(input_mic_path, output_mic_path))
						duplicated_in_output_dir_mic_id_substr_list.append(mic_id_substr)
					# else:
					# 	# This is normal case!
					# 	assert (not subkey_output_mic_path in mic_id_entry)
			
				if len(warinnig_messages) > 0:
					print("WARNING!!! Micrograph ID %s have inconsistency amoung provided information:"%(mic_id_substr))
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
			assert (os.path.exists(output_dir))
			
			# Open the consistency check file
			inconsist_mic_list_path = os.path.join(output_dir,"inconsist_mic_id_file.txt")
			print(" ")
			print("Generating the provided information consistency report in %s..."%(inconsist_mic_list_path))
			inconsist_mic_list_file = open(inconsist_mic_list_path, "w")
			inconsist_mic_list_file.write("# The information about inconsistent micrograph IDs\n")
			# Loop over substring id list
			for mic_id_substr in global_entry_dict:
				mic_id_entry = global_entry_dict[mic_id_substr]
				
				consistency_messages = []
				# Check if associated input micrograph path exists
				if not subkey_input_mic_path in mic_id_entry:
					input_mic_path = input_mic_pattern.replace("*", mic_id_substr)
					consistency_messages.append("    associated micrograph %s does not exist in input directory."%(input_mic_path))
			
				if subkey_output_mic_path in mic_id_entry:
					output_mic_path = output_mic_pattern.replace("*", mic_id_substr)
					consistency_messages.append("    associated micrograph %s already exist in output directory."%(output_mic_path))
			
				# Check if associated selected micrograph basename exists
				if not subkey_selected_mic_basename in mic_id_entry:
					input_mic_path = input_mic_pattern.replace("*", mic_id_substr)
					consistency_messages.append("    associated selected micrograph %s is missing. This must be discared."%(input_mic_path))
			
				if len(consistency_messages) > 0:
					inconsist_mic_list_file.write("Micrograph ID %s have inconsistency amoung provided information:\n"%(mic_id_substr))
					for consistency_message in consistency_messages:
						inconsist_mic_list_file.write(consistency_message)
						inconsist_mic_list_file.write("\n")
		
			# Close the consistency check file, if necessary
			inconsist_mic_list_file.flush()
			inconsist_mic_list_file.close()
		
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
		print("Detected                           : %6d"%(len(global_entry_dict)))
		print("Valid                              : %6d"%(len(valid_mic_id_substr_list)))
		print("Rejected by not found at all       : %6d"%(len(no_mic_in_both_dir_id_substr_list)))
		print("Rejected by already in output dir  : %6d"%(len(already_in_output_dir_mic_id_substr_list)))
		print("Rejected by duplicated             : %6d"%(len(duplicated_in_output_dir_mic_id_substr_list)))
		
		# --------------------------------------------------------------------------------
		# Clean up variables related to tracking of invalid (rejected) micrographs 
		# --------------------------------------------------------------------------------
		del no_mic_in_both_dir_id_substr_list
		del already_in_output_dir_mic_id_substr_list
		del duplicated_in_output_dir_mic_id_substr_list

"""
	# --------------------------------------------------------------------------------
	# Create output directory
	# --------------------------------------------------------------------------------
	print_progress("Creating output directory...")
	os.mkdir(args.output_directory)
	print(" ")
	
	# --------------------------------------------------------------------------------
	# Load micrograph list form the file
	# --------------------------------------------------------------------------------
	print_progress("Loading micrograph name list from the selection list file...")
	assert (os.path.exists(args.selection_list))
	mic_path_list = read_text_file(args.selection_list)
	
	print_progress("# Found %d entries in %s"%(len(mic_path_list), args.selection_list))
	print(" ")

	# move micrographs in the list to output directory
	for mic_path in mic_path_list:
		if os.path.exists(mic_path):
			print_progress("Moving %s..."%(mic_path))
			shutil.move(mic_path, args.output_directory)
		else:
			assert (os.path.exists(mic_path) == False)
			print_progress "WARNINNG: Cannot find %s!!! Skipping this micrograph..."%(mic_path)
			continue
	
	# Print summary of processing
	print(" ")
	print_progress("Summary of processing...")
	print_progress("Detected  : %6d"%(n_img_detected))
	print_progress("Processed : %6d"%(n_img_processed))
	print(" ")
"""
	
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
	parser_isac_subset.add_argument("input_bdb_stack_path",          type=str,                             help="Input BDB image stack: Specify the same BDB image stack used for the associated ISAC run. (default required string)")
	parser_isac_subset.add_argument("input_isac_class_avgs_path",    type=str,                             help="ISAC class average file path: Input ISAC class average file path. (default required string)")
	parser_isac_subset.add_argument("output_directory",              type=str,                             help="Output directory: The results will be written here. This directory will be created automatically and it must not exist previously. (default required string)")
	### 
	### NOTE: Toshio Moriya 2017/11/16
	### The following options are not implemented yet.
	### parser_isac_subset.add_argument("--isac_class_id",               type=int,             default=-1,     help="ISAC class average ID: Retrieve only particle members of the specifed ISAC class. By default, retrieve from all classes. (default -1)")
	### parser_isac_subset.add_argument("--no_virtual_stack",            action="store_true",  default=False,  help="Do not create virtual stack: Use this option to create only the particle ID list text file associated with the ISAC class averages. (default False)")
	parser_isac_subset.set_defaults(func=isac_substack)
	
	# create the parser for the "move_micrographs" command
	parser_move_micrographs = subparsers.add_parser("organize_micrographs", help="Organize micrographs: Organize micrographs by moving micrographs in a selecting file from specifed input directory (given by input micrographs pattern) to specified output directory.")
	parser_move_micrographs.add_argument("input_micrograph_pattern", type=str,                             help="Input micrograph path pattern: Specify path pattern of input micrographs with a wild card (*). Use the wild card to indicate the place of variable part of the file names (e.g. serial number, time stamp, and etc). The path pattern must be enclosed by single quotes (\') or double quotes (\"). (Note: sxgui.py automatically adds single quotes (\')). The substring at the variable part must be same between the associated pair of input micrograph and coordinates file. bdb files can not be selected as input micrographs. (default required string)")
	parser_move_micrographs.add_argument("output_directory",         type=str,                             help="Output directory: The results will be written here. This directory will be created automatically and it must not exist previously. (default required string)")
	parser_move_micrographs.add_argument("selection_list",           type=str,                             help="Micrograph selecting list: Specify a name of text file containing a list of selected micrograph names or paths. The directory path of each entry will be ignored if there is. (default required string)")
	parser_move_micrographs.add_argument("--check_consistency",      action="store_true",  default=False,  help="Check consistency of inputs: Create a text file containing the list of inconsistent Micrograph ID entries (i.e. inconsist_mic_list_file.txt). (default False)")
	parser_move_micrographs.set_defaults(func=organize_micrographs)
	
	# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	# Run specified subcommand
	# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	args = parser.parse_args() # Get namespace object from parser
	# args_dict = vars(parser.parse_args()) # convert it to dictionary object
	# print (args_dict)
	
	print_progress(get_cmd_line())
	print(" ")
	
	# Call the associated function of the specified subcommand
	args.func(args)

	print_progress("DONE!!!")
	print(" ")

# ----------------------------------------------------------------------------------------
if __name__ == "__main__":
	main()

# ========================================================================================
# END OF SCRIPT
# ========================================================================================
