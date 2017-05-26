#!/usr/bin/env python
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
import	global_def
from	global_def 	import *

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
	subcommand_name = 'isac_substack'
	if not os.path.exists(args.input_isac_class_avgs_path):
		ERROR('Input ISAC class average stack file does not exist. Please check the file path and restart the program.', subcommand_name) # action=1 - fatal error, exit
	if os.path.exists(args.output_directory):
		ERROR('Output directory exists. Please change the name and restart the program.', subcommand_name) # action=1 - fatal error, exit
	
	assert(os.path.exists(args.input_isac_class_avgs_path))
	assert(not os.path.exists(args.output_directory))
	
	# Create output directory
	os.mkdir(args.output_directory)
	
	# Retrieve original particle IDs of member particles listed in ISAC class average stack
	n_img_processed = EMUtil.get_image_count(args.input_isac_class_avgs_path)
	isac_substack_particle_id_list = []
	for i_img in xrange(n_img_processed):
		isac_substack_particle_id_list += get_im(args.input_isac_class_avgs_path, i_img).get_attr('members')
	isac_substack_particle_id_list.sort()
	
	# Save the substack particle id list
	isac_substack_particle_id_list_file_path = os.path.join(args.output_directory, 'isac_substack_particle_id_list.txt')
	write_text_file(isac_substack_particle_id_list, isac_substack_particle_id_list_file_path)
	
	# Open the output BDB dictionary
	assert(args.output_directory != '')
	output_virtual_bdb_stack_real_path = 'bdb:%s#isac_substack' % args.output_directory
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
	for i_img_detected, isac_substack_particle_id in enumerate(isac_substack_particle_id_list):
		# Print progress
		if i_img_detected % 1000 == 0:
			try:
				print_progress("Progress %5.2f%%: Processing %6dth entry (Particle ID %6d)."%(float(i_img_detected)/n_img_detected*100.0, i_img_detected, isac_substack_particle_id))
				sys.stdout.flush()
			except:
				pass
		
		# Read a particle image header from input bdb stack
		try: 
			img_header = input_bdb_stack.get(isac_substack_particle_id, nodata=1).get_attr_dict() # Need only header information
		except:
			ERROR('Failed to read image header of particle #%d from %s. Skipping this image...' % (isac_substack_particle_id, args.input_bdb_stack_path), subcommand_name, action = 0) # action = 0 - non-fatal, print a warning;
			continue
		
		# Convert an absolute path to the actual input data to a relative path by eliminating any symbolic links 
		try:
			input_bdb_stack_real_path = os.path.realpath(input_bdb_stack.get_data_path(isac_substack_particle_id))
			# Conver the path to OS specific format
			if os.name == 'nt':
				output_virtual_bdb_stack_real_path = output_virtual_bdb_stack_real_path.replace("\\", '/')
				input_bdb_stack_real_path = input_bdb_stack_real_path.replace('\\', '/')
			# Takes a pair of paths /a/b/c/d and /a/b/e/f/g and returns a relative path to b from a, ../../e/f/g
			common_relative_path = makerelpath(output_virtual_bdb_stack_real_path, input_bdb_stack_real_path)
		except:
			ERROR('Failure to find common relative data path for particle image #%d. Skipping this image...' % (isac_substack_particle_id), subcommand_name, action = 0) # action = 0 - non-fatal, print a warning;
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

# ========================================================================================
# Main function
# ========================================================================================
def main():
	# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	# Set up argument parser (supports subcommand)
	# ><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><
	parser = argparse.ArgumentParser(description='The collection of SPHIRE small pipleline tools.')
	parser.add_argument('--version', action='version', version=SPARXVERSION)
	# subparsers = parser.add_subparsers(title='subcommands', description='valid subcommands', help='additional help')	
	subparsers = parser.add_subparsers(help='sub-command help')	

	# create the parser for the 'isac_substack' command
	parser_isac_subset = subparsers.add_parser('isac_substack', help='Create Stack Subset: Create virtual subset stack consisting from ISAC accounted particles by retrieving particle numbers associated with the class averages. The command also saves a list text file containing the retrieved original image numbers.')
	parser_isac_subset.add_argument('input_bdb_stack_path',        type=str,                             help='Input BDB image stack: Specify the same BDB image stack used for the associated ISAC run. (default required string)')
	parser_isac_subset.add_argument('input_isac_class_avgs_path',  type=str,                             help='ISAC class average file path: Input ISAC class average file path. (default required string)')
	parser_isac_subset.add_argument('output_directory',            type=str,                             help='Output directory: The results will be written here. This directory will be created automatically and it must not exist previously. (default required string)')
	parser_isac_subset.add_argument('--isac_class_id',             type=int,             default=-1,     help='ISAC class average ID: Retrieve only particle members of the specifed ISAC class. By default, retrieve from all classes. (default -1)')
	parser_isac_subset.add_argument('--no_virtual_stack',          action="store_true",  default=False,  help='Do not create virtual stack: Use this option to create only the particle ID list text file associated with the ISAC class averages. (default False)')
	parser_isac_subset.set_defaults(func=isac_substack)
	
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
if __name__ == '__main__':
	main()

# ========================================================================================
# END OF SCRIPT
# ========================================================================================
