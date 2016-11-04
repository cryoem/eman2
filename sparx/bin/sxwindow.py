#!/usr/bin/env python
#
# Author: Toshio Moriya 10/20/2016 (toshio.moriya@mpi-dortmund.mpg.de)
# Author: T. Durmaz 08/29/2014 (tunay.durmaz@uth.tmc.edu)
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPHIRE software packages have some GPL dependencies,
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

from __future__ import print_function
import os, sys
from optparse import OptionParser, SUPPRESS_HELP
import glob
import json

from EMAN2 import *
from EMAN2db import *
from EMAN2jsondb import *
from sparx import *
from applications import MPI_start_end
from inspect import currentframe, getframeinfo

# ========================================================================================
# Define functions for reading coordinates files of different formats.
# One of these will be used in main() through a first-class data type variable of Python 
# (like function pointer in C/C++)
# This way, switch statement is unnecessary inside of the coordinates loop.
# ========================================================================================
def read_sparx_coords_file(coords_path):
	coords_list = read_text_row(coords_path)
	return coords_list

def read_eman1_coords_file(coords_path):
	coords_list = read_text_row(coords_path)
	for i in xrange(len(coords_list)):
		coords_list[i] = [(coords_list[i][0] + coords_list[i][2] // 2), (coords_list[i][1] + coords_list[i][3] // 2)]
	return coords_list

def read_eman2_coords_file(coords_path):
	coords_list = js_open_dict(coords_path)["boxes"]
	for i in xrange(len(coords_list)):
		coords_list[i] = [coords_list[i][0], coords_list[i][1]]
	return coords_list

def read_spider_coords_file(coords_path):
	coords_list = read_text_row(coords_path)
	for i in xrange(len(coords_list)):
		coords_list[i] = [coords_list[i][2], coords_list[i][3]]
	return coords_list

# ========================================================================================
#  Helper functions
# ========================================================================================
def check_options(options, program_name):

	error_status = None
	# not a real while, an if with the opportunity to use break when errors need to be reported
	while True:
		if options.selection_list != None:
			if not os.path.exists(options.selection_list): 
				error_status = ("File specified by selection_list option does not exists. Please check selection_list option. Run %s -h for help." % (program_name), getframeinfo(currentframe()))
				break
		
		if options.coordinates_format.lower() not in ["sparx", "eman1", "eman2", "spider"]:
			error_status = ("Invalid option value: --coordinates_format=%s. Please run %s -h for help." % (options.coordinates_format, program_name), getframeinfo(currentframe()))
			break
		
		if options.import_ctf:
			if os.path.exists(options.import_ctf) == False:
				error_status = ("Specified CTER CTF file is not found. Please check --import_ctf option. Run %s -h for help." % (program_name), getframeinfo(currentframe()))
				break
		else:
			assert (not options.import_ctf)
			if options.limit_ctf:
				error_status = ("--limit_ctf option requires valid CTER CTF File (--import_ctf). Please run %s -h for help." % (program_name), getframeinfo(currentframe()))
				break
		
		if (options.resample_ratio <= 0.0 or options.resample_ratio > 1.0):
			error_status = ("Invalid option value: --resample_ratio=%s. Please run %s -h for help." % (options.resample_ratio, program_name), getframeinfo(currentframe()))
			break
		
		break
	if_error_then_all_processes_exit_program(error_status)
	
# ========================================================================================
#  Main function
# ========================================================================================
def main():
	program_name = os.path.basename(sys.argv[0])
	usage = program_name + """  input_micrograph_pattern  input_coordinates_pattern  output_directory  --selection_list=selection_list  --coordinates_format  --box_size=box_size  --skip_invert  --import_ctf=ctf_file  --limit_ctf  --astigmatism_error=astigmatism_error  --resample_ratio=resample_ratio  --check_consistency
	
Window particles from micrographs using the particles coordinates.

All Micrographs Mode - Process all micrographs in a directory:
	Specify path pattern of input micrographs and coordinates files with a wild card (*). 
	Use the wild card to indicate the place of variable part of the file names (e.g. serial number, time stamp, and etc). 
	The path pattern must be enclosed by single quotes (') or double quotes ("). (Note: sxgui.py automatically adds single quotes (')). 
	The substring at the variable part must be same between the associated pair of input micrograph and coordinates file.
	bdb files can not be selected as input micrographs.
	In this mode, all micrographs matching the path pattern will be processed.
	
	mpirun  -np  32  sxwindow.py  './mic*.hdf'  'info/mic*_info.json'  particles  --coordinates_format=eman2  --box_size=64  --import_ctf=outdir_cter/partres/partres.txt

Selected Micrographs Mode - Process all micrographs in a selection list file:
	In addition to path pattern of input micrographs and coordinates files, specify a name of micrograph selection list text file using --selection_list option.
	In this mode, only micrographs in the selection list which matches the file name part of the pattern (i.e. ignores the directory paths) will be processed.
	If a micrograph name in the selection list does not exists in the directory specified by the micrograph path pattern, processing of the micrograph will be skipped.

	mpirun  -np  32  sxwindow.py  './mic*.hdf'  'info/mic*_info.json'  particles  --selection_list=mic_list.txt  --coordinates_format=eman2  --box_size=64  --import_ctf=outdir_cter/partres/partres.txt

Single Micrograph Mode - Process a single micrograph:
	In addition to path pattern of input micrographs and coordinates files, specify a single micrograph name using --selection_list option.
	In this mode, only the specified single micrograph will be processed.
	If this micrograph name does not matches the file name part of the pattern (i.e. ignores the directory paths), the process will exit without processing it.
	If this micrograph name matches the file name part of the pattern but does not exists in the directory which specified by the micrograph path pattern, again the process will exit without processing it.

	sxwindow.py  './mic*.hdf'  'info/mic*_info.json'  particles  --selection_list=mic0.hdf  --coordinates_format=eman2  --box_size=64  --import_ctf=outdir_cter/partres/partres.txt

"""
	parser = OptionParser(usage, version=SPARXVERSION)
	parser.add_option("--selection_list",      type="string",        default=None,      help="Micrograph selecting list: Specify a name of micrograph selection list text file for Selected Micrographs Mode. The file extension must be \'.txt\'. Alternatively, the file name of a single micrograph can be specified for Single Micrograph Mode. (default none)")
	parser.add_option("--coordinates_format",  type="string",        default="eman1",   help="Coordinate file format: Allowed values are \'sparx\', \'eman1\', \'eman2\', or \'spider\'. The sparx, eman2, and spider formats use the particle center as coordinates. The eman1 format uses the lower left corner of the box as coordinates. (default eman1)")
	parser.add_option("--box_size",            type="int",           default=256,       help="Particle box size [Pixels]: The x and y dimensions of square area to be windowed. The box size after resampling is assumed when resample_ratio < 1.0. (default 256)")
	parser.add_option("--skip_invert",         action="store_true",  default=False,     help="Skip invert image contrast: Use this option for negative staining data. By default, the image contrast is inverted for cryo data. (default False)")
	parser.add_option("--import_ctf",          type="string",        default="",        help="CTER CTF parameter file: The file produced by sxcter and normally called partres.txt. (default none)")
	parser.add_option("--limit_ctf",           action="store_true",  default=False,     help="Use CTF limit filter: Frequencies whose oscillation can not be properly modeled at the current pixel size are discarded in the images with the appropriate low-pass filter. This option requires --import_ctf. (default False)")
	parser.add_option("--astigmatism_error",   type="float",         default=360.0,     help="Astigmatism error limit [Degrees]: Set astigmatism to zero for all micrographs where the angular error computed by sxcter is larger than the desired value. (default 360.0)")
	parser.add_option("--resample_ratio",      type="float",         default=1.0,       help="Ratio between new and original pixel size: Use values between 0.0 and 1.0. The new pixel size is automatically recalculated when resample_ratio < 1.0 is used. (default 1.0)")
	parser.add_option("--check_consistency",   action="store_true",  default=False,     help="Check consistency of inputs: Create the text file containing the list of inconsistent Micrograph ID entries (i.e. inconsist_mic_list_file.txt). (default False)")
	
	# ====================================================================================
	# Prepare processing
	# ====================================================================================
	# ------------------------------------------------------------------------------------
	# Set up MPI related variables
	# ------------------------------------------------------------------------------------
	# Detect if program is running under MPI
	RUNNING_UNDER_MPI = "OMPI_COMM_WORLD_SIZE" in os.environ
	
	main_mpi_proc = 0
	if RUNNING_UNDER_MPI:
		from mpi import mpi_init
		from mpi import MPI_COMM_WORLD, mpi_comm_rank, mpi_comm_size, mpi_barrier, mpi_reduce, MPI_INT, MPI_SUM
		
		mpi_init(0, [])
		my_mpi_proc_id = mpi_comm_rank(MPI_COMM_WORLD)
		n_mpi_procs = mpi_comm_size(MPI_COMM_WORLD)
	else:
		my_mpi_proc_id = 0
		n_mpi_procs = 1
	
	(options, args) = parser.parse_args(sys.argv[1:])
	
	# ------------------------------------------------------------------------------------
	# Check error conditions of arguments and prepare variables for them
	# ------------------------------------------------------------------------------------
	mic_pattern = None
	coords_pattern = None
	out_dir = None
	error_status = None
	# not a real while, an if with the opportunity to use break when errors need to be reported
	while True:
		if len(args) != 3:
			error_status = ("Please check usage for number of arguments.\n Usage: " + usage + "\n" + "Please run %s -h for help." % (program_name), getframeinfo(currentframe()))
			break
		else: 
			assert (len(args) == 3)
			mic_pattern = args[0]
			coords_pattern = args[1]
			out_dir = args[2]
		
		if mic_pattern[:len("bdb:")].lower() == "bdb":
			error_status = ("BDB file can not be selected as input micrographs. Please convert the format, and restart the program. Run %s -h for help." % (program_name), getframeinfo(currentframe()))
			break
		
		if mic_pattern.find("*") == -1:
			error_status = ("Input micrograph file name pattern must contain wild card (*). Please check input_micrograph_pattern argument. Run %s -h for help." % (program_name), getframeinfo(currentframe()))
			break
		
		if coords_pattern.find("*") == -1:
			error_status = ("Input coordinates file name pattern must contain wild card (*). Please check input_coordinates_pattern argument. Run %s -h for help." % (program_name), getframeinfo(currentframe()))
			break
		
		if my_mpi_proc_id == main_mpi_proc:
			if os.path.exists(out_dir):
				error_status = ("Output directory exists. Please change the name and restart the program.", getframeinfo(currentframe()))
				break
		break
	if_error_then_all_processes_exit_program(error_status)
	assert(coords_pattern != None)
	assert(mic_pattern != None)
	assert(out_dir != None)
	
	# ------------------------------------------------------------------------------------
	# Check error conditions of options
	# ------------------------------------------------------------------------------------
	check_options(options, program_name)
	
	# ------------------------------------------------------------------------------------
	# Prepare the variables for all sections
	# ------------------------------------------------------------------------------------
	# Micrograph basename pattern (directory path is removed from micrograph path pattern)
	mic_basename_pattern = os.path.basename(mic_pattern)
	
	# Global entry dictionary (all possible entries from all lists) for all mic id substring
	global_entry_dict = {} # mic id substring is the key
	subkey_input_mic_path = "Input Micrograph Path"
	subkey_selected_mic_basename = "Selected Micrograph Basename"
	subkey_coords_path = "Input Coordinates File Path"
	subkey_cter_entry = "CTER Entry"
	
	# List keeps only id substrings of micrographs whose all necessary information are available
	valid_mic_id_substr_list = [] 
	
	# Indices of CTER Parameters
	# All mpi processes must have access to indices
	if options.import_ctf:
		i_enum = -1
		i_enum += 1; idx_cter_def          = i_enum # defocus [um]; index must be same as ctf object format
		i_enum += 1; idx_cter_cs           = i_enum # Cs [mm]; index must be same as ctf object format
		i_enum += 1; idx_cter_vol          = i_enum # voltage[kV]; index must be same as ctf object format
		i_enum += 1; idx_cter_apix         = i_enum # pixel size [A]; index must be same as ctf object format
		i_enum += 1; idx_cter_bfactor      = i_enum # B-factor [A^2]; index must be same as ctf object format
		i_enum += 1; idx_cter_ac           = i_enum # amplitude contrast [%]; index must be same as ctf object format
		i_enum += 1; idx_cter_astig_amp    = i_enum # astigmatism amplitude [um]; index must be same as ctf object format
		i_enum += 1; idx_cter_astig_ang    = i_enum # astigmatism angle [degree]; index must be same as ctf object format
		i_enum += 1; idx_cter_sd_def       = i_enum # std dev of defocus [um]
		i_enum += 1; idx_cter_sd_astig_amp = i_enum # std dev of ast amp [A]
		i_enum += 1; idx_cter_sd_astig_ang = i_enum # std dev of ast angle [degree]
		i_enum += 1; idx_cter_cv_def       = i_enum # coefficient of variation of defocus [%]
		i_enum += 1; idx_cter_cv_astig_amp = i_enum # coefficient of variation of ast amp [%]
		i_enum += 1; idx_cter_spectra_diff = i_enum # average of differences between with- and without-astig. experimental 1D spectra at extrema
		i_enum += 1; idx_cter_error_def    = i_enum # frequency at which signal drops by 50% due to estimated error of defocus alone [1/A]
		i_enum += 1; idx_cter_error_astig  = i_enum # frequency at which signal drops by 50% due to estimated error of defocus and astigmatism [1/A]
		i_enum += 1; idx_cter_error_ctf    = i_enum # limit frequency by CTF error [1/A]
		i_enum += 1; idx_cter_mic_path     = i_enum # micrograph file path
		i_enum += 1; n_idx_cter            = i_enum
	
	# ====================================================================================
	# Obtain the list of micrograph id sustrings using a single CPU (i.e. main mpi process)
	# ====================================================================================
	# NOTE: Toshio Moriya 2016/10/24
	# The below is not a real while.  
	# It gives if-statements an opportunity to use break when errors need to be reported
	# However, more elegant way is to use 'raise' statement of exception mechanism...
	# 
	error_status = None
	while my_mpi_proc_id == main_mpi_proc:
		# --------------------------------------------------------------------------------
		# Prepare variables for this section
		# --------------------------------------------------------------------------------
		# Prefix and suffix of micrograph basename pattern 
		# to find the head/tail indices of micrograph id substring
		mic_basename_tokens = mic_basename_pattern.split('*')
		assert (len(mic_basename_tokens) == 2)
		# Find head index of micrograph id substring
		mic_id_substr_head_idx = len(mic_basename_tokens[0])
		
		# Prefix and suffix of coordinates file path pattern
		# to find the head/tail indices of coordinates file id substring
		coords_pattern_tokens = coords_pattern.split('*') 
		assert (len(coords_pattern_tokens) == 2)
		# Find head index of coordinates id substring
		coords_id_substr_head_idx = len(coords_pattern_tokens[0])
		
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
		if len(input_mic_path_list) == 0:
			error_status = ("No micrograph files are found in the directory specified by micrograph path pattern (%s). Please check input_micrograph_pattern argument. Run %s -h for help." % (os.path.dirname(mic_pattern), program_name), getframeinfo(currentframe()))
			break
		assert(len(input_mic_path_list) > 0)
		
		# Register micrograph id substrings to the global entry dictionary
		for input_mic_path in input_mic_path_list:
			# Find tail index of micrograph id substring and extract the substring from the micrograph name
			input_mic_basename = os.path.basename(input_mic_path)
			mic_id_substr_tail_idx = input_mic_basename.index(mic_basename_tokens[1])
			mic_id_substr = input_mic_basename[mic_id_substr_head_idx:mic_id_substr_tail_idx]
			assert(input_mic_path == mic_pattern.replace("*", mic_id_substr))
			if not mic_id_substr in global_entry_dict:
				# print("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from input_mic_path_list " % (mic_id_substr))
				global_entry_dict[mic_id_substr] = {}
			assert(mic_id_substr in global_entry_dict)
			global_entry_dict[mic_id_substr][subkey_input_mic_path] = input_mic_path
		assert(len(global_entry_dict) > 0)
		
		# --------------------------------------------------------------------------------
		# Register micrograph id substrings found in the selection list
		# to the global entry dictionary
		# --------------------------------------------------------------------------------
		# Generate the list of selected micrograph paths in the selection file
		selected_mic_path_list = []
		# Generate micrograph lists according to the execution mode
		if options.selection_list == None:
			print(" ")
			print("----- Running with All Micrographs Mode -----")
			# Treat all micrographs in the input directory as selected ones
			selected_mic_path_list = input_mic_path_list
		else:
			assert(options.selection_list != None and os.path.exists(options.selection_list))
			if os.path.splitext(options.selection_list)[1] == ".txt":
				print(" ")
				print("----- Running with Selected Micrographs Mode -----")
				print(" ")
				print("Checking the selection list...")
				selected_mic_path_list = read_text_file(options.selection_list)
				
				# Check error condition of micrograph entry lists
				print("Found %d microgarph entries in %s." % (len(selected_mic_path_list), options.selection_list))
				if len(selected_mic_path_list) == 0:
					error_status = ("No micrograph entries are found in the selection list file. Please check selection_list option. Run %s -h for help." % (program_name), getframeinfo(currentframe()))
					break
			else:
				print(" ")
				print("----- Running with Single Micrograph Mode -----")
				print(" ")
				print("Processing a single micorgprah: %s..." % (options.selection_list))
				selected_mic_path_list = [options.selection_list]
			assert(len(selected_mic_path_list) > 0)
			
			selected_mic_directory = os.path.dirname(selected_mic_path_list[0])
			if selected_mic_directory != "":
				print("    NOTE: Program disregards the directory paths in the selection list (%s)." % (selected_mic_directory))
			
		assert(len(selected_mic_path_list) > 0)
		
		# Register micrograph id substrings to the global entry dictionary
		for selected_mic_path in selected_mic_path_list:
			# Find tail index of micrograph id substring and extract the substring from the micrograph name
			selected_mic_basename = os.path.basename(selected_mic_path)
			mic_id_substr_tail_idx = selected_mic_basename.index(mic_basename_tokens[1])
			mic_id_substr = selected_mic_basename[mic_id_substr_head_idx:mic_id_substr_tail_idx]
			assert(selected_mic_basename == mic_basename_pattern.replace("*", mic_id_substr))
			if not mic_id_substr in global_entry_dict:
				# print("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from selected_mic_path_list " % (mic_id_substr))
				global_entry_dict[mic_id_substr] = {}
			assert(mic_id_substr in global_entry_dict)
			global_entry_dict[mic_id_substr][subkey_selected_mic_basename] = selected_mic_basename
		assert(len(global_entry_dict) > 0)
		
		del selected_mic_path_list # Do not need this anymore
		del input_mic_path_list # Do not need this anymore
		
		# --------------------------------------------------------------------------------
		# Register coordinates id substrings in coordinate path list to the global entry dictionary.
		# coordinates id substring (coords_id_substr) and micrograph id substring (mic_id_substr)
		# should be the same for the associated pair of micrograph and coordnates file.
		# --------------------------------------------------------------------------------
		print(" ")
		print("Checking the coordinates files...")
		coords_path_list = glob.glob(coords_pattern)
		
		# Check error condition of coordinates file path list
		print("Found %d coordinates files in %s directory." % (len(coords_path_list), os.path.dirname(coords_pattern)))
		if len(coords_path_list) == 0:
			error_status = ("No coordinates files are found in the directory specified by coordinates file path pattern (%s). Please check input_coordinates_pattern argument. Run %s -h for help." % (os.path.dirname(coords_pattern), program_name), getframeinfo(currentframe()))
			break
		assert(len(coords_path_list) > 0)
		
		for coords_path in coords_path_list:
			# Find tail index of coordinates id substring and extract the substring from the coordinates file path
			coords_id_substr_tail_idx = coords_path.index(coords_pattern_tokens[1])
			coords_id_substr = coords_path[coords_id_substr_head_idx:coords_id_substr_tail_idx]
			assert(coords_path == coords_pattern.replace("*", coords_id_substr))
			if not coords_id_substr in global_entry_dict:
				# print("MRK_DEBUG: Added new coords_id_substr (%s) to global_entry_dict from coords_path_list " % (coords_id_substr))
				global_entry_dict[coords_id_substr] = {}
			assert(coords_id_substr in global_entry_dict)
			global_entry_dict[coords_id_substr][subkey_coords_path] = coords_path
		assert(len(global_entry_dict) > 0)
		
		del coords_path_list # Do not need this anymore
		
		# --------------------------------------------------------------------------------
		# If necessary, register micrograph id substrings of CTER entries to the global entry dictionary
		# --------------------------------------------------------------------------------
		if options.import_ctf:
			print(" ")
			print("Checking the CTER CTF file...")
			assert(os.path.exists(options.import_ctf))
			cter_entry_list = read_text_row(options.import_ctf)
			
			# Check error condition of CTER entry list
			print("Found %d CTER entries in %s." % (len(cter_entry_list), options.import_ctf))
			if len(cter_entry_list) == 0:
				error_status = ("No CTER entries are found in %s. Please check --import_ctf option. Run %s -h for help." % (options.import_ctf, program_name), getframeinfo(currentframe()))
				break
			assert(len(cter_entry_list) > 0)
			
			if (len(cter_entry_list[0]) != n_idx_cter):
				error_status = ("Number of columns (%d) must be %d in %s. The format might be old. Please run sxcter.py again." % (len(cter_entry_list[0]), n_idx_cter, options.import_ctf), getframeinfo(currentframe()))
				break
			assert(len(cter_entry_list[0]) == n_idx_cter)
			
			cter_mic_directory = os.path.dirname(cter_entry_list[0][idx_cter_mic_path])
			if cter_mic_directory != "":
				print("    NOTE: Program disregards the directory paths in the CTER CTF file (%s)." % (cter_mic_directory))
			
			for cter_entry in cter_entry_list:
				assert(len(cter_entry) == n_idx_cter)
				# Find tail index of micrograph id substring and extract the substring from the micrograph path of CTER entry
				cter_mic_path = cter_entry[idx_cter_mic_path]
				cter_mic_basename = os.path.basename(cter_mic_path)
				mic_id_substr_tail_idx = cter_mic_basename.index(mic_basename_tokens[1])
				mic_id_substr = cter_mic_basename[mic_id_substr_head_idx:mic_id_substr_tail_idx]
				# Between cter_mic_path and mic_path, directory paths might be different but the basenames should be same!
				assert(cter_mic_basename == mic_basename_pattern.replace("*", mic_id_substr))
				
				if(cter_entry[idx_cter_sd_astig_ang] > options.astigmatism_error):
					print("    NOTE: Astigmatism angular SD of %s (%f degree) exceeds specified limit (%f degree). Resetting astigmatism parameters to zeros..." % (cter_mic_basename, cter_entry[idx_cter_sd_astig_ang], options.astigmatism_error))
					cter_entry[idx_cter_astig_amp] = 0.0
					cter_entry[idx_cter_astig_ang] = 0.0
				
				if not mic_id_substr in global_entry_dict:
					# print("MRK_DEBUG: Added new mic_id_substr (%s) to global_entry_dict from cter_entry_list " % (mic_id_substr))
					global_entry_dict[mic_id_substr] = {}
				assert(mic_id_substr in global_entry_dict)
				global_entry_dict[mic_id_substr][subkey_cter_entry] = cter_entry
			assert(len(global_entry_dict) > 0)
			
			del cter_entry_list # Do not need this anymore
		
		# --------------------------------------------------------------------------------
		# Clean up variables related to registration to the global entry dictionary
		# --------------------------------------------------------------------------------
		del coords_id_substr_head_idx
		del coords_pattern_tokens
		del mic_id_substr_head_idx
		del mic_basename_tokens
		
		# --------------------------------------------------------------------------------
		# Create the list containing only valid micrograph id substrings
		# --------------------------------------------------------------------------------
		# Prepare lists to keep track of invalid (rejected) micrographs 
		no_input_mic_id_substr_list = []
		no_coords_mic_id_substr_list = []
		no_cter_entry_mic_id_substr_list = []
		
		print(" ")
		print("Checking the input datasets consistency...")
		
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
				
				# Check if associated coordinate file exists
				if not subkey_coords_path in mic_id_entry:
					coords_path = coords_pattern.replace("*", mic_id_substr)
					warinnig_messages.append("    associated coordinates file %s." % (coords_path))
					no_coords_mic_id_substr_list.append(mic_id_substr)
				
				if options.import_ctf:
					# Check if associated CTER entry exists
					if not subkey_cter_entry in mic_id_entry:
						mic_basename = mic_basename_pattern.replace("*", mic_id_substr)
						warinnig_messages.append("    associated CTF parameter entry with %s in the CTER paramter file %s." % (mic_basename, options.import_ctf))
						no_cter_entry_mic_id_substr_list.append(mic_id_substr)
				
				if len(warinnig_messages) > 0:
					print("WARNING!!! Micrograph ID %s does not have:" % (mic_id_substr))
					for warinnig_message in warinnig_messages:
						print(warinnig_message)
					print("    Ignores this as an invalid entry.")
				else:
					# print("MRK_DEBUG: adding mic_id_substr := ", mic_id_substr)
					valid_mic_id_substr_list.append(mic_id_substr)
			# else:
			# 	assert(not subkey_selected_mic_basename in mic_id_entry)
			# 	# This entry is not in the selection list. Do nothing
			
		# Check the input dataset consistency and save the result to a text file, if necessary.
		if options.check_consistency:
			# Create output directory
			assert(not os.path.exists(out_dir))
			os.mkdir(out_dir)
			
			# Open the consistency check file
			inconsist_mic_list_path = os.path.join(out_dir,"inconsist_mic_id_file.txt")
			print(" ")
			print("Generating the input datasets consistency report in %s..." % (inconsist_mic_list_path))
			inconsist_mic_list_file = open(inconsist_mic_list_path, "w")
			inconsist_mic_list_file.write("# The information about inconsistent micrograph IDs\n")
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
				
				# Check if associated coordinate file exists
				if not subkey_coords_path in mic_id_entry:
					coords_path = coords_pattern.replace("*", mic_id_substr)
					consistency_messages.append("    associated coordinates file %s." % (coords_path))
				
				if options.import_ctf:
					# Check if associated CTER entry exists
					if not subkey_cter_entry in mic_id_entry:
						mic_basename = mic_basename_pattern.replace("*", mic_id_substr)
						consistency_messages.append("    associated CTF parameter entry with %s in the CTER paramter file %s." % (mic_basename, options.import_ctf))
				
				if len(consistency_messages) > 0:
					inconsist_mic_list_file.write("Micrograph ID %s does not have:\n" % (mic_id_substr))
					for consistency_message in consistency_messages:
						inconsist_mic_list_file.write(consistency_message)
						inconsist_mic_list_file.write("\n")
			
			# Close the consistency check file, if necessary
			inconsist_mic_list_file.flush()
			inconsist_mic_list_file.close()
			
		# --------------------------------------------------------------------------------
		# Print out the summary of input consistency
		# --------------------------------------------------------------------------------
		print(" ")
		print("Summary of dataset consistency check...")
		print("Detected                           : %6d" % (len(global_entry_dict)))
		print("Valid                              : %6d" % (len(valid_mic_id_substr_list)))
		print("Rejected by no coordinates file    : %6d" % (len(no_coords_mic_id_substr_list)))
		if options.import_ctf:
			print("Rejected by no CTER entry          : %6d" % (len(no_cter_entry_mic_id_substr_list)))
		
		# Since mic_id_substr is once stored as the key of global_entry_dict and extracted with the key order
		# we need sort the valid_mic_id_substr_list here
		# print("MRK_DEBUG: before sort, no_cter_entry_mic_id_substr_list := ", valid_mic_id_substr_list)
		valid_mic_id_substr_list.sort()
		# print("MRK_DEBUG: after sort, no_cter_entry_mic_id_substr_list := ", valid_mic_id_substr_list)
		
		# --------------------------------------------------------------------------------
		# Check MPI error condition
		# --------------------------------------------------------------------------------
		if len(valid_mic_id_substr_list) < n_mpi_procs:
			error_status = ("Number of MPI processes (%d) supplied by --np in mpirun cannot be greater than %d (number of valid micrographs that satisfy all criteria to be processed)." % (n_mpi_procs, len(valid_mic_id_substr_list)), getframeinfo(currentframe()))
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
	box_size = options.box_size
	box_half = box_size // 2
	mask2d = model_circle(box_size//2, box_size, box_size) # Create circular 2D mask to Util.infomask of particle images
	resample_ratio = options.resample_ratio
	
	# Prepare the function for reading coordinates files with the specified format.
	# This way, the following switch statement is unnecessary inside of the coordinates loop.
	coords_format = options.coordinates_format.lower()
	read_coords_file = None
	if coords_format == "sparx":
		read_coords_file = read_sparx_coords_file
	elif coords_format == "eman1":
		read_coords_file = read_eman1_coords_file
	elif coords_format == "eman2":
		read_coords_file = read_eman2_coords_file
	elif coords_format == "spider":
		read_coords_file = read_spider_coords_file
	else: 
		assert (False) # Unreachable code
	assert(read_coords_file != None)
	
	# Preapre variables related to CTF limit option
	if options.limit_ctf:
		assert (options.import_ctf)
		abs_ctf_limit_histogram = []  # compute the histogram for micrographs cut of by cter_entry limit.
			
	# Micrograph baseroot pattern (extension are removed from micrograph basename pattern)
	# for substack file names
	mic_baseroot_pattern = os.path.splitext(mic_basename_pattern)[0]  
	
	# Prepare the variables for the global summary of micrographs
	n_mic_process = 0
	n_mic_reject_no_coords_entry = 0
	n_global_coords_detect = 0
	n_global_coords_process = 0
	n_global_coords_reject_out_of_boundary = 0
	
	# keep a copy of the root output directory where the final bdb will be created
	root_out_dir = out_dir
	unsliced_valid_serial_id_list = valid_mic_id_substr_list
	if RUNNING_UNDER_MPI:
		mpi_barrier(MPI_COMM_WORLD)
		# All mpi processes should know global entry directory and valid micrograph id substring list
		global_entry_dict = wrap_mpi_bcast(global_entry_dict, main_mpi_proc)
		valid_mic_id_substr_list = wrap_mpi_bcast(valid_mic_id_substr_list, main_mpi_proc)
		
		# Slice the list of valid micrograph id substrings for this mpi process
		mic_start, mic_end = MPI_start_end(len(valid_mic_id_substr_list), n_mpi_procs, my_mpi_proc_id)
		valid_mic_id_substr_list = valid_mic_id_substr_list[mic_start:mic_end]
		
		# generate subdirectories of out_dir, one for each process
		out_dir = os.path.join(out_dir, "mpi_proc_%03d" % my_mpi_proc_id)
	
	# Set up progress message
	if my_mpi_proc_id == main_mpi_proc:
		print(" ")
		print("Micrographs processed by main process (including percent of progress):")
		progress_percent_step = len(valid_mic_id_substr_list)/100.0 # the number of micrograms for main node divided by 100
	
	# ------------------------------------------------------------------------------------
	# Starting main parallel execution
	# ------------------------------------------------------------------------------------
	for mic_id_substr_idx, mic_id_substr in enumerate(valid_mic_id_substr_list):
		
		# --------------------------------------------------------------------------------
		# Print out progress if necessary
		# --------------------------------------------------------------------------------
		mic_basename = global_entry_dict[mic_id_substr][subkey_selected_mic_basename]
		assert(mic_basename == mic_basename_pattern.replace("*", mic_id_substr))
		if my_mpi_proc_id == main_mpi_proc:
			print("%s ---> % 2.2f%%" % (mic_basename, mic_id_substr_idx / progress_percent_step))
		
		# --------------------------------------------------------------------------------
		# Read the associated coordinates according to the specified format and 
		# make the coordinates the center of particle image if necessary
		# Do this first because error might happen 
		# --------------------------------------------------------------------------------
		coords_path = global_entry_dict[mic_id_substr][subkey_coords_path]
		assert(os.path.exists(coords_path))
		assert(read_coords_file != None)
		coords_list = read_coords_file(coords_path)
		if (len(coords_list) == 0):
			print("For %s, the associate coordinates file %s does not contain any entries. Skipping..." % (mic_basename, coords_path))
			n_mic_reject_no_coords_entry += 1
			continue
		
		# --------------------------------------------------------------------------------
		# Get CTF parameter if necessary
		# Calculate the resampled pixel size and store it to the cter_entry if necessary
		# Do before expensive micrograph processing
		# --------------------------------------------------------------------------------
		if options.import_ctf:
			cter_entry = global_entry_dict[mic_id_substr][subkey_cter_entry] # <<== MUST BE FIXED!!!
			pixel_size = cter_entry[idx_cter_apix]
			
			if resample_ratio < 1.0:
				assert (resample_ratio > 0.0)
				resampled_pixel_size = pixel_size / resample_ratio
				if my_mpi_proc_id == main_mpi_proc:
					print("Resample micrograph to pixel size %6.4f and window segments from resampled micrograph." % resampled_pixel_size)
			else:
				assert (resample_ratio == 1.0)
				resampled_pixel_size = pixel_size
			
			# store the resampled pixel size to the cter_entry
			cter_entry[idx_cter_apix] = resampled_pixel_size
			
			# Generate CTF object of this micrograph
			from utilities import generate_ctf
			ctf_obj = generate_ctf(cter_entry) # indexes 0 to 7 (idx_cter_def to idx_cter_astig_ang) must be same in cter format & ctf object format.
		
		else:
			assert (not options.import_ctf)
			if resample_ratio < 1.0:
				assert (resample_ratio > 0.0)
				if my_mpi_proc_id == main_mpi_proc:
					print("Resample micrograph with ratio %6.4f and window segments from resampled micrograph." % resample_ratio)
			else:
				assert (resample_ratio == 1.0)
			# Note resampled_pixel_size is not declared if not options.import_ctf!
		
		# --------------------------------------------------------------------------------
		# Read micrograph
		# --------------------------------------------------------------------------------
		mic_path = global_entry_dict[mic_id_substr][subkey_input_mic_path]
		assert(mic_path == mic_pattern.replace("*", mic_id_substr))
		mic_img = get_im(mic_path)
		
		# --------------------------------------------------------------------------------
		# Move to the Fourier spsace processing
		# --------------------------------------------------------------------------------
		fftip(mic_img) # In-place fft
		
		# --------------------------------------------------------------------------------
		# If necessary, apply the hyperbolic tangent low-pass Fourier filter based on the (resampled) CTF limit;
		# Cut off frequency components higher than the (resampled) CTF limit.
		# --------------------------------------------------------------------------------
		if options.limit_ctf:
			assert (options.import_ctf)
			
			# Comput absolute frequency of CTF limit (abs_ctf_limit) with the original pixel size
			abs_ctf_limit, angstrom_ctf_limit = ctflimit(box_size, cter_entry[idx_cter_def], cter_entry[idx_cter_cs], cter_entry[idx_cter_vol], resampled_pixel_size)
			
			# Adjust the CTF limit according to the resampling ratio and box size
			if resample_ratio < 1.0:
				assert (resample_ratio > 0.0)
				abs_ctf_limit = resample_ratio * abs_ctf_limit / float(box_size)
			else:
				assert (resample_ratio == 1.0) # -> pixel_size == resampled_pixel_size -> pixel_size / resampled_pixel_size == 1.0
				abs_ctf_limit = abs_ctf_limit / float(box_size)
			
			# If ctf limit is lower than Nyquist frequency, apply the low pass filter with the cutoff at CTF limit frequencye.
			if abs_ctf_limit < 0.5:
				mic_img = filt_tanl(mic_img, abs_ctf_limit, 0.01)
				abs_ctf_limit_histogram.append(abs_ctf_limit)
		
		# --------------------------------------------------------------------------------
		# Apply the Gaussian high-pass Fourier filter to micrograph based on the (resampled) box size;
		# Cut off frequency components lower than one that the (resampled) box size can express.
		# Then, move back to the real space processing
		# --------------------------------------------------------------------------------
		mic_img = fft(filt_gaussh(mic_img, resample_ratio / box_size))
		
		# --------------------------------------------------------------------------------
		# Resample micrograph, map coordinates, and window segments from resampled micrograph using new coordinates
		# after resampling by resample_ratio, new pixel size will be pixel_size/resample_ratio = resampled_pixel_size
		# --------------------------------------------------------------------------------
		# NOTE: 2015/04/13 Toshio Moriya
		# resample() efficiently takes care of the case resample_ratio = 1.0 but
		# it does not set apix_*. Even though it sets apix_* when resample_ratio < 1.0...
		mic_img = resample(mic_img, resample_ratio)
		
		# --------------------------------------------------------------------------------
		# If necessary, invert image contrast of this micrograph 
		# --------------------------------------------------------------------------------
		if not options.skip_invert:
			mic_stats = Util.infomask(mic_img, None, True) # mic_stat[0:mean, 1:SD, 2:min, 3:max]
			Util.mul_scalar(mic_img, -1.0)
			mic_img += 2 * mic_stats[0]
		
		# --------------------------------------------------------------------------------
		# Generate the output file path of particle stack for this mpi process
		# --------------------------------------------------------------------------------
		mic_baseroot = mic_baseroot_pattern.replace("*", mic_id_substr)
		local_stack_path  = "bdb:%s#" % out_dir + mic_baseroot + "_ptcls"
		
		# --------------------------------------------------------------------------------
		# Prepare coordinates loop variables
		# --------------------------------------------------------------------------------
		nx = mic_img.get_xsize() 
		ny = mic_img.get_ysize()
		x0 = nx//2
		y0 = ny//2
		
		local_particle_id = 0 # can be different from coordinates_id
		coords_reject_out_of_boundary_messages = []
		
		# Loop through coordinates
		for coords_id in xrange(len(coords_list)):
			
			x = int(coords_list[coords_id][0])
			y = int(coords_list[coords_id][1])
			
			if resample_ratio < 1.0:
				assert (resample_ratio > 0.0)
				x = int(x * resample_ratio)
				y = int(y * resample_ratio)
			else:
				assert(resample_ratio == 1.0)
			
			if( (0 <= x - box_half) and ( x + box_half <= nx ) and (0 <= y - box_half) and ( y + box_half <= ny ) ):
				particle_img = Util.window(mic_img, box_size, box_size, 1, x-x0, y-y0)
			else:
				coords_reject_out_of_boundary_messages.append("coordinates ID = %04d: x = %4d, y = %4d, box_size = %4d " % (coords_id, x, y, box_size))
				# print("MRK_DEBUG: coords_reject_out_of_boundary_messages[-1] := %s" % coords_reject_out_of_boundary_messages[-1])
				continue
			
			particle_img = ramp(particle_img)
			particle_stats = Util.infomask(particle_img, mask2d, False) # particle_stats[0:mean, 1:SD, 2:min, 3:max]
			particle_img -= particle_stats[0]
			particle_img /= particle_stats[1]
			
			# NOTE: 2015/04/09 Toshio Moriya
			# ptcl_source_image might be redundant information...
			# Consider re-organizing header entries...
			# 
			particle_img.set_attr("ptcl_source_image", mic_path)
			particle_img.set_attr("ptcl_source_coord_id", coords_id)
			particle_img.set_attr("ptcl_source_coord", [int(coords_list[coords_id][0]), int(coords_list[coords_id][1])])
			particle_img.set_attr("resample_ratio", resample_ratio)
			
			# NOTE: 2015/04/13 Toshio Moriya
			# apix_* attributes are updated by resample() only when resample_ratio != 1.0
			# Let's make sure header info is consistent by setting apix_* = 1.0 
			# regardless of options, so it is not passed down the processing line
			# 
			particle_img.set_attr("apix_x", 1.0)
			particle_img.set_attr("apix_y", 1.0)
			particle_img.set_attr("apix_z", 1.0)
			if options.import_ctf:
				particle_img.set_attr("ctf",ctf_obj)
				particle_img.set_attr("ctf_applied", 0)
				particle_img.set_attr("pixel_size", pixel_size)
				# particle_img.set_attr("apix_x", resampled_pixel_size)
				# particle_img.set_attr("apix_y", resampled_pixel_size)
				# particle_img.set_attr("apix_z", resampled_pixel_size)
			# NOTE: 2015/04/13 Toshio Moriya 
			# Pawel Comment: Micrograph is not supposed to have CTF header info.
			# So, let's assume it does not exist & ignore its presence.
			# Note that resample() "correctly" updates pixel size of CTF header info if it exists
			# 
			# elif (particle_img.has_ctff()):
			# 	assert(not options.import_ctf)
			# 	ctf_origin = particle_img.get_attr("ctf_obj")
			# 	pixel_size = round(ctf_origin.apix, 5) # Because SXCTER ouputs up to 5 digits 
			# 	particle_img.set_attr("apix_x",pixel_size)
			# 	particle_img.set_attr("apix_y",pixel_size)
			# 	particle_img.set_attr("apix_z",pixel_size)	
			
			# print("MRK_DEBUG: local_stack_path, local_particle_id", local_stack_path, local_particle_id)
			particle_img.write_image(local_stack_path, local_particle_id)
			local_particle_id += 1
		
		# Save the message list of rejected coordinates because of out-of-boundary
		# print("MRK_DEBUG: len(coords_reject_out_of_boundary_messages) := %d" % len(coords_reject_out_of_boundary_messages))
		if len(coords_reject_out_of_boundary_messages) > 0:
			# Open file path to save the message list
			coords_reject_out_of_boundary_path = os.path.join(root_out_dir, os.path.splitext(os.path.basename(coords_path))[0] + "_reject_out_of_boundary.txt")
			# print("MRK_DEBUG: coords_reject_out_of_boundary_path := %s" % coords_reject_out_of_boundary_path)
			coords_reject_out_of_boundary_file = open(coords_reject_out_of_boundary_path, "w")
			
			for coords_reject_out_of_boundary_message in coords_reject_out_of_boundary_messages:
				coords_reject_out_of_boundary_file.write(coords_reject_out_of_boundary_message)
				coords_reject_out_of_boundary_file.write("\n")
			
			# Close the consistency check file, if necessary
			coords_reject_out_of_boundary_file.flush()
			coords_reject_out_of_boundary_file.close()
		
		n_mic_process += 1
		n_global_coords_detect += len(coords_list)
		n_global_coords_process += local_particle_id
		n_global_coords_reject_out_of_boundary += len(coords_reject_out_of_boundary_messages)
		
		# Release the data base of local stack from this process
		# so that the subprocess can access to the data base
		db_close_dict(local_stack_path)
	
	# ------------------------------------------------------------------------------------
	# Print out CTF limit information
	# ------------------------------------------------------------------------------------
	if options.limit_ctf:
		assert (options.import_ctf)
		if RUNNING_UNDER_MPI:
			abs_ctf_limit_histogram = wrap_mpi_gatherv(abs_ctf_limit_histogram, main_mpi_proc)
		
		if my_mpi_proc_id == main_mpi_proc:
			# Print out the summary of CTF limit absolute frequency
			print(" ")
			print("Global summary of CTF limit absolute frequency (--limit_ctf)...")
			print("Percentage of filtered micrographs: %8.2f\n" % (len(abs_ctf_limit_histogram) * 100.0 / len(unsliced_valid_serial_id_list)))
			
			n_bins = 10
			if len(abs_ctf_limit_histogram) >= n_bins:
				from statistics import hist_list
				cutoff_region, cutoff_counts = hist_list(abs_ctf_limit_histogram, n_bins)
				print("Histogram of CTF limit absolute frequency used for the filtering:")
				print("      CTF limit       counts")
				for bin_id in xrange(n_bins):
					print(" %14.7f     %7d" % (cutoff_region[bin_id], cutoff_counts[bin_id]))
			else:
				print("The number of filtered micrographs (%d) is less than the number of bins (%d). No histogram is produced." % (len(abs_ctf_limit_histogram), n_bins))
	
	# ------------------------------------------------------------------------------------
	# Print out summary of processing
	# ------------------------------------------------------------------------------------
	if RUNNING_UNDER_MPI:
		n_mic_process = mpi_reduce(n_mic_process, 1, MPI_INT, MPI_SUM, main_mpi_proc, MPI_COMM_WORLD)
		n_mic_reject_no_coords_entry = mpi_reduce(n_mic_reject_no_coords_entry, 1, MPI_INT, MPI_SUM, main_mpi_proc, MPI_COMM_WORLD)
		n_global_coords_detect = mpi_reduce(n_global_coords_detect, 1, MPI_INT, MPI_SUM, main_mpi_proc, MPI_COMM_WORLD)
		n_global_coords_process = mpi_reduce(n_global_coords_process, 1, MPI_INT, MPI_SUM, main_mpi_proc, MPI_COMM_WORLD)
		n_global_coords_reject_out_of_boundary = mpi_reduce(n_global_coords_reject_out_of_boundary, 1, MPI_INT, MPI_SUM, main_mpi_proc, MPI_COMM_WORLD)
	
	# Print out the summary of all micrographs
	if main_mpi_proc == my_mpi_proc_id:
		print(" ")
		print("Summary of micrograph level processing...")
		print("Valid                              : %6d" % (len(unsliced_valid_serial_id_list)))
		print("Processed                          : %6d" % (n_mic_process))
		print("Rejected by no coordinates entries : %6d" % (n_mic_reject_no_coords_entry))
		print(" ")
		print("Global summary of coordinates level processing...")
		print("Detected                           : %6d" % (n_global_coords_detect))
		print("Processed                          : %6d" % (n_global_coords_process))
		print("Rejected by out of boundary        : %6d" % (n_global_coords_reject_out_of_boundary))
		if n_global_coords_reject_out_of_boundary > 0:
			print("    NOTE: Information of rejected coordinates by out of boundary are saved in %s files." % (os.path.join(root_out_dir, os.path.splitext(os.path.basename(coords_pattern))[0] + "_reject_out_of_boundary.txt" )))
	
	if RUNNING_UNDER_MPI:
		mpi_barrier(MPI_COMM_WORLD)
	
	if main_mpi_proc == my_mpi_proc_id:
		# NOTE: Toshio Moriya 2016/10/27
		# Pawel commented out the command execution below because of an MPI issue.
		# When the program is running, all CPUs are writing to the disk. 
		# However, at the end of it, there is no guarantee all files will be closed 
		# as seen from the main node. 
		# The only way to prevent the crash is to execute e2bdb after the program finished, 
		# and even with that one is well advised to wait. 
		"""
		print("\n Creating bdb:%s/data\n"%root_out_dir)
		for proc_i in range(n_mpi_procs):
			mic_start, mic_end = MPI_start_end(len(unsliced_valid_serial_id_list), n_mpi_procs, proc_i)
			for mic_id_substr in unsliced_valid_serial_id_list[mic_start:mic_end]:
				e2bdb_command = "e2bdb.py "
				mic_baseroot = mic_baseroot_pattern.replace("*", mic_id_substr)
				if RUNNING_UNDER_MPI:
					e2bdb_command += "bdb:" + os.path.join(root_out_dir,"%03d/"%proc_i) + mic_baseroot + "_ptcls "
				else:
					e2bdb_command += "bdb:" + os.path.join(root_out_dir, mic_baseroot + "_ptcls ") 
				
				e2bdb_command += " --appendvstack=bdb:%s/data  1>/dev/null"%root_out_dir
		"""
		if RUNNING_UNDER_MPI:
			e2bdb_command = "e2bdb.py  " + root_out_dir + "/mpi_proc_*  --makevstack=bdb:" + root_out_dir + "/data"
			print(" ")
			print("Please execute from the command line :  ", e2bdb_command)
		else:
			e2bdb_command = "e2bdb.py  " + root_out_dir + "  --makevstack=bdb:" + root_out_dir + "/data"
			cmdexecute(e2bdb_command, printing_on_success = False)
		
		print(" ")
		print("DONE!!!\n")
	
	if RUNNING_UNDER_MPI:
		mpi_barrier(MPI_COMM_WORLD)
		from mpi import mpi_finalize
		mpi_finalize()
	
	sys.stdout.flush()
	sys.exit(0)
	
# ========================================================================================
# Define main function for command line execution
# ========================================================================================
if __name__=="__main__":
	main()

# ========================================================================================
#  END OF FILE
# ========================================================================================
