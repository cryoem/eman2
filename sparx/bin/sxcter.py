#!/usr/bin/env python
#
# Author: Pawel A.Penczek and Edward H. Egelman 05/27/2009 (Pawel.A.Penczek@uth.tmc.edu)
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
# Copyright (c) 2008-Forever The University of Virginia
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
from optparse import OptionParser

import global_def
from global_def import *
from inspect import currentframe, getframeinfo
from utilities import if_error_then_all_processes_exit_program

def main():
	program_name = os.path.basename(sys.argv[0])
	usage = program_name + """  input_image_path  output_directory  --selection_list=selection_list  --wn=CTF_WINDOW_SIZE --apix=PIXEL_SIZE  --Cs=CS  --voltage=VOLTAGE  --ac=AMP_CONTRAST  --f_start=FREA_START  --f_stop=FREQ_STOP  --kboot=KBOOT  --overlap_x=OVERLAP_X  --overlap_y=OVERLAP_Y  --edge_x=EDGE_X  --edge_y=EDGE_Y  --set_ctf_header  --check_consistency  --stack_mode  --debug_mode

Automated estimation of CTF parameters with error assessment.

All Micrographs Mode - Process all micrographs in a directory: 
	Specify a list of input micrographs using a wild card (*), called here input micrographs path pattern. 
	Use the wild card to indicate the place of variable part of the file names (e.g. serial number, time stamp, and etc). 
	Running from the command line requires enclosing the string by single quotes (') or double quotes ("). 
	sxgui.py will automatically adds single quotes to the string. 
	BDB files can not be selected as input micrographs. 
	Then, specify output directory where all outputs should be saved. 
	In this mode, all micrographs matching the path pattern will be processed.

	mpirun -np 16 sxcter.py './mic*.hdf' outdir_cter --wn=512 --apix=2.29 --Cs=2.0 --voltage=300 --ac=10.0

Selected Micrographs Mode - Process all micrographs in a selection list file:
	In addition to input micrographs path pattern and output directry arguments, 
	specify a name of micrograph selection list text file using --selection_list option 
	(e.g. output of sxgui_unblur.py or sxgui_cter.py). The file extension must be ".txt". 
	In this mode, only micrographs in the selection list which matches the file name part of the pattern (ignoring the directory paths) will be processed. 
	If a micrograph name in the selection list does not exists in the directory specified by the micrograph path pattern, processing of the micrograph will be skipped.

	mpirun -np 16 sxcter.py './mic*.hdf' outdir_cter --selection_list=mic_list.txt --wn=512 --apix=2.29 --Cs=2.0 --voltage=300 --ac=10.0

Single Micrograph Mode - Process a single micrograph: 
	In addition to input micrographs path pattern and output directry arguments, 
	specify a single micrograph name using --selection_list option. 
	In this mode, only the specified single micrograph will be processed. 
	If this micrograph name does not matches the file name part of the pattern (ignoring the directory paths), the process will exit without processing it. 
	If this micrograph name matches the file name part of the pattern but does not exists in the directory which specified by the micrograph path pattern, again the process will exit without processing it. 
	Use single processor for this mode.

	sxcter.py './mic*.hdf' outdir_cter --selection_list=mic0.hdf --wn=512 --apix=2.29 --Cs=2.0 --voltage=300 --ac=10.0

Stack Mode - Process a particle stack (Not supported by SPHIRE GUI)):: 
	Use --stack_mode option, then specify the path of particle stack file (without wild card "*") and output directory as arguments. 
	This mode ignores --selection_list, --wn --overlap_x, --overlap_y, --edge_x, and --edge_y options. 
	Use single processor for this mode. Not supported by SPHIRE GUI (sxgui.py). 

	sxcter.py bdb:stack outdir_cter --apix=2.29 --Cs=2.0 --voltage=300 --ac=10.0 --stack_mode 

"""
	parser = OptionParser(usage, version=SPARXVERSION)
	parser.add_option("--selection_list",    type="string",        default=None,   help="Micrograph selecting list: Specify path of a micrograph selection list text file for Selected Micrographs Mode. The file extension must be \'.txt\'. Alternatively, the file name of a single micrograph can be specified for Single Micrograph Mode. (default none)")
	parser.add_option("--wn",                type="int",           default=512,    help="CTF window size [pixels]: The size should be slightly larger than particle box size. This will be ignored in Stack Mode. (default 512)")
	parser.add_option("--apix",              type="float",         default=-1.0,   help="Pixel size [A/Pixels]: The pixel size of input micrograph(s) or images in input particle stack. (default -1.0)")
	parser.add_option("--Cs",                type="float",         default=2.0,    help="Microscope spherical aberration (Cs) [mm]: The spherical aberration (Cs) of microscope used for imaging. (default 2.0)")
	parser.add_option("--voltage",           type="float",         default=300.0,  help="Microscope voltage [kV]: The acceleration voltage of microscope used for imaging. (default 300.0)")
	parser.add_option("--ac",                type="float",         default=10.0,   help="Amplitude contrast [%]: The typical amplitude contrast is in the range of 7% - 14%. The value mainly depends on the thickness of the ice embedding the particles. (default 10.0)")
	parser.add_option("--f_start",           type="float",         default=-1.0,   help="Lowest frequency [1/A]: Lowest frequency to be considered in the CTF estimation. Determined automatically by default. (default -1.0)")
	parser.add_option("--f_stop",            type="float",         default=-1.0,   help="Highest frequency [1/A]: Highest frequency to be considered in the CTF estimation. Determined automatically by default. (default -1.0)")
	parser.add_option("--kboot",             type="int",           default=16,     help="Number of CTF estimates per micrograph: Used for error assessment. (default 16)")
	parser.add_option("--overlap_x",         type="int",           default=50,     help="X overlap [%]: Overlap between the windows in the x direction. This will be ignored in Stack Mode. (default 50)")
	parser.add_option("--overlap_y",         type="int",           default=50,     help="Y overlap [%]: Overlap between the windows in the y direction. This will be ignored in Stack Mode. (default 50)")
	parser.add_option("--edge_x",            type="int",           default=0,      help="Edge x [pixels]: Defines the edge of the tiling area in the x direction. Normally it does not need to be modified. This will be ignored in Stack Mode. (default 0)")
	parser.add_option("--edge_y",            type="int",           default=0,      help="Edge y [pixels]: Defines the edge of the tiling area in the y direction. Normally it does not need to be modified. This will be ignored in Stack Mode. (default 0)")
	parser.add_option("--set_ctf_header",    action="store_true",  default=False,  help="Export CTF parameters to header: Exports the estimated CTF parameters to the image header. (default False)")
	parser.add_option("--check_consistency", action="store_true",  default=False,  help="Check consistency of inputs: Create a text file containing the list of inconsistent Micrograph ID entries (i.e. inconsist_mic_list_file.txt). (default False)")
	parser.add_option("--stack_mode",        action="store_true",  default=False,  help="Use stack mode: Use a stack as the input. Please set the file path of a stack as the first argument and output directory for the second argument. This is advanced option. Not supported by sxgui. (default False)")
	parser.add_option("--debug_mode",        action="store_true",  default=False,  help="Enable debug mode: Print out debug information. (default False)")
	
	(options, args) = parser.parse_args(sys.argv[1:])
	
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
		from mpi import mpi_init, mpi_comm_rank, mpi_comm_size, mpi_barrier, MPI_COMM_WORLD
		
		sys.argv = mpi_init(len(sys.argv), sys.argv)
		my_mpi_proc_id = mpi_comm_rank(MPI_COMM_WORLD)
		n_mpi_procs = mpi_comm_size(MPI_COMM_WORLD)
	else:
		my_mpi_proc_id = 0
		n_mpi_procs = 1
	
	# ------------------------------------------------------------------------------------
	# Set up SPHIRE global definitions
	# ------------------------------------------------------------------------------------
	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()
	
	# Change the name log file for error message
	original_logfilename = global_def.LOGFILE
	global_def.LOGFILE = os.path.splitext(program_name)[0] + '_' + original_logfilename + '.txt'
	
	# ------------------------------------------------------------------------------------
	# Check error conditions of arguments and options, then prepare variables for arguments
	# ------------------------------------------------------------------------------------
	input_image_path = None
	output_directory = None
	# not a real while, an if with the opportunity to use break when errors need to be reported
	error_status = None
	while True:
		# --------------------------------------------------------------------------------
		# Check the number of arguments. If OK, then prepare variables for them
		# --------------------------------------------------------------------------------
		if len(args) != 2:
			error_status = ("Please check usage for number of arguments.\n Usage: " + usage + "\n" + "Please run %s -h for help." % (program_name), getframeinfo(currentframe()))
			break
		assert (len(args) == 2)
		
		# NOTE: 2015/11/27 Toshio Moriya
		# Require single quotes (') or double quotes (") when input micrograph pattern is give for input_image_path
		#  so that sys.argv does not automatically expand wild card and create a list of file names
		#
		input_image_path = args[0]
		output_directory = args[1]
		
		# --------------------------------------------------------------------------------
		# NOTE: 2016/03/17 Toshio Moriya
		# cter_mrk() will take care of all the error conditions 
		# --------------------------------------------------------------------------------
		
		break
	if_error_then_all_processes_exit_program(error_status)
	assert (input_image_path != None)
	assert (output_directory != None)
	
	if my_mpi_proc_id == main_mpi_proc:
		command_line = ""
		for command_token in sys.argv:
			command_line += command_token + "  "
		print(" ")
		print("Shell line command:")
		print(command_line)
	
	from morphology import cter_mrk
	result = cter_mrk(input_image_path, output_directory, options.selection_list, options.wn, options.apix, options.Cs, options.voltage, options.ac, options.f_start, options.f_stop, options.kboot, options.overlap_x, options.overlap_y, options.edge_x, options.edge_y, options.set_ctf_header, options.check_consistency, options.stack_mode, options.debug_mode, program_name, RUNNING_UNDER_MPI, main_mpi_proc, my_mpi_proc_id, n_mpi_procs)
	
	if RUNNING_UNDER_MPI:
		mpi_barrier(MPI_COMM_WORLD)
	
	if main_mpi_proc == my_mpi_proc_id:
		if options.debug_mode:
			print "Returned value from cter_mrk() := ", result
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
	
	# ------------------------------------------------------------------------------------
	# Clean up MPI related variables
	# ------------------------------------------------------------------------------------
	if RUNNING_UNDER_MPI:
		mpi_barrier(MPI_COMM_WORLD)
		from mpi import mpi_finalize
		mpi_finalize()
	
	sys.stdout.flush()
	sys.exit(0)
	
if __name__ == "__main__":
	main()
