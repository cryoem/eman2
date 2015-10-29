#!/usr/bin/env python

import global_def
from   global_def import *

import os
import sys
from optparse import OptionParser

def main(sys_arg_list):
	progname = os.path.basename(sys_arg_list[0])
	usage = progname + """ stack outdir1 outdir2 --indir --nameroot --micsuffix --wn --apix --Cs --voltage --ac --kboot --MPI

	Process micrographs:
	
		Specify directory and prefix and suffix (type) of micrographs to process through --indir, --nameroot, and --micsuffix
		Specify output directories pwrot and partres as arguments.
		
			mpirun -np 16 sxcter.py pwrot partres --indir=. --nameroot=micrograph_PSC23_A8A_1GD_11112_135 --micsuffix=mrc --wn=512 --apix=2.29 --Cs=2.0 --voltage=300 --ac=10.0 --MPI
		After the program stops, it is advisable to concatenate all output files in partres directory:
		cd partres
		cat */* >>allctfs.txt

	Process stack:
	
		Specify name of stack and output directories as arguments.
			sxcter.py bdb:stack pwrot partres --apix=2.29 --Cs=2.0 --voltage=300 --ac=10.0

	
	"""
	parser = OptionParser(usage,version=SPARXVERSION)
	
	parser.add_option("--indir",            	  type="string",		    default=".",     				 help="<Micrograph Directory> Directory containing micrographs to be processed. (Default: current directory)")
	parser.add_option("--nameroot",          	  type="string",		    default="",     				 help="<Micrograph Root Name> Root name (Prefix) of micrographs to be processed.")
	parser.add_option("--micsuffix",              type=str,                 default="",                      help="<Micrograph Extension> A string denoting micrograph type. For example 'mrc', 'hdf', 'ser' ... ")
	parser.add_option("--wn",				  	  type="int",				default=512, 					 help="<Window Size> Size of window to use (should be slightly larger than particle box size). (Default 512)")
	
	parser.add_option("--apix",               	  type="float",			 	default=-1,               	     help="<Pixel Size> pixel size in Angstroms. (Default: Invalid Value)")   
	parser.add_option("--Cs",               	  type="float",			 	default=2.0,               	     help="<Microscope Cs> Microscope Cs (spherical aberation) in mm. (Default: 2.0)")
	parser.add_option("--voltage",				  type="float",				default=300.0, 					 help="<Microscope Voltage> Microscope voltage in kV. (Default: 300.0)")
	parser.add_option("--ac",					  type="float",				default=10.0, 					 help="<Amplitude contrast> Amplitude contrast in percentage. (Default: 10)")
	parser.add_option("--kboot",				  type="int",				default=16, 					 help="<kboot> Number of defocus estimates for micrograph to assess error. (Default: 16) (advanced)")
	parser.add_option("--MPI",               	  action="store_true",   	default=False,              	 help="use MPI version. (Default: not use MPI)")
	parser.add_option("--debug",               	  action="store_true",   	default=False,              	 help="<Enable Debug> debug. (Default: disable) (advanced)")
	parser.add_option("--overlap_x",			  type="int",				default=50, 					 help="<Overlap X> overlap x (Default: 50%) (advanced)")
	parser.add_option("--overlap_y",			  type="int",				default=50, 					 help="<Overlap Y> overlap y (Default: 50%) (advanced)")
	parser.add_option("--edge_x",			  	  type="int",				default=0, 					     help="<Edge X> edge x (Default: 0) (advanced)")
	parser.add_option("--edge_y",			      type="int",				default=0, 					     help="<Edge Y> edge y (Default: 0) (advanced)")
	parser.add_option("--f_start",                type="float",			 	default=-1.0,               	 help="<Start Frequency> starting frequency, units [1/A]. (Default: determined automatically) (advanced)")   
	parser.add_option("--f_stop",                 type="float",			 	default=-1.0,               	 help="<Stop Frequency> stop frequency, units [1/A]. (Default: determined automatically) (advanced)")

	# Only for GUI support
	from optparse import OptionGroup
	# Add argument group to parser
	arg_group = OptionGroup(parser, "Arguments", "These options are interpretated as arguments by GUI.")

	arg_group.add_option("--sxgui_arguments",                                                                   help="GUI uses this option to get argument group.")
	# arg_group.add_option("--stack",                  type="string",		       default=".",                     help="<Particle Stack> Set of 2-D images in a stack file. ")
	arg_group.add_option("--outdir1",                type="string",		       default="",                      help="<Output Directory 1> Output Directory to store a file containing the list of estimated CTF parameters for all image (partres).")
	arg_group.add_option("--outdir2",                type="string",            default="",                      help="<Output Directory 2> Output Directory to store files containg experimental and fitted rotational averaged power spectrum for images (rotinf****).")

	parser.add_option_group(arg_group)

	
	# NOTE: 2015/10/22 Toshio Moriya
	# The followings are necessary because
	# some scripts support only MPI version but does not have --MPI option, and
	# some other scripts does not support MPI and does not have --MPI option.
			
	# Add MPI related option group to parser
	mpi_group = OptionGroup(parser, "MPI Options", "These options are used only by GUI.")
			
	mpi_group.add_option("--sxgui_mpi_options",                                       help="GUI uses this option to get MPI option group.")
	mpi_group.add_option("--MPI_support",        action="store_true",  default=True,  help="No --MPI option doesn't always mean that script does not support MPI.")
	mpi_group.add_option("--MPI_add_flag",       action="store_true",  default=True,  help="Need to add --MPI in command line.")
	
	parser.add_option_group(mpi_group)
	
	return parser
