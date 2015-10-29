#!/usr/bin/env python

import os, sys
import json

# from optparse import *
from optparse import OptionParser
# from argparse import *
from EMAN2 import *
from EMAN2db import *
from EMAN2jsondb import *
# from emboxerbase import *
from sparx import *

def main(sys_arg_list):
	progname = os.path.basename(sys_arg_list[0])
	usage = progname + " micrographs_list  --coords_dir=coords_dir  --coords_suffix=coords_suffix" + \
	                                          "  --coords_extension=coords_extension  --coords_format=coords_format" + \
	                                          "  --indir=input_dir  --importctf=ctf_file  --limitctf" + \
	                                          "  --resample_ratio=resample_ratio  --box_size=box_size" + \
	                                          "  --outdir=outdir  --outsuffix=outsuffix  --micsuffix=micsuffix" + \
	                                          "  --nameroot=nameroot  --invert" + \
	                                          "  --defocuserror=defocuserror  --astigmatismerror=astigmatismerror"
	parser = OptionParser(usage, version=SPARXVERSION)

	parser.add_option("--coords_dir",       type="string",        default=".",       help="<Coordinates Directory> Directory containing files with particle coordinates. (Default: current directory)")
	parser.add_option("--coords_suffix",    type="string",        default="",        help="<Coordinates File Suffix> Suffix of coordinate files. For example '_ptcls'. ")
	parser.add_option("--coords_extension", type="string",        default="box",     help="<Coordinates File Extension> File extension of coordinate files. e.g 'box' for eman1, 'json' for eman2, ...")
	parser.add_option("--coords_format",    type="string",        default="eman1",   help="<Coordinates File Format> Format of coordinates file: 'sparx', 'eman1', 'eman2', or 'spider'. The coordinates of sparx, eman2, and spider format is particle center. The coordinates of eman1 format is particle box conner associated with the original box size.")
	parser.add_option("--indir",            type="string",        default=".",       help="<Micrograph Directory> Directory containing micrographs to be processed. (Default: current directory)")
	parser.add_option("--nameroot",         type="string",        default="",        help="<Micrograph Root Name> Root name (Prefix) of micrographs to be processed.")
	parser.add_option("--micsuffix",        type="string",        default="hdf",     help="<Micrograph Extension> A string denoting micrograph type. (Default 'hdf')")
	parser.add_option("--outdir",           type="string",        default=".",       help="<Output Directory> Output directory (Default: current directory)")
	parser.add_option("--outsuffix",        type="string",        default="_ptcls",  help="<Output File Suffix> Suffix for output stack. (Default '_ptcls')")
	parser.add_option("--importctf",        type="string",        default="",        help="<CTER CTF File> File name with CTF parameters produced by sxcter.") 
	parser.add_option("--box_size",         type="int",           default=256,       help="<Box Size> x and y dimension in pixels of square area to be windowed. Pixel size after resampling is assumed when resample_ratio < 1.0 (Default 256)")
	parser.add_option("--invert",           action="store_true",  default=False,     help="<Invert Contrast> Invert image contrast (recommended for cryo data) (Default, no contrast inversion)")
	parser.add_option("--resample_ratio",   type="float",         default=1.0,       help="<Resample Ratio> Ratio of new to old image size (or old to new pixel size) for resampling. Valid range is 0.0 < resample_ratio <= 1.0. (Default: 1.0)  (advanced)")
	parser.add_option("--limitctf",         action="store_true",  default=False,     help="<Apply CTF-Limit Filter> Filter micrographs based on the CTF limit. It requires --importctf. (Default: no filter) (advanced)")	
	parser.add_option("--defocuserror",     type="float",         default=1000000.0, help="<Defocus Error Limit> Exclude micrographs whose relative defocus error as estimated by sxcter is larger than defocuserror percent.  The error is computed as (std dev defocus)/defocus*100%. (Default: include all irrespective of error values.) (advanced)" )
	parser.add_option("--astigmatismerror", type="float",         default=360.0,     help="<Astigmatism Error Limit> Set to zero astigmatism for micrographs whose astigmatism angular error as estimated by sxcter is larger than astigmatismerror degrees. (Default: include all irrespective of error values.)  (advanced)")

	# Only for GUI support
	from optparse import OptionGroup
	# Add argument group to parser
	arg_group = OptionGroup(parser, "Arguments", "These options are interpretated as arguments by GUI.")

	arg_group.add_option("--sxgui_arguments",                                                                   help="GUI uses this option to get argument group.")
	# arg_group.add_option("--micrographs_list",                  type="string",		       default=".",                     help="<Micrograph List> List of input micrographs to be processed. ")

	parser.add_option_group(arg_group)
			
	# NOTE: 2015/10/22 Toshio Moriya
	# The followings are necessary because
	# some scripts support only MPI version but does not have --MPI option, and
	# some other scripts does not support MPI and does not have --MPI option.
			
	# Add MPI related option group to parser
	mpi_group = OptionGroup(parser, "MPI Options", "These options are used only by GUI.")
			
	mpi_group.add_option("--sxgui_mpi_options",                                        help="GUI uses this option to get MPI option group.")
	mpi_group.add_option("--MPI_support",        action="store_true",  default=False,  help="No --MPI option doesn't always mean that script does not support MPI.")
	mpi_group.add_option("--MPI_add_flag",       action="store_true",  default=False,  help="Need to add '--MPI' in command line.")
	
	parser.add_option_group(mpi_group)
	
	return parser