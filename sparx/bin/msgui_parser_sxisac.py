#!/usr/bin/env python

import os
import sys
import random

import global_def
from   global_def import *
# from optparse import OptionParser, SUPPRESS_HELP
from optparse import OptionParser
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

NAME_OF_JSON_STATE_FILE = "my_state.json"
NAME_OF_ORIGINAL_IMAGE_INDEX = "originalid"
NAME_OF_RUN_DIR = "run"
NAME_OF_MAIN_DIR = "generation_"
DIR_DELIM = os.sep

def main(sys_arg_list):
	progname = os.path.basename(sys_arg_list[0])
	usage = ( progname + " stack_file  output_directory --radius=particle_radius --img_per_grp=img_per_grp --CTF --restart_section<The remaining parameters are optional --ir=ir --rs=rs --xr=xr --yr=yr --ts=ts --maxit=maxit --dst=dst --FL=FL --FH=FH --FF=FF --init_iter=init_iter --main_maxit=main_iter" +
			" --iter_reali=iter_reali --match_first=match_first --max_round=max_round --match_second=match_second --stab_ali=stab_ali --thld_err=thld_err --indep_run=indep_run --thld_grp=thld_grp" +
			"  --generation=generation  --rand_seed=rand_seed>" )
	
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--radius",         type="int",          default=-1,      help="<Particle Radius> It has to be provided.")
	parser.add_option("--img_per_grp",    type="int",          default=100,     help="<Number of Images Per Group> in the ideal case (essentially maximum size of class) (100)")
	parser.add_option("--CTF",            action="store_true", default=False,   help="<CTF Flag>, if set the data will be phase-flipped")
	parser.add_option("--ir",             type="int",          default=1,       help="<Inner Ring> of the resampling to polar coordinates (1)")
	parser.add_option("--rs",             type="int",          default=1,       help="<Ring Step> of the resampling to polar coordinates (1)")
	parser.add_option("--xr",             type="int",          default=-1,      help="<X Range> of translational search (By default set by the program) (advanced)")
	parser.add_option("--yr",             type="int",          default=-1,      help="<Y Range> of translational search (same as xr) (advanced)")
	parser.add_option("--ts",             type="float",        default=1.0,     help="<Search Step> of translational search (1.0)")
	parser.add_option("--maxit",          type="int",          default=30,      help="number of iterations for reference-free alignment (30)")
	#parser.add_option("--snr",            type="float",        default=1.0,     help="signal-to-noise ratio (only meaningful when CTF is enabled, currently not supported)")
	parser.add_option("--center_method",  type="int",          default=7,       help="<Centering Method> of global 2D average during initial prealignment of data (default : 7; 0 : no centering; -1 : average shift method; please see center_2D in utilities.py for methods 1-7)")
	parser.add_option("--dst",            type="float",        default=90.0,    help="discrete angle used in within group alignment ")
	parser.add_option("--FL",             type="float",        default=0.2,     help="<Lowest Stopband> frequency used in the tangent filter (0.2)")
	parser.add_option("--FH",             type="float",        default=0.3,     help="<Highest Stopband> frequency used in the tangent filter (0.3)")
	parser.add_option("--FF",             type="float",        default=0.2,     help="<Tangent Filter Falloff> falloff of tangent filter (0.2)")
	parser.add_option("--init_iter",      type="int",          default=3,       help="<Init Iterations> number of iterations of ISAC program in initialization (3)")
	parser.add_option("--main_iter",      type="int",          default=3,       help="<Main Iterations> number of iterations of ISAC program in main part (3)")
	parser.add_option("--iter_reali",     type="int",          default=1,       help="<Realignment Iterations> number of iterations in ISAC before checking stability (1)")
	parser.add_option("--match_first",    type="int",          default=1,       help="number of iterations to run 2-way matching in the first phase (1)")
	parser.add_option("--max_round",      type="int",          default=20,      help="maximum rounds of generating candidate averages in the first phase (20)")
	parser.add_option("--match_second",   type="int",          default=5,       help="number of iterations to run 2-way (or 3-way) matching in the second phase (5)")
	parser.add_option("--stab_ali",       type="int",          default=5,       help="number of alignments when checking stability (5)")
	parser.add_option("--thld_err",       type="float",        default=0.7,     help="the threshold of pixel error when checking stability (0.7)")
	parser.add_option("--indep_run",      type="int",          default=4,       help="number of independent runs for reproducibility (default=4, only values 2, 3 and 4 are supported (4)")
	parser.add_option("--thld_grp",       type="int",          default=10,      help="minimum size of class (10)")
	parser.add_option("--n_generations",     type="int",          default=100,       help="<Number of Generations> program stops when reaching this total number of generations (advanced)")
	#parser.add_option("--candidatesexist",action="store_true", default=False,   help="Candidate class averages exist use them (default False)")
	parser.add_option("--rand_seed",      type="int",          default=None,    help="random seed set before calculations, useful for testing purposes (default None - total randomness)")
	parser.add_option("--new",            action="store_true", default=False,   help="use new code (default = False)")
	parser.add_option("--debug",          action="store_true", default=False,   help="debug info printout (default = False)")

	# must be switched off in production
	parser.add_option("--use_latest_master_directory", action="store_true", dest="use_latest_master_directory", default=False)
	
	parser.add_option("--restart_section", type="string", default="", help="<Restart Section Name> (no spaces) followed immediately by comma, followed immediately by generation to restart, example: \n--restart_section=candidate_class_averages,1         (Sections: restart, candidate_class_averages, reproducible_class_averages)")
	parser.add_option("--stop_after_candidates",          action="store_true", default=False,   help="<Stop After Candidates> stops after the 'candidate_class_averages' section")
	
	# Only for GUI support
	from optparse import OptionGroup
	# Add argument group to parser
	arg_group = OptionGroup(parser, "Arguments", "These options are interpretated as arguments by GUI.")
	
	arg_group.add_option("--sxgui_arguments",                              help="GUI uses this option to get argument group.")
	arg_group.add_option("--stack_file",       type="string", default=".", help="<Particle Stack> Set of 2-D images in a stack file (format must be bdb), images have to be square (nx=ny).")
	arg_group.add_option("--output_directory", type="string", default="",  help="<Output Directory> Directory name into which the results will be written (if it does not exist, it will be created, if it does exist, the results will be written possibly overwriting previous results).")

	parser.add_option_group(arg_group)
	
	# NOTE: 2015/10/22 Toshio Moriya
	# The followings are necessary because
	# some scripts support only MPI version but does not have --MPI option, and
	# some other scripts does not support MPI and does not have --MPI option.
			
	# Add MPI related option group to parser
	mpi_group = OptionGroup(parser, "MPI Options", "These options are used only by GUI.")
			
	mpi_group.add_option("--sxgui_mpi_options",                                        help="GUI uses this option to get MPI option group.")
	mpi_group.add_option("--MPI_support",        action="store_true",  default=True,   help="No --MPI option doesn't always mean that script does not support MPI.")
	mpi_group.add_option("--MPI_add_flag",       action="store_true",  default=False,  help="Need to add '--MPI' in command line.")
	
	parser.add_option_group(mpi_group)
	
	return parser



