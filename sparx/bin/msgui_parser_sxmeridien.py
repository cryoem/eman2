#!/usr/bin/env python

from __future__ import print_function
from EMAN2 import *
from sparx import *
from logger import Logger, BaseLogger_Files
import global_def

from mpi   import  *
from math  import  *

import os
import sys
import subprocess
import time
import string
from   sys import exit
from   time import localtime, strftime

def main(sys_arg_list):

	from utilities import write_text_row, drop_image, model_gauss_noise, get_im, set_params_proj, wrap_mpi_bcast, model_circle, get_shrink_data
	import user_functions
	from applications import MPI_start_end
	from optparse import OptionParser
	from global_def import SPARXVERSION
	from EMAN2 import EMData
	from multi_shc import multi_shc
	#from development import do_volume_mrk01
	from logger import Logger, BaseLogger_Files
	import sys
	import os
	import time
	import socket
	
	# ------------------------------------------------------------------------------------
	# PARSE COMMAND OPTIONS
	progname = os.path.basename(sys_arg_list[0])
	usage = progname + " stack  [output_directory]  initial_volume  --radius=particle_radius --ref_a=S --sym=c1 --startangles --inires  --mask3D --CTF --function=user_function"
	parser = OptionParser(usage,version=SPARXVERSION)
	#parser.add_option("--ir",      		type= "int",   default= 1,			help="inner radius for rotational correlation > 0 (set to 1)")
	parser.add_option("--radius",      		type= "int",   default= -1,			help="<Particle Radius> Outer radius [in pixels] for rotational correlation < int(nx/2)-1 (Please set to the radius of the particle)")
	##parser.add_option("--rs",      		type= "int",   default= 1,			help="step between rings in rotational correlation >0  (set to 1)" ) 
	#parser.add_option("--xr",      		type="string", default= "-1",		help="range for translation search in x direction, search is +/xr (default 0)")
	#parser.add_option("--yr",      		type="string", default= "-1",		help="range for translation search in y direction, search is +/yr (default = same as xr)")
	#parser.add_option("--ts",      		type="string", default= "1",		help="step size of the translation search in both directions, search is -xr, -xr+ts, 0, xr-ts, xr, can be fractional")
	#parser.add_option("--delta",   		type="string", default= "-1",		help="angular step of reference projections during initialization step (default automatically selected based on radius of the structure.)")
	#parser.add_option("--an",      		type="string", default= "-1",		help="angular neighborhood for local searches (phi and theta) (Default exhaustive searches)")
	#parser.add_option("--center",  		type="int",  default= 0,			help="-1: average shift method; 0: no centering; 1: center of gravity (default=0)")
	#parser.add_option("--maxit",   		type="int",  	default= 400,		help="maximum number of iterations performed for the GA part (set to 400) ")
	parser.add_option("--outlier_percentile",type="float",    default= 95,	help="percentile above which outliers are removed every iteration")
	#parser.add_option("--iteration_start",   type="int",    default= 0,		help="starting iteration for rviper, 0 means go to the most recent one (default).")
	parser.add_option("--CTF",     		     action="store_true", default=False,	help="<Use CTF Correction> Use CTF (Default no CTF correction)")
	#parser.add_option("--snr",     		type="float",  default= 1.0,		help="Signal-to-Noise Ratio of the data (default 1.0)")
	parser.add_option("--ref_a",   		    type="string", default= "S",		help="method for generating the quasi-uniformly distributed projection directions (default S)")
	parser.add_option("--sym",     		    type="string", default= "c1",		help="<Point-Group Symmetry> Point-group symmetry of the refined structure")
	#parser.add_option("--npad",    		type="int",    default= 2,			help="padding size for 3D reconstruction (default=2)")
	#parser.add_option("--nsoft",    	     type="int",    default= 0,			help="Use SHC in first phase of refinement iteration (default=0, to turn it on set to 1)")
	parser.add_option("--startangles",      action="store_true", default=False,	help="Use orientation parameters in the input file header to jumpstart the procedure")
	parser.add_option("--restrict_shifts",  type="int",    default= -1,			help="Restrict initial searches for translation [unit - original size pixel] (default=-1, no restriction)")
	parser.add_option("--local_filter",     action="store_true", default=False,	help="Use local filtration (Default generic tangent filter)")
	parser.add_option("--smear",            action="store_true", default=False,	help="Use rotational smear")
	parser.add_option("--sausage",          action="store_true", default=False,	help="Sausage-making filter")

	#options introduced for the do_volume function
	#parser.add_option("--fl",			type="float",	default=0.12,		help="cut-off frequency of hyperbolic tangent low-pass Fourier filter (default 0.12)")
	#parser.add_option("--aa",			type="float",	default=0.1,		help="fall-off of hyperbolic tangent low-pass Fourier filter (default 0.1)")
	parser.add_option("--inires",		     type="float",	default=25.,		help="<Initial Resolution> Resolution of the initial_volume volume (default 25A)")
	parser.add_option("--pwreference",	     type="string",	default="",			help="<Power Spectrum> text file with a reference power spectrum (default no power spectrum adjustment)")
	parser.add_option("--mask3D",		     type="string",	default=None,		help="<3D Mask> 3D mask file (default a sphere with radius (nx/2)-1)")
	parser.add_option("--function",          type="string", default="do_volume_mrk02",  help="<User Function Name> name of the reference preparation function (default do_volume_mrk02) (advanced)")

	# Only for GUI support
	from optparse import OptionGroup
	# Add argument group to parser
	arg_group = OptionGroup(parser, "Arguments", "These options are interpretated as arguments by GUI.")

	arg_group.add_option("--sxgui_arguments",                                                                   help="GUI uses this option to get argument group.")
	arg_group.add_option("--stack",                  type="string",		       default=".",                     help="<Particle Stack> Set of 2-D images in a stack file. ")
	arg_group.add_option("--output_directory",       type="string",		       default="",                      help="<Output Directory> Output Directory to store all outputs.")
	arg_group.add_option("--initial_volume",         type="string",		       default="",                      help="<Initial 3D Structure> 3D Structure.")

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
