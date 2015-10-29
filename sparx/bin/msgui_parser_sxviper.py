#!/usr/bin/env python

import sys
import os

import global_def
from global_def import *

def main(sys_arg_list):
	from utilities import write_text_row, drop_image, model_gauss_noise, get_im, set_params_proj, wrap_mpi_bcast, model_circle
	from logger import Logger, BaseLogger_Files
	from mpi import mpi_init, mpi_finalize, MPI_COMM_WORLD, mpi_comm_rank, mpi_comm_size, mpi_barrier
	import user_functions
	import sys
	import os
	from applications import MPI_start_end
	from optparse import OptionParser, SUPPRESS_HELP
	from global_def import SPARXVERSION
	from EMAN2 import EMData
	from multi_shc import multi_shc

	progname = os.path.basename(sys_arg_list[0])
	usage = progname + " stack  output_directory  --ir=inner_radius --ou=outer_radius --rs=ring_step --xr=x_range --yr=y_range  --ts=translational_search_step  --delta=angular_step --center=center_type --maxit1=max_iter1 --maxit2=max_iter2 --L2threshold=0.1 --ref_a=S --sym=c1"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--ir",       type= "int",   default= 1,                  help="<Inner Radius> for rotational correlation > 0 (set to 1)")
	parser.add_option("--ou",       type= "int",   default= -1,                 help="<Outer Radius> for rotational correlation < int(nx/2)-1 (set to the radius of the particle)")
	parser.add_option("--rs",       type= "int",   default= 1,                  help="<Step Between> rings in rotational correlation >0  (set to 1)" ) 
	parser.add_option("--xr",       type="string", default= "0",                help="<X Search Range> for translation search in x direction, search is +/xr (default 0)")
	parser.add_option("--yr",       type="string", default= "-1",               help="<Y Search Range> for translation search in y direction, search is +/yr (default = same as xr)")
	parser.add_option("--ts",       type="string", default= "1",                help="<Search Step Size> of the translation search in both directions, search is -xr, -xr+ts, 0, xr-ts, xr, can be fractional")
	parser.add_option("--delta",    type="string", default= "2",                help="<Angular Step> of reference projections (default 2)")
	parser.add_option("--center",   type="float",  default= -1,                 help="-1: average shift method; 0: no centering; 1: center of gravity (default=-1)")
	parser.add_option("--maxit1",   type="float",  default= 400,                help="maximum number of iterations performed for the GA part (set to 400) ")
	parser.add_option("--maxit2",   type="float",  default= 50,                 help="maximum number of iterations performed for the finishing up part (set to 50) ")
	parser.add_option("--L2threshold", type="float",  default= 0.03,            help="Stopping criterion of GA given as a maximum relative dispersion of L2 norms (set to 0.03) ")
	parser.add_option("--ref_a",    type="string", default= "S",                help="method for generating the quasi-uniformly distributed projection directions (default S)")
	parser.add_option("--sym",      type="string", default= "c1",               help="<Symmetry> of the refined structure")
	
	# parser.add_option("--function", type="string", default="ref_ali3d",         help="name of the reference preparation function (ref_ali3d by default)")
	parser.add_option("--function", type="string", default="ref_ali3d",         help= SUPPRESS_HELP)
	
	parser.add_option("--nruns",    type="int",    default= 6,                  help="number of quasi-independent runs (default=6)")
	parser.add_option("--doga",     type="float",  default= 0.1,                help="do GA when fraction of orientation changes less than 1.0 degrees is at least doga (default=0.1)")
	parser.add_option("--npad",     type="int",    default= 2,                  help="padding size for 3D reconstruction (default=2)")
	parser.add_option("--fl",       type="float",  default=0.25,                help="<Cutoff Frequency> of hyperbolic tangent low-pass Fourier filter (default 0.25)")
	parser.add_option("--aa",       type="float",  default=0.1,                 help="<Falloff Frequency> of hyperbolic tangent low-pass Fourier filter (default 0.1)")
	parser.add_option("--pwreference",      type="string",  default="",         help="<Power Spectrum> reference text file (default no power spectrum adjustment) (advanced)")
	parser.add_option("--mask3D",      type="string",  default=None,            help="3D mask file (default a sphere)")
	parser.add_option("--moon_elimination",      type="string",  default="",    help="<Moon Elimination> mass in KDa and resolution in px/A separated by comma, no space (advanced)")
	parser.add_option("--debug",          action="store_true", default=False,   help="<Debug> info printout (default = False)")
		
	#parser.add_option("--an",       type="string", default= "-1",               help="NOT USED angular neighborhood for local searches (phi and theta)")
	#parser.add_option("--CTF",      action="store_true", default=False,         help="NOT USED Consider CTF correction during the alignment ")
	#parser.add_option("--snr",      type="float",  default= 1.0,                help="NOT USED Signal-to-Noise Ratio of the data (default 1.0)")

	# Only for GUI support
	from optparse import OptionGroup
	# Add argument group to parser
	arg_group = OptionGroup(parser, "Arguments", "These options are interpretated as arguments by GUI.")
	
	arg_group.add_option("--sxgui_arguments",                              help="GUI uses this option to get argument group.")
	arg_group.add_option("--stack",            type="string", default=".", help="<Projection Stack File> Set of 2-D images in a stack file (format must be bdb), images have to be square (nx=ny).")
	arg_group.add_option("--output_directory", type="string", default="",  help="<Output Directory> Directory name into which the results will be written (if it does not exist, it will be created, if it does exist, the results will be written possibly overwriting previous results)")

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

