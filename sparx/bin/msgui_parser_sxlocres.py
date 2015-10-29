#!/usr/bin/env python

import global_def
from   global_def import *
from EMAN2 import *
from sparx import *
from optparse import OptionParser

def main(sys_arg_list):
	progname = os.path.basename(sys_arg_list[0])
	usage = progname + """ firstvolume  secondvolume maskfile outputfile --wn --step --cutoff  --radius  --fsc --MPI

	Compute local resolution in real space within area outlined by the maskfile and within regions wn x wn x wn
	"""
	parser = OptionParser(usage,version=SPARXVERSION)
	
	parser.add_option("--wn",		type="int",		default=7, 			help="<Window Size in Voxels> Size of window within which local real-space FSC is computed (default 7")
	parser.add_option("--step",     type="float",	default= 1.0,       help="<Resolution Step in Voxels> Shell step in Fourier size in pixels (default 1.0)")   
	parser.add_option("--cutoff",   type="float",	default= 0.5,       help="<Resolution Cutoff> resolution cut-off for FSC (default 0.5)")
	parser.add_option("--radius",	type="int",		default=-1, 		help="<Particule Radius> if there is no maskfile, sphere with r=radius will be used, by default the radius is nx/2-wn")
	parser.add_option("--fsc",      type="string",	default= None,      help="<FSC File Name> overall FSC curve (might be truncated) (default no curve)")
	parser.add_option("--MPI",      action="store_true",   	default=False,  help="use MPI version")

	# Only for GUI support
	from optparse import OptionGroup
	# Add argument group to parser
	arg_group = OptionGroup(parser, "Arguments", "These options are interpretated as arguments by GUI.")

	arg_group.add_option("--sxgui_arguments",                            help="GUI uses this option to get argument group.")
	arg_group.add_option("--firstvolume",   type="string",  default="",  help="<First Half-set Volume> input file name of first half-set volume.")
	arg_group.add_option("--secondvolume",  type="string",	default="",  help="<Second Half-set Volume> input file name of second half-set volume.")
	arg_group.add_option("--maskfile",      type="string",  default="",  help="<Mask Volume> input file name of 3D mask.")
	arg_group.add_option("--outputfile",    type="string",  default="",  help="<Local Resolution Volume> output file name of local resolution volume.")

	parser.add_option_group(arg_group)

	
	# NOTE: 2015/10/22 Toshio Moriya
	# The followings are necessary because
	# some scripts support only MPI version but does not have --MPI option, and
	# some other scripts does not support MPI and does not have --MPI option.
			
	# Add MPI related option group to parser
	mpi_group = OptionGroup(parser, "MPI Options", "These options are used only by GUI.")
			
	mpi_group.add_option("--sxgui_mpi_options",                                       help="GUI uses this option to get MPI option group.")
	mpi_group.add_option("--MPI_support",        action="store_true",  default=True,  help="No --MPI option doesn't always mean that script does not support MPI.")
	mpi_group.add_option("--MPI_add_flag",       action="store_true",  default=True,  help="Need to add '--MPI' in command line.")
	
	parser.add_option_group(mpi_group)
	
	return parser
