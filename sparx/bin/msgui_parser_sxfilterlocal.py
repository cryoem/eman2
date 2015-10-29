#!/usr/bin/env python

import global_def
from   global_def import *
from   EMAN2 import *
from   sparx import *
from optparse import OptionParser

def main(sys_arg_list):
	progname = os.path.basename(sys_arg_list[0])
	usage = progname + """ inputvolume  locresvolume maskfile outputfile   --radius --falloff  --MPI

	    Locally filer a volume based on local resolution volume (sxlocres.py) within area outlined by the maskfile
	"""
	parser = OptionParser(usage,version=SPARXVERSION)

	parser.add_option("--radius",	type="int",		        default=-1, 	help="<Particle Radius> if there is no maskfile, sphere with r=radius will be used, by default the radius is nx/2-1")
	parser.add_option("--falloff",	type="float",		    default=0.1,    help="<Falloff of Tangent Filter> falloff of tangent low-pass filter (tanl) (default 0.1)")
	parser.add_option("--MPI",      action="store_true",   	default=False,  help="use MPI version")

	# Only for GUI support
	from optparse import OptionGroup
	# Add argument group to parser
	arg_group = OptionGroup(parser, "Arguments", "These options are interpretated as arguments by GUI.")

	arg_group.add_option("--sxgui_arguments",                            help="GUI uses this option to get argument group.")
	arg_group.add_option("--inputvolume",   type="string",  default="",  help="<First Half-set Volume> input file name of first half-set volume.")
	arg_group.add_option("--locresvolume",  type="string",  default="",  help="<Local Resolution Volume> input file name of local resolution volume (as produced by sxlocres.py).")
	arg_group.add_option("--maskfile",      type="string",  default="",  help="<Mask Volume> input file name of mask volume outlining the region within which local filtration will be applied. (optional)")
	arg_group.add_option("--outputfile",    type="string",  default="",  help="<Locally Filtered Volume> output file name of locally filtered volume.")

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
