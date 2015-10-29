#!/usr/bin/env python

from EMAN2 import *
from sparx import *


import	global_def
from	global_def 	import *
from	optparse 	import OptionParser
from	EMAN2 		import EMUtil
import	os
import	sys
from 	time		import	time


def main(sys_arg_list):

	def params_3D_2D_NEW(phi, theta, psi, s2x, s2y, mirror):
		if mirror:
			m = 1
			alpha, sx, sy, scalen = compose_transform2(0, s2x, s2y, 1.0, 540.0-psi, 0, 0, 1.0)
		else:
			m = 0
			alpha, sx, sy, scalen = compose_transform2(0, s2x, s2y, 1.0, 360.0-psi, 0, 0, 1.0)
		return  alpha, sx, sy, m
	
	progname = os.path.basename(sys_arg_list[0])
	usage = progname + " prj_stack  --ave2D= --var2D=  --ave3D= --var3D= --img_per_grp= --fl=0.2 --aa=0.1  --sym=symmetry --CTF"
	parser = OptionParser(usage, version=SPARXVERSION)

	parser.add_option("--ave2D",		type="string"	   ,	default=False,				help="<Save 2D Averages> write to the disk a stack of 2D averages")
	parser.add_option("--var2D",		type="string"	   ,	default=False,				help="<Save 2D Variances> write to the disk a stack of 2D variances")
	parser.add_option("--ave3D",		type="string"	   ,	default=False,				help="<Save 3D Average> write to the disk reconstructed 3D average")
	parser.add_option("--var3D",		type="string"	   ,	default=False,				help="<Save 3D Variance> compute 3D variability (time consuming!)")
	parser.add_option("--img_per_grp",	type="int"         ,	default=10   ,				help="<Number of Neighbour Projections> number of neighbouring projections")
	parser.add_option("--no_norm",		action="store_true",	default=False,				help="<Disable Normalization> do not use normalization")
	parser.add_option("--radiusvar", 	type="int"         ,	default=-1   ,				help="<Radius of 3D Variance> radius for 3D var" )
	parser.add_option("--npad",			type="int"         ,	default=2    ,				help="<Padding Size> number of time to pad the original images (advanced)")
	parser.add_option("--sym" , 		type="string"      ,	default="c1" ,				help="<Point-Group Symmetry> Point-group symmetry")
	parser.add_option("--fl",			type="float"       ,	default=0.0  ,				help="<Stop-Band Frequency of Filter> stop-band frequency of low pass filter (Default - no filtration) (advanced)")
	parser.add_option("--aa",			type="float"       ,	default=0.0  ,				help="<Falloff of Filter> fall off of the local pass filter (Default - no filtration) (advanced)")
	parser.add_option("--CTF",			action="store_true",	default=False,				help="<Use CTF Correction> use CFT correction")
	parser.add_option("--VERBOSE",		action="store_true",	default=False,				help="Long output for debugging")
	#parser.add_option("--MPI" , 		action="store_true",	default=False,				help="use MPI version")

	#parser.add_option("--radiuspca", 	type="int"         ,	default=-1   ,				help="radius for PCA" )
	#parser.add_option("--iter", 		type="int"         ,	default=40   ,				help="maximum number of iterations (stop criterion of reconstruction process)" )
	#parser.add_option("--abs", 			type="float"       ,	default=0.0  ,				help="minimum average absolute change of voxels' values (stop criterion of reconstruction process)" )
	#parser.add_option("--squ", 			type="float"       ,	default=0.0  ,				help="minimum average squared change of voxels' values (stop criterion of reconstruction process)" )
	parser.add_option("--VAR" , 		action="store_true",	default=False,				help="<VAR> stack on input consists of 2D variances (Default False) (advanced)")
	parser.add_option("--SND",			action="store_true",	default=False,				help="<SND> compute squared normalized differences (Default False) (advanced)")
	#parser.add_option("--nvec",			type="int"         ,	default=0    ,				help="number of eigenvectors, default = 0 meaning no PCA calculated")
	parser.add_option("--symmetrize",	action="store_true",	default=False,				help="<Symmetrize> Prepare input stack for handling symmetry (Default False) (advanced)")

	# Only for GUI support
	from optparse import OptionGroup
	# Add argument group to parser
	arg_group = OptionGroup(parser, "Arguments", "These options are interpretated as arguments by GUI.")

	arg_group.add_option("--sxgui_arguments",                               help="GUI uses this option to get argument group.")
	arg_group.add_option("--prj_stack",        type="string",  default="",  help="<Stack File> Stack of 2D images with 3D orientation parameters in header and (optionally) CTF information.")

	parser.add_option_group(arg_group)

	
	# NOTE: 2015/10/22 Toshio Moriya
	# The followings are necessary because
	# some scripts support only MPI version but does not have --MPI option, and
	# some other scripts does not support MPI and does not have --MPI option.
			
	# Add MPI related option group to parser
	mpi_group = OptionGroup(parser, "MPI Options", "These options are used only by GUI.")
			
	mpi_group.add_option("--sxgui_mpi_options",                                        help="GUI uses this option to get MPI option group.")
	mpi_group.add_option("--MPI_support",        action="store_true",  default=True,   help="No --MPI option doesn't always mean that script does not support MPI.")
	mpi_group.add_option("--MPI_add_flag",       action="store_true",  default=False,  help="Need to add --MPI in command line.")
	
	parser.add_option_group(mpi_group)
	
	return parser

