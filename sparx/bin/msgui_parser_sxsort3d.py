#!/usr/bin/env python

from sparx import *
import os
import global_def
from   global_def import *
from   optparse  import OptionParser
import sys
from   numpy     import array
import types
from   logger    import Logger, BaseLogger_Files
				
def main(sys_arg_list):
	from logger import Logger, BaseLogger_Files
        arglist = []
        i = 0
        while( i < len(sys_arg_list) ):
            if sys_arg_list[i]=='-p4pg':
                i = i+2
            elif sys_arg_list[i]=='-p4wd':
                i = i+2
            else:
                arglist.append( sys_arg_list[i] )
                i = i+1
	progname = os.path.basename(arglist[0])
	usage = progname + " stack  outdir  <mask> --focus=3Dmask --ir=inner_radius --radius=outer_radius --rs=ring_step --xr=x_range --yr=y_range  --ts=translational_searching_step " +\
	" --delta=angular_step --an=angular_neighborhood --center=1 --nassign=reassignment_number --nrefine=alignment_number --maxit=max_iter --stoprnct=percentage_to_stop " + \
	" --CTF --snr=1.0 --ref_a=S --sym=c1 --function=user_function --independent=indenpendent_runs  --number_of_images_per_group=number_of_images_per_group  --resolution  "
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--focus",    type="string",       default=None,             help="<Focuse Clustering Mask Volume> 3D mask for focused clustering ")
	parser.add_option("--ir",       type= "int",         default=1, 	       help="<Inner Radius> inner radius for rotational correlation > 0 (set to 1). (advanced)")
	parser.add_option("--radius",   type= "int",         default="-1",	       help="<Outer Radius> outer radius for rotational correlation <nx-1 (set to the radius of the particle). (advanced)")
	parser.add_option("--maxit",	type= "int",         default=5, 	       help="<Iterations> maximum number of iteration")
	parser.add_option("--rs",       type= "int",         default="1",	       help="<Ring Step> step between rings in rotational correlation >0 (set to 1) (advanced)" ) 
	parser.add_option("--xr",       type="string",       default="4 2 1 1 1",      help="<X Search Range> range for translation search in x direction, search is +/-xr ")
	parser.add_option("--yr",       type="string",       default="-1",	       help="<Y Search Range> range for translation search in y direction, search is +/-yr (default = same as xr)")
	parser.add_option("--ts",       type="string",       default="0.25",           help="<Search Step Size> step size of the translation search in both directions direction, search is -xr, -xr+ts, 0, xr-ts, xr ")
	parser.add_option("--delta",    type="string",       default="10 6 4  3   2",  help="<Angular Step> angular step of reference projections")
	parser.add_option("--an",       type="string",       default="-1",	       help="<Angular Neighbors> angular neighborhood for local searches")
	parser.add_option("--center",   type="int",          default=0,	               help="<Centering Method> 0 - if you do not want the volume to be centered, 1 - center the volume using cog (default=0) (advanced)")
	parser.add_option("--nassign",  type="int",          default=0, 	       help="<Reassignment Iterations> number of reassignment iterations performed for each angular step (set to 3). (advanced)")
	parser.add_option("--nrefine",  type="int",          default=1, 	       help="<Alingment Iterations> number of alignment iterations performed for each angular step (set to 1). (advanced)")
	parser.add_option("--CTF",      action="store_true", default=False,            help="<Use CTF Correction> Consider CTF correction during the alignment ")
	parser.add_option("--snr",      type="float",        default=1.0,              help="<Signal to Noize Ratio> Signal-to-Noise Ratio of the data")
	parser.add_option("--stoprnct", type="float",        default=0.0,              help="<Termination Change Percentage> Minimum percentage of assignment change to stop the program")   
	parser.add_option("--ref_a",    type="string",       default="S",              help="<Angular Sampling Method> method for generating the quasi-uniformly distributed projection directions (default S) (advanced)")
	parser.add_option("--sym",      type="string",       default="c1",             help="<Point-Group Symmetry> symmetry of the structure ")
	parser.add_option("--function", type="string",       default="ref_ali3dm",     help="<User Function Name> name of the reference preparation function (advanced)")
	parser.add_option("--npad",     type="int",          default= 2,               help="<Padding Size> padding size for 3D reconstruction (advanced)")
	parser.add_option("--independent", type="int",       default= 3,               help="<Number of Indepndent Run> number of independent run (advanced)")
	parser.add_option("--number_of_images_per_group",    type="int",               default=-1,               help="<Number of Images per Groups> number of groups")
	parser.add_option("--resolution",  type="float",     default= .40,             help="<Resolution Limit> structure is low-pass-filtered to this resolution for clustering (advanced)" )
	parser.add_option("--mode",        type="string",    default="EK_only",        help="<Clusterning Mode> mode options: EK_only, Kgroup_guess, auto_search, lpf_search,large_number_run,isac. (advanced)" )
	#parser.add_option("--importali3d", type="string",    default="",               help="import the xform.projection parameters as the initial configuration for 3-D reconstruction" )
	#parser.add_option("--Kgroup_guess",  action="store_true",default=False,        help="Guess the possible number of groups existing in one dataset" )
	#parser.add_option("--frequency_start_search",  type="float",default=.10,       help="start frequency for low pass filter search")
	#parser.add_option("--frequency_stop_search",   type="float",default=.40,       help="stop frequency for low pass filter search")
	#parser.add_option("--frequency_search_step",   type="float",default=.02,       help="frequency step for low pass filter search")
	parser.add_option("--scale_of_number", type="float",default=1.,                 help="<Scale number> scale number to control particle number per group")
	
	# Only for GUI support
	from optparse import OptionGroup
	# Add argument group to parser
	arg_group = OptionGroup(parser, "Arguments", "These options are interpretated as arguments by GUI.")
	
	arg_group.add_option("--sxgui_arguments",                     help="GUI uses this option to get argument group.")
	arg_group.add_option("--stack",   type="string", default="",  help="<Particle Stack> Set of 2-D images in a stack file.")
	arg_group.add_option("--outdir",  type="string", default="",  help="<Output Directory> Directory name into which the results will be written (if it does not exist, it will be created, if it does exist, the results will be written possibly overwriting previous results).")
	arg_group.add_option("--mask",    type="string", default="",  help="<Mask Volume> input file name of 3D mask.")

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
