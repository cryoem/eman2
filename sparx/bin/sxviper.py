#!/usr/bin/env python

import sys
import os

import global_def
from global_def import *

def main(args):
	from utilities import if_error_all_processes_quit_program, write_text_row, drop_image, model_gauss_noise, get_im, set_params_proj, wrap_mpi_bcast, model_circle
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

	progname = os.path.basename(sys.argv[0])
	usage = progname + " stack  [output_directory] --ir=inner_radius --rs=ring_step --xr=x_range --yr=y_range  --ts=translational_search_step  --delta=angular_step --center=center_type --maxit1=max_iter1 --maxit2=max_iter2 --L2threshold=0.1 --ref_a=S --sym=c1"
	usage += """

stack			2D images in a stack file: (default required string)
directory		output directory name: into which the results will be written (if it does not exist, it will be created, if it does exist, the results will be written possibly overwriting previous results) (default required string)
"""
	
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--radius", type="int",            help="radius of the particle: has to be less than < int(nx/2)-1 (default required int)")

	parser.add_option("--xr", type="string",  default= '0',            help="range for translation search in x direction: search is +/xr in pixels (default '0')")
	parser.add_option("--yr", type="string",  default= '0',            help="range for translation search in y direction: if omitted will be set to xr, search is +/yr in pixels (default '0')")
	parser.add_option("--mask3D", type="string",  default= None,            help="3D mask file: (default sphere)")
	parser.add_option("--moon_elimination", type="string",  default= "",            help="elimination of disconnected pieces: two arguments: mass in KDa and pixel size in px/A separated by comma, no space (default none)")
	parser.add_option("--ir", type="int",  default= 1,            help="inner radius for rotational search: > 0 (default 1)")
	
	# 'radius' and 'ou' are the same as per Pawel's request; 'ou' is hidden from the user
	# the 'ou' variable is not changed to 'radius' in the 'sparx' program. This change is at interface level only for sxviper.
	##### XXXXXXXXXXXXXXXXXXXXXX option does not exist in docs XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	parser.add_option("--ou",       type= "int",   default= -1,                 help= SUPPRESS_HELP)
	parser.add_option("--rs", type="int",  default= 1,            help="step between rings in rotational search: >0 (default 1)")
	parser.add_option("--ts", type="string",  default= '1.0',            help="step size of the translation search in x-y directions: search is -xr, -xr+ts, 0, xr-ts, xr, can be fractional (default '1.0')")
	parser.add_option("--delta", type="string",  default= '2.0',            help="angular step of reference projections: (default '2.0')")
	parser.add_option("--center", type="float",  default= -1.0,            help="centering of 3D template: average shift method; 0: no centering; 1: center of gravity (default -1.0)")
	parser.add_option("--maxit1", type="float",  default= 400.0,            help="maximum number of iterations performed for the GA part: (default 400.0)")
	parser.add_option("--maxit2", type="float",  default= 50.0,            help="maximum number of iterations performed for the finishing up part: (default 50.0)")
	parser.add_option("--L2threshold", type="float",  default= 0.03,            help="stopping criterion of GA: given as a maximum relative dispersion of volumes' L2 norms: (default 0.03)")
	parser.add_option("--ref_a", type="string",  default= "S",            help="method for generating the quasi-uniformly distributed projection directions: (default S)")
	parser.add_option("--sym", type="string",  default= "c1",            help="point-group symmetry of the structure: (default c1)")
	
	# parser.add_option("--function", type="string", default="ref_ali3d",         help="name of the reference preparation function (ref_ali3d by default)")
	##### XXXXXXXXXXXXXXXXXXXXXX option does not exist in docs XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	parser.add_option("--function", type="string", default="ref_ali3d",         help= SUPPRESS_HELP)
	
	parser.add_option("--nruns", type="int",  default= 6,            help="GA population: aka number of quasi-independent volumes (default 6)")
	parser.add_option("--doga", type="float",  default= 0.1,            help="do GA when fraction of orientation changes less than 1.0 degrees is at least doga: (default 0.1)")
	##### XXXXXXXXXXXXXXXXXXXXXX option does not exist in docs XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	parser.add_option("--npad",     type="int",    default= 2,                  help="padding size for 3D reconstruction (default=2)")
	parser.add_option("--fl", type="float",  default= 0.25,            help="cut-off frequency applied to the template volume: using a hyperbolic tangent low-pass filter (default 0.25)")
	parser.add_option("--aa", type="float",  default= 0.1,            help="fall-off of hyperbolic tangent low-pass filter: (default 0.1)")
	parser.add_option("--pwreference", type="string",  default= "",            help="text file with a reference power spectrum: (default none)")
	parser.add_option("--debug", action="store_true",  default= False,            help="debug info printout: (default False)")
	
	##### XXXXXXXXXXXXXXXXXXXXXX option does not exist in docs XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	parser.add_option("--return_options", action="store_true", dest="return_options", default=False, help = SUPPRESS_HELP)	
	
	#parser.add_option("--an",       type="string", default= "-1",               help="NOT USED angular neighborhood for local searches (phi and theta)")
	#parser.add_option("--CTF",      action="store_true", default=False,         help="NOT USED Consider CTF correction during the alignment ")
	#parser.add_option("--snr",      type="float",  default= 1.0,                help="NOT USED Signal-to-Noise Ratio of the data (default 1.0)")
	# (options, args) = parser.parse_args(sys.argv[1:])

	required_option_list = ['radius']
	(options, args) = parser.parse_args(args)
	# option_dict = vars(options)
	# print parser
	
	if options.return_options:
		return parser
	
	if options.moon_elimination == "":
		options.moon_elimination = []
	else:
		options.moon_elimination = map(float, options.moon_elimination.split(","))

	# Making sure all required options appeared.
	for required_option in required_option_list:
		if not options.__dict__[required_option]:
			print "\n ==%s== mandatory option is missing.\n"%required_option
			print "Please run '" + progname + " -h' for detailed options"
			return 1



	if len(args) < 2 or len(args) > 3:
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
		return 1

	mpi_init(0, [])

	log = Logger(BaseLogger_Files())

	# 'radius' and 'ou' are the same as per Pawel's request; 'ou' is hidden from the user
	# the 'ou' variable is not changed to 'radius' in the 'sparx' program. This change is at interface level only for sxviper.
	options.ou = options.radius 
	runs_count = options.nruns
	mpi_rank = mpi_comm_rank(MPI_COMM_WORLD)
	mpi_size = mpi_comm_size(MPI_COMM_WORLD)	# Total number of processes, passed by --np option.

	if mpi_rank == 0:
		all_projs = EMData.read_images(args[0])
		subset = range(len(all_projs))
		# if mpi_size > len(all_projs):
		# 	ERROR('Number of processes supplied by --np needs to be less than or equal to %d (total number of images) ' % len(all_projs), 'sxviper', 1)
		# 	mpi_finalize()
		# 	return
	else:
		all_projs = None
		subset = None

	outdir = args[1]
	if mpi_rank == 0:
		if mpi_size % options.nruns != 0:
			ERROR('Number of processes needs to be a multiple of total number of runs. Total runs by default are 3, you can change it by specifying --nruns option.', 'sxviper', 1)
			mpi_finalize()
			return

		if os.path.exists(outdir):
			ERROR('Output directory exists, please change the name and restart the program', "sxviper", 1)
			mpi_finalize()
			return

		os.mkdir(outdir)
		import global_def
		global_def.LOGFILE =  os.path.join(outdir, global_def.LOGFILE)

	mpi_barrier(MPI_COMM_WORLD)

	if outdir[-1] != "/":
		outdir += "/"
	log.prefix = outdir
	
	# if len(args) > 2:
	# 	ref_vol = get_im(args[2])
	# else:
	ref_vol = None

	options.user_func = user_functions.factory[options.function]

	options.CTF = False
	options.snr = 1.0
	options.an  = -1.0

	out_params, out_vol, out_peaks = multi_shc(all_projs, subset, runs_count, options, mpi_comm=MPI_COMM_WORLD, log=log, ref_vol=ref_vol)

	mpi_finalize()

if __name__=="__main__":
	main(sys.argv[1:])

