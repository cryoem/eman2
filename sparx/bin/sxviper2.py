#!/usr/bin/env python

import global_def
from global_def import *

def main():
	from utilities import write_text_row, drop_image, model_gauss_noise, get_im, set_params_proj, wrap_mpi_bcast, model_circle
	from logger import Logger, BaseLogger_Files
	from mpi import mpi_init, mpi_finalize, MPI_COMM_WORLD, mpi_comm_rank, mpi_comm_size, mpi_barrier
	import sys
	import os
	import user_functions
	from applications import MPI_start_end
	from optparse import OptionParser
	from global_def import SPARXVERSION
	from EMAN2 import EMData
	from multi_shc import multi_shc

	progname = os.path.basename(sys.argv[0])
	usage = progname + " stack  output_directory  [initial_volume]  --ir=inner_radius --ou=outer_radius --rs=ring_step --xr=x_range --yr=y_range  --ts=translational_search_step  --delta=angular_step --an=angular_neighborhood  --center=center_type --maxit1=max_iter1 --maxit2=max_iter2 --L2threshold=0.1  --CTF --snr=SNR  --ref_a=S --sym=c1 --function=user_function"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--ir",       type= "int",   default= 1,                  help="inner radius for rotational correlation > 0 (set to 1)")
	parser.add_option("--ou",       type= "int",   default= -1,                 help="outer radius for rotational correlation < int(nx/2)-1 (set to the radius of the particle)")
	parser.add_option("--rs",       type= "int",   default= 1,                  help="step between rings in rotational correlation >0  (set to 1)" ) 
	parser.add_option("--xr",       type="string", default= "0",                help="range for translation search in x direction, search is +/xr (default 0)")
	parser.add_option("--yr",       type="string", default= "-1",               help="range for translation search in y direction, search is +/yr (default = same as xr)")
	parser.add_option("--ts",       type="string", default= "1",                help="step size of the translation search in both directions, search is -xr, -xr+ts, 0, xr-ts, xr, can be fractional")
	parser.add_option("--delta",    type="string", default= "2",                help="angular step of reference projections (default 2)")
	parser.add_option("--an",       type="string", default= "-1",               help="angular neighborhood for local searches (phi and theta)")
	parser.add_option("--center",   type="float",  default= -1,                 help="-1: average shift method; 0: no centering; 1: center of gravity (default=-1)")
	parser.add_option("--maxit1",    type="float",  default= 400,               help="maximum number of iterations performed for the GA part (set to 400) ")
	parser.add_option("--maxit2",    type="float",  default= 30,                help="maximum number of iterations performed for the finishing up part (set to 30) ")
	parser.add_option("--L2threshold", type="float",  default= 0.05,            help="Stopping criterion of GA given as a maximum relative dispersion of L2 norms (set to 0.05) ")
	parser.add_option("--CTF",      action="store_true", default=False,         help="Consider CTF correction during the alignment ")
	parser.add_option("--snr",      type="float",  default= 1.0,                help="Signal-to-Noise Ratio of the data (default 1.0)")
	parser.add_option("--ref_a",    type="string", default= "S",                help="method for generating the quasi-uniformly distributed projection directions (default S)")
	parser.add_option("--sym",      type="string", default= "c1",               help="symmetry of the refined structure")
	parser.add_option("--function", type="string", default="ref_ali3d",         help="name of the reference preparation function (ref_ali3d by default)")
	parser.add_option("--nruns",    type="int",    default= 3,                  help="number of quasi-independent runs (default=3)")
	parser.add_option("--niruns",       type= "int",   default= 1,                  help="number of independent runs")
	parser.add_option("--doga",     type="float",  default= 0.3,                help="do GA when fraction of orientation changes less than 1.0 degrees is at least doga (default=0.3)")
	parser.add_option("--npad",     type="int",    default= 2,                  help="padding size for 3D reconstruction (default=2)")
	#parser.add_option("--MPI",      action="store_true", default=True,          help="whether to use MPI version - this is always set to True")
	(options, args) = parser.parse_args(sys.argv[1:])
	if len(args) < 2 or len(args) > 3:
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
		return 1


	number_of_runs = options.niruns 
	number_of_subruns = options.nruns 



	mpi_init(0, [])

	mpi_comm = MPI_COMM_WORLD

	log = Logger(BaseLogger_Files())

	mpi_rank = mpi_comm_rank(MPI_COMM_WORLD)
	mpi_size = mpi_comm_size(MPI_COMM_WORLD)	# Total number of processes, passed by --np option.


	if mpi_rank == 0:
		all_projs = EMData.read_images(args[0])
		subse_ = range(len(all_projs))
		if mpi_size > len(all_projs):
			ERROR('Number of processes supplied by --np needs to be less than or equal to %d (total number of images) ' % len(all_projs), 'sxviper', 1)
			mpi_finalize()
			return
	else:
		all_projs = None
		subset = None

	mpi_subcomm = wrap_mpi_split(mpi_comm, number_of_runs)

	for runs_iter in range(0, number_of_runs):
		"""
		Iteration over the independent runs of viper function
		""" 
		
		"""
		Generate new directories for processing data
		"""
		postfix_out_dir_name = time.strftime("_%Y_%m_%d__%H_%M_%S__", time.localtime())
		outdir = args[1] + postfix_out_dir_name  + '%03d' % (runs_iter)
		os.mkdir(outdir)

		if mpi_rank == 0:
			if mpi_size % number_of_subruns != 0:
				ERROR('Number of processes needs to be a multiple of total number of runs. Total runs by default are 3, you can change it by specifying --nruns option.', 'sxviper', 1)
				mpi_finalize()
				return

		import global_def
		global_def.LOGFILE =  os.path.join(outdir, global_def.LOGFILE)

		mpi_barrier(MPI_COMM_WORLD)

		if outdir[-1] != "/":
			outdir += "/"
		log.prefix = outdir
		
		if len(args) > 2:
			ref_vol = get_im(args[2])
		else:
			ref_vol = None

		options.user_func = user_functions.factory[options.function]

		#out_params, out_vol, out_peaks = multi_shc(all_projs, subset, number_of_subruns, options, mpi_comm=MPI_COMM_WORLD, log=log, ref_vol=ref_vol)
		out_params, out_vol, out_peaks = multi_shc(all_projs, subset, number_of_subruns, options, mpi_comm=mpi_subcomm,D log=log, ref_vol=ref_vol)

	mpi_finalize()

if __name__=="__main__":
	main()

