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
	usage = progname + " stack  output_directory  [initial_volume]  --ir=inner_radius --ou=outer_radius --rs=ring_step --xr=x_range --yr=y_range  --ts=translational_search_step  --delta=angular_step --an=angular_neighborhood  --center=center_type --maxit=max_iter --CTF --snr=SNR  --ref_a=S --sym=c1 --function=user_function"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--ir",       type= "int",   default= 1,                  help="inner radius for rotational correlation > 0 (set to 1)")
	parser.add_option("--ou",       type= "int",   default= -1,                 help="outer radius for rotational correlation < int(nx/2)-1 (set to the radius of the particle)")
	parser.add_option("--rs",       type= "int",   default= 1,                  help="step between rings in rotational correlation >0  (set to 1)" ) 
	parser.add_option("--xr",       type="string", default= "0",                help="range for translation search in x direction, search is +/xr")
	parser.add_option("--yr",       type="string", default= "-1",               help="range for translation search in y direction, search is +/yr (default = same as xr)")
	parser.add_option("--ts",       type="string", default= "1",                help="step size of the translation search in both directions, search is -xr, -xr+ts, 0, xr-ts, xr, can be fractional")
	parser.add_option("--delta",    type="string", default= "2",                help="angular step of reference projections")
	parser.add_option("--an",       type="string", default= "-1",               help="angular neighborhood for local searches (phi and theta)")
	parser.add_option("--center",   type="float",  default= -1,                 help="-1: average shift method; 0: no centering; 1: center of gravity (default=-1)")
	parser.add_option("--maxit",    type="float",  default= 100,                help="maximum number of iterations performed for each angular step (set to 100) ")
	parser.add_option("--CTF",      action="store_true", default=False,         help="Consider CTF correction during the alignment ")
	parser.add_option("--snr",      type="float",  default= 1.0,                help="Signal-to-Noise Ratio of the data")	
	parser.add_option("--ref_a",    type="string", default= "S",                help="method for generating the quasi-uniformly distributed projection directions (default S)")
	parser.add_option("--sym",      type="string", default= "c1",               help="symmetry of the refined structure")
	parser.add_option("--function", type="string", default="ref_ali3d",         help="name of the reference preparation function (ref_ali3d by default)")
	parser.add_option("--npad",     type="int",    default= 2,                  help="padding size for 3D reconstruction (default=2)")
	#parser.add_option("--MPI",      action="store_true", default=True,          help="whether to use MPI version - this is always set to True")
	(options, args) = parser.parse_args(sys.argv[1:])
	if len(args) < 2 or len(args) > 3:
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
		return 1

	mpi_init(0, [])

	log = Logger(BaseLogger_Files())

	runs_count = 3
	mpi_rank = mpi_comm_rank(MPI_COMM_WORLD)

	if mpi_rank == 0:
		all_projs = EMData.read_images(args[0])
		subset = range(len(all_projs))
	else:
		all_projs = None
		subset = None

	outdir = args[1]
	if mpi_rank == 0:
		if os.path.exists(outdir):
			ERROR('Output directory exists, please change the name and restart the program', "sxmulti_shc", 1)
			mpi_finalize()
			return
		os.mkdir(outdir)
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

	out_params, out_vol, out_peaks = multi_shc(all_projs, subset, runs_count, options, mpi_comm=MPI_COMM_WORLD, log=log, ref_vol=ref_vol)

	mpi_finalize()

if __name__=="__main__":
	main()

