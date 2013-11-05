#!/usr/bin/env python



def wrap_mpi_split(mpi_comm, number_of_subcomm):
	from mpi import mpi_comm_rank, mpi_comm_size, mpi_comm_split
	from air import mpi_env_type
	
	main_size = mpi_comm_size(mpi_comm)
	if number_of_subcomm > main_size:
		raise RuntimeError("number_of_subcomm > main_size")
	
	me = mpi_env_type()
	me.main_comm = mpi_comm
	me.main_rank = mpi_comm_rank(mpi_comm)
	me.subcomm_id = me.main_rank % number_of_subcomm
	me.sub_rank = me.main_rank / number_of_subcomm
	me.sub_comm = mpi_comm_split(mpi_comm, me.subcomm_id, me.sub_rank)
	me.subcomms_count = number_of_subcomm
	me.subcomms_roots = range(number_of_subcomm)

	return me


def main():
	from EMAN2 import EMData
	from utilities import write_text_file
	from mpi import mpi_init, mpi_finalize, MPI_COMM_WORLD, mpi_comm_rank, mpi_comm_size, mpi_comm_split, mpi_barrier
	from logger import Logger, BaseLogger_Files
	from air import air
	import sys
	import os
	import user_functions
	from optparse import OptionParser
	from global_def import SPARXVERSION
	
	progname = os.path.basename(sys.argv[0])
	usage = progname + " projections  minimal_subset_size  target_threshold  output_directory --ir=inner_radius --ou=outer_radius --rs=ring_step --xr=x_range --yr=y_range  --ts=translational_search_step  --delta=angular_step --an=angular_neighborhood  --center=center_type --maxit=max_iter --CTF --snr=SNR  --ref_a=S --sym=c1 --function=user_function --MPI"
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
	parser.add_option("--maxit",    type="float",  default= 50,                 help="maximum number of iterations performed for each angular step (set to 50) ")
	parser.add_option("--CTF",      action="store_true", default=False,         help="Consider CTF correction during the alignment ")
	parser.add_option("--snr",      type="float",  default= 1.0,                help="Signal-to-Noise Ratio of the data")	
	parser.add_option("--ref_a",    type="string", default= "S",                help="method for generating the quasi-uniformly distributed projection directions (default S)")
	parser.add_option("--sym",      type="string", default= "c1",               help="symmetry of the refined structure")
	parser.add_option("--function", type="string", default="ref_ali3d",         help="name of the reference preparation function (ref_ali3d by default)")
	parser.add_option("--npad",     type="int",    default= 2,                  help="padding size for 3D reconstruction (default=2)")
	parser.add_option("--MPI",      action="store_true", default=True,          help="whether to use MPI version - this is always set to True")
	parser.add_option("--proc_mshc",type="int",    default=3,                   help="number of MPI processes per multiSHC, 3 is minimum (default=3)")
	(options, args) = parser.parse_args(sys.argv[1:])
	
	if len(args) < 4:
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
		return 1
	
	
	mpi_init(0, [])
	
	mpi_size = mpi_comm_size(MPI_COMM_WORLD)
	mpi_rank = mpi_comm_rank(MPI_COMM_WORLD)
	
	proc_per_mshc = int(options.proc_mshc)
	
	if mpi_size < proc_per_mshc:
		print "Number of processes can't be smaller than value given as the parameter --proc_mshc"
		mpi_finalize()
		return
	
	log = Logger(BaseLogger_Files())
	
	projs = EMData.read_images(args[0])
	minimal_subset_size = int(args[1])
	target_threshold = float(args[2])
	outdir = args[3]
	
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
	
	me = wrap_mpi_split(MPI_COMM_WORLD, mpi_size / proc_per_mshc )

	options.user_func = user_functions.factory[options.function]

	new_subset, new_threshold = air(projs, minimal_subset_size, target_threshold, options, number_of_runs=6, number_of_winners=3, mpi_env=me, log=log)

	if mpi_rank == 0:
		log.add("Output threshold =", new_threshold)
		log.add("Output subset: ", len(new_subset), new_subset)
		write_text_file(new_subset, log.prefix + "final_subset.txt")

	mpi_finalize()


if __name__=="__main__":
	main()
