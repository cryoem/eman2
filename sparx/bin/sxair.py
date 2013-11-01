#!/usr/bin/env python


# returns: subcommunicator, subcommunicator id (from 0..number_of_subcomm-1) and list of roots of subcommunicators (against given mpi_comm)
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
	from mpi import mpi_init, mpi_finalize, MPI_COMM_WORLD, mpi_comm_rank, mpi_comm_size, mpi_comm_split
	from logger import Logger, BaseLogger_ManyFiles
	import sys
	import user_functions
	from optparse import OptionParser
	from global_def import SPARXVERSION
	
	progname = "xxx" #os.path.basename(sys.argv)
	usage = progname + " projections  minimal_subset_size  target_threshold  output_prefix  --ir=inner_radius --ou=outer_radius --rs=ring_step --xr=x_range --yr=y_range  --ts=translational_search_step  --delta=angular_step --an=angular_neighborhood --deltapsi=Delta_psi --startpsi=Start_psi --center=center_type --maxit=max_iter --stoprnct=percentage_to_stop --CTF --snr=SNR  --ref_a=S --sym=c1 --function=user_function --Fourvar=Fourier_variance --debug --MPI"
	parser = OptionParser(usage,version=SPARXVERSION)
#	parser.add_option("--ir",       type= "int",   default= 1,                  help="inner radius for rotational correlation > 0 (set to 1)")
	parser.add_option("--ou",       type= "int",   default= -1,                 help="outer radius for rotational correlation < int(nx/2)-1 (set to the radius of the particle)")
#	parser.add_option("--rs",       type= "int",   default= 1,                  help="step between rings in rotational correlation >0  (set to 1)" ) 
	parser.add_option("--xr",       type="string", default= "0",        help="range for translation search in x direction, search is +/xr")
#	parser.add_option("--yr",       type="string", default= "-1",               help="range for translation search in y direction, search is +/yr (default = same as xr)")
	parser.add_option("--ts",       type="string", default= "1",   help="step size of the translation search in both directions, search is -xr, -xr+ts, 0, xr-ts, xr, can be fractional")
	parser.add_option("--delta",    type="string", default= "2",       help="angular step of reference projections")
#	parser.add_option("--an",       type="string", default= "-1",               help="angular neighborhood for local searches (phi and theta)")
#	parser.add_option("--apsi",     type="string", default= "-1",               help="angular neighborhood for local searches (psi)")
#	parser.add_option("--deltapsi", type="string", default= "-1",               help="Delta psi for coarse search")
#	parser.add_option("--startpsi", type="string", default= "-1",               help="Start psi for coarse search")
	parser.add_option("--center",   type="float",  default= -1,                 help="-1: average shift method; 0: no centering; 1: center of gravity (default=-1)")
	parser.add_option("--maxit",    type="float",  default= 50,                  help="maximum number of iterations performed for each angular step (set to 100) ")
#	parser.add_option("--stoprnct", type="float",  default=0.0,                 help="Minimum percentage of particles that change orientation to stop the program")   
	parser.add_option("--CTF",      action="store_true", default=False,         help="Consider CTF correction during the alignment ")
	parser.add_option("--snr",      type="float",  default= 1.0,                help="Signal-to-Noise Ratio of the data")	
#	parser.add_option("--ref_a",    type="string", default= "S",                help="method for generating the quasi-uniformly distributed projection directions (default S)")
	parser.add_option("--symmetry", type="string", default= "c1",               help="symmetry of the refined structure")
	parser.add_option("--function", type="string", default="[./,paus_302,ali3d_reference3]",         help="name of the reference preparation function ([./,paus_302,ali3d_reference3])")
#	parser.add_option("--MPI",      action="store_true", default=False,         help="whether to use MPI version")
#	parser.add_option("--Fourvar",  action="store_true", default=False,         help="compute Fourier variance")
	parser.add_option("--npad",     type="int",    default= 2,                  help="padding size for 3D reconstruction (default=2)")
#	parser.add_option("--debug",    action="store_true", default=False,         help="debug")
#	parser.add_option("--shc",      action="store_true", default=False,         help="use SHC algorithm")
#	parser.add_option("--nh2",      action="store_true", default=False,         help="new - SHC2")
#	parser.add_option("--ns",       action="store_true", default=False,         help="new - saturn")
#	parser.add_option("--ns2",      action="store_true", default=False,         help="new - saturn2")
#	parser.add_option("--chunk",    type="float",  default= 0.2,                help="percentage of data used for alignment")
#	parser.add_option("--rantest",  action="store_true", default=False,         help="rantest")
#	parser.add_option("--searchpsi",action="store_true", default= False,        help="psi refinement")
	(options, args) = parser.parse_args(sys.argv[1:])
	
	if len(args) < 4:
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
		return 1
	
	
	mpi_init(0, [])
	
	mpi_size = mpi_comm_size(MPI_COMM_WORLD)
	
	if mpi_size < 12:
		print "Number of processes can't be smaller than 12"
		mpi_finalize()
		return
	
	prefix = sys.argv[4]
	
	log=Logger(prefix + "_full", BaseLogger_ManyFiles("log_"))
	
	projs = EMData.read_images(sys.argv[1])
	minimal_subset_size = int(sys.argv[2])
	target_threshold = float(sys.argv[3])
	
	me = wrap_mpi_split(MPI_COMM_WORLD, 4)

	if options.ou < 1:
		options.ou = (projs[0].get_xsize() - 1) / 2

# 	options.symmetry = "c1"
# 	options.ou=30
# 	options.xr="0"
# 	options.yr="0"
# 	options.ts="1"
# 	options.delta="2"
# 	options.center=0
# 	options.maxit=100
# 	options.CTF=False
# 	options.snr=1.0
# 	options.function="[./,paus_302,ali3d_reference3]"
# 	options.npad=2
	options.user_func = user_functions.factory[options.function]

	new_subset, new_threshold = full_proc(projs, minimal_subset_size, target_threshold, options, number_of_runs=4, number_of_winners=3, mpi_env=me, log=log, prefix=prefix)

	log.add("Output threshold =", new_threshold)
	log.add("Output subset =", new_subset)

	write_text_file(new_subset, prefix + "_out.txt")

	mpi_finalize()


if __name__=="__main__":
	main()
