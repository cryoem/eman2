#!/usr/bin/env python
import sys
import os
from global_def import SPARXVERSION

def main():
	from   optparse       import OptionParser
	progname = os.path.basename(sys.argv[0])
	usage = progname + " filelist outdir  --fl=flit_low_value --aa=filt_fall_off --radccc=radius_ccc --writelp --writestack --MPI"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--fl",        type="float",  default=0.0,    help="cut-off frequency of hyperbolic tangent low-pass Fourier filter")
	parser.add_option("--aa",        type="float",  default=0.0,    help="fall-off of hyperbolic tangent low-pass Fourier filter")
	parser.add_option("--radccc",    type="int",    default=-1,     help="radius for ccc caclualtion")
	parser.add_option("--writelp",   action="store_true", default=False, help="write the low pass filtered volume to disk (default is False)" )
	parser.add_option("--writestack", action="store_true", default=False, help="write the stack contain all variance map" )
	parser.add_option("--MPI", action="store_true", default=False, help="use MPI version" )
	parser.add_option("--pca", action="store_true", default=False, help="run pca" )
	parser.add_option("--pcamask", type="string", help="mask for pca" )
	parser.add_option("--pcanvec", type="int", help="number of eigvectors for pca")
	parser.add_option("--method", type="string", default="inc", help="calculation method: def: calculate by definition, inc: use incremental. default is inc")
	(options, args) = parser.parse_args(sys.argv[1:])

	if len(args)<2 :
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
	else:
		files = args[0:-1]
		outdir = args[-1]

		if options.MPI:
			from mpi import mpi_init
			sys.argv = mpi_init( len(sys.argv), sys.argv )


			arglist = []
			for arg in sys.argv:
				arglist.append( arg )
			from applications import var_mpi
			from utilities import init_mpi_bdb
			init_mpi_bdb()


			var_mpi( files, outdir, options.fl, options.aa, options.radccc, options.writelp, options.writestack, options.method, options.pca, options.pcamask, options.pcanvec)
		else:
			from applications import defvar
			defvar( files, outdir, options.fl, options.aa, options.radccc, options.writelp, options.writestack)


if __name__ == "__main__":
	main()
