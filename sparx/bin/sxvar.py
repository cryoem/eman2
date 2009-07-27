#!/usr/bin/env python
import sys
import os
import global_def
from global_def import *

def main():
	from   optparse       import OptionParser
	progname = os.path.basename(sys.argv[0])
	usage = progname + " filelist outdir  --fl=filter_low_value --aa=filter_fall_off --radccc=radius_ccc --overwrite --filtered --repair=repairfile --pca --pcamask --pcanvec --MPI"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--fl",         type="float",  default=0.0,    help="cut-off frequency of hyperbolic tangent low-pass Fourier filter")
	parser.add_option("--aa",         type="float",  default=0.0,    help="fall-off of hyperbolic tangent low-pass Fourier filter")
	parser.add_option("--radccc",     type="int",    default=-1,     help="radius for ccc calculation")
	parser.add_option("--overwrite",  action="store_true", default=False, help="write repaired bootstrap volumes to original files (default is False)" )
	parser.add_option("--filtered",   action="store_true", default=False, help="write the stack containing all low-pass filtered volumes to disk (default is False)" )
	parser.add_option("--MPI",        action="store_true", default=False, help="use MPI version" )
	parser.add_option("--repair",     type="string", help="repair file for original bootstrap volumes")
	parser.add_option("--pca",        action="store_true", default=False, help="run pca" )
	parser.add_option("--pcamask",    type="string", help="mask for pca" )
	parser.add_option("--pcanvec",    type="int", help="number of eigvectors computed in PCA")
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


			global_def.BATCH = True
			var_mpi( files, outdir, options.fl, options.aa, options.radccc, options.overwrite, options.filtered, options.repair, options.pca, options.pcamask, options.pcanvec)
			global_def.BATCH = False
		else:
			from applications import defvar
			global_def.BATCH = True
			defvar(  files, outdir, options.fl, options.aa, options.radccc, options.overwrite, options.filtered, options.repair, options.pca, options.pcamask, options.pcanvec)
			global_def.BATCH = False


if __name__ == "__main__":
	main()
