#!/usr/bin/env python
import sys
import os

import global_def
from global_def import *

def main():
	from mpi import mpi_init
	sys.argv = mpi_init( len(sys.argv), sys.argv )

	arglist = []
	for arg in sys.argv:
		arglist.append( arg )


	from   optparse       import OptionParser
	progname = os.path.basename(arglist[0])
	usage = progname + " filelist outdir  --fl=flit_low_value --fh=filt_high_value --radccc=radius_ccc --writelp --writestack --MPI"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--fl",        type="float",  default=0.2,    help="cut-off frequency of hyperbolic tangent low-pass Fourier filter")
	parser.add_option("--fh",        type="float",  default=0.4,    help="fall-off of hyperbolic tangent low-pass Fourier filter")
	parser.add_option("--radccc",    type="int",    default=-1,     help="radius for ccc caclualtion")
	parser.add_option("--writelp",   action="store_true", default=False, help="write the low pass filtered volume to disk (default is False)" )
	parser.add_option("--writestack", action="store_true", default=False, help="write the stack contain all variance map" )
	parser.add_option("--MPI", action="store_true", default=False, help="use MPI version" )
	(options, args) = parser.parse_args(arglist[1:])

	if len(args)<2 :
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
	else:
		files = args[0:-1]
		outdir = args[-1]

		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()

		if options.MPI:
			from applications import defvar_mpi
			defvar_mpi( files, outdir, options.fl, options.fh, options.radccc, options.writelp, options.writestack)
			from mpi import mpi_finalize
			mpi_finalize()
		else:
			from applications import defvar
			defvar( files, outdir, options.fl, options.fh, options.radccc, options.writelp, options.writestack)


if __name__ == "__main__":
	main()
