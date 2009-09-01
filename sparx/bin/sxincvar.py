#!/usr/bin/env python
import sys
import os
from global_def import SPARXVERSION
import global_def

def main():
	from mpi import mpi_init
	sys.argv = mpi_init( len(sys.argv), sys.argv )


	arglist = []
	for arg in sys.argv:
		arglist.append( arg )


	from   optparse       import OptionParser
	progname = os.path.basename(arglist[0])
	usage = progname + " prefix nfile nprj outputfile  --fl=flit_low_value --fh=filt_high_value --radccc=radius_ccc --writelp --writestack"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--fl",        type="float",  default=0.2,    help="first parameter for low pass filter")
	parser.add_option("--fh",        type="float",  default=0.4,    help="second parameter for low pass filter")
	parser.add_option("--radccc",    type="int",    default=-1,     help="radius for ccc caclualtion")
	parser.add_option("--writelp",   action="store_true", default=False, help="if write the low pass filtered volume to disk (default is False)" )
	parser.add_option("--writestack", action="store_true", default=False, help="if write the stack contain all variance map" )
	parser.add_option("--MPI", action="store_true", default=False, help="if use MPI version" )
	(options, args) = parser.parse_args(arglist[1:])

	if len(args)!=4 :
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
	else:
		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()
		prefix = args[0]
		nfile  = int(args[1])
		nprj   = int(args[2])
		output = args[3]
		if options.MPI:
			print 'mpi version of incvar'
			from applications import incvar_mpi
			incvar_mpi( prefix, nfile, nprj, output, options.fl, options.fh, options.radccc, options.writelp, options.writestack)
		else:
			from applications import incvar
			incvar( prefix, nfile, nprj, output, options.fl, options.fh, options.radccc, options.writelp, options.writestack)


if __name__ == "__main__":
	main()
