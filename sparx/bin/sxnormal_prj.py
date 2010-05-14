#!/usr/bin/env python


import os
import global_def
from   global_def import *
from   optparse import OptionParser
import sys


def main():
        arglist = []
        for arg in sys.argv:
        	arglist.append( arg )
	progname = os.path.basename(sys.argv[0])
	usage = progname + " stack outdir <refvol> --apply_weights=weights  --r=radius --niter=number_of_iteration --snr=Signal-to-Noise Ratio --sym=symmetry --CTF --MPI"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--apply_weights",     type="string",   default=None,  help="text file with projections' weights in the first column" )
	parser.add_option("--r",       type="int",      default=-1,    help="radius < int(nx/2)-1 (set to int(nx/2)-1)")
	parser.add_option("--niter",   type="int",	default=1,     help="number of iteration" )
	parser.add_option("--snr",     type="float",    default=1.0,   help="Signal-to-Noise ratio" )
	parser.add_option("--sym",     type="string",   default="c1",  help="symmetry" )
	parser.add_option("--CTF", action="store_true", default=False, help="consider CTF")
	parser.add_option("--MPI", action="store_true", default=False, help="use MPI version")
	parser.add_option("--verbose", type="int",      default=0,     help="verbose level: 0 no, 1 yes" )

	(options, args) = parser.parse_args( arglist[1:] )

	if( len(args)!=2 and len(args)!=3 ):
		print "usage: " + usage
		print "Please run'" + progname + " -h' for detailed options"
		from sys import exit
		exit(-1)

	if len(args)==3 :
		from utilities import get_im
		refvol = get_im( args[2] )
	else:
		refvol = None

	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()
	from applications import normal_prj
	if options.MPI:
		from mpi import mpi_init
		sys.argv = mpi_init(len(sys.argv), sys.argv)
	normal_prj( args[0], args[1], refvol, options.apply_weights, options.r, options.niter, options.snr, options.sym, options.verbose, options.CTF, options.MPI )


if __name__ == "__main__":
	main()
