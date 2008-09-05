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
	usage = progname + " stack outdir <refvol> --r=radius --niter=number_of_iteration --snr=Signal-to-Noise Ratio --sym=symmetry --MPI"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--r",     type="int",        default=-1,    help="  radius < int(nx/2)-1 (set to int(nx/2)-1)")
	parser.add_option("--niter", type="int",	default=1,     help="  number of iteration" )
	parser.add_option("--snr",   type="float",      default=1.0,   help="  Signal-to-Noise ratio" )
	parser.add_option("--sym",   type="string",     default="c1",  help="  Symmetry" )
	parser.add_option("--MPI", action="store_true", default=False,     help="  whether using MPI version ")

	(options, args) = parser.parse_args( arglist[1:] )

	if( len(args)!=2 and len(args)!=3 ):
		print "usage: " + usage
		print "Please run'" + progname + " -h' for detailed options"
		from sys import exit
		exit(-1)

	if len(args)==3 :
		from utilities import getImage
		refvol = getImage( args[2] )
	else:
		refvol = None

	from applications import normal_prj
	if options.MPI:
		from mpi import mpi_init
		sys.argv = mpi_init(len(sys.argv), sys.argv)
	normal_prj( args[0], args[1], refvol, options.r, options.niter, options.snr, options.sym, options.MPI )


if __name__ == "__main__":
	main()
