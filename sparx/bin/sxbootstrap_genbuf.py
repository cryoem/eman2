#!/usr/bin/env python

import global_def
from global_def import *
from optparse import OptionParser

import os
import sys
import string


def main():

	arglist = []
	for arg in sys.argv:
		arglist.append( arg )

	progname = os.path.basename( arglist[0] )
	usage = progname + " prj_stack outdir --verbose=(0|1) --MPI"
	parser = OptionParser(usage, version=SPARXVERSION)
	parser.add_option("--verbose", type="int", default=0, help="verbose level: 0 no verbose, 1 verbose" )
	parser.add_option("--MPI", action="store_true", default=False,     help="  whether using MPI version ")
	parser.add_option("--CTF", action="store_true", default=False,     help="  whether consider CTF" )
	(options,args) = parser.parse_args(arglist[1:])

	if len(args) != 2 :
		print usage
		sys.exit(-1)

	proj_stack = args[0]
	outdir = args[1]

	from applications import bootstrap_genbuf
	global_def.BATCH = True
	if options.MPI: 
		from mpi import mpi_init
		sys.argv = mpi_init(len(sys.argv), sys.argv)
	bootstrap_genbuf(proj_stack, outdir, options.verbose, options.CTF, options.MPI) 
	global_def.BATCH = False

if __name__ == "__main__":
	main()


