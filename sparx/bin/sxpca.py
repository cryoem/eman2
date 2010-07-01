#!/usr/bin/env python

import global_def
from global_def import *
from optparse import OptionParser
from EMAN2_cppwrap import *

import os
import sys

      
def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " input_stack1 ... output_stack --subavg=average_image --rad=mask_radius --nvec=number_of_eigenvectors --mask=maskfile"
	parser = OptionParser(usage, version=SPARXVERSION)
	parser.add_option("--subavg",  type="string",       default="",    help="subtract average")
	parser.add_option("--rad",     type="int",          default=-1,    help="radius of mask")
	parser.add_option("--nvec",    type="int",          default=1,     help="number of eigenvectors")
	parser.add_option("--mask",    type="string",       default="",    help="mask file" )
	parser.add_option("--sdir",    type="string",       default=".",   help="scratch directory")
	parser.add_option("--usebuf",  action="store_true", default=False, help="use existing buffer")
	parser.add_option("--MPI",     action="store_true", default=False, help="run mpi version" )
	parser.add_option("--shuffle", action="store_true", default=False, help="use shuffle")

	(options, args) = parser.parse_args()

	input_stacks = args[0:-1]
	output_stack = args[-1]

	if options.nvec is None:
		print "Error: number of components is not given"
		sys.exit(-2) 

	if options.MPI:
		from mpi import mpi_init
		sys.argv = mpi_init( len(sys.argv), sys.argv )

	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()
	from applications import pca
	global_def.BATCH = True
	pca(input_stacks, output_stack, options.subavg, options.rad, options.sdir, options.nvec, options.shuffle, not(options.usebuf), options.mask, options.MPI)
	global_def.BATCH = False
        if options.MPI:
		from mpi import mpi_finalize
		mpi_finalize()


if __name__ == "__main__":
	main()
