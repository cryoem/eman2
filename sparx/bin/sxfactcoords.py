#!/usr/bin/env python

import global_def
from global_def import *
from optparse import OptionParser
from EMAN2_cppwrap import *

import os
import sys

      
def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " prj_stack .. average eigvol output_factcoords --rad=radius --neigvol=number_of_eigvol"
	parser = OptionParser(usage, version=SPARXVERSION)
	parser.add_option("--rad",       type="int",    default=-1, help="radius of mask")
	parser.add_option("--neigvol",   type="int",    default=-1, help="number of eigvenvectors to use (default all)")
	parser.add_option("--fl",        type="float",  default=0.0,    help="cut-off frequency of hyperbolic tangent low-pass Fourier filter")
	parser.add_option("--aa",        type="float",  default=0.0,    help="fall-off of hyperbolic tangent low-pass Fourier filter")
	parser.add_option("--MPI",       action="store_true", help="flag for MPI version")

	(options, args) = parser.parse_args()

	if( len(args) < 4 ):
		print "usage: " + usage
		print "Please run '" + progname + " -h' for details"
	else:
		stacks = args[0:-3]
		avgvol = args[-3]
		eigvol = args[-2]
		output = args[-1]
		
		if options.rad < 0:
			print "Error: mask radius is not given"
			sys.exit(-1)
		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()

		from utilities import get_im
		global_def.BATCH = True
		if( get_im( stacks[0]).get_zsize() == 1):
			from applications import factcoords_prj
			factcoords_prj(stacks, avgvol, eigvol, output, options.rad, options.neigvol, options.fl, options.aa, options.MPI)
		else:
			from applications import factcoords_vol
			factcoords_vol(stacks, avgvol, eigvol, output, options.rad, options.neigvol, options.fl, options.aa, options.MPI)
		global_def.BATCH = False
		

if __name__ == "__main__":
	main()
