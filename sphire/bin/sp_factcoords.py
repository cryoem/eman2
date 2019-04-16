#!/usr/bin/env python
from __future__ import print_function

import sp_global_def
from sp_global_def import sxprint, ERROR

from sp_global_def import *
from optparse import OptionParser
from EMAN2_cppwrap import *

import os
import sys

import mpi

mpi.mpi_init( 0, [] )
      
def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " prj_stack .. average eigvol output_factcoords --rad=radius --neigvol=number_of_eigvol  --CTF"
	parser = OptionParser(usage, version=SPARXVERSION)
	parser.add_option("--rad",       type="int",    default=-1,     help="radius of mask")
	parser.add_option("--neigvol",   type="int",    default=-1,     help="number of eigvenvectors to use (default all)")
	parser.add_option("--fl",        type="float",  default=0.0,    help="cut-off frequency of hyperbolic tangent low-pass Fourier filter")
	parser.add_option("--aa",        type="float",  default=0.0,    help="fall-off of hyperbolic tangent low-pass Fourier filter")
	parser.add_option("--CTF",       action="store_true", default=False,  help="Use CTF")
	parser.add_option("--MPI",       action="store_true",           help="use MPI")

	(options, args) = parser.parse_args()

	if( len(args) < 4 ):
		sxprint( "Usage: " + usage )
		sxprint( "Please run \'" + progname + " -h\' for detailed options" )
		ERROR( "Invalid number of parameters used. Please see usage information above." )
		return

	else:
		stacks = args[0:-3]
		avgvol = args[-3]
		eigvol = args[-2]
		output = args[-1]
		
		if options.rad < 0:
			ERROR( "Mask radius is not given" )
			return

		if sp_global_def.CACHE_DISABLE:
			from sp_utilities import disable_bdb_cache
			disable_bdb_cache()

		from sp_utilities import get_im
		sp_global_def.BATCH = True
		
		if( get_im( stacks[0]).get_zsize() == 1 and get_im( eigvol).get_zsize() > 1):
			from sp_applications import factcoords_prj
			factcoords_prj(stacks, avgvol, eigvol, output, options.rad, options.neigvol, options.fl, options.aa, options.CTF, options.MPI)
		else:
			from sp_applications import factcoords_vol
			factcoords_vol(stacks, avgvol, eigvol, output, options.rad, options.neigvol, options.fl, options.aa, options.MPI)
		sp_global_def.BATCH = False
		

if __name__ == "__main__":
	sp_global_def.print_timestamp( "Start" )
	sp_global_def.write_command()
	main()
	sp_global_def.print_timestamp( "Finish" )
	mpi.mpi_finalize()