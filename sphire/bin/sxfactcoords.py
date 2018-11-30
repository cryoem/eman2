#!/usr/bin/env python
from __future__ import print_function

pass#IMPORTIMPORTIMPORT import global_def
pass#IMPORTIMPORTIMPORT from global_def import *
pass#IMPORTIMPORTIMPORT from optparse import OptionParser
pass#IMPORTIMPORTIMPORT from EMAN2_cppwrap import *

pass#IMPORTIMPORTIMPORT import os
pass#IMPORTIMPORTIMPORT import sys

      
pass#IMPORTIMPORTIMPORT import EMAN2_cppwrap
pass#IMPORTIMPORTIMPORT import applications
pass#IMPORTIMPORTIMPORT import global_def
pass#IMPORTIMPORTIMPORT import mpi
pass#IMPORTIMPORTIMPORT import optparse
pass#IMPORTIMPORTIMPORT import os
pass#IMPORTIMPORTIMPORT import sys
pass#IMPORTIMPORTIMPORT import utilities
import applications
import global_def
import mpi
import optparse
import os
import sys
import utilities
def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " prj_stack .. average eigvol output_factcoords --rad=radius --neigvol=number_of_eigvol  --CTF"
	parser = optparse.OptionParser(usage, version=global_def.SPARXVERSION)
	parser.add_option("--rad",       type="int",    default=-1,     help="radius of mask")
	parser.add_option("--neigvol",   type="int",    default=-1,     help="number of eigvenvectors to use (default all)")
	parser.add_option("--fl",        type="float",  default=0.0,    help="cut-off frequency of hyperbolic tangent low-pass Fourier filter")
	parser.add_option("--aa",        type="float",  default=0.0,    help="fall-off of hyperbolic tangent low-pass Fourier filter")
	parser.add_option("--CTF",       action="store_true", default=False,  help="Use CTF")
	parser.add_option("--MPI",       action="store_true",           help="use MPI")

	(options, args) = parser.parse_args()

	if( len(args) < 4 ):
		print("usage: " + usage)
		print("Please run '" + progname + " -h' for details")
	else:
		stacks = args[0:-3]
		avgvol = args[-3]
		eigvol = args[-2]
		output = args[-1]
		
		if options.rad < 0:
			print("Error: mask radius is not given")
			sys.exit(-1)
		if global_def.CACHE_DISABLE:
			pass#IMPORTIMPORTIMPORT from utilities import disable_bdb_cache
			utilities.disable_bdb_cache()
		if options.MPI:
			pass#IMPORTIMPORTIMPORT from mpi import mpi_init
			sys.argv = mpi.mpi_init(len(sys.argv), sys.argv)

		pass#IMPORTIMPORTIMPORT from utilities import get_im
		global_def.BATCH = True
		if( utilities.get_im( stacks[0]).get_zsize() == 1 and utilities.get_im( eigvol).get_zsize() > 1):
			pass#IMPORTIMPORTIMPORT from applications import factcoords_prj
			applications.factcoords_prj(stacks, avgvol, eigvol, output, options.rad, options.neigvol, options.fl, options.aa, options.CTF, options.MPI)
		else:
			pass#IMPORTIMPORTIMPORT from applications import factcoords_vol
			applications.factcoords_vol(stacks, avgvol, eigvol, output, options.rad, options.neigvol, options.fl, options.aa, options.MPI)
		global_def.BATCH = False
		

if __name__ == "__main__":
	main()
