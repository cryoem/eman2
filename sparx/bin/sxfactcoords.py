#!/usr/bin/env python

import global_def
from global_def import *
from optparse import OptionParser
from EMAN2_cppwrap import *

import os
import sys

      
def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " prj_stack average eigvol output_factcoords --rad=radius --neigvol=number_of_eigvol --of=output_format"
	parser = OptionParser(usage, version=SPARXVERSION)
	parser.add_option("--rad", type="int",  default=-1, help="radius of mask")
	parser.add_option("--neigvol", type="int", default=-1, help="number of eigvenvectors to use (default all)")
	parser.add_option("--of", type="string", default="hdf", help="output format: hdf or txt (default is hdf)")

	(options, args) = parser.parse_args()

	if( len(args) > 4 or len(args) < 3):
		print "usage: " + usage
		print "Please run '" + progname + " -h' for details"
	else:
		prj_stack  = args[0]
		if len(args) == 4:
			avgvol = args[1]
			eigvol = args[2]
			output = args[3]
		else:
			avgvol = None
			eigvol = args[1]
			output = args[2]

		if options.rad < 0:
			print "Error: mask radius is not given"
			sys.exit(-1)

		global_def.BATCH = True
		v = EMData()
		v.read_image(prj_stack, True)
		if(v.get_zsize() > 1):
			from applications import factcoords3D
			factcoords3D(prj_stack, avgvol, eigvol, output, options.rad, options.neigvol, options.of)
		else:
			from applications import factcoords2D
			factcoords2D(prj_stack, avgvol, eigvol, output, options.rad, options.neigvol, options.of)
		global_def.BATCH = False

if __name__ == "__main__":
	main()
