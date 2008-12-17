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
	parser.add_option("--subavg", type="string", default="", help="subtract average")
	parser.add_option("--rad",    type="int",  default=-1, help="radius of mask")
	parser.add_option("--nvec",   type="int", help="number of eigenvectors")
	parser.add_option("--type",   type="string", default="out_of_core", help="out_of_core calculations")
	parser.add_option("--mask",   type="string", default="", help="mask file" )

	(options, args) = parser.parse_args()

	input_stacks = args[0:-1]
	output_stack = args[1-]

	if options.nvec is None:
		print "Error: number of components is not given"
		sys.exit(-2) 

	from applications import pca
	global_def.BATCH = True
	pca(input_stacks, output_stack, options.subavg, options.rad, options.nvec, options.type, options.mask)
	global_def.BATCH = False

if __name__ == "__main__":
	main()
