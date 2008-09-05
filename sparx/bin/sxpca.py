#!/usr/bin/env python

import global_def
from global_def import *
from optparse import OptionParser
from EMAN2_cppwrap import *

import os
import sys

      
def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " input_stack output_stack [start end step] --nfile=number_of_input --subavg=average_image --rad=mask_radius --nvec=number_of_eigenvectors --mask=maskfile"
	parser = OptionParser(usage, version=SPARXVERSION)
	parser.add_option("--nfile",  type="int", default=1, help="number of files")
	parser.add_option("--subavg", type="string", default="", help="subtract average")
	parser.add_option("--rad",    type="int",  default=-1, help="radius of mask")
	parser.add_option("--nvec",   type="int", help="number of eigenvectors")
	parser.add_option("--type",   type="string", default="out_of_core", help="out_of_core calculations")
	parser.add_option("--mask",   type="string", default="", help="mask file" )

	(options, args) = parser.parse_args()

	if len(args) != 2  and len(args) != 5:
		print "usage: " + usage
		print "Please run '" + progname + " -h' for details"
	else:
		input_stack  = args[0]
		output_stack = args[1]

		if options.nfile==1:
			if len(args)==2:
				imgstart = 0
				imgend = EMUtil.get_image_count( input_stack )
				imgstep = 1
				imglist = range(imgstart,imgend,imgstep)
		    	else:
				from string import atoi
				imgstart = atoi( args[2] )
				imgend   = atoi( args[3] )
				imgstep  = atoi( args[4] )
				imglist = range(imgstart, imgend, imgstep)
		else:
			imglist = None

		if options.rad is None:
			print "Error: mask radius is not given"
			sys.exit(-1)

		if options.nvec is None:
			print "Error: number of components is not given"
			sys.exit(-2) 

		from applications import pca
		global_def.BATCH = True
		pca(input_stack, output_stack, imglist, options.nfile, options.subavg, options.rad, options.nvec, options.type, options.mask)
		global_def.BATCH = False

if __name__ == "__main__":
	main()
