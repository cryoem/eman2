#!/usr/bin/env python

import global_def
from global_def import *
from optparse import OptionParser

import os
import sys

      
def main():
    progname = os.path.basename(sys.argv[0])
    usage = progname + " input_stack start end output_stack --rad=mask_radius"
    parser = OptionParser(usage, version=SPARXVERSION)
    parser.add_option("--rad", type="int",  help="radius of mask")
    parser.add_option("--verbose", type="int", default=0,  help="verbose level (0|1)")


    (options, args) = parser.parse_args()

    

    if len(args) !=  4:
        print "usage: " + usage
	print "Please run '" + progname + " -h' for details"
    else:
        from string import atoi
        input_stack  = args[0]
        imgstart     = atoi( args[1] )
        imgend       = atoi( args[2] ) +1
        output_stack = args[3]

        if options.rad is None:
            print "Error: mask radius is not given"
            sys.exit(-1)

        from applications import varimax
	global_def.BATCH = True
        varimax(input_stack, range(imgstart, imgend), output_stack, options.rad, options.verbose)
	global_def.BATCH = False

if __name__ == "__main__":
	main()     


