#!/usr/bin/env python

import global_def
from global_def import *
from optparse import OptionParser

import os
import sys

      
def main():
    progname = os.path.basename(sys.argv[0])
    usage = progname + " input_stack start end output_stack  <mask> --rad=mask_radius"
    parser = OptionParser(usage, version=SPARXVERSION)
    parser.add_option("--rad",     type="int", default=-1, help="radius of mask")
    parser.add_option("--verbose", type="int", default=0,  help="verbose level (0|1)")


    (options, args) = parser.parse_args()

    

    if len(args) < 4:
        print "usage: " + usage
	print "Please run '" + progname + " -h' for details"
    else:
        from string import atoi
        input_stack  = args[0]
        imgstart     = atoi( args[1] )
        imgend       = atoi( args[2] ) +1
        output_stack = args[3]
	if(len(args) == 5):  mask = args[4]
	else:               mask = None

	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()
        from applications import varimax
	global_def.BATCH = True
        varimax(input_stack, range(imgstart, imgend), output_stack, mask, options.rad, options.verbose)
	global_def.BATCH = False

if __name__ == "__main__":
	main()     


