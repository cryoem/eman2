#!/bin/env python
# e2boxer.py  07/27/3004  Steven Ludtke
# This program is used to box out particles from micrographs/CCD frames

from EMAN2 import *
from optparse import OptionParser
from math import *
import os
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage: %prog [options] <image>
	
Automatic and manual particle selection. This version is specifically aimed at square boxes
for single particle analysis."""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--gui",action="store_true",help="Start the GUI for interactive boxing",default=false)
	parser.add_option("--box","-B",type="int",help="Box size in pixels",default=-1)
	parser.add_option("--ptclsize","-P",type="int",help="Approximate size (diameter) of the particle in pixels. Not required if reference particles are provided.",default=-1)
	parser.add_option("--refptcl","-R",type="string",help="A stack of reference images. Must have the same scale as the image being boxed.",default=None)
	parser.add_option("--auto","-A",type="string",action="append",help="Autobox using specified method: circle")
			
	(options, args) = parser.parse_args()
	if len(args)<1 : parser.error("Input image required")

	image=EMData()
	image.read_image(args[0])
	
	if "circle" in options.auto:
		
	
	
	
	
if __name__ == "__main__":
	main()
