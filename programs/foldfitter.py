#!/bin/env python

from EMAN2 import *
from optparse import OptionParser

def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " options inputfile outputfile"

	parser = OptionParser(usage)

	parser.add_option("--apix", type="float", help="the Angstrom/pixel for S scaling")
	parser.add_option("--average", action="store_true", help="Averages all input images (without alignment) and writes a single (normalized) output image")
	parser.add_option("--calcsf", type="string", nargs=2, help="calculate a radial structure factor for the image and write it to the output file, must specify apix. divide into <n> angular bins")    
	parser.add_option("--clip", type="float", nargs=2, action="append", help="Define the output image size")
	
	(options, args) = parser.parse_args()

	
	
if __name__ == "__main__":
    main()
