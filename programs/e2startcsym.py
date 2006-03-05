#!/usr/bin/env python

# $Id$


from EMAN2 import *

from optparse import OptionParser
import os.path
import sys
import pyemtbx.files
import e2sysbest

def main():
	progname = os.path.basename(sys.argv[0])
	
	usage = progname + " options inputfile"
	parser = OptionParser(usage,version=EMANVERSION)

	parser.add_option("--nkeep", metavar="N", type="int", help="Number of particles to use for each view")
	parser.add_option("--out", metavar="outputfile", type="string", help="output filename")
	parser.add_option("--mode", metavar="n", type="int", help="")
	parser.add_option("--sym", metavar="Cn", type="string", help="Symmetry of model")

	parser.add_option("--imask", metavar="rad", type="int", help="Inside mask uesd to exclude inside regions")
	parser.add_option("--ccl", action="store_true", help="")

	parser.add_option("--nosym", action="store_true", help="Skips the initial symmetry search.")
	parser.add_option("--fixrot", metavar="side-angle", type="float", help="Used to correct side view orientation when startcsym makes a mistake.")
	
	(options, args) = parser.parse_args()

	if len(args) != 1:
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
        sys.exit(1)

		
	inputfile = args[0]
	
	nk = options.nkeep
	if options.nkeep < 1:
		print "Error: Number to keep 'nkeep' must be > 0"
		sys.exit(1)

	nimg = EMUtil.get_image_count(inputfile)
	if nimg < 3:
		print "Error: input file too small!"
		sys.exit(1)

	if n/3 < options.nkeep:
		print "Error: Number to keep 'nkeep' should be less than 1/3 number of particles"
		sys.exit(1)

	if not options.nosym:
		files.remove_files("cls*lst")
		os.remove("sym.hed")
		os.remove("sym.img")

		symbest_options = "--sym=" + options.sym + " --mirror=cls0001.lst --nkeep=" + string(options.nkeep)
		e2sysbest.main("e2symbest.py " + symbest_options)
		
		# next: classalignall

if __name__ == "__main__":
    main()
    
