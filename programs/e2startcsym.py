#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
#
# This software is issued under a joint BSD/GNU license. You may use the
# source code in this file under either license. However, note that the
# complete EMAN2 and SPARX software packages have some GPL dependencies,
# so you are responsible for compliance with the licenses of these packages
# if you opt to use BSD licensing. The warranty disclaimer below holds
# in either instance.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

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
	
	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n", type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
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
