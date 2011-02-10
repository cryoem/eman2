#!/usr/bin/env python

#
# Author: Steven Ludtke, 02/10/2011 (sludtke@bcm.edu)
# Copyright (c) 2011 Baylor College of Medicine
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
import sys
from math import *
import os.path
import pyemtbx.options
from pyemtbx.options import intvararg_callback
from pyemtbx.options import floatvararg_callback

def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + """ [options] <inputfile>
	This is a specialized version of e2proc3d.py targeted at performing a limited set of operations on
very large volumes in-place (such as tomograms) which may not readily fit into system memory. Operations are 
performed by reading portions of the image, processing, then writing the portion back to disk. Unlike e2proc3d.py
you may pass only a single operation to the program for each invocation, or behavior will be undefined. It will 
process a single volume in a single file in-place.

"""
	parser = OptionParser(usage)
	

	parser.add_option("--streaksubtract",type="string",help="This will subtract the histogram peak value along a single axis in the volume.",default=None)

	parser.add_option("--process", metavar="processor_name:param1=value1:param2=value2", type="string",
								action="append", help="apply a processor named 'processorname' with all its parameters/values. WARNING: this works by operating on fragments of the overall image at a time, and some processors won't work properly this way.")

	parser.add_option("--mult", metavar="f", type="float", 
								help="Scales the densities by a fixed number in the output")
	
	parser.add_option("--multfile", type="string", action="append",
								help="Multiplies the volume by another volume of identical size. This can be used to apply masks, etc.")
		
	parser.add_option("--add", metavar="f", type="float", 
								help="Adds a constant 'f' to the densities")

	parser.add_option("--trans", metavar="dx,dy,dz", type="string", default=0, help="Translate map by dx,dy,dz ")

	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n", type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")
		
	(options, args) = parser.parse_args()

	print "Sorry, this program still under development. Not functional yet."
	sys.exit(1)

	try:
		hdr=EMData(args[1],0,1)
	except:
		print "ERROR: Can't read input file header"
		sys.exit(1)


	if options.mediansubtract!=None :
		
		

def findmode(img) :
	"""This computes something akin to the mode


if __name__ == "__main__":
	main()
