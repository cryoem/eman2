#!/usr/bin/env python

#
# Author: Steven Ludtke, 02/03/2007 (sludtke@bcm.edu)
# Copyright (c) 2000-2007 Baylor College of Medicine
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

# e2classesbysim.py  03/02/2007	Steven Ludtke
# This program takes a similarity matrix and classifies particles

from EMAN2 import *
from optparse import OptionParser
from math import *
import os
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <sim mx>
	Classifies particles based on a similarity matrix"""
	parser = OptionParser(usage=usage,version=EMANVERSION)

	#parser.add_option("--apix", "-A", type="float", help="A/voxel", default=1.0)
	#parser.add_option("--box", "-B", type="string", help="Box size in pixels, <xyz> or <x>,<y>,<z>")
	#parser.add_option("--het", action="store_true", help="Include HET atoms in the map", default=False)
	parser.add_option("--first", "-f", type="int", help="First particle to classify", default=-1)
	parser.add_option("--last", "-l", type="int", help="Last particle to classify", default=-1)
	
#	parser.add_option("--align",type="string",help="The name of an 'aligner' to use prior to comparing the images", default=None)
#	parser.add_option("--saveali",action="store_true",help="Save alignment values, output is c x r x 4 instead of c x r x 1",default=False)
	(options, args) = parser.parse_args()
	
	if len(args)<1 : parser.error("Input and output files required")
	
	E2n=E2init(sys.argv)
	
	
	E2end(E2n)

def classify(simmx,neach):
	"""Classifies particles based on a similarity matrix. References in columns, particles in rows.
	neach is the number of classes to put each particle in."""
	
	
if __name__ == "__main__":
    main()
