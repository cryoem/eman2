#!/usr/bin/env python

#
# Author: Steven Ludtke, 10/12/2009 (sludtke@bcm.edu)
# Copyright (c) 2000-2009 Baylor College of Medicine
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

# e2refmerge.py  10/12/2009	Steven Ludtke
# This program computes a similarity matrix between two sets of images

from EMAN2 import *
from math import *
import os
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <ref imgs> <ref self-simmx> <merged>
	
	THIS PROGRAM IS OBSOLETE

	This program is used by e2simmx2stage.py

	Given a stack of reference projections and a self similarity matrix, this program will find similar subsets of the references, align them
	to each other and average them together, producing a new, smaller, set of averaged projections suitable for preliminary classification of
	images. Classification of images is made based on image similarity, not projection geometry, so this should be a fairly reliable scheme for
	two-stage orientation determination. Each averaged projection will also contain metadata identifying which projections were used to produce
	it, and thus which should be checked in second-stage classification. This takes what is ordinarily an order n calculation and makes it an
	3.5*sqrt(3n) calculation."""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--nmerged",type=int,help="Number of merged references to generate. Default = sqrt(#proj*nbest)",default=0)
	parser.add_argument("--nbest",type=int,help="This will associate each projection with the best N merged references",default=3)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	#parser.add_argument("--apix", "-A", type="float", help="A/voxel", default=1.0)
	#parser.add_argument("--box", "-B", type="string", help="Box size in pixels, <xyz> or <x>,<y>,<z>")
	#parser.add_argument("--het", action="store_true", help="Include HET atoms in the map", default=False)

	(options, args) = parser.parse_args()
	
	if len(args)<3 : parser.error("Please specify <ref img file> <self simmx file> <merged output file>")
	
	
	
if __name__ == "__main__":
    main()

