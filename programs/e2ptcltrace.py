#!/usr/bin/env python

#
# Author: Steven Ludtke, 1/17/2011 (sludtke@bcm.edu)  (rewrote older broken program)
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

# e2ptcltrace.py  01/17/11  Steven Ludtke
# This program will follow particles through refinement and assess how self consistent particle orientation assignments are


from EMAN2 import *
from optparse import OptionParser
from math import *
import os
import sys
import copy

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options]

	This program traces the orientation of particles through multiple iterations.
"""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--sym",type="string",help="The symmetry to be used to nearness testing and, if it is specified, reduction", default="c1")
	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n", type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	print "WARNING: this program may currently be broken. It is on a TODO list to fix..."
	
	(options, args) = parser.parse_args()
	
	E2n=E2init(sys.argv)
	
	
	E2end(E2n)



if __name__ == "__main__":
    main()
