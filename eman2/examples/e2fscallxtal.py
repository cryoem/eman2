#!/usr/bin/env python

# Author: Steven Ludtke, 05/18/14 (sludtke@bcm.edu)
# Copyright (c) 2000- Baylor College of Medicine
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

import os
from EMAN2 import *

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <xtalfile>
	Will take a 3-D density computed from a crystal structure or other source and compute an FSC curve between it and the final
	iteration from all refine_XX folders. The provided <xtalfile> MUST have the same sampling, box size and orientation as was
	used for the refinement.
	"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()

	# find all directories named refine_*
	ds=[i for i in os.listdir(".") if "refine_" in i and os.path.isdir(i)]
	
	# iterate over the refine_xx directories
	for d in ds:
		tds=[i for i in os.listdir(d) if len(i)==13 and "threed_" in i and ".hdf" in i]		
		big=max(tds)
		nd=d[7:9]
		ntd=big[7:9]
		
		com="e2proc3d.py {} fsc_{}_{}.txt --calcfsc {}/{}".format(args[0],nd,ntd,d,big)
		print com
		os.system(com)


if __name__ == "__main__":
	main()
