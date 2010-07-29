#!/usr/bin/env python
#
# Author: Pawel A.Penczek and Edward H. Egelman 05/27/2009 (Pawel.A.Penczek@uth.tmc.edu)
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
# Copyright (c) 2008-Forever The University of Virginia
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#

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
