#! /usr/bin/env python

#
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
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

import os
import global_def
from global_def import *
from applications import ali2d_rac
from optparse import OptionParser
import sys
def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " stack <maskfile> --ir=inner_radius --ou=outer_radius --rs=ring_step --nclass=number_of_classes --maxit=max_iter --maxin=max_internal --mirror=1 --rand_alpha --rand_seed=1000 --MPI"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--ir", type="float", default=1, help="  inner radius for rotational correlation > 0 (set to 1)")
	parser.add_option("--ou", type="float", default=-1, help="  outer radius for rotational correlation <nx-1 (set to the radius of the particle)")
	parser.add_option("--rs", type="float", default=1, help="  step between rings in rotational correlation >0 (set to 1)" )	
	parser.add_option("--nclass", type="float", default=2, help=" number of classes (set to 2) ")
	parser.add_option("--maxit", type="float", default=10, help="  maximum number of iterations (set to 10) ")
	parser.add_option("--maxin", type="float", default=10, help="  maximum number of interal (per group) iterations (set to 10) ")
	parser.add_option("--check_mirror", action="store_true", default=False, help="  whether to check mirror (set to False) ")
	parser.add_option("--rand_alpha", action="store_true", default=False, help=" start with random alpha")
	parser.add_option("--rand_seed", type="int", default=1000, help=" random seed of initial (set to 1000)" )
	parser.add_option("--MPI", action="store_true", default=False,     help="  whether using MPI version ")
	(options, args) = parser.parse_args()
	if len(args) < 1 or len(args)> 2:
    		print "usage: " + usage
    		print "Please run '" + progname + " -h' for detailed options"
	else:
		if len(args)==1:
			mask = None
		else:
			mask = args[1]

	        if options.MPI:
			from mpi import mpi_init
   			sys.argv = mpi_init(len(sys.argv), sys.argv)

		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()

		global_def.BATCH = True
		ali2d_rac(args[0], mask, options.ir, options.ou, options.rs, options.nclass, options.maxit, options.maxin, options.check_mirror, options.rand_seed, options.rand_alpha, options.MPI)
		global_def.BATCH = False
		
if __name__ == "__main__":
	main()
