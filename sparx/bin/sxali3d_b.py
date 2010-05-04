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
from applications import ali3d_b
from optparse import OptionParser
import sys


def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " stack ref_vol outdir --ir=inner_radius --ou=outer_radius --rs=ring_step --xr=x_range --yr=y_range  --ts=translation_step --delta=angular_step --maxit=max_iter --ref_a=S --sym=c1"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--ir", type="float", default=1, help="  inner radius for rotational correlation > 0 (set to 1)")
	parser.add_option("--ou", type="float", default=-1, help="  outer radius for rotational correlation <nx-1 (set to the radius of the particle)")
	parser.add_option("--rs", type="float", default=1, help="  step between rings in rotational correlation >0 (set to 1)" )	
	parser.add_option("--xr", type="float", default=0, help="  range for translation search in x direction, search is +/xr ")
	parser.add_option("--yr", type="float", default=0, help="  range for translation search in y direction, search is +/yr ")
	parser.add_option("--ts", type="float", default=1, help="  step of translation search in both directions direction, search is -xr, -xr+ts, 0, xr-ts, xr ")
	parser.add_option("--delta", type="float", default=15, help="  initial angular step ")
	parser.add_option("--maxit", type="float", default=10, help="  maximum number of iterations (set to 10) ")
	parser.add_option("--ref_a", type="string", default="S", help="  method for generating the quasi-uniformly distributed projection directions (default S) ")
	parser.add_option("--sym", type="string", default="c1", help=" The symmetry of the 3D structure to be constructed, it can be any of Cn, or Dn")
	(options, args) = parser.parse_args()
	if len(args) != 3:
    		print "usage: " + usage
    		print "Please run '" + progname + " -h' for detailed options"
	else:

		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()

		ali3d_b(args[0], args[1], args[2], options.ir, options.ou, options.rs, options.xr, options.yr, options.ts, options.delta, options.maxit, options.ref_a, options.sym)


if __name__ == "__main__":
	main()
