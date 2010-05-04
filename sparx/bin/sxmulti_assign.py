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
from   global_def import *
from   optparse import OptionParser
import sys, ConfigParser

def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " configure_file.cfg"
	
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--ir", type="float", default=1, help="  inner radius for rotational correlation (set to 1)")
	parser.add_option("--ou", type="float", default=-1, help="  outer radius for rotational correlation (set to the radius of the particle)")
	parser.add_option("--rs", type="float", default=1, help="  step between rings in rotational correlation (set to 1)" )
	parser.add_option("--xr", type="float", default=0, help="  range for translation search in x direction, search is +/-xr ")
	parser.add_option("--yr", type="float", default=0, help="  range for translation search in y direction, search is +/-yr ")
	parser.add_option("--ts", type="float", default=1, help="  step of translation search in both directions")
	parser.add_option("--maxit", type="float", default=10, help="  maximum number of iterations (set to 10) ")
	parser.add_option("--min_class_size", type="float", default=1, help="  mininum size of class (set to 1) ")
	parser.add_option("--max_class_size", type="float", default=-1, help="  maximum size of class (set to -1, means can be any number) ")
	parser.add_option("--split_group", action="store_true", default=False, help=" whether to split group is allowed")
	parser.add_option("--combine_group", action="store_true", default=False, help=" whether to combine group is allowed")
	parser.add_option("--CTF", action="store_true", default=False, help=" Consider CTF correction during multiple reference assignment")
	parser.add_option("--MPI", action="store_true", default=False, help="  whether to use MPI version ")

	(options, args) = parser.parse_args()
	if len(args) < 2 or len(args) > 3:
    		print "usage: " + usage
    		print "Please run '" + progname + " -h' for detailed options"
		sys.exit()

	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()
		
	from development import multi_assign
	global_def.BATCH = True
	multi_assign(args[0], args[1], args[2], mask, options.ir, options.ou, options.rs, options.xr, options.yr, options.ts, options.maxit, 
			options.min_class_size, options.max_class_size, options.split_group, options.combine_group, options.CTF, options.MPI)
	global_def.BATCH = False


if __name__ == "__main__":
	main()
