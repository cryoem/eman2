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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#


import os
import global_def
from   global_def import *
from   optparse import OptionParser
import sys
def main():
	
	progname = os.path.basename(sys.argv[0])
	usage = progname + " stack dendoname <maskfile> --link=kind_of_link --dist=kind_of_dist --dissimilar"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--link",       type='string',       default="single",     help="Kind of linkage: single, complete, average (default single)")
	parser.add_option("--dist",       type='string',       default="sim_SqEuc",  help="Kind of distance: SqEuc, CCC (default SqEuc)")
	parser.add_option("--dissimilar", action='store_true', default=False,        help="Change the distance to the negative value (default False)")

	chk_link = ['single', 'complete', 'average']
	chk_dist = ['SqEuc', 'CCC']
	
	(options, args) = parser.parse_args()
    	if len(args) < 2 or len(args) > 3:
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
	elif options.link not in chk_link:
		sys.stderr.write('ERROR: Kind of linkage unknown.\n\n')
		sys.exit()
	elif options.dist not in chk_dist:
		sys.stderr.write('ERROR: Kind of distance unknown.\n\n')
		sys.exit()
	else:
		if len(args) == 2: maskname = None
		else:              maskname = args[2]

		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()
		from  applications  import  HAC_clustering
		global_def.BATCH = True
		HAC_clustering(args[0], args[1], maskname, options.link, options.dist, options.dissimilar)
		global_def.BATCH = False

if __name__ == "__main__":
	        main()
