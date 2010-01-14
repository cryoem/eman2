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
import sys, ConfigParser
def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " config_file.cfg --ite=1 --ali --clus --reali"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option('--ite',      type='int',          default=-1,    help='Start at the iteration number defined')
	parser.add_option('--ali',      action='store_true', default=False, help='Perform global alignement')
	parser.add_option('--clus',     action='store_true', default=False, help='Perform clustering')
	parser.add_option('--reali',    action='store_true', default=False, help='Perform realignment')
	parser.add_option('--decimate', action='store_true', default=False, help='Used decimated data')
	(options, args) = parser.parse_args()

    	if len(args) != 1:
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
		sys.exit()

	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()

	from applications import isc
	global_def.BATCH = True
	isc(args[0], options.ite, options.ali, options.clus, options.reali, options.decimate)
	global_def.BATCH = False
	

if __name__ == "__main__":
	        main()
