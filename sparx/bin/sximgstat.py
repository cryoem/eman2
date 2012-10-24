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
from   global_def     import *

def main():
	import os
	import sys
	from optparse    import OptionParser

	arglist = []
	for arg in sys.argv:
		arglist.append( arg )

	progname = os.path.basename( arglist[0] )
	usage = progname + "  stack1 <stack2> <mask> --ccc --fsc file --inf --rad=r"
	parser = OptionParser(usage, version=SPARXVERSION)

	parser.add_option( "--ccc", action="store_true", default=False, help="print cross corelation coefficient" )
	parser.add_option( "--fsc", type="string",       default="",    help="calculate resolution curve" )
	parser.add_option( "--inf", action="store_true", default=False, help="print basic infomation of the img" )
	parser.add_option( "--rad", type="int",          default=-1,    help="radius of operation" )

        (options,args) = parser.parse_args( arglist[1:] )
     
	if len(args)<1 or len(args)>3:
    		print "usage: " + usage
    		print "Please run '" + progname + " -h' for detailed options"
		sys.exit(-1)


	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()
	from applications import imgstat
	global_def.BATCH = True
	imgstat( args, options.ccc, options.fsc, options.inf, options.rad )

if __name__=="__main__":
	main()
