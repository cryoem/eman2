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
from   global_def import *

from   optparse import OptionParser

def genbuf( prjfile, bufprefix, beg, end, CTF, npad, verbose = 0 ):
	from EMAN2  import newfile_store
	from utilities import get_im
	from time import time
	import os
	if(verbose == 1):  finfo=open( os.path.join(outdir, "progress.txt"), "w" )
	else:              finfo = None
	start_time = time()
	istore = newfile_store( bufprefix, npad, CTF )
	for i in xrange( beg, end ):
		prj = get_im( prjfile, i )
		istore.add_image( prj, prj.get_attr("xform.projection") )
		if( not(finfo is None) and ((i%100==99 or i==end-1))):
			finfo.write( "%6d buffered, time: %10.3f\n" % (i+1, time()-start_time) )
			finfo.flush()

def main():

	import sys
	import os

        arglist = []
        for arg in sys.argv:
		arglist.append( arg )

	progname = os.path.basename(arglist[0])
	usage = progname + " prjstack bufprefix --npad --CTF --verbose"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--CTF",      action="store_true", default=False, help="use CTF")
	parser.add_option("--npad",     type="int",          default=2,     help="times of padding")
	parser.add_option("--verbose",  type="int",          default=0,     help="verbose level: 0 no, 1 yes")

	(options, args) = parser.parse_args( arglist[1:] )

	if( len(args) != 2):
		print "usage: " + usage
		return None

	prjfile = args[0]

	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()

	bufprefix = args[1]
	nprj = EMUtil.get_image_count( prjfile )
	genbuf( prjfile, bufprefix, 0, nprj, options.CTF, options.npad, options.verbose )


if __name__ == "__main__":
	main()
