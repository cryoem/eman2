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
def main():

	import os,sys

	arglist = []
	for arg in sys.argv:
		arglist.append( arg )

	progname = os.path.basename(arglist[0])
	usage = progname + " volume binarymask smoothmask --variance --repair"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--variance", type="string",    default=None,         help="variance map")
	parser.add_option("--repair",   type="string",    default="repair.hdf", help="repair map")

	(options, args) = parser.parse_args( arglist[1:] )

	if( len(args) != 3):
		print "usage: " + usage
		return None


	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()
	from utilities import get_im
	from morphology import adaptive_mask, binarize
	v = get_im( args[0] )
	m = adaptive_mask( v , 2.0, 2)
	bm = binarize( m, 0.5 )
	bm.write_image( args[1] )

	adaptive_mask( v , 2.0, 2, 9, 3).write_image( args[2] )
	if(options.variance != None):
		from fundamentals import rot_avg_image
		from morphology import square_root
		trovc = rot_avg_image(get_im(options.variance))
		nc = trovc.get_xsize()//2
		trovc /= trovc[nc,nc,nc]
		square_root(trovc).write_image( options.repair )


if __name__ == "__main__":
	main()

