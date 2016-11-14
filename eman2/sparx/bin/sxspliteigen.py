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


from EMAN2 import *


def main():

	import sys

        arglist = []
        for arg in sys.argv:
		arglist.append( arg )

	progname = os.path.basename(arglist[0])
	usage = progname + " eigvol EIG_prefix (output volumes multiplied by sqrt(eigval), if set"
	parser = OptionParser(usage,version=SPARXVERSION)

	(options, args) = parser.parse_args( arglist[1:] )

	if( len(args) != 2):
		print "usage: " + usage
		return None

	from math import sqrt
	nimage = EMUtil.get_image_count( args[0] )

	for i in xrange(nimage) :
	        data = EMData()
	        data.read_image( args[0], i )

	        eigval = data.get_attr_default( "eigval", 1.0 )
	        Util.mul_scalar(data , sqrt(eigval) )

	        fname = args[1] + ("%04d_pos.hdf" % (i+1) )
	        data.write_image( fname )

	        fname = args[1] + ("%04d_neg.hdf" % (i+1) )
	        Util.mul_scalar(data , -1 )
	        data.write_image( fname )


if __name__ == "__main__":
	main()
