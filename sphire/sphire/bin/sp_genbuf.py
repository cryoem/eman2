#!/usr/bin/env python
#
# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Tapu Shaikh 2019 (tapu.shaikh@mpi-dortmund.mpg.de)
# Author: Adnan Ali 2019 (adnan.ali@mpi-dortmund.mpg.de)
# Author: Luca Lusnig 2019 (luca.lusnig@mpi-dortmund.mpg.de)
# Author: Toshio Moriya 2019 (toshio.moriya@kek.jp)
# Author: Pawel A.Penczek 05/27/2009 (Pawel.A.Penczek@uth.tmc.edu)
#
# Copyright (c) 2019 Max Planck Institute of Molecular Physiology
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

from builtins import range
from ..libpy import sp_global_def
from ..libpy.sp_global_def import sxprint, ERROR

from ..libpy.sp_global_def import *

from   optparse import OptionParser

from EMAN2  import newfile_store
from ..libpy.sp_utilities import get_im
from time import time
import os

import sys

def genbuf( prjfile, bufprefix, beg, end, CTF, npad, verbose = 0 ):
	
	if(verbose == 1):  
		finfo=open( os.path.join(outdir, "progress.txt"), "w" )
	else:
		finfo = None
	
	start_time = time()
	istore = newfile_store( bufprefix, npad, CTF )
	for i in range( beg, end ):
		prj = get_im( prjfile, i )
		istore.add_image( prj, prj.get_attr("xform.projection") )
		if( not(finfo is None) and ((i%100==99 or i==end-1))):
			finfo.write( "%6d buffered, time: %10.3f\n" % (i+1, time()-start_time) )
			finfo.flush()

def run():
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
		sxprint( "Usage: " + usage )
		sxprint( "Please run \'" + progname + " -h\' for detailed options" )
		ERROR( "Invalid number of parameters used. Please see usage information above." )
		return

	prjfile = args[0]

	if sp_global_def.CACHE_DISABLE:
		from ..libpy.sp_utilities import disable_bdb_cache
		disable_bdb_cache()

	bufprefix = args[1]
	nprj = EMUtil.get_image_count( prjfile )
	genbuf( prjfile, bufprefix, 0, nprj, options.CTF, options.npad, options.verbose )

def main():
	sp_global_def.print_timestamp( "Start" )
	sp_global_def.write_command()
	run()
	sp_global_def.print_timestamp( "Finish" )
if __name__ == "__main__":
	main()
