#! /usr/bin/env python
#
# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Tapu Shaikh 2019 (tapu.shaikh@mpi-dortmund.mpg.de)
# Author: Adnan Ali 2019 (adnan.ali@mpi-dortmund.mpg.de)
# Author: Luca Lusnig 2019 (luca.lusnig@mpi-dortmund.mpg.de)
# Author: Toshio Moriya 2019 (toshio.moriya@kek.jp)
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
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

import os
from ..libpy import sp_global_def
from ..libpy.sp_global_def import sxprint, ERROR

from ..libpy.sp_global_def import *
from   optparse import OptionParser
import sys

def run():
	
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
		sxprint( "Usage: " + usage )
		sxprint( "Please run '" + progname + " -h' for detailed options" )
		ERROR( "Invalid number of parameters. Please see usage information above." )
		return

	elif options.link not in chk_link:
		ERROR( "Kind of linkage unknown" )
		return

	elif options.dist not in chk_dist:
		ERROR( "Kind of distance unknown" )
		return

	else:
		if len(args) == 2: maskname = None
		else:              maskname = args[2]

		if sp_global_def.CACHE_DISABLE:
			from ..libpy.sp_utilities import disable_bdb_cache
			disable_bdb_cache()
		from ..libpy.sp_applications  import  HAC_clustering
		sp_global_def.BATCH = True
		HAC_clustering(args[0], args[1], maskname, options.link, options.dist, options.dissimilar)
		sp_global_def.BATCH = False

def main():
	sp_global_def.print_timestamp( "Start" )
	sp_global_def.write_command()
	run()
	sp_global_def.print_timestamp( "Finish" )

if __name__ == "__main__":
	main()
