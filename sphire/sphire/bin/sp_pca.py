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
from optparse import OptionParser
from EMAN2_cppwrap import *

import os
import sys
import mpi

mpi.mpi_init( 0, [] )


      
def run():
	progname = os.path.basename(sys.argv[0])
	usage = progname + "  input_stack output_stack --subavg=average_image --rad=mask_radius --nvec=number_of_eigenvectors --incore --mask=maskfile --shuffle --usebuf --MPI"
	parser = OptionParser(usage, version=SPARXVERSION)
	parser.add_option("--subavg",  type="string",       default="",    help="subtract average")
	parser.add_option("--rad",     type="int",          default=-1,    help="radius of mask")
	parser.add_option("--nvec",    type="int",          default=1,     help="number of eigenvectors")
	parser.add_option("--mask",    type="string",       default="",    help="mask file" )
	parser.add_option("--genbuf",  action="store_true", default=False, help="use existing buffer")
	parser.add_option("--shuffle", action="store_true", default=False, help="use shuffle")
	parser.add_option("--incore",  action="store_true", default=False, help="no buffer on a disk" )
	parser.add_option("--MPI",     action="store_true", default=False, help="run mpi version" )

	(options, args) = parser.parse_args()

	input_stacks = args[0:-1]
	output_stack = args[-1]

	if options.nvec is None:
		ERROR( "Error: number of components is not given" )
		return

	isRoot = True
	if options.MPI:
		isRoot = ( mpi.mpi_comm_rank(mpi.MPI_COMM_WORLD) == 0 )
		
	if sp_global_def.CACHE_DISABLE:
		from ..libpy.sp_utilities import disable_bdb_cache
		disable_bdb_cache()

	from ..libpy.sp_applications import pca
	sp_global_def.BATCH = True
	vecs = []
	vecs = pca(input_stacks, options.subavg, options.rad, options.nvec, options.incore, options.shuffle, not(options.genbuf), options.mask, options.MPI)

	if isRoot:
		for i in range(len(vecs)):
			vecs[i].write_image(output_stack, i)
	
	sp_global_def.BATCH = False

def main():
	sp_global_def.print_timestamp( "Start" )
	sp_global_def.write_command()
	run()
	sp_global_def.print_timestamp( "Finish" )
	mpi.mpi_finalize()

if __name__ == "__main__":
	main()
	
