#!/usr/bin/env python
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
# Copyright (c) 2012 The University of Texas - Houston Medical School
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
	usage = progname + " prj_stack .. average eigvol output_factcoords --rad=radius --neigvol=number_of_eigvol  --CTF"
	parser = OptionParser(usage, version=SPARXVERSION)
	parser.add_option("--rad",       type="int",    default=-1,     help="radius of mask")
	parser.add_option("--neigvol",   type="int",    default=-1,     help="number of eigvenvectors to use (default all)")
	parser.add_option("--fl",        type="float",  default=0.0,    help="cut-off frequency of hyperbolic tangent low-pass Fourier filter")
	parser.add_option("--aa",        type="float",  default=0.0,    help="fall-off of hyperbolic tangent low-pass Fourier filter")
	parser.add_option("--CTF",       action="store_true", default=False,  help="Use CTF")
	parser.add_option("--MPI",       action="store_true",           help="use MPI")

	(options, args) = parser.parse_args()

	if( len(args) < 4 ):
		sxprint( "Usage: " + usage )
		sxprint( "Please run \'" + progname + " -h\' for detailed options" )
		ERROR( "Invalid number of parameters used. Please see usage information above." )
		return

	else:
		stacks = args[0:-3]
		avgvol = args[-3]
		eigvol = args[-2]
		output = args[-1]
		
		if options.rad < 0:
			ERROR( "Mask radius is not given" )
			return

		if sp_global_def.CACHE_DISABLE:
			from ..libpy.sp_utilities import disable_bdb_cache
			disable_bdb_cache()

		from ..libpy.sp_utilities import get_im
		sp_global_def.BATCH = True
		
		if( get_im( stacks[0]).get_zsize() == 1 and get_im( eigvol).get_zsize() > 1):
			from ..libpy.sp_applications import factcoords_prj
			factcoords_prj(stacks, avgvol, eigvol, output, options.rad, options.neigvol, options.fl, options.aa, options.CTF, options.MPI)
		else:
			from ..libpy.sp_applications import factcoords_vol
			factcoords_vol(stacks, avgvol, eigvol, output, options.rad, options.neigvol, options.fl, options.aa, options.MPI)
		sp_global_def.BATCH = False
		
def main():
	sp_global_def.print_timestamp( "Start" )
	sp_global_def.write_command()
	run()
	sp_global_def.print_timestamp( "Finish" )
	mpi.mpi_finalize()

if __name__ == "__main__":
	main()
