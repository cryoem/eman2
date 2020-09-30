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
#


import os
from ..libpy import sp_global_def
from ..libpy.sp_global_def import sxprint, ERROR
from ..libpy.sp_global_def import *
from   optparse import OptionParser
import sys

import mpi

mpi.mpi_init( 0, [] )


def run():
	arglist = []
	for arg in sys.argv:
		arglist.append( arg )
	progname = os.path.basename(sys.argv[0])
	usage = progname + " stack outdir <maskfile> --ou=outer_radius --delta=angular_bracket --maxit=max_iter --chunk=data_chunk_for_update --center --CTF --snr=SNR --sym=symmetry  --function=user_function --MPI"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--ou",       type="float",        default=-1,      help="outer radius of a circular mask that should encompass the particle< int(nx/2)-1 (set to int(nx/2)-1)")
	parser.add_option("--delta",    type="float",        default=2,       help="angular bracket (set to 2)")
	parser.add_option("--ts",       type="float",        default=2,       help="shift bracket (set to 2)")
	parser.add_option("--center",   type="float",        default=0,       help="-1 - average centering method; 0 - no cetnering of template volume (default), 1 - center the volume using center of gravity")
	parser.add_option("--maxit",    type="int",          default=10,      help="maximum number of iterations (set to 10)")
	parser.add_option("--chunk",    type="float",        default=1.0,     help="chunk of data after which the 3-D structure will be updated 0<chunk<=1.0 (set to 1.0)")
	parser.add_option("--CTF",      action="store_true", default=False,   help="Consider CTF correction during the alignments")
	parser.add_option("--snr",      type="float", 	     default=1,       help="SNR > 0.0 (set to 1.0)")
	parser.add_option("--sym",      type="string",       default="c1",    help="symmetry group (set to c1)")
	parser.add_option("--function", type="string",       default="ref_ali3d", help="name of the user-supplied reference preparation function")
	parser.add_option("--npad",     type="int",          default= 2,      help="padding size for 3D reconstruction")
	parser.add_option("--debug",    action="store_true", default=False,   help="Debug printout")
	parser.add_option("--MPI",      action="store_true", default=False,   help="use MPI version")
	parser.add_option("--fourvar",  action="store_true", default=False,   help="compute Fourier variance")
	parser.add_option("--scipy_minimization",  action="store_true", default=False,   help="use scipy minimization instead of amoeba")
	(options, args) = parser.parse_args(arglist[1:])
	if(len(args) < 2 or len(args) > 3):
		sxprint("usage: " + usage)
		sxprint("Please run '" + progname + " -h' for detailed options")
		ERROR( "Invalid number of parameters used. Please see usage information above." )
		return

	else:

		if(len(args) == 2):
			mask = None
		else:
			mask = args[2]

		if sp_global_def.CACHE_DISABLE:
			from ..libpy.sp_utilities import disable_bdb_cache
			disable_bdb_cache()

		
		sp_global_def.BATCH = True
		if options.fourvar:
			from ..libpy.sp_development import nlocal_ali3d_MPI
			nlocal_ali3d_MPI(args[0], args[1], mask, options.ou, options.delta, options.ts, options.center, options.maxit,
			options.CTF, options.snr, options.sym, options.chunk, options.function, options.fourvar,
			options.npad, options.debug)
		elif options.scipy_minimization:
			if options.MPI:
				from ..libpy.sp_applications import local_ali3d_MPI_scipy_minimization
				local_ali3d_MPI_scipy_minimization(args[0], args[1], mask, options.ou, options.delta, options.ts, options.center, options.maxit,
				options.CTF, options.snr, options.sym, options.chunk, options.function, options.fourvar,
				options.npad, options.debug)
		else:
			from ..libpy.sp_applications import local_ali3d
			local_ali3d(args[0], args[1], mask, options.ou, options.delta, options.ts, options.center, options.maxit,
			options.CTF, options.snr, options.sym, options.chunk, options.function, options.fourvar,
			options.npad, options.debug, options.MPI)
		sp_global_def.BATCH = False

def main():
	sp_global_def.print_timestamp( "Start" )
	sp_global_def.write_command()
	run()
	sp_global_def.print_timestamp( "Finish" )
	mpi.mpi_finalize()

if __name__ == "__main__":
	main()
