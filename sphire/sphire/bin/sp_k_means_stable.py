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
	
	progname = os.path.basename(sys.argv[0])
	usage = progname + " stack outdir <maskfile> --K=2 --nb_part=5  --th_nobj=10 --rand_seed=10 --opt_method=SSE --maxit=1000 --normalize --CTF  --MPI"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--K",              type="int",          default=2,         help="Number of classes for K-means (default 2)")
	parser.add_option("--nb_part",        type="int",          default=5,         help="Number of partitions used to calculate the stability (default 5)")
	#parser.add_option("--F",              type="float",        default=0.0,       help="Cooling factor in simulated annealing, <1.0")
	#parser.add_option("--T0",             type="float",        default=0.0,       help="Simulated annealing first temperature")
	parser.add_option("--th_nobj",        type="int",          default=1,         help="Cleanning threshold, classes with number of images < th_nobj are removed (default 1)")
	parser.add_option("--rand_seed",      type="int",          default=0,         help="Random seed")
	#parser.add_option("--opt_method",     type='string',       default='SSE',     help="K-means method: SSE (default), cla")
	#parser.add_option("--match",          type='string',       default='bbenum',     help='Algorithm to match partitions: pwa, pair-wise agreement (default), or hh, hierarchical Hungarian algorithm, or bbenum')
	parser.add_option("--maxit",          type="int",          default=1e9,       help="Maximum number of iterations for k-means")
	parser.add_option("--normalize",      action="store_true", default=False,     help="Normalize images under the mask")
	parser.add_option("--CTF",            action="store_true", default=False,     help="Perform classification using CTF information")
	#parser.add_option("--CUDA",           action="store_true", default=False,     help="CUDA version")
	parser.add_option("--MPI",            action="store_true", default=False,     help="Use MPI version ")	
	(options, args) = parser.parse_args()
	if len(args) < 2 or len(args) > 3:
		sxprint("usage: " + usage)
		sxprint("Please run '" + progname + " -h' for detailed options")
		ERROR( "Invalid number of parameters used. Please see usage information above." )
		return
	else:
		if len(args) == 2:
			mask = None
		else:
			mask = args[2]

		if options.K < 2:
			ERROR( "K must be > 1 group" )
			return

		if options.nb_part < 2:
			ERROR( "nb_part must be > 1 partition" )
			return

		if sp_global_def.CACHE_DISABLE:
			from ..libpy.sp_utilities import disable_bdb_cache
			disable_bdb_cache()

		sp_global_def.BATCH = True
		if options.MPI:
			from ..libpy.sp_statistics import  k_means_stab_MPI_stream
			k_means_stab_MPI_stream(args[0], args[1], mask, options.K, options.nb_part, 0.0, 0.0, options.th_nobj, options.rand_seed, "SSE", options.CTF, options.maxit)
		else:
			from ..libpy.sp_statistics  import  k_means_stab_stream
			k_means_stab_stream(args[0], args[1], mask, options.K, options.nb_part, 0.0, 0.0, options.th_nobj, options.rand_seed, "SSE", options.CTF, options.maxit)
		sp_global_def.BATCH = False


def main():
	sp_global_def.print_timestamp( "Start" )
	sp_global_def.write_command()
	run()
	sp_global_def.print_timestamp( "Finish" )
	mpi.mpi_finalize()

if __name__ == "__main__":
	main()
