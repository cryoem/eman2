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

def run():

	progname = os.path.basename(sys.argv[0])
	usage = progname + " stackfile outdir  <maskfile> --K1=Min_number_of_Cluster --K2=Max_number_of_Clusters --opt_method=K-means_method --trials=Number_of_trials_of_K-means --CTF --rand_seed=1000 --maxit=Maximum_number_of_iterations --F=simulated_annealing --T0=simulated_annealing --MPI --CUDA --debug"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--K1",          type="int",          default=2,          help="Mimimum number of clusters")
	parser.add_option("--K2",          type="int",          default=3,          help="Maximum number of clusters")
	parser.add_option("--trials",      type="int",          default=1,          help="Number of trials in K-means (default 1)")
	parser.add_option("--CTF",         action="store_true", default=False,      help="Perform clustering using CTF information")
	parser.add_option("--rand_seed",   type="int",          default=-1,         help="Random seed of initial (default random)" )
	parser.add_option("--maxit",       type="int",          default=100,        help="Mimimum number of iterations within K-means")
	#parser.add_option("--F",           type="float",        default=0.0,        help="Factor to decrease temperature in simulated annealing, ex.: 0.9")
	#parser.add_option("--T0",          type="float",        default=0.0,        help="Initial temperature in simulated annealing, ex: 100")
	parser.add_option("--MPI",         action="store_true", default=False,      help="Use MPI version")
	parser.add_option("--debug",       action="store_true", default=False,      help="Debug output")

	(options, args) = parser.parse_args()
	if len(args) < 2 or len(args) > 3:
		sxprint("usage: " + usage)
		sxprint("Please run '" + progname + " -h' for detailed options")
		ERROR( "Invalid number of parameters used. Please see usage information above." )
		return

	elif options.trials < 1:
		ERROR( "Number of trials should be at least 1" )
		return

	else: 
		if len(args)==2: mask = None
		else:            mask = args[2]

		if options.K1 < 2:
			ERROR( "K1 must be > 1 group" )
			return

		if sp_global_def.CACHE_DISABLE:
			from ..libpy.sp_utilities import disable_bdb_cache
			disable_bdb_cache()
			
		from ..libpy.sp_applications import k_means_groups
		sp_global_def.BATCH = True
		k_means_groups(args[0], args[1], mask, "SSE", options.K1, options.K2, options.rand_seed, options.maxit, options.trials, options.CTF, 0.0, 0.0, options.MPI, False, options.debug)
		sp_global_def.BATCH = False
		
		if options.MPI:
			from mpi import mpi_finalize
			mpi_finalize()

def main():
	sp_global_def.print_timestamp( "Start" )
	sp_global_def.write_command()
	run()
	sp_global_def.print_timestamp( "Finish" )

if __name__ == "__main__":
	main()
