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
from   optparse   import OptionParser

import sys
import mpi

mpi.mpi_init( 0, [] )

def run():
	progname = os.path.basename(sys.argv[0])
	usage    = progname + " stack outdir --ir --ou --delta --dpsi --lf --hf --rand_seed --maxit --debug --noweights --trials --given --first_zero --weights --MPIGA--pcross --pmut --maxgen --MPI --trials"
	parser   = OptionParser(usage, version = SPARXVERSION)
	parser.add_option("--ir",         type="float",        default=-1,       help="Inner radius of particle (set to 1)")
	parser.add_option("--ou",         type="float",        default=-1,       help="Outer radius of particle < int(nx/2)-1")
	parser.add_option("--delta",      type="float",        default=5.0,      help="Angle step" )
	parser.add_option("--dpsi",       type="int",          default=1,        help="Angle accuracy for sinogram (set to 1)")
	parser.add_option("--lf",         type="float",        default=0.0,      help="Filter, minimum frequency (set to 0.0)")
	parser.add_option("--hf",         type="float",        default=0.5,      help="Filter, maximum frequency (set to 0.5)")
	parser.add_option("--given",      action="store_true", default=False,    help="Start from given projections orientation (set to False, means start with randomize orientations)")
	parser.add_option("--rand_seed",  type="int",          default=-1,       help="Random seed of initial orientations (if set to randomly)")
	parser.add_option("--maxit",      type="int",          default=100,      help="Maximum number of iterations ")
	parser.add_option("--debug",      action="store_true", default=False,    help="Help to debug")
	parser.add_option("--first_zero", action="store_true", default=False,    help="Assign the first projection orientation to 0")
	parser.add_option("--noweights",  action="store_true", default=False,    help="Use Voronoi weighting (by default use weights)")
	parser.add_option("--MPI",        action="store_true", default=False,    help="MPI version")
	parser.add_option("--trials",     type="int",          default=1,        help="Number of trials for the MPI version")
	parser.add_option("--MPIGA",      action="store_true", default=False,    help="MPI version (Genetic algorithm)")
	parser.add_option("--pcross",     type="float",        default=0.95,     help="Cross-over probability (set to 0.95)")
	parser.add_option("--pmut",       type="float",        default=0.05,     help="Mutation probability (set to 0.05)")
	parser.add_option("--maxgen",     type="int",          default=10,       help="Maximum number of generations (set to 10)")

	(options, args) = parser.parse_args()

	if len(args) != 2:
		sxprint( "usage: " + usage )
		sxprint( "Please run '" + progname + " -h' for detailed options" )
		ERROR( "Invalid number of parameters used. Please see usage information above." )

	else:
		if options.maxit < 1: 
			options.maxit = 1
		
		if options.noweights: 
			weights = False
		else:
			weights = True

		if sp_global_def.CACHE_DISABLE:
			from ..libpy.sp_utilities import disable_bdb_cache
			disable_bdb_cache()

		if options.MPIGA:
			from ..libpy.sp_development import cml2_main_mpi
			sp_global_def.BATCH = True
			cml2_main_mpi(args[0], args[1], options.ir, options.ou, options.delta, options.dpsi, 
				      options.lf, options.hf, options.rand_seed, options.maxit, options.given, options.first_zero, 
				      weights, options.debug, options.maxgen, options.pcross, options.pmut)
			sp_global_def.BATCH = False

		elif options.MPI:

			from ..libpy.sp_applications import cml_find_structure_MPI2
			sp_global_def.BATCH = True
			cml_find_structure_MPI2(args[0], args[1], options.ir, options.ou, options.delta, options.dpsi, 
				    options.lf, options.hf, options.rand_seed, options.maxit, options.given, options.first_zero, 
				    weights, options.debug, options.trials)
			sp_global_def.BATCH = False
			
		else:
			from ..libpy.sp_applications import cml_find_structure_main
			sp_global_def.BATCH = True
			cml_find_structure_main(args[0], args[1], options.ir, options.ou, options.delta, options.dpsi, 
				    options.lf, options.hf, options.rand_seed, options.maxit, options.given, options.first_zero, 
				    weights, options.debug, options.trials)
			sp_global_def.BATCH = False

def main():
	sp_global_def.print_timestamp( "Start" )
	sp_global_def.write_command()
	run()
	sp_global_def.print_timestamp( "Finish" )
	mpi.mpi_finalize()

if __name__ == "__main__":
	main()
