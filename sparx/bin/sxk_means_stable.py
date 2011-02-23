#! /usr/bin/env python

#
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
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
import global_def
from   global_def import *
from   optparse import OptionParser
import sys
def main():
	
	progname = os.path.basename(sys.argv[0])
	usage = progname + " stack outdir <maskfile> --K=2 --nb_part=5 --F=0.9 --T0=5.0 --th_nobj=10 --rand_seed=10 --opt_method=SSE --maxit=1000 --normalize --CTF --CUDA --MPI"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--K",              type="int",          default=2,         help="Number of classes for K-means (default 2)")
	parser.add_option("--nb_part",        type="int",          default=5,         help="Number of partitions used to calculate the stability (default 5)")
	parser.add_option("--F",              type="float",        default=0.0,       help="Cooling factor in simulated annealing, <1.0")
	parser.add_option("--T0",             type="float",        default=0.0,       help="Simulated annealing first temperature")
	parser.add_option("--th_nobj",        type="int",          default=1,         help="Cleanning threshold, classes with number of images < th_nobj are removed (default 1)")
	parser.add_option("--rand_seed",      type="int",          default=0,         help="Random seed")
	parser.add_option("--opt_method",     type='string',       default='SSE',     help="K-means method: SSE (default), cla")
	#parser.add_option("--match",          type='string',       default='bbenum',     help='Algorithm to match partitions: pwa, pair-wise agreement (default), or hh, hierarchical Hungarian algorithm, or bbenum')
	parser.add_option("--maxit",          type="int",          default=1e9,       help="Maximum number of iterations for k-means")
	parser.add_option("--normalize",      action="store_true", default=False,     help="Normalize images under the mask")
	parser.add_option("--CTF",            action="store_true", default=False,     help="Perform classification using CTF information")
	parser.add_option("--CUDA",           action="store_true", default=False,     help="CUDA version")
	parser.add_option("--MPI",            action="store_true", default=False,     help="Use MPI version ")	
	(options, args) = parser.parse_args()
    	if len(args) < 2 or len(args) > 3:
				print "usage: " + usage
        			print "Please run '" + progname + " -h' for detailed options"
	else:
		if len(args) == 2: mask = None
		else:              mask = args[2]

		if options.K < 2:
			sys.stderr.write('ERROR: K must be > 1 group\n\n')
			sys.exit()

		if options.nb_part < 2:
			sys.stderr.write('ERROR: nb_part must be > 1 partition\n\n')
			sys.exit()

		if options.F != 0 and options.T0 == 0:
			sys.stderr.write('ERROR: is F is set you must also set T0\n\n')
			sys.exit()

		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()

		global_def.BATCH = True
		if options.MPI:
			from mpi import mpi_init
			sys.argv = mpi_init(len(sys.argv), sys.argv)
			if options.CUDA:
				from  development import  k_means_stab_MPICUDA_stream_YANG
				k_means_stab_MPICUDA_stream_YANG(args[0], args[1], mask, options.K, options.nb_part, options.F, options.T0, options.th_nobj, options.rand_seed, options.maxit)
			else:
				from  development import  k_means_stab_MPI_stream
				k_means_stab_MPI_stream(args[0], args[1], mask, options.K, options.nb_part, options.F, options.T0, options.th_nobj, options.rand_seed, options.opt_method, options.CTF, options.maxit)
		else:
			if options.CUDA:
				from  development  import  k_means_stab_CUDA_stream
				k_means_stab_CUDA_stream(args[0], args[1], mask, options.K, options.nb_part, options.F, options.T0, options.th_nobj, options.rand_seed, options.maxit)
			else:
				from  development  import  k_means_stab_stream
				k_means_stab_stream(args[0], args[1], mask, options.K, options.nb_part, options.F, options.T0, options.th_nobj, options.rand_seed, options.opt_method, options.CTF, options.maxit)
		global_def.BATCH = False

		if options.MPI:
			from mpi import mpi_finalize
			mpi_finalize()
if __name__ == "__main__":
	        main()
