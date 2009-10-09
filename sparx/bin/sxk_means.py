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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#


import os
import global_def
from   global_def import *
from   optparse import OptionParser
import sys
def main():
	
	progname = os.path.basename(sys.argv[0])
	usage = progname + " stack outdir <maskfile> --K=number_of_classes --trials=num_trials --opt_method=optimization_method --maxit=max_iter --CTF --rand_seed=1000 --crit=criterion_name --F=factor_temperature --T0=init_temperature --normalize --MPI --CUDA"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--K",          type="int",          default=2,         help="Number of classes (default 2)")
	parser.add_option("--trials",     type="int",          default=1,         help="Number of trials of K-means (default 1)")
	parser.add_option("--opt_method", type='string',       default="cla",     help="K-means method: SSE (default), cla")
	parser.add_option("--maxit",      type="int",          default=100,       help="Maximum number of iterations within K-means")
	parser.add_option("--CTF",        action="store_true", default=False,     help="Perform classification using CTF information")
	parser.add_option("--rand_seed",  type="int",          default=-1,        help="random seed of initial (default random)" )
	parser.add_option("--crit",       type="string",       default="D",       help="Name of criterions: Coleman [C], Harabasz[H], Davies-Bouldin[D], All [all]")
	parser.add_option("--F",          type="float",        default=0.0,       help="Cooling in simulate annealing, ex.: 0.9")
	parser.add_option("--T0",         type="float",        default=0.0,       help="Initial temperature in simulate annealing, ex: 100")
	parser.add_option("--MPI",        action="store_true", default=False,     help="Use MPI version")
	parser.add_option("--CUDA",       action="store_true", default=False,     help="Use CUDA version")
	parser.add_option("--debug",      action="store_true", default=False,     help="")
	parser.add_option("--normalize",  action="store_true", default=False,     help="Normalize images under the mask")
	
	(options, args) = parser.parse_args()
    	if len(args) < 2 or len(args) > 3:
				print "usage: " + usage
        			print "Please run '" + progname + " -h' for detailed options"
	elif options.trials < 1:
			sys.stderr.write("ERROR: Number of trials should be at least 1.\n\n")
			sys.exit()
	elif(options.opt_method != "cla"  and options.opt_method != "SSE"):
			sys.stderr.write("ERROR: unknown method\n\n")
			sys.exit()
	else:
		if len(args)==2: mask = None
		else:            mask = args[2]

		if options.K < 2:
			sys.stderr.write('ERROR: K must be > 1 group\n\n')
			sys.exit()

		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()
		from  applications  import  k_means_main
		global_def.BATCH = True
		k_means_main(args[0], args[1], mask, options.opt_method, options.K, options.rand_seed, options.maxit, options.trials, options.crit, options.CTF, options.F, options.T0, options.MPI, options.CUDA, options.debug, options.normalize)
		global_def.BATCH = False

if __name__ == "__main__":
	        main()
