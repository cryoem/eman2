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
	usage = progname + " stack outdir <maskfile> --K=2 --nb_part=5 --opt_method='SSE' --CTF --F=0.9 --max_run=10 --th_nobj=10 --th_stab_max=6.0 --th_stab_min=3.0 --min_dec_K=5 --FK=0 --rand_seed=0 --MPI --CUDA"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--K",              type="int",          default=2,         help="Number of classes for k-means (default 2)")
	parser.add_option("--nb_part",        type="int",          default=5,         help="Number of partitions used to calculate the stability (default 5)")
	parser.add_option("--opt_method",     type='string',       default="SSE",     help="K-means method: SSE (default), cla")
	parser.add_option("--CTF",            action="store_true", default=False,     help="Perform classification using CTF information")
	parser.add_option("--F",              type="float",        default=0.0,       help="Factor to decrease temperature in simulate annealing, ex.: 0.95")
	parser.add_option("--FK",             type="float",        default=0,         help="Cooling factor to decrease K according the number of runs")
	parser.add_option("--max_run",        type="int",          default=50,        help="Maximum number of runs (default 50)")
	parser.add_option("--th_nobj",        type="int",          default=10,        help="Cleanning threshold, classes with number of images < th_nobj are removed (default 10)")
	parser.add_option("--th_stab_max",    type="float",        default=6.0,       help="Stability threshold where K stop to decrease (default 6.0)")
	parser.add_option("--th_stab_min",    type="float",        default=0.0,       help="Stability threshold required to start the next run in percent (default 3.0)")
	parser.add_option("--min_dec_K",      type="int",          default=5,         help="Minimum decrement value to K when the stability is to low (default 5)")
	parser.add_option("--rand_seed",      type="int",          default=0,         help="Random seed")
	parser.add_option("--MPI",            action="store_true", default=False,     help="MPI version")
	parser.add_option("--CUDA",           action="store_true", default=False,     help="CUDA version")
	
	(options, args) = parser.parse_args()
    	if len(args) < 2 or len(args) > 3:
				print "usage: " + usage
        			print "Please run '" + progname + " -h' for detailed options"
	elif(options.opt_method != "cla"  and options.opt_method != "SSE"):
			sys.stderr.write("ERROR: unknown method\n\n")
			sys.exit()
	else:
		if len(args) == 2: mask = None
		else:              mask = args[2]

		if options.K < 2:
			sys.stderr.write('ERROR: K must be > 1 group\n\n')
			sys.exit()

		if options.nb_part < 2:
			sys.stderr.write('ERROR: nb_part must be > 1 partition\n\n')
			sys.exit()

		from  applications  import  k_means_stab
		global_def.BATCH = True
		k_means_stab(args[0], args[1], mask, options.opt_method, options.K, options.nb_part, options.CTF, options.F, options.FK, options.max_run, options.th_nobj, options.th_stab_max, options.th_stab_min, options.min_dec_K, options.rand_seed, options.MPI, options.CUDA)
		global_def.BATCH = False

if __name__ == "__main__":
	        main()
