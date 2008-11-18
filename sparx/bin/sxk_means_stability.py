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
	usage = progname + " stack <maskfile> --K=number_of_classes --nb_part=number_of_partitions --opt_method=optimization_method --CTF --F=factor_temperature"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--K",          type="int",          default=2,         help=" Number of classes (default 2)")
	parser.add_option("--nb_part",    type="int",          default=5,         help=" Number of partitions used to calculate the stability (default 5)")
	parser.add_option("--opt_method", type='string',       default="SSE",     help=" K-means method: SSE (default), cla")
	parser.add_option("--CTF",        action="store_true", default=False,     help=" Perform classification using CTF information")
	parser.add_option("--F",          type="float",        default=0.0,       help=" Factor to decrease temperature in simulate annealing, ex.: 0.9")
	parser.add_option("--debug",      action="store_true", default=False,     help="")
	
	(options, args) = parser.parse_args()
    	if len(args) < 1 or len(args) > 2:
				print "usage: " + usage
        			print "Please run '" + progname + " -h' for detailed options"
	elif(options.opt_method != "cla"  and options.opt_method != "SSE"):
			sys.stderr.write("ERROR: unknown method\n\n")
			sys.exit()
	else:
		if len(args) == 1: mask = None
		else:              mask = args[1]

		if options.K < 2:
			sys.stderr.write('ERROR: K must be > 1 group\n\n')
			sys.exit()

		if options.nb_part < 2:
			sys.stderr.write('ERROR: nb_part must be > 1 partition\n\n')
			sys.exit()

		from  development  import  k_means_stab
		global_def.BATCH = True
		k_means_stab(args[0], mask, options.opt_method, options.K, options.nb_part, options.CTF, options.F, options.debug)
		global_def.BATCH = False

if __name__ == "__main__":
	        main()
