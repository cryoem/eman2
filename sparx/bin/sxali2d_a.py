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
from   global_def     import *
from   user_functions import *
from   optparse       import OptionParser
import sys
def main():
	progname = os.path.basename(sys.argv[0])
	#usage = progname + " stack outdir <maskfile> --ir=inner_radius --ou=outer_radius --rs=ring_step --xr=x_range --yr=y_range --ts=translation_step --center=center --maxit=max_iter --CTF --function=user_function_name --randomize --T0=T0 --F=F --MPI"
	usage = progname + " stack outdir <maskfile> --ir=inner_radius --ou=outer_radius --rs=ring_step --xr=x_range --yr=y_range --ts=translation_step --center=center --maxit=max_iter --CTF --function=user_function_name --restart=0 --MPI"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--ir",    type="float",  default=1,             help="inner radius for rotational correlation > 0 (set to 1)")
	parser.add_option("--ou",    type="float",  default=-1,            help="outer radius for rotational correlation < int(nx/2)-1 (set to the radius of the particle)")
	parser.add_option("--rs",    type="float",  default=1,             help="step between rings in rotational correlation > 0 (set to 1)" ) 
	parser.add_option("--xr",    type="string", default="4 2 1 1",     help="range for translation search in x direction, search is +/xr ")
	parser.add_option("--yr",    type="string", default="-1",          help="range for translation search in y direction, search is +/yr ")
	parser.add_option("--ts",    type="string", default="2 1 0.5 0.25",help="step of translation search in both directions direction, search is -xr, -xr+ts, 0, xr-ts, xr ")
	parser.add_option("--center",type="float",  default=-1,            help="-1. average center method 0.not centered 1.phase approximation 2.cc with Gaussian function 3.cc with donut-shaped image 4.cc with user-defined reference 5.cc with self-rotated average")
	parser.add_option("--maxit", type="float",  default=0,             help="maximum number of iterations (0 means the maximum iterations is 10, but it will automatically stop should the criterion falls")
	parser.add_option("--CTF", action="store_true", default=False,     help="Consider CTF correction during the alignment ")
	parser.add_option("--function", type="string", default="ref_ali2d",help="name of the reference preparation function")
	#parser.add_option("--random_method", type="string", default="",    help="SA: Simulated Annealing   ML: Maximum Likelihood")
	#parser.add_option("--T0",    type="float",  default=1.0,           help="initial temperature for simulated annealing")
	#parser.add_option("--F",     type="float",  default=0.996,         help="cooling rate for simulated annealing")
	parser.add_option("--restart", type="int", default=-1,             help="name of the reference preparation function")
	parser.add_option("--MPI", action="store_true", default=False,     help="use MPI version ")
	(options, args) = parser.parse_args()
	if len(args) < 2 or len(args) >3:
    		print "usage: " + usage
    		print "Please run '" + progname + " -h' for detailed options"
	else:
		
		if len(args) == 2: mask = None
		else:              mask = args[2]

		from applications import ali2d_a
		if options.MPI:
			from mpi import mpi_init
			sys.argv = mpi_init(len(sys.argv),sys.argv)		

			from utilities import init_mpi_bdb
			init_mpi_bdb()

		global_def.BATCH = True
		ali2d_a(args[0], args[1], mask, options.ir, options.ou, options.rs, options.xr, options.yr, options.ts, options.center, options.maxit, options.CTF, options.function, options.restart, options.MPI)
		global_def.BATCH = False

if __name__ == "__main__":
	main()
