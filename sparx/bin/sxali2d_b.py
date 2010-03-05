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
from   global_def     import *
from   user_functions import *
from   optparse       import OptionParser
import sys
def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " stack outdir nodelist <maskfile> --ir=inner_radius --ou=outer_radius --rs=ring_step --xr=x_range --yr=y_range --ts=translation_step --maxit=max_iter --num_ali=number_of_alignment --CTF --snr=SNR --function=user_function_name --MPI"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--ir",    type="int",  default=1,             help="inner radius for rotational correlation > 0 (set to 1)")
	parser.add_option("--ou",    type="int",  default=-1,            help="outer radius for rotational correlation < int(nx/2)-1 (set to the radius of the particle)")
	parser.add_option("--rs",    type="int",  default=1,             help="step between rings in rotational correlation > 0 (set to 1)" ) 
	parser.add_option("--xr",    type="string", default="4 2 1 1",     help="range for translation search in x direction, search is +/xr ")
	parser.add_option("--yr",    type="string", default="-1",          help="range for translation search in y direction, search is +/yr ")
	parser.add_option("--ts",    type="string", default="2 1 0.5 0.25",help="step of translation search in both directions direction, search is -xr, -xr+ts, 0, xr-ts, xr ")
	parser.add_option("--maxit", type="int",  default=0,             help="maximum number of iterations (0 means the maximum iterations is 10, but it will automatically stop should the criterion falls")
	parser.add_option("--num_ali", type="int",    default=4,           help="Number of alignments performed")
	parser.add_option("--max_merge", type="int",   default=50,          help="The maximum merge allowed")
	parser.add_option("--CTF", action="store_true", default=False,     help="Consider CTF correction during the alignment ")
	parser.add_option("--Fourvar", action="store_true", default=False,     help="Whether to divided by variance")
	parser.add_option("--adw", action="store_true", default=False,     help="Whether to use new CTF correction")
	parser.add_option("--Ng", type="int", default=1,                   help="Number of groups presumably in the dataset")
	parser.add_option("--snr",   type="float",  default=1.0,           help="Signal-to-noise ratio of the dataset")
	parser.add_option("--thr",   type="float",  default=5.0,           help="threshold")
	parser.add_option("--function", type="string", default="ref_ali2d",help="name of the reference preparation function")
	parser.add_option("--CUDA", action="store_true", default=False,    help="use CUDA program")
	parser.add_option("--GPU", type="int", default=0,                  help="number of GPUs available")
	parser.add_option("--proc", type="int", default=4,                 help="number of processors available")
	(options, args) = parser.parse_args()
	
	if len(args) < 3 or len(args) > 4:
    		print "usage: " + usage
    		print "Please run '" + progname + " -h' for detailed options"
	else:
		
		if len(args) == 3: mask = None
		else:              mask = args[3]

		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()

		from development import ali2d_b
		if options.proc > 1:
			from mpi import mpi_init
			sys.argv = mpi_init(len(sys.argv),sys.argv)		

		global_def.BATCH = True
		ali2d_b(args[0], args[1], args[2], mask, options.ir, options.ou, options.rs, options.xr, options.yr, options.ts, options.maxit, options.num_ali, options.max_merge, \
			options.CTF, options.Fourvar, options.adw, options.Ng, options.snr, options.function, options.CUDA, options.GPU, options.proc, options.thr)
		global_def.BATCH = False

if __name__ == "__main__":
	main()
