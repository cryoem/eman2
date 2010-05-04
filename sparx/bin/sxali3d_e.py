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
from   utilities import parse_user_function
import sys


def main():
	arglist = []
	for arg in sys.argv:
		arglist.append( arg )
	progname = os.path.basename(sys.argv[0])
	usage = progname + " stack outdir <maskfile> --ou=outer_radius --delta=angular_bracket --maxit=max_iter --chunk=data_chunk_for_update --center --CTF --snr=SNR --sym=symmetry  --function=user_function --MPI"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--ou",     type="float",        default=-1,      help="outer radius of a circular mask that should encompass the particle< int(nx/2)-1 (set to int(nx/2)-1)")
	parser.add_option("--delta",  type="float",        default=2,       help="angular bracket (set to 2)")
	parser.add_option("--ts",     type="float",        default=0.25,    help="shift bracket (set to 2)")
	parser.add_option("--center", type="float",        default=-1,      help="-1 - average centering method (default) 0 - if you do not want the volume to be centered, 1 - center the volume using cog")
	parser.add_option("--maxit",  type="int",          default=10,      help="maximum number of iterations (set to 10)")
	parser.add_option("--chunk",  type="float",        default=1.0,     help="chunk of data after which the 3-D structure will be updated 0<chunk<=1.0 (set to 1.0)")
	parser.add_option("--CTF",    action="store_true", default=False,   help="Consider CTF correction during the alignments")
	parser.add_option("--snr",    type="float", 	   default=1,       help="SNR > 0.0 (set to 1.0)")
	parser.add_option("--sym",    type="string",       default="c1",    help="symmetry group (set to c1)")
	parser.add_option("--function", type="string", default="ref_ali3d", help="name of the user-supplied reference preparation function")
	parser.add_option("--debug",  action="store_true", default=False,   help="Debug printout")
	parser.add_option("--MPI",    action="store_true", default=False,   help="use MPI version")
	parser.add_option("--fourvar",action="store_true", default=False,   help="compute Fourier variance")
	(options, args) = parser.parse_args(arglist[1:])
	if(len(args) < 2 or len(args) > 3):
    		print "usage: " + usage
    		print "Please run '" + progname + " -h' for detailed options"
	else:
	
		if(len(args) == 2):
			mask = None
		else:
			mask = args[2]

		if options.MPI:
			from mpi import mpi_init
	   		sys.argv = mpi_init(len(sys.argv), sys.argv)

		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()

		from applications import ali3d_e
		global_def.BATCH = True
		ali3d_e(args[0], args[1], mask, options.ou, options.delta, options.ts, options.center, options.maxit, options.CTF, options.snr, options.sym, options.chunk, parse_user_function(options.function), options.fourvar, options.debug, options.MPI)
		global_def.BATCH = False

if __name__ == "__main__":
	main()
