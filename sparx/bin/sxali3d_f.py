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
from   global_def import *
from   optparse import OptionParser
import sys


def main():
        arglist = []
        for arg in sys.argv:
	    arglist.append( arg )

	progname = os.path.basename(arglist[0])
	usage = progname + " stack ref_vol outdir <maskfile> --r=radius --CTF --snr=SNR --dtheta=angular_bracket --maxit=max_iter --chunk=data_chunk_for_update --MPI"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--r", type="float", default=-1, help="  radius  <nx-1")
	parser.add_option("--snr", type="float", default=1, help="  SNR  >0.0")
	parser.add_option("--dtheta", type="float", default=2, help="  angular bracket ")
	parser.add_option("--maxit", type="int", default=10, help="  maximum number of iterations (set to 10) ")
	parser.add_option("--chunk", type="float", default=1.0, help="  chunk of data after which the 3D will be updated 0<chunk<=1.0 (default 1.0) ")
	parser.add_option("--CTF", action="store_true", default=False, help="  Consider CTF correction during the alignment ")
	parser.add_option("--symmetry", type="string", default="c1", help="  symmetry group ")
	parser.add_option("--MPI", action="store_true", default=False,     help="  whether using MPI version ")
	(options, args) = parser.parse_args(arglist[1:])
	if(len(args) < 3 or len(args) > 4):
    		print "usage: " + usage
    		print "Please run '" + progname + " -h' for detailed options"

        prj_stack = args[0]
        ref_vol   = args[1]
        outdir    = args[2]
        if( len(args)==3 ) :
        	maskfile = None
        else:
        	maskfile = args[3]

	from applications import ali3d_f

	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()

        if options.MPI:
        	from mpi import mpi_init
        	sys.argv = mpi_init( len(sys.argv), sys.argv )
	ali3d_f(prj_stack, ref_vol, outdir, maskfile, options.r, options.snr, options.dtheta, options.maxit, options.symmetry, options.CTF, options.chunk, options.MPI)

if __name__ == "__main__":
	main()

