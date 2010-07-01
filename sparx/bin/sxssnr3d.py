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
	arglist = []
    	for arg in sys.argv: arglist.append( arg )
	progname = os.path.basename(arglist[0])
	usage = progname + " stack <output_volume> <ssnr_text_file> <reference_structure> <2Dmaskfile> --ou=outer_radius --rw=ring_width --npad=padding_times --CTF --MPI --sign=CTF_sign --sym=symmetry --random_angles=0"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--ou",       type= "int",            default=  -1,      help="  radius of particle (set to int(nx/2)-1)")
	parser.add_option("--rw",       type= "float",	        default=  1.0,     help="  ring width for calculating Fourier shell/ring correlation (set to 1)")
	parser.add_option("--npad",     type= "int",            default=  1,       help="  image padding for 3D reconstruction (set to 1)")
	parser.add_option("--CTF",      action="store_true",    default=  False,   help="  Consider CTF correction during the reconstruction (set to False)")
	parser.add_option("--sign",     type= "int",            default=  1,       help="  sign of the CTF (set to 1)")
	parser.add_option("--sym",      type= "string",         default= "c1",     help="  symmetry of the structure (set to c1)")
	parser.add_option("--MPI",      action="store_true",    default=False,     help="  use MPI version ")
	parser.add_option("--random_angles",  type= "int",      default= "0",      help="  randomize Euler angles: 0 - no, 1 - only psi, 2 - all three")
	(options, args) = parser.parse_args(arglist[1:])
	if len(args) < 1 or len(args) > 4:
    		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
	else:
		stack = args[0]
		if   len(args) == 1:
			out_vol	  = "SSNR.spi"
			ssnr_file = "ssnr"
			reference = None
			mask      = None
		elif len(args) == 2:
			out_vol   = args[1]
			ssnr_file = "ssnr"
			reference = None
			mask = None
		elif len(args) == 3:
			out_vol   = args[1]
			ssnr_file = args[2]
			reference = None
			mask = None
		elif len(args) == 4:
			out_vol   = args[1]
			ssnr_file = args[2]
			reference = args[3]
			mask = None
		elif len(args) == 5:
			out_vol   = args[1]
			ssnr_file = args[2]
			reference = args[3]
			mask      = args[4]

		if options.MPI:
			from mpi import mpi_init
			sys.argv = mpi_init(len(sys.argv), sys.argv)
			from mpi import mpi_finalize
			mpi_finalize()

		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()

		from applications import ssnr3d		
		global_def.BATCH = True
		ssnr3d(stack, out_vol, ssnr_file, mask, reference, options.ou, options.rw, options.npad, options.CTF, options.sign, options.sym, options.MPI, options.random_angles)
		global_def.BATCH = False

if __name__ == "__main__":
	main()
