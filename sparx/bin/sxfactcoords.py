#!/usr/bin/env python
from __future__ import print_function
#
# Author: Pawel A.Penczek 05/27/2009 (Pawel.A.Penczek@uth.tmc.edu)
# Copyright (c) 2000-2019 The University of Texas - Houston Medical School
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

import global_def
from global_def import *
from optparse import OptionParser
from EMAN2_cppwrap import *

import os
import sys

      
def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " prj_stack .. average eigvol output_factcoords --rad=radius --neigvol=number_of_eigvol  --CTF"
	parser = OptionParser(usage, version=SPARXVERSION)
	parser.add_option("--rad",       type="int",    default=-1,     help="radius of mask")
	parser.add_option("--neigvol",   type="int",    default=-1,     help="number of eigvenvectors to use (default all)")
	parser.add_option("--fl",        type="float",  default=0.0,    help="cut-off frequency of hyperbolic tangent low-pass Fourier filter")
	parser.add_option("--aa",        type="float",  default=0.0,    help="fall-off of hyperbolic tangent low-pass Fourier filter")
	parser.add_option("--CTF",       action="store_true", default=False,  help="Use CTF")
	parser.add_option("--MPI",       action="store_true",           help="use MPI")

	(options, args) = parser.parse_args()

	if( len(args) < 4 ):
		print("usage: " + usage)
		print("Please run '" + progname + " -h' for details")
	else:
		stacks = args[0:-3]
		avgvol = args[-3]
		eigvol = args[-2]
		output = args[-1]
		
		if options.rad < 0:
			print("Error: mask radius is not given")
			sys.exit(-1)
		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()
		if options.MPI:
			from mpi import mpi_init
			sys.argv = mpi_init(len(sys.argv), sys.argv)

		from utilities import get_im
		global_def.BATCH = True
		if( get_im( stacks[0]).get_zsize() == 1 and get_im( eigvol).get_zsize() > 1):
			from applications import factcoords_prj
			factcoords_prj(stacks, avgvol, eigvol, output, options.rad, options.neigvol, options.fl, options.aa, options.CTF, options.MPI)
		else:
			from applications import factcoords_vol
			factcoords_vol(stacks, avgvol, eigvol, output, options.rad, options.neigvol, options.fl, options.aa, options.MPI)
		global_def.BATCH = False
		

if __name__ == "__main__":
	main()
