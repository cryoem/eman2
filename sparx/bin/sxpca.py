#!/usr/bin/env python
#
# Author: Pawel A.Penczek and Edward H. Egelman 05/27/2009 (Pawel.A.Penczek@uth.tmc.edu)
# Copyright (c) 2000-2006 The University of Texas - Houston Medical School
# Copyright (c) 2008-Forever The University of Virginia
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
	usage = progname + " input_stack1 ... output_stack --sdir --usebuf --MPI --shuffle --subavg=average_image --rad=mask_radius --nvec=number_of_eigenvectors --mask=maskfile"
	parser = OptionParser(usage, version=SPARXVERSION)
	parser.add_option("--subavg",  type="string",       default="",    help="subtract average")
	parser.add_option("--rad",     type="int",          default=-1,    help="radius of mask")
	parser.add_option("--nvec",    type="int",          default=1,     help="number of eigenvectors")
	parser.add_option("--mask",    type="string",       default="",    help="mask file" )
	parser.add_option("--sdir",    type="string",       default=".",   help="scratch directory")
	parser.add_option("--usebuf",  action="store_true", default=False, help="use existing buffer")
	parser.add_option("--MPI",     action="store_true", default=False, help="run mpi version" )
	parser.add_option("--shuffle", action="store_true", default=False, help="use shuffle")

	(options, args) = parser.parse_args()

	input_stacks = args[0:-1]
	output_stack = args[-1]

	if options.nvec is None:
		print "Error: number of components is not given"
		sys.exit(-2) 

	if options.MPI:
		from mpi import mpi_init
		sys.argv = mpi_init( len(sys.argv), sys.argv )

	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()
	from applications import pca
	global_def.BATCH = True
	pca(input_stacks, output_stack, options.subavg, options.rad, options.sdir, options.nvec, options.shuffle, not(options.usebuf), options.mask, options.MPI)
	global_def.BATCH = False
        if options.MPI:
		from mpi import mpi_finalize
		mpi_finalize()


if __name__ == "__main__":
	main()
