#!/usr/bin/env python
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
from global_def import *
from optparse import OptionParser
import sys
def main():
	arglist = []
	for arg in sys.argv:
		arglist.append( arg )
	progname = os.path.basename(arglist[0])
	usage = progname + " stack ref_vol outdir --dp=rise --dphi=rotation --apix=pixel_size --phistep=phi_step --zstep=z_step --fract=helicising_fraction --rmax=maximum_radius --rmin=min_radius --CTF --sym=c1 --function=user_function --maxit=max_iter --MPI"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--dp",       type="float",        default= 1.0,                help="delta z - translation in Angstroms")   
	parser.add_option("--dphi",     type="float",        default= 1.0,                help="delta phi - rotation in degrees")  
	parser.add_option("--apix",     type="float",        default= 1.84,               help="pixel size in Angstroms")
	parser.add_option("--rmin",     type="int",          default= 0,                  help="minimal radial extent of structure")   
	parser.add_option("--rmax",     type="int",          default= 70,                 help="maximal radial extent of structure")
	parser.add_option("--fract",    type="float",        default= 0.66,               help="fraction of the volume used for helical search")
	parser.add_option("--sym",      type="string",       default="c1",                help="symmetry of the structure")
	parser.add_option("--function", type="string",       default="helical",  	      help="name of the reference preparation function")
	parser.add_option("--zstep",    type="int",          default= 1,                  help="Step size for translational search along z")   
	parser.add_option("--CTF",      action="store_true", default=False,               help="CTF correction")
	parser.add_option("--maxit",    type="int",          default=5,                   help="maximum number of iterations performed")
	parser.add_option("--MPI",      action="store_true", default=False,               help="use MPI version")
	(options, args) = parser.parse_args(arglist[1:])
	if len(args) != 3:
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
	else:
		if options.MPI:
			from mpi import mpi_init, mpi_finalize
			sys.argv = mpi_init(len(sys.argv), sys.argv)
		else:
			print "There is only MPI version of sxfilrecons3d.py. See SPARX wiki page for downloading MyMPI details."
			sys.exit()
			
		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()

		from development import filrecons3D_MPI
		global_def.BATCH = True
		filrecons3D_MPI(args[0], args[1], args[2], options.dp, options.dphi, options.apix, options.function, options.zstep, options.fract, options.rmax, options.rmin,
		                options.CTF, options.maxit, options.sym)
		
		global_def.BATCH = False

		if options.MPI:  mpi_finalize()

if __name__ == "__main__":
	main()
