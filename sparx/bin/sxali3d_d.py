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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
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
	usage = progname + " stack ref_vol outdir <maskfile> --ir=inner_radius --ou=outer_radius --rs=ring_step --xr=x_range --yr=y_range  --ts=translational_search_step  --delta=angular_step --an=angular_neighborhood --center=center_type --maxit=max_iter --CTF --snr=SNR  --ref_a=S --sym=c1 --function=user_function --debug --MPI"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--ir",       type= "int",   default= 1,                  help="inner radius for rotational correlation > 0 (set to 1)")
	parser.add_option("--ou",       type= "int",   default= -1,                 help="outer radius for rotational correlation < int(nx/2)-1 (set to the radius of the particle)")
	parser.add_option("--rs",       type= "int",   default= 1,                  help="step between rings in rotational correlation >0  (set to 1)" ) 
	parser.add_option("--xr",       type="string", default= "4 2 1 1 1",        help="range for translation search in x direction, search is +/xr")
	parser.add_option("--yr",       type="string", default= "-1",               help="range for translation search in y direction, search is +/yr (default = same as xr)")
	parser.add_option("--ts", 	type="string", default= "1 1 1 0.5 0.25",   help="step size of the translation search in both directions, search is -xr, -xr+ts, 0, xr-ts, xr, can be fractional")
	parser.add_option("--delta",    type="string", default= "10 6 4 3 2",       help="angular step of reference projections")
	parser.add_option("--an",       type="string", default= "-1",               help="angular neighborhood for local searches")
	parser.add_option("--center",   type="float",  default= -1,                 help="-1: average shift method; 0: no centering; 1: center of gravity (default=-1)")
	parser.add_option("--maxit",    type="float",  default= 5,                  help="maximum number of iterations performed for each angular step (set to 5) ")
	parser.add_option("--CTF",      action="store_true", default=False,         help="Consider CTF correction during the alignment ")
	parser.add_option("--snr",      type="float",  default= 1.0,                help="Signal-to-Noise Ratio of the data")	
	parser.add_option("--ref_a",    type="string", default= "S",                help="method for generating the quasi-uniformly distributed projection directions (default S)")
	parser.add_option("--sym",      type="string", default= "c1",               help="symmetry of the refined structure")
	parser.add_option("--function", type="string", default="ref_ali3d",         help="name of the reference preparation function (ref_ali3d)")
	parser.add_option("--MPI",      action="store_true", default=False,         help="whether to use MPI version")
	parser.add_option("--Fourvar",  action="store_true", default=False,         help="compute Fourier variance")
	parser.add_option("--debug",    action="store_true", default=False,         help="debug")
	parser.add_option("--n",        action="store_true", default=False,         help="new")
	(options, args) = parser.parse_args(arglist[1:])
	if len(args) < 3 or len(args) > 4:
    		print "usage: " + usage
    		print "Please run '" + progname + " -h' for detailed options"
	else:
		if len(args) == 3 :
			mask = None
		else:
			mask = args[3]
		if options.MPI:
			from mpi import mpi_init
			sys.argv = mpi_init(len(sys.argv), sys.argv)

		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()

		if(options.n):
			from development import ali3d_n
			global_def.BATCH = True
			ali3d_n(args[0], args[1], args[2], mask, options.ir, options.ou, options.rs, options.xr,
				options.yr, options.ts, options.delta, options.an,
				options.center, options.maxit, options.CTF, options.snr, options.ref_a, options.sym,
				options.function, options.Fourvar, options.debug, options.MPI)
			global_def.BATCH = False
			sys.exit()
		from applications import ali3d_d
		global_def.BATCH = True
		ali3d_d(args[0], args[1], args[2], mask, options.ir, options.ou, options.rs, options.xr,
			options.yr, options.ts, options.delta, options.an,
			options.center, options.maxit, options.CTF, options.snr, options.ref_a, options.sym,
			options.function, options.Fourvar, options.debug, options.MPI)
		global_def.BATCH = False

if __name__ == "__main__":
	main()
