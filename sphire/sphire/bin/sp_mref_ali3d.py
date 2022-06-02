#! /usr/bin/env python
#
# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Tapu Shaikh 2019 (tapu.shaikh@mpi-dortmund.mpg.de)
# Author: Adnan Ali 2019 (adnan.ali@mpi-dortmund.mpg.de)
# Author: Luca Lusnig 2019 (luca.lusnig@mpi-dortmund.mpg.de)
# Author: Toshio Moriya 2019 (toshio.moriya@kek.jp)
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
#
# Copyright (c) 2019 Max Planck Institute of Molecular Physiology
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
from ..libpy import sp_global_def
from ..libpy.sp_global_def import sxprint, ERROR

from ..libpy.sp_global_def import *
from optparse import OptionParser
import sys

import mpi

mpi.mpi_init( 0, [] )


def run():
	arglist = []
	i = 0
	while( i < len(sys.argv) ):
		if sys.argv[i]=='-p4pg':
			i = i+2
		elif sys.argv[i]=='-p4wd':
			i = i+2
		else:
			arglist.append( sys.argv[i] )
			i = i+1
	progname = os.path.basename(arglist[0])
	usage = progname + " stack ref_vols outdir <mask> --focus=3Dmask --ir=inner_radius --ou=outer_radius --rs=ring_step --xr=x_range --yr=y_range  --ts=translational_searching_step " +\
	" --delta=angular_step --an=angular_neighborhood --center=1 --nassign=reassignment_number --nrefine=alignment_number --maxit=max_iter --stoprnct=percentage_to_stop " + \
	" --debug --fourvar=fourier_variance --CTF --snr=1.0 --ref_a=S --sym=c1 --function=user_function --MPI --kmeans"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--focus",    type="string",       default=None,             help="3D mask for focused clustering ")
	parser.add_option("--ir",       type= "int",         default=1, 	           help="inner radius for rotational correlation > 0 (set to 1)")
	parser.add_option("--ou",       type= "int",         default="-1",	           help="outer radius for rotational correlation <nx-1 (set to the radius of the particle)")
	parser.add_option("--maxit",	type= "int",         default=5, 	           help="maximum number of iteration")
	parser.add_option("--rs",       type= "int",         default="1",	           help="step between rings in rotational correlation >0 (set to 1)" ) 
	parser.add_option("--xr",       type="string",       default="4 2 1 1 1",      help="range for translation search in x direction, search is +/-xr ")
	parser.add_option("--yr",       type="string",       default="-1",	           help="range for translation search in y direction, search is +/-yr (default = same as xr)")
	parser.add_option("--ts",       type="string",       default="0.25",           help="step size of the translation search in both directions direction, search is -xr, -xr+ts, 0, xr-ts, xr ")
	parser.add_option("--delta",    type="string",       default="10 6 4  3   2",  help="angular step of reference projections")
	parser.add_option("--an",       type="string",       default="-1",	           help="angular neighborhood for local searches")
	parser.add_option("--center",   type="int",          default=0,	               help="0 - if you do not want the volume to be centered, 1 - center the volume using cog (default=0)")
	parser.add_option("--nassign",  type="int",          default=0, 	           help="number of reassignment iterations performed for each angular step (set to 3) ")
	parser.add_option("--nrefine",  type="int",          default=1, 	           help="number of alignment iterations performed for each angular step (set to 1) ")
	parser.add_option("--CTF",      action="store_true", default=False,            help="Consider CTF correction during the alignment ")
	parser.add_option("--snr",      type="float",        default=1.0,              help="Signal-to-Noise Ratio of the data")   
	parser.add_option("--stoprnct", type="float",        default=0.0,              help="Minimum percentage of assignment change to stop the program")   
	parser.add_option("--ref_a",    type="string",       default="S",              help="method for generating the quasi-uniformly distributed projection directions (default S) ")
	parser.add_option("--sym",      type="string",       default="c1",             help="symmetry of the structure ")
	parser.add_option("--function", type="string",       default="ref_ali3dm",     help="name of the reference preparation function")
	parser.add_option("--MPI",      action="store_true", default=False,            help="Use MPI version ")
	parser.add_option("--npad",     type="int",          default= 2,               help="padding size for 3D reconstruction")
	parser.add_option("--debug",    action="store_true", default=False,            help="debug ")
	parser.add_option("--fourvar",  action="store_true", default=False,            help="compute and use fourier variance")
	parser.add_option("--kmeans",   action="store_true", default=False,            help="use kmeansmref instead of equalmref")
	parser.add_option("--kmeans2",  action="store_true", default=False,            help="use kmeansmref2 (no assignment step!) instead of equalmref")

	
	(options, args) = parser.parse_args(arglist[1:])
	if len(args) < 3 or len(args) > 4:
		sxprint( "Usage: " + usage )
		sxprint( "Please run '" + progname + " -h' for detailed options" )
		ERROR( "Invalid number of parameters used. Please see usage information above." )
		return
		
	else:

		if len(args) == 3 :
			maskfile = None
		else:
			maskfile = args[3]

		if sp_global_def.CACHE_DISABLE:
			from ..libpy.sp_utilities import disable_bdb_cache
			disable_bdb_cache()
		
		sp_global_def.BATCH = True
		if options.MPI:

			if options.kmeans:
				from ..libpy.sp_applications import Kmref_ali3d_MPI
				Kmref_ali3d_MPI(args[0], args[1], args[2], maskfile, options.focus, options.maxit, options.ir, options.ou, options.rs, \
				options.xr, options.yr, options.ts, options.delta, options.an, options.center, \
				options.nassign, options.nrefine, options.CTF, options.snr, options.ref_a, options.sym, \
				options.function,  options.npad, options.debug, options.fourvar, options.stoprnct, mpi_comm=None, log=None)

			elif options.kmeans2:
				if( options.nassign != 0):
					sxprint("  Setting nassign to zero")
					options.nassign = 0
				from ..libpy.sp_applications import Kmref2_ali3d_MPI
				Kmref2_ali3d_MPI(args[0], args[1], args[2], maskfile, options.focus, options.maxit, options.ir, options.ou, options.rs, \
				options.xr, options.yr, options.ts, options.delta, options.an, options.center, \
				options.nassign, options.nrefine, options.CTF, options.snr, options.ref_a, options.sym, \
				options.function,  options.npad, options.debug, options.fourvar, options.stoprnct, mpi_comm=None, log=None)

			else:
				from ..libpy.sp_applications import mref_ali3d_MPI
				mref_ali3d_MPI(args[0], args[1], args[2], maskfile, options.focus, options.maxit, options.ir, options.ou, options.rs, \
				options.xr, options.yr, options.ts, options.delta, options.an, options.center, \
				options.nassign, options.nrefine, options.CTF, options.snr, options.ref_a, options.sym, \
				options.function,  options.npad, options.debug, options.fourvar, options.stoprnct, mpi_comm = None, log = None)
		else:

			from ..libpy.sp_applications import mref_ali3d
			mref_ali3d(args[0], args[1], args[2], maskfile, options.focus, options.maxit, options.ir, options.ou, options.rs, 
			options.xr, options.yr, options.ts, options.delta, options.an, options.center,
			options.nassign, options.nrefine, options.CTF, options.snr, options.ref_a, options.sym,
			options.function,  options.npad, options.debug, options.fourvar, options.stoprnct)
		sp_global_def.BATCH = False
		
def main():
	sp_global_def.print_timestamp( "Start" )
	sp_global_def.write_command()
	run()
	sp_global_def.print_timestamp( "Finish" )
	mpi.mpi_finalize()

if __name__ == "__main__":
	main()
