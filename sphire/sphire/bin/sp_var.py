#!/usr/bin/env python
#
# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Tapu Shaikh 2019 (tapu.shaikh@mpi-dortmund.mpg.de)
# Author: Adnan Ali 2019 (adnan.ali@mpi-dortmund.mpg.de)
# Author: Luca Lusnig 2019 (luca.lusnig@mpi-dortmund.mpg.de)
# Author: Toshio Moriya 2019 (toshio.moriya@kek.jp)
# Author: Pawel A.Penczek 05/27/2009 (Pawel.A.Penczek@uth.tmc.edu)
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



import sys
import os
from ..libpy import sp_global_def
from ..libpy.sp_global_def import sxprint, ERROR
from ..libpy.sp_global_def import *

import mpi

mpi.mpi_init( 0, [] )


def run():
	from   optparse       import OptionParser
	progname = os.path.basename(sys.argv[0])
	usage = progname + " filelist outdir  --fl=filter_low_value --aa=filter_fall_off --radccc=radius_ccc  -repair=repairfile --pca --pcamask --pcanvec --MPI"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--fl",             type="float",        default=0.0,       help="cut-off frequency of hyperbolic tangent low-pass Fourier filter")
	parser.add_option("--aa",             type="float",        default=0.0,       help="fall-off of hyperbolic tangent low-pass Fourier filter")
	parser.add_option("--radccc",         type="int",          default=-1,        help="radius for ccc calculation")
	parser.add_option("--MPI",            action="store_true", default=False,     help="use MPI version" )
	parser.add_option("--repair",         type="string",       default="default", help="repair original bootstrap volumes: None or repair file name")
	parser.add_option("--pca",            action="store_true", default=False,     help="run pca" )
	parser.add_option("--pcamask",        type="string",       default=None,      help="mask for pca" )
	parser.add_option("--pcanvec",        type="int",          default=2,         help="number of eigvectors computed in PCA")
	parser.add_option("--n",              action="store_true", default=False,     help="new")
	parser.add_option("--scratch",        type="string",       default="./",      help="scratch directory")
	(options, args) = parser.parse_args(sys.argv[1:])

	if len(args)<2 :
		sxprint("Usage: " + usage)
		sxprint("Please run \'" + progname + " -h\' for detailed options")
		ERROR( "Invalid number of parameters used. Please see usage information above." )
		return
		
	else:
		files = args[0:-1]
		outdir = args[-1]

		if sp_global_def.CACHE_DISABLE:
			from ..libpy.sp_utilities import disable_bdb_cache
			disable_bdb_cache()

		if options.MPI:

			arglist = []
			for arg in sys.argv:
				arglist.append( arg )

			sp_global_def.BATCH = True
			
			if(options.n):
				from ..libpy.sp_development import var_mpi_new
				var_mpi_new( files[0], outdir, options.scratch, options.fl, options.aa, options.radccc, False, False, options.repair, options.pca, options.pcamask, options.pcanvec)
			else:
				from ..libpy.sp_applications import var_mpi
				var_mpi( files, outdir, options.fl, options.aa, options.radccc, options.repair, options.pca, options.pcamask, options.pcanvec)

			sp_global_def.BATCH = False
		else:
			sp_global_def.BATCH = True
			ERROR( "Please use MPI version" )
			from ..libpy.sp_applications import defvar
			defvar(  files, outdir, options.fl, options.aa, options.radccc, options.repair, options.pca, options.pcamask, options.pcanvec)
			sp_global_def.BATCH = False

def main():
	sp_global_def.print_timestamp( "Start" )
	sp_global_def.write_command()
	run()
	sp_global_def.print_timestamp( "Finish" )
	mpi.mpi_finalize()

if __name__ == "__main__":
	main()
