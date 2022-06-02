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
from   optparse import OptionParser
import sys

def run():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " volume stack  <maskfile> --delta=angular_step --method=S --phiEqpsi=Minus --symmetry=c1 --angles=angles.txt --CTF=ctf.txt --noise=s"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--delta",    type="float",   default=2,       help="angular step ")
	parser.add_option("--phiEqpsi", type="string",  default="Minus", help="if Minus, psi is set to minus phi (default), if Zero, set to zero ")
	parser.add_option("--method",   type="string",  default="S",     help="method of quasi-uniformly distributed Eulerian angles S (default) or P")
	parser.add_option("--symmetry", type="string",  default="c1",    help="symmetry group")
	parser.add_option("--angles",   type="string",  default=None,    help="List of angles (phi, theta, psi) or with shifts (phi, theta, psi, tx, ty)")
	parser.add_option("--noise",    type="float",   default=None,    help="add Gaussian noise with standard deviation s and zero mean")
	parser.add_option("--CTF",      type="string",  default=None,    help="list of CTF parameters")
	parser.add_option("--realspace",action="store_true", default=False,   help="real space projection")
	parser.add_option("--tril",     action="store_true", default=False,   help="trilinear interpolation projection")

	(options, args) = parser.parse_args()

	if(len(args) < 2 or len(args) > 3):
		sxprint("Usage: " + usage)
		sxprint("Please run \'" + progname + " -h\' for detailed options")
		ERROR( "Invalid number of parameters used. Please see usage information above." )
		return
		
	else:
		if len(args) == 2:
			mask = None
		else:
			mask = args[2]
			
		if sp_global_def.CACHE_DISABLE:
			from ..libpy.sp_utilities import disable_bdb_cache
			disable_bdb_cache()
		from ..libpy.sp_applications import project3d
		sp_global_def.BATCH = True
		project3d(args[0], args[1], mask, options.delta, options.method, options.phiEqpsi, options.symmetry, options.angles, \
		  listctfs=options.CTF, noise=options.noise, realsp=options.realspace, trillinear=options.tril)
		sp_global_def.BATCH = False

def main():
	sp_global_def.print_timestamp( "Start" )
	sp_global_def.write_command()
	run()
	sp_global_def.print_timestamp( "Finish" )

if __name__ == "__main__":
	main()
