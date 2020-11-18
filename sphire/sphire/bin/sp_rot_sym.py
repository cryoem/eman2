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

def run():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " input_filename output_filename --sym=Symmetry group --phi --theta --psi=The 3 Eulerian angles in degrees --r=Radius of mask --phirange --thetarange --psirange=A search scale for each angle --ftol --xtol = convergence criterion the function and angles values"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--sym",        type="string", default="c1", help="  String that specifies the point group symmetry. default = 'c1'")
	parser.add_option("--phi",        type='float',default=0.0, help="  phi angle, default = 0")
	parser.add_option("--theta",      type='float',default=0.0, help=" theta angle, default=0")
	parser.add_option("--psi",        type='float',default=0.0, help=" phi angle, default=0")
	parser.add_option("--r",          type='float',default=None,help=" Input the radius of the mask. default=None")
	parser.add_option("--phirange",   type='float',default=20.0,help=" The search scale for phi angle...default=20")
	parser.add_option("--thetarange", type='float',default=20.0,help=" The search scale for theta angle...default=20")
	parser.add_option("--psirange",   type='float',default=20.0,help=" The search scale for psi angle...default=20")
	parser.add_option("--ftol",       type='float',default=1.e-4,help=" convergence criterion on the function values...default = 1.e-4")
	parser.add_option("--xtol",       type='float',default=1.e-4,help=" convergence criterion on the variable values...default = 1.e-4")
	
	(options, args) = parser.parse_args()
	
	if len(args) != 2:
		sxprint("Usage: " + usage)
		sxprint("Please run \'" + progname + " -h\' for detailed options")
		ERROR( "Invalid number of parameters used. Please see usage information above." )
		return
	
	else:
		if sp_global_def.CACHE_DISABLE:
			from ..libpy.sp_utilities import disable_bdb_cache
			disable_bdb_cache()
		from ..libpy.sp_applications  import  rot_sym
		sp_global_def.BATCH = True
		rot_sym(args[0],args[1],options.sym,options.r,options.phi,options.theta,options.psi,options.phirange,options.thetarange,options.psirange,options.ftol,options.xtol)
		sp_global_def.write_command('.')
		sp_global_def.BATCH = False

def main():
	sp_global_def.print_timestamp( "Start" )
	run()
	sp_global_def.print_timestamp( "Finish" )

if __name__ == "__main__":
	main()
