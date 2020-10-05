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
from ..libpy.sp_global_def import sxprint

from ..libpy.sp_global_def import *
from optparse import OptionParser
import sys
def run():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " volume ref_volume --discrepancy=ccc --ang_scale=angular range  --shift_scale=shift range  --mag_scale=magnification range --r=radius of a mask"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--discrepancy", type="string", default="ccc", help="  Discrepancy measure used: ccc - crosscorrelation coefficient (default), SqEuclidean - Euclidean squared) ")
	parser.add_option("--ang_scale",    type='float', default=None, help="  Correct angles will be within +/- ang_scale of initial values")
	parser.add_option("--shift_scale",  type='float', default=None, help="  Correct shifts will be within +/- shift_scale of initial values")
	parser.add_option("--mag_scale",    type='float', default=None, help="  Correct magnification will be within +/- mag_scale of initial value")
	parser.add_option("--r",            type='float', default=None, help="  Radius of a spherical mask (nx/2-2)")
	(options, args) = parser.parse_args()    	

	if sp_global_def.CACHE_DISABLE:
		from ..libpy.sp_utilities import disable_bdb_cache
		disable_bdb_cache()

	if len(args) != 2:
		sxprint( "Usage: " + usage )
		sxprint( "Please run \'" + progname + " -h\' for detailed options" )
		sp_global_def.ERROR( "Invalid number of parameters used. Please see usage information above." )
		return
		
	elif(options.ang_scale != None and options.shift_scale != None and options.mag_scale != None):
		from ..libpy.sp_applications  import ali_vol_scale
		sp_global_def.BATCH = True
		ali_vol_scale(args[0], args[1], options.ang_scale, options.shift_scale, options.mag_scale, options.r, options.discrepancy)
		sp_global_def.BATCH = False
	elif(options.ang_scale is None and options.shift_scale is None and options.mag_scale != None):
		from ..libpy.sp_applications  import ali_vol_only_scale
		sp_global_def.BATCH = True
		ali_vol_only_scale(args[0], args[1], options.mag_scale, options.r, options.discrepancy)
		sp_global_def.BATCH = False
	elif(options.ang_scale is None and options.shift_scale != None and options.mag_scale is None):
		from ..libpy.sp_applications  import ali_vol_shift
		sp_global_def.BATCH = True
		ali_vol_shift(args[0], args[1], options.shift_scale, options.r, options.discrepancy)
		sp_global_def.BATCH = False
	elif(options.ang_scale != None and options.shift_scale != None and options.mag_scale is None):
		from ..libpy.sp_applications  import ali_vol
		sp_global_def.BATCH = True
		ali_vol(args[0], args[1], options.ang_scale, options.shift_scale, options.r, options.discrepancy)
		sp_global_def.BATCH = False
	elif(options.ang_scale != None and options.shift_scale is None and options.mag_scale is None):
		from ..libpy.sp_applications  import ali_vol_rotate
		sp_global_def.BATCH = True
		ali_vol_rotate(args[0], args[1], options.ang_scale, options.r, options.discrepancy)
		sp_global_def.BATCH = False

def main():
	sp_global_def.print_timestamp( "Start" )
	sp_global_def.write_command()
	run()
	sp_global_def.print_timestamp( "Finish" )

if __name__ == "__main__":
	main()
