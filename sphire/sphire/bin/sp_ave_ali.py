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
from   optparse import OptionParser
import sys


def run():
	
	progname = os.path.basename(sys.argv[0])
	# # horatio active_refactoring Jy51i1EwmLD4tWZ9_00003_1
	# usage = progname + " stack <name_output> --ali --active --set_size=param_name_with_size --set_members=param_name_with_id"
	# # horatio active_refactoring Jy51i1EwmLD4tWZ9_00003_2	
	usage = progname + " stack <name_output> --ali --set_size=param_name_with_size --set_members=param_name_with_id"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--ali"        , action = "store_true", default=False, help="Perform average using alignment parameters")

	# # horatio active_refactoring Jy51i1EwmLD4tWZ9_00004_1	
	# parser.add_option("--active"     , action = "store_true", default=False, help="Perform average only for active images")
	
	parser.add_option("--set_size"   , type   = "string"    , default=None , help="Save number of input images to parameter with given name")
	parser.add_option("--set_members", type   = "string"    , default=None , help="Save list of id of input images to parameter named \"members\", id of input images are taken from parameter with given name")
	parser.add_option("--filament"        , action = "store_true", default=False, help="Calculate stack of averages according to filament membership")
	(options, args) = parser.parse_args()
	if len(args) < 1 or len(args) > 2:
		sxprint( "Usage: " + usage )
		sxprint( "Please run \'" + progname + " -h\' for detailed options" )
		sp_global_def.ERROR( "Invalid number of parameters used. Please see usage information above." )
		return
		
	else: 
		if len(args) == 1: 
			name_output = None
		else:
			name_output = args[1]
		
		if options.filament:
			from ..libpy.sp_development import ave_ali_filament
	
			if sp_global_def.CACHE_DISABLE:
				from ..libpy.sp_utilities import disable_bdb_cache
				disable_bdb_cache()
	
			sp_global_def.BATCH = True
			ave_ali_filament(args[0], name_output, options.ali)
			sp_global_def.BATCH = False
		
		else:
			from ..libpy.sp_applications import ave_ali
	
			if sp_global_def.CACHE_DISABLE:
				from ..libpy.sp_utilities import disable_bdb_cache
				disable_bdb_cache()
	
			sp_global_def.BATCH = True
			ave_ali(args[0], name_output, options.ali, options.set_size, options.set_members)
			sp_global_def.BATCH = False

def main():
	sp_global_def.print_timestamp( "Start" )
	sp_global_def.write_command()
	main()
	sp_global_def.print_timestamp( "Finish" )

if __name__ == "__main__":
	main()
