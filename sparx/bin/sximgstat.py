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
from   global_def     import *

def main():
	import os
	import sys
	from optparse    import OptionParser

	arglist = []
	for arg in sys.argv:
		arglist.append( arg )

	progname = os.path.basename( arglist[0] )
	usage = progname + "  stack1 <stack2> <mask> --ccc --fsc file --inf --rad=r"
	parser = OptionParser(usage, version=SPARXVERSION)

	parser.add_option( "--ccc", action="store_true", default=False, help="print cross corelation coefficient" )
	parser.add_option( "--fsc", type="string",       default="",    help="calculate resolution curve" )
	parser.add_option( "--inf", action="store_true", default=False, help="print basic infomation of the img" )
	parser.add_option( "--rad", type="int",          default=-1,    help="radius of operation" )
	parser.add_option( "--hfsc", type="string",       default="",    help="calculate resolution curve for helical structures" )
	parser.add_option("--MPI",      action="store_true", default=False,         help="use MPI version")
	parser.add_option("--dp",       type="float",  default= 1.0,                help="delta z - translation in Angstroms")   
	parser.add_option("--dphi",     type="float",  default= 1.0,                help="delta phi - rotation in degrees")  
	parser.add_option("--rmin",     type="int",  default= 0,                help="minimal radial extent of structure")   
	parser.add_option("--rmax",     type="int",  default= 70,               help="maximal radial extent of structure")
	parser.add_option("--fract",    type="float",  default= 0.66,                help="fraction of the volume used for helical search")
	parser.add_option("--sym",      type="string", default="c1",               help="symmetry of the structure")
	parser.add_option("--CTF",      action="store_true", default=False,         help="CTF correction")
	parser.add_option("--fil_attr",      type="string", default="",               help="attribute with filament information")
	parser.add_option("--xysize",     type="int",  default=1,               help="xy dimensions of rectangular reconstruction for helical fsc (--hfsc option)")
	parser.add_option("--snr",      type="float",  default= 1.0,                help="Signal-to-Noise Ratio of the data")	
	parser.add_option("--npad",     type="int",    default=2,                  help="padding size for 3D reconstruction")
	(options,args) = parser.parse_args( arglist[1:] )
	
	if len(args)<1 or len(args)>3:
		print "usage: " + usage
		print "Please run '" + progname + " -h' for detailed options"
		sys.exit(-1)


	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()
	
	if len(options.hfsc) > 0:
		if not(options.MPI):
			print "The single processor version of hfsc option is not supported, please run the mpi version instead with --MPI. To install MPI (if it is not already installed), see http://sparx-em.org/sparxwiki/MPI-installation" 
			sys.exit()
		from mpi import mpi_init
		sys.argv = mpi_init(len(sys.argv), sys.argv)
		from development import imgstat_hfsc_MPI
		global_def.BATCH = True
		imgstat_hfsc_MPI( args[0], options.hfsc, options.sym, options.dp, options.dphi, options.fract, options.CTF, options.rmax, options.rmin, options.xysize, options.fil_attr, options.snr)
		global_def.BATCH = False
		
		from mpi import mpi_finalize
		mpi_finalize()
	else:
		from applications import imgstat
		global_def.BATCH = True
		imgstat( args, options.ccc, options.fsc, options.inf, options.rad )
		
if __name__=="__main__":
	main()
