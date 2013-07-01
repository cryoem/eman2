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


def main():
	import os
	import sys
	from optparse import OptionParser
	from global_def import SPARXVERSION
	import global_def
        arglist = []
        for arg in sys.argv:
        	arglist.append( arg )
	progname = os.path.basename(arglist[0])
	usage = progname + " <stack> outdir1 outdir2 --indir --nameroot --nx --apix --Cs --voltage --ac --kboot --MPI"
	'''
	mpirun -np 1 sxcter.py pwrot partres --indir=. --nameroot=micrograph_PSC23_A8A_1GD_11112_135 --nx=256 --apix=2.29 --Cs=2.0 --voltage=300 --ac=10.0 --kboot=16 --MPI
	'''
	parser = OptionParser(usage,version=SPARXVERSION)
	
	parser.add_option("--indir",            	  type="string",		    default= ".",     				 help="Directory containing micrographs to be processed.")
	parser.add_option("--nameroot",          	  type="string",		    default= "",     				 help="Prefix of micrographs to be processed.")
	parser.add_option("--nx",				  	  type="int",				default=256, 					 help="Size of window to use (should be slightly larger than particle box size)")
	
	parser.add_option("--apix",               	  type="float",			 	default= -1,               	     help="pixel size in Angstroms")   
	parser.add_option("--Cs",               	  type="float",			 	default= 2.0,               	 help="Microscope Cs (spherical aberation)")
	parser.add_option("--voltage",				  type="float",				default=300.0, 					 help="Microscope voltage in KV")
	parser.add_option("--ac",					  type="float",				default=10.0, 					 help="Amplitude contrast (percentage, default=10)")
	parser.add_option("--kboot",				  type="int",				default=16, 					 help="kboot")
	parser.add_option("--MPI",               	  action="store_true",   	default=False,              	 help="use MPI version")
	parser.add_option("--debug",               	  action="store_true",   	default=False,              	 help="debug")
	parser.add_option("--overlap_x",			  type="int",				default=50, 					 help="overlap x")
	parser.add_option("--overlap_y",			  type="int",				default=50, 					 help="overlap y")
	parser.add_option("--edge_x",			  	  type="int",				default=0, 					     help="edge x")
	parser.add_option("--edge_y",			      type="int",				default=0, 					     help="edge y")
	parser.add_option("--f_start",                type="float",			 	default=-1.0,               	 help="starting frequency")   
	parser.add_option("--f_stop",                 type="float",			 	default=-1.0,               	 help="stop frequency")

	(options, args) = parser.parse_args(arglist[1:])
	
	if len(args) <2 or len(args) > 3:
		print "see usage"
		sys.exit()
	
	stack = None
	
	if len(args) == 3:
		if options.MPI:
			print "Please use single processor version if specifying a stack."
			sys.exit()
		stack = args[0]
		out1 = args[1]
		out2 = args[2]
		
	if len(args) == 2:
		out1 = args[0]
		out2 = args[1]
	
	if options.apix < 0:
		print "Enter pixel size"
		sys.exit()
	
	if options.MPI:
		from mpi import mpi_init, mpi_finalize
		sys.argv = mpi_init(len(sys.argv), sys.argv)

	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()

	from morphology import cter
	global_def.BATCH = True
	cter(stack, out1, out2, options.indir, options.nameroot, options.nx, f_start=options.f_start, f_stop=options.f_stop, voltage=300.0, Pixel_size=options.apix, Cs = options.Cs, wgh=options.ac, kboot=options.kboot, MPI=options.MPI, DEBug = options.debug, overlap_x = options.overlap_x, overlap_y = options.overlap_y, edge_x = options.edge_x, edge_y = options.edge_y, guimic=None)
	global_def.BATCH = False

	if options.MPI:
		from mpi import mpi_finalize
		mpi_finalize()

if __name__ == "__main__":
	main()
