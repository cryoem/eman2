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
from   global_def import *


def main():
	import os
	import sys
	from optparse import OptionParser
        arglist = []
        for arg in sys.argv:
        	arglist.append( arg )
	progname = os.path.basename(arglist[0])
	usage = progname + """ stack outdir1 outdir2 --indir --nameroot --micsuffix --wn --apix --Cs --voltage --ac --kboot --MPI

	Process micrographs:
	
		Specify directory and prefix and suffix (type) of micrographs to process through --indir, --nameroot, and --micsuffix
		Specify output directories pwrot and partres as arguments.
		
			mpirun -np 16 sxcter.py pwrot partres --indir=. --nameroot=micrograph_PSC23_A8A_1GD_11112_135 --micsuffix=mrc --wn=512 --apix=2.29 --Cs=2.0 --voltage=300 --ac=10.0 --MPI
		After the program stops, it is advisable to concatenate all output files in partres directory:
		cd partres
		cat */* >>allctfs.txt

	Process stack:
	
		Specify name of stack and output directories as arguments.
			sxcter.py bdb:stack pwrot partres --apix=2.29 --Cs=2.0 --voltage=300 --ac=10.0

	
	"""
	parser = OptionParser(usage,version=SPARXVERSION)
	
	parser.add_option("--indir",            	  type="string",		    default= ".",     				 help="Directory containing micrographs to be processed.")
	parser.add_option("--nameroot",          	  type="string",		    default="",     				 help="Prefix of micrographs to be processed.")
	parser.add_option("--micsuffix",              type=str,                 default="",                      help="A string denoting micrograph type. For example 'mrc', 'hdf', 'ser' ...")
	parser.add_option("--wn",				  	  type="int",				default=512, 					 help="Size of window to use (should be slightly larger than particle box size, default 512)")
	
	parser.add_option("--apix",               	  type="float",			 	default= -1,               	     help="pixel size in Angstroms (default 1.0)")   
	parser.add_option("--Cs",               	  type="float",			 	default= 2.0,               	 help="Microscope Cs (spherical aberation, default 2.0)")
	parser.add_option("--voltage",				  type="float",				default=300.0, 					 help="Microscope voltage in KV (default 300.0)")
	parser.add_option("--ac",					  type="float",				default=10.0, 					 help="Amplitude contrast (percentage, default=10)")
	parser.add_option("--kboot",				  type="int",				default=16, 					 help="kboot (default 16)")
	parser.add_option("--MPI",               	  action="store_true",   	default=False,              	 help="use MPI version")
	parser.add_option("--debug",               	  action="store_true",   	default=False,              	 help="debug")
	parser.add_option("--overlap_x",			  type="int",				default=50, 					 help="overlap x (default 50%)")
	parser.add_option("--overlap_y",			  type="int",				default=50, 					 help="overlap y (default 50%)")
	parser.add_option("--edge_x",			  	  type="int",				default=0, 					     help="edge x (default 0)")
	parser.add_option("--edge_y",			      type="int",				default=0, 					     help="edge y (default 0)")
	parser.add_option("--f_start",                type="float",			 	default=-1.0,               	 help="starting frequency, units [1/A], (by default determined automatically)")   
	parser.add_option("--f_stop",                 type="float",			 	default=-1.0,               	 help="stop frequency, units [1/A], (by default determined automatically)")

	(options, args) = parser.parse_args(arglist[1:])
	
	if len(args) <2 or len(args) > 3:
		print "see usage " + usage
		sys.exit()
	
	stack = None
	
	if len(args) == 3:
		if options.MPI:
			ERROR("Please use single processor version if specifying a stack", "sxcter", 1)
			sys.exit()
		stack = args[0]
		out1 = args[1]
		out2 = args[2]
		
	elif len(args) == 2:
		out1 = args[0]
		out2 = args[1]
		if options.micsuffix == "" or options.nameroot == "":
			ERROR("Micrograph prefix and suffix (type) have to be specified", "sxcter", 1)
			sys.exit()
	else:
		ERROR("Incorrect number of parameters","sxcter",1)
	
	if options.apix < 0:
		ERROR("Pixel size has to be specified", "sxcter", 1)
		sys.exit()
	
	if options.MPI:
		from mpi import mpi_init, mpi_finalize
		sys.argv = mpi_init(len(sys.argv), sys.argv)

	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()

	from morphology import cter
	global_def.BATCH = True

	cter(stack, out1, out2, options.indir, options.nameroot, options.micsuffix, options.wn, \
		f_start=options.f_start, f_stop=options.f_stop, voltage=options.voltage, Pixel_size=options.apix, \
		Cs = options.Cs, wgh=options.ac, kboot=options.kboot, MPI=options.MPI, DEBug = options.debug, \
		overlap_x = options.overlap_x, overlap_y = options.overlap_y, edge_x = options.edge_x, \
		edge_y = options.edge_y, guimic=None)

	global_def.BATCH = False

	if options.MPI:
		from mpi import mpi_finalize
		mpi_finalize()

if __name__ == "__main__":
	main()
