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
	usage = progname + """  input_image  output_directory  --wn  --apix  --Cs  --voltage  --ac  --f_start  --f_stop  --kboot  --overlap_x  --overlap_y  --edge_x  --edge_y  --MPI  --debug
	
	Micrograph Mode - Process a set of micrographs:
	
		Specify micrograph name with wild card (*) enclosed by single quotes (') or double quotes (") (Note: sxgui.py automatically adds single quotes (')). 
		Specify output directory as an argument.
		
			mpirun -np 16 sxcter.py 'Micrographs/mic*.mrc' outdir_cter --wn=512 --apix=2.29 --Cs=2.0 --voltage=300 --ac=10.0 --MPI
			
	Stack Mode - Process a stack:
	
		Specify name of stack (without wild card "*") and output directory as arguments. 
			sxcter.py bdb:stack outdir_cter --apix=2.29 --Cs=2.0 --voltage=300 --ac=10.0
	
	"""
	parser = OptionParser(usage, version=SPARXVERSION)
	parser.add_option("--wn",         type="int",           default=512,    help="size of window to use: should be slightly larger than particle box size (default 512)")
	parser.add_option("--apix",       type="float",         default= -1,    help="pixel size in angstroms: (default 1.0)")
	parser.add_option("--Cs",         type="float",         default= 2.0,   help="microscope Cs (spherical aberration): (default 2.0)")
	parser.add_option("--voltage",    type="float",         default=300.0,  help="microscope voltage in KV: (default 300.0)")
	parser.add_option("--ac",         type="float",         default=10.0,   help="amplitude contrast in percentage: (default 10.0)")
	parser.add_option("--f_start",    type="float",         default=-1.0,   help="starting frequency in 1/A: by default determined automatically (default -1.0)")
	parser.add_option("--f_stop",     type="float",         default=-1.0,   help="stop frequency in 1/A: by default determined automatically (default -1.0)")
	parser.add_option("--kboot",      type="int",           default=16,     help="number of defocus estimates for micrograph: used for error assessment (default 16)")
	parser.add_option("--overlap_x",  type="int",           default=50,     help="overlap x in percentage: (default 50)")
	parser.add_option("--overlap_y",  type="int",           default=50,     help="overlap y in percentage: (default 50)")
	parser.add_option("--edge_x",     type="int",           default=0,      help="edge x in pixels: (default 0)")
	parser.add_option("--edge_y",     type="int",           default=0,      help="edge y in pixels: (default 0)")
	parser.add_option("--MPI",        action="store_true",  default=False,  help="use MPI version (default False)")
	parser.add_option("--debug",      action="store_true",  default=False,  help="debug info printout: (default False)")

	(options, args) = parser.parse_args(arglist[1:])
	
	if len(args) != 2:
		print "see usage " + usage
		sys.exit()
	
	# NOTE: 2015/11/27 Toshio Moriya
	# Require single quotes (') or double quotes (") for input_image so that
	# sys.argv does not automatically expand wild card and create a list of file names
	#
	input_image = args[0]
	output_directory = args[1]
	# 
	# NOTE: 2016/03/17 Toshio Moriya
	# cter_mrk() will take care of all error conditions 
	
	if options.MPI:
		from mpi import mpi_init
		sys.argv = mpi_init(len(sys.argv), sys.argv)
		
	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()
	
	from morphology import cter_mrk
	# global_def.BATCH = True
	
	cter_mrk(input_image, output_directory, options.wn, pixel_size=options.apix, \
			Cs = options.Cs, voltage=options.voltage, wgh=options.ac, \
			f_start=options.f_start, f_stop=options.f_stop, kboot=options.kboot, \
			overlap_x = options.overlap_x, overlap_y = options.overlap_y, \
			edge_x = options.edge_x, edge_y = options.edge_y, \
			MPI=options.MPI, debug_mode = options.debug)
	
	# global_def.BATCH = False
	
	if options.MPI:
		from mpi import mpi_comm_rank, MPI_COMM_WORLD, mpi_finalize
		if mpi_comm_rank(MPI_COMM_WORLD) == 0:
			print "DONE!!!"
		mpi_finalize()
	else:
		print "DONE!!!"
	
if __name__ == "__main__":
	main()
