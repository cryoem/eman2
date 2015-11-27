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
	from optparse import OptionParser, SUPPRESS_HELP
	from mpi import mpi_init, mpi_finalize, MPI_COMM_WORLD, mpi_comm_size, mpi_comm_rank, mpi_barrier
	
	progname = os.path.basename(sys.argv[0])
	usage = progname + """  input_image  output_directory  --wn  --apix  --Cs  --voltage  --ac  --kboot  --overlap_x  --overlap_y  --edge_x  --edge_y  --f_start  --f_stop  --debug
	
	Process a set of micrographs:
	
		Specify micrograph name with wild card (*) enclosed by single quotes (') or double quotes (") (Note: sxgui.py automatically adds single quotes (')). 
		The wild card (*) has to be in front of the extension. The extension must be 3 letter long excluding dot (.).
		Specify output directory as an argument.
		
			mpirun -np 16 sxcter.py 'Micrographs/mic*.mrc' outdir_cter --wn=512 --apix=2.29 --Cs=2.0 --voltage=300 --ac=10.0
			
	Process a stack:
	
		Specify name of stack (without wild card "*") and output directory as arguments. 
			sxcter.py bdb:stack outdir_cter --apix=2.29 --Cs=2.0 --voltage=300 --ac=10.0
	
	"""
	parser = OptionParser(usage, version=SPARXVERSION)
	parser.add_option("--wn",         type="int",           default=512,    help="size of window to use: should be slightly larger than particle box size (default 512)")
	parser.add_option("--apix",       type="float",         default= -1,    help="pixel size in angstroms: (default 1.0)")
	parser.add_option("--Cs",         type="float",         default= 2.0,   help="microscope Cs (spherical aberration): (default 2.0)")
	parser.add_option("--voltage",    type="float",         default=300.0,  help="microscope voltage in KV: (default 300.0)")
	parser.add_option("--ac",         type="float",         default=10.0,   help="amplitude contrast in percentage: (default 10.0)")
	parser.add_option("--kboot",      type="int",           default=16,     help="number of defocus estimates for micrograph: used for error assessment (default 16)")
	parser.add_option("--overlap_x",  type="int",           default=50,     help="overlap x in percentage: (default 50)")
	parser.add_option("--overlap_y",  type="int",           default=50,     help="overlap y in percentage: (default 50)")
	parser.add_option("--edge_x",     type="int",           default=0,      help="edge x in pixels: (default 0)")
	parser.add_option("--edge_y",     type="int",           default=0,      help="edge y in pixels: (default 0)")
	parser.add_option("--f_start",    type="float",         default=-1.0,   help="starting frequency in 1/A: by default determined automatically (default -1.0)")
	parser.add_option("--f_stop",     type="float",         default=-1.0,   help="stop frequency in 1/A: by default determined automatically (default -1.0)")
	parser.add_option("--debug",      action="store_true",  default=False,  help="debug info printout: (default False)")

	(options, args) = parser.parse_args(sys.argv[1:])
	
	# Initialize MPI related variables
	mpi_init(0, [])
	nproc     = mpi_comm_size(MPI_COMM_WORLD)
	myid      = mpi_comm_rank(MPI_COMM_WORLD)
	main_node = 0

	if len(args) != 2:
		if myid == main_node: print "see usage " + usage
		mpi_finalize()
		exit()
	
	if options.apix < 0:
		if myid == main_node: ERROR("Pixel size has to be specified", "sxcter", 1, myid)
		mpi_finalize()
		exit()
	
	input_image = args[0]
	# NOTE: 2015/11/27 Toshio Moriya
	# Require single quotes (') or double quotes (") for input_image so that
	# sys.argv does not automatically expand wild card and create a list of file names
	if input_image.find("*") != -1:
		# This is a micrograph file name pattern because the string contains wild card "*"
		stack = None
		MPI_support = True # Always use MPI version of cter() when it is micrographs
		micrograph_name = input_image
		indir, basename = os.path.split(input_image)
		nameroot, micsuffix = os.path.splitext(basename)
		
		if nameroot[-1] != "*":
			if myid == main_node: ERROR("input image file name for micrograph name (%s) must contain wild card * in front of the extension." % micrograph_name, "sxcter", 1, myid)
			mpi_finalize()
			exit()
		
		if micsuffix[0] != ".":
			if myid == main_node: ERROR("input image file name for micrograph name (%s) must contain extension." % micrograph_name, "sxcter", 1, myid)
			mpi_finalize()
			exit()
		
		# cter() will take care of the other error case of image image
		
		if not indir:
			# For input directory path, interpretate empty string as a current directory
			# Necessary to avoid error of os.listdir("") called by cter() in morphology.py
			indir = '.'
		
		nameroot = nameroot[:-1]
		
	else: 
		# This is a stack file name because the string does NOT contains wild card "*"
		stack = input_image
		MPI_support = False # Always use non MPI version of cter() when it is stack
		indir = "."
		nameroot = ""
		micsuffix = ""
		
		if nproc > 1: 
			if myid == main_node: ERROR("Please use '-np 1'. Stack mode supports only single processor version"  "sxcter", 1, myid)
			mpi_finalize()
			exit()
	
	output_directory = args[1]
	if os.path.exists(output_directory):
		if myid == main_node: ERROR('Output directory exists, please change the name and restart the program', "sxcter", 1, myid)
		mpi_finalize()
		exit()
	mpi_barrier(MPI_COMM_WORLD)
	
	if myid == main_node: 
		os.mkdir(output_directory)
	mpi_barrier(MPI_COMM_WORLD)
			
	out1 = "%s/pwrot" % (output_directory)
	out2 = "%s/partres" % (output_directory)
	# cter() will take care of the error case of output directory
	
	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()
	
	from morphology import cter
	global_def.BATCH = True
	
	cter(stack, out1, out2, indir, nameroot, micsuffix, options.wn, \
		f_start=options.f_start, f_stop=options.f_stop, voltage=options.voltage, Pixel_size=options.apix, \
		Cs = options.Cs, wgh=options.ac, kboot=options.kboot, MPI=MPI_support, DEBug = options.debug, \
		overlap_x = options.overlap_x, overlap_y = options.overlap_y, edge_x = options.edge_x, \
		edge_y = options.edge_y, guimic=None)
	
	global_def.BATCH = False
	
	mpi_finalize()
	
if __name__ == "__main__":
	main()
