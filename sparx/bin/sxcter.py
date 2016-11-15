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
from global_def import *

def main():
	import os
	import sys
	from optparse import OptionParser
	arglist = []
	for arg in sys.argv:
		arglist.append( arg )
	
	progname = os.path.basename(arglist[0])
	usage = progname + """  input_image  output_directory  --wn=CTF_WINDOW_SIZE --apix=PIXEL_SIZE  --Cs=CS  --voltage=VOLATEGE  --ac=AMP_CONTRAST  --f_start=FREA_START  --f_stop=FREQ_STOP  --kboot=KBOOT  --overlap_x=OVERLAP_X  --overlap_y=OVERLAP_Y  --edge_x=EDGE_X  --edge_y=EDGE_Y  --set_ctf_header  --MPI  --stack_mode  --debug  
sxcter exists in for both MPI and non-MPI versions.

Milti-Micrograph Mode - Process a set of micrographs in a list file or in a directory:

	Specify a micrograph list file name (e.g. output of sxgui_unblur.py or sxgui_cter.py) and output directory as arguments. The file extension must be ".txt".
	
	mpirun -np 16 sxcter.py mic_list.txt outdir_cter --wn=512 --apix=2.29 --Cs=2.0 --voltage=300 --ac=10.0 --MPI
	
	Alternativel, specify micrograph name pattern with wild card [*] enclosed by single quotes ['] or double quotes ["] 
	(Note: sxgui.py automatically adds single quotes [']) and output directory as arguments. 
	BDB files can not be selected as input micrographs.
	
	mpirun -np 16 sxcter.py 'Micrographs/mic*.mrc' outdir_cter --wn=512 --apix=2.29 --Cs=2.0 --voltage=300 --ac=10.0 --MPI
	
Single-Micrograph Mode - Process a single micrograph:

	Specify micrograph name (without wild card "*") and output directory as arguments.
	BDB file can not be selected as input micrograph. Use single processor for this mode.
	
	sxcter.py Micrographs/mic0.mrc outdir_cter --wn=512 --apix=2.29 --Cs=2.0 --voltage=300 --ac=10.0
		
Stack Mode - Process a stack (Advanced Option):

	Use option --stack_mode, then specify name of stack (without wild card "*") and output directory as arguments. 
	--wn will be not used with this mode. Use single processor for this mode. 
	
	sxcter.py bdb:stack outdir_cter --apix=2.29 --Cs=2.0 --voltage=300 --ac=10.0 --stack_mode

"""
	parser = OptionParser(usage, version=SPARXVERSION)
	parser.add_option("--wn",              type="int",           default=512,    help="CTF window size [pixels]: should be slightly larger than particle box size. used only in micrograph modes. (default 512)")
	parser.add_option("--apix",            type="float",         default=-1.0,   help="pixel size [A]: (default -1.0)")
	parser.add_option("--Cs",              type="float",         default=2.0,    help="microscope spherical aberration (Cs) [mm]: (default 2.0)")
	parser.add_option("--voltage",         type="float",         default=300.0,  help="microscope voltage [kV]: (default 300.0)")
	parser.add_option("--ac",              type="float",         default=10.0,   help="amplitude contrast [%]: (default 10.0)")
	parser.add_option("--f_start",         type="float",         default=-1.0,   help="starting frequency [1/A]: by default determined automatically (default -1.0)")
	parser.add_option("--f_stop",          type="float",         default=-1.0,   help="stop frequency [1/A]: by default determined automatically (default -1.0)")
	parser.add_option("--kboot",           type="int",           default=16,     help="number of defocus estimates for micrograph: used for error assessment (default 16)")
	parser.add_option("--overlap_x",       type="int",           default=50,     help="overlap x [%]: (default 50)")
	parser.add_option("--overlap_y",       type="int",           default=50,     help="overlap y [%]: (default 50)")
	parser.add_option("--edge_x",          type="int",           default=0,      help="edge x [pixels]: (default 0)")
	parser.add_option("--edge_y",          type="int",           default=0,      help="edge y [pixels]: (default 0)")
	parser.add_option("--set_ctf_header",  action="store_true",  default=False,  help="set estimated CTF parameters to image header: used only in micrograph modes. (default False)")
	parser.add_option("--MPI",             action="store_true",  default=False,  help="use MPI version: (default False)")
	parser.add_option("--stack_mode",      action="store_true",  default=False,  help="use stack mode: also set a stack name to input image. this is advanced option. (default False)")
	parser.add_option("--debug",           action="store_true",  default=False,  help="print out debug info: (default False)")

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
	# NOTE: 2016/11/15 Toshio Moriya
	# Disabled global_def.BATCH so that the error message will be also written to standard output.
	# In addition, change the name log file
	# 
	# global_def.BATCH = True
	original_logfilename = global_def.LOGFILE
	global_def.LOGFILE = 'sxcter_' + original_logfilename + '.txt'
	
	result = cter_mrk(input_image, output_directory, options.wn, pixel_size=options.apix, \
					Cs = options.Cs, voltage=options.voltage, wgh=options.ac, \
					f_start=options.f_start, f_stop=options.f_stop, kboot=options.kboot, \
					overlap_x = options.overlap_x, overlap_y = options.overlap_y, \
					edge_x = options.edge_x, edge_y = options.edge_y, \
					set_ctf_header = options.set_ctf_header, MPI=options.MPI, \
					stack_mode = options.stack_mode, debug_mode = options.debug)
	
	# NOTE: 2016/11/15 Toshio Moriya
	# Disabled global_def.BATCH so that the error message will be also written to standard output
	# In addition, resotre the name log file
	# 
	# global_def.BATCH = False
	global_def.LOGFILE = original_logfilename
	
	if options.MPI:
		from mpi import mpi_comm_rank, MPI_COMM_WORLD, mpi_finalize
		if mpi_comm_rank(MPI_COMM_WORLD) == 0:
			if options.debug:
				print "returned value from cter_mrk() := ", result
			print "DONE!!!"
		mpi_finalize()
	else:
		if options.debug:
			print "returned value from cter_mrk() := ", result
		print "DONE!!!"
	
if __name__ == "__main__":
	main()
