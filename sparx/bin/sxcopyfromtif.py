#! /usr/bin/env python

#
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
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


import global_def
from global_def 	import *
from optparse 		import OptionParser
import sys
import os
def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " tifdir <micdir> --inx=tif --foc=f --ext=spi --cst=1 pixel_size=2 --sca_a=1 --sca_b=1 --step=63.5 --mag=40 --MPI"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--inx",        type = "string", default="tif", help =" input extension ")
	parser.add_option("--foc",        type = "string", default="f",   help =" film or CCD frames ")
	parser.add_option("--ext",        type = "string", default="spi", help =" extenstion of output file")
	parser.add_option("--cst",        type = "float",  default=1,     help =" contrast invert or not, -1=invert ")
	parser.add_option("--pixel_size", type = "float",  default=1,     help =" the dimension adjusted output image pixel size")
	parser.add_option("--sca_a",      type = "float",  default=1,     help =" scanner OD converting parameter a, check manual of the scanner ")
	parser.add_option("--sca_b",      type = "float",  default=1,     help =" scanner OD converting parameter b, check manual of the scanner ")
	parser.add_option("--step",       type = "float",  default=63.5,  help =" scan step size of scanner or CCD camera ")
 	parser.add_option("--mag",        type = "float",  default=40,    help =" magnification at which the images are taken ")		
	parser.add_option("--MPI", action="store_true", default=False,     help="  whether using MPI version ")
	(options, args) = parser.parse_args()    	
    	if len(args) < 1:
        	print "usage: "      + usage
        	print "Please run '" + progname + " -h' for detailed options"
	else:
	
		if len(args) == 1: 
			outdir = None
		else:
			outdir = args[1]

		from applications import copyfromtif

		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()

		if options.MPI:
			from mpi import mpi_init
			sys.argv = mpi_init(len(sys.argv),sys.argv)		

		global_def.BATCH = True

		copyfromtif(args[0], outdir, options.inx, options.foc, options.ext, options.cst, options.pixel_size, options.sca_a, options.sca_b, options.step, options.mag, options.MPI)
		global_def.BATCH = False


if __name__ == "__main__":
	        main()
