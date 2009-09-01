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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  2111-1307 USA
#
#

import os
import global_def
from global_def import *
from applications import autowin
from optparse import OptionParser
import sys
def main():
	arglist = []
	for arg in sys.argv:
		arglist.append( arg )	
	progname = os.path.basename(arglist[0])
	usage = progname + " mics ptl noisedoc noisemic template --deci=decimation --pck=1 --p_size=image size --sigma=1 --hf_p=half size --n_tem= number of templates --n_ptl=number of particles --cst=1 --prm=micrograph --MPI"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--deci",   type = "float",   default= 1,               help="  decimation ratio of the micrograph in case of particle is too large ")
	parser.add_option("--pck" ,   type = "float",   default= 1,               help="  method used in particle detection pck=1 Gaussian high pass filter pck=2 fast local normalization")
	parser.add_option("--p_size", type = "float",   default= 128,             help="  window size")
	parser.add_option("--sigma",  type = "float",   default= 1,               help="  sigma of Gaussian high pass filter, set as 1, or 2")
	parser.add_option("--hf_p",   type = "float",   default= 64,              help="  approximate particle diameter")
	parser.add_option("--n_ptl",  type = "float",   default= 1000,            help="  number of particles to be picked")		
	parser.add_option("--cst",    type = "float",   default= 1,               help="  invert image contrast (=-1) or not")
	parser.add_option("--CTF", action="store_true", default=False,            help="  if true, ctf paramters would be copied from micrograph headers into particle image headers")
	parser.add_option("--prm",    type = "string",  default= "micrograph",    help="  invert image contrast (=-1) or not")
	parser.add_option("--MPI", action="store_true", default=False,     help="  whether using MPI version ")
	(options, args) = parser.parse_args(arglist[1:]) 
    	if len(args)  != 5:
        	print "usage: " + usage
        	print "Please run '" + progname + " -h' for detailed options"
	else: 
		from applications import autowin

		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()

		if options.MPI: 
			from mpi import mpi_init
			sys.argv = mpi_init(len(sys.argv), sys.argv)
		global_def.BATCH = True
		autowin(args[0], args[1], args[2], args[3], args[4], options.deci, options.pck, options.p_size, options.sigma, options.hf_p, options.n_ptl, options.cst, options.MPI)
		global_def.BATCH = False

if __name__ == "__main__":
	        main()
