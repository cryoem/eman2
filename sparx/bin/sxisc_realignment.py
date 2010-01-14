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


import os, sys
import global_def
from   global_def     import *
from   user_functions import *
from   optparse       import OptionParser

def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " stack averages out_averages outdir --ou=outer_radius --xr=x_range --ts=translation_step --maxit=max_iteration --CTF --snr=SNR --function=user_function_name --Fourvar --th_err=threshold_cutoff --ncpu=number_of_cpu --nodelist=list_of_nodes --ali=kind_of_alignment"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--ou",       type="int",        default=-1,             help="outer radius for rotational correlation < nx/2-1 (set to the radius of the particle)")
	parser.add_option("--xr",       type="string",       default="4 2 1 1",      help="range for translation search in x direction, search is +/xr ")
	parser.add_option("--ts",       type="string",       default="2 1 0.5 0.25", help="step of translation search in both directions")
	parser.add_option("--maxit",    type="float",        default=0,              help="maximum number of iterations (0 means the maximum iterations is 10, but it will automatically stop should the criterion falls")
	parser.add_option("--CTF",      action="store_true", default=False,          help="Consider CTF correction during the alignment ")
	parser.add_option("--snr",      type="float",        default=1.0,            help="signal-to-noise ratio of the data (set to 1.0)")
	parser.add_option("--Fourvar",  action="store_true", default=False,          help="compute Fourier variance")
	parser.add_option("--function", type="string",       default="ref_ali2d",    help="name of the reference preparation function")
	parser.add_option('--th_err',   type='float',        default=5.0,            help='cutoff threshold to pixel error')
	parser.add_option('--ncpu',     type='int',          default=1,              help='number of cpus required in mpirun, if only one the program run without mpirun')
	parser.add_option('--nodelist', type='string',       default='local',        help='list of nodes to mpirun, if "local" run mpirun on local machine')
	parser.add_option('--ali',      type='string',       default='ali2d_c',      help='name of the function to aligned, ali2d_c or ali2d_a')
	(options, args) = parser.parse_args()
	if len(args) != 4:
    		print "usage: " + usage
    		print "Please run '" + progname + " -h' for detailed options"
	else:
		from applications import isc_realignment

		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()

		global_def.BATCH = True
		isc_realignment(args[0], args[1], args[2], args[3], options.ou, options.xr, options.ts, 
				options.maxit, options.function, options.th_err, options.snr, options.CTF,
				options.Fourvar, options.ncpu, options.nodelist, options.ali)
		global_def.BATCH = False

if __name__ == "__main__":
	main()
