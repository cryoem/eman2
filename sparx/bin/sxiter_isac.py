#!/usr/bin/env python

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

import os
import global_def
from   global_def import *
from   optparse import OptionParser
import sys, ConfigParser

def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " stack_file --n_run=n_run --img_per_grp=img_per_grps"

	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--ir",             type="int",          default=1,       help="inner ring of the resampling to polar coordinates ")
	parser.add_option("--ou",             type="int",          default=-1,      help="outer ring of the resampling to polar coordinates ")
	parser.add_option("--rs",             type="int",          default=1,       help="ring step of the resampling to polar coordinates ")
	parser.add_option("--xr",             type="float",        default=1.0,     help="x range of translational search ")
	parser.add_option("--yr",             type="float",        default=1.0,     help="y range of translational search ")
	parser.add_option("--ts",             type="float",        default=1.0,     help="search step of translational search ")
	parser.add_option("--init_maxit",     type="int",          default=2,       help="maximum iteration of ISAC program in the initialization phase ")
	parser.add_option("--main_maxit",     type="int",          default=3,       help="maximum iteration of ISAC program in the main phase ")
	parser.add_option("--match_loop_first",     type="int",          default=3,       help="number of iterations to run 2-way matching in the first phase ")
	parser.add_option("--match_loop_second",    type="int",          default=5,       help="number of iterations to run 2-way (or 3-way) matching in the second phase ")
	parser.add_option("--CTF",            action="store_true", default=False,   help="whether to use CTF information ")
	parser.add_option("--snr",            type="float",        default=1.0,     help="signal-to-noise ratio ")
	parser.add_option("--num_ali",        type="int",          default=5,       help="number of alignments when checking for stability ")
	parser.add_option("--loops_reali",    type="int",          default=1,       help="number of iterations in ISAC before checking stability ")
	parser.add_option("--th_err",         type="float",        default=1.0,     help="the threshold of stability ")
	parser.add_option("--max_Iter",       type="int",          default=10,      help="maximum overall iterations in the first phase ")
	parser.add_option("--n_run",          type="int",          default=4,       help="number of indepentdent runs ")
	parser.add_option("--th_grp",         type="int",          default=10,      help="the threshold of size of reproducible group ")
	parser.add_option("--img_per_grp",    type="int",          default=100,     help="number of images per group ")
	parser.add_option("--FL",             type="float",        default=0.1,     help="number of images per group ")
	parser.add_option("--FH",             type="float",        default=0.3,     help="number of images per group ")
	parser.add_option("--MPI",            action="store_true", default=True,    help="use MPI version ")
	(options, args) = parser.parse_args()

	if len(args) != 1:
    		print "usage: " + usage
    		print "Please run '" + progname + " -h' for detailed options"
		sys.exit()

	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()

	if options.MPI:
		from mpi import mpi_init
		sys.argv = mpi_init(len(sys.argv),sys.argv)

	from development import iter_isac
	global_def.BATCH = True
	iter_isac(args[0], options.ir, options.ou, options.rs, options.xr, options.yr, options.ts, options.init_maxit, options.main_maxit, \
		    options.match_loop_first, options.match_loop_second, options.CTF, options.snr, options.num_ali, options.loops_reali, \
		    options.th_err, options.max_Iter, options.n_run, options.th_grp, options.img_per_grp, options.FL, options.FH)
	global_def.BATCH = False


if __name__ == "__main__":
	main()
