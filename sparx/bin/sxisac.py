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
	usage = progname + " stack_file --ir=ir --ou=ou --rs=rs --xr=xr --yr=yr --ts=ts --maxit=maxit --CTF --snr=snr --dst=dst --FL=FL --FH=FH --FF=FF --init_iter=init_iter --main_maxit=main_iter --iter_reali=iter_reali --match_first=match_first --max_round=max_round --match_second=match_second --stab_ali=stab_ali --thld_err=thld_err --indep_run=indep_run --thld_grp=thld_grp --img_per_grp=img_per_grp --generation=generation --MPI"

	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--ir",             type="int",          default=1,       help="inner ring of the resampling to polar coordinates ")
	parser.add_option("--ou",             type="int",          default=-1,      help="outer ring of the resampling to polar coordinates ")
	parser.add_option("--rs",             type="int",          default=1,       help="ring step of the resampling to polar coordinates ")
	parser.add_option("--xr",             type="float",        default=1.0,     help="x range of translational search ")
	parser.add_option("--yr",             type="float",        default=1.0,     help="y range of translational search ")
	parser.add_option("--ts",             type="float",        default=1.0,     help="search step of translational search ")
	parser.add_option("--maxit",          type="int",          default=30,      help="number of iterations for reference-free alignment ")
	parser.add_option("--CTF",            action="store_true", default=False,   help="whether to use CTF information (default=False, currently True is not supported)")
	parser.add_option("--snr",            type="float",        default=1.0,     help="signal-to-noise ratio (only meaningful when CTF is enabled, currently not supported)")
	parser.add_option("--dst",            type="float",        default=90.0,    help="discrete angle used in within group alignment ")
	parser.add_option("--FL",             type="float",        default=0.1,     help="lowest stopband frequency used in the tangent filter")
	parser.add_option("--FH",             type="float",        default=0.3,     help="highest stopband frequency used in the tangent filter ")
	parser.add_option("--FF",             type="float",        default=0.2,     help="fall-off of the tangent filter ")
	parser.add_option("--init_iter",      type="int",          default=3,       help="number of iterations of ISAC program in initialization ")
	parser.add_option("--main_iter",      type="int",          default=3,       help="number of iterations of ISAC program in main part ")
	parser.add_option("--iter_reali",     type="int",          default=1,       help="number of iterations in ISAC before checking stability ")
	parser.add_option("--match_first",    type="int",          default=1,       help="number of iterations to run 2-way matching in the first phase ")
	parser.add_option("--max_round",      type="int",          default=20,      help="maximum rounds of generating candidate averages in the first phase ")
	parser.add_option("--match_second",   type="int",          default=5,       help="number of iterations to run 2-way (or 3-way) matching in the second phase ")
	parser.add_option("--stab_ali",       type="int",          default=5,       help="number of alignments when checking stability ")
	parser.add_option("--thld_err",       type="float",        default=1.732,   help="the threshold of pixel error when checking stability ")
	parser.add_option("--indep_run",      type="int",          default=4,       help="number of indepentdent runs for reproducibility (default=4, currently other values not supported")
	parser.add_option("--thld_grp",       type="int",          default=10,      help="the threshold of size of reproducible class (essentially minimum size of class)")
	parser.add_option("--img_per_grp",    type="int",          default=100,     help="number of images per group in the ideal case (essentially maximum size of class)")
	parser.add_option("--generation",     type="int",          default=1,       help="the n-th approach on the dataset ")
	parser.add_option("--MPI",            action="store_true", default=True,    help="use MPI version (default=True, currently False is not supported)")
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

	from applications import iter_isac
	global_def.BATCH = True
	iter_isac(args[0], options.ir, options.ou, options.rs, options.xr, options.yr, options.ts, options.maxit, options.CTF, options.snr, \
		options.dst, options.FL, options.FH, options.FF, options.init_iter, options.main_iter, options.iter_reali, options.match_first, \
		options.max_round, options.match_second, options.stab_ali, options.thld_err, options.indep_run, options.thld_grp, \
		options.img_per_grp, options.generation)
	global_def.BATCH = False


if __name__ == "__main__":
	main()
