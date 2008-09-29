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
from   global_def import *
from   optparse   import OptionParser
import sys
def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " stack outdir --delta=angular_delta --rand_seed=random_seed --trials=number_of_trials --refine --MPI --ite --given --ou --ir"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--delta",     type="float",        default=10.0,     help=" Angle step " )
	parser.add_option("--ir",        type="float",        default=-1,       help=" Inner radius of particle (set to )")
	parser.add_option("--ou",        type="float",        default=-1,       help=" Outer radius of particle < int(nx/2)-1")
	parser.add_option("--ite",       type="int",          default=10,       help=" Maximum number of iterations (10)")
	parser.add_option("--rand_seed", type="int",          default=-1,       help=" Random seed of initial" )
	parser.add_option("--trials",    type="int",          default=1,        help=" Number of trials")
	parser.add_option("--refine",    action="store_true", default=False,    help=" Angle refinement with simplex method")
	parser.add_option("--MPI",       action="store_true", default=False,    help=" MPI version")
	parser.add_option("--debug",     action="store_true", default=False,    help=" Activate debug")
	parser.add_option("--given",     action="store_true", default=False,    help=" Start with given angles in the header")
	(options, args) = parser.parse_args()
    	if len(args) != 2:
		print "usage: " + usage
        	print "Please run '" + progname + " -h' for detailed options"
	else:
		if options.trials <= 0: ERROR('Number of trials must be > 0: %d' % options.trials, 'find_struct', 1)
		from  applications  import  find_struct
		global_def.BATCH = True
		find_struct(args[0], args[1], options.delta, options.ir, options.ou, options.ite, options.rand_seed, options.trials, options.refine, options.MPI, options.debug, options.given)
		global_def.BATCH = False

if __name__ == "__main__":
	main()
