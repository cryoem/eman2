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


import os
import global_def
from global_def import *
from optparse import OptionParser
import sys
def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " stack outdir <maskfile> --ou=outer_radius --br=brackets --center=center_type --eps=epsilon --maxit=max_iter --CTF --snr=SNR --function=user_function_name"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--ou", type="float", default=-1, help="  outer radius for a 2-D mask within which the alignment is performed")
	parser.add_option("--br", type="float", default=1.75, help="  brackets for the search of orientation parameters (each parameter will be checked +/- bracket (set to 1.75)")
	parser.add_option("--center", type="float", default=1, help="  0 - if you do not want the average to be centered, 1 - center the average (default=1)")
	parser.add_option("--eps", type="float", default=0.001, help="  stopping criterion, program will terminate when the relative increase of the criterion is less than epsilon  ")
	parser.add_option("--maxit", type="float", default=10, help="  maximum number of iterations (set to 10) ")
	parser.add_option("--CTF", action="store_true", default=False, help="  Consider CTF correction during the alignment ")
	parser.add_option("--snr", type="float", default=1.0, help="  Signal-to-Noise Ratio of the data")	
	parser.add_option("--function", type="string", default="ref_ali2d", help="  name of the reference preparation function")
	(options, args) = parser.parse_args()
	if len(args) < 2 or len(args) >3:
    		print "usage: " + usage
    		print "Please run '" + progname + " -h' for detailed options"
	else:
		if len(args) == 2:
			mask = None
		else:
			mask = args[2]
		
		from applications import local_ali2d

		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()
		
		global_def.BATCH = True	
		local_ali2d(args[0], args[1], mask, options.ou, options.br, options.center, options.eps, options.maxit, options.CTF, options.snr, options.function)
		global_def.BATCH = False	
if __name__ == "__main__":
	main()
