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
from   global_def     import *
from   optparse       import OptionParser
import sys
def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " averages1 averages2 --th_grp"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--T",           type="int",     default=0,        help=" Threshold for matching")
	parser.add_option("--verbose",     action="store_true",     default=False,        help=" Threshold for matching")
	(options, args) = parser.parse_args()
	if len(args) != 2:
    		print "usage: " + usage
    		print "Please run '" + progname + " -h' for detailed options"
	else:
		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()

		global_def.BATCH = True

		from numpy import array
		from statistics import k_means_stab_bbenum

		data1 = EMData.read_images(args[0])
		data2 = EMData.read_images(args[1])
		
		Parts = []
	        
		part = []
		mem1 = []
        	for k in xrange(len(data1)):
                	lid = data1[k].get_attr('members') 
			mem1.extend(lid)
                        lid = array(lid, 'int32') 
                        lid.sort() 
                        part.append(lid.copy())
		Parts.append(part)

	        part = []
		mem2 = []
        	for k in xrange(len(data2)):
                	lid = data2[k].get_attr('members')
			mem2.extend(lid)
                        lid = array(lid, 'int32') 
                        lid.sort() 
                        part.append(lid.copy())
		Parts.append(part)

	        MATCH, STB_PART, CT_s, CT_t, ST, st = k_means_stab_bbenum(Parts, T=options.T, J=50, max_branching=40, stmult=0.1, branchfunc=2)

		print MATCH
		print STB_PART
		print CT_s
		print CT_t
		print ST
		print st

		for i in xrange(len(MATCH)):
			assert len(STB_PART[i]) == CT_s[i]
			print "Group %3d matchs Group %3d : group size = %3d %3d    matched size = %3d"%(MATCH[i][0], MATCH[i][1], len(Parts[0][MATCH[i][0]]), len(Parts[1][MATCH[i][1]]), CT_s[i]),
			if options.verbose:
				print "   matched group = %s"%(STB_PART[i])
			else: print ""

		print "Total number of particles = %5d %5d     number of matched particles = %5d"%(len(mem1), len(mem2), sum(CT_s))

		global_def.BATCH = False



if __name__ == "__main__":
	main()
