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

	if global_def.CACHE_DISABLE:
		from utilities import disable_bdb_cache
		disable_bdb_cache()

	global_def.BATCH = True

	from numpy import array
	from statistics import k_means_stab_bbenum

	R = len(args)
	Parts = []
	mem = [0]*R
	avg = [0]*R
	for r in xrange(R):
		data = EMData.read_images(args[r])
		avg[r] = len(data)
		
		part = []
        	for k in xrange(len(data)):
                	lid = data[k].get_attr('members') 
			mem[r] += len(lid)
                        lid = array(lid, 'int32') 
                        lid.sort() 
                        part.append(lid.copy())
		Parts.append(part)

	MATCH, STB_PART, CT_s, CT_t, ST, st = k_means_stab_bbenum(Parts, T=options.T, J=50, max_branching=40, stmult=0.1, branchfunc=2)

	if options.verbose:
		print MATCH
		print STB_PART
		print CT_s
		print CT_t
		print ST
		print st
		print " "

	for i in xrange(len(STB_PART)-1, -1, -1):
		if STB_PART[i] == []: del STB_PART[i]
	for i in xrange(len(CT_s)-1, -1, -1):
		if CT_s[i] == 0: del CT_s[i]

	for i in xrange(len(MATCH)):
		assert len(STB_PART[i]) == CT_s[i]
		print "Group %3d matches Group %3d "%(MATCH[i][0], MATCH[i][1]),
		for r in xrange(2, R):
			print " Group %3d"%(MATCH[i][r]),
		print ":    group size = ",
		for r in xrange(R):
			print " %3d"%len(Parts[r][MATCH[i][r]]), 
		print "     matched size = %3d"%(CT_s[i]),
		if options.verbose:
			print "   matched group = %s"%(STB_PART[i])
		else: print ""

	print "\nNumber of averages = ",
	for r in xrange(R):
		print "%3d"%(avg[r]),
	print "\nTotal number of particles = ",
	for r in xrange(R):
		print "%3d"%(mem[r]), 
	print "     number of matched particles = %5d"%(sum(CT_s))

	global_def.BATCH = False

if __name__ == "__main__":
	main()
