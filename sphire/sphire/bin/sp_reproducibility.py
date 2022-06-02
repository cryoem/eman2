#!/usr/bin/env python
#
# Author: Markus Stabrin 2019 (markus.stabrin@mpi-dortmund.mpg.de)
# Author: Fabian Schoenfeld 2019 (fabian.schoenfeld@mpi-dortmund.mpg.de)
# Author: Thorsten Wagner 2019 (thorsten.wagner@mpi-dortmund.mpg.de)
# Author: Tapu Shaikh 2019 (tapu.shaikh@mpi-dortmund.mpg.de)
# Author: Adnan Ali 2019 (adnan.ali@mpi-dortmund.mpg.de)
# Author: Luca Lusnig 2019 (luca.lusnig@mpi-dortmund.mpg.de)
# Author: Toshio Moriya 2019 (toshio.moriya@kek.jp)
# Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
#
# Copyright (c) 2019 Max Planck Institute of Molecular Physiology
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
from builtins import range
import os
from ..libpy import sp_global_def
from ..libpy.sp_global_def import sxprint
from ..libpy.sp_global_def import *
from   optparse   import OptionParser
import sys

def run():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " averages1 averages2 --th_grp"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option( "--T",             type="int",          default=0,     help=" Threshold for matching" )
	parser.add_option( "--J",             type="int",          default=50,    help=" J" )
	parser.add_option( "--max_branching", type="int",     	  default=40,     help=" maximum branching" )
	parser.add_option( "--verbose",       action="store_true", default=False, help=" Threshold for matching" )
	parser.add_option( "--timing",        action="store_true", default=False, help=" Get the timing" )

	(options, args) = parser.parse_args()

	if sp_global_def.CACHE_DISABLE:
		from ..libpy.sp_utilities import disable_bdb_cache
		disable_bdb_cache()

	sp_global_def.BATCH = True

	from numpy import array
	from ..libpy.sp_statistics import k_means_stab_bbenum

	R = len(args)
	Parts = []
	mem = [0]*R
	avg = [0]*R
	for r in range(R):
		data = EMData.read_images(args[r])
		avg[r] = len(data)
		
		part = []
		for k in range(len(data)):
			lid = data[k].get_attr('members') 
			mem[r] += len(lid)
			lid = array(lid, 'int32') 
			lid.sort() 
			part.append(lid.copy())
		Parts.append(part)

	if options.timing:
		from time import time
		time1 = time()

	MATCH, STB_PART, CT_s, CT_t, ST, st = k_means_stab_bbenum(Parts, T=options.T, J=options.J, max_branching=options.max_branching, stmult=0.1, branchfunc=2)

	if options.verbose:
		sxprint(MATCH)
		sxprint(STB_PART)
		sxprint(CT_s)
		sxprint(CT_t)
		sxprint(ST)
		sxprint(st)
		sxprint(" ")

	for i in range(len(MATCH)):
		u = MATCH[i][0]  # u is the group in question in partition 1
		assert len(STB_PART[u]) == CT_s[u]
		sxprint("Group ", end=' ') 
		for r in range(R):
			sxprint("%3d "%(MATCH[i][r]), end=' ')
		sxprint(" matches:   group size = ", end=' ')
		for r in range(R):
			sxprint(" %3d"%len(Parts[r][MATCH[i][r]]), end=' ') 
		sxprint("     matched size = %3d"%(CT_s[u]), end=' ')
		if options.verbose:
			sxprint("   matched group = %s"%(STB_PART[u]))
		else: 
			sxprint("")

	sxprint("\nNumber of averages = ", end=' ')
	for r in range(R):
		sxprint("%3d"%(avg[r]), end=' ')
	sxprint("\nTotal number of particles = ", end=' ')
	for r in range(R):
		sxprint("%3d"%(mem[r]), end=' ') 
	sxprint("     number of matched particles = %5d"%(sum(CT_s)))

	if options.timing:
		sxprint("Elapsed time = ", time() - time1)

	sp_global_def.BATCH = False

def main():
	sp_global_def.print_timestamp( "Start" )
	sp_global_def.write_command()
	run()
	sp_global_def.print_timestamp( "Finish" )

if __name__ == "__main__":
	main()

