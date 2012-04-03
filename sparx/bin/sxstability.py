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
	usage = progname + " stack averages --ou=ou --th_grp=th_grp --thld_err=thld_err --num_ali=num_ali --verbose"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--ou",           type="int",     default=-1,        help=" outer radius for alignment")
	parser.add_option("--thld_grp",     type="int",     default=5,         help=" mininum number of objects to consider for stability (default = 5)")
	parser.add_option("--thld_err",     type="float",   default=1.732,     help=" threshld of pixel error (default = 1.732)")
	parser.add_option("--num_ali",      type="int",     default=5,         help=" number of alignments performed for stability (default = 5)")
	parser.add_option("--verbose",      action="store_true",     default=False,       help=" whether to print individual pixel error (default = False)")
	parser.add_option("--stab_part",	action="store_true",	 default=False,	      help=" whether to output the stable particles number in file")
	(options, args) = parser.parse_args()
	if len(args) != 2:
    		print "usage: " + usage
    		print "Please run '" + progname + " -h' for detailed options"
	else:
		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()

		from mpi import mpi_init, mpi_comm_size, mpi_comm_rank, MPI_COMM_WORLD
		sys.argv = mpi_init(len(sys.argv),sys.argv)
		number_of_proc = mpi_comm_size(MPI_COMM_WORLD)
		myid = mpi_comm_rank(MPI_COMM_WORLD)

		from applications import within_group_refinement
		from pixel_error import multi_align_stability
		from utilities import write_text_file	

		global_def.BATCH = True

		data = EMData.read_images(args[0])
		averages = EMData.read_images(args[1])

		nx = data[0].get_xsize()
		ou = options.ou
		num_ali = options.num_ali
		if ou == -1: ou = nx/2-2
		from utilities import model_circle, get_params2D, set_params2D
		mask = model_circle(ou, nx, nx)

		if myid == 0:
			print "%14s %20s %20s %20s %20s"%("", "Mirror stab rate", "Pixel error", "Size of stable set", "Size of set")
		for i in xrange(len(averages)):
			if i%number_of_proc == myid:
				mem = averages[i].get_attr('members')
				mem = map(int, mem)
				if len(mem) < options.thld_grp:
					print "Average %4d: Group size too small to consider for stability."%i
				else:
					class_data = [data[im] for im in mem]
					for im in class_data: set_params2D(im, [0.0, 0.0, 0.0, 0, 1.0])
					all_ali_params = []
					for ii in xrange(num_ali):
						ali_params = []
						if options.verbose:
							ALPHA = []
							SX = []
							SY = []
							MIRROR = []
							SCALE = []
						dummy = within_group_refinement(class_data, mask, True, 1, ou, 1, [2, 1], [2, 1], [1, 0.5], 90.0, 30, 0.3, 0.2)
						for im in class_data:
							alpha, sx, sy, mirror, scale = get_params2D(im)
							ali_params.extend([alpha, sx, sy, mirror])
							if options.verbose:
								ALPHA.append(alpha)
								SX.append(sx)
								SY.append(sy)
								MIRROR.append(mirror)
								SCALE.append(scale)
						all_ali_params.append(ali_params)
						if options.verbose:
							write_text_file([ALPHA, SX, SY, MIRROR, SCALE], "ali_params_grp_%03d_run_%d"%(i, ii)) 
					stable_set, mir_stab_rate, pix_err = multi_align_stability(all_ali_params, 0.0, 10000.0, options.thld_err, options.verbose)
					print "Average %4d : %20.3f %20.3f %20d %20d"%(i, mir_stab_rate, pix_err, len(stable_set), len(mem))
					if stab_part and len(stable_set) >= options.thld_grp:
						stab_mem = [0]*len(stable_set)
						for j in xrange(len(stable_set)): stab_mem[j] = mem[stable_set[j]]
						write_text_file(mem, "stab_part_%03d.txt"%i)

		global_def.BATCH = False

		from mpi import mpi_finalize
		mpi_finalize()



if __name__ == "__main__":
	main()
