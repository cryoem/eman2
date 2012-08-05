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
	from utilities import get_input_from_string
	progname = os.path.basename(sys.argv[0])
	usage = progname + " stack <averages> --ou=ou --xr=xr --yr=yr --ts=ts --thld_grp=thld_grp --thld_err=thld_err --num_ali=num_ali --fl=fl --aa=aa --CTF --verbose --stab_part"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--ou",           type="int",              default=-1,          help=" outer radius for alignment")
	parser.add_option("--xr",           type="string"      ,     default="2 1",       help="range for translation search in x direction, search is +/xr")
	parser.add_option("--yr",           type="string"      ,     default="-1",        help="range for translation search in y direction, search is +/yr (default = same as xr)")
	parser.add_option("--ts",           type="string"      ,     default="1 0.5",     help="step size of the translation search in both directions, search is -xr, -xr+ts, 0, xr-ts, xr, can be fractional")
	parser.add_option("--thld_grp",     type="int",              default=5,           help="mininum number of objects to consider for stability (default = 5)")
	parser.add_option("--thld_err",     type="float",            default=1.732,       help="threshld of pixel error (default = 1.732)")
	parser.add_option("--num_ali",      type="int",              default=5,           help="number of alignments performed for stability (default = 5)")
	parser.add_option("--fl",           type="float"       ,     default=0.3,         help="cut-off frequency of hyperbolic tangent low-pass Fourier filter")
	parser.add_option("--aa",           type="float"       ,     default=0.2,         help="fall-off of hyperbolic tangent low-pass Fourier filter")
	parser.add_option("--CTF",          action="store_true",     default=False,       help="Consider CTF correction during the alignment ")
	parser.add_option("--verbose",      action="store_true",     default=False,       help="print individual pixel error (default = False)")
	parser.add_option("--stab_part",	action="store_true",	 default=False,	      help="output the stable particles number in file")
	(options, args) = parser.parse_args()
	if len(args) != 1 and len(args) != 2:
    		print "usage: " + usage
    		print "Please run '" + progname + " -h' for detailed options"
	else:
		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()

		from applications import within_group_refinement
		from pixel_error import multi_align_stability_new
		from utilities import write_text_file, write_text_row

		global_def.BATCH = True

		xrng        = get_input_from_string(options.xr)
		if  options.yr == "-1":  yrng = xrng
		else          :  yrng = get_input_from_string(options.yr)
		step        = get_input_from_string(options.ts)

		data = EMData.read_images(args[0])
		if len(args) > 1:
			averages = EMData.read_images(args[1])
		else:
			averages = [ EMData() ]
			averages[0].set_attr( "members", range(len(data)) )

		nx = data[0].get_xsize()
		ou = options.ou
		num_ali = options.num_ali
		if ou == -1: ou = nx/2-2
		from utilities import model_circle, get_params2D, set_params2D
		mask = model_circle(ou, nx, nx)

		print "%14s %20s %20s %20s %20s"%("", "Mirror stab rate", "Pixel error", "Size of stable set", "Size of set")
		for i in xrange(len(averages)):
			mem = averages[i].get_attr('members')
			mem = map(int, mem)
			if len(mem) < options.thld_grp:
				print "Average %4d: Group size too small to consider for stability."%i
			else:
				class_data = [data[im] for im in mem]
				if options.CTF :
					from filter import filt_ctf
					for im in xrange(len(class_data)):
						class_data[im] = filt_ctf(class_data[im], class_data[im].get_attr("ctf"), binary=1)
				for im in class_data:
					try:
						t = im.get_attr("xform.align2d") # if they are there, no need to set them!
					except:
						try:
							t = im.get_attr("xform.projection")
							d = t.get_params("spider")
							set_params2D(im, [0.0,-d["tx"],-d["ty"],0,1.0])
						except:
							set_params2D(im, [0.0, 0.0, 0.0, 0, 1.0])
				all_ali_params = []
				for ii in xrange(num_ali):
					ali_params = []
					if options.verbose:
						ALPHA = []
						SX = []
						SY = []
						MIRROR = []
						SCALE = []
					dummy = within_group_refinement(class_data, mask, True, 1, ou, 1, xrng, yrng, step, 90.0, 30, options.fl, options.aa)
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
				stable_set, mir_stab_rate, pix_err = multi_align_stability_new(all_ali_params, 0.0, 10000.0, options.thld_err, options.verbose, 2*ou+1)
				print "Average %4d : %20.3f %20.3f %20d %20d"%(i, mir_stab_rate, pix_err, len(stable_set), len(mem))
				if options.stab_part and len(stable_set) >= options.thld_grp:
					stab_mem = [[0,0.0,0] for j in xrange(len(stable_set))]
					for j in xrange(len(stable_set)): stab_mem[j] = [mem[int(stable_set[j][1])], stable_set[j][0], j]
					write_text_row(stab_mem, "stab_part_%06d.txt"%i)

		global_def.BATCH = False

if __name__ == "__main__":
	main()
