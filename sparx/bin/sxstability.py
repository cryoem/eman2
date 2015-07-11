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
	usage = progname + " stack output_average --radius=particle_radius --xr=xr --yr=yr --ts=ts --thld_err=thld_err --num_ali=num_ali --fl=fl --aa=aa --CTF --verbose --stables"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--radius",       type="int",              default=-1,          help=" particle radius for alignment")
	parser.add_option("--xr",           type="string"      ,     default="2 1",       help="range for translation search in x direction, search is +/xr (default 2,1)")
	parser.add_option("--yr",           type="string"      ,     default="-1",        help="range for translation search in y direction, search is +/yr (default = same as xr)")
	parser.add_option("--ts",           type="string"      ,     default="1 0.5",     help="step size of the translation search in both directions, search is -xr, -xr+ts, 0, xr-ts, xr, can be fractional (default: 1,0.5)")
	parser.add_option("--thld_err",     type="float",            default=0.75,        help="threshld of pixel error (default = 0.75)")
	parser.add_option("--num_ali",      type="int",              default=5,           help="number of alignments performed for stability (default = 5)")
	parser.add_option("--maxit",        type="int",              default=30,          help="number of iterations for each xr (default = 30)")
	parser.add_option("--fl",           type="float"       ,     default=0.3,         help="cut-off frequency of hyperbolic tangent low-pass Fourier filter (default = 0.3)")
	parser.add_option("--aa",           type="float"       ,     default=0.2,         help="fall-off of hyperbolic tangent low-pass Fourier filter (default = 0.2)")
	parser.add_option("--CTF",          action="store_true",     default=False,       help="Use CTF correction during the alignment ")
	parser.add_option("--verbose",      action="store_true",     default=False,       help="print individual pixel error (default = False)")
	parser.add_option("--stables",		action="store_true",	 default=False,	      help="output the stable particles number in file (default = False)")
	parser.add_option("--method",		type="string"      ,	 default=" ",	      help="SHC (standard method is default when flag is ommitted)")
	(options, args) = parser.parse_args()
	if len(args) != 1 and len(args) != 2:
    		print "usage: " + usage
    		print "Please run '" + progname + " -h' for detailed options"
	else:
		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()

		from applications   import within_group_refinement, ali2d_ras
		from pixel_error    import multi_align_stability
		from utilities      import write_text_file, write_text_row

		global_def.BATCH = True

		xrng        = get_input_from_string(options.xr)
		if  options.yr == "-1":  yrng = xrng
		else          :  yrng = get_input_from_string(options.yr)
		step        = get_input_from_string(options.ts)

		class_data = EMData.read_images(args[0])

		nx = class_data[0].get_xsize()
		ou = options.radius
		num_ali = options.num_ali
		if ou == -1: ou = nx/2-2
		from utilities import model_circle, get_params2D, set_params2D
		mask = model_circle(ou, nx, nx)

		if options.CTF :
			from filter import filt_ctf
			for im in xrange(len(class_data)):
				#  Flip phases
				class_data[im] = filt_ctf(class_data[im], class_data[im].get_attr("ctf"), binary=1)
		for im in class_data:
			im.set_attr("previousmax", -1.0e10)
			try:
				t = im.get_attr("xform.align2d") # if they are there, no need to set them!
			except:
				try:
					t = im.get_attr("xform.projection")
					d = t.get_params("spider")
					set_params2D(im, [0.0, -d["tx"], -d["ty"], 0, 1.0])
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
			if( xrng[0] == 0.0 and yrng[0] == 0.0 ):
				avet = ali2d_ras(class_data, randomize = True, ir = 1, ou = ou, rs = 1, step = 1.0, dst = 90.0, \
						maxit = options.maxit, check_mirror = True, FH=options.fl, FF=options.aa)
			else:
				avet = within_group_refinement(class_data, mask, True, 1, ou, 1, xrng, yrng, step, 90.0, \
						maxit = options.maxit, FH=options.fl, FF=options.aa, method = options.method)
				from utilities import info
				#print "  avet  ",info(avet)
			for im in class_data:
				alpha, sx, sy, mirror, scale = get_params2D(im)
				ali_params.extend([alpha, sx, sy, mirror])
				if options.verbose:
					ALPHA.append(alpha)
					SX.append(sx)
					SY.append(sy)
					MIRROR.append(mirror)
			all_ali_params.append(ali_params)
			if options.verbose:
				write_text_file([ALPHA, SX, SY, MIRROR], "ali_params_run_%d"%ii)
		"""
		avet = class_data[0]
		from utilities import read_text_file
		all_ali_params = []
		for ii in xrange(5):
			temp = read_text_file( "ali_params_run_%d"%ii,-1)
			uuu = []
			for k in xrange(len(temp[0])):
				uuu.extend([temp[0][k],temp[1][k],temp[2][k],temp[3][k]])
			all_ali_params.append(uuu)


		"""

		stable_set, mir_stab_rate, pix_err = multi_align_stability(all_ali_params, 0.0, 10000.0, options.thld_err, options.verbose, 2*ou+1)
		print "%4s %20s %20s %20s %30s %6.2f"%("", "Size of set", "Size of stable set", "Mirror stab rate", "Pixel error prior to pruning the set above threshold of",options.thld_err)
		print "Average stat: %10d %20d %20.2f   %15.2f"%( len(class_data), len(stable_set), mir_stab_rate, pix_err)
		if( len(stable_set) > 0):
			if options.stables:
				stab_mem = [[0,0.0,0] for j in xrange(len(stable_set))]
				for j in xrange(len(stable_set)): stab_mem[j] = [int(stable_set[j][1]), stable_set[j][0], j]
				write_text_row(stab_mem, "stable_particles.txt")

			stable_set_id = []
			particle_pixerr = []
			for s in stable_set:
				stable_set_id.append(s[1])
				particle_pixerr.append(s[0])
			from fundamentals import rot_shift2D
			avet.to_zero()
			l = -1
			print "average parameters:  angle, x-shift, y-shift, mirror"
			for j in stable_set_id:
				l += 1
				print " %4d  %4d  %12.2f %12.2f %12.2f        %1d"%(l,j, stable_set[l][2][0], stable_set[l][2][1], stable_set[l][2][2], int(stable_set[l][2][3]))
				avet += rot_shift2D(class_data[j], stable_set[l][2][0], stable_set[l][2][1], stable_set[l][2][2], stable_set[l][2][3] )
			avet /= (l+1)
			avet.set_attr('members', stable_set_id)
			avet.set_attr('pix_err', pix_err)
			avet.set_attr('pixerr', particle_pixerr)
			avet.write_image(args[1])



		global_def.BATCH = False

if __name__ == "__main__":
	main()
