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
from   user_functions import *
from   optparse       import OptionParser
import sys

def write_montage_file(stack, montage_file, N, gx, gy, bg, scale):

	from utilities import model_blank

	if gy == -1: gy = gx

	data = EMData.read_images(stack)
	if scale:
		for im in data:
			st = Util.infomask(im, None, True)
			im -= st[0]
			im /= st[1]
			
	nx = data[0].get_xsize()
	
	K = len(data)
	M = (K-1)/N+1
	
	NX = (nx+gx)*N
	NY = (nx+gy)*M
	
	maxn = -1.e-20
	minn = 1.e+20
	avgn = 0
	for im in data:
		st = Util.infomask(im, None, True)
		avgn += st[0]
		if st[3] > maxn: maxn = st[3]
		if st[2] < minn: minn = st[2]
	avgn /= K		
	
	if bg == 0:
		bk = minn
	elif bg == 1:
		bk = maxn
	elif bg == 2:
		bk = 0
	else:
		bk = avgn
	
	montage = model_blank(NX, NY, 1, bk)
	
	for i in xrange(K):
		col = i%N
		row = M-1-i/N
		for s in xrange(nx):
			for t in xrange(nx):
				v = data[i].get_value_at(s, t)
				montage.set_value_at(col*(nx+gx)+s, row*(nx+gy)+t, v)
	montage.write_image(montage_file)

def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " stack montagefile --N=image_per_row --gx=gap_on_x_axis --gy=gap_on_y_axis --bg=background_option --scale"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--N",       type="int",  default=10,           help="number of images per row (default = 10)")
	parser.add_option("--gx",      type="int",  default=0,            help="number of pixels of the gap on x-axis (default = 0)")
	parser.add_option("--gy",      type="int",  default=-1,           help="number of pixels of the gap on y-axis (default = gx)")
	parser.add_option("--bg",      type="int",  default=0,            help="background option: 0. use minimum as background; 1. use maximum as background; 2. use 0 as background; 3. use average as background (default = 0)") 
	parser.add_option("--scale",   action="store_true",  default=False,   help="whether to scale all particles to (0, 1) ")
	(options, args) = parser.parse_args()
	if len(args) != 2:
    		print "usage: " + usage
    		print "Please run '" + progname + " -h' for detailed options"
	else:
		if global_def.CACHE_DISABLE:
			from utilities import disable_bdb_cache
			disable_bdb_cache()
		
		global_def.BATCH = True
		write_montage_file(args[0], args[1], options.N, options.gx, options.gy, options.bg, options.scale)
		global_def.BATCH = False


if __name__ == "__main__":
	main()
