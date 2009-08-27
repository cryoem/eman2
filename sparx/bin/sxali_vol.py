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
from global_def import *
from optparse import OptionParser
import sys
def main():
	progname = os.path.basename(sys.argv[0])
	usage = progname + " volume ref_volume --discrepancy=ccc --ang_scale=angular range  --shift_scale=shift range  --mag_scale=magnification range --r=radius of a mask"
	parser = OptionParser(usage,version=SPARXVERSION)
	parser.add_option("--discrepancy", type="string", default="ccc", help="  Discrepancy measure used: ccc - crosscorrelation coefficient (default), SqEuclidean - Euclidean squared) ")
	parser.add_option("--ang_scale",    type='float', default=None, help="  Correct angles will be within +/- ang_scale of initial values")
	parser.add_option("--shift_scale",  type='float', default=None, help="  Correct shifts will be within +/- shift_scale of initial values")
	parser.add_option("--mag_scale",    type='float', default=None, help="  Correct magnification will be within +/- mag_scale of initial value")
	parser.add_option("--r",            type='float', default=None, help="  Radius of a spherical mask (nx/2-2)")
	(options, args) = parser.parse_args()    	
    	if len(args) != 2:
		print "usage: " + usage
        	print "Please run '" + progname + " -h' for detailed options"
		exit(1)
	elif(options.ang_scale != None and options.shift_scale != None and options.mag_scale != None):
		from applications  import ali_vol_scale
		global_def.BATCH = True
		ali_vol_scale(args[0], args[1], options.ang_scale, options.shift_scale, options.mag_scale, options.r, options.discrepancy)
		global_def.BATCH = False
	elif(options.ang_scale is None and options.shift_scale is None and options.mag_scale != None):
		from applications  import ali_vol_only_scale
		global_def.BATCH = True
		ali_vol_only_scale(args[0], args[1], options.mag_scale, options.r, options.discrepancy)
		global_def.BATCH = False
	elif(options.ang_scale is None and options.shift_scale != None and options.mag_scale is None):
		from applications  import ali_vol_shift
		global_def.BATCH = True
		ali_vol_shift(args[0], args[1], options.shift_scale, options.r, options.discrepancy)
		global_def.BATCH = False
	elif(options.ang_scale != None and options.shift_scale != None and options.mag_scale is None):
		from applications  import ali_vol
		global_def.BATCH = True
		ali_vol(args[0], args[1], options.ang_scale, options.shift_scale, options.r, options.discrepancy)
		global_def.BATCH = False
	elif(options.ang_scale != None and options.shift_scale is None and options.mag_scale is None):
		from applications  import ali_vol_rotate
		global_def.BATCH = True
		ali_vol_rotate(args[0], args[1], options.ang_scale, options.r, options.discrepancy)
		global_def.BATCH = False
		
if __name__ == "__main__":
		main()
