#!/usr/bin/env python

#
# Author: Steven Ludtke, 02/12/2013 (sludtke@bcm.edu)
# Copyright (c) 2000-2013 Baylor College of Medicine
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

import pprint
from EMAN2 import *


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <ddd_movie_stack>

	This program will take an image stack from movie mode on a DirectElectron DDD camera and process it in various ways.
	The input stack should be <dark ref> <gain ref> <img 1> ...
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--darkgain", action="store_true",help="Performs dark and gain correction on each frame as a preprocessing step. Default = false, but should usually be performed",default=False)
	parser.add_argument("--align_frames", action="store_true",help="Perform whole-frame alignment of the stack",default=False)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	
	(options, args) = parser.parse_args()
	
	if len(args)<1:
		print usage
		parser.error("Specify input DDD stack")

	# the user may provide multiple movies to process at once
	for fsp in args:
		outname=fsp.rsplit(".",1)[0]+"_proc.hdf"		# always output to an HDF file. Output contents vary with options
		
		n=EMUtil.get_image_count(fsp)
		if n<3 : 
			print "ERROR: {} has only {} images. Min 3 required.".format(fsp,n)
			continue
		
		dark=EMData(fsp,0)
		gain=EMData(fsp,0)
		gain.sub(dark)												# dark correct the gain-reference
		gain.process_inplace("math.reciprocal",{"zero_to":1.0})		# so we can multiply by the gain reference
		gain.mult(1.0/gain["mean"])									# normalize so gain reference on average multiplies by 1.0
		
		for ii in xrange(2,n):
			im=EMData(fsp,ii)
			
			if options.darkgain:
				im.sub(dark)
				im.mult(gain)
				
				im.write_image(outname,ii-2)

if __name__ == "__main__":
	main()
