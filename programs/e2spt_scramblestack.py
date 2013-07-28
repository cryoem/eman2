#!/usr/bin/env python

#
# Author: Jesus Galaz  20/3/2012
# Copyright (c) 2011- Baylor College of Medicine
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA

from sys import argv
import os
from EMAN2 import *
import random

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <Volume file>

	WARNING: This program is still under development.
	
	It takes in an HDF stack and randomizes the position of images in it."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--input", type=str, help="The name of the volumes stack that HAVE BEEN ALIGNED to a common reference", default=None)
	parser.add_argument("--output", type=str, help="The name of the volumes stack that HAVE BEEN ALIGNED to a common reference", default=None)

	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")				
	
	(options, args) = parser.parse_args()	
	
	n = EMUtil.get_image_count(options.input)
	indexes = list(xrange(n))
	
	if not options.output:
		options.output = options.input.replace('.hdf','_scrambled.hdf')
		
	for i in range(n):
		num = random.choice(indexes)		
		a = EMData(options.input,num)

		a.write_image(options.output,i)
		indexes.remove(num)
		print "I have chosen taken particle %d from the original stack" %num
		print "And have put it into index %d in the randomized stack" %i
	print "DONE!"	
	return()
	
if __name__ == "__main__":
    main()

