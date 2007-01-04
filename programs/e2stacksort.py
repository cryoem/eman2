#!/usr/bin/env python

#
# Author: Steven Ludtke, 01/03/07 (sludtke@bcm.edu)
# Copyright (c) 2000-2007 Baylor College of Medicine
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

# e2stacksort.py  01/03/07  Steven Ludtke
# This program will sort a stack of images based on some similarity criterion


from EMAN2 import *
from optparse import OptionParser
from math import *
import os
import sys

def sortstack(stack,cmptype,cmpopts):
	"""Sorts a list of images based on a standard 'cmp' metric. cmptype is the name
	of a valid cmp type. cmpopts is a dictionary. Returns a new (sorted) stack.
	The original stack is destroyed."""
	
	ret=[stack[0]]
	del stack[0]
	while (len(stack)>0) :
		best=(1.0e38,-1)
		for i,im in enumerate(stack):
			c=ret[-1].cmp(cmptype,im,cmpopts)+im.cmp(cmptype,ret[-1],cmpopts)	# symmetrize results
			if c<best[0] or best[1]<0 : best=(c,i)
		ret.append(stack[best[1]])
		del stack[best[1]]
		print "%d.\t%d  (%1.4f)"%(len(ret)-1,best[1],best[0])

	return ret

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <input_stack> <output_stack>
	
This program will sort a stack of images based on some similarity criterion"""

	parser = OptionParser(usage=usage,version=EMANVERSION)

#	parser.add_option("--tilt", "-T", type="float", help="Angular spacing between tilts (fixed)",default=0.0)
#	parser.add_option("--maxshift","-M", type="int", help="Maximum translational error between images (pixels), default=64",default=64.0)
#	parser.add_option("--mode",type="string",help="centering mode 'modeshift', 'censym' or 'region,<x>,<y>,<clipsize>,<alisize>",default="censym")
#	parser.add_option("--twopass",action="store_true",default=False,help="Skip automatic tilt axis location, use fixed angle from x")
	
	(options, args) = parser.parse_args()
	if len(args)<2 : parser.error("Input and output files required")
	
	a=EMData.read_images(args[0])
	b=sortstack(a,"optvariance",{"keepzero":1,"matchfilt":1})
	for i,im in enumerate(b): im.write_image(args[1],i)

if __name__ == "__main__":
    main()
