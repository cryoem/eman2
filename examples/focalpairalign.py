#!/usr/bin/env python

#
# Author: Steven Ludtke, 01/18/2011 (sludtke@bcm.edu), 
# Copyright (c) 2000-2011 Baylor College of Medicine
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

from EMAN2 import *
from optparse import OptionParser
import os
import sys
import math
import random
import traceback

def get_usage():
	progname = os.path.basename(sys.argv[0])
	usage = progname + """ <close> <far> <close> <far> ... [options]
	This will take a set of focal pair images and align the second (far from focus) image to the 
	first (close to focus) image.
	"""
	return usage

def print_usage():
	
	usage = get_usage()
	print "usage " + usage;
	print "Please run '" + progname + " -h' for detailed options"

def main():
	parser=OptionParser(usage=get_usage())
	parser.add_option("--align",type="string",help="This is the aligner used to align particles to the previous class average. Default is None.", default="translational:nozero=1")
	parser.add_option("--verbose", "-v", dest="verbose", action="store", metavar="n", type="int", default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()
	align=parsemodopt(options.align)

	if len(args)<2 or len(args)%2!=0 :
		print "Arguments must be pairs of near/far focus images"
		sys.exit(1)
	

	for i in range(0,len(args),2):
		if options.verbose : print "%s -> %s"%(args[i+1],args[i])
		im1=EMData(args[i],0)
		im2=EMData(args[i+1],0)
		
		if im1["nx"]!=im2["nx"] or im1["ny"]!=im2["ny"] :
			print "Second image not same size as first, clipping"
			im2=im2.get_clip(Region(0,0,im1["nx"],im1["ny"]))
		
		im1.process_inplace("normalize.edgemean")
		im1.process_inplace("filter.highpass.gauss",{"cutoff_abs":.01})
		im2.process_inplace("normalize.edgemean")
		im2.process_inplace("filter.highpass.gauss",{"cutoff_abs":.01})
		
		ali=im2.align(align[0],im1,align[1],"ccc",{})

		if options.verbose>1 : print "Refineing"
		ali=im2.align("refine",im1,{"xform.align2d":ali.get_attr("xform.align2d")},"ccc",{})
		xform=ali["xform.align2d"]
		ali=EMData(args[i+1],0)
		ali.process_inplace("normalize.edgemean")
		ali.process_inplace("xform",{"transform":xform})

		a=args[i+1].rsplit(".",1)
		a[0]=a[0]+"_ali"
		ali.write_image(".".join(a),0)
		

if __name__ == "__main__":
    main()
