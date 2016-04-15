#!/usr/bin/env python

#
# Author: Steve Ludtke, 7/5/14 
# Copyright (c) 2000- Baylor College of Medicine
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

from math import *
import os
import sys
from EMAN2db import db_check_dict
from EMAN2 import *

def main():
	
	progname = os.path.basename(sys.argv[0])
	usage = """prog <references> <particles> <classmx> [options]
	
	*** THIS PROGRAM IS NOT YET FUNCTIONAL ***

	This program classifies a set of particles based on a set of references (usually projections). Historically this
	was done in two steps, first by e2simmx.py or e2simmx2stage.py then by e2classify. However, this approach produced
	very large intermediate files (simmx file), and was inefficient in a number of other ways. While the similartiy matrix
	can sometimes be useful for other purposes, this is not true in the vast majority of situtations.
	
	This program performs multi-stage classification, all within this single program. It first classifies the references
	themselves against one another. This produces first-stage references, each of which is associated with several
	(similar) individual references. Each particle is first classified against the first stage references, then a different
	algorithm is used to discriminate among the more similar subset of references.

	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--sep", type=int, help="The number of classes a particle can contribute towards (default is 1)", default=1)
	parser.add_argument("--force", "-f",dest="force",default=False, action="store_true",help="Force overwrite the output file if it exists")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	
	if options.nofilecheck: options.check = True

	if (len(args)<3 ): parser.error("Please specify <references> <particles> <classmx file>")
	
	E2n=E2init(sys.argv, options.ppid)

	if os.path.exists(args[2]):
		remove_file(args[2])

	nref=EMUtil.get_image_count(args[0])
	nptcl=EMUtil.get_image_count(args[1])

	tmp=EMData(args[0],0,True)
	nx,ny=tmp["nx"],tmp["ny"]


	E2end(E2n)

	print "Classification complete, writing classmx"
	clsmx[0].write_image(args[1],0)
	clsmx[1].write_image(args[1],1)
	if num_sim>=5 :
		clsmx[2].write_image(args[1],2)
		clsmx[3].write_image(args[1],3)
		clsmx[4].write_image(args[1],4)
		clsmx[5].write_image(args[1],5)
		if num_sim>5 : clsmx[6].write_image(args[1],6)
		
	

def footprint(image):
	"""Computes a "footprint" for an image. Note that it normalizes the input as a side-effect !"""
	image.process_inplace("normalize.edgemean")
	fp=image.window_center(image["nx"]*2)
	fp=fp.calc_ccf(fp)
	fp.process_inplace("xform.phaseorigin.tocenter")
	fp.process_inplace("math.rotationalsubtract")
	fp=fp.unwrap(4,image["nx"]/2)
	fp=fp.do_fft()
	fp.ri2inten()
#	return fp.calc_ccfx(fp,0,fp["ny"],1)
	

if __name__ == "__main__":
    main()
