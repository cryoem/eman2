#!/usr/bin/env python

#
# Author: Steve Ludtke, 8/24/2010
# Copyright (c) 2010- Baylor College of Medicine
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
import math
from copy import deepcopy
import os
import sys
from random import choice

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog <input> [options]

	Aligns a stack of particle images to a reference. If output is not specified, input file
	is overwritten with the aligned particles.
	"""
		
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--output", type=str, help="The name of the output class stack", default=None)
	parser.add_argument("--ref", type=str, help="Reference imageimage to use for alignment. Required", default=None)
	parser.add_argument("--refn", type=int, help="Number of the reference image in 'ref'. Default=0", default=0)
	parser.add_argument("--align",type=str,help="This is the aligner used to align particles to the previous class average. Default is None.", default="rotate_translate_flip")
	parser.add_argument("--aligncmp",type=str,help="The comparitor used for the --align aligner. Default is dot.",default="ccc")
	parser.add_argument("--ralign",type=str,help="This is the second stage aligner used to refine the first alignment. This is usually the \'refine\' aligner.", default=None)
	parser.add_argument("--raligncmp",type=str,help="The comparitor used by the second stage aligner.",default="ccc")
	parser.add_argument("--cmp",type=str,help="The comparitor used to generate quality scores for the purpose of particle exclusion in classes, strongly linked to the keep argument.", default=None)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	
	(options, args) = parser.parse_args()

	if options.output==None : options.output=args[0]
	if options.ref==None : 
		print "Must specify reference image"
		sys.exit(1)
	options.align=parsemodopt(options.align)
	options.aligncmp=parsemodopt(options.aligncmp)
	options.ralign=parsemodopt(options.ralign)
	options.raligncmp=parsemodopt(options.raligncmp)
	options.cmp=parsemodopt(options.cmp)

	logger=E2init(sys.argv, options.ppid)
	
	reference=EMData(options.ref,options.refn)
	alignstack(args[0],options.output,reference,options.align,options.aligncmp,options.ralign,options.raligncmp,options.cmp,options.verbose)
	
	E2end(logger)

def alignstack(inp,out,ref,align,aligncmp,ralign=None,raligncmp=None,fincmp=None,verbose=0):
	"""Aligns a stack of particles (from a file on disk) to a single reference image (already in memory). Writes
	results to disk. Ok for inp and out to be the same. align,aligncmp,ralign, raligncmp and fincmp are usually
	produced by parsemodopt(), producing a 2 tuple with the name of the modular class and a dictionary of options."""
	n=EMUtil.get_image_count(inp)
	
	for i in range(n):
		im=EMData(inp,i)
		aim=im.align(align[0],ref,align[1],aligncmp[0],aligncmp[1])
		if verbose : print "%d. "%i,aim["xform.align2d"],
		if ralign!=None and ralign[0]!=None:
			rparms=ralign[1]
			rparms["xform.align2d"]=aim["xform.align2d"]
			aim=im.align(ralign[0],ref,rparms,raligncmp[0],raligncmp[1])
			if verbose : print aim["xform.align2d"],
		elif verbose: print ""
		
		if fincmp!=None and fincmp[0]!=None :
			aim["match_qual"]=aim.cmp(fincmp[0],ref,fincmp[1])
			
		aim.write_image(out,i)

if __name__ == "__main__":
    main()
