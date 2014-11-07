#!/usr/bin/env python

#
# Author: Steven Ludtke, 10/04/2013 (sludtke@bcm.edu)
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

# e2classvsproj.py  Steven Ludtke

from EMAN2 import *
from math import *
import os
import sys
import traceback

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <classes> <projection stack or 3Dmap> <output>
	Compares each class-average (or particle) in the classes input stack to each projection in the 'projection stack'. If a
	3-D map is provided as the second input, then projections are generated using the specified angular spacing.

	e2ptclvsmap.py is similar, but also sorts the results """


	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--ang",type=float,help="Angle between projections if 3Dmap provided",default=10.0)
	parser.add_argument("--sym", dest="sym", default="c1", help="Set the symmetry; if no value is given then the model is assumed to have no symmetry.\nChoices are: i, c, d, tet, icos, or oct.")
	parser.add_argument("--align",type=str,help="The name of an 'aligner' to use prior to comparing the images", default="rotate_translate_flip")
	parser.add_argument("--aligncmp",type=str,help="Name of the aligner along with its construction arguments",default="ccc")
	parser.add_argument("--ralign",type=str,help="The name and parameters of the second stage aligner which refines the results of the first alignment", default="refine")
	parser.add_argument("--raligncmp",type=str,help="The name and parameters of the comparitor used by the second stage aligner. Default is dot.",default="ccc")
	parser.add_argument("--cmp",type=str,help="The name of a 'cmp' to be used in comparing the aligned images", default="ccc")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()

	if len(args)<3 : parser.error("Input and output files required")

	options.align=parsemodopt(options.align)
	options.aligncmp=parsemodopt(options.aligncmp)
	options.ralign=parsemodopt(options.ralign)
	options.raligncmp=parsemodopt(options.raligncmp)
	options.cmp=parsemodopt(options.cmp)

	# initialize projections either a stack from a file, or by making projections of a 3-D map
	projs=EMData.read_images(args[1])
	if len(projs)==1 :
		mdl=projs[0]
		projs=[]
		if mdl["nz"]==1 : 
			print "Error, you must provide either a stack of 2-D images or a 3-D map as the second argument"
			sys.exit(1)

		print "Generating projections with an angular step of {} and symmetry {}".format(options.ang,options.sym)
		
		E2n=E2init(sys.argv, options.ppid)
		mdl.process_inplace("normalize.edgemean")
		sym_object = parsesym(options.sym)
		eulers = sym_object.gen_orientations("eman", {"delta":options.ang,"inc_mirror":0,"perturb":0})
		for i,euler in enumerate(eulers):
			p=mdl.project("standard",euler)
			p.set_attr("xform.projection",euler)
			p.set_attr("ptcl_repr",0)
			p.process_inplace("normalize.edgemean")
			projs.append(p)
			if options.verbose : print i,euler
	else:
		E2n=E2init(sys.argv, options.ppid)

	# Now find the best match for each particle. We could use e2simmx, but more efficient to just do it in place (though not parallel this way)
	nptcl=EMUtil.get_image_count(args[0])
	for i in xrange(nptcl):
		best=None
		ptcl=EMData(args[0],i)
		ptcl.process_inplace("normalize.edgemean")
		for proj in projs:
			aligned=proj.align(options.align[0],ptcl,options.align[1],options.aligncmp[0],options.aligncmp[1])
			if options.ralign != None: # potentially employ refine alignment
				refine_parms=options.ralign[1]
				refine_parms["xform.align2d"] = aligned.get_attr("xform.align2d")
				proj.del_attr("xform.align2d")
				aligned = proj.align(options.ralign[0],ptcl,refine_parms,options.raligncmp[0],options.raligncmp[1])
			
			c=ptcl.cmp(options.cmp[0],aligned,options.cmp[1])
			if best==None or c<best[0] : best=(c,aligned)
			if options.verbose>1 : print i,projs.index(proj),c
		
		best[1].write_image(args[2],i*2)
		ptcl.write_image(args[2],i*2+1)
		if options.verbose: print "Class-average {} with projection {}".format(i,str(best[1]["xform.projection"]))
		

	E2end(E2n)


if __name__ == "__main__":
    main()
