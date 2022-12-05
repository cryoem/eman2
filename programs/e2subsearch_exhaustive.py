#!/usr/bin/env python
#
# Author: Steven Ludtke, 12/3/2022
# Copyright (c) 2022 Baylor College of Medicine
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


from EMAN2 import *
from math import *
import numpy as np
import queue
from time import time

def makeproj(jsd,mdl,eul,newsize):
	x0=(mdl["nx"]-newsize)//2
	for e in eul:
		prj=mdl.project("standard",e)
		prj=prj.get_clip(Region(x0,x0,newsize,newsize))
		prj.process_inplace("xform.phaseorigin.tocorner")
		jsd.put(prj.do_fft())

def compute_ccfs(jsd,

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """e2subsearch_exhaustive <particles> <template3d> [options]
	This program will do an exhaustive search for 3-D subunit projections across individual 2-D images. The template 3-D map will be exhaustively rotated in 3-D and
	the corresponding projection cross-correlated with each 2-D image to find the best match(es), typically for extraction of "subparticles".

	Please insure that A/pix is the same for images and the template, and that the particles referenced in the .lst input have been phase-flipped.
	"""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_pos_argument(name="particles",help=".lst file containing particle stack to be searched", default="", guitype='filebox', browser="EMRefine2dTable(withmodal=True,multiselect=False)",  filecheck=False, row=0, col=0,rowspan=1, colspan=2)
	parser.add_pos_argument(name="template",help="3D volume template to search for", default="", guitype='filebox', browser="EMRefine2dTable(withmodal=True,multiselect=False)",  filecheck=False, row=1, col=0,rowspan=1, colspan=2)
	parser.add_argument("--angstep", type=float, default=-1.0, help="Angular step for search in degrees, default=auto (-1)")
	parser.add_argument("--extractbox", type=int, default=-1, help="Box size to use for extracted matches, default=auto (-1)")
	parser.add_argument("--maxmatches", type=int, default=-1, help="Maximum number of matches to find in each particle image")
	parser.add_argument("--sym", dest = "sym", help = "Symmetry of the template, if any. NOT the symmetry of the target particles.",default="c1")
	parser.add_argument("--threads", default=4,type=int,help="Number of threads to run in parallel on a single computer",guitype='intbox', row=2, col=0, rowspan=1, colspan=1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbosity [0-9], higher number means more detailed output")


	# parser.add_pos_argument(name="refine_xx",help="Name of a completed refine_xx folder.", default="", guitype='filebox', browser="EMRefine2dTable(withmodal=True,multiselect=False)",  filecheck=False, row=0, col=0,rowspan=1, colspan=2, mode='evalptcl')
	# parser.add_argument("--timing", default=False, action="store_true", help="Report on the time required for each step of each refinement run")
	# parser.add_argument("--timingbypath", default=False, action="store_true", help="Report on the CPU time required in each refine_xx folder")
	# parser.add_argument("--resolution", default=False, action="store_true", help="generates a resolution and convergence plot for a single refinement run.")
	# parser.add_argument("--resolution_all", default=False, action="store_true", help="generates resolution plot with the last iteration of all refine_xx directories")
	# parser.add_argument("--resolution_vsref", type=str, default=None, help="Computes the FSC between the last iteration of each refine_xx directory and a specified reference map. Map must be aligned, but will be rescaled if necessary.")
	# parser.add_argument("--evalptclqual", default=False, action="store_true", help="Evaluates the particle-map agreement using the refine_xx folder name. This may be used to identify bad particles.",guitype='boolbox', row=8, col=1, rowspan=1, colspan=1, mode='evalptcl[True]')
	# parser.add_argument("--evalclassqual", default=False, action="store_true", help="Evaluates the class-average-projection agreement using the refine_xx folder name.",guitype='boolbox', row=8, col=2, rowspan=1, colspan=1, mode='evalptcl[False]')
	# parser.add_argument("--extractorientptcl", default=None, type=str, help="Extracts the particles from a refinement with per-particle orientation information. If HDF output, will store as xform_align3d in header.")
	# parser.add_argument("--anisotropy", type=int, default=-1, help="Specify a class-number (more particles better). Will use that class to evaluate magnification anisotropy in the data. ")
	# parser.add_argument("--evalclassdetail", default=False, action="store_true", help="If specified with evalclassqual, will generate individual FRC curves for each class average in the even subset")
	# parser.add_argument("--includeprojs", default=False, action="store_true", help="If specified with --evalptclqual, projections will be written to disk for easy comparison.",guitype='boolbox', row=8, col=0, rowspan=1, colspan=1, mode='evalptcl[True]')
	# parser.add_argument("--iter", type=int, default=None, help="If a refine_XX folder is being used, this selects a particular refinement iteration. Otherwise the last complete iteration is used.")
	# parser.add_argument("--mask",type=str,help="Mask to be used to focus --evalptclqual and other options. May be useful for separating heterogeneous data.", default=None)
	# parser.add_argument("--sym",type=str,help="Symmetry to be used in searching adjacent unit cells, default from refine_xx parms", default=None)
	# #parser.add_argument("--parmcmp",  default=False, action="store_true",help="Compare parameters used in different refinement rounds")
	# #parser.add_argument("--parmpair",default=None,type=str,help="Specify iter,iter to compare the parameters used between 2 itertions.")
	# parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	# #options associated with e2refine.py
	# #parser.add_argument("--iter", dest = "iter", type = int, default=0, help = "The total number of refinement iterations to perform")
	# #parser.add_argument("--check", "-c", dest="check", default=False, action="store_true",help="Checks the contents of the current directory to verify that e2refine.py command will work - checks for the existence of the necessary starting files and checks their dimensions. Performs no work ")
	# #parser.add_argument("--input", dest="input", default=None,type=str, help="The name of the image containing the particle data")

	(options, args) = parser.parse_args()

	template=EMData(args[1])
	if template["nz"]==1 : error_exit("template must be a 3-D volume")

	nthr=options.threads
	if nthr<4: nthr=4	# using 4 threads minimum

	nx=EMData(args[0],0,True)["nx"]

	if options.angstep<=0 :
		angstep=atan(4.0/template["nx"])*180.0/pi		# 2 pixels at edge of template image
		if options.verbose: print(f"Using angular step {angstep:1.2f}")
	else: angstep=options.angstep

	if angstep<2.0 : print(f"WARNING: {angstep:1.3g} is a very small angular step and may take prohibitively long to run")

	sym_object = parsesym(options.sym)
	eulers = sym_object.gen_orientations("eman", {"delta":angstep,"inc_mirror":1,"perturb":0})

	jsd=queue.Queue(0)
	thrds=[threading.Thread(target=makeproj,args=(jsd,template,eulers[i::nthr]),nx) for i in range(nthr)]

	t0=time()
	for t in thrds: t.start()
	for t in thrds: t.join()

	projs=[]
	while not jsd.empty(): projs.append(jsd.get())

	if options.verbose: print(f"{len(projs)} projections generated in {time()-t0}s")





if __name__ == "__main__":
	main()
