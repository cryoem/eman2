#!/usr/bin/env python
from __future__ import print_function
# This program performs simple processing of .LST files

# Author: Steven Ludtke, 10/06/14 (sludtke@bcm.edu), modified: May 15, 2017 (Jesus GalazMontoya)
# Copyright (c) 2014- Baylor College of Medicine
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

from EMAN2 import *
from math import *
import os
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """Usage:\ne2spt_buildsets.py <stack 1> <stack 2> ... [options]\nCombine 3D particles from multiple tomograms using a LST file."""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)	
	parser.add_pos_argument(name="particle_stacks",help="List the SPT particle stacks in the particles3d directory to include in the labeled set.", default="", guitype='filebox', browser="EMSPTParticleTable(withmodal=True,multiselect=True)", row=0, col=0,rowspan=1, colspan=3, mode="sets")
	parser.add_argument("--name", type=str, default=None, help="Specify the name of the set to be generated. Example: 'ribosome' will create the file 'sets/ribosome.lst'", guitype='strbox',row=2, col=0,rowspan=1, colspan=2, mode="sets")
	parser.add_argument("--inplace", action="store_true", default=False, help="If the set indicated in --name already exists, this will prevent appending to it. rather, the file will be modified in place. Use this only when an existing set is present.",guitype='boolbox',  row=2, col=2,rowspan=1, colspan=1, mode="sets[False]")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-2)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, help="verbose level [0-9], higner number means higher level of verboseness",default=1)
	(options, args) = parser.parse_args()

	if len(args)<1 :
		parser.error("At least one file required")
		sys.exit(1)

	# check box size of all inputs

	logid=E2init(sys.argv,options.ppid)

	if '.lst' not in options.name and '.lsx' not in options.name:
		#if options.verbose:
		#	print("Note: A .lst/.lsx extension was not included in --name. Adding correct extension automatically.")
		options.name = options.name.split(".")[0] + ".lst"

	if "sets/" not in options.name : options.name = "sets/{}".format(options.name)
	
	lst=LSXFile(options.name,False)

	for f in args:
		n=EMUtil.get_image_count(f)
		if f.endswith(".lst"):
			lstin=LSXFile(f,True)
			fromlst=True
		else : fromlst=False
		
		indxsinclude = xrange(n) #by default, assume all particles in input file will be part of output lsx; otherwise, modify indexes to include according to options

		if options.verbose : print("Processing {} images in {}".format(len(indxsinclude),f))
		
		kk=0
		for i in indxsinclude:
		
			if fromlst:
				ln = lstin.read(i)
				if options.inplace : lst.write(kk,ln[0],ln[1],ln[2])
				else : lst.write(-1,ln[0],ln[1],ln[2])
			else:
				if options.inplace : lst.write(kk,i,f)
				else : lst.write(-1,i,f)
			kk+=1

	E2end(logid)




if __name__ == "__main__":
	main()
