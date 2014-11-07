#!/usr/bin/env python

#
# Author: Steven Ludtke, 1/30/2013 (sludtke@bcm.edu)  (rewrote older broken program)
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

# e2ptcltrace.py  1/30/2013 (rewrite)  Steven Ludtke
# This program will follow particles through refinement and assess how self consistent particle orientation assignments are


from EMAN2 import *
from math import *
import os
import sys
import copy

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <classify_xx> <classify_yy> ...

	WARNING: experimental program

	This program traces the orientation of particles through multiple iterations. Specify a list of classify_xx files for the comparison. 
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--sym",type=str,help="Symmetry to be used in searching adjacent unit cells", default="c1")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()
	
	# read all classification matrix data into a list of lists
	cls=[EMData.read_images(i) for i in args]			
	nptcl=cls[0][0]["ny"]								# particles are along the y axis, so the y height of the first image is the number of particles
	
	for i in cls[1:]: 
		if i[0]["ny"]!=nptcl:
			print "ERROR: classify files must have exactly the same number of particles"
			sys.exit(1)

	# wait until after error checking
	E2n=E2init(sys.argv,options.ppid)

	# This block creates a list of lists of Transforms representing the orientations of the reference projections for each classify file
	clsort=[]
	for i in args:
		# Now we try to get projection orientation information for each class for each classify file
		projfile=i.replace("classify","projections")
		ncls=EMUtil.get_image_count(projfile)			# may be different for different classify files
		
		orts=[]
		for i in xrange(ncls):
			orts.append(EMData(projfile,i,True)["xform.projection"])
		
		clsort.append(orts)
	
	syms=parsesym(options.sym).get_syms()			# this gets a list of Transform objects to move to each other asymmetric unit in the symmetry group
	
	for p in xrange(nptcl):
		print "%d. "%p,
		for i in xrange(1,len(args)):
			ort1=clsort[i-1][int(cls[i-1][0][0,p])]	# orientation of particle in first classify
			ort2=clsort[i][int(cls[i][0][0,p])]		# orientation of particle in second classify
			diffs=[]
			for t in syms: 
				ort2p=ort2*t
				diffs.append((ort1*ort2p.inverse()).get_rotation("spin")["omega"])		# make a list of the rotation angle to each other symmetry-related point

			diff=min(diffs)			# The angular error for the best-agreeing orientation
				
			
			if options.verbose>0 : print "%1.1f,%1.1f (%d) -> %1.1f,%1.1f (%d)\t%1.2f"%(ort1.get_rotation("eman")["alt"],ort1.get_rotation("eman")["az"],cls[i-1][0][0,p],ort2.get_rotation("eman")["alt"],ort2.get_rotation("eman")["az"],cls[i][0][0,p],diff),
		print ""
		
	
	E2end(E2n)



if __name__ == "__main__":
    main()
