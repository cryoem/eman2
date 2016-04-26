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
	usage = """prog [options] <classmx_xx> <classmx_yy> ...

	WARNING: experimental program

	This program traces the orientation of particles through multiple iterations. Specify a list of classify_xx files for the comparison. 
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	
	parser.add_argument("--trace",type=str,help="Name of output file.", default="ptcltrace.txt")
	parser.add_argument("--sym",type=str,help="Symmetry to be used in searching adjacent unit cells", default="c1")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")

	(options, args) = parser.parse_args()
	
	cmx = []
	proj = []
	for c in args:
		if "classmx" in c:
			if os.path.isfile(c): cmx.append(c)
			p = c.replace("classmx","projections")
			if os.path.isfile(p): proj.append(p)
		else:
			print("{} is not a classmx file. Will not process.".format(f))
	
	if len(cmx) < 2: 
		print("ERROR: You must specify at least two classmx files.")
		sys.exit(1)
	
	if len(cmx) != len(proj):
		print("ERROR: Could not find matching projection files for your input classmx files.")
		sys.exit(1)
		
	cls = [EMData.read_images(c) for c in cmx] # read all classification matrix data into a list of lists
	nptcl = cls[0][0]['ny'] # particles are along the y axis
	
	for i in cls[1:]:
		if i[0]["ny"]!=nptcl:
			print "ERROR: classmx files must have exactly the same number of particles"
			sys.exit(1)

	# wait until after error checking
	E2n=E2init(sys.argv,options.ppid)

	# Create a list of lists of Transforms representing the orientations of the reference projections 
	# for each classmx file and try to get projection orientation information for each class
	
	clsort=[]
	for c,p in zip(cmx,proj): 
		ncls=EMUtil.get_image_count(p)
		orts = []
		for i in xrange(ncls):
			orts.append( EMData(p,i,True)["xform.projection"] )
		clsort.append(orts)
	
	# Get a list of Transform objects to move to each other asymmetric unit in the symmetry group
	syms=parsesym( str(options.sym) ).get_syms()
	
	with open(options.trace,"w") as outf:
		set = "placeholder.lst" # for use with 2D plot ptcl viewing
		fmt = "{:.3f}\t{:.3f}\t{:.0f}\t{:.3f}\t{:.3f}\t{:.0f}\t{:.3f} # {};{}\n"
		for p in xrange(nptcl):
			if options.verbose: 
				sys.stdout.write('\rparticle: {0:.0f} / {1:.0f}'.format(p+1,nptcl))
			outf.write("{}".format(p))
			for i in xrange(1,len(cmx)):
				ort1=clsort[i-1][int(cls[i-1][0][0,p])]	# orientation of particle in first classmx
				ort2=clsort[i][int(cls[i][0][0,p])]		# orientation of particle in second classmx
				diffs=[] # make a list of the rotation angle to each other symmetry-related point
				for t in syms:
					ort2p=ort2*t
					diffs.append((ort1*ort2p.inverse()).get_rotation("spin")["omega"])
				diff=min(diffs)	# The angular error for the best-agreeing orientation
				e1 = ort1.get_rotation("eman")
				e2 = ort2.get_rotation("eman")
				outf.write(fmt.format(e1["alt"],e1["az"],cls[i-1][0][0,p],e2["alt"],e2["az"],cls[i][0][0,p],diff,p,set))
	
	print("\nSUMMARY:\n")
	print("Argument\tMean\tConfidence\tShiftMag\tDispersion\t")
	for i,(c,p) in enumerate(zip(cmx,proj)):
		print("{}: {},{}".format(i,c,p))
	
	print("\nComparisons and assigned euler angles stored in {}.".format(options.trace))
	
	E2end(E2n)


if __name__ == "__main__":
    main()
