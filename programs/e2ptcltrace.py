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
	classes = []
	for c in args:
		if "classmx" in c:
			if os.path.isfile(c): cmx.append(c)
			p = c.replace("classmx","projections")
			cs = c.replace("classmx","classes")
			if os.path.isfile(p): proj.append(p)
			if os.path.isfile(cs): classes.append(cs)
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

	if options.verbose: print("Parsing assigned projection orientations")
	clsort=[]
	for x,p,c in zip(cmx,proj,classes):
		ncls=EMUtil.get_image_count(p)
		orts = []
		for i in xrange(ncls):
			if options.verbose:
				sys.stdout.write('\r{}, {}\t{}/{}\t'.format(x,p,i+1,ncls))
			orts.append( EMData(p,i,True)["xform.projection"] )
		clsort.append(orts)
		if options.verbose: print("")

	# Get a list of Transform objects to move to each other asymmetric unit in the symmetry group
	syms=parsesym( str(options.sym) ).get_syms()

	if options.verbose: print("Tracing particles through input classmx files")
	with open(options.trace,"w") as outf:
		for p in xrange(nptcl):
			if options.verbose:
				sys.stdout.write('\r{0:.0f} / {1:.0f}\t'.format(p+1,nptcl))
			for i in xrange(1,len(cmx)):
				ort1=clsort[i-1][int(cls[i-1][0][0,p])]	# orientation of particle in first classmx
				ort2=clsort[i][int(cls[i][0][0,p])]		# orientation of particle in second classmx

				diffs=[] # make a list of the rotation angle to each other symmetry-related point
				for t in syms:
					ort2p=ort2*t
					diffs.append((ort1*ort2p.inverse()).get_rotation("spin")["omega"])
				diff=min(diffs)	# The angular error for the best-agreeing orientation

				e1 = ort1.get_rotation("eman")
				alt1 = e1["alt"]
				az1 = e1["az"]
				cls1 = int(cls[i-1][0][0,p])

				e2 =  ort2.get_rotation("eman")
				alt2 = e2["alt"]
				az2 = e2["az"]
				cls2 = int(cls[i][0][0,p])

				clsdiff = abs(cls2-cls1)

				try:
					classes2 = cmx[i].replace("classmx","classes")
					hdr2 = EMData(classes2,cls2,True)
					idx2 = hdr2["projection_image_idx"]
					proj2 = hdr2["projection_image"]

					classes1 = cmx[i-1].replace("classmx","classes")
					hdr1 = EMData(classes1,cls1,True)
					idx1 = hdr1["projection_image_idx"]
					proj1 = hdr1["projection_image"]

					ptcl_src = hdr2["class_ptcl_src"]
					ptcl_idx = [i for i in hdr2["class_ptcl_idxs"] if i == p][0]

					hdr1 = None
					hdr2 = None

					dat = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}"
					outf.write(dat.format(alt1,az1,cls1,alt2,az2,cls2,diff,clsdiff))

					cmt = " # {};{};{};{};{};{};{};{};{};{}\n"
					outf.write(cmt.format(cls2,classes2,cls1,classes1,idx2,proj2,idx1,proj1,ptcl_idx,ptcl_src))

				except:
					pass # no data in class corresponding to projection

	if ".txt" in options.trace: kf = options.trace.replace(".txt",".key")
	else: kf = options.trace + ".key"

	with open(kf,"w") as keyfile:
		ctr = 0
		for i,c in enumerate(cmx):
			if i > 0:
				k = []
				k.append("{}:\talt from {} (input {})".format(ctr,c,i-1))
				k.append("{}\taz from {} (input {})".format(ctr+1,c,i-1))
				k.append("{}:\trotation of class {} to axis of symmetry (input {})".format(ctr+2,c,i-1))
				k.append("{}:\talt from {} (input {})".format(ctr+3,c,i))
				k.append("{}\taz from {} (input {})".format(ctr+4,c,i))
				k.append("{}:\trotation of class {} to axis of symmetry (input {})".format(ctr+5,c,i))
				k.append("{}:\tangular error for best agreeing orientation (difference between input {} and {})".format(ctr+6,i,i-1))
				k.append("{}:\tabsolute difference between class assignment {} and {}".format(ctr+7,i,i-1))
				keyfile.write("\n".join([x for x in k])+"\n")
				ctr+=len(k)

	#print("\n\nSUMMARY:")
	#print("UNDER CONSTRUCTION!")
	#print("DISPERSION")
	#print("COVERAGE")
	#print("mean/median diff")
	#print("
	#for i,(c,p) in enumerate(zip(cmx,proj)):
	#	print("{}: {},{}".format(i,c,p))

	print("Particle trace results stored in {}.\nThe file {} describes the contents of each column.".format(options.trace,kf))

	E2end(E2n)


if __name__ == "__main__":
    main()
