#!/usr/bin/env python

#
# Author: Steven Ludtke, 3/29/15 (sludtke@bcm.edu)
# Copyright (c) 2000-2015 Baylor College of Medicine
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
import time
import os
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <sets/file.lst> <classmx_xx.hdf> <box3d 1> <box3d 2> ...

This program is part of a pipeline to improve centering of subtomograms prior to alignment and averaging.
This is a typical workflow:

1) select 3-D particles and store locations in .box3d files with names agreeing with the name of the tomogram, eg
  rawtomograms/file0001.hdf -> rawtomograms/file0001.box3d

2) extract the 3-D particles from the tomograms, and generate Z-projections corresponding to each particle (this may be done automatically by e2spt_boxer)

3) insure that the Z-projections are normalized, filtered and have postive (white) contrast against the background, eg:
  e2proc2d.py spt_particles/file0001_prjsz.hdf spt_particles/file0001_prjsz.hdf --inplace --process filter.lowpass.gauss:cutoff_abs=0.1 --process filter.highpass.gauss:cutoff_abs=0.02 --mult -1 --process normalize.edgemean

4) build a set containing all of the z projection stacks, eg
  e2proclst.py spt_particles/file0001_prjsz.hdf spt_particles/file0002_prjsz.hdf --create sets/all_prjsz.lst

5) run e2refine2d.py

6) run this program, eg
  e2boxadjust.py sets/all_prjsz.lst r2d_01/classmx_08.hdf rawtomograms/*box3d -v 2

7) re-extract the particles using the new _cen.box3d files

"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	#parser.add_argument("--normproj",action="store_true",help="Normalize the projections resulting from 'project', such that the length of each vector is 1",default=False)
	#parser.add_argument("--normalize",type=str,help="Normalize the input images using the named processor. Specify None to disable.",default="normalize.unitlen")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	
	(options, args) = parser.parse_args()
#	if len(args)>0 : parser.error("e2basis.py takes no arguments, only options")

	logid=E2init(sys.argv,options.ppid)
	
	# The classmx file contains particle alignment information
	classmx=EMData(args[1],0)
	nptcl=classmx["ny"]
	cmxtx=EMData(args[1],2)
	cmxty=EMData(args[1],3)
	cmxalpha=EMData(args[1],4)
	cmxmirror=EMData(args[1],5)
	
	print "Classmx has info on ",nptcl," particles"
	
	# The files containing the particle locations
	boxfiles=[base_name(args[i],nodir=True) for i in xrange(2,len(args))]
	
	# The .lst file allowing us to reference original files from the information in cls files 
	lsx=LSXFile(args[0])
	
	lpfile=None
	skipfile=True
	for p in xrange(nptcl):
		# The number and file of particle N
		pn,pfile,com = lsx[p]
		
		if pfile!=lpfile:
			# write the boxes from the last file
			if not skipfile: write_boxes(curboxfile,curboxes)

			skipfile=False
			pfileb=base_name(pfile,nodir=True)
			if not pfileb in boxfiles :
				print "No box file found for: ",pfileb
				lpfile=pfile
				skipfile=True
				continue
			
			
			# This is the file containing the box locations for this range of particles
			curboxfile=args[boxfiles.index(pfileb)+2]
			p0=p
			if options.verbose: print pfileb,"->",curboxfile
			
			# These are the box locations within that file
			curboxes=[[int(j) for j in i.split()] for i in file(curboxfile,"r") if i[0]!="#"]
			lpfile=pfile
		else:
			if skipfile : continue		# we've already identified this as a file we don't have box locations for
		
		ptclxf=Transform({"type":"2d","alpha":cmxalpha[0,p],"mirror":int(cmxmirror[0,p]),"tx":cmxtx[0,p],"ty":cmxty[0,p]})
		pt2d=ptclxf.get_pre_trans_2d()
		curboxes[p-p0][0]-=pt2d[0]
		curboxes[p-p0][1]-=pt2d[1]
		
		
		if options.verbose>1: 
			try: print "{}) {}: {}\t {:d},{:d}".format(p,p-p0,pfileb,int(pt2d[0]),int(pt2d[1]))
			except: pass
	
	if not skipfile: write_boxes(curboxfile,curboxes)
	
	E2end(logid)

def write_boxes(curboxfile,curboxes):
	print "Writing updated boxes for: ",curboxfile
	out=file(curboxfile.split(".")[0]+"_cen.box3d","w")
	for b in curboxes: out.write("{:d}\t{:d}\t{:d}\n".format(int(b[0]),int(b[1]),int(b[2])))

if __name__== "__main__":
	main()
	
