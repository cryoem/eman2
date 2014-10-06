#!/usr/bin/env python
# This program performs simple processing of .LST files

# Author: Steven Ludtke, 10/06/14 (sludtke@bcm.edu)
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
	usage = """Usage:\nproclst.py [options] <lst 1> <lst 2> ... \nSimple manipulations of LST files. If your goal is to produce an actual image file rather than the
sort of virtual stack represented by .lst files, use e2proc2d.py or e2proc3d.py instead. Those programs will treat LST files as normal image files for input.\n."""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	####################
#	parser.add_argument("--average", action="store_true", help="Averages all input images (without alignment) and writes a single output image")
	parser.add_argument("--merge",type=str,help="Specify the output name here. This will concatenate all of the input .lst files into a single output",default=None)
	parser.add_argument("--mergesort",type=str,help="Specify the output name here. This will merge all of the input .lst files into a single (resorted) output",default=None)
	parser.add_argument("--retype",type=str,help="If a lst file is referencing a set of particles from particles/imgname__oldtype.hdf, this will change oldtype to the specified string in-place (modifies input files)",default=None)
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n", type=int, help="verbose level [0-9], higner number means higher level of verboseness",default=1)

	(options, args) = parser.parse_args()
	
	if len(args)<1 : 
		parser.error("At least one lst file required")
		sys.eixt(1)

	if options.retype != None:
		
		# if the user provided the leading __ for us, we strip it off and add it back later
		if options.retype[:2]=="__" : 
			options.retype=options.retype[2:]
		
		for f in args:
			if options.verbose : print "Processing ",f
			lst=LSXFile(f,True)
			
			a=lst.read(0)
			if a[1][:10]!="particles/" :
				print "To use the --retype option, the .lst file must reference image files in particles/*"
				
			if options.verbose>1 : 
				b=base_name(a[1])
				print "{} -> {}".format(a[1],b+"__"+options.retype+".hdf")
			
			# loop over the images in the lst file
			for i in xrange(len(lst)):
				im=lst.read(i)
				outname="particles/{}__{}.hdf".format(base_name(im[1]),options.retype)
				lst.write(i,im[0],outname,im[2])
		
			lst.normalize()			# clean up at the end

			if options.verbose>1 : print len(lst)," particles adjusted"

		if options.verbose : print "Done processing {} files".format(len(args))

	if options.merge!=None:
		
		# create/update output lst
		lsto=LSXFile(options.merge)
		ntot=0
		
		# loop over input files
		for f in args:
			lst=LSXFile(f,True)
			ntot+=len(lst)
			
			for i in xrange(len(lst)):
				im=lst.read(i)
				lsto.write(-1,im[0],im[1],im[2])

		if options.verbose : print "{} particles added to {}".format(ntot,options.merge)

	if options.mergesort!=None:
		# create/update output lst
		lsto=LSXFile(options.mergesort)
		ntot=0
		
		# loop over input files
		ptcls=[]
		for f in args:
			lst=LSXFile(f,True)
			ntot+=len(lst)
			
			for i in xrange(len(lst)):
				im=lst.read(i)
				ptcls.append((im[1],im[0],im[2]))
			ptcls.sort()
			
		for i in ptcls:
			lsto.write(-1,i[1],i[0],i[2])

		if options.verbose : print "{} particles in {}".format(ntot,options.mergesort)


if __name__ == "__main__":
	main()
