#!/usr/bin/env python

#
# Author: Steven Ludtke, 02/03/2007 (sludtke@bcm.edu)
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

# e2simmx.py  02/03/2007	Steven Ludtke
# This program computes a similarity matrix between two sets of images

from EMAN2 import *
from optparse import OptionParser
from math import *
import os
import sys

def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <c input> <r input> <output>
	Computes a similarity matrix between c-input (col) and r-input (row) stacks of 2-D images. Images may
	optionally be aligned before comparison. Output is a matrix stored as an image with similarity value
	pairs"""
	parser = OptionParser(usage=usage,version=EMANVERSION)

	#parser.add_option("--apix", "-A", type="float", help="A/voxel", default=1.0)
	#parser.add_option("--box", "-B", type="string", help="Box size in pixels, <xyz> or <x>,<y>,<z>")
	#parser.add_option("--het", action="store_true", help="Include HET atoms in the map", default=False)
	parser.add_option("--align",type="string",help="The name of an 'aligner' to use prior to comparing the images", default=None)
	parser.add_option("--aligncmp",type="string",help="Name of a 'cmp' to be used in the aligner",default="dot")
	parser.add_option("--cmp",type="string",help="The name of a 'cmp' to in comparing the aligned images", default="dot(normalize=1)")
	parser.add_option("--range",type="string",help="Range of images to process (c0,r0,c1,r1) c0,r0 inclusive c1,r1 exclusive", default=None)
	parser.add_option("--saveali",action="store_true",help="Save alignment values, output is c x r x 4 instead of c x r x 1",default=False)
	parser.add_option("--init",action="store_true",help="Initialize the output matrix file before performing 'range' calculations",default=False)
	(options, args) = parser.parse_args()
	
	if len(args)<3 : parser.error("Input and output files required")
	
	E2n=E2init(sys.argv)
	
	options.align=parsemodopt(options.align)
	options.aligncmp=parsemodopt(options.aligncmp)
	options.cmp=parsemodopt(options.cmp)
	
	clen=EMUtil.get_image_count(args[0])
	rlen=EMUtil.get_image_count(args[1])
	
	if options.init:
		a=EMData()
		if options.saveali : a.set_size(clen,rlen,4)
		else : a.set_size(clen,rlen,1)
		a.to_zero()
		a.write_image(arg[2])
		E2end(E2n)
		sys.exit(0)
	
	# Compute range in c and r
	if options.range :
		crange=options.range.split(",")[0:4:2]
		rrange=options.range.split(",")[1:4:2]
		crange[0]=int(crange[0])
		crange[1]=int(crange[1])
		rrange[0]=int(rrange[0])
		rrange[1]=int(rrange[1])
	else:
		crange=[0,clen]
		rrange=[0,rlen]

	# initialize output array
	mxout=EMData()
	if options.saveali : mxout.set_size(crange[1]-crange[0],rrange[1]-rrange[0],4)
	else : mxout.set_size(crange[1]-crange[0],rrange[1]-rrange[0],1)
	mxout.to_zero()

	# Read all c images, then read and compare one r image at a time
	cimgs=EMData.read_images(args[0],range(*crange))
	rimg=EMData()
	for r in range(*rrange):
		rimg.read_image(args[1],r)
		row=cmponetomany(cimgs,rimg,options.align,options.aligncmp,options.cmp)
		for c,v in enumerate(row):
			mxout.set_value_at(c,r,0,v[0])
		
		if options.saveali :
			for c,v in enumerate(row):
				mxout.set_value_at(c,r,1,v[1])
				mxout.set_value_at(c,r,2,v[2])
				mxout.set_value_at(c,r,3,v[3])
	
	# write the results into the full-sized matrix
	if options.saveali : mxout.write_image(args[2],0,IMAGE_UNKNOWN,Region(crange[0],rrange[0],0,crange[1]-crange[0],rrange[1]-rrange[0],4))
	else : mxout.write_image(args[2],0,IMAGE_UNKNOWN,Region(crange[0],rrange[0],crange[1]-crange[0],rrange[1]-rrange[0]))
	
	E2end(E2n)
	
def cmponetomany(reflist,target,align=None,alicmp=("dot",{}),cmp=("dot",{})):
	"""Compares one image (target) to a list of many images (reflist). Returns """
	
	ret=[None for i in reflist]
	for i,r in enumerate(reflist):
		if align[0] : 
			ta=target.align(align[0],r,align[1],alicmp[0],alicmp[1])
			ret[i]=(ta.cmp(cmp[0],r,cmp[1]),ta.get_attr_default("translational.dx",0),ta.get_attr_default("translational.dy",0),ta.get_attr_default("rotational",0))
		else : 
			ret[i]=(target.cmp(cmp[0],r,cmp[1]),0,0,0)
		
	return ret

				
if __name__ == "__main__":
    main()
