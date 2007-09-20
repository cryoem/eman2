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
	Computes a similarity matrix between c-input (col - projections) and r-input (row - particles) stacks of 2-D images. Images may
	optionally be aligned before comparison. Output is a matrix stored as an image with similarity value
	pairs. When used for classifiaction, c input is the references and r input are the particles."""
	parser = OptionParser(usage=usage,version=EMANVERSION)

	#parser.add_option("--apix", "-A", type="float", help="A/voxel", default=1.0)
	#parser.add_option("--box", "-B", type="string", help="Box size in pixels, <xyz> or <x>,<y>,<z>")
	#parser.add_option("--het", action="store_true", help="Include HET atoms in the map", default=False)
	parser.add_option("--align",type="string",help="The name of an 'aligner' to use prior to comparing the images", default=None)
	parser.add_option("--aligncmp",type="string",help="Name of the aligner along with its construction arguments",default="dot")
	parser.add_option("--alignr",type="string",help="The name and parameters of the second stage aligner which refines the results of the first alignment", default=None)
	parser.add_option("--alignrcmp",type="string",help="The name and parameters of the comparitor used by the second stage aligner. Default is dot.",default="dot")
	parser.add_option("--cmp",type="string",help="The name of a 'cmp' to be used in comparing the aligned images", default="dot:normalize=1")
	parser.add_option("--range",type="string",help="Range of images to process (c0,r0,c1,r1) c0,r0 inclusive c1,r1 exclusive", default=None)
	parser.add_option("--saveali",action="store_true",help="Save alignment values, output is c x r x 4 instead of c x r x 1",default=False)
	parser.add_option("--verbose","-v",action="store_true",help="Verbose display during run",default=False)
	parser.add_option("--init",action="store_true",help="Initialize the output matrix file before performing 'range' calculations",default=False)
	(options, args) = parser.parse_args()
	
	if len(args)<3 : parser.error("Input and output files required")
	
	if os.path.exists(args[2]):
		parser.error("File %s exists, will not write over, exiting" %args[2])
	
	E2n=E2init(sys.argv)
	
	options.align=parsemodopt(options.align)
	options.aligncmp=parsemodopt(options.aligncmp)
	options.alignr=parsemodopt(options.alignr)
	options.alignrcmp=parsemodopt(options.alignrcmp)
	options.cmp=parsemodopt(options.cmp)

	clen=EMUtil.get_image_count(args[0])
	rlen=EMUtil.get_image_count(args[1])
	
	if options.init:
		a=EMData()
		if options.saveali : a.set_size(clen,rlen,4)
		else : a.set_size(clen,rlen,1)
		a.to_zero()
		a.write_image(args[2])
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
	if options.verbose: print "Computing Similarities"

	# Read all c images, then read and compare one r image at a time
	cimgs=EMData.read_images(args[0],range(*crange))
	rimg=EMData()
	for r in range(*rrange):
		if options.verbose: 
			print "%d/%d\r"%(r,rrange[1]),
			sys.stdout.flush()
		rimg.read_image(args[1],r)
		row=cmponetomany(cimgs,rimg,options.align,options.aligncmp,options.cmp, options.alignr, options.alignrcmp)
		for c,v in enumerate(row):
			mxout.set_value_at(c,r,0,v[0])
		
		if options.saveali :
			for c,v in enumerate(row):
				#print "%f %f %f " %(v[1],v[2],v[3])
				mxout.set_value_at(c,r,1,v[1])
				mxout.set_value_at(c,r,2,v[2])
				mxout.set_value_at(c,r,3,v[3])
	
	if options.verbose : print"\nSimilarity computation complete"
	
	# write the results into the full-sized matrix
	if crange==[0,clen] and rrange==[0,rlen] :
		mxout.write_image(args[2],0)
	else :
		if options.saveali : mxout.write_image(args[2],0,IMAGE_UNKNOWN,0,Region(crange[0],rrange[0],0,crange[1]-crange[0],rrange[1]-rrange[0],4))
		else : mxout.write_image(args[2],0,IMAGE_UNKNOWN,0,Region(crange[0],rrange[0],crange[1]-crange[0],rrange[1]-rrange[0]))
	
	E2end(E2n)
	
def cmponetomany(reflist,target,align=None,alicmp=("dot",{}),cmp=("dot",{}), alignr=None, alircmp=("dot",{})):
	"""Compares one image (target) to a list of many images (reflist). Returns """
	
	ret=[None for i in reflist]
	for i,r in enumerate(reflist):
		if align[0] :
			ta=target.align(align[0],r,align[1],alicmp[0],alicmp[1])
			#ta.debug_print_params()
			
			if alignr[0]:
				alignr[1]["az"] = ta.get_attr_default("align.az",0)-1
				alignr[1]["dx"] = ta.get_attr_default("align.dx",0)-1
				alignr[1]["dy"] = ta.get_attr_default("align.dy",0)-1
				ta = target.align(alignr[0],r,alignr[1],alircmp[0],alircmp[1])
				
			ret[i]=(ta.cmp(cmp[0],r,cmp[1]),ta.get_attr_default("align.dx",0),ta.get_attr_default("align.dy",0),ta.get_attr_default("align.az",0))
		else :
			ret[i]=(target.cmp(cmp[0],r,cmp[1]),0,0,0)
		
	return ret

				
if __name__ == "__main__":
    main()
