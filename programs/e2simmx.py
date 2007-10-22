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
	parser.add_option("--ralign",type="string",help="The name and parameters of the second stage aligner which refines the results of the first alignment", default=None)
	parser.add_option("--raligncmp",type="string",help="The name and parameters of the comparitor used by the second stage aligner. Default is dot.",default="dot")
	parser.add_option("--cmp",type="string",help="The name of a 'cmp' to be used in comparing the aligned images", default="dot:normalize=1")
	parser.add_option("--range",type="string",help="Range of images to process (c0,r0,c1,r1) c0,r0 inclusive c1,r1 exclusive", default=None)
	parser.add_option("--saveali",action="store_true",help="Save alignment values, output is c x r x 4 instead of c x r x 1",default=False)
	parser.add_option("--verbose","-v",action="store_true",help="Verbose display during run",default=False)
	parser.add_option("--lowmem",action="store_true",help="prevent the bulk reading of the reference images - this will save memory but potentially increase CPU time",default=False)
	parser.add_option("--init",action="store_true",help="Initialize the output matrix file before performing 'range' calculations",default=False)
	parser.add_option("--force", "-f",dest="force",default=False, action="store_true",help="Force overwrite the output file if it exists")
	parser.add_option("--nofilecheck",action="store_true",help="Turns file checking off in the check functionality - used by e2refine.py.",default=False)
	parser.add_option("--check","-c",action="store_true",help="Performs a command line argument check only.",default=False)
	
	(options, args) = parser.parse_args()
	
	if len(args)<3 : parser.error("Input and output files required")
	
	if (options.check): options.verbose = True # turn verbose on if the user is only checking...
	
	if (options.verbose):
		print ""
		print "### Testing to see if I can run e2simmx.py"
		
	if (options.nofilecheck == False):
		options.reffile=args[0]
		options.datafile = args[1]
		options.outfile = args[2]
	
	error = check(options,True)
	
	if (options.verbose):
		if (error):
			print "e2simmx.py test.... FAILED"
		else:
			print "e2simmx.py test.... PASSED"
		
	if ( options.check or error ) : exit(1)
	
	E2n=E2init(sys.argv)
	
	# just remove the file - if the user didn't specify force then the error should have been found in the check function
	if os.path.exists(options.outfile):
		if (options.force):
			remove_file(options.outfile)
	
	options.align=parsemodopt(options.align)
	options.aligncmp=parsemodopt(options.aligncmp)
	options.ralign=parsemodopt(options.ralign)
	options.raligncmp=parsemodopt(options.raligncmp)
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
	if (options.lowmem):
		rimg=EMData()
	else:
		rimages = EMData.read_images(args[1],range(*rrange))
	
	for r in range(*rrange):
		if options.verbose: 
			print "%d/%d\r"%(r,rrange[1]),
			sys.stdout.flush()
			
		if ( options.lowmem ):
			rimg.read_image(args[1],r)
		else:
			rimg = rimages[r]
		
		row=cmponetomany(cimgs,rimg,options.align,options.aligncmp,options.cmp, options.ralign, options.raligncmp)
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
	
def cmponetomany(reflist,target,align=None,alicmp=("dot",{}),cmp=("dot",{}), ralign=None, alircmp=("dot",{})):
	"""Compares one image (target) to a list of many images (reflist). Returns """
	
	ret=[None for i in reflist]
	for i,r in enumerate(reflist):
		if align[0] :
			ta=target.align(align[0],r,align[1],alicmp[0],alicmp[1])
			#ta.debug_print_params()
			
			if ralign[0]:
				ralign[1]["az"] = ta.get_attr_default("align.az",0)-1
				ralign[1]["dx"] = ta.get_attr_default("align.dx",0)-1
				ralign[1]["dy"] = ta.get_attr_default("align.dy",0)-1
				refineparms[1]["mode"] = 0
				refineparms[1]["stepx"] = 2
				refineparms[1]["stepy"] = 2
				refineparms[1]["stepaz"] = 5
				ta = target.align(ralign[0],r,ralign[1],alircmp[0],alircmp[1])
				
			ret[i]=(ta.cmp(cmp[0],r,cmp[1]),ta.get_attr_default("align.dx",0),ta.get_attr_default("align.dy",0),ta.get_attr_default("align.az",0))
		else :
			ret[i]=(target.cmp(cmp[0],r,cmp[1]),0,0,0)
		
	return ret

def check(options,verbose):
	
	error = False
	if ( options.nofilecheck == False ):
		if not os.path.exists(options.datafile):
			if verbose:
				print "Error: the file expected to contain the particle images (%s) does not exist." %(options.reffile)
			error = True
		if not os.path.exists(options.reffile):
			if verbose:
				print "Error: the file expected to contain the projection images (%s) does not exist." %(options.reffile)
			error = True
		
		if ( os.path.exists(options.datafile) and os.path.exists(options.reffile) ):
			(xsize, ysize ) = gimme_image_dimensions2D(options.datafile);
			(pxsize, pysize ) = gimme_image_dimensions2D(options.reffile);
			if ( xsize != pxsize ):
				if verbose:
					print "Error - the (x) dimension of the reference images %d does not match that of the particle data %d" %(xsize,pxsize)
				error = True
			elif ( ysize != pysize ):
				if verbose:
					print "Error - the (y) dimension of the reference images %d does not match that of the particle data %d" %(ysize,pysize)
				error = True
		
		if os.path.exists(options.outfile):
			if ( not options.force):
				if verbose:
					print "Error: File %s exists, will not write over - specify the force option" %options.outfile
				error = True
	
	if (options.cmp == None or options.cmp == ""):
		if verbose:
			print "Error: the --cmp argument must be specified"
		error = True
	else:
		if ( check_eman2_type(options.cmp,Cmps,"Comparitor") == False ):
			error = True
	
	if (options.saveali):
		if   (options.align == None or options.align == ""):
			if verbose:
				print "Error: the --align argument must be specified if --saveali is specificed"
			error = True
		else:
			if ( check_eman2_type(options.align, Aligners,"Aligner") == False ):
				error = True
		
		if ( (options.aligncmp == None or options.aligncmp == "") and options.saveali):
			if verbose:
				print "Error: the --aligncmp argument must be specified if --saveali is specificed"
			error = True
		else:
			if ( check_eman2_type(options.aligncmp,Cmps,"Comparitor") == False ):
				error = True
		
	
		if ( options.ralign != None and options.ralign != ""):
			
			if ( check_eman2_type(options.ralign,Aligners,"Aligner") == False ):
				error = True
			
			if ( options.raligncmp == None or options.raligncmp == ""):
				if verbose:
					print "Error: the --raligncmp argument must be specified if --ralign is specificed"
			else:
				if ( check_eman2_type(options.raligncmp,Cmps,"Comparitor") == False ):
					error = True
	
	return error
	
	
	
if __name__ == "__main__":
    main()
