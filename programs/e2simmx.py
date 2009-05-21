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


def opt_rectangular_subdivision(x,y,n):
	xparts = 1 # x partitions
	yparts = 1 # y partitions
	
	width = x
	height = y
	part_width = width
	part_height = height
	
	
	candidates = []
	
	min = None
	best = None
	for cn in xrange(1,n/2+1):
		rn = n/cn
		
		f = rn*x + cn*y
		print f,rn,cn
		if min == None or f < min:
			best = [rn,cn,y/cn,x/rn]
			min = f
	return best
#	for xdiv in xrange(1,width/2+1):
#		ydiv = n/xdiv
#		if ydiv == 0: ydiv = 1
#		total = xdiv*ydiv
#		if total > n-5 and total <= n:candidates.append([xdiv,ydiv])
#		
#	for ydiv in xrange(1,height/2+1):
#		xdiv = y/ydiv
#		if xdiv == 0: xdiv = 1
#		if total > n-5 and total <= n:candidates.append([xdiv,ydiv])
#		
#	
#	if len(candidates) == 0: raise RuntimeError
#	
#	min = None
#	solution = []
#	
#	for xdiv,ydiv in candidates:
#		nx = width/xdiv
#		ny = height/ydiv
#		total = nx + ny
#		if min == None or total < min:
#			solution = [nx,ny,xdiv*ydiv,total,width-nx*xdiv,height-ny*ydiv]
#			min = total
#			
#	return solution
	
	
#	while xparts*yparts < n:
#		if part_height > part_width:
#			yparts += 1 
#			part_height = height/yparts
#		else:
#			xparts += 1
#			part_width = width/xparts
#			
#		solutions.append([part_width,part_height,xparts*yparts])
			
	
	
	


class EMParallelSimMX:
	def __init__(self,options,args,logger=None):
		'''
		@param options the options produced by (options, args) = parser.parse_args()
		@param args the options produced by (options, args) = parser.parse_args()
		@param logger and EMAN2 logger, i.e. logger=E2init(sys.argv)
		assumes you have already called the check function.
		'''
		self.options = options
		self.args = args
		self.logger = logger
		
	
		from EMAN2PAR import EMTaskCustomer
		etc=EMTaskCustomer("dc:localhost:9990")
		self.num_cpus = etc.cpu_est()
		self.num_cpus = 32
	
	
	def __init_memory(self):
		self.clen=EMUtil.get_image_count(self.args[0])
		self.rlen=EMUtil.get_image_count(self.args[1])
	
	def __get_blocks(self):
		'''
		Gets the blocks that will be processed in parallel, these are essentially ranges
		'''
		total_jobs = self.num_cpus
		block_c = self.clen/total_jobs
		block_r = self.rlen/total_jobs
		
		residual_c = self.clen-block_c*total_jobs # residual left over by integer division
		residual_r = self.rlen-block_r*total_jobs # residual left over by integer division

		current_c = 0
		current_r = 0

		ranges = []
		for i in xrange(0,total_jobs):
			last_c = current_c + block_c
			if residual_c > 0:
				last_c += 1
				residual_c -= 1
			
			last_r = current_r + block_r
			if residual_r > 0:
				last_r += 1
				residual_r -= 1
			
			
			ranges.append([current_c,last_c,current_r,last_r])
			current_c = last_c
			current_r = last_r
		
		return ranges
	
	
	def execute(self):
		'''
		The main function to be called
		'''
		self.__init_memory()
		print self.__get_blocks()
		pass
	
		
		
	

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
	parser.add_option("--saveali",action="store_true",help="Save alignment values, output is c x r x 4 instead of c x r x 1",default=True)
	parser.add_option("--verbose","-v",type="int",help="Verbose display during run",default=0)
#	parser.add_option("--lowmem",action="store_true",help="prevent the bulk reading of the reference images - this will save memory but potentially increase CPU time",default=False)
	parser.add_option("--init",action="store_true",help="Initialize the output matrix file before performing 'range' calculations",default=False)
	parser.add_option("--force", "-f",dest="force",default=False, action="store_true",help="Force overwrite the output file if it exists")
	parser.add_option("--exclude", type="string",default=None,help="The named file should contain a set of integers, each representing an image from the input file to exclude. Matrix elements will still be created, but will be zeroed.")
	parser.add_option("--shrink", type="int",default=None,help="Optionally shrink the input particles by an integer amount prior to computing similarity scores. This will speed the process up.")
	parser.add_option("--nofilecheck",action="store_true",help="Turns file checking off in the check functionality - used by e2refine.py.",default=False)
	parser.add_option("--check","-c",action="store_true",help="Performs a command line argument check only.",default=False)
	parser.add_option("--parallel",type="string",help="Parallelism string",default=None)

	(options, args) = parser.parse_args()
	
	if len(args)<3 : parser.error("Input and output files required")
	
	if (options.check): options.verbose = True # turn verbose on if the user is only checking...
		
	if (options.nofilecheck == False):
		options.reffile=args[0]
		options.datafile = args[1]
		options.outfile = args[2]
	
	error = check(options,True)
	
	if (options.verbose):
		if (error):
			print "e2simmx.py command line arguments test.... FAILED"
		else:
			print "e2simmx.py command line arguments test.... PASSED"
		
	if error : exit(1)
	if options.check: exit(0)
	
	E2n=E2init(sys.argv)
	
	if options.parallel:
		parsimmx = EMParallelSimMX(options,args,E2n)
		parsimmx.execute()
		E2end(E2n)
		sys.exit(0)
		
	
	# just remove the file - if the user didn't specify force then the error should have been found in the check function
	if os.path.exists(options.outfile):
		if (options.force):
			remove_file(options.outfile)
	
	options.align=parsemodopt(options.align)
	options.aligncmp=parsemodopt(options.aligncmp)
	options.ralign=parsemodopt(options.ralign)
	options.raligncmp=parsemodopt(options.raligncmp)
	options.cmp=parsemodopt(options.cmp)

	if options.exclude: 
		try:
			excl=file(options.exclude,"r").readlines()
			excl=[int(i) for i in excl]
			excl=set(excl)
		except: print "Warning: exclude file failed"		# it's ok if this fails
		

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
	mxout=[EMData()]
	mxout[0].set_size(crange[1]-crange[0],rrange[1]-rrange[0],1)
	mxout[0].to_zero()
	if options.saveali : 
		mxout.append(mxout[0].copy()) # dx
		mxout.append(mxout[0].copy()) # dy
		mxout.append(mxout[0].copy()) # alpha (angle)
		mxout.append(mxout[0].copy()) # mirror
	if options.verbose: print "Computing Similarities"

	# Read all c images, then read and compare one r image at a time
	cimgs=EMData.read_images(args[0],range(*crange))
	if options.shrink != None: # the check function guarantees that shrink is an integer greater than 1
		#d = [ image.process("math.meanshrink",{"n":options.shrink}) for image in cimgs]
		#cimgs = d
		for image in cimgs:
			image.process_inplace("math.meanshrink",{"n":options.shrink})
	
#	if (options.lowmem):
	rimg=EMData()
#	else:
#		rimages = EMData.read_images(args[1],range(*rrange))
#		if options.shrink != None: # the check function guarantees that shrink is an integer greater than 1
#			#d = [ image.process("math.meanshrink",{"n":options.shrink}) for image in rimages]
#			#rimages = d
#			# I chose this way in the end for memory efficiency. There's probably a better way to do it
#			for image in rimages:
#				image.process_inplace("math.meanshrink",{"n":options.shrink})
		
	#dimages =  EMData.read_images(args[1],range(*rrange))
	#d = [ image.process_inplace("math.meanshrink",{"n":options.shrink}) for image in dimages]
	
	for r in range(*rrange):
		if options.exclude and r in excl : continue
		
		if options.verbose: 
			print "%d/%d\r"%(r,rrange[1]),
			sys.stdout.flush()
			
#		if ( options.lowmem ):
		rimg.read_image(args[1],r)
		if options.shrink != None: # the check function guarantees that shrink is an integer greater than 
			rimg.process_inplace("math.meanshrink",{"n":options.shrink})
#		else:
#			rimg = rimages[r]
		
		E2progress(E2n,float(r-rrange[0])/(rrange[1]-rrange[0]))
		row=cmponetomany(cimgs,rimg,options.align,options.aligncmp,options.cmp, options.ralign, options.raligncmp)
		for c,v in enumerate(row):
			mxout[0].set_value_at(c,r,0,v[0])
		
		if options.saveali :
			for c,v in enumerate(row):
#				print row
				#print "%f %f %f " %(v[1],v[2],v[3])
				mxout[1].set_value_at(c,r,0,v[1])
				mxout[2].set_value_at(c,r,0,v[2])
				mxout[3].set_value_at(c,r,0,v[3])
				mxout[4].set_value_at(c,r,0,v[4])
	
	if options.verbose : print"\nSimilarity computation complete"
	
	# write the results into the full-sized matrix
	if crange==[0,clen] and rrange==[0,rlen] :
		for i,j in enumerate(mxout) : j.write_image(args[2],i)
	else :
		for i,j in enumerate(mxout) : j.write_image(args[2],i,IMAGE_UNKNOWN,0,Region(crange[0],rrange[0],0,crange[1]-crange[0],rrange[1]-rrange[0],1))
	
	E2end(E2n)
	
def cmponetomany(reflist,target,align=None,alicmp=("dot",{}),cmp=("dot",{}), ralign=None, alircmp=("dot",{})):
	"""Compares one image (target) to a list of many images (reflist). Returns """
	
	ret=[None for i in reflist]
	for i,r in enumerate(reflist):
		if align[0] :
			ta=target.align(align[0],r,align[1],alicmp[0],alicmp[1])
			#ta.debug_print_params()
			
			if ralign and ralign[0]:
				ralign[1]["xform.align2d"] = ta.get_attr("xform.align2d")
				ta = target.align(ralign[0],r,ralign[1],alircmp[0],alircmp[1])
			
			t = ta.get_attr("xform.align2d")
			p = t.get_params("2d")
			
			ret[i]=(ta.cmp(cmp[0],r,cmp[1]),p["tx"],p["ty"],p["alpha"],p["mirror"])
		else :
			ret[i]=(target.cmp(cmp[0],r,cmp[1]),0,0,0,False)
		
	return ret

def check(options,verbose):
	
	error = False
	if ( options.nofilecheck == False ):
		if not os.path.exists(options.datafile) and not db_check_dict(options.datafile):
			if verbose:
				print "Error: the file expected to contain the particle images (%s) does not exist." %(options.reffile)
			error = True
		if not os.path.exists(options.reffile) and not db_check_dict(options.reffile):
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
	
	if (options.shrink != None):
		if options.shrink == 1:
			options.shrink = None # just leave it as None please
			print "Warning, setting shrink to 1 does nothing. If you don't want shrinking to occur just forget the shrink argument"
			
		if options.shrink <= 1:
			print "Error: shrink must be greater than 1"
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
