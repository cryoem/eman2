#!/usr/bin/env python

#
# Author: Philip Baldwin, 9/12/2007 (woolford@bcm.edu)
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

from math import *
import os
import sys
from EMAN2db import db_check_dict
from EMAN2 import *

def main():
	
	progname = os.path.basename(sys.argv[0])
	usage = """prog <simMatrixIn> <classMatrixOut> [options]

	This program analyzes a similarity matrix as produced by e2simmx.py and produces a classification matrix mapping class
	membership for each particle, which can in-turn be used with e2classaverage.py.

	Takes a similarity matrix that has been created between reprojections-input (col) and particles-input (row) stacks of 2-D images.
	Typically these have been created via e2simmx such as

	e2simmx.py proj.hed part.hed simMatrix.hed --saveali --align=rotate_translate:maxshift=5

	The similarity matrix is a stack of 1 or 5 images: similarity, dx,dy,dalpha,mirror.
	The output is 6 images: class, weight, dx, dy, dalpha, mirror).
	See the wiki for more complete documentation of the files.
	"""
	
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)
	parser.add_argument("--sep", type=int, help="The number of classes a particle can contribute towards (default is 1)", default=1)
	parser.add_argument("--simvec",action="store_true",help="Instead of using the class for the peak value, uses the pattern of similarities for each orientation for assignment.",default=False)
	parser.add_argument("--force", "-f",dest="force",default=False, action="store_true",help="Force overwrite the output file if it exists")
	parser.add_argument("--verbose", "-v", dest="verbose", action="store", metavar="n",type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
	parser.add_argument("--nofilecheck",action="store_true",help="Turns file checking off in the check functionality - used by e2refine.py.",default=False)
	parser.add_argument("--check","-c",action="store_true",help="Performs a command line argument check only.",default=False)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)

	(options, args) = parser.parse_args()
	
	if options.nofilecheck: options.check = True

	if (len(args)<2 ): parser.error("Input and output files required")
	
	if (options.check): 
		options.verbose = 9 # turn verbose on if the user is only checking...
		
	if ( options.nofilecheck == False ):
		options.simmxfile = args[0]
		options.outfile = args[1]
		
	error = check(options,True)
	
	if (options.verbose > 0):
		if (error):
			print "e2classify.py command line arguments test.... FAILED"
		else:
			print "e2classify.py command line arguments test.... PASSED"
		
	if error : exit(1)
	if options.check: exit(0)

	
	E2n=E2init(sys.argv, options.ppid)

	if os.path.exists(args[1]):
		if (options.force):
			remove_file(args[1])

	num_sim =  EMUtil.get_image_count(args[0])
	if (num_sim < 5):
		print "Warning, similarity matrix did not contain alignments, neither will the classification matrix"

	tmp=EMData(args[0],0,True)
	nptcl=tmp["ny"]
	nref=tmp["nx"]
	if num_sim==5:
		clsmx=(EMData(options.sep,nptcl),EMData(options.sep,nptcl),EMData(options.sep,nptcl),EMData(options.sep,nptcl),EMData(options.sep,nptcl),EMData(options.sep,nptcl))
	elif num_sim==6:
		clsmx=(EMData(options.sep,nptcl),EMData(options.sep,nptcl),EMData(options.sep,nptcl),EMData(options.sep,nptcl),EMData(options.sep,nptcl),EMData(options.sep,nptcl),EMData(options.sep,nptcl))		
	else :
		clsmx=(EMData(options.sep,nptcl),EMData(options.sep,nptcl))

	# preparation for the simvec option. This finds all particles with a peak value corresponding to a particular orientation
	# and then computes an average similarity vector (across all references).
	if options.simvec:
		print "Computing average unit vectors"
		
		bvecs={}
		# compile vector sums for each class
		for y in range(nptcl):
			im=EMData(args[0],0,False,Region(0,y,nref,1))
			N=im.calc_min_index()
			im.process_inplace("normalize")
			try: bvecs[N].add(im)
			except: bvecs[N]=im

		# normalize all vector sums
		for im in bvecs.values(): 
			im.process_inplace("normalize.unitlen")

		## Make an output image of vectors
		#mx=EMData(nx,nx,1)
		#mx.to_zero()
		#for i in range(nx):
			#try: mx.insert_clip(bvecs[i],(0,i))
			#except: pass
			
		#mx.write_image("simvec.hdf",0)
		
		
		
	simmx=(EMData(),EMData(),EMData(),EMData(),EMData(),EMData())
	for iptcl in xrange(nptcl):
		simmx[0].read_image(args[0],0,False,Region(0,iptcl,nref,1))
		if num_sim>=5 :
			simmx[1].read_image(args[0],1,False,Region(0,iptcl,nref,1))		#tx
			simmx[2].read_image(args[0],2,False,Region(0,iptcl,nref,1))		#ty
			simmx[3].read_image(args[0],3,False,Region(0,iptcl,nref,1))		#dalpha
			simmx[4].read_image(args[0],4,False,Region(0,iptcl,nref,1))		#mirror
			try:
				simmx[5].read_image(args[0],5,False,Region(0,iptcl,nref,1))		#scale
			except: pass
		
		# We replace simmx[0] with a new version computed via average vectors
		if options.simvec:
			newmx=simmx[0].copy()
			for i in xrange(newmx["nx"]):
				try: 
#					simmx[0][i]=newmx.cmp("ccc",bvecs[i])
					simmx[0][i]=newmx.cmp("sqeuclidean",bvecs[i],{"normto":1})
				except: simmx[0][i]=100000.0		# bad value if we have no reference
		
		#hmmm, this code is pretty silly, but harmless, I suppose...
		maximum=simmx[0]["maximum"]
		for ic in xrange(options.sep):
			cls=simmx[0].calc_min_index()
			clsmx[0][ic,iptcl]=cls
			clsmx[1][ic,iptcl]=1.0		# weight
			if num_sim>=5:
				clsmx[2][ic,iptcl]=simmx[1][cls]
				clsmx[3][ic,iptcl]=simmx[2][cls]
				clsmx[4][ic,iptcl]=simmx[3][cls]
				clsmx[5][ic,iptcl]=simmx[4][cls]
				if num_sim>5 : clsmx[6][ic,iptcl]=simmx[5][cls]
			simmx[0][cls]=maximum		# this is so the next minimum search gets the next highest value
			
		E2end(E2n)

	print "Classification complete, writing classmx"
	clsmx[0].write_image(args[1],0)
	clsmx[1].write_image(args[1],1)
	if num_sim>=5 :
		clsmx[2].write_image(args[1],2)
		clsmx[3].write_image(args[1],3)
		clsmx[4].write_image(args[1],4)
		clsmx[5].write_image(args[1],5)
		if num_sim>5 : clsmx[6].write_image(args[1],6)
		
	

def check(options,verbose):
	error = False
	
	if (options.sep < 1):
		if verbose>0:
			print "Error: the --sep argument must be greater than zero, currently it is %d" %(options.sep)
		error = True
	

	
	if ( options.nofilecheck == False ):
		if os.path.exists(options.outfile):
			if (not options.force):
				if verbose>0:
					print "File %s exists, will not write over, exiting" %options.outfile
				error = True
		
		if not os.path.exists(options.simmxfile) and not db_check_dict(options.simmxfile):
			if verbose>0:
				print "Error: the similarity matrix file (%s) was not found, cannot run e2classify.py" %(options.simmxfile)
			error = True
		else:
			num_sim =  EMUtil.get_image_count(options.simmxfile)
			if (num_sim<5):
				if verbose>0:
					print "Error, the similarity matrix did not contain 5 images - be sure to use the --saveali argument when running e2simmx.py"
				error = True
	
	return error

if __name__ == "__main__":
    main()
