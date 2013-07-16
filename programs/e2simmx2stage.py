#!/usr/bin/env python

#
# Author: Steven Ludtke, 12/01/2009 (sludtke@bcm.edu)
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

# e2simmx2stage.py  12/01/2009	Steven Ludtke
# This program computes a similarity matrix between two sets of images

from EMAN2 import *
from math import *
import os
import sys

#a = EMUtil.ImageType.IMAGE_UNKNOWN

PROJ_FILE_ATTR = "projection_file" # this attribute important to e2simmxxplor
PART_FILE_ATTR = "particle_file" # this attribute important to e2simmxxplor


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """prog [options] <c input> <r input> <output> <ref simmx> <stg1 refs> <stg1 simmx>
	Computes a similarity matrix between c-input (col - projections) and r-input (row - particles) stacks of 2-D images. Unlike
	e2simmx.py, this will perform classification in two stages, first a coarse classification, then a local classification. Particle
	association for coarse classification, however, is not assigned based on Euler angle, but rather on mutual similarity in a subsampled
	reference-image self-classification. Output is the same as e2simmx, with coarsely sampled results inserted in unsampled local regions.
	When used for classification, c input is the references and r input are the particles."""
	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

	parser.add_argument("--align",type=str,help="The name of an 'aligner' to use prior to comparing the images", default=None)
	parser.add_argument("--aligncmp",type=str,help="Name of the aligner along with its construction arguments",default="dot")
	parser.add_argument("--ralign",type=str,help="The name and parameters of the second stage aligner which refines the results of the first alignment", default=None)
	parser.add_argument("--raligncmp",type=str,help="The name and parameters of the comparitor used by the second stage aligner. Default is dot.",default="dot")
	parser.add_argument("--cmp",type=str,help="The name of a 'cmp' to be used in comparing the aligned images", default="dot:normalize=1")
	parser.add_argument("--mask",type=str,help="File containing a single mask image to apply before similarity comparison",default=None)
	parser.add_argument("--colmasks",type=str,help="File containing one mask for each column (projection) image, to be used when refining row (particle) image alignments.",default=None)
	parser.add_argument("--saveali",action="store_true",help="Save alignment values, output is c x r x 4 instead of c x r x 1",default=False)
	parser.add_argument("--prefilt",action="store_true",help="Filter each reference (c) to match the power spectrum of each particle (r) before alignment and comparison",default=False)
	parser.add_argument("--verbose", "-v", dest="verbose", type=int, default=0, help="verbose level [0-9], higner number means higher level of verboseness")
#	parser.add_argument("--lowmem",action="store_true",help="prevent the bulk reading of the reference images - this will save meclen,mory but potentially increase CPU time",default=False)
	parser.add_argument("--exclude", type=str,default=None,help="The named file should contain a set of integers, each representing an image from the input file to exclude. Matrix elements will still be created, but will be zeroed.")
	parser.add_argument("--shrink", type=int,default=None,help="Optionally shrink the input particles by an integer amount prior to computing similarity scores. This will speed the process up but may change classifications.")
	parser.add_argument("--shrinks1", type=int,help="Shrinking performed for first stage classification, default=2",default=2)
	parser.add_argument("--finalstage",action="store_true",help="Assume that existing preliminary particle classifications are correct, and only recompute final local orientations",default=False)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--parallel",type=str,help="Parallelism string",default="thread:1")
	parser.add_argument("--force", "-f",dest="force",default=True, action="store_true",help="Deprecated. Value ignored")

	(options, args) = parser.parse_args()

	if len(args)<6 : parser.error("Please specify all filenames : <c input> <r input> <output> <ref simmx> <stg1 refs> <stg1 simmx>")

	E2n=E2init(sys.argv,options.ppid)


	clen=EMUtil.get_image_count(args[0])
	rlen=EMUtil.get_image_count(args[1])
	clen_stg1=3*int(sqrt(clen))

	print "%d references, using %d stage 1 averaged references"%(clen,clen_stg1)

	if not options.finalstage :
		############### Step 1 - classify the reference images

		E2progress(E2n,0.001)
		# compute the reference self-similarity matrix
		cmd="e2simmx.py %s %s %s --shrink=%d --align=rotate_translate_flip --aligncmp=sqeuclidean:normto=1 --cmp=sqeuclidean --saveali --force --verbose=%d"%(args[0],args[0],args[3],options.shrinks1, options.verbose-1)
		if options.prefilt : cmd+=" --prefilt"
		if options.parallel!=None : cmd+=" --parallel="+options.parallel
		print "executing ",cmd
		launch_childprocess(cmd)

		E2progress(E2n,0.10)
		# Go through the reference self-simmx and determine the most self-dissimilar set of references
		print "Finding %d dissimilar classification centers"%clen_stg1
		ref_simmx=EMData(args[3],0)
		ref_orts=EMData.read_images(args[3],(1,2,3,4))
		centers=[0]		# start with the first (generally a top view) image regardless
		for i in xrange(clen_stg1-1) :
			best=(0,-1)
			for j in xrange(clen):
				if j in centers : continue

				#here we find the sum of the distances to this point
				#simsum=0
				#for k in centers:
					#simsum+=ref_simmx[j,k]

				# here we find the distance to the closest matching existing point
				# so we find the point with the largest minimum distance to its nearest neighbor
				simsum=1.0e30
				for k in centers:
					if ref_simmx[j,k]<simsum: simsum=ref_simmx[j,k]

				if best[1]<0 or simsum>best[0] : best=(simsum,j)
			centers.append(best[1])

	#	print centers
		# Re-sort references by similarity
		print "Sort references"
		for i in range(1,clen_stg1-1):
			for j in range(i+1,clen_stg1):
				if ref_simmx[centers[i-1],centers[i]]>ref_simmx[centers[i-1],centers[j]] : centers[i],centers[j]=centers[j],centers[i]

		# now associate each reference with the closest center
		print "Associating references with centers"
		classes=[[] for i in centers]	# each center becomes a list to start the process
		for i in xrange(clen):
			quals=[(ref_simmx[i,k],j) for j,k in enumerate(centers)]
			quals.sort()
#			for j in xrange(4): classes[quals[j][1]].append(i)		# we used to associate each reference with 3 closest centers
#			classes[quals[0][1]].append(i)							# now we just associate it with the closest one, but use multiple centers when searching
			for j in xrange(2): classes[quals[j][1]].append(i)		# bring this idea back again with 2 closest centers

		# now generate an averaged reference for each center
		print "Averaging each center"
		for ii,i in enumerate(classes):
	#		print "%d.  %d"%(ii,len(i)),i[:6]
			avg=EMData(args[0],i[0])
			for j in i[1:]:
				tmp=EMData(args[0],j)
	#			print ref_orts[0][i[0],j],"\t",ref_orts[1][i[0],j],"\t",ref_orts[2][i[0],j],"\t",int(ref_orts[3][i[0],j])
				xf=Transform({"type":"2d","tx":ref_orts[0][i[0],j],"ty":ref_orts[1][i[0],j],"alpha":ref_orts[2][i[0],j],"mirror":bool(ref_orts[3][i[0],j])})
				tmp.process_inplace("xform",{"transform":xf})
	#			tmp.write_image("testing/rcls.%03d.hdf"%ii,-1)
				avg.add(tmp)

			avg.mult(1.0/len(i))
			avg["class_ptcl_idxs"]=i
			avg["class_ptcl_src"]=args[0]
			avg.write_image(args[4],ii)

		E2progress(E2n,0.15)
		############### Step 2 - classify the particles against the averaged references
		print "First stage particle classification"
		cmd="e2simmx.py %s %s %s --shrink=%d --align=%s --aligncmp=%s  --cmp=%s  --saveali --force --verbose=%d"%(args[4],args[1],args[5],options.shrinks1,
			options.align,options.aligncmp,options.cmp, options.verbose-1)
		if options.prefilt : cmd+=" --prefilt"
		if options.ralign!=None : cmd+=" --ralign=%s --raligncmp=%s"%(options.ralign,options.raligncmp)
		if options.parallel!=None : cmd+=" --parallel="+options.parallel
		if options.exclude!=None : cmd+=" --exclude="+options.exclude
		print "executing ",cmd
		launch_childprocess(cmd)
	else :
		# reread classification info
		classes=[i["class_ptcl_idxs"] for i in EMData.read_images(args[4],None,True)]


	E2progress(E2n,0.60)
	############### Step 3 - classify particles against subset of original projections
	# Now we need to convert this small classification into a 'seed' for the large classification matrix for simplicity
	print "Seeding full classification matrix (%d x %d) -> %s"%(clen,rlen,args[2])
	mxstg1=EMData(args[5],0)
	mx=EMData(clen,rlen,1)
	mx.add(-1.0e38)	# a large negative value to be replaced later

	for ptcl in xrange(rlen):
		# find the best classes from the coarse search
		vals=[(mxstg1[cls1,ptcl],cls1) for cls1 in range(clen_stg1)]
		vals.sort()

		# then set the corresponding values in the best 4 stage 1 classes in the full matrix to 0
		for j in xrange(4):
			for i in classes[vals[j][1]]: mx[i,ptcl]=0.0

	mx.update()
	mx.write_image(args[2],0)
#	mx.write_image("bdb:refine_02#simmx_00_x",0)

	mx.to_zero()
	if options.saveali:
		for i in xrange(1,6):
			mx.write_image(args[2],i)		# seed alignment data with nothing
#			mx.write_image("bdb:refine_02#simmx_00_x",i)


	# the actual final classification
	cmd = "e2simmx.py %s %s %s -f --saveali --cmp=%s --align=%s --aligncmp=%s --fillzero --nofilecheck --force --verbose=%d"  %(args[0],args[1],args[2],options.cmp,options.align,options.aligncmp, options.verbose-1)
	if options.mask!=None : cmd += " --mask=%s"%options.mask
	if options.colmasks!=None : cmd += " --colmasks=%s"%options.colmasks

	if ( options.ralign != None ):
		cmd += " --ralign=%s --raligncmp=%s" %(options.ralign,options.raligncmp)

	if options.verbose>0:
		cmd += " --verbose=%d"%(options.verbose-1)

	if options.prefilt : cmd+=" --prefilt"
	if options.parallel: cmd += " --parallel=%s" %options.parallel

	#if (options.lowmem): e2simmxcmd += " --lowmem"

	if (options.shrink):
		cmd += " --shrink="+str(options.shrink)

	print "executing ",cmd
	launch_childprocess(cmd)

	E2progress(E2n,1.0)

#	E2progress(E2n,float(r-rrange[0])/(rrange[1]-rrange[0]))

	E2end(E2n)


if __name__ == "__main__":
    main()
