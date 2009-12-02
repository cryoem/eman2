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
from optparse import OptionParser
from math import *
import os
import sys

#a = EMUtil.ImageType.IMAGE_UNKNOWN

PROJ_FILE_ATTR = "projection_file" # this attribute important to e2simmxxplor
PART_FILE_ATTR = "particle_file" # this attribute important to e2simmxxplor


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] <c input> <r input> <output> <ref simmx> <stg1 refs> <stg1 simmx>
	Computes a similarity matrix between c-input (col - projections) and r-input (row - particles) stacks of 2-D images. Unlike
	e2simmx.py, this will perform classification in two stages, first a coarse classification, then a local classification. Particle
	association for coarse classification, however, is not assigned based on Euler angle, but rather on mutual similarity in a subsampled
	reference-image self-classification. Output is the same as e2simmx, with coarsely sampled results inserted in unsampled local regions.
	When used for classification, c input is the references and r input are the particles."""
	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--align",type="string",help="The name of an 'aligner' to use prior to comparing the images", default=None)
	parser.add_option("--aligncmp",type="string",help="Name of the aligner along with its construction arguments",default="dot")
	parser.add_option("--ralign",type="string",help="The name and parameters of the second stage aligner which refines the results of the first alignment", default=None)
	parser.add_option("--raligncmp",type="string",help="The name and parameters of the comparitor used by the second stage aligner. Default is dot.",default="dot")
	parser.add_option("--cmp",type="string",help="The name of a 'cmp' to be used in comparing the aligned images", default="dot:normalize=1")
	parser.add_option("--mask",type="string",help="File containing a single mask image to apply before similarity comparison",default=None)
	parser.add_option("--saveali",action="store_true",help="Save alignment values, output is c x r x 4 instead of c x r x 1",default=False)
	parser.add_option("--verbose","-v",type="int",help="Verbose display during run",default=0)
#	parser.add_option("--lowmem",action="store_true",help="prevent the bulk reading of the reference images - this will save memory but potentially increase CPU time",default=False)
	parser.add_option("--exclude", type="string",default=None,help="The named file should contain a set of integers, each representing an image from the input file to exclude. Matrix elements will still be created, but will be zeroed.")
	parser.add_option("--shrink", type="int",default=None,help="Optionally shrink the input particles by an integer amount prior to computing similarity scores. This will speed the process up but may change classifications.")
	parser.add_option("--shrinks1", type="int",default=None,help="Shrinking performed for first stage classification, default=2",default=2)
	parser.add_option("--parallel",type="string",help="Parallelism string",default=None)

	(options, args) = parser.parse_args()
	
	if len(args)<6 : parser.error("Please specify all filenames : <c input> <r input> <output> <ref simmx> <stg1 refs> <stg1 simmx>")
	
	E2n=E2init(sys.argv)
	

	clen=EMUtil.get_image_count(args[0])
	rlen=EMUtil.get_image_count(args[1])
	clen_stg1=sqrt(clen*3)

	print "%d references, using %d stage 1 averaged references"%(clen,clen_stg1)

	############### Step 1 - classify the reference images
	cmd="e2simmx.py %s %s %s --shrink=%d --align=rotate_translate_flip --aligncmp=dot --cmp=phase --saveali"%(args[0],args[0],args[3],options.shrinks1)
	if options.parallel!=None : cmd+=" --parallel="+options.parallel
	print "executing ",cmd
	os.system(cmd)
	
	# Go through the reference self-simmx and determine the most self-dissimilar set of references
	ref_simmx=EMData(args[3],0)
	centers=[0]		# start with the first (generally a top view) image regardless 
	for i in xrange(clen_stg1-1) :
		best=(0,-1)
		for j in xrange(clen):
			if j in centers : continue
			simsum=0
			for k in refs: 
				simsum+=ref_simmx[j,k]
			if best[1]<0 or simsum>best[0] : best=(simsum,j)
		centers.append(best[1])

	# now associate each reference with the closest 3 centers
	centers=[[i] for i in centers]		# each center becomes a list to start the process
	for i in xrange(clen):
		quals=[(ref_simmx[i,k[0]],j) for j,k in enumerate(centers)]
		quals.sort()
		for j in xrange(3):

	E2progress(E2n,float(r-rrange[0])/(rrange[1]-rrange[0]))
	
	E2end(E2n)
	

if __name__ == "__main__":
    main()
