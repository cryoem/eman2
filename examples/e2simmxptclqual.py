#!/usr/bin/env python

#
# Author: Steven Ludtke, 03/16/2012 (sludtke@bcm.edu)
# Copyright (c) 2000- Baylor College of Medicine
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
import os
import time
from EMAN2 import *
from numpy import *

def main():
        progname = os.path.basename(sys.argv[0])
        usage = """e2simmxptclqual.py [options] <simmx file in> 
	Computes the average simmx score vector for each orientation, normalizes it, then uses it to compute per-particle projections which are hopefully representative of particle quality.
        
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

#       parser.add_option("--input",type=str,help="Similarity matrix to analyze",default=None)
#        parser.add_argument("--refine",type=str,default=None,help="Automatically get parameters for a refine directory")
        #parser.add_argument("--output",type=str,help="Output text file",default="zvssim.txt")
        #parser.add_argument("--refs",type=str,help="Reference images from the similarity matrix (projections)",default=None)
        #parser.add_argument("--inimgs",type=str,help="Input image file",default=None)
        #parser.add_argument("--outimgs",type=str,help="Output image file",default="imgs.hdf")
        #parser.add_argument("--filtimgs",type=str,help="A python expression using Z[n], Q[n] and N[n] for selecting specific particles to output. n is the 0 indexed number of the input file",default=None)
        #parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	parser.add_argument("--refs",type=str,help="Reference images from the similarity matrix (projections)",default=None)
        
	(options, args) = parser.parse_args()

	mx=EMData(args[0],0,True)
	nx=mx["nx"]
	ny=mx["ny"]

	bvecs={}

# read in projection Euler angles
	if options.refs:
		ORTs=[]
		for i in xrange(nx):	
			# this reads the header, gets the orientation, and reads it out EMAN style
			ort=EMData(options.refs,i,True)["xform.projection"]
			o=ort.get_rotation("eman")
			o["phi"]=0
			ort.set_rotation(o)
			ORTs.append(ort)
		print nx," projections read"


	print "Computing average unit vectors"
	# compile vector sums for each class
	for y in range(ny):
		im=EMData(args[0],0,False,Region(0,y,nx,1))
		N=im.calc_min_index()
		im.process_inplace("normalize")
		try: bvecs[N].add(im)
		except: bvecs[N]=im

	# normalize all vector sums
	for im in bvecs.values(): 
		im.process_inplace("normalize.unitlen")

	# Make an output image of vectors
	mx=EMData(nx,nx,1)
	mx.to_zero()
	for i in range(nx):
		try: mx.insert_clip(bvecs[i],(0,i))
		except: pass
		
	mx.write_image("simvec.hdf",0)
	print "Output mean quality vector per class in simvec.hdf"

	print "Particle quality file"
	# Output particle quality file
	out=file("simqual.txt","w")
	t=time.time()
	for y in xrange(ny):
		if time.time()-t>.2 :
			print " %d\r"%y,
			sys.stdout.flush()
			t=time.time()

		im=EMData(args[0],0,False,Region(0,y,nx,1))
		N=im.calc_min_index()
		best=-10.0,-10.0,-1
		for r in xrange(nx):
			try:
				c=-im.cmp("ccc",bvecs[r])
				if c>best[0]: best=[c,0,r]
			except: pass
		best[1]=-im.cmp("dot",bvecs[best[2]])
		
		Nq=-im.cmp("ccc",bvecs[N])
		Nqd=-im.cmp("dot",bvecs[N])

		out.write("%d\t%d\t%1.4f\t%1.4g\t%d\t%1.4f\t%1.4g"%(y,N,Nq,Nqd,best[2],best[0],best[1]))
	
		if options.refs:
			# This is the amount of rotation required to move from one best Euler to the other, ignoring in-plane
			angdiff=fabs((ORTs[N]*ORTs[best[2]].inverse()).get_rotation("spin")["Omega"])
			out.write("\t%1.3g"%angdiff)
		out.write("\n")
	
	print "Output particle quality file simqual.txt"

if __name__ == "__main__":  main()

