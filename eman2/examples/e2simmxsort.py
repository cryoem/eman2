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
from EMAN2 import *
from numpy import *

def main():
        progname = os.path.basename(sys.argv[0])
        usage = """e2simmxsort.py [options] <simmx file in> <simmx file out>
Sorts a similarity matrix by classified particle to compare patterns.
        
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
        
	(options, args) = parser.parse_args()

	mx=EMData(args[0],0) 		# read in whole simmx
	nx=mx["nx"]
	ny=mx["ny"]

	# Make a list of all of the best class numbers for the images in the simmx
	Ns=[]
	for y in range(ny):
#		im=EMData(args[0],0,False,Region(0,y,nx,1))
		im=mx.get_clip(Region(0,y,nx,1))
		N=im.calc_min_index()
		Ns.append(N)

	Ns=array(Ns)
	Nsord=Ns.argsort()

	# generate the new sorted matrix
	last=0
	for i,j in enumerate(Nsord):
		im=EMData(args[0],0,False,Region(0,int(j),nx,1))
#		im.process_inplace("normalize")
		mx.insert_clip(im,(0,i))

	# chop the matrix up into per-class pieces
	CNs=bincount(Ns)	# number of particles in each class
	
	sm=0
	for i,j in enumerate(CNs):
		if j==0 : continue
		mx2=mx.get_clip(Region(0,int(sm),mx["nx"],int(j)))
		sm+=j
		mx2.write_image(args[1],i)
		print i,sm,j

if __name__ == "__main__":  main()

