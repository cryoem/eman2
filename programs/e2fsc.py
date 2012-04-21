#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/20/2012 (sludtke@bcm.edu)
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



from EMAN2 import *
from optparse import OptionParser
from math import *
import os
import sys
import time


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """e2fsc.py [options] input1 input2

Simple 2 volume FSCs can be computed with e2proc3d.py. This program is designed for more esoteric tasks, such as localized FSC calculations.

Two volumes with identical dimensions, and presumably identical orientations, should be provided at the command-line.
"""

	parser = EMArgumentParser(usage=usage,version=EMANVERSION)

#	parser.add_option("--input",type=str,help="Similarity matrix to analyze",default=None)
#	parser.add_argument("--refine",type=str,default=None,help="Automatically get parameters for a refine directory")
	parser.add_argument("--output",type=str,help="Output text file",default="zvssim.txt")
	#parser.add_argument("--refs",type=str,help="Reference images from the similarity matrix (projections)",default=None)
	#parser.add_argument("--inimgs",type=str,help="Input image file",default=None)
	#parser.add_argument("--outimgs",type=str,help="Output image file",default="imgs.hdf")
	#parser.add_argument("--filtimgs",type=str,help="A python expression using Z[n], Q[n] and N[n] for selecting specific particles to output. n is the 0 indexed number of the input file",default=None)
	parser.add_argument("--ppid", type=int, help="Set the PID of the parent process, used for cross platform PPID",default=-1)
	
	(options, args) = parser.parse_args()
		
	if len(args)<2 : 
		print "Please specify 2 input files"
		sys.exit(1)

	logid=E2init(sys.argv,options.ppid)

	v1=EMData(args[0],0)
	v2=EMData(args[1],0)
	
	nx,ny,nz=v1["nx"],v1["ny"],v1["nz"]
	print "%d x %d x %d"%(nx,ny,nz)
	
	if nx!=ny or nx!=nz : print "Warning: non-cubic volumes may produce unexpected results"
	
	# Create a centered Gaussian mask with a size ~1/10th of the box size
	cenmask=EMData(nx,ny,nz)
	cenmask.to_one()
	cenmask.process_inplace("mask.gaussian",{"outer_radius":nx/20})
	
	fys=[]
	for z in range(nz/20,nz,nz/10):
		for y in range(nz/20,ny,ny/10):
			for x in range(nz/20,nx,nx/10):
				print "%d, %d, %d :"%(x,y,z),
				
				# make a copy of the mask centered at the current position
				mask=cenmask.process("xform",{"transform":Transform({"tx":x-nx/2,"ty":y-ny/2,"tz":z-nz/2})})		

				v1m=v1*mask
				v2m=v2*mask
				if v1m["sigma_nonzero"]<.001 or v2m["sigma_nonzero"]<.001 :
					print "0 sigma"
					continue

				print v1m["sigma_nonzero"], v2m["sigma_nonzero"]
				
				fsc=v1m.calc_fourier_shell_correlation(v2m)
				fx=fsc[0:len(fsc)/3]
				fy=fsc[len(fsc)/3:len(fsc)*2/3]
				
				fys.append(fy)
				if isnan(fy[0]): print "NAN"
				else: 
					pass
					#plot((fx,fy))
	
	out=file("rslt.txt","w")
	for i,x in enumerate(fx):
		out.write( "%f\t"%x)
		for j in range(len(fy)):
			out.write( "%f\t"%fys[j][i])
			
		out.write("\n")

	E2end(logid)

if __name__ == "__main__":
        main()

