#!/usr/bin/env python

#
# Author: Steven Ludtke, 10/21/2011 (sludtke@bcm.edu)
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


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] input1 input2 ...
	
Reads a full similarity matrix. Computes a Z score for each particle and plots vs actual
similarity score."""

	parser = OptionParser(usage=usage,version=EMANVERSION)

#	parser.add_option("--input",type="string",help="Similarity matrix to analyze",default=None)
	parser.add_option("--output",type="string",help="Output text file",default="zvssim.txt")
	parser.add_option("--refs",type="string",help="Reference images from the similarity matrix (projections)",default=None)
	parser.add_option("--inimgs",type="string",help="Input image file",default=None)
	parser.add_option("--outimgs",type="string",help="Output image file",default="imgs.hdf")
	parser.add_option("--filtimgs",type="string",help="A python expression using Z[n], Q[n] and N[n] for selecting specific particles to output. n is the 0 indexed number of the input file",default=None)
	
	(options, args) = parser.parse_args()
	if len(args)<1 : 
		print "Please specify input files"
		sys.exit(1)

	tmp=EMData(args[0],0,True)		# header info from first input, assume others are same
	nx=tmp["nx"]
	ny=tmp["ny"]

	# read in projection Euler angles
	if options.refs:
		ALTs=[]
		AZs=[]
		for i in xrange(nx):	
			# this reads the header, gets the orientation, and reads it out EMAN style
			ort=EMData(options.refs,i,True)["xform.projection"].get_rotation("eman")
			ALTs.append(ort["alt"])
			AZs.append(ort["az"])
		print nx," projections read"

	out=file(options.output,"w")

	# We read one line of the simmx at a time. The line represents all values
	# for a single particle
	for y in xrange(ny):
		Qs=[]		# Quality
		Zs=[]		# Z score for classification
		Ns=[]		# classified best class
		for cm in args:
			im=EMData(cm,0,False,Region(0,y,nx,1))
	
			Z=(im["mean"]-im["minimum"])/im["sigma"]
			Zs.append(Z)
			Q=im["minimum"]
			Qs.append(Q)
			N=im.calc_min_index()
			Ns.append(N)

		for q in Qs : out.write("%1.4g\t"%q)
		for z in Zs : out.write("%1.4g\t"%z)
		for n in Ns : out.write("%d\t"%n)

		# if refs were provided we also write out Euler angle columns
		if options.refs :
			for n in Ns : out.write("%1.5g\t"%ALTs[n])
			for n in Ns : out.write("%1.5g\t"%AZs[n])

		out.write("\n")

		if options.filtimgs!=None and options.inimgs!=None and eval(options.filtimgs):
			EMData(options.inimgs,y).write_image(options.outimgs,-1)
	
if __name__ == "__main__":  main()
