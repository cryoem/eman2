#!/usr/bin/env python

#
# Author: Steven Ludtke, 04/10/2003 (sludtke@bcm.edu)
# Copyright (c) 2000-2006 Baylor College of Medicine
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

# e2aligntest.py  09/21/2004  Steven Ludtke
# This program is used to generate various alignment test images


from EMAN2 import *
from optparse import OptionParser
from math import *
import os
import sys
#import pdb
from bisect import insort


def main():
	progname = os.path.basename(sys.argv[0])
	usage = """%prog [options] input1 input2
	
Locates the best 'docking' locations for a small probe in a large target map."""

	parser = OptionParser(usage=usage,version=EMANVERSION)

	parser.add_option("--axes",type="string",help="String list 3 axes from xyzaqp",default="xaq")
	
	(options, args) = parser.parse_args()
	if len(args)<2 : parser.error("Input1 and Input2 required")

	axes=options.axes
	
	
#	pdb.set_trace()
	i1=EMData()
	i2=EMData()
	
	i1.read_image(args[0],0)
	i2.read_image(args[1],0)
	
	result=EMData()
	result.set_size(64,64,64)
	
	alt=0.
	az=0.
	phi=0.
	dx=0.
	dy=0.
	dz=0.
	v=[0,0,0]
	for v[2] in range(32):
		print v[2],"/32"
		for v[1] in range(32):
			for v[0] in range(32):
				if "x" in axes: dx=(v[axes.find("x")]-16)/2.0
				if "y" in axes: dy=(v[axes.find("y")]-16)/2.0
				if "z" in axes: dz=(v[axes.find("z")]-16)/2.0
				if "a" in axes: alt=(v[axes.find("a")]-16)*2*pi/180.0
				if "q" in axes: az =(v[axes.find("q")]-16)*2*pi/180.0
				if "p" in axes: phi=(v[axes.find("p")]-16)*2*pi/180.0
				
				i2a=i2.copy(0)
				i2a.rotate_translate(alt,az,phi,dx,dy,dz)
				
				dot=i1.cmp("dot",{"with":EMObject(i2a)})
				result.set_value_at(v[0],v[1],v[2],dot)
	
	result.update()
	result.write_image("result.mrc")
	
	
if __name__ == "__main__":  main()
