#!/bin/env python
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
	usage = """Usage: %prog [options] input1 input2
	
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
				
				dot=i1.cmp("Dot",{"with":EMObject(i2a)})
				result.set_value_at(v[0],v[1],v[2],dot)
	
	result.update()
	result.write_image("result.mrc")
	
	
if __name__ == "__main__":  main()
