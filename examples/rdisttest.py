#!/usr/bin/env python

### This program computes the maximum radial distribution plot and characterizes it

from EMAN2 import *
from sys import argv

im=EMData(argv[1]+"/threed_even_unmasked.hdf",0)

for f in xrange(15):
	aa=im.process("filter.lowpass.gauss",{"cutoff_freq":.08})
	if f : aa=aa.process("math.gausskernelfix",{"gauss_width":float(f)})
	md=aa.calc_radial_dist(aa["nx"]/2,0,1,3)
	out=open("cor{:02d}.txt".format(f),"w")
	for x,y in enumerate(md):
		out.write("{}\t{}\n".format(x,y))

#im=EMData(argv[1]+"/threed_odd_unmasked.hdf",0)
#aa=im.process("filter.lowpass.gauss",{"cutoff_freq":.08})
#md2=aa.calc_radial_dist(aa["nx"]/2,0,1,3)

#plot(md,md2)
