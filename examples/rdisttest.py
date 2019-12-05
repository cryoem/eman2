#!/usr/bin/env python
from __future__ import division

### This program computes the maximum radial distribution plot and characterizes it

from past.utils import old_div
from builtins import range
from EMAN2 import *
from sys import argv

im=EMData(argv[1]+"/threed_even_unmasked.hdf",0)

for f in range(6):
	aa=im.process("filter.lowpass.gauss",{"cutoff_freq":.08})
	if f : aa=aa.process("math.gausskernelfix",{"gauss_width":float(f)})
	md=aa.calc_radial_dist(old_div(aa["nx"],2),0,1,1)
	out=open("corm{:02d}.txt".format(f),"w")
	for x,y in enumerate(md):
		out.write("{}\t{}\n".format(x,y))

#im=EMData(argv[1]+"/threed_odd_unmasked.hdf",0)
#aa=im.process("filter.lowpass.gauss",{"cutoff_freq":.08})
#md2=aa.calc_radial_dist(aa["nx"]/2,0,1,3)

#plot(md,md2)
