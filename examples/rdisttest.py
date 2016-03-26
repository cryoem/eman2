#!/usr/bin/env python

### This program computes the maximum radial distribution plot and characterizes it

from EMAN2 import *
from sys import argv

im=EMData(argv[1],0)

aa=im.process("filter.lowpass.gauss",{"cutoff_freq":.05})
md=aa.calc_radial_dist(aa["nx"]/2,0,1,3)
plot(md)
