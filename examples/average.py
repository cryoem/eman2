#!/usr/bin/env python

### This program will average a list of images or volumes passed in on the command-line
### average.py <vol1.mrc> <vol2.mrc> ...
### will produce average.hdf
### All inputs must have the same dimensions

from EMAN2 import *
from sys import argv

av=EMData(argv[1],0)
for i in argv[2:]:
	im=EMData(i,0)
	av.add(im)

av.mult(1.0/(len(argv)-1))
av.write_image("average.hdf",0)

