#!/bin/env python

from EMAN2 import *
import os

e = EMData()
e.read_image(os.environ['HOME'] + "/images/search.dm3")
pixels = e.calc_highest_locations(1200)

print "length=",len(pixels)
print "first one = ", pixels[0]


