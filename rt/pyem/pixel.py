#!/bin/env python

from EMAN2 import *


e = EMData()
e.read_image(Util.get_debug_image("search.dm3")
pixels = e.calc_highest_locations(1200)

print "length=",len(pixels)
print "first one = ", pixels[0]


