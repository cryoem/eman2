#!/bin/env python

from EMAN2 import *


e = EMData()
e.read_image(TestUtil.get_debug_image("search.dm3"))
pixels = e.calc_highest_locations(1200)
assert(len(pixels) == 612502)

p = pixels[0]
assert (p == Pixel(776,677,0, 1201))

