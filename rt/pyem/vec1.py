#!/bin/env python

from EMAN2 import *

v = Vec3i(1,2,3)
vlist = v.as_list()
assert (vlist == [1,2,3])

e = EMData()
e.set_size(32,32,1)
intp = e.calc_min_location()

assert (intp == (0,0,0))
assert (type(intp) == type((1,2)))

e2 = EMData()
e2.set_size(12,12,1)
e.insert_clip(e2, [1,1,0])
e.insert_clip(e2, (1,1,0))

e2.translate(Vec3f(1,2,3))
e2.translate((1,2,3))
e2.translate([1,2,3])
