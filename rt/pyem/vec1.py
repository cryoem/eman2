#!/bin/env python

from EMAN2 import *

v = Vec3i(1,2,3)
print v.get_as_list()

m = Matrix3f()
print m.get_as_list()

e = EMData()
e.set_size(32,32,1)
intp = e.calc_min_location()

print intp
print type(intp)

e2 = EMData()
e2.set_size(12,12,1)
e.insert_clip(e2, [1,1,0])
