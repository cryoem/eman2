#!/bin/env python

from EMAN2 import *

v1 = Vec3f(1,2,3)
v2 = Vec3f(2,4,6)

v1 = v1 + v2
v1list = v1.as_list()
assert(v1list == [3, 6, 9])
