#!/bin/env python

from EMAN2 import *

m1 = Transform()
v1 = Vec3f(1.1, 2.2, 3.3)
v2 = v1 * m1
print v2.as_list()

