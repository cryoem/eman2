#!/bin/env python

from EMAN2 import *

m1 = Matrix3f()
v1 = Vec3f(1.1, 2.2, 3.3)
v2 = m1 * v1
print v1.get_as_list(), " * ", m1.get_as_list(), " = ", v2.get_as_list()

