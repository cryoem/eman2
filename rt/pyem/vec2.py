#!/bin/env python

from EMAN2 import *

m1 = Transform()
v1 = Vec3f(1.5, 2.5, 3.5)
v2 = v1 * m1

assert (v2 == v1)

