#!/bin/env python

from EMAN2 import *

x = 1
y = 2
z = 3
v = 2.5

p1 = Pixel(x,y,z,v)
p2 = Pixel(x,y,z,v)
assert(p1 == p2)

p3 = Pixel(1,2,3,2.3)
p4 = Pixel(11,2,3,2.3)
assert(p3 != p4)

assert(p1.x == x and p1.y == y and p1.z == z and p1.value == v)


