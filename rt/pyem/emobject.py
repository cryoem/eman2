#!/bin/env python

from EMAN2 import *

n = 2
f = 3.5

e1 = EMObject(n)
e2 = EMObject(n)

assert(e1 == e2)

e3 = EMObject(f)
e4 = EMObject(f)

assert(e3 == e4)

assert(int(e1) == n)
assert(float(e3) == f)

n2 = 2.0
e5 = EMObject(float(n2))

assert(e1 == e5)
