#!/bin/env python

from EMAN2 import *
import os
import Numeric

n = 100

l = range(2*n*n)

a = Numeric.reshape(Numeric.array(l, Numeric.Float32), (2*n, n))

print a.shape
print a.typecode()

e = EMData()
Wrapper.numpy2em(a, e)
e.write_image("numpy.mrc")

for i in range(n):
	assert(e.get_value_at(i, 0) == i)
	
	
