#!/bin/env python

from EMAN2 import *
import os
import numarray

n = 100

l = range(2*n*n)

a = numarray.reshape(numarray.array(l, numarray.Float32), (2*n, n))

print a.shape
print a.typecode()

e = EMData()
Wrapper.numpy2em(a, e)
e.write_image("numpy.mrc")

for i in range(n):
	assert(e.get_value_at(i, 0) == i)
	
	
