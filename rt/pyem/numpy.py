#!/bin/env python

from EMAN2 import *

e = EMData()
e.read_image(Util.get_debug_image("tablet.mrc"))
nx = e.get_xsize()
ny = e.get_ysize()

a = Wrapper.em2numpy(e)
n = ny/2

for i in range(nx):
	assert(e.get_value_at(i, n) == a[n][i])

print a.shape
print a.typecode()

for x in range(nx):
	for y in range(n):
		a[y][x] = 0

e.write_image("test.mrc")

e2 = EMData()
Wrapper.numpy2em(a, e2)
e2.write_image("test2.mrc")
