#!/bin/env python

from EMAN2 import *
import os

e = EMData()
e.read_image(os.environ['HOME'] + "/images/tablet.mrc")
nx = e.get_xsize()
ny = e.get_ysize()

a = Wrapper.em2numpy(e)
n = ny/2

for i in range(nx):
	print e.get_value_at(i, n), a[n][i]


print a.shape

for x in range(nx):
	for y in range(n):
		a[y][x] = 0

e.write_image("test.mrc")

