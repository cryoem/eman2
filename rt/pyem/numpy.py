#!/bin/env python

from EMAN2 import *
import os

e = EMData()
e.read_image(os.environ['HOME'] + "/images/tablet.mrc")

n = 500
ny = e.get_ysize()

a = Wrapper.em2numpy2(e)

for i in range(500):
	print e.get_value_at(n, i), a[i][n]


print a.shape
