#!/bin/env python

from EMAN2 import *

# 479x479x240
img1 = Util.get_debug_image("3f-avg.mrc")
# 100x100x100
img2 = Util.get_debug_image("3d.mrc")

a=EMData()
b=EMData()

a.read_image(img1,0)
b.read_image(img2,0)
b.read_image(img1,0)
