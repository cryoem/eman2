#!/bin/env python

from EMAN2 import *

# 479x479x240
img1 = TestUtil.get_debug_image("groel2d.mrc")
# 100x100x100
img2 = TestUtil.get_debug_image("3d.mrc")

a=EMData()
b=EMData()

a.read_image(img1,0)
b.read_image(img2,0)
b.read_image(img1,0)
