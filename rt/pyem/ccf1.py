#!/bin/env python

from EMAN import *


img1 = Util.get_debug_image("3d86_1.mrc")
img2 = Util.get_debug_image("3d86_2.mrc")

a=EMData()
b=EMData()
a.readImage(img1)
b.readImage(img2)

c=a.calcCCF(b,1)
c.writeImage("test_ccf1.mrc")
