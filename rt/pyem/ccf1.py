#!/bin/env python

from EMAN import *
import EMAN2

img1 = EMAN2.Util.get_debug_image("3d86_1.mrc")
img2 = EMAN2.Util.get_debug_image("3d86_2.mrc")

a=EMData()
b=EMData()
a.readImage(img1)
b.readImage(img2)

c=a.calcCCF(b,1)
c.writeImage("test_ccf1.mrc")
