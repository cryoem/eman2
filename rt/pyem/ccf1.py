#!/bin/env python

from EMAN import *
import EMAN2
import testlib
import sys

img1 = EMAN2.TestUtil.get_debug_image("3d86_1.mrc")
img2 = EMAN2.TestUtil.get_debug_image("3d86_2.mrc")

a=EMData()
b=EMData()
a.readImage(img1)
b.readImage(img2)

c=a.calcCCF(b,1)
testlib.check_emdata(c, sys.argv[0])

