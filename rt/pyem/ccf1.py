#!/bin/env python

from EMAN import *
import os

img1 = os.environ['HOME'] + "/images/3d86_1.mrc"
img2 = os.environ['HOME'] + "/images/3d86_2.mrc"

a=EMData()
b=EMData()
a.readImage(img1)
b.readImage(img2)

c=a.calcCCF(b,1)
c.writeImage("test_ccf1.mrc")
