#!/bin/env python

from EMAN2 import *
import os

testimg = os.environ['HOME'] + "/images/groel3d.mrc"

img1 = EMData()
img1.read_image(testimg)

img2 = EMData()
img2.read_image(testimg)

cmpscore = img1.cmp("Dot", {"with": EMObject(img2)})

print "cmp score = ", cmpscore

