#!/bin/env python

from EMAN2 import *
import os

imgpath = os.environ["HOME"] + "/images/"
# 479x479x240
img1 = imgpath + "3f-avg.mrc"
# 100x100x100
img2 = imgpath + "3d.mrc"

a=EMData()
b=EMData()

a.read_image(img1,0)
b.read_image(img2,0)
b.read_image(img1,0)
