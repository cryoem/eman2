#!/bin/env python

from EMAN2 import *
import os

img1 = os.environ['HOME'] + "/images/3d86_1.mrc"
img2 = os.environ['HOME'] + "/images/3d86_2.mrc"

a=EMData()
b=EMData()
a.read_image(img1)
b.read_image(img2)
c=a.calc_ccf(b,1)
c.write_image("test_ccf2.mrc")
