#!/bin/env python

from EMAN2 import *
import os

image1 = EMData()
image1.read_image(os.environ['HOME'] + "/images/samesize1.mrc")

image2 = EMData()
image2.read_image(os.environ['HOME'] + "/images/samesize2.mrc")

image3 = image1.calc_ccf(image2)
image3.write_image("result.mrc")
