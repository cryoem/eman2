#!/bin/env python

from EMAN2 import *


image1 = EMData()
image1.read_image(Util.get_debug_image("samesize1.mrc"))

image2 = EMData()
image2.read_image(Util.get_debug_image("samesize2.mrc"))

image3 = image1.calc_ccf(image2)
image3.write_image("result.mrc")
