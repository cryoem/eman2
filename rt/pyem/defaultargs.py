#!/bin/env python

from EMAN2 import *
import sys
from testlib import *

image1 = EMData()
image1.read_image(TestUtil.get_debug_image("samesize1.mrc"))

image2 = EMData()
image2.read_image(TestUtil.get_debug_image("samesize2.mrc"))

image3 = image1.calc_ccf(image2)
check_emdata(image3, sys.argv[0])

