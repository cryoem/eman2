#!/bin/env python

from EMAN2 import *
import math
import sys
import testlib

volume = EMData()
volume.read_image(TestUtil.get_debug_image("groel3d.mrc"))
pi = math.pi

proj = volume.project("Standard", { "alt" : pi/3, "az" : pi/5, "phi" : 1})
assert(proj.get_xsize() == 100)
assert(proj.get_ysize() == 100)
assert(proj.get_zsize() == 1)

testlib.check_emdata(proj, sys.argv[0])

