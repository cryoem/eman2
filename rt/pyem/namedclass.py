#!/bin/env python

from EMAN2 import *
import testlib
import sys

imagename = TestUtil.get_debug_image("monomer.mrc")
x = EMData()
x.read_image(imagename)
y=x.get_rotated_clip(Transform([16,16,16],EULER_EMAN,0,0,0),[16,16,16],1.0)
y2=x.get_rotated_clip(Transform((16,16,16),EULER_EMAN,0,0,0),(16,16,16),1.0)
testlib.check_emdata(y, sys.argv[0])
testlib.check_emdata(y2, sys.argv[0])

