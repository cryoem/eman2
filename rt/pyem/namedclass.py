#!/bin/env python

from EMAN2 import *


imagename = TestUtil.get_debug_image("monomer.mrc")
x = EMData()
x.read_image(imagename)
y=x.get_rotated_clip([16,16,16],Rotation(0,0,0,Rotation.Type.EMAN),[16,16,16],1.0)
y2=x.get_rotated_clip((16,16,16),Rotation(0,0,0,Rotation.Type.EMAN),(16,16,16),1.0)
