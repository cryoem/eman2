#!/bin/env python

import os
from EMAN2 import *

imagename = os.environ['HOME'] + "/images/monomer.mrc"

a=EMData()
a.read_image(imagename)
b=a.get_rotated_clip([24,24,24],Rotation(0,0,0,Rotation.Type.EMAN),[32,32,32],1.0)
b.write_image("z.mrc")
