#!/bin/env python

from EMAN2 import *
from testlib import *

az = -0.60170830102
alt = 1.45232928554
phi = 0

t = Transform(EULER_EMAN, alt, az, phi)
rot = t.get_rotation(EULER_EMAN)

assertfloat(az, float(rot["az"]))
assertfloat(alt, float(rot["alt"]))
assertfloat(phi, float(rot["phi"]))
