#!/bin/env python

from EMAN2 import *


az = -0.60170830102
alt = 1.45232928554
phi = 0

rot = Rotation(alt, az, phi, Rotation.EulerType.EMAN)
print az, rot.eman_az()
print alt, rot.eman_alt()
print phi, rot.eman_phi()
