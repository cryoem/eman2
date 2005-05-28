#!/usr/bin/python

from EMAN2 import *
e = EMData()
e.set_size(100,100,1)
e.to_zero()
e.make_rotational_footprint()

