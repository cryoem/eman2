#!/bin/env python

# to test boost python's argument overloading
from EMAN2 import *

e = EMData()
e.set_size(10,20,2)
e.make_rotational_footprint()
