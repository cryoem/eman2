#!/usr/bin/python
# usage: filter1.py inputfile outputfile

from EMAN2 import *
import os

testimg = os.environ['HOME'] + "/images/groel3d.mrc"

e = EMData()
e.read_image(testimg)
e.filter("ValueSqrt")
e.filter("RangeThreshold", {"low" : 5.2, "high" : 10})
e.write_image("groel3df.mrc", 0, MRC)
