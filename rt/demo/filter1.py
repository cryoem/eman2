#!/usr/bin/python
# usage: filter1.py inputfile outputfile

from EMAN2 import *
import sys

e = EMData()
e.read_image(sys.argv[1])
e.filter("ValueSqrt")
e.filter("RangeMask", {"low" : 5.2, "high" : 10})
e.write_image(sys.argv[2], 0, MRC)
