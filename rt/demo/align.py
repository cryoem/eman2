#!/usr/bin/python
# usage: align.py referencefile inputfile outputfile

from EMAN2 import *
import sys

ref = EMData()
ref.read_image(sys.argv[1])

e = EMData()
e.read_image(sys.argv[2])

e.align("RTFBest", Dict("with", ref, "maxshift", 8))
e.write_image(sys.argv[3], 0, MRC)
