#!/usr/bin/python
# usage: filter2.py inputfile outputfile

from EMAN2 import *
import sys

e = EMData()
e.read_image(sys.argv[1])

m = EMData()
m = e.copy(False, False)

e.filter("eman1.MaskNormalize", {"mask" : EMObject(m), "no_sigma" : 1})
e.write_image(sys.argv[2], 0, MRC)
