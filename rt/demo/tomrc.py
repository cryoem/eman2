#!/usr/bin/python
# usage: tomrc.py inputfile outputfile

from EMAN2 import *
import sys

e = EMData()
e.read_image(sys.argv[1])
e.write_image(sys.argv[2], 0, MRC)
