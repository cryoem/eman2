#!/usr/bin/python
# usage: project.py inputfile outputfile

from EMAN2 import *
import sys

e = EMData()
e.read_image(sys.argv[1])

p = e.project("Pawel")
p.write_image(sys.argv[2], 0, IMAGIC)
