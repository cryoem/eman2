#!/usr/bin/python
# usage: align.py referencefile [inputfile outputfile]

from EMAN2 import *
import sys

if (len(sys.argv) != 1 and len(sys.argv) != 4):
	print "usage: align.py [referencefile inputfile outputfile]"
	sys.exit(1)

if len(sys.argv) == 1:
	reffile = TestUtil.get_debug_image("samesize1.mrc")
	inputfile = TestUtil.get_debug_image("samesize2.mrc")
	outputfile = "align.mrc"
	
else:
	reffile = sys.argv[1]
	inputfile = sys.argv[2]
	outputfile = sys.argv[3]
	
ref = EMData()
ref.read_image(reffile)

e = EMData()
e.read_image(inputfile)

e.align("rtf_best", {"to": ref, "maxshift": 8})

e.write_image(outputfile)

