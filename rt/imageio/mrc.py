#!/usr/bin/python

from EMAN2 import *
import sys

infile = sys.argv[1]

e = EMData()
e.read_image(infile)

d = e.get_attr_dict()
nlabels = int(d["MRC.nlabels"])
print "nlabels = ", nlabels

label = "MRC.label" + str(nlabels)
e.set_attr_dict(label, EMObject("liwei peng"))
nlabels = nlabels + 1
e.set_attr_dict("MRC.nlabels", EMObject(nlabels))

e.write_image("a.mrc", 0, MRC)
