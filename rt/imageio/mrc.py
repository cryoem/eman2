#!/usr/bin/python

from EMAN2 import *
import sys

infile = sys.argv[1]

e = EMData()
e.read_image(infile)
e.set_attr_dict("MRC.label3", EMObject("liwei peng"))
e.write_image("a.mrc", 0, MRC)

'''
e2 = EMData()
e2.read_image("a.mrc")
d = e2.get_attr_dict()
nlabels = d["MRC.nlabels"]
print "nlabels = ", nlabels

for i in range(int(nlabels)):
    label = "MRC.nlabel" + str(i)
    print d[label]
    
'''
