#!/usr/bin/env python

from EMAN2 import *
import testlib
import sys

filternames = Filters.get_list()
assert(len(filternames) > 0)
print "number of filters = ", len(filternames)

e = EMData()
e.read_image(TestUtil.get_debug_image("search.dm3"))

fnum = 1000
f1 = Filters.get("Binarize", {'value': fnum})

new_params = f1.get_params()
assert(float(new_params["value"]) == fnum)

f1.process(e)

testlib.check_emdata(e, sys.argv[0])

