#!/usr/bin/env python

from EMAN2 import *

filternames = Filters.get_list()
print "all filters: ", filternames

e = EMData()
e.read_image(Util.get_debug_image("search.dm3")

f1 = Filters.get("Binarize", {'value': 1000})

new_params = f1.get_params()
print "Params in C++: value = ", float(new_params["value"])

f1.process(e)

e.write_image("search_f1.mrc", 0, MRC)
