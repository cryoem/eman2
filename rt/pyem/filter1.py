#!/usr/bin/env python

from EMAN2 import *
import os

filters = FilterFactory.instance()
filternames = filters.get_list()
print "all filters: ", filternames

e = EMData()
e.read_image(os.environ['HOME'] + "/images/search.dm3")

params = {'value': EMObject(int(1000))}
f1 = filters.get("Binarize", params)

new_params = f1.get_params()
print "Params in C++: ", new_params

f1.process(e)

e.write_image("search_f1.mrc", 0, MRC)
