#!/usr/bin/env python

from libpyEM import *

e = EMData()
e.read_image("/home/lpeng/images/tablet.mrc")
d = e.get_attr_dict()
print type(d)
print d

