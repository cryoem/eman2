#!/bin/env python

from EMAN2 import *

e = EMData()
e.read_image(TestUtil.get_debug_image("tablet.mrc"))
d = e.get_attr_dict()

for k in d.keys():
    print k, d[k].to_str()

EMUtil.dump_dict(d)

nums = [1, 3]
images = EMData.read_images(TestUtil.get_debug_image("ali.hed"), nums)

nimg = len(images)
print "nimg = ", nimg

for d1 in images:
    EMUtil.dump_dict(d1.get_attr_dict())
    
