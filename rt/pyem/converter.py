#!/usr/bin/env python

from EMAN2 import *

e = EMData()
e.read_image("/home/lpeng/images/tablet.mrc")
d = e.get_attr_dict()

for k in d.keys():
    print k, d[k].to_str()

EMUtil.dump_dict(d)

    
nums = [1, 3]
images = EMData.read_images_by_index("/home/lpeng/images/ali.hed", nums)

nimg = len(images)
print "nimg = ", nimg

for d1 in images:
    EMUtil.dump_dict(d1.get_attr_dict())
    
