#!/bin/env python

from EMAN2 import *

e = EMData()
all_imgs = e.read_images(TestUtil.get_debug_image("start.hed"))
for img in all_imgs:
	img.append_image("qq.hed")
	
