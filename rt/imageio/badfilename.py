#!/usr/bin/python

from EMAN2 import *

img = EMData()
imagename = Util.get_debug_image("square.mrc")
img.read_image(imagename)
img.write_image("out1", 0, IMAGIC)
