#!/bin/env python

from EMAN2 import *
import math

volume = EMData()
volume.read_image(Util.get_debug_image("groel3d.mrc"))
pi = math.pi

proj = volume.project("Standard", { "alt" : pi/3, "az" : pi/5, "phi" : pi/7,})
print proj.get_xsize(), proj.get_ysize(), proj.get_zsize()
proj.write_image("proj.hed", 0, IMAGIC)
