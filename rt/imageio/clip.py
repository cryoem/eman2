#!/usr/bin/python

from EMAN2 import *

d = EMData()

d.read_image("GroEL-20031117173916.img",0)
d.write_image("x.mrc", 0, MRC)

oversample = 2
d2 = d.get_clip(Region(d.get_xsize()/2*(1-oversample),
                       d.get_ysize()/2*(1-oversample),
                       d.get_xsize()*oversample,
                       d.get_ysize()*oversample) )

d2.write_image("x2.mrc", 0, MRC)
