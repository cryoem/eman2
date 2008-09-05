#!/usr/bin/env python

from EMAN2 import *
from sparx import *
from math import sqrt
from sys import argv
from sys import exit

if len(argv) != 4:
    print "filt_by_rops input target output"
    exit(-1)

img_stack = argv[1]
img_dst = getImage( argv[2] )
out_stack = argv[2]

nimage = EMUtil.get_image_count( img_stack )

if nimage < 100:
    navg = nimage
else:
    navg = 100

img_src = getImage(img_stack, 0)
for i in xrange(1,navg):
    img_tmp = get_im(img_stack, i)
    img_src += img_tmp
img_src /= navg


rops_src = rops_table(img_src)
rops_dst = rops_table(img_dst)

assert len(rops_dst) == len(rops_src)

table = [0.0]*len(rops_dst)
for i in xrange( len(rops_dst) ):
    table[i] = sqrt( rops_dst[i]/rops_src[i] )


for i in xrange(nimage):
    img = get_im(img_stack, i)
    img_out = filt_table(img, table)
    img_out.write_image(out_stack, i)
    print i, ' done'

