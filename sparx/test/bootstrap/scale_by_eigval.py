#!/usr/bin/env python

from EMAN2 import *
from sparx import *

from sys import argv
from sys import exit

if len(argv) != 3:
    print "scale_by_eigval input_stack output_stack"
    exit(-1)

nimage = EMUtil.get_image_count( argv[1] )

for i in xrange(nimage):
    from math import sqrt
    data = get_im( argv[1], i )
    eigval = data.get_attr( 'eigval' )
    data *= sqrt(eigval)
    data.write_image( argv[2], i )
    print i, ' done'
