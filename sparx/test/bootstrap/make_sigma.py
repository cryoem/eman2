#!/usr/bin/env python

from EMAN2 import *
from sparx import *
from sys import argv, exit

if len(argv)!=3:
    print "Usage: makevar.py prj_stack var_stack"
    exit(-1)


prj_stack = argv[1]
var_stack = argv[2]

nimage = EMUtil.get_image_count( prj_stack )

var = model_blank( nimage )

for i in xrange(nimage):
    prj = get_im( prj_stack, i )
    [favg, fsigma, fmin, fmax] = Util.infomask( prj, None, True )

    var.set_value_at( i, 0, 0, fsigma )

    print i, ' done'

dropImage( var, var_stack )
