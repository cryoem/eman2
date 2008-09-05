#!/usr/bin/env python

from EMAN2 import *
from sparx import *
from sys import argv, exit
from random import random

if len(argv)!=3:
	print "change_norm.py stack_in stack_ot"
	exit(-1)

stack_in = argv[1]
stack_ot = argv[2]


nimage = EMUtil.get_image_count( argv[1] )

for i in xrange( nimage ):
	img = get_im( stack_in, i )
	r = random() - 0.5
	r = r*0.3 + 1.0
	img *= r
	img.write_image( stack_ot, i )
	print 'i,r: %6d %10.8f' % (i,r)


