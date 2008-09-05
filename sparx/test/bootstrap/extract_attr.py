#!/usr/bin/env python

from EMAN2 import *
from sparx import *

from sys import argv
from sys import exit
from sys import stdout

if len(argv) < 3:
    print "Usage: extract_attr img_stack attr1 attr2 ... "
    exit(-1)

nimage = EMUtil.get_image_count(argv[1])

for i in xrange(nimage):
    data = get_im( argv[1], i )
    for j in xrange(2, len(argv) ):
        a = data.get_attr( argv[j] )
        stdout.write( "%10.3e " % a )

    stdout.write("\n")
