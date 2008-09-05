#!/usr/bin/env python

from EMAN2 import *
from sparx import *

from sys import argv
from sys import exit

from string import atoi

if len(argv) !=4 and len(argv) != 5:
    print "Usage subtrace input1 input2 ouput [nimg]"
    exit(-1)


input1 = argv[1]
input2 = argv[2]
output = argv[3]

if len(argv) == 5 :
    nimage = atoi( argv[4] )
else :
    nimage = EMUtil.get_image_count( input1 )

nimage2 = EMUtil.get_image_count( input2 )

if( nimage2 != 1 and nimage2 != nimage ) :
    print "Error: input2 should single or input1's size"
    exit(-2)

for i in xrange(nimage): 
    a = EMData()
    a.read_image( input1, i )

    b = EMData()
    if nimage2 == 1 :
        b.read_image( input2, 0 )
    else:
        b.read_image( input2, i )

    a = a - b
    a.write_image( output, i )

    print "image %4d" % i


