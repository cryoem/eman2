#!/usr/bin/env python

from sparx import *
from EMAN2 import *
from sys import argv
from sys import exit
from string import atoi

if(len(argv) != 3 and len(argv) != 4) :
    print "average input output [nimg]"
    exit(-1)

input = argv[1]
output = argv[2]

sum = None

if len(argv) == 4 :
    nimage = atoi( argv[3] )
else :
    nimage = EMUtil.get_image_count( input )

a = EMData()
for i in xrange( nimage ) :
    print "processing %4d" % i
    a.read_image( input, i )
    if( sum == None ) :
        sum = a.copy()
    else :
        sum = sum + a
    info(a)
    info(sum)

ave = sum / nimage

ave.write_image( output )


