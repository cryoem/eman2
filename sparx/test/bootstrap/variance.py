#!/usr/bin/env python

from sparx import *
from EMAN2 import *
from sys import argv
from sys import exit, stdout
from string import atoi

if( len(argv) != 3 and len(argv) != 4 ) :
    print "variance input output [nimage]"
    exit(-1)

input = argv[1]
output = argv[2]

if len(argv) == 4 :
    nimage = atoi( argv[3] )
else :
    nimage = EMUtil.get_image_count( input )

var = var_bydef( input, range(nimage), stdout )

dropImage( var, output )
