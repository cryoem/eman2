#!/usr/bin/env python

from EMAN2 import *
from sparx import *

from sys import argv
from sys import exit
from sys import stdout

from string import atof

if len(argv) != 5 : 
    print "btwl input_stack low_range high_range output_stack"
    exit(-1)

input = argv[1]
lower = atof( argv[2] )
higher = atof( argv[3] )
output = argv[4]

nimage = EMUtil.get_image_count( input )

info = stdout

for i in xrange(nimage) :
    if( i % 60 == 0 ) :
        info.write( " %4d " % i )
        info.flush()

    data = EMData()
    data.read_image( input, i )
    btwl = filt_btwl( data, lower, higher )
    btwl.write_image( output, i )

    info.write( "." )
    info.flush()

    if( i % 60 == 59 ) :
        info.write( "\n" )

info.write( "\n" )
  
