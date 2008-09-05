#!/usr/bin/env python

from EMAN2 import *
from sparx import *
from sys   import argv
from sys   import exit
from sys   import stdout
from string import atoi
from string import replace

#info = open( "pack.inf", "w" )

info = stdout

if( len( argv ) < 4 ):
    print "Usage: pack.py filelist output"
    exit(-1)



output = argv[-1]

total = 0
for input in argv[1:-1] :
    info.write("packing file " + input + "\n" )
    nimage = EMUtil.get_image_count( input )
    for j in xrange(nimage) :
        if( j%60==0 ) :
            info.write( "%4d " %j )
            info.flush()

        a = EMData()
        a.read_image( input, j )
        a.write_image(output, total+j)
        info.write( "." )
        info.flush()
        if( (j+1)%60==0 ):
            info.write( "\n" )

    if( nimage %60 != 0 ) :
        info.write( "\n" )

    info.write( " %4d volumes packed\n" % nimage )
    total = total + nimage

info.write( "totally %4d volumes packed\n" % total )
