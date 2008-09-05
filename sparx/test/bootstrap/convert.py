#!/usr/bin/env python

from EMAN2 import *
from sparx import *

from sys import argv
from sys import exit

if( len(argv) != 3 ) :
    print "Usage: convert hdf_in spi_out"
    exit(-1)


a = getImage( argv[1] )
dropImage( a, argv[2] )
