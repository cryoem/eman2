#!/usr/bin/env python

from sparx import *
from EMAN2 import *
from sys import argv
from sys import exit
from string import atoi

def reverse_slice( slice ) :
    nx = slice.get_xsize()
    ny = slice.get_ysize()

    reverse = slice.copy()

    for iy in xrange( ny ) :
        for ix in xrange( nx ) :
	    value = slice.get_value_at( ix, iy, 0 )
	    reverse.set_value_at( ix, ny-1-iy, 0, value)
    
    return reverse


if( len(argv) != 4 ) :
    print "zslice vol_input iz slice_output"
    exit(-1)

a = getImage( argv[1] )

slice = get_z_slice( a, atoi( argv[2] ) )

nx = slice.get_xsize()
ny = slice.get_ysize()

mask = model_circle( nx/2 -3, nx, ny )

slice = slice*mask

reverse = reverse_slice( slice ) 

dropImage( reverse, argv[3] )

