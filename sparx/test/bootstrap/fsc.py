#!/usr/bin/env python

from EMAN2 import *
from sparx import *

from sys import argv
from sys import exit

from string import atoi

if( len(argv) != 4  and len(argv) != 5 ) :
    print "Usage: fsc img1 img2 output [radius]"
    exit(-1)

a = getImage( argv[1] )
b = getImage( argv[2] )

nx = a.get_xsize()
ny = a.get_ysize()
nz = a.get_zsize()

if( len(argv) == 5 ) :
    radius = atoi( argv[4] )
    mask = model_circle( radius, nx, ny, nz )
    [mean_a,sigma,imin,imax]=Util.infomask(a,mask,False)
    [mean_b,sigma,imin,imax]=Util.infomask(b,mask,False)

    a = (a-mean_a)*mask
    b = (b-mean_b)*mask

    print "before mask, mean_a is ", mean_a 
    print "before mask, mean_b is ", mean_b

    [mean_a,sigma,imin,imax]=Util.infomask(a,mask,False)
    [mean_b,sigma,imin,imax]=Util.infomask(b,mask,False)

    if( mean_a != 0.0 ) :
        print "Waring: average of image1 under mask is ", mean_a, ", should be zero"

    if( mean_b != 0.0 ) :
        print "Waring: average of image2 under mask is ", mean_b, ", should be zero"

data = fsc( a, b, filename=argv[3] )


