#!/usr/bin/env python

from EMAN2 import *
from sparx import *

from sys import argv
from sys import exit

from string import atoi

if len(argv) != 2 and len(argv)!=3 :
    print "Usage: info img [radius]"
    exit(-1)

imgstack = argv[1]

a = getImage( imgstack )

xsize = a.get_xsize()
ysize = a.get_ysize()
zsize = a.get_zsize()

if len(argv)==2:
    mask = None
else:
    try:
        radius = atoi( argv[2] )
	mask = model_circle( radius, xsize, ysize, zsize )
    except:
    	mask = getImage( argv[2] )

nimage = EMUtil.get_image_count( imgstack )

for i in xrange(nimage):
    a = EMData()
    a.read_image( imgstack, i)
    [mean,sigma,min,max] = Util.infomask( a, mask, True )
    print "%8d %15.8f %15.8f %15.8f %15.8f" % ( i, mean, sigma, min, max )

