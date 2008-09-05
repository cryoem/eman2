#!/usr/bin/env python

from EMAN2 import *
from sparx import *

from sys import argv
from sys import exit

from string import atoi

if len(argv) !=3 and len(argv) != 4:
    print "Usage ccc img1 img2 [radius]"
    exit(-1)

nimage1 = EMUtil.get_image_count( argv[1] )
nimage2 = EMUtil.get_image_count( argv[2] )

img0 = getImage( argv[1] )

xsize = img0.get_xsize()
ysize = img0.get_ysize()
zsize = img0.get_zsize()


if len(argv) == 4 :
    try:
        radius = atoi( argv[3] )
        mask = model_circle( radius, xsize, ysize, zsize )
    except:
        mask = getImage( argv[3] )
else:
    mask = None

if( nimage2>nimage1 or nimage2==1 ) :
    nimage = nimage1
else:
    nimage = nimage2


for i in xrange(nimage) :
    img1 = EMData()
    img1.read_image( argv[1], i )
    img2 = EMData()

    if( nimage2 == 1 ) :
        img2.read_image( argv[2], 0 )
    else:
        img2.read_image( argv[2], i )

    if mask is None :
        result = ccc(img1,img2)
    else :
        result = ccc(img1, img2, mask )

    print "Image %5d %10.7f" % (i,result)
