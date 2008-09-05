#!/usr/bin/env python

from EMAN2 import *
from sparx import *

from sys import argv
from sys import exit

from string import atoi

if( len(argv) != 5 ) :
    print "Usage: mask img1 inner_radius outer_radius img_out"
    exit(-1)

file_in = argv[1]
inner_radius = atoi( argv[2] )
outer_radius = atoi( argv[3] )
file_out = argv[4]

a = getImage(file_in)

nx = a.get_xsize()
ny = a.get_ysize()
nz = a.get_zsize()

outer_sphere = model_circle( outer_radius, nx, ny, nz )
inner_sphere = model_circle( inner_radius, nx, ny, nz )

shell = outer_sphere - inner_sphere


nimage = EMUtil.get_image_count(file_in)

for i in xrange(nimage) :
    a = EMData()
    a.read_image( file_in, i )
    [mean_a,sigma,imin,imax]=Util.infomask(a,shell,False)
    
    va = threshold_to_minval( a, mean_a )

    print 'mean_a:', mean_a

    if nimage==1:
        dropImage( va, file_out, 's' )
    else:
        va.write_image( file_out, i )

    print "processed %4d" % i



