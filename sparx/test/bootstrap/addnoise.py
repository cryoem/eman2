#!/usr/bin/env python

from EMAN2 import *
from sparx import *
from sys import argv
from sys import exit
from string import atof

if len(argv) != 4 :
    print "Usage add_noise input_proj sigma output_proj"
    exit(-1)

nproj = EMUtil.get_image_count( argv[1] )
sigma = atof( argv[2] )

for i in xrange(nproj) :
    print "Proj # %4d" % i
    a = EMData()
    a.read_image( argv[1], i )

    nx = a.get_xsize()
    ny = a.get_ysize()
    nz = a.get_zsize()
    b = model_gauss_noise(sigma,nx,ny,nz)

    a = a + b
    a.write_image( argv[3], i )

