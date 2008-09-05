#!/usr/bin/env python

from EMAN2 import *
from sparx import *
from sys import argv
from sys import exit
from string import atof
from string import atoi
from random import normalvariate

if len(argv) != 6 :
    print "Usage gennoise model sigma low high output_stack"
    exit(-1)

nproj = EMUtil.get_image_count( argv[1] )
sigma = atof( argv[2] )
low   = atof( argv[3] )
high  = atof( argv[4] )

a = getImage( argv[1] )
size = a.get_xsize()

step = (high-low)/(size-1)

table =[]
for i in xrange(size):
    table.append( low + step*i )

print table

for i in xrange(nproj) :
    print "Proj # %4d" % i

    a =  EMData()
    a.read_image( argv[1], i )
    nx = a.get_xsize()
    ny = a.get_ysize()
    nz = a.get_zsize()

    b = model_gauss_noise(sigma,nx,ny,nz)
    
    c = filt_table( b, table )

    a += c

    a.write_image( argv[5], i )
