#!/usr/bin/env python

from EMAN2 import *
from sparx import *
from sys import argv
from string import split, atof
from math import cos, sin, exp

angs = open( argv[1], "r" )
hist = model_blank( 200, 200 )

s2 = 0.1
line  = angs.readline()
iang = 0
while len(line) > 0 :

    items = split( line )

    phi = atof( items[0] )
    tht = atof( items[1] )
    
    if tht < 0: 
        tht = tht + 360.0

    if tht > 270.0 :
        tht = 360.0 - tht
    if tht > 180.0 :
        tht = tht - 180.0
        phi = phi + 180.0
    elif tht > 90.0 :
        tht = 180.0 - tht
        phi = phi + 180.0
 
    if phi > 360.0:
        phi -= 360.0
 
    assert tht >= 0.0 and tht <= 90.0

    x = tht * cos( phi*3.1415926/180.0 ) + 100.0
    y = tht * sin( phi*3.1415926/180.0 ) + 100.0

    assert x > 0 and y > 0

    ix = int(x)
    iy = int(y)

    dx1 = x - ix
    dy1 = y - iy
   
    dx2 = ix + 1 - x
    dy2 = iy + 1 - y

    v = hist.get_value_at( ix, iy )
    v += exp( -(dx1*dx1+dy1*dy1)/s2 )
    hist.set_value_at( ix, iy, 0, v )

    v = hist.get_value_at( ix+1, iy )
    v += exp( -(dx2*dx2+dy1*dy1)/s2 )
    hist.set_value_at( ix+1, iy, 0, v )

    v = hist.get_value_at( ix, iy+1 )
    v += exp( -(dx1*dx1+dy2*dy2)/s2 )
    hist.set_value_at( ix, iy+1, 0, v )

    v = hist.get_value_at( ix+1, iy+1 )
    v += exp( -(dx2*dx2+dy2*dy2)/s2 )
    hist.set_value_at( ix+1, iy+1, 0, v )


    line = angs.readline()


    print iang, ' done'

    iang += 1


dropImage( hist, argv[2] )


