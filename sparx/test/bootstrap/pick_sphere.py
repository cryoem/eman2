#!/usr/bin/env python

from sparx import *
from EMAN2 import *
from sys import argv
from sys import exit
from string import atoi
from string import split

def calc_avg(img, sx, sy, sz, r) :
    nx = img.get_xsize()
    ny = img.get_ysize()
    nz = img.get_zsize()
    n  = 0 
    sum = 0.0
    for ix in xrange(nx):
        for iy in xrange(ny):
            for iz in xrange(nz):
	        dx = ix - sx
    	        dy = iy - sy
	        dz = iz - sz
	        if( dx*dx + dy*dy + dz*dz <= r*r ) :
                    n = n+1
                    sum = sum+img.get_value_at(ix,iy,iz)

    return sum/n

if( len(argv) != 6 and len(argv) != 4 ) :
    print "Usage: pick_sphere imgfile x y z r | pick_sphere imgfile -i points.in"
    exit(-1)

input = argv[1]

xl = []
yl = []
zl = []
rl = []

if( argv[2] == "-i" ) :
    fpoints = open( argv[3], "r" )
    line = fpoints.readline()
    while line:
        items = split(line)
	xl.append( atoi(items[0]) )
	yl.append( atoi(items[1]) )
	zl.append( atoi(items[2]) )
	rl.append( atoi(items[3]) )
	line = fpoints.readline()
else :
    xl.append( atoi(argv[2]) )
    yl.append( atoi(argv[3]) )
    zl.append( atoi(argv[4]) )
    rl.append( atoi(argv[5]) )


nimg = EMUtil.get_image_count( input ) 

for i in xrange(nimg):
    img = EMData()
    img.read_image(input, i)
    for j in xrange( len(xl) ):
        sx = xl[j]
	sy = yl[j]
	sz = zl[j]
	r  = rl[j]
        avg = calc_avg( img, sx, sy, sz, r )
        print "x,y,z,r,avg: %4d %4d %4d %4d %8.5f" % (sx,sy,sz,r,avg)
    print "\n"
