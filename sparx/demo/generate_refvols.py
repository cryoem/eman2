#!/usr/bin/env python

from EMAN2 import *
from sparx import *
from sys import argv
from string import atof
from math import sqrt

vol = get_im( argv[1] )
var = get_im( argv[2] )
thr = atof( argv[3] )

var = threshold( var - thr )

nx = var.get_xsize()
ny = var.get_ysize()
nz = var.get_zsize()

stderr = model_blank(nx, ny, nz)
for ix in xrange(nx):
    for iy in xrange(ny):
        for iz in xrange(nz):
            v = var.get_value_at(ix, iy, iz)*1440
            stderr.set_value_at(ix, iy, iz, sqrt(v) )

v0 = vol - 3.0*stderr
v1 = vol + 3.0*stderr

v0.write_image( "rem.hdf", 0 )
v1.write_image( "rem.hdf", 1 )

