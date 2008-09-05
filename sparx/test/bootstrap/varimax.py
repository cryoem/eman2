#!/usr/bin/env python

from EMAN2 import *
from sparx import *

from sys import argv
from sys import exit

from string import atoi

def varimax( input_stack, imglist, mask ) :
    ana = Analyzers.get( "varimax", {"mask":mask} )

    for i in imglist:
        data = EMData()
        data.read_image( input_stack, i)
        ana.insert_image( data )
        print "Inserting image %4d" % i

    return ana.analyze()
       

if len(argv) != 6:
    print "Usage: varimax input_stack mask_radius imgstart imgend output_stack"
    exit(-1)
  
input_stack = argv[1]
mask_radius = atoi( argv[2] )
imgstart = atoi( argv[3] )
imgend   = atoi( argv[4] )
output_stack = argv[5]

data = getImage( argv[1] )
mask = model_circle( mask_radius, data.get_xsize(), data.get_ysize(), data.get_zsize() )

totsize = data.get_xsize()*data.get_ysize()*data.get_zsize()

print "before rotation:"
for i in xrange(imgstart,imgend):
    data = EMData()
    data.read_image(input_stack, i)
    data *= mask
    print "i,dot: %4d %10.3f" % (i, -data.cmp("dot", data)*totsize )

vecs = varimax(input_stack, range(imgstart,imgend), mask)

print "after rotation:"
for i in xrange( len(vecs) ):
    vecs[i].write_image(output_stack, i)
    print "i,dot: %4d %10.3f" % (i, -vecs[i].cmp("dot", vecs[i])*totsize )

