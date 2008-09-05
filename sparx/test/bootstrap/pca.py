#!/usr/bin/env python

from EMAN2 import *
from sparx import *

from sys import argv
from sys import exit

from string import atoi

def pca( input_stack, imglist, mask, nvec, type="out_of_core" ) :
    if type == "out_of_core" :
        ana = Analyzers.get( "pca_large", {"mask":mask, "nvec":nvec} )
    else :
        ana = Analyzers.get( "pca", {"mask":mask, "nvec":nvec} )

    for i in imglist:
        data = EMData()
        data.read_image( input_stack, i)
        ana.insert_image( data )
        print "Inserting image %4d" % i

    return ana.analyze()
       

if len(argv) != 5 and len(argv) !=6 :
    print "Usage: pca input_stack mask_radius nvec output_stack [nimg]"
    exit(-1)
  
input_stack = argv[1]
mask_radius = atoi( argv[2] )
nvec = atoi( argv[3] )
output_stack = argv[4]

if len(argv) == 6 :
    nimage = atoi(argv[5] )
else :
    nimage = EMUtil.get_image_count(input_stack); 

print "nimage:", nimage

data = getImage( argv[1] )
mask = model_circle( mask_radius, data.get_xsize(), data.get_ysize(), data.get_zsize() )

vecs = pca( input_stack, range(nimage), mask, nvec )

iout = 0
for vec in vecs:
    vec.write_image( output_stack, iout)
    iout = iout + 1

