#!/usr/bin/env python

from EMAN2 import *
from sparx import *

from sys import argv
from sys import exit
from string import replace

if len(argv) != 3 :
    print "Usage: split input output_prefix"
    exit(-1)

file_i = argv[1]
prefix = argv[2]

nimage = EMUtil.get_image_count( file_i )

for i in xrange(nimage) :
    data = EMData()
    data.read_image( file_i, i )
    file_o = prefix + ( "%4d.spi" % (i+1) )
    file_o = replace( file_o, ' ', '0' )
    data.write_image( file_o, 0, EMUtil.ImageType.IMAGE_SINGLE_SPIDER )
    print "file ", file_o, " written!"


