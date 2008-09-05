#!/usr/bin/env python

from EMAN2 import *
from sparx import *

from sys import exit
from sys import argv

from string import replace

if len(argv) != 4 :
    print "apply_ang.py stack_in prefix stack_ot"
    exit(-1)

stack_in = argv[1]
prefix   = argv[2]
stack_ot = argv[3]

already  = EMUtil.get_image_count( argv[3] )


nimage = EMUtil.get_image_count( stack_in )
image_per_node = nimage/64
imgid = 0
for inode in xrange(64):

        image_start = inode*image_per_node
        image_end = image_start + image_per_node
        if inode==63:
            image_end = nimage
 
        if( already > image_end ):
            continue

        imgid = image_start
      
        filename = "%s%4d.txt" % (prefix, inode)
        filename = replace( filename, ' ', '0' )
        print 'loading ', filename
        param = readSpiderDoc( filename )

        for i in xrange( len(param) ):
            if( imgid < already ):
                imgid += 1
                continue

            img = get_im( stack_in, imgid )

            newscale = param[i][0]
    
            img *= newscale

            img.write_image( stack_ot, imgid )

            imgid += 1

            print 'img %6d done' % imgid

        assert imgid==image_end

