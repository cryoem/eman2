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

ncpu = 60
nimage = EMUtil.get_image_count( stack_in )
image_per_node = nimage/ncpu
imgid = 0

allscale = []
for inode in xrange(ncpu):
        image_start = inode*image_per_node
        image_end = image_start + image_per_node
        if inode==(ncpu-1):
            image_end = nimage
       
        filename = "%s%4d.txt" % (prefix, inode+1)
        filename = replace( filename, ' ', '0' )
        print 'loading ', filename
        param = readSpiderDoc( filename )

        for i in xrange( len(param) ):
            newscale = param[i][0]
            allscale.append( newscale )   
            imgid += 1
            print 'img %6d read' % imgid

        while imgid < image_end:
            allscale.append( 1.0 )
            imgid += 1
            print 'img %6d unchanged' % imgid

        assert imgid==image_end

avg = sum(allscale)/len(allscale)

for i in xrange(nimage):
    img = get_im( stack_in, i )
    img *= allscale[i]/avg
    img.write_image( stack_ot, i )
    print i, ' wrote'


