#!/usr/bin/env python

from EMAN2  import *
from sparx  import *
from random import randint
from string import atoi
from string import atof
from mpi import *

import sys
import string

sys.argv = mpi_init( len(sys.argv), sys.argv )
size = mpi_comm_size(MPI_COMM_WORLD)
myid = mpi_comm_rank(MPI_COMM_WORLD)


#if len(argv) != 4 :
#    print "usage: bootstrap_genbuf proj_file npad bufprefix"
#    exit(-1)

proj_stack = sys.argv[1]
npad = atoi( sys.argv[2] )
buff_prefix = sys.argv[3]

store = file_store(buff_prefix, npad, 1)

mystatus = "genbuf%4d.txt" % ( myid )
mystatus = string.replace( mystatus, ' ', '0' )

output = open( mystatus, "w" )

nimage = EMUtil.get_image_count( proj_stack )
for i in xrange(nimage):
    proj = EMData()
    proj.read_image( proj_stack, i )
    store.add_image( proj )

    if( (i+1) % 100 == 0 ) :
        output.write( "proj %4d done\n" % (i+1) )
        output.flush()
   
output.write( "proj %4d done\n" % nimage )
output.flush
   
