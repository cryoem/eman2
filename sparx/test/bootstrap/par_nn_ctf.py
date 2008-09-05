#!/usr/bin/env python

from EMAN2 import *
from sparx import *
from mpi   import *
from string import *
from time  import time
from Numeric import reshape

import sys



sys.argv = mpi_init(len(sys.argv),sys.argv)
nproc = mpi_comm_size(MPI_COMM_WORLD)
myid = mpi_comm_rank(MPI_COMM_WORLD)



infofile = "progress%4d.txt" % (myid+1)
infofile = replace(infofile, ' ', '0')
info = open( infofile, 'w' )

if len(sys.argv) != 3 and len(sys.argv) != 6:
    print "usage: par_nnctf.py input_stack output_stack [start end step]"

prj_stack = sys.argv[1]
vol_stack = sys.argv[2]

time_start = time()

if len(sys.argv) == 3 :
    nimage = EMUtil.get_image_count(prj_stack)
    list_prj = range(0,nimage)
else :
    start = atoi( sys.argv[3] )
    end   = atoi( sys.argv[4] )
    step  = atoi( sys.argv[5] )
    list_prj = range(start, end, step)

     
nimage = len(list_prj)
nimage_per_node = nimage/nproc
image_start = myid * nimage_per_node

if myid == nproc-1 :
    image_end = nimage
else:
    image_end = image_start + nimage_per_node

npad = 4
symmetry = "c1"
snr = 10.0
sign = 1
Ttype = Transform3D.EulerType.SPIDER

prjlist = []
for i in range(image_start,image_end):
    prj = EMData()
    prj.read_image( prj_stack,list_prj[i] )
    prjlist.append( prj )
    info.write( "%4d read\n" % i )

vol = recons3d_4nn_ctf_MPI(myid, prjlist, snr, sign, symmetry, info)

if myid == 0 :
    dropImage( vol, vol_stack )
    info.write( "result wrote to " + vol_stack + "\n")
    info.write( "Total time: %10.3f\n" % (time()-time_start) )
    info.flush()

