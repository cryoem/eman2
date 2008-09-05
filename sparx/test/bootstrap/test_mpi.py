#!/usr/bin/env python

import sys
from mpi import *
from numpy import array

sys.argv = mpi_init(len(sys.argv),sys.argv)

myid = mpi_comm_rank(MPI_COMM_WORLD)
size = mpi_comm_size(MPI_COMM_WORLD)

print "Hello from ",myid
print "Numprocs is ",size

vec1d1 = array( ([11.0*myid, 22.0*myid, 33.0*myid, 44.0*myid]), "f")
vec1d2 = array( ([12.0*myid, 24.0*myid, 36.0*myid, 48.0*myid]), "f")

vec2d  = array( ([vec1d1,vec1d2]), "f" )

print "myid,vec2d:", myid, vec2d

sum2d  = mpi_reduce(vec2d, 8, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD)

if myid==0 : print "myid,sum  :", myid, sum2d


from EMAN2 import *
from sparx import *
img = model_blank( 10, 1, 1, float(myid+1) )
img = fft(img)

print 'before reduce, myid, img: ', myid, img.get_value_at(0, 0, 0), img.get_value_at(1, 0, 0)
reduce_EMData_to_root( img, myid, 0 )
print ' after reduce, myid, img: ', myid, img.get_value_at(0, 0, 0), img.get_value_at(1, 0, 0)


