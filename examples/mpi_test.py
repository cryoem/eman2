#!/usr/bin/env python

from sys import argv
from EMAN2 import *
from mpi import *
from mpi_eman import *

mpi_init(0,[])
proc=mpi_comm_rank(MPI_COMM_WORLD)
nproc=mpi_comm_size(MPI_COMM_WORLD)

a=test_image_3d(type=7,size=(256,256,256))

print "Running on %d/%d"%(proc,nproc)

# stage 1, synchronous send/recv
if proc==0 :
	a.write_image("final.hdf",0)

	for i in range(1,nproc):
		mpi_eman2_send(a,i,0)
		a=mpi_eman2_recv(i,0)[0]
		a.write_image("test_mpi_1.hdf",i)

else :
	b=mpi_eman2_recv(0,0)[0]
	b.mult(-1.5)
	mpi_eman2_send(b,0,0)

# stage 2, broadcast EMData
if proc==0:
#	a=test_image(1)
	a=test_image_3d(type=7,size=(195,195,195))
	mpi_bcast_send(a)
	
	allsrc=set(range(1,nproc))	# we use this to make sure we get a reply from all nodes
	while (1):
		if len(allsrc)==0 : break
		mpi_probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD)
		starr = mpi_status()
		src = (int)(starr[0])
		tag = (int)(starr[1])
		l = mpi_get_count(MPI_CHAR)

		print "%d replied"%src
		
		b=mpi_eman2_recv(src,tag)[0]
		b.write_image("test_mpi_2.hdf",src)
		allsrc.remove(src)

else:
	a=mpi_bcast_recv(0)
	a.mult(proc)
	mpi_eman2_send(a,0,1)

mpi_barrier(MPI_COMM_WORLD)

# stage 3, broadcast a string
if proc==0:
	a="You should see nodes-1 of these lines"
	mpi_bcast_send(a)
	
	allsrc=set(range(1,nproc))	# we use this to make sure we get a reply from all nodes
	while (1):
		if len(allsrc)==0 : break
		mpi_probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD)
		starr = mpi_status()
		src = (int)(starr[0])
		tag = (int)(starr[1])
		l = mpi_get_count(MPI_CHAR)
		print "%d replied"%src
		
		b=mpi_eman2_recv(src,tag)[0]
		print b
		allsrc.remove(src)

else:
	x=mpi_bcast_recv(0)
	mpi_eman2_send(x,0,1)

mpi_finalize()

if proc==0:
	print "done"

	print "\nIf you didn't see any errors above, then the test was a success."

#print "running on CPU ",proc
#if proc==0 : mpi_bcast("testing")
#else : print mpi_bcast(8,0)


