#!/usr/bin/env python

from sys import argv,stdout
from EMAN2 import *
from mpi import *
from mpi_eman import *

mpi_init(0,[])
mpi_barrier(MPI_COMM_WORLD)
proc=mpi_comm_rank(MPI_COMM_WORLD)
nproc=mpi_comm_size(MPI_COMM_WORLD)

a=test_image_3d(type=7,size=(64,64,64))

print "Running on %d/%d"%(proc,nproc)

# stage 1, synchronous send/recv
if proc==0 :
#	a.write_image("final.hdf",0)

	print "Stage 1, synchronous send/receive"
	print "Rank ",
	for i in range(1,nproc):
		mpi_eman2_send("DATA",a,i)
		com,data,src=mpi_eman2_recv(i)
#		a.write_image("test_mpi_1.hdf",i)
		print i,
		stdout.flush()
	print "\nStage 1 complete, all responses in"

else :
	com,data,src=mpi_eman2_recv(0)
	data.mult(-1.5)
	mpi_eman2_send("DATA",data,0)

# stage 2, broadcast EMData
if proc==0:
#	a=test_image(1)
	a=test_image_3d(type=7,size=(256,256,256))
	print "Stage 2, broadcast test"
	print "Rank ",
	mpi_bcast_send(a)
	
	allsrc=set(range(1,nproc))	# we use this to make sure we get a reply from all nodes
	while (1):
		if len(allsrc)==0 : break
		mpi_probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD)
		
		com,data,src=mpi_eman2_recv(MPI_ANY_SOURCE)
		print src,
		stdout.flush()

		allsrc.remove(src)

else:
	a=mpi_bcast_recv(0)
	a.mult(proc)
	mpi_eman2_send("DATA",a,0)

mpi_barrier(MPI_COMM_WORLD)

mpi_finalize()

if proc==0:
	print "\nStage 2, broadcast test complete"
	print "done"

	print "\nIf you didn't see any errors above, then the test was a success."

#print "running on CPU ",proc
#if proc==0 : mpi_bcast("testing")
#else : print mpi_bcast(8,0)


