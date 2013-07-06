#!/usr/bin/env python

from EMAN2 import *
from mpi_eman import *

proc,nproc=mpi_init()

a=test_image()

print "Running on %d/%d"%(proc,nproc)

# stage 1, synchronous send/recv
if proc==0 :
	a.write_image("final.hdf",0)

	for i in range(1,nproc):
		mpi_send(a,i,0)
		a=mpi_recv(i,0)[0]
		a.write_image("test_mpi_1.hdf",i)

else :
	b=mpi_recv(0,0)[0]
	b.mult(-1.5)
	mpi_send(b,0,0)

# stage 2, broadcast EMData
if proc==0:
	a=test_image(1)
	mpi_bcast_send(a)
	
	allsrc=set(range(1,nproc))	# we use this to make sure we get a reply from all nodes
	while (1):
		if len(allsrc)==0 : break
		l,src,tag = mpi_probe(-1,-1)
		print "%d replied"%src
		
		b=mpi_recv(src,tag)[0]
		b.write_image("test_mpi_2.hdf",src)
		allsrc.remove(src)

else:
	a=mpi_bcast_recv(0)
	a.mult(proc)
	mpi_send(a,0,1)

mpi_barrier()

# stage 3, broadcast a string
if proc==0:
	a="You should see nodes-1 of these lines"
	mpi_bcast_send(a)
	
	allsrc=set(range(1,nproc))	# we use this to make sure we get a reply from all nodes
	while (1):
		if len(allsrc)==0 : break
		l,src,tag = mpi_probe(-1,-1)
		print "%d replied"%src
		
		b=mpi_recv(src,tag)[0]
		print b
		allsrc.remove(src)

else:
	x=mpi_bcast_recv(0)
	mpi_send(x,0,1)

mpi_finalize()

if proc==0:
	print "done"

	print "\nIf you didn't see any errors above, then the test was a success."

#print "running on CPU ",proc
#if proc==0 : mpi_bcast("testing")
#else : print mpi_bcast(8,0)


