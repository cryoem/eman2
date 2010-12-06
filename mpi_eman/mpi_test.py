from EMAN2 import *
from mpi_eman import *

proc,nproc=mpi_init()

a=test_image()

print "Running on %d/%d"%(proc,nproc)
if proc==0 :
	mpi_send(a,1,0)
	b=mpi_recv(1,0)[0]
	mpi_send(b,2,0)
	c=mpi_recv(2,0)[0]
	
	b.write_image("final.hdf",0)
	c.write_image("final.hdf",1)
else :
	b=mpi_recv(0,0)[0]
	b.mult(-1.0)
	mpi_send(b,0,0)

#print "running on CPU ",proc
#if proc==0 : mpi_bcast("testing")
#else : print mpi_bcast(8,0)

mpi_finalize()

