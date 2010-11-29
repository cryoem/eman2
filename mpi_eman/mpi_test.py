from EMAN2 import *
from mpi_eman import *
from cPickle import *

proc,nproc=mpi_init()

a=test_image()

print "Running on %d/%d"%(proc,nproc)
if proc==0 :
	mpi_send(dumps(a,-1),1,0)
	b=loads(mpi_recv(1,0)[0])
	mpi_send(dumps(b,-1),2,0)
	c=loads(mpi_recv(2,0)[0])
	
	b.write_image("final.hdf",0)
	c.write_image("final.hdf",1)
else :
	b=loads(mpi_recv(0,0)[0])
	b.mult(-1.0)
	mpi_send(dumps(b,-1),0,0)

#print "running on CPU ",proc
#if proc==0 : mpi_bcast("testing")
#else : print mpi_bcast(8,0)

mpi_finalize()

