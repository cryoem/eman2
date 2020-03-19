#!/usr/bin/env python
from builtins import range
from sys import argv,stdout
from EMAN2 import *
from mpi import *
from mpi_eman import *

mpi_init(0,[])
mpi_barrier(MPI_COMM_WORLD)
proc=mpi_comm_rank(MPI_COMM_WORLD)
nproc=mpi_comm_size(MPI_COMM_WORLD)

a=test_image_3d(type=7,size=(64,64,64))

print("Running on %d/%d"%(proc,nproc))

# stage 1, synchronous send/recv
if proc==0 :
#	a.write_image("final.hdf",0)

	print("Stage 1, synchronous send/receive")
	print("Rank ", end=' ')
	for i in range(1,nproc):
		mpi_eman2_send("DATA",a,i)
		com,data,src=mpi_eman2_recv(i)
#		a.write_image("test_mpi_1.hdf",i)
		print(i, end=' ')
		stdout.flush()
	print("\nStage 1 complete, all responses in")

else :
	com,data,src=mpi_eman2_recv(0)
	data.mult(-1.5)
	mpi_eman2_send("DATA",data,0)

# stage 2, broadcast EMData
if proc==0:
#	a=test_image(1)
	a=test_image_3d(type=7,size=(256,256,256))
	print("Stage 2, broadcast test")
	print("Rank ", end=' ')
	mpi_bcast_send(a)
	
	allsrc=set(range(1,nproc))	# we use this to make sure we get a reply from all nodes
	while (1):
		if len(allsrc)==0 : break
		mpi_probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD)
		
		com,data,src=mpi_eman2_recv(MPI_ANY_SOURCE)
		print(src, end=' ')
		stdout.flush()

		allsrc.remove(src)
	print("\nStage 2, broadcast test complete")

else:
	a=mpi_bcast_recv(0)
	a.mult(proc)
	mpi_eman2_send("DATA",a,0)

mpi_barrier(MPI_COMM_WORLD)

mpi_finalize()
from pickle import dump,load
import select

fname='mpisock'
if proc==0:
	print("\nStage 3, test socket channel...")
	mpisock=socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
	
	try:  os.remove(fname)
	except:  pass
	mpisock.bind(fname)
	mpisock.listen(1)
	mpiconn, mpiaddr = mpisock.accept()
	print('r0: testtest')
	mpifile=mpiconn.makefile(mode='wb')
	mpifile.write(b"HELO")
	mpifile.flush()
	#s=mpifile.read(4)
	#print('r0:',s)
	dump(('DATA','testtest'),mpifile,-1)
	mpifile.flush()

elif proc==1:
	mpisock=socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
	time.sleep(2)
	mpisock.connect(fname)
	mpifile=mpisock.makefile(mode='rb')
	print("r1: Connected to Controller")

	# Initial handshake to make sure we're both here
	s=mpifile.read(4)
	print('r1:',s)

	if select.select([mpifile],[],[],0)[0]:
		com,data=load(mpifile)
		print('r1:',com,data)
	#print('done')


if proc==0:
	time.sleep(5)
	
	try:  os.remove(fname)
	except:  pass
	print("done")

	print("\nIf you didn't see any errors above, then the test was a success.")
	

#print "running on CPU ",proc
#if proc==0 : mpi_bcast("testing")
#else : print mpi_bcast(8,0)


