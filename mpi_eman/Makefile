mpi_eman.so: mpi_eman.o
mpi_eman.o: mpi_eman.c
	mpicc -fno-strict-aliasing -DNDEBUG -g -fwrapv -O2 -Wall -Wstrict-prototypes -fPIC -I/usr/include/python2.6 -c mpi_eman.c -o mpi_eman.o
	mpicc -pthread -shared -Wl,-O1 -Wl,-Bsymbolic-functions build/temp.linux-x86_64-2.6/mpi_eman.o -o mpi_eman.so

install: mpi_eman.so
	cp mpi_eman.so ${EMAN2DIR}/lib

