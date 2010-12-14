mpi_eman_c.so: mpi_eman_c.o
mpi_eman_c.o: mpi_eman.c
	mpicc -fno-strict-aliasing -DNDEBUG -g -fwrapv -O2 -Wall -Wstrict-prototypes -fPIC -I/usr/include/python2.6 -c mpi_eman.c -o mpi_eman_c.o
	mpicc -pthread -shared -Wl,-O1 -Wl,-Bsymbolic-functions mpi_eman_c.o -o mpi_eman_c.so

install: mpi_eman_c.so mpi_eman.py
	cp mpi_eman_c.so mpi_eman.py ${EMAN2DIR}/lib

clean:
	rm mpi_eman_c.o mpi_eman_c.so
