MPI_INC=/usr/include
MPI_LIB=/usr/lib/libmpi.so
PYTHON=/usr/include/python2.4
NUMERIC=/Users/zhengfanyang/sparxroot/include/python2.5/
SITE=/usr/lib/python2.5/site-packages/mpi
SITE2=/users/tkaiser/bin/python/lib/python2.4/site-packages/mpi

PIC=-fPIC
default: mpi.so


mpi.so: mympimodule.c
	mkdir -p build/temp.linux-ia64-2.2
	mkdir -p build/lib.linux-ia64-2.2
	gcc -DNDEBUG -O2 -g -D_GNU_SOURCE -DNUMERIC $(PIC) -I$(MPI_INC) -I$(PYTHON)   -I$(NUMERIC)  -c mympimodule.c -o mympimodule.o
	gcc -shared -o mpi.so  mympimodule.o /usr/lib64/libpython2.4.so /usr/lib/libmpi.so


install:mpi.so
	mkdir -p $(SITE)
	cp mpi.so $(SITE)

clean:
	rm -rf mpi.so
	rm -rf mpi_tests/mpi.so
	rm -rf build
#	rm -rf $(SITE)

tar: clean
	tar -cf mympimod.tar makefile* mympimodule.c setupmpimod*
