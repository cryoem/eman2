CC = mpicxx

CFLAGS = -fPIC -ftemplate-depth-120 -w -DBOOST_DISABLE_THREADS  -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGE_FILES -DHAVE_FSEEKO -DHAVE_FTELLO -DFFTW2 -DEM_TIFF -DEM_PNG -DEM_BOOST -DUSE_MPI -DEM2_EXPORTS -Wall -O3 

LIB    = -L../../../../install/lib -lEM2 -lEMext -L../../../../../usr/lib \
         -lgsl -lgslcblas -lhdf5 -lm

EMAN2INCLUDE = -I. -I../ -I../../ -I../../../ -I../../../../../usr/include 

CCTBXINCLUDE = -I../../../../../cctbx/cctbx_sources/scitbx/include \
               -I../../../../../cctbx/install/include 

.cpp.o:
	$(CC) $(CFLAGS) $(EMAN2INCLUDE) $(CCTBXINCLUDE) -c $<

OBJ = ali3d_d_mpi.o ali3d_unified_mpi.o alignoptions.o utilcomm.o sirt.o

all: runhybrid rununified runali3d_d

runhybrid: runhybrid.o $(OBJ)
	$(CC) runhybrid.o $(OBJ) $(LIB) -o runhybrid

rununified: rununified.o ali3d_unified_mpi.o utilcomm.o
	$(CC) rununified.o ali3d_unified_mpi.o utilcomm.o $(LIB) -o rununified

runali3d_d: runali3d_d.o $(OBJ)
	$(CC) runali3d_d.o $(OBJ) $(LIB) -o runali3d_d

clean:
	rm -f *.o runhybrid rununified runali3d_d
