CC     = mpicxx

CFLAGS = -fPIC -ftemplate-depth-120 -w -DBOOST_DISABLE_THREADS  -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGE_FILES -DHAVE_FSEEKO -DHAVE_FTELLO -DFFTW2 -DEM_TIFF -DEM_PNG -DEM_BOOST -DUSE_MPI -DEM2_EXPORTS -Wall 

LIB    = -L../../../../install/lib -lEM2 -lEMext -L../../../../../usr/lib \
         -lgsl -lgslcblas -lhdf5 -lm

EMAN2INCLUDE = -I. -I../ -I../../ -I../../../ -I../../../../../usr/include 

CCTBXINCLUDE = -I../../../../../cctbx/cctbx_sources/scitbx/include \
               -I../../../../../cctbx/install/include 

.cpp.o:
	$(CC) $(CFLAGS) $(EMAN2INCLUDE) $(CCTBXINCLUDE) -c $<

OBJ = runrefine.o ali3d_d_mpi.o ali3d_unified_mpi.o alignoptions.o utilcomm.o

runrefine: $(OBJ)
	$(CC) $(OBJ) $(LIB) -o runrefine

clean:
	rm -f *.o runrefine
