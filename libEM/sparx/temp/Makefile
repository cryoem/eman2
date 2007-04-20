CC = mpicxx
CFLAGS = -fPIC -ftemplate-depth-120 -w -DBOOST_DISABLE_THREADS  -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGE_FILES -DHAVE_FSEEKO -DHAVE_FTELLO -DFFTW2 -DEM_TIFF -DEM_PNG -DEM_BOOST -DUSE_MPI -DEM2_EXPORTS -Wall 

#EMAN2INCLUDE = -I. -I/scratch/scratchdirs/cyang/Sparx/eman2/eman2 -I/scratch/scratchdirs/cyang/Sparx/eman2/eman2/libEM -I/scratch/scratchdirs/cyang/Sparx/eman2/eman2/libEM/sparx -I/scratch/scratchdirs/cyang/Sparx/usr/include 
#CCTBXINCLUDE = -I/scratch/scratchdirs/cyang/Sparx/cctbx/cctbx_sources/boost -Iinclude -I/scratch/scratchdirs/cyang/Sparx/cctbx/cctbx_sources/include -I/scratch/scratchdirs/cyang/Sparx/cctbx/cctbx_sources/scitbx/include/scitbx -I/scratch/scratchdirs/cyang/Sparx/cctbx/cctbx_sources/scitbx/include -I/scratch/scratchdirs/cyang/Sparx/cctbx/cctbx_sources/scitbx -I/scratch/scratchdirs/cyang/Sparx/cctbx/cctbx_build/include 

#CCTBXINCLUDE = -I/scratch/scratchdirs/cyang/Sparx/cctbx/cctbx_sources/scitbx/include -I/scratch/scratchdirs/cyang/Sparx/cctbx/install/include -I./
#EMAN2INCLUDE = -I. -I/home/psternberg/spx/eman2/install/include -I/home/psternberg/spx/eman2/eman2/libEM
EMAN2INCLUDE = -I. -I../../../install/include -I../..

#CCTBXINCLUDE = -I/home/psternberg/spx/cctbx/cctbx_sources/scitbx/include -I/home/psternberg/spx/cctbx/install/include
CCTBXINCLUDE = -I../../../../../cctbx/cctbx_sources/scitbx/include -I../../../../../cctbx/install/include

#LIB = -L/scratch/scratchdirs/cyang/Sparx/eman2/install/lib -lEM2 -lEMext -L/scratch/scratchdirs/cyang/Sparx/usr/lib -lgsl -lgslcblas -lhdf5 -lm
#LIB = -L/home/psternberg/spx/eman2/install/lib -L/home/psternberg/spx/usr/lib -lEM2 -lm
LIB = -L../../../../../eman2/install/lib -L../../../../../usr/lib -lEM2 -lm

.cpp.o:
	$(CC) $(CFLAGS) $(EMAN2INCLUDE) $(CCTBXINCLUDE) -c $<

OBJ = runrefine.o ali3d_d_mpi.o ali3d_unified_mpi.o alignoptions.o utilcomm.o

runrefine: $(OBJ)
	$(CC) $(OBJ) $(LIB) -o runrefine

clean:
	rm -f *.o runrefine
