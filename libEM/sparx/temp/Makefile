CC = mpicxx

CFLAGS = -fPIC -ftemplate-depth-120 -w -DBOOST_DISABLE_THREADS  -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGE_FILES -DHAVE_FSEEKO -DHAVE_FTELLO -DFFTW2 -DEM_TIFF -DEM_PNG -DEM_BOOST -DUSE_MPI -DEM2_EXPORTS -Wall -O3 

LIBS    = -L../../../../install/lib -lEM2 -L../../../../../usr/lib \
         -lgsl -lgslcblas -lhdf5 -lm

EMAN2INCLUDE = -I. -I../ -I../../ -I../../../ -I../../../../../usr/include 

CCTBXINCLUDE = -I../../../../../cctbx/cctbx_sources/scitbx/include \
               -I../../../../../cctbx/install/include 

.cpp.o:
	$(CC) $(CFLAGS) $(EMAN2INCLUDE) $(CCTBXINCLUDE) -c $<

OBJS =   alignoptions.o utilcomm.o project3d.o fgcalc.o utilnum.o utilparse.o

OBJ_CART = utilcomm_Cart.o project3d_Cart.o project3d.o utilnum.o

all: runhybrid rununified runali3d_d runsirt runsirt_cart

runhybrid: runhybrid.o ali3d_d_mpi.o ali3d_unified_mpi.o sirt.o $(OBJS)
	$(CC) runhybrid.o ali3d_d_mpi.o ali3d_unified_mpi.o sirt.o $(OBJS) $(LIBS) -o runhybrid

rununified: rununified.o ali3d_unified_mpi.o $(OBJS) 
	$(CC) rununified.o ali3d_unified_mpi.o $(OBJS) $(LIBS) -o rununified

runali3d_d: runali3d_d.o ali3d_d_mpi.o sirt.o $(OBJS)
	$(CC) runali3d_d.o ali3d_d_mpi.o sirt.o $(OBJS) $(LIBS) -o runali3d_d

runsirt: runsirt.o sirt.o $(OBJS) 
	$(CC) runsirt.o sirt.o $(OBJS) $(LIBS) -o runsirt

runsirt_cart: runsirt_Cart.o sirt_Cart.o HyBR_Cart.o $(OBJ_CART) 
	$(CC) runsirt_Cart.o sirt_Cart.o HyBR_Cart.o $(OBJ_CART) $(LIBS) -o runsirt_cart

clean:
	rm -f *.o runhybrid rununified runali3d_d runsirt runsirt_cart
