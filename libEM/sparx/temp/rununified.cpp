/*
 * Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
 * Please do not copy or modify this file without written consent of the author.
 * Copyright (c) 2000-2019 The University of Texas - Houston Medical School
 *
 * This software is issued under a joint BSD/GNU license. You may use the
 * source code in this file under either license. However, note that the
 * complete EMAN2 and SPARX software packages have some GPL dependencies,
 * so you are responsible for compliance with the licenses of these packages
 * if you opt to use BSD licensing. The warranty disclaimer below holds
 * in either instance.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 */
#include "mpi.h"
#include "stdlib.h"

// include EMAN2
#include "emdata.h"

#include "ali3d_d_mpi.h"
#include "ali3d_unified_mpi.h"
#include "alignoptions.h"
#include "utilcomm.h"

int main(int argc, char *argv[])
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int ncpus, mypid, ierr, mpierr=0;
    int nloc; 
    double t0;

    MPI_Status mpistatus;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(comm,&ncpus);
    MPI_Comm_rank(comm,&mypid);
    printf("mypid = %d, ncpus = %d\n", mypid, ncpus);
    char  volfname[100], paramfname[100], stackfname[100],voutfname[100];
    EMData **expimages;
    int maxiter=0;

    // parse the command line and set filenames	
    if (argc < 4) {
      if (mypid == 0) {
          printf("Not enough arguments to the command...\n");
          printf("Usage: rununified -data=<imagestack> ");
          printf("-model=<initial 3D volume filename> ");
          printf("-param=<initial angles&shifts> ");
          printf("-out=<output volume filename>\n"); 
          printf("[-maxit=<maximum number of iterations>]\n"); 
          printf("[-sym=<symmtry type>]\n"); 
          printf("[-CTF]\n"); 
      }
      ierr = MPI_Finalize();
      exit(1);
    }

    int ia=1;
    while (ia < argc) {
       if ( !strncmp(argv[ia],"-data",5) ) {
          strcpy(stackfname,&argv[ia][6]);
          if (mypid == 0)
             printf("image stack        = %s\n", stackfname);
       }
       else if ( !strncmp(argv[ia],"-model",6) ) {
          strcpy(volfname,&argv[ia][7]);
          if (mypid == 0)
             printf("initial model      = %s\n", volfname);
       }
       else if ( !strncmp(argv[ia],"-param",6) ) {
          strcpy(paramfname,&argv[ia][7]);
          if (mypid == 0)
             printf("initial parameters = %s\n", paramfname);
       }
       else if ( !strncmp(argv[ia],"-maxit",6) ) {
	  maxiter = atoi(&(argv[ia][7]));
          if (mypid == 0)
             printf("maxiter = %d\n", maxiter);
       }
       else if ( !strncmp(argv[ia],"-out",4) ) {
          strcpy(voutfname,&argv[ia][5]);
          if (mypid == 0)
             printf("output volume file = %s\n", voutfname);
       }
       ia++;
    }

    // read and broadcast the initial model
    t0 = MPI_Wtime();
    EMData *volume = new EMData();
    ierr = ReadVandBcast(comm, volume, volfname);
    if (ierr == 0) {
	if (mypid == 0) {
	   printf("Finished reading and replicating volume\n");
	   printf("I/O time for reading volume = %11.3e\n",
		  MPI_Wtime() - t0);
	}
    }
    else {
       if (mypid == 0) printf("Failed to read the model volume %s! exit...\n", volfname);
       ierr = MPI_Finalize();
       exit(1);
    }
    
    // read and distribute a stack of experimental images
    t0 = MPI_Wtime();
    ierr = ReadStackandDist(comm, &expimages, stackfname,&nloc);
    if (ierr == 0) {
	if (mypid == 0) {
	   printf("Finished reading and distributing image stack\n");
	   printf("I/O time for reading image stack = %11.3e\n",
		  MPI_Wtime() - t0);
	}
    }
    else {
       if (mypid == 0) printf("Failed to read the image stack %s! exit...\n", stackfname);
       ierr = MPI_Finalize();
       exit(1);
    }

    Vec3i volsize;
    Vec3i origin;
    volsize[0] = volume->get_xsize();
    volsize[1] = volume->get_ysize();
    volsize[2] = volume->get_zsize();
    origin[0] = volume->get_xsize()/2 + 1;
    origin[1] = volume->get_ysize()/2 + 1;
    origin[2] = volume->get_zsize()/2 + 1;
    int    ri = volume->get_xsize()/2 - 2;
    ierr = CleanStack(comm, expimages, nloc, ri, volsize, origin);

    // read angles and shifts from the parameter file
    float * angleshift = new float[5*nloc];
    ierr = ReadAngTrandDist(comm, angleshift, paramfname, nloc);
    if (ierr!=0) { 
       mpierr = MPI_Finalize();
       return 1;
    }

    if (maxiter <= 0) maxiter = 10;  
    try {
	unified(comm, volume, expimages, angleshift, nloc, 
                maxiter, voutfname);
    }
    catch (std::exception const& e) {
	printf("%s\n", e.what());
    }

#ifdef DEBUG
    for (int i = 0; i < nloc; i++)
        printf("%11.3e   %11.3e   %11.3e   %11.3e   %11.3e\n", 
           angleshift[5*i+0], 
           angleshift[5*i+1],
           angleshift[5*i+2],
           angleshift[5*i+3],
           angleshift[5*i+4]);
#endif

    EMDeletePtr(volume);
    for ( int i = 0 ; i < nloc; ++i ) {
	EMDeletePtr(expimages[i]);
    }
    EMDeleteArray(expimages);
    EMDeleteArray(angleshift);
    ierr = MPI_Finalize();
    return 0;
}

