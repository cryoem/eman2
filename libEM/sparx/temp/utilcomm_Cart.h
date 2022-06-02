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
#ifndef UTILCOMM_CART_H
#define UTILCOMM_CART_H

#include "mpi.h"
#include "stdlib.h"
#include "emdata.h"
#include "alignoptions.h"

using namespace EMAN;

int ReadStackandDist_Cart(MPI_Comm comm, EMData ***images2D, char *stackfname, int *nloc);
int CleanStack_Cart(MPI_Comm comm, EMData ** image_stack, int nloc, int ri, Vec3i volsize, Vec3i origin);
int ReadAngTrandDist_Cart(MPI_Comm comm_2d, MPI_Comm comm_row, int *dim, float *angleshift, char *angfname, int nloc);
int setpart_gc1(MPI_Comm comm_2d, int nangs, int *psize, int *nbase);
int setpart_gr1(MPI_Comm comm_2d, int nnz, int *nnzpart, int *nnzbase);

int sphpart(MPI_Comm comm_2d, int nrays, int *ptrs, int *nnzbase, int *ptrstart);
int getcb2sph(Vec3i volsize, int ri, Vec3i origin, int nnz0, int *ptrs, int *cord);

int fwdpj3_Cart(Vec3i volsize, int nraysloc, int nnzloc, float *dm, Vec3i origin, int ri, int *ptrs, int *cord, int myptrstart, float *x, float *y);

int bckpj3_Cart(Vec3i volsize, int nraysloc, int nnzloc, float *dm, Vec3i origin, int ri, int *ptrs, int *cord, int myptrstart, float *x, float *y);


#endif // UTILCOMM_CART_H

