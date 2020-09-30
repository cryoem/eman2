/*
 * Copyright (c) 2019- Baylor College of Medicine
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
 * */

/** This file is a part of util.h,
 * To use this file's functions, you should include "util.h"
 * NEVER directly include this file */

#ifndef util__sphire_h__
#define util__sphire_h__

// utility function used in (GPU) ISAC
static vector<int> sp_assign_groups(std::string matrix_address, int nref, int nima);

static vector<float> pw_extract_sphire(vector<float>pw, int n, int iswi,float ps);
static vector<float> call_cl1_sphire(long int *k,long int *n, float *ps, long int *iswi, float *pw, float *q2, double *q, double *x, double *res, double *cu, double *s, long int *iu);
static vector<float> lsfit_sphire(long int *ks,long int *n, long int *klm2d, long int *iswi, float *q1,double *q, double *x, double *res, double *cu, double *s,long int *iu);
static void cl1_sphire(long int *k, long int *l, long int *m, long int *n, long int *klm2d,double *q, double *x, double *res, double *cu, long
int *iu, double *s);

static void set_freq_sphire(EMData* freqvol, EMData* temp, EMData* mask, float cutoff, float freq);

#endif	//util__sphire_h__
