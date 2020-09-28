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

#ifdef _WIN32
#pragma warning(disable:4819)
#include <malloc.h>
#endif	//_WIN32

#include <cstring>
#include <string>
#include <ctime>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <sstream>
#include "emdata.h"
#include "util.h"
using namespace EMAN;
#include "emassert.h"
#include "randnum.h"
#include <algorithm>
#include <vector>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_bessel.h>
#include <cmath>
#include "emfft.h"
#include <fftw3.h>
#include <limits.h>
#include <float.h>

using namespace std;
using std::complex;


#define    QUADPI      		    3.141592653589793238462643383279502884197
#define    PI2                  QUADPI/2.0
#define    TWOPI                2*QUADPI

#define deg_rad  QUADPI/180.0
#define rad_deg  180.0/QUADPI

///*
// * This file is based largely on the following software distribution:
// *
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// *
// *                              FFTPACK
// *
// * Reference
// *    P.N. Swarztrauber, Vectorizing the FFTs, in Parallel Computations
// *    (G. Rodrigue, ed.), Academic Press, 1982, pp. 51--83.
// *
// *     http://www.netlib.org/fftpack/
// *
// * Updated to single, double, and extended precision,
// * and translated to ISO-Standard C/C++ (without aliasing)
// * on 10 October 2005 by Andrew Fernandes <andrew_AT_fernandes.org>
// *
// *                   Version 4  April 1985
// *
// *      A Package of Fortran Subprograms for the Fast Fourier
// *       Transform of Periodic and other Symmetric Sequences
// *
// *                          by
// *
// *                   Paul N Swarztrauber
// *
// *   National Center for Atmospheric Research, Boulder, Colorado 80307,
// *
// *    which is sponsored by the National Science Foundation
// *
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// *
// * There appears to be no explicit license for FFTPACK. However, the
// * package has been incorporated verbatim into a large number of software
// * systems over the years with numerous types of license without complaint
// * from the original author; therefore it would appear
// * that the code is effectively public domain. If you are in doubt,
// * however, you will need to contact the author or the  National Center
// * for Atmospheric Research to be sure.
// *
// * All the changes from the original FFTPACK to the current file
// * fall under the following BSD-style open-source license:
// *
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
// *
// * Copyright (c) 2005, Andrew Fernandes (andrew@fernandes.org);
// * All rights reserved.
// *
// * Redistribution and use in source and binary forms, with or without
// * modification, are permitted provided that the following conditions
// * are met:
// *
// * - Redistributions of source code must retain the above copyright
// * notice, this list of conditions and the following disclaimer.
// *
// * - Redistributions in binary form must reproduce the above copyright
// * notice, this list of conditions and the following disclaimer in the
// * documentation and/or other materials provided with the distribution.
// *
// * - Neither the name of the North Carolina State University nor the
// * names of its contributors may be used to endorse or promote products
// * derived from this software without specific prior written permission.
// *
// * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
// * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
// * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// * POSSIBILITY OF SUCH DAMAGE.
// *
// */
//#include "fftpack4.h"

#ifdef __cplusplus

#include <cmath> /* the correct precision will be automatically selected */
using std::cos;
using std::sin;

#else /* ! __cplusplus */

#include <math.h> /* you must define/typedef the functions 'cos/cosf/cosl' and 'sin/sinf/sinl' as appropriate */
/* real_t cos(real_t); */
/* real_t sin(real_t); */
#endif

// utility struct for sp_assign_groups()
struct assign_groups_comparator {
	const float * values;
	bool operator() (int i,int j) { return (values[i] > values[j]); }
	assign_groups_comparator(const float *v) : values(v) {}
};

// utility function used in (GPU) ISAC
vector<int> Util::sp_assign_groups(std::string matrix_address, int nref, int nima)
{
	// convert memory address sent as string to float pointer
	const float * matrix;
	size_t addr = 0;
	for ( std::string::const_iterator i = matrix_address.begin();  i != matrix_address.end();  ++i ) {
		int digit = *i - '0';
		addr *= 10;
		addr += digit;
	}
	matrix = reinterpret_cast<float*>(addr);

	int kt = nref;
	unsigned int maxasi = (unsigned int)(nima/nref);
	vector< vector<int> > id_list(nref);
	int group, ima;

	// allocate and sort vector of indices
	size_t count = (size_t)nref * (size_t)nima; // danger: nref*nima will not fit into int
	std::vector<int> dd(count) ;
	for (size_t i=0; i<count; i++){
		dd[i] = i;
	}

	assign_groups_comparator comparator(matrix);
	sort(dd.begin(), dd.end(), comparator);
	
	// main loop
	size_t begin = 0;
	std::vector<bool> del_row(nref, false);
	std::vector<bool> del_column(nima, false);
	while (kt > 0) {
		bool flag = true;
		while (flag) {
			int l = dd[begin];
			group = l/nima;
			ima = l%nima;
			if (del_column[ima] || del_row[group])
				begin++;
			else
				flag = false;
		}

		id_list[group].push_back(ima);
		if (kt > 1) {
			if (id_list[group].size() < maxasi)
				group = -1;
			else
				kt -= 1;
		} else {
			if (id_list[group].size() < maxasi+nima%nref)
				group = -1;
			else
				kt -= 1;
		}
		del_column[ima] = true;
		if (group != -1) {
			del_row[group] = true;
		}
	}
	
	vector<int> id_list_1;

	for (int iref=0; iref<nref; iref++){
		for (unsigned int im=0; im<maxasi; im++){
			id_list_1.push_back(id_list[iref][im]);
		}
	}
	for (unsigned int im=maxasi; im<maxasi+nima%nref; im++){
		id_list_1.push_back(id_list[group][im]);
	}

	id_list_1.push_back(group);

	return id_list_1;
}

vector<float> Util::pw_extract_sphire(vector<float>pw, int n, int iswi, float ps)
{
	int k,m,n1,klmd,klm2d,nklmd,n2d,n_larg,l, n2;

	k=(int)pw.size();
	l=0;
	m=k;
	n2=n+2;
	n1=n+1;
	klmd=k+l+m;
	klm2d= k+l+m+2;
	nklmd=k+l+m+n;
	n2d=n+2;
	/*size has to be increased when N is large*/
	n_larg=klmd*2;
	klm2d=n_larg+klm2d;
	klmd=n_larg+klmd;
	nklmd=n_larg+nklmd;
	int size_q=klm2d*n2d;
	int size_cu=nklmd*2;
	static int i__;

	double *q ;
	double *x ;
	double *res;
	double *cu;
	float  *q2;
	float  *pw_;
	long int *iu;
	double *s;
	q   = (double*)calloc(size_q,sizeof(double));
	x   = (double*)calloc(n2d,sizeof(double));
	res = (double*)calloc(klmd,sizeof(double));
	cu  = (double*)calloc(size_cu,sizeof(double));
	s   = (double*)calloc(klmd,sizeof(double));
	q2  = (float*)calloc(size_q,sizeof(float));
	iu  = (long int*)calloc(size_cu,sizeof(long int));
	pw_ = (float*)calloc(k,sizeof(float));

	for( i__ =0; i__<k; ++i__) pw_[i__] = log(pw[i__]);
	long int l_k=k;
	long int l_n=n;
	long int l_iswi=iswi;
	vector<float> cl1_res;
	cl1_res = Util::call_cl1_sphire(&l_k, &l_n, &ps, &l_iswi, pw_, q2, q, x, res, cu, s, iu);
	free(q);
	free(x);
	free(res);
	free(s);
	free(cu);
	free(q2);
	free(iu);
	free(pw_);
	return cl1_res;
}



vector<float> Util::call_cl1_sphire(long int *k, long int *n, float *ps, long int *iswi, float *pw, float *q2,double *q, double *x, double *res, double *cu, double *s, long int *iu)
{
    long int q2_dim1, q2_offset, q_dim1, q_offset, i__1, i__2;
    float r__1;
    int tmp__i;
    long int i__, j;
    --s;
    --res;
    iu -= 3;
    cu -= 3;
    --x;
    long int klm2d;
    klm2d= *k+*k+2;
    klm2d=klm2d+klm2d;
    q_dim1 = klm2d;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    q2_dim1 = klm2d;
    q2_offset = 1 + q2_dim1;
    q2 -= q2_offset;
    i__2=0;
    i__1 = *n - 1;
    tmp__i=0;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *k;
	tmp__i+=1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    r__1 = float(i__ - 1) /(float) *k / (*ps * 2);
	    q2[i__ + j * q2_dim1] = pow(r__1, tmp__i);
	    }
    }
    for  (i__ = 1; i__ <= i__2; ++i__)
      { q2[i__ + *n * q2_dim1] = 1.f;
	    q2[i__ + (*n + 1) * q2_dim1] = pw[i__-1];
	}
   vector<float> fit_res;
   fit_res=Util::lsfit_sphire(k, n, &klm2d, iswi, &q2[q2_offset], &q[q_offset], &x[1], &res[1], &cu[3], &s[1], &iu[3]);
   return fit_res;
}
vector<float> Util::lsfit_sphire(long int *ks,long int *n, long int *klm2d, long int *iswi,float *q1,double *q, double *x, double *res, double *cu, double *s, long int *iu)
{
    /* System generated locals */
    long int q_dim1, q_offset, q1_dim1, q1_offset, i__1, i__2;

    /* Local variables */
    long int i__, j, m, n1, ii, jj;
    double tmp;
    vector<float> p;
    --x;
    q_dim1 = *klm2d;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    q1_dim1 = *klm2d;
    q1_offset = 1 + q1_dim1;
    q1 -= q1_offset;
    --s;
    --res;
    iu -= 3;
    cu -= 3;

    /* Function Body */
    long int l = 0;

/* C==ZHONG HUANG,JULY,12,02;L=0,1,2,3,4,5,6 correspond to different equality constraints */
    m = *ks;
    n1 = *n + 1;
    if (*iswi == 1) {
	i__1 = n1;
	for (jj = 1; jj <= i__1; ++jj) {
	    i__2 = *ks;
	    for (ii = 1; ii <= i__2; ++ii) {
	/*	q[ii + jj * q_dim1] = (double) q1[ii + jj * q1_dim1];*/

		q[*ks + ii + jj * q_dim1] = (double) q1[ii + jj * q1_dim1]
			;
	    }
	}
    } else if (*iswi == 2) {
	i__1 = *ks;
	for (ii = 1; ii <= i__1; ++ii) {
	    i__2 = n1;
	    for (jj = 1; jj <= i__2; ++jj) {
		q[ii + jj * q_dim1] = (double) q1[ii + jj * q1_dim1];
		q[*ks + ii + jj * q_dim1] = -((double) q1[ii + jj * q1_dim1]);
	    }
	}
    } else if (*iswi == 3) {
	l = 2;
	i__1 = n1;
	for (jj = 1; jj <= i__1; ++jj) {
	    i__2 = *ks + 2;
	    for (ii = 1; ii <= i__2; ++ii) {
		q[ii + jj * q_dim1] = (double) q1[ii + jj * q1_dim1];
	    }
	    i__2 = *ks;
	    for (ii = 1; ii <= i__2; ++ii) {
		q[*ks + 2 + ii + jj * q_dim1] = (double) q1[ii + jj * q1_dim1];
	    }
	}
    } else if (*iswi == 4) {
	l = 2;
	i__1 = n1;
	for (jj = 1; jj <= i__1; ++jj) {
	    i__2 = *ks + 2;
	    for (ii = 1; ii <= i__2; ++ii) {
		q[ii + jj * q_dim1] = (double) q1[ii + jj * q1_dim1];
	    }
	    i__2 = *ks;
	    for (ii = 1; ii <= i__2; ++ii) {
		q[*ks + 2 + ii + jj * q_dim1] = -((double) q1[ii + jj * q1_dim1]);
	    }
	}
    } else if (*iswi == 5) {
	l = 1;
	i__1 = n1;
	for (jj = 1; jj <= i__1; ++jj) {
	    i__2 = *ks + 1;
	    for (ii = 1; ii <= i__2; ++ii) {
		q[ii + jj * q_dim1] = (double) q1[ii + jj * q1_dim1];
	    }
	    i__2 = *ks;
	    for (ii = 1; ii <= i__2; ++ii) {
		q[*ks + 1 + ii + jj * q_dim1] = -((double) q1[ii + jj * q1_dim1]);
	    }
	}
    } else if (*iswi == 6) {
	l = 1;
	i__1 = n1;
	for (jj = 1; jj <= i__1; ++jj) {
	    i__2 = *ks + 1;
	    for (ii = 1; ii <= i__2; ++ii) {
		q[ii + jj * q_dim1] = (double) q1[ii + jj * q1_dim1];
	    }
	    i__2 = *ks;
	    for (ii = 1; ii <= i__2; ++ii) {
		q[*ks + 1 + ii + jj * q_dim1] = (double) q1[ii + jj * q1_dim1];
	    }
	}
    } else if (*iswi == 7) {
	l = 3;
	i__1 = n1;
	for (jj = 1; jj <= i__1; ++jj) {
	    i__2 = *ks + 3;
	    for (ii = 1; ii <= i__2; ++ii) {
		q[ii + jj * q_dim1] = (double) q1[ii + jj * q1_dim1];
	    }
	    i__2 = *ks;
	    for (ii = 1; ii <= i__2; ++ii) {
		q[*ks + 3 + ii + jj * q_dim1] = -((double) q1[ii + jj * q1_dim1]);
	    }
	}
    } else if (*iswi == 8) {
	l = 4;
	i__1 = n1;
	for (jj = 1; jj <= i__1; ++jj) {
	    i__2 = *ks + 4;
	    for (ii = 1; ii <= i__2; ++ii) {
		q[ii + jj * q_dim1] = (double) q1[ii + jj * q1_dim1];
	    }
	    i__2 = *ks;
	    for (ii = 1; ii <= i__2; ++ii) {
		q[*ks + 4 + ii + jj * q_dim1] = -((double) q1[ii + jj * q1_dim1]);
	    }
	}
    }

    Util::cl1_sphire(ks, &l, &m, n, klm2d, &q[q_offset], &x[1], &res[1], &cu[3], &iu[3], &s[1]);
    i__1 = *ks;
    int tmp__j=0;
    for (i__ = 1; i__ <= i__1; ++i__) {
	tmp = 0.f;
	i__2 = *n - 1;
	for (j = 1; j <= i__2; ++j) {
	    tmp__j=j;
	    tmp += pow(q1[i__ + q1_dim1], tmp__j) * x[j];
	    }
	tmp += x[*n];
	p.push_back(static_cast<float>(exp(tmp)));
	p.push_back(q1[i__ + q1_dim1]);
    }
    i__2=*n;
    for (i__=1;i__<=i__2;++i__)
    	{ p.push_back(static_cast<float>(x[i__]));}
    return p;
}
void Util::cl1_sphire(long int *k, long int *l, long int *m, long int *n, long int *klm2d,
	double *q, double *x, double *res, double *cu, long int *iu, double *s)
{

    long int q_dim1, q_offset, i__1, i__2;
    double d__1;

    static long int i__, j;
    static double z__;
    static long int n1, n2, ia, ii, kk, in, nk, js;
    static double sn, zu, zv;
    static long int nk1, klm, nkl, jmn, jpn;
    static double cuv;
    static long int klm1, nkl1, klm2, kode, iimn, nklm, iter;
    static float xmin;
    static double xmax;
    static long int iout;
    static double xsum;
    static long int iineg, maxit;
    static double toler;
    static float error;
    static double pivot;
    static long int kforce, iphase;
    static double tpivot;

    --s;
    --res;
    iu -= 3;
    cu -= 3;
    --x;
    q_dim1 = *klm2d;
    q_offset = 1 + q_dim1;
    q -= q_offset;

    /* Function Body */
    maxit = 500;
    kode = 0;
    toler = 1e-4f;
    iter = 0;
    n1 = *n + 1;
    n2 = *n + 2;
    nk = *n + *k;
    nk1 = nk + 1;
    nkl = nk + *l;
    nkl1 = nkl + 1;
    klm = *k + *l + *m;
    klm1 = klm + 1;
    klm2 = klm + 2;
    nklm = *n + klm;
    kforce = 1;
    iter = 0;
    js = 1;
    ia = 0;
/* SET UP LABELS IN Q. */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	q[klm2 + j * q_dim1] = (double) j;
/* L10: */
    }
    i__1 = klm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	q[i__ + n2 * q_dim1] = (double) (*n + i__);
	if (q[i__ + n1 * q_dim1] >= 0.f) {
	    goto L30;
	}
	i__2 = n2;
	for (j = 1; j <= i__2; ++j) {
	    q[i__ + j * q_dim1] = -q[i__ + j * q_dim1];
/* L20: */
	}
L30:
	;
    }
/* SET UP PHASE 1 COSTS. */
    iphase = 2;
    i__1 = nklm;
    for (j = 1; j <= i__1; ++j) {
	cu[(j << 1) + 1] = 0.f;
	cu[(j << 1) + 2] = 0.f;
	iu[(j << 1) + 1] = 0;
	iu[(j << 1) + 2] = 0;
/* L40: */
    }
    if (*l == 0) {
	goto L60;
    }
    i__1 = nkl;
    for (j = nk1; j <= i__1; ++j) {
	cu[(j << 1) + 1] = 1.f;
	cu[(j << 1) + 2] = 1.f;
	iu[(j << 1) + 1] = 1;
	iu[(j << 1) + 2] = 1;
/* L50: */
    }
    iphase = 1;
L60:
    if (*m == 0) {
	goto L80;
    }
    i__1 = nklm;
    for (j = nkl1; j <= i__1; ++j) {
	cu[(j << 1) + 2] = 1.f;
	iu[(j << 1) + 2] = 1;
	jmn = j - *n;
	if (q[jmn + n2 * q_dim1] < 0.f) {
	    iphase = 1;
	}
/* L70: */
    }
L80:
    if (kode == 0) {
	goto L150;
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	if ((d__1 = x[j]) < 0.) {
	    goto L90;
	} else if (d__1 == 0) {
	    goto L110;
	} else {
	    goto L100;
	}
L90:
	cu[(j << 1) + 1] = 1.f;
	iu[(j << 1) + 1] = 1;
	goto L110;
L100:
	cu[(j << 1) + 2] = 1.f;
	iu[(j << 1) + 2] = 1;
L110:
	;
    }
    i__1 = *k;
    for (j = 1; j <= i__1; ++j) {
	jpn = j + *n;
	if ((d__1 = res[j]) < 0.) {
	    goto L120;
	} else if (d__1 == 0) {
	    goto L140;
	} else {
	    goto L130;
	}
L120:
	cu[(jpn << 1) + 1] = 1.f;
	iu[(jpn << 1) + 1] = 1;
	if (q[j + n2 * q_dim1] > 0.f) {
	    iphase = 1;
	}
	goto L140;
L130:
	cu[(jpn << 1) + 2] = 1.f;
	iu[(jpn << 1) + 2] = 1;
	if (q[j + n2 * q_dim1] < 0.f) {
	    iphase = 1;
	}
L140:
	;
    }
L150:
    if (iphase == 2) {
	goto L500;
    }
/* COMPUTE THE MARGINAL COSTS. */
L160:
    i__1 = n1;
    for (j = js; j <= i__1; ++j) {
	xsum = 0.;
	i__2 = klm;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    ii = (long int) q[i__ + n2 * q_dim1];
	    if (ii < 0) {
		goto L170;
	    }
	    z__ = cu[(ii << 1) + 1];
	    goto L180;
L170:
	    iineg = -ii;
	    z__ = cu[(iineg << 1) + 2];
L180:
	    xsum += q[i__ + j * q_dim1] * z__;
/*  180       XSUM = XSUM + Q(I,J)*Z */
/* L190: */
	}
	q[klm1 + j * q_dim1] = xsum;
/* L200: */
    }
    i__1 = *n;
    for (j = js; j <= i__1; ++j) {
	ii = (long int) q[klm2 + j * q_dim1];
	if (ii < 0) {
	    goto L210;
	}
	z__ = cu[(ii << 1) + 1];
	goto L220;
L210:
	iineg = -ii;
	z__ = cu[(iineg << 1) + 2];
L220:
	q[klm1 + j * q_dim1] -= z__;
/* L230: */
    }
/* DETERMINE THE VECTOR TO ENTER THE BASIS. */
L240:
    xmax = 0.f;
    if (js > *n) {
	goto L490;
    }
    i__1 = *n;
    for (j = js; j <= i__1; ++j) {
	zu = q[klm1 + j * q_dim1];
	ii = (long int) q[klm2 + j * q_dim1];
	if (ii > 0) {
	    goto L250;
	}
	ii = -ii;
	zv = zu;
	zu = -zu - cu[(ii << 1) + 1] - cu[(ii << 1) + 2];
	goto L260;
L250:
	zv = -zu - cu[(ii << 1) + 1] - cu[(ii << 1) + 2];
L260:
	if (kforce == 1 && ii > *n) {
	    goto L280;
	}
	if (iu[(ii << 1) + 1] == 1) {
	    goto L270;
	}
	if (zu <= xmax) {
	    goto L270;
	}
	xmax = zu;
	in = j;
L270:
	if (iu[(ii << 1) + 2] == 1) {
	    goto L280;
	}
	if (zv <= xmax) {
	    goto L280;
	}
	xmax = zv;
	in = j;
L280:
	;
    }
    if (xmax <= toler) {
	goto L490;
    }
    if (q[klm1 + in * q_dim1] == xmax) {
	goto L300;
    }
    i__1 = klm2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	q[i__ + in * q_dim1] = -q[i__ + in * q_dim1];
/* L290: */
    }
    q[klm1 + in * q_dim1] = xmax;
/* DETERMINE THE VECTOR TO LEAVE THE BASIS. */
L300:
    if (iphase == 1 || ia == 0) {
	goto L330;
    }
    xmax = 0.f;
    i__1 = ia;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__ = (d__1 = q[i__ + in * q_dim1], abs(d__1));
	if (z__ <= xmax) {
	    goto L310;
	}
	xmax = z__;
	iout = i__;
L310:
	;
    }
    if (xmax <= toler) {
	goto L330;
    }
    i__1 = n2;
    for (j = 1; j <= i__1; ++j) {
	z__ = q[ia + j * q_dim1];
	q[ia + j * q_dim1] = q[iout + j * q_dim1];
	q[iout + j * q_dim1] = z__;
/* L320: */
    }
    iout = ia;
    --ia;
    pivot = q[iout + in * q_dim1];
    goto L420;
L330:
    kk = 0;
    i__1 = klm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__ = q[i__ + in * q_dim1];
	if (z__ <= toler) {
	    goto L340;
	}
	++kk;
	res[kk] = q[i__ + n1 * q_dim1] / z__;
	s[kk] = (double) i__;
L340:
	;
    }
L350:
    if (kk > 0) {
	goto L360;
    }
    kode = 2;
    goto L590;
L360:
    xmin = static_cast<float>( res[1] );
    iout = (long int) s[1];
    j = 1;
    if (kk == 1) {
	goto L380;
    }
    i__1 = kk;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (res[i__] >= xmin) {
	    goto L370;
	}
	j = i__;
	xmin = static_cast<float>( res[i__] );
	iout = (long int) s[i__];
L370:
	;
    }
    res[j] = res[kk];
    s[j] = s[kk];
L380:
    --kk;
    pivot = q[iout + in * q_dim1];
    ii = (long int) q[iout + n2 * q_dim1];
    if (iphase == 1) {
	goto L400;
    }
    if (ii < 0) {
	goto L390;
    }
    if (iu[(ii << 1) + 2] == 1) {
	goto L420;
    }
    goto L400;
L390:
    iineg = -ii;
    if (iu[(iineg << 1) + 1] == 1) {
	goto L420;
    }
/* 400 II = IABS(II) */
L400:
    ii = abs(ii);
    cuv = cu[(ii << 1) + 1] + cu[(ii << 1) + 2];
    if (q[klm1 + in * q_dim1] - pivot * cuv <= toler) {
	goto L420;
    }
/* BYPASS INTERMEDIATE VERTICES. */
    i__1 = n1;
    for (j = js; j <= i__1; ++j) {
	z__ = q[iout + j * q_dim1];
	q[klm1 + j * q_dim1] -= z__ * cuv;
	q[iout + j * q_dim1] = -z__;
/* L410: */
    }
    q[iout + n2 * q_dim1] = -q[iout + n2 * q_dim1];
    goto L350;
/* GAUSS-JORDAN ELIMINATION. */
L420:
    if (iter < maxit) {
	goto L430;
    }
    kode = 3;
    goto L590;
L430:
    ++iter;
    i__1 = n1;
    for (j = js; j <= i__1; ++j) {
	if (j != in) {
	    q[iout + j * q_dim1] /= pivot;
	}
/* L440: */
    }
/* IF PERMITTED, USE SUBROUTINE COL OF THE DESCRIPTION */
/* SECTION AND REPLACE THE FOLLOWING SEVEN STATEMENTS DOWN */
/* TO AND INCLUDING STATEMENT NUMBER 460 BY.. */
/*     DO 460 J=JS,N1 */
/*        IF(J .EQ. IN) GO TO 460 */
/*        Z = -Q(IOUT,J) */
/*        CALL COL(Q(1,J), Q(1,IN), Z, IOUT, KLM1) */
/* 460 CONTINUE */
    i__1 = n1;
    for (j = js; j <= i__1; ++j) {
	if (j == in) {
	    goto L460;
	}
	z__ = -q[iout + j * q_dim1];
	i__2 = klm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    if (i__ != iout) {
		q[i__ + j * q_dim1] += z__ * q[i__ + in * q_dim1];
	    }
/* L450: */
	}
L460:
	;
    }
    tpivot = -pivot;
    i__1 = klm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ != iout) {
	    q[i__ + in * q_dim1] /= tpivot;
	}
/* L470: */
    }
    q[iout + in * q_dim1] = 1.f / pivot;
    z__ = q[iout + n2 * q_dim1];
    q[iout + n2 * q_dim1] = q[klm2 + in * q_dim1];
    q[klm2 + in * q_dim1] = z__;
    ii = (long int) abs(z__);
    if (iu[(ii << 1) + 1] == 0 || iu[(ii << 1) + 2] == 0) {
	goto L240;
    }
    i__1 = klm2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__ = q[i__ + in * q_dim1];
	q[i__ + in * q_dim1] = q[i__ + js * q_dim1];
	q[i__ + js * q_dim1] = z__;
/* L480: */
    }
    ++js;
    goto L240;
/* TEST FOR OPTIMALITY. */
L490:
    if (kforce == 0) {
	goto L580;
    }
    if (iphase == 1 && q[klm1 + n1 * q_dim1] <= toler) {
	goto L500;
    }
    kforce = 0;
    goto L240;
/* SET UP PHASE 2 COSTS. */
L500:
    iphase = 2;
    i__1 = nklm;
    for (j = 1; j <= i__1; ++j) {
	cu[(j << 1) + 1] = 0.f;
	cu[(j << 1) + 2] = 0.f;
/* L510: */
    }
    i__1 = nk;
    for (j = n1; j <= i__1; ++j) {
	cu[(j << 1) + 1] = 1.f;
	cu[(j << 1) + 2] = 1.f;
/* L520: */
    }
    i__1 = klm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = (long int) q[i__ + n2 * q_dim1];
	if (ii > 0) {
	    goto L530;
	}
	ii = -ii;
	if (iu[(ii << 1) + 2] == 0) {
	    goto L560;
	}
	cu[(ii << 1) + 2] = 0.f;
	goto L540;
L530:
	if (iu[(ii << 1) + 1] == 0) {
	    goto L560;
	}
	cu[(ii << 1) + 1] = 0.f;
L540:
	++ia;
	i__2 = n2;
	for (j = 1; j <= i__2; ++j) {
	    z__ = q[ia + j * q_dim1];
	    q[ia + j * q_dim1] = q[i__ + j * q_dim1];
	    q[i__ + j * q_dim1] = z__;
/* L550: */
	}
L560:
	;
    }
    goto L160;
L570:
    if (q[klm1 + n1 * q_dim1] <= toler) {
	goto L500;
    }
    kode = 1;
    goto L590;
L580:
    if (iphase == 1) {
	goto L570;
    }
/* PREPARE OUTPUT. */
    kode = 0;
L590:
    xsum = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	x[j] = 0.f;
/* L600: */
    }
    i__1 = klm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	res[i__] = 0.f;
/* L610: */
    }
    i__1 = klm;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = (long int) q[i__ + n2 * q_dim1];
	sn = 1.f;
	if (ii > 0) {
	    goto L620;
	}
	ii = -ii;
	sn = -1.f;
L620:
	if (ii > *n) {
	    goto L630;
	}
	x[ii] = sn * q[i__ + n1 * q_dim1];
	goto L640;
L630:
	iimn = ii - *n;
	res[iimn] = sn * q[i__ + n1 * q_dim1];
	if (ii >= n1 && ii <= nk) {
	    xsum += q[i__ + n1 * q_dim1];
	}
L640:
	;
    }
    error = (float)xsum;
    return;
}

void Util::set_freq_sphire(EMData* freqvol, EMData* temp, EMData* mask, float cutoff, float freq)
{
	ENTERFUNC;
	/* Exception Handle */
	if (!freqvol) {
		throw NullPointerException("NULL input image");
	}

	int nx=freqvol->get_xsize(),ny=freqvol->get_ysize(),nz=freqvol->get_zsize();
	size_t size = (size_t)nx*ny*nz;
	float *freqvol_ptr  = freqvol->get_data();
	float *temp_ptr = temp->get_data();
	float *mask_ptr = mask->get_data();

	for (size_t i=0; i<size; ++i) {
		if(mask_ptr[i] >0.5f) {
			if(freqvol_ptr[i]  == 0.0f) {
				if(temp_ptr[i] < cutoff) freqvol_ptr[i] = freq;
			}
		}
	}

	EXITFUNC;
	freqvol->update();
}
