/**
 * $Id$
 */

/*
 * Author: Pawel A.Penczek, 09/09/2006 (Pawel.A.Penczek@uth.tmc.edu)
 * Copyright (c) 2000-2006 The University of Texas - Houston Medical School
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

#ifdef NATIVE_FFT

#include <cmath>
#include "native_fft.h" 
#include "log.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "util.h"

using namespace EMAN;

int Nativefft::fftmcf_(float *a, float *b, integer *ntot, 
            integer *n, integer *nspan, integer *isn)
{
    /* System generated locals */
    integer i__1;
    float r__1, r__2;
    int   status = -1;

    /* Builtin functions */
    // double atan(), cos(), sin(), sqrt();

    /* Local variables */
    static integer nfac[11];
    static float radf;
    static integer maxf, maxp, i__, j, k, m, kspan;
    static float c1, c2, c3;
    static integer kspnn, k1, k2, k3, k4;
    static float s1, s2, s3, aa, bb, cd, aj, c72;
    static integer jc;
    static float ck[23], ak;
    static integer jf;
    static float bk, bj;
    static integer jj;
    static float at[23], bt[23], sd;
    static integer kk;
    static float s72;
    static integer nn, np[209];
    static float sk[23];
    static integer ks, kt, nt;
    static float s120, rad, ajm, akm;
    static integer inc;
    static float ajp, akp, bkp, bkm, bjp, bjm;


/* ------- ARRAY STORAGE FOR MAXIMUM PRIME FACTOR OF 23 */

/*      EQUIVALENCE (I,II) */

/* ----- THE FOLLOWING TWO CONSTANTS SHOULD AGREE WITH */
/* ----- THE ARRAY DIMENSIONS. */

    /* Parameter adjustments */
    --b;
    --a;

    /* Function Body */
    maxf = 23;
    maxp = 209;
    if (*n < 2) {
	return 0;
    }
    inc = *isn;
    rad = atan((float)1.) * (float)8.;
    s72 = rad / (float)5.;
    c72 = cos(s72);
    s72 = sin(s72);
    s120 = sqrt((float).75);
    if (*isn >= 0) {
	goto L10;
    }
    s72 = -s72;
    s120 = -s120;
    rad = -rad;
    inc = -inc;
L10:
    nt = inc * *ntot;
    ks = inc * *nspan;
    kspan = ks;
    nn = nt - inc;
    jc = ks / *n;
    radf = rad * (float) jc * (float).5;
    i__ = 0;
    jf = 0;

/* -------DETERMINE THE FACTORS OF N------------------------- */

    m = 0;
    k = *n;
    goto L20;
L15:
    ++m;
    nfac[m - 1] = 4;
    k /= 16;
L20:
    if (k - (k / 16 << 4) == 0) {
	goto L15;
    }
    j = 3;
    jj = 9;
    goto L30;
L25:
    ++m;
    nfac[m - 1] = j;
    k /= jj;
L30:
    if (k % jj == 0) {
	goto L25;
    }
    j += 2;
/* Computing 2nd power */
    i__1 = j;
    jj = i__1 * i__1;
    if (jj <= k) {
	goto L30;
    }
    if (k > 4) {
	goto L40;
    }
    kt = m;
    nfac[m] = k;
    if (k != 1) {
	++m;
    }
    goto L80;
L40:
    if (k - (k / 4 << 2) != 0) {
	goto L50;
    }
    ++m;
    nfac[m - 1] = 2;
    k /= 4;
L50:
    kt = m;
    j = 2;
L60:
    if (k % j != 0) {
	goto L70;
    }
    ++m;
    nfac[m - 1] = j;
    k /= j;
L70:
    j = ((j + 1) / 2 << 1) + 1;
    if (j <= k) {
	goto L60;
    }
L80:
    if (kt == 0) {
	goto L100;
    }
    j = kt;
L90:
    ++m;
    nfac[m - 1] = nfac[j - 1];
    --j;
    if (j != 0) {
	goto L90;
    }

/* -------COMPUTE FOURIER TRANSFORM----------------------------- */

L100:
    sd = radf / (float) kspan;
/* Computing 2nd power */
    r__1 = sin(sd);
    cd = r__1 * r__1 * (float)2.;
    sd = sin(sd + sd);
    kk = 1;
    ++i__;
    if (nfac[i__ - 1] != 2) {
	goto L400;
    }

/* -------TRANSFORM FOR FACTOR OF 2 (INCLUDING ROTATION FACTOR)- */

    kspan /= 2;
    k1 = kspan + 2;
L210:
    k2 = kk + kspan;
    ak = a[k2];
    bk = b[k2];
    a[k2] = a[kk] - ak;
    b[k2] = b[kk] - bk;
    a[kk] += ak;
    b[kk] += bk;
    kk = k2 + kspan;
    if (kk <= nn) {
	goto L210;
    }
    kk -= nn;
    if (kk <= jc) {
	goto L210;
    }
    if (kk > kspan) {
	goto L800;
    }
L220:
    c1 = (float)1. - cd;
    s1 = sd;
L230:
    k2 = kk + kspan;
    ak = a[kk] - a[k2];
    bk = b[kk] - b[k2];
    a[kk] += a[k2];
    b[kk] += b[k2];
    a[k2] = c1 * ak - s1 * bk;
    b[k2] = s1 * ak + c1 * bk;
    kk = k2 + kspan;
    if (kk < nt) {
	goto L230;
    }
    k2 = kk - nt;
    c1 = -c1;
    kk = k1 - k2;
    if (kk > k2) {
	goto L230;
    }
    ak = c1 - (cd * c1 + sd * s1);
    s1 = sd * c1 - cd * s1 + s1;

/* ------THE FOLLOWING THREE STATEMENTS COMPENSATE FOR TRUNCATION-- */
/*      ERROR. IF ROUNDED ARITHMETIC IS USED, SUBSTITUTE */
/*      C1 = AK */

/* Computing 2nd power */
    r__1 = ak;
/* Computing 2nd power */
    r__2 = s1;
    c1 = (float).5 / (r__1 * r__1 + r__2 * r__2) + (float).5;
    s1 = c1 * s1;
    c1 *= ak;
    kk += jc;
    if (kk < k2) {
	goto L230;
    }
    k1 = k1 + inc + inc;
    kk = (k1 - kspan) / 2 + jc;
    if (kk <= jc + jc) {
	goto L220;
    }
    goto L100;

/* -------TRANSFORM FOR FACTOR OF 3 (OPTIONAL CODE)---------- */

L320:
    k1 = kk + kspan;
    k2 = k1 + kspan;
    ak = a[kk];
    bk = b[kk];
    aj = a[k1] + a[k2];
    bj = b[k1] + b[k2];
    a[kk] = ak + aj;
    b[kk] = bk + bj;
    ak = aj * (float)-.5 + ak;
    bk = bj * (float)-.5 + bk;
    aj = (a[k1] - a[k2]) * s120;
    bj = (b[k1] - b[k2]) * s120;
    a[k1] = ak - bj;
    b[k1] = bk + aj;
    a[k2] = ak + bj;
    b[k2] = bk - aj;
    kk = k2 + kspan;
    if (kk < nn) {
	goto L320;
    }
    kk -= nn;
    if (kk <= kspan) {
	goto L320;
    }
    goto L700;

/* -------TRANSFORM FOR FACTOR OF 4--------------- */

L400:
    if (nfac[i__ - 1] != 4) {
	goto L600;
    }
    kspnn = kspan;
    kspan /= 4;
L410:
    c1 = (float)1.;
    s1 = (float)0.;
L420:
    k1 = kk + kspan;
    k2 = k1 + kspan;
    k3 = k2 + kspan;
    akp = a[kk] + a[k2];
    akm = a[kk] - a[k2];
    ajp = a[k1] + a[k3];
    ajm = a[k1] - a[k3];
    a[kk] = akp + ajp;
    ajp = akp - ajp;
    bkp = b[kk] + b[k2];
    bkm = b[kk] - b[k2];
    bjp = b[k1] + b[k3];
    bjm = b[k1] - b[k3];
    b[kk] = bkp + bjp;
    bjp = bkp - bjp;
    if (*isn < 0) {
	goto L450;
    }
    akp = akm - bjm;
    akm += bjm;
    bkp = bkm + ajm;
    bkm -= ajm;
    if (s1 == (float)0.) {
	goto L460;
    }
L430:
    a[k1] = akp * c1 - bkp * s1;
    b[k1] = akp * s1 + bkp * c1;
    a[k2] = ajp * c2 - bjp * s2;
    b[k2] = ajp * s2 + bjp * c2;
    a[k3] = akm * c3 - bkm * s3;
    b[k3] = akm * s3 + bkm * c3;
    kk = k3 + kspan;
    if (kk <= nt) {
	goto L420;
    }
L440:
    c2 = c1 - (cd * c1 + sd * s1);
    s1 = sd * c1 - cd * s1 + s1;

/* -------THE FOLLOWING THREE STATEMENTS COMPENSATE FOR TRUNCATION */
/*       ERROR. IF ROUNDED ARITHMETIC IS USED, SUBSTITUTE */
/*       C1 = C2 */

/* Computing 2nd power */
    r__1 = c2;
/* Computing 2nd power */
    r__2 = s1;
    c1 = (float).5 / (r__1 * r__1 + r__2 * r__2) + (float).5;
    s1 = c1 * s1;
    c1 *= c2;
/* Computing 2nd power */
    r__1 = c1;
/* Computing 2nd power */
    r__2 = s1;
    c2 = r__1 * r__1 - r__2 * r__2;
    s2 = c1 * (float)2. * s1;
    c3 = c2 * c1 - s2 * s1;
    s3 = c2 * s1 + s2 * c1;
    kk = kk - nt + jc;
    if (kk <= kspan) {
	goto L420;
    }
    kk = kk - kspan + inc;
    if (kk <= jc) {
	goto L410;
    }
    if (kspan == jc) {
	goto L800;
    }
    goto L100;
L450:
    akp = akm + bjm;
    akm -= bjm;
    bkp = bkm - ajm;
    bkm += ajm;
    if (s1 != (float)0.) {
	goto L430;
    }
L460:
    a[k1] = akp;
    b[k1] = bkp;
    a[k2] = ajp;
    b[k2] = bjp;
    a[k3] = akm;
    b[k3] = bkm;
    kk = k3 + kspan;
    if (kk <= nt) {
	goto L420;
    }
    goto L440;

/* -------TRANSFORM FOR FACTOR OF 5 (OPTIONAL CODE)-------- */

L510:
/* Computing 2nd power */
    r__1 = c72;
/* Computing 2nd power */
    r__2 = s72;
    c2 = r__1 * r__1 - r__2 * r__2;
    s2 = c72 * (float)2. * s72;
L520:
    k1 = kk + kspan;
    k2 = k1 + kspan;
    k3 = k2 + kspan;
    k4 = k3 + kspan;
    akp = a[k1] + a[k4];
    akm = a[k1] - a[k4];
    bkp = b[k1] + b[k4];
    bkm = b[k1] - b[k4];
    ajp = a[k2] + a[k3];
    ajm = a[k2] - a[k3];
    bjp = b[k2] + b[k3];
    bjm = b[k2] - b[k3];
    aa = a[kk];
    bb = b[kk];
    a[kk] = aa + akp + ajp;
    b[kk] = bb + bkp + bjp;
    ak = akp * c72 + ajp * c2 + aa;
    bk = bkp * c72 + bjp * c2 + bb;
    aj = akm * s72 + ajm * s2;
    bj = bkm * s72 + bjm * s2;
    a[k1] = ak - bj;
    a[k4] = ak + bj;
    b[k1] = bk + aj;
    b[k4] = bk - aj;
    ak = akp * c2 + ajp * c72 + aa;
    bk = bkp * c2 + bjp * c72 + bb;
    aj = akm * s2 - ajm * s72;
    bj = bkm * s2 - bjm * s72;
    a[k2] = ak - bj;
    a[k3] = ak + bj;
    b[k2] = bk + aj;
    b[k3] = bk - aj;
    kk = k4 + kspan;
    if (kk < nn) {
	goto L520;
    }
    kk -= nn;
    if (kk <= kspan) {
	goto L520;
    }
    goto L700;

/* -------TRANSFORM FOR ODD FACTORS-------------------- */

L600:
    k = nfac[i__ - 1];
    kspnn = kspan;
    kspan /= k;
    if (k == 3) {
	goto L320;
    }
    if (k == 5) {
	goto L510;
    }
    if (k == jf) {
	goto L640;
    }
    jf = k;
    s1 = rad / (float) k;
    c1 = cos(s1);
    s1 = sin(s1);
    if (jf > maxf) {
	goto L998;
    }
    ck[jf - 1] = (float)1.;
    sk[jf - 1] = (float)0.;
    j = 1;
L630:
    ck[j - 1] = ck[k - 1] * c1 + sk[k - 1] * s1;
    sk[j - 1] = ck[k - 1] * s1 - sk[k - 1] * c1;
    --k;
    ck[k - 1] = ck[j - 1];
    sk[k - 1] = -sk[j - 1];
    ++j;
    if (j < k) {
	goto L630;
    }
L640:
    k1 = kk;
    k2 = kk + kspnn;
    aa = a[kk];
    bb = b[kk];
    ak = aa;
    bk = bb;
    j = 1;
    k1 += kspan;
L650:
    k2 -= kspan;
    ++j;
    at[j - 1] = a[k1] + a[k2];
    ak = at[j - 1] + ak;
    bt[j - 1] = b[k1] + b[k2];
    bk = bt[j - 1] + bk;
    ++j;
    at[j - 1] = a[k1] - a[k2];
    bt[j - 1] = b[k1] - b[k2];
    k1 += kspan;
    if (k1 < k2) {
	goto L650;
    }
    a[kk] = ak;
    b[kk] = bk;
    k1 = kk;
    k2 = kk + kspnn;
    j = 1;
L660:
    k1 += kspan;
    k2 -= kspan;
    jj = j;
    ak = aa;
    bk = bb;
    aj = (float)0.;
    bj = (float)0.;
    k = 1;
L670:
    ++k;
    ak = at[k - 1] * ck[jj - 1] + ak;
    bk = bt[k - 1] * ck[jj - 1] + bk;
    ++k;
    aj = at[k - 1] * sk[jj - 1] + aj;
    bj = bt[k - 1] * sk[jj - 1] + bj;
    jj += j;
    if (jj > jf) {
	jj -= jf;
    }
    if (k < jf) {
	goto L670;
    }
    k = jf - j;
    a[k1] = ak - bj;
    b[k1] = bk + aj;
    a[k2] = ak + bj;
    b[k2] = bk - aj;
    ++j;
    if (j < k) {
	goto L660;
    }
    kk += kspnn;
    if (kk <= nn) {
	goto L640;
    }
    kk -= nn;
    if (kk <= kspan) {
	goto L640;
    }

/* ------- MULTIPLY BY ROTATION FACTOR (EXCEPT FOR FACTORS OF 2 AND 4) */

L700:
    if (i__ == m) {
	goto L800;
    }
    kk = jc + 1;
L710:
    c2 = (float)1. - cd;
    s1 = sd;
L720:
    c1 = c2;
    s2 = s1;
    kk += kspan;
L730:
    ak = a[kk];
    a[kk] = c2 * ak - s2 * b[kk];
    b[kk] = s2 * ak + c2 * b[kk];
    kk += kspnn;
    if (kk <= nt) {
	goto L730;
    }
    ak = s1 * s2;
    s2 = s1 * c2 + c1 * s2;
    c2 = c1 * c2 - ak;
    kk = kk - nt + kspan;
    if (kk <= kspnn) {
	goto L730;
    }
    c2 = c1 - (cd * c1 + sd * s1);
    s1 += sd * c1 - cd * s1;

/* -------THE FOLLOWING THREE STATEMENTS COMPENSATE FOR TRUNCATION */
/*       ERROR. IF ROUNDED ARITHMETIC IS USED, THEY MAY */
/*       BE DELETED. */

/* Computing 2nd power */
    r__1 = c2;
/* Computing 2nd power */
    r__2 = s1;
    c1 = (float).5 / (r__1 * r__1 + r__2 * r__2) + (float).5;
    s1 = c1 * s1;
    c2 = c1 * c2;
    kk = kk - kspnn + jc;
    if (kk <= kspan) {
	goto L720;
    }
    kk = kk - kspan + jc + inc;
    if (kk <= jc + jc) {
	goto L710;
    }
    goto L100;

/* -------PERMUTE THE RESULTS TO NORMAL ORDER---DONE IN TWO STAGES- */
/*       PERMUTATION FOR SQUARE FACTORS OF N */

L800:
    np[0] = ks;
    if (kt == 0) {
	goto L890;
    }
    k = kt + kt + 1;
    if (m < k) {
	--k;
    }
    j = 1;
    np[k] = jc;
L810:
    np[j] = np[j - 1] / nfac[j - 1];
    np[k - 1] = np[k] * nfac[j - 1];
    ++j;
    --k;
    if (j < k) {
	goto L810;
    }
    k3 = np[k];
    kspan = np[1];
    kk = jc + 1;
    j = 1;
    k2 = kspan + 1;
    if (*n != *ntot) {
	goto L850;
    }

/* -------PERMUTATION FOR SINGLE-VARIATE TRANSFORM (OPTIONAL CODE) */

L820:
    ak = a[kk];
    a[kk] = a[k2];
    a[k2] = ak;
    bk = b[kk];
    b[kk] = b[k2];
    b[k2] = bk;
    kk += inc;
    k2 = kspan + k2;
    if (k2 < ks) {
	goto L820;
    }
L830:
    k2 -= np[j - 1];
    ++j;
    k2 = np[j] + k2;
    if (k2 > np[j - 1]) {
	goto L830;
    }
    j = 1;
L840:
    if (kk < k2) {
	goto L820;
    }
    kk += inc;
    k2 = kspan + k2;
    if (k2 < ks) {
	goto L840;
    }
    if (kk < ks) {
	goto L830;
    }
    jc = k3;
    goto L890;

/* -------PERMUTATION FOR MULTIVARIATE TRANSFORM------ */

L850:
    k = kk + jc;
L860:
    ak = a[kk];
    a[kk] = a[k2];
    a[k2] = ak;
    bk = b[kk];
    b[kk] = b[k2];
    b[k2] = bk;
    kk += inc;
    k2 += inc;
    if (kk < k) {
	goto L860;
    }
    kk = kk + ks - jc;
    k2 = k2 + ks - jc;
    if (kk < nt) {
	goto L850;
    }
    k2 = k2 - nt + kspan;
    kk = kk - nt + jc;
    if (k2 < ks) {
	goto L850;
    }
L870:
    k2 -= np[j - 1];
    ++j;
    k2 = np[j] + k2;
    if (k2 > np[j - 1]) {
	goto L870;
    }
    j = 1;
L880:
    if (kk < k2) {
	goto L850;
    }
    kk += jc;
    k2 = kspan + k2;
    if (k2 < ks) {
	goto L880;
    }
    if (kk < ks) {
	goto L870;
    }
    jc = k3;
L890:
    if ((kt << 1) + 1 >= m) {
	return 0;
    }
    kspnn = np[kt];

/* -------PERMUTATION FOR SQUARE-FREE FACTORS OF N----- */

    j = m - kt;
    nfac[j] = 1;
L900:
    nfac[j - 1] *= nfac[j];
    --j;
    if (j != kt) {
	goto L900;
    }
    ++kt;
    nn = nfac[kt - 1] - 1;
    if (nn > maxp) {
	goto L998;
    }
    jj = 0;
    j = 0;
    goto L906;
L902:
    jj -= k2;
    k2 = kk;
    ++k;
    kk = nfac[k - 1];
L904:
    jj = kk + jj;
    if (jj >= k2) {
	goto L902;
    }
    np[j - 1] = jj;
L906:
    k2 = nfac[kt - 1];
    k = kt + 1;
    kk = nfac[k - 1];
    ++j;
    if (j <= nn) {
	goto L904;
    }

/* -------DETERMINE THE PERMUTATION CYCLES OF LENGTH GREATER THAN 1 */

    j = 0;
    goto L914;
L910:
    k = kk;
    kk = np[k - 1];
    np[k - 1] = -kk;
    if (kk != j) {
	goto L910;
    }
    k3 = kk;
L914:
    ++j;
    kk = np[j - 1];
    if (kk < 0) {
	goto L914;
    }
    if (kk != j) {
	goto L910;
    }
    np[j - 1] = -j;
    if (j != nn) {
	goto L914;
    }
    maxf = inc * maxf;

/* -------REORDER A AND B, FOLLOWING THE PERMUTATION CYCLES- */

    goto L950;
L924:
    --j;
    if (np[j - 1] < 0) {
	goto L924;
    }
    jj = jc;
L926:
    kspan = jj;
    if (jj > maxf) {
	kspan = maxf;
    }
    jj -= kspan;
    k = np[j - 1];
    kk = jc * k + i__ + jj;
    k1 = kk + kspan;
    k2 = 0;
L928:
    ++k2;
    at[k2 - 1] = a[k1];
    bt[k2 - 1] = b[k1];
    k1 -= inc;
    if (k1 != kk) {
	goto L928;
    }
L932:
    k1 = kk + kspan;
    k2 = k1 - jc * (k + np[k - 1]);
    k = -np[k - 1];
L936:
    a[k1] = a[k2];
    b[k1] = b[k2];
    k1 -= inc;
    k2 -= inc;
    if (k1 != kk) {
	goto L936;
    }
    kk = k2;
    if (k != j) {
	goto L932;
    }
    k1 = kk + kspan;
    k2 = 0;
L940:
    ++k2;
    a[k1] = at[k2 - 1];
    b[k1] = bt[k2 - 1];
    k1 -= inc;
    if (k1 != kk) {
	goto L940;
    }
    if (jj != 0) {
	goto L926;
    }
    if (j != 1) {
	goto L924;
    }
L950:
    j = k3 + 1;
    nt -= kspnn;
    i__ = nt - inc + 1;
    if (nt >= 0) {
	goto L924;
    }
    return 0;

/* --------ERROR FINISH, INSUFFICIENT ARRAY STORAGE-- */

L998:
    *isn = 0;
    return status;
} /* fftmcf_ */

#define work(i)  work[(i)-1]
#define x(i)     x[(i)-1]
#define y(i)     y[(i)-1]

// 1d inplace FFT
int Nativefft::fmrs_1rf(float *x, float *work, int nsam)
{
	// input : x(nsam+2-nsam%2), work(nsam+2-nsam%2)
	// output: x(nsam+2-nsam%2), overwritten by Fourier coefficients

	int i, n, status;
	int inv=1;

	n = nsam;

	for (i=1; i<=n; i++)  work(i) = 0.0;
	status = fftmcf_(x,work,&n,&n,&n,&inv);
	// should put in appropriate exception handling here
	if (status == -1) {
		fprintf(stderr, "Native FFT cannot be performed on images with one of ");
		fprintf(stderr, "the dimensions set to n = %d\n", nsam);  
		exit(1); 
	}

	// check even/odd
	if (n%2 == 0) {
		for (i=n+1;i>=3;i=i-2) {
			x(i)   = x((i+1)/2);
			x(i+1) = work((i+1)/2);
		}
		x(2)   = 0.0;
		x(n+2) = 0.0;
	}
	else {
		for (i=n;i>=3;i=i-2) {
			x(i)   = x(i/2+1);
			x(i+1) = work(i/2+1);
		}
		x(2)=0.0;
	}
	return status;
}

// 1d inplace IFFT
int Nativefft::fmrs_1rb(float *x, float *work, int nsam)
{
	// input:  x(nsam+2-nsam%2), work(nsam+2-nsam%2)
	// output: x(nsam+2-nsam%2), overwritten with Fourier coefficients

	int i, n, status;
	int inv=-1;

	n = nsam; 

	for (i=2;i<=n/2+1;i++) {
		work(i)     = x(2*i)/n;
		work(n-i+2) = -work(i);
	}
	work(1) = 0.0;

	for (i=1;i<=n/2+1;i++)  x(i) = x(2*i-1)/n;
	for (i=n;i>=n/2+2;i--)  x(i) = x(n-i+2);

	status = fftmcf_(x,work,&n,&n,&n,&inv);
	// should put in appropriate exception handling here
	if (status == -1) {
		fprintf(stderr, "Native IFT cannot be performed on images with one of ");
		fprintf(stderr, "the dimensions set to n = %d\n", nsam);  
		exit(1); 
	}

	return status;
}

#undef x
#undef y
#undef work


//----------2D FFT INPLACE--------------------

#define x(i,j) x[((j)-1)*lda + (i)-1]
#define y(i,j) y[((j)-1)*lda + (i)-1]

// 2d inplace fft
int Nativefft::fmrs_2rf(float *x, float *work, int lda, int nsam, int nrow)
{
	// input:  x(lda,nrow), work(lda)
	// output: x(lda,nrow) overwritten with Fourier coefficients

	int   ins, status=0, i, j;

	ins = lda;
	/*
	int  l;
	l=(int)(log2(nsam));
	for (j=1;j<=nrow;j++) {
		Util::fftr_q(&x(1,j),l);
		status = fmrs_1rf(&x(1,j), work, nsam);
	}
	*/

	for (i=1;i<=lda;i=i+2){
		status = fftmcf_(&x(i,1),&x(i+1,1),&nrow,&nrow,&nrow,&ins);
		if (status == -1) {
			fprintf(stderr, "Native FFT cannot be performed on images with one of ");
			fprintf(stderr, "the dimensions set to n = %d\n", nrow);
			exit(1); 
		}
	}

	return status;
}


// 2d inplace IFFT
int Nativefft::fmrs_2rb(float *y, float *work, int lda, int nsam, int nrow)
{
	// input:  y(lda,nrow), work(lda)
	// output: y(lda,nrow), overwritten with Fourier coefficients

	int   ins, status=0, i, j;
	float q;

	ins=-lda;

	for (i=1; i<=lda; i=i+2) {
		fftmcf_(&y(i,1),&y(i+1,1),&nrow,&nrow,&nrow,&ins);
		if (status == -1) {
			fprintf(stderr, "Native IFT cannot be performed on images with one of ");
			fprintf(stderr, "the dimensions set to n = %d\n", nrow);
			exit(1); 
		}
	}
	 
	// normalize for inverse
	q = 1.0/(float)(nrow);
	for (j=1; j<=nrow; j++)  for (i=1; i<=lda; i++) y(i,j)*=q;

	// need an inplace 1d ifft 
	for (j=1; j<=nrow; j++) status = fmrs_1rb(&y(1,j), work, nsam);

	return status;
}
#undef x
#undef y

//----------2D FFT out of place--------------------

#define complexd(i,j) complexd[(j-1)*lda + i-1]
#define reald(i,j)    reald[(j-1)*nsam + i-1]

// 2D out of place fft
int Nativefft::ftp_2rf(float *reald, float *complexd, int lda, int nsam, int nrow)
{
	// input:  reald(nsam,nrow), 
	//  work(2*nsam+15)
	// output: complexd(lda,nrow) Fourier coefficients

	int   ins, status=0, i, j;

	float * work = (float*) malloc((2*nsam+15)*sizeof(float));
	if (!work) {
	 	 fprintf(stderr,"real_to_complex_nd(2df): failed to allocate work\n");
	 	 LOGERR("real_to_complex_nd(2df): failed to allocate work\n");
	}
	//  Do columns with fftpack	
	rffti(nsam, work);

	for (j=1; j<=nrow; j++)  {
		memcpy(&complexd(2,j), &reald(1,j), nsam * sizeof(float));
		rfftf(nsam, &complexd(2,j), work);
		complexd(1,j) = complexd(2,j) ;
		complexd(2,j) = 0.0f ;
		if (nsam%2 == 0)  complexd(nsam+2,j) = 0.0f ;
	}

	free(work);

	ins = lda;

	for (i=1; i<=lda; i=i+2) {
	   status = fftmcf_(&complexd(i,1),&complexd(i+1,1),&nrow,&nrow,&nrow,&ins);
	   if (status == -1) {
		  fprintf(stderr, "Native FFT cannot be performed on images with one of ");
		  fprintf(stderr, "the dimensions set to n = %d\n", nrow);
		  exit(1); 
	   }
	}

   return status;
}


int Nativefft::ftp_2rb(float *complexd, int lda, int nsam, int nrow)
{
	// input:  complexd(lsd,nrow), 
	//  work(2*nsam+15)
	// output: complexd(lda,nrow) Fourier coefficients

	int   ins, status=0, i, j;

	ins = -lda;

	for (i=1; i<=lda; i=i+2) {
		status = fftmcf_(&complexd(i,1),&complexd(i+1,1),&nrow,&nrow,&nrow,&ins);
		if (status == -1) {
			fprintf(stderr, "Native FFT cannot be performed on images with one of ");
			fprintf(stderr, "the dimensions set to n = %d\n", nrow);
			exit(1); 
		}
	}


	float * work = (float*) malloc((2*nsam+15)*sizeof(float));
	if (!work) {
		fprintf(stderr,"real_to_complex_nd(2df): failed to allocate work\n");
		LOGERR("real_to_complex_nd(2df): failed to allocate work\n");
	}
	//  Do inv columns with fftpack	
	rffti(nsam, work);

	for (j=1; j<=nrow; j++) {
		memmove(&complexd(2,j), &complexd(3,j), (nsam-1) * sizeof(float));
		rfftb(nsam, &complexd(1,j), work);
	}
	free(work);
	//  Normalize
	float nrm = 1.0f/float(nsam*nrow);
	for (int j = 1; j<=nrow; j++) for (int i = 1; i<=nsam; i++) complexd(i,j) *= nrm;

   return status;
}
#undef complexd
#undef reald

#define  reald(i,j,k)     reald[i-1 + (j-1+(k-1)*nrow)*nsam]
#define  complexd(i,j,k)  complexd[i-1 + (j-1+(k-1)*nrow)*lda]
// 3D out of place FFT
int Nativefft::ftp_3rf(float *reald, float *complexd, int lda, int nsam, int nrow, int nslice)
{
	// input :  reald(nsam,nrow,nslice)
	// output:  complexd(lda,nrow,nslice), overwritten with Fourier coefficients

	int ndr, i, j, k, status=0;

	for (k=1; k<=nslice;k ++) status = ftp_2rf(&reald(1,1,k), &complexd(1,1,k), lda, nsam, nrow);

	ndr=lda*nrow;

	for (j=1; j<=nrow; j++) {
		for (i=1; i<=lda-1; i=i+2) {
			status = fftmcf_(&complexd(i,j,1),&complexd(i+1,j,1),&nslice,&nslice,&nslice,&ndr);
			if (status == -1) {
				fprintf(stderr, "Native FFT cannot be performed on images with one of ");
				fprintf(stderr, "the dimensions set to n = %d\n", nslice);
				exit(1); 
			}
		}
	}

	return status;
}
#undef reald

// 3d out of place IFFT
int Nativefft::ftp_3rb(float *complexd, int lda, int nsam, int nrow, int nslice)
{
	// input:  complexd(lda,nrow,nslice)
	// output: complexd(lda,nrow,nslice), with Fourier coefficients 

	int ndr, i, j, k, status=0;
	float q;

	ndr=-lda*nrow;

	for (j=1; j<=nrow; j++) {
		for (i=1; i<=lda-1; i=i+2) {
			status = fftmcf_(&complexd(i,j,1),&complexd(i+1,j,1),&nslice,&nslice,&nslice,&ndr);
			if (status == -1) {
				fprintf(stderr, "Native IFT cannot be performed on images with one of ");
				fprintf(stderr, "the dimensions set to n = %d\n", nslice);
				exit(1); 
			}
		}
	}

	// normalize for inverse
	q=1.0/(float)(nslice);
	for (k=1; k<=nslice; k++) for (j=1;j<=nrow;j++) for (i=1;i<=lda;i++) complexd(i,j,k) *= q;

	for (k=1; k<=nslice; k++) status = ftp_2rb(&complexd(1,1,k), lda, nsam, nrow);

	return status;
}
#undef complexd

//---------------3D INPLACE FFT----------------------
#define b(i,j,k) b[((k)-1)*lda*nrow + ((j)-1)*lda + (i) - 1]

// 3d inplace FFT
int Nativefft::fmrs_3rf(float *b, float *work, int lda, int nsam, int nrow, int nslice)
{
	// input :  b(lda,nrow,nslice), work(lda)
	// output:  b(lda,nrow,nslice), overwritten with Fourier coefficients

	int ndr, i, j, k, status=0;

	ndr=lda*nrow;

	for (k=1; k<=nslice; k++) status = fmrs_2rf(&b(1,1,k), work, lda, nsam, nrow);

	for (j=1; j<=nrow; j++) {
		for (i=1; i<=lda-1; i=i+2) {
			status = fftmcf_(&b(i,j,1),&b(i+1,j,1),&nslice,&nslice,&nslice,&ndr);
			if (status == -1) {
				fprintf(stderr, "Native FFT cannot be performed on images with one of ");
				fprintf(stderr, "the dimensions set to n = %d\n", nslice);
				exit(1); 
			}
		}
	}
	// should check for error here by looking at ndrt

	return status;
}

// 3d inplace IFFT
int Nativefft::fmrs_3rb(float *b, float *work, int lda, int nsam, int nrow, int nslice)
{
	// input:  b(lda,nrow,nslice), work(lda)
	// output: b(lda,nrow,nslice), overwritten with Fourier coefficients 

	int ndr, i, j, k, status=0;
	float q;

	ndr=-lda*nrow;

	for (j=1;j<=nrow;j++) {
		for (i=1;i<=lda-1;i=i+2) {
			status = fftmcf_(&b(i,j,1),&b(i+1,j,1),&nslice,&nslice,&nslice,&ndr);
			if (status == -1) {
				fprintf(stderr, "Native IFT cannot be performed on images with one of ");
				fprintf(stderr, "the dimensions set to n = %d\n", nslice);
				exit(1); 
			}
		}
	}

	// should check for error here by looking at ndrt

	// normalize for inverse
	q=1.0/(float)(nslice);
	for (k=1;k<=nslice;k++) for (j=1;j<=nrow;j++) for (i=1;i<=lda;i++) b(i,j,k)*=q;

	for (k=1; k<=nslice; k++) status = fmrs_2rb(&b(1,1,k), work, lda, nsam, nrow);
	return status;
}
#undef b



/*
fftpack.c : A set of FFT routines in C.
Algorithmically based on Fortran-77 FFTPACK by Paul N. Swarztrauber (Version 4, 1985).

*/

/* isign is +1 for backward and -1 for forward transforms */

#include <math.h>
#include <stdio.h>
//#define DOUBLE

#ifdef DOUBLE
#define Treal double
#else
#define Treal float
#endif


#define ref(u,a) u[a]

#define MAXFAC 13    /* maximum number of factors in factorization of n */
#define NSPECIAL 4   /* number of factors for which we have special-case routines */

//#define __cplusplus
//#ifdef __cplusplus
//extern "C" {
//#endif


/* ----------------------------------------------------------------------
   passf2, passf3, passf4, passf5, passf. Complex FFT passes fwd and bwd.
---------------------------------------------------------------------- */

static void passf2(int ido, int l1, const Treal cc[], Treal ch[], const Treal wa1[], int isign)
  /* isign==+1 for backward transform */
  {
	int i, k, ah, ac;
	Treal ti2, tr2;
	if (ido <= 2) {
		for (k=0; k<l1; k++) {
			ah = k*ido;
			ac = 2*k*ido;
			ch[ah]  	    = ref(cc,ac) + ref(cc,ac + ido);
			ch[ah + ido*l1]     = ref(cc,ac) - ref(cc,ac + ido);
			ch[ah+1]	    = ref(cc,ac+1) + ref(cc,ac + ido + 1);
			ch[ah + ido*l1 + 1] = ref(cc,ac+1) - ref(cc,ac + ido + 1);
		}
	} else {
		for (k=0; k<l1; k++) {
			for (i=0; i<ido-1; i+=2) {
				ah = i + k*ido;
				ac = i + 2*k*ido;
				ch[ah]   = ref(cc,ac) + ref(cc,ac + ido);
				tr2	 = ref(cc,ac) - ref(cc,ac + ido);
				ch[ah+1] = ref(cc,ac+1) + ref(cc,ac + 1 + ido);
				ti2	 = ref(cc,ac+1) - ref(cc,ac + 1 + ido);
				ch[ah+l1*ido+1] = wa1[i]*ti2 + isign*wa1[i+1]*tr2;
				ch[ah+l1*ido]	= wa1[i]*tr2 - isign*wa1[i+1]*ti2;
			}
		}
	}
  } /* passf2 */


static void passf3(int ido, int l1, const Treal cc[], Treal ch[],
      const Treal wa1[], const Treal wa2[], int isign)
  /* isign==+1 for backward transform */
  {
    static const Treal taur = -0.5;
    static const Treal taui = 0.866025403784439;
    int i, k, ac, ah;
    Treal ci2, ci3, di2, di3, cr2, cr3, dr2, dr3, ti2, tr2;
    if (ido == 2) {
      for (k=1; k<=l1; k++) {
        ac = (3*k - 2)*ido;
        tr2 = ref(cc,ac) + ref(cc,ac + ido);
        cr2 = ref(cc,ac - ido) + taur*tr2;
        ah = (k - 1)*ido;
        ch[ah] = ref(cc,ac - ido) + tr2;

        ti2 = ref(cc,ac + 1) + ref(cc,ac + ido + 1);
        ci2 = ref(cc,ac - ido + 1) + taur*ti2;
        ch[ah + 1] = ref(cc,ac - ido + 1) + ti2;

        cr3 = isign*taui*(ref(cc,ac) - ref(cc,ac + ido));
        ci3 = isign*taui*(ref(cc,ac + 1) - ref(cc,ac + ido + 1));
        ch[ah + l1*ido] = cr2 - ci3;
        ch[ah + 2*l1*ido] = cr2 + ci3;
        ch[ah + l1*ido + 1] = ci2 + cr3;
        ch[ah + 2*l1*ido + 1] = ci2 - cr3;
      }
    } else {
      for (k=1; k<=l1; k++) {
        for (i=0; i<ido-1; i+=2) {
          ac = i + (3*k - 2)*ido;
          tr2 = ref(cc,ac) + ref(cc,ac + ido);
          cr2 = ref(cc,ac - ido) + taur*tr2;
          ah = i + (k-1)*ido;
          ch[ah] = ref(cc,ac - ido) + tr2;
          ti2 = ref(cc,ac + 1) + ref(cc,ac + ido + 1);
          ci2 = ref(cc,ac - ido + 1) + taur*ti2;
          ch[ah + 1] = ref(cc,ac - ido + 1) + ti2;
          cr3 = isign*taui*(ref(cc,ac) - ref(cc,ac + ido));
          ci3 = isign*taui*(ref(cc,ac + 1) - ref(cc,ac + ido + 1));
          dr2 = cr2 - ci3;
          dr3 = cr2 + ci3;
          di2 = ci2 + cr3;
          di3 = ci2 - cr3;
          ch[ah + l1*ido + 1] = wa1[i]*di2 + isign*wa1[i+1]*dr2;
          ch[ah + l1*ido] = wa1[i]*dr2 - isign*wa1[i+1]*di2;
          ch[ah + 2*l1*ido + 1] = wa2[i]*di3 + isign*wa2[i+1]*dr3;
          ch[ah + 2*l1*ido] = wa2[i]*dr3 - isign*wa2[i+1]*di3;
        }
      }
    }
  } /* passf3 */


static void passf4(int ido, int l1, const Treal cc[], Treal ch[],
      const Treal wa1[], const Treal wa2[], const Treal wa3[], int isign)
  /* isign == -1 for forward transform and +1 for backward transform */
  {
    int i, k, ac, ah;
    Treal ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
    if (ido == 2) {
      for (k=0; k<l1; k++) {
        ac = 4*k*ido + 1;
        ti1 = ref(cc,ac) - ref(cc,ac + 2*ido);
        ti2 = ref(cc,ac) + ref(cc,ac + 2*ido);
        tr4 = ref(cc,ac + 3*ido) - ref(cc,ac + ido);
        ti3 = ref(cc,ac + ido) + ref(cc,ac + 3*ido);
        tr1 = ref(cc,ac - 1) - ref(cc,ac + 2*ido - 1);
        tr2 = ref(cc,ac - 1) + ref(cc,ac + 2*ido - 1);
        ti4 = ref(cc,ac + ido - 1) - ref(cc,ac + 3*ido - 1);
        tr3 = ref(cc,ac + ido - 1) + ref(cc,ac + 3*ido - 1);
        ah = k*ido;
        ch[ah] = tr2 + tr3;
        ch[ah + 2*l1*ido] = tr2 - tr3;
        ch[ah + 1] = ti2 + ti3;
        ch[ah + 2*l1*ido + 1] = ti2 - ti3;
        ch[ah + l1*ido] = tr1 + isign*tr4;
        ch[ah + 3*l1*ido] = tr1 - isign*tr4;
        ch[ah + l1*ido + 1] = ti1 + isign*ti4;
        ch[ah + 3*l1*ido + 1] = ti1 - isign*ti4;
      }
    } else {
      for (k=0; k<l1; k++) {
        for (i=0; i<ido-1; i+=2) {
          ac = i + 1 + 4*k*ido;
          ti1 = ref(cc,ac) - ref(cc,ac + 2*ido);
          ti2 = ref(cc,ac) + ref(cc,ac + 2*ido);
          ti3 = ref(cc,ac + ido) + ref(cc,ac + 3*ido);
          tr4 = ref(cc,ac + 3*ido) - ref(cc,ac + ido);
          tr1 = ref(cc,ac - 1) - ref(cc,ac + 2*ido - 1);
          tr2 = ref(cc,ac - 1) + ref(cc,ac + 2*ido - 1);
          ti4 = ref(cc,ac + ido - 1) - ref(cc,ac + 3*ido - 1);
          tr3 = ref(cc,ac + ido - 1) + ref(cc,ac + 3*ido - 1);
          ah = i + k*ido;
          ch[ah] = tr2 + tr3;
          cr3 = tr2 - tr3;
          ch[ah + 1] = ti2 + ti3;
          ci3 = ti2 - ti3;
          cr2 = tr1 + isign*tr4;
          cr4 = tr1 - isign*tr4;
          ci2 = ti1 + isign*ti4;
          ci4 = ti1 - isign*ti4;
          ch[ah + l1*ido] = wa1[i]*cr2 - isign*wa1[i + 1]*ci2;
          ch[ah + l1*ido + 1] = wa1[i]*ci2 + isign*wa1[i + 1]*cr2;
          ch[ah + 2*l1*ido] = wa2[i]*cr3 - isign*wa2[i + 1]*ci3;
          ch[ah + 2*l1*ido + 1] = wa2[i]*ci3 + isign*wa2[i + 1]*cr3;
          ch[ah + 3*l1*ido] = wa3[i]*cr4 -isign*wa3[i + 1]*ci4;
          ch[ah + 3*l1*ido + 1] = wa3[i]*ci4 + isign*wa3[i + 1]*cr4;
        }
      }
    }
  } /* passf4 */


static void passf5(int ido, int l1, const Treal cc[], Treal ch[],
      const Treal wa1[], const Treal wa2[], const Treal wa3[], const Treal wa4[], int isign)
  /* isign == -1 for forward transform and +1 for backward transform */
  {
    static const Treal tr11 = 0.309016994374947;
    static const Treal ti11 = 0.951056516295154;
    static const Treal tr12 = -0.809016994374947;
    static const Treal ti12 = 0.587785252292473;
    int i, k, ac, ah;
    Treal ci2, ci3, ci4, ci5, di3, di4, di5, di2, cr2, cr3, cr5, cr4, ti2, ti3,
        ti4, ti5, dr3, dr4, dr5, dr2, tr2, tr3, tr4, tr5;
    if (ido == 2) {
      for (k = 1; k <= l1; ++k) {
        ac = (5*k - 4)*ido + 1;
        ti5 = ref(cc,ac) - ref(cc,ac + 3*ido);
        ti2 = ref(cc,ac) + ref(cc,ac + 3*ido);
        ti4 = ref(cc,ac + ido) - ref(cc,ac + 2*ido);
        ti3 = ref(cc,ac + ido) + ref(cc,ac + 2*ido);
        tr5 = ref(cc,ac - 1) - ref(cc,ac + 3*ido - 1);
        tr2 = ref(cc,ac - 1) + ref(cc,ac + 3*ido - 1);
        tr4 = ref(cc,ac + ido - 1) - ref(cc,ac + 2*ido - 1);
        tr3 = ref(cc,ac + ido - 1) + ref(cc,ac + 2*ido - 1);
        ah = (k - 1)*ido;
        ch[ah] = ref(cc,ac - ido - 1) + tr2 + tr3;
        ch[ah + 1] = ref(cc,ac - ido) + ti2 + ti3;
        cr2 = ref(cc,ac - ido - 1) + tr11*tr2 + tr12*tr3;
        ci2 = ref(cc,ac - ido) + tr11*ti2 + tr12*ti3;
        cr3 = ref(cc,ac - ido - 1) + tr12*tr2 + tr11*tr3;
        ci3 = ref(cc,ac - ido) + tr12*ti2 + tr11*ti3;
        cr5 = isign*(ti11*tr5 + ti12*tr4);
        ci5 = isign*(ti11*ti5 + ti12*ti4);
        cr4 = isign*(ti12*tr5 - ti11*tr4);
        ci4 = isign*(ti12*ti5 - ti11*ti4);
        ch[ah + l1*ido] = cr2 - ci5;
        ch[ah + 4*l1*ido] = cr2 + ci5;
        ch[ah + l1*ido + 1] = ci2 + cr5;
        ch[ah + 2*l1*ido + 1] = ci3 + cr4;
        ch[ah + 2*l1*ido] = cr3 - ci4;
        ch[ah + 3*l1*ido] = cr3 + ci4;
        ch[ah + 3*l1*ido + 1] = ci3 - cr4;
        ch[ah + 4*l1*ido + 1] = ci2 - cr5;
      }
    } else {
      for (k=1; k<=l1; k++) {
        for (i=0; i<ido-1; i+=2) {
          ac = i + 1 + (k*5 - 4)*ido;
          ti5 = ref(cc,ac) - ref(cc,ac + 3*ido);
          ti2 = ref(cc,ac) + ref(cc,ac + 3*ido);
          ti4 = ref(cc,ac + ido) - ref(cc,ac + 2*ido);
          ti3 = ref(cc,ac + ido) + ref(cc,ac + 2*ido);
          tr5 = ref(cc,ac - 1) - ref(cc,ac + 3*ido - 1);
          tr2 = ref(cc,ac - 1) + ref(cc,ac + 3*ido - 1);
          tr4 = ref(cc,ac + ido - 1) - ref(cc,ac + 2*ido - 1);
          tr3 = ref(cc,ac + ido - 1) + ref(cc,ac + 2*ido - 1);
          ah = i + (k - 1)*ido;
          ch[ah] = ref(cc,ac - ido - 1) + tr2 + tr3;
          ch[ah + 1] = ref(cc,ac - ido) + ti2 + ti3;
          cr2 = ref(cc,ac - ido - 1) + tr11*tr2 + tr12*tr3;

          ci2 = ref(cc,ac - ido) + tr11*ti2 + tr12*ti3;
          cr3 = ref(cc,ac - ido - 1) + tr12*tr2 + tr11*tr3;

          ci3 = ref(cc,ac - ido) + tr12*ti2 + tr11*ti3;
          cr5 = isign*(ti11*tr5 + ti12*tr4);
          ci5 = isign*(ti11*ti5 + ti12*ti4);
          cr4 = isign*(ti12*tr5 - ti11*tr4);
          ci4 = isign*(ti12*ti5 - ti11*ti4);
          dr3 = cr3 - ci4;
          dr4 = cr3 + ci4;
          di3 = ci3 + cr4;
          di4 = ci3 - cr4;
          dr5 = cr2 + ci5;
          dr2 = cr2 - ci5;
          di5 = ci2 - cr5;
          di2 = ci2 + cr5;
          ch[ah + l1*ido] = wa1[i]*dr2 - isign*wa1[i+1]*di2;
          ch[ah + l1*ido + 1] = wa1[i]*di2 + isign*wa1[i+1]*dr2;
          ch[ah + 2*l1*ido] = wa2[i]*dr3 - isign*wa2[i+1]*di3;
          ch[ah + 2*l1*ido + 1] = wa2[i]*di3 + isign*wa2[i+1]*dr3;
          ch[ah + 3*l1*ido] = wa3[i]*dr4 - isign*wa3[i+1]*di4;
          ch[ah + 3*l1*ido + 1] = wa3[i]*di4 + isign*wa3[i+1]*dr4;
          ch[ah + 4*l1*ido] = wa4[i]*dr5 - isign*wa4[i+1]*di5;
          ch[ah + 4*l1*ido + 1] = wa4[i]*di5 + isign*wa4[i+1]*dr5;
        }
      }
    }
  } /* passf5 */


static void passf(int *nac, int ido, int ip, int l1, int idl1,
      Treal cc[], Treal ch[],
      const Treal wa[], int isign)
  /* isign is -1 for forward transform and +1 for backward transform */
  {
    int idij, idlj, idot, ipph, i, j, k, l, jc, lc, ik, nt, idj, idl, inc,idp;
    Treal wai, war;

    idot = ido / 2;
    nt = ip*idl1;
    ipph = (ip + 1) / 2;
    idp = ip*ido;
    if (ido >= l1) {
      for (j=1; j<ipph; j++) {
        jc = ip - j;
        for (k=0; k<l1; k++) {
          for (i=0; i<ido; i++) {
            ch[i + (k + j*l1)*ido] =
                ref(cc,i + (j + k*ip)*ido) + ref(cc,i + (jc + k*ip)*ido);
            ch[i + (k + jc*l1)*ido] =
                ref(cc,i + (j + k*ip)*ido) - ref(cc,i + (jc + k*ip)*ido);
          }
        }
      }
      for (k=0; k<l1; k++)
        for (i=0; i<ido; i++)
          ch[i + k*ido] = ref(cc,i + k*ip*ido);
    } else {
      for (j=1; j<ipph; j++) {
        jc = ip - j;
        for (i=0; i<ido; i++) {
          for (k=0; k<l1; k++) {
            ch[i + (k + j*l1)*ido] = ref(cc,i + (j + k*ip)*ido) + ref(cc,i + (jc + k*
                ip)*ido);
            ch[i + (k + jc*l1)*ido] = ref(cc,i + (j + k*ip)*ido) - ref(cc,i + (jc + k*
                ip)*ido);
          }
        }
      }
      for (i=0; i<ido; i++)
        for (k=0; k<l1; k++)
          ch[i + k*ido] = ref(cc,i + k*ip*ido);
    }

    idl = 2 - ido;
    inc = 0;
    for (l=1; l<ipph; l++) {
      lc = ip - l;
      idl += ido;
      for (ik=0; ik<idl1; ik++) {
        cc[ik + l*idl1] = ch[ik] + wa[idl - 2]*ch[ik + idl1];
        cc[ik + lc*idl1] = isign*wa[idl-1]*ch[ik + (ip-1)*idl1];
      }
      idlj = idl;
      inc += ido;
      for (j=2; j<ipph; j++) {
        jc = ip - j;
        idlj += inc;
        if (idlj > idp) idlj -= idp;
        war = wa[idlj - 2];
        wai = wa[idlj-1];
        for (ik=0; ik<idl1; ik++) {
          cc[ik + l*idl1] += war*ch[ik + j*idl1];
          cc[ik + lc*idl1] += isign*wai*ch[ik + jc*idl1];
        }
      }
    }
    for (j=1; j<ipph; j++)
      for (ik=0; ik<idl1; ik++)
        ch[ik] += ch[ik + j*idl1];
    for (j=1; j<ipph; j++) {
      jc = ip - j;
      for (ik=1; ik<idl1; ik+=2) {
        ch[ik - 1 + j*idl1] = cc[ik - 1 + j*idl1] - cc[ik + jc*idl1];
        ch[ik - 1 + jc*idl1] = cc[ik - 1 + j*idl1] + cc[ik + jc*idl1];
        ch[ik + j*idl1] = cc[ik + j*idl1] + cc[ik - 1 + jc*idl1];
        ch[ik + jc*idl1] = cc[ik + j*idl1] - cc[ik - 1 + jc*idl1];
      }
    }
    *nac = 1;
    if (ido == 2) return;
    *nac = 0;
    for (ik=0; ik<idl1; ik++)
      cc[ik] = ch[ik];
    for (j=1; j<ip; j++) {
      for (k=0; k<l1; k++) {
        cc[(k + j*l1)*ido + 0] = ch[(k + j*l1)*ido + 0];
        cc[(k + j*l1)*ido + 1] = ch[(k + j*l1)*ido + 1];
      }
    }
    if (idot <= l1) {
      idij = 0;
      for (j=1; j<ip; j++) {
        idij += 2;
        for (i=3; i<ido; i+=2) {
          idij += 2;
          for (k=0; k<l1; k++) {
            cc[i - 1 + (k + j*l1)*ido] =
                wa[idij - 2]*ch[i - 1 + (k + j*l1)*ido] -
                isign*wa[idij-1]*ch[i + (k + j*l1)*ido];
            cc[i + (k + j*l1)*ido] =
                wa[idij - 2]*ch[i + (k + j*l1)*ido] +
                isign*wa[idij-1]*ch[i - 1 + (k + j*l1)*ido];
          }
        }
      }
    } else {
      idj = 2 - ido;
      for (j=1; j<ip; j++) {
        idj += ido;
        for (k = 0; k < l1; k++) {
          idij = idj;
          for (i=3; i<ido; i+=2) {
            idij += 2;
            cc[i - 1 + (k + j*l1)*ido] =
                wa[idij - 2]*ch[i - 1 + (k + j*l1)*ido] -
                isign*wa[idij-1]*ch[i + (k + j*l1)*ido];
            cc[i + (k + j*l1)*ido] =
                wa[idij - 2]*ch[i + (k + j*l1)*ido] +
                isign*wa[idij-1]*ch[i - 1 + (k + j*l1)*ido];
          }
        }
      }
    }
  } /* passf */


  /* ----------------------------------------------------------------------
radf2,radb2, radf3,radb3, radf4,radb4, radf5,radb5, radfg,radbg.
Treal FFT passes fwd and bwd.
---------------------------------------------------------------------- */

static void radf2(int ido, int l1, const Treal cc[], Treal ch[], const Treal wa1[])
  {
    int i, k, ic;
    Treal ti2, tr2;
    for (k=0; k<l1; k++) {
      ch[2*k*ido] =
          ref(cc,k*ido) + ref(cc,(k + l1)*ido);
      ch[(2*k+1)*ido + ido-1] =
          ref(cc,k*ido) - ref(cc,(k + l1)*ido);
    }
    if (ido < 2) return;
    if (ido != 2) {
      for (k=0; k<l1; k++) {
        for (i=2; i<ido; i+=2) {
          ic = ido - i;
          tr2 = wa1[i - 2]*ref(cc, i-1 + (k + l1)*ido) + wa1[i - 1]*ref(cc, i + (k + l1)*ido);
          ti2 = wa1[i - 2]*ref(cc, i + (k + l1)*ido) - wa1[i - 1]*ref(cc, i-1 + (k + l1)*ido);
          ch[i + 2*k*ido] = ref(cc,i + k*ido) + ti2;
          ch[ic + (2*k+1)*ido] = ti2 - ref(cc,i + k*ido);
          ch[i - 1 + 2*k*ido] = ref(cc,i - 1 + k*ido) + tr2;
          ch[ic - 1 + (2*k+1)*ido] = ref(cc,i - 1 + k*ido) - tr2;
        }
      }
      if (ido % 2 == 1) return;
    }
    for (k=0; k<l1; k++) {
      ch[(2*k+1)*ido] = -ref(cc,ido-1 + (k + l1)*ido);
      ch[ido-1 + 2*k*ido] = ref(cc,ido-1 + k*ido);
    }
  } /* radf2 */


static void radb2(int ido, int l1, const Treal cc[], Treal ch[], const Treal wa1[])
  {
    int i, k, ic;
    Treal ti2, tr2;
    for (k=0; k<l1; k++) {
      ch[k*ido] =
          ref(cc,2*k*ido) + ref(cc,ido-1 + (2*k+1)*ido);
      ch[(k + l1)*ido] =
          ref(cc,2*k*ido) - ref(cc,ido-1 + (2*k+1)*ido);
    }
    if (ido < 2) return;
    if (ido != 2) {
      for (k = 0; k < l1; ++k) {
        for (i = 2; i < ido; i += 2) {
          ic = ido - i;
          ch[i-1 + k*ido] =
              ref(cc,i-1 + 2*k*ido) + ref(cc,ic-1 + (2*k+1)*ido);
          tr2 = ref(cc,i-1 + 2*k*ido) - ref(cc,ic-1 + (2*k+1)*ido);
          ch[i + k*ido] =
              ref(cc,i + 2*k*ido) - ref(cc,ic + (2*k+1)*ido);
          ti2 = ref(cc,i + (2*k)*ido) + ref(cc,ic + (2*k+1)*ido);
          ch[i-1 + (k + l1)*ido] =
              wa1[i - 2]*tr2 - wa1[i - 1]*ti2;
          ch[i + (k + l1)*ido] =
              wa1[i - 2]*ti2 + wa1[i - 1]*tr2;
        }
      }
      if (ido % 2 == 1) return;
    }
    for (k = 0; k < l1; k++) {
      ch[ido-1 + k*ido] = 2*ref(cc,ido-1 + 2*k*ido);
      ch[ido-1 + (k + l1)*ido] = -2*ref(cc,(2*k+1)*ido);
    }
  } /* radb2 */


static void radf3(int ido, int l1, const Treal cc[], Treal ch[],
      const Treal wa1[], const Treal wa2[])
  {
    static const Treal taur = -0.5;
    static const Treal taui = 0.866025403784439;
    int i, k, ic;
    Treal ci2, di2, di3, cr2, dr2, dr3, ti2, ti3, tr2, tr3;
    for (k=0; k<l1; k++) {
      cr2 = ref(cc,(k + l1)*ido) + ref(cc,(k + 2*l1)*ido);
      ch[3*k*ido] = ref(cc,k*ido) + cr2;
      ch[(3*k+2)*ido] = taui*(ref(cc,(k + l1*2)*ido) - ref(cc,(k + l1)*ido));
      ch[ido-1 + (3*k + 1)*ido] = ref(cc,k*ido) + taur*cr2;
    }
    if (ido == 1) return;
    for (k=0; k<l1; k++) {
      for (i=2; i<ido; i+=2) {
        ic = ido - i;
        dr2 = wa1[i - 2]*ref(cc,i - 1 + (k + l1)*ido) +
            wa1[i - 1]*ref(cc,i + (k + l1)*ido);
        di2 = wa1[i - 2]*ref(cc,i + (k + l1)*ido) - wa1[i - 1]*ref(cc,i - 1 + (k + l1)*ido);
        dr3 = wa2[i - 2]*ref(cc,i - 1 + (k + l1*2)*ido) + wa2[i - 1]*ref(cc,i + (k + l1*2)*ido);
        di3 = wa2[i - 2]*ref(cc,i + (k + l1*2)*ido) - wa2[i - 1]*ref(cc,i - 1 + (k + l1*2)*ido);
        cr2 = dr2 + dr3;
        ci2 = di2 + di3;
        ch[i - 1 + 3*k*ido] = ref(cc,i - 1 + k*ido) + cr2;
        ch[i + 3*k*ido] = ref(cc,i + k*ido) + ci2;
        tr2 = ref(cc,i - 1 + k*ido) + taur*cr2;
        ti2 = ref(cc,i + k*ido) + taur*ci2;
        tr3 = taui*(di2 - di3);
        ti3 = taui*(dr3 - dr2);
        ch[i - 1 + (3*k + 2)*ido] = tr2 + tr3;
        ch[ic - 1 + (3*k + 1)*ido] = tr2 - tr3;
        ch[i + (3*k + 2)*ido] = ti2 + ti3;
        ch[ic + (3*k + 1)*ido] = ti3 - ti2;
      }
    }
  } /* radf3 */


static void radb3(int ido, int l1, const Treal cc[], Treal ch[],
      const Treal wa1[], const Treal wa2[])
  {
    static const Treal taur = -0.5;
    static const Treal taui = 0.866025403784439;
    int i, k, ic;
    Treal ci2, ci3, di2, di3, cr2, cr3, dr2, dr3, ti2, tr2;
    for (k=0; k<l1; k++) {
      tr2 = 2*ref(cc,ido-1 + (3*k + 1)*ido);
      cr2 = ref(cc,3*k*ido) + taur*tr2;
      ch[k*ido] = ref(cc,3*k*ido) + tr2;
      ci3 = 2*taui*ref(cc,(3*k + 2)*ido);
      ch[(k + l1)*ido] = cr2 - ci3;
      ch[(k + 2*l1)*ido] = cr2 + ci3;
    }
    if (ido == 1) return;
    for (k=0; k<l1; k++) {
      for (i=2; i<ido; i+=2) {
        ic = ido - i;
        tr2 = ref(cc,i - 1 + (3*k + 2)*ido) + ref(cc,ic - 1 + (3*k + 1)*ido);
        cr2 = ref(cc,i - 1 + 3*k*ido) + taur*tr2;
        ch[i - 1 + k*ido] = ref(cc,i - 1 + 3*k*ido) + tr2;
        ti2 = ref(cc,i + (3*k + 2)*ido) - ref(cc,ic + (3*k + 1)*ido);
        ci2 = ref(cc,i + 3*k*ido) + taur*ti2;
        ch[i + k*ido] = ref(cc,i + 3*k*ido) + ti2;
        cr3 = taui*(ref(cc,i - 1 + (3*k + 2)*ido) - ref(cc,ic - 1 + (3*k + 1)*ido));
        ci3 = taui*(ref(cc,i + (3*k + 2)*ido) + ref(cc,ic + (3*k + 1)*ido));
        dr2 = cr2 - ci3;
        dr3 = cr2 + ci3;
        di2 = ci2 + cr3;
        di3 = ci2 - cr3;
        ch[i - 1 + (k + l1)*ido] = wa1[i - 2]*dr2 - wa1[i - 1]*di2;
        ch[i + (k + l1)*ido] = wa1[i - 2]*di2 + wa1[i - 1]*dr2;
        ch[i - 1 + (k + 2*l1)*ido] = wa2[i - 2]*dr3 - wa2[i - 1]*di3;
        ch[i + (k + 2*l1)*ido] = wa2[i - 2]*di3 + wa2[i - 1]*dr3;
      }
    }
  } /* radb3 */


static void radf4(int ido, int l1, const Treal cc[], Treal ch[],
      const Treal wa1[], const Treal wa2[], const Treal wa3[])
  {
    static const Treal hsqt2 = 0.7071067811865475;
    int i, k, ic;
    Treal ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
    for (k=0; k<l1; k++) {
      tr1 = ref(cc,(k + l1)*ido) + ref(cc,(k + 3*l1)*ido);
      tr2 = ref(cc,k*ido) + ref(cc,(k + 2*l1)*ido);
      ch[4*k*ido] = tr1 + tr2;
      ch[ido-1 + (4*k + 3)*ido] = tr2 - tr1;
      ch[ido-1 + (4*k + 1)*ido] = ref(cc,k*ido) - ref(cc,(k + 2*l1)*ido);
      ch[(4*k + 2)*ido] = ref(cc,(k + 3*l1)*ido) - ref(cc,(k + l1)*ido);
    }
    if (ido < 2) return;
    if (ido != 2) {
      for (k=0; k<l1; k++) {
        for (i=2; i<ido; i += 2) {
          ic = ido - i;
          cr2 = wa1[i - 2]*ref(cc,i - 1 + (k + l1)*ido) + wa1[i - 1]*ref(cc,i + (k + l1)*ido);
          ci2 = wa1[i - 2]*ref(cc,i + (k + l1)*ido) - wa1[i - 1]*ref(cc,i - 1 + (k + l1)*ido);
          cr3 = wa2[i - 2]*ref(cc,i - 1 + (k + 2*l1)*ido) + wa2[i - 1]*ref(cc,i + (k + 2*l1)*
              ido);
          ci3 = wa2[i - 2]*ref(cc,i + (k + 2*l1)*ido) - wa2[i - 1]*ref(cc,i - 1 + (k + 2*l1)*
              ido);
          cr4 = wa3[i - 2]*ref(cc,i - 1 + (k + 3*l1)*ido) + wa3[i - 1]*ref(cc,i + (k + 3*l1)*
              ido);
          ci4 = wa3[i - 2]*ref(cc,i + (k + 3*l1)*ido) - wa3[i - 1]*ref(cc,i - 1 + (k + 3*l1)*
              ido);
          tr1 = cr2 + cr4;
          tr4 = cr4 - cr2;
          ti1 = ci2 + ci4;
          ti4 = ci2 - ci4;
          ti2 = ref(cc,i + k*ido) + ci3;
          ti3 = ref(cc,i + k*ido) - ci3;
          tr2 = ref(cc,i - 1 + k*ido) + cr3;
          tr3 = ref(cc,i - 1 + k*ido) - cr3;
          ch[i - 1 + 4*k*ido] = tr1 + tr2;
          ch[ic - 1 + (4*k + 3)*ido] = tr2 - tr1;
          ch[i + 4*k*ido] = ti1 + ti2;
          ch[ic + (4*k + 3)*ido] = ti1 - ti2;
          ch[i - 1 + (4*k + 2)*ido] = ti4 + tr3;
          ch[ic - 1 + (4*k + 1)*ido] = tr3 - ti4;
          ch[i + (4*k + 2)*ido] = tr4 + ti3;
          ch[ic + (4*k + 1)*ido] = tr4 - ti3;
        }
      }
      if (ido % 2 == 1) return;
    }
    for (k=0; k<l1; k++) {
      ti1 = -hsqt2*(ref(cc,ido-1 + (k + l1)*ido) + ref(cc,ido-1 + (k + 3*l1)*ido));
      tr1 = hsqt2*(ref(cc,ido-1 + (k + l1)*ido) - ref(cc,ido-1 + (k + 3*l1)*ido));
      ch[ido-1 + 4*k*ido] = tr1 + ref(cc,ido-1 + k*ido);
      ch[ido-1 + (4*k + 2)*ido] = ref(cc,ido-1 + k*ido) - tr1;
      ch[(4*k + 1)*ido] = ti1 - ref(cc,ido-1 + (k + 2*l1)*ido);
      ch[(4*k + 3)*ido] = ti1 + ref(cc,ido-1 + (k + 2*l1)*ido);
    }
  } /* radf4 */


static void radb4(int ido, int l1, const Treal cc[], Treal ch[],
      const Treal wa1[], const Treal wa2[], const Treal wa3[])
  {
    static const Treal sqrt2 = 1.414213562373095;
    int i, k, ic;
    Treal ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
    for (k = 0; k < l1; k++) {
      tr1 = ref(cc,4*k*ido) - ref(cc,ido-1 + (4*k + 3)*ido);
      tr2 = ref(cc,4*k*ido) + ref(cc,ido-1 + (4*k + 3)*ido);
      tr3 = ref(cc,ido-1 + (4*k + 1)*ido) + ref(cc,ido-1 + (4*k + 1)*ido);
      tr4 = ref(cc,(4*k + 2)*ido) + ref(cc,(4*k + 2)*ido);
      ch[k*ido] = tr2 + tr3;
      ch[(k + l1)*ido] = tr1 - tr4;
      ch[(k + 2*l1)*ido] = tr2 - tr3;
      ch[(k + 3*l1)*ido] = tr1 + tr4;
    }
    if (ido < 2) return;
    if (ido != 2) {
      for (k = 0; k < l1; ++k) {
        for (i = 2; i < ido; i += 2) {
          ic = ido - i;
          ti1 = ref(cc,i + 4*k*ido) + ref(cc,ic + (4*k + 3)*ido);
          ti2 = ref(cc,i + 4*k*ido) - ref(cc,ic + (4*k + 3)*ido);
          ti3 = ref(cc,i + (4*k + 2)*ido) - ref(cc,ic + (4*k + 1)*ido);
          tr4 = ref(cc,i + (4*k + 2)*ido) + ref(cc,ic + (4*k + 1)*ido);
          tr1 = ref(cc,i - 1 + 4*k*ido) - ref(cc,ic - 1 + (4*k + 3)*ido);
          tr2 = ref(cc,i - 1 + 4*k*ido) + ref(cc,ic - 1 + (4*k + 3)*ido);
          ti4 = ref(cc,i - 1 + (4*k + 2)*ido) - ref(cc,ic - 1 + (4*k + 1)*ido);
          tr3 = ref(cc,i - 1 + (4*k + 2)*ido) + ref(cc,ic - 1 + (4*k + 1)*ido);
          ch[i - 1 + k*ido] = tr2 + tr3;
          cr3 = tr2 - tr3;
          ch[i + k*ido] = ti2 + ti3;
          ci3 = ti2 - ti3;
          cr2 = tr1 - tr4;
          cr4 = tr1 + tr4;
          ci2 = ti1 + ti4;
          ci4 = ti1 - ti4;
          ch[i - 1 + (k + l1)*ido] = wa1[i - 2]*cr2 - wa1[i - 1]*ci2;
          ch[i + (k + l1)*ido] = wa1[i - 2]*ci2 + wa1[i - 1]*cr2;
          ch[i - 1 + (k + 2*l1)*ido] = wa2[i - 2]*cr3 - wa2[i - 1]*ci3;
          ch[i + (k + 2*l1)*ido] = wa2[i - 2]*ci3 + wa2[i - 1]*cr3;
          ch[i - 1 + (k + 3*l1)*ido] = wa3[i - 2]*cr4 - wa3[i - 1]*ci4;
          ch[i + (k + 3*l1)*ido] = wa3[i - 2]*ci4 + wa3[i - 1]*cr4;
        }
      }
      if (ido % 2 == 1) return;
    }
    for (k = 0; k < l1; k++) {
      ti1 = ref(cc,(4*k + 1)*ido) + ref(cc,(4*k + 3)*ido);
      ti2 = ref(cc,(4*k + 3)*ido) - ref(cc,(4*k + 1)*ido);
      tr1 = ref(cc,ido-1 + 4*k*ido) - ref(cc,ido-1 + (4*k + 2)*ido);
      tr2 = ref(cc,ido-1 + 4*k*ido) + ref(cc,ido-1 + (4*k + 2)*ido);
      ch[ido-1 + k*ido] = tr2 + tr2;
      ch[ido-1 + (k + l1)*ido] = sqrt2*(tr1 - ti1);
      ch[ido-1 + (k + 2*l1)*ido] = ti2 + ti2;
      ch[ido-1 + (k + 3*l1)*ido] = -sqrt2*(tr1 + ti1);
    }
  } /* radb4 */


static void radf5(int ido, int l1, const Treal cc[], Treal ch[],
      const Treal wa1[], const Treal wa2[], const Treal wa3[], const Treal wa4[])
  {
    static const Treal tr11 = 0.309016994374947;
    static const Treal ti11 = 0.951056516295154;
    static const Treal tr12 = -0.809016994374947;
    static const Treal ti12 = 0.587785252292473;
    int i, k, ic;
    Treal ci2, di2, ci4, ci5, di3, di4, di5, ci3, cr2, cr3, dr2, dr3, dr4, dr5,
        cr5, cr4, ti2, ti3, ti5, ti4, tr2, tr3, tr4, tr5;
    for (k = 0; k < l1; k++) {
      cr2 = ref(cc,(k + 4*l1)*ido) + ref(cc,(k + l1)*ido);
      ci5 = ref(cc,(k + 4*l1)*ido) - ref(cc,(k + l1)*ido);
      cr3 = ref(cc,(k + 3*l1)*ido) + ref(cc,(k + 2*l1)*ido);
      ci4 = ref(cc,(k + 3*l1)*ido) - ref(cc,(k + 2*l1)*ido);
      ch[5*k*ido] = ref(cc,k*ido) + cr2 + cr3;
      ch[ido-1 + (5*k + 1)*ido] = ref(cc,k*ido) + tr11*cr2 + tr12*cr3;
      ch[(5*k + 2)*ido] = ti11*ci5 + ti12*ci4;
      ch[ido-1 + (5*k + 3)*ido] = ref(cc,k*ido) + tr12*cr2 + tr11*cr3;
      ch[(5*k + 4)*ido] = ti12*ci5 - ti11*ci4;
    }
    if (ido == 1) return;
    for (k = 0; k < l1; ++k) {
      for (i = 2; i < ido; i += 2) {
        ic = ido - i;
        dr2 = wa1[i - 2]*ref(cc,i - 1 + (k + l1)*ido) + wa1[i - 1]*ref(cc,i + (k + l1)*ido);
        di2 = wa1[i - 2]*ref(cc,i + (k + l1)*ido) - wa1[i - 1]*ref(cc,i - 1 + (k + l1)*ido);
        dr3 = wa2[i - 2]*ref(cc,i - 1 + (k + 2*l1)*ido) + wa2[i - 1]*ref(cc,i + (k + 2*l1)*ido);
        di3 = wa2[i - 2]*ref(cc,i + (k + 2*l1)*ido) - wa2[i - 1]*ref(cc,i - 1 + (k + 2*l1)*ido);
        dr4 = wa3[i - 2]*ref(cc,i - 1 + (k + 3*l1)*ido) + wa3[i - 1]*ref(cc,i + (k + 3*l1)*ido);
        di4 = wa3[i - 2]*ref(cc,i + (k + 3*l1)*ido) - wa3[i - 1]*ref(cc,i - 1 + (k + 3*l1)*ido);
        dr5 = wa4[i - 2]*ref(cc,i - 1 + (k + 4*l1)*ido) + wa4[i - 1]*ref(cc,i + (k + 4*l1)*ido);
        di5 = wa4[i - 2]*ref(cc,i + (k + 4*l1)*ido) - wa4[i - 1]*ref(cc,i - 1 + (k + 4*l1)*ido);
        cr2 = dr2 + dr5;
        ci5 = dr5 - dr2;
        cr5 = di2 - di5;
        ci2 = di2 + di5;
        cr3 = dr3 + dr4;
        ci4 = dr4 - dr3;
        cr4 = di3 - di4;
        ci3 = di3 + di4;
        ch[i - 1 + 5*k*ido] = ref(cc,i - 1 + k*ido) + cr2 + cr3;
        ch[i + 5*k*ido] = ref(cc,i + k*ido) + ci2 + ci3;
        tr2 = ref(cc,i - 1 + k*ido) + tr11*cr2 + tr12*cr3;
        ti2 = ref(cc,i + k*ido) + tr11*ci2 + tr12*ci3;
        tr3 = ref(cc,i - 1 + k*ido) + tr12*cr2 + tr11*cr3;
        ti3 = ref(cc,i + k*ido) + tr12*ci2 + tr11*ci3;
        tr5 = ti11*cr5 + ti12*cr4;
        ti5 = ti11*ci5 + ti12*ci4;
        tr4 = ti12*cr5 - ti11*cr4;
        ti4 = ti12*ci5 - ti11*ci4;
        ch[i - 1 + (5*k + 2)*ido] = tr2 + tr5;
        ch[ic - 1 + (5*k + 1)*ido] = tr2 - tr5;
        ch[i + (5*k + 2)*ido] = ti2 + ti5;
        ch[ic + (5*k + 1)*ido] = ti5 - ti2;
        ch[i - 1 + (5*k + 4)*ido] = tr3 + tr4;
        ch[ic - 1 + (5*k + 3)*ido] = tr3 - tr4;
        ch[i + (5*k + 4)*ido] = ti3 + ti4;
        ch[ic + (5*k + 3)*ido] = ti4 - ti3;
      }
    }
  } /* radf5 */


static void radb5(int ido, int l1, const Treal cc[], Treal ch[],
      const Treal wa1[], const Treal wa2[], const Treal wa3[], const Treal wa4[])
  {
    static const Treal tr11 = 0.309016994374947;
    static const Treal ti11 = 0.951056516295154;
    static const Treal tr12 = -0.809016994374947;
    static const Treal ti12 = 0.587785252292473;
    int i, k, ic;
    Treal ci2, ci3, ci4, ci5, di3, di4, di5, di2, cr2, cr3, cr5, cr4, ti2, ti3,
        ti4, ti5, dr3, dr4, dr5, dr2, tr2, tr3, tr4, tr5;
    for (k = 0; k < l1; k++) {
      ti5 = 2*ref(cc,(5*k + 2)*ido);
      ti4 = 2*ref(cc,(5*k + 4)*ido);
      tr2 = 2*ref(cc,ido-1 + (5*k + 1)*ido);
      tr3 = 2*ref(cc,ido-1 + (5*k + 3)*ido);
      ch[k*ido] = ref(cc,5*k*ido) + tr2 + tr3;
      cr2 = ref(cc,5*k*ido) + tr11*tr2 + tr12*tr3;
      cr3 = ref(cc,5*k*ido) + tr12*tr2 + tr11*tr3;
      ci5 = ti11*ti5 + ti12*ti4;
      ci4 = ti12*ti5 - ti11*ti4;
      ch[(k + l1)*ido] = cr2 - ci5;
      ch[(k + 2*l1)*ido] = cr3 - ci4;
      ch[(k + 3*l1)*ido] = cr3 + ci4;
      ch[(k + 4*l1)*ido] = cr2 + ci5;
    }
    if (ido == 1) return;
    for (k = 0; k < l1; ++k) {
      for (i = 2; i < ido; i += 2) {
        ic = ido - i;
        ti5 = ref(cc,i + (5*k + 2)*ido) + ref(cc,ic + (5*k + 1)*ido);
        ti2 = ref(cc,i + (5*k + 2)*ido) - ref(cc,ic + (5*k + 1)*ido);
        ti4 = ref(cc,i + (5*k + 4)*ido) + ref(cc,ic + (5*k + 3)*ido);
        ti3 = ref(cc,i + (5*k + 4)*ido) - ref(cc,ic + (5*k + 3)*ido);
        tr5 = ref(cc,i - 1 + (5*k + 2)*ido) - ref(cc,ic - 1 + (5*k + 1)*ido);
        tr2 = ref(cc,i - 1 + (5*k + 2)*ido) + ref(cc,ic - 1 + (5*k + 1)*ido);
        tr4 = ref(cc,i - 1 + (5*k + 4)*ido) - ref(cc,ic - 1 + (5*k + 3)*ido);
        tr3 = ref(cc,i - 1 + (5*k + 4)*ido) + ref(cc,ic - 1 + (5*k + 3)*ido);
        ch[i - 1 + k*ido] = ref(cc,i - 1 + 5*k*ido) + tr2 + tr3;
        ch[i + k*ido] = ref(cc,i + 5*k*ido) + ti2 + ti3;
        cr2 = ref(cc,i - 1 + 5*k*ido) + tr11*tr2 + tr12*tr3;

        ci2 = ref(cc,i + 5*k*ido) + tr11*ti2 + tr12*ti3;
        cr3 = ref(cc,i - 1 + 5*k*ido) + tr12*tr2 + tr11*tr3;

        ci3 = ref(cc,i + 5*k*ido) + tr12*ti2 + tr11*ti3;
        cr5 = ti11*tr5 + ti12*tr4;
        ci5 = ti11*ti5 + ti12*ti4;
        cr4 = ti12*tr5 - ti11*tr4;
        ci4 = ti12*ti5 - ti11*ti4;
        dr3 = cr3 - ci4;
        dr4 = cr3 + ci4;
        di3 = ci3 + cr4;
        di4 = ci3 - cr4;
        dr5 = cr2 + ci5;
        dr2 = cr2 - ci5;
        di5 = ci2 - cr5;
        di2 = ci2 + cr5;
        ch[i - 1 + (k + l1)*ido] = wa1[i - 2]*dr2 - wa1[i - 1]*di2;
        ch[i + (k + l1)*ido] = wa1[i - 2]*di2 + wa1[i - 1]*dr2;
        ch[i - 1 + (k + 2*l1)*ido] = wa2[i - 2]*dr3 - wa2[i - 1]*di3;
        ch[i + (k + 2*l1)*ido] = wa2[i - 2]*di3 + wa2[i - 1]*dr3;
        ch[i - 1 + (k + 3*l1)*ido] = wa3[i - 2]*dr4 - wa3[i - 1]*di4;
        ch[i + (k + 3*l1)*ido] = wa3[i - 2]*di4 + wa3[i - 1]*dr4;
        ch[i - 1 + (k + 4*l1)*ido] = wa4[i - 2]*dr5 - wa4[i - 1]*di5;
        ch[i + (k + 4*l1)*ido] = wa4[i - 2]*di5 + wa4[i - 1]*dr5;
      }
    }
  } /* radb5 */


static void radfg(int ido, int ip, int l1, int idl1,
      Treal cc[], Treal ch[], const Treal wa[])
  {
    static const Treal twopi = 6.28318530717959;
    int idij, ipph, i, j, k, l, j2, ic, jc, lc, ik, is, nbd;
    Treal dc2, ai1, ai2, ar1, ar2, ds2, dcp, arg, dsp, ar1h, ar2h;
    arg = twopi / ip;
    dcp = cos(arg);
    dsp = sin(arg);
    ipph = (ip + 1) / 2;
    nbd = (ido - 1) / 2;
    if (ido != 1) {
      for (ik=0; ik<idl1; ik++) ch[ik] = cc[ik];
      for (j=1; j<ip; j++)
        for (k=0; k<l1; k++)
          ch[(k + j*l1)*ido] = cc[(k + j*l1)*ido];
      if (nbd <= l1) {
        is = -ido;
        for (j=1; j<ip; j++) {
          is += ido;
          idij = is-1;
          for (i=2; i<ido; i+=2) {
            idij += 2;
            for (k=0; k<l1; k++) {
              ch[i - 1 + (k + j*l1)*ido] =
                  wa[idij - 1]*cc[i - 1 + (k + j*l1)*ido] + wa[idij]*cc[i + (k + j*l1)*ido];
              ch[i + (k + j*l1)*ido] =
                  wa[idij - 1]*cc[i + (k + j*l1)*ido] - wa[idij]*cc[i - 1 + (k + j*l1)*ido];
            }
          }
        }
      } else {
        is = -ido;
        for (j=1; j<ip; j++) {
          is += ido;
          for (k=0; k<l1; k++) {
            idij = is-1;
            for (i=2; i<ido; i+=2) {
              idij += 2;
              ch[i - 1 + (k + j*l1)*ido] =
                  wa[idij - 1]*cc[i - 1 + (k + j*l1)*ido] + wa[idij]*cc[i + (k + j*l1)*ido];
              ch[i + (k + j*l1)*ido] =
                  wa[idij - 1]*cc[i + (k + j*l1)*ido] - wa[idij]*cc[i - 1 + (k + j*l1)*ido];
            }
          }
        }
      }
      if (nbd >= l1) {
        for (j=1; j<ipph; j++) {
          jc = ip - j;
          for (k=0; k<l1; k++) {
            for (i=2; i<ido; i+=2) {
              cc[i - 1 + (k + j*l1)*ido] = ch[i - 1 + (k + j*l1)*ido] + ch[i - 1 + (k + jc*l1)*ido];
              cc[i - 1 + (k + jc*l1)*ido] = ch[i + (k + j*l1)*ido] - ch[i + (k + jc*l1)*ido];
              cc[i + (k + j*l1)*ido] = ch[i + (k + j*l1)*ido] + ch[i + (k + jc*l1)*ido];
              cc[i + (k + jc*l1)*ido] = ch[i - 1 + (k + jc*l1)*ido] - ch[i - 1 + (k + j*l1)*ido];
            }
          }
        }
      } else {
        for (j=1; j<ipph; j++) {
          jc = ip - j;
          for (i=2; i<ido; i+=2) {
            for (k=0; k<l1; k++) {
              cc[i - 1 + (k + j*l1)*ido] =
                  ch[i - 1 + (k + j*l1)*ido] + ch[i - 1 + (k + jc*l1)*ido];
              cc[i - 1 + (k + jc*l1)*ido] = ch[i + (k + j*l1)*ido] - ch[i + (k + jc*l1)*ido];
              cc[i + (k + j*l1)*ido] = ch[i + (k + j*l1)*ido] + ch[i + (k + jc*l1)*ido];
              cc[i + (k + jc*l1)*ido] = ch[i - 1 + (k + jc*l1)*ido] - ch[i - 1 + (k + j*l1)*ido];
            }
          }
        }
      }
    } else {  /* now ido == 1 */
      for (ik=0; ik<idl1; ik++) cc[ik] = ch[ik];
    }
    for (j=1; j<ipph; j++) {
      jc = ip - j;
      for (k=0; k<l1; k++) {
        cc[(k + j*l1)*ido] = ch[(k + j*l1)*ido] + ch[(k + jc*l1)*ido];
        cc[(k + jc*l1)*ido] = ch[(k + jc*l1)*ido] - ch[(k + j*l1)*ido];
      }
    }

    ar1 = 1;
    ai1 = 0;
    for (l=1; l<ipph; l++) {
      lc = ip - l;
      ar1h = dcp*ar1 - dsp*ai1;
      ai1 = dcp*ai1 + dsp*ar1;
      ar1 = ar1h;
      for (ik=0; ik<idl1; ik++) {
        ch[ik + l*idl1] = cc[ik] + ar1*cc[ik + idl1];
        ch[ik + lc*idl1] = ai1*cc[ik + (ip-1)*idl1];
      }
      dc2 = ar1;
      ds2 = ai1;
      ar2 = ar1;
      ai2 = ai1;
      for (j=2; j<ipph; j++) {
        jc = ip - j;
        ar2h = dc2*ar2 - ds2*ai2;
        ai2 = dc2*ai2 + ds2*ar2;
        ar2 = ar2h;
        for (ik=0; ik<idl1; ik++) {
          ch[ik + l*idl1] += ar2*cc[ik + j*idl1];
          ch[ik + lc*idl1] += ai2*cc[ik + jc*idl1];
        }
      }
    }
    for (j=1; j<ipph; j++)
      for (ik=0; ik<idl1; ik++)
        ch[ik] += cc[ik + j*idl1];

    if (ido >= l1) {
      for (k=0; k<l1; k++) {
        for (i=0; i<ido; i++) {
          ref(cc,i + k*ip*ido) = ch[i + k*ido];
        }
      }
    } else {
      for (i=0; i<ido; i++) {
        for (k=0; k<l1; k++) {
          ref(cc,i + k*ip*ido) = ch[i + k*ido];
        }
      }
    }
    for (j=1; j<ipph; j++) {
      jc = ip - j;
      j2 = 2*j;
      for (k=0; k<l1; k++) {
        ref(cc,ido-1 + (j2 - 1 + k*ip)*ido) =
            ch[(k + j*l1)*ido];
        ref(cc,(j2 + k*ip)*ido) =
            ch[(k + jc*l1)*ido];
      }
    }
    if (ido == 1) return;
    if (nbd >= l1) {
      for (j=1; j<ipph; j++) {
        jc = ip - j;
        j2 = 2*j;
        for (k=0; k<l1; k++) {
          for (i=2; i<ido; i+=2) {
            ic = ido - i;
            ref(cc,i - 1 + (j2 + k*ip)*ido) = ch[i - 1 + (k + j*l1)*ido] + ch[i - 1 + (k + jc*l1)*ido];
            ref(cc,ic - 1 + (j2 - 1 + k*ip)*ido) = ch[i - 1 + (k + j*l1)*ido] - ch[i - 1 + (k + jc*l1)*ido];
            ref(cc,i + (j2 + k*ip)*ido) = ch[i + (k + j*l1)*ido] + ch[i + (k + jc*l1)*ido];
            ref(cc,ic + (j2 - 1 + k*ip)*ido) = ch[i + (k + jc*l1)*ido] - ch[i + (k + j*l1)*ido];
          }
        }
      }
    } else {
      for (j=1; j<ipph; j++) {
        jc = ip - j;
        j2 = 2*j;
        for (i=2; i<ido; i+=2) {
          ic = ido - i;
          for (k=0; k<l1; k++) {
            ref(cc,i - 1 + (j2 + k*ip)*ido) = ch[i - 1 + (k + j*l1)*ido] + ch[i - 1 + (k + jc*l1)*ido];
            ref(cc,ic - 1 + (j2 - 1 + k*ip)*ido) = ch[i - 1 + (k + j*l1)*ido] - ch[i - 1 + (k + jc*l1)*ido];
            ref(cc,i + (j2 + k*ip)*ido) = ch[i + (k + j*l1)*ido] + ch[i + (k + jc*l1)*ido];
            ref(cc,ic + (j2 - 1 + k*ip)*ido) = ch[i + (k + jc*l1)*ido] - ch[i + (k + j*l1)*ido];
          }
        }
      }
    }
  } /* radfg */


static void radbg(int ido, int ip, int l1, int idl1,
      Treal cc[], Treal ch[], const Treal wa[])
  {
    static const Treal twopi = 6.28318530717959;
    int idij, ipph, i, j, k, l, j2, ic, jc, lc, ik, is;
    Treal dc2, ai1, ai2, ar1, ar2, ds2;
    int nbd;
    Treal dcp, arg, dsp, ar1h, ar2h;
    arg = twopi / ip;
    dcp = cos(arg);
    dsp = sin(arg);
    nbd = (ido - 1) / 2;
    ipph = (ip + 1) / 2;
    if (ido >= l1) {
      for (k=0; k<l1; k++) {
        for (i=0; i<ido; i++) {
          ch[i + k*ido] = ref(cc,i + k*ip*ido);
        }
      }
    } else {
      for (i=0; i<ido; i++) {
        for (k=0; k<l1; k++) {
          ch[i + k*ido] = ref(cc,i + k*ip*ido);
        }
      }
    }
    for (j=1; j<ipph; j++) {
      jc = ip - j;
      j2 = 2*j;
      for (k=0; k<l1; k++) {
        ch[(k + j*l1)*ido] = ref(cc,ido-1 + (j2 - 1 + k*ip)*ido) + ref(cc,ido-1 + (j2 - 1 + k*ip)*
            ido);
        ch[(k + jc*l1)*ido] = ref(cc,(j2 + k*ip)*ido) + ref(cc,(j2 + k*ip)*ido);
      }
    }

    if (ido != 1) {
      if (nbd >= l1) {
        for (j=1; j<ipph; j++) {
          jc = ip - j;
          for (k=0; k<l1; k++) {
            for (i=2; i<ido; i+=2) {
              ic = ido - i;
              ch[i - 1 + (k + j*l1)*ido] = ref(cc,i - 1 + (2*j + k*ip)*ido) + ref(cc,
                  ic - 1 + (2*j - 1 + k*ip)*ido);
              ch[i - 1 + (k + jc*l1)*ido] = ref(cc,i - 1 + (2*j + k*ip)*ido) -
                  ref(cc,ic - 1 + (2*j - 1 + k*ip)*ido);
              ch[i + (k + j*l1)*ido] = ref(cc,i + (2*j + k*ip)*ido) - ref(cc,ic
                  + (2*j - 1 + k*ip)*ido);
              ch[i + (k + jc*l1)*ido] = ref(cc,i + (2*j + k*ip)*ido) + ref(cc,ic
                  + (2*j - 1 + k*ip)*ido);
            }
          }
        }
      } else {
        for (j=1; j<ipph; j++) {
          jc = ip - j;
          for (i=2; i<ido; i+=2) {
            ic = ido - i;
            for (k=0; k<l1; k++) {
              ch[i - 1 + (k + j*l1)*ido] = ref(cc,i - 1 + (2*j + k*ip)*ido) + ref(cc,
                  ic - 1 + (2*j - 1 + k*ip)*ido);
              ch[i - 1 + (k + jc*l1)*ido] = ref(cc,i - 1 + (2*j + k*ip)*ido) -
                  ref(cc,ic - 1 + (2*j - 1 + k*ip)*ido);
              ch[i + (k + j*l1)*ido] = ref(cc,i + (2*j + k*ip)*ido) - ref(cc,ic
                  + (2*j - 1 + k*ip)*ido);
              ch[i + (k + jc*l1)*ido] = ref(cc,i + (2*j + k*ip)*ido) + ref(cc,ic
                  + (2*j - 1 + k*ip)*ido);
            }
          }
        }
      }
    }

    ar1 = 1;
    ai1 = 0;
    for (l=1; l<ipph; l++) {
      lc = ip - l;
      ar1h = dcp*ar1 - dsp*ai1;
      ai1 = dcp*ai1 + dsp*ar1;
      ar1 = ar1h;
      for (ik=0; ik<idl1; ik++) {
        cc[ik + l*idl1] = ch[ik] + ar1*ch[ik + idl1];
        cc[ik + lc*idl1] = ai1*ch[ik + (ip-1)*idl1];
      }
      dc2 = ar1;
      ds2 = ai1;
      ar2 = ar1;
      ai2 = ai1;
      for (j=2; j<ipph; j++) {
        jc = ip - j;
        ar2h = dc2*ar2 - ds2*ai2;
        ai2 = dc2*ai2 + ds2*ar2;
        ar2 = ar2h;
        for (ik=0; ik<idl1; ik++) {
          cc[ik + l*idl1] += ar2*ch[ik + j*idl1];
          cc[ik + lc*idl1] += ai2*ch[ik + jc*idl1];
        }
      }
    }
    for (j=1; j<ipph; j++) {
      for (ik=0; ik<idl1; ik++) {
        ch[ik] += ch[ik + j*idl1];
      }
    }
    for (j=1; j<ipph; j++) {
      jc = ip - j;
      for (k=0; k<l1; k++) {
        ch[(k + j*l1)*ido] = cc[(k + j*l1)*ido] - cc[(k + jc*l1)*ido];
        ch[(k + jc*l1)*ido] = cc[(k + j*l1)*ido] + cc[(k + jc*l1)*ido];
      }
    }

    if (ido == 1) return;
    if (nbd >= l1) {
      for (j=1; j<ipph; j++) {
        jc = ip - j;
        for (k=0; k<l1; k++) {
          for (i=2; i<ido; i+=2) {
            ch[i - 1 + (k + j*l1)*ido] = cc[i - 1 + (k + j*l1)*ido] - cc[i + (k + jc*l1)*ido];
            ch[i - 1 + (k + jc*l1)*ido] = cc[i - 1 + (k + j*l1)*ido] + cc[i + (k + jc*l1)*ido];
            ch[i + (k + j*l1)*ido] = cc[i + (k + j*l1)*ido] + cc[i - 1 + (k + jc*l1)*ido];
            ch[i + (k + jc*l1)*ido] = cc[i + (k + j*l1)*ido] - cc[i - 1 + (k + jc*l1)*ido];
          }
        }
      }
    } else {
      for (j=1; j<ipph; j++) {
        jc = ip - j;
        for (i=2; i<ido; i+=2) {
          for (k=0; k<l1; k++) {
            ch[i - 1 + (k + j*l1)*ido] = cc[i - 1 + (k + j*l1)*ido] - cc[i + (k + jc*l1)*ido];
            ch[i - 1 + (k + jc*l1)*ido] = cc[i - 1 + (k + j *l1)*ido] + cc[i + (k + jc*l1)*ido];
            ch[i + (k + j*l1)*ido] = cc[i + (k + j*l1)*ido] + cc[i - 1 + (k + jc*l1)*ido];
            ch[i + (k + jc*l1)*ido] = cc[i + (k + j*l1)*ido] - cc[i - 1 + (k + jc*l1)*ido];
          }
        }
      }
    }
    for (ik=0; ik<idl1; ik++) cc[ik] = ch[ik];
    for (j=1; j<ip; j++)
      for (k=0; k<l1; k++)
        cc[(k + j*l1)*ido] = ch[(k + j*l1)*ido];
    if (nbd <= l1) {
      is = -ido;
      for (j=1; j<ip; j++) {
        is += ido;
        idij = is-1;
        for (i=2; i<ido; i+=2) {
          idij += 2;
          for (k=0; k<l1; k++) {
            cc[i - 1 + (k + j*l1)*ido] = wa[idij - 1]*ch[i - 1 + (k + j*l1)*ido] - wa[idij]*
                ch[i + (k + j*l1)*ido];
            cc[i + (k + j*l1)*ido] = wa[idij - 1]*ch[i + (k + j*l1)*ido] + wa[idij]*ch[i - 1 + (k + j*l1)*ido];
          }
        }
      }
    } else {
      is = -ido;
      for (j=1; j<ip; j++) {
        is += ido;
        for (k=0; k<l1; k++) {
          idij = is - 1;
          for (i=2; i<ido; i+=2) {
            idij += 2;
            cc[i - 1 + (k + j*l1)*ido] = wa[idij-1]*ch[i - 1 + (k + j*l1)*ido] - wa[idij]*
                ch[i + (k + j*l1)*ido];
            cc[i + (k + j*l1)*ido] = wa[idij-1]*ch[i + (k + j*l1)*ido] + wa[idij]*ch[i - 1 + (k + j*l1)*ido];
          }
        }
      }
    }
  } /* radbg */

  /* ----------------------------------------------------------------------
cfftf1, cfftf, cfftb, cffti1, cffti. Complex FFTs.
---------------------------------------------------------------------- */

static void cfftf1(int n, Treal c[], Treal ch[], const Treal wa[], const int ifac[MAXFAC+2], int isign)
  {
    int idot, i;
    int k1, l1, l2;
    int na, nf, ip, iw, ix2, ix3, ix4, nac, ido, idl1;
    Treal *cinput, *coutput;
    nf = ifac[1];
    na = 0;
    l1 = 1;
    iw = 0;
    for (k1=2; k1<=nf+1; k1++) {
      ip = ifac[k1];
      l2 = ip*l1;
      ido = n / l2;
      idot = ido + ido;
      idl1 = idot*l1;
      if (na) {
        cinput = ch;
        coutput = c;
      } else {
        cinput = c;
        coutput = ch;
      }
      switch (ip) {
      case 4:
        ix2 = iw + idot;
        ix3 = ix2 + idot;
        passf4(idot, l1, cinput, coutput, &wa[iw], &wa[ix2], &wa[ix3], isign);
        na = !na;
        break;
      case 2:
        passf2(idot, l1, cinput, coutput, &wa[iw], isign);
        na = !na;
        break;
      case 3:
        ix2 = iw + idot;
        passf3(idot, l1, cinput, coutput, &wa[iw], &wa[ix2], isign);
        na = !na;
        break;
      case 5:
        ix2 = iw + idot;
        ix3 = ix2 + idot;
        ix4 = ix3 + idot;
        passf5(idot, l1, cinput, coutput, &wa[iw], &wa[ix2], &wa[ix3], &wa[ix4], isign);
        na = !na;
        break;
      default:
        passf(&nac, idot, ip, l1, idl1, cinput, coutput, &wa[iw], isign);
        if (nac != 0) na = !na;
      }
      l1 = l2;
      iw += (ip - 1)*idot;
    }
    if (na == 0) return;
    for (i=0; i<2*n; i++) c[i] = ch[i];
  } /* cfftf1 */


void cfftf(int n, Treal c[], Treal wsave[])
  {
    int iw1, iw2;
    if (n == 1) return;
    iw1 = 2*n;
    iw2 = iw1 + 2*n;
    cfftf1(n, c, wsave, wsave+iw1, (int*)(wsave+iw2), -1);
  } /* cfftf */


void cfftb(int n, Treal c[], Treal wsave[])
  {
    int iw1, iw2;
    if (n == 1) return;
    iw1 = 2*n;
    iw2 = iw1 + 2*n;
    cfftf1(n, c, wsave, wsave+iw1, (int*)(wsave+iw2), +1);
  } /* cfftb */


static void factorize(int n, int ifac[MAXFAC+2], const int ntryh[NSPECIAL])
  /* Factorize n in factors in ntryh and rest. On exit,
ifac[0] contains n and ifac[1] contains number of factors,
the factors start from ifac[2]. */
  {
    int ntry=3, i, j=0, ib, nf=0, nl=n, nq, nr;
startloop:
    if (j < NSPECIAL)
      ntry = ntryh[j];
    else
      ntry+= 2;
    j++;
    do {
      nq = nl / ntry;
      nr = nl - ntry*nq;
      if (nr != 0) goto startloop;
      nf++;
      ifac[nf + 1] = ntry;
      nl = nq;
      if (ntry == 2 && nf != 1) {
        for (i=2; i<=nf; i++) {
          ib = nf - i + 2;
          ifac[ib + 1] = ifac[ib];
        }
        ifac[2] = 2;
      }
    } while (nl != 1);
    ifac[0] = n;
    ifac[1] = nf;
  }


static void cffti1(int n, Treal wa[], int ifac[MAXFAC+2])
  {
    static const Treal twopi = 6.28318530717959;
    Treal arg, argh, argld, fi;
    int idot, i, j;
    int i1, k1, l1, l2;
    int ld, ii, nf, ip;
    int ido, ipm;

    static const int ntryh[NSPECIAL] = {
      3,4,2,5    }; /* Do not change the order of these. */

    factorize(n,ifac,ntryh);
    nf = ifac[1];
    argh = twopi/(Treal)n;
    i = 1;
    l1 = 1;
    for (k1=1; k1<=nf; k1++) {
      ip = ifac[k1+1];
      ld = 0;
      l2 = l1*ip;
      ido = n / l2;
      idot = ido + ido + 2;
      ipm = ip - 1;
      for (j=1; j<=ipm; j++) {
        i1 = i;
        wa[i-1] = 1;
        wa[i] = 0;
        ld += l1;
        fi = 0;
        argld = ld*argh;
        for (ii=4; ii<=idot; ii+=2) {
          i+= 2;
          fi+= 1;
          arg = fi*argld;
          wa[i-1] = cos(arg);
          wa[i] = sin(arg);
        }
        if (ip > 5) {
          wa[i1-1] = wa[i-1];
          wa[i1] = wa[i];
        }
      }
      l1 = l2;
    }
  } /* cffti1 */


void cffti(int n, Treal wsave[])
 {
    int iw1, iw2;
    if (n == 1) return;
    iw1 = 2*n;
    iw2 = iw1 + 2*n;
    cffti1(n, wsave+iw1, (int*)(wsave+iw2));
  } /* cffti */

  /* ----------------------------------------------------------------------
rfftf1, rfftb1, rfftf, rfftb, rffti1, rffti. Treal FFTs.
---------------------------------------------------------------------- */

static void rfftf1(int n, Treal c[], Treal ch[], const Treal wa[], const int ifac[MAXFAC+2])
  {
    int i;
    int k1, l1, l2, na, kh, nf, ip, iw, ix2, ix3, ix4, ido, idl1;
    Treal *cinput, *coutput;
    nf = ifac[1];
    na = 1;
    l2 = n;
    iw = n-1;
    for (k1 = 1; k1 <= nf; ++k1) {
      kh = nf - k1;
      ip = ifac[kh + 2];
      l1 = l2 / ip;
      ido = n / l2;
      idl1 = ido*l1;
      iw -= (ip - 1)*ido;
      na = !na;
      if (na) {
        cinput = ch;
        coutput = c;
      } else {
        cinput = c;
        coutput = ch;
      }
      switch (ip) {
      case 4:
        ix2 = iw + ido;
        ix3 = ix2 + ido;
        radf4(ido, l1, cinput, coutput, &wa[iw], &wa[ix2], &wa[ix3]);
        break;
      case 2:
        radf2(ido, l1, cinput, coutput, &wa[iw]);
        break;
      case 3:
        ix2 = iw + ido;
        radf3(ido, l1, cinput, coutput, &wa[iw], &wa[ix2]);
        break;
      case 5:
        ix2 = iw + ido;
        ix3 = ix2 + ido;
        ix4 = ix3 + ido;
        radf5(ido, l1, cinput, coutput, &wa[iw], &wa[ix2], &wa[ix3], &wa[ix4]);
        break;
      default:
        if (ido == 1)
          na = !na;
        if (na == 0) {
          radfg(ido, ip, l1, idl1, c, ch, &wa[iw]);
          na = 1;
        } else {
          radfg(ido, ip, l1, idl1, ch, c, &wa[iw]);
          na = 0;
        }
      }
      l2 = l1;
    }
    if (na == 1) return;
    for (i = 0; i < n; i++) c[i] = ch[i];
  } /* rfftf1 */


void rfftb1(int n, Treal c[], Treal ch[], const Treal wa[], const int ifac[MAXFAC+2])
  {
    int i;
    int k1, l1, l2, na, nf, ip, iw, ix2, ix3, ix4, ido, idl1;
    Treal *cinput, *coutput;
    nf = ifac[1];
    na = 0;
    l1 = 1;
    iw = 0;
    for (k1=1; k1<=nf; k1++) {
      ip = ifac[k1 + 1];
      l2 = ip*l1;
      ido = n / l2;
      idl1 = ido*l1;
      if (na) {
        cinput = ch;
        coutput = c;
      } else {
        cinput = c;
        coutput = ch;
      }
      switch (ip) {
      case 4:
        ix2 = iw + ido;
        ix3 = ix2 + ido;
        radb4(ido, l1, cinput, coutput, &wa[iw], &wa[ix2], &wa[ix3]);
        na = !na;
        break;
      case 2:
        radb2(ido, l1, cinput, coutput, &wa[iw]);
        na = !na;
        break;
      case 3:
        ix2 = iw + ido;
        radb3(ido, l1, cinput, coutput, &wa[iw], &wa[ix2]);
        na = !na;
        break;
      case 5:
        ix2 = iw + ido;
        ix3 = ix2 + ido;
        ix4 = ix3 + ido;
        radb5(ido, l1, cinput, coutput, &wa[iw], &wa[ix2], &wa[ix3], &wa[ix4]);
        na = !na;
        break;
      default:
        radbg(ido, ip, l1, idl1, cinput, coutput, &wa[iw]);
        if (ido == 1) na = !na;
      }
      l1 = l2;
      iw += (ip - 1)*ido;
    }
    if (na == 0) return;
    for (i=0; i<n; i++) c[i] = ch[i];
  } /* rfftb1 */


void EMAN::rfftf(int n, Treal r[], Treal wsave[])
  {
    if (n == 1) return;
    rfftf1(n, r, wsave, wsave+n, (int*)(wsave+2*n));
  } /* rfftf */


void EMAN::rfftb(int n, Treal r[], Treal wsave[])
  {
    if (n == 1) return;
    rfftb1(n, r, wsave, wsave+n, (int*)(wsave+2*n));
  } /* rfftb */


static void rffti1(int n, Treal wa[], int ifac[MAXFAC+2])
  {
    static const Treal twopi = 6.28318530717959;
    Treal arg, argh, argld, fi;
    int i, j;
    int k1, l1, l2;
    int ld, ii, nf, ip, is;
    int ido, ipm, nfm1;
    static const int ntryh[NSPECIAL] = {
      4,2,3,5    }; /* Do not change the order of these. */
    factorize(n,ifac,ntryh);
    nf = ifac[1];
    argh = twopi / n;
    is = 0;
    nfm1 = nf - 1;
    l1 = 1;
    if (nfm1 == 0) return;
    for (k1 = 1; k1 <= nfm1; k1++) {
      ip = ifac[k1 + 1];
      ld = 0;
      l2 = l1*ip;
      ido = n / l2;
      ipm = ip - 1;
      for (j = 1; j <= ipm; ++j) {
        ld += l1;
        i = is;
        argld = (Treal) ld*argh;
        fi = 0;
        for (ii = 3; ii <= ido; ii += 2) {
          i += 2;
          fi += 1;
          arg = fi*argld;
          wa[i - 2] = cos(arg);
          wa[i - 1] = sin(arg);
        }
        is += ido;
      }
      l1 = l2;
    }
  } /* rffti1 */


void EMAN::rffti(int n, Treal wsave[])
  {
	if (n == 1) return;
	rffti1(n, wsave+n, (int*)(wsave+2*n));
  } /* rffti */

//#ifdef __cplusplus
//}
//#endif

#endif	//NATIVE_FFT
