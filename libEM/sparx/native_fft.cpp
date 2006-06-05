/**
 * $Id$
 */
#ifdef NATIVE_FFT

#include <cmath>
#include "native_fft.h" 

#include <stdio.h>
#include <stdlib.h>

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

   for (i=1; i<=n; i++) {
      work(i) = 0.0;
   }
   status = fftmcf_(x,work,&n,&n,&n,&inv);
   // should put in appropriate exception handling here
   if (status == -1) {
       fprintf(stderr, "Native FFT cannot be performed on images with the leading ");
       fprintf(stderr, "dimension set to = %d\n", nsam);  
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
       fprintf(stderr, "Native IFT cannot be performed on images with the leading ");
       fprintf(stderr, "dimension set to = %d\n", nsam);  
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

   int   ins, inv, invt, status=0, i, j;
 
   inv = 1;
   ins = inv*lda;

   for (j=1;j<=nrow;j++) {
       status = fmrs_1rf(&x(1,j), work, nsam);
   }

   for (i=1;i<=lda;i=i+2){
      invt = ins;
      status = fftmcf_(&x(i,1),&x(i+1,1),&nrow,&nrow,&nrow,&invt);
   }

   return status;
}

// 2d inplace IFFT
int Nativefft::fmrs_2rb(float *y, float *work, int lda, int nsam, int nrow)
{
   // input:  y(lda,nrow), work(lda)
   // output: y(lda,nrow), overwritten with Fourier coefficients

   int   ins, inv, invt, status=0, i, j;
   float q;

   inv = -1;  
   ins=inv*lda;

   for (i=1;i<=lda;i=i+2){
      invt = ins;
      fftmcf_(&y(i,1),&y(i+1,1),&nrow,&nrow,&nrow,&invt);
   }
    
   // normalize for inverse
   q = 1.0/(float)(nrow);
   for (j=1;j<=nrow;j++) {
       for (i=1;i<=lda;i++) y(i,j)=y(i,j)*q;
   }
 
   for (j=1;j<=nrow;j++) {
      // need an inplace 1d ifft 
      status = fmrs_1rb(&y(1,j), work, nsam);
   }

   return status;
}
#undef x
#undef y

//---------------3D INPLACE FFT----------------------
#define b(i,j,k) b[((k)-1)*lda*nrow + ((j)-1)*lda + (i) - 1]

// 3d inplace FFT
int Nativefft::fmrs_3rf(float *b, float *work, int lda, int nsam, int nrow, int nslice)
{
    // input :  b(lda,nrow,nslice), work(lda)
    // output:  b(lda,nrow,nslice), overwritten with Fourier coefficients

    int ndr, ndrt, inv, i, j, k, status=0;

    inv = 1;
    ndr=inv*lda*nrow;

    for (k=1;k<=nslice;k++) {
        status = fmrs_2rf(&b(1,1,k), work, lda, nsam, nrow);
    } 

    for (j=1;j<=nrow;j++) {
       for (i=1;i<=lda-1;i=i+2) {
	  ndrt=ndr;
          status = fftmcf_(&b(i,j,1),&b(i+1,j,1),&nslice,&nslice,&nslice,&ndrt);
       }
    }
    // should check for error here by looking at ndrt

    return status;
}

// 3d inplace IFFT
int Nativefft::fmrs_3rb(float *b, float *work, 
             int lda, int nsam, int nrow, int nslice)
{
    // input:  b(lda,nrow,nslice), work(lda)
    // output: b(lda,nrow,nslice), overwritten with Fourier coefficients 

    int ndr, ndrt, inv, i, j, k, status=0;
    float q;

    inv = -1;
    ndr=inv*lda*nrow;

    for (j=1;j<=nrow;j++) {
       for (i=1;i<=lda-1;i=i+2) {
	  ndrt=ndr;
          status = fftmcf_(&b(i,j,1),&b(i+1,j,1),&nslice,&nslice,&nslice,&ndrt);
       }
    }

    // should check for error here by looking at ndrt

    // normalize for inverse
    q=1.0/(float)(nslice);
    for (k=1;k<=nslice;k++) 
	for (j=1;j<=nrow;j++)
	    for (i=1;i<=lda;i++)
  	        b(i,j,k)=b(i,j,k)*q;

    for (k=1;k<=nslice;k++) {
       status = fmrs_2rb(&b(1,1,k), work, lda, nsam, nrow);
    }
    return status;
}
#undef b

#endif	//NATIVE_FFT
