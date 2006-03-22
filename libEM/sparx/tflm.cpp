/* TFLM.f -- translated by f2c (version 20020621).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <cerrno>

using namespace std;

namespace {
    /* Common Block Declarations */

    struct commonblock1 {
        float ps, xlambda, cs, contrast, defocus, quadpi, xalf1, xbeta1, alf1, 
              alf2, beta1, beta2, v_ratio__,aa,bb,cc;
        long iswi;
    } ctf1_;

    struct commonblock1& ctf1_1 = ctf1_;
    const bool TRUE_ = 1;
    const bool FALSE_ = 0;

//#define ctf1_1 ctf1_
//#define TRUE_ (1)
//#define FALSE_ (1)

    typedef int /* Unknown procedure type */ (*U_fp)(...);
    typedef /* Subroutine */ int (*S_fp)(...);

    inline double d_sign(double *a, double *b) {
        double x;
        x = (*a >= 0 ? *a : - *a);
        return (*b >= 0 ? x : -x);
    }

    /* Table of constant values */

    static long c__1 = 1;
    static long c__2 = 2;
    static long c__64 = 64;
    static long c__6 = 6;
    static long c__4 = 4;
    static long c__0 = 0;
    static long c__3 = 3;

    int ada_(long *l1, long *nl, long *n, long *nmax,
            long *lpp2, long *iv, double *a, long *inc, double *
            t, double *alf, long *isel);
    int varpro_(long *l, long *nl, long *n, long *
            nmax, long *lpp2, long *iv, double *t, double *y, 
            double *w, U_fp ada, double *a, long *iprint, double *
            alf, double *beta, long *ierr);


    /* CC============================================= */
    /* Subroutine */ int call_lm2__()
    {
        /* System generated locals */
        long i__1, i__2;
        float r__1, r__2, r__3;

        /* Local variables */
        static long l, n;
        static float base_line__, ctf_noise__, aa, bb, cc;
        static long ii, ip, nl, iv;
        static float pu, tx;
        extern /* Subroutine */ int lm2_(float *, long *, long *);
        static float ctf;
        static long lpp2;
        static float tmp3;
        static long nmax;
        static float plot[64];
        static long unit;
        static float xalf2;
        static long unit1;
        extern /* Subroutine */ int sub_base_line__(float *, long *, float *, 
                float *, float *), sub_ctf_noise__(float *, long *);
        static float xbeta2, dlist1, dlist2, dlist3, dlist4, dlist5, dlist6, 
                     dlist7;
        static long lm_ierr__;
        static float xsignal;

        /* C==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==*==* */
        unit = 90;
        unit1 = 91;
        aa = (float)-3.;
        bb = (float)10.8959393613873;
        cc = (float).043;
        ctf1_1.iswi = 2;
        /* L100: */
        /* L101: */
        /* L102: */
        unit = 90;
        ctf1_1.quadpi = (float)3.14159265;
        ctf1_1.iswi = 4;
        ctf1_1.xlambda = (float).017;
        ctf1_1.cs = (float)2e7;
        ctf1_1.defocus = (float)2.4e4;
        ctf1_1.ps = (float)2.3644;
        ctf1_1.contrast = (float).1;
        ctf1_1.v_ratio__ = (float).05;
        n = 64;
        ifstream in("tflm2.dat");
        for (ii = 1; ii <= 64; ++ii) {
            in >> plot[ii-1];
        }
        in.close();
        /*  Number of Linear seperatable parameters Beta */
        l = 1;
        /*  Number of non-linear parameters Alf */
        nl = 0;
        /*  observations */
        /* number of imdependent variables T */
        iv = 1;

        ip = 0;
        lpp2 = l + ip + 2;
        /* Computing MAX */
        i__1 = n, i__2 = (nl << 1) + 3;
        nmax = max(i__1,i__2);
        /* 	32.1086655 -23.4352131  7.52145386  7.72586966 */
        xalf2 = (float)-23.4352131;
        xbeta2 = (float)7.72586966;
        ctf1_1.xalf1 = (float)32.1086655;
        ctf1_1.xbeta1 = (float)7.52145386;
        sub_base_line__(plot, &c__64, &xalf2, &xbeta2, &ctf1_1.ps);
        sub_ctf_noise__(plot, &c__64);
        lm2_(plot, &n, &lm_ierr__);
        for (ii = 1; ii <= 64; ++ii) {
            dlist1 = (float) ii;
            tx = (float) (ii - 1) / (float)64. / (float)2. / ctf1_1.ps;
            /* Computing 2nd power */
            r__1 = tx;
            /* Computing 3rd power */
            r__2 = ctf1_1.xlambda;
            /* Computing 4th power */
            r__3 = tx, r__3 *= r__3;
            ctf = sin(-(ctf1_1.quadpi * (ctf1_1.defocus * ctf1_1.xlambda * (r__1 *
                                r__1) - ctf1_1.cs * (r__2 * (r__2 * r__2)) * (r__3 * r__3) / 
                            (float)2.)) - atan(ctf1_1.contrast / ((float)1. - 
                                                                  ctf1_1.contrast)));
            /* Computing 2nd power */
            r__1 = tx / cc + (float)1.;
            pu = exp(aa + bb / (r__1 * r__1));
            base_line__ = xbeta2 * exp(xalf2 * tx);
            if (ctf1_1.v_ratio__ > (float)1.) {
                /* Computing 2nd power */
                r__1 = tx * ctf1_1.xalf1;
                /* Computing 2nd power */
                r__2 = ctf;
                ctf_noise__ = ctf1_1.beta2 * ctf1_1.xbeta1 * exp(-(r__1 * r__1)) *
                    (r__2 * r__2);
            } else {
                /* Computing 2nd power */
                r__1 = tx * ctf1_1.xalf1;
                /* Computing 2nd power */
                r__2 = ctf;
                ctf_noise__ = ctf1_1.v_ratio__ * ctf1_1.xbeta1 * exp(-(r__1 * 
                            r__1)) * (r__2 * r__2);
            }
            /* Computing 2nd power */
            r__1 = tx * ctf1_1.xalf1;
            /* Computing 2nd power */
            r__2 = ctf;
            xsignal = ctf1_1.beta1 * pu * exp(-(r__1 * r__1)) * (r__2 * r__2);
            dlist2 = tx;
            dlist3 = base_line__;
            dlist4 = ctf_noise__;
            dlist5 = xsignal;
            dlist6 = xsignal + ctf_noise__;
            dlist7 = plot[ii - 1];
            ofstream out91("fort.91");
            out91 << dlist1 << " " << dlist2 << " " << dlist3 << " "
                << dlist4 << " " << dlist5 << " " << dlist6 << " "
                << dlist7 << endl;
        }
        tmp3 = (float) lm_ierr__;
        /* c  XALF1  B-factor */
        /* c  Xbeta1 scaling factor of envelope function */
        /* c  XALF2 base line noise ax */
        /* c  Xbeta2  base line b=log(Xbeta2) */
        /* c  Beta1 scaling factor of the model Pu */
        /* Computing 2nd power */
        r__1 = ctf1_1.xalf1;
        ctf1_1.xalf1 = r__1 * r__1;
        cout << ctf1_1.xalf1 << " " << xalf2 << " "
            << ctf1_1.xbeta1 << " " << xbeta2 << " "
            << ctf1_1.beta1 << endl;
        return 0;
    } /* call_lm2__ */

    /* C====================================================================== */
    /* Subroutine */ int sub_base_line__(float *plot, long *n, float *xalf2, 
            float *xbeta2, float *xps)
    {
        /* System generated locals */
        long i__1;

        /* Local variables */
        static long ii;
        static float tmpx, tmpy;

        /* Parameter adjustments */
        --plot;

        /* Function Body */
        i__1 = *n;
        for (ii = 1; ii <= i__1; ++ii) {
            tmpx = (float) (ii - 1) / (float)2. / (float) (*n) / *xps;
            tmpy = *xbeta2 * exp(*xalf2 * tmpx);
            plot[ii] -= tmpy;
        }
        return 0;
    } /* sub_base_line__ */

    /* c================================================================= */
    /* Subroutine */ int sub_ctf_noise__(float *plot, long *n)
    {
        /* System generated locals */
        long i__1;
        float r__1, r__2, r__3;

        /* Local variables */
        static float ctf_noise__;
        static long ii;
        static float tx, ctf;

        /* Parameter adjustments */
        --plot;

        /* Function Body */
        i__1 = *n;
        for (ii = 1; ii <= i__1; ++ii) {
            tx = (float) (ii - 1) / (float) (*n) / (float)2. / ctf1_1.ps;
            /* Computing 2nd power */
            r__1 = tx;
            /* Computing 3rd power */
            r__2 = ctf1_1.xlambda;
            /* Computing 4th power */
            r__3 = tx, r__3 *= r__3;
            ctf = sin(-(ctf1_1.quadpi * (ctf1_1.defocus * ctf1_1.xlambda * (r__1 *
                                r__1) - ctf1_1.cs * (r__2 * (r__2 * r__2)) * (r__3 * r__3) / 
                            (float)2.)) - atan(ctf1_1.contrast / ((float)1. - 
                                                                  ctf1_1.contrast)));
            /* Computing 2nd power */
            r__1 = tx * ctf1_1.xalf1;
            /* Computing 2nd power */
            r__2 = ctf;
            ctf_noise__ = ctf1_1.xbeta1 * exp(-(r__1 * r__1)) * (r__2 * r__2);
            plot[ii] -= ctf_noise__ * ctf1_1.v_ratio__;
            if (plot[ii] < (float)0.) {
                plot[ii] = (float)0.;
            }
        }
        return 0;
    } /* sub_ctf_noise__ */

    /* C================================================================= */
    /* Test levenberg-Marqardt optimization algorithm */
    /* 			Zhong Huang, March, 3, 2005 */
    /* Subroutine */ int lm2_(float *plot, long *n, long *lm_ierr__)
    {
        /* System generated locals */
        long i__1;
        float r__1, r__2, r__3;
        double d__1, d__2;

        /* Local variables */
        static double a[192]	/* was [64][3] */, t[64]	/* was [64][1]
                                                                 */, w[64], y[64];
            static long ii;
        static double alf;
        static double beta[1];
        static long ierr;
        static float tmpx;
        static long iprint;

        /*  Number of Linear seperatable parameters Beta */
        /* 	L=1 */
        /*  Number of non-linear parameters Alf */
        /* 	NL=0 */
        /*  observations */
        /* number of imdependent variables T */
        /* 	IV=1 */

        /* 	IP=0 */
        /* 	LPP2=L+IP+2 */
        /* 	NMAX=MAX(N,2*NL+3) */
        /* Parameter adjustments */
        --plot;

        /* Function Body */
        iprint = -1;
        /* simple test, we assume the observation is from one simple gaussian fucniton */
        /* c define weighting function */
        /* cWe always discuss functions within [0,.2] only! */
        i__1 = *n;
        for (ii = 1; ii <= i__1; ++ii) {
            t[ii - 1] = (double) (ii - 1) / (double) (*n) / (float)2. / (
                    double) ctf1_1.ps;
        }
        i__1 = *n;
        for (ii = 1; ii <= i__1; ++ii) {
            y[ii - 1] = (double) plot[ii];
        }
        i__1 = *n;
        for (ii = 1; ii <= i__1; ++ii) {
            /* Computing 2nd power */
            r__1 = (float) t[ii - 1];
            /* Computing 3rd power */
            r__2 = ctf1_1.xlambda;
            /* Computing 4th power */
            r__3 = (float) t[ii - 1], r__3 *= r__3;
            tmpx = sin(-(ctf1_1.quadpi * (ctf1_1.defocus * ctf1_1.xlambda * (r__1 
                                * r__1) - ctf1_1.cs * (r__2 * (r__2 * r__2)) * (r__3 * r__3) /
                            (float)2.)) - atan(ctf1_1.contrast / ((float)1. - 
                                                                  ctf1_1.contrast)));
            /* Computing 2nd power */
            d__1 = (double) tmpx;
            /* Computing 2nd power */
            d__2 = t[ii - 1] * ctf1_1.xalf1;
            w[ii - 1] = d__1 * d__1 * exp(-(d__2 * d__2));
        }
        varpro_(&c__1, &c__0, n, &c__64, &c__3, &c__1, t, y, w, (U_fp)ada_, a, &
                iprint, &alf, beta, &ierr);
        *lm_ierr__ = ierr;
        ctf1_1.beta1 = (float) beta[0];
        return 0;
    } /* lm2_ */

    /* CC============================================== */
    /* Subroutine */ int ada_(long *l1, long *nl, long *n, long *nmax,
            long *lpp2, long *iv, double *a, long *inc, double *
            t, double *alf, long *isel)
    {
        /* System generated locals */
        long a_dim1, a_offset, t_dim1, t_offset, i__1;
        float r__1, r__2, r__3;
        double d__1, d__2, d__3;

        /* Local variables */
        double* b = NULL;
        double* ctf = NULL;
        static double aa, bb, cc;
        static long ii, jj;
        static double pu;
        static float tmpx;

        /* Parameter adjustments */
        --alf;
        a_dim1 = *nmax;
        a_offset = 1 + a_dim1 * 1;
        a -= a_offset;
        t_dim1 = *nmax;
        t_offset = 1 + t_dim1 * 1;
        t -= t_offset;
        inc -= 13;
        b = new double [(*nmax)*(*l1)];  // Alloca
        ctf = new double [*n];
        /* Function Body */
        if (ctf1_1.iswi == 4) {
            i__1 = *n;
            for (ii = 1; ii <= i__1; ++ii) {
                /* Computing 2nd power */
                r__1 = (float) t[ii + t_dim1];
                /* Computing 3rd power */
                r__2 = ctf1_1.xlambda;
                /* Computing 4th power */
                r__3 = (float) t[ii + t_dim1], r__3 *= r__3;
                tmpx = sin(-(ctf1_1.quadpi * (ctf1_1.defocus * ctf1_1.xlambda * (
                                    r__1 * r__1) - ctf1_1.cs * (r__2 * (r__2 * r__2)) * (r__3 
                                        * r__3) / (float)2.)) - atan(ctf1_1.contrast / ((float)1. 
                                        - ctf1_1.contrast)));
                ctf[ii - 1] = (double) tmpx;
            }
            if (*isel == 1) {
                for (ii = 1; ii <= 12; ++ii) {
                    for (jj = 1; jj <= 8; ++jj) {
                        inc[ii + jj * 12] = 0;
                    }
                }
                inc[13] = 1;
                inc[26] = 1;
                /* Initilization WE set BB=6, CC=1. and start optimization */
                i__1 = *n;
                for (ii = 1; ii <= i__1; ++ii) {
                    a[ii + (a_dim1 << 1)] = exp(alf[2] * t[ii + t_dim1]);
                    /* Computing 2nd power */
                    d__1 = t[ii + t_dim1] * alf[1];
                    /* Computing 2nd power */
                    d__2 = ctf[ii - 1];
                    a[ii + a_dim1] = exp(-(d__1 * d__1)) * (d__2 * d__2);
                    a[ii + a_dim1 * 5] = exp(alf[2] * t[ii + t_dim1]) * t[ii + 
                        t_dim1];
                    /* Computing 2nd power */
                    d__1 = t[ii + t_dim1] * alf[1];
                    /* Computing 2nd power */
                    d__2 = t[ii + t_dim1];
                    /* Computing 2nd power */
                    d__3 = ctf[ii - 1];
                    a[ii + (a_dim1 << 2)] = exp(-(d__1 * d__1)) * (float)-2. * 
                        alf[1] * (d__2 * d__2) * (d__3 * d__3);
                }
            }
            if (*isel == 2) {
                i__1 = *n;
                for (ii = 1; ii <= i__1; ++ii) {
                    a[ii + (a_dim1 << 1)] = exp(alf[2] * t[ii + t_dim1]);
                    /* Computing 2nd power */
                    d__1 = t[ii + t_dim1] * alf[1];
                    /* Computing 2nd power */
                    d__2 = ctf[ii - 1];
                    a[ii + a_dim1] = exp(-(d__1 * d__1)) * (d__2 * d__2);
                    /* Computing 2nd power */
                    d__1 = t[ii + t_dim1] * alf[1];
                    b[ii - 1] = exp(-(d__1 * d__1));
                    b[ii + (*n)-1] = exp(alf[2] * t[ii + t_dim1]);
                }
            }
            if (*isel == 3) {
                i__1 = *n;
                for (ii = 1; ii <= i__1; ++ii) {
                    a[ii + a_dim1 * 5] = b[ii +(*n)-1 ] * t[ii + t_dim1];
                    /* Computing 2nd power */
                    d__1 = t[ii + t_dim1];
                    a[ii + (a_dim1 << 2)] = b[ii - 1] * (float)-2. * alf[1] * (
                            d__1 * d__1);
                }
            }
        } else {
            aa = ctf1_1.aa ;    /*(float)-3.;*/;
            bb = ctf1_1.bb ;   /* (float)10.8959393613873; */
            cc = ctf1_1.cc ; /*(float).043; */
            i__1 = *n;
            for (ii = 1; ii <= i__1; ++ii) {
                /* Computing 2nd power */
                r__1 = (float) t[ii + t_dim1];
                /* Computing 3rd power */
                r__2 = ctf1_1.xlambda;
                /* Computing 4th power */
                r__3 = (float) t[ii + t_dim1], r__3 *= r__3;
                tmpx = sin(-(ctf1_1.quadpi * (ctf1_1.defocus * ctf1_1.xlambda * (
                                    r__1 * r__1) - ctf1_1.cs * (r__2 * (r__2 * r__2)) * (r__3 
                                        * r__3) / (float)2.)) - atan(ctf1_1.contrast / ((float)1. 
                                        - ctf1_1.contrast)));
                ctf[ii - 1] = (double) tmpx;
            }
            if (*isel == 1) {
                for (ii = 1; ii <= 12; ++ii) {
                    for (jj = 1; jj <= 8; ++jj) {
                        inc[ii + jj * 12] = 0;
                    }
                }
                /* Initilization WE set BB=6, CC=1. and start optimization */
                i__1 = *n;
                for (ii = 1; ii <= i__1; ++ii) {
                    /* Computing 2nd power */
                    d__1 = t[ii + t_dim1] / cc + (float)1.;
                    pu = exp(aa + bb / (d__1 * d__1));
                    /* Computing 2nd power */
                    d__1 = t[ii + t_dim1] * ctf1_1.xalf1;
                    /* Computing 2nd power */
                    d__2 = ctf[ii - 1];
                    a[ii + a_dim1] = pu * exp(-(d__1 * d__1)) * (d__2 * d__2);
                    /* Computing 2nd power */
                    d__1 = t[ii + t_dim1] * ctf1_1.xalf1;
                    /* Computing 2nd power */
                    d__2 = ctf[ii - 1];
                    a[ii + (a_dim1 << 1)] = ctf1_1.xbeta1 * exp(-(d__1 * d__1)) * 
                        (d__2 * d__2);
                }
            }
        }
        delete [] b;
        delete [] ctf;
        b = NULL;
        ctf =NULL;
        return 0;
    } /* ada_ */

    /* cc================================================================= */
    /* Subroutine */ int varpro_(long *l, long *nl, long *n, long *
            nmax, long *lpp2, long *iv, double *t, double *y, 
            double *w, U_fp ada, double *a, long *iprint, double *
            alf, double *beta, long *ierr)
    {
        /* Initialized data */

        static double eps1 = 1e-8;
        static long itmax = 50;

        /* System generated locals */
        long a_dim1, a_offset, t_dim1, t_offset, i__1, i__2;


        /* Local variables */
        static long j, k;
        static double r__;
        static long b1;
        static double nu;
        static long lp1, lnl2, nlp1;
        static double acum;
        static long isub, iter, ksub;
        static bool skip;
        static long jsub;
        static double rnew;
        extern /* Subroutine */ int vpdpa_(long *, long *, long *, 
                long *, long *, long *, double *, double *, 
                double *, double *, U_fp, long *, long *, 
                double *, double *, double *, double *);
        static long modit;
        extern /* Subroutine */ int vperr_(long *, long *, long *), 
               vpfac1_(long *, long *, long *, long *, long *, 
                       double *, double *, long *), vpfac2_(long *, 
                           long *, double *, double *);
        static long iterin;
        static double gnstep, prjres;
        extern /* Subroutine */ int vpbsol_(long *, long *, double *, 
                double *);
        extern double vpnorm_(long *, double *);
        extern /* Subroutine */ int vppost_(long *, long *, long *, 
                long *, long *, double *, double *, long *, 
                double *, double *, double *, double *, 
                double *, long *);

        /*           GIVEN A SET OF N OBSERVATIONS, CONSISTING OF VALUES Y(1), */
        /*        Y(2), ..., Y(N) OF A DEPENDENT VARIABLE Y, WHERE Y(I) */
        /*        CORRESPONDS TO THE IV INDEPENDENT VARIABLE(S) T(I,1), T(I,2), */
        /*        ..., T(I,IV), VARPRO ATTEMPTS TO COMPUTE A WEIGHTED LEAST */
        /*        SQUARES FIT TO A FUNCTION ETA (THE 'MODEL') WHICH IS A LINEAR */
        /*        COMBINATION */
        /*                              L                                        CTF1 */
        /*       ETA(ALF, BETA; T)  =  SUM  BETA  * PHI (ALF; T) + PHI   (ALF; T) */
        /*                             J=1      J      J              L+1 */

        /*        OF NONLINEAR FUNCTIONS PHI(J) (E.G., A SUM OF EXPONENTIALS AND/ */
        /*        OR GAUSSIANS).  THAT IS, DETERMINE THE LINEAR PARAMETERS */
        /*        BETA(J) AND THE VECTOR OF NONLINEAR PARAMETERS ALF BY MINIMIZ- */
        /*        ING */

        /*                         2     N                                 2 */
        /*           NORM(RESIDUAL)  =  SUM  W  * (Y  - ETA(ALF, BETA; T )) . */
        /*                              I=1   I     I                   I */

        /*        THE (L+1)-ST TERM IS OPTIONAL, AND IS USED WHEN IT IS DESIRED */
        /*        TO FIX ONE OR MORE OF THE BETA'S (RATHER THAN LET THEM BE */
        /*        DETERMINED).  VARPRO REQUIRES FIRST DERIVATIVES OF THE PHI'S. */

        /*                                NOTES: */

        /*        A)  THE ABOVE PROBLEM IS ALSO REFERRED TO AS 'MULTIPLE */
        /*        NONLINEAR REGRESSION'.  FOR USE IN STATISTICAL ESTIMATION, */
        /*        VARPRO RETURNS THE RESIDUALS, THE COVARIANCE MATRIX OF THE */
        /*        LINEAR AND NONLINEAR PARAMETERS, AND THE ESTIMATED VARIANCE OF */
        /*        THE OBSERVATIONS. */

        /*        B) AN ETA OF THE ABOVE FORM IS CALLED 'SEPARABLE'.  THE */
        /*        CASE OF A NONSEPARABLE ETA CAN BE HANDLED BY SETTING L = 0 */
        /*        AND USING PHI(L+1). */

        /*        C) VARPRO MAY ALSO BE USED TO SOLVE LINEAR LEAST SQUARES */
        /*        PROBLEMS (IN THAT CASE NO ITERATIONS ARE PERFORMED).  SET */
        /*        NL = 0. */

        /*        D)  THE MAIN ADVANTAGE OF VARPRO OVER OTHER LEAST SQUARES */
        /*        PROGRAMS IS THAT NO INITIAL GUESSES ARE NEEDED FOR THE LINEAR */
        /*        PARAMETERS.  NOT ONLY DOES THIS MAKE IT EASIER TO USE, BUT IT */
        /*        OFTEN LEADS TO FASTER CONVERGENCE. */


        /*     DESCRIPTION OF PARAMETERS */

        /*        L       NUMBER OF LINEAR PARAMETERS BETA (MUST BE .GE. 0). */
        /*        NL      NUMBER OF NONLINEAR PARAMETERS ALF (MUST BE .GE. 0). */
        /*        N       NUMBER OF OBSERVATIONS.  N MUST BE GREATER THAN L + NL */
        /*                (I.E., THE NUMBER OF OBSERVATIONS MUST EXCEED THE */
        /*                NUMBER OF PARAMETERS). */
        /*        IV      NUMBER OF INDEPENDENT VARIABLES T. */
        /*        T       float N BY IV MATRIX OF INDEPENDENT VARIABLES.  T(I, J) */
        /*                CONTAINS THE VALUE OF THE I-TH OBSERVATION OF THE J-TH */
        /*                INDEPENDENT VARIABLE. */
        /*        Y       N-VECTOR OF OBSERVATIONS, ONE FOR EACH ROW OF T. */
        /*        W       N-VECTOR OF NONNEGATIVE WEIGHTS.  SHOULD BE SET TO 1'S */
        /*                IF WEIGHTS ARE NOT DESIRED.  IF VARIANCES OF THE */
        /*                INDIVIDUAL OBSERVATIONS ARE KNOWN, W(I) SHOULD BE SET */
        /*                TO 1./VARIANCE(I). */
        /*        INC     NL X (L+1) long INCIDENCE MATRIX.  INC(K, J) = 1 IF */
        /*                NON-LINEAR PARAMETER ALF(K) APPEARS IN THE J-TH */
        /*                FUNCTION PHI(J).  (THE PROGRAM SETS ALL OTHER INC(K, J) */
        /*                TO ZERO.)  IF PHI(L+1) IS INCLUDED IN THE MODEL, */
        /*                THE APPROPRIATE ELEMENTS OF THE (L+1)-ST COLUMN SHOULD */
        /*                BE SET TO 1'S.  INC IS NOT NEEDED WHEN L = 0 OR NL = 0. */
        /*                CAUTION:  THE DECLARED ROW DIMENSION OF INC (IN ADA) */
        /*                MUST CURRENTLY BE SET TO 12.  SEE 'RESTRICTIONS' BELOW. */
        /*        NMAX    THE DECLARED ROW DIMENSION OF THE MATRICES A AND T. */
        /*                IT MUST BE AT LEAST MAX(N, 2*NL+3). */
        /*        LPP2    L+P+2, WHERE P IS THE NUMBER OF ONES IN THE MATRIX INC. */
        /*                THE DECLARED COLUMN DIMENSION OF A MUST BE AT LEAST */
        /*                LPP2.  (IF L = 0, SET LPP2 = NL+2. IF NL = 0, SET LPP2 */
        /*                L+2.) */
        /*        A       float MATRIX OF SIZE MAX(N, 2*NL+3) BY L+P+2.  ON INPUT */
        /*                IT CONTAINS THE PHI(J)'S AND THEIR DERIVATIVES (SEE */
        /*                BELOW).  ON OUTPUT, THE FIRST L+NL ROWS AND COLUMNS OF */
        /*                A WILL CONTAIN AN APPROXIMATION TO THE (WEIGHTED) */
        /*                COVARIANCE MATRIX AT THE SOLUTION (THE FIRST L ROWS */
        /*                CORRESPOND TO THE LINEAR PARAMETERS, THE LAST NL TO THE */
        /*                NONLINEAR ONES), COLUMN L+NL+1 WILL CONTAIN THE */
        /*                WEIGHTED RESIDUALS (Y - ETA), A(1, L+NL+2) WILL CONTAIN */
        /*                THE (EUCLIDEAN) NORM OF THE WEIGHTED RESIDUAL, AND */
        /*                A(2, L+NL+2) WILL CONTAIN AN ESTIMATE OF THE (WEIGHTED) */
        /*                VARIANCE OF THE OBSERVATIONS, NORM(RESIDUAL)**2/ */
        /*                (N - L - NL). */
        /*        IPRINT  INPUT long CONTROLLING PRINTED OUTPUT.  IF IPRINT IS */
        /*                POSITIVE, THE NONLINEAR PARAMETERS, THE NORM OF THE */
        /*                RESIDUAL, AND THE MARQUARDT PARAMETER WILL BE OUTPUT */
        /*                EVERY IPRINT-TH ITERATION (AND INITIALLY, AND AT THE */
        /*                FINAL ITERATION).  THE LINEAR PARAMETERS WILL BE */
        /*                PRINTED AT THE FINAL ITERATION.  ANY ERROR MESSAGES */
        /*                WILL ALSO BE PRINTED.  (IPRINT = 1 IS RECOMMENDED AT */
        /*                FIRST.) IF IPRINT = 0, ONLY THE FINAL QUANTITIES WILL */
        /*                BE PRINTED, AS WELL AS ANY ERROR MESSAGES.  IF IPRINT = */
        /*                -1, NO PRINTING WILL BE DONE.  THE USER IS THEN */
        /*                RESPONSIBLE FOR CHECKING THE PARAMETER IERR FOR ERRORS. */
        /*        ALF     NL-VECTOR OF ESTIMATES OF NONLINEAR PARAMETERS */
        /*                (INPUT).  ON OUTPUT IT WILL CONTAIN OPTIMAL VALUES OF */
        /*                THE NONLINEAR PARAMETERS. */
        /*        BETA    L-VECTOR OF LINEAR PARAMETERS (OUTPUT ONLY). */
        /*        IERR    long ERROR FLAG (OUTPUT): */
        /*                .GT. 0 - SUCCESSFUL CONVERGENCE, IERR IS THE NUMBER OF */
        /*                    ITERATIONS TAKEN. */
        /*                -1  TERMINATED FOR TOO MANY ITERATIONS. */
        /*                -2  TERMINATED FOR ILL-CONDITIONING (MARQUARDT */
        /*                    PARAMETER TOO LARGE.)  ALSO SEE IERR = -8 BELOW. */
        /*                -4  INPUT ERROR IN PARAMETER N, L, NL, LPP2, OR NMAX. */
        /*                -5  INC MATRIX IMPROPERLY SPECIFIED, OR P DISAGREES */
        /*                    WITH LPP2. */
        /*                -6  A WEIGHT WAS NEGATIVE. */
        /*                -7  'CONSTANT' COLUMN WAS COMPUTED MORE THAN ONCE. */
        /*                -8  CATASTROPHIC FAILURE - A COLUMN OF THE A MATRIX HAS */
        /*                    BECOME ZERO.  SEE 'CONVERGENCE FAILURES' BELOW. */

        /*                (IF IERR .LE. -4, THE LINEAR PARAMETERS, COVARIANCE */
        /*                MATRIX, ETC. ARE NOT RETURNED.) */

        /*     SUBROUTINES REQUIRED */

        /*           NINE SUBROUTINES, VPDPA, VPFAC1, VPFAC2, VPBSOL, VPPOST, */
        /*        VPCOV ,VPNORM, VPINIT, AND VPERR ARE PROVIDED.  IN ADDITION, */
        /*        THE USER MUST PROVIDE A SUBROUTINE (CORRESPONDING TO THE ARGU- */
        /*        MENT ADA) WHICH, GIVEN ALF, WILL EVALUATE THE FUNCTIONS PHI(J) */
        /*        AND THEIR PARTIAL DERIVATIVES D PHI(J)/D ALF(K), AT THE SAMPLE */
        /*        POINTS T(I).  THIS ROUTINE MUST BE DECLARED 'EXTERNAL' IN THE */
        /*        CALLING PROGRAM.  ITS CALLING SEQUENCE IS */

        /*        SUBROUTINE ADA (L+1, NL, N, NMAX, LPP2, IV, A, INC, T, ALF, */
        /*        ISEL) */

        /*           THE USER SHOULD MODIFY THE EXAMPLE SUBROUTINE 'ADA' (GIVEN */
        /*        ELSEWHERE) FOR HIS OWN FUNCTIONS. */

        /*           THE VECTOR SAMPLED FUNCTIONS PHI(J) SHOULD BE STORED IN THE */
        /*        FIRST N ROWS AND FIRST L+1 COLUMNS OF THE MATRIX A, I.E., */
        /*        A(I, J) SHOULD CONTAIN PHI(J, ALF; T(I,1), T(I,2), ..., */
        /*        T(I,IV)), I = 1, ..., N; J = 1, ..., L (OR L+1).  THE (L+1)-ST */
        /*        COLUMN OF A CONTAINS PHI(L+1) IF PHI(L+1) IS IN THE MODEL, */
        /*        OTHERWISE IT IS RESERVED FOR WORKSPACE.  THE 'CONSTANT' FUNC- */
        /*        TIONS (THESE ARE FUNCTIONS PHI(J) WHICH DO NOT DEPEND UPON ANY */
        /*        NONLINEAR PARAMETERS ALF, E.G., T(I)**J) (IF ANY) MUST APPEAR */
        /*        FIRST, STARTING IN COLUMN 1.  THE COLUMN N-VECTORS OF NONZERO */
        /*        PARTIAL DERIVATIVES D PHI(J) / D ALF(K) SHOULD BE STORED */
        /*        SEQUENTIALLY IN THE MATRIX A IN COLUMNS L+2 THROUGH L+P+1. */
        /*        THE ORDER IS */

        /*          D PHI(1)  D PHI(2)        D PHI(L)  D PHI(L+1)  D PHI(1) */
        /*          --------, --------, ...,  --------, ----------, --------, */
        /*          D ALF(1)  D ALF(1)        D ALF(1)   D ALF(1)   D ALF(2) */

        /*          D PHI(2)       D PHI(L+1)       D PHI(1)        D PHI(L+1) */
        /*          --------, ..., ----------, ..., ---------, ..., ----------, */
        /*          D ALF(2)        D ALF(2)        D ALF(NL)       D ALF(NL) */

        /*        OMITTING COLUMNS OF DERIVATIVES WHICH ARE ZERO, AND OMITTING */
        /*        PHI(L+1) COLUMNS IF PHI(L+1) IS NOT IN THE MODEL.  NOTE THAT */
        /*        THE LINEAR PARAMETERS BETA ARE NOT USED IN THE MATRIX A. */
        /*        COLUMN L+P+2 IS RESERVED FOR WORKSPACE. */

        /*        THE CODING OF ADA SHOULD BE ARRANGED SO THAT: */

        /*        ISEL = 1  (WHICH OCCURS THE FIRST TIME ADA IS CALLED) MEANS: */
        /*                  A.  FILL IN THE INCIDENCE MATRIX INC */
        /*                  B.  STORE ANY CONSTANT PHI'S IN A. */
        /*                  C.  COMPUTE NONCONSTANT PHI'S AND PARTIAL DERIVA- */
        /*                      TIVES. */
        /*             = 2  MEANS COMPUTE ONLY THE NONCONSTANT FUNCTIONS PHI */
        /*             = 3  MEANS COMPUTE ONLY THE DERIVATIVES */

        /*        (WHEN THE PROBLEM IS LINEAR (NL = 0) ONLY ISEL = 1 IS USED, AND */
        /*        DERIVATIVES ARE NOT NEEDED.) */

        /*     RESTRICTIONS */

        /*           THE SUBROUTINES VPDPA, VPINIT (AND ADA) CONTAIN THE LOCALLY */
        /*        DIMENSIONED MATRIX INC, WHOSE DIMENSIONS ARE CURRENTLY SET FOR */
        /*        MAXIMA OF L+1 = 8, NL = 12.  THEY MUST BE CHANGED FOR LARGER */
        /*        PROBLEMS.  DATA PLACED IN ARRAY A IS OVERWRITTEN ('DESTROYED'). */
        /*        DATA PLACED IN ARRAYS T, Y AND INC IS LEFT INTACT.  THE PROGRAM */
        /*        RUNS IN WATFIV, EXCEPT WHEN L = 0 OR NL = 0. */

        /*           IT IS ASSUMED THAT THE MATRIX PHI(J, ALF; T(I)) HAS FULL */
        /*        COLUMN RANK.  THIS MEANS THAT THE FIRST L COLUMNS OF THE MATRIX */
        /*        A MUST BE LINEARLY INDEPENDENT. */

        /*           OPTIONAL NOTE:  AS WILL BE NOTED FROM THE SAMPLE SUBPROGRAM */
        /*        ADA, THE DERIVATIVES D PHI(J)/D ALF(K) (ISEL = 3) MUST BE */
        /*        COMPUTED INDEPENDENTLY OF THE FUNCTIONS PHI(J) (ISEL = 2), */
        /*        SINCE THE FUNCTION VALUES ARE OVERWRITTEN AFTER ADA IS CALLED */
        /*        WITH ISEL = 2.  THIS IS DONE TO MINIMIZE STORAGE, AT THE POS- */
        /*        SIBLE EXPENSE OF SOME RECOMPUTATION (SINCE THE FUNCTIONS AND */
        /*        DERIVATIVES FREQUENTLY HAVE SOME COMMON SUBEXPRESSIONS).  TO */
        /*        REDUCE THE AMOUNT OF COMPUTATION AT THE EXPENSE OF SOME */
        /*        STORAGE, CREATE A MATRIX B OF DIMENSION NMAX BY L+1 IN ADA, AND */
        /*        AFTER THE COMPUTATION OF THE PHI'S (ISEL = 2), COPY THE VALUES */
        /*        INTO B.  THESE VALUES CAN THEN BE USED TO CALCULATE THE DERIV- */
        /*        ATIVES (ISEL = 3).  (THIS MAKES USE OF THE FACT THAT WHEN A */
        /*        CALL TO ADA WITH ISEL = 3 FOLLOWS A CALL WITH ISEL = 2, THE */
        /*        ALFS ARE THE SAME.) */

        /*           TO CONVERT TO OTHER MACHINES, CHANGE THE OUTPUT UNIT IN THE */
        /*        DATA STATEMENTS IN VARPRO, VPDPA, VPPOST, AND VPERR.  THE */
        /*        PROGRAM HAS BEEN CHECKED FOR PORTABILITY BY THE BELL LABS PFORT */
        /*        VERIFIER.  FOR MACHINES WITHOUT DOUBLE PRECISION HARDWARE, IT */
        /*        MAY BE DESIRABLE TO CONVERT TO SINGLE PRECISION.  THIS CAN BE */
        /*        DONE BY CHANGING (A) THE DECLARATIONS 'DOUBLE PRECISION' TO */
        /*        'float', (B) THE PATTERN '.D' TO '.E' IN THE 'DATA' STATEMENT IN */
        /*        VARPRO, (C) DSIGN, DSQRT AND DABS TO SIGN, SQRT AND ABS, */
        /*        RESPECTIVELY, AND (D) DEXP TO EXP IN THE SAMPLE PROGRAMS ONLY. */

        /*     NOTE ON INTERPRETATION OF COVARIANCE MATRIX */

        /*           FOR USE IN STATISTICAL ESTIMATION (MULTIPLE NONLINEAR */
        /*        REGRESSION) VARPRO RETURNS THE COVARIANCE MATRIX OF THE LINEAR */
        /*        AND NONLINEAR PARAMETERS.  THIS MATRIX WILL BE USEFUL ONLY IF */
        /*        THE USUAL STATISTICAL ASSUMPTIONS HOLD:  AFTER WEIGHTING, THE */
        /*        ERRORS IN THE OBSERVATIONS ARE INDEPENDENT AND NORMALLY DISTRI- */
        /*        BUTED, WITH MEAN ZERO AND THE SAME VARIANCE.  IF THE ERRORS DO */
        /*        NOT HAVE MEAN ZERO (OR ARE UNKNOWN), THE PROGRAM WILL ISSUE A */
        /*        WARNING MESSAGE (UNLESS IPRINT .LT. 0) AND THE COVARIANCE */
        /*        MATRIX WILL NOT BE VALID.  IN THAT CASE, THE MODEL SHOULD BE */
        /*        ALTERED TO INCLUDE A CONSTANT TERM (SET PHI(1) = 1.). */

        /*           NOTE ALSO THAT, IN ORDER FOR THE USUAL ASSUMPTIONS TO HOLD, */
        /*        THE OBSERVATIONS MUST ALL BE OF APPROXIMATELY THE SAME */
        /*        MAGNITUDE (IN THE ABSENCE OF INFORMATION ABOUT THE ERROR OF */
        /*        EACH OBSERVATION), OTHERWISE THE VARIANCES WILL NOT BE THE */
        /*        SAME.  IF THE OBSERVATIONS ARE NOT THE SAME SIZE, THIS CAN BE */
        /*        CURED BY WEIGHTING. */

        /*           IF THE USUAL ASSUMPTIONS HOLD, THE SQUARE ROOTS OF THE */
        /*        DIAGONALS OF THE COVARIANCE MATRIX A GIVE THE STANDARD ERROR */
        /*        S(I) OF EACH PARAMETER.  DIVIDING A(I,J) BY S(I)*S(J) YIELDS */
        /*        THE CORRELATION MATRIX OF THE PARAMETERS.  PRINCIPAL AXES AND */
        /*        CONFIDENCE ELLIPSOIDS CAN BE OBTAINED BY PERFORMING AN EIGEN- */
        /*        VALUE/EIGENVECTOR ANALYSIS ON A.  ONE SHOULD CALL THE EISPACK */
        /*        PROGRAM TRED2, FOLLOWED BY TQL2 (OR USE THE EISPAC CONTROL */
        /*        PROGRAM). */

        /*     CONVERGENCE FAILURES */

        /*           IF CONVERGENCE FAILURES OCCUR, FIRST CHECK FOR INCORRECT */
        /*        CODING OF THE SUBROUTINE ADA.  CHECK ESPECIALLY THE ACTION OF */
        /*        ISEL, AND THE COMPUTATION OF THE PARTIAL DERIVATIVES.  IF THESE */
        /*        ARE CORRECT, TRY SEVERAL STARTING GUESSES FOR ALF.  IF ADA */
        /*        IS CODED CORRECTLY, AND IF ERROR RETURNS IERR = -2 OR -8 */
        /*        PERSISTENTLY OCCUR, THIS IS A SIGN OF ILL-CONDITIONING, WHICH */
        /*        MAY BE CAUSED BY SEVERAL THINGS.  ONE IS POOR SCALING OF THE */
        /*        PARAMETERS; ANOTHER IS AN UNFORTUNATE INITIAL GUESS FOR THE */
        /*        PARAMETERS, STILL ANOTHER IS A POOR CHOICE OF THE MODEL. */

        /*     ALGORITHM */

        /*           THE RESIDUAL R IS MODIFIED TO INCORPORATE, FOR ANY FIXED */
        /*        ALF, THE OPTIMAL LINEAR PARAMETERS FOR THAT ALF.  IT IS THEN */
        /*        POSSIBLE TO MINIMIZE ONLY ON THE NONLINEAR PARAMETERS.  AFTER */
        /*        THE OPTIMAL VALUES OF THE NONLINEAR PARAMETERS HAVE BEEN DETER- */
        /*        MINED, THE LINEAR PARAMETERS CAN BE RECOVERED BY LINEAR LEAST */
        /*        SQUARES TECHNIQUES (SEE REF. 1). */

        /*           THE MINIMIZATION IS BY A MODIFICATION OF OSBORNE'S (REF. 3) */
        /*        MODIFICATION OF THE LEVENBERG-MARQUARDT ALGORITHM.  INSTEAD OF */
        /*        SOLVING THE NORMAL EQUATIONS WITH MATRIX */

        /*                 T      2 */
        /*               (J J + NU  * D),      WHERE  J = D(ETA)/D(ALF), */

        /*        STABLE ORTHOGONAL (HOUSEHOLDER) REFLECTIONS ARE USED ON A */
        /*        MODIFICATION OF THE MATRIX */
        /*                                   (   J  ) */
        /*                                   (------) , */
        /*                                   ( NU*D ) */

        /*        WHERE D IS A DIAGONAL MATRIX CONSISTING OF THE LENGTHS OF THE */
        /*        COLUMNS OF J.  THIS MARQUARDT STABILIZATION ALLOWS THE ROUTINE */
        /*        TO RECOVER FROM SOME RANK DEFICIENCIES IN THE JACOBIAN. */
        /*        OSBORNE'S EMPIRICAL STRATEGY FOR CHOOSING THE MARQUARDT PARAM- */
        /*        ETER HAS PROVEN REASONABLY SUCCESSFUL IN PRACTICE.  (GAUSS- */
        /*        NEWTON WITH STEP CONTROL CAN BE OBTAINED BY MAKING THE CHANGE */
        /*        INDICATED BEFORE THE INSTRUCTION LABELED 5).  A DESCRIPTION CAN */
        /*        BE FOUND IN REF. (3), AND A FLOW CHART IN (2), P. 22. */

        /*        FOR REFERENCE, SEE */

        /*        1.  GENE H. GOLUB AND V. PEREYRA, 'THE DIFFERENTIATION OF */
        /*            PSEUDO-INVERSES AND NONLINEAR LEAST SQUARES PROBLEMS WHOSE */
        /*            VARIABLES SEPARATE,' SIAM J. NUMER. ANAL. 10, 413-432 */
        /*            (1973). */
        /*        2.  ------, SAME TITLE, STANFORD C.S. REPORT 72-261, FEB. 1972. */
        /*        3.  OSBORNE, MICHAEL R., 'SOME ASPECTS OF NON-LINEAR LEAST */
        /*            SQUARES CALCULATIONS,' IN LOOTSMA, ED., 'NUMERICAL METHODS */
        /*            FOR NON-LINEAR OPTIMIZATION,' ACADEMIC PRESS, LONDON, 1972. */
        /*        4.  KROGH, FRED, 'EFFICIENT IMPLEMENTATION OF A VARIABLE PRO- */
        /*            JECTION ALGORITHM FOR NONLINEAR LEAST SQUARES PROBLEMS,' */
        /*            COMM. ACM 17, PP. 167-169 (MARCH, 1974). */
        /*        5.  KAUFMAN, LINDA, 'A VARIABLE PROJECTION METHOD FOR SOLVING */
        /*            SEPARABLE NONLINEAR LEAST SQUARES PROBLEMS', B.I.T. 15, */
        /*            49-57 (1975). */
        /*        6.  DRAPER, N., AND SMITH, H., APPLIED REGRESSION ANALYSIS, */
        /*            WILEY, N.Y., 1966 (FOR STATISTICAL INFORMATION ONLY). */
        /*        7.  C. LAWSON AND R. HANSON, SOLVING LEAST SQUARES PROBLEMS, */
        /*            PRENTICE-HALL, ENGLEWOOD CLIFFS, N. J., 1974. */

        /*                      JOHN BOLSTAD */
        /*                      COMPUTER SCIENCE DEPT., SERRA HOUSE */
        /*                      STANFORD UNIVERSITY */
        /*                      JANUARY, 1977 */

        /*     .................................................................. */

        /* Parameter adjustments */
        --beta;
        --alf;
        --w;
        --y;
        a_dim1 = *nmax;
        a_offset = 1 + a_dim1 * 1;
        a -= a_offset;
        t_dim1 = *nmax;
        t_offset = 1 + t_dim1 * 1;
        t -= t_offset;

        /* Function Body */

        /*           THE FOLLOWING TWO PARAMETERS ARE USED IN THE CONVERGENCE */
        /*           TEST:  EPS1 IS AN ABSOLUTE AND RELATIVE TOLERANCE FOR THE */
        /*           NORM OF THE PROJECTION OF THE RESIDUAL ONTO THE RANGE OF THE */
        /*           JACOBIAN OF THE VARIABLE PROJECTION FUNCTIONAL. */
        /*           ITMAX IS THE MAXIMUM NUMBER OF FUNCTION AND DERIVATIVE */
        /*           EVALUATIONS ALLOWED.  CAUTION:  EPS1 MUST NOT BE */
        /*           SET SMALLER THAN 10 TIMES THE UNIT ROUND-OFF OF THE MACHINE. */

        /* ----------------------------------------------------------------- */
        /* ALL LIB MONITOR FROM VARPRO, MAINTENANCE NUMBER 509, DATE 77178 */
        /*        CALL LIBMON(509) */
        /* ***PLEASE DON'T REMOVE OR CHANGE THE ABOVE CALL. IT IS YOUR ONLY */
        /* ***PROTECTION AGAINST YOUR USING AN OUT-OF-DATE OR INCORRECT */
        /* ***VERSION OF THE ROUTINE. THE LIBRARY MONITOR REMOVES THIS CALL, */
        /* ***SO IT ONLY OCCURS ONCE, ON THE FIRST ENTRY TO THIS ROUTINE. */
        /* ----------------------------------------------------------------- */
        *ierr = 1;
        iter = 0;
        lp1 = *l + 1;
        b1 = *l + 2;
        lnl2 = *l + *nl + 2;
        nlp1 = *nl + 1;
        skip = FALSE_;
        modit = *iprint;
        if (*iprint <= 0) {
            modit = itmax + 2;
        }
        nu = (float)0.;
        /*              IF GAUSS-NEWTON IS DESIRED REMOVE THE NEXT STATEMENT. */
        nu = (float)1.;

        /*              BEGIN OUTER ITERATION LOOP TO UPDATE ALF. */
        /*              CALCULATE THE NORM OF THE RESIDUAL AND THE DERIVATIVE OF */
        /*              THE MODIFIED RESIDUAL THE FIRST TIME, BUT ONLY THE */
        /*              DERIVATIVE IN SUBSEQUENT ITERATIONS. */

L5:
        vpdpa_(l, nl, n, nmax, lpp2, iv, &t[t_offset], &y[1], &w[1], &alf[1], (
                    U_fp)ada, ierr, iprint, &a[a_offset], &beta[1], &a[lp1 * a_dim1 + 
                1], &r__);
        gnstep = (float)1.;
        iterin = 0;
        if (iter > 0) {
            goto L10;
        }
        if (*nl == 0) {
            goto L90;
        }
        if (*ierr != 1) {
            goto L99;
        }

        if (*iprint <= 0) {
            goto L10;
        }
        cout << " " << iterin << " Norm of residual = " << r__ << endl;
        cout << "    Nu = " << nu << endl;
        /*                              BEGIN TWO-STAGE ORTHOGONAL FACTORIZATION */
L10:
        vpfac1_(&nlp1, nmax, n, l, iprint, &a[b1 * a_dim1 + 1], &prjres, ierr);
        if (*ierr < 0) {
            goto L99;
        }
        *ierr = 2;
        if (nu == (float)0.) {
            goto L30;
        }

        /*              BEGIN INNER ITERATION LOOP FOR GENERATING NEW ALF AND */
        /*              TESTING IT FOR ACCEPTANCE. */

L25:
        vpfac2_(&nlp1, nmax, &nu, &a[b1 * a_dim1 + 1]);

        /*              SOLVE A NL X NL UPPER TRIANGULAR SYSTEM FOR DELTA-ALF. */
        /*              THE TRANSFORMED RESIDUAL (IN COL. LNL2 OF A) IS OVER- */
        /*              WRITTEN BY THE RESULT DELTA-ALF. */

L30:
        vpbsol_(nmax, nl, &a[b1 * a_dim1 + 1], &a[lnl2 * a_dim1 + 1]);
        i__1 = *nl;
        for (k = 1; k <= i__1; ++k) {
            /* L35: */
            a[k + b1 * a_dim1] = alf[k] + a[k + lnl2 * a_dim1];
        }
        /*           NEW ALF(K) = ALF(K) + DELTA ALF(K) */

        /*              STEP TO THE NEW POINT NEW ALF, AND COMPUTE THE NEW */
        /*              NORM OF RESIDUAL.  NEW ALF IS STORED IN COLUMN B1 OF A. */

L40:
        vpdpa_(l, nl, n, nmax, lpp2, iv, &t[t_offset], &y[1], &w[1], &a[b1 * 
                a_dim1 + 1], (U_fp)ada, ierr, iprint, &a[a_offset], &beta[1], &a[
                lp1 * a_dim1 + 1], &rnew);
        if (*ierr != 2) {
            goto L99;
        }
        ++iter;
        ++iterin;
        skip = iter % modit != 0;
        if (skip) {
            goto L45;
        }
        cout << "  Iteration " << iter << "   Nonlinear parameters" << endl;
        i__1 = *nl;
        for (k = 1; k <= i__1; ++k) {
            cout << " " << scientific << setw(15) << setprecision(7)
                << a[k+b1*a_dim1];
            if (0 == k%8) cout << endl;
        }
        cout << endl;
        cout << " " << iterin << " Norm of residual = " << rnew << endl;

L45:
        if (iter < itmax) {
            goto L50;
        }
        *ierr = -1;
        vperr_(iprint, ierr, &c__1);
        goto L95;
L50:
        if (rnew - r__ < eps1 * (r__ + 1.)) {
            goto L75;
        }

        /*              RETRACT THE STEP JUST TAKEN */

        if (nu != (float)0.) {
            goto L60;
        }
        /*                                             GAUSS-NEWTON OPTION ONLY */
        gnstep *= (float).5;
        if (gnstep < eps1) {
            goto L95;
        }
        i__1 = *nl;
        for (k = 1; k <= i__1; ++k) {
            /* L55: */
            a[k + b1 * a_dim1] = alf[k] + gnstep * a[k + lnl2 * a_dim1];
        }
        goto L40;
        /*                                        ENLARGE THE MARQUARDT PARAMETER */
L60:
        nu *= (float)1.5;
        if (! skip) {
            cout << "    Step retracted, nu = "
                << scientific << setprecision(7) << setw(15) << nu << endl;
        }
        if (nu <= (float)100.) {
            goto L65;
        }
        *ierr = -2;
        vperr_(iprint, ierr, &c__1);
        goto L95;
        /*                                        RETRIEVE UPPER TRIANGULAR FORM */
        /*                                        AND RESIDUAL OF FIRST STAGE. */
L65:
        i__1 = *nl;
        for (k = 1; k <= i__1; ++k) {
            ksub = lp1 + k;
            i__2 = nlp1;
            for (j = k; j <= i__2; ++j) {
                jsub = lp1 + j;
                isub = nlp1 + j;
                /* L70: */
                a[k + jsub * a_dim1] = a[isub + ksub * a_dim1];
            }
        }
        goto L25;
        /*                                        END OF INNER ITERATION LOOP */
        /*           ACCEPT THE STEP JUST TAKEN */

L75:
        r__ = rnew;
        i__2 = *nl;
        for (k = 1; k <= i__2; ++k) {
            /* L80: */
            alf[k] = a[k + b1 * a_dim1];
        }
        /*                                        CALC. NORM(DELTA ALF)/NORM(ALF) */
        acum = gnstep * vpnorm_(nl, &a[lnl2 * a_dim1 + 1]) / vpnorm_(nl, &alf[1]);

        /*           IF ITERIN IS GREATER THAN 1, A STEP WAS RETRACTED DURING */
        /*           THIS OUTER ITERATION. */

        if (iterin == 1) {
            nu *= (float).5;
        }
        if (skip) {
            goto L85;
        }
        cout << "     nu = " << scientific << setprecision(7) << setw(15)
            << nu << endl;
        cout << "    Norm(delta-alf) / norm(alf) = " << setprecision(3)
            << setw(12) << acum << endl;
L85:
        *ierr = 3;
        if (prjres > eps1 * (r__ + 1.)) {
            goto L5;
        }
        /*           END OF OUTER ITERATION LOOP */

        /*           CALCULATE FINAL QUANTITIES -- LINEAR PARAMETERS, RESIDUALS, */
        /*           COVARIANCE MATRIX, ETC. */

L90:
        *ierr = iter;
L95:
        if (*nl > 0) {
            vpdpa_(l, nl, n, nmax, lpp2, iv, &t[t_offset], &y[1], &w[1], &alf[1], 
                    (U_fp)ada, &c__4, iprint, &a[a_offset], &beta[1], &a[lp1 * 
                    a_dim1 + 1], &r__);
        }
        vppost_(l, nl, n, nmax, &lnl2, &eps1, &r__, iprint, &alf[1], &w[1], &a[
                a_offset], &a[lp1 * a_dim1 + 1], &beta[1], ierr);
L99:
        return 0;

    } /* varpro_ */


    /* Subroutine */ int vpfac1_(long *nlp1, long *nmax, long *n, 
            long *l, long *iprint, double *b, double *prjres, 
            long *ierr)
    {
        /* System generated locals */
        long b_dim1, b_offset, i__1, i__2, i__3;
        double d__1;

        /* Local variables */
        static long i__, j, k;
        static double u;
        static long nl, kp1, lp1, nl23, lpk;
        static double beta, acum;
        static long jsub;
        static double alpha;
        extern /* Subroutine */ int vperr_(long *, long *, long *);
        extern double vpnorm_(long *, double *);


        /*            STAGE 1:  HOUSEHOLDER REDUCTION OF */

        /*                   (    .    )      ( DR'. R3 )    NL */
        /*                   ( DR . R2 )  TO  (----. -- ), */
        /*                   (    .    )      (  0 . R4 )  N-L-NL */

        /*                     NL    1          NL   1 */

        /*         WHERE DR = -D(Q2)*Y IS THE DERIVATIVE OF THE MODIFIED RESIDUAL */
        /*         PRODUCED BY VPDPA, R2 IS THE TRANSFORMED RESIDUAL FROM DPA, */
        /*         AND DR' IS IN UPPER TRIANGULAR FORM (AS IN REF. (2), P. 18). */
        /*         DR IS STORED IN ROWS L+1 TO N AND COLUMNS L+2 TO L + NL + 1 OF */
        /*         THE MATRIX A (I.E., COLUMNS 1 TO NL OF THE MATRIX B).  R2 IS */
        /*         STORED IN COLUMN L + NL + 2 OF THE MATRIX A (COLUMN NL + 1 OF */
        /*         B).  FOR K = 1, 2, ..., NL, FIND REFLECTION I - U * U' / BETA */
        /*         WHICH ZEROES B(I, K), I = L+K+1, ..., N. */

        /*     .................................................................. */


        /* Parameter adjustments */
        b_dim1 = *nmax;
        b_offset = 1 + b_dim1 * 1;
        b -= b_offset;

        /* Function Body */
        nl = *nlp1 - 1;
        nl23 = (nl << 1) + 3;
        lp1 = *l + 1;

        i__1 = nl;
        for (k = 1; k <= i__1; ++k) {
            lpk = *l + k;
            i__2 = *n + 1 - lpk;
            d__1 = vpnorm_(&i__2, &b[lpk + k * b_dim1]);
            alpha = d_sign(&d__1, &b[lpk + k * b_dim1]);
            u = b[lpk + k * b_dim1] + alpha;
            b[lpk + k * b_dim1] = u;
            beta = alpha * u;
            if (alpha != (float)0.) {
                goto L13;
            }
            /*                                                   COLUMN WAS ZERO */
            *ierr = -8;
            i__2 = lp1 + k;
            vperr_(iprint, ierr, &i__2);
            goto L99;
            /*                                APPLY REFLECTIONS TO REMAINING COLUMNS */
            /*                                OF B AND TO RESIDUAL VECTOR. */
L13:
            kp1 = k + 1;
            i__2 = *nlp1;
            for (j = kp1; j <= i__2; ++j) {
                acum = (float)0.;
                i__3 = *n;
                for (i__ = lpk; i__ <= i__3; ++i__) {
                    /* L20: */
                    acum += b[i__ + k * b_dim1] * b[i__ + j * b_dim1];
                }
                acum /= beta;
                i__3 = *n;
                for (i__ = lpk; i__ <= i__3; ++i__) {
                    /* L25: */
                    b[i__ + j * b_dim1] -= b[i__ + k * b_dim1] * acum;
                }
            }
            /* L30: */
            b[lpk + k * b_dim1] = -alpha;
        }

        *prjres = vpnorm_(&nl, &b[lp1 + *nlp1 * b_dim1]);

        /*           SAVE UPPER TRIANGULAR FORM AND TRANSFORMED RESIDUAL, FOR USE */
        /*           IN CASE A STEP IS RETRACTED.  ALSO COMPUTE COLUMN LENGTHS. */

        if (*ierr == 4) {
            goto L99;
        }
        i__1 = nl;
        for (k = 1; k <= i__1; ++k) {
            lpk = *l + k;
            i__3 = *nlp1;
            for (j = k; j <= i__3; ++j) {
                jsub = *nlp1 + j;
                b[k + j * b_dim1] = b[lpk + j * b_dim1];
                /* L40: */
                b[jsub + k * b_dim1] = b[lpk + j * b_dim1];
            }
            /* L50: */
            b[nl23 + k * b_dim1] = vpnorm_(&k, &b[lp1 + k * b_dim1]);
        }

L99:
        return 0;
    } /* vpfac1_ */


    /* Subroutine */ int vpfac2_(long *nlp1, long *nmax, double *nu, 
            double *b)
    {
        /* System generated locals */
        long b_dim1, b_offset, i__1, i__2, i__3;
        double d__1;

        /* Local variables */
        static long i__, j, k;
        static double u;
        static long nl, nl2, kp1, nl23;
        static double beta, acum;
        static long nlpk;
        static double alpha;
        static long nlpkm1;
        extern double vpnorm_(long *, double *);


        /*        STAGE 2:  SPECIAL HOUSEHOLDER REDUCTION OF */

        /*                      NL     ( DR' . R3 )      (DR'' . R5 ) */
        /*                             (-----. -- )      (-----. -- ) */
        /*                  N-L-NL     (  0  . R4 )  TO  (  0  . R4 ) */
        /*                             (-----. -- )      (-----. -- ) */
        /*                      NL     (NU*D . 0  )      (  0  . R6 ) */

        /*                                NL    1          NL    1 */

        /*        WHERE DR', R3, AND R4 ARE AS IN VPFAC1, NU IS THE MARQUARDT */
        /*        PARAMETER, D IS A DIAGONAL MATRIX CONSISTING OF THE LENGTHS OF */
        /*        THE COLUMNS OF DR', AND DR'' IS IN UPPER TRIANGULAR FORM. */
        /*        DETAILS IN (1), PP. 423-424.  NOTE THAT THE (N-L-NL) BAND OF */
        /*        ZEROES, AND R4, ARE OMITTED IN STORAGE. */

        /*     .................................................................. */


        /* Parameter adjustments */
        b_dim1 = *nmax;
        b_offset = 1 + b_dim1 * 1;
        b -= b_offset;

        /* Function Body */
        nl = *nlp1 - 1;
        nl2 = nl << 1;
        nl23 = nl2 + 3;
        i__1 = nl;
        for (k = 1; k <= i__1; ++k) {
            kp1 = k + 1;
            nlpk = nl + k;
            nlpkm1 = nlpk - 1;
            b[nlpk + k * b_dim1] = *nu * b[nl23 + k * b_dim1];
            b[nl + k * b_dim1] = b[k + k * b_dim1];
            i__2 = k + 1;
            d__1 = vpnorm_(&i__2, &b[nl + k * b_dim1]);
            alpha = d_sign(&d__1, &b[k + k * b_dim1]);
            u = b[k + k * b_dim1] + alpha;
            beta = alpha * u;
            b[k + k * b_dim1] = -alpha;
            /*                        THE K-TH REFLECTION MODIFIES ONLY ROWS K, */
            /*                        NL+1, NL+2, ..., NL+K, AND COLUMNS K TO NL+1. */
            i__2 = *nlp1;
            for (j = kp1; j <= i__2; ++j) {
                b[nlpk + j * b_dim1] = (float)0.;
                acum = u * b[k + j * b_dim1];
                i__3 = nlpkm1;
                for (i__ = *nlp1; i__ <= i__3; ++i__) {
                    /* L20: */
                    acum += b[i__ + k * b_dim1] * b[i__ + j * b_dim1];
                }
                acum /= beta;
                b[k + j * b_dim1] -= u * acum;
                i__3 = nlpk;
                for (i__ = *nlp1; i__ <= i__3; ++i__) {
                    /* L30: */
                    b[i__ + j * b_dim1] -= b[i__ + k * b_dim1] * acum;
                }
            }
        }

        return 0;
    } /* vpfac2_ */


    /* Subroutine */ int vpdpa_(long *l, long *nl, long *n, long *
            nmax, long *lpp2, long *iv, double *t, double *y, 
            double *w, double *alf, S_fp ada, long *isel, long *
            iprint, double *a, double *u, double *r__, double *
            rnorm)
    {
        /* System generated locals */
        long a_dim1, a_offset, t_dim1, t_offset, i__1, i__2, i__3;
        double d__1;

        /* Local variables */
        static long i__, j, k, m, kp1, lp1, lp2, inc[96]	/* was [12][8] */, 
                    lnl2, lpp1;
        static double beta, acum;
        static long ncon;
        static double save;
        static long ksub;
        static double alpha;
        static long lastc;
        extern /* Subroutine */ int vperr_(long *, long *, long *);
        static bool philp1;
        static long nconp1, firstc;
        static bool nowate;
        extern /* Subroutine */ int vpbsol_(long *, long *, double *, 
                double *);
        static long firstr;
        extern /* Subroutine */ int vpinit_(long *, long *, long *, 
                long *, long *, long *, double *, double *, 
                double *, S_fp, long *, long *, double *, long *,
                long *, long *, bool *, bool *);
        extern double vpnorm_(long *, double *);


        /*        COMPUTE THE NORM OF THE RESIDUAL (IF ISEL = 1 OR 2), OR THE */
        /*        (N-L) X NL DERIVATIVE OF THE MODIFIED RESIDUAL (N-L) VECTOR */
        /*        Q2*Y (IF ISEL = 1 OR 3).  HERE Q * PHI = S, I.E., */

        /*         L     ( Q1 ) (     .   .        )   ( S  . R1 .  F1  ) */
        /*               (----) ( PHI . Y . D(PHI) ) = (--- . -- . ---- ) */
        /*         N-L   ( Q2 ) (     .   .        )   ( 0  . R2 .  F2  ) */

        /*                 N       L    1      P         L     1     P */

        /*        WHERE Q IS N X N ORTHOGONAL, AND S IS L X L UPPER TRIANGULAR. */
        /*        THE NORM OF THE RESIDUAL = NORM(R2), AND THE DESIRED DERIVATIVE */
        /*        ACCORDING TO REF. (5), IS */
        /*                                               -1 */
        /*                    D(Q2 * Y) = -Q2 * D(PHI)* S  * Q1* Y. */

        /*     .................................................................. */


        /* Parameter adjustments */
        --u;
        --alf;
        --r__;
        --w;
        --y;
        a_dim1 = *nmax;
        a_offset = 1 + a_dim1 * 1;
        a -= a_offset;
        t_dim1 = *nmax;
        t_offset = 1 + t_dim1 * 1;
        t -= t_offset;

        /* Function Body */
        if (*isel != 1) {
            goto L3;
        }
        lp1 = *l + 1;
        lnl2 = *l + 2 + *nl;
        lp2 = *l + 2;
        lpp1 = *lpp2 - 1;
        firstc = 1;
        lastc = lpp1;
        firstr = lp1;
        vpinit_(l, nl, n, nmax, lpp2, iv, &t[t_offset], &w[1], &alf[1], (S_fp)ada,
                isel, iprint, &a[a_offset], inc, &ncon, &nconp1, &philp1, &
                nowate);
        if (*isel != 1) {
            goto L99;
        }
        goto L30;

L3:
        i__1 = min(*isel,3L);
        (*ada)(&lp1, nl, n, nmax, lpp2, iv, &a[a_offset], inc, &t[t_offset], &alf[
               1], &i__1);
        if (*isel == 2) {
            goto L6;
        }
        /*                                                 ISEL = 3 OR 4 */
        firstc = lp2;
        lastc = lpp1;
        firstr = (4 - *isel) * *l + 1;
        goto L50;
        /*                                                 ISEL = 2 */
L6:
        firstc = nconp1;
        lastc = lp1;
        if (ncon == 0) {
            goto L30;
        }
        if (a[ncon * a_dim1 + 1] == save) {
            goto L30;
        }
        cout << "NCON " << isel << endl;
        *isel = -7;
        vperr_(iprint, isel, &ncon);
        goto L99;
        /*                                                  ISEL = 1 OR 2 */
L30:
        if (philp1) {
            goto L40;
        }
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            /* L35: */
            r__[i__] = y[i__];
        }
        goto L50;
L40:
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            /* L45: */
            r__[i__] = y[i__] - r__[i__];
        }
        /*                                             WEIGHT APPROPRIATE COLUMNS */
L50:
        if (nowate) {
            goto L58;
        }
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            acum = w[i__];
            i__2 = lastc;
            for (j = firstc; j <= i__2; ++j) {
                /* L55: */
                a[i__ + j * a_dim1] *= acum;
            }
        }

        /*           COMPUTE ORTHOGONAL FACTORIZATIONS BY HOUSEHOLDER */
        /*           REFLECTIONS.  IF ISEL = 1 OR 2, REDUCE PHI (STORED IN THE */
        /*           FIRST L COLUMNS OF THE MATRIX A) TO UPPER TRIANGULAR FORM, */
        /*           (Q*PHI = S), AND TRANSFORM Y (STORED IN COLUMN L+1), GETTING */
        /*           Q*Y = R.  IF ISEL = 1, ALSO TRANSFORM J = D PHI (STORED IN */
        /*           COLUMNS L+2 THROUGH L+P+1 OF THE MATRIX A), GETTING Q*J = F. */
        /*           IF ISEL = 3 OR 4, PHI HAS ALREADY BEEN REDUCED, TRANSFORM */
        /*           ONLY J.  S, R, AND F OVERWRITE PHI, Y, AND J, RESPECTIVELY, */
        /*           AND A FACTORED FORM OF Q IS SAVED IN U AND THE LOWER */
        /*           TRIANGLE OF PHI. */

L58:
        if (*l == 0) {
            goto L75;
        }
        i__2 = *l;
        for (k = 1; k <= i__2; ++k) {
            kp1 = k + 1;
            if (*isel >= 3 || *isel == 2 && k < nconp1) {
                goto L66;
            }
            i__1 = *n + 1 - k;
            d__1 = vpnorm_(&i__1, &a[k + k * a_dim1]);
            alpha = d_sign(&d__1, &a[k + k * a_dim1]);
            u[k] = a[k + k * a_dim1] + alpha;
            a[k + k * a_dim1] = -alpha;
            firstc = kp1;
            if (alpha != (float)0.) {
                goto L66;
            }
            *isel = -8;
            vperr_(iprint, isel, &k);
            goto L99;
            /*                                        APPLY REFLECTIONS TO COLUMNS */
            /*                                        FIRSTC TO LASTC. */
L66:
            beta = -a[k + k * a_dim1] * u[k];
            i__1 = lastc;
            for (j = firstc; j <= i__1; ++j) {
                acum = u[k] * a[k + j * a_dim1];
                i__3 = *n;
                for (i__ = kp1; i__ <= i__3; ++i__) {
                    /* L68: */
                    acum += a[i__ + k * a_dim1] * a[i__ + j * a_dim1];
                }
                acum /= beta;
                a[k + j * a_dim1] -= u[k] * acum;
                i__3 = *n;
                for (i__ = kp1; i__ <= i__3; ++i__) {
                    /* L70: */
                    a[i__ + j * a_dim1] -= a[i__ + k * a_dim1] * acum;
                }
            }
        }

L75:
        if (*isel >= 3) {
            goto L85;
        }
        i__3 = *n - *l;
        *rnorm = vpnorm_(&i__3, &r__[lp1]);
        if (*isel == 2) {
            goto L99;
        }
        if (ncon > 0) {
            save = a[ncon * a_dim1 + 1];
        }

        /*           F2 IS NOW CONTAINED IN ROWS L+1 TO N AND COLUMNS L+2 TO */
        /*           L+P+1 OF THE MATRIX A.  NOW SOLVE THE L X L UPPER TRIANGULAR */
        /*           SYSTEM S*BETA = R1 FOR THE LINEAR PARAMETERS BETA.  BETA */
        /*           OVERWRITES R1. */

L85:
        if (*l > 0) {
            vpbsol_(nmax, l, &a[a_offset], &r__[1]);
        }

        /*           MAJOR PART OF KAUFMAN'S SIMPLIFICATION OCCURS HERE.  COMPUTE */
        /*           THE DERIVATIVE OF ETA WITH RESPECT TO THE NONLINEAR */
        /*           PARAMETERS */

        /*   T   D ETA        T    L          D PHI(J)    D PHI(L+1) */
        /*  Q * --------  =  Q * (SUM BETA(J) --------  + ----------)  =  F2*BETA */
        /*      D ALF(K)          J=1         D ALF(K)     D ALF(K) */

        /*           AND STORE THE RESULT IN COLUMNS L+2 TO L+NL+1.  IF ISEL NOT */
        /*           = 4, THE FIRST L ROWS ARE OMITTED.  THIS IS -D(Q2)*Y.  IF */
        /*           ISEL NOT = 4 THE RESIDUAL R2 = Q2*Y (IN COL. L+1) IS COPIED */
        /*           TO COLUMN L+NL+2.  OTHERWISE ALL OF COLUMN L+1 IS COPIED. */

        i__3 = *n;
        for (i__ = firstr; i__ <= i__3; ++i__) {
            if (*l == ncon) {
                goto L95;
            }
            m = lp1;
            i__1 = *nl;
            for (k = 1; k <= i__1; ++k) {
                acum = (float)0.;
                i__2 = *l;
                for (j = nconp1; j <= i__2; ++j) {
                    if (inc[k + j * 12 - 13] == 0) {
                        goto L88;
                    }
                    ++m;
                    acum += a[i__ + m * a_dim1] * r__[j];
L88:
                    ;
                }
                ksub = lp1 + k;
                if (inc[k + lp1 * 12 - 13] == 0) {
                    goto L90;
                }
                ++m;
                acum += a[i__ + m * a_dim1];
L90:
                a[i__ + ksub * a_dim1] = acum;
            }
L95:
            a[i__ + lnl2 * a_dim1] = r__[i__];
        }

L99:
        return 0;
    } /* vpdpa_ */

    /* CC======================================================================== */
    /* Subroutine */ int vpinit_(long *l, long *nl, long *n, long *
            nmax, long *lpp2, long *iv, double *t, double *w, 
            double *alf, S_fp ada, long *isel, long *iprint, double 
            *a, long *inc, long *ncon, long *nconp1, bool *philp1, 
            bool *nowate)
    {
        /* Initialized data */

        /* System generated locals */
        long a_dim1, a_offset, t_dim1, t_offset, i__1, i__2;

        /* Local variables */
        static long i__, j, k, p, lp1, lnl2, inckj;
        extern /* Subroutine */ int vperr_(long *, long *, long *);


        /*        CHECK VALIDITY OF INPUT PARAMETERS, AND DETERMINE NUMBER OF */
        /*        CONSTANT FUNCTIONS. */

        /*     .................................................................. */

        /* Parameter adjustments */
        --alf;
        --w;
        a_dim1 = *nmax;
        a_offset = 1 + a_dim1 * 1;
        a -= a_offset;
        t_dim1 = *nmax;
        t_offset = 1 + t_dim1 * 1;
        t -= t_offset;
        inc -= 13;

        /* Function Body */

        lp1 = *l + 1;
        lnl2 = *l + 2 + *nl;
        /*                                          CHECK FOR VALID INPUT */
        if (*l >= 0 && *nl >= 0 && *l + *nl < *n && lnl2 <= *lpp2 && (*nl << 1) + 
                3 <= *nmax && *n <= *nmax && *iv > 0 && ! (*nl == 0 && *l == 0)) {
            goto L1;
        }
        *isel = -4;
        vperr_(iprint, isel, &c__1);
        goto L99;

L1:
        if (*l == 0 || *nl == 0) {
            goto L3;
        }
        i__1 = lp1;
        for (j = 1; j <= i__1; ++j) {
            i__2 = *nl;
            for (k = 1; k <= i__2; ++k) {
                /* L2: */
                inc[k + j * 12] = 0;
            }
        }

L3:
        (*ada)(&lp1, nl, n, nmax, lpp2, iv, &a[a_offset], &inc[13], &t[t_offset], 
               &alf[1], isel);

        *nowate = TRUE_;
        i__2 = *n;
        for (i__ = 1; i__ <= i__2; ++i__) {
            *nowate = *nowate && w[i__] == (float)1.;
            if (w[i__] >= (float)0.) {
                goto L9;
            }
            /*                                                ERROR IN WEIGHTS */
            *isel = -6;
            vperr_(iprint, isel, &i__);
            goto L99;
L9:
            w[i__] = sqrt(w[i__]);
        }

        *ncon = *l;
        *nconp1 = lp1;
        *philp1 = *l == 0;
        if (*philp1 || *nl == 0) {
            goto L99;
        }
        /*                                   CHECK INC MATRIX FOR VALID INPUT AND */
        /*                                   DETERMINE NUMBER OF CONSTANT FCNS. */
        p = 0;
        i__2 = lp1;
        for (j = 1; j <= i__2; ++j) {
            if (p == 0) {
                *nconp1 = j;
            }
            i__1 = *nl;
            for (k = 1; k <= i__1; ++k) {
                inckj = inc[k + j * 12];
                if (inckj != 0 && inckj != 1) {
                    goto L15;
                }
                if (inckj == 1) {
                    ++p;
                }
                /* L11: */
            }
        }

        *ncon = *nconp1 - 1;
        if (*iprint >= 0) {
            cout << "  Number of constant functions = " << *ncon << endl;
        }
        if (*l + p + 2 == *lpp2) {
            goto L20;
        }
        /*                                              INPUT ERROR IN INC MATRIX */
L15:
        *isel = -5;
        vperr_(iprint, isel, &c__1);
        goto L99;
        /*                                 DETERMINE IF PHI(L+1) IS IN THE MODEL. */
L20:
        i__1 = *nl;
        for (k = 1; k <= i__1; ++k) {
            /* L25: */
            if (inc[k + lp1 * 12] == 1) {
                *philp1 = TRUE_;
            }
        }

L99:
        return 0;
    } /* vpinit_ */

    /* ccc=========================================================================================== */
    /* Subroutine */ int vpbsol_(long *nmax, long *n, double *a, 
            double *x)
    {
        /* System generated locals */
        long a_dim1, a_offset, i__1, i__2;

        /* Local variables */
        static long i__, j, ip1, np1;
        static double acum;
        static long iback;


        /*        BACKSOLVE THE N X N UPPER TRIANGULAR SYSTEM A*X = B. */
        /*        THE SOLUTION X OVERWRITES THE RIGHT SIDE B. */


        /* Parameter adjustments */
        --x;
        a_dim1 = *nmax;
        a_offset = 1 + a_dim1 * 1;
        a -= a_offset;

        /* Function Body */
        x[*n] /= a[*n + *n * a_dim1];
        if (*n == 1) {
            goto L30;
        }
        np1 = *n + 1;
        i__1 = *n;
        for (iback = 2; iback <= i__1; ++iback) {
            i__ = np1 - iback;
            /*           I = N-1, N-2, ..., 2, 1 */
            ip1 = i__ + 1;
            acum = x[i__];
            i__2 = *n;
            for (j = ip1; j <= i__2; ++j) {
                /* L10: */
                acum -= a[i__ + j * a_dim1] * x[j];
            }
            /* L20: */
            x[i__] = acum / a[i__ + i__ * a_dim1];
        }

L30:
        return 0;
    } /* vpbsol_ */

    /* Subroutine */ int vppost_(long *l, long *nl, long *n, long *
            nmax, long *lnl2, double *eps, double *rnorm, long *
            iprint, double *alf, double *w, double *a, double *
            r__, double *u, long *ierr)
    {
        /* System generated locals */
        long a_dim1, a_offset, i__1, i__2;
        double d__1;

        /* Local variables */
        static long i__, j, k, kp1, lp1, lnl1;
        static double acum, save;
        static long lpnl, kback;
        extern /* Subroutine */ int vpcov_(long *, long *, double *, 
                double *), vpfac1_(long *, long *, long *, long *,
                    long *, double *, double *, long *);
        static double prjres;


        /*        CALCULATE RESIDUALS, SAMPLE VARIANCE, AND COVARIANCE MATRIX. */
        /*        ON INPUT, U CONTAINS INFORMATION ABOUT HOUSEHOLDER REFLECTIONS */
        /*        FROM VPDPA.  ON OUTPUT, IT CONTAINS THE LINEAR PARAMETERS. */

        /* Parameter adjustments */
        --u;
        --alf;
        --r__;
        --w;
        a_dim1 = *nmax;
        a_offset = 1 + a_dim1 * 1;
        a -= a_offset;

        /* Function Body */

        lp1 = *l + 1;
        lpnl = *lnl2 - 2;
        lnl1 = lpnl + 1;
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            /* L10: */
            /* Computing 2nd power */
            d__1 = w[i__];
            w[i__] = d__1 * d__1;
        }

        /*              UNWIND HOUSEHOLDER TRANSFORMATIONS TO GET RESIDUALS, */
        /*              AND MOVE THE LINEAR PARAMETERS FROM R TO U. */

        if (*l == 0) {
            goto L30;
        }
        i__1 = *l;
        for (kback = 1; kback <= i__1; ++kback) {
            k = lp1 - kback;
            kp1 = k + 1;
            acum = (float)0.;
            i__2 = *n;
            for (i__ = kp1; i__ <= i__2; ++i__) {
                /* L20: */
                acum += a[i__ + k * a_dim1] * r__[i__];
            }
            save = r__[k];
            r__[k] = acum / a[k + k * a_dim1];
            acum = -acum / (u[k] * a[k + k * a_dim1]);
            u[k] = save;
            i__2 = *n;
            for (i__ = kp1; i__ <= i__2; ++i__) {
                /* L25: */
                r__[i__] -= a[i__ + k * a_dim1] * acum;
            }
        }
        /*                                              COMPUTE MEAN ERROR */
L30:
        acum = (float)0.;
        i__2 = *n;
        for (i__ = 1; i__ <= i__2; ++i__) {
            /* L35: */
            acum += r__[i__];
        }
        save = acum / *n;

        /*              THE FIRST L COLUMNS OF THE MATRIX HAVE BEEN REDUCED TO */
        /*              UPPER TRIANGULAR FORM IN VPDPA.  FINISH BY REDUCING ROWS */
        /*              L+1 TO N AND COLUMNS L+2 THROUGH L+NL+1 TO TRIANGULAR */
        /*              FORM.  THEN SHIFT COLUMNS OF DERIVATIVE MATRIX OVER ONE */
        /*              TO THE LEFT TO BE ADJACENT TO THE FIRST L COLUMNS. */

        if (*nl == 0) {
            goto L45;
        }
        i__2 = *nl + 1;
        vpfac1_(&i__2, nmax, n, l, iprint, &a[(*l + 2) * a_dim1 + 1], &prjres, &
                c__4);
        i__2 = *n;
        for (i__ = 1; i__ <= i__2; ++i__) {
            a[i__ + *lnl2 * a_dim1] = r__[i__];
            i__1 = lnl1;
            for (k = lp1; k <= i__1; ++k) {
                /* L40: */
                a[i__ + k * a_dim1] = a[i__ + (k + 1) * a_dim1];
            }
        }
        /*                                              COMPUTE COVARIANCE MATRIX */
L45:
        a[*lnl2 * a_dim1 + 1] = *rnorm;
        acum = *rnorm * *rnorm / (*n - *l - *nl);
        a[*lnl2 * a_dim1 + 2] = acum;
        vpcov_(nmax, &lpnl, &acum, &a[a_offset]);

        if (*iprint < 0) {
            goto L99;
        }
        cout << " ''''''''''''''''''''''''''''''''''''''''''''''''''" << endl;
        if (*l > 0) {
            i__1 = *l;
            cout << "  Linear parameters // ";
            for (j = 1; j <= i__1; ++j) {
                cout << scientific << setprecision(7) << setw(15)
                    << u[j];
                if (0 == j%8) cout << endl;
            }
        }
        if (*nl > 0) {
            cout << "  Nonlinear parameters // ";
            i__1 = *nl;
            for (k = 1; k <= i__1; ++k) {
                cout << scientific << setprecision(7) << setw(15)
                    << alf[k];
                if (0 == k%8) cout << endl;
            }
        }
        cout << "  Norm of residual = " << scientific << setprecision(7)
            << setw(15) << *rnorm << ", expected error of observations = "
            << save << "\n  estimated variance of observations = "
            << acum << endl;
        if (abs(save) > *eps) {
            cout << " WARNING -- expected error of observations is not zero."
                << " Covariance matrix may be meaningless" << endl;
        }
        cout << " ''''''''''''''''''''''''''''''''''''''''''''''''''" << endl;
L99:
        return 0;

    } /* vppost_ */

    /* C=============================================================================== */
    /* Subroutine */ int vpcov_(long *nmax, long *n, double *sigma2, 
            double *a)
    {
        /* System generated locals */
        long a_dim1, a_offset, i__1, i__2, i__3;

        /* Local variables */
        static long i__, j, m, jm1, ip1, nm1;
        static double sum;


        /*           COMPUTE THE SCALED COVARIANCE MATRIX OF THE L + NL */
        /*        PARAMETERS.  THIS INVOLVES COMPUTING */

        /*                               2     -1    -T */
        /*                          SIGMA  *  T   * T */

        /*        WHERE THE (L+NL) X (L+NL) UPPER TRIANGULAR MATRIX T IS */
        /*        DESCRIBED IN SUBROUTINE VPPOST.  THE RESULT OVERWRITES THE */
        /*        FIRST L+NL ROWS AND COLUMNS OF THE MATRIX A.  THE RESULTING */
        /*        MATRIX IS SYMMETRIC.  SEE REF. 7, PP. 67-70, 281. */

        /*     .................................................................. */


        /* Parameter adjustments */
        a_dim1 = *nmax;
        a_offset = 1 + a_dim1 * 1;
        a -= a_offset;

        /* Function Body */
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            /* L10: */
            a[j + j * a_dim1] = (float)1. / a[j + j * a_dim1];
        }

        /*                 INVERT T UPON ITSELF */

        if (*n == 1) {
            goto L70;
        }
        nm1 = *n - 1;
        i__1 = nm1;
        for (i__ = 1; i__ <= i__1; ++i__) {
            ip1 = i__ + 1;
            i__2 = *n;
            for (j = ip1; j <= i__2; ++j) {
                jm1 = j - 1;
                sum = (float)0.;
                i__3 = jm1;
                for (m = i__; m <= i__3; ++m) {
                    /* L50: */
                    sum += a[i__ + m * a_dim1] * a[m + j * a_dim1];
                }
                /* L60: */
                a[i__ + j * a_dim1] = -sum * a[j + j * a_dim1];
            }
        }

        /*                 NOW FORM THE MATRIX PRODUCT */

L70:
        i__2 = *n;
        for (i__ = 1; i__ <= i__2; ++i__) {
            i__1 = *n;
            for (j = i__; j <= i__1; ++j) {
                sum = (float)0.;
                i__3 = *n;
                for (m = j; m <= i__3; ++m) {
                    /* L80: */
                    sum += a[i__ + m * a_dim1] * a[j + m * a_dim1];
                }
                sum *= *sigma2;
                a[i__ + j * a_dim1] = sum;
                /* L90: */
                a[j + i__ * a_dim1] = sum;
            }
        }

        return 0;
    } /* vpcov_ */

    /* CC======================================================== */
    /* Subroutine */ int vperr_(long *iprint, long *ierr, long *k)
    {
        /* Local variables */
        //static long errno;


        /*        PRINT ERROR MESSAGES */


        if (*iprint < 0) {
            goto L99;
        }
        errno = abs(*ierr);
        switch (errno) {
            case 1:  goto L1;
            case 2:  goto L2;
            case 3:  goto L99;
            case 4:  goto L4;
            case 5:  goto L5;
            case 6:  goto L6;
            case 7:  goto L7;
            case 8:  goto L8;
        }

L1:
        cout << "Problem terminated for excessive iterations\n\n" << endl;
        goto L99;
L2:
        cout << "Problem terminated because of ill-conditioning\n\n" << endl;
        goto L99;
L4:
        cout << "Input error in parameter L, NL, N. Lpp2, or NMAX.\n" << endl;
        goto L99;
L5:
        cout << "Error -- inc matrix improperly specified, "
            << "or disagrees with lpp2.\n" << endl;
        goto L99;
L6:
        cout << "Error -- weight(" << *k << ") is negative.\n" << endl;
        goto L99;
L7:
        cout << "Error -- constant column " << *k 
            << "must be computed only when isel = 1.\n" << endl;
        goto L99;
L8:
        cout << "Catastrophic failure -- column " << *k 
            << "is zero, see documentation.\n" << endl;

L99:
        return 0;
    } /* vperr_ */

    /* CC================================================================================= */
    double vpnorm_(long *n, double *x)
    {
        /* System generated locals */
        long i__1;
        double ret_val, d__1, d__2;

        /* Local variables */
        static long i__;
        static double sum, rmax, term;


        /*        COMPUTE THE L2 (EUCLIDEAN) NORM OF A VECTOR, MAKING SURE TO */
        /*        AVOID UNNECESSARY UNDERFLOWS.  NO ATTEMPT IS MADE TO SUPPRESS */
        /*        OVERFLOWS. */


        /*           FIND LARGEST (IN ABSOLUTE VALUE) ELEMENT */
        /* Parameter adjustments */
        --x;

        /* Function Body */
        rmax = (float)0.;
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            if ((d__1 = x[i__], abs(d__1)) > rmax) {
                rmax = (d__2 = x[i__], abs(d__2));
            }
            /* L10: */
        }

        sum = (float)0.;
        if (rmax == (float)0.) {
            goto L30;
        }
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
            term = (float)0.;
            if (rmax + (d__1 = x[i__], abs(d__1)) != rmax) {
                term = x[i__] / rmax;
            }
            /* L20: */
            sum += term * term;
        }

L30:
        ret_val = rmax * sqrt(sum);
        /* L99: */
        return ret_val;
    } /* vpnorm_ */
}

namespace EMAN {
    long tflm(long l, long nl, long n, long nmax, long lpp2, long iv,
            double t, double y, double w, double a, long iprint,
            double alf, double beta) {
        long ierr;
        varpro_(&l, &nl, &n, &nmax, &lpp2, &iv, &t, &y, &w,
                (U_fp)ada_, &a, &iprint, &alf, &beta, &ierr);
        return ierr;
    }
}

